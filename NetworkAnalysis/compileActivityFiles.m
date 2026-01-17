function [] = compileActivityFiles(data)
% Compiles distributions and summary metrics into a single CSV file
% Processes HD-MEA ActivityScan .h5 files, extracts spike-based features, and saves summary montage


% Setup and extract variables
mfilename('fullpath');
tic;
fig = data.fig;
logFile = data.logFile;
plotFig = data.plotFig;
refDir = data.refDir;
opDir = data.opDir;
div0 = data.div0Date;
parentFolderPath = data.parentFolderPath;
d = uiprogressdlg(fig,'Title','Compiling files','Message','Start','Cancelable','on');

refTable = readtable(refDir);
run_ids = unique(refTable.Run_);
div0_date = datetime(div0, "InputFormat",'MM/dd/yyyy');

% Output folder
if ~isfolder(fullfile(opDir, 'ActivityScan_outputs'))
    mkdir(fullfile(opDir, 'ActivityScan_outputs'));
end

% Gather files
filePattern = fullfile(parentFolderPath, '**/ActivityScan/**/*raw.h5'); 
theFiles = dir(filePattern);
numFiles = length(theFiles);
results = cell(numFiles, 1);
skippedFiles = cell(numFiles, 1); 

% Use parallel pool
numCores = 6;
currentPool = gcp('nocreate');
if isempty(currentPool) || currentPool.NumWorkers ~= numCores
    if ~isempty(currentPool), delete(currentPool); end
    poolObj = parpool(numCores);
end

for k = 1:numFiles 
   % if k > 1, continue;end    %for debugging
    if d.CancelRequested, d.close(); error("User interruption"); end
    d.Value = k/numFiles;
    pathFileActivityScan = fullfile(theFiles(k).folder);
    fileDirParts = strsplit(pathFileActivityScan, filesep);
    scan_runID = str2double(fileDirParts{end});
    scan_chipID = fileDirParts{end-2};
    scan_runID_text = fileDirParts{end};

    if ~ismember(scan_runID, run_ids), continue; end
    idx = refTable.Run_ == scan_runID;
    wellsIDs = str2double(strsplit(strtrim(string(refTable.Wells_Recorded(idx))), ','));
    neuronTypes = strtrim(strsplit(string(refTable.NeuronSource(idx)), ','));
    assayColumn = refTable.Assay(idx);

    if length(wellsIDs) ~= length(neuronTypes)
        error('Mismatch: Each well must have a corresponding neuron type.');
    end

    numWells = length(wellsIDs);
    fileResults = cell(numWells, 1);
    skippedWells = cell(numWells, 1);

    parfor z = 1:numWells
        %if z>1,continue;end      %for testing purpose
        wellID = wellsIDs(z);
        neuronType = neuronTypes(z);
        try
            activityScanData = mxw.fileManager(pathFileActivityScan, wellID);
        catch
            skippedWells{z} = pathFileActivityScan;
            continue;
        end
        msg = sprintf('[Run %d | Well %d] Processing...\n', scan_runID, wellID);
        %fprintf(logFile, '%s', msg);   % Write to log file
        fprintf(1, '%s', msg);     
        
        try
                     
            hd5_time = activityScanData.fileObj.stopTime;
            try
                hd5Date = datetime(hd5_time,'InputFormat', 'yyyy-MM-dd HH:mm:ss');
            catch
                hd5Date = datetime(hd5_time,'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
            end
            div0_datetime = datetime(div0_date, 'InputFormat', 'yourInputFormatHere');
            scan_div = floor(days(hd5Date - div0_date));
            meanFiringRate = mxw.activityMap.computeSpikeRate(activityScanData);
            amplitude90perc = abs(mxw.activityMap.computeAmplitude90percentile(activityScanData));

            % ISI metrics
            nElectrodes = length(activityScanData.extractedSpikes(1).frameno);
            nRecs = length(activityScanData.extractedSpikes);
            spikeTimesCell = cell(nElectrodes, 1);
            
            for e = 1:nElectrodes
                % Collect all frame numbers from all recordings for this electrode
                allFrames = cell(nRecs, 1);
                for r = 1:nRecs
                    allFrames{r} = activityScanData.extractedSpikes(r).frameno{e};
                end
                % Concatenate once and convert to seconds
                spikeTimesCell{e} = vertcat(allFrames{:}) / activityScanData.fileObj(1).samplingFreq;
            end

            % Preallocate result arrays
            meanISI = NaN(nElectrodes,1);
            stdISI = NaN(nElectrodes,1);
            isiCV = NaN(nElectrodes,1);
            fano = NaN(nElectrodes,1); % Assuming fano factor is calculated elsewhere
            
            % Calculate Inter-Spike Intervals (ISIs)
            isiCell = cellfun(@diff, spikeTimesCell, 'UniformOutput', false);
            isiCell = cellfun(@double, isiCell, 'UniformOutput', false);
            % --- Calculate Mean ISI ---
            % A mean can be calculated on one or more ISIs.
            validMean = cellfun(@(x) isnumeric(x) && numel(x) >= 1, isiCell);
            meanISI(validMean) = cellfun(@mean, isiCell(validMean));
            
            % --- Calculate Standard Deviation of ISI ---
            % A standard deviation requires at least two ISIs to be meaningful.
            % This new filter prevents the error and ensures statistical validity.
            validStd = cellfun(@(x) numel(x) > 1, isiCell);
            stdISI(validStd) = cellfun(@std, isiCell(validStd));
            
            % --- Calculate Coefficient of Variation (CV) of ISI ---
            % The CV should only be calculated where a valid standard deviation and mean exist.
            % We use validStd to ensure we don't perform NaN/value or 0/0.
            isiCV(validStd) = stdISI(validStd) ./ meanISI(validStd);


            validSpikes = cellfun(@(x) numel(x) >= 1, spikeTimesCell);
            maxTime = max(cellfun(@(x) max(x), spikeTimesCell(validSpikes)));
            binEdges = 0:1:ceil(maxTime);
            for e = find(validSpikes)'
                st = spikeTimesCell{e};
                sc = histcounts(st, binEdges);
                fano(e) = var(double(sc)) / max(eps, mean(sc));
            end

            % Summary Metrics
            thrFiringRate = 0.1;
            thrAmp = 20;
            scan_meanFiringRate = mean(meanFiringRate(meanFiringRate > thrFiringRate));
            scan_meanSpikeAmplitude = mean(amplitude90perc(amplitude90perc > thrAmp));
            idx = (meanFiringRate > thrFiringRate & amplitude90perc > thrAmp);
            Active_area = sum(idx) / length(idx) * 100;
            xpos_elec = activityScanData.processedMap.xpos(idx);
            ypos_elec= activityScanData.processedMap.ypos(idx);

            summaryRow = table(scan_runID, assayColumn, scan_div, wellID, {neuronType}, hd5Date, {string(scan_chipID)}, ...
                scan_meanFiringRate, scan_meanSpikeAmplitude, Active_area, ...
                nanmean(meanISI), nanmean(stdISI), nanmean(isiCV), nanmean(fano), ...
                'VariableNames', {'Run_ID','AssayType','DIV','Well','NeuronType','Time','Chip_ID', ...
                                  'Mean_FiringRate','Mean_SpikeAmplitude','Active_area', ...
                                  'Mean_ISI','Std_ISI','Mean_ISI_CV','Mean_Fano'});

         
            % This contains the large 26400x1 vectors inside cells
            nestedDistTable = table(scan_runID,assayColumn, scan_div, wellID, {neuronType}, {string(scan_chipID)}, ... % Add identifiers
                                    {meanFiringRate(:)}, {amplitude90perc(:)}, ...
                                    {meanISI(:)}, {stdISI(:)}, {isiCV(:)}, {fano(:)}, ...
                                    'VariableNames', {'Run_ID','AssayType','DIV','Well','NeuronType','Chip_ID', ...
                                                      'FR_Hz','Amplitude_uV','ISI_s','ISI_STD','ISI_CV','Fano'});

            %fullTable = [summaryRow array2table(nan(height(summaryRow), width(distTable) - width(summaryRow))) ; distTable];
           
            fileResults{z} =struct('summary', summaryRow, 'nested', nestedDistTable);


            % --- Plotting Montage ---
            if plotFig
                %frLim = [0.1, mxw.util.percentile(meanFiringRate(meanFiringRate>0),99)];
                %ampLim = [10, mxw.util.percentile(amplitude90perc(amplitude90perc>0),99)];
                frLim = [0.1, 15];
                ampLim = [10, 250];  %hard coding
                aaLim = [0, 5];
                % -- Add plot title-- (SS)
                plotTitle = sprintf('Activity_Run%s_Chip%s_Well%d_DIV%d_%s', ...
                    scan_runID_text, scan_chipID, wellID, scan_div, strrep(neuronType, ' ', '_'));
                plotSummaryMontage(activityScanData, meanFiringRate, amplitude90perc, xpos_elec, ypos_elec, ...
                    meanISI, scan_runID_text, scan_chipID, wellID, scan_div, neuronType, ...
                    fullfile(opDir, 'ActivityScan_outputs'), frLim, ampLim, aaLim, plotTitle);
            end

        catch ME
                % Create a detailed error message
                errorReport = sprintf('ERROR in Run %d, Well %d:\n', scan_runID, wellID);
                errorReport = [errorReport, sprintf('Message: %s\n', ME.message)];
                
                % Include the function and line number where the error occurred
                if ~isempty(ME.stack)
                    errorReport = [errorReport, sprintf('File: %s\nFunction: %s\nLine: %d\n', ...
                        ME.stack(1).file, ME.stack(1).name, ME.stack(1).line)];
                end
                
                % Print the detailed report to the command window
                fprintf(1, '%s\n', errorReport); skippedWells{z} = pathFileActivityScan; continue;
        end
        msg = sprintf('[Run %d | Well %d] Finished successfully...\n', scan_runID, wellID);
        %fprintf(logFile, '%s', msg);   % Write to log file
        fprintf(1, '%s', msg); 
    end

    results{k} = vertcat(fileResults{:});
    skippedFiles{k} = vertcat(skippedWells{:});
end

% Consolidate all the structs from the nested cell arrays
allDataStructs = vertcat(results{:});

% --- Separate the two types of data ---
% Use vertcat to stack all the summary tables into one master table
allSummaries = vertcat(allDataStructs.summary);

% Do the same for all the nested distribution tables
allNestedDistributions = vertcat(allDataStructs.nested);

% --- Save the final files in their respective formats ---
% Save the summary data as a human-readable CSV file
writetable(allSummaries, fullfile(opDir, 'ActivityScan_outputs', 'Compiled_Activity_Summary.csv'));

% Save nested distributions to an HDF5 file that Python can read with h5py
h5file = fullfile(opDir, 'ActivityScan_outputs', 'Compiled_Distributions.h5');
if exist(h5file, 'file'), delete(h5file); end  % remove if exists

nRows = height(allNestedDistributions);
for colIdx = 1:width(allNestedDistributions)
    colName = allNestedDistributions.Properties.VariableNames{colIdx};
    colData = allNestedDistributions.(colName);

    % Handle cell array of vectors (e.g. ISI_s, Fano) separately
    if iscell(colData)
        % Get flattened sizes and max length
        maxLen = max(cellfun(@(x) numel(x), colData));
        dataMatrix = NaN(nRows, maxLen);

        for i = 1:nRows
            thisVec = colData{i};
            if ~isempty(thisVec)
                dataMatrix(i, 1:numel(thisVec)) = thisVec(:)';
            end
        end

        % Save as dataset: shape [nRows, maxLen]
        h5create(h5file, ['/' colName], size(dataMatrix), 'Datatype', 'double');
        h5write(h5file, ['/' colName], dataMatrix);

    else
        % Numeric column (e.g., Run_ID, DIV)
        h5create(h5file, ['/' colName], size(colData), 'Datatype', class(colData));
        h5write(h5file, ['/' colName], colData);
    end
end


% This part for skipped files remains the same
skippedFiles = vertcat(skippedFiles{~cellfun(@isempty, skippedFiles )});
if ~isempty(skippedFiles)
    fprintf('Skipped %d files due to read errors.\n', numel(skippedFiles));
end

fprintf('Activity Scan analysis successfully compiled.\n');
fprintf('Total time: %.2f seconds\n', toc);
end



function plotSummaryMontage(activityScanData, meanFiringRate, amplitude90perc, xpos, ypos, ...
    meanISI, runID, chipID, wellID, div, neuronType, ...
    outputDir, fixedFR_Caxis, fixedAmp_Caxis, fixedAA_Caxis, plotTitle)

% Create a summary montage for Firing Rate, Amplitude, Active Area, and ISI metrics

fig = figure('Visible','off','Color','w','Position',[100 100 1800 1200]);
t = tiledlayout(3,3,'TileSpacing','tight','Padding','compact');

% -- Add Plot Title --- 
title(t, plotTitle, 'Interpreter', 'none', 'FontSize', 12, 'FontWeight', 'bold')
% --- 1. Firing Rate Map ---
nexttile(t,1);
customCmap = customDivergingColorMap(256);
mxw.plot.activityMap(activityScanData, meanFiringRate,'ColorMap',customCmap, 'Ylabel', '[Hz]',...
    'CaxisLim', fixedFR_Caxis,'Interpolate',true,'Figure',false,'Title','Firing Rate Map');
title('Firing Rate Map');

% --- 2. Firing Rate Distribution ---
nexttile(t,4);
histogram(meanFiringRate(meanFiringRate>0.1), 0:.1:fixedFR_Caxis(2), 'FaceColor','#135ba3');
xlabel('Firing Rate [Hz]'); ylabel('Count'); title('Firing Rate Distribution');

% --- 3. Amplitude Map ---
nexttile(t,2);
mxw.plot.activityMap(activityScanData, amplitude90perc,'ColorMap',customCmap, 'Ylabel', '[uV]',...
    'CaxisLim', fixedAmp_Caxis,'Interpolate',true,'Figure',false,'Title','Spike Amplitude Map');
title('Amplitude Map');

% --- 4. Amplitude Distribution ---
nexttile(t,5);
histogram(amplitude90perc(amplitude90perc>10), 0:1:fixedAmp_Caxis(2), 'FaceColor','#135ba3');
xlabel('Amplitude [uV]'); ylabel('Count'); title('Amplitude Distribution');
% 
% % --- 5. Active Area Map (density) ---
% nexttile(t,3);
% xlim_chip = [0 3850]; ylim_chip = [0 2100]; BinSize = 10;
% xedges = xlim_chip(1):BinSize:xlim_chip(2);
% yedges = ylim_chip(1):BinSize:ylim_chip(2);
% N = histcounts2(ypos, xpos, yedges, xedges);
% N_smooth = imgaussfilt(N,2);
% imagesc(xedges(1:end-1), yedges(1:end-1), N_smooth);
% axis xy; axis equal; axis off;
% colormap('gray'); colorbar; caxis(fixedAA_Caxis);
% title('Active Area Density');


% --- 7. ISI Distribution ---
% Capture the handle for the specific axes you want to modify
ax_isi = nexttile(t,3); 

validISI = meanISI(~isnan(meanISI));

% Check if there is data to plot to avoid errors with empty histograms
if ~isempty(validISI)
    histogram(ax_isi, validISI, 0:0.01:2, 'FaceColor','#135ba3');
    
    % Use the captured handle 'ax_isi' instead of 'gca'
    set(ax_isi, 'Yscale', 'log'); 
else
    % If there's no data, you might want to turn off the axes ticks or label it
    set(ax_isi, 'Yscale', 'linear'); % Set to linear if log is not needed
end

% Also apply the handle to other commands for robustness
xlabel(ax_isi, 'ISI [s]'); ylabel(ax_isi, 'Count'); title(ax_isi, 'ISI Distribution');


% Save figure
figname = sprintf('%s/Montage_Run%s_Chip%s_Well%d_DIV%d_%s.svg', outputDir, runID, chipID, wellID, div, strrep(neuronType,' ','_'));
print(fig, figname, '-dsvg');close(fig);
end

