% This script goes through the h5 files in a parent folder,
% plot the Rasters, FiringRate Maps, and, Amplitude Maps
% and compile a csv with the mean Firing Rate and Amplitude over days


function [] = compileActivityFiles(data)

customAmpOutsideLoop =[];  % or customAmpOutsideLoop = [] write the value to fix the color maps if needed

tic;
mfilename('fullpath')
fileDir = [pwd,'/',mfilename];

plotFig =data.plotFig;


%ui components.
fig =data.fig;
d = uiprogressdlg(fig,'Title','Compiling files',...
    'Message','Start','Cancelable','on');
drawnow
logFile = data.logFile;
% Unpack the data structure

div0 = data.div0Date;
parentFolderPath = data.parentFolderPath;
refDir = data.refDir;
opDir = data.opDir;

    
% make output folder
outputFolders = {'AmplitudeMap', 'FiringRateMap','Active_Area','ElectrodeAmp'};


for i = 1:length(outputFolders)
    folderPath = fullfile(opDir, 'ActivityScan_outputs/', outputFolders{i});
    if ~isfolder(folderPath)
        mkdir(folderPath);
    end
end

% extract runID info from reference excel sheet
refTable = readtable(refDir);

% % Convert the 'Assay' column to lowercase and trim any leading/trailing whitespace
%assayColumn = strtrim(lower(refTable.Assay));
% 
% % Find rows that contain 'network today' or 'network'
% containsNetworkToday = contains(assayColumn, 'network today');
% containsNetwork = contains(assayColumn, 'network');
% 
% % Combine the conditions using logical OR
% combinedCondition = containsNetworkToday | containsNetwork;
% 
% % Extract unique run_ids based on the combined condition
% run_ids = unique(refTable.Run_(combinedCondition)); 
run_ids = unique(refTable.Run_);

% defines
% convert div 0 to date datatype
div0_date = datetime(div0, "InputFormat",'MM/dd/yyyy');





%% iterate through ActivityScans
% get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(parentFolderPath, '**/ActivityScan/**/*raw.h5'); 
theFiles = dir(filePattern);
numFiles = length(theFiles);
results = cell(numFiles, 1);
skippedFiles = cell(numFiles, 1); 

% Specify the exact number of cores to use
numCores = 6;  % Adjust this number based on your needs and resource availability

% Initialize or modify the existing parallel pool
currentPool = gcp('nocreate');  % Check for existing parallel pool

if isempty(currentPool)
    poolObj = parpool(numCores);  % Create a new pool if none exists
elseif currentPool.NumWorkers ~= numCores
    delete(currentPool);  % Close the existing pool if it does not match the desired size
    poolObj = parpool(numCores);  % Create a new pool with the correct number of workers
else
    poolObj = currentPool;  % Use the existing pool if it already matches the desired number of workers
end

if ~plotFig
    fprintf(1, 'skipping plotting\n');
end

for k = 1 : numFiles
    % Check for Cancel button press
    if d.CancelRequested
       d.close()
       error("User interruption")
    end
    % Update progress, report current estimate
    d.Value = k/length(theFiles);
    % reset recording info
    scan_runID = nan;
    scan_chipID = nan;
    scan_meanFiringRate = nan;
    scan_meanSpikeAmplitude = nan;
    hd5Date = nan; 
    scan_div = nan;
    
    pathFileActivityScan = fullfile(theFiles(k).folder);
    % extract dir information
    fileDirParts = strsplit(pathFileActivityScan, filesep); % split dir into elements
    scan_runID = str2double(fileDirParts{end-0}); % extract runID
    scan_runID_text = fileDirParts{end-0};
    scan_chipID = fileDirParts{end-2}; % extract chipID
    % folderDate = datetime(fileDirParts{end-3},'InputFormat','yyMMdd','Format','MM-dd-yyyy'); % extract date and convert to date object
    
    if ismember(scan_runID,run_ids)
        fprintf(1, 'Now reading %s\n', pathFileActivityScan);
        fprintf(logFile, 'Now reading %s\n', pathFileActivityScan);
        % create fileManager object for the Activity Scan
            
        idx = refTable.Run_ == scan_runID;

        % --- Wells ---
        if iscell(refTable.Wells_Recorded(idx))
            wells_str = refTable.Wells_Recorded{idx};
        else
            wells_str = refTable.Wells_Recorded(idx);
        end
        wellsIDs = str2double(strsplit(strtrim(wells_str), ','));
        
        % --- Neuron Types ---
        if iscell(refTable.NeuronSource(idx))
            neuron_str = refTable.NeuronSource{idx};
        else
            neuron_str = refTable.NeuronSource(idx);
        end
        neuronTypes = strtrim(strsplit(neuron_str, ','));
        
        % --- Enforce correspondence ---
        if length(neuronTypes) ~= length(wellsIDs)
            error('Mismatch: Each well must have a corresponding neuron type.');
        end
        assayColumn = refTable.Assay(idx);
        numWells = length(wellsIDs);
        fileResults = cell(numWells, 1);
        skippedWells = cell(numWells,1);
       
            for z = 1:numWells
            customAmp = customAmpOutsideLoop;
            wellID=wellsIDs(z);
            fprintf(1, 'Processing Well %d\n', wellID);
            neuronSourceType = neuronTypes(z);
            try

                activityScanData = mxw.fileManager(pathFileActivityScan,wellID);
            catch
                skippedWells{z} = [pathFileActivityScan,num2str(wellID)];
                continue
            end
    
            % get the startTime of the recordings
            hd5_time = activityScanData.fileObj.stopTime;
            try
                hd5Date = datetime(hd5_time,'InputFormat', 'yyyy-MM-dd HH:mm:ss');
            catch
                hd5Date = datetime(hd5_time,'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
            end

            try
            div0_datetime = datetime(div0_date, 'InputFormat', 'yourInputFormatHere');
            scan_div = floor(days(hd5Date - div0_datetime));
    
            % compute means

            % get the mean firing rate for each electrode
            meanFiringRate = mxw.activityMap.computeSpikeRate(activityScanData);
            % calculate the overall mean firing rate of the scan
            thrFiringRate = 0.1; % set a minimum spike rate threshold (Hz)
            scan_meanFiringRate = mean(meanFiringRate(meanFiringRate>thrFiringRate));
            % get th 90th percentile spike amplitude value for each electrode
            amplitude90perc = abs(mxw.activityMap.computeAmplitude90percentile(activityScanData));
            % calculate the overall mean spike amplitude of the scan
            thrAmp = 20;  % set a minimum spike amplitude threshold (uV)
            scan_meanSpikeAmplitude = mean(amplitude90perc(amplitude90perc>thrAmp));
            idx = (meanFiringRate>thrFiringRate & amplitude90perc>thrAmp);
            Active_area = sum(idx)/length(idx)*100;

            %xpos_elec = activityScanData.processedMap.xpos(idx);
            %ypos_elec= activityScanData.processedMap.ypos(idx);

            %compute ISI metrics
            nRecordings = length(activityScanData.extractedSpikes);
            nElectrodes = length(activityScanData.extractedSpikes(1).frameno);
            spikeTimesCell = cell(nElectrodes, 1);
            for r = 1:nRecordings
                for e = 1:nElectrodes
                    spikeTimesCell{e} = [spikeTimesCell{e}; activityScanData.extractedSpikes(r).frameno{e}./activityScanData.fileObj(1).samplingFreq];
                end
            end
            % Preallocate
            meanISI = NaN(nElectrodes, 1);
            stdISI  = NaN(nElectrodes, 1);
            isiCV   = NaN(nElectrodes, 1);
            fano    = NaN(nElectrodes, 1);
            
            % Compute ISI stats in vectorized fashion using cellfun
            isiCell = cellfun(@diff, spikeTimesCell, 'UniformOutput', false);
            validISI = cellfun(@(x) numel(x) >= 1 && isnumeric(x), isiCell);
            % Vectorized mean & std of ISIs
            % Compute mean/std
            
            meanISI(validISI) = cellfun(@(x) mean(double(x)), isiCell(validISI), 'UniformOutput', true);
            stdISI(validISI) = cellfun(@(x) std(double(x)), isiCell(validISI), 'UniformOutput', true);
            isiCV(validISI)   = stdISI(validISI) ./ meanISI(validISI);
            
            % --- Compute Fano factor ---
            validSpikes = cellfun(@(x) numel(x) >= 1, spikeTimesCell);
            
            % Estimate max duration for binning (in seconds)
            maxTime = max(cellfun(@(x) max(x), spikeTimesCell(validSpikes)));
            binEdges = 0:1:ceil(maxTime);  % 1-second bins
            
            for e = find(validSpikes)'
                st = spikeTimesCell{e};
                spikeCounts = histcounts(st, binEdges);
                if mean(spikeCounts) > 0
                    fano(e) = var(double(spikeCounts)) / mean(spikeCounts);
                else
                    fano(e) = NaN;
                end
            end

            % Save only summary metrics to main table
            fileResults{z} = table(scan_runID, assayColumn, scan_div, wellID, {neuronSourceType{1}}, hd5Date, {scan_chipID}, ...
                scan_meanFiringRate, scan_meanSpikeAmplitude, Active_area, ...
                mean(meanISI),mean(stdISI),mean(isiCV),mean(fano),...
                'VariableNames', {
                    'Run_ID', 'AssayType', 'DIV', 'Well', 'NeuronType', 'Time', 'Chip_ID', ...
                    'Mean_FiringRate', 'Mean_SpikeAmplitude', 'Active_area',...
                    'Mean_ISI','Std_ISI','Mean_ISI_CV','Mean_Fano'...
            });
            catch ME
                disp(ME.message)
                
                skippedWells{z} = [pathFileActivityScan,num2str(wellID)];
                
                continue
            end
            % plot amplitude and firing rate maps
            if plotFig

            try
                % get the mean firing rate for each electrode
                meanFiringRate = mxw.activityMap.computeSpikeRate(activityScanData);
                
                % define maximum firing rate used in plots as 99th percentile of firing
                % rate values
                maxFiringRate = mxw.util.percentile(meanFiringRate(meanFiringRate~=0),99);
                    
                % plot the list of mean firing rates as a map, given the electrode 
                % x and y coordinates in 'fileManagerObj.processedMap'
                customCmap = customDivergingColorMap(256); 
                figure('Color','w','visible','off');
                subplot(2,1,1);
                mxw.plot.activityMap(activityScanData, meanFiringRate,'ColorMap',customCmap, 'Ylabel', '[Hz]',...
                    'CaxisLim', [0.1 max(meanFiringRate)/5],'Interpolate',true,'Figure',false,'Title','Firing Rate Activity Map');
                % run the line above several times, experimenting with different 
                % [min max] range of the color gradient 'CaxisLim'
                
                % add scale bar 
                line([300 800],[2000+400 2000+400],'Color','k','LineWidth',4);
                axis off;
                text(340,2100+500,'0.5 mm','color','k');
                xlim([200 3750]);ylim([150 2500])
                
                % plot the distribution of firing rates across all of the electrodes
                subplot(2,1,2);
                thrFiringRate = 0.1; % set a minimum spike rate threshold (Hz)
                histogram(meanFiringRate(meanFiringRate>thrFiringRate),0:.1:ceil(maxFiringRate), ...
                    'FaceColor','#135ba3', "EdgeColor",'#414042')
                xlim([thrFiringRate maxFiringRate])
                ylabel('Counts');xlabel('Firing Rate [Hz]');
                box off;
                legend(['Mean Firing Rate = ',num2str(mean(meanFiringRate(meanFiringRate>thrFiringRate)),'%.2f'),...
                    ' Hz,  sd = ',num2str(std(meanFiringRate(meanFiringRate>thrFiringRate)),'%.2f')])
            
                print(gcf,append(opDir,'ActivityScan_outputs/FiringRateMap/FiringRateMap',scan_runID_text,'_WellID_',num2str(wellID),'_',num2str(scan_chipID),'_DIV',num2str(scan_div),'_',strrep(neuronSourceType{1},' ',''),'.svg'),'-dsvg');
           
                %exportgraphics(fig,fullfile(append(opDir,'ActivityScan_outputs/FiringRateMap/FiringRateMap',scan_runID_text,'_WellID_',num2str(wellID),'_',num2str(scan_chipID),'_DIV',num2str(scan_div),'_',strrep(neuronSourceType{1},' ',''),'.svg')),'ContentType','vector','BackgroundColor','none');
            catch ME
                fprintf('Unable to plot Firing Rate Map for run s%/n',scan_runID_text)
                fprintf('%s\n',ME.message)
                %fprintf(logFile,'Unable to plot Firing Rate Map for run s%/n',scan_runID_text);
            end
            
            % III. Spike Amplitude Activity Map and Distribution
            try
                % get th 90th percentile spike amplitude value for each electrode
                amplitude90perc = abs(mxw.activityMap.computeAmplitude90percentile(activityScanData));
                
                % define maximum amplitude used in plots as 99th percentile of amplitude values
                maxAmp = mxw.util.percentile(amplitude90perc(amplitude90perc~=0),99);
                
                % plot the mean firing rate vector as a map, given the electrode 
                % x and y coordinates in 'fileManagerObj.processedMap'
                figure('Color','w','visible','off');
                subplot(2,1,1);
                if isempty(customAmp)
                    customAmp = maxAmp;
                end
                mxw.plot.activityMap(activityScanData, amplitude90perc,'ColorMap',customCmap,'Ylabel', '[\muV]',...
                    'CaxisLim', [10 customAmp], 'Interpolate',true,'Figure',false,'Title','Spike Amplitude Activity Map');
                % run the line above several times, experimenting with different 
                % [min max] range of the color gradient 'CaxisLim'
                
                % add scale bar 
                line([300 800],[2000+400 2000+400],'Color','k','LineWidth',4);
                axis off;
                text(340,2100+500,'0.5 mm','color','k');
                xlim([200 3750])
                ylim([150 2500])
                
                % plot the distribution of spike amplitudes across all of the electrodes
                subplot(2,1,2);
                thrAmp = 10;  % set a minimum spike amplitude threshold (uV)
                histogram(amplitude90perc(amplitude90perc>thrAmp),ceil(0:1:maxAmp), ...
                    'FaceColor','#135ba3', "EdgeColor",'#414042')
                xlim([thrAmp maxAmp])
                ylabel('Counts');xlabel('Spike Amplitude [\muV]');
                box off;
                legend(['Mean Spike Amplitude = ',num2str(mean(amplitude90perc(amplitude90perc>thrAmp)),'%.2f'),...
                    ' \muV,  sd = ',num2str(std(amplitude90perc(amplitude90perc>thrAmp)),'%.2f')])
            
                print(gcf,append(opDir,'ActivityScan_outputs/AmplitudeMap/AmplitudeMap',scan_runID_text,'_WellID_',num2str(wellID),'_',num2str(scan_chipID),'_DIV',num2str(scan_div),'_',strrep(neuronSourceType{1},' ',''),'.svg'),'-dsvg');
                %exportgraphics(fig,append(opDir,'ActivityScan_outputs/AmplitudeMap/AmplitudeMap',scan_runID_text,'_WellID_',num2str(wellID),'_',num2str(scan_chipID),'_DIV',num2str(scan_div),'_',strrep(neuronSourceType{1},' ',''),'.svg'),"ContentType",'vector','BackgroundColor','none');
            catch ME
                fprintf('Unable to plot Amplitude Map for run %s\n', scan_runID_text)
                fprintf('%s\n',ME.message)
                %fprintf(logFile,'Unable to plot Amplitude Map for run %s/n', scan_runID_text);
            end
          try
            figure('Color', 'w', 'Visible', 'off');  % create invisible figure
        
            % Correct Maxwell Chip Physical Area
            xlim_chip = [0 3850];  % 3.85 mm = 3850 microns
            ylim_chip = [0 2100];  % 2.10 mm = 2100 microns
        
            % === [NEW] Use finer bin size ===
            BinSize = 10;  % Higher resolution than previous BinSize = 30
        
            % Define bin edges
            xedges = xlim_chip(1):BinSize:xlim_chip(2);
            yedges = ylim_chip(1):BinSize:ylim_chip(2);
        
            % Compute 2D histogram (density)
            N = histcounts2(ypos_elec, xpos_elec, yedges, xedges);  % Note (y,x) order!
        
            % === [NEW] Apply Gaussian smoothing to reduce graininess ===
           % Define a Gaussian kernel manually
            sigma = 2;
            kernelSize = 5;
            [xg, yg] = meshgrid(-kernelSize:kernelSize, -kernelSize:kernelSize);
            gKernel = exp(-(xg.^2 + yg.^2) / (2 * sigma^2));
            gKernel = gKernel / sum(gKernel(:));  % Normalize
            
            % Apply smoothing using convolution
            N_smooth = conv2(N, gKernel, 'same');
        
            % Plot the heatmap
            imagesc(xedges(1:end-1), yedges(1:end-1), N_smooth);  % Align with bin centers
            axis xy;
            axis equal;
            hold on;
        
            % Set colormap and colorbar
            colormap("gray");
            colorbar;
        
            % Set color limits dynamically or manually
            CaxisLim = [];
            if isempty(CaxisLim)
                caxis([0 max(N_smooth(:))]);
            else
                caxis(CaxisLim);
            end
        
            % Scalebar
            line([300 800], [1800 1800], 'Color', 'k', 'LineWidth', 4);  % Adjusted position
            text(340, 1850, '0.5 mm', 'color', 'k', 'FontSize', 8);
        
            % Set axis limits correctly
            xlim(xlim_chip);
            ylim(ylim_chip);
        
            axis off;
            box off;
        
            % Add title with active area percentage
            title(append('Active Area ', num2str(Active_area, '%.1f'), '%'), 'FontSize', 12);
        
            % Save figure
            activeHeatmapPath = append(opDir, 'ActivityScan_outputs/Active_Area/ActiveAreaDensity_', ...
                scan_runID_text, '_WellID_', num2str(wellID), '_', num2str(scan_chipID), ...
                '_DIV', num2str(scan_div), '_', strrep(neuronSourceType{1}, ' ', ''), '.svg');
            print(gcf, activeHeatmapPath, '-dsvg');
        
        catch ME
            fprintf('Unable to plot Active Area Density Heatmap for run %s\n', scan_runID_text)
            fprintf('%s\n', ME.message)
        end
            % % III. Electrodes 
            % try
            %     % get th 90th percentile spike amplitude value for each electrode
            %     amplitude90perc = abs(mxw.activityMap.computeAmplitude90percentile(activityScanData));
            % 
            %     % define maximum amplitude used in plots as 99th percentile of amplitude values
            %     maxAmp = mxw.util.percentile(amplitude90perc(amplitude90perc~=0),99);
            % 
            %     % plot the mean firing rate vector as a map, given the electrode 
            %     % x and y coordinates in 'fileManagerObj.processedMap'
            %     figure('Color','w','visible','off');
            %     subplot(2,1,1);
            %     if isempty(customAmp)
            %         customAmp = maxAmp;
            %     end
            %     mxw.plot.activityMap(activityScanData, amplitude90perc,'ColorMap',customCmap,'Ylabel', '[\muV]',...
            %         'CaxisLim', [10 customAmp], 'Interpolate',true,'Figure',false,'Title','Spike Amplitude Activity Map');
            %     % run the line above several times, experimenting with different 
            %     % [min max] range of the color gradient 'CaxisLim'
            %     selectNetworkElectrodes = mxw.util.electrodeSelection.networkRec(activityScanData,'ActivityThreshold',thrFiringRate,'AmplitudeThreshold',thrAmp);
            %     hold
            %     %scatter(activityScanData.processedMap.xpos(selectNetworkElectrodes.electrodes), activityScanData.processedMap.ypos(selectNetworkElectrodes.electrodes), 15, 'r','s', 'filled');
            %     scatter(selectNetworkElectrodes.xpos, selectNetworkElectrodes.ypos, 15, 'r','s', 'filled');
            % 
            %     % add scale bar 
            %     line([300 800],[2000+400 2000+400],'Color','k','LineWidth',4);
            %     axis off;
            %     text(340,2100+500,'0.5 mm','color','k');
            %     xlim([200 3750])
            %     ylim([150 2500])
            % 
            %     % plot the distribution of spike amplitudes across all of the electrodes
            %     subplot(2,1,2);
            %     thrAmp = 10;  % set a minimum spike amplitude threshold (uV)
            %     histogram(amplitude90perc(amplitude90perc>thrAmp),ceil(0:1:maxAmp), ...
            %         'FaceColor','#135ba3', "EdgeColor",'#414042')
            %     xlim([thrAmp maxAmp])
            %     ylabel('Counts');xlabel('Spike Amplitude [\muV]');
            %     box off;
            %     legend(['Mean Spike Amplitude = ',num2str(mean(amplitude90perc(amplitude90perc>thrAmp)),'%.2f'),...
            %         ' \muV,  sd = ',num2str(std(amplitude90perc(amplitude90perc>thrAmp)),'%.2f')])
            % 
            %     print(gcf,append(opDir,'ActivityScan_outputs/ElectrodeAmp/ElectrodeAmp',scan_runID_text,'_WellID_',num2str(wellID),'_',num2str(scan_chipID),'_DIV',num2str(scan_div),'_',strrep(neuronSourceType{1},' ',''),'.svg'),'-dsvg');
            %     %exportgraphics(fig,append(opDir,'ActivityScan_outputs/AmplitudeMap/AmplitudeMap',scan_runID_text,'_WellID_',num2str(wellID),'_',num2str(scan_chipID),'_DIV',num2str(scan_div),'_',strrep(neuronSourceType{1},' ',''),'.svg'),"ContentType",'vector','BackgroundColor','none');
            % catch ME
            %     fprintf('Unable to plot Amplitude Map for run %s\n', scan_runID_text)
            %     fprintf('%s\n',ME.message)
            %     %fprintf(logFile,'Unable to plot Amplitude Map for run %s/n', scan_runID_text);
            % end

            try

                matFileName = sprintf('Distributions_Run%s_Chip%s_Well%d_DIV%d.mat', ...
                    scan_runID_text, scan_chipID, wellID, scan_div);
                matFilePath = fullfile(opDir, 'ActivityScan_outputs', 'Distributions', matFileName);
                
                % Ensure folder exists
                if ~exist(fullfile(opDir, 'ActivityScan_outputs', 'Distributions'), 'dir')
                    mkdir(fullfile(opDir, 'ActivityScan_outputs', 'Distributions'));
                end
                
                % Save variables to MAT file
               % parsave(matFilePath, 'meanFiringRate', 'amplitude90perc', 'scan_runID', ...
               %      'scan_chipID', 'wellID', 'scan_div');
            catch ME
                fprintf('%s\n',ME.message)
            end 
          % active area
        end
        end
    else
     continue
    end



results{k} = vertcat(fileResults{~cellfun(@isempty, fileResults)});
skippedFiles{k} = vertcat(skippedWells{~cellfun(@isempty, skippedWells)});
end
    delete(gcp('nocreate'));
   % profile viewer;
    
    skippedFiles = vertcat(skippedFiles{~cellfun(@isempty, skippedFiles )});
    % Concatenate all tables from each primary file into the final table
    finalWriteTable = vertcat(results{~cellfun(@isempty, results)});
    
    writetable(finalWriteTable, fullfile(opDir,'ActivityScan_outputs/Compiled_ActivityScan.csv'));
    

    if ~isempty(skippedFiles)
        
        fprintf(1,'Unable to read file with runID: %s, file(s) skipped.\n',skippedFiles);
        %fprintf(1,'Unable to read file with runID: %s, file(s) skipped.\n',skippedFiles);
    end

fprintf('Activity Scan analysis successfully compiled.\n')
fprintf(' Total elapsed time for execution: %f seconds\n', toc);
fprintf(logFile,'Activity Scan analysis successfully compiled.\n');
end
