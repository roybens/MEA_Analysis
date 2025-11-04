function [] = compileNetworkFiles(data)
    % compileNetworkFiles  Process HD-MEA “Network” .h5 files in parallel,
    % detect bursts, plot rasters, and save compiled CSV and extended metrics.
    %
    % Usage:
    %   compileNetworkFiles(data)
    %
    % Inputs (fields of data struct):
    %   .fig                : UI figure handle (for progress dialog)
    %   .div0Date           : DIV-0 date string, e.g. '01/15/2025'
    %   .parentFolderPath   : root path containing **/Network/**/*raw.h5
    %   .refDir             : Excel/CSV “ref” file with Run_, Assay, Wells_Recorded, NeuronSource
    %   .opDir              : output base folder
    %   .gaussianSigma      : Gaussian kernel σ for firing rate
    %   .binSize            : bin size for firing rate
    %   .minPeakDistance    : min peak distance for burst detection
    %   .minProminence      : min peak prominence for burst detection
    %   .thresholdBurst     : threshold for burst detection
    %   .thresholdStartStop : fraction for burst start/stop edges
    %   .thresholdMethod    : threshold method name
    %   .xlim               : x-axis limit for plots
    %   .ylim               : y-axis limit for plots
    %   .plotFig            : true/false to generate raster plots
    %   .extMetricsFlag     : true/false to save extended metrics
    %   .logFile            : path to log file (optional)
    %
    % Example:
    %   data.parentFolderPath = '/data/HDMEA';
    %   data.refDir           = '/data/ref.xlsx';
    %   data.opDir            = '/data/output';
  



    %   compileNetworkFiles(data);

    required = {'fig','div0Date','parentFolderPath','refDir','opDir', ...    % --- Validate required fields ---
                'gaussianSigma','binSize','minPeakDistance','minPeakProminence', ...
                'thresholdBurst','thresholdStartStop','thresholdMethod', ...
                'xlim','ylim','plotFig','extMetricsFlag'};
    for f = required
        assert(isfield(data,f{1}), 'Missing data.%s', f{1});
    end

    
    fig                = data.fig;% --- Unpack data struct ---
    div0_date          = datetime(data.div0Date, "InputFormat",'MM/dd/yyyy');
    parentFolderPath   = data.parentFolderPath;
    refDir             = data.refDir;
    opDir              = data.opDir;
    gaussianSigma      = data.gaussianSigma;
    binSize            = data.binSize;
    minPeakDistance    = data.minPeakDistance;
    minPeakProminence      = data.minPeakProminence;          
    thresholdBurst     = data.thresholdBurst;
    thresholdStartStop = data.thresholdStartStop;      
    thresholdMethod    = data.thresholdMethod;
    xlimNetwork        = data.xlim;
    ylimNetwork        = data.ylim;
    plotFig            = data.plotFig;
    extMetricsFlag     = data.extMetricsFlag;
    logFile            = data.logFile;
   
    % --- Read reference table and identify runs of interest ---
    refTable = readtable(refDir);
    assayList = strtrim(lower(refTable.Assay));
    selNetwork = contains(assayList,'network') | contains(assayList,'network today');
    run_ids    = unique(refTable.Run_(selNetwork));

    % --- Discover .h5 files ---
    patterns = {fullfile(parentFolderPath, '**','Network','**','*raw.h5'),
                    fullfile(parentFolderPath,"data.raw.h5"),
                    fullfile(parentFolderPath),
                    };
    theFiles =[];
    for p = 1:length(patterns)
        candidateFiles = dir(patterns{p});
        if ~isempty(candidateFiles)
            theFiles = candidateFiles;
            break;
        end
    end
    if isempty(theFiles)
        error('No .h5 files found recheck your data path')
    end
    numFiles =numel(theFiles);

    % --- Prepare parallel pool ---
    numCores = 6;
    poolObj = gcp('nocreate');
    if isempty(poolObj) || poolObj.NumWorkers~=numCores
        if ~isempty(poolObj), delete(poolObj); end
        parpool(numCores);
    end

    networkResults = cell(numFiles,1);
    skippedFiles   = cell(numFiles,1);
    allExtMetrics  = cell(numFiles,1);
    d = uiprogressdlg(fig, 'Title', 'Compiling files', ...
        'Message', 'Start', 'Cancelable', 'on');

    info = struct( ...
            'div0_date', datetime(div0_date), ...
            'parentFolderPath', parentFolderPath, ...
            'refDir', refDir, ...
            'opDir', opDir, ...
            'gaussianSigma', gaussianSigma, ...
            'binSize', binSize, ...
            'minPeakDistance', minPeakDistance, ...
            'minPeakProminence', minPeakProminence, ...
            'thresholdBurst', thresholdBurst, ...
            'thresholdStartStop', thresholdStartStop, ...
            'thresholdMethod', thresholdMethod, ...
            'xlimNetwork', xlimNetwork, ...
            'ylimNetwork', ylimNetwork, ...
            'plotFig', plotFig, ...
            'extMetricsFlag', extMetricsFlag, ...
            'logFile', logFile ...
        );
        
        % --- Convert to JSON string ---
        jsonText = jsonencode(info, 'PrettyPrint', true);
        
        % --- Write to file ---
        jsonFile = fullfile(opDir, "network_param_data.json");
        fid = fopen(jsonFile, 'w');
        fprintf(fid, '%s', jsonText);
        fclose(fid);
    % --- Loop over files ---
    for k = 1:numFiles
        %if k > 1, continue;end    %for debugging
        if d.CancelRequested, error('User interruption'); end
        d.Value = k/numFiles;
        baseName = theFiles(k).name;
        fullPath = fullfile(theFiles(k).folder,baseName);

        % Extract runID and chipID from path
        parts = strsplit(theFiles(k).folder, filesep);
        scan_runID      = str2double(parts{end});
        scan_runID_text = parts{end};
        scan_chipID     = parts{end-2};

        % Skip runs not in ref table
        if ~ismember(scan_runID, run_ids), msg =sprintf('[Run %d not in Ref File...\n', scan_runID);fprintf(1, '%s', msg);continue;end
        idx         = refTable.Run_ == scan_runID;
        wellsIDs    = str2double( strsplit(string(refTable.Wells_Recorded(idx)),',') );
        neuronTypes = strtrim( strsplit(string(refTable.NeuronSource(idx)),',') );
        thisAssay   = refTable.Assay(idx);   

        if numel(wellsIDs)~=numel(neuronTypes)
            error('Mismatch between Wells_Recorded and NeuronSource counts for run %d.', scan_runID);
        end

        numWells    = numel(wellsIDs);
        fileResults = cell(numWells,1);
        skippedWells= cell(numWells,1);
        if extMetricsFlag, extData = cell(numWells,1); end

        % successFile = false; attempt = 0; maxRetries = 3;
        % while ~successFile && attempt<maxRetries
        %     attempt = attempt + 1;

            try
                durations    = {'Plot60s','Plot120s','Plot300s'};
                formats      = {'svg','png'};
                wellStatus = strings(numWells,1);
                parfor z = 1:numWells
                  try
                    %if z>1,continue;end      %for testing purpose
                    wellID = wellsIDs(z);
                    neuronSourceType = neuronTypes{z};

                    % Prepare output folders for rasters/plots %not every
                    % iteration it is required
                    folderName   = sprintf('%s_Well%d', scan_chipID, wellID);
                    chipWellDir  = fullfile(opDir,'Network_outputs','Raster_BurstActivity',folderName);
                       for d_index = 1:numel(durations)
                            for f = 1:numel(formats)
                                outDir = fullfile(chipWellDir, durations{d_index}, formats{f});
                                if ~isfolder(outDir)
                                    mkdir(outDir);  % safe because this is outside parfor
                                end
                            end
                       end
                    msg = sprintf('[Run %d | Well %d] Processing...\n', scan_runID, wellID);
                        %fprintf(logFile, '%s', msg);   % Write to log file
                    fprintf(1, '%s', msg);  
                    % Load networkData
                    networkData = mxw.fileManager(fullPath, wellID);

                    % Compute DIV relative to div0_date
                    try
                        hd5Date = datetime(networkData.fileObj.stopTime,'InputFormat','yyyy-MM-dd HH:mm:ss');
                    catch
                        hd5Date = datetime(networkData.fileObj.stopTime,'InputFormat','dd-MMM-yyyy HH:mm:ss');
                    end
                    scan_div = floor(days(hd5Date - div0_date));

                    % Compute firing rate & basic network stats
                    relativeSpikeTimes = computeRelativeSpikeTimes(networkData);
                    networkAct = mxw.networkActivity.computeNetworkAct( networkData,'BinSize',binSize,'GaussianSigma',gaussianSigma);
                    networkStats = computeNetworkStatsModified( networkAct.firingRate, networkAct.time, 'ThresholdMethod', thresholdMethod, ...
                        'Threshold', thresholdBurst, 'MinPeakProminence', minPeakProminence,  'MinPeakDistance', minPeakDistance);

                    AbsBP = networkAct.absfiringRate( ismember(networkAct.time, networkStats.maxAmplitudesTimes) );

                    % --- Burst detection & extended metrics ---
                    [burstsSummaryRow, burstMetricsVariables] = gaussianFiringRateBurstDetector(networkData,networkAct, networkStats, AbsBP, thresholdStartStop);

                    baselineFiringRate = burstsSummaryRow.BaselineFiringRate;
                    % --- 2) Original log‐ISI burst detection ---
                    % For each channel:
                   
                    unique_channels = unique(relativeSpikeTimes.channel);
                    chanBursts = cell(numel(unique_channels), 1);
                    
                    for i = 1:numel(unique_channels)
                        ch = unique_channels(i);
                        st = relativeSpikeTimes.time(relativeSpikeTimes.channel == ch);
                        chanBursts{i} = detectBursts_logISI(st, 10, 0.100, []); %Chiappalone et al., 2005 10 spikes 
                    end

                    % --- Step 1: Collect all bursts from each electrode ---
                    allBurstEvents = [];  % Will store: [start_time, end_time, channel_index]
                    
                    for i = 1:length(chanBursts)
                        bursts = chanBursts{i};
                        for b = 1:numel(bursts)
                            % Each row = [start_time, end_time, channel_index]
                            allBurstEvents(end+1, :) = [bursts(b).start_time, bursts(b).end_time, i]; %#ok<AGROW>
                        end
                    end
                    
                   
              
                    % Build metadata + results table
                    metaTbl = table( scan_runID, scan_div, thisAssay, wellID, string(neuronSourceType), ...
                                     hd5Date, string(scan_chipID),'VariableNames',{'Run_ID','DIV','Assay','Well','NeuronType','Time','Chip_ID'} );
                    
                    fileResults{z} = [metaTbl, burstsSummaryRow];%,logSummary];

                    % Save extended metrics if desired
                    if extMetricsFlag
                        ts = ((double(networkData.fileObj.spikes.frameno)...
                            - double(networkData.fileObj.firstFrameNum))/networkData.fileObj.samplingFreq)';
                        ts = ts(ts>0);
                        extData{z} = [metaTbl,computeExtendedBurstMetrics( burstMetricsVariables, relativeSpikeTimes, ts)];
                    end
                    
                    % --- Plot rasters & network traces ---
                    if plotFig
                        textStr = sprintf('#Bursts:%d | MeanDur:%.1fs | SpB:%.1f | IBI:%.1fs | Peak:%.1f', ...
                                          burstsSummaryRow.Number_Bursts, burstsSummaryRow.mean_BurstDuration,  burstsSummaryRow.mean_Spike_per_Burst, ...
                                          burstsSummaryRow.mean_IBI,burstsSummaryRow.mean_Burst_Peak);
                        
                        plotFileBase = sprintf('Raster_%s_Well%d_%s_DIV%d_%s',scan_runID_text, wellID, scan_chipID, scan_div, neuronSourceType);
                         % CREATE TITLE STRING:
                        plotTitle = sprintf('Raster_%s_Well%d_%s_DIV%d_%s', scan_runID_text, wellID, ...
                        scan_chipID, scan_div, strrep(neuronSourceType, ' ', ''));
                        locData = struct()
                        locData.channel = networkData.fileObj.map.channel;
                        locData.xpos = networkData.fileObj.map.x;
                        locData.ypos = networkData.fileObj.map.y;
                        plotRasterNetwork( networkAct, networkStats, ...
                                           relativeSpikeTimes,locData, binSize, opDir, chipWellDir, ...
                                           xlimNetwork, ylimNetwork, textStr, plotFileBase, ...
                                            baselineFiringRate);
     
                    end

                wellStatus(z) ='SUCCESS';
                catch ME
                % Restart pool on parallel errors
                if contains(ME.message,'parallel')
                    delete(gcp('nocreate'));
                    parpool(numCores);
                else

                       wellStatus(z) = "FAILED";
                       skippedWells{z} = sprintf('[Run %d | Well %d] ERROR: %s\n', scan_runID, wellID, ME.message);

                       % Log specific error per well
                       fprintf(2, '[Run %d | Well %d] FAILED: %s\n', scan_runID, wellID, ME.message);
                    if ~isempty(ME.stack)
                        fprintf(2, 'In function: %s\n', ME.stack(1).name);
                        fprintf(2, 'At line: %d\n', ME.stack(1).line);
                        fprintf(2, 'Full Call Stack:\n');
                        for si = 1:length(ME.stack)
                            fprintf(2, '   [%02d] %s (line %d)\n', si, ME.stack(si).name, ME.stack(si).line);
                        end
                    end

                    continue;
                end
            end
                end  % parfor

                % Collect results
                % only keep the non‐empty table cells
                isTbl = cellfun(@(c) istable(c) && ~isempty(c), fileResults);
                if any(isTbl)
                    networkResults{k} = vertcat( fileResults{isTbl} );
                else
                    networkResults{k} = table();   % empty table if nothing found
                end
                skippedFiles{k}   = vertcat(skippedWells{~cellfun(@isempty,skippedWells)});
                if extMetricsFlag
                    allExtMetrics{k} = vertcat(extData{:});
                end

                successFile = true;
                catch ME
                    % Structured error logging
                    fprintf(2, '\n============================================\n');
                    fprintf(2, 'Error in file: %s\n', fullPath);
                    fprintf(2, 'Time: %s\n', datestr(now));
                    fprintf(2, 'Error Message: %s\n', ME.message);
                
                    % Show where the error originated
                    if ~isempty(ME.stack)
                        fprintf(2, 'In function: %s\n', ME.stack(1).name);
                        fprintf(2, 'At line: %d\n', ME.stack(1).line);
                        fprintf(2, 'Full Call Stack:\n');
                        for si = 1:length(ME.stack)
                            fprintf(2, '   [%02d] %s (line %d)\n', si, ME.stack(si).name, ME.stack(si).line);
                        end
                    end
                
                    % Code excerpt around the error line
                    try
                        fileText = fileread(which(ME.stack(1).file));
                        fileLines = strsplit(fileText, '\n');
                        errLine = ME.stack(1).line;
                        range = max(1, errLine-2):min(length(fileLines), errLine+2);
                        fprintf(2, '\nCode excerpt:\n');
                        for li = range
                            mark = '';
                            if li == errLine, mark = '>> '; else, mark = '   '; end
                            fprintf(2, '%s%4d: %s\n', mark, li, fileLines{li});
                        end
                    catch
                        fprintf(2, 'Could not load code excerpt.\n');
                    end
                
                    fprintf(2, '============================================\n\n');
                
                    % Handle parallel pool recovery if needed
                    if contains(ME.message, 'parallel')
                        delete(gcp('nocreate'));
                        parpool(numCores);
                    else
                        continue;
                    end
        end
      end

    %     if ~successFile
    %         warning('Failed processing %s after %d attempts', fullPath, maxRetries);
    %     end
    % end  % for files

    delete(gcp('nocreate'));
    delete(d);

    % --- Write compiled CSV ---
    tblOK    = cellfun(@(c) istable(c) && ~isempty(c), networkResults);
    allResults = vertcat( networkResults{tblOK} );

    csvPath = fullfile(opDir,'Network_outputs','Compiled_Networks.csv');
    if isfile(csvPath), delete(csvPath); end
    writetable(allResults, csvPath);

    % --- Save extended metrics ---
    if extMetricsFlag
        save(fullfile(opDir,'Network_outputs','extendedMetrics.mat'),'allExtMetrics');
    end

    fprintf('Compilation complete. %d files processed in %.1fs\n', ...
            numel(networkResults), toc);


end



    % function out = plotFRMatrices(data)
    %             frbinSizes = [0.01,0.1,1,10];
    %             if plotFRmatrix
    %             for frbinSize = frbinSizes
    % 
    %                     % Ensure output directories exist for this bin size
    %                     % Define directory paths
    %                     firingRateMatrixFolder = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', sprintf('FiringRateMatrix_BinSize_%2f', frbinSize));
    %                     heatmapFolder = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', sprintf('Heatmaps_BinSize_%2f', frbinSize));
    %                     svdFolder = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', sprintf('SVD_BinSize_%2f', frbinSize));
    %                     logheatmapFolder = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', sprintf('Log_Heatmaps_BinSize_%2f', frbinSize));
    %                     % Create directories if they don't exist
    %                     if ~isfolder(firingRateMatrixFolder)
    %                         mkdir(firingRateMatrixFolder);
    %                     end
    %                     if ~isfolder(logheatmapFolder)
    %                         mkdir(logheatmapFolder);
    %                     end
    %                     if ~isfolder(heatmapFolder)
    %                         mkdir(heatmapFolder);
    %                     end
    %                     if ~isfolder(svdFolder)
    %                         mkdir(svdFolder);
    %                     end
    %                     % Compute firing rate matrix
    %                     binEdges = 0:frbinSize:networkData.fileObj.dataLenTime;
    %                     firingRateMatrix = zeros(max(ch), length(binEdges) - 1);
    %                     for chIdx = 1:max(ch)
    %                         spikeTimes = ts(ch == chIdx);
    %                         firingRateMatrix(chIdx, :) = histcounts(spikeTimes, binEdges);
    %                     end
    %                                             %%SVD computation
    %                     [U, S, V] = svd(firingRateMatrix, 'econ');
    % 
    %                     % Extract the singular values from the diagonal matrix S
    %                     singularValues = diag(S);
    % 
    %                     % Number of singular vectors (columns) to visualize
    %                     %numSingularValues = length(singularValues);
    %                     numSingularValues = 20;
    % 
    %                     % Create a single figure for combined visualization
    %                     figure('Visible', 'off');
    % 
    %                     % Subplot 1: Singular values as a line plot
    %                     subplot(3, 1, 1);
    %                     plot(singularValues, '-o', 'LineWidth', 1.5);
    %                     xlabel('Component number');
    %                     ylabel('Singular value (variance explained)');
    %                     title(sprintf('Singular values from SVD (Bin Size: %f)', frbinSize));
    %                     grid on;
    % 
    %                     % Subplot 2: Heatmap of all left singular vectors (U)
    %                     subplot(3, 1, 2);
    %                     imagesc(U(:, 1:numSingularValues));  % Use all columns of U
    %                     colorbar;
    %                     xlabel('Singular vector index');
    %                     ylabel('Channel');
    %                     title('All left singular vectors (U)');
    % 
    %                     % Subplot 3: Heatmap of all right singular vectors (V)
    %                     subplot(3, 1, 3);
    %                     imagesc(V(:, 1:numSingularValues));  % Use all columns of V
    %                     colorbar;
    %                     xlabel('Singular vector index');
    %                     ylabel('Time bin');
    %                     title('All right singular vectors (V)');
    % 
    %                     % Save the combined figure
    %                     combinedPlotPath = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity',...
    %                         sprintf('SVD_BinSize_%2f', frbinSize),...
    %                         sprintf('CombinedSVD_RunID%d_%s_WellID%d_DIV%d.png', scan_runID,scan_chipID ,wellID, scan_div));
    %                     saveas(gcf, combinedPlotPath);
    %                     close(gcf);
    %                     % % Save firing rate matrix
    %                     % firingRateMatrixPath = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', ...
    %                     %     sprintf('FiringRateMatrix_BinSize_%.2f', frbinSize), ...
    %                     %     sprintf('FiringRateMatrix_RunID%d_WellID%d.mat', scan_runID, wellID));
    %                     % save(firingRateMatrixPath, 'firingRateMatrix');
    % 
    %                     % Generate and save heatmaps for 60s, 120s, and 300s
    %                     timeIntervals = [60, 120, 300];
    %                     for tIdx = 1:length(timeIntervals)
    %                         tLimit = timeIntervals(tIdx);
    %                         tBinIdx = find(binEdges <= tLimit, 1, 'last')-1;
    %                         figure('Visible', 'off');
    %                         imagesc(firingRateMatrix(:, 1:tBinIdx));
    %                         title(sprintf('Firing Rate Heatmap (First %ds, Bin Size %2f)', tLimit, frbinSize));
    %                         xlabel('Time Bins');
    %                         ylabel('Channel');
    %                         colorbar;
    %                         colormap('hot');
    %                         heatmapPath = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', ...
    %                             sprintf('Heatmaps_BinSize_%2f', frbinSize), ...
    %                             sprintf('Heatmap_RunID%d_%s_WellID%d_DIV%d_%ds.png', scan_runID,scan_chipID ,wellID, scan_div,tLimit));
    %                         saveas(gcf, heatmapPath);
    %                         close(gcf);
    % 
    %                         % Log-scaled heatmap
    %                         figure('Visible', 'off');
    %                         % Add 1 to avoid log(0) issues
    %                         imagesc(log10(firingRateMatrix(:, 1:tBinIdx) + 1));
    %                         title(sprintf('Log-Scaled Firing Rate Heatmap (First %ds, Bin Size %2f)', tLimit, frbinSize));
    %                         xlabel('Time Bins');
    %                         ylabel('Channel');
    %                         colorbar;
    %                         colormap('hot');
    %                         logHeatmapPath = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', ...
    %                             sprintf('Log_Heatmaps_BinSize_%2f', frbinSize), ...
    %                             sprintf('Log_Heatmap_RunID%d_%s_WellID%d_DIV%d_%ds.png', scan_runID,scan_chipID, wellID, scan_div, tLimit));
    %                         saveas(gcf, logHeatmapPath);
    %                         close(gcf);
    %                     end
    % 
    %               end
    %             end
    % end
    % 
    % 
    % 
   



