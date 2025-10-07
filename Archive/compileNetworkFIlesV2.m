function [] = compileNetworkFIlesV2(data)
% Main function to compile MEA network recording files, calculate metrics,
% and generate plots. This version is modularized for better readability
% and maintainability.

    % --- 1. UI Components and Initialization ---
    fig = data.fig;
    fprintf(1, 'Inside compileNetworkFiles (Modularized)\n');
    d = uiprogressdlg(fig, 'Title', 'Compiling files', ...
        'Message', 'Start', 'Cancelable', 'on');
    drawnow;
    tic;

    % --- 2. Unpack Configuration and Parameters ---
    plotFig = data.plotFig;
    extMetricsFlag = data.extMetricsFlag;
    logFile = data.logFile;

    % Analysis parameters
    analysisParams = struct();
    analysisParams.div0_date = datetime(data.div0Date, "InputFormat", 'MM/dd/yyyy');
    analysisParams.gaussianSigma = data.gaussianSigma;
    analysisParams.binSize = data.binSize;
    analysisParams.minPeakDistance = data.minPeakDistance;
    analysisParams.minProminence = data.minProminience; % Corrected typo from minProminience
    analysisParams.thresholdBurst = data.thresholdBurst;
    analysisParams.thresholdStartStop = data.thresholdStartStop;
    analysisParams.thresholdFunction = data.thresholdMethod;
    analysisParams.frBinSizesSVD = [0.01, 0.1, 1, 10]; % For SVD/Heatmap
    analysisParams.plotFRmatrix = false; % Set this to true if SVD/heatmaps are needed, pass from data struct if dynamic

    % LogISI specific parameters (can also be moved into data struct)
    analysisParams.logISI_min_spikes = 5;
    analysisParams.logISI_min_duration_s = 0.05; % 50 ms
    analysisParams.logISI_tuning_params = struct('log_isi_bins', 50, ...
                                               'smooth_span_hist', 5, ...
                                               'max_intra_burst_isi_s', 0.15); % 150ms

    % Plotting parameters
    plottingParams = struct();
    plottingParams.opDir = data.opDir;
    plottingParams.xlimNetwork = data.xlim;
    plottingParams.ylimNetwork = data.ylim;
    plottingParams.epsPlot = false; % Set true to save EPS plots
    plottingParams.logFile = logFile; % Pass logFile for error logging in plots if needed

    % Ensure base output directories for aggregated results exist
    if ~isfolder(fullfile(plottingParams.opDir, 'Network_outputs'))
        mkdir(fullfile(plottingParams.opDir, 'Network_outputs'));
    end
    if plotFig && ~isfolder(fullfile(plottingParams.opDir, 'Network_outputs', 'Raster_BurstActivity'))
         mkdir(fullfile(plottingParams.opDir, 'Network_outputs', 'Raster_BurstActivity'))
    end
    if plotFig && analysisParams.plotFRmatrix && ~isfolder(fullfile(plottingParams.opDir, 'Network_outputs', 'FiringRateData'))
        mkdir(fullfile(plottingParams.opDir, 'Network_outputs', 'FiringRateData'))
    end


    % --- 3. Load Reference Data and Discover Files ---
    [run_ids_to_process, refTable] = loadAndFilterReferenceData(data.refDir, {'network today', 'network'});
    theFiles = findMEADataFiles(data.parentFolderPath, '**/Network/**/*raw.h5');
    
    numFiles = length(theFiles);
    if numFiles == 0
        d.Message = 'No raw.h5 files found in specified path.';
        pause(2); close(d);
        fprintf(1,'No raw.h5 files found. Exiting.\n');
        return;
    end

    networkResultsAggregated = cell(numFiles, 1);
    allExtMetricsAggregated = cell(numFiles, 1);
    skippedFilesLog = {};

    % --- 4. Parallel Pool Setup ---
    numCores = 6; % Or data.numCores if passed
    currentPool = gcp('nocreate');
    if isempty(currentPool) || currentPool.NumWorkers ~= numCores
        if ~isempty(currentPool), delete(currentPool); end
        try
            parpool(numCores);
        catch parE
            fprintf(logFile, 'Warning: Could not start parallel pool with %d cores. Using default. Error: %s\n', numCores, parE.message);
            parpool(); % Try default
        end
    end
    
    fprintf(1, 'Plotting Rasters: %d, Extended Metrics: %d\n', plotFig, extMetricsFlag);
    fprintf(logFile, 'Plotting Rasters: %d, Extended Metrics: %d\n', plotFig, extMetricsFlag);

    % --- 5. Main File Processing Loop ---
    for k = 1:numFiles
        if d.CancelRequested
            fprintf(logFile, 'User cancelled processing.\n');
            break;
        end
        d.Value = k / numFiles;
        d.Message = sprintf('Processing file %d/%d: %s', k, numFiles, theFiles(k).name);
        drawnow;

        currentFilePath = fullfile(theFiles(k).folder, theFiles(k).name);
        [file_scan_runID, file_scan_chipID, file_wellsToProcess, file_neuronTypes, file_assayType] = ...
            getFileMetadata(currentFilePath, refTable, run_ids_to_process);

        if isempty(file_wellsToProcess)
            fprintf(1, 'Skipping file (no relevant RunID or wells): %s\n', currentFilePath);
            fprintf(logFile, 'Skipping file (no relevant RunID or wells): %s\n', currentFilePath);
            skippedFilesLog{end+1} = currentFilePath; %#ok<AGROW>
            continue;
        end
        
        numWellsInFile = length(file_wellsToProcess);
        wellResultsForFile = cell(numWellsInFile, 1);
        extMetricsForFile = cell(numWellsInFile, 1);
        
        % Retry logic for each file (especially for parfor issues)
        successFileProcess = false; attempt = 0; maxRetries = 3;
        while ~successFileProcess && attempt < maxRetries
            attempt = attempt + 1;
            try
                parfor z = 1:numWellsInFile
                    well_ID = file_wellsToProcess(z);
                    neuron_Type = file_neuronTypes{z}; % Ensure this indexing is robust
                    assay_Type = file_assayType{1}; % Assuming assay type is same for all wells in a run
                    
                    fprintf(1, 'File %d/%d, Well %d: Starting...\n', k, numFiles, well_ID);

                    % Create dynamic output folders for this well's plots
                    if plotFig
                        createWellPlotFolders(plottingParams.opDir, file_scan_chipID, well_ID, analysisParams.plotFRmatrix, analysisParams.frBinSizesSVD);
                    end
                    
                    % --- Core Well Processing ---
                    [wellData, relSpikeTimes, ts, ch, recordingTime_s, scan_div] = ...
                        loadWellData(currentFilePath, well_ID, analysisParams.div0_date);

                    if isempty(ts) || recordingTime_s <= 0
                        fprintf(1, 'Well %d: No spike data or zero recording time. Skipping metrics.\n', well_ID);
                        wellResultsForFile{z} = createEmptyMetricsRow(file_scan_runID, scan_div, assay_Type, well_ID, neuron_Type, now, file_scan_chipID, "No spike data");
                        continue;
                    end

                    networkActStruct = calculateSmoothedNetworkActivity(wellData, analysisParams.binSize, analysisParams.gaussianSigma);
                    
                    networkBurstPeakStats = detectNetworkBursts_PeakProminence(networkActStruct.firingRate, networkActStruct.time, analysisParams);
                    
                    [burstSummaryMetrics, derivedBurstChars] = calculateDerivedNetworkBurstStats(networkActStruct, networkBurstPeakStats, ts, ch, analysisParams.thresholdStartStop, recordingTime_s);
                    
                    % Combine basic metadata with calculated burst metrics
                    wellRow = table(file_scan_runID, scan_div, {assay_Type}, well_ID, {strtrim(neuron_Type)}, wellData.fileObj.hd5Date, {file_scan_chipID}, ...
                        'VariableNames', {'Run_ID', 'DIV', 'Assay', 'Well', 'NeuronType', 'Time', 'Chip_ID'});
                    wellRow = [wellRow, burstSummaryMetrics]; %#ok<AGROW>
                    
                    if extMetricsFlag
                        extMetricsData = calculateExtendedWellMetrics(ts, ch, derivedBurstChars.burstEdges, recordingTime_s);
                        % --- Add Per-Channel LogISI Burst Analysis ---
                        extMetricsData.perChannelLogISI = analyzePerChannelBursts_logISI(ts, ch, analysisParams.logISI_min_spikes, analysisParams.logISI_min_duration_s, analysisParams.logISI_tuning_params);
                        extMetricsForFile{z} = extMetricsData;
                    end
                    
                    if plotFig
                        plottingParamsCurrentWell = plottingParams; % Pass a copy
                        plottingParamsCurrentWell.scan_chipID = file_scan_chipID;
                        plottingParamsCurrentWell.wellID = well_ID;
                        plottingParamsCurrentWell.scan_div = scan_div;
                        plottingParamsCurrentWell.neuronSourceType = neuron_Type;
                        plottingParamsCurrentWell.scan_runID_text = sprintf('%d',file_scan_runID); % For filenames

                        % Add derivedBurstChars to plottingParams if plotNetworkActivityModified needs them
                        plottingParamsCurrentWell.nBursts = burstSummaryMetrics.Number_Bursts;
                        plottingParamsCurrentWell.meanBurstDuration = burstSummaryMetrics.mean_BurstDuration;
                        plottingParamsCurrentWell.meanSpikesPerBurst = burstSummaryMetrics.mean_Spike_per_Burst; % Check if this is directly used in plot annotation
                        plottingParamsCurrentWell.meanIBI = burstSummaryMetrics.mean_IBI;
                        plottingParamsCurrentWell.meanBurstPeak = burstSummaryMetrics.mean_Burst_Peak;
                        plottingParamsCurrentWell.baselineFiringRate = burstSummaryMetrics.Baseline;


                        plotAllWellFigures(wellData, relSpikeTimes, networkActStruct, networkBurstPeakStats, plottingParamsCurrentWell, analysisParams, derivedBurstChars.baselineFiringRate);
                        
                        if analysisParams.plotFRmatrix
                             plotFiringRateSVD(ts, ch, recordingTime_s, analysisParams.frBinSizesSVD, plottingParamsCurrentWell);
                        end
                    end
                    wellResultsForFile{z} = wellRow;
                    fprintf(1, 'File %d/%d, Well %d: Completed.\n', k, numFiles, well_ID);
                end % end parfor z
                successFileProcess = true; % If parfor completes without throwing an error to this level
            catch ME_parfor
                fprintf(logFile, 'ERROR during parfor processing of file %s (attempt %d/%d): %s\n', currentFilePath, attempt, maxRetries, ME_parfor.message);
                fprintf(1, 'ERROR during parfor processing of file %s (attempt %d/%d): %s. Check log.\n', currentFilePath, attempt, maxRetries, ME_parfor.message);
                for i_stack = 1:length(ME_parfor.stack)
                    fprintf(logFile, '  at %s (line %d)\n', ME_parfor.stack(i_stack).name, ME_parfor.stack(i_stack).line);
                end

                if contains(ME_parfor.message, 'parallel') && attempt < maxRetries
                    fprintf(logFile,'Restarting parallel pool due to error...\n');
                    delete(gcp('nocreate')); parpool(numCores); 
                elseif attempt == maxRetries
                    fprintf(logFile,'Max retries reached for file %s. Skipping this file.\n', currentFilePath);
                    skippedFilesLog{end+1} = [currentFilePath, ' (parfor error)']; %#ok<AGROW>
                end
            end % end try-catch for parfor
        end % end while retry loop for file

        if successFileProcess && ~isempty(wellResultsForFile)
            networkResultsAggregated{k} = vertcat(wellResultsForFile{~cellfun(@isempty, wellResultsForFile)});
            if extMetricsFlag && ~isempty(extMetricsForFile)
                allExtMetricsAggregated{k} = [extMetricsForFile{~cellfun(@isempty, extMetricsForFile)}];
            end
        elseif ~successFileProcess
             skippedFilesLog{end+1} = [currentFilePath, ' (failed all retries)']; %#ok<AGROW>
        end

    end % end for k (files)

    % --- 6. Finalize and Save Results ---
    if ~isempty(gcp('nocreate')), delete(gcp('nocreate')); end % Close parallel pool

    finalTable = vertcat(networkResultsAggregated{~cellfun(@isempty, networkResultsAggregated)});
    if extMetricsFlag
        finalExtMetrics = [allExtMetricsAggregated{~cellfun(@isempty, allExtMetricsAggregated)}];
         % The above concatenation might need adjustment depending on how extMetricsForFile are structured.
         % If extMetricsForFile{z} is a struct, then finalExtMetrics will be an array of structs.
    end

    saveCompiledResults(finalTable, ...
        fullfile(plottingParams.opDir, 'Network_outputs/Compiled_Networks.csv'), ...
        extMetricsFlag, ...
        finalExtMetrics, ... % Corrected variable name
        fullfile(plottingParams.opDir, 'Network_outputs/extendedMetrics.mat'));

    if ~isempty(skippedFilesLog)
        fprintf(1, 'Warning: Some files/wells were skipped. See log file and below:\n');
        fprintf(logFile, '--- Skipped Files/Wells Summary ---\n');
        for i=1:length(skippedFilesLog)
            fprintf(1, '  %s\n', skippedFilesLog{i});
            fprintf(logFile, '  %s\n', skippedFilesLog{i});
        end
    end
    
    close(d);
    fprintf(logFile, 'Network analysis successfully compiled.\n');
    fprintf(1, 'Total elapsed time for execution: %f seconds\n', toc);
    fprintf(1, 'Network analysis successfully compiled.\n');

end

% =========================================================================
% LOCAL HELPER FUNCTIONS (Modularized Parts)
% These can also be moved to separate .m files if preferred
% =========================================================================

function [run_ids, refTable] = loadAndFilterReferenceData(refDir, assayKeywords)
    % Loads reference table and filters run_ids based on assay keywords.
    refTable = readtable(refDir);
    assayColumn = strtrim(lower(refTable.Assay));
    combinedCondition = false(size(assayColumn));
    for i = 1:length(assayKeywords)
        combinedCondition = combinedCondition | contains(assayColumn, assayKeywords{i});
    end
    run_ids = unique(refTable.Run_(combinedCondition));
end

function [theFiles] = findMEADataFiles(parentFolderPath, filePattern)
    % Finds MEA data files matching the pattern.
    theFiles = dir(fullfile(parentFolderPath, filePattern));
end

function [scan_runID, scan_chipID, wellsToProcess, neuronTypes, assayType] = getFileMetadata(filePath, refTable, valid_run_ids)
    % Extracts metadata for a given file from the reference table.
    fileDirParts = strsplit(filePath, filesep);
    scan_runID_str = fileDirParts{end-1};
    scan_runID = str2double(scan_runID_str);
    scan_chipID = fileDirParts{end-3};
    
    wellsToProcess = [];
    neuronTypes = {};
    assayType = {};

    if ismember(scan_runID, valid_run_ids)
        idx_run = refTable.Run_ == scan_runID;
        if ~any(idx_run)
             fprintf(1, 'Warning: Run ID %d not found in reference table after initial filter. Skipping.\n', scan_runID);
             scan_runID = []; % Mark to skip
             return;
        end
        
        wellsIDsRaw = refTable.Wells_Recorded(idx_run);
        if iscell(wellsIDsRaw) && ~isempty(wellsIDsRaw) % Ensure it's a cell and not empty
            wellsCell = wellsIDsRaw{1}; % Assuming one entry per Run_ID
            if ischar(wellsCell) && contains(wellsCell, ',')
                wellsToProcess = str2double(strsplit(wellsCell, ','));
            elseif isnumeric(wellsCell)
                 wellsToProcess = wellsCell;
            elseif ischar(wellsCell) && ~isempty(str2num(wellsCell)) %#ok<ST2NM>
                wellsToProcess = str2double(wellsCell); % single numeric string
            else
                 fprintf(1, 'Warning: Wells_Recorded for RunID %d has unexpected format: %s. Skipping wells for this file.\n', scan_runID, disp(wellsCell));
                 wellsToProcess = [];
            end
        elseif isnumeric(wellsIDsRaw) % If it's directly a numeric array not in a cell
             wellsToProcess = wellsIDsRaw;
        else
            fprintf(1, 'Warning: Wells_Recorded for RunID %d is empty or has unexpected format. Skipping wells for this file.\n', scan_runID);
            wellsToProcess = [];
        end

        neuronTypesRaw = refTable.NeuronSource(idx_run);
        if iscell(neuronTypesRaw) && ~isempty(neuronTypesRaw)
            neuronStr = neuronTypesRaw{1};
             if ischar(neuronStr) && contains(neuronStr, ',')
                neuronTypes = strtrim(strsplit(neuronStr, ','));
            elseif ischar(neuronStr)
                neuronTypes = {strtrim(neuronStr)}; % Ensure cell array
             else
                neuronTypes = {''}; % Default if format is unexpected
             end
        else
            neuronTypes = {''}; % Default
        end
        
        assayType = refTable.Assay(idx_run); % This will be a cell
        if isempty(assayType), assayType = {''}; end


        % Ensure neuronTypes matches length of wellsToProcess
        if ~isempty(wellsToProcess) && length(neuronTypes) ~= length(wellsToProcess)
            if length(neuronTypes) == 1
                neuronTypes = repmat(neuronTypes, 1, length(wellsToProcess));
            else
                fprintf(1, 'Warning: Mismatch wells and neuron types for RunID %d. Using first type or empty.\n', scan_runID);
                if ~isempty(neuronTypes)
                    neuronTypes = repmat(neuronTypes(1), 1, length(wellsToProcess));
                else
                    neuronTypes = repmat({''}, 1, length(wellsToProcess));
                end
            end
        end
    else
        scan_runID = []; % Mark to skip by emptying runID
    end
end

function [] = createWellPlotFolders(baseOpDir, chipID, wellID, createFRMatrixFolders, frBinSizes)
    % Creates the nested directory structure for saving plots for a specific well.
    folderName = sprintf('%s_Well%d', chipID, wellID);
    chipWellFolder = fullfile(baseOpDir, 'Network_outputs', 'Raster_BurstActivity', folderName);

    outputSubFolders = {'Plot60s', 'Plot120s', 'Plot300s', 'Plot600s'};
    formatSubFolders = {'png', 'eps'}; % 'eps' might be conditional via plottingParams.epsPlot

    if ~isfolder(chipWellFolder), mkdir(chipWellFolder); end

    for i = 1:length(outputSubFolders)
        durationFolder = fullfile(chipWellFolder, outputSubFolders{i});
        if ~isfolder(durationFolder), mkdir(durationFolder); end
        for j = 1:length(formatSubFolders)
            formatPath = fullfile(durationFolder, formatSubFolders{j});
            if ~isfolder(formatPath), mkdir(formatPath); end
        end
    end
    
    % Folders for Firing Rate Matrix and SVD plots (if applicable)
    if createFRMatrixFolders
        baseFRPath = fullfile(baseOpDir, 'Network_outputs', 'FiringRateData', folderName);
        if ~isfolder(baseFRPath), mkdir(baseFRPath); end

        for frBin = frBinSizes
            subNames = {sprintf('FiringRateMatrix_BinSize_%.2f',frBin), ...
                        sprintf('Heatmaps_BinSize_%.2f',frBin), ...
                        sprintf('Log_Heatmaps_BinSize_%.2f',frBin), ...
                        sprintf('SVD_BinSize_%.2f',frBin)};
            for i=1:length(subNames)
                if ~isfolder(fullfile(baseFRPath,subNames{i})), mkdir(fullfile(baseFRPath,subNames{i}));end
            end
        end
    end
end


function [wellData, relSpikeTimes, ts, ch, recordingTime_s, scan_div] = loadWellData(filePath, wellID, div0_date)
    % Loads data for a single well using mxw.fileManager and computes DIV.
    wellData = mxw.fileManager(filePath, wellID);
    relSpikeTimes = mxw.util.computeRelativeSpikeTimes(wellData); % For raster plot
    
    % Get raw spike times and channels for other calculations
    ts = ((double(wellData.fileObj.spikes.frameno) - double(wellData.fileObj.firstFrameNum)) / wellData.fileObj.samplingFreq)';
    ch = wellData.fileObj.spikes.channel;
    
    % Filter out potential negative timestamps or spikes before recording start
    valid_idx = ts >= 0;
    ts = ts(valid_idx);
    ch = ch(valid_idx);

    recordingTime_s = wellData.fileObj.dataLenTime;
    
    hd5_time_str = wellData.fileObj.stopTime; % Or startTime if more appropriate
    try
        hd5Date = datetime(hd5_time_str, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    catch
        hd5Date = datetime(hd5_time_str, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss'); % Fallback
    end
    scan_div = floor(days(hd5Date - div0_date));
end

function [networkActivityStruct] = calculateSmoothedNetworkActivity(wellData, binSize, gaussianSigma)
    % Computes smoothed network activity.
    networkActivityStruct = mxw.networkActivity.computeNetworkAct(wellData, 'BinSize', binSize, 'GaussianSigma', gaussianSigma);
end

function [networkBurstPeakStats] = detectNetworkBursts_PeakProminence(firingRate, timeVector, analysisParams)
    % Detects network bursts using peak prominence. This is essentially your computeNetworkStatsModified.
    % Assuming computeNetworkStatsModified is available and works as in original code:
    networkBurstPeakStats = computeNetworkStatsModified(firingRate, timeVector, ...
        'ThresholdMethod', analysisParams.thresholdFunction, ...
        'Threshold', analysisParams.thresholdBurst, ...
        'MinPeakProminence', analysisParams.minProminence, ...
        'MinPeakDistance', analysisParams.minPeakDistance);
end

function [burstSummaryMetrics, derivedBurstChars] = calculateDerivedNetworkBurstStats(networkActStruct, networkBurstPeakStats, ts, ch, thresholdStartStop, recordingTime_s)
    % Calculates various metrics from detected network burst peaks and raw spikes.
    
    % Initialize outputs
    derivedBurstChars = struct();
    meanIBI = NaN; covIBI = NaN; meanBurstPeak = NaN; covBurstPeak = NaN;
    nBursts = 0; meanSpikesPerBurst = NaN; covSpikesPerBurst = NaN;
    meanAbsBP = NaN; covAbsBP = NaN; meanBurstDuration = NaN; covBurstDuration = NaN;
    baselineFiringRate = mean(networkActStruct.firingRate); % Default
    
    IBI_List_str = {''}; Burst_Peak_List_str = {''}; Abs_Burst_Peak_List_str = {''};
    Burst_Duration_List_str = {''}; Spikes_per_Burst_List_str = {''};
    derivedBurstChars.burstEdges = []; % Initialize to empty

    if isempty(networkBurstPeakStats) || isempty(networkBurstPeakStats.maxAmplitudesTimes) || length(networkBurstPeakStats.maxAmplitudesTimes) < 2 % Need at least 2 bursts for IBI
        nBursts = length(networkBurstPeakStats.maxAmplitudesTimes); % Could be 0 or 1
        if nBursts == 1
            meanBurstPeak = networkBurstPeakStats.maxAmplitudesValues(1);
            % Potentially calculate single burst duration and spikes
        end
        % Keep other metrics as NaN or default
    else
        peakAmps = networkBurstPeakStats.maxAmplitudesValues';
        peakTimes = networkBurstPeakStats.maxAmplitudesTimes;
        nBursts = length(peakTimes);
        burstRate = nBursts / recordingTime_s;

        absBPs = networkActStruct.absfiringRate(ismember(networkActStruct.time, peakTimes));
        
        meanIBI = mean(diff(peakTimes)); % diff(peakTimes) gives IBIs
        stdIBI = std(diff(peakTimes));
        covIBI = (stdIBI / meanIBI) * 100;
        
        meanBurstPeak = mean(peakAmps);
        stdBurstPeak = std(peakAmps);
        covBurstPeak = (stdBurstPeak / meanBurstPeak) * 100;
        
        meanAbsBP = mean(absBPs);
        stdAbsBP = std(absBPs);
        covAbsBP = (stdAbsBP / meanAbsBP) * 100;

        % Determine burst edges and spikes per burst
        burstEdges = zeros(nBursts, 2);
        spikesPerBurstVec = zeros(nBursts, 1);
        burstDurationsVec = zeros(nBursts, 1);

        for i = 1:nBursts
            idx_act = networkActStruct.time > (peakTimes(i) - 6) & networkActStruct.time < (peakTimes(i) + 6);
            t1_act = networkActStruct.time(idx_act);
            a1_act = networkActStruct.firingRate(idx_act)';
            
            peakWidthAmp = peakAmps(i) * (1 - thresholdStartStop); % Threshold relative to peak height from zero
                                                                 % Or: baseline + (peakAmps(i)-baseline)*(1-thresholdStartStop) if baseline is non-zero
            
            idx1 = find(a1_act < peakWidthAmp & t1_act < peakTimes(i), 1, 'last');
            idx2 = find(a1_act < peakWidthAmp & t1_act > peakTimes(i), 1, 'first');
            
            if ~isempty(idx1) && ~isempty(idx2)
                tBefore = t1_act(idx1);
                tAfter = t1_act(idx2);
                burstEdges(i,:) = [tBefore, tAfter];
                burstDurationsVec(i) = tAfter - tBefore;
                spikesPerBurstVec(i) = sum(ts >= tBefore & ts <= tAfter);
            else % Could not define edges for this burst
                burstEdges(i,:) = [NaN, NaN];
                burstDurationsVec(i) = NaN;
                spikesPerBurstVec(i) = NaN;
            end
        end
        derivedBurstChars.burstEdges = burstEdges;

        validBurstIdx = ~isnan(burstDurationsVec);
        spikesPerBurstVec = spikesPerBurstVec(validBurstIdx);
        burstDurationsVec = burstDurationsVec(validBurstIdx);
        
        if ~isempty(spikesPerBurstVec)
            meanSpikesPerBurst = mean(spikesPerBurstVec);
            stdSpikesPerBurst = std(spikesPerBurstVec);
            covSpikesPerBurst = (stdSpikesPerBurst / meanSpikesPerBurst) * 100;
        end
        if ~isempty(burstDurationsVec)
            meanBurstDuration = mean(burstDurationsVec);
            stdBurstDuration = std(burstDurationsVec);
            covBurstDuration = (stdBurstDuration / meanBurstDuration) * 100;
        end

        % Recalculate baseline firing rate excluding burst periods
        if ~isempty(burstEdges) && any(validBurstIdx)
            excludeTimeMask = false(size(networkActStruct.time));
            validEdges = burstEdges(validBurstIdx,:);
            for i = 1:size(validEdges, 1)
                excludeTimeMask = excludeTimeMask | (networkActStruct.time >= validEdges(i,1) & networkActStruct.time <= validEdges(i,2));
            end
            includeTimeMask = ~excludeTimeMask;
            if any(includeTimeMask)
                baselineFiringRate = mean(networkActStruct.firingRate(includeTimeMask));
            else % All time is within bursts, baseline is effectively zero or undefined for this calc
                baselineFiringRate = 0; 
            end
        end
        
        % String versions for table
        IBI_List_str = {strjoin(arrayfun(@num2str, diff(peakTimes), 'UniformOutput', false), ',')};
        Burst_Peak_List_str = {strjoin(arrayfun(@num2str, peakAmps, 'UniformOutput', false), ',')};
        Abs_Burst_Peak_List_str = {strjoin(arrayfun(@num2str, absBPs, 'UniformOutput', false), ',')};
        Burst_Duration_List_str = {strjoin(arrayfun(@num2str, burstDurationsVec', 'UniformOutput', false), ',')}; % Ensure row vector for strjoin
        Spikes_per_Burst_List_str = {strjoin(arrayfun(@num2str, spikesPerBurstVec', 'UniformOutput', false), ',')};
    end

    burstSummaryMetrics = table(meanIBI, covIBI, meanBurstPeak, covBurstPeak, ...
        nBursts, burstRate, meanSpikesPerBurst, covSpikesPerBurst, meanAbsBP, covAbsBP, meanBurstDuration, covBurstDuration, baselineFiringRate, ...
        IBI_List_str, Burst_Peak_List_str, Abs_Burst_Peak_List_str, Burst_Duration_List_str, Spikes_per_Burst_List_str, ...
        'VariableNames', {
        'mean_IBI', 'cov_IBI', 'mean_Burst_Peak', 'cov_Burst_Peak', ...
        'Number_Bursts', 'BurstRate', 'mean_Spike_per_Burst', 'cov_Spike_per_Burst', 'mean_Burst_Peak_Abs', 'cov_Burst_Peak_Abs', 'mean_BurstDuration', 'cov_BurstDuration','Baseline', ...
        'IBI_List', 'Burst_Peak_List', 'Abs_Burst_Peak_List', 'Burst_Times_List', 'SpikesPerBurst_List' % Note: Burst_Times_List was BurstDuration, SpikesPerBurst_List was spikesPerBurst
        });
end


function [extMetrics] = calculateExtendedWellMetrics(ts, ch, burstEdges, recordingTime_s) % Added recordingTime_s
    % Calculates extended metrics like ISI distributions, Fano factor.
    % Note: Fano factor typically requires windowed spike counts.
    % The original Fano calculation seemed to be on the entire recording's spike times binned at 100ms.
    
    extMetrics = struct();

    % Overall ISI stats
    if length(ts) > 1
        all_ISIs_flat = [];
        unique_ch = unique(ch);
        for i_ch = 1:length(unique_ch)
            ch_ts = sort(ts(ch == unique_ch(i_ch)));
            if length(ch_ts) > 1
                all_ISIs_flat = [all_ISIs_flat; diff(ch_ts)]; %#ok<AGROW>
            end
        end
        all_ISIs_flat = all_ISIs_flat(all_ISIs_flat > 0 & all_ISIs_flat < 1); % Filter, e.g. < 1s
        
        if ~isempty(all_ISIs_flat)
            extMetrics.MeanNetworkISI = mean(all_ISIs_flat);
            extMetrics.StdNetworkISI = std(all_ISIs_flat);
            extMetrics.CoVNetworkISI = (extMetrics.StdNetworkISI / extMetrics.MeanNetworkISI) * 100;
            
            instFreq = 1./all_ISIs_flat;
            minFreq = min(instFreq(instFreq > 0 & isfinite(instFreq)));
            maxFreq = max(instFreq(isfinite(instFreq)));
            if ~isempty(minFreq) && ~isempty(maxFreq) && minFreq < maxFreq
                iqrFreq = iqr(instFreq(isfinite(instFreq)));
                 binWidth = 2 * iqrFreq / (length(instFreq(isfinite(instFreq)))^(1/3));
                 if binWidth > 0
                    numBins = ceil((maxFreq-minFreq) / binWidth);
                    numBins = max(numBins,5); % ensure some bins
                    binedgesFreq = logspace(log10(minFreq), log10(min(maxFreq,500)), numBins); % Cap at 500Hz
                    [extMetrics.networkAPFreqBins, extMetrics.networkAPFreqEdges] = histcounts(instFreq, binedgesFreq);
                 else
                    extMetrics.networkAPFreqBins = []; extMetrics.networkAPFreqEdges = [];
                 end
            else
                 extMetrics.networkAPFreqBins = []; extMetrics.networkAPFreqEdges = [];
            end
        else
            extMetrics.MeanNetworkISI = NaN; extMetrics.CoVNetworkISI = NaN;
            extMetrics.networkAPFreqBins = []; extMetrics.networkAPFreqEdges = [];
        end
    else
        extMetrics.MeanNetworkISI = NaN; extMetrics.CoVNetworkISI = NaN;
        extMetrics.networkAPFreqBins = []; extMetrics.networkAPFreqEdges = [];
    end

    % Within-burst and Outside-burst ISIs
    tsWithinBurst = []; chWithinBurst = [];
    tsOutsideBurst = ts; chOutsideBurst = ch; % Start with all, then remove those in bursts

    if ~isempty(burstEdges) && size(burstEdges,2)==2
        isSpikeInAnyBurst = false(size(ts));
        for i_burst = 1:size(burstEdges, 1)
            if ~any(isnan(burstEdges(i_burst,:)))
                 isSpikeInAnyBurst(ts >= burstEdges(i_burst,1) & ts <= burstEdges(i_burst,2)) = true;
            end
        end
        tsWithinBurst = ts(isSpikeInAnyBurst);
        chWithinBurst = ch(isSpikeInAnyBurst);
        tsOutsideBurst = ts(~isSpikeInAnyBurst);
        chOutsideBurst = ch(~isSpikeInAnyBurst);
    end
    
    % Calculate for tsWithinBurst
    if length(tsWithinBurst) > 1
        within_burst_ISIs_flat = [];
        unique_ch_wb = unique(chWithinBurst);
        for i_ch = 1:length(unique_ch_wb)
            ch_ts_wb = sort(tsWithinBurst(chWithinBurst == unique_ch_wb(i_ch)));
            if length(ch_ts_wb) > 1
                within_burst_ISIs_flat = [within_burst_ISIs_flat; diff(ch_ts_wb)]; %#ok<AGROW>
            end
        end
        within_burst_ISIs_flat = within_burst_ISIs_flat(within_burst_ISIs_flat > 0 & within_burst_ISIs_flat < 1);

        if ~isempty(within_burst_ISIs_flat)
            extMetrics.MeanWithinBurstISI = mean(within_burst_ISIs_flat);
            extMetrics.StdWithinBurstISI = std(within_burst_ISIs_flat);
            extMetrics.CoVWithinBurstISI = (extMetrics.StdWithinBurstISI / extMetrics.MeanWithinBurstISI) * 100;

            instAPFreqWithinBurst = 1./within_burst_ISIs_flat;
            minFreq = min(instAPFreqWithinBurst(instAPFreqWithinBurst > 0 & isfinite(instAPFreqWithinBurst)));
            maxFreq = max(instAPFreqWithinBurst(isfinite(instAPFreqWithinBurst)));
             if ~isempty(minFreq) && ~isempty(maxFreq) && minFreq < maxFreq
                iqrFreq = iqr(instAPFreqWithinBurst(isfinite(instAPFreqWithinBurst)));
                binWidth = 2 * iqrFreq / (length(instAPFreqWithinBurst(isfinite(instAPFreqWithinBurst)))^(1/3));
                 if binWidth > 0
                    numBins = ceil((maxFreq-minFreq) / binWidth);
                    numBins = max(numBins,5);
                    binedgesAP = logspace(log10(minFreq), log10(min(maxFreq,500)), numBins);
                    [extMetrics.burstAPFreqBins, extMetrics.burstAPFreqEdges] = histcounts(instAPFreqWithinBurst, binedgesAP);
                else
                    extMetrics.burstAPFreqBins = []; extMetrics.burstAPFreqEdges = [];
                end
            else
                extMetrics.burstAPFreqBins = []; extMetrics.burstAPFreqEdges = [];
            end
        else
            extMetrics.MeanWithinBurstISI = NaN; extMetrics.CoVWithinBurstISI = NaN;
            extMetrics.burstAPFreqBins = []; extMetrics.burstAPFreqEdges = [];
        end
    else
        extMetrics.MeanWithinBurstISI = NaN; extMetrics.CoVWithinBurstISI = NaN;
        extMetrics.burstAPFreqBins = []; extMetrics.burstAPFreqEdges = [];
    end

    % Calculate for tsOutsideBurst
    if length(tsOutsideBurst) > 1
        outside_burst_ISIs_flat = [];
        unique_ch_ob = unique(chOutsideBurst);
        for i_ch = 1:length(unique_ch_ob)
            ch_ts_ob = sort(tsOutsideBurst(chOutsideBurst == unique_ch_ob(i_ch)));
            if length(ch_ts_ob) > 1
                outside_burst_ISIs_flat = [outside_burst_ISIs_flat; diff(ch_ts_ob)]; %#ok<AGROW>
            end
        end
        outside_burst_ISIs_flat = outside_burst_ISIs_flat(outside_burst_ISIs_flat > 0 & outside_burst_ISIs_flat < 1);

        if ~isempty(outside_burst_ISIs_flat)
            extMetrics.MeanOutsideBurstISI = mean(outside_burst_ISIs_flat);
            extMetrics.StdOutsideBurstISI = std(outside_burst_ISIs_flat);
            extMetrics.CoVOutsideBurstISI = (extMetrics.StdOutsideBurstISI / extMetrics.MeanOutsideBurstISI) * 100;
            
            instAPFreqOutsideBurst = 1./outside_burst_ISIs_flat;
            minFreq = min(instAPFreqOutsideBurst(instAPFreqOutsideBurst > 0 & isfinite(instAPFreqOutsideBurst)));
            maxFreq = max(instAPFreqOutsideBurst(isfinite(instAPFreqOutsideBurst)));
            if ~isempty(minFreq) && ~isempty(maxFreq) && minFreq < maxFreq
                iqrFreq = iqr(instAPFreqOutsideBurst(isfinite(instAPFreqOutsideBurst)));
                binWidth = 2 * iqrFreq / (length(instAPFreqOutsideBurst(isfinite(instAPFreqOutsideBurst)))^(1/3));
                if binWidth > 0
                    numBins = ceil((maxFreq-minFreq) / binWidth);
                    numBins = max(numBins,5);
                    binedgesAPOutside = logspace(log10(minFreq), log10(min(maxFreq,500)), numBins);
                    [extMetrics.nonburstAPFreqBins, extMetrics.nonburstAPFreqEdges] = histcounts(instAPFreqOutsideBurst, binedgesAPOutside);
                else
                    extMetrics.nonburstAPFreqBins = []; extMetrics.nonburstAPFreqEdges = [];
                end
            else
                extMetrics.nonburstAPFreqBins = []; extMetrics.nonburstAPFreqEdges = [];
            end
        else
            extMetrics.MeanOutsideBurstISI = NaN; extMetrics.CoVOutsideBurstISI = NaN;
            extMetrics.nonburstAPFreqBins = []; extMetrics.nonburstAPFreqEdges = [];
        end
    else
        extMetrics.MeanOutsideBurstISI = NaN; extMetrics.CoVOutsideBurstISI = NaN;
        extMetrics.nonburstAPFreqBins = []; extMetrics.nonburstAPFreqEdges = [];
    end

    % Fano Factor (using 100ms bins across the whole recording for all spikes)
    if ~isempty(ts) && recordingTime_s > 0.1
        fanoBinSize = 0.1; % 100 ms
        binedgesFano = 0:fanoBinSize:recordingTime_s;
        if length(binedgesFano) > 1
            spikeCountsFano = histcounts(ts, binedgesFano);
            meanSpikeCounts = mean(spikeCountsFano);
            varSpikeCounts = var(spikeCountsFano);
            if meanSpikeCounts > 0
                extMetrics.Fanofactor = varSpikeCounts / meanSpikeCounts;
            else
                extMetrics.Fanofactor = NaN;
            end
        else
            extMetrics.Fanofactor = NaN;
        end
    else
        extMetrics.Fanofactor = NaN;
    end
end

function [perChannelBurstData] = analyzePerChannelBursts_logISI(ts, ch, min_spikes, min_duration, logISI_tune_params)
    % Analyzes bursts per channel using the logISI method.
    unique_channels = unique(ch);
    num_unique_channels = length(unique_channels);
    perChannelBurstData = cell(1, num_unique_channels);

    for ch_idx = 1:num_unique_channels
        current_channel_id = unique_channels(ch_idx);
        channel_spike_times = sort(ts(ch == current_channel_id)); % Ensure sorted
        
        if length(channel_spike_times) >= min_spikes
            channel_bursts_struct = detectBursts_logISI(channel_spike_times, min_spikes, min_duration, logISI_tune_params);
            
            if ~isempty(channel_bursts_struct)
                perChannelBurstData{ch_idx} = struct(...
                    'channel_id', current_channel_id, ...
                    'num_bursts', length(channel_bursts_struct), ...
                    'mean_duration_s', mean([channel_bursts_struct.duration_s]), ...
                    'mean_spikes_per_burst', mean([channel_bursts_struct.num_spikes]), ...
                    'mean_peak_isi_s', mean([channel_bursts_struct.peak_isi_s]), ... % Can be NaN if no clear peak
                    'mean_threshold_isi_s', mean([channel_bursts_struct.threshold_isi_s]), ...
                    'all_bursts_details', channel_bursts_struct);
            else
                perChannelBurstData{ch_idx} = struct('channel_id', current_channel_id, 'num_bursts', 0, 'all_bursts_details', []);
            end
        else
            perChannelBurstData{ch_idx} = struct('channel_id', current_channel_id, 'num_bursts', 0, 'all_bursts_details', []);
        end
    end
    perChannelBurstData = [perChannelBurstData{:}]; % Convert cell array of structs to struct array
end


function [bursts] = detectBursts_logISI(spike_times, min_spikes_in_burst, min_burst_duration_s, isi_threshold_method_params)
% detectBursts_logISI: Detects bursts in a spike train using the logISI histogram method.
% (This function was provided in the previous response, ensure it's included here or in a separate file on the path)
% ... (implementation from previous response) ...
    bursts = struct('start_time', {}, 'end_time', {}, 'duration_s', {}, ...
                    'num_spikes', {}, 'spike_times', {}, ...
                    'peak_isi_s', {}, 'threshold_isi_s', {});

    if nargin < 4 || isempty(isi_threshold_method_params)
        isi_threshold_method_params = struct();
    end
    % Set defaults if not provided
    if ~isfield(isi_threshold_method_params, 'log_isi_bins'), isi_threshold_method_params.log_isi_bins = 50; end
    if ~isfield(isi_threshold_method_params, 'smooth_span_hist'), isi_threshold_method_params.smooth_span_hist = 5; end
    if ~isfield(isi_threshold_method_params, 'max_intra_burst_isi_s'), isi_threshold_method_params.max_intra_burst_isi_s = 0.2; end

    if length(spike_times) < min_spikes_in_burst
        return; 
    end

    isis = diff(spike_times);
    if isempty(isis) || length(isis) < min_spikes_in_burst -1 % Need at least min_spikes_in_burst-1 ISIs for that many spikes
        % Check if the few spikes present could form a single short burst
        if length(spike_times) >= min_spikes_in_burst
            duration = spike_times(end) - spike_times(1);
            if duration >= min_burst_duration_s && all(isis <= isi_threshold_method_params.max_intra_burst_isi_s)
                 bursts(1).start_time = spike_times(1);
                 bursts(1).end_time = spike_times(end);
                 bursts(1).duration_s = duration;
                 bursts(1).num_spikes = length(spike_times);
                 bursts(1).spike_times = spike_times;
                 bursts(1).peak_isi_s = ifelse(isempty(isis), NaN, mean(isis)); 
                 bursts(1).threshold_isi_s = isi_threshold_method_params.max_intra_burst_isi_s;
            end
        end
        return;
    end
    
    positive_isis = isis(isis > 1e-9); % Avoid log(0) or very small numbers; 1ns precision assumed
    if isempty(positive_isis) || length(positive_isis) < min_spikes_in_burst -1
        return;
    end
    log_isis = log10(positive_isis); 
    
    min_log_isi = min(log_isis);
    max_log_isi = max(log_isis);

    if isinf(min_log_isi) || isinf(max_log_isi) || (max_log_isi - min_log_isi < 1e-6) % Effectively all ISIs are same or problematic
        if length(isis) >= min_spikes_in_burst-1 && all(isis <= isi_threshold_method_params.max_intra_burst_isi_s)
            duration = spike_times(end) - spike_times(1);
            if duration >= min_burst_duration_s
                bursts(1).start_time = spike_times(1);
                bursts(1).end_time = spike_times(end);
                bursts(1).duration_s = duration;
                bursts(1).num_spikes = length(spike_times);
                bursts(1).spike_times = spike_times;
                bursts(1).peak_isi_s = mean(isis);
                bursts(1).threshold_isi_s = isi_threshold_method_params.max_intra_burst_isi_s;
            end
        end
        return;
    end

    hist_edges = linspace(min_log_isi, max_log_isi, isi_threshold_method_params.log_isi_bins + 1);
    [counts, edges_hist] = histcounts(log_isis, hist_edges); % Changed variable name
    centers = (edges_hist(1:end-1) + edges_hist(2:end)) / 2;

    smoothed_counts = smoothdata(counts, 'movmean', isi_threshold_method_params.smooth_span_hist);
    
    final_peak_isi = NaN;
    final_isi_threshold = isi_threshold_method_params.max_intra_burst_isi_s; % Default

    [peak_values, peak_locs_indices] = findpeaks(smoothed_counts);
    if ~isempty(peak_locs_indices)
        % Find first peak that corresponds to ISI <= max_intra_burst_isi_s
        candidate_intra_burst_peaks_idx = peak_locs_indices(centers(peak_locs_indices) <= log10(isi_threshold_method_params.max_intra_burst_isi_s));
        if ~isempty(candidate_intra_burst_peaks_idx)
            % Select the one with highest count among these candidates
            [~, max_val_idx] = max(peak_values(ismember(peak_locs_indices, candidate_intra_burst_peaks_idx)));
            main_peak_hist_idx = candidate_intra_burst_peaks_idx(max_val_idx);
            final_peak_isi = 10^centers(main_peak_hist_idx);

            % Find first valley *after* this main_peak_hist_idx
            temp_counts_for_valley = smoothed_counts;
            temp_counts_for_valley(1:main_peak_hist_idx) = Inf; % Ignore valleys before or at the main peak
            [~, valley_locs_indices] = findpeaks(-temp_counts_for_valley); % Find peaks in inverted signal

            if ~isempty(valley_locs_indices)
                first_valley_hist_idx = valley_locs_indices(1);
                final_isi_threshold = 10^centers(first_valley_hist_idx);
                 % Heuristic: if threshold is too large, cap it or use fallback.
                if final_isi_threshold > isi_threshold_method_params.max_intra_burst_isi_s * 2 % e.g. 2x the expected max
                    final_isi_threshold = isi_threshold_method_params.max_intra_burst_isi_s;
                end
            else
                 % Fallback: if no valley, threshold might be to the right of the peak
                 % This part is tricky and often where logISI methods differ in robustness
                 idx_half_max_right = find(smoothed_counts(main_peak_hist_idx:end) < smoothed_counts(main_peak_hist_idx)/2, 1, 'first');
                 if ~isempty(idx_half_max_right) && (main_peak_hist_idx + idx_half_max_right -1) <= length(centers)
                     final_isi_threshold = 10^(centers(main_peak_hist_idx + idx_half_max_right -1));
                 else % Fallback to user-defined max or a slightly wider version
                     final_isi_threshold = isi_threshold_method_params.max_intra_burst_isi_s;
                 end
            end
        end % else, final_isi_threshold remains default
    end % else, final_isi_threshold remains default

    % --- Identify bursts based on the final_isi_threshold ---
    in_burst_flag = false;
    current_burst_start_spike_idx = 0; % Index in spike_times
    burst_counter = 0;

    for i = 1:length(isis) % Iterate through ISIs
        if isis(i) <= final_isi_threshold
            if ~in_burst_flag
                in_burst_flag = true;
                current_burst_start_spike_idx = i; % The burst starts with spike_times(i)
            end
        else % ISI is greater than threshold, potential end of burst
            if in_burst_flag
                % Burst ended with spike_times(i). Spikes from current_burst_start_spike_idx to i.
                num_s = (i - current_burst_start_spike_idx) + 1; % Number of spikes in this potential burst
                
                if num_s >= min_spikes_in_burst
                    b_start_time = spike_times(current_burst_start_spike_idx);
                    b_end_time = spike_times(i); % ISI(i) was between spike(i) and spike(i+1)
                    b_duration = b_end_time - b_start_time;

                    if b_duration >= min_burst_duration_s
                        burst_counter = burst_counter + 1;
                        bursts(burst_counter).start_time = b_start_time;
                        bursts(burst_counter).end_time = b_end_time;
                        bursts(burst_counter).duration_s = b_duration;
                        bursts(burst_counter).num_spikes = num_s;
                        bursts(burst_counter).spike_times = spike_times(current_burst_start_spike_idx:i);
                        bursts(burst_counter).peak_isi_s = final_peak_isi; % Store determined peak ISI
                        bursts(burst_counter).threshold_isi_s = final_isi_threshold;
                    end
                end
                in_burst_flag = false;
            end
        end
    end

    % Check for a burst that extends to the end of the spike train
    if in_burst_flag
        % Burst includes spikes from current_burst_start_spike_idx to end of spike_times
        num_s = (length(spike_times) - current_burst_start_spike_idx) + 1;
        if num_s >= min_spikes_in_burst
            b_start_time = spike_times(current_burst_start_spike_idx);
            b_end_time = spike_times(end);
            b_duration = b_end_time - b_start_time;

            if b_duration >= min_burst_duration_s
                burst_counter = burst_counter + 1;
                bursts(burst_counter).start_time = b_start_time;
                bursts(burst_counter).end_time = b_end_time;
                bursts(burst_counter).duration_s = b_duration;
                bursts(burst_counter).num_spikes = num_s;
                bursts(burst_counter).spike_times = spike_times(current_burst_start_spike_idx:end);
                bursts(burst_counter).peak_isi_s = final_peak_isi;
                bursts(burst_counter).threshold_isi_s = final_isi_threshold;
            end
        end
    end
end


function [] = plotAllWellFigures(wellData, relSpikeTimes, networkActStruct, networkBurstPeakStats, plotParams, analysisParams, baselineFiringRate)
    % Master plotting function for a single well.
    % plotParams should contain opDir, scan_chipID, wellID, scan_div, neuronSourceType, xlimNetwork, ylimNetwork etc.
    
    basePlotDir = fullfile(plotParams.opDir, 'Network_outputs', 'Raster_BurstActivity', ...
                           sprintf('%s_Well%d', plotParams.scan_chipID, plotParams.wellID));

    timeWindows_s = [60, 120, 300, 600];
    plotFileNameBase = sprintf('Raster_BurstActivity_%s_WellID_%d_%s_DIV%d_%s', ...
        plotParams.scan_runID_text, plotParams.wellID, plotParams.scan_chipID, plotParams.scan_div, strrep(plotParams.neuronSourceType, ' ', ''));

    % General Raster & Network Activity Plot
    fig_main = figure('Color','w', 'Position', [0 0 400 800], 'Visible', 'off');
    ax1 = subplot(2,1,1);
    mxw.plot.rasterPlot(wellData, 'Figure', false, 'Axes', ax1); % Pass axes
    box off;
    ylabel(ax1,'Channel ID');
    
    ax2 = subplot(2,1,2);
    plotNetworkActivityModified(networkActStruct, ... % Assuming this is your plotting function
        'ThresholdFunction', networkBurstPeakStats.thresholdFunction, ...
        'Threshold', networkBurstPeakStats.threshold, 'Figure', false, 'Axes', ax2); % Pass axes
    hold(ax2, 'on');
    if ~isempty(networkBurstPeakStats.maxAmplitudesTimes)
        plot(ax2, networkBurstPeakStats.maxAmplitudesTimes, networkBurstPeakStats.maxAmplitudesValues, 'or', 'MarkerFaceColor', 'r');
    end
    % Plot baseline firing rate
    if ~isnan(baselineFiringRate) && baselineFiringRate > 0
        line(ax2, xlim(ax2), [baselineFiringRate baselineFiringRate], 'Color', 'g', 'LineStyle', '--', 'DisplayName', 'Baseline FR');
    end
    hold(ax2, 'off');
    box off; 
    xlabel(ax2,'Time (s)'); 
    ylabel(ax2,'Network Firing Rate (Hz)');
    
    % Save for different time windows
    for tWin = timeWindows_s
        xlim(ax1, [0 tWin]);
        if ~isempty(relSpikeTimes.channel) % only set ylim if channels exist
            ylim(ax1, [0.5 max(relSpikeTimes.channel)+0.5]); % Adjusted ylim
        else
            ylim(ax1, [0.5 1.5]); % Default if no spikes
        end

        xlim(ax2, [0 tWin]);
        ylim(ax2, [0 plotParams.ylimNetwork]); % Use passed ylim
        
        if tWin == 300 % Add annotation only for 300s plot
            textString = sprintf('#NBursts: %d | Mean NBDur: %.2fs | Mean Sp/NB: %.2f | Mean IBI: %.2fs | Mean NBPk: %.2fHz', ...
                plotParams.nBursts, plotParams.meanBurstDuration, plotParams.meanSpikesPerBurst, plotParams.meanIBI, plotParams.meanBurstPeak);
            annotation(fig_main,'textbox', [0.1, 0.47, 0.8, 0.05], 'String', textString, ... % Adjusted position slightly
                       'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FitBoxToText','on');
        end
        
        % PNG
        pngDir = fullfile(basePlotDir, sprintf('Plot%ds', tWin), 'png');
        print(fig_main, fullfile(pngDir, [plotFileNameBase '.png']), '-dpng', '-r300');
        % EPS (optional)
        if plotParams.epsPlot
            epsDir = fullfile(basePlotDir, sprintf('Plot%ds', tWin), 'eps');
            print(fig_main, fullfile(epsDir, [plotFileNameBase '.eps']), '-depsc', '-r300');
        end
    end
    
    % Save the .fig file once (e.g., with the longest view or a default view)
    xlim(ax1, [0 plotParams.xlimNetwork]); % Reset to overall xlim for .fig
    xlim(ax2, [0 plotParams.xlimNetwork]);
    figSavePath = fullfile(plotParams.opDir, 'Network_outputs', 'Raster_BurstActivity', [plotFileNameBase '.fig']);
    savefig(fig_main, figSavePath);
    close(fig_main);

    % Amplitude Map Plot
    try
        amplitude90perc = abs(mxw.activityMap.computeAmplitude90percentile(wellData));
        if isempty(amplitude90perc) || all(amplitude90perc == 0)
             error('Amplitude data is empty or all zeros.');
        end
        maxAmp = mxw.util.percentile(amplitude90perc(amplitude90perc~=0 & isfinite(amplitude90perc)),99);
        if isempty(maxAmp) || maxAmp == 0, maxAmp = 100; end % Fallback

        fig_amp = figure('Color','w','visible','off', 'Position', [0 0 400 600]);
        subplot(2,1,1);
        % Ensure map coordinates are present
        if isfield(wellData, 'processedMap') && isfield(wellData.processedMap,'xpos') && ~isempty(wellData.processedMap.xpos)
            scatter(wellData.processedMap.xpos, wellData.processedMap.ypos, 20, amplitude90perc ,'s', 'filled');
            colormap(gca, 'hot'); % Or 'parula'
            cb = colorbar; ylabel(cb, '90th Pctl Amp (\muV)');
            caxis([0 maxAmp]); % Cap color axis
        else
            text(0.5,0.5, 'No map data for scatter', 'HorizontalAlignment','center');
        end
        
        chipX=min(wellData.electrodeMap.x); chipY=min(wellData.electrodeMap.y); 
        chipW=max(wellData.electrodeMap.x)-min(wellData.electrodeMap.x); 
        chipH=max(wellData.electrodeMap.y)-min(wellData.electrodeMap.y);
        
        hold on;
        rectangle('Position', [chipX, chipY, chipW, chipH], 'EdgeColor', 'k', 'LineWidth', 0.5);
        axis tight; axis equal; axis off;
        title(sprintf('Amplitude Map: Well %d, DIV %d', plotParams.wellID, plotParams.scan_div));

        subplot(2,1,2);
        thrAmp = 5; % Min amplitude to include in histogram
        validAmps = amplitude90perc(amplitude90perc > thrAmp & isfinite(amplitude90perc));
        if ~isempty(validAmps)
            histogram(validAmps, ceil(0:max(1,maxAmp/50):maxAmp)); % Adjust binning
            xlim([thrAmp maxAmp]);
            meanValidAmp = mean(validAmps); stdValidAmp = std(validAmps);
            legend(sprintf('Mean Amp = %.2f \\muV, SD = %.2f', meanValidAmp, stdValidAmp), 'Location', 'best');
        else
            text(0.5,0.5, 'No valid amplitudes for histogram', 'HorizontalAlignment','center');
        end
        ylabel('Counts'); xlabel('90th Pctl Spike Amplitude (\muV)');
        box off;
        
        ampMapBaseDir = fullfile(plotParams.opDir, 'Network_outputs', 'ElectrodesAmplitudeMap', ...
                                 sprintf('%s_Well%d', plotParams.scan_chipID, plotParams.wellID));
        if ~isfolder(ampMapBaseDir), mkdir(ampMapBaseDir); end
        ampPlotFileName = sprintf('AmpMap_%s_WellID_%d_DIV%d_%s.png', ...
            plotParams.scan_runID_text, plotParams.wellID, plotParams.scan_div, strrep(plotParams.neuronSourceType, ' ', ''));
        print(fig_amp, fullfile(ampMapBaseDir, ampPlotFileName), '-dpng', '-r300');
        close(fig_amp);

    catch ME_amp
        fprintf(1, 'Unable to plot Amplitude Map for Run %s, Well %d: %s\n', plotParams.scan_runID_text, plotParams.wellID, ME_amp.message);
        fprintf(plotParams.logFile, 'Unable to plot Amplitude Map for Run %s, Well %d: %s\n', plotParams.scan_runID_text, plotParams.wellID, ME_amp.message);
    end
end


function [] = plotFiringRateSVD(ts, ch, recordingTime_s, frBinSizesSVD, plotParams)
    % Calculates and plots firing rate matrices and their SVD.
    % plotParams should contain opDir, scan_chipID, wellID, scan_div, neuronSourceType etc.
    baseOutputDir = fullfile(plotParams.opDir, 'Network_outputs', 'FiringRateData', ...
                             sprintf('%s_Well%d', plotParams.scan_chipID, plotParams.wellID));

    for frBinSize = frBinSizesSVD
        binEdges = 0:frBinSize:recordingTime_s;
        if length(binEdges) < 2, continue; end % Not enough bins

        % Initialize firingRateMatrix, ensure max(ch) is valid
        maxChannel = max(ch);
        if isempty(maxChannel) || maxChannel == 0, continue; end % No channels or spikes
        firingRateMatrix = zeros(maxChannel, length(binEdges) - 1);

        for chIdx = 1:maxChannel
            spikeTimesThisChannel = ts(ch == chIdx);
            if ~isempty(spikeTimesThisChannel)
                firingRateMatrix(chIdx, :) = histcounts(spikeTimesThisChannel, binEdges);
            end
        end
        if isempty(firingRateMatrix), continue; end

        % --- SVD ---
        try
            [U, S, V] = svd(firingRateMatrix, 'econ');
            singularValues = diag(S);
            numSingularValuesPlot = min(20, length(singularValues));

            fig_svd = figure('Visible', 'off', 'Position', [0 0 600 900]);
            subplot(3, 1, 1);
            plot(singularValues(1:numSingularValuesPlot), '-o', 'LineWidth', 1.5);
            xlabel('Component Number'); ylabel('Singular Value');
            title(sprintf('SVD - Singular Values (Bin Size: %.2fs)', frBinSize)); grid on;

            subplot(3, 1, 2);
            imagesc(U(:, 1:numSingularValuesPlot)); colorbar;
            xlabel('Singular Vector Index'); ylabel('Channel'); title('Left Singular Vectors (U)');

            subplot(3, 1, 3);
            imagesc(V(:, 1:numSingularValuesPlot)'); colorbar; % Transpose V for channels x components like U
            xlabel('Time Bin'); ylabel('Singular Vector Index'); title('Right Singular Vectors (V)'); % Corrected labels if V is time x components

            svdPlotDir = fullfile(baseOutputDir, sprintf('SVD_BinSize_%.2f', frBinSize));
            svdPlotFileName = sprintf('SVD_RunID%s_DIV%d.png', plotParams.scan_runID_text, plotParams.scan_div);
            print(fig_svd, fullfile(svdPlotDir, svdPlotFileName), '-dpng', '-r300');
            close(fig_svd);
        catch ME_svd
             fprintf(1,'SVD Error for Well %d, BinSize %.2f: %s\n', plotParams.wellID, frBinSize, ME_svd.message);
             fprintf(plotParams.logFile,'SVD Error for Well %d, BinSize %.2f: %s\n', plotParams.wellID, frBinSize, ME_svd.message);
        end

        % --- Heatmaps ---
        timeIntervals = [60, 120, 300]; % Up to 300s for heatmaps
        for tLimit = timeIntervals
            if recordingTime_s < tLimit && tLimit ~= timeIntervals(end) % if rec time less than current tlimit, just use full recording for last interval
                 tActualLimit = recordingTime_s;
            else
                 tActualLimit = min(tLimit, recordingTime_s);
            end
            if tActualLimit < frBinSize, continue; end % Not enough time for even one bin
            
            tBinIdx = find(binEdges <= tActualLimit, 1, 'last') - 1;
            if isempty(tBinIdx) || tBinIdx == 0, continue; end
            
            currentFRMatrixPortion = firingRateMatrix(:, 1:tBinIdx);

            % Normal Heatmap
            fig_hm = figure('Visible', 'off');
            imagesc(currentFRMatrixPortion); colormap('hot'); colorbar;
            title(sprintf('Firing Rate (First %ds, Bin %.2fs)', round(tActualLimit), frBinSize));
            xlabel('Time Bins'); ylabel('Channel');
            hmDir = fullfile(baseOutputDir, sprintf('Heatmaps_BinSize_%.2f', frBinSize));
            hmFileName = sprintf('Heatmap_RunID%s_DIV%d_%ds.png', plotParams.scan_runID_text, plotParams.scan_div, round(tActualLimit));
            print(fig_hm, fullfile(hmDir, hmFileName), '-dpng', '-r300');
            close(fig_hm);

            % Log-Scaled Heatmap
            fig_loghm = figure('Visible', 'off');
            imagesc(log10(currentFRMatrixPortion + 1)); colormap('hot'); colorbar;
            title(sprintf('Log FR (First %ds, Bin %.2fs)', round(tActualLimit), frBinSize));
            xlabel('Time Bins'); ylabel('Channel');
            logHmDir = fullfile(baseOutputDir, sprintf('Log_Heatmaps_BinSize_%.2f', frBinSize));
            logHmFileName = sprintf('LogHeatmap_RunID%s_DIV%d_%ds.png',plotParams.scan_runID_text, plotParams.scan_div, round(tActualLimit));
            print(fig_loghm, fullfile(logHmDir, logHmFileName), '-dpng', '-r300');
            close(fig_loghm);
        end
    end
end

function [] = saveCompiledResults(finalTable, csvFilePath, extMetricsFlag, finalExtMetrics, matFilePath)
    % Saves the compiled results table to CSV and extended metrics to .mat.
    if ~isempty(finalTable)
        if isfile(csvFilePath), delete(csvFilePath); end % Overwrite if exists
        writetable(finalTable, csvFilePath);
        fprintf(1, 'Successfully wrote compiled network data to %s\n', csvFilePath);
    else
        fprintf(1, 'Warning: Final compiled table is empty. Nothing written to CSV.\n');
    end

    if extMetricsFlag && ~isempty(finalExtMetrics)
        allExtMetrics = finalExtMetrics; % Ensure variable name matches save command
        save(matFilePath, 'allExtMetrics');
        fprintf(1, 'Successfully wrote extended metrics to %s\n', matFilePath);
    elseif extMetricsFlag && isempty(finalExtMetrics)
        fprintf(1, 'Warning: Extended metrics flag was true, but no extended metrics were aggregated to save.\n');
    end
end

function [emptyRow] = createEmptyMetricsRow(run_id, div, assay, well_id, neuron_type, timestamp, chip_id, reason)
    % Creates a table row with NaNs for a skipped or empty well.
    % This function needs to match the columns of your 'burstSummaryMetrics' + metadata
    % For simplicity, I'm creating a basic one. You'll need to expand it.
    
    varNames = {'Run_ID', 'DIV', 'Assay', 'Well', 'NeuronType', 'Time', 'Chip_ID', ...
                'mean_IBI', 'cov_IBI', 'mean_Burst_Peak', 'cov_Burst_Peak', ...
                'Number_Bursts', 'BurstRate', 'mean_Spike_per_Burst', 'cov_Spike_per_Burst', ...
                'mean_Burst_Peak_Abs', 'cov_Burst_Peak_Abs', 'mean_BurstDuration', 'cov_BurstDuration', 'Baseline', ...
                'IBI_List', 'Burst_Peak_List', 'Abs_Burst_Peak_List', 'Burst_Times_List', 'SpikesPerBurst_List'};
    
    emptyRow = table(run_id, div, {assay}, well_id, {neuron_type}, timestamp, {chip_id}, ...
                     NaN, NaN, NaN, NaN, ...
                     0, NaN, NaN, NaN, ...
                     NaN, NaN, NaN, NaN, NaN, ...
                     {reason}, {''}, {''}, {''}, {''}, ...
                     'VariableNames', varNames);
end

% Ensure ifelse is available or define it
function out = ifelse(condition, true_val, false_val)
    if condition
        out = true_val;
    else
        out = false_val;
    end
end 