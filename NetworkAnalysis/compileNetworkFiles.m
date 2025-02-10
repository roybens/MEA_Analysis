function [] = compileNetworkFiles(data)

   
    %ui components.
    fig = data.fig;
    fprintf(1,"inside the compilenetwork")
    d = uiprogressdlg(fig,'Title','Compiling files',...
        'Message','Start','Cancelable','on');
    drawnow
    tic;
    plotFig=data.plotFig;
    extMetricsFlag = data.extMetricsFlag;
    epsPlot= false
    % Unpack the data structure
   
    div0_date =  datetime(data.div0Date, "InputFormat",'MM/dd/yyyy');
    parentFolderPath = data.parentFolderPath;
    refDir = data.refDir;
    opDir = data.opDir;
    gaussianSigma = data.gaussianSigma;
    binSize = data.binSize;
    minPeakDistance = data.minPeakDistance;
    minProminence = data.minProminience ;
    thresholdBurst = data.thresholdBurst;
   
    thresholdStartStop = data.thresholdStartStop;
    thresholdFunction = data.thresholdMethod;
   
    xlimNetwork = data.xlim;
    ylimNetwork = data.ylim;


    % Ensure output directories exist

    
    % outputFolders = {'Plot60s', 'Plot120s', 'Plot300s', 'Plot600s', 'Figureformat'};
    % for i = 1:length(outputFolders)
    %     folderPath = fullfile(opDir, 'Network_outputs/Raster_BurstActivity', outputFolders{i});
    %     if ~isfolder(folderPath)
    %         mkdir(folderPath);
    %     end
    % end
    
    

    logFile = data.logFile;

    % extract runID info from reference excel sheet
    refTable = readtable(refDir);
    % Convert the 'Assay' column to lowercase and trim any leading/trailing whitespace
    assayColumn = strtrim(lower(refTable.Assay));
    
    % Find rows that contain 'network today' or 'network'
    containsNetworkToday = contains(assayColumn, 'network today');
    containsNetwork = contains(assayColumn, 'network');
    
    % Combine the conditions using logical OR
    combinedCondition = containsNetworkToday | containsNetwork;
    
    % Extract unique run_ids based on the combined condition
    run_ids = unique(refTable.Run_(combinedCondition)); 


    

    % Get a list of all files in the folder with the desired file name pattern.
    filePattern = fullfile(parentFolderPath, '**/Network/**/*raw.h5'); 
    theFiles = dir(filePattern);
    numFiles = length(theFiles);
    networkResults = cell(numFiles, 1);
 
    
    skippedFiles = cell(numFiles, 1); 
      
    allExtMetrics = {};
  

    %%parallelizing the files processsing

    % %%Specify the exact number of cores to use
    numCores = 3;  % Adjust this number based on your needs and resource availability

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


    fprintf(1, 'Plotraster %d , Extended Metrics %d\n ', plotFig,extMetricsFlag);
 
    for k = 1 : numFiles
        if d.CancelRequested
            break
        end
        % Update progress, report current estimate
        d.Value = k/length(theFiles);
        % extract dir informationfileNames
        baseFileName = theFiles(k).name;
        pathFileNetwork = fullfile(theFiles(k).folder, baseFileName);
        
        %this is where the fiel is processed
        fileDirParts = strsplit(pathFileNetwork, filesep); % split dir into elements
        scan_runID = str2double(fileDirParts{end-1}); % extract runID
        scan_runID_text = fileDirParts{end-1};
        scan_chipID = fileDirParts{end-3}; % extract chipID
    
        if ismember(scan_runID,run_ids)
            fprintf(1, 'Now reading %s\n', pathFileNetwork);
            fprintf(logFile, 'Now reading %s\n', pathFileNetwork);

            % create fileManager object for the Network recording
            idx = refTable.Run_ == scan_runID;   % todo each loop takes a copy of this which is an overhead.

            wellsIDs = refTable.Wells_Recorded(idx);
            if iscell(wellsIDs)
            if ismember(',', wellsIDs{1})
                wellsIDs = strsplit(wellsIDs{1}, ',');
                wellsIDs = cellfun(@str2double,wellsIDs);
            else
                error('wellsIDs are not comma separated correctly');
            end
            end
           
            neuronTypes = refTable.NeuronSource(idx);
            if ismember(',', neuronTypes{1})
                neuronTypes = strsplit(neuronTypes{1}, ',');
            end
            
            numWells = length(wellsIDs);
            fileResults = cell(numWells, 1);
            skippedWells = cell(numWells,1);
            if extMetricsFlag
                extData = cell(numWells,1);
            end
            % for loop for processing each well

            sucessFile = false;
            attempt = 0;
            maxRetries = 3;
            while ~sucessFile && attempt < maxRetries
            attempt = attempt + 1;
            try
            for z = 1:numWells
                try
                wellID=wellsIDs(z);
                fprintf(1, 'Processing Well %d\n', wellID);
                %tic;
                neuronSourceType = neuronTypes(z);
                        % Define ChipID_WellID folder
                    folderName = sprintf('%s_Well%d', scan_chipID, wellID); 
                    chipWellFolder = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', folderName);
            
                    % Ensure ChipID_WellID folder exists
                    if ~isfolder(chipWellFolder)
                        mkdir(chipWellFolder);
                    end
            
                    % Define time durations and file formats
                    outputFolders = {'Plot60s', 'Plot120s', 'Plot300s','Plot600s'};
                    formatFolders = {'eps', 'png'};
            
                    % Ensure each duration folder exists inside ChipID_WellID
                    for i = 1:length(outputFolders)
                        durationFolder = fullfile(chipWellFolder, outputFolders{i});
                        if ~isfolder(durationFolder)
                            mkdir(durationFolder);
                        end
            
                        % Create eps and png subfolders
                        for j = 1:length(formatFolders)
                            formatPath = fullfile(durationFolder, formatFolders{j});
                            if ~isfolder(formatPath)
                                mkdir(formatPath);
                            end
                        end
                    end
                networkData = mxw.fileManager(pathFileNetwork,wellID);
                
               
                % get the startTime of the recordings
                hd5_time = networkData.fileObj.stopTime;
                try
                    hd5Date = datetime(hd5_time,'InputFormat', 'yyyy-MM-dd HH:mm:ss');
                catch
                    hd5Date = datetime(hd5_time,'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
                end
                div0_datetime = datetime(div0_date, 'InputFormat', 'yourInputFormatHere');
                scan_div = floor(days(hd5Date - div0_datetime));
                
        
                % compute Network Activity and detect bursts
                relativeSpikeTimes = mxw.util.computeRelativeSpikeTimes(networkData);
                networkAct = mxw.networkActivity.computeNetworkAct(networkData, 'BinSize', binSize,'GaussianSigma', gaussianSigma);
                
                %for 600s second manually reducing to 300
                %networkAct.time = networkAct.time(1:3000);
                %networkAct.firingRate= networkAct.firingRate(1:3000);
                networkStats = computeNetworkStatsModified(networkAct.firingRate,networkAct.time, 'ThresholdMethod',thresholdFunction,'Threshold', thresholdBurst, 'MinPeakProminence', minProminence,'MinPeakDistance', minPeakDistance);
                %networkStatsNorm = computeNetworkStatsNew_JL(networkAct.firingRateNorm,networkAct.time, threshold_fn, thresholdBurst, 'MinPeakDistance', minPeakDistance);
                %networkStatsAbs = computeNetworkStatsModified(networkAct.absfiringRate,networkAct.time, 'ThresholdMethod',thresholdFunction,'Threshold', thresholdBurst, 'MinPeakProminence', minProminence,'MinPeakDistance', minPeakDistance);
                AbsBP = networkAct.absfiringRate(ismember(networkAct.time,networkStats.maxAmplitudesTimes));
                %tally the peaks with normal 

                %% Tim's code for averaging and aggregating mean spiking data (IBI, Burst peaks, Spikes within Bursts, # of Bursts etc.)
                %average IBI
                meanIBI = mean(networkStats.maxAmplitudeTimeDiff);
                %average Burst peak (burst firing rate y-value)
                meanBurstPeak = mean(networkStats.maxAmplitudesValues);
                %meanBPNorm = mean(networkStatsNorm.maxAmplitudesValues);
                meanAbsBP = mean(AbsBP);
                %Number of bursts
                nBursts = length(networkStats.maxAmplitudesTimes);

                %data length time
                recordingTime= networkData.fileObj.dataLenTime;

                burstRate =  nBursts / recordingTime;



               % Calculate standard deviations
                stdIBI = std(networkStats.maxAmplitudeTimeDiff);
                stdBurstPeak = std(networkStats.maxAmplitudesValues);
                %stdBPNorm = std(networkStatsNorm.maxAmplitudesValues);
                stdAbsBP = std(AbsBP);
                % Calculate coefficients of variation
                covIBI = (stdIBI / meanIBI) * 100;  % CV as a percentage
                covBurstPeak = (stdBurstPeak / meanBurstPeak) * 100;  % CV as a percentage
                %covBPNorm = (stdBPNorm / meanBPNorm) * 100;  % CV as a percentage
                covAbsBP = (stdAbsBP/meanAbsBP)* 100;

                
                ts = ((double(networkData.fileObj.spikes.frameno)...
                    - double(networkData.fileObj.firstFrameNum))/networkData.fileObj.samplingFreq)';
                
                ch = networkData.fileObj.spikes.channel;
                ch = ch(ts>0);
                ts = ts(ts>0);
                %%firing rate
              % Loop through each bin size
                   
                frbinSizes = [0.01,0.1,1,10];
            
                for frbinSize = frbinSizes

                        % Ensure output directories exist for this bin size
                        % Define directory paths
                        firingRateMatrixFolder = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', sprintf('FiringRateMatrix_BinSize_%2f', frbinSize));
                        heatmapFolder = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', sprintf('Heatmaps_BinSize_%2f', frbinSize));
                        svdFolder = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', sprintf('SVD_BinSize_%2f', frbinSize));
                        logheatmapFolder = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', sprintf('Log_Heatmaps_BinSize_%2f', frbinSize));
                        % Create directories if they don't exist
                        if ~isfolder(firingRateMatrixFolder)
                            mkdir(firingRateMatrixFolder);
                        end
                        if ~isfolder(logheatmapFolder)
                            mkdir(logheatmapFolder);
                        end
                        if ~isfolder(heatmapFolder)
                            mkdir(heatmapFolder);
                        end
                        if ~isfolder(svdFolder)
                            mkdir(svdFolder);
                        end
                        % Compute firing rate matrix
                        binEdges = 0:frbinSize:networkData.fileObj.dataLenTime;
                        firingRateMatrix = zeros(max(ch), length(binEdges) - 1);
                        for chIdx = 1:max(ch)
                            spikeTimes = ts(ch == chIdx);
                            firingRateMatrix(chIdx, :) = histcounts(spikeTimes, binEdges);
                        end
                                                %%SVD computation
                        [U, S, V] = svd(firingRateMatrix, 'econ');

                        % Extract the singular values from the diagonal matrix S
                        singularValues = diag(S);
                    
                        % Number of singular vectors (columns) to visualize
                        %numSingularValues = length(singularValues);
                        numSingularValues = 20;
                    
                        % Create a single figure for combined visualization
                        figure('Visible', 'off');
                    
                        % Subplot 1: Singular values as a line plot
                        subplot(3, 1, 1);
                        plot(singularValues, '-o', 'LineWidth', 1.5);
                        xlabel('Component number');
                        ylabel('Singular value (variance explained)');
                        title(sprintf('Singular values from SVD (Bin Size: %f)', frbinSize));
                        grid on;
                    
                        % Subplot 2: Heatmap of all left singular vectors (U)
                        subplot(3, 1, 2);
                        imagesc(U(:, 1:numSingularValues));  % Use all columns of U
                        colorbar;
                        xlabel('Singular vector index');
                        ylabel('Channel');
                        title('All left singular vectors (U)');
                    
                        % Subplot 3: Heatmap of all right singular vectors (V)
                        subplot(3, 1, 3);
                        imagesc(V(:, 1:numSingularValues));  % Use all columns of V
                        colorbar;
                        xlabel('Singular vector index');
                        ylabel('Time bin');
                        title('All right singular vectors (V)');
                    
                        % Save the combined figure
                        combinedPlotPath = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity',...
                            sprintf('SVD_BinSize_%2f', frbinSize),...
                            sprintf('CombinedSVD_RunID%d_%s_WellID%d_DIV%d.png', scan_runID,scan_chipID ,wellID, scan_div));
                        saveas(gcf, combinedPlotPath);
                        close(gcf);
                        % % Save firing rate matrix
                        % firingRateMatrixPath = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', ...
                        %     sprintf('FiringRateMatrix_BinSize_%.2f', frbinSize), ...
                        %     sprintf('FiringRateMatrix_RunID%d_WellID%d.mat', scan_runID, wellID));
                        % save(firingRateMatrixPath, 'firingRateMatrix');

                        % Generate and save heatmaps for 60s, 120s, and 300s


                        timeIntervals = [60, 120, 300];
                        for tIdx = 1:length(timeIntervals)
                            tLimit = timeIntervals(tIdx);
                            tBinIdx = find(binEdges <= tLimit, 1, 'last')-1;
                            figure('Visible', 'off');
                            imagesc(firingRateMatrix(:, 1:tBinIdx));
                            title(sprintf('Firing Rate Heatmap (First %ds, Bin Size %2f)', tLimit, frbinSize));
                            xlabel('Time Bins');
                            ylabel('Channel');
                            colorbar;
                            colormap('hot');
                            heatmapPath = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', ...
                                sprintf('Heatmaps_BinSize_%2f', frbinSize), ...
                                sprintf('Heatmap_RunID%d_%s_WellID%d_DIV%d_%ds.png', scan_runID,scan_chipID ,wellID, scan_div,tLimit));
                            saveas(gcf, heatmapPath);
                            close(gcf);
                        
                            % Log-scaled heatmap
                            figure('Visible', 'off');
                            % Add 1 to avoid log(0) issues
                            imagesc(log10(firingRateMatrix(:, 1:tBinIdx) + 1));
                            title(sprintf('Log-Scaled Firing Rate Heatmap (First %ds, Bin Size %2f)', tLimit, frbinSize));
                            xlabel('Time Bins');
                            ylabel('Channel');
                            colorbar;
                            colormap('hot');
                            logHeatmapPath = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', ...
                                sprintf('Log_Heatmaps_BinSize_%2f', frbinSize), ...
                                sprintf('Log_Heatmap_RunID%d_%s_WellID%d_DIV%d_%ds.png', scan_runID,scan_chipID, wellID, scan_div, tLimit));
                            saveas(gcf, logHeatmapPath);
                            close(gcf);
                        end

                  end

                %average spikesPerBurst        
                if length(networkStats.maxAmplitudesTimes)>3
                    peakAmps = networkStats.maxAmplitudesValues';
                    peakTimes = networkStats.maxAmplitudesTimes;
                    
                    % get the times of the burst start and stop edges
                    edges = double.empty(length(peakAmps),0);
                    for i = 1:length(peakAmps)
                       % take a sizeable (�6 s) chunk of the network activity curve 
                       % around each burst peak point
                       %MANDAR I THink using this there is a p
                       idx = networkAct.time>(peakTimes(i)-6) & networkAct.time<(peakTimes(i)+6);
                       t1 = networkAct.time(idx);
                       a1 = networkAct.firingRate(idx)';
                      
                       % get the arunIDstemp = run_id_and_type(:,1);mplitude at the desired peak width
                       peakWidthAmp = (peakAmps(i)-peakAmps(i)*thresholdStartStop);
                        %the fraction of the total peak height measured
                        %from the top at hich the burst start stop are
                        %defined.
                       
                       % get the indices of the peak edges
                       idx1 = find(a1<peakWidthAmp & t1<peakTimes(i));
                       idx2 = find(a1<peakWidthAmp & t1>peakTimes(i));
                       
                       if ~isempty(idx1)&&~isempty(idx2)       
                           tBefore = t1(idx1(end));
                           tAfter = t1(idx2(1));
                           edges(i,[1 2]) = [tBefore tAfter];
                       end
                    end
                    
                   % identify spikes that fall within the bursts


                    if isempty(edges)

                    spikesPerBurst = NaN;
                    tsWithinBurst = NaN;
                    chWithinBurst = NaN;
                    meanSpikesPerBurst=NaN;
                    meanBurstDuration =NaN;
                    covSpikesPerBurst =NaN ;
                    covBurstDuration = NaN;
                    IBIStrings = {''};
                    BurstPeakStrings = {''};
                    AbsBurstPeakStrings = {''};
                    BurstDurationStrings = {''};
                    spikesPerBurstStrings = {''};

                    else
                    % Initialize the spikesPerBurst array to zeros 
                    spikesPerBurst = zeros(length(edges), 1);
                    
                    % Create an array where each element is the index of the timestamp
                    % This is used for broadcasting comparisons against edges
                    tsIndices = repmat(1:length(ts), size(edges, 1), 1);
                    
                    % Create matrices of start and end times for easy comparison
                    startTimes = repmat(edges(:,1), 1, length(ts));
                    endTimes = repmat(edges(:,2), 1, length(ts));
                    
                    % Logical matrices where true indicates a timestamp within the burst interval
                    inBurstWindows = (ts > startTimes) & (ts < endTimes);
                    
                    % Sum each row to get the count of true values (spikes within each burst)
                    spikesPerBurst = sum(inBurstWindows, 2);

           
                    
                    
                    % For tsWithinBurst and chWithinBurst, you can use logical indexing
                    tsWithinBurst = ts(any(inBurstWindows, 1));
                    chWithinBurst = ch(any(inBurstWindows, 1));
                    tsOutsideBurst = ts(~any(inBurstWindows, 1));
                    chOutsideBurst = ch(~any(inBurstWindows, 1));
                    meanSpikesPerBurst = mean(spikesPerBurst);
                    stdSpikesPerBurst = std(spikesPerBurst);
                    bursts = abs(edges(:,1) - edges(:,2));
                    meanBurstDuration = mean(bursts);
                    stdBurstDuration = std(meanBurstDuration);
                    covSpikesPerBurst = (stdSpikesPerBurst/meanSpikesPerBurst)*100;
                    covBurstDuration = (stdBurstDuration/meanBurstDuration) * 100;
                    if isempty(edges)
                        % If edges are empty, use all time points
                        baselineFiringRate = mean(networkAct.firingrate);
                    else
                        % Otherwise, exclude the times within the burst intervals
                        excludeIdx = false(size(networkAct.time));
                    
                        for i = 1:size(edges, 1)
                            excludeIdx = excludeIdx | (networkAct.time >= edges(i,1) & networkAct.time <= edges(i,2));
                        end
                    
                        includeIdx = ~excludeIdx;
                        baselineFiringRate = mean(networkAct.firingrate(includeIdx));
                    end
                                        
                    % Convert spikesPerBurst to a comma-separated string
                    IBIStrings = {strjoin(arrayfun(@num2str, networkStats.maxAmplitudeTimeDiff, 'UniformOutput', false), ',')};
                    BurstPeakStrings = {strjoin(arrayfun(@num2str, networkStats.maxAmplitudesValues, 'UniformOutput', false), ',')};
                    AbsBurstPeakStrings = {strjoin(arrayfun(@num2str, AbsBP, 'UniformOutput', false), ',')};
                    BurstDurationStrings = {strjoin(arrayfun(@num2str, bursts, 'UniformOutput', false), ',')};
                    spikesPerBurstStrings = {strjoin(arrayfun(@num2str, spikesPerBurst, 'UniformOutput', false), ',')};
                                        

                    end
                else
                    spikesPerBurst = NaN;
                    tsWithinBurst = NaN;
                    chWithinBurst = NaN;
                    meanSpikesPerBurst=NaN;
                    meanBurstDuration =NaN;
                    covSpikesPerBurst =NaN ;
                    covBurstDuration = NaN;
                    IBIStrings = {''};
                    BurstPeakStrings = {''};
                    AbsBurstPeakStrings = {''};
                    BurstDurationStrings = {''};
                    spikesPerBurstStrings = {''};

                end
                %totalTime = toc;  % Measure the total elapsed time after the loop
    
                 % Display total elapsed time
                %fprintf('in Well %d , Total elapsed time for Basic calc: %f seconds\n', wellID,totalTime);  
                %tic;

                %%  calculating the Interspikeiinterval observed in each channel
                if  extMetricsFlag
                
                if isnan(tsWithinBurst)

                    % Set all variables to empty or default values because tsWithinBurst is empty
                    burstISIs = {};
                    burstISIsCombined = [];
                    meanBurstISI = NaN;
                    IQRABurstISI = NaN;
                    stdBurstISI = NaN;
                    covBurstISI = NaN;
                    instAPFreqWithinBurst =NaN;
                    iqrAPFreqWithinBurst = NaN;
                    rangeAPFreq = NaN;
                    binWidth = NaN;
                    numBins = NaN;
                    binedgesAP =NaN;
                    binAPcountsWithinBurst = NaN;
                    binedgesAPWithinBurst =NaN;

                     meanBurstISIOutside =NaN;
                     covBurstISIOutside = NaN;
                     binAPcountsOutsideBurst =NaN;
                     binedgesAPOutsideBurst =NaN;

                else
                burstISIs = accumarray(chWithinBurst,tsWithinBurst,[],@(x){diff(x)});
                burstISIsCombined = vertcat(burstISIs{:});
                burstISIsCombined = burstISIsCombined(burstISIsCombined<1);
                meanBurstISI = mean(burstISIsCombined);
                IQRABurstISI = iqr(burstISIsCombined);
                stdBurstISI = std(burstISIsCombined);
                covBurstISI = (stdBurstISI/meanBurstISI)*100;
                instAPFreqWithinBurst = (1./burstISIsCombined);
                iqrAPFreqWithinBurst = iqr(instAPFreqWithinBurst);
                rangeAPFreq = max(instAPFreqWithinBurst)- min(instAPFreqWithinBurst);
                binWidth = 2 * iqrAPFreqWithinBurst /(length(IQRABurstISI)^(1/3)) ;%Freedman-Diaconis Rule
                numBins = ceil(rangeAPFreq/binWidth);
                binedgesAP = logspace(log10(min(instAPFreqWithinBurst)),log10(500),numBins);   % hardcoding 500 hz here ,, as anything above is not physiological ( but What about MUA)
              
                [binAPcountsWithinBurst,binedgesAPWithinBurst]= histcounts(instAPFreqWithinBurst,binedgesAP);
                
                % Calculate inter-spike intervals (ISIs) for data outside bursts
                burstISIsOutside = accumarray(chOutsideBurst, tsOutsideBurst, [], @(x) {diff(sort(x))});
                
                % Combine ISIs from all channels
                burstISIsCombinedOutside = vertcat(burstISIsOutside{:});
                burstISIsCombinedOutside = burstISIsCombinedOutside(burstISIsCombinedOutside < 1);  % Filter ISIs less than 1 second
                
                % Calculate mean, IQR, standard deviation, and coefficient of variation for the ISIs
                meanBurstISIOutside = mean(burstISIsCombinedOutside);
                IQRBurstISIOutside = iqr(burstISIsCombinedOutside);
                stdBurstISIOutside = std(burstISIsCombinedOutside);
                covBurstISIOutside = (stdBurstISIOutside / meanBurstISIOutside) * 100;
                
                % Calculate instantaneous AP frequency from ISIs
                instAPFreqOutsideBurst = 1 ./ burstISIsCombinedOutside;
                
                % Compute IQR and range for AP frequencies
                iqrAPFreqOutsideBurst = iqr(instAPFreqOutsideBurst);
                rangeAPFreqOutside = max(instAPFreqOutsideBurst) - min(instAPFreqOutsideBurst);
                
                % Determine histogram bins using the Freedman-Diaconis Rule
                binWidthOutside = 2 * iqrAPFreqOutsideBurst / (length(instAPFreqOutsideBurst)^(1/3));
                numBinsOutside = ceil(rangeAPFreqOutside / binWidthOutside);
                
                % Define bin edges for histogram - assuming physiological limit as 500 Hz
                binedgesAPOutside = logspace(log10(min(instAPFreqOutsideBurst)), log10(500), numBinsOutside);
                
                % Histogram counts using the defined bin edges
                [binAPcountsOutsideBurst, binedgesAPOutsideBurst] = histcounts(instAPFreqOutsideBurst, binedgesAPOutside);
                
                end
                % Combine ISIs from all channels into a single array
                ISIs = accumarray(relativeSpikeTimes.channel, relativeSpikeTimes.time, [], @(x){diff(x)});
                ISIsCombined = vertcat(ISIs{:});
                ISIsCombined = ISIsCombined(ISIsCombined > 0);  % Optional: filter out non-positive values if any
                
                % Calculate mean, standard deviation, and IQR
                meanISI = mean(ISIsCombined);
                stdISI = std(ISIsCombined);
                IQRISI = iqr(ISIsCombined);
                
                % Calculate coefficient of variation (CV)
                covISI = (stdISI / meanISI) * 100;
                
                % Calculate instantaneous frequencies from ISIs
                instFreq = 1 ./ ISIsCombined;
                
                % Compute IQR and range for instantaneous frequencies
                iqrFreq = iqr(instFreq);
                rangeFreq = max(instFreq) - min(instFreq);
                
                % Determine histogram bins using the Freedman-Diaconis Rule
                binWidth = 2 * iqrFreq / (length(instFreq)^(1/3));
                numBins = ceil(rangeFreq / binWidth);
                minFreq = min(instFreq(instFreq > 0));  % Avoid log(0) by ensuring positive frequencies
                
                % Define bin edges for histogram - setting physiological reasonable limits
                binedgesFreq = logspace(log10(minFreq), log10(max(instFreq)), numBins);
                
                % Histogram counts using the defined bin edges
                [binFreqCounts, binedgesFreq] = histcounts(instFreq, binedgesFreq);
                

                %%%now calculating the Fanoi factor
                % ts is the spike times, let me set 100 ms bin with
                binedgesfano = min(ts): 0.1 : max(ts);
                spikeCountsFano = histcounts(ts,binedgesfano);
                meanSpikeCounts = mean(spikeCountsFano);
                varSpikeCounts = var(spikeCountsFano);
                fanoFactor = varSpikeCounts / meanSpikeCounts;


                %channelISIStrings = cellfun(@(x) strjoin(arrayfun(@num2str, x, 'UniformOutput', false), ','), ISIs, 'UniformOutput', false);
                %combinedISIString = strjoin(channelISIStrings, ';');   
            % Add spikesPerBurst to the table
                  fileResults{z} = table(scan_runID, scan_div, wellID, strtrim({neuronSourceType{1}}), hd5Date, {scan_chipID}, ...
                meanIBI, covIBI, meanBurstPeak, covBurstPeak, ...
                nBursts, meanSpikesPerBurst, covSpikesPerBurst, meanAbsBP, covAbsBP, meanBurstDuration, covBurstDuration, ...
                meanISI, covISI, ...
                meanBurstISI, covBurstISI, ...
                meanBurstISIOutside, covBurstISIOutside, ...
                fanoFactor, burstRate,baselineFiringRate, IBIStrings, BurstPeakStrings, AbsBurstPeakStrings, BurstDurationStrings, spikesPerBurstStrings, ...
                'VariableNames', {
                'Run_ID', 'DIV', 'Well', 'NeuronType', 'Time', 'Chip_ID', ...
                'mean_IBI', 'cov_IBI', 'mean_Burst_Peak', 'cov_Burst_Peak', ...
                'Number_Bursts', 'mean_Spike_per_Burst', 'cov_Spike_per_Burst', 'mean_Burst_Peak_Abs', 'cov_Burst_Peak_Abs', 'mean_BurstDuration', 'cov_BurstDuration', ...
                'MeanNetworkISI', 'CoVNetworkISI', ...
                'MeanWithinBurstISI', 'CoVWithinBurstISI', ...
                'MeanOutsideBurstISI', 'CoVOutsideBurstISI', ...
                'Fanofactor', 'BurstRate','Baseline' ...
                'IBI_List', 'Burst_Peak_List', 'Abs_Burst_Peak_List', 'Burst_Times_List', 'SpikesPerBurst_List'
                });

                extMetrics = struct('Run_ID',scan_runID, 'DIV', scan_div,'Well' ,wellID,'Chip_ID',{scan_chipID},'NeuronType',strtrim({neuronSourceType{1}}),...
                    'networkAPFreqBins',binFreqCounts,'networkAPFreqEdges',binedgesFreq,...
                    'burstAPFreqBins',binAPcountsWithinBurst,'burstAPFreqEdges',binedgesAPWithinBurst,...
                    'nonburstAPFreqBins',binAPcountsOutsideBurst,'nonburstAPFreqEdges',binedgesAPOutsideBurst);

                extData{z} = extMetrics;
                else
                %%
              % Create a new row for the table, now including combinedISIString
                fileResults{z} = table(scan_runID, scan_div, wellID, strtrim({neuronSourceType{1}}), hd5Date, {scan_chipID}, ...
                       meanIBI,covIBI, meanBurstPeak,covBurstPeak, meanBPNorm,covBPNorm, meanAbsBP,covAbsBP, ...
                       nBursts, meanSpikesPerBurst,covSpikesPerBurst, meanBurstDuration,covBurstDuration,...
                       'VariableNames', {
                'Run_ID', 'DIV', 'Well', 'NeuronType', 'Time', 'Chip_ID', ...
                'mean_IBI','cov_IBI', 'mean_Burst_Peak','cov_Burst_Peak',...
                'Number_Bursts', 'mean_Spike_per_Burst','cov_Spike_per_Burst','mean_Burst_Peak_Abs','cov_Burst_Peak_Abs', 'mean_BurstDuration','cov_BurstDuration',...             
                });
                end
         
               %totalTime = toc;  % Measure the total elapsed time after the loop
    
                 % Display total elapsed time
                %fprintf('in Well %d Total elapsed time for ISI calc: %f seconds\n', wellID,totalTime); 




                %%
                if plotFig
                 % Define time limits and corresponding directories for saving the plots
                    %tic;
                    % Create figure in the background
                    f = figure('Color','w', 'Position', [0 0 400 800], 'Visible', 'off');
                       
                    % Plot raster and network activity
                    subplot(2,1,1);
                    mxw.plot.rasterPlot(networkData, 'Figure', false);
                    box off;
                    xlim([0 60]);
                    ylim([1 max(relativeSpikeTimes.channel)]);
                
                    subplot(2,1,2);
                    plotNetworkActivityModified(networkAct, 'ThresholdFunction',networkStats.thresholdFunction, 'Threshold',networkStats.threshold, 'Figure', false);
                    box off;
                    hold on;    
                    plot(networkStats.maxAmplitudesTimes, networkStats.maxAmplitudesValues, 'or');
                    plotFileName =  sprintf('Raster_BurstActivity_%s_WellID_%d_%s_DIV%d_%s', scan_runID_text, wellID, scan_chipID, scan_div, strrep(neuronSourceType{1}, ' ', ''));
                    %fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', 'Figureformat',plotFileName);
                    fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity',plotFileName);
                    savefig(f, [fileNameBase '.fig']);
                    
                    ylim([0 ylimNetwork]);
                    xlim([0 xlimNetwork])
                    %rohan made changes here
                    fileNameBase = fullfile(chipWellFolder, 'Plot60s', 'png', plotFileName);
                    %fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', chipWellFolder, 'Plot60s', 'png', plotFileName);
                    print(f, [fileNameBase '.png'], '-dpng', '-r300');
                    if epsPlot
                    %fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', chipWellFolder, 'Plot60s', 'eps', plotFileName);
                    fileNameBase = fullfile(chipWellFolder, 'Plot60s', 'eps', plotFileName);
                    print(f, [fileNameBase '.eps'], '-depsc', '-r300');
                    end

                    %totalTime = toc;  % Measure the total elapsed time after the loop
    
                 % Display total elapsed time
                    %fprintf('in Well %d Total elapsed time for plot and save for 60s: %f seconds\n', wellID,totalTime); 
                    subplot(2,1,1);
                    xlim([0 120])
                    ylim([1 max(relativeSpikeTimes.channel)])
                    subplot(2,1,2);
                    xlim([0 120])
                    ylim([0 ylimNetwork])
                    %fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', 'Plot120s',plotFileName); 
                    %fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', chipWellFolder, 'Plot120s', 'png', plotFileName);
                    fileNameBase = fullfile(chipWellFolder, 'Plot120s', 'png', plotFileName);
                    print(f, [fileNameBase '.png'], '-dpng', '-r300');

                    if epsPlot
                    %fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', chipWellFolder, 'Plot120s', 'eps', plotFileName);
                    fileNameBase = fullfile(chipWellFolder, 'Plot120s', 'eps', plotFileName);
                    print(f, [fileNameBase '.eps'], '-depsc', '-r300');
                    end
                    
                    % Save the figure in EPS format
                   % print(f, [fileNameBase '.eps'], '-depsc', '-r300');


    
                 
                    subplot(2,1,1);
                    xlim([0 300])
                    ylim([0 max(relativeSpikeTimes.channel)])
                    subplot(2,1,2);
                    xlim([0 300])
                    ylim([0 ylimNetwork])
                     % Assuming you want to display metrics above the raster plot
                    textString = sprintf('# Bursts: %d | Mean Burst Duration: %.2fs | mean SpB: %.2f | mean IBI: %.2fs | Mean BP: %.2f', nBursts,meanBurstDuration, meanSpikesPerBurst, meanIBI, meanBurstPeak);
                    timeVector = 0:binSize:max(relativeSpikeTimes.time);   %need to optimize
                    plot(baselineFiringRate*ones(ceil(timeVector(end)),1))

                    % Create one annotation box containing all the text entries
                    % Adjust the position vector [x y width height] as needed
                    annotation('textbox', [0.1, 0.425, 0.9, 0.1], 'String', textString, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
                    %fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', 'Plot300s',plotFileName);

                    %fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', chipWellFolder, 'Plot300s', 'png', plotFileName);
                    fileNameBase = fullfile(chipWellFolder, 'Plot300s', 'png', plotFileName);
                    print(f, [fileNameBase '.png'], '-dpng', '-r300');

                    if epsPlot
                   % fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', chipWellFolder, 'Plot300s', 'eps', plotFileName);
                    fileNameBase = fullfile(chipWellFolder, 'Plot300s', 'eps', plotFileName);
                    print(f, [fileNameBase '.eps'], '-depsc', '-r300');
                    end    

                    
    
    
                    subplot(2,1,1);
                    xlim([0 600])
                    ylim([0 max(relativeSpikeTimes.channel)])
                    subplot(2,1,2);
                    xlim([0 600])
                    ylim([0 ylimNetwork])
                    %fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', 'Plot600s',plotFileName);
                    fileNameBase = fullfile(chipWellFolder, 'Plot600s', 'png', plotFileName);
                        % Save the figure in PNG format
                    print(f, [fileNameBase '.png'], '-dpng', '-r300');
                    
                    if epsPlot
                    fileNameBase = fullfile(chipWellFolder, 'Plot600s', 'eps', plotFileName);
                    print(f, [fileNameBase '.eps'], '-depsc', '-r300');
                    end
                    

                    close(f);
                    %totalTime = toc;  % Measure the total elapsed time after the loop
    
                 % Display total elapsed time
                     %fprintf('in Well %d Total elapsed time for plot and save all graphs: %f seconds\n', wellID,totalTime); 
            end
                catch ME
                fprintf('Error: %s\n', ME.message);
                skippedWells{z} = [pathFileNetwork,num2str(wellID)];
                    continue
                end

            end
               networkResults{k} = vertcat(fileResults{~cellfun(@isempty, fileResults)});
               skippedFiles{k} = vertcat(skippedWells{~cellfun(@isempty, skippedWells)});
               if extMetricsFlag
               allExtMetrics{k} = extData;
               end
               sucessFile = true;
            
            catch ME
                % Handle specific errors
                if contains(ME.message, 'parallel')
                    fprintf('Error: %s\n', ME.message);
                    fprintf('Attempt %d/%d: Restarting parallel pool for file %s.\n', attempt, maxRetries, pathFileNetwork);
                    delete(gcp('nocreate')); % Shut down any existing parallel pool
                    parpool(numCores);    % Restart the parallel pool
                else
                    fprintf('Unexpected error: %s\n', ME.message);
                    break; % Break the loop for unexpected errors
                end
            end
            end
            if ~sucessFile

             fprintf('Failed to process file %s after %d attempts. Skipping this file.\n', pathFileNetwork, maxRetries);
             delete(gcp('nocreate')); % Shut down any existing parallel pool
             parpool(numCores); 
             continue
            end

        end 
    end
    delete(gcp('nocreate'));

   % profile viewer;
    
    skippedFiles = vertcat(skippedFiles{~cellfun(@isempty, skippedFiles )});
    % Concatenate all tables from each primary file into the final table
    finalWriteTable = vertcat(networkResults{~cellfun(@isempty, networkResults)});
    
    csvfilepath =  fullfile(opDir,'Network_outputs/Compiled_Networks.csv');
    if isfile(csvfilepath)
        delete(csvfilepath);
    end
    writetable(finalWriteTable,csvfilepath );
    if extMetricsFlag
       save(fullfile(opDir,'Network_outputs/extendedMetrics.mat'),'allExtMetrics');
    end


    if ~isempty(skippedFiles)
        
        fprintf(1,'Unable to read file with runID: %s, file(s) skipped.\n',skippedFiles);
        %fprintf(1,'Unable to read file with runID: %s, file(s) skipped.\n',skippedFiles);
    end
    %printf(logFile,'Network analysis successfully compiled.\n');
    fprintf(' Total elapsed time for execution: %f seconds\n', toc);
    fprintf(1,'Network analysis successfully comp   iled.\n');

end
