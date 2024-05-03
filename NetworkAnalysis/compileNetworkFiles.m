function [] = compileNetworkFiles(data)

   
    %ui components.
    fig = data.fig;
    d = uiprogressdlg(fig,'Title','Compiling files',...
        'Message','Start','Cancelable','on');
    drawnow

    tic;
    plotFig=data.plotFig;
    extMetricsFlag = data.extMetricsFlag;
    % Unpack the data structure
    %projectName = data.projectName;
    div0_date =  datetime(data.div0Date, "InputFormat",'MM/dd/yyyy');
    parentFolderPath = data.parentFolderPath;
    refDir = data.refDir;
    opDir = data.opDir;
    gaussianSigma = data.gaussianSigma;
    binSize = data.binSize;
    minPeakDistance = data.minPeakDistance;
    thresholdBurst = data.thresholdBurst;
   
    thresholdStartStop = data.thresholdStartStop;
    % Set Threshold function for later use
    use_fix_threshold = false;
    if use_fix_threshold
    threshold_fn = 'FixedThreshold' ;
    else
    threshold_fn = 'Threshold';
    end
    xlimNetwork = data.xlim;
    ylimNetwork = data.ylim;


    % Ensure output directories exist
    outputFolders = {'Plot60s', 'Plot120s', 'Plot300s', 'Plot600s', 'Figureformat'};
    for i = 1:length(outputFolders)
        folderPath = fullfile(opDir, 'Network_outputs/Raster_BurstActivity', outputFolders{i});
        if ~isfolder(folderPath)
            mkdir(folderPath);
        end
    end


    logFile = data.logFile;

    % extract runID info from reference excel sheet
    refTable = readtable(refDir);
    run_ids = unique(refTable.Run_(strcmp(strtrim(lower(refTable.Assay)), 'network today')));  

    % Initialize an empty table with the correct types for each column,
    
    finalWriteTable = table([], [], [], cell(0,1), [], cell(0,1), [], [], [], [], [], [], [], cell(0,1), ...
        'VariableNames', {
            'Run_ID', 'DIV', 'Well', 'NeuronType', 'Time', 'Chip_ID', ...
            'IBI', 'Burst_Peak', 'Burst_Peak_Normalized', 'Burst_Peak_Abs', ...
            'Number_Bursts', 'Spike_per_Burst', 'BurstDuration', 'ISIString'...
    });
    

    % Get a list of all files in the folder with the desired file name pattern.
    filePattern = fullfile(parentFolderPath, '**/Network/**/*raw.h5'); 
    theFiles = dir(filePattern);
    numFiles = length(theFiles);
    networkResults = cell(numFiles, 1);
 
    
    skippedFiles = cell(numFiles, 1); 
       
    if extMetricsFlag
        extendendMetricsTable = table([], [], [], cell(0,1), [], cell(0,1),[],  ...
        'VariableNames', {
                        'Run_ID', 'DIV', 'Well', 'NeuronType', 'Time', 'Chip_ID','ISIString'...
        });
        extNetworkResults = cell(numFiles,1);
    end

    %parallelizing the files processsing

    % Specify the exact number of cores to use
    numCores = 10;  % Adjust this number based on your needs and resource availability

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
 %profile on;
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
                extFileResults = cell(numWells,1);
            end
            % for loop for processing each well
            parfor z = 1:numWells

                wellID=wellsIDs(z);
                fprintf(1, 'Processing Well %d\n', wellID);
                %tic;
                neuronSourceType = neuronTypes(z);
                try
                    networkData = mxw.fileManager(pathFileNetwork,wellID);
                catch
                    skippedWells{z} = [pathFileNetwork,num2str(wellID)];
                    continue
                end
               
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
                networkStats = computeNetworkStatsNew_JL(networkAct.firingRate,networkAct.time, threshold_fn, thresholdBurst, 'MinPeakDistance', minPeakDistance);
                networkStatsNorm = computeNetworkStatsNew_JL(networkAct.firingRateNorm,networkAct.time, threshold_fn, thresholdBurst, 'MinPeakDistance', minPeakDistance);
                networkStatsAbs = computeNetworkStatsNew_JL(networkAct.absfiringRate,networkAct.time, threshold_fn, thresholdBurst, 'MinPeakDistance', minPeakDistance);
                
                %% Tim's code for averaging and aggregating mean spiking data (IBI, Burst peaks, Spikes within Bursts, # of Bursts etc.)
                %average IBI
                meanIBI = mean(networkStats.maxAmplitudeTimeDiff);
                %average Burst peak (burst firing rate y-value)
                meanBurstPeak = mean(networkStats.maxAmplitudesValues);
                meanBPNorm = mean(networkStatsNorm.maxAmplitudesValues);
                meanAbsBP = mean(networkStatsAbs.maxAmplitudesValues);
                %Number of bursts
                nBursts = length(networkStats.maxAmplitudesTimes);
               
                %average spikesPerBurst        
                if length(networkStats.maxAmplitudesTimes)>3
                    peakAmps = networkStats.maxAmplitudesValues';
                    peakTimes = networkStats.maxAmplitudesTimes;
                    
                    % get the times of the burst start and stop edges
                    edges = double.empty(length(peakAmps),0);
                    for i = 1:length(peakAmps)
                       % take a sizeable (±6 s) chunk of the network activity curve 
                       % around each burst peak point
                       idx = networkAct.time>(peakTimes(i)-6) & networkAct.time<(peakTimes(i)+6);
                       t1 = networkAct.time(idx);
                       a1 = networkAct.firingRate(idx)';
                      
                       % get the arunIDstemp = run_id_and_type(:,1);mplitude at the desired peak width
                       peakWidthAmp = (peakAmps(i)-round(peakAmps(i)*thresholdStartStop));
                       
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
                    ts = ((double(networkData.fileObj.spikes.frameno)...
                        - double(networkData.fileObj.firstFrameNum))/networkData.fileObj.samplingFreq)';
                    ch = networkData.fileObj.spikes.channel;
                    
                    spikesPerBurst = double.empty(length(edges),0);
                    tsWithinBurst = [];
                    chWithinBurst = [];
                    for i = 1:length(edges)
                       idx = (ts>edges(i,1) & ts<edges(i,2));
                       spikesPerBurst(i) = sum(idx); 
                       tsWithinBurst = [tsWithinBurst ts(idx)];
                       chWithinBurst = [chWithinBurst ch(idx)'];
                    end
                    meanSpikesPerBurst = mean(spikesPerBurst);
                    meanBurstDuration = mean(abs(edges(:,1) - edges(:,2)));
                end
                %totalTime = toc;  % Measure the total elapsed time after the loop
    
                 % Display total elapsed time
                %fprintf('in Well %d , Total elapsed time for Basic calc: %f seconds\n', wellID,totalTime);  
                %tic;

                %%  calculating the Interspikeiinterval observed in each channel
                if  extMetricsFlag
                ISIs = accumarray(relativeSpikeTimes.channel, relativeSpikeTimes.time, [], @(x) {diff(x)});
                channelISIStrings = cellfun(@(x) strjoin(arrayfun(@num2str, x, 'UniformOutput', false), ','), ISIs, 'UniformOutput', false);
                combinedISIString = strjoin(channelISIStrings, ';');   
                extFileResults{z} = table(scan_runID, scan_div, wellID, {neuronSourceType{1}}, hd5Date, {scan_chipID}, ...
                        {combinedISIString},...
                       'VariableNames', {
                'Run_ID', 'DIV', 'Well', 'NeuronType', 'Time', 'Chip_ID', ...
                 'ISIString'...
                });
                end
                %%
              % Create a new row for the table, now including combinedISIString
                fileResults{z} = table(scan_runID, scan_div, wellID, {neuronSourceType{1}}, hd5Date, {scan_chipID}, ...
                       meanIBI, meanBurstPeak, meanBPNorm, meanAbsBP, ...
                       nBursts, meanSpikesPerBurst, meanBurstDuration,...
                       'VariableNames', {
                'Run_ID', 'DIV', 'Well', 'NeuronType', 'Time', 'Chip_ID', ...
                'IBI', 'Burst_Peak', 'Burst_Peak_Normalized', 'Burst_Peak_Abs', ...
                'Number_Bursts', 'Spike_per_Burst', 'BurstDuration'...
                });

         
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
                    mxw.plot.networkActivity(networkAct, 'Threshold', thresholdBurst, 'Figure', false);
                    box off;
                    hold on;    
                    plot(networkStats.maxAmplitudesTimes, networkStats.maxAmplitudesValues, 'or');
                    plotFileName =  sprintf('Raster_BurstActivity_%s_WellID_%d_%s_DIV%d_%s', scan_runID_text, wellID, scan_chipID, scan_div, strrep(neuronSourceType{1}, ' ', ''));
                    fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', 'Figureformat',plotFileName);
                    savefig(f, [fileNameBase '.fig']);
                    
                    ylim([0 ylimNetwork]);
                    xlim([0 xlimNetwork])
                    fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', 'Plot60s',plotFileName);
                           
                           
                    % Save in different formats
                    print(f, [fileNameBase '.png'], '-dpng', '-r300');
                    
                    % Save the figure in EPS format
                    print(f, [fileNameBase '.eps'], '-depsc', '-r300');
                    

                    %totalTime = toc;  % Measure the total elapsed time after the loop
    
                 % Display total elapsed time
                    %fprintf('in Well %d Total elapsed time for plot and save for 60s: %f seconds\n', wellID,totalTime); 
                    subplot(2,1,1);
                    xlim([0 120])
                    ylim([1 max(relativeSpikeTimes.channel)])
                    subplot(2,1,2);
                    xlim([0 120])
                    ylim([0 ylimNetwork])
                    fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', 'Plot120s',plotFileName); 
                  
                    % Save in different formats
                    print(f, [fileNameBase '.png'], '-dpng', '-r300');
                    
                    % Save the figure in EPS format
                    print(f, [fileNameBase '.eps'], '-depsc', '-r300');

    
                 
                    subplot(2,1,1);
                    xlim([0 300])
                    ylim([0 max(relativeSpikeTimes.channel)])
                    subplot(2,1,2);
                    xlim([0 300])
                    ylim([0 ylimNetwork])
                     % Assuming you want to display metrics above the raster plot
                    textString = sprintf('# Bursts: %d | Mean Burst Duration: %.2fs | mean SpB: %.2f | mean IBI: %.2fs | Mean BP: %.2fs', nBursts,meanBurstDuration, meanSpikesPerBurst, meanIBI, meanBurstPeak);
                    
                    % Create one annotation box containing all the text entries
                    % Adjust the position vector [x y width height] as needed
                    annotation('textbox', [0.1, 0.425, 0.9, 0.1], 'String', textString, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
                    fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', 'Plot300s',plotFileName);

                    % Save in different formats
                    print(f, [fileNameBase '.png'], '-dpng', '-r300');
                    
                    % Save the figure in EPS format
                    print(f, [fileNameBase '.eps'], '-depsc', '-r300');
                    
    
    
                    subplot(2,1,1);
                    xlim([0 600])
                    ylim([0 max(relativeSpikeTimes.channel)])
                    subplot(2,1,2);
                    xlim([0 600])
                    ylim([0 ylimNetwork])
                    fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', 'Plot600s',plotFileName);
                    
                        % Save the figure in PNG format
                    print(f, [fileNameBase '.png'], '-dpng', '-r300');
                    
                    % Save the figure in EPS format
                    print(f, [fileNameBase '.eps'], '-depsc', '-r300');
                    

                    close(f);
                    %totalTime = toc;  % Measure the total elapsed time after the loop
    
                 % Display total elapsed time
                     %fprintf('in Well %d Total elapsed time for plot and save all graphs: %f seconds\n', wellID,totalTime); 
            end

            end

        end
       networkResults{k} = vertcat(fileResults{~cellfun(@isempty, fileResults)});
       skippedFiles{k} = vertcat(skippedWells{~cellfun(@isempty, skippedWells)});
       if extMetricsFlag
           extNetworkResults{k} = vertcat(extFileResults{~cellfun(@isempty, extFileResults)});
       end
    end
    delete(gcp('nocreate'));
   % profile viewer;
    
    skippedFiles = vertcat(skippedFiles{~cellfun(@isempty, skippedFiles )});
    % Concatenate all tables from each primary file into the final table
    finalWriteTable = vertcat(networkResults{~cellfun(@isempty, networkResults)});
    
    writetable(finalWriteTable, fullfile(opDir,'Network_outputs/Compiled_Networks.csv'));
    
    if extMetricsFlag

        extendendMetricsTable =vertcat(extNetworkResults{~cellfun(@isempty, extNetworkResults)});
        writetable(extendendMetricsTable, fullfile(opDir,'Network_outputs/Extended_Metrics.csv'));
    end

    if ~isempty(skippedFiles)
        
        fprintf(1,'Unable to read file with runID: %s, file(s) skipped.\n',skippedFiles);
        %fprintf(1,'Unable to read file with runID: %s, file(s) skipped.\n',skippedFiles);
    end
    %printf(logFile,'Network analysis successfully compiled.\n');
    fprintf(' Total elapsed time for execution: %f seconds\n', toc);
    fprintf(1,'Network analysis successfully compiled.\n');

    end





