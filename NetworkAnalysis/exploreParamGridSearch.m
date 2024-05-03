

gridSearchWrapper();

function gridSearchWrapper()
    % Setup the parallel pool with a predefined number of workers (e.g., 10)
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

    % Define directories
    fileDir ='/mnt/disk20tb/PrimaryNeuronData/Maxtwo/SHANK3_1/SHANK3_1/231205/';
    refDir = '/mnt/disk15tb/paula/Main_DA_Projects/Ref_Files/SHANK3_1 data/SHANK3_1_ref.xlsx';
    outDir = '/mnt/disk15tb/paula/Main_DA_Projects/data_analysis_output/Primary Neurons/shank3_test/';

    % Define parameter ranges in the form [start, step, end]
    paramRanges = {
        [0.1, 0.02, 2],  ... Range for GaussianSigma
        [0.01, 0.01, 0.2], ...  Range for BinSize
        [0.8, 0.1, 1.5],   ... Range for Threshold
        [0.5, 0.5, 2.0] ,  ... Range for MinPeakDistance
        [0.2 0.1 0.7] ... Range for threshold start-stop
    };

    % Execute the grid search
    Compare_NetworkParametersFunction(fileDir, refDir, outDir, paramRanges);
end

function Compare_NetworkParametersFunction(fileDir, refDir, outDir, paramRanges)
    % Prepare output directory
    opDir = fullfile(outDir, 'ParameterComparisons_Allcomb/');
    if ~exist(opDir, 'dir')
        mkdir(opDir);
    end

    % Generate all combinations of parameter values
    [gaussians, bins, thresholds, peaks,thresholdStartStop] = ndgrid(paramRanges{1}(1):paramRanges{1}(2):paramRanges{1}(3), ...
                                                  paramRanges{2}(1):paramRanges{2}(2):paramRanges{2}(3), ...
                                                  paramRanges{3}(1):paramRanges{3}(2):paramRanges{3}(3), ...
                                                  paramRanges{4}(1):paramRanges{4}(2):paramRanges{4}(3),...
                                                  paramRanges{5}(1):paramRanges{5}(2):paramRanges{5}(3));
    allCombos = [gaussians(:) bins(:) thresholds(:) peaks(:) thresholdStartStop(:)];
      % Parallel processing over all parameter combinations
    numCombos = size(allCombos, 1);
    resultsArray = cell(numCombos+1, 1);  % Preallocate cell array for results
    resultsArray{1} = {'ChipID','WellID','NeuronType','gaussian','binSize','Threshold','MinPeakDist','ThresholdStartStop','IBI','Burst Peak','# Bursts','Spikes per Burst'};
  
    parfor idx = 1:numCombos
        gaussianSigma = allCombos(idx, 1);
        binSize = allCombos(idx, 2);
        threshold = allCombos(idx, 3);
        minPeakDistance = allCombos(idx, 4);
        thresholdStartStop = allCombos(idx,5)
        try
            % Simulate analysis function
            
            resultsArray{idx+1} = simulateAnalysis(fileDir,refDir ,gaussianSigma, binSize, threshold, minPeakDistance,thresholdStartStop);
            % Save results functionally, consider using a more detailed  save function to handle large dat
        catch ME
            fprintf('Error processing index %d: %s\n', idx, ME.message);
            % Optionally log errors to a file or error handling structure
        end

    end
    % Combine results into a single table
    allResults = vertcat(resultsArray{:});
    writetable(allResults, fullfile(opDir, 'combined_results.csv'));
    delete(gcp('nocreate'))
    print("Exploration done")
end

function Averages_chipwise = simulateAnalysis(parentFolderPath,refDir, gaussianSigma, binSize, threshold, minPeakDistance,thresholdStartStop)
    % Placeholder for actual analysis function
    fprintf('Running analysis with Gaussian=%f, BinSize=%f, Threshold=%f, MinPeak=%f\n', ...
            gaussianSigma, binSize, threshold, minPeakDistance);
    % Set Threshold function for later use
    use_fixed_threshold =0;
    threshold_fn = 'Threshold';
    if use_fixed_threshold
    threshold_fn = 'FixedThreshold';
    end
    T = readtable(refDir);
        % create a list to catch error runIDs
    error_l = [];
    
    Assay ="today";
    % extract run ids based on the desired assay type
    %to do: check if only for network today/best it needs to be done.
    assay_T = T(contains(T.(3),'network today/best',IgnoreCase=true) & contains(T.(3),Assay, IgnoreCase=true),:);
    asssy_runIDs = unique(assay_T.(4)).';
    
    % Get a list of all files in the folder with the desired file name pattern.
    filePattern = fullfile(parentFolderPath, '**/Network/**/*raw.h5'); 
    theFiles = dir(filePattern);
    Averages_chipwise ={};
    for f = 1 : length(theFiles)
        baseFileName = theFiles(f).name;
        %fullFileName = fullfile(theFiles(k).folder, baseFileName);
        pathFileNetwork = fullfile(theFiles(f).folder, baseFileName);
        % extract dir information
        fileDirParts = strsplit(pathFileNetwork, filesep); % split dir into elements
        runID = str2double(fileDirParts{end-1}); % extract runID
        if ismember(runID,asssy_runIDs)
           
            fprintf(1, 'Now reading %s\n', pathFileNetwork);
            idx = T.Run_ == runID;
    
            wellsIDs = T.Wells_Recorded(idx);
            if iscell(wellsIDs)
            if ismember(',', wellsIDs{1})
            wellsIDs = strsplit(wellsIDs{1}, ',');
            wellsIDs = cellfun(@str2double,wellsIDs);
            else
                error('wellsIDs are not comma separated correclty');
            end
            end
         
            neuronTypes = T.NeuronSource(idx);
            
            if ismember(',', neuronTypes{1})
            neuronTypes = strsplit(neuronTypes{1}, ',');
            end
            Averages_wellwise ={};
            for z = 1:length(wellsIDs)
            
            wellID=wellsIDs(z);
            fprintf(1, 'Processing Well %d\n', wellID);
            %fprintf(logFile, 'Processing Well %d\n', wellID);
            neuronSourceType = neuronTypes(z);
            % create fileManager object for the Network recording
            try
                networkData = mxw.fileManager(pathFileNetwork,wellID);
            catch
                error_l = [error_l string(runID)];
                continue
            end
            

          
        

    
            % reset values
            meanSpikesPerBurst = nan;
            meanIBI = nan;        
            meanBurstPeak = nan;       
            nBursts = nan;
            spikesPerBurst = NaN;
            networkAct_opt = mxw.networkActivity.computeNetworkAct(networkData, 'BinSize', binSize,'GaussianSigma', gaussianSigma);
            networkStats_opt = computeNetworkStats_JL(networkAct_opt, threshold_fn, threshold, 'MinPeakDistance', minPeakDistance);

    
            meanIBI_opt = mean(networkStats_opt.maxAmplitudeTimeDiff);
            meanBurstPeak_opt = mean(networkStats_opt.maxAmplitudesValues);
            nBursts_opt = length(networkStats_opt.maxAmplitudesTimes);
            
            %{
            % Set the threshold to find the start and stop time of the bursts.
            %This is the Start-Stop threshold from the Scope software
            %thresholdStartStop = 0.3;
            % 0.3 means 30% value of the burst peak. Note that by raising
            % the value, the percentage of spikes within bursts and the burst duration 
            % increase, since the bursts are considered wider. 
            %}
    
            if length(networkStats_opt.maxAmplitudesTimes)>3
                peakAmps = networkStats_opt.maxAmplitudesValues';
                peakTimes = networkStats_opt.maxAmplitudesTimes;
                
                % get the times of the burst start and stop edges
                edges = double.empty(length(peakAmps),0);
                for i = 1:length(peakAmps)
               % take a sizeable (±6 s) chunk of the network activity curve 
               % around each burst peak point
                   idx = networkAct_opt.time>(peakTimes(i)-6) & networkAct_opt.time<(peakTimes(i)+6);
                   t1 = networkAct_opt.time(idx);
                   a1 = networkAct_opt.firingRate(idx)';
                  
                   % get the amplitude at the desired peak width
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
                   meanSpikesPerBurst = mean(spikesPerBurst);
                end
                end
                
            %%Tim's code for averaging and aggregating mean spiking data (IBI, Burst
            %%peaks, Spikes within bursts, # of Bursts etc.)
                
            %chipID =  str2num( regexprep( pathFileNetwork, {'\D*([\d\.]+\d)[^\d]*', '[^\d\.]*'}, {'$1 ', ' '} ) )
            extractChipIDWellID = regexp(pathFileNetwork,'\d{5}\.?','match');
            %chipID = regexp(pathFileNetwork,'\d{5}\.?\d*','match')
            chipID = extractChipIDWellID(:,2);
            
            
            chipAverages = [];
            
            %average spikesPerBurst
            meanSpikesPerBurst = mean(spikesPerBurst);
            
            %average IBI
            meanIBI = mean(networkStats_opt.maxAmplitudeTimeDiff);
            
            %average Burst peak (burst firing rate y-value)
            meanBurstPeak = mean(networkStats_opt.maxAmplitudesValues);
            
            %Number of bursts
            nBursts = length(networkStats_opt.maxAmplitudesTimes);
            
            chipAverages = [meanSpikesPerBurst, meanIBI, meanBurstPeak, nBursts];
            Averages_wellwise = [Averages_wellwise; {chipID,wellID,strrep(neuronSourceType{1},' ',''),binSize,threshold,minPeakDistance,thresholdStartStop,meanIBI_opt, meanBurstPeak_opt, nBursts_opt, meanSpikesPerBurst}];
               
            
            end
            % T1 = cell2table(Averages_opt(2:end,:),'VariableNames',Averages_opt(1,:));
            % IDstring = string(chipID);
            % WellString = string(wellID);
            % writetable(T1,opDir + IDstring +WellString+ '.csv');
        end
        Averages_chipwise =[Averages_chipwise;Averages_wellwise];
        
    end

    end

 
function saveResults(results, opDir, idx)
    % Function to save results
    filename = sprintf('results_%d.mat', idx);
    save(fullfile(opDir, filename), 'results');
end


