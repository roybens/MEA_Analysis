clear
close all

%Read the project names from the files and the DIV start date of the
%project
helperFile= '/home/mmp/Documents/wt_network_sup.xlsx';


% extract runID info from reference excel sheet
T = readtable(helperFile);
ProjectFilePaths = T.(1);
DIV0 = T.(2);
CHIPS = T.(3);
DIVOPT = T.(4);
DIVOPT_formatted = datestr(DIVOPT,'yymmdd');
refFiles = T.(5);
%select the DIV to be scanned. eg : 25 or 27
fileNames= {};
fileFolders ={};

%Iterate through the different projects.
for k = 1 : length(ProjectFilePaths)
    % Iterate over the CHIPID values
    chip_names  = char(CHIPS(k));
    chips = strsplit(chip_names,',');
    div0_date = datetime(DIV0(k), "InputFormat",'MM/dd/yyyy');
    
    % set Gaussian kernel standard deviation [s] (smoothing window)
    gaussianSigma = 0.18; %0.18
    % set histogram bin size [s]
    binSize = 0.3;
    % set minimum peak distance [s]
    minPeakDistance = 0.025;
    % set burst detection threshold [rms / fixed]
    thresholdBurst =1.2; %1.2
    % set fixed threshold;
    use_fixed_threshold = false;
    % Set the threshold to find the start and stop time of the bursts. (start-stop threshold)
    thresholdStartStop = 0.3; %0.3
    % Set Threshold function for later use
    threshold_fn = 'Threshold';
    if use_fixed_threshold
        threshold_fn = 'FixedThreshold';
    end

    
    folderstmp = strsplit(ProjectFilePaths{1},filesep);
    projectName = folderstmp{end-1};
    opDir = '/home/mmp/Documents/script_output/';
    opDir = strcat(opDir,projectName,'/');

    % make output folder
    if not(isfolder(append(opDir,'Network_outputs/Raster_BurstActivity/')))
        mkdir(append(opDir,'Network_outputs/Raster_BurstActivity/'));
    end
    for chip = chips
        % Generate the file pattern with the updated chip ID
        filePattern = fullfile(ProjectFilePaths(k), sprintf('**/%s/%s/Network/**/*raw.h5', DIVOPT_formatted(k,:), char(chip)));
        filePattern = char(filePattern);
        
        % Get the files that match the pattern
        files = dir(filePattern);
        
        % Append the file names to the fileNames cell array
        fileNames = [fileNames;{files.name}'];
        fileFolders = [fileFolders;{files.folder}'];
    end
    refDir = refFiles{k};
    % extract runID info from reference file
    T2 = readtable(refDir);
    folderTable = table();
    folderTable.fileFolders = fileFolders;
    folderTable.fileNames = fileNames;

    %Here I need to find the div opt run ids

    % Extract data where Div based on DIV opt
    divValues = DIVOPT(k);
    filteredTable = T2(ismember(T2.Date, divValues)&ismember(T2.NeuronSource,'WT cortex')&ismember(T2.Assay,'Network today'), :);
    run_ids = unique(filteredTable.(4));
    run_id_and_type = [T2(:,[4,7])];
    
  
    % Filter paths based on the matched numbers
    filteredPathFolders = folderTable(endsWith(fileFolders,cellstr(num2str(run_ids))),:);
        % create table elements
    Run_ID = [];
    DIV = [];
    Time = [];
    Chip_ID = [];
    IBI = [];
    Burst_Peak = [];
    Number_Bursts = [];
    Spike_per_Burst = [];
    Project = [];

    % create a list to catch error runIDs
    error_l = [];
    %for all the paths that needs to be analysed
    for j = 1 : length(filteredPathFolders.fileFolders)
        % reset recording info
        scan_runID = nan;
        scan_chipID = nan;
        meanIBI = nan;
        meanBurstPeak = nan;
        nBursts = nan;
        meanSpikesPerBurst = nan;
        spikesPerBurst = nan;
        hd5Date = nan; 
        scan_div = nan;
    
        baseFileName = filteredPathFolders.fileNames(j);
        pathFileNetwork = fullfile(filteredPathFolders.fileFolders(j), baseFileName);
        % extract dir information
        fileDirParts = strsplit(pathFileNetwork{1}, filesep); % split dir into elements
        scan_runID = str2double(fileDirParts{end-1}); % extract runID
        scan_runID_text = fileDirParts{end-1};
        scan_chipID = str2double(fileDirParts{end-3}); % extract chipID

        if ismember(scan_runID,run_ids)
            fprintf(1, 'Now reading %s\n', pathFileNetwork{1});

        % create fileManager object for the Network recording
        try
            networkData = mxw.fileManager(pathFileNetwork);
        catch
            error_l = [error_l string(scan_runID_text)];
            continue
        end
       
        % get the startTime of the recordings
        hd5_time = networkData.fileObj.stopTime;
        try
            hd5Date = datetime(hd5_time,'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        catch
            hd5Date = datetime(hd5_time,'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
        end
        scan_div = fix(daysact(div0_date , hd5Date));
    

        % compute Network Activity and detect bursts
        relativeSpikeTimes = mxw.util.computeRelativeSpikeTimes(networkData);
        networkAct = mxw.networkActivity.computeNetworkAct(networkData, 'BinSize', binSize,'GaussianSigma', gaussianSigma);
        networkStats = computeNetworkStats_JL(networkAct, threshold_fn, thresholdBurst, 'MinPeakDistance', minPeakDistance);
    
        
        %% Tim's code for averaging and aggregating mean spiking data (IBI, Burst peaks, Spikes within Bursts, # of Bursts etc.)
        %average IBI
        meanIBI = mean(networkStats.maxAmplitudeTimeDiff);
        %average Burst peak (burst firing rate y-value)
        meanBurstPeak = mean(networkStats.maxAmplitudesValues);
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
            end
            meanSpikesPerBurst = mean(spikesPerBurst);
        end
   
        % append information to table elements
        Run_ID = [Run_ID scan_runID];
        DIV = [DIV scan_div];
        Time = [Time hd5Date];
        Chip_ID = [Chip_ID scan_chipID];
        IBI = [IBI meanIBI];
        Burst_Peak = [Burst_Peak meanBurstPeak];
        Number_Bursts = [Number_Bursts nBursts];
        Spike_per_Burst = [Spike_per_Burst meanSpikesPerBurst];
        Project = [Project projectName];
        runIDstemp = run_id_and_type.Run_;
        types = run_id_and_type.NeuronSource;
        index = find(runIDstemp == scan_runID);
        targetType = types{index};


        % plot results
            figure('Color','w','Position',[0 0 400 800],'Visible','off');
            subplot(2,1,1);
            mxw.plot.rasterPlot(networkData,'Figure',false);
            box off;
            %xlim([0 round(max(relativeSpikeTimes.time)/4)])
            xlim([0 120])
            ylim([1 max(relativeSpikeTimes.channel)])
            
            subplot(2,1,2);
            mxw.plot.networkActivity(networkAct,'Threshold',thresholdBurst,'Figure',false);
            box off;
            hold on;
            plot(networkStats.maxAmplitudesTimes,networkStats.maxAmplitudesValues,'or')
            %xlim([0 round(max(relativeSpikeTimes.time)/4)])
            xlim([0 120])
            ylim([0 20])
            saveas(gcf,append(opDir,'Network_outputs/Raster_BurstActivity/Raster_BurstActivity',scan_runID_text,'_',num2str(scan_chipID),'_DIV',num2str(scan_div),'_',targetType,'.png'))
    
   
        end


    end


%% construct table
% convert row list to columns
Run_ID = Run_ID';
DIV = DIV';
Time = Time';
Chip_ID = Chip_ID';
IBI = IBI';
Burst_Peak = Burst_Peak';
Number_Bursts = Number_Bursts';
Spike_per_Burst = Spike_per_Burst';
Project = Project';
% make table
T = table(Run_ID,DIV,Time,Chip_ID,IBI,Burst_Peak,Number_Bursts,Spike_per_Burst,Project);
T = sortrows(T,"Run_ID","ascend");
writetable(T, fullfile(opDir,'Network_outputs/Compiled_Networks_Whole.csv'));

if ~isempty(error_l)% create a list to catch error runIDs
    error_l = [];
    error_str = strjoin(error_l,', ');
    fprintf('Unable to read file with runID: %s, file(s) skipped.\n',error_str);
end

fprintf('Network analysis for %s successfully compiled.\n',ProjectFilePaths{k})


end

%CHIP IDS of interest that too from the file.

