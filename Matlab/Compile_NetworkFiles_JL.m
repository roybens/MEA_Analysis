% This script goes through the h5 files in a parent folder,
% plot the Raster and Burst Activities,
% and compile a csv with the 4 mean burst properties over days

clear
close all

% setting starts here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% settings
% set path to folder containing subfolders that contain h5 files
parentFolderPath = '/mnt/harddrive-2/ADNP/';
% set path to excel file that has the reference note
refDir = '/home/jonathan/Documents/Scripts/Python/ADNP_Notes.xlsx';
% set output folder
opDir = '/home/jonathan/Documents/Scripts/Matlab/scripts_output/ADNP/';

% set DIV 0 date
div0 = '05/05/2023'; % format: MM/DD/YYYY

% set Gaussian kernel standard deviation [s] (smoothing window)
gaussianSigma = 0.18;
% set histogram bin size [s]
binSize = 0.1;
% set minimum peak distance [s]
minPeakDistance = 1.0;
% set burst detection threshold [rms firing rate]
thresholdBurst = 0.85;
% set fixed threshold;
use_fixed_threshold = false;
% Set the threshold to find the start and stop time of the bursts. (start-stop threshold)
thresholdStartStop = 0.3;

% setting ends here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set Threshold function for later use
threshold_fn = 'Threshold';
if use_fixed_threshold
    threshold_fn = 'FixedThreshold';
end

% make output folder
if not(isfolder(append(opDir,'Network_outputs/Raster_BurstActivity/')))
    mkdir(append(opDir,'Network_outputs/Raster_BurstActivity/'));
end

% extract runID info from reference excel sheet
T = readtable(refDir);
run_ids = unique(T.(4));


% defines
% convert div 0 to date datatype
div0_date = datetime(div0, "InputFormat",'MM/dd/yyyy');
% create table elements
Run_ID = [];
DIV = [];
Time = [];
Chip_ID = [];
IBI = [];
Burst_Peak = [];
Number_Bursts = [];
Spike_per_Burst = [];

% create a list to catch error runIDs
error_l = [];

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(parentFolderPath, '**/Network/**/*raw.h5'); 
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
    % reset recording info
    scan_runID = 0;
    scan_chipID = 0;
    meanIBI = 0;
    meanBurstPeak = 0;
    nBursts = 0;
    meanSpikesPerBurst = 0;
    hd5Date = 0; 
    scan_div = 0;

    baseFileName = theFiles(k).name;
    pathFileNetwork = fullfile(theFiles(k).folder, baseFileName);
    % extract dir information
    fileDirParts = strsplit(pathFileNetwork, filesep); % split dir into elements
    scan_runID = str2double(fileDirParts{end-1}); % extract runID
    scan_runID_text = fileDirParts{end-1};
    scan_chipID = str2double(fileDirParts{end-3}); % extract chipID

    if ismember(scan_runID,run_ids)
        fprintf(1, 'Now reading %s\n', pathFileNetwork);

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
        networkStats = mxw.networkActivity.computeNetworkStats_JL(networkAct, threshold_fn, thresholdBurst, 'MinPeakDistance', minPeakDistance);
    
        
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
                   meanSpikesPerBurst = mean(spikesPerBurst);
                end
            end
        meanSpikesPerBurst = mean(spikesPerBurst);
    
        % append information to table elements
        Run_ID = [Run_ID scan_runID];
        DIV = [DIV scan_div];
        Time = [Time hd5Date];
        Chip_ID = [Chip_ID scan_chipID];
        IBI = [IBI meanIBI];
        Burst_Peak = [Burst_Peak meanBurstPeak];
        Number_Bursts = [Number_Bursts nBursts];
        Spike_per_Burst = [Spike_per_Burst meanSpikesPerBurst];
      
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

            saveas(gcf,append(opDir,'Network_outputs/Raster_BurstActivity/',scan_runID_text,'.png'))
            %savefig(append(opDir,'Network_outputs/Raster_BurstActivity/',scan_runID_text,'.fig'))
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
% make table
T = table(Run_ID,DIV,Time,Chip_ID,IBI,Burst_Peak,Number_Bursts,Spike_per_Burst);
T = sortrows(T,"Run_ID","ascend");
writetable(T, fullfile(opDir,'Network_outputs/Compiled_Networks.csv'));

if ~isempty(error_l)
    error_str = strjoin(error_l,', ');
    fprintf('Unable to read file with runID: %s, file skipped.\n',error_str);
end

fprintf('Network analysis successfully compiled.\n')
