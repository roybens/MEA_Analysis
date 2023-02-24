%This script will go through the output folders from Scope and pull the
%network .h5 files and put them in a single folder.



clear
close all

%% settings
% set DIV 0 date
div0 = '11/09/2022'; % format: MM/DD/YYYY

%set path to folder containing subfolders that contain h5 files
%parentFolderPath = 'C:/Users/Tim/Documents/Silverman/Syngap/MEA/Syngap 2/Raw Data h5 files/testParentFolder/';
parentFolderPath = '/home/jonathan/Documents/Data/B6J_Drugs';
%make output folder
opDir = '/home/jonathan/Documents/Scripts/Matlab/scrpits_output/B6J_Drugs/network_outputs/';

%% defines
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


% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(parentFolderPath, '**/Network/**/*raw.h5'); 
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    %fullFileName = fullfile(theFiles(k).folder, baseFileName);
    pathFileNetwork = fullfile(theFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', pathFileNetwork);
    
    % reset recording info
    scan_runID = 0;
    scan_chipID = 0;
    meanIBI = 0;
    meanBurstPeak = 0;
    nBursts = 0;
    meanSpikesPerBurst = 0;
    hd5Date = 0; 
    scan_div = 0;

    % extract dir information
    fileDirParts = strsplit(pathFileNetwork, filesep); % split dir into elements
    scan_runID = str2double(fileDirParts{end-1}); % extract runID
    scan_runID_text = fileDirParts{end-1};
    scan_chipID = str2double(fileDirParts{end-3}); % extract chipID

    if strcmp(scan_runID_text,'000189')
        continue
    end

    % create fileManager object for the Network recording
    networkData = mxw.fileManager(pathFileNetwork);
   
    % get the startTime of the recordings
    hd5_time = networkData.fileObj.stopTime;
    try
        hd5Date = datetime(hd5_time,'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    catch
        hd5Date = datetime(hd5_time,'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
    end
    scan_div = fix(daysact(div0_date , hd5Date));


    %% Parameter Adjustment
    relativeSpikeTimes = mxw.util.computeRelativeSpikeTimes(networkData);
    % set histogram bin size [s]
    binSize = 0.1;
    % set Gaussian kernel standard deviation [s]
    gaussianSigma = 0.14;
    % set burst detection threshold [rms firing rate]
    thresholdBurst = 1.3;
    % compute Network Activity and detect bursts
    networkAct = mxw.networkActivity.computeNetworkAct(networkData, 'BinSize', binSize,'GaussianSigma', gaussianSigma);
    networkStats = mxw.networkActivity.computeNetworkStats(networkAct, 'Threshold', thresholdBurst);
    

    
    %% Tim's code for averaging and aggregating mean spiking data (IBI, Burst peaks, Spikes within Bursts, # of Bursts etc.)
    %average IBI
    meanIBI = mean(networkStats.maxAmplitudeTimeDiff);
    %average Burst peak (burst firing rate y-value)
    meanBurstPeak = mean(networkStats.maxAmplitudesValues);
    %Number of bursts
    nBursts = length(networkStats.maxAmplitudesTimes);
    %average spikesPerBurst
    % Set the threshold to find the start and stop time of the bursts.
        %This is the Start-Stop threshold from the Scope software
        thresholdStartStop = 0.3;
        % 0.3 means 30% value of the burst peak. Note that by raising
        % the value, the percentage of spikes within bursts and the burst duration 
        % increase, since the bursts are considered wider. 
    
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
        fig = figure('Color', 'w','Position',[0 0 1600 800],'Visible','off');
        %subplot(1,1,1);
        yyaxis right
        mxw.plot.networkActivity(networkAct, 'Threshold', thresholdBurst, 'Figure', false,'Color','#000000')
        %box off;
        %hold on;
        plot(networkStats.maxAmplitudesTimes,networkStats.maxAmplitudesValues,'or');
        ylim([0 5])
        %xlim([0 round(max(relativeSpikeTimes.time)/4)])
        %xlim([0 60])
        hold on;
    
        %subplot(1,1,1);
        yyaxis left
        hold on;
        mxw.plot.rasterPlot(networkData, 'Figure', false, 'Color','#c6c6c6');
        box off;
        %xlim([0 round(max(relativeSpikeTimes.time)/4)])
        xlim([0 120])
        ylim([1 max(relativeSpikeTimes.channel)])
        hold on;
        title("DIV" + scan_div + ' Chip' + scan_chipID)
           
        saveas(gcf,append(opDir,'/network_graphs/',scan_runID_text,'.png'))
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
writetable(T, fullfile(opDir,'Compiled_Networks.csv'));

