%This script will go through the output folders from Scope and pull the
%network .h5 files and put them in a single folder.



clear
close all


%set path to folder containing subfolders that contain h5 files
%parentFolderPath = 'C:/Users/Tim/Documents/Silverman/Syngap/MEA/Syngap 2/Raw Data h5 files/testParentFolder/';
parentFolderPath = 'D:/MEA Raw Data Backup/Syngap3/230103/';

%make output folder
opDir = char('C:/Users/Tim/Documents/MATLAB/Maxwell Analysis Stuff/MEA_DataOutput_syngap3_DIV21/');
mkdir(opDir);

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(parentFolderPath, '**/Network/**/*raw.h5'); 
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    %fullFileName = fullfile(theFiles(k).folder, baseFileName);
    pathFileNetwork = fullfile(theFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', pathFileNetwork);
    tic
   
    
    % create fileManager object for the Network recording
    networkData = mxw.fileManager(pathFileNetwork);
    
    % %% Parameter Adjustment
    
    relativeSpikeTimes = mxw.util.computeRelativeSpikeTimes(networkData);
    
    % set histogram bin size [s]
    binSize = 0.1;
    
    % set Gaussian kernel standard deviation [s]
    gaussianSigma = 0.1;
    
    % set burst detection threshold [rms firing rate]
    thresholdBurst = 0.4;
    
    % compute Network Activity and detect bursts
    
    networkAct = mxw.networkActivity.computeNetworkAct(networkData, 'BinSize', binSize,'GaussianSigma', gaussianSigma);
    networkStats = mxw.networkActivity.computeNetworkStats(networkAct, 'Threshold', thresholdBurst);
    
    %% Loop through different parameters to find best fitting one. _opt suffix on variable to denote optimization.
    
    
    
    relativeSpikeTimes_opt = mxw.util.computeRelativeSpikeTimes(networkData);
    
    Averages_opt = [{'ChipID','Gaussian','IBI','Burst Peak','# Bursts','Spikes per Burst' }];
    
    for k=.02:.02:1
    % set Gaussian kernel standard deviation [s]
    %This is the smoothing window from Scope software
        gaussianSigma_opt = k;
    % set histogram bin size [s]
        binSize_opt = 0.1;
    % set burst detection threshold [rms firing rate]
        thresholdBurst_opt = 0.4;
    
        networkAct_opt = mxw.networkActivity.computeNetworkAct(networkData, 'BinSize', binSize_opt,'GaussianSigma', gaussianSigma_opt);
        networkStats_opt = mxw.networkActivity.computeNetworkStats(networkAct_opt, 'Threshold', thresholdBurst_opt);
    
        meanIBI_opt = mean(networkStats_opt.maxAmplitudeTimeDiff);
        meanBurstPeak_opt = mean(networkStats_opt.maxAmplitudesValues);
        nBursts_opt = length(networkStats_opt.maxAmplitudesTimes);
    
        % Set the threshold to find the start and stop time of the bursts.
        %This is the Start-Stop threshold from the Scope software
        thresholdStartStop = 0.3;
        % 0.3 means 30% value of the burst peak. Note that by raising
        % the value, the percentage of spikes within bursts and the burst duration 
        % increase, since the bursts are considered wider. 
    
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
    extractChipID = regexp(pathFileNetwork,'\d{5}\.?','match');
    %chipID = regexp(pathFileNetwork,'\d{5}\.?\d*','match')
    chipID = extractChipID(:,2);
    
    
    chipAverages = [];
    
    %average spikesPerBurst
    meanSpikesPerBurst = mean(spikesPerBurst);
    
    %average IBI
    meanIBI = mean(networkStats.maxAmplitudeTimeDiff);
    
    %average Burst peak (burst firing rate y-value)
    meanBurstPeak = mean(networkStats.maxAmplitudesValues);
    
    %Number of bursts
    nBursts = length(networkStats.maxAmplitudesTimes);
    
    chipAverages = [meanSpikesPerBurst, meanIBI, meanBurstPeak, nBursts];
    Averages_opt = [Averages_opt; {chipID,k,meanIBI_opt, meanBurstPeak_opt, nBursts_opt, meanSpikesPerBurst}];
       
    
    end
    T = cell2table(Averages_opt(2:end,:),'VariableNames',Averages_opt(1,:));
    IDstring = string(chipID);
    writetable(T,'C:/Users/Tim/Documents/MATLAB/Maxwell Analysis Stuff/MEA_DataOutput_syngap3_DIV21/' + IDstring + '.csv');
    
    
    %% Tim's code for averaging and aggregating mean spiking data (IBI, Burst peaks, Spikes within Bursts, # of Bursts etc.)
    
    %chipID =  str2num( regexprep( pathFileNetwork, {'\D*([\d\.]+\d)[^\d]*', '[^\d\.]*'}, {'$1 ', ' '} ) )
    extractChipID = regexp(pathFileNetwork,'\d{5}\.?','match');
    %chipID = regexp(pathFileNetwork,'\d{5}\.?\d*','match')
    chipID = extractChipID(:,2);
    
    
    chipAverages = [];
    
    %average spikesPerBurst
    meanSpikesPerBurst = mean(spikesPerBurst);
    
    %average IBI
    meanIBI = mean(networkStats.maxAmplitudeTimeDiff);
    
    %average Burst peak (burst firing rate y-value)
    meanBurstPeak = mean(networkStats.maxAmplitudesValues);
    
    %Number of bursts
    nBursts = length(networkStats.maxAmplitudesTimes);
    
    chipAverages = [meanSpikesPerBurst, meanIBI, meanBurstPeak, nBursts];
   
    %% Make Optimization Plots
    
    relativeSpikeTimes_opt = mxw.util.computeRelativeSpikeTimes(networkData);
    
    for j=.02:.06:.2
        
    % set Gaussian kernel standard deviation [s]
    %This is the smoothing window from Scope software
        gaussianSigma_opt = j;
    
    % set histogram bin size [s]
        binSize_opt = 0.1;
    % set burst detection threshold [rms firing rate]
        thresholdBurst_opt = 0.4;
    
        networkAct_opt = mxw.networkActivity.computeNetworkAct(networkData, 'BinSize', binSize_opt,'GaussianSigma', gaussianSigma_opt);
        networkStats_opt = mxw.networkActivity.computeNetworkStats(networkAct_opt, 'Threshold', thresholdBurst_opt);
    
    % plot results
        fig = figure('Color', 'w','Position',[0 0 1600 800],'Visible','off');
        
    
        %subplot(1,1,1);
        yyaxis right
        mxw.plot.networkActivity(networkAct_opt, 'Threshold', thresholdBurst, 'Figure', false,'Color','#000000')
        %box off;
        %hold on;
        plot(networkStats_opt.maxAmplitudesTimes,networkStats_opt.maxAmplitudesValues,'or');
        ylim([0 30])
        %xlim([0 round(max(relativeSpikeTimes.time)/4)])
        %xlim([0 60])
        hold on;
    
        %subplot(1,1,1);
        yyaxis left
        hold on;
        mxw.plot.rasterPlot(networkData, 'Figure', false, 'Color','#c6c6c6') %
        box off;
        %xlim([0 round(max(relativeSpikeTimes.time)/4)])
        xlim([0 120])
        ylim([1 max(relativeSpikeTimes.channel)])
        hold on;
        title(IDstring +' Gaussian = ' + string(j))
    
        exportFile = [opDir 'figcompare.pdf']; %folder and filename for raster figures
        exportgraphics(fig, exportFile ,'Append',true,'Resolution',150)
        toc    
    end
end
toc
%% PSD code
%signal = reshape(networkStats_opt.maxAmplitudesTimes,1,[]);
% signal = networkStats_opt(:).';
% n = length(signal);
% x = signal;
% 
% pxx = pwelch(x);
% pwelch(x)

%https://stackoverflow.com/questions/21750075/matlab-python-power-spectral-density-of-non-uniform-time-series


