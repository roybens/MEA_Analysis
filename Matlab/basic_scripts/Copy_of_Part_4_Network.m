clear
close all

% set paths to the Network recording file
pathFileNetwork =  '/mnt/disk15tb/jonathan/Syngap3/Syngap3/230109/16657/Network/000178/data.raw.h5';

% select which well to analyze (an integer between 1 and 6; leave as 1 for MaxOne data)
wellID = 1;

% create fileManager object for the Network recording
networkData = mxw.fileManager(pathFileNetwork,wellID);

electrodesInterest = [9703, 13043, 15450, 13186, 13879, 16705, 12284, 9744, 9767, 12357, 15153, 12357, 12357, 8422, 13385, 13254, 13277, 7291, 11634, 5065, 5055, 5292, 19424, 5519, 6337, 8759, 11950, 16913, 8376, 18524, 16785, 16124, 18087, 23103, 8576, 17640, 21399, 688, 24879, 9416, 22207, 24019, 5279, 2637, 11748, 13687, 16435, 1134, 2838, 10051, 19387, 2844, 17369, 22076, 13259, 5261, 24037, 9515, 22938, 26229, 12177, 25991, 20523, 8269, 10406, 25925, 76, 605, 8957, 2380, 4695, 22919, 21601, 8730, 9979, 4580, 10654, 20343, 21288, 20164, 3342, 24681, 18797, 26330, 13339, 25029, 11064, 3907, 12622, 21417, 1467, 21300, 1050, 7846, 20546, 7370, 7370, 19671, 19671, 18159, 26324, 11101, 21865, 2340, 25401, 21004, 6477, 2553, 5404, 24738, 9568, 12001, 769, 14174, 9572, 16613, 7815, 6494, 12666, 9781, 13964, 4728, 5187, 11767, 2730, 4065];

matchingIndices =ismember(networkData.rawMap.map.electrode,electrodesInterest);

%updating networkData.rawMap

networkData.rawMap.map.channel=networkData.rawMap.map.channel(matchingIndices);
networkData.rawMap.map.electrode=networkData.rawMap.map.electrode(matchingIndices);
networkData.rawMap.map.x=networkData.rawMap.map.x(matchingIndices);
networkData.rawMap.map.y=networkData.rawMap.map.y(matchingIndices);


%updating the extracted spikes.
networkData.extractedSpikes.frameno = networkData.extractedSpikes.frameno(matchingIndices);
networkData.extractedSpikes.amplitude = networkData.extractedSpikes.amplitude(matchingIndices);

%now get channels mapping.
matchingChannelIndices = ismember(networkData.rawMap.spikes.channel,networkData.rawMap.map.channel);
networkData.rawMap.spikes.channel = networkData.rawMap.spikes.channel(matchingChannelIndices);
networkData.rawMap.spikes.frameno = networkData.rawMap.spikes.frameno(matchingChannelIndices);
networkData.rawMap.spikes.amplitude = networkData.rawMap.spikes.amplitude(matchingChannelIndices);

%updating the networkData.fileObj
networkData.fileObj.map = networkData.rawMap.map;
networkData.fileObj.spikes = networkData.rawMap.spikes;


%updating networkData.processsedMap
matchingIndices =ismember(networkData.processedMap.electrode,electrodesInterest);
networkData.processedMap.electrode = networkData.processedMap.electrode(matchingIndices);
networkData.processedMap.xpos = networkData.processedMap.xpos(matchingIndices);
networkData.processedMap.ypos = networkData.processedMap.ypos(matchingIndices);
networkData.processedMap.recordingIndex = ones(length(networkData.processedMap.electrode),1);
networkData.processedMap.nonRoutedElec= (0:26399)';
indicesToRemove = ismember(networkData.processedMap.nonRoutedElec, networkData.processedMap.electrode);
networkData.processedMap.nonRoutedElec(indicesToRemove)=[];

%% Network Activity and Burst Detection
% Next, we use a function that converts the spike times from frame numbers into 
% seconds, referenced to the start time of the recording. The spike times from 
% all of the channels are cancatenated into a single long vector and the channel 
% numbers are listed in a corresponding vector. Both variables are stored in the 
% structure relativeSpikeTimes.

relativeSpikeTimes = mxw.util.computeRelativeSpikeTimes(networkData);

% Let's use the time and channel vectors to visualise our data in a raster plot.

figure('Color','w','position',[0 0 400 800]);
subplot(2,1,1);
plot(relativeSpikeTimes.time,relativeSpikeTimes.channel,'.','MarkerSize',2,'Color','#135ba3')
ylabel('Channel')
xlabel('Time [s]')
title('Raster Plot','fontsize',11)
xlim([0 round(max(relativeSpikeTimes.time)/4)])
ylim([1 max(relativeSpikeTimes.channel)])
box off;

% In a well-interconnected neuronal culture, bursts of activity will often be 
% visible by eye. In order to detect them automatically and to quantify their 
% amplitude, we take a three-step approach. First, we bin all the spike times 
% into small time windows (size adjusted by the parameter _binSize_). Note that 
% from this point on, we disregard channel number, treating all the spikes together 
% as network activity. 

binSize = 0.02;
timeVector = 0:binSize:max(relativeSpikeTimes.time);
[binnedTimes, ~] = histcounts(relativeSpikeTimes.time,timeVector);
binnedTimes(end+1) = 0;
binnedTimes = binnedTimes.';

% Second, we convolve the resulting histogram with a Gaussian kernel to produce 
% a smoothed curve of network activity. The parameter gaussianSigma, which is 
% the kernel standard deviation in seconds, determines how much the curve will 
% be smoothed; with smaller values, finer details of the burst dynamics will be 
% resolved (ie, onset and offset peaks in activity bursts), although noise may 
% be captured as well. We normalise network activity by the number of active electrodes, 
% which may vary widely from culture to culture. 

gaussianSigma = 0.14;
kernel = mxw.util.normpdf(-3*gaussianSigma:binSize:3*gaussianSigma,0,gaussianSigma); 
kernel = kernel*binSize;
firingRate = conv(binnedTimes,kernel,'same');
firingRate = firingRate/binSize;
firingRateNorm = firingRate/length(unique(networkData.rawMap.spikes.channel)); 
%save('NetworkAct16795.mat',"firingRateNorm")
% Let's plot the Network Activity below the raster. 

subplot(2,1,2);
plot(timeVector,firingRateNorm,'Color','#135ba3')
xlim([0 round(max(relativeSpikeTimes.time)/4)])
ylabel('Firing Rate [Hz]')
xlabel('Time [s]')
title('Network Activity','fontsize',11)
hold on;

electrodesInterest = [9703, 13043, 15450, 13186, 13879, 16705, 12284, 9744, 9767, 12357, 15153, 12357, 12357, 8422, 13385, 13254, 13277, 7291, 11634, 5065, 5055, 5292, 19424, 5519, 6337, 8759, 11950, 16913, 8376, 18524, 16785, 16124, 18087, 23103, 8576, 17640, 21399, 688, 24879, 9416, 22207, 24019, 5279, 2637, 11748, 13687, 16435, 1134, 2838, 10051, 19387, 2844, 17369, 22076, 13259, 5261, 24037, 9515, 22938, 26229, 12177, 25991, 20523, 8269, 10406, 25925, 76, 605, 8957, 2380, 4695, 22919, 21601, 8730, 9979, 4580, 10654, 20343, 21288, 20164, 3342, 24681, 18797, 26330, 13339, 25029, 11064, 3907, 12622, 21417, 1467, 21300, 1050, 7846, 20546, 7370, 7370, 19671, 19671, 18159, 26324, 11101, 21865, 2340, 25401, 21004, 6477, 2553, 5404, 24738, 9568, 12001, 769, 14174, 9572, 16613, 7815, 6494, 12666, 9781, 13964, 4728, 5187, 11767, 2730, 4065];
keepChannels = networkData.rawMap.map.channel(ismember(networkData.rawMap.map.electrode,electrodesInterest));

channelsToKeep = ismember(relativeSpikeTimes.channel,keepChannels);

% Extract the filtered time and channel data
filteredTime = relativeSpikeTimes.time(channelsToKeep);
filteredChannels = relativeSpikeTimes.channel(channelsToKeep);

figure('Color','w','position',[0 0 400 800]);
subplot(3,1,1);
plot(filteredTime,filteredChannels,'.','MarkerSize',2,'Color','#135ba3')
ylabel('Channel')
xlabel('Time [s]')
title('Raster Plot','fontsize',11)
xlim([0 round(max(relativeSpikeTimes.time)/4)])
ylim([1 max(relativeSpikeTimes.channel)])
box off;

binSize = 0.02;
timeVector = 0:binSize:max(filteredTime);
[binnedTimes, ~] = histcounts(filteredTime,timeVector);
binnedTimes(end+1) = 0;
binnedTimes = binnedTimes.';

% Second, we convolve the resulting histogram with a Gaussian kernel to produce 
% a smoothed curve of network activity. The parameter gaussianSigma, which is 
% the kernel standard deviation in seconds, determines how much the curve will 
% be smoothed; with smaller values, finer details of the burst dynamics will be 
% resolved (ie, onset and offset peaks in activity bursts), although noise may 
% be captured as well. We normalise network activity by the number of active electrodes, 
% which may vary widely from culture to culture. 

gaussianSigma = 0.14;
kernel = mxw.util.normpdf(-3*gaussianSigma:binSize:3*gaussianSigma,0,gaussianSigma); 
kernel = kernel*binSize;
firingRate = conv(binnedTimes,kernel,'same');
firingRate = firingRate/binSize;
firingRateNorm = firingRate/length(unique(filteredChannels)); 

subplot(3,1,2);
plot(timeVector,firingRateNorm,'Color','#135ba3')
xlim([0 round(max(filteredTime)/4)])
ylabel('Firing Rate [Hz]')
xlabel('Time [s]')
title('Network Activity','fontsize',11)
hold on;


% As a third step, we carry out peak detection on the network activity curve. 
% Here we must set a third parameter, _thresholdBurst:_ the threshold, in root-mean-square 
% of the network firing rate, above which the bursts will be detected.  

thresholdBurst = 1; % in rms of the firing rate
rmsFiringRate = mxw.util.rms(firingRateNorm);

[tmpTimes, burstPeakValues] = mxw.util.findPeaks(firingRateNorm,...
    'PositiveThreshold', thresholdBurst*rmsFiringRate);
burstPeakTimes = timeVector(tmpTimes);
plot(thresholdBurst*rmsFiringRate*ones(ceil(timeVector(end)),1))
plot(burstPeakTimes,burstPeakValues,'or')

%% Parameter Adjustment

% With the parameters chosen, you may find that some smaller-amplitude bursts 
% may have been missed. Below, play around with the parameters to find an optimal 
% set of values that best capture the network bursting activity. In order to speed 
% up testing, we will use two dedicated functions of the MaxWell Matalab Toolbox 
% that perform all of the steps we have covered above, as well as two dedicated 
% plotting functions. 

% set histogram bin size [s]
binSize = 0.1;

% set Gaussian kernel standard deviation [s]
gaussianSigma = 0.1;

% set burst detection threshold [rms firing rate]
thresholdBurst = 0.4;

% compute network activity and detect bursts
networkAct = mxw.networkActivity.computeNetworkAct(networkData,'BinSize',binSize,'GaussianSigma',gaussianSigma);
networkStats = mxw.networkActivity.computeNetworkStats(networkAct,'Threshold',thresholdBurst);

% plot results
figure('Color','w','Position',[0 0 400 800]);
subplot(2,1,1);
mxw.plot.rasterPlot(networkData,'Figure',false);
box off;
xlim([0 round(max(relativeSpikeTimes.time)/4)])
ylim([1 max(relativeSpikeTimes.channel)])

subplot(2,1,2);
mxw.plot.networkActivity(networkAct,'Threshold',thresholdBurst,'Figure',false);
box off;
hold on;
plot(networkStats.maxAmplitudesTimes,networkStats.maxAmplitudesValues,'or')
xlim([0 round(max(relativeSpikeTimes.time)/4)])

%% Burst Peak Amplitude and Interburst Interval
% Let's take a look at the distributions of the burst peak values and the interburst 
% intervals. 

% make sure at least three bursts were detected
if length(networkStats.maxAmplitudesTimes)>3
    % Burst Peak
    mxw.plot.networkStats(networkStats,'Option','maxAmplitude',...
        'Figure',true,'Ylabel','Counts','Xlabel','Burst Peak [Hz]',...
        'Title','Burst Peak Distribution','Bins',20);
    box off;
    legend(['Mean Burst Peak = ',num2str(mean(networkStats.maxAmplitudesValues),'%.2f Hz'),...
        ', sd = ',num2str(std(networkStats.maxAmplitudesValues),'%.2f')])
    
    % IBI
    mxw.plot.networkStats(networkStats,'Option','maxAmplitudeTimeDiff',...
        'Figure',true,'Ylabel','Counts','Xlabel','Interburst Interval [s]',...
        'Title','Interburst Interval Distribution','Bins',20);
    box off;
    legend(['Mean Interburst Interval = ',num2str(mean(networkStats.maxAmplitudeTimeDiff),'%.2f s'),...
        ', sd = ',num2str(std(networkStats.maxAmplitudeTimeDiff),'%.2f')])
end
%% Burst Duration and the Percentage of Spikes within Bursts
% Finally, we can identify the start and stop times of the bursts and look at 
% the distribution of the burst duration as well as the number of spikes per burst. 

thresholdStartStop = 0.4;
% threshold to find the start and stop time of the bursts
% 0.3 means 30% value of the burst peak. Note that by raising
% the value, the percentage of spikes within bursts increases,
% since the bursts are considered wider. 

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
    
    % plot the distribution of the number of spikes per burst
    figure('Color','w');
    h = histogram(spikesPerBurst,20);
    h.FaceColor = '#719cc7';
    h.EdgeColor = '#414042';
    h.FaceAlpha = 1;
    box off
    ylabel('Counts')
    xlabel('Number of Spikes per Burst')
    title(['Spikes within Bursts = ', num2str(sum(spikesPerBurst/length(ts))*100,'%.1f'),' %'],...
        'FontSize',11)
    legend(['Mean Spikes per Burst = ',num2str(mean(spikesPerBurst),'%.2f'),...
        ', sd = ',num2str(std(spikesPerBurst),'%.2f')])
    
    % Burst Duration
    figure('Color','w');
    h = histogram(abs(edges(:,1) - edges(:,2)),20);
    h.FaceColor = '#719cc7';
    h.EdgeColor = '#414042';
    h.FaceAlpha = 1;
    box off
    ylabel('Counts')
    xlabel('Time [s]')
    title('Burst Duration', 'FontSize',11)
    legend(['Mean Burst Duration = ',num2str(mean(abs(edges(:,1)...
        - edges(:,2))),'%.2f'), ' s, sd = ',num2str(std(abs(edges(:,1)...
        - edges(:,2))),'%.2f')])
end