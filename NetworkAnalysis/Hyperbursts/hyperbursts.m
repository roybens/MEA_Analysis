%Hyperbusrts adpated from network activity originally  by Elena Gronskaya [updated 25.11.2020]
%modified by Mandar Patil (github : mandarmp)



clear
close all

% set paths to the Network recording file
pathFileNetwork =  '/mnt/disk20tb/PrimaryNeuronData/Maxtwo/ADNP_2/ADNP_2/ADNP_2/231106/M06792/Network/000018/data.raw.h5';

% select which well to analyze (an integer between 1 and 6; leave as 1 for MaxOne data)
wellID = 1;

% create fileManager object for the Network recording
networkData = mxw.fileManager(pathFileNetwork,wellID);

%% Visualise the rasters to find the duration of hyperbursts.


relativeSpikeTimes = mxw.util.computeRelativeSpikeTimes(networkData);

% Let's use the time and channel vectors to visualise our data in a raster plot.
starttime = 100;
%endtime = round(max(relativeSpikeTimes.time)/4);
endtime = 160;

figure('Color','w','position',[0 0 400 400]);

plot(relativeSpikeTimes.time,relativeSpikeTimes.channel,'.','MarkerSize',2,'Color','#135ba3')
ylabel('Channel')
xlabel('Time [s]')
title('Raster Plot','fontsize',11)
xlim([starttime endtime])
ylim([1 max(relativeSpikeTimes.channel)])
box off;

% In a well-interconnected neuronal culture, bursts of activity will often be 
% visible by eye. In order to detect them automatically and to quantify their 
% amplitude, we take a three-step approach. First, we bin all the spike times 
% into small time windows (size adjusted by the parameter _binSize_). Note that 
% from this point on, we disregard channel number, treating all the spikes together 
% as network activity. 

%% Now choose the hyperburst intervals and adjust parameters

starttime = 100;
endtime= 160;
binSize = 0.1;
gaussianSigma = 0.5;
thresholdBurst = 1; % in rms of the firing rate
MinPeakDistance = 2;


thres_fn = 'dynamic';

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
val = -3*gaussianSigma:binSize:3*gaussianSigma;

kernel = mxw.util.normpdf(-3*gaussianSigma:binSize:3*gaussianSigma,0,gaussianSigma); 
kernel = kernel*binSize;
firingRate = conv(binnedTimes,kernel,'same');
firingRate = firingRate/binSize;
firingRateNorm = firingRate/length(unique(networkData.rawMap.spikes.channel)); 

% Let's plot the Network Activity below the raster. 


figure('Color','w','position',[0 0 400 400]);
subplot(2,1,1);
plot(relativeSpikeTimes.time,relativeSpikeTimes.channel,'.','MarkerSize',2,'Color','#135ba3')
ylabel('Channel')
xlabel('Time [s]')
title('Raster Plot','fontsize',11)
xlim([starttime endtime])
ylim([1 max(relativeSpikeTimes.channel)])
box off;
subplot(2,1,2);
plot(timeVector,firingRateNorm,'Color','#135ba3')
xlim([starttime endtime])
ylabel('Firing Rate [Hz]')
xlabel('Time [s]')
title('Network Activity','fontsize',11)
hold on;

% As a third step, we carry out peak detection on the network activity curve. 
% Here we must set a third parameter, _thresholdBurst:_ the threshold, in root-mean-square 
% of the network firing rate, above which the bursts will be detected.  


rmsFiringRate = mxw.util.rms(firingRateNorm);

if strcmp(thres_fn ,'dynamic')
[tmpTimes, burstPeakValues] = findPeaks(firingRateNorm,...
    'MinPeakHeight',thresholdBurst * rmsFiringRate, 'MinPeakDistance', MinPeakDistance);
burstPeakTimes = timeVector(tmpTimes);
plot(thresholdBurst*rmsFiringRate*ones(ceil(timeVector(end)),1))
else
[tmpTimes, burstPeakValues] = findPeaks(firingRateNorm,...
    'PositiveThreshold', thresholdBurst, 'MinPeakDistance', MinPeakDistance);
burstPeakTimes = timeVector(tmpTimes);
plot(thresholdBurst*ones(ceil(timeVector(end)),1))
end


plot(burstPeakTimes,burstPeakValues,'or')





%% These are calculations of the Network metrics.


% Burst Peak Amplitude and Interburst Interval
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