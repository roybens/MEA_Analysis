clear all
close all


matFile = '/home/mmp/disktb/mmpatil/MEA_Analysis/IPNAnalysis/Organoid.mat';

load(matFile);

load

% Burst Parameters.


% set Gaussian kernel standard deviation [s] (smoothing window)
gaussianSigma = 0.16; %0.18
% set histogram bin size [s]
binSize = 0.02;
% set minimum peak distance [s]
minPeakDistance = 3.0;
% set burst detection threshold [rms / fixed]
thresholdBurst =1.5; %1.2
% set fixed threshold;
use_fixed_threshold = false;
% Set the threshold to find the start and stop time of the bursts. (start-stop threshold)
thresholdStartStop = 1; %0.3

dataFile = '/mnt/harddrive-2/Organoids_Mandeep_Fink_Lab/Cdkl5_Organoids_Mano0855-D/221118/16719/Network/000032/data.raw.h5';

% Set Threshold function for later use
threshold_fn = 'Threshold';
if use_fixed_threshold
    threshold_fn = 'FixedThreshold';
end

error_l =[];
% create fileManager object for the Network recording
relativeSpikeTimes.time = [];
relativeSpikeTimes.channel = [];
relativeSpikeTimes.mappedChannel =[];
try
    networkData = mxw.fileManager(dataFile);
catch
    error_l = [error_l string(scan_runID_text)];
    
end

spikeTimes = double(spiking_data.spike_frames)/ double(networkData.fileObj.samplingFreq);
firstFrame = min(spikeTimes);

relativeSpikeTimes.time= spikeTimes - firstFrame;
relativeSpikeTimes.channel= spiking_data.spike_new_channels;

unique_values = unique(relativeSpikeTimes.channel);
unique_values =sort(unique_values);
value_mapping=containers.Map(unique_values,1:numel(unique_values));
relativeSpikeTimes.mappedChannel = cellfun(@(x) value_mapping(x),num2cell(relativeSpikeTimes.channel));



networkAct = mxw.networkActivity.computeNetworkAct(relativeSpikeTimes, 'BinSize', binSize,'GaussianSigma', gaussianSigma);
networkStats = computeNetworkStats_JL(networkAct, threshold_fn, thresholdBurst, 'MinPeakDistance', minPeakDistance);
    
% Let's use the time and channel vectors to visualise our data in a raster plot.

figure('Color','w','position',[0 0 400 800]);
subplot(2,1,1);
plot(relativeSpikeTimes.time,relativeSpikeTimes.mappedChannel,'.','MarkerSize',2,'Color','#135ba3')
ylabel('Channel')
xlabel('Time [s]')
title('Raster Plot','fontsize',11)
xlim([0 round(max(relativeSpikeTimes.time)/4)])
ylim([1 max(relativeSpikeTimes.mappedChannel)])
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
ylim([0 2])
ylabel('Firing Rate [Hz]')
xlabel('Time [s]')
title('Network Activity','fontsize',11)
hold on;
thresholdBurst = 1; % in rms of the firing rate
minPeakDist =3.0;
rmsFiringRate = mxw.util.rms(firingRateNorm);

[peaks, peakIndices] = findpeaks(firingRateNorm,timeVector,'MinPeakHeight',...
    thresholdBurst*rmsFiringRate, 'MinPeakDistance', minPeakDist);
% 
% [tmpTimes, burstPeakValues] = mxw.util.findPeaks(firingRateNorm,...
%     'PositiveThreshold', thresholdBurst*rmsFiringRate);
% burstPeakTimes = timeVector(peakIndices);
% burstPeakDiffTimes = diff(burstPeakTimes);
% burstPeakDiffAmpValues =  diff(burstPeakValues);
% indices = find(burstPeakDiff<minPeakDist);
% for i = 1:length(indices)
%     idx = indices(i);
%     if burstPeakDiffAmpValues(i) < 0
%         
%     end
% 
% end
% burstPeakValues(indices)=[];
% burstPeakTimes(indices)=[];
plot(thresholdBurst*rmsFiringRate*ones(ceil(timeVector(end)),1))
%plot(burstPeakTimes,burstPeakValues,'or')
plot(peakIndices,peaks,'or')

saveas(gcf,'~/Documents/Images_25sep/SpikeSortedNetwork.pdf','pdf')


% 
% %% Tim's code for averaging and aggregating mean spiking data (IBI, Burst peaks, Spikes within Bursts, # of Bursts etc.)
% %average IBI
% meanIBI = mean(networkStats.maxAmplitudeTimeDiff);
% %average Burst peak (burst firing rate y-value)
% meanBurstPeak = mean(networkStats.maxAmplitudesValues);
% %Number of bursts
% nBursts = length(networkStats.maxAmplitudesTimes);
% 
% 
% 
% %average spikesPerBurst        
% if length(networkStats.maxAmplitudesTimes)>3
%     peakAmps = networkStats.maxAmplitudesValues';
%     peakTimes = networkStats.maxAmplitudesTimes;
%     
%     % get the times of the burst start and stop edges
%     edges = double.empty(length(peakAmps),0);
%     for i = 1:length(peakAmps)
%        % take a sizeable (±6 s) chunk of the network activity curve 
%        % around each burst peak point
%        idx = networkAct.time>(peakTimes(i)-6) & networkAct.time<(peakTimes(i)+6);
%        t1 = networkAct.time(idx);
%        a1 = networkAct.firingRate(idx)';
%       
%        % get the amplitude at the desired peak width
%        peakWidthAmp = (peakAmps(i)-round(peakAmps(i)*thresholdStartStop));
%        
%        % get the indices of the peak edges
%        idx1 = find(a1<peakWidthAmp & t1<peakTimes(i));
%        idx2 = find(a1<peakWidthAmp & t1>peakTimes(i));
%        
%        if ~isempty(idx1)&&~isempty(idx2)       
%            tBefore = t1(idx1(end));
%            tAfter = t1(idx2(1));
%            edges(i,[1 2]) = [tBefore tAfter];
%        end
%     end
%     if spikeSorted
%     
%     ts = ((double(spikeTimes)...
%         - double(firstFrame))/networkData.fileObj.samplingFreq)';
%     ch = spiking_data.units;
% 
%     else
%    % identify spikes that fall within the bursts
%     ts = ((double(networkData.fileObj.spikes.frameno)...
%         - double(networkData.fileObj.firstFrameNum))/networkData.fileObj.samplingFreq)';
%     ch = networkData.fileObj.spikes.channel;
%     end
%     spikesPerBurst = double.empty(length(edges),0);
%     tsWithinBurst = [];
%     chWithinBurst = [];
%     for i = 1:length(edges)
%        idx = (ts>edges(i,1) & ts<edges(i,2));
%        spikesPerBurst(i) = sum(idx); 
%        tsWithinBurst = [tsWithinBurst ts(idx)];
%        chWithinBurst = [chWithinBurst ch(idx)'];
%     end
%     meanSpikesPerBurst = mean(spikesPerBurst);
% end
%    
% % append information to table elements
% Run_ID = [Run_ID scan_runID];
% DIV = [DIV scan_div];
% Time = [Time hd5Date];
% Chip_ID = [Chip_ID scan_chipID];
% IBI = [IBI meanIBI];
% Burst_Peak = [Burst_Peak meanBurstPeak];
% Number_Bursts = [Number_Bursts nBursts];
% Spike_per_Burst = [Spike_per_Burst meanSpikesPerBurst];
% runIDstemp = run_id_and_type.Run_;
% types = run_id_and_type.NeuronSource;
% index = find(runIDstemp == scan_runID);
% targetType = types{index};
% % plot results
%     figure('Color','w','Position',[0 0 400 800],'Visible','off');
%     subplot(2,1,1);
%     mxw.plot.rasterPlot(relativeSpikeTimes,'Figure',false);
%     box off;
%     %xlim([0 round(max(relativeSpikeTimes.time)/4)])
%     xlim([0 120])
%     ylim([1 max(relativeSpikeTimes.channel)])
%     
%     subplot(2,1,2);
%     mxw.plot.networkActivity(networkAct,'Threshold',thresholdBurst,'Figure',false);
%     box off;
%     hold on;
%     plot(networkStats.maxAmplitudesTimes,networkStats.maxAmplitudesValues,'or')
%     %xlim([0 round(max(relativeSpikeTimes.time)/4)])
%     xlim([0 120])
%     ylim([0 20])
%            % saveas(gcf,append(opDir,'Network_outputs/Raster_BurstActivity/Raster_BurstActivity',scan_runID_text,'_',num2str(scan_chipID),'_DIV',num2str(scan_div),'_',targetType,'.png'))
