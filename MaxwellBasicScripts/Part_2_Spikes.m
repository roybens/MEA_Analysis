%% MaxWell Matlab Toolboxes Tutorials
% by Elena Gronskaya [updated 25.11.2020]
% 
%% PART II: Spikes
% In this tutorial, we will see how to plot the spike data in time (raster plots) 
% and space (activity maps). We will also look at some basic spike stats: firing 
% rate and spike amplitude distributions for all the recorded electrodes.
% 
% Sample data can be downloaded from here: <https://share.mxwbio.com/d/dfcccbd3e5a74d2a8ae8/ 
% https://share.mxwbio.com/d/dfcccbd3e5a74d2a8ae8/> (the folders 'scan' and 'network')
% 
% For this tutorial, please make sure to use the scan_data_spike_only file. 
% (For the Network recording, the same file is always used).  
% 
% Let's start with loading the data and creating fileManager objects.  

clear
close all

% set paths to the Activity Scan  and the Network recording files
pathFileActivityScan = '/home/mmpatil/Downloads/scan_data_full.raw.h5'; 
pathFileNetwork = '/home/mmpatil/Downloads/network_data_full.raw.h5'; 

% select which well to analyze (an integer between 1 and 6; leave as 1 for MaxOne data)
wellID = 1;

% create fileManager object for the Activity Scan
activityScanData = mxw.fileManager(pathFileActivityScan,wellID);
% create fileManager object for the Network recording
networkData = mxw.fileManager(pathFileNetwork,wellID);
%% Plotting the Spike Data
% I. Raster Plots
% 
% To visualise the activity of the neural population over time, we can plot 
% a raster plot. In fileManager.fileObj, we can find all the information needed 
% for a raster:

% - the sampling rate
% - the spike times (in frame numbers)
% - the first frame number of the recording, relative to which all the spike 
% frame numbers are reported
% - the channels on which the spikes were recorded 

% Let's compare a raster plot from one of the Activity Scan recording files 
% with a raster plot from a Network recording, which uses the most active (highest 
% firing-rate) electrodes. 

% The recording number should be set according to its order in the 'fileObj" array. 
recNr = 3;
% sampling frequency
fsActivity = activityScanData.fileObj(recNr).samplingFreq;
% get spike time stamps 
tsActivity = double(activityScanData.fileObj(recNr).spikes.frameno - activityScanData.fileObj(recNr).firstFrameNum)/fsActivity;
% channel list, where spike time stamps where detected
chActivity = activityScanData.fileObj(recNr).spikes.channel;

% plot raster
figure('Color','w');
subplot(2,1,1)
plot(tsActivity, chActivity,'.','Color','#135ba3','MarkerSize',2)
box off 
xlabel('Time [s]') 
ylabel('Channel')
title('Activity Scan Raster Plot')
xmaxActivity = max(tsActivity);
xlim([0 xmaxActivity])
ylim([0 max(chActivity)])

% sampling frequency
fsNetwork = networkData.fileObj.samplingFreq;
% get spike time stamps 
tsNetwork = double(networkData.fileObj.spikes.frameno - networkData.fileObj.firstFrameNum)/fsNetwork;
% channel list, where spike time stamps where detected
chNetwork = networkData.fileObj.spikes.channel;

% plot raster
subplot(2,1,2)
plot(tsNetwork, chNetwork,'.','Color','#135ba3','MarkerSize',2)
box off 
xlabel('Time [s]') 
ylabel('Channel')
title('Network Recording Raster Plot')
xmaxNetwork = max(tsNetwork);
xlim([0 xmaxNetwork]) 
ylim([0 max(chNetwork)])
%% 
% II. Spike Rate Activity Map and Distribution
% 
% From the spike times, we can also calculate the mean spike rate of each electrode, 
% and see how the spike rate is distributed over the area of the chip (an "activity 
% map"). Typically, data from the Activity Scan is used to visualise the activity 
% map. 

% get the mean firing rate for each electrode
meanFiringRate = mxw.activityMap.computeSpikeRate(activityScanData);

% define maximum firing rate used in plots as 99th percentile of firing
% rate values
maxFiringRate = mxw.util.percentile(meanFiringRate(meanFiringRate~=0),99);

% plot the list of mean firing rates as a map, given the electrode 
% x and y coordinates in 'fileManagerObj.processedMap'
figure('Color','w');
subplot(2,1,1);
mxw.plot.activityMap(activityScanData, meanFiringRate, 'Ylabel', '[Hz]',...
    'CaxisLim', [0.1 max(meanFiringRate)/5],'Figure',false,'Title','Firing Rate Activity Map');
% run the line above several times, experimenting with different 
% [min max] range of the color gradient 'CaxisLim'

% add scale bar 
line([300 800],[2000+400 2000+400],'Color','k','LineWidth',4);
axis off;
text(340,2100+500,'0.5 mm','color','k');
xlim([200 3750]);ylim([150 2500])

% plot the distribution of firing rates across all of the electrodes
subplot(2,1,2);
thrFiringRate = 0.1; % set a minimum spike rate threshold (Hz)
histogram(meanFiringRate(meanFiringRate>thrFiringRate),0:.1:ceil(maxFiringRate), ...
    'FaceColor','#135ba3', "EdgeColor",'#414042')
xlim([thrFiringRate maxFiringRate])
ylabel('Counts');xlabel('Firing Rate [Hz]');
box off;
legend(['Mean Firing Rate = ',num2str(mean(meanFiringRate(meanFiringRate>thrFiringRate)),'%.2f'),...
    ' Hz,  sd = ',num2str(std(meanFiringRate(meanFiringRate>thrFiringRate)),'%.2f')])
%% 
% III. Spike Amplitude Activity Map and Distribution
% 
% An activity map can also be plotted on the basis of spike amplitude. We use 
% the 90th percentile of all the spike amplitudes recorded on each electrode (a 
% more appropriate measure than the mean amplitude, since the multi-unit activity 
% picked up by the electrodes is dominated by smaller spikes). 

% get th 90th percentile spike amplitude value for each electrode
amplitude90perc = abs(mxw.activityMap.computeAmplitude90percentile(activityScanData));

% define maximum amplitude used in plots as 99th percentile of amplitude values
maxAmp = mxw.util.percentile(amplitude90perc(amplitude90perc~=0),99);

% plot the mean firing rate vector as a map, given the electrode 
% x and y coordinates in 'fileManagerObj.processedMap'
figure('color','w');
subplot(2,1,1);
mxw.plot.activityMap(activityScanData, amplitude90perc,'Ylabel', '[\muV]',...
    'CaxisLim', [10 maxAmp], 'Figure',false,'Title','Spike Amplitude Activity Map');
% run the line above several times, experimenting with different 
% [min max] range of the color gradient 'CaxisLim'

% add scale bar 
line([300 800],[2000+400 2000+400],'Color','k','LineWidth',4);
axis off;
text(340,2100+500,'0.5 mm','color','k');
xlim([200 3750])
ylim([150 2500])

% plot the distribution of spike amplitudes across all of the electrodes
subplot(2,1,2);
thrAmp = 10;  % set a minimum spike amplitude threshold (uV)
histogram(amplitude90perc(amplitude90perc>thrAmp),ceil(0:1:maxAmp), ...
    'FaceColor','#135ba3', "EdgeColor",'#414042')
xlim([thrAmp maxAmp])
ylabel('Counts');xlabel('Spike Amplitude [\muV]');
box off;
legend(['Mean Spike Amplitude = ',num2str(mean(amplitude90perc(amplitude90perc>thrAmp)),'%.2f'),...
    ' \muV,  sd = ',num2str(std(amplitude90perc(amplitude90perc>thrAmp)),'%.2f')])