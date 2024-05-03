clear
close all



% set paths to the Network recording file
pathFileNetwork =  '~/Documents/MATLAB/mmpatil/network_data_full.raw.h5';
pathFileActivityScan = '~/Documents/MATLAB/mmpatil/scan_data_full.raw.h5';

% select which well to analyze (an integer between 1 and 6; leave as 1 for MaxOne data)
wellID = 1;

% create fileManager object for the Network recording
networkData = mxw.fileManager(pathFileNetwork,wellID);
activityScanData = mxw.fileManager(pathFileActivityScan,wellID);

startSample = 1; % sample number from which to start data extraction
dataChunkSize = 5; % length of the data (s) you want to extract

%  get the sampling rate
fsNetwork = networkData.fileObj.samplingFreq;


xpos = networkData.processedMap.xpos;
ypos = networkData.processedMap.ypos;
electrodes = networkData.processedMap.electrode;
spikeCount = mxw.activityMap.computeSpikeCount(networkData);
spikeRate = mxw.activityMap.computeSpikeRate(networkData);
%% get spike time stamps
tsNetwork = double(networkData.fileObj.spikes.frameno - networkData.fileObj.firstFrameNum)/fsNetwork;
% channel list, where spike time stamps where detected
chNetwork = networkData.fileObj.spikes.channel;

%relative spike time wrt to first spike
relativetsNetwork = round(double((networkData.fileObj.spikes.frameno-networkData.fileObj.spikes.frameno(1))));

deltaSamples = 5;
a = [inf; diff(relativetsNetwork)];
b = find(a<=deltaSamples);   %% kinda has the index
chs = chNetwork(b);
d = find(a>deltaSamples);

% extract unfiltered Network data
% [ntwTracesRaw, ~, ~] = networkData.extractRawData(startSample,dataChunkSize*fsNetwork);
% [electrode_groups, channel_groups] = circus.electrodeGroups(networkData.rawMap.map, 100);
% 
% %get 10% high spiking electrodes



% extract unfiltered Network data
% [ntwTracesRaw, ~, ~] = networkData.extractRawData(startSample,dataChunkSize*fsNetwork);
% [electrode_groups, channel_groups] = circus.electrodeGroups(networkData.rawMap.map, 100);
% 
% %get 10% high spiking electrodes
% spikeCount = mxw.activityMap.computeSpikeCount(networkData);
% spikeRate = mxw.activityMap.computeSpikeRate(activityScanData);
% activityThreshold = 10;
% activityThrValue = mxw.util.percentile(spikeCount, 100 - activityThreshold);
% selectedElectrodeTruthValues = (spikeCount > activityThrValue);
% electrodes = networkData.processedMap.electrode(selectedElectrodeTruthValues);
% xpos = networkData.processedMap.xpos(selectedElectrodeTruthValues);
% ypos = networkData.processedMap.ypos(selectedElectrodeTruthValues);
% 
% %need estimation of max how many electrodes may be present.???
% 
% for i=1:length(electrodes)
%     
%     [neighEls, neighInd] = mxw.util.get_neighbouring_electrodes(electrodes(i), networkData, 10*10);
%     selectedElectrodes(i).elec = electrodes(i);
%     selectedElectrodes(i).xpos = xpos(i);
%     selectedElectrodes(i).ypos = ypos(i);
%     selectedElectrodes(i).neighEls = neighEls';
%     selectedElectrodes(i).neighInd = neighInd;
% 
% end  

% relativeSpikeValues  = mxw.util.computeRelativeSpikeTimes( networkData );

%selected electrodes have 100 high activity electrode.
% each electrode we have their 100 neighbours.

% new_selected = mxw.util.electrodeSelection.networkRec( networkData );

[grp grp2] = detRedspike(networkData,5,15);






