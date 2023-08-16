
clear
close all


% %pathFileNetwork = '/mnt/disk15tb/mmpatil/1000_Electrodes_1hz/Trace_20230427_16_22_05_18908.raw.h5';
%pathFileNetwork = '/mnt/disk15tb/mmpatil/Spikesorting/Data/230327fourblocks/Trace_20230327_16_35_34.raw.h5';
%pathFileNetwork = '/mnt/ben-shalom_nas/irc/home/mxwbio/Desktop/Mandar/april 28/230428/16757/Network/000027 firing rate/data.raw.h5';
pathFileNetwork = '/mnt/harddrive-2/ADNP/ADNP/230516/16848/Network/000060/data.raw.h5';
networkData = mxw.fileManager(pathFileNetwork);
spikeRate = mxw.activityMap.computeSpikeRate(networkData);
amplitude90perc = abs(mxw.activityMap.computeAmplitude90percentile(networkData));
% rawFrames = double(networkData.rawMap.spikes.frameno);
% rawSpikes = double(networkData.rawMap.spikes.amplitude);
srthreshold = prctile(spikeRate,90);
indicesSR = find(spikeRate>=srthreshold);
requiredElectrodesSR = networkData.processedMap.electrode(indices);
ampThreshold = prctile(amplitude90perc,90);
indices = find(amplitude90perc >= ampThreshold);
requiredElectrodesAmp = networkData.processedMap.electrode(indices);
    % sampling frequency
fsNetwork = networkData.fileObj.samplingFreq;

relativeSpiketime = mxw.util.computeRelativeSpikeTimes(networkData);
N =10;
downsample_spiketimes  = downsample(relativeSpiketime.time,N);
spiketimes = downsample_spiketimes * 1e4;
spiketimes = round(spiketimes);
max_time = max(spiketimes);
max_channel = max(relativeSpiketime.channel);
channel = downsample(relativeSpiketime.channel,N);
spike_matrix = zeros(max_channel, max_time);

for i  = 1:length(channel)
    chan = channel(i);
    time = spiketimes(i);
    spike_matrix(chan,time) = 1;
end

save('spiketimes_firingrate.mat','spike_matrix','-v7.3');


