

close all;
clear all;

pathFileNetwork =  '/mnt/disk15tb/paula/PAULA/Trace_20230516_16_52_27.raw.h5';

wellID = 1;
networkData = mxw.fileManager(pathFileNetwork,wellID);
relativeSpikeTimes = mxw.util.computeRelativeSpikeTimes(networkData);


% select channels of interest
InterestChannels = [875, 876, 877];

% select time period of interest (stimulus onset,stimulus offset)
periodStart = 30;
periodEnd = 90;

% binning
binSize = 0.1;  % The bin size in seconds
binEdgesTime = periodStart:binSize:periodEnd; % Bin Edges

% % determine spike counts that are  binned
% spikeCounts_IC = histcounts(relativeSpikeTimes.time(ismember(relativeSpikeTimes.channel, InterestChannels) & ...
%     relativeSpikeTimes.time >= periodStart & relativeSpikeTimes.time < periodEnd), binEdgesTime);
spikeCounts = zeros(numel(InterestChannels),numel(binEdgesTime)-1);
for i = 1:numel(InterestChannels)
    spikeCounts(i,:) = histcounts(relativeSpikeTimes.time(ismember(relativeSpikeTimes.channel, InterestChannels(i)) & ...
     relativeSpikeTimes.time >= periodStart & relativeSpikeTimes.time < periodEnd), binEdgesTime);
    fr(i,:)=spikeCounts(i,:)/binSize;
end
% Calculate firing rate 
%firingRate = spikeCounts_IC / binSize;
% determine channels from time range
uniqueChannels = unique(relativeSpikeTimes.channel(ismember(relativeSpikeTimes.channel, InterestChannels) & ...
    relativeSpikeTimes.time >= periodStart & relativeSpikeTimes.time < periodEnd));
% Define colors for the bars
numChannels = numel(uniqueChannels);
colorMap = prism(numChannels);  % Or specify your own custom color map
barColors = colorMap(1:numChannels, :);
% Create a bar graph of the firing rate
bar(binEdgesTime(1:end-1), [fr(1,:);fr(2,:);fr(3,:)])% 'hist')
colormap(barColors)
xlabel('Time (s)')
ylabel('Average Spikes per Binsize (spikes/s)')
title('Firing Rate for Each Bin')


% Add channel legends to the graph
legendLabels = cell(1, numel(uniqueChannels));
for i = 1:numel(uniqueChannels)
    channel = uniqueChannels(i);
    legendLabels{i} = sprintf('Channel %d', channel);
end
legend(legendLabels, 'Location', 'northeast')

% Adjust the appearance of the plot
xlim([binEdgesTime(1), binEdgesTime(end-1)])
ylim([0, max(fr(:))*1.1])
%grid on