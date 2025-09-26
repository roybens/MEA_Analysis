
clear all
close all


% %pathFileNetwork = '/mnt/disk15tb/mmpatil/1000_Electrodes_1hz/Trace_20230427_16_22_05_18908.raw.h5';
%pathFileNetwork = '/mnt/disk15tb/mmpatil/Spikesorting/Data/230327fourblocks/Trace_20230327_16_35_34.raw.h5';
%pathFileNetwork = '/mnt/ben-shalom_nas/irc/home/mxwbio/Desktop/Mandar/april 28/230428/16757/Network/000027 firing rate/data.raw.h5';
pathFileNetwork = '/mnt/disk20tb/PrimaryNeuronData/Maxtwo/CDKL5-R59X_SingleBrainPerChip_11252024_PS/CDKL5-R59X_SingleBrainPerChip_11252024_PS/241208/M07038/Network/000027/data.raw.h5';
% Load Maxwell data
networkData = mxw.fileManager(pathFileNetwork);

% Get spike times (ts) in seconds, and corresponding channels (ch)
tsAll = (double(networkData.fileObj.spikes.frameno) - ...
         double(networkData.fileObj.firstFrameNum)) / networkData.fileObj.samplingFreq;
chAll = networkData.fileObj.spikes.channel;

% Exclude any spikes before t = 0 (if needed)
validIdx = tsAll > 0;
tsAll   = tsAll(validIdx);
chAll   = chAll(validIdx);

%% --------------------- Group Spikes by Channel ---------------------------
[uniqueChannels, ~, channelIndices] = unique(chAll);
spikeTimesByChannel = accumarray(channelIndices, tsAll, [], @(x){sort(x)});

%% --------------------- Parameters for n-step ISI -------------------------
nValues = [1, 10];  % We want n=1 and n=10

% Storage for intervals across all channels for n=1 and n=10
ISI_n_all = cell(numel(nValues), 1);

%% --------------------- Compute n-step intervals --------------------------
for chanIdx = 1:length(spikeTimesByChannel)
    spikeTimes = spikeTimesByChannel{chanIdx};
    
    % Only proceed if we have at least 2 spikes for n=1,
    % and at least 11 spikes for n=10, etc.
    numSpikes = length(spikeTimes);
    
    % For each requested n, compute intervals for this channel
    for nv = 1:numel(nValues)
        nStep = nValues(nv);
        
        % Check if we have enough spikes to compute n-step intervals
        if numSpikes >= (nStep + 1)
            % n-step intervals: time difference between spike(i) and spike(i+nStep)
            intervals = spikeTimes((1 + nStep):end) - spikeTimes(1:(end - nStep));
            
            % Store intervals in a common cell array for all channels
            ISI_n_all{nv} = [ISI_n_all{nv}; intervals]; %#ok<AGROW>
        end
    end
end

%% --------------------- Log Transform and Filter --------------------------
% For each n, convert intervals to log scale, filter out any zero/negative
logISI = cell(size(ISI_n_all));
for nv = 1:numel(nValues)
    validIntervals = ISI_n_all{nv}(ISI_n_all{nv} > 0);  % remove non-positive intervals if any
    logISI{nv} = log(validIntervals);  % natural log
    % Alternatively, use log10(...) if you prefer base 10
end

%% --------------------- Plot Distributions on One Figure ------------------
figure('Name','n-step ISI Distributions','Color','w');
hold on;

% Number of bins for the histogram
numBins = 50;

% Plot each distribution with a different color/opacity
colors = lines(numel(nValues));  % default color set
legendEntries = cell(numel(nValues),1);

for nv = 1:numel(nValues)
    histogram(logISI{nv}, numBins, ...
        'Normalization', 'probability', ...
        'FaceAlpha', 0.4, ...
        'EdgeColor', 'none', ...
        'FaceColor', colors(nv,:));
    legendEntries{nv} = ['n = ', num2str(nValues(nv))];
end

legend(legendEntries, 'Location', 'best');
xlabel('log(ISI) [natural log of seconds]');
ylabel('Probability');
title('Comparison of n-step ISI Distributions (n=1 vs. n=10) Across Channels');
grid on;

%% --------------------- Export Figure to PDF ------------------------------
outputPDF = 'nStep_ISI_Distributions.pdf';
exportgraphics(gcf, outputPDF, 'ContentType','vector');
disp(['Figure exported to ', outputPDF]);

% % Sort spikes by channel and time
% [sortedCh, sortIdx] = sort(ch);
% sortedTs = ts(sortIdx);
% 
% % Get channel transitions
% chTransitions = find(diff(sortedCh));
% chTransitions = [0; chTransitions; length(sortedCh)];
% 
% % Preallocate result array
% logIsiN10 = zeros(size(ts));
% 
% % Calculate log ISI N-10 using vectorized operations
% for i = 1:length(chTransitions)-1
%     startIdx = chTransitions(i) + 1;
%     endIdx = chTransitions(i+1);
% 
%     channelTs = sortedTs(startIdx:endIdx);
%     if length(channelTs) >= 11
%         isiN10_temp = channelTs(11:end) - channelTs(1:end-10);
%         % Calculate log of ISI N-10, add small constant to avoid log(0)
%         logIsiN10(sortIdx(startIdx+10:endIdx)) = log10(isiN10_temp + eps);
%     end
% end
% 
% % Restore original order
% logIsiN10 = logIsiN10(sortIdx);
% 
% % Remove -Inf values (if any)
% logIsiN10(isinf(logIsiN10)) = NaN;
% 
% % Create figure
% figure('Position', [100 100 800 600]);
% 
% % Plot histogram
% histogram(logIsiN10, 100, 'Normalization', 'probability', 'EdgeColor', 'none');
% xlabel('log_{10}(ISI N-10) [s]');
% ylabel('Probability');
% title('Distribution of log ISI N-10');
% grid on;
% 
% % Add mean and median lines
% hold on;
% xline(mean(logIsiN10, 'omitnan'), 'r--', 'LineWidth', 2, 'Label', sprintf('Mean: %.2f', mean(logIsiN10, 'omitnan')));
% xline(median(logIsiN10, 'omitnan'), 'g--', 'LineWidth', 2, 'Label', sprintf('Median: %.2f', median(logIsiN10, 'omitnan')));
% hold off;
% 
% % Make it pretty
% set(gca, 'FontSize', 12);
% box on;
% 
% % Add text with basic statistics
% stats = sprintf('n = %d\nMean = %.2f\nMedian = %.2f\nStd = %.2f', ...
%     sum(~isnan(logIsiN10)), ...
%     mean(logIsiN10, 'omitnan'), ...
%     median(logIsiN10, 'omitnan'), ...
%     std(logIsiN10, 'omitnan'));
% annotation('textbox', [0.7 0.7 0.2 0.2], 'String', stats, 'EdgeColor', 'none', 'FontSize', 10);
% 
% 
% 
% % amplitude90perc = abs(mxw.activityMap.computeAmplitude90percentile(networkData));
% % % rawFrames = double(networkData.rawMap.spikes.frameno);
% % % rawSpikes = double(networkData.rawMap.spikes.amplitude);
% % srthreshold = prctile(spikeRate,90);
% % indicesSR = find(spikeRate>=srthreshold);
% % requiredElectrodesSR = networkData.processedMap.electrode(indices);
% % ampThreshold = prctile(amplitude90perc,90);
% % indices = find(amplitude90perc >= ampThreshold);
% % requiredElectrodesAmp = networkData.processedMap.electrode(indices);
% %     % sampling frequency
% % fsNetwork = networkData.fileObj.samplingFreq;
% % 
% % relativeSpiketime = mxw.util.computeRelativeSpikeTimes(networkData);
% % % N =10;
% % downsample_spiketimes  = downsample(relativeSpiketime.time,N);
% % spiketimes = downsample_spiketimes * 1e4;
% % spiketimes = round(spiketimes);
% % max_time = max(spiketimes);
% % max_channel = max(relativeSpiketime.channel);
% % channel = downsample(relativeSpiketime.channel,N);
% % spike_matrix = zeros(max_channel, max_time);
% % 
% % for i  = 1:length(channel)
% %     chan = channel(i);
% %     time = spiketimes(i);
% %     spike_matrix(chan,time) = 1;
% % end

%save('spiketimes_firingrate.mat','spike_matrix','-v7.3');


