clear all
close all

% Enable offscreen rendering
set(0, 'DefaultFigureVisible', 'off');
opengl('save', 'software');

% Define path for data file
pathFileNetwork = '/mnt/Vol20tb1/PrimaryNeuronData/MaxTwo/KCNT1_T5_C1_04132024/KCNT1_T5_C1_04132024/240430/M08029/Network/000039/data.raw.h5';

% Load Maxwell data
networkData = mxw.fileManager(pathFileNetwork, 5);

% Compute relative spike times
relativeSpikeTimes = mxw.util.computeRelativeSpikeTimes(networkData);

% Debug: Plot raster plot
figure('Color', 'w', 'Position', [100, 100, 800, 400]);
scatter(relativeSpikeTimes.time, relativeSpikeTimes.channel, 5, 'filled');
ylabel('Channel');
xlabel('Time [s]');
title('Raster Plot');
grid on;

% Save raster plot
print(gcf, 'RasterPlot.png', '-dpng', '-r300');
close;

% ISI Calculation for 10th Spike Ahead
ts = (double(networkData.fileObj.spikes.frameno) - double(networkData.fileObj.firstFrameNum)) / networkData.fileObj.samplingFreq;
ch = networkData.fileObj.spikes.channel;

% Validate time and channel data
validIdx = ts > 0;
ts = ts(validIdx);
ch = ch(validIdx);

% Get unique channels
uniqueChannels = unique(ch);
ISI_perElectrode = struct();
ISI_sec_perElectrode = struct();
spikeStep = 10; % Step for the 10th spike ahead
outputPDF = 'ISI_Distributions_with_Bursts.pdf'; % PDF to store histograms

% Loop over each channel to calculate ISIs
for i = 1:length(uniqueChannels)
    currentChannel = uniqueChannels(i);
    spikeTimes = sort(ts(ch == currentChannel));
    
    % Calculate ISI for the 10th spike ahead
    if length(spikeTimes) > spikeStep
        ISIs = spikeTimes(1+spikeStep:end) - spikeTimes(1:end-spikeStep);
        logISIs = log10(ISIs);
        ISI_perElectrode.(sprintf('Channel_%d', currentChannel)) = logISIs;

        % Step 1: Calculate Histogram
        [counts, edges] = histcounts(logISIs, 'Normalization', 'probability', 'BinWidth', 0.1);
        centers = edges(1:end-1) + diff(edges) / 2;

        % Step 2: Smooth the Histogram
        if length(counts) >= 3
            smoothedCounts = smoothdata(counts, 'gaussian', 1);

            % Plot the Smoothed Histogram
            figure('Color', 'w', 'Position', [100, 100, 800, 400]);
            plot(centers, smoothedCounts, 'LineWidth', 2);
            xlabel('log(ISI)');
            ylabel('Smoothed Probability');
            title(['Smoothed ISI Histogram for Channel ', num2str(currentChannel)]);
            grid on;

            % Step 3: Detect Peaks
            if length(smoothedCounts) >= 3
                [peakValues, peakIndices] = findpeaks(smoothedCounts);

                % Plot Peaks
                hold on;
                plot(centers(peakIndices), smoothedCounts(peakIndices), 'or', 'MarkerSize', 8, 'DisplayName', 'Peaks');
                legend('Smoothed Histogram', 'Detected Peaks');
                hold off;

                % Step 4: Find Valleys Between Peaks
                if length(peakIndices) >= 2
                    valleys = []; % To store valley indices
                    for j = 1:length(peakIndices) - 1
                        % Range between two peaks
                        rangeStart = peakIndices(j);
                        rangeEnd = peakIndices(j + 1);

                        % Find the minimum in this range
                        [valleyValue, valleyIdx] = min(smoothedCounts(rangeStart:rangeEnd));
                        valleys = [valleys; rangeStart + valleyIdx - 1]; % Adjust index for full array
                    end

                    % Plot Valleys
                    hold on;
                    plot(centers(valleys), smoothedCounts(valleys), 'xg', 'MarkerSize', 10, 'DisplayName', 'Valleys');
                    legend('Smoothed Histogram', 'Detected Peaks', 'Detected Valleys');
                    hold off;

                    % Use the first valley as the burst threshold
                    burstThreshold = 10^centers(valleys(1));
                    disp(['Channel ', num2str(currentChannel), ' Burst Threshold (sec): ', num2str(burstThreshold)]);

                    % Separate ISIs into bursts and non-bursts
                    intraBurst_ISIs = ISIs(ISIs <= burstThreshold);
                    interBurst_ISIs = ISIs(ISIs > burstThreshold);

                    % Display counts
                    disp(['Intra-Burst ISIs: ', num2str(length(intraBurst_ISIs))]);
                    disp(['Inter-Burst ISIs: ', num2str(length(interBurst_ISIs))]);

                    % Save burst threshold
                    ISI_sec_perElectrode.(sprintf('Channel_%d', currentChannel)) = burstThreshold;
                else
                    disp(['Channel ', num2str(currentChannel), ' does not have enough peaks to determine valleys.']);
                    ISI_sec_perElectrode.(sprintf('Channel_%d', currentChannel)) = [];
                end
            else
                disp(['Channel ', num2str(currentChannel), ' skipped: insufficient data for peak detection.']);
                ISI_sec_perElectrode.(sprintf('Channel_%d', currentChannel)) = [];
            end
        else
            disp(['Channel ', num2str(currentChannel), ' skipped: insufficient data for histogram smoothing.']);
            ISI_sec_perElectrode.(sprintf('Channel_%d', currentChannel)) = [];
        end

        % Save each plot as an image file
        print(gcf, ['Channel_', num2str(currentChannel), '_ISI_Histogram.png'], '-dpng', '-r300');
        close; % Close the figure after saving
    else
        ISI_perElectrode.(sprintf('Channel_%d', currentChannel)) = [];
        ISI_sec_perElectrode.(sprintf('Channel_%d', currentChannel)) = [];
        disp(['Channel ', num2str(currentChannel), ' skipped: insufficient spikes for 10th spike calculation.']);
    end
end

% Summary of Results
disp('Burst Thresholds (seconds) for all channels:');
disp(ISI_sec_perElectrode);