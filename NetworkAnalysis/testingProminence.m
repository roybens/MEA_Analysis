
% Parameters
numNeurons = 100;  % Total number of neurons
totalTime = 100;  % Total time in seconds
baselineRate = 0.1;  % Baseline firing rate per second
burstRate = 5;  % Firing rate during bursts per second
burstTimes = [20, 40, 60, 80];  % Times at which bursts occur
burstDuration = 2;  % Duration of each burst in seconds

% Initialize spike times
spikeTimes = [];

% Generate random spikes based on baseline activity
for i = 1:numNeurons
    numSpikes = poissrnd(baselineRate * totalTime);
    spikes = rand(1, numSpikes) * totalTime;
    spikeTimes = [spikeTimes spikes];
end

% Add bursts
for t = burstTimes
    for i = 1:numNeurons
        numSpikes = poissrnd(burstRate * burstDuration);
        spikes = t + rand(1, numSpikes) * burstDuration;
        spikeTimes = [spikeTimes spikes];
    end
end

% Sorting spike times (optional, for neatness)
spikeTimes = sort(spikeTimes);

% Plot histogram of spike times to visualize
figure;
histogram(spikeTimes, 100);
xlabel('Time (s)');
ylabel('Number of spikes');
title('Histogram of Spike Times');


% Assume spikeTimes and neuronIndices are vectors containing spike times and their corresponding neuron indices
timeVector = linspace(min(spikeTimes), max(spikeTimes), 1000);  % Define a time vector for smoothing
spikeCount = histcounts(spikeTimes, timeVector);  % Count spikes in each time bin

% Apply Gaussian smoothing
sigma = 1;  % Standard deviation of the Gaussian kernel
windowSize = 5 * sigma;  % Define the size of the Gaussian window
gaussianWindow = normpdf(-windowSize:windowSize, 0, sigma);
smoothedRate = conv(spikeCount, gaussianWindow, 'same');  % Convolve spike count with the Gaussian window

% Plot smoothed firing rate
figure;
plot(timeVector(1:end-1), smoothedRate);
xlabel('Time');
ylabel('Smoothed Firing Rate');
title('Gaussian Smoothed Spike Rate');


% Parameters for peak detection
minProminence = 0.15 * max(smoothedRate);  % Set minimum prominence as a fraction of the maximum rate

[pks, locs, widths, proms] = findpeaks(smoothedRate, timeVector(1:end-1), ...
                                      'MinPeakProminence', 0);

% Plot peaks on the smoothed data
hold on;
plot(locs, pks, 'ro');
legend('Smoothed Rate', 'Significant Bursts');
title('Identified Bursts with Prominence Criterion');
