%%
% developed by Elena, Jonathan and Mandar.
%%


function [values] = computeNetworkStatsModified(networkActivityVectorFR, networkActivityVectortime, varargin)

p = inputParser;

p.addRequired('networkActivityVectorFR');
p.addRequired('networkActivityVectortime');
p.addParameter('ThresholdMethod', 'RMS'); % New parameter to specify the threshold method

p.addParameter('Threshold', 1.2);
p.addParameter('MinPeakDistance', 0.8);
p.addParameter('MinPeakProminence', 1.5); % New parameter for minimum peak prominence

p.parse(networkActivityVectorFR, networkActivityVectortime, varargin{:});
args = p.Results;

peakParams =[];
if ~isempty(args.MinPeakDistance)
    peakParams = [peakParams, {'MinPeakDistance', args.MinPeakDistance}];
end
if ~isempty(args.MinPeakProminence)
    peakParams = [peakParams, {'MinPeakProminence', args.MinPeakProminence}];
end

% Determine the appropriate threshold based on the specified method

switch args.ThresholdMethod
    case 'Fixed'
        if isempty(args.FixedThreshold)
            error('FixedThreshold must be provided for fixed threshold method');
        end
        threshold = args.FixedThreshold;
        [maxAmplitudesValues, maxAmplitudesTimes] = findpeaks(args.networkActivityVectorFR, ...
            args.networkActivityVectortime,'MinPeakHeight',threshold, peakParams{:});
    case 'RMS'
        if isempty(args.Threshold)
            error('Threshold must be provided for RMS threshold method');
        end
        rmsFiringRate = mxw.util.rms(args.networkActivityVectorFR);
        threshold = args.Threshold * rmsFiringRate;
        [maxAmplitudesValues, maxAmplitudesTimes] = findpeaks(args.networkActivityVectorFR, ...
            args.networkActivityVectortime,'MinPeakHeight',threshold, peakParams{:});
    case 'Adaptive'
        % Use findpeaks to get initial peaks
        [initialPeaks, initialLocs] = findpeaks(args.networkActivityVectorFR, ...
            args.networkActivityVectortime, peakParams{:});
        % threshold = movmean(args.networkActivityVectorFR, 30) + 2 * movstd(args.networkActivityVectorFR, 30);   % 10 is hardcoded
        % % Filter peaks using the adaptive threshold
        % [~, locIndices] = ismember(initialLocs, networkActivityVectortime);
        % validPeaks = initialPeaks > threshold(locIndices);
        % maxAmplitudesValues = initialPeaks(validPeaks);
        % maxAmplitudesTimes = initialLocs(validPeaks);
        maxAmplitudesValues = initialPeaks;
        maxAmplitudesTimes = initialLocs;
        threshold =zeros(1,3000);
    
    otherwise
        error('Unknown ThresholdMethod. Choose from "fixed", "rms", or "adaptive".');
end


values.maxAmplitudesValues = maxAmplitudesValues;
maxAmplitudeTimeDiff = diff(maxAmplitudesTimes);
values.maxAmplitudeTimeDiff = maxAmplitudeTimeDiff;
values.maxAmplitudesTimes = maxAmplitudesTimes;
values.maxAmplitudesValues = maxAmplitudesValues;
values.thresholdFunction = args.ThresholdMethod;
values.threshold = threshold;

end