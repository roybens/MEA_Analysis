function [ values ] = computeNetworkStats_JL( networkActivityVector, varargin )

p = inputParser;

p.addRequired('struct', @(x) isstruct(x));
p.addParameter('Threshold', []);
p.addParameter('FixedThreshold', []);
p.addParameter('MinPeakDistance', []);


p.parse(networkActivityVector, varargin{:});
args = p.Results;

rmsFiringRate = mxw.util.rms(networkActivityVector.firingRate);

if ~isempty(args.MinPeakDistance)
    if ~isempty(args.Threshold)
        [maxAmplitudesValues, maxAmplitudesTimes] = findpeaks(networkActivityVector.firingRate,...
            networkActivityVector.time, 'MinPeakHeight', args.Threshold * rmsFiringRate, 'MinPeakDistance', args.MinPeakDistance);
    elseif ~isempty(args.FixedThreshold)
        [maxAmplitudesValues, maxAmplitudesTimes] = findpeaks(networkActivityVector.firingRate,...
            networkActivityVector.time, 'MinPeakHeight', args.FixedThreshold, 'MinPeakDistance', args.MinPeakDistance);
    else
        error('Please set a threshold')
    end
else
    if ~isempty(args.Threshold)
        [maxAmplitudesValues, maxAmplitudesTimes] = findpeaks(networkActivityVector.firingRate,...
            networkActivityVector.time, 'MinPeakHeight', args.Threshold * rmsFiringRate);
    elseif ~isempty(args.FixedThreshold)
        [maxAmplitudesValues, maxAmplitudesTimes] = findpeaks(networkActivityVector.firingRate,...
            networkActivityVector.time, 'MinPeakHeight', args.FixedThreshold);
    else
        error('Please set a threshold')
    end
end

values.maxAmplitudesValues = maxAmplitudesValues;

maxAmplitudeTimeDiff = diff(maxAmplitudesTimes);
values.maxAmplitudeTimeDiff = maxAmplitudeTimeDiff;
values.maxAmplitudesTimes = maxAmplitudesTimes;
values.maxAmplitudesValues = maxAmplitudesValues;


end