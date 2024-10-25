function [ values ] = computeNetworkStatsNew_JL( networkActivityVectorFR,networkActivityVectortime, varargin )

p = inputParser;

p.addRequired('networkActivityVectorFR');
p.addRequired('networkActivityVectortime');
p.addParameter('ThresholdMethod', 'adaptive');
p.addParameter('Threshold', []);
p.addParameter('MinPeakDistance', []);


p.parse(networkActivityVectorFR, networkActivityVectortime,varargin{:});
args = p.Results;

rmsFiringRate = mxw.util.rms(args.networkActivityVectorFR);

if ~isempty(args.MinPeakDistance)
    if ~isempty(args.Threshold)
        % [maxAmplitudesValues, maxAmplitudesTimes] = findpeaks(args.networkActivityVectorFR,...
        %     args.networkActivityVectortime, 'MinPeakHeight', args.Threshold * rmsFiringRate, 'MinPeakDistance', args.MinPeakDistance);
        
        [maxAmplitudesValues, maxAmplitudesTimes] = findpeaks(args.networkActivityVectorFR,...
             args.networkActivityVectortime, 'MinPeakProminence',1.5 , 'MinPeakDistance', args.MinPeakDistance,'MinPeakHeight', args.Threshold * rmsFiringRate);
        
    elseif ~isempty(args.FixedThreshold)
        [maxAmplitudesValues, maxAmplitudesTimes] = findpeaks(args.networkActivityVectorFR,...
            args.networkActivityVectortime, 'MinPeakHeight', args.FixedThreshold, 'MinPeakDistance', args.MinPeakDistance);
    else
        error('Please set a threshold')
    end
else
    if ~isempty(args.Threshold)
        [maxAmplitudesValues, maxAmplitudesTimes] = findpeaks(args.networkActivityVectorFR,...
            args.networkActivityVectortime, 'MinPeakHeight', args.Threshold * rmsFiringRate);
    elseif ~isempty(args.FixedThreshold)
        [maxAmplitudesValues, maxAmplitudesTimes] = findpeaks(args.networkActivityVectorFR,...
            args.networkActivityVectortime, 'MinPeakHeight', args.FixedThreshold);
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