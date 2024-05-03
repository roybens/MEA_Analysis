clear
close all




project_name = 'shank3';

% set DIV 0 date
div0 = '11/15/2023';

dataDir = '/mnt/disk20tb/PrimaryNeuronData/Maxtwo/SHANK3_1/SHANK3_1/231205/';
refDir ='/mnt/disk15tb/paula/Main_DA_Projects/Ref_Files/SHANK3_1 data/SHANK3_1_ref.xlsx';

opDir =  

% set Gaussian kernel standard deviation [s] (smoothing window)
gaussianSigma = 0.14; %0.18
% set histogram bin size [s]
binSize = 0.1;
% set minimum peak distance [s]
minPeakDistance = 1.0;
% set burst detection threshold [rms / fixed]
thresholdBurst =1.0; %1.2

% Set the threshold to find the start and stop time of the bursts. (start-stop threshold)
thresholdStartStop = 0.3; %0.3
% set fixed threshold;
use_fixed_threshold = 0;
% Set Threshold function for later use
threshold_fn = 'Threshold';

if use_fixed_threshold
    threshold_fn = 'FixedThreshold';
end

% condense the parameter set into an array to input into the the function below
base_parameters = [gaussianSigma,binSize,minPeakDistance,thresholdBurst,use_fixed_threshold,thresholdStartStop];
outDir = append(opDir, '/Network_outputs/');

% Define the parameters and their corresponding values
parameters = {'Gaussian', 'BinSize', 'Threshold', 'StartStopThreshold', 'MinPeakDistance'};
parameterValues = {gaussianSigma, binSize, thresholdBurst, thresholdStartStop, minPeakDistance};

% Define variable parameters for each main parameter
variableParams = {
    [0, 0.02, 1],  ... Gaussian sigma variations
    [0, 0.02, 1],  ... Bin size variations
    [0, 0.02, 2],   ...Threshold burst variations
    [0, 0.02, 1],    ... Start-stop threshold variations
    [0, 0.1, 1]     ... Min peak distance variations
};
% Specify the exact number of cores to use
numCores = 10;  % Adjust this number based on your needs and resource availability

% Initialize or modify the existing parallel pool
currentPool = gcp('nocreate');  % Check for existing parallel pool

if isempty(currentPool)
    poolObj = parpool(numCores);  % Create a new pool if none exists
elseif currentPool.NumWorkers ~= numCores
    delete(currentPool);  % Close the existing pool if it does not match the desired size
    poolObj = parpool(numCores);  % Create a new pool with the correct number of workers
else
    poolObj = currentPool;  % Use the existing pool if it already matches the desired number of workers
end
% Loop through the parameters and call Compare_NetworkParameters function
parfor i = 1:numel(parameters)


    % Call Compare_NetworkParameters function with the current parameter
    Compare_NetworkParameters(dataDir, refDir, outDir, parameters{i}, parameterValues{i}, 'VarParameter',variableParams{i}, 'BaseParameters', base_parameters);
    

end

delete(gcp('nocreate'));







