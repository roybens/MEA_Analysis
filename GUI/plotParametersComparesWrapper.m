% This script goes thruough the h5 files in a DIV,
% uses a base set of parameters and change one parameter at a time
% to plot parameter compare figures
% Developers: Jonathan Lin, Mandar Patil and Tim Fenton


function [success] = plotParametersComparesWrapper(data)

success = false;


% Unpack the data structure
dataDir = data.dataDir;
refDir = data.refDir;
opDir = data.opDir;
gaussianSigma = data.gaussianSigma;
binSize = data.binSize;
minPeakDistance = data.minPeakDistance;
thresholdBurst = data.thresholdBurst;
use_fixed_threshold = false;
thresholdStartStop = data.thresholdStartStop; 




% condense the parameter set into an array to input into the the function below
base_parameters = [gaussianSigma,binSize,minPeakDistance,thresholdBurst,use_fixed_threshold,thresholdStartStop];
outDir = append(opDir, '/Network_outputs/');

% Define the parameters and their corresponding values
parameters = {'Gaussian', 'BinSize', 'Threshold', 'StartStopThreshold', 'MinPeakDistance'};
parameterValues = {gaussianSigma, binSize, thresholdBurst, thresholdStartStop, minPeakDistance};

fig = data.fig;
% Create the progress dialog
d = uiprogressdlg(fig, 'Title', 'Compiling files', 'Message', '1', 'Cancelable', 'on');
drawnow;

% Loop through the parameters and call Compare_NetworkParameters function
for i = 1:numel(parameters)
    % Update the progress dialog message
    d.Message = parameters{i};
    drawnow;
    % Check if the Cancel button is pressed
    if d.CancelRequested
        disp('Operation canceled by user.');
        break;
    end
    % Call Compare_NetworkParameters function with the current parameter
    Compare_NetworkParameters(dataDir, refDir, outDir, parameters{i}, parameterValues{i}, 'VarParameter', [0, 0.04, 1], 'BaseParameters', base_parameters);
    
    % Update the progress value
    d.Value = i / numel(parameters);
    drawnow;
end
% Check if the Cancel button is pressed
if d.CancelRequested
    success = false;
    % Close the progress dialog
    delete(d);
    error('Operation canceled by user.');
end
% Close the progress dialog
delete(d);
success = true;



