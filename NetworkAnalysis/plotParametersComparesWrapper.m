% This script goes thruough the h5 files in a DIV,
% uses a base set of parameters and change one parameter at a time
% to plot parameter compare figures
% Developers: Jonathan Lin, Mandar Patil and Tim Fenton


function [success] = plotParametersComparesWrapper(data)

success = false;

% Validate that each field is a comma-separated array with three elements

fields = {'gaussianSigma', 'binSize', 'minPeakDistance', 'thresholdBurst','minPeakProminence'};
% if any(~isfield(data, fields) | cellfun(@(f) length(data.(f)), fields) ~= 3)
%     errordlg('All fields in data must be arrays with (default, start, end) separated by commas.');
%     return;
% end
% Unpack the data structure
dataDir = data.dataDir;
refDir = data.refDir;
opDir = data.opDir;
gaussianSigma = data.gaussianSigma;
binSize = data.binSize;
minPeakDistance = data.minPeakDistance;
minPeakProminence = data.minProminience;
thresholdBurst = data.thresholdBurst;
thresholdMethod = data.thresholdMethod;
thresholdStartStop = data.thresholdStartStop; 
logFile= data.logFile;



% condense the parameter set into an array to input into the the function below
base_parameters = {gaussianSigma,binSize,minPeakDistance,thresholdBurst,thresholdMethod,thresholdStartStop,minPeakProminence};
outDir = append(opDir, '/Network_outputs/');

% Define the parameters and their corresponding values
parameters = {'Gaussian', 'BinSize', 'Threshold', 'StartStopThreshold', 'MinPeakDistance','MinPeakProminence'};
parameterValues = {gaussianSigma, binSize, thresholdBurst, thresholdStartStop, minPeakDistance,minPeakProminence};

fig = data.fig;
% Create the progress dialog
d = uiprogressdlg(fig, 'Title', 'Compiling files', 'Message', '1', 'Cancelable', 'on');
drawnow;


try
    % Pre-calculate total number of parameters
    numParameters = numel(parameters);
    
    % Read JSON file
    fid = fopen('expParameterSetting.json');
    if fid == -1
        error('Cannot open params.json file');
    end
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    
     % Decode JSON file
    params = jsondecode(str);

    % Specify the exact number of cores to use
    numCores =  numParameters;  % Adjust this number based on your needs and resource availability
    
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
    % Execute computation in parallel
    parfor i = 1:numParameters

            if isfield(params, parameters{i})
                varParams = params.(parameters{i})';
            else
                error('Parameter name "%s" not found in the JSON file', parameters{i});
            end
        % Call the Compare_NetworkParameters function with current parameter
        Compare_NetworkParameters(dataDir, refDir, outDir, parameters{i}, parameterValues{i}, 'VarParameter', varParams, 'BaseParameters', base_parameters);
    end

    % If completed successfully
    d.Message = 'Completion successful!';
    d.Value = 1.0;
    pause(1); % Display completion message briefly
catch ME
    if strcmp(ME.identifier, 'MATLAB:OperationCancelled')
        disp('Operation canceled by user.');
        d.Message = 'Operation canceled by user.';
        d.Value = 0; % Reflect the cancellation in the progress dialog
        pause(1); % Display cancel message briefly
        delete(d);
        error('Operation canceled by user.');
    else
        fprintf(1,"%s",ME.message);
        rethrow(ME);
    end
end

% Close the progress dialog after a successful operation
delete(d);
success = true;


delete(gcp('nocreate'));
%Find all existing jobs
allJobs = findJob(parcluster);

% Delete each job
for i = 1:length(allJobs)
    delete(allJobs(i));
end

fprintf(1,'Completion successful!');
