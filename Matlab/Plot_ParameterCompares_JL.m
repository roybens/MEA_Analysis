% This script goes thruough the h5 files in a DIV,
% uses a base set of parameters and change one parameter at a time
% to plot parameter compare figures

clear
close all

% set path to folder containing subfolders that contain h5 files
dataDir = '/mnt/harddrive-2/ADNP/ADNP/230510';
% set path to excel file that has the reference note
refDir = '/home/jonathan/Documents/Scripts/Python/ADNP_Notes.xlsx';
% set path to an output folder
opDir = '/home/jonathan/Documents/Scripts/Matlab/scripts_output/ADNP/';

% set Gaussian kernel standard deviation [s] (smoothing window)
gaussianSigma = 0.18;
% set histogram bin size [s]
binSize = 0.1;
% set minimum peak distance [s]
minPeakDistance = 1.0;
% set burst detection threshold [rms firing rate]
thresholdBurst = 0.85;
% set fixed threshold;
use_fixed_threshold = false;
% Set the threshold to find the start and stop time of the bursts. (start-stop threshold)
thresholdStartStop = 0.3;

% condense the parameter set into an array to input into the the function below
base_parameters = [gaussianSigma,binSize,minPeakDistance,thresholdBurst,use_fixed_threshold,thresholdStartStop];
outDir = append(opDir, 'Network_outputs/');

Compare_NetworkParameters(dataDir,refDir, outDir, 'Gaussian', 'VarParameter', [0,0.04,1], 'BaseParameters',base_parameters);
Compare_NetworkParameters(dataDir,refDir, outDir, 'BinSize', 'VarParameter', [0,0.04,1], 'BaseParameters',base_parameters);
Compare_NetworkParameters(dataDir,refDir, outDir, 'Threshold', 'VarParameter', [0,0.08,2], 'BaseParameters',base_parameters);
Compare_NetworkParameters(dataDir,refDir, outDir, 'FixedThreshold','VarParameter', [0,1,20], 'BaseParameters',base_parameters);
Compare_NetworkParameters(dataDir,refDir, outDir, 'StartStopThreshold','VarParameter', [0,0.04,1], 'BaseParameters',base_parameters);
Compare_NetworkParameters(dataDir,refDir, outDir, 'MinPeakDistance','VarParameter', [0,0.08,2], 'BaseParameters',base_parameters);