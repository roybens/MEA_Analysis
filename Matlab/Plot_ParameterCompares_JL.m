% This script goes thruough the h5 files in a DIV,
% uses a base set of parameters and change one parameter at a time
% to plot parameter compare figures

clear
close all

% manually set path to folder containing subfolders that contain h5 files
%parentFolderPath = '/mnt/harddrive-2/ADNP/';
dataDir = '/mnt/disk15tb/jonathan/Syngap3/Syngap3/230109/';
% set path to excel file that has the reference note
%refDir = '/home/jonathan/Documents/Scripts/Python/ADNP_Notes.xlsx';
refDir = '/home/mmp/Documents/Syngap3_Notes.xlsx';
% set output folder
%opDir = '/home/jonathan/Documents/Scripts/Matlab/scripts_output/ADNP/';
opDir = '/home/mmp/Documents/script_output/Syngap3/';

% set Gaussian kernel standard deviation [s] (smoothing window)
gaussianSigma = 0.18; %0.18
% set histogram bin size [s]
binSize = 0.3;
% set minimum peak distance [s]
minPeakDistance = 0.025;
% set burst detection threshold [rms / fixed]
thresholdBurst =1.2; %1.2
% set fixed threshold;
use_fixed_threshold = false;
% Set the threshold to find the start and stop time of the bursts. (start-stop threshold)
thresholdStartStop = 0.3; %0.3


% condense the parameter set into an array to input into the the function below
base_parameters = [gaussianSigma,binSize,minPeakDistance,thresholdBurst,use_fixed_threshold,thresholdStartStop];
outDir = append(opDir, 'Network_outputs/');

Compare_NetworkParameters(dataDir,refDir, outDir, 'Gaussian', gaussianSigma,'VarParameter', [0,0.04,1], 'BaseParameters',base_parameters);
Compare_NetworkParameters(dataDir,refDir, outDir, 'BinSize', binSize,'VarParameter', [0,0.04,1], 'BaseParameters',base_parameters);
Compare_NetworkParameters(dataDir,refDir, outDir, 'Threshold', thresholdBurst, 'VarParameter', [0,0.1,3], 'BaseParameters',base_parameters);
%Compare_NetworkParameters(dataDir,refDir, outDir, 'FixedThreshold','VarParameter', [0,1,20], 'BaseParameters',base_parameters);
Compare_NetworkParameters(dataDir,refDir, outDir, 'StartStopThreshold',thresholdStartStop,'VarParameter', [0,0.04,1], 'BaseParameters',base_parameters);
Compare_NetworkParameters(dataDir,refDir, outDir, 'MinPeakDistance',minPeakDistance,'VarParameter', [0,0.08,2], 'BaseParameters',base_parameters);