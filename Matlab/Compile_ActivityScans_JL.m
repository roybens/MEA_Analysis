clear
close all

%% settings
% set DIV 0 date
div0 = '08/05/2022'; % format: MM/DD/YYYY
% set path to folder containing subfolders that contain h5 files
parentFolderPath = '/mnt/harddrive-2/Organoids_Mandeep_Fink_Lab';
% set output folder
outputFolderPath = '/home/jonathan/Documents/Scripts/Matlab/scrpits_output/Organoids/activityscan_outputs';

%% defines
% convert div 0 to date datatype
div0_date = datetime(div0, "InputFormat",'MM/dd/yyyy');
% create table elements
Run_ID = [];
DIV = [];
Time = [];
Chip_ID = [];
Mean_FiringRate = [];
Mean_SpikeAmplitude = [];

%% interate through ActivityScans
% get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(parentFolderPath, '**/ActivityScan/**/*.mxassay'); 
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
    pathFileActivityScan = fullfile(theFiles(k).folder);
    fprintf(1, 'Now reading %s\n',  pathFileActivityScan);
    tic
    
    %% reset recording info
    scan_runID = 0;
    scan_chipID = 0;
    scan_meanFiringRate = 0;
    scan_meanSpikeAmplitude = 0;
    hd5Date = 0; 
    scan_div = 0;
    

    %% extract dir information
    fileDirParts = strsplit(pathFileActivityScan, filesep); % split dir into elements
    scan_runID = str2double(fileDirParts{end-0}); % extract runID
    scan_chipID = str2double(fileDirParts{end-2}); % extract chipID
    % folderDate = datetime(fileDirParts{end-3},'InputFormat','yyMMdd','Format','MM-dd-yyyy'); % extract date and convert to date object
    
    %% extract hd5 information
    wellID = 1; % select which well to analyze
    % create fileManager object for the Activity Scan
    activityScanData = mxw.fileManager(pathFileActivityScan,wellID);
    
    % get the mean firing rate for each electrode
    meanFiringRate = mxw.activityMap.computeSpikeRate(activityScanData);
    % calculate the overall mean firing rate of the scan
    thrFiringRate = 0.1; % set a minimum spike rate threshold (Hz)
    scan_meanFiringRate = mean(meanFiringRate(meanFiringRate>thrFiringRate));
    
    % get th 90th percentile spike amplitude value for each electrode
    amplitude90perc = abs(mxw.activityMap.computeAmplitude90percentile(activityScanData));
    % calculate the overall mean spike amplitude of the scan
    thrAmp = 10;  % set a minimum spike amplitude threshold (uV)
    scan_meanSpikeAmplitude = mean(amplitude90perc(amplitude90perc>thrAmp));


    % get the startTime of the recordings
    hd5_time = activityScanData.fileObj.stopTime;
    try
        hd5Date = datetime(hd5_time,'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    catch
        hd5Date = datetime(hd5_time,'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
    end
    scan_div = fix(daysact(div0_date , hd5Date));
    
    %% append information to table elements
    Run_ID = [Run_ID scan_runID];
    DIV = [DIV scan_div];
    Time = [Time hd5Date];
    Chip_ID = [Chip_ID scan_chipID];
    Mean_FiringRate = [Mean_FiringRate scan_meanFiringRate];
    Mean_SpikeAmplitude = [Mean_SpikeAmplitude scan_meanSpikeAmplitude];
end

%% construct table
% convert row list to columns
Run_ID = Run_ID';
DIV = DIV';
Time = Time';
Chip_ID = Chip_ID';
Mean_FiringRate = Mean_FiringRate';
Mean_SpikeAmplitude = Mean_SpikeAmplitude';
% make table
T = table(Run_ID,DIV,Time,Chip_ID,Mean_FiringRate,Mean_SpikeAmplitude);
T = sortrows(T,"Run_ID","ascend");
writetable(T, fullfile(outputFolderPath,'Compiled_ActivityScan.csv'));
