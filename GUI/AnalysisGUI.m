function mainFunction
   % Create a uifigure window
fig = uifigure('Name', 'Process Activity /Network', 'Position', [100, 100, 700, 500], 'SizeChangedFcn', @updateLayout);

% Create a project name label and edit field
projectNameLabel = uilabel(fig, 'Text', 'Project Name:', 'Position', [50, 350, 100, 20]);
projectNameEdit = uieditfield(fig, 'Position', [150, 350, 200, 25]);

% Create a DIV 0 date label and edit field
div0DateLabel = uilabel(fig, 'Text', 'DIV 0 Date:', 'Position', [50, 310, 100, 20]);
div0DateEdit = uieditfield(fig, 'Position', [150, 310, 200, 25]);

% Create a parent folder path label and edit field
parentFolderLabel = uilabel(fig, 'Text', 'Parent Folder Path:', 'Position', [50, 270, 120, 20]);
parentFolderEdit = uieditfield(fig, 'Position', [180, 270, 170, 25]);

% Create a reference directory label and edit field
refDirLabel = uilabel(fig, 'Text', 'Reference Directory:', 'Position', [50, 230, 120, 20]);
refDirEdit = uieditfield(fig, 'Position', [180, 230, 170, 25]);

% Create an output directory label and edit field
opDirLabel = uilabel(fig, 'Text', 'Output Directory:', 'Position', [50, 190, 120, 20]);
opDirEdit = uieditfield(fig, 'Position', [180, 190, 170, 25]);

% Create labels and edit fields in a column on the right side
columnX = 400;
spacingY = 40;

gaussianSigmaLabel = uilabel(fig, 'Text', 'Gaussian Sigma:', 'Position', [columnX, 300, 100, 20]);
gaussianSigmaEdit = uieditfield(fig, 'numeric', 'Position', [columnX + 100, 300, 100, 20], 'Value', 0.18);

binSizeLabel = uilabel(fig, 'Text', 'Bin Size:', 'Position', [columnX, 300 - spacingY, 100, 20]);
binSizeEdit = uieditfield(fig, 'numeric', 'Position', [columnX + 100, 300 - spacingY, 100, 20], 'Value', 0.3);

minPeakDistanceLabel = uilabel(fig, 'Text', 'Min Peak Distance:', 'Position', [columnX, 300 - 2 * spacingY, 120, 20]);
minPeakDistanceEdit = uieditfield(fig, 'numeric', 'Position', [columnX + 120, 300 - 2 * spacingY, 100, 20], 'Value', 0.025);

thresholdBurstLabel = uilabel(fig, 'Text', 'Threshold Burst:', 'Position', [columnX, 300 - 3 * spacingY, 120, 20]);
thresholdBurstEdit = uieditfield(fig, 'numeric', 'Position', [columnX + 120, 300 - 3 * spacingY, 100, 20], 'Value', 1.2);

thresholdStartStopLabel = uilabel(fig, 'Text', 'Threshold Start/Stop:', 'Position', [columnX, 300 - 4 * spacingY, 120, 20]);
thresholdStartStopEdit = uieditfield(fig, 'numeric', 'Position', [columnX + 120, 300 - 4 * spacingY, 100, 20], 'Value', 0.3);

% Create buttons at the bottom
buttonY = 50;

processActivity = uibutton(fig, 'Text', 'Process Activity', 'Position', [45, buttonY, 90, 25], 'ButtonPushedFcn', @processActivityButtonCallback);
processNetwork = uibutton(fig, 'Text', 'Process Network', 'Position', [155, buttonY, 90, 25], 'ButtonPushedFcn', @processNetworkButtonCallback);
clearButton = uibutton(fig, 'Text', 'Clear', 'Position', [265, buttonY, 80, 25], 'ButtonPushedFcn', @clearButtonCallback);
quitButton = uibutton(fig, 'Text', 'Quit', 'Position', [365, buttonY, 80, 25], 'ButtonPushedFcn', @quitButtonCallback);



    
    % Callback for processing data button
    function processActivityButtonCallback(~, ~)

    % Get the values from the edit fields
    projectName = projectNameEdit.Value;
    div0Date = div0DateEdit.Value;
    parentFolderPath = parentFolderEdit.Value;
    refDir = refDirEdit.Value;
    % Get the values from the sliders
    gaussianSigma = gaussianSigmaEdit.Value;
    binSize = binSizeEdit.Value;
    minPeakDistance = minPeakDistanceEdit.Value;
    thresholdBurst = thresholdBurstEdit.Value;
    thresholdStartStop = thresholdStartStopEdit.Value;
    opDir = opDirEdit.Value;


    % Create a struct and store the values
    data = struct();
    data.projectName = projectName;
    data.div0Date = div0Date;
    data.parentFolderPath = parentFolderPath;
    data.refDir = refDir;
    data.gaussianSigma = gaussianSigma;
    data.binSize = binSize;
    data.minPeakDistance = minPeakDistance;
    data.thresholdBurst = thresholdBurst;
    data.thresholdStartStop = thresholdStartStop;
    data.opDir = opDir;

    data.fig = fig;
    % Log the data into a log file
    logFileName = './network_log_file.txt';
    logFile = fopen(logFileName, 'a'); % 'a' for append mode
        
    if logFile == -1
        error('Error opening the log file.');
    end
    % Get the current timestamp
    currentTimestamp = datetime('now');
    fprintf(logFile, 'Timestamp: %s\n', char(currentTimestamp));
    data.logFile = logFile;
    fprintf(logFile,"Processing activity");
    disp("Activity Processing")
    fprintf(logFile, 'Project Name: %s\n', data.projectName);
    fprintf(logFile, 'DIV 0 Date: %s\n', data.div0Date);
    fprintf(logFile, 'Parent Folder Path: %s\n', data.parentFolderPath);
    fprintf(logFile, 'Reference Directory: %s\n', data.refDir);
    fprintf(logFile, 'Gaussian Sigma: %f\n', data.gaussianSigma);
    fprintf(logFile, 'Bin Size: %f\n', data.binSize);
    fprintf(logFile, 'Min Peak Distance: %f\n', data.minPeakDistance);
    fprintf(logFile, 'Threshold Burst: %f\n', data.thresholdBurst);
    fprintf(logFile, 'Threshold Start-Stop: %f\n', data.thresholdStartStop);
    fprintf(logFile, 'Output Directory: %s\n', data.opDir);


  

    % Now you have all the values stored in the "data" struct.
    % You can use it as required in your program logic.

    disp('Process Data button clicked.');
    disp('Activity File');
    compileActivityFiles(data);

    % Close the log file
    fclose(logFile);
    end
        % Callback for processing data button
    function processNetworkButtonCallback(~, ~)

    % Get the values from the edit fields
    projectName = projectNameEdit.Value;
    div0Date = div0DateEdit.Value;
    parentFolderPath = parentFolderEdit.Value;
    refDir = refDirEdit.Value;
      % Get the values from the sliders
    gaussianSigma = gaussianSigmaEdit.Value;
    binSize = binSizeEdit.Value;
    minPeakDistance = minPeakDistanceEdit.Value;
    thresholdBurst = thresholdBurstEdit.Value;
    thresholdStartStop = thresholdStartStopEdit.Value;
    opDir = opDirEdit.Value;


    % Create a struct and store the values
    data = struct();
    data.projectName = projectName;
    data.div0Date = div0Date;
    data.parentFolderPath = parentFolderPath;
    data.refDir = refDir;
    data.gaussianSigma = gaussianSigma;
    data.binSize = binSize;
    data.minPeakDistance = minPeakDistance;
    data.thresholdBurst = thresholdBurst;
    data.thresholdStartStop = thresholdStartStop;
    data.opDir = opDir;


    % Log the data into a log file
    logFileName = './activity_log_file.txt';
    logFile = fopen(logFileName, 'a'); % 'a' for append mode

    if logFile == -1
        error('Error opening the log file.');
    end
    data.logFile = logFile;
    data.fig = fig;
    currentTimestamp = datetime('now');
    fprintf(logFile, 'Timestamp: %s\n', char(currentTimestamp));
    fprintf(logFile,"Processing Network");
    disp("Network Processing")
    fprintf(logFile, 'Project Name: %s\n', data.projectName);
    fprintf(logFile, 'DIV 0 Date: %s\n', data.div0Date);
    fprintf(logFile, 'Parent Folder Path: %s\n', data.parentFolderPath);
    fprintf(logFile, 'Reference Directory: %s\n', data.refDir);
    fprintf(logFile, 'Gaussian Sigma: %f\n', data.gaussianSigma);
    fprintf(logFile, 'Bin Size: %f\n', data.binSize);
    fprintf(logFile, 'Min Peak Distance: %f\n', data.minPeakDistance);
    fprintf(logFile, 'Threshold Burst: %f\n', data.thresholdBurst);
    fprintf(logFile, 'Threshold Start-Stop: %f\n', data.thresholdStartStop);
    fprintf(logFile, 'Output Directory: %s\n', data.opDir);


  

    % Now you have all the values stored in the "data" struct.
    % You can use it as required in your program logic.

    disp('Process Network Data button clicked.');
    disp('Network File');
    compileNetworkFiles(data);

    % Close the log file
    fclose(logFile);
    end
    % Callback for clearing input fields button
    function clearButtonCallback(~, ~)
        % Clear the input fields
        projectNameEdit.Value = '';
        div0DateEdit.Value = '';
        parentFolderEdit.Value = '';
        refDirEdit.Value = '';
      % Reset the numeric fields to default values
        gaussianSigmaEdit.Value = 0.18;
        binSizeEdit.Value = 0.3;
        minPeakDistanceEdit.Value = 0.025;
        thresholdBurstEdit.Value = 1.2;
        thresholdStartStopEdit.Value = 0.3;
        opDirEdit.Value = '';

    end

    % Callback for quitting the GUI button
    function quitButtonCallback(~, ~)
        % Close the figure window
        delete(fig);
    end

    % Callback for Gaussian Sigma slider
    function gaussianSigmaCallback(slider, ~)
        value = slider.Value;
        gaussianSigmaValueLabel.Text = num2str(value);
    end

    % Callback for Bin Size slider
    function binSizeCallback(slider, ~)
        value = slider.Value;
        binSizeValueLabel.Text = num2str(value);
    end

    % Callback for Min Peak Distance slider
    function minPeakDistanceCallback(slider, ~)
        value = slider.Value;
        minPeakDistanceValueLabel.Text = num2str(value);
    end

    % Callback for Threshold Burst slider
    function thresholdBurstCallback(slider, ~)
        value = slider.Value;
        thresholdBurstValueLabel.Text = num2str(value);
    end

    % Callback for Threshold Start-Stop slider
    function thresholdStartStopCallback(slider, ~)
        value = slider.Value;
        thresholdStartStopValueLabel.Text = num2str(value);
    end

   

end
