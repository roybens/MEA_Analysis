function mainFunction
    % Create a uifigure window with specified position and callback for dynamic layout
    fig = uifigure('Name', 'Process Activity /Network', 'Position', [100, 100, 700, 500], 'SizeChangedFcn', @updateLayout);
    fig.Color = 'white'; % Change the background color to white for better contrast

    % Improved alignment and consistent spacing
    labelWidth = 110;
    editWidth = 200;
    editX = 170;
    labelX = editX - labelWidth - 10;
    startY = 420;
    spacingY = 30;

    % Create a project name label and edit field
    projectNameLabel = uilabel(fig, 'Text', 'Project Name:', 'Position', [labelX, startY, labelWidth, 22]);
    projectNameEdit = uieditfield(fig, 'Position', [editX, startY, editWidth, 22]);

    % Same spacing and positioning scheme applied to all elements for consistency
    div0DateLabel = uilabel(fig, 'Text', 'DIV 0 Date:', 'Position', [labelX, startY - spacingY, labelWidth, 22]);
    div0DateEdit = uieditfield(fig, 'Position', [editX, startY - spacingY, editWidth, 22]);

    parentFolderLabel = uilabel(fig, 'Text', 'Data Path:', 'Position', [labelX, startY - 2 * spacingY, labelWidth, 22]);
    parentFolderEdit = uieditfield(fig, 'Position', [editX, startY - 2 * spacingY, editWidth, 22]);

    refDirLabel = uilabel(fig, 'Text', 'Reference File:', 'Position', [labelX, startY - 3 * spacingY, labelWidth, 22]);
    refDirEdit = uieditfield(fig, 'Position', [editX, startY - 3 * spacingY, editWidth, 22]);

    opDirLabel = uilabel(fig, 'Text', 'Output Directory:', 'Position', [labelX, startY - 4 * spacingY, labelWidth, 22]);
    opDirEdit = uieditfield(fig, 'Position', [editX, startY - 4 * spacingY, editWidth, 22]);

    % Right column for numerical inputs, aligned with labels on the left for visual balance
    columnX = fig.Position(3) - editWidth - labelWidth - 20; % Automatically aligns to the right
    gaussianSigmaLabel = uilabel(fig, 'Text', 'Gaussian Sigma:', 'Position', [columnX, startY, labelWidth, 22]);
    gaussianSigmaEdit = uieditfield(fig, 'numeric', 'Position', [columnX + labelWidth + 10, startY, editWidth, 22], 'Value', 0.14);

    binSizeLabel = uilabel(fig, 'Text', 'Bin Size:', 'Position', [columnX, startY - spacingY, labelWidth, 22]);
    binSizeEdit = uieditfield(fig, 'numeric', 'Position', [columnX + labelWidth + 10, startY - spacingY, editWidth, 22], 'Value', 0.1);

    minPeakDistanceLabel = uilabel(fig, 'Text', 'Min Peak Distance:', 'Position', [columnX, startY - 2 * spacingY, labelWidth, 22]);
    minPeakDistanceEdit = uieditfield(fig, 'numeric', 'Position', [columnX + labelWidth + 10, startY - 2 * spacingY, editWidth, 22], 'Value', 1.0);

    thresholdBurstLabel = uilabel(fig, 'Text', 'Threshold Burst:', 'Position', [columnX, startY - 3 * spacingY, labelWidth, 22]);
    thresholdBurstEdit = uieditfield(fig, 'numeric', 'Position', [columnX + labelWidth + 10, startY - 3 * spacingY, editWidth, 22], 'Value', 1.0);

    thresholdStartStopLabel = uilabel(fig, 'Text', 'Threshold Start/Stop:', 'Position', [columnX, startY - 4 * spacingY, labelWidth, 22]);
    thresholdStartStopEdit = uieditfield(fig, 'numeric', 'Position', [columnX + labelWidth + 10, startY - 4 * spacingY, editWidth, 22], 'Value', 0.3);

    % Centered buttons with consistent sizing
    buttonY = 50;
    buttonWidth = 125;
    buttonHeight = 30;
    spacingX = 10;
    buttonStartX = (fig.Position(3) - 5 * buttonWidth - 4 * spacingX) / 2; % Centering buttons

    % Updated buttons with consistent sizes and callback functions defined below
    processActivityButton = uibutton(fig, 'Text', 'Process Activity', 'Position', [buttonStartX, buttonY, buttonWidth, buttonHeight], 'ButtonPushedFcn', @processActivityButtonCallback);
    processNetworkButton = uibutton(fig, 'Text', 'Process Network', 'Position', [buttonStartX + buttonWidth + spacingX, buttonY, buttonWidth, buttonHeight], 'ButtonPushedFcn', @processNetworkButtonCallback);
    exploreParamsButton = uibutton(fig, 'Text', 'Explore Params', 'Position', [buttonStartX + 2 * buttonWidth + 2 * spacingX, buttonY, buttonWidth, buttonHeight], 'ButtonPushedFcn', @exploreParamsButtonCallback);
    clearButton = uibutton(fig, 'Text', 'Clear', 'Position', [buttonStartX + 3 * buttonWidth + 3 * spacingX, buttonY, buttonWidth, buttonHeight], 'ButtonPushedFcn', @clearButtonCallback);
    quitButton = uibutton(fig, 'Text', 'Quit', 'Position', [buttonStartX + 4 * buttonWidth + 4 * spacingX, buttonY, buttonWidth, buttonHeight], 'ButtonPushedFcn', @quitButtonCallback);
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

    % Event handler for the process button
    function exploreParamsButtonCallback(~, ~)
        % Get the values from the edit fields
        dataDirPath = parentFolderEdit.Value;
        refDir = refDirEdit.Value;
        opDir = opDirEdit.Value;

        % Get the values from the sliders
        gaussianSigma = gaussianSigmaEdit.Value;
        binSize = binSizeEdit.Value;
        minPeakDistance = minPeakDistanceEdit.Value;
        thresholdBurst = thresholdBurstEdit.Value;
        thresholdStartStop = thresholdStartStopEdit.Value;

        % Validate the fields
        if isempty(dataDirPath) || isempty(refDir) || isempty(opDir)
            disp('Error: Please fill in all the required fields.');
            return;
        end
         % Specify the log file name and open it in append mode
        logFileName = './explore_parm_log_file.txt';
        logFile = fopen(logFileName, 'a'); % 'a' for append mode
    
        % Check if the file was opened successfully
        if logFile == -1
            error('Error opening the log file.');
        end
        % Create the data structure
        data.dataDir = dataDirPath;
        data.refDir = refDir;
        data.gaussianSigma = gaussianSigma;
        data.binSize = binSize;
        data.minPeakDistance = minPeakDistance;
        data.thresholdBurst = thresholdBurst;
        data.thresholdStartStop = thresholdStartStop;
        data.opDir = opDir;
        data.fig = fig;
        data.logFile = logFile;
        fprintf(logFile, 'data.dataDir = %s\n', data.dataDir);
        fprintf(logFile, 'data.refDir = %s\n', data.refDir);
        fprintf(logFile, 'data.gaussianSigma = %f\n', data.gaussianSigma);
        fprintf(logFile, 'data.binSize = %f\n', data.binSize);
        fprintf(logFile, 'data.minPeakDistance = %f\n', data.minPeakDistance);
        fprintf(logFile, 'data.thresholdBurst = %f\n', data.thresholdBurst);
        fprintf(logFile, 'data.thresholdStartStop = %f\n', data.thresholdStartStop);
        fprintf(logFile, 'data.opDir = %s\n', data.opDir);
        

        % Call the processing function
        plotParametersComparesWrapper(data);

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
