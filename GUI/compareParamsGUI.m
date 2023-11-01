function mainFunction
    % Create a uifigure window
    addpath('../../MEA_Analysis/');
    addpath('../../MxW_Matlab_22.2/');
    fig = uifigure('Name', 'Plot Parameters', 'Position', [100, 100, 700, 400]);

    % Create UI components
  % Create labels and edit fields instead of sliders
    dataDirLabel = uilabel(fig, 'Text', 'Specific DIV Path:', 'Position', [50, 300, 120, 20]);
    dataDirEdit = uieditfield(fig, 'Position', [180, 300, 170, 25]);

    refDirLabel = uilabel(fig, 'Text', 'Reference Directory:', 'Position', [50, 250, 120, 20]);
    refDirEdit = uieditfield(fig, 'Position', [180, 250, 170, 25]);

    opDirLabel = uilabel(fig, 'Text', 'Output Directory:', 'Position', [50, 200, 120, 20]);
    opDirEdit = uieditfield(fig, 'Position', [180, 200, 170, 25]);

    gaussianSigmaLabel = uilabel(fig, 'Text', 'Gaussian Sigma:', 'Position', [400, 300, 100, 20]);
    gaussianSigmaEdit = uieditfield(fig, 'numeric', 'Position', [500, 300, 100, 20], 'Value', 0.18);

    binSizeLabel = uilabel(fig, 'Text', 'Bin Size:', 'Position', [400, 250, 100, 20]);
    binSizeEdit = uieditfield(fig, 'numeric', 'Position', [500, 250, 100, 20], 'Value', 0.3);

    minPeakDistanceLabel = uilabel(fig, 'Text', 'Min Peak Distance:', 'Position', [400, 200, 120, 20]);
    minPeakDistanceEdit = uieditfield(fig, 'numeric', 'Position', [500, 200, 100, 20], 'Value', 0.025);

    thresholdBurstLabel = uilabel(fig, 'Text', 'Threshold Burst:', 'Position', [400, 150, 120, 20]);
    thresholdBurstEdit = uieditfield(fig, 'numeric', 'Position', [500, 150, 100, 20], 'Value', 1.2);

    thresholdStartStopLabel = uilabel(fig, 'Text', 'Threshold Start/Stop:', 'Position', [400, 100, 120, 20]);
    thresholdStartStopEdit = uieditfield(fig, 'numeric', 'Position', [500, 100, 100, 20], 'Value', 0.3);

    processButton = uibutton(fig, 'Text', 'Process Data', 'Position', [50, 150, 100, 30],'ButtonPushedFcn',@processButtonCallback);     % Add your callback function
    clearButton = uibutton(fig, 'Text', 'Clear', 'Position', [175, 150, 100, 30],'ButtonPushedFcn',@clearButtonCallback);  % Add your callback function
    quitButton = uibutton(fig, 'Text', 'Quit', 'Position', [300, 150, 100, 30],'ButtonPushedFcn',@quitButtonCallback);  % Add your callback function

    % Event handler for the process button
    function processButtonCallback(~, ~)
        % Get the values from the edit fields
        dataDirPath = dataDirEdit.Value;
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
        logFileName = './parm_log_file.txt';
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

    % Event handler for the clear button
    function clearButtonCallback(~, ~)
        % Clear the input fields
        dataDirEdit.Value = '';
        refDirEdit.Value = '';
        opDirEdit.Value = '';
        % Reset the numeric fields to default values
        gaussianSigmaEdit.Value = 0.18;
        binSizeEdit.Value = 0.3;
        minPeakDistanceEdit.Value = 0.025;
        thresholdBurstEdit.Value = 1.2;
        thresholdStartStopEdit.Value = 0.3;

    end

    % Event handler for the quit button
    function quitButtonCallback(~, ~)
        % Close the figure window
        close(fig);
    end

    end
