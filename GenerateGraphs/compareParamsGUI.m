function mainFunction
    % Create a uifigure window
    fig = uifigure('Name', 'Plot Parameters', 'Position', [100, 100, 700, 400]);

    % Create UI components
    dataDirLabel = uilabel(fig, 'Text', 'Specific DIV Path:', 'Position', [50, 300, 120, 20]);
    dataDirEdit = uieditfield(fig, 'Position', [180, 300, 170, 25]);

    refDirLabel = uilabel(fig, 'Text', 'Reference Directory:', 'Position', [50, 250, 120, 20]);
    refDirEdit = uieditfield(fig, 'Position', [180, 250, 170, 25]);

    opDirLabel = uilabel(fig, 'Text', 'Output Directory:', 'Position', [50, 200, 120, 20]);
    opDirEdit = uieditfield(fig, 'Position', [180, 200, 170, 25]);

    gaussianSigmaLabel = uilabel(fig, 'Text', 'Gaussian Sigma:', 'Position', [400, 300, 100, 20]);
    gaussianSigmaSlider = uislider(fig, 'Position', [500, 300, 100, 20], 'Limits', [0, 1], 'Value', 0.18, 'ValueChangedFcn', @gaussianSigmaCallback);
    gaussianSigmaValueLabel = uilabel(fig, 'Text', '0.18', 'Position', [610, 300, 30, 20]);

    binSizeLabel = uilabel(fig, 'Text', 'Bin Size:', 'Position', [400, 250, 100, 20]);
    binSizeSlider = uislider(fig, 'Position', [500, 250, 100, 20], 'Limits', [0, 1], 'Value', 0.3, 'ValueChangedFcn', @binSizeCallback);
    binSizeValueLabel = uilabel(fig, 'Text', '0.3', 'Position', [610, 250, 30, 20]);

    minPeakDistanceLabel = uilabel(fig, 'Text', 'Min Peak Distance:', 'Position', [400, 200, 120, 20]);
    minPeakDistanceSlider = uislider(fig, 'Position', [500, 200, 100, 20], 'Limits', [0, 1], 'Value', 0.025, 'ValueChangedFcn', @minPeakDistanceCallback);
    minPeakDistanceValueLabel = uilabel(fig, 'Text', '0.025', 'Position', [610, 200, 30, 20]);

    thresholdBurstLabel = uilabel(fig, 'Text', 'Threshold Burst:', 'Position', [400, 150, 120, 20]);
    thresholdBurstSlider = uislider(fig, 'Position', [500, 150, 100, 20], 'Limits', [0, 2], 'Value', 1.2, 'ValueChangedFcn', @thresholdBurstCallback);
    thresholdBurstValueLabel = uilabel(fig, 'Text', '1.2', 'Position', [610, 150, 30, 20]);

    thresholdStartStopLabel = uilabel(fig, 'Text', 'Threshold Start/Stop:', 'Position', [400, 100, 120, 20]);
    thresholdStartStopSlider = uislider(fig, 'Position', [500, 100, 100, 20], 'Limits', [0, 1], 'Value', 0.3, 'ValueChangedFcn', @thresholdStartStopCallback);
    thresholdStartStopValueLabel = uilabel(fig, 'Text', '0.3', 'Position', [610, 100, 30, 20]);

    processButton = uibutton(fig, 'Text', 'Process Data', 'Position', [50, 150, 100, 30], 'ButtonPushedFcn', @processButtonCallback);
    clearButton = uibutton(fig, 'Text', 'Clear', 'Position', [175, 150, 100, 30], 'ButtonPushedFcn', @clearButtonCallback);
    quitButton = uibutton(fig, 'Text', 'Quit', 'Position', [300, 150, 100, 30], 'ButtonPushedFcn', @quitButtonCallback);

    % Event handler for the process button
    function processButtonCallback(~, ~)
        % Get the values from the edit fields
        dataDirPath = dataDirEdit.Value;
        refDir = refDirEdit.Value;
        opDir = opDirEdit.Value;

        % Get the values from the sliders
        gaussianSigma = gaussianSigmaSlider.Value;
        binSize = binSizeSlider.Value;
        minPeakDistance = minPeakDistanceSlider.Value;
        thresholdBurst = thresholdBurstSlider.Value;
        thresholdStartStop = thresholdStartStopSlider.Value;

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
        fprintf(logFile, 'data.fig = %d\n', data.fig);

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
    end

    % Event handler for the quit button
    function quitButtonCallback(~, ~)
        % Close the figure window
        close(fig);
    end

    % Callback for Gaussian Sigma slider
    function gaussianSigmaCallback(~, event)
        value = event.Value;
        gaussianSigmaValueLabel.Text = num2str(value);
    end

    % Callback for Bin Size slider
    function binSizeCallback(~, event)
        value = event.Value;
        binSizeValueLabel.Text = num2str(value);
    end

    % Callback for Min Peak Distance slider
    function minPeakDistanceCallback(~, event)
        value = event.Value;
        minPeakDistanceValueLabel.Text = num2str(value);
    end

    % Callback for Threshold Burst slider
    function thresholdBurstCallback(~, event)
        value = event.Value;
        thresholdBurstValueLabel.Text = num2str(value);
    end
    
    % Callback for Threshold Start/Stop slider
    function thresholdStartStopCallback(~, event)
        value = event.Value;
        thresholdStartStopValueLabel.Text = num2str(value);
    end
end
