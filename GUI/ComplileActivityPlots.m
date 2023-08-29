function mainFunction
    % Create a uifigure window
    fig = uifigure('Name', 'GUI Example', 'Position', [100, 100, 700, 500], 'SizeChangedFcn', @updateLayout);

    % Create a project name label
    projectNameLabel = uilabel(fig, 'Text', 'Project Name:', 'Position', [50, 350, 100, 20]);

    % Create an edit field for project name
    projectNameEdit = uieditfield(fig, 'Position', [150, 350, 200, 25]);

    % Create a DIV 0 date label
    div0DateLabel = uilabel(fig, 'Text', 'DIV 0 Date:', 'Position', [50, 310, 100, 20]);

    % Create an edit field for DIV 0 date
    div0DateEdit = uieditfield(fig, 'Position', [150, 310, 200, 25]);

    % Create a parent folder path label
    parentFolderLabel = uilabel(fig, 'Text', 'Parent Folder Path:', 'Position', [50, 270, 120, 20]);

    % Create an edit field for parent folder path
    parentFolderEdit = uieditfield(fig, 'Position', [180, 270, 170, 25]);

    % Create a reference directory label
    refDirLabel = uilabel(fig, 'Text', 'Reference Directory:', 'Position', [50, 230, 120, 20]);

    % Create an edit field for reference directory
    refDirEdit = uieditfield(fig, 'Position', [180, 230, 170, 25]);

    % Create an output directory label
    opDirLabel = uilabel(fig, 'Text', 'Output Directory:', 'Position', [50, 190, 120, 20]);

    % Create an edit field for output directory
    opDirEdit = uieditfield(fig, 'Position', [180, 190, 170, 25]);

    % Create a button for processing data
    processActivity = uibutton(fig, 'Text', 'Process Activity', 'Position', [45, 150, 90, 25], 'ButtonPushedFcn', @processActivityButtonCallback);
    % Create a button for processing data
    processNetwork = uibutton(fig, 'Text', 'Process Network', 'Position', [155, 150, 90, 25], 'ButtonPushedFcn', @processNetworkButtonCallback);
    % Create a button for clearing input fields
    clearButton = uibutton(fig, 'Text', 'Clear', 'Position', [265, 150, 80, 25], 'ButtonPushedFcn', @clearButtonCallback);

    % Create a button for quitting the GUI
    quitButton = uibutton(fig, 'Text', 'Quit', 'Position', [365, 150,80,25], 'ButtonPushedFcn', @quitButtonCallback);

    % Create a slider for Gaussian kernel standard deviation
    gaussianSigmaLabel = uilabel(fig, 'Text', 'Gaussian Sigma:', 'Position', [350, 350, 100, 20]);
    gaussianSigmaSlider = uislider(fig, 'Position', [450, 350, 100, 20], 'Limits', [0, 1], 'Value', 0.18, 'ValueChangedFcn', @gaussianSigmaCallback);
    gaussianSigmaValueLabel = uilabel(fig, 'Text', '0.18', 'Position', [560, 350, 30, 20]);

    % Create a slider for histogram bin size
    binSizeLabel = uilabel(fig, 'Text', 'Bin Size:', 'Position', [350, 310, 100, 20]);
    binSizeSlider = uislider(fig, 'Position', [450, 310, 100, 20], 'Limits', [0, 1], 'Value', 0.3, 'ValueChangedFcn', @binSizeCallback);
    binSizeValueLabel = uilabel(fig, 'Text', '0.3', 'Position', [560, 310, 30, 20]);

    % Create a slider for minimum peak distance
    minPeakDistanceLabel = uilabel(fig, 'Text', 'Min Peak Distance:', 'Position', [350, 270, 120, 20]);
    minPeakDistanceSlider = uislider(fig, 'Position', [450, 270, 100, 20], 'Limits', [0, 1], 'Value', 0.025, 'ValueChangedFcn', @minPeakDistanceCallback);
    minPeakDistanceValueLabel = uilabel(fig, 'Text', '0.025', 'Position', [560, 270, 30, 20]);

    % Create a slider for burst detection threshold
    thresholdBurstLabel = uilabel(fig, 'Text', 'Threshold Burst:', 'Position', [350, 230, 120, 20]);
    thresholdBurstSlider = uislider(fig, 'Position', [450, 230, 100, 20], 'Limits', [0, 2], 'Value', 1.2, 'ValueChangedFcn', @thresholdBurstCallback);
    thresholdBurstValueLabel = uilabel(fig, 'Text', '1.2', 'Position', [560, 230, 30, 20]);

    % Create a slider for threshold start-stop
    thresholdStartStopLabel = uilabel(fig, 'Text', 'Threshold Start-Stop:', 'Position', [350, 190, 120, 20]);
    thresholdStartStopSlider = uislider(fig, 'Position', [450, 190, 100, 20], 'Limits', [0, 1], 'Value', 0.3, 'ValueChangedFcn', @thresholdStartStopCallback);
    thresholdStartStopValueLabel = uilabel(fig, 'Text', '0.3', 'Position', [560, 190, 30, 20]);

    % Create checkboxes for assay types
    assayTypesLabel = uilabel(fig, 'Text', 'Net Assays:', 'Position', [50, 130, 100, 20]);
    todayCheckbox = uicheckbox(fig, 'Text', 'Today', 'Position', [150, 130, 60, 20]);
    lastCheckbox = uicheckbox(fig, 'Text', 'Last', 'Position', [220, 130, 60, 20]);
    bestCheckbox = uicheckbox(fig, 'Text', 'Best', 'Position', [290, 130, 60, 20]);

    % Create checkboxes for assay types
    activityAssayTypesLabel = uilabel(fig, 'Text', 'Activity Assays:', 'Position', [350, 130, 100, 20]);
    sparse7xCheckbox = uicheckbox(fig, 'Text', 'Sparse 7x', 'Position', [450, 130, 60, 20]);
    
    % Create exclude_chips label and edit field
    excludeChipsLabel = uilabel(fig, 'Text', 'Exclude Chips:', 'Position', [50, 90, 120, 20]);
    excludeChipsEdit = uieditfield(fig, 'Position', [180, 90, 170, 25]);

    % Create exclude_runids label and edit field
    excludeRunIDsLabel = uilabel(fig, 'Text', 'Exclude Run IDs:', 'Position', [50, 50, 120, 20]);
    excludeRunIDsEdit = uieditfield(fig, 'Position', [180, 50, 170, 25]);

    % Create a button for plotting graph
    plotActButton = uibutton(fig, 'Text', 'Act settings', 'Position', [360, 50, 100, 30], 'ButtonPushedFcn', @plotActButtonCallback);
    % Create a button for plotting graph
    plotNetButton = uibutton(fig, 'Text', 'Net settings', 'Position', [460, 50, 100, 30], 'ButtonPushedFcn', @plotNetButtonCallback);
    
%     % Callback for updating the layout when figure is resized
%     function updateLayout(~, ~)
%         % Get the current figure size
%         figSize = fig.Position(3:4);
%         
%         % Calculate the new positions of the UI components based on the figure size
%         
%         % Project Name components
%         projectNameLabel.Position(2) = figSize(2) - 50;
%         projectNameEdit.Position(2) = figSize(2) - 50;
%         
%         % DIV 0 Date components
%         div0DateLabel.Position(2) = figSize(2) - 90;
%         div0DateEdit.Position(2) = figSize(2) - 90;
%         
%         % Parent Folder Path components
%         parentFolderLabel.Position(2) = figSize(2) - 130;
%         parentFolderEdit.Position(2) = figSize(2) - 130;
%         
%         % Reference Directory components
%         refDirLabel.Position(2) = figSize(2) - 170;
%         refDirEdit.Position(2) = figSize(2) - 170;
%         
%         % Output Directory components
%         opDirLabel.Position(2) = figSize(2) - 210;
%         opDirEdit.Position(2) = figSize(2) - 210;
%         
%         % Process Data button
%         processButton.Position(2) = figSize(2) - 250;
%         
%         % Clear button
%         clearButton.Position(2) = figSize(2) - 250;
%         clearButton.Position(1) = 175 + (figSize(1) - 600);
%         
%         % Quit button
%         quitButton.Position(2) = figSize(2) - 250;
%         quitButton.Position(1) = 300 + (figSize(1) - 600);
%         
%         % Gaussian Sigma components
%         gaussianSigmaLabel.Position(2) = figSize(2) - 50;
%         gaussianSigmaSlider.Position(2) = figSize(2) - 50;
%         gaussianSigmaValueLabel.Position(2) = figSize(2) - 50;
%         
%         % Bin Size components
%         binSizeLabel.Position(2) = figSize(2) - 90;
%         binSizeSlider.Position(2) = figSize(2) - 90;
%         binSizeValueLabel.Position(2) = figSize(2) - 90;
%         
%         % Min Peak Distance components
%         minPeakDistanceLabel.Position(2) = figSize(2) - 130;
%         minPeakDistanceSlider.Position(2) = figSize(2) - 130;
%         minPeakDistanceValueLabel.Position(2) = figSize(2) - 130;
%         
%         % Threshold Burst components
%         thresholdBurstLabel.Position(2) = figSize(2) - 170;
%         thresholdBurstSlider.Position(2) = figSize(2) - 170;
%         thresholdBurstValueLabel.Position(2) = figSize(2) - 170;
%         
%         % Threshold Start-Stop components
%         thresholdStartStopLabel.Position(2) = figSize(2) - 210;
%         thresholdStartStopSlider.Position(2) = figSize(2) - 210;
%         thresholdStartStopValueLabel.Position(2) = figSize(2) - 210;
%         
%         % Assay Types components
%         assayTypesLabel.Position(2) = figSize(2) - 250;
%         todayCheckbox.Position(2) = figSize(2) - 250;
%         lastCheckbox.Position(2) = figSize(2) - 250;
%         bestCheckbox.Position(2) = figSize(2) - 250;
%         
%         % Exclude Chips components
%         excludeChipsLabel.Position(2) = figSize(2) - 290;
%         excludeChipsEdit.Position(2) = figSize(2) - 290;
%         
%         % Exclude Run IDs components
%         excludeRunIDsLabel.Position(2) = figSize(2) - 330;
%         excludeRunIDsEdit.Position(2) = figSize(2) - 330;
%         
%         % Plot Graph button
%         plotButton.Position(2) = figSize(2) - 370;
%     end
    
    % Callback for processing data button
    function processActivityButtonCallback(~, ~)

    % Get the values from the edit fields
    projectName = projectNameEdit.Value;
    div0Date = div0DateEdit.Value;
    parentFolderPath = parentFolderEdit.Value;
    refDir = refDirEdit.Value;
    gaussianSigma = gaussianSigmaSlider.Value;
    binSize = binSizeSlider.Value;
    minPeakDistance = minPeakDistanceSlider.Value;
    thresholdBurst = thresholdBurstSlider.Value;
    thresholdStartStop = thresholdStartStopSlider.Value;
    opDir = opDirEdit.Value;
    assayTypes = [todayCheckbox.Value, lastCheckbox.Value, bestCheckbox.Value];
    excludeChips = excludeChipsEdit.Value;
    excludeRunIDs = excludeRunIDsEdit.Value;

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
    data.assayTypes = assayTypes;
    data.excludeChips = excludeChips;
    data.excludeRunIDs = excludeRunIDs;
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
    fprintf(logFile, 'Assay Types: %d %d %d\n', data.assayTypes);
    fprintf(logFile, 'Exclude Chips: %s\n', data.excludeChips);
    fprintf(logFile, 'Exclude Run IDs: %s\n', data.excludeRunIDs);

  

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
    gaussianSigma = gaussianSigmaSlider.Value;
    binSize = binSizeSlider.Value;
    minPeakDistance = minPeakDistanceSlider.Value;
    thresholdBurst = thresholdBurstSlider.Value;
    thresholdStartStop = thresholdStartStopSlider.Value;
    opDir = opDirEdit.Value;
    assayTypes = [todayCheckbox.Value, lastCheckbox.Value, bestCheckbox.Value];
    excludeChips = excludeChipsEdit.Value;
    excludeRunIDs = excludeRunIDsEdit.Value;

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
    data.assayTypes = assayTypes;
    data.excludeChips = excludeChips;
    data.excludeRunIDs = excludeRunIDs;

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
    fprintf(logFile, 'Assay Types: %d %d %d\n', data.assayTypes);
    fprintf(logFile, 'Exclude Chips: %s\n', data.excludeChips);
    fprintf(logFile, 'Exclude Run IDs: %s\n', data.excludeRunIDs);

  

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
        gaussianSigmaSlider.Value = 0.18;
        binSizeSlider.Value = 0.3;
        minPeakDistanceSlider.Value = 0.025;
        thresholdBurstSlider.Value = 1.2;
        thresholdStartStopSlider.Value = 0.3;
        opDirEdit.Value = '';
        todayCheckbox.Value = false;
        lastCheckbox.Value = false;
        bestCheckbox.Value = false;
        sparse7xCheckbox.Value = false;
        excludeChipsEdit.Value = '';
        excludeRunIDsEdit.Value = '';
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

    % Callback for plot button
    function plotActButtonCallback(~, ~)
        % Implement your logic here
        disp('Plot Act Graph button clicked.');
        projectName = projectNameEdit.Value;
        refDir = refDirEdit.Value;
        opDir = opDirEdit.Value;
        assayTypes = {'sparse 7x'}; %need to remove the hardcoding
        excludeChips = excludeChipsEdit.Value;
        excludeRunIDs = excludeRunIDsEdit.Value;
        % Create a struct to hold the settings
        settings = struct();
        settings.projectName = projectName;
        settings.refDir = refDir;
        settings.opDir = opDir;
        settings.assayTypes = assayTypes;
        settings.excludeChips = excludeChips;
        settings.excludeRunIDs = excludeRunIDs;
    
        % Convert the settings struct to a JSON string
        jsonString = jsonencode(settings);
    
        % Save the JSON string to a file
        settingsFileName = 'act_plt_settings.json';
        fileID = fopen(settingsFileName, 'w');
        fprintf(fileID, jsonString);
        fclose(fileID);
        msgbox('settings saved')


    end
    % Callback for plot button
    function plotNetButtonCallback(~, ~)
        % Implement your logic here
        disp('Plot Network Graph button clicked.');
        projectName = projectNameEdit.Value;
        refDir = refDirEdit.Value;
        opDir = opDirEdit.Value;
        assayTypes = {};
        % Check the checkboxes and add assay types to the cell array if checked
        if todayCheckbox.Value
            assayTypes{end+1} = 'today';
        end
        
        if lastCheckbox.Value
            assayTypes{end+1} = 'last';
        end
        
        if bestCheckbox.Value
            assayTypes{end+1} = 'best';
        end
        %assayTypes = [todayCheckbox.Value, lastCheckbox.Value, bestCheckbox.Value];
        excludeChips = excludeChipsEdit.Value;
        excludeRunIDs = excludeRunIDsEdit.Value;
        % Create a struct to hold the settings
        settings = struct();
        settings.projectName = projectName;
        settings.refDir = refDir;
        settings.opDir = opDir;
        settings.assayTypes = assayTypes;
        settings.excludeChips = excludeChips;
        settings.excludeRunIDs = excludeRunIDs;
    
        % Convert the settings struct to a JSON string
        jsonString = jsonencode(settings);
    
        % Save the JSON string to a file
        settingsFileName = 'net_plt_settings.json';
        fileID = fopen(settingsFileName, 'w');
        fprintf(fileID, jsonString);
        fclose(fileID);
        msgbox('Settings saved')

    end

end
