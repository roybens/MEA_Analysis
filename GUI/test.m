%% Plot all lines on one plot
% Get a list of all files in the folder with the desired file name pattern.
IBI_max = 0;
BurstPeak_max = 0;
nBursts_max = 0;
spikePerBurst_max = 0;

filePattern = fullfile(parentFolderPath, '*.csv'); 
theFiles = dir(filePattern);
for f = 1 : length(theFiles)
    baseFileName = theFiles(f).name;
    pathFileNetwork = fullfile(theFiles(f).folder, baseFileName);
    datafile = (pathFileNetwork);
    data = readtable(datafile,'PreserveVariableNames',true);
    IBI_max = max([IBI_max, max(data.("IBI"))]);
    BurstPeak_max = max([BurstPeak_max, max(data.("Burst Peak"))]);
    nBursts_max = max([nBursts_max, max(data.("# Bursts"))]);
    spikePerBurst_max = max([spikePerBurst_max, max(data.("Spikes per Burst"))]);
end

fig2 =figure('color','w','Position',[0 0 800 800],'Visible','off');
% Define a map to hold genoStr-color pairs
genoColorMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
legendHandlesMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

% Define a list of colors to use for plots
colorList = {'k', 'r', 'b', 'c', 'm', 'y', 'g'}; % Add more colors as needed
colorIndex = 1;

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(parentFolderPath, '*.csv'); 
theFiles = dir(filePattern);
iterations = length(theFiles);
newGenoStrs = cell(iterations,1);
newIBIs = cell(iterations,1);
newBPs = cell(iterations,1);
newSPBs=cell(iterations,1);
newNBs=cell(iterations,1);
for k = 1 : iterations
    baseFileName = theFiles(k).name;
    %fullFileName = fullfile(theFiles(k).folder, baseFileName);
    pathFileNetwork = fullfile(theFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', pathFileNetwork);
    
    datafile = (pathFileNetwork);
    data = readtable(datafile,'PreserveVariableNames',true);
    
    extractChipID=data.("ChipID")(1);
    WellID=data.("WellID")(1);
    %extractChipID = regexp(pathFileNetwork,'\d{5}?','match');

    %idDouble = str2double(extractChipIDWellID);
%     if ismember(idDouble,wt)
%         geno = 'WT';
%         genoStr = string(geno);
%     end
%     
%     if ismember(idDouble,het)
%         geno = 'HET';
%         genoStr = string(geno);
%     end    
    genoStr=data.("NeuronType"){1};
    % Check if the genoStr is already in the map
    if ~isKey(genoColorMap, genoStr)
        % Assign the next color from colorList to this genoStr
        genoColorMap(genoStr) = colorList{colorIndex};
        
        % Increment the colorIndex, and reset if it exceeds the length of colorList
        colorIndex = mod(colorIndex, length(colorList)) + 1;
    end
    newGenoStrs{k} = genoStr;
    newIBIs{k} = data.("IBI");
    color = genoColorMap(genoStr);
    subplot(3,2,1);
    plot(data.(parameter),data.("IBI"),color);
    title('IBI')
    xlabel(plot_x_title)
    %xticks(parameter_start:plot_inc:parameter_end)
    xticks('auto')
    ylim([0 IBI_max*10/8])
    xline(param_val,'b--', 'LineWidth', 0.5)
    % Add a text annotation for the xline value
    text(param_val, 0, ['x = ', num2str(param_val)], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')
    grid on
    hold on
    newBPs{k} = data.("Burst Peak");
    subplot(3,2,2);
    plot(data.(parameter),data.("Burst Peak"),color);
    title('Burst Peak')
    xlabel(plot_x_title)
    %xticks(parameter_start:plot_inc:parameter_end)
    xticks('auto')
    ylim([0 BurstPeak_max*10/8])
    xline(param_val,'b--', 'LineWidth', 0.5)
    % Add a text annotation for the xline value
    text(param_val, 0, ['x = ', num2str(param_val)], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')
    
    grid on
    hold on
    newNBs{k} = data.("# Bursts");
    subplot(3,2,3);
    plot(data.(parameter),data.("# Bursts"),color);
    title('# of Bursts')
    xlabel(plot_x_title)
    %xticks(parameter_start:plot_inc:parameter_end)
    xticks('auto')
    ylim([0 nBursts_max*10/8])
    xline(param_val,'b--', 'LineWidth', 0.5)
    % Add a text annotation for the xline value
    text(param_val, 0, ['x = ', num2str(param_val)], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')
    
    grid on
    hold on
    newSPBs{k} = data.("Spikes per Burst");
    subplot(3,2,4);
    p=plot(data.(parameter),data.("Spikes per Burst"),color);
    title('Spikes per Burst')
    xlabel(plot_x_title)
    %xticks(parameter_start:plot_inc:parameter_end)
    xticks('auto')
    ylim([0 spikePerBurst_max*10/8])
    xline(param_val,'b--', 'LineWidth', 0.5)
    % Add a text annotation for the xline value
    text(param_val, 0, ['x = ', num2str(param_val)], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')
    
    grid on
    hold on
    % If this is the first time this genoStr is plotted, add its handle to legendHandles
    if  ~isKey(legendHandlesMap, genoStr)
        legendHandlesMap(genoStr) = p; % Store the plot handle for the legend
    end
end

legendHandles = values(legendHandlesMap, keys(genoColorMap));
legendLabels = keys(genoColorMap);
legend([legendHandles{:}], legendLabels);
newTable = table(newGenoStrs, newIBIs,newBPs,newNBs,newSPBs, 'VariableNames', {'GenoStr', 'IBI','BP','NB','SPB'});
% uniqueGenoTypes = unique(newTable.GenoStr);
% for i = 1:length(uniqueGenoTypes)
%     genoStr = uniqueGenoTypes{i};
%     genotypeFilter = strcmp(newTable.GenoStr, genoStr);
%     colDataIBI = horzcat(newTable{genotypeFilter, 'IBI'}{:});
%     colDataBP = horzcat(newTable{genotypeFilter, 'BP'}{:});
%     colDataNB = horzcat(newTable{genotypeFilter, 'NB'}{:});
%     colDataSPB = horzcat(newTable{genotypeFilter, 'SPB'}{:});
%     meanIBI = mean(colDataIBI,2);
%     stdIBI=std(colDataIBI,0,2);
%     % Calculate mean and standard deviation for BP
%     meanBP = mean(colDataBP, 2);
%     stdBP = std(colDataBP, 0, 2);
%     
%     % Calculate mean and standard deviation for NB
%     meanNB = mean(colDataNB, 2);
%     stdNB = std(colDataNB, 0, 2);
%     
%     % Calculate mean and standard deviation for SPB
%     meanSPB = mean(colDataSPB, 2);
%     stdSPB = std(colDataSPB, 0, 2);
% 
%     grid on;
%     hold on;
%     subplot(3,2,1);
%     % Prepare X and Y values for the shaded area (std deviation)
%     xFill = [data.(parameter);flipud(data.(parameter))];
%     yFill = [meanIBI - 2* stdIBI; flipud(meanIBI + 2* stdIBI)];
%     grid on;
%     hold on;
%     % Plot the standard deviation as a shaded area
%     fillHandle = fill(xFill, yFill, genoColorMap(genoStr), 'LineStyle', 'none', 'FaceAlpha', 0.1);
%     subplot(3,2,2);
%     xFillBP = [data.(parameter); flipud(data.(parameter))];
%     yFillBP = [meanBP - 2*stdBP; flipud(meanBP + 2*stdBP)];
%     fill(xFillBP, yFillBP, genoColorMap(genoStr), 'LineStyle', 'none', 'FaceAlpha', 0.1);
%     grid on;
%     hold on;
%     subplot(3,2,3);
%     xFillNB = [data.(parameter); flipud(data.(parameter))];
%     yFillNB = [meanNB - 2*stdNB; flipud(meanNB +2* stdNB)];
%     fill(xFillNB, yFillNB, genoColorMap(genoStr), 'LineStyle', 'none', 'FaceAlpha', 0.1);
%     
%     subplot(3,2,4);
%     xFillSPB = [data.(parameter); flipud(data.(parameter))];
%     yFillSPB = [meanSPB - 2*stdSPB; flipud(meanSPB + 2*stdSPB)];
%     fill(xFillSPB, yFillSPB, genoColorMap(genoStr), 'LineStyle', 'none', 'FaceAlpha', 0.1);
% 
% end

uniqueGenoTypes = unique(newTable.GenoStr);
for i = 1:length(uniqueGenoTypes)
    genoStr = uniqueGenoTypes{i};
    genotypeFilter = strcmp(newTable.GenoStr, genoStr);
    colDataIBI = horzcat(newTable{genotypeFilter, 'IBI'}{:});
    colDataBP = horzcat(newTable{genotypeFilter, 'BP'}{:});
    colDataNB = horzcat(newTable{genotypeFilter, 'NB'}{:});
    colDataSPB = horzcat(newTable{genotypeFilter, 'SPB'}{:});

    % Calculate mean for each parameter
    meanIBI = mean(colDataIBI, 2);
    meanBP = mean(colDataBP, 2);
    meanNB = mean(colDataNB, 2);
    meanSPB = mean(colDataSPB, 2);

    % Calculate IQR for each parameter
    iqrIBI = iqr(colDataIBI,2);
    iqrBP = iqr(colDataBP,2);
    iqrNB = iqr(colDataNB,2);
    iqrSPB = iqr(colDataSPB,2);

    % Calculate 95th percentile of IQR
    upperLimitIBI = prctile(iqrIBI, 95,2);
    upperLimitBP = prctile(iqrBP, 95,2);
    upperLimitNB = prctile(iqrNB, 95,2);
    upperLimitSPB = prctile(iqrSPB, 95,2);
    % Calculate the 5th percentile of IQR as the lower limit for IBI
    lowerLimitIBI = prctile(iqrIBI, 5, 2);
    
    % Calculate the 5th percentile of IQR as the lower limit for BP
    lowerLimitBP = prctile(iqrBP, 5, 2);
    
    % Calculate the 5th percentile of IQR as the lower limit for NB
    lowerLimitNB = prctile(iqrNB, 5, 2);
    
    % Calculate the 5th percentile of IQR as the lower limit for SPB
    lowerLimitSPB = prctile(iqrSPB, 5, 2);

    

    grid on;
    hold on;

    % Plot for IBI
    subplot(3,2,1);
    xFillIBI = [data.(parameter); flipud(data.(parameter))];
    yFillIBI = [meanIBI - lowerLimitIBI; flipud(meanIBI + upperLimitIBI)];
    fill(xFillIBI, yFillIBI, genoColorMap(genoStr), 'LineStyle', 'none', 'FaceAlpha', 0.25);

    % Plot for BP
    subplot(3,2,2);
    xFillBP = [data.(parameter); flipud(data.(parameter))];
    yFillBP = [meanBP - lowerLimitBP; flipud(meanBP + upperLimitBP)];
    fill(xFillBP, yFillBP, genoColorMap(genoStr), 'LineStyle', 'none', 'FaceAlpha', 0.25);

    % Plot for NB
    subplot(3,2,3);
    xFillNB = [data.(parameter); flipud(data.(parameter))];
    yFillNB = [meanNB - lowerLimitNB; flipud(meanNB + upperLimitNB)];
    fill(xFillNB, yFillNB, genoColorMap(genoStr), 'LineStyle', 'none', 'FaceAlpha', 0.25);

    % Plot for SPB
    subplot(3,2,4);
    xFillSPB = [data.(parameter); flipud(data.(parameter))];
    yFillSPB = [meanSPB - lowerLimitSPB; flipud(meanSPB + upperLimitSPB)];
    fill(xFillSPB, yFillSPB, genoColorMap(genoStr), 'LineStyle', 'none', 'FaceAlpha', 0.25);
end


%exportFile = [opDir 'paramsCompare_singlePlot.pdf']; %folder and filename for raster figures
%exportgraphics(fig2, exportFile ,'Resolution',150)

%fprintf("%s compare saved at %s", parameter, opDir);
%output_arg = append(parameter, " compare saved at ", opDir);