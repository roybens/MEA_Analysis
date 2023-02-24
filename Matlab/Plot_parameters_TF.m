%plot different parameter values based on gaussian


parentFolderPath = 'C:/Users/Tim/Documents/MATLAB/Maxwell Analysis Stuff/MEA_DataOutput_syngap3_DIV6/';


wt = [16657,16665,16856,16795,16666,16861,16867,16873];
het = [16742,16821,16850,16787,16792,16940,16869,16880];

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(parentFolderPath, '*.csv'); 
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    %fullFileName = fullfile(theFiles(k).folder, baseFileName);
    pathFileNetwork = fullfile(theFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', pathFileNetwork);
    
    datafile = (pathFileNetwork);
    data = readtable(datafile,'PreserveVariableNames',true);
    
    extractChipID = regexp(pathFileNetwork,'\d{5}?','match');

    idDouble = str2double(extractChipID);
    if ismember(idDouble,wt)
        geno = 'WT';
        genoStr = string(geno);
    else 
        geno = 'het';
        genoStr = string(geno);
    end    

    
    fig = figure('color','w','Position',[0 0 1600 800],'Visible','off');
    subplot(2,2,1);
    plot(data.Gaussian,data.IBI);
    title(genoStr + ' IBI')
    xlabel('Gaussian')
    xticks(0:.04:1)
    ylim([0 10])
    grid on

    subplot(2,2,2);
    plot(data.Gaussian,data.(4));
    title(genoStr +' Burst Peak')
    xlabel('Gaussian')
    xticks(0:.04:1)
    ylim([0 20])
    grid on

    subplot(2,2,3);
    plot(data.Gaussian,data.(5));
    title(genoStr + ' # of Bursts')
    xlabel('Gaussian')
    xticks(0:.04:1)
    ylim([0 800])
    grid on

    subplot(2,2,4);
    plot(data.Gaussian,data.(6));
    title(genoStr + ' Spikes per Burst')
    xlabel('Gaussian')
    xticks(0:.04:1)
    ylim([0 150000])
    grid on
    
    exportFile = [parentFolderPath 'paramsCompare.pdf']; %folder and filename for raster figures
    exportgraphics(fig, exportFile ,'Append',true,'Resolution',150)
end
%% Plot all lines on one plot


fig2 = figure('color','w','Position',[0 0 1600 800],'Visible','off');


for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    %fullFileName = fullfile(theFiles(k).folder, baseFileName);
    pathFileNetwork = fullfile(theFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', pathFileNetwork);
    
    datafile = (pathFileNetwork);
    data = readtable(datafile,'PreserveVariableNames',true);
    
    extractChipID = regexp(pathFileNetwork,'\d{5}?','match');

    idDouble = str2double(extractChipID);
    if ismember(idDouble,wt)
        geno = 'WT';
        genoStr = string(geno);
        color = '-k';
    else 
        geno = 'het';
        genoStr = string(geno);
        color = '--r';
    end    

    
    subplot(2,2,1);
    plot(data.Gaussian,data.IBI,color);
    title('IBI')
    xlabel('Gaussian')
    xticks(0:.04:1)
    ylim([0 10])
    grid on
    hold on

    subplot(2,2,2);
    plot(data.Gaussian,data.(4),color);
    title('Burst Peak')
    xlabel('Gaussian')
    xticks(0:.04:1)
    ylim([0 20])
    grid on
    hold on

    subplot(2,2,3);
    plot(data.Gaussian,data.(5),color);
    title('# of Bursts')
    xlabel('Gaussian')
    xticks(0:.04:1)
    ylim([0 800])
    grid on
    hold on

    subplot(2,2,4);
    plot(data.Gaussian,data.(6),color);
    title('Spikes per Burst')
    xlabel('Gaussian')
    xticks(0:.04:1)
    ylim([0 150000])
    grid on
    hold on
    
end
exportFile = [parentFolderPath 'paramsCompare_singlePlot.pdf']; %folder and filename for raster figures
exportgraphics(fig2, exportFile ,'Resolution',150)

%% Read all CSVs and add means of specified gaussian values to single table

gaussian = .16;

for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    %fullFileName = fullfile(theFiles(k).folder, baseFileName);
    pathFileNetwork = fullfile(theFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', pathFileNetwork);
    
    datafile = (pathFileNetwork);
    data = readtable(datafile,'PreserveVariableNames',true);
    
    extractChipID = regexp(pathFileNetwork,'\d{5}?','match');

    idDouble = str2double(extractChipID);
    if ismember(idDouble,wt)
        geno = 'WT';
        genoStr = string(geno);
    else 
        geno = 'het';
        genoStr = string(geno);
    end 

    selectedParam = ismember(guassian,data.Gaussian,'rows');

%     if data(:,2) == gaussian
%         %if any row in column 2 = preset gaussian, then append that row to
%         %new table 
%     end

end


