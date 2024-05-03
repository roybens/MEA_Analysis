close all;
clear all;

figDirectory= '/mnt/disk15tb/paula/Main_DA_Projects/data_analysis_output/SPTAN1_2_T2/Network_outputs/Raster_BurstActivity/300s';

opDirectory='/mnt/disk15tb/paula/Main_DA_Projects/data_analysis_output/SPTAN1_2_T2/Network_outputs/Raster_BurstActivity/new/'

ymin = 0;
ymax= 20;
newYLim = [ymin, ymax]; % Replace with your desired y-axis limits

% Create the destination directory if it doesn't exist

mkdir(opDirectory); 
mkdir(strcat(opDirectory,'60s/png'));
mkdir(strcat(opDirectory,'60s/svg'));
mkdir(strcat(opDirectory,'120s/png'));
mkdir(strcat(opDirectory,'120s/svg'));
mkdir(strcat(opDirectory,'300s/png'));
mkdir(strcat(opDirectory,'300s/svg'));



% Get a list of .fig files in the source directory
figFiles = dir(fullfile(figDirectory, '*.fig'));

% Process each .fig file
for i = 1:length(figFiles)
    % Load the figure
    fig = openfig(fullfile(figDirectory, figFiles(i).name), 'invisible');
    fprintf('Processing file: %s\n', figFiles(i).name);
    % Modify the y-axis limits
    ax = findall(fig, 'Type', 'axes');
   
    set(ax(1), 'YLim', newYLim);
    

    % Save the figure as .png and .svg in the destination directory
    [~, name, ~] = fileparts(figFiles(i).name);
    set(ax(1), 'XLim', [0, 60]);
    pngFile = fullfile(strcat(opDirectory,'60s/png'), [name, '.png']);
    svgFile = fullfile(strcat(opDirectory,'60s/svg'), [name, '.svg']);
    
    saveas(fig, pngFile);
    saveas(fig, svgFile);
    
    set(ax(1), 'XLim', [0, 120]);
    pngFile = fullfile(strcat(opDirectory,'120s/png'), [name, '.png']);
    svgFile = fullfile(strcat(opDirectory,'120s/svg'), [name, '.svg']);
    saveas(fig, pngFile);
    saveas(fig, svgFile);

    set(ax(1), 'XLim', [0 ,300]);
    pngFile = fullfile(strcat(opDirectory,'300s/png'), [name, '.png']);
    svgFile = fullfile(strcat(opDirectory,'300s/svg'), [name, '.svg']);
    saveas(fig, pngFile);
    saveas(fig, svgFile);
    
    % Close the figure
    close(fig);
end
