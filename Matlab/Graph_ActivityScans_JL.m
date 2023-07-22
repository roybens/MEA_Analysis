clear
close all

%% settings
% set path to folder containing subfolders that contain h5 files
parentFolderPath = '/mnt/harddrive-2/Organoids_Mandeep_Fink_Lab';
% set output folder
outputFolderPath = '/home/jonathan/Documents/Scripts/Matlab/scrpits_output/Organoids/activityscan_outputs';

%% interate through ActivityScans
% get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(parentFolderPath, '**/ActivityScan/**/*.mxassay'); 
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
    pathFileActivityScan = fullfile(theFiles(k).folder);
    fprintf(1, 'Now reading %s\n',  pathFileActivityScan);
    tic

    %% extract dir information
    fileDirParts = strsplit(pathFileActivityScan, filesep); % split dir into elements
    scan_chipID = fileDirParts{end-2}; % extract chipID
    scan_runID = fileDirParts{end-0}; % extract runID
    
    %% extract hd5 information
    wellID = 1; % select which well to analyze
    % create fileManager object for the Activity Scan
    activityScanData = mxw.fileManager(pathFileActivityScan,wellID);
    
    %% Plotting the Spike Data
    %% I. Raster Plots
    recNr = 3;
    % sampling frequency
    fsActivity = activityScanData.fileObj(recNr).samplingFreq;
    % get spike time stamps 
    tsActivity = double(activityScanData.fileObj(recNr).spikes.frameno - activityScanData.fileObj(recNr).firstFrameNum)/fsActivity;
    % channel list, where spike time stamps where detected
    chActivity = activityScanData.fileObj(recNr).spikes.channel;
    
    % plot raster
    figure('Color','w','visible','off');
    subplot(2,1,1)
    plot(tsActivity, chActivity,'.','Color','#135ba3','MarkerSize',2)
    box off 
    xlabel('Time [s]') 
    ylabel('Channel')
    title('Activity Scan Raster Plot')
    xmaxActivity = max(tsActivity);
    xlim([0 xmaxActivity])
    ylim([0 max(chActivity)])
    
    saveas(gcf,append(outputFolderPath,'/activityscan_Raster',scan_runID,'.png'))

    %% II. Spike Rate Activity Map and Distribution
    % get the mean firing rate for each electrode
    meanFiringRate = mxw.activityMap.computeSpikeRate(activityScanData);
    
    % define maximum firing rate used in plots as 99th percentile of firing
    % rate values
    maxFiringRate = mxw.util.percentile(meanFiringRate(meanFiringRate~=0),99);
    
    % plot the list of mean firing rates as a map, given the electrode 
    % x and y coordinates in 'fileManagerObj.processedMap'
    figure('Color','w','visible','off');
    subplot(2,1,1);
    mxw.plot.activityMap(activityScanData, meanFiringRate, 'Ylabel', '[Hz]',...
        'CaxisLim', [0.1 max(meanFiringRate)/5],'Figure',false,'Title','Firing Rate Activity Map');
    % run the line above several times, experimenting with different 
    % [min max] range of the color gradient 'CaxisLim'
    
    % add scale bar 
    line([300 800],[2000+400 2000+400],'Color','k','LineWidth',4);
    axis off;
    text(340,2100+500,'0.5 mm','color','k');
    xlim([200 3750]);ylim([150 2500])
    
    % plot the distribution of firing rates across all of the electrodes
    subplot(2,1,2);
    thrFiringRate = 0.1; % set a minimum spike rate threshold (Hz)
    histogram(meanFiringRate(meanFiringRate>thrFiringRate),0:.1:ceil(maxFiringRate), ...
        'FaceColor','#135ba3', "EdgeColor",'#414042')
    xlim([thrFiringRate maxFiringRate])
    ylabel('Counts');xlabel('Firing Rate [Hz]');
    box off;
    legend(['Mean Firing Rate = ',num2str(mean(meanFiringRate(meanFiringRate>thrFiringRate)),'%.2f'),...
        ' Hz,  sd = ',num2str(std(meanFiringRate(meanFiringRate>thrFiringRate)),'%.2f')])

    saveas(gcf,append(outputFolderPath,'/activityscan_FiringRateMap',scan_runID,'.png'))

    %% III. Spike Amplitude Activity Map and Distribution
    % get th 90th percentile spike amplitude value for each electrode
    amplitude90perc = abs(mxw.activityMap.computeAmplitude90percentile(activityScanData));
    
    % define maximum amplitude used in plots as 99th percentile of amplitude values
    maxAmp = mxw.util.percentile(amplitude90perc(amplitude90perc~=0),99);
    
    % plot the mean firing rate vector as a map, given the electrode 
    % x and y coordinates in 'fileManagerObj.processedMap'
    figure('Color','w','visible','off');
    subplot(2,1,1);
    mxw.plot.activityMap(activityScanData, amplitude90perc,'Ylabel', '[\muV]',...
        'CaxisLim', [10 maxAmp], 'Figure',false,'Title','Spike Amplitude Activity Map');
    % run the line above several times, experimenting with different 
    % [min max] range of the color gradient 'CaxisLim'
    
    % add scale bar 
    line([300 800],[2000+400 2000+400],'Color','k','LineWidth',4);
    axis off;
    text(340,2100+500,'0.5 mm','color','k');
    xlim([200 3750])
    ylim([150 2500])
    
    % plot the distribution of spike amplitudes across all of the electrodes
    subplot(2,1,2);
    thrAmp = 10;  % set a minimum spike amplitude threshold (uV)
    histogram(amplitude90perc(amplitude90perc>thrAmp),ceil(0:1:maxAmp), ...
        'FaceColor','#135ba3', "EdgeColor",'#414042')
    xlim([thrAmp maxAmp])
    ylabel('Counts');xlabel('Spike Amplitude [\muV]');
    box off;
    legend(['Mean Spike Amplitude = ',num2str(mean(amplitude90perc(amplitude90perc>thrAmp)),'%.2f'),...
        ' \muV,  sd = ',num2str(std(amplitude90perc(amplitude90perc>thrAmp)),'%.2f')])

    saveas(gcf,append(outputFolderPath,'/activityscan_AmplitudeMap',scan_runID,'.png'))
end
