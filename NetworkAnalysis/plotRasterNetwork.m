function plotRasterNetwork(networkData,networkAct,networkStats,relativeSpikeTimes,binSize, ...
    opDir,chipWellFolder,xlimNetwork,ylimNetwork,textString,plotFileName,epsPlot,baselineFiringRate,logBursts)
    % Create figure in the background
    f = figure('Color','w', 'Position', [0 0 400 800], 'Visible', 'off');  
    % Plot raster and network activity
    subplot(2,1,1);
    mxw.plot.rasterPlot(networkData, 'Figure', false);
    box off;
    xlim([0 60]);
    ylim([1 max(relativeSpikeTimes.channel)]);
    subplot(2,1,2);
    plotNetworkActivityModified(networkAct, 'ThresholdFunction',networkStats.thresholdFunction, 'Threshold',networkStats.threshold, 'Figure', false);
    box off;hold on;    
    plot(networkStats.maxAmplitudesTimes, networkStats.maxAmplitudesValues, 'or');
    %fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', 'Figureformat',plotFileName);
    fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity',plotFileName);
    savefig(f, [fileNameBase '.fig']);
    ylim([0 ylimNetwork]);
    xlim([0 xlimNetwork])
    %rohan made changes here
    fileNameBase = fullfile(chipWellFolder, 'Plot60s', 'png', plotFileName);
    %fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', chipWellFolder, 'Plot60s', 'png', plotFileName);
    print(f, [fileNameBase '.png'], '-dpng', '-r300');
    if epsPlot
    %fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', chipWellFolder, 'Plot60s', 'eps', plotFileName);
    fileNameBase = fullfile(chipWellFolder, 'Plot60s', 'eps', plotFileName);
    print(f, [fileNameBase '.eps'], '-depsc', '-r300');
    end
    %totalTime = toc;  % Measure the total elapsed time after the loop
    %fprintf('in Well %d Total elapsed time for plot and save for 60s: %f seconds\n', wellID,totalTime); 
    subplot(2,1,1);
    xlim([0 120])
    ylim([1 max(relativeSpikeTimes.channel)])
    subplot(2,1,2);
    xlim([0 120])
    ylim([0 ylimNetwork])
    %fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', 'Plot120s',plotFileName); 
    %fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', chipWellFolder, 'Plot120s', 'png', plotFileName);
    fileNameBase = fullfile(chipWellFolder, 'Plot120s', 'png', plotFileName);
    print(f, [fileNameBase '.png'], '-dpng', '-r300');
    if epsPlot
    %fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', chipWellFolder, 'Plot120s', 'eps', plotFileName);
    fileNameBase = fullfile(chipWellFolder, 'Plot120s', 'eps', plotFileName);
    print(f, [fileNameBase '.eps'], '-depsc', '-r300');
    end
    subplot(2,1,1);
    xlim([0 300])
    ylim([0 max(relativeSpikeTimes.channel)])
    subplot(2,1,2);
    xlim([0 300])
    ylim([0 ylimNetwork])
     % Assuming you want to display metrics above the raster plot
    timeVector = 0:binSize:max(relativeSpikeTimes.time);   %need to optimize
    plot(baselineFiringRate*ones(ceil(timeVector(end)),1))
    % Create one annotation box containing all the text entries
    % Adjust the position vector [x y width height] as needed
    annotation('textbox', [0.1, 0.425, 0.9, 0.1], 'String', textString, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    %fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', 'Plot300s',plotFileName);
        % --- after drawing the 300s raster (subplot 2,1,1) and before saving ---
   
    %fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', chipWellFolder, 'Plot300s', 'png', plotFileName);
    fileNameBase = fullfile(chipWellFolder, 'Plot300s', 'png', plotFileName);
    print(f, [fileNameBase '.png'], '-dpng', '-r300');

    if epsPlot
   % fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', chipWellFolder, 'Plot300s', 'eps', plotFileName);
    fileNameBase = fullfile(chipWellFolder, 'Plot300s', 'eps', plotFileName);
    print(f, [fileNameBase '.eps'], '-depsc', '-r300');
    end    

    

    set(f,'Position',[0 0 1000 800])
    subplot(2,1,1);
    xlim([0 30])
    ylim([0 max(relativeSpikeTimes.channel)])
    subplot(2,1,2);
    xlim([0 30])
    ylim([0 ylimNetwork])
    if exist('logBursts','var') && ~isempty(logBursts)
        ax = subplot(2,1,1);
        hold(ax,'on');
        yMax = max(relativeSpikeTimes.channel)+1;
        for i = 1:numel(logBursts.start_time)
            rectangle(ax, ...
                'Position', [logBursts.start_time(i), 0.5, logBursts.duration_s(i), yMax], ...
                'FaceColor', [0.6 0 0 0.4], ...
                'EdgeColor', 'none');
        end
        hold(ax,'off');
    end
    %fileNameBase = fullfile(opDir, 'Network_outputs', 'Raster_BurstActivity', 'Plot600s',plotFileName);
    fileNameBase = fullfile(chipWellFolder, 'Plot600s', 'png', plotFileName);
        % Save the figure in PNG format
    print(f, [fileNameBase '.png'], '-dpng', '-r300');
    
    if epsPlot
    fileNameBase = fullfile(chipWellFolder, 'Plot600s', 'eps', plotFileName);
    print(f, [fileNameBase '.eps'], '-depsc', '-r300');
    end
    

    close(f);
end