function plotRasterNetwork(networkAct,networkStats,relativeSpikeTimes,locData,binSize, ...
    opDir,chipWellFolder,xlimNetwork,ylimNetwork,textString,plotFileName,baselineFiringRate, plotTitle)
    
    %% NO INPUT PARSER - Simpler for PARFOR
    % Direct color settings
    RasterColor = '#000000';
    NetworkColor = '#0066CC';
    BaseLineColor = '#FF6600';
    ThresholdColor = '#C0392B';
    PeakColor = '#FF3333';
    LineWidth = 1.2;
    SpikeHeight = 0.9;
    SvgPlot = false;
    
    % Create figure
    f = figure('Color','white', 'Position', [0 0 800 600], 'Visible', 'off');
    set(f, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
    
    %% RASTER PLOT
    subplot(2,1,1);
    
    spikeTimes = single(relativeSpikeTimes.time);
    channels = single(relativeSpikeTimes.channel);
    nSpikes = length(spikeTimes);
    
    if nSpikes > 0
        x_plot = nan(3*nSpikes, 1, 'single');
        y_plot = nan(3*nSpikes, 1, 'single');
        indices = 1:3:3*nSpikes;
        x_plot(indices) = spikeTimes;
        x_plot(indices+1) = spikeTimes;
        half_height = SpikeHeight / 2;
        y_plot(indices) = channels - half_height;
        y_plot(indices+1) = channels + half_height;
        plot(x_plot, y_plot, 'Color', RasterColor, 'LineWidth', LineWidth);
    end
    
    ylabel('Channel', 'FontSize', 10, 'Color', [0.2 0.2 0.2]);

    % Title Text Added (SS) --- 
    titleText = {plotTitle, 'Neural Activity Raster'};  % Two-line title
    title(titleText, 'Interpreter', 'none', 'FontSize', 10, 'FontWeight', 'bold');

    %title('Neural Activity Raster', 'FontSize', 11, 'FontWeight', 'bold');
    box off; set(gca, 'YDir', 'normal', 'FontSize', 9, 'TickDir', 'out');
    set(gca, 'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2]);
    set(gca, 'XTickLabel', []);
    
    if ~isempty(channels)
        ylim([min(channels) - 1, max(channels) + 1]);
    end
    
    %% NETWORK ACTIVITY
    subplot(2,1,2);
    
    plot(networkAct.time, networkAct.firingRate, 'Color', NetworkColor, 'LineWidth', LineWidth*1.2);
    hold on;
    
    % Threshold
    if ~isempty(networkStats.threshold)
        switch upper(networkStats.thresholdFunction)
            case {'FIXED', 'RMS'}
                t_limits = xlim;
                plot(t_limits, [networkStats.threshold networkStats.threshold], ...
                     '--', 'Color', ThresholdColor, 'LineWidth', LineWidth);
        end
    end
    
    % Peaks
    plot(networkStats.maxAmplitudesTimes, networkStats.maxAmplitudesValues, 'o', ...
         'Color', PeakColor, 'MarkerSize', 4, 'MarkerFaceColor', 'white', 'LineWidth', 1);
    
    % Baseline
    plot([0 300], [baselineFiringRate baselineFiringRate], '--', 'Color', BaseLineColor, 'LineWidth', LineWidth);
    
    xlabel('Time (s)', 'FontSize', 10, 'Color', [0.2 0.2 0.2]);
    ylabel('Firing Rate (Hz)', 'FontSize', 10, 'Color', [0.2 0.2 0.2]);
    title('Network Activity', 'FontSize', 11, 'FontWeight', 'bold');
    box off; set(gca, 'FontSize', 9, 'TickDir', 'out');
    set(gca, 'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2]);
    
    if ~isempty(networkAct.firingRate)
        y_max = max([max(networkAct.firingRate), max(networkStats.maxAmplitudesValues)]);
        ylim([0, y_max * 1.1]);
    end
    hold off;
    
    %% PARFOR-SAFE SAVE OPERATIONS
    try
        % Create local copy of location data to avoid broadcast variable issues
        locDataLocal = struct();
        if isstruct(locData)
            locDataLocal = locData;  % Copy to make it loop-local
        end
        
        % Save spike times - should work (loop variable)
        fileNameBase = fullfile(chipWellFolder, append(plotFileName,'_spiketimes'));
        save(fileNameBase, '-struct', 'relativeSpikeTimes', '-v7.3');
        
        % Save location data - use local copy
        locFileName = fullfile(chipWellFolder, append(plotFileName,'_locs'));
        save(locFileName, '-struct', 'locDataLocal', '-v7.3');
        
    catch ME
        warning('Save failed: %s', ME.message);
    end
    
    %% EXPORT TIME WINDOWS
    % 60s window
    subplot(2,1,1); xlim([0 60]);
    if ~isempty(channels)
        ylim([min(channels) - 1, max(channels) + 1]);
    end
    subplot(2,1,2); xlim([0 60]); ylim([0 ylimNetwork]);
    
    fileNameBase = fullfile(chipWellFolder, 'Plot60s', 'png', plotFileName);
    print(f, [fileNameBase '.png'], '-dpng', '-r300');
    if SvgPlot
        fileNameBase = fullfile(chipWellFolder, 'Plot60s', 'svg', plotFileName);
        print(f, [fileNameBase '.svg'], '-dsvg', '-r300');
    end
    
    % 120s window
    subplot(2,1,1); xlim([0 120]);
    if ~isempty(channels)
        ylim([min(channels) - 1, max(channels) + 1]);
    end
    subplot(2,1,2); xlim([0 120]); ylim([0 ylimNetwork]);
    
    fileNameBase = fullfile(chipWellFolder, 'Plot120s', 'png', plotFileName);
    print(f, [fileNameBase '.png'], '-dpng', '-r300');
    if SvgPlot
        fileNameBase = fullfile(chipWellFolder, 'Plot120s', 'svg', plotFileName);
        print(f, [fileNameBase '.svg'], '-dsvg', '-r300');
    end
    
    % 300s window
    subplot(2,1,1); xlim([0 300]);
    if ~isempty(channels)
        ylim([min(channels) - 1, max(channels) + 1]);
    end
    subplot(2,1,2); xlim([0 300]); ylim([0 ylimNetwork]);
    
    annotation('textbox', [0.15, 0.48, 0.7, 0.06], 'String', textString, ...
        'EdgeColor', 'none', 'BackgroundColor', 'white', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 8, 'Color', [0.2 0.2 0.2], 'Margin', 2);
    
    fileNameBase = fullfile(chipWellFolder, 'Plot300s', 'png', plotFileName);
    print(f, [fileNameBase '.png'], '-dpng', '-r300');
    if SvgPlot
        fileNameBase = fullfile(chipWellFolder, 'Plot300s', 'svg', plotFileName);
        print(f, [fileNameBase '.svg'], '-dsvg', '-r300');
    end
    
    close(f);
end

