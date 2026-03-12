function plotRasterNetwork(networkAct,networkStats,relativeSpikeTimes,locData,binSize, ...
    opDir,chipWellFolder,xlimNetwork,ylimNetwork,textString,plotFileName,baselineFiringRate, plotTitle, plotMode)
    
    %% NO INPUT PARSER - Simpler for PARFOR
    % Direct color settings
    RasterColor = '#000000';
    NetworkColor = '#B22222';     % firebrick, matching IPNAnalysis plot_clean_network
    BaseLineColor = '#FF6600';
    ThresholdColor = '#C0392B';
    PeakColor = 'red';
    LineWidth = 1.2;
    SpikeHeight = 0.9;
    SvgPlot = true;
    if nargin < 14 || isempty(plotMode)
        plotMode = 'separate';
    end
    mergedMode = strcmpi(plotMode, 'merged');
    
    % Create figure
    f = figure('Color','white', 'Position', [0 0 800 600], 'Visible', 'off');
    set(f, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
    if mergedMode
        ax1 = axes(f, 'Position', [0.10 0.12 0.78 0.78]);
        ax2 = axes(f, 'Position', ax1.Position, 'Color', 'none', ...
            'YAxisLocation', 'right', ...
            'Box', 'off', 'FontSize', 9, 'TickDir', 'out');
        linkaxes([ax1, ax2], 'x');
    else
        tl = tiledlayout(f, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        ax1 = nexttile(tl);
        ax2 = nexttile(tl);
    end

    %% RASTER PLOT
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
        plot(ax1, x_plot, y_plot, 'Color', RasterColor, 'LineWidth', LineWidth);
    end
    
    ylabel(ax1, 'Channel', 'FontSize', 10, 'Color', [0.2 0.2 0.2]);

    % Title Text Added (SS) --- 
    if mergedMode
        titleText = {plotTitle, 'Neural Activity Raster + Network Activity'};
    else
        titleText = {plotTitle, 'Neural Activity Raster'};  % Two-line title
    end
    title(ax1, titleText, 'Interpreter', 'none', 'FontSize', 10, 'FontWeight', 'bold');

    %title('Neural Activity Raster', 'FontSize', 11, 'FontWeight', 'bold');
    box off; set(ax1, 'YDir', 'normal', 'FontSize', 9, 'TickDir', 'out');
    set(ax1, 'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2]);
    if ~mergedMode
        set(ax1, 'XTickLabel', []);
    end
    
    if ~isempty(channels)
        ylim(ax1, [min(channels) - 1, max(channels) + 1]);
    end
    
    %% NETWORK ACTIVITY
    plot(ax2, networkAct.time, networkAct.firingRate, 'Color', NetworkColor, 'LineWidth', LineWidth*1.2);
    hold(ax2, 'on');
    
    % Threshold
    if ~isempty(networkStats.threshold)
        switch upper(networkStats.thresholdFunction)
            case {'FIXED', 'RMS'}
                t_limits = xlim(ax2);
                plot(ax2, t_limits, [networkStats.threshold networkStats.threshold], ...
                     '--', 'Color', ThresholdColor, 'LineWidth', LineWidth);
        end
    end
    
    % Peaks
    plot(ax2, networkStats.maxAmplitudesTimes, networkStats.maxAmplitudesValues, 'o', ...
         'Color', PeakColor, 'MarkerSize', 4, 'MarkerFaceColor', 'white', 'LineWidth', 1);
    
    % Baseline
    plot(ax2, [0 300], [baselineFiringRate baselineFiringRate], '--', 'Color', BaseLineColor, 'LineWidth', LineWidth);
    
    if mergedMode
        xlabel(ax1, 'Time (s)', 'FontSize', 10, 'Color', [0.2 0.2 0.2]);
        ylabel(ax2, 'Firing Rate (Hz)', 'FontSize', 10, 'Color', [0.2 0.2 0.2]);
        set(ax2, 'XTick', [], 'XColor', 'none', 'YColor', [0.2 0.2 0.2]);
    else
        xlabel(ax2, 'Time (s)', 'FontSize', 10, 'Color', [0.2 0.2 0.2]);
        ylabel(ax2, 'Firing Rate (Hz)', 'FontSize', 10, 'Color', [0.2 0.2 0.2]);
        title(ax2, 'Network Activity', 'FontSize', 11, 'FontWeight', 'bold');
        set(ax2, 'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2]);
    end
    box off; set(ax2, 'FontSize', 9, 'TickDir', 'out');
    
    if ~isempty(networkAct.firingRate)
        y_max = max([max(networkAct.firingRate), max(networkStats.maxAmplitudesValues)]);
        ylim(ax2, [0, y_max * 1.1]);
    end
    hold(ax2, 'off');
    
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
    xlim(ax1, [0 60]);
    if ~isempty(channels)
        ylim(ax1, [min(channels) - 1, max(channels) + 1]);
    end
    xlim(ax2, [0 60]); ylim(ax2, [0 ylimNetwork]);
    
    fileNameBase = fullfile(chipWellFolder, 'Plot60s', 'png', plotFileName);
    print(f, [fileNameBase '.png'], '-dpng', '-r300');
    if SvgPlot
        fileNameBase = fullfile(chipWellFolder, 'Plot60s', 'svg', plotFileName);
        print(f, [fileNameBase '.svg'], '-dsvg', '-painters');
    end
    
    % 120s window
    xlim(ax1, [0 120]);
    if ~isempty(channels)
        ylim(ax1, [min(channels) - 1, max(channels) + 1]);
    end
    xlim(ax2, [0 120]); ylim(ax2, [0 ylimNetwork]);
    
    fileNameBase = fullfile(chipWellFolder, 'Plot120s', 'png', plotFileName);
    print(f, [fileNameBase '.png'], '-dpng', '-r300');
    if SvgPlot
        fileNameBase = fullfile(chipWellFolder, 'Plot120s', 'svg', plotFileName);
        print(f, [fileNameBase '.svg'], '-dsvg', '-painters');
    end
    
    % 300s window
    xlim(ax1, [0 300]);
    if ~isempty(channels)
        ylim(ax1, [min(channels) - 1, max(channels) + 1]);
    end
    xlim(ax2, [0 300]); ylim(ax2, [0 ylimNetwork]);
    
    if mergedMode
        textY = 0.02;
    else
        textY = 0.48;
    end
    annotation('textbox', [0.15, textY, 0.7, 0.06], 'String', textString, ...
        'EdgeColor', 'none', 'BackgroundColor', 'white', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 8, 'Color', [0.2 0.2 0.2], 'Margin', 2);
    
    fileNameBase = fullfile(chipWellFolder, 'Plot300s', 'png', plotFileName);
    print(f, [fileNameBase '.png'], '-dpng', '-r300');
    if SvgPlot
        fileNameBase = fullfile(chipWellFolder, 'Plot300s', 'svg', plotFileName);
        print(f, [fileNameBase '.svg'], '-dsvg', '-painters');
    end
    
    close(f);
end
