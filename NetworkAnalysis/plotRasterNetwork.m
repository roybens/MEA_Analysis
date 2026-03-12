function plotRasterNetwork(networkAct,networkStats,relativeSpikeTimes,locData,binSize, ...
    opDir,chipWellFolder,xlimNetwork,ylimNetwork,textString,plotFileName,baselineFiringRate, plotTitle, plotMode)

    %% SIMPLIFIED / ROBUST VERSION
    % API unchanged

    % ----------------------------
    % Visual settings
    % ----------------------------
    RasterColor    = [0 0 0];
    NetworkColor   = [178 34 34] / 255;   % firebrick
    BaseLineColor  = [1 102 0] / 255;
    ThresholdColor = [192 57 43] / 255;
    PeakEdgeColor  = [1 0 0];
    PeakFaceColor  = [1 1 1];

    RasterLineWidth  = 0.6;
    NetworkLineWidth = 1.8;
    SpikeHeight      = 0.7;

    SvgPlot = false;

    if nargin < 14 || isempty(plotMode)
        plotMode = 'separate';
    end
    mergedMode = strcmpi(plotMode, 'merged');

    % ----------------------------
    % Prepare spike data
    % ----------------------------
    if isfield(relativeSpikeTimes, 'time') && ~isempty(relativeSpikeTimes.time)
        spikeTimes = double(relativeSpikeTimes.time(:));
    else
        spikeTimes = [];
    end

    if isfield(relativeSpikeTimes, 'channel') && ~isempty(relativeSpikeTimes.channel)
        channelsRaw = double(relativeSpikeTimes.channel(:));
    else
        channelsRaw = [];
    end

    nSpikes = numel(spikeTimes);

    % ----------------------------
    % Remap channels to dense 1:N for display
    % This is for visualization only.
    % ----------------------------
    if ~isempty(channelsRaw)
        [uniqueChannels, ~, channelMapped] = unique(channelsRaw, 'stable');
        nDisplayChannels = numel(uniqueChannels);
    else
        uniqueChannels   = [];
        channelMapped    = [];
        nDisplayChannels = 0;
    end

    % ----------------------------
    % Create figure
    % ----------------------------
    f = figure('Color','white', 'Position', [100 100 900 650], 'Visible', 'off');
    set(f, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');

    if mergedMode
        ax = axes(f, 'Position', [0.10 0.12 0.82 0.78], ...
            'FontSize', 9, 'TickDir', 'out', 'Box', 'off');
    else
        tl = tiledlayout(f, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        axRaster = nexttile(tl);
        axNet    = nexttile(tl);
    end

    %% =========================================================
    % MERGED MODE  (single axes + yyaxis)
    % ==========================================================
    if mergedMode

        % ----------------------------
        % Left axis: raster
        % ----------------------------
        yyaxis(ax, 'left');
        hold(ax, 'on');

        if nSpikes > 0
            x_plot = nan(3*nSpikes, 1);
            y_plot = nan(3*nSpikes, 1);

            idx = 1:3:(3*nSpikes);
            half_height = SpikeHeight / 2;

            x_plot(idx)   = spikeTimes;
            x_plot(idx+1) = spikeTimes;

            y_plot(idx)   = channelMapped - half_height;
            y_plot(idx+1) = channelMapped + half_height;

            plot(ax, x_plot, y_plot, '-', 'Color', RasterColor, 'LineWidth', RasterLineWidth);
        end

        ylabel(ax, 'Active channel index', 'FontSize', 10, 'Color', [0.15 0.15 0.15]);
        ax.YColor = [0.15 0.15 0.15];

        if nDisplayChannels > 0
            ylim(ax, [0.5, nDisplayChannels + 0.5]);
        else
            ylim(ax, [0 1]);
        end

        % ----------------------------
        % Right axis: network activity
        % ----------------------------
        yyaxis(ax, 'right');
        hold(ax, 'on');

        if isfield(networkAct, 'time') && isfield(networkAct, 'firingRate') && ~isempty(networkAct.time)
            plot(ax, networkAct.time, networkAct.firingRate, ...
                'Color', NetworkColor, 'LineWidth', NetworkLineWidth);
        end

        % Threshold
        if isfield(networkStats, 'threshold') && ~isempty(networkStats.threshold) && ...
           isfield(networkStats, 'thresholdFunction') && ~isempty(networkStats.thresholdFunction)

            switch upper(string(networkStats.thresholdFunction))
                case {"FIXED","RMS"}
                    xLimitsNow = [0 300];
                    if ~isempty(networkAct.time)
                        xLimitsNow = [min(networkAct.time), max(networkAct.time)];
                    end
                    plot(ax, xLimitsNow, [networkStats.threshold networkStats.threshold], '--', ...
                        'Color', ThresholdColor, 'LineWidth', 1.1);
            end
        end

        % Peaks
        if isfield(networkStats, 'maxAmplitudesTimes') && isfield(networkStats, 'maxAmplitudesValues') && ...
           ~isempty(networkStats.maxAmplitudesTimes)
            plot(ax, networkStats.maxAmplitudesTimes, networkStats.maxAmplitudesValues, 'o', ...
                'MarkerEdgeColor', PeakEdgeColor, ...
                'MarkerFaceColor', PeakFaceColor, ...
                'MarkerSize', 4.5, 'LineWidth', 1);
        end

        % Baseline
        xBase = [0 300];
        if ~isempty(networkAct.time)
            xBase = [min(networkAct.time), max(networkAct.time)];
        end
        plot(ax, xBase, [baselineFiringRate baselineFiringRate], '--', ...
            'Color', BaseLineColor, 'LineWidth', 1.1);

        ylabel(ax, 'Firing Rate (Hz)', 'FontSize', 10, 'Color', [0.15 0.15 0.15]);
        ax.YColor = [0.15 0.15 0.15];

        if ~isempty(networkAct.firingRate)
            peakVals = [];
            if isfield(networkStats, 'maxAmplitudesValues') && ~isempty(networkStats.maxAmplitudesValues)
                peakVals = networkStats.maxAmplitudesValues(:);
            end
            yMaxNet = max([networkAct.firingRate(:); peakVals; baselineFiringRate; 1]);
            ylim(ax, [0, 1.1*yMaxNet]);
        else
            ylim(ax, [0, max(1, ylimNetwork)]);
        end

        xlabel(ax, 'Time (s)', 'FontSize', 10, 'Color', [0.15 0.15 0.15]);
        title(ax, {plotTitle, 'Neural Activity Raster + Network Activity'}, ...
            'Interpreter', 'none', 'FontSize', 11, 'FontWeight', 'bold');

        xlim(ax, [0 60]);
        grid(ax, 'off');

    %% =========================================================
    % SEPARATE MODE
    % ==========================================================
    else
        % ----------------------------
        % Raster subplot
        % ----------------------------
        hold(axRaster, 'on');

        if nSpikes > 0
            x_plot = nan(3*nSpikes, 1);
            y_plot = nan(3*nSpikes, 1);

            idx = 1:3:(3*nSpikes);
            half_height = SpikeHeight / 2;

            x_plot(idx)   = spikeTimes;
            x_plot(idx+1) = spikeTimes;

            y_plot(idx)   = channelMapped - half_height;
            y_plot(idx+1) = channelMapped + half_height;

            plot(axRaster, x_plot, y_plot, '-', 'Color', RasterColor, 'LineWidth', RasterLineWidth);
        end

        ylabel(axRaster, 'Active channel index', 'FontSize', 10, 'Color', [0.2 0.2 0.2]);
        title(axRaster, {plotTitle, 'Neural Activity Raster'}, ...
            'Interpreter', 'none', 'FontSize', 10, 'FontWeight', 'bold');

        box(axRaster, 'off');
        set(axRaster, 'YDir', 'normal', 'FontSize', 9, 'TickDir', 'out');
        set(axRaster, 'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2], 'XTickLabel', []);

        if nDisplayChannels > 0
            ylim(axRaster, [0.5, nDisplayChannels + 0.5]);
        else
            ylim(axRaster, [0 1]);
        end

        % ----------------------------
        % Network subplot
        % ----------------------------
        hold(axNet, 'on');

        if ~isempty(networkAct.time)
            plot(axNet, networkAct.time, networkAct.firingRate, ...
                'Color', NetworkColor, 'LineWidth', NetworkLineWidth);
        end

        % Threshold
        if isfield(networkStats, 'threshold') && ~isempty(networkStats.threshold) && ...
           isfield(networkStats, 'thresholdFunction') && ~isempty(networkStats.thresholdFunction)

            switch upper(string(networkStats.thresholdFunction))
                case {"FIXED","RMS"}
                    xLimitsNow = [0 300];
                    if ~isempty(networkAct.time)
                        xLimitsNow = [min(networkAct.time), max(networkAct.time)];
                    end
                    plot(axNet, xLimitsNow, [networkStats.threshold networkStats.threshold], '--', ...
                        'Color', ThresholdColor, 'LineWidth', 1.1);
            end
        end

        % Peaks
        if isfield(networkStats, 'maxAmplitudesTimes') && isfield(networkStats, 'maxAmplitudesValues') && ...
           ~isempty(networkStats.maxAmplitudesTimes)
            plot(axNet, networkStats.maxAmplitudesTimes, networkStats.maxAmplitudesValues, 'o', ...
                'MarkerEdgeColor', PeakEdgeColor, ...
                'MarkerFaceColor', PeakFaceColor, ...
                'MarkerSize', 4.5, 'LineWidth', 1);
        end

        % Baseline
        xBase = [0 300];
        if ~isempty(networkAct.time)
            xBase = [min(networkAct.time), max(networkAct.time)];
        end
        plot(axNet, xBase, [baselineFiringRate baselineFiringRate], '--', ...
            'Color', BaseLineColor, 'LineWidth', 1.1);

        xlabel(axNet, 'Time (s)', 'FontSize', 10, 'Color', [0.2 0.2 0.2]);
        ylabel(axNet, 'Firing Rate (Hz)', 'FontSize', 10, 'Color', [0.2 0.2 0.2]);
        title(axNet, 'Network Activity', 'FontSize', 11, 'FontWeight', 'bold');

        box(axNet, 'off');
        set(axNet, 'FontSize', 9, 'TickDir', 'out');
        set(axNet, 'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2]);

        if ~isempty(networkAct.firingRate)
            peakVals = [];
            if isfield(networkStats, 'maxAmplitudesValues') && ~isempty(networkStats.maxAmplitudesValues)
                peakVals = networkStats.maxAmplitudesValues(:);
            end
            yMaxNet = max([networkAct.firingRate(:); peakVals; baselineFiringRate; 1]);
            ylim(axNet, [0, 1.1*yMaxNet]);
        else
            ylim(axNet, [0, max(1, ylimNetwork)]);
        end
    end

    %% =========================================================
    % SAVE MAT FILES
    % ==========================================================
    try
        locDataLocal = struct();
        if isstruct(locData)
            locDataLocal = locData;
        end

        fileNameBase = fullfile(chipWellFolder, append(plotFileName,'_spiketimes'));
        save(fileNameBase, '-struct', 'relativeSpikeTimes', '-v7.3');

        locFileName = fullfile(chipWellFolder, append(plotFileName,'_locs'));
        save(locFileName, '-struct', 'locDataLocal', '-v7.3');

    catch ME
        warning('Save failed: %s', ME.message);
    end

    %% =========================================================
    % EXPORT WINDOWS
    % ==========================================================
    exportWindows = [60, 120, 300];

    for k = 1:numel(exportWindows)
        win = exportWindows(k);

        if mergedMode
            xlim(ax, [0 win]);

            yyaxis(ax, 'left');
            if nDisplayChannels > 0
                ylim(ax, [0.5, nDisplayChannels + 0.5]);
            else
                ylim(ax, [0 1]);
            end

            yyaxis(ax, 'right');
            ylim(ax, [0 ylimNetwork]);

        else
            xlim(axRaster, [0 win]);
            if nDisplayChannels > 0
                ylim(axRaster, [0.5, nDisplayChannels + 0.5]);
            else
                ylim(axRaster, [0 1]);
            end

            xlim(axNet, [0 win]);
            ylim(axNet, [0 ylimNetwork]);
        end

        % Add annotation only for 300 s export
        annHandle = [];
        if win == 300
            if mergedMode
                textY = 0.03;
            else
                textY = 0.48;
            end

            annHandle = annotation(f, 'textbox', [0.15, textY, 0.7, 0.06], ...
                'String', textString, ...
                'EdgeColor', 'none', ...
                'BackgroundColor', 'white', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 8, ...
                'Color', [0.2 0.2 0.2], ...
                'Margin', 2);
        end

        outFolderPng = fullfile(chipWellFolder, sprintf('Plot%ds', win), 'png');
        if ~exist(outFolderPng, 'dir')
            mkdir(outFolderPng);
        end

        fileNameBase = fullfile(outFolderPng, plotFileName);
        print(f, [fileNameBase '.png'], '-dpng', '-r300');

        if SvgPlot
            outFolderSvg = fullfile(chipWellFolder, sprintf('Plot%ds', win), 'svg');
            if ~exist(outFolderSvg, 'dir')
                mkdir(outFolderSvg);
            end
            fileNameBase = fullfile(outFolderSvg, plotFileName);
            print(f, [fileNameBase '.svg'], '-dsvg', '-painters');
        end

        if ~isempty(annHandle) && isvalid(annHandle)
            delete(annHandle);
        end
    end

    close(f);
end