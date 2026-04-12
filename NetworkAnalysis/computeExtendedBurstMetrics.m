function extMetrics = computeExtendedBurstMetrics(burstMetricsVariables, relativeSpikeTimes, ts)

    % -------------------------------
    % SAFETY HELPERS
    % -------------------------------
    safeCov = @(x) (mean(x) > 0) * (std(x)/mean(x)*100) + (mean(x) == 0) * NaN;

    % -------------------------------
    % HANDLE EMPTY BURST CASE
    % -------------------------------
    noBurst = isempty(burstMetricsVariables.tsWithinBurst) || ...
              (isscalar(burstMetricsVariables.tsWithinBurst) && isnan(burstMetricsVariables.tsWithinBurst));

    % -------------------------------
    % PER-CHANNEL ISI (DO NOT COLLAPSE)
    % -------------------------------
    ISIs = accumarray(relativeSpikeTimes.channel, relativeSpikeTimes.time, [], @(x){diff(x)});
    ISIs = ISIs(~cellfun(@isempty, ISIs));

    % Per-electrode stats
    meanISI_e = cellfun(@mean, ISIs);
    medianISI_e = cellfun(@median, ISIs);
    stdISI_e = cellfun(@std, ISIs);
    iqrISI_e = cellfun(@iqr, ISIs);
    covISI_e = arrayfun(@(i) safeCov(ISIs{i}), 1:numel(ISIs));

    % Aggregate electrode stats (ROBUST)
    meanISI = mean(meanISI_e);
    covISI = mean(covISI_e, 'omitnan');
    medianISI = median(medianISI_e);
    iqrISI = median(iqrISI_e);

    % -------------------------------
    % INSTANTANEOUS FREQUENCY (NETWORK)
    % -------------------------------
    ISIsCombined = vertcat(ISIs{:});
    ISIsCombined = ISIsCombined(ISIsCombined > 0);

    instFreq = 1 ./ ISIsCombined;
    instFreq = instFreq(isfinite(instFreq) & instFreq > 0);

    % Fixed log bins (robust across datasets)
    freqEdges = logspace(log10(0.1), log10(500), 50);
    [binFreqCounts, ~] = histcounts(instFreq, freqEdges);

    % -------------------------------
    % BURST vs NON-BURST
    % -------------------------------
    if noBurst
        meanBurstISI = NaN; covBurstISI = NaN;
        meanBurstISIOutside = NaN; covBurstISIOutside = NaN;

        binAPcountsWithinBurst = NaN; binedgesAPWithinBurst = NaN;
        binAPcountsOutsideBurst = NaN; binedgesAPOutsideBurst = NaN;

        burstFraction_e = NaN;
        participation_e = NaN;

    else
        % --- WITHIN BURST ---
        burstISIs = accumarray(burstMetricsVariables.chWithinBurst, ...
                              burstMetricsVariables.tsWithinBurst, [], @(x){diff(x)});
        burstISIs = burstISIs(~cellfun(@isempty, burstISIs));

        burstISIsCombined = vertcat(burstISIs{:});
        burstISIsCombined = burstISIsCombined(burstISIsCombined > 0);

        meanBurstISI = mean(burstISIsCombined);
        covBurstISI = safeCov(burstISIsCombined);

        instAPFreqWithinBurst = 1 ./ burstISIsCombined;
        instAPFreqWithinBurst = instAPFreqWithinBurst(instAPFreqWithinBurst > 0);

        edges = logspace(log10(0.1), log10(500), 50);
        [binAPcountsWithinBurst, binedgesAPWithinBurst] = histcounts(instAPFreqWithinBurst, edges);

        % --- OUTSIDE BURST ---
        burstISIsOutside = accumarray(burstMetricsVariables.chOutsideBurst, ...
                                     burstMetricsVariables.tsOutsideBurst, [], @(x){diff(sort(x))});
        burstISIsOutside = burstISIsOutside(~cellfun(@isempty, burstISIsOutside));

        burstISIsCombinedOutside = vertcat(burstISIsOutside{:});
        burstISIsCombinedOutside = burstISIsCombinedOutside(burstISIsCombinedOutside > 0);

        meanBurstISIOutside = mean(burstISIsCombinedOutside);
        covBurstISIOutside = safeCov(burstISIsCombinedOutside);

        instAPFreqOutsideBurst = 1 ./ burstISIsCombinedOutside;
        instAPFreqOutsideBurst = instAPFreqOutsideBurst(instAPFreqOutsideBurst > 0);

        [binAPcountsOutsideBurst, binedgesAPOutsideBurst] = histcounts(instAPFreqOutsideBurst, edges);

        % -------------------------------
        % ELECTRODE-LEVEL METRICS (NEW)
        % -------------------------------
        % Spike counts per electrode
        spikeCounts = cellfun(@length, ISIs);

        % Burst spike counts per electrode
        burstSpikeCounts = accumarray(burstMetricsVariables.chWithinBurst, 1);
        totalSpikeCounts = accumarray(relativeSpikeTimes.channel, 1);

        % Align sizes
        maxCh = max([length(totalSpikeCounts), length(burstSpikeCounts)]);
        totalSpikeCounts(end+1:maxCh) = 0;
        burstSpikeCounts(end+1:maxCh) = 0;

        % Burst fraction per electrode
        burstFraction_e = burstSpikeCounts ./ max(totalSpikeCounts,1);

        % Participation (binary: did electrode participate at least once)
        participation_e = burstSpikeCounts > 0;

    end

    % -------------------------------
    % ELECTRODE AGGREGATES (ROBUST)
    % -------------------------------
    meanBurstFraction = mean(burstFraction_e, 'omitnan');
    medianBurstFraction = median(burstFraction_e, 'omitnan');

    participationRate = mean(participation_e, 'omitnan'); % fraction of active electrodes

    % -------------------------------
    % FANO FACTOR (GLOBAL)
    % -------------------------------
    binEdgesFano = min(ts):0.1:max(ts);
    spikeCountsFano = histcounts(ts, binEdgesFano);
    meanSpikeCounts = mean(spikeCountsFano);
    varSpikeCounts = var(spikeCountsFano);

    if meanSpikeCounts > 0
        fanoFactor = varSpikeCounts / meanSpikeCounts;
    else
        fanoFactor = NaN;
    end

    % -------------------------------
    % FINAL TABLE (API PRESERVED + EXTENDED)
    % -------------------------------
    extMetrics = table( ...
        {binFreqCounts}, {freqEdges}, ...
        {binAPcountsWithinBurst}, {binedgesAPWithinBurst}, ...
        {binAPcountsOutsideBurst}, {binedgesAPOutsideBurst}, ...
        meanISI, covISI, fanoFactor, ...
        medianISI, iqrISI, ...                      % NEW robust ISI
        meanBurstFraction, medianBurstFraction, ... % NEW electrode burstiness
        participationRate, ...                      % NEW recruitment
        'VariableNames', { ...
        'networkAPFreqBins', 'networkAPFreqEdges', ...
        'burstAPFreqBins', 'burstAPFreqEdges', ...
        'nonburstAPFreqBins', 'nonburstAPFreqEdges', ...
        'meanISI','covISI','fanoFactor', ...
        'medianISI','iqrISI', ...
        'meanBurstFraction','medianBurstFraction', ...
        'participationRate' ...
        });

end
