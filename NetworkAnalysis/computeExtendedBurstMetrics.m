 function extMetrics = computeExtendedBurstMetrics(burstMetricsVariables, relativeSpikeTimes, ts)

    if isnan(burstMetricsVariables.tsWithinBurst)
        burstISIs = {};
        burstISIsCombined = [];
        meanBurstISI = NaN;
        IQRABurstISI = NaN;
        stdBurstISI = NaN;
        covBurstISI = NaN;
        instAPFreqWithinBurst = NaN;
        iqrAPFreqWithinBurst = NaN;
        rangeAPFreq = NaN;
        binWidth = NaN;
        numBins = NaN;
        binedgesAP = NaN;
        binAPcountsWithinBurst = NaN;
        binedgesAPWithinBurst = NaN;

        meanBurstISIOutside = NaN;
        covBurstISIOutside = NaN;
        binAPcountsOutsideBurst = NaN;
        binedgesAPOutsideBurst = NaN;
    else
        burstISIs = accumarray(burstMetricsVariables.chWithinBurst, burstMetricsVariables.tsWithinBurst, [], @(x){diff(x)});
        burstISIsCombined = vertcat(burstISIs{:});
        burstISIsCombined = burstISIsCombined(burstISIsCombined < 1);
        meanBurstISI = mean(burstISIsCombined);
        IQRABurstISI = iqr(burstISIsCombined);
        stdBurstISI = std(burstISIsCombined);
        covBurstISI = (stdBurstISI / meanBurstISI) * 100;

        instAPFreqWithinBurst = 1 ./ burstISIsCombined;
        iqrAPFreqWithinBurst = iqr(instAPFreqWithinBurst);
        rangeAPFreq = max(instAPFreqWithinBurst) - min(instAPFreqWithinBurst);
        binWidth = 2 * iqrAPFreqWithinBurst / (length(instAPFreqWithinBurst)^(1/3));
        numBins = ceil(rangeAPFreq / binWidth);
        binedgesAP = logspace(log10(min(instAPFreqWithinBurst)), log10(500), numBins);
        [binAPcountsWithinBurst, binedgesAPWithinBurst] = histcounts(instAPFreqWithinBurst, binedgesAP);

        burstISIsOutside = accumarray(burstMetricsVariables.chOutsideBurst, burstMetricsVariables.tsOutsideBurst, [], @(x) {diff(sort(x))});
        burstISIsCombinedOutside = vertcat(burstISIsOutside{:});
        burstISIsCombinedOutside = burstISIsCombinedOutside(burstISIsCombinedOutside < 1);
        meanBurstISIOutside = mean(burstISIsCombinedOutside);
        stdBurstISIOutside = std(burstISIsCombinedOutside);
        covBurstISIOutside = (stdBurstISIOutside / meanBurstISIOutside) * 100;

        instAPFreqOutsideBurst = 1 ./ burstISIsCombinedOutside;
        iqrAPFreqOutsideBurst = iqr(instAPFreqOutsideBurst);
        rangeAPFreqOutside = max(instAPFreqOutsideBurst) - min(instAPFreqOutsideBurst);
        binWidthOutside = 2 * iqrAPFreqOutsideBurst / (length(instAPFreqOutsideBurst)^(1/3));
        numBinsOutside = ceil(rangeAPFreqOutside / binWidthOutside);
        binedgesAPOutside = logspace(log10(min(instAPFreqOutsideBurst)), log10(500), numBinsOutside);
        [binAPcountsOutsideBurst, binedgesAPOutsideBurst] = histcounts(instAPFreqOutsideBurst, binedgesAPOutside);
    end

    ISIs = accumarray(relativeSpikeTimes.channel, relativeSpikeTimes.time, [], @(x){diff(x)});
    ISIsCombined = vertcat(ISIs{:});
    ISIsCombined = ISIsCombined(ISIsCombined > 0);
    meanISI = mean(ISIsCombined);
    stdISI = std(ISIsCombined);
    covISI = (stdISI / meanISI) * 100;

    instFreq = 1 ./ ISIsCombined;
    iqrFreq = iqr(instFreq);
    rangeFreq = max(instFreq) - min(instFreq);
    binWidth = 2 * iqrFreq / (length(instFreq)^(1/3));
    numBins = ceil(rangeFreq / binWidth);
    minFreq = min(instFreq(instFreq > 0));
    binedgesFreq = logspace(log10(minFreq), log10(max(instFreq)), numBins);
    [binFreqCounts, binedgesFreq] = histcounts(instFreq, binedgesFreq);

    binedgesfano = min(ts):0.1:max(ts);
    spikeCountsFano = histcounts(ts, binedgesfano);
    meanSpikeCounts = mean(spikeCountsFano);
    varSpikeCounts = var(spikeCountsFano);
    fanoFactor = varSpikeCounts / meanSpikeCounts;

   
    extMetrics = table( ...
    {binFreqCounts},      {binedgesFreq}, {binAPcountsWithinBurst}, {binedgesAPWithinBurst}, ...
    {binAPcountsOutsideBurst}, {binedgesAPOutsideBurst},meanISI,covISI,fanoFactor,...
    'VariableNames', { 'networkAPFreqBins', 'networkAPFreqEdges', ...
      'burstAPFreqBins',   'burstAPFreqEdges', 'nonburstAPFreqBins','nonburstAPFreqEdges', 'meanISI','covISI','fanoFactor' } );

end
