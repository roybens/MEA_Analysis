function [bursts,burstTimeVariables] = gaussianFiringRateBurstDetector(networkData,networkAct,networkStats,AbsBP,thresholdStartStop)

                %% Tim's code for averaging and aggregating mean spiking data (IBI, Burst peaks, Spikes within Bursts, # of Bursts etc.)
                %average IBI
                meanIBI = mean(networkStats.maxAmplitudeTimeDiff);
                %average Burst peak (burst firing rate y-value)
                meanBurstPeak = mean(networkStats.maxAmplitudesValues);
                %meanBPNorm = mean(networkStatsNorm.maxAmplitudesValues);
                meanAbsBP = mean(AbsBP);
                %Number of bursts
                nBursts = length(networkStats.maxAmplitudesTimes);
                recordingTime= networkData.fileObj.dataLenTime;
                burstRate =  nBursts / recordingTime;
                % Calculate standard deviations
                stdIBI = std(networkStats.maxAmplitudeTimeDiff);
                stdBurstPeak = std(networkStats.maxAmplitudesValues);
                %stdBPNorm = std(networkStatsNorm.maxAmplitudesValues);
                stdAbsBP = std(AbsBP);
                % Calculate coefficients of variation
                covIBI = (stdIBI / meanIBI) * 100;  % CV as a percentage
                covBurstPeak = (stdBurstPeak / meanBurstPeak) * 100;  % CV as a percentage
                %covBPNorm = (stdBPNorm / meanBPNorm) * 100;  % CV as a percentage
                covAbsBP = (stdAbsBP/meanAbsBP)* 100;
                ts = ((double(networkData.fileObj.spikes.frameno)...
                    - double(networkData.fileObj.firstFrameNum))/networkData.fileObj.samplingFreq)';
                ch = networkData.fileObj.spikes.channel;
                ch = ch(ts>0);
                ts = ts(ts>0);
                %%firing rate matrix diagrams
                %if plotFRmatrix, plotFRMatrices(variables);end
                %average spikesPerBurst        
                if ~isempty(networkStats.maxAmplitudesTimes)
                    peakAmps = networkStats.maxAmplitudesValues';
                    peakTimes = networkStats.maxAmplitudesTimes;
                    % get the times of the burst start and stop edges
                    edges = double.empty(length(peakAmps),0);
                    for i = 1:length(peakAmps)
                       % take a sizeable (�6 s) chunk of the network activity curve 
                       % around each burst peak point
                       %MANDAR I THink using this there is a p
                       idx = networkAct.time>(peakTimes(i)-6) & networkAct.time<(peakTimes(i)+6);
                       t1 = networkAct.time(idx);
                       a1 = networkAct.firingRate(idx)';
                       % get the arunIDstemp = run_id_and_type(:,1);mplitude at the desired peak width
                       peakWidthAmp = (peakAmps(i)-peakAmps(i)*thresholdStartStop);
                       %the fraction of the total peak height measuredfrom the top at hich the burst start stop are %defined.
                       % get the indices of the peak edges
                       idx1 = find(a1<peakWidthAmp & t1<peakTimes(i));
                       idx2 = find(a1<peakWidthAmp & t1>peakTimes(i));
                       
                       if ~isempty(idx1)&&~isempty(idx2)       
                           tBefore = t1(idx1(end));
                           tAfter = t1(idx2(1));
                           edges(i,[1 2]) = [tBefore tAfter];
                       end
                    end
                    
                   % identify spikes that fall within the bursts
                    if isempty(edges)
                    spikesPerBurst = NaN;
                    tsWithinBurst = NaN;
                    chWithinBurst = NaN;
                    tsOutsideBurst = NaN;
                    chOutsideBurst = NaN;
                    meanSpikesPerBurst=NaN;
                    meanBurstDuration =NaN;
                    covSpikesPerBurst =NaN ;
                    covBurstDuration = NaN;
                    IBIStrings = {''};
                    BurstPeakStrings = {''};
                    AbsBurstPeakStrings = {''};
                    BurstDurationStrings = {''};
                    spikesPerBurstStrings = {''};
                    baselineFiringRate = mean(networkAct.firingRate);
                    else
                    % Initialize the spikesPerBurst array to zeros 
                    spikesPerBurst = zeros(length(edges), 1);
                    % Create an array where each element is the index of the timestamp
                    % This is used for broadcasting comparisons against edges
                    tsIndices = repmat(1:length(ts), size(edges, 1), 1);
                    % Create matrices of start and end times for easy comparison
                    startTimes = repmat(edges(:,1), 1, length(ts));
                    endTimes = repmat(edges(:,2), 1, length(ts));
                    % Logical matrices where true indicates a timestamp within the burst interval
                    inBurstWindows = (ts > startTimes) & (ts < endTimes);
                    % Sum each row to get the count of true values (spikes within each burst)
                    spikesPerBurst = sum(inBurstWindows, 2);
                    % For tsWithinBurst and chWithinBurst, you can use logical indexing
                    tsWithinBurst = ts(any(inBurstWindows, 1));
                    chWithinBurst = ch(any(inBurstWindows, 1));
                    tsOutsideBurst = ts(~any(inBurstWindows, 1));
                    chOutsideBurst = ch(~any(inBurstWindows, 1));
                    meanSpikesPerBurst = mean(spikesPerBurst);
                    stdSpikesPerBurst = std(spikesPerBurst);
                    bursts = abs(edges(:,1) - edges(:,2));
                    meanBurstDuration = mean(bursts);
                    stdBurstDuration = std(meanBurstDuration);
                    covSpikesPerBurst = (stdSpikesPerBurst/meanSpikesPerBurst)*100;
                    covBurstDuration = (stdBurstDuration/meanBurstDuration) * 100;
                    % Default value in case no bursts are detected
                    baselineFiringRate = mean(networkAct.firingRate);
                    if ~isempty(edges)
                        excludeIdx = false(size(networkAct.time));
                        for i = 1:size(edges, 1)
                            excludeIdx = excludeIdx | (networkAct.time >= edges(i,1) & networkAct.time <= edges(i,2));
                        end
                        includeIdx = ~excludeIdx;
                        if any(includeIdx)  % Avoid empty selection leading to NaN values
                            baselineFiringRate = mean(networkAct.firingRate(includeIdx));
                        end
                    end                                          
                    % Convert spikesPerBurst to a comma-separated string
                    IBIStrings = {strjoin(arrayfun(@num2str, networkStats.maxAmplitudeTimeDiff, 'UniformOutput', false), ',')};
                    BurstPeakStrings = {strjoin(arrayfun(@num2str, networkStats.maxAmplitudesValues, 'UniformOutput', false), ',')};
                    AbsBurstPeakStrings = {strjoin(arrayfun(@num2str, AbsBP, 'UniformOutput', false), ',')};
                    BurstDurationStrings = {strjoin(arrayfun(@num2str, bursts, 'UniformOutput', false), ',')};
                    spikesPerBurstStrings = {strjoin(arrayfun(@num2str, spikesPerBurst, 'UniformOutput', false), ',')};
                    end
                    else
                    spikesPerBurst = NaN;
                    tsWithinBurst = NaN;
                    chWithinBurst = NaN;
                    tsOutsideBurst = NaN;
                    chOutsideBurst = NaN;
                    meanSpikesPerBurst=NaN;
                    meanBurstDuration =NaN;
                    covSpikesPerBurst =NaN ;
                    covBurstDuration = NaN;
                    IBIStrings = {''};
                    BurstPeakStrings = {''};
                    AbsBurstPeakStrings = {''};
                    BurstDurationStrings = {''};
                    spikesPerBurstStrings = {''};
                    baselineFiringRate = mean(networkAct.firingRate);
                end
               bursts = table(...
                    meanIBI, covIBI, ...
                    meanBurstPeak, covBurstPeak, ...
                    nBursts, burstRate,...
                    meanSpikesPerBurst, covSpikesPerBurst, ...
                    meanAbsBP, covAbsBP, ...
                    meanBurstDuration, covBurstDuration, ...
                    baselineFiringRate, ...
                    IBIStrings, BurstPeakStrings, AbsBurstPeakStrings, ...
                    BurstDurationStrings, spikesPerBurstStrings, ...
                    'VariableNames', { ...
                        'mean_IBI', 'cov_IBI', ...
                        'mean_Burst_Peak', 'cov_Burst_Peak', ...
                        'Number_Bursts','burstRate' ...
                        'mean_Spike_per_Burst', 'cov_Spike_per_Burst', ...
                        'mean_Burst_Peak_Abs', 'cov_Burst_Peak_Abs', ...
                        'mean_BurstDuration', 'cov_BurstDuration', ...
                        'BaselineFiringRate', ...
                        'IBI_Distribution', 'Burst_Peak_Distribution', 'Abs_Burst_Peak_Distribution', ...
                        'Burst_Duration_Distribution', 'SpikesPerBurst_Distribution' ...
                    });
               burstTimeVariables =struct();
               burstTimeVariables.chWithinBurst = chWithinBurst;
               burstTimeVariables.tsWithinBurst = tsWithinBurst;
               burstTimeVariables.tsOutsideBurst = tsOutsideBurst;
               burstTimeVariables.chOutsideBurst = chOutsideBurst;

    end