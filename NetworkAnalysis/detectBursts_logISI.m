function [bursts] = detectBursts_logISI(spike_times, min_spikes_in_burst, min_burst_duration_s, isi_threshold_method_params)
% detectBursts_logISI: Detects bursts in a spike train using the logISI histogram method.
% (Vectorized burst identification step)
%
% Inputs:
%   spike_times:            Vector of spike times (in seconds). Must be sorted.
%   min_spikes_in_burst:    Minimum number of spikes to constitute a burst.
%   min_burst_duration_s:   Minimum duration (in seconds) for a burst to be considered valid.
%   isi_threshold_method_params: (Optional) Structure with parameters for ISI threshold detection:
%       .log_isi_bins:       Number of bins for the log(ISI) histogram (e.g., 50 or 100). Default: 50.
%       .smooth_span_hist:   Span for moving average smoothing of log(ISI) histogram (e.g., 5). Default: 5.
%       .max_intra_burst_isi_s: Initial estimate/upper bound for intra-burst ISI (e.g., 0.1 for 100ms). Default: 0.2.
%
% Outputs:
%   bursts: A structure array, where each element represents a detected burst and contains:
%       .start_time:    Start time of the burst (s).
%       .end_time:      End time of the burst (s).
%       .duration_s:    Duration of the burst (s).
%       .num_spikes:    Number of spikes in the burst.
%       .spike_times:   Array of spike times within the burst (s).
%       .peak_isi_s:    The ISI value corresponding to the main peak in the log(ISI) histogram (intra-burst ISIs)
%       .threshold_isi_s: The adaptively determined ISI threshold used to define bursts.

    bursts = struct('start_time', {}, 'end_time', {}, 'duration_s', {}, ...
                    'num_spikes', {}, 'spike_times', {}, ...
                    'peak_isi_s', {}, 'threshold_isi_s', {});

    if nargin < 4 || isempty(isi_threshold_method_params)
        isi_threshold_method_params = struct();
    end
    if ~isfield(isi_threshold_method_params, 'log_isi_bins'), isi_threshold_method_params.log_isi_bins = 50; end
    if ~isfield(isi_threshold_method_params, 'smooth_span_hist'), isi_threshold_method_params.smooth_span_hist = 5; end
    if ~isfield(isi_threshold_method_params, 'max_intra_burst_isi_s'), isi_threshold_method_params.max_intra_burst_isi_s = 0.1; end  %Chiappalone et al. 2005

    if length(spike_times) < min_spikes_in_burst
        return; 
    end

    isis = diff(spike_times);
    if isempty(isis) || length(isis) < min_spikes_in_burst -1
        if length(spike_times) >= min_spikes_in_burst
            duration = spike_times(end) - spike_times(1);
            if duration >= min_burst_duration_s && (isempty(isis) || all(isis <= isi_threshold_method_params.max_intra_burst_isi_s))
                 bursts(1).start_time = spike_times(1);
                 bursts(1).end_time = spike_times(end);
                 bursts(1).duration_s = duration;
                 bursts(1).num_spikes = length(spike_times);
                 bursts(1).spike_times = spike_times;
                 bursts(1).peak_isi_s = ifelse(isempty(isis), NaN, mean(isis)); 
                 bursts(1).threshold_isi_s = isi_threshold_method_params.max_intra_burst_isi_s;
            end
        end
        return;
    end
    
    positive_isis = isis(isis > 1e-9); 
    if isempty(positive_isis) || length(positive_isis) < min_spikes_in_burst -1
        return;
    end
    log_isis = log10(positive_isis); 
    
    min_log_isi = min(log_isis);
    max_log_isi = max(log_isis);

    if isinf(min_log_isi) || isinf(max_log_isi) || (max_log_isi - min_log_isi < 1e-6) 
        if length(spike_times) >= min_spikes_in_burst && all(positive_isis <= isi_threshold_method_params.max_intra_burst_isi_s)
            duration = spike_times(end) - spike_times(1);
            if duration >= min_burst_duration_s
                bursts(1).start_time = spike_times(1);
                bursts(1).end_time = spike_times(end);
                bursts(1).duration_s = duration;
                bursts(1).num_spikes = length(spike_times);
                bursts(1).spike_times = spike_times;
                bursts(1).peak_isi_s = mean(positive_isis);
                bursts(1).threshold_isi_s = isi_threshold_method_params.max_intra_burst_isi_s;
            end
        end
        return;
    end

    hist_edges = linspace(min_log_isi, max_log_isi, isi_threshold_method_params.log_isi_bins + 1);
    [counts, edges_hist] = histcounts(log_isis, hist_edges);
    centers = (edges_hist(1:end-1) + edges_hist(2:end)) / 2;
    smoothed_counts = smoothdata(counts, 'movmean', isi_threshold_method_params.smooth_span_hist);
    
    final_peak_isi = NaN;
    final_isi_threshold = isi_threshold_method_params.max_intra_burst_isi_s; %fallback threshold
    [peak_values_hist, peak_locs_indices_hist] = findpeaks(smoothed_counts);
    if ~isempty(peak_locs_indices_hist)
        candidate_intra_burst_peaks_log_centers = centers(peak_locs_indices_hist);
        valid_candidate_indices = find(candidate_intra_burst_peaks_log_centers <= log10(isi_threshold_method_params.max_intra_burst_isi_s));
        
        if ~isempty(valid_candidate_indices)
            peak_counts_of_candidates = peak_values_hist(ismember(peak_locs_indices_hist, peak_locs_indices_hist(valid_candidate_indices)));
            [~, max_val_idx_within_candidates] = max(peak_counts_of_candidates);
            
            main_peak_hist_idx = peak_locs_indices_hist(valid_candidate_indices(max_val_idx_within_candidates));
            final_peak_isi = 10^centers(main_peak_hist_idx);

            temp_counts_for_valley = smoothed_counts;
            temp_counts_for_valley(1:main_peak_hist_idx) = Inf; 
            [~, valley_locs_indices_hist] = findpeaks(-temp_counts_for_valley); 

            if ~isempty(valley_locs_indices_hist)
                first_valley_hist_idx = valley_locs_indices_hist(1);
                threshold_from_valley = 10^centers(first_valley_hist_idx);
                if threshold_from_valley < final_isi_threshold * 2.5 % Heuristic check
                    final_isi_threshold = threshold_from_valley;
                end
            else
                 idx_half_max_right = find(smoothed_counts(main_peak_hist_idx:end) < smoothed_counts(main_peak_hist_idx)/2, 1, 'first');
                 if ~isempty(idx_half_max_right) && (main_peak_hist_idx + idx_half_max_right - 1) <= length(centers)
                     final_isi_threshold = 10^(centers(main_peak_hist_idx + idx_half_max_right - 1));
                 end
            end
        end
    end

    % --- Step 3: Identify Bursts (Vectorized Approach) ---
    is_intra_burst_isi = (isis <= final_isi_threshold);

    % Find start and end indices of consecutive true values in is_intra_burst_isi
    % These indices refer to the 'isis' array.
    diff_is_intra = diff([0; is_intra_burst_isi; 0]); % Pad to detect bursts at start/end
    
    % Start indices of potential bursts (where sequence of short ISIs begins)
    % These are indices IN THE 'isis' ARRAY.
    start_isi_indices = find(diff_is_intra == 1);
    % End indices of potential bursts (where sequence of short ISIs ends)
    % These are also indices IN THE 'isis' ARRAY.
    end_isi_indices = find(diff_is_intra == -1) - 1;

    burst_counter = 0;
    if ~isempty(start_isi_indices) && ~isempty(end_isi_indices)
        % Ensure consistent lengths if file ends with a burst sequence
        if length(start_isi_indices) > length(end_isi_indices)
             end_isi_indices(end+1) = length(is_intra_burst_isi); % Assume burst goes to end
        end
        
        for i = 1:length(start_isi_indices)
            current_start_isi_idx = start_isi_indices(i);
            current_end_isi_idx = end_isi_indices(i);

            % Number of ISIs in this burst candidate
            num_isis_in_candidate = current_end_isi_idx - current_start_isi_idx + 1;
            % Number of spikes in this burst candidate
            num_spikes_current_burst = num_isis_in_candidate + 1;

            if num_spikes_current_burst >= min_spikes_in_burst
                % Spike times:
                % The first spike of the burst is spike_times(current_start_isi_idx)
                % The last spike of the burst is spike_times(current_end_isi_idx + 1)
                burst_start_time = spike_times(current_start_isi_idx);
                burst_end_time = spike_times(current_end_isi_idx + 1);
                duration = burst_end_time - burst_start_time;

                if duration >= min_burst_duration_s
                    burst_counter = burst_counter + 1;
                    bursts(burst_counter).start_time = burst_start_time;
                    bursts(burst_counter).end_time = burst_end_time;
                    bursts(burst_counter).duration_s = duration;
                    bursts(burst_counter).num_spikes = num_spikes_current_burst;
                    bursts(burst_counter).spike_times = spike_times(current_start_isi_idx : current_end_isi_idx + 1);
                    bursts(burst_counter).peak_isi_s = final_peak_isi;
                    bursts(burst_counter).threshold_isi_s = final_isi_threshold;
                end
            end
        end
    end
end

% Helper function (if not already defined elsewhere in your path)
function out = ifelse(condition, true_val, false_val)
    if condition
        out = true_val;
    else
        out = false_val;
    end
end