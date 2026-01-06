import numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt

def logISI_burst_detection(spike_times, 
                           max_isi_cutoff=0.1,  # 100 ms default cutoff
                           void_threshold=0.7,  # Void parameter threshold
                           min_spikes_in_burst=3,
                           hist_bins=100,
                           smooth_sigma=2.0):
    """
    Detect bursts using the LogISI method (Pasquale et al. 2010).
    
    This method analyzes the histogram of log-transformed interspike intervals (ISIs)
    to automatically determine the optimal ISI threshold for burst detection.
    
    Algorithm:
    1. Calculate all ISIs and transform to log scale
    2. Create histogram of log(ISI) values
    3. Find peaks in the histogram
    4. Identify the "intraburst peak" (largest peak with ISI ≤ max_isi_cutoff)
    5. Calculate "void parameter" between peaks to assess separation
    6. Set ISI cutoff where void parameter exceeds threshold
    7. Detect bursts as sequences of ≥3 spikes with ISI < cutoff
    
    Parameters:
    -----------
    spike_times : array-like
        Array of spike times in seconds
    max_isi_cutoff : float
        Maximum allowed ISI cutoff value in seconds (default: 0.1 = 100 ms)
    void_threshold : float
        Threshold for void parameter (default: 0.7, range: 0-1)
        Higher values = more stringent separation requirement
    min_spikes_in_burst : int
        Minimum number of spikes required to define a burst (default: 3)
    hist_bins : int
        Number of bins for log(ISI) histogram (default: 100)
    smooth_sigma : float
        Sigma for Gaussian smoothing of histogram (default: 2.0)
    
    Returns:
    --------
    bursts : list of tuples
        List of (burst_start_time, burst_end_time, num_spikes) for each detected burst
    isi_cutoff : float
        The ISI cutoff value determined by the algorithm (in seconds)
    burst_info : dict
        Dictionary containing detailed burst statistics
    """
    
    spike_times = np.array(spike_times)
    if len(spike_times) < min_spikes_in_burst:
        return [], np.nan, {}
    
    # Step 1: Calculate interspike intervals
    isis = np.diff(spike_times)
    
    if len(isis) == 0:
        return [], np.nan, {}
    
    # Step 2: Log-transform ISIs (add small value to avoid log(0))
    log_isis = np.log10(isis + 1e-10)
    
    # Step 3: Create histogram of log(ISI) values
    hist_counts, hist_edges = np.histogram(log_isis, bins=hist_bins)
    hist_centers = (hist_edges[:-1] + hist_edges[1:]) / 2
    
    # Smooth the histogram for better peak detection
    hist_counts_smooth = gaussian_filter1d(hist_counts.astype(float), sigma=smooth_sigma)
    
    # Step 4: Find peaks in the log(ISI) histogram
    peaks, properties = find_peaks(hist_counts_smooth, prominence=np.max(hist_counts_smooth) * 0.1)
    
    if len(peaks) == 0:
        # No clear peaks found - use default cutoff
        isi_cutoff = max_isi_cutoff
    else:
        # Find the "intraburst peak" - largest peak corresponding to ISI ≤ max_isi_cutoff
        log_max_cutoff = np.log10(max_isi_cutoff)
        intraburst_candidates = peaks[hist_centers[peaks] <= log_max_cutoff]
        
        if len(intraburst_candidates) == 0:
            # No peak found below max_isi_cutoff - use default
            isi_cutoff = max_isi_cutoff
        else:
            # Select the largest peak (highest histogram count) as intraburst peak
            intraburst_peak_idx = intraburst_candidates[np.argmax(hist_counts_smooth[intraburst_candidates])]
            
            # Step 5: Calculate void parameter between intraburst peak and subsequent peaks
            # The void parameter quantifies how well-separated two peaks are
            subsequent_peaks = peaks[peaks > intraburst_peak_idx]
            
            if len(subsequent_peaks) == 0:
                # No subsequent peaks - use default cutoff
                isi_cutoff = max_isi_cutoff
            else:
                # Find minima between intraburst peak and each subsequent peak
                void_params = []
                cutoff_positions = []
                
                for next_peak in subsequent_peaks:
                    # Find minimum between the two peaks
                    search_region = hist_counts_smooth[intraburst_peak_idx:next_peak+1]
                    if len(search_region) > 0:
                        min_idx = np.argmin(search_region) + intraburst_peak_idx
                        
                        # Calculate void parameter
                        # Void = (height_peak1 + height_peak2 - 2*height_min) / (height_peak1 + height_peak2)
                        h1 = hist_counts_smooth[intraburst_peak_idx]
                        h2 = hist_counts_smooth[next_peak]
                        h_min = hist_counts_smooth[min_idx]
                        
                        if (h1 + h2) > 0:
                            void_param = (h1 + h2 - 2*h_min) / (h1 + h2)
                            void_params.append(void_param)
                            cutoff_positions.append(min_idx)
                
                # Step 6: Set ISI cutoff at first minimum where void parameter exceeds threshold
                if len(void_params) > 0:
                    valid_cutoffs = [cutoff_positions[i] for i, vp in enumerate(void_params) 
                                   if vp >= void_threshold]
                    
                    if len(valid_cutoffs) > 0:
                        cutoff_idx = valid_cutoffs[0]
                        isi_cutoff = 10 ** hist_centers[cutoff_idx]
                        
                        # Ensure cutoff doesn't exceed maximum
                        if isi_cutoff > max_isi_cutoff:
                            isi_cutoff = max_isi_cutoff
                    else:
                        # No void parameter exceeded threshold - use default
                        isi_cutoff = max_isi_cutoff
                else:
                    isi_cutoff = max_isi_cutoff
    
    # Step 7: Detect bursts using the determined ISI cutoff
    bursts = []
    current_burst_spikes = [0]  # Start with first spike index
    
    for i in range(len(isis)):
        if isis[i] <= isi_cutoff:
            # Continue current burst
            current_burst_spikes.append(i + 1)
        else:
            # End current burst and check if valid
            if len(current_burst_spikes) >= min_spikes_in_burst:
                burst_start = spike_times[current_burst_spikes[0]]
                burst_end = spike_times[current_burst_spikes[-1]]
                num_spikes = len(current_burst_spikes)
                bursts.append((burst_start, burst_end, num_spikes))
            
            # Start new burst
            current_burst_spikes = [i + 1]
    
    # Check last burst
    if len(current_burst_spikes) >= min_spikes_in_burst:
        burst_start = spike_times[current_burst_spikes[0]]
        burst_end = spike_times[current_burst_spikes[-1]]
        num_spikes = len(current_burst_spikes)
        bursts.append((burst_start, burst_end, num_spikes))
    
    # Calculate burst statistics
    if len(bursts) > 0:
        burst_durations = [b[1] - b[0] for b in bursts]
        spikes_per_burst = [b[2] for b in bursts]
        
        # Calculate inter-burst intervals (offset to onset)
        if len(bursts) > 1:
            ibis = [bursts[i+1][0] - bursts[i][1] for i in range(len(bursts)-1)]
            mean_ibi = np.mean(ibis)
            cov_ibi = np.std(ibis) / mean_ibi if mean_ibi > 0 else np.nan
        else:
            mean_ibi = np.nan
            cov_ibi = np.nan
        
        # Count spikes in bursts
        total_burst_spikes = sum(spikes_per_burst)
        burst_percentage = (total_burst_spikes / len(spike_times)) * 100
        
        burst_info = {
            'num_bursts': len(bursts),
            'mean_burst_duration': np.mean(burst_durations),
            'std_burst_duration': np.std(burst_durations),
            'cov_burst_duration': np.std(burst_durations) / np.mean(burst_durations) if np.mean(burst_durations) > 0 else np.nan,
            'mean_spikes_per_burst': np.mean(spikes_per_burst),
            'std_spikes_per_burst': np.std(spikes_per_burst),
            'mean_ibi': mean_ibi,
            'cov_ibi': cov_ibi,
            'burst_percentage': burst_percentage,
            'isi_cutoff_used': isi_cutoff
        }
    else:
        burst_info = {
            'num_bursts': 0,
            'mean_burst_duration': np.nan,
            'std_burst_duration': np.nan,
            'cov_burst_duration': np.nan,
            'mean_spikes_per_burst': np.nan,
            'std_spikes_per_burst': np.nan,
            'mean_ibi': np.nan,
            'cov_ibi': np.nan,
            'burst_percentage': 0.0,
            'isi_cutoff_used': isi_cutoff
        }
    
    return bursts, isi_cutoff, burst_info


def ISI_N_burst_detection(spike_times,
                         N=3,  # Minimum number of consecutive spikes
                         ISI_threshold=0.1,  # Maximum ISI in seconds
                         min_spikes_in_burst=3):
    """
    Detect bursts using the ISI_N method (simple threshold-based).
    
    A burst is defined as N or more consecutive spikes with ISI < ISI_threshold.
    This is a simpler, more deterministic alternative to logISI.
    
    Parameters:
    -----------
    spike_times : array-like
        Array of spike times in seconds
    N : int
        Minimum number of consecutive spikes required (default: 3)
    ISI_threshold : float
        Maximum ISI for spikes to be considered in a burst (seconds, default: 0.1)
    min_spikes_in_burst : int
        Minimum spikes in a burst (usually same as N, default: 3)
    
    Returns:
    --------
    bursts : list of tuples
        List of (burst_start_time, burst_end_time, num_spikes) for each detected burst
    burst_info : dict
        Dictionary containing detailed burst statistics
    """
    
    spike_times = np.array(spike_times)
    if len(spike_times) < min_spikes_in_burst:
        return [], {}
    
    # Calculate ISIs
    isis = np.diff(spike_times)
    
    if len(isis) == 0:
        return [], {}
    
    # Detect bursts
    bursts = []
    current_burst_spikes = [0]  # Start with first spike
    
    for i in range(len(isis)):
        if isis[i] <= ISI_threshold:
            # Continue current burst
            current_burst_spikes.append(i + 1)
        else:
            # End current burst and check if valid
            if len(current_burst_spikes) >= min_spikes_in_burst:
                burst_start = spike_times[current_burst_spikes[0]]
                burst_end = spike_times[current_burst_spikes[-1]]
                num_spikes = len(current_burst_spikes)
                bursts.append((burst_start, burst_end, num_spikes))
            
            # Start new burst
            current_burst_spikes = [i + 1]
    
    # Check last burst
    if len(current_burst_spikes) >= min_spikes_in_burst:
        burst_start = spike_times[current_burst_spikes[0]]
        burst_end = spike_times[current_burst_spikes[-1]]
        num_spikes = len(current_burst_spikes)
        bursts.append((burst_start, burst_end, num_spikes))
    
    # Calculate statistics
    if len(bursts) > 0:
        burst_durations = [b[1] - b[0] for b in bursts]
        spikes_per_burst = [b[2] for b in bursts]
        
        if len(bursts) > 1:
            ibis = [bursts[i+1][0] - bursts[i][1] for i in range(len(bursts)-1)]
            mean_ibi = np.mean(ibis)
            cov_ibi = np.std(ibis) / mean_ibi if mean_ibi > 0 else np.nan
        else:
            mean_ibi = np.nan
            cov_ibi = np.nan
        
        total_burst_spikes = sum(spikes_per_burst)
        burst_percentage = (total_burst_spikes / len(spike_times)) * 100
        
        burst_info = {
            'num_bursts': len(bursts),
            'mean_burst_duration': np.mean(burst_durations),
            'std_burst_duration': np.std(burst_durations),
            'cov_burst_duration': np.std(burst_durations) / np.mean(burst_durations) if np.mean(burst_durations) > 0 else np.nan,
            'mean_spikes_per_burst': np.mean(spikes_per_burst),
            'std_spikes_per_burst': np.std(spikes_per_burst),
            'mean_ibi': mean_ibi,
            'cov_ibi': cov_ibi,
            'burst_percentage': burst_percentage
        }
    else:
        burst_info = {
            'num_bursts': 0,
            'mean_burst_duration': np.nan,
            'std_burst_duration': np.nan,
            'cov_burst_duration': np.nan,
            'mean_spikes_per_burst': np.nan,
            'std_spikes_per_burst': np.nan,
            'mean_ibi': np.nan,
            'cov_ibi': np.nan,
            'burst_percentage': 0.0
        }
    
    return bursts, burst_info



import numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt

def plot_network_bursts_logISI(ax, SpikeTimes,
                                max_isi_cutoff=0.1,
                                void_threshold=0.7,
                                min_spikes_in_burst=3,
                                min_active_electrodes=0.25,  # Fraction of electrodes that must burst
                                time_window=0.05,  # 50 ms window for synchrony detection
                                min_active_threshold=0.1):  # Hz for active electrode
    """
    Detect NETWORK bursts using LogISI method applied to individual electrodes,
    then aggregating across the network.
    
    Network burst = when a sufficient fraction of electrodes burst within a time window.
    
    Algorithm:
    1. Apply LogISI burst detection to each electrode independently
    2. Slide a time window through the recording
    3. Count how many electrodes are bursting at each time point
    4. Network burst = when ≥ threshold electrodes burst simultaneously
    
    Parameters:
    -----------
    ax : matplotlib axis object
        Axis to plot on
    SpikeTimes : dict
        Dictionary where keys are unit/electrode IDs and values are arrays of spike times (in seconds)
    max_isi_cutoff : float
        Maximum allowed ISI cutoff for LogISI (default: 0.1 = 100 ms)
    void_threshold : float
        Void parameter threshold for LogISI (default: 0.7)
    min_spikes_in_burst : int
        Minimum spikes per single-electrode burst (default: 3)
    min_active_electrodes : float or int
        If float < 1: fraction of active electrodes that must participate
        If int >= 1: minimum number of electrodes that must participate
        (default: 0.25 = 25% of electrodes)
    time_window : float
        Time window for detecting synchronous bursting (default: 0.05 = 50 ms)
    min_active_threshold : float
        Minimum firing rate (Hz) for an electrode to be considered active
    
    Returns:
    --------
    network_data : dict
        Dictionary containing network burst metrics including:
        - num_network_bursts: Number of detected network bursts
        - network_burst_rate: Network bursts per second
        - mean_network_burst_duration: Mean network burst duration (seconds)
        - cov_network_burst_duration: Coefficient of variation of duration
        - mean_ibi: Mean inter-burst interval (seconds)
        - cov_ibi: Coefficient of variation of IBI
        - mean_spikes_per_network_burst: Average spikes per network burst
        - mean_electrodes_per_network_burst: Average electrode participation
        - network_burst_percentage: Percentage of spikes in network bursts
        - active_electrodes: Number of active electrodes
        - bursting_electrodes: Number of electrodes with detected bursts
    """
    
    # Step 1: Detect bursts on each electrode using LogISI
    electrode_bursts = {}
    active_electrodes = []
    all_spike_times = []
    
    for electrode_id, spike_times in SpikeTimes.items():
        if len(spike_times) == 0:
            continue
            
        all_spike_times.extend(spike_times)
        
        # Check if electrode is active
        recording_duration = max(spike_times) - min(spike_times) if len(spike_times) > 1 else 1.0
        firing_rate = len(spike_times) / recording_duration
        
        if firing_rate < min_active_threshold:
            continue
        
        active_electrodes.append(electrode_id)
        
        # Detect bursts on this electrode
        bursts, isi_cutoff, burst_info = logISI_burst_detection(
            spike_times, max_isi_cutoff, void_threshold, min_spikes_in_burst
        )
        
        if len(bursts) > 0:
            electrode_bursts[electrode_id] = {
                'bursts': bursts,
                'isi_cutoff': isi_cutoff,
                'info': burst_info
            }
    
    if len(active_electrodes) == 0:
        print("Warning: No active electrodes detected!")
        return None
    
    if len(electrode_bursts) == 0:
        print("Warning: No bursts detected on any electrode!")
        return None
    
    # Determine minimum number of electrodes needed for network burst
    if isinstance(min_active_electrodes, float) and min_active_electrodes < 1:
        min_electrodes_threshold = int(np.ceil(min_active_electrodes * len(active_electrodes)))
    else:
        min_electrodes_threshold = int(min_active_electrodes)
    
    # Step 2: Find network bursts by detecting temporal overlap
    all_spike_times = np.array(sorted(all_spike_times))
    if len(all_spike_times) == 0:
        return None
    
    recording_start = all_spike_times[0]
    recording_end = all_spike_times[-1]
    recording_duration = recording_end - recording_start
    
    # Slide a window through time and count how many electrodes are bursting
    time_step = time_window / 2  # Overlap windows by 50%
    time_points = np.arange(recording_start, recording_end, time_step)
    
    electrodes_bursting_count = []
    
    for t in time_points:
        count = 0
        for electrode_id, burst_data in electrode_bursts.items():
            for burst_start, burst_end, num_spikes in burst_data['bursts']:
                # Check if this burst overlaps with current time window
                if burst_start <= t + time_window and burst_end >= t:
                    count += 1
                    break  # Count each electrode only once per time point
        electrodes_bursting_count.append(count)
    
    electrodes_bursting_count = np.array(electrodes_bursting_count)
    
    # Step 3: Identify network bursts where enough electrodes burst simultaneously
    network_burst_active = electrodes_bursting_count >= min_electrodes_threshold
    
    # Find continuous regions of network bursting
    network_bursts = []
    in_network_burst = False
    current_burst_start = None
    current_burst_max_electrodes = 0
    
    for i, is_active in enumerate(network_burst_active):
        if is_active and not in_network_burst:
            # Start of network burst
            in_network_burst = True
            current_burst_start = time_points[i]
            current_burst_max_electrodes = electrodes_bursting_count[i]
        elif is_active and in_network_burst:
            # Continue network burst
            current_burst_max_electrodes = max(current_burst_max_electrodes, 
                                               electrodes_bursting_count[i])
        elif not is_active and in_network_burst:
            # End of network burst
            in_network_burst = False
            current_burst_end = time_points[i-1] + time_window
            
            # Count total spikes in this network burst
            spikes_in_burst = np.sum((all_spike_times >= current_burst_start) & 
                                    (all_spike_times <= current_burst_end))
            
            network_bursts.append((current_burst_start, current_burst_end, 
                                  spikes_in_burst, current_burst_max_electrodes))
    
    # Check if recording ended during a network burst
    if in_network_burst:
        current_burst_end = time_points[-1] + time_window
        spikes_in_burst = np.sum((all_spike_times >= current_burst_start) & 
                                (all_spike_times <= current_burst_end))
        network_bursts.append((current_burst_start, current_burst_end, 
                              spikes_in_burst, current_burst_max_electrodes))
    
    # Step 4: Calculate network burst statistics
    if len(network_bursts) > 0:
        burst_durations = [nb[1] - nb[0] for nb in network_bursts]
        spikes_per_burst = [nb[2] for nb in network_bursts]
        electrodes_per_burst = [nb[3] for nb in network_bursts]
        
        if len(network_bursts) > 1:
            ibis = [network_bursts[i+1][0] - network_bursts[i][1] 
                   for i in range(len(network_bursts)-1)]
            mean_ibi = np.mean(ibis)
            cov_ibi = np.std(ibis) / mean_ibi if mean_ibi > 0 else np.nan
        else:
            mean_ibi = np.nan
            cov_ibi = np.nan
        
        total_burst_spikes = sum(spikes_per_burst)
        network_burst_percentage = (total_burst_spikes / len(all_spike_times)) * 100
        
        network_data = {
            'num_network_bursts': len(network_bursts),
            'network_burst_rate': len(network_bursts) / recording_duration,
            'mean_network_burst_duration': np.mean(burst_durations),
            'std_network_burst_duration': np.std(burst_durations),
            'cov_network_burst_duration': np.std(burst_durations) / np.mean(burst_durations) if np.mean(burst_durations) > 0 else np.nan,
            'mean_ibi': mean_ibi,
            'cov_ibi': cov_ibi,
            'mean_spikes_per_network_burst': np.mean(spikes_per_burst),
            'mean_electrodes_per_network_burst': np.mean(electrodes_per_burst),
            'network_burst_percentage': network_burst_percentage,
            'active_electrodes': len(active_electrodes),
            'bursting_electrodes': len(electrode_bursts),
            'total_spikes': len(all_spike_times),
            'recording_duration': recording_duration
        }
    else:
        network_data = {
            'num_network_bursts': 0,
            'network_burst_rate': 0.0,
            'mean_network_burst_duration': np.nan,
            'std_network_burst_duration': np.nan,
            'cov_network_burst_duration': np.nan,
            'mean_ibi': np.nan,
            'cov_ibi': np.nan,
            'mean_spikes_per_network_burst': np.nan,
            'mean_electrodes_per_network_burst': np.nan,
            'network_burst_percentage': 0.0,
            'active_electrodes': len(active_electrodes),
            'bursting_electrodes': len(electrode_bursts),
            'total_spikes': len(all_spike_times),
            'recording_duration': recording_duration
        }
    
    # Step 5: Plot network burst activity
    # Plot electrode participation over time
    ax.plot(time_points, electrodes_bursting_count, 'b-', linewidth=1.5, 
           label='Electrodes Bursting')
    
    # Mark threshold
    ax.axhline(y=min_electrodes_threshold, color='gray', linestyle='--', 
              linewidth=1.5, alpha=0.7, 
              label=f'Threshold ({min_electrodes_threshold} electrodes)')
    
    # Highlight network bursts
    for nb_start, nb_end, nb_spikes, nb_electrodes in network_bursts:
        ax.axvspan(nb_start, nb_end, alpha=0.3, color='gray')
    
    # Mark network burst peaks
    if len(network_bursts) > 0:
        burst_centers = [(nb[0] + nb[1]) / 2 for nb in network_bursts]
        burst_peaks = [nb[3] for nb in network_bursts]
        ax.plot(burst_centers, burst_peaks, 'r*', markersize=3, 
               label=f'Network Bursts (n={len(network_bursts)})')
    
    # Set axis properties
    ax.set_xlim([recording_start, recording_end])
    ax.set_ylim([0, len(active_electrodes) * 1.1])
    ax.set_ylabel('Number of Bursting Electrodes', fontsize=11)
    ax.set_xlabel('Time (s)', fontsize=11)
    ax.set_title(f'Network Burst Detection (LogISI) - {len(network_bursts)} Network Bursts', 
                fontsize=12, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)
    
    return network_data
