import numpy as np
from scipy.stats import norm
from scipy.signal import find_peaks, convolve
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def plot_network_activity(
                         ax=None,
                         SpikeTimes=None,
                         binSize=0.01,  # 10 ms bins (improved from 0.1)
                         gaussianSigma=0.1,  # 100 ms sigma (improved from 0.16)
                         thresholdBurst=3.0,  # Multiple of std for adaptive threshold
                         min_peak_distance=1.0,  # seconds
                         min_active_threshold=0.1,  # Hz for active electrode definition
                         onset_offset_threshold_factor=0.3,  # Fraction of peak for burst boundaries
                         figSize=(10, 6),
                         do_plot=True):
    """
    Detect and analyze network bursts in MEA data using Gaussian smoothing and adaptive peak detection.
    
    Parameters:
    -----------
    ax : matplotlib axis object
        Axis to plot on
    SpikeTimes : dict
        Dictionary where keys are unit/electrode IDs and values are arrays of spike times (in seconds)
    binSize : float
        Bin size for spike histogram in seconds (default: 0.01 = 10 ms)
    gaussianSigma : float
        Standard deviation of Gaussian smoothing kernel in seconds (default: 0.1 = 100 ms)
    thresholdBurst : float
        Threshold multiplier for adaptive burst detection (mean + threshold*std)
    min_peak_distance : float
        Minimum time between detected peaks in seconds
    min_active_threshold : float
        Minimum firing rate (Hz) for an electrode to be considered active
    onset_offset_threshold_factor : float
        Threshold level as fraction of peak height for determining burst boundaries
    
    Returns:
    --------
    network_data : dict
        Dictionary containing network burst metrics including:
        - num_network_bursts: Number of detected network bursts
        - network_burst_rate: Bursts per second
        - mean_IBI: Mean inter-burst interval (seconds)
        - CoV_IBI: Coefficient of variation of IBI (regularity measure)
        - mean_network_burst_peak: Mean peak firing rate during bursts (Hz)
        - mean_network_burst_duration: Mean burst duration (seconds)
        - CoV_burst_duration: Coefficient of variation of burst duration
        - CoV_burst_peak: Coefficient of variation of burst peaks
        - mean_spikes_per_burst: Average number of spikes per burst
        - network_burst_percentage: Percentage of spikes in bursts
        - active_electrodes: Number of active electrodes
        - baseline_firing_rate: Baseline firing rate (Hz)
        - baseline_std: Standard deviation of baseline firing rate
    """
    
    # Step 1: Collect all spike times and count active units
    relativeSpikeTimes = []
    active_units = 0
    
    if SpikeTimes is None:
        SpikeTimes = {}

    for unit_id, spike_times in SpikeTimes.items():
        if len(spike_times) == 0:
            continue
        relativeSpikeTimes.extend(spike_times)
        
        # Calculate firing rate for this unit
        recording_duration = max(spike_times) - min(spike_times) if len(spike_times) > 1 else 1.0
        firing_rate = len(spike_times) / recording_duration
        
        if firing_rate >= min_active_threshold:
            active_units += 1
    
    if len(relativeSpikeTimes) == 0 or active_units == 0:
        print("Warning: No active units detected!")
        return None
    
    relativeSpikeTimes = np.array(relativeSpikeTimes)
    relativeSpikeTimes.sort()
    
    # Step 2: Bin spike times with appropriate resolution
    timeVector = np.arange(min(relativeSpikeTimes), max(relativeSpikeTimes) + binSize, binSize)
    binnedTimes, _ = np.histogram(relativeSpikeTimes, bins=timeVector)
    
    # Step 3: Create Gaussian kernel with proper sampling
    # Extend to ±4 sigma for better approximation (covers 99.99% of distribution)
    kernelRange = np.arange(-4*gaussianSigma, 4*gaussianSigma + binSize, binSize)
    kernel = norm.pdf(kernelRange, 0, gaussianSigma)
    kernel = kernel / np.sum(kernel)  # Proper normalization to preserve spike counts
    
    # Step 4: Smooth the binned spike times
    smoothedSpikes = convolve(binnedTimes, kernel, mode='same')
    
    # Convert to firing rate (Hz)
    # This is the network-wide firing rate (total spikes/sec across all active electrodes)
    firingRate = smoothedSpikes / binSize
    
    # Alternative: Mean firing rate across active electrodes
    # firingRate_mean = firingRate / active_units
    
    # Ensure timeVector matches firingRate length
    timeVector = timeVector[:-1]  # Remove last bin edge to match histogram output
    
    # Step 5: Calculate adaptive threshold for peak detection
    baseline_mean = np.mean(firingRate)
    baseline_std = np.std(firingRate)
    baseline_median = np.median(firingRate)
    rmsFiringRate = np.sqrt(np.mean(firingRate**2))
    
    # Adaptive peak threshold
    peak_threshold = baseline_mean + thresholdBurst * baseline_std
    
    # Step 6: Detect peaks with adaptive threshold
    min_distance_samples = int(min_peak_distance / binSize)
    
    peaks, properties = find_peaks(
        firingRate, 
        height=peak_threshold,  # Minimum height threshold
        prominence=baseline_std,  # Prominence relative to noise level
        distance=min_distance_samples  # Minimum distance in samples
    )
    
    if len(peaks) == 0:
        print("Warning: No network bursts detected!")
        network_data = {
            "num_network_bursts": 0,
            "network_burst_rate": 0.0,
            "mean_IBI": np.nan,
            "CoV_IBI": np.nan,
            "mean_network_burst_peak": np.nan,
            "mean_network_burst_duration": np.nan,
            "CoV_burst_duration": np.nan,
            "CoV_burst_peak": np.nan,
            "mean_spikes_per_burst": np.nan,
            "network_burst_percentage": 0.0,
            "active_electrodes": active_units,
            "baseline_firing_rate": baseline_mean,
            "baseline_std": baseline_std,
            "burst_onset_times": np.array([]),
            "burst_offset_times": np.array([]),
            "burst_peak_times": np.array([]),
            "burst_peak_values": np.array([]),
            "peak_threshold": float(peak_threshold),
        }
        return network_data
    
    burstPeakTimes = timeVector[peaks]
    burstPeakValues = firingRate[peaks]
    
    # Step 7: Calculate burst onset and offset times using threshold crossing
    # Define threshold as fraction of peak height (or alternatively: baseline_mean + 1*baseline_std)
    burst_onset_times = []
    burst_offset_times = []
    burst_durations = []
    spikes_per_burst = []
    
    for peak_idx, peak_time in zip(peaks, burstPeakTimes):
        peak_value = firingRate[peak_idx]
        threshold_level = baseline_mean + 1.0 * baseline_std  # Alternative: peak_value * onset_offset_threshold_factor
        
        # Find onset (search backward from peak)
        onset_idx = peak_idx
        while onset_idx > 0 and firingRate[onset_idx] > threshold_level:
            onset_idx -= 1
        
        # Find offset (search forward from peak)
        offset_idx = peak_idx
        while offset_idx < len(firingRate) - 1 and firingRate[offset_idx] > threshold_level:
            offset_idx += 1
        
        onset_time = timeVector[onset_idx]
        offset_time = timeVector[offset_idx]
        duration = offset_time - onset_time
        
        burst_onset_times.append(onset_time)
        burst_offset_times.append(offset_time)
        burst_durations.append(duration)
        
        # Count spikes in this burst
        burst_spikes = np.sum((relativeSpikeTimes >= onset_time) & 
                             (relativeSpikeTimes <= offset_time))
        spikes_per_burst.append(burst_spikes)
    
    burst_onset_times = np.array(burst_onset_times)
    burst_offset_times = np.array(burst_offset_times)
    burst_durations = np.array(burst_durations)
    spikes_per_burst = np.array(spikes_per_burst)
    
    # Step 8: Calculate Inter-Burst Intervals (IBI)
    # IBI = time from offset of one burst to onset of next burst
    if len(burst_offset_times) > 1:
        ibi_intervals = burst_onset_times[1:] - burst_offset_times[:-1]
        mean_ibi = np.mean(ibi_intervals)
        std_ibi = np.std(ibi_intervals)
        cov_ibi = std_ibi / mean_ibi if mean_ibi > 0 else np.nan
    else:
        ibi_intervals = np.array([])
        mean_ibi = np.nan
        cov_ibi = np.nan
    
    # Step 9: Calculate Coefficient of Variation for burst properties
    mean_burst_duration = np.mean(burst_durations)
    std_burst_duration = np.std(burst_durations)
    cov_burst_duration = std_burst_duration / mean_burst_duration if mean_burst_duration > 0 else np.nan
    
    mean_peak_height = np.mean(burstPeakValues)
    std_peak_height = np.std(burstPeakValues)
    cov_burst_peak = std_peak_height / mean_peak_height if mean_peak_height > 0 else np.nan
    
    # Step 10: Calculate network burst percentage
    total_spikes = len(relativeSpikeTimes)
    spikes_in_bursts = np.sum(spikes_per_burst)
    network_burst_percentage = (spikes_in_bursts / total_spikes * 100) if total_spikes > 0 else 0.0
    
    # Step 11: Calculate network burst rate
    recording_duration = timeVector[-1] - timeVector[0]
    network_burst_rate = len(peaks) / recording_duration  # bursts per second
    
    if do_plot and ax is not None:
        # Step 12: Plot the network activity
        ax.plot(timeVector, firingRate, color='royalblue', linewidth=1.2, label='Network Firing Rate')

        # Plot detected peaks
        ax.plot(burstPeakTimes, burstPeakValues, 'r*', markersize=3,
                label=f'Detected Bursts (n={len(peaks)})')

        # Plot threshold line
        ax.axhline(y=peak_threshold, color='brown', linestyle='--', linewidth=1,
                   alpha=0.7, label=f'Threshold ({thresholdBurst}×σ)')

        # Plot baseline
        ax.axhline(y=baseline_mean, color='green', linestyle=':', linewidth=1,
                   alpha=0.5, label='Baseline Mean')

        # Shade burst regions
        for onset, offset in zip(burst_onset_times, burst_offset_times):
            ax.axvspan(onset, offset, alpha=0.2, color='gray')

        # Set axis properties
        ax.set_xlim([timeVector[0], timeVector[-1]])
        ax.set_ylim([0, max(firingRate) * 1.1])
        ax.set_ylabel('Network Firing Rate (Hz)', fontsize=11)
        ax.set_xlabel('Time (s)', fontsize=11)
        ax.set_title(f'Network Activity - {len(peaks)} Network Bursts Detected', fontsize=12, fontweight='bold')
        ax.legend(loc='upper right', fontsize=9)
        ax.grid(True, alpha=0.3)
    
    # Step 13: Compile network metrics
    network_data = {
        "num_network_bursts": len(peaks),
        "network_burst_rate": network_burst_rate,  # Hz
        "mean_IBI": mean_ibi,  # seconds
        "CoV_IBI": cov_ibi,  # Coefficient of Variation (not covariance!)
        "mean_network_burst_peak": mean_peak_height,  # Hz
        "mean_network_burst_duration": mean_burst_duration,  # seconds
        "CoV_burst_duration": cov_burst_duration,
        "CoV_burst_peak": cov_burst_peak,
        "mean_spikes_per_burst": np.mean(spikes_per_burst),
        "network_burst_percentage": network_burst_percentage,  # percentage
        "active_electrodes": active_units,
        "baseline_firing_rate": baseline_mean,  # Hz
        "baseline_std": baseline_std,  # Hz
        "burst_onset_times": burst_onset_times,
        "burst_offset_times": burst_offset_times,
        "burst_peak_times": burstPeakTimes,
        "burst_peak_values": burstPeakValues,
        "peak_threshold": float(peak_threshold),
    }

    if do_plot and ax is not None:
        ax.plot(burstPeakTimes, burstPeakValues, 'or')  # Plot burst peaks as red circles
    return network_data