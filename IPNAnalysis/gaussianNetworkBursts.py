import numpy as np
from scipy.stats import norm
from scipy.signal import find_peaks, convolve


def plot_network_activity(
    ax=None,
    SpikeTimes=None,
    binSize=0.01,
    gaussianSigma=0.1,
    thresholdBurst=3.0,
    min_peak_distance=1.0,
    min_active_threshold=0.1,
    onset_offset_threshold_factor=0.3,
    figSize=(10, 6),
    do_plot=True
):

    # -------------------------------
    # STEP 1 — Flatten spikes + keep channels
    # -------------------------------
    relativeSpikeTimes = []
    spike_channels = []

    active_units = 0

    if SpikeTimes is None:
        SpikeTimes = {}

    for unit_id, spike_times in SpikeTimes.items():
        if len(spike_times) == 0:
            continue

        spike_times = np.array(spike_times)
        relativeSpikeTimes.extend(spike_times)
        spike_channels.extend([unit_id] * len(spike_times))

        duration = spike_times.max() - spike_times.min() if len(spike_times) > 1 else 1.0
        fr = len(spike_times) / duration

        if fr >= min_active_threshold:
            active_units += 1

    if len(relativeSpikeTimes) == 0 or active_units == 0:
        return None

    relativeSpikeTimes = np.array(relativeSpikeTimes)
    spike_channels = np.array(spike_channels)

    sort_idx = np.argsort(relativeSpikeTimes)
    relativeSpikeTimes = relativeSpikeTimes[sort_idx]
    spike_channels = spike_channels[sort_idx]

    total_channels = len(np.unique(spike_channels))

    # -------------------------------
    # STEP 2 — Histogram + smoothing
    # -------------------------------
    timeVector = np.arange(relativeSpikeTimes.min(),
                           relativeSpikeTimes.max() + binSize,
                           binSize)

    binnedTimes, _ = np.histogram(relativeSpikeTimes, bins=timeVector)

    kernelRange = np.arange(-4*gaussianSigma, 4*gaussianSigma + binSize, binSize)
    kernel = norm.pdf(kernelRange, 0, gaussianSigma)
    kernel /= np.sum(kernel)

    smoothedSpikes = convolve(binnedTimes, kernel, mode='same')
    firingRate = smoothedSpikes / binSize
    timeVector = timeVector[:-1]

    # -------------------------------
    # STEP 3 — Adaptive threshold
    # -------------------------------
    baseline_mean = np.mean(firingRate)
    baseline_std = np.std(firingRate)

    peak_threshold = baseline_mean + thresholdBurst * baseline_std
    min_distance_samples = int(min_peak_distance / binSize)

    peaks, _ = find_peaks(
        firingRate,
        height=peak_threshold,
        prominence=baseline_std,
        distance=min_distance_samples
    )

    if len(peaks) == 0:
        return None

    peak_times = timeVector[peaks]
    peak_values = firingRate[peaks]

    # -------------------------------
    # STEP 4 — Burst boundaries (MATLAB-style)
    # -------------------------------
    burst_onsets = []
    burst_offsets = []

    for peak_idx in peaks:
        peak_val = firingRate[peak_idx]

        # MATLAB-style: relative to peak
        threshold_level = peak_val * (1 - onset_offset_threshold_factor)

        i = peak_idx
        while i > 0 and firingRate[i] > threshold_level:
            i -= 1
        onset = timeVector[i]

        j = peak_idx
        while j < len(firingRate)-1 and firingRate[j] > threshold_level:
            j += 1
        offset = timeVector[j]

        burst_onsets.append(onset)
        burst_offsets.append(offset)

    burst_onsets = np.array(burst_onsets)
    burst_offsets = np.array(burst_offsets)

    burst_durations = burst_offsets - burst_onsets

    # -------------------------------
    # STEP 5 — Assign spikes to bursts
    # -------------------------------
    spikes_per_burst = []
    participation_per_burst = []

    in_burst_mask = np.zeros_like(relativeSpikeTimes, dtype=bool)

    for onset, offset in zip(burst_onsets, burst_offsets):
        mask = (relativeSpikeTimes >= onset) & (relativeSpikeTimes <= offset)

        spikes_per_burst.append(np.sum(mask))

        ch_in_burst = spike_channels[mask]
        if len(ch_in_burst) > 0:
            participation = len(np.unique(ch_in_burst)) / total_channels
        else:
            participation = np.nan

        participation_per_burst.append(participation)

        in_burst_mask |= mask

    spikes_per_burst = np.array(spikes_per_burst)
    participation_per_burst = np.array(participation_per_burst)

    tsWithinBurst = relativeSpikeTimes[in_burst_mask]
    chWithinBurst = spike_channels[in_burst_mask]

    tsOutsideBurst = relativeSpikeTimes[~in_burst_mask]
    chOutsideBurst = spike_channels[~in_burst_mask]

    # -------------------------------
    # STEP 6 — Proper baseline (outside bursts)
    # -------------------------------
    outside_mask_time = np.ones_like(timeVector, dtype=bool)

    for onset, offset in zip(burst_onsets, burst_offsets):
        outside_mask_time &= ~((timeVector >= onset) & (timeVector <= offset))

    if np.any(outside_mask_time):
        baseline_outside = np.mean(firingRate[outside_mask_time])
    else:
        baseline_outside = baseline_mean

    # -------------------------------
    # STEP 7 — IBI
    # -------------------------------
    if len(burst_onsets) > 1:
        ibis = burst_onsets[1:] - burst_offsets[:-1]
    else:
        ibis = np.array([])

    # -------------------------------
    # STEP 8 — CoV helper
    # -------------------------------
    def cov(x):
        return np.std(x) / np.mean(x) if len(x) > 0 and np.mean(x) > 0 else np.nan

    # -------------------------------
    # STEP 9 — Metrics (MATLAB-aligned)
    # -------------------------------
    network_data = {
        # Core
        "num_network_bursts": len(peaks),
        "network_burst_rate": len(peaks) / (timeVector[-1] - timeVector[0]),

        # IBI
        "mean_IBI": np.mean(ibis) if len(ibis) else np.nan,
        "CoV_IBI": cov(ibis),

        # Peaks
        "mean_network_burst_peak": np.mean(peak_values),
        "CoV_burst_peak": cov(peak_values),

        # Duration
        "mean_network_burst_duration": np.mean(burst_durations),
        "CoV_burst_duration": cov(burst_durations),

        # Spikes
        "mean_spikes_per_burst": np.mean(spikes_per_burst),
        "CoV_spikes_per_burst": cov(spikes_per_burst),

        # Participation
        "mean_participation_ratio": np.mean(participation_per_burst),
        "CoV_participation_ratio": cov(participation_per_burst),

        # Baseline
        "baseline_firing_rate": baseline_outside,

        # Distributions (CRITICAL)
        "IBI_distribution": ibis,
        "burst_duration_distribution": burst_durations,
        "spikes_per_burst_distribution": spikes_per_burst,
        "participation_distribution": participation_per_burst,

        # Time variables (MATLAB equivalent)
        "tsWithinBurst": tsWithinBurst,
        "chWithinBurst": chWithinBurst,
        "tsOutsideBurst": tsOutsideBurst,
        "chOutsideBurst": chOutsideBurst,

        # Events
        "burst_onset_times": burst_onsets,
        "burst_offset_times": burst_offsets,
        "burst_peak_times": peak_times,
        "burst_peak_values": peak_values,
    }

    return network_data
