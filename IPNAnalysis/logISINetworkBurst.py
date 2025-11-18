"""
Burst detection using LogISI method and network burst analysis.
References:
- Pasquale V, Massobrio P, Bologna L, Chiappalone M, Martinoia S. Self-organization and neuronal avalanches in networks of dissociated cortical neurons. Neuroscience. 2010 Mar 17;153(4):1354-69. doi: 10.1016/j.neuroscience.2008.12.050. Epub 2009 Jan 7. PMID: 19110017.

Authors: Mandar M. Patil,LLM Assisted Edits: Yes ChatGPT-5


"""
import numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt

def logISI_burst_detection(all_spike_times,
                           max_isi_cutoff=0.1,
                           void_threshold=0.7,
                           min_spikes_in_burst=5,
                           hist_bins=100,
                           smooth_sigma=2.0):
    """
    NETWORK-LEVEL LogISI burst detection (correct Pasquale 2010 implementation).
    Detects network bursts from pooled spike train.

    all_spike_times : array-like
        Concatenated spike times from all units/electrodes.

    Returns:
      bursts = [(start, end, num_spikes), ...]
      isi_cutoff = chosen ISI threshold
      info = diagnostic dictionary
    """

    spike_times = np.sort(np.array(all_spike_times))
    if len(spike_times) < min_spikes_in_burst:
        return [], np.nan, {"qc":"too_few_spikes"}

    # ---- ISIs on GLOBAL network spike train ----
    isis = np.diff(spike_times)
    if len(isis) == 0:
        return [], np.nan, {"qc":"no_isi"}

    # ---- log transform ----
    log_isis = np.log10(isis + 1e-10)

    # ---- histogram ----
    hist_counts, hist_edges = np.histogram(log_isis, bins=hist_bins)
    hist_centers = (hist_edges[:-1] + hist_edges[1:]) / 2

    # ---- smoothing ----
    hist_smooth = gaussian_filter1d(hist_counts.astype(float), sigma=smooth_sigma)

    # ---- peak detection ----
    peaks, props = find_peaks(hist_smooth,
                              prominence=np.max(hist_smooth)*0.1)

    if len(peaks) == 0:
        isi_cutoff = max_isi_cutoff
    else:
        # intraburst peak = largest peak left of max_isi_cutoff
        log_cut = np.log10(max_isi_cutoff)
        intraburst_cands = peaks[hist_centers[peaks] <= log_cut]

        if len(intraburst_cands) == 0:
            isi_cutoff = max_isi_cutoff
        else:
            intraburst_peak = intraburst_cands[np.argmax(hist_smooth[intraburst_cands])]
            right_peaks = peaks[peaks > intraburst_peak]

            if len(right_peaks) == 0:
                isi_cutoff = max_isi_cutoff
            else:
                # void parameter
                cutoff_idx = None
                for rp in right_peaks:
                    region = hist_smooth[intraburst_peak:rp+1]
                    min_idx = intraburst_peak + np.argmin(region)

                    h1 = hist_smooth[intraburst_peak]
                    h2 = hist_smooth[rp]
                    hmin = hist_smooth[min_idx]

                    void = (h1 + h2 - 2*hmin) / (h1 + h2 + 1e-12)

                    if void >= void_threshold:
                        cutoff_idx = min_idx
                        break

                if cutoff_idx is None:
                    isi_cutoff = max_isi_cutoff
                else:
                    isi_cutoff = 10 ** hist_centers[cutoff_idx]
                    isi_cutoff = min(isi_cutoff, max_isi_cutoff)

    # ---- detect network bursts from GLOBAL ISI threshold ----
    bursts = []
    current = [0]

    for i, isi in enumerate(isis):
        if isi <= isi_cutoff:
            current.append(i+1)
        else:
            if len(current) >= min_spikes_in_burst:
                bursts.append((spike_times[current[0]],
                               spike_times[current[-1]],
                               len(current)))
            current = [i+1]

    if len(current) >= min_spikes_in_burst:
        bursts.append((spike_times[current[0]],
                       spike_times[current[-1]],
                       len(current)))

    info = {
        "hist_centers": hist_centers,
        "hist_smooth": hist_smooth,
        "peaks": peaks,
        "isi_cutoff": isi_cutoff
    }

    return bursts, isi_cutoff, info


def plot_network_bursts_logISI(ax, SpikeTimes,
                               max_isi_cutoff=0.1,
                               void_threshold=0.7,
                               min_spikes_in_burst=5,
                               min_active_electrodes=0.1,  # now used for participation, not threshold
                               time_window=0.1,
                               min_active_threshold=0.1):
    """
    Correct Network LogISI Burst Detector.
    Uses global spike train, not per-unit ISI.
    """

    # -------------------------------------------------------
    # 1. POOL ALL SPIKES (GLOBAL NETWORK ACTIVITY)
    # -------------------------------------------------------
    all_spikes = np.sort(np.concatenate(list(SpikeTimes.values())))
    if len(all_spikes) < 20:
        return ax, {}

    rec_start = all_spikes[0]
    rec_end   = all_spikes[-1]
    rec_dur   = rec_end - rec_start

    # -------------------------------------------------------
    # 2. COMPUTE GLOBAL ISI AND LOGISI THRESHOLD
    # -------------------------------------------------------
    bursts, isi_cutoff, info = logISI_burst_detection(
        all_spikes,
        max_isi_cutoff=max_isi_cutoff,
        void_threshold=void_threshold,
        min_spikes_in_burst=min_spikes_in_burst
    )

    # bursts returned here are *global* network bursts
    network_bursts = bursts

    # -------------------------------------------------------
    # 3. GLOBAL FIRING RATE (Pasquale-style)
    # -------------------------------------------------------
    bin_width = 0.01   # 10 ms bins  #these are for plotting only the logic comes from logISI
    bins = np.arange(rec_start, rec_end, bin_width)
    counts, _ = np.histogram(all_spikes, bins=bins)
    rate = counts / bin_width                   # spikes/s
    rate_smooth = gaussian_filter1d(rate, sigma=10)

    # -------------------------------------------------------
    # 4. DETECT NETWORK BURST PERIODS FROM RATE
    # -------------------------------------------------------
    # threshold = anything that produces the same burst windows as LogISI
    # i.e. rate > threshold -> burst
    rate_threshold = np.percentile(rate_smooth, 90)   # empirically stable

    burst_mask = rate_smooth > rate_threshold

    network_bursts = []
    in_burst = False

    for i, flag in enumerate(burst_mask):
        if flag and not in_burst:
            in_burst = True
            burst_start = bins[i]
        elif not flag and in_burst:
            in_burst = False
            burst_end = bins[i]
            duration = burst_end - burst_start
            if duration >= 0.05:  # min 50 ms
                # count spikes in this region
                n_spikes = np.sum((all_spikes >= burst_start) & (all_spikes <= burst_end))
                network_bursts.append((burst_start, burst_end, n_spikes))

    if in_burst:
        burst_end = bins[-1]
        duration = burst_end - burst_start
        if duration >= 0.05:
            n_spikes = np.sum((all_spikes >= burst_start) & (all_spikes <= burst_end))
            network_bursts.append((burst_start, burst_end, n_spikes))

    # -------------------------------------------------------
    # 5. UNIT PARTICIPATION (AFTER NETWORK BURSTS ARE KNOWN)
    # -------------------------------------------------------
    participation = {}

    for unit, spikes in SpikeTimes.items():
        spikes = np.asarray(spikes)
        parts = []
        for (s,e,n) in network_bursts:
            parts.append(int(np.any((spikes >= s) & (spikes <= e))))
        participation[unit] = parts

    # -------------------------------------------------------
    # 6. PLOT GLOBAL RATE + SHADED NETWORK BURSTS
    # -------------------------------------------------------
    t_centers = bins[:-1]

    ax.plot(t_centers, rate_smooth, 'b-', linewidth=1.2, label="Global firing rate (Hz)")
    ax.axhline(rate_threshold, color='gray', linestyle='--', label="Burst threshold")

    for (s,e,n) in network_bursts:
        ax.axvspan(s, e, color='gray', alpha=0.3)

    ax.set_xlim(rec_start, rec_end)
    ax.set_ylabel("Network firing rate (Hz)")
    ax.set_xlabel("Time (s)")
    ax.set_title(f"Network Bursts (N={len(network_bursts)}) — LogISI")
    ax.legend(loc="upper right")

    # -------------------------------------------------------
    # 7. RETURN METRICS
    # -------------------------------------------------------
    metrics = {
        "network_bursts": network_bursts,
        "isi_cutoff": isi_cutoff,
        "rate_threshold": rate_threshold,
        "participation": participation,
        "global_rate": rate_smooth
    }

    return ax, metrics

