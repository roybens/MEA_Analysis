import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, peak_widths
import matplotlib.pyplot as plt
from . import helper_functions as helper

def compute_network_bursts(
    ax_raster=None,
    ax_macro=None,
    SpikeTimes=None,

    # ---------------- core params ----------------
    isi_threshold=0.1,
    bin_ms=10,
    gamma=1.0,

    # ---------------- smoothing ----------------
    smoothing_min_ms=20,

    # ---------------- burstlet filters ----------------
    min_burstlet_participation=0.20,
    min_absolute_rate_Hz=0.5,    
    min_burst_density_Hz=1.0,    
    min_relative_height=0.1,
    extent_frac=0.30,

    # ---------------- merging ----------------
    burstlet_merge_gap_s=0.1,
    network_merge_gap_s=1.0,
    superburst_min_duration_s=2.5,

    plot=False,
    verbose=True
):
    # ... [Sanity, ISI-N, and Signal Generation stay EXACTLY the same] ...
    # ---------------------------
    # 0. Sanity
    # ---------------------------
    units = list(SpikeTimes.keys())
    if not units: return {"error": "no_units"}
    all_spikes = np.sort(np.concatenate([SpikeTimes[u] for u in units]))
    if all_spikes.size == 0: return {"error": "no_spikes"}
    rec_start, rec_end = float(all_spikes[0]), float(all_spikes[-1])
    total_dur = rec_end - rec_start

    # ---------------------------
    # 1. Calibration
    all_log_isis = []
    N = 2

    for u in units:
        t = np.unique(np.sort(SpikeTimes[u]))
        if len(t) >= N:
            isin = t[N-1:] - t[:-N+1]
            isin = isin[isin > 0]
            all_log_isis.extend(np.log10(isin))

    biological_isi_s = (
        10 ** np.median(all_log_isis)
        if len(all_log_isis) > 0
        else 0.1
    )
    adaptive_bin_ms = np.clip(
    biological_isi_s * 1000.0,   # not 0.5 × ISI
    10.0,
    30.0
    )
    bin_size = adaptive_bin_ms / 1000.0
    bins = np.arange(rec_start, rec_end + bin_size, bin_size)
    t_centers = (bins[:-1] + bins[1:]) / 2.0

    # ---------------------------
    # 2. Network Synchrony (P and W)
    # ---------------------------
    P = np.zeros(len(t_centers))
    for u in units:
        counts, _ = np.histogram(SpikeTimes[u], bins=bins)
        P += (counts > 0).astype(float)
    W_gamma = P ** gamma
    W = P.astype(float)
    # Population firing rate (spikes, not electrodes)
    PFR = np.zeros(len(t_centers))
    for u in units:
        counts, _ = np.histogram(SpikeTimes[u], bins=bins)
        PFR += counts
    PFR = PFR / bin_size  # Hz
    # ---------------------------
    # 3. Dual Smoothing
    # ---------------------------
    # assume bin_size already clipped to 10–30 ms

    isi_bins = biological_isi_s / bin_size

    sigma_fast_bins = np.clip(1.0 * isi_bins, 1.0, 2.0)
    sigma_slow_bins = np.clip(5.0 * isi_bins, 3.0, 8.0)

    ws_sharp = gaussian_filter1d(W, sigma_fast_bins)
    ws_smooth  = gaussian_filter1d(W_gamma, sigma_slow_bins)
    burstlet_merge_gap_s = 3 * biological_isi_s
    network_merge_gap_s  = 10* biological_isi_s

    # ---------------------------------------------------------
    # 4. Burstlets via Peak Finding (NEW)
    # ---------------------------------------------------------
    # Calculate Threshold
    #global_max_height = np.percentile(ws_onset, 99) if len(ws_onset) > 0 else 0
    # global_max_height = np.max(ws_onset) if len(ws_onset) > 0 else 0
    #relative_threshold_val = global_max_height * min_relative_height

    baseline = np.median(ws_sharp)
    spread   = np.median(np.abs(ws_sharp - baseline))  # MAD

    relative_threshold_val = baseline + 3 * spread



    # 1. Find Peaks (The "Summits")
    # We use a tiny height threshold just to avoid finding peaks in absolute zero noise
    # The real filtering happens in the GATES below.
    #min_peak_height = np.max(ws_onset) * 0.01 
    min_peak_height = baseline + 0.5 * spread
    peaks, _ = find_peaks(ws_sharp, height=min_peak_height)
    peak_filtered =[]
    peak_filtered_times =[]
    # 2. Find Boundaries (The "Base" of the mountain)
    # rel_height=0.95 means "Measure width at 95% down from the peak" (near the bottom)
    starts_idx = []
    ends_idx   = []

    valid_peaks = []
    starts_idx  = []
    ends_idx    = []

    n = len(ws_sharp)

    for p in peaks:

        # ---- expand LEFT while signal is rising ----
        s = p
        while s > 1 and ws_sharp[s - 1] <= ws_sharp[s]:
            s -= 1

        # ---- expand RIGHT while signal is falling ----
        e = p
        while e < n - 2 and ws_sharp[e + 1] <= ws_sharp[e]:
            e += 1

        if e > s:
            valid_peaks.append(p)
            starts_idx.append(s)
            ends_idx.append(e)

    valid_peaks = np.asarray(valid_peaks, dtype=int)
    starts_idx  = np.asarray(starts_idx, dtype=int)
    ends_idx    = np.asarray(ends_idx, dtype=int)

    burstlets = []
    candidate_log = []

    # Loop through the detected peaks
    for i in range(len(valid_peaks)):
        p,s, e = valid_peaks[i], starts_idx[i], ends_idx[i]
        
        # If peak is just 1 bin wide, expand it slightly to catch spikes
        if s == e: e += 1

        start_t = t_centers[s]
        end_t   = t_centers[min(e, len(t_centers)-1)]
        duration_s = end_t - start_t
        
        if duration_s <= 0: continue

        # --- Metrics & Gating (SAME AS BEFORE) ---
        participating = 0
        total_spikes = 0
        for u in units:
            ts = SpikeTimes[u]
            c = np.sum((ts >= start_t) & (ts <= end_t))
            total_spikes += c
            if c > 0: participating += 1
        
        participation_frac = participating / len(units)
        denom = (duration_s * max(1, participating))
        burst_density = total_spikes / denom if denom > 0 else 0
        raw_peak_counts = np.max(P[s:e]) if (e > s) else 0
        peak_rate_hz_per_unit = raw_peak_counts / bin_size / len(units)
        
        peak_fast = ws_sharp[p]
        peak_fast_time= t_centers[p]
        local_slow_baseline = np.max(ws_smooth[s:e]) if (e > s) else 0

        burst_peak_pfr = np.max(PFR[s:e]) if (e > s) else 0.0

        # --- Diagnostics Log ---
        candidate_info = {
            "start_t": round(float(start_t), 3),
            "end_t": round(float(end_t), 3),
            "dur_s": round(float(duration_s), 4),
            "partic": round(float(participation_frac), 3),
            "density": round(float(burst_density), 3),
            "peak_hz": round(float(peak_rate_hz_per_unit), 4),
            "status": "pending"
        }

        # --- Gates ---
        if peak_fast < relative_threshold_val:
            candidate_info["status"] = f"REJECT_HEIGHT ({peak_fast:.2f} < {relative_threshold_val:.2f})"
            candidate_log.append(candidate_info)
            continue

        if participation_frac < min_burstlet_participation:
            candidate_info["status"] = f"REJECT_PARTICIPATION ({participation_frac:.2f} < {min_burstlet_participation})"
            candidate_log.append(candidate_info)
            continue

        if burst_density < min_burst_density_Hz:
            candidate_info["status"] = f"REJECT_DENSITY ({burst_density:.2f} < {min_burst_density_Hz})"
            candidate_log.append(candidate_info)
            continue

        if peak_rate_hz_per_unit < min_absolute_rate_Hz:
            candidate_info["status"] = f"REJECT_FLOOR ({peak_rate_hz_per_unit:.2f} < {min_absolute_rate_Hz})"
            candidate_log.append(candidate_info)
            continue

        # Optional Amplitude Gate
        # if peak_fast <= local_slow_baseline:
        #      candidate_info["status"] = "REJECT_AMPLITUDE"
        #      candidate_log.append(candidate_info)
        #      continue

        candidate_info["status"] = "ACCEPTED"
        candidate_log.append(candidate_info)

        burstlets.append({
            "start": float(start_t),
            "end": float(end_t),
            "duration_s": float(duration_s),
            "peak_synchrony": float(peak_fast),
            "synchrony_energy": float(np.sum(ws_smooth[s:e]) * bin_size),
            "participation": participation_frac,
            "total_spikes": int(total_spikes),
            "peak_time": float(peak_fast_time),
            "burst_peak": float(burst_peak_pfr),

        })

        peak_filtered.append(peak_fast)
        peak_filtered_times.append(peak_fast_time)
       

    # ---------------------------
    # 5. Merging (Critical for overlaps)
    # ---------------------------
    # Since find_peaks can find two peaks on one mountain, we rely on this
    # merge function to combine overlapping burstlets.
    def merge(events, gap, min_dur=0):
        if not events: return []
        # Sort by start time just in case peaks are out of order
        events = sorted(events, key=lambda x: x['start'])
        
        merged = []
        curr = [events[0]]
        s, e = events[0]["start"], events[0]["end"]
        for nxt in events[1:]:
            # Overlap logic: if next starts before current ends + gap
            if nxt["start"] <= e + gap:
                curr.append(nxt)
                e = max(e, nxt["end"])
            else:
                if (e - s) >= min_dur: merged.append(finalize(curr, s, e))
                curr = [nxt]
                s, e = nxt["start"], nxt["end"]
        if (e - s) >= min_dur: merged.append(finalize(curr, s, e))
        return merged

    def finalize(events, s, e):
        ev_max = max(events, key=lambda x: x["peak_synchrony"])

        return {
            "start": s,
            "end": e,
            "duration_s": e - s,
            "peak_synchrony": ev_max["peak_synchrony"],
            "peak_time": ev_max["peak_time"], 
            "synchrony_energy": sum(ev["synchrony_energy"] for ev in events),
            "fragment_count": len(events),
            "total_spikes": sum(ev["total_spikes"] for ev in events),
            "participation": float(np.mean([ev["participation"] for ev in events])),
            "burst_peak": max(ev["burst_peak"] for ev in events),
        }

    network_bursts = merge(burstlets, burstlet_merge_gap_s)
    superbursts = merge(network_bursts, network_merge_gap_s, superburst_min_duration_s)

    # ... [Metrics and Plotting remain the same] ...
    def stats(x):
        x = np.asarray(x)
        if x.size < 2: return {"mean": float(x.mean()) if x.size else 0.0, "std": 0.0, "cv": 0.0}
        return {"mean": float(x.mean()), "std": float(x.std()), "cv": float(x.std() / x.mean())}

    def level_metrics(events):
        if not events: return {}
        starts = [e["start"] for e in events]
        return {
            "count": len(events),
            "rate": len(events) / total_dur,
            "duration": stats([e["duration_s"] for e in events]),
            "inter_event_interval": stats(np.diff(starts)) if len(starts) > 1 else stats([]),
            "intensity": stats([e["synchrony_energy"] for e in events]),
            "participation": stats([e["participation"] for e in events]),
            "spikes_per_burst": stats([e["total_spikes"] for e in events]),
            "burst_peak": stats([e["burst_peak"] for e in events]),

        }

    if plot and ax_macro is not None:
        ax_macro.plot(t_centers, ws_smooth, color="tab:blue", lw=2)
        ax_macro.plot(t_centers, ws_sharp, color="tab:orange", lw=1)
        # Highlight peaks found
        #ax_macro.plot(t_centers[peaks], ws_onset[peaks], "x", color="red", label="Peaks")
        
        for b in network_bursts:
            ax_macro.axvspan(b["start"], b["end"], color="tab:blue", alpha=0.25)
        for s in superbursts:
            ax_macro.axvspan(s["start"], s["end"], color="purple", alpha=0.3)

    return {
        "burstlets": {"events": burstlets, "metrics": level_metrics(burstlets)},
        "network_bursts": {"events": network_bursts, "metrics": level_metrics(network_bursts)},
        "superbursts": {"events": superbursts, "metrics": level_metrics(superbursts)},
        "diagnostics": {
            "adaptive_bin_ms": adaptive_bin_ms,
            "biological_isi_s": biological_isi_s,
            "smoothing_sigma_fast_bins": sigma_fast_bins,
            "smoothing_sigma_slow_bins": sigma_slow_bins,
            "threshold_value": relative_threshold_val,
            "burstlet_merge_gap_s": burstlet_merge_gap_s,
            "network_merge_gap_s": network_merge_gap_s
            #"candidate_event_log": candidate_log
        },
        "plot_data": {
            "t": t_centers,
            "signal": ws_sharp,
            "signal_smooth": ws_smooth,
            "burst_peak_times": np.array([b["peak_time"] for b in network_bursts]),
            "burst_peak_values": np.array([b["peak_synchrony"] for b in network_bursts]),
            "baseline": np.percentile(ws_sharp, 10),
            "threshold": relative_threshold_val
        }

    }