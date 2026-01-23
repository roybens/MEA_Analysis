import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, peak_widths
import matplotlib.pyplot as plt

def compute_network_bursts(
    ax_raster,
    ax_macro,
    SpikeTimes,

    # ---------------- core params ----------------
    isi_threshold=0.1,
    bin_ms=10,
    gamma=2.0,

    # ---------------- smoothing ----------------
    smoothing_min_ms=20,

    # ---------------- burstlet filters ----------------
    min_burstlet_participation=0.20,
    min_absolute_rate_Hz=0.5,    
    min_burst_density_Hz=1.0,    
    min_relative_height=0.25,

    # ---------------- merging ----------------
    burstlet_merge_gap_s=0.1,
    network_merge_gap_s=1.0,
    superburst_min_duration_s=2.5,

    plot=True,
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
    # ---------------------------
    all_isis = []
    N = 2
    for u in units:
        t = np.unique(np.sort(SpikeTimes[u]))
        if len(t) >= N:
            isin = t[N-1:] - t[:-N+1]
            isin = isin[isin > 0]
            all_isis.extend(isin.tolist())
    biological_isi_s = np.median(all_isis) if all_isis else 0.1
    adaptive_bin_ms = np.clip(0.5 * biological_isi_s * 1000.0, 1.0, bin_ms)
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
    W = P ** gamma

    # ---------------------------
    # 3. Dual Smoothing
    # ---------------------------
    sigma_min_bins = max(1.0, smoothing_min_ms / adaptive_bin_ms)
    sigma_slow_bins = max(sigma_min_bins, 5.0 * biological_isi_s / bin_size)
    sigma_fast_bins = np.clip(sigma_slow_bins / 5.0, 1.0, 5.0)
    
    ws_onset = gaussian_filter1d(W, sigma_fast_bins)
    ws_peak  = gaussian_filter1d(W, sigma_slow_bins)
    
    burstlet_merge_gap_s = 2.0 * biological_isi_s
    network_merge_gap_s  = 5.0 * biological_isi_s

    # ---------------------------------------------------------
    # 4. Burstlets via Peak Finding (NEW)
    # ---------------------------------------------------------
    # Calculate Threshold
    #global_max_height = np.percentile(ws_onset, 99) if len(ws_onset) > 0 else 0
    global_max_height = np.max(ws_onset) if len(ws_onset) > 0 else 0
    relative_threshold_val = global_max_height * min_relative_height

    # 1. Find Peaks (The "Summits")
    # We use a tiny height threshold just to avoid finding peaks in absolute zero noise
    # The real filtering happens in the GATES below.
    min_peak_height = np.max(ws_onset) * 0.01 
    peaks, properties = find_peaks(ws_onset, height=min_peak_height)
    
    # 2. Find Boundaries (The "Base" of the mountain)
    # rel_height=0.95 means "Measure width at 95% down from the peak" (near the bottom)
    results_width = peak_widths(ws_onset, peaks, rel_height=0.80)
    
    # These indices are interpolated (floats), so we convert to int
    starts_idx = np.floor(results_width[2]).astype(int)
    ends_idx   = np.ceil(results_width[3]).astype(int)
    
    # Clip to array bounds
    starts_idx = np.clip(starts_idx, 0, len(ws_onset)-1)
    ends_idx   = np.clip(ends_idx, 0, len(ws_onset)-1)

    burstlets = []
    candidate_log = []

    # Loop through the detected peaks
    for i in range(len(peaks)):
        s, e = starts_idx[i], ends_idx[i]
        
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
        
        peak_fast = np.max(ws_onset[s:e]) if (e > s) else 0
        local_slow_baseline = np.max(ws_peak[s:e]) if (e > s) else 0

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
            "synchrony_energy": float(np.sum(ws_peak[s:e]) * bin_size),
            "participation": participation_frac,
            "total_spikes": int(total_spikes)
        })

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
        return {
            "start": s, "end": e, "duration_s": e - s,
            "peak_synchrony": max(ev["peak_synchrony"] for ev in events),
            "synchrony_energy": sum(ev["synchrony_energy"] for ev in events),
            "fragment_count": len(events),
            "total_spikes": sum(ev["total_spikes"] for ev in events),
            "participation": float(np.mean([ev["participation"] for ev in events]))
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
            "rate_bpm": len(events) / total_dur * 60.0,
            "duration": stats([e["duration_s"] for e in events]),
            "inter_event_interval": stats(np.diff(starts)) if len(starts) > 1 else stats([]),
            "intensity": stats([e["synchrony_energy"] for e in events]),
            "participation": stats([e["participation"] for e in events])
        }

    if plot and ax_macro is not None:
        ax_macro.plot(t_centers, ws_peak, color="tab:blue", lw=2)
        ax_macro.plot(t_centers, ws_onset, color="tab:orange", lw=1)
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
            #"candidate_event_log": candidate_log
        }
    }