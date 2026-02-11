import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

def compute_network_bursts(
    ax_raster=None, ax_macro=None, SpikeTimes=None,
    isi_threshold=0.1, bin_ms=10, gamma=1.0,
    smoothing_min_ms=20, min_burstlet_participation=0.20,
    min_absolute_rate_Hz=0.5, min_burst_density_Hz=1.0,
    min_relative_height=0.1, extent_frac=0.30,
    burstlet_merge_gap_s=0.1, network_merge_gap_s=1.0,
    superburst_min_duration_s=2.5, plot=False, verbose=True
):
    # ---------------------------------------------------------
    # 0. Sanity & Signal Initialization
    # ---------------------------------------------------------
    units = list(SpikeTimes.keys())
    if not units: return {"error": "no_units"}
    all_spikes = np.sort(np.concatenate([SpikeTimes[u] for u in units]))
    if all_spikes.size == 0: return {"error": "no_spikes"}
    rec_start, rec_end = float(all_spikes[0]), float(all_spikes[-1])
    total_dur = rec_end - rec_start

    # 1. Biological Calibration (Adaptive Binning)
    all_log_isis = []
    for u in units:
        t = np.unique(np.sort(SpikeTimes[u]))
        if len(t) >= 2:
            isin = np.diff(t)
            all_log_isis.extend(np.log10(isin[isin > 0]))

    biological_isi_s = 10 ** np.median(all_log_isis) if all_log_isis else 0.1
    adaptive_bin_ms = np.clip(biological_isi_s * 1000.0, 10.0, 30.0)
    bin_size = adaptive_bin_ms / 1000.0
    bins = np.arange(rec_start, rec_end + bin_size, bin_size)
    t_centers = (bins[:-1] + bins[1:]) / 2.0

    # 2. Synchrony Signals
    P = np.zeros(len(t_centers))
    PFR = np.zeros(len(t_centers))
    for u in units:
        counts, _ = np.histogram(SpikeTimes[u], bins=bins)
        P += (counts > 0).astype(float)
        PFR += counts
    PFR /= bin_size 
    W, W_gamma = P.astype(float), P ** gamma

    # 3. Dual Smoothing
    isi_bins = biological_isi_s / bin_size
    sigma_fast = np.clip(1.0 * isi_bins, 1.0, 2.0)
    sigma_slow = np.clip(5.0 * isi_bins, 3.0, 8.0)
    ws_sharp = gaussian_filter1d(W, sigma_fast)
    ws_smooth = gaussian_filter1d(W_gamma, sigma_slow)
    
    # Merge gap adaptation
    burstlet_merge_gap_s = 3 * biological_isi_s
    network_merge_gap_s  = 10 * biological_isi_s

    # ---------------------------------------------------------
    # 4. Permissible Detection Threshold & Quorum
    # ---------------------------------------------------------
    # Participation Quorum: Scaled for sparse vs dense populations
    participation_floor = max(5, 0.15 * len(units)) if len(units) < 50 else max(10, 0.05 * len(units))
    baseline_val = np.median(ws_sharp)
    
    # Detection Threshold (1.2x Baseline)
    permissible_detect_thresh = baseline_val * 1.2
    relative_threshold_val = max(participation_floor, permissible_detect_thresh)

    # ---------------------------------------------------------
    # 5. Peak Finding & Burstlet Extraction
    # ---------------------------------------------------------
    n = len(ws_sharp) 
    # Use standard MAD for finding candidate summits
    spread_mad = np.median(np.abs(ws_sharp - baseline_val))
    min_peak_height = baseline_val + (0.5 * spread_mad)
    peaks, _ = find_peaks(ws_sharp, height=min_peak_height)
    
    burstlets = []
    for p in peaks:
        # Summit-to-Base boundary expansion
        s, e = p, p
        while s > 1 and ws_sharp[s - 1] <= ws_sharp[s]: s -= 1
        while e < n - 2 and ws_sharp[e + 1] <= ws_sharp[e]: e += 1
        
        start_idx, end_idx = s, min(e, n - 1)
        start_t, end_t = t_centers[start_idx], t_centers[end_idx]
        duration_s = end_t - start_t
        if duration_s <= 0: continue

        participating = sum(1 for u in units if np.any((SpikeTimes[u] >= start_t) & (SpikeTimes[u] <= end_t)))
        participation_frac = participating / len(units)
        
        # Detection Gating
        if ws_sharp[p] < relative_threshold_val: continue
        if participation_frac < min_burstlet_participation: continue
        
        # Metric Calculations
        total_spikes = int(np.sum(PFR[start_idx:end_idx+1]) * bin_size)
        denom = (duration_s * max(1, participating))
        burst_density = total_spikes / denom if denom > 0 else 0
        peak_rate_hz_per_unit = np.max(P[start_idx:end_idx+1]) / bin_size / len(units)

        if burst_density < min_burst_density_Hz or peak_rate_hz_per_unit < min_absolute_rate_Hz: continue

        burstlets.append({
            "start": float(start_t), "end": float(end_t), "duration_s": float(duration_s),
            "peak_synchrony": float(ws_sharp[p]), "peak_time": float(t_centers[p]),
            "synchrony_energy": float(np.sum(ws_smooth[start_idx:end_idx+1]) * bin_size),
            "participation": participation_frac, "total_spikes": total_spikes,
            "burst_peak": float(np.max(PFR[start_idx:end_idx+1]))
        })

    # ---------------------------------------------------------
    # 6. Permissive "Fall to Baseline" Merging
    # ---------------------------------------------------------
    # Merge Floor: 1.05x noise floor to prevent fragmentation
    merge_floor = baseline_val * 1.05 

    def finalize(evs, s, e):
        best = max(evs, key=lambda x: x["peak_synchrony"])
        return {
            "start": s, "end": e, "duration_s": e - s,
            "peak_synchrony": best["peak_synchrony"], "peak_time": best["peak_time"], 
            "synchrony_energy": sum(ev["synchrony_energy"] for ev in evs),
            "fragment_count": len(evs), "total_spikes": sum(ev["total_spikes"] for ev in evs),
            "participation": float(np.mean([ev["participation"] for ev in evs])),
            "burst_peak": max(ev["burst_peak"] for ev in evs)
        }

    def merge(events, gap, floor_val, min_dur=0):
        if not events: return []
        events = sorted(events, key=lambda x: x['start'])
        merged, curr = [], [events[0]]
        s, e = events[0]["start"], events[0]["end"]

        for nxt in events[1:]:
            mid_start_idx, mid_end_idx = int(e / bin_size), int(nxt["start"] / bin_size)
            valley_min = np.min(ws_sharp[mid_start_idx:mid_end_idx]) if mid_end_idx > mid_start_idx else ws_sharp[mid_start_idx]
            
            # Join if within temporal gap AND quorum stayed above merge floor
            if nxt["start"] <= e + gap and valley_min >= floor_val:
                curr.append(nxt)
                e = max(e, nxt["end"])
            else:
                merged.append(finalize(curr, s, e))
                curr, s, e = [nxt], nxt["start"], nxt["end"]
        merged.append(finalize(curr, s, e))
        return [m for m in merged if m["duration_s"] >= min_dur]

    network_bursts = merge(burstlets, burstlet_merge_gap_s, merge_floor)
    superbursts = merge(network_bursts, network_merge_gap_s, merge_floor, superburst_min_duration_s)

    # ---------------------------------------------------------
    # 7. Metrics Summary & Return
    # ---------------------------------------------------------
    def stats(x):
        x = np.asarray(x)
        if x.size < 2: return {"mean": float(x.mean()) if x.size else 0.0, "std": 0.0, "cv": 0.0}
        return {"mean": float(x.mean()), "std": float(x.std()), "cv": float(x.std() / x.mean())}

    def level_metrics(events):
        if not events: return {}
        starts = [ev["start"] for ev in events]
        return {
            "count": len(events), "rate": len(events) / total_dur,
            "duration": stats([ev["duration_s"] for ev in events]),
            "inter_event_interval": stats(np.diff(starts)) if len(starts) > 1 else stats([]),
            "intensity": stats([ev["synchrony_energy"] for ev in events]),
            "participation": stats([ev["participation"] for ev in events]),
            "spikes_per_burst": stats([ev["total_spikes"] for ev in events]),
            "burst_peak": stats([ev["burst_peak"] for ev in events]),
            "peak_synchrony": stats([ev["peak_synchrony"] for ev in events])
        }

    if plot and ax_macro is not None:
        ax_macro.plot(t_centers, ws_smooth, color="tab:blue", lw=2)
        ax_macro.plot(t_centers, ws_sharp, color="tab:orange", lw=1)
        ax_macro.axhline(relative_threshold_val, color="red", ls="--", alpha=0.6, label="Detect Thresh")
        ax_macro.axhline(merge_floor, color="green", ls=":", alpha=0.6, label="Merge Floor")
        for b in network_bursts: ax_macro.axvspan(b["start"], b["end"], color="tab:blue", alpha=0.25)
        for s in superbursts: ax_macro.axvspan(s["start"], s["end"], color="purple", alpha=0.3)

    # ---------------------------------------------------------
    # 8. Registered Diagnostics
    # ---------------------------------------------------------
    return {
        "burstlets": {"events": burstlets, "metrics": level_metrics(burstlets)},
        "network_bursts": {"events": network_bursts, "metrics": level_metrics(network_bursts)},
        "superbursts": {"events": superbursts, "metrics": level_metrics(superbursts)},
        "diagnostics": {
            "adaptive_bin_ms": adaptive_bin_ms,
            "biological_isi_s": biological_isi_s,
            "baseline_value": float(baseline_val),
            "spread_mad": float(spread_mad),
            "participation_quorum_floor": float(participation_floor),
            "threshold_detect_1.2x": float(relative_threshold_val),
            "threshold_merge_1.05x": float(merge_floor),
            "burstlet_merge_gap_s": float(burstlet_merge_gap_s),
            "network_merge_gap_s": float(network_merge_gap_s),
            "merge_floor": float(merge_floor)
        },
        "plot_data": {
            "t": t_centers, "signal": ws_sharp, "signal_smooth": ws_smooth,
            "burst_peak_times": np.array([b["peak_time"] for b in network_bursts]),
            "burst_peak_values": np.array([b["peak_synchrony"] for b in network_bursts]),
            "baseline": baseline_val, "threshold": relative_threshold_val
            
        }
    }