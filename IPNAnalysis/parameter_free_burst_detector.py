import numpy as np
import json
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d, percentile_filter

def compute_network_bursts(
    ax_raster,
    ax_macro,
    SpikeTimes,
    # ---------------- core detection params ----------------
    isi_threshold=0.1,             # seconds for unit burst detection
    bin_ms=10,                     # bin width (hardcoded to 10ms for biophysics)
    gamma=2.0,                     # Contrast amplifier (Signal^2)
    
    # ---------------- adaptive smoothing params ----------------
    smoothing_min_ms=20,           # absolute minimum smoothing
    
    # ---------------- thresholding params ----------------
    baseline_pct=50.0,             # Median baseline (handles active cultures)
    local_baseline_win_s=5.0,      # Rolling baseline window
    delta_thresh=0.25,             # Relative threshold (25% above baseline)
    
    # ---------------- peak / noise params ----------------
    max_back_ms=500.0,             # Backtracking window
    slope_k=2.0,
    min_burstlet_duration_ms=20,   # CHANGED: 20ms (2 bins) to catch sharp spikes
    
    # ---------------- merging params ----------------
    burstlet_merge_gap_s=0.1,      # CHANGED: 0.1s to separate fast stuttering
    network_merge_gap_s=1.0,       # 1.0s to group superbursts
    superburst_min_duration_s=5.0, 
    
    # ---------------- plotting / verbosity ----------------
    plot=True,
    verbose=False
):
    """
    Final Robust Burst Detection Pipeline.
    Integrates Adaptive Smoothing, Anchored Squelch, Fast-Rescue, and Literature Metrics.
    """

    # ---------------------------
    # 0) Sanity Checks & Setup
    # ---------------------------
    units = list(SpikeTimes.keys())
    if len(units) == 0:
        return {"error": "no_units"}, None

    # Flatten spikes for global stats
    all_spikes_flat = np.sort(np.concatenate([np.asarray(SpikeTimes[u]) for u in units]))
    if all_spikes_flat.size == 0:
        return {"error": "no_spikes"}, None

    rec_start, rec_end = float(all_spikes_flat[0]), float(all_spikes_flat[-1])
    total_rec_duration = rec_end - rec_start
    
    if verbose:
        print(f"Recording: {rec_start:.2f} -> {rec_end:.2f} s | Units: {len(units)}")

    # ---------------------------
    # Phase 1: Unit-level ISI Bursts (Biological Calibration)
    # ---------------------------
    def detect_unit_bursts_logic(spike_times_dict, thresh):
        out = {}
        for u, times in spike_times_dict.items():
            times = np.asarray(times)
            if times.size < 2:
                out[u] = {"bursts": [], "isis_all": []}
                continue
            
            isis = np.diff(times)
            mask = isis < thresh
            
            # Find runs of short ISIs
            starts = np.where(np.diff(np.insert(mask.astype(int), 0, 0)) == 1)[0]
            ends   = np.where(np.diff(np.append(mask.astype(int), 0)) == -1)[0]
            
            bursts = [list(times[s:e+1]) for s,e in zip(starts, ends)]
            out[u] = {
                "bursts": bursts,
                "isis_all": isis.tolist() 
            }
        return out

    unit_bursts = detect_unit_bursts_logic(SpikeTimes, isi_threshold)

    # ---------------------------
    # Phase 2: Construct "Network Synchrony" Signal
    # ---------------------------
    bin_size = float(bin_ms) / 1000.0
    bins = np.arange(rec_start, rec_end + bin_size, bin_size)
    t_centers = (bins[:-1] + bins[1:]) / 2.0

    # Calculate Population Firing Rate (Participation)
    P = np.zeros(len(t_centers), dtype=float)
    for u in units:
        counts, _ = np.histogram(SpikeTimes[u], bins=bins)
        P += (counts > 0).astype(float) # Binary participation per unit per bin

    # Contrast Enhancement (Synchrony Amplifier)
    W = (P ** gamma) if gamma != 1.0 else P.copy()

    # ---------------------------
    # Phase 3: Adaptive Dual-Stream Smoothing
    # ---------------------------
    
    # A. Calculate "Biological ISI" (Tempo of the culture)
    all_unit_isis = []
    for u_data in unit_bursts.values():
        all_unit_isis.extend(u_data['isis_all'])
    
    biological_isi_s = np.median(all_unit_isis) if len(all_unit_isis) > 0 else 0.1
    
    # B. Define Kernels
    target_sigma_s = 5.0 * biological_isi_s
    sigma_min_bins = max(1.0, smoothing_min_ms / bin_ms)
    target_sigma_bins = float(target_sigma_s / bin_size)
    
    sigma_slow_bins = max(sigma_min_bins, target_sigma_bins)       # The "Integrator"
    sigma_fast_bins = max(1.0, sigma_slow_bins / 5.0)              # The "Detector"
    sigma_fast_bins = min(5.0, sigma_fast_bins)                    # Cap fast stream for sharpness

    if verbose:
        print(f"Adaptive Smoothing | Bio ISI: {biological_isi_s*1000:.1f}ms -> Sigma Slow: {sigma_slow_bins*bin_ms:.1f}ms")

    # C. Apply Filters
    ws_onset = gaussian_filter1d(W, sigma=sigma_fast_bins) # Fast Stream (Orange)
    ws_peak  = gaussian_filter1d(W, sigma=sigma_slow_bins) # Slow Stream (Blue)

    # ---------------------------
    # Phase 4: Anchored Squelch & Detection
    # ---------------------------
    
    # A. Calculate Squelch (Noise Floor)
    # Anchor to 5th percentile to find true silence even in active cultures
    noise_anchor = np.percentile(ws_peak, 5.0)
    below_median = ws_peak[ws_peak < np.median(ws_peak)]
    noise_mad = np.median(np.abs(below_median - noise_anchor)) if len(below_median) > 0 else 1.0
    
    squelch_threshold = noise_anchor + (5.0 * noise_mad)
    
    # Safety: Ensure squelch is at least 5% of max peak (prevents 0-threshold in silent files)
    min_safety = np.max(ws_peak) * 0.05
    squelch_threshold = max(squelch_threshold, min_safety)

    if verbose:
        print(f"Squelch: {squelch_threshold:.2f} (Anchor: {noise_anchor:.2f})")

    # B. Find Candidates (Relative Thresholding)
    win_bins = max(3, int(round(local_baseline_win_s / bin_size)))
    baseline_local = percentile_filter(ws_onset, percentile=baseline_pct, size=win_bins, mode='reflect')
    eps = 1e-9
    delta = (ws_onset - baseline_local) / (baseline_local + eps)
    
    above = delta > delta_thresh
    actives_idx = []
    i = 0
    while i < len(above):
        if above[i]:
            s = i
            while i < len(above) and above[i]: i += 1
            actives_idx.append((s, i-1))
        else:
            i += 1

    # C. Validate Candidates (Dual-Gate Logic)
    peaks_idx, _ = find_peaks(ws_peak, prominence=(0.1 * np.max(ws_peak)))
    peaks_idx = set(peaks_idx)

    burstlets = []
    max_back_bins = int(max_back_ms / bin_ms)
    slope_noise = 1.4826 * np.median(np.abs(np.diff(ws_onset))) if len(ws_onset)>1 else 0

    for (s_idx, e_idx) in actives_idx:
        # Define Rescue Threshold (2x Squelch for sharp events)
        fast_rescue_thresh = squelch_threshold * 2.0

        # Check 1: Slow Stream Peak (Standard)
        p_inside = [p for p in peaks_idx if s_idx <= p <= e_idx]
        
        # Check 2: Fast Stream Peak (Rescue)
        segment_fast = ws_onset[s_idx:e_idx+1]
        fast_peak_val = np.max(segment_fast) if len(segment_fast) > 0 else 0
        
        pk_idx = None
        pk_val = 0.0
        is_rescued = False

        if p_inside:
            # Path A: Found slow peak
            pk_idx = max(p_inside, key=lambda p: ws_peak[p])
            pk_val = float(ws_peak[pk_idx])
            
            # If slow peak is too weak, try fast rescue
            if pk_val < squelch_threshold:
                if fast_peak_val > fast_rescue_thresh:
                    pk_idx = s_idx + np.argmax(segment_fast) 
                    is_rescued = True
                else:
                    continue # Reject
        else:
            # Path B: No slow peak, try fast rescue immediately
            if fast_peak_val > fast_rescue_thresh:
                pk_idx = s_idx + np.argmax(segment_fast)
                pk_val = float(ws_peak[pk_idx]) 
                is_rescued = True
            else:
                continue # Reject

        # ---------------------------------------------------------
        # D. Refine Boundaries
        # ---------------------------------------------------------
        # Backtrack to Onset (Valley finding)
        search_start = max(0, pk_idx - max_back_bins)
        segment = ws_onset[search_start:pk_idx+1]
        onset_idx = search_start + np.argmin(segment)
        
        # Slope Check
        if slope_k > 0:
            rise_dur = (pk_idx - onset_idx) * bin_size
            rise_amp = ws_onset[pk_idx] - ws_onset[onset_idx]
            slope = rise_amp / rise_dur if rise_dur > 0 else 0
            if slope < (slope_k * slope_noise):
                continue 

        onset_t = float(t_centers[onset_idx])
        end_t = float(t_centers[e_idx])
        
        # Duration Check (20ms)
        if (end_t - onset_t)*1000 < min_burstlet_duration_ms:
            continue

        # E. Calculate Energy (Area under curve)
        # Sum the Slow Stream (ws_peak) for the duration of the burst
        # Indices: onset_idx to e_idx
        burst_energy = np.sum(ws_peak[onset_idx:e_idx+1]) * bin_size

        # F. Collect Stats
        total_spikes = 0
        participating_units = 0
        for u in units:
            ts = np.asarray(SpikeTimes[u])
            c = np.sum((ts >= onset_t) & (ts <= end_t))
            total_spikes += c
            if c > 0: participating_units += 1

        burstlets.append({
            "start": onset_t, "end": end_t,
            "duration_s": end_t - onset_t,
            "total_spikes": int(total_spikes),
            "participating_units": int(participating_units),
            "peak_synchrony": pk_val,       # Renamed from peak_amp
            "synchrony_energy": burst_energy, # Calculated Area
            "rescued": is_rescued
        })

    burstlets.sort(key=lambda x: x["start"])

    # ---------------------------
    # Phase 5: Hierarchical Merging
    # ---------------------------
    def merge_events(events, gap_s, min_dur_s=0):
        if not events: return []
        merged = []
        curr = events[0]
        c_start, c_end = curr["start"], curr["end"]
        sub_events = [curr]
        
        for nxt in events[1:]:
            if (nxt["start"] - c_end) < gap_s:
                # Merge
                c_end = max(c_end, nxt["end"])
                sub_events.append(nxt)
            else:
                # Finalize
                dur = c_end - c_start
                if dur >= min_dur_s:
                    # Sum energies and spikes from sub-events
                    total_energy = sum(e["synchrony_energy"] for e in sub_events)
                    total_s = sum(e["total_spikes"] for e in sub_events)
                    peak_s = max(e["peak_synchrony"] for e in sub_events)
                    
                    merged.append({
                        "start": c_start, "end": c_end, 
                        "duration_s": dur, 
                        "total_spikes": total_s,
                        "peak_synchrony": peak_s,
                        "synchrony_energy": total_energy,
                        "fragment_count": len(sub_events)
                    })
                # Start new
                c_start, c_end = nxt["start"], nxt["end"]
                sub_events = [nxt]
        
        # Finalize last
        if (c_end - c_start) >= min_dur_s:
            total_energy = sum(e["synchrony_energy"] for e in sub_events)
            total_s = sum(e["total_spikes"] for e in sub_events)
            peak_s = max(e["peak_synchrony"] for e in sub_events)
            merged.append({
                "start": c_start, "end": c_end, 
                "duration_s": c_end - c_start, 
                "total_spikes": total_s,
                "peak_synchrony": peak_s,
                "synchrony_energy": total_energy,
                "fragment_count": len(sub_events)
            })
        return merged

    network_bursts = merge_events(burstlets, gap_s=burstlet_merge_gap_s)
    superbursts = merge_events(network_bursts, gap_s=network_merge_gap_s, min_dur_s=superburst_min_duration_s)

    # ---------------------------
    # Phase 6: Plotting
    # ---------------------------
    if plot and ax_macro is not None:
        # Plot Signal Traces
        ax_macro.plot(t_centers, ws_peak, color='tab:blue', lw=1.5, label='Synchrony Envelope (Slow)')
        ax_macro.plot(t_centers, ws_onset, color='tab:orange', lw=1, alpha=0.6, label='Instantaneous Synchrony (Fast)')
        ax_macro.plot(t_centers, baseline_local, ':', color='tab:red', label='Baseline')
        
        # Visualize Squelch
        ax_macro.axhline(squelch_threshold, color='green', linestyle='--', alpha=0.5, label='Squelch (Noise Floor)')
        
        # Visualize Network Bursts (Blue Shading)
        for nb in network_bursts:
            ax_macro.axvspan(nb["start"], nb["end"], color='tab:blue', alpha=0.15)
        
        # Visualize Superbursts (Black Bar)
        for sb in superbursts:
            ax_macro.hlines(y=np.max(ws_peak)*1.05, xmin=sb["start"], xmax=sb["end"], color='k', lw=3)

        # Labels
        ax_macro.set_ylabel('Network Synchrony (a.u.)')
        ax_macro.set_xlabel('Time (s)')
        ax_macro.legend(loc='upper right', fontsize='small')
        #ax_macro.set_title(f"Network Bursts: {len(network_bursts)} | Superbursts: {len(superbursts)}")
        #have it as a annotation at right instead of title
        #ax_macro.an
        ax_macro.annotate(f"Network Bursts: {len(network_bursts)} | Superbursts: {len(superbursts)}", xy=(0.5, 1.05), xycoords='axes fraction', ha='center', fontsize='medium')

    # ---------------------------
    # Phase 7: Calculate Literature Metrics (The Final Data)
    # ---------------------------
    def calculate_lit_metrics(nbs, total_dur):
        if not nbs:
            return {"count": 0, "rate_hz": 0.0, "duration": {}, "ibi": {}, "intensity": {}}
        
        # 1. Sort
        nbs.sort(key=lambda x: x["start"])
        
        # 2. Calculate IBI (Inter-Burst Interval)
        #    IBI = Time from End of Previous to Start of Current
        ibis = []
        for i in range(1, len(nbs)):
            ibi_val = nbs[i]["start"] - nbs[i-1]["end"]
            ibis.append(ibi_val)
            nbs[i]["ibi_pre_s"] = float(ibi_val)
        if len(nbs) > 0: nbs[0]["ibi_pre_s"] = None

        # 3. Aggregators
        def get_stats(data):
            arr = np.array(data)
            if len(arr) < 2: return {"mean": float(np.mean(arr)) if len(arr)==1 else 0, "std": 0, "cv": 0}
            mean = float(np.mean(arr))
            std = float(np.std(arr))
            return {"mean": mean, "std": std, "cv": std/mean if mean > 0 else 0}

        # 4. Global Stats
        stats_dur = get_stats([n["duration_s"] for n in nbs])
        stats_ibi = get_stats(ibis)
        stats_energy = get_stats([n["synchrony_energy"] for n in nbs])
        stats_spikes = get_stats([n["total_spikes"] for n in nbs])
        
        burst_spikes = sum(n["total_spikes"] for n in nbs)
        total_spikes = len(all_spikes_flat)
        psib = (burst_spikes / total_spikes * 100.0) if total_spikes > 0 else 0

        return {
            "count": len(nbs),
            "rate_hz": len(nbs) / total_dur,
            "rate_bpm": (len(nbs) / total_dur) * 60.0,
            "duration": stats_dur,
            "inter_burst_interval": stats_ibi,
            "intensity": {
                "energy": stats_energy,
                "spikes": stats_spikes,
                "psib": psib
            }
        }

    lit_metrics = calculate_lit_metrics(network_bursts, total_rec_duration)

    # ---------------------------
    # Phase 8: Clean Output
    # ---------------------------
    def clean(obj):
        if isinstance(obj, dict): return {k: clean(v) for k,v in obj.items()}
        if isinstance(obj, list): return [clean(v) for v in obj]
        if isinstance(obj, np.integer): return int(obj)
        if isinstance(obj, np.floating): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        return obj

    return clean({
        "unit_bursts": unit_bursts,
        "burstlets": burstlets,
        "network_bursts": network_bursts,
        "superbursts": superbursts,
        "network_stats": lit_metrics,  # <--- New Standardized Metrics
        "diagnostics": {
            "sigma_slow_ms": sigma_slow_bins * bin_ms,
            "biological_isi_ms": biological_isi_s * 1000,
            "squelch_threshold": squelch_threshold
        }
    })