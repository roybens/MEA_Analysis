import numpy as np
import math
import json
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d, percentile_filter, median_filter
from collections import defaultdict

def compute_network_bursts(
    ax_raster,
    ax_macro,
    SpikeTimes,
    # ---------------- core detection params ----------------
    isi_threshold=0.1,             # seconds for unit burst detection
    bin_ms=10,                     # bin width for participation histogram (ms)
    gamma=2.0,                     # exponent for weighted participation
    smoothing_min_ms=20,           # min smoothing (ms)
    baseline_pct=20.0,             # rolling baseline percentile
    local_baseline_win_s=5.0,      # rolling baseline window (s)
    delta_thresh=0.30,             # ΔP/P threshold to start active candidate
    peak_top_percent=5.0,          # for dyn estimate
    noise_k=4.0,                   # noise floor multiplier for dyn/prom
    prom_frac=0.15,                # fraction of dyn used for prominence
    max_back_ms=500.0,             # backtracking window for onset (ms)
    slope_k=2.0,                   # slope multiplier for slope-threshold
    min_burstlet_duration_ms=40,   # minimum burstlet duration (ms)
    # ---------------- merging params ----------------
    burstlet_merge_gap_s=0.05,     # small-gap to merge burstlets -> network bursts (s)
    network_merge_gap_s=0.1,       # long-gap to merge network bursts -> superbursts (s)
    superburst_min_duration_s=5.0, # minimum duration for a superburst (s)
    # ---------------- plotting / verbosity ----------------
    plot=True,
    verbose=False
):
    """
    Full 4-level burst pipeline:
      - unit_bursts: per-unit ISI-bursts
      - burstlets: ΔP/P active windows validated by micro-peaks (with backtracking onset)
      - network_bursts: merge burstlets closer than burstlet_merge_gap_s
      - superbursts: merge network_bursts closer than network_merge_gap_s and keep those longer than superburst_min_duration_s

    Returns:
      results (dict, JSON-serializable) and also draws onto ax_raster and ax_macro if provided.
    """

    # ---------------------------
    # Sanity checks / quick returns
    # ---------------------------
    units = list(SpikeTimes.keys())
    if len(units) == 0:
        return {"error": "no_units"}, None

    # Concatenate and sort all spikes (for recording bounds)
    all_spikes = np.sort(np.concatenate([np.asarray(SpikeTimes[u]) for u in units]))
    if all_spikes.size == 0:
        return {"error": "no_spikes"}, None

    rec_start, rec_end = float(all_spikes[0]), float(all_spikes[-1])
    if verbose:
        print(f"Recording: {rec_start:.3f} → {rec_end:.3f} s, Units: {len(units)}, Total spikes: {len(all_spikes)}")

    # ---------------------------
    # 1) Unit-level ISI bursts
    # ---------------------------
    def detect_unit_bursts(spike_times_dict, isi_thresh):
        out = {}
        for unit, times in spike_times_dict.items():
            times = np.asarray(times)
            if times.size < 2:
                out[unit] = {
                    "bursts": [],
                    "isis_within": [],
                    "isis_outside": [],
                    "isis_all": []
                }
                continue
            isis = np.diff(times)
            mask = isis < isi_thresh
            starts = np.where(np.diff(np.insert(mask.astype(int), 0, 0)) == 1)[0]
            ends   = np.where(np.diff(np.append(mask.astype(int), 0)) == -1)[0]
            bursts = [list(times[s:e+1]) for s,e in zip(starts, ends)]
            isis_within = np.concatenate([isis[s:e] for s,e in zip(starts, ends)]) if len(starts)>0 else np.array([])
            isis_outside = isis[~mask] if isis.size>0 else np.array([])
            out[unit] = {
                "bursts": bursts,
                "isis_within": isis_within.tolist(),
                "isis_outside": isis_outside.tolist(),
                "isis_all": isis.tolist()
            }
        return out

    unit_bursts = detect_unit_bursts(SpikeTimes, isi_threshold)

    # ---------------------------
    # 2) Participation histogram / weighted participation
    # ---------------------------
    bin_size = float(bin_ms) / 1000.0
    bins = np.arange(rec_start, rec_end + bin_size, bin_size)
    t_centers = (bins[:-1] + bins[1:]) / 2.0

    P = np.zeros(len(t_centers), dtype=float)
    for u in units:
        counts, _ = np.histogram(SpikeTimes[u], bins=bins)
        P += (counts > 0).astype(float)

    W = (P ** gamma) if gamma != 1.0 else P.copy()

    # ---------------------------
    # 3) Two-stream smoothing (fast for onsets, slow for peaks) - recommended
    # ---------------------------
    # convert smoothing minima
    sigma_min_bins = max(1.0, smoothing_min_ms / bin_ms)

    # adapt sigma using median ISI if available
    diffs = np.diff(all_spikes)
    diffs = diffs[diffs > 0]
    median_isi = np.median(diffs) if len(diffs) > 0 else 0.001
    target_sigma_s = 5.0 * median_isi
    # convert to bins
    target_sigma_bins = float(np.clip(target_sigma_s / bin_size, 1.0, 30.0))

    # choose fast and slow sigmas (fast for onset, slow for peak)
    sigma_fast_bins = max(1.0, min(5.0, sigma_min_bins))         # small smoothing for onset
    sigma_slow_bins = max(sigma_min_bins, target_sigma_bins)     # slower for peak detection

    W_smooth_fast = gaussian_filter1d(W, sigma=sigma_fast_bins)
    W_smooth_slow = gaussian_filter1d(W, sigma=sigma_slow_bins)

    # default ws for different steps
    ws_onset = W_smooth_fast   # used for delta, onset timing, slope
    ws_peak  = W_smooth_slow   # used for peak detection and dyn/prom

    # ---------------------------
    # 4) Micro-peaks on slow trace & dyn/prom
    # ---------------------------
    ws_p = ws_peak
    n_top = max(1, int(len(ws_p) * (peak_top_percent / 100.0)))
    peak_mean = float(np.mean(np.sort(ws_p)[-n_top:]))
    baseline_global = float(np.percentile(ws_onset, baseline_pct))

    dws_peak = np.abs(np.diff(ws_p))
    noise_est = float(1.4826 * np.median(dws_peak)) if len(dws_peak) > 0 else 0.0

    dyn = max(peak_mean - baseline_global, noise_k * noise_est, 1e-9)
    prom = max(prom_frac * dyn, noise_k * noise_est, 1e-9)

    peaks_idx, peak_props = find_peaks(ws_p, prominence=prom)
    # ensure numpy ints -> python ints for json
    peaks_idx = peaks_idx.astype(int).tolist()

    # ---------------------------
    # 5) Rolling baseline (fast trace) and ΔP/P
    # ---------------------------
    win_bins = max(3, int(round(local_baseline_win_s / bin_size)))
    try:
        baseline_local = percentile_filter(ws_onset.astype(float), percentile=baseline_pct, size=win_bins, mode='reflect')
    except Exception:
        baseline_local = np.full_like(ws_onset, baseline_global)

    eps = 1e-9
    delta = (ws_onset - baseline_local) / (baseline_local + eps)

    # ---------------------------
    # 6) Active candidate windows using ΔP/P
    # ---------------------------
    above = delta > delta_thresh
    actives_idx = []
    i = 0
    n = len(above)
    while i < n:
        if above[i]:
            s = i
            while i < n and above[i]:
                i += 1
            e = i - 1
            actives_idx.append((s, e))
        else:
            i += 1

    # ---------------------------
    # 7) Validate active windows -> burstlets (micro-bursts)
    #     require micropeak inside, backtrack to onset (local min on ws_onset),
    #     apply slope check, enforce min duration
    # ---------------------------
    burstlets = []   # list of dicts
    max_back_bins = max(1, int(round((max_back_ms / 1000.0) / bin_size)))
    slope_noise = 1.4826 * np.median(np.abs(np.diff(ws_onset))) if len(ws_onset) > 1 else 0.0
    slope_thresh = slope_k * slope_noise

    for (s_idx, e_idx) in actives_idx:
        s_t = float(t_centers[s_idx]); e_t = float(t_centers[e_idx])
        dur_ms = (e_t - s_t) * 1000.0
        if dur_ms < min_burstlet_duration_ms:
            # too short -> skip
            continue

        # find peaks in slow trace that fall inside candidate
        p_inside = [p for p in peaks_idx if (p >= s_idx and p <= e_idx)]
        if len(p_inside) == 0:
            continue

        # pick primary peak (max amplitude on ws_p)
        p_vals = [ws_p[p] for p in p_inside]
        pk_idx = int(p_inside[int(np.argmax(p_vals))])
        pk_t = float(t_centers[pk_idx])

        # backtrack onset: local minimum before pk within max_back_bins
        start_search = max(0, pk_idx - max_back_bins)
        seg = ws_onset[start_search:pk_idx+1]
        if seg.size == 0:
            onset_idx = s_idx
        else:
            rel_min = int(np.argmin(seg))
            onset_idx = start_search + rel_min

        # slope check on ws_onset near onset
        post_idx = min(len(ws_onset)-1, onset_idx + 3)
        denom = ((post_idx - onset_idx + 1) * bin_size)
        avg_slope = (ws_onset[post_idx] - ws_onset[onset_idx]) / denom if denom > 0 else 0.0
        if slope_k > 0 and avg_slope < slope_thresh:
            # try to find earliest index with sufficient slope between onset_idx and pk_idx
            found = False
            for kk in range(onset_idx, pk_idx):
                pw_end = min(len(ws_onset)-1, kk + 3)
                denom2 = ((pw_end - kk + 1) * bin_size)
                avg_s = (ws_onset[pw_end] - ws_onset[kk]) / denom2 if denom2 > 0 else 0.0
                if avg_s >= slope_thresh:
                    onset_idx = kk
                    found = True
                    break
            if not found:
                # treat as valid if the peak is very prominent relative to baseline (defensive)
                if (ws_p[pk_idx] - baseline_global) < 0.5 * dyn:
                    continue

        onset_t = float(t_centers[onset_idx])
        # final duration check (from onset to e_t)
        final_dur_ms = (e_t - onset_t) * 1000.0
        if final_dur_ms < min_burstlet_duration_ms:
            continue

        # compute spike counts and participating units inside [onset_t, e_t]
        total_spikes = 0
        participating_units = 0
        for u in units:
            times = np.asarray(SpikeTimes[u])
            if times.size == 0:
                continue
            mask = (times >= onset_t) & (times <= e_t)
            cnt = int(np.sum(mask))
            total_spikes += cnt
            if cnt > 0:
                participating_units += 1

        # peak rate and energy inside
        mask_bins = (t_centers >= onset_t) & (t_centers <= e_t)
        peak_rate_inside = float(np.max(ws_p[mask_bins])) if np.any(mask_bins) else None
        energy = float(np.sum(ws_p[mask_bins]) * bin_size) if np.any(mask_bins) else 0.0

        burstlets.append({
            "start": onset_t,
            "end": e_t,
            "duration_s": float((e_t - onset_t)),
            "total_spikes": int(total_spikes),
            "participating_units": int(participating_units),
            "primary_peak_index": int(pk_idx),
            "primary_peak_time": pk_t,
            "primary_peak_amplitude": float(ws_p[pk_idx]),
            "peak_rate_inside": peak_rate_inside,
            "burst_energy": energy
        })

    # sort burstlets by start
    burstlets = sorted(burstlets, key=lambda x: x["start"])

    # ---------------------------
    # 8) Merge burstlets -> network_bursts (short-gap merging)
    # ---------------------------
    def merge_windows(windows, gap_s, min_duration_s=0.0, return_fragments=False):
        """
        windows: list of (start, end) or dicts with "start","end"
        gap_s: maximum gap to merge
        returns list of dicts with start,end,duration,inner_windows
        """
        if len(windows) == 0:
            return []

        # canonicalize to (s,e)
        se = [(w["start"], w["end"]) if isinstance(w, dict) else (float(w[0]), float(w[1])) for w in windows]
        se_sorted = sorted(se, key=lambda x: x[0])

        merged = []
        cs, ce = se_sorted[0]
        fragments = [(cs, ce)]
        for (s, e) in se_sorted[1:]:
            if s - ce <= gap_s:
                ce = max(ce, e)
                fragments.append((s, e))
            else:
                if (ce - cs) >= min_duration_s:
                    merged.append({"start": float(cs), "end": float(ce), "duration_s": float(ce - cs), "fragments": list(fragments)})
                cs, ce = s, e
                fragments = [(s, e)]
        # last
        if (ce - cs) >= min_duration_s:
            merged.append({"start": float(cs), "end": float(ce), "duration_s": float(ce - cs), "fragments": list(fragments)})
        return merged

    # prepare list-of-(s,e) for merging
    burstlet_windows = [(b["start"], b["end"]) for b in burstlets]
    network_bursts = merge_windows(burstlet_windows, gap_s=burstlet_merge_gap_s, min_duration_s=0.0)

    # enrich network bursts with stats (counts, IBIs inside, peak rate, energy, fragment count)
    enriched_network_bursts = []
    for nb in network_bursts:
        ms, me = nb["start"], nb["end"]
        # inner burstlets
        inner = [b for b in burstlets if (b["start"] >= ms and b["end"] <= me)]
        frag_count = len(inner)
        total_spikes = sum(int(b["total_spikes"]) for b in inner)
        participating_units = len(set([u for b in inner for u in units if np.any((np.asarray(SpikeTimes[u]) >= b["start"]) & (np.asarray(SpikeTimes[u]) <= b["end"]))]))
        # IBIs between burstlet starts
        if frag_count > 1:
            starts = np.array([b["start"] for b in inner])
            IBIs = np.diff(starts)
            mean_IBI_inside = float(np.mean(IBIs))
        else:
            mean_IBI_inside = None
        mask_bins = (t_centers >= ms) & (t_centers <= me)
        peak_inside = float(np.max(ws_p[mask_bins])) if np.any(mask_bins) else None
        energy = float(np.sum(ws_p[mask_bins]) * bin_size) if np.any(mask_bins) else 0.0

        enriched_network_bursts.append({
            "start": ms,
            "end": me,
            "duration_s": me - ms,
            "fragment_count": frag_count,
            "total_spikes": int(total_spikes),
            "participating_units": int(participating_units),
            "mean_IBI_inside": mean_IBI_inside,
            "peak_rate_inside": peak_inside,
            "burst_energy": energy,
            "inner_burstlets": inner
        })

    # ---------------------------
    # 9) Merge network_bursts -> superbursts (long-gap merging)
    # ---------------------------
    net_windows = [(nb["start"], nb["end"]) for nb in enriched_network_bursts]
    superbursts_raw = merge_windows(net_windows, gap_s=network_merge_gap_s, min_duration_s=superburst_min_duration_s)
    enriched_superbursts = []
    for sb in superbursts_raw:
        ms, me = sb["start"], sb["end"]
        inner_nbs = [nb for nb in enriched_network_bursts if (nb["start"] >= ms and nb["end"] <= me)]
        frag_count = len(inner_nbs)
        total_spikes = sum(nb["total_spikes"] for nb in inner_nbs)
        participating_units = len(set([u for nb in inner_nbs for u in units if np.any((np.asarray(SpikeTimes[u]) >= nb["start"]) & (np.asarray(SpikeTimes[u]) <= nb["end"]))]))
        if frag_count > 1:
            starts = np.array([nb["start"] for nb in inner_nbs])
            IBIs = np.diff(starts)
            mean_IBI_inside = float(np.mean(IBIs))
        else:
            mean_IBI_inside = None
        mask_bins = (t_centers >= ms) & (t_centers <= me)
        peak_inside = float(np.max(ws_p[mask_bins])) if np.any(mask_bins) else None
        energy = float(np.sum(ws_p[mask_bins]) * bin_size) if np.any(mask_bins) else 0.0

        enriched_superbursts.append({
            "start": ms,
            "end": me,
            "duration_s": me - ms,
            "fragment_count": frag_count,
            "total_spikes": int(total_spikes),
            "participating_units": int(participating_units),
            "mean_IBI_inside": mean_IBI_inside,
            "peak_rate_inside": peak_inside,
            "burst_energy": energy,
            "inner_network_bursts": inner_nbs
        })

    # ---------------------------
    # 10) Aggregate statistics at each level
    # ---------------------------
    def summarize_events(events, duration_key="duration_s"):
        if len(events) == 0:
            return {
                "count": 0,
                "mean_duration": None,
                "median_duration": None,
                "std_duration": None,
                "mean_participating_units": None,
                "mean_total_spikes": None
            }
        durations = np.array([e[duration_key] for e in events if e[duration_key] is not None])
        participating = np.array([e.get("participating_units", np.nan) or np.nan for e in events])
        spikes = np.array([e.get("total_spikes", np.nan) or np.nan for e in events])
        return {
            "count": int(len(events)),
            "mean_duration": float(np.nanmean(durations)) if durations.size>0 else None,
            "median_duration": float(np.nanmedian(durations)) if durations.size>0 else None,
            "std_duration": float(np.nanstd(durations)) if durations.size>0 else None,
            "mean_participating_units": float(np.nanmean(participating)) if participating.size>0 else None,
            "mean_total_spikes": float(np.nanmean(spikes)) if spikes.size>0 else None
        }

    agg_unit = {"count": len(unit_bursts)}
    agg_burstlets = summarize_events(burstlets, duration_key="duration_s")
    agg_network = summarize_events(enriched_network_bursts, duration_key="duration_s")
    agg_super = summarize_events(enriched_superbursts, duration_key="duration_s")

    # ---------------------------
    # 11) Plotting
    # ---------------------------
    if plot and ax_raster is not None:
        # raster sorted by spike count (low -> high)
        spike_counts = {u: len(SpikeTimes[u]) for u in units}
        sorted_units = sorted(spike_counts, key=spike_counts.get)
        y = 0
        for unit in sorted_units:
            times = np.asarray(SpikeTimes[unit])
            ax_raster.plot(times, np.ones_like(times) + y, '|', color='royalblue', markersize=1)
            # highlight unit bursts
            for ub in unit_bursts[unit]["bursts"]:
                ax_raster.plot(ub, np.ones_like(ub) + y, '|', color='black', markersize=1)
            y += 1
        ax_raster.set_ylabel('Units')
        ax_raster.set_title('Raster with Unit Bursts (black)')

    if plot and ax_macro is not None:
        # plot slow participation trace and decorations
        ax_macro.plot(t_centers, ws_peak, lw=1.2, label='W_smooth (peak stream)')
        ax_macro.plot(t_centers, ws_onset, lw=1.0, alpha=0.6, label='W_smooth (onset stream)')
        # micro-peaks
        ax_macro.plot([t_centers[int(p)] for p in peaks_idx], [float(ws_p[int(p)]) for p in peaks_idx], 'x', label='micro-peaks')
        # local baseline
        ax_macro.plot(t_centers, baseline_local, ls=':', lw=1.0, label=f'baseline_local (p{baseline_pct})')
        # # burstlets
        # for b in burstlets:
        #     ax_macro.axvspan(b["start"], b["end"], color='C1', alpha=0.25)
        #     ax_macro.plot([b["start"]], [b["primary_peak_amplitude"]], 'o', markersize=3, color='C1')
        # network bursts (merge of burstlets)
        for nb in enriched_network_bursts:
            ax_macro.axvspan(nb["start"], nb["end"], color='C0', alpha=0.2)
        # superbursts
        for sb in enriched_superbursts:
            #ax_macro.axvspan(sb["start"], sb["end"], color='gray', alpha=0.25)
            ax_macro.hlines(np.max(ws_p) * 1.04, sb["start"], sb["end"], colors='k', linewidth=2)
        ax_macro.set_xlabel('Time (s)')
        ax_macro.set_ylabel('Weighted participation')
        ax_macro.set_title(f'Network bursts: {len(enriched_network_bursts)} | Superbursts: {len(enriched_superbursts)}')
        ax_macro.legend(fontsize='small')

    # ---------------------------
    # 12) Package JSON-serializable results
    # ---------------------------
    def safe_numpy(obj):
        if isinstance(obj, np.generic):
            return obj.item()
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    results = {
        "recording_start": float(rec_start),
        "recording_end": float(rec_end),
        "params": {
            "isi_threshold": isi_threshold,
            "bin_ms": bin_ms,
            "gamma": gamma,
            "smoothing_min_ms": smoothing_min_ms,
            "baseline_pct": baseline_pct,
            "local_baseline_win_s": local_baseline_win_s,
            "delta_thresh": delta_thresh,
            "peak_top_percent": peak_top_percent,
            "noise_k": noise_k,
            "prom_frac": prom_frac,
            "max_back_ms": max_back_ms,
            "min_burstlet_duration_ms": min_burstlet_duration_ms,
            "burstlet_merge_gap_s": burstlet_merge_gap_s,
            "network_merge_gap_s": network_merge_gap_s,
            "superburst_min_duration_s": superburst_min_duration_s
        },
        "diagnostics": {
            "time_centers": t_centers.tolist(),
            "participation": P.tolist(),
            "weighted": W.tolist(),
            "W_smooth_onset": ws_onset.tolist(),
            "W_smooth_peak": ws_peak.tolist(),
            "baseline_local": baseline_local.tolist(),
            "delta": delta.tolist(),
            "micro_peaks_idx": [int(x) for x in peaks_idx]
        },
        "unit_bursts": unit_bursts,
        "burstlets": burstlets,
        "network_bursts": enriched_network_bursts,
        "superbursts": enriched_superbursts,
        "aggregates": {
            "unit_level": {"n_units": len(units)},
            "burstlet_level": agg_burstlets,
            "network_level": agg_network,
            "superburst_level": agg_super
        }
    }

    # convert any numpy types inside nested dicts to Python natives where necessary
    # The structure above uses only native types and lists, but ensure safety by json.dumps roundtrip
    try:
        _ = json.dumps(results)  # validate serializability (will raise if not)
    except TypeError:
        # fallback: do lightweight conversion
        def convert(o):
            if isinstance(o, dict):
                return {str(k): convert(v) for k, v in o.items()}
            if isinstance(o, (list, tuple)):
                return [convert(v) for v in o]
            if isinstance(o, np.ndarray):
                return o.tolist()
            if isinstance(o, np.generic):
                return o.item()
            return o
        results = convert(results)

    return results