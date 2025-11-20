import numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
from math import exp
from typing import Dict, Iterable


def plot_network_bursts(
    ax,
    SpikeTimes,
    min_duration_ms_floor=10.0,
    participation_frac_floor=0.15):


    # ------------------------------------------------------------
    # PREP SPIKES
    # ------------------------------------------------------------

    units = list(SpikeTimes.keys())
    if len(units) == 0:
        return ax, {"qc": "no_units"}

    all_spikes = np.sort(np.concatenate([np.asarray(SpikeTimes[u]) for u in units]))
    if all_spikes.size == 0:
        return ax, {"qc": "no_spikes"}

    # Adaptive dedupe
    def _dedupe_spikes_adaptive(spike_times):
        if spike_times.size == 0:
            return spike_times, 0.0
        spike_times = np.sort(spike_times)
        diffs = np.diff(spike_times)
        pos = diffs[diffs > 0]
        eps = 1e-6 if pos.size == 0 else max(1e-6, 0.5 * np.median(pos))
        groups = []
        cur = [spike_times[0]]
        for t in spike_times[1:]:
            if t - cur[-1] <= eps:
                cur.append(t)
            else:
                groups.append(np.mean(cur))
                cur = [t]
        groups.append(np.mean(cur))
        return np.asarray(groups), eps

    pooled, dedupe_eps = _dedupe_spikes_adaptive(all_spikes)

    info = {
        "n_units": len(units),
        "n_spikes_raw": int(all_spikes.size),
        "n_spikes_dedup": int(pooled.size),
        "dedupe_eps": float(dedupe_eps),
    }

    if pooled.size < 5:
        info["qc"] = "too_few_spikes_after_dedupe"
        return ax, info

    rec_start, rec_end = pooled[0], pooled[-1]
    rec_dur = rec_end - rec_start

    # ------------------------------------------------------------
    # TIMESCALE FROM AUTOCORR
    # ------------------------------------------------------------
    def _autotimescale_from_autocorr(spike_times, max_ms=500.0):
        if spike_times.size < 10:
            return 100.0
        rec_start, rec_end = spike_times[0], spike_times[-1]
        dur = rec_end - rec_start
        if dur <= 0:
            return 100.0

        n_bins = int(min(dur * 1000 + 1, 200000))
        bins = np.linspace(rec_start, rec_end, n_bins)
        counts, _ = np.histogram(spike_times, bins=bins)

        max_lag = int(min(max_ms, n_bins - 2))
        ac = np.correlate(counts - np.mean(counts),
                          counts - np.mean(counts),
                          mode='full')
        c = len(ac) // 2
        ac = ac[c:c + max_lag + 1]
        if np.all(ac == 0):
            return 100.0

        ac = ac / (ac[0] + 1e-12)
        below = np.where(ac < np.exp(-1))[0]
        if below.size > 1:
            tau = below[1]
        elif below.size == 1:
            tau = below[0]
        else:
            pks, _ = find_peaks(ac)
            tau = pks[0] if len(pks) else min(100, max_lag)

        return float(np.clip(tau, 20.0, 500.0))

    timescale_ms = _autotimescale_from_autocorr(pooled)
    info["timescale_ms"] = timescale_ms

    # ------------------------------------------------------------
    # BUILD POPULATION RATE
    # ------------------------------------------------------------
    bin_ms = 10.0
    bin_size = bin_ms / 1000.0
    bins = np.arange(rec_start, rec_end + bin_size, bin_size)
    counts, _ = np.histogram(pooled, bins=bins)
    t_centers = bins[:-1] + bin_size / 2
    rate = counts / bin_size

    # SMOOTHING
    sigma_from_timescale = max(1.0, timescale_ms / bin_ms)
    rate_std = np.std(rate)
    rate_med = np.median(rate) + 1e-12
    noise_ratio = rate_std / rate_med
    sigma_min = np.clip(5 + 10 * noise_ratio, 5, 30)
    sigma_bins = max(sigma_from_timescale, sigma_min)
    rate_smooth = gaussian_filter1d(rate.astype(float), sigma=sigma_bins)

    # ------------------------------------------------------------
    # PEAK / VALLEY DETECTION
    # ------------------------------------------------------------
    dyn = np.max(rate_smooth) - np.median(rate_smooth)
    prominence = max(0.15 * dyn, 1e-6)
    peaks, _ = find_peaks(rate_smooth, prominence=prominence)
    valley_prominence = max(0.08 * dyn, 1e-6)   # you can tune 0.08 → 0.10–0.15
    valleys, valley_props = find_peaks(-rate_smooth, prominence=valley_prominence)
    valleys= np.sort(valleys)

    if len(peaks) == 0 or len(valleys) < 2:
        info["qc"] = "no_peaks_or_valleys"
        ax.plot(t_centers, rate_smooth)
        return ax, info

    # ------------------------------------------------------------
    # CANDIDATE WINDOWS
    # ------------------------------------------------------------
    candidates = []
    for pk in peaks:
        L = valleys[valleys < pk]
        R = valleys[valleys > pk]
        if len(L) == 0 or len(R) == 0:
            continue
        vL = L[-1]
        vR = R[0]
        s = t_centers[vL]
        e = t_centers[vR]
        if e <= s:
            continue
        nspk = int(np.sum((pooled >= s) & (pooled <= e)))
        peakval = rate_smooth[pk]
        valleyval = min(rate_smooth[vL], rate_smooth[vR])
        vfrac = valleyval / (peakval + 1e-12)
        candidates.append({
            "s": s, "e": e, "n_spk": nspk,
            "pk_idx": pk,
            "peak_val": peakval,
            "valley_frac": vfrac
        })

    if len(candidates) == 0:
        info["qc"] = "no_candidate_windows"
        return ax, info

    # ------------------------------------------------------------
    # VALIDATION (PARTICIPATION + PVAL)
    # ------------------------------------------------------------
    baseline_rate = float(np.percentile(rate_smooth, 25) + 1e-12)
    unit_spikes = {u: np.sort(np.asarray(SpikeTimes[u])) for u in units}
    U = len(units)

    kept = []
    for w in candidates:
        s, e = w["s"], w["e"]
        dur = max(1e-6, e - s)

        # participation
        pcount = 0
        for sps in unit_spikes.values():
            if sps.size == 0:
                continue
            i0 = np.searchsorted(sps, s)
            if i0 < sps.size and sps[i0] <= e:
                pcount += 1

        w["part_frac"] = pcount / U

        # Poisson tail
        expected = baseline_rate * dur

        # compute poisson tail p-value
        def _poisson_tail_p_value(k, lam):
            if lam <= 0:
                return 0.0 if k > 0 else 1.0
            term = np.exp(-lam)
            cdf = term
            for i in range(1, k):
                term = term * lam / i
                cdf += term
            return max(1e-300, 1.0 - cdf)

        w["pval"] = _poisson_tail_p_value(w["n_spk"], expected)

    # adaptive thresholds
    part_list = [w["part_frac"] for w in candidates]
    adaptive_part_floor = max(participation_frac_floor,
                              0.8 * np.median(part_list))
    alpha = 0.01 / max(1, len(candidates))

    # select final bursts
    for w in candidates:
        if (w["e"] - w["s"]) * 1000 < min_duration_ms_floor:
            continue
        if w["part_frac"] >= adaptive_part_floor or w["pval"] < alpha:
            kept.append(w)

    if len(kept) == 0:
        info["qc"] = "no_windows_passed_validation"
        return ax, info

    # ------------------------------------------------------------
    # COMPUTE 90% MASS BURST CORES
    # ------------------------------------------------------------

    burst_list = []
    core_list = []   # (core_start, core_end)

    for w in kept:
        s, e = w["s"], w["e"]
        sp = pooled[(pooled >= s) & (pooled <= e)]

        if len(sp) < 5:
            # trivial burst
            burst_list.append((s, e, w["n_spk"]))
            core_list.append((s, e))
            continue

        N = len(sp)
        K = int(0.9 * N)

        widths = sp[K:] - sp[:N - K]
        min_idx = np.argmin(widths)
        core_start = sp[min_idx]
        core_end = sp[min_idx + K]
        burst_list.append((core_start, core_end, w["n_spk"]))
        core_list.append((core_start, core_end))

    # ------------------------------------------------------------
    # PARTICIPATION + METRICS WITH 90% DURATION
    # ------------------------------------------------------------

    participation = {}
    for u, sps in unit_spikes.items():
        arr = []
        for (cs, ce, _) in burst_list:
            if sps.size == 0:
                arr.append(0)
                continue
            i0 = np.searchsorted(sps, cs)
            arr.append(int(i0 < sps.size and sps[i0] <= ce))
        participation[u] = arr

    burst_metrics = []
    for idx, (cs, ce, n_spk) in enumerate(burst_list):
        dur = ce - cs
        mask = (t_centers >= cs) & (t_centers <= ce)
        peak_rate = np.max(rate_smooth[mask]) if np.any(mask) else np.nan
        parts = [participation[u][idx] for u in participation]
        participation_fraction = np.mean(parts)

        burst_metrics.append({
            "burst_id": idx,
            "start": cs, "end": ce,
            "n_spikes": n_spk,
            "peak_rate": peak_rate,
            "participation_fraction": participation_fraction,
            "duration": dur})
    
    info["n_bursts"] = len(burst_metrics)
    # ------------------------------------------------------------
    #  MEGA-BURST MERGING  (coarse UP-state detection)
    # ------------------------------------------------------------

    MEGA_GAP = 0.2   # 200 ms threshold for UP-state gap

    # burst_list = [(core_start, core_end, n_spikes)]
    burst_list_sorted = sorted(burst_list, key=lambda x: x[0])

    megabursts = []
    cur_s, cur_e, cur_n = burst_list_sorted[0]
    cur_fragments = 1     # count microbursts merged

    for (s, e, n) in burst_list_sorted[1:]:
        gap = s - cur_e
        if gap < MEGA_GAP:
            # extend megaburst
            cur_e = max(cur_e, e)
            cur_n += n
            cur_fragments += 1
        else:
            megabursts.append((cur_s, cur_e, cur_n, cur_fragments))
            cur_s, cur_e, cur_n = s, e, n
            cur_fragments = 1

    megabursts.append((cur_s, cur_e, cur_n, cur_fragments))

    # Minimum duration threshold for a biologically valid mega-burst
    MEGA_MIN_DURATION = 2.0   # seconds, based on MEA superburst literature

    # ------------------------------------------------------------
    #  ENHANCED MEGABURST METRICS
    # ------------------------------------------------------------

    megabursts_out = []

    for (ms, me, mn, frag_count) in megabursts:

        duration = me - ms
        if duration < MEGA_MIN_DURATION:
            continue

        # microbursts inside this megaburst
        inner = [(cs, ce, n)
                 for (cs, ce, n) in burst_list_sorted
                 if cs >= ms and ce <= me]

        # --- A) IBIs inside megaburst ---
        if len(inner) > 1:
            IBIs = np.diff([b[0] for b in inner])
            mean_IBI_inside = float(np.mean(IBIs))
        else:
            mean_IBI_inside = np.nan

        # --- B) Peak rate inside megaburst ---
        mask = (t_centers >= ms) & (t_centers <= me)
        if np.any(mask):
            peak_rate_inside = float(np.max(rate_smooth[mask]))
        else:
            peak_rate_inside = np.nan

        # --- C) Burst energy (integrated rate) ---
        dt = bin_size
        if np.any(mask):
            burst_energy = float(np.sum(rate_smooth[mask]) * dt)
        else:
            burst_energy = 0.0

        # --- assemble dictionary ---
        megabursts_out.append({
            "start": float(ms),
            "end": float(me),
            "duration": float(duration),
            "total_spikes": int(mn),
            "fragment_count": int(frag_count),
            "mean_IBI_inside": mean_IBI_inside,
            "peak_rate_inside": peak_rate_inside,
            "burst_energy": burst_energy
        })

    # store enhanced megabursts
    info["megaburst_count"] = len(megabursts_out)
    info["megabursts"] = megabursts_out



    # ------------------------------------------------------------
    # PLOT
    # ------------------------------------------------------------

    ax.plot(t_centers, rate_smooth, '-', lw=1.2, label="Network rate (Hz)")

    # accepted peaks only
    kept_peak_idx = [w["pk_idx"] for w in kept]
    ax.plot(t_centers[kept_peak_idx],
            rate_smooth[kept_peak_idx], 'x', label="peaks")

    # valleys for reference
    ax.plot(t_centers[valleys], rate_smooth[valleys], 'o', ms=4, label="valleys")

    # SHADE 90% MASS CORES
    for (cs, ce) in core_list:
        ax.axvspan(cs, ce, color='gray', alpha=0.25)

        # thin bar on top for publication
        ymax = np.max(rate_smooth)
        ax.hlines(ymax * 1.02, cs, ce,
                  colors='blue', linewidth=2)
        
    for megaburst in megabursts_out:
        ax.hlines(ymax * 1.05, megaburst["start"], megaburst["end"], colors='black', linewidth=2)

    ax.set_xlim(rec_start, rec_end)
    ax.set_ylabel("Rate (Hz)")
    ax.set_xlabel("Time (s)")
    ax.set_title(f"Unified Parameter-free Bursts (N={len(burst_list)})")
    ax.legend()

    return ax, {
        "network_bursts": burst_list,
        "participation": participation,
        "burst_metrics": burst_metrics,
        "info": info,
        "rate_smooth": rate_smooth
    }