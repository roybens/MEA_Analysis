import numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
from math import exp


# ------------------------------------------------------------
# UNIT BURSTINESS METRIC (CV2)
# ------------------------------------------------------------
def compute_cv2(spike_train):
    """
    Computes CV2: 2|ISI(i+1)-ISI(i)| / (ISI(i+1)+ISI(i)).
    Good for discriminating tonic vs bursty units.
    """
    if len(spike_train) < 3:
        return np.nan
    isi = np.diff(spike_train)
    if len(isi) < 2:
        return np.nan
    return np.mean(2 * np.abs(np.diff(isi)) / (isi[:-1] + isi[1:] + 1e-12))


def plot_network_bursts(
    ax,
    SpikeTimes,
    min_duration_ms_floor=10.0,
    participation_frac_floor=0.15):


    # ------------------------------------------------------------
    # STEP 0: REMOVE TONIC / NON-BURSTY UNITS
    # ------------------------------------------------------------
    units = list(SpikeTimes.keys())

    unit_cv2 = {}
    for u in units:
        unit_cv2[u] = compute_cv2(np.asarray(SpikeTimes[u]))

    # CV2 < 0.7 -> strongly regular -> tonic -> suppress
    tonic_units = [u for u,cv in unit_cv2.items()
                   if not np.isnan(cv) and cv < 0.7]

    print(f"[Tonic filter] Removing {len(tonic_units)} tonic units:", tonic_units[:10])

    good_units = [u for u in units if u not in tonic_units]

    # fail-safe fallback
    if len(good_units) < 5:
        print("[Tonic filter] Too aggressive → using all units.")
        good_units = units

    # prune SpikeTimes
    SpikeTimes = {u: SpikeTimes[u] for u in good_units}
    units = good_units

    # ------------------------------------------------------------
    # PREP SPIKES
    # ------------------------------------------------------------
    if len(units) == 0:
        return ax, {"qc": "no_units"}

    all_spikes = np.sort(np.concatenate([np.asarray(SpikeTimes[u]) for u in units]))
    if all_spikes.size == 0:
        return ax, {"qc": "no_spikes"}

    # Adaptive dedupe ------------------------------------------
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
        "tonic_units_removed": len(tonic_units)
    }

    if pooled.size < 5:
        info["qc"] = "too_few_spikes_after_dedupe"
        return ax, info

    rec_start, rec_end = pooled[0], pooled[-1]


    # ------------------------------------------------------------
    # MEDIAN-ISI TIMESCALE
    # ------------------------------------------------------------
    ISI = np.diff(pooled)
    ISI = ISI[ISI > 0]
    median_isi = np.median(ISI) if ISI.size > 0 else 0.005
    target_timescale = 5 * median_isi


    # ------------------------------------------------------------
    # BUILD POPULATION RATE
    # ------------------------------------------------------------
    bin_ms = 10
    bin_size = bin_ms / 1000.0
    bins = np.arange(rec_start, rec_end + bin_size, bin_size)
    counts, _ = np.histogram(pooled, bins=bins)
    t_centers = bins[:-1] + bin_size / 2
    rate = counts / bin_size


    # ------------------------------------------------------------
    # ADAPTIVE SMOOTHING (still parameter-free)
    # ------------------------------------------------------------
    sigma_bins = target_timescale / bin_size
    sigma_bins = float(np.clip(sigma_bins, 2, 10))

    print(f"[Adaptive σ] median ISI = {median_isi*1000:.2f} ms → σ = {sigma_bins:.2f} bins")

    rate_smooth = gaussian_filter1d(rate.astype(float), sigma=sigma_bins)


    # ------------------------------------------------------------
    # PEAKS / VALLEYS
    # ------------------------------------------------------------
    dyn = np.max(rate_smooth) - np.median(rate_smooth)
    prominence = max(0.15 * dyn, 1e-6)
    peaks, _ = find_peaks(rate_smooth, prominence=prominence)

    valley_prominence = max(0.08 * dyn, 1e-6)
    valleys, _ = find_peaks(-rate_smooth, prominence=valley_prominence)
    valleys = np.sort(valleys)

    if len(peaks) == 0 or len(valleys) < 2:
        info["qc"] = "no_peaks_or_valleys"
        ax.plot(t_centers, rate_smooth)
        return ax, info


    # ------------------------------------------------------------
    # CANDIDATE WINDOWS VIA PK–VALLEY PAIRS
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
            "s": s,
            "e": e,
            "n_spk": nspk,
            "pk_idx": pk,
            "peak_val": peakval,
            "valley_frac": vfrac
        })

    if len(candidates) == 0:
        info["qc"] = "no_candidate_windows"
        return ax, info


    # ------------------------------------------------------------
    # VALIDATION: participation & Poisson p-value
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
            idx = np.searchsorted(sps, s)
            if idx < sps.size and sps[idx] <= e:
                pcount += 1
        w["part_frac"] = pcount / U

        # poisson tail
        expected = baseline_rate * dur

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

    # adaptive participation
    part_list = [w["part_frac"] for w in candidates]
    adaptive_part_floor = max(participation_frac_floor,
                              0.8 * np.median(part_list))
    alpha = 0.01 / max(1, len(candidates))

    for w in candidates:
        if (w["e"] - w["s"]) * 1000 < min_duration_ms_floor:
            continue
        if w["part_frac"] >= adaptive_part_floor or w["pval"] < alpha:
            kept.append(w)

    if len(kept) == 0:
        info["qc"] = "no_windows_passed_validation"
        return ax, info


    # ------------------------------------------------------------
    # 90% MASS CORES
    # ------------------------------------------------------------
    burst_list = []
    core_list = []

    for w in kept:
        s, e = w["s"], w["e"]
        sp = pooled[(pooled >= s) & (pooled <= e)]

        if len(sp) < 5:
            burst_list.append((s, e, w["n_spk"]))
            core_list.append((s, e))
            continue

        N = len(sp)
        K = int(0.9 * N)
        widths = sp[K:] - sp[:N - K]
        i = np.argmin(widths)
        core_start = sp[i]
        core_end = sp[i + K]

        burst_list.append((core_start, core_end, w["n_spk"]))
        core_list.append((core_start, core_end))


    # ------------------------------------------------------------
    # PARTICIPATION METRICS
    # ------------------------------------------------------------
    participation = {}
    for u, sps in unit_spikes.items():
        arr = []
        for (cs, ce, _) in burst_list:
            idx = np.searchsorted(sps, cs)
            arr.append(int(idx < sps.size and sps[idx] <= ce))
        participation[u] = arr

    burst_metrics = []
    for j, (cs, ce, n_spk) in enumerate(burst_list):
        dur = ce - cs
        mask = (t_centers >= cs) & (t_centers <= ce)
        peak_rate = np.max(rate_smooth[mask]) if np.any(mask) else np.nan
        parts = [participation[u][j] for u in participation]
        pf = np.mean(parts)

        burst_metrics.append({
            "burst_id": j,
            "start": cs,
            "end": ce,
            "duration": dur,
            "peak_rate": peak_rate,
            "n_spikes": n_spk,
            "participation_fraction": pf,
        })


    info["n_bursts"] = len(burst_metrics)


    # ------------------------------------------------------------
    # MEGABURSTS
    # ------------------------------------------------------------
    MEGA_GAP = 0.2
    burst_sorted = sorted(burst_list, key=lambda x: x[0])

    megabursts = []
    cur_s, cur_e, cur_n = burst_sorted[0]
    cur_frag = 1

    for (s, e, n) in burst_sorted[1:]:
        if s - cur_e < MEGA_GAP:
            cur_e = max(cur_e, e)
            cur_n += n
            cur_frag += 1
        else:
            megabursts.append((cur_s, cur_e, cur_n, cur_frag))
            cur_s, cur_e, cur_n, cur_frag = s, e, n, 1

    megabursts.append((cur_s, cur_e, cur_n, cur_frag))

    MEGA_MIN_DURATION = 2.0
    megabursts_out = []

    for (ms, me, mn, frag) in megabursts:
        dur = me - ms
        if dur < MEGA_MIN_DURATION:
            continue

        inner = [(cs, ce, n) for (cs, ce, n) in burst_sorted
                 if cs >= ms and ce <= me]

        if len(inner) > 1:
            IBIs = np.diff([b[0] for b in inner])
            mean_IBI = float(np.mean(IBIs))
        else:
            mean_IBI = np.nan

        mask = (t_centers >= ms) & (t_centers <= me)
        peak_inside = float(np.max(rate_smooth[mask])) if np.any(mask) else np.nan
        energy = float(np.sum(rate_smooth[mask]) * bin_size) if np.any(mask) else 0.0

        megabursts_out.append({
            "start": float(ms),
            "end": float(me),
            "duration": float(dur),
            "total_spikes": int(mn),
            "fragment_count": int(frag),
            "mean_IBI_inside": mean_IBI,
            "peak_rate_inside": peak_inside,
            "burst_energy": energy,
        })

    info["megabursts"] = megabursts_out
    info["megaburst_count"] = len(megabursts_out)


    # ------------------------------------------------------------
    # PLOTTING
    # ------------------------------------------------------------
    ax.plot(t_centers, rate_smooth, '-', lw=1.2, label="Network rate (Hz)")

    kept_peak_idx = [w["pk_idx"] for w in kept]
    ax.plot(t_centers[kept_peak_idx],
            rate_smooth[kept_peak_idx], 'x', label="peaks")

    ax.plot(t_centers[valleys], rate_smooth[valleys], 'o', ms=4, label="valleys")

    ymax = np.max(rate_smooth)

    for (cs, ce) in core_list:
        ax.axvspan(cs, ce, color='gray', alpha=0.25)
        ax.hlines(ymax * 1.02, cs, ce, colors='blue', linewidth=2)

    for meg in megabursts_out:
        ax.hlines(ymax * 1.05, meg["start"], meg["end"], colors='black', linewidth=2)

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