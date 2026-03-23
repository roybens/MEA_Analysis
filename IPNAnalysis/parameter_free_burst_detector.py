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
    # 0. Sanity checks
    # ---------------------------------------------------------
    units = list(SpikeTimes.keys())
    if not units:
        return {"error": "no_units"}

    all_spikes = np.sort(np.concatenate([SpikeTimes[u] for u in units if len(SpikeTimes[u]) > 0]))
    if all_spikes.size == 0:
        return {"error": "no_spikes"}

    rec_start = float(all_spikes[0])
    rec_end = float(all_spikes[-1])
    total_dur = rec_end - rec_start

    # ---------------------------------------------------------
    # 1. Biological calibration
    # ---------------------------------------------------------
    all_log_isis = []

    for u in units:
        t = np.unique(np.sort(SpikeTimes[u]))
        if len(t) >= 2:
            isi = np.diff(t)
            isi = isi[isi > 0]
            if isi.size > 0:
                all_log_isis.extend(np.log10(isi))

    biological_isi_s = 10 ** np.median(all_log_isis) if all_log_isis else 0.1

    adaptive_bin_ms = np.clip(biological_isi_s * 1000, 10, 30)
    bin_size = adaptive_bin_ms / 1000.0

    bins = np.arange(rec_start, rec_end + bin_size, bin_size)
    t_centers = (bins[:-1] + bins[1:]) / 2

    # ---------------------------------------------------------
    # 2. Population signals
    # ---------------------------------------------------------
    n_bins = len(t_centers)
    n_units = len(units)

    active_unit_counts = np.zeros(n_bins)
    spike_counts_total = np.zeros(n_bins)

    for u in units:

        spk = np.asarray(SpikeTimes[u])

        if spk.size == 0:
            continue

        counts, _ = np.histogram(spk, bins=bins)

        active_unit_counts += (counts > 0)
        spike_counts_total += counts

    participation_signal_raw = active_unit_counts / max(1, n_units)

    rate_signal_raw = spike_counts_total / bin_size / max(1, n_units)

    # compatibility
    P = active_unit_counts
    PFR = spike_counts_total / bin_size

    # ---------------------------------------------------------
    # 3. Smoothing
    # ---------------------------------------------------------
    isi_bins = biological_isi_s / bin_size

    sigma_fast = np.clip(1.0 * isi_bins, 1, 2)
    sigma_slow = np.clip(5.0 * isi_bins, 3, 8)

    participation_signal_smooth = gaussian_filter1d(participation_signal_raw, sigma_fast)
    rate_signal_smooth = gaussian_filter1d(rate_signal_raw, sigma_slow)

    ws_sharp = participation_signal_smooth
    ws_smooth = rate_signal_smooth

    # adaptive merge gaps
    burstlet_merge_gap_s = 3 * biological_isi_s
    network_merge_gap_s = 10 * biological_isi_s

    # ---------------------------------------------------------
    # 4. Detection thresholds
    # ---------------------------------------------------------
    participation_floor_count = max(5, 0.15 * n_units) if n_units < 50 else max(10, 0.05 * n_units)
    participation_floor = participation_floor_count / max(1, n_units)

    baseline_val = np.median(ws_sharp)
    spread_mad = np.median(np.abs(ws_sharp - baseline_val))

    relative_threshold_val = max(participation_floor, baseline_val + 0.5 * spread_mad)

    # ---------------------------------------------------------
    # 5. Burstlet detection
    # ---------------------------------------------------------
    min_peak_height = baseline_val + 0.5 * spread_mad
    peaks, _ = find_peaks(ws_sharp, height=min_peak_height)

    burstlets = []

    for p in peaks:

        s = p
        e = p

        while s > 1 and ws_sharp[s - 1] <= ws_sharp[s]:
            s -= 1

        while e < n_bins - 2 and ws_sharp[e + 1] <= ws_sharp[e]:
            e += 1

        start_idx = s
        end_idx = min(e, n_bins - 1)

        start_t = t_centers[start_idx]
        end_t = t_centers[end_idx]

        duration_s = end_t - start_t

        if duration_s <= 0:
            continue

        participating = sum(
            1 for u in units
            if np.any((SpikeTimes[u] >= start_t) & (SpikeTimes[u] <= end_t))
        )

        participation_frac = participating / n_units

        if ws_sharp[p] < relative_threshold_val:
            continue

        if participation_frac < min_burstlet_participation:
            continue

        total_spikes = int(np.sum(spike_counts_total[start_idx:end_idx + 1]))

        denom = duration_s * max(1, participating)
        burst_density = total_spikes / denom if denom > 0 else 0

        peak_drive_rate = np.max(rate_signal_raw[start_idx:end_idx + 1])

        if burst_density < min_burst_density_Hz:
            continue

        if peak_drive_rate < min_absolute_rate_Hz:
            continue

        burstlets.append({
            "start": float(start_t),
            "end": float(end_t),
            "duration_s": float(duration_s),
            "peak_synchrony": float(ws_sharp[p]),
            "peak_time": float(t_centers[p]),
            "synchrony_energy": float(np.sum(ws_smooth[start_idx:end_idx + 1]) * bin_size),
            "participation": participation_frac,
            "total_spikes": total_spikes,
            "burst_peak": float(np.max(PFR[start_idx:end_idx + 1]))
        })

    # ---------------------------------------------------------
    # 6. Merge logic
    # ---------------------------------------------------------

    merge_floor = baseline_val + 0.5 * spread_mad
    drive_floor = np.median(ws_smooth)
    max_valley_duration = 2 * biological_isi_s

    def finalize(evs, s, e):

        best = max(evs, key=lambda x: x["peak_synchrony"])

        participating_units = sum(
            1 for u in units
            if np.any((SpikeTimes[u] >= s) & (SpikeTimes[u] <= e))
        )

        participation_frac = participating_units / n_units

        return {
            "start": s,
            "end": e,
            "duration_s": e - s,
            "peak_synchrony": best["peak_synchrony"],
            "peak_time": best["peak_time"],
            "synchrony_energy": sum(ev["synchrony_energy"] for ev in evs),
            "fragment_count": len(evs),
            "total_spikes": sum(ev["total_spikes"] for ev in evs),
            "participation": participation_frac,
            "burst_peak": max(ev["burst_peak"] for ev in evs)
        }

    def merge(events, gap, floor_val, min_dur=0):

        if not events:
            return []

        events = sorted(events, key=lambda x: x["start"])

        merged = []
        curr = [events[0]]

        s = events[0]["start"]
        e = events[0]["end"]

        for nxt in events[1:]:

            mid_start_idx = np.searchsorted(t_centers, e)
            mid_end_idx = np.searchsorted(t_centers, nxt["start"])

            mid_start_idx = int(np.clip(mid_start_idx, 0, len(ws_sharp) - 1))
            mid_end_idx = int(np.clip(mid_end_idx, 0, len(ws_sharp)))

            if mid_end_idx > mid_start_idx:
                valley_part = np.min(ws_sharp[mid_start_idx:mid_end_idx])
                valley_drive = np.min(ws_smooth[mid_start_idx:mid_end_idx])
            else:
                valley_part = ws_sharp[mid_start_idx]
                valley_drive = ws_smooth[mid_start_idx]

            valley_duration = nxt["start"] - e
            peak_ratio = nxt["peak_synchrony"] / max(curr[-1]["peak_synchrony"], 1e-9)

            merge_condition = (
                (nxt["start"] <= e + gap)
                and (valley_part >= floor_val)
                and (valley_drive >= drive_floor)
                and (valley_duration <= max_valley_duration)
                and (peak_ratio > 0.2)
            )

            if merge_condition:
                curr.append(nxt)
                e = max(e, nxt["end"])
            else:
                merged.append(finalize(curr, s, e))
                curr = [nxt]
                s = nxt["start"]
                e = nxt["end"]

        merged.append(finalize(curr, s, e))

        return [m for m in merged if m["duration_s"] >= min_dur]

    network_bursts = merge(burstlets, burstlet_merge_gap_s, merge_floor)

    superbursts = merge(
        network_bursts,
        network_merge_gap_s,
        merge_floor,
        superburst_min_duration_s
    )

    # ---------------------------------------------------------
    # 7. Metrics
    # ---------------------------------------------------------
    def stats(x):

        x = np.asarray(x)

        if x.size < 2:
            return {"mean": float(x.mean()) if x.size else 0.0, "std": 0.0, "cv": 0.0}

        return {
            "mean": float(x.mean()),
            "std": float(x.std()),
            "cv": float(x.std() / x.mean())
        }

    def level_metrics(events):

        if not events:
            return {}

        starts = [ev["start"] for ev in events]

        return {
            "count": len(events),
            "rate": len(events) / total_dur,
            "duration": stats([ev["duration_s"] for ev in events]),
            "inter_event_interval": stats(np.diff(starts)) if len(starts) > 1 else stats([]),
            "intensity": stats([ev["synchrony_energy"] for ev in events]),
            "participation": stats([ev["participation"] for ev in events]),
            "spikes_per_burst": stats([ev["total_spikes"] for ev in events]),
            "burst_peak": stats([ev["burst_peak"] for ev in events]),
            "peak_synchrony": stats([ev["peak_synchrony"] for ev in events])
        }

    # ---------------------------------------------------------
    # 8. Plot
    # ---------------------------------------------------------
    if plot and ax_macro is not None:

        ax_macro.plot(t_centers, ws_smooth, lw=2)
        ax_macro.plot(t_centers, ws_sharp, lw=1)

        ax_macro.axhline(relative_threshold_val, ls="--")

        for b in network_bursts:
            ax_macro.axvspan(b["start"], b["end"], alpha=0.25)

        for s in superbursts:
            ax_macro.axvspan(s["start"], s["end"], alpha=0.3)

    # ---------------------------------------------------------
    # 9. Return
    # ---------------------------------------------------------
    return {

        "burstlets": {"events": burstlets, "metrics": level_metrics(burstlets)},

        "network_bursts": {"events": network_bursts, "metrics": level_metrics(network_bursts)},

        "superbursts": {"events": superbursts, "metrics": level_metrics(superbursts)},

        "diagnostics": {
            "adaptive_bin_ms": adaptive_bin_ms,
            "biological_isi_s": biological_isi_s,
            "baseline_value": baseline_val,
            "spread_mad": spread_mad,
            "merge_floor": merge_floor,
            "burstlet_merge_gap_s": burstlet_merge_gap_s,
            "network_merge_gap_s": network_merge_gap_s,
            "n_units": n_units,
            "sigma_fast_bins": sigma_fast,
            "sigma_slow_bins": sigma_slow
        },

        "plot_data": {
            "t": t_centers,
            "participation_signal": ws_sharp,          # dimensionless, fraction-like
            "rate_signal": ws_smooth,                  # Hz / unit
            "burst_peak_times": np.array([b["peak_time"] for b in network_bursts]),
            "burst_peak_values": np.array([b["peak_synchrony"] for b in network_bursts]),
            "participation_baseline": baseline_val,
            "participation_threshold": relative_threshold_val
        }
    }