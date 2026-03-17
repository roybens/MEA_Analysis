from __future__ import annotations

import configparser
import json
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

from ...h5_utils import (
    ensure_maxwell_hdf5_plugin_path,
    list_stream_recording_names,
    read_well_rec_frame_nos_and_trigger_settings,
)


@dataclass(frozen=True)
class RawPreprocessPlan:
    h5_path: Path
    stream_id: str
    cfg_files: tuple[Path, ...]


def discover_cfg_files(h5_path: Path) -> list[Path]:
    return sorted(Path(h5_path).parent.glob("*.cfg"))


def parse_cfg_channel_locations(cfg_path: Path) -> dict:
    cfg_path = Path(cfg_path)
    raw = cfg_path.read_text(errors="replace")
    parser = configparser.ConfigParser()
    sections: dict[str, dict[str, str]] = {}
    try:
        parser.read_string(raw)
        for section in parser.sections():
            sections[section] = dict(parser.items(section))
    except configparser.Error:
        sections = {}
    return {"path": str(cfg_path), "sections": sections, "raw": raw}


def build_preprocess_plan(*, h5_path: Path, stream_id: str, cfg_files: Optional[Iterable[Path]] = None) -> RawPreprocessPlan:
    h5_path = Path(h5_path).expanduser().resolve()
    if cfg_files is None:
        cfg_files = discover_cfg_files(h5_path)
    return RawPreprocessPlan(
        h5_path=h5_path,
        stream_id=str(stream_id),
        cfg_files=tuple(Path(p).expanduser().resolve() for p in cfg_files),
    )


def _extract_xy_from_contact_vector(contact_vector) -> tuple[list[float], list[float]]:
    x = [float(v) for v in contact_vector["x"]]
    y = [float(v) for v in contact_vector["y"]]
    return x, y


def _safe_artifact_token(token: object) -> str:
    text = str(token)
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)


def _trace_matrix_for_plot(recording, *, max_channels: int = 32, max_samples: int = 2_000):
    import numpy as np

    try:
        n_frames = int(recording.get_num_samples(segment_index=0))
    except Exception:
        n_frames = int(recording.get_num_samples())
    if n_frames <= 0:
        return np.empty((0, 0), dtype=float), np.asarray([], dtype=object)

    end_frame = min(int(max_samples), n_frames)
    channel_ids = np.asarray(recording.get_channel_ids(), dtype=object)
    if channel_ids.size == 0:
        return np.empty((0, 0), dtype=float), channel_ids
    selected_ids = channel_ids[: min(int(max_channels), int(channel_ids.size))]

    traces = recording.get_traces(
        segment_index=0,
        start_frame=0,
        end_frame=int(end_frame),
        channel_ids=list(selected_ids),
    )
    traces = np.asarray(traces, dtype=float)
    if traces.ndim != 2:
        return np.empty((0, 0), dtype=float), selected_ids
    return traces, selected_ids


def _choose_plot_channels(*, channel_ids, n_trace_channels: int, rng_seed: int = 0):
    import numpy as np

    ids = np.asarray(channel_ids, dtype=object)
    if ids.size == 0:
        return ids
    n_select = min(max(1, int(n_trace_channels)), int(ids.size))
    rng = np.random.default_rng(int(rng_seed))
    idx = rng.choice(int(ids.size), size=int(n_select), replace=False)
    idx = np.asarray(idx, dtype=int)
    return ids[idx]


def _extract_selected_traces(
    *,
    recording,
    selected_channel_ids,
    max_samples: Optional[int],
):
    import numpy as np

    available = np.asarray(recording.get_channel_ids(), dtype=object)
    selected = [ch for ch in selected_channel_ids if ch in set(available.tolist())]
    if not selected:
        selected = list(available[: min(1, int(available.size))])
    if not selected:
        return np.empty((0, 0), dtype=float), np.asarray([], dtype=object)

    total_samples = int(recording.get_num_samples())
    if max_samples is None:
        n_samples = total_samples
    else:
        n_samples = min(int(max_samples), total_samples)
    if n_samples <= 0:
        return np.empty((0, 0), dtype=float), np.asarray([], dtype=object)

    traces = recording.get_traces(
        segment_index=0,
        start_frame=0,
        end_frame=int(n_samples),
        channel_ids=list(selected),
    )
    traces = np.asarray(traces, dtype=float)
    if traces.ndim != 2:
        return np.empty((0, 0), dtype=float), np.asarray([], dtype=object)
    return traces, np.asarray(selected, dtype=object)


def _save_channel_layout_plots(
    *,
    h5_path: Path,
    stream_id: str,
    rec_names: list[str],
    common_electrodes: list[int],
    hdmea_geometry: Optional[dict[str, object]],
    out_dir: Path,
    layout_dir: Optional[Path] = None,
    logger: logging.Logger,
) -> None:
    log = logger.getChild("_save_channel_layout_plots")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        import spikeinterface.extractors as se
    except Exception as e:
        log.warning("Skipping channel layout plots: plotting dependencies unavailable (%s)", e)
        return

    _ensure_maxwell_hdf5_plugin_path()
    layout_dir = Path(layout_dir) if layout_dir is not None else (Path(out_dir) / "channel_layouts")
    layout_dir.mkdir(parents=True, exist_ok=True)
    log.info("Generating channel layout plots under: %s", layout_dir)
    common = {int(el) for el in common_electrodes}

    pitch_um = None
    electrode_size_x_um = None
    electrode_size_y_um = None
    electrodes_per_well = None
    physical_w_um = None
    physical_h_um = None
    if isinstance(hdmea_geometry, dict):
        try:
            if hdmea_geometry.get("pitch_um") is not None:
                pitch_um = float(hdmea_geometry.get("pitch_um"))
        except Exception:
            pitch_um = None
        try:
            if hdmea_geometry.get("physical_width_mm") is not None:
                physical_w_um = float(hdmea_geometry.get("physical_width_mm")) * 1000.0
        except Exception:
            physical_w_um = None
        try:
            if hdmea_geometry.get("physical_height_mm") is not None:
                physical_h_um = float(hdmea_geometry.get("physical_height_mm")) * 1000.0
        except Exception:
            physical_h_um = None

        try:
            if hdmea_geometry.get("active_width_mm") is not None:
                physical_w_um = float(hdmea_geometry.get("active_width_mm")) * 1000.0
        except Exception:
            pass
        try:
            if hdmea_geometry.get("active_height_mm") is not None:
                physical_h_um = float(hdmea_geometry.get("active_height_mm")) * 1000.0
        except Exception:
            pass

        try:
            if hdmea_geometry.get("electrode_size_x_um") is not None:
                electrode_size_x_um = float(hdmea_geometry.get("electrode_size_x_um"))
        except Exception:
            electrode_size_x_um = None
        try:
            if hdmea_geometry.get("electrode_size_y_um") is not None:
                electrode_size_y_um = float(hdmea_geometry.get("electrode_size_y_um"))
        except Exception:
            electrode_size_y_um = None
        try:
            if hdmea_geometry.get("electrodes_per_well") is not None:
                electrodes_per_well = int(hdmea_geometry.get("electrodes_per_well"))
            elif hdmea_geometry.get("electrode_count") is not None:
                electrodes_per_well = int(hdmea_geometry.get("electrode_count"))
        except Exception:
            electrodes_per_well = None

    def _square_marker_area_from_geometry(x_span_um: float, y_span_um: float, fig_w_in: float, fig_h_in: float) -> tuple[float, float]:
        ax_w_in = max(1.0, float(fig_w_in) * 0.82)
        ax_h_in = max(1.0, float(fig_h_in) * 0.80)
        x_scale_pt_per_um = (ax_w_in * 72.0) / max(float(x_span_um), 1.0)
        y_scale_pt_per_um = (ax_h_in * 72.0) / max(float(y_span_um), 1.0)

        if electrode_size_x_um is not None and electrode_size_y_um is not None and electrode_size_x_um > 0 and electrode_size_y_um > 0:
            electrode_side_um = (float(electrode_size_x_um) * float(electrode_size_y_um)) ** 0.5
        elif electrodes_per_well is not None and electrodes_per_well > 0:
            electrode_area_um2 = (float(x_span_um) * float(y_span_um)) / float(electrodes_per_well)
            electrode_side_um = max(1.0, electrode_area_um2 ** 0.5)
        elif pitch_um is not None and pitch_um > 0:
            electrode_side_um = float(pitch_um)
        else:
            electrode_side_um = 12.0

        side_pt = max(0.8, min(8.0, float(electrode_side_um) * min(x_scale_pt_per_um, y_scale_pt_per_um)))
        return float(side_pt * side_pt), float(side_pt)

    for rec_name in rec_names:
        try:
            rec = se.read_maxwell(file_path=str(h5_path), stream_id=stream_id, rec_name=rec_name)
            cv = rec.get_property("contact_vector")
            electrodes = np.asarray(cv["electrode"], dtype=int)
            x, y = _extract_xy_from_contact_vector(cv)
            x = np.asarray(x, dtype=float)
            y = np.asarray(y, dtype=float)
            if electrodes.size == 0:
                continue
            is_common = np.asarray([int(el) in common for el in electrodes], dtype=bool)

            # Keep native contact_vector coordinates unchanged for plotting.
            x_plot = x
            y_plot = y
            x0 = float(np.nanmin(x_plot)) if x_plot.size > 0 else 0.0
            y0 = float(np.nanmin(y_plot)) if y_plot.size > 0 else 0.0

            x_span_data = float(np.nanmax(x_plot) - np.nanmin(x_plot)) if x_plot.size > 0 else 1.0
            y_span_data = float(np.nanmax(y_plot) - np.nanmin(y_plot)) if y_plot.size > 0 else 1.0
            plot_w_um = float(physical_w_um) if physical_w_um is not None else max(x_span_data, 1.0)
            plot_h_um = float(physical_h_um) if physical_h_um is not None else max(y_span_data, 1.0)
            aspect = max(0.25, min(6.0, plot_w_um / max(plot_h_um, 1.0)))
            fig_h = 6.2
            fig_w = max(6.0, min(14.0, fig_h * aspect))
            marker_segment, marker_side_pt = _square_marker_area_from_geometry(
                x_span_um=plot_w_um,
                y_span_um=plot_h_um,
                fig_w_in=fig_w,
                fig_h_in=fig_h,
            )
            marker_common = max(1.0, marker_segment * 1.15)

            fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=150)
            ax.scatter(x_plot[~is_common], y_plot[~is_common], s=marker_segment, marker="s", c="#b9b9b9", alpha=0.55, label="segment-only")
            ax.scatter(x_plot[is_common], y_plot[is_common], s=marker_common, marker="s", c="#d62728", alpha=0.9, label="common")
            ax.set_title(f"{stream_id} {rec_name} channel layout")
            ax.set_xlabel("x (um)")
            ax.set_ylabel("y (um)")
            ax.legend(loc="best", frameon=False)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlim(0.0, float(plot_w_um))
            ax.set_ylim(0.0, float(plot_h_um))
            ax.grid(alpha=0.15)
            fig.tight_layout()

            if isinstance(hdmea_geometry, dict):
                log.debug(
                    "HDMEA layout style rec=%s fig=(%.2f,%.2f) axis_origin_um=(0.0,0.0) axis_span_um=(%.1f,%.1f) data_min_um=(%.1f,%.1f) aspect=%.3f marker_side_pt=%.2f pitch_um=%s active_mm=(%s,%s)",
                    rec_name,
                    float(fig_w),
                    float(fig_h),
                    float(plot_w_um),
                    float(plot_h_um),
                    float(x0),
                    float(y0),
                    float(aspect),
                    float(marker_side_pt),
                    str(hdmea_geometry.get("pitch_um")),
                    str(hdmea_geometry.get("active_width_mm")),
                    str(hdmea_geometry.get("active_height_mm")),
                )

            out_path = layout_dir / f"channel_layout_{_safe_artifact_token(stream_id)}_{_safe_artifact_token(rec_name)}.png"
            fig.savefig(out_path)
            plt.close(fig)
            log.debug("Wrote channel layout plot: %s", out_path)
        except Exception as e:
            log.warning("Failed channel layout plot for rec=%s: %s", rec_name, e)


def _save_channel_layout_overlap_heatmap(
    *,
    h5_path: Path,
    stream_id: str,
    rec_names: list[str],
    hdmea_geometry: Optional[dict[str, object]],
    out_path: Path,
    logger: logging.Logger,
) -> None:
    log = logger.getChild("_save_channel_layout_overlap_heatmap")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import BoundaryNorm, ListedColormap
        import numpy as np
        import spikeinterface.extractors as se
    except Exception as e:
        log.warning("Skipping channel overlap heatmap: plotting dependencies unavailable (%s)", e)
        return

    if not rec_names:
        log.warning("Skipping channel overlap heatmap: no segment names available")
        return

    _ensure_maxwell_hdf5_plugin_path()

    active_w_um = None
    active_h_um = None
    pitch_um = None
    electrode_size_x_um = None
    electrode_size_y_um = None
    electrodes_per_well = None
    if isinstance(hdmea_geometry, dict):
        try:
            if hdmea_geometry.get("active_width_mm") is not None:
                active_w_um = float(hdmea_geometry.get("active_width_mm")) * 1000.0
        except Exception:
            active_w_um = None
        try:
            if hdmea_geometry.get("active_height_mm") is not None:
                active_h_um = float(hdmea_geometry.get("active_height_mm")) * 1000.0
        except Exception:
            active_h_um = None
        try:
            if hdmea_geometry.get("pitch_um") is not None:
                pitch_um = float(hdmea_geometry.get("pitch_um"))
        except Exception:
            pitch_um = None
        try:
            if hdmea_geometry.get("electrode_size_x_um") is not None:
                electrode_size_x_um = float(hdmea_geometry.get("electrode_size_x_um"))
        except Exception:
            electrode_size_x_um = None
        try:
            if hdmea_geometry.get("electrode_size_y_um") is not None:
                electrode_size_y_um = float(hdmea_geometry.get("electrode_size_y_um"))
        except Exception:
            electrode_size_y_um = None
        try:
            if hdmea_geometry.get("electrodes_per_well") is not None:
                electrodes_per_well = int(hdmea_geometry.get("electrodes_per_well"))
            elif hdmea_geometry.get("electrode_count") is not None:
                electrodes_per_well = int(hdmea_geometry.get("electrode_count"))
        except Exception:
            electrodes_per_well = None

    participation: dict[int, int] = {}
    xy_by_electrode: dict[int, tuple[float, float]] = {}

    for rec_name in rec_names:
        rec = se.read_maxwell(file_path=str(h5_path), stream_id=stream_id, rec_name=str(rec_name))
        cv = rec.get_property("contact_vector")
        electrodes = np.asarray(cv["electrode"], dtype=int)
        x = np.asarray(cv["x"], dtype=float)
        y = np.asarray(cv["y"], dtype=float)
        for el, xx, yy in zip(electrodes.tolist(), x.tolist(), y.tolist(), strict=False):
            key = int(el)
            participation[key] = int(participation.get(key, 0)) + 1
            if key not in xy_by_electrode:
                xy_by_electrode[key] = (float(xx), float(yy))

    if not participation:
        log.warning("Skipping channel overlap heatmap: no channel participation data computed")
        return

    items = sorted(
        ((el, xy_by_electrode[el], int(cnt)) for el, cnt in participation.items() if el in xy_by_electrode),
        key=lambda it: (float(it[1][1]), float(it[1][0]), int(it[0])),
    )
    x = np.asarray([float(it[1][0]) for it in items], dtype=float)
    y = np.asarray([float(it[1][1]) for it in items], dtype=float)
    counts = np.asarray([int(it[2]) for it in items], dtype=float)
    n_segments = max(1, int(len(rec_names)))

    x_span_data = float(np.nanmax(x) - np.nanmin(x)) if x.size > 0 else 1.0
    y_span_data = float(np.nanmax(y) - np.nanmin(y)) if y.size > 0 else 1.0
    plot_w_um = float(active_w_um) if active_w_um is not None else max(x_span_data, 1.0)
    plot_h_um = float(active_h_um) if active_h_um is not None else max(y_span_data, 1.0)
    aspect = max(0.25, min(6.0, float(plot_w_um) / max(float(plot_h_um), 1.0)))
    fig_h = 6.2
    fig_w = max(6.0, min(14.0, fig_h * aspect))

    ax_w_in = max(1.0, float(fig_w) * 0.82)
    ax_h_in = max(1.0, float(fig_h) * 0.80)
    x_scale_pt_per_um = (ax_w_in * 72.0) / max(float(plot_w_um), 1.0)
    y_scale_pt_per_um = (ax_h_in * 72.0) / max(float(plot_h_um), 1.0)
    if electrode_size_x_um is not None and electrode_size_y_um is not None and electrode_size_x_um > 0 and electrode_size_y_um > 0:
        electrode_side_um = (float(electrode_size_x_um) * float(electrode_size_y_um)) ** 0.5
    elif electrodes_per_well is not None and electrodes_per_well > 0:
        electrode_area_um2 = (float(plot_w_um) * float(plot_h_um)) / float(electrodes_per_well)
        electrode_side_um = max(1.0, electrode_area_um2 ** 0.5)
    elif pitch_um is not None and pitch_um > 0:
        electrode_side_um = float(pitch_um)
    else:
        electrode_side_um = 12.0
    marker_side_pt = max(0.8, min(8.0, float(electrode_side_um) * min(x_scale_pt_per_um, y_scale_pt_per_um)))
    marker_area = float(marker_side_pt * marker_side_pt)

    if int(n_segments) <= 1:
        color_steps = ["#b9b9b9"]
    else:
        # Keep step=1 truly gray, then use a high-contrast red ramp for steps >=2.
        reds = plt.get_cmap("Reds")
        red_steps = [reds(float(v)) for v in np.linspace(0.55, 0.98, int(n_segments - 1))]
        color_steps = ["#b9b9b9", *red_steps]
    cmap = ListedColormap(color_steps, name="seg_overlap_steps")
    boundaries = np.arange(0.5, float(n_segments) + 1.5, 1.0)
    norm = BoundaryNorm(boundaries=boundaries, ncolors=cmap.N, clip=True)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=150)
    sc = ax.scatter(x, y, c=counts, cmap=cmap, norm=norm, marker="s", s=marker_area, alpha=0.95)
    ax.set_title(f"{stream_id} channel overlap across segments")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(0.0, float(plot_w_um))
    ax.set_ylim(0.0, float(plot_h_um))
    ax.grid(alpha=0.15)

    cbar = fig.colorbar(sc, ax=ax, boundaries=boundaries, spacing="proportional")
    cbar.set_label("segment participation count")
    ticks = list(range(1, int(max(1, n_segments)) + 1))
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([str(t) for t in ticks])

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
    log.debug(
        "Overlap heatmap marker sizing: marker_side_pt=%.3f marker_area=%.3f electrode_side_um=%.3f",
        float(marker_side_pt),
        float(marker_area),
        float(electrode_side_um),
    )
    log.info("Wrote channel overlap heatmap: %s", out_path)


def _save_single_trace_panel(*, traces, channel_ids, title: str, out_path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    traces = np.asarray(traces, dtype=float)
    if traces.ndim != 2 or traces.shape[0] <= 0 or traces.shape[1] <= 0:
        return

    n_samples = int(traces.shape[0])
    t = np.arange(n_samples, dtype=float)

    centered = traces - np.nanmedian(traces, axis=0, keepdims=True)
    per_ch_scale = np.nanpercentile(np.abs(centered), 95, axis=0)
    scale = float(np.nanmedian(per_ch_scale[per_ch_scale > 0])) if np.any(per_ch_scale > 0) else 1.0
    scale = max(scale, 1e-9)
    norm = centered / scale

    n_ch = int(norm.shape[1])
    offsets = np.arange(n_ch, dtype=float) * 4.0

    fig_h = min(16.0, max(6.0, 0.28 * n_ch))
    fig, ax = plt.subplots(figsize=(14, fig_h), dpi=150)
    for idx in range(n_ch):
        ax.plot(t, norm[:, idx] + offsets[idx], lw=0.65, color="#1f77b4", alpha=0.92)

    tick_stride = max(1, n_ch // 20)
    tick_idx = np.arange(0, n_ch, tick_stride, dtype=int)
    ax.set_yticks(offsets[tick_idx])
    ax.set_yticklabels([str(channel_ids[i]) for i in tick_idx], fontsize=7)
    ax.set_xlabel("samples")
    ax.set_ylabel("channel id")
    ax.set_title(title)
    ax.grid(alpha=0.2, axis="x")
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def _save_trace_panel_with_axis(
    *,
    x,
    traces,
    channel_ids,
    title: str,
    x_label: str,
    out_path: Path,
    vlines: Optional[list[float]] = None,
    break_on_gaps: bool = False,
    gap_threshold: float = 1.0,
    logger: Optional[logging.Logger] = None,
) -> None:
    log = logger.getChild("_save_trace_panel_with_axis") if logger is not None else None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    x = np.asarray(x, dtype=float)
    traces = np.asarray(traces, dtype=float)
    if traces.ndim != 2 or traces.shape[0] <= 0 or traces.shape[1] <= 0:
        return
    if x.ndim != 1 or int(x.size) != int(traces.shape[0]):
        return

    centered = traces - np.nanmedian(traces, axis=0, keepdims=True)
    per_ch_scale = np.nanpercentile(np.abs(centered), 95, axis=0)
    scale = float(np.nanmedian(per_ch_scale[per_ch_scale > 0])) if np.any(per_ch_scale > 0) else 1.0
    scale = max(scale, 1e-9)
    norm = centered / scale

    n_ch = int(norm.shape[1])
    offsets = np.arange(n_ch, dtype=float) * 4.0

    if bool(break_on_gaps) and int(x.size) > 1:
        gap_idx = np.flatnonzero(np.diff(x) > float(gap_threshold))
        if log is not None:
            log.debug(
                "Trace gap detection: threshold=%.3f n_gaps=%d n_points=%d",
                float(gap_threshold),
                int(gap_idx.size),
                int(x.size),
            )
        if int(gap_idx.size) > 0:
            new_len = int(x.size) + int(gap_idx.size)
            x_plot = np.full(new_len, np.nan, dtype=float)
            norm_plot = np.full((new_len, n_ch), np.nan, dtype=float)
            src_i = 0
            dst_i = 0
            gap_set = {int(i) for i in gap_idx.tolist()}
            while src_i < int(x.size):
                x_plot[dst_i] = float(x[src_i])
                norm_plot[dst_i, :] = norm[src_i, :]
                if src_i in gap_set:
                    dst_i += 1
                src_i += 1
                dst_i += 1
            if log is not None:
                log.debug("Inserted NaN plot breaks across %d gaps (new_points=%d)", int(gap_idx.size), int(new_len))
        else:
            x_plot = x
            norm_plot = norm
    else:
        x_plot = x
        norm_plot = norm

    fig_h = min(16.0, max(6.0, 0.28 * n_ch))
    fig, ax = plt.subplots(figsize=(14, fig_h), dpi=150)
    for idx in range(n_ch):
        ax.plot(x_plot, norm_plot[:, idx] + offsets[idx], lw=0.65, color="#1f77b4", alpha=0.92)

    if vlines:
        for xline in vlines:
            ax.axvline(float(xline), color="red", linestyle=":", linewidth=1.2, alpha=0.9)

    tick_stride = max(1, n_ch // 20)
    tick_idx = np.arange(0, n_ch, tick_stride, dtype=int)
    ax.set_yticks(offsets[tick_idx])
    ax.set_yticklabels([str(channel_ids[i]) for i in tick_idx], fontsize=7)
    ax.set_xlabel(x_label)
    ax.set_ylabel("channel id")
    ax.set_title(title)
    ax.grid(alpha=0.2, axis="x")
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
    if log is not None:
        log.info("Wrote trace plot: %s", out_path)


def _save_trace_diagnostics(
    *,
    h5_path: Path,
    stream_id: str,
    rec_names: list[str],
    rec_list: list[object],
    multirecording,
    out_dir: Path,
    segment_traces_dir: Optional[Path] = None,
    concat_cluster_reps_path: Optional[Path] = None,
    segment_filename_prefix: str = "segment_trace",
    concat_title_label: str = "concatenated representative traces",
    n_trace_channels: int = 24,
    logger: logging.Logger,
) -> None:
    log = logger.getChild("_save_trace_diagnostics")
    try:
        import numpy as np
    except Exception as e:
        log.warning("Skipping trace diagnostics plots: plotting dependencies unavailable (%s)", e)
        return

    safe_stream = _safe_artifact_token(stream_id)
    segment_dir = Path(segment_traces_dir) if segment_traces_dir is not None else (Path(out_dir) / "segment_traces")
    if not rec_list:
        return

    selected_channels = _choose_plot_channels(
        channel_ids=rec_list[0].get_channel_ids(),
        n_trace_channels=int(n_trace_channels),
        rng_seed=0,
    )
    log.info(
        "Trace diagnostics settings: n_trace_channels=%d selected=%d",
        int(n_trace_channels),
        int(selected_channels.size),
    )
    log.debug("Selected trace channels: %s", [str(ch) for ch in selected_channels.tolist()])
    if selected_channels.size == 0:
        log.warning("Skipping trace diagnostics plots: no channels available")
        return

    concat_x_segments_abs_sec: list[np.ndarray] = []
    concat_traces_segments: list[np.ndarray] = []
    concat_segment_starts_abs_sec: list[float] = []
    segment_timing_summary: list[tuple[str, float, float, int]] = []

    start_ms_values: list[int] = []
    for rec_name in rec_names:
        try:
            _info = _read_well_rec_frame_nos_and_trigger_settings(
                h5_path=h5_path,
                stream_id=stream_id,
                rec_name=rec_name,
            )
            start_ms_values.append(int(_info.get("start_ms", 0)))
        except Exception:
            continue
    t0_ms = int(min(start_ms_values)) if start_ms_values else 0

    for rec_name, rec in zip(rec_names, rec_list, strict=False):
        try:
            fs_extractor_hz = float(rec.get_sampling_frequency())
            frame_info = _read_well_rec_frame_nos_and_trigger_settings(
                h5_path=h5_path,
                stream_id=stream_id,
                rec_name=rec_name,
            )
            start_ms = int(frame_info.get("start_ms", 0))
            stop_ms = int(frame_info.get("stop_ms", start_ms))
            fs_h5_hz = float(frame_info.get("sampling_hz", fs_extractor_hz))
            fs_hz = float(fs_h5_hz if fs_h5_hz > 0 else fs_extractor_hz)
            frame_nos = np.asarray(frame_info.get("frame_nos", []), dtype=float)
            log.debug("Segment %s frame_nos loaded: n=%d", rec_name, int(frame_nos.size))
            if abs(float(fs_extractor_hz) - float(fs_hz)) > 1e-6:
                log.warning(
                    "Sampling-rate mismatch for %s: h5=%.6fHz extractor=%.6fHz (using h5)",
                    rec_name,
                    float(fs_hz),
                    float(fs_extractor_hz),
                )
            else:
                log.debug("Segment %s sampling_hz from h5 metadata: %.6f", rec_name, float(fs_hz))

            traces, channel_ids = _extract_selected_traces(
                recording=rec,
                selected_channel_ids=selected_channels,
                max_samples=None,
            )
            if traces.size == 0:
                continue
            n = min(int(traces.shape[0]), int(frame_nos.size))
            if n <= 0:
                continue
            traces = traces[:n, :]
            seg_anchor_sec = (float(start_ms) - float(t0_ms)) / 1000.0
            frame_rel_sec = (frame_nos[:n] - float(frame_nos[0])) / fs_hz
            x_segment_seconds_abs = seg_anchor_sec + frame_rel_sec

            segment_timing_summary.append((str(rec_name), float(start_ms), float(stop_ms), int(n)))
            seg_out = segment_dir / (
                f"{_safe_artifact_token(segment_filename_prefix)}_{safe_stream}_{_safe_artifact_token(rec_name)}.png"
            )
            _save_trace_panel_with_axis(
                x=x_segment_seconds_abs,
                traces=traces,
                channel_ids=channel_ids,
                title=(
                    f"{stream_id} {rec_name} {_safe_artifact_token(segment_filename_prefix)} sample "
                    f"(n_ch={int(channel_ids.size)}, frame gaps visible)"
                ),
                x_label="time (s)",
                out_path=seg_out,
                break_on_gaps=True,
                gap_threshold=(1.5 / fs_hz),
                logger=log,
            )

            concat_segment_starts_abs_sec.append(float(seg_anchor_sec))
            log.debug("Concat segment start marker (abs sec): rec=%s x=%.3f", rec_name, float(seg_anchor_sec))
            concat_x_segments_abs_sec.append(x_segment_seconds_abs)
            concat_traces_segments.append(traces)
        except Exception as e:
            log.warning("Failed segment trace plot for rec=%s: %s", rec_name, e)

    if concat_x_segments_abs_sec and concat_traces_segments:
        try:
            x_concat_all_seconds_abs = np.concatenate(concat_x_segments_abs_sec)
            traces_concat_all = np.concatenate(concat_traces_segments, axis=0)
            x_concat_all_minutes = x_concat_all_seconds_abs / 60.0
            fs_first_hz = 1.0
            try:
                first_info = _read_well_rec_frame_nos_and_trigger_settings(
                    h5_path=h5_path,
                    stream_id=stream_id,
                    rec_name=str(rec_names[0]),
                )
                fs_first_hz = max(1.0, float(first_info.get("sampling_hz", 1.0)))
            except Exception:
                fs_first_hz = max(1.0, float(rec_list[0].get_sampling_frequency()))
            concat_out = (
                Path(concat_cluster_reps_path)
                if concat_cluster_reps_path is not None
                else (Path(out_dir) / f"concat_cluster_reps_{safe_stream}.png")
            )
            _save_trace_panel_with_axis(
                x=x_concat_all_minutes,
                traces=traces_concat_all,
                channel_ids=selected_channels,
                title=(
                    f"{stream_id} {concat_title_label} "
                    f"(n_ch={int(selected_channels.size)}, segment starts marked)"
                ),
                x_label="time (min)",
                out_path=concat_out,
                vlines=[float(v) / 60.0 for v in concat_segment_starts_abs_sec],
                break_on_gaps=True,
                gap_threshold=(1.5 / (60.0 * max(1.0, float(fs_first_hz)))),
                logger=log,
            )
            log.debug("Concat segment starts (abs sec): %s", [float(v) for v in concat_segment_starts_abs_sec])

            if segment_timing_summary:
                earliest_start_s = (min(float(s) for _, s, _, _ in segment_timing_summary) - float(t0_ms)) / 1000.0
                latest_stop_s = (max(float(e) for _, _, e, _ in segment_timing_summary) - float(t0_ms)) / 1000.0
                log.info(
                    "Timing summary: segments=%d elapsed_sec=%.3f start_sec=%.3f stop_sec=%.3f",
                    int(len(segment_timing_summary)),
                    float(max(0.0, latest_stop_s - earliest_start_s)),
                    float(earliest_start_s),
                    float(latest_stop_s),
                )
        except Exception as e:
            log.warning("Failed concatenated representative trace plot: %s", e)


def find_common_electrodes_from_segments(*, h5_path: Path, stream_id: str) -> tuple[list[str], list[int]]:
    try:
        import spikeinterface.extractors as se
    except Exception as e:  # pragma: no cover
        raise RuntimeError("raw preprocessing requires `h5py` and `spikeinterface` installed") from e

    _ensure_maxwell_hdf5_plugin_path()

    rec_names = list_stream_recording_names(h5_path=h5_path, stream_id=stream_id)

    common: Optional[set[int]] = None
    for rec_name in rec_names:
        rec = se.read_maxwell(file_path=str(h5_path), stream_id=stream_id, rec_name=rec_name)
        electrodes = rec.get_property("contact_vector")["electrode"]
        electrode_set = {int(x) for x in electrodes}
        common = electrode_set if common is None else (common & electrode_set)

    return rec_names, sorted(common or set())


def _get_runtime_logger(*, logger: Optional[logging.Logger] = None) -> logging.Logger:
    runtime_logger = logging.getLogger(__name__)
    if logger is not None:
        runtime_logger.setLevel(logging.DEBUG)
        runtime_logger.propagate = False
        for handler in logger.handlers:
            if handler not in runtime_logger.handlers:
                runtime_logger.addHandler(handler)
    return runtime_logger


def _process_rec_segment_for_concatenation(
    *,
    h5_path: Path,
    stream_id: str,
    rec_name: str,
    common_el: list[int],
    center_chunk_size: int,
    expected_xy_by_electrode: Optional[dict[int, tuple[float, float]]] = None,
    expected_xy_atol: float = 0.0,
):
    try:
        import numpy as np
        import spikeinterface.extractors as se
        import spikeinterface.full as si
    except Exception as e:  # pragma: no cover
        raise RuntimeError("raw preprocessing requires `numpy` and `spikeinterface` installed") from e

    _ensure_maxwell_hdf5_plugin_path()

    rec = se.read_maxwell(file_path=str(h5_path), stream_id=stream_id, rec_name=rec_name)
    fs = float(rec.get_sampling_frequency())
    n_samples = int(rec.get_num_samples())

    chunk = min(center_chunk_size, n_samples) - 100
    chunk = max(chunk, 100)
    rec_centered = si.center(rec, chunk_size=chunk)

    rec_el = np.asarray(rec.get_property("contact_vector")["electrode"], dtype=int)
    if int(np.unique(rec_el).size) != int(rec_el.size):
        raise RuntimeError(f"Duplicate electrode ids found in contact_vector for segment {rec_name}")

    el_to_idx = {int(el): int(i) for i, el in enumerate(rec_el)}
    chan_idx = [el_to_idx[int(el)] for el in common_el]
    sel_channels = np.asarray(rec.get_channel_ids(), dtype=object)[chan_idx]
    processed = rec_centered.select_channels(list(sel_channels))

    expected_el = np.asarray(common_el, dtype=int)
    processed_el = np.asarray(processed.get_property("contact_vector")["electrode"], dtype=int)
    if processed_el.shape != expected_el.shape or not np.array_equal(processed_el, expected_el):
        raise RuntimeError(f"Selected electrodes mismatch for segment {rec_name}")

    if expected_xy_by_electrode is not None:
        cv = processed.get_property("contact_vector")
        x, y = _extract_xy_from_contact_vector(cv)
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        expected_x = np.asarray([expected_xy_by_electrode[int(el)][0] for el in expected_el], dtype=float)
        expected_y = np.asarray([expected_xy_by_electrode[int(el)][1] for el in expected_el], dtype=float)
        if not (
            np.allclose(x, expected_x, atol=float(expected_xy_atol))
            and np.allclose(y, expected_y, atol=float(expected_xy_atol))
        ):
            raise RuntimeError(f"Electrode x/y locations differ from reference for segment {rec_name}")

    processed = processed.rename_channels([int(el) for el in expected_el])
    return processed, {"rec_name": rec_name, "fs": fs, "n_samples": n_samples}


def _read_well_rec_frame_nos_and_trigger_settings(*, h5_path: Path, stream_id: str, rec_name: str):
    return read_well_rec_frame_nos_and_trigger_settings(
        h5_path=h5_path,
        stream_id=stream_id,
        rec_name=rec_name,
    )


def _ensure_maxwell_hdf5_plugin_path(*, prefix: str = "[mea_analysis]") -> None:
    ensure_maxwell_hdf5_plugin_path(prefix=prefix)


def _apply_standard_preprocessing(*, recording, logger: logging.Logger):
    log = logger.getChild("_apply_standard_preprocessing")
    try:
        import spikeinterface.preprocessing as spre
    except Exception as e:  # pragma: no cover
        raise RuntimeError("preprocessing requires `spikeinterface.preprocessing` installed") from e

    rec = recording
    try:
        dtype_str = str(rec.get_dtype())
    except Exception:
        dtype_str = ""
    if dtype_str.startswith("uint"):
        rec = spre.unsigned_to_signed(rec)

    rec = spre.highpass_filter(rec, freq_min=300.0)
    try:
        rec = spre.common_reference(rec, reference="local", operator="median", local_radius=(250, 250))
    except Exception as e:
        log.warning("Local common_reference failed; falling back to global median reference (%s)", e)
        rec = spre.common_reference(rec, reference="global", operator="median")

    try:
        rec.annotate(is_filtered=True)
    except Exception:
        pass

    try:
        dtype_after = str(rec.get_dtype())
    except Exception:
        dtype_after = ""
    if dtype_after != "float32":
        rec = spre.astype(rec, "float32")
    return rec


def build_concatenated_recording(
    *,
    h5_path: Path,
    stream_id: str,
    n_jobs: int = 8,
    center_chunk_size: int = 10_000,
    temporal_resample_factor: Optional[int] = None,
    temporal_resample_rate_hz: Optional[int] = None,
    temporal_resample_margin_ms: float = 100.0,
    temporal_resample_dtype: Optional[str] = None,
    plot_output_dir: Optional[Path] = None,
    channel_layouts_output_dir: Optional[Path] = None,
    channel_layout_heatmap_output_path: Optional[Path] = None,
    channel_layout_heatmap_enabled: bool = False,
    plot_channel_layouts: bool = True,
    segment_traces_output_dir: Optional[Path] = None,
    preprocessed_segment_traces_output_dir: Optional[Path] = None,
    concat_cluster_reps_output_path: Optional[Path] = None,
    preprocessed_concat_cluster_reps_output_path: Optional[Path] = None,
    hdmea_geometry: Optional[dict[str, object]] = None,
    n_trace_channels: int = 24,
    plot_centered_traces: bool = True,
    plot_preprocessed_traces: bool = True,
    epoch_markers_output_dir: Optional[Path] = None,
    return_artifacts: bool = False,
    logger: Optional[logging.Logger] = None,
) -> tuple[object, list[int]] | tuple[object, list[int], dict[str, object]]:
    """Minimal stg1 runtime: segment load/center/common-electrode concat + epoch markers."""

    try:
        import numpy as np
        import spikeinterface.extractors as se
        import spikeinterface.full as si
        import spikeinterface.preprocessing as spre
    except Exception as e:  # pragma: no cover
        raise RuntimeError(
            "raw preprocessing requires `numpy`, `h5py`, and `spikeinterface` installed"
        ) from e

    log = _get_runtime_logger(logger=logger)

    _ensure_maxwell_hdf5_plugin_path()
    log.debug(
        "Starting build_concatenated_recording: h5_path=%s stream_id=%s n_jobs=%s",
        h5_path,
        stream_id,
        int(n_jobs),
    )

    rec_names = list_stream_recording_names(h5_path=h5_path, stream_id=stream_id)
    if not rec_names:
        raise RuntimeError(f"No recording segments found under /wells/{stream_id} in {h5_path}")
    log.info("Discovered %d Maxwell recording segment(s)", int(len(rec_names)))

    is_single_segment = len(rec_names) == 1
    if is_single_segment:
        rec = se.read_maxwell(file_path=str(h5_path), stream_id=stream_id, rec_name=rec_names[0])
        chunk = max(min(center_chunk_size, int(rec.get_num_samples())) - 100, 100)
        rec_centered = si.center(rec, chunk_size=chunk)
        rec_el = np.asarray(rec.get_property("contact_vector")["electrode"], dtype=int)
        common_el = [int(e) for e in rec_el.tolist()]
        rec_list = [rec_centered.rename_channels(common_el)]
    else:
        rec_names, common_el = find_common_electrodes_from_segments(h5_path=h5_path, stream_id=stream_id)
        rec0 = se.read_maxwell(file_path=str(h5_path), stream_id=stream_id, rec_name=rec_names[0])
        cv0 = rec0.get_property("contact_vector")
        el0 = np.asarray(cv0["electrode"], dtype=int)
        x0, y0 = _extract_xy_from_contact_vector(cv0)
        expected_xy_by_electrode = {
            int(e): (float(x), float(y))
            for e, x, y in zip(el0, np.asarray(x0, dtype=float), np.asarray(y0, dtype=float), strict=False)
        }

        from concurrent.futures import ThreadPoolExecutor

        max_workers = min(len(rec_names), max(1, int(n_jobs)))
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            results = list(
                ex.map(
                    lambda rn: _process_rec_segment_for_concatenation(
                        h5_path=h5_path,
                        stream_id=stream_id,
                        rec_name=rn,
                        common_el=common_el,
                        center_chunk_size=center_chunk_size,
                        expected_xy_by_electrode=expected_xy_by_electrode,
                        expected_xy_atol=0.0,
                    ),
                    rec_names,
                )
            )

        rec_list = [r for r, _ in results]

    centered_rec_list = rec_list
    centered_concat = centered_rec_list[0] if is_single_segment else si.concatenate_recordings(centered_rec_list)
    preprocessed_rec_list = [_apply_standard_preprocessing(recording=rec, logger=log) for rec in centered_rec_list]
    multirecording = preprocessed_rec_list[0] if is_single_segment else si.concatenate_recordings(preprocessed_rec_list)

    seg_lengths = [int(r.get_num_samples()) for r in preprocessed_rec_list]
    seg_offsets: list[int] = []
    acc = 0
    for n_frames in seg_lengths:
        seg_offsets.append(acc)
        acc += int(n_frames)

    concat_epochs: list[dict] = []
    for i, (rn, n_frames, start) in enumerate(zip(rec_names, seg_lengths, seg_offsets, strict=False)):
        concat_epochs.append(
            {
                "segment_index": int(i),
                "rec_name": str(rn),
                "start_sample": int(start),
                "end_sample": int(start + int(n_frames)),
                "n_samples": int(n_frames),
            }
        )

    maxwell_epochs: list[dict] = []
    for seg_index, rn in enumerate(rec_names):
        info = _read_well_rec_frame_nos_and_trigger_settings(h5_path=h5_path, stream_id=stream_id, rec_name=rn)
        frame_nos = info["frame_nos"]
        diffs = np.diff(frame_nos)
        split_points = np.flatnonzero(diffs != 1) + 1
        run_starts = np.concatenate(([0], split_points))
        run_ends = np.concatenate((split_points, [frame_nos.size]))
        seg_offset = int(seg_offsets[seg_index]) if seg_index < len(seg_offsets) else 0
        for rs, re in zip(run_starts, run_ends, strict=False):
            rs_i = int(rs)
            re_i = int(re)
            if re_i <= rs_i:
                continue
            maxwell_epochs.append(
                {
                    "segment_index": int(seg_index),
                    "rec_name": str(rn),
                    "start_sample": int(seg_offset + rs_i),
                    "end_sample": int(seg_offset + re_i),
                    "segment_start_sample": int(rs_i),
                    "segment_end_sample": int(re_i),
                    "frame_no_start": int(frame_nos[rs_i]),
                    "frame_no_end": int(frame_nos[re_i - 1]),
                }
            )

    if temporal_resample_rate_hz is not None or temporal_resample_factor is not None:
        old_fs = float(multirecording.get_sampling_frequency())
        if temporal_resample_rate_hz is not None:
            new_fs_int = int(temporal_resample_rate_hz)
        else:
            factor = int(temporal_resample_factor or 0)
            if factor < 2:
                raise ValueError(f"temporal_resample_factor must be >=2, got {factor}")
            new_fs_int = int(round(old_fs * factor))

        ratio = float(new_fs_int) / float(old_fs)
        dtype = None if temporal_resample_dtype is None else temporal_resample_dtype
        multirecording = spre.resample(
            multirecording,
            resample_rate=int(new_fs_int),
            margin_ms=float(temporal_resample_margin_ms),
            dtype=dtype,
            skip_checks=False,
        )

        def _scale_epoch(ep: dict) -> dict:
            out = dict(ep)
            for k in ("start_sample", "end_sample", "segment_start_sample", "segment_end_sample", "n_samples"):
                if k in out and out[k] is not None:
                    out[k] = int(round(float(out[k]) * ratio))
            return out

        concat_epochs = [_scale_epoch(ep) for ep in concat_epochs]
        maxwell_epochs = [_scale_epoch(ep) for ep in maxwell_epochs]

    if plot_output_dir is not None:
        out_dir = Path(plot_output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        log.info("Generating preprocess diagnostics plots under: %s", out_dir)
        if bool(plot_channel_layouts):
            _save_channel_layout_plots(
                h5_path=h5_path,
                stream_id=stream_id,
                rec_names=[str(rn) for rn in rec_names],
                common_electrodes=[int(el) for el in common_el],
                hdmea_geometry=hdmea_geometry,
                out_dir=out_dir,
                layout_dir=(Path(channel_layouts_output_dir) if channel_layouts_output_dir is not None else None),
                logger=log,
            )
        if bool(channel_layout_heatmap_enabled) and channel_layout_heatmap_output_path is not None:
            _save_channel_layout_overlap_heatmap(
                h5_path=h5_path,
                stream_id=stream_id,
                rec_names=[str(rn) for rn in rec_names],
                hdmea_geometry=hdmea_geometry,
                out_path=Path(channel_layout_heatmap_output_path),
                logger=log,
            )
        if bool(plot_centered_traces):
            _save_trace_diagnostics(
                h5_path=h5_path,
                stream_id=stream_id,
                rec_names=[str(rn) for rn in rec_names],
                rec_list=centered_rec_list,
                multirecording=centered_concat,
                out_dir=out_dir,
                segment_traces_dir=(Path(segment_traces_output_dir) if segment_traces_output_dir is not None else None),
                concat_cluster_reps_path=(
                    Path(concat_cluster_reps_output_path) if concat_cluster_reps_output_path is not None else None
                ),
                segment_filename_prefix="segment_trace",
                concat_title_label="concatenated representative traces (centered)",
                n_trace_channels=int(n_trace_channels),
                logger=log,
            )

        if bool(plot_preprocessed_traces):
            _save_trace_diagnostics(
                h5_path=h5_path,
                stream_id=stream_id,
                rec_names=[str(rn) for rn in rec_names],
                rec_list=preprocessed_rec_list,
                multirecording=multirecording,
                out_dir=out_dir,
                segment_traces_dir=(
                    Path(preprocessed_segment_traces_output_dir)
                    if preprocessed_segment_traces_output_dir is not None
                    else None
                ),
                concat_cluster_reps_path=(
                    Path(preprocessed_concat_cluster_reps_output_path)
                    if preprocessed_concat_cluster_reps_output_path is not None
                    else None
                ),
                segment_filename_prefix="preprocessed_segment_trace",
                concat_title_label="concatenated representative traces (preprocessed)",
                n_trace_channels=int(n_trace_channels),
                logger=log,
            )

    if epoch_markers_output_dir is not None:
        out_dir = Path(epoch_markers_output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        maxwell_path = out_dir / f"maxwell_contiguous_epochs_{stream_id}.json"
        concat_path = out_dir / f"concatenation_stitch_epochs_{stream_id}.json"

        with open(maxwell_path, "w", encoding="utf-8") as f:
            json.dump(list(maxwell_epochs), f, indent=2)
        log.info("Wrote maxwell epochs marker file: %s (n=%d)", maxwell_path, int(len(maxwell_epochs)))
        if not is_single_segment:
            with open(concat_path, "w", encoding="utf-8") as f:
                json.dump(list(concat_epochs), f, indent=2)
            log.info("Wrote concat epochs marker file: %s (n=%d)", concat_path, int(len(concat_epochs)))

    log.debug(
        "Completed build_concatenated_recording: segments=%d common_electrodes=%d",
        int(len(rec_names)),
        int(len(common_el)),
    )

    if bool(return_artifacts):
        return (
            multirecording,
            common_el,
            {
                "rec_names": [str(rn) for rn in rec_names],
                "centered_rec_list": centered_rec_list,
                "preprocessed_rec_list": preprocessed_rec_list,
                "centered_concat_recording": centered_concat,
                "preprocessed_concat_recording": multirecording,
            },
        )

    return multirecording, common_el


__all__ = [
    "RawPreprocessPlan",
    "discover_cfg_files",
    "parse_cfg_channel_locations",
    "build_preprocess_plan",
    "find_common_electrodes_from_segments",
    "_process_rec_segment_for_concatenation",
    "_read_well_rec_frame_nos_and_trigger_settings",
    "_ensure_maxwell_hdf5_plugin_path",
    "build_concatenated_recording",
]
