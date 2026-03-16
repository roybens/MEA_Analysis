from __future__ import annotations

import configparser
import json
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
    epoch_markers_output_dir: Optional[Path] = None,
) -> tuple[object, list[int]]:
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

    _ensure_maxwell_hdf5_plugin_path()

    rec_names = list_stream_recording_names(h5_path=h5_path, stream_id=stream_id)
    if not rec_names:
        raise RuntimeError(f"No recording segments found under /wells/{stream_id} in {h5_path}")

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

    multirecording = rec_list[0] if is_single_segment else si.concatenate_recordings(rec_list)

    seg_lengths = [int(r.get_num_samples()) for r in rec_list]
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

    if epoch_markers_output_dir is not None:
        out_dir = Path(epoch_markers_output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        maxwell_path = out_dir / f"maxwell_contiguous_epochs_{stream_id}.json"
        concat_path = out_dir / f"concatenation_stitch_epochs_{stream_id}.json"

        with open(maxwell_path, "w", encoding="utf-8") as f:
            json.dump(list(maxwell_epochs), f, indent=2)
        if not is_single_segment:
            with open(concat_path, "w", encoding="utf-8") as f:
                json.dump(list(concat_epochs), f, indent=2)

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
