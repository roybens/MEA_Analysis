from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class WaveformMultisegPlan:
    """Planning payload for optional multisegment waveform extraction flow."""

    enabled: bool
    segment_count: int
    mode: str


@dataclass(frozen=True)
class WaveformSortingSource:
    source_dir: Path
    source_kind: str


def _segment_count(recording: Any) -> int:
    try:
        return int(recording.get_num_segments())
    except Exception:
        return 1


def build_waveform_multiseg_plan(
    recording: Any,
    *,
    mode: str = "none",
) -> WaveformMultisegPlan:
    """Return an initial optional multiseg waveform extraction plan.

    This is a non-invasive scaffold to support phased migration of stage-3
    mechanics into stage-1 analyzer hooks.
    """

    segments = _segment_count(recording)
    token = str(mode or "none").strip().lower()
    enabled = token not in {"", "none", "off", "disabled"} and segments > 1
    return WaveformMultisegPlan(enabled=enabled, segment_count=segments, mode=token or "none")


def resolve_waveform_sorting_source_dir(
    *,
    output_dir: Path,
    merged_sorting_dir: str | Path | None = None,
    prefer_merged_sorting: bool = False,
) -> WaveformSortingSource:
    """Resolve canonical sorting source directory for waveform extraction.

    Supports both legacy stage-3 layout and stage-1-collapsed layout.
    """

    base = Path(output_dir)
    if merged_sorting_dir is not None:
        explicit = Path(merged_sorting_dir)
        if explicit.exists():
            return WaveformSortingSource(source_dir=explicit, source_kind="explicit_merged")

    candidates = [
        (base / "unitmatch_outputs" / "final_merged_sorting", "local_merged"),
        (base / "stg2_spikesorting_outputs" / "unitmatch_outputs" / "final_merged_sorting", "nested_merged"),
    ]
    if bool(prefer_merged_sorting):
        for path, kind in candidates:
            if path.exists():
                return WaveformSortingSource(source_dir=path, source_kind=kind)

    non_merged_candidates = [
        (base / "stg2_spikesorting_outputs" / "sorter_output", "nested_sorter_output"),
        (base / "sorter_output", "local_sorter_output"),
    ]
    for path, kind in non_merged_candidates:
        if path.exists():
            return WaveformSortingSource(source_dir=path, source_kind=kind)

    # Default path for clear error messaging in downstream loaders.
    return WaveformSortingSource(source_dir=base / "sorter_output", source_kind="default_missing")


def load_sorting_from_output_dir(*, sorter_output_dir: Path, sorter: str) -> Any:
    """Load a sorting extractor from a sorter output directory."""

    import spikeinterface.full as si  # type: ignore[import-not-found]

    if hasattr(si, "read_sorter_folder"):
        try:
            return si.read_sorter_folder(sorter_output_dir, sorter_name=sorter)
        except TypeError:
            try:
                return si.read_sorter_folder(sorter_output_dir, sorter)
            except Exception:
                pass
        except Exception:
            pass

    try:
        return si.load_extractor(sorter_output_dir)
    except Exception:
        pass

    raise RuntimeError(f"Could not load sorting from sorter_output_dir={sorter_output_dir} (sorter={sorter}).")


def load_preprocessed_recording_from_output_dir(
    *,
    output_dir: Path,
    preprocess_outputs_dirname: str = "stg1_preprocess_outputs",
) -> Any:
    """Load stage-1 preprocessed recording from a canonical output directory."""

    recording_dir = Path(output_dir) / str(preprocess_outputs_dirname) / "preprocessed_recording"
    if not recording_dir.exists():
        raise FileNotFoundError(f"preprocessed_recording not found: {recording_dir}")

    import spikeinterface.full as si  # type: ignore[import-not-found]

    try:
        return si.load(recording_dir)
    except Exception:
        return si.load_extractor(recording_dir)


def resolve_waveform_epoch_marker_paths(
    *,
    output_dir: Path,
    stream_id: str,
    preprocess_outputs_dirname: str = "stg1_preprocess_outputs",
) -> dict[str, Path]:
    """Resolve best-effort locations for waveform epoch marker JSON files."""

    base = Path(output_dir)
    stream = str(stream_id)

    local_maxwell = base / f"maxwell_contiguous_epochs_{stream}.json"
    local_concat = base / f"concatenation_stitch_epochs_{stream}.json"

    preprocess_dir = base / str(preprocess_outputs_dirname)
    nested_maxwell = preprocess_dir / f"maxwell_contiguous_epochs_{stream}.json"
    nested_concat = preprocess_dir / f"concatenation_stitch_epochs_{stream}.json"

    return {
        "maxwell": local_maxwell if local_maxwell.exists() else nested_maxwell,
        "concat": local_concat if local_concat.exists() else nested_concat,
    }


def harmonize_waveform_sorting(
    *,
    sorting: Any,
    recording: Any,
    debug_max_units: int | None = None,
    logger: Any | None = None,
) -> Any:
    """Apply best-effort sorting cleanup and optional debug unit limiting."""

    out = sorting

    try:
        import spikeinterface.full as si  # type: ignore[import-not-found]

        out = si.remove_excess_spikes(out, recording)
        out = out.remove_empty_units()
    except Exception:
        if logger is not None:
            try:
                logger.debug(
                    "Sorting cleanup (remove_excess_spikes/remove_empty_units) failed; continuing without cleanup.",
                    exc_info=True,
                )
            except Exception:
                pass

    try:
        if debug_max_units is not None:
            max_units = int(debug_max_units)
            if max_units > 0:
                try:
                    unit_ids = list(out.get_unit_ids())
                except Exception:
                    unit_ids = list(getattr(out, "unit_ids", []))

                if unit_ids:
                    try:
                        unit_ids_sorted = sorted(unit_ids)
                    except Exception:
                        unit_ids_sorted = sorted(unit_ids, key=lambda x: str(x))

                    keep = unit_ids_sorted[:max_units]
                    try:
                        out = out.select_units(unit_ids=keep)
                    except Exception:
                        try:
                            out = out.select_units(keep)
                        except Exception:
                            pass

                    if logger is not None:
                        try:
                            logger.warning(
                                "DEBUG: limiting waveforms to first %d units (of %d)",
                                int(len(keep)),
                                int(len(unit_ids_sorted)),
                            )
                        except Exception:
                            pass
    except Exception:
        pass

    return out


def epochs_to_intervals(epochs: list[dict]) -> list[tuple[int, int]]:
    intervals: list[tuple[int, int]] = []
    for e in epochs:
        try:
            start = int(e["start_sample"])
            end = int(e["end_sample"])
        except Exception:
            continue
        if end > start:
            intervals.append((start, end))
    intervals.sort()
    return intervals


def maxwell_epochs_to_segment_local_intervals(*, maxwell_epochs: list[dict], segment_index: int) -> list[tuple[int, int]]:
    intervals: list[tuple[int, int]] = []
    for e in maxwell_epochs:
        try:
            if int(e.get("segment_index")) != int(segment_index):
                continue
            start = int(e["segment_start_sample"])
            end = int(e["segment_end_sample"])
        except Exception:
            continue
        if end > start:
            intervals.append((start, end))

    intervals.sort()
    return intervals


def filter_spike_train_by_intervals(
    *,
    spike_train: list[int],
    intervals: list[tuple[int, int]],
    pre_samples: int,
    post_samples: int,
) -> tuple[list[int], int, int, list[int], list[int]]:
    if not intervals:
        return spike_train, 0, 0, [], []

    kept: list[int] = []
    removed_outside = 0
    removed_edge = 0
    removed_outside_spikes: list[int] = []
    removed_edge_spikes: list[int] = []

    i = 0
    for t in spike_train:
        t_int = int(t)
        t0 = t_int - int(pre_samples)
        t1 = t_int + int(post_samples)

        while i < len(intervals) and intervals[i][1] <= t_int:
            i += 1

        in_interval = False
        ok = False
        if i < len(intervals):
            start, end = intervals[i]
            if start <= t_int < end:
                in_interval = True
                if t0 >= start and t1 < end:
                    ok = True

        if ok:
            kept.append(t_int)
        else:
            if in_interval:
                removed_edge += 1
                removed_edge_spikes.append(t_int)
            else:
                removed_outside += 1
                removed_outside_spikes.append(t_int)

    return kept, removed_outside, removed_edge, removed_outside_spikes, removed_edge_spikes


__all__ = [
    "WaveformMultisegPlan",
    "WaveformSortingSource",
    "build_waveform_multiseg_plan",
    "resolve_waveform_sorting_source_dir",
    "load_sorting_from_output_dir",
    "load_preprocessed_recording_from_output_dir",
    "resolve_waveform_epoch_marker_paths",
    "harmonize_waveform_sorting",
    "epochs_to_intervals",
    "maxwell_epochs_to_segment_local_intervals",
    "filter_spike_train_by_intervals",
]
