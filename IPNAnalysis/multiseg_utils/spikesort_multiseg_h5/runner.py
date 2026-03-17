from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class SpikeSortMultisegPlan:
    """Planning payload for optional multisegment sorting flow.

    This module is intentionally lightweight in the first implementation slice.
    Stage-1 orchestration can consume this plan as a no-op baseline and evolve it
    as mechanical logic is migrated from axon_reconstructor stage-2.
    """

    enabled: bool
    segment_count: int
    mode: str


@dataclass(frozen=True)
class UnitMatchArtifactPaths:
    output_subdir: Path
    throughput_subdir: Path
    report_subdir: Path
    final_merged_sorting: Path


def _segment_count(recording: Any) -> int:
    try:
        return int(recording.get_num_segments())
    except Exception:
        return 1


def build_spikesort_multiseg_plan(
    recording: Any,
    *,
    mode: str = "none",
) -> SpikeSortMultisegPlan:
    """Return an initial optional multiseg sorting plan.

    Current behavior is deliberately conservative: this function only reports
    whether multisegment-specific handling should be considered. Execution
    remains owned by existing stage-1 code paths until migration is completed.
    """

    segments = _segment_count(recording)
    token = str(mode or "none").strip().lower()
    enabled = token not in {"", "none", "off", "disabled"} and segments > 1
    return SpikeSortMultisegPlan(enabled=enabled, segment_count=segments, mode=token or "none")


def resolve_unitmatch_artifact_paths(
    *,
    output_dir: Path,
    output_subdir_name: str = "unitmatch_outputs",
    throughput_subdir_name: str = "unitmatch_throughput",
    report_subdir_name: str = "unitmatch_reports",
) -> UnitMatchArtifactPaths:
    """Resolve canonical UnitMatch artifact locations under a stage output root."""

    return UnitMatchArtifactPaths(
        output_subdir=Path(output_dir) / str(output_subdir_name),
        throughput_subdir=Path(output_dir) / str(throughput_subdir_name),
        report_subdir=Path(output_dir) / str(report_subdir_name),
        final_merged_sorting=Path(output_dir) / str(output_subdir_name) / "final_merged_sorting",
    )


def sync_unitmatch_summary_json(*, summary: dict[str, Any], summary_path: str | Path | None) -> bool:
    """Best-effort summary file sync after artifact/path updates.

    Returns True when a summary path was provided and write succeeded.
    """

    if summary_path is None:
        return False

    path = Path(summary_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return True
