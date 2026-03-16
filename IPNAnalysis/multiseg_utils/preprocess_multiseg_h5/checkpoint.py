from __future__ import annotations

import json
import os
import re
import tempfile
import traceback
from dataclasses import dataclass
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any, Optional


class ProcessingStage(Enum):
    """Checkpoint stage values aligned with MEA_Analysis pipeline stages."""

    NOT_STARTED = 0
    PREPROCESSING = 1
    PREPROCESSING_COMPLETE = 2
    SORTING = 3
    SORTING_COMPLETE = 4
    ANALYZER = 5
    ANALYZER_COMPLETE = 6
    REPORTS = 7
    REPORTS_COMPLETE = 8


def _now_str() -> str:
    return datetime.now(timezone.utc).isoformat()


def _safe_write_json(path: Path, payload: dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        Path(tmp_name).replace(path)
    finally:
        try:
            Path(tmp_name).unlink(missing_ok=True)
        except Exception:
            pass


def parse_mea_style_metadata(file_path: Path, *, stream_id: str) -> dict[str, Optional[str]]:
    """Best-effort metadata parsing compatible with MEA_Analysis path conventions."""

    file_path = Path(file_path).expanduser().resolve()
    meta: dict[str, Optional[str]] = {
        "run_id": None,
        "chip_id": None,
        "project": None,
        "relative_pattern": f"{file_path.parent.parent.name}/{file_path.parent.name}/{file_path.name}",
        "date": None,
        "well": stream_id,
    }

    try:
        path_str = str(file_path)
        match = re.search(r"/(\d+)/data\.raw\.h5", path_str)
        if match:
            meta["run_id"] = match.group(1)

        parts = path_str.split(os.sep)
        if len(parts) > 5:
            meta["relative_pattern"] = os.path.join(*parts[-6:-1])
            meta["project"] = parts[-6]
            meta["date"] = parts[-5]
            meta["chip_id"] = parts[-4]
    except Exception:
        pass

    meta_file = file_path.parent / ".metadata"
    if meta_file.exists():
        try:
            import configparser

            cfg = configparser.ConfigParser()
            cfg.read(meta_file, encoding="utf-8")
            if "properties" in cfg:
                meta["run_id"] = cfg["properties"].get("runid", meta.get("run_id"))
                meta["project"] = cfg["properties"].get("project_title", meta.get("project"))
            if "runtime" in cfg:
                meta["chip_id"] = cfg["runtime"].get("chipid", meta.get("chip_id"))
        except Exception:
            pass

    return meta


def compute_checkpoint_file(
    *,
    output_dir: Path,
    file_path: Path,
    stream_id: str,
    checkpoint_root: Optional[Path] = None,
) -> Path:
    """Compute the canonical stage checkpoint file under the well output root."""

    meta = parse_mea_style_metadata(file_path, stream_id=stream_id)
    project = meta.get("project") or "UnknownProject"
    run_id = meta.get("run_id") or "UnknownRun"

    ckpt_root = Path(checkpoint_root) if checkpoint_root else Path(output_dir) / "checkpoints"
    ckpt_root.mkdir(parents=True, exist_ok=True)
    return ckpt_root / f"{project}_{run_id}_{stream_id}_checkpoint.json"


@dataclass
class CheckpointState:
    stage: int
    failed_stage: Optional[str]
    last_updated: Optional[str]
    run_id: str
    chip_id: str
    well: str
    project: str
    date: str
    output_dir: str
    error: Optional[dict[str, Any]]
    extras: dict[str, Any]

    @classmethod
    def fresh(cls, *, output_dir: Path, file_path: Path, stream_id: str) -> "CheckpointState":
        meta = parse_mea_style_metadata(file_path, stream_id=stream_id)
        return cls(
            stage=ProcessingStage.NOT_STARTED.value,
            failed_stage=None,
            last_updated=None,
            run_id=(meta.get("run_id") or "UnknownRun"),
            chip_id=(meta.get("chip_id") or "UnknownChip"),
            well=(meta.get("well") or stream_id),
            project=(meta.get("project") or "UnknownProject"),
            date=(meta.get("date") or "UnknownDate"),
            output_dir=str(Path(output_dir)),
            error=None,
            extras={},
        )

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "CheckpointState":
        known = {
            "stage",
            "failed_stage",
            "last_updated",
            "run_id",
            "chip_id",
            "well",
            "project",
            "date",
            "output_dir",
            "error",
        }
        return cls(
            stage=int(payload.get("stage", ProcessingStage.NOT_STARTED.value)),
            failed_stage=payload.get("failed_stage"),
            last_updated=payload.get("last_updated"),
            run_id=str(payload.get("run_id", "UnknownRun")),
            chip_id=str(payload.get("chip_id", "UnknownChip")),
            well=str(payload.get("well", "UnknownWell")),
            project=str(payload.get("project", "UnknownProject")),
            date=str(payload.get("date", "UnknownDate")),
            output_dir=str(payload.get("output_dir", "")),
            error=payload.get("error"),
            extras={k: v for k, v in payload.items() if k not in known},
        )

    def to_dict(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "stage": int(self.stage),
            "failed_stage": self.failed_stage,
            "last_updated": self.last_updated,
            "run_id": self.run_id,
            "chip_id": self.chip_id,
            "well": self.well,
            "project": self.project,
            "date": self.date,
            "output_dir": self.output_dir,
            "error": self.error,
        }
        payload.update(self.extras)
        return payload


def load_checkpoint(
    *,
    checkpoint_file: Path,
    force_restart: bool,
    output_dir: Path,
    file_path: Path,
    stream_id: str,
) -> CheckpointState:
    checkpoint_file = Path(checkpoint_file)
    if checkpoint_file.exists() and not force_restart:
        with open(checkpoint_file, "r", encoding="utf-8") as f:
            payload = json.load(f)
        return CheckpointState.from_dict(payload)

    return CheckpointState.fresh(output_dir=output_dir, file_path=file_path, stream_id=stream_id)


def save_checkpoint(
    *,
    checkpoint_file: Path,
    state: CheckpointState,
    stage: ProcessingStage,
    failed_stage: Optional[str] = None,
    error: Optional[dict[str, Any]] = None,
    extra_fields: Optional[dict[str, Any]] = None,
) -> CheckpointState:
    payload = state.to_dict()
    payload["stage"] = stage.value
    payload["failed_stage"] = failed_stage
    payload["error"] = error
    payload["last_updated"] = _now_str()
    if extra_fields:
        payload.update(extra_fields)

    _safe_write_json(Path(checkpoint_file), payload)
    return CheckpointState.from_dict(payload)


def exception_to_error_dict(exc: BaseException) -> dict[str, Any]:
    return {
        "type": type(exc).__name__,
        "message": str(exc),
        "traceback": "".join(traceback.format_exception(type(exc), exc, exc.__traceback__)),
    }


__all__ = [
    "ProcessingStage",
    "CheckpointState",
    "parse_mea_style_metadata",
    "compute_checkpoint_file",
    "load_checkpoint",
    "save_checkpoint",
    "exception_to_error_dict",
]
