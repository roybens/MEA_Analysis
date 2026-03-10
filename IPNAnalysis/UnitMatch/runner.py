from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any


@dataclass
class UnitMatchConfig:
    enabled: bool = False
    dry_run: bool = True
    fail_open: bool = True


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def run_unitmatch_merge_if_enabled(
    *,
    sorting: Any,
    output_dir: str | Path,
    logger: Any,
    config: UnitMatchConfig,
) -> tuple[Any, dict[str, Any]]:
    """Phase-0 UnitMatch integration hook.

    This hook keeps pipeline wiring stable while we implement full DeepUnitMatch
    adapters and merge operations in later phases.
    """
    out_dir = Path(output_dir) / "unitmatch"
    summary_path = out_dir / "unitmatch_summary.json"

    summary: dict[str, Any] = {
        "timestamp": str(datetime.now()),
        "config": asdict(config),
        "status": "skipped",
        "message": "UnitMatch disabled",
        "n_units_before": None,
        "n_units_after": None,
        "n_candidate_pairs": 0,
        "n_accepted_pairs": 0,
    }

    if not bool(config.enabled):
        _write_json(summary_path, summary)
        return sorting, summary

    try:
        n_before = int(sorting.get_num_units()) if hasattr(sorting, "get_num_units") else None
        summary["n_units_before"] = n_before

        # Phase-0 behavior: wiring + reporting only.
        # We intentionally avoid applying merges until DeepUnitMatch input/output
        # conversion and validation are implemented.
        summary["status"] = "dry_run_only" if bool(config.dry_run) else "not_implemented"
        summary["message"] = (
            "UnitMatch integration wired at analyzer stage; merge execution pending implementation"
        )
        summary["n_units_after"] = n_before

        logger.info("UnitMatch hook active (%s); no merges applied in this phase", summary["status"])
        _write_json(summary_path, summary)
        return sorting, summary
    except Exception as exc:
        summary["status"] = "error"
        summary["message"] = str(exc)
        _write_json(summary_path, summary)
        if bool(config.fail_open):
            logger.warning("UnitMatch failed-open: %s", exc)
            return sorting, summary
        raise
