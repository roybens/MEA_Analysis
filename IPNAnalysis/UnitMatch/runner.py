from __future__ import annotations

import csv
import importlib.util
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
    max_candidate_pairs: int = 20000


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def _normalize_unit_id(uid: Any) -> str:
    try:
        return str(int(uid))
    except Exception:
        return str(uid)


def _extract_unit_ids(sorting: Any) -> list[str]:
    if not hasattr(sorting, "get_unit_ids"):
        return []
    raw = list(sorting.get_unit_ids())
    return [_normalize_unit_id(u) for u in raw]


def _build_candidate_pairs(unit_ids: list[str], max_pairs: int) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    n = len(unit_ids)
    for i in range(n):
        for j in range(i + 1, n):
            pairs.append((unit_ids[i], unit_ids[j]))
            if len(pairs) >= max_pairs:
                return pairs
    return pairs


def _write_candidates_csv(path: Path, pairs: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "unit_id_a",
            "unit_id_b",
            "candidate_rank",
            "placeholder_score",
            "source",
            "note",
        ])
        for idx, (a, b) in enumerate(pairs, start=1):
            w.writerow([a, b, idx, "", "sorting_all_pairs", "phase1_dry_run_placeholder"])


def run_unitmatch_merge_if_enabled(
    *,
    sorting: Any,
    output_dir: str | Path,
    logger: Any,
    config: UnitMatchConfig,
) -> tuple[Any, dict[str, Any]]:
    """UnitMatch integration hook.

    Phase-1 behavior supports dry-run candidate generation directly from an
    existing sorting object and records DeepUnitMatch import readiness.
    """
    out_dir = Path(output_dir) / "unitmatch"
    summary_path = out_dir / "unitmatch_summary.json"
    config_path = out_dir / "unitmatch_config.json"
    candidates_path = out_dir / "match_candidates.csv"

    summary: dict[str, Any] = {
        "timestamp": str(datetime.now()),
        "config": asdict(config),
        "status": "skipped",
        "message": "UnitMatch disabled",
        "n_units_before": None,
        "n_units_after": None,
        "n_candidate_pairs": 0,
        "n_accepted_pairs": 0,
        "deepunitmatch_import_ok": False,
        "deepunitmatch_module": None,
        "artifacts": {
            "summary_json": str(summary_path),
            "config_json": str(config_path),
            "candidates_csv": None,
        },
    }

    _write_json(config_path, asdict(config))

    if not bool(config.enabled):
        _write_json(summary_path, summary)
        return sorting, summary

    try:
        n_before = int(sorting.get_num_units()) if hasattr(sorting, "get_num_units") else None
        summary["n_units_before"] = n_before
        unit_ids = _extract_unit_ids(sorting)
        summary["n_units_detected"] = len(unit_ids)

        # Probe import availability so dry-runs reveal whether DeepUnitMatch is
        # actually importable in the active environment.
        dm_spec = importlib.util.find_spec("DeepUnitMatch")
        if dm_spec is not None:
            summary["deepunitmatch_import_ok"] = True
            summary["deepunitmatch_module"] = "DeepUnitMatch"
        else:
            dm_spec = importlib.util.find_spec("UnitMatchPy")
            summary["deepunitmatch_import_ok"] = False
            summary["deepunitmatch_module"] = ("UnitMatchPy" if dm_spec is not None else None)

        if bool(config.dry_run):
            pairs = _build_candidate_pairs(unit_ids, int(max(0, config.max_candidate_pairs)))
            _write_candidates_csv(candidates_path, pairs)
            summary["artifacts"]["candidates_csv"] = str(candidates_path)
            summary["n_candidate_pairs"] = len(pairs)
            summary["status"] = "dry_run_candidates"
            summary["message"] = (
                "Dry-run generated candidate pairs from sorting unit IDs; DeepUnitMatch scoring adapter pending"
            )
            logger.info(
                "UnitMatch dry-run: units=%d candidate_pairs=%d deepunitmatch_import_ok=%s",
                len(unit_ids),
                len(pairs),
                summary["deepunitmatch_import_ok"],
            )
        else:
            summary["status"] = "not_implemented"
            summary["message"] = (
                "Merge execution not implemented yet; enable dry-run to generate candidate artifacts"
            )

        summary["n_units_after"] = n_before
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
