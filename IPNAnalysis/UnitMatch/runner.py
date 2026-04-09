from __future__ import annotations

import csv
import json
import math
import re
import traceback
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from .backends.deepunitmatch_adapter import (
    check_backend_availability,
    infer_similarity_matrix,
    preprocess_waveforms,
)


UNITMATCH_CLONE_ROOT = Path("/home/adamm/dev/pkgs/UnitMatch/UnitMatchPy")
DEEPUNITMATCH_CLONE_DIR = UNITMATCH_CLONE_ROOT / "DeepUnitMatch"


@dataclass
class UnitMatchConfig:
    enabled: bool = False
    dry_run: bool = True
    scored_dry_run: bool = True
    fail_open: bool = True
    max_candidate_pairs: int = 20000
    output_subdir_name: str = "unitmatch_outputs"
    throughput_subdir_name: str = "unitmatch_throughput"
    oversplit_min_probability: float = 0.90
    oversplit_max_suggestions: int = 2000
    apply_merges: bool = False
    recursive: bool = False
    max_iterations: int = 5
    max_spikes_per_unit: int = 100
    keep_all_iterations: bool = True
    backend_name: str = "deepunitmatch_clone"
    backend_device: str = "cpu"
    backend_clone_root: str | None = None


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
    if int(max_pairs) == 0:
        return []

    unlimited = int(max_pairs) < 0
    pairs: list[tuple[str, str]] = []
    n = len(unit_ids)
    for i in range(n):
        for j in range(i + 1, n):
            pairs.append((unit_ids[i], unit_ids[j]))
            if (not unlimited) and len(pairs) >= int(max_pairs):
                return pairs
    return pairs


def _limit_to_summary(limit_value: int) -> dict[str, Any]:
    if int(limit_value) < 0:
        return {"mode": "unlimited", "value": int(limit_value)}
    if int(limit_value) == 0:
        return {"mode": "none", "value": 0}
    return {"mode": "bounded", "value": int(limit_value)}


def _write_candidates_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "unit_id_a",
            "unit_id_b",
            "candidate_rank",
            "score",
            "reverse_score",
            "reciprocal_min_score",
            "score_kind",
            "source",
            "score_error",
            "note",
        ])
        for row in rows:
            w.writerow([
                str(row.get("unit_id_a", "")),
                str(row.get("unit_id_b", "")),
                int(row.get("candidate_rank", 0)),
                row.get("score", ""),
                row.get("reverse_score", ""),
                row.get("reciprocal_min_score", ""),
                str(row.get("score_kind", "")),
                str(row.get("source", "")),
                str(row.get("score_error", "")),
                str(row.get("note", "")),
            ])


def _write_thresholded_csv(path: Path, rows: list[dict[str, Any]], threshold: float) -> dict[str, int]:
    path.parent.mkdir(parents=True, exist_ok=True)
    n_good = 0
    n_bad = 0
    n_unscored = 0
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "unit_id_a",
            "unit_id_b",
            "candidate_rank",
            "score",
            "reverse_score",
            "reciprocal_min_score",
            "threshold",
            "is_good_match",
            "source",
            "score_kind",
            "note",
        ])
        for row in rows:
            score = row.get("score", "")
            if isinstance(score, (int, float)) and not math.isnan(float(score)):
                is_good_bool = float(score) >= float(threshold)
                is_good = "1" if is_good_bool else "0"
                if is_good_bool:
                    n_good += 1
                else:
                    n_bad += 1
            else:
                is_good = ""
                n_unscored += 1

            w.writerow([
                str(row.get("unit_id_a", "")),
                str(row.get("unit_id_b", "")),
                int(row.get("candidate_rank", 0)),
                score,
                row.get("reverse_score", ""),
                row.get("reciprocal_min_score", ""),
                float(threshold),
                is_good,
                str(row.get("source", "")),
                str(row.get("score_kind", "")),
                str(row.get("note", "")),
            ])

    return {
        "n_good_matches": int(n_good),
        "n_bad_matches": int(n_bad),
        "n_unscored": int(n_unscored),
    }


def _write_oversplit_suggestions_csv(path: Path, suggestions: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "unit_id_a",
            "unit_id_b",
            "score",
            "reverse_score",
            "reciprocal_min_score",
            "score_kind",
            "source",
            "recommended_action",
            "candidate_rank",
        ])
        for s in suggestions:
            w.writerow([
                str(s.get("unit_id_a", "")),
                str(s.get("unit_id_b", "")),
                float(s.get("score", float("nan"))),
                s.get("reverse_score", ""),
                s.get("reciprocal_min_score", ""),
                str(s.get("score_kind", "")),
                str(s.get("source", "")),
                str(s.get("recommended_action", "suggest_merge")),
                int(s.get("candidate_rank", 0)),
            ])


def _placeholder_rows_from_pairs(pairs: list[tuple[str, str]], note: str, source: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for idx, (a, b) in enumerate(pairs, start=1):
        rows.append(
            {
                "unit_id_a": a,
                "unit_id_b": b,
                "candidate_rank": idx,
                "score": "",
                "score_kind": "",
                "source": source,
                "score_error": "",
                "note": note,
            }
        )
    return rows


def _rows_from_pair_scores(
    *,
    pairs: list[tuple[str, str]],
    scores: list[float],
    reverse_scores: list[float] | None,
    source: str,
    score_kind: str,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    reverse_list = reverse_scores if reverse_scores is not None else [float("nan")] * len(scores)
    for idx, ((a, b), score, reverse_score) in enumerate(zip(pairs, scores, reverse_list), start=1):
        reciprocal_min = float("nan")
        if isinstance(score, (int, float)) and isinstance(reverse_score, (int, float)):
            score_f = float(score)
            reverse_f = float(reverse_score)
            if (not math.isnan(score_f)) and (not math.isnan(reverse_f)):
                reciprocal_min = min(score_f, reverse_f)
        rows.append(
            {
                "unit_id_a": a,
                "unit_id_b": b,
                "candidate_rank": idx,
                "score": float(score),
                "reverse_score": float(reverse_score) if isinstance(reverse_score, (int, float)) else "",
                "reciprocal_min_score": reciprocal_min,
                "score_kind": score_kind,
                "source": source,
                "score_error": "",
                "note": "deepunitmatch_scored_dry_run",
            }
        )
    return rows


def _candidate_ks_dirs(output_dir: str | Path) -> list[Path]:
    base = Path(output_dir)
    return [
        base / "sorter_output",
        base / "sorter_output" / "sorter_output",
        base / "sorter_output" / "in_container_sorting",
    ]


def _raw_waveforms_dir_for_ks(ks_dir: Path) -> Path | None:
    candidates = [
        ks_dir / "RawWaveforms",
        ks_dir / "qMetrics" / "RawWaveforms",
        ks_dir / "bombcell" / "RawWaveforms",
    ]
    for c in candidates:
        if c.exists() and c.is_dir():
            return c
    return None


def _choose_ks_dir(output_dir: str | Path) -> tuple[Path | None, dict[str, Any]]:
    candidates = _candidate_ks_dirs(output_dir)
    searched = [str(c) for c in candidates]

    for c in candidates:
        if c.exists() and c.is_dir() and (c / "channel_positions.npy").exists():
            return c, {"searched_ks_candidates": searched, "selected_ks_dir": str(c)}

    return None, {"searched_ks_candidates": searched, "selected_ks_dir": None}


def _ensure_unitmatch_rawwaveforms_from_si(
    *,
    ks_dir: Path,
    sorting: Any,
    recording: Any,
    raw_waveforms_dir: Path,
    max_spikes_per_unit: int,
    logger: Any,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "attempted": False,
        "ok": False,
        "raw_waveforms_dir": None,
        "message": None,
    }

    if raw_waveforms_dir.exists() and any(raw_waveforms_dir.glob("Unit*_RawSpikes.npy")):
        result["ok"] = True
        result["raw_waveforms_dir"] = str(raw_waveforms_dir)
        result["message"] = "already_present"
        return result

    if sorting is None or recording is None:
        result["message"] = "sorting_or_recording_missing"
        return result

    if not (ks_dir / "channel_positions.npy").exists():
        result["message"] = "channel_positions_missing"
        return result

    result["attempted"] = True
    raw_dir = Path(raw_waveforms_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)

    try:
        import numpy as np
        import spikeinterface.full as si

        # Design decision:
        # DeepUnitMatch has historically assumed Neuropixels-like template shapes,
        # including a time axis convention of 82 samples. We patched DeepUnitMatch
        # for more variable Maxwell HDMEA inputs, but this exporter still targets
        # 82 time points for compatibility and to reduce shape-related regressions.
        # Impact on matching performance is currently unknown and should be
        # validated empirically on held-out data.
        target_fs_hz = 30000.0
        target_n_time = 82
        target_window_ms = float((target_n_time - 1) * 1000.0 / target_fs_hz)
        ms_before = target_window_ms / 3.0
        ms_after = target_window_ms - ms_before

        def _resample_to_target_window(template_tc: Any, source_fs_hz: float) -> Any:
            if template_tc.ndim != 2:
                return template_tc

            n_time, n_channels = int(template_tc.shape[0]), int(template_tc.shape[1])
            if n_time <= 1:
                return np.repeat(template_tc, target_n_time, axis=0)

            # Resample each channel over the full captured duration to exactly 82 samples.
            duration_s = float((n_time - 1) / max(source_fs_hz, 1.0))
            t_old = np.linspace(0.0, duration_s, n_time, dtype=np.float64)
            t_new = np.linspace(0.0, duration_s, target_n_time, dtype=np.float64)

            resampled = np.empty((target_n_time, n_channels), dtype=np.float32)
            for ch in range(n_channels):
                resampled[:, ch] = np.interp(t_new, t_old, template_tc[:, ch]).astype(np.float32)
            return resampled

        n_samples = int(recording.get_num_samples())
        source_fs_hz = float(recording.get_sampling_frequency())
        split_idx = max(1, n_samples // 2)

        sorting_halves = [
            sorting.frame_slice(start_frame=0, end_frame=split_idx),
            sorting.frame_slice(start_frame=split_idx, end_frame=n_samples),
        ]
        recording_halves = [
            recording.frame_slice(start_frame=0, end_frame=split_idx),
            recording.frame_slice(start_frame=split_idx, end_frame=n_samples),
        ]

        template_by_half: list[dict[int, Any]] = []
        for half_idx in range(2):
            ana = si.create_sorting_analyzer(
                sorting_halves[half_idx],
                recording_halves[half_idx],
                sparse=False,
            )
            random_spikes_kwargs: dict[str, Any] = {"method": "uniform"}
            if int(max_spikes_per_unit) >= 0:
                random_spikes_kwargs["max_spikes_per_unit"] = int(max_spikes_per_unit)
            ana.compute("random_spikes", **random_spikes_kwargs)
            ana.compute("waveforms", ms_before=ms_before, ms_after=ms_after, dtype="float32", save=False)
            ana.compute("templates")

            t_ext = ana.get_extension("templates")
            t_data = t_ext.get_data()
            t_unit_ids = [int(u) for u in list(ana.sorting.get_unit_ids())]
            template_by_half.append({uid: t_data[i] for i, uid in enumerate(t_unit_ids)})

        unit_ids_0 = set(template_by_half[0].keys())
        unit_ids_1 = set(template_by_half[1].keys())
        common_ids = sorted(list(unit_ids_0.intersection(unit_ids_1)))
        if len(common_ids) == 0:
            result["message"] = "no_common_units_between_halves"
            return result

        n_saved = 0
        for uid in common_ids:
            cv0 = template_by_half[0][uid]
            cv1 = template_by_half[1][uid]
            cv0 = _resample_to_target_window(cv0, source_fs_hz)
            cv1 = _resample_to_target_window(cv1, source_fs_hz)
            avg = np.stack((cv0, cv1), axis=-1)
            np.save(raw_dir / f"Unit{int(uid)}_RawSpikes.npy", avg)
            n_saved += 1

        if n_saved == 0:
            result["message"] = "no_waveforms_saved"
            return result

        result["ok"] = True
        result["raw_waveforms_dir"] = str(raw_dir)
        result["message"] = f"generated_{n_saved}_units"
        result["generation_mode"] = "split_half_sortinganalyzer_templates"
        result["n_units_common_halves"] = int(len(common_ids))
        result["source_sampling_frequency_hz"] = float(source_fs_hz)
        result["target_sampling_frequency_hz"] = float(target_fs_hz)
        result["target_waveform_n_time"] = int(target_n_time)
        result["waveform_ms_before"] = float(ms_before)
        result["waveform_ms_after"] = float(ms_after)
        result["waveform_window_ms"] = float(target_window_ms)
        return result
    except Exception as exc:
        logger.warning("Failed to auto-generate UnitMatch RawWaveforms in %s: %s", ks_dir, exc)
        result["message"] = str(exc)
        result["traceback"] = traceback.format_exc()
        return result


def _parse_unit_id_from_filename(name: str) -> int | None:
    m = re.match(r"^Unit(\d+).*_RawSpikes\.npy$", name)
    if not m:
        return None
    try:
        return int(m.group(1))
    except Exception:
        return None


def _load_waveforms_from_raw_dir(raw_dir: Path) -> tuple[Any, list[str], dict[str, Any]]:
    import numpy as np

    files = sorted([p for p in raw_dir.glob("Unit*_RawSpikes.npy") if p.is_file()])
    if not files:
        return None, [], {"reason": "no_rawwaveform_files"}

    unit_ids: list[str] = []
    waves: list[Any] = []
    bad_shapes: list[dict[str, Any]] = []
    accepted_shape_histogram: dict[str, int] = {}
    rejection_counts: dict[str, int] = {
        "bad_rank": 0,
        "too_short_time": 0,
        "too_few_channels": 0,
        "bad_repeat_axis": 0,
    }

    min_time = 20
    min_channels = 8

    for fp in files:
        uid = _parse_unit_id_from_filename(fp.name)
        if uid is None:
            continue
        arr = np.load(fp)
        reject_reason = None
        wave = None

        if arr.ndim == 2:
            t, c = arr.shape
            if t < min_time:
                reject_reason = "too_short_time"
            elif c < min_channels:
                reject_reason = "too_few_channels"
            else:
                wave = np.stack((arr, arr), axis=-1)
        elif arr.ndim == 3:
            t, c, cv = arr.shape
            if t < min_time:
                reject_reason = "too_short_time"
            elif c < min_channels:
                reject_reason = "too_few_channels"
            elif cv < 1:
                reject_reason = "bad_repeat_axis"
            elif cv == 1:
                wave = np.stack((arr[:, :, 0], arr[:, :, 0]), axis=-1)
            else:
                wave = arr[:, :, :2]
        else:
            reject_reason = "bad_rank"

        if wave is None:
            if reject_reason is None:
                reject_reason = "bad_rank"
            rejection_counts[reject_reason] = int(rejection_counts.get(reject_reason, 0) + 1)
            bad_shapes.append({"file": fp.name, "shape": list(arr.shape), "reason": reject_reason})
            continue

        unit_ids.append(str(uid))
        waves.append(wave)
        shape_key = str(tuple(int(x) for x in wave.shape))
        accepted_shape_histogram[shape_key] = int(accepted_shape_histogram.get(shape_key, 0) + 1)

    if not waves:
        return None, [], {
            "reason": "no_valid_rawwaveforms",
            "bad_shapes": bad_shapes,
            "rejection_counts": rejection_counts,
            "accepted_shape_histogram": accepted_shape_histogram,
        }

    waveforms = np.stack(waves, axis=0)
    return waveforms, unit_ids, {
        "reason": "ok",
        "n_loaded": len(unit_ids),
        "bad_shapes": bad_shapes,
        "rejection_counts": rejection_counts,
        "accepted_shape_histogram": accepted_shape_histogram,
        "accepted_shape_examples": [list(w.shape) for w in waves[:3]],
    }


def _write_score_matrix_artifacts(
    *,
    sim_matrix: Any,
    matrix_unit_ids: list[str],
    out_dir: Path,
) -> dict[str, str]:
    import numpy as np

    out_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = out_dir / "deepunitmatch_similarity_matrix.npy"
    unit_ids_path = out_dir / "deepunitmatch_similarity_matrix_unit_ids.json"
    np.save(matrix_path, sim_matrix)
    with unit_ids_path.open("w", encoding="utf-8") as f:
        json.dump({"unit_ids": [str(u) for u in matrix_unit_ids]}, f, indent=2)
    return {
        "score_matrix_npy": str(matrix_path),
        "score_matrix_unit_ids_json": str(unit_ids_path),
    }


def _score_stats(scores: list[float]) -> dict[str, Any]:
    import numpy as np

    valid = [float(s) for s in scores if isinstance(s, (int, float)) and not math.isnan(float(s))]
    n_total = int(len(scores))
    n_valid = int(len(valid))
    n_nan = int(n_total - n_valid)
    if n_valid == 0:
        return {
            "n_total": n_total,
            "n_valid": n_valid,
            "n_nan": n_nan,
            "has_valid_scores": False,
        }

    arr = np.asarray(valid, dtype=np.float64)
    return {
        "n_total": n_total,
        "n_valid": n_valid,
        "n_nan": n_nan,
        "has_valid_scores": True,
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "percentiles": {
            "p50": float(np.percentile(arr, 50.0)),
            "p75": float(np.percentile(arr, 75.0)),
            "p90": float(np.percentile(arr, 90.0)),
            "p95": float(np.percentile(arr, 95.0)),
            "p99": float(np.percentile(arr, 99.0)),
        },
    }


def _matrix_directional_scores_for_pairs(
    *,
    pairs: list[tuple[str, str]],
    matrix: Any,
    matrix_unit_ids: list[str],
) -> tuple[list[float], list[float], int]:
    id_to_idx = {uid: idx for idx, uid in enumerate(matrix_unit_ids)}
    scores: list[float] = []
    reverse_scores: list[float] = []
    missing = 0
    for a, b in pairs:
        ia = id_to_idx.get(str(a))
        ib = id_to_idx.get(str(b))
        if ia is None or ib is None:
            missing += 1
            scores.append(float("nan"))
            reverse_scores.append(float("nan"))
            continue
        scores.append(float(matrix[ia, ib]))
        reverse_scores.append(float(matrix[ib, ia]))
    return scores, reverse_scores, missing


def _score_pairs_deepunitmatch_clone_only(
    *,
    pairs: list[tuple[str, str]],
    sorting: Any,
    recording: Any,
    output_dir: str | Path,
    output_subdir_name: str,
    throughput_iteration_dir: Path,
    logger: Any,
    backend_device: str,
    backend_clone_root: str | None,
    max_spikes_per_unit: int,
) -> tuple[list[float] | None, list[float] | None, str, dict[str, Any]]:
    clone_root = Path(backend_clone_root) if backend_clone_root else UNITMATCH_CLONE_ROOT
    clone_dir = clone_root / "DeepUnitMatch"
    availability = check_backend_availability(clone_root)
    probe: dict[str, Any] = {
        "module": "DeepUnitMatch",
        "strategy": "deepunitmatch_clone_inference",
        "clone_root": str(clone_root),
        "clone_dir": str(clone_dir),
        "clone_present": bool(availability.diagnostics.get("clone_present", False)),
        "import_ok": bool(availability.available),
        "reason": None,
    }
    if not availability.available:
        probe["reason"] = str(availability.reason_code)
        probe.update(availability.diagnostics)
        return None, None, "", probe

    ks_dir, ks_diag = _choose_ks_dir(output_dir)
    probe.update(ks_diag)
    if ks_dir is None:
        probe["reason"] = "no_ks_dir_with_channel_positions"
        return None, None, "", probe

    raw_dir = throughput_iteration_dir / "RawWaveforms"
    if not raw_dir.exists():
        gen = _ensure_unitmatch_rawwaveforms_from_si(
            ks_dir=ks_dir,
            sorting=sorting,
            recording=recording,
            raw_waveforms_dir=raw_dir,
            max_spikes_per_unit=int(max_spikes_per_unit),
            logger=logger,
        )
        probe["rawwaveforms_generation"] = gen
        if not bool(gen.get("ok")):
            probe["reason"] = "rawwaveforms_generation_failed"
            return None, None, "", probe

    if not raw_dir.exists():
        probe["reason"] = "throughput_rawwaveforms_missing"
        return None, None, "", probe

    probe["throughput_iteration_dir"] = str(throughput_iteration_dir)
    probe["throughput_rawwaveforms_dir"] = str(raw_dir)

    try:
        import numpy as np

        channel_pos = np.load(ks_dir / "channel_positions.npy")
        waveforms, matrix_unit_ids, load_diag = _load_waveforms_from_raw_dir(raw_dir)
        probe["waveform_load"] = load_diag
        if waveforms is None or len(matrix_unit_ids) == 0:
            probe["reason"] = str(load_diag.get("reason") or "waveform_load_failed")
            return None, None, "", probe

        # Build one-session demo-style snippets into output-local temp folder.
        prep_dir = throughput_iteration_dir / "deepunitmatch_preprocessed"
        processed_dir, prep_diag = preprocess_waveforms(
            waveforms=waveforms,
            channel_pos=channel_pos,
            matrix_unit_ids=matrix_unit_ids,
            prep_dir=prep_dir,
            clone_root=clone_root,
        )
        probe["preprocess_diagnostics"] = prep_diag
        if processed_dir is None:
            probe["reason"] = str(prep_diag.get("reason") or "preprocessing_failed")
            return None, None, "", probe

        infer_result = infer_similarity_matrix(
            processed_dir=processed_dir,
            clone_root=clone_root,
            device=str(backend_device),
        )
        probe["inference_diagnostics"] = infer_result.diagnostics
        probe["import_ok"] = bool(infer_result.import_ok)
        if not infer_result.success or infer_result.sim_matrix is None:
            probe["reason"] = str(infer_result.reason_code)
            return None, None, "", probe

        sim_matrix = infer_result.sim_matrix
        probe["artifacts"] = _write_score_matrix_artifacts(
            sim_matrix=sim_matrix,
            matrix_unit_ids=matrix_unit_ids,
            out_dir=throughput_iteration_dir,
        )

        scores, reverse_scores, missing = _matrix_directional_scores_for_pairs(
            pairs=pairs,
            matrix=sim_matrix,
            matrix_unit_ids=matrix_unit_ids,
        )

        probe["n_matrix_units"] = int(len(matrix_unit_ids))
        probe["n_missing_pair_mappings"] = int(missing)
        asym_count = 0
        for fwd, rev in zip(scores, reverse_scores):
            if math.isnan(float(fwd)) or math.isnan(float(rev)):
                continue
            if not math.isclose(float(fwd), float(rev), rel_tol=1e-6, abs_tol=1e-8):
                asym_count += 1
        probe["n_asymmetric_pairs"] = int(asym_count)
        probe["reason"] = "ok" if missing < len(pairs) else "no_pair_id_overlap"
        if missing >= len(pairs):
            return None, None, "", probe
        return scores, reverse_scores, "DeepUnitMatch.adapter_clone_inference", probe
    except Exception as exc:
        probe["reason"] = "deepunitmatch_runtime_error"
        probe["message"] = str(exc)
        probe["traceback"] = traceback.format_exc()
        return None, None, "", probe


def _build_oversplit_suggestions(
    *,
    rows: list[dict[str, Any]],
    min_probability: float,
    max_suggestions: int,
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    suggestions: list[dict[str, Any]] = []
    n_rejected_non_numeric = 0
    n_rejected_forward_below = 0
    n_rejected_reciprocal_below = 0
    for row in rows:
        score = row.get("score", "")
        if not isinstance(score, (int, float)):
            n_rejected_non_numeric += 1
            continue
        score_f = float(score)
        if math.isnan(score_f):
            n_rejected_non_numeric += 1
            continue
        if score_f < float(min_probability):
            n_rejected_forward_below += 1
            continue

        reverse = row.get("reverse_score", "")
        if not isinstance(reverse, (int, float)):
            n_rejected_non_numeric += 1
            continue
        reverse_f = float(reverse)
        if math.isnan(reverse_f):
            n_rejected_non_numeric += 1
            continue

        reciprocal_min = min(score_f, reverse_f)
        if reciprocal_min < float(min_probability):
            n_rejected_reciprocal_below += 1
            continue

        a = str(row.get("unit_id_a", ""))
        b = str(row.get("unit_id_b", ""))
        if not a or not b or a == b:
            continue

        suggestions.append(
            {
                "unit_id_a": a,
                "unit_id_b": b,
                "score": score_f,
                "reverse_score": reverse_f,
                "reciprocal_min_score": reciprocal_min,
                "score_kind": str(row.get("score_kind", "")),
                "source": str(row.get("source", "")),
                "recommended_action": "suggest_merge",
                "candidate_rank": int(row.get("candidate_rank", 0)),
            }
        )

    suggestions.sort(key=lambda x: float(x["score"]), reverse=True)
    max_s = int(max_suggestions)
    if max_s < 0:
        return suggestions, {
            "n_rejected_non_numeric": int(n_rejected_non_numeric),
            "n_rejected_forward_below_threshold": int(n_rejected_forward_below),
            "n_rejected_reciprocal_below_threshold": int(n_rejected_reciprocal_below),
        }
    if max_s == 0:
        return [], {
            "n_rejected_non_numeric": int(n_rejected_non_numeric),
            "n_rejected_forward_below_threshold": int(n_rejected_forward_below),
            "n_rejected_reciprocal_below_threshold": int(n_rejected_reciprocal_below),
        }
    return suggestions[:max_s], {
        "n_rejected_non_numeric": int(n_rejected_non_numeric),
        "n_rejected_forward_below_threshold": int(n_rejected_forward_below),
        "n_rejected_reciprocal_below_threshold": int(n_rejected_reciprocal_below),
    }


def _select_conflict_free_top_pairs(suggestions: list[dict[str, Any]]) -> list[dict[str, Any]]:
    ordered = sorted(
        suggestions,
        key=lambda x: (float(x.get("score", float("nan"))), -int(x.get("candidate_rank", 0))),
        reverse=True,
    )
    selected: list[dict[str, Any]] = []
    used_units: set[str] = set()
    for s in ordered:
        a = str(s.get("unit_id_a", ""))
        b = str(s.get("unit_id_b", ""))
        if not a or not b or a == b:
            continue
        if a in used_units or b in used_units:
            continue
        selected.append(s)
        used_units.add(a)
        used_units.add(b)
    return selected


def _suggestion_partner_summary(suggestions: list[dict[str, Any]]) -> dict[str, Any]:
    by_unit: dict[str, list[dict[str, Any]]] = {}
    for s in suggestions:
        a = str(s.get("unit_id_a", ""))
        b = str(s.get("unit_id_b", ""))
        if not a or not b:
            continue
        by_unit.setdefault(a, []).append({"partner": b, "score": float(s.get("score", float("nan")))})
        by_unit.setdefault(b, []).append({"partner": a, "score": float(s.get("score", float("nan")))})

    counts = {u: int(len(v)) for u, v in by_unit.items()}
    top_partner: dict[str, Any] = {}
    for u, partners in by_unit.items():
        partners_sorted = sorted(partners, key=lambda x: float(x["score"]), reverse=True)
        top_partner[u] = partners_sorted[0] if partners_sorted else None

    return {
        "n_units_with_suggestions": int(len(by_unit)),
        "partner_count_by_unit": counts,
        "top_partner_by_unit": top_partner,
    }


def _write_per_unit_partner_summary(path: Path, summary: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)


def _read_selected_pairs_csv(path: Path) -> list[tuple[str, str]]:
    if not path.exists() or not path.is_file():
        return []

    pairs: list[tuple[str, str]] = []
    try:
        with path.open("r", encoding="utf-8", newline="") as f:
            r = csv.DictReader(f)
            for row in r:
                a = str(row.get("unit_id_a", ""))
                b = str(row.get("unit_id_b", ""))
                if a and b and a != b:
                    pairs.append((a, b))
    except Exception:
        return []
    return pairs


def _apply_selected_merges(
    *,
    sorting: Any,
    selected_suggestions: list[dict[str, Any]],
) -> tuple[Any, dict[str, Any]]:
    try:
        from spikeinterface.core import apply_merges_to_sorting
    except Exception as exc:
        return sorting, {
            "attempted": False,
            "applied": False,
            "reason": f"spikeinterface_merge_import_failed: {exc}",
            "n_requested_groups": 0,
            "n_applied_groups": 0,
        }

    sorting_ids = list(getattr(sorting, "get_unit_ids")()) if hasattr(sorting, "get_unit_ids") else []
    id_lookup = {_normalize_unit_id(uid): uid for uid in sorting_ids}
    merge_groups: list[list[Any]] = []
    skipped = 0
    for s in selected_suggestions:
        a_norm = str(s.get("unit_id_a", ""))
        b_norm = str(s.get("unit_id_b", ""))
        a = id_lookup.get(a_norm)
        b = id_lookup.get(b_norm)
        if a is None or b is None or a == b:
            skipped += 1
            continue
        merge_groups.append([a, b])

    if not merge_groups:
        return sorting, {
            "attempted": True,
            "applied": False,
            "reason": "no_valid_merge_groups",
            "n_requested_groups": int(len(selected_suggestions)),
            "n_applied_groups": 0,
            "n_skipped_groups": int(skipped),
        }

    try:
        merged = apply_merges_to_sorting(sorting, merge_groups, censor_ms=None)
        return merged, {
            "attempted": True,
            "applied": True,
            "reason": "ok",
            "n_requested_groups": int(len(selected_suggestions)),
            "n_applied_groups": int(len(merge_groups)),
            "n_skipped_groups": int(skipped),
        }
    except Exception as exc:
        return sorting, {
            "attempted": True,
            "applied": False,
            "reason": f"apply_merges_failed: {exc}",
            "n_requested_groups": int(len(selected_suggestions)),
            "n_applied_groups": 0,
            "n_skipped_groups": int(skipped),
        }


def run_unitmatch_merge_if_enabled(
    *,
    sorting: Any,
    recording: Any = None,
    output_dir: str | Path,
    logger: Any,
    config: UnitMatchConfig,
    iteration_index: int = 0,
) -> tuple[Any, dict[str, Any]]:
    clone_root_cfg = Path(str(config.backend_clone_root)) if config.backend_clone_root else UNITMATCH_CLONE_ROOT
    clone_dir_cfg = clone_root_cfg / "DeepUnitMatch"
    out_dir = Path(output_dir) / str(config.output_subdir_name)
    throughput_root = Path(output_dir) / str(config.throughput_subdir_name)
    iteration_dir = throughput_root / f"iteration_{int(max(0, iteration_index))}"
    summary_path = out_dir / "unitmatch_summary.json"
    config_path = out_dir / "unitmatch_config.json"
    candidates_path = out_dir / "match_candidates.csv"
    thresholded_path = out_dir / "match_candidates_thresholded.csv"
    oversplit_path = out_dir / "oversplit_suggestions.csv"

    summary: dict[str, Any] = {
        "timestamp": str(datetime.now()),
        "config": asdict(config),
        "status": "disabled",
        "message": "UnitMatch disabled",
        "n_units_before": None,
        "n_units_after": None,
        "n_candidate_pairs": 0,
        "n_accepted_pairs": 0,
        "n_good_matches": 0,
        "n_bad_matches": 0,
        "n_unscored_matches": 0,
        "n_oversplit_candidates": 0,
        "n_oversplit_suggestions": 0,
        "n_scored_pairs": 0,
        "n_selected_non_conflicting_pairs": 0,
        "clone_import": {
            "clone_root": str(clone_root_cfg),
            "clone_dir": str(clone_dir_cfg),
            "clone_present": bool(clone_root_cfg.exists() and clone_dir_cfg.exists()),
            "import_ok": False,
            "failure_reason": None,
        },
        "artifacts": {
            "summary_json": str(summary_path),
            "config_json": str(config_path),
            "candidates_csv": None,
            "thresholded_csv": None,
            "oversplit_suggestions_csv": None,
            "throughput_iteration_dir": str(iteration_dir),
            "per_unit_partner_summary_json": None,
            "final_merged_sorting_folder": None,
        },
        "limit_policy": {
            "max_candidate_pairs": _limit_to_summary(int(config.max_candidate_pairs)),
            "oversplit_max_suggestions": _limit_to_summary(int(config.oversplit_max_suggestions)),
        },
        "iteration": {
            "index": int(max(0, iteration_index)),
            "throughput_root": str(throughput_root),
            "throughput_iteration_dir": str(iteration_dir),
            "convergence": {
                "status": "single_pass",
                "stop_reason": "not_recursive_in_runner",
            },
        },
    }

    _write_json(config_path, asdict(config))

    if not bool(config.enabled):
        _write_json(summary_path, summary)
        return sorting, summary

    iteration_dir.mkdir(parents=True, exist_ok=True)

    try:
        n_before = int(sorting.get_num_units()) if hasattr(sorting, "get_num_units") else None
        summary["n_units_before"] = n_before
        unit_ids = _extract_unit_ids(sorting)
        summary["n_units_detected"] = len(unit_ids)

        pairs = _build_candidate_pairs(unit_ids, int(config.max_candidate_pairs))
        rows: list[dict[str, Any]]

        if not bool(config.scored_dry_run):
            rows = _placeholder_rows_from_pairs(
                pairs,
                note="dry_run_candidates_only",
                source="sorting_all_pairs",
            )
            summary["status"] = "candidates_only"
            summary["message"] = "Scoring disabled; wrote unscored candidate pairs"
        else:
            scores, reverse_scores, backend, probe = _score_pairs_deepunitmatch_clone_only(
                pairs=pairs,
                sorting=sorting,
                recording=recording,
                output_dir=output_dir,
                output_subdir_name=str(config.output_subdir_name),
                throughput_iteration_dir=iteration_dir,
                logger=logger,
                backend_device=str(config.backend_device),
                backend_clone_root=config.backend_clone_root,
                max_spikes_per_unit=int(config.max_spikes_per_unit),
            )
            summary["clone_import"]["import_ok"] = bool(probe.get("import_ok", False))
            summary["clone_import"]["failure_reason"] = probe.get("reason")
            summary["deepunitmatch_scoring_diagnostics"] = probe

            if scores is None:
                rows = _placeholder_rows_from_pairs(
                    pairs,
                    note="scored_dry_run_fallback_unscored",
                    source="deepunitmatch_unavailable",
                )
                summary["status"] = "scored_fallback_unscored"
                summary["message"] = (
                    "Scored dry-run failed; wrote unscored candidate pairs "
                    f"(reason={probe.get('reason')})"
                )
            else:
                if not isinstance(reverse_scores, list) or len(reverse_scores) != len(scores):
                    reverse_scores = [float("nan")] * len(scores)
                rows = _rows_from_pair_scores(
                    pairs=pairs,
                    scores=scores,
                    reverse_scores=reverse_scores,
                    source="DeepUnitMatch",
                    score_kind="similarity_score",
                )
                summary["n_scored_pairs"] = int(len(rows))
                summary["score_statistics"] = _score_stats(scores)
                summary["status"] = "scored_ok"
                summary["message"] = "Scored dry-run completed with DeepUnitMatch clone inference"

        _write_candidates_csv(candidates_path, rows)
        summary["artifacts"]["candidates_csv"] = str(candidates_path)
        summary["n_candidate_pairs"] = len(rows)

        threshold = float(config.oversplit_min_probability)
        thresholded_counts = _write_thresholded_csv(thresholded_path, rows, threshold)
        summary["artifacts"]["thresholded_csv"] = str(thresholded_path)
        summary["n_good_matches"] = int(thresholded_counts.get("n_good_matches", 0))
        summary["n_bad_matches"] = int(thresholded_counts.get("n_bad_matches", 0))
        summary["n_unscored_matches"] = int(thresholded_counts.get("n_unscored", 0))
        summary["n_accepted_pairs"] = int(summary["n_good_matches"])

        suggestions, reciprocal_gate_stats = _build_oversplit_suggestions(
            rows=rows,
            min_probability=float(config.oversplit_min_probability),
            max_suggestions=int(config.oversplit_max_suggestions),
        )
        selected_for_merge = _select_conflict_free_top_pairs(suggestions)
        summary["n_selected_non_conflicting_pairs"] = int(len(selected_for_merge))

        _write_oversplit_suggestions_csv(oversplit_path, suggestions)
        summary["artifacts"]["oversplit_suggestions_csv"] = str(oversplit_path)
        summary["n_oversplit_candidates"] = int(summary["n_good_matches"])
        summary["n_oversplit_suggestions"] = int(len(suggestions))

        selected_merge_path = iteration_dir / "selected_non_conflicting_pairs.csv"
        _write_oversplit_suggestions_csv(selected_merge_path, selected_for_merge)
        summary["artifacts"]["selected_non_conflicting_pairs_csv"] = str(selected_merge_path)

        partner_summary = _suggestion_partner_summary(suggestions)
        partner_summary_path = iteration_dir / "per_unit_partner_summary.json"
        _write_per_unit_partner_summary(partner_summary_path, partner_summary)
        summary["artifacts"]["per_unit_partner_summary_json"] = str(partner_summary_path)
        summary["per_unit_suggestion_summary"] = {
            "n_units_with_suggestions": int(partner_summary.get("n_units_with_suggestions", 0))
        }
        summary["reciprocal_gate"] = {
            "enabled": True,
            "min_probability": float(config.oversplit_min_probability),
            "stats": reciprocal_gate_stats,
        }

        summary["oversplit_policy"] = {
            "min_probability": float(config.oversplit_min_probability),
            "reciprocal_min_probability": float(config.oversplit_min_probability),
            "max_suggestions": int(config.oversplit_max_suggestions),
            "strategy": "threshold_plus_reciprocal_min",
        }

        merge_result = {
            "attempted": False,
            "applied": False,
            "reason": "dry_run_or_disabled",
            "n_requested_groups": int(len(selected_for_merge)),
            "n_applied_groups": 0,
        }
        merged_sorting = sorting
        if bool(config.apply_merges) and (not bool(config.dry_run)):
            merged_sorting, merge_result = _apply_selected_merges(
                sorting=sorting,
                selected_suggestions=selected_for_merge,
            )

        merge_result["persisted_final_merged_sorting"] = False
        if bool(merge_result.get("applied", False)):
            final_merged_sorting_dir = out_dir / "final_merged_sorting"
            try:
                import shutil

                if final_merged_sorting_dir.exists():
                    shutil.rmtree(final_merged_sorting_dir, ignore_errors=True)

                if not hasattr(merged_sorting, "save"):
                    raise RuntimeError("Merged sorting object has no save(folder=...) method")

                merged_sorting.save(folder=str(final_merged_sorting_dir))
                merge_result["persisted_final_merged_sorting"] = bool(final_merged_sorting_dir.exists())
                if bool(merge_result["persisted_final_merged_sorting"]):
                    summary["artifacts"]["final_merged_sorting_folder"] = str(final_merged_sorting_dir)
                else:
                    merge_result["persist_final_merged_sorting_error"] = (
                        f"Final merged sorting folder not found after save: {final_merged_sorting_dir}"
                    )
            except Exception as persist_exc:
                merge_result["persisted_final_merged_sorting"] = False
                merge_result["persist_final_merged_sorting_error"] = str(persist_exc)
                logger.warning(
                    "Could not persist final merged sorting artifact to %s (%s)",
                    final_merged_sorting_dir,
                    persist_exc,
                    exc_info=True,
                )

        summary["merge_application"] = merge_result

        logger.info(
            "UnitMatch dry-run: units=%d candidate_pairs=%d scored=%s clone_import_ok=%s reason=%s",
            len(unit_ids),
            len(rows),
            bool(summary["status"] == "scored_ok"),
            bool(summary["clone_import"].get("import_ok")),
            summary["clone_import"].get("failure_reason"),
        )

        if hasattr(merged_sorting, "get_num_units"):
            summary["n_units_after"] = int(merged_sorting.get_num_units())
        else:
            summary["n_units_after"] = n_before
        _write_json(summary_path, summary)
        return merged_sorting, summary
    except Exception as exc:
        summary["status"] = "error"
        summary["message"] = str(exc)
        _write_json(summary_path, summary)
        if bool(config.fail_open):
            logger.warning("UnitMatch failed-open: %s", exc)
            return sorting, summary
        raise


def run_unitmatch_merge_with_recursion(
    *,
    sorting: Any,
    recording: Any = None,
    output_dir: str | Path,
    logger: Any,
    config: UnitMatchConfig,
) -> tuple[Any, dict[str, Any]]:
    throughput_root = Path(output_dir) / str(config.throughput_subdir_name)
    throughput_root.mkdir(parents=True, exist_ok=True)

    # Always start recursive runs from a clean iteration workspace.
    # This avoids stale artifacts when rerunning with the same throughput path,
    # even when keep_all_iterations is enabled.
    for child in throughput_root.glob("iteration_*"):
        if child.is_dir():
            import shutil

            shutil.rmtree(child, ignore_errors=True)

    if not bool(config.recursive):
        return run_unitmatch_merge_if_enabled(
            sorting=sorting,
            recording=recording,
            output_dir=output_dir,
            logger=logger,
            config=config,
            iteration_index=0,
        )

    current_sorting = sorting
    iteration = 0
    max_iterations_raw = int(config.max_iterations)
    max_iterations = max(0, int(max_iterations_raw))
    uncapped = bool(max_iterations_raw < 0)
    history: list[dict[str, Any]] = []
    last_signature: tuple[tuple[str, str], ...] | None = None
    final_summary: dict[str, Any] = {}

    while True:
        iter_cfg = UnitMatchConfig(**asdict(config))
        iter_cfg.recursive = False

        current_sorting, iter_summary = run_unitmatch_merge_if_enabled(
            sorting=current_sorting,
            recording=recording,
            output_dir=output_dir,
            logger=logger,
            config=iter_cfg,
            iteration_index=iteration,
        )
        final_summary = iter_summary

        selected_csv_s = str(iter_summary.get("artifacts", {}).get("selected_non_conflicting_pairs_csv") or "")
        selected_pairs = _read_selected_pairs_csv(Path(selected_csv_s)) if selected_csv_s else []
        signature = tuple(sorted(selected_pairs))

        stop_reason = None
        if not bool(config.apply_merges):
            stop_reason = "merge_disabled"
        elif len(selected_pairs) == 0:
            stop_reason = "no_selected_pairs"
        elif last_signature is not None and signature == last_signature:
            stop_reason = "repeated_selected_pairs"
        elif not bool(iter_summary.get("merge_application", {}).get("applied")):
            stop_reason = "merge_not_applied"
        elif (not uncapped) and (iteration + 1 >= max_iterations):
            stop_reason = "max_iterations_reached"

        history.append(
            {
                "iteration": int(iteration),
                "status": str(iter_summary.get("status", "unknown")),
                "n_selected_pairs": int(len(selected_pairs)),
                "n_oversplit_suggestions": int(iter_summary.get("n_oversplit_suggestions", 0) or 0),
                "stop_reason": stop_reason,
            }
        )

        if stop_reason is not None:
            break

        last_signature = signature
        iteration += 1

    convergence = {
        "recursive_enabled": True,
        "uncapped_iterations": bool(uncapped),
        "max_iterations": int(max_iterations_raw),
        "keep_all_iterations": bool(config.keep_all_iterations),
        "n_iterations_executed": int(len(history)),
        "history": history,
    }
    final_summary["iteration_convergence_report"] = convergence

    merge_history_path = throughput_root / "merge_history.json"
    _write_json(merge_history_path, convergence)
    final_summary.setdefault("artifacts", {})["merge_history_json"] = str(merge_history_path)

    if not bool(config.keep_all_iterations):
        keep_idx = int(history[-1]["iteration"]) if history else 0
        for child in throughput_root.glob("iteration_*"):
            try:
                idx = int(str(child.name).split("_")[-1])
            except Exception:
                idx = -1
            if idx >= 0 and idx != keep_idx and child.is_dir():
                import shutil

                shutil.rmtree(child, ignore_errors=True)

    summary_path_s = str(final_summary.get("artifacts", {}).get("summary_json") or "")
    if summary_path_s:
        _write_json(Path(summary_path_s), final_summary)

    return current_sorting, final_summary
