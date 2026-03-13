from __future__ import annotations

import json
import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@dataclass(frozen=True)
class UnitMatchReportConfig:
    throughput_subdir_name: str = "unitmatch_throughput"
    output_subdir_name: str = "unitmatch_outputs"
    report_subdir_name: str = "unitmatch_reports"
    max_heatmap_units: int = 200


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _safe_float(value: Any) -> float:
    try:
        out = float(value)
        if math.isnan(out):
            return float("nan")
        return out
    except Exception:
        return float("nan")


def _ensure_logger(logger: Any):
    if logger is not None:
        return logger

    class _FallbackLogger:
        def info(self, msg: str, *args: Any) -> None:
            if args:
                msg = msg % args
            print(msg)

        def warning(self, msg: str, *args: Any) -> None:
            if args:
                msg = msg % args
            print(msg)

    return _FallbackLogger()


def _latest_iteration_dir(throughput_root: Path) -> Path | None:
    candidates = [p for p in throughput_root.glob("iteration_*") if p.is_dir()]
    if not candidates:
        return None

    def _idx(p: Path) -> int:
        try:
            return int(p.name.split("_")[-1])
        except Exception:
            return -1

    candidates.sort(key=_idx)
    return candidates[-1]


def _render_score_plots(scores: list[float], out_dir: Path) -> dict[str, Any]:
    valid = [s for s in scores if not math.isnan(s)]
    if not valid:
        return {"ok": False, "reason": "no_valid_scores"}

    arr = np.asarray(valid, dtype=np.float64)
    hist_png = out_dir / "score_histogram.png"
    cdf_png = out_dir / "score_cdf.png"

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.hist(arr, bins=40, color="#1f77b4", alpha=0.85)
    ax.set_title("UnitMatch Score Distribution")
    ax.set_xlabel("score")
    ax.set_ylabel("count")
    fig.tight_layout()
    fig.savefig(hist_png, dpi=180)
    plt.close(fig)

    xs = np.sort(arr)
    ys = np.arange(1, len(xs) + 1, dtype=np.float64) / float(len(xs))
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(xs, ys, color="#2ca02c", linewidth=2)
    ax.set_title("UnitMatch Score CDF")
    ax.set_xlabel("score")
    ax.set_ylabel("cumulative fraction")
    ax.set_ylim(0.0, 1.0)
    fig.tight_layout()
    fig.savefig(cdf_png, dpi=180)
    plt.close(fig)

    stats_csv = out_dir / "score_statistics_table.csv"
    stats = {
        "n": int(len(arr)),
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "p50": float(np.percentile(arr, 50.0)),
        "p75": float(np.percentile(arr, 75.0)),
        "p90": float(np.percentile(arr, 90.0)),
        "p95": float(np.percentile(arr, 95.0)),
        "p99": float(np.percentile(arr, 99.0)),
    }
    pd.DataFrame([stats]).to_csv(stats_csv, index=False)

    return {
        "ok": True,
        "histogram_png": str(hist_png),
        "cdf_png": str(cdf_png),
        "stats_csv": str(stats_csv),
    }


def _render_matrix_heatmap(matrix_path: Path, out_dir: Path, max_units: int) -> dict[str, Any]:
    if not matrix_path.exists():
        return {"ok": False, "reason": "matrix_missing"}

    mat = np.load(matrix_path)
    if mat.ndim != 2:
        return {"ok": False, "reason": "matrix_not_2d"}

    n = int(mat.shape[0])
    n_show = min(n, int(max_units))
    view = mat[:n_show, :n_show]

    heatmap_png = out_dir / "similarity_matrix_heatmap.png"
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(view, aspect="auto", interpolation="nearest", cmap="viridis")
    ax.set_title(f"Similarity Matrix Heatmap (first {n_show}/{n} units)")
    ax.set_xlabel("unit index")
    ax.set_ylabel("unit index")
    fig.colorbar(im, ax=ax, shrink=0.9)
    fig.tight_layout()
    fig.savefig(heatmap_png, dpi=180)
    plt.close(fig)

    return {
        "ok": True,
        "heatmap_png": str(heatmap_png),
        "n_matrix_units": int(n),
        "n_displayed_units": int(n_show),
    }


def _render_partner_reports(summary_path: Path, out_dir: Path) -> dict[str, Any]:
    if not summary_path.exists():
        return {"ok": False, "reason": "partner_summary_missing"}

    payload = _read_json(summary_path)
    counts = payload.get("partner_count_by_unit", {})
    tops = payload.get("top_partner_by_unit", {})

    rows: list[dict[str, Any]] = []
    for unit_id, n_partners in counts.items():
        top = tops.get(str(unit_id)) or {}
        rows.append(
            {
                "unit_id": str(unit_id),
                "n_suggested_partners": int(n_partners),
                "top_partner": str(top.get("partner", "")),
                "top_partner_score": _safe_float(top.get("score", "")),
            }
        )

    if not rows:
        return {"ok": False, "reason": "partner_rows_empty"}

    df = pd.DataFrame(rows)
    df = df.sort_values(by=["n_suggested_partners", "unit_id"], ascending=[False, True])

    csv_path = out_dir / "per_unit_partner_table.csv"
    df.to_csv(csv_path, index=False)

    top_n = min(len(df), 40)
    plot_df = df.head(top_n)
    plot_png = out_dir / "per_unit_partner_counts.png"
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(plot_df["unit_id"], plot_df["n_suggested_partners"], color="#ff7f0e")
    ax.set_title(f"Top {top_n} Units by Suggested Partners")
    ax.set_xlabel("unit_id")
    ax.set_ylabel("n_suggested_partners")
    ax.tick_params(axis="x", labelrotation=90)
    fig.tight_layout()
    fig.savefig(plot_png, dpi=180)
    plt.close(fig)

    return {
        "ok": True,
        "table_csv": str(csv_path),
        "counts_png": str(plot_png),
        "n_units": int(len(df)),
    }


def _render_merge_diagnostics(
    *,
    suggestions_csv: Path | None,
    selected_csv: Path | None,
    out_dir: Path,
) -> dict[str, Any]:
    n_suggestions = 0
    n_selected = 0

    if suggestions_csv is not None and suggestions_csv.exists():
        try:
            n_suggestions = int(len(pd.read_csv(suggestions_csv)))
        except Exception:
            n_suggestions = 0
    if selected_csv is not None and selected_csv.exists():
        try:
            n_selected = int(len(pd.read_csv(selected_csv)))
        except Exception:
            n_selected = 0

    n_rejected = max(0, int(n_suggestions - n_selected))
    csv_path = out_dir / "merge_diagnostics_counts.csv"
    pd.DataFrame(
        [
            {
                "n_suggestions": int(n_suggestions),
                "n_selected_conflict_free": int(n_selected),
                "n_not_selected": int(n_rejected),
            }
        ]
    ).to_csv(csv_path, index=False)

    plot_png = out_dir / "merge_diagnostics_counts.png"
    fig, ax = plt.subplots(figsize=(7, 4.5))
    labels = ["suggestions", "selected", "not_selected"]
    values = [n_suggestions, n_selected, n_rejected]
    ax.bar(labels, values, color=["#1f77b4", "#2ca02c", "#d62728"])
    ax.set_title("Merge Diagnostics")
    ax.set_ylabel("count")
    fig.tight_layout()
    fig.savefig(plot_png, dpi=180)
    plt.close(fig)

    return {
        "ok": True,
        "counts_csv": str(csv_path),
        "counts_png": str(plot_png),
    }


def _render_convergence(convergence_payload: dict[str, Any], out_dir: Path) -> dict[str, Any]:
    history = convergence_payload.get("history", [])
    if not isinstance(history, list) or len(history) == 0:
        csv_path = out_dir / "iteration_convergence.csv"
        pd.DataFrame(
            [
                {
                    "iteration": 0,
                    "n_selected_pairs": 0,
                    "n_oversplit_suggestions": 0,
                    "stop_reason": "not_available",
                }
            ]
        ).to_csv(csv_path, index=False)
        return {"ok": False, "reason": "no_history", "table_csv": str(csv_path)}

    rows: list[dict[str, Any]] = []
    for step in history:
        rows.append(
            {
                "iteration": int(step.get("iteration", 0)),
                "n_selected_pairs": int(step.get("n_selected_pairs", 0)),
                "n_oversplit_suggestions": int(step.get("n_oversplit_suggestions", 0)),
                "stop_reason": str(step.get("stop_reason", "")),
            }
        )

    df = pd.DataFrame(rows).sort_values(by="iteration")
    csv_path = out_dir / "iteration_convergence.csv"
    df.to_csv(csv_path, index=False)

    plot_png = out_dir / "iteration_convergence.png"
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(df["iteration"], df["n_selected_pairs"], marker="o", label="n_selected_pairs")
    ax.plot(df["iteration"], df["n_oversplit_suggestions"], marker="o", label="n_oversplit_suggestions")
    ax.set_title("UnitMatch Iteration Convergence")
    ax.set_xlabel("iteration")
    ax.set_ylabel("count")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(plot_png, dpi=180)
    plt.close(fig)

    return {
        "ok": True,
        "table_csv": str(csv_path),
        "plot_png": str(plot_png),
    }


def generate_unitmatch_static_report_pack(
    *,
    output_dir: str | Path,
    logger: Any = None,
    config: UnitMatchReportConfig = UnitMatchReportConfig(),
) -> dict[str, Any]:
    log = _ensure_logger(logger)
    base = Path(output_dir)

    throughput_root = base / str(config.throughput_subdir_name)
    output_root = base / str(config.output_subdir_name)
    report_root = base / str(config.report_subdir_name)
    report_root.mkdir(parents=True, exist_ok=True)

    summary_path = output_root / "unitmatch_summary.json"
    summary: dict[str, Any] = _read_json(summary_path) if summary_path.exists() else {}

    artifacts = summary.get("artifacts", {}) if isinstance(summary, dict) else {}
    iteration_dir = _latest_iteration_dir(throughput_root)

    matrix_path = None
    partner_summary_path = None
    selected_pairs_path = None
    suggestions_path = None
    candidates_path = None

    if isinstance(artifacts, dict):
        matrix_path_s = artifacts.get("score_matrix_npy")
        if matrix_path_s:
            matrix_path = Path(str(matrix_path_s))

        partner_path_s = artifacts.get("per_unit_partner_summary_json")
        if partner_path_s:
            partner_summary_path = Path(str(partner_path_s))

        selected_path_s = artifacts.get("selected_non_conflicting_pairs_csv")
        if selected_path_s:
            selected_pairs_path = Path(str(selected_path_s))

        suggestions_path_s = artifacts.get("oversplit_suggestions_csv")
        if suggestions_path_s:
            suggestions_path = Path(str(suggestions_path_s))

        candidates_path_s = artifacts.get("candidates_csv")
        if candidates_path_s:
            candidates_path = Path(str(candidates_path_s))

    if iteration_dir is not None:
        if matrix_path is None:
            matrix_path = iteration_dir / "deepunitmatch_similarity_matrix.npy"
        if partner_summary_path is None:
            partner_summary_path = iteration_dir / "per_unit_partner_summary.json"
        if selected_pairs_path is None:
            selected_pairs_path = iteration_dir / "selected_non_conflicting_pairs.csv"

    if suggestions_path is None:
        suggestions_path = output_root / "oversplit_suggestions.csv"
    if candidates_path is None:
        candidates_path = output_root / "match_candidates.csv"

    report: dict[str, Any] = {
        "timestamp": str(datetime.now()),
        "output_dir": str(base),
        "report_root": str(report_root),
        "inputs": {
            "summary_json": str(summary_path),
            "iteration_dir": (str(iteration_dir) if iteration_dir is not None else None),
            "matrix_npy": (str(matrix_path) if matrix_path is not None else None),
            "partner_summary_json": (str(partner_summary_path) if partner_summary_path is not None else None),
            "candidates_csv": (str(candidates_path) if candidates_path is not None else None),
            "suggestions_csv": (str(suggestions_path) if suggestions_path is not None else None),
            "selected_pairs_csv": (str(selected_pairs_path) if selected_pairs_path is not None else None),
        },
        "outputs": {},
    }

    scores: list[float] = []
    if candidates_path is not None and candidates_path.exists():
        try:
            df = pd.read_csv(candidates_path)
            scores = [_safe_float(v) for v in list(df.get("score", []))]
        except Exception as exc:
            log.warning("UnitMatch report: failed reading candidates CSV %s (%s)", candidates_path, exc)
    report["outputs"]["score_plots"] = _render_score_plots(scores, report_root)

    if matrix_path is None:
        report["outputs"]["matrix_heatmap"] = {"ok": False, "reason": "matrix_path_missing"}
    else:
        report["outputs"]["matrix_heatmap"] = _render_matrix_heatmap(
            matrix_path=matrix_path,
            out_dir=report_root,
            max_units=int(config.max_heatmap_units),
        )

    if partner_summary_path is None:
        report["outputs"]["partner_reports"] = {"ok": False, "reason": "partner_summary_path_missing"}
    else:
        report["outputs"]["partner_reports"] = _render_partner_reports(partner_summary_path, report_root)

    report["outputs"]["merge_diagnostics"] = _render_merge_diagnostics(
        suggestions_csv=suggestions_path,
        selected_csv=selected_pairs_path,
        out_dir=report_root,
    )

    merge_history_path = throughput_root / "merge_history.json"
    convergence_payload: dict[str, Any] = {}
    if merge_history_path.exists():
        try:
            convergence_payload = _read_json(merge_history_path)
        except Exception:
            convergence_payload = {}
    elif isinstance(summary, dict):
        convergence_payload = summary.get("iteration_convergence_report", {})

    report["outputs"]["convergence"] = _render_convergence(convergence_payload, report_root)

    manifest_path = report_root / "report_manifest.json"
    with manifest_path.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    log.info("UnitMatch static report pack written to %s", report_root)
    return report
