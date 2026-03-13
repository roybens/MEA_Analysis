# ==========================================================
# pipeline/cli.py
# CLI entry point for mea_analysis_routine (single-well processing).
# ==========================================================

import argparse
import sys
import traceback

from .options import MEARunOptions
from .runner import run_mea_pipeline

try:
    from ..config_loader import load_config, resolve_args
except ImportError:
    try:
        from MEA_Analysis.IPNAnalysis.config_loader import load_config, resolve_args
    except ImportError:
        from config_loader import load_config, resolve_args


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="MEA Analysis Routine — processes a single well from an MEA recording",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # --- Positional ---
    parser.add_argument(
        "file_path",
        help="Path to .h5, .nwb, or .raw MEA recording file",
    )

    # --- Input / Output ---
    io_group = parser.add_argument_group("input/output")
    io_group.add_argument("--config", type=str, default=None,
        help="Path to config JSON file (CLI flags always override config)")
    io_group.add_argument("--well", required=True,
        help="Well ID to process (e.g. well000)")
    io_group.add_argument("--rec", type=str, default=None,
        help="Recording name inside HDF5 file (default: rec0000)")
    io_group.add_argument("--output-dir", type=str, default=None,
        help="Output directory for results")
    io_group.add_argument("--checkpoint-dir", type=str, default=None,
        help="Checkpoint directory (default: <output-dir>/checkpoints)")
    io_group.add_argument("--output-subdir-after-well", type=str, default=None,
        help="Optional single subdirectory appended under the resolved well output directory")
    io_group.add_argument("--export-to-phy", action="store_true",
        help="Export results to Phy format")
    io_group.add_argument("--clean-up", action="store_true",
        help="Remove intermediate files after processing")

    # --- Sorting ---
    sort_group = parser.add_argument_group("sorting")
    sort_group.add_argument("--sorter", type=str, default=None,
        help="Spike sorter to use (default: kilosort4)")
    sort_group.add_argument("--docker", type=str, default=None,
        help="Docker image name for containerized sorting")
    sort_group.add_argument("--skip-spikesorting", action="store_true",
        help="Run spike detection only, skip full sorting")

    # --- Plotting ---
    plot_group = parser.add_argument_group("plotting")
    plot_group.add_argument("--plot-mode", choices=["separate", "merged"], default=None,
        help="Plot raster and network on separate axes or merged twin-axis\n(default: separate)")
    plot_group.add_argument("--raster-sort",
        choices=["none", "firing_rate", "location_y", "unit_id"], default=None,
        help="How to sort units on raster y-axis (default: none)")
    plot_group.add_argument("--plot-debug", action="store_true",
        help="Overlay burst and superburst intervals on raster plot")
    plot_group.add_argument("--fixed-y", action="store_true",
        help="Use fixed y-axis limits for raster plots — run once without it first to generate summary")

    # --- Curation ---
    cur_group = parser.add_argument_group("curation")
    cur_group.add_argument("--no-curation", action="store_true",
        help="Skip automatic unit curation")
    cur_group.add_argument("--params", type=str, default=None,
        help="JSON string or file path with quality thresholds")

    # --- Run Control ---
    ctrl_group = parser.add_argument_group("run control")
    ctrl_group.add_argument("--force-restart", action="store_true",
        help="Ignore checkpoint and restart from scratch")
    ctrl_group.add_argument(
        "--resume-from", "--resume_from", dest="resume_from", type=str, default=None,
        choices=["preprocessing", "sorting", "merge", "analyzer", "reports"],
        help="Resume by rewinding checkpoint to just before this stage and rerunning from there",
    )
    ctrl_group.add_argument("--reanalyze-bursts", action="store_true",
        help="Re-run burst analysis on existing spike times only")
    ctrl_group.add_argument("--debug", action="store_true",
        help="Enable verbose logging")

    # --- UnitMatch ---
    parser.add_argument("--unitmatch-merge-units", action='store_true',
        help="Optional (default off): run UnitMatch integration as an alternative to auto_merge_units.")
    parser.add_argument("--unitmatch-dry-run", action='store_true',
        help="When UnitMatch is enabled, produce UnitMatch reports without applying merges.")
    parser.add_argument("--unitmatch-scored-dry-run",
        action=argparse.BooleanOptionalAction, default=None,
        help="When UnitMatch dry-run is enabled, attempt backend scoring (default: enabled).")
    parser.add_argument("--unitmatch-output-subdir-name", type=str, default=None,
        help="UnitMatch artifact subdirectory under output_dir (default: unitmatch_outputs).")
    parser.add_argument("--unitmatch-throughput-subdir-name", type=str, default=None,
        help="UnitMatch throughput subdirectory under output_dir (default: unitmatch_throughput).")
    parser.add_argument("--unitmatch-max-candidate-pairs", type=int, default=None,
        help="Maximum UnitMatch candidate pairs (-1 unlimited, 0 none, default: 20000).")
    parser.add_argument("--unitmatch-oversplit-min-probability", type=float, default=None,
        help="Minimum UnitMatch probability for oversplit suggestions (default: 0.80).")
    parser.add_argument("--unitmatch-oversplit-max-suggestions", type=int, default=None,
        help="Maximum oversplit suggestions (-1 unlimited, 0 none, default: 2000).")
    parser.add_argument("--unitmatch-apply-merges", action='store_true',
        help="Apply conflict-free top UnitMatch suggestions to mutate sorting.")
    parser.add_argument("--unitmatch-recursive", action='store_true',
        help="Recursively run UnitMatch merge iterations until convergence or cap.")
    parser.add_argument("--unitmatch-max-iterations", type=int, default=None,
        help="Maximum recursive UnitMatch iterations (-1 uncapped, default: 5).")
    parser.add_argument("--unitmatch-max-spikes-per-unit", type=int, default=None,
        help="Max spikes per unit for UnitMatch raw-waveform generation (-1 uncapped, default: 100).")
    parser.add_argument("--unitmatch-keep-all-iterations",
        action=argparse.BooleanOptionalAction, default=None,
        help="Keep all unitmatch_throughput iteration folders (default: enabled).")
    parser.add_argument("--unitmatch-generate-reports",
        action=argparse.BooleanOptionalAction, default=True,
        help="Generate static UnitMatch report pack from existing artifacts (default: enabled).")
    parser.add_argument("--unitmatch-report-subdir-name", type=str, default="unitmatch_reports",
        help="UnitMatch report output subdirectory under output_dir (default: unitmatch_reports).")
    parser.add_argument("--unitmatch-report-max-heatmap-units", type=int, default=200,
        help="Maximum units rendered in UnitMatch similarity heatmap (default: 200).")

    # --- Auto-merge ---
    parser.add_argument("--auto-merge-units", action='store_true',
        help="Optional (default off): run SpikeInterface auto_merge_units during analyzer stage.")
    parser.add_argument("--auto-merge-template-diff-thresh", default="0.05,0.15,0.25",
        help="Comma-separated template_diff_thresh values for auto-merge "
             "(used with preset x_contaminations).")
    parser.add_argument("--rerun-analyzer", action='store_true',
        help="Recompute analyzer_output even if checkpoint says complete "
             "(does not rerun spikesorting).")

    return parser


def main():
    parser = _build_parser()
    args = parser.parse_args()
    config = load_config(args.config)
    resolved = resolve_args(args, config)

    plot_mode = resolved["plot_mode"]
    plot_debug = resolved["plot_debug"]
    raster_sort = resolved["raster_sort"]
    fixed_y = resolved["fixed_y"]
    sorter = resolved["sorter"]
    thresholds = resolved["quality_thresholds"]
    rec = args.rec or "rec0000"

    try:
        if args.reanalyze_bursts:
            print("Re-analyzing bursts only on existing spike times...")

        result = run_mea_pipeline(
            MEARunOptions(
                file_path=args.file_path,
                stream_id=args.well,
                recording_num=rec,
                output_root=resolved["output_dir"],
                checkpoint_root=resolved["checkpoint_dir"],
                sorter=sorter,
                docker_image=resolved["docker_image"],
                verbose=bool(args.debug),
                cleanup=bool(resolved["clean_up"]),
                force_restart=bool(args.force_restart),
                resume_from=args.resume_from,
                um_kwargs={
                    "merge_units": bool(args.unitmatch_merge_units),
                    "dry_run": bool(args.unitmatch_dry_run),
                    "scored_dry_run": bool(resolved["unitmatch_scored_dry_run"]),
                    "output_subdir_name": str(resolved["unitmatch_output_subdir_name"]),
                    "throughput_subdir_name": str(resolved["unitmatch_throughput_subdir_name"]),
                    "max_candidate_pairs": int(resolved["unitmatch_max_candidate_pairs"]),
                    "oversplit_min_probability": float(
                        resolved["unitmatch_oversplit_min_probability"]
                    ),
                    "oversplit_max_suggestions": int(
                        resolved["unitmatch_oversplit_max_suggestions"]
                    ),
                    "apply_merges": bool(resolved["unitmatch_apply_merges"]),
                    "recursive": bool(resolved["unitmatch_recursive"]),
                    "max_iterations": int(resolved["unitmatch_max_iterations"]),
                    "max_spikes_per_unit": int(resolved["unitmatch_max_spikes_per_unit"]),
                    "keep_all_iterations": bool(resolved["unitmatch_keep_all_iterations"]),
                    "generate_reports": bool(args.unitmatch_generate_reports),
                    "report_subdir_name": str(args.unitmatch_report_subdir_name),
                    "report_max_heatmap_units": int(args.unitmatch_report_max_heatmap_units),
                },
                am_kwargs={
                    "enabled": bool(args.auto_merge_units),
                    "template_diff_thresh": str(args.auto_merge_template_diff_thresh),
                },
                option_kwargs={
                    "force_rerun_analyzer": bool(args.rerun_analyzer),
                    "output_subdir_after_well": resolved.get("output_subdir_after_well"),
                },
                reanalyze_bursts=bool(args.reanalyze_bursts),
                skip_spikesorting=bool(args.skip_spikesorting),
                run_analyzer=True,
                run_reports=True,
                thresholds=thresholds,
                no_curation=bool(resolved["no_curation"]),
                export_to_phy=bool(resolved["export_to_phy"]),
                plot_mode=plot_mode,
                plot_debug=bool(plot_debug),
                raster_sort=raster_sort,
                fixed_y=bool(fixed_y),
                auto_merge_template_diff_thresh=str(args.auto_merge_template_diff_thresh),
            )
        )

        if result.reanalyzed_bursts:
            print("Burst Re-analysis Complete.")
            sys.exit(0)

        if result.skipped:
            sys.exit(0)

        print(f"Processing Complete for {args.well}")

    except Exception as e:
        print(f"CRITICAL FAILURE in {args.well}: {e}")
        traceback.print_exc()
        sys.exit(1)
