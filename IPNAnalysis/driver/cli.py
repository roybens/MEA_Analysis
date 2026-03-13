# ==========================================================
# driver/cli.py
# CLI entry point for run_pipeline_driver (batch well processing).
# ==========================================================

import argparse
import os
import sys
import time
from pathlib import Path

import pandas as pd

from .file_scanner import find_data_files, read_h5_recording_map
from .launcher import launch_sorting_subprocess, setup_driver_logger

try:
    from ..config_loader import load_config, resolve_args, build_extra_args
except ImportError:
    try:
        from MEA_Analysis.IPNAnalysis.config_loader import load_config, resolve_args, build_extra_args
    except ImportError:
        from config_loader import load_config, resolve_args, build_extra_args


# ------------------------------------------------------------------
# Argument parser
# ------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="MEA Pipeline Driver — batch processes MEA data files",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # --- Positional ---
    parser.add_argument("path", type=str,
        help="File or directory path containing MEA data")

    # --- Input / Output ---
    io_group = parser.add_argument_group("input/output")
    io_group.add_argument("--config", type=str, default=None,
        help="Path to config JSON file (CLI flags always override config)")
    io_group.add_argument("--output-dir", type=str, default=None,
        help="Output directory for all results")
    io_group.add_argument("--checkpoint-dir", type=str, default=None,
        help="Checkpoint directory (default: <output-dir>/checkpoints)")
    io_group.add_argument("--output-subdir-after-well", type=str, default=None,
        help="Optional single subdirectory appended under each resolved well output directory")
    io_group.add_argument("--export-to-phy", action="store_true",
        help="Export results to Phy format")
    io_group.add_argument("--clean-up", action="store_true",
        help="Remove intermediate files after processing")

    # --- Filtering (directory mode only) ---
    filter_group = parser.add_argument_group("filtering (directory mode only)")
    filter_group.add_argument("--reference", type=str, default=None,
        help="Excel file to filter runs by assay type\n"
             "(must have 'Run #' and 'Assay' columns)")
    filter_group.add_argument(
        "--type", nargs="+",
        default=['network today', 'network today/best'],
        help="Assay types to include\n"
             "(default: 'network today' 'network today/best')",
    )

    # --- Sorting (passed to each well) ---
    sort_group = parser.add_argument_group("sorting (passed to each well)")
    sort_group.add_argument("--sorter", type=str, default=None,
        help="Spike sorter to use (default: kilosort4)")
    sort_group.add_argument("--docker", type=str, default=None,
        help="Docker image name for containerized sorting")
    sort_group.add_argument("--skip-spikesorting", action="store_true",
        help="Run spike detection only, skip full sorting")
    sort_group.add_argument("--unitmatch-merge-units", action="store_true",
        help="Run UnitMatch dry-run/merge phase (optional, passed to each well)")
    sort_group.add_argument("--unitmatch-dry-run", action="store_true",
        help="When UnitMatch is enabled, run in dry-run mode (passed to each well)")
    sort_group.add_argument("--unitmatch-scored-dry-run",
        action=argparse.BooleanOptionalAction, default=None,
        help="When UnitMatch dry-run is enabled, attempt backend scoring (default: enabled)")
    sort_group.add_argument("--unitmatch-output-subdir-name", type=str, default=None,
        help="UnitMatch artifact subdirectory under each well output (default: unitmatch_outputs)")
    sort_group.add_argument("--unitmatch-throughput-subdir-name", type=str, default=None,
        help="UnitMatch throughput subdirectory under each well output (default: unitmatch_throughput)")
    sort_group.add_argument("--unitmatch-max-candidate-pairs", type=int, default=None,
        help="Max candidate pairs for UnitMatch scoring (-1 unlimited, 0 none, default: 20000)")
    sort_group.add_argument("--unitmatch-oversplit-min-probability", type=float, default=None,
        help="Minimum UnitMatch probability for oversplit suggestions (default: 0.80)")
    sort_group.add_argument("--unitmatch-oversplit-max-suggestions", type=int, default=None,
        help="Max oversplit suggestions to emit (-1 unlimited, 0 none, default: 2000)")
    sort_group.add_argument("--unitmatch-apply-merges", action="store_true",
        help="Apply UnitMatch-selected merges instead of report-only mode")
    sort_group.add_argument("--unitmatch-recursive", action="store_true",
        help="Recursively rerun UnitMatch merge iterations until convergence or iteration cap")
    sort_group.add_argument("--unitmatch-max-iterations", type=int, default=None,
        help="Maximum recursive UnitMatch iterations (-1 uncapped, default: 5)")
    sort_group.add_argument("--unitmatch-max-spikes-per-unit", type=int, default=None,
        help="Max spikes per unit for UnitMatch raw-waveform generation (-1 uncapped, default: 100)")
    sort_group.add_argument("--unitmatch-keep-all-iterations",
        action=argparse.BooleanOptionalAction, default=None,
        help="Keep all per-iteration throughput artifacts (default: enabled)")
    sort_group.add_argument("--unitmatch-generate-reports",
        action=argparse.BooleanOptionalAction, default=None,
        help="Generate static UnitMatch report pack from existing artifacts (default: enabled)")
    sort_group.add_argument("--unitmatch-report-subdir-name", type=str, default=None,
        help="UnitMatch report output subdirectory under each well output (default: unitmatch_reports)")
    sort_group.add_argument("--unitmatch-report-max-heatmap-units", type=int, default=None,
        help="Maximum units rendered in UnitMatch similarity heatmap (default: 200)")

    # --- Plotting (passed to each well) ---
    plot_group = parser.add_argument_group("plotting (passed to each well)")
    plot_group.add_argument("--plot-mode", choices=["separate", "merged"], default=None,
        help="Plot raster and network on separate axes or merged twin-axis\n(default: separate)")
    plot_group.add_argument(
        "--raster-sort",
        choices=["none", "firing_rate", "location_y", "unit_id"],
        default=None,
        help="How to sort units on raster y-axis (default: none)",
    )
    plot_group.add_argument("--plot-debug", action="store_true",
        help="Overlay burst and superburst intervals on raster plot")
    plot_group.add_argument("--fixed-y", action="store_true",
        help="Use fixed y-axis limits for raster plots — run once without it first to generate summary")

    # --- Curation (passed to each well) ---
    cur_group = parser.add_argument_group("curation (passed to each well)")
    cur_group.add_argument("--no-curation", action="store_true",
        help="Skip automatic unit curation")
    cur_group.add_argument("--params", type=str, default=None,
        help="JSON string or file path with quality thresholds")

    # --- Run Control ---
    ctrl_group = parser.add_argument_group("run control")
    ctrl_group.add_argument("--force-restart", action="store_true",
        help="Ignore checkpoints and restart from scratch")
    ctrl_group.add_argument("--reanalyze-bursts", action="store_true",
        help="Re-run burst analysis on existing spike times")
    ctrl_group.add_argument(
        "--resume-from", "--resume_from", dest="resume_from",
        type=str, default=None,
        choices=["preprocessing", "sorting", "merge", "analyzer", "reports"],
        help="Rewind checkpoint and resume each well from this stage",
    )
    ctrl_group.add_argument("--dry", action="store_true",
        help="Print what would run without any processing")
    ctrl_group.add_argument("--debug", action="store_true",
        help="Enable verbose logging")

    return parser


# ------------------------------------------------------------------
# Directory/file processing helpers
# ------------------------------------------------------------------

def _process_h5_file(
    file_path: str,
    args: argparse.Namespace,
    extra_arg_string: str,
    valid_runs: set | None,
    logger,
) -> None:
    """Process a single .h5 file: read recording map and launch per-well subprocesses."""
    try:
        run_id = int(Path(file_path).parent.name)
        if valid_runs is not None and run_id not in valid_runs:
            logger.info(f"[SKIP] Run {run_id} not in reference list")
            return
    except Exception as e:
        logger.warning(f"⚠️ Could not extract run ID from {file_path}: {e}")

    try:
        recording_map = read_h5_recording_map(file_path)
    except Exception as e:
        logger.error(f"Error opening HDF5 file {file_path}: {e}")
        return

    for recording, wells in recording_map.items():
        for well in wells:
            if args.dry:
                logger.info(f"[DRY-RUN] Would process {file_path} / {well}")
            else:
                logger.info(
                    f"Processing : {file_path} recording : {recording} well_id : {well}"
                )
                launch_sorting_subprocess(str(file_path), recording, well, extra_arg_string)


def _load_valid_runs(resolved: dict, logger) -> set | None:
    """Load the set of valid run IDs from a reference Excel file if provided."""
    if not resolved.get('reference_file'):
        return None
    try:
        df = pd.read_excel(resolved['reference_file'])
        filtered = df[
            df['Assay'].str.lower().isin(
                [t.lower() for t in resolved['assay_types']]
            )
        ]
        valid_runs = set(filtered['Run #'].astype(int).tolist())
        logger.info(
            f"Reference filter applied: {len(valid_runs)} valid runs from "
            f"'{resolved['reference_file']}' for assay types {resolved['assay_types']}"
        )
        return valid_runs
    except Exception as e:
        logger.error(f"Failed to load reference: {e}")
        sys.exit(1)


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

def main():
    parser = _build_parser()
    args = parser.parse_args()
    config = load_config(args.config)
    resolved = resolve_args(args, config)

    # Normalise output / checkpoint directories before logging setup.
    if args.output_dir:
        args.output_dir = os.path.normpath(os.path.abspath(args.output_dir))
        os.makedirs(args.output_dir, exist_ok=True)
        if args.checkpoint_dir:
            args.checkpoint_dir = os.path.normpath(os.path.abspath(args.checkpoint_dir))
        else:
            args.checkpoint_dir = os.path.join(args.output_dir, "checkpoints")
    else:
        if args.checkpoint_dir:
            args.checkpoint_dir = os.path.normpath(os.path.abspath(args.checkpoint_dir))
        else:
            args.checkpoint_dir = None

    logger = setup_driver_logger(log_path=args.output_dir)
    logger.info(f"Starting pipeline driver on path: {args.path}")
    logger.debug(f"Parsed arguments:\n{args}")

    path = Path(args.path)
    if not path.exists():
        logger.error(f"❌ Path does not exist: {path}")
        sys.exit(1)

    extra_arg_string = build_extra_args(resolved, args)
    logger.debug(f"Constructed extra argument string for subprocesses: {extra_arg_string}")

    logger.info(f"start time : {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    start_clock = time.time()

    # ------------------------------------------------------------------
    # Directory Mode
    # ------------------------------------------------------------------
    if path.is_dir():
        logger.info(f"[DIR MODE] Scanning directory: {path}")
        valid_runs = _load_valid_runs(resolved, logger)

        all_files = find_data_files(path, logger=logger)
        if not all_files:
            logger.critical("❌ No data files found (.h5, .nwb, .raw)")
            sys.exit(1)

        for fp in all_files:
            suffix = fp.suffix
            if suffix == ".h5":
                _process_h5_file(str(fp), args, extra_arg_string, valid_runs, logger)
            elif suffix == ".nwb":
                logger.info(f"[PLACEHOLDER] NWB file support not implemented yet: {fp}")
            elif suffix == ".raw":
                logger.info(f"[PLACEHOLDER] Raw binary file support not implemented yet: {fp}")
            else:
                logger.warning(f"[SKIP] Unknown file type: {fp}")

    # ------------------------------------------------------------------
    # Single File Mode
    # ------------------------------------------------------------------
    elif path.is_file():
        if path.suffix == ".h5":
            _process_h5_file(str(path), args, extra_arg_string, None, logger)
        elif path.suffix == ".nwb":
            logger.info(f"[PLACEHOLDER] NWB file support not implemented yet: {path}")
            sys.exit(1)
        elif path.suffix == ".raw":
            logger.info(f"[PLACEHOLDER] Raw file support not implemented yet: {path}")
            sys.exit(1)
        else:
            logger.error(f"Unsupported file extension: {path.suffix}")
            sys.exit(1)

    else:
        logger.error(f"Invalid path type: {path}")
        sys.exit(1)

    end_clock = time.time()
    logger.info(f"Pipeline driver completed in {end_clock - start_clock:.2f} seconds")
