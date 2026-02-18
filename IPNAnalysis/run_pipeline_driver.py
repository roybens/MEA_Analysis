# ==========================================================
# run_pipeline_driver.py
# Author: Mandar Patil
# LLM Assisted Edits: Yes ChatGPT-4

# ==========================================================

import time
import os
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"
import sys
import argparse
import subprocess
import h5py
import json
from pathlib import Path
import shlex
import traceback
import pandas as pd
import logging
from datetime import datetime

# This ensures Python looks in the same folder as this script for modules
script_dir = Path(__file__).resolve().parent
if str(script_dir) not in sys.path:
    sys.path.append(str(script_dir))

# Add the parent directory of MEA_Analysis to path so absolute imports work
root_dir = script_dir.parent.parent 
if str(root_dir) not in sys.path:
    sys.path.append(str(root_dir))

# Import your custom modules
try:
    from config_loader import load_config, resolve_args, build_extra_args
except ImportError as e:
    print(f"CRITICAL ERROR: Could not import helper modules. {e}")
    print(f"Current sys.path: {sys.path}")
    sys.exit(1)
BASE_FILE_PATH = str(Path(__file__).resolve().parent)

# ----------------------------------------------------------
# Logger setup
# ----------------------------------------------------------
def setup_driver_logger(log_path=None):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    if log_path is None:
        log_file= os.path.join(BASE_FILE_PATH, f"run_pipeline_driver_{timestamp}.log")
    else:
        log_file = os.path.join(log_path, f"run_pipeline_driver_{timestamp}.log")
    open(log_file, 'w').close()  # Clear log file on each run

    logger = logging.getLogger("driver")
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s')

    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    logger.info(f"Logging to {log_file}")
    return logger


# ----------------------------------------------------------
# Subprocess launcher
# ----------------------------------------------------------
def launch_sorting_subprocess(file_path, rec_name, stream_id, extra_args=""):
    
    # Calculate rec_name for 24-well MaxTwo plates
    well_id = int(stream_id[-3:])
    #rec_name = 'rec000' + str(well_id // 6) # each rec contain 6 wells (1 row on 24-well plate)

    command = f"python3 {BASE_FILE_PATH}/mea_analysis_routine.py '{file_path}' --rec {rec_name} --well {stream_id} {extra_args}"
    logger = logging.getLogger("driver")
    logger.info(f"[DRIVER] Launching: {command}")
    
    try:
        subprocess.run(shlex.split(command), check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Subprocess failed for {stream_id} in {file_path}")
        logger.error(f"Exit Code: {e.returncode}")
        logger.error(traceback.format_exc())


# ----------------------------------------------------------
# Main
# ----------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="MEA Pipeline Driver — batch processes MEA data files",
        formatter_class=argparse.RawTextHelpFormatter
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
    io_group.add_argument("--export-to-phy", action="store_true",
        help="Export results to Phy format")
    io_group.add_argument("--clean-up", action="store_true",
        help="Remove intermediate files after processing")

    # --- Filtering (directory mode only) ---
    filter_group = parser.add_argument_group("filtering (directory mode only)")
    filter_group.add_argument("--reference", type=str, default=None,
        help="Excel file to filter runs by assay type\n(must have 'Run #' and 'Assay' columns)")
    filter_group.add_argument("--type", nargs="+", default=['network today', 'network today/best'],
        help="Assay types to include\n(default: 'network today' 'network today/best')")

    # --- Sorting (passed to each well) ---
    sort_group = parser.add_argument_group("sorting (passed to each well)")
    sort_group.add_argument("--sorter", type=str, default=None,
        help="Spike sorter to use (default: kilosort4)")
    sort_group.add_argument("--docker", type=str, default=None,
        help="Docker image name for containerized sorting")
    sort_group.add_argument("--skip-spikesorting", action="store_true",
        help="Run spike detection only, skip full sorting")

    # --- Plotting (passed to each well) ---
    plot_group = parser.add_argument_group("plotting (passed to each well)")
    plot_group.add_argument("--plot-mode", choices=["separate", "merged"], default=None,
        help="Plot raster and network on separate axes or merged twin-axis\n(default: separate)")
    plot_group.add_argument("--raster-sort", choices=["none", "firing_rate", "location_y", "unit_id"], default=None,
        help="How to sort units on raster y-axis (default: none)")
    plot_group.add_argument("--plot-debug", action="store_true",
    help="Overlay burst and superburst intervals on raster plot")
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
    ctrl_group.add_argument("--dry", action="store_true",
        help="Print what would run without any processing")
    ctrl_group.add_argument("--debug", action="store_true",
        help="Enable verbose logging")

    args = parser.parse_args()
    config   = load_config(args.config)
    resolved = resolve_args(args, config)
    logger = None

    # -------------------------------------------------------------------
    # Normalize output + checkpoint directory before logging setup
    # -------------------------------------------------------------------

    if args.output_dir:
        # Make absolute normalized path
        args.output_dir = os.path.normpath(os.path.abspath(args.output_dir))
        os.makedirs(args.output_dir, exist_ok=True)

        # Default checkpoint is inside output
        if args.checkpoint_dir:
            # Respect user-provided checkpoint dir (normalize)
            args.checkpoint_dir = os.path.normpath(os.path.abspath(args.checkpoint_dir))
        else:
            # Set checkpoint inside output directory
            args.checkpoint_dir = os.path.join(args.output_dir, "checkpoints")

    else:
        # No output dir ⇒ keep legacy behavior
        if args.checkpoint_dir:
            args.checkpoint_dir = os.path.normpath(os.path.abspath(args.checkpoint_dir))
        else:
            args.checkpoint_dir = None  # handled later

    # -------------------------------------------------------------------
    # Logger setup (driver)
    # -------------------------------------------------------------------

    if args.output_dir:
        logger = setup_driver_logger(log_path=args.output_dir)
    else:
        logger = setup_driver_logger()

    logger.info(f"Starting pipeline driver on path: {args.path}")
    logger.debug(f"Parsed arguments:\n{args}")


    path = Path(args.path)
    if not path.exists():
        logger.error(f"❌ Path does not exist: {path}")
        sys.exit(1)

    # ------------------------------------------------------
    # Build extra argument string for subprocess
    # ------------------------------------------------------
    extra_arg_string = build_extra_args(resolved, args)
    logger.debug(f"Constructed extra argument string for subprocesses: {extra_arg_string}")

    logger.info(f"start time : {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    #ticker
    start_clock = time.time()
    # ------------------------------------------------------
    # Handle Directory Mode
    # ------------------------------------------------------
    if path.is_dir():
        logger.info(f"[DIR MODE] Scanning directory: {path}")

        valid_runs = None
        if args.reference:
            try:
                df = pd.read_excel(args.reference)
                filtered = df[df['Assay'].str.lower().isin([t.lower() for t in args.type])]
                valid_runs = set(filtered['Run #'].astype(int).tolist())
                logger.info(f"Reference filter applied: {len(valid_runs)} valid runs")
            except Exception as e:
                logger.error(f"Failed to load reference: {e}")
                sys.exit(1)
        h5_files,nwb_files,raw_files = [],[],[]
        try:
            
            file_name_pattern = "data.raw.h5"
            subfolder_name = "Network"
            import helper_functions as helper
            h5_files = helper.find_files_with_subfolder(path, file_name_pattern, subfolder_name)
        except Exception as e:
            logger.warning(f"Error finding HDF5 files: {e}")
        

        nwb_files = list(path.rglob("*.nwb"))
        raw_files = list(path.rglob("*.raw"))
        all_files = h5_files + nwb_files + raw_files

        if not all_files:
            logger.critical("❌ No data files found (.h5, .nwb, .raw)")
            sys.exit(1)

        for file_path in all_files:
            file_path = str(file_path)
            suffix = Path(file_path).suffix

            if suffix == ".h5":
                try:
                    run_id = int(Path(file_path).parent.name)
                    if valid_runs and run_id not in valid_runs:
                        logger.info(f"[SKIP] Run {run_id} not in reference list")
                        continue
                except Exception as e:
                    logger.warning(f"⚠️ Could not extract run ID from {file_path}: {e}")
                try:
                    with h5py.File(file_path, "r") as h5f:
                        recording_map = {
                            recording: list(h5f["recordings"][recording].keys())
                            for recording in h5f["recordings"].keys()
                        }
                except Exception as e:
                    logger.error(f"Error opening HDF5 file {file_path}: {e}")
                    sys.exit(1)
                for recording, wells in recording_map.items():
                    for well in wells:
                        if args.dry:
                            logger.info(f"[DRY-RUN] Would process {file_path} / {well}")
                        else:
                            logger.info(
                                f"Processing : {file_path} recording : {recording} well_id : {well}"
                            )
                            launch_sorting_subprocess(
                                str(file_path), recording, well, extra_arg_string
                            )
                            
            elif suffix == ".nwb":
                logger.info(f"[PLACEHOLDER] NWB file support not implemented yet: {file_path}")
                continue

            elif suffix == ".raw":
                logger.info(f"[PLACEHOLDER] Raw binary file support not implemented yet: {file_path}")
                continue

            else:
                logger.warning(f"[SKIP] Unknown file type: {file_path}")
                continue

    # ------------------------------------------------------
    # Handle Single File Mode
    # ------------------------------------------------------
    elif path.is_file():
        if path.suffix == ".h5":
            try:
                with h5py.File(path, "r") as h5f:
                    recording_map = {
                        recording: list(h5f["recordings"][recording].keys())
                        for recording in h5f["recordings"].keys()
                    }
            except Exception as e:
                logger.error(f"Error opening HDF5 file {path}: {e}")
                sys.exit(1)
            for recording, wells in recording_map.items():
                for well in wells:
                    if args.dry:
                        logger.info(f"[DRY-RUN] Would process {path} / {well}")
                    else:
                        logger.info(
                            f"Processing : {path} recording : {recording} well_id : {well}"
                        )
                        launch_sorting_subprocess(
                            str(path), recording, well, extra_arg_string
                        )
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

    # Stop the timer and log the duration
    end_clock = time.time()
    duration = end_clock - start_clock
    logger.info(f"Pipeline driver completed in {duration:.2f} seconds")
    
# ----------------------------------------------------------
# Entry
# ----------------------------------------------------------
if __name__ == "__main__":
    main()