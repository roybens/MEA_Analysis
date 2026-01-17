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
def launch_sorting_subprocess(file_path, stream_id, extra_args=""):
    command = f"python3 {BASE_FILE_PATH}/mea_analysis_routine.py '{file_path}' --well {stream_id} {extra_args}"
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
    parser = argparse.ArgumentParser(description="MEA Pipeline Driver")

    parser.add_argument("path", type=str, help="File or directory path containing MEA data")
    parser.add_argument("--reference", type=str, help="Optional Excel file to filter runs")
    parser.add_argument("--type", nargs="+", default=['network today', 'network today/best'], help="Assay types to include")
    parser.add_argument("--params", type=str, help="JSON file or string with quality thresholds")
    parser.add_argument("--docker", type=str, help="Docker image name if using containerized sorting")
    #parser.add_argument("--resume", action="store_true", help="Resume from checkpoint"). #auto resumes
    parser.add_argument("--skip-spikesorting", action="store_true", help="Skip spike sorting stage")
    parser.add_argument("--force-restart", action="store_true", help="Restart even if checkpoint complete")
    parser.add_argument("--sorter", default="kilosort4", help="Sorter to use")
    parser.add_argument("--debug", action="store_true", help="Debug mode")
    parser.add_argument("--dry", action="store_true", help="Dry run only (no processing)")
    parser.add_argument("--clean-up", action="store_true", help="Clear sorter outputs etc.")
    parser.add_argument("--export-to-phy", action="store_true", help="Export results to Phy format")
    parser.add_argument("--checkpoint-dir", type=str, default=None, help="Checkpoint directory")
    parser.add_argument("--no-curation", action="store_true", help="Skip automatic curation")
    parser.add_argument("--output-dir", type=str, help="Output directory for results")

    # -------------------------------------------------------------------
    # HPC/Slurm-friendly knobs (opt-in)
    # -------------------------------------------------------------------
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=None,
        help="Number of parallel workers used inside each well subprocess (defaults to 16, or SLURM_CPUS_PER_TASK if set).",
    )
    parser.add_argument(
        "--chunk-duration",
        type=str,
        default=None,
        help="Chunk duration for SpikeInterface jobs (e.g. '1s', '500ms').",
    )

    parser.add_argument(
        "--scratch-dir",
        type=str,
        default=None,
        help=(
            "Optional fast local scratch base directory (e.g. $SLURM_TMPDIR). "
            "If set, heavy intermediate folders (binary/sorter/analyzer) are written there per well."
        ),
    )
    parser.add_argument(
        "--stage-back",
        choices=["none", "sorter", "all"],
        default="sorter",
        help=(
            "When using --scratch-dir, controls what is copied back into --output-dir at the end. "
            "'sorter' stages back sorter_output only (default). 'all' stages binary, sorter_output, analyzer_output."
        ),
    )
    parser.add_argument(
        "--stage-back-mode",
        choices=["copy", "move"],
        default="copy",
        help=(
            "When using --scratch-dir and stage-back is enabled, choose whether to copy or move results back. "
            "'move' reclaims scratch space but removes the staged folders from scratch."
        ),
    )

    # -------------------------------------------------------------------
    # GPU-node carveouts (opt-in)
    # -------------------------------------------------------------------
    parser.add_argument(
        "--cuda-visible-devices",
        type=str,
        default=None,
        help="Value to set CUDA_VISIBLE_DEVICES for each well subprocess (e.g. '0' or '0,1').",
    )
    parser.add_argument(
        "--require-gpu",
        action="store_true",
        help="Fail fast if --docker is not used and no CUDA GPU is visible.",
    )

    args = parser.parse_args()

    logger = None

    # -------------------------------------------------------------------
    # Normalize output + checkpoint directory before logging setup
    # -------------------------------------------------------------------

    # If output-dir is missing, pick a sensible default.
    # This avoids failing later because mea_analysis_routine.py requires --output-dir.
    if not args.output_dir:
        env_out = os.environ.get("MEA_ANALYSIS_OUTPUT_DIR")
        slurm_submit = os.environ.get("SLURM_SUBMIT_DIR")
        slurm_job = os.environ.get("SLURM_JOB_ID")
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        if env_out:
            args.output_dir = env_out
        elif slurm_submit:
            job_tag = slurm_job if slurm_job else f"manual_{timestamp}"
            args.output_dir = os.path.join(slurm_submit, "mea_analysis_outputs", job_tag)
        else:
            args.output_dir = os.path.join(os.getcwd(), "mea_analysis_outputs", f"manual_{timestamp}")

    # Normalize output + checkpoint dir
    args.output_dir = os.path.normpath(os.path.abspath(args.output_dir))
    os.makedirs(args.output_dir, exist_ok=True)

    if args.checkpoint_dir:
        args.checkpoint_dir = os.path.normpath(os.path.abspath(args.checkpoint_dir))
    else:
        args.checkpoint_dir = os.path.join(args.output_dir, "checkpoints")

    # -------------------------------------------------------------------
    # Logger setup (driver)
    # -------------------------------------------------------------------

    logger = setup_driver_logger(log_path=args.output_dir)

    logger.info(f"Starting pipeline driver on path: {args.path}")
    logger.debug(f"Parsed arguments:\n{args}")


    path = Path(args.path)
    if not path.exists():
        logger.error(f"❌ Path does not exist: {path}")
        sys.exit(1)

    # ------------------------------------------------------
    # Build extra argument string for subprocess
    # ------------------------------------------------------
    extra_args = []
    #if args.resume: extra_args.append("--resume")
    if args.force_restart: extra_args.append("--force-restart")
    if args.docker: extra_args.append(f"--docker {args.docker}")
    if args.sorter: extra_args.append(f"--sorter {args.sorter}")
    if args.debug: extra_args.append("--debug")
    if args.params: extra_args.append(f"--params '{args.params}'")
    if args.clean_up: extra_args.append("--clean-up")
    if args.checkpoint_dir: extra_args.append(f"--checkpoint-dir '{args.checkpoint_dir}'")
    if args.export_to_phy: extra_args.append("--export-to-phy")
    if args.no_curation: extra_args.append("--no-curation")
    if args.skip_spikesorting: extra_args.append("--skip-spikesorting")
    if args.output_dir: extra_args.append(f"--output-dir '{args.output_dir}'")

    # HPC knobs
    if args.n_jobs is not None: extra_args.append(f"--n-jobs {int(args.n_jobs)}")
    if args.chunk_duration is not None: extra_args.append(f"--chunk-duration '{args.chunk_duration}'")

    # GPU-node carveouts
    if args.cuda_visible_devices is not None:
        extra_args.append(f"--cuda-visible-devices '{args.cuda_visible_devices}'")
    if args.require_gpu:
        extra_args.append("--require-gpu")

    # Scratch carveouts
    if args.scratch_dir is not None:
        extra_args.append(f"--scratch-dir '{args.scratch_dir}'")
        if args.stage_back:
            extra_args.append(f"--stage-back {args.stage_back}")
        if args.stage_back_mode:
            extra_args.append(f"--stage-back-mode {args.stage_back_mode}")

    extra_arg_string = " ".join(extra_args)
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
                        wells = list(h5f["wells"].keys())
                except Exception as e:
                    logger.error(f"Failed to open HDF5: {e}")
                    continue

                for well in wells:
                    if args.dry:
                        logger.info(f"[DRY-RUN] Would process {file_path} / {well}")
                    else:
                        launch_sorting_subprocess(file_path, well, extra_arg_string)

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
                    wells = list(h5f["wells"].keys())
            except Exception as e:
                logger.error(f"Could not open file: {e}")
                sys.exit(1)

            for well in wells:
                if args.dry:
                    logger.info(f"[DRY-RUN] Would process {path} / {well}")
                else:
                    launch_sorting_subprocess(str(path), well, extra_arg_string)

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