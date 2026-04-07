# ==========================================================
# driver/launcher.py
# Subprocess launcher and driver-level logger setup.
# ==========================================================

import logging
import os
import shlex
import subprocess
import sys
import traceback
from datetime import datetime
from pathlib import Path

# Absolute path to the IPNAnalysis package directory.
_IPNANALYSIS_DIR = str(Path(__file__).resolve().parent.parent)


# ------------------------------------------------------------------
# Driver logger
# ------------------------------------------------------------------

def setup_driver_logger(log_path: str | None = None) -> logging.Logger:
    """Create (or retrieve) a ``"driver"`` logger writing to a timestamped file.

    Parameters
    ----------
    log_path:
        Directory in which to create the log file.  Defaults to
        the IPNAnalysis package directory.

    Returns
    -------
    logging.Logger
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base = log_path if log_path else _IPNANALYSIS_DIR
    log_file = os.path.join(base, f"run_pipeline_driver_{timestamp}.log")
    open(log_file, 'w').close()  # Clear / create log file

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


# ------------------------------------------------------------------
# Per-well subprocess launcher
# ------------------------------------------------------------------

def launch_sorting_subprocess(
    file_path: str,
    rec_name: str,
    stream_id: str,
    extra_args: str = "",
) -> None:
    """Launch ``mea_analysis_routine.py`` for a single well in a subprocess.

    Parameters
    ----------
    file_path:
        Absolute path to the .h5 data file.
    rec_name:
        Recording name (e.g. ``rec0000``).
    stream_id:
        Well ID (e.g. ``well000``).
    extra_args:
        Additional CLI flags built by :func:`config_loader.build_extra_args`.
    """
    routine_path = os.path.join(_IPNANALYSIS_DIR, "mea_analysis_routine.py")
    command = (
        f"python3 {routine_path} '{file_path}' "
        f"--rec {rec_name} --well {stream_id} {extra_args}"
    )
    logger = logging.getLogger("driver")
    logger.info(f"[DRIVER] Launching: {command}")

    try:
        subprocess.run(shlex.split(command), check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Subprocess failed for {stream_id} in {file_path}")
        logger.error(f"Exit Code: {e.returncode}")
        logger.error(traceback.format_exc())
