# ==========================================================
# run_pipeline_driver.py
# Author: Mandar Patil
# LLM Assisted Edits: Yes ChatGPT-4
#
# REFACTORED: The implementation now lives in the `driver/` sub-package.
# This file is kept as a thin backward-compatible facade so that existing
# callers (CLI invocations, imports) continue to work without any changes.
# ==========================================================

import os
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"

import sys
from pathlib import Path

# Allow direct-script execution (`python run_pipeline_driver.py ...`) by
# adding the IPNAnalysis folder to sys.path before importing the driver package.
_script_dir = Path(__file__).resolve().parent
if str(_script_dir) not in sys.path:
    sys.path.insert(0, str(_script_dir))
_root_dir = _script_dir.parent.parent
if str(_root_dir) not in sys.path:
    sys.path.insert(0, str(_root_dir))


# Re-export the public symbols for existing import callers.
try:
    if __package__:
        from .driver import main, launch_sorting_subprocess, setup_driver_logger
        from .driver.file_scanner import find_data_files, read_h5_recording_map
    else:
        from driver import main, launch_sorting_subprocess, setup_driver_logger
        from driver.file_scanner import find_data_files, read_h5_recording_map
except ImportError:
    from MEA_Analysis.IPNAnalysis.driver import (
        main,
        launch_sorting_subprocess,
        setup_driver_logger,
    )
    from MEA_Analysis.IPNAnalysis.driver.file_scanner import find_data_files, read_h5_recording_map

if __name__ == "__main__":
    main()
