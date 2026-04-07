# ==========================================================
# mea_analysis_routine.py
# Author: Mandar Patil
# Contributors: Yuxin Ren, Shruti Shah, Adam Weiner
# LLM Assisted Edits: Yes ChatGPT-4, Claude sonnet 4.6, ChatGPT-5.3-Codex
#
# REFACTORED: The implementation now lives in the `pipeline/` sub-package.
# This file is kept as a thin backward-compatible facade so that existing
# callers (import statements, CLI invocations, test suites) continue to work
# without any changes.
# ==========================================================

import os
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"
os.environ["QT_QPA_PLATFORM"] = "offscreen"  # For headless environments

import sys
from pathlib import Path

# Allow direct-script execution (`python mea_analysis_routine.py ...`) by
# adding the IPNAnalysis folder (and repo root) to sys.path before importing
# the pipeline sub-package.
if not __package__:
    _script_dir = Path(__file__).resolve().parent
    if str(_script_dir) not in sys.path:
        sys.path.insert(0, str(_script_dir))
    _root_dir = _script_dir.parent.parent
    if str(_root_dir) not in sys.path:
        sys.path.insert(0, str(_root_dir))


# Re-export the entire public API so existing `from mea_analysis_routine import X`
# statements keep working unchanged.
try:
    if __package__:
        from .pipeline import (
            ProcessingStage,
            CHECKPOINT_SCHEMA_VERSION,
            NpEncoder,
            _default_um_kwargs,
            _default_am_kwargs,
            _default_option_kwargs,
            MEAPipeline,
            MEARunOptions,
            MEARunResult,
            run_mea_pipeline,
        )
        from .pipeline.runner import _normalize_resume_from_stage, _apply_resume_from_stage
        from .pipeline.cli import main
    else:
        from pipeline import (
            ProcessingStage,
            CHECKPOINT_SCHEMA_VERSION,
            NpEncoder,
            _default_um_kwargs,
            _default_am_kwargs,
            _default_option_kwargs,
            MEAPipeline,
            MEARunOptions,
            MEARunResult,
            run_mea_pipeline,
        )
        from pipeline.runner import _normalize_resume_from_stage, _apply_resume_from_stage
        from pipeline.cli import main
except ImportError:
    from MEA_Analysis.IPNAnalysis.pipeline import (
        ProcessingStage,
        CHECKPOINT_SCHEMA_VERSION,
        NpEncoder,
        _default_um_kwargs,
        _default_am_kwargs,
        _default_option_kwargs,
        MEAPipeline,
        MEARunOptions,
        MEARunResult,
        run_mea_pipeline,
    )
    from MEA_Analysis.IPNAnalysis.pipeline.runner import (
        _normalize_resume_from_stage,
        _apply_resume_from_stage,
    )
    from MEA_Analysis.IPNAnalysis.pipeline.cli import main

if __name__ == "__main__":
    main()
