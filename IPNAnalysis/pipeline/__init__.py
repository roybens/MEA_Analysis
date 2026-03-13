# ==========================================================
# pipeline/__init__.py
# Public API for the IPNAnalysis pipeline package.
# ==========================================================

# Set environment variables early so they take effect before heavy imports
# (PyTorch/CUDA, matplotlib headless backend).
import os as _os
_os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
_os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from .stages import ProcessingStage, CHECKPOINT_SCHEMA_VERSION, NpEncoder
from .defaults import _default_um_kwargs, _default_am_kwargs, _default_option_kwargs
from .core import MEAPipeline
from .options import MEARunOptions, MEARunResult
from .runner import run_mea_pipeline

__all__ = [
    "ProcessingStage",
    "CHECKPOINT_SCHEMA_VERSION",
    "NpEncoder",
    "_default_um_kwargs",
    "_default_am_kwargs",
    "_default_option_kwargs",
    "MEAPipeline",
    "MEARunOptions",
    "MEARunResult",
    "run_mea_pipeline",
]
