"""Backend adapters for UnitMatch integrations."""

from .deepunitmatch_adapter import (
    BackendCheckResult,
    BackendInferenceResult,
    check_backend_availability,
    infer_similarity_matrix,
    preprocess_waveforms,
)

__all__ = [
    "BackendCheckResult",
    "BackendInferenceResult",
    "check_backend_availability",
    "infer_similarity_matrix",
    "preprocess_waveforms",
]
