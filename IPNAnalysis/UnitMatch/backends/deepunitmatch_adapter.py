from __future__ import annotations

import importlib
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class BackendCheckResult:
    available: bool
    reason_code: str
    import_target: str
    diagnostics: dict[str, Any]


@dataclass(frozen=True)
class BackendInferenceResult:
    success: bool
    reason_code: str
    import_ok: bool
    sim_matrix: Any | None
    diagnostics: dict[str, Any]


def _import_installed(module_name: str) -> Any:
    return importlib.import_module(module_name)


def check_backend_availability(clone_root: str | Path | None = None) -> BackendCheckResult:
    _ = clone_root

    try:
        _import_installed("DeepUnitMatch.testing.test")
        _import_installed("DeepUnitMatch.utils.param_fun")
    except Exception as exc:
        return BackendCheckResult(
            available=False,
            reason_code="import_failed",
            import_target="DeepUnitMatch",
            diagnostics={
                "message": str(exc),
                "traceback": traceback.format_exc(),
            },
        )

    return BackendCheckResult(
        available=True,
        reason_code="ok",
        import_target="DeepUnitMatch",
        diagnostics={
            "import_strategy": "installed_package",
        },
    )


def preprocess_waveforms(
    *,
    waveforms: Any,
    channel_pos: Any,
    matrix_unit_ids: list[str],
    prep_dir: str | Path,
    clone_root: str | Path | None = None,
) -> tuple[Path | None, dict[str, Any]]:
    _ = clone_root
    diag: dict[str, Any] = {
        "import_target": "DeepUnitMatch",
        "import_strategy": "installed_package",
        "import_ok": False,
        "reason": None,
    }

    try:
        import numpy as np

        param_mod = _import_installed("DeepUnitMatch.utils.param_fun")
        diag["import_ok"] = True

        prep_path = Path(prep_dir)
        session_id = np.zeros(len(matrix_unit_ids), dtype=int)
        param_mod.get_snippets(
            waveforms,
            [channel_pos],
            session_id,
            save_path=str(prep_path),
            unit_ids=np.array([int(u) for u in matrix_unit_ids]),
        )

        processed_dir = prep_path / "processed_waveforms"
        diag["processed_waveforms_dir"] = str(processed_dir)
        if not processed_dir.exists() or not processed_dir.is_dir():
            diag["reason"] = "processed_waveforms_missing"
            return None, diag

        session_dirs = [p for p in processed_dir.iterdir() if p.is_dir()]
        n_processed_files = sum(len([x for x in sd.glob("Unit*_RawSpikes.npy") if x.is_file()]) for sd in session_dirs)
        diag["processed_waveforms_n_sessions"] = int(len(session_dirs))
        diag["processed_waveforms_n_files"] = int(n_processed_files)
        if len(session_dirs) == 0 or n_processed_files == 0:
            diag["reason"] = "processed_waveforms_empty"
            return None, diag

        diag["reason"] = "ok"
        return processed_dir, diag
    except Exception as exc:
        diag["reason"] = "preprocessing_failed"
        diag["message"] = str(exc)
        diag["traceback"] = traceback.format_exc()
        return None, diag


def infer_similarity_matrix(
    *,
    processed_dir: str | Path,
    clone_root: str | Path | None = None,
    device: str = "cpu",
) -> BackendInferenceResult:
    _ = clone_root
    diag: dict[str, Any] = {
        "import_target": "DeepUnitMatch",
        "import_strategy": "installed_package",
        "import_ok": False,
        "reason": None,
    }

    try:
        test_mod = _import_installed("DeepUnitMatch.testing.test")
        diag["import_ok"] = True

        model = test_mod.load_trained_model(device=device)
        model = model.float()
        try:
            sim_matrix = test_mod.inference(model, str(processed_dir))
        except RuntimeError as exc:
            msg = str(exc)
            if "Input type (torch.FloatTensor) and weight type (torch.DoubleTensor)" not in msg:
                raise
            sim_matrix = test_mod.inference(model.float(), str(processed_dir))

        diag["reason"] = "ok"
        return BackendInferenceResult(
            success=True,
            reason_code="ok",
            import_ok=True,
            sim_matrix=sim_matrix,
            diagnostics=diag,
        )
    except Exception as exc:
        diag["reason"] = "inference_failed"
        diag["message"] = str(exc)
        diag["traceback"] = traceback.format_exc()
        return BackendInferenceResult(
            success=False,
            reason_code="inference_failed",
            import_ok=bool(diag.get("import_ok", False)),
            sim_matrix=None,
            diagnostics=diag,
        )
