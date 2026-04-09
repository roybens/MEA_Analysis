from __future__ import annotations

import importlib
import sys
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any


DEFAULT_CLONE_ROOT = Path("/home/adamm/dev/pkgs/UnitMatch/UnitMatchPy")
DEFAULT_DEEPUNITMATCH_DIR = DEFAULT_CLONE_ROOT / "DeepUnitMatch"


@dataclass(frozen=True)
class BackendCheckResult:
    available: bool
    reason_code: str
    clone_root: str
    clone_dir: str
    diagnostics: dict[str, Any]


@dataclass(frozen=True)
class BackendInferenceResult:
    success: bool
    reason_code: str
    import_ok: bool
    sim_matrix: Any | None
    diagnostics: dict[str, Any]


def _import_clone_only(module_name: str, clone_root: Path, clone_dir: Path) -> Any:
    if not clone_root.exists() or not clone_dir.exists():
        raise RuntimeError(
            f"UnitMatch clone not found. Expected: {clone_root} and {clone_dir}"
        )

    added: list[str] = []
    clone_paths = [str(clone_dir), str(clone_root)]
    try:
        for p in clone_paths:
            if p not in sys.path:
                sys.path.insert(0, p)
                added.append(p)
        mod = importlib.import_module(module_name)
    finally:
        for p in added:
            try:
                sys.path.remove(p)
            except ValueError:
                pass

    module_file = str(getattr(mod, "__file__", "") or "")
    if module_file and not module_file.startswith(str(clone_root)):
        raise RuntimeError(
            f"Module '{module_name}' was not imported from clone path. Imported from: {module_file}"
        )
    return mod


def check_backend_availability(clone_root: str | Path | None = None) -> BackendCheckResult:
    root = Path(clone_root) if clone_root is not None else DEFAULT_CLONE_ROOT
    deep_dir = root / "DeepUnitMatch"
    present = bool(root.exists() and deep_dir.exists())

    if not present:
        return BackendCheckResult(
            available=False,
            reason_code="clone_not_found",
            clone_root=str(root),
            clone_dir=str(deep_dir),
            diagnostics={
                "clone_present": False,
            },
        )

    try:
        _import_clone_only("DeepUnitMatch.testing.test", root, deep_dir)
        _import_clone_only("DeepUnitMatch.utils.param_fun", root, deep_dir)
    except Exception as exc:
        return BackendCheckResult(
            available=False,
            reason_code="import_failed",
            clone_root=str(root),
            clone_dir=str(deep_dir),
            diagnostics={
                "clone_present": True,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            },
        )

    return BackendCheckResult(
        available=True,
        reason_code="ok",
        clone_root=str(root),
        clone_dir=str(deep_dir),
        diagnostics={
            "clone_present": True,
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
    root = Path(clone_root) if clone_root is not None else DEFAULT_CLONE_ROOT
    deep_dir = root / "DeepUnitMatch"
    diag: dict[str, Any] = {
        "clone_root": str(root),
        "clone_dir": str(deep_dir),
        "import_ok": False,
        "reason": None,
    }

    try:
        import numpy as np

        param_mod = _import_clone_only("DeepUnitMatch.utils.param_fun", root, deep_dir)
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
    root = Path(clone_root) if clone_root is not None else DEFAULT_CLONE_ROOT
    deep_dir = root / "DeepUnitMatch"
    diag: dict[str, Any] = {
        "clone_root": str(root),
        "clone_dir": str(deep_dir),
        "import_ok": False,
        "reason": None,
    }

    try:
        test_mod = _import_clone_only("DeepUnitMatch.testing.test", root, deep_dir)
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
