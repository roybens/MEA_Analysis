from __future__ import annotations

import os
from pathlib import Path


def ensure_maxwell_hdf5_plugin_path(*, prefix: str = "[mea_analysis]") -> None:
    """Best-effort discovery of Maxwell HDF5 plugin directory."""

    env = os.environ.get("HDF5_PLUGIN_PATH")
    if env:
        try:
            if not Path(env).expanduser().exists():
                os.environ.pop("HDF5_PLUGIN_PATH", None)
        except Exception:
            pass

    if os.environ.get("HDF5_PLUGIN_PATH"):
        return

    here = Path(__file__).resolve()
    candidates: list[Path] = []
    for parent in [here] + list(here.parents):
        candidates.append(parent / "vendor" / "maxwell_hdf5_plugin" / "Linux")
        candidates.append(parent / "axon_reconstructor" / "vendor" / "maxwell_hdf5_plugin" / "Linux")

    for cand in candidates:
        if (cand / "libcompression.so").exists():
            os.environ["HDF5_PLUGIN_PATH"] = str(cand)
            print(f"{prefix}[DEBUG] set HDF5_PLUGIN_PATH={cand}", flush=True)
            return


__all__ = ["ensure_maxwell_hdf5_plugin_path"]
