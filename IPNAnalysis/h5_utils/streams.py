from __future__ import annotations

from pathlib import Path


def list_stream_recording_names(*, h5_path: Path, stream_id: str) -> list[str]:
    """Return recording segment names under /wells/<stream_id>."""

    try:
        import h5py
    except Exception as e:  # pragma: no cover
        raise RuntimeError("listing stream recording names requires h5py") from e

    with h5py.File(h5_path, "r") as h5:
        return list(h5["wells"][stream_id].keys())


__all__ = ["list_stream_recording_names"]
