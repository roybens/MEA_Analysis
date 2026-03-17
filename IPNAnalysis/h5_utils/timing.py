from __future__ import annotations

from pathlib import Path


def read_well_rec_frame_nos_and_trigger_settings(*, h5_path: Path, stream_id: str, rec_name: str) -> dict:
    """Read per-recording timing data for Maxwell /wells/<stream>/<rec>."""

    try:
        import h5py
        import numpy as np
    except Exception as e:  # pragma: no cover
        raise RuntimeError("reading well timing requires h5py/numpy") from e

    with h5py.File(h5_path, "r") as h5:
        rec = h5["wells"][str(stream_id)][str(rec_name)]
        start_ms = int(np.asarray(rec["start_time"][()]).ravel()[0])
        stop_ms = int(np.asarray(rec["stop_time"][()]).ravel()[0])
        frame_nos = np.asarray(rec["groups"]["routed"]["frame_nos"], dtype=np.int64)
        sampling_raw = rec["settings"]["sampling"][()]
        sampling_hz = float(np.asarray(sampling_raw).ravel()[0])

    return {
        "start_ms": start_ms,
        "stop_ms": stop_ms,
        "frame_nos": frame_nos,
        "sampling_hz": sampling_hz,
    }


__all__ = ["read_well_rec_frame_nos_and_trigger_settings"]
