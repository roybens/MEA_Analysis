# ==========================================================
# driver/file_scanner.py
# Utilities for discovering MEA data files under a directory tree.
# ==========================================================

import logging
from pathlib import Path

try:
    import h5py
    _HAS_H5PY = True
except ImportError:
    _HAS_H5PY = False


def find_data_files(path: Path, logger: logging.Logger | None = None) -> list[Path]:
    """Return a list of MEA data files (.h5/.nwb/.raw) under *path*.

    H5 files are discovered using the ``helper_functions`` module when
    available, so only files located inside a ``Network/`` subfolder are
    returned (matching the Maxwell export convention).  NWB and RAW files
    are discovered via a simple ``rglob``.

    Parameters
    ----------
    path:
        Root directory to scan.
    logger:
        Optional logger for warnings.

    Returns
    -------
    list[Path]
        Deduplicated, sorted list of discovered data files.
    """
    h5_files: list[Path] = []
    nwb_files: list[Path] = []
    raw_files: list[Path] = []

    try:
        from .. import helper_functions as helper
        h5_files = helper.find_files_with_subfolder(path, "data.raw.h5", "Network")
    except ImportError:
        try:
            import helper_functions as helper
            h5_files = helper.find_files_with_subfolder(path, "data.raw.h5", "Network")
        except Exception as e:
            if logger:
                logger.warning(f"Error finding HDF5 files via helper_functions: {e}")
    except Exception as e:
        if logger:
            logger.warning(f"Error finding HDF5 files: {e}")

    nwb_files = list(path.rglob("*.nwb"))
    raw_files = list(path.rglob("*.raw"))

    return h5_files + nwb_files + raw_files


def read_h5_recording_map(file_path: str | Path) -> dict[str, list[str]]:
    """Read the ``recordings/{rec}/{well}`` structure from an HDF5 file.

    Returns
    -------
    dict[str, list[str]]
        Mapping of recording name → list of well IDs.

    Raises
    ------
    ImportError
        If ``h5py`` is not installed.
    """
    if not _HAS_H5PY:
        raise ImportError("h5py is required to read HDF5 files")

    with h5py.File(str(file_path), "r") as h5f:
        return {
            recording: list(h5f["recordings"][recording].keys())
            for recording in h5f["recordings"].keys()
        }
