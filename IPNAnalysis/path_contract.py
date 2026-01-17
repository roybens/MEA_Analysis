"""Stable, dependency-free path contract for MEA_Analysis outputs.

Author: Adam Weiner
Datetime: 2026-01-16
Assisted by: ChatGPT-5.2

Purpose
-------
This module intentionally contains *only* lightweight path/metadata logic (no spikeinterface,
no torch, no matplotlib). It exists so that:

- The main pipeline can compute output locations deterministically.
- External callers (e.g., axon_reconstructor, Slurm wrappers) can compute/validate output
  paths without importing the heavy scientific stack.

If this contract ever changes, downstream consumers should only need to track changes here.
"""

from __future__ import annotations

import os
from pathlib import Path


def compute_relative_pattern(data_file: os.PathLike[str] | str) -> str:
    """Compute the per-run relative pattern used in MEA_Analysis outputs.

    Mirrors the logic historically implemented in
    `IPNAnalysis/mea_analysis_routine.py::MEAPipeline._parse_metadata`.

    The returned value does not include the well/stream id.
    """

    file_path = Path(data_file).resolve()

    # Legacy/initial default used by the pipeline
    try:
        relative_pattern = f"{file_path.parent.parent.name}/{file_path.parent.name}/{file_path.name}"
    except Exception:
        relative_pattern = str(file_path.name)

    # Path slicing override when there is enough depth.
    parts = str(file_path).split(os.sep)
    if len(parts) > 5:
        relative_pattern = os.path.join(*parts[-6:-1])

    return relative_pattern


def compute_output_dir(
    *,
    output_root: os.PathLike[str] | str,
    data_file: os.PathLike[str] | str,
    well: str,
) -> Path:
    """Compute per-well output directory under output_root."""

    output_root = Path(output_root).resolve()
    relative_pattern = compute_relative_pattern(data_file)
    return output_root / relative_pattern / well


def compute_sorter_output_dir(
    *,
    output_root: os.PathLike[str] | str,
    data_file: os.PathLike[str] | str,
    well: str,
) -> Path:
    """Compute the folder that contains sorter artifacts."""

    return compute_output_dir(output_root=output_root, data_file=data_file, well=well) / "sorter_output"
