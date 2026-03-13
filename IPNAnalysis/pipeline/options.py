# ==========================================================
# pipeline/options.py
# MEARunOptions and MEARunResult dataclasses.
# ==========================================================

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class MEARunOptions:
    """Immutable configuration bundle for a single pipeline run."""

    file_path: str | Path
    stream_id: str
    recording_num: str = "rec0000"
    output_root: str | Path | None = None
    checkpoint_root: str | Path | None = None
    output_subdir_after_well: str | None = None
    sorter: str = "kilosort4"
    docker_image: str | None = None
    verbose: bool = False
    cleanup: bool = False
    force_restart: bool = False
    resume_from: str | None = None
    n_jobs: int | None = None
    chunk_duration: str | None = None
    sorter_kwargs: dict | None = None
    um_kwargs: dict | None = None
    am_kwargs: dict | None = None
    option_kwargs: dict | None = None
    reanalyze_bursts: bool = False
    skip_spikesorting: bool = False
    run_analyzer: bool = True
    run_reports: bool = True
    thresholds: dict | None = None
    no_curation: bool = False
    export_to_phy: bool = False
    plot_mode: str = "separate"
    plot_debug: bool = False
    raster_sort: str | None = None
    fixed_y: bool = False
    auto_merge_template_diff_thresh: str = "0.05,0.15,0.25"


@dataclass(frozen=True)
class MEARunResult:
    """Result returned by :func:`run_mea_pipeline`."""

    pipeline: object  # MEAPipeline — kept as object to avoid circular imports
    skipped: bool = False
    reanalyzed_bursts: bool = False
