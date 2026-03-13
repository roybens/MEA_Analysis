# ==========================================================
# pipeline/core.py
# MEAPipeline — the main pipeline class, assembled from mixins.
# ==========================================================

import gc
import shutil
from pathlib import Path

from .stages import ProcessingStage
from .defaults import _default_um_kwargs, _default_am_kwargs, _default_option_kwargs
from ._setup_mixin import SetupMixin
from ._checkpoint_mixin import CheckpointMixin
from ._preprocessing_mixin import PreprocessingMixin
from ._sorting_mixin import SortingMixin
from ._merge_mixin import MergeMixin
from ._analyzer_mixin import AnalyzerMixin
from ._reports_mixin import ReportsMixin
from ._burst_mixin import BurstMixin


class MEAPipeline(
    SetupMixin,
    CheckpointMixin,
    PreprocessingMixin,
    SortingMixin,
    MergeMixin,
    AnalyzerMixin,
    ReportsMixin,
    BurstMixin,
):
    """
    SOTA Pipeline for 2D MEA Analysis (Maxwell Biosystems).

    Encapsulates Preprocessing, Sorting, Merge, Analysis, and Curation in a
    resumable, checkpoint-aware pipeline.  Each phase is implemented as a
    dedicated mixin so the codebase stays navigable as the pipeline grows.

    Phases
    ------
    1. Preprocessing  (``run_preprocessing``)
    2. Spike Sorting  (``run_sorting``)
    2.5 Merge         (``run_optional_merge_phase``)
    3. Analyzer       (``run_analyzer``)
    4. Reports        (``generate_reports``)
    """

    def __init__(
        self,
        file_path,
        stream_id='well000',
        recording_num='rec0000',
        output_root=None,
        checkpoint_root=None,
        sorter='kilosort4',
        docker_image=None,
        verbose=True,
        cleanup=False,
        force_restart=False,
        n_jobs=None,
        chunk_duration=None,
        sorter_kwargs=None,
        um_kwargs=None,
        am_kwargs=None,
        option_kwargs=None,
    ):
        self.file_path = Path(file_path).resolve()
        self.stream_id = stream_id
        self.recording_num = recording_num
        self.sorter = sorter
        self.docker_image = docker_image
        self.verbose = verbose
        self.cleanup_flag = cleanup
        self.force_restart = force_restart

        # Merge / option kwargs with defaults
        self.um_kwargs = _default_um_kwargs()
        if isinstance(um_kwargs, dict):
            self.um_kwargs.update(um_kwargs)
        self.am_kwargs = _default_am_kwargs()
        if isinstance(am_kwargs, dict):
            self.am_kwargs.update(am_kwargs)
        self.option_kwargs = _default_option_kwargs()
        if isinstance(option_kwargs, dict):
            self.option_kwargs.update(option_kwargs)

        # Resource hints
        self.n_jobs = n_jobs
        self.chunk_duration = chunk_duration
        self.output_subdir_after_well = self._validate_output_subdir_after_well(
            self.option_kwargs.get("output_subdir_after_well")
        )
        self.sorter_kwargs = sorter_kwargs

        # UnitMatch settings
        self.unitmatch_merge_units = bool(self.um_kwargs.get("merge_units"))
        self.unitmatch_dry_run = bool(self.um_kwargs.get("dry_run"))
        self.unitmatch_scored_dry_run = bool(self.um_kwargs.get("scored_dry_run"))
        self.unitmatch_output_subdir_name = (
            self._validate_output_subdir_after_well(self.um_kwargs.get("output_subdir_name"))
            or "unitmatch_outputs"
        )
        self.unitmatch_throughput_subdir_name = (
            self._validate_output_subdir_after_well(self.um_kwargs.get("throughput_subdir_name"))
            or "unitmatch_throughput"
        )
        self.unitmatch_max_candidate_pairs = int(self.um_kwargs.get("max_candidate_pairs"))
        self.unitmatch_oversplit_min_probability = float(
            self.um_kwargs.get("oversplit_min_probability")
        )
        self.unitmatch_oversplit_max_suggestions = int(
            self.um_kwargs.get("oversplit_max_suggestions")
        )
        self.unitmatch_apply_merges = bool(self.um_kwargs.get("apply_merges"))
        self.unitmatch_recursive = bool(self.um_kwargs.get("recursive"))
        self.unitmatch_max_iterations = int(self.um_kwargs.get("max_iterations"))
        self.unitmatch_max_spikes_per_unit = int(self.um_kwargs.get("max_spikes_per_unit"))
        self.unitmatch_keep_all_iterations = bool(self.um_kwargs.get("keep_all_iterations"))
        self.unitmatch_generate_reports = bool(self.um_kwargs.get("generate_reports"))
        self.unitmatch_report_subdir_name = (
            self._validate_output_subdir_after_well(self.um_kwargs.get("report_subdir_name"))
            or "unitmatch_reports"
        )
        self.unitmatch_report_max_heatmap_units = int(
            self.um_kwargs.get("report_max_heatmap_units")
        )

        # Auto-merge settings
        self.auto_merge_units = bool(self.am_kwargs.get("enabled"))
        self.auto_merge_presets = self.am_kwargs.get("presets")
        self.auto_merge_steps_params = self.am_kwargs.get("steps_params")

        # Miscellaneous options
        self.force_rerun_analyzer = bool(self.option_kwargs.get("force_rerun_analyzer"))
        self.preprocessed_recording = self.option_kwargs.get("preprocessed_recording")
        self.skip_preprocessing = bool(self.option_kwargs.get("skip_preprocessing"))
        self.cuda_visible_devices = self.option_kwargs.get("cuda_visible_devices")

        # Parse metadata and build directory structure.
        self.metadata = self._parse_metadata()
        self.run_id = self.metadata.get('run_id', 'UnknownRun')
        self.project_name = self.metadata.get('project', 'UnknownProject')
        self.well = self.metadata.get('well', 'UnknownWell')
        self.chip_id = self.metadata.get('chip_id', 'UnknownChip')
        self.date = self.metadata.get('date', 'UnknownDate')
        self.relative_pattern = self.metadata.get('relative_pattern', 'UnknownPattern')

        self.output_root = Path(output_root)
        base_output_dir = Path(output_root) / self.relative_pattern / self.stream_id
        self.output_dir = (
            base_output_dir / self.output_subdir_after_well
            if self.output_subdir_after_well is not None
            else base_output_dir
        )
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Logger
        log_file = self.output_dir / f"{self.run_id}_{self.stream_id}_pipeline.log"
        self.logger = self._setup_logger(log_file)

        # Checkpointing
        ckpt_root = Path(checkpoint_root) if checkpoint_root else self.output_dir / "checkpoints"
        ckpt_root.mkdir(parents=True, exist_ok=True)
        self.checkpoint_file = (
            ckpt_root / f"{self.project_name}_{self.run_id}_{self.stream_id}_checkpoint.json"
        )
        self.state = self._load_checkpoint()

        # Apply runtime controls after logger exists so errors are visible.
        self._apply_runtime_controls()
        self._log_runtime_controls()

        # Data placeholders
        self.recording = None
        self.sorting = None
        self.analyzer = None

    # ------------------------------------------------------------------
    # Cleanup
    # ------------------------------------------------------------------

    def cleanup(self):
        if self.cleanup_flag:
            self.logger.info("Cleaning up temp files...")
            shutil.rmtree(self.output_dir / "binary", ignore_errors=True)
            shutil.rmtree(self.output_dir / "sorter_output", ignore_errors=True)
        self.recording = None
        self.sorting = None
        self.analyzer = None
        gc.collect()
