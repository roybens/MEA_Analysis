# ==========================================================
# mea_analysis_routine.py
# Author: Mandar Patil
# Contributors: Yuxin Ren, Shruti Shah, Adam Weiner
# LLM Assisted Edits: Yes ChatGPT-4, Claude sonnet 4.6, ChatGPT-5.3-Codex
# ==========================================================

import os
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"
os.environ["QT_QPA_PLATFORM"] = "offscreen"  # For headless environments
import sys
import shutil
import gc
import json
import re
import traceback
import logging
import argparse
import configparser
from dataclasses import dataclass
from pathlib import Path
from datetime import datetime
from timeit import default_timer as timer
from enum import Enum
from math import floor
from typing import Any
import psutil
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg') # Non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf

from spikeinterface.sortingcomponents.peak_detection import detect_peaks

# Scientific Libraries
import spikeinterface.full as si
import spikeinterface.preprocessing as spre
import spikeinterface.curation as sic

# If this file is executed as a *script* (not imported as a package module),
# allow local imports by adding the script folder and repo root to sys.path.
# When imported as `MEA_Analysis.IPNAnalysis.mea_analysis_routine`, we avoid
# mutating sys.path and rely on proper package-relative imports.
if not __package__:
    script_dir = Path(__file__).resolve().parent
    if str(script_dir) not in sys.path:
        sys.path.append(str(script_dir))

    root_dir = script_dir.parent.parent
    if str(root_dir) not in sys.path:
        sys.path.append(str(root_dir))


# Import helper modules.
# - When imported as a package module, prefer relative imports.
# - When executed as a script, prefer local (same-folder) imports.
try:
    if __package__:
        from .parameter_free_burst_detector import compute_network_bursts
        from . import helper_functions as helper
        from .scalebury import add_scalebar
        from .config_loader import load_config, resolve_args
        from .UnitMatch.runner import run_unitmatch_merge_with_recursion, UnitMatchConfig
        from .UnitMatch.reporting import UnitMatchReportConfig, generate_unitmatch_static_report_pack
    else:
        from parameter_free_burst_detector import compute_network_bursts
        import helper_functions as helper
        from scalebury import add_scalebar
        from config_loader import load_config, resolve_args
        from UnitMatch.runner import run_unitmatch_merge_with_recursion, UnitMatchConfig
        from UnitMatch.reporting import UnitMatchReportConfig, generate_unitmatch_static_report_pack
except ImportError:
    try:
        from MEA_Analysis.IPNAnalysis.parameter_free_burst_detector import compute_network_bursts
        from MEA_Analysis.IPNAnalysis import helper_functions as helper
        from MEA_Analysis.IPNAnalysis.scalebury import add_scalebar
        from MEA_Analysis.IPNAnalysis.config_loader import load_config, resolve_args
        from MEA_Analysis.IPNAnalysis.UnitMatch.runner import run_unitmatch_merge_with_recursion, UnitMatchConfig
        from MEA_Analysis.IPNAnalysis.UnitMatch.reporting import (
            UnitMatchReportConfig,
            generate_unitmatch_static_report_pack,
        )
    except ImportError as e:
        raise ImportError(
            "Could not import MEA_Analysis.IPNAnalysis helper modules "
            "(parameter_free_burst_detector/helper_functions/scalebury/config_loader). "
            "If you are running this file directly, prefer: "
            "`python -m MEA_Analysis.IPNAnalysis.mea_analysis_routine ...`"
        ) from e



class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer): return int(obj)
        if isinstance(obj, np.floating): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        return super(NpEncoder, self).default(obj)
    
# --- Enums for Checkpointing ---
class ProcessingStage(Enum):
    NOT_STARTED = 0
    PREPROCESSING = 1
    PREPROCESSING_COMPLETE = 2
    SORTING = 3
    SORTING_COMPLETE = 4
    MERGE = 5
    MERGE_COMPLETE = 6
    ANALYZER = 7
    ANALYZER_COMPLETE = 8
    REPORTS = 9
    REPORTS_COMPLETE = 10


CHECKPOINT_SCHEMA_VERSION = 2



# --- The Main Pipeline Class ---
def _default_um_kwargs() -> dict[str, Any]:
    return {
        "merge_units": False,
        "dry_run": True,
        "scored_dry_run": True,
        "output_subdir_name": "unitmatch_outputs",
        "throughput_subdir_name": "unitmatch_throughput",
        "max_candidate_pairs": 20000,
        "oversplit_min_probability": 0.99,
        "oversplit_max_suggestions": 2000,
        "apply_merges": False,
        "recursive": False,
        "max_iterations": 5,
        "max_spikes_per_unit": 100,
        "keep_all_iterations": True,
        "generate_reports": True,
        "report_subdir_name": "unitmatch_reports",
        "report_max_heatmap_units": 200,
    }


def _default_am_kwargs() -> dict[str, Any]:
    return {
        "enabled": False,
        "presets": None,
        "steps_params": None,
        "template_diff_thresh": "0.05,0.15,0.25",
    }


def _default_option_kwargs() -> dict[str, Any]:
    return {
        "force_rerun_analyzer": False,
        "preprocessed_recording": None,
        "skip_preprocessing": False,
        "cuda_visible_devices": None,
        "output_subdir_after_well": None,
    }


class MEAPipeline:
    """
    SOTA Pipeline for 2D MEA Analysis (Maxwell Biosystems).
    Encapsulates Preprocessing, Sorting, Analysis, and Curation.
    """

    def __init__(self, file_path, stream_id='well000', recording_num='rec0000', output_root=None, checkpoint_root=None, 
                 sorter='kilosort4', docker_image=None, verbose=True, 
                 cleanup=False, force_restart=False,
                 n_jobs: int | None = None,
                 chunk_duration: str | None = None,
                 sorter_kwargs: dict | None = None,
                 um_kwargs: dict | None = None,
                 am_kwargs: dict | None = None,
                 option_kwargs: dict | None = None):
        
        self.file_path = Path(file_path).resolve()
        self.stream_id = stream_id
        self.recording_num = recording_num
        self.sorter = sorter
        self.docker_image = docker_image
        self.verbose = verbose
        self.cleanup_flag = cleanup
        self.force_restart = force_restart

        self.um_kwargs = _default_um_kwargs()
        if isinstance(um_kwargs, dict):
            self.um_kwargs.update(um_kwargs)
        self.am_kwargs = _default_am_kwargs()
        if isinstance(am_kwargs, dict):
            self.am_kwargs.update(am_kwargs)
        self.option_kwargs = _default_option_kwargs()
        if isinstance(option_kwargs, dict):
            self.option_kwargs.update(option_kwargs)

        # Optional resource hints (used by some SpikeInterface steps).
        self.n_jobs = n_jobs
        self.chunk_duration = chunk_duration
        self.output_subdir_after_well = self._validate_output_subdir_after_well(self.option_kwargs.get("output_subdir_after_well"))

        # Optional sorter kwargs override (e.g. Kilosort4 batch_size tuning).
        self.sorter_kwargs = sorter_kwargs

        # Optional UnitMatch alternative merge path (default off).
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
        self.unitmatch_oversplit_min_probability = float(self.um_kwargs.get("oversplit_min_probability"))
        self.unitmatch_oversplit_max_suggestions = int(self.um_kwargs.get("oversplit_max_suggestions"))
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
        self.unitmatch_report_max_heatmap_units = int(self.um_kwargs.get("report_max_heatmap_units"))

        # Optional post-sorting unit merge (default-off).
        self.auto_merge_units = bool(self.am_kwargs.get("enabled"))
        self.auto_merge_presets = self.am_kwargs.get("presets")
        self.auto_merge_steps_params = self.am_kwargs.get("steps_params")

        # Allow re-running analyzer without re-running spikesorting.
        self.force_rerun_analyzer = bool(self.option_kwargs.get("force_rerun_analyzer"))

        # Optional preprocessed recording injection path.
        self.preprocessed_recording = self.option_kwargs.get("preprocessed_recording")
        self.skip_preprocessing = bool(self.option_kwargs.get("skip_preprocessing"))

        # Optional runtime/resource controls (default behavior when unset).
        self.cuda_visible_devices = self.option_kwargs.get("cuda_visible_devices")

        # 1. Parse Metadata & Paths
        self.metadata = self._parse_metadata()
        self.run_id = self.metadata.get('run_id', 'UnknownRun')
        self.project_name = self.metadata.get('project', 'UnknownProject')
        self.well = self.metadata.get('well', 'UnknownWell')
        self.chip_id = self.metadata.get('chip_id', 'UnknownChip')
        self.date = self.metadata.get('date', 'UnknownDate')

        
        # Define Directory Structure: Output / Pattern / Well
        self.relative_pattern = self.metadata.get('relative_pattern', 'UnknownPattern')
        # Default output_root: "AnalyzedData" sibling of the recording file directory
        effective_output_root = Path(output_root) if output_root is not None else (self.file_path.parent / "AnalyzedData")
        self.output_root = effective_output_root
        base_output_dir = effective_output_root / self.relative_pattern / self.stream_id
        self.output_dir = (
            base_output_dir / self.output_subdir_after_well
            if self.output_subdir_after_well is not None
            else base_output_dir
        )
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 2. Logger Setup
        log_file = self.output_dir / f"{self.run_id}_{self.stream_id}_pipeline.log"
        self.logger = self._setup_logger(log_file)
        
        # 3. Checkpointing Setup
        ckpt_root = Path(checkpoint_root) if checkpoint_root else self.output_dir / "checkpoints"
        ckpt_root.mkdir(parents=True, exist_ok=True)
        self.checkpoint_file = ckpt_root / f"{self.project_name}_{self.run_id}_{self.stream_id}_checkpoint.json"
        self.state = self._load_checkpoint()

        # Apply optional runtime controls after logger exists so failures are visible.
        self._apply_runtime_controls()
        self._log_runtime_controls()

        # Data placeholders
        self.recording = None
        self.sorting = None
        self.analyzer = None

    # --- Setup Methods ---
    def _setup_logger(self, log_file):
        logger = logging.getLogger(f"mea_{self.stream_id}")
        logger.setLevel(logging.DEBUG if self.verbose else logging.INFO)
        # Prevent duplicate handlers if running in loop
        if not logger.handlers:
            formatter = logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s')
            fh = logging.FileHandler(log_file, mode='a') # Append to log per run
            #add a big header
            fh.stream.write("\n" + "="*80 + "\n")
            fh.setFormatter(formatter)
            logger.addHandler(fh)
            ch = logging.StreamHandler(sys.stdout)
            ch.setFormatter(formatter)
            logger.addHandler(ch)
        return logger

    def _apply_runtime_controls(self):
        if self.cuda_visible_devices is not None:
            try:
                os.environ["CUDA_VISIBLE_DEVICES"] = str(self.cuda_visible_devices)
            except Exception:
                pass

    def _log_runtime_controls(self):
        def _env_or_none(name):
            value = os.environ.get(name)
            if value is None:
                return None
            token = str(value).strip()
            return token if token else None

        self.logger.info(
            "Runtime snapshot: pid=%s cpu_count=%s n_jobs=%s chunk_duration=%s",
            os.getpid(),
            os.cpu_count(),
            self.n_jobs,
            self.chunk_duration,
        )
        self.logger.info(
            "Runtime controls: cuda_visible_devices=%s",
            self.cuda_visible_devices,
        )
        self.logger.info(
            "Runtime env effective: CUDA_VISIBLE_DEVICES=%s",
            _env_or_none("CUDA_VISIBLE_DEVICES"),
        )

    def _parse_metadata(self):
    
        meta = {
            'run_id': None, 'chip_id': None, 'project': None, 
            'relative_pattern': f"{self.file_path.parent.parent.name}/{self.file_path.parent.name}/{self.file_path.name}",
            'date': None, 'well': None
        }
        
        # Strategy A: Regex on path (Fallback)
        try:
            path_str = str(self.file_path)
            match = re.search(r"/(\d+)/data.raw.h5", path_str)
            if match: meta['run_id'] = match.group(1)
            parts = path_str.split(os.sep)
            if len(parts) > 5:
                meta['relative_pattern'] = os.path.join(*parts[-6:-1])
                meta['project'] = parts[-6]
                meta['date'] = parts[-5]
                meta['chip_id'] = parts[-4]
                meta['well'] = self.stream_id
        except Exception: pass

        # Strategy B: .metadata file (Overrides regex)
        meta_file = self.file_path.parent / ".metadata"
        if meta_file.exists():
            try:
                cfg = configparser.ConfigParser()
                cfg.read(meta_file, encoding='utf-8') 
                if 'properties' in cfg:
                    meta['run_id'] = cfg['properties'].get('runid', meta.get('run_id'))
                    meta['project'] = cfg['properties'].get('project_title', meta.get('project'))
                if 'runtime' in cfg:
                    meta['chip_id'] = cfg['runtime'].get('chipid', meta.get('chip_id'))
            except: pass
        return meta

    def _validate_output_subdir_after_well(self, value):
        if value is None:
            return None

        token = str(value).strip()
        if not token:
            return None

        if "/" in token or "\\" in token:
            raise ValueError(
                "output_subdir_after_well must be a single directory name, not a path"
            )

        candidate = Path(token)
        if candidate.is_absolute() or token in (".", ".."):
            raise ValueError(
                "output_subdir_after_well must be a relative single directory name"
            )

        return token

    # --- Checkpointing Methods ---
    def _load_checkpoint(self):
        if self.checkpoint_file.exists() and not self.force_restart:
            with open(self.checkpoint_file, 'r') as f:
                state = json.load(f)

            # Migrate legacy stage numbering (schema v1 had no merge stage).
            try:
                schema_version = int(state.get("checkpoint_schema_version", 1))
            except Exception:
                schema_version = 1
            if schema_version < CHECKPOINT_SCHEMA_VERSION:
                try:
                    old_stage = int(state.get("stage", ProcessingStage.NOT_STARTED.value))
                except Exception:
                    old_stage = ProcessingStage.NOT_STARTED.value
                # v1 mapping: ANALYZER/ANALYZER_COMPLETE/REPORTS/REPORTS_COMPLETE
                # were 5/6/7/8; shift by +2 to account for MERGE stages.
                if old_stage >= 5:
                    state["stage"] = old_stage + 2
                state["checkpoint_schema_version"] = CHECKPOINT_SCHEMA_VERSION

            return state
        return {'stage': ProcessingStage.NOT_STARTED.value, #last stage completed
                'checkpoint_schema_version': CHECKPOINT_SCHEMA_VERSION,
                'failed_stage': None,
                'last_updated': None,
                'run_id': self.run_id,
                'chip_id': self.chip_id,
                'well': self.well,
                'project': self.project_name,
                'date': self.date,
                'output_dir': str(self.output_dir),
                'error': None,
                }

    def _save_checkpoint(self, stage, **kwargs):
        self.state['stage'] = stage.value
        self.state['checkpoint_schema_version'] = CHECKPOINT_SCHEMA_VERSION
        self.state['last_updated'] = str(datetime.now())
        self.state.update(kwargs)
        with open(self.checkpoint_file, 'w') as f: json.dump(self.state, f, indent=2)
        self.logger.info(f"Checkpoint Saved: {stage.name}")

    def should_skip(self):
        if self.state['stage'] == ProcessingStage.REPORTS_COMPLETE.value and not self.force_restart:
            self.logger.info("Pipeline already completed. Skipping.")
            return True
        return False

    def _load_recording_file(self):
        fpath = str(self.file_path)
        if fpath.endswith(".h5"):
            # Maxwell loader does NOT accept load_if_exists
            return si.read_maxwell(fpath, stream_id=self.stream_id, rec_name=self.recording_num)
        elif fpath.endswith(".nwb"):
            return si.read_nwb(fpath, load_if_exists=True)
        elif self.file_path.is_dir():
            return si.load_extractor(self.file_path)
        raise ValueError(f"Unknown format: {fpath}")

    # --- Phase 1: Preprocessing ---
    def run_preprocessing(self):
        binary_folder = self.output_dir / "binary"

        # Optional shortcut for callers that provide a preprocessed recording.
        if self.preprocessed_recording is not None:
            self.logger.info(
                "Using injected preprocessed recording (skip_preprocessing=%s)",
                bool(self.skip_preprocessing),
            )
            self.recording = self.preprocessed_recording
            # Do not rewind checkpoint progress when resuming from later stages.
            if int(self.state.get('stage', 0)) < ProcessingStage.PREPROCESSING_COMPLETE.value:
                self._save_checkpoint(
                    ProcessingStage.PREPROCESSING_COMPLETE,
                    failed_stage=None,
                    error=None,
                    preprocessing_source="injected_preprocessed_recording",
                )
            return

        # 1. Resume Check
        if self.state['stage'] >= ProcessingStage.PREPROCESSING_COMPLETE.value and binary_folder.exists():
            self.logger.info("Resuming: Loading preprocessed data from binary cache.")
            try: self.recording = si.load(binary_folder)
            except: self.recording = si.load_extractor(binary_folder)
            return 
        self._save_checkpoint(ProcessingStage.PREPROCESSING)
        self.logger.info("--- [Phase 1] Preprocessing ---")
        
        # 2. Universal Loading
        rec = self._load_recording_file()

        # 3. Time Slicing (Remove last 1s to avoid end-of-file artifacts)
        fs = rec.get_sampling_frequency()
        self.metadata['fs'] = fs
        total_frames = rec.get_num_frames()
        # Calculate new end frame flooring to nearest integer second
        end_frame = floor(total_frames)
        #end_frame = total_frames - int(1.0 * fs)
        if end_frame > 0:
            self.logger.info(f"Trimming recording: {total_frames} -> {end_frame} frames (removed last 1s).")
            rec = rec.frame_slice(start_frame=0, end_frame=end_frame)

        if rec.get_dtype().kind == 'u':
            rec = spre.unsigned_to_signed(rec)


        rec = spre.highpass_filter(rec, freq_min=300)

        # Local Common Median Reference (Preserves Network Bursts)
        # NOTE: I'm not 100% sure because its a bit annoying to test, but runing this with local_radius (250, 250) may be failing every time.
        # If our intent is to get all channels within 250um, we should set local_radius=(0, 250). According to my cursory understanding of the docs, (250, 250) would create
        # an annulus excluding the inner 250um radius. 
        # -- aw 2026-02-01 20:44:37
        try:
            rec = spre.common_reference(rec, reference='local', operator='median', local_radius=(250, 250))
        except:
            self.logger.warning("Local CMR failed (missing locations?), using Global CMR.")
            rec = spre.common_reference(rec, reference='global', operator='median')
        
        rec.annotate(is_filtered=True)

        # F. Force Float32 (Prevents Signal Crushing)
        # Int16 scaling caused Kilosort crash on artifacts. Float32 is safer.
        if rec.get_dtype() != 'float32':
            self.logger.info("Converting to float32 to preserve signal fidelity...")
            rec = spre.astype(rec, 'float32')

        # 5. Save to BINARY
        if binary_folder.exists():
            shutil.rmtree(binary_folder)
        
        self.logger.info(f"Saving binary recording to {binary_folder}...")
        
        rec.save(
            folder=binary_folder, 
            format='binary', 
            overwrite=True, 
            n_jobs=(int(self.n_jobs) if self.n_jobs is not None else 16), 
            chunk_duration=(str(self.chunk_duration) if self.chunk_duration is not None else '1s'), 
            progress_bar=self.verbose
        )
        
        # 6. Reload from Binary (Memory Map)
        self.recording = si.load(binary_folder)
        self._save_checkpoint(ProcessingStage.PREPROCESSING_COMPLETE)

    # --- Phase 2: Sorting ---
    def run_sorting(self):
        sorter_folder = self.output_dir / "sorter_output"
        if self.state['stage'] >= ProcessingStage.SORTING_COMPLETE.value:
            self.logger.info("Resuming: Loading existing sorting.")
            try: self.sorting = si.read_sorter_folder(sorter_folder)
            except: self.sorting = si.read_kilosort(sorter_folder)
            return
        self._save_checkpoint(ProcessingStage.SORTING)
        self.logger.info(f"--- [Phase 2] Spike Sorting ({self.sorter}) ---")
        
        import torch
        if torch.cuda.is_available():
            total_vram = torch.cuda.get_device_properties(0).total_memory / (1024**3) # Convert to GB
            self.logger.info(f"GPU Detected: {torch.cuda.get_device_name(0)} with {total_vram:.2f} GB VRAM")
        else:
            self.logger.warning("No GPU detected! Kilosort4 will likely fail or run extremely slowly on CPU.")
            total_vram = 0

        fs = int(self.recording.get_sampling_frequency())

        # Parameters for 2D Static Cultures (Kilosort4)
        # Th_universal lowered to 8 to improve detection on organoids
        ks_params_high_vram ={
                        'batch_size': int(self.recording.get_sampling_frequency()) * 2,
                        'clear_cache': True,
                        'invert_sign': True,
                        'cluster_downsampling': 20,
                        'max_cluster_subset': None,
                        'nblocks':0,
                        'dmin':17,
                        'do_correction': False,
                    }
        ks_params_low_vram = {
                        'batch_size': int(self.recording.get_sampling_frequency() * 0.5),
                        'clear_cache': True,
                        'invert_sign': True,
                        'cluster_downsampling': 30, 
                        'max_cluster_subset': 50000,
                        'nblocks':0,
                        'do_correction': False,
        }

        if total_vram >= 14:
            ks_params = ks_params_high_vram
        else:
            ks_params = ks_params_low_vram

        # Optional override from caller (e.g. debug harness). This lets us tune
        # Kilosort parameters (like batch_size) without forking MEA_Analysis.
        if getattr(self, "sorter_kwargs", None):
            try:
                ks_params = dict(ks_params)
                ks_params.update(dict(self.sorter_kwargs))
            except Exception:
                pass

        start = timer()
        try:
            self.logger.info(f"Running Kilosort4 with parameters: {ks_params}")
            torch.cuda.empty_cache()
            self.sorting = si.run_sorter(
                sorter_name=self.sorter,
                recording=self.recording,
                folder=sorter_folder,
                delete_output_folder=False,
                remove_existing_folder=True, # Safety flag
                verbose=self.verbose,
                docker_image=self.docker_image, 
                **ks_params
            )
            self.logger.info(f"Sorting finished in {timer()-start:.2f}s")

            # [FIX] CRITICAL: Remove spikes that extend beyond the recording length
            # This prevents the IndexError during Analyzer computation
            self.logger.info("Cleaning sorting (removing excess spikes)...")
            self.sorting = si.remove_excess_spikes(self.sorting, self.recording)
            self.sorting = self.sorting.remove_empty_units()
            self._save_checkpoint(ProcessingStage.SORTING_COMPLETE, failed_stage=None, error=None)
        except Exception as e:
            err = {
                "failed_stage": ProcessingStage.SORTING.name,
                "exception": type(e).__name__,
                "message": str(e),
                "traceback": traceback.format_exc(),
                "time": str(datetime.now())
            }
            self.logger.error(err["traceback"])
            self._save_checkpoint(ProcessingStage.PREPROCESSING_COMPLETE, error=err)
            raise

    def _spike_detection_only(self):
        """Detects spikes (threshold crossings) without sorting."""
        self.logger.info("--- [Phase 2-Alt] Spike Detection (No Sorting) ---")
        job_kwargs = {
            'n_jobs': (int(self.n_jobs) if self.n_jobs is not None else 16),
            'chunk_duration': (str(self.chunk_duration) if self.chunk_duration is not None else '1s'),
            'progress_bar': self.verbose,
        }
        
        peaks = detect_peaks(
            self.recording, method='by_channel', 
            detect_threshold=5, peak_sign='neg', exclude_sweep_ms=0.1, **job_kwargs
        )

        
        self.logger.info(f"Detected {len(peaks)} total spikes.")
        fs = self.recording.get_sampling_frequency()
        self.metadata['fs'] = fs
        channel_ids = self.recording.get_channel_ids()
        spike_times = {}

        for ch_index, ch_id in enumerate(channel_ids):
            mask = peaks['channel_index'] == ch_index
            spike_times[ch_id] = (
                peaks['sample_index'][mask] / fs
                if np.any(mask)
                else np.array([])
            )
        
        np.save(self.output_dir / "spike_times.npy", spike_times)

        return list(spike_times.keys())

    # --- Phase 2.5: Merge (Optional) ---
    def run_optional_merge_phase(self):
        self.logger.info("--- [Phase 2.5] Merge (Optional) ---")

        try:
            current_stage = ProcessingStage(int(self.state.get('stage', ProcessingStage.SORTING_COMPLETE.value)))
        except Exception:
            current_stage = ProcessingStage.SORTING_COMPLETE

        analyzer_folder = self.output_dir / "analyzer_output"
        if (
            (not self.force_rerun_analyzer)
            and int(self.state.get('stage', 0)) >= ProcessingStage.ANALYZER_COMPLETE.value
            and analyzer_folder.exists()
        ):
            self.logger.info("Merge step skipped: analyzer already complete and rerun not requested")
            self._save_checkpoint(
                current_stage,
                merge_phase={
                    "status": "skipped_analyzer_complete",
                    "mode": None,
                    "applied": False,
                },
            )
            return

            self._save_checkpoint(ProcessingStage.MERGE)

        if self.unitmatch_merge_units:
            try:
                if self.auto_merge_units:
                    self.logger.warning(
                        "Both UnitMatch and auto_merge_units requested; using UnitMatch path and skipping auto_merge_units"
                    )
                merged_sorting, um_summary = run_unitmatch_merge_with_recursion(
                    sorting=self.sorting,
                    recording=self.recording,
                    output_dir=self.output_dir,
                    logger=self.logger,
                    config=UnitMatchConfig(
                        enabled=True,
                        dry_run=self.unitmatch_dry_run,
                        scored_dry_run=self.unitmatch_scored_dry_run,
                        fail_open=True,
                        max_candidate_pairs=self.unitmatch_max_candidate_pairs,
                        output_subdir_name=self.unitmatch_output_subdir_name,
                        throughput_subdir_name=self.unitmatch_throughput_subdir_name,
                        oversplit_min_probability=self.unitmatch_oversplit_min_probability,
                        oversplit_max_suggestions=self.unitmatch_oversplit_max_suggestions,
                        apply_merges=self.unitmatch_apply_merges,
                        recursive=self.unitmatch_recursive,
                        max_iterations=self.unitmatch_max_iterations,
                        max_spikes_per_unit=self.unitmatch_max_spikes_per_unit,
                        keep_all_iterations=self.unitmatch_keep_all_iterations,
                    ),
                )
                self.sorting = merged_sorting
                self.logger.info("UnitMatch summary: status=%s", um_summary.get("status"))

                artifacts = um_summary.setdefault("artifacts", {})
                merge_application = um_summary.setdefault("merge_application", {})
                final_merged_sorting_folder = artifacts.get("final_merged_sorting_folder")
                persisted_final_merged_sorting = bool(merge_application.get("persisted_final_merged_sorting", False))

                summary_path_s = str(artifacts.get("summary_json") or "")
                if summary_path_s:
                    try:
                        summary_path = Path(summary_path_s)
                        summary_path.parent.mkdir(parents=True, exist_ok=True)
                        summary_path.write_text(json.dumps(um_summary, indent=2), encoding="utf-8")
                    except Exception as summary_sync_exc:
                        self.logger.warning(
                            "Failed to sync UnitMatch summary JSON after merge artifact update (%s)",
                            summary_sync_exc,
                        )

                self._save_checkpoint(
                    ProcessingStage.MERGE_COMPLETE,
                    merge_phase={
                        "status": str(um_summary.get("status", "unknown")),
                        "mode": "unitmatch",
                        "applied": bool(um_summary.get("merge_application", {}).get("applied", False)),
                        "final_merged_sorting_folder": (
                            str(final_merged_sorting_folder) if final_merged_sorting_folder is not None else None
                        ),
                        "persisted_final_merged_sorting": persisted_final_merged_sorting,
                        "recursive": bool(self.unitmatch_recursive),
                        "n_iterations": int(
                            um_summary.get("iteration_convergence_report", {}).get("n_iterations_executed", 1)
                        ),
                        "summary": um_summary,
                    },
                )
                if self.unitmatch_generate_reports:
                    try:
                        generate_unitmatch_static_report_pack(
                            output_dir=self.output_dir,
                            logger=self.logger,
                            config=UnitMatchReportConfig(
                                throughput_subdir_name=str(self.unitmatch_throughput_subdir_name),
                                output_subdir_name=str(self.unitmatch_output_subdir_name),
                                report_subdir_name=str(self.unitmatch_report_subdir_name),
                                max_heatmap_units=int(self.unitmatch_report_max_heatmap_units),
                            ),
                        )
                    except Exception as report_exc:
                        self.logger.warning(
                            "UnitMatch static report generation failed-open; continuing (%s)",
                            report_exc,
                        )
            except Exception as e:
                self.logger.warning("UnitMatch merge path failed; continuing without UnitMatch merge (%s)", e)
                self._save_checkpoint(
                    ProcessingStage.MERGE_COMPLETE,
                    merge_phase={
                        "status": "error",
                        "mode": "unitmatch",
                        "applied": False,
                        "error": str(e),
                    },
                )
            return

        if self.auto_merge_units:
            merge_analyzer_folder = self.output_dir / "analyzer_output_merge_tmp"
            try:
                presets = self.auto_merge_presets
                steps_params = self.auto_merge_steps_params

                if presets is None or steps_params is None:
                    template_diff_thresh = [0.05, 0.15, 0.25]
                    presets = ["x_contaminations"] * len(template_diff_thresh)
                    steps_params = [
                        {"template_similarity": {"template_diff_thresh": float(t)}}
                        for t in template_diff_thresh
                    ]

                if merge_analyzer_folder.exists():
                    shutil.rmtree(merge_analyzer_folder)
                merge_sparsity = si.estimate_sparsity(
                    self.sorting,
                    self.recording,
                    method="radius",
                    radius_um=50,
                    peak_sign='neg',
                )
                merge_analyzer = si.create_sorting_analyzer(
                    self.sorting,
                    self.recording,
                    format="binary_folder",
                    folder=merge_analyzer_folder,
                    sparsity=merge_sparsity,
                    return_in_uV=True,
                )

                merge_exts = ["random_spikes", "templates", "template_similarity", "correlograms"]
                merge_kwargs = {'verbose': self.verbose}
                if self.n_jobs is not None:
                    merge_kwargs['n_jobs'] = int(self.n_jobs)
                merge_analyzer.compute(merge_exts, **merge_kwargs)

                n_before = int(merge_analyzer.sorting.get_num_units())
                try:
                    merged = sic.auto_merge_units(
                        merge_analyzer,
                        presets=presets,
                        steps_params=steps_params,
                        recursive=True,
                        n_jobs=(int(self.n_jobs) if self.n_jobs is not None else 1),
                    )
                except TypeError:
                    merged = sic.auto_merge_units(
                        merge_analyzer,
                        presets=presets,
                        steps_params=steps_params,
                        recursive=True,
                    )

                if isinstance(merged, tuple) and len(merged) >= 1:
                    merged = merged[0]

                if hasattr(merged, "sorting"):
                    merged_sorting = merged.sorting
                else:
                    merged_sorting = merged

                n_after = int(merged_sorting.get_num_units())
                self.logger.info("Auto-merge units: %d -> %d", n_before, n_after)
                self.sorting = merged_sorting
                self._save_checkpoint(
                    ProcessingStage.MERGE_COMPLETE,
                    merge_phase={
                        "status": "ok",
                        "mode": "auto_merge",
                        "applied": True,
                        "n_before": n_before,
                        "n_after": n_after,
                    },
                )
            except Exception as e:
                self.logger.warning("Auto-merge units failed; continuing without merging (%s)", e)
                self._save_checkpoint(
                    ProcessingStage.MERGE_COMPLETE,
                    merge_phase={
                        "status": "error",
                        "mode": "auto_merge",
                        "applied": False,
                        "error": str(e),
                    },
                )
            finally:
                try:
                    if merge_analyzer_folder.exists():
                        shutil.rmtree(merge_analyzer_folder)
                except Exception:
                    pass
            return

        self.logger.info("Merge step disabled; continuing without unit merging")
        self._save_checkpoint(
            ProcessingStage.MERGE_COMPLETE,
            merge_phase={
                "status": "skipped_disabled",
                "mode": None,
                "applied": False,
            },
        )

    # --- Phase 3: Analyzer ---
    def run_analyzer(self):
        analyzer_folder = self.output_dir / "analyzer_output"
        if (
            (not self.force_rerun_analyzer)
            and self.state['stage'] >= ProcessingStage.ANALYZER_COMPLETE.value
            and analyzer_folder.exists()
        ):
            self.logger.info("Resuming: Loading Sorting Analyzer.")
            self.analyzer = si.load_sorting_analyzer(analyzer_folder)
            # Keep sorting/analyzer unit IDs consistent for downstream reports.
            self.sorting = self.analyzer.sorting
            return
        
        self._save_checkpoint(ProcessingStage.ANALYZER)
        self.logger.info("--- [Phase 3] Computing Sorting Analyzer ---")
        try:
            if analyzer_folder.exists():
                shutil.rmtree(analyzer_folder)

            sparsity = si.estimate_sparsity(
                self.sorting, self.recording,
                method="radius", radius_um=50, peak_sign='neg'
            )

            # Build final analyzer once on current canonical sorting.
            self.analyzer = si.create_sorting_analyzer(
                self.sorting,
                self.recording,
                format="binary_folder",
                folder=analyzer_folder,
                sparsity=sparsity,
                return_in_uV=True,
            )

            ext_list = ["random_spikes","spike_amplitudes", "waveforms", "templates", "noise_levels", 
                        "quality_metrics", "template_metrics", "unit_locations"]
            
            ext_params = {
                "waveforms": {"ms_before": 1.0, "ms_after": 2.0},
                "unit_locations": {"method": "monopolar_triangulation"}
            }
            
            compute_kwargs = {'verbose': self.verbose}
            if self.n_jobs is not None:
                compute_kwargs['n_jobs'] = int(self.n_jobs)
            self.analyzer.compute(ext_list, extension_params=ext_params, **compute_kwargs)

            # Always use analyzer sorting as the source of truth going forward.
            self.sorting = self.analyzer.sorting
            self._save_checkpoint(ProcessingStage.ANALYZER_COMPLETE,failed_stage=None, error=None)
        except Exception as e:
            err = {
                "failed_stage": ProcessingStage.ANALYZER.name,
                "exception": type(e).__name__,
                "message": str(e),
                "traceback": traceback.format_exc(),
                "time": str(datetime.now())
            }
            self.logger.error(err["traceback"])
            try:
                fallback_stage = ProcessingStage(int(self.state.get('stage', ProcessingStage.MERGE_COMPLETE.value)))
            except Exception:
                fallback_stage = ProcessingStage.MERGE_COMPLETE
            self._save_checkpoint(fallback_stage, error=err)
            raise

    # --- Phase 4: Reports & Curation ---
    def generate_reports(self, thresholds=None, no_curation=False, export_phy=False,plot_mode="separate", plot_debug=False, raster_sort=None, fixed_y=False):
        if self.state['stage'] == ProcessingStage.REPORTS_COMPLETE.value: return

        self.logger.info("--- [Phase 4] Reports & Curation ---")
        try:
            q_metrics = self.analyzer.get_extension("quality_metrics").get_data()
            t_metrics = self.analyzer.get_extension("template_metrics").get_data()
            locations = self.analyzer.get_extension("unit_locations").get_data()
            
            q_metrics['loc_x'] = locations[:, 0]
            q_metrics['loc_y'] = locations[:, 1]
            
            q_metrics.to_excel(self.output_dir / "qm_unfiltered.xlsx")
            t_metrics.to_excel(self.output_dir / "tm_unfiltered.xlsx")
            self._plot_probe_locations(q_metrics.index.values, locations, "locations_unfiltered.pdf")

            if no_curation:
                self.logger.info("Skipping curation.")
                clean_units = q_metrics.index.values
            else:
                self.logger.info("Applying curation.")
                clean_metrics, rejection_log = self._apply_curation_logic(q_metrics, thresholds)
                clean_units = clean_metrics.index.values
                clean_metrics.to_excel(self.output_dir / "metrics_curated.xlsx")
                rejection_log.to_excel(self.output_dir / "rejection_log.xlsx")
                t_metrics.loc[clean_units].to_excel(self.output_dir / "tm_curated.xlsx")

            if len(clean_units) == 0:
                self.logger.warning("No units passed curation.")
                self._save_checkpoint(ProcessingStage.REPORTS_COMPLETE, n_units=0)
                return

            mask = np.isin(self.analyzer.unit_ids, clean_units)
            self._plot_probe_locations(clean_units, locations[mask], f"locations_{len(clean_units)}_units.pdf")
            self._plot_waveforms_grid(clean_units)
            self._run_burst_analysis(clean_units,plot_mode=plot_mode, plot_debug=plot_debug, raster_sort=raster_sort, fixed_y=fixed_y)

            if export_phy:
                si.export_to_phy(self.analyzer.select_units(clean_units), 
                            output_folder=self.output_dir/"phy_output", remove_if_exists=True, copy_binary=False)

            self._save_checkpoint(ProcessingStage.REPORTS_COMPLETE, n_units=len(clean_units),failed_stage=None, error=None)
        except Exception as e:
            err = {
                "failed_stage": ProcessingStage.REPORTS.name,
                "exception": type(e).__name__,
                "message": str(e),
                "traceback": traceback.format_exc(),
                "time": str(datetime.now())
            }
            self.logger.error(err["traceback"])
            self._save_checkpoint(ProcessingStage.ANALYZER_COMPLETE, error=err)
            raise

    def _apply_curation_logic(self, metrics, user_thresholds):
        defaults = {'presence_ratio': 0.75, 'rp_contamination': 0.15, 'firing_rate': 0.05, 
                    'amplitude_median': -20, 'amplitude_cv_median': 0.5}
        if user_thresholds: defaults.update(user_thresholds)
        
        keep_mask = np.ones(len(metrics), dtype=bool)
        rejections = []
        for idx, row in metrics.iterrows():
            reasons = []
            if row.get('presence_ratio', 1) < defaults['presence_ratio']: reasons.append("Low Presence")
            if row.get('rp_contamination', 0) > defaults['rp_contamination']: reasons.append("High Contam")
            if row.get('firing_rate', 0) < defaults['firing_rate']: reasons.append("Low FR")
            if row.get('amplitude_median', -100) > defaults['amplitude_median']: reasons.append("Low Amp")
            #TODO : need to add cv_median logic after checking if metric exists in current version of SI
            
            if reasons:
                keep_mask[metrics.index.get_loc(row.name)] = False
                rejections.append({"unit_id": row.name, "reasons": "; ".join(reasons)})
        
        return metrics[keep_mask], pd.DataFrame(rejections)

    def _plot_probe_locations(self, unit_ids, locations, filename):
        fig, ax = plt.subplots(figsize=(10.5, 6.5))
        si.plot_probe_map(self.recording, ax=ax, with_channel_ids=False)
        ax.scatter(locations[:, 0], locations[:, 1], s=10, c='blue', alpha=0.6)
        ax.invert_yaxis()
        fig.savefig(self.output_dir / filename)
        plt.close(fig)


    def _plot_waveforms_grid(self, unit_ids):
        pdf_path = self.output_dir / "waveforms_grid.pdf"
        self.logger.info(f"Generating PDF: {pdf_path}")
        
        wf_ext = self.analyzer.get_extension("waveforms")
        
        # [1] Get sampling freq to calculate milliseconds
        fs = self.recording.get_sampling_frequency()
        
        with pdf.PdfPages(pdf_path) as pdf_doc:
            units_per_page = 12 
            for i in range(0, len(unit_ids), units_per_page):
                batch = unit_ids[i : i + units_per_page]
                fig, axes = plt.subplots(3, 4, figsize=(12, 9))
                axes = axes.flatten()
                
                for ax, uid in zip(axes, batch):
                    wf = wf_ext.get_waveforms_one_unit(uid)
                    mean_wf = np.mean(wf, axis=0)
                    best_ch = np.argmin(np.min(mean_wf, axis=0))
                    
                    # [2] Create time vector in milliseconds
                    time_ms = np.arange(wf.shape[1]) / fs * 1000

                    # Random subset logic
                    n_spikes = wf.shape[0]
                    if n_spikes > 100:
                        indices = np.random.choice(n_spikes, 100, replace=False)
                        spikes_to_plot = wf[indices, :, best_ch]
                    else:
                        spikes_to_plot = wf[:, :, best_ch]
                    
                    # [3] Plot against time_ms (X-axis)
                    ax.plot(time_ms, spikes_to_plot.T, c='gray', lw=0.5, alpha=0.3)
                    ax.plot(time_ms, mean_wf[:, best_ch], c='red', lw=1.5)
                    
                    ax.set_title(f"Unit {uid} | Ch {best_ch}", fontsize=10)
                    
                    # [4] Add Scale Bar (1 ms x 50 uV)
                    try:
                        add_scalebar(ax, 
                                    matchx=False, matchy=False, 
                                    sizex=1.0, labelx='1 ms', 
                                    sizey=50, labely='50 µV', 
                                    loc='lower right', 
                                    hidex=True, hidey=True)
                    except Exception:
                        ax.spines['top'].set_visible(False)
                        ax.spines['right'].set_visible(False)

                for j in range(len(batch), len(axes)): 
                    axes[j].axis('off')
                
                pdf_doc.savefig(fig)
                plt.close(fig)

    def _run_burst_analysis(self, ids_list=None, plot_mode='separate', plot_debug=False, raster_sort='none', fixed_y = False):
        self.logger.info("Running Network Burst Analysis...")

        spike_times = {}

        # 1. Load Spike Times
        if self.sorting:
            fs = self.recording.get_sampling_frequency()
            if ids_list is None:
                ids_list = self.analyzer.unit_ids

            missing_unit_ids = []
            for uid in ids_list:
                try:
                    spike_times[uid] = self.sorting.get_unit_spike_train(uid) / fs
                except KeyError:
                    missing_unit_ids.append(uid)

            if missing_unit_ids:
                self.logger.warning(
                    "Skipping %d unit(s) not present in active sorting during burst analysis: %s",
                    len(missing_unit_ids),
                    missing_unit_ids[:20],
                )

            if not spike_times:
                self.logger.error(
                    "No valid units left for burst analysis after filtering missing unit IDs."
                )
                return

            np.save(self.output_dir / "spike_times.npy", spike_times)
        else:
            # No spike sorter available: load from the spike_times.npy file saved by _spike_detection_only
            spike_times_file = self.output_dir / "spike_times.npy"
            if spike_times_file.exists():
                spike_times = np.load(spike_times_file, allow_pickle=True).item()
                if not spike_times:
                    self.logger.error("No spike times found in saved spike_times.npy file.")
                    return
                if ids_list is not None:
                    spike_times = {uid: spike_times[uid] for uid in ids_list if uid in spike_times}
                if not spike_times:
                    self.logger.error("No spike times found for burst analysis after filtering by ids_list.")
                    return
            else:
                self.logger.error("No spike times found for burst analysis.")
                return

        if not spike_times:
            self.logger.warning("Spike times dictionary is empty. Skipping burst analysis.")
            return

        # ---------------------------------------------------------
        # 2. Main analysis block
        # ---------------------------------------------------------
        try:
            # A. Run network burst detector
            network_data = compute_network_bursts(
                SpikeTimes=spike_times,
            )

            if isinstance(network_data, dict) and "error" in network_data:
                self.logger.error(f"Burst detector returned error: {network_data['error']}")
                return

            # B. Save clean JSON
            network_data_clean = helper.recursive_clean(network_data)
            network_data_clean["n_units"] = len(spike_times)

            temp_file = self.output_dir / "network_results.tmp.json"
            final_file = self.output_dir / "network_results.json"

            with open(temp_file, "w") as f:
                json.dump(network_data_clean, f, indent=2)

            if temp_file.exists():
                os.replace(temp_file, final_file)
                self.logger.info(f"Successfully saved: {final_file}")

            # C. Sort units for raster
            sorted_units = self._sort_units_for_raster(spike_times, raster_sort)

            # D. Build figure
            ax_network_red = None

            if plot_mode == "separate":
                fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
                ax_raster, ax_network = axs

                helper.plot_clean_raster(
                    ax_raster,
                    spike_times,
                    sorted_units,
                    color="gray",
                    markersize=4,
                    markeredgewidth=0.5,
                    alpha=1.0
                )

                ax_network, ax_network_red = helper.plot_clean_network(
                    ax_network,
                    **network_data["plot_data"],
                    use_twinx=True
                )

            elif plot_mode == "merged":
                fig, ax_raster = plt.subplots(figsize=(12, 5))

                helper.plot_clean_raster(
                    ax_raster,
                    spike_times,
                    sorted_units,
                    color="gray",
                    markersize=4,
                    markeredgewidth=0.5,
                    alpha=1.0
                )

                ax_network = ax_raster.twinx()

                ax_network, ax_network_red = helper.plot_clean_network(
                    ax_network,
                    **network_data["plot_data"],
                    use_twinx=False
                )

                ax_raster.spines["right"].set_visible(False)
                ax_network.spines["right"].set_visible(True)

            else:
                self.logger.warning(f"Unknown plot mode: {plot_mode}")
                return

            # ---------------------------------------------------------
            # 3. Overlay burst hierarchy
            # ---------------------------------------------------------
            burstlet_events = network_data["burstlets"]["events"]
            network_burst_events = network_data["network_bursts"]["events"]
            superburst_events = network_data["superbursts"]["events"]

            helper.mark_burst_hierarchy(
                ax_raster=ax_raster,
                ax_network=ax_network,
                burstlets=burstlet_events,
                network_bursts=network_burst_events,
                superbursts=superburst_events,
                show_raster_spans=False,      # no raster shading for now
                show_burstlet_ticks=True,     # black top ticks
                show_network_ticks=True,      # blue top ticks
                show_superburst_bars=True,    # purple top bars
                min_superburst_duration_s=2.5
            )

            # ---------------------------------------------------------
            # 4. Legend
            # ---------------------------------------------------------
            from matplotlib.lines import Line2D

            hierarchy_handles = [
                Line2D([0], [0], color="black", lw=1.2, label="Burstlet ticks"),
                Line2D([0], [0], color="steelblue", lw=2.0, label="Network burst ticks"),
                Line2D([0], [0], color="mediumpurple", lw=2.2, label="Superbursts"),
                Line2D([0], [0], marker='o', color='red', lw=0, markersize=5, label="Network burst centers"),
            ]
            ax_raster.legend(handles=hierarchy_handles, loc="upper right", frameon=False, fontsize=8)

            # ---------------------------------------------------------
            # 5. Layout and save
            # ---------------------------------------------------------
            plt.tight_layout()
            if plot_mode == "separate":
                plt.subplots_adjust(hspace=0.05)

            full_svg = self.output_dir / "raster_burst_plot.svg"
            full_png = self.output_dir / "raster_burst_plot.png"
            zoom60_svg = self.output_dir / "raster_burst_plot_60s.svg"
            zoom30_svg = self.output_dir / "raster_burst_plot_30s.svg"

            plt.savefig(full_svg)

            # 60 s zoom
            ax_raster.set_xlim(0, 60)
            ax_network.set_xlim(0, 60)
            if ax_network_red is not None and ax_network_red is not ax_network:
                ax_network_red.set_xlim(0, 60)
            plt.savefig(zoom60_svg)

            # 30 s zoom
            ax_raster.set_xlim(0, 30)
            ax_network.set_xlim(0, 30)
            if ax_network_red is not None and ax_network_red is not ax_network:
                ax_network_red.set_xlim(0, 30)
            ax_network.set_xlabel("Time (s)")
            plt.savefig(zoom30_svg)

            plt.savefig(full_png, dpi=300)
            plt.close(fig)

            self.logger.info("Burst analysis plots saved successfully.")

            # Fixed Y Logic
            if fixed_y:
                summary_file = self.output_root / self.project_name / f"{self.project_name}_y_max_summary.json"
                if not summary_file.exists():
                    self.logger.error(f"No y-max summary found at {summary_file}. Run without --fixed-y first.")
                else:
                    with open(summary_file, 'r') as f:
                        summary = json.load(f)
                    # Flatten all values and find global max
                    all_maxima = [
                        v for date in summary.values()
                        for chip in date.values()
                        for v in chip.values()
                    ]
                    global_max = max(all_maxima)
                    self.logger.info(f"Applying fixed y-max: {global_max:.4f}")

                    # Replot with fixed y
                    fig2, axs2 = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
                    ax_raster2, ax_network2 = axs2
                    helper.plot_clean_raster(ax_raster2, spike_times, color='gray', markersize=4, markeredgewidth=0.5, alpha=1.0)
                    helper.plot_clean_network(ax_network2, **network_data["plot_data"])
                    ax_network2.set_ylim(0, global_max)
                    plt.tight_layout()
                    plt.subplots_adjust(hspace=0.05)
                    for start, end in sb_intervals:
                        ax_network2.axvspan(start, end, color='gray', alpha=0.3)
                    plt.savefig(self.output_dir / "fixed_y_raster_burst_plot.svg")
                    plt.savefig(self.output_dir / "fixed_y_raster_burst_plot.png", dpi=300)
                    ax_raster2.set_xlim(0, 60)
                    ax_network2.set_xlim(0, 60)
                    plt.savefig(self.output_dir / "fixed_y_raster_burst_plot_60s.svg")
                    ax_raster2.set_xlim(0, 30)
                    ax_network2.set_xlim(0, 30)
                    ax_network2.set_xlabel("Time (s)")
                    plt.savefig(self.output_dir / "fixed_y_raster_burst_plot_30s.svg")
                    plt.savefig(self.output_dir / "fixed_y_raster_burst_plot_30s.png", dpi=300)

                    plt.close(fig2)

        except Exception as e:
            self.logger.error(f"Burst analysis error: {e}")
            traceback.print_exc()
            raise e

    def _sort_units_for_raster(self, spike_times, raster_sort):
        """Returns ordered list of unit keys for raster y-axis."""
        if raster_sort == 'none':
            return None  # plot_clean_raster handles default ordering itself

        if raster_sort == 'firing_rate':
            return sorted(spike_times.keys(), key=lambda uid: len(spike_times[uid]))

        elif raster_sort == 'unit_id':
            return sorted(spike_times.keys())

        self.logger.warning(f"Unknown raster_sort: {raster_sort}. Falling back to none.")
        return None

    # --- Cleanup ---
    def cleanup(self):
        if self.cleanup_flag:
            self.logger.info("Cleaning up temp files...")
            shutil.rmtree(self.output_dir / "binary", ignore_errors=True)
            shutil.rmtree(self.output_dir / "sorter_output", ignore_errors=True)
        self.recording = None
        self.sorting = None
        self.analyzer = None
        gc.collect()


@dataclass(frozen=True)
class MEARunOptions:
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
    pipeline: MEAPipeline
    skipped: bool = False
    reanalyzed_bursts: bool = False


def _normalize_resume_from_stage(resume_from: str | None) -> str | None:
    if resume_from is None:
        return None

    token = str(resume_from).strip().lower().replace("-", "_")
    if not token:
        return None

    aliases = {
        "preprocess": "preprocessing",
        "preprocessing": "preprocessing",
        "sorting": "sorting",
        "sort": "sorting",
        "merge": "merge",
        "analyzer": "analyzer",
        "analyse": "analyzer",
        "analysis": "analyzer",
        "report": "reports",
        "reports": "reports",
    }
    normalized = aliases.get(token)
    if normalized is None:
        valid = ", ".join(["preprocessing", "sorting", "merge", "analyzer", "reports"])
        raise ValueError(f"Invalid resume_from stage '{resume_from}'. Valid stages: {valid}")
    return normalized


def _apply_resume_from_stage(pipeline: MEAPipeline, resume_from: str | None) -> None:
    stage_name = _normalize_resume_from_stage(resume_from)
    if stage_name is None:
        return

    # Resume semantics: set checkpoint to the stage immediately before target execution.
    resume_checkpoint_stage = {
        "preprocessing": ProcessingStage.NOT_STARTED,
        "sorting": ProcessingStage.PREPROCESSING_COMPLETE,
        "merge": ProcessingStage.SORTING_COMPLETE,
        "analyzer": ProcessingStage.MERGE_COMPLETE,
        "reports": ProcessingStage.ANALYZER_COMPLETE,
    }[stage_name]

    if stage_name in {"merge", "analyzer", "reports"}:
        pipeline.force_rerun_analyzer = True

    pipeline._save_checkpoint(
        resume_checkpoint_stage,
        failed_stage=None,
        error=None,
        resume_from=stage_name,
        resume_forced_rerun_analyzer=bool(pipeline.force_rerun_analyzer),
    )
    pipeline.logger.info(
        "Resume-from requested: %s (checkpoint set to %s)",
        stage_name,
        resume_checkpoint_stage.name,
    )


def run_mea_pipeline(options: MEARunOptions) -> MEARunResult:
    um_kwargs = _default_um_kwargs()
    if isinstance(options.um_kwargs, dict):
        um_kwargs.update(options.um_kwargs)

    am_kwargs = _default_am_kwargs()
    if isinstance(options.am_kwargs, dict):
        am_kwargs.update(options.am_kwargs)

    option_kwargs = _default_option_kwargs()
    if isinstance(options.option_kwargs, dict):
        option_kwargs.update(options.option_kwargs)

    auto_merge_presets = am_kwargs.get("presets")
    auto_merge_steps_params = am_kwargs.get("steps_params")

    if bool(am_kwargs.get("enabled")) and (auto_merge_presets is None or auto_merge_steps_params is None):
        try:
            diffs = [
                float(x.strip())
                for x in str(am_kwargs.get("template_diff_thresh", "0.05,0.15,0.25")).split(",")
                if x.strip()
            ]
            auto_merge_presets = ["x_contaminations"] * len(diffs)
            auto_merge_steps_params = [
                {"template_similarity": {"template_diff_thresh": float(t)}}
                for t in diffs
            ]
        except Exception:
            auto_merge_presets = None
            auto_merge_steps_params = None

    am_kwargs["presets"] = auto_merge_presets
    am_kwargs["steps_params"] = auto_merge_steps_params

    pipeline = MEAPipeline(
        file_path=options.file_path,
        stream_id=options.stream_id,
        recording_num=options.recording_num,
        output_root=options.output_root,
        checkpoint_root=options.checkpoint_root,
        sorter=options.sorter,
        docker_image=options.docker_image,
        verbose=bool(options.verbose),
        cleanup=bool(options.cleanup),
        force_restart=bool(options.force_restart),
        n_jobs=options.n_jobs,
        chunk_duration=options.chunk_duration,
        sorter_kwargs=options.sorter_kwargs,
        um_kwargs=um_kwargs,
        am_kwargs=am_kwargs,
        option_kwargs=option_kwargs,
    )

    _apply_resume_from_stage(pipeline, options.resume_from)

    if bool(options.reanalyze_bursts):
        pipeline._run_burst_analysis(
            plot_mode=options.plot_mode,
            plot_debug=bool(options.plot_debug),
            raster_sort=options.raster_sort,
            fixed_y=bool(options.fixed_y),
        )
        current_stage = ProcessingStage(pipeline.state['stage'])
        pipeline._save_checkpoint(
            current_stage,
            note="Burst Re-analysis Performed",
            last_updated=str(datetime.now()),
        )
        return MEARunResult(pipeline=pipeline, skipped=False, reanalyzed_bursts=True)

    if pipeline.should_skip():
        return MEARunResult(pipeline=pipeline, skipped=True, reanalyzed_bursts=False)

    pipeline.run_preprocessing()

    if not bool(options.skip_spikesorting):
        pipeline.run_sorting()
        pipeline.run_optional_merge_phase()

        if bool(options.run_analyzer):
            pipeline.run_analyzer()

        if bool(options.run_reports):
            if (not bool(options.run_analyzer)) and pipeline.analyzer is None:
                raise RuntimeError(
                    "Requested report generation but analyzer was not run and no existing analyzer was loaded. "
                    "Set run_analyzer=True or run once to populate analyzer_output."
                )
            pipeline.generate_reports(
                options.thresholds,
                bool(options.no_curation),
                bool(options.export_to_phy),
                plot_mode=options.plot_mode,
                plot_debug=bool(options.plot_debug),
                raster_sort=options.raster_sort,
                fixed_y=bool(options.fixed_y),
            )
    else:
        ids = pipeline._spike_detection_only()
        pipeline._run_burst_analysis(
            ids,
            plot_mode=options.plot_mode,
            plot_debug=bool(options.plot_debug),
            raster_sort=options.raster_sort,
            fixed_y=bool(options.fixed_y),
        )

    if bool(options.cleanup):
        pipeline.cleanup()

    return MEARunResult(pipeline=pipeline, skipped=False, reanalyzed_bursts=False)

# --- CLI Entry Point ---
def main():
    parser = argparse.ArgumentParser(
        description="MEA Analysis Routine — processes a single well from an MEA recording",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # --- Positional ---
    parser.add_argument("file_path",
        help="Path to .h5, .nwb, or .raw MEA recording file")

    # --- Input / Output ---
    io_group = parser.add_argument_group("input/output")
    io_group.add_argument("--config", type=str, default=None,
        help="Path to config JSON file (CLI flags always override config)")
    io_group.add_argument("--well", required=True,
        help="Well ID to process (e.g. well000)")
    io_group.add_argument("--rec", type=str, default=None,
        help="Recording name inside HDF5 file (default: rec0000)")
    io_group.add_argument("--output-dir", type=str, default=None,
        help="Output directory for results")
    io_group.add_argument("--checkpoint-dir", type=str, default=None,
        help="Checkpoint directory (default: <output-dir>/checkpoints)")
    io_group.add_argument("--output-subdir-after-well", type=str, default=None,
        help="Optional single subdirectory appended under the resolved well output directory")
    io_group.add_argument("--export-to-phy", action="store_true",
        help="Export results to Phy format")
    io_group.add_argument("--clean-up", action="store_true",
        help="Remove intermediate files after processing")

    # --- Sorting ---
    sort_group = parser.add_argument_group("sorting")
    sort_group.add_argument("--sorter", type=str, default=None,
        help="Spike sorter to use (default: kilosort4)")
    sort_group.add_argument("--docker", type=str, default=None,
        help="Docker image name for containerized sorting")
    sort_group.add_argument("--skip-spikesorting", action="store_true",
        help="Run spike detection only, skip full sorting")

    # --- Plotting ---
    plot_group = parser.add_argument_group("plotting")
    plot_group.add_argument("--plot-mode", choices=["separate", "merged"], default=None,
        help="Plot raster and network on separate axes or merged twin-axis\n(default: separate)")
    plot_group.add_argument("--raster-sort", choices=["none", "firing_rate", "location_y", "unit_id"], default=None,
        help="How to sort units on raster y-axis (default: none)")
    plot_group.add_argument("--plot-debug", action="store_true",
        help="Overlay burst and superburst intervals on raster plot")
    plot_group.add_argument("--fixed-y", action="store_true",
        help="Use fixed y-axis limits for raster plots — run once without it first to generate summary")
    # --- Curation ---
    cur_group = parser.add_argument_group("curation")
    cur_group.add_argument("--no-curation", action="store_true",
        help="Skip automatic unit curation")
    cur_group.add_argument("--params", type=str, default=None,
        help="JSON string or file path with quality thresholds")

    # --- Run Control ---
    ctrl_group = parser.add_argument_group("run control")
    ctrl_group.add_argument("--force-restart", action="store_true",
        help="Ignore checkpoint and restart from scratch")
    ctrl_group.add_argument("--resume-from", "--resume_from", dest="resume_from", type=str, default=None,
        choices=["preprocessing", "sorting", "merge", "analyzer", "reports"],
        help="Resume by rewinding checkpoint to just before this stage and rerunning from there")
    ctrl_group.add_argument("--reanalyze-bursts", action="store_true",
        help="Re-run burst analysis on existing spike times only")
    ctrl_group.add_argument("--debug", action="store_true",
        help="Enable verbose logging")

    # Optional post-spikesort step(s)
    parser.add_argument(
        "--unitmatch-merge-units",
        action='store_true',
        help="Optional (default off): run UnitMatch integration as an alternative to auto_merge_units.",
    )
    parser.add_argument(
        "--unitmatch-dry-run",
        action='store_true',
        help="When UnitMatch is enabled, produce UnitMatch reports without applying merges.",
    )
    parser.add_argument(
        "--unitmatch-scored-dry-run",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="When UnitMatch dry-run is enabled, attempt backend scoring (default: enabled).",
    )
    parser.add_argument(
        "--unitmatch-output-subdir-name",
        type=str,
        default=None,
        help="UnitMatch artifact subdirectory under output_dir (default: unitmatch_outputs).",
    )
    parser.add_argument(
        "--unitmatch-throughput-subdir-name",
        type=str,
        default=None,
        help="UnitMatch throughput subdirectory under output_dir (default: unitmatch_throughput).",
    )
    parser.add_argument(
        "--unitmatch-max-candidate-pairs",
        type=int,
        default=None,
        help="Maximum UnitMatch candidate pairs (-1 unlimited, 0 none, default: 20000).",
    )
    parser.add_argument(
        "--unitmatch-oversplit-min-probability",
        type=float,
        default=None,
        help="Minimum UnitMatch probability for oversplit suggestions (default: 0.80).",
    )
    parser.add_argument(
        "--unitmatch-oversplit-max-suggestions",
        type=int,
        default=None,
        help="Maximum oversplit suggestions (-1 unlimited, 0 none, default: 2000).",
    )
    parser.add_argument(
        "--unitmatch-apply-merges",
        action='store_true',
        help="Apply conflict-free top UnitMatch suggestions to mutate sorting.",
    )
    parser.add_argument(
        "--unitmatch-recursive",
        action='store_true',
        help="Recursively run UnitMatch merge iterations until convergence or cap.",
    )
    parser.add_argument(
        "--unitmatch-max-iterations",
        type=int,
        default=None,
        help="Maximum recursive UnitMatch iterations (-1 uncapped, default: 5).",
    )
    parser.add_argument(
        "--unitmatch-max-spikes-per-unit",
        type=int,
        default=None,
        help="Max spikes per unit for UnitMatch raw-waveform generation (-1 uncapped, default: 100).",
    )
    parser.add_argument(
        "--unitmatch-keep-all-iterations",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Keep all unitmatch_throughput iteration folders (default: enabled).",
    )
    parser.add_argument(
        "--unitmatch-generate-reports",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Generate static UnitMatch report pack from existing artifacts (default: enabled).",
    )
    parser.add_argument(
        "--unitmatch-report-subdir-name",
        type=str,
        default="unitmatch_reports",
        help="UnitMatch report output subdirectory under output_dir (default: unitmatch_reports).",
    )
    parser.add_argument(
        "--unitmatch-report-max-heatmap-units",
        type=int,
        default=200,
        help="Maximum units rendered in UnitMatch similarity heatmap (default: 200).",
    )
    parser.add_argument(
        "--auto-merge-units",
        action='store_true',
        help="Optional (default off): run SpikeInterface auto_merge_units during analyzer stage.",
    )
    parser.add_argument(
        "--auto-merge-template-diff-thresh",
        default="0.05,0.15,0.25",
        help="Comma-separated template_diff_thresh values for auto-merge (used with preset x_contaminations).",
    )
    parser.add_argument(
        "--rerun-analyzer",
        action='store_true',
        help="Recompute analyzer_output even if checkpoint says complete (does not rerun spikesorting).",
    )

    args = parser.parse_args()
    config = load_config(args.config)
    resolved = resolve_args(args, config)
    fixed_y = resolved["fixed_y"]

   
    plot_mode = resolved["plot_mode"]
    plot_debug = resolved["plot_debug"]
    raster_sort = resolved["raster_sort"]
    sorter = resolved["sorter"]
    thresholds = resolved["quality_thresholds"]
    rec = args.rec or "rec0000"

    try:
        if args.reanalyze_bursts:
            print("Re-analyzing bursts only on existing spike times...")

        result = run_mea_pipeline(
            MEARunOptions(
                file_path=args.file_path,
                stream_id=args.well,
                recording_num=rec,
                output_root=resolved["output_dir"],
                checkpoint_root=resolved["checkpoint_dir"],
                sorter=sorter,
                docker_image=resolved["docker_image"],
                verbose=bool(args.debug),
                cleanup=bool(resolved["clean_up"]),
                force_restart=bool(args.force_restart),
                resume_from=args.resume_from,
                um_kwargs={
                    "merge_units": bool(args.unitmatch_merge_units),
                    "dry_run": bool(args.unitmatch_dry_run),
                    "scored_dry_run": bool(resolved["unitmatch_scored_dry_run"]),
                    "output_subdir_name": str(resolved["unitmatch_output_subdir_name"]),
                    "throughput_subdir_name": str(resolved["unitmatch_throughput_subdir_name"]),
                    "max_candidate_pairs": int(resolved["unitmatch_max_candidate_pairs"]),
                    "oversplit_min_probability": float(resolved["unitmatch_oversplit_min_probability"]),
                    "oversplit_max_suggestions": int(resolved["unitmatch_oversplit_max_suggestions"]),
                    "apply_merges": bool(resolved["unitmatch_apply_merges"]),
                    "recursive": bool(resolved["unitmatch_recursive"]),
                    "max_iterations": int(resolved["unitmatch_max_iterations"]),
                    "max_spikes_per_unit": int(resolved["unitmatch_max_spikes_per_unit"]),
                    "keep_all_iterations": bool(resolved["unitmatch_keep_all_iterations"]),
                    "generate_reports": bool(args.unitmatch_generate_reports),
                    "report_subdir_name": str(args.unitmatch_report_subdir_name),
                    "report_max_heatmap_units": int(args.unitmatch_report_max_heatmap_units),
                },
                am_kwargs={
                    "enabled": bool(args.auto_merge_units),
                    "template_diff_thresh": str(args.auto_merge_template_diff_thresh),
                },
                option_kwargs={
                    "force_rerun_analyzer": bool(args.rerun_analyzer),
                    "output_subdir_after_well": resolved.get("output_subdir_after_well"),
                },
                reanalyze_bursts=bool(args.reanalyze_bursts),
                skip_spikesorting=bool(args.skip_spikesorting),
                run_analyzer=True,
                run_reports=True,
                thresholds=thresholds,
                no_curation=bool(resolved["no_curation"]),
                export_to_phy=bool(resolved["export_to_phy"]),
                plot_mode=plot_mode,
                plot_debug=bool(plot_debug),
                raster_sort=raster_sort,
                fixed_y=bool(fixed_y),
                auto_merge_template_diff_thresh=str(args.auto_merge_template_diff_thresh),
            )
        )

        if result.reanalyzed_bursts:
            print("Burst Re-analysis Complete.")
            sys.exit(0)

        if result.skipped:
            sys.exit(0)

        print(f"Processing Complete for {args.well}")

    except Exception as e:
        print(f"CRITICAL FAILURE in {args.well}: {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()