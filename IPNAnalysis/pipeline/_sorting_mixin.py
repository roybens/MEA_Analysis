# ==========================================================
# pipeline/_sorting_mixin.py
# SortingMixin — Phase 2: spike sorting and spike-detection-only path.
# ==========================================================

import traceback
from datetime import datetime
from timeit import default_timer as timer

import numpy as np
import spikeinterface.full as si
from spikeinterface.sortingcomponents.peak_detection import detect_peaks

from .stages import ProcessingStage


class SortingMixin:
    """Encapsulates the spike-sorting and spike-detection phases."""

    # ------------------------------------------------------------------
    # Phase 2: Spike Sorting
    # ------------------------------------------------------------------

    def run_sorting(self):
        sorter_folder = self.output_dir / "sorter_output"

        # Resume check — reload from existing sorter output.
        if self.state['stage'] >= ProcessingStage.SORTING_COMPLETE.value:
            self.logger.info("Resuming: Loading existing sorting.")
            try:
                self.sorting = si.read_sorter_folder(sorter_folder)
            except Exception:
                self.sorting = si.read_kilosort(sorter_folder)
            return

        self._save_checkpoint(ProcessingStage.SORTING)
        self.logger.info(f"--- [Phase 2] Spike Sorting ({self.sorter}) ---")

        import torch
        if torch.cuda.is_available():
            total_vram = torch.cuda.get_device_properties(0).total_memory / (1024 ** 3)
            self.logger.info(
                f"GPU Detected: {torch.cuda.get_device_name(0)} with {total_vram:.2f} GB VRAM"
            )
        else:
            self.logger.warning(
                "No GPU detected! Kilosort4 will likely fail or run extremely slowly on CPU."
            )
            total_vram = 0

        # Parameters for 2D static cultures (Kilosort4).
        ks_params_high_vram = {
            'batch_size': int(self.recording.get_sampling_frequency()) * 2,
            'clear_cache': True,
            'invert_sign': True,
            'cluster_downsampling': 20,
            'max_cluster_subset': None,
            'nblocks': 0,
            'dmin': 17,
            'do_correction': False,
        }
        ks_params_low_vram = {
            'batch_size': int(self.recording.get_sampling_frequency() * 0.5),
            'clear_cache': True,
            'invert_sign': True,
            'cluster_downsampling': 30,
            'max_cluster_subset': 50000,
            'nblocks': 0,
            'do_correction': False,
        }

        ks_params = ks_params_high_vram if total_vram >= 14 else ks_params_low_vram

        # Optional override from caller (e.g. debug harness / batch tuning).
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
                remove_existing_folder=True,
                verbose=self.verbose,
                docker_image=self.docker_image,
                **ks_params,
            )
            self.logger.info(f"Sorting finished in {timer() - start:.2f}s")

            # Remove spikes that extend beyond the recording length to prevent
            # IndexError during analyzer computation.
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
                "time": str(datetime.now()),
            }
            self.logger.error(err["traceback"])
            self._save_checkpoint(ProcessingStage.PREPROCESSING_COMPLETE, error=err)
            raise

    # ------------------------------------------------------------------
    # Phase 2-Alt: Spike Detection Only (no sorting)
    # ------------------------------------------------------------------

    def _spike_detection_only(self):
        """Detects spikes (threshold crossings) without running a full sorter."""
        self.logger.info("--- [Phase 2-Alt] Spike Detection (No Sorting) ---")

        job_kwargs = {
            'n_jobs': (int(self.n_jobs) if self.n_jobs is not None else 16),
            'chunk_duration': (str(self.chunk_duration) if self.chunk_duration is not None else '1s'),
            'progress_bar': self.verbose,
        }

        peaks = detect_peaks(
            self.recording,
            method='by_channel',
            detect_threshold=5,
            peak_sign='neg',
            exclude_sweep_ms=0.1,
            **job_kwargs,
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
