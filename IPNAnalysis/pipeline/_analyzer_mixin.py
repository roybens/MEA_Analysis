# ==========================================================
# pipeline/_analyzer_mixin.py
# AnalyzerMixin — Phase 3: SpikeInterface sorting analyzer.
# ==========================================================

import shutil
import traceback
from datetime import datetime

import spikeinterface.full as si

from .stages import ProcessingStage


class AnalyzerMixin:
    """Computes the SpikeInterface SortingAnalyzer and its extensions."""

    def run_analyzer(self):
        analyzer_folder = self.output_dir / "analyzer_output"

        # Resume check — reload existing analyzer if already computed.
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
                self.sorting,
                self.recording,
                method="radius",
                radius_um=50,
                peak_sign='neg',
            )

            self.analyzer = si.create_sorting_analyzer(
                self.sorting,
                self.recording,
                format="binary_folder",
                folder=analyzer_folder,
                sparsity=sparsity,
                return_in_uV=True,
            )

            ext_list = [
                "random_spikes",
                "spike_amplitudes",
                "waveforms",
                "templates",
                "noise_levels",
                "quality_metrics",
                "template_metrics",
                "unit_locations",
            ]
            ext_params = {
                "waveforms": {"ms_before": 1.0, "ms_after": 2.0},
                "unit_locations": {"method": "monopolar_triangulation"},
            }

            compute_kwargs = {'verbose': self.verbose}
            if self.n_jobs is not None:
                compute_kwargs['n_jobs'] = int(self.n_jobs)
            self.analyzer.compute(ext_list, extension_params=ext_params, **compute_kwargs)

            # Always use analyzer sorting as the canonical source of truth.
            self.sorting = self.analyzer.sorting
            self._save_checkpoint(ProcessingStage.ANALYZER_COMPLETE, failed_stage=None, error=None)
        except Exception as e:
            err = {
                "failed_stage": ProcessingStage.ANALYZER.name,
                "exception": type(e).__name__,
                "message": str(e),
                "traceback": traceback.format_exc(),
                "time": str(datetime.now()),
            }
            self.logger.error(err["traceback"])
            try:
                fallback_stage = ProcessingStage(
                    int(self.state.get('stage', ProcessingStage.MERGE_COMPLETE.value))
                )
            except Exception:
                fallback_stage = ProcessingStage.MERGE_COMPLETE
            self._save_checkpoint(fallback_stage, error=err)
            raise
