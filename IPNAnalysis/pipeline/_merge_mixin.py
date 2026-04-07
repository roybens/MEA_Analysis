# ==========================================================
# pipeline/_merge_mixin.py
# MergeMixin — Phase 2.5: UnitMatch and SpikeInterface auto-merge.
# ==========================================================

import json
import shutil

import spikeinterface.full as si
import spikeinterface.curation as sic

from .stages import ProcessingStage

try:
    from ..UnitMatch.runner import run_unitmatch_merge_with_recursion, UnitMatchConfig
    from ..UnitMatch.reporting import UnitMatchReportConfig, generate_unitmatch_static_report_pack
except ImportError:
    try:
        from MEA_Analysis.IPNAnalysis.UnitMatch.runner import (
            run_unitmatch_merge_with_recursion,
            UnitMatchConfig,
        )
        from MEA_Analysis.IPNAnalysis.UnitMatch.reporting import (
            UnitMatchReportConfig,
            generate_unitmatch_static_report_pack,
        )
    except ImportError:
        from UnitMatch.runner import run_unitmatch_merge_with_recursion, UnitMatchConfig
        from UnitMatch.reporting import UnitMatchReportConfig, generate_unitmatch_static_report_pack


class MergeMixin:
    """Handles the optional post-sorting unit-merge phase."""

    def run_optional_merge_phase(self):
        self.logger.info("--- [Phase 2.5] Merge (Optional) ---")

        try:
            current_stage = ProcessingStage(
                int(self.state.get('stage', ProcessingStage.SORTING_COMPLETE.value))
            )
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
            self._run_unitmatch_merge(current_stage)
            return

        if self.auto_merge_units:
            self._run_auto_merge()
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

    # ------------------------------------------------------------------
    # UnitMatch path
    # ------------------------------------------------------------------

    def _run_unitmatch_merge(self, current_stage):
        try:
            if self.auto_merge_units:
                self.logger.warning(
                    "Both UnitMatch and auto_merge_units requested; "
                    "using UnitMatch path and skipping auto_merge_units"
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
            persisted_final_merged_sorting = bool(
                merge_application.get("persisted_final_merged_sorting", False)
            )

            summary_path_s = str(artifacts.get("summary_json") or "")
            if summary_path_s:
                try:
                    from pathlib import Path as _Path
                    summary_path = _Path(summary_path_s)
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
                        str(final_merged_sorting_folder)
                        if final_merged_sorting_folder is not None
                        else None
                    ),
                    "persisted_final_merged_sorting": persisted_final_merged_sorting,
                    "recursive": bool(self.unitmatch_recursive),
                    "n_iterations": int(
                        um_summary.get("iteration_convergence_report", {}).get(
                            "n_iterations_executed", 1
                        )
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
            self.logger.warning(
                "UnitMatch merge path failed; continuing without UnitMatch merge (%s)", e
            )
            self._save_checkpoint(
                ProcessingStage.MERGE_COMPLETE,
                merge_phase={
                    "status": "error",
                    "mode": "unitmatch",
                    "applied": False,
                    "error": str(e),
                },
            )

    # ------------------------------------------------------------------
    # SpikeInterface auto_merge_units path
    # ------------------------------------------------------------------

    def _run_auto_merge(self):
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

            merged_sorting = merged.sorting if hasattr(merged, "sorting") else merged

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
