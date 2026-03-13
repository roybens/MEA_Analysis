# ==========================================================
# pipeline/runner.py
# High-level pipeline runner: run_mea_pipeline and stage helpers.
# ==========================================================

from datetime import datetime

from .core import MEAPipeline
from .stages import ProcessingStage
from .defaults import _default_um_kwargs, _default_am_kwargs, _default_option_kwargs
from .options import MEARunOptions, MEARunResult


# ------------------------------------------------------------------
# Stage-name helpers
# ------------------------------------------------------------------

def _normalize_resume_from_stage(resume_from: str | None) -> str | None:
    """Normalise user-supplied stage name to a canonical token."""
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
        raise ValueError(
            f"Invalid resume_from stage '{resume_from}'. Valid stages: {valid}"
        )
    return normalized


def _apply_resume_from_stage(pipeline: MEAPipeline, resume_from: str | None) -> None:
    """Rewind the pipeline checkpoint so execution restarts at *resume_from*."""
    stage_name = _normalize_resume_from_stage(resume_from)
    if stage_name is None:
        return

    # Set checkpoint to the stage immediately *before* the target stage.
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


# ------------------------------------------------------------------
# Main pipeline runner
# ------------------------------------------------------------------

def run_mea_pipeline(options: MEARunOptions) -> MEARunResult:
    """Run the full MEA analysis pipeline for a single well.

    Parameters
    ----------
    options:
        Immutable configuration bundle (see :class:`MEARunOptions`).

    Returns
    -------
    MEARunResult
        Contains the pipeline instance and run-status flags.
    """
    um_kwargs = _default_um_kwargs()
    if isinstance(options.um_kwargs, dict):
        um_kwargs.update(options.um_kwargs)

    am_kwargs = _default_am_kwargs()
    if isinstance(options.am_kwargs, dict):
        am_kwargs.update(options.am_kwargs)

    option_kwargs = _default_option_kwargs()
    if isinstance(options.option_kwargs, dict):
        option_kwargs.update(options.option_kwargs)

    # Pre-resolve auto-merge presets if enabled but not yet fully specified.
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

    # Burst re-analysis shortcut (skips the full pipeline).
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
                    "Requested report generation but analyzer was not run and no existing "
                    "analyzer was loaded.  Set run_analyzer=True or run once to populate "
                    "analyzer_output."
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
