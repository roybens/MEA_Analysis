# ==========================================================
# _mea_runner.py
# Author: Mandar Patil
# Contributors: Yuxin Ren, Shruti Shah, Adam Weiner
# LLM Assisted Edits: Yes ChatGPT-4, Claude sonnet 4.6, ChatGPT-5.3-Codex
#
# Pipeline runner: MEARunOptions, MEARunResult, and run_mea_pipeline.
# MEAPipeline class lives in _mea_pipeline.py.
# CLI entry point (main) lives in mea_analysis_routine.py.
# ==========================================================

import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

# If this file is loaded directly (not imported as a package module),
# allow local imports by adding the script folder and repo root to sys.path.
if not __package__:
    script_dir = Path(__file__).resolve().parent
    if str(script_dir) not in sys.path:
        sys.path.append(str(script_dir))

    root_dir = script_dir.parent.parent
    if str(root_dir) not in sys.path:
        sys.path.append(str(root_dir))

# Import MEAPipeline and supporting symbols from the class module.
try:
    if __package__:
        from ._mea_pipeline import (
            MEAPipeline,
            ProcessingStage,
            _default_um_kwargs,
            _default_am_kwargs,
            _default_option_kwargs,
        )
    else:
        from _mea_pipeline import (
            MEAPipeline,
            ProcessingStage,
            _default_um_kwargs,
            _default_am_kwargs,
            _default_option_kwargs,
        )
except ImportError:
    try:
        from MEA_Analysis.IPNAnalysis._mea_pipeline import (
            MEAPipeline,
            ProcessingStage,
            _default_um_kwargs,
            _default_am_kwargs,
            _default_option_kwargs,
        )
    except ImportError as e:
        raise ImportError(
            "Could not import MEA_Analysis.IPNAnalysis._mea_pipeline. "
            "If you are running this file directly, prefer: "
            "`python -m MEA_Analysis.IPNAnalysis.mea_analysis_routine ...`"
        ) from e


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
