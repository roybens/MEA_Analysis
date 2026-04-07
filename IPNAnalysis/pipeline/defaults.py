# ==========================================================
# pipeline/defaults.py
# Factory functions for default keyword-argument dictionaries.
# ==========================================================

from typing import Any


def _default_um_kwargs() -> dict[str, Any]:
    """Default settings for the UnitMatch merge phase."""
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
    """Default settings for the SpikeInterface auto-merge phase."""
    return {
        "enabled": False,
        "presets": None,
        "steps_params": None,
        "template_diff_thresh": "0.05,0.15,0.25",
    }


def _default_option_kwargs() -> dict[str, Any]:
    """Default settings for miscellaneous pipeline options."""
    return {
        "force_rerun_analyzer": False,
        "preprocessed_recording": None,
        "skip_preprocessing": False,
        "cuda_visible_devices": None,
        "output_subdir_after_well": None,
    }
