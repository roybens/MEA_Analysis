# ==========================================================
# config_loader.py
# Shared config loader for run_pipeline_driver.py and mea_analysis_routine.py
# Priority: CLI arg -> config file -> hardcoded default
# ==========================================================

import os
import json
from pathlib import Path

# ----------------------------------------------------------
# Hardcoded Defaults
# ----------------------------------------------------------
DEFAULTS = {
    "io": {
        "output_dir":       None,
        "checkpoint_dir":   None,
        "export_to_phy":    False,
        "clean_up":         False,
    },
    "sorting": {
        "sorter":           "kilosort4",
        "docker_image":     None,
    },
    "filtering": {
        "reference_file":   None,
        "assay_types":      ["network today", "network today/best"],
    },
    "plotting": {
        "plot_mode":        "separate",
        "raster_sort":      "none",
        "plot_debug":       False,
        "fixed_y":       False,
    },
    "curation": {
        "no_curation":      False,
        "quality_thresholds": {
            "presence_ratio":       0.75,
            "rp_contamination":     0.15,
            "firing_rate":          0.05,
            "amplitude_median":     -20,
            "amplitude_cv_median":  0.5,
        }
    }
}


# ----------------------------------------------------------
# Load config file
# ----------------------------------------------------------
def load_config(config_path=None):
    if config_path is None:
        return {}
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    with open(config_path) as f:
        return json.load(f)


# ----------------------------------------------------------
# Resolve a single value: CLI -> config -> default
# ----------------------------------------------------------
def _resolve(cli_val, cfg_val, default):
    if cli_val is not None:
        return cli_val
    if cfg_val is not None:
        return cfg_val
    return default


def _cfg(config, section, key):
    return config.get(section, {}).get(key, None)


# ----------------------------------------------------------
# Thresholds need deep merge: defaults -> config -> CLI --params
# ----------------------------------------------------------
def _resolve_thresholds(args, config):
    thresholds = DEFAULTS["curation"]["quality_thresholds"].copy()
    thresholds.update(config.get("curation", {}).get("quality_thresholds", {}))
    cli_params = getattr(args, "params", None)
    if cli_params:
        if os.path.exists(cli_params):
            with open(cli_params) as f:
                thresholds.update(json.load(f))
        else:
            thresholds.update(json.loads(cli_params))
    return thresholds


# ----------------------------------------------------------
# Bool store_true flags: only True if explicitly passed
# False from argparse means "not passed" so return None
# ----------------------------------------------------------
def _bool(args, attr):
    return True if getattr(args, attr, None) is True else None


# ----------------------------------------------------------
# Main resolver - works for both scripts
# ----------------------------------------------------------
def resolve_args(args, config):
    return {
        # io
        "output_dir":       _resolve(getattr(args, "output_dir", None),     _cfg(config, "io", "output_dir"),            DEFAULTS["io"]["output_dir"]),
        "checkpoint_dir":   _resolve(getattr(args, "checkpoint_dir", None),  _cfg(config, "io", "checkpoint_dir"),        DEFAULTS["io"]["checkpoint_dir"]),
        "export_to_phy":    _resolve(_bool(args, "export_to_phy"),           _cfg(config, "io", "export_to_phy"),         DEFAULTS["io"]["export_to_phy"]),
        "clean_up":         _resolve(_bool(args, "clean_up"),                _cfg(config, "io", "clean_up"),              DEFAULTS["io"]["clean_up"]),
        # sorting
        "sorter":           _resolve(getattr(args, "sorter", None),          _cfg(config, "sorting", "sorter"),           DEFAULTS["sorting"]["sorter"]),
        "docker_image":     _resolve(getattr(args, "docker", None),          _cfg(config, "sorting", "docker_image"),     DEFAULTS["sorting"]["docker_image"]),
        # filtering (driver only - ignored silently in routine)
        "reference_file":   _resolve(getattr(args, "reference", None),       _cfg(config, "filtering", "reference_file"), DEFAULTS["filtering"]["reference_file"]),
        "assay_types":      _resolve(getattr(args, "type", None),            _cfg(config, "filtering", "assay_types"),    DEFAULTS["filtering"]["assay_types"]),
        # plotting
        "plot_mode":        _resolve(getattr(args, "plot_mode", None),       _cfg(config, "plotting", "plot_mode"),       DEFAULTS["plotting"]["plot_mode"]),
        "raster_sort":      _resolve(getattr(args, "raster_sort", None),     _cfg(config, "plotting", "raster_sort"),     DEFAULTS["plotting"]["raster_sort"]),
        "plot_debug":       _resolve(_bool(args, "plot_debug"),              _cfg(config, "plotting", "plot_debug"),      DEFAULTS["plotting"]["plot_debug"]),
        "fixed_y":        _resolve(_bool(args, "fixed_y"),               _cfg(config, "plotting", "fixed_y"),         DEFAULTS["plotting"].get("fixed_y", False)),
        # curation
        "no_curation":      _resolve(_bool(args, "no_curation"),             _cfg(config, "curation", "no_curation"),     DEFAULTS["curation"]["no_curation"]),
        "quality_thresholds": _resolve_thresholds(args, config),
    }


# ----------------------------------------------------------
# Build extra_arg_string for driver subprocess
# Receives resolved dict + raw args for CLI-only pass-throughs
# ----------------------------------------------------------
def build_extra_args(resolved, cli_args):
    extra = []

    # Pass config downstream so subprocess inherits it
    if getattr(cli_args, "config", None):
        extra.append(f"--config '{cli_args.config}'")

    # io
    if resolved["output_dir"]:      extra.append(f"--output-dir '{resolved['output_dir']}'")
    if resolved["checkpoint_dir"]:  extra.append(f"--checkpoint-dir '{resolved['checkpoint_dir']}'")
    if resolved["export_to_phy"]:   extra.append("--export-to-phy")
    if resolved["clean_up"]:        extra.append("--clean-up")

    # sorting
    if resolved["sorter"]:          extra.append(f"--sorter {resolved['sorter']}")
    if resolved["docker_image"]:    extra.append(f"--docker {resolved['docker_image']}")

    # plotting
    if resolved["plot_mode"]:       extra.append(f"--plot-mode {resolved['plot_mode']}")
    if resolved["raster_sort"]:     extra.append(f"--raster-sort {resolved['raster_sort']}")
    if resolved["plot_debug"]:      extra.append("--plot-debug")

    # curation
    if resolved["no_curation"]:     extra.append("--no-curation")
    if resolved["fixed_y"]:        extra.append("--fixed-y")

    # CLI-only flags - passed through directly, never in config
    if getattr(cli_args, "force_restart", False):     extra.append("--force-restart")
    if getattr(cli_args, "debug", False):             extra.append("--debug")
    if getattr(cli_args, "skip_spikesorting", False): extra.append("--skip-spikesorting")
    if getattr(cli_args, "reanalyze_bursts", False):  extra.append("--reanalyze-bursts")

    return " ".join(extra)


# ----------------------------------------------------------
# Template generator - run directly to get a starter config
# python config_loader.py [optional_output_path]
# ----------------------------------------------------------
if __name__ == "__main__":
    import sys
    output = sys.argv[1] if len(sys.argv) > 1 else "mea_config_template.json"
    with open(output, "w") as f:
        json.dump(DEFAULTS, f, indent=2)
    print(f"Template written to: {output}")