from __future__ import annotations

from pathlib import Path
from typing import Any


_MULTISEG_OUTPUT_NAMES = {
    "common_electrodes",
    "maxwell_epochs",
    "concat_epochs",
    "final_merged_sorting",
    "concat_waveforms",
    "segment_waveforms",
    "wf_rejection_log",
    "curated_sorting",
}


def _is_multiseg_enabled(multiseg_mode: bool | None) -> bool:
    return bool(multiseg_mode)


def _get_output_override(
    *,
    phase_output_paths: dict[str, Any],
    phase_name: str,
    output_name: str,
) -> str | None:
    phase_cfg = phase_output_paths.get(str(phase_name), {})
    if not isinstance(phase_cfg, dict):
        return None

    outputs_cfg = phase_cfg.get("outputs", phase_cfg)
    if not isinstance(outputs_cfg, dict):
        return None

    if str(output_name) in _MULTISEG_OUTPUT_NAMES:
        multiseg_cfg = outputs_cfg.get("multiseg", {})
        if not isinstance(multiseg_cfg, dict):
            return None
        raw = multiseg_cfg.get(str(output_name), None)
    else:
        raw = outputs_cfg.get(str(output_name), None)

    if raw is None:
        return None

    candidate = Path(str(raw).strip())
    if not str(candidate).strip() or candidate.is_absolute() or ".." in candidate.parts:
        return None

    return str(candidate)


def coerce_phase_output_paths(value: Any) -> dict[str, dict[str, dict[str, str]]]:
    if not isinstance(value, dict) or not value:
        return {}

    cleaned: dict[str, dict[str, dict[str, str]]] = {}
    for phase_name, phase_cfg in value.items():
        if not isinstance(phase_cfg, dict):
            continue

        outputs_cfg = phase_cfg.get("outputs", phase_cfg)
        if not isinstance(outputs_cfg, dict):
            continue

        clean_outputs: dict[str, str] = {}
        clean_multiseg_outputs: dict[str, str] = {}
        for output_name, raw in outputs_cfg.items():
            if output_name == "multiseg":
                if not isinstance(raw, dict):
                    continue
                for ms_output_name, ms_raw in raw.items():
                    if ms_raw is None:
                        continue
                    ms_token = str(ms_raw).strip()
                    if not ms_token:
                        continue
                    ms_candidate = Path(ms_token)
                    if ms_candidate.is_absolute() or ".." in ms_candidate.parts:
                        continue
                    clean_multiseg_outputs[str(ms_output_name)] = ms_token
                continue

            if raw is None:
                continue
            token = str(raw).strip()
            if not token:
                continue
            candidate = Path(token)
            if candidate.is_absolute() or ".." in candidate.parts:
                continue
            clean_outputs[str(output_name)] = token

        if clean_outputs:
            cleaned[str(phase_name)] = {"outputs": clean_outputs}
        if clean_multiseg_outputs:
            phase_payload = cleaned.setdefault(str(phase_name), {"outputs": {}})
            phase_payload["outputs"]["multiseg"] = clean_multiseg_outputs

    return cleaned


def resolve_phase_output_path(
    *,
    output_dir: Path,
    phase_output_paths: dict[str, Any],
    phase_name: str,
    output_name: str,
    default_rel_path: str | None,
    logger: Any | None = None,
) -> Path:
    rel_path = _get_output_override(
        phase_output_paths=phase_output_paths,
        phase_name=phase_name,
        output_name=output_name,
    )

    if rel_path is None:
        if default_rel_path is None:
            raise ValueError(
                f"Missing required phase output path for {phase_name}.{output_name}; "
                "set option_kwargs.phase_output_paths accordingly."
            )
        rel_path = str(default_rel_path)

    resolved = Path(output_dir) / rel_path
    try:
        resolved.parent.mkdir(parents=True, exist_ok=True)
    except Exception:
        if logger is not None:
            try:
                logger.debug("Could not create parent dir for resolved phase output path %s", resolved, exc_info=True)
            except Exception:
                pass
    return resolved


def build_resolved_phase_output_paths(
    *,
    output_dir: Path,
    phase_output_paths: dict[str, Any],
    unitmatch_output_subdir_name: str,
    unitmatch_throughput_subdir_name: str,
    multiseg_mode: bool | None,
    stream_id: str | None = None,
    logger: Any | None = None,
) -> dict[str, dict[str, str]]:
    stream_token = str(stream_id).strip() if stream_id is not None else ""
    maxwell_default = (
        f"stg1_preprocess_outputs/maxwell_contiguous_epochs_{stream_token}.json"
        if stream_token
        else "stg1_preprocess_outputs/maxwell_contiguous_epochs.json"
    )
    concat_default = (
        f"stg1_preprocess_outputs/concatenation_stitch_epochs_{stream_token}.json"
        if stream_token
        else "stg1_preprocess_outputs/concatenation_stitch_epochs.json"
    )

    default_map: dict[str, dict[str, str]] = {
        "preprocessing": {
            "binary": "binary",
            "preprocessed_recording": "stg1_preprocess_outputs/preprocessed_recording",
            "maxwell_epochs": maxwell_default,
            "concat_epochs": concat_default,
            "preprocess_config": "stg1_preprocess_outputs/preprocess_config.json",
        },
        "spikesorting": {
            "sorter_output": "sorter_output",
        },
        "merge": {},
        "analyzer": {
            "analyzer_output": "analyzer_output",
        },
        "reports": {
            "qm_unfiltered": "qm_unfiltered.xlsx",
            "tm_unfiltered": "tm_unfiltered.xlsx",
            "metrics_curated": "metrics_curated.xlsx",
            "tm_curated": "tm_curated.xlsx",
        },
        "curation": {
            "rejection_log": "rejection_log.xlsx",
        },
    }

    required_multiseg_outputs = [
        ("preprocessing", "common_electrodes"),
        ("preprocessing", "maxwell_epochs"),
        ("preprocessing", "concat_epochs"),
        ("merge", "final_merged_sorting"),
        ("analyzer", "concat_waveforms"),
        ("analyzer", "segment_waveforms"),
        ("reports", "wf_rejection_log"),
        ("curation", "curated_sorting"),
    ]

    if _is_multiseg_enabled(multiseg_mode):
        missing: list[str] = []
        for phase_name, output_name in required_multiseg_outputs:
            if _get_output_override(
                phase_output_paths=phase_output_paths,
                phase_name=phase_name,
                output_name=output_name,
            ) is None:
                missing.append(f"{phase_name}.{output_name}")

        if missing:
            raise ValueError(
                "multiseg_mode requires explicit phase_output_paths for: " + ", ".join(missing)
            )

    resolved: dict[str, dict[str, str]] = {}
    for phase_name, outputs in default_map.items():
        phase_payload: dict[str, str] = {}
        for output_name, rel in outputs.items():
            phase_payload[str(output_name)] = str(
                resolve_phase_output_path(
                    output_dir=output_dir,
                    phase_output_paths=phase_output_paths,
                    phase_name=phase_name,
                    output_name=output_name,
                    default_rel_path=rel,
                    logger=logger,
                )
            )
        resolved[str(phase_name)] = phase_payload

    # Optional outputs are included only when explicitly configured.
    optional_outputs = [
        ("preprocessing", "common_electrodes"),
        ("preprocessing", "maxwell_epochs"),
        ("preprocessing", "concat_epochs"),
        ("merge", "final_merged_sorting"),
        ("analyzer", "concat_waveforms"),
        ("analyzer", "segment_waveforms"),
        ("reports", "wf_rejection_log"),
        ("curation", "curated_sorting"),
    ]
    for phase_name, output_name in optional_outputs:
        override = _get_output_override(
            phase_output_paths=phase_output_paths,
            phase_name=phase_name,
            output_name=output_name,
        )
        if override is None:
            continue

        phase_payload = resolved.setdefault(str(phase_name), {})
        phase_payload[str(output_name)] = str(
            resolve_phase_output_path(
                output_dir=output_dir,
                phase_output_paths=phase_output_paths,
                phase_name=phase_name,
                output_name=output_name,
                default_rel_path=None,
                logger=logger,
            )
        )

    return resolved
