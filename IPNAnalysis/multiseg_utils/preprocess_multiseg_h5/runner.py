from __future__ import annotations

import logging
from typing import Any

try:
    import spikeinterface.full as si
except Exception:  # pragma: no cover - environment dependent
    si = None


def _segment_count(recording: Any) -> int:
    try:
        return int(recording.get_num_segments())
    except Exception:
        return 1


def _as_optional_bool(value: Any) -> bool | None:
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    token = str(value).strip().lower()
    if token in {"", "none", "null", "auto"}:
        return None
    if token in {"true", "1", "yes", "y", "multi", "multiseg", "multisegment"}:
        return True
    if token in {"false", "0", "no", "n", "single", "single_segment"}:
        return False
    raise ValueError(
        "Invalid expect_multisegment value. Use true/false or single/multi tokens."
    )


def _coerce_multiseg_mode(mode: Any) -> bool:
    if mode is None:
        return False
    if isinstance(mode, bool):
        return mode

    token = str(mode).strip().lower()
    if token in {"", "none", "off", "disabled", "false", "0", "noop"}:
        return False
    if token in {"true", "1", "yes", "concat", "concatenate", "merge_segments"}:
        return True

    raise ValueError(
        f"Invalid multiseg_mode '{mode}'. Valid values: true, false"
    )


def _concatenate_segments(recording: Any) -> Any:
    n_segments = _segment_count(recording)
    if n_segments <= 1:
        return recording

    if si is None:
        raise RuntimeError(
            "multiseg_mode=concatenate requires spikeinterface to be installed"
        )

    segment_recordings = [si.select_segment_recording(recording, i) for i in range(n_segments)]
    return si.concatenate_recordings(segment_recordings)


def prepare_multisegment_recording(
    recording: Any,
    *,
    expect_multisegment: Any = None,
    mode: Any = False,
    logger: logging.Logger | None = None,
) -> tuple[Any, dict[str, Any]]:
    """Validate multisegment topology and optionally transform it.

    This is a strict policy gate used by stage-1 preprocessing. When
    ``expect_multisegment`` is provided, the function raises immediately if
    recording topology disagrees with the expectation.
    """

    log = logger or logging.getLogger(__name__)

    expected_multi = _as_optional_bool(expect_multisegment)
    mode_enabled = _coerce_multiseg_mode(mode)

    segments_input = _segment_count(recording)
    actual_multi = segments_input > 1

    if expected_multi is True and not actual_multi:
        raise RuntimeError(
            "Multisegment policy mismatch: expect_multisegment=true but input recording has a single segment."
        )
    if expected_multi is False and actual_multi:
        raise RuntimeError(
            "Multisegment policy mismatch: expect_multisegment=false but input recording has multiple segments."
        )

    out_recording = recording
    if mode_enabled:
        if segments_input <= 1:
            log.info("multiseg_mode=true requested but input has <=1 segment; no-op")
        else:
            log.info("Applying multiseg_mode=true across %d segments", segments_input)
            out_recording = _concatenate_segments(recording)

    segments_output = _segment_count(out_recording)
    info = {
        "mode": bool(mode_enabled),
        "segments_input": int(segments_input),
        "segments_output": int(segments_output),
        "expected_multi": expected_multi,
    }
    log.info(
        "Multisegment prep summary: mode=%s expected_multi=%s input_segments=%d output_segments=%d",
        bool(mode_enabled),
        expected_multi,
        int(segments_input),
        int(segments_output),
    )

    return out_recording, info


__all__ = ["prepare_multisegment_recording"]
