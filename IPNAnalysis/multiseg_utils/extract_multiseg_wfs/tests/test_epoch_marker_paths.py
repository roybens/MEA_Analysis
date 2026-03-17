from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from IPNAnalysis.multiseg_utils.extract_multiseg_wfs import resolve_waveform_epoch_marker_paths


class TestEpochMarkerPaths(unittest.TestCase):
    def test_prefers_local_epoch_markers_when_present(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            base = Path(td)
            stream_id = "well001"

            local_maxwell = base / f"maxwell_contiguous_epochs_{stream_id}.json"
            local_concat = base / f"concatenation_stitch_epochs_{stream_id}.json"
            local_maxwell.write_text("[]", encoding="utf-8")
            local_concat.write_text("[]", encoding="utf-8")

            resolved = resolve_waveform_epoch_marker_paths(
                output_dir=base,
                stream_id=stream_id,
            )

            self.assertEqual(resolved["maxwell"], local_maxwell)
            self.assertEqual(resolved["concat"], local_concat)

    def test_falls_back_to_preprocess_dir_when_local_missing(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            base = Path(td)
            stream_id = "well002"

            preprocess_dir = base / "stg1_preprocess_outputs"
            preprocess_dir.mkdir(parents=True, exist_ok=True)
            nested_maxwell = preprocess_dir / f"maxwell_contiguous_epochs_{stream_id}.json"
            nested_concat = preprocess_dir / f"concatenation_stitch_epochs_{stream_id}.json"

            resolved = resolve_waveform_epoch_marker_paths(
                output_dir=base,
                stream_id=stream_id,
            )

            self.assertEqual(resolved["maxwell"], nested_maxwell)
            self.assertEqual(resolved["concat"], nested_concat)


if __name__ == "__main__":
    unittest.main()