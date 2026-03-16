from __future__ import annotations

from dataclasses import dataclass
import unittest
from unittest import mock

from IPNAnalysis.multiseg_utils.preprocess_multiseg_h5 import runner as m


@dataclass
class DummyRecording:
    n_segments: int

    def get_num_segments(self) -> int:
        return int(self.n_segments)


class DummySI:
    @staticmethod
    def select_segment_recording(recording: DummyRecording, i: int) -> DummyRecording:
        return DummyRecording(n_segments=1)

    @staticmethod
    def concatenate_recordings(recordings: list[DummyRecording]) -> DummyRecording:
        return DummyRecording(n_segments=1)


class TestPreprocessMultisegH5(unittest.TestCase):
    def test_expect_multi_true_with_single_segment_raises(self) -> None:
        rec = DummyRecording(n_segments=1)
        with self.assertRaises(RuntimeError) as ctx:
            m.prepare_multisegment_recording(rec, expect_multisegment=True, mode="none")
        self.assertIn("expect_multisegment=true", str(ctx.exception))

    def test_expect_multi_false_with_multi_segment_raises(self) -> None:
        rec = DummyRecording(n_segments=3)
        with self.assertRaises(RuntimeError) as ctx:
            m.prepare_multisegment_recording(rec, expect_multisegment=False, mode="none")
        self.assertIn("expect_multisegment=false", str(ctx.exception))

    def test_mode_none_preserves_segment_count(self) -> None:
        rec = DummyRecording(n_segments=4)
        out, info = m.prepare_multisegment_recording(rec, expect_multisegment=None, mode="none")

        self.assertIs(out, rec)
        self.assertEqual(info["segments_input"], 4)
        self.assertEqual(info["segments_output"], 4)
        self.assertEqual(info["mode"], "none")

    def test_mode_concatenate_collapses_to_single_segment(self) -> None:
        rec = DummyRecording(n_segments=5)
        with mock.patch.object(m, "si", DummySI()):
            out, info = m.prepare_multisegment_recording(rec, expect_multisegment=True, mode="concatenate")

        self.assertIsInstance(out, DummyRecording)
        self.assertEqual(out.get_num_segments(), 1)
        self.assertEqual(info["segments_input"], 5)
        self.assertEqual(info["segments_output"], 1)
        self.assertEqual(info["mode"], "concatenate")


if __name__ == "__main__":
    unittest.main()
