from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import h5py

from IPNAnalysis.h5_utils.streams import list_stream_recording_names


class TestH5Streams(unittest.TestCase):
    def test_list_stream_recording_names_single_segment(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            h5_path = Path(td) / "single.raw.h5"
            with h5py.File(h5_path, "w") as h5:
                h5.create_group("wells").create_group("well000").create_group("rec0000")

            names = list_stream_recording_names(h5_path=h5_path, stream_id="well000")
            self.assertListEqual(names, ["rec0000"])

    def test_list_stream_recording_names_multisegment(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            h5_path = Path(td) / "multi.raw.h5"
            with h5py.File(h5_path, "w") as h5:
                well = h5.create_group("wells").create_group("well001")
                well.create_group("rec0000")
                well.create_group("rec0001")

            names = list_stream_recording_names(h5_path=h5_path, stream_id="well001")
            self.assertSetEqual(set(names), {"rec0000", "rec0001"})


if __name__ == "__main__":
    unittest.main()
