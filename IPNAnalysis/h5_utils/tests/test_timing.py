from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

from IPNAnalysis.h5_utils.timing import read_well_rec_frame_nos_and_trigger_settings


class TestH5Timing(unittest.TestCase):
    def test_read_well_rec_frame_nos_and_trigger_settings(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            h5_path = Path(td) / "test.raw.h5"
            with h5py.File(h5_path, "w") as h5:
                wells = h5.create_group("wells")
                well = wells.create_group("well001")
                rec = well.create_group("rec0000")
                rec.create_dataset("start_time", data=np.array([1000], dtype=np.int64))
                rec.create_dataset("stop_time", data=np.array([2000], dtype=np.int64))
                routed = rec.create_group("groups").create_group("routed")
                routed.create_dataset("frame_nos", data=np.array([10, 11, 12, 20, 21], dtype=np.int64))

            info = read_well_rec_frame_nos_and_trigger_settings(
                h5_path=h5_path,
                stream_id="well001",
                rec_name="rec0000",
            )

            self.assertEqual(info["start_ms"], 1000)
            self.assertEqual(info["stop_ms"], 2000)
            self.assertListEqual(info["frame_nos"].tolist(), [10, 11, 12, 20, 21])


if __name__ == "__main__":
    unittest.main()
