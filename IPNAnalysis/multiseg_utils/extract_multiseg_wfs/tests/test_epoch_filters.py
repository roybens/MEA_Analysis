from __future__ import annotations

import unittest

from IPNAnalysis.multiseg_utils.extract_multiseg_wfs import (
    epochs_to_intervals,
    filter_spike_train_by_intervals,
    maxwell_epochs_to_segment_local_intervals,
)


class TestEpochFilters(unittest.TestCase):
    def test_epochs_to_intervals_sorts_and_filters_invalid(self) -> None:
        epochs = [
            {"start_sample": 30, "end_sample": 40},
            {"start_sample": 10, "end_sample": 20},
            {"start_sample": 50, "end_sample": 50},
            {"foo": 1},
        ]
        self.assertEqual(epochs_to_intervals(epochs), [(10, 20), (30, 40)])

    def test_maxwell_epochs_to_segment_local_intervals_filters_segment(self) -> None:
        maxwell_epochs = [
            {"segment_index": 0, "segment_start_sample": 5, "segment_end_sample": 10},
            {"segment_index": 1, "segment_start_sample": 1, "segment_end_sample": 2},
            {"segment_index": 0, "segment_start_sample": 20, "segment_end_sample": 25},
        ]
        self.assertEqual(
            maxwell_epochs_to_segment_local_intervals(maxwell_epochs=maxwell_epochs, segment_index=0),
            [(5, 10), (20, 25)],
        )

    def test_filter_spike_train_by_intervals(self) -> None:
        spike_train = [8, 12, 18, 24, 31]
        intervals = [(10, 20), (30, 40)]
        kept, removed_outside, removed_edge, removed_outside_spikes, removed_edge_spikes = (
            filter_spike_train_by_intervals(
                spike_train=spike_train,
                intervals=intervals,
                pre_samples=2,
                post_samples=2,
            )
        )

        self.assertEqual(kept, [12])
        self.assertEqual(removed_outside, 2)
        self.assertEqual(removed_edge, 2)
        self.assertEqual(removed_outside_spikes, [8, 24])
        self.assertEqual(removed_edge_spikes, [18, 31])


if __name__ == "__main__":
    unittest.main()
