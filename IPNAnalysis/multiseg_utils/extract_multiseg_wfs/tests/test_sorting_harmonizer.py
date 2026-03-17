from __future__ import annotations

import unittest

from IPNAnalysis.multiseg_utils.extract_multiseg_wfs import harmonize_waveform_sorting


class _DummySorting:
    def __init__(self, unit_ids):
        self._unit_ids = list(unit_ids)

    def get_unit_ids(self):
        return list(self._unit_ids)

    @property
    def unit_ids(self):
        return list(self._unit_ids)

    def select_units(self, unit_ids=None):
        if unit_ids is None:
            return _DummySorting(self._unit_ids)
        return _DummySorting(list(unit_ids))


class _DummyRecording:
    pass


class TestSortingHarmonizer(unittest.TestCase):
    def test_debug_limit_keeps_sorted_first_n_units(self) -> None:
        sorting = _DummySorting([7, 2, 9, 1])
        recording = _DummyRecording()

        out = harmonize_waveform_sorting(
            sorting=sorting,
            recording=recording,
            debug_max_units=2,
            logger=None,
        )

        self.assertEqual(out.get_unit_ids(), [1, 2])

    def test_none_debug_limit_leaves_units_unchanged(self) -> None:
        sorting = _DummySorting([3, 1, 2])
        recording = _DummyRecording()

        out = harmonize_waveform_sorting(
            sorting=sorting,
            recording=recording,
            debug_max_units=None,
            logger=None,
        )

        self.assertEqual(out.get_unit_ids(), [3, 1, 2])


if __name__ == "__main__":
    unittest.main()