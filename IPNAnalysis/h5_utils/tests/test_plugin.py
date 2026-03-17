from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path

from IPNAnalysis.h5_utils.plugin import ensure_maxwell_hdf5_plugin_path


class TestH5Plugin(unittest.TestCase):
    def test_invalid_env_path_is_cleared(self) -> None:
        old = os.environ.get("HDF5_PLUGIN_PATH")
        try:
            os.environ["HDF5_PLUGIN_PATH"] = "/definitely/not/a/real/path"
            ensure_maxwell_hdf5_plugin_path(prefix="[test]")
            current = os.environ.get("HDF5_PLUGIN_PATH")
            self.assertTrue(current is None or Path(current).exists())
        finally:
            if old is None:
                os.environ.pop("HDF5_PLUGIN_PATH", None)
            else:
                os.environ["HDF5_PLUGIN_PATH"] = old

    def test_valid_env_path_is_kept(self) -> None:
        old = os.environ.get("HDF5_PLUGIN_PATH")
        with tempfile.TemporaryDirectory() as td:
            try:
                os.environ["HDF5_PLUGIN_PATH"] = td
                ensure_maxwell_hdf5_plugin_path(prefix="[test]")
                self.assertEqual(os.environ.get("HDF5_PLUGIN_PATH"), td)
            finally:
                if old is None:
                    os.environ.pop("HDF5_PLUGIN_PATH", None)
                else:
                    os.environ["HDF5_PLUGIN_PATH"] = old


if __name__ == "__main__":
    unittest.main()
