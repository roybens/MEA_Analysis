from __future__ import annotations

from setuptools import setup


# This repository is both the VCS root and the top-level Python package directory.
# We provide explicit package mappings so `pip install -e .` supports imports like:
#   import MEA_Analysis
#   from MEA_Analysis.IPNAnalysis import mea_analysis_routine
setup(
    name="MEA_Analysis",
    version="0.0.0",
    packages=[
        "MEA_Analysis",
        "MEA_Analysis.IPNAnalysis",
    ],
    package_dir={
        "MEA_Analysis": ".",
        "MEA_Analysis.IPNAnalysis": "IPNAnalysis",
    },
)
