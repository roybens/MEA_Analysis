# UnitMatch Setup

## Purpose
This folder contains the MEA_Analysis integration layer for UnitMatch/DeepUnitMatch.

## Install (current project practice)
The to use DeepUnitMatch we install from GitHub directly. This will install the latest code from the `main` branch of the UnitMatch repository, specifically targeting the `UnitMatchPy` subdirectory.

```bash
python -m pip install "git+https://github.com/EnnyvanBeest/UnitMatch.git@main#subdirectory=UnitMatchPy"
```

## Triggering dry-run from driver
From `run_pipeline_driver.py`:

```bash
python IPNAnalysis/run_pipeline_driver.py <path_to_data_or_dir> \
  --unitmatch-merge-units \
  --unitmatch-dry-run
```

This passes through to `mea_analysis_routine.py` and executes:
- Phase 2 sorting (or resume)
- Phase 2.5 merge optional stage in UnitMatch dry-run mode
- Analyzer/reports continue as usual

## Dry-run artifacts
Under each well output folder:
- `unitmatch/unitmatch_config.json`
- `unitmatch/unitmatch_summary.json`
- `unitmatch/match_candidates.csv`
