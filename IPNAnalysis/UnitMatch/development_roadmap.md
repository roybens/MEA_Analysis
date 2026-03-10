# UnitMatch + DeepUnitMatch Integration Roadmap

## Goal
Integrate UnitMatch/DeepUnitMatch into the MEA pipeline with minimal edits to `IPNAnalysis/mea_analysis_routine.py` by placing implementation details in a dedicated module under `IPNAnalysis/UnitMatch/`.

## Constraints From Project Direction
- Keep `mea_analysis_routine.py` edits minimal.
- Insert matching+merge flow after spikesorting and before curation/report logic.
- Accept a SpikeInterface sorting object from the pipeline.
- Run DeepUnitMatch-based matching logic.
- Report candidate/matched units clearly in logs + output artifacts.
- Merge matched units using SpikeInterface-supported methods.
- Continue downstream routine using updated sorting object.

## Proposed Architecture

### 1. New UnitMatch integration layer (this subdir)
Create focused modules under `IPNAnalysis/UnitMatch/`:
- `config.py`: integration config dataclass + defaults.
- `runner.py`: top-level orchestration entrypoint called by pipeline.
- `adapters.py`: conversion between SpikeInterface analyzer/sorting structures and UnitMatch/DeepUnitMatch inputs.
- `matching.py`: UnitMatch/DeepUnitMatch invocation and match scoring.
- `merge.py`: apply accepted matches through SpikeInterface merge APIs.
- `reporting.py`: save JSON/CSV summaries and human-readable logs.
- `exceptions.py`: explicit exception types to allow safe fallback behavior.

This keeps all UnitMatch concerns isolated and easy to disable.

### 2. Minimal touchpoint in `mea_analysis_routine.py`
Add one integration hook near the end of spikesorting / before curation steps:
- Call a single function, e.g. `run_unitmatch_merge_if_enabled(...)`.
- Pass required context: `sorting`, optional `recording`/`analyzer`, output dir, logger, and config flags.
- Receive updated `sorting` (merged) + summary metadata.
- Store summary in checkpoint/state for reproducibility.

## Execution Flow (Predicted)
1. Pipeline completes spikesorting and has `self.sorting` (+ optional `self.analyzer`).
2. If UnitMatch integration is enabled:
   - Build UnitMatch input payload from SpikeInterface objects.
   - Run DeepUnitMatch scoring/inference.
   - Convert pair scores into merge candidates using explicit thresholds.
   - Validate merge set (no cycles/conflicts, quality checks).
   - Apply merges via SpikeInterface merge utilities.
   - Return merged sorting object and summary.
3. Pipeline continues analyzer/curation/reports using merged sorting.
4. Write artifacts under a dedicated folder, e.g. `output_dir/unitmatch/`.

## Proposed Configuration Surface
Initial flags (CLI/config, all default-off):
- `unitmatch_enabled: bool = False`
- `unitmatch_mode: str = "deepunitmatch"` (future extensibility)
- `unitmatch_score_threshold: float`
- `unitmatch_max_pairs_per_unit: int | None`
- `unitmatch_dry_run: bool` (report only, no merge)
- `unitmatch_fail_open: bool = True` (on error: log + continue without merge)

Rationale: this protects collaborative stability while enabling controlled rollout.

## Artifact Plan
Under `.../unitmatch/` write:
- `unitmatch_config.json` (resolved config)
- `match_candidates.csv` (all evaluated pairs + scores)
- `accepted_matches.csv` (pairs selected for merge)
- `merge_actions.json` (exact merge operations applied)
- `unitmatch_summary.json` (counts, thresholds, timing, warnings)

## Validation Strategy

### Phase 0: Wiring only
- Hook call from pipeline with `enabled=False` default.
- Confirm zero behavior change when disabled.

### Phase 1: Dry-run matching
- Run DeepUnitMatch and generate reports only.
- No merges applied.
- Validate score distributions and candidate plausibility.

### Phase 2: Controlled merge enablement
- Enable merge for small test wells / subset.
- Validate unit count transitions and downstream analyzer compatibility.
- Confirm curation/report stages run unchanged with merged sorting.

### Phase 3: Wider rollout
- Add guardrails + stronger tests.
- Optional tuning presets per data modality.

## Testing Plan
- Unit tests for:
  - adapter conversion correctness,
  - merge-candidate conflict resolution,
  - fail-open behavior.
- Integration tests for:
  - end-to-end dry-run on a small fixture,
  - merge-on path updates sorting and preserves downstream execution.
- Regression check:
  - pipeline output unchanged when UnitMatch disabled.

## Risks and Mitigations
- DeepUnitMatch API/package instability:
  - isolate imports in `matching.py` and fail-open with actionable logs.
- Over-merging risk:
  - start dry-run default, require explicit thresholding.
- SpikeInterface merge API changes:
  - isolate SI-specific calls in `merge.py` adapter layer.
- Runtime/memory overhead:
  - cache intermediates and expose optional limits.

## First Implementation Milestones
1. Add module skeleton + dataclass config + reporting scaffolding.
2. Add one pipeline hook call with default-off behavior.
3. Implement dry-run matching path + artifacts.
4. Implement merge application path behind explicit enable flag.
5. Add tests + docs section for new flags and outputs.

## Definition of Done (Initial Feature)
- A single optional call in `mea_analysis_routine.py` triggers UnitMatch integration.
- When disabled, behavior is unchanged.
- When enabled in dry-run, candidate reports are produced.
- When merge enabled, sorting object is updated and downstream pipeline completes.
- Outputs are reproducible and documented.

## Progress Notes

### 2026-03-10 (Phase 0 wiring)
- Completed: created `IPNAnalysis/UnitMatch/` integration module scaffold.
- Completed: added `IPNAnalysis/UnitMatch/runner.py` with `run_unitmatch_merge_if_enabled(...)` and `UnitMatchConfig`.
- Completed: added `IPNAnalysis/UnitMatch/__init__.py` export surface.
- Completed: wired `run_analyzer()` in `mea_analysis_routine.py` to call UnitMatch path as an alternative to `auto_merge_units`.
- Completed: added CLI flags `--unitmatch-merge-units` and `--unitmatch-dry-run`.
- Completed: conflict handling when both UnitMatch and auto-merge are requested (UnitMatch takes precedence, warning logged).
- Completed: added per-run UnitMatch artifact output at `output_dir/unitmatch/unitmatch_summary.json`.
- Completed: removed duplicate `argparse` definitions for `--auto-merge-units` that were causing CLI conflict errors.
- Completed: extracted merge logic to explicit `Phase 2.5` method (`run_optional_merge_phase`) above analyzer stage.
- Completed: main pipeline order updated to `run_sorting()` -> `run_optional_merge_phase()` -> `run_analyzer()`.
- Completed: analyzer stage now computes the final analyzer once on canonical post-merge sorting.
- Completed: wired Phase 2.5 into checkpoint metadata via `self._save_checkpoint(..., merge_phase={...})`.
- Completed: default skipped behavior now records explicit merge status (`skipped_disabled` or `skipped_analyzer_complete`).
- Completed: added explicit checkpoint enum stages for merge (`MERGE`, `MERGE_COMPLETE`) in `ProcessingStage`.
- Completed: added checkpoint schema migration for legacy runs (v1 stage values shifted to account for merge stage insertion).

### Next
- Implement Phase 1 in `UnitMatch/runner.py`: DeepUnitMatch-backed scoring and candidate export (`match_candidates.csv`) while keeping dry-run/no-merge behavior.
- Define first-pass pair acceptance policy and conflict pruning before enabling real merges.
