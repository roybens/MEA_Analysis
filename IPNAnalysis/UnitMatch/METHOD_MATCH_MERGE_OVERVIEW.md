# UnitMatch Matching And Merge Method (Current Implementation)

## Purpose
This module performs optional post-sorting unit matching and merge suggestion/application for a single well. It can run in report-only mode, scored dry-run mode, or apply selected merges directly to the sorting object.

The implementation is centered in runner.py and uses a DeepUnitMatch clone adapter backend for similarity scoring.

## High-Level Pipeline
1. Build candidate unit pairs from current sorting unit IDs.
2. Optionally score pairs with DeepUnitMatch-derived similarity matrix.
3. Produce candidate and thresholded match tables.
4. Build oversplit merge suggestions using forward and reciprocal score gating.
5. Select conflict-free top suggestions.
6. Optionally apply merges to sorting.
7. Persist artifacts and summary JSON.
8. If recursive mode is enabled, repeat iterations until convergence/stop conditions.

## Inputs And Configuration
Primary config object is UnitMatchConfig with controls for:
- Enable/disable mode.
- Dry-run versus merge-apply behavior.
- Candidate pair and suggestion limits, including unlimited mode via negative values.
- Recursive iteration behavior and max iterations.
- Spike sampling cap for RawWaveforms generation (negative means uncapped).
- Backend selection and clone root.

## Step-By-Step Method

### 1) Candidate Pair Construction
- Unit IDs are normalized to string form.
- All unique unordered pairs are built by nested traversal.
- Pair generation stops early when max_candidate_pairs is reached.
- Limit semantics:
  - max < 0 means unlimited.
  - max = 0 means no pairs.
  - max > 0 means bounded.

### 2) Scoring Backend Path
If scored dry-run is enabled, the backend path is attempted:

1. Verify DeepUnitMatch clone availability and clone-only import validity.
2. Locate valid sorter output directory containing channel_positions.npy.
3. Ensure RawWaveforms are available for this iteration:
   - If missing, generate from SpikeInterface sorting/recording using split-halves.
   - For each half, compute random spikes, waveforms, templates.
   - Resample template time axis to 82 samples (compatibility target) and save UnitX_RawSpikes.npy.
   - max_spikes_per_unit is applied when non-negative; negative means uncapped sampling.
4. Load RawWaveforms files and validate shapes.
5. Preprocess waveforms through DeepUnitMatch preprocessing utility.
6. Run DeepUnitMatch inference to obtain similarity matrix.
7. Convert candidate pairs into directional scores from matrix lookup:
   - score(a,b)
   - reverse_score(b,a)
   - reciprocal_min_score = min(score, reverse_score)
8. Persist matrix artifacts for diagnostics.

If backend scoring fails, flow degrades to unscored placeholders and continues in fail-open style unless configured otherwise.

### 3) Candidate Tables And Thresholding
- match_candidates.csv is written for all candidate pairs.
- match_candidates_thresholded.csv is written with threshold comparison against oversplit_min_probability.
- Summary counts track good, bad, and unscored matches.

### 4) Oversplit Suggestion Gate
Suggestions are produced only when all gates pass:
1. score is numeric and not NaN.
2. forward score >= threshold.
3. reverse score is numeric and not NaN.
4. reciprocal minimum score >= threshold.
5. unit IDs are distinct and valid.

Suggestions are sorted by score descending and optionally truncated by oversplit_max_suggestions:
- max < 0 unlimited.
- max = 0 none.
- max > 0 bounded.

### 5) Conflict-Free Pair Selection
A second pass selects non-conflicting pairs greedily:
- Suggestions are ordered by score (with rank tie handling).
- A pair is selected only if neither unit has already been used.
- This produces selected_non_conflicting_pairs.csv.

This ensures one merge partner per unit within an iteration.

### 6) Merge Application
If apply_merges is true and dry_run is false:
1. Convert selected unit IDs back to live sorting unit IDs.
2. Build merge groups as pairs.
3. Apply merges with SpikeInterface apply_merges_to_sorting.
4. On success, persist final merged sorting artifact to:
   - unitmatch_outputs/final_merged_sorting

If merge application cannot proceed, merge_application reason is captured and flow remains fail-open by default.

### 7) Iterative Recursive Mode
run_unitmatch_merge_with_recursion executes repeated single-iteration runs when recursive is enabled.

Per run:
- Starts from clean iteration workspace (iteration_N folders).
- Runs one iteration using non-recursive config clone.
- Reads selected pairs signature.

Stop conditions:
1. apply_merges is disabled.
2. No selected pairs.
3. Selected pair signature repeats previous iteration.
4. Merge not applied.
5. Max iterations reached (unless max_iterations < 0, which is uncapped).

Outputs include merge_history.json and an iteration convergence report embedded in summary.

## Artifacts Produced
Main outputs under output_subdir_name and throughput_subdir_name:
- unitmatch_summary.json
- unitmatch_config.json
- match_candidates.csv
- match_candidates_thresholded.csv
- oversplit_suggestions.csv
- selected_non_conflicting_pairs.csv
- per_unit_partner_summary.json
- deepunitmatch_similarity_matrix.npy
- deepunitmatch_similarity_matrix_unit_ids.json
- final_merged_sorting folder (when merges applied)
- merge_history.json (recursive mode)

## Reporting Layer
reporting.py creates a static report pack by consuming summary and artifacts:
- Score histogram and CDF.
- Similarity matrix heatmap.
- Partner suggestion tables/plots.
- Merge diagnostics counts.
- Iteration convergence plots.
- report_manifest.json.

## Operational Characteristics
- Fail-open behavior is intentional in many error paths so downstream pipeline can continue.
- Backend path enforces clone-only import source checks to avoid accidental environment mismatches.
- Reciprocal gating is conservative and reduces one-direction-only false positives.
- Conflict-free greedy selection favors top confidence while preventing multi-merge conflicts in one pass.

## Diagram-Friendly Flow Skeleton
Input sorting + recording
-> candidate pair generation
-> scored backend path (or unscored fallback)
-> thresholding + reciprocal gate
-> oversplit suggestions
-> conflict-free selection
-> optional merge application
-> iteration convergence check
-> artifacts + report generation
