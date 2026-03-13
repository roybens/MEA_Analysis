# UnitMatch Development Decisions

## 82 timepoint convention for template waveforms

DeepUnitMatch has historically assumed Neuropixels-like waveform/template shapes, including a time-axis convention of 82 samples.

We patched DeepUnitMatch to better tolerate variable shapes for Maxwell HDMEA support, but in this pipeline we still generate template waveforms with 82 time points for compatibility and to reduce shape-related regressions.

Current implementation details:
- Target sampling frequency: 30000 Hz
- Target timepoints: 82
- Window length: (82 - 1) / 30000 s = 2.7 ms
- Extraction window split: 0.9 ms before spike, 1.8 ms after spike

Performance impact is currently unknown. This compatibility choice should be validated empirically on held-out data.
