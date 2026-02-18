# MEA_Analysis

End-to-end pipeline for neuronal spike sorting and network burst analysis on **Maxwell Biosystems MEA** recordings.

Built on [SpikeInterface](https://github.com/SpikeInterface/spikeinterface) with Kilosort4 as the default sorter.

---

## Repository Layout

```
MEA_Analysis/
├── IPNAnalysis/          ← Main production pipeline (start here)
├── NetworkAnalysis/      ← MATLAB-based network analysis tools
├── NeuronClassification/ ← UMAP-based neuron classification
├── Archive/              ← Legacy scripts (not maintained)
├── MEAProcessingLibrary/ ← Installable Python processing library
├── MaxwellBiosystemsDeviceInterface/  ← Device control scripts
├── StimulationAnalysis/  ← Stimulation experiment analysis
└── dockers/              ← Docker build configs for spike sorters
```

## Getting Started

The active pipeline lives in [`IPNAnalysis/`](IPNAnalysis/README.md).  
See the **[IPNAnalysis README](IPNAnalysis/README.md)** for:

- Architecture overview (driver → routine → config_loader)
- Pipeline stages (preprocessing, sorting, analyzer, reports)
- Configuration system (`mea_config.json`)
- Typical workflows and full CLI reference

### Quick Start

```bash
# Generate a config template
python IPNAnalysis/config_loader.py mea_config.json

# Dry run on a directory
python IPNAnalysis/run_pipeline_driver.py /data/experiment --config mea_config.json --dry

# Full batch run
python IPNAnalysis/run_pipeline_driver.py /data/experiment --config mea_config.json

# Single well
python IPNAnalysis/mea_analysis_routine.py /data/exp/run_001/Network/data.raw.h5 \
  --well well000 --rec rec0001 --config mea_config.json
```

## Dependencies

- Python ≥ 3.9
- [SpikeInterface](https://github.com/SpikeInterface/spikeinterface)
- Kilosort4 (GPU recommended, ≥ 8 GB VRAM)
- See `requirements.txt` for full list

## License

Pending

