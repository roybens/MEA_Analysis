# NetworkAnalysis

MATLAB-based tools for network burst analysis on Maxwell Biosystems MEA recordings.

---

## MXW Toolbox Functions Used in the NetworkAnalysis Workflow

The table below lists every function from the MaxWell MATLAB Toolbox (`mxw`) called within the `NetworkAnalysis/` scripts.

### Core

| Function | Description | Used in |
|---|---|---|
| `mxw.fileManager(path, wellID)` | Creates a `fileManager` object that loads an HDF5 recording file and exposes spike data, electrode maps, and metadata. | `NetworkMetricsSingleFolder.m`, `ActivityMetricsSingleFolder.m`, `Part_4_Network.m`, `compileNetworkFiles.m`, `compileActivityFiles.m`, and others |

### Network Activity

| Function | Description | Used in |
|---|---|---|
| `mxw.networkActivity.computeNetworkAct(networkData, 'BinSize', b, 'GaussianSigma', s)` | Bins all spike times into a histogram, convolves with a Gaussian kernel to produce a smoothed network firing-rate curve, and returns a struct with `time` and `firingRate` fields. | `NetworkMetricsSingleFolder.m`, `Part_4_Network.m`, `compileNetworkFiles.m`, `Compare_NetworkParameters.m` |
| `mxw.networkActivity.computeNetworkStats(networkAct, 'Threshold', t)` | Detects network bursts in the firing-rate curve and returns burst statistics: peak amplitudes, peak times, interburst intervals, etc. | `Part_4_Network.m`, `Copy_of_Part_4_Network.m`, `Compare_NetworkParameters.m` |

### Activity Map

| Function | Description | Used in |
|---|---|---|
| `mxw.activityMap.computeSpikeRate(fileManagerObj)` | Returns the mean per-electrode firing rate (Hz) across the recording. | `ActivityMetricsSingleFolder.m`, `Part_2_Spikes.m`, `compileActivityFiles.m`, `RastorPlotSegregation.m` |
| `mxw.activityMap.computeSpikeCount(fileManagerObj)` | Returns the total spike count per electrode across the recording. | `RastorPlotSegregation.m` |
| `mxw.activityMap.computeAmplitude90percentile(fileManagerObj)` | Returns the 90th-percentile spike amplitude (µV) per electrode. | `ActivityMetricsSingleFolder.m`, `Part_2_Spikes.m`, `compileActivityFiles.m` |

### Plotting

| Function | Description | Used in |
|---|---|---|
| `mxw.plot.rasterPlot(networkData, ...)` | Plots a spike raster (channel vs. time) for a network recording. | `NetworkMetricsSingleFolder.m`, `Part_4_Network.m`, `Copy_of_Part_4_Network.m` |
| `mxw.plot.networkActivity(networkAct, 'Threshold', t, ...)` | Plots the smoothed network firing-rate curve together with the burst-detection threshold line. | `NetworkMetricsSingleFolder.m`, `Part_4_Network.m`, `Copy_of_Part_4_Network.m` |
| `mxw.plot.networkStats(networkStats, 'Option', opt, ...)` | Plots histograms of burst statistics (e.g. peak amplitude or interburst interval). | `Part_4_Network.m`, `Copy_of_Part_4_Network.m`, `Hyperbursts/hyperbursts.m` |
| `mxw.plot.activityMap(fileManagerObj, metric, ...)` | Plots a 2-D spatial heatmap of a per-electrode metric (e.g. firing rate or amplitude) overlaid on the electrode array layout. | `ActivityMetricsSingleFolder.m`, `Part_2_Spikes.m`, `compileActivityFiles.m` |
| `mxw.plot.axonTraces(x, y, waveformCutOuts, ...)` | Plots single-neuron "footprint" traces from spike waveform cut-outs across all electrodes. | `MaxwellAnalysisBasicScripts/Part_3_Traces.m` |

### Utilities

| Function | Description | Used in |
|---|---|---|
| `mxw.util.computeRelativeSpikeTimes(fileManagerObj)` | Converts raw spike frame numbers to timestamps in seconds relative to the start of the recording; returns a struct with `time` and `channel` vectors. | `NetworkMetricsSingleFolder.m`, `Part_4_Network.m`, `workspace.m`, `hyperbursts.m`, and others |
| `mxw.util.normpdf(x, mu, sigma)` | Evaluates the normal (Gaussian) probability density function; used to construct the convolution kernel for firing-rate smoothing. | `Part_4_Network.m`, `SpikeSortedNetworkActivitySingleFile.m`, `hyperbursts.m`, and others |
| `mxw.util.rms(signal)` | Computes the root-mean-square (RMS) of a signal vector; used to set the adaptive burst-detection threshold. | `Part_4_Network.m`, `SpikeSortedNetworkActivitySingleFile.m`, `computeNetworkStatsModified.m`, and others |
| `mxw.util.findPeaks(signal, 'PositiveThreshold', t)` | Detects positive peaks above a threshold in a signal vector; returns peak indices and values. | `Part_4_Network.m`, `Copy_of_Part_4_Network.m` |
| `mxw.util.percentile(data, p)` | Returns the *p*-th percentile of a data vector; used to cap color scales in activity-map plots. | `ActivityMetricsSingleFolder.m`, `Part_2_Spikes.m` |
| `mxw.util.groupElectrodes(positions, deltaDist)` | Groups electrodes whose spatial positions are within `deltaDist` of each other; used for redundant-spike detection. | `MaxwellAnalysisBasicScripts/detRedspike.m` |

---

## Summary

In total, **17 distinct `mxw` toolbox functions** are actively used across the NetworkAnalysis workflow:

```
mxw.fileManager
mxw.networkActivity.computeNetworkAct
mxw.networkActivity.computeNetworkStats
mxw.activityMap.computeSpikeRate
mxw.activityMap.computeSpikeCount
mxw.activityMap.computeAmplitude90percentile
mxw.plot.rasterPlot
mxw.plot.networkActivity
mxw.plot.networkStats
mxw.plot.activityMap
mxw.plot.axonTraces
mxw.util.computeRelativeSpikeTimes
mxw.util.normpdf
mxw.util.rms
mxw.util.findPeaks
mxw.util.percentile
mxw.util.groupElectrodes
```
