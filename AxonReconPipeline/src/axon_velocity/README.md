[![PyPI version](https://badge.fury.io/py/axon-velocity.svg)](https://badge.fury.io/py/axon-velocity) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4896745.svg)](https://doi.org/10.5281/zenodo.4896745)

# Axon velocity

Graph-based algorithm to reconstruct axonal branches and compute velocities from extracellular 
multi-electrode array (MEA) recordings in Python.

## Installation 

To install the `axon_velocity` package, you can clone the repo and install using pip:

```bash
pip install axon_velocity
```

To install from sources:

```bash
git clone https://github.com/alejoe91/axon_velocity.git
cd axon_velocity
python setup.py install (or develop)
```

### Requirements

`axon_velocity` depends on the following packages, which are automatically installed

- numpy
- matplotlib
- scipy
- networkx
- sklearn
- MEAutility
- probeinterface

For the simulation notebooks in the `simulation_notebooks` folder, additional requirements are needed:

- [NEURON](https://www.neuron.yale.edu/neuron/)
- [LFPy](https://lfpy.readthedocs.io/en/latest/)
- [neuroplotlib](https://github.com/alejoe91/neuroplotlib)

All additional requirements can be installed with: `pip install -r requirements_fill.txt`

## Usage

The inputs to the tracking algorithm are:

- templates: mean extracellular waveforms (n_channels x n_samples)
- locations: x-y position of electrodes (n_channels x 2)
- fs: sampling frequency (float)

The graph-based method can be run as follows:

```python
import axon_velocity as av

gtr = av.compute_graph_propagation_velocity(template=your_template, locations=your_locations, fs=fs)
```

To inspect available arguments, you can use `av.compute_graph_propagation_velocity?`. 

The output `gtr` is an object of a class called `GraphAxonTracking`. 
It contains the following fields:

- `branches`: List of dictionaries containing the following fields:
    - 'selected_channels': selected channels in the path
    - 'velocity': velocity estimate in mm/s (if locations in um and fs in Hz)
    - 'offset': offset of velocity estimate (same units as locations)
    - 'r2': r-squared of the velocoty fit
    - 'error': standard error of the linear fit
    - 'pval': p_value of the fit
    - 'distances': array with distances computed along the branch
    - 'peak_times': array with peak differences with initial channel
    - 'init_channel': channel used as initial channel
- `selected_channels`: List of selected channels
- `graph`: NetworkX directed graph 

The `GraphAxonTracking` also implements useful methods for plotting the selected channels 
(`gtr.plot_channel_selection()`), plot the underlying graph (`gtr.plot_graph()`), plot the selected axonal branches 
(`gtr.plot_branches()`), and plot the estimated velocities for each branch (`gtr.plot_velocities()`).
    
  
## Contribute

Contributions are welcome! Before pushing, make sure to clean up all notebooks with `nbconvert`:

`pip install nbconvert` (just once)

`jupyter nbconvert --to notebook --ClearOutputPreprocessor.enabled=True --ClearMetadataPreprocessor.enabled=True  --inplace **/*.ipynb` (before committing) 
