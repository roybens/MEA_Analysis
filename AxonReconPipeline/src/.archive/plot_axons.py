

from IPython.core.debugger import set_trace
from pathlib import Path

### Initiate Axon trace algo
import numpy as np
import matplotlib.pylab as plt
import MEAutility as mu
from scipy.signal import resample_poly
from scipy.stats import kurtosis, linregress
from matplotlib import gridspec
from scipy import io
import numpy as np
import networkx as nx
from pathlib import Path
import sys
from pprint import pprint
from axon_velocity import *
from axon_velocity.models import load_cell
from axon_velocity.evaluation import *
import os
from pathlib import Path
import pickle
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from probeinterface import plotting as plotting_probe
import json

import concurrent.futures
from tqdm import tqdm
from multiprocessing import Manager

#for testing
import numpy as np
from shapely.geometry import Polygon
import itertools
from shapely.geometry import Point
from shapely.geometry import Polygon
from collections import Counter
import math
from collections import defaultdict
import mea_processing_library as MPL

#import all functions in main.py
from main import *

#Define Algorithm Params
params = get_default_graph_velocity_params()
# change params (fig6 Buccino axon_velocity paper example)
params['detect_threshold'] = 0.01
params['kurt_threshold'] = 0.1
params['peak_std_threshold'] = 0.8
params['upsample'] = 5
params['neighbor_radius'] = 100
params['r2_threshold'] = 0.8

# AW invention. I dont remember why I did this. 4Jan2023
#params["theilsen_maxiter"] = 10000
#pprint(params)
#

#debug
debug_mode = True
pre_selected_folders = [        
    #"/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/Tests_Adam/Tests_Adam/240118/M06844/ActivityScan/000023",
    #"/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/Tests_Adam/Tests_Adam/240118/M06844/Network/000024",
    "/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/Tests_Adam/Tests_Adam/240118/M06844/AxonTracking/000026",
    #testing continuity test with non-continuous data
    #"/mnt/disk20tb/PrimaryNeuronData/Maxone/SPTAN1_1/230725/16657/ActivityScan",

    #new comparison 09Feb24
    #"/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/FolicAcid_AxonTracking_05022024/FolicAcid_AxonTrack_Comparison/240205/M05506/AxonTracking/000003"
    #"/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/FolicAcid_AxonTracking_05022024/FolicAcid_AxonTrack_Comparison/240205/M07038/AxonTracking"
    ]
continuous_h5_dirs = select_folders_and_get_continuous_h5_dirs(pre_selected_folders = pre_selected_folders, debug_mode = debug_mode)


stream_select = 1
templates_dir='./AxonReconPipeline/data/temp_data/templates'
for i, continuous_h5_dir in continuous_h5_dirs.items():
    recording_details = MPL.extract_recording_details(continuous_h5_dir['h5_file_path'])
    date = recording_details[0]['date']
    chip_id = recording_details[0]['chipID']
    scanType = recording_details[0]['scanType']
    run_id = recording_details[0]['runID']
    h5 = h5py.File(continuous_h5_dir['h5_file_path'], 'r')
    stream_ids = list(h5['wells'].keys())
    for stream_id in stream_ids:
        if stream_select is not None:
            if stream_id != f'well{stream_select:03}':
                continue
            else:
                print(f'Processing stream {stream_id}')
                break

    temp_dir = os.path.dirname(templates_dir)
    sorting_dir = Path(temp_dir) /'sortings'/ date / chip_id / scanType / run_id / stream_id / 'sorter_output'
    sorting = ss.Kilosort2_5Sorter._get_result_from_folder(sorting_dir)
    recon_dir = Path(temp_dir) /'reconstructions'/ date / chip_id / scanType / run_id / stream_id  
    unit_ids = sorting.get_unit_ids()
    for unit_id in unit_ids:

        ## debug
        if unit_id != 203: continue
        ##
                
        merged_template_dir = Path(templates_dir) / date / chip_id / scanType / run_id / stream_id / 'templates' / f'{unit_id}.npy'
        merged_channel_dir = Path(templates_dir) / date / chip_id / scanType / run_id / stream_id / 'templates' / f'{unit_id}_channels.npy'
        try:              

            merged_template = np.load(merged_template_dir)
            merged_channel_loc = np.load(merged_channel_dir)

            #merged_template = np.load('/mnt/disk15tb/adam/git_workspace/MEA_Analysis/AxonReconPipeline/data/temp_data/templates/240118/M06844/AxonTracking/000026/well001/partial_templates/seg3_unit17.npy')
            #merged_channel_loc = np.load('/mnt/disk15tb/adam/git_workspace/MEA_Analysis/AxonReconPipeline/data/temp_data/templates/240118/M06844/AxonTracking/000026/well001/partial_templates/seg3_unit17_channels.npy')
        except:
            #print(f'No merged template found for {unit_id}')
            continue

        template = merged_template[-1]
        # Iterate over the merged templates
        #Initialize list to store gtr objects
        gtrs = {}
        #for i, template in enumerate(templates):
            
        # Extract the sampling frequency for the current unit
        #hardcoded for now, fix later
        fs = 10000      

        transformed_template = template.T
        #transformed_template = np.array(template[0]).T
        #transformed_template = np.array(template)
        #trans_loc = merged_channel_loc[0]
        trans_loc = merged_channel_loc[-1]
        plot_dir = Path(recon_dir) / f'unit_{unit_id}' / 'axon_tracking_plots'
        key = 0
        
        # Generate amplitude_map figure
        try:
            fig_amp = plt.figure(figsize=(10, 5))
            ax_amp = fig_amp.add_subplot(111)
            #set_trace()
            ax_amp = plot_amplitude_map(transformed_template, 
                                        trans_loc, 
                                        log=False, 
                                        ax=ax_amp, 
                                        cmap="PRGn", 
                                        colorbar=True, 
                                        colorbar_orientation="horizontal")
            ax_amp.set_title(f"Amplitude", fontsize=20)

            #Save figure
            fig_name = "amplitude_map.png"
            fig_path = plot_dir / fig_name
            if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)
            plt.savefig(fig_path, dpi=300, bbox_inches='tight')
            plt.close()

            #accept_units_amp_map.append(key)
        except Exception as e:
            #reject_units_amp_map.append(key)
            print(f"unit {key} failed to generate amplitude_map, error: {e}")

        # Generate peak latency map figure
        try:
            fig_peaks = plt.figure(figsize=(10, 5))
            ax_peaks = fig_peaks.add_subplot(111)
            ax_peaks = plot_peak_latency_map(transformed_template, merged_channel_loc[i], fs=fs, log=False, ax=ax_peaks, colorbar=True, colorbar_orientation="horizontal")
            ax_peaks.set_title(f"Peak latency", fontsize=20)

            #Save figure
            fig_name = "peak_latency_map.png"
            fig_path = plot_dir / fig_name
            plt.savefig(fig_path, dpi=300, bbox_inches='tight')
            plt.close()

            #accept_units_peak_latency_map.append(key)
        except Exception as e:
            #reject_units_peak_latency_map.append(key)
            print(f"unit {key} failed to generate peak_latency_map, error: {e}")

        # Select channels and generate selected channels figures
        def plot_selected_channels(fig_title, selected_channels, locations, filename):
            fig = plt.figure(figsize=(10, 5))
            ax = fig.add_subplot(111)
            ax.set_title(fig_title, fontsize=20)
            ax.plot(locations[:, 0], locations[:, 1], marker=".", color="grey", ls="", alpha=0.2)
            ax.plot(locations[selected_channels, 0], 
                    locations[selected_channels, 1], marker=".", color="k", ls="", alpha=0.5)
            ax.axis("off")
            fig_name = filename
            fig_path = plot_dir / fig_name
            plt.savefig(fig_path, dpi=300, bbox_inches='tight')
            plt.close()
        try:
            # Build axon tracking object
            verbose_bool = True
            gtr0 = GraphAxonTracking(transformed_template, merged_channel_loc[i], fs, verbose=verbose_bool, **params)
            

            # Select channels
            gtr0.select_channels()                

            # Use the function for each block
            plot_selected_channels(f"Selected after detection threshold: {gtr0._detect_threshold} ms", 
                                np.array(list(gtr0._selected_channels_detect)), gtr0.locations, "detection_threshold.png")
            plot_selected_channels(f"Selected after peak std threshold: {gtr0._peak_std_threhsold} ms", 
                                    np.array(list(gtr0._selected_channels_peakstd)), gtr0.locations, "peak_std_threshold.png")
            plot_selected_channels(f"Selected after init delay threshold: {gtr0._init_delay} ms", 
                                    np.array(list(gtr0._selected_channels_init)), gtr0.locations, "init_delay_threshold.png")
            plot_selected_channels("Selected after all thresholds", gtr0.selected_channels, gtr0.locations, "all_thresholds.png")
            #accept_units_selected_channels.append(key)
        except Exception as e:
            #reject_units_selected_channels.append(key)
            print(f"unit {key} failed to generate selected_channels, error: {e}")              

        # Generate axon_reconstruction_heuristics figure
        try:
            gtr0.build_graph()
            gtr0.find_paths()
            gtr0._verbose = 1

            fig_graph = plt.figure(figsize=(10, 7))
            fig_graph = gtr0.plot_graph(node_search_labels=False, fig=fig_graph, cmap_nodes="viridis", cmap_edges="YlGn")

            fig_name = "axon_reconstruction_heuristics.png"
            fig_path = plot_dir / fig_name
            plt.savefig(fig_path, dpi=300, bbox_inches='tight')
            plt.close()

            #accept_units_axon_reconstruction_heuristics.append(key)
        except Exception as e:
            #reject_units_axon_reconstruction_heuristics.append(key)
            print(f"unit {key} failed to generate axon_reconstruction_heuristics, error: {e}")

        # Generate axon_reconstruction figure
        try:
            fig_graph = plt.figure(figsize=(12, 7))
            ax_raw = fig_graph.add_subplot(111)                
            axpaths_raw = ax_raw
            axpaths_raw = gtr0.plot_raw_branches(cmap="tab20", plot_bp=True, plot_neighbors=True, 
                                                plot_full_template=True, ax=axpaths_raw)
            axpaths_raw.legend(fontsize=6)
            axpaths_raw.set_title("Raw Branches")

            fig_name = "axon_reconstruction.png"
            fig_path = plot_dir / fig_name
            plt.savefig(fig_path, dpi=600, bbox_inches='tight')
            plt.close()
            # gtrs[key] = compute_graph_propagation_velocity(template.T, merged_channel_loc[i], fs, 
            #                                         verbose=False, **params)

            #accept_units_axon_reconstruction.append(key)
            #flag = True
            #break
        except Exception as e:
            #reject_units_axon_reconstruction.append(key)
            print(f"unit {key} failed to generate axon_reconstruction, error: {e}")
            
        # Generate axon_reconstruction_velocities figure
        try:
            fvel, ax_vel = plt.subplots(figsize=(9, 15))
            ax_vel.axis("off")
            fvel = gtr0.plot_velocities(cmap="tab20", plot_outliers=True, fig=fvel, markersize=12, markersize_out=18, fs=20)
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)

            fig_name = "axon_reconstruction_velocities.png"
            fig_path = plot_dir / fig_name
            plt.savefig(fig_path, dpi=600, bbox_inches='tight')
            plt.close()

            #accept_units_axon_reconstruction_velocities.append(key)
        except Exception as e:
            #reject_units_axon_reconstruction_velocities.append(key)
            print(f"unit {key} failed to generate axon_reconstruction_velocities, error: {e}")