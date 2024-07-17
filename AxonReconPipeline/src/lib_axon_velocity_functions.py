''' This file contains plotting and analysis functions that use the axon_velocity package in main pipeline functions. '''

'''imports'''
import sys
import os
import h5py
import spikeinterface.sorters as ss
import shutil
import numpy as np
from tqdm import tqdm
from pathlib import Path
import matplotlib
matplotlib.use('Agg')  # Use the Anti-Grain Geometry (Agg) backend
import matplotlib.pyplot as plt
import json
import pandas as pd

#parallel processing
import concurrent.futures
import threading
lock = threading.Lock()
from multiprocessing import Pool, Manager

# Axon trace algo imports
import MEAutility as mu
import axon_velocity as av
from axon_velocity import *
from axon_velocity.models import load_cell
from axon_velocity.evaluation import *
from probeinterface import plotting as plotting_probe

''' Local imports '''
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from MEAProcessingLibrary import mea_processing_library as MPL

'''functions'''
def get_extremum(template, locations):
    """
    Get the extremum of the template.
    """
    # Get the extremum index along the first dim of the template, aligning with channel locations
    extremum_idx = np.unravel_index(np.argmax(np.abs(template)), template.shape)[0]
    extremum = locations[extremum_idx]
    return tuple(extremum)

def extract_recording_info(continuous_h5_dir):
    recording_details = MPL.extract_recording_details(continuous_h5_dir['h5_file_path'])
    date = recording_details[0]['date']
    chip_id = recording_details[0]['chipID']
    scanType = recording_details[0]['scanType']
    run_id = recording_details[0]['runID']
    return date, chip_id, scanType, run_id

def get_stream_ids(h5_file_path):
    h5 = h5py.File(h5_file_path, 'r')
    return list(h5['wells'].keys())

def transform_data(merged_template, merged_channel_loc, merged_template_filled, merged_channel_filled_loc):
    transformed_template = merged_template.T
    transformed_template_filled = merged_template_filled.T
    trans_loc = np.array([[loc[0], loc[1]*-1] for loc in merged_channel_loc])
    trans_loc_filled = np.array([[loc[0], loc[1]*-1] for loc in merged_channel_filled_loc])
    return transformed_template, transformed_template_filled, trans_loc, trans_loc_filled

def create_plot_dir(recon_dir, unit_id):
    plot_dir = Path(recon_dir) / f'unit_{unit_id}' / 'axon_tracking_plots'
    # if os.path.exists(plot_dir):
    #     shutil.rmtree(plot_dir)
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    return plot_dir

'''ploltting and analysis functions'''
def save_figure(fig, fig_path):
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()

def generate_amplitude_map(transformed_template_filled, trans_loc_filled, plot_dir, fresh_plots=False):
    fig_amp = plt.figure(figsize=(10, 5))
    ax_amp = fig_amp.add_subplot(111)
    fig_path = plot_dir / "amplitude_map.png"
    if os.path.exists(fig_path) and fresh_plots == False: return
    ax_amp = plot_amplitude_map(transformed_template_filled, 
                                trans_loc_filled, 
                                log=False, 
                                ax=ax_amp, 
                                cmap="PRGn", 
                                colorbar=True, 
                                colorbar_orientation="horizontal")
    ax_amp.set_title(f"Amplitude", fontsize=20)
    save_figure(fig_amp, fig_path)

def generate_peak_latency_map(transformed_template_filled, trans_loc_filled, plot_dir, fresh_plots=False):
    fig_peaks = plt.figure(figsize=(10, 5))
    ax_peaks = fig_peaks.add_subplot(111)
    fig_path = plot_dir / "peak_latency_map.png"
    if os.path.exists(fig_path) and fresh_plots == False: return
    ax_peaks = plot_peak_latency_map(transformed_template_filled, 
                                    trans_loc_filled, fs=10000, log=False, 
                                    ax=ax_peaks, colorbar=True, 
                                    colorbar_orientation="horizontal")
    ax_peaks.set_title(f"Peak latency", fontsize=20)
    save_figure(fig_peaks, fig_path)

def plot_selected_channels(fig_title, selected_channels, locations, filename, plot_dir, fresh_plots=False):
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111)
    ax.set_title(fig_title, fontsize=20)
    fig_path = plot_dir / filename
    if os.path.exists(fig_path) and fresh_plots == False: return
    ax.plot(locations[:, 0], locations[:, 1], marker=".", color="grey", ls="", alpha=0.2)
    ax.plot(locations[selected_channels, 0], 
            locations[selected_channels, 1], marker=".", color="k", ls="", alpha=0.5)
    ax.axis("off")
    save_figure(fig, fig_path)

def select_and_plot_channels(gtr0, plot_dir):
    gtr0.select_channels()
    plot_selected_channels(f"Selected after detection threshold: {gtr0._detect_threshold} ms", 
                        np.array(list(gtr0._selected_channels_detect)), gtr0.locations, "detection_threshold.png", plot_dir)
    plot_selected_channels(f"Selected after detection threshold: {gtr0._kurt_threshold} ms", 
                        np.array(list(gtr0._selected_channels_kurt)), gtr0.locations, "kurt_threshold.png", plot_dir)
    plot_selected_channels(f"Selected after peak std threshold: {gtr0._peak_std_threhsold} ms", 
                            np.array(list(gtr0._selected_channels_peakstd)), gtr0.locations, "peak_std_threshold.png", plot_dir)
    plot_selected_channels(f"Selected after init delay threshold: {gtr0._init_delay} ms", 
                            np.array(list(gtr0._selected_channels_init)), gtr0.locations, "init_delay_threshold.png", plot_dir)
    plot_selected_channels("Selected after all thresholds", gtr0.selected_channels, gtr0.locations, "all_thresholds.png", plot_dir)

def generate_axon_analytics(gtr0, units, branch_ids, velocities, path_lengths, r2s, extremums, unit_id, transformed_template, trans_loc):
    units.append(unit_id)
    for i, br in enumerate(gtr0.branches):
        path = br["channels"]
        velocity = br["velocity"]
        r2 = br["r2"]
        length = gtr0.compute_path_length(path)
        branch_ids.append(i)
        velocities.append(velocity)
        path_lengths.append(length)
        r2s.append(r2)
    extremums.append(get_extremum(transformed_template, trans_loc))

def generate_axon_reconstruction_heuristics(gtr0, plot_dir, unit_id, fresh_plots=False):
    gtr0.build_graph()
    gtr0.find_paths()
    gtr0._verbose = 1
    fig_graph = plt.figure(figsize=(10, 7))
    fig_path = plot_dir / "axon_reconstruction_heuristics.png"
    if os.path.exists(fig_path) and fresh_plots == False: return
    fig_graph = gtr0.plot_graph(node_search_labels=False, fig=fig_graph, cmap_nodes="viridis", cmap_edges="YlGn")
    save_figure(fig_graph, fig_path)
    #logger.info(f"unit {unit_id} successfully generated axon_reconstruction_heuristics")

def generate_axon_reconstruction_raw(gtr0, plot_dir, unit_id, recon_dir, successful_recons, fresh_plots=False):
    fig_graph = plt.figure(figsize=(12, 7))
    ax_raw = fig_graph.add_subplot(111)
    axpaths_raw = ax_raw
    fig_path = plot_dir / "axon_reconstruction_raw.png"
    if os.path.exists(fig_path) and fresh_plots == False: return
    axpaths_raw = gtr0.plot_raw_branches(cmap="tab20", plot_bp=True, plot_neighbors=True, plot_full_template=True, ax=axpaths_raw)
    axpaths_raw.legend(fontsize=6)
    axpaths_raw.set_title("Raw Branches")
    plt.savefig(fig_path, dpi=600, bbox_inches='tight')
    plt.close()
    successful_recons[str(recon_dir)]["successful_units"][str(unit_id)] = {}
    #logger.info(f"unit {unit_id} successfully generated raw axon_reconstruction")

def generate_axon_reconstruction_clean(gtr0, plot_dir, unit_id, recon_dir, successful_recons, fresh_plots=False):
    fig_graph = plt.figure(figsize=(12, 7))
    ax_raw = fig_graph.add_subplot(111)
    axpaths_raw = ax_raw
    fig_path = plot_dir / "axon_reconstruction_clean.png"
    if os.path.exists(fig_path) and fresh_plots == False: return
    axpaths_raw = gtr0.plot_clean_branches(cmap="tab20", plot_bp=True, plot_full_template=True, ax=axpaths_raw)
    axpaths_raw.legend(fontsize=6)
    axpaths_raw.set_title("Clean Branches")
    save_figure(fig_graph, fig_path)
    successful_recons[str(recon_dir)]["successful_units"][str(unit_id)] = {}
    #logger.info(f"unit {unit_id} successfully generated clean axon_reconstruction")

def generate_axon_reconstruction_velocities(gtr0, plot_dir, unit_id, fresh_plots=False):
    fvel, ax_vel = plt.subplots(figsize=(9, 15))
    ax_vel.axis("off")
    fig_path = plot_dir / "axon_reconstruction_velocities.png"
    if os.path.exists(fig_path) and fresh_plots == False: return
    fvel = gtr0.plot_velocities(cmap="tab20", plot_outliers=True, fig=fvel, markersize=12, markersize_out=18, fs=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    save_figure(fvel, fig_path)
    #logger.info(f"unit {unit_id} successfully generated axon_reconstruction_velocities")

def save_axon_analytics(stream_id, units, extremums, branch_ids, velocities, path_lengths, r2s, num_channels_included, channel_density, recon_dir):
    df_mea1k = pd.DataFrame({"unit_ids": units, "unit location": extremums, "branch_id": branch_ids, "velocity": velocities, "length": path_lengths, "r2": r2s, "num_channels_included": num_channels_included, "channel_density": channel_density})
    recon_dir_parent = os.path.dirname(recon_dir)
    if not os.path.exists(recon_dir_parent): os.makedirs(recon_dir_parent)
    df_mea1k.to_csv(Path(recon_dir_parent) / f"{stream_id}_axon_analytics.csv", index=False)

def process_unit(unit_id, unit_templates, recon_dir, params, successful_recons, logger=None):
    if logger is not None: logger.info(f'Processing unit {unit_id}')
    else: print(f'Processing unit {unit_id}')
    #merged_template_dir, merged_channel_dir, merged_template_filled_dir, merged_channel_filled_dir = get_paths(templates_dir, date, chip_id, scanType, run_id, stream_id, unit_id)
    #merged_template, merged_channel_loc, merged_template_filled, merged_channel_filled_loc = load_data(merged_template_dir, merged_channel_dir, merged_template_filled_dir, merged_channel_filled_dir, unit_id)
    
    merged_template = unit_templates['merged_template']
    merged_channel_loc = unit_templates['merged_channel_loc']
    merged_template_filled = unit_templates['merged_template_filled']
    merged_channel_filled_loc = unit_templates['merged_channel_locs_filled']

    if merged_template is None: return None  # skip if no merged template found
    
    transformed_template, transformed_template_filled, trans_loc, trans_loc_filled = transform_data(merged_template, merged_channel_loc, merged_template_filled, merged_channel_filled_loc)
    plot_dir = create_plot_dir(recon_dir, unit_id)
    units, branch_ids, velocities, path_lengths, r2s, extremums, num_channels_included, channel_density = [], [], [], [], [], [], [], []

    # Assuming merged_channel_loc and merged_channel_filled_loc are numpy arrays of shape (n, 2) where n is the number of channels
    #frac_chans_included.append(len(merged_channel_loc) / len(merged_channel_filled_loc))
    num_channels_included.append(len(merged_channel_loc))

    # Calculate the area covered by the channels
    # Get the bounding box of the channels to estimate the area
    x_coords = merged_channel_loc[:, 0]
    y_coords = merged_channel_loc[:, 1]

    # Calculate the area of the bounding box
    width = np.max(x_coords) - np.min(x_coords)
    height = np.max(y_coords) - np.min(y_coords)
    area = width * height

    # Calculate the channel density as the number of channels per unit area
    channel_density.append(len(merged_channel_loc) / area) #channels / um^2
    
    try:
        #with lock: 
        generate_amplitude_map(transformed_template_filled, trans_loc_filled, plot_dir)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate amplitude_map, error: {e}")
        plt.close()
    
    try:
        #with lock: 
        generate_peak_latency_map(transformed_template_filled, trans_loc_filled, plot_dir)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate peak_latency_map, error: {e}")
        plt.close()
    
    try:
        #with lock: 
        gtr0 = GraphAxonTracking(transformed_template, trans_loc, 10000, **params)
        select_and_plot_channels(gtr0, plot_dir)
        gtr0 = compute_graph_propagation_velocity(transformed_template, trans_loc, 10000, **params)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to select and plot channels or compute_graph_propagation_velocity, error: {e}")
        plt.close()
        return None
    
    try:
        #with lock:
        generate_axon_analytics(gtr0, units, branch_ids, velocities, path_lengths, r2s, extremums, unit_id, transformed_template, trans_loc)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate axon_analytics, error: {e}")
        #plt.close()

    try:
        #with lock: 
        generate_axon_reconstruction_heuristics(gtr0, plot_dir, unit_id)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate axon_reconstruction_heuristics, error: {e}")
        plt.close() #close the figure, avoid memory leak
    
    try:
        #with lock:
        generate_axon_reconstruction_raw(gtr0, plot_dir, unit_id, recon_dir, successful_recons)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate raw axon_reconstruction, error: {e}")
        plt.close() #close the figure, avoid memory leak
    
    try:
        #with lock: 
        generate_axon_reconstruction_clean(gtr0, plot_dir, unit_id, recon_dir, successful_recons)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate clean axon_reconstruction, error: {e}")
        plt.close() #close the figure, avoid memory leak
    
    try:
        #with lock: 
        generate_axon_reconstruction_velocities(gtr0, plot_dir, unit_id)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate axon_reconstruction_velocities, error: {e}")
        plt.close() #close the figure, avoid memory leak
    
    return units, extremums, branch_ids, velocities, path_lengths, r2s, num_channels_included, channel_density

def analyze_and_reconstruct(templates, params=av.get_default_graph_velocity_params(), recon_dir = None, stream_select=None, n_jobs = 8, logger=None):
    for key, tmps in templates.items():        
        for stream_id, stream_templates in tmps['streams'].items():
            if stream_select is not None: 
                if stream_id != f'well{stream_select:03}': continue
            unit_templates = stream_templates['units']
            date = tmps['date']
            chip_id = tmps['chip_id']
            scanType = tmps['scanType']
            run_id = tmps['run_id']
            recon_dir = Path(recon_dir) / date / chip_id / scanType / run_id / stream_id
            successful_recons = {}
            unit_ids = [unit_id for unit_id in unit_templates.keys()]
            successful_recons[str(recon_dir)] = {}
            successful_recons[str(recon_dir)]["successful_units"] = {}
            all_units, all_branch_ids, all_velocities, all_path_lengths, all_r2s, all_extremums, all_num_channels_included, all_channel_densities = [], [], [], [], [], [], [], []
            logger.info(f'Processing {len(unit_ids)} units, with {n_jobs} workers')
            
            with concurrent.futures.ProcessPoolExecutor(max_workers=n_jobs) as executor:
                futures = [executor.submit(process_unit, unit_id, unit_templates, recon_dir, params, successful_recons, logger = logger) for unit_id, unit_templates in unit_templates.items()]
                for future in concurrent.futures.as_completed(futures):
                    try:
                        result = future.result()
                        if result:
                            units, extremums, branch_ids, velocities, path_lengths, r2s, channels_included, channel_density= result
                            all_units.append(units)
                            all_branch_ids.append(branch_ids)
                            all_velocities.append(velocities)
                            all_path_lengths.append(path_lengths)
                            all_r2s.append(r2s)
                            all_extremums.append(extremums)
                            all_num_channels_included.append(channels_included)
                            all_channel_densities.append(channel_density)
                    except Exception as exc:
                        logger.info(f'Unit generated an exception: {exc}')

            save_axon_analytics(stream_id, all_units, all_extremums, all_branch_ids, all_velocities, all_path_lengths, all_r2s, all_num_channels_included, all_channel_densities, recon_dir)
            
            # with open(f'{recon_dir}/successful_recons.json', 'w') as f:
            #     json.dump(successful_recons, f)