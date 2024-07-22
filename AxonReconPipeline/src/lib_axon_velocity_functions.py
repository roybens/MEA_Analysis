''' This file contains plotting and analysis functions that use the axon_velocity package in main pipeline functions. '''

'''imports'''
import os
import h5py
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use('Agg')  # Use the Anti-Grain Geometry (Agg) backend
import matplotlib.pyplot as plt

# Local imports
# import sys
# sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from MEAProcessingLibrary import mea_processing_library as MPL

# Axon trace algorithm imports
#import MEAutility as mu
from AxonReconPipeline.axon_velocity import axon_velocity as av  # Import submodule
#from AxonReconPipeline.axon_velocity.axon_velocity.models import load_cell
from AxonReconPipeline.axon_velocity.axon_velocity.evaluation import *

'''functions'''
def get_extremum(template, locations):
    """
    Get the extremum of the template.
    """
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
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    return plot_dir

'''plotting and analysis functions'''
def save_figure(fig, fig_path):
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()

def generate_amplitude_map(transformed_template_filled, trans_loc_filled, plot_dir, title, fresh_plots=False):
    fig_amp = plt.figure(figsize=(10, 5))
    ax_amp = fig_amp.add_subplot(111)
    fig_path = plot_dir / f"{title}_amplitude_map.png"
    if os.path.exists(fig_path) and fresh_plots == False: return
    ax_amp = av.plot_amplitude_map(transformed_template_filled, 
                                   trans_loc_filled, 
                                   log=False, 
                                   ax=ax_amp, 
                                   cmap="PRGn", 
                                   colorbar=True, 
                                   colorbar_orientation="horizontal")
    ax_amp.set_title(f"Amplitude {title}", fontsize=20)
    save_figure(fig_amp, fig_path)

def generate_peak_latency_map(transformed_template_filled, trans_loc_filled, plot_dir, title, fresh_plots=False):
    fig_peaks = plt.figure(figsize=(10, 5))
    ax_peaks = fig_peaks.add_subplot(111)
    fig_path = plot_dir / f"{title}_peak_latency_map.png"
    if os.path.exists(fig_path) and fresh_plots == False: return
    ax_peaks = av.plot_peak_latency_map(transformed_template_filled, 
                                        trans_loc_filled, fs=10000, log=False, 
                                        ax=ax_peaks, colorbar=True, 
                                        colorbar_orientation="horizontal")
    ax_peaks.set_title(f"Peak latency {title}", fontsize=20)
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

def plot_template_propagation_wrapper(template, locations, selected_channels, plot_dir, unit_id, title, fresh_plots=False):
    fig = plt.figure(figsize=(8, 40))
    ax = fig.add_subplot(111)
    fig_path = plot_dir / f"{title}_template_propagation_unit_{unit_id}.png"
    if os.path.exists(fig_path) and fresh_plots == False: return

    ax = av.plot_template_propagation(template, locations, selected_channels, ax=ax, sort_templates=True)
    ax.set_title(f'Template Propagation for Unit {unit_id} {title}', fontsize=16)
    save_figure(fig, fig_path)

def select_and_plot_channels(gtr0, plot_dir, suffix=""):
    gtr0.select_channels()
    plot_selected_channels(f"Selected after detection threshold: {gtr0._detect_threshold} ms", 
                        np.array(list(gtr0._selected_channels_detect)), gtr0.locations, f"detection_threshold{suffix}.png", plot_dir)
    plot_selected_channels(f"Selected after detection threshold: {gtr0._kurt_threshold} ms", 
                        np.array(list(gtr0._selected_channels_kurt)), gtr0.locations, f"kurt_threshold{suffix}.png", plot_dir)
    plot_selected_channels(f"Selected after peak std threshold: {gtr0._peak_std_threhsold} ms", 
                            np.array(list(gtr0._selected_channels_peakstd)), gtr0.locations, f"peak_std_threshold{suffix}.png", plot_dir)
    plot_selected_channels(f"Selected after init delay threshold: {gtr0._init_delay} ms", 
                            np.array(list(gtr0._selected_channels_init)), gtr0.locations, f"init_delay_threshold{suffix}.png", plot_dir)
    plot_selected_channels("Selected after all thresholds", gtr0.selected_channels, gtr0.locations, f"all_thresholds{suffix}.png", plot_dir)

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

def generate_axon_analytics_dvdt(gtr, units, branch_ids, velocities, path_lengths, r2s, extremums, unit_id, transformed_template, trans_loc):
    units.append(unit_id)
    for i, br in enumerate(gtr.branches):
        path = br["channels"]
        velocity = br["velocity"]
        r2 = br["r2"]
        length = gtr.compute_path_length(path)
        branch_ids.append(i)
        velocities.append(velocity)
        path_lengths.append(length)
        r2s.append(r2)
    extremums.append(get_extremum(transformed_template, trans_loc))

def generate_axon_reconstruction_heuristics(gtr0, plot_dir, unit_id, fresh_plots=False, suffix=""):
    gtr0.build_graph()
    gtr0.find_paths()
    gtr0._verbose = 1
    fig_graph = plt.figure(figsize=(10, 7))
    fig_path = plot_dir / f"axon_reconstruction_heuristics{suffix}.png"
    if os.path.exists(fig_path) and fresh_plots == False: return
    fig_graph = gtr0.plot_graph(node_search_labels=False, fig=fig_graph, cmap_nodes="viridis", cmap_edges="YlGn")
    save_figure(fig_graph, fig_path)

def generate_axon_reconstruction_raw(gtr0, plot_dir, unit_id, recon_dir, successful_recons, fresh_plots=False, suffix=""):
    fig_graph = plt.figure(figsize=(12, 7))
    ax_raw = fig_graph.add_subplot(111)
    axpaths_raw = ax_raw
    fig_path = plot_dir / f"axon_reconstruction_raw{suffix}.png"
    if os.path.exists(fig_path) and fresh_plots == False: return
    axpaths_raw = gtr0.plot_raw_branches(cmap="tab20", plot_bp=True, plot_neighbors=True, plot_full_template=True, ax=axpaths_raw)
    axpaths_raw.legend(fontsize=6)
    axpaths_raw.set_title("Raw Branches")
    plt.savefig(fig_path, dpi=600, bbox_inches='tight')
    plt.close()
    successful_recons[str(recon_dir)]["successful_units"][str(unit_id)] = {}

def generate_axon_reconstruction_clean(gtr0, plot_dir, unit_id, recon_dir, successful_recons, fresh_plots=False, suffix=""):
    fig_graph = plt.figure(figsize=(12, 7))
    ax_raw = fig_graph.add_subplot(111)
    axpaths_raw = ax_raw
    fig_path = plot_dir / f"axon_reconstruction_clean{suffix}.png"
    if os.path.exists(fig_path) and fresh_plots == False: return
    axpaths_raw = gtr0.plot_clean_branches(cmap="tab20", plot_bp=True, plot_full_template=True, ax=axpaths_raw)
    axpaths_raw.legend(fontsize=6)
    axpaths_raw.set_title("Clean Branches")
    save_figure(fig_graph, fig_path)
    successful_recons[str(recon_dir)]["successful_units"][str(unit_id)] = {}

def generate_axon_reconstruction_velocities(gtr0, plot_dir, unit_id, fresh_plots=False, suffix=""):
    fvel, ax_vel = plt.subplots(figsize=(9, 15))
    ax_vel.axis("off")
    fig_path = plot_dir / f"axon_reconstruction_velocities{suffix}.png"
    if os.path.exists(fig_path) and fresh_plots == False: return
    fvel = gtr0.plot_velocities(cmap="tab20", plot_outliers=True, fig=fvel, markersize=12, markersize_out=18, fs=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    save_figure(fvel, fig_path)

def save_axon_analytics(stream_id, units, extremums, branch_ids, velocities, path_lengths, r2s, num_channels_included, channel_density, recon_dir,  suffix=""):
    df_mea1k = pd.DataFrame({"unit_ids": units, "unit location": extremums, "branch_id": branch_ids, "velocity": velocities, "length": path_lengths, "r2": r2s, "num_channels_included": num_channels_included, "channel_density": channel_density})
    recon_dir_parent = os.path.dirname(recon_dir)
    if not os.path.exists(recon_dir_parent): os.makedirs(recon_dir_parent)
    df_mea1k.to_csv(Path(recon_dir_parent) / f"{stream_id}_axon_analytics{suffix}.csv", index=False)

def plot_voltage_signals(template, locations, selected_channels, plot_dir, unit_id, fresh_plots=False):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111)
    fig_path = plot_dir / f"voltage_signals_unit_{unit_id}.png"
    if os.path.exists(fig_path) and fresh_plots == False: return

    template_selected = template[selected_channels]
    locs = locations[selected_channels]
    
    peaks = np.argmin(template_selected, 1)
    sorted_idxs = np.argsort(peaks)
    template_sorted = template_selected[sorted_idxs]
    locs_sorted = locs[sorted_idxs]
    
    dist_peaks_sorted = np.array([np.linalg.norm(loc - locs_sorted[0]) for loc in locs_sorted])
    dist_peaks_sorted /= np.max(dist_peaks_sorted)
    dist_peaks_sorted *= len(template_sorted)
    
    ptp_glob = np.max(np.ptp(template_sorted, 1))
    for i, temp in enumerate(template_sorted):
        temp_shifted = temp + i * 1.5 * ptp_glob
        min_t = np.min(temp_shifted)
        min_idx = np.argmin(temp_shifted)
        ax.plot(temp_shifted, color='C0')
        ax.plot(min_idx, min_t, marker='o', color='r')
    
    ax.axis('off')
    ax.set_title(f'Voltage Signals for Unit {unit_id}', fontsize=16)
    save_figure(fig, fig_path)

def plot_all_traces(template, locations, plot_dir, unit_id, fresh_plots=False):
    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(111)
    fig_path = plot_dir / f"all_traces_unit_{unit_id}.png"
    if os.path.exists(fig_path) and fresh_plots == False: return

    ptp_glob = np.max(np.ptp(template, 1))
    for i, temp in enumerate(template):
        temp_shifted = temp + i * 1.5 * ptp_glob
        min_t = np.min(temp_shifted)
        min_idx = np.argmin(temp_shifted)
        ax.plot(temp_shifted, color='C0')
        ax.plot(min_idx, min_t, marker='o', color='r')
    
    ax.axis('off')
    ax.set_title(f'All Voltage Traces for Unit {unit_id}', fontsize=16)
    save_figure(fig, fig_path)

def plot_template_propagation_wrapper(template, locations, selected_channels, plot_dir, unit_id, title, fresh_plots=False):
    fig = plt.figure(figsize=(8, 40))
    ax = fig.add_subplot(111)
    fig_path = plot_dir / f"{title}_template_propagation_unit_{unit_id}.png"
    if os.path.exists(fig_path) and fresh_plots == False: return

    ax = av.plot_template_propagation(template, locations, selected_channels, ax=ax, sort_templates=True)
    ax.set_title(f'Template Propagation for Unit {unit_id} {title}', fontsize=16)
    save_figure(fig, fig_path)