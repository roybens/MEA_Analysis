''' This file contains plotting and analysis functions that use the axon_velocity package in main pipeline functions. '''

import os
import h5py
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

from MEAProcessingLibrary import mea_processing_library as MPL
from AxonReconPipeline.axon_velocity.axon_velocity import GraphAxonTracking
import AxonReconPipeline.axon_velocity.axon_velocity as av
import AxonReconPipeline.axon_velocity.axon_velocity.plotting as av_plotting


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

def transform_data(merged_template, merged_channel_loc, merged_template_filled=None, merged_channel_filled_loc=None):
    transformed_template = merged_template.T  # array with time samples and amplitude values (voltage or v/s)
    if merged_template_filled is not None: 
        transformed_template_filled = merged_template_filled.T  # array with time samples and amplitude values (voltage or v/s)
    else: 
        transformed_template_filled = None
    trans_loc = np.array([[loc[0], loc[1] * -1] for loc in merged_channel_loc])
    if merged_channel_filled_loc is not None: 
        trans_loc_filled = np.array([[loc[0], loc[1] * -1] for loc in merged_channel_filled_loc])
    else: 
        trans_loc_filled = None
    return transformed_template, transformed_template_filled, trans_loc, trans_loc_filled

def create_plot_dir(recon_dir, unit_id):
    plot_dir = Path(recon_dir) / f'unit_{unit_id}' / 'axon_tracking_plots'
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    return plot_dir

def save_figure(fig, fig_path, dpi=300):
    plt.savefig(fig_path, dpi=dpi, bbox_inches='tight')
    plt.close()

def generate_amplitude_map(template, locations, plot_dir, title, fresh_plots=False, cmap='viridis', log=False):
    fig_amp = plt.figure(figsize=(10, 5))
    ax_amp = fig_amp.add_subplot(111)
    fig_path = plot_dir / f"{title}_amplitude_map.png"
    if os.path.exists(fig_path) and not fresh_plots:
        return
    ax_amp = av.plotting.plot_amplitude_map(template, locations, cmap=cmap, log=log,
                                            elec_size=8, alpha=0.9, ax=None, colorbar=True,
                                            colorbar_orientation="vertical", colorbar_shrink=0.5,
                                            plot_image=True)
    ax_amp.set_title(f"Amplitude {title}", fontsize=20)
    save_figure(fig_amp, fig_path)

def generate_peak_latency_map(template, locations, plot_dir, title, fresh_plots=False, cmap='viridis', fs=10000, log=False):
    fig_peaks = plt.figure(figsize=(10, 5))
    ax_peaks = fig_peaks.add_subplot(111)
    fig_path = plot_dir / f"{title}_peak_latency_map.png"
    if os.path.exists(fig_path) and not fresh_plots:
        return
    ax_peaks = av.plotting.plot_peak_latency_map(template, locations, fs, cmap=cmap, log=log,
                                                 elec_size=8, alpha=0.9, ax=None, colorbar=True,
                                                 colorbar_orientation="vertical", colorbar_shrink=0.5,
                                                 plot_image=True)
    ax_peaks.set_title(f"Peak latency {title}", fontsize=20)
    save_figure(fig_peaks, fig_path)

def plot_selected_channels(gtr, plot_dir, suffix="", fresh_plots=False):
    plot_selected_channels_helper(f"Selected after detection threshold: {gtr._detect_threshold} uV", 
                                  np.array(list(gtr._selected_channels_detect)), gtr.locations, 
                                  f"detection_threshold{suffix}.png", plot_dir, fresh_plots=fresh_plots)
    plot_selected_channels_helper(f"Selected after kurtosis threshold: {gtr._kurt_threshold}", 
                                  np.array(list(gtr._selected_channels_kurt)), gtr.locations, 
                                  f"kurt_threshold{suffix}.png", plot_dir, fresh_plots=fresh_plots)
    plot_selected_channels_helper(f"Selected after peak std threshold: {gtr._peak_std_threhsold} ms", 
                                  np.array(list(gtr._selected_channels_peakstd)), gtr.locations, 
                                  f"peak_std_threshold{suffix}.png", plot_dir, fresh_plots=fresh_plots)
    plot_selected_channels_helper(f"Selected after init delay threshold: {gtr._init_delay} ms", 
                                  np.array(list(gtr._selected_channels_init)), gtr.locations, 
                                  f"init_delay_threshold{suffix}.png", plot_dir, fresh_plots=fresh_plots)
    plot_selected_channels_helper("Selected after all thresholds", gtr.selected_channels, gtr.locations, 
                                  f"all_thresholds{suffix}.png", plot_dir, fresh_plots=fresh_plots)

def plot_selected_channels_helper(fig_title, selected_channels, locations, filename, plot_dir, fresh_plots=False):
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111)
    ax.set_title(fig_title, fontsize=20)
    fig_path = plot_dir / filename
    if os.path.exists(fig_path) and not fresh_plots:
        return
    ax.plot(locations[:, 0], locations[:, 1], marker=".", color="grey", ls="", alpha=0.2)
    ax.plot(locations[selected_channels, 0], locations[selected_channels, 1], marker=".", color="k", ls="", alpha=0.5)
    ax.axis("off")
    save_figure(fig, fig_path)

def plot_template_propagation(gtr, fig_dir, unit_id, title, fresh_plots=False, figsize=(8, 40), 
                              linewidth=1, markerfacecolor='r', markersize=5, dpi=300, 
                              max_overlapped_lines=1, color_marker='none', marker='.'):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    fig_path = fig_dir
    if os.path.exists(fig_path) and not fresh_plots:
        return
    template = gtr.template
    selected_channels = gtr.selected_channels
    locations = gtr.locations
    ax = av.plotting.plot_template_propagation(
        template, locations, selected_channels, sort_templates=True,
        color=None, color_marker=color_marker, ax=None, linewidth=linewidth, 
        markerfacecolor=markerfacecolor, markersize=markersize, 
        max_overlapped_lines=max_overlapped_lines, marker=marker
    )
    ax.set_title(f'Template Propagation for Unit {unit_id} {title}', fontsize=16)
    save_figure(fig, fig_path, dpi=dpi)

def generate_axon_analytics(gtr, units, branch_ids, velocities, path_lengths, r2s, extremums, init_chans, num_channels_included, 
                            channel_density, unit_id, transformed_template, trans_loc):
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
    num_channels_included.append(len(trans_loc))
    x_coords = trans_loc[:, 0]
    y_coords = trans_loc[:, 1]
    width = np.max(x_coords) - np.min(x_coords)
    height = np.max(y_coords) - np.min(y_coords)
    area = width * height
    channel_density_value = len(trans_loc) / area  # channels / um^2
    channel_density.append(channel_density_value)
    init_chans.append(gtr.init_channel)
    extremums.append(get_extremum(transformed_template, trans_loc))

def generate_axon_reconstruction_heuristics(gtr, plot_dir, unit_id, fresh_plots=False, suffix="", figsize=(10, 7)):
    gtr.build_graph()
    gtr.find_paths()
    gtr._verbose = 1
    fig_graph = plt.figure(figsize=figsize)
    fig_path = plot_dir / f"axon_reconstruction_heuristics{suffix}.png"
    if os.path.exists(fig_path) and not fresh_plots:
        return
    fig_graph = gtr.plot_graph(node_search_labels=False, fig=fig_graph, cmap_nodes="viridis", cmap_edges="YlGn")
    save_figure(fig_graph, fig_path)

def generate_axon_reconstruction_raw(gtr, plot_dir, unit_id, recon_dir, successful_recons, fresh_plots=False, suffix="", plot_full_template=False):
    fig_graph = plt.figure(figsize=(12, 7))
    ax_raw = fig_graph.add_subplot(111)
    axpaths_raw = ax_raw
    fig_path = plot_dir / f"axon_reconstruction_raw{suffix}.png"
    if os.path.exists(fig_path) and not fresh_plots:
        return
    axpaths_raw = gtr.plot_raw_branches(
        plot_full_template=plot_full_template, ax=axpaths_raw, cmap="rainbow",
        plot_labels=True, plot_bp=True, plot_neighbors=True
    )
    axpaths_raw.legend(fontsize=6)
    axpaths_raw.set_title("Raw Branches")
    plt.savefig(fig_path, dpi=600, bbox_inches='tight')
    plt.close()
    successful_recons[str(recon_dir)]["successful_units"][str(unit_id)] = {}

def generate_axon_reconstruction_clean(gtr, plot_dir, unit_id, recon_dir, successful_recons, fresh_plots=False, suffix="", plot_full_template=False):
    fig_graph = plt.figure(figsize=(12, 7))
    ax_raw = fig_graph.add_subplot(111)
    axpaths_raw = ax_raw
    fig_path = plot_dir / f"axon_reconstruction_clean{suffix}.png"
    if os.path.exists(fig_path) and not fresh_plots:
        return
    axpaths_raw = gtr.plot_clean_branches(
        plot_full_template=plot_full_template, ax=axpaths_raw, cmap="rainbow",
        plot_bp=True, branch_colors=None
    )
    axpaths_raw.legend(fontsize=6)
    axpaths_raw.set_title("Clean Branches")
    save_figure(fig_graph, fig_path)
    successful_recons[str(recon_dir)]["successful_units"][str(unit_id)] = {}

def generate_axon_reconstruction_velocities(gtr, plot_dir, unit_id, fresh_plots=False, suffix=""):
    fvel, ax_vel = plt.subplots(figsize=(9, 15))
    ax_vel.axis("off")
    fig_path = plot_dir / f"axon_reconstruction_velocities{suffix}.png"
    if os.path.exists(fig_path) and not fresh_plots:
        return
    fvel = gtr.plot_velocities(cmap="tab20", plot_outliers=True, fig=fvel, markersize=12, markersize_out=18, fs=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    save_figure(fvel, fig_path)

def save_axon_analytics(stream_id, units, extremums, branch_ids, velocities, path_lengths, r2s, num_channels_included, 
                        channel_density, init_chans, recon_dir, suffix=""):
    df_mea1k = pd.DataFrame({"unit_ids": units, "unit location": extremums, "branch_id": branch_ids, 
                             "velocity": velocities, "length": path_lengths, "r2": r2s, 
                             "num_channels_included": num_channels_included, "channel_density": channel_density, 
                             "init_chan": init_chans})
    recon_dir_parent = os.path.dirname(recon_dir)
    if not os.path.exists(recon_dir_parent):
        os.makedirs(recon_dir_parent)
    df_mea1k.to_csv(Path(recon_dir_parent) / f"{stream_id}_axon_analytics{suffix}.csv", index=False)

def plot_branch_neurites_wrapper(data, save_path, **kwargs):
    """
    Wrapper for axon_velocity.plotting.plot_branch_neurites.

    Parameters:
    data : Data required for plotting branch neurites.
    save_path : Path where the plot will be saved.
    **kwargs : Additional arguments for plot customization.
    """
    av_plotting.plot_branch_neurites(data, save_path, **kwargs)

def play_template_map_wrapper(save_path, **kwargs):
    """
    Wrapper for axon_velocity.plotting.play_template_map.

    Parameters:
    data : Data required for playing the template map.
    save_path : Path where the plot will be saved.
    **kwargs : Additional arguments for plot customization.
    """
    log = kwargs.pop('log', False)
    template = kwargs.pop('template', None)
    locations = kwargs.pop('locations', None)
    gtr = kwargs.pop('gtr', None)
    elec_size = kwargs.pop('elec_size', 8)
    cmap = kwargs.pop('cmap', 'viridis')
    ax = kwargs.pop('ax', None)
    skip_frames = kwargs.pop('skip_frames', 1)
    interval = kwargs.pop('interval', 10)
    ani = av_plotting.play_template_map(
        template, locations, gtr=gtr, elec_size=elec_size, cmap=cmap, 
        log=log, ax=ax, skip_frames=skip_frames, interval=interval
    )
    ani.save(save_path, writer='ffmpeg')

def plot_template_wrapper(save_path, title, fresh_plots=False, **kwargs):
    """
    Wrapper for axon_velocity.plotting.plot_template.

    Parameters:
    data : Data required for plotting the template.
    save_path : Path where the plot will be saved.
    **kwargs : Additional arguments for plot customization.
    """
    if os.path.exists(save_path) and not fresh_plots:
        return
    fig = plt.figure(figsize=(100, 200))
    ax = fig.add_subplot(111)
    ax = av_plotting.plot_template(**kwargs)
    ax.set_title(title, fontsize=20)
    if save_path is not None: 
        save_figure(ax, save_path, dpi=2400)

def plot_axon_summary_wrapper(save_path, title, fresh_plots, **kwargs):
    """
    Wrapper for axon_velocity.plotting.plot_axon_summary.

    Parameters:
    data : Data required for plotting the axon summary.
    save_path : Path where the plot will be saved.
    **kwargs : Additional arguments for plot customization.
    """
    if os.path.exists(save_path) and not fresh_plots:
        return
    fig, axes = av_plotting.plot_axon_summary(**kwargs)
    fig.set_title(title, fontsize=20)
    if save_path is not None: 
        save_figure(fig, save_path, dpi=600)