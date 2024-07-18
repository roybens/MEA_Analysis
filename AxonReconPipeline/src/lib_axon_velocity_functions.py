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
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from MEAProcessingLibrary import mea_processing_library as MPL
from AxonReconPipeline.axon_velocity.axon_velocity import GraphAxonTracking
import AxonReconPipeline.axon_velocity.axon_velocity as av


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
    transformed_template = merged_template.T #array with time samples and amplitude values (voltage or v/s)
    if merged_template_filled is not None: transformed_template_filled = merged_template_filled.T #array with time samples and amplitude values (voltage or v/s)
    else: transformed_template_filled = None
    trans_loc = np.array([[loc[0], loc[1]*-1] for loc in merged_channel_loc])
    if merged_channel_filled_loc is not None: trans_loc_filled = np.array([[loc[0], loc[1]*-1] for loc in merged_channel_filled_loc])
    else: trans_loc_filled = None
    return transformed_template, transformed_template_filled, trans_loc, trans_loc_filled

def create_plot_dir(recon_dir, unit_id):
    plot_dir = Path(recon_dir) / f'unit_{unit_id}' / 'axon_tracking_plots'
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    return plot_dir

<<<<<<< HEAD
def save_figure(fig, fig_path, dpi=300):
    plt.savefig(fig_path, dpi=dpi, bbox_inches='tight')
=======
'''plotting and analysis functions'''
def save_figure(fig, fig_path):
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
>>>>>>> 82357b0 (Pipeline is functional. Modificaitons to plot generated underway.)
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

def generate_peak_latency_map(template, locations, plot_dir, title, fresh_plots=False, cmap='viridis', fs = 10000, log=False):
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
    plot_selected_channels_helper(f"Selected after detection threshold: {gtr._detect_threshold} uV", np.array(list(gtr._selected_channels_detect)), gtr.locations, f"detection_threshold{suffix}.png", plot_dir, fresh_plots=fresh_plots)
    plot_selected_channels_helper(f"Selected after kurtosis threshold: {gtr._kurt_threshold}", np.array(list(gtr._selected_channels_kurt)), gtr.locations, f"kurt_threshold{suffix}.png", plot_dir, fresh_plots=fresh_plots)
    plot_selected_channels_helper(f"Selected after peak std threshold: {gtr._peak_std_threhsold} ms", np.array(list(gtr._selected_channels_peakstd)), gtr.locations, f"peak_std_threshold{suffix}.png", plot_dir, fresh_plots=fresh_plots)
    plot_selected_channels_helper(f"Selected after init delay threshold: {gtr._init_delay} ms", np.array(list(gtr._selected_channels_init)), gtr.locations, f"init_delay_threshold{suffix}.png", plot_dir, fresh_plots=fresh_plots)
    plot_selected_channels_helper("Selected after all thresholds", gtr.selected_channels, gtr.locations, f"all_thresholds{suffix}.png", plot_dir, fresh_plots=fresh_plots)
    #return gtr

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

def plot_template_propagation(gtr, plot_dir, unit_id, title, fresh_plots=False):
    fig = plt.figure(figsize=(8, 40))
    ax = fig.add_subplot(111)
    fig_path = plot_dir / f"{title}_template_propagation_unit_{unit_id}.png"
    if os.path.exists(fig_path) and not fresh_plots:
        return
    template = gtr.template
    selected_channels = gtr.selected_channels
    locations = gtr.locations
    ax = av.plotting.plot_template_propagation(
                                                #ax=ax, sort_templates=True
                                                template, locations, selected_channels, sort_templates=True,
                                                color=None, color_marker=None, ax=None
                              )
    ax.set_title(f'Template Propagation for Unit {unit_id} {title}', fontsize=16)
    save_figure(fig, fig_path)

def generate_axon_analytics(gtr, units, branch_ids, velocities, path_lengths, r2s, extremums, init_chans, num_channels_included, 
                            channel_density, 
                            unit_id, transformed_template, trans_loc):
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
                                        #cmap="tab20", plot_bp=True, plot_neighbors=True, plot_full_template=True, ax=axpaths_raw
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
                                        #cmap="tab20", plot_bp=True, plot_full_template=True, ax=axpaths_raw
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

def save_axon_analytics(stream_id, units, extremums, branch_ids, velocities, path_lengths, r2s, num_channels_included, channel_density, init_chans, recon_dir, suffix=""):
    df_mea1k = pd.DataFrame({"unit_ids": units, "unit location": extremums, 
                            "branch_id": branch_ids, "velocity": velocities, "length": path_lengths, "r2": r2s, "num_channels_included": num_channels_included, 
                            "channel_density": channel_density, "init_chan": init_chans})
    
    # df_mea1k = pd.DataFrame({"unit_ids": [units], "unit location": [extremums], 
    #                          "branch_id": [branch_ids], "velocity": [velocities], "length": [path_lengths], "r2": [r2s], "num_channels_included": [num_channels_included], 
    #                          "channel_density": [channel_density], "init_chan": [init_chans]})

                            #  "length": str(path_lengths), "r2": str(r2s), "num_channels_included": str(num_channels_included), 
                            #  "channel_density": str(channel_density), "init_chan": str(init_chans)})
    recon_dir_parent = os.path.dirname(recon_dir)
<<<<<<< HEAD
    if not os.path.exists(recon_dir_parent):
        os.makedirs(recon_dir_parent)
    df_mea1k.to_csv(Path(recon_dir_parent) / f"{stream_id}_axon_analytics{suffix}.csv", index=False)
=======
    if not os.path.exists(recon_dir_parent): os.makedirs(recon_dir_parent)
    df_mea1k.to_csv(Path(recon_dir_parent) / f"{stream_id}_axon_analytics.csv", index=False)
>>>>>>> 1f4fae2 (Major changes to pipeline logic + axon_velocity submod for TK project.)

<<<<<<< HEAD
# lib_plotting_and_analysis.py
=======
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

def plot_template_propagation_wrapper(template, locations, selected_channels, plot_dir, unit_id, fresh_plots=False):
    fig = plt.figure(figsize=(8, 40))
    ax = fig.add_subplot(111)
    fig_path = plot_dir / f"template_propagation_unit_{unit_id}.png"
    if os.path.exists(fig_path) and fresh_plots == False: return

    ax = plot_template_propagation(template, locations, selected_channels, ax=ax, sort_templates=True)
    ax.set_title(f'Template Propagation for Unit {unit_id}', fontsize=16)
    save_figure(fig, fig_path)

def process_unit(unit_id, unit_templates, recon_dir, params, successful_recons, logger=None):
    if logger is not None: logger.info(f'Processing unit {unit_id}')
    else: print(f'Processing unit {unit_id}')
    
    merged_template = unit_templates['merged_template']
    merged_channel_loc = unit_templates['merged_channel_loc']
    merged_template_filled = unit_templates['merged_template_filled']
    merged_channel_filled_loc = unit_templates['merged_channel_locs_filled']
>>>>>>> 82357b0 (Pipeline is functional. Modificaitons to plot generated underway.)

import AxonReconPipeline.axon_velocity.axon_velocity.plotting as av_plotting

<<<<<<< HEAD
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
    #log = kwargs.pop('log', None)
    log = kwargs.pop('log', False)
    template = kwargs.pop('template', None)
    locations = kwargs.pop('locations', None)
    gtr = kwargs.pop('gtr', None)
    elec_size = kwargs.pop('elec_size', 8)
    cmap = kwargs.pop('cmap', 'viridis')
    ax = kwargs.pop('ax', None)
    skip_frames = kwargs.pop('skip_frames', 1)
    interval = kwargs.pop('interval', 10)
    skip_frames = kwargs.pop('skip_frames', 1)
    ani = av_plotting.play_template_map(
        template, locations, gtr=gtr, elec_size=elec_size, cmap=cmap, 
        log=log, ax=ax, skip_frames=skip_frames, interval=interval, 
        #**kwargs
        )
    ani.save(save_path, writer='ffmpeg')

def plot_template_wrapper(save_path, title, fresh_plots = False, **kwargs):
    """
    Wrapper for axon_velocity.plotting.plot_template.

    Parameters:
    data : Data required for plotting the template.
    save_path : Path where the plot will be saved.
    **kwargs : Additional arguments for plot customization.
    """
    if os.path.exists(save_path) and not fresh_plots: return
    else:
        fig = plt.figure(figsize=(100, 200))
        ax = fig.add_subplot(111)
        ax = av_plotting.plot_template(**kwargs)
        ax.set_title(title, fontsize=20)
        if save_path is not None: save_figure(ax, save_path, dpi=2400)

def plot_axon_summary_wrapper(save_path, title, fresh_plots, **kwargs):
    """
    Wrapper for axon_velocity.plotting.plot_axon_summary.

    Parameters:
    data : Data required for plotting the axon summary.
    save_path : Path where the plot will be saved.
    **kwargs : Additional arguments for plot customization.
    """
    if os.path.exists(save_path) and not fresh_plots: return
    else:
        fig, axes = av_plotting.plot_axon_summary(**kwargs)
        fig.set_title(title, fontsize=20)
        if save_path is not None: save_figure(fig, save_path, dpi=600)
=======
    num_channels_included.append(len(merged_channel_loc))

    x_coords = merged_channel_loc[:, 0]
    y_coords = merged_channel_loc[:, 1]
    width = np.max(x_coords) - np.min(x_coords)
    height = np.max(y_coords) - np.min(y_coords)
    area = width * height
    channel_density.append(len(merged_channel_loc) / area) #channels / um^2
    
    try:
        generate_amplitude_map(transformed_template_filled, trans_loc_filled, plot_dir)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate amplitude_map, error: {e}")
        plt.close()
    
    try:
        generate_peak_latency_map(transformed_template_filled, trans_loc_filled, plot_dir)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate peak_latency_map, error: {e}")
        plt.close()
    
    try:
        gtr0 = GraphAxonTracking(transformed_template, trans_loc, 10000, **params)
        select_and_plot_channels(gtr0, plot_dir)
        gtr0 = compute_graph_propagation_velocity(transformed_template, trans_loc, 10000, **params)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to select and plot channels or compute_graph_propagation_velocity, error: {e}")
        plt.close()
        return None
    
    try:
        generate_axon_analytics(gtr0, units, branch_ids, velocities, path_lengths, r2s, extremums, unit_id, transformed_template, trans_loc)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate axon_analytics, error: {e}")
    
    try:
        generate_axon_reconstruction_heuristics(gtr0, plot_dir, unit_id)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate axon_reconstruction_heuristics, error: {e}")
        plt.close() #close the figure, avoid memory leak
    
    try:
        generate_axon_reconstruction_raw(gtr0, plot_dir, unit_id, recon_dir, successful_recons)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate raw axon_reconstruction, error: {e}")
        plt.close() #close the figure, avoid memory leak
    
    try:
        generate_axon_reconstruction_clean(gtr0, plot_dir, unit_id, recon_dir, successful_recons)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate clean axon_reconstruction, error: {e}")
        plt.close() #close the figure, avoid memory leak
    
    try:
        generate_axon_reconstruction_velocities(gtr0, plot_dir, unit_id)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate axon_reconstruction_velocities, error: {e}")
        plt.close() #close the figure, avoid memory leak
    
    try:
        plot_voltage_signals(transformed_template, trans_loc, gtr0.selected_channels, plot_dir, unit_id)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to plot voltage signals, error: {e}")
        plt.close() #close the figure, avoid memory leak

    try:
        plot_all_traces(transformed_template, trans_loc, plot_dir, unit_id)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to plot all voltage traces, error: {e}")
        plt.close() #close the figure, avoid memory leak

    try:
        plot_template_propagation_wrapper(transformed_template, trans_loc, gtr0.selected_channels, plot_dir, unit_id)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to plot template propagation, error: {e}")
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
>>>>>>> 82357b0 (Pipeline is functional. Modificaitons to plot generated underway.)
