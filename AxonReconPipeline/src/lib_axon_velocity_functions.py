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

'''logging setup'''
import logging
logger = logging.getLogger(__name__) #Create a logger
logger.setLevel(logging.DEBUG)
stream_handler = logging.StreamHandler()  # Create handlers, logs to console
stream_handler.setLevel(logging.DEBUG) # Set level of handlers
logger.addHandler(stream_handler) # Add handlers to the logger
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s - %(module)s.%(funcName)s') # Create formatters and add it to handlers
stream_handler.setFormatter(formatter)

'''default paths'''
default_templates_dir='./AxonReconPipeline/data/temp_data/templates'

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
def should_process_stream(stream_id, stream_select):
    if stream_select is not None:
        if stream_id != f'well{stream_select:03}':
            return False
    logger.info(f'Processing stream {stream_id}')
    return True
def get_paths(templates_dir, date, chip_id, scanType, run_id, stream_id, unit_id):
    temp_dir = os.path.dirname(templates_dir)
    merged_template_dir = Path(templates_dir) / date / chip_id / scanType / run_id / stream_id / 'templates' / f'{unit_id}.npy'
    merged_channel_dir = Path(templates_dir) / date / chip_id / scanType / run_id / stream_id / 'templates' / f'{unit_id}_channels.npy'
    merged_template_filled_dir = Path(templates_dir) / date / chip_id / scanType / run_id / stream_id / 'templates' / f'{unit_id}_filled.npy'
    merged_channel_filled_dir = Path(templates_dir) / date / chip_id / scanType / run_id / stream_id / 'templates' / f'{unit_id}_channels_filled.npy'
    return merged_template_dir, merged_channel_dir, merged_template_filled_dir, merged_channel_filled_dir
def load_data(merged_template_dir, merged_channel_dir, merged_template_filled_dir, merged_channel_filled_dir, unit_id):
    try:              
        merged_template = np.load(merged_template_dir)
        merged_channel_loc = np.load(merged_channel_dir)
        merged_template_filled = np.load(merged_template_filled_dir)
        merged_channel_filled_loc = np.load(merged_channel_filled_dir)
    except:
        logger.info(f'No merged template found for {unit_id}')
        return None, None, None, None
    return merged_template, merged_channel_loc, merged_template_filled, merged_channel_filled_loc
def transform_data(merged_template, merged_channel_loc, merged_template_filled, merged_channel_filled_loc):
    transformed_template = merged_template[-1].T
    transformed_template_filled = merged_template_filled[-1].T
    trans_loc = np.array([[loc[0], loc[1]*-1] for loc in merged_channel_loc[-1]])
    trans_loc_filled = np.array([[loc[0], loc[1]*-1] for loc in merged_channel_filled_loc[-1]])
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
    logger.info(f"unit {unit_id} successfully generated axon_reconstruction_heuristics")
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
    logger.info(f"unit {unit_id} successfully generated raw axon_reconstruction")
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
    logger.info(f"unit {unit_id} successfully generated clean axon_reconstruction")
def generate_axon_reconstruction_velocities(gtr0, plot_dir, unit_id, fresh_plots=False):
    fvel, ax_vel = plt.subplots(figsize=(9, 15))
    ax_vel.axis("off")
    fig_path = plot_dir / "axon_reconstruction_velocities.png"
    if os.path.exists(fig_path) and fresh_plots == False: return
    fvel = gtr0.plot_velocities(cmap="tab20", plot_outliers=True, fig=fvel, markersize=12, markersize_out=18, fs=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    save_figure(fvel, fig_path)
    logger.info(f"unit {unit_id} successfully generated axon_reconstruction_velocities")
def save_axon_analytics(stream_id, units, extremums, branch_ids, velocities, path_lengths, r2s, num_channels_included, channel_density, recon_dir):
    df_mea1k = pd.DataFrame({"unit_ids": units, "unit location": extremums, "branch_id": branch_ids, "velocity": velocities, "length": path_lengths, "r2": r2s, "num_channels_included": num_channels_included, "channel_density": channel_density})
    recon_dir_parent = os.path.dirname(recon_dir)
    df_mea1k.to_csv(Path(recon_dir_parent) / f"{stream_id}_axon_analytics.csv", index=False)
def process_unit(unit_id, templates_dir, date, chip_id, scanType, run_id, stream_id, recon_dir, params, successful_recons):
    logger.info(f'Processing unit {unit_id}')
    merged_template_dir, merged_channel_dir, merged_template_filled_dir, merged_channel_filled_dir = get_paths(templates_dir, date, chip_id, scanType, run_id, stream_id, unit_id)
    merged_template, merged_channel_loc, merged_template_filled, merged_channel_filled_loc = load_data(merged_template_dir, merged_channel_dir, merged_template_filled_dir, merged_channel_filled_dir, unit_id)
    
    if merged_template is None:
        return None  # skip if no merged template found
    
    transformed_template, transformed_template_filled, trans_loc, trans_loc_filled = transform_data(merged_template, merged_channel_loc, merged_template_filled, merged_channel_filled_loc)
    plot_dir = create_plot_dir(recon_dir, unit_id)
    units, branch_ids, velocities, path_lengths, r2s, extremums, num_channels_included, channel_density = [], [], [], [], [], [], [], []

    # Assuming merged_channel_loc and merged_channel_filled_loc are numpy arrays of shape (n, 2) where n is the number of channels
    #frac_chans_included.append(len(merged_channel_loc) / len(merged_channel_filled_loc))
    num_channels_included.append(len(merged_channel_loc[0]))

    # Calculate the area covered by the channels
    # Get the bounding box of the channels to estimate the area
    x_coords = merged_channel_loc[0][:, 0]
    y_coords = merged_channel_loc[0][:, 1]

    # Calculate the area of the bounding box
    width = np.max(x_coords) - np.min(x_coords)
    height = np.max(y_coords) - np.min(y_coords)
    area = width * height

    # Calculate the channel density as the number of channels per unit area
    channel_density.append(len(merged_channel_loc[0]) / area) #channels / um^2
    
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
def analyze_and_reconstruct(continuous_h5_file_info, templates_dir=default_templates_dir, params=av.get_default_graph_velocity_params(), stream_select=None, n_jobs = 8):
    for i, continuous_h5_dir in continuous_h5_file_info.items():
        date, chip_id, scanType, run_id = extract_recording_info(continuous_h5_dir)
        stream_ids = get_stream_ids(continuous_h5_dir['h5_file_path'])
        
        for stream_id in stream_ids:
            if not should_process_stream(stream_id, stream_select):
                continue

            temp_dir = os.path.dirname(templates_dir)
            sorting_dir = Path(temp_dir) / 'sortings' / date / chip_id / scanType / run_id / stream_id / 'sorter_output'
            sorting = ss.Kilosort2_5Sorter._get_result_from_folder(sorting_dir)
            recon_dir = Path(temp_dir) / 'reconstructions' / date / chip_id / scanType / run_id / stream_id
            successful_recons = {}
            unit_ids = sorting.get_unit_ids()
            successful_recons[str(recon_dir)] = {}
            successful_recons[str(recon_dir)]["successful_units"] = {}
            all_units, all_branch_ids, all_velocities, all_path_lengths, all_r2s, all_extremums, all_num_channels_included, all_channel_densities = [], [], [], [], [], [], [], []
            logger.info(f'Processing {len(unit_ids)} units, with {n_jobs} workers')
            
            with concurrent.futures.ProcessPoolExecutor(max_workers=n_jobs) as executor:
                futures = [executor.submit(process_unit, unit_id, templates_dir, date, chip_id, scanType, run_id, stream_id, recon_dir, params, successful_recons) for unit_id in unit_ids]
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


def analyze_and_reconstruct_v2(continuous_h5_file_info, templates_dir=default_templates_dir, allowed_scan_types=['AxonTracking'], params=av.get_default_graph_velocity_params(), stream_select=None):
    for i, continuous_h5_dir in continuous_h5_file_info.items():
        date, chip_id, scanType, run_id = extract_recording_info(continuous_h5_dir)
        stream_ids = get_stream_ids(continuous_h5_dir['h5_file_path'])
        
        for stream_id in stream_ids:
            if not should_process_stream(stream_id, stream_select):
                continue

            temp_dir = os.path.dirname(templates_dir)
            sorting_dir = Path(temp_dir) / 'sortings' / date / chip_id / scanType / run_id / stream_id / 'sorter_output'
            sorting = ss.Kilosort2Sorter._get_result_from_folder(sorting_dir)
            recon_dir = Path(temp_dir) / 'reconstructions' / date / chip_id / scanType / run_id / stream_id
            successful_recons = {}
            unit_ids = sorting.get_unit_ids()
            successful_recons[str(recon_dir)] = {}
            successful_recons[str(recon_dir)]["successful_units"] = {}
            units, branch_ids, velocities, path_lengths, r2s, extremums = [], [], [], [], [], []
            
            for unit_id in unit_ids:
                logger.info(f'Processing unit {unit_id}')
                merged_template_dir, merged_channel_dir, merged_template_filled_dir, merged_channel_filled_dir = get_paths(templates_dir, date, chip_id, scanType, run_id, stream_id, unit_id)
                merged_template, merged_channel_loc, merged_template_filled, merged_channel_filled_loc = load_data(merged_template_dir, merged_channel_dir, merged_template_filled_dir, merged_channel_filled_dir, unit_id)
                
                if merged_template is None: continue #skip if no merged template found
                
                transformed_template, transformed_template_filled, trans_loc, trans_loc_filled = transform_data(merged_template, merged_channel_loc, merged_template_filled, merged_channel_filled_loc)
                plot_dir = create_plot_dir(recon_dir, unit_id)

                #if os.path.exists(plot_dir): continue
                
                try: generate_amplitude_map(transformed_template_filled, trans_loc_filled, plot_dir)
                except Exception as e: logger.info(f"unit {unit_id} failed to generate amplitude_map, error: {e}")
                
                try: generate_peak_latency_map(transformed_template_filled, trans_loc_filled, plot_dir)
                except Exception as e: logger.info(f"unit {unit_id} failed to generate peak_latency_map, error: {e}")
                
                try:
                    gtr0 = GraphAxonTracking(transformed_template, trans_loc, 10000, **params)
                    select_and_plot_channels(gtr0, plot_dir)
                    gtr0 = compute_graph_propagation_velocity(transformed_template, trans_loc, 10000, **params)
                except Exception as e: logger.info(f"unit {unit_id} failed to select and plot channels or compute_graph_propagation_velocity, error: {e}")
                
                try: generate_axon_analytics(gtr0, units, branch_ids, velocities, path_lengths, r2s, extremums, unit_id, transformed_template, trans_loc)
                except Exception as e: logger.info(f"unit {unit_id} failed to generate axon_analytics, error: {e}")
                
                try: generate_axon_reconstruction_heuristics(gtr0, plot_dir, unit_id)
                except Exception as e: logger.info(f"unit {unit_id} failed to generate axon_reconstruction_heuristics, error: {e}")
                
                try: generate_axon_reconstruction_raw(gtr0, plot_dir, unit_id, recon_dir, successful_recons)
                except Exception as e: logger.info(f"unit {unit_id} failed to generate raw axon_reconstruction, error: {e}")
                
                try: generate_axon_reconstruction_clean(gtr0, plot_dir, unit_id, recon_dir, successful_recons)
                except Exception as e: logger.info(f"unit {unit_id} failed to generate clean axon_reconstruction, error: {e}")
                
                try: generate_axon_reconstruction_velocities(gtr0, plot_dir, unit_id)
                except Exception as e: logger.info(f"unit {unit_id} failed to generate axon_reconstruction_velocities, error: {e}")
            
            save_axon_analytics(units, extremums, branch_ids, velocities, path_lengths, r2s, recon_dir)
            
            # with open(f'{recon_dir}/successful_recons.json', 'w') as f:
            #     json.dump(successful_recons, f)

'''for reference'''
def analyze_and_reconstruct_og(continuous_h5_file_info, templates_dir = default_templates_dir, params = av.get_default_graph_velocity_params(),  stream_select=None):
    for i, continuous_h5_dir in continuous_h5_file_info.items():
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
                    logger.info(f'Processing stream {stream_id}')
                    break

        temp_dir = os.path.dirname(templates_dir)
        sorting_dir = Path(temp_dir) /'sortings'/ date / chip_id / scanType / run_id / stream_id / 'sorter_output'
        sorting = ss.Kilosort2_5Sorter._get_result_from_folder(sorting_dir)
        recon_dir = Path(temp_dir) /'reconstructions'/ date / chip_id / scanType / run_id / stream_id
        successful_recons = {}  
        unit_ids = sorting.get_unit_ids()
        successful_recons[str(recon_dir)] = {}
        successful_recons[str(recon_dir)]["successful_units"] = {}
        units = []
        branch_ids = []
        velocities = []
        path_lengths = []
        r2s = []
        extremums =[]
        for unit_id in unit_ids:

            print(f'\nProcessing unit {unit_id}')
            # debug
            #if unit_id != 51: continue
            #
                    
            merged_template_dir = Path(templates_dir) / date / chip_id / scanType / run_id / stream_id / 'templates' / f'{unit_id}.npy'
            merged_channel_dir = Path(templates_dir) / date / chip_id / scanType / run_id / stream_id / 'templates' / f'{unit_id}_channels.npy'
            merged_template_filled_dir = Path(templates_dir) / date / chip_id / scanType / run_id / stream_id / 'templates' / f'{unit_id}_filled.npy'
            merged_channel_filled_dir = Path(templates_dir) / date / chip_id / scanType / run_id / stream_id / 'templates' / f'{unit_id}_channels_filled.npy'
            try:              

                merged_template = np.load(merged_template_dir)
                merged_channel_loc = np.load(merged_channel_dir)
                merged_template_filled = np.load(merged_template_filled_dir)
                merged_channel_filled_loc = np.load(merged_channel_filled_dir)

                #merged_template = np.load('/mnt/disk15tb/adam/git_workspace/MEA_Analysis/AxonReconPipeline/data/temp_data/templates/240118/M06844/AxonTracking/000026/well001/partial_templates/seg3_unit17.npy')
                #merged_channel_loc = np.load('/mnt/disk15tb/adam/git_workspace/MEA_Analysis/AxonReconPipeline/data/temp_data/templates/240118/M06844/AxonTracking/000026/well001/partial_templates/seg3_unit17_channels.npy')
            except:
                print(f'No merged template found for {unit_id}')
                continue

            ## Channels
            template = merged_template[-1]
            transformed_template = template.T
            #mirror matrix about x-axis to match maxwell orientation
            #transformed_template = np.flip(transformed_template, axis=0)
            
            template_filled = merged_template_filled[-1]
            transformed_template_filled = template_filled.T
            #mirror matrix about x-axis to match maxwell orientation\
            #transformed_template_filled = np.flip(transformed_template_filled, axis=0)
            
            ## Locations
            trans_loc = merged_channel_loc[-1]
            #multiply all y values by -1 to match maxwell orientation
            trans_loc = np.array([[loc[0], loc[1]*-1] for loc in trans_loc])

            trans_loc_filled = merged_channel_filled_loc[-1]
            #multiply all y values by -1 to match maxwell orientation
            trans_loc_filled = np.array([[loc[0], loc[1]*-1] for loc in trans_loc_filled])
            #mirror matrix about x-axis to match maxwell orientation
            #trans_loc_filled = np.flip(trans_loc_filled, axis = 0)

            ##
            gtrs = {}
            fs = 10000
            plot_dir = Path(recon_dir) / f'unit_{unit_id}' / 'axon_tracking_plots'
            if os.path.exists(plot_dir):
                #print(f'Plot directory already exists for unit {unit_id}, skipping...')
                #continue
                shutil.rmtree(plot_dir)
            key = 0

            # Generate amplitude_map figure
            try:
                fig_amp = plt.figure(figsize=(10, 5))
                ax_amp = fig_amp.add_subplot(111)
                #set_trace()
                ax_amp = plot_amplitude_map(transformed_template_filled, 
                                            trans_loc_filled, 
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
                plt.close()
                print(f"unit {unit_id} failed to generate amplitude_map, error: {e}")

            # Generate peak latency map figure
            try:
                fig_peaks = plt.figure(figsize=(10, 5))
                ax_peaks = fig_peaks.add_subplot(111)
                ax_peaks = plot_peak_latency_map(transformed_template_filled, 
                                                trans_loc_filled, fs=fs, log=False, 
                                                ax=ax_peaks, colorbar=True, 
                                                colorbar_orientation="horizontal")
                ax_peaks.set_title(f"Peak latency", fontsize=20)

                #Save figure
                fig_name = "peak_latency_map.png"
                fig_path = plot_dir / fig_name
                plt.savefig(fig_path, dpi=300, bbox_inches='tight')
                plt.close()

                #accept_units_peak_latency_map.append(key)
            except Exception as e:
                #reject_units_peak_latency_map.append(key)
                plt.close()
                print(f"unit {unit_id} failed to generate peak_latency_map, error: {e}")

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
                gtr0 = None #initialize gtr0
                gtr0 = GraphAxonTracking(transformed_template, trans_loc, fs, 
                                        #verbose=verbose_bool, 
                                        **params)           

                # Select channels
                gtr0.select_channels()                

                # Use the function for each block
                plot_selected_channels(f"Selected after detection threshold: {gtr0._detect_threshold} ms", 
                                    np.array(list(gtr0._selected_channels_detect)), gtr0.locations, "detection_threshold.png")
                plot_selected_channels(f"Selected after detection threshold: {gtr0._kurt_threshold} ms", 
                                    np.array(list(gtr0._selected_channels_kurt)), gtr0.locations, "kurt_threshold.png")
                plot_selected_channels(f"Selected after peak std threshold: {gtr0._peak_std_threhsold} ms", 
                                        np.array(list(gtr0._selected_channels_peakstd)), gtr0.locations, "peak_std_threshold.png")
                plot_selected_channels(f"Selected after init delay threshold: {gtr0._init_delay} ms", 
                                        np.array(list(gtr0._selected_channels_init)), gtr0.locations, "init_delay_threshold.png")
                plot_selected_channels("Selected after all thresholds", gtr0.selected_channels, gtr0.locations, "all_thresholds.png")

                #if all thresholds are passed, compute graph propagation velocity
                try:
                    gtr0 = compute_graph_propagation_velocity(transformed_template, trans_loc, fs, 
                        #verbose=False, 
                        **params)
                except Exception as e:
                    print(f"unit {unit_id} failed to compute_graph_propagation_velocity, error: {e}")
                #accept_units_selected_channels.append(key)
            except Exception as e:
                #reject_units_selected_channels.append(key)
                plt.close()
                print(f"unit {unit_id} failed to generate selected_channels, error: {e}")              
            
            # Generate axon analytics
            try:

                for i, br in enumerate(gtr0.branches):
                    path = br["channels"]
                    velocity = br["velocity"]
                    r2 = br["r2"]
                    length = gtr0.compute_path_length(path)
                    units.append(unit_id)
                    branch_ids.append(i)
                    velocities.append(velocity)
                    path_lengths.append(length)
                    r2s.append(r2)
                    extremums.append(get_extremum(transformed_template, trans_loc))

                #accept_units_axon_reconstruction_heuristics.append(key)
            except Exception as e:
                #reject_units_axon_reconstruction_heuristics.append(key)
                plt.close()
                print(f"unit {unit_id} failed to generate axon_analytics, error: {e}")

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
                print(f"unit {unit_id} successfully generated axon_reconstruction_heuristics")

                #accept_units_axon_reconstruction_heuristics.append(key)
            except Exception as e:
                #reject_units_axon_reconstruction_heuristics.append(key)
                plt.close()
                print(f"unit {unit_id} failed to generate axon_reconstruction_heuristics, error: {e}")

            # # Generate axon_reconstruction figures
            # try:
            #     fig_graph = plt.figure(figsize=(12, 7))
            #     #ax_raw = fig_graph.add_subplot(111)                
            #     #axpaths_raw = ax_raw
            #     axpaths_raw = gtr0.plot_branches(cmap="tab20", fig = fig_graph,
            #                                      #plot_bp=True, plot_neighbors=True, 
            #                                         #plot_full_template=True, 
            #                                      #ax=axpaths_raw
            #                                      )
            #     #axpaths_raw.legend(fontsize=6)
            #     #axpaths_raw.set_title("All Branches")

            #     fig_name = "axon_reconstruction_all.png"
            #     fig_path = plot_dir / fig_name
            #     plt.savefig(fig_path, dpi=600, bbox_inches='tight')
            #     plt.close()
                
            #     successful_recons[str(recon_dir)]["successful_units"][str(unit_id)] = {}
            #     #successful_recons[str(recon_dir)]["successful_units"][str(unit_id)]["stream_id"] = stream_id
            #     #successful_recons[str(recon_dir)]["successful_units"][str(unit_id)]["fig_path"] = fig_path
            #     print(f"unit {unit_id} successfully generated axon_reconstruction")
            #     #successful_recons[continuous_h5_dir]["template"] = template
            #     #successful_recons[continuous_h5_dir]["merged_channel_loc"] = merged_channel_loc
            #     # gtrs[key] = compute_graph_propagation_velocity(template.T, merged_channel_loc[i], fs, 
            #     #                                         verbose=False, **params)

            #     #accept_units_axon_reconstruction.append(key)
            #     #flag = True
            #     #break
            # except Exception as e:
            #     #reject_units_axon_reconstruction.append(key)
            #     plt.close()
            #     print(f"unit {unit_id} failed to generate all axon_reconstruction, error: {e}")
            
            # Generate axon_reconstruction figures
            try:
                fig_graph = plt.figure(figsize=(12, 7))
                ax_raw = fig_graph.add_subplot(111)                
                axpaths_raw = ax_raw
                axpaths_raw = gtr0.plot_clean_branches(cmap="tab20", plot_bp=True, 
                                                    #plot_neighbors=True, 
                                                    plot_full_template=True, ax=axpaths_raw)
                axpaths_raw.legend(fontsize=6)
                axpaths_raw.set_title("Clean Branches")

                fig_name = "axon_reconstruction_clean.png"
                fig_path = plot_dir / fig_name
                plt.savefig(fig_path, dpi=600, bbox_inches='tight')
                plt.close()
                
                successful_recons[str(recon_dir)]["successful_units"][str(unit_id)] = {}
                #successful_recons[str(recon_dir)]["successful_units"][str(unit_id)]["stream_id"] = stream_id
                #successful_recons[str(recon_dir)]["successful_units"][str(unit_id)]["fig_path"] = fig_path
                print(f"unit {unit_id} successfully generated clean axon_reconstruction")
                #successful_recons[continuous_h5_dir]["template"] = template
                #successful_recons[continuous_h5_dir]["merged_channel_loc"] = merged_channel_loc
                # gtrs[key] = compute_graph_propagation_velocity(template.T, merged_channel_loc[i], fs, 
                #                                         verbose=False, **params)

                #accept_units_axon_reconstruction.append(key)
                #flag = True
                #break
            except Exception as e:
                #reject_units_axon_reconstruction.append(key)
                plt.close()
                print(f"unit {unit_id} failed to generate clean axon_reconstruction, error: {e}")

            # Generate axon_reconstruction figures
            try:
                fig_graph = plt.figure(figsize=(12, 7))
                ax_raw = fig_graph.add_subplot(111)                
                axpaths_raw = ax_raw
                axpaths_raw = gtr0.plot_raw_branches(cmap="tab20", plot_bp=True, plot_neighbors=True, 
                                                    plot_full_template=True, ax=axpaths_raw)
                axpaths_raw.legend(fontsize=6)
                axpaths_raw.set_title("Raw Branches")

                fig_name = "axon_reconstruction_raw.png"
                fig_path = plot_dir / fig_name
                plt.savefig(fig_path, dpi=600, bbox_inches='tight')
                plt.close()
                
                successful_recons[str(recon_dir)]["successful_units"][str(unit_id)] = {}
                #successful_recons[str(recon_dir)]["successful_units"][str(unit_id)]["stream_id"] = stream_id
                #successful_recons[str(recon_dir)]["successful_units"][str(unit_id)]["fig_path"] = fig_path
                print(f"unit {unit_id} successfully generated axon_reconstruction")
                #successful_recons[continuous_h5_dir]["template"] = template
                #successful_recons[continuous_h5_dir]["merged_channel_loc"] = merged_channel_loc
                # gtrs[key] = compute_graph_propagation_velocity(template.T, merged_channel_loc[i], fs, 
                #                                         verbose=False, **params)

                #accept_units_axon_reconstruction.append(key)
                #flag = True
                #break
            except Exception as e:
                #reject_units_axon_reconstruction.append(key)
                plt.close()
                print(f"unit {unit_id} failed to generate raw axon_reconstruction, error: {e}")
                
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
                print(f"unit {unit_id} successfully generated axon_reconstruction_velocities")

                #accept_units_axon_reconstruction_velocities.append(key)
            except Exception as e:
                #reject_units_axon_reconstruction_velocities.append(key)
                plt.close()
                print(f"unit {unit_id} failed to generate axon_reconstruction_velocities, error: {e}")

        import pandas as pd
        df_mea1k = pd.DataFrame({"unit_ids": units, "unit location": extremums, "branch_id": branch_ids, "velocity": velocities,
                            "length": path_lengths, "r2": r2s})
        df_mea1k.to_csv(Path(recon_dir) / "axon_analytics.csv", index=False)
        
        # def convert_keys_to_int(d):
        #     if isinstance(d, dict):
        #         return {int(k) if isinstance(k, np.int64) else k: convert_keys_to_int(v) for k, v in d.items()}
        #     elif isinstance(d, list):
        #         return [convert_keys_to_int(v) for v in d]
        #     else:
        #         return d
        with open(f'{recon_dir}/successful_recons.json', 'w') as f:
            json.dump(successful_recons, f)