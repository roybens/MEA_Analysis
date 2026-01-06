'''process each unit template, generate axon reconstructions and analytics'''

from venv import logger
from lib_axon_velocity_functions import *
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
from axon_velocity import GraphAxonTracking
from axon_velocity import get_default_graph_velocity_params
from copy import deepcopy
import shutil
from scipy.signal import find_peaks
from scipy.spatial import KDTree


def sinc_interp(x, s, u):
    """
    Interpolates x at the indices u, using the indices s. This function assumes x is sampled
    at regular intervals given by s, and interpolates it to find the values at the positions u.

    Parameters:
    - x: The original signal.
    - s: The original sample indices.
    - u: The desired interpolated indices.

    Returns:
    - The interpolated signal.
    """
    sinc_matrix = np.tile(u, (len(s), 1)) - np.tile(s[:, np.newaxis], (1, len(u)))
    return np.dot(x, np.sinc(sinc_matrix))

def upsample_whitaker_shannon(template, target_fs=200000, original_fs=10000):
    """
    Upsamples the data using Whitaker-Shannon interpolation.
    
    Parameters:
    - template: 2D array where each row corresponds to a channel and each column is a time point.
    - target_fs: The target sampling rate in Hz (default is 200 kHz).
    - original_fs: The original sampling rate in Hz (default is 10 kHz).
    
    Returns:
    - upsampled_template: The upsampled template.
    """
    num_channels, num_samples = template.shape
    
    # Original sample indices
    original_indices = np.arange(num_samples)
    
    # Desired upsampled indices
    upsample_factor = target_fs / original_fs
    new_num_samples = int(num_samples * upsample_factor)
    upsampled_indices = np.linspace(0, num_samples - 1, new_num_samples)
    
    # Perform the interpolation
    upsampled_template = np.zeros((num_channels, new_num_samples))
    
    for channel in range(num_channels):
        upsampled_template[channel, :] = sinc_interp(template[channel, :], original_indices, upsampled_indices)
    
    return upsampled_template

def estimate_background_noise(template, peak_threshold=0.01):
    """
    Cuts out windows around each spike in the template and estimates the background noise for each channel.
    
    Parameters:
    - template: 2D array where each row corresponds to a channel and each column is a time point.
    - peak_threshold: Threshold multiplier for peak detection based on the maximum amplitude across all channels.
    
    Returns:
    - noise_std_per_channel: A list of estimated standard deviations of the background noise for each channel.
    """
    
    num_channels, num_timepoints = template.shape
    
    # Calculate the global maximum amplitude across all channels
    global_max_amplitude = np.max(np.abs(template))
    
    # Set the threshold for peak detection
    threshold = peak_threshold * global_max_amplitude
    
    # Initialize a list to store the noise standard deviations for each channel
    noise_std_per_channel = []
    
    for channel in range(num_channels):
        # Get the template for the current channel
        signal = template[channel, :]
        
        # Calculate the baseline (average of the signal)
        baseline = np.mean(signal)
        
        # Find peaks above the threshold
        peaks, _ = find_peaks(signal, height=threshold)
        
        # Initialize a mask to include all data points
        mask = np.ones(num_timepoints, dtype=bool)
        
        # Exclude windows around each detected peak
        for peak in peaks:
            # Find the start of the peak (when it leaves the baseline)
            start = peak
            while start > 0 and signal[start] > baseline:
                start -= 1
            
            # Find the end of the peak (when it returns to the baseline)
            end = peak
            while end < num_timepoints and signal[end] > baseline:
                end += 1
            
            # Update the mask to exclude the peak window
            mask[start:end] = False
        
        # Collect noise data excluding windows around peaks
        noise_data = signal[mask]

        # Calculate the standard deviation of the noise for this channel
        noise_std_per_channel.append(np.std(noise_data))
    
    return noise_std_per_channel


def divide_into_frames(template, frame_interval=50e-6, original_fs=10000):
    target_fs = 200000  # 200 kHz for finer temporal resolution
    upsampled_template = upsample_whitaker_shannon(template, target_fs=target_fs, original_fs=original_fs)
    
    num_channels, num_samples_upsampled = upsampled_template.shape
    
    # Calculate the number of samples corresponding to the frame interval
    samples_per_frame = int(frame_interval * target_fs)
    
    # Generate the frame indices
    frame_indices = np.arange(0, num_samples_upsampled, samples_per_frame)
    
    # Extract frames
    frames = [upsampled_template[:, i:i+1] for i in frame_indices]
    
    # Find the homing frame and shift
    homing_frame_idx = find_homing_frame(frames)
    shifted_frames = shift_frames_to_include_homing(frames, homing_frame_idx)
    
    return shifted_frames

def find_homing_frame(frames):
    """
    Find the homing frame, which contains the absolute peak value in time.
    Args:
        frames (list of numpy.ndarray): List of frames from `divide_into_frames`.
    
    Returns:
        int: The index of the homing frame.
    """
    max_peak = -np.inf
    homing_frame_idx = 0
    for idx, frame in enumerate(frames):
        peak_in_frame = np.max(np.abs(frame))
        if peak_in_frame > max_peak:
            max_peak = peak_in_frame
            homing_frame_idx = idx
    return homing_frame_idx

def shift_frames_to_include_homing(frames, homing_frame_idx):
    """
    Shift timing of frames so that the homing frame is included.
    Args:
        frames (list of numpy.ndarray): List of frames from `divide_into_frames`.
        homing_frame_idx (int): The index of the homing frame.
    
    Returns:
        list: The shifted frames.
    """
    return frames[homing_frame_idx:] + frames[:homing_frame_idx]

from scipy.spatial import KDTree

def channel_selector(time_frames, noise_thresholds, radius_thresholds, locations, previously_selected_channels=None):
    selected_channels = []

    for pass_idx in range(3):  # Three passes
        noise_threshold = noise_thresholds[pass_idx]
        radius_threshold = radius_thresholds[pass_idx]
        
        for t, frame in enumerate(time_frames):
            peak_channels = np.where(frame > noise_threshold)[0]
            
            groups = group_contiguous_channels(peak_channels, locations, radius_threshold)
            
            if pass_idx == 0:
                selected_in_frame = [max(group, key=lambda ch: frame[ch]) for group in groups]
            else:
                if previously_selected_channels is None or len(previously_selected_channels) == 0:
                    selected_in_frame = peak_channels  # If there are no previous channels, just select from peak_channels
                else:
                    kd_tree = KDTree(locations[previously_selected_channels])  # Use locations of previously selected channels
                    selected_in_frame = []
                    for ch in peak_channels:
                        dist, idx = kd_tree.query([locations[ch]], distance_upper_bound=radius_threshold)
                        if dist[0] < np.inf:  # If within radius, select the channel
                            selected_in_frame.append(ch)

            selected_channels.extend(selected_in_frame)
            previously_selected_channels = np.array(selected_in_frame, dtype=int)

    return selected_channels




from scipy.spatial import KDTree

def group_contiguous_channels(peak_channels, locations, radius_threshold):
    # Convert radius threshold from microns to index distance if needed
    radius_threshold = radius_threshold / 17.5  # assuming 17.5 microns between channels
    
    # Create a KDTree using the locations of the peak channels
    kdtree = KDTree(locations[peak_channels])
    
    groups = []
    visited = set()
    
    for i, ch in enumerate(peak_channels):
        if i in visited:
            continue
        # Query the KDTree to find all neighbors within the radius threshold
        indices = kdtree.query_ball_point(locations[ch], radius_threshold)
        group = set(indices)
        
        # Mark these indices as visited
        visited.update(group)
        
        # Convert the set of indices back to the actual channel numbers
        grouped_channels = peak_channels[list(group)]
        groups.append(grouped_channels)
    
    return groups



def skeletonize(template, locations, params, plot_dir, fresh_plots, unit_id):
    # Estimate the noise standard deviation
    noise_std = estimate_background_noise(template)

    # Ensure the plot directory exists
    frame_dir = os.path.join(plot_dir, f"frames_unit_{unit_id}")
    os.makedirs(frame_dir, exist_ok=True)

    # Template Propagation - show signal propagation through channels selected
    skeleton_params = params.copy()
    skeleton_params.update({
        'detect_threshold': 0.05,
        'remove_isolated': False,
        'peak_std_threshold': None, 
        'init_delay': 0.1,
    })
    
    # Generate GraphAxonTracking object
    gtr = GraphAxonTracking(template, locations, 10000, **skeleton_params)
    gtr.noise_std = noise_std

    gtr.select_channels()
    try: 
        plot_template_propagation(gtr, plot_dir, unit_id, title='test_prop', fresh_plots=fresh_plots)
    except Exception as e: 
        plt.close()

    # Use Channel_Selector to perform the multi-pass selection
    time_frames = divide_into_frames(template)
    homing_frame_idx = find_homing_frame(time_frames)
    time_frames = shift_frames_to_include_homing(time_frames, homing_frame_idx)
    
    # Set thresholds for passes
    noise_thresholds = [9 * gtr.noise_std, 2 * gtr.noise_std, 1 * gtr.noise_std]
    radius_thresholds = [4000, 50, 100]

    selected_channels = channel_selector(time_frames, noise_thresholds, radius_thresholds)

    if selected_channels is None:
        print("Error: channel_selector returned None")
        return []
    
    # Now use these selected channels for further processing
    gtr.selected_channels = selected_channels

    try: 
        plot_template_propagation(gtr, plot_dir, unit_id, title='test_prop', fresh_plots=fresh_plots)
    except Exception as e: 
        plt.close()

    # Template Plot
    try:
        kwargs = {
            'save_path': f"{plot_dir}/template_ploted_{unit_id}.png",
            'title': f'Template {unit_id}',
            'fresh_plots': True,
            'template': gtr.template,
            'locations': gtr.locations,
            'lw': 0.1, # line width
        }
        plot_template_wrapper(**kwargs)
    except Exception as e:
        plt.close()

    direct_links = []
    selected_channels = gtr.selected_channels

    # Loop through selected channels from last to first
    for i in range(len(selected_channels) - 1, 0, -1):
        current_peak = selected_channels[i]
        next_peak = selected_channels[i - 1]
        distance = np.linalg.norm(gtr.locations[current_peak] - gtr.locations[next_peak])
        if distance <= 100:
            direct_links.append((current_peak, next_peak))

            # Plot the direct links incrementally
            plt.figure(figsize=(10, 8))
            plt.scatter(gtr.locations[:, 0], gtr.locations[:, 1], c='gray', label='All Channels')
            plt.scatter(gtr.locations[selected_channels, 0], gtr.locations[selected_channels, 1], c='blue', label='Selected Channels')
            
            for link in direct_links:
                start, end = link
                plt.plot([gtr.locations[start, 0], gtr.locations[end, 0]], 
                         [gtr.locations[start, 1], gtr.locations[end, 1]], 
                         c='red', lw=2)

            plt.xlabel('X Coordinate')
            plt.ylabel('Y Coordinate (Inverted)')
            plt.title(f'Direct Links between Channels for Unit {unit_id} - Frame {len(direct_links)}')
            plt.legend()
            plt.savefig(os.path.join(frame_dir, f"frame_{len(direct_links)}.png"))
            plt.close()

    return direct_links

def approx_milos_tracking(template_dvdt, locations, params, plot_dir, fresh_plots, **kwargs):
    #skeletonize
    #unit_id = 0
    #skeleton = skeletonize(template_dvdt, locations, params, plot_dir, fresh_plots, unit_id)
    
    #simulated_background_signals = estimate_background_noise(template_dvdt, peak_threshold=params['detect_threshold'])
    #noise = np.array(simulated_background_signals)
    #noise = np.nanmean(noise)
    gtr = GraphAxonTracking(template_dvdt, locations, 10000, **params)

    successful_recons = {str(plot_dir): {"successful_units": {}}}
    analysis_options = {'generate_animation': False}
    
    # Calculate noise standard deviation
    noise_std = np.std(template_dvdt)  
    
    time_frames = divide_into_frames(template_dvdt)
    noise_thresholds = [9 * noise_std, 2 * noise_std, 1 * noise_std]
    radius_thresholds = [4000, 50, 100]
    
    # Pass the locations to the channel_selector function
    custom_selected_channels = channel_selector(time_frames, noise_thresholds, radius_thresholds, locations)
    
    #skeletonized_pathways = skeletonize(template_dvdt, locations, params, plot_dir, fresh_plots, unit_id=0, selected_channels=custom_selected_channels)
    
    gtr.selected_channels = custom_selected_channels
    gtr.track_axons()
    
    paths_files_generated = []
    unit_id = 485 
    suffix = ""  

    try:
        # Generate Axon Reconstruction Heuristics
        generate_axon_reconstruction_heuristics(gtr, plot_dir, unit_id, suffix=f"_{suffix}", fresh_plots=True, figsize=(20, 10))
        paths_files_generated.append(plot_dir / f"axon_reconstruction_heuristics{suffix}.png")
    except Exception as e:
        if logger: 
            logger.info(f"unit {unit_id}_{suffix} failed to generate axon_reconstruction_heuristics, error: {e}")
        plt.close()

    # Generate Axon Reconstruction Raw
    try:
        generate_axon_reconstruction_raw(gtr, plot_dir, unit_id, plot_dir, successful_recons, suffix=f"_{suffix}", fresh_plots=True, plot_full_template=False)
        generate_axon_reconstruction_raw(gtr, plot_dir, unit_id, plot_dir, successful_recons, suffix=f"_{suffix}_full", fresh_plots=True, plot_full_template=True)
        paths_files_generated.append(plot_dir / f"axon_reconstruction_raw{suffix}.png")
        paths_files_generated.append(plot_dir / f"axon_reconstruction_raw{suffix}_full.png")
    except Exception as e:
        if logger: 
            logger.info(f"unit {unit_id}_{suffix} failed to generate axon_reconstruction_raw, error: {e}")
        plt.close()

    # Generate Axon Reconstruction Clean
    try:
        generate_axon_reconstruction_clean(gtr, plot_dir, unit_id, plot_dir, successful_recons, suffix=f"_{suffix}", fresh_plots=True, plot_full_template=False)
        generate_axon_reconstruction_clean(gtr, plot_dir, unit_id, plot_dir, successful_recons, suffix=f"_{suffix}_full", fresh_plots=True, plot_full_template=True)
        paths_files_generated.append(plot_dir / f"axon_reconstruction_clean{suffix}.png")
        paths_files_generated.append(plot_dir / f"axon_reconstruction_clean{suffix}_full.png")
    except Exception as e:
        if logger: 
            logger.info(f"unit {unit_id}_{suffix} failed to generate axon_reconstruction_clean, error: {e}")
        plt.close()

    # Generate Axon Reconstruction Velocities
    try:
        generate_axon_reconstruction_velocities(gtr, plot_dir, unit_id, suffix=f"_{suffix}")
        paths_files_generated.append(plot_dir / f"axon_reconstruction_velocities{suffix}.png")
    except Exception as e:
        if logger: 
            logger.info(f"unit {unit_id}_{suffix} failed to generate axon_reconstruction_velocities, error: {e}")
        plt.close()

    # Optional: Animation generation if switch is turned on
    switch = False  
    if switch:
        try:
            assert analysis_options['generate_animation'] == True, "generate_animation set to False. Skipping template map animation generation."
            kwargs = {
                'template': template_dvdt,
                'locations': locations,
                'gtr': gtr,
                'elec_size': 8,
                'cmap': 'viridis',
                'log': False,
                'save_path': f"{plot_dir}/template_map_{unit_id}_{suffix}.gif",
            }
            play_template_map_wrapper(**kwargs)
            paths_files_generated.append(kwargs['save_path'])
        except Exception as e:
            if logger: 
                logger.info(f"unit {unit_id}_{suffix} failed to play template map, error: {e}")
            plt.close()

    # #First pass
    # std_level = 9
    # params9 = params.copy()
    # params9.update({
    #     'detection_type': 'absolute',
    #     'detect_threshold': noise*std_level,
    #     'kurt_threshold': params['kurt_threshold']*std_level,
    #     'remove_isolated': False,
    #     'peak_std_threshold': None,
    #     'init_delay': 0.05,
    #     'max_distance_to_init': 4000.0,
    #     'max_distance_for_edge': 4000.0,
    #     'edge_dist_amp_ratio': 1,
    # })
    # milosh_dir = Path(os.path.join(plot_dir, 'milos_frames'))
    # if os.path.exists(milosh_dir) is False: os.makedirs(milosh_dir)
    # suffix = 'frame1'
    # gtr9 = GraphAxonTracking(template_dvdt, locations, 10000, **params9)
    # gtr9.select_channels()
    # plot_selected_channels_helper(f"Selected after detection threshold: {gtr9._detect_threshold} uV", np.array(list(gtr9._selected_channels_detect)), gtr9.locations, f"detection_threshold{suffix}.png", milosh_dir, fresh_plots=fresh_plots)
    # plot_selected_channels_helper(f"Selected after kurtosis threshold: {gtr9._kurt_threshold}", np.array(list(gtr9._selected_channels_kurt)), gtr9.locations, f"kurt_threshold{suffix}.png", milosh_dir, fresh_plots=fresh_plots)
    # plot_selected_channels_helper(f"Selected after init delay threshold: {gtr9._init_delay} ms", np.array(list(gtr9._selected_channels_init)), gtr9.locations, f"init_delay_threshold{suffix}.png", milosh_dir, fresh_plots=fresh_plots)
    # plot_selected_channels_helper("Selected after all thresholds", gtr9.selected_channels, gtr9.locations, f"all_thresholds{suffix}.png", milosh_dir, fresh_plots=fresh_plots)
    # print()

def process_unit(unit_id, unit_templates, recon_dir, params, analysis_options, successful_recons, failed_recons, logger=None, density_dict=None):
    if logger is not None: 
        logger.info(f'Processing unit {unit_id}')
    else: 
        print(f'Processing unit {unit_id}')

    templates = {
        'merged_template': unit_templates['merged_template'],
        'dvdt_merged_template': unit_templates['dvdt_merged_template'],
    }

    transformed_templates = {}
    for key, template in templates.items():
        if template is not None:
            transformed_template, transformed_template_filled, trans_loc, trans_loc_filled = transform_data(
                template, 
                unit_templates['merged_channel_loc'], 
                unit_templates['merged_template_filled'], 
                unit_templates['merged_channel_locs_filled']
            )
            transformed_templates[key] = (transformed_template, transformed_template_filled, trans_loc, trans_loc_filled)

    plot_dir = create_plot_dir(recon_dir, unit_id)
    analytics_data = {key: {
        'units': [], 'branch_ids': [], 'velocities': [], 'path_lengths': [], 'r2s': [], 'extremums': [], 
        'num_channels_included': [], 'channel_density': [], 'init_chans': []
    } for key in templates.keys()}

    x_coords = unit_templates['merged_channel_loc'][:, 0]
    y_coords = unit_templates['merged_channel_loc'][:, 1]
    width = np.max(x_coords) - np.min(x_coords)
    height = np.max(y_coords) - np.min(y_coords)
    area = width * height
    channel_density_value = len(unit_templates['merged_channel_loc']) / area  # channels / um^2

    # Store the channel density value in the dictionary
    if density_dict is not None:
        density_dict[unit_id] = channel_density_value
        print(f"Stored channel density for unit {unit_id}: {channel_density_value}")
    
    failed_recons[unit_id] = []

    for key, (transformed_template, transformed_template_filled, trans_loc, trans_loc_filled) in transformed_templates.items():
        paths_files_generated = []
        suffix = key.split('_')[-1] if '_' in key else 'template'
        if 'dvdt' in key: suffix = 'dvdt'
        gtr = None
        # Generate Gtr object
        try: gtr = GraphAxonTracking(transformed_template, trans_loc, 10000, **params)
        except Exception as e:
            if logger: logger.info(f"unit {unit_id}_{suffix} failed to initialize GraphAxonTracking for {key}, error: {e}")

        if gtr is not None:
            analytics_data[key]['num_channels_included'].append(len(unit_templates['merged_channel_loc']))
            analytics_data[key]['channel_density'].append(channel_density_value)

            # Template Plot
            try:
                kwargs = {
                    'save_path': f"{plot_dir}/template_ploted_{unit_id}_{suffix}.png",
                    'title': f'Template {unit_id} {suffix}',
                    'fresh_plots': True,
                    'template': gtr.template,
                    'locations': gtr.locations,
                    'lw': 0.1, # line width
                }
                plot_template_wrapper(**kwargs)
                paths_files_generated.append(kwargs['save_path'])
            except Exception as e:
                if logger: 
                    logger.info(f"unit {unit_id}_{suffix} failed to plot template for {key}, error: {e}")
                plt.close()
            
            # Amplitude Maps
            try:
                title = f'{key}_{suffix}'
                generate_amplitude_map(transformed_template_filled, trans_loc_filled, plot_dir, title=f'{key}_{suffix}', fresh_plots=True, log=False, 
                                       cmap='terrain',
                                       )
                paths_files_generated.append(plot_dir / f"{title}_amplitude_map.png")
            except Exception as e:
                if logger: 
                    logger.info(f"unit {unit_id}_{suffix} failed to generate amplitude_map for {key}, error: {e}")
                plt.close()

            # Peak Latency Map
            try:
                generate_peak_latency_map(transformed_template_filled, trans_loc_filled, plot_dir, title=f'{key}_{suffix}', fresh_plots=True, log=False, 
                                       cmap='terrain',
                                       )
                paths_files_generated.append(f"{plot_dir}/{key}_{suffix}_peak_latency_map.png")
            except Exception as e:
                if logger: 
                    logger.info(f"unit {unit_id}_{suffix} failed to generate peak_latency_map for {key}, error: {e}")
                plt.close()
            
            # Select Channels
            try:
                fresh_plots = True
                gtr.select_channels()
                plot_selected_channels_helper(f"Selected after detection threshold: {gtr._detect_threshold} uV", np.array(list(gtr._selected_channels_detect)), gtr.locations, f"detection_threshold{suffix}.png", plot_dir, fresh_plots=fresh_plots)
                plot_selected_channels_helper(f"Selected after kurtosis threshold: {gtr._kurt_threshold}", np.array(list(gtr._selected_channels_kurt)), gtr.locations, f"kurt_threshold{suffix}.png", plot_dir, fresh_plots=fresh_plots)
                plot_selected_channels_helper(f"Selected after peak std threshold: {gtr._peak_std_threhsold} ms", np.array(list(gtr._selected_channels_peakstd)), gtr.locations, f"peak_std_threshold{suffix}.png", plot_dir, fresh_plots=fresh_plots)
                plot_selected_channels_helper(f"Selected after init delay threshold: {gtr._init_delay} ms", np.array(list(gtr._selected_channels_init)), gtr.locations, f"init_delay_threshold{suffix}.png", plot_dir, fresh_plots=fresh_plots)
                plot_selected_channels_helper("Selected after all thresholds", gtr.selected_channels, gtr.locations, f"all_thresholds{suffix}.png", plot_dir, fresh_plots=fresh_plots)
                paths_files_generated.append(f"{plot_dir}/detection_threshold{suffix}.png")
                paths_files_generated.append(f"{plot_dir}/kurt_threshold{suffix}.png")
                paths_files_generated.append(f"{plot_dir}/peak_std_threshold{suffix}.png")
                paths_files_generated.append(f"{plot_dir}/init_delay_threshold{suffix}.png")
                paths_files_generated.append(f"{plot_dir}/all_thresholds{suffix}.png")
            except Exception as e:
                if logger: 
                    logger.info(f"unit {unit_id}_{suffix} failed to select and plot channels for {key}, error: {e}")
                plt.close()

            switch = False
            if switch:
                gtr_milos = None
                if 'dvdt' in key:
                    kwargs = {
                        'template': transformed_template,
                        'locations': trans_loc,
                        'params': params,
                        'plot_dir': plot_dir,
                        'fresh_plots': True,
                    }
                    try: gtr_milos = approx_milos_tracking(**kwargs)
                    except Exception as e:
                        if logger: logger.info(f"unit {unit_id}_{suffix} failed to track axons in milos style for {key}, error: {e}")

            # Track Axons
            try: 
                gtr.track_axons()
            except Exception as e:
                if logger: logger.info(f"unit {unit_id}_{suffix} failed to track axons for {key}, error: {e}")

            # Template Propagation - show signal propagation through channels selected
            try:
                title = key
                fig_dir = plot_dir / f"{title}_template_propagation_unit_{unit_id}.png"
                plot_template_propagation(gtr, fig_dir, unit_id, title=key, figsize = (200, 50), 
                                          fresh_plots=True, linewidth=0.1, markersize=0.1, color_marker = 'r', markerfacecolor='none', 
                                          dpi=1000, max_overlapped_lines=10)
                paths_files_generated.append(fig_dir)
            except Exception as e:
                if logger: 
                    logger.info(f"unit {unit_id}_{suffix} failed to plot template propagation for {key}, error: {e}")
                plt.close()

            # Axon Reconstruction Heuristics
            try:
                generate_axon_reconstruction_heuristics(gtr, plot_dir, unit_id, suffix=f"_{suffix}", fresh_plots=True, figsize=(20, 10)) #Sravya
                paths_files_generated.append(plot_dir / f"axon_reconstruction_heuristics{suffix}.png")
            except Exception as e:
                if logger: 
                    logger.info(f"unit {unit_id}_{suffix} failed to generate axon_reconstruction_heuristics for {key}, error: {e}")
                plt.close()

            # Axon Reconstruction Raw
            try:
                generate_axon_reconstruction_raw(gtr, plot_dir, unit_id, recon_dir, successful_recons, suffix=f"_{suffix}", fresh_plots=True, plot_full_template=False) #Sravya
                generate_axon_reconstruction_raw(gtr, plot_dir, unit_id, recon_dir, successful_recons, suffix=f"_{suffix}_full", fresh_plots=True, plot_full_template=True) #Sravya
                paths_files_generated.append(plot_dir / f"axon_reconstruction_raw{suffix}.png")
                paths_files_generated.append(plot_dir / f"axon_reconstruction_raw{suffix}_full.png")
            except Exception as e:
                if logger: 
                    logger.info(f"unit {unit_id}_{suffix} failed to generate raw axon_reconstruction for {key}, error: {e}")
                plt.close()

            # Axon Reconstruction Clean
            try:
                generate_axon_reconstruction_clean(gtr, plot_dir, unit_id, recon_dir, successful_recons, suffix=f"_{suffix}", fresh_plots=True, plot_full_template=False) #Sravya
                generate_axon_reconstruction_clean(gtr, plot_dir, unit_id, recon_dir, successful_recons, suffix=f"_{suffix}_full", fresh_plots=True, plot_full_template=True) #Sravya
                paths_files_generated.append(plot_dir / f"axon_reconstruction_clean{suffix}.png")
                paths_files_generated.append(plot_dir / f"axon_reconstruction_clean{suffix}_full.png")
            except Exception as e:
                if logger: 
                    logger.info(f"unit {unit_id}_{suffix} failed to generate clean axon_reconstruction for {key}, error: {e}")
                plt.close()

            # Axon Reconstruction Velocities
            try:
                generate_axon_reconstruction_velocities(gtr, plot_dir, unit_id, suffix=f"_{suffix}") #@Sravya want
                paths_files_generated.append(plot_dir / f"axon_reconstruction_velocities{suffix}.png")
            except Exception as e:
                if logger: 
                    logger.info(f"unit {unit_id}_{suffix} failed to generate axon_reconstruction_velocities for {key}, error: {e}")
                plt.close()

            switch = False
            if switch:
                try:
                    assert analysis_options['generate_animation'] == True, "generate_animation set to False. Skipping template map animation generation."
                    kwargs = {
                        'template': transformed_template_filled,
                        'locations': trans_loc_filled,
                        'gtr': gtr,
                        'elec_size':8, 
                        'cmap':'viridis',
                        'log': False,
                        'save_path': f"{plot_dir}/template_map_{unit_id}_{suffix}.gif",
                    }
                    play_template_map_wrapper(**kwargs)
                    paths_files_generated.append(kwargs['save_path'])
                except Exception as e:
                    if logger: 
                        logger.info(f"unit {unit_id}_{suffix} failed to play template map for {key}, error: {e}")
                    plt.close()

            try:
                assert len(gtr.branches)>0, f"No axon branches found, deleting all generated plots for unit {unit_id}_{suffix}"
                generate_axon_analytics(gtr, 
                                        analytics_data[key]['units'], 
                                        analytics_data[key]['branch_ids'], 
                                        analytics_data[key]['velocities'], 
                                        analytics_data[key]['path_lengths'], 
                                        analytics_data[key]['r2s'], 
                                        analytics_data[key]['extremums'],
                                        analytics_data[key]['init_chans'],
                                        analytics_data[key]['num_channels_included'],
                                        analytics_data[key]['channel_density'], 
                                        unit_id, 
                                        transformed_template, 
                                        trans_loc)
            except Exception as e:
                if logger: 
                    logger.info(f"unit {unit_id}_{suffix} failed to generate axon_analytics for {key}, error: {e}")
                for file in paths_files_generated:
                    try: Path(file).unlink()
                    except: pass
                analytics_data[key] = e

    return analytics_data

def analyze_and_reconstruct(templates, params=None, analysis_options=None, recon_dir=None, stream_select=None, n_jobs=8, logger=None, unit_select=None, **recon_kwargs):
    if params is None:
        params = get_default_graph_velocity_params()

    def process_result(result, all_analytics):
        for key, data in result.items():
            if isinstance(data, Exception): continue
            for metric, values in data.items():
                all_analytics[key][metric].append(values)

    def save_results(stream_id, all_analytics, recon_dir, logger=None):
        for key, data in all_analytics.items():
            suffix = key.split('_')[-1] if '_' in key else 'template'
            if 'dvdt' in key: suffix = 'dvdt'
            save_axon_analytics(
                stream_id, data['units'], data['extremums'], data['branch_ids'], data['velocities'], 
                data['path_lengths'], data['r2s'], data['num_channels_included'], data['channel_density'], data['init_chans'],
                recon_dir, suffix=f"_{suffix}"
            )

    def mark_for_deletion(unit_id, plot_dir, failed_recons):
        for suffix in failed_recons[unit_id]:
            for ext in ['png', 'svg', 'jpg']:
                for file in plot_dir.glob(f"*{suffix}*.{ext}"):
                    file.unlink()

    # Accessing the 'streams' dictionary directly
    streams = templates['streams']
    
    for stream_id, stream_templates in streams.items():
        if stream_select is not None and stream_id != f'well{stream_select:03}':
            continue

        unit_templates = stream_templates['units']
        date = templates['date']
        chip_id = templates['chip_id']
        scanType = templates['scanType']
        run_id = templates['run_id']
        recon_dir = Path(recon_dir) / date / chip_id / scanType / run_id / stream_id
        successful_recons = {str(recon_dir): {"successful_units": {}}}
        failed_recons = {}

        if logger:
            logger.info(f'Processing {len(unit_templates)} units, with {n_jobs} workers')

        all_analytics = {key: {metric: [] for metric in [
            'units', 'branch_ids', 'velocities', 'path_lengths', 'r2s', 'extremums', 
            'num_channels_included', 'channel_density', 'init_chans'
        ]} for key in [
            'merged_template', 'dvdt_merged_template', 
        ]}

        if unit_select is not None: n_jobs = 1
        if n_jobs > 1:
            with ProcessPoolExecutor(max_workers=n_jobs) as executor:
                futures = [
                    executor.submit(process_unit, unit_id, unit_templates, recon_dir, params, analysis_options, successful_recons, failed_recons, logger)
                    for unit_id, unit_templates in unit_templates.items()
                ]
                for future in as_completed(futures):
                    try:
                        result = future.result()
                        if result:
                            process_result(result, all_analytics)
                    except Exception as exc:
                        if logger:
                            logger.info(f'Unit generated an exception: {exc}')
                        if 'A process in the process pool was terminated abruptly' in str(exc):
                            raise exc
        else:
            for unit_id, unit_templates in unit_templates.items():
                if unit_select is not None and unit_id != unit_select: continue
                result = process_unit(unit_id, unit_templates, recon_dir, params, analysis_options, successful_recons, failed_recons, logger=logger)
                if result:
                    process_result(result, all_analytics)

        save_results(stream_id, all_analytics, recon_dir, logger=logger)
