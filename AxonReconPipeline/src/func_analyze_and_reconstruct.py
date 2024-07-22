'''process each unit template, generate axon reconstructions and analytics'''

from lib_axon_velocity_functions import *
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
from axon_velocity import GraphAxonTracking
from axon_velocity import get_default_graph_velocity_params

def process_unit(unit_id, unit_templates, recon_dir, params, successful_recons, failed_recons, logger=None):
    if logger is not None: 
        logger.info(f'Processing unit {unit_id}')
    else: 
        print(f'Processing unit {unit_id}')

    templates = {
        'merged_template': unit_templates['merged_template'],
        'dvdt_merged_template': unit_templates['dvdt_merged_template'],
        #'merged_template_filled': unit_templates['merged_template_filled'],
        #'dvdt_merged_template_filled': unit_templates['dvdt_merged_template_filled']
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
        'num_channels_included': [], 'channel_density': []
    } for key in templates.keys()}

    x_coords = unit_templates['merged_channel_loc'][:, 0]
    y_coords = unit_templates['merged_channel_loc'][:, 1]
    width = np.max(x_coords) - np.min(x_coords)
    height = np.max(y_coords) - np.min(y_coords)
    area = width * height
    channel_density_value = len(unit_templates['merged_channel_loc']) / area  # channels / um^2

    failed_recons[unit_id] = []

    for key, (transformed_template, transformed_template_filled, trans_loc, trans_loc_filled) in transformed_templates.items():
        suffix = key.split('_')[-1] if '_' in key else 'template'
        try:
            gtr = GraphAxonTracking(transformed_template, trans_loc, 10000, **params)
            gtr_filled = GraphAxonTracking(transformed_template_filled, trans_loc_filled, 10000, **params)
        except Exception as e:
            if logger: 
                logger.info(f"unit {unit_id} failed to initialize GraphAxonTracking for {key}, error: {e}")
            gtr = None
            gtr_filled = None

        if gtr:
            analytics_data[key]['num_channels_included'].append(len(unit_templates['merged_channel_loc']))
            analytics_data[key]['channel_density'].append(channel_density_value)

            try:
                #generate_amplitude_map(gtr, plot_dir, title=key)
                generate_amplitude_map(transformed_template_filled, trans_loc_filled, plot_dir, title=f'{key}_filled')
            except Exception as e:
                if logger: 
                    logger.info(f"unit {unit_id} failed to generate amplitude_map for {key}, error: {e}")
                plt.close()

            try:
                #generate_peak_latency_map(gtr_filled, plot_dir, title=key)
                #generate_peak_latency_map(gtr, plot_dir, title=key)
                generate_peak_latency_map(gtr_filled, plot_dir, title=f'{key}_filled')
            except Exception as e:
                if logger: 
                    logger.info(f"unit {unit_id} failed to generate peak_latency_map for {key}, error: {e}")
                plt.close()

            try:
                plot_selected_channels(gtr_filled, plot_dir, suffix=f"_{suffix}")
                gtr.track_axons()
            except Exception as e:
                if logger: 
                    logger.info(f"unit {unit_id} failed to select and plot channels or track axons for {key}, error: {e}")
                plt.close()

            try:
                generate_axon_analytics(gtr, 
                                        analytics_data[key]['units'], 
                                        analytics_data[key]['branch_ids'], 
                                        analytics_data[key]['velocities'], 
                                        analytics_data[key]['path_lengths'], 
                                        analytics_data[key]['r2s'], 
                                        analytics_data[key]['extremums'], 
                                        unit_id, transformed_template, trans_loc)
            except Exception as e:
                if logger: 
                    logger.info(f"unit {unit_id} failed to generate axon_analytics for {key}, error: {e}")

            try:
                generate_axon_reconstruction_heuristics(gtr, plot_dir, unit_id, suffix=f"_{suffix}")
                generate_axon_reconstruction_heuristics(gtr_filled, plot_dir, unit_id, suffix=f"_{suffix}_filled")
            except Exception as e:
                if logger: 
                    logger.info(f"unit {unit_id} failed to generate axon_reconstruction_heuristics for {key}, error: {e}")
                plt.close()

            try:
                generate_axon_reconstruction_raw(gtr, plot_dir, unit_id, recon_dir, successful_recons, suffix=f"_{suffix}")
            except Exception as e:
                if logger: 
                    logger.info(f"unit {unit_id} failed to generate raw axon_reconstruction for {key}, error: {e}")
                failed_recons[unit_id].append(f"raw_{suffix}")
                plt.close()

            try:
                generate_axon_reconstruction_clean(gtr, plot_dir, unit_id, recon_dir, successful_recons, suffix=f"_{suffix}")
            except Exception as e:
                if logger: 
                    logger.info(f"unit {unit_id} failed to generate clean axon_reconstruction for {key}, error: {e}")
                failed_recons[unit_id].append(f"clean_{suffix}")
                plt.close()

            try:
                generate_axon_reconstruction_velocities(gtr, plot_dir, unit_id, suffix=f"_{suffix}")
            except Exception as e:
                if logger: 
                    logger.info(f"unit {unit_id} failed to generate axon_reconstruction_velocities for {key}, error: {e}")
                plt.close()

            try:
                plot_template_propagation(gtr, plot_dir, unit_id, title=key)
            except Exception as e:
                if logger: 
                    logger.info(f"unit {unit_id} failed to plot template propagation for {key}, error: {e}")
                plt.close()

    return analytics_data

def analyze_and_reconstruct(templates, params=None, recon_dir=None, stream_select=None, n_jobs=8, logger=None):
    if params is None:
        params = get_default_graph_velocity_params()

    def process_result(result, all_analytics):
        for key, data in result.items():
            for metric, values in data.items():
                all_analytics[key][metric].extend(values)

    def save_results(stream_id, all_analytics, recon_dir):
        for key, data in all_analytics.items():
            suffix = key.split('_')[-1] if '_' in key else 'template'
            save_axon_analytics(
                stream_id, data['units'], data['extremums'], data['branch_ids'], data['velocities'], 
                data['path_lengths'], data['r2s'], data['num_channels_included'], data['channel_density'], 
                recon_dir, suffix=f"_{suffix}"
            )

    def mark_for_deletion(unit_id, plot_dir, failed_recons):
        for suffix in failed_recons[unit_id]:
            for ext in ['png', 'svg', 'jpg']:
                for file in plot_dir.glob(f"*{suffix}*.{ext}"):
                    file.unlink()
    
    for key, tmps in templates.items():
        for stream_id, stream_templates in tmps['streams'].items():
            if stream_select is not None and stream_id != f'well{stream_select:03}':
                continue
            
            unit_templates = stream_templates['units']
            date = tmps['date']
            chip_id = tmps['chip_id']
            scanType = tmps['scanType']
            run_id = tmps['run_id']
            recon_dir = Path(recon_dir) / date / chip_id / scanType / run_id / stream_id
            successful_recons = {str(recon_dir): {"successful_units": {}}}
            failed_recons = {}
            
            logger.info(f'Processing {len(unit_templates)} units, with {n_jobs} workers')

            # Initialize data containers
            all_analytics = {key: {metric: [] for metric in [
                'units', 'branch_ids', 'velocities', 'path_lengths', 'r2s', 'extremums', 
                'num_channels_included', 'channel_density'
            ]} for key in [
                'merged_template', 'dvdt_merged_template', 
                #'merged_template_filled', 'dvdt_merged_template_filled'
            ]}
            
            n_jobs = 1 #debug setting            
            if n_jobs > 1:
                with ProcessPoolExecutor(max_workers=n_jobs) as executor:
                    futures = [
                        executor.submit(process_unit, unit_id, unit_templates, recon_dir, params, successful_recons, failed_recons, logger)
                        for unit_id, unit_templates in unit_templates.items()
                    ]
                    for future in as_completed(futures):
                        try:
                            result = future.result()
                            if result:
                                process_result(result, all_analytics)
                        except Exception as exc:
                            logger.info(f'Unit generated an exception: {exc}')
                            if 'A process in the process pool was terminated abruptly' in str(exc):
                                raise exc
            else:
                for unit_id, unit_templates in unit_templates.items():
                    result = process_unit(unit_id, unit_templates, recon_dir, params, successful_recons, failed_recons, logger=logger)
                    if result:
                        process_result(result, all_analytics)
                    if unit_id in failed_recons and failed_recons[unit_id]:
                        mark_for_deletion(unit_id, recon_dir, failed_recons)

            save_results(stream_id, all_analytics, recon_dir)
