'''process each unit template, generate axon reconstructions and analytics'''

from lib_axon_velocity_functions import *

#parallel processing
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

def process_unit(unit_id, unit_templates, recon_dir, params, successful_recons, logger=None):
    if logger is not None:
        logger.info(f'Processing unit {unit_id}')
    else:
        print(f'Processing unit {unit_id}')
    
    merged_template = unit_templates['merged_template']
    dvdt_merged_template = unit_templates['dvdt_merged_template']
    merged_channel_loc = unit_templates['merged_channel_loc']
    merged_template_filled = unit_templates['merged_template_filled']
    dvdt_merged_template_filled = unit_templates['dvdt_merged_template_filled']
    merged_channel_filled_loc = unit_templates['merged_channel_locs_filled']

    if merged_template is None:
        return None  # skip if no merged template found
    
    transformed_template, transformed_template_filled, trans_loc, trans_loc_filled = transform_data(merged_template, merged_channel_loc, merged_template_filled, merged_channel_filled_loc)
    transformed_dvdt_template, transformed_dvdt_template_filled, trans_loc_dvdt, trans_loc_filled_dvdt = transform_data(dvdt_merged_template, merged_channel_loc, dvdt_merged_template_filled, merged_channel_filled_loc)
    
    plot_dir = create_plot_dir(recon_dir, unit_id)
    units, branch_ids, velocities, path_lengths, r2s, extremums, num_channels_included, channel_density = [], [], [], [], [], [], [], []
    dvdt_units, dvdt_branch_ids, dvdt_velocities, dvdt_path_lengths, dvdt_r2s, dvdt_extremums, dvdt_num_channels_included, dvdt_channel_density = [], [], [], [], [], [], [], []

    num_channels_included.append(len(merged_channel_loc))
    dvdt_num_channels_included.append(len(merged_channel_loc))

    x_coords = merged_channel_loc[:, 0]
    y_coords = merged_channel_loc[:, 1]
    width = np.max(x_coords) - np.min(x_coords)
    height = np.max(y_coords) - np.min(y_coords)
    area = width * height
    channel_density.append(len(merged_channel_loc) / area)  # channels / um^2
    dvdt_channel_density.append(len(merged_channel_loc) / area)  # channels / um^2
    
    try:
        generate_amplitude_map(transformed_template_filled, trans_loc_filled, plot_dir, title="merged_template")
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate amplitude_map for merged_template, error: {e}")
        plt.close()

    try:
        generate_amplitude_map(transformed_dvdt_template_filled, trans_loc_filled_dvdt, plot_dir, title="dvdt_merged_template")
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate amplitude_map for dvdt_merged_template, error: {e}")
        plt.close()
    
    try:
        generate_peak_latency_map(transformed_template_filled, trans_loc_filled, plot_dir, title="merged_template")
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate peak_latency_map for merged_template, error: {e}")
        plt.close()

    try:
        generate_peak_latency_map(transformed_dvdt_template_filled, trans_loc_filled_dvdt, plot_dir, title="dvdt_merged_template")
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate peak_latency_map for dvdt_merged_template, error: {e}")
        plt.close()
    
    try:
        gtr0 = av.GraphAxonTracking(transformed_template, trans_loc, 10000, **params)
        select_and_plot_channels(gtr0, plot_dir)
        gtr0.track_axons()
        #gtr0 = av.compute_graph_propagation_velocity(transformed_template, trans_loc, 10000, **params)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to select and plot channels or compute_graph_propagation_velocity for merged_template, error: {e}")
        plt.close()
        return None
    
    try:
        gtr1 = av.GraphAxonTracking(transformed_dvdt_template, trans_loc_dvdt, 10000, **params)
        select_and_plot_channels(gtr1, plot_dir, suffix="_dvdt")
        gtr1.track_axons()
        #gtr1 = av.compute_graph_propagation_velocity(transformed_dvdt_template, trans_loc_dvdt, 10000, **params)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to select and plot channels or compute_graph_propagation_velocity for dvdt_merged_template, error: {e}")
        plt.close()
        return None
    
    try:
        generate_axon_analytics(gtr0, units, branch_ids, velocities, path_lengths, r2s, extremums, unit_id, transformed_template, trans_loc)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate axon_analytics for merged_template, error: {e}")

    try:
        generate_axon_analytics_dvdt(gtr1, dvdt_units, dvdt_branch_ids, dvdt_velocities, dvdt_path_lengths, dvdt_r2s, dvdt_extremums, unit_id, transformed_dvdt_template, trans_loc_dvdt)
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate axon_analytics for dvdt_merged_template, error: {e}")
    
    try:
        generate_axon_reconstruction_heuristics(gtr0, plot_dir, unit_id, suffix="_merged_template")
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate axon_reconstruction_heuristics for merged_template, error: {e}")
        plt.close()  # close the figure, avoid memory leak

    try:
        generate_axon_reconstruction_heuristics(gtr1, plot_dir, unit_id, suffix="_dvdt_merged_template")
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate axon_reconstruction_heuristics for dvdt_merged_template, error: {e}")
        plt.close()  # close the figure, avoid memory leak
    
    try:
        generate_axon_reconstruction_raw(gtr0, plot_dir, unit_id, recon_dir, successful_recons, suffix="_merged_template")
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate raw axon_reconstruction for merged_template, error: {e}")
        plt.close()  # close the figure, avoid memory leak

    try:
        generate_axon_reconstruction_raw(gtr1, plot_dir, unit_id, recon_dir, successful_recons, suffix="_dvdt_merged_template")
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate raw axon_reconstruction for dvdt_merged_template, error: {e}")
        plt.close()  # close the figure, avoid memory leak
    
    try:
        generate_axon_reconstruction_clean(gtr0, plot_dir, unit_id, recon_dir, successful_recons, suffix="_merged_template")
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate clean axon_reconstruction for merged_template, error: {e}")
        plt.close()  # close the figure, avoid memory leak

    try:
        generate_axon_reconstruction_clean(gtr1, plot_dir, unit_id, recon_dir, successful_recons, suffix="_dvdt_merged_template")
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate clean axon_reconstruction for dvdt_merged_template, error: {e}")
        plt.close()  # close the figure, avoid memory leak
    
    try:
        generate_axon_reconstruction_velocities(gtr0, plot_dir, unit_id, suffix="_merged_template")
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate axon_reconstruction_velocities for merged_template, error: {e}")
        plt.close()  # close the figure, avoid memory leak

    try:
        generate_axon_reconstruction_velocities(gtr1, plot_dir, unit_id, suffix="_dvdt_merged_template")
    except Exception as e:
        logger.info(f"unit {unit_id} failed to generate axon_reconstruction_velocities for dvdt_merged_template, error: {e}")
        plt.close()  # close the figure, avoid memory leak
    
    try:
        plot_template_propagation_wrapper(transformed_template, trans_loc, gtr0.selected_channels, plot_dir, unit_id, title="merged_template")
    except Exception as e:
        logger.info(f"unit {unit_id} failed to plot template propagation for merged_template, error: {e}")
        plt.close()  # close the figure, avoid memory leak

    try:
        plot_template_propagation_wrapper(transformed_dvdt_template, trans_loc_dvdt, gtr1.selected_channels, plot_dir, unit_id, title="dvdt_merged_template")
    except Exception as e:
        logger.info(f"unit {unit_id} failed to plot template propagation for dvdt_merged_template, error: {e}")
        plt.close()  # close the figure, avoid memory leak

    return (units, extremums, branch_ids, velocities, path_lengths, 
            r2s, num_channels_included, channel_density, 
            dvdt_units, dvdt_extremums, dvdt_branch_ids, dvdt_velocities, 
            dvdt_path_lengths, dvdt_r2s, dvdt_num_channels_included, 
            dvdt_channel_density)



def analyze_and_reconstruct(templates, params=av.get_default_graph_velocity_params(), recon_dir=None, stream_select=None, n_jobs=8, logger=None):
    #aw21July2024 - TODO: need to verify units of dv/dt plots. Probably important.
    
    def process_result(result, all_data, all_dvdt_data):
        (
            units, extremums, branch_ids, velocities, path_lengths, r2s, num_channels_included, channel_density, 
            dvdt_units, dvdt_extremums, dvdt_branch_ids, dvdt_velocities, dvdt_path_lengths, dvdt_r2s, 
            dvdt_num_channels_included, dvdt_channel_density
        ) = result
        
        all_data['units'].append(units)
        all_data['branch_ids'].append(branch_ids)
        all_data['velocities'].append(velocities)
        all_data['path_lengths'].append(path_lengths)
        all_data['r2s'].append(r2s)
        all_data['extremums'].append(extremums)
        all_data['num_channels_included'].append(num_channels_included)
        all_data['channel_densities'].append(channel_density)
        
        all_dvdt_data['units'].append(dvdt_units)
        all_dvdt_data['branch_ids'].append(dvdt_branch_ids)
        all_dvdt_data['velocities'].append(dvdt_velocities)
        all_dvdt_data['path_lengths'].append(dvdt_path_lengths)
        all_dvdt_data['r2s'].append(dvdt_r2s)
        all_dvdt_data['extremums'].append(dvdt_extremums)
        all_dvdt_data['num_channels_included'].append(dvdt_num_channels_included)
        all_dvdt_data['channel_densities'].append(dvdt_channel_density)

    def save_results(stream_id, all_data, all_dvdt_data, recon_dir):
        save_axon_analytics(stream_id, all_data['units'], all_data['extremums'], all_data['branch_ids'], all_data['velocities'], all_data['path_lengths'], all_data['r2s'], all_data['num_channels_included'], all_data['channel_densities'], recon_dir)
        save_axon_analytics(stream_id, all_dvdt_data['units'], all_dvdt_data['extremums'], all_dvdt_data['branch_ids'], all_dvdt_data['velocities'], all_dvdt_data['path_lengths'], all_dvdt_data['r2s'], all_dvdt_data['num_channels_included'], all_dvdt_data['channel_densities'], recon_dir, suffix="_dvdt")

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
            
            logger.info(f'Processing {len(unit_templates)} units, with {n_jobs} workers')

            # Initialize data containers
            all_data = {key: [] for key in ['units', 'branch_ids', 'velocities', 'path_lengths', 'r2s', 'extremums', 'num_channels_included', 'channel_densities']}
            all_dvdt_data = {key: [] for key in ['units', 'branch_ids', 'velocities', 'path_lengths', 'r2s', 'extremums', 'num_channels_included', 'channel_densities']}
            
            concurrent_bool = False  # Debug setting
            n_jobs = 1  # Debug setting
            
            if concurrent_bool:
                with ProcessPoolExecutor(max_workers=n_jobs) as executor:
                    futures = [
                        executor.submit(process_unit, unit_id, unit_templates, recon_dir, params, successful_recons, logger)
                        for unit_id, unit_templates in unit_templates.items()
                    ]
                    for future in as_completed(futures):
                        try:
                            result = future.result()
                            if result:
                                process_result(result, all_data, all_dvdt_data)
                        except Exception as exc:
                            logger.info(f'Unit generated an exception: {exc}')
                            if 'A process in the process pool was terminated abruptly' in str(exc):
                                raise exc
            else:
                for unit_id, unit_templates in unit_templates.items():
                    result = process_unit(unit_id, unit_templates, recon_dir, params, successful_recons, logger=logger)
                    if result:
                        process_result(result, all_data, all_dvdt_data)

            save_results(stream_id, all_data, all_dvdt_data, recon_dir)