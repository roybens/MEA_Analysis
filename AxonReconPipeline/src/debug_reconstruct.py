from pathlib import Path
import numpy as np
import os

toy_data_dir = r"C:\Users\sravy\OneDrive\Documents\MEA_Analysis\MEA_Analysis\AxonReconPipeline\toy_data_for_development\240402_M07037_well003_HOM_KCNT1_Templates"

# Output directory
output_dir = r'C:\Users\sravy\OneDrive\Documents\MEA_Analysis\MEA_Analysis\AxonReconPipeline\toy_data_for_development\testing_results\axon_tracking_plots'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
output_dir = Path(output_dir)

input_params = {
    'upsample': 2, 
    'min_selected_points': 20, 
    'verbose': False, 
    'detect_threshold': 0.01, 
    'detection_type': 'relative', 
    'kurt_threshold': 0.2, 
    'peak_std_threshold': 0.5, 
    'init_delay': 0.05, 
    'peak_std_distance': 20.0, 
    'remove_isolated': False, 
    'init_amp_peak_ratio': 0.2, 
    'max_distance_for_edge': 150.0, 
    'max_distance_to_init': 300.0, 
    'n_neighbors': 5, 
    'distance_exp': 1.5, 
    'edge_dist_amp_ratio': 0.2, 
    'min_path_length': 80.0, 
    'min_path_points': 3, 
    'neighbor_radius': 80.0, 
    'min_points_after_branching': 2, 
    'mad_threshold': 10.0, 
    'split_paths': True, 
    'max_peak_latency_for_splitting': 0.7, 
    'r2_threshold': 0.85, 
    'r2_threshold_for_outliers': 0.95, 
    'min_outlier_tracking_error': 40.0
}

# List of specific units to process, add more later 
units_to_process = [2, 3, 8, 14]

# Prepare the templates dictionary for the selected units
templates = {
    'date': '240322',  # Date of the experiment
    'chip_id': 'M07037',  # Chip ID
    'scanType': 'AxonTracking',  # Scan type
    'run_id': '000025',  # Run ID
    'streams': {
        'well000': {  # Well ID
            'units': {}
        }
    }
}

# Load data for each selected unit and update the templates dictionary
for unit in units_to_process:
    try:
        template = np.load(toy_data_dir + f'/{unit}.npy')
        template_dvdt = np.load(toy_data_dir + f'/{unit}_dvdt.npy')
        locations = np.load(toy_data_dir + f'/{unit}_channels.npy')
        template_filled = np.load(toy_data_dir + f'/{unit}_filled.npy')
        locations_filled = np.load(toy_data_dir + f'/{unit}_channels_filled.npy')

        # Update the templates dictionary with the loaded data
        templates['streams']['well000']['units'][unit] = {
            'merged_template': template,
            'dvdt_merged_template': template_dvdt,
            'merged_channel_loc': locations,
            'merged_template_filled': template_filled,
            'merged_channel_locs_filled': locations_filled
        }
    except FileNotFoundError:
        print(f"Files for unit {unit} not found. Skipping this unit.")

from func_analyze_and_reconstruct import analyze_and_reconstruct

if __name__ == '__main__':
    # Update kwargs as needed
    kwargs = {
        'n_jobs': 4,
    }
    analyze_and_reconstruct(templates, params=input_params, recon_dir=output_dir, **kwargs)
