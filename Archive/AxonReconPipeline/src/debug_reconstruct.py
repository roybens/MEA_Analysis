from pathlib import Path
import json
import numpy
import os

#TODO: Use this script to debug and develop approximation of milos tracking methods using axon_velocity

toy_data_dir = "/home/adam/workspace/git_workspace/MEA_Analysis/AxonReconPipeline/toy_data_for_development/" #Note: data comes from 240322/M07037/AxonTracking/000025/well000/unit_0/
#toy_dvdt_template_dir = "/home/adam/workspace/git_workspace/MEA_Analysis/AxonReconPipeline/toy_data_for_development/"
output_dir = 'AxonReconPipeline/toy_data_for_development/testing_results/axon_tracking_plots'
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
    'distance_exp': 1.5, 'edge_dist_amp_ratio': 0.2, 'min_path_length': 80.0, 'min_path_points': 3, 'neighbor_radius': 80.0, 
    'min_points_after_branching': 2, 'mad_threshold': 10.0, 'split_paths': True, 'max_peak_latency_for_splitting': 0.7, 'r2_threshold': 0.85, 
    'r2_threshold_for_outliers': 0.95, 'min_outlier_tracking_error': 40.0}
toy_template = numpy.load(toy_data_dir + '/template.npy')
toy_template_dvdt = numpy.load(toy_data_dir + '/template_dvdt.npy')
toy_locs = numpy.load(toy_data_dir + '/locations.npy')

#Milos recon debugging
from func_analyze_and_reconstruct import approx_milos_tracking, transform_data
transformed_template, _, trans_loc, _ = transform_data(toy_template, toy_locs)
transformed_template_dvdt, _, trans_loc, _ = transform_data(toy_template_dvdt, toy_locs)
approx_milos_tracking(transformed_template_dvdt, trans_loc, input_params, output_dir, fresh_plots=True)

#Channel density analysis
#basic objecitves:
#1. Reconstruct axons and plot them
#2. Collect a few examples of low and high channel density reconstructions, put them on slides to qualitatively compare

from func_analyze_and_reconstruct import analyze_and_reconstruct
#update kwargs as needed
kwargs = {
    'n_jobs': 4,
}
#this probably wont work right away, need to debug
analyze_and_reconstruct(toy_template, toy_locs, output_dir, **kwargs)