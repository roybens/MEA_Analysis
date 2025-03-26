import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from scipy.spatial.distance import mahalanobis

from sklearn.experimental import enable_iterative_imputer  # Required for IterativeImputer
from sklearn.impute import IterativeImputer
import traceback

import spikeinterface.postprocessing as spost

def check_features(feature_dict):
    """Check feature values for issues and report problematic feature names."""
    for feature_name, feature_value in feature_dict.items():
        if feature_value is None:
            print(f"WARNING: Feature '{feature_name}' is None")
        elif feature_value is np.nan:
            print(f"WARNING: Feature '{feature_name}' is np.nan")
        elif isinstance(feature_value, list) and len(feature_value) == 0:
            print(f"WARNING: Feature '{feature_name}' is an empty list")
        elif isinstance(feature_value, dict) and len(feature_value) == 0:
            print(f"WARNING: Feature '{feature_name}' is an empty dict")
        elif isinstance(feature_value, np.ndarray) and len(feature_value) == 0:
            print(f"WARNING: Feature '{feature_name}' is an empty ndarray")
        elif isinstance(feature_value, np.ndarray) and np.isnan(feature_value).all():
            print(f"WARNING: Feature '{feature_name}' is all NaNs")
        elif isinstance(feature_value, np.ndarray) and np.isinf(feature_value).all():
            print(f"WARNING: Feature '{feature_name}' is all Infs")
        elif isinstance(feature_value, np.ndarray) and np.isnan(feature_value).any():
            print(f"WARNING: Feature '{feature_name}' contains NaNs")
        elif isinstance(feature_value, list) and len(feature_value) > 0 and any(x is None for x in feature_value):
            print(f"WARNING: Feature '{feature_name}' contains None values")
        elif isinstance(feature_value, list) and len(feature_value) > 0 and any(x is np.nan for x in feature_value):
            print(f"WARNING: Feature '{feature_name}' contains np.nan values")
        elif isinstance(feature_value, list) and len(feature_value) > 0 and any(x == '' for x in feature_value):
            print(f"WARNING: Feature '{feature_name}' contains empty strings")

def extract_features_v2(we, unit_ids, network_data, sorting_object, sampling_frequency=None):
    """
    Extracts neuronal features related to bursting and firing patterns.
    """
    # init
    if sampling_frequency is None:
        sampling_frequency = sorting_object.get_sampling_frequency()
    
    # include/exclude units
    include_units = []
    exclude_units = {}
    exclusion_criteria = {
        # very low firing rate, trying 0.1 Hz for now # aw 2025-02-18 20:00:37
        # aw 2025-02-18 20:23:03 - firing rate alone is very good at excluding units and make the clustering more clean and as expected
        'fr': 0.1,
        'num_spikes': 50 # also just trying a minimum number of spikes and seeing how it goes. 
    }
    features = []
    feature_dicts = {}
    
    for unit_id in unit_ids:      
        
        try:
            
            ## ** spiking features **
            spiking_data = network_data['spiking_data']['spiking_metrics_by_unit'].get(unit_id, {})
            if len(spiking_data) == 0:
                exclude_units[unit_id] = {
                    'reason': 'no_spiking_data'
                    }
                continue
            
            num_spikes = spiking_data['num_spikes']
            fr = spiking_data.get('fr', np.nan)
            isi_mean = spiking_data['isi'].get('mean', np.nan)
            isi_std = spiking_data['isi'].get('std', np.nan)
            isi_cov = spiking_data['isi'].get('cov', np.nan)
            isi_median = spiking_data['isi'].get('median', np.nan)
            
            ### exclude very quiet units and/or too few spikes
            if fr < exclusion_criteria['fr']:
                exclude_units[unit_id] = {
                    'reason': 'low_fr',
                    'fr': fr,
                    'num_spikes': num_spikes
                    }
                continue
            
            if num_spikes < exclusion_criteria['num_spikes']:
                exclude_units[unit_id] = {
                    'reason': 'low_num_spikes',
                    'num_spikes': num_spikes,
                    'fr': fr
                    }
                continue
            
            ## ** wf features **
            wf_metrics = network_data['spiking_data']['spiking_metrics_by_unit'][unit_id]['wf_metrics']
            # weighted_var = wf_metrics.get('weighted_variability', np.nan)
            # weighted_cov = wf_metrics.get('weighted_coefficient_of_variation', np.nan)
            # aw 2025-02-18 15:06:39 - excluding above metrics because I need single values, not lists
            mean_weighted_var = wf_metrics['bio_variability_metrics'].get('mean_variance', np.nan)
            std_weighted_var = wf_metrics['bio_variability_metrics'].get('std_variance', np.nan)
            cov_weighted_var = wf_metrics['bio_variability_metrics'].get('cov_variance', np.nan)
            mean_weighted_cov = wf_metrics['bio_variability_metrics'].get('mean_cv', np.nan)
            std_weighted_cov = wf_metrics['bio_variability_metrics'].get('std_cv', np.nan)
            cov_weighted_cov = wf_metrics['bio_variability_metrics'].get('cov_cv', np.nan)
            
            #amplitude
            peak_amplitude = wf_metrics['amplitude_metrics'].get('peak_amplitude', np.nan)
            trough_amplitude = wf_metrics['amplitude_metrics'].get('trough_amplitude', np.nan)
            peak_to_trough_amplitude = wf_metrics['amplitude_metrics'].get('peak_to_trough_amplitude', np.nan)
            
            #temporal
            trough_time_ms = wf_metrics['temporal_metrics'].get('trough_time_ms', np.nan)
            peak_time_ms = wf_metrics['temporal_metrics'].get('peak_time_ms', np.nan)
            peak_to_trough_time_ms = wf_metrics['temporal_metrics'].get('peak_to_trough_time_ms', np.nan)
            ap_phase_duration_ms = wf_metrics['temporal_metrics'].get('ap_phase_duration_ms', np.nan)
            ap_start_ms = wf_metrics['temporal_metrics'].get('ap_start_ms', np.nan)
            ap_end_ms = wf_metrics['temporal_metrics'].get('ap_end_ms', np.nan)
            refractory_phase_duration_ms = wf_metrics['temporal_metrics'].get('refractory_phase_duration_ms', np.nan)
            refractory_end_ms = wf_metrics['temporal_metrics'].get('refractory_end_ms', np.nan)
            spike_width_half_max_ms = wf_metrics['temporal_metrics'].get('spike_width_half_max_ms', np.nan)
            
            #slope
            max_depolarization_slope = wf_metrics['slope_metrics'].get('max_depolarization_slope', np.nan)  
            max_repolarization_slope = wf_metrics['slope_metrics'].get('max_repolarization_slope', np.nan)
            slope_ratio = wf_metrics['slope_metrics'].get('slope_ratio', np.nan)
            
            #asymmetry
            trough_to_peak_ratio = wf_metrics['waveform_asymmetry'].get('trough_to_peak_ratio', np.nan)
            waveform_asymmetry_index = wf_metrics['waveform_asymmetry'].get('waveform_asymmetry_index', np.nan)
            
            #energy
            ap_phase_power_uv2 = wf_metrics['energy_metrics'].get('ap_phase_power_uv2', np.nan)
            refractory_phase_power_uv2 = wf_metrics['energy_metrics'].get('refractory_phase_power_uv2', np.nan)
            total_spike_power_uv2 = wf_metrics['energy_metrics'].get('total_spike_power_uv2', np.nan)
            
            #shape
            waveform_skewness = wf_metrics['waveform_shape'].get('waveform_skewness', np.nan)
            waveform_kurtosis = wf_metrics['waveform_shape'].get('waveform_kurtosis', np.nan)
            
            # vertical smearing
            #std_across_time = wf_metrics['vertical_smearing'].get('std_across_time', np.nan)
            cv_peak_amplitude = wf_metrics['vertical_smearing'].get('cv_peak_amplitude', np.nan)
            cv_trough_amplitude = wf_metrics['vertical_smearing'].get('cv_trough_amplitude', np.nan)
            peak_to_trough_variance = wf_metrics['vertical_smearing'].get('peak_to_trough_variance', np.nan)
            
            # horizontal smearing
            trough_time_std_ms = wf_metrics['horizontal_smearing'].get('trough_time_std_ms', np.nan)
            peak_time_std_ms = wf_metrics['horizontal_smearing'].get('peak_time_std_ms', np.nan)
            trough_to_peak_jitter_ms = wf_metrics['horizontal_smearing'].get('trough_to_peak_jitter_ms', np.nan)        
            
            ## ** burst features **
            burst_data = network_data['bursting_data']['unit_metrics'].get(unit_id, {})
            burst_part_rate = burst_data.get('burst_part_rate', np.nan)
            quiet_part_rate = burst_data.get('quiet_part_rate', np.nan)
            burst_part_perc = burst_data.get('burst_part_perc', np.nan)
            in_burst_fr_mean = burst_data['fr']['in_burst'].get('mean', np.nan)
            in_burst_fr_std = burst_data['fr']['in_burst'].get('std', np.nan)
            in_burst_fr_cov = burst_data['fr']['in_burst'].get('cov', np.nan)
            out_burst_fr_mean = burst_data['fr']['out_burst'].get('mean', np.nan)
            out_burst_fr_std = burst_data['fr']['out_burst'].get('std', np.nan)
            out_burst_fr_cov = burst_data['fr']['out_burst'].get('cov', np.nan)
            in_burst_isi_mean = burst_data['isi']['in_burst'].get('mean', np.nan)
            in_burst_isi_std = burst_data['isi']['in_burst'].get('std', np.nan)
            in_burst_isi_cov = burst_data['isi']['in_burst'].get('cov', np.nan)
            out_burst_isi_mean = burst_data['isi']['out_burst'].get('mean', np.nan)
            out_burst_isi_std = burst_data['isi']['out_burst'].get('std', np.nan)
            out_burst_isi_cov = burst_data['isi']['out_burst'].get('cov', np.nan)
            spike_count_in_burst_mean = burst_data['spike_counts']['in_burst'].get('mean', np.nan)
            spike_count_in_burst_std = burst_data['spike_counts']['in_burst'].get('std', np.nan)
            spike_count_in_burst_cov = burst_data['spike_counts']['in_burst'].get('cov', np.nan)
            spike_count_out_burst_mean = burst_data['spike_counts']['out_burst'].get('mean', np.nan)
            spike_count_out_burst_std = burst_data['spike_counts']['out_burst'].get('std', np.nan)
            spike_count_out_burst_cov = burst_data['spike_counts']['out_burst'].get('cov', np.nan)
            ff_in_burst = burst_data['fano_factor']['in_burst']
            ff_out_burst = burst_data['fano_factor']['out_burst']       
            
            ## ** mega burst features **
            mega_burst_data = network_data['mega_bursting_data']['unit_metrics'].get(unit_id, {})
            mega_burst_part_rate = mega_burst_data.get('burst_part_rate', np.nan)
            mega_quiet_part_rate = mega_burst_data.get('quiet_part_rate', np.nan)
            mega_burst_part_perc = mega_burst_data.get('burst_part_perc', np.nan)
            mega_in_burst_fr_mean = mega_burst_data['fr']['in_burst'].get('mean', np.nan)
            mega_in_burst_fr_std = mega_burst_data['fr']['in_burst'].get('std', np.nan)
            mega_in_burst_fr_cov = mega_burst_data['fr']['in_burst'].get('cov', np.nan)
            mega_out_burst_fr_mean = mega_burst_data['fr']['out_burst'].get('mean', np.nan)
            mega_out_burst_fr_std = mega_burst_data['fr']['out_burst'].get('std', np.nan)
            mega_out_burst_fr_cov = mega_burst_data['fr']['out_burst'].get('cov', np.nan)
            mega_in_burst_isi_mean = mega_burst_data['isi']['in_burst'].get('mean', np.nan)
            mega_in_burst_isi_std = mega_burst_data['isi']['in_burst'].get('std', np.nan)
            mega_in_burst_isi_cov = mega_burst_data['isi']['in_burst'].get('cov', np.nan)
            mega_out_burst_isi_mean = mega_burst_data['isi']['out_burst'].get('mean', np.nan)
            mega_out_burst_isi_std = mega_burst_data['isi']['out_burst'].get('std', np.nan)
            mega_out_burst_isi_cov = mega_burst_data['isi']['out_burst'].get('cov', np.nan)
            mega_spike_count_in_burst_mean = mega_burst_data['spike_counts']['in_burst'].get('mean', np.nan)
            mega_spike_count_in_burst_std = mega_burst_data['spike_counts']['in_burst'].get('std', np.nan)
            mega_spike_count_in_burst_cov = mega_burst_data['spike_counts']['in_burst'].get('cov', np.nan)
            mega_spike_count_out_burst_mean = mega_burst_data['spike_counts']['out_burst'].get('mean', np.nan)
            mega_spike_count_out_burst_std = mega_burst_data['spike_counts']['out_burst'].get('std', np.nan)
            mega_spike_count_out_burst_cov = mega_burst_data['spike_counts']['out_burst'].get('cov', np.nan)
            mega_ff_in_burst = mega_burst_data['fano_factor']['in_burst']
            mega_ff_out_burst = mega_burst_data['fano_factor']['out_burst']        
            
            # ** Store features as a dictionary **
            feature_dict = {
                # Spikes
                "fr": fr,
                "isi_mean": isi_mean,
                "isi_std": isi_std,
                "isi_cov": isi_cov,
                "isi_median": isi_median,

                # Waveforms
                # "weighted_var": weighted_var,
                # "weighted_cov": weighted_cov,
                "mean_weighted_var": mean_weighted_var,
                "std_weighted_var": std_weighted_var,
                "cov_weighted_var": cov_weighted_var,
                "mean_weighted_cov": mean_weighted_cov,
                "std_weighted_cov": std_weighted_cov,
                "cov_weighted_cov": cov_weighted_cov,
                "peak_amplitude": peak_amplitude,
                "trough_amplitude": trough_amplitude,
                "peak_to_trough_amplitude": peak_to_trough_amplitude,
                "trough_time_ms": trough_time_ms,
                "peak_time_ms": peak_time_ms,
                "peak_to_trough_time_ms": peak_to_trough_time_ms,
                "ap_phase_duration_ms": ap_phase_duration_ms,
                #"ap_start_ms": ap_start_ms,
                #"ap_end_ms": ap_end_ms,
                "refractory_phase_duration_ms": refractory_phase_duration_ms,
                #"refractory_end_ms": refractory_end_ms,
                "spike_width_half_max_ms": spike_width_half_max_ms,
                "max_depolarization_slope": max_depolarization_slope,
                "max_repolarization_slope": max_repolarization_slope,
                "slope_ratio": slope_ratio,
                "trough_to_peak_ratio": trough_to_peak_ratio,
                "waveform_asymmetry_index": waveform_asymmetry_index,
                "ap_phase_power_uv2": ap_phase_power_uv2,
                "refractory_phase_power_uv2": refractory_phase_power_uv2,
                "total_spike_power_uv2": total_spike_power_uv2,
                "waveform_skewness": waveform_skewness,
                "waveform_kurtosis": waveform_kurtosis,
                #"std_across_time": std_across_time,
                "cv_peak_amplitude": cv_peak_amplitude,
                "cv_trough_amplitude": cv_trough_amplitude,
                "peak_to_trough_variance": peak_to_trough_variance,
                "trough_time_std_ms": trough_time_std_ms,
                "peak_time_std_ms": peak_time_std_ms,
                "trough_to_peak_jitter_ms": trough_to_peak_jitter_ms,

                # Bursts
                "burst_part_rate": burst_part_rate,
                "quiet_part_rate": quiet_part_rate,
                "burst_part_perc": burst_part_perc,
                "in_burst_fr_mean": in_burst_fr_mean,
                "in_burst_fr_std": in_burst_fr_std,
                "in_burst_fr_cov": in_burst_fr_cov,
                "out_burst_fr_mean": out_burst_fr_mean,
                "out_burst_fr_std": out_burst_fr_std,
                "out_burst_fr_cov": out_burst_fr_cov,
                "in_burst_isi_mean": in_burst_isi_mean,
                "in_burst_isi_std": in_burst_isi_std,
                "in_burst_isi_cov": in_burst_isi_cov,
                "out_burst_isi_mean": out_burst_isi_mean,
                "out_burst_isi_std": out_burst_isi_std,
                "out_burst_isi_cov": out_burst_isi_cov,
                "spike_count_in_burst_mean": spike_count_in_burst_mean,
                "spike_count_in_burst_std": spike_count_in_burst_std,
                "spike_count_in_burst_cov": spike_count_in_burst_cov,
                "spike_count_out_burst_mean": spike_count_out_burst_mean,
                "spike_count_out_burst_std": spike_count_out_burst_std,
                "spike_count_out_burst_cov": spike_count_out_burst_cov,
                "ff_in_burst": ff_in_burst,
                "ff_out_burst": ff_out_burst,

                # Mega Bursts
                "mega_burst_part_rate": mega_burst_part_rate,
                "mega_quiet_part_rate": mega_quiet_part_rate,
                "mega_burst_part_perc": mega_burst_part_perc,
                "mega_in_burst_fr_mean": mega_in_burst_fr_mean,
                "mega_in_burst_fr_std": mega_in_burst_fr_std,
                "mega_in_burst_fr_cov": mega_in_burst_fr_cov,
                "mega_out_burst_fr_mean": mega_out_burst_fr_mean,
                "mega_out_burst_fr_std": mega_out_burst_fr_std,
                "mega_out_burst_fr_cov": mega_out_burst_fr_cov,
                "mega_in_burst_isi_mean": mega_in_burst_isi_mean,
                "mega_in_burst_isi_std": mega_in_burst_isi_std,
                "mega_in_burst_isi_cov": mega_in_burst_isi_cov,
                "mega_out_burst_isi_mean": mega_out_burst_isi_mean,
                "mega_out_burst_isi_std": mega_out_burst_isi_std,
                "mega_out_burst_isi_cov": mega_out_burst_isi_cov,
                "mega_spike_count_in_burst_mean": mega_spike_count_in_burst_mean,
                "mega_spike_count_in_burst_std": mega_spike_count_in_burst_std,
                "mega_spike_count_in_burst_cov": mega_spike_count_in_burst_cov,
                "mega_spike_count_out_burst_mean": mega_spike_count_out_burst_mean,
                "mega_spike_count_out_burst_std": mega_spike_count_out_burst_std,
                "mega_spike_count_out_burst_cov": mega_spike_count_out_burst_cov,
                "mega_ff_in_burst": mega_ff_in_burst,
                "mega_ff_out_burst": mega_ff_out_burst,
            }
        
            # ** Run feature check **
            check_features(feature_dict)
            
            # include/exclude units
            #exclude_bool = exclude_check(exclusion_criteria, feature_dict)
            exclude_bool = False #HACK
            if exclude_bool:
                # exclude_units.append(unit_id)
                exclude_units[unit_id] = {
                    'reason': 'exclusion_criteria',
                    'feature_dict': feature_dict
                    }
                continue
            else:
                include_units.append(unit_id)
                
            # append to feature_dicts
            feature_dicts[unit_id] = feature_dict
            
            # Convert to a list after validation
            feature_vector = list(feature_dict.values())
            
            # Store extracted features
            features.append(feature_vector)
        except Exception as e:
            print(f"Error extracting features for unit {unit_id}: {e}")
            traceback.print_exc()
            #exclude_units.append(unit_id)
            exclude_units[unit_id] = {
                'reason': 'error',
                'error': str(e)
                }
            continue
        
    print(f'{len(exclude_units)} units excluded based on feature extraction')
    print(f'{len(include_units)} units included based on feature extraction')
    return features, feature_dicts, include_units, exclude_units
    
def classify(network_data, plot_wfs=False, **kwargs):
    '''
    Classifies neurons into clusters based on multiple waveform and burst activity features.
    '''
    
    ## init
    wfs_output_path = network_data['waveform_output']
    sorting_object = kwargs['sorting_object']
    sorting_object.remove_empty_units()
    sampling_frequency = sorting_object.get_sampling_frequency()
    we = kwargs['wf_extractor']
    
    ## ** Extract Features from Spiking and Bursting Analysis **
    unit_ids = sorting_object.get_unit_ids()
    
    features, feature_dicts, include_units, exclude_units = extract_features_v2(we, unit_ids, network_data, sorting_object, sampling_frequency)    
    
    ## ** Feature Scaling and Dimensionality Reduction **
    
    # Convert to numpy array
    features = np.array(features)
    
    # change Nones to nans
    features = np.where(features == None, np.nan, features)
    
    # Create NaN indicators to preserve missing information
    nan_mask = np.isnan(features).astype(float)  # 1 for NaN, 0 for non-NaN
    features_with_nan_info = np.hstack((features, nan_mask))  # Append mask as extra features
    
    # Replace all nans with 0
    features = np.nan_to_num(features, nan=0.0)
    
    # Normalize features
    scaler = StandardScaler()
    features_scaled = scaler.fit_transform(features)
    
    # Reduce dimensionality using PCA (keep 2 components)
    explained_variance = None
    n_components = 2
    while explained_variance is None or explained_variance < 0.9:

        pca = PCA(n_components=n_components)
        features_pca = pca.fit_transform(features_scaled)
        explained_variance = np.sum(pca.explained_variance_ratio_)
        
        print(f'{n_components} components capture {explained_variance:.2f} variance')
        n_components += 1
        
    print(f'{n_components-1} components required to capture {explained_variance:.2f} variance (>= 0.9)')
    #print(f'{len(features_scaled)} features reduced to {features_pca.shape[1]} components')
    
    #check how well n_components=2 captures variance
    print(f"Explained variance ratio: {pca.explained_variance_ratio_}")
    
    #TODO: loop here, add PCA components until explained variance is > 0.9
    
    # Apply K-means clustering
    kmeans = KMeans(n_clusters=2)
    cluster_labels = kmeans.fit_predict(features_pca)
    
    # Outlier Removal Using Mahalanobis Distance
    remove_outliers = False
    if remove_outliers:
        cluster_centroid = np.mean(features_pca, axis=0)
        cov_matrix = np.cov(features_pca, rowvar=False)
        inv_cov_matrix = np.linalg.inv(cov_matrix)

        distances = np.array([mahalanobis(point, cluster_centroid, inv_cov_matrix) for point in features_pca])
        outlier_threshold = 1.4
        
        valid_units = distances < outlier_threshold
        features_pca_filtered = features_pca[valid_units]
        cluster_labels_filtered = cluster_labels[valid_units]
        
        #redefine
        cluster_labels = cluster_labels_filtered
        features_pca = features_pca_filtered
    
    # check which cluster has higher mean firing rate - assume this is the inhibitory cluster
    cluster_1_frs = []
    cluster_2_frs = []
    cluster_1_ffs = []
    cluster_2_ffs = []
    for i in range(2):
        #frs.append(np.mean(features[cluster_labels == i, 2]))
        for unit_idx, unit_id in enumerate(feature_dicts.keys()):
            designated_cluster = cluster_labels[unit_idx]
            fr = feature_dicts[unit_id]['fr']
            ff = feature_dicts[unit_id]['ff_in_burst']
            if designated_cluster == 0:
                #cluster_1_frs.append(fr)
                cluster_1_frs.append(fr)
                cluster_1_ffs.append(ff)
            elif designated_cluster == 1:
                #cluster_2_frs.append(fr)
                cluster_2_frs.append(fr)
                cluster_2_ffs.append(ff)
    
    verbose = True
    if verbose:
        
        #
        mean_fr_cluster_1 = np.mean(cluster_1_frs)
        max_fr_cluster_1 = np.max(cluster_1_frs)
        min_fr_cluster_1 = np.min(cluster_1_frs)
        cov_fr_cluster_1 = np.std(cluster_1_frs) / mean_fr_cluster_1
        # print in multiple lines
        print(f"Cluster 1 FR: {mean_fr_cluster_1:.2f} Hz")
        print(f"Max FR: {max_fr_cluster_1:.2f} Hz")
        print(f"Min FR: {min_fr_cluster_1:.2f} Hz")
        print(f"CoV FR: {cov_fr_cluster_1:.2f}")
        
        mean_ff_cluster_1 = np.mean(cluster_1_ffs)
        max_ff_cluster_1 = np.max(cluster_1_ffs)
        min_ff_cluster_1 = np.min(cluster_1_ffs)
        cov_ff_cluster_1 = np.std(cluster_1_ffs) / mean_ff_cluster_1
        # print in multiple lines
        print(f"Cluster 1 FF: {mean_ff_cluster_1:.2f}")
        print(f"Max FF: {max_ff_cluster_1:.2f}")
        print(f"Min FF: {min_ff_cluster_1:.2f}")
        print(f"CoV FF: {cov_ff_cluster_1:.2f}")
        
        mean_fr_cluster_2 = np.mean(cluster_2_frs)
        max_fr_cluster_2 = np.max(cluster_2_frs)
        min_fr_cluster_2 = np.min(cluster_2_frs)
        cov_fr_cluster_2 = np.std(cluster_2_frs) / mean_fr_cluster_2
        # print in multiple lines
        print(f"Cluster 2 FR: {mean_fr_cluster_2:.2f} Hz")
        print(f"Max FR: {max_fr_cluster_2:.2f} Hz")
        print(f"Min FR: {min_fr_cluster_2:.2f} Hz")
        print(f"CoV FR: {cov_fr_cluster_2:.2f}")
        
        mean_ff_cluster_2 = np.mean(cluster_2_ffs)
        max_ff_cluster_2 = np.max(cluster_2_ffs)
        min_ff_cluster_2 = np.min(cluster_2_ffs)
        cov_ff_cluster_2 = np.std(cluster_2_ffs) / mean_ff_cluster_2
        # print in multiple lines
        print(f"Cluster 2 FF: {mean_ff_cluster_2:.2f}")
        print(f"Max FF: {max_ff_cluster_2:.2f}")
        print(f"Min FF: {min_ff_cluster_2:.2f}")
        print(f"CoV FF: {cov_ff_cluster_2:.2f}")
        
        # assign cluster labels
        cluster_info_dict = {}
        cluster_info_dict['cluster_metrics'] = {}
        cluster_metrics = cluster_info_dict['cluster_metrics']
        unique_cluster_labels = np.unique(cluster_labels)
        centroids = kmeans.cluster_centers_
        for cluster_label in unique_cluster_labels:
            cluster_metrics[cluster_label] = {
                'centroid': centroids[cluster_label],
                'mean_fr': np.mean(cluster_1_frs) if cluster_label == 0 else np.mean(cluster_2_frs),
                'mean_ff': np.mean(cluster_1_ffs) if cluster_label == 0 else np.mean(cluster_2_ffs),
                'max_fr': np.max(cluster_1_frs) if cluster_label == 0 else np.max(cluster_2_frs),
                'max_ff': np.max(cluster_1_ffs) if cluster_label == 0 else np.max(cluster_2_ffs),
                'min_fr': np.min(cluster_1_frs) if cluster_label == 0 else np.min(cluster_2_frs),
                'min_ff': np.min(cluster_1_ffs) if cluster_label == 0 else np.min(cluster_2_ffs),
                'cov_fr': np.std(cluster_1_frs) / np.mean(cluster_1_frs) if cluster_label == 0 else np.std(cluster_2_frs) / np.mean(cluster_2_frs),
                'cov_ff': np.std(cluster_1_ffs) / np.mean(cluster_1_ffs) if cluster_label == 0 else np.std(cluster_2_ffs) / np.mean(cluster_2_ffs)
            }
            
        mean_fr_all_clusters = [cluster_metrics[cluster_label]['mean_fr'] for cluster_label in unique_cluster_labels]
        
        # inhib cluster will be the one with faster mean fr
        inhibitory_cluster = np.argmax(mean_fr_all_clusters)
        excitatory_cluster = np.argmin(mean_fr_all_clusters)
        
        # assign cluster descriptions
        for cluster_id, cluster_metric in cluster_metrics.items():
            cluster_metric['description'] = 'inhib' if cluster_id == inhibitory_cluster else 'excit'
            
        # save cluster percentages
        n_units = len(cluster_labels)
        n_inhib = np.sum(cluster_labels == inhibitory_cluster)
        n_excit = n_units - n_inhib
        print(f"Cluster 0 (Inhibitory) contains {n_inhib} units ({n_inhib/n_units*100:.2f}%)")
        print(f"Cluster 1 (Excitatory) contains {n_excit} units ({n_excit/n_units*100:.2f}%)")
        cluster_info_dict['n_units'] = n_units
        cluster_info_dict['n_inhib'] = n_inhib
        cluster_info_dict['n_excit'] = n_excit
        cluster_info_dict['perc_inhib'] = n_inhib/n_units*100
        cluster_info_dict['perc_excit'] = n_excit/n_units*100
        
        #print
        print(f'Inhibitory cluster - contains {n_inhib} units ({n_inhib/n_units*100:.2f}%)')
        print(f'Excitatory cluster - contains {n_excit} units ({n_excit/n_units*100:.2f}%)')
        
        # inhibitory_cluster = np.argmax([mean_fr_cluster_1, mean_fr_cluster_2])
        # #inhibitory_cluster = np.argmax([mean_ff_cluster_1, mean_ff_cluster_2])
        # excitatory_cluster = np.argmin([mean_fr_cluster_1, mean_fr_cluster_2])
        # cluster_descriptions = np.where(cluster_labels == inhibitory_cluster, 'inhib', 'excit')
        
        # #print percent of units in each cluster
        # n_units = len(cluster_labels)
        # n_inhib = np.sum(cluster_labels == inhibitory_cluster)
        # n_excit = n_units - n_inhib
        # print(f"Cluster 0 (Inhibitory) contains {n_inhib} units ({n_inhib/n_units*100:.2f}%)")
        # print(f"Cluster 1 (Excitatory) contains {n_excit} units ({n_excit/n_units*100:.2f}%)")
    
    # Plot the clusters
    plot_clusters = True
    if plot_clusters:
        plt.figure(figsize=(8, 6))
        #plt.scatter(features_pca[:, 0], features_pca[:, 1], c=cluster_labels, s=50, alpha=0.5, cmap='bwr')
        
        #scatter inhib and excit neurons separately
        # excitatory_units = features_pca[cluster_labels == 1]
        # inhibitory_units = features_pca[cluster_labels == 0]
        excitatory_units = features_pca[cluster_labels == excitatory_cluster]
        inhibitory_units = features_pca[cluster_labels == inhibitory_cluster]
        plt.scatter(excitatory_units[:, 0], excitatory_units[:, 1], c='blue', s=100, alpha=0.5, label='Excitatory')
        plt.scatter(inhibitory_units[:, 0], inhibitory_units[:, 1], c='red', s=100, alpha=0.5, label='Inhibitory')
        
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        plt.title('K-Means Clustering of Neuronal Features')
        
        # add centroids - NOTE: these aren't classified as inhib/excit at this step
        centroids = kmeans.cluster_centers_
        
        # Plot inhibitory and excitatory centroids with different colors and larger empty circles
        plt.scatter(centroids[inhibitory_cluster, 0], centroids[inhibitory_cluster, 1], edgecolor='red', facecolor='none', s=300, label='Inhib. Centroid')
        plt.scatter(centroids[excitatory_cluster, 0], centroids[excitatory_cluster, 1], edgecolor='blue', facecolor='none', s=300, label='Excit. Centroid')
        
        # save
        plt.tight_layout()
        plt.legend()
        #plt.savefig(os.path.join(wfs_output_path, 'cluster_plot.png'))
        #replace 'waveforms' w/ 'class_cluster_plots'
        #cluster_output_path = os.path.join(wfs_output_path, 'class_cluster_plots')
        cluster_output_path = wfs_output_path.replace('waveforms', 'class_cluster_plots')
        if not os.path.exists(cluster_output_path):
            os.makedirs(cluster_output_path)
        #plt.savefig(os.path.join(cluster_output_path, 'cluster_plot.png'))
        png_path = os.path.join(cluster_output_path, 'cluster_plot.png')
        pdf_path = png_path.replace('png', 'pdf')
        plt.savefig(png_path, dpi=300)
        plt.savefig(pdf_path)
        print(f"Saved cluster plot to {png_path}")
        print(f"Saved cluster plot to {pdf_path}")
        plt.close()
        
        # import sys
        # sys.exit()

        #print(f"Saved cluster plot to {os.path.join(cluster_output_path, 'cluster_plot.png')}")
    
    # ** Store Classification Results **
    # init dict
    #classified_units = {unit_id: {'desc': cluster_descriptions[i]} for i, unit_id in enumerate(feature_dicts.keys())}
    classified_units = {unit_id: {'desc': cluster_metrics[cluster_labels[unit_idx]]['description']} for unit_idx, unit_id in enumerate(feature_dicts.keys())}
    
    # add data to dict
    for unit in classified_units.keys():
        classified_units[unit]['features'] = feature_dicts[unit]
                
        # include pcs as a list
        classified_units[unit]['pca_features'] = features_pca[i].tolist() #NOTE: use to plot later if needed
        
        
        #classified_units[unit]['cluster_label'] = cluster_descriptions[i]
        
        
    return classified_units, include_units, exclude_units, features_pca, cluster_info_dict
            
def classify_neurons_v2(network_data, plot_wfs=False, **kwargs):
    
    classified_units, include_units, exclude_units, features_pca, cluster_info_dict = classify(network_data, plot_wfs=plot_wfs, **kwargs)
    
    # ** Store Classification Results **
    classification_dict = {
        'classified_units': classified_units,
        'include_units': include_units,
        'exclude_units': exclude_units,
        'features_pca': features_pca,
        'cluster_info': cluster_info_dict   
    }
    
    # save as json
    # import json
    # classification_output_path = os.path.join(network_data['waveform_output'], 'classified_units.json')
    # with open(classification_output_path, 'w') as f:
    #     json.dump(classified_units, f, indent=4)
    # print(f"Saved classified_units to {classification_output_path}")
    
    #print(classified_units)
    
    network_data['classification_output'] = classification_dict
    
    return network_data   
    
    # ## Extract Unit Locations
    # channel_locations = we.get_channel_locations()
    # unit_locations = spost.compute_unit_locations(we)
    # unit_locations_dict = {unit_id: unit_locations[i] for i, unit_id in enumerate(unit_ids)}
    # unit_locations = unit_locations_dict

def classify_neurons(network_data, plot_wfs=False, **kwargs):
    """
    Classifies neurons into clusters based on multiple waveform and burst activity features.
    Uses PCA for dimensionality reduction, K-means clustering, and Mahalanobis distance-based
    outlier removal.

    Parameters:
    - network_data (dict): Dictionary containing waveform, spiking, and bursting data.
    - sorting_object: SpikeInterface sorting object with unit IDs.

    Returns:
    - classification_dict (dict): Dictionary mapping unit_id -> 'Cluster 1' or 'Cluster 2'.
    """
    classification_data = {}
    classification_dict = {}

    # Load waveforms
    wfs_output_path = network_data['waveform_output']
    sorting_object = kwargs['sorting_object']
    #recording_object = kwargs['recording_object']
    sorting_object.remove_empty_units()
    unit_ids = sorting_object.get_unit_ids()
    we = mea.load_waveforms(wfs_output_path)  # Waveform extractor
    channel_locations = we.get_channel_locations()
    unit_locations = spost.compute_unit_locations(we)
    unit_locations_dict = {unit_id: unit_locations[i] for i, unit_id in enumerate(unit_ids)}
    unit_locations = unit_locations_dict
    sampling_frequency = we.sampling_frequency

    # Feature storage
    features = []
    unit_id_list = []
    
    # data to keep
    waveform_durations = {} 
    should_exclude = []    #units excluded based on bad waveforms

    def extract_features(we, unit_ids, network_data, sorting_object):
        """
        Extracts neuronal features related to bursting and firing patterns.
        """
        for unit_id in unit_ids:
            
            try:
                # Retrieve waveforms for a specific unit
                unit_wfs = we.get_waveforms(unit_id)  # Shape: (n_spikes, n_samples, n_channels)

                # Compute mean waveform across all spikes
                avg_waveform = np.nanmean(unit_wfs, axis=0)  # Shape: (n_samples, n_channels)

                # Find the best channel (largest absolute amplitude)
                best_channel_idx = np.argmax(np.max(np.abs(avg_waveform), axis=0))  # Index of best channel

                # Extract waveform for the best channel
                best_avg_waveform = avg_waveform[:, best_channel_idx]  # Shape: (n_samples,)

                # Plot all waveforms for the best channel in grey, with mean in red
                #plot = True
                #plot = False
                plot = plot_wfs
                if plot:
                    plt.figure(figsize=(10, 4))
                    for i in range(unit_wfs.shape[0]):
                        try:
                            wf_i = unit_wfs[i, :, best_channel_idx]
                            # Determine the half-amplitude points
                            # First half-amplitude point (rising phase)
                            # Find negative peak (trough) and subsequent positive peak
                            trough_idx = np.argmin(wf_i)
                            peak_idx = np.argmax(wf_i[trough_idx:]) + trough_idx
                            first_half_idx = np.where(wf_i[:trough_idx] <= wf_i[trough_idx] + half_amplitude)[0][-1]                        
                            plt.plot(unit_wfs[i, :, best_channel_idx], color='grey', alpha=0.1)  # Individual waveforms
                        except:
                            continue
                    plt.plot(best_avg_waveform, color='red', linewidth=2, label="Mean Waveform")  # Mean waveform

                    plt.title(f'Unit {unit_id} Waveforms (Best Channel: {best_channel_idx})')
                    plt.xlabel('Time (samples)')
                    plt.ylabel('Amplitude (uV)')
                    plt.legend()
                    plt.tight_layout()

                    # Save figure
                    wfplot_output_path = os.path.join(wfs_output_path, 'waveform_plots', f'unit_{unit_id}_waveforms.png')
                    os.makedirs(os.path.dirname(wfplot_output_path), exist_ok=True)
                    plt.savefig(wfplot_output_path)
                    print(f"Saved waveforms plot for unit {unit_id} to {wfplot_output_path}")
                
                wf_durations = []
                for i in range(unit_wfs.shape[0]):
                    try:
                        wf_i = unit_wfs[i, :, best_channel_idx]
                    
                        # Find negative peak (trough) and subsequent positive peak
                        trough_idx = np.argmin(wf_i)
                        peak_idx = np.argmax(wf_i[trough_idx:]) + trough_idx

                        # Calculate half of the peak-to-trough amplitude
                        half_amplitude = (wf_i[peak_idx] - wf_i[trough_idx]) / 2
                        
                        # Determine the half-amplitude points
                        # First half-amplitude point (rising phase)
                        first_half_idx = np.where(wf_i[:trough_idx] <= wf_i[trough_idx] + half_amplitude)[0][-1]

                        # Second half-amplitude point (falling phase)
                        second_half_idx = np.where(wf_i[trough_idx:] >= wf_i[trough_idx] + half_amplitude)[0][0] + trough_idx

                        # Compute action potential width at half-maximum (FWHM) in milliseconds
                        ap_half_width = (second_half_idx - first_half_idx) / sampling_frequency * 1000  # in ms
                        
                        wf_durations.append(ap_half_width)
                    except:
                        continue

                # Store the result
                unit_waveform_durations = wf_durations
                mean_waveform_duration = np.nanmean(unit_waveform_durations)
                waveform_durations[unit_id] = {
                    'all_durations' : unit_waveform_durations,
                    'mean': np.nanmean(unit_waveform_durations),
                    'var': np.nanvar(unit_waveform_durations),
                }
            except:
                should_exclude.append(unit_id)
                pass        
            
            # Firing rate
            firing_rate = network_data['spiking_data']['spiking_data_by_unit'].get(unit_id, {}).get('FireRate', np.nan)

            # Mean ISI values
            mean_isi_outside_burst = network_data['bursting_data']['bursting_data_by_unit'].get(unit_id, {}).get('mean_isi_outside', np.nan)
            mean_isi_within_burst = network_data['bursting_data']['bursting_data_by_unit'].get(unit_id, {}).get('mean_isi_within', np.nan)

            approximate_firing_rate_outside_burst = 1 / mean_isi_outside_burst if mean_isi_outside_burst > 0 else np.nan
            approximate_firing_rate_within_burst = 1 / mean_isi_within_burst if mean_isi_within_burst > 0 else np.nan

            # Burst-related features
            burst_duration = network_data['bursting_data']['bursting_data_by_unit'].get(unit_id, {}).get('burst_duration', np.nan)
            spike_number = network_data['bursting_data']['bursting_data_by_unit'].get(unit_id, {}).get('spike_number', np.nan)
            intra_burst_spike_rate = network_data['bursting_data']['bursting_data_by_unit'].get(unit_id, {}).get('intra_burst_spike_rate', np.nan)
            intra_burst_interval = network_data['bursting_data']['bursting_data_by_unit'].get(unit_id, {}).get('intra_burst_interval', np.nan)

            # Fano factor (time window = 6s)
            fano_factor_all = network_data['bursting_data']['bursting_data_by_unit'].get(unit_id, {}).get('fano_factor_all', np.nan)
            fano_factor_within = network_data['bursting_data']['bursting_data_by_unit'].get(unit_id, {}).get('fano_factor_within', np.nan)
            fano_factor_outside = network_data['bursting_data']['bursting_data_by_unit'].get(unit_id, {}).get('fano_factor_outside', np.nan)
            fano_factor = network_data['spiking_data']['spiking_data_by_unit'].get(unit_id, {}).get('fano_factor', np.nan)

            # ISI histogram-based CV2 measure
            cv2 = network_data['bursting_data']['bursting_data_by_unit'].get(unit_id, {}).get('cv2', np.nan)

            #HACK
            if fano_factor is None:
                fano_factor = 0.0
            
            # Store extracted features
            feature_vector = [
                mean_waveform_duration,
                fano_factor,
                firing_rate,
                approximate_firing_rate_outside_burst,
                #approximate_firing_rate_within_burst,
                #burst_duration,
                #spike_number,
                #intra_burst_spike_rate,
                #intra_burst_interval,
                #fano_factor_all,
                #fano_factor_within,
                fano_factor_outside,
                #cv2,
            ]

            features.append(feature_vector)
            unit_id_list.append(unit_id)

    extract_features(we, unit_ids, network_data, sorting_object)

    # Convert to numpy array
    features = np.array(features)
    
    # change Nones to nans
    #features = np.where(features == None, np.nan, features)

    # keep nans?
    # # Create NaN indicators to preserve missing information
    # nan_mask = np.isnan(features).astype(float)  # 1 for NaN, 0 for non-NaN
    # features_with_nan_info = np.hstack((features, nan_mask))  # Append mask as extra features

    # # Impute missing values using IterativeImputer
    # imputer = IterativeImputer(
    #     #random_state=42, 
    #     max_iter=10, initial_strategy="mean")
    # features_imputed = imputer.fit_transform(features_with_nan_info)

    # # Normalize features
    # scaler = StandardScaler()
    # features_scaled = scaler.fit_transform(features_imputed)
    
    #get rid of nans?
    # replace all nans with 0
    features = np.nan_to_num(features, nan=0.0)
    scaler = StandardScaler()
    features_scaled = scaler.fit_transform(features)

    # Reduce dimensionality using PCA (keep 2 components)
    #pca = PCA(n_components=2)
    pca = PCA(n_components=3)  # Keep 95% of variance
    features_pca = pca.fit_transform(features_scaled)

    # Apply K-means clustering
    kmeans = KMeans(n_clusters=2, 
                    #random_state=42
                    )
    cluster_labels = kmeans.fit_predict(features_pca)

    # Outlier Removal Using Mahalanobis Distance
    remove_outliers = False
    if remove_outliers:
        cluster_centroid = np.mean(features_pca, axis=0)
        cov_matrix = np.cov(features_pca, rowvar=False)
        inv_cov_matrix = np.linalg.inv(cov_matrix)

        distances = np.array([mahalanobis(point, cluster_centroid, inv_cov_matrix) for point in features_pca])
        outlier_threshold = 1.4  # Fixed threshold as per original method

        valid_units = distances < outlier_threshold
        features_pca_filtered = features_pca[valid_units]
        cluster_labels_filtered = cluster_labels[valid_units]
        unit_id_list_filtered = [unit_id_list[i] for i in range(len(unit_id_list)) if valid_units[i]]
        
        #redefine
        unit_id_list = unit_id_list_filtered
        cluster_labels = cluster_labels_filtered
        features_pca = features_pca_filtered

    # # Plot PCA-clustered data
    # plt.figure(figsize=(10, 8))
    # plt.scatter(features_pca_filtered[:, 0], features_pca_filtered[:, 1], c=cluster_labels_filtered, cmap='viridis', marker='o')
    # plt.title('PCA of Neuron Features with K-means Clustering')
    # plt.xlabel('Principal Component 1')
    # plt.ylabel('Principal Component 2')
    # plt.colorbar(label='Cluster Label')
    # plt.grid()
    # well_path = os.path.dirname(wfs_output_path)
    # plt.savefig(os.path.join(os.path.dirname(well_path), 'PCA_KMeans_Clustering.png'))
    # print(f"Saved PCA and K-means clustering plot to {os.path.join(os.path.dirname(well_path), 'PCA_KMeans_Clustering.png')}")
    
    # Assign classifications
    for i, unit_id in enumerate(unit_id_list):
        classification_dict[unit_id] = f"Cluster {cluster_labels[i] + 1}"
        
    # Average firing rates for each cluster
    cluster_firing_rates = {}
    for cluster in np.unique(cluster_labels):
        cluster_units = [unit_id_list[i] for i in range(len(unit_id_list)) if cluster_labels[i] == cluster]
        cluster_firing_rates[cluster] = np.mean([network_data['spiking_data']['spiking_data_by_unit'][unit_id]['FireRate'] for unit_id in cluster_units])
    print(f"Cluster firing rates: {cluster_firing_rates}")
    
    # assume faster firing rate is inhibitory cluster
    inhibitory_cluster = 1 if cluster_firing_rates[0] < cluster_firing_rates[1] else 0

    # Compute and print neuron ratios
    # cluster_1_count = sum(1 for unit_id in classification_dict if classification_dict[unit_id] == "Cluster 1")
    # cluster_2_count = sum(1 for unit_id in classification_dict if classification_dict[unit_id] == "Cluster 2")
    # inhibitory_cluster_count = sum(1 for unit_id in classification_dict if classification_dict[unit_id] == f"Cluster {inhibitory_cluster + 1}")
    # excitatory_cluster_count = sum(1 for unit_id in classification_dict if classification_dict[unit_id] == f"Cluster {1 - inhibitory_cluster + 1}")
    # print(f"Putative Inhibitory Neurons: {inhibitory_cluster_count}")
    # print(f"Putative Excitatory Neurons: {excitatory_cluster_count}")
    # print(f'Percent Inhibitory: {inhibitory_cluster_count / (inhibitory_cluster_count + excitatory_cluster_count) * 100:.2f}%')
    # print(f'Percent Excitatory: {excitatory_cluster_count / (inhibitory_cluster_count + excitatory_cluster_count) * 100:.2f}%')
    
    plot=True
    if plot:
    
        def plot_pca_clusters_2d(features_pca, cluster_labels, kmeans, inhibitory_cluster, well_path):
            plt.figure(figsize=(10, 8))
            plt.scatter(features_pca[:, 0], features_pca[:, 1], c=cluster_labels, cmap='viridis', marker='o')
            plt.scatter(kmeans.cluster_centers_[:, 0], kmeans.cluster_centers_[:, 1], c='red', marker='o', s=200, label='Centroids')
            
            plt.title('PCA of Neuron Features with K-means Clustering')
            plt.xlabel('Principal Component 1')
            plt.ylabel('Principal Component 2')
            
            for i, center in enumerate(kmeans.cluster_centers_):
                cluster_id = i + 1
                if i == inhibitory_cluster:
                    plt.text(center[0], center[1]+0.1, 'Inhibitory Cluster', fontsize=12, ha='center', va='center', color='red')
                else:
                    plt.text(center[0], center[1]+0.1, 'Excitatory Cluster', fontsize=12, ha='center', va='center', color='red')
            
            plt.legend()
            plt.grid()
            plt.tight_layout()
            plt.savefig(os.path.join(well_path, 'PCA_KMeans_Clustering_with_Centroids_2D.png'))
            print(f"Saved PCA and K-means clustering plot to {os.path.join(well_path, 'PCA_KMeans_Clustering_with_Centroids_2D.png')}")
        
        def plot_pca_clusters_3d(features_pca, cluster_labels, kmeans, inhibitory_cluster, well_path):
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter3D(features_pca[:, 0], features_pca[:, 1], features_pca[:, 2], c=cluster_labels, cmap='viridis', marker='o')
            ax.scatter3D(kmeans.cluster_centers_[:, 0], kmeans.cluster_centers_[:, 1], kmeans.cluster_centers_[:, 2], c='red', marker='o', s=200, label='Centroids')
            
            ax.set_title('PCA of Neuron Features with K-means Clustering')
            ax.set_xlabel('Principal Component 1')
            ax.set_ylabel('Principal Component 2')
            ax.set_zlabel('Principal Component 3')
            
            for i, center in enumerate(kmeans.cluster_centers_):
                cluster_id = i + 1
                if i == inhibitory_cluster:
                    ax.text(center[0], center[1]+0.1, center[2], 'Inhibitory Cluster', fontsize=12, ha='center', va='center', color='red')
                else:
                    ax.text(center[0], center[1]+0.1, center[2], 'Excitatory Cluster', fontsize=12, ha='center', va='center', color='red')
            
            plt.legend()
            plt.grid()
            plt.tight_layout()
            plt.savefig(os.path.join(well_path, 'PCA_KMeans_Clustering_with_Centroids_3D.png'))
            print(f"Saved PCA and K-means clustering plot to {os.path.join(well_path, 'PCA_KMeans_Clustering_with_Centroids_3D.png')}")
        
        # Example usage
        well_path = os.path.dirname(wfs_output_path)
        #plot_pca_clusters_2d(features_pca, cluster_labels, kmeans, inhibitory_cluster, well_path)
        plot_pca_clusters_3d(features_pca, cluster_labels, kmeans, inhibitory_cluster, well_path)
        # # Print cluster statistics
        
        #remove waveform exclusions from unit_id_list
        unit_id_list = [unit_id for unit_id in unit_id_list if unit_id not in should_exclude]
        
        # plot histogram of waveform durations
        plt.figure(figsize=(10, 6))
        durations_inhib = [waveform_durations[unit_id]['mean'] for unit_id in unit_id_list if classification_dict[unit_id] == f"Cluster {inhibitory_cluster + 1}"]
        durations_excit = [waveform_durations[unit_id]['mean'] for unit_id in unit_id_list if classification_dict[unit_id] == f"Cluster {1 - inhibitory_cluster + 1}"]
        #remove nans
        durations_inhib = [x for x in durations_inhib if not np.isnan(x)]
        durations_excit = [x for x in durations_excit if not np.isnan(x)]
        
        bins = np.histogram_bin_edges(durations_inhib + durations_excit, bins=100)
        plt.hist(durations_inhib, bins=bins, alpha=0.5, label='Inhibitory Cluster')
        plt.hist(durations_excit, bins=bins, alpha=0.5, label='Excitatory Cluster')
        plt.title('Waveform Duration Distribution by Cluster')
        plt.xlabel('Waveform Duration (ms)')
        plt.ylabel('Number of Neurons')
        plt.legend()
        #plt.savefig(os.path.join(os.path.dirname(well_path), 'Waveform_Duration_Distribution_by_Cluster.png'))
        #print(f"Saved waveform duration distribution plot to {os.path.join(os.path.dirname(well_path), 'Waveform_Duration_Distribution_by_Cluster.png')}")
        plt.savefig(os.path.join(well_path, 'Waveform_Duration_Distribution_by_Cluster.png'))
        print(f"Saved waveform duration distribution plot to {os.path.join(well_path, 'Waveform_Duration_Distribution_by_Cluster.png')}")
                
        # plot histogram of firing rates, visualize separation of clusters
        plt.figure(figsize=(10, 6))
        # plt.hist([network_data['spiking_data']['spiking_data_by_unit'][unit_id]['FireRate'] for unit_id in unit_id_list if classification_dict[unit_id] == f"Cluster {inhibitory_cluster + 1}"], bins=20, alpha=0.5, label='Inhibitory Cluster')
        # plt.hist([network_data['spiking_data']['spiking_data_by_unit'][unit_id]['FireRate'] for unit_id in unit_id_list if classification_dict[unit_id] == f"Cluster {1 - inhibitory_cluster + 1}"], bins=20, alpha=0.5, label='Excitatory Cluster')
        firing_rates_inhib = [network_data['spiking_data']['spiking_data_by_unit'][unit_id]['FireRate'] for unit_id in unit_id_list if classification_dict[unit_id] == f"Cluster {inhibitory_cluster + 1}"]
        firing_rates_excit = [network_data['spiking_data']['spiking_data_by_unit'][unit_id]['FireRate'] for unit_id in unit_id_list if classification_dict[unit_id] == f"Cluster {1 - inhibitory_cluster + 1}"]
        bins = np.histogram_bin_edges(firing_rates_inhib + firing_rates_excit, bins=20)
        plt.hist(firing_rates_inhib, bins=bins, alpha=0.5, label='Inhibitory Cluster')
        plt.hist(firing_rates_excit, bins=bins, alpha=0.5, label='Excitatory Cluster')
        plt.title('Firing Rate Distribution by Cluster')
        plt.xlabel('Firing Rate (Hz)')
        plt.ylabel('Number of Neurons')
        plt.legend()
        #plt.savefig(os.path.join(os.path.dirname(well_path), 'Firing_Rate_Distribution_by_Cluster.png'))
        # print(f"Cluster 1 Neurons: {cluster_1_count}")
        #print(f"Saved firing rate distribution plot to {os.path.join(os.path.dirname(well_path), 'Firing_Rate_Distribution_by_Cluster.png')}")
        plt.savefig(os.path.join(well_path, 'Firing_Rate_Distribution_by_Cluster.png'))
        print(f"Saved firing rate distribution plot to {os.path.join(well_path, 'Firing_Rate_Distribution_by_Cluster.png')}")
        
        # plot histogram firing rates vs within_burst firig rates, visualize separation of clusters
        plt.figure(figsize=(10, 6))
        within_bursts_ISI_inhib = [network_data['bursting_data']['bursting_data_by_unit'][unit_id]['mean_isi_within'] for unit_id in unit_id_list if classification_dict[unit_id] == f"Cluster {inhibitory_cluster + 1}"]
        within_bursts_ISI_excit = [network_data['bursting_data']['bursting_data_by_unit'][unit_id]['mean_isi_within'] for unit_id in unit_id_list if classification_dict[unit_id] == f"Cluster {1 - inhibitory_cluster + 1}"]
        firing_rates_inhib = [1/x for x in within_bursts_ISI_inhib if not np.isnan(x)]
        firing_rates_excit = [1/x for x in within_bursts_ISI_excit if not np.isnan(x)]
        bins = np.histogram_bin_edges(firing_rates_inhib + firing_rates_excit, bins=20)
        plt.hist(firing_rates_inhib, bins=bins, alpha=0.5, label='Inhibitory Cluster')
        plt.hist(firing_rates_excit, bins=bins, alpha=0.5, label='Excitatory Cluster')
        plt.title('Firing Rate Distribution by Cluster')
        plt.xlabel('Firing Rate (Hz)')
        plt.ylabel('Number of Neurons')
        plt.legend()
        #plt.savefig(os.path.join(os.path.dirname(well_path), 'Firing_Rate_Distribution_by_Cluster_within_burst.png'))
        #print(f"Saved firing rate within burst distribution plot to {os.path.join(os.path.dirname(well_path), 'Firing_Rate_Distribution_by_Cluster_within_burst.png')}")
        plt.savefig(os.path.join(well_path, 'Firing_Rate_Distribution_by_Cluster_within_burst.png'))
        print(f"Saved firing rate within burst distribution plot to {os.path.join(well_path, 'Firing_Rate_Distribution_by_Cluster_within_burst.png')}")
        
        # plot histogram of firing rates outside bursts, visualize separation of clusters
        plt.figure(figsize=(10, 6))
        outside_bursts_ISI_inhib = [network_data['bursting_data']['bursting_data_by_unit'][unit_id]['mean_isi_outside'] for unit_id in unit_id_list if classification_dict[unit_id] == f"Cluster {inhibitory_cluster + 1}"]
        outside_bursts_ISI_excit = [network_data['bursting_data']['bursting_data_by_unit'][unit_id]['mean_isi_outside'] for unit_id in unit_id_list if classification_dict[unit_id] == f"Cluster {1 - inhibitory_cluster + 1}"]
        firing_rates_inhib = [1/x for x in outside_bursts_ISI_inhib if not np.isnan(x)]
        firing_rates_excit = [1/x for x in outside_bursts_ISI_excit if not np.isnan(x)]
        bins = np.histogram_bin_edges(firing_rates_inhib + firing_rates_excit, bins=20)
        plt.hist(firing_rates_inhib, bins=bins, alpha=0.5, label='Inhibitory Cluster')
        plt.hist(firing_rates_excit, bins=bins, alpha=0.5, label='Excitatory Cluster')
        plt.title('Firing Rate Distribution by Cluster')
        plt.xlabel('Firing Rate (Hz)')
        plt.ylabel('Number of Neurons')
        plt.legend()
        #plt.savefig(os.path.join(os.path.dirname(well_path), 'Firing_Rate_Distribution_by_Cluster_outside_burst.png'))
        #print(f"Saved firing rate outside burst distribution plot to {os.path.join(os.path.dirname(well_path), 'Firing_Rate_Distribution_by_Cluster_outside_burst.png')}")
        plt.savefig(os.path.join(well_path, 'Firing_Rate_Distribution_by_Cluster_outside_burst.png'))
        print(f"Saved firing rate outside burst distribution plot to {os.path.join(well_path, 'Firing_Rate_Distribution_by_Cluster_outside_burst.png')}")
            
        # plot histogram of fano factors, visualize separation of clusters
        plt.figure(figsize=(10, 6))    
        fano_factor_hist_inhib_data = [network_data['spiking_data']['spiking_data_by_unit'][unit_id]['fano_factor'] for unit_id in unit_id_list if classification_dict[unit_id] == f"Cluster {inhibitory_cluster + 1}" and network_data['spiking_data']['spiking_data_by_unit'][unit_id]['fano_factor'] is not None]
        fano_factor_hist_excit_data = [network_data['spiking_data']['spiking_data_by_unit'][unit_id]['fano_factor'] for unit_id in unit_id_list if classification_dict[unit_id] == f"Cluster {1 - inhibitory_cluster + 1}" and network_data['spiking_data']['spiking_data_by_unit'][unit_id]['fano_factor'] is not None]
        #remove nans if any
        fano_factor_hist_excit_data = [x for x in fano_factor_hist_excit_data if not np.isnan(x)]
        fano_factor_hist_inhib_data = [x for x in fano_factor_hist_inhib_data if not np.isnan(x)]
        bins = np.histogram_bin_edges(fano_factor_hist_inhib_data + fano_factor_hist_excit_data, bins=20)
        plt.hist(fano_factor_hist_inhib_data, bins=bins, alpha=0.5, label='Inhibitory Cluster')
        plt.hist(fano_factor_hist_excit_data, bins=bins, alpha=0.5, label='Excitatory Cluster')
        plt.title('Fano Factor Distribution by Cluster')
        plt.xlabel('Fano Factor')
        plt.ylabel('Number of Neurons')
        plt.legend()
        # plt.savefig(os.path.join(os.path.dirname(well_path), 'Fano_Factor_Distribution_by_Cluster.png'))
        # print(f"Saved fano factor distribution plot to {os.path.join(os.path.dirname(well_path), 'Fano_Factor_Distribution_by_Cluster.png')}")
        plt.savefig(os.path.join(well_path, 'Fano_Factor_Distribution_by_Cluster.png'))
        print(f"Saved fano factor distribution plot to {os.path.join(well_path, 'Fano_Factor_Distribution_by_Cluster.png')}")
        
        # repeat for fano factor within burst and fano factor outside burst
        plt.figure(figsize=(10, 6))
        fano_factor_hist_inhib_data = [network_data['bursting_data']['bursting_data_by_unit'][unit_id]['fano_factor_within'] for unit_id in unit_id_list if classification_dict[unit_id] == f"Cluster {inhibitory_cluster + 1}" and network_data['bursting_data']['bursting_data_by_unit'][unit_id]['fano_factor_within'] is not None]
        fano_factor_hist_excit_data = [network_data['bursting_data']['bursting_data_by_unit'][unit_id]['fano_factor_within'] for unit_id in unit_id_list if classification_dict[unit_id] == f"Cluster {1 - inhibitory_cluster + 1}" and network_data['bursting_data']['bursting_data_by_unit'][unit_id]['fano_factor_within'] is not None]
        #remove nans if any
        fano_factor_hist_excit_data = [x for x in fano_factor_hist_excit_data if not np.isnan(x)]
        fano_factor_hist_inhib_data = [x for x in fano_factor_hist_inhib_data if not np.isnan(x)]
        bins = np.histogram_bin_edges(fano_factor_hist_inhib_data + fano_factor_hist_excit_data, bins=20)
        plt.hist(fano_factor_hist_inhib_data, bins=bins, alpha=0.5, label='Inhibitory Cluster')
        plt.hist(fano_factor_hist_excit_data, bins=bins, alpha=0.5, label='Excitatory Cluster')
        plt.title('Fano Factor Distribution by Cluster')
        plt.xlabel('Fano Factor')
        plt.ylabel('Number of Neurons')
        plt.legend()
        # plt.savefig(os.path.join(os.path.dirname(well_path), 'Fano_Factor_Distribution_by_Cluster_within_burst.png'))
        # print(f"Saved fano factor within burst distribution plot to {os.path.join(os.path.dirname(well_path), 'Fano_Factor_Distribution_by_Cluster_within_burst.png')}")
        plt.savefig(os.path.join(well_path, 'Fano_Factor_Distribution_by_Cluster_within_burst.png'))
        print(f"Saved fano factor within burst distribution plot to {os.path.join(well_path, 'Fano_Factor_Distribution_by_Cluster_within_burst.png')}")
        
        plt.figure(figsize=(10, 6))
        fano_factor_hist_inhib_data = [network_data['bursting_data']['bursting_data_by_unit'][unit_id]['fano_factor_outside'] for unit_id in unit_id_list if classification_dict[unit_id] == f"Cluster {inhibitory_cluster + 1}" and network_data['bursting_data']['bursting_data_by_unit'][unit_id]['fano_factor_outside'] is not None]
        fano_factor_hist_excit_data = [network_data['bursting_data']['bursting_data_by_unit'][unit_id]['fano_factor_outside'] for unit_id in unit_id_list if classification_dict[unit_id] == f"Cluster {1 - inhibitory_cluster + 1}" and network_data['bursting_data']['bursting_data_by_unit'][unit_id]['fano_factor_outside'] is not None]
        #remove nans if any
        fano_factor_hist_excit_data = [x for x in fano_factor_hist_excit_data if not np.isnan(x)]
        fano_factor_hist_inhib_data = [x for x in fano_factor_hist_inhib_data if not np.isnan(x)]
        bins = np.histogram_bin_edges(fano_factor_hist_inhib_data + fano_factor_hist_excit_data, bins=20)
        plt.hist(fano_factor_hist_inhib_data, bins=bins, alpha=0.5, label='Inhibitory Cluster')
        plt.hist(fano_factor_hist_excit_data, bins=bins, alpha=0.5, label='Excitatory Cluster')
        plt.title('Fano Factor Distribution by Cluster')
        plt.xlabel('Fano Factor')
        plt.ylabel('Number of Neurons')
        plt.legend()
        # plt.savefig(os.path.join(os.path.dirname(well_path), 'Fano_Factor_Distribution_by_Cluster_outside_burst.png'))
        # print(f"Saved fano factor outside burst distribution plot to {os.path.join(os.path.dirname(well_path), 'Fano_Factor_Distribution_by_Cluster_outside_burst.png')}")
        plt.savefig(os.path.join(well_path, 'Fano_Factor_Distribution_by_Cluster_outside_burst.png'))
        print(f"Saved fano factor outside burst distribution plot to {os.path.join(well_path, 'Fano_Factor_Distribution_by_Cluster_outside_burst.png')}")
        
        #all
        plt.figure(figsize=(10, 6))
        fano_factor_hist_inhib_data = [network_data['bursting_data']['bursting_data_by_unit'][unit_id]['fano_factor_all'] for unit_id in unit_id_list if classification_dict[unit_id] == f"Cluster {inhibitory_cluster + 1}" and network_data['bursting_data']['bursting_data_by_unit'][unit_id]['fano_factor_all'] is not None]
        fano_factor_hist_excit_data = [network_data['bursting_data']['bursting_data_by_unit'][unit_id]['fano_factor_all'] for unit_id in unit_id_list if classification_dict[unit_id] == f"Cluster {1 - inhibitory_cluster + 1}" and network_data['bursting_data']['bursting_data_by_unit'][unit_id]['fano_factor_all'] is not None]
        #remove nans if any
        fano_factor_hist_excit_data = [x for x in fano_factor_hist_excit_data if not np.isnan(x)]
        fano_factor_hist_inhib_data = [x for x in fano_factor_hist_inhib_data if not np.isnan(x)]    
        bins = np.histogram_bin_edges(fano_factor_hist_inhib_data + fano_factor_hist_excit_data, bins=20)
        plt.hist(fano_factor_hist_inhib_data, bins=bins, alpha=0.5, label='Inhibitory Cluster')
        plt.hist(fano_factor_hist_excit_data, bins=bins, alpha=0.5, label='Excitatory Cluster')
        plt.title('Fano Factor Distribution by Cluster')
        plt.xlabel('Fano Factor')
        plt.ylabel('Number of Neurons')
        plt.legend()
        # plt.savefig(os.path.join(os.path.dirname(well_path), 'Fano_Factor_Distribution_by_Cluster_all.png'))
        # print(f"Saved fano factor all distribution plot to {os.path.join(os.path.dirname(well_path), 'Fano_Factor_Distribution_by_Cluster_all.png')}")  
        plt.savefig(os.path.join(well_path, 'Fano_Factor_Distribution_by_Cluster_all.png'))
        print(f"Saved fano factor all distribution plot to {os.path.join(well_path, 'Fano_Factor_Distribution_by_Cluster_all.png')}")
        
        # plot locations of neurons on probe 2000 x 4000 um rectangle. ignore z coordinate in unit_locations
        #unit_locations_2d = np.array([[loc[0], loc[1]] for loc in unit_locations])
        unit_locations_2d = {unit_id: np.array([loc[0], loc[1]]) for unit_id,loc in unit_locations.items()}
        #unit_locations_2d = unit_locations_2d.T # transpose to get x and y coordinates
        # label inhibitory and excitatory neurons
        # inhib_neuron_locs = unit_locations_2d[[i for i, unit_id in enumerate(unit_id_list) if classification_dict[unit_id] == f"Cluster {inhibitory_cluster + 1}"]]
        # excit_neuron_locs = unit_locations_2d[[i for i, unit_id in enumerate(unit_id_list) if classification_dict[unit_id] == f"Cluster {1 - inhibitory_cluster + 1}"]]
        inhib_neuron_locs = np.array([unit_locations_2d[unit_id] for unit_id in unit_id_list if classification_dict[unit_id] == f"Cluster {inhibitory_cluster + 1}"])
        excit_neuron_locs = np.array([unit_locations_2d[unit_id] for unit_id in unit_id_list if classification_dict[unit_id] == f"Cluster {1 - inhibitory_cluster + 1}"])
        
        min_x = min(np.min(inhib_neuron_locs[:, 0]), np.min(excit_neuron_locs[:, 0]))
        min_y = min(np.min(inhib_neuron_locs[:, 1]), np.min(excit_neuron_locs[:, 1]))
        max_x = max(np.max(inhib_neuron_locs[:, 0]), np.max(excit_neuron_locs[:, 0]))
        max_y = max(np.max(inhib_neuron_locs[:, 1]), np.max(excit_neuron_locs[:, 1]))
        
        if min_x > 0: min_x = 0
        if min_y > 0: min_y = 0
        if max_x < 4000: max_x = 4000
        if max_y < 2000: max_y = 2000
        
        
        plt.figure(figsize=(10, 6))
        plt.scatter(excit_neuron_locs[:, 0], excit_neuron_locs[:, 1], c='orange', label='Excitatory Neurons', s=10)
        #plt.scatter(inhib_neuron_locs[:, 0], inhib_neuron_locs[:, 1], c='blue', label='Inhibitory Neurons', s=10)
        #plot inhib with open circles in case of overlap
        plt.scatter(inhib_neuron_locs[:, 0], inhib_neuron_locs[:, 1], facecolors='none', edgecolors='blue', label='Inhibitory Neurons', s=10)
        
        
        # make the dots smaller
        plt.title('Neuron Locations on MEA')
        plt.xlabel('X Coordinate (um)')
        plt.ylabel('Y Coordinate (um)')
        # plt.xlim(0, 2000)
        # plt.xlim(0, 4000)
        plt.xlim(min_x, max_x)
        plt.ylim(min_y, max_y)
        plt.legend()
        plt.savefig(os.path.join(os.path.dirname(wfs_output_path), 'Neuron_Locations_on_MEA.png'))
        print(f"Saved neuron locations plot to {os.path.join(os.path.dirname(wfs_output_path), 'Neuron_Locations_on_MEA.png')}")
        
    
    classification_data = {
        
        'waveform_durations': waveform_durations,
        
        'no. of units': len(unit_id_list),
        'no. of units classified:': len([unit_id for unit_id, classification in classification_dict.items() if classification in [f"Cluster {inhibitory_cluster + 1}", f"Cluster {1 - inhibitory_cluster + 1}"]]),
        'inhibitory_neurons': [unit_id for unit_id, classification in classification_dict.items() if classification == f"Cluster {inhibitory_cluster + 1}"],
        'excitatory_neurons': [unit_id for unit_id, classification in classification_dict.items() if classification == f"Cluster {1 - inhibitory_cluster + 1}"],
        'unclassified_neurons': [unit_id for unit_id, classification in classification_dict.items() if classification not in [f"Cluster {inhibitory_cluster + 1}", f"Cluster {1 - inhibitory_cluster + 1}"]],
        'unit_locations': unit_locations,
        # 'percent_inhibitory': inhibitory_cluster_count / (inhibitory_cluster_count + excitatory_cluster_count) * 100,
        # 'percent_excitatory': excitatory_cluster_count / (inhibitory_cluster_count + excitatory_cluster_count) * 100,
        'percent_inhibitory': len([unit_id for unit_id, classification in classification_dict.items() if classification == f"Cluster {inhibitory_cluster + 1}"]) / len(unit_id_list) * 100,
        'percent_excitatory': len([unit_id for unit_id, classification in classification_dict.items() if classification == f"Cluster {1 - inhibitory_cluster + 1}"]) / len(unit_id_list) * 100,
        'percent_unclassified': len([unit_id for unit_id, classification in classification_dict.items() if classification not in [f"Cluster {inhibitory_cluster + 1}", f"Cluster {1 - inhibitory_cluster + 1}"]]) / len(unit_id_list) * 100,
    }
    
    # #save classification data as json to os.path.dirname(wfs_output_path)
    # import numpy as np
    # import json
    # # Convert numpy data types to native Python types
    # def convert_to_native_types(obj):
    #     if isinstance(obj, np.ndarray):
    #         return obj.tolist()
    #     elif isinstance(obj, np.generic):
    #         return obj.item()
    #     elif isinstance(obj, dict):
    #         return {k: convert_to_native_types(v) for k, v in obj.items()}
    #     elif isinstance(obj, list):
    #         return [convert_to_native_types(i) for i in obj]
    #     else:
    #         return obj
    
    network_data['classification_data'] = classification_data
    
    return network_data
    
    
    
    #return classification_data
