import spikeinterface.full as si
import os
from spikeinterface.extractors import MaxwellRecordingExtractor
from MEA_Analysis.MEAProcessingLibrary import mea_processing_library as mea

def classify_unit_FR(
    unit_data, 
    max_E_assumption, 
    max_I_assumption, 
    min_E_assumption, 
    min_I_assumption, 
    max_E_assumption_inBurst, 
    max_I_assumption_inBurst, 
    min_E_assumption_inBurst, 
    min_I_assumption_inBurst):
    
    # Extract firing rates and ISIs
    mean_FR = unit_data['spikes']['FireRate']
    mean_isi_outside_burst = unit_data['bursts']['mean_isi_outside']
    fr_outside_burst = 1 / mean_isi_outside_burst if mean_isi_outside_burst > 0 else 0
    mean_isi_within_burst = unit_data['bursts']['mean_isi_within']
    fr_within_burst = 1 / mean_isi_within_burst if mean_isi_within_burst > 0 else 0

    # Hierarchical logic
    # Step 1: High-confidence conditions
    if fr_within_burst > 50:  # Strong indication of an inhibitory neuron
        return {
            'type_guess': 'I',
            'confidence_E': 0.0,
            'confidence_I': 1.0,
            'reasoning': 'Firing rate within burst exceeds 50 Hz.'
        }
    if mean_FR < 1 and fr_outside_burst < 1:  # Likely a quiescent excitatory neuron
        return {
            'type_guess': 'E',
            'confidence_E': 0.9,
            'confidence_I': 0.1,
            'reasoning': 'Mean and outside burst firing rates are very low.'
        }
    
    # Step 2: Voting-based classification
    E_conditions = [
        mean_FR <= max_E_assumption,
        mean_FR >= min_E_assumption,
        fr_outside_burst <= max_E_assumption,
        fr_outside_burst >= min_E_assumption,
        fr_within_burst <= max_E_assumption_inBurst,
        fr_within_burst >= min_E_assumption_inBurst,
    ]
    
    I_conditions = [
        mean_FR <= max_I_assumption,
        mean_FR >= min_I_assumption,
        fr_outside_burst <= max_I_assumption,
        fr_outside_burst >= min_I_assumption,
        fr_within_burst <= max_I_assumption_inBurst,
        fr_within_burst >= min_I_assumption_inBurst,
    ]
    
    votes_E = sum(E_conditions)
    votes_I = sum(I_conditions)

    # Step 3: Final classification based on votes
    if votes_E < votes_I:
        neuron_classification = {
            'type_guess': 'I',
            'confidence_E': votes_E / len(E_conditions),
            'confidence_I': votes_I / len(I_conditions),
            'reasoning': 'More conditions satisfied for inhibitory neuron classification.'
        }
    elif votes_E > votes_I:
        neuron_classification = {
            'type_guess': 'E',
            'confidence_E': votes_E / len(E_conditions),
            'confidence_I': votes_I / len(I_conditions),
            'reasoning': 'More conditions satisfied for excitatory neuron classification.'
        }
    else:
        neuron_classification = {
            'type_guess': 'U',
            'confidence_E': votes_E / len(E_conditions),
            'confidence_I': votes_I / len(I_conditions),
            'reasoning': 'Equally likely to be excitatory or inhibitory based on current thresholds.'
        }
    
    return neuron_classification

# #def classify_neurons(sorting_object, recording_object, sampling_rate=10000, output_path=None):
# import numpy as np
# # import spikeinterface.extractors as mea

# def classify_neurons_simple(network_data, **kwargs):
#     """
#     Classifies neurons into excitatory or inhibitory based on waveform features.

#     Parameters:
#     - network_data (dict): Dictionary containing 'waveform_output' path.
#     - sorting_object: SpikeInterface sorting object with unit IDs.

#     Returns:
#     - classification_dict (dict): Dictionary mapping unit_id -> 'excitatory' or 'inhibitory'.
#     """
#     classification_dict = {}
    
#     # Load waveforms
#     wfs_output_path = network_data['waveform_output']
#     sorting_object = kwargs['sorting_object']
#     unit_ids = sorting_object.get_unit_ids()
#     we = mea.load_waveforms(wfs_output_path) #waveform extractor
    
#     def cluster_by_peak_to_trough(we, unit_ids):
#         #Clustering logic:
#         for unit_id in unit_ids:
#             # unit_wfs is a np.memmap of np.memmap
#             unit_wfs = we.get_waveforms(unit_id)  # Shape: (n_spikes, n_samples, n_channels)
            
#             # Compute mean waveform across all spikes and channels
#             avg_waveform = np.mean(unit_wfs, axis=(0, 2))  # Shape: (n_samples,)
            
#             # Find the negative peak (trough) and the subsequent positive peak
#             trough_idx = np.argmin(avg_waveform)  # Index of minimum (negative peak)
#             peak_idx = np.argmax(avg_waveform[trough_idx:]) + trough_idx  # Index of max after trough

#             # Compute peak-to-trough duration (in samples)
#             peak_to_trough_samples = peak_idx - trough_idx
#             sampling_rate = we.sampling_frequency  # Get sampling rate in Hz
#             peak_to_trough_duration = peak_to_trough_samples / sampling_rate * 1000  # Convert to ms

#             # Classify neuron based on peak-to-trough duration
#             if peak_to_trough_duration > 0.4:
#                 classification_dict[unit_id] = "excitatory"
#             else:
#                 classification_dict[unit_id] = "inhibitory"
#         #return {unit_id: [unit_id] for unit_id in unit_ids}
#         return classification_dict
#     classification_dict = cluster_by_peak_to_trough(we, unit_ids)
            
#     # print ration excitatory to inhibitory neurons
#     excitatory_count = sum(1 for unit_id in classification_dict if classification_dict[unit_id] == "excitatory")
#     inhibitory_count = sum(1 for unit_id in classification_dict if classification_dict[unit_id] == "inhibitory")
#     print(f"Excitatory neurons: {excitatory_count}")
#     print(f"Inhibitory neurons: {inhibitory_count}")
#     ration = excitatory_count / inhibitory_count if inhibitory_count > 0 else float('inf')
#     print(f"Ratio excitatory to inhibitory neurons: {ration}")

#     return classification_dict

# import numpy as np
# import os
# import matplotlib.pyplot as plt
# from sklearn.mixture import GaussianMixture
# from sklearn.decomposition import PCA
# #from sklearn.impute import IterativeImputer
# from sklearn.preprocessing import StandardScaler

# from sklearn.experimental import enable_iterative_imputer # NOTE: from sklearn.experimental import enable_iterative_imputer
#                                                             #from sklearn.impute import IterativeImputer

# from sklearn.impute import IterativeImputer


# def classify_neurons_more_complex(network_data, **kwargs):
#     """
#     Classifies neurons into excitatory or inhibitory based on multiple waveform features 
#     using PCA for dimensionality reduction and Gaussian Mixture Models (GMM) for clustering.

#     Handles NaNs by keeping them informative with IterativeImputer and NaN indicator variables.

#     Parameters:
#     - network_data (dict): Dictionary containing 'waveform_output' path.
#     - sorting_object: SpikeInterface sorting object with unit IDs.

#     Returns:
#     - classification_dict (dict): Dictionary mapping unit_id -> 'excitatory' or 'inhibitory'.
#     """
#     classification_dict = {}

#     # Load waveforms
#     wfs_output_path = network_data['waveform_output']
#     sorting_object = kwargs['sorting_object']
#     unit_ids = sorting_object.get_unit_ids()
#     we = mea.load_waveforms(wfs_output_path)  # Waveform extractor
    
#     # Feature storage
#     features = []
#     unit_id_list = []

#     def extract_features(we, unit_ids, network_data):
#         """
#         Extracts key features: Peak-to-Trough ratio, Firing Rate, AP Width.
#         """
#         for unit_id in unit_ids:
#             # Estimate firing rate
#             firing_rate = network_data['spiking_data']['spiking_data_by_unit'].get(unit_id, {}).get('FireRate', np.nan)
            
#             # Mean ISI outside burst
#             mean_isi_outside_burst = network_data['bursting_data']['bursting_data_by_unit'].get(unit_id, {}).get('mean_isi_outside', np.nan)
#             approximate_firing_rate_outside_burst = 1 / mean_isi_outside_burst if mean_isi_outside_burst > 0 else np.nan
            
#             # Mean ISI within burst
#             mean_isi_within_burst = network_data['bursting_data']['bursting_data_by_unit'].get(unit_id, {}).get('mean_isi_within', np.nan)
#             approximate_firing_rate_within_burst = 1 / mean_isi_within_burst if mean_isi_within_burst > 0 else np.nan

#             # Fano factors
#             fano_factor_all = network_data['bursting_data']['bursting_data_by_unit'].get(unit_id, {}).get('fano_factor_all', np.nan)
#             fano_factor_within = network_data['bursting_data']['bursting_data_by_unit'].get(unit_id, {}).get('fano_factor_within', np.nan)
#             fano_factor_outside = network_data['bursting_data']['bursting_data_by_unit'].get(unit_id, {}).get('fano_factor_outside', np.nan)
            
#             # Store extracted features
#             feature_vector = [
#                 firing_rate,
#                 approximate_firing_rate_outside_burst, 
#                 approximate_firing_rate_within_burst, 
#                 fano_factor_all,
#                 fano_factor_within,
#                 fano_factor_outside,
#             ]

#             features.append(feature_vector)
#             unit_id_list.append(unit_id)

#     extract_features(we, unit_ids, network_data)
    
#     # Convert to numpy array
#     features = np.array(features)

#     # Create a mask for NaNs (to preserve their presence as a feature)
#     nan_mask = np.isnan(features).astype(float)  # 1 for NaN, 0 for non-NaN
#     features_with_nan_info = np.hstack((features, nan_mask))  # Append mask as additional features

#     # Impute missing values while preserving their importance
#     imputer = IterativeImputer(random_state=42, max_iter=10, initial_strategy="mean")
#     features_imputed = imputer.fit_transform(features_with_nan_info)

#     # Normalize features
#     scaler = StandardScaler()
#     features_scaled = scaler.fit_transform(features_imputed)

#     # Reduce dimensionality using PCA (keep 2 components)
#     pca = PCA(n_components=2)
#     features_pca = pca.fit_transform(features_scaled)

#     # Cluster using Gaussian Mixture Model (GMM)
#     gmm = GaussianMixture(n_components=2, covariance_type="full", random_state=42)
#     cluster_labels = gmm.fit_predict(features_pca)

#     # Determine which cluster is excitatory and which is inhibitory
#     cluster_mean_firing_rates = [np.nanmean(features[:, 0][cluster_labels == i]) for i in range(2)]
#     excitatory_cluster = np.argmin(cluster_mean_firing_rates)  # Lower FR → Excitatory
#     inhibitory_cluster = 1 - excitatory_cluster  # Other cluster is Inhibitory
    
#     # Plot PCA-clustered data
#     plt.figure(figsize=(10, 8))
#     plt.scatter(features_pca[:, 0], features_pca[:, 1], c=cluster_labels, cmap='viridis', marker='o')
#     plt.title('PCA of Neuron Waveform Features with GMM Clustering')
#     plt.xlabel('Principal Component 1')
#     plt.ylabel('Principal Component 2')
#     plt.colorbar(label='Cluster Label')
#     plt.grid()
#     well_path = os.path.dirname(wfs_output_path)
#     plt.savefig(os.path.join(os.path.dirname(well_path), 'PCA_GMM_Clustering.png'))
#     print(f"Saved PCA and GMM clustering plot to {os.path.join(os.path.dirname(well_path), 'PCA_GMM_Clustering.png')}")

#     # Assign classifications
#     for i, unit_id in enumerate(unit_id_list):
#         classification_dict[unit_id] = "excitatory" if cluster_labels[i] == excitatory_cluster else "inhibitory"

#     # Compute and print neuron ratios
#     excitatory_count = sum(1 for unit_id in classification_dict if classification_dict[unit_id] == "excitatory")
#     inhibitory_count = sum(1 for unit_id in classification_dict if classification_dict[unit_id] == "inhibitory")
#     print(f"Excitatory neurons: {excitatory_count}")
#     print(f"Inhibitory neurons: {inhibitory_count}")
#     print(f'Percent excitatory: {excitatory_count / (excitatory_count + inhibitory_count) * 100:.2f}%')
#     print(f'Percent inhibitory: {inhibitory_count / (excitatory_count + inhibitory_count) * 100:.2f}%')

#     return classification_dict

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

import spikeinterface.postprocessing as spost

def classify_neurons(network_data, **kwargs):
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
                plot = True
                #plot = False
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
