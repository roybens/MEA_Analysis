import json
import numpy as np
import math
import os
import fnmatch
import matplotlib.pyplot as plt
from scipy.signal import convolve, find_peaks
from scipy.stats import norm
from scipy.interpolate import interp1d

# --- Utility Functions ---

def detect_peaks(trace, peak_sign, abs_threshold):
    peaks_sample_inds = []
    peaks_chan_inds = []

    if peak_sign in ("pos", "both"):
        peaks, _ = find_peaks(trace, height=abs_threshold)
        peaks_sample_inds.extend(peaks)
        peaks_chan_inds.extend([0] * len(peaks))

    if peak_sign in ("neg", "both"):
        peaks, _ = find_peaks(-trace, height=abs_threshold)
        peaks_sample_inds.extend(peaks)
        peaks_chan_inds.extend([0] * len(peaks))

    return np.array(peaks_sample_inds), np.array(peaks_chan_inds)

def get_surrounding_coordinates(x, y, number):
    surrounding_coordinates = []
    number = math.sqrt(number)
    for i in range(-2, 3):
        for j in range(-2, 3):
            surrounding_coordinates.append((x + i, y + j))
    
    return surrounding_coordinates


def convert_int64_keys_to_ints(d):
    new_dict = {}
    for k, v in d.items():
        if isinstance(k, np.int64):
            k = int(k)
        if isinstance(v, dict):
            v = convert_int64_keys_to_ints(v)
        new_dict[k] = v
    return new_dict


def load_json(filename):
    d = {}
    with open(filename,'rb') as f:
        d = json.load(f)
    return d

def get_templates_with_same_channels(electrode_file):

    my_dict = load_json(electrode_file)

    templates_with_channel = {}

    for template_name, template_data in my_dict.items():
        for channel_name in template_data:
            if channel_name in templates_with_channel:
                templates_with_channel[channel_name].append(template_name)
            else:
                templates_with_channel[channel_name] = [template_name]
    
    same_channel_templates = [templates for channel, templates in templates_with_channel.items() if len(templates) > 1]

    return same_channel_templates

def empty_directory(directory_path):
    if os.path.exists(directory_path) and os.path.isdir(directory_path):
        # List all files and subdirectories in the directory
        file_list = os.listdir(directory_path)
        
        for filename in file_list:
            file_path = os.path.join(directory_path, filename)
            if os.path.isfile(file_path):
                # Remove files
                os.remove(file_path)
            elif os.path.isdir(file_path):
                # Remove subdirectories and their contents recursively
                empty_directory(file_path)
                os.rmdir(file_path)
        print(f"The directory '{directory_path}' has been emptied.")
    else:
        print(f"The directory '{directory_path}' does not exist.")



def find_files_with_subfolder(root_dir, file_name_pattern, subfolder_name):
    file_paths = []
    for dirpath, _, filenames in os.walk(root_dir):
        if subfolder_name in dirpath.split(os.path.sep):
            for filename in fnmatch.filter(filenames, file_name_pattern):
                file_paths.append(os.path.join(dirpath, filename))
    return file_paths

def isexists_folder_not_empty(folder_path):
    # Check if the folder exists
    if not os.path.exists(folder_path):
        print("sorting dir doesnt exist: making new one")
        os.makedirs(folder_path,0o777)
        return 0
    # Get the list of items (files and folders) in the specified directory
    items = os.listdir(folder_path)

    # If the folder is not empty, return True
    return len(items) > 0

def dumpdicttofile(data,filename):
    data = convert_int64_keys_to_ints(data)
    json_data = json.dumps(data,indent=4)

    with open(filename,'w') as fileptr:
        fileptr.write(json_data)

def save_json(file_path, data):
    def default_converter(o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return float(o)
        raise TypeError(f'Object of type {o.__class__.__name__} is not JSON serializable')

    print(f"[HELPER] Saving JSON to {file_path}...")
    with open(file_path, 'w') as f:
        json.dump(data, f, indent=4, default=default_converter)

# --- Burst Detection & Statistics ---

def detect_bursts_statistics(spike_times, isi_threshold):
    """
    Detect bursts and calculate stats with DEBUG prints and SAFE covariance.
    """
    print(f"[HELPER] Detecting bursts for {len(spike_times)} units (ISI thresh={isi_threshold}s)...")
    results = {}

    for unit, times in spike_times.items():
        # Handle empty spike trains
        if len(times) < 2:
            print(f"[HELPER WARNING] Unit {unit} has < 2 spikes. Skipping stats.")
            results[unit] = {
                "bursts": [], "mean_isi_within": np.nan, "cov_isi_within": np.nan,
                "mean_isi_outside": np.nan, "cov_isi_outside": np.nan,
                "mean_isi_all": np.nan, "cov_isi_all": np.nan,
                "isis_within_bursts": [], "isis_outside_bursts": [], "isis_all": []
            }
            continue

        # Step 1: Calculate ISIs
        isis = np.diff(times)
        
        # Step 2: Identify bursts
        burst_mask = isis < isi_threshold
        
        # Step 3: Find indices
        # Pad with False to handle edge cases
        padded_mask = np.hstack(([False], burst_mask, [False]))
        diffs = np.diff(padded_mask.astype(int))
        burst_starts = np.where(diffs == 1)[0]
        burst_ends = np.where(diffs == -1)[0]
        
        # Step 4: Group spikes
        bursts = [times[start:end + 1] for start, end in zip(burst_starts, burst_ends)]
        
        # Extract ISIs
        if len(burst_starts) > 0:
            isis_within_bursts = np.concatenate([isis[start:end] for start, end in zip(burst_starts, burst_ends)])
        else:
            isis_within_bursts = np.array([])
            
        isis_outside_bursts = isis[~burst_mask]

        # --- SAFE STATISTICS (Fixes RuntimeWarning) ---
        def safe_mean(arr):
            return np.mean(arr) if arr.size > 0 else np.nan

        def safe_cov(arr):
            # Covariance requires at least 2 data points to be valid
            # If size is 0 or 1, return NaN to avoid 'Degrees of freedom <= 0' warning
            return np.cov(arr) if arr.size > 1 else np.nan

        results[unit] = {
            "bursts": bursts,
            "mean_isi_within": safe_mean(isis_within_bursts),
            "cov_isi_within": safe_cov(isis_within_bursts),
            "mean_isi_outside": safe_mean(isis_outside_bursts),
            "cov_isi_outside": safe_cov(isis_outside_bursts),
            "mean_isi_all": safe_mean(isis),
            "cov_isi_all": safe_cov(isis),
            "isis_within_bursts": isis_within_bursts,
            "isis_outside_bursts": isis_outside_bursts,
            "isis_all": isis
        }

    return results

def plot_raster_with_bursts(ax, spike_times, bursts, sorted_units=None, title_suffix=""):
    """
    Plots raster with extensive debug checks for type mismatches.
    """
    print(f"[HELPER] Plotting raster. SpikeTimes keys: {len(spike_times)}. Bursts type: {type(bursts)}")
    
    # Debug: Check structure of 'bursts'
    if isinstance(bursts, list):
        print(f"[HELPER WARNING] 'bursts' passed as LIST (len={len(bursts)}). Converting to DICT based on spike_times keys.")
        # Attempt to map list back to keys if lengths match
        keys = list(spike_times.keys())
        if len(keys) == len(bursts):
            bursts = dict(zip(keys, bursts))
        else:
            print("[HELPER ERROR] Cannot map bursts list to units. Length mismatch!")
            return ax

    y_offset = 0
    units_to_plot = sorted_units if sorted_units else sorted(list(spike_times.keys()))
    
    for unit in units_to_plot:
        # Check integrity
        if unit not in spike_times:
            print(f"[HELPER SKIP] Unit {unit} missing from spike_times")
            continue
        if unit not in bursts:
            print(f"[HELPER SKIP] Unit {unit} missing from bursts dict")
            continue

        times = spike_times[unit]
        unit_bursts = bursts[unit]
        
        # Plot spikes
        ax.plot(times, np.ones_like(times) * y_offset, '|', color='royalblue', markersize=3, rasterized=True)
        
        # Plot bursts
        for i, burst in enumerate(unit_bursts):
            if len(burst) > 0:
                ax.plot(burst, np.ones_like(burst) * y_offset, '|', color='black', markersize=3, rasterized=True)
        
        y_offset += 1
    
    #ax.set_xlabel('Time (s)')
    ax.set_ylabel('Unit Index')
    ax.set_yticks([1, y_offset//2, y_offset-1]) 
    ax.set_title(f'Raster Plot {title_suffix}')
    
    print(f"[HELPER] Raster plot complete. Plotted {y_offset} units.")
    return ax

def plot_network_activity(ax, SpikeTimes, min_peak_distance=1.0, binSize=0.1, gaussianSigma=0.16, thresholdBurst=1.2):
    print("[HELPER] Calculating Network Activity...")
    
    # Flatten all spikes
    all_spikes = []
    for times in SpikeTimes.values():
        all_spikes.extend(times)
    
    if not all_spikes:
        print("[HELPER WARNING] No spikes found in any unit. Skipping network plot.")
        return ax, {}

    all_spikes = np.array(all_spikes)
    all_spikes.sort()
    
    start_time = np.floor(all_spikes[0])
    end_time = np.ceil(all_spikes[-1])
    
    if end_time - start_time < binSize:
        print("[HELPER WARNING] Recording duration too short for binning.")
        return ax, {}

    # Binning
    bins = np.arange(start_time, end_time + binSize, binSize)
    binned_counts, _ = np.histogram(all_spikes, bins=bins)
    
    # Convert to Firing Rate (Hz)
    # Rate = Count / BinSize / NumUnits
    num_units = len(SpikeTimes)
    firing_rate_raw = binned_counts / binSize / num_units

    # Gaussian Smoothing
    sigma_bins = gaussianSigma / binSize
    window_size = int(6 * sigma_bins)
    if window_size % 2 == 0: window_size += 1
    
    # Create Gaussian window
    x = np.linspace(-3*sigma_bins, 3*sigma_bins, window_size)
    kernel = np.exp(-0.5 * x**2 / sigma_bins**2)
    kernel /= kernel.sum() # Normalize
    
    firing_rate_smooth = convolve(firing_rate_raw, kernel, mode='same')
    
    # Time vector for plotting (centers of bins)
    time_vector = bins[:-1] + binSize/2

    # Plot
    ax.plot(time_vector, firing_rate_smooth, color='royalblue')
    ax.set_ylabel('Avg Firing Rate [Hz]')
    #ax.set_xlabel('Time [s]')
    ax.set_title('Population Firing Rate')
    ax.set_xlim([start_time, end_time])

    # Peak Detection
    peaks, _ = find_peaks(firing_rate_smooth, prominence=0.5, distance=int(min_peak_distance/binSize))
    
    burst_peaks_t = time_vector[peaks]
    burst_peaks_val = firing_rate_smooth[peaks]
    
    ax.plot(burst_peaks_t, burst_peaks_val, 'or', label='Bursts')
    if len(peaks) > 0: ax.legend()

    # Metrics
    intervals = np.diff(burst_peaks_t)
    
    stats = {
        "Number_Bursts": len(peaks),
        "mean_IBI": np.mean(intervals) if len(intervals) > 0 else np.nan,
        "cov_IBI": np.cov(intervals) if len(intervals) > 1 else np.nan,
        "mean_Burst_Peak": np.mean(burst_peaks_val) if len(burst_peaks_val) > 0 else np.nan,
        "cov_Burst_Peak": np.cov(burst_peaks_val) if len(burst_peaks_val) > 1 else np.nan
    }
    
    print(f"[HELPER] Network Activity: Found {len(peaks)} bursts.")
    return ax, stats