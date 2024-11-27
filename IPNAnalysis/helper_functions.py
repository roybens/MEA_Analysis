import json
import numpy as np
import math
import os
import fnmatch
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve, find_peaks
from scipy.stats import norm


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
    """
    Recursively convert all numpy.int64 keys in a dictionary to integers.
    """
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

def dumpdicttofile(data,filename):
    data = convert_int64_keys_to_ints(data)
    json_data = json.dumps(data,indent=4)

    with open(filename,'w') as fileptr:
        fileptr.write(json_data)

def get_key_by_value(d, value):
    """
    Useful to get the templates associated with a extreme electrode
    """
    keys = [k for k, v in d.items() if v == value]
    return keys if keys else None



def inarrxnotinarry(arr1,arr2):
    print(f" in array1 not in array2 :{[x for x in arr1 if x not in arr2]} ")
    print(f" in array2 not in array1 :{[x for x in arr2 if x not in arr1]} ")


def flatten(lst):
    flat_list = []
    for item in lst:
        if isinstance(item, list):
            flat_list.extend(flatten(item))
        else:
            flat_list.append(item)
    return flat_list

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

def detect_bursts_statistics(spike_times, isi_threshold):
    """
    Detect bursts in a spike train and calculate burst statistics.

    Parameters  
    ----------
    spike_times : dict  
        A dictionary of spike times for each unit.
    isi_threshold : float
        The threshold for detecting bursts (in seconds).

    Returns

    -------
    bursts : dict   
        A dictionary of burst times for each unit.
    
    """
# Dictionary to store results
    results = {}

    for unit, times in spike_times.items():
        
        # Step 1: Calculate the ISIs
        isis = np.diff(times)
        
        # Step 2: Identify where ISIs are below the threshold
        burst_mask = isis < isi_threshold
        
        # Step 3: Find the start and end indices of bursts using diff and boolean indexing
        burst_starts = np.where(np.diff(np.insert(burst_mask, 0, False).astype(int)) == 1)[0]
        burst_ends = np.where(np.diff(np.append(burst_mask, False).astype(int)) == -1)[0]
        
        # Step 4: Group spikes into bursts
        bursts = [times[start:end + 1] for start, end in zip(burst_starts, burst_ends)]
        
        # Extract ISIs within bursts using list comprehension and vectorized operations
        isis_within_bursts = np.concatenate([isis[start:end] for start, end in zip(burst_starts, burst_ends)]) if len(burst_starts) > 0 else np.array([])
        isis_outside_bursts = isis[~burst_mask]

        # Calculate statistics using vectorized numpy operations
        mean_isi_within = np.mean(isis_within_bursts) if isis_within_bursts.size > 0 else np.nan
        cov_isi_within = np.cov(isis_within_bursts) if isis_within_bursts.size > 0 else np.nan
        
        mean_isi_outside = np.mean(isis_outside_bursts) if isis_outside_bursts.size > 0 else np.nan
        cov_isi_outside = np.cov(isis_outside_bursts) if isis_outside_bursts.size > 0 else np.nan
        
        mean_isi_all = np.mean(isis) if isis.size > 0 else np.nan
        cov_isi_all = np.cov(isis) if isis.size > 0 else np.nan

        results[unit] = {
            "bursts": bursts,
            "mean_isi_within": mean_isi_within,
            "cov_isi_within": cov_isi_within,
            "mean_isi_outside": mean_isi_outside,
            "cov_isi_outside": cov_isi_outside,
            "mean_isi_all": mean_isi_all,
            "cov_isi_all": cov_isi_all,
            "isis_within_bursts": isis_within_bursts,
            "isis_outside_bursts": isis_outside_bursts,
            "isis_all": isis
        }

    return results
    






# Plot raster plot with bursts highlighted for all units
def plot_raster_with_bursts(ax, spike_times, bursts, sorted_units=None, title_suffix=""):
    y_offset = 0
    units_to_plot = sorted_units if sorted_units else list(spike_times.keys())
    
    for unit in units_to_plot:
        times = spike_times[unit]
        unit_bursts = bursts[unit]
        
        # Plot all spike times for this unit
        ax.plot(times, np.ones_like(times) + y_offset, '|', color='royalblue', markersize=1, rasterized=True)
        
        # Highlight bursts for this unit
        for burst in unit_bursts:
            ax.plot(burst, np.ones_like(burst) + y_offset, '|', color='black', markersize=1, rasterized=True)
        
        y_offset += 1  # Move to the next unit
    
    ax.set_xlabel('Time (s)')
    num_units = len(units_to_plot)
    yticks = [1, num_units//3, 2*num_units//3, num_units]
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks)
    ax.set_ylabel('Units')
    ax.set_title(f'Raster Plot with Bursts Highlighted {title_suffix}')

    return ax

def plot_network_activity(ax,SpikeTimes, min_peak_distance=1.0, binSize=0.1, gaussianSigma=0.16, thresholdBurst=1.2, figSize=(10, 6)):
    relativeSpikeTimes = []
    units = 0
    for unit_id, spike_times in SpikeTimes.items():
        relativeSpikeTimes.extend(spike_times) 
        units += 1 # Set the first spike time to 0
    
    relativeSpikeTimes = np.array(relativeSpikeTimes)
    relativeSpikeTimes.sort() # Convert to NumPy array

    # Step 1: Bin all spike times into small time windows
    timeVector = np.arange(min(relativeSpikeTimes), max(relativeSpikeTimes), binSize)  # Time vector for binning
    binnedTimes, _ = np.histogram(relativeSpikeTimes, bins=timeVector)  # Bin spike times
    binnedTimes = np.append(binnedTimes, 0)  # Append 0 to match MATLAB's binnedTimes length

    # Step 2: Smooth the binned spike times with a Gaussian kernel
    kernelRange = np.arange(-3*gaussianSigma, 3*gaussianSigma + binSize, binSize)  # Range for Gaussian kernel
    kernel = norm.pdf(kernelRange, 0, gaussianSigma)  # Gaussian kernel
    kernel *= binSize  # Normalize kernel by bin size
    firingRate = convolve(binnedTimes, kernel, mode='same') / binSize  # Convolve and normalize by bin size
    firingRate = firingRate / units  # Convert to Hz

    # Create a new figure with a specified size (width, height)
    #fig, ax = plt.subplots(figsize=figSize)

    # Plot the smoothed network activity
    ax.plot(timeVector, firingRate, color='royalblue')
    # Restrict the plot to the first and last 100 ms
    ax.set_xlim([min(relativeSpikeTimes), max(relativeSpikeTimes)])
    ax.set_ylim([min(firingRate)*0.8, max(firingRate)*1.2])  # Set y-axis limits to min and max of firingRate
    ax.set_ylabel('Firing Rate [Hz]')
    ax.set_xlabel('Time [ms]')
    ax.set_title('Network Activity', fontsize=11)

    # Step 3: Peak detection on the smoothed firing rate curve
    rmsFiringRate = np.sqrt(np.mean(firingRate**2))  # Calculate RMS of the firing rate

    peaks, properties = find_peaks(firingRate, prominence=1, distance=min_peak_distance)  # Find peaks above the threshold
    print(properties.keys())
    burstPeakTimes = timeVector[peaks]  # Convert peak indices to times
    #burstPeakValues = properties['prominences']  # Get the peak values
    burstPeakValues = firingRate[peaks]  # Get the peak values

    # Calculate the derivative
    derivative = np.diff(firingRate)
    derivative = np.append(derivative, 0)  # Match length with the signal

    # Identify zero crossings
    zero_crossings = np.where(np.diff(np.sign(derivative)))[0]

    # Calculate widths using zero crossings
    widths = []
    for peak in peaks:
        # Find the nearest zero crossings before and after the peak
        left_zero = zero_crossings[zero_crossings < peak][-1]
        right_zero = zero_crossings[zero_crossings > peak][0]
        width = right_zero - left_zero
        widths.append(width)

    mean_burst_duration = np.mean(widths) if widths else np.nan
    #Calculate the ISIs between spiketimes
    # Calculate the intervals between consecutive peaks
    intervals = np.diff(burstPeakTimes)

    # Calculate the mean interval between peaks
    mean_interburstinterval = np.mean(intervals)
    
    covariance_interburstinterval = np.cov(intervals)
    # Count the number of peaks
    num_peaks = len(peaks)

    # Calculate the mean peak height
    mean_peak_height = np.mean(burstPeakValues)
    cov_peak_height = np.cov(burstPeakValues)

    # # Print the results
    # print("Mean Interval between Peaks:", mean_interburstinterval)
    # print("Number of Peaks:", num_peaks)
    # print("Mean Peak Height:", mean_peak_height)

    # #calculate the mean isi and the variance of the isi
    # mean_isi =np.mean([SpikeTimes.values()])

    spikecounts= [len(x) for x in SpikeTimes.values()]
    var_spikecounts = np.var(spikecounts)
    mean_spikecounts = np.mean(spikecounts)
    fanofact = var_spikecounts/mean_spikecounts

    network_data = {
        "Number_Bursts": num_peaks,
        "mean_IBI": mean_interburstinterval,
        "cov_IBI": covariance_interburstinterval,
        "mean_Burst_Peak": mean_peak_height,
        "mean_Burst_Duration" : mean_burst_duration,
        "cov_Burst_Peak": cov_peak_height,
        "fano_factor": fanofact
    }   


    # Plot the threshold line and burst peaks
    #ax.plot(np.arange(timeVector[-1]), thresholdBurst * rmsFiringRate * np.ones(np.ceil(timeVector[-1]).astype(int)), color='gray')
    ax.plot(burstPeakTimes, burstPeakValues, 'or')  # Plot burst peaks as red circles    

    return ax,network_data


def save_json(file_path, data):
    def default_converter(o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(f'Object of type {o.__class__.__name__} is not JSON serializable')

    with open(file_path, 'w') as f:
        json.dump(data, f, indent=4, default=default_converter)