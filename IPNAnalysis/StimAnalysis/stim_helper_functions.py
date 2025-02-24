import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import h5py
from sklearn.manifold import TSNE

def filter_spikes_by_isi(spike_times, isi_threshold=0.01):
    """
    Filter spikes based on inter-spike-interval threshold
    
    Params:
        spike_times : array of spike timestamps in seconds
        isi_threshold : minimum allowed time between spikes in seconds (default: 10ms)
        
    Returns:
        filtered_spikes : array of spike times with bursts filtered out
    """
    if len(spike_times) < 2:
        return spike_times
        
    spike_times = np.sort(spike_times)
    isis = np.diff(spike_times)
    
    filtered_spikes = [spike_times[0]]
    
    # Add spikes only if they are separated by at least isi_threshold
    for i in range(1, len(spike_times)):
        if isis[i-1] >= isi_threshold:
            filtered_spikes.append(spike_times[i])
            
    return np.array(filtered_spikes)

def plot_spike_waveforms_subplots(waveforms, fs=10000, num_spikes=10):
    """
    Plots multiple spike waveforms in subplots within a single figure.

    Parameters:
    - waveforms: Dictionary where keys are spike times (in seconds) and values are waveform arrays.
    - fs: Sampling frequency (Hz).
    - num_spikes: Number of spikes to plot (default: 10).
    """

    # Limit number of spikes plotted
    spike_times = list(waveforms.keys())[:num_spikes]
    num_spikes = len(spike_times)

    fig, axes = plt.subplots(num_spikes, 1, figsize=(6, num_spikes * 2), sharex=True)
    
    if num_spikes == 1:
        axes = [axes]  # Ensure iterable when only one subplot

    for i, spike_time in enumerate(spike_times):
        waveform = waveforms[spike_time]
        time_axis = np.linspace(-len(waveform) / (2 * fs) * 1000, len(waveform) / (2 * fs) * 1000, len(waveform))  # Time in ms

        axes[i].plot(time_axis, waveform, color='black')
        axes[i].set_ylabel("Amplitude (µV)")
        axes[i].set_title(f"Spike at {spike_time:.4f} sec")
        axes[i].grid()

    plt.xlabel("Time (ms)")
    plt.show()

def plot_cluster_examples(spike_waveforms_dict, tsne_results, n_examples=3):
    """
    Plot example waveforms from different clusters to compare their shapes
    Params: 
        spike_waveforms_dict: dictionary of spike waveforms (from extract_spike_waveforms())
        tsne_results: t-SNE results (from tsne_spike_visualization())
        n_examples: number of examples to plot
    """
    # Convert dictionary to list for indexing
    waveforms = list(spike_waveforms_dict.values())
    
    # Define regions of interest based on your t-SNE plot
    regions = {
        'Main Cluster': (tsne_results[:, 0] < 20) & (tsne_results[:, 0] > -40) & (tsne_results[:, 1] < 30),
        'Vertical Chain': (tsne_results[:, 0] > 30) & (tsne_results[:, 0] < 50),
        'Top Cluster': (tsne_results[:, 1] > 30),
        'Far Right Outliers': (tsne_results[:, 0] > 50)
    }
    
    fig, axes = plt.subplots(len(regions), n_examples, figsize=(15, 10))
    
    for i, (name, mask) in enumerate(regions.items()):
        # Get indices for this region
        indices = np.where(mask)[0]
        
        # Randomly sample n_examples from this region
        sample_indices = np.random.choice(indices, size=min(n_examples, len(indices)), replace=False)
        
        for j, idx in enumerate(sample_indices):
            waveform = waveforms[idx]
            axes[i, j].plot(waveform)
            axes[i, j].set_title(f'{name}\nExample {j+1}')
            axes[i, 0].set_ylabel('Amplitude')
    
    plt.tight_layout()
    plt.show()

def tsne_spike_visualization(spike_waveforms_dict):
    """
    Visualize spike waveforms using t-SNE dimensionality reduction.
    
    Params:
    spike_waveforms_dict : dict Dictionary mapping spike times to waveform arrays
    """
    # convert dictionary of waveforms to a matrix
    waveforms_matrix = np.array(list(spike_waveforms_dict.values()))
    
    # t-SNE
    tsne = TSNE(n_components=2, perplexity=min(30, len(waveforms_matrix)-1), random_state=42)
    tsne_results = tsne.fit_transform(waveforms_matrix)
    
    # scatter plot
    plt.figure(figsize=(8,6))
    plt.scatter(tsne_results[:, 0], tsne_results[:, 1], alpha=0.5)
    plt.title('t-SNE Visualization of Spike Waveforms')
    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')
    plt.show()
    
    return tsne_results

def get_spike_cluster_times(tsne_results, spike_waveforms_dict):
    """
    Get spike times separated by cluster type from t-SNE analysis.
    Params:
        tsne_results: t-SNE results
        spike_waveform_dict: dictionary of spike waveforms
    Returns:
        dict : Dictionary containing lists of spike times (in seconds) for each cluster type
    """
    # Define regions
    regions = {
        'Main Cluster': (tsne_results[:, 0] < 20) & (tsne_results[:, 0] > -40) & (tsne_results[:, 1] < 30),
        'Vertical Chain': (tsne_results[:, 0] > 30) & (tsne_results[:, 0] < 50),
        'Top Cluster': (tsne_results[:, 1] > 30),
        'Far Right Outliers': (tsne_results[:, 0] > 50)
    }
    
    spike_times = np.array(list(spike_waveforms_dict.keys()))
    
    cluster_times = {}
    
    # Get times for each cluster type
    for region_name, region_mask in regions.items():
        cluster_times[region_name] = spike_times[region_mask]
    
    return cluster_times


def print_file_structure(self):
        with h5py.File(self.file_path, 'r') as h5file:
            def print_hdf5_structure(name, obj):
                print(name)

            h5file.visititems(print_hdf5_structure)

def plot_frequency_spectrum(trace, fs):
    """
    Plots the frequency spectrum of a given trace.
    trace: The signal (array of amplitudes).
    fs: Sampling frequency (in Hz).
    """
    # Calculate the FFT of the signal
    fft_result = np.fft.fft(trace)
    
    # Compute the corresponding frequencies
    freqs = np.fft.fftfreq(len(trace), d=1/fs)
    
    # Take the magnitude of the FFT (absolute value)
    fft_magnitude = np.abs(fft_result)
    
    # Only consider positive frequencies
    positive_freqs = freqs[freqs >= 0]
    positive_magnitude = fft_magnitude[freqs >= 0]
    
    # Plot the frequency spectrum
    plt.figure(figsize=(10, 6))
    plt.plot(positive_freqs, positive_magnitude)
    plt.title("Frequency Spectrum of the Signal")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Amplitude")
    plt.xlim([0, fs/2])  # Nyquist limit
    plt.grid()
    plt.show()


def bandpass_filter(trace, fs, lowcut, highcut, order=3):
    nyquist = fs / 2
    low = lowcut / nyquist
    high = highcut / nyquist
    b, a = butter(order, [low, high], btype='band')
    filtered_trace = filtfilt(b, a, trace)
    return filtered_trace

def polynomial_interpolation(trace, start_idx, end_idx, degree=3):
    # Define x and y points outside the artifact region
    x = np.concatenate((
        np.arange(0, start_idx),  # Points before
        np.arange(end_idx + 1, len(trace) - 1)   # Points after
    ))
    y = np.concatenate((
        trace[0: start_idx],
        trace[end_idx + 1 : len(trace) - 1]
    ))

    # Fit polynomial to surrounding points
    # poly_coeffs = np.polyfit(x, y, degree)

    
    spline = CubicSpline(x, y)


    # Interpolate the artifact region
    artifact_indices = np.arange(start_idx, end_idx + 1)
    # trace[artifact_indices] = np.polyval(poly_coeffs, artifact_indices)
    trace[artifact_indices] = spline(artifact_indices)

    # plt.plot(x, y, 'o', label="Surrounding Points")
    # plt.plot(artifact_indices, trace[artifact_indices], 'r-', label="Interpolated Region")
    # plt.legend()
    # plt.show()

    return trace

def compare_early_late_waveforms(spike_waveforms_dict, spike_times, early_threshold=60, late_threshold=90):
    """
    Compare spike waveforms from early (<60s) and late (>90s) time periods
    
    Parameters:
    -----------
    spike_waveforms_dict : dict
        Dictionary with spike times as keys and waveform arrays as values
    spike_times : array-like
        Array of spike times in seconds
    early_threshold : float
        Time threshold for early spikes (default: 60s)
    late_threshold : float
        Time threshold for late spikes (default: 90s)
        
    Returns:
    --------
    None (displays plot)
    """
    
    # Get early and late spike times
    early_spikes = spike_times[spike_times < early_threshold]
    late_spikes = spike_times[spike_times > late_threshold]
    
    # Get corresponding waveforms
    early_waveforms = [spike_waveforms_dict[t] for t in early_spikes]
    late_waveforms = [spike_waveforms_dict[t] for t in late_spikes]
    
    # Convert to numpy arrays
    early_waveforms = np.array(early_waveforms)
    late_waveforms = np.array(late_waveforms)
    
    # Calculate means and standard deviations
    early_mean = np.mean(early_waveforms, axis=0)
    early_std = np.std(early_waveforms, axis=0)
    late_mean = np.mean(late_waveforms, axis=0)
    late_std = np.std(late_waveforms, axis=0)
    
    # Create time axis (assuming waveform length matches your data)
    time_ms = np.linspace(-1, 1, len(early_mean))  # Adjust range as needed
    
    # Plot
    plt.figure(figsize=(12, 8))
    
    # Plot early spikes
    plt.fill_between(time_ms, 
                     early_mean - early_std, 
                     early_mean + early_std, 
                     alpha=0.3, 
                     color='blue',
                     label='Early Std')
    plt.plot(time_ms, early_mean, 'b-', label=f'Early Mean (n={len(early_spikes)})')
    
    # Plot late spikes
    plt.fill_between(time_ms, 
                     late_mean - late_std, 
                     late_mean + late_std, 
                     alpha=0.3, 
                     color='red',
                     label='Late Std')
    plt.plot(time_ms, late_mean, 'r-', label=f'Late Mean (n={len(late_spikes)})')
    
    plt.title('Comparison of Early (<60s) vs Late (>90s) Spike Waveforms')
    plt.xlabel('Time (ms)')
    plt.ylabel('Amplitude (µV)')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    # Print some basic statistics
    print(f"Number of early spikes: {len(early_spikes)}")
    print(f"Number of late spikes: {len(late_spikes)}")
    print(f"Early mean amplitude: {np.max(early_mean)-np.min(early_mean):.2f} µV")
    print(f"Late mean amplitude: {np.max(late_mean)-np.min(late_mean):.2f} µV")

