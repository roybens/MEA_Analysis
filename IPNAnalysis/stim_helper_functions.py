import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import h5py
from sklearn.manifold import TSNE

def plot_spike_waveforms(spike_waveforms):
    """
    Plots all extracted spike waveforms to verify the alignment and window size.

    Parameters:
    - spike_waveforms: (num_spikes, window_size) array of extracted spike waveforms.
    """

    plt.figure(figsize=(10, 6))

    for waveform in spike_waveforms:
        plt.plot(waveform, alpha=0.3, color="black")  # Low opacity for overlapping visualization

    plt.xlabel("Time (samples)")
    plt.ylabel("Amplitude (mV)")
    plt.title("Extracted Spike Waveforms")
    plt.show()

def extract_spike_waveforms(recording_trace, spike_indices, window_size=50):
    """
    Extracts spike waveforms from a continuous recording trace.

    Parameters:
    - recording_trace: 1D array of the voltage signal.
    - spike_indices: List or array of indices where spikes occur.
    - window_size: Total number of samples per waveform.

    Returns:
    - waveforms: 2D array (num_spikes, window_size) of extracted waveforms.
    """

    half_window = window_size // 2
    waveforms = []

    for idx in spike_indices:
        if idx - half_window < 0 or idx + half_window >= len(recording_trace):
            continue  # Skip spikes near the edges

        waveform = recording_trace[idx - half_window : idx + half_window]
        waveforms.append(waveform)

    return np.array(waveforms)

def tsne_spike_visualization(spike_waveforms):
    """
    Uses t-SNE to visualize spike waveform clustering.

    Parameters:
    - spike_waveforms: (num_spikes, waveform_length) array.
    """

    tsne = TSNE(n_components=2, perplexity=30, random_state=42)
    tsne_results = tsne.fit_transform(spike_waveforms)

    plt.figure(figsize=(8,6))
    plt.scatter(tsne_results[:, 0], tsne_results[:, 1], alpha=0.5)
    plt.xlabel("t-SNE Component 1")
    plt.ylabel("t-SNE Component 2")
    plt.title("t-SNE Visualization of Spike Waveforms")
    plt.show()


def print_file_structure(self):
        with h5py.File(self.file_path, 'r') as h5file:
            def print_hdf5_structure(name, obj):
                print(name)

            h5file.visititems(print_hdf5_structure)

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