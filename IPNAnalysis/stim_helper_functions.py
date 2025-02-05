import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import h5py
from sklearn.manifold import TSNE

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
from scipy.signal import butter, filtfilt

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