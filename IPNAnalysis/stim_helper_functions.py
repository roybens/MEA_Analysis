import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
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