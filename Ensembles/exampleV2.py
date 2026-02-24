import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from concurrent.futures import ProcessPoolExecutor



def compute_histogram(spikes, bins):
    """
    Compute histogram for a given neuron's spike times.
    Parameters:
        spikes (array): Spike times for a neuron.
        bins (array): Bin edges for the histogram.
    Returns:
        spike_counts (array): Histogram spike counts.
    """
    spike_counts, _ = np.histogram(spikes, bins=bins)
    return spike_counts

def gaussian(x, amp, mean, std):
    """Gaussian function."""
    return amp * np.exp(-((x - mean) ** 2) / (2 * std ** 2))

def multi_gaussian(x, *params):
    """
    Sum of multiple Gaussian functions.
    params: [amp1, mean1, std1, amp2, mean2, std2, ...]
    """
    n_gaussians = len(params) // 3
    result = np.zeros_like(x)
    for i in range(n_gaussians):
        amp, mean, std = params[i * 3:(i + 1) * 3]
        result += gaussian(x, amp, mean, std)
    return result

def fit_gaussians(bin_centers, spike_counts):
    """
    Fit Gaussians to the histogram data.
    """
    # Detect peaks in the histogram
    peaks, _ = find_peaks(spike_counts, height=1, distance=5)
    
    # Prepare initial guesses for Gaussian fitting
    initial_params = []
    for peak in peaks:
        amp = spike_counts[peak]
        mean = bin_centers[peak]
        std = 2  # Initial guess for standard deviation
        initial_params.extend([amp, mean, std])
    
    # Fit multi-Gaussian model
    try:
        popt, _ = curve_fit(multi_gaussian, bin_centers, spike_counts, p0=initial_params)
    except RuntimeError:
        popt = []
    
    return popt, peaks

def plot_histograms_to_pdf(spike_times, bin_size=0.1, duration=300, output_pdf="neuron_histograms_gaussians.pdf"):
    """
    Generate histograms for all neurons in spike_times, fit Gaussians, and save to a PDF.
    
    Parameters:
        spike_times (dict): Keys are neuron IDs, values are arrays of spike times (seconds).
        bin_size (float): Width of bins for the histogram in seconds.
        duration (float): Total duration of the recording in seconds.
        output_pdf (str): Output PDF file name.
    """
    bins = np.arange(0, duration + bin_size, bin_size)
    neuron_ids = list(spike_times.keys())
    neuron_spikes = [spike_times[nid] for nid in neuron_ids]

    # Compute histograms in parallel
    print("Computing histograms...")
    with ProcessPoolExecutor() as executor:
        histograms = list(tqdm(executor.map(compute_histogram, neuron_spikes, [bins]*len(neuron_ids)), 
                               desc="Processing neurons", unit="neuron", total=len(neuron_ids)))

    # Save histograms and Gaussian fits to PDF
    print("Generating PDF...")
    with PdfPages(output_pdf) as pdf:
        for neuron_id, spike_counts in zip(neuron_ids, histograms):
            bin_centers = bins[:-1] + bin_size / 2

            # Fit Gaussians to the histogram
            popt, peaks = fit_gaussians(bin_centers, spike_counts)

            # Plot histogram
            plt.figure(figsize=(8, 4))  # Adjusted for Gaussian fits
            plt.bar(bin_centers, spike_counts, width=bin_size, alpha=0.6, color="blue", label="Spike Histogram")
            
            # Plot Gaussian fits
            if popt:
                fit_curve = multi_gaussian(bin_centers, *popt)
                plt.plot(bin_centers, fit_curve, color="red", label="Gaussian Fit")
                for i in range(len(popt) // 3):
                    amp, mean, std = popt[i * 3:(i + 1) * 3]
                    plt.plot(bin_centers, gaussian(bin_centers, amp, mean, std), linestyle="--", label=f"Gaussian {i+1}")
            
            # Plot peak markers
            #plt.scatter(bin_centers[peaks], spike_counts[peaks], color="black", zorder=5, label="Peaks")
            
            # Customize plot
            plt.title(f"Neuron {neuron_id} Spike Histogram with Gaussian Fits")
            plt.xlabel("Time (s)")
            plt.ylabel("Spike Count")
            plt.xlim(0, duration)
            plt.ylim(0, max(spike_counts) + 1)
            plt.legend()
            plt.grid(False)

            # Save the figure to the PDF
            pdf.savefig()
            plt.close()

    print(f"Histograms with Gaussian fits saved to {output_pdf}")


# Load the .npz file
data = np.load('../AnalyzedData/CDKL5_R59X_SingleBrainPerChip/CDKL5-R59X_SingleBrainPerChip_11252024_PS/241208/M07038/Network/000027/well000/spike_times.npy',allow_pickle=True).item()



# Define the total recording duration in seconds
total_recording_duration = 300  # Example: 5 minutes (300 seconds)

# Calculate firing rates for all neurons
firing_rates = {}
for neuron_id in range(len(data)):
    spike_times = data[neuron_id]  # Get spike times for this neuron
    num_spikes = len(spike_times)  # Number of spikes
    firing_rate = num_spikes / total_recording_duration  # Firing rate in Hz
    firing_rates[neuron_id] = firing_rate  # Store firing rate

spike_times=data



plot_histograms_to_pdf(spike_times, bin_size=0.01, duration=300, output_pdf="neuron_histograms_gaussians_V2.pdf")