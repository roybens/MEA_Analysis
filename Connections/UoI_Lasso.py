import numpy as np
from pyuoi.linear_model import UoI_Lasso
import matplotlib.pyplot as plt
from multiprocessing import Pool
###############################################################################
# Fit a UoI-Lasso model
# ---------------------
#
# UoI-Lasso can fit low bias model parameters with feature selectivity. We can
# evaluate the predictions of the model, compare the fit :math:`\beta`, and look
# at the fraction of false positive and false negatives.
data = np.load('spikearray.npz')


spike_matrix = data['data']

# Define the lag (number of previous time steps you want to use for prediction)
lag = 1  # This will use the previous time point to predict the next

# Extract the total number of neurons and time points
num_neurons, num_timepoints = spike_matrix.shape

# Function to run UoI-Lasso for a given neuron
def run_uoi_lasso(neuron_of_interest):
    print(f"Running UoI-Lasso for neuron {neuron_of_interest+1}/{num_neurons}...")

    # Define the response `y` for the neuron of interest (shifted by 1 time step to the right)
    y = spike_matrix[neuron_of_interest, lag:]  # Response for neuron_of_interest

    # Define the input feature matrix `X` (activity of all neurons at previous time step)
    X = spike_matrix[:, :-lag].T  # Spike activity of all neurons, lagged by 1 (transposed to shape (timepoints, neurons))

    # Remove the neuron_of_interest from the features matrix (we don't want it to predict itself)
    X = np.delete(X, neuron_of_interest, axis=1)  # Remove neuron_of_interest's activity from the predictors

    # Initialize UoI-Lasso model
    uoi_lasso = UoI_Lasso()

    # Fit the UoI-Lasso model to the data
    uoi_lasso.fit(X, y)

    # Get the coefficients (connectivity to neuron_of_interest from other neurons)
    coefficients = uoi_lasso.coef_

    # Insert a 0 for the neuron_of_interest itself (since we removed it earlier)
    return np.insert(coefficients, neuron_of_interest, 0)

# Function to handle parallel execution
def parallel_uoi_lasso():
    # Create a pool of workers, using as many CPU cores as available
    with Pool() as pool:
        # Run UoI-Lasso in parallel for all neurons
        connectivity_matrix = pool.map(run_uoi_lasso, range(num_neurons))
        
    # Convert the result into a NumPy array
    return np.array(connectivity_matrix)

# Run the parallel UoI-Lasso process
connectivity_matrix = parallel_uoi_lasso()

# Check the shape of the connectivity matrix
print("Connectivity Matrix Shape:", connectivity_matrix.shape)


# Plot the functional connectivity matrix
plt.figure(figsize=(10, 8))
plt.imshow(connectivity_matrix, cmap='coolwarm', interpolation='none')
plt.colorbar(label="Connection Strength")
plt.title('Functional Connectivity Matrix (UoI-Lasso)')
plt.xlabel('Predictor Neuron')
plt.ylabel('Target Neuron')

# Invert the y-axis so the first neuron is at the top
plt.gca().invert_yaxis()

# Save the figure as a PDF
plt.savefig("functional_connectivity_matrix.pdf", format="pdf", bbox_inches="tight")

# No plt.show() to avoid displaying the plot




