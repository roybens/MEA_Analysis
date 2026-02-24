import os
import json
import logging
from datetime import datetime
import numpy as np
from pyuoi.linear_model import UoI_Lasso
import matplotlib.pyplot as plt
from multiprocessing import Pool
import time
import pandas as pd
import networkx as nx
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from functools import reduce
from sklearn.linear_model import Lasso
import argparse 
#from UoIVarLassoMicehal import UoI_VAR_LASSO

 # For handling arguments
###############################################################################



# Define output folder paths
base_output_dir = os.path.join(os.getenv("SCRATCH", "."), "MEA_Analysis/Connections/output")
job_id = os.getenv("SLURM_JOB_ID", "local_run")
output_dir = os.path.join(base_output_dir, job_id)
os.makedirs(output_dir, exist_ok=True)

# Set up logging
log_file = os.path.join(output_dir, f"log_{job_id}.txt")
logging.basicConfig(
    filename=log_file,
    filemode="w",
    format="%(asctime)s - %(levelname)s - %(message)s",
    level=logging.INFO
)

# Log start
logging.info("Job started.")



# Parse command-line arguments
parser = argparse.ArgumentParser(description="Run UoI Lasso or VAR Lasso on input data.")
parser.add_argument("--data_path", required=True, help="Path to the directory containing the input data.")
args = parser.parse_args()


try:

     # Define parameters
    start_index = 0
    end_index = 6000
    lag = 1  # Number of previous time steps to use for prediction
    method_choice = "LASSO"  # Options: "LASSO" or "VAR_LASSO"

    # UoI-Lasso parameters (important parameters to log)
    uoi_params = {
        "n_boots_sel": 10,
        "n_boots_est": 10,
        "selection_frac": 0.9,
        "estimation_frac": 0.9,
        "n_lambdas": 48,
        "stability_selection": 1.0
    }
    logging.info(f"Method: {method_choice}")
    logging.info(f"UoI-Lasso parameters: {json.dumps(uoi_params, indent=4)}")

    start_time = time.time()


    # Load data
    data_path = args.data_path
    logging.info(f"Loading data from {data_path}.")
    data = np.load(data_path,allow_pickle=True).item()
    
    logging.info("Data loaded successfully.")

    num_neurons = len(data)

    duration = 300
    target_fs =1000
        # Calculate the number of bins at the target sampling frequency
    num_bins = int(duration * target_fs)
    
    # Create an empty binary matrix for the downsampled data
    spike_matrix = np.zeros((num_neurons, num_bins), dtype=int)
    
    # Define the bin edges based on the target sampling frequency
    bin_edges = np.linspace(0, duration, num_bins + 1)  # num_bins + 1 to include the last edge
    
    # Iterate over each neuron and populate the matrix
    for neuron_id, spikes in data.items():
        # Digitize the spike times into bins
        bin_indices = np.digitize(spikes, bin_edges) - 1  # Convert to 0-based index
        bin_indices = bin_indices[(bin_indices >= 0) & (bin_indices < num_bins)]  # Remove out-of-range indices
        
        # Mark the bins where spikes occurred
        spike_matrix[neuron_id, bin_indices] = 1


    spike_matrix = spike_matrix[:,:6000]
    #spike_matrix = data['data'][:,:6000]

    # Define the lag (number of previous time steps you want to use for prediction)
    lag = 1  # This will use the previous time point to predict the next

    # Extract the total number of neurons and time points
    num_neurons, num_timepoints = spike_matrix.shape
    logging.info(f"Spike matrix dimensions: {num_neurons} neurons, {num_timepoints} timepoints.")

    if method_choice == "LASSO":
        # UoI-Lasso function
        def run_uoi_lasso(neuron_of_interest):
            try:
                logging.info(f"Running UoI-Lasso for neuron {neuron_of_interest+1}/{num_neurons}...")

                # Define response and predictors
                y = spike_matrix[neuron_of_interest, lag:]
                X = spike_matrix[:, :-lag].T
                X = np.delete(X, neuron_of_interest, axis=1)

                # Initialize and fit UoI-Lasso
                uoi_lasso = UoI_Lasso(
                    n_boots_sel=uoi_params["n_boots_sel"],
                    n_boots_est=uoi_params["n_boots_est"],
                    selection_frac=uoi_params["selection_frac"],
                    estimation_frac=uoi_params["estimation_frac"],
                    n_lambdas=uoi_params["n_lambdas"],
                    stability_selection=uoi_params["stability_selection"],
                    estimation_score='r2',
                    max_iter=1000,
                    tol=1e-4
                )
                uoi_lasso.fit(X, y)

                # Get coefficients
                coefficients = uoi_lasso.coef_
                return np.insert(coefficients, neuron_of_interest, 0)
            except Exception as e:
                logging.error(f"Error running UoI-Lasso for neuron {neuron_of_interest}: {e}")
                raise

        # Parallel execution
        def parallel_uoi_lasso():
            with Pool() as pool:
                return np.array(pool.map(run_uoi_lasso, range(num_neurons)))

        # Run UoI-Lasso
        connectivity_matrix = parallel_uoi_lasso()
        logging.info(f"Connectivity matrix computed. Shape: {connectivity_matrix.shape}")
    """
    elif method_choice == "VAR_LASSO":
        # Define the VAR-LASSO function
        def create_lagged_data(df, lag_order):
            n_obs, n_var = df.shape
            X = np.zeros((n_obs - lag_order, n_var * lag_order))
            y = df.iloc[lag_order:].values  # Response matrix
            for t in range(lag_order, n_obs):
                for lag in range(lag_order):
                    X[t - lag_order, lag * n_var:(lag + 1) * n_var] = df.iloc[t - lag - 1].values
            return X, y

        def run_uoi_var_lasso():
            df_spikes = pd.DataFrame(spike_matrix.T)
            lag_order = 1
            X, y = create_lagged_data(df_spikes, lag_order)

            lambda_grid = np.linspace(0.0001, 0.01, 50)
            results = UoI_VAR_LASSO(
                df=df_spikes,
                reg_path=lambda_grid,
                B1=50,  # Bootstraps for selection
                d=lag_order,
                intercept=False,
                method='OLS',
                n_blocks=10,
                block_length=100,
                B2=50,  # Bootstraps for estimation
                gamma=1
            )
            n_neurons = df_spikes.shape[1]
            return results['Beta_hat'].reshape((n_neurons, n_neurons))

        # Run UoI-VAR_LASSO
        connectivity_matrix = run_uoi_var_lasso()
        logging.info(f"Connectivity matrix computed. Shape: {connectivity_matrix.shape}")
    """
    # Save connectivity matrix
    connectivity_matrix_file = os.path.join(
        output_dir, f'connectivity_matrix_{method_choice}_nbs{uoi_params["n_boots_sel"]}_nbe{uoi_params["n_boots_est"]}_nl{uoi_params["n_lambdas"]}_ss{uoi_params["stability_selection"]}.npy'
    )
    np.save(connectivity_matrix_file, connectivity_matrix)
    logging.info(f"Connectivity matrix saved to {connectivity_matrix_file}.")

    # Plot the functional connectivity matrix
    plt.figure(figsize=(10, 8))
    plt.imshow(connectivity_matrix, cmap='coolwarm', interpolation='none')
    plt.colorbar(label="Connection Strength")
    plt.title(f'Functional Connectivity Matrix ({method_choice})')
    plt.xlabel('Predictor Neuron')
    plt.ylabel('Target Neuron')
    plt.gca().invert_yaxis()
    plot_file = os.path.join(
        output_dir, f"functional_connectivity_matrix_{method_choice}_nbs{uoi_params['n_boots_sel']}_nbe{uoi_params['n_boots_est']}_nl{uoi_params['n_lambdas']}_ss{uoi_params['stability_selection']}.pdf"
    )
    plt.savefig(plot_file, format="pdf", bbox_inches="tight")
    plt.close()
    logging.info(f"Functional connectivity matrix plot saved to {plot_file}.")

    elapsed_time = time.time() - start_time
    logging.info(f"Job completed successfully in {elapsed_time:.2f} seconds.")

except Exception as e:
    logging.error("Job failed.", exc_info=True)
    raise

