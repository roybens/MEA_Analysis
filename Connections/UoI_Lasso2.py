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
try:
    # Define parameters
    start_index = 0
    end_index = 6000
    lag = 1  # Number of previous time steps to use for prediction

    # UoI-Lasso parameters (important parameters to log)
    uoi_params = {
        "n_boots_sel": 10,
        "n_boots_est": 10,
        "selection_frac": 0.9,
        "estimation_frac": 0.9,
        "n_lambdas": 48,
        "stability_selection": 1.0
    }
    logging.info(f"UoI-Lasso parameters: {json.dumps(uoi_params, indent=4)}")

    start_time = time.time()

    # Load data
    data_path = '/pscratch/sd/m/mpatil1/MEA_Analysis/Connections/data/downsampled_spike_data.npz'
    logging.info(f"Loading data from {data_path}.")
    data = np.load(data_path)
    spike_matrix = data['spikes'][:, start_index:end_index]
    logging.info("Data loaded successfully.")

    # Get number of neurons and time points
    num_neurons, num_timepoints = spike_matrix.shape
    logging.info(f"Spike matrix dimensions: {num_neurons} neurons, {num_timepoints} timepoints.")

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

    # Save connectivity matrix
    connectivity_matrix_file = os.path.join(
        output_dir, f'connectivity_matrix_nbs{uoi_params["n_boots_sel"]}_nbe{uoi_params["n_boots_est"]}_nl{uoi_params["n_lambdas"]}_ss{uoi_params["stability_selection"]}.npy'
    )
    np.save(connectivity_matrix_file, connectivity_matrix)
    logging.info(f"Connectivity matrix saved to {connectivity_matrix_file}.")

    # Plot the functional connectivity matrix
    plt.figure(figsize=(10, 8))
    plt.imshow(connectivity_matrix, cmap='coolwarm', interpolation='none')
    plt.colorbar(label="Connection Strength")
    plt.title('Functional Connectivity Matrix (UoI-Lasso)')
    plt.xlabel('Predictor Neuron')
    plt.ylabel('Target Neuron')
    plt.gca().invert_yaxis()
    plot_file = os.path.join(
        output_dir, f"functional_connectivity_matrix_nbs{uoi_params['n_boots_sel']}_nbe{uoi_params['n_boots_est']}_nl{uoi_params['n_lambdas']}_ss{uoi_params['stability_selection']}.pdf"
    )
    plt.savefig(plot_file, format="pdf", bbox_inches="tight")
    plt.close()
    logging.info(f"Functional connectivity matrix plot saved to {plot_file}.")

    # # Additional analysis (Neuron connection visualization and counts)
    # df = pd.read_excel('./data/unitinfos.xlsx')
    # x_coords = df['LocX']
    # y_coords = df['LocY']
    with open('./data/gen_1_cand_29_data.json', 'r') as file:
        data = json.load(file)
    
    # Extract neuron locations from `net.cells`
    def extract_locations_from_cells(data):
        if "cells" in data["net"]:
            cells = data["net"]["cells"]
            x_coords = []
            y_coords=[]
            for cell in cells:
                x_coords.append(cell.get("tags", {}).get("x"))
                y_coords.append(cell.get("tags", {}).get("y"))
            return x_coords, y_coords
        else:
            print("No `cells` section found in the JSON data.")
            return []
   
    x_coords, y_coords = extract_locations_from_cells(data)
    G = nx.DiGraph()
    for i, (x, y) in enumerate(zip(x_coords, y_coords)):
        G.add_node(i, pos=(x, y))
    for i in range(num_neurons):
        for j in range(num_neurons):
            weight = connectivity_matrix[i, j]
            if weight > 0:
                G.add_edge(j, i, weight=weight, color='blue')
            elif weight < 0:
                G.add_edge(j, i, weight=weight, color='red')
    pos = {i: (x, y) for i, (x, y) in enumerate(zip(x_coords, y_coords))}
    input_excitatory = {i: 0 for i in range(num_neurons)}
    output_excitatory = {i: 0 for i in range(num_neurons)}
    input_inhibitory = {i: 0 for i in range(num_neurons)}
    output_inhibitory = {i: 0 for i in range(num_neurons)}
    for i, j, data in G.edges(data=True):
        weight = data['weight']
        if weight > 0:
            input_excitatory[j] += 1
            output_excitatory[i] += 1
        elif weight < 0:
            input_inhibitory[j] += 1
            output_inhibitory[i] += 1
    df_counts = pd.DataFrame({
        "Neuron": list(range(num_neurons)),
        "Incoming Excitatory": [input_excitatory[i] for i in range(num_neurons)],
        "Outgoing Excitatory": [output_excitatory[i] for i in range(num_neurons)],
        "Incoming Inhibitory": [input_inhibitory[i] for i in range(num_neurons)],
        "Outgoing Inhibitory": [output_inhibitory[i] for i in range(num_neurons)]
    })
    excel_file = os.path.join(
        output_dir, f"neuron_connection_counts_nbs{uoi_params['n_boots_sel']}_nbe{uoi_params['n_boots_est']}_nl{uoi_params['n_lambdas']}_ss{uoi_params['stability_selection']}.xlsx"
    )
    df_counts.to_excel(excel_file, index=False)
    logging.info(f"Neuron connection counts saved to {excel_file}.")

    viz_file = os.path.join(
        output_dir, f"neuron_connections_visualizations_nbs{uoi_params['n_boots_sel']}_nbe{uoi_params['n_boots_est']}_nl{uoi_params['n_lambdas']}_ss{uoi_params['stability_selection']}.pdf"
    )
    with PdfPages(viz_file) as pdf:
        plt.figure(figsize=(10, 8))
        nx.draw(G, pos, with_labels=False, node_size=50, edge_color=[data['color'] for _, _, data in G.edges(data=True)], width=[abs(data['weight']) for _, _, data in G.edges(data=True)], arrows=True)
        plt.title('Neuron Connections with Strengths')
        pdf.savefig()
        plt.close()

    elapsed_time = time.time() - start_time
    logging.info(f"Job completed successfully in {elapsed_time:.2f} seconds.")

except Exception as e:
    logging.error("Job failed.", exc_info=True)
    raise