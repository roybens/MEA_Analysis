import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import networkx as nx
import numpy as np
import seaborn as sns
number = 12
# Load the .xlsx file containing neuron locations
df = pd.read_excel('./data/unitinfos.xlsx')  # Replace with your actual filename

# Load the connectivity matrix
connectivity_matrix = np.load(f'/pscratch/sd/m/mpatil1/MEA_Analysis/Connections/output/32956845/connectivity_matrix_nbs10_nbe10_nl48_ss1.0.npy')

# Extract neuron locations from the DataFrame
x_coords = df['LocX']
y_coords = df['LocY']

# Create a directed graph for connections using NetworkX
G = nx.DiGraph()  # Use a directed graph (since the connections are directional)

# Add nodes with positions to the graph
for i, (x, y) in enumerate(zip(x_coords, y_coords)):
    G.add_node(i, pos=(x, y))  # Add a node with the neuron's position

# Use the connectivity_matrix to add edges
num_neurons = len(df)  # Assuming each row represents a neuron
for i in range(num_neurons):
    for j in range(num_neurons):
        weight = connectivity_matrix[i, j]  # Extract the connection strength
        if weight > 0:  # Excitatory connection (positive weight)
            G.add_edge(j, i, weight=weight, color='blue')  # j influences i
        elif weight < 0:  # Inhibitory connection (negative weight)
            G.add_edge(j, i, weight=weight, color='red')  # j inhibits i

# Get the positions of the nodes (neurons) from the graph
pos = {i: (x, y) for i, (x, y) in enumerate(zip(x_coords, y_coords))}

# Extract edge attributes (colors and weights) for visualization
edges = G.edges(data=True)
colors = [edge[2]['color'] for edge in edges]  # Red for inhibitory, blue for excitatory
weights = [abs(edge[2]['weight']) for edge in edges]  # Thickness proportional to weight

# Initialize dictionaries to store counts for each neuron
input_excitatory = {i: 0 for i in range(num_neurons)}
output_excitatory = {i: 0 for i in range(num_neurons)}
input_inhibitory = {i: 0 for i in range(num_neurons)}
output_inhibitory = {i: 0 for i in range(num_neurons)}

# Count each neuron's incoming and outgoing excitatory and inhibitory connections
for i, j, data in G.edges(data=True):
    weight = data['weight']
    if weight > 0:  # Excitatory connection
        input_excitatory[j] += 1  # Count as incoming to j
        output_excitatory[i] += 1  # Count as outgoing from i
    elif weight < 0:  # Inhibitory connection
        input_inhibitory[j] += 1  # Count as incoming to j
        output_inhibitory[i] += 1  # Count as outgoing from i

# Create a DataFrame with the counts
df_counts = pd.DataFrame({
    "Neuron": list(range(num_neurons)),
    "Incoming Excitatory": [input_excitatory[i] for i in range(num_neurons)],
    "Outgoing Excitatory": [output_excitatory[i] for i in range(num_neurons)],
    "Incoming Inhibitory": [input_inhibitory[i] for i in range(num_neurons)],
    "Outgoing Inhibitory": [output_inhibitory[i] for i in range(num_neurons)]
})

# Save the connection counts to an Excel file
df_counts.to_excel(f"neuron_connection_counts{number}.xlsx", index=False)
print(f"Neuron connection counts saved to 'neuron_connection_counts{number}.xlsx'")

# Open a PDF file to save all plots
with PdfPages(f'neuron_connections_visualizations{number}.pdf') as pdf:

    # 1. Save the original network plot with different colors and weights
    plt.figure(figsize=(10, 8))
    nx.draw(G, pos, with_labels=False, node_size=50, edge_color=colors, width=weights, arrows=True)
    plt.title('Neuron Connections with Different Strengths')
    pdf.savefig()  # Save to PDF
    plt.close()

    # 2. Histogram of Connection Strengths
    connection_strengths = [data['weight'] for _, _, data in G.edges(data=True)]
    plt.figure(figsize=(8, 6))
    plt.hist(connection_strengths, bins=20, color='skyblue', edgecolor='black')
    plt.title('Distribution of Connection Strengths')
    plt.xlabel('Connection Strength (Weight)')
    plt.ylabel('Frequency')
    pdf.savefig()  # Save to PDF
    plt.close()

    # 3. Scatter Plot of Neuron Locations with Connection Strength as a Third Dimension
    plt.figure(figsize=(10, 8))
    for i, (x, y) in enumerate(zip(x_coords, y_coords)):
        for j in range(num_neurons):
            weight = connectivity_matrix[i, j]
            if weight != 0:  # Only plot connections that exist
                plt.plot([x, x_coords[j]], [y, y_coords[j]], 'k-', lw=abs(weight) / max(connection_strengths) * 5, alpha=0.5)
                plt.scatter([x, x_coords[j]], [y, y_coords[j]], s=abs(weight) * 10, c='blue' if weight > 0 else 'red')
    plt.xlabel('LocX')
    plt.ylabel('LocY')
    plt.title('Neuron Locations with Connection Strengths')
    pdf.savefig()  # Save to PDF
    plt.close()

    # 4. Box Plot of Connection Counts with Quartiles
    fig, ax = plt.subplots(2, 2, figsize=(12, 10))

    sns.boxplot(data=df_counts, y="Incoming Excitatory", ax=ax[0, 0])
    ax[0, 0].set_title("Incoming Excitatory Connections Distribution")

    sns.boxplot(data=df_counts, y="Incoming Inhibitory", ax=ax[0, 1])
    ax[0, 1].set_title("Incoming Inhibitory Connections Distribution")

    sns.boxplot(data=df_counts, y="Outgoing Excitatory", ax=ax[1, 0])
    ax[1, 0].set_title("Outgoing Excitatory Connections Distribution")

    sns.boxplot(data=df_counts, y="Outgoing Inhibitory", ax=ax[1, 1])
    ax[1, 1].set_title("Outgoing Inhibitory Connections Distribution")

    plt.tight_layout()
    pdf.savefig(fig)  # Save the box plot figure to PDF
    plt.close(fig)

print(f"All plots saved to 'neuron_connections_visualizations{number}.pdf'")