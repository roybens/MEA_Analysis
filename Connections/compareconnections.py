import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

# List of file paths and their respective parameter set names
file_paths = [
    ('/pscratch/sd/m/mpatil1/MEA_Analysis/Connections/neuron_connection_counts4.xlsx', 'Default'),
    ('/pscratch/sd/m/mpatil1/MEA_Analysis/Connections/neuron_connection_counts7.xlsx', 'N-boot_sel: 10, stability_selection = .7'),
    ('/pscratch/sd/m/mpatil1/MEA_Analysis/Connections/neuron_connection_counts10.xlsx', 'N-boot_sel: 24, stability_selection = .7'),
    ('/pscratch/sd/m/mpatil1/MEA_Analysis/Connections/neuron_connection_counts9.xlsx', 'N-boot_sel:10, stability_selection = 1')
]

# Initialize an empty list to store DataFrames
df_list = []

# Read each file and add a 'Parameter Set' column
for file_path, param_name in file_paths:
    df = pd.read_excel(file_path)
    df['Parameter Set'] = param_name
    df_list.append(df)

# Concatenate all DataFrames into a single DataFrame
df_combined = pd.concat(df_list, ignore_index=True)

# Open a PDF file to save the combined figure with subplots
with PdfPages('neuron_connections_comparison_dynamic.pdf') as pdf:
    # Create a figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))  # 2x2 grid for line plots
    
    # Plot Incoming Excitatory Connections
    sns.lineplot(data=df_combined, x='Neuron', y='Incoming Excitatory', hue='Parameter Set', ax=axes[0, 0])
    axes[0, 0].set_title('Incoming Excitatory Connections per Neuron')
    axes[0, 0].set_xlabel('Neuron')
    axes[0, 0].set_ylabel('Incoming Excitatory Connections')

    # Plot Incoming Inhibitory Connections
    sns.lineplot(data=df_combined, x='Neuron', y='Incoming Inhibitory', hue='Parameter Set', ax=axes[0, 1])
    axes[0, 1].set_title('Incoming Inhibitory Connections per Neuron')
    axes[0, 1].set_xlabel('Neuron')
    axes[0, 1].set_ylabel('Incoming Inhibitory Connections')

    # Plot Outgoing Excitatory Connections
    sns.lineplot(data=df_combined, x='Neuron', y='Outgoing Excitatory', hue='Parameter Set', ax=axes[1, 0])
    axes[1, 0].set_title('Outgoing Excitatory Connections per Neuron')
    axes[1, 0].set_xlabel('Neuron')
    axes[1, 0].set_ylabel('Outgoing Excitatory Connections')

    # Plot Outgoing Inhibitory Connections
    sns.lineplot(data=df_combined, x='Neuron', y='Outgoing Inhibitory', hue='Parameter Set', ax=axes[1, 1])
    axes[1, 1].set_title('Outgoing Inhibitory Connections per Neuron')
    axes[1, 1].set_xlabel('Neuron')
    axes[1, 1].set_ylabel('Outgoing Inhibitory Connections')

    # Adjust layout for the line plots
    plt.tight_layout()
    
    # Save the line plot figure as a single PDF page
    pdf.savefig(fig)
    plt.close(fig)

    # Create a new figure for the violin plot
    fig_violin, ax = plt.subplots(figsize=(10, 8))
    
    # Melt the DataFrame for violin plotting (long format needed)
    df_melted = df_combined.melt(id_vars=["Neuron", "Parameter Set"],
                                 value_vars=["Incoming Excitatory", "Incoming Inhibitory", 
                                             "Outgoing Excitatory", "Outgoing Inhibitory"],
                                 var_name="Connection Type", value_name="Connection Count")
    
    # Plot the violin plot
    sns.violinplot(data=df_melted, x="Connection Type", y="Connection Count", hue="Parameter Set", split=True, ax=ax)
    ax.set_title('Distribution of Connection Counts by Type and Parameter Set')
    ax.set_xlabel('Connection Type')
    ax.set_ylabel('Connection Count')
    
    # Save the violin plot figure as a PDF page
    pdf.savefig(fig_violin)
    plt.close(fig_violin)

print("Subplots with line plots and violin plot saved to 'neuron_connections_comparison_dynamic.pdf'")