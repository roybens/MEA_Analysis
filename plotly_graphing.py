import pandas as pd
import numpy as np
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import plotly.express as px
import plotly.graph_objects as go
import scipy.io
import csv
import sys
import socket
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from scipy.stats import ttest_ind, t
import matplotlib.pyplot as plt
import io
import base64

# Set Matplotlib backend to 'Agg' for non-GUI rendering
plt.switch_backend('Agg')

# Increase the maximum field size allowed
csv.field_size_limit(sys.maxsize)

# Initialize the Dash app
app = dash.Dash(__name__)

# App layout
app.layout = html.Div([
    dcc.Dropdown(
        id='dataset-dropdown',
        options=[
            {'label': 'Compiled Networks 4', 'value': 'Compiled_Networks_4.csv'},
            {'label': 'Extended Metrics MAT', 'value': 'extendedMetrics.mat'},
            {'label': 'Compiled Networks 5', 'value': 'Compiled_Networks_5.csv'},
            {'label': 'Compiled Networks SHANK3', 'value': 'Compiled_Networks_shank3_1.csv'}
        ],
        value='Compiled_Networks_4.csv'  # Default dataset
    ),
    dcc.Checklist(
        id='div-checklist',
        options=[],  # Options will be populated based on the selected dataset
        value=[],
        inline=True
    ),
    dcc.Dropdown(
        id='plot-type-dropdown',
        options=[
            {'label': 'Bar Graph', 'value': 'bar'},
            {'label': 'Histogram', 'value': 'hist'},
            {'label': 'Line Plot', 'value': 'line'},
            {'label': 'LDA', 'value': 'lda'}
        ],
        value='bar'
    ),
    dcc.Dropdown(
        id='metric-dropdown',
        options=[],  # Options will be populated based on the selected dataset
        value=None
    ),
    html.Div(id='graphs-container')
])

# Callback to update DIV checklist and metric dropdown based on the selected dataset
@app.callback(
    [Output('div-checklist', 'options'),
     Output('div-checklist', 'value'),
     Output('metric-dropdown', 'options'),
     Output('metric-dropdown', 'value')],
    [Input('dataset-dropdown', 'value')]
)
def update_div_and_metric_options(selected_dataset):
    if selected_dataset == 'Compiled_Networks_4.csv':
        data = pd.read_csv(selected_dataset)
        unique_divs = sorted(data['DIV'].unique())
        metrics = ['IBI', 'Burst_Peak', 'Number_Bursts', 'Spike_per_Burst', 'BurstDuration']
    elif selected_dataset == 'Compiled_Networks_5.csv':
        data = pd.read_csv(selected_dataset)
        unique_divs = sorted(data['Run_ID'].unique())
        # Exclude 'Chip_ID' from the metrics
        metrics = data.columns.difference(['Run_ID', 'NeuronType', 'Chip_ID']).tolist() + ['LDA']
    elif selected_dataset == 'Compiled_Networks_shank3_1.csv':
        data = pd.read_csv(selected_dataset)
        unique_divs = sorted(data['DIV'].unique())
        exclude_columns = ['Run_ID', 'Chip_ID', 'Time', 'NeuronType', 'DIV', 'Well']
        metrics = [col for col in data.columns if col not in exclude_columns]
    else:  # extendedMetrics.mat
        mat_data = scipy.io.loadmat(selected_dataset)
        allExtMetrics = mat_data['allExtMetrics']
        divs = []
        for i in range(allExtMetrics.shape[1]):
            for j in range(allExtMetrics[0, i].shape[0]):
                metric = allExtMetrics[0, i][j, 0]
                divs.append(metric['DIV'][0][0])
        unique_divs = sorted(set(int(div) for div in divs))  # Convert to integers for hashing
        metrics = ['networkAPFreqBins', 'burstAPFreqBins', 'nonburstAPFreqBins']

    div_options = [{'label': str(div), 'value': div} for div in unique_divs]
    return div_options, unique_divs, [{'label': metric, 'value': metric} for metric in metrics], metrics[0]


def perform_lda(file_path, run_ids):
    # Load the dataset
    df = pd.read_csv(file_path)

    # Filter the DataFrame for specified Run_IDs
    filtered_df = df[df['Run_ID'].isin(run_ids)]

    # Columns to consider for LDA
    columns_to_consider = ['NeuronType', 'mean_IBI', 'cov_IBI', 'mean_Burst_Peak', 'cov_Burst_Peak',
                           'mean_Burst_Peak_Abs', 'cov_Burst_Peak_Abs', 'Number_Bursts', 'mean_Spike_per_Burst',
                           'cov_Spike_per_Burst', 'mean_BurstDuration', 'MeanNetworkISI', 'CoVNetworkISI',
                           'MeanWithinBurstISI', 'CoVWithinBurstISI', 'MeanOutsideBurstISI', 'CoVOutsideBurstISI',
                           'Fanofactor']

    # Ensure the columns to consider exist in the filtered dataframe
    columns_to_consider = [col for col in columns_to_consider if col in filtered_df.columns]

    filtered_df = filtered_df[columns_to_consider].dropna()

    # Separate features and target variable
    X = filtered_df.drop(columns=['NeuronType'])
    y = filtered_df['NeuronType']

    # Apply LDA
    lda = LinearDiscriminantAnalysis()
    try:
        X_lda = lda.fit_transform(X, y)
    except Exception as e:
        return go.Figure()

    # Create a DataFrame for LDA results
    lda_df = pd.DataFrame(X_lda, columns=['LD1', 'LD2'])
    lda_df['NeuronType'] = y.values

    # Plot the LDA results
    fig = px.scatter(lda_df, x='LD1', y='LD2', color='NeuronType', title='LDA of NeuronType for selected Run_IDs')
    fig.update_layout(plot_bgcolor='white', paper_bgcolor='white')
    return fig


def plot_bar_with_p_values(data, divs, metrics):
    unique_genotypes = data['NeuronType'].unique()
    colors = ['blue', 'red']
    markers = ['o', 's']

    images = []
    for output_type in metrics:
        total_genotypes = len(unique_genotypes)
        print(f"Number of unique Genotypes: {total_genotypes}")

        output_arrays = {genotype: [] for genotype in unique_genotypes}
        chip_arrays = {genotype: [] for genotype in unique_genotypes}
        well_arrays = {genotype: [] for genotype in unique_genotypes}
        
        for i in divs:
            for genotype in unique_genotypes:
                temp_df = data.loc[(data['DIV'] == i) & (data['NeuronType'].str.strip() == genotype)]
                output_arrays[genotype].append(np.array(temp_df[output_type]))
                chip_arrays[genotype].append(np.array(temp_df['Chip_ID']))
                well_arrays[genotype].append(np.array(temp_df['Well']))
        
        bar_width = 0.25  # Slimmer bar width
        gap_between_bars = 0.0  # Small gap between bars within a group
        gap_between_groups = bar_width  # Gap between groups of bars

        total_bar_group_width = total_genotypes * bar_width + (total_genotypes - 1) * gap_between_bars
        total_plot_width = len(divs) * (total_bar_group_width + gap_between_groups) + gap_between_groups

        x_genotype = {genotype: [] for genotype in unique_genotypes}
        base_x_coordinate = np.arange(len(divs)) * (total_bar_group_width + gap_between_groups) + gap_between_groups + bar_width
        offset = (total_genotypes * bar_width + (total_genotypes - 1) * gap_between_bars) / 2
        centered_x = base_x_coordinate - offset + bar_width / 2
        for i, genotype in enumerate(unique_genotypes):
            x_genotype[genotype] = centered_x + (i) * (bar_width + gap_between_bars)

        fig, ax = plt.subplots()
        mean_data_all = {}
        yerr_data_all = {}
        n_data_all = {}

        for i, (genotype, color, marker) in enumerate(zip(unique_genotypes, colors, markers)):
            y_data = output_arrays[genotype]
            chipy_data = chip_arrays[genotype]
            welly_data = well_arrays[genotype]
            mean_data = [np.mean([n for n in yi if np.isfinite(n)]) for yi in y_data]
            yerr_data = [np.std([n for n in yi if np.isfinite(n)], ddof=1) / np.sqrt(np.size(yi)) for yi in y_data]
            n_data = [len(yi) for yi in y_data]
            mean_data_all[genotype] = mean_data
            yerr_data_all[genotype] = yerr_data
            n_data_all[genotype] = n_data
            
            alpha_value = 0.5
            ax.bar(x_genotype[genotype], mean_data, yerr=yerr_data, capsize=3, width=bar_width, color=color, label=genotype, alpha=0.7)

            jitter_amount = 0.07
            for j in range(len(x_genotype[genotype])):
                combined_data = [str(chip) + str(well) for chip, well in zip(chipy_data[j], welly_data[j])]
                for k in range(len(y_data[j])):
                    ax.scatter(
                        x_genotype[genotype][j] + np.random.uniform(-jitter_amount, jitter_amount, 1),
                        y_data[j][k],
                        s=20,
                        color=color,
                        marker=marker
                    )

        for i in range(len(base_x_coordinate)):
            maxim = max(max(array) for genotype_arrays in output_arrays.values() for array in genotype_arrays)
            count = 1
            p_values = []
            for j, genotype1 in enumerate(unique_genotypes):
                for k, genotype2 in enumerate(unique_genotypes):
                    if j < k:
                        mean1, sem1, n1 = mean_data_all[genotype1][i], yerr_data_all[genotype1][i], n_data_all[genotype1][i]
                        mean2, sem2, n2 = mean_data_all[genotype2][i], yerr_data_all[genotype2][i], n_data_all[genotype2][i]
                        sed = np.sqrt(sem1**2.0 + sem2**2.0)
                        t_stat = (mean1 - mean2) / sed
                        degreef = n1 + n2 - 2
                        p_value = (1.0 - t.cdf(abs(t_stat), degreef)) * 2.0
                        p_values.append([mean1, sem1, mean2, sem2, p_value])

                        x1, x2 = x_genotype[genotype1][i], x_genotype[genotype2][i]
                        sign = "***" if p_value <= 0.001 else "**" if p_value <= 0.01 else "*" if p_value <= 0.05 else "ns"
                        if sign != 'ns':
                            ax.plot([x1, x2], [maxim + 0.05 * maxim * (count)] * 2, 'k', linewidth=1.5)
                            ax.text((x1 + x2) / 2, maxim + 0.05 * maxim * (count), f"{sign}", ha='center', va='bottom', fontsize=7)
                        else:
                            ax.plot([x1, x2], [maxim + 0.05 * maxim * (count)] * 2, 'k', linewidth=1.5)
                            ax.text((x1 + x2) / 2, maxim + 0.05 * maxim * (count), f"{sign}", ha='center', va='bottom', fontsize=7)
                        count += 1
                        y_ticks = ax.get_yticks()
                        y_ticks_selected = [y_ticks[0], y_ticks[len(y_ticks) // 2], y_ticks[-1]]
                        ax.set_yticks(y_ticks_selected)
                        ax.spines['top'].set_visible(False)
                        ax.spines['right'].set_visible(False)
        plt.title(f"{output_type}", fontsize=14)
        plt.xlabel('DIV', fontsize=12)
        plt.ylabel(f"{output_type}", fontsize=12)
        plt.xticks(base_x_coordinate, divs, fontsize=10)
        plt.legend(title='NeuronType', loc='upper left', bbox_to_anchor=(1, 1), fontsize='small')
        plt.tight_layout()
        plt.xlim([0, total_plot_width])
        
        # Save the plot to a BytesIO object
        img_bytes = io.BytesIO()
        plt.savefig(img_bytes, format='png')
        img_bytes.seek(0)
        img_data = base64.b64encode(img_bytes.read()).decode('utf-8')
        img_element = html.Img(src='data:image/png;base64,{}'.format(img_data))
        images.append(img_element)
        plt.close(fig)
    return images


# Callback to update the graph based on the selected dataset, DIVs, plot type, and metric
@app.callback(
    Output('graphs-container', 'children'),
    [Input('dataset-dropdown', 'value'),
     Input('div-checklist', 'value'),
     Input('plot-type-dropdown', 'value'),
     Input('metric-dropdown', 'value')]
)
def update_graph(selected_dataset, selected_divs, plot_type, selected_metric):
    if selected_dataset == 'Compiled_Networks_4.csv':
        data = pd.read_csv(selected_dataset)
        filtered_data = data[data['DIV'].isin(selected_divs) & data['NeuronType'].str.contains('WT|HET', regex=True)]
        
        if plot_type == 'bar':
            return plot_bar_with_p_values(filtered_data, selected_divs, [selected_metric])
        else:
            fig = go.Figure()
            means = filtered_data.groupby(['DIV', 'NeuronType'])[selected_metric].mean().unstack()
            if plot_type == 'hist':
                for neuron_type in ['WT cortex', 'HET cortex']:
                    for div in selected_divs:
                        div_data = filtered_data[filtered_data['DIV'] == div][selected_metric]
                        fig.add_trace(go.Histogram(
                            x=div_data,
                            name=f'{neuron_type} DIV {div}',
                            opacity=0.75
                        ))
                fig.update_layout(barmode='overlay')

            elif plot_type == 'line':
                for neuron_type in ['WT cortex', 'HET cortex']:
                    fig.add_trace(go.Scatter(
                        x=selected_divs,
                        y=means.get(neuron_type, pd.Series(index=selected_divs, data=np.nan)).values,
                        name=neuron_type,
                        mode='lines+markers'
                    ))

            fig.update_layout(
                title=f'{selected_metric} by DIV',
                xaxis_title='DIV',
                yaxis_title=f'Mean {selected_metric}',
                xaxis={'type': 'category'},
                plot_bgcolor='white',
                paper_bgcolor='white'
            )
            return [dcc.Graph(figure=fig)]
    elif selected_dataset == 'Compiled_Networks_shank3_1.csv':
        data = pd.read_csv(selected_dataset)
        filtered_data = data[data['DIV'].isin(selected_divs) & data['NeuronType'].str.contains('WT|HET', regex=True)]

        if plot_type == 'bar':
            return plot_bar_with_p_values(filtered_data, selected_divs, [selected_metric])
        elif plot_type == 'hist':
                for neuron_type in ['WT', 'HET']:
                    for div in selected_divs:
                        div_data = filtered_data[filtered_data['DIV'] == div][selected_metric]
                        fig.add_trace(go.Histogram(
                            x=div_data,
                            name=f'{neuron_type} DIV {div}',
                            opacity=0.75
                        ))
                fig.update_layout(barmode='overlay')

        elif plot_type == 'line':
                for neuron_type in ['WT', 'HET']:
                    fig.add_trace(go.Scatter(
                        x=selected_divs,
                        y=means.get(neuron_type, pd.Series(index=selected_divs, data=np.nan)).values,
                        name=neuron_type,
                        mode='lines+markers'
                    ))

        fig.update_layout(
                title=f'{selected_metric} by DIV',
                xaxis_title='DIV',
                yaxis_title=f'Mean {selected_metric}',
                xaxis={'type': 'category'},
                plot_bgcolor='white',
                paper_bgcolor='white'
            )
        return [dcc.Graph(figure=fig)]
            
    elif selected_dataset == 'Compiled_Networks_5.csv':
        if plot_type == 'lda':
            fig = perform_lda(selected_dataset, selected_divs)
            return [dcc.Graph(figure=fig)]
        else:
            data = pd.read_csv(selected_dataset)
            filtered_data = data[data['Run_ID'].isin(selected_divs) & data['NeuronType'].str.contains('WT|HET|HOM', regex=True)]
            means = filtered_data.groupby(['Run_ID', 'NeuronType'])[selected_metric].mean().unstack()

            fig = go.Figure()
            if plot_type == 'bar':
                return plot_bar_with_p_values(filtered_data, selected_divs, [selected_metric])
            elif plot_type == 'hist':
                for neuron_type in ['WT', 'HET', 'HOM']:
                    for div in selected_divs:
                        div_data = filtered_data[filtered_data['Run_ID'] == div][selected_metric]
                        fig.add_trace(go.Histogram(
                            x=div_data,
                            name=f'{neuron_type} Run_ID {div}',
                            opacity=0.75
                        ))
                fig.update_layout(barmode='overlay')

            elif plot_type == 'line':
                for neuron_type in ['WT', 'HET', 'HOM']:
                    fig.add_trace(go.Scatter(
                        x=selected_divs,
                        y=means.get(neuron_type, pd.Series(index=selected_divs, data=np.nan)).values,
                        name=neuron_type,
                        mode='lines+markers'
                    ))

            fig.update_layout(
                title=f'{selected_metric} by Run_ID',
                xaxis_title='Run_ID',
                yaxis_title=f'Mean {selected_metric}',
                xaxis={'type': 'category'},
                plot_bgcolor='white',
                paper_bgcolor='white'
            )
            return [dcc.Graph(figure=fig)]

    else:  # extendedMetrics.mat
        mat_data = scipy.io.loadmat(selected_dataset)
        allExtMetrics = mat_data['allExtMetrics']
        neuron_types = []
        metric_bins = []
        metric_edges = []
        run_ids = []

        for i in range(allExtMetrics.shape[1]):
            for j in range(allExtMetrics[0, i].shape[0]):
                metric = allExtMetrics[0, i][j, 0]
                if metric['DIV'][0][0] in selected_divs:
                    run_id = metric['Run_ID'][0][0]
                    if run_id in [57, 59]:
                        neuron_type = metric['NeuronType'][0][0]
                        bins = metric[selected_metric][0][0]
                        edges = metric[selected_metric.replace('Bins', 'Edges')][0][0]

                        if bins.size > 0 and edges.size > 0:
                            neuron_types.append(neuron_type)
                            metric_bins.append(bins)
                            metric_edges.append(edges)
                            run_ids.append(run_id)

        neuron_types = [str(nt[0]) for nt in neuron_types]
        run_ids = [int(rid[0]) for rid in run_ids]
        metric_bins = [bins.flatten() for bins in metric_bins]
        metric_edges = [edges.flatten() for edges in metric_edges]

        def normalize_bins(bins):
            total_count = np.sum(bins)
            return bins / total_count if total_count != 0 else bins

        def interpolate_bins(x, y, common_x):
            return np.interp(common_x, x, y)

        fig = go.Figure()
        common_x = np.logspace(-1, 3, 100)
        colors = {'HOM': 'red', 'WT': 'blue', 'HET': 'green'}

        unique_neuron_types = set(neuron_types)
        for neuron_type in unique_neuron_types:
            indices = [i for i, x in enumerate(neuron_types) if x == neuron_type]

            all_bins = []
            for i, idx in enumerate(indices):
                x = metric_edges[idx][:-1]
                y = normalize_bins(metric_bins[idx])
                if len(x) == len(y):
                    y_interp = interpolate_bins(x, y, common_x)
                    all_bins.append(y_interp)

            if all_bins:
                all_bins = np.array(all_bins)
                mean_bins = np.mean(all_bins, axis=0)
                sem_bins = np.std(all_bins, axis=0) / np.sqrt(all_bins.shape[0])

                fig.add_trace(go.Scatter(
                    x=common_x,
                    y=mean_bins,
                    mode='lines+markers',
                    name=f'{neuron_type} Neurons',
                    line=dict(color=colors.get(neuron_type, 'black')),
                    error_y=dict(type='data', array=sem_bins, visible=True)
                ))

        fig.update_layout(
            title=f'Normalized Line Plot of {selected_metric} for Neuron Types',
            xaxis_title=f'{selected_metric}',
            yaxis_title='Normalized Count',
            xaxis=dict(type='log'),
            yaxis=dict(showgrid=False),
            plot_bgcolor='white',
            paper_bgcolor='white'
        )
        return [dcc.Graph(figure=fig)]

def find_free_port():
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(('', 0))
    port = s.getsockname()[1]
    s.close()
    return port


# Run the app on an available port
if __name__ == '__main__':
    free_port = find_free_port()
    app.run_server(debug=True, port=free_port)
