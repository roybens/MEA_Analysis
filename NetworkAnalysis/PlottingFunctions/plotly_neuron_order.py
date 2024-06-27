import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State, ALL
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import csv
import sys
import socket
from scipy.stats import ttest_ind, t
import matplotlib.pyplot as plt
import io
import base64
import zipfile
from dash import callback_context

# Set Matplotlib backend to 'Agg' for non-GUI rendering
plt.switch_backend('Agg')

# Increase the maximum field size allowed
csv.field_size_limit(sys.maxsize)

# Initialize the Dash app
app = dash.Dash(__name__, suppress_callback_exceptions=True)

# Define the list of colors for selection
default_colors = ['blue', 'green', 'red', 'purple', 'orange', 'brown', 'pink', 'gray', 'olive', 'cyan']

app.layout = html.Div([
    dcc.Loading(
        id='loading-upload',
        type='circle',
        children=dcc.Upload(
            id='upload-data',
            children=html.Div([
                'Drag and Drop or ',
                html.A('Select Files')
            ]),
            style={
                'width': '100%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px'
            },
            multiple=False
        )
    ),
    html.Div(id='output-data-upload'),
    html.Div(id='neuron-options-container'),
    dcc.Checklist(
        id='div-checklist',
        options=[],  # Options will be populated based on the selected dataset
        value=[],
        inline=True
    ),
    dcc.Checklist(
        id='chip-well-checklist',
        options=[],  # Options will be populated based on the selected dataset
        value=[],
        inline=True
    ),
    html.Div(id='color-inputs-container'),  # Color inputs
    html.Button('Submit Colors', id='submit-colors-button', n_clicks=0),
    html.Button('Submit Neuron Order', id='submit-neuron-order-button', n_clicks=0),  
    dcc.Store(id='hidden-hex-codes'),  # Store for updated hex codes
    dcc.Store(id='neuron-order-store'),  # Store for neuron order
    html.Button("Download as SVG", id="download-svg-button"),
    html.Button("Download as PNG", id="download-png-button"),
    dcc.Download(id="download-svg"),
    dcc.Download(id="download-png"),
    html.Div(id='color-dropdown-container'),
    dcc.Loading(
        id='loading-graphs',
        type='circle',
        children=html.Div(id='graphs-container')  # Wrap graphs-container with dcc.Loading
    )  # Graph container
])





def parse_contents(contents, filename):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if 'csv' in filename:
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
            df.columns = df.columns.str.strip()  # Strip whitespace from column names
            df['Chip_Well'] = df['Chip_ID'].astype(str) + '_' + df['Well'].astype(str)
        else:
            return html.Div([
                'Unsupported file format'
            ])
    except Exception as e:
        return html.Div([
            'There was an error processing this file.'
        ])

    return df

def plot_activity_graphs(data, div, output_types, ordered_genotypes, selected_colors):
    markers = ['o', 's', 'D', '^', 'v', 'x', '*', 'H', '8', 'p']
    images = []
    svg_bytes_list = []
    png_bytes_list = []
    for output_type in output_types:
        total_genotypes = len(ordered_genotypes)
        if total_genotypes == 0:
            return [html.Div('No data to display.')], [], []

        output_arrays = {genotype: [] for genotype in ordered_genotypes}
        chip_arrays = {genotype: [] for genotype in ordered_genotypes}
        well_arrays = {genotype: [] for genotype in ordered_genotypes}

        for i in div:
            for genotype in ordered_genotypes:
                temp_df = data.loc[(data['DIV'] == i) & (data['NeuronType'] == genotype)]
                if output_type not in temp_df.columns:
                    return [html.Div(f'Output type {output_type} not found in data.')], [], []
                output_arrays[genotype].append(np.array(temp_df[output_type]))
                chip_arrays[genotype].append(np.array(temp_df['Chip_ID']))
                well_arrays[genotype].append(np.array(temp_df['Well']))

        bar_width = 0.25
        gap_between_bars = 0

        total_bar_group_width = total_genotypes * bar_width + (total_genotypes - 1) * gap_between_bars

        x_genotype = {genotype: [] for genotype in ordered_genotypes}
        base_x_coordinate = np.arange(len(div))
        offset = (total_genotypes * bar_width + (total_genotypes - 1) * gap_between_bars) / 2
        centered_x = base_x_coordinate - offset + bar_width / 2
        for i, genotype in enumerate(ordered_genotypes):
            x_genotype[genotype] = centered_x + i * (bar_width + gap_between_bars)

        fig, ax = plt.subplots()

        mean_data_all = {}
        yerr_data_all = {}
        n_data_all = {}

        for i, genotype in enumerate(ordered_genotypes):
            y_data = output_arrays[genotype]
            chipy_data = chip_arrays[genotype]
            welly_data = well_arrays[genotype]

            mean_data = [np.nanmean(yi) if len(yi) > 0 else np.nan for yi in y_data]
            yerr_data = [np.nanstd(yi, ddof=1) / np.sqrt(len(yi)) if len(yi) > 1 else np.nan for yi in y_data]
            n_data = [len(yi) for yi in y_data]

            mean_data_all[genotype] = mean_data
            yerr_data_all[genotype] = yerr_data
            n_data_all[genotype] = n_data

            if len(x_genotype[genotype]) != len(mean_data):
                print(f"Length mismatch for genotype {genotype}: {len(x_genotype[genotype])} vs {len(mean_data)}")

            ax.bar(x_genotype[genotype], mean_data, yerr=yerr_data, capsize=3, width=bar_width, color=selected_colors.get(genotype, 'blue'), edgecolor='black', ecolor='black', label=genotype, alpha=0.6)

            jitter_amount = 0.05
            for j in range(len(x_genotype[genotype])):
                combined_data = [str(chip) + str(well) for chip, well in zip(chipy_data[j], welly_data[j])]
                for k in range(len(y_data[j])):
                    ax.scatter(x_genotype[genotype][j] + np.random.uniform(-jitter_amount, jitter_amount, 1), y_data[j][k], s=10, color=selected_colors.get(genotype, 'blue'), marker=markers[i % len(markers)])

        valid_arrays = [array for genotype_arrays in output_arrays.values() for array in genotype_arrays if len(array) > 0]
        if valid_arrays:
            maxim = max(max(array[np.isfinite(array)]) for array in valid_arrays)
        else:
            maxim = 0

        for i in range(len(base_x_coordinate)):
            count = 1
            for j, genotype1 in enumerate(ordered_genotypes):
                for k, genotype2 in enumerate(ordered_genotypes):
                    if j < k:
                        mean1, sem1, n1 = mean_data_all[genotype1][i], yerr_data_all[genotype1][i], n_data_all[genotype1][i]
                        mean2, sem2, n2 = mean_data_all[genotype2][i], yerr_data_all[genotype2][i], n_data_all[genotype2][i]
                        if np.isnan(mean1) or np.isnan(mean2):
                            continue  # Skip t-test if any mean is NaN
                        sed = np.sqrt(sem1 ** 2.0 + sem2 ** 2.0)
                        t_stat = (mean1 - mean2) / sed
                        degreef = n1 + n2 - 2
                        alpha = 0.05
                        cv = t.ppf(1.0 - alpha, degreef)
                        p_value = (1.0 - t.cdf(abs(t_stat), degreef)) * 2.0
                        sign = "***" if p_value <= 0.001 else "**" if p_value <= 0.01 else "*" if p_value <= 0.05 else "ns"
                        if sign != 'ns':
                            ax.plot([x_genotype[genotype1][i], x_genotype[genotype2][i]], [maxim + 0.05 * maxim * count] * 2, 'k', linewidth=1.3)
                            ax.text((x_genotype[genotype1][i] + x_genotype[genotype2][i]) / 2, maxim + 0.05 * maxim * count, sign, ha='center', va='bottom', fontsize=7)
                            ax.axvline(x_genotype[genotype1][i], color='black', linestyle=':', linewidth=0.5)
                        count += 1

        plt.title(f"{output_type}", fontsize=14)
        plt.xlabel('DIV', fontsize=12)
        plt.ylabel(f"{output_type}", fontsize=12)
        plt.xticks(base_x_coordinate, div, fontsize=10)
        plt.legend(title='NeuronType', loc='upper left', bbox_to_anchor=(1, 1), fontsize='small')
        plt.tight_layout()

        img_bytes = io.BytesIO()
        plt.savefig(img_bytes, format='png')
        img_bytes.seek(0)
        img_data = base64.b64encode(img_bytes.read()).decode('utf-8')
        img_element = html.Img(src='data:image/png;base64,{}'.format(img_data))
        images.append(img_element)
        png_bytes_list.append(img_bytes.getvalue())

        svg_bytes = io.BytesIO()
        plt.savefig(svg_bytes, format='svg')
        svg_bytes.seek(0)
        svg_data = svg_bytes.getvalue()
        svg_bytes_list.append(svg_data)

        plt.close(fig)

    return images, svg_bytes_list, png_bytes_list

def plot_bar_with_p_values(data, divs, metrics, ordered_genotypes, selected_colors):
    markers = ['o', 's', 'D', '^', 'v', 'x', '*', 'H', '8', 'p']

    images = []
    svg_bytes_list = []
    png_bytes_list = []

    for output_type in metrics:
        total_genotypes = len(ordered_genotypes)
        print(f"Number of unique Genotypes: {total_genotypes}")

        output_arrays = {genotype: [] for genotype in ordered_genotypes}
        chip_arrays = {genotype: [] for genotype in ordered_genotypes}
        well_arrays = {genotype: [] for genotype in ordered_genotypes}

        for i in divs:
            for genotype in ordered_genotypes:
                temp_df = data.loc[(data['DIV'] == i) & (data['NeuronType'].str.strip() == genotype)]
                output_arrays[genotype].append(np.array(temp_df[output_type]))
                chip_arrays[genotype].append(np.array(temp_df['Chip_ID']))
                well_arrays[genotype].append(np.array(temp_df['Well']))

        bar_width = 0.20
        gap_between_bars = 0.0
        gap_between_groups = bar_width

        total_bar_group_width = total_genotypes * bar_width + (total_genotypes - 1) * gap_between_bars
        total_plot_width = len(divs) * (total_bar_group_width + gap_between_groups) + gap_between_groups

        x_genotype = {genotype: [] for genotype in ordered_genotypes}
        base_x_coordinate = np.arange(len(divs)) * (total_bar_group_width + gap_between_groups) + gap_between_groups + bar_width
        offset = (total_genotypes * bar_width + (total_genotypes - 1) * gap_between_bars) / 2
        centered_x = base_x_coordinate - offset + bar_width / 2
        for i, genotype in enumerate(ordered_genotypes):
            x_genotype[genotype] = centered_x + (i) * (bar_width + gap_between_bars)

        fig, ax = plt.subplots()
        mean_data_all = {}
        yerr_data_all = {}
        n_data_all = {}

        for i, genotype in enumerate(ordered_genotypes):
            y_data = output_arrays[genotype]
            chipy_data = chip_arrays[genotype]
            welly_data = well_arrays[genotype]
            mean_data = [np.mean(yi[np.isfinite(yi)]) if len(yi[np.isfinite(yi)]) > 0 else np.nan for yi in y_data]
            yerr_data = [np.std(yi[np.isfinite(yi)], ddof=1) / np.sqrt(np.size(yi[np.isfinite(yi)])) if len(yi[np.isfinite(yi)]) > 1 else np.nan for yi in y_data]
            n_data = [len(yi[np.isfinite(yi)]) for yi in y_data]
            mean_data_all[genotype] = mean_data
            yerr_data_all[genotype] = yerr_data
            n_data_all[genotype] = n_data

            ax.bar(x_genotype[genotype], mean_data, yerr=yerr_data, capsize=3, width=bar_width, color=selected_colors.get(genotype, 'blue'), label=genotype, alpha=0.7)

            jitter_amount = 0.05
            for j in range(len(x_genotype[genotype])):
                combined_data = [str(chip) + str(well) for chip, well in zip(chipy_data[j], welly_data[j])]
                for k in range(len(y_data[j])):
                    ax.scatter(
                        x_genotype[genotype][j] + np.random.uniform(-jitter_amount, jitter_amount, 1),
                        y_data[j][k],
                        s=10,
                        color=selected_colors.get(genotype, 'blue'),
                        marker=markers[i % len(markers)]
                    )

        valid_arrays = [array for genotype_arrays in output_arrays.values() for array in genotype_arrays if len(array) > 0]
        if valid_arrays:
            maxim = max(max(array) for array in valid_arrays)
        else:
            maxim = 0

        for i in range(len(base_x_coordinate)):
            count = 1
            p_values = []
            for j, genotype1 in enumerate(ordered_genotypes):
                for k, genotype2 in enumerate(ordered_genotypes):
                    if j < k:
                        mean1, sem1, n1 = mean_data_all[genotype1][i], yerr_data_all[genotype1][i], n_data_all[genotype1][i]
                        mean2, sem2, n2 = mean_data_all[genotype2][i], yerr_data_all[genotype2][i], n_data_all[genotype2][i]
                        if n1 > 1 and n2 > 1:
                            sed = np.sqrt(sem1 ** 2.0 + sem2 ** 2.0)
                            t_stat = (mean1 - mean2) / sed
                            degreef = n1 + n2 - 2

                            p_value = (1.0 - t.cdf(abs(t_stat), degreef)) * 2.0
                            p_values.append([mean1, sem1, mean2, sem2, p_value])

                            x1, x2 = x_genotype[genotype1][i], x_genotype[genotype2][i]
                            sign = "***" if p_value <= 0.001 else "**" if p_value <= 0.01 else "*" if p_value <= 0.05 else "ns"
                            if sign != 'ns':
                                ax.plot([x1, x2], [maxim + 0.05 * maxim * (count)] * 2, 'k', linewidth=1.5)
                                ax.text((x1 + x2) / 2, maxim + 0.05 * maxim * (count), f"{sign}", ha='center', va='bottom', fontsize=7)
                                ax.axvline(x1, color='black', linestyle=':', linewidth=0.5)
                                ax.axvline(x2, color='black', linestyle=':', linewidth=0.5)
                            count += 1
            y_ticks = ax.get_yticks()
            y_ticks_selected = [y_ticks[0], y_ticks[len(y_ticks) // 2], y_ticks[-1]]
            ax.set_yticks(y_ticks_selected)

        plt.title(f"{output_type}", fontsize=14)
        plt.xlabel('DIV', fontsize=12)
        plt.ylabel(f"{output_type}", fontsize=12)
        plt.xticks(base_x_coordinate, divs, fontsize=10)
        plt.legend(title='NeuronType', loc='upper left', bbox_to_anchor=(1, 1), fontsize='small')
        plt.tight_layout()
        plt.xlim([0, total_plot_width])

        img_bytes = io.BytesIO()
        plt.savefig(img_bytes, format='png')
        img_bytes.seek(0)
        img_data = base64.b64encode(img_bytes.read()).decode('utf-8')
        img_element = html.Img(src='data:image/png;base64,{}'.format(img_data))
        images.append(img_element)
        png_bytes_list.append(img_bytes.getvalue())

        svg_bytes = io.BytesIO()
        plt.savefig(svg_bytes, format='svg')
        svg_bytes.seek(0)
        svg_data = svg_bytes.getvalue()
        svg_bytes_list.append(svg_data)

        plt.close(fig)

    return images, svg_bytes_list, png_bytes_list

@app.callback(
    [Output('div-checklist', 'options'),
     Output('div-checklist', 'value'),
     Output('chip-well-checklist', 'options'),
     Output('chip-well-checklist', 'value'),
     Output('color-inputs-container', 'children'),  # Updated to generate color inputs dynamically
     Output('neuron-options-container', 'children')],
    [Input('upload-data', 'contents')],
    [State('upload-data', 'filename')]
)
def update_options(contents, filename):
    if contents is None:
        return [], [], [], [], [], []

    df = parse_contents(contents, filename)
    if isinstance(df, html.Div):
        return [], [], [], [], [], []

    unique_divs = sorted(df['DIV'].unique())
    df['NeuronType'] = df['NeuronType'].str.strip()
    unique_genotypes = df['NeuronType'].unique()
    unique_chip_wells = df['Chip_Well'].unique()

    div_options = [{'label': str(div), 'value': div} for div in unique_divs]
    chip_well_options = [{'label': str(chip_well), 'value': chip_well} for chip_well in unique_chip_wells]

    default_color_values = {}
    for i, neuron in enumerate(unique_genotypes):
        default_color_values[neuron] = default_colors[i % len(default_colors)]

    color_inputs = []
    neuron_checklist_items = []
    for i, (neuron, default_color) in enumerate(default_color_values.items()):
        color_inputs.append(html.Label(f'Color for {neuron}'))
        color_inputs.append(dcc.Input(
            id={'type': 'color-hex', 'index': neuron},
            type='text',
            placeholder='Enter hex code',
            value=default_color  # Set the initial value to the default color
        ))
        neuron_checklist_items.append(html.Div([
            dcc.Checklist(
                id={'type': 'neuron-checklist', 'index': neuron},
                options=[{'label': neuron, 'value': neuron}],
                value=[neuron],
                inline=True
            ),
            dcc.Dropdown(
                id={'type': 'neuron-order', 'index': neuron},
                options=[{'label': str(idx+1), 'value': idx+1} for idx in range(len(unique_genotypes))],
                value=i+1,
                clearable=False
            )
        ], style={'display': 'flex', 'alignItems': 'center'}))

    return div_options, unique_divs, chip_well_options, unique_chip_wells, color_inputs, neuron_checklist_items




@app.callback(
    Output('hidden-hex-codes', 'data'),
    Input('submit-colors-button', 'n_clicks'),
    State({'type': 'color-hex', 'index': ALL}, 'value'),
    State('upload-data', 'contents'),
    State('upload-data', 'filename'),
    prevent_initial_call=True
)
def update_colors_on_submit(n_clicks, hex_codes, contents, filename):
    if contents is None:
        return {}
    df = parse_contents(contents, filename)
    if isinstance(df, html.Div):
        return {}
    unique_genotypes = df['NeuronType'].str.strip().unique()
    selected_colors_dict = {genotype: hex_codes[i] for i, genotype in enumerate(unique_genotypes)}
    return selected_colors_dict

@app.callback(
    Output('neuron-order-store', 'data'),
    Input('submit-neuron-order-button', 'n_clicks'),
    Input({'type': 'neuron-order', 'index': ALL}, 'value'),
    State('upload-data', 'contents'),
    State('upload-data', 'filename'),
    prevent_initial_call=True
)
def update_neuron_order(n_clicks, selected_orders, contents, filename):
    if contents is None:
        return []

    if selected_orders is None or len(selected_orders) == 0:
        return []

    if len(selected_orders) != len(set(selected_orders)):
        raise ValueError("Duplicate neuron orders found. Please assign unique orders to each neuron.")

    df = parse_contents(contents, filename)
    if isinstance(df, html.Div):
        return []

    unique_genotypes = df['NeuronType'].str.strip().unique()
    genotype_order_dict = {genotype: selected_orders[i] for i, genotype in enumerate(unique_genotypes)}
    sorted_genotypes = [genotype for genotype, order in sorted(genotype_order_dict.items(), key=lambda item: item[1])]
    print("sorted genotypes",sorted_genotypes)
    return sorted_genotypes


@app.callback(
    [Output('graphs-container', 'children'),
     Output('download-svg', 'data'),
     Output('download-png', 'data')],
    [Input('hidden-hex-codes', 'data'),
     Input('div-checklist', 'value'),
     Input('chip-well-checklist', 'value'),
     Input('submit-neuron-order-button', 'n_clicks'),
     Input('submit-colors-button', 'n_clicks'),
     Input('download-svg-button', 'n_clicks'),
     Input('download-png-button', 'n_clicks')],
    [State('upload-data', 'contents'),
     State('upload-data', 'filename'),
     State('neuron-order-store', 'data'),
     State({'type': 'neuron-checklist', 'index': ALL}, 'value')]
)
def update_graph(hex_codes, selected_divs, selected_chip_wells, n_clicks_order, n_clicks_colors, n_clicks_svg, n_clicks_png, contents, filename, neuron_order, selected_neuron_types_lists):
    triggered_id = callback_context.triggered[0]['prop_id'].split('.')[0]

    if contents is None:
        return html.Div(['No data uploaded yet.']), None, None

    df = parse_contents(contents, filename)
    if isinstance(df, html.Div):
        return df, None, None

    selected_divs = sorted(selected_divs)
    df['NeuronType'] = df['NeuronType'].str.strip()
    selected_neuron_types = [item for sublist in selected_neuron_types_lists for item in sublist]  # Flatten list of lists
    df = df[df['NeuronType'].isin(selected_neuron_types)]
    df = df[df['Chip_Well'].isin(selected_chip_wells)]
    unique_genotypes = df['NeuronType'].unique()

    hex_codes = hex_codes or {}
    selected_colors_dict = {genotype: hex_codes.get(genotype, default_colors[i % len(default_colors)]) for i, genotype in enumerate(unique_genotypes)}
    ordered_genotypes = neuron_order if neuron_order else unique_genotypes

    ordered_genotypes = [genotype for genotype in ordered_genotypes if genotype in selected_neuron_types]
    if not ordered_genotypes:
        ordered_genotypes = unique_genotypes

    graphs = []
    excluded_columns = ['DIV', 'NeuronType', 'Chip_ID', 'Well', 'Chip_Well', 'Run_ID', 'Time']
    selected_metrics = [col for col in df.columns if col not in excluded_columns]

    all_svg_bytes_list = []
    all_png_bytes_list = []


    for metric in selected_metrics:
        if 'activity' in filename.lower():
            images, svg_bytes_list, png_bytes_list = plot_activity_graphs(df, selected_divs, [metric], ordered_genotypes, selected_colors_dict)
        else:
            images, svg_bytes_list, png_bytes_list = plot_bar_with_p_values(df, selected_divs, [metric], ordered_genotypes, selected_colors_dict)
        
        graphs.extend(images)
        for i,svg_data in enumerate(svg_bytes_list):
            all_svg_bytes_list.append((f"{metric}_plot.svg", svg_data))
        for i, png_data in enumerate(png_bytes_list):
            all_png_bytes_list.append((f"{metric}_plot.png", png_data))

    table_rows = []
    row = []
    for i, graph in enumerate(graphs):
        row.append(html.Td(graph, style={'width': '33%', 'padding': '10px'}))
        if (i + 1) % 3 == 0:
            table_rows.append(html.Tr(row))
            row = []
    if row:
        table_rows.append(html.Tr(row))

    table = html.Table(table_rows, style={'width': '100%', 'borderCollapse': 'collapse'})

    if triggered_id == 'download-svg-button':
        if all_svg_bytes_list:
            # Create a ZIP file in memory
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, 'w') as zf:
                for filename, svg_data in all_svg_bytes_list:
                    zf.writestr(filename, svg_data)
            zip_buffer.seek(0)
            zip_base64 = base64.b64encode(zip_buffer.read()).decode('utf-8')
            return table, dict(content=zip_base64, base64=True, filename="plots.zip", type="application/zip"), None

            # svg_data = svg_bytes_list[0]
            # svg_base64 = base64.b64encode(svg_data).decode('utf-8')
            # return table, dict(content=svg_base64, base64=True, filename="plot.svg", type="image/svg+xml"), None

    if triggered_id == 'download-png-button':
        if all_png_bytes_list:
            # Create a ZIP file in memory
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, 'w') as zf:
                for filename, png_data in all_png_bytes_list:
                    zf.writestr(filename, png_data)
            zip_buffer.seek(0)
            zip_base64 = base64.b64encode(zip_buffer.read()).decode('utf-8')
            return table, None, dict(content=zip_base64, base64=True, filename="plots.zip", type="application/zip")
            # png_data = png_bytes_list[0]
            # png_base64 = base64.b64encode(png_data).decode('utf-8')
            # return table, None, dict(content=png_base64, base64=True, filename="plot.png", type="image/png")

    return table, None, None





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
