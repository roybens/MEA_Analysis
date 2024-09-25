from dash import dcc, html, callback_context
from dash.dependencies import Input, Output, State, ALL
from data_processing import parse_contents, parse_mat
from plot_functions import plot_activity_graphs, plot_isi_graph, plot_bar_with_p_values
import io, zipfile, base64

default_colors = ['blue', 'green', 'red', 'purple', 'orange', 'brown', 'pink', 'gray', 'olive', 'cyan']

def register_callbacks(app):
    @app.callback(
        Output('tabs-content-example', 'children'),
        Input('tabs-example', 'value')
    )
    def render_content(tab):
        if tab == 'tab-upload':
            return html.Div([
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
                dcc.Loading(
                    id='loading-output',
                    type='circle',
                    children=html.Div(id='output-data-upload')
                ),
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
                dcc.Loading(
                    id='submit-button-loading',
                    type='circle',
                    children=html.Button('Submit', id='submit-button', n_clicks=0)
                ),
                dcc.Store(id='hidden-hex-codes'),  # Store for updated hex codes
                dcc.Store(id='neuron-order-store'),  # Store for neuron order
                dcc.Loading(
                    id='downloading-graphs-svg',
                    type='circle',
                    children=html.Button("Download as SVG", id="download-svg-button")
                ),
                dcc.Loading(
                    id='downloading-graphs-png',
                    type='circle',
                    children=html.Button("Download as PNG", id="download-png-button")
                ),
                dcc.Download(id="download-svg"),
                dcc.Download(id="download-png"),
                html.Div(id='color-dropdown-container'),
                dcc.Loading(
                    id='loading-graphs',
                    type='circle',
                    children=html.Div(id='graphs-container')  # Wrap graphs-container with dcc.Loading
                )  # Graph container
            ])
        elif tab == 'tab-isi':
            return html.Div([
                dcc.Upload(
                    id='upload-mat-file',
                    children=html.Div([
                        'Drag and Drop or ',
                        html.A('Select a MATLAB File')
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
                ),
                html.Div(id='mat-upload-output'),
                dcc.Checklist(
                    id='div-checklist-isi',
                    options=[],  # Options will be populated based on the uploaded data
                    value=[],
                    inline=True
                ),
                html.Div(id='color-inputs-container-isi'),  # Color inputs for ISI
                dcc.Dropdown(
                    id='metric-dropdown',
                    options=[],  # No default options
                    value=None,
                    clearable=False
                ),
                html.Button('Submit', id='submit-isi', n_clicks=0),
                html.Button("Download All as SVG", id="download-all-isi-svg-button"),
                html.Button("Download All as PNG", id="download-all-isi-png-button"),
                dcc.Store(id='hidden-hex-codes-isi'),  # Store for updated hex codes for ISI
                dcc.Loading(id='loading-isi', type='circle', children=html.Div(id='isi-graph-container')),

                dcc.Download(id="download-isi-svg"),
                dcc.Download(id="download-isi-png")
            ])
        elif tab == 'tab-lda':
            return html.Div([
                dcc.Upload(
                    id='upload-csv-lda',
                    children=html.Div([
                        'Drag and Drop or ',
                        html.A('Select CSV File for LDA')
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
                ),
                html.Div(id='output-csv-upload'),
                dcc.Checklist(
                    id='div-checklist-lda',
                    options=[],  # Options will be populated based on the selected dataset
                    value=[],
                    inline=True
                ),
                html.Button('Submit', id='submit-lda', n_clicks=0),
                dcc.Loading(id='loading-lda', type='circle', children=html.Div(id='lda-graph-container')),
                html.Button("Download as SVG", id="download-lda-svg-button"),
                html.Button("Download as PNG", id="download-lda-png-button"),
                dcc.Download(id="download-lda-svg"),
                dcc.Download(id="download-lda-png")
            ])

    @app.callback(
        [Output('div-checklist', 'options'),
         Output('div-checklist', 'value'),
         Output('chip-well-checklist', 'options'),
         Output('chip-well-checklist', 'value'),
         Output('color-inputs-container', 'children'),
         Output('neuron-options-container', 'children'),
         Output('graphs-container', 'children', allow_duplicate=True)],
        [Input('upload-data', 'contents')],
        [State('upload-data', 'filename')],
        prevent_initial_call=True  # Add this to prevent initial call
    )
    def update_options(contents, filename):
        if contents is None:
            return [], [], [], [], [], html.Div([])

        df = parse_contents(contents, filename)
        if isinstance(df, html.Div):
            return [], [], [], [], [], df

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

        return div_options, unique_divs, chip_well_options, unique_chip_wells, color_inputs, neuron_checklist_items, html.Div([])

    @app.callback(
        [Output('graphs-container', 'children', allow_duplicate=True),
         Output('download-svg', 'data', allow_duplicate=True),
         Output('download-png', 'data', allow_duplicate=True)],
        [Input('submit-button', 'n_clicks'),
         Input('download-svg-button', 'n_clicks'),
         Input('download-png-button', 'n_clicks')],
        [State('upload-data', 'contents'),
         State('upload-data', 'filename'),
         State('div-checklist', 'value'),
         State('chip-well-checklist', 'value'),
         State({'type': 'color-hex', 'index': ALL}, 'value'),
         State({'type': 'color-hex', 'index': ALL}, 'id'),
         State('neuron-order-store', 'data'),
         State({'type': 'neuron-checklist', 'index': ALL}, 'value')],
        prevent_initial_call=True  # Add this to prevent initial call
    )
    def update_graph(n_clicks, n_clicks_svg, n_clicks_png, contents, filename, selected_divs, selected_chip_wells, hex_codes, hex_codes_ids, neuron_order, selected_neuron_types_lists):

        triggered_id = callback_context.triggered[0]['prop_id'].split('.')[0]

        if contents is None or n_clicks == 0:
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
        # Build the color mapping dictionary using IDs
        selected_colors_dict = {hex_code_id['index']: hex_code for hex_code_id, hex_code in zip(hex_codes_ids, hex_codes)}
        hex_codes = hex_codes or {}
        selected_colors_dict = {genotype: selected_colors_dict[genotype] for genotype in selected_neuron_types if genotype in selected_colors_dict}
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
            elif 'networks' in filename.lower():
                images, svg_bytes_list, png_bytes_list = plot_bar_with_p_values(df, selected_divs, [metric], ordered_genotypes, selected_colors_dict)
            else:
                raise(TypeError("No file of name activity or network found"))
            graphs.extend(images)
            for i, svg_data in enumerate(svg_bytes_list):
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

        if triggered_id == 'download-svg-button' and all_svg_bytes_list:
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, 'w') as zf:
                for filename, svg_data in all_svg_bytes_list:
                    zf.writestr(filename, svg_data)
            zip_buffer.seek(0)
            zip_base64 = base64.b64encode(zip_buffer.read()).decode('utf-8')
            return table, dict(content=zip_base64, base64=True, filename="plots.zip", type="application/zip"), None

        if triggered_id == 'download-png-button' and all_png_bytes_list:
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, 'w') as zf:
                for filename, png_data in all_png_bytes_list:
                    zf.writestr(filename, png_data)
            zip_buffer.seek(0)
            zip_base64 = base64.b64encode(zip_buffer.read()).decode('utf-8')
            return table, None, dict(content=zip_base64, base64=True, filename="plots.zip", type="application/zip")

        return table, None, None

    @app.callback(
        [Output('div-checklist-isi', 'options'),
         Output('div-checklist-isi', 'value'),
         Output('color-inputs-container-isi', 'children'),
         Output('metric-dropdown', 'options'),
         Output('isi-graph-container', 'children', allow_duplicate=True)],
        Input('upload-mat-file', 'contents'),
        State('upload-mat-file', 'filename'),
        prevent_initial_call=True  # Add this to prevent initial call
    )
    def update_div_options_isi(contents, filename):
        if contents is None:
            return [], [], [], [], []
        
        data_df, available_metrics = parse_mat(contents)
        
        divs = sorted(data_df['DIV'].unique())
        div_options = [{'label': f'DIV {div}', 'value': div} for div in divs]

        default_color_values = {}
        unique_genotypes = data_df['NeuronType'].str.strip().unique()
        for i, neuron in enumerate(unique_genotypes):
            default_color_values[neuron] = default_colors[i % len(default_colors)]

        color_inputs = []
        for i, (neuron, default_color) in enumerate(default_color_values.items()):
            color_inputs.append(html.Label(f'Color for {neuron}'))
            color_inputs.append(dcc.Input(
                id={'type': 'color-hex-isi', 'index': neuron},
                type='text',
                placeholder='Enter hex code',
                value=default_color  # Set the initial value to the default color
            ))

        metric_options = [{'label': metric.replace('APFreqBins', ' AP Frequency'), 'value': metric} for metric in available_metrics]
        
        return div_options, [], color_inputs, metric_options, []  # Do not select any DIV initially

    @app.callback(
        [Output('isi-graph-container', 'children'),
         Output('download-isi-svg', 'data'),
         Output('download-isi-png', 'data')],
        [Input('submit-isi', 'n_clicks'),
         Input('download-all-isi-svg-button', 'n_clicks'),
         Input('download-all-isi-png-button', 'n_clicks')],
        [State('upload-mat-file', 'contents'),
         State('upload-mat-file', 'filename'),
         State('div-checklist-isi', 'value'),
         State('metric-dropdown', 'value'),
         State({'type': 'color-hex-isi', 'index': ALL}, 'value')],
        prevent_initial_call=True  # Add this to prevent initial call
    )
    def plot_isi_graph_callback(n_clicks, n_clicks_svg, n_clicks_png, contents, filename, selected_divs, selected_metric, hex_codes):
        triggered_id = callback_context.triggered[0]['prop_id'].split('.')[0]

        if contents is None or not selected_divs or not selected_metric:
            return html.Div('Please upload a file and select DIVs and a metric.'), None, None

        data_df, available_metrics = parse_mat(contents)
        hex_codes = hex_codes or {}
        unique_genotypes = data_df['NeuronType'].str.strip().unique()
        selected_colors_dict = {genotype: hex_codes[i] for i, genotype in enumerate(unique_genotypes)}

        img_element, svg_bytes, png_bytes, all_svg_bytes = plot_isi_graph(data_df, selected_divs, selected_metric, selected_colors_dict)
        svg_bytes.seek(0)
        
        if triggered_id == 'download-all-isi-svg-button' and all_svg_bytes:
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, 'w') as zf:
                zf.writestr(f"{selected_metric}_plot.svg", svg_bytes.read())
            zip_buffer.seek(0)
            zip_base64 = base64.b64encode(zip_buffer.read()).decode('utf-8')
            return None, dict(content=zip_base64, base64=True, filename="isi_plots.zip", type="application/zip"), None

        if triggered_id == 'download-all-isi-png-button':
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, 'w') as zf:
                zf.writestr(f"{selected_metric}_plot.png", png_bytes)
            zip_buffer.seek(0)
            zip_base64 = base64.b64encode(zip_buffer.read()).decode('utf-8')
            return None, None, dict(content=zip_base64, base64=True, filename="isi_plots.zip", type="application/zip")

        return img_element, None, None

    @app.callback(
        [Output('div-checklist-lda', 'options'),
         Output('output-csv-upload', 'children'),
         Output('lda-graph-container', 'children', allow_duplicate=True)],
        Input('upload-csv-lda', 'contents'),
        State('upload-csv-lda', 'filename'),
        prevent_initial_call=True  # Add this to prevent initial call
    )
    def update_lda_options(contents, filename):
        if contents is None:
            return [], 'No file uploaded yet.', []

        df = parse_contents(contents, filename)
        if isinstance(df, html.Div):
            return [], df, []

        unique_divs = sorted(df['DIV'].unique())
        div_options = [{'label': str(div), 'value': div} for div in unique_divs]
        return div_options, f'File {filename} successfully uploaded.', []

    @app.callback(
        [Output('lda-graph-container', 'children'),
         Output('download-lda-svg', 'data'),
         Output('download-lda-png', 'data')],
        [Input('submit-lda', 'n_clicks'),
         Input('download-lda-svg-button', 'n_clicks'),
         Input('download-lda-png-button', 'n_clicks')],
        [State('upload-csv-lda', 'contents'),
         State('upload-csv-lda', 'filename'),
         State('div-checklist-lda', 'value')],
        prevent_initial_call=True  # Add this to prevent initial call
    )
    def update_lda_graph(n_clicks, n_clicks_svg, n_clicks_png, contents, filename, selected_divs):
        triggered_id = callback_context.triggered[0]['prop_id'].split('.')[0]

        if contents is None or not selected_divs:
            return html.Div('Please upload a CSV file and select DIVs.'), None, None

        df = parse_contents(contents, filename)
        if isinstance(df, html.Div):
            return df, None, None

        filtered_df = df[df['DIV'].isin(selected_divs)]
        columns_to_remove = ['cov_BurstDuration', 'mean_Burst_Peak_Abs', 'cov_Burst_Peak_Abs',
                             'mean_Burst_Peak', 'cov_Burst_Peak', 'Time', 'Chip_ID', 'Well', 'Run_ID', 'DIV']
        columns_to_remove = [col for col in columns_to_remove if col in filtered_df.columns]
        filtered_df = filtered_df.drop(columns=columns_to_remove)
        filtered_df = filtered_df.dropna()

        X = filtered_df.drop(columns=['NeuronType'])
        y = filtered_df['NeuronType']
        
        # Get the unique neuron types
        neuron_types = filtered_df['NeuronType'].unique()
        palette = {neuron_type: f'C{i}' for i, neuron_type in enumerate(neuron_types)}
        lda = LinearDiscriminantAnalysis()
        X_lda = lda.fit_transform(X, y)

        column_names = [f'LD{i+1}' for i in range(X_lda.shape[1])]
        lda_df = pd.DataFrame(X_lda, columns=column_names)
        lda_df['NeuronType'] = y.values

        fig, ax = plt.subplots()
        sns.scatterplot(x=lda_df.columns[0], y=lda_df.columns[1] if len(lda_df.columns) > 1 else lda_df.columns[0], hue='NeuronType', data=lda_df, ax=ax)
        plt.xlabel('Linear Discriminant 1')
        plt.ylabel('Linear Discriminant 2' if len(lda_df.columns) > 1 else 'Linear Discriminant 1')
        plt.legend(loc='best')
        
        svg_bytes = io.BytesIO()
        plt.savefig(svg_bytes, format='svg')
        svg_bytes.seek(0)
        svg_data = svg_bytes.getvalue()

        png_bytes = io.BytesIO()
        plt.savefig(png_bytes, format='png')
        png_bytes.seek(0)
        png_data = png_bytes.getvalue()

        img_bytes = io.BytesIO()
        plt.savefig(img_bytes, format='png')
        img_bytes.seek(0)
        img_data = base64.b64encode(img_bytes.read()).decode('utf-8')
        img_element = html.Img(src='data:image/png;base64,{}'.format(img_data))

        plt.close(fig)

        if triggered_id == 'download-lda-svg-button':
            return img_element, dict(content=svg_data, base64=True, filename="lda_plot.svg"), None

        if triggered_id == 'download-lda-png-button':
            return img_element, None, dict(content=png_data, base64=True, filename="lda_plot.png")

        return img_element, None, None