from dash import dcc, html

def create_layout():
    return html.Div([
        dcc.Tabs(id="tabs-example", value='tab-upload', children=[
            dcc.Tab(label='Upload Data', value='tab-upload'),
            dcc.Tab(label='ISI', value='tab-isi'),
            dcc.Tab(label='LDA', value='tab-lda'),
        ]),
        html.Div(id='tabs-content-example')
    ])
