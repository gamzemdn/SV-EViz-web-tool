from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc

def get_visualize_layout():
    return dbc.Container([
        html.Div([
            html.Div(id="visual_output"),
            dcc.Store(id='circos-data-store'),
            dcc.Store(id='svtype-colors-store'),          
            html.Div(id="circos-container"),
            html.Div(id="default-circos-output"),
            dcc.Store(id='manhattan-store'),
            html.Div(id="manhattan-controls-container", style={"display": "none"}, children=[
                html.Label("Manhattan Threshold:", style={'marginTop': '10px'}),
                dcc.Slider(
                    id='manhattan-threshold-slider',
                    min=0, max=10, step=0.1, value=6,
                    marks={i: str(i) for i in range(0, 11)},
                    tooltip={"placement": "bottom", "always_visible": False}
                ),
                html.Br(),
                html.Label("Select SVTYPEs:", style={'marginTop': '10px'}),
                dcc.Dropdown(
                    id='manhattan-svtype-selector',
                    options=[],  # to be filled later
                    value=[],
                    multi=True,
                    placeholder="Select one or more SV types",
                    style={"width": "300px"}
                )
            ]),

            html.Div(id="default-dashbio-manhattanplot")  # ➕ Display plot here
        ])
    ], fluid=True)
