from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc

def get_visualize_layout():
    return dbc.Container([
        html.Div([
            dcc.Store(id='circos-data-store'),
            dcc.Store(id='svtype-colors-store'),
            dcc.Store(id='manhattan-store'),
            html.Div(id="visual_output"),
            html.Div(id="default-circos-output")
        ]),
        
        
        dcc.Store(id='selected-input-source', data='caller'),
       # html.Div(id='vcf-files-list', style={'marginTop': '20px'}),
        dcc.Store(id='uploaded-file-path', data=""),
        dcc.Store(id='filtered-data', data=""),

       # html.Div(id='visualize-status', style={'marginTop': '20px'}),
        dbc.Row([
            dbc.Col(html.Div(id='vcf-table-container', style={'marginTop': '30px'}), width=12)
        ]),

        dbc.Row([
            dbc.Col(html.Div(id='vcf-figures-container', style={'marginTop': '30px'}), width=12)
        ])

        
    ], fluid=True)
