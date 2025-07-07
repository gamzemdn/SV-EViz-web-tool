from dash import html, dcc
import dash_bootstrap_components as dbc

def get_survivor_layout():
    shared_input_style = {
        'fontFamily': '"Times New Roman", Times, serif',
        'width': '200px',
        'marginLeft': '10px'
    }

    shared_label_style = {
        'fontFamily': '"Times New Roman", Times, serif',
        'width': '250px',
        'display': 'inline-block'
    }

    shared_row_style = {
        'marginBottom': '8px'
    }

    return html.Div([
        html.H3("SURVIVOR", style={'padding': '1rem', 'fontSize': '24px', 'fontFamily': '"Times New Roman", Times, serif', 'fontWeight': 'bold'}),
        html.H3("Upload Caller VCF File", style={'marginTop': '0px', 'fontFamily': '"Times New Roman", Times, serif', 'fontSize': '20px'}),
        # Upload caller area
        dcc.Upload(
            id='survivor-upload',
            children=html.Div(['📎 Drag and Drop or Select File']),
            style={
                'width': '460px', 'height': '60px', 'lineHeight': '60px',
                'borderWidth': '1px', 'borderStyle': 'dashed', 'borderRadius': '5px',
                'textAlign': 'center', 'margin': '10px', 'fontFamily': '"Times New Roman", Times, serif'
            },
            multiple=False
        ),
        
        html.H3("Upload Reference VCF File", style={'marginTop': '0px', 'fontFamily': '"Times New Roman", Times, serif', 'fontSize': '20px'}),
        # Upload reference area
        dcc.Upload(
            id='survivor-upload-ref',
            children=html.Div(['📎 Drag and Drop or Select File']),
            style={
                'width': '460px', 'height': '60px', 'lineHeight': '60px',
                'borderWidth': '1px', 'borderStyle': 'dashed', 'borderRadius': '5px',
                'textAlign': 'center', 'margin': '10px', 'fontFamily': '"Times New Roman", Times, serif'
            },
            multiple=False
        ),
        
        html.Br(),

        # Parameter options
        html.Div([

            html.Div([
                html.Label('Use Univar SV Catalog [1:Yes 0:No] (default: No)', style=shared_label_style),
                dcc.RadioItems(
                    id='param_univar',
                    options=[{'label': 'Yes', 'value': '1'}, {'label': 'No', 'value': '0'}],
                    value='0',
                    inline=True,
                    labelStyle={'marginLeft': '10px', 'fontFamily': '"Times New Roman", Times, serif'}
                )
            ], style=shared_row_style),

            html.Div([
                html.Label('Maximum allowed distance (bp):', style=shared_label_style),
                dcc.Input(id='param_dist', type='number', min=0, placeholder='default: 1000', style=shared_input_style)
            ], style=shared_row_style),

            html.Div([
                html.Label('Number of callers supporting the call:', style=shared_label_style),
                dcc.Input(id='param_callers', type='number', min=0, placeholder='default: 2', style=shared_input_style)
            ], style=shared_row_style),

            html.Div([
                html.Label('Type matching [1:Yes 0:No] (default: Yes)', style=shared_label_style),
                dcc.RadioItems(
                    id='param_type',
                    options=[{'label': 'Yes', 'value': '1'}, {'label': 'No', 'value': '0'}],
                    value='1',
                    inline=True,
                    labelStyle={'marginLeft': '10px', 'fontFamily': '"Times New Roman", Times, serif'}
                )
            ], style=shared_row_style),

            html.Div([
                html.Label('Strand matching [1:Yes 0:No] (default: Yes)', style=shared_label_style),
                dcc.RadioItems(
                    id='param_strand',
                    options=[{'label': 'Yes', 'value': '1'}, {'label': 'No', 'value': '0'}],
                    value='1',
                    inline=True,
                    labelStyle={'marginLeft': '10px', 'fontFamily': '"Times New Roman", Times, serif'}
                )
            ], style=shared_row_style),

            html.Div([
                html.Label('Allow duplicate entry [1:Yes 0:No] (default: No)', style=shared_label_style),
                dcc.RadioItems(
                    id='param_dup',
                    options=[{'label': 'Yes', 'value': '1'}, {'label': 'No', 'value': '0'}],
                    value='0',
                    inline=True,
                    labelStyle={'marginLeft': '10px', 'fontFamily': '"Times New Roman", Times, serif'}
                )
            ], style=shared_row_style),

            html.Div([
                html.Label('Minimum SV size (bp):', style=shared_label_style),
                dcc.Input(id='param_sv_size', type='number', min=0, placeholder='default: 30', style=shared_input_style)
            ], style=shared_row_style),

            dbc.Button('Merge Files', id='merge-button', n_clicks=0, color="primary", className="mt-2", style={'marginTop': '15px', 'fontFamily': '"Times New Roman", Times, serif'}),
        ], style={'margin': '10px'}),
        
        dcc.Store(id='survivor-uploaded-paths', data=[]),
        # Output area
        html.Div(id='survivor-output', style={'marginTop': '20px'})
    ])
