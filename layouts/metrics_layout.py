from dash import html, dcc
import dash_bootstrap_components as dbc


def get_metrics_layout():
    return html.Div([
        html.Div(id='output-metrics'),  # For Truvari results
        html.Div(id='summary-visualization-output')  # For EvalSVCallers results
    ])

'''
def get_metricsss_layout():
    return html.Div([
        # --- SURVIVOR Upload Section ---
        html.Div(id='survivor-upload-section', children=[
            html.Label("Select Reference File Type:",style={'fontFamily': '"Times New Roman", Times, serif'}),
            dcc.Dropdown(
                id='reference-type',
                options=[
                    {'label': 'Standard VCF', 'value': 'normal'},
                    {'label': 'NA12878 EvalSVCallers', 'value': 'na12878'}
                ],
                value='normal',
                style={'width': '50%', 'fontFamily': '"Times New Roman", Times, serif'}
            ),  
            dcc.Upload(id='upload-caller', children=html.Button('Upload Caller File'), style={'marginTop': '20px', 'fontFamily': '"Times New Roman", Times, serif'},multiple=False),
            dcc.Upload(id='upload-reference', children=html.Button('Upload Reference File'), style={'fontFamily': '"Times New Roman", Times, serif'}, multiple=False),
            dcc.Upload(id='upload-intersection', children=html.Button('Upload Intersection File'), style={'fontFamily': '"Times New Roman", Times, serif'}, multiple=False),
            dbc.Button('Process and Calculate Metrics', id='process-button', n_clicks=0, color="primary", className="mt-2", style={'marginTop': '20px', 'fontFamily': '"Times New Roman", Times, serif'}),
            html.Div(id='output-metrics')
        ], style={'display': 'block'}),  # shown only in SURVIVOR mode

        # --- EvalSVcallers Upload Section ---
        html.Div(id='evalsvcallers-upload-section', children=[
            html.Label("Upload EvalSVCallers .eval.txt Output File:",style={'fontFamily': '"Times New Roman", Times, serif'}),
            dcc.Upload(
                id='upload-eval-file',
                children=html.Div(['📎 Drag and Drop or ', html.A('Select *.eval.txt File')]),
                style={
                    'width': '400px',
                    'height': '60px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin': '10px',
                    'fontFamily': '"Times New Roman", Times, serif'
                },
                multiple=False
            ), 
            dbc.Button('Process EvalSVCallers File', id='process-eval-button', n_clicks=0, color="primary", className="mt-2", style={'marginTop': '10px','fontFamily': '"Times New Roman", Times, serif'}),
            html.Div(id="summary-visualization-output")
        ], style={'display': 'none'})  # shown only in EvalSVcallers mode
    ])
'''
