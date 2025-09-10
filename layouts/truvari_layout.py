from dash import html, dcc
import dash_bootstrap_components as dbc
import os

def get_truvari_layout():
    shared_input_style = {
        'fontFamily': '"Times New Roman", Times, serif',
        'width': '200px',
        'marginLeft': '0px'
    }

    shared_label_style = {
        'fontFamily': '"Times New Roman", Times, serif',
        'width': '400px',
        'display': 'inline-block'
    }

    shared_row_style = {
        'marginBottom': '8px'
    }

    # Read .fa files from uploaded_files for reference selector
    ref_options = [
        {"label": "GRCh38", "value": os.path.join("uploaded_files", f)}
        for f in os.listdir("uploaded_files") if "GRCh38" in f and f.endswith(".fa")
    ] + [
        {"label": "HG19", "value": os.path.join("uploaded_files", f)}
        for f in os.listdir("uploaded_files") if ("HG19" in f or "hs37d5" in f) and f.endswith(".fa")
    ]

    return html.Div([
        html.H3("Select Base VCF File (.vcf)", style={'marginTop': '20px', 'fontFamily': '"Times New Roman", Times, serif', 'fontSize': '20px'}),

        dcc.RadioItems(
            id="base_choice",
            options=[
                {"label": "Univar SV Katalog", "value": "univar"},
                {"label": "Upload Custom Reference", "value": "custom"}
                
               
            ],
            value="custom",
            inline=True,
            labelStyle={'display': 'block', 'fontFamily': '"Times New Roman", Times, serif'}
        ),

        
        html.Div(id='tru-base-upload-container', children=[

            dcc.Upload(
                id='tru-base-upload',
                children=html.Div(['📎 Drag and Drop or Select File']),
                style={
                    'width': '100%', 'height': '60px', 'lineHeight': '60px',
                    'borderWidth': '1px', 'borderStyle': 'dashed', 'borderRadius': '5px',
                    'textAlign': 'center', 'margin': '10px', 'fontFamily': '"Times New Roman", Times, serif'
                },
                multiple=False
            )
        ], style={'width': '100%'}),
        
        html.H3("Upload Comparison VCF File (.vcf)", style={'marginTop': '20px', 'fontFamily': '"Times New Roman", Times, serif', 'fontSize': '20px','whiteSpace': 'nowrap' }),
        html.Div([
        
            dcc.Upload(
                id='tru-comp-upload',
                children=html.Div(['📎 Drag and Drop or Select File']),
                style={
                    'width': '100%', 'height': '60px', 'lineHeight': '60px',
                    'borderWidth': '1px', 'borderStyle': 'dashed', 'borderRadius': '5px',
                    'textAlign': 'center', 'margin': '10px', 'fontFamily': '"Times New Roman", Times, serif'
                },
                multiple=False
            )
        ], style={'width': '100%'}),

        html.Div([
            html.Label("Select Reference Genome:", style=shared_label_style),
            dcc.Dropdown(
                id='tru-ref-selector',
                options=ref_options,
                placeholder="Select reference (.fa) from uploaded_files",
                style={'width': '100%'}
            )
        ], style={'width': '100%','marginBottom': '20px'}),

     #   html.H4("Truvari Parameters", style={'fontFamily': '"Times New Roman", Times, serif', 'marginTop': '20px'}),

        html.Div([
            html.Div([
                html.Label("Max reference location distance (--refdist):", style=shared_label_style),
                dcc.Input(id='tru-param-refdist', type='number', min=0, placeholder='default: 500', style=shared_input_style)
            ], style=shared_row_style),

            html.Div([
                html.Label("Min percent allele seq similarity (--pctsim):", style=shared_label_style),
                dcc.Input(id='tru-param-pctsim', type='number', min=0, max=1, step=0.01, placeholder='default: 0.7', style=shared_input_style)
            ], style=shared_row_style),

            html.Div([
                html.Label("Min allele size similarity (--pctsize):", style=shared_label_style),
                dcc.Input(id='tru-param-pctsize', type='number', min=0, max=1, step=0.01, placeholder='default: 0.7', style=shared_input_style)
            ], style=shared_row_style),

            html.Div([
                html.Label("Min reciprocal overlap (--pctovl):", style=shared_label_style),
                dcc.Input(id='tru-param-pctovl', type='number', min=0, max=1, step=0.01, placeholder='default: 0.0', style=shared_input_style)
            ], style=shared_row_style),

            html.Div([
                html.Label("Min variant size to compare (--sizemin):", style=shared_label_style),
                dcc.Input(id='tru-param-sizemin', type='number', placeholder='default: 50', style=shared_input_style)
            ], style=shared_row_style),

            html.Div([
                html.Label("Min variant size to load (--sizefilt):", style=shared_label_style),
                dcc.Input(id='tru-param-sizefilt', type='number', placeholder='default: 30', style=shared_input_style)
            ], style=shared_row_style),

            html.Div([
                html.Label("Max variant size to compare (--sizemax):", style=shared_label_style),
                dcc.Input(id='tru-param-sizemax', type='number', placeholder='default: 50000', style=shared_input_style)
            ], style=shared_row_style),

            html.Div([
                html.Label("Variant type matching (--typeignore)", style=shared_label_style),
                dcc.RadioItems(
                    id='tru-flag-typeignore',
                    options=[{'label': 'Ignore types', 'value': True}, {'label': 'Match types', 'value': False}],
                    value=False, inline=True,labelStyle={'display': 'inline-block', 'marginRight': '10px', 'fontFamily': '"Times New Roman", Times, serif'}
                )
            ], style=shared_row_style),

            html.Div([
                html.Label("Use Levenshtein distance (--use-lev)", style=shared_label_style),
                dcc.RadioItems(
                    id='tru-flag-uselev',
                    options=[{'label': 'Yes', 'value': True}, {'label': 'No', 'value': False}],
                    value=False, inline=True,labelStyle={'display': 'inline-block', 'marginRight': '10px', 'fontFamily': '"Times New Roman", Times, serif'}
                )
            ], style=shared_row_style),

            html.Div([
                html.Label("Compare genotypes (--gtcomp)", style=shared_label_style),
                dcc.RadioItems(
                    id='tru-flag-gtcomp',
                    options=[{'label': 'Yes', 'value': True}, {'label': 'No', 'value': False}],
                    value=False, inline=True,labelStyle={'display': 'inline-block', 'marginRight': '10px', 'fontFamily': '"Times New Roman", Times, serif'}
                )
            ], style=shared_row_style),

            html.Div([
                html.Label("Include only FILTER=PASS (--passonly)", style=shared_label_style),
                dcc.RadioItems(
                    id='tru-flag-passonly',
                    options=[{'label': 'Yes', 'value': True}, {'label': 'No', 'value': False}],
                    value=False, inline=True,labelStyle={'display': 'inline-block', 'marginRight': '10px', 'fontFamily': '"Times New Roman", Times, serif'}
                )
            ], style=shared_row_style),

            html.Div([
                html.Label("Allow multiple matches (--multimatch)", style=shared_label_style),
                dcc.RadioItems(
                    id='tru-flag-multimatch',
                    options=[{'label': 'Yes', 'value': True}, {'label': 'No', 'value': False}],
                    value=False, inline=True,labelStyle={'display': 'inline-block', 'marginRight': '10px', 'fontFamily': '"Times New Roman", Times, serif'}
                )
            ], style=shared_row_style),
        ]),

        dbc.Button('Run Truvari bench', id='tru-run-btn', color='primary', className='mt-3', style={'fontFamily': '"Times New Roman", Times, serif'})

    ])