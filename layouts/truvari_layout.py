from dash import html, dcc
import dash_bootstrap_components as dbc
import os
from layouts.truvari_functions import get_reference_options, REFERENCE_FILES_DIR

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
    ref_options = get_reference_options()

    return html.Div([
        html.H3("Select Base VCF File (.vcf)", style={'marginTop': '20px', 'fontFamily': '"Times New Roman", Times, serif', 'fontSize': '20px'}),

        dcc.RadioItems(
            id="base_choice",
            options=[
                {"label": "Univar SV Katalog", "value": "univar"},
                {"label": "Upload Custom Truth Set", "value": "custom"}
                
               
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
            html.Label("Reference genome input:", style=shared_label_style),
        
            dcc.RadioItems(
                id="tru-ref-input-mode",
                options=[
                    {"label": "Select from reference folder", "value": "folder"},
                    {"label": "Upload custom reference FASTA", "value": "upload"}
                ],
                value="folder",
                inline=False,
                labelStyle={
                    'display': 'block',
                    'fontFamily': '"Times New Roman", Times, serif',
                    'marginBottom': '5px'
                }
            ),
        
            html.Div(
                id="tru-ref-dropdown-container",
                children=[
                    html.Label("Select reference genome:", style=shared_label_style),
        
                    dcc.Dropdown(
                        id='tru-ref-selector',
                        options=ref_options,
                        placeholder="Select reference FASTA from reference_files folder",
                        style={'width': '100%'}
                    ),
        
                    html.Div(
                        id='tru-ref-selected-info',
                        style={
                            'marginTop': '5px',
                            'fontFamily': '"Times New Roman", Times, serif'
                        }
                    ),
        
                    html.Small(
                        "If the required reference genome is not listed, copy the FASTA file into "
                        "`uploaded_files/reference_files/` in a local deployment, or mount this directory in Docker. "
                        "After refreshing the page, the reference FASTA will appear in this dropdown.",
                        style={
                            'display': 'block',
                            'color': 'gray',
                            'fontFamily': '"Times New Roman", Times, serif',
                            'marginTop': '5px',
                            'marginBottom': '10px'
                        }
                    )
                                    ],
                style={'display': 'block', 'width': '100%'}
            ),
        
            html.Div(
                id="tru-ref-upload-container",
                children=[
                    dbc.Alert(
                        [
                            html.Div("Upload mode selected.", style={"fontWeight": "bold"}),
                            html.Div(
                                "After selecting a reference FASTA file, please wait until the upload-completed message appears. "
                                "Large reference files may take several minutes and the browser may appear inactive during upload."
                            ),
                            html.Div(
                                "For large reference genomes, copying the FASTA file into the reference_files folder or mounting it with Docker is recommended."
                            )
                        ],
                        color="warning",
                        style={
                            "fontFamily": '"Times New Roman", Times, serif',
                            "fontSize": "13px",
                            "marginTop": "10px"
                        }
                    ),
            
                    dcc.Upload(
                        id='tru-ref-upload',
                        children=html.Div(['📎 Upload custom reference FASTA (.fa/.fasta/.fna)']),
                        style={
                            'width': '100%',
                            'height': '60px',
                            'lineHeight': '60px',
                            'borderWidth': '1px',
                            'borderStyle': 'dashed',
                            'borderRadius': '5px',
                            'textAlign': 'center',
                            'margin': '10px 0',
                            'fontFamily': '"Times New Roman", Times, serif',
                            'backgroundColor': '#f8f9fa'
                        },
                        multiple=False
                    )
            
                ],
                style={'display': 'none', 'width': '100%'}
            )
        ], style={'width': '100%', 'marginBottom': '20px'}),
                
            html.Div([
                html.Label("High-confident region BED restriction (--includebed):", style=shared_label_style),
            
                dcc.RadioItems(
                    id="tru-bed-mode",
                    options=[
                        {"label": "Do not use BED file", "value": "no_bed"},
                        {"label": "Use BED file", "value": "use_bed"}
                    ],
                    value="no_bed",
                    inline=False,
                    labelStyle={
                        'display': 'block',
                        'fontFamily': '"Times New Roman", Times, serif',
                        'marginBottom': '5px'
                    }
                ),
            
                html.Div(
                    id="tru-bed-upload-container",
                    children=[
                        dcc.Upload(
                            id='tru-bed-upload',
                            children=html.Div(['📎 Drag and Drop or Select BED File']),
                            style={
                                'width': '100%',
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
            
                        html.Div(
                            id='tru-bed-upload-info',
                            style={
                                'marginTop': '5px',
                                'fontFamily': '"Times New Roman", Times, serif'
                            }
                        ),
            
                        html.Small(
                            "If provided, Truvari benchmarking will be restricted to the confident genomic regions in this BED file.",
                            style={
                                'display': 'block',
                                'color': 'gray',
                                'fontFamily': '"Times New Roman", Times, serif',
                                'marginLeft': '10px',
                                'marginBottom': '15px'
                            }
                        )
                    ],
                    style={'display': 'none', 'width': '100%'}
                )
            ], style={'width': '100%', 'marginBottom': '20px'}),

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