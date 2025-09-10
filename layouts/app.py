from dash import Dash, html, dcc, Input, Output, State, ctx, get_asset_url, callback_context
import dash_bio as dashbio
import datetime
import dash_bootstrap_components as dbc
from dash_bootstrap_templates import load_figure_template
from dash.dependencies import MATCH
from dash import no_update
import dash
import flask
import time 
import plotly.express as px
import tempfile
import io
import uuid
import os
import re
import base64
import json
import plotly.graph_objects as go  # ← EK
import pandas as pd
from dash_bio import Clustergram
from layouts.survivor_layout import get_survivor_layout
from layouts.evalsvcallers_layout import get_evalsvcallers_layout
from layouts.visualize_layout import get_visualize_layout
from layouts.metrics_layout import get_metrics_layout
from layouts.truvari_layout import get_truvari_layout
from layouts.truvari_functions import save_uploaded_file, run_truvari_pipeline,_open_any,_parse_vcf_record,_hex_to_rgba,vcf_to_circos_json_truvari,update_tracks,plot_manhattan_truvari
from layouts.survivor_functions import parse_uploaded_files as parse_survivor_uploaded_files, prepare_vcf_files_for_merge, run_survivor_merge, get_merge_preview
from layouts.evalsvcallers_functions import (
    save_file, save_custom_reference_file, parse_uploaded_files as parse_eval_uploaded_files,
    run_conversion, run_evaluation
)
from layouts.visualize_functions import (plot_circos, plot_manhattan, plot_manhattan_svlen, plot_clustergram, plot_sankey,load_vcf_dataframe, update_tracks, vcf_to_circos_json, detect_genome_version, save_file, get_variant_extraction_status, parse_uploaded_vcf, extract_variant_types, plot_vcf_data, serve_file_for_download)
from layouts.metrics_functions import (
    parse_contents, parse_reference_normal, parse_reference_na12878,
    parse_caller_vcf, parse_intersection_vcf, calculate_metrics_explicit,get_survivor_upload_section,get_truvari_upload_section, get_evalsvcallers_upload_section,parse_truvari_file, generate_truvari_visuals, parse_evalsvcallers_file, generate_evalsvcallers_visuals,process_survivor_metrics_from_content, generate_survivor_visuals)
import plotly.express as px
from dash import dash_table

UPLOAD_DIRECTORY = "./uploaded_files/"
SURVIVOR_OUTPUT_DIR= os.path.join(UPLOAD_DIRECTORY, "survivor_output")
EVAL_OUTPUT_DIR = os.path.join(UPLOAD_DIRECTORY, "evalsvcallers_output")
VISUALIZATION_OUTPUT_DIR = os.path.join(UPLOAD_DIRECTORY, "visualization_output")
TRUVARI_OUTPUT_DIR = os.path.join(UPLOAD_DIRECTORY, "truvari_output")
os.makedirs(TRUVARI_OUTPUT_DIR, exist_ok=True)
app_directory = os.path.abspath(os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        os.pardir,
        os.pardir
    ))

assets_directory = os.path.join(app_directory, 'assets')
load_figure_template('LITERA')
app = Dash(__name__, suppress_callback_exceptions=True, external_stylesheets=[dbc.themes.BOOTSTRAP],title="SV Eval & Viz")
active_tab_style = {
    'background': 'linear-gradient(90deg, #0984e3 0%, #74b9ff 100%)',
    'color': 'white',
    'borderRadius': '10px',
    'border': '2px solid #0984e3',
    'fontWeight': '500',
 #   'padding': '0.4rem 0.8rem',
    'boxShadow': '0 2px 6px rgba(0,0,0,0.15)'
}
tab_style = {
    'padding': '0.5rem 1rem',
    'fontSize': '16px',
    'fontFamily': '"Times New Roman", Times, serif'
}
app.layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            dbc.Row([
                dbc.Col(
                    html.Div(
                        html.Img(
                            src=dash.get_asset_url('sv_icon.png'),
                            height='100px',
                            style={
                                'paddingLeft': '4rem'
                            }
                        ),
                        style={
                            'display': 'flex',
                            'alignItems': 'center',
                            'height': '100%'
                        }
                    ),
                    width='auto'
                ),
                dbc.Col(
                    html.Div(
                        html.H2(
                            "SV Eval & Viz",
                            style={
                                'textAlign': 'center',
                                'fontSize': '24px',
                                'color': '#2c3e50',
                                'fontFamily': '\"Times New Roman\", Times, serif',
                                'fontWeight': 'bold',
                                'margin': '0'
                            }
                        ),
                        style={
                            'display': 'flex',
                            'alignItems': 'center',
                            'height': '100%'
                           
                        }
                    )
                ), 
            ], align='center', justify='center')
        ], width='auto')
    ]),
    dbc.Row([
        dbc.Col([
            dbc.Tabs(
                [
                    dbc.Tab(label="Welcome", tab_id="tab-welcome", active_tab_style=active_tab_style),
                    dbc.Tab(label="Comparison", tab_id="tab-comparison", active_tab_style=active_tab_style),
                    dbc.Tab(label="Visualization", tab_id="tab-visualization", active_tab_style=active_tab_style),
                    dbc.Tab(label="Metrics", tab_id="tab-metrics",active_tab_style=active_tab_style),
                ],
                id="main-tabs",
                active_tab="tab-welcome",
                style={
                    "width": "100%",
                    "minWidth": "0",
                    "flexWrap": "wrap",
                    'fontFamily': '\"Times New Roman\", Times, serif',
                    'justifyContent': 'center'
                }
            ),
            html.Div(id='left-extra', style={'marginTop': '1rem'})
        ], width=4, style={
            'backgroundColor': '#f8f9fa',
            'height': 'auto',
            'minHeight': '100vh',
            'padding': '1rem',
            'overflowY': 'auto'
        }),
        dbc.Col(
            html.Div(id='tabs-content', style={'padding': '1rem'}),
            width=8
        )
    ])
], fluid=True)

################### MAIN LAYOUT ROUTING ########################

@app.callback(
    [Output('left-extra', 'children'),
     Output('tabs-content', 'children')],
    Input('main-tabs', 'active_tab'),
)
def update_layout(active_tab):
    if active_tab == 'tab-welcome':
        return html.Div(""), html.Div([
            html.H2("Welcome to the Structural Variant Comparison and Visualization Tool!",
                    style={
                        'padding': '1rem',
                        'fontSize': '20px',
                       # 'fontFamily': '"Times New Roman", Times, serif',
                        'fontWeight': 'bold'
                    }),
        
            html.P("This tool provides a unified platform to compare, visualize, and evaluate Structural Variant (SV) calling outputs using evaluation tools and advanced visualizations.",
                   style={
                       'padding': '0 1rem',
                       'fontSize': '15px'
                      # 'fontFamily': '"Times New Roman", Times, serif'
                   }),
        
            html.H3("🔍 Comparison Tab", style={'padding': '0.75rem 1rem 0 1rem', 'fontSize': '16px', 'fontWeight': 'bold'}),
            html.P("Use this tab to compare SV output with a reference file",
                   style={'padding': '0 1rem', 'fontSize': '15px'}),
            html.Ul([
                html.Li("SURVIVOR: Upload maximum 3 SV caller vcf files to perform merge operation.", style={'padding': '0.25rem'}),
                html.Li("EvalSVcallers: Upload SV caller vcf file for convertion before evaluation and upload/choose a reference file to perform evaluation.", style={'padding': '0.25rem'}),
                 html.Li("Truvari: Upload SV caller vcf file and upload/choose a reference file to perform evaluation.", style={'padding': '0.25rem'}),
            ], style={'padding': '0 2rem', 'fontSize': '15px'}),
        
            html.H3("📊 Visualization Tab", style={'padding': '0.75rem 1rem 0 1rem', 'fontSize': '16px', 'fontWeight': 'bold'}),
            html.P("Create basic or advanced visualizations of output SV data based on caller, SURVIVOR or Evalsvcallers:",
                   style={'padding': '0 1rem', 'fontSize': '15px'}),
            html.Ul([
                html.Li("Basic charts: Data header table, SVTYPE counts, SVLEN distributions, Spyder plot, chromosome-level plots.", style={'padding': '0.25rem'}),
                html.Li("Advanced plots: Sankey, Circos, Clustergram,  and Manhattan visualizations.", style={'padding': '0.25rem'}),
                html.Li("Compatible input types: Caller VCFs (converted via Evalsvcallers convert function), SURVIVOR merged VCF, EvalSVcallers TF VCF.", style={'padding': '0.25rem'}),
            ], style={'padding': '0 2rem', 'fontSize': '15px'}),
        
            html.H3("📈 Metrics Tab", style={'padding': '0.75rem 1rem 0 1rem', 'fontSize': '16px', 'fontWeight': 'bold'}),
            html.P("Calculate and visualize SV calling performance metrics (Precision, Recall, F1-Score):",
                   style={'padding': '0 1rem', 'fontSize': '15px'}),
            html.Ul([
                
                html.Li("EvalSVcallers: Upload the `.eval.txt` metrics file to instantly view tabular and graphical results.", style={'padding': '0.25rem'}),
                html.Li("Truvari: Upload the `summary.txt` metrics file to instantly view tabular and graphical results.", style={'padding': '0.25rem'}),
            ], style={'padding': '0 2rem', 'fontSize': '15px'}),
                
            html.P("Navigate through the tabs to begin your analysis. Hover over graphs for more details, and ensure required files are uploaded in correct formats.",
                   style={'padding': '0 1rem 1rem 1rem', 'fontSize': '15px'
                         # 'fontFamily': '"Times New Roman", Times, serif'
                         })
        ])
    elif active_tab == 'tab-visualization':
        controls = html.Div([
            dcc.RadioItems(
                id='visualize-input-source',
                options=[
                    {'label': 'Caller / Truvari', 'value': 'caller'},
                    {'label': 'SURVIVOR', 'value': 'survivor'},
                    {'label': 'EvalSVcallers', 'value': 'evalsvcallers'}
                   # {'label': 'Truvari', 'value': 'survivor'}
                ],
                value='caller',
                labelStyle={'display': 'block', 'fontFamily': '\"Times New Roman\", Times, serif'}
            ),
            dcc.Upload(
                id='upload-data',
                children=html.Div(['📎 Drag and Drop or ', html.A('Select *vcf File')]),
                style={
                    'width': '100%', 'height': '60px', 'lineHeight': '60px',
                    'borderWidth': '1px', 'borderStyle': 'dashed', 'borderRadius': '5px',
                    'textAlign': 'center', 'margin': '10px 0', 'fontFamily': '\"Times New Roman\", Times, serif'
                },
                multiple=False
            ),
            html.Div(id='vcf-files-list', style={'marginTop': '10px'}),
            dbc.Spinner([dbc.Button("Extract Variant Types", id='visualize-button', n_clicks=0, color="primary", className="mt-2", style={'fontFamily': '"Times New Roman", Times, serif'})], size="sm", color="primary", type="border", fullscreen=False),           
            html.Div(id='visualize-status', style={'marginTop': '10px'}),
            html.Label("Filter by Variant Type:", style={'marginTop': '10px','fontFamily': '"Times New Roman", Times, serif'}),
            dcc.Dropdown(id="variant-type-filter", options=[], multi=False, placeholder="Select a variant type", style={'fontFamily': '"Times New Roman", Times, serif'}),
            html.H2("Select Advanced Visualizations", style={'marginTop': '10px','fontSize': '24px', 'fontFamily': '"Times New Roman", Times, serif'}),
            dcc.RadioItems(
                id="viz_selector",
                options=[
                    {"label": " Data Table & Basic Visuals", "value": "basic"},
                    {"label": " Sankey", "value": "sankey"},
                    {"label": " Circos", "value": "circos"},
                    {"label": " Clustergram", "value": "clustergram"},
                    {"label": " Manhattan", "value": "manhattan"},
                ],
                value=None,
                labelStyle={"display": "block", 'fontFamily': '"Times New Roman", Times, serif'}
            ),
            dcc.Store(id='selected-input-source', data=""),
            dcc.Store(id='uploaded-file-path', data=""),
            dcc.Store(id='filtered-data', data="")
           
        ])
        return controls, get_visualize_layout()

    elif active_tab == 'tab-metrics':
        radio = html.Div([
            html.P("Select a tool from below.",style={'padding': '0.5rem', 'fontSize': '15px', 'fontFamily': '\"Times New Roman\", Times, serif'}),
            dcc.RadioItems(
                id='metrics-input-source',
                options=[
                  #  {'label': 'SURVIVOR', 'value': 'survivor'},
                    {'label': 'EvalSVcallers', 'value': 'evalsvcallers'},
                    {'label': 'Truvari', 'value': 'truvari'},
                ],
                value=None,
                labelStyle={'display': 'block', 'fontFamily': '\"Times New Roman\", Times, serif'}
            ),
            html.Div(id='metrics-upload-section') ])
        return radio, get_metrics_layout()

    elif active_tab == 'tab-comparison':
        left_column = html.Div([
            html.P("Select a tool from below.",
                   style={'padding': '0.5rem', 'fontSize': '15px', 'fontFamily': '"Times New Roman", Times, serif'}),
            
            html.Div([
                dcc.RadioItems(
                    id='comparison-type',
                    options=[
                        {'label': 'SURVIVOR', 'value': 'survivor'},
                        {'label': 'EvalSVcallers', 'value': 'evalsvcallers'},
                        {'label': 'Truvari', 'value': 'truvari'} 
                    ],
                    value=None,
                    inline=True,
                    labelStyle={'marginRight': '20px', 'marginLeft': '20px','fontFamily': '"Times New Roman", Times, serif'}
                )
            ], style={'whiteSpace': 'nowrap'}),
    
            html.Div(id="survivor-params-container", children=get_survivor_layout(), style={'display': 'none'}),
            html.Div(id="evalsvcallers-params-container", children=get_evalsvcallers_layout(), style={'display': 'none'}),
            html.Div(id="truvari-params-container", children=get_truvari_layout(), style={'display': 'none'})
        ], style={'width': '100%', 'display': 'inline-block', 'verticalAlign': 'top'})

    
        right_column = html.Div([
            dcc.Store(id='survivor-results-store', data=None),
            dcc.Store(id='evalsvcallers-results-store', data=None),

            dcc.Store(id="merged-path-store"), #survivor
            dcc.Store(id='converted-file-store'), #evalsvcallers
            dcc.Store(id="tp-fp-file-store"),
            dcc.Store(id="metrics-file-store"),
            dcc.Store(id='survivor-uploaded-paths', data=[]), #uploaded-files-list #converted-files-list
           # dcc.Store(id='survivor-uploaded-paths-ref', data=[]),#reference_upload_section
       
            html.Div(id='tru-base-upload-info', style={'marginTop': '0px'}),
            html.Div(id='tru-comp-upload-info', style={'marginTop': '0px'}),
            html.Div(id='tru-output-list', style={'marginTop': '0px'}),
            dcc.Store(id='truvari-uploaded-vcf-paths'),

            
            html.Div(id='eval-output-ref', style={'marginTop': '20px'}),
            # Output area
            dcc.Loading(
                id="loading-survivor-output",
                type="default",  # or "circle", "dot", "graph"
                fullscreen=False,
                children=html.Div([
                    html.Div(id='survivor-output', style={'marginTop': '20px'}), #convert-status
                    html.Div(id='eval-output-list', style={'marginTop': '10px', 'fontFamily': '"Times New Roman", Times, serif'}),
                    html.Div(id='tru-status', style={'marginTop': '10px'}),
                    html.Div(id="truvari-metrics-preview", style={"marginTop": "10px"}),
                    html.Div(id='truvari-visual-output', style={"marginTop": "10px"})

                ])
                
            ),
            # Separate visualization containers
            html.Div(id='survivor-visuals-container', children=[
                html.Div(id='comparison-extra-output')
            ], style={'display': 'block'}),


            html.Div(id='evalsvcallers-visuals-container', children=[
                
            ], style={'display': 'none'}),
            html.Div(id='truvari-visuals-container', children=[
                html.Div(id='truvari-output-list'),  # veya başka görsel alan
            ], style={'display': 'none'})

        ])
        return left_column, right_column

    else:
        return html.Div(""), html.Div("Unknown tab.")
# --- Place the toggle callback after layout ---
@app.callback(
    Output('survivor-visuals-container', 'style'),
    Output('evalsvcallers-visuals-container', 'style'),
    Output('truvari-visuals-container', 'style'),
    Input('comparison-type', 'value')
)
def toggle_visuals(selected_source):
    if selected_source == 'survivor':
        return {'display': 'block'}, {'display': 'none'}, {'display': 'none'}
    elif selected_source == 'evalsvcallers':
        return {'display': 'none'}, {'display': 'block'}, {'display': 'none'}
    elif selected_source == 'truvari':
        return {'display': 'none'}, {'display': 'none'}, {'display': 'block'}
    return {'display': 'none'}, {'display': 'none'}, {'display': 'none'}
@app.callback(
    Output('survivor-params-container', 'style'),
    Output('evalsvcallers-params-container', 'style'),
    Output('truvari-params-container', 'style'),  # ✅ yeni çıktı
    Input('comparison-type', 'value'),
)
def toggle_param_layout(selected_tool):
    if selected_tool == 'survivor':
        return {'display': 'block'}, {'display': 'none'},{'display': 'none'}
    elif selected_tool == 'evalsvcallers':
        return {'display': 'none'}, {'display': 'block'}, {'display': 'none'}
    elif selected_tool == 'truvari':
        return {'display': 'none'}, {'display': 'none'}, {'display': 'block'}
    return {'display': 'none'}, {'display': 'none'}, {'display': 'none'}
@app.callback(
    Output('converted-file-store', 'data', allow_duplicate=True),
    Output('tp-fp-file-store', 'data', allow_duplicate=True),
    Output('metrics-file-store', 'data', allow_duplicate=True),
    Output('survivor-uploaded-paths', 'data', allow_duplicate=True),
    Output('merged-path-store', 'data', allow_duplicate=True),
    Output('truvari-output-list', 'children', allow_duplicate=True),
    Output('tru-status', 'children', allow_duplicate=True),
 #   Output('eval-preview-card', 'children', allow_duplicate=True),
 #   Output('eval-visuals-output', 'children', allow_duplicate=True),
 #   Output('eval-clustergram-output', 'children', allow_duplicate=True),
 #   Output('eval-circos-output', 'children', allow_duplicate=True),
 #   Output('eval-manhattan-output', 'children', allow_duplicate=True),
 #   Output('preview-card', 'children', allow_duplicate=True),
 #   Output('visuals-output', 'children', allow_duplicate=True),
  #  Output('clustergram-output', 'children', allow_duplicate=True),
  #  Output('circos-output', 'children', allow_duplicate=True),
  #  Output('manhattan-output', 'children', allow_duplicate=True),
    Input('comparison-type', 'value'),
    prevent_initial_call=True
)
def reset_all_outputs(selected_tool):
    return (
        {}, {}, {},
        [], None,
        html.Div(),  # truvari-output-list
        html.Div()   # tru-status
    )       
################### SURVIVOR FUNCTIONS ########################
@app.callback(
    Output('ref-upload-section', 'style'),
    Input('ref_source_selector', 'value')
)
def toggle_ref_upload(ref_source):
    if ref_source == 'upload':
        return {'display': 'block'}
    else:
        return {'display': 'none'}
@app.callback(
    Output('survivor-output', 'children'),
    Input('comparison-type', 'value'),
    prevent_initial_call=True
)
def clear_survivor_status_on_switch(tool):
    if tool == 'evalsvcallers':
        return ""  # clears previous eval messages
    return dash.no_update    

@app.callback(
    Output('preview-card', 'children'),
    Output('visuals-output', 'children'),
    Output('clustergram-output', 'children'),
    Output('circos-output', 'children'),
    Output('manhattan-output', 'children'),
    [
        Input('comparison-svtype-filter', 'value'),
        Input('clustergram-chromosome-selector', 'value'),
        Input('graph-type-dropdown', 'value'),
        Input('svtype-filter-dropdown', 'value'),
        Input('manhattan-svtype-selector', 'value'),
        Input('manhattan-threshold-slider', 'value'),
    ],
    State("survivor-uploaded-paths", "data"),
    State("merged-path-store", "data"),
    prevent_initial_call=True
)

def update_survivor_comparison_visuals(selected_svtype, selected_chroms, circos_graph_type, circos_svtypes, manhattan_svtypes, manhattan_threshold, uploaded_paths,merged_path):

    
    
    # 🚀 Early exit block
    if not uploaded_paths or not merged_path:
        return html.Div("⚠️ No merged file available."), html.Div(), html.Div(), html.Div(), html.Div() 
        
    vcf_filename = os.path.basename(merged_path)
    json_filename = f"{vcf_filename.replace('.vcf', '')}_circos.json"
    json_path = os.path.join(SURVIVOR_OUTPUT_DIR, json_filename)

    with open(json_path, "r") as f:
        circos_data = json.load(f)

    svtype_colors_circos = {track["name"]: track.get("color", "#95a5a6") for track in circos_data.get("tracks", [])}
    
    df = load_vcf_dataframe(merged_path, "survivor")

    if df is None or df.empty:
        return html.Div("⚠️ No data for preview.")

    # Ensure SVTYPE is correctly extracted from INFO
    df['SVTYPE'] = df['INFO'].str.extract(r'SVTYPE=([^;]+)')

    # Filter by selected SVTYPE if not "ALL"
    df_filtered = df.copy()
    if selected_svtype and selected_svtype != "ALL":
        df_filtered = df_filtered[df_filtered["SVTYPE"] == selected_svtype]

    # Table preview
    table_card = dbc.Card([
        dbc.CardHeader(html.H4(f"VCF File Preview ({selected_svtype})")),
        dbc.CardBody([
            dash_table.DataTable(
                data=df_filtered.head(5).to_dict('records'),
                columns=[{"name": col, "id": col} for col in df.columns],
                style_table={'overflowX': 'auto'},
                style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'},
                style_header={'fontWeight': 'bold', 'backgroundColor': 'lightgrey', 'border': '1px solid black'},
                style_data={'border': '1px solid black'}
            )
        ])
    ], className="mt-4")


    # Basic visuals plot
    visuals_html = plot_vcf_data(df_filtered)

    
    # Clustergram: ensure selected chromosomes exist
    valid_chroms = [f"chr{i}" for i in range(1,23)] + ["chrX","chrY"]
    df["CHROM"] = df["CHROM"].astype(str).apply(lambda x: f"chr{x}" if not x.startswith("chr") else x)
    available_chroms = [c for c in valid_chroms if c in df["CHROM"].unique()]

    selected_existing_chroms = [chrom for chrom in selected_chroms if chrom in available_chroms]
    if not selected_existing_chroms:
        clustergram_html = html.Div("⚠️ No selected chromosomes found in data.")
    else:
        clustergram_html = plot_clustergram(df, selected_existing_chroms)


    # Circos plot
    circos_html = plot_circos(circos_graph_type, circos_svtypes, circos_data, svtype_colors_circos)
    if not manhattan_svtypes:
        manhattan_svtypes = df["SVTYPE"].dropna().unique().tolist()
    
    if not manhattan_threshold:
        manhattan_threshold = 6
    
    # Manhattan plot
    manhattan_html = plot_manhattan(df, manhattan_svtypes, manhattan_threshold)

    return table_card, visuals_html, clustergram_html, circos_html, manhattan_html
    
@app.callback(
    Output('survivor-output', 'children', allow_duplicate=True),
    Output('survivor-uploaded-paths', 'data'),
  #  Output('survivor-uploaded-paths-ref', 'data'),
    Output('comparison-extra-output', 'children', allow_duplicate=True),
    Output('merged-path-store', 'data'),
    [Input('survivor-upload', 'contents'),
#     Input('survivor-upload-ref', 'contents'),
     Input('merge-button', 'n_clicks')],
    [State('survivor-upload', 'filename'),
     State('survivor-upload', 'last_modified'),
 #    State('survivor-upload-ref', 'filename'),
 #    State('survivor-upload-ref', 'last_modified'),
  #   State('ref_source_selector', 'value'),
     State('param_dist', 'value'),
     State('param_callers', 'value'),
     State('param_type', 'value'),
     State('param_strand', 'value'),
     State('param_dup', 'value'),
     State('param_sv_size', 'value'),
     State('survivor-uploaded-paths', 'data'),
   #  State('survivor-uploaded-paths-ref', 'data'),
     State('main-tabs', 'active_tab')],
    prevent_initial_call=True
)
def handle_survivor_actions(contents, n_clicks,
                            names, dates,
                            dist, callers, sv_type, strand, dup, sv_size,
                            uploaded_paths, active_tab):
    comparison_output = dash.no_update
    merged_path= dash.no_update
    triggered_id = ctx.triggered_id
    
    # Handle uploads
    if triggered_id == 'survivor-upload' and contents:
        if isinstance(contents, str):
            contents, names, dates = [contents], [names], [dates]
        div, new_paths = parse_survivor_uploaded_files(contents, names, dates)
     #   div = html.Div([
    #        div, 
    #        html.Br(),
    #        html.H5("Uploaded Files:"),
     #       html.Ul([html.Li(os.path.basename(p)) for p in uploaded_paths + new_paths])
     #   ])
        
        return div, uploaded_paths + new_paths, dash.no_update, dash.no_update

 #   if triggered_id == 'survivor-upload-ref' and ref_contents:
  #      if isinstance(ref_contents, str):
  #          ref_contents, ref_filename, ref_dates = [ref_contents], [ref_filename], [ref_dates]
  #      div, new_paths = parse_survivor_uploaded_files(ref_contents, ref_filename, ref_dates)
  #      return div, dash.no_update, new_paths[-1:],dash.no_update, dash.no_update

    # Handle merge
    if triggered_id == 'merge-button':
        
        params = {
#            'use_univar': 1 if ref_source == 'univar' else 0,
            'max_distance': int(dist) if dist else 1000,
            'min_callers': callers,
            'type_match': int(sv_type) if sv_type else 1,
            'strand_match': int(strand) if strand else 1,
            'allow_duplicates': int(dup) if dup else 0,
            'min_sv_size': int(sv_size) if sv_size else 30,
        }

        vcf_files=[]
        if uploaded_paths:
            vcf_files = [os.path.abspath(p) for p in uploaded_paths]
        else:
            vcf_files = []
  #      if params['use_univar'] == 1:
 #           univar_path = os.path.abspath(os.path.join(UPLOAD_DIRECTORY, "univar.vcf"))
 #           if os.path.exists(univar_path):
  #              vcf_files.append(univar_path)
  #      elif uploaded_paths_ref:
  #          vcf_files.append(os.path.abspath(uploaded_paths_ref[-1]))

     #   if len(set(vcf_files)) > 3:
      #      return html.Div("⚠️ More than 3 VCF files uploaded, last 3 VCF files are selected for analysis."), uploaded_paths, uploaded_paths_ref, dash.no_update, dash.no_update

        
        # Run SURVIVOR merge
        status, merged_path = run_survivor_merge(params, vcf_files)

        if not merged_path:
            return html.Div([status]), uploaded_paths, dash.no_update, None

        # ✅ Now load to get SVTYPEs for dropdowns
        df = load_vcf_dataframe(merged_path, "survivor")
        if df is None or df.empty:
            return html.Div([status, html.Br(), html.Div("⚠️ Merged file loaded but no variants found. Make sure reference version is hg38")]), uploaded_paths, dash.no_update, dash.no_update

        df["SVTYPE"] = df["INFO"].str.extract(r"SVTYPE=([^;]+)").fillna("UNK")
        unique_svtypes = sorted(df["SVTYPE"].dropna().unique())

        
        if active_tab == 'tab-comparison' and merged_path:
            try:
                # Helper to convert file -> base64
                def read_file_as_base64(file_path):
                    with open(file_path, "rb") as f:
                        encoded = base64.b64encode(f.read()).decode('utf-8')
                    return "data:text/plain;base64," + encoded

                inter_content_b64 = read_file_as_base64(merged_path)
                # Read files for evaluation
         #       if ref_source == 'univar':
          #          univar_path = os.path.abspath(os.path.join(UPLOAD_DIRECTORY, "univar.vcf"))
         #           if not os.path.exists(univar_path):
         #               raise FileNotFoundError(f"Univar VCF not found at {univar_path}")
         #           ref_content_b64 = read_file_as_base64(univar_path)
         #       else:
         #           ref_content_b64 = read_file_as_base64(uploaded_paths_ref[-1]) if uploaded_paths_ref else None
                
                caller_content_b64 = read_file_as_base64(uploaded_paths[-1]) if uploaded_paths else None
                
         #       if not (ref_content_b64 and caller_content_b64):
         #           raise ValueError("Uploaded caller or reference file content missing.")

                # Metrics
         #       metrics_df, error_msg = process_survivor_metrics_from_content(
          #          ref_content_b64, caller_content_b64, inter_content_b64, 'normal'
         #       )
        #        metrics_html = (
          #          generate_survivor_visuals(metrics_df)
        #            if not error_msg else html.Div(error_msg, style={'color': 'red'})
         #       )

                # Basic visuals
                df = load_vcf_dataframe(merged_path, "survivor")
                # Correct SVTYPE
                df['SVTYPE'] = df['INFO'].str.extract(r'SVTYPE=([^;]+)')
    
                if df is not None:
                
                    # Combine preview table + your visuals into one html.Div
                    visuals_html = html.Div([
                  #      preview_table,
                        html.Div(plot_vcf_data(df))
                    ])
                else:
                    visuals_html = html.Div("⚠️ No data for visualization.")

                # Advanced plots
                try:
                    sankey_html = plot_sankey(df)
                except Exception as e:
                    sankey_html = html.Div(f"❌ Sankey error: {str(e)}")

                try:
                    valid_chroms = [f"chr{i}" for i in range(1,23)] + ["chrX","chrY"]
                    df["CHROM"] = df["CHROM"].astype(str).apply(
                        lambda x: f"chr{x}" if not x.startswith("chr") else x
                    )
                    available_chroms = [c for c in valid_chroms if c in df["CHROM"].unique()]
                    clustergram_html = plot_clustergram(df, available_chroms) if available_chroms else html.Div("⚠️ No chromosomes for Clustergram.")
                except Exception as e:
                    clustergram_html = html.Div(f"❌ Clustergram error: {str(e)}")

                try:
                    vcf_filename = os.path.basename(merged_path)
                    timestamp = time.strftime("%Y%m%d_%H%M%S")
                    json_filename = f"{vcf_filename.replace('.vcf', '')}_circos.json"
                    json_path = os.path.join(SURVIVOR_OUTPUT_DIR, json_filename)
                
                    vcf_to_circos_json(merged_path, json_path, "survivor")
                    with open(json_path, "r") as f:
                        circos_data = json.load(f)
                    available_svtypes_circos = [track["name"] for track in circos_data.get("tracks", [])]
                    svtype_colors_circos = {track["name"]: track.get("color", "#95a5a6") for track in circos_data.get("tracks", [])}
                
                    circos_html = plot_circos("histogram", available_svtypes_circos, circos_data, svtype_colors_circos)
                except Exception as e:
                    circos_html = html.Div(f"❌ Circos error: {str(e)}")
 
                try:
                    all_svtypes = df["SVTYPE"].dropna().unique().tolist()
                    manhattan_html = plot_manhattan(df, all_svtypes, 6)
                except Exception as e:
                    manhattan_html = html.Div(f"❌ Manhattan error: {str(e)}")

                # Build combined output

                comparison_output = html.Div([
                    html.H3("📊 Visualizations", style={'marginTop': '20px'}),
                #    metrics_html,
                    html.Br(),
                    # General filters for all visuals
                    html.Div([
                        html.Label("Filter by Variant Type:"),
                        dcc.Dropdown(
                            id="comparison-svtype-filter",
                            options=[{"label": "ALL", "value": "ALL"}] + [{"label": sv, "value": sv} for sv in unique_svtypes],
                            value="ALL",
                            clearable=False,
                            style={"width": "300px"}
                        )
                    ], style={"marginBottom": "20px"}),
                    html.Div([
                        dbc.Card([
                            dbc.CardHeader(html.H4(f"VCF File Preview (ALL)")),
                            dbc.CardBody([
                                dash_table.DataTable(
                                    data=df.head(5).to_dict('records'),
                                    columns=[{"name": col, "id": col} for col in df.columns],
                                    style_table={'overflowX': 'auto'},
                                    style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'},
                                    style_header={'fontWeight': 'bold', 'backgroundColor': 'lightgrey', 'border': '1px solid black'},
                                    style_data={'border': '1px solid black'}
                                )
                            ])
                        ], className="mt-4")
                    ], id="preview-card"),
                    html.Br(),
                    # Basic visualizations block
                    html.H4("Basic Visualizations"),
                    html.Div(visuals_html, id="visuals-output"),
                    html.Br(),
                
                    # Advanced Visuals
                    html.H4("Advanced Visualizations"),
                    sankey_html,
                    html.Br(),
                
                    # Clustergram filters + Clustergram
                    html.H4("Clustergram Visualization"),
                    html.Div([
                        html.Label("Select Chromosomes for Clustergram:"),
                        dcc.Dropdown(
                            id="clustergram-chromosome-selector",
                            options=[{"label": c, "value": c} for c in available_chroms],
                            value=available_chroms,
                            multi=True,
                            style={"width": "600px"}
                        )
                    ], style={"marginBottom": "20px"}),
                    html.Div(clustergram_html, id="clustergram-output"),
                    html.Br(),
                
                    # Circos filters + Circos
                    html.H4("Circos Visualization"),
                    html.Div([
                        html.Label("Graph Type:"),
                        dcc.Dropdown(
                            id="graph-type-dropdown",
                            options=[{"label": "Histogram", "value": "histogram"}],
                            value="histogram",
                            clearable=False,
                            style={"width": "200px"}
                        ),
                        html.Br(),
                        html.Label("Select SV Types:"),
                        dcc.Dropdown(
                            id="svtype-filter-dropdown",
                            options=[{"label": sv, "value": sv} for sv in unique_svtypes],
                            value=unique_svtypes,
                            multi=True,
                            style={"width": "400px"}
                        )
                    ], style={"marginBottom": "20px"}),
                    html.Div(circos_html, id="circos-output"),
                    html.Br(),
                
                    # Manhattan filters + Manhattan
                    html.H4("Manhattan Visualization"),
                    html.Div([
                        html.Label("Select SV Types for Manhattan Plot:"),
                        dcc.Dropdown(
                            id="manhattan-svtype-selector",
                            options=[{"label": sv, "value": sv} for sv in unique_svtypes],
                            value=unique_svtypes,
                            multi=True,
                            style={"width": "400px"}
                        ),
                        html.Label("Manhattan Threshold:"),
                        dcc.Slider(
                            id="manhattan-threshold-slider",
                            min=0,
                            max=10,
                            step=1,
                            value=6,
                            marks={i: str(i) for i in range(11)},
                            tooltip={"always_visible": False, "placement": "bottom"}
                            
                        )
                    ], style={"marginBottom": "20px"}),
                    html.Div(manhattan_html, id="manhattan-output")
                    
                ])

            except Exception as e:
                comparison_output = html.Div(f"❌ Error generating visuals: {str(e)}")

        
        # ✅ Return merge status + filtering controls (visuals filled by next callback)
        return html.Div([status, html.Br()]), uploaded_paths, comparison_output, merged_path

    return dash.no_update, uploaded_paths,comparison_output, merged_path
################### EVALSVCALLERS FUNCTIONS ########################
@app.callback(
    Output('conversion-section', 'style', allow_duplicate=True),
    Input('comparison-type', 'value'),
    prevent_initial_call=True
)
def show_conversion_if_evalsvcallers(selected_tool):
    if selected_tool == 'evalsvcallers':
        return {'display': 'block'}
    return {'display': 'none'}

@app.callback(
    Output('eval-output-ref', 'children'),
    Input('dynamic-upload', 'contents'),
    State('dynamic-upload', 'filename'),
    prevent_initial_call=True
)
def handle_reference_upload(contents, filename):
    if not contents:
        return "⚠️ Please upload reference file."

    result_msg = save_custom_reference_file(contents, filename)
    return html.Div(result_msg, style={'color': 'green' if '✅' in result_msg else 'red'})


# ========== Upload + Convert ==========
@app.callback(
    Output('eval-output-list', 'children', allow_duplicate=True),
    Output('converted-file-store', 'data'),
    [Input('eval-upload', 'contents'),
     Input('upload-converted-data', 'contents'),
     Input('convert-button', 'n_clicks')],
    [State('eval-upload', 'filename'),
     State('eval-upload', 'last_modified'),
     State('upload-converted-data', 'filename'),
     State('upload-converted-data', 'last_modified'),
     State('caller_tool', 'value')],
    prevent_initial_call=True
)
def handle_upload_or_convert(
    eval_contents, uploaded_contents, convert_clicks,
    eval_name, eval_date, uploaded_name, uploaded_date, caller_tool
):
    triggered_id = ctx.triggered_id

    # Uploading already converted file
    if 'upload-converted-data' in str(triggered_id) and uploaded_contents:
        parse_eval_uploaded_files([uploaded_contents], [uploaded_name], [uploaded_date])
        filename = uploaded_name
        return html.Div(f"✅ Uploaded converted file: {filename}"), {"filename": filename, "timestamp": None}


    # Uploading raw file in convert mode
    if 'eval-upload' in str(triggered_id) and eval_contents:
        parse_eval_uploaded_files([eval_contents], [eval_name], [eval_date])
        return html.Div(f"✅ Uploaded raw file: {eval_name}"), html.Div()

    # Clicking Convert button
    if triggered_id == 'convert-button':
        if eval_name:
            filename = eval_name
            result, converted_filename, timestamp = run_conversion(caller_tool, filename)
            return html.Div([result]), {"filename": converted_filename, "timestamp": timestamp}
        else:
            return html.Div("❌ Please upload a file before converting."), html.Div()

    return html.Div(), html.Div()
    
@app.callback(
    Output('advanced_params', 'style'),
    Input('param_mode', 'value')
)
def toggle_param_mode(mode):
    if mode == 'advanced':
        return {'display': 'block'}
    return {'display': 'none'}
@app.callback(
    Output('eval-output-list', 'children'),
    Input('comparison-type', 'value'),
    prevent_initial_call=True
)
def clear_eval_status_on_switch(tool):
    if tool == 'survivor':
        return html.Div()   # clears previous eval messages
    return html.Div()
@app.callback(
    Output('eval-output-list', 'children', allow_duplicate=True),
    Input('upload-converted-data', 'filename'),
    prevent_initial_call=True
)
def clear_eval_message_on_file_upload(filename):
    if filename:
        return html.Div()  # clear message area
    return html.Div()
@app.callback(
    Output('eval-manhattan-selector', 'options'),
    [Input('tp-fp-file-store', 'data')],
    [State('comparison-type', 'value')]
)
def update_manhattan_svtypes(data, comparison_type):
    if comparison_type != 'evalsvcallers':
        raise dash.exceptions.PreventUpdate

    if not data or data in ["null", "None", "", {}]:
        raise dash.exceptions.PreventUpdate

    try:
        df = pd.read_json(data, orient="split")
    except (ValueError, TypeError):
        raise dash.exceptions.PreventUpdate

    if "SVTYPE" not in df.columns:
        return []

    svtypes = df["SVTYPE"].dropna().unique().tolist()
    return [{"label": sv, "value": sv} for sv in svtypes]
    
    
@app.callback(
    Output('eval-manhattan-slider', 'max'),
    Output('eval-manhattan-slider', 'marks'),
    [Input('tp-fp-file-store', 'data')],
    [State('comparison-type', 'value')]
)
def update_eval_threshold_slider_range(data, comparison_type):
    if comparison_type != 'evalsvcallers':
        raise dash.exceptions.PreventUpdate

    if not data or data in ["null", "None", "", {}]:
        raise dash.exceptions.PreventUpdate

    try:
        df = pd.read_json(data, orient="split")
    except (ValueError, TypeError):
        raise dash.exceptions.PreventUpdate

    return 10, {i: str(i) for i in range(11)}
    

@app.callback(
    Output('eval-preview-card', 'children'),
    Output('eval-visuals-output', 'children'),
    Output('eval-clustergram-output', 'children'),
    Output('eval-circos-output', 'children'),
    Output('eval-manhattan-output', 'children'),
    [
        Input('eval-svtype-filter', 'value'),
        Input('eval-clustergram-selector', 'value'),
        Input('eval-graph-type-dropdown', 'value'),
        Input('eval-svtype-dropdown', 'value'),
        Input('eval-manhattan-selector', 'value'),
        Input('eval-manhattan-slider', 'value'),
    ],
    State("tp-fp-file-store", "data"),
    prevent_initial_call=True
)
def update_eval_comparison_visuals(selected_svtype, selected_chroms, circos_graph_type,
                                   circos_svtypes, manhattan_svtypes, manhattan_threshold,
                                   tp_fp):
    if not tp_fp or isinstance(tp_fp, dict):
        return html.Div("⚠️ No evaluation file loaded."), html.Div(), html.Div(), html.Div(), html.Div()

    # Load Circos JSON for eval
    vcf_filename = os.path.basename(tp_fp)
    json_filename = f"{vcf_filename}_circos.json"
    json_path = os.path.join(EVAL_OUTPUT_DIR, json_filename)
    with open(json_path, "r") as f:
        circos_data = json.load(f)
    svtype_colors_circos = {track["name"]: track.get("color", "#95a5a6") for track in circos_data.get("tracks", [])}

    # Load DataFrame
    df = load_vcf_dataframe(tp_fp, "evalsvcallers")
    if df is None or df.empty:
        return html.Div("⚠️ No data available."), html.Div(), html.Div(), html.Div(), html.Div()
     
    # Ensure SVTYPE is correctly extracted from INFO
    df['SVTYPE'] = df['INFO'].str.extract(r'SVTYPE=([^;]+)')
    # Filter by SVTYPE
    df_filtered = df.copy()
    if selected_svtype and selected_svtype != "ALL":
        df_filtered = df_filtered[df_filtered["SVTYPE"] == selected_svtype]

    # Table preview
    table_card = dbc.Card([
        dbc.CardHeader(html.H4(f"VCF File Preview ({selected_svtype})")),
        dbc.CardBody([
            dash_table.DataTable(
                data=df_filtered.head(5).to_dict('records'),
                columns=[{"name": col, "id": col} for col in df.columns],
                style_table={'overflowX': 'auto'},
                style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'},
                style_header={'fontWeight': 'bold', 'backgroundColor': 'lightgrey', 'border': '1px solid black'},
                style_data={'border': '1px solid black'}
            )
        ])
    ], className="mt-4")

    # Basic visuals
    visuals_html = plot_vcf_data(df_filtered)

    # Clustergram
    valid_chroms = [f"chr{i}" for i in range(1,23)] + ["chrX","chrY"]
    df["CHROM"] = df["CHROM"].astype(str).apply(lambda x: f"chr{x}" if not x.startswith("chr") else x)
    available_chroms = [c for c in valid_chroms if c in df["CHROM"].unique()]
    selected_existing_chroms = [chrom for chrom in selected_chroms if chrom in available_chroms]

    if not selected_existing_chroms:
        clustergram_html = html.Div("⚠️ No selected chromosomes found in data.")
    else:
        clustergram_html = plot_clustergram(df, selected_existing_chroms)

    # Circos
    circos_html = plot_circos(circos_graph_type, circos_svtypes, circos_data, svtype_colors_circos)
    if not manhattan_svtypes:
        manhattan_svtypes = df["SVTYPE"].dropna().unique().tolist()
    
    if not manhattan_threshold:
        manhattan_threshold = 6
    # Manhattan
    manhattan_html = plot_manhattan(df, manhattan_svtypes, manhattan_threshold)

    return table_card, visuals_html, clustergram_html, circos_html, manhattan_html
@app.callback(
    Output('eval-output-list', 'children', allow_duplicate=True),
    Output('tp-fp-file-store', 'data'),
    Output('metrics-file-store', 'data'),
    Output('evalsvcallers-visuals-container', 'children'),
    Input('eval-button', 'n_clicks'),
    State('converted-file-store', 'data'),
    State("reference_choice", "value"),
    State('caller_tool', 'value'),

    # BASIC param states
    State('param_st_basic', 'value'),
    State('param_l_basic', 'value'),
    State('param_xl_basic', 'value'),
    State('param_rl_basic', 'value'),
    State('param_rxl_basic', 'value'),
    State('param_mo_basic', 'value'),
    State('param_of_basic', 'value'),

    # ADVANCED param states
    State('param_parent1', 'contents'),
    State('param_parent1', 'filename'),
    State('param_parent2', 'contents'),
    State('param_parent2', 'filename'),
    State('param_rb', 'contents'),
    State('param_rb', 'filename'),
    State('param_c', 'value'),
    State('param_mr', 'value'),
    State('param_rs', 'value'),
    State('param_mins', 'value'),
    State('param_eg', 'value'),
    State('param_eb', 'value'),
    State('param_i', 'value'),
    State('param_y', 'value'),
    State('param_sm', 'value'),
    prevent_initial_call=True
)
def run_evalsvcallers_combined(
    n_clicks, selected_converted_file, reference_choice, caller_tool,
    st_basic, l_basic, xl_basic, rl_basic, rxl_basic, mo_basic, of_basic,
    parent1_content, parent1, parent2_content, parent2, rb_content, rb,
    c, mr, rs, mins, eg, eb, i, y, sm
):
    try:
        if not n_clicks or not selected_converted_file:
            return html.Div(), [], [], html.Div()

        # Run the evaluation pipeline
        #selected_converted_file_path = os.path.join(EVAL_OUTPUT_DIR, selected_converted_file)
         # unroll dict
        if isinstance(selected_converted_file, dict):
            converted_filename = selected_converted_file.get("filename")
            timestamp = selected_converted_file.get("timestamp")
        else:
            converted_filename = selected_converted_file
            timestamp = None

        converted_file_path = os.path.join(EVAL_OUTPUT_DIR, converted_filename)
        
        status_html, tp_fp_path, metrics_path = run_evaluation(
            converted_file_path, reference_choice, caller_tool,timestamp,
            st_basic, l_basic, xl_basic, rl_basic, rxl_basic, mo_basic, of_basic,
            parent1_content, parent1, parent2_content, parent2, rb_content, rb,
            c, mr, rs, mins, eg, eb, i, y, sm
        )

        # Load metrics (same as your build_evalsvcallers_visuals)
        with open(metrics_path, "r") as f:
            encoded_contents = "data:text/plain;base64," + base64.b64encode(f.read().encode()).decode()
        pivot_ref, df_long, df_block = parse_evalsvcallers_file(encoded_contents)
        metrics_html = generate_evalsvcallers_visuals(pivot_ref, df_long, df_block)

        # Load the VCF dataframe
        df = load_vcf_dataframe(tp_fp_path, "evalsvcallers")
        unique_svtypes = sorted(df["SVTYPE"].dropna().unique())
        valid_chroms = [f"chr{i}" for i in range(1,23)] + ["chrX","chrY"]
        df["CHROM"] = df["CHROM"].astype(str).apply(lambda x: f"chr{x}" if not x.startswith("chr") else x)
        available_chroms = [c for c in valid_chroms if c in df["CHROM"].unique()]

        # Generate all visuals
        visuals_html = plot_vcf_data(df)
        sankey_html = plot_sankey(df)
        clustergram_html = plot_clustergram(df, available_chroms) if available_chroms else html.Div("⚠️ No chromosomes.")
        
        # Circos JSON
        json_path = os.path.join(EVAL_OUTPUT_DIR, f"{os.path.basename(tp_fp_path)}_circos.json")
        vcf_to_circos_json(tp_fp_path, json_path, "evalsvcallers")
        with open(json_path, "r") as f:
            circos_data = json.load(f)
        svtype_colors = {track["name"]: track.get("color", "#95a5a6") for track in circos_data.get("tracks", [])}
        circos_html = plot_circos("histogram", unique_svtypes, circos_data, svtype_colors)

        try:
            all_svtypes = df["SVTYPE"].dropna().unique().tolist()
            manhattan_html = plot_manhattan(df, all_svtypes, 6)
        except Exception as e:
            manhattan_html = html.Div(f"❌ Manhattan error: {str(e)}")        # Manhattan

        # Build final UI block, exactly like your working visual callback
        visuals_ui = html.Div([
            html.H3("📊 Automatic Metrics & Visuals", style={'marginTop': '20px'}),
            metrics_html,
            html.Br(),

            # Basic visuals
            html.Label("Filter by Variant Type:"),
            dcc.Dropdown(
                id="eval-svtype-filter",
                options=[{"label": "ALL", "value": "ALL"}] + [{"label": sv, "value": sv} for sv in unique_svtypes],
                value="ALL",
                clearable=False,
                style={"width": "300px"}
            ),
            html.Div([
                dbc.Card([
                    dbc.CardHeader(html.H4(f"VCF File Preview (ALL)")),
                    dbc.CardBody([
                        dash_table.DataTable(
                                 data=df.head(5).to_dict('records'),
                                 columns=[{"name": col, "id": col} for col in df.columns],
                                 style_table={'overflowX': 'auto'},
                                 style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'},
                                  style_header={'fontWeight': 'bold', 'backgroundColor': 'lightgrey', 'border': '1px solid black'},
                                 style_data={'border': '1px solid black'}
                            )
                        ])
                    ], className="mt-4")
                 ], id="eval-preview-card"),

            # Basic visualizations block
            html.H4("Basic Visualizations"),
            html.Div(visuals_html, id="eval-visuals-output"),
            html.Br(),
            sankey_html,
            html.Br(),
            # Clustergram
            html.H4("Clustergram Visualization"),
            html.Div([
                html.Label("Select Chromosomes:"),
                dcc.Dropdown(
                    id="eval-clustergram-selector",
                    options=[{"label": c, "value": c} for c in available_chroms],
                    value=available_chroms,
                    multi=True,
                    style={"width": "600px"}
                )
            ]),
            html.Div(clustergram_html, id="eval-clustergram-output"),
            html.Br(),

            # Circos
            html.H4("Circos Visualization"),
            html.Div([
                html.Label("Graph Type:"),
                dcc.Dropdown(
                    id="eval-graph-type-dropdown",
                    options=[{"label": "Histogram", "value": "histogram"}],
                    value="histogram",
                    clearable=False,
                    style={"width": "200px"}
                ),
                html.Br(),
                html.Label("Select SV Types:"),
                dcc.Dropdown(
                    id="eval-svtype-dropdown",
                    options=[{"label": sv, "value": sv} for sv in unique_svtypes],
                    value=unique_svtypes,
                    multi=True,
                    style={"width": "400px"}
                )
            ], style={"marginBottom": "20px"}),
            html.Div(circos_html, id="eval-circos-output"),
            html.Br(),

            # Manhattan
            html.H4("Manhattan Visualization"),
            html.Div([
                html.Label("Select SV Types for Manhattan Plot:"),
                dcc.Dropdown(
                    id="eval-manhattan-selector",
                    options=[{"label": sv, "value": sv} for sv in unique_svtypes],
                    value=unique_svtypes,
                    multi=True,
                    style={"width": "400px"}
                ),
                html.Label("Manhattan Threshold:"),
                dcc.Slider(
                    id="eval-manhattan-slider",
                    min=0, max=10, step=1, value=6,
                    marks={i: str(i) for i in range(11)},
                    tooltip={"always_visible": False, "placement": "bottom"}
                )
            ], style={"marginBottom": "20px"}),
            html.Div(manhattan_html, id="eval-manhattan-output")
        ])

        return html.Div([status_html]), tp_fp_path, metrics_path, visuals_ui

    except Exception as e:
        return html.Div(f"❌ Evaluation failed: {str(e)}"), [], [], html.Div()
        
@app.callback(
    Output("reference_upload_section", "style"),
    Input("reference_choice", "value")
)
def toggle_reference_upload(reference_choice):
    if reference_choice == "custom":
        return {'width': '100%'}
    else:
        return {'display': 'none'}
        
   #     return html.Div([
   #         dcc.Upload(
   #             id={'type': 'dynamic-upload', 'index': 0},
   #             children=html.Div(['📎 Drag and Drop or ', html.A('Select VCF File')]),
   #             style={
   #                 'width': '400px',
   #                 'height': '60px',
   #                 'lineHeight': '60px',
   #                 'borderWidth': '1px',
   #                 'borderStyle': 'dashed',
   #                 'borderRadius': '5px',
   #                 'textAlign': 'center',
    #                'margin': '10px',
   #                 'fontFamily': '"Times New Roman", Times, serif'
    #            },
   #             multiple=False
   #         )
         #   html.Div(id={'type': 'dynamic-upload-status', 'index': 0})
    #    ])
   # return ""


@app.callback(
    Output('converted-upload-section', 'style'),
    Input('use_conversion_toggle', 'value')
)
def toggle_converted_upload_section(selected):
    return {'display': 'none'}

################### TRUVARI FUNCTIONS ########################
from dash import Output, Input
@app.callback(
    Output("tru-base-upload-info", "children"),
    Input("tru-base-upload", "filename"),
    Input("tru-base-upload", "contents"),
    prevent_initial_call=True
)
def show_base_upload_status(filename, contents):
    if filename and contents:
        return html.Div(f"✅ BASE file uploaded: {filename}", style={"color": "green"})
    else:
        return html.Div("⚠️ Waiting for BASE file upload...", style={"color": "gray"})

@app.callback(
    Output("tru-comp-upload-info", "children"),
    Input("tru-comp-upload", "filename"),
    Input("tru-comp-upload", "contents"),
    prevent_initial_call=True
)
def show_comp_upload_status(filename, contents):
    if filename and contents:
        return html.Div(f"✅ COMP file uploaded: {filename}", style={"color": "green"})
    else:
        return html.Div("⚠️ Waiting for COMP file upload...", style={"color": "gray"})



@app.callback(
    Output('tru-base-upload-container', 'style'),
    Input('base_choice', 'value')
)
def toggle_base_upload(choice):
    if choice == 'univar':
        return {'display': 'none'}
    else:
        return {'width': '100%'}
# ==================== TRUVARI: Run bench ====================
import gzip, shutil, os  # os zaten varsa tekrar ekleme

def _ensure_plain_vcf_for_circos(vcf_path: str) -> str:
    """
    visualization.vcf_to_circos_json() plain-text VCF beklediği için
    .vcf.gz dosyasını aynı klasöre .decompressed.vcf olarak açar ve o yolu döndürür.
    .gz değilse olduğu gibi döner.
    """
    if not vcf_path or not os.path.exists(vcf_path):
        return vcf_path
    if not vcf_path.endswith(".gz"):
        return vcf_path

    plain_path = vcf_path + ".decompressed.vcf"
    if not os.path.exists(plain_path):  # gereksiz tekrar yazmayı önle
        with gzip.open(vcf_path, "rb") as fin, open(plain_path, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    return plain_path
    
@app.callback(
    [Output("tru-status", "children"),
     Output("truvari-output-list", "children"),
     Output("truvari-metrics-preview", "children"),
     Output("truvari-visual-output", "children")
    ],  
    Input("tru-run-btn", "n_clicks"),
    State("base_choice", "value"),  # <-- yeni eklendi
    State("tru-base-upload", "filename"),
    State("tru-base-upload", "contents"),
    State("tru-comp-upload", "filename"),
    State("tru-comp-upload", "contents"),
    State("tru-ref-selector", "value"),
    State("tru-param-refdist", "value"),
    State("tru-param-pctsim", "value"),
    State("tru-param-pctsize", "value"),
    State("tru-param-pctovl", "value"),
    State("tru-param-sizemin", "value"),
    State("tru-param-sizefilt", "value"),
    State("tru-param-sizemax", "value"),
    State("tru-flag-typeignore", "value"),
    State("tru-flag-uselev", "value"),
    State("tru-flag-gtcomp", "value"),
    State("tru-flag-passonly", "value"),
    State("tru-flag-multimatch", "value"),
    prevent_initial_call=True
)
def run_truvari_bench(n_clicks,base_choice,
                      base_fname, base_contents,
                      comp_fname, comp_contents,
                      reference_value,
                      refdist, pctsim, pctsize, pctovl, sizemin, sizefilt, sizemax,
                      typeignore, uselev, gtcomp, passonly, multimatch):
    if not n_clicks:
        return dash.no_update, dash.no_update, dash.no_update, dash.no_update

    # --- Validate inputs
    # --- Handle base_vcf based on base_choice
    if base_choice == "univar":
        base_vcf_path = os.path.join(UPLOAD_DIRECTORY, "univar.vcf")
        if not os.path.exists(base_vcf_path):
            return html.Div("❌ Univar SV catalog (univar.vcf) not found."), html.Ul([html.Li("No output.")]), dash.no_update, dash.no_update
    elif base_choice == "custom":
        if not base_contents or not base_fname:
            return html.Div("❌ Please upload the BASE VCF."), html.Ul([html.Li("No output.")]), dash.no_update, dash.no_update
        try:
            base_vcf_path = save_uploaded_file(base_fname, base_contents)
        except Exception as e:
            return html.Div(f"❌ BASE file save error: {e}"), html.Ul([html.Li("No output.")]), dash.no_update, dash.no_update
    else:
        return html.Div("❌ Invalid base choice."), html.Ul([html.Li("No output.")]), dash.no_update, dash.no_update
        
    if not comp_contents or not comp_fname:
        return html.Div("❌ Please upload the COMPARISON VCF."), html.Ul([html.Li("No output.")]), dash.no_update, dash.no_update
    if not reference_value:
        return html.Div("❌ Please select a reference (.fa)."), html.Ul([html.Li("No output.")]), dash.no_update, dash.no_update

    # --- Save uploaded VCFs
    try:
        base_vcf_path = save_uploaded_file(base_fname, base_contents)
        comp_vcf_path = save_uploaded_file(comp_fname, comp_contents)
    except Exception as e:
        return html.Div(f"❌ File save error: {e}"), html.Ul([html.Li()]), dash.no_update, dash.no_update

    # --- Absolute reference path (layout stores relative 'uploaded_files/..')
    reference_fa = os.path.abspath(reference_value)

    # --- Build Truvari params
    def _maybe(v):  # pass only if user provided
        return v if v not in (None, "") else None

    params = {}
    # flags
    if bool(typeignore): params["--typeignore"] = True
    if bool(uselev):     params["--use-lev"]    = True
    if bool(gtcomp):     params["--gtcomp"]     = True
    if bool(passonly):   params["--passonly"]   = True
    if bool(multimatch): params["--multimatch"] = True
    # keyed values
    if _maybe(refdist)  is not None: params["--refdist"]  = refdist
    if _maybe(pctsim)   is not None: params["--pctsim"]   = pctsim
    if _maybe(pctsize)  is not None: params["--pctsize"]  = pctsize
    if _maybe(pctovl)   is not None: params["--pctovl"]   = pctovl
    if _maybe(sizemin)  is not None: params["--sizemin"]  = sizemin
    if _maybe(sizefilt) is not None: params["--sizefilt"] = sizefilt
    if _maybe(sizemax)  is not None: params["--sizemax"]  = sizemax

    # --- Run pipeline
    status_div, tru_paths = run_truvari_pipeline(
        base_choice,
        base_vcf_uploaded=base_vcf_path,
        comp_vcf=comp_vcf_path,
        reference_fa=reference_fa,
        params=params
    )
    # Pipeline başarısızsa
    if not isinstance(tru_paths, dict):
        return status_div, html.Ul([]), dash.no_update, dash.no_update

    # --- summary.json'u oku ve metrik görsellerini oluştur
    metrics_preview = dash.no_update
    try:
        summary_path = tru_paths.get("summary")
        if summary_path and os.path.exists(summary_path):
            with open(summary_path) as f:
                summary_data = json.load(f)
            # gt_matrix gibi nested yapıları dışla
            df_summary = pd.DataFrame(list(summary_data.items()), columns=["Metric", "Value"])
            df_summary = df_summary[~df_summary["Value"].apply(lambda x: isinstance(x, (dict, list)) or x is None)].copy()
            metrics_preview = generate_truvari_visuals(df_summary)
        else:
            metrics_preview = html.Div("❌ summary.json not found.", style={"color": "red"})
    except Exception as e:
        metrics_preview = html.Div(f"❌ Metrics render error: {e}", style={"color": "red"})

  # --- TRUVARI GÖRSELLERİ (tp-comp, tp-base, fp, fn)
    import plotly.express as px
    from dash import dcc
    
    def _load_and_tag(vcf_path, tag):
        if not vcf_path or not os.path.exists(vcf_path):
            return None
    
        df = load_vcf_dataframe(vcf_path, "caller")  # Truvari çıktı VCF'leri
        if df is None or df.empty:
            return None
    
        # ✅ CHROM normalize
        df["CHROM"] = df["CHROM"].astype(str).str.strip().str.replace("^chr", "", regex=True)
    
        # ✅ SVTYPE'i HER ZAMAN INFO'dan yeniden çıkar (ID/SMAP karışıklığını engeller)
        df["SVTYPE"] = df["INFO"].astype(str).str.extract(r'(?:^|;)SVTYPE=([^;]+)', expand=False)
    
        # (İstersen BND'leri tamamen dışla)
        # df = df[df["SVTYPE"] != "BND"]
    
        # ✅ SVLEN yoksa INFO'dan çıkar, sonra sayısala çevir
        if "SVLEN" not in df or df["SVLEN"].isna().all():
            df["SVLEN"] = df["INFO"].astype(str).str.extract(r'(?:^|;)SVLEN=(-?\d+)', expand=False)
        df["SVLEN"] = pd.to_numeric(df["SVLEN"], errors="coerce")
    
        # ✅ QUAL sayısal
        df["QUAL"] = pd.to_numeric(df.get("QUAL"), errors="coerce")
    
        # Kaynak etiketle
        df["__SOURCE__"] = tag
    
        # SVTYPE olmayanları at (x ekseninde yanlışları temizler)
        df = df.dropna(subset=["SVTYPE"])
    
        return df
        
    def _grouped_basic_section(title, pairs):
        """
        pairs: list of tuples (vcf_path, 'LABEL'), örn:
          [ (tp_base_path, 'TP-BASE'), (tp_comp_path, 'TP-COMP') ]
          [ (fp_path, 'FP'), (fn_path, 'FN') ]
        """
        dfs = []
        for p, lab in pairs:
            d = _load_and_tag(p, lab)
            if d is not None and not d.empty:
                dfs.append(d)
        if not dfs:
            return html.Div(f"⚠️ {title}: no data.")
    
        data = pd.concat(dfs, ignore_index=True)
    
        # 1) SVTYPE dağılımı: grouped bar
        sv_df = (
            data.dropna(subset=["SVTYPE"])
                .groupby(["__SOURCE__", "SVTYPE"])
                .size()
                .reset_index(name="Count")
        )
        fig_svtype = px.bar(
            sv_df, x="SVTYPE", y="Count", color="__SOURCE__",  color_discrete_sequence=px.colors.qualitative.Set2,barmode="group",
            title=f"{title} — SVTYPE Distribution (Grouped)"
        )
        fig_svtype.update_layout(margin=dict(l=20, r=20, t=60, b=30), height=420)
    
        # 2) Kromozom dağılımı: grouped bar (chr1..22, X, Y ile sınırlayalım)
        chrom_order = [str(i) for i in range(1,23)] + ["X","Y"]
        chrom_df = data[data["CHROM"].isin(chrom_order)].copy()
        chrom_df["CHROM"] = pd.Categorical(chrom_df["CHROM"], categories=chrom_order, ordered=True)
        chrom_df = (
            chrom_df.groupby(["__SOURCE__", "CHROM"])
                    .size()
                    .reset_index(name="Variant Count")
                    .sort_values(["CHROM", "__SOURCE__"])
        )
        fig_chr = px.bar(
            chrom_df, x="CHROM", y="Variant Count", color="__SOURCE__", color_discrete_sequence=px.colors.qualitative.Set2,barmode="group",
            title=f"{title} — Chromosome-wise Distribution (Grouped)"
        )
        fig_chr.update_layout(margin=dict(l=20, r=20, t=60, b=30), height=420)
    
        # 3) SVLEN dağılımı: histogram (kaynak renkle), violin (yan yana)
        svlen_df = data.dropna(subset=["SVTYPE", "SVLEN"]).copy()
        # histogram (log gösterim istersen log_y veya x ekseninde transform ekleyebiliriz)
        fig_svlen = px.histogram(
            svlen_df, x="SVLEN", color="__SOURCE__", color_discrete_sequence=px.colors.qualitative.Set2,nbins=40,
            title=f"{title} — SVLEN Histogram (by Source)"
        )
        fig_svlen.update_layout(margin=dict(l=20, r=20, t=60, b=30), height=420)
    
        # violin: SVTYPE’a göre, source yan yana
        violin_df = svlen_df.copy()
        # log kullanmak istersen:
        # violin_df["logSVLEN"] = np.log10(violin_df["SVLEN"].abs() + 1)
        fig_violin = px.violin(
            violin_df, x="SVTYPE", y="SVLEN", color="__SOURCE__", color_discrete_sequence=px.colors.qualitative.Set2,box=True, points="all",
            title=f"{title} — SVLEN by SVTYPE (Grouped Violin)"
        )
        fig_violin.update_layout(margin=dict(l=20, r=20, t=60, b=30), height=480)


        all_types = sorted(sv_df["SVTYPE"].unique().tolist())
        preferred = ["DEL", "INS", "DUP", "INV", "BND", "CNV"]
        ordered_types = [t for t in preferred if t in all_types] + [t for t in all_types if t not in preferred]
        
        fig_spider = go.Figure()
        if ordered_types:  # veri varsa çiz
            # Kaynaklara (TP-BASE, TP-COMP, FP, FN) göre pivotla
            pivot = (
                sv_df.pivot_table(index="SVTYPE", columns="__SOURCE__", values="Count", aggfunc="sum")
                    .reindex(index=ordered_types)
                    .fillna(0)
                    .astype(int)
            )
            thetas = ordered_types + [ordered_types[0]]  # poligonu kapatmak için ilkini sona ekle
            palette = px.colors.qualitative.Set2
        
            for idx, lab in enumerate(sv_df["__SOURCE__"].drop_duplicates().tolist()):
                r_vals = pivot[lab].tolist() if lab in pivot.columns else [0] * len(ordered_types)
                r_closed = r_vals + [r_vals[0]]
                fig_spider.add_trace(go.Scatterpolar(
                    r=r_closed,
                    theta=thetas,
                    fill='toself',
                    name=lab,
                    line=dict(width=2, color=palette[idx % len(palette)]),
                    opacity=0.6
                ))
        
        fig_spider.update_layout(
            polar=dict(radialaxis=dict(visible=True)),
            title=f"{title} — Variant Type Radar Chart",
            showlegend=True,
            margin=dict(l=20, r=20, t=60, b=30),
            height=480
                )
        
        return html.Div([
            html.H4(title, style={'marginTop': '15px'}),
            dcc.Graph(figure=fig_svtype),
            dcc.Graph(figure=fig_spider),   # ← EK
            dcc.Graph(figure=fig_chr),
            dcc.Graph(figure=fig_svlen),
            dcc.Graph(figure=fig_violin)
        ])



    def advanced_truvari_section(title, pairs):
        import plotly.express as px 
        from dash import html
        import os
        dfs = []
        for p, lab in pairs:
            d = _load_and_tag(p, lab)
            if d is not None and not d.empty:
                dfs.append(d)
        if not dfs:
            return html.Div(f"⚠️ {title}: no data.")
        
        data = pd.concat(dfs, ignore_index=True)
        # --- SANKEY: TP-BASE ➝ SVTYPE ➝ TP-COMP ---
        sankey_df = data.dropna(subset=["CHROM", "SVTYPE", "__SOURCE__"]).copy()
        sankey_df["CHROM"] = "chr" + sankey_df["CHROM"].astype(str)
        
        # CHROM etiketlerini BASE ve COMP'e göre ayır
        sankey_df["CHROM_SRC"] = sankey_df.apply(
            lambda row: f"{row['CHROM']} (BASE)" if row["__SOURCE__"] == "TP-BASE" else f"{row['CHROM']} (COMP)", axis=1
        )
        
        # BASE için: CHROM_SRC → SVTYPE
        base_df = sankey_df[sankey_df["__SOURCE__"] == "TP-BASE"].copy()
        base_df_grouped = base_df.groupby(["CHROM_SRC", "SVTYPE"]).size().reset_index(name="Count")
        
        # COMP için: SVTYPE → CHROM_SRC
        comp_df = sankey_df[sankey_df["__SOURCE__"] == "TP-COMP"].copy()
        comp_df_grouped = comp_df.groupby(["SVTYPE", "CHROM_SRC"]).size().reset_index(name="Count")
        
        # Etiketler
        base_chroms = sorted(base_df_grouped["CHROM_SRC"].unique())
        svtypes = sorted(set(base_df_grouped["SVTYPE"]).union(set(comp_df_grouped["SVTYPE"])))
        comp_chroms = sorted(comp_df_grouped["CHROM_SRC"].unique())
        all_labels = base_chroms + svtypes + comp_chroms
        label_map = {label: i for i, label in enumerate(all_labels)}
        
        # BASE → SVTYPE bağlantıları
        base_source = base_df_grouped["CHROM_SRC"].map(label_map)
        base_target = base_df_grouped["SVTYPE"].map(label_map)
        base_value = base_df_grouped["Count"]
        
        # SVTYPE → COMP bağlantıları
        comp_source = comp_df_grouped["SVTYPE"].map(label_map)
        comp_target = comp_df_grouped["CHROM_SRC"].map(label_map)
        comp_value = comp_df_grouped["Count"]
        
        # SVTYPE'a özel renkler
        svtype_color_map = {sv: px.colors.qualitative.Set2[i % 10] for i, sv in enumerate(svtypes)}
        base_color = base_df_grouped["SVTYPE"].map(svtype_color_map)
        comp_color = comp_df_grouped["SVTYPE"].map(svtype_color_map)
        
        # Sankey oluştur
        sankey_fig = go.Figure(go.Sankey(
            arrangement="snap",
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color="black", width=0.5),
                label=all_labels,
                color=["#DC143C"] * len(all_labels)
            ),
            link=dict(
                source=pd.concat([base_source, comp_source]),
                target=pd.concat([base_target, comp_target]),
                value=pd.concat([base_value, comp_value]),
                color=pd.concat([base_color, comp_color])
            )
        ))
        
        sankey_fig.update_layout(
            title=f"{title} — TP-BASE → SVTYPE → TP-COMP Sankey",
            font_size=12,
            height=600,
            margin=dict(l=20, r=20, t=60, b=30)
        )

        # ---------- CLUSTERGRAM ----------
        clust_df = data.copy()
        
        # Sadece TP-BASE / TP-COMP kullan (FP/FN vs. gelirse karışmasın)
        clust_df = clust_df[clust_df["__SOURCE__"].isin(["TP-BASE", "TP-COMP"])].copy()
        
        # CHROM "chr" prefiksi (görsel tutarlılık)
        clust_df["CHROM"] = clust_df["CHROM"].astype(str).str.strip()
        clust_df["CHROM"] = clust_df["CHROM"].apply(lambda x: x if x.startswith("chr") else f"chr{x}")
        
        # Çoklu sütun: (__SOURCE__, SVTYPE)
        pivot = clust_df.pivot_table(
            index="CHROM",
            columns=["__SOURCE__", "SVTYPE"],
            aggfunc="size",
            fill_value=0
        )
        
        # Kromozom sıralaması (chr1..22, chrX, chrY) + sadece mevcut olanları al
        chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        pivot = pivot.reindex([c for c in chrom_order if c in pivot.index])
        
        # Bazı SVTYPE/SOURCE kombinasyonları yoksa kolonları tamlamak isteyebilirsin (opsiyonel):
        # desired_sv = ["DEL", "INS", "DUP", "INV", "BND", "CNV"]
        # desired_sources = ["TP-BASE", "TP-COMP"]
        # import itertools
        # full_cols = pd.MultiIndex.from_product([desired_sources, desired_sv], names=["__SOURCE__", "SVTYPE"])
        # pivot = pivot.reindex(columns=full_cols, fill_value=0)
        
        # Kolon etiketlerini okunur yap: "TP-BASE:DEL" gibi düz string'e çevir
        if isinstance(pivot.columns, pd.MultiIndex):
            nice_cols = [f"{src}:{sv}" for (src, sv) in pivot.columns]
        else:
            # Çok nadir bir durumda MultiIndex düşerse
            nice_cols = pivot.columns.astype(str).tolist()
        
        clustergram_fig = Clustergram(
            data=pivot.values,
            row_labels=pivot.index.tolist(),     # chr1..chrY
            column_labels=nice_cols,             # "TP-BASE:DEL", "TP-COMP:DEL", ...
            color_map="RdBu",
            height=800,
            width=1000,
            display_ratio=[0.85, 0.15],
            # hidden_labels=[],
            color_threshold={"row": 0.5, "col": 0.5}
        )
        
        # --- Tüm x/y eksenlerinde tick label'ları aç ---
        layout_dict = clustergram_fig.layout.to_plotly_json()
        axis_updates = {}
        for key in layout_dict.keys():
            if str(key).startswith("xaxis") or str(key).startswith("yaxis"):
                axis_updates[key] = dict(showticklabels=True)
        
        clustergram_fig.update_layout(
            **axis_updates,
            margin=dict(l=110, r=60, t=50, b=110)  # etiketlere alan
        )
        
    # ------------------ CIRCOS (TP-BASE & TP-COMP) ENTEGRASYONU --------------------
        # TP-BASE / TP-COMP VCF yollarını 'pairs'ten çek
        # --- Genome length dictionaries ---
        # --- Genome lengths ---
        hg19_lengths = {
            "chr1": 249250621, "chr2": 243199373, "chr3": 198022430,
            "chr4": 191154276, "chr5": 180915260, "chr6": 171115067,
            "chr7": 159138663, "chr8": 146364022, "chr9": 141213431,
            "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
            "chr13": 115169878, "chr14": 107349540, "chr15": 102531392,
            "chr16": 90354753,  "chr17": 81195210,  "chr18": 78077248,
            "chr19": 59128983,  "chr20": 63025520,  "chr21": 48129895,
            "chr22": 51304566,  "chrX": 155270560,  "chrY": 59373566
        }
        hg38_lengths = {
            "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
            "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
            "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
            "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
            "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
            "chr16": 90338345,  "chr17": 83257441,  "chr18": 80373285,
            "chr19": 58617616,  "chr20": 64444167,  "chr21": 46709983,
            "chr22": 50818468,  "chrX": 156040895,  "chrY": 57227415
        }

        # --- Chromosome colors and order ---
        chrom_colors = [
            "#f28e2b", "#4e79a7", "#e15759", "#76b7b2", "#59a14f", "#edc948",
            "#b07aa1", "#ff9da7", "#9c755f", "#bab0ab", "#8cd17d", "#b6992d",
            "#499894", "#d37295", "#fabfd2", "#d4a6c8", "#9f9f9f", "#bcbd22",
            "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#c5b0d5", "#c49c94"
        ]
        chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        chrom_color_map = dict(zip(chrom_order, chrom_colors))
                       
        # ------------------ CIRCOS (TP-BASE & TP-COMP) ENTEGRASYONU --------------------
        tp_base_vcf = next((p for p, lab in pairs if lab == "TP-BASE"), None)
        tp_comp_vcf = next((p for p, lab in pairs if lab == "TP-COMP"), None)
        
        circos_div = html.Div("⚠️ Circos: TP-BASE / TP-COMP not found.")
        if tp_base_vcf and tp_comp_vcf and os.path.exists(tp_base_vcf) and os.path.exists(tp_comp_vcf):
            import tempfile, uuid, re, json, math
            import plotly.express as px
            import time
            from pathlib import Path
    #        tmp_json = os.path.join(tempfile.gettempdir(), f"circos_truvari_{uuid.uuid4().hex}.json")
            run_dir = Path(summary_path).resolve().parents[0]

            ts = time.strftime("%Y%m%d_%H%M%S")
            base_tag = Path(tp_base_vcf).stem
            comp_tag = Path(tp_comp_vcf).stem
         #   dest_json = os.path.join(TRUVARI_OUTPUT_DIR, f"circos_{ts}_{base_tag}_vs_{comp_tag}.json")
            dest_json = run_dir / f"circos_{ts}_{base_tag}_vs_{comp_tag}.json"   # ← DOSYA YOLU!

            # 1) JSON üret
            vcf_to_circos_json_truvari(
                tp_base_vcf, tp_comp_vcf,
                dest_json,
                svtype_color_map=svtype_color_map,   # layout/ideogram için kullanılabilir
                chrom_color_map=chrom_color_map,
                hg38_lengths=hg38_lengths,
                hg19_lengths=hg19_lengths,
                bin_size=10_000_000
            )
        
            # 2) JSON oku
            with open(dest_json, "r", encoding="utf-8") as f:
                circos_data = json.load(f)
                
            saved_note = html.Div([
                html.Span("💾 Circos JSON kaydedildi: "),
                html.Code(dest_json, style={"fontSize": "0.9em"})
            ], style={"margin": "6px 0 12px"})

            layout = circos_data.get("layout", [])
            if not layout:
                circos_div = html.Div("⚠️ Circos: invalid layout.")
            else:
                layout_ids = {str(b.get("id")) for b in layout}
                raw_tracks = circos_data.get("tracks", [])
                all_tracks = [t for t in raw_tracks if isinstance(t, dict) and isinstance(t.get("name"), str)]
                if not all_tracks:
                    circos_div = html.Div("⚠️ Circos: no tracks in JSON.")
                else:
                    # --- Yardımcılar
                    def parse_name(name: str):
                        # "DEL (TP-BASE)" / "INS (TP-COMP)" gibi adları yakala
                        m = re.match(r"\s*([A-Za-z0-9_+-]+)\s*\(\s*(TP-BASE|TP-COMP)\s*\)\s*$", name, flags=re.I)
                        if m:
                            return m.group(1).upper(), m.group(2).upper()
                        low = name.lower()
                        grp = "TP-BASE" if "tp-base" in low else ("TP-COMP" if "tp-comp" in low else None)
                        # SVTYPE fallback
                        if ":" in name:
                            sv = name.split(":")[-1].strip().upper()
                        else:
                            sv = name.split("(")[0].split("-")[0].strip().upper()
                        return sv, grp
        
                    def _colorize_bins(bins, col):
                        # her item’a da renk yaz (sürüm uyumluluğu)
                        return [{**it, "color": col} for it in bins]
        
                    # 3) HISTOGRAM verilerini topla (sv, grp, data)
                    hist_elems = []  # list[(sv, grp, data)]
                    seen_sv_order, seen_sv = [], set()
                    for t in all_tracks:
                        if str(t.get("type", "")).lower() != "histogram":
                            continue
                        sv, grp = parse_name(t.get("name", ""))
                        if grp not in ("TP-BASE", "TP-COMP") or not sv:
                            continue
        
                        clean = []
                        for item in t.get("data", []):
                            blk = str(item.get("block_id"))
                            if blk not in layout_ids:
                                continue
                            try:
                                s = int(item.get("start", 0)); e = int(item.get("end", 0)); v = int(item.get("value", 0))
                            except Exception:
                                continue
                            if e <= s or v < 0:
                                continue
                            clean.append({"block_id": blk, "start": s, "end": e, "value": v})
                        if not clean:
                            continue
        
                        hist_elems.append((sv, grp, clean))
                        if sv not in seen_sv:
                            seen_sv.add(sv); seen_sv_order.append(sv)
        
                    if not hist_elems:
                        circos_div = html.Div("⚠️ Circos: no histogram data.")
                    else:
                        # 4) SVTYPE bazında BASE/COMP var mı yok mu (sadece varsa yerleştirilecek)
                        base_bins_by_sv = {sv: [] for sv in seen_sv_order}
                        comp_bins_by_sv = {sv: [] for sv in seen_sv_order}
                        for sv, grp, data in hist_elems:
                            if grp == "TP-BASE":
                                base_bins_by_sv[sv].extend(data)
                            elif grp == "TP-COMP":
                                comp_bins_by_sv[sv].extend(data)
        
                        # 5) KOMBO-BAZLI RENK: her (grp, sv) kombinasyonu için benzersiz renk
                        #    - önce geniş bir palette döndür
                        #    - palet yetmezse HSL hash ile HEX üret
                        big_palette = (
                            px.colors.qualitative.Vivid
                            + px.colors.qualitative.Set2
                            + px.colors.qualitative.Plotly
                            + px.colors.qualitative.Safe
                            + px.colors.qualitative.Bold
                            + px.colors.qualitative.Pastel
                        )
                        combo_keys = []  # sıralı
                        for grp in ("TP-BASE", "TP-COMP"):
                            for sv in seen_sv_order:
                                if grp == "TP-BASE" and base_bins_by_sv[sv]:
                                    combo_keys.append(f"{grp}:{sv}")
                                if grp == "TP-COMP" and comp_bins_by_sv[sv]:
                                    combo_keys.append(f"{grp}:{sv}")
        
                        def hsl_hash_hex(key: str):
                            # stabil, geniş renk alanı için hash tabanlı HSL → HEX
                            h = abs(hash(key))
                            hue = (h % 360)            # 0..359
                            sat = 65 + (h % 25)        # 65..89
                            lig = 45 + (h % 10)        # 45..54
                            # HSL -> RGB
                            c = (1 - abs(2*lig/100 - 1)) * (sat/100)
                            x = c * (1 - abs(((hue/60) % 2) - 1))
                            m = lig/100 - c/2
                            if   0 <= hue < 60:   r,g,b = c,x,0
                            elif 60 <= hue <120:  r,g,b = x,c,0
                            elif 120<= hue <180:  r,g,b = 0,c,x
                            elif 180<= hue <240:  r,g,b = 0,x,c
                            elif 240<= hue <300:  r,g,b = x,0,c
                            else:                 r,g,b = c,0,x
                            R = int((r+m)*255); G = int((g+m)*255); B = int((b+m)*255)
                            return f"#{R:02x}{G:02x}{B:02x}"
        
                        combo_color_map = {}
                        for i, key in enumerate(combo_keys):
                            if i < len(big_palette):
                                combo_color_map[key] = big_palette[i]
                            else:
                                combo_color_map[key] = hsl_hash_hex(key)
        
                        # 6) TRACKLER: yalnızca ilgili SVTYPE o grupta varsa ekle; renk = combo_color_map
                        tracks = []
        
                        # İç halka: TP-BASE
                        for sv in seen_sv_order:
                            bins = base_bins_by_sv.get(sv, [])
                            if not bins:
                                continue
                            key = f"TP-BASE:{sv}"
                            col = combo_color_map[key]
                            bins_colored = _colorize_bins(bins, col)
                            tracks.append({
                                "type": "HISTOGRAM",
                                "name": key,
                                "data": bins_colored,
                                "innerRadius": 160,
                                "outerRadius": 195,
                                "color": col,
                                "config": {
                                    "innerRadius": 160,
                                    "outerRadius": 195,
                                    "opacity": 1.0,
                                    "strokeWidth": 0,
                                    "color": col
                                }
                            })
        
                        # Dış halka: TP-COMP
                        for sv in seen_sv_order:
                            bins = comp_bins_by_sv.get(sv, [])
                            if not bins:
                                continue
                            key = f"TP-COMP:{sv}"
                            col = combo_color_map[key]
                            bins_colored = _colorize_bins(bins, col)
                            tracks.append({
                                "type": "HISTOGRAM",
                                "name": key,
                                "data": bins_colored,
                                "innerRadius": 200,
                                "outerRadius": 235,
                                "color": col,
                                "config": {
                                    "innerRadius": 200,
                                    "outerRadius": 235,
                                    "opacity": 1.0,
                                    "strokeWidth": 0,
                                    "color": col
                                }
                            })
        
                        if not tracks:
                            circos_div = html.Div("⚠️ Circos: no histogram data after filtering.")
                        else:
                            # 7) LEGEND: sadece var olan kombinasyonlar ve birebir renkleri
                            legend_items = []
                            for key in combo_keys:
                                legend_items.append(
                                    html.Div([
                                        html.Span(style={
                                            'backgroundColor': combo_color_map[key],
                                            'display': 'inline-block','width': 12,'height': 12,
                                            'marginRight': 6,'borderRadius': 2
                                        }),
                                        html.Span(key)
                                    ], style={'display': 'inline-block','marginRight': 14,'marginBottom': 6})
                                )
                            legend = html.Div(legend_items, style={'padding': '6px 0'})
        
                            # 8) Render
                            circos_div = html.Div([
                                html.Hr(),
                                html.H4("Circos — TP-BASE & TP-COMP by SVTYPE (unique colors per combo)"),
                                legend,
                                dashbio.Circos(
                                    id="circos-truvari",
                                    layout=layout,
                                    tracks=tracks,
                                    config={
                                        "innerRadius": 260,
                                        "outerRadius": 320,
                                        "ticks": {"display": True, "spacing": 10_000_000, "labelSpacing": 5},
                                        "labelLayout": {"spacing": 60, "radialOffset": 90}
                                    },
                                    size=800
                                )
                            ])
        manhattan_div = plot_manhattan_truvari(
            pairs,                # [(tp_base, "TP-BASE"), (tp_comp, "TP-COMP")]
            selected_svtypes=None,  # istersen bir UI state geçir
            threshold=6
        )

        # -------------------------------------------------------------------------------

        return html.Div([
            html.H4(f"{title} — Advanced Graphs", style={'marginTop': '15px'}),
            dcc.Graph(figure=sankey_fig),
            dcc.Graph(figure=clustergram_fig),
            circos_div,  # Circos block
            manhattan_div,
        ])
    # ✅ TP-BASE vs TP-COMP (yan yana)
    tp_block_basic = _grouped_basic_section(
        "TP (TP-BASE vs TP-COMP)",
        [(tru_paths.get("tp_base"), "TP-BASE"), (tru_paths.get("tp_comp"), "TP-COMP")]
    )
    # ✅ Advanced TP grafikleri (sadece TP için)
    tp_block_advanced = advanced_truvari_section(
        "TP (TP-BASE vs TP-COMP)",
        [(tru_paths.get("tp_base"), "TP-BASE"), (tru_paths.get("tp_comp"), "TP-COMP")]
    )    
    # ❌ FP vs FN (yan yana)
    err_block = _grouped_basic_section(
        "FP vs FN",
        [(tru_paths.get("fp"), "FP"), (tru_paths.get("fn"), "FN")]
    )

    visuals_block = html.Div([
        tp_block_basic,
        html.Hr(),
        tp_block_advanced,
        html.Hr(),
        err_block  # FP/FN basic only
    ])    
    
    # Dosya listesi istenmiyor -> boş
    return status_div, html.Ul([]), metrics_preview, visuals_block

################### TRUVARI FUNCTIONS END ########################

################### VISUALIZATION FUNCTIONS ########################
@app.callback(
    Output('manhattan-threshold-slider', 'max'),
    Output('manhattan-threshold-slider', 'marks'),
    [Input('filtered-data', 'data')],
    [State('selected-input-source', 'data'),
     State('main-tabs', 'active_tab')])
def update_threshold_slider_range(json_data, source_type, active_tab):
    if active_tab != 'tab-visualization':
        raise dash.exceptions.PreventUpdate

    if not json_data:
        return 10, {i: str(i) for i in range(11)}

    try:
        df = pd.read_json(io.StringIO(json_data))
        if df is None or "QUAL" not in df.columns:
            return 10, {i: str(i) for i in range(11)}
        max_threshold = min(30, df["QUAL"].max() / 200)
        return max_threshold, {i: str(i) for i in range(0, int(max_threshold) + 1)}
    except Exception:
        return 10, {i: str(i) for i in range(11)}


@app.callback(
    Output("circos-container", "children", allow_duplicate=True),
    Output("default-circos-output", "children", allow_duplicate=True),
    Input("viz_selector", "value"),
    prevent_initial_call=True
)
def clear_circos_outputs(viz_selected):
    if "circos" not in viz_selected:
        return "", ""
    raise dash.exceptions.PreventUpdate

@app.callback(
    Output('manhattan-controls-container', 'style'),
    Input('viz_selector', 'value')
)
def toggle_manhattan_controls(selected):
    if selected == "manhattan":
        return {"display": "block", "marginTop": "20px"}
    return {"display": "none"}
@app.callback(
    Output('manhattan-svtype-selector', 'options'),
    [Input('filtered-data', 'data')],
    [State('main-tabs', 'active_tab')]
)
def update_manhattan_svtypes(json_data, active_tab):
    if active_tab != 'tab-visualization':
        raise dash.exceptions.PreventUpdate

    if not json_data:
        return []

    df = pd.read_json(io.StringIO(json_data))
    if "SVTYPE" not in df.columns:
        return []

    svtypes = df["SVTYPE"].dropna().unique().tolist()
    return [{"label": sv, "value": sv} for sv in svtypes]
    
@app.callback(
    Output('visualize-status', 'children'),
    Input('visualize-button', 'n_clicks'),
    prevent_initial_call=True
)
def show_status(n_clicks):
    return get_variant_extraction_status()
    
@app.callback(
    Output('selected-input-source', 'data'),
    Input('visualize-input-source', 'value')
)
def store_selected_source(value):
    return value

@app.callback(
    [Output('vcf-files-list', 'children'),
     Output('uploaded-file-path', 'data')],
    Input('upload-data', 'contents'),
    State('upload-data', 'filename'),
    State('upload-data', 'last_modified')
)
def handle_uploaded_files(list_of_contents, list_of_names, list_of_dates):
    if list_of_contents is None:
        return html.P("No file uploaded yet.", style={"color": "gray", 'fontFamily': '"Times New Roman", Times, serif'}), ""

    if isinstance(list_of_contents, str):
        list_of_contents = [list_of_contents]
        list_of_names = [list_of_names]
        list_of_dates = [list_of_dates]

    file_path = ""
    children = []
    for contents, filename, date in zip(list_of_contents, list_of_names, list_of_dates):
        item, file_path = parse_uploaded_vcf(contents, filename, date)
        children.append(item)

    return html.Ul(children), file_path

@app.callback(
    [Output('variant-type-filter', 'options'),
     Output('filtered-data', 'data')],
    Input('visualize-button', 'n_clicks'),
    State('uploaded-file-path', 'data'),
    State('selected-input-source', 'data'),
    prevent_initial_call=True
)
def extract_variant_types_callback(n_clicks, file_path, source_type):
    return extract_variant_types(file_path, source_type)

@app.callback(
    Output("visual_output", "children"),
    Input("viz_selector", "value"),
    Input('uploaded-file-path', 'data'),
    Input('variant-type-filter', 'value'),           # ✅ INPUT: to re-trigger on change
    Input('filtered-data', 'data'),                  # ✅ INPUT: filtered JSON data
    Input('manhattan-svtype-selector', 'value'),       # ✅ NEW
    Input('manhattan-threshold-slider', 'value'),
    State('selected-input-source', 'data'),
    prevent_initial_call=True
)
def render_selected_visuals(selected, file_path, selected_svtype, json_data, svtypes_manhattan, threshold_value, source_type):
    visuals = []

    # ✅ Load filtered data from JSON
    if not json_data:
        return [html.Div("")]

    df = pd.read_json(io.StringIO(json_data))
    if df is None or df.empty:
        return [html.Div("⚠️ No valid structural variant data to display.")]

    # ✅ BASIC VISUALS SECTION
    if selected == "basic":
        try:
            df_basic = df.copy()
            if selected_svtype and selected_svtype != "ALL":
                df_basic = df_basic[df_basic["SVTYPE"] == selected_svtype]

            figures = plot_vcf_data(df_basic)

            visuals.append(
                dbc.Card([
                    dbc.CardHeader(html.H4(f"Filtered VCF File Preview ({selected_svtype})")),
                    dbc.CardBody([
                        dcc.Loading(
                            dash.dash_table.DataTable(
                                data=df_basic.head(5).to_dict('records'),
                                columns=[{"name": col, "id": col} for col in df_basic.columns],
                                style_table={'overflowX': 'auto'},
                                style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'},
                                style_header={'fontWeight': 'bold', 'backgroundColor': 'lightgrey', 'border': '1px solid black'},
                                style_data={'border': '1px solid black'}
                            )
                        )
                    ])
                ], className="mt-4")
            )

            visuals.append(
                dbc.Card([
                    dbc.CardHeader(html.H4("Filtered Visualizations")),
                    dbc.CardBody(figures)
                ], className="mt-4")
            )
            visuals.append(html.Hr(style={"marginTop": "30px", "marginBottom": "30px"}))


        except Exception as e:
            visuals.append(html.Div(f"❌ Basic visuals error: {e}"))
        
    # ✅ Sankey
    if selected == "sankey":
        try:
            visuals.append(plot_sankey(df))
            visuals.append(html.Hr(style={"marginTop": "30px", "marginBottom": "30px"}))

        except Exception as e:
            visuals.append(html.Div(f"❌ Sankey error: {e}"))

    # ✅ Clustergram
    if selected == "clustergram":
        try:
            df_cluster = df.copy()
            df_cluster["CHROM"] = df_cluster["CHROM"].astype(str).str.strip()
            valid_chroms = [str(i) for i in range(1, 23)] + ["X", "Y"]
            df_cluster["CHROM"] = df_cluster["CHROM"].apply(
                lambda x: x if x.startswith("chr") and x[3:] in valid_chroms
                else f"chr{x}" if x in valid_chroms
                else None
            )
            df_cluster = df_cluster.dropna(subset=["CHROM"])
    
            pivot_df = df_cluster.groupby(["CHROM", "SVTYPE"]).size().unstack(fill_value=0)
            chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
            available_chroms = [c for c in chrom_order if c in pivot_df.index]


            if not available_chroms:
                visuals.append(html.Div("⚠️ No available chromosomes for clustergram."))
            else:
                visuals.append(html.Div([
                    html.Label("Rows to display (chromosomes)", style={"marginTop": "10px"}),
                    dcc.Dropdown(
                        id='clustergram-chromosome-selector',
                        options=[{"label": c, "value": c} for c in available_chroms],
                        value=available_chroms,
                        multi=True,
                        style={"width": "600px"}
                    ),
                    html.Div(id='clustergram-display-area')
                ]))
                visuals.append(html.Hr(style={"marginTop": "30px", "marginBottom": "30px"}))

        except Exception as e:
            import traceback
            error_details = traceback.format_exc()
            visuals.append(html.Div([
                html.Div("❌ Clustergram error:"),
                html.Pre(error_details)
            ], style={"color": "red"}))
                
    # ✅ Manhattan (you must get 'svtypes' and 'threshold' from UI inputs later)
    if selected == "manhattan":
        try:
            svtypes_to_plot = svtypes_manhattan or df["SVTYPE"].dropna().unique().tolist()

    
            if not svtypes_to_plot:
                svtypes_to_plot = df["SVTYPE"].dropna().unique().tolist()
            if "SVLEN" not in df.columns or df["SVLEN"].isnull().all():
                df["SVLEN"] = df["INFO"].apply(
                    lambda x: int(re.search(r"SVLEN=-?\d+", x).group().split("=")[1])
                    if pd.notnull(x) and isinstance(x, str) and re.search(r"SVLEN=-?\d+", x)
                    else None
                )
            visuals.append(plot_manhattan(df, svtypes_to_plot, threshold_value))
            visuals.append(html.Hr(style={"marginTop": "30px", "marginBottom": "30px"}))

        except Exception as e:
            visuals.append(html.Div(f"❌ Manhattan error: {e}"))

    # ✅ Circos
    if selected == "circos":
        try:
            if file_path.endswith(".vcf"):
                vcf_filename = os.path.basename(file_path)
                timestamp = time.strftime("%Y%m%d_%H%M%S")
                json_filename = f"{vcf_filename.replace('.vcf', '')}_{timestamp}_circos.json"
                json_path = os.path.join(VISUALIZATION_OUTPUT_DIR, json_filename)

                vcf_to_circos_json(file_path, json_path, source_type)
                visuals.append(html.Div(f"✅ Circos JSON saved: {os.path.basename(json_path)} (Full path: {json_path})"))

            elif file_path.endswith(".json"):
                json_path = file_path
                visuals.append(html.Div(f"✅ Using uploaded Circos JSON file: {os.path.basename(json_path)} (Full path: {json_path})"))
            else:
                visuals.append(html.Div("❌ Invalid file. Please upload a `.vcf` or `.json` file."))
                return visuals

            with open(json_path, "r") as f:
                circos_data = json.load(f)
                available_svtypes = [track["name"] for track in circos_data.get("tracks", [])]
                svtype_colors = {track["name"]: track.get("color", "#95a5a6") for track in circos_data.get("tracks", [])}

            visuals.extend([
                dcc.Store(id='circos-data-store', data=circos_data),
                dcc.Store(id='svtype-colors-store', data=svtype_colors),
                html.Div([
                    html.Label("Graph Type:"),
                    dcc.Dropdown(
                        id="graph-type-dropdown",
                        options=[#{"label": "Line", "value": "chords"},
                                 {"label": "Histogram", "value": "histogram"}],
                        value="histogram", clearable=False, style={"width": "200px"}
                    )
                ], style={"marginTop": "20px"}),
                html.Div([
                    html.Label("Select SV Types:"),
                    dcc.Dropdown(
                        id="svtype-filter-dropdown",
                        options=[{"label": sv, "value": sv} for sv in available_svtypes],
                        value=available_svtypes, multi=True, style={"width": "400px"}
                    )
                ], id="svtype-dropdown-container", style={"marginTop": "10px"})
                
         #       html.Div(id="circos-container", style={"marginTop": "20px"})
        #        html.Div(id="default-circos-output", style={"marginTop": "20px"})
            ])
            visuals.append(html.Hr(style={"marginTop": "30px", "marginBottom": "30px"}))


        except Exception as e:
            visuals.append(html.Div(f"❌ Circos error: {e}"))

    return visuals

@app.callback(
    Output("clustergram-display-area", "children"),
    Input("clustergram-chromosome-selector", "value"),
    State("uploaded-file-path", "data"),
    State("selected-input-source", "data")
)
def update_clustergram(selected_chroms, file_path, source_type):
    if not selected_chroms:
        return html.Div("⚠️ Please select at least one chromosome.")
    
    df = load_vcf_dataframe(file_path, source_type)
    if df is None or df.empty:
        return html.Div("❌ No valid data for Clustergram.")

    try:
        return plot_clustergram(df, selected_chroms)
    except Exception as e:
        return html.Div(f"❌ Failed to render Clustergram: {e}")
from dash.exceptions import PreventUpdate

@app.callback(
    Output("circos-container", "children"),
    Input("graph-type-dropdown", "value"),
    Input("svtype-filter-dropdown", "value"),
    State("circos-data-store", "data"),
    State("svtype-colors-store", "data"),
    prevent_initial_call=True
)
def update_circos_plot(graph_type, selected_svtypes, circos_data, svtype_colors):
    if not circos_data:
        raise PreventUpdate

    return plot_circos(graph_type, selected_svtypes, circos_data, svtype_colors)
@app.callback(
    Output("svtype-colors-store", "style"),
    Input("graph-type-dropdown", "value")
)
def toggle_svtype_dropdown(graph_type):
    return {"display": "block", "marginTop": "10px"} if graph_type == "histogram" else {"display": "none"}

# Hover feedback
@app.callback(
    Output("default-circos-output", "children"),
    Input("circos-plot", "eventDatum")
)
def show_event_data(value):
    if value:
        return [html.Div(f"{k.title()}: {v}") for k, v in value.items()]
    return "Hover over a data point to see more info."
################### METRICS FUNCTIONS ########################
@app.callback(
    [Output('output-metrics', 'children', allow_duplicate=True),
     Output('summary-visualization-output', 'children', allow_duplicate=True)],
    Input('metrics-input-source', 'value'),
    prevent_initial_call=True
)
def clear_other_output(comp_type):
    if comp_type == 'truvari':
        return dash.no_update, ""  # ✅ Keep SURVIVOR, clear Eval
    elif comp_type == 'evalsvcallers':
        return "", dash.no_update  # ✅ Clear SURVIVOR, keep Eval
    return "", ""


@app.callback(
    Output('metrics-upload-section', 'children'),
    Input('metrics-input-source', 'value')
)
def render_upload_section(source):
    if source == 'truvari':
        return get_truvari_upload_section()
    elif source == 'evalsvcallers':
        return get_evalsvcallers_upload_section()
    #return html.Div("Please select a tool.", style={'fontFamily': '"Times New Roman", Times, serif'})
@app.callback(
    Output('caller-upload-status', 'children'),
    Output('reference-upload-status', 'children'),
    Output('intersection-upload-status', 'children'),
    Input('upload-caller', 'contents'),
    Input('upload-reference', 'contents'),
    Input('upload-intersection', 'contents'),
    State('upload-caller', 'filename'),
    State('upload-reference', 'filename'),
    State('upload-intersection', 'filename'),
    prevent_initial_call=True
)
def show_survivor_upload_status(caller_contents, ref_contents, inter_contents,
                                 caller_name, ref_name, inter_name):
    def msg(contents, name):
        if contents and name:
            return html.Div(f"✅ '{name}' uploaded successfully.", style={'color': 'green','fontFamily': '\"Times New Roman\", Times, serif'})
        return html.Div("❌ No file uploaded.", style={'color': 'red','fontFamily': '\"Times New Roman\", Times, serif'})

    return (
        msg(caller_contents, caller_name),
        msg(ref_contents, ref_name),
        msg(inter_contents, inter_name)
    )    
#@app.callback(
#    Output('output-metrics', 'children'),
#    Input('process-button', 'n_clicks'),
#    State('upload-reference', 'contents'),
#    State('upload-caller', 'contents'),
#    State('upload-intersection', 'contents'),
#    State('reference-type', 'value'),
#    prevent_initial_call=True
#)
#def process_survivor_metrics(n_clicks, ref_content, caller_content, inter_content, reference_type):
#    if not (ref_content and caller_content and inter_content):
#        return html.Div("⚠️ Please upload a valid summary.txt file.", style={'fontFamily': '\"Times New Roman\", Times, serif'})

#    metrics_df, error_msg = process_survivor_metrics_from_content(ref_content, caller_content, inter_content, reference_type)

#    if error_msg:
#        return html.Div(error_msg, style={'color': 'red'})

#    return generate_survivor_visuals(metrics_df)


@app.callback(
    Output('output-metrics', 'children'),
    Input("process-tru-button", "n_clicks"),
    State("upload-tru-file", "contents"),
    State("upload-tru-file", "filename"),
    prevent_initial_call=True
)
def process_truvari_file(n_clicks, contents, filename):
    if not contents:
        return html.Div(
            "⚠️ Please upload a valid summary file (.json).",
            style={'fontFamily': '"Times New Roman", Times, serif'}
        )
    # Accept both .txt and .json; skip strict extension checks
    try:
        df = parse_truvari_file(contents)  # returns LONG df
        if df is None or df.empty:
            return html.Div(
                "❌ File is empty or not valid JSON.",
                style={'fontFamily': '"Times New Roman", Times, serif'}
            )
        return generate_truvari_visuals(df)
    except Exception as e:
        return html.Div(
            f"❌ Error processing file: {e}",
            style={'fontFamily': '"Times New Roman", Times, serif'}
        )

@app.callback(
    Output('upload-tru-txt-status', 'children'),
    Input('upload-tru-file', 'contents'),
    State('upload-tru-file', 'filename'),
    prevent_initial_call=True
)
def show_tru_file_upload_status(contents, filename):
    if contents and filename:
        return html.Div(f"✅ '{filename}' uploaded successfully.", style={'color': 'green','fontFamily': '\"Times New Roman\", Times, serif'})
    return html.Div("❌ Upload failed.", style={'color': 'red','fontFamily': '\"Times New Roman\", Times, serif'})
   
@app.callback(
    Output("summary-visualization-output", "children"),
    Input("process-eval-button", "n_clicks"),
    State("upload-eval-file", "contents"),
    State("upload-eval-file", "filename"),
    prevent_initial_call=True
)
def process_evalsvcallers_file(n_clicks, contents, filename):
    if not contents or not filename.endswith(".eval.txt"):
        return html.Div("⚠️ Please upload a valid .eval.txt file.",style={'fontFamily': '\"Times New Roman\", Times, serif'})

    try:
        pivot_ref, df_long, df_block = parse_evalsvcallers_file(contents)
        return generate_evalsvcallers_visuals(pivot_ref, df_long, df_block)
    except Exception as e:
        return html.Div(f"❌ Error processing file: {e}",style={'fontFamily': '\"Times New Roman\", Times, serif'})

@app.callback(
    Output('upload-eval-txt-status', 'children'),
    Input('upload-eval-file', 'contents'),
    State('upload-eval-file', 'filename'),
    prevent_initial_call=True
)
def show_eval_file_upload_status(contents, filename):
    if contents and filename:
        return html.Div(f"✅ '{filename}' uploaded successfully.", style={'color': 'green','fontFamily': '\"Times New Roman\", Times, serif'})
    return html.Div("❌ Upload failed.", style={'color': 'red','fontFamily': '\"Times New Roman\", Times, serif'})

################### DOWNLOAD ROUTES ########################

@app.server.route('/download/survivor/<path:filename>')
def download_survivor_file(filename):
    survivor_dir = os.path.join(app.server.root_path, 'uploaded_files', 'survivor_output')
    return flask.send_from_directory(survivor_dir, filename, as_attachment=True)

@app.server.route('/download/evalsvcallers/<path:filename>')
def download_eval_file(filename):
    return flask.send_from_directory('./uploaded_files/evalsvcallers_output', filename, as_attachment=True)

@app.server.route('/download/visualization/<path:filename>')
def download_visualization_file(filename):
    return serve_file_for_download(filename)

################### RUN APP ########################

if __name__ == "__main__":
    app.run(debug=True, port=8040)
    #app.run_server(debug=True, port=8040)
