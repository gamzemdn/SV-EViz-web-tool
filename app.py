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
import io
import os
import re
import base64
import json
import pandas as pd
from layouts.survivor_layout import get_survivor_layout
from layouts.evalsvcallers_layout import get_evalsvcallers_layout
from layouts.visualize_layout import get_visualize_layout
from layouts.metrics_layout import get_metrics_layout
from layouts.survivor_functions import parse_uploaded_files as parse_survivor_uploaded_files, prepare_vcf_files_for_merge, run_survivor_merge, get_merge_preview
from layouts.evalsvcallers_functions import (
    save_file, save_custom_reference_file, parse_uploaded_files as parse_eval_uploaded_files,
    run_conversion, run_evaluation
)
from layouts.visualize_functions import (plot_circos, plot_manhattan, plot_manhattan_svlen, plot_clustergram, plot_sankey,load_vcf_dataframe, update_tracks, vcf_to_circos_json, detect_genome_version, save_file, get_variant_extraction_status, parse_uploaded_vcf, extract_variant_types, plot_vcf_data, serve_file_for_download)
from layouts.metrics_functions import (
    parse_contents, parse_reference_normal, parse_reference_na12878,
    parse_caller_vcf, parse_intersection_vcf, calculate_metrics_explicit,get_survivor_upload_section,get_evalsvcallers_upload_section,parse_evalsvcallers_file, generate_evalsvcallers_visuals,process_survivor_metrics_from_content, generate_survivor_visuals)
import plotly.express as px
from dash import dash_table

UPLOAD_DIRECTORY = "./uploaded_files/"
SURVIVOR_OUTPUT_DIR= os.path.join(UPLOAD_DIRECTORY, "survivor_output")
EVAL_OUTPUT_DIR = os.path.join(UPLOAD_DIRECTORY, "evalsvcallers_output")
VISUALIZATION_OUTPUT_DIR = os.path.join(UPLOAD_DIRECTORY, "visualization_output")
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
        ], width=3, style={
            'backgroundColor': '#f8f9fa',
            'height': 'auto',
            'minHeight': '100vh',
            'padding': '1rem',
            'overflowY': 'auto'
        }),
        dbc.Col(
            html.Div(id='tabs-content', style={'padding': '1rem'}),
            width=6
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
                html.Li("SURVIVOR: Upload SV caller vcf file along with a reference file or choose Univar reference file to perform merge operation.", style={'padding': '0.25rem'}),
                html.Li("EvalSVcallers: Upload SV caller vcf file for convertion before evaluation and upload/choose a reference file to perform evaluation.", style={'padding': '0.25rem'}),
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
                html.Li("SURVIVOR: Upload the caller VCF, SURVIVOR result VCF, and reference VCF to compute metrics.", style={'padding': '0.25rem'}),
                html.Li("EvalSVcallers: Upload the `.eval.txt` metrics file to instantly view tabular and graphical results.", style={'padding': '0.25rem'}),
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
                    {'label': 'Caller', 'value': 'caller'},
                    {'label': 'SURVIVOR', 'value': 'survivor'},
                    {'label': 'EvalSVcallers', 'value': 'evalsvcallers'}
                ],
                value='caller',
                labelStyle={'display': 'block', 'fontFamily': '\"Times New Roman\", Times, serif'}
            ),
            dcc.Upload(
                id='upload-data',
                children=html.Div(['📎 Drag and Drop or ', html.A('Select *vcf File Here')]),
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
                    {'label': 'SURVIVOR', 'value': 'survivor'},
                    {'label': 'EvalSVcallers', 'value': 'evalsvcallers'}
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
                        {'label': 'EvalSVcallers', 'value': 'evalsvcallers'}
                    ],
                    value=None,
                    inline=True,
                    labelStyle={'marginRight': '20px', 'marginLeft': '20px','fontFamily': '"Times New Roman", Times, serif'}
                )
            ], style={'whiteSpace': 'nowrap'}),
    
            html.Div(id="survivor-params-container", children=get_survivor_layout(), style={'display': 'none'}),
            html.Div(id="evalsvcallers-params-container", children=get_evalsvcallers_layout(), style={'display': 'none'})
        ], style={'width': '50%', 'display': 'inline-block', 'verticalAlign': 'top'})

    
        right_column = html.Div([
            dcc.Store(id="merged-path-store"), #survivor
            dcc.Store(id='converted-file-store'), #evalsvcallers
            dcc.Store(id="tp-fp-file-store"),
            dcc.Store(id="metrics-file-store"),
            dcc.Store(id='survivor-uploaded-paths', data=[]), #uploaded-files-list #converted-files-list
            dcc.Store(id='survivor-uploaded-paths-ref', data=[]),#reference_upload_section
            html.Div(id={'type': 'eval-output-ref', 'index': 0}),
            # Output area
            dcc.Loading(
                id="loading-survivor-output",
                type="default",  # or "circle", "dot", "graph"
                fullscreen=False,
                children=html.Div([
                    html.Div(id='survivor-output', style={'marginTop': '20px'}), #convert-status
                    html.Div(id='eval-output-list', style={'marginTop': '10px', 'fontFamily': '"Times New Roman", Times, serif'})
                ])
                
            ),
            # Separate visualization containers
            html.Div(id='survivor-visuals-container', children=[
                html.Div(id='comparison-extra-output')
            ], style={'display': 'block'}),


            html.Div(id='evalsvcallers-visuals-container', children=[
                
            ], style={'display': 'none'})

        ])
        return left_column, right_column

    else:
        return html.Div(""), html.Div("Unknown tab.")
# --- Place the toggle callback after layout ---
@app.callback(
    Output('survivor-visuals-container', 'style'),
    Output('evalsvcallers-visuals-container', 'style'),
    Input('comparison-type', 'value')
)
def toggle_visuals(selected_source):
    if selected_source == 'survivor':
        return {'display': 'block'}, {'display': 'none'}
    elif selected_source == 'evalsvcallers':
        return {'display': 'none'}, {'display': 'block'}
    return {'display': 'none'}, {'display': 'none'}
@app.callback(
    Output('survivor-params-container', 'style'),
    Output('evalsvcallers-params-container', 'style'),
    Input('comparison-type', 'value'),
)
def toggle_param_layout(selected_tool):
    if selected_tool == 'survivor':
        return {'display': 'block'}, {'display': 'none'}
    elif selected_tool == 'evalsvcallers':
        return {'display': 'none'}, {'display': 'block'}
    return {'display': 'none'}, {'display': 'none'}
@app.callback(
    Output('converted-file-store', 'data', allow_duplicate=True),
    Output('tp-fp-file-store', 'data', allow_duplicate=True),
    Output('metrics-file-store', 'data', allow_duplicate=True),
    Output('survivor-uploaded-paths', 'data', allow_duplicate=True),
    Output('merged-path-store', 'data', allow_duplicate=True),
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
        [], None
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

    # Manhattan plot
    manhattan_html = plot_manhattan(df, manhattan_svtypes, manhattan_threshold)

    return table_card, visuals_html, clustergram_html, circos_html, manhattan_html
    
@app.callback(
    Output('survivor-output', 'children', allow_duplicate=True),
    Output('survivor-uploaded-paths', 'data'),
    Output('survivor-uploaded-paths-ref', 'data'),
    Output('comparison-extra-output', 'children', allow_duplicate=True),
    Output('merged-path-store', 'data'),
    [Input('survivor-upload', 'contents'),
     Input('survivor-upload-ref', 'contents'),
     Input('merge-button', 'n_clicks')],
    [State('survivor-upload', 'filename'),
     State('survivor-upload', 'last_modified'),
     State('survivor-upload-ref', 'filename'),
     State('survivor-upload-ref', 'last_modified'),
     State('ref_source_selector', 'value'),
     State('param_dist', 'value'),
     State('param_callers', 'value'),
     State('param_type', 'value'),
     State('param_strand', 'value'),
     State('param_dup', 'value'),
     State('param_sv_size', 'value'),
     State('survivor-uploaded-paths', 'data'),
     State('survivor-uploaded-paths-ref', 'data'),
     State('main-tabs', 'active_tab')],
    prevent_initial_call=True
)
def handle_survivor_actions(contents, ref_contents, n_clicks,
                            names, dates, ref_filename, ref_dates,
                            ref_source, dist, callers, sv_type, strand, dup, sv_size,
                            uploaded_paths, uploaded_paths_ref, active_tab):
    comparison_output = dash.no_update
    merged_path= dash.no_update
    triggered_id = ctx.triggered_id

    # Handle uploads
    if triggered_id == 'survivor-upload' and contents:
        if isinstance(contents, str):
            contents, names, dates = [contents], [names], [dates]
        div, new_paths = parse_survivor_uploaded_files(contents, names, dates)
        return div, new_paths[-1:], dash.no_update,dash.no_update, dash.no_update

    if triggered_id == 'survivor-upload-ref' and ref_contents:
        if isinstance(ref_contents, str):
            ref_contents, ref_filename, ref_dates = [ref_contents], [ref_filename], [ref_dates]
        div, new_paths = parse_survivor_uploaded_files(ref_contents, ref_filename, ref_dates)
        return div, dash.no_update, new_paths[-1:],dash.no_update, dash.no_update

    # Handle merge
    if triggered_id == 'merge-button':
        
        params = {
            'use_univar': 1 if ref_source == 'univar' else 0,
            'max_distance': int(dist) if dist else 1000,
            'min_callers': callers,
            'type_match': int(sv_type) if sv_type else 1,
            'strand_match': int(strand) if strand else 1,
            'allow_duplicates': int(dup) if dup else 0,
            'min_sv_size': int(sv_size) if sv_size else 30,
        }

        vcf_files = []
        if uploaded_paths:
            vcf_files.append(os.path.abspath(uploaded_paths[-1]))
        if params['use_univar'] == 1:
            univar_path = os.path.abspath(os.path.join(UPLOAD_DIRECTORY, "univar.vcf"))
            if os.path.exists(univar_path):
                vcf_files.append(univar_path)
        elif uploaded_paths_ref:
            vcf_files.append(os.path.abspath(uploaded_paths_ref[-1]))

        if len(set(vcf_files)) != 2:
            return html.Div("⚠️ Please upload at least two different VCF files."), uploaded_paths, uploaded_paths_ref, dash.no_update, dash.no_update

        # Run SURVIVOR merge
        status, merged_path = run_survivor_merge(params, vcf_files)

        if not merged_path:
            return html.Div([status]), uploaded_paths, uploaded_paths_ref

        # ✅ Now load to get SVTYPEs for dropdowns
        df = load_vcf_dataframe(merged_path, "survivor")
        if df is None or df.empty:
            return html.Div([status, html.Br(), html.Div("⚠️ Merged file loaded but no variants found. Make sure reference version is hg38")]), uploaded_paths, uploaded_paths_ref, dash.no_update, dash.no_update

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
                if ref_source == 'univar':
                    univar_path = os.path.abspath(os.path.join(UPLOAD_DIRECTORY, "univar.vcf"))
                    if not os.path.exists(univar_path):
                        raise FileNotFoundError(f"Univar VCF not found at {univar_path}")
                    ref_content_b64 = read_file_as_base64(univar_path)
                else:
                    ref_content_b64 = read_file_as_base64(uploaded_paths_ref[-1]) if uploaded_paths_ref else None
                
                caller_content_b64 = read_file_as_base64(uploaded_paths[-1]) if uploaded_paths else None
                
                if not (ref_content_b64 and caller_content_b64):
                    raise ValueError("Uploaded caller or reference file content missing.")

                # Metrics
                metrics_df, error_msg = process_survivor_metrics_from_content(
                    ref_content_b64, caller_content_b64, inter_content_b64, 'normal'
                )
                metrics_html = (
                    generate_survivor_visuals(metrics_df)
                    if not error_msg else html.Div(error_msg, style={'color': 'red'})
                )

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
                    html.H3("📊 Automatic Metrics & Visuals", style={'marginTop': '20px'}),
                    metrics_html,
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
                comparison_output = html.Div(f"❌ Error generating metrics/visuals: {str(e)}")

        
        # ✅ Return merge status + filtering controls (visuals filled by next callback)
        return html.Div([status, html.Br()]), uploaded_paths, uploaded_paths_ref, comparison_output, merged_path

    return dash.no_update, uploaded_paths, uploaded_paths_ref,comparison_output, merged_path
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
    Output({'type': 'eval-output-ref', 'index': MATCH}, 'children'),
    Input({'type': 'dynamic-upload', 'index': MATCH}, 'contents'),
    State({'type': 'dynamic-upload', 'index': MATCH}, 'filename'),
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
    Input('tp-fp-file-store', 'data')
)
def update_manhattan_svtypes(tp_fp_path):
    if not tp_fp_path or not isinstance(tp_fp_path, str):
        return [], {i: str(i) for i in range(11)}
    df = load_vcf_dataframe(tp_fp_path, "evalsvcallers")
    if "SVTYPE" not in df.columns:
        return []
    svtypes = df["SVTYPE"].dropna().unique().tolist()
    return [{"label": sv, "value": sv} for sv in svtypes]
@app.callback(
    Output('eval-manhattan-slider', 'max'),
    Output('eval-manhattan-slider', 'marks'),
    Input('tp-fp-file-store', 'data')
)
def update_eval_threshold_slider_range(tp_fp_path):
    if not tp_fp_path or not isinstance(tp_fp_path, str):
        return [], {i: str(i) for i in range(11)}
    # Always keep slider clean and small
    marks = {i: str(i) for i in range(0, 11)}
    return 10, marks
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
    Output("reference_upload_section", "children"),
    Input("reference_choice", "value")
)
def toggle_reference_upload(reference_choice):
    if reference_choice == "custom":
        return html.Div([
            dcc.Upload(
                id={'type': 'dynamic-upload', 'index': 0},
                children=html.Div(['📎 Drag and Drop or ', html.A('Select VCF File')]),
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
            )
        #    html.Div(id={'type': 'dynamic-upload-status', 'index': 0})
        ])
    return ""


@app.callback(
    Output('converted-upload-section', 'style'),
    Input('use_conversion_toggle', 'value')
)
def toggle_converted_upload_section(selected):
    return {'display': 'none'}
################### VISUALIZATION FUNCTIONS ########################
@app.callback(
    Output('manhattan-threshold-slider', 'max'),
    Output('manhattan-threshold-slider', 'marks'),
    Input('filtered-data', 'data')
)
def update_threshold_slider_range(json_data):
    if not json_data:
        return 10, {i: str(i) for i in range(0, 11)}
    df = pd.read_json(io.StringIO(json_data))
    print(df["QUAL"])
    max_threshold = min(30, df["QUAL"].max() / 200 if "QUAL" in df.columns else 10)
    return max_threshold, {i: str(i) for i in range(0, int(max_threshold) + 1)}

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
    Input('filtered-data', 'data')
)
def update_manhattan_svtypes(json_data):
    if not json_data:
        return []
    df = pd.read_json(io.StringIO(json_data))
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
    if comp_type == 'survivor':
        return dash.no_update, ""  # ✅ Keep SURVIVOR, clear Eval
    elif comp_type == 'evalsvcallers':
        return "", dash.no_update  # ✅ Clear SURVIVOR, keep Eval
    return "", ""


@app.callback(
    Output('metrics-upload-section', 'children'),
    Input('metrics-input-source', 'value')
)
def render_upload_section(source):
    if source == 'survivor':
        return get_survivor_upload_section()
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
@app.callback(
    Output('output-metrics', 'children'),
    Input('process-button', 'n_clicks'),
    State('upload-reference', 'contents'),
    State('upload-caller', 'contents'),
    State('upload-intersection', 'contents'),
    State('reference-type', 'value'),
    prevent_initial_call=True
)
def process_survivor_metrics(n_clicks, ref_content, caller_content, inter_content, reference_type):
    if not (ref_content and caller_content and inter_content):
        return html.Div("⚠️ All three files must be uploaded.", style={'color': 'red','fontFamily': '\"Times New Roman\", Times, serif'})

    metrics_df, error_msg = process_survivor_metrics_from_content(ref_content, caller_content, inter_content, reference_type)

    if error_msg:
        return html.Div(error_msg, style={'color': 'red'})

    return generate_survivor_visuals(metrics_df)
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
