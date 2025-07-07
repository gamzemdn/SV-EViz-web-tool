from dash import Dash, html, dcc, Input, Output, State, ctx, get_asset_url, callback_context
import dash_bio as dashbio
import datetime
import dash_bootstrap_components as dbc
from dash_bootstrap_templates import load_figure_template
from dash.dependencies import MATCH
from dash import no_update
import dash
import flask
import io
import os
import re
import json
import pandas as pd
from survivor_layout import get_survivor_layout
from evalsvcallers_layout import get_evalsvcallers_layout
from visualize_layout import get_visualize_layout
from metrics_layout import get_metrics_layout
from survivor_functions import parse_uploaded_files, run_survivor_merge, get_merge_preview
from evalsvcallers_functions import ( save_file, save_custom_reference_file, parse_uploaded_files, run_conversion, run_evaluation )
from visualize_functions import (plot_circos, plot_manhattan, plot_clustergram, plot_sankey,load_vcf_dataframe, update_tracks, vcf_to_circos_json, detect_genome_version, save_file, get_variant_extraction_status, parse_uploaded_vcf, extract_variant_types, plot_vcf_data, serve_file_for_download)
from metrics_functions import (
    parse_contents, parse_reference_normal, parse_reference_na12878,
    parse_caller_vcf, parse_intersection_vcf, calculate_metrics_explicit,get_survivor_upload_section,get_evalsvcallers_upload_section,parse_evalsvcallers_file, generate_evalsvcallers_visuals,process_survivor_metrics_from_content, generate_survivor_visuals)
import plotly.express as px
from dash import dash_table

UPLOAD_DIRECTORY = "./uploaded_files/"
VISUALIZATION_OUTPUT_DIR = os.path.join(UPLOAD_DIRECTORY, "visualization_output")
load_figure_template('LITERA')
app = Dash(__name__, suppress_callback_exceptions=True, assets_folder='assets', external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            dbc.Row([
                dbc.Col(html.Img(
                    src=dash.get_asset_url('sv_icon.svg'),
                    height='40px',
                    style={
                        'marginRight': '8px',
                        'marginLeft': '4px',
                        'borderRadius': '4px',
                        'border': '1px solid #ccc'
                    }
                ), width='auto'),
                dbc.Col(html.H2(
                    "SV Comparison and Visualization",
                    style={'padding': '1rem', 'fontSize': '24px', 'color': '#2c3e50', 'fontFamily': '\"Times New Roman\", Times, serif', 'fontWeight': 'bold'}
                ))
            ], align='center')
        ], width=12)
    ]),
    dbc.Row([
        dbc.Col([
            dbc.Tabs(
                [
                    dbc.Tab(label="Welcome", tab_id="tab-welcome"),
                    dbc.Tab(label="Comparison", tab_id="tab-comparison"),
                    dbc.Tab(label="Visualization", tab_id="tab-visualization"),
                    dbc.Tab(label="Metrics", tab_id="tab-metrics"),
                ],
                id="main-tabs",
                active_tab="tab-welcome",
                style={"width": "100%", "minWidth": "0", "flexWrap": "wrap", 'fontFamily': '\"Times New Roman\", Times, serif'}
            ),
            html.Div(id='left-extra', style={'marginTop': '1rem'})
        ], width=3, style={'backgroundColor': '#f8f9fa', 'height': '100vh', 'padding': '1rem'}),
        dbc.Col(
            html.Div(id='tabs-content', style={'padding': '1rem'}),
            width=9
        )
    ])
], fluid=True)

################### SURVIVOR FUNCTIONS ########################

@app.callback(
    Output('survivor-output', 'children'),
    [Input('survivor-upload', 'contents'),
     Input('merge-button', 'n_clicks')],
    [State('survivor-upload', 'filename'),
     State('survivor-upload', 'last_modified'),
     State('param_univar', 'value'),
     State('param_dist', 'value'),
     State('param_callers', 'value'),
     State('param_type', 'value'),
     State('param_strand', 'value'),
     State('param_dup', 'value'),
     State('param_sv_size', 'value')],
    prevent_initial_call=True
)
def handle_survivor_actions(contents, n_clicks, names, dates,
                            univar, dist, callers, sv_type, strand, dup, sv_size):
    triggered_id = ctx.triggered_id

    if triggered_id == 'survivor-upload' and contents:
        return parse_uploaded_files(contents, names, dates)

    if triggered_id == 'merge-button':
        params = {
            'param_univar': int(univar) if univar else 0,
            'max_distance': int(dist) if dist else 1000,
            'min_callers': int(callers) if callers else 2,
            'type_match': int(sv_type) if sv_type else 1,
            'strand_match': int(strand) if strand else 1,
            'allow_duplicates': int(dup) if dup else 0,
            'min_sv_size': int(sv_size) if sv_size else 30,
        }
        status = run_survivor_merge(params)
        preview = get_merge_preview()
        download_link = html.A("⬇️ Download merged file",
                               href='/download/survivor/merged_output.vcf',
                               target='_blank',
                               style={'display': 'block', 'marginTop': '10px'})
        return html.Div([html.Div(status), preview, download_link])

    return html.Div("No action detected.")

################### EVALSVCALLERS FUNCTIONS ########################

@app.callback(
    Output({'type': 'dynamic-upload-status', 'index': MATCH}, 'children'),
    Input({'type': 'dynamic-upload', 'index': MATCH}, 'contents'),
    State({'type': 'dynamic-upload', 'index': MATCH}, 'filename'),
    prevent_initial_call=True
)
def handle_reference_upload(contents, filename):
    if not contents:
        return "❌ No file uploaded."

    result_msg = save_custom_reference_file(contents, filename)
    return html.Div(result_msg, style={'color': 'green' if '✅' in result_msg else 'red'})


# ========== Upload + Convert ==========
@app.callback(
    Output('eval-output', 'children'),
    [Input('eval-upload', 'contents'),
     Input('convert-button', 'n_clicks')],
    [State('eval-upload', 'filename'),
     State('eval-upload', 'last_modified'),
     State('caller_tool', 'value')],
    prevent_initial_call=True
)
def handle_evalsvcallers_actions(contents, convert_clicks, names, dates, caller_tool):
    triggered_id = ctx.triggered_id

    if 'eval-upload' in str(triggered_id) and contents:
        return parse_uploaded_files([contents], [names], [dates])

    if triggered_id == 'convert-button':
        if names:
            filenames = [names] if isinstance(names, str) else names
            result, file_list = run_conversion(caller_tool, filenames)
            return html.Div([result, file_list])
        else:
            return html.Div("❌ Please upload a file before converting.")

    return html.Div("No action detected.")

@app.callback(
    Output('converted-files-list', 'children'),
    Input('upload-converted-data', 'contents'),
    State('upload-converted-data', 'filename'),
    State('upload-converted-data', 'last_modified'),
    prevent_initial_call=True
)
def update_converted_files(list_of_contents, list_of_names, list_of_dates):
    if not list_of_contents:
        return html.Div("No converted files uploaded yet.")

    uploaded_files = []
    for content, name, date in zip(list_of_contents, list_of_names, list_of_dates):
        formatted_date = datetime.datetime.fromtimestamp(date).strftime("%Y-%m-%d %H:%M:%S")
        uploaded_files.append(html.Li(f"{name} (Saved at: {formatted_date})"))

    return html.Div(f"Converted file '{list_of_names}' uploaded successfully!")

# ========== EvalSVcallers Run ==========
@app.callback(
    Output('eval-status', 'children'),
    Input('eval-button', 'n_clicks'),
    prevent_initial_call=True
)
def show_running_status(n):
    return "⏳ Running EvalSVcallers..."
    
@app.callback(
    [Output('basic_params', 'style'),
     Output('advanced_params', 'style')],
    Input('param_mode', 'value')
)
def toggle_param_mode(mode):
    if mode == 'advanced':
        return {'display': 'block', 'marginTop': '10px'}, {'display': 'block', 'marginTop': '10px'}
    else:
        return {'display': 'block', 'marginTop': '10px'}, {'display': 'none'}
        
@app.callback(
    Output('eval-output-list', 'children'),
    Input('eval-button', 'n_clicks'),
    State("upload-converted-data", "filename"),
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
    State('param_rl', 'value'),
    State('param_rxl', 'value'),
    State('param_mr', 'value'),
    State('param_rs', 'value'),
    State('param_mins', 'value'),
    State('param_eg', 'value'),
    State('param_eb', 'value'),
    State('param_i', 'value'),
    State('param_y', 'value'),
    State('param_of', 'value'),
    State('param_sm', 'value'),
    State('param_st_advanced', 'value'),
    State('param_l_advanced', 'value'),
    State('param_xl_advanced', 'value'),
    State('param_mo_advanced', 'value'),
    prevent_initial_call=True
)
def run_evalsvcallers(n_clicks, selected_converted_file, reference_choice, caller_tool,
                      st_basic, l_basic, xl_basic, rl_basic, rxl_basic, mo_basic, of_basic,
                      parent1_content, parent1, parent2_content, parent2, rb_content, rb,
                      c, rl, rxl, mr, rs, mins, eg, eb, i, y, of, sm, st_advanced, l_advanced, xl_advanced, mo_advanced):

    return run_evaluation(
        selected_converted_file, reference_choice, caller_tool,
        st_basic, l_basic, xl_basic, rl_basic, rxl_basic, mo_basic, of_basic,
        parent1_content, parent1, parent2_content, parent2, rb_content, rb,
        c, rl, rxl, mr, rs, mins, eg, eb, i, y, of, sm, st_advanced, l_advanced, xl_advanced, mo_advanced
    )
@app.callback(
    Output("reference_upload_section", "children"),
    Input("reference_choice", "value")
)
def toggle_reference_upload(reference_choice):
    if reference_choice == "custom":
        return html.Div([
            dcc.Upload(
                id={'type': 'dynamic-upload', 'index': 0},
                children=html.Div(['Drag and Drop or ', html.A('Select VCF File')]),
                style={
                    'width': '460px',
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
            html.Div(id={'type': 'dynamic-upload-status', 'index': 0})
        ])
    return ""

@app.callback(
    Output('conversion-section', 'style'),
    Input('use_conversion_toggle', 'value')
)
def toggle_conversion_section(selected):
    if 'use' in selected:
        return {'display': 'block'}
    return {'display': 'none'}

################### VISUALIZATION FUNCTIONS ########################
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
    if "manhattan" in selected:
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
"""
@app.callback(
    [Output('vcf-table-container', 'children'),
     Output('vcf-figures-container', 'children')],
    Input('variant-type-filter', 'value'),
    State('filtered-data', 'data')
)
def update_visualizations(selected_svtype, json_data):
    if not json_data:
        return html.Div(""), ""

    df = pd.read_json(io.StringIO(json_data))

    if selected_svtype and selected_svtype != "ALL":
        df = df[df["SVTYPE"] == selected_svtype]

    figures = plot_vcf_data(df)

    return (
        dbc.Card([
            dbc.CardHeader(html.H4(f"Filtered VCF File Preview ({selected_svtype})")),
            dbc.CardBody([
                dcc.Loading(
                    dash.dash_table.DataTable(
                        data=df.head(5).to_dict('records'),
                        columns=[{"name": col, "id": col} for col in df.columns],
                        style_table={'overflowX': 'auto'},
                        style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'},
                        style_header={'fontWeight': 'bold', 'backgroundColor': 'lightgrey', 'border': '1px solid black'},
                        style_data={'border': '1px solid black'}
                    )
                )
            ])
        ], className="mt-4"),
        dbc.Card([
            dbc.CardHeader(html.H4("Filtered Visualizations")),
            dbc.CardBody(figures)
        ], className="mt-4")
    )

"""
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
    if "basic" in selected:
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

        except Exception as e:
            visuals.append(html.Div(f"❌ Basic visuals error: {e}"))
        
    # ✅ Sankey
    if "sankey" in selected:
        try:
            visuals.append(plot_sankey(df))
        except Exception as e:
            visuals.append(html.Div(f"❌ Sankey error: {e}"))

    # ✅ Clustergram
    if "clustergram" in selected:
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
        except Exception as e:
            visuals.append(html.Div(f"❌ Clustergram error: {e}"))
    # ✅ Manhattan (you must get 'svtypes' and 'threshold' from UI inputs later)
    if "manhattan" in selected:
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
        except Exception as e:
            visuals.append(html.Div(f"❌ Manhattan error: {e}"))

    # ✅ Circos
    if "circos" in selected:
        try:
            if file_path.endswith(".vcf"):
                vcf_filename = os.path.basename(file_path)
                json_filename = vcf_filename.replace(".vcf", "_circos.json")
                json_path = os.path.join(VISUALIZATION_OUTPUT_DIR, json_filename)

                if not os.path.exists(json_path):
                    vcf_to_circos_json(file_path, json_path,source_type)
                    visuals.append(html.Div(f"✅ Circos JSON saved: {os.path.basename(json_path)} (Full path: {json_path})"))
                else:
                    visuals.append(html.Div(f"✅ Using existing Circos JSON file: {os.path.basename(json_path)} (Full path: {json_path})"))

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
    from visualize_functions import load_vcf_dataframe, plot_clustergram

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
################### MAIN LAYOUT ROUTING ########################

@app.callback(
    [Output('left-extra', 'children'),
     Output('tabs-content', 'children')],
    Input('main-tabs', 'active_tab'),
)
def update_layout(active_tab):
    if active_tab == 'tab-welcome':
        return html.Div(""), html.Div([
            html.H2("Welcome to Structural Variant Comparison and Visualization Tool!",
                    style={'padding': '1rem', 'fontSize': '20px', 'fontFamily': '\"Times New Roman\", Times, serif', 'fontWeight': 'bold'}),
            html.P("Select an option from the tabs on the left.",
                   style={'padding': '1rem', 'fontSize': '15px', 'fontFamily': '\"Times New Roman\", Times, serif'}),
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
            dcc.Checklist(
                id="viz_selector",
                options=[
                    {"label": " Data Table & Basic Visuals", "value": "basic"},
                    {"label": " Sankey", "value": "sankey"},
                    {"label": " Circos", "value": "circos"},
                    {"label": " Clustergram", "value": "clustergram"},
                    {"label": " Manhattan", "value": "manhattan"},
                ],
                value=[],
                labelStyle={"display": "block", 'fontFamily': '"Times New Roman", Times, serif'}
            ),
            dcc.Store(id='selected-input-source', data='caller'),
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
                   style={'padding': '0.5rem', 'fontSize': '15px', 'fontFamily': '\"Times New Roman\", Times, serif'}),
            dcc.RadioItems(
                id='comparison-type',
                options=[
                    {'label': 'SURVIVOR', 'value': 'survivor'},
                    {'label': 'EvalSVcallers', 'value': 'evalsvcallers'}
                ],
                value=None,
                labelStyle={'display': 'block', 'fontFamily': '\"Times New Roman\", Times, serif'}
            )
        ])
        right_column = html.Div([
            html.Div(id='comparison-content')
        ])
        return left_column, right_column

    else:
        return html.Div(""), html.Div("Unknown tab.")

@app.callback(
    Output('comparison-content', 'children'),
    Input('comparison-type', 'value'),
    prevent_initial_call=True
)
def update_comparison_section(selection):
    if selection == 'evalsvcallers':
        return get_evalsvcallers_layout()
    else:
        return get_survivor_layout()

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
    app.run_server(debug=True, port=8040)
