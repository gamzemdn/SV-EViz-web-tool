from pathlib import Path
import os
import shutil
import base64
import datetime
import json
import pandas as pd

import dash
from dash import html, dcc, dash_table, Input, Output, State
import dash_bootstrap_components as dbc


def _example_root(deps):
    return Path(deps["APP_ROOT"]).resolve() / "example_data"

def file_to_dash_contents(file_path, mime="text/plain"):
    with open(file_path, "rb") as f:
        encoded = base64.b64encode(f.read()).decode("utf-8")
    return f"data:{mime};base64,{encoded}"

def get_comparison_example_panel():
    return html.Div([
        html.Hr(),
        dbc.Alert(
            "Run the selected Comparison module using the built-in example dataset.",
            color="secondary",
            style={
                "fontSize": "13px",
                "fontFamily": '"Times New Roman", Times, serif',
                "padding": "8px"
            }
        ),
        dbc.Button(
            "Run selected example workflow",
            id="run-comparison-example",
            n_clicks=0,
            color="secondary",
            className="mt-1",
            style={'fontFamily': '"Times New Roman", Times, serif'}
        ),
        html.Div(
            id="comparison-example-status",
            style={
                "marginTop": "8px",
                "fontSize": "13px",
                "fontFamily": '"Times New Roman", Times, serif'
            }
        ),
        html.Hr()
    ])


def get_visualization_example_panel():
    return html.Div([
        html.Hr(),
        dbc.Alert(
            "Load a built-in example file for the selected Visualization source.",
            color="secondary",
            style={
                "fontSize": "13px",
                "fontFamily": '"Times New Roman", Times, serif',
                "padding": "8px"
            }
        ),
        dbc.Button(
            "Load visualization example",
            id="load-visualization-example",
            n_clicks=0,
            color="secondary",
            className="mt-1",
            style={'fontFamily': '"Times New Roman", Times, serif'}
        ),
        html.Div(
            id="visualization-example-status",
            style={
                "marginTop": "8px",
                "fontSize": "13px",
                "fontFamily": '"Times New Roman", Times, serif'
            }
        ),
        html.Hr()
    ])


def get_metrics_example_panel():
    return html.Div([
        html.Hr(),
        dbc.Alert(
            "Load a built-in example file for the selected Metrics module.",
            color="secondary",
            style={
                "fontSize": "13px",
                "fontFamily": '"Times New Roman", Times, serif',
                "padding": "8px"
            }
        ),
        dbc.Button(
            "Load metrics example",
            id="load-metrics-example",
            n_clicks=0,
            color="secondary",
            className="mt-1",
            style={'fontFamily': '"Times New Roman", Times, serif'}
        ),
        html.Div(
            id="metrics-example-status",
            style={
                "marginTop": "8px",
                "fontSize": "13px",
                "fontFamily": '"Times New Roman", Times, serif'
            }
        ),
        html.Hr()
    ])


def _make_run_dir(base_dir, prefix):
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = Path(base_dir).resolve() / f"{prefix}_{timestamp}"
    (run_dir / "inputs").mkdir(parents=True, exist_ok=True)
    (run_dir / "tables").mkdir(parents=True, exist_ok=True)
    (run_dir / "figures").mkdir(parents=True, exist_ok=True)
    (run_dir / "downloads").mkdir(parents=True, exist_ok=True)
    return run_dir


def _copy_to_inputs(src_path, run_dir):
    src_path = Path(src_path).resolve()
    dst_path = Path(run_dir) / "inputs" / src_path.name
    shutil.copy2(src_path, dst_path)
    return str(dst_path)


def _file_to_upload_contents(file_path):
    file_path = Path(file_path).resolve()
    encoded = base64.b64encode(file_path.read_bytes()).decode("utf-8")
    return "data:text/plain;base64," + encoded


def _missing_file_box(path):
    return html.Div(
        f"❌ Example file not found: {path}",
        style={"color": "red", "fontFamily": '"Times New Roman", Times, serif'}
    )


def register_example_workflows(app, deps):

    # ============================================================
    # 1) COMPARISON TAB — SURVIVOR EXAMPLE WORKFLOW
    # ============================================================
    @app.callback(
        Output("comparison-example-status", "children"),
        Output("survivor-uploaded-paths", "data", allow_duplicate=True),
        Output("merged-path-store", "data", allow_duplicate=True),
        Output("comparison-extra-output", "children", allow_duplicate=True),
        Input("run-comparison-example", "n_clicks"),
        State("comparison-type", "value"),
        prevent_initial_call=True
    )
    def run_comparison_example(n_clicks, selected_tool):
        if not n_clicks:
            raise dash.exceptions.PreventUpdate

        if selected_tool != "survivor":
            return (
                html.Div(
                    "⚠️ This example auto-run is currently connected to the SURVIVOR example. "
                    "Please select SURVIVOR to run the example workflow.",
                    style={"color": "orange"}
                ),
                dash.no_update,
                dash.no_update,
                dash.no_update
            )

        ex_root = _example_root(deps)
        survivor_dir = ex_root / "comparison" / "survivor"

        example_files = [
            survivor_dir / "Manta.vcf",
            survivor_dir / "Delly.vcf",
            survivor_dir / "Lumpy.vcf"
        ]

        for p in example_files:
            if not p.exists():
                return _missing_file_box(p), dash.no_update, dash.no_update, dash.no_update

        params = {
            "max_distance": 1000,
            "min_callers": 1,
            "type_match": 1,
            "strand_match": 1,
            "allow_duplicates": 0,
            "min_sv_size": 30,
        }

        status, merged_path = deps["run_survivor_merge"](
            params,
            [str(p.resolve()) for p in example_files]
        )

        if not merged_path or not os.path.exists(merged_path):
            return (
                html.Div([status, html.Br(), html.Div("❌ SURVIVOR example did not generate merged output.")]),
                [str(p.resolve()) for p in example_files],
                None,
                html.Div()
            )

        df = deps["load_vcf_dataframe"](merged_path, "survivor")
        if df is None or df.empty:
            return (
                html.Div([status, html.Br(), html.Div("⚠️ Merged example file was created, but no variants were loaded.")]),
                [str(p.resolve()) for p in example_files],
                merged_path,
                html.Div()
            )

        df["SVTYPE"] = df["INFO"].astype(str).str.extract(r"SVTYPE=([^;]+)")
        unique_svtypes = sorted(df["SVTYPE"].dropna().unique().tolist())

        # Basic visuals
        visuals_html = html.Div(deps["plot_vcf_data"](df))

        # Sankey
        try:
            sankey_html = deps["plot_sankey"](df.copy())
        except Exception as e:
            sankey_html = html.Div(f"❌ Sankey error: {e}", style={"color": "red"})

        # Clustergram
        try:
            valid_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
            df_cluster = df.copy()
            df_cluster["CHROM"] = df_cluster["CHROM"].astype(str).apply(
                lambda x: f"chr{x}" if not x.startswith("chr") else x
            )
            available_chroms = [c for c in valid_chroms if c in df_cluster["CHROM"].unique()]
            clustergram_html = (
                deps["plot_clustergram"](df_cluster, available_chroms)
                if available_chroms else html.Div("⚠️ No chromosomes for Clustergram.")
            )
        except Exception as e:
            clustergram_html = html.Div(f"❌ Clustergram error: {e}", style={"color": "red"})
            available_chroms = []

        # Circos
        try:
            json_path = os.path.join(
                os.path.dirname(merged_path),
                f"{os.path.basename(merged_path).replace('.vcf', '')}_circos.json"
            )
            deps["vcf_to_circos_json"](merged_path, json_path, "survivor")

            with open(json_path, "r") as f:
                circos_data = json.load(f)

            available_svtypes_circos = [track["name"] for track in circos_data.get("tracks", [])]
            svtype_colors_circos = {
                track["name"]: track.get("color", "#95a5a6")
                for track in circos_data.get("tracks", [])
            }

            circos_html = deps["plot_circos"](
                "histogram",
                available_svtypes_circos,
                circos_data,
                svtype_colors_circos
            )
        except Exception as e:
            circos_html = html.Div(f"❌ Circos error: {e}", style={"color": "red"})

        # Manhattan
        try:
            manhattan_html = deps["plot_manhattan"](df.copy(), unique_svtypes, 6)
        except Exception as e:
            manhattan_html = html.Div(f"❌ Manhattan error: {e}", style={"color": "red"})

        comparison_output = html.Div([
            html.H3("📊 Example SURVIVOR Workflow Output", style={"marginTop": "20px"}),

            dbc.Alert(
                "This example was generated automatically from the example_data/comparison/survivor files.",
                color="info"
            ),

            html.Div([
                html.Label("Filter by Variant Type:"),
                dcc.Dropdown(
                    id="comparison-svtype-filter",
                    options=[{"label": "ALL", "value": "ALL"}] +
                            [{"label": sv, "value": sv} for sv in unique_svtypes],
                    value="ALL",
                    clearable=False,
                    style={"width": "300px"}
                )
            ], style={"marginBottom": "20px"}),

            html.Div([
                dbc.Card([
                    dbc.CardHeader(html.H4("VCF File Preview (ALL)")),
                    dbc.CardBody([
                        dash_table.DataTable(
                            data=df.head(5).to_dict("records"),
                            columns=[{"name": col, "id": col} for col in df.columns],
                            style_table={"overflowX": "auto"},
                            style_cell={"textAlign": "left", "whiteSpace": "normal", "height": "auto"},
                            style_header={
                                "fontWeight": "bold",
                                "backgroundColor": "lightgrey",
                                "border": "1px solid black"
                            },
                            style_data={"border": "1px solid black"}
                        )
                    ])
                ], className="mt-4")
            ], id="preview-card"),

            html.Br(),
            html.H4("Basic Visualizations"),
            html.Div(visuals_html, id="visuals-output"),

            html.Br(),
            html.H4("Sankey Visualization"),
            sankey_html,

            html.Br(),
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
            html.H4("Circos Visualization"),
            circos_html,

            html.Br(),
            html.H4("Manhattan Visualization"),
            manhattan_html
        ])

        try:
            survivor_run_dir = Path(merged_path).resolve().parents[1]

            zip_path, saved_files = deps["create_run_export_package"](
                run_dir=survivor_run_dir,
                components=comparison_output,
                zip_prefix="survivor_example_exports"
            )

            download_card = deps["make_zip_download_card"](
                tool_name="survivor",
                run_dir=survivor_run_dir,
                zip_path=zip_path,
                saved_files=saved_files
            )

            comparison_output = html.Div([
                comparison_output,
                download_card
            ])

        except Exception as e:
            comparison_output = html.Div([
                comparison_output,
                html.Div(
                    f"⚠️ Export package could not be created: {e}",
                    style={"color": "orange", "marginTop": "10px"}
                )
            ])

        return (
            html.Div("✅ SURVIVOR example workflow completed.", style={"color": "green"}),
            [str(p.resolve()) for p in example_files],
            merged_path,
            comparison_output
        )

    # ============================================================
    # 2) VISUALIZATION TAB — LOAD EXAMPLE FILE AND TRIGGER BASIC
    # ============================================================
    @app.callback(
        Output("visualization-example-status", "children"),
        Output("vcf-files-list", "children", allow_duplicate=True),
        Output("uploaded-file-path", "data", allow_duplicate=True),
        Output("visualization-run-dir", "data", allow_duplicate=True),
        Output("variant-type-filter", "options", allow_duplicate=True),
        Output("variant-type-filter", "value", allow_duplicate=True),
        Output("filtered-data", "data", allow_duplicate=True),
        Output("viz_selector", "value"),
        Input("load-visualization-example", "n_clicks"),
        State("visualize-input-source", "value"),
        prevent_initial_call=True
    )
    def load_visualization_example(n_clicks, selected_source):
        if not n_clicks:
            raise dash.exceptions.PreventUpdate

        ex_root = _example_root(deps)

        if selected_source == "survivor":
            src = ex_root / "visualization" / "survivor" / "survivor_manta_delly_lumpy_merged_output.vcf"
            source_type = "survivor"

        elif selected_source == "evalsvcallers":
            src = ex_root / "visualization" / "evalsvcallers" / "NA12878_svaba.vcf.gz"
            source_type = "evalsvcallers"

        else:
            src = ex_root / "visualization" / "caller-truvari" / "caller_lumpy.vcf"
            source_type = "caller"

        if not src.exists():
            return _missing_file_box(src), dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        run_dir = _make_run_dir(deps["VISUALIZATION_OUTPUT_DIR"], "example_visualization_run")
        copied_path = _copy_to_inputs(src, run_dir)

        options, json_data = deps["extract_variant_types"](copied_path, source_type)

        file_list = html.Ul([
            html.Li(
                f"{src.name} loaded from example_data.",
                style={"fontFamily": '"Times New Roman", Times, serif'}
            )
        ])

        return (
            html.Div("✅ Visualization example loaded. Basic visuals will be displayed.", style={"color": "green"}),
            file_list,
            copied_path,
            str(run_dir),
            options,
            "ALL",
            json_data,
            "basic"
        )

    # ============================================================
    # 3) METRICS TAB — LOAD EXAMPLE METRICS FILE
    # ============================================================
    @app.callback(
        Output("metrics-example-status", "children"),
        Output("output-metrics", "children", allow_duplicate=True),
        Output("summary-visualization-output", "children", allow_duplicate=True),
        Input("load-metrics-example", "n_clicks"),
        State("metrics-input-source", "value"),
        prevent_initial_call=True
    )
    def load_metrics_example(n_clicks, selected_source):
        if not n_clicks:
            raise dash.exceptions.PreventUpdate

        if not selected_source:
            return (
                html.Div("⚠️ Please select Truvari or EvalSVcallers first.", style={"color": "orange"}),
                dash.no_update,
                dash.no_update
            )

        ex_root = _example_root(deps)
        run_dir = _make_run_dir(deps["METRICS_OUTPUT_DIR"], f"example_metrics_{selected_source}")

        if selected_source == "truvari":
            src = ex_root / "metrics" / "truvari_summary.json"

            if not src.exists():
                return _missing_file_box(src), dash.no_update, dash.no_update

            contents = _file_to_upload_contents(src)
            _copy_to_inputs(src, run_dir)

            df = deps["parse_truvari_file"](contents)
            metrics_components = deps["generate_truvari_visuals"](df)

            zip_path, saved_files = deps["create_run_export_package"](
                run_dir=run_dir,
                components=metrics_components,
                zip_prefix="truvari_metrics_example_exports"
            )

            download_card = deps["make_zip_download_card"](
                tool_name="metrics",
                run_dir=run_dir,
                zip_path=zip_path,
                saved_files=saved_files
            )

            return (
                html.Div("✅ Truvari metrics example loaded.", style={"color": "green"}),
                html.Div([metrics_components, download_card]),
                html.Div()
            )

        if selected_source == "evalsvcallers":
            src = ex_root / "metrics" / "converted_Manta.eval.txt"

            if not src.exists():
                return _missing_file_box(src), dash.no_update, dash.no_update

            contents = _file_to_upload_contents(src)
            _copy_to_inputs(src, run_dir)

            pivot_ref, df_long, df_block = deps["parse_evalsvcallers_file"](contents)
            metrics_components = deps["generate_evalsvcallers_visuals"](pivot_ref, df_long, df_block)

            zip_path, saved_files = deps["create_run_export_package"](
                run_dir=run_dir,
                components=metrics_components,
                zip_prefix="evalsvcallers_metrics_example_exports"
            )

            download_card = deps["make_zip_download_card"](
                tool_name="metrics",
                run_dir=run_dir,
                zip_path=zip_path,
                saved_files=saved_files
            )

            return (
                html.Div("✅ EvalSVcallers metrics example loaded.", style={"color": "green"}),
                html.Div(),
                html.Div([metrics_components, download_card])
            )

        return html.Div(), html.Div(), html.Div()