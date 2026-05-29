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
EXAMPLE_CREATE_EXPORTS = False

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
        Output("survivor-output", "children", allow_duplicate=True),
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
            # Other tools are handled by their own dedicated callbacks below.
            raise dash.exceptions.PreventUpdate

        ex_root = _example_root(deps)
        survivor_dir = ex_root / "comparison" / "survivor"

        example_files = [
            survivor_dir / "Manta.vcf",
            survivor_dir / "Delly.vcf",
            survivor_dir / "Lumpy.vcf"
        ]

        for p in example_files:
            if not p.exists():
                return _missing_file_box(p), dash.no_update, dash.no_update, dash.no_update, dash.no_update

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

        uploaded_example_paths = [str(p.resolve()) for p in example_files]

        if not merged_path or not os.path.exists(merged_path):
            msg = html.Div([status, html.Br(), html.Div("❌ SURVIVOR example did not generate merged output.")])
            return msg, msg, uploaded_example_paths, None, html.Div()

        df = deps["load_vcf_dataframe"](merged_path, "survivor")
        if df is None or df.empty:
            msg = html.Div([status, html.Br(), html.Div("⚠️ Merged example file was created, but no variants were loaded.")])
            return msg, msg, uploaded_example_paths, merged_path, html.Div()

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
            html.H3("📊 Example SURVIVOR Workflow Output", style={}),

#            dbc.Alert(
#                "This example was generated automatically from the example_data/comparison/survivor files.",
#                color="info"
#            ),

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

        if EXAMPLE_CREATE_EXPORTS:
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
            html.Div(),
            html.Div([status, html.Br()]),
            uploaded_example_paths,
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
            src = ex_root / "visualization" / "evalsvcallers" / "NA12878_manta.TF.vcf"
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
                html.Div([metrics_components]),
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
                html.Div(),
                html.Div([metrics_components])
            )

        return html.Div(), html.Div(), html.Div()

    # ============================================================
    # 4) COMPARISON TAB — TRUVARI EXAMPLE WORKFLOW
    # ============================================================
    # Reference FASTAs (~3 GB each) are not bundled with the repo or Docker
    # image, so we cannot run `truvari bench` live in the example. Instead,
    # the example uses a pre-computed Truvari bench output (tp-base/tp-comp/
    # fp/fn/summary.json) generated from a real run, and routes those files
    # through the SAME metrics + visualization logic as the live pipeline
    # (verbatim copy of app.py post-processing) so the reviewer sees an
    # output identical to a manual `Run Truvari bench`.
    @app.callback(
        Output("tru-status", "children", allow_duplicate=True),
        Output("truvari-output-list", "children", allow_duplicate=True),
        Output("truvari-metrics-preview", "children", allow_duplicate=True),
        Output("truvari-visual-output", "children", allow_duplicate=True),
        Input("run-comparison-example", "n_clicks"),
        State("comparison-type", "value"),
        prevent_initial_call=True
    )
    def run_truvari_example(n_clicks, selected_tool):
        if not n_clicks:
            raise dash.exceptions.PreventUpdate
        if selected_tool != "truvari":
            raise dash.exceptions.PreventUpdate

        # ------ Local imports required by the verbatim app.py block below ------
        import plotly.express as px
        import plotly.graph_objects as go
        import numpy as np
        import re
        import dash_bio as dashbio
        from dash_bio import Clustergram

        # Aliases pulled from deps so the verbatim block can use the same names
        load_vcf_dataframe        = deps["load_vcf_dataframe"]
        generate_truvari_visuals  = deps["generate_truvari_visuals"]
        plot_manhattan_truvari    = deps["plot_manhattan_truvari"]
        vcf_to_circos_json_truvari = deps["vcf_to_circos_json_truvari"]
        create_run_export_package = deps["create_run_export_package"]
        make_zip_download_card    = deps["make_zip_download_card"]
        get_manhattan_controls    = deps["get_manhattan_controls"]

        # ------ Locate and stage precomputed Truvari bench output ------
        ex_root = _example_root(deps)
        precomp_dir = ex_root / "comparison" / "truvari" / "precomputed_output_for_example_workflow"

        required = {
            "summary": precomp_dir / "summary.json",
            "tp_base": precomp_dir / "tp-base.vcf.gz",
            "tp_comp": precomp_dir / "tp-comp.vcf.gz",
            "fp":      precomp_dir / "fp.vcf.gz",
            "fn":      precomp_dir / "fn.vcf.gz",
        }
        missing = [str(p) for p in required.values() if not p.exists()]
        if missing:
            err = html.Div([
                html.Div("❌ Truvari example output files not found:"),
                html.Ul([html.Li(m) for m in missing]),
                html.Div(
                    "Please place a pre-computed Truvari bench output in "
                    f"{precomp_dir}/ (summary.json, tp-base.vcf.gz, "
                    "tp-comp.vcf.gz, fp.vcf.gz, fn.vcf.gz)."
                )
            ], style={"color": "red"})
            return err, html.Ul([]), html.Div(), html.Div()

        run_dir = _make_run_dir(
            Path(deps["APP_ROOT"]) / "uploaded_files" / "truvari_output",
            "example_truvari_run"
        )
        output_dir = run_dir / "truvari_output"
        output_dir.mkdir(parents=True, exist_ok=True)

        tru_paths = {}
        for key, src in required.items():
            dst = output_dir / src.name
            shutil.copy2(src, dst)
            tru_paths[key] = str(dst.resolve())

        # `pairs` is referenced inside the verbatim advanced_truvari_section
        # below; it is built from tru_paths exactly as the live callback does.
        pairs = [(tru_paths.get("tp_base"), "TP-BASE"),
                 (tru_paths.get("tp_comp"), "TP-COMP")]

        # Status banner for the top of the page
        status_div = html.Div(
            "✅ Truvari example workflow loaded (pre-computed bench output).",
            style={"color": "green"}
        )

        # ============================================================
        # VERBATIM COPY OF app.py TRUVARI POST-PROCESSING (lines 1976-2628)
        # Only the run_truvari_pipeline() call at the top is skipped — we
        # already have tru_paths from the bundled example output above.
        # ============================================================
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
                html.H4("Manhattan Visualization"),
                get_manhattan_controls("tru-manhattan-slider", "tru-manhattan-selector", prefix="tru-"),
                html.Div(manhattan_div, id="tru-manhattan-output")
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
        if EXAMPLE_CREATE_EXPORTS:
            try:
                summary_path = tru_paths.get("summary")
                truvari_run_dir = Path(summary_path).resolve().parents[1]
            
                export_components = html.Div([
                    metrics_preview,
                    visuals_block
                ])
            
                zip_path, saved_files = create_run_export_package(
                    run_dir=truvari_run_dir,
                    components=export_components,
                    zip_prefix="truvari_exports"
                )
            
                download_card = make_zip_download_card(
                    tool_name="truvari",
                    run_dir=truvari_run_dir,
                    zip_path=zip_path,
                    saved_files=saved_files
                )
            
                visuals_block = html.Div([
                    visuals_block,
                    download_card
                ])
            
            except Exception as e:
                visuals_block = html.Div([
                    visuals_block,
                    html.Div(
                        f"⚠️ Export package could not be created: {e}",
                        style={"color": "orange", "marginTop": "10px"}
                    )
                ])
        # Dosya listesi istenmiyor -> boş
        return status_div, html.Ul([]), metrics_preview, visuals_block

    # ============================================================
    # 5) COMPARISON TAB — EVALSVCALLERS EXAMPLE WORKFLOW
    # ============================================================
    # EvalSVcallers is light-weight (no reference FASTA needed), so the
    # example runs the actual convert + evaluate pipeline on the bundled
    # HG00514 Manta calls + GRCh38 truthset. The post-processing block
    # below is a verbatim copy of run_evalsvcallers_combined in app.py so
    # the produced visuals are IDENTICAL to a manual `Run EvalSVcallers`.
    @app.callback(
        Output("eval-output-list", "children", allow_duplicate=True),
        Output("tp-fp-file-store", "data", allow_duplicate=True),
        Output("metrics-file-store", "data", allow_duplicate=True),
        Output("evalsvcallers-visuals-container", "children", allow_duplicate=True),
        Input("run-comparison-example", "n_clicks"),
        State("comparison-type", "value"),
        prevent_initial_call=True
    )
    def run_evalsvcallers_example(n_clicks, selected_tool):
        if not n_clicks:
            raise dash.exceptions.PreventUpdate
        if selected_tool != "evalsvcallers":
            raise dash.exceptions.PreventUpdate

        # ------ Local imports + aliases for the verbatim app.py block ------
        import base64
        load_vcf_dataframe              = deps["load_vcf_dataframe"]
        plot_vcf_data                   = deps["plot_vcf_data"]
        plot_sankey                     = deps["plot_sankey"]
        plot_clustergram                = deps["plot_clustergram"]
        plot_circos                     = deps["plot_circos"]
        plot_manhattan                  = deps["plot_manhattan"]
        vcf_to_circos_json              = deps["vcf_to_circos_json"]
        parse_evalsvcallers_file        = deps["parse_evalsvcallers_file"]
        generate_evalsvcallers_visuals  = deps["generate_evalsvcallers_visuals"]
    #    create_run_export_package       = deps["create_run_export_package"]
    #    make_zip_download_card          = deps["make_zip_download_card"]
        get_manhattan_controls          = deps["get_manhattan_controls"]

        # ------ Stage bundled example files into the eval upload cache ------
        ex_root = _example_root(deps)
        eval_dir = ex_root / "comparison" / "evalsvcallers"
        input_vcf = eval_dir / "HG00514_Manta.vcf"
        truthset_vcf = eval_dir / "HG00514_GRCh38_truthset.vcf"

        missing = [p for p in [input_vcf, truthset_vcf] if not p.exists()]
        if missing:
            err = html.Div([
                html.Div("❌ EvalSVcallers example files not found:"),
                html.Ul([html.Li(str(p)) for p in missing])
            ], style={"color": "red"})
            return err, dash.no_update, dash.no_update, dash.no_update

        try:
            from layouts.evalsvcallers_functions import EVAL_UPLOAD_CACHE_DIR
        except Exception:
            EVAL_UPLOAD_CACHE_DIR = None

        if not EVAL_UPLOAD_CACHE_DIR:
            err = html.Div(
                "❌ Could not access EvalSVcallers upload cache directory.",
                style={"color": "red"}
            )
            return err, dash.no_update, dash.no_update, dash.no_update

        os.makedirs(EVAL_UPLOAD_CACHE_DIR, exist_ok=True)
        staged_input = Path(EVAL_UPLOAD_CACHE_DIR) / input_vcf.name
        shutil.copy2(input_vcf, staged_input)
        staged_ref = Path(EVAL_UPLOAD_CACHE_DIR) / "reference.vcf"
        shutil.copy2(truthset_vcf, staged_ref)
        uploaded_filename = input_vcf.name

        # ------ 1) CONVERSION ------
        try:
            conv_result, converted_file, timestamp = deps["run_eval_conversion"](
                "Manta", uploaded_filename
            )
        except Exception as e:
            err = html.Div(
                f"❌ EvalSVcallers conversion error: {e}",
                style={"color": "red"}
            )
            return err, dash.no_update, dash.no_update, dash.no_update

        if not converted_file or not os.path.exists(converted_file):
            err = html.Div([
                html.Div("❌ EvalSVcallers conversion did not produce an output file."),
                conv_result if conv_result else html.Div()
            ], style={"color": "red"})
            return err, dash.no_update, dash.no_update, dash.no_update

        # ------ 2) EVALUATION ------
        # Pass None for every UI parameter so evaluate_SV_callers.pl uses
        # its own built-in defaults (mirroring an empty manual form).
        try:
            status_html, tp_fp_path, metrics_path = deps["run_eval_evaluation"](
                converted_file,   # selected_converted_file
                "custom",         # reference_choice -> staged truthset VCF
                "Manta",          # caller_tool
                timestamp,
                None, None, None, None, None, None, None,
                None, None, None, None, None, None,
                None, None, None, None, None, None, None, None, None
            )
        except Exception as e:
            err = html.Div([
                html.Div(f"❌ EvalSVcallers evaluation error: {e}"),
                conv_result
            ], style={"color": "red"})
            return err, dash.no_update, dash.no_update, dash.no_update

        if not tp_fp_path or not os.path.exists(tp_fp_path) \
                or not metrics_path or not os.path.exists(metrics_path):
            err = html.Div([
                html.Div("❌ EvalSVcallers evaluation did not produce expected output files."),
                status_html if status_html else html.Div()
            ], style={"color": "red"})
            return err, dash.no_update, dash.no_update, dash.no_update

        # ============================================================
        # VERBATIM COPY OF app.py EVALSVCALLERS POST-PROCESSING
        # (lines 1440-1587 of run_evalsvcallers_combined) — wrapped in
        # try/except for safety. Variable names match the source 1:1 so
        # the block can be pasted without any modification.
        # ============================================================
        try:
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
            eval_run_output_dir = os.path.dirname(tp_fp_path)
            json_path = os.path.join(
                eval_run_output_dir,
                f"{os.path.basename(tp_fp_path)}_circos.json"
            )
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
                get_manhattan_controls("eval-manhattan-slider", "eval-manhattan-selector", prefix="eval-"),
                html.Div(manhattan_html, id="eval-manhattan-output")
    
            ])
            if EXAMPLE_CREATE_EXPORTS:
                try:
                    eval_run_dir = Path(metrics_path).resolve().parents[1]
            
                    zip_path, saved_files = create_run_export_package(
                        run_dir=eval_run_dir,
                        components=visuals_ui,
                        zip_prefix="evalsvcallers_example_exports"
                    )
            
                    download_card = make_zip_download_card(
                        tool_name="evalsvcallers",
                        run_dir=eval_run_dir,
                        zip_path=zip_path,
                        saved_files=saved_files
                    )
            
                    visuals_ui = html.Div([
                        visuals_ui,
                        download_card
                    ])
            
                except Exception as e:
                    visuals_ui = html.Div([
                        visuals_ui,
                        html.Div(
                            f"⚠️ Export package could not be created: {e}",
                            style={"color": "orange", "marginTop": "10px"}
                        )
                    ])
            
            return html.Div([status_html]), tp_fp_path, metrics_path, visuals_ui
        except Exception as e:
            # Post-processing failed; show the error but still expose the
            # tp-fp and metrics files so the user can pick up from them.
            err_block = html.Div([
                html.Div(f"❌ EvalSVcallers post-processing error: {e}"),
                status_html if status_html else html.Div()
            ], style={"color": "red"})
            return err_block, tp_fp_path, metrics_path, html.Div()
