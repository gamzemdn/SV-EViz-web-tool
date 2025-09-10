import os
import io
import base64
import pandas as pd
import plotly.express as px
from collections import defaultdict
from dash import html, dcc
import json
from dash import dash_table
import numpy as np
import dash_bootstrap_components as dbc
UPLOAD_DIRECTORY = "./uploaded_files/"
METRIC_OUTPUT_DIR = os.path.join(UPLOAD_DIRECTORY, "metrics_output")
os.makedirs(UPLOAD_DIRECTORY, exist_ok=True)
os.makedirs(METRIC_OUTPUT_DIR, exist_ok=True)
svtypes = ['DEL', 'DUP', 'INS', 'INV', 'BND', 'CNV']  # or whatever types you have
base_colors = px.colors.qualitative.Vivid.copy()  # or px.colors.qualitative.Bold
# Optionally remove yellowish tone (like '#FECB52')
exclude_colors = ['rgb(218, 165, 27)']  # Plotly yellow
color_palette = [c for c in base_colors if c not in exclude_colors]
svtype_color_map = {
    sv: color_palette[i % len(color_palette)] for i, sv in enumerate(svtypes)
}



#size_classes = ["A", "SS", "S", "M", "L"]
#size_color_discrete_map = {
#    size_class: color_palette[i % len(color_palette)]
#    for i, size_class in enumerate(size_classes)
#}
#tsml_order      = ["Total", "tiny (<= 100 bp)", "short (<= 1.0 kp)", "middle (<= 100 kb)", "large (> 100 kb)"]  # bar plot

#size_color_discrete_map = {size_classes[i]: color_palette[i] for i in range(len(size_classes))}
#tsml_color_discrete_map = {tsml_order[i]: color_palette[i] for i in range(len(tsml_order))}


metric_classes = ["Precision", "Recall", "F1-Score"]
metric_color_map = {
    metric_class: color_palette[i % len(color_palette)]
    for i, metric_class in enumerate(metric_classes)
}

# --- Decode uploaded content ---
def parse_contents(contents):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    return io.StringIO(decoded.decode('utf-8'))

# --- Reference parsers ---
def parse_reference_normal(file_like_obj):
    lines = [line for line in file_like_obj if not line.startswith('##')]
    vcf_rows = [line.strip().split('\t') for line in lines if not line.startswith('#')]

    if vcf_rows:
        vcf_data = pd.DataFrame(vcf_rows, columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'])
        vcf_data['SVTYPE'] = vcf_data['INFO'].str.extract(r'SVTYPE=([^;]+)')
        for svtype in vcf_data['SVTYPE'].dropna().unique():
            subset = vcf_data[vcf_data['SVTYPE'] == svtype]
            file_path = os.path.join(METRIC_OUTPUT_DIR, f'{svtype}_variants.csv')
            subset.to_csv(file_path, index=False)
        return vcf_data
    else:
        return pd.DataFrame()

def parse_reference_na12878(file_like_obj):
    lines = [line.strip() for line in file_like_obj if line.strip() and not line.startswith('##')]
    vcf_rows = [line.split('\t') for line in lines if not line.startswith('#')]

    if not vcf_rows:
        return pd.DataFrame()

    expected_columns = ['CHROM', 'POS', 'SVTYPE', 'ID', 'REF', 'QUAL', 'FILTER', 'INFO']
    vcf_data = pd.DataFrame(vcf_rows, columns=expected_columns[:len(vcf_rows[0])])
    vcf_data['SVTYPE'] = vcf_data['INFO'].str.extract(r'SVTYPE=([^;]+)')
    for svtype in vcf_data['SVTYPE'].dropna().unique():
        subset = vcf_data[vcf_data['SVTYPE'] == svtype]
        file_path = os.path.join(METRIC_OUTPUT_DIR, f'{svtype}_variants.csv')
        subset.to_csv(file_path, index=False)
    return vcf_data

def parse_caller_vcf(file_like_obj):
    lines = [line for line in file_like_obj if not line.startswith('##')]
    vcf_rows = [line.strip().split('\t') for line in lines if not line.startswith('#')]

    if not vcf_rows:
        return pd.DataFrame()

    column_count = len(vcf_rows[0])
    if column_count == 9:
        vcf_data = pd.DataFrame(vcf_rows, columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'SAMPLE'])
    elif column_count == 10:
        vcf_data = pd.DataFrame(vcf_rows, columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'])
    elif column_count == 8:
        vcf_data = pd.DataFrame(vcf_rows, columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])
    else:
        raise ValueError(f"Unexpected number of columns ({column_count}) in Caller VCF")

    vcf_data['SVTYPE'] = vcf_data['INFO'].str.extract(r'SVTYPE=([^;]+)')
    for svtype in vcf_data['SVTYPE'].dropna().unique():
        subset = vcf_data[vcf_data['SVTYPE'] == svtype]
        file_path = os.path.join(METRIC_OUTPUT_DIR, f'{svtype}_caller_variants_output.csv')
        subset.to_csv(file_path, index=False)
    return vcf_data

def parse_intersection_vcf(file_like_obj):
    lines = [line for line in file_like_obj if not line.startswith('##')]
    vcf_rows = [line.strip().split('\t') for line in lines if not line.startswith('#')]

    if vcf_rows:
        column_count = len(vcf_rows[0])
        if column_count == 11:
            vcf_data = pd.DataFrame(vcf_rows, columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'CALLER', 'SAMPLE'])
        else:
            raise ValueError(f"Unexpected number of columns ({column_count}) in Intersection VCF")

        vcf_data['SVTYPE'] = vcf_data['INFO'].str.extract(r'SVTYPE=([^;]+)')
        for svtype in vcf_data['SVTYPE'].dropna().unique():
            subset = vcf_data[vcf_data['SVTYPE'] == svtype]
            file_path = os.path.join(METRIC_OUTPUT_DIR, f'{svtype}_survivor_variants.csv')
            subset.to_csv(file_path, index=False)
        return vcf_data
    else:
        return pd.DataFrame()

# --- Calculate Metrics ---
def calculate_metrics_explicit(caller_df, intersection_df, ground_truth_df):
    tp_df = intersection_df[['CHROM', 'POS']].dropna()
    caller_df = caller_df[['CHROM', 'POS']].dropna()
    ground_truth_df = ground_truth_df[['CHROM', 'POS']].dropna()

    tp_df['CHROM'] = tp_df['CHROM'].astype(str)
    caller_df['CHROM'] = caller_df['CHROM'].astype(str)
    ground_truth_df['CHROM'] = ground_truth_df['CHROM'].astype(str)

    tp_df['POS'] = tp_df['POS'].astype(int)
    caller_df['POS'] = caller_df['POS'].astype(int)
    ground_truth_df['POS'] = ground_truth_df['POS'].astype(int)

    tp_count = tp_df.shape[0]
    fp_df = caller_df.merge(ground_truth_df, on=["CHROM", "POS"], how='left', indicator=True)
    fp_count = fp_df[fp_df['_merge'] == 'left_only'].shape[0]
    fn_df = ground_truth_df.merge(tp_df, on=["CHROM", "POS"], how='left', indicator=True)
    fn_count = fn_df[fn_df['_merge'] == 'left_only'].shape[0]

    precision = tp_count / (tp_count + fp_count) if (tp_count + fp_count) > 0 else 0
    recall = tp_count / (tp_count + fn_count) if (tp_count + fn_count) > 0 else 0
    f1_score = (2 * precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    return precision, recall, f1_score, tp_count, fp_count, fn_count
def get_survivor_upload_section():
    return html.Div([
        html.Label("Select Reference File Type:", style={'fontFamily': '"Times New Roman", Times, serif'}),
        dcc.Dropdown(
            id='reference-type',
            options=[
                {'label': 'Standard VCF', 'value': 'normal'},
                {'label': 'NA12878 EvalSVCallers', 'value': 'na12878'}
            ],
            value='normal',
            style={'width': '100%', 'fontFamily': '"Times New Roman", Times, serif'}
        ),  
        dcc.Upload(
            id='upload-caller',
            children=html.Div(['📎 Drag and Drop or ', html.A('Select Caller VCF File')]),
            style={
                'width': '400px',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'marginTop': '20px',
                'margin': '10px',
                'fontFamily': '"Times New Roman", Times, serif'
            },
            multiple=False
        ),
        html.Div(id='caller-upload-status', style={'marginTop': '5px'}),
        dcc.Upload(
            id='upload-reference',
            children=html.Div(['📎 Drag and Drop or ', html.A('Select Reference VCF File')]),
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
        html.Div(id='reference-upload-status', style={'marginTop': '5px'}),
        dcc.Upload(
            id='upload-intersection',
            children=html.Div(['📎 Drag and Drop or ', html.A('Select Intersection VCF File')]),
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
        html.Div(id='intersection-upload-status', style={'marginTop': '5px'}),
        dbc.Button('Process and Calculate Metrics', id='process-button', n_clicks=0, color="primary", className="mt-2", style={'marginTop': '20px', 'fontFamily': '"Times New Roman", Times, serif'})
    ])


def get_truvari_upload_section():
    return html.Div([
        html.Label("Upload Truvari Summary Output File:", style={'fontFamily': '"Times New Roman", Times, serif'}),
        dcc.Upload(
            id='upload-tru-file',
            children=html.Div(['📎 Drag and Drop or ', html.A('Select *summary File')]),
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
        html.Div(id='upload-tru-txt-status', style={'marginLeft': '10px', 'fontFamily': '"Times New Roman", Times, serif'}),
        dbc.Button('Process Truvari File', id='process-tru-button', n_clicks=0, color="primary", className="mt-2", style={'marginTop': '10px', 'fontFamily': '"Times New Roman", Times, serif'})
    ])





def get_evalsvcallers_upload_section():
    return html.Div([
        html.Label("Upload EvalSVCallers .eval.txt Output File:", style={'fontFamily': '"Times New Roman", Times, serif'}),
        dcc.Upload(
            id='upload-eval-file',
            children=html.Div(['📎 Drag and Drop or ', html.A('Select *.eval.txt File')]),
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
        html.Div(id='upload-eval-txt-status', style={'marginLeft': '10px', 'fontFamily': '"Times New Roman", Times, serif'}),
        dbc.Button('Process EvalSVCallers File', id='process-eval-button', n_clicks=0, color="primary", className="mt-2", style={'marginTop': '10px', 'fontFamily': '"Times New Roman", Times, serif'})
    ])

READ_COUNTS = ["2", "3", "4", "5", "6", "7", "8", "9", "10", "12"]



def parse_truvari_file(contents):
    """
    Accepts a Dash upload 'contents' string for a Truvari summary file
    (JSON-formatted text, regardless of .txt or .json extension).
    Returns a DataFrame in LONG form with columns: ['Metric', 'Value'].
    """
    content_type, content_string = contents.split(',', 1)
    decoded = base64.b64decode(content_string).decode('utf-8', errors='replace').strip()

    # Robust JSON load (handles txt/json the same)
    data = json.loads(decoded)  # will raise if not valid JSON

    # Convert to long form table for easy plotting
    df_long = pd.DataFrame(list(data.items()), columns=["Metric", "Value"])
    return df_long

def generate_truvari_visuals(df):
    """
    Generates:
      1) Summary table
      2) Metric bar: precision, recall, f1 (+ gt_concordance çizgisi)
      3) Count bar: FP, FN, base count, comp count
      4) TP breakdown: TP-base, TP-comp, TP-comp_TP-gt, TP-comp_FP-gt, TP-base_TP-gt, TP-base_FP-gt
    """
    import numpy as np
    import pandas as pd
    import plotly.express as px
    from dash import html, dcc
    import dash_table

    # Nested yapıları at
    df = df[~df["Value"].apply(lambda x: isinstance(x, (dict, list)) or x is None)].copy()

    # Long forma normalize et
    if {"Metric", "Value"}.issubset(df.columns):
        df_long = df.copy()
    else:
        df_long = df.T.reset_index()
        df_long.columns = ["Metric", "Value"]

    # ✅ İSİM NORMALİZASYONU (JSON farklarını gider)
    # - call/comp farkları
    # - TP-call -> TP-comp eşlemesi
    replace_map = {
        "base cnt": "base count",
        "call cnt": "comp count",
        "comp cnt": "comp count",
        "TP-call": "TP-comp",
        "TP-call_TP-gt": "TP-comp_TP-gt",
        "TP-call_FP-gt": "TP-comp_FP-gt",
    }
    df_long["Metric"] = df_long["Metric"].replace(replace_map)

    # Wide tablo üretmek için sözlük
    wide = pd.DataFrame([dict(df_long.values)])

    # Sayısallaştır
    df_long["Value"] = pd.to_numeric(df_long["Value"], errors="coerce")

    # === Summary Table ===
    rate_keys = {"precision", "recall", "f1", "gt_concordance"}
    table_row = {}
    for col, val in wide.iloc[0].items():
        if col in rate_keys and pd.notna(val):
            table_row[col] = f"{float(val):.4f}"
        else:
            if pd.notna(val) and (isinstance(val, (int, np.integer)) or (isinstance(val, float) and float(val).is_integer())):
                table_row[col] = int(val)
            else:
                table_row[col] = val

    table_df = pd.DataFrame([table_row])
    table = dash_table.DataTable(
        data=table_df.to_dict("records"),
        columns=[{"name": c, "id": c} for c in table_df.columns],
        style_table={'overflowX': 'auto'},
        style_cell={'textAlign': 'center', 'fontFamily': '"Times New Roman", Times, serif'},
        style_header={'backgroundColor': 'lightgrey', 'fontWeight': 'bold'},
        page_size=15
    )

    # === Figure 1: precision, recall, f1 (+ gt_concordance çizgisi) ===
    fig1_keys = ["precision", "recall", "f1"]
    fig1_df = df_long[df_long["Metric"].isin(fig1_keys)].copy()
    fig1_df["Value_fmt"] = fig1_df["Value"].map(lambda v: f"{v:.4f}" if pd.notna(v) else "")

    gt_concordance_val = df_long.loc[df_long["Metric"].eq("gt_concordance"), "Value"]
    gt_concordance_val = float(gt_concordance_val.iloc[0]) if len(gt_concordance_val) else None

    fig1 = px.bar(
        fig1_df, x="Metric", y="Value", text="Value_fmt",
        title="Truvari Performance Metrics (precision/recall/f1) + GT Concordance",
        color="Metric", color_discrete_sequence=px.colors.qualitative.Pastel,
        range_y=[0, 1]
    )
    fig1.update_traces(textposition='outside')

    if gt_concordance_val is not None:
        fig1.add_shape(
            type="line",
            x0=-0.5, x1=2.5, y0=gt_concordance_val, y1=gt_concordance_val,
            line=dict(color="black", width=2, dash="dash"),
            name="GT Concordance"
        )
        fig1.add_annotation(
            x=1.5, y=gt_concordance_val,
            text=f"GT Concordance = {gt_concordance_val:.4f}",
            showarrow=False, yshift=10,
            font=dict(size=12, family="Times New Roman")
        )
    fig1.update_layout(height=420, margin=dict(l=20, r=20, t=60, b=40))

    # === Figure 2: Sayımlar — FP, FN, base count, comp count ===
    wanted_order = ["FP", "FN", "base count", "comp count", "call count"]  # call count varsa da göster (eski veri için)
    present = [k for k in wanted_order if k in set(df_long["Metric"])]
    fig2_df = df_long[df_long["Metric"].isin(present)].copy()
    fig2_df["Value_fmt"] = fig2_df["Value"].map(lambda v: f"{v:.0f}" if pd.notna(v) else "")
    fig2_df["Metric"] = pd.Categorical(fig2_df["Metric"], categories=present, ordered=True)
    fig2_df = fig2_df.sort_values("Metric")

    fig2 = px.bar(
        fig2_df, x="Metric", y="Value", text="Value_fmt",
        title="Callset & Error Counts",
        color="Metric", color_discrete_sequence=px.colors.qualitative.Vivid
    )
    fig2.update_traces(textposition='outside')
    fig2.update_layout(height=420, margin=dict(l=20, r=20, t=60, b=40))

    # === Figure 3: TP ayrımı — base/comp ve GT breakdown ===
    fig3_keys_ordered = [
        "TP-base", "TP-comp",
        "TP-base_TP-gt", "TP-comp_TP-gt",
        "TP-base_FP-gt", "TP-comp_FP-gt"
    ]
    fig3_df = df_long[df_long["Metric"].isin(fig3_keys_ordered)].copy()
    fig3_df["Metric"] = pd.Categorical(fig3_df["Metric"], categories=fig3_keys_ordered, ordered=True)
    fig3_df = fig3_df.sort_values("Metric")
    fig3_df["Value_fmt"] = fig3_df["Value"].map(lambda v: f"{v:.0f}" if pd.notna(v) else "")

    fig3 = px.bar(
        fig3_df, x="Metric", y="Value", text="Value_fmt",
        title="True Positive Breakdown (base/comp + GT)",
        color="Metric", color_discrete_sequence=px.colors.qualitative.Set2
    )
    fig3.update_traces(textposition='outside')
    fig3.update_layout(height=420, margin=dict(l=20, r=20, t=60, b=40))

    # === Final Output ===
    return html.Div([
        html.H4("Truvari Summary Metrics", style={'marginTop': '0px'}),
        table,
        html.Br(),
        dcc.Graph(figure=fig1),
        dcc.Graph(figure=fig2),
        dcc.Graph(figure=fig3)
    ])
    
def parse_evalsvcallers_file(contents):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    lines = io.StringIO(decoded.decode("utf-8")).readlines()

    # --- Part 1: Summary Table ---
    summary_lines = lines[:4]
    ref_data = []
    for line in summary_lines:
        parts = line.strip().split("\t")
        ref_type = parts[0].split(":")[0]
        total = parts[0].split(":")[1].strip()
        for part in parts[1:]:
            if ": " in part:
                k, v = part.split(": ")
                ref_data.append({
                    "Ref Variant": ref_type,
                    "Total": total,
                    "Size Category": k.strip(),
                    "Count": v.strip()
                })
    df_ref = pd.DataFrame(ref_data)
    pivot_ref = df_ref.pivot(index=["Ref Variant", "Total"], columns="Size Category", values="Count").reset_index()
    pivot_ref = pivot_ref.apply(pd.to_numeric, errors="ignore")
    if "short (<= 1.0 kp)" in pivot_ref.columns and "short (<=1.0 kp)" in pivot_ref.columns:
        pivot_ref["short (<=1.0 kp)"] = pivot_ref["short (<= 1.0 kp)"].fillna(0).astype(int) + pivot_ref["short (<=1.0 kp)"].fillna(0).astype(int)
        pivot_ref.drop(columns=["short (<= 1.0 kp)"], inplace=True)
    pivot_ref = pivot_ref.rename(columns=lambda c: c.replace(" (", "\n("))

    # --- Part 2: Metric Tables and Long Format ---
    sv_sections = {}
    current_type = None
    for line in lines:
        if line.startswith("##"):
            current_type = line.strip().replace("##", "").strip()
            sv_sections[current_type] = []
        elif current_type and line.strip() and not line.startswith("<<") and not line.startswith("Ref-") and "<Number of" not in line:
            sv_sections[current_type].append(line.strip())

    records = []
    df_blocks = []
    precision_dict, recall_dict = {}, {}

    for svtype, block_lines in sv_sections.items():
        for l in block_lines:
            parts = l.strip().split("\t")
            if not parts or parts[0].strip().replace('.', '', 1).isdigit():
                continue
            metric_raw = parts[0].strip()
            values = parts[1:11]

            if "(" in metric_raw:
                metric = metric_raw.split("(")[0].strip()
                size_class = metric_raw.split("(")[1].replace(")", "").strip()
            else:
                metric, size_class = metric_raw, "All"

            row = {"SVTYPE": svtype, "Metric": metric, "SizeClass": size_class}
            row.update({
                rc: round(pd.to_numeric(val, errors='coerce') / 100, 4) if metric.lower() in ["precis", "recall"]
                else pd.to_numeric(val, errors='coerce')
                for rc, val in zip(READ_COUNTS, values)
            })
            df_blocks.append(row)

            # Record for long-form
            for rc, val in zip(READ_COUNTS, values):
                numeric_val = pd.to_numeric(val, errors='coerce')
                value = numeric_val / 100 if metric.lower() in ["precis", "recall"] else numeric_val
                records.append({
                    "SVTYPE": svtype,
                    "Read Count": rc,
                    "Metric": metric,
                    "SizeClass": size_class,
                    "Value": value
                })
                if metric.lower() == "precis":
                    precision_dict[(svtype, size_class, rc)] = value
                elif metric.lower() == "recall":
                    recall_dict[(svtype, size_class, rc)] = value
    # Step 1: Accumulate all read counts in one row per SVTYPE and SizeClass
    f1_rows = {}
    # F1 calculation
    for key in precision_dict:
        if key in recall_dict:
            p, r = precision_dict[key], recall_dict[key]
            f1 = round(2 * p * r / (p + r), 3) if p + r > 0 else 0
            svtype, size_class, rc = key
            records.append({
                "SVTYPE": svtype,
                "Read Count": rc,
                "Metric": "F1-Score",
                "SizeClass": size_class,
                "Value": f1
            })
        # Wide-form: accumulate
        row_key = (svtype, size_class)
        if row_key not in f1_rows:
            f1_rows[row_key] = {
                "SVTYPE": svtype,
                "Metric": "F1-Score",
                "SizeClass": size_class
            }
        f1_rows[row_key][rc] = f1

    # Step 2: Add combined rows to df_blocks
    for row in f1_rows.values():
        df_blocks.append(row)
    
    return pivot_ref, pd.DataFrame(records), pd.DataFrame(df_blocks)

def generate_evalsvcallers_visuals(pivot_ref, df_long, df_block):
    import dash_table
    import plotly.express as px
    from dash import html, dcc
    

    # Define metric order for main plot (excluding 'Call')
    metric_order = ["Recall", "Precis", "F1-Score"]

    # === Filter dataframes ===
    df_metrics_main = df_block[df_block["Metric"] != "Call"].copy()
    df_metrics_main["Metric_Order"] = df_metrics_main["Metric"].apply(
        lambda m: metric_order.index(m) if m in metric_order else 999
    )
    df_metrics_main = df_metrics_main.sort_values(by=["SVTYPE", "SizeClass", "Metric_Order"]).drop(columns=["Metric_Order"])

    df_long_main = df_long[df_long["Metric"] != "Call"]
    df_long_call = df_long[df_long["Metric"] == "Call"]


     # ✅ Fix newline and typo in column names
    pivot_ref.columns = (
        pivot_ref.columns
        .str.replace("\n", " ", regex=False)
        .str.strip()
        .str.replace("(<=1.0", "(<= 1.0", regex=False)  # ensure space after <=
    )

    pivot_ref = pivot_ref.rename(columns={"tinny (<= 100 bp)": "tiny (<= 100 bp)"})

    # Melt after renaming
    df_melted = pivot_ref.melt(
        id_vars=["Ref Variant"],
        var_name="Size Category",
        value_name="Count"
    )

    static_colors = {
        "A": "rgb(204, 97, 176)",           # orange
        "SS": "rgb(82, 188, 163)",         # blue
        "S": "rgb(47, 138, 196)",          # teal
        "M": "rgb(237, 100, 90)",          # green
        "L": "rgb(93, 105, 177)",          # pink
    
        "Total": "rgb(204, 97, 176)",       # matches A
        "tiny (<= 100 bp)": "rgb(82, 188, 163)",  # matches SS
        "short (<= 1.0 kp)": "rgb(47, 138, 196)", # matches S
        "middle (<= 100 kb)": "rgb(237, 100, 90)",# matches M
        "large (> 100 kb)": "rgb(93, 105, 177)"   # matches L
    }
    
    # Use in plots:
    size_color_discrete_map = {k: static_colors[k] for k in ["A", "SS", "S", "M", "L"]}
    tsml_color_discrete_map = {k: static_colors[k] for k in ["Total", "tiny (<= 100 bp)", "short (<= 1.0 kp)", "middle (<= 100 kb)", "large (> 100 kb)"]}
    
    
    # === Reference stacked bar plot ===
    bar_fig = px.bar(
        df_melted,
        x="Ref Variant", y="Count", color="Size Category", color_discrete_map=tsml_color_discrete_map, barmode="stack",
        title="Reference Variant Distribution"
    )
    # === Line plot for main metrics ===
    line_fig_main = px.line(
        df_long_main,
        x="Read Count",
        y="Value",
        color="SizeClass",             # now coloring by SizeClass
        color_discrete_map=size_color_discrete_map, 
        line_dash="SizeClass",         # still keep line style
        facet_row="SVTYPE",            # separate by SV type (DEL, DUP, INS)
        facet_col="Metric",            # separate by metric (Precision, Recall, F1)
        markers=True,
        title="Metrics by SV Type, Size Class and Read Count",
        category_orders={
            "SizeClass": ["A", "SS", "S", "M", "L"],
            "Metric": ["Precis","Recall","F1-Score"]
        }
      #  color_discrete_map={
      #      "A": "#1f77b4", "SS": "#ff7f0e", "S": "#2ca02c", "M": "#d62728", "L": "#9467bd"
      #  }
    )

    #line_fig_main.update_traces(marker=dict(size=8))  # make markers larger
    line_fig_main.update_layout(height=800)
    line_fig_main.update_yaxes(range=[0, 1])
    # === Separate line plot for 'Call' ===
    line_fig_call = px.line(
        df_long_call, x="Read Count", y="Value",color="SizeClass", color_discrete_map=size_color_discrete_map, 
        line_dash="SizeClass",  facet_row="SVTYPE",
        markers=True, title="Variant Count (Call) by SV Type and Size Class"
    )
    line_fig_call.update_layout(height=800)
    # Correct the column name typo
    
    # === Return all visuals ===
    return html.Div([
        html.H4("Reference Summary Table", style={'marginTop': '0px', 'fontFamily': '"Times New Roman", Times, serif', 'fontSize': '20px', 'fontWeight': 'bold'}),
        dash_table.DataTable(
            data=pivot_ref.to_dict("records"),
            columns=[{"name": col, "id": col} for col in pivot_ref.columns],
            style_table={'overflowX': 'auto'},
            page_size=10,
            style_cell={'textAlign': 'center'}
        ),
        html.Br(),
        dcc.Graph(figure=bar_fig),
        html.H4("Read Support Metrics Table (All SV Types)", style={'marginTop': '0px', 'fontFamily': '"Times New Roman", Times, serif', 'fontSize': '20px', 'fontWeight': 'bold'}),
        dash_table.DataTable(
            data=df_metrics_main.to_dict("records"),
            columns=[{"name": col, "id": col} for col in df_metrics_main.columns],
            style_table={'overflowX': 'auto'},
            page_size=20,
            style_cell={'textAlign': 'center'}
        ),
        html.Br(),
        dcc.Graph(figure=line_fig_main),
    #    html.H4("Variant Count (Call) Plot", style={'marginTop': '20px', 'fontFamily': '"Times New Roman", Times, serif', 'fontSize': '20px', 'fontWeight': 'bold'}),
        dcc.Graph(figure=line_fig_call)
    ])
def process_survivor_metrics_from_content(ref_content, caller_content, inter_content, reference_type): 

    # Parse reference
    if reference_type == 'normal':
        ref_vcf = parse_reference_normal(parse_contents(ref_content))
    else:
        ref_vcf = parse_reference_na12878(parse_contents(ref_content))

    # Parse caller and intersection
    caller_vcf = parse_caller_vcf(parse_contents(caller_content))
    inter_vcf = parse_intersection_vcf(parse_contents(inter_content))

    # Group by SVTYPE
    ref_svtypes = {svtype: ref_vcf[ref_vcf['SVTYPE'] == svtype] for svtype in ref_vcf['SVTYPE'].dropna().unique()}
    caller_svtypes = {svtype: caller_vcf[caller_vcf['SVTYPE'] == svtype] for svtype in caller_vcf['SVTYPE'].dropna().unique()}
    inter_svtypes = {svtype: inter_vcf[inter_vcf['SVTYPE'] == svtype] for svtype in inter_vcf['SVTYPE'].dropna().unique()}

    common_svtypes = set(ref_svtypes.keys()) & set(caller_svtypes.keys()) & set(inter_svtypes.keys())

    if not common_svtypes:
        return None, "⚠️ No common SVTYPEs found across the three files."

    # Compute metrics
    metrics = []
    for svtype in sorted(common_svtypes):
        ref_df = ref_svtypes[svtype]
        caller_df = caller_svtypes[svtype]
        inter_df = inter_svtypes[svtype]

        precision, recall, f1, tp, fp, fn = calculate_metrics_explicit(caller_df, inter_df, ref_df)
        metrics.append({
            'SVTYPE': svtype,
            'Precision': round(precision, 4),
            'Recall': round(recall, 4),
            'F1-Score': round(f1, 4),
            'TP': tp,
            'FP': fp,
            'FN': fn
        })

    return pd.DataFrame(metrics), None

def generate_survivor_visuals(metrics_df):
    import plotly.express as px
    from dash import html, dcc
    import dash_table

    return html.Div([
        dash_table.DataTable(
            columns=[{"name": col, "id": col} for col in metrics_df.columns],
            data=metrics_df.to_dict('records'),
            style_cell={'textAlign': 'center'},
            style_header={'backgroundColor': 'lightgrey', 'fontWeight': 'bold'},
            page_size=20
        ),
        html.Br(),
        dcc.Graph(
            figure=px.bar(
                metrics_df.melt(id_vars='SVTYPE', value_vars=['Precision', 'Recall', 'F1-Score'],
                                var_name='Metric', value_name='Score'),
                x='SVTYPE',
                y='Score',
                color='Metric', 
                color_discrete_map=metric_color_map, 
                barmode='group',
                title='Precision, Recall, and F1 Score per SVTYPE',
                height=500
            ).update_yaxes(range=[0, 1])
        )
    ])