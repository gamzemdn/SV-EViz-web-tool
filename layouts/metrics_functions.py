import os
import io
import base64
import pandas as pd
import plotly.express as px
from collections import defaultdict
from dash import html, dcc
import dash_bootstrap_components as dbc
UPLOAD_DIRECTORY = "./uploaded_files/"
METRIC_OUTPUT_DIR = os.path.join(UPLOAD_DIRECTORY, "metrics_output")
os.makedirs(UPLOAD_DIRECTORY, exist_ok=True)
os.makedirs(METRIC_OUTPUT_DIR, exist_ok=True)
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

def get_evalsvcallers_upload_section():
    return html.Div([
        html.Label("Upload EvalSVCallers .eval.txt Output File:", style={'fontFamily': '"Times New Roman", Times, serif'}),
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
        html.Div(id='upload-eval-txt-status', style={'marginLeft': '10px', 'fontFamily': '"Times New Roman", Times, serif'}),
        dbc.Button('Process EvalSVCallers File', id='process-eval-button', n_clicks=0, color="primary", className="mt-2", style={'marginTop': '10px', 'fontFamily': '"Times New Roman", Times, serif'})
    ])

READ_COUNTS = ["2", "3", "4", "5", "6", "7", "8", "9", "10", "12"]

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

    # === Reference stacked bar plot ===
    bar_fig = px.bar(
        pivot_ref.melt(id_vars=["Ref Variant", "Total"], var_name="Size Category", value_name="Count"),
        x="Ref Variant", y="Count", color="Size Category", barmode="stack",
        title="Reference Variant Distribution"
    )

    # === Line plot for main metrics ===
    line_fig_main = px.line(
        df_long_main,
        x="Read Count",
        y="Value",
        color="SizeClass",             # now coloring by SizeClass
        line_dash="SizeClass",         # still keep line style
        facet_row="SVTYPE",            # separate by SV type (DEL, DUP, INS)
        facet_col="Metric",            # separate by metric (Precision, Recall, F1)
        markers=True,
        title="Metrics by SV Type, Size Class and Read Count",
        category_orders={
            "SizeClass": ["A", "SS", "S", "M", "L"],
            "Metric": ["Precis","Recall","F1-Score"]
        },
        color_discrete_map={
            "A": "#1f77b4", "SS": "#ff7f0e", "S": "#2ca02c", "M": "#d62728", "L": "#9467bd"
        }
    )

    #line_fig_main.update_traces(marker=dict(size=8))  # make markers larger
    line_fig_main.update_layout(height=800)
    line_fig_main.update_yaxes(range=[0, 1])
    # === Separate line plot for 'Call' ===
    line_fig_call = px.line(
        df_long_call, x="Read Count", y="Value",color="SizeClass", 
        line_dash="SizeClass",  facet_row="SVTYPE",
        markers=True, title="Variant Count (Call) by SV Type and Size Class"
    )
    line_fig_call.update_layout(height=800)

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
                barmode='group',
                title='Precision, Recall, and F1 Score per SVTYPE',
                height=500
            ).update_yaxes(range=[0, 1])
        )
    ])