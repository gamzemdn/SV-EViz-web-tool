import os
import base64
import datetime
import json
import io
from collections import defaultdict
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import dash_bio as dashbio
from dash import dcc, html, dash_table
import flask
import re
import numpy as np
from dash_bio import Clustergram
# Directories
UPLOAD_DIRECTORY = "./uploaded_files/"
VISUALIZATION_OUTPUT_DIR = os.path.join(UPLOAD_DIRECTORY, "visualization_output")
os.makedirs(VISUALIZATION_OUTPUT_DIR, exist_ok=True)
svtypes = ['DEL', 'DUP', 'INS', 'INV', 'BND', 'CNV']  # or whatever types you have
base_colors = px.colors.qualitative.Vivid.copy()  # or px.colors.qualitative.Bold
# Optionally remove yellowish tone (like '#FECB52')
exclude_colors = ['rgb(218, 165, 27)']  # Plotly yellow
color_palette = [c for c in base_colors if c not in exclude_colors]
svtype_color_map = {
    sv: color_palette[i % len(color_palette)] for i, sv in enumerate(svtypes)
}

# --- Genome length dictionaries ---
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
# Colors
chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
chrom_color_map = dict(zip(chrom_order, chrom_colors))

svtype_colors = {
    "DEL": "#e74c3c", "INS": "#3498db", "DUP": "#2ecc71",
    "INV": "#f1c40f", "BND": "#9b59b6"
}

# ---- File Saving ----
def save_file(name, content):
    data = content.encode("utf8").split(b";base64,")[1]
    file_path = os.path.join(VISUALIZATION_OUTPUT_DIR, name)
    with open(file_path, "wb") as fp:
        fp.write(base64.b64decode(data))
    return file_path
def get_variant_extraction_status():
    return html.Li(
        "Variant types extracted and dropdown updated successfully.",
        style={'fontFamily': '"Times New Roman", Times, serif', 'marginRight': '20px'}
    )    
# ---- Data Parsing ----
def parse_uploaded_vcf(contents, filename, date):
    file_path = save_file(filename, contents)
    return html.Li(
    f"{filename} (Uploaded at {datetime.datetime.fromtimestamp(date)})",
    style={'fontFamily': '"Times New Roman", Times, serif', 'marginRight': '20px'}
), file_path

# ---- Data Extraction ----
def extract_variant_types(file_path, source_type):
    if not file_path or not os.path.exists(file_path):
        return [], ""

    try:
        if source_type == 'caller':
            try:
                df = pd.read_csv(file_path, sep="\t", comment="#", header=None)
                if df.shape[1] >= 8:
                    df = df.iloc[:, :8]
                    df.columns = ["CHROM", "POS", "SVTYPE", "col4", "col5", "QUAL", "FILTER", "INFO"]
                else:
                    raise ValueError("Caller file must have at least 8 columns.")
            except Exception as e:
                print(f"❌ Error loading caller file: {e}")
                return [], ""

            df["SVTYPE"] = df["INFO"].str.extract(r"SVTYPE=([^;]+)")
            df = df[df["SVTYPE"] != "BND"]  # ✅ Filter out BND
            df["CHROM"] = df["CHROM"].astype(str).str.replace("^chr", "", regex=True)
            svtypes = df["SVTYPE"].dropna().unique().tolist()

        elif source_type == 'survivor':
            df = pd.read_csv(file_path, sep="\t", comment="#", header=None, usecols=range(8))
            df.columns = ["CHROM", "POS", "SVTYPE", "col4", "col5", "QUAL", "FILTER", "INFO"]
            df["SVTYPE"] = df["INFO"].str.extract(r"SVTYPE=([^;]+)")
            df["CHROM"] = df["CHROM"].astype(str).str.replace("^chr", "", regex=True)
            svtypes = df["SVTYPE"].dropna().unique().tolist()

        elif source_type == 'evalsvcallers':
            df = pd.read_csv(file_path, sep="\t", header=None)
            tf_chrom = df[0].astype(str).str.strip().str.split(" ", expand=True)
            df.insert(0, "TF", tf_chrom[0])
            df.insert(1, "CHROM", tf_chrom[1])
            df = df.drop(columns=[0])
            df = df[df["TF"] == "TP"]
            df.columns = ["TF", "CHROM", "POS", "SVTYPE", "col4", "col5", "QUAL", "FILTER", "INFO"]
            df["CHROM"] = df["CHROM"].astype(str).str.replace("^chr", "", regex=True)
            svtypes = df["SVTYPE"].dropna().unique().tolist()
            df = df.drop(columns=["TF"])
        else:
            return [], ""

        svtypes = sorted(set(svtypes))
        svtypes.insert(0, "ALL")
        return [{"label": sv, "value": sv} for sv in svtypes], df.to_json(date_format='iso')

    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        return [], ""

# ---- Plotting ----
def plot_vcf_data(df):
    figures = []
    if "INFO" in df.columns:
        df["SVTYPE"] = df["INFO"].str.extract(r'SVTYPE=([^;]+)')
        if df["SVTYPE"].notna().any():
            sv_counts = df["SVTYPE"].value_counts().reset_index()
            sv_counts.columns = ["Variant Type", "Count"]
            fig_svtype = px.bar(
                sv_counts, x="Variant Type", y="Count", color="Variant Type", color_discrete_map=svtype_color_map,
                title="Variant Type Distribution"
            )
            figures.append(dcc.Graph(figure=fig_svtype))

            fig_spyder = go.Figure()
            fig_spyder.add_trace(go.Scatterpolar(
                r=sv_counts["Count"].tolist() + [sv_counts["Count"].tolist()[0]],
                theta=sv_counts["Variant Type"].tolist() + [sv_counts["Variant Type"].tolist()[0]],
                fill='toself', name='Variant Types'
            ))
            fig_spyder.update_layout(
                polar=dict(radialaxis=dict(visible=True)),
                title="Variant Type Radar Chart",
                showlegend=False
            )
            figures.append(dcc.Graph(figure=fig_spyder))

    if "#CHROM" in df.columns:
        chrom_counts = df["#CHROM"].value_counts().reset_index()
        chrom_counts.columns = ["Chromosome", "Variant Count"]
        fig_chr = px.bar(
            chrom_counts, x="Chromosome", y="Variant Count",
            title="Chromosome-Wise Variant Distribution"
        )
        figures.append(dcc.Graph(figure=fig_chr))

    if "INFO" in df.columns:
        df["SVLEN"] = df["INFO"].str.extract(r'SVLEN=(-?\d+)')
        df["SVLEN"] = pd.to_numeric(df["SVLEN"], errors='coerce')
        if df["SVLEN"].notna().any():
            fig_svlen = px.histogram(
                df, x="SVLEN", nbins=30,
                title="Variant Length Distribution"
            )
            figures.append(dcc.Graph(figure=fig_svlen))

          # ✅ Violin plot (log-transformed absolute SVLEN)
            df["SVLEN_abs"] = df["SVLEN"].abs()
            df["logSVLEN"] = np.log10(df["SVLEN_abs"] + 1)
            violin_df = df[df["SVTYPE"].notna() & df["logSVLEN"].notna()]

            if not violin_df.empty:
                fig_violin = px.violin(
                    violin_df, x="SVTYPE", y="logSVLEN", color="SVTYPE", color_discrete_map=svtype_color_map, box=True, points="all",
                    title="Log-transformed SV Length by Variant Type",
                    labels={"logSVLEN": "log10(SVLEN + 1)"}
                )
                figures.append(dcc.Graph(figure=fig_violin))
    return figures

# ---- Circos JSON Builder ----
def detect_genome_version(vcf_path):
    with open(vcf_path) as f:
        for line in f:
            if not line.startswith("##") and not line.startswith("#CHROM"):
                chrom = line.split("\t")[0]
                return "hg38" if chrom.startswith("chr") else "hg19"
    return "hg38"
def vcf_to_circos_json(vcf_file_path, output_json_path, source_type, bin_size=10_000_000):
    genome = detect_genome_version(vcf_file_path)
    chromosome_lengths = hg38_lengths if genome == "hg38" else hg19_lengths
    chromosome_lengths = {chrom: max(1, length) for chrom, length in chromosome_lengths.items()}
    valid_chroms = set(chromosome_lengths.keys())

    layout_chromosomes = set()
    raw_chords = []
    chord_metadata = []
    tracks = []
    histograms = defaultdict(lambda: defaultdict(int))
    svtypes_found = set()
    
    with open(vcf_file_path) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue

            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue

            chrom = parts[0]
            try:
                pos = int(parts[1])
            except:
                continue

            # Defaults
            svtype = "UNK"
            end = None
            svlen = None
            qual = 0
            reads = 0
            gt = ""

            if source_type in ["caller", "survivor"]:
                info_field = parts[7]
                info = {}
                for kv in info_field.split(";"):
                    if "=" in kv:
                        k, v = kv.split("=", 1)
                        info[k] = v
                svtype = info.get("SVTYPE", "UNK").strip().upper()
                if svtype == "BND":
                    continue

                try:
                    qual = int(float(parts[5])) if parts[5].replace('.', '', 1).isdigit() else 0
                except:
                    qual = 0
                    
                end_str = info.get("END", "")
                svlen_str = info.get("SVLEN", "")

                if svlen_str:
                    try:
                        svlen = abs(int(svlen_str))
                    except:
                        svlen = None

                if end_str:
                    try:
                        end = int(end_str)
                    except:
                        end = None

            elif source_type == "evalsvcallers":
                tf_chrom = parts[0].strip().split(" ")
                if len(tf_chrom) != 2 or tf_chrom[0] != "TP":
                    continue
                chrom = tf_chrom[1]

                try:
                    pos = int(parts[1])
                    svtype = parts[2].strip().upper()
                    if svtype == "BND":
                        continue

                    end = int(parts[3]) if parts[3].isdigit() else None
                    svlen = abs(int(parts[4])) if parts[4].lstrip("-").isdigit() else None
                    try:
                        qual = int(float(parts[5])) if parts[5].replace('.', '', 1).isdigit() else 0
                    except:
                        qual = 0
                except Exception as ex:
                    print(f"❌ Exception parsing line: {parts} → {ex}")
                    continue

            # Final fallback logic
            if end is None and svlen is not None:
                end = pos + abs(svlen)
            if svlen is None and end is not None:
                svlen = abs(end - pos)
            if end is None and svlen is None:
                end = pos + 1000
                svlen = 1000

            chrom = f"chr{chrom}" if not chrom.startswith("chr") else chrom
            if chrom not in valid_chroms:
                continue

            layout_chromosomes.add(chrom)
            svtypes_found.add(svtype)

            bin_start = (pos // bin_size) * bin_size
            histograms[svtype][(chrom, bin_start)] += 1

            if end != pos:
                source = {"id": chrom, "start": pos, "end": end}
                target = {"id": chrom, "start": end, "end": end + 1}
                color = svtype_color_map.get(svtype, "#95a5a6")

                raw_chords.append({
                    "source": source,
                    "target": target,
                    "color": color,
                    "value": svlen
                })

                chord_metadata.append({
                    "source": source,
                    "target": target,
                    "svtype": svtype,
                    "gt": gt,
                    "reads": reads,
                    "qual": qual,
                    "svlen": svlen
                })

    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    layout = [
        {
            "id": chrom,
            "label": chrom,
            "color": chrom_color_map.get(chrom, "#bbbbbb"),
            "len": chromosome_lengths[chrom],
        }
        for chrom in chrom_order if chrom in layout_chromosomes
    ]

    for svtype, bins in histograms.items():
        data = [
            {
                "block_id": chrom,
                "start": start,
                "end": start + bin_size,
                "value": count
            }
            for (chrom, start), count in bins.items() if count > 0
        ]
        if data:
            tracks.append({
                "type": "histogram",
                "data": data,
                "name": svtype,
                "color": svtype_color_map.get(svtype, "#95a5a6"),
                "height": 0.08
            })

    circos_json = {
        "layout": layout,
        "chords": raw_chords,
        "chord_metadata": chord_metadata,
        "tracks": tracks,
        "detected_svtypes": sorted(svtypes_found)
    }

    with open(output_json_path, "w") as f:
        json.dump(circos_json, f, indent=2)

def plot_sankey(df):
 
    df["CHROM"] = df["CHROM"].astype(str).str.replace("^chr", "", regex=True)
    chrom_svtype_counts = df.groupby(["CHROM", "SVTYPE"]).size().reset_index(name="count")

    chrom_order = [str(i) for i in range(1, 23)] + ["X", "Y"]
    chrom_svtype_counts = chrom_svtype_counts[chrom_svtype_counts["CHROM"].isin(chrom_order)]
    chrom_svtype_counts["Chromosome_Sort"] = chrom_svtype_counts["CHROM"].apply(lambda v: chrom_order.index(v))
    chrom_svtype_counts = chrom_svtype_counts.sort_values("Chromosome_Sort")
    chrom_svtype_counts["CHROM"] = "chr" + chrom_svtype_counts["CHROM"].astype(str)

    chrom_order_sorted = chrom_svtype_counts["CHROM"].unique().tolist()
    svtypes_sorted = chrom_svtype_counts["SVTYPE"].unique().tolist()
    labels = chrom_order_sorted + svtypes_sorted
    label_to_idx = {label: idx for idx, label in enumerate(labels)}

    sources = chrom_svtype_counts["CHROM"].map(label_to_idx)
    targets = chrom_svtype_counts["SVTYPE"].map(label_to_idx)
    values = chrom_svtype_counts["count"]
    
    node_colors = (
        ['#DC143C'] * len(chrom_order_sorted) +  # chromosomes (default light blue)
        [svtype_color_map.get(sv, '#999999') for sv in svtypes_sorted]  # svtypes
    )
    # 🎨 Define link colors based on target SVTYPE
    target_colors = chrom_svtype_counts["SVTYPE"].map(
        lambda sv: svtype_color_map.get(sv, 'rgba(160,160,160,0.4)')
    )
    link_colors = target_colors.apply(lambda c: c.replace("rgb(", "rgba(").replace(")", ",0.5)") if c.startswith("rgb(") else c)
    
    fig = go.Figure(go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labels,
            color=node_colors
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color=link_colors
        )
    ))
    fig.update_layout(title_text="Chromosome to SVTYPE Sankey Plot", font_size=14)
    legend = html.Div([
            html.Div([
                html.Span(style={
                    'backgroundColor': svtype_color_map.get(sv, '#999'),
                    'display': 'inline-block',
                    'width': '15px',
                    'height': '15px',
                    'marginRight': '5px',
                    'borderRadius': '3px'
                }),
                html.Span(f"{sv}", style={'marginRight': '15px', 'fontSize': '14px'})
            ], style={'display': 'inline-block'}) for sv in svtypes_sorted
        ], style={'padding': '10px 0', 'textAlign': 'left'})


    return html.Div([
        legend,
        dcc.Graph(figure=fig, style={"height": "700px", "width": "100%"})
    ])
def plot_clustergram(df, selected_chroms=None):
    # ✅ Always re-extract SVTYPE from INFO
    df["SVTYPE"] = df["INFO"].str.extract(r"SVTYPE=([^;]+)")
    df = df[df["SVTYPE"].notnull()]
    df = df[df["SVTYPE"] != "BND"]

    # ✅ Normalize CHROM
    df["CHROM"] = df["CHROM"].astype(str).str.strip()
    valid_chroms = [str(i) for i in range(1, 23)] + ["X", "Y"]
    df["CHROM"] = df["CHROM"].apply(
        lambda x: f"chr{x}" if not x.startswith("chr") and x in valid_chroms
        else x if x.startswith("chr") and x[3:] in valid_chroms
        else None
    )
    df = df.dropna(subset=["CHROM", "SVTYPE"])

    # ✅ Pivot only on valid CHROM–SVTYPE
    pivot_df = df.groupby(["CHROM", "SVTYPE"]).size().unstack(fill_value=0)
    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    pivot_df = pivot_df.reindex(chrom_order).dropna(how="all")

    if pivot_df.empty:
        return html.Div("⚠️ No data available for Clustergram.")

    if selected_chroms is None:
        selected_chroms = [c for c in chrom_order if c in pivot_df.index]
    pivot_df = pivot_df.loc[selected_chroms]

    # ✅ Remove any zero-only rows or cols
    pivot_df = pivot_df.loc[(pivot_df != 0).any(axis=1)]
    pivot_df = pivot_df.loc[:, (pivot_df != 0).any(axis=0)]

    if pivot_df.shape[0] < 2 or pivot_df.shape[1] < 2:
        return html.Div("⚠️ Not enough informative data to generate Clustergram.")

    # 1) Figure'ı üret
    fig = Clustergram(
        data=pivot_df.values,
        row_labels=pivot_df.index.tolist(),
        column_labels=pivot_df.columns.tolist(),
        color_map="RdBu",
        color_threshold={'row': 250, 'col': 700},
        height=800,
        width=900,
        display_ratio=[0.85, 0.15],
        # hidden_labels=[]  # gerekirse aç
    )
    
    # 2) Tüm x/y eksenlerinde tick label’ları AÇ
    layout_dict = fig.layout.to_plotly_json()
    
    axis_updates = {}
    for key in layout_dict.keys():
        if str(key).startswith("xaxis") or str(key).startswith("yaxis"):
            axis_updates[key] = dict(showticklabels=True)
    
    fig.update_layout(
        **axis_updates,
        margin=dict(l=110, r=60, t=50, b=110)  # etiketler için alan
    )
    
    return dcc.Graph(figure=fig)


def plot_manhattan(df, selected_svtypes, threshold):
    import numpy as np
    import plotly.graph_objects as go
    from dash import dcc

    highlight_df_list = []

    df = df[df["SVTYPE"].isin(selected_svtypes)]

    # ✅ Normalize CHROM and map only valid chromosomes
    df["CHROM"] = df["CHROM"].astype(str).str.strip().str.replace("^chr", "", regex=True)
    valid_chr_map = {str(i): i for i in range(1, 23)}
    valid_chr_map.update({"X": 23, "Y": 24})
    df = df[df["CHROM"].isin(valid_chr_map.keys())]
    df["CHR"] = df["CHROM"].map(valid_chr_map).astype(int)

    df["BP"] = df["POS"]
    df["P"] = 10 ** (-df["QUAL"].astype(float) / 200)
    df = df.dropna(subset=["CHR", "BP", "P", "SVTYPE"])

    offset = 0
    tick_positions = []
    tick_labels = []
    fig = go.Figure()

    for chrom in sorted(df["CHR"].unique()):
        chr_df = df[df["CHR"] == chrom].copy()
        chr_df["cumulative_bp"] = chr_df["BP"] + offset

        midpoint = chr_df["cumulative_bp"].median()
        tick_positions.append(midpoint)
        tick_labels.append(f"chr{chrom if chrom <= 22 else ('X' if chrom == 23 else 'Y')}")

        for svtype in chr_df["SVTYPE"].unique():
            sv_df = chr_df[chr_df["SVTYPE"] == svtype]
            fig.add_trace(go.Scattergl(
                x=sv_df["cumulative_bp"],
                y=-np.log10(sv_df["P"]),
                mode='markers',
                name=f"chr{chrom}-{svtype}",
                marker=dict(color=svtype_color_map.get(svtype, "#999999"),size=4),
                customdata=sv_df[["CHROM", "POS", "SVTYPE", "SVLEN", "QUAL"]],
                hovertemplate=(
                    "CHROM: %{customdata[0]}<br>"
                    "POS: %{customdata[1]}<br>"
                    "SVTYPE: %{customdata[2]}<br>"
                    "SVLEN: %{customdata[3]}<br>"
                    "QUAL: %{customdata[4]}<br>"
                    "-log10(P): %{y}<extra></extra>"
                )
            ))

        highlight_chr_df = chr_df[-np.log10(chr_df["P"]) > threshold]
        if not highlight_chr_df.empty:
            highlight_df_list.append(highlight_chr_df)

        offset = chr_df["cumulative_bp"].max() + 1_000_000

    if highlight_df_list:
        all_high = pd.concat(highlight_df_list)
        fig.add_trace(go.Scattergl(
            x=all_high["cumulative_bp"],
            y=-np.log10(all_high["P"]),
            mode='markers',
            name="Points of Interest",
            marker=dict(color="red", size=6),
            showlegend=True
        ))

    fig.add_shape(
        type="line",
        x0=0, x1=offset,
        y0=threshold, y1=threshold,
        line=dict(color="red", dash="dash")
    )

    fig.update_layout(
            title="Interactive Manhattan Plot for Structural Variants",
            xaxis=dict(tickmode='array', tickvals=tick_positions, ticktext=tick_labels),
            yaxis_title="-log10(P) (simulated from QUAL)",
            height=600,
            margin=dict(t=50, l=60, r=30, b=60)
        )

    # === SVTYPE Legend ===
    svtypes = df["SVTYPE"].unique().tolist()
    legend = html.Div([
        html.Div([
            html.Span(style={
                'backgroundColor': svtype_color_map.get(sv, '#999999'),
                'display': 'inline-block',
                'width': '15px',
                'height': '15px',
                'marginRight': '5px',
                'borderRadius': '3px'
            }),
            html.Span(f"{sv}", style={'marginRight': '15px', 'fontSize': '14px'})
        ], style={'display': 'inline-block'}) for sv in svtypes
    ], style={'padding': '10px 0', 'textAlign': 'left'})


    
    return html.Div([
        html.Div(
            "⚠️ If plot is empty, check QUAL column.",
            style={
                'color': 'orange',
                'fontWeight': 'bold',
                'textAlign': 'left',
                'marginTop': '0px',
                'marginBottom': '20px'
            }
        ),
        legend,
        dcc.Graph(figure=fig)

    ])

def plot_manhattan_svlen(df, selected_svtypes, threshold):
    import numpy as np
    import plotly.graph_objects as go
    from dash import dcc

    highlight_df_list = []

    df = df[df["SVTYPE"].isin(selected_svtypes)]
    df["CHROM"] = df["CHROM"].astype(str).str.strip().str.replace("^chr", "", regex=True)
    valid_chr_map = {str(i): i for i in range(1, 23)}
    valid_chr_map.update({"X": 23, "Y": 24})
    df = df[df["CHROM"].isin(valid_chr_map.keys())]
    df["CHR"] = df["CHROM"].map(valid_chr_map).astype(int)

    df["BP"] = pd.to_numeric(df["POS"], errors="coerce")
    df["SVLEN"] = pd.to_numeric(df["SVLEN"], errors="coerce")
    df = df.dropna(subset=["CHR", "BP", "SVLEN", "SVTYPE"])
    df = df[df["SVLEN"] > 0]  # remove zeros since log scale can't handle them

    offset = 0
    tick_positions = []
    tick_labels = []
    fig = go.Figure()

    for chrom in sorted(df["CHR"].unique()):
        chr_df = df[df["CHR"] == chrom].copy()
        chr_df["cumulative_bp"] = chr_df["BP"] + offset
        midpoint = chr_df["cumulative_bp"].median()
        tick_positions.append(midpoint)
        tick_labels.append(f"chr{chrom if chrom <= 22 else ('X' if chrom == 23 else 'Y')}")

        for svtype in chr_df["SVTYPE"].unique():
            sv_df = chr_df[chr_df["SVTYPE"] == svtype]
            fig.add_trace(go.Scattergl(
                x=sv_df["cumulative_bp"],
                y=sv_df["SVLEN"],
                mode='markers',
                name=f"chr{chrom}-{svtype}",
                marker=dict(size=4),
                customdata=sv_df[["CHROM", "POS", "SVTYPE", "SVLEN"]],
                hovertemplate=(
                    "CHROM: %{customdata[0]}<br>"
                    "POS: %{customdata[1]}<br>"
                    "SVTYPE: %{customdata[2]}<br>"
                    "SVLEN: %{customdata[3]}<br>"
                    "SVLEN: %{y}<extra></extra>"
                )
            ))

        highlight_chr_df = chr_df[chr_df["SVLEN"] > threshold]
        if not highlight_chr_df.empty:
            highlight_df_list.append(highlight_chr_df)

        offset = chr_df["cumulative_bp"].max() + 1_000_000

    if highlight_df_list:
        all_high = pd.concat(highlight_df_list)
        fig.add_trace(go.Scattergl(
            x=all_high["cumulative_bp"],
            y=all_high["SVLEN"],
            mode='markers',
            name="Points of Interest",
            marker=dict(color="red", size=6),
            showlegend=True
        ))

    fig.add_shape(
        type="line",
        x0=0, x1=offset,
        y0=threshold, y1=threshold,
        line=dict(color="red", dash="dash")
    )

    fig.update_layout(
        title="Interactive Manhattan Plot for SVLEN (log scale)",
        xaxis=dict(tickmode='array', tickvals=tick_positions, ticktext=tick_labels),
        yaxis=dict(title="SVLEN", type="log"),
        height=600,
        margin=dict(t=50, l=60, r=30, b=60)
    )

    return dcc.Graph(figure=fig)

def plot_circos(graph_type, selected_svtypes, circos_data, svtype_color_map):
    if not circos_data:
        return html.Div("No Circos data available.")

    tracks = update_tracks(graph_type, selected_svtypes, circos_data, svtype_color_map)
 # 🔶 Create a color legend based on selected SVTYPEs
    legend = html.Div([
        html.Div([
            html.Span(style={
                'backgroundColor': svtype_color_map.get(sv, '#999'),
                'display': 'inline-block',
                'width': '15px',
                'height': '15px',
                'marginRight': '5px',
                'borderRadius': '3px'
            }),
            html.Span(f"{sv}", style={'marginRight': '15px', 'fontSize': '14px'})
        ], style={'display': 'inline-block'}) for sv in selected_svtypes
    ], style={'padding': '10px 0', 'textAlign': 'left'})

    
    return html.Div([
        legend,
        dashbio.Circos(
        id='circos-plot',
        layout=circos_data["layout"],
        selectEvent={"0": "hover"},
        tracks=tracks,
        config={
            "innerRadius": 250,
            "outerRadius": 300,
            "ticks": {
                "display": True,
                "spacing": 10000000,
                "labelSpacing": 5
            },
            "labelLayout": {
                "spacing": 60,
                "radialOffset": 90
            },
            "tooltipContent": {
                "source": "source.id + ':' + source.start + '-' + source.end",
                "target": "target.id + ':' + target.start + '-' + target.end"
            }
        },
        size=800
    )
])

def load_vcf_dataframe(file_path, source_type):

    try:
        if not file_path or not os.path.exists(file_path):
            return None

        if source_type == 'evalsvcallers':
            df = pd.read_csv(file_path, sep="\t", header=None)
            # ✅ Step 1: Split TF + CHROM from column 0 if needed
            
            tf_chrom = df[0].astype(str).str.strip().str.split(" ", expand=True)
            df.insert(0, "TF", tf_chrom[0])
            df.insert(1, "CHROM", tf_chrom[1])
            df = df.drop(columns=[0])
            df = df[df["TF"] == "TP"]
            df = df.drop(columns=["TF"])
            # ✅ Rename columns if not already done
            if df.shape[1] == 8 and not all(col in df.columns for col in ["CHROM", "POS", "SVTYPE", "QUAL", "INFO"]):
                df.columns = ["CHROM", "POS", "SVTYPE", "col4", "col5", "QUAL", "FILTER", "INFO"]
        else:
            try:
                df = pd.read_csv(file_path, sep="\t", comment="#", header=None)
                if df.shape[1] >= 8:
                    df = df.iloc[:, :8]
                    df.columns = ["CHROM", "POS", "SVTYPE", "col4", "col5", "QUAL", "FILTER", "INFO"]
                else:
                    raise ValueError("VCF/Caller file must have at least 8 columns.")
            except Exception as e:
                print(f"❌ Failed to load VCF/Caller file: {e}")
                return None

        # Clean CHROM
        df["CHROM"] = df["CHROM"].astype(str).str.strip().str.replace("^chr", "", regex=True)

        # Extract SVTYPE from INFO if missing or empty
        if df["SVTYPE"].isnull().all() or df["SVTYPE"].astype(str).str.strip().eq("").all():
            df["SVTYPE"] = df["INFO"].str.extract(r"SVTYPE=([^;]+)")

        # Ensure QUAL is numeric
        df["QUAL"] = pd.to_numeric(df["QUAL"], errors='coerce')

        # Extract SVLEN from INFO
        df["SVLEN"] = df["INFO"].apply(
            lambda x: int(re.search(r"SVLEN=-?\d+", x).group().split("=")[1])
            if pd.notnull(x) and isinstance(x, str) and re.search(r"SVLEN=-?\d+", x)
            else None
        )

        return df

    except Exception as e:
        print(f"❌ Unexpected error in load_vcf_dataframe(): {e}")
        return None
def update_tracks(graph_type, selected_svtypes, circos_data, svtype_color_map):
    if graph_type == "histogram":
        tracks = []
        for track in circos_data.get("tracks", []):
            if track["name"] in selected_svtypes:
                tracks.append({
                    "type": "HISTOGRAM",
                    "data": track["data"],
                    "config": {
                        "color": svtype_color_map.get(track["name"], "#95a5a6"),
                        "innerRadius": 200,
                        "outerRadius": 250,
                        "opacity": 0.9,
                        "tooltipContent": {"name": "block_id", "value": "value"}
                    }
                })
        return tracks

    return []

# ---- File Download ----
def serve_file_for_download(filename):
    return flask.send_from_directory(VISUALIZATION_OUTPUT_DIR, filename, as_attachment=True)
