import os
import base64
import datetime
import subprocess
import json
import sys
from dash import html
import shutil
import plotly.express as px
from layouts.visualize_functions import load_vcf_dataframe
APP_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(APP_ROOT, "..")))
UPLOAD_DIRECTORY = os.path.join(APP_ROOT, "..", "uploaded_files")
TRUVARI_OUTPUT_DIR = os.path.join(UPLOAD_DIRECTORY, "truvari_output")
os.makedirs(TRUVARI_OUTPUT_DIR, exist_ok=True)
# Son oluşturulan run klasörü (global pointer)


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

def save_uploaded_file(name, content):
    data = content.encode("utf8").split(b";base64,")[1]
    filepath = os.path.join(TRUVARI_OUTPUT_DIR, name)
    with open(filepath, "wb") as f:
        f.write(base64.b64decode(data))
    return os.path.abspath(filepath)

def _run(cmd, logf, env=None):
    logf.write(" ".join(cmd) + "\n")
    logf.flush()
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)
    logf.write(r.stdout + "\n")
    logf.flush()
    r.check_returncode()
    return r

def run_truvari_pipeline(base_choice, base_vcf_uploaded, comp_vcf, reference_fa, params=None):
    """
    Uçtan uca: sort/index -> ref.bed -> comp temizle (INFO/END>POS) -> onref -> GT header -> reheader -> sort -> chr filtresi -> truvari bench
    base_vcf: .vcf veya .vcf.gz
    comp_vcf: .vcf veya .vcf.gz
    reference_fa: .fa
    """
    if base_choice == 'univar':
        base_vcf = os.path.join(UPLOAD_DIRECTORY, "univar.vcf")
        if not os.path.exists(base_vcf):
            return html.Div("❌ Univar SV catalog (univar.vcf) not found."), None
    else:
        if not base_vcf_uploaded:
            return html.Div("❌ No custom base VCF uploaded."), None
        base_vcf = os.path.abspath(base_vcf_uploaded)



    
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = os.path.join(TRUVARI_OUTPUT_DIR, f"run_{ts}")
    os.makedirs(run_dir, exist_ok=True)
    output_dir = os.path.join(run_dir, "truvari_output")
    # Truvari kendi oluşturmak ister; varsa sil, asla önceden yaratma
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    log_path = os.path.join(run_dir, "truvari_run.log")
    with open(log_path, "w") as logf:
        try:
            # 1) BASE sort/index
            logf.write("🔧 Sorting and indexing BASE VCF\n")
            base_sorted = os.path.join(run_dir, "base.sorted.vcf.gz")
            _run(["bcftools", "sort", "-Oz", "-o", base_sorted, base_vcf], logf)
            _run(["tabix", "-f", "-p", "vcf", base_sorted], logf)

            # 2) COMP sort/index
            logf.write("🔧 Sorting and indexing COMP VCF\n")
            comp_sorted = os.path.join(run_dir, "comp.sorted.vcf.gz")
            _run(["bcftools", "sort", comp_vcf, "-Oz", "-o", comp_sorted], logf)
            _run(["tabix", "-f", "-p", "vcf", comp_sorted], logf)

            # 3) REF .fai + ref.bed(.gz)
            logf.write("🔧 Creating .fai from reference\n")
            _run(["samtools", "faidx", reference_fa], logf)
            fai_path = reference_fa + ".fai"

            logf.write("🔧 Creating ref.bed from .fai\n")
            bed_path = os.path.join(run_dir, "ref.bed")
            with open(fai_path) as fi, open(bed_path, "w") as fo:
                for line in fi:
                    chrom, length = line.strip().split("\t")[:2]
                    fo.write(f"{chrom}\t0\t{length}\n")
            _run(["bgzip", "-f", bed_path], logf)
            bed_gz = bed_path + ".gz"
            _run(["tabix", "-f", "-p", "bed", bed_gz], logf)

            # 4) COMP temizliği (INFO/END>POS)
            logf.write("🔧 Cleaning COMP: filter INFO/END>POS\n")
            comp_clean = os.path.join(run_dir, "comp.cleaned.vcf.gz")
            _run(["bcftools", "view", "-i", "INFO/END>POS", comp_sorted, "-Oz", "-o", comp_clean], logf)
            _run(["tabix", "-f", "-p", "vcf", comp_clean], logf)

            # 5) COMP -> onref
            logf.write("🔧 Filtering COMP to regions in ref.bed.gz\n")
            onref_vcf = os.path.join(run_dir, "comp.onref.vcf.gz")
            _run(["bcftools", "view", "-R", bed_gz, "-Oz", "-o", onref_vcf, comp_clean], logf)

            # 6) GT header ekle
            logf.write("🔧 Adding GT header to COMP VCF\n")
            hdr_file = os.path.join(run_dir, "add.hdr")
            with open(hdr_file, "w") as f:
                f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            hdrfix_vcf = os.path.join(run_dir, "comp.onref.hdrfix.vcf.gz")
            _run(["bcftools", "annotate", "-h", hdr_file, "-Oz", "-o", hdrfix_vcf, onref_vcf], logf)
            _run(["tabix", "-f", "-p", "vcf", hdrfix_vcf], logf)

            # 7) Reheader
            logf.write("🔧 Reheader with .fai chrom names\n")
            rehead_vcf = os.path.join(run_dir, "comp.onref.rehead.vcf.gz")
            _run(["bcftools", "reheader", "-f", fai_path, hdrfix_vcf, "-o", rehead_vcf], logf)

            # 8) Son sort/index
            logf.write("🔧 Sorting and indexing reheaded COMP\n")
            final_sorted = os.path.join(run_dir, "comp.onref.sorted.vcf.gz")
            _run(["bcftools", "sort", "-T", os.path.join(run_dir, "tmp"), "-Oz", "-o", final_sorted, rehead_vcf], logf)
            _run(["tabix", "-f", "-p", "vcf", final_sorted], logf)

            # 9) chr1..22,X,Y
            logf.write("🔧 Filtering COMP to chr1..22,X,Y\n")
            final_main = os.path.join(run_dir, "comp.main.vcf.gz")
            chr_list = ",".join([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])
            _run(["bcftools", "view", "-r", chr_list, "-Oz", "-o", final_main, final_sorted], logf)
            _run(["tabix", "-f", "-p", "vcf", final_main], logf)

            
            # 🔧 Reheader BASE with .fai  (contig length lines için)
            logf.write("🔧 Reheader BASE with .fai\n")
            base_rehead = os.path.join(run_dir, "base.rehead.vcf.gz")
            _run(["bcftools", "reheader", "-f", fai_path, base_sorted, "-o", base_rehead], logf)
            _run(["tabix", "-f", "-p", "vcf", base_rehead], logf)
            base_sorted = base_rehead  # Bundan sonra BASE için bunu kullan
            
            # 🔧 Reheader COMP with .fai  (contig length lines için)
            logf.write("🔧 Reheader COMP with .fai\n")
            comp_rehead = os.path.join(run_dir, "comp.main.rehead.vcf.gz")
            _run(["bcftools", "reheader", "-f", fai_path, final_main, "-o", comp_rehead], logf)
            _run(["tabix", "-f", "-p", "vcf", comp_rehead], logf)
            final_main = comp_rehead  # Bundan sonra COMP için bunu kullan


            # 10) truvari bench
            logf.write("🚀 Running truvari bench\n")
            reference_fa_abs = os.path.abspath(reference_fa)
            base_vcf_final   = os.path.abspath(base_sorted)
            comp_vcf_final   = os.path.abspath(final_main)
            output_dir_abs   = os.path.abspath(output_dir)

            # Önce 'truvari' binary denenir, yoksa -m truvari.bench
            try_cmd = ["truvari", "bench",
                       "-b", base_vcf_final,
                       "-c", comp_vcf_final,
                       "-f", reference_fa_abs,
                       "-o", output_dir_abs]

            if params:
                for flag in ["--typeignore", "--use-lev", "--gtcomp", "--passonly", "--multimatch", "--giabreport", "--debug", "--prog"]:
                    if params.get(flag, False):
                        try_cmd.append(flag)
                for key in ["--refdist","--pctsim","--pctsize","--pctovl","--sizemin","--sizefilt","--sizemax","--no-ref","--bSample","--cSample"]:
                    if key in params and params[key] not in [None, ""]:
                        try_cmd += [key, str(params[key])]

            logf.write("Running:\n" + " ".join(try_cmd) + "\n")
            try:
                _run(try_cmd, logf)
            except subprocess.CalledProcessError:
                # fallback: python -m truvari.bench
                env = os.environ.copy()
                app_parent = os.path.abspath(os.path.join(APP_ROOT, ".."))
                env["PYTHONPATH"] = (env.get("PYTHONPATH", "") + os.pathsep + app_parent).strip(os.pathsep)
                py_cmd = [sys.executable, "-m", "truvari.bench",
                          "-b", base_vcf_final,
                          "-c", comp_vcf_final,
                          "-f", reference_fa_abs,
                          "-o", output_dir_abs]
                if params:
                    for flag in ["--typeignore","--use-lev","--gtcomp","--passonly","--multimatch","--giabreport","--debug","--prog"]:
                        if params.get(flag, False):
                            py_cmd.append(flag)
                    for key in ["--refdist","--pctsim","--pctsize","--pctovl","--sizemin","--sizefilt","--sizemax","--no-ref","--bSample","--cSample"]:
                        if key in params and params[key] not in [None, ""]:
                            py_cmd += [key, str(params[key])]
                logf.write("Fallback:\n" + " ".join(py_cmd) + "\n")
                _run(py_cmd, logf, env=env)

            # summary.json
            summary_json = os.path.join(output_dir_abs, "summary.json")
            summary_data = {}
            if os.path.exists(summary_json):
                with open(summary_json) as f:
                    summary_data = json.load(f)

            return html.Div([
                html.Div("✅ Truvari bench completed."),
                html.Div(f"📁 Output directory: {output_dir_abs}"),
                html.Br()
            ]), {
                "tp_comp": os.path.join(output_dir_abs, "tp-comp.vcf.gz"),
                "tp_base": os.path.join(output_dir_abs, "tp-base.vcf.gz"),
                "fp": os.path.join(output_dir_abs, "fp.vcf.gz"),
                "fn": os.path.join(output_dir_abs, "fn.vcf.gz"),
                "summary": os.path.join(output_dir_abs, "summary.json")
            }

        except Exception as e:
            return html.Div(f"❌ Truvari error: {str(e)}"), None
import gzip
from collections import defaultdict

def _open_any(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")

def _parse_vcf_record(parts, source_kind):
    """
    parts: VCF split line
    source_kind: 'TP-BASE' veya 'TP-COMP'
    """
    chrom = parts[0].strip()
    try:
        pos = int(parts[1])
    except:
        return None

    # INFO alanını ayıkla
    info = {}
    if len(parts) >= 8:
        for kv in parts[7].split(";"):
            if "=" in kv:
                k, v = kv.split("=", 1)
                info[k] = v

    svtype = info.get("SVTYPE", "UNK").strip().upper()
    if svtype == "BND":
        return None

    # QUAL
    try:
        qual = int(float(parts[5])) if parts[5].replace('.', '', 1).isdigit() else 0
    except:
        qual = 0

    # END ve SVLEN
    end = None
    svlen = None
    if "END" in info:
        try:
            end = int(info["END"])
        except:
            end = None
    if "SVLEN" in info:
        try:
            svlen = abs(int(info["SVLEN"]))
        except:
            svlen = None

    # Fallback hesaplar
    if end is None and svlen is not None:
        end = pos + abs(svlen)
    if svlen is None and end is not None:
        svlen = abs(end - pos)
    if end is None and svlen is None:
        end = pos + 1000
        svlen = 1000

    chrom = f"chr{chrom}" if not chrom.startswith("chr") else chrom

    return {
        "chrom": chrom,
        "pos": pos,
        "end": end,
        "svtype": svtype,
        "svlen": svlen,
        "qual": qual,
        "source": source_kind  # TP-BASE / TP-COMP
    }

def _hex_to_rgba(color_str, alpha=1.0):
    """
    Accepts both hex ("#RRGGBB") and rgb/rgba("rgb(r,g,b)") formats
    Returns rgba string with given alpha
    """
    color_str = color_str.strip()

    if color_str.startswith("#"):
        # HEX → RGBA
        hex_color = color_str.lstrip("#")
        if len(hex_color) == 3:
            hex_color = "".join([c*2 for c in hex_color])
        r = int(hex_color[0:2], 16)
        g = int(hex_color[2:4], 16)
        b = int(hex_color[4:6], 16)
        return f"rgba({r},{g},{b},{alpha})"

    elif color_str.startswith("rgb"):
        # Already rgb/rgba → sadece alpha güncelle
        nums = color_str[color_str.find("(")+1:color_str.find(")")].split(",")
        nums = [n.strip() for n in nums]
        if len(nums) >= 3:
            r, g, b = nums[:3]
            return f"rgba({r},{g},{b},{alpha})"
    
    # fallback
    return f"rgba(128,128,128,{alpha})"
def vcf_to_circos_json_truvari(tp_base_vcf, tp_comp_vcf, output_json_path,
                               svtype_color_map, chrom_color_map,
                               hg38_lengths, hg19_lengths, bin_size=10_000_000):
    """
    İki Truvari VCF’ini (TP-BASE, TP-COMP) okuyup tek bir Circos JSON üretir.
    - Histogramları BASE ve COMP için ayrı ayrı yazar (aynı SVTYPE için iki track).
    - Chord’ları BASE ve COMP için farklı opaklıkta renklendirir (aynı SVTYPE rengi, BASE=1.0, COMP=0.45).
    """
    # Hangi genom? Basit tespit: ilk vcf’in ilk veri satırı
    def detect_genome_from_vcf(vcf):
        with _open_any(vcf) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                chrom = line.split("\t")[0]
                return "hg38" if chrom.startswith("chr") else "hg19"
        return "hg38"

    genome = detect_genome_from_vcf(tp_base_vcf)
    chromosome_lengths = hg38_lengths if genome == "hg38" else hg19_lengths
    chromosome_lengths = {chrom: max(1, length) for chrom, length in chromosome_lengths.items()}
    valid_chroms = set(chromosome_lengths.keys())

    # Toplayıcılar
    layout_chromosomes = set()
    raw_chords = []
    chord_metadata = []

    # BASE/COMP & SVTYPE’a göre histogram sayacı
    # histograms[source][svtype][(chrom, bin_start)] = count
    histograms = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    svtypes_found = set()
    sources_found = set()

    def ingest_vcf(vcf_path, source_kind):
        with _open_any(vcf_path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    continue
                rec = _parse_vcf_record(parts, source_kind)
                if rec is None:
                    continue
                chrom = rec["chrom"]
                if chrom not in valid_chroms:
                    continue

                layout_chromosomes.add(chrom)
                svtypes_found.add(rec["svtype"])
                sources_found.add(source_kind)

                # Histogram bin sayımı
                bin_start = (rec["pos"] // bin_size) * bin_size
                histograms[source_kind][rec["svtype"]][(chrom, bin_start)] += 1

                # Chord (aynı krom içinde pos→end tek çizgi)
                source = {"id": chrom, "start": rec["pos"], "end": rec["end"]}
                target = {"id": chrom, "start": rec["end"], "end": rec["end"] + 1}
                base_color = svtype_color_map.get(rec["svtype"], "#95a5a6")
                color = _hex_to_rgba(base_color, 1.0 if source_kind == "TP-BASE" else 0.45)

                raw_chords.append({
                    "source": source,
                    "target": target,
                    "color": color,
                    "value": rec["svlen"]
                })
                chord_metadata.append({
                    "source": source,
                    "target": target,
                    "svtype": rec["svtype"],
                    "source_kind": source_kind,
                    "qual": rec["qual"],
                    "svlen": rec["svlen"]
                })

    # İki VCF’i içeri al
    ingest_vcf(tp_base_vcf, "TP-BASE")
    ingest_vcf(tp_comp_vcf, "TP-COMP")

    # Layout (görünen kromlar)
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

    # Tracks: BASE ve COMP için ayrı histogram track’ları
    tracks = []
    # Sıra sabit olsun diye SVTYPE’ları sıralayalım
    for svtype in sorted(svtypes_found):
        for src in ["TP-BASE", "TP-COMP"]:
            bins = histograms[src][svtype]
            if not bins:
                continue
            data = [
                {
                    "block_id": chrom,
                    "start": start,
                    "end": start + bin_size,
                    "value": count,
                    "name": f"{src}:{svtype}"  # ← buraya 'name' ekleniyor!
                }
                for (chrom, start), count in bins.items()
                if count > 0
            ]
            if data:
                base_hex = svtype_color_map.get(svtype, "#95a5a6")
                color = _hex_to_rgba(base_hex, 1.0 if src == "TP-BASE" else 0.55)
                tracks.append({
                    "type": "histogram",
                    "data": data,
                    "name": f"{src}:{svtype}",
                    "color": color,
                    "height": 0.08
                })

    circos_json = {
        "layout": layout,
        "chords": raw_chords,
        "chord_metadata": chord_metadata,
        "tracks": tracks,
        "detected_svtypes": sorted(svtypes_found),
        "detected_sources": sorted(sources_found),  # ["TP-BASE", "TP-COMP"]
        "bin_size": bin_size,
        "genome": genome
    }

    with open(output_json_path, "w") as f:
        json.dump(circos_json, f, indent=2)

def update_tracks(graph_type, selected_svtypes, circos_data, svtype_color_map):
    """
    graph_type: 'hist' veya 'chord' (sen nasıl kullanıyorsan ona uy)
    selected_svtypes: ["DEL","INS",...] gibi seçim
    circos_data: vcf_to_circos_json_truvari çıktısı
    """
    selected = set(selected_svtypes) if selected_svtypes else set(circos_data.get("detected_svtypes", []))

    # Histogram track'larını SVTYPE filtresine göre süz
    if graph_type == "hist":
        tracks = []
        for tr in circos_data.get("tracks", []):
            # isim formatı "TP-BASE:DEL" gibi
            name = tr.get("name", "")
            try:
                src, sv = name.split(":", 1)
            except ValueError:
                # Eski format desteklemek istersen:
                sv = name
            if sv in selected:
                tracks.append(tr)
        return tracks

    # Chord'lar için "highlight" track
    if graph_type == "chord":
        chords = []
        for ch, meta in zip(circos_data.get("chords", []), circos_data.get("chord_metadata", [])):
            if meta.get("svtype") in selected:
                # chord tek başına döndürülmüyor; Circos'ta chord'lar ayrı argüman
                chords.append(ch)
        # dashbio.Circos chord’ları tracks listesi yerine component argümanı ile kullanılıyor.
        # Ama sen plot_circos'ta chords'i doğrudan vermiyorsun; istersen burada chords'u
        # "chords" isimli track gibi döndürüp plot_circos’ta yakalayabilirsin.
        # Mevcut kodunda tracks dışı chords’u kullanmıyorsan hist’e devam edebilirsin.
        # Örnek olarak chord yoksa fallback:
        return circos_data.get("tracks", [])

    # Varsayılan: tüm histogramlar
    return circos_data.get("tracks", [])


def plot_manhattan_truvari(pairs, selected_svtypes=None, threshold=6):
    """
    Manhattan plot for Truvari (TP-BASE vs TP-COMP) with legend grouped by (Chromosome × Source).
    Legend entry text also lists which SVTYPEs exist in that group, e.g. "chr1TPBASE (DEL, INS)".
    Clicking a legend entry toggles ALL traces for that chromosome+source (groupclick=togglegroup).
    Points still colored uniquely per (SOURCE, SVTYPE) combo.
    """
    import os
    import numpy as np
    import pandas as pd
    import plotly.graph_objects as go
    from dash import html, dcc
    import plotly.express as px

    # ---- helper to load and tag ----
    def _load_and_tag(vcf_path, tag):
        if not vcf_path or not os.path.exists(vcf_path):
            return None
        df = load_vcf_dataframe(vcf_path, "caller")
        if df is None or df.empty:
            return None
        df = df.copy()
        # CHROM normalize to 1..22,X,Y
        df["CHROM"] = df["CHROM"].astype(str).str.strip().str.replace("^chr", "", regex=True)
        # Robust SVTYPE
        df["SVTYPE"] = df["INFO"].astype(str).str.extract(r'(?:^|;)SVTYPE=([^;]+)', expand=False)
        # SVLEN
        if "SVLEN" not in df or df["SVLEN"].isna().all():
            df["SVLEN"] = df["INFO"].astype(str).str.extract(r'(?:^|;)SVLEN=(-?\d+)', expand=False)
        df["SVLEN"] = pd.to_numeric(df["SVLEN"], errors="coerce")
        # POS/QUAL numeric
        df["POS"]  = pd.to_numeric(df.get("POS"), errors="coerce")
        df["QUAL"] = pd.to_numeric(df.get("QUAL"), errors="coerce")
        df["__SOURCE__"] = tag
        # drop bad rows
        df = df.dropna(subset=["SVTYPE", "POS", "QUAL"])
        return df

    # ---- collect TP-BASE/TP-COMP ----
    dfs = []
    for p, lab in (pairs or []):
        if lab not in ("TP-BASE", "TP-COMP"):
            continue
        d = _load_and_tag(p, lab)
        if d is not None and not d.empty:
            dfs.append(d)
    if not dfs:
        return html.Div("⚠️ Manhattan: no TP-BASE/TP-COMP data.")
    df = pd.concat(dfs, ignore_index=True)

    # ---- optional SVTYPE filter ----
    if selected_svtypes:
        df = df[df["SVTYPE"].isin(selected_svtypes)]
    if df.empty:
        return html.Div("⚠️ Manhattan: no data after SVTYPE filter.")

    # ---- keep valid chromosomes only ----
    valid_chr_map = {str(i): i for i in range(1, 23)}
    valid_chr_map.update({"X": 23, "Y": 24})
    df = df[df["CHROM"].isin(valid_chr_map.keys())].copy()
    if df.empty:
        return html.Div("⚠️ Manhattan: no valid chromosomes (1..22, X, Y).")

    df["CHR"] = df["CHROM"].map(valid_chr_map).astype(int)

    # ---- QUAL surrogate: -log10(P) = QUAL/200 ----
    df["neglog10P"] = df["QUAL"] / 200.0
    df["BP"] = df["POS"]

    # ---- build cumulative genome axis ----
    offsets = {}
    tick_positions, tick_labels = [], []
    offset = 0.0
    for chrom in sorted(df["CHR"].unique()):
        chr_df = df[df["CHR"] == chrom]
        if chr_df.empty:
            continue
        offsets[chrom] = offset
        span_max = float(chr_df["BP"].max())
        midpoint = offset + (float(chr_df["BP"].median()) if not np.isnan(chr_df["BP"].median()) else span_max / 2.0)
        tick_positions.append(midpoint)
        tick_labels.append(f"chr{chrom if chrom <= 22 else ('X' if chrom == 23 else 'Y')}")
        offset += span_max + 1_000_000  # gap

    df["cumulative_bp"] = df.apply(lambda r: r["BP"] + offsets.get(r["CHR"], 0.0), axis=1)

    # ---- unique colors per (SOURCE, SVTYPE) combo (kept for point-level distinction) ----
    big_palette = (
        px.colors.qualitative.Vivid
        + px.colors.qualitative.Set2
        + px.colors.qualitative.Plotly
        + px.colors.qualitative.Safe
        + px.colors.qualitative.Bold
        + px.colors.qualitative.Pastel
    )
    # stable svtype order
    sv_order = list(pd.unique(df["SVTYPE"]))
    combo_keys = []
    for grp in ("TP-BASE", "TP-COMP"):
        present = df[df["__SOURCE__"] == grp]["SVTYPE"].dropna().unique().tolist()
        ordered = [sv for sv in sv_order if sv in present] + [sv for sv in present if sv not in sv_order]
        combo_keys += [f"{grp}:{sv}" for sv in ordered]

    def hsl_hash_hex(key: str):
        h = abs(hash(key))
        hue = h % 360
        sat = 65 + (h % 25)     # 65..89
        lig = 45 + (h % 10)     # 45..54
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

    combo_color = {}
    for i, key in enumerate(combo_keys):
        combo_color[key] = big_palette[i] if i < len(big_palette) else hsl_hash_hex(key)

    # ---- build figure; legend grouped by (Chromosome × Source) ----
    fig = go.Figure()

    for chrom in sorted(df["CHR"].unique()):
        chr_label = f"chr{chrom if chrom <= 22 else ('X' if chrom == 23 else 'Y')}"
        for grp in ("TP-BASE", "TP-COMP"):
            sub = df[(df["CHR"] == chrom) & (df["__SOURCE__"] == grp)]
            if sub.empty:
                continue

            # Legend group label like "chr1TPBASE"
            group_label = f"{chr_label}-{grp}"

            # SVTYPE list for legend text (ordered)
            present_types = sub["SVTYPE"].dropna().unique().tolist()
            present_types = [sv for sv in sv_order if sv in present_types] + \
                            [sv for sv in present_types if sv not in sv_order]
            types_text = ", ".join(present_types)
            legend_name = f"{group_label} ({types_text})" if types_text else group_label

            first_in_group = True
            # stable per-group svtype order
            for sv in present_types:
                sv_df = sub[sub["SVTYPE"] == sv]
                if sv_df.empty:
                    continue
                col = combo_color.get(f"{grp}:{sv}", "#999999")

                fig.add_trace(go.Scattergl(
                    x=sv_df["cumulative_bp"],
                    y=sv_df["neglog10P"],
                    mode="markers",
                    name=legend_name if first_in_group else f"{group_label}:{sv}",  # sadece ilk trace legend'de
                    showlegend=first_in_group,
                    legendgroup=group_label,
                    marker=dict(color=col, size=4),
                    customdata=sv_df[["CHROM", "POS", "SVTYPE", "SVLEN", "QUAL", "__SOURCE__"]],
                    hovertemplate=(
                        "LegendGroup: " + group_label + "<br>" +
                        "SRC: %{customdata[5]}<br>" +
                        "CHROM: %{customdata[0]}<br>" +
                        "POS: %{customdata[1]}<br>" +
                        "SVTYPE: %{customdata[2]}<br>" +
                        "SVLEN: %{customdata[3]}<br>" +
                        "QUAL: %{customdata[4]}<br>" +
                        "-log10(P): %{y}<extra></extra>"
                    )
                ))
                first_in_group = False

    # Highlight points above threshold
    high_df = df[df["neglog10P"] > float(threshold)]
    if not high_df.empty:
        fig.add_trace(go.Scattergl(
            x=high_df["cumulative_bp"],
            y=high_df["neglog10P"],
            mode="markers",
            name="Points of Interest",
            marker=dict(color="#d62728", size=6),
            showlegend=True
        ))

    # Threshold line
    fig.add_shape(
        type="line",
        x0=0, x1=float(df["cumulative_bp"].max()) + 1_000_000,
        y0=float(threshold), y1=float(threshold),
        line=dict(color="#d62728", dash="dash")
    )

    fig.update_layout(
        title="Truvari Manhattan — TP-BASE vs TP-COMP (legend = Chrom×Source + SVTYPE list)",
        xaxis=dict(tickmode='array', tickvals=tick_positions, ticktext=tick_labels),
        yaxis_title="-log10(P) (≈ QUAL/200)",
        height=600,
        margin=dict(t=50, l=60, r=30, b=60),
        legend=dict(groupclick="togglegroup")   # legend'den grup halinde aç/kapa
    )

    return html.Div([
        dcc.Graph(figure=fig)
    ])
