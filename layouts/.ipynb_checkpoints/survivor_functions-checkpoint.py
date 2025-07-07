import os
import subprocess
import base64
import datetime
import pandas as pd
from dash import html, dash_table

# Define paths
UPLOAD_DIRECTORY = "./uploaded_files/"
SURVIVOR_OUTPUT_DIR = os.path.join(UPLOAD_DIRECTORY, "survivor")
os.makedirs(SURVIVOR_OUTPUT_DIR, exist_ok=True)

def save_file(name, content):
    """Decode and save uploaded file to survivor folder."""
    data = content.encode("utf8").split(b";base64,")[1]
    with open(os.path.join(SURVIVOR_OUTPUT_DIR, name), "wb") as fp:
        fp.write(base64.b64decode(data))

def parse_uploaded_files(contents_list, names_list, dates_list):
    """Save and return file info as list items."""
    items = []
    for c, n, d in zip(contents_list, names_list, dates_list):
        save_file(n, c)
        timestamp = datetime.datetime.fromtimestamp(d)
        items.append(html.Li(f"{n} (Uploaded at {timestamp})"))
    return html.Ul(items)

def run_survivor_merge(params):
    """Run SURVIVOR merge command with provided parameters."""
    sample_files_path = os.path.join(SURVIVOR_OUTPUT_DIR, "sample_files")
    merged_file_path = os.path.join(SURVIVOR_OUTPUT_DIR, "merged_output.vcf")
    univar_path = os.path.join(UPLOAD_DIRECTORY, "univar.vcf")

    vcf_files = [
        os.path.join(SURVIVOR_OUTPUT_DIR, f)
        for f in os.listdir(SURVIVOR_OUTPUT_DIR)
        if f.endswith(".vcf") and f != os.path.basename(merged_file_path)
    ]

    if params['param_univar'] == 1 and os.path.exists(univar_path):
        vcf_files.append(univar_path)

    if len(vcf_files) < 2:
        return "Error: At least two VCF files are required to merge."

    with open(sample_files_path, "w") as sf:
        for vcf in vcf_files:
            sf.write(vcf + "\n")

    survivor_executable = "./SURVIVOR/Debug/SURVIVOR"
    if not os.access(survivor_executable, os.X_OK):
        os.chmod(survivor_executable, 0o755)

    command = [
        survivor_executable, "merge", sample_files_path,
        str(params['max_distance']),
        str(params['min_callers']),
        str(params['type_match']),
        str(params['strand_match']),
        str(params['allow_duplicates']),
        str(params['min_sv_size']),
        merged_file_path
    ]

    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
        return f"✅ SURVIVOR merge completed! Output saved to {merged_file_path}"
    except Exception as e:
        return f"❌ Error during SURVIVOR merge: {str(e)}"

def get_merge_preview():
    """Return a DataTable preview of merged VCF."""
    merged_file_path = os.path.join(SURVIVOR_OUTPUT_DIR, "merged_output.vcf")
    try:
        with open(merged_file_path, "r") as file:
            for line in file:
                if line.startswith("#CHROM"):
                    columns = line.strip().split("\t")
                    break
            else:
                return html.Div("❌ Error: #CHROM header not found.")

        data_rows = []
        with open(merged_file_path, "r") as file:
            for line in file:
                if not line.startswith("#"):
                    data_rows.append(line.strip().split("\t"))
                if len(data_rows) == 5:
                    break

        df = pd.DataFrame(data_rows, columns=columns)
        return dash_table.DataTable(
            data=df.to_dict('records'),
            columns=[{"name": col, "id": col} for col in df.columns],
            style_table={'overflowX': 'auto'},
            style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'},
            style_header={'fontWeight': 'bold', 'backgroundColor': 'lightgrey'},
            style_data={'border': '1px solid black'},
        )
    except Exception as e:
        return html.Div(f"❌ Error reading merged file: {str(e)}")
