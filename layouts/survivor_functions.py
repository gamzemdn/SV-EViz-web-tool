import os
import subprocess
import base64
import datetime
import pandas as pd
from dash import html, dash_table
import time
# Define paths
UPLOAD_DIRECTORY = "./uploaded_files/"
SURVIVOR_OUTPUT_DIR = os.path.join(UPLOAD_DIRECTORY, "survivor_output")
SURVIVOR_UPLOAD_CACHE_DIR = os.path.join(SURVIVOR_OUTPUT_DIR, "upload_cache")

os.makedirs(SURVIVOR_OUTPUT_DIR, exist_ok=True)
os.makedirs(SURVIVOR_UPLOAD_CACHE_DIR, exist_ok=True)

def _decode_upload_content(content):
    content_type, content_string = content.split(",", 1)
    return base64.b64decode(content_string)


def _detect_survivor_extension(path):
    lower = str(path).lower()

    if lower.endswith(".vcf.gz"):
        return ".vcf.gz"
    if lower.endswith(".vcf"):
        return ".vcf"

    return os.path.splitext(path)[1]

def save_file(name, content):
    """
    Save uploaded SURVIVOR input files into temporary upload cache.
    Files are copied into run-specific inputs/ folder during merge.
    """
    os.makedirs(SURVIVOR_UPLOAD_CACHE_DIR, exist_ok=True)

    if not name or not content:
        raise ValueError("Uploaded file name or content is missing.")

    data = _decode_upload_content(content)

    base = os.path.basename(name)
    root, ext = os.path.splitext(base)

    save_path = os.path.abspath(os.path.join(SURVIVOR_UPLOAD_CACHE_DIR, base))

    if os.path.exists(save_path):
        stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        save_path = os.path.abspath(
            os.path.join(SURVIVOR_UPLOAD_CACHE_DIR, f"{root}_{stamp}{ext}")
        )

    with open(save_path, "wb") as fp:
        fp.write(data)

    return save_path

def _copy_input_to_survivor_run(src_path, inputs_dir, index):
    src_path = os.path.abspath(src_path)

    if not os.path.exists(src_path):
        raise FileNotFoundError(f"Input file not found: {src_path}")

    ext = _detect_survivor_extension(src_path)
    dst_path = os.path.abspath(os.path.join(inputs_dir, f"caller_{index}{ext}"))

    import shutil
    if os.path.abspath(src_path) != os.path.abspath(dst_path):
        shutil.copy2(src_path, dst_path)

    return dst_path


def _write_survivor_input_manifest(manifest_path, original_files, run_files, params):
    with open(manifest_path, "w") as f:
        f.write("SURVIVOR Run Input Manifest\n")
        f.write("===========================\n\n")

        f.write("[INPUT CALLER VCF FILES]\n")
        for idx, (orig, run_copy) in enumerate(zip(original_files, run_files), start=1):
            f.write(f"Caller {idx} original path: {orig}\n")
            f.write(f"Caller {idx} run copy:      {run_copy}\n\n")

        f.write("[SURVIVOR PARAMETERS]\n")
        for key, value in params.items():
            f.write(f"{key}: {value}\n")


def prepare_vcf_files_for_merge(caller_paths, use_univar, ref_filename=None):
    vcf_files = [os.path.abspath(p) for p in caller_paths if p]
    vcf_files = list(dict.fromkeys(vcf_files))[:3]
    return vcf_files

    #univar_path = os.path.abspath(os.path.join(UPLOAD_DIRECTORY, "univar.vcf"))
   # ref_path = os.path.abspath(os.path.join(SURVIVOR_OUTPUT_DIR, ref_filename)) if ref_filename else None

   # if use_univar == 1 and os.path.exists(univar_path):
   #     if univar_path not in vcf_files:
   #         vcf_files.append(univar_path)

  #  elif use_univar == 0 and ref_path and os.path.exists(ref_path):
        # Only add if it's not the exact same file
   #     if not any(os.path.samefile(ref_path, p) for p in vcf_files):
   #         vcf_files.append(ref_path)

   # return vcf_files

def parse_uploaded_files(contents, names, dates):
    saved_paths = []
    messages = []

    for content, name in zip(contents, names):
        save_path = save_file(name, content)

        saved_paths.append(save_path)
        messages.append(
            html.P([
                f"✅ Uploaded: {os.path.basename(save_path)}",
                html.Br(),
                html.Small(f"Temporary path: {save_path}", style={"color": "gray"})
            ])
        )

    return html.Div(messages), saved_paths

def run_survivor_merge(params, vcf_files):
    """
    Run SURVIVOR merge in a run-specific folder.
    Inputs are copied to:
        survivor_output/run_YYYYMMDD_HHMMSS/inputs/
    Output is written to:
        survivor_output/run_YYYYMMDD_HHMMSS/survivor_output/
    """
    os.makedirs(SURVIVOR_OUTPUT_DIR, exist_ok=True)

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    run_dir = os.path.abspath(os.path.join(SURVIVOR_OUTPUT_DIR, f"run_{timestamp}"))
    inputs_dir = os.path.join(run_dir, "inputs")
    output_dir = os.path.join(run_dir, "survivor_output")

    os.makedirs(inputs_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    log_file = os.path.join(run_dir, "survivor_run.log")

    # Normalize paths and keep unique files
    original_files = [os.path.abspath(p) for p in (vcf_files or []) if p]
    original_files = list(dict.fromkeys(original_files))

    if len(original_files) < 2:
        return html.Div([
            "❌ Error: At least two caller VCF files are required for SURVIVOR merge.",
            html.Br(),
            html.Div("Current files:"),
            html.Ul([html.Li(os.path.basename(p)) for p in original_files]) if original_files else html.Div("— none —")
        ]), None

    warning_msg = None

    if len(original_files) > 3:
        warning_msg = html.Div(
            f"⚠️ Only the first 3 of {len(original_files)} uploaded files will be merged.",
            style={"color": "orange"}
        )

    selected_original_files = original_files[:3]

    try:
        run_input_files = []

        for idx, src in enumerate(selected_original_files, start=1):
            copied_path = _copy_input_to_survivor_run(src, inputs_dir, idx)
            run_input_files.append(copied_path)

        sample_files_path = os.path.join(inputs_dir, "sample_files.txt")

        with open(sample_files_path, "w") as sf:
            for vcf in run_input_files:
                sf.write(vcf + "\n")

        manifest_path = os.path.join(inputs_dir, "input_manifest.txt")
        _write_survivor_input_manifest(
            manifest_path=manifest_path,
            original_files=selected_original_files,
            run_files=run_input_files,
            params=params
        )

    except Exception as e:
        return html.Div(f"❌ Failed to prepare SURVIVOR run inputs: {str(e)}"), None

    merged_filename = f"merged_{timestamp}_output.vcf"
    merged_file_path = os.path.abspath(os.path.join(output_dir, merged_filename))

    survivor_executable = os.path.abspath("./SURVIVOR/Debug/SURVIVOR")

    if not os.access(survivor_executable, os.X_OK):
        try:
            os.chmod(survivor_executable, 0o755)
        except Exception as e:
            return html.Div(f"❌ Cannot make SURVIVOR executable: {str(e)}"), None

    command = [
        survivor_executable,
        "merge",
        sample_files_path,
        str(params["max_distance"]),
        str(params["min_callers"]),
        str(params["type_match"]),
        str(params["strand_match"]),
        str(params["allow_duplicates"]),
        str(params["min_sv_size"]),
        merged_file_path
    ]

    try:
        with open(log_file, "w") as log:
            log.write(f"Run directory: {run_dir}\n")
            log.write(f"Input directory: {inputs_dir}\n")
            log.write(f"Output directory: {output_dir}\n")
            log.write(f"Sample file list: {sample_files_path}\n")
            log.write(f"Running: {' '.join(command)}\n\n")

            result = subprocess.run(
                command,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True
            )

            log.write(result.stdout)

        status_div = html.Div([
            html.Div("✅ SURVIVOR merge completed!"),
            html.Div(f"📁 Run directory: {run_dir}"),
            html.Div(f"📁 Output directory: {output_dir}"),
            html.Div(f"Output saved to: {merged_file_path}"),
            html.Br(),
            html.Div("🗂️ Files merged:"),
            html.Ul([html.Li(os.path.basename(file)) for file in run_input_files]),
            html.Br(),
            warning_msg if warning_msg else ""
        ])

        return status_div, merged_file_path

    except subprocess.CalledProcessError as e:
        return html.Div(f"❌ SURVIVOR failed: {str(e)}"), None

    except Exception as e:
        return html.Div(f"❌ Unexpected error during SURVIVOR merge: {str(e)}"), None
def get_merge_preview(merged_path):
    """Return a DataTable preview of merged VCF."""
    try:
        with open(merged_path, "r") as file:
            for line in file:
                if line.startswith("#CHROM"):
                    columns = line.strip().split("\t")
                    break
            else:
                return html.Div("❌ Error: #CHROM header not found.")

        data_rows = []
        with open(merged_path, "r") as file:
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