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
os.makedirs(SURVIVOR_OUTPUT_DIR, exist_ok=True)

def save_file(name, content):
    """Decode and save uploaded file to survivor folder."""
    data = content.encode("utf8").split(b";base64,")[1]
    with open(os.path.join(SURVIVOR_OUTPUT_DIR, name), "wb") as fp:
        fp.write(base64.b64decode(data))
        
        
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
        content_type, content_string = content.split(',')
        decoded = base64.b64decode(content_string)

        base = os.path.basename(name)
        root, ext = os.path.splitext(base)
        save_path = os.path.join(SURVIVOR_OUTPUT_DIR, base)

        # Ensure unique filename
        if os.path.exists(save_path):
            stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
            save_path = os.path.join(SURVIVOR_OUTPUT_DIR, f"{root}_{stamp}{ext}")

        with open(save_path, 'wb') as f:
            f.write(decoded)

        saved_paths.append(save_path)
        messages.append(html.P(f"✅ Uploaded: {os.path.basename(save_path)}", style={}))

    return html.Div(messages), saved_paths

def run_survivor_merge(params, vcf_files):
    import subprocess
    import datetime
    import os

    # Ensure output directory and filenames
    os.makedirs(SURVIVOR_OUTPUT_DIR, exist_ok=True)
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    merged_filename = f"merged_{timestamp}_output.vcf"
    merged_file_path = os.path.join(SURVIVOR_OUTPUT_DIR, merged_filename)
    sample_files_path = os.path.join(SURVIVOR_OUTPUT_DIR, "sample_files")

   # --- Normalize & keep last 3 unique (by most recent in list) ---
 #   vcf_files = [os.path.abspath(p) for p in (vcf_files or [])]
 #   unique_files = list(dict.fromkeys(reversed(vcf_files)))[::-1][-3:]

    # --- Enforce minimum 2 distinct files ---
#    if len(unique_files) < 2:
#        msg = html.Div([
 #           "❌ Error: At least two caller VCF files are required.",
 #           html.Br(),
  #          html.Div("Current files:"),
#            html.Ul([html.Li(os.path.basename(p)) for p in unique_files]) if unique_files else html.Div("— none —")
  #      ])
 #       return msg, None

            
    # Ensure minimum input requirement
 #   if len(vcf_files) < 2:
  #      return "❌ Error: At least two VCF files are required to merge.", None

    # Write file paths to sample list
    try:
       
        with open(sample_files_path, "w") as sf:
            for vcf in vcf_files[:3]:
                sf.write(vcf + "\n")
    except Exception as e:
        return f"❌ Failed to write sample file list: {str(e)}", None

    # Ensure SURVIVOR executable has permission
    survivor_executable = "./SURVIVOR/Debug/SURVIVOR"
    if not os.access(survivor_executable, os.X_OK):
        try:
            os.chmod(survivor_executable, 0o755)
        except Exception as e:
            return f"❌ Cannot make SURVIVOR executable: {str(e)}", None

    # Construct merge command
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
        warning_msg = None
        if len(vcf_files) > 3:
            warning_msg = html.Div(f"⚠️ Only the first 3 of {len(vcf_files)} uploaded files will be merged.")
    
        with open(sample_files_path, "w") as sf:
            for vcf in vcf_files[:3]:
                sf.write(vcf + "\n")
    except Exception as e:
        return f"❌ Failed to write sample file list: {str(e)}", None
    
    # Run SURVIVOR
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
        #return f"✅ SURVIVOR merge completed! Output saved to {merged_file_path}", merged_file_path
        
        with open(sample_files_path, "r") as f:
            file_lines = [line.strip() for line in f if line.strip()]

        status_div = html.Div([
            html.Div("✅ SURVIVOR merge completed!"),
            html.Div(f"Output saved to: {merged_file_path}"),
            html.Br(),
            html.Div("🗂️ Files merged:"),
            html.Ul([html.Li(os.path.basename(file)) for file in vcf_files[:3]]),
            html.Br(),
            warning_msg if warning_msg else ""
        ])

        return status_div, merged_file_path
    
    except subprocess.CalledProcessError as e:
        return f"❌ SURVIVOR failed:\n{e.stderr.strip()}", None
    except Exception as e:
        return f"❌ Unexpected error during SURVIVOR merge: {str(e)}", None
    
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