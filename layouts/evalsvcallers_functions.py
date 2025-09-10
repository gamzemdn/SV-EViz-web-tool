import os
import subprocess
import base64
import datetime
from dash import html
import time

UPLOAD_DIRECTORY = "./uploaded_files/"
EVAL_OUTPUT_DIR = os.path.join(UPLOAD_DIRECTORY, "evalsvcallers_output")
os.makedirs(EVAL_OUTPUT_DIR, exist_ok=True)

def save_file(name, content):
    data = content.encode("utf8").split(b";base64,")[1]
    with open(os.path.join(EVAL_OUTPUT_DIR, name), "wb") as fp:
        fp.write(base64.b64decode(data))
def save_custom_reference_file(content, original_filename=None):
    """
    Save uploaded custom reference file as 'reference.vcf' regardless of original name.
    """
    try:
        if content:
            content_type, content_string = content.split(',')
            decoded = base64.b64decode(content_string)
            ref_path = os.path.join(EVAL_OUTPUT_DIR, "reference.vcf")
            with open(ref_path, "wb") as f:
                f.write(decoded)
            return f"✅ Reference file uploaded successfully as reference.vcf."
        else:
            return "❌ No content to save."
    except Exception as e:
        return f"❌ Failed to save reference file: {str(e)}"
def parse_uploaded_files(contents_list, names_list, dates_list):
    items = []
    for c, n, d in zip(contents_list, names_list, dates_list):
        save_file(n, c)
        timestamp = datetime.datetime.fromtimestamp(d)
        items.append(html.P(f"✅ Uploaded {n}"))
    return html.Ul(items)

def run_conversion(caller_tool, uploaded_filename):
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    if not caller_tool or not uploaded_filename:
        return "❌ Caller tool or filename missing.", "", ""

    conversion_script = os.path.expanduser("~/EvalSVcallers-master/scripts/convert_SV_callers_vcf.pl")
    log_file = os.path.join(EVAL_OUTPUT_DIR, "conversion_log.txt")
    input_file = os.path.abspath(os.path.join(EVAL_OUTPUT_DIR, uploaded_filename))
    converted_filename = f"converted_{caller_tool}_{timestamp}.vcf"
    converted_file = os.path.abspath(os.path.join(EVAL_OUTPUT_DIR, converted_filename))

    command = f"{conversion_script} -t {caller_tool} '{input_file}' > '{converted_file}'"
    try:
        with open(log_file, "a") as log:
            log.write(f"Command: {command}\n")
            subprocess.run(command, shell=True, check=True)
        if os.path.exists(converted_file) and os.path.getsize(converted_file) > 0:
            return f"✅ Conversion completed. Please set the parameters and click Run Evalsvcallers button! Converted file: {converted_file}", converted_filename,timestamp
        else:
            return f"❌ Conversion failed. Output file is empty.", "", ""
    except Exception as e:
        return f"❌ Error during conversion: {str(e)}", "", ""

def run_evaluation(
    selected_converted_file, reference_choice, caller_tool,timestamp,
    st_basic, l_basic, xl_basic, rl_basic, rxl_basic, mo_basic, of_basic,
    parent1_content, parent1, parent2_content, parent2, rb_content, rb,
    c, mr, rs, mins, eg, eb, i, y, sm
):
   # timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    evalsvcallers_script = os.path.expanduser("~/EvalSVcallers-master/scripts/evaluate_SV_callers.pl")
    log_file = os.path.join(EVAL_OUTPUT_DIR, "evaluation_log.txt")

    if not selected_converted_file:
        return dash.no_update, []

    if isinstance(selected_converted_file, list):
        selected_converted_file = selected_converted_file[0]

    converted_file = os.path.abspath(selected_converted_file)

    if not caller_tool:
        caller_tool = os.path.basename(selected_converted_file).split("_")[1]

    if timestamp:
        output_vcf_tp_fp = os.path.join(EVAL_OUTPUT_DIR, f"converted_{caller_tool}_{timestamp}.TF.vcf")
        output_txt_eval = os.path.join(EVAL_OUTPUT_DIR, f"converted_{caller_tool}_{timestamp}.eval.txt")
    else:
        # fallback for manually uploaded converted file
        base_name = os.path.splitext(os.path.basename(selected_converted_file))[0]
        output_vcf_tp_fp = os.path.join(EVAL_OUTPUT_DIR, f"{base_name}.TF.vcf")
        output_txt_eval = os.path.join(EVAL_OUTPUT_DIR, f"{base_name}.eval.txt")

    command = [evalsvcallers_script]

    # Reference handling
    if reference_choice in ["N", "A"]:
        command.extend(["-r", reference_choice, converted_file])
    elif reference_choice == "univar":
        ref_file = os.path.abspath(os.path.join(UPLOAD_DIRECTORY, "univar.vcf"))
        if not os.path.exists(ref_file):
            return html.Div("❌ Univar SV katalog not found."), []
        command.extend(["-r2", ref_file, converted_file])
    else:
        ref_file = os.path.abspath(os.path.join(EVAL_OUTPUT_DIR, "reference.vcf"))
        if not os.path.exists(ref_file):
            return html.Div("❌ Reference file not uploaded."), []
        command.extend(["-r2", ref_file, converted_file])

    # Basic params
    if st_basic: command.extend(["-st", st_basic])
    if l_basic is not None: command.extend(["-l", str(l_basic)])
    if xl_basic is not None: command.extend(["-xl", str(xl_basic)])
    if rl_basic is not None: command.extend(["-rl", str(rl_basic)])
    if rxl_basic is not None: command.extend(["-rxl", str(rxl_basic)])
    if mo_basic is not None: command.extend(["-mo", str(mo_basic)])

    # Advanced param uploads
    def decode_and_save(content, filename):
        if content and filename:
            content_type, content_string = content.split(',')
            decoded = base64.b64decode(content_string)
            path = os.path.join(EVAL_OUTPUT_DIR, filename)
            with open(path, "wb") as f:
                f.write(decoded)
            return path
        return None

    p1_path = decode_and_save(parent1_content, parent1)
    p2_path = decode_and_save(parent2_content, parent2)
    rb_path = decode_and_save(rb_content, rb)

    if p1_path: command.extend(["-p1", p1_path])
    if p2_path: command.extend(["-p2", p2_path])
    if rb_path: command.extend(["-rb", rb_path])

    # Advanced numeric/string flags
    if c: command.extend(["-c", c])
    if mr is not None: command.extend(["-mr", str(mr)])
    if rs is not None: command.extend(["-rs", str(rs)])
    if mins is not None: command.extend(["-mins", str(mins)])
    if eg: command.append("-eg")
    if eb: command.append("-eb")
    if i:  command.append("-i")
    if y:  command.append("-y")
    if sm: command.append("-sm")
    
    # Output file type (-of) only uses of_basic now
    if of_basic is not None:
        command.extend(["-of", str(of_basic)])
    else:
        command.extend(["-of", "3"])

    # Final log and run
    with open(log_file, "a") as log:
        log.write(f"Running: {' '.join(command)}\n")
        try:            
            subprocess.run(command, stderr=log, cwd=EVAL_OUTPUT_DIR, check=True)
            return html.Div([
                "✅ EvalSVcallers evaluation completed.",
                html.Br(),
                html.Div(f"TP/FP File: {output_vcf_tp_fp}"),
                html.Div(f"Metrics File: {output_txt_eval}")
            ]),  output_vcf_tp_fp, output_txt_eval
        except Exception as e:
            return html.Div([f"❌ Error: {str(e)}"]), []