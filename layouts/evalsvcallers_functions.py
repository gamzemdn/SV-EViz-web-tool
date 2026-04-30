import os
import subprocess
import base64
import datetime
from dash import html
import time

UPLOAD_DIRECTORY = "./uploaded_files/"
EVAL_OUTPUT_DIR = os.path.join(UPLOAD_DIRECTORY, "evalsvcallers_output")
os.makedirs(EVAL_OUTPUT_DIR, exist_ok=True)

EVAL_UPLOAD_CACHE_DIR = os.path.join(EVAL_OUTPUT_DIR, "upload_cache")
os.makedirs(EVAL_UPLOAD_CACHE_DIR, exist_ok=True)

APP_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
EVALSVCALLERS_DIR = os.path.join(APP_ROOT, "EvalSVcallers-master")

def _decode_upload_content(content):
    content_type, content_string = content.split(",", 1)
    return base64.b64decode(content_string)


def _safe_name(name):
    return os.path.basename(name)


def _detect_extension(path):
    lower = str(path).lower()

    if lower.endswith(".vcf.gz"):
        return ".vcf.gz"
    if lower.endswith(".vcf"):
        return ".vcf"
    if lower.endswith(".bed.gz"):
        return ".bed.gz"
    if lower.endswith(".bed"):
        return ".bed"
    if lower.endswith(".txt"):
        return ".txt"

    return os.path.splitext(path)[1]
def _remove_from_eval_upload_cache(path):
    """
    Remove temporary EvalSVcallers uploaded file after it has been copied
    into the run-specific inputs folder.
    Only files inside evalsvcallers_output/upload_cache are deleted.
    """
    try:
        if not path:
            return

        path = os.path.abspath(path)
        cache_dir = os.path.abspath(EVAL_UPLOAD_CACHE_DIR)

        if os.path.exists(path) and os.path.commonpath([path, cache_dir]) == cache_dir:
            os.remove(path)

    except Exception:
        pass
def save_file(name, content):
    """
    Save uploaded EvalSVcallers files into a temporary upload cache.
    The files are copied into the specific run folder when evaluation starts.
    """
    os.makedirs(EVAL_UPLOAD_CACHE_DIR, exist_ok=True)

    if not name or not content:
        raise ValueError("Uploaded file name or content is missing.")

    data = _decode_upload_content(content)
    safe_name = _safe_name(name)

    path = os.path.abspath(os.path.join(EVAL_UPLOAD_CACHE_DIR, safe_name))

    with open(path, "wb") as fp:
        fp.write(data)

    return path

def save_custom_reference_file(content, original_filename=None):
    """
    Save uploaded custom reference file into EvalSVcallers upload cache.
    It is stored as reference.vcf so the later evaluation step can locate it.
    """
    try:
        os.makedirs(EVAL_UPLOAD_CACHE_DIR, exist_ok=True)

        if not content:
            return "❌ No content to save."

        decoded = _decode_upload_content(content)

        ref_path = os.path.abspath(os.path.join(EVAL_UPLOAD_CACHE_DIR, "reference.vcf"))

        with open(ref_path, "wb") as f:
            f.write(decoded)

        shown_name = original_filename if original_filename else "reference.vcf"

        return (
            f"✅ Reference file uploaded successfully: {shown_name}. "
            f"It will be copied into the run-specific inputs folder during evaluation."
        )

    except Exception as e:
        return f"❌ Failed to save reference file: {str(e)}"
def parse_uploaded_files(contents_list, names_list, dates_list):
    items = []

    for c, n, d in zip(contents_list, names_list, dates_list):
        saved_path = save_file(n, c)
        timestamp = datetime.datetime.fromtimestamp(d)

        items.append(
            html.P([
                f"✅ Uploaded {n}",
                html.Br(),
                html.Small(f"Temporary path: {saved_path}", style={"color": "gray"})
            ])
        )

    return html.Ul(items)
def _resolve_eval_input_file(path_or_name):
    """
    Resolve uploaded/cache/root paths safely.
    """
    if not path_or_name:
        return None

    # absolute path
    if os.path.isabs(path_or_name) and os.path.exists(path_or_name):
        return os.path.abspath(path_or_name)

    # cache path
    cache_path = os.path.abspath(os.path.join(EVAL_UPLOAD_CACHE_DIR, os.path.basename(path_or_name)))
    if os.path.exists(cache_path):
        return cache_path

    # old root path fallback
    root_path = os.path.abspath(os.path.join(EVAL_OUTPUT_DIR, os.path.basename(path_or_name)))
    if os.path.exists(root_path):
        return root_path

    return None


def _copy_input_to_eval_run(src_path, inputs_dir, target_prefix):
    src_path = os.path.abspath(src_path)

    if not os.path.exists(src_path):
        raise FileNotFoundError(f"Input file not found: {src_path}")

    ext = _detect_extension(src_path)
    dst_path = os.path.abspath(os.path.join(inputs_dir, f"{target_prefix}{ext}"))

    if os.path.abspath(src_path) != os.path.abspath(dst_path):
        import shutil
        shutil.copy2(src_path, dst_path)

    return dst_path


def _write_eval_input_manifest(
    manifest_path,
    converted_original,
    converted_run,
    reference_choice,
    reference_original=None,
    reference_run=None,
    parent1_path=None,
    parent2_path=None,
    region_bed_path=None
):
    with open(manifest_path, "w") as f:
        f.write("EvalSVcallers Run Input Manifest\n")
        f.write("================================\n\n")

        f.write("[CONVERTED CALLER VCF]\n")
        f.write(f"Original path: {converted_original}\n")
        f.write(f"Run copy:      {converted_run}\n\n")

        f.write("[REFERENCE]\n")
        f.write(f"Reference choice: {reference_choice}\n")

        if reference_original and reference_run:
            f.write(f"Original path:    {reference_original}\n")
            f.write(f"Run copy:         {reference_run}\n")
        else:
            f.write("No external reference VCF was copied. Embedded reference mode was used.\n")

        f.write("\n[ADVANCED INPUTS]\n")
        f.write(f"Parent 1:   {parent1_path if parent1_path else 'Not used'}\n")
        f.write(f"Parent 2:   {parent2_path if parent2_path else 'Not used'}\n")
        f.write(f"Region BED: {region_bed_path if region_bed_path else 'Not used'}\n")


def _get_run_dir_from_converted_file(converted_file):
    """
    If converted_file is already inside run_xxx/inputs, reuse that run.
    Otherwise return None.
    """
    converted_file = os.path.abspath(converted_file)
    parent_dir = os.path.dirname(converted_file)

    if os.path.basename(parent_dir) == "inputs":
        run_dir = os.path.dirname(parent_dir)
        if os.path.basename(run_dir).startswith("run_"):
            return run_dir

    return None

def run_conversion(caller_tool, uploaded_filename):
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    if not caller_tool or not uploaded_filename:
        return "❌ Caller tool or filename missing.", "", ""

    conversion_script = os.path.join(
        EVALSVCALLERS_DIR,
        "scripts",
        "convert_SV_callers_vcf.pl"
    )

    input_file = _resolve_eval_input_file(uploaded_filename)
    if not input_file:
        return f"❌ Uploaded input file not found: {uploaded_filename}", "", ""

    run_dir = os.path.abspath(os.path.join(EVAL_OUTPUT_DIR, f"run_{timestamp}"))
    inputs_dir = os.path.join(run_dir, "inputs")
    os.makedirs(inputs_dir, exist_ok=True)

    log_file = os.path.join(run_dir, "conversion_log.txt")

    raw_input_run = _copy_input_to_eval_run(
        input_file,
        inputs_dir,
        "raw_caller_input"
    )

    converted_filename = f"converted_{caller_tool}_{timestamp}.vcf"
    converted_file = os.path.abspath(os.path.join(inputs_dir, converted_filename))

    command = f"{conversion_script} -t {caller_tool} '{raw_input_run}' > '{converted_file}'"

    try:
        with open(log_file, "a") as log:
            log.write(f"Command: {command}\n")
            subprocess.run(command, shell=True, check=True, stderr=log)
        
        if os.path.exists(converted_file) and os.path.getsize(converted_file) > 0:
            _remove_from_eval_upload_cache(input_file)
            return (
                html.Div([
                    html.Div("✅ Conversion completed."),
                    html.Div(f"Run directory: {run_dir}"),
                    html.Div(f"Converted file: {converted_file}"),
                    html.Div("Please set the parameters and click Run EvalSVcallers button.")
                ]),
                converted_file,
                timestamp
            )

        return "❌ Conversion failed. Output file is empty.", "", ""

    except Exception as e:
        return f"❌ Error during conversion: {str(e)}", "", ""


def run_evaluation(
    selected_converted_file, reference_choice, caller_tool, timestamp,
    st_basic, l_basic, xl_basic, rl_basic, rxl_basic, mo_basic, of_basic,
    parent1_content, parent1, parent2_content, parent2, rb_content, rb,
    c, mr, rs, mins, eg, eb, i, y, sm
):
    evalsvcallers_script = os.path.join(
        EVALSVCALLERS_DIR,
        "scripts",
        "evaluate_SV_callers.pl"
    )
    if not selected_converted_file:
        return html.Div("❌ No converted file selected."), None, None

    if isinstance(selected_converted_file, list):
        selected_converted_file = selected_converted_file[0]

    converted_file_original = _resolve_eval_input_file(selected_converted_file)

    if not converted_file_original:
        return html.Div(f"❌ Converted file not found: {selected_converted_file}"), None, None

    existing_run_dir = _get_run_dir_from_converted_file(converted_file_original)

    if existing_run_dir:
        run_dir = existing_run_dir
    else:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        run_dir = os.path.abspath(os.path.join(EVAL_OUTPUT_DIR, f"run_{timestamp}"))

    inputs_dir = os.path.join(run_dir, "inputs")
    output_dir = os.path.join(run_dir, "evalsvcallers_output")

    os.makedirs(inputs_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    log_file = os.path.join(run_dir, "evaluation_log.txt")

    # Copy converted file into this run inputs folder if it is not already there
    if os.path.abspath(os.path.dirname(converted_file_original)) == os.path.abspath(inputs_dir):
        converted_file = converted_file_original
    else:
        converted_file = _copy_input_to_eval_run(
            converted_file_original,
            inputs_dir,
            "converted_input"
        )

    if not caller_tool:
        try:
            caller_tool = os.path.basename(converted_file).split("_")[1]
        except Exception:
            caller_tool = "caller"

    base_name = os.path.splitext(os.path.basename(converted_file))[0]

    output_vcf_tp_fp = os.path.abspath(os.path.join(output_dir, f"{base_name}.TF.vcf"))
    output_txt_eval = os.path.abspath(os.path.join(output_dir, f"{base_name}.eval.txt"))

    command = [evalsvcallers_script]

    reference_original = None
    reference_run = None

    # Reference handling
    if reference_choice in ["N", "A"]:
        command.extend(["-r", reference_choice, converted_file])

    elif reference_choice == "univar":
        reference_original = os.path.abspath(os.path.join(UPLOAD_DIRECTORY, "univar.vcf"))

        if not os.path.exists(reference_original):
            return html.Div("❌ Univar SV catalog not found."), None, None

        reference_run = _copy_input_to_eval_run(
            reference_original,
            inputs_dir,
            "reference_univar"
        )

        command.extend(["-r2", reference_run, converted_file])

    else:
        reference_original = os.path.abspath(os.path.join(EVAL_UPLOAD_CACHE_DIR, "reference.vcf"))

        if not os.path.exists(reference_original):
            return html.Div("❌ Reference file not uploaded."), None, None

        reference_run = _copy_input_to_eval_run(
            reference_original,
            inputs_dir,
            "reference_input"
        )

        command.extend(["-r2", reference_run, converted_file])

    # Basic params
    if st_basic:
        command.extend(["-st", st_basic])
    if l_basic is not None:
        command.extend(["-l", str(l_basic)])
    if xl_basic is not None:
        command.extend(["-xl", str(xl_basic)])
    if rl_basic is not None:
        command.extend(["-rl", str(rl_basic)])
    if rxl_basic is not None:
        command.extend(["-rxl", str(rxl_basic)])
    if mo_basic is not None:
        command.extend(["-mo", str(mo_basic)])

    # Advanced uploads are saved directly into run inputs folder
    def decode_and_save_to_inputs(content, filename, prefix):
        if content and filename:
            decoded = _decode_upload_content(content)
            ext = _detect_extension(filename)
            path = os.path.abspath(os.path.join(inputs_dir, f"{prefix}{ext}"))

            with open(path, "wb") as f:
                f.write(decoded)

            return path

        return None

    p1_path = decode_and_save_to_inputs(parent1_content, parent1, "parent1_input")
    p2_path = decode_and_save_to_inputs(parent2_content, parent2, "parent2_input")
    rb_path = decode_and_save_to_inputs(rb_content, rb, "region_bed_input")

    if p1_path:
        command.extend(["-p1", p1_path])
    if p2_path:
        command.extend(["-p2", p2_path])
    if rb_path:
        command.extend(["-rb", rb_path])

    # Advanced numeric/string flags
    if c:
        command.extend(["-c", c])
    if mr is not None:
        command.extend(["-mr", str(mr)])
    if rs is not None:
        command.extend(["-rs", str(rs)])
    if mins is not None:
        command.extend(["-mins", str(mins)])
    if eg:
        command.append("-eg")
    if eb:
        command.append("-eb")
    if i:
        command.append("-i")
    if y:
        command.append("-y")
    if sm:
        command.append("-sm")

    # Output file type
    if of_basic is not None:
        command.extend(["-of", str(of_basic)])
    else:
        command.extend(["-of", "3"])

    manifest_path = os.path.join(inputs_dir, "input_manifest.txt")
    _write_eval_input_manifest(
        manifest_path=manifest_path,
        converted_original=converted_file_original,
        converted_run=converted_file,
        reference_choice=reference_choice,
        reference_original=reference_original,
        reference_run=reference_run,
        parent1_path=p1_path,
        parent2_path=p2_path,
        region_bed_path=rb_path
    )
    _remove_from_eval_upload_cache(converted_file_original)
    
    if reference_original:
        _remove_from_eval_upload_cache(reference_original)
    with open(log_file, "a") as log:
        log.write(f"Run directory: {run_dir}\n")
        log.write(f"Output directory: {output_dir}\n")
        log.write(f"Running: {' '.join(command)}\n")

        try:
            subprocess.run(
                command,
                stderr=log,
                stdout=log,
                cwd=output_dir,
                check=True
            )

            return html.Div([
                html.Div("✅ EvalSVcallers evaluation completed."),
                html.Br(),
                html.Div(f"📁 Run directory: {run_dir}"),
                html.Div(f"📁 Output directory: {output_dir}"),
                html.Div(f"TP/FP File: {output_vcf_tp_fp}"),
                html.Div(f"Metrics File: {output_txt_eval}")
            ]), output_vcf_tp_fp, output_txt_eval

        except Exception as e:
            return html.Div([f"❌ Error: {str(e)}"]), None, None