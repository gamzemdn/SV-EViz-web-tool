import os
import re
import zipfile
import pandas as pd
import plotly.graph_objects as go
from dash import html
from dash import dcc
from dash import dash_table


def _safe_filename(text, fallback):
    if not text:
        text = fallback

    text = str(text)
    text = re.sub(r"<[^>]+>", "", text)
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    text = text.strip("_")

    if not text:
        text = fallback

    return text[:120]


def _iter_dash_components(obj):
    """
    Recursively walk through Dash components/lists/tuples.
    """
    if obj is None:
        return

    if isinstance(obj, (list, tuple)):
        for item in obj:
            yield from _iter_dash_components(item)
        return

    yield obj

    children = getattr(obj, "children", None)

    if children is not None:
        yield from _iter_dash_components(children)


def _extract_plotly_title(fig):
    try:
        title = fig.layout.title.text
        if title:
            return title
    except Exception:
        pass

    return None
def _extract_graph_name(comp, fig, graph_counter):
    """
    Create a meaningful filename for exported Plotly/Dash Bio figures.
    Priority:
    1. Plot title
    2. dcc.Graph id
    3. fallback figure number
    """
    title = _extract_plotly_title(fig)
    if title:
        return title

    graph_id = getattr(comp, "id", None)
    if graph_id:
        return str(graph_id)

    return f"figure_{graph_counter:03d}"

def _save_graphs_from_components(components, figures_dir):
    os.makedirs(figures_dir, exist_ok=True)

    saved = []
    graph_counter = 1

    for comp in _iter_dash_components(components):
        if isinstance(comp, dcc.Graph):
            fig_obj = getattr(comp, "figure", None)

            if fig_obj is None:
                continue

            try:
                fig = go.Figure(fig_obj)
            except Exception:
                continue

            graph_name = _extract_graph_name(comp, fig, graph_counter)
            base_name = _safe_filename(graph_name, f"figure_{graph_counter:03d}")

            html_path = os.path.join(figures_dir, f"{graph_counter:03d}_{base_name}.html")
            fig.write_html(html_path)

            saved.append(html_path)

            # PNG requires kaleido. If unavailable, do not break the app.
            png_path = os.path.join(figures_dir, f"{graph_counter:03d}_{base_name}.png")
            try:
                fig.write_image(png_path, scale=2)
                saved.append(png_path)
            except Exception:
                pass

            graph_counter += 1

    return saved


def _save_tables_from_components(components, tables_dir):
    os.makedirs(tables_dir, exist_ok=True)

    saved = []
    table_counter = 1
    excel_path = os.path.join(tables_dir, "tables.xlsx")

    table_dfs = []

    for comp in _iter_dash_components(components):
        if isinstance(comp, dash_table.DataTable):
            data = getattr(comp, "data", None)
            columns = getattr(comp, "columns", None)

            if not data:
                continue

            try:
                df = pd.DataFrame(data)

                if columns:
                    ordered_cols = [
                        col.get("id") for col in columns
                        if isinstance(col, dict) and col.get("id") in df.columns
                    ]
                    if ordered_cols:
                        df = df[ordered_cols]

                table_name = f"table_{table_counter:03d}"
                csv_path = os.path.join(tables_dir, f"{table_name}.csv")
                df.to_csv(csv_path, index=False)

                saved.append(csv_path)
                table_dfs.append((table_name, df))

                table_counter += 1

            except Exception:
                continue

    if table_dfs:
        try:
            with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
                for sheet_name, df in table_dfs:
                    df.to_excel(writer, sheet_name=sheet_name[:31], index=False)

            saved.append(excel_path)
        except Exception:
            pass

    return saved


def _zip_folder(source_dir, zip_path):
    os.makedirs(os.path.dirname(zip_path), exist_ok=True)

    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(source_dir):
            for file in files:
                abs_path = os.path.join(root, file)

                # Do not put the zip inside itself
                if os.path.abspath(abs_path) == os.path.abspath(zip_path):
                    continue

                rel_path = os.path.relpath(abs_path, source_dir)
                zipf.write(abs_path, rel_path)

    return zip_path


def create_run_export_package(run_dir, components, zip_prefix="exports"):
    """
    Save all visible Dash DataTables and Plotly Graphs from a result block
    into the run folder, then create a zip package.

    Returns:
        zip_path, saved_files
    """
    run_dir = os.path.abspath(run_dir)

    tables_dir = os.path.join(run_dir, "tables")
    figures_dir = os.path.join(run_dir, "figures")
    downloads_dir = os.path.join(run_dir, "downloads")

    os.makedirs(tables_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(downloads_dir, exist_ok=True)

    saved_files = []

    saved_files.extend(_save_tables_from_components(components, tables_dir))
    saved_files.extend(_save_graphs_from_components(components, figures_dir))

    run_name = os.path.basename(run_dir)
    zip_name = f"{run_name}_{zip_prefix}.zip"
    zip_path = os.path.join(downloads_dir, zip_name)

    _zip_folder(run_dir, zip_path)

    return zip_path, saved_files


def make_zip_download_card(tool_name, run_dir, zip_path, saved_files=None):
    """
    Create a small UI card containing the zip download link.
    The Flask route must be added in app.py.
    """
    run_name = os.path.basename(os.path.abspath(run_dir))
    zip_name = os.path.basename(zip_path)

    href = f"/download/run_zip/{tool_name}/{run_name}/{zip_name}"

    saved_count = len(saved_files) if saved_files else 0

    return html.Div([
        html.H5("Downloadable Outputs", style={"marginTop": "15px"}),
        html.Div(
            f"Saved files: {saved_count}. Tables and figures were saved inside the run folder.",
            style={"fontSize": "13px", "color": "gray"}
        ),
        html.A(
            "⬇ Download ZIP package",
            href=href,
            target="_blank",
            style={
                "display": "inline-block",
                "marginTop": "8px",
                "padding": "8px 12px",
                "backgroundColor": "#0984e3",
                "color": "white",
                "borderRadius": "6px",
                "textDecoration": "none",
                "fontFamily": '"Times New Roman", Times, serif'
            }
        )
    ], style={
        "border": "1px solid #ddd",
        "padding": "10px",
        "borderRadius": "8px",
        "marginTop": "15px",
        "marginBottom": "15px",
        "backgroundColor": "#f8f9fa"
    })