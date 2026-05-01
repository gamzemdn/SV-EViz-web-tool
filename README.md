# SV-EViz: Structural Variant Evaluation and Visualization Web Tool

SV-EViz is a web-based platform for structural variant (SV) comparison, benchmarking, metric inspection, and visualization. It provides a graphical interface for commonly used SV analysis workflows and reduces the need to run multiple command-line tools manually.

SV-EViz integrates:

- **SURVIVOR** for multi-caller SV consensus/merge generation
- **Truvari** for query-versus-truth SV benchmarking
- **EvalSVcallers** for caller-specific SV conversion and evaluation
- Built-in visualization modules for SV summaries and genome-wide inspection

A lightweight public version is available for interface testing:

**Public Render-hosted web version:**  
https://sv-eviz-web-tool.onrender.com/

> Note: The Render-hosted version is intended for lightweight testing and user-interface exploration. For large whole-genome VCF, BED, or FASTA files, local Docker or server-based deployment is recommended.

---

## Main Features

- Web-based graphical interface for SV comparison and visualization
- Support for three major SV workflows:
  - SURVIVOR merge
  - Truvari bench
  - EvalSVcallers conversion and evaluation
- Module-specific input handling
- Example datasets for quick testing
- Interactive visualizations:
  - SVTYPE distributions
  - SV length distributions
  - Manhattan plots
  - Clustergrams
  - Circos plots
  - Sankey diagrams
- Adjustable parameters for benchmarking modules
- Run-specific output folders
- Downloadable ZIP archives containing generated outputs
- Exported tables in CSV and Excel formats
- Exported Plotly figures as interactive HTML and, when supported, static PNG files

---

## Deployment Options

| Deployment option | Description | Recommended use |
|---|---|---|
| Public Render-hosted version | Browser-accessible version without local installation | Lightweight testing and interface exploration |
| Docker local deployment | Runs SV-EViz locally using Docker | Recommended for most users and reproducible testing |
| Server-based deployment | Runs SV-EViz on an institutional or cloud server | Recommended for larger datasets and shared use |

---

## Technical Requirements

For local or server-based use, the following are recommended:

| Requirement | Recommended configuration |
|---|---|
| Operating system | Linux or Docker-supported environment |
| Python | Python 3.11 inside Docker |
| Memory | At least 8 GB RAM for small to medium examples; more for large WGS files |
| Disk space | Depends on uploaded VCF, BED, FASTA, and generated output files |
| Browser | Modern browser such as Chrome, Firefox, or Edge |
| Reference FASTA | Required for Truvari benchmarking |

Large reference FASTA files are not included in the repository because of their size. For Truvari analyses, users should place reference FASTA files in:

```bash
uploaded_files/reference_files/
