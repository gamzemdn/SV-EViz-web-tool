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
```

For Docker use, this folder is mounted into the container through `docker-compose.yml`.

---

## Reference Genome Setup (Required for Truvari)

Truvari benchmarking requires a reference FASTA file matching the genome
build of your VCF files. Due to the file size (~3 GB per genome build),
reference genomes are **not included** in this repository or the Docker
image. You must download them locally before running Truvari.

### Where to place reference files

Reference FASTA files (and their `.fai` indices) must be placed in:

```bash
uploaded_files/reference_files/
```

This directory is bind-mounted into the Docker container by
`docker-compose.yml`, so any FASTA you place here becomes immediately
available to Truvari without rebuilding the image. The application
auto-detects all `.fa` files in this directory and lists them in the
Truvari "Reference Genome" dropdown.

### Downloading GRCh38 (hg38)

```bash
cd uploaded_files/reference_files/
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
mv hg38.fa GRCh38.fa
samtools faidx GRCh38.fa
```

### Downloading GRCh37 / hs37d5

```bash
cd uploaded_files/reference_files/
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz
samtools faidx hs37d5.fa
```

### Verifying the setup

After downloading, you should see:

```bash
$ ls uploaded_files/reference_files/
GRCh38.fa  GRCh38.fa.fai
# or
hs37d5.fa  hs37d5.fa.fai
```

When you open the Truvari tab in the web interface, the reference dropdown
should now list the available FASTA file(s). If the dropdown is empty,
verify the files are placed directly in `uploaded_files/reference_files/`
(not in a subdirectory) and have the `.fa` extension.

### Why is the reference required?

The Truvari pipeline uses `bcftools reheader -f <ref.fai>` to attach
sequence dictionaries to base and comparison VCFs before sorting and
benchmarking. Without a valid reference FASTA index, the reheader step
fails and the sort step exits with code 255. Other modules of the
application (SURVIVOR merging, EvalSVcallers, visualization) do **not**
require the reference and will work without it.

### Note on Docker deployment

The reference FASTA files are not stored inside the Docker image.

For Docker Compose, the included `docker-compose.yml` mounts
`./uploaded_files` from the host into the container. Place the reference files in:

```text
uploaded_files/reference_files/

---
Then start the application with:

```bash
sudo docker compose up
```

The included `docker-compose.yml` mounts the local `uploaded_files` directory
into the container.

When starting the published Docker image directly, mount the local reference
directory into the container:

```bash
sudo docker run -d \
  --name sv-eviz \
  -p 8040:8040 \
  -v "$(pwd)/uploaded_files/reference_files:/app/uploaded_files/reference_files:ro" \
  maden21/sv-eviz:latest
```

Without this volume mount, the reference directory inside the container will
be empty and no reference genome will appear in the Truvari dropdown.

## Installation with Docker Compose

Clone the repository:

```bash
git clone https://github.com/gamzemdn/SV-EViz-web-tool.git
cd SV-EViz-web-tool
```

Build the Docker image:

```bash
sudo docker compose build --no-cache
```

Run the application:

```bash
sudo docker compose up
```

Open the application in a browser:

```text
http://localhost:8040
```

To stop the application:

```bash
sudo docker compose down
```

---

## Manual Installation

Manual installation is possible, but Docker is recommended because SV-EViz depends on external tools and system-level packages.

```bash
git clone https://github.com/gamzemdn/SV-EViz-web-tool.git
cd SV-EViz-web-tool
pip install -r requirements.txt
python app.py
```

---

## Example Data

Example input files are provided under the `example_data` directory to allow quick testing after installation.

The example data are organized by module:

```text
example_data/
├── comparison/
│   ├── survivor/
│   ├── truvari/
│   └── evalsvcallers/
├── visualization/
│   ├── caller-truvari/
│   ├── survivor/
│   └── evalsvcallers/
└── metrics/
```

These files allow users to test:

- SURVIVOR merge using multiple caller VCF files
- Truvari benchmarking using query/truth VCF files and optional BED input
- EvalSVcallers conversion and evaluation
- Visualization of caller-like VCF, SURVIVOR output, Truvari-compatible VCF, and EvalSVcallers-style output
- Metric inspection using Truvari `summary.json` , `summary.txt` and EvalSVcallers `.eval.txt` outputs

Depending on the module, users can either upload their own files or use the provided example datasets to test the workflow.

---

## Supported Inputs by Module

| Module | Required input | Minimal required fields | Output |
|---|---|---|---|
| SURVIVOR | Multiple caller VCF files | CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO with SVTYPE | Merged/consensus VCF |
| Truvari | BASE/truth VCF, COMP/query VCF, reference FASTA, optional BED | Sorted/indexed VCF-compatible SV records; reference FASTA | TP, FP, FN VCF outputs and `summary.json` |
| EvalSVcallers | Caller VCF and selected reference truth set | Caller-compatible VCF with SVTYPE and SV coordinates | Converted VCF, TP/FP output, `.eval.txt` metrics |
| Visualization | Caller-like, SURVIVOR, Truvari-compatible, or EvalSVcallers-style VCF files | SVTYPE in INFO for type-based plots; SVLEN or END for length-based plots | Interactive visual summaries |
| Metrics | Truvari `summary.json` or EvalSVcallers `.eval.txt` | Tool-specific metric output file | Parsed tables and metric plots |

---

## Important Input Notes

- SVTYPE should be present in the INFO field for full type-based visualization.
- SVLEN or END is recommended for length-based plots.
- Generic VCF files may work in the Visualization tab but may not be valid for Comparison or Metrics workflows.
- The Comparison tab requires workflow-specific inputs, such as:
  - multiple caller VCFs for SURVIVOR,
  - query/truth VCF pairs for Truvari or EvalSVcallers.
- The Metrics tab requires pre-computed metric files, not generic VCF files.

---

## Output Organization

SV-EViz stores generated files in run-specific output folders under `uploaded_files/`.

A typical run folder may contain:

```text
inputs/
tables/
figures/
downloads/
runtime_memory_info.txt
```

Generated outputs include:

- Uploaded input files used for the run
- Tool-generated output files
- Parsed tables
- CSV and Excel exports
- Interactive HTML figures
- Static PNG figures when supported
- ZIP archive containing the run outputs

This structure prevents outputs from different analyses from being overwritten and improves reproducibility.

---

## Visualization Modules

SV-EViz provides the following visualization options:

- Data table and basic visual summaries
- SVTYPE distribution plots
- SV length distribution plots
- Chromosome-wise SV distributions
- Sankey diagrams
- Clustergrams
- Circos plots
- Manhattan plots

The Manhattan plot supports QUAL-based genome-wide visualization. The displayed score is calculated as:

```text
QUAL / 200
```

which is equivalent to:

```text
−log10(P), where P = 10^(−QUAL / 200)
```

Highlighted points indicate variants exceeding the selected threshold.

---

## Workflow Comparison: CLI vs SV-EViz

| Task | Native CLI workflow | SV-EViz workflow |
|---|---|---|
| SURVIVOR merge | Prepare file list, define parameters, run command, locate output manually | Upload VCFs, select parameters, run from interface |
| Truvari benchmarking | Prepare query/truth VCFs, provide reference FASTA, define parameters, run command manually | Select files, provide reference FASTA, optional BED, run from interface |
| EvalSVcallers benchmarking | Select caller-specific converter, run conversion script, run evaluation, parse metrics manually | Select caller type, upload files, run conversion/evaluation, view parsed metrics |
| Visualization | Write or reuse separate plotting scripts | Use built-in visualization tab |
| Metric summary | Parse tool-specific output files manually | View harmonized tables and plots |

---

## Case Study and Validation

The interface examples were generated using representative germline WGS structural variant callsets. SV-EViz was tested across multiple datasets and workflows, including HG00514, HG002, and NA12878 examples.

The validation included:

- File upload
- Backend workflow execution
- Output generation
- Metric parsing
- Visualization generation
- Downloadable output packaging

These tests were used to confirm that the workflow is not limited to a single dataset or caller representation.

---

## Repository Structure

```text
SV-EViz-web-tool/
├── app.py
├── Dockerfile
├── docker-compose.yml
├── requirements.txt
├── assets/
├── layouts/
├── example_data/
├── uploaded_files/
│   ├── reference_files/
│   ├── survivor_output/
│   ├── evalsvcallers_output/
│   ├── truvari_output/
│   ├── visualization_output/
│   └── metrics_output/
├── EvalSVcallers-master/
└── SURVIVOR/
```

---

## Notes on Large Files

Large FASTA, BAM, FASTQ, and WGS-scale VCF files should not be committed to the GitHub repository. For local or Docker-based analyses, place reference FASTA files inside:

```bash
uploaded_files/reference_files/
```

For example:

```text
uploaded_files/reference_files/GRCh38.fa
uploaded_files/reference_files/GRCh38.fa.fai
uploaded_files/reference_files/hs37d5.fa
uploaded_files/reference_files/hs37d5.fa.fai
```

If the FASTA index file is not available, SV-EViz attempts to create it during Truvari execution when the required tools are available.

---

## Troubleshooting

### Port 8040 is already in use

If Docker cannot start because port 8040 is already occupied, stop the existing process or change the exposed port in `docker-compose.yml`.

To identify the process using port 8040:

```bash
sudo lsof -i :8040
```

To release the port:

```bash
sudo fuser -k 8040/tcp
```

Then restart:

```bash
sudo docker compose up
```

### Reference FASTA not listed in Truvari

Make sure the FASTA file and its `.fai` index are placed directly inside:

```bash
uploaded_files/reference_files/
```

If running with Docker Compose, this folder is mounted into the container.

### PNG export does not work

Static PNG export requires Plotly/Kaleido support and a Chromium-compatible browser inside the environment. The Docker image installs Chromium and sets the required environment variables.

---

## Citations for Integrated Frameworks

SV-EViz integrates existing SV tools and provides a graphical workflow around them. Please cite the original tools when using their corresponding modules.

### SURVIVOR

Jeffares, D.C., Jolly, C., Hoti, M., Speed, D., Shaw, L., Rallis, C., Balloux, F., Dessimoz, C., Bähler, J., and Sedlazeck, F.J. (2017).  
“Transient structural variations have strong effects on quantitative traits and reproductive isolation in fission yeast.”  
*Nature Communications*, 8, 14061.  
https://doi.org/10.1038/ncomms14061  

GitHub: https://github.com/fritzsedlazeck/SURVIVOR

### EvalSVcallers

Kosugi, S., Momozawa, Y., Liu, X., Terao, C., Kubo, M., and Kamatani, Y. (2019).  
“Comprehensive evaluation of structural variation detection algorithms for whole genome sequencing.”  
*Genome Biology*, 20, 117.  
https://doi.org/10.1186/s13059-019-1720-5  

GitHub: https://github.com/stat-lab/EvalSVcallers

### Truvari

English, A.C., Menon, V.K., Gibbs, R.A., Metcalf, G.A., Sedlazeck, F.J., and others. (2022).  
“Truvari: refined structural variant comparison preserves allelic diversity.”  
*Genome Biology*, 23, 271.  
https://doi.org/10.1186/s13059-022-02840-6  

GitHub: https://github.com/ACEnglish/truvari

---

## License

Please see the repository license file for usage terms.
