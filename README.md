# INCI: Insecticidal Cross-kingdom RNA Interference pipeline

INCI is a research pipeline for identifying candidate double-stranded RNAs, examining small RNA coverage, predicting RNA targets, and evaluating degradome support. It combines Python analysis modules with a browser interface that runs locally on your computer.

This repository contains a research software snapshot. Interfaces, methods, and output formats may change. Predicted dsRNA loci, target interactions, and adaptation scores require biological interpretation and experimental validation.

## Branches

- **`_cur`** is the default branch and presents the curated set of working tools.
- **`_dev`** preserves the broader development snapshot, including experimental tools, integrations, and unfinished sections. Availability in this branch does not mean a tool has been fully validated.

Clone the default branch using the instructions below. To work on the development snapshot, run `git switch _dev` after cloning.

## Development tool coverage

- **Read preprocessing:** batch adapter and quality trimming with Trim Galore or Cutadapt.
- **dsRNA analysis:** bidirectional RNA-seq coverage, candidate locus ranking, and grouped replicate coverage plots.
- **Small RNA analysis:** Bowtie 1 mapping, coverage and length distributions, candidate dsRNA regions, and control-reference filtering.
- **Target analysis:** RNAplex-based pairing and MFE-ratio filtering, CleaveLand-style degradome analysis, and siRNA adaptation comparisons across transcriptomes.
- **Supporting code and tools:** FASTA deduplication, BLAST helper functions, siRNA generation and mismatch scoring, and optional dsRIP integration.

The list above describes `_dev`. The BLAST Tool, RNAi Susceptibility Prediction, and Orthology Inference pages are placeholders without runnable interface workflows. The curated branch excludes those pages, siRNA adaptation, dsRIP enhancement, and other unfinished interface tools.

The degradome implementation is INCI's own CleaveLand-style workflow; the CleaveLand4 distribution is not required or bundled.

## Install and run

Use Python 3.10 or later. Run the following from a terminal on macOS or Linux:

```sh
git clone https://github.com/dogapm123/INCI-Insecticidal-Cross-kingdom-RNA-Interference-pipeline.git
cd INCI-Insecticidal-Cross-kingdom-RNA-Interference-pipeline
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python inci_main.py
```

The launcher opens your browser and prints a local address such as `http://127.0.0.1:<port>/`. It selects an available port each time. Keep the terminal process running while using the interface; press `Ctrl+C` to stop it.

Create or select a project, provide your own input files and references, and check the external tools in **Settings** before starting an analysis. File and folder dialogs currently use macOS integration; paths can also be entered in the interface.

Python requirements cover the core pipeline and its supporting scripts. They are not a locked, fully reproduced analysis environment. External executables and optional dsRIP dependencies must be installed separately.

## External executables

Install the tools needed by the workflows you use and make their executables available on `PATH`. The interface also checks several local installation locations.

| Tool | Used for |
| --- | --- |
| Trim Galore or Cutadapt | Sequencing-read preprocessing |
| Bowtie 1 (`bowtie`, `bowtie-build`) | Small RNA mapping, control filtering, and mapping-based adaptation |
| minimap2, SAMtools, and SeqKit | RNA-seq mapping, genomic reference binning, and alignment processing |
| ViennaRNA (`RNAplex`, `RNAfold`) | Target interaction energies and RNA structure calculations |
| NCBI BLAST+ (`blastn`, `makeblastdb`) | Local sequence similarity searches |
| CD-HIT (`cd-hit-est`) | Nucleotide FASTA deduplication |
| BEDTools | Supporting interval and sequence-extraction workflows |

Bowtie 2 does not replace Bowtie 1 for these workflows. The Python dependency file does not install the programs in this table.

## Optional dsRIP integration in `_dev`

The adapter in `dsrip_sirna_api.py` expects a separate dsRIP installation at:

```text
dsRIP/dsRIP_web/main_site/
```

The dsRIP source, reference databases, and constant files are not included in this repository. The full efficiency workflow needs that installation and its Python dependencies, including ViennaRNA bindings (`RNA`), Biopython, RNAtweaks, and orffinder, as well as its relevant external tools and reference assets.

The single-siRNA feature adapter also uses the ViennaRNA Python bindings and dsRIP's thermodynamic lookup database (`constant_files/lookup.db`). The current adapter can omit thermodynamic lookup contributions when that database is absent, so provide the appropriate dsRIP assets before interpreting those scores. Installing the core requirements alone does not enable dsRIP-based analyses.

## Inputs and outputs

Supply your own FASTA references and sequencing files. Supported inputs depend on the module and include FASTQ, gzipped FASTQ, FASTA, SAM/BAM, and CSV/TSV sample manifests. Consult the selected tool or its `--help` output for the exact format and parameters.

For example, the dsRNA plotter accepts a paired-end sample CSV with these columns:

```csv
sample_id,group,replicate,read1,read2
sample_1,condition_1,1,/path/to/sample_1_R1.fastq.gz,/path/to/sample_1_R2.fastq.gz
```

Project outputs are saved under `outputs/<project>/` by default or in the project directory you select. Local paths and interface preferences are stored in `.inci_pipeline_paths.json`.

Raw sequencing data, reference collections, generated results, alignment indexes, local environments, and machine-specific settings are excluded from this repository. The interface's predefined control-reference shortcuts require local reference assets; use your own reference FASTA files in a fresh installation.

## Command-line use

Several analyses can run independently of the browser interface. Explore their options from the repository root:

```sh
python dsRNA_identification.py --help
python dsRNA_plotter.py --help
python sRNA_identification.py --help
python MFE_ratio.py --help
```

Additional development commands include:

```sh
python dsrna_adaptation.py --help
python bowtie_multitranscriptome_adaptation.py --help
```

## Development checks

Run the Python regression checks in the activated environment:

```sh
python -m pip install pytest
python -m pytest tests
```

The interface polling check additionally requires Node.js and uses `python3` from your active environment:

```sh
node tests/interface_polling.cjs
```

These checks cover selected analysis and background-job behavior. They do not constitute end-to-end validation of every workflow or external tool installation.
