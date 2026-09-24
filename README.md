# INCI: Insecticidal Cross-kingdom RNA Interference pipeline

INCI is a research pipeline for identifying candidate double-stranded RNAs, examining small RNA coverage, predicting RNA targets, and evaluating degradome support. It combines Python analysis modules with a browser interface that runs locally on your computer.

This repository contains a research software snapshot. Interfaces, methods, and output formats may change. Predicted dsRNA loci and target interactions require biological interpretation and experimental validation.

## Branches

- **`_cur`** is the default branch and contains the nine tools listed below.
- **`_dev`** preserves the broader development snapshot, including siRNA adaptation, dsRNA enhancement with dsRIP, and placeholder pages for BLAST, RNAi susceptibility prediction, and orthology inference. Availability in that branch does not mean a tool has been fully validated.

Clone the default branch using the instructions below. To work on the development snapshot, run `git switch _dev` after cloning.

## Included tools

| Tool | Purpose |
| --- | --- |
| RNA-seq Preprocessing | Batch adapter and quality trimming with Trim Galore or Cutadapt |
| dsRNA Identification | Rank candidate loci by bidirectional RNA-seq coverage |
| dsRNA Plotter | Plot directional RNA-seq coverage across grouped biological replicates |
| sRNA-based dsRNA Identification | Rank reference regions by bidirectional small RNA coverage |
| Small RNA Mapping | Map small RNAs and summarize coverage, read lengths, and unique sequences |
| Small RNA Control Mapping & Filtering | Screen against control references and flag low-complexity sequences |
| Target Prediction | Predict small RNA target sites using RNAplex pairing and MFE-ratio filtering |
| Degradome Analysis | Evaluate predicted target sites using degradome cleavage evidence |
| FASTA Deduplication | Cluster nucleotide FASTA records with CD-HIT-EST |

The degradome implementation is INCI's own CleaveLand-style workflow. The CleaveLand4/GSTAr distribution and dsRIP are not required or bundled for these tools.

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

Create or select a project, provide your own input files and references, and check the external tools in **Settings** before starting an analysis. File and folder dialogs use macOS integration. On other operating systems, type or paste file and folder paths into the interface.

Python requirements cover the core pipeline and its supporting scripts. They are not a locked, fully reproduced analysis environment. External executables must be installed separately.

## External executables

Install the tools needed by the workflows you use and make their executables available on `PATH`. The interface also checks several local installation locations.

| Tool | Used for |
| --- | --- |
| Trim Galore or Cutadapt | Sequencing-read preprocessing |
| Bowtie 1 (`bowtie`, `bowtie-build`) | Small RNA mapping, control filtering, and target/degradome workflows |
| minimap2, SAMtools, and SeqKit | RNA-seq mapping, genomic reference binning, and alignment processing |
| ViennaRNA (`RNAplex`) | Target interaction energies and pairing |
| CD-HIT (`cd-hit-est`) | Nucleotide FASTA deduplication |

Bowtie 2 does not replace Bowtie 1 for these workflows. The Python dependency file does not install the programs in this table.

## Inputs and outputs

Supply your own FASTA references and sequencing files. Supported inputs depend on the module and include FASTQ, gzipped FASTQ, FASTA, SAM/BAM, and CSV/TSV sample manifests. Consult the selected tool or its `--help` output for the exact format and parameters.

For example, the dsRNA plotter accepts a paired-end sample CSV with these columns:

```csv
sample_id,group,replicate,read1,read2
sample_1,condition_1,1,/path/to/sample_1_R1.fastq.gz,/path/to/sample_1_R2.fastq.gz
```

Project outputs are saved under `outputs/<project>/` by default or in the project directory you select. Local paths and interface preferences are stored in `.inci_pipeline_paths.json`.

Raw sequencing data, reference collections, generated results, alignment indexes, local environments, and machine-specific settings are excluded from this repository. Supply your own genome or control small RNA FASTA for Small RNA Control Mapping & Filtering; the development machine's default control reference is not included.

## Command-line use

Several analyses can run independently of the browser interface. Explore their options from the repository root:

```sh
python dsRNA_identification.py --help
python dsRNA_plotter.py --help
python sRNA_identification.py --help
python MFE_ratio.py --help
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
