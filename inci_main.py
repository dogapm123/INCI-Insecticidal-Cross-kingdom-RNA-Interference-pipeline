"""Offline browser interface for the INCI RNA interference pipeline.

Run this file from PyCharm or a terminal. It starts a local-only Python server,
opens a dashboard in your browser, and records paths to local input files.
"""

from __future__ import annotations

import html
import csv
import datetime
import json
import math
import os
import re
import shutil
import signal
import socket
import subprocess
import sys
import threading
import time
import traceback
import urllib.parse
import uuid
import webbrowser
from dataclasses import dataclass, field
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any

from degradome_analysis import DegradomeConfig, DegradomeSample, normalize_pvalue_method, run_degradome_analysis
from sirna_annotations import draw_sirna_annotation_track, find_sirna_annotations, write_sirna_annotation_table
from target_prediction import TargetPredictionConfig, run_target_prediction


APP_DIR = Path(__file__).resolve().parent
STATE_FILE = APP_DIR / ".inci_pipeline_paths.json"
OUTPUT_DIR = APP_DIR / "outputs"
DEFAULT_CONTROL_SRNA_FASTA = APP_DIR / "data" / "Tcastaneum_sRNA_control.fasta"
TRIM_GALORE_ENV_DIR = APP_DIR / "external_tools" / "trim_galore_env"
GENERAL_TOOLS_ENV_DIR = APP_DIR / "external_tools" / "inci_tools_env"
HOST = "127.0.0.1"
FAVICON_SVG = """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64">
<rect width="64" height="64" rx="13" fill="#0f5f5a"/>
<path d="M17 15c22 0 8 34 30 34M47 15c-22 0-8 34-30 34" fill="none" stroke="#d7f3ee" stroke-width="6" stroke-linecap="round"/>
<circle cx="22" cy="22" r="3" fill="#f4b740"/><circle cx="42" cy="42" r="3" fill="#f4b740"/>
</svg>"""
os.environ.setdefault("MPLCONFIGDIR", str(Path(os.environ.get("TMPDIR", "/tmp")) / "matplotlib-inci"))
RUNNING_PROCESSES: dict[int, tuple[subprocess.Popen[str], str]] = {}
RUNNING_PROCESSES_LOCK = threading.Lock()
PIPELINE_JOBS: dict[str, dict[str, Any]] = {}
PIPELINE_JOBS_LOCK = threading.Lock()


@dataclass(frozen=True)
class InputSpec:
    key: str
    label: str
    example: str
    extensions: str


@dataclass(frozen=True)
class OutputSpec:
    label: str
    filename: str


@dataclass(frozen=True)
class ModuleSpec:
    key: str
    title: str
    subtitle: str
    inputs: tuple[str, ...]
    outputs: tuple[OutputSpec, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class ModuleGroup:
    title: str
    keys: tuple[str, ...]


@dataclass(frozen=True)
class ExternalToolSpec:
    key: str
    title: str
    executables: tuple[str, ...]
    conda_package: str
    brew_package: str
    used_by: str
    description: str


INPUTS: tuple[InputSpec, ...] = (
    InputSpec("rnaseq", "RNA-seq", "RNAseq.fastq / RNAseq.fastq.gz / sample_manifest.csv", ".fastq, .fq, .fastq.gz, .fq.gz, .bam, .sam, .csv, .tsv"),
    InputSpec("transcriptome", "Transcriptome", "Transcriptome.fasta", ".fasta, .fa, .fas, .fna"),
    InputSpec("srna_fasta", "sRNA", "sRNA.fasta", ".fasta, .fa, .fas, .fna"),
    InputSpec("degradome", "Degradome-seq", "Degradome.fastq / Degradome.fasta", ".fastq, .fq, .fasta, .fa"),
    InputSpec("srnaseq", "sRNA-seq", "sRNAseq.fastq / sRNAseq.fastq.gz", ".fastq, .fq, .fastq.gz, .fq.gz, .bam, .sam"),
)


MODULES: tuple[ModuleSpec, ...] = (
    ModuleSpec(
        "preprocessing",
        "RNA-seq Preprocessing",
        "Trim raw sequencing files in reusable batches.",
        ("rnaseq",),
        (
            OutputSpec("Preprocessed RNA-seq", "rnaseq_preprocessed.fastq"),
            OutputSpec("Preprocessing report", "preprocessing_report.html"),
        ),
    ),
    ModuleSpec(
        "dsrna-identification",
        "dsRNA Identification",
        "Screen genomic or multi-contig references for locally balanced, bidirectionally covered dsRNA loci.",
        ("rnaseq", "transcriptome"),
        (
            OutputSpec("All scored bins", "tables/dsrna_loci_all_scored_bins.csv"),
            OutputSpec("Top-hit table", "tables/dsrna_loci_top_hits.csv"),
            OutputSpec("Summary plots", "phase1/summary_plots"),
            OutputSpec("Top-hit coverage plots", "top_hits/plots"),
            OutputSpec("Run manifest", "run_manifest.json"),
        ),
    ),
    ModuleSpec(
        "dsrna-plotter",
        "dsRNA Plotter",
        "Map grouped paired-end RNA-seq biological replicates to template sequences and plot bidirectional coverage per contig.",
        ("rnaseq", "transcriptome"),
        (
            OutputSpec("Directional coverage plots", "plots"),
            OutputSpec("Per-base global CPM tables", "tables"),
            OutputSpec("Sample mapping summary", "tables/sample_summary.csv"),
            OutputSpec("Contig summary", "tables/contig_summary.csv"),
            OutputSpec("Run manifest", "run_manifest.json"),
        ),
    ),
    ModuleSpec(
        "srna-dsrna-identification",
        "sRNA-based dsRNA Identification",
        "Rank reference bins by high bidirectional sRNA coverage and plot top contexts.",
        ("srnaseq", "transcriptome"),
        (
            OutputSpec("All scored sRNA bins", "tables/srna_dsrna_all_scored_bins.tsv"),
            OutputSpec("Top-hit table", "tables/srna_dsrna_top_hits.tsv"),
            OutputSpec("Top-hit context plots", "top_hits/plots"),
            OutputSpec("Length distributions", "length_distributions"),
            OutputSpec("Top-hit context FASTA", "top_hits/srna_dsrna_top_hit_contexts.fasta"),
            OutputSpec("Sample filter summary TSV", "tables/sample_summary.tsv"),
            OutputSpec("Run manifest", "run_manifest.json"),
        ),
    ),
    ModuleSpec(
        "srna-mapping",
        "Small RNA Mapping",
        "Map grouped sRNA samples to pasted or FASTA references with Bowtie1.",
        ("srnaseq",),
        (
            OutputSpec("Coverage plots", "plots"),
            OutputSpec("Length distributions", "length_distributions"),
            OutputSpec("Mapping summary TSV", "tables/contig_summary.tsv"),
            OutputSpec("Sample filter summary TSV", "tables/sample_summary.tsv"),
            OutputSpec("Unique siRNA table TSV", "tables/unique_sRNAs.tsv"),
            OutputSpec("Top exported unique siRNA table TSV", "tables/top_unique_sRNAs.tsv"),
            OutputSpec("Per-contig length table TSV", "tables/mapped_length_distribution_by_contig.tsv"),
            OutputSpec("Filtered mapped sRNA FASTA", "filtered_mapped_sRNAs.fasta"),
        ),
    ),
    ModuleSpec(
        "srna-control-filtering",
        "Small RNA Control Mapping & Filtering",
        "Screen mapped unique sRNAs against control references and flag low-complexity reads.",
        ("srna_fasta",),
        (
            OutputSpec("Control mapping summary TSV", "tables/control_mapping_summary.tsv"),
            OutputSpec("Complexity summary TSV", "tables/complexity_summary.tsv"),
            OutputSpec("Clean unique sRNA FASTA", "filtered_unique_sRNAs.clean.fasta"),
            OutputSpec("Short-contained collapse summary", "tables/short_contained_collapse_summary.tsv"),
            OutputSpec("Short-contained removed reads", "tables/short_contained_removed.tsv"),
        ),
    ),
    ModuleSpec(
        "degradome-analysis",
        "Degradome Analysis",
        "Analyze degradome evidence for cleavage-supported targets.",
        ("degradome", "transcriptome", "srna_fasta"),
        (
            OutputSpec("Cleavage candidates", "degradome_cleavage_sites.tsv"),
            OutputSpec("Target predictions", "tables/target_predictions.tsv"),
            OutputSpec("Transcript summary", "tables/transcript_summary.tsv"),
            OutputSpec("sRNA summary", "tables/srna_summary.tsv"),
            OutputSpec("Overall summary", "tables/degradome_overall_summary.tsv"),
            OutputSpec("Degradome density", "tables/degradome_5p_density_by_sample.tsv"),
            OutputSpec("PNG plots", "plots"),
            OutputSpec("Run report", "degradome_analysis_report.html"),
        ),
    ),
    ModuleSpec(
        "target-prediction",
        "Target Prediction",
        "Predict CleaveLand/GSTAr-style sRNA-transcript pairs without degradome evidence.",
        ("transcriptome", "srna_fasta"),
        (
            OutputSpec("Retained sRNA-transcript pairs", "target_prediction_pairs.tsv"),
            OutputSpec("Transcript summary", "tables/transcript_target_summary.tsv"),
            OutputSpec("sRNA summary", "tables/srna_target_summary.tsv"),
            OutputSpec("Summary plot", "plots/target_prediction_summary.png"),
            OutputSpec("Transcript prediction plots", "plots"),
            OutputSpec("Run report", "target_prediction_report.html"),
        ),
    ),
    ModuleSpec(
        "fasta-deduplication",
        "FASTA Deduplication",
        "Cluster nucleotide FASTA records with CD-HIT-EST and keep nonredundant representatives.",
        ("transcriptome",),
        (
            OutputSpec("Deduplicated FASTA", "deduplicated.fasta"),
            OutputSpec("Cluster summary", "clusters.csv"),
            OutputSpec("CD-HIT cluster file", "deduplicated.fasta.clstr"),
            OutputSpec("Run manifest", "run_manifest.json"),
        ),
    ),
)


MODULE_GROUPS: tuple[ModuleGroup, ...] = (
    ModuleGroup("RNA-seq Preprocessing", ("preprocessing",)),
    ModuleGroup("dsRNA Analysis", ("dsrna-identification", "dsrna-plotter")),
    ModuleGroup("Small RNA Analysis", ("srna-mapping", "srna-dsrna-identification", "srna-control-filtering")),
    ModuleGroup("Target Analysis", ("target-prediction", "degradome-analysis")),
    ModuleGroup("Supporting Tools", ("fasta-deduplication",)),
)


TOOL_PURPOSES: dict[str, str] = {
    "preprocessing": "Trim raw sequencing reads into reusable analysis batches.",
    "dsrna-identification": "Find bidirectionally covered dsRNA loci in a reference.",
    "dsrna-plotter": "Plot directional RNA-seq coverage across dsRNA templates.",
    "srna-mapping": "Map small RNAs and compare coverage across groups.",
    "srna-dsrna-identification": "Find dsRNA loci supported by bidirectional small RNAs.",
    "srna-control-filtering": "Remove control-mapping and low-complexity small RNAs.",
    "target-prediction": "Predict small RNA target sites in transcripts.",
    "degradome-analysis": "Test predicted targets for degradome cleavage support.",
    "fasta-deduplication": "Collapse redundant FASTA records into representatives.",
}


SEQUENCING_DATASETS: tuple[tuple[str, str, str], ...] = (
    ("rnaseq", "dsRNA-seq data", "Paired- or single-end dsRNA/RNA-seq reads or a sample manifest"),
    ("srnaseq", "sRNA-seq data", "Small RNA reads or a sample manifest"),
    ("degradome", "Degradome-seq data", "Degradome/PARE reads or a sample manifest"),
)


MODULE_DATASETS: dict[str, tuple[str, ...]] = {
    "dsrna-identification": ("rnaseq",),
    "dsrna-plotter": ("rnaseq",),
    "srna-mapping": ("srnaseq", "rnaseq"),
    "srna-dsrna-identification": ("srnaseq",),
    "degradome-analysis": ("degradome",),
}


EXTERNAL_TOOLS: tuple[ExternalToolSpec, ...] = (
    ExternalToolSpec("trim-galore", "Trim Galore", ("trim_galore",), "trim-galore", "trim-galore", "RNA-seq Preprocessing", "Automated adapter and quality trimming for sequencing reads."),
    ExternalToolSpec("cutadapt", "Cutadapt", ("cutadapt",), "cutadapt", "cutadapt", "RNA-seq Preprocessing", "Configurable adapter and quality trimming for sequencing reads."),
    ExternalToolSpec("bowtie", "Bowtie 1", ("bowtie", "bowtie-build"), "bowtie", "bowtie", "Small RNA and Target Analysis", "Short-read alignment and reference indexing."),
    ExternalToolSpec("samtools", "SAMtools", ("samtools",), "samtools", "samtools", "dsRNA Analysis", "SAM/BAM conversion, sorting, indexing, and filtering."),
    ExternalToolSpec("minimap2", "minimap2", ("minimap2",), "minimap2", "minimap2", "dsRNA Analysis", "Fast sequence mapping for transcript and genome-scale references."),
    ExternalToolSpec("seqkit", "SeqKit", ("seqkit",), "seqkit", "seqkit", "dsRNA Identification", "Read subsampling for the initial genomic dsRNA screen."),
    ExternalToolSpec("viennarna", "ViennaRNA", ("RNAplex",), "viennarna", "viennarna", "Target Analysis", "RNA duplex energy and structural scoring with RNAplex."),
    ExternalToolSpec("cd-hit", "CD-HIT", ("cd-hit-est",), "cd-hit", "cd-hit", "FASTA Deduplication", "Nucleotide FASTA clustering and nonredundant representative selection with CD-HIT-EST."),
)


DEFAULT_PLOT_SETTINGS: dict[str, Any] = {
    "font_family": "DejaVu Sans",
    "font_size": 11.0,
    "title_size": 14.0,
    "line_width": 2.0,
    "grid_width": 0.7,
    "marker_size": 6.0,
    "dpi": 220,
    "figure_width": 12.0,
    "figure_height": 5.0,
}


def default_state() -> dict[str, Any]:
    return {
        "project": {
            "name": "",
            "dir": "",
        },
        "paths": {},
        "preprocessed_paths": {},
        "preprocessing_status": {},
        "raw_overrides": {},
        "input_sources": {},
        "module_paths": {},
        "module_samples": {},
        "trimmed_batches": {},
        "plot_settings": dict(DEFAULT_PLOT_SETTINGS),
        "trimming": {
            "batch_name": "",
            "data_type": "rnaseq",
            "tool": "trim_galore",
            "mode": "paired",
            "quality": 20,
            "length": 18,
            "adapter1": "",
            "adapter2": "",
            "cutadapt_error_rate": 0.1,
            "cutadapt_overlap": 3,
            "cutadapt_max_n": "",
            "cutadapt_trim_n": True,
            "cutadapt_cores": 1,
            "cutadapt_pair_filter": "any",
            "samples": [],
            "manifest": "",
            "last_run": {},
        },
    }


def normalize_state(data: Any) -> dict[str, Any]:
    state = default_state()
    if not isinstance(data, dict):
        return state

    legacy_paths = {key: value for key, value in data.items() if isinstance(value, str)}
    if legacy_paths and "paths" not in data:
        state["paths"] = legacy_paths
        return state

    for section in state:
        value = data.get(section, {})
        if isinstance(value, dict):
            if section == "trimming":
                merged = dict(state["trimming"])
                merged.update(value)
                if not isinstance(merged.get("samples"), list):
                    merged["samples"] = []
                state[section] = merged
            else:
                state[section] = dict(value)
    return state


def load_state() -> dict[str, Any]:
    if not STATE_FILE.exists():
        return default_state()
    try:
        data = json.loads(STATE_FILE.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return default_state()
    return normalize_state(data)


def save_state(state: dict[str, Any]) -> None:
    STATE_FILE.write_text(json.dumps(normalize_state(state), indent=2, sort_keys=True), encoding="utf-8")


def normalized_plot_settings(value: Any) -> dict[str, Any]:
    incoming = value if isinstance(value, dict) else {}
    settings = dict(DEFAULT_PLOT_SETTINGS)
    family = str(incoming.get("font_family", settings["font_family"])).strip()
    settings["font_family"] = family or DEFAULT_PLOT_SETTINGS["font_family"]
    limits = {
        "font_size": (6.0, 40.0),
        "title_size": (6.0, 48.0),
        "line_width": (0.25, 10.0),
        "grid_width": (0.1, 5.0),
        "marker_size": (1.0, 30.0),
        "figure_width": (4.0, 30.0),
        "figure_height": (3.0, 20.0),
    }
    for key, (minimum, maximum) in limits.items():
        try:
            number = float(incoming.get(key, settings[key]))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{key.replace('_', ' ').title()} must be a number.") from exc
        if not minimum <= number <= maximum:
            raise ValueError(f"{key.replace('_', ' ').title()} must be between {minimum:g} and {maximum:g}.")
        settings[key] = number
    try:
        dpi = int(incoming.get("dpi", settings["dpi"]))
    except (TypeError, ValueError) as exc:
        raise ValueError("DPI must be a whole number.") from exc
    if not 72 <= dpi <= 600:
        raise ValueError("DPI must be between 72 and 600.")
    settings["dpi"] = dpi
    return settings


def plot_settings(state: dict[str, Any] | None = None) -> dict[str, Any]:
    if state is None:
        state = load_state()
    return normalized_plot_settings(state.get("plot_settings", {}))


def apply_global_plot_settings(state: dict[str, Any] | None = None) -> dict[str, Any]:
    settings = plot_settings(state)
    os.environ["INCI_PLOT_SETTINGS"] = json.dumps(settings)
    try:
        import matplotlib

        matplotlib.rcParams.update(
            {
                "font.family": settings["font_family"],
                "font.size": settings["font_size"],
                "axes.titlesize": settings["title_size"],
                "axes.labelsize": settings["font_size"],
                "lines.linewidth": settings["line_width"],
                "lines.markersize": settings["marker_size"],
                "grid.linewidth": settings["grid_width"],
                "legend.fontsize": max(6.0, settings["font_size"] - 1.0),
                "figure.dpi": settings["dpi"],
                "savefig.dpi": settings["dpi"],
            }
        )
    except ImportError:
        pass
    return settings


def project_name(state: dict[str, Any]) -> str:
    project = state.get("project", {})
    name = str(project.get("name", "") if isinstance(project, dict) else "").strip()
    return name or "Untitled project"


def project_slug(state: dict[str, Any] | None = None) -> str:
    if state is None:
        state = load_state()
    return slugify_sample_token(project_name(state), "untitled_project").lower()


def project_custom_dir(state: dict[str, Any]) -> str:
    project = state.get("project", {})
    directory = str(project.get("dir", "") if isinstance(project, dict) else "").strip()
    return directory


def project_output_root(state: dict[str, Any] | None = None) -> Path:
    if state is None:
        state = load_state()
    custom_dir = project_custom_dir(state)
    if custom_dir:
        return Path(custom_dir).expanduser()
    return OUTPUT_DIR / project_slug(state)


def raw_paths(state: dict[str, Any]) -> dict[str, str]:
    return {key: value for key, value in state.get("paths", {}).items() if isinstance(value, str)}


def preprocessed_paths(state: dict[str, Any]) -> dict[str, str]:
    return {key: value for key, value in state.get("preprocessed_paths", {}).items() if isinstance(value, str)}


def module_uses_raw(state: dict[str, Any], module_key: str, input_key: str | None = None) -> bool:
    if input_key:
        module_sources = state.get("input_sources", {}).get(module_key, {})
        if isinstance(module_sources, dict) and input_key in module_sources:
            return module_sources[input_key] == "loaded"
    return bool(state.get("raw_overrides", {}).get(module_key, False))


def effective_path(state: dict[str, Any], key: str, module_key: str | None = None) -> tuple[str, str]:
    if module_key:
        module_paths = state.get("module_paths", {}).get(module_key, {})
        if isinstance(module_paths, dict):
            selected = str(module_paths.get(key, "")).strip()
            if selected:
                return selected, "Selected for this tool"
        if key in MODULE_DATASETS.get(module_key, ()):
            return "", "Not selected"
    paths = raw_paths(state)
    preprocessed = preprocessed_paths(state)
    if module_key and module_uses_raw(state, module_key, key):
        return paths.get(key, ""), "Raw selected"
    if preprocessed.get(key):
        return preprocessed[key], "Preprocessed"
    if paths.get(key):
        return paths[key], "Raw"
    return "", "Missing"


def valid_module_dataset(module_key: str, input_key: str) -> bool:
    return module_for_key(module_key) is not None and input_key in MODULE_DATASETS.get(module_key, ())


def normalized_tool_samples(samples: Any) -> list[dict[str, Any]]:
    if not isinstance(samples, list):
        return []
    normalized: list[dict[str, Any]] = []
    group_counts: dict[str, int] = {}
    for index, value in enumerate(samples, start=1):
        if not isinstance(value, dict):
            continue
        read1 = str(value.get("trimmed_read1", value.get("read1", ""))).strip()
        if not read1:
            continue
        sample_id = slugify_sample_token(str(value.get("sample_id", "")).strip() or Path(read1).name.split(".")[0], f"sample_{index}")
        group = str(value.get("group", "")).strip() or sample_id
        included_value = value.get("included", True)
        included = included_value if isinstance(included_value, bool) else str(included_value).strip().lower() not in {"0", "false", "no", "off"}
        replicate = ""
        if included:
            group_key = group.casefold()
            group_counts[group_key] = group_counts.get(group_key, 0) + 1
            replicate = str(group_counts[group_key])
        read2 = str(value.get("trimmed_read2", value.get("read2", ""))).strip()
        normalized.append(
            {
                "id": str(value.get("id", "")).strip() or uuid.uuid4().hex,
                "sample_id": sample_id,
                "group": group,
                "included": included,
                "replicate": replicate,
                "mode": "paired" if read2 else "single",
                "trimmed_read1": read1,
                "trimmed_read2": read2,
                "status": str(value.get("status", "user_processed")) or "user_processed",
                "source": str(value.get("source", "local")) or "local",
            }
        )
    return normalized


def stored_tool_samples(state: dict[str, Any], module_key: str, input_key: str) -> list[dict[str, Any]]:
    module_samples = state.get("module_samples", {}).get(module_key, {})
    if not isinstance(module_samples, dict):
        return []
    return normalized_tool_samples(module_samples.get(input_key, []))


def set_tool_samples(state: dict[str, Any], module_key: str, input_key: str, samples: Any) -> list[dict[str, Any]]:
    normalized = normalized_tool_samples(samples)
    state.setdefault("module_samples", {}).setdefault(module_key, {})[input_key] = normalized
    return normalized


def write_tool_sample_manifest(state: dict[str, Any], module_key: str, input_key: str) -> str:
    samples = [sample for sample in stored_tool_samples(state, module_key, input_key) if sample.get("included", True)]
    module_paths = state.setdefault("module_paths", {}).setdefault(module_key, {})
    if not samples:
        module_paths.pop(input_key, None)
        return ""
    outdir = project_output_root(state) / "selections" / module_key
    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = outdir / f"{input_key}_samples.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample_id", "group", "replicate", "mode", "trimmed_read1", "trimmed_read2", "status"],
        )
        writer.writeheader()
        for sample in samples:
            writer.writerow({key: sample[key] for key in writer.fieldnames})
    module_paths[input_key] = str(manifest_path)
    return str(manifest_path)


def rows_from_sample_manifest(path: str) -> list[dict[str, str]]:
    dataset_path = Path(path).expanduser()
    if dataset_path.suffix.lower() not in {".csv", ".tsv"} or not dataset_path.exists():
        return []
    delimiter = "\t" if dataset_path.suffix.lower() == ".tsv" else ","
    with dataset_path.open(newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter=delimiter)]


def append_tool_samples(
    state: dict[str, Any], module_key: str, input_key: str, additions: list[dict[str, Any]]
) -> tuple[list[dict[str, Any]], int]:
    samples = stored_tool_samples(state, module_key, input_key)
    known_paths = {(sample["trimmed_read1"], sample.get("trimmed_read2", "")) for sample in samples}
    known_sample_ids = {sample["sample_id"] for sample in samples}
    added = 0
    for value in additions:
        read1 = str(value.get("trimmed_read1", value.get("read1", ""))).strip()
        read2 = str(value.get("trimmed_read2", value.get("read2", ""))).strip()
        if not read1 or (read1, read2) in known_paths:
            continue
        base_sample_id = slugify_sample_token(
            str(value.get("sample_id", "")).strip() or Path(read1).name.split(".")[0], "sample"
        )
        sample_id = base_sample_id
        suffix = 2
        while sample_id in known_sample_ids:
            sample_id = f"{base_sample_id}_{suffix}"
            suffix += 1
        samples.append(
            {
                "id": uuid.uuid4().hex,
                "sample_id": sample_id,
                "group": str(value.get("group", "")).strip(),
                "trimmed_read1": read1,
                "trimmed_read2": read2,
                "status": str(value.get("status", "user_processed")),
                "source": str(value.get("source", "local")),
                "included": True,
            }
        )
        known_paths.add((read1, read2))
        known_sample_ids.add(sample_id)
        added += 1
    return set_tool_samples(state, module_key, input_key, samples), added


def add_batch_to_tool_samples(state: dict[str, Any], module_key: str, input_key: str, batch_id: str) -> int:
    if not valid_module_dataset(module_key, input_key):
        raise ValueError("Unknown tool input.")
    batch = state.get("trimmed_batches", {}).get(batch_id, {})
    if not isinstance(batch, dict) or batch.get("data_type") != input_key:
        raise ValueError("That preprocessing batch is not available for this input.")
    rows = rows_from_sample_manifest(str(batch.get("manifest", "")))
    if not rows:
        raise ValueError("The preprocessing batch has no readable samples.")
    additions = [{**row, "source": f"batch:{batch_id}"} for row in rows]
    _, added = append_tool_samples(state, module_key, input_key, additions)
    write_tool_sample_manifest(state, module_key, input_key)
    return added


def add_local_tool_sample(
    state: dict[str, Any], module_key: str, input_key: str, read1: str, read2: str
) -> int:
    if not valid_module_dataset(module_key, input_key):
        raise ValueError("Unknown tool input.")
    first = Path(read1).expanduser()
    second = Path(read2).expanduser() if read2 else None
    if not read1 or not first.is_file():
        raise ValueError("Choose an existing preprocessed sequence file.")
    if second is not None and not second.is_file():
        raise ValueError("Read 2 does not exist.")
    if input_key != "rnaseq" and second is not None:
        raise ValueError("Read 2 is only supported for RNA-seq samples.")
    _, added = append_tool_samples(
        state,
        module_key,
        input_key,
        [
            {
                "sample_id": first.name.split(".")[0],
                "group": "",
                "trimmed_read1": str(first),
                "trimmed_read2": str(second) if second is not None else "",
                "status": "user_processed",
                "source": "local",
            }
        ],
    )
    if not added:
        raise ValueError("That sample is already in the selected sample set.")
    write_tool_sample_manifest(state, module_key, input_key)
    return added


def update_tool_sample_settings(
    state: dict[str, Any], module_key: str, input_key: str, settings: dict[str, Any]
) -> list[dict[str, Any]]:
    if not valid_module_dataset(module_key, input_key):
        raise ValueError("Unknown tool input.")
    samples = stored_tool_samples(state, module_key, input_key)
    for sample in samples:
        value = settings.get(sample["id"])
        if isinstance(value, dict):
            sample["group"] = str(value.get("group", "")).strip() or sample["sample_id"]
            sample["included"] = bool(value.get("included", False))
        elif value is not None:
            sample["group"] = str(value).strip() or sample["sample_id"]
    normalized = set_tool_samples(state, module_key, input_key, samples)
    write_tool_sample_manifest(state, module_key, input_key)
    return normalized


def remove_tool_sample(state: dict[str, Any], module_key: str, input_key: str, sample_id: str) -> bool:
    if not valid_module_dataset(module_key, input_key):
        raise ValueError("Unknown tool input.")
    samples = stored_tool_samples(state, module_key, input_key)
    remaining = [sample for sample in samples if sample["id"] != sample_id]
    if len(remaining) == len(samples):
        return False
    set_tool_samples(state, module_key, input_key, remaining)
    write_tool_sample_manifest(state, module_key, input_key)
    return True


def preprocessing_text(state: dict[str, Any], key: str) -> str:
    preprocessed = preprocessed_paths(state).get(key, "")
    if preprocessed:
        return f"Preprocessed > {Path(preprocessed).name}"
    if raw_paths(state).get(key):
        return "Preprocessing pending"
    return "No raw path yet"


def planned_preprocessed_path(key: str) -> Path:
    filenames = {
        "rnaseq": "rnaseq_preprocessed.fastq",
        "transcriptome": "transcriptome_preprocessed.fasta",
        "srna_fasta": "srna_preprocessed.fasta",
        "degradome": "degradome_preprocessed.fastq",
        "srnaseq": "srnaseq_preprocessed.fastq",
    }
    return project_output_root() / "preprocessing" / filenames[key]


def trimming_output_dir(batch_id: str = "latest") -> Path:
    return project_output_root() / "trimming" / slugify_sample_token(batch_id, "latest")


def trimming_manifest_path(batch_id: str = "latest") -> Path:
    return trimming_output_dir(batch_id) / "trimmed_samples_manifest.csv"


def slugify_sample_token(value: str, fallback: str = "group") -> str:
    token = re.sub(r"[^A-Za-z0-9]+", "_", value.strip()).strip("_")
    return token or fallback


def trimmed_fastq_candidates(sample_id: str, paired: bool, mate: int = 1, outdir: Path | None = None) -> list[Path]:
    outdir = outdir or trimming_output_dir()
    if paired:
        suffixes = [f"_val_{mate}.fq.gz", f"_val_{mate}.fastq.gz", f"_val_{mate}.fq", f"_val_{mate}.fastq"]
    else:
        suffixes = ["_trimmed.fq.gz", "_trimmed.fastq.gz", "_trimmed.fq", "_trimmed.fastq"]
    return [outdir / f"{sample_id}{suffix}" for suffix in suffixes]


def first_existing_path(paths: list[Path]) -> Path:
    return next((path for path in paths if path.exists()), paths[0])


def cutadapt_fastq_path(sample_id: str, paired: bool, mate: int, outdir: Path) -> Path:
    if paired:
        return outdir / f"{sample_id}_R{mate}.trimmed.fq.gz"
    return outdir / f"{sample_id}.trimmed.fq.gz"


def trimming_tool_title(tool: str) -> str:
    return "Cutadapt" if tool == "cutadapt" else "Trim Galore"


def now_text() -> str:
    return datetime.datetime.now().isoformat(timespec="seconds")


def create_pipeline_job(label: str) -> str:
    job_id = uuid.uuid4().hex
    with PIPELINE_JOBS_LOCK:
        PIPELINE_JOBS[job_id] = {
            "id": job_id,
            "label": label,
            "status": "running",
            "started": now_text(),
            "finished": "",
            "current_message": f"Started {label}.",
            "logs": [{"time": now_text(), "level": "info", "message": f"Started {label}."}],
            "result": {},
        }
    return job_id


def append_job_log(job_id: str | None, message: str, level: str = "info") -> None:
    if not job_id:
        return
    with PIPELINE_JOBS_LOCK:
        job = PIPELINE_JOBS.get(job_id)
        if not job:
            return
        logs = job.setdefault("logs", [])
        if isinstance(logs, list):
            logs.append({"time": now_text(), "level": level, "message": message})
            del logs[:-120]
        job["current_message"] = message


def finish_pipeline_job(job_id: str, status: str, result: dict[str, Any] | None = None, message: str = "") -> None:
    with PIPELINE_JOBS_LOCK:
        job = PIPELINE_JOBS.get(job_id)
        if not job:
            return
        if job.get('cancel_requested'):
            status = 'stopped'
            message = 'Analysis stopped by the user.'
        job["status"] = status
        job["finished"] = now_text()
        if result is not None:
            job["result"] = result
        if message:
            job["current_message"] = message
            logs = job.setdefault("logs", [])
            if isinstance(logs, list):
                logs.append({"time": now_text(), "level": "error" if status == "failed" else "info", "message": message})


def pipeline_job_snapshot(job_id: str) -> dict[str, Any] | None:
    with PIPELINE_JOBS_LOCK:
        job = PIPELINE_JOBS.get(job_id)
        if not job:
            return None
        if job['status'] == 'running' and not any(
            thread.is_alive() and thread.name.endswith(job_id[:8])
            for thread in threading.enumerate()
        ):
            job.update(status='failed', finished=now_text(),
                       current_message='The analysis worker stopped unexpectedly.',
                       result={'message': 'The analysis worker stopped unexpectedly. Check the process console and run again.'})
        return json.loads(json.dumps(job))


def running_process_labels() -> list[str]:
    with RUNNING_PROCESSES_LOCK:
        return [label for process, label in RUNNING_PROCESSES.values() if process.poll() is None]


def run_tracked_command(command: list[str], label: str, job_id: str | None = None) -> subprocess.CompletedProcess[str]:
    if job_id is None:
        worker_name = threading.current_thread().name
        with PIPELINE_JOBS_LOCK:
            job_id = next((key for key in PIPELINE_JOBS if worker_name.endswith(key[:8])), None)
    with PIPELINE_JOBS_LOCK:
        if job_id and PIPELINE_JOBS.get(job_id, {}).get('cancel_requested'):
            raise RuntimeError('Analysis stopped by the user.')
    append_job_log(job_id, f"Running {label}.")
    popen_kwargs: dict[str, Any] = {
        "cwd": APP_DIR,
        "text": True,
        "stdout": subprocess.PIPE,
        "stderr": subprocess.PIPE,
    }
    if os.name != "nt":
        popen_kwargs["start_new_session"] = True
    process = subprocess.Popen(command, **popen_kwargs)
    with RUNNING_PROCESSES_LOCK:
        RUNNING_PROCESSES[process.pid] = (process, label)
    try:
        stdout, stderr = process.communicate()
        if process.returncode == 0:
            append_job_log(job_id, f"Finished {label}.")
        else:
            append_job_log(job_id, f"{label} exited with code {process.returncode}.", "error")
        return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)
    finally:
        with RUNNING_PROCESSES_LOCK:
            RUNNING_PROCESSES.pop(process.pid, None)


def stop_process_group(process: subprocess.Popen[str]) -> None:
    if process.poll() is not None:
        return
    try:
        if os.name == "nt":
            process.terminate()
        else:
            os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    except OSError:
        process.terminate()


def kill_process_group(process: subprocess.Popen[str]) -> None:
    if process.poll() is not None:
        return
    try:
        if os.name == "nt":
            process.kill()
        else:
            os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        return
    except OSError:
        process.kill()


def stop_running_processes() -> dict[str, Any]:
    with PIPELINE_JOBS_LOCK:
        jobs = [job for job in PIPELINE_JOBS.values() if job['status'] == 'running']
        for job in jobs:
            job['cancel_requested'] = True
            job['current_message'] = 'Stopping analysis...'
    with RUNNING_PROCESSES_LOCK:
        running = list(RUNNING_PROCESSES.items())
    if not running:
        return {"count": len(jobs), "labels": [job['label'] for job in jobs]}

    labels: list[str] = []
    for _, (process, label) in running:
        labels.append(label)
        stop_process_group(process)

    time.sleep(0.8)
    for _, (process, _) in running:
        kill_process_group(process)

    return {"count": len(running), "labels": labels}


RUN_PRESERVED_INPUT_DIRS = {"pasted_inputs"}
RUN_MARKER_NAME = ".inci_run.json"


def tool_output_root(module_key: str, state: dict[str, Any] | None = None) -> Path:
    folder = "trimming" if module_key == "preprocessing" else module_key
    return project_output_root(state) / folder


def prepare_tool_output_dir(module_key: str, state: dict[str, Any] | None = None) -> Path:
    tool_root = tool_output_root(module_key, state)
    tool_root.mkdir(parents=True, exist_ok=True)
    for path in list(tool_root.iterdir()):
        if path.name in RUN_PRESERVED_INPUT_DIRS:
            continue
        if path.is_dir() and (path / RUN_MARKER_NAME).exists():
            continue
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()
    return tool_root


def reset_output_paths_for_page(module_key: str, state: dict[str, Any]) -> dict[str, Any]:
    root = project_output_root(state)
    if module_key == "preprocessing":
        targets = [root / "preprocessing", root / "trimming"]
    elif module_for_key(module_key) is not None:
        targets = [root / module_key]
    else:
        raise ValueError("Unknown pipeline page.")

    removed: list[str] = []
    for target in targets:
        try:
            target.relative_to(root)
        except ValueError as exc:
            raise ValueError("Refusing to reset a path outside the active project output folder.") from exc
        if target.exists():
            if target.is_dir():
                shutil.rmtree(target)
            else:
                target.unlink()
            removed.append(str(target))

    if module_key == "preprocessing":
        manifest = str(state.get("trimming", {}).get("manifest", "") or "")
        state.setdefault("trimming", {})["manifest"] = ""
        state.setdefault("trimming", {})["last_run"] = {}
        state["trimmed_batches"] = {}
        for module_key_value, paths in list(state.get("module_paths", {}).items()):
            if not isinstance(paths, dict):
                continue
            for input_key, value in list(paths.items()):
                try:
                    Path(str(value)).relative_to(root / "trimming")
                except ValueError:
                    continue
                paths.pop(input_key, None)
        for key, value in list(state.get("preprocessed_paths", {}).items()):
            value_path = Path(str(value))
            if value == manifest:
                state["preprocessed_paths"].pop(key, None)
                state.get("preprocessing_status", {}).pop(key, None)
                continue
            try:
                value_path.relative_to(root)
            except ValueError:
                continue
            state["preprocessed_paths"].pop(key, None)
            state.get("preprocessing_status", {}).pop(key, None)
        save_state(state)

    return {"removed": removed, "root": str(root)}


def executable_candidates(name: str) -> list[Path]:
    home = Path.home()
    return [
        GENERAL_TOOLS_ENV_DIR / "bin" / name,
        TRIM_GALORE_ENV_DIR / "bin" / name,
        home / "miniconda3" / "bin" / name,
        home / "mambaforge" / "bin" / name,
        home / "micromamba" / "bin" / name,
        home / "anaconda3" / "bin" / name,
        Path("/opt/homebrew/bin") / name,
        Path("/usr/local/bin") / name,
        Path("/usr/bin") / name,
    ]


def resolve_executable(name: str) -> str:
    found = shutil.which(name)
    if found:
        return found
    for candidate in executable_candidates(name):
        if candidate.exists() and os.access(candidate, os.X_OK):
            return str(candidate)
    return ""


def resolve_trim_galore() -> str:
    return resolve_executable("trim_galore")


def resolve_cutadapt() -> str:
    return resolve_executable("cutadapt")


def external_tool_for_key(key: str) -> ExternalToolSpec | None:
    return next((tool for tool in EXTERNAL_TOOLS if tool.key == key), None)


def external_tool_status(tool: ExternalToolSpec) -> dict[str, Any]:
    paths = {name: resolve_executable(name) for name in tool.executables}
    return {"ready": all(paths.values()), "paths": paths}


def resolve_package_manager() -> tuple[str, str]:
    for name in ("mamba", "conda", "micromamba"):
        found = resolve_executable(name)
        if found:
            return name, found
    brew = resolve_executable("brew")
    if brew:
        return "brew", brew
    return "", ""


def install_trim_galore(job_id: str | None = None) -> dict[str, str]:
    existing = resolve_trim_galore()
    if existing:
        append_job_log(job_id, f"Trim Galore already available at {existing}.")
        return {"status": "ready", "path": existing, "message": "Trim Galore is already available."}

    manager_name, manager_path = resolve_package_manager()
    if not manager_path:
        raise ValueError(
            "Could not find conda, mamba, micromamba, or brew. Install one package manager or add trim_galore to PATH."
        )

    TRIM_GALORE_ENV_DIR.parent.mkdir(parents=True, exist_ok=True)
    append_job_log(job_id, f"Installing Trim Galore with {manager_name} into {TRIM_GALORE_ENV_DIR}.")
    if manager_name in {"conda", "mamba", "micromamba"}:
        command = [
            manager_path,
            "create",
            "-y",
            "-p",
            str(TRIM_GALORE_ENV_DIR),
            "-c",
            "bioconda",
            "-c",
            "conda-forge",
            "trim-galore",
            "cutadapt",
            "fastqc",
        ]
    else:
        command = [manager_path, "install", "trim-galore"]

    completed = run_tracked_command(command, f"Install Trim Galore with {manager_name}", job_id)
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout or "").strip()
        raise RuntimeError(f"Automatic Trim Galore installation failed: {detail}")

    installed = resolve_trim_galore()
    if not installed:
        env_executable = TRIM_GALORE_ENV_DIR / "bin" / "trim_galore"
        installed = str(env_executable) if env_executable.exists() else ""
    if not installed:
        raise RuntimeError("Installation finished, but trim_galore was still not found.")

    append_job_log(job_id, f"Trim Galore ready at {installed}.")
    return {
        "status": "installed",
        "path": installed,
        "message": f"Trim Galore installed with {manager_name}.",
    }


def install_cutadapt(job_id: str | None = None) -> dict[str, str]:
    existing = resolve_cutadapt()
    if existing:
        append_job_log(job_id, f"Cutadapt already available at {existing}.")
        return {"status": "ready", "path": existing, "message": "Cutadapt is already available."}

    manager_name, manager_path = resolve_package_manager()
    if not manager_path:
        raise ValueError("Could not find conda, mamba, micromamba, or brew. Install one package manager or add cutadapt to PATH.")

    append_job_log(job_id, f"Installing Cutadapt with {manager_name}.")
    if manager_name in {"conda", "mamba", "micromamba"}:
        GENERAL_TOOLS_ENV_DIR.parent.mkdir(parents=True, exist_ok=True)
        action = "install" if GENERAL_TOOLS_ENV_DIR.exists() else "create"
        command = [
            manager_path,
            action,
            "-y",
            "-p",
            str(GENERAL_TOOLS_ENV_DIR),
            "-c",
            "bioconda",
            "-c",
            "conda-forge",
            "cutadapt",
        ]
    else:
        command = [manager_path, "install", "cutadapt"]

    completed = run_tracked_command(command, f"Install Cutadapt with {manager_name}", job_id)
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout or "").strip()
        raise RuntimeError(f"Automatic Cutadapt installation failed: {detail}")

    installed = resolve_cutadapt()
    if not installed:
        env_executable = GENERAL_TOOLS_ENV_DIR / "bin" / "cutadapt"
        installed = str(env_executable) if env_executable.exists() else ""
    if not installed:
        raise RuntimeError("Installation finished, but cutadapt was still not found.")

    append_job_log(job_id, f"Cutadapt ready at {installed}.")
    return {
        "status": "installed",
        "path": installed,
        "message": f"Cutadapt installed with {manager_name}.",
    }


def install_external_tool(tool: ExternalToolSpec, job_id: str | None = None) -> dict[str, Any]:
    status = external_tool_status(tool)
    if status["ready"]:
        append_job_log(job_id, f"{tool.title} is already available.")
        return {"status": "ready", "tool": tool.key, "paths": status["paths"]}
    if tool.key == "trim-galore":
        result = install_trim_galore(job_id)
        return {"tool": tool.key, **result}
    if tool.key == "cutadapt":
        result = install_cutadapt(job_id)
        return {"tool": tool.key, **result}

    manager_name, manager_path = resolve_package_manager()
    if not manager_path:
        raise ValueError("Could not find conda, mamba, micromamba, or Homebrew for external-tool installation.")

    if manager_name in {"conda", "mamba", "micromamba"}:
        GENERAL_TOOLS_ENV_DIR.parent.mkdir(parents=True, exist_ok=True)
        action = "install" if GENERAL_TOOLS_ENV_DIR.exists() else "create"
        command = [
            manager_path,
            action,
            "-y",
            "-p",
            str(GENERAL_TOOLS_ENV_DIR),
            "-c",
            "bioconda",
            "-c",
            "conda-forge",
            tool.conda_package,
        ]
    else:
        command = [manager_path, "install", tool.brew_package]

    append_job_log(job_id, f"Installing {tool.title} with {manager_name}.")
    completed = run_tracked_command(command, f"Install {tool.title}", job_id)
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout or "").strip()
        raise RuntimeError(f"Could not install {tool.title}: {detail}")
    status = external_tool_status(tool)
    if not status["ready"]:
        raise RuntimeError(f"{tool.title} installation completed, but its executable was not found.")
    return {"status": "installed", "tool": tool.key, "paths": status["paths"]}


def start_external_tools_install_job(tool_keys: list[str]) -> str:
    selected: list[ExternalToolSpec] = []
    for key in tool_keys:
        tool = external_tool_for_key(key)
        if tool and tool not in selected:
            selected.append(tool)
    if not selected:
        raise ValueError("Choose at least one recognized external tool.")

    job_id = create_pipeline_job("external tool installation")

    def worker() -> None:
        results: list[dict[str, Any]] = []
        try:
            for tool in selected:
                results.append(install_external_tool(tool, job_id))
        except Exception as exc:
            finish_pipeline_job(job_id, "failed", {"message": str(exc), "tools": results}, f"External tool installation failed: {exc}")
            return
        finish_pipeline_job(job_id, "finished", {"tools": results}, f"Finished checking/installing {len(results)} external tool(s).")

    threading.Thread(target=worker, name=f"inci-install-tools-{job_id[:8]}", daemon=True).start()
    return job_id


def input_for_key(key: str) -> InputSpec:
    return next(spec for spec in INPUTS if spec.key == key)


def module_for_key(key: str) -> ModuleSpec | None:
    return next((module for module in MODULES if module.key == key), None)


def reveal_path(path: str | Path) -> None:
    target = Path(path).expanduser()
    if not target.exists():
        target = target.parent if target.parent.exists() else APP_DIR

    if sys.platform == "darwin":
        if target.is_file():
            subprocess.run(["open", "-R", str(target)], check=False)
        else:
            subprocess.run(["open", str(target)], check=False)
    elif os.name == "nt":
        if target.is_file():
            subprocess.run(["explorer", "/select,", str(target)], check=False)
        else:
            os.startfile(str(target))  # type: ignore[attr-defined]
    else:
        opener = "xdg-open"
        subprocess.run([opener, str(target if target.is_dir() else target.parent)], check=False)


def choose_local_file(title: str) -> str | None:
    if sys.platform == "darwin":
        script = f'POSIX path of (choose file with prompt "{title}")'
        result = subprocess.run(["osascript", "-e", script], text=True, capture_output=True, check=False)
        chosen = result.stdout.strip()
        return chosen or None
    return None


def choose_local_dir(title: str) -> str | None:
    if sys.platform == "darwin":
        script = f'POSIX path of (choose folder with prompt "{title}")'
        result = subprocess.run(["osascript", "-e", script], text=True, capture_output=True, check=False)
        chosen = result.stdout.strip()
        return chosen.rstrip("/") if chosen else None
    return None


def find_free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind((HOST, 0))
        return int(sock.getsockname()[1])


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def format_preview_value(value: object) -> str:
    text = str(value).strip()
    if not text:
        return ""
    try:
        numeric = float(text)
    except ValueError:
        return text
    if not math.isfinite(numeric):
        return text
    rounded = f"{numeric:.3f}".rstrip("0").rstrip(".")
    return rounded if rounded else "0"


def url_for(path: str | Path) -> str:
    return urllib.parse.quote(str(path), safe="")


def output_file_url(path: str | Path) -> str:
    target = Path(path)
    try:
        version = target.stat().st_mtime_ns
    except OSError:
        version = 0
    return f"/output-file?path={url_for(target)}&v={version}"


def render_file_cards(state: dict[str, Any]) -> str:
    cards: list[str] = []
    paths = raw_paths(state)
    for spec in INPUTS:
        path = paths.get(spec.key, "")
        status = "Path added" if path else "Waiting for file"
        status_class = "ready" if path else "missing"
        name = Path(path).name if path else "No path selected"
        prep = preprocessing_text(state, spec.key)
        reveal = (
            f'<a class="ghost" href="/reveal?path={url_for(path)}" onclick="writeTerminal(\'Revealing {esc(Path(path).name)} in Finder.\')">Reveal</a>'
            if path
            else '<span class="ghost disabled">Reveal</span>'
        )
        cards.append(
            f"""
            <section class="card input-card">
                <div>
                    <p class="eyebrow">{esc(spec.extensions)}</p>
                    <h3>{esc(spec.label)}</h3>
                    <p class="{status_class}">{esc(status)}</p>
                </div>
                <p class="filename">{esc(name)}</p>
                <p class="prep-status">{esc(prep)}</p>
                <input id="{esc(spec.key)}" value="{esc(path)}" placeholder="{esc(spec.example)}">
                <div class="button-row">
                    <button onclick="browseFile('{esc(spec.key)}')">Browse</button>
                    <button class="secondary" onclick="savePath('{esc(spec.key)}')">Save Path</button>
                    {reveal}
                </div>
            </section>
            """
        )
    return "\n".join(cards)


def render_summary(state: dict[str, Any]) -> str:
    rows: list[str] = []
    paths = raw_paths(state)
    for spec in INPUTS:
        path = paths.get(spec.key, "")
        label = f"{spec.label} > {Path(path).name}" if path else f"{spec.label} > {spec.example}"
        prep = preprocessing_text(state, spec.key)
        detail = esc(path) if path else "Add a local path when this input is available."
        status = "OK" if path else "--"
        action = (
            f'<a href="/reveal?path={url_for(path)}" onclick="writeTerminal(\'Revealing {esc(Path(path).name)} in Finder.\')">Open in Finder</a>'
            if path
            else ""
        )
        rows.append(
            f"""
            <div class="path-row">
                <span class="badge">{esc(status)}</span>
                <div>
                    <strong>{esc(label)}</strong>
                    <small>{detail}</small>
                    <small>{esc(prep)}</small>
                </div>
                {action}
            </div>
            """
        )
    return "\n".join(rows)


def sequencing_dataset_spec(key: str) -> tuple[str, str]:
    for dataset_key, label, description in SEQUENCING_DATASETS:
        if dataset_key == key:
            return label, description
    return key, "Shared sequencing data"


def render_dataset_selector(state: dict[str, Any], module_key: str, key: str) -> str:
    label, description = sequencing_dataset_spec(key)
    batches = state.get("trimmed_batches", {})
    batch_options: list[str] = []
    if isinstance(batches, dict):
        ordered = sorted(
            (
                batch
                for batch in batches.values()
                if isinstance(batch, dict)
                and batch.get("data_type") == key
                and str(batch.get("project_root", project_output_root(state))) == str(project_output_root(state))
            ),
            key=lambda batch: str(batch.get("finished", "")),
            reverse=True,
        )
        for batch in ordered:
            path = str(batch.get("manifest", ""))
            if not path or not Path(path).exists():
                continue
            batch_options.append(
                f"""
                <div class="dataset-option batch-option">
                    <span>
                        <strong>{esc(batch.get('name', 'Trimmed batch'))}</strong>
                        <small>{esc(batch.get('sample_count', 0))} sample(s) · {esc(batch.get('finished', ''))}</small>
                        <small>{esc(Path(path).name)}</small>
                    </span>
                    <button type="button" class="secondary" onclick="addToolBatch('{esc(module_key)}', '{esc(key)}', '{esc(batch.get('id', ''))}')">Add Batch</button>
                </div>
                """
            )
    batch_body = "".join(batch_options) or '<div class="empty-state compact-empty">No completed preprocessing batch for this data type.</div>'
    samples = stored_tool_samples(state, module_key, key)
    included_count = sum(1 for sample in samples if sample.get("included", True))
    sample_rows: list[str] = []
    for sample in samples:
        files = Path(sample["trimmed_read1"]).name
        if sample.get("trimmed_read2"):
            files += f" + {Path(sample['trimmed_read2']).name}"
        sample_rows.append(
            f"""
            <div class="tool-sample-row {'is-excluded' if not sample.get('included', True) else ''}" data-tool-sample data-sample-id="{esc(sample['id'])}" data-sample-name="{esc(sample['sample_id'])}" data-read1="{esc(sample['trimmed_read1'])}">
                <label class="sample-include"><input data-field="included" type="checkbox" {'checked' if sample.get('included', True) else ''} onchange="refreshToolReplicates('{esc(module_key)}', '{esc(key)}')"><span>Include</span></label>
                <div class="sample-file"><strong>{esc(sample['sample_id'])}</strong><small>{esc(files)}</small></div>
                <label class="field compact-field"><span>Group</span><input data-field="group" value="{esc(sample['group'])}" onchange="refreshToolReplicates('{esc(module_key)}', '{esc(key)}')"></label>
                <div class="replicate-value"><small>Biological replicate</small><strong data-replicate>{esc(sample['replicate'] or 'Excluded')}</strong></div>
                <button type="button" class="ghost icon-button" title="Remove sample" onclick="removeToolSample('{esc(module_key)}', '{esc(key)}', '{esc(sample['id'])}')">x</button>
            </div>
            """
        )
    selected_body = "".join(sample_rows) or '<div class="empty-state compact-empty">No samples added yet.</div>'
    read2 = ""
    if key == "rnaseq":
        read2 = f"""
            <div class="file-picker">
                <input id="tool-read2-{esc(module_key)}-{esc(key)}" placeholder="Read 2 (optional for single-end data)">
                <button type="button" class="secondary" onclick="browseToolSample('{esc(module_key)}', '{esc(key)}', 'read2')">Browse</button>
            </div>
        """
    return f"""
        <section class="panel dataset-selector">
            <div class="panel-heading">
                <div>
                    <p class="eyebrow">Dataset selection</p>
                    <h3>{esc(label)}</h3>
                    <p>{esc(description)}</p>
                </div>
            </div>
            <div class="dataset-source-grid">
                <div>
                    <h3>Completed preprocessing batches</h3>
                    <div class="dataset-options stacked">{batch_body}</div>
                </div>
                <div>
                    <h3>Already-preprocessed file</h3>
                    <p class="muted">Add one trimmed sample at a time, then assign its group in Selected Samples below.</p>
                    <div class="file-picker">
                        <input id="tool-read1-{esc(module_key)}-{esc(key)}" placeholder="{'Read 1' if key == 'rnaseq' else 'Trimmed sequence file'}">
                        <button type="button" class="secondary" onclick="browseToolSample('{esc(module_key)}', '{esc(key)}', 'read1')">Browse</button>
                    </div>
                    {read2}
                    <div class="button-row">
                        <button type="button" class="secondary" onclick="addToolSample('{esc(module_key)}', '{esc(key)}')">Add Sample</button>
                    </div>
                </div>
            </div>
            <div class="selected-samples">
                <div class="panel-heading sample-heading">
                    <div><p class="eyebrow">Analysis sample set</p><h3>{esc(label)} Samples</h3><p><span data-included-count>{included_count}</span> of {len(samples)} included</p></div>
                    <button type="button" class="secondary" {'disabled' if not samples else ''} onclick="saveToolSampleSettings('{esc(module_key)}', '{esc(key)}')">Save Sample Set</button>
                </div>
                <div id="tool-samples-{esc(module_key)}-{esc(key)}" class="tool-sample-list">{selected_body}</div>
            </div>
        </section>
    """


def render_module_dataset_selectors(state: dict[str, Any], module_key: str) -> str:
    return "".join(render_dataset_selector(state, module_key, key) for key in MODULE_DATASETS.get(module_key, ()))


def render_analysis_action(
    module_key: str,
    run_js: str,
    run_label: str,
    *,
    disabled: bool = False,
    extra_html: str = "",
) -> str:
    disabled_attr = " disabled" if disabled else ""
    return f"""
        <div class="analysis-action-group">
            <div class="analysis-action-bar">
                <div class="analysis-action-copy">
                    <p class="eyebrow">Ready to run</p>
                    <strong>Start Analysis</strong>
                    <small>Results appear below when the analysis completes.</small>
                </div>
                <div class="analysis-action-buttons">
                    {extra_html}
                    <a class="ghost" href="/reveal?path={url_for(tool_output_root(module_key))}">Output Folder</a>
                    <button type="button" class="clear-outputs" onclick="resetCurrentPage('{esc(module_key)}')">Clear Previous Outputs</button>
                    <button type="button" class="run-primary" onclick="{run_js}"{disabled_attr}>{esc(run_label)}</button>
                </div>
            </div>
            <div class="analysis-action-status-slot" data-analysis-status-slot></div>
        </div>
    """


def render_modules(state: dict[str, Any]) -> str:
    groups: list[str] = []
    for group in MODULE_GROUPS:
        cards: list[str] = []
        for key in group.keys:
            module = module_for_key(key)
            if module is None:
                continue
            cards.append(
                f"""
                <a class="module-card" href="/module/{esc(module.key)}">
                    <div>
                        <h3>{esc(module.title)}</h3>
                        <p>{esc(TOOL_PURPOSES.get(module.key, module.subtitle))}</p>
                    </div>
                    <span>Open</span>
                </a>
                """
            )
        groups.append(
            f"""
            <section class="tool-menu-group">
                <div class="tool-menu-heading"><h2>{esc(group.title)}</h2><span>{len(cards)} tool{'s' if len(cards) != 1 else ''}</span></div>
                <div class="tool-menu-grid">
                    {''.join(cards)}
                </div>
            </section>
            """
        )
    return "\n".join(groups)


def render_trimming_sample_rows(state: dict[str, Any]) -> str:
    trimming = state.get("trimming", {})
    samples = trimming.get("samples", [])
    if not isinstance(samples, list) or not samples:
        samples = [{"group": "", "read1": "", "read2": ""}]

    rows: list[str] = []
    for index, sample in enumerate(samples):
        if not isinstance(sample, dict):
            sample = {}
        rows.append(
            f"""
            <div class="trim-row" data-trim-row>
                <div class="file-picker">
                    <input data-field="read1" value="{esc(sample.get('read1', ''))}" placeholder="Read 1 FASTQ">
                    <button type="button" class="secondary" onclick="browseTrimFile(this, 'read1')">Browse</button>
                </div>
                <div class="file-picker paired-only">
                    <input data-field="read2" value="{esc(sample.get('read2', ''))}" placeholder="Read 2 FASTQ">
                    <button type="button" class="secondary" onclick="browseTrimFile(this, 'read2')">Browse</button>
                </div>
                <button type="button" class="ghost icon-button" title="Remove row" onclick="removeTrimRow(this)">x</button>
            </div>
            """
        )
    return "\n".join(rows)


def render_trimming_content(state: dict[str, Any]) -> str:
    trimming = state.get("trimming", {})
    trim_galore_path = resolve_trim_galore()
    cutadapt_path = resolve_cutadapt()
    batch_name = str(trimming.get("batch_name", "") or "")
    data_type = str(trimming.get("data_type", "rnaseq") or "rnaseq")
    tool = str(trimming.get("tool", "trim_galore") or "trim_galore")
    mode = trimming.get("mode", "paired")
    quality = trimming.get("quality", 20)
    length = trimming.get("length", 18)
    adapter1 = str(trimming.get("adapter1", "") or "")
    adapter2 = str(trimming.get("adapter2", "") or "")
    cutadapt_error_rate = trimming.get("cutadapt_error_rate", 0.1)
    cutadapt_overlap = trimming.get("cutadapt_overlap", 3)
    cutadapt_max_n = str(trimming.get("cutadapt_max_n", "") or "")
    cutadapt_trim_n = bool(trimming.get("cutadapt_trim_n", True))
    cutadapt_cores = trimming.get("cutadapt_cores", 1)
    cutadapt_pair_filter = str(trimming.get("cutadapt_pair_filter", "any") or "any")
    last_run = trimming.get("last_run", {})
    run_detail = ""
    if isinstance(last_run, dict) and last_run:
        run_detail = (
            f"<p class=\"muted\">Last completed: {esc(last_run.get('name', 'Batch'))} · "
            f"{esc(last_run.get('sample_count', 0))} sample(s) · {esc(last_run.get('finished', ''))}</p>"
        )
    tool_status_parts = [
        f"Trim Galore: {'ready' if trim_galore_path else 'will install when needed'}",
        f"Cutadapt: {'ready' if cutadapt_path else 'will install when needed'}",
    ]
    tool_status = " · ".join(tool_status_parts)

    return f"""
        <section class="panel trim-panel">
            <div class="panel-heading">
                <div>
                    <p class="eyebrow">Batch preprocessing</p>
                    <h3>Load and Trim Sequencing Files</h3>
                    <p>Add raw FASTQ files, choose the batch type, and run trimming. Completed batches become selectable in the relevant analysis tools.</p>
                    {run_detail}
                </div>
                <button type="button" class="secondary" onclick="addTrimRow()">Add Files</button>
            </div>

            <div class="trim-controls">
                <label class="field">
                    <span>Batch name</span>
                    <input id="trim-batch-name" value="{esc(batch_name)}" placeholder="e.g. leaf sRNA batch">
                </label>
                <label class="field">
                    <span>Preprocessing tool</span>
                    <select id="trim-tool" onchange="syncTrimTool()">
                        <option value="trim_galore" {'selected' if tool == 'trim_galore' else ''}>Trim Galore</option>
                        <option value="cutadapt" {'selected' if tool == 'cutadapt' else ''}>Cutadapt</option>
                    </select>
                </label>
                <label class="field">
                    <span>Data type</span>
                    <select id="trim-data-type">
                        <option value="rnaseq" {'selected' if data_type == 'rnaseq' else ''}>dsRNA-seq</option>
                        <option value="srnaseq" {'selected' if data_type == 'srnaseq' else ''}>sRNA-seq</option>
                        <option value="degradome" {'selected' if data_type == 'degradome' else ''}>Degradome-seq</option>
                    </select>
                </label>
                <label class="field">
                    <span>Read layout</span>
                    <select id="trim-mode" onchange="syncTrimMode()">
                        <option value="paired" {'selected' if mode == 'paired' else ''}>Paired-end</option>
                        <option value="single" {'selected' if mode == 'single' else ''}>Single-end</option>
                    </select>
                </label>
                <label class="field">
                    <span>Minimum quality</span>
                    <input id="trim-quality" type="number" min="0" max="40" value="{esc(quality)}">
                </label>
                <label class="field">
                    <span>Minimum length</span>
                    <input id="trim-length" type="number" min="1" value="{esc(length)}">
                </label>
            </div>
            <div class="cutadapt-controls" id="cutadapt-controls">
                <label class="field">
                    <span>Read 1 adapter</span>
                    <input id="cutadapt-adapter1" value="{esc(adapter1)}" placeholder="optional adapter sequence">
                </label>
                <label class="field paired-only">
                    <span>Read 2 adapter</span>
                    <input id="cutadapt-adapter2" value="{esc(adapter2)}" placeholder="optional adapter sequence">
                </label>
                <label class="field">
                    <span>Error rate</span>
                    <input id="cutadapt-error-rate" type="number" min="0" max="0.5" step="0.01" value="{esc(cutadapt_error_rate)}">
                </label>
                <label class="field">
                    <span>Minimum overlap</span>
                    <input id="cutadapt-overlap" type="number" min="1" step="1" value="{esc(cutadapt_overlap)}">
                </label>
                <label class="field">
                    <span>Max N</span>
                    <input id="cutadapt-max-n" type="number" min="0" step="1" value="{esc(cutadapt_max_n)}" placeholder="blank">
                </label>
                <label class="field">
                    <span>Cores</span>
                    <input id="cutadapt-cores" type="number" min="1" step="1" value="{esc(cutadapt_cores)}">
                </label>
                <label class="field paired-only">
                    <span>Pair filter</span>
                    <select id="cutadapt-pair-filter">
                        <option value="any" {'selected' if cutadapt_pair_filter == 'any' else ''}>Discard if either read fails</option>
                        <option value="both" {'selected' if cutadapt_pair_filter == 'both' else ''}>Discard only if both fail</option>
                        <option value="first" {'selected' if cutadapt_pair_filter == 'first' else ''}>Use Read 1 filter only</option>
                    </select>
                </label>
                <label class="toggle-row cutadapt-toggle">
                    <input id="cutadapt-trim-n" type="checkbox" {'checked' if cutadapt_trim_n else ''}>
                    <span>Trim terminal Ns</span>
                </label>
            </div>
            <div id="trim-rows" class="trim-rows">
                {render_trimming_sample_rows(state)}
            </div>

            {render_analysis_action(
                'preprocessing',
                'runTrimming()',
                'Run Batch Trimming',
                extra_html=f'<a class="ghost" href="/reveal?path={url_for(project_output_root(state) / "trimming")}">Open Trimmed Files</a>',
            )}
            <p class="muted">{esc(tool_status)}</p>
        </section>
    """


def srna_output_dir() -> Path:
    return project_output_root() / "srna-mapping"


def srna_dsrna_output_dir() -> Path:
    return project_output_root() / "srna-dsrna-identification"


def srna_control_output_dir() -> Path:
    return project_output_root() / "srna-control-filtering"


def sequencing_sample_rows(state: dict[str, Any], key: str, module_key: str) -> list[dict[str, str]]:
    path, source = effective_path(state, key, module_key)
    if not path:
        return []
    dataset_path = Path(path).expanduser()
    if dataset_path.suffix.lower() in {".csv", ".tsv"} and dataset_path.exists():
        delimiter = "\t" if dataset_path.suffix.lower() == ".tsv" else ","
        with dataset_path.open(newline="", encoding="utf-8") as handle:
            return [dict(row) for row in csv.DictReader(handle, delimiter=delimiter)]
    return [
        {
            "sample_id": slugify_sample_token(dataset_path.name.split(".")[0], "sample"),
            "group": "",
            "replicate": "1",
            "mode": "single",
            "trimmed_read1": str(dataset_path),
            "trimmed_read2": "",
            "status": "preprocessed" if source == "Preprocessed" else "loaded",
        }
    ]


def selected_shared_samples(sample_ids: list[str], state: dict[str, Any]) -> list[dict[str, str]]:
    selected = set(sample_ids)
    rows = [row for row in sequencing_sample_rows(state, "srnaseq", "srna-mapping") if row.get("sample_id") in selected]
    if not rows:
        raise ValueError("Select at least one preprocessed sRNA sample.")
    for row in rows:
        path = row.get("trimmed_read1", "")
        if not path or not Path(path).exists():
            raise ValueError(f"Sample {row.get('sample_id', '')} does not have a readable processed FASTQ.")
    return rows


def selected_dsrna_overlay_samples(sample_ids: list[str], state: dict[str, Any]) -> list[dict[str, str]]:
    selected = set(sample_ids)
    rows = [row for row in sequencing_sample_rows(state, "rnaseq", "srna-mapping") if row.get("sample_id") in selected]
    for row in rows:
        read1 = row.get("trimmed_read1", "")
        read2 = row.get("trimmed_read2", "")
        if not read1 or not Path(read1).exists():
            raise ValueError(f"dsRNA-seq sample {row.get('sample_id', '')} does not have a readable processed Read 1 FASTQ.")
        if not read2 or not Path(read2).exists():
            raise ValueError(f"dsRNA-seq overlay requires paired reads; sample {row.get('sample_id', '')} is missing Read 2.")
    return rows


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


PLOT_DISPLAY_LIMIT = 4
SUMMARY_PLOT_TERMS = (
    "summary",
    "overview",
    "heatmap",
    "heat_map",
    "distribution",
    "barplot",
    "balance",
    "metric",
    "ranking",
)


def table_distinct_value_count(path: Path, column: str) -> int:
    if not path.exists():
        return 0
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            return len({str(row.get(column, "")).strip() for row in csv.DictReader(handle, delimiter=delimiter) if str(row.get(column, "")).strip()})
    except (OSError, csv.Error):
        return 0


def prioritized_plot_paths(paths: Any, *, single_contig: bool = False) -> list[Path]:
    unique: dict[str, Path] = {}
    for value in paths:
        path = Path(value)
        if path.exists() and path.suffix.lower() in {".png", ".svg"}:
            unique[str(path)] = path

    def is_summary(path: Path) -> bool:
        token = path.stem.lower()
        return any(term in token for term in SUMMARY_PLOT_TERMS)

    ordered = sorted(unique.values(), key=lambda path: path.name.lower())
    summaries = [path for path in ordered if is_summary(path)]
    contig_plots = [path for path in ordered if not is_summary(path)]
    preferred = [contig_plots[0], *summaries, *contig_plots[1:]] if single_contig and contig_plots else summaries + contig_plots
    return preferred[:PLOT_DISPLAY_LIMIT]


def render_srna_mapping_content(state: dict[str, Any]) -> str:
    output_dir = tool_output_root("srna-mapping", state)
    plot_links: list[str] = []
    summary_probe = output_dir / "tables" / "contig_summary.tsv"
    contig_count = table_distinct_value_count(summary_probe, "contig")
    all_plots = sorted((output_dir / "plots").glob("*.png"))
    plots_to_show = prioritized_plot_paths(all_plots, single_contig=contig_count == 1)
    for plot in plots_to_show:
        plot_links.append(
            f"""
            <a class="plot-card" href="/reveal?path={url_for(plot)}">
                <strong>{esc(plot.stem)}</strong>
                <img src="{output_file_url(plot)}" alt="{esc(plot.name)}">
            </a>
            """
        )
    plots = "\n".join(plot_links) if plot_links else '<div class="empty-state">Run small RNA mapping to generate plots.</div>'
    summary_tsv = output_dir / "tables" / "contig_summary.tsv"
    alignments_tsv = output_dir / "tables" / "alignments.tsv"
    length_tsv = output_dir / "tables" / "length_distribution.tsv"
    top_unique_tsv = output_dir / "tables" / "top_unique_sRNAs.tsv"
    annotation_tsv = output_dir / "tables" / "siRNA_annotation_matches.tsv"
    fasta = output_dir / "filtered_mapped_sRNAs.fasta"
    return f"""
        <section class="hero">
            <div>
                <p class="eyebrow">Analysis tool</p>
                <h2>Small RNA Mapping</h2>
                <p class="muted">Map grouped small RNA samples to pasted or FASTA references with Bowtie1.</p>
            </div>
            <div class="button-row top-actions">
                <a class="button-link" href="/reveal?path={url_for(output_dir)}">Open Output Folder</a>
            </div>
        </section>
        <div class="notice" id="status">Choose a completed preprocessing batch or an already-preprocessed file, then set mapping options.</div>
        {render_module_dataset_selectors(state, 'srna-mapping')}
        <section class="panel">
            <div>
                <p class="eyebrow">Reference</p>
                <h3>Reference Sequence</h3>
                <textarea id="srna-reference-text" placeholder="Paste one sequence or FASTA records here"></textarea>
                <div class="file-picker">
                    <input id="srna-reference-fasta" placeholder="Or select reference FASTA">
                    <button type="button" class="secondary" onclick="browseSrnaReference()">Browse</button>
                </div>
            </div>
        </section>
        {render_sirna_annotation_selector(state, 'srna-mapping', 'srna-sirna-annotations-fasta')}
        <section class="panel">
            <p class="eyebrow">Mapping options</p>
            <div class="tool-grid">
                <label class="field"><span>Mismatches</span><select id="srna-mismatches"><option selected>0</option><option>1</option><option>2</option><option>3</option></select></label>
                <label class="field"><span>Focus length</span><input id="srna-focus-length" type="number" min="0" placeholder="e.g. 21"></label>
                <label class="field"><span>Plot CPM threshold</span><input id="srna-cpm-threshold" type="number" min="0" step="0.1" value="0"></label>
                <label class="field"><span>Export min CPM</span><input id="srna-export-min-cpm" type="number" min="0" step="0.1" value="0"></label>
                <label class="field"><span>Export length</span><input id="srna-export-length" type="number" min="0" placeholder="optional"></label>
                <label class="field"><span>Top exported unique sRNAs</span><input id="srna-export-top-n" type="number" min="1" step="1" value="20"></label>
                <label class="field"><span>Threads</span><input id="srna-threads" type="number" min="1" value="4"></label>
            </div>
            <div class="button-row">
                <label class="inline-check"><input id="srna-report-all" type="checkbox"><span>Report all valid mappings</span></label>
                <label class="inline-check"><input id="srna-filter-simple" type="checkbox" checked><span>Filter simple mapped reads</span></label>
            </div>
            {render_analysis_action('srna-mapping', 'runSrnaMapping()', 'Run Small RNA Mapping')}
        </section>
        <section class="panel">
            <div class="panel-heading">
                <div><p class="eyebrow">Outputs</p><h3>Plots and Tables</h3></div>
                <div class="button-row">
                    <a class="ghost" href="/reveal?path={url_for(summary_tsv)}">Summary TSV</a>
                    <a class="ghost" href="/reveal?path={url_for(alignments_tsv)}">Alignments TSV</a>
                    <a class="ghost" href="/reveal?path={url_for(length_tsv)}">Lengths TSV</a>
                    <a class="ghost" href="/reveal?path={url_for(top_unique_tsv)}">Top sRNAs TSV</a>
                    <a class="ghost" href="/reveal?path={url_for(annotation_tsv)}">siRNA annotations TSV</a>
                    <a class="ghost" href="/reveal?path={url_for(fasta)}">Filtered FASTA</a>
                </div>
            </div>
            <div class="plot-grid">{plots}</div>
        </section>
    """


def render_srna_control_filtering_content(state: dict[str, Any]) -> str:
    output_dir = tool_output_root("srna-control-filtering", state)
    mapped_unique = tool_output_root("srna-mapping", state) / "filtered_mapped_sRNAs.fasta"
    default_unique = mapped_unique
    default_control_srna = str(DEFAULT_CONTROL_SRNA_FASTA) if DEFAULT_CONTROL_SRNA_FASTA.exists() else ""
    summary_tsv = output_dir / "tables" / "control_mapping_summary.tsv"
    complexity_tsv = output_dir / "tables" / "complexity_summary.tsv"
    collapse_summary_tsv = output_dir / "tables" / "short_contained_collapse_summary.tsv"
    collapse_removed_tsv = output_dir / "tables" / "short_contained_removed.tsv"
    clean_fasta = output_dir / "filtered_unique_sRNAs.clean.fasta"
    dsrna_sam = output_dir / "tables" / "dsrna_reference.mapped.sam"
    genome_sam = output_dir / "tables" / "control_genome.mapped.sam"
    control_srna_sam = output_dir / "tables" / "control_srna.mapped.sam"
    plot_links: list[str] = []
    all_plots = sorted((output_dir / "plots").glob("*.png"))
    single_contig = len([path for path in all_plots if "single_overlay_control_hits" in path.stem]) == 1
    for plot in prioritized_plot_paths(all_plots, single_contig=single_contig):
        plot_links.append(
            f"""
            <a class="plot-card" href="/reveal?path={url_for(plot)}">
                <strong>{esc(plot.stem)}</strong>
                <img src="{output_file_url(plot)}" alt="{esc(plot.name)}">
            </a>
            """
        )
    plots = "\n".join(plot_links) if plot_links else '<div class="empty-state">Run control mapping to generate dsRNA locus plots.</div>'
    return f"""
        <section class="hero">
            <div>
                <p class="eyebrow">Analysis tool</p>
                <h2>Small RNA Control Mapping & Filtering</h2>
                <p class="muted">Screen unique mapped sRNAs against control references and flag sequence simplicity before downstream use.</p>
            </div>
            <div class="button-row top-actions">
                <a class="ghost" href="/module/srna-mapping">Small RNA Mapping</a>
                <a class="button-link" href="/reveal?path={url_for(output_dir)}">Open Output Folder</a>
            </div>
        </section>
        <div class="notice" id="status">Use the unique FASTA from Small RNA Mapping, map it to the dsRNA reference locus, then mark reads that also hit control references.</div>
        <section class="panel">
            <div class="two-col">
                <div>
                    <p class="eyebrow">Input</p>
                    <h3>Unique sRNA FASTA</h3>
                    <textarea id="control-unique-text" placeholder="Paste unique sRNA FASTA records here"></textarea>
                    <div class="file-picker">
                        <input id="control-unique-fasta" value="{esc(str(default_unique))}" placeholder="filtered_mapped_sRNAs.fasta">
                        <button type="button" class="secondary" onclick="browseSrnaControlFile('control-unique-fasta')">Browse</button>
                    </div>
                    <textarea id="control-dsrna-text" placeholder="Paste dsRNA reference sequence or FASTA records here"></textarea>
                    <div class="file-picker">
                        <input id="control-dsrna-fasta" placeholder="dsRNA reference FASTA for locus plots">
                        <button type="button" class="secondary" onclick="browseSrnaControlFile('control-dsrna-fasta')">Browse</button>
                    </div>
                </div>
                <div>
                    <p class="eyebrow">Control references</p>
                    <h3>Genome and sRNA FASTA</h3>
                    <textarea id="control-genome-text" placeholder="Paste control genome FASTA records here"></textarea>
                    <div class="file-picker">
                        <input id="control-genome-fasta" placeholder="Control reference genome FASTA">
                        <button type="button" class="secondary" onclick="browseSrnaControlFile('control-genome-fasta')">Browse</button>
                    </div>
                    <textarea id="control-srna-text" placeholder="Paste control sRNA FASTA records here"></textarea>
                    <div class="file-picker">
                        <input id="control-srna-fasta" value="{esc(default_control_srna)}" placeholder="Control sRNA FASTA">
                        <button type="button" class="secondary" onclick="browseSrnaControlFile('control-srna-fasta')">Browse</button>
                    </div>
                </div>
            </div>
        </section>
        <section class="panel">
            <p class="eyebrow">Filtering options</p>
            <div class="tool-grid">
                <label class="field"><span>dsRNA mapping mismatches</span><select id="control-dsrna-mismatches"><option selected>0</option><option>1</option><option>2</option><option>3</option></select></label>
                <label class="field"><span>Control search mismatches</span><select id="control-reference-mismatches"><option selected>0</option><option>1</option><option>2</option><option>3</option></select></label>
                <label class="field"><span>sRNA length filter</span><input id="control-length-filter" type="number" min="0" placeholder="optional, e.g. 21"></label>
                <label class="field"><span>Threads</span><input id="control-threads" type="number" min="1" value="4"></label>
                <label class="field"><span>Index build mode</span><select id="control-index-mode"><option value="auto" selected>Automatic / balanced</option><option value="fast">Faster build, more memory</option><option value="lowmem">Low memory, slower build</option></select></label>
            </div>
            <div class="button-row">
                <label class="inline-check"><input id="control-remove-low-complexity" type="checkbox" checked><span>Remove low-complexity reads from clean FASTA</span></label>
                <label class="inline-check"><input id="control-remove-control-mappers" type="checkbox" checked><span>Remove reads mapping to controls from clean FASTA</span></label>
                <label class="inline-check"><input id="control-collapse-contained" type="checkbox"><span>Collapse shorter reads contained in longer reads</span></label>
            </div>
            {render_analysis_action('srna-control-filtering', 'runSrnaControlFiltering()', 'Run Control Mapping & Filtering')}
        </section>
        <section class="panel">
            <div class="panel-heading">
                <div><p class="eyebrow">Outputs</p><h3>Tables and Clean FASTA</h3></div>
                <div class="button-row">
                    <a class="ghost" href="/reveal?path={url_for(summary_tsv)}">Summary TSV</a>
                    <a class="ghost" href="/reveal?path={url_for(complexity_tsv)}">Complexity TSV</a>
                    <a class="ghost" href="/reveal?path={url_for(collapse_summary_tsv)}">Collapse Summary</a>
                    <a class="ghost" href="/reveal?path={url_for(collapse_removed_tsv)}">Removed Contained</a>
                    <a class="ghost" href="/reveal?path={url_for(clean_fasta)}">Clean FASTA</a>
                    <a class="ghost" href="/reveal?path={url_for(dsrna_sam)}">dsRNA SAM</a>
                    <a class="ghost" href="/reveal?path={url_for(genome_sam)}">Genome SAM</a>
                    <a class="ghost" href="/reveal?path={url_for(control_srna_sam)}">Control sRNA SAM</a>
                </div>
            </div>
            <div class="plot-grid">{plots}</div>
        </section>
    """


def render_preprocessing_content(state: dict[str, Any]) -> str:
    module = module_for_key("preprocessing")
    if module is None:
        return ""

    output_dir = project_output_root(state) / "trimming"

    return f"""
        <section class="hero">
            <div>
                <p class="eyebrow">Pipeline section</p>
                <h2>{esc(module.title)}</h2>
                <p class="muted">{esc(module.subtitle)}</p>
            </div>
            <div class="button-row top-actions">
                <a class="button-link" href="/reveal?path={url_for(output_dir)}">Open Trimmed Files</a>
            </div>
        </section>
        <div class="notice" id="status">Add raw files and run one batch. Already-trimmed files are selected directly from the analysis tool that uses them.</div>
        {render_trimming_content(state)}
    """


def normalize_trimming_payload(payload: dict[str, Any]) -> dict[str, Any]:
    batch_name = str(payload.get("batch_name", "")).strip()
    if not batch_name:
        raise ValueError("Enter a batch name.")
    data_type = str(payload.get("data_type", "rnaseq")).strip()
    if data_type not in {key for key, _, _ in SEQUENCING_DATASETS}:
        raise ValueError("Choose dsRNA-seq, sRNA-seq, or degradome-seq as the data type.")
    tool = str(payload.get("tool", "trim_galore")).strip().lower().replace("-", "_")
    if tool not in {"trim_galore", "cutadapt"}:
        raise ValueError("Choose Trim Galore or Cutadapt as the preprocessing tool.")
    mode = str(payload.get("mode", "paired")).strip().lower()
    if mode not in {"paired", "single"}:
        raise ValueError("Choose paired-end or single-end trimming.")

    try:
        quality = int(payload.get("quality", 20))
        length = int(payload.get("length", 18))
        cutadapt_overlap = int(payload.get("cutadapt_overlap", 3))
        cutadapt_cores = int(payload.get("cutadapt_cores", 1))
    except (TypeError, ValueError) as exc:
        raise ValueError("Quality, length, overlap, and cores must be whole numbers.") from exc

    try:
        cutadapt_error_rate = float(payload.get("cutadapt_error_rate", 0.1))
    except (TypeError, ValueError) as exc:
        raise ValueError("Cutadapt error rate must be a number.") from exc

    if quality < 0 or quality > 40:
        raise ValueError("Minimum quality must be between 0 and 40.")
    if length < 1:
        raise ValueError("Minimum length must be at least 1.")
    if not 0 <= cutadapt_error_rate <= 0.5:
        raise ValueError("Cutadapt error rate must be between 0 and 0.5.")
    if cutadapt_overlap < 1:
        raise ValueError("Cutadapt minimum overlap must be at least 1.")
    if cutadapt_cores < 1:
        raise ValueError("Cutadapt cores must be at least 1.")

    cutadapt_max_n = str(payload.get("cutadapt_max_n", "")).strip()
    if cutadapt_max_n:
        try:
            max_n_number = float(cutadapt_max_n)
        except ValueError as exc:
            raise ValueError("Cutadapt max N must be a number or blank.") from exc
        if max_n_number < 0:
            raise ValueError("Cutadapt max N must be 0 or greater.")

    cutadapt_pair_filter = str(payload.get("cutadapt_pair_filter", "any")).strip().lower()
    if cutadapt_pair_filter not in {"any", "both", "first"}:
        raise ValueError("Cutadapt pair filter must be any, both, or first.")

    samples = payload.get("samples", [])
    if not isinstance(samples, list) or not samples:
        raise ValueError("Add at least one sample row before trimming.")

    clean_samples: list[dict[str, str]] = []
    for index, sample in enumerate(samples, start=1):
        if not isinstance(sample, dict):
            continue
        read1 = str(sample.get("read1", "")).strip()
        read2 = str(sample.get("read2", "")).strip()
        if not read1:
            raise ValueError(f"Sample row {index} is missing Read 1.")
        if mode == "paired" and not read2:
            raise ValueError(f"Sample row {index} is missing Read 2 for paired-end trimming.")
        if not Path(read1).expanduser().exists():
            raise ValueError(f"Read 1 file does not exist for row {index}: {read1}")
        if mode == "paired" and not Path(read2).expanduser().exists():
            raise ValueError(f"Read 2 file does not exist for row {index}: {read2}")
        clean_samples.append({"group": "", "read1": read1, "read2": read2, "skip_preprocessing": False})

    if not clean_samples:
        raise ValueError("Add at least one complete sample row before trimming.")

    return {
        "batch_name": batch_name,
        "data_type": data_type,
        "tool": tool,
        "mode": mode,
        "quality": quality,
        "length": length,
        "adapter1": str(payload.get("adapter1", "")).strip(),
        "adapter2": str(payload.get("adapter2", "")).strip(),
        "cutadapt_error_rate": cutadapt_error_rate,
        "cutadapt_overlap": cutadapt_overlap,
        "cutadapt_max_n": cutadapt_max_n,
        "cutadapt_trim_n": bool(payload.get("cutadapt_trim_n", True)),
        "cutadapt_cores": cutadapt_cores,
        "cutadapt_pair_filter": cutadapt_pair_filter,
        "samples": clean_samples,
    }


def save_trimming_settings(state: dict[str, Any], settings: dict[str, Any]) -> None:
    current = dict(state.get("trimming", {}))
    current.update({key: value for key, value in settings.items() if not key.startswith("_")})
    state["trimming"] = current
    save_state(state)


def run_trimming(settings: dict[str, Any], state: dict[str, Any], job_id: str | None = None) -> dict[str, Any]:
    needs_trimming = any(not sample.get("skip_preprocessing") for sample in settings["samples"])
    tool = str(settings.get("tool", "trim_galore"))
    tool_title = trimming_tool_title(tool)
    executable = ""
    if needs_trimming:
        if tool == "cutadapt":
            executable = resolve_cutadapt()
            if not executable:
                append_job_log(job_id, "Cutadapt not found; starting automatic installation.")
                executable = install_cutadapt(job_id)["path"]
        else:
            executable = resolve_trim_galore()
            if not executable:
                append_job_log(job_id, "Trim Galore not found; starting automatic installation.")
                executable = install_trim_galore(job_id)["path"]
        if executable:
            append_job_log(job_id, f"Using {tool_title} at {executable}.")
        else:
            raise RuntimeError(f"{tool_title} is not available.")

    batch_id = "latest"
    outdir = trimming_output_dir(batch_id)
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = trimming_manifest_path(batch_id)
    paired = settings["mode"] == "paired"
    sample_counts: dict[str, int] = {}
    rows: list[dict[str, str]] = []

    for index, sample in enumerate(settings["samples"], start=1):
        filename = Path(sample["read1"]).name
        sample_label = re.sub(r"(?:_R?1)?\.(?:fastq|fq)(?:\.gz)?$", "", filename, flags=re.IGNORECASE)
        sample_slug = slugify_sample_token(sample_label, f"sample_{index}")
        sample_counts[sample_slug] = sample_counts.get(sample_slug, 0) + 1
        sample_id = sample_slug if sample_counts[sample_slug] == 1 else f"{sample_slug}_{sample_counts[sample_slug]}"
        append_job_log(job_id, f"Preparing {sample_id}.")

        if tool == "cutadapt":
            trimmed_read1 = cutadapt_fastq_path(sample_id, paired, 1, outdir)
            trimmed_read2 = cutadapt_fastq_path(sample_id, True, 2, outdir) if paired else Path("")
            command = [
                executable,
                "--quality-cutoff",
                str(settings["quality"]),
                "--minimum-length",
                str(settings["length"]),
                "--cores",
                str(settings.get("cutadapt_cores", 1)),
                "--error-rate",
                str(settings.get("cutadapt_error_rate", 0.1)),
                "--overlap",
                str(settings.get("cutadapt_overlap", 3)),
                "--output",
                str(trimmed_read1),
            ]
            if paired:
                command.extend(["--paired-output", str(trimmed_read2)])
            if settings.get("adapter1"):
                command.extend(["--adapter", settings["adapter1"]])
            if paired and settings.get("adapter2"):
                command.extend(["--adapter2", settings["adapter2"]])
            if settings.get("cutadapt_max_n"):
                command.extend(["--max-n", str(settings["cutadapt_max_n"])])
            if settings.get("cutadapt_trim_n", True):
                command.append("--trim-n")
            if paired:
                command.extend(["--pair-filter", str(settings.get("cutadapt_pair_filter", "any"))])
                command.extend([sample["read1"], sample["read2"]])
            else:
                command.append(sample["read1"])
        else:
            command = [
                executable,
                "--quality",
                str(settings["quality"]),
                "--length",
                str(settings["length"]),
                "--gzip",
                "--output_dir",
                str(outdir),
                "--basename",
                sample_id,
            ]
            if settings.get("adapter1"):
                command.extend(["--adapter", settings["adapter1"]])
            if paired:
                if settings.get("adapter2"):
                    command.extend(["--adapter2", settings["adapter2"]])
                command.extend(["--paired", sample["read1"], sample["read2"]])
            else:
                command.append(sample["read1"])

            trimmed_read1 = first_existing_path(trimmed_fastq_candidates(sample_id, paired, 1, outdir))
            trimmed_read2 = first_existing_path(trimmed_fastq_candidates(sample_id, True, 2, outdir)) if paired else Path("")

        completed = run_tracked_command(command, f"{tool_title} {sample_id}", job_id)
        if completed.returncode != 0:
            detail = (completed.stderr or completed.stdout or "").strip()
            raise RuntimeError(f"{tool_title} failed for {sample_id}: {detail}")

        rows.append(
            {
                "sample_id": sample_id,
                "group": "",
                "replicate": str(index),
                "batch_id": batch_id,
                "tool": tool_title,
                "data_type": settings["data_type"],
                "mode": settings["mode"],
                "raw_read1": sample["read1"],
                "raw_read2": sample["read2"],
                "trimmed_read1": str(trimmed_read1),
                "trimmed_read2": str(trimmed_read2) if paired else "",
                "command": " ".join(command),
                "status": "trimmed",
            }
        )
        append_job_log(job_id, f"Wrote {tool_title} outputs for {sample_id}.")

    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sample_id",
                "group",
                "replicate",
                "batch_id",
                "tool",
                "data_type",
                "mode",
                "raw_read1",
                "raw_read2",
                "trimmed_read1",
                "trimmed_read2",
                "command",
                "status",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    settings["manifest"] = str(manifest_path)
    settings["last_run"] = {
        "status": "finished",
        "name": settings["batch_name"],
        "batch_id": batch_id,
        "sample_count": len(rows),
        "finished": now_text(),
    }
    state["trimmed_batches"] = {batch_id: {
        "id": batch_id,
        "name": settings["batch_name"],
        "tool": tool,
        "data_type": settings["data_type"],
        "mode": settings["mode"],
        "manifest": str(manifest_path),
        "outdir": str(outdir),
        "sample_count": len(rows),
        "finished": settings["last_run"]["finished"],
        "project_root": str(project_output_root(state)),
    }}
    save_trimming_settings(state, settings)
    append_job_log(job_id, f"Wrote trimmed sample manifest to {manifest_path}.")
    return {"batch_id": batch_id, "manifest": str(manifest_path), "sample_count": len(rows), "outdir": str(outdir)}


def start_trimming_job(settings: dict[str, Any]) -> str:
    tool_title = trimming_tool_title(str(settings.get("tool", "trim_galore")))
    job_id = create_pipeline_job(f"{tool_title} batch")

    def worker() -> None:
        try:
            result = run_trimming(settings, load_state(), job_id)
        except Exception as exc:
            finish_pipeline_job(job_id, "failed", {"message": str(exc)}, f"{tool_title} batch failed: {exc}")
            return
        finish_pipeline_job(
            job_id,
            "finished",
            result,
            f"{tool_title} finished for {result.get('sample_count', 0)} sample(s).",
        )

    threading.Thread(target=worker, name=f"inci-trimming-{job_id[:8]}", daemon=True).start()
    return job_id


def start_trim_galore_install_job() -> str:
    job_id = create_pipeline_job("Trim Galore install")

    def worker() -> None:
        try:
            result = install_trim_galore(job_id)
        except Exception as exc:
            finish_pipeline_job(job_id, "failed", {"message": str(exc)}, f"Trim Galore installation failed: {exc}")
            return
        finish_pipeline_job(job_id, "finished", result, result.get("message", "Trim Galore is ready."))

    threading.Thread(target=worker, name=f"inci-install-trimgalore-{job_id[:8]}", daemon=True).start()
    return job_id


def start_srna_mapping_job(settings: dict[str, Any]) -> str:
    job_id = create_pipeline_job("sRNA mapping")
    prepare_tool_output_dir("srna-mapping")

    def worker() -> None:
        try:
            result = run_srna_mapping(settings, load_state(), job_id)
        except Exception as exc:
            append_job_log(job_id, traceback.format_exc(), "error")
            finish_pipeline_job(job_id, "failed", {"message": str(exc)}, f"sRNA mapping failed: {exc}")
            return
        finish_pipeline_job(job_id, "finished", result, f"sRNA mapping finished with {result.get('contigs', 0)} output record(s).")

    threading.Thread(target=worker, name=f"inci-srna-mapping-{job_id[:8]}", daemon=True).start()
    return job_id


def parse_srna_mapping_payload(payload: dict[str, Any]) -> dict[str, Any]:
    sample_ids = payload.get("sample_ids", [])
    if not isinstance(sample_ids, list) or not sample_ids:
        raise ValueError("Select at least one sRNA sample.")
    dsrna_sample_ids = payload.get("dsrna_sample_ids", [])
    if not isinstance(dsrna_sample_ids, list):
        dsrna_sample_ids = []
    mismatches = int(payload.get("mismatches", 0))
    if mismatches not in {0, 1, 2, 3}:
        raise ValueError("Mismatches must be 0, 1, 2, or 3.")
    return {
        "sample_ids": [str(item) for item in sample_ids],
        "dsrna_sample_ids": [str(item) for item in dsrna_sample_ids],
        "include_dsrna_overlay": bool(payload.get("include_dsrna_overlay", False)),
        "reference_text": str(payload.get("reference_text", "")).strip(),
        "reference_fasta": str(payload.get("reference_fasta", "")).strip(),
        "sirna_annotations_fasta": str(payload.get("sirna_annotations_fasta", "")).strip(),
        "mismatches": mismatches,
        "report_all": bool(payload.get("report_all", False)),
        "filter_simple": bool(payload.get("filter_simple", True)),
        "focus_length": int(payload.get("focus_length") or 0),
        "cpm_threshold": float(payload.get("cpm_threshold") or 0),
        "export_min_cpm": float(payload.get("export_min_cpm") or 0),
        "export_length": int(payload.get("export_length") or 0),
        "export_top_n": max(1, int(payload.get("export_top_n") or 20)),
        "threads": max(1, int(payload.get("threads") or 4)),
    }


def write_srna_reference(settings: dict[str, Any], outdir: Path) -> Path:
    if settings["reference_fasta"]:
        path = Path(settings["reference_fasta"]).expanduser()
        if not path.exists():
            raise ValueError(f"Reference FASTA does not exist: {path}")
        return path
    text = settings["reference_text"].strip()
    if not text:
        raise ValueError("Paste a reference sequence or choose a FASTA file.")
    ref_path = outdir / "reference.fasta"
    if text.startswith(">"):
        ref_path.write_text(text + "\n", encoding="utf-8")
    else:
        seq = re.sub(r"[^A-Za-z]", "", text).upper().replace("U", "T")
        if not seq:
            raise ValueError("Pasted reference sequence has no nucleotide bases.")
        ref_path.write_text(f">pasted_reference\n{seq}\n", encoding="utf-8")
    return ref_path


def fastq_records(path: Path):
    import sRNA_identification as srna_core

    yield from srna_core.iter_fastq(path)


def fast_simple_srna_reason(sequence: str) -> str | None:
    seq = sequence.upper().replace("U", "T")
    length = len(seq)
    if not length:
        return "empty"
    counts = {base: seq.count(base) for base in "ACGT"}
    acgt_total = sum(counts.values())
    if acgt_total == 0:
        return "no_acgt"
    if seq.count("N") / length > 0.10:
        return "too_many_n"
    nonzero = [count for count in counts.values() if count]
    if len(nonzero) <= 2:
        return "one_or_two_base_alphabet"
    if max(nonzero) / acgt_total >= 0.80:
        return "dominant_single_base"
    longest = 1
    run = 1
    for prev, current in zip(seq, seq[1:]):
        run = run + 1 if current == prev else 1
        longest = max(longest, run)
    if longest >= 8 or longest / length >= 0.60:
        return "long_homopolymer"
    for k, cutoff in ((2, 0.70), (3, 0.70)):
        if length >= k:
            kmers: dict[str, int] = {}
            total = length - k + 1
            for idx in range(total):
                kmer = seq[idx : idx + k]
                kmers[kmer] = kmers.get(kmer, 0) + 1
            if max(kmers.values()) / total >= cutoff:
                return f"dominant_{k}mer"
    return None


def filter_srna_fastq(source: Path, dest: Path, filter_simple: bool, job_id: str | None = None, label: str = "") -> tuple[int, int]:
    total = 0
    kept = 0
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", encoding="utf-8") as out:
        for record in fastq_records(source):
            total += 1
            seq = record.sequence.upper().replace("U", "T")
            if filter_simple:
                reason = fast_simple_srna_reason(seq)
                if reason:
                    continue
            kept += 1
            out.write(f"{record.name}\n{seq}\n{record.plus}\n{record.quality}\n")
            if total % 500000 == 0:
                append_job_log(job_id, f"Filtered {total:,} reads for {label or source.name}; kept {kept:,}.")
    if total == 0:
        raise ValueError(f"No reads found in {source}")
    append_job_log(job_id, f"Finished read filtering for {label or source.name}: kept {kept:,} of {total:,} reads.")
    return total, kept


def count_srna_fastq(source: Path, job_id: str | None = None, label: str = "") -> int:
    total = 0
    for _ in fastq_records(source):
        total += 1
        if total % 1000000 == 0:
            append_job_log(job_id, f"Counted {total:,} reads for {label or source.name}.")
    if total == 0:
        raise ValueError(f"No reads found in {source}")
    append_job_log(job_id, f"Using {total:,} preprocessed reads for CPM normalization in {label or source.name}.")
    return total


def read_reference_lengths(fasta: Path) -> dict[str, int]:
    import sRNA_identification as srna_core

    return {contig.name: len(contig.sequence) for contig in srna_core.read_fasta(fasta)}


def read_reference_sequences(fasta: Path) -> dict[str, str]:
    import sRNA_identification as srna_core

    return {contig.name: contig.sequence.upper().replace("U", "T") for contig in srna_core.read_fasta(fasta)}


def parse_srna_sam(
    sam_path: Path,
    sample: dict[str, str],
    denominator: int,
    focus_length: int,
    filter_simple: bool = True,
    job_id: str | None = None,
    label: str = "",
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    import sRNA_identification as srna_core

    alignments: list[dict[str, Any]] = []
    length_counts: dict[tuple[int, str], int] = {}
    raw_mapped_alignments = 0
    simple_filtered_alignments = 0
    focus_filtered_alignments = 0
    with sam_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            try:
                flag = int(fields[1])
                start = int(fields[3])
            except ValueError:
                continue
            if flag & 4:
                continue
            seq = fields[9].upper().replace("U", "T")
            length = len(seq)
            strand = "antisense" if flag & 16 else "sense"
            nm = srna_core.mismatch_count(fields) or 0
            raw_mapped_alignments += 1
            if raw_mapped_alignments % 1_000_000 == 0:
                append_job_log(job_id, f"Parsed {raw_mapped_alignments:,} mapped SAM alignment(s) for {label or sample['sample_id']}; retained {len(alignments):,}.")
            if filter_simple and fast_simple_srna_reason(seq):
                simple_filtered_alignments += 1
                continue
            length_counts[(length, strand)] = length_counts.get((length, strand), 0) + 1
            if focus_length and length != focus_length:
                focus_filtered_alignments += 1
                continue
            ref_len = srna_core.cigar_reference_length(fields[5], seq)
            alignments.append(
                {
                    "sample_id": sample["sample_id"],
                    "group": sample["group"],
                    "replicate": sample["replicate"],
                    "contig": fields[2],
                    "start": start,
                    "end": start + ref_len - 1,
                    "length": length,
                    "strand": strand,
                    "mismatches": nm,
                    "sequence": seq,
                    "cpm": 1_000_000.0 / denominator if denominator else 0.0,
                }
            )
    stats = {
        "sample_id": sample["sample_id"],
        "group": sample["group"],
        "replicate": sample["replicate"],
        "denominator_reads": denominator,
        "mapped_alignments": len(alignments),
        "raw_mapped_alignments": raw_mapped_alignments,
        "simple_filtered_alignments": simple_filtered_alignments,
        "focus_filtered_alignments": focus_filtered_alignments,
        "length_counts": length_counts,
    }
    return alignments, stats


def group_mean_sem(values: list[float]) -> tuple[float, float]:
    if not values:
        return 0.0, 0.0
    import numpy as np

    arr = np.array(values, dtype=float)
    mean = float(np.mean(arr))
    sem = float(np.std(arr, ddof=1) / np.sqrt(len(arr))) if len(arr) > 1 else 0.0
    return mean, sem


def run_srna_mapping(settings: dict[str, Any], state: dict[str, Any], job_id: str | None = None) -> dict[str, Any]:
    apply_global_plot_settings(state)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    samples = selected_shared_samples(settings["sample_ids"], state)
    dsrna_samples = (
        selected_dsrna_overlay_samples(settings.get("dsrna_sample_ids", []), state)
        if settings.get("include_dsrna_overlay")
        else []
    )
    outdir = srna_output_dir()
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "plots").mkdir(exist_ok=True)
    (outdir / "tables").mkdir(exist_ok=True)
    (outdir / "fastq").mkdir(exist_ok=True)

    all_alignments: list[dict[str, Any]] = []
    sample_stats: list[dict[str, Any]] = []

    reference = write_srna_reference(settings, outdir)
    ref_lengths = read_reference_lengths(reference)
    reference_sequences = read_reference_sequences(reference)
    annotation_path = optional_sirna_annotation_path(settings.get("sirna_annotations_fasta", ""))
    sirna_annotations = find_sirna_annotations(annotation_path, reference_sequences) if annotation_path else []
    if annotation_path:
        remember_sirna_annotation_path(state, "srna-mapping", annotation_path)
        write_sirna_annotation_table(outdir / "tables" / "siRNA_annotation_matches.tsv", sirna_annotations)
        append_job_log(job_id, f"Found {len(sirna_annotations):,} exact siRNA annotation match(es) for positional plots.")
    index_prefix = outdir / "bowtie_index" / reference.stem
    index_prefix.parent.mkdir(parents=True, exist_ok=True)
    expected_index = [Path(f"{index_prefix}.{suffix}.ebwt") for suffix in ("1", "2", "3", "4", "rev.1", "rev.2")]
    if not all(path.exists() for path in expected_index):
        completed = run_tracked_command(["bowtie-build", "--threads", str(settings["threads"]), str(reference), str(index_prefix)], "Build sRNA Bowtie index", job_id)
        if completed.returncode != 0:
            raise RuntimeError(completed.stderr or completed.stdout or "Bowtie index build failed.")

    bowtie = resolve_executable("bowtie-align-s") or resolve_executable("bowtie")
    if not bowtie:
        raise ValueError("Bowtie1 was not found on PATH.")

    for sample in samples:
        append_job_log(job_id, f"Preparing sRNA sample {sample['sample_id']} ({sample.get('group', 'ungrouped')} replicate {sample.get('replicate', '')}).")
        source_fastq = Path(sample["trimmed_read1"])
        filtered = source_fastq
        denominator = count_srna_fastq(source_fastq, job_id, sample["sample_id"])
        command = [
            bowtie,
            "-S",
            "--no-unal",
            "-v",
            str(settings["mismatches"]),
            "-p",
            str(settings["threads"]),
            "-x",
            str(index_prefix),
            str(filtered),
            str(outdir / "tables" / f"{sample['sample_id']}.sam"),
        ]
        if settings["report_all"]:
            command.insert(2, "-a")
            append_job_log(job_id, f"Report-all mapping is enabled for {sample['sample_id']}; this can be slow and produce large SAM files.")
        else:
            append_job_log(job_id, f"Mapping sRNA sample {sample['sample_id']} with Bowtie1 best/valid default reporting.")
        completed = run_tracked_command(command, f"sRNA Bowtie {sample['sample_id']}", job_id)
        if completed.returncode != 0:
            raise RuntimeError(completed.stderr or completed.stdout or f"Bowtie failed for {sample['sample_id']}")
        append_job_log(job_id, f"Parsing Bowtie SAM for {sample['sample_id']}.")
        alignments, stats = parse_srna_sam(
            outdir / "tables" / f"{sample['sample_id']}.sam",
            sample,
            denominator,
            settings["focus_length"],
            settings.get("filter_simple", True),
            job_id,
            sample["sample_id"],
        )
        all_alignments.extend(alignments)
        sample_stats.append(stats)
        append_job_log(
            job_id,
            f"Done sRNA sample {sample['sample_id']}: {len(alignments):,} alignment(s) retained; "
            f"{stats['simple_filtered_alignments']:,} simple mapped alignment(s) filtered.",
        )

    append_job_log(job_id, "Summarizing sRNA mapping coverage across contigs and biological groups.")
    summary_rows, coverage_by_group, length_rows = summarize_srna_mapping(all_alignments, sample_stats, ref_lengths)
    dsrna_coverage_by_group, dsrna_stats_rows = load_dsrna_overlay_coverage(
        dsrna_samples,
        reference,
        ref_lengths,
        outdir,
        threads=settings["threads"],
        job_id=job_id,
    )
    append_job_log(job_id, "Rendering sRNA positional coverage plots.")
    plot_srna_outputs(outdir, summary_rows, coverage_by_group, length_rows, settings, dsrna_coverage_by_group, sirna_annotations)
    length_by_contig_rows = plot_mapped_srna_length_distributions(outdir, all_alignments)
    append_job_log(job_id, "Writing sRNA mapping tables and filtered FASTA.")
    export_records = srna_export_records(all_alignments, settings, sample_stats)
    fasta_path = write_filtered_srna_fasta(outdir / "filtered_mapped_sRNAs.fasta", export_records)
    write_tsv(outdir / "tables" / "contig_summary.tsv", summary_rows)
    write_tsv(outdir / "tables" / "alignments.tsv", all_alignments)
    write_tsv(outdir / "tables" / "length_distribution.tsv", length_rows)
    write_tsv(outdir / "tables" / "mapped_length_distribution_by_contig.tsv", length_by_contig_rows)
    write_tsv(outdir / "tables" / "sample_summary.tsv", srna_sample_summary_rows(sample_stats))
    write_tsv(outdir / "tables" / "unique_sRNAs.tsv", unique_srna_table_rows(all_alignments, settings, sample_stats))
    write_tsv(outdir / "tables" / "top_unique_sRNAs.tsv", export_records)
    if dsrna_stats_rows:
        write_tsv(outdir / "tables" / "dsrna_overlay_sample_summary.tsv", dsrna_stats_rows)
    append_job_log(job_id, f"Wrote sRNA mapping outputs to {outdir}.")
    return {
        "outdir": str(outdir),
        "summary_tsv": str(outdir / "tables" / "contig_summary.tsv"),
        "fasta": str(fasta_path),
        "plots": str(outdir / "plots"),
        "contigs": len(summary_rows),
        "dsrna_overlay_samples": len(dsrna_stats_rows),
        "sirna_annotation_matches": len(sirna_annotations),
        "sirna_annotation_tsv": str(outdir / "tables" / "siRNA_annotation_matches.tsv") if annotation_path else "",
    }


def parse_srna_dsrna_payload(payload: dict[str, Any]) -> dict[str, Any]:
    sample_ids = payload.get("sample_ids", [])
    if not isinstance(sample_ids, list):
        sample_ids = []
    if not sample_ids:
        raise ValueError("Select at least one sRNA-seq sample.")

    def integer(name: str, default: int, minimum: int) -> int:
        try:
            value = int(payload.get(name, default) or default)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{name.replace('_', ' ').capitalize()} must be an integer.") from exc
        if value < minimum:
            raise ValueError(f"{name.replace('_', ' ').capitalize()} must be at least {minimum}.")
        return value

    mismatches = integer("mismatches", 0, 0)
    if mismatches > 3:
        raise ValueError("Mismatches must be 0, 1, 2, or 3.")
    return {
        "sample_ids": [str(item) for item in sample_ids],
        "scoring_group": str(payload.get("scoring_group", "")).strip(),
        "bin_size": integer("bin_size", 250, 25),
        "plot_context_bp": integer("plot_context_bp", 500, 0),
        "top_n": integer("top_n", 20, 1),
        "mismatches": mismatches,
        "focus_length": integer("focus_length", 0, 0),
        "filter_simple": bool(payload.get("filter_simple", True)),
        "max_multimappers": integer("max_multimappers", 50, 1),
        "report_all": bool(payload.get("report_all", False)),
        "threads": integer("threads", 4, 1),
    }


def selected_srna_dsrna_samples(sample_ids: list[str], state: dict[str, Any]) -> list[dict[str, str]]:
    selected = set(sample_ids)
    rows = [row for row in stored_tool_samples(state, "srna-dsrna-identification", "srnaseq") if row.get("sample_id") in selected]
    if not rows:
        rows = [row for row in sequencing_sample_rows(state, "srnaseq", "srna-dsrna-identification") if row.get("sample_id") in selected]
    if not rows:
        raise ValueError("Select at least one preprocessed sRNA sample.")
    for row in rows:
        path = row.get("trimmed_read1", "")
        if not path or not Path(path).exists():
            raise ValueError(f"Sample {row.get('sample_id', '')} does not have a readable processed FASTQ.")
    return rows


def read_fasta_sequences(path: Path) -> dict[str, str]:
    import sRNA_identification as srna_core

    return {record.name: record.sequence.upper().replace("U", "T") for record in srna_core.read_fasta(path)}


def score_srna_dsrna_bins(
    coverage_by_group: dict[str, Any],
    ref_lengths: dict[str, int],
    settings: dict[str, Any],
    scoring_group: str,
) -> list[dict[str, Any]]:
    import numpy as np

    bin_size = int(settings["bin_size"])
    rows: list[dict[str, Any]] = []
    for contig, group_data in coverage_by_group.items():
        sample_values = group_data.get(scoring_group, {})
        if not sample_values:
            continue
        length = int(ref_lengths.get(contig, 0))
        if length <= 0:
            length = max((len(values.get("sense", [])) for values in sample_values.values()), default=0)
        if length <= 0:
            continue
        sense_arrays = [np.array(values["sense"], dtype=float)[:length] for values in sample_values.values()]
        antisense_arrays = [np.array(values["antisense"], dtype=float)[:length] for values in sample_values.values()]
        if not sense_arrays:
            continue
        for start0 in range(0, length, bin_size):
            end0 = min(start0 + bin_size, length)
            bin_len = end0 - start0
            if bin_len <= 0:
                continue
            sense_areas = np.array([float(np.sum(arr[start0:end0])) for arr in sense_arrays], dtype=float)
            antisense_areas = np.array([float(np.sum(arr[start0:end0])) for arr in antisense_arrays], dtype=float)
            mean_sense_area = float(np.mean(sense_areas))
            mean_antisense_area = float(np.mean(antisense_areas))
            if mean_sense_area <= 0 and mean_antisense_area <= 0:
                continue
            sense_depth = mean_sense_area / bin_len
            antisense_depth = mean_antisense_area / bin_len
            duplex_depth = min(sense_depth, antisense_depth)
            product_score = mean_sense_area * mean_antisense_area
            total_depth = sense_depth + antisense_depth
            balance = 100.0 * (2.0 * duplex_depth / total_depth) if total_depth else 0.0
            rows.append(
                {
                    "target_bin": f"{contig}_{start0 + 1}_{end0}",
                    "source_contig": contig,
                    "scoring_group": scoring_group,
                    "start_1based": start0 + 1,
                    "end_1based": end0,
                    "bin_length_nt": bin_len,
                    "sense_area_cpm_bp": mean_sense_area,
                    "antisense_area_cpm_bp": mean_antisense_area,
                    "sense_depth_cpm": sense_depth,
                    "antisense_depth_cpm": antisense_depth,
                    "bidirectional_srna_product_score": product_score,
                    "bidirectional_srna_score_cpm": duplex_depth,
                    "total_srna_depth_cpm": total_depth,
                    "strand_balance_percent": balance,
                    "replicates": len(sense_arrays),
                }
            )
    rows.sort(key=lambda row: (float(row["bidirectional_srna_product_score"]), float(row["bidirectional_srna_score_cpm"]), float(row["total_srna_depth_cpm"])), reverse=True)
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
    return rows


def write_srna_dsrna_context_fasta(path: Path, ranked: list[dict[str, Any]], sequences: dict[str, str], context_bp: int) -> list[dict[str, Any]]:
    path.parent.mkdir(parents=True, exist_ok=True)
    context_rows: list[dict[str, Any]] = []
    with path.open("w", encoding="utf-8") as handle:
        for row in ranked:
            contig = str(row["source_contig"])
            sequence = sequences.get(contig, "")
            if not sequence:
                continue
            start0 = int(row["start_1based"]) - 1
            end0 = int(row["end_1based"])
            context_start0 = max(0, start0 - context_bp)
            context_end0 = min(len(sequence), end0 + context_bp)
            header = (
                f"rank_{int(row['rank']):03d}|{contig}:{context_start0 + 1}-{context_end0}|"
                f"scored_bin={start0 + 1}-{end0}|"
                f"bidirectional_srna_product_score={float(row['bidirectional_srna_product_score']):.6g}"
            )
            handle.write(f">{header}\n")
            context_seq = sequence[context_start0:context_end0]
            for idx in range(0, len(context_seq), 80):
                handle.write(context_seq[idx : idx + 80] + "\n")
            context_rows.append({**row, "context_start_1based": context_start0 + 1, "context_end_1based": context_end0})
    return context_rows


def plot_srna_dsrna_top_contexts(
    outdir: Path,
    ranked: list[dict[str, Any]],
    coverage_by_group: dict[str, Any],
    ref_lengths: dict[str, int],
    context_bp: int,
) -> None:
    style_settings = apply_global_plot_settings()
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import numpy as np

    plots_dir = outdir / "top_hits" / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    colors = ["#0f766e", "#7c3aed", "#ea580c", "#2563eb", "#be123c", "#16a34a"]
    for row in ranked:
        contig = str(row["source_contig"])
        length = int(ref_lengths.get(contig, 0))
        if length <= 0:
            continue
        start0 = int(row["start_1based"]) - 1
        end0 = int(row["end_1based"])
        context_start0 = max(0, start0 - context_bp)
        context_end0 = min(length, end0 + context_bp)
        x = np.arange(context_start0 + 1, context_end0 + 1)
        fig, ax = plt.subplots(figsize=(style_settings["figure_width"], style_settings["figure_height"]), constrained_layout=True)
        local_abs_max = 0.0
        group_handles: list[Any] = []
        for idx, (group, sample_values) in enumerate(sorted(coverage_by_group.get(contig, {}).items())):
            sense_arrays = [np.array(values["sense"], dtype=float)[context_start0:context_end0] for values in sample_values.values()]
            antisense_arrays = [np.array(values["antisense"], dtype=float)[context_start0:context_end0] for values in sample_values.values()]
            if not sense_arrays:
                continue
            sense = np.vstack(sense_arrays)
            antisense = np.vstack(antisense_arrays)
            sense_mean = np.mean(sense, axis=0)
            antisense_mean = np.mean(antisense, axis=0)
            sense_sem = np.std(sense, axis=0, ddof=1) / np.sqrt(sense.shape[0]) if sense.shape[0] > 1 else np.zeros_like(sense_mean)
            antisense_sem = np.std(antisense, axis=0, ddof=1) / np.sqrt(antisense.shape[0]) if antisense.shape[0] > 1 else np.zeros_like(antisense_mean)
            antisense_signed = -antisense_mean
            color = colors[idx % len(colors)]
            local_abs_max = max(
                local_abs_max,
                float(np.max(np.abs(sense_mean + sense_sem))) if sense_mean.size else 0.0,
                float(np.max(np.abs(antisense_signed - antisense_sem))) if antisense_signed.size else 0.0,
            )
            group_handles.append(Line2D([0], [0], color=color, linewidth=style_settings["line_width"], label=str(group)))
            ax.plot(x, sense_mean, color=color, linewidth=style_settings["line_width"])
            ax.fill_between(x, sense_mean - sense_sem, sense_mean + sense_sem, color=color, alpha=0.15)
            ax.plot(x, antisense_signed, color=color, linewidth=style_settings["line_width"], linestyle="--")
            ax.fill_between(x, antisense_signed - antisense_sem, antisense_signed + antisense_sem, color=color, alpha=0.10)
        ax.axhline(0, color="#303030", linewidth=max(0.5, style_settings["line_width"] * 0.45))
        ax.axvspan(start0 + 1, end0, color="#64748b", alpha=0.14)
        if local_abs_max > 0:
            ax.set_ylim(-local_abs_max * 1.15, local_abs_max * 1.15)
        rank = int(row["rank"])
        ax.set_title(f"rank {rank}: {contig}:{start0 + 1}-{end0}")
        ax.set_xlabel("Position")
        ax.set_ylabel("sRNA CPM")
        ax.grid(axis="x", color="#e5e7eb", linewidth=style_settings["grid_width"])
        style_handles = [
            Line2D([0], [0], color="#303030", linewidth=style_settings["line_width"], linestyle="-", label="sense +"),
            Line2D([0], [0], color="#303030", linewidth=style_settings["line_width"], linestyle="--", label="antisense -"),
        ]
        style_legend = ax.legend(handles=style_handles, loc="upper left", frameon=True, framealpha=0.82, fontsize=9, borderpad=0.35)
        ax.add_artist(style_legend)
        if group_handles:
            ax.legend(
                handles=group_handles,
                loc="upper right",
                frameon=True,
                framealpha=0.82,
                fontsize=9,
                borderpad=0.35,
                ncol=1 if len(group_handles) <= 4 else 2,
            )
        fig.savefig(plots_dir / f"rank_{rank:03d}_{re.sub(r'[^A-Za-z0-9_.-]+', '_', contig)[:80]}_{start0 + 1}_{end0}.png", dpi=style_settings["dpi"])
        plt.close(fig)


def parse_srna_dsrna_sam_for_bins(
    sam_path: Path,
    sample: dict[str, str],
    denominator: int,
    settings: dict[str, Any],
    ref_lengths: dict[str, int],
    job_id: str | None = None,
) -> tuple[dict[tuple[str, str, str, int, str], float], list[dict[str, Any]], dict[str, Any]]:
    import sRNA_identification as srna_core

    bin_size = int(settings["bin_size"])
    cpm = 1_000_000.0 / denominator if denominator else 0.0
    bin_counts: dict[tuple[str, str, str, int, str], float] = {}
    length_totals: dict[tuple[str, str, str, str, int, str], float] = {}
    length_counts: dict[tuple[int, str], int] = {}
    raw_mapped = 0
    retained = 0
    simple_filtered = 0
    focus_filtered = 0
    sample_id = str(sample["sample_id"])
    group = str(sample["group"])

    with sam_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            try:
                flag = int(fields[1])
                start = int(fields[3])
            except ValueError:
                continue
            if flag & 4:
                continue
            raw_mapped += 1
            if raw_mapped % 1_000_000 == 0:
                append_job_log(job_id, f"Scanned {raw_mapped:,} mapped SAM alignment(s) for {sample_id}; retained {retained:,}.")
            seq = fields[9].upper().replace("U", "T")
            length = len(seq)
            strand = "antisense" if flag & 16 else "sense"
            if settings.get("filter_simple", True) and fast_simple_srna_reason(seq):
                simple_filtered += 1
                continue
            length_counts[(length, strand)] = length_counts.get((length, strand), 0) + 1
            if settings["focus_length"] and length != int(settings["focus_length"]):
                focus_filtered += 1
                continue
            contig = fields[2]
            ref_len = srna_core.cigar_reference_length(fields[5], seq)
            if ref_len <= 0:
                continue
            contig_len = int(ref_lengths.get(contig, start + ref_len - 1))
            start0 = max(0, start - 1)
            end0 = min(contig_len, start0 + ref_len)
            if end0 <= start0:
                continue
            retained += 1
            for scope, scope_contig in (("overall", "all_contigs"), ("contig", contig)):
                key = (scope, scope_contig, group, sample_id, length, strand)
                length_totals[key] = length_totals.get(key, 0.0) + cpm
            first_bin = start0 // bin_size
            last_bin = (end0 - 1) // bin_size
            for bin_index in range(first_bin, last_bin + 1):
                bin_start = bin_index * bin_size
                bin_end = min(bin_start + bin_size, contig_len)
                overlap = min(end0, bin_end) - max(start0, bin_start)
                if overlap > 0:
                    key = (group, sample_id, contig, bin_index, strand)
                    bin_counts[key] = bin_counts.get(key, 0.0) + cpm * overlap

    length_rows = [
        {
            "scope": scope,
            "contig": contig,
            "group": group_name,
            "sample_id": sample_name,
            "length": length,
            "strand": strand,
            "cpm": value,
        }
        for (scope, contig, group_name, sample_name, length, strand), value in sorted(length_totals.items())
    ]
    stats = {
        "sample_id": sample_id,
        "group": group,
        "replicate": sample.get("replicate", ""),
        "denominator_reads": denominator,
        "mapped_alignments": retained,
        "raw_mapped_alignments": raw_mapped,
        "simple_filtered_alignments": simple_filtered,
        "focus_filtered_alignments": focus_filtered,
        "length_counts": length_counts,
    }
    return bin_counts, length_rows, stats


def score_srna_dsrna_bin_counts(
    bin_counts: dict[tuple[str, str, str, int, str], float],
    samples: list[dict[str, str]],
    ref_lengths: dict[str, int],
    settings: dict[str, Any],
    scoring_group: str,
) -> list[dict[str, Any]]:
    bin_size = int(settings["bin_size"])
    scoring_samples = [str(sample["sample_id"]) for sample in samples if str(sample["group"]) == scoring_group]
    bins = sorted({(contig, bin_index) for group, _sample, contig, bin_index, _strand in bin_counts if group == scoring_group})
    rows: list[dict[str, Any]] = []
    for contig, bin_index in bins:
        contig_len = int(ref_lengths.get(contig, 0))
        start0 = bin_index * bin_size
        end0 = min(start0 + bin_size, contig_len) if contig_len else start0 + bin_size
        bin_len = end0 - start0
        if bin_len <= 0:
            continue
        sense_values = [bin_counts.get((scoring_group, sample_id, contig, bin_index, "sense"), 0.0) for sample_id in scoring_samples]
        antisense_values = [bin_counts.get((scoring_group, sample_id, contig, bin_index, "antisense"), 0.0) for sample_id in scoring_samples]
        if not sense_values:
            continue
        mean_sense_area = sum(sense_values) / len(sense_values)
        mean_antisense_area = sum(antisense_values) / len(antisense_values)
        if mean_sense_area <= 0 and mean_antisense_area <= 0:
            continue
        sense_depth = mean_sense_area / bin_len
        antisense_depth = mean_antisense_area / bin_len
        duplex_depth = min(sense_depth, antisense_depth)
        total_depth = sense_depth + antisense_depth
        rows.append(
            {
                "target_bin": f"{contig}_{start0 + 1}_{end0}",
                "source_contig": contig,
                "scoring_group": scoring_group,
                "start_1based": start0 + 1,
                "end_1based": end0,
                "bin_length_nt": bin_len,
                "sense_area_cpm_bp": mean_sense_area,
                "antisense_area_cpm_bp": mean_antisense_area,
                "sense_depth_cpm": sense_depth,
                "antisense_depth_cpm": antisense_depth,
                "bidirectional_srna_product_score": mean_sense_area * mean_antisense_area,
                "bidirectional_srna_score_cpm": duplex_depth,
                "total_srna_depth_cpm": total_depth,
                "strand_balance_percent": 100.0 * (2.0 * duplex_depth / total_depth) if total_depth else 0.0,
                "replicates": len(scoring_samples),
            }
        )
    rows.sort(key=lambda row: (float(row["bidirectional_srna_product_score"]), float(row["bidirectional_srna_score_cpm"]), float(row["total_srna_depth_cpm"])), reverse=True)
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
    return rows


def build_srna_dsrna_context_coverage(
    sam_paths: dict[str, Path],
    samples: list[dict[str, str]],
    denominators: dict[str, int],
    settings: dict[str, Any],
    ranked: list[dict[str, Any]],
    ref_lengths: dict[str, int],
    context_bp: int,
    job_id: str | None = None,
) -> dict[int, dict[str, Any]]:
    import sRNA_identification as srna_core
    import numpy as np

    contexts: dict[int, dict[str, Any]] = {}
    by_contig: dict[str, list[dict[str, Any]]] = {}
    for row in ranked:
        contig = str(row["source_contig"])
        contig_len = int(ref_lengths.get(contig, 0))
        start0 = int(row["start_1based"]) - 1
        end0 = int(row["end_1based"])
        context_start0 = max(0, start0 - context_bp)
        context_end0 = min(contig_len, end0 + context_bp)
        rank = int(row["rank"])
        context = {
            "row": row,
            "contig": contig,
            "context_start0": context_start0,
            "context_end0": context_end0,
            "coverage": {},
        }
        contexts[rank] = context
        by_contig.setdefault(contig, []).append(context)
    for sample in samples:
        sample_id = str(sample["sample_id"])
        group = str(sample["group"])
        denominator = denominators.get(sample_id, 0)
        cpm = 1_000_000.0 / denominator if denominator else 0.0
        for context in contexts.values():
            length = int(context["context_end0"]) - int(context["context_start0"])
            context["coverage"].setdefault(group, {})[sample_id] = {
                "sense": np.zeros(length, dtype=float),
                "antisense": np.zeros(length, dtype=float),
            }
        scanned = 0
        append_job_log(job_id, f"Building top-hit coverage from SAM for {sample_id}.")
        with sam_paths[sample_id].open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("@"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 11:
                    continue
                try:
                    flag = int(fields[1])
                    start = int(fields[3])
                except ValueError:
                    continue
                if flag & 4:
                    continue
                scanned += 1
                if scanned % 1_000_000 == 0:
                    append_job_log(job_id, f"Checked {scanned:,} SAM alignment(s) for top-hit plotting in {sample_id}.")
                seq = fields[9].upper().replace("U", "T")
                if settings.get("filter_simple", True) and fast_simple_srna_reason(seq):
                    continue
                if settings["focus_length"] and len(seq) != int(settings["focus_length"]):
                    continue
                contig = fields[2]
                if contig not in by_contig:
                    continue
                strand = "antisense" if flag & 16 else "sense"
                ref_len = srna_core.cigar_reference_length(fields[5], seq)
                start0 = max(0, start - 1)
                end0 = start0 + ref_len
                for context in by_contig[contig]:
                    overlap_start0 = max(start0, int(context["context_start0"]))
                    overlap_end0 = min(end0, int(context["context_end0"]))
                    if overlap_end0 <= overlap_start0:
                        continue
                    local_start = overlap_start0 - int(context["context_start0"])
                    local_end = overlap_end0 - int(context["context_start0"])
                    context["coverage"][group][sample_id][strand][local_start:local_end] += cpm
    return contexts


def plot_srna_dsrna_streamed_contexts(outdir: Path, contexts: dict[int, dict[str, Any]]) -> None:
    style_settings = apply_global_plot_settings()
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import numpy as np

    plots_dir = outdir / "top_hits" / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    colors = ["#0f766e", "#7c3aed", "#ea580c", "#2563eb", "#be123c", "#16a34a"]
    for rank, context in sorted(contexts.items()):
        row = context["row"]
        contig = str(context["contig"])
        context_start0 = int(context["context_start0"])
        context_end0 = int(context["context_end0"])
        start0 = int(row["start_1based"]) - 1
        end0 = int(row["end_1based"])
        x = np.arange(context_start0 + 1, context_end0 + 1)
        fig, ax = plt.subplots(figsize=(style_settings["figure_width"], style_settings["figure_height"]), constrained_layout=True)
        local_abs_max = 0.0
        group_handles: list[Any] = []
        for idx, (group, sample_values) in enumerate(sorted(context["coverage"].items())):
            sense_arrays = [values["sense"] for values in sample_values.values()]
            antisense_arrays = [values["antisense"] for values in sample_values.values()]
            if not sense_arrays:
                continue
            sense = np.vstack(sense_arrays)
            antisense = np.vstack(antisense_arrays)
            sense_mean = np.mean(sense, axis=0)
            antisense_mean = np.mean(antisense, axis=0)
            sense_sem = np.std(sense, axis=0, ddof=1) / np.sqrt(sense.shape[0]) if sense.shape[0] > 1 else np.zeros_like(sense_mean)
            antisense_sem = np.std(antisense, axis=0, ddof=1) / np.sqrt(antisense.shape[0]) if antisense.shape[0] > 1 else np.zeros_like(antisense_mean)
            antisense_signed = -antisense_mean
            color = colors[idx % len(colors)]
            local_abs_max = max(
                local_abs_max,
                float(np.max(np.abs(sense_mean + sense_sem))) if sense_mean.size else 0.0,
                float(np.max(np.abs(antisense_signed - antisense_sem))) if antisense_signed.size else 0.0,
            )
            group_handles.append(Line2D([0], [0], color=color, linewidth=style_settings["line_width"], label=str(group)))
            ax.plot(x, sense_mean, color=color, linewidth=style_settings["line_width"])
            ax.fill_between(x, sense_mean - sense_sem, sense_mean + sense_sem, color=color, alpha=0.15)
            ax.plot(x, antisense_signed, color=color, linewidth=style_settings["line_width"], linestyle="--")
            ax.fill_between(x, antisense_signed - antisense_sem, antisense_signed + antisense_sem, color=color, alpha=0.10)
        ax.axhline(0, color="#303030", linewidth=max(0.5, style_settings["line_width"] * 0.45))
        ax.axvspan(start0 + 1, end0, color="#64748b", alpha=0.14)
        if local_abs_max > 0:
            ax.set_ylim(-local_abs_max * 1.15, local_abs_max * 1.15)
        ax.set_title(f"rank {rank}: {contig}:{start0 + 1}-{end0}")
        ax.set_xlabel("Position")
        ax.set_ylabel("sRNA CPM")
        ax.grid(axis="x", color="#e5e7eb", linewidth=style_settings["grid_width"])
        style_handles = [
            Line2D([0], [0], color="#303030", linewidth=style_settings["line_width"], linestyle="-", label="sense +"),
            Line2D([0], [0], color="#303030", linewidth=style_settings["line_width"], linestyle="--", label="antisense -"),
        ]
        style_legend = ax.legend(handles=style_handles, loc="upper left", frameon=True, framealpha=0.82, fontsize=9, borderpad=0.35)
        ax.add_artist(style_legend)
        if group_handles:
            ax.legend(handles=group_handles, loc="upper right", frameon=True, framealpha=0.82, fontsize=9, borderpad=0.35, ncol=1 if len(group_handles) <= 4 else 2)
        fig.savefig(plots_dir / f"rank_{rank:03d}_{re.sub(r'[^A-Za-z0-9_.-]+', '_', contig)[:80]}_{start0 + 1}_{end0}.png", dpi=style_settings["dpi"])
        plt.close(fig)


def run_srna_dsrna_identification(settings: dict[str, Any], state: dict[str, Any], job_id: str | None = None) -> dict[str, Any]:
    apply_global_plot_settings(state)
    reference_path = dsrna_reference_path(state, "srna-dsrna-identification")
    if not reference_path or not Path(reference_path).exists():
        raise ValueError("Paste or select a reference FASTA/FNA before running sRNA-based dsRNA Identification.")
    reference = Path(reference_path)
    samples = selected_srna_dsrna_samples(settings["sample_ids"], state)
    outdir = srna_dsrna_output_dir()
    for subdir in ("tables", "top_hits", "top_hits/plots", "bowtie_index", "fastq"):
        (outdir / subdir).mkdir(parents=True, exist_ok=True)

    sequences = read_fasta_sequences(reference)
    ref_lengths = {name: len(seq) for name, seq in sequences.items()}
    index_prefix = outdir / "bowtie_index" / reference.stem
    expected_index = [Path(f"{index_prefix}.{suffix}.ebwt") for suffix in ("1", "2", "3", "4", "rev.1", "rev.2")]
    if not all(path.exists() for path in expected_index):
        completed = run_tracked_command(["bowtie-build", "--threads", str(settings["threads"]), str(reference), str(index_prefix)], "Build sRNA dsRNA-identification Bowtie index", job_id)
        if completed.returncode != 0:
            raise RuntimeError(completed.stderr or completed.stdout or "Bowtie index build failed.")

    bowtie = resolve_executable("bowtie-align-s") or resolve_executable("bowtie")
    if not bowtie:
        raise ValueError("Bowtie1 was not found on PATH.")

    all_bin_counts: dict[tuple[str, str, str, int, str], float] = {}
    all_length_rows: list[dict[str, Any]] = []
    sample_stats: list[dict[str, Any]] = []
    denominators: dict[str, int] = {}
    sam_paths: dict[str, Path] = {}
    for sample in samples:
        sample_id = sample["sample_id"]
        append_job_log(job_id, f"Preparing sRNA dsRNA-identification sample {sample_id}.")
        source_fastq = Path(sample["trimmed_read1"])
        denominator = count_srna_fastq(source_fastq, job_id, sample_id)
        denominators[sample_id] = denominator
        sam_path = outdir / "tables" / f"{sample_id}.sam"
        sam_paths[sample_id] = sam_path
        command = [
            bowtie,
            "-S",
            "--no-unal",
            "-v",
            str(settings["mismatches"]),
            "-m",
            str(settings.get("max_multimappers", 50)),
            "-p",
            str(settings["threads"]),
            "-x",
            str(index_prefix),
            str(source_fastq),
            str(sam_path),
        ]
        if settings["report_all"]:
            command.insert(2, "-a")
            append_job_log(
                job_id,
                f"Report-all mapping is enabled for {sample_id}; reads with more than {settings.get('max_multimappers', 50)} placements are suppressed.",
                "warn",
            )
        else:
            append_job_log(job_id, f"Mapping sRNA sample {sample_id} to dsRNA-identification reference; suppressing >{settings.get('max_multimappers', 50)} multi-mappers.")
        completed = run_tracked_command(command, f"sRNA dsRNA-id Bowtie {sample_id}", job_id)
        if completed.returncode != 0:
            raise RuntimeError(completed.stderr or completed.stdout or f"Bowtie failed for {sample_id}")
        append_job_log(job_id, f"Streaming Bowtie SAM into bin scores for {sample_id}.")
        bin_counts, length_rows, stats = parse_srna_dsrna_sam_for_bins(
            sam_path,
            sample,
            denominator,
            settings,
            ref_lengths,
            job_id,
        )
        for key, value in bin_counts.items():
            all_bin_counts[key] = all_bin_counts.get(key, 0.0) + value
        all_length_rows.extend(length_rows)
        sample_stats.append(stats)
        append_job_log(
            job_id,
            f"Done {sample_id}: {stats['mapped_alignments']:,} retained mapped alignment(s); {stats['simple_filtered_alignments']:,} simple alignment(s) filtered.",
        )

    group_order = list(dict.fromkeys(sample["group"] for sample in samples))
    scoring_group = settings["scoring_group"] or (group_order[0] if group_order else "")
    if scoring_group not in group_order:
        scoring_group = group_order[0] if group_order else ""
    append_job_log(job_id, f"Scoring {settings['bin_size']:,} bp reference bins using group {scoring_group}.")
    scored_rows = score_srna_dsrna_bin_counts(all_bin_counts, samples, ref_lengths, settings, scoring_group)
    top_rows = [row for row in scored_rows if float(row["bidirectional_srna_product_score"]) > 0][: int(settings["top_n"])]
    if not top_rows:
        append_job_log(job_id, "No bins had both sense and antisense sRNA coverage.", "warn")
    context_rows = write_srna_dsrna_context_fasta(outdir / "top_hits" / "srna_dsrna_top_hit_contexts.fasta", top_rows, sequences, int(settings["plot_context_bp"]))
    append_job_log(job_id, "Building coverage only for top-hit context plots.")
    context_coverage = build_srna_dsrna_context_coverage(sam_paths, samples, denominators, settings, top_rows, ref_lengths, int(settings["plot_context_bp"]), job_id)
    plot_srna_dsrna_streamed_contexts(outdir, context_coverage)
    length_by_contig_rows = plot_mapped_srna_length_distribution_rows(outdir, all_length_rows)
    write_tsv(outdir / "tables" / "srna_dsrna_all_scored_bins.tsv", scored_rows)
    write_tsv(outdir / "tables" / "srna_dsrna_top_hits.tsv", top_rows)
    write_tsv(outdir / "tables" / "srna_dsrna_top_hit_contexts.tsv", context_rows)
    write_tsv(outdir / "tables" / "length_distribution.tsv", [row for row in all_length_rows if row.get("scope") == "overall"])
    write_tsv(outdir / "tables" / "mapped_length_distribution_by_contig.tsv", length_by_contig_rows)
    write_tsv(outdir / "tables" / "sample_summary.tsv", srna_sample_summary_rows(sample_stats))
    (outdir / "run_manifest.json").write_text(
        json.dumps(
            {
                "reference_fasta": str(reference),
                "sample_ids": [sample["sample_id"] for sample in samples],
                "groups": group_order,
                "scoring_group": scoring_group,
                "bin_size": settings["bin_size"],
                "plot_context_bp": settings["plot_context_bp"],
                "top_n": settings["top_n"],
                "mismatches": settings["mismatches"],
                "focus_length": settings["focus_length"],
                "filter_simple_mapped_reads": settings.get("filter_simple", True),
                "max_reportable_alignments_per_read": settings.get("max_multimappers", 50),
                "ranking": "Bins are ranked by bidirectional_srna_product_score, defined as mean total sense CPM coverage multiplied by mean total antisense CPM coverage in the bin.",
                "normalization": "Each sample uses raw retained mapped alignment counts normalized to all reads in that sample FASTQ.",
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    append_job_log(job_id, f"Wrote sRNA-based dsRNA identification outputs to {outdir}.")
    return {
        "outdir": str(outdir),
        "all_bins": str(outdir / "tables" / "srna_dsrna_all_scored_bins.tsv"),
        "top_hits": str(outdir / "tables" / "srna_dsrna_top_hits.tsv"),
        "plots": str(outdir / "top_hits" / "plots"),
        "contexts": str(outdir / "top_hits" / "srna_dsrna_top_hit_contexts.fasta"),
        "scored_bins": len(scored_rows),
        "top_bins": len(top_rows),
    }


def start_srna_dsrna_identification_job(settings: dict[str, Any]) -> str:
    job_id = create_pipeline_job("sRNA-based dsRNA identification")
    prepare_tool_output_dir("srna-dsrna-identification")

    def worker() -> None:
        try:
            result = run_srna_dsrna_identification(settings, load_state(), job_id)
        except Exception as exc:
            finish_pipeline_job(job_id, "failed", {"message": str(exc)}, f"sRNA-based dsRNA identification failed: {exc}")
            return
        finish_pipeline_job(
            job_id,
            "finished",
            result,
            f"sRNA-based dsRNA identification finished with {result['scored_bins']:,} scored bin(s) and {result['top_bins']:,} top plot(s).",
        )

    threading.Thread(target=worker, name=f"inci-srna-dsrna-id-{job_id[:8]}", daemon=True).start()
    return job_id


def summarize_srna_mapping(alignments: list[dict[str, Any]], sample_stats: list[dict[str, Any]], ref_lengths: dict[str, int]) -> tuple[list[dict[str, Any]], dict[str, Any], list[dict[str, Any]]]:
    per_sample_contig: dict[tuple[str, str, str], float] = {}
    covered: dict[tuple[str, str, str], set[int]] = {}
    coverage_by_group: dict[str, Any] = {}
    inferred_lengths = dict(ref_lengths)
    for aln in alignments:
        contig = str(aln["contig"])
        inferred_lengths[contig] = max(int(inferred_lengths.get(contig, 0)), int(aln["end"]))

    for aln in alignments:
        key = (aln["sample_id"], aln["contig"], aln["strand"])
        per_sample_contig[key] = per_sample_contig.get(key, 0.0) + float(aln["cpm"])
        covered.setdefault(key, set()).update(range(int(aln["start"]), int(aln["end"]) + 1))
        group = aln["group"]
        contig = aln["contig"]
        strand = aln["strand"]
        sample_id = aln["sample_id"]
        length = int(inferred_lengths.get(contig, int(aln["end"])))
        arr = coverage_by_group.setdefault(contig, {}).setdefault(group, {}).setdefault(sample_id, {"sense": [0.0] * length, "antisense": [0.0] * length})
        if len(arr["sense"]) < length:
            arr["sense"].extend([0.0] * (length - len(arr["sense"])))
            arr["antisense"].extend([0.0] * (length - len(arr["antisense"])))
        signed = float(aln["cpm"])
        for pos in range(max(1, int(aln["start"])), min(length, int(aln["end"])) + 1):
            arr[strand][pos - 1] += signed

    rows: list[dict[str, Any]] = []
    contigs = sorted({aln["contig"] for aln in alignments})
    for contig in contigs:
        groups = sorted({aln["group"] for aln in alignments if aln["contig"] == contig})
        for group in groups:
            sample_ids = sorted({aln["sample_id"] for aln in alignments if aln["contig"] == contig and aln["group"] == group})
            sense_values = [per_sample_contig.get((sample, contig, "sense"), 0.0) for sample in sample_ids]
            antisense_values = [per_sample_contig.get((sample, contig, "antisense"), 0.0) for sample in sample_ids]
            sense_mean, sense_sem = group_mean_sem(sense_values)
            antisense_mean, antisense_sem = group_mean_sem(antisense_values)
            rows.append(
                {
                    "contig": contig,
                    "group": group,
                    "replicates": len(sample_ids),
                    "length": inferred_lengths.get(contig, ""),
                    "sense_cpm_mean": sense_mean,
                    "sense_cpm_sem": sense_sem,
                    "antisense_cpm_mean": antisense_mean,
                    "antisense_cpm_sem": antisense_sem,
                    "total_cpm_mean": sense_mean + antisense_mean,
                    "sense_positions_covered": len(set().union(*(covered.get((sample, contig, "sense"), set()) for sample in sample_ids))) if sample_ids else 0,
                    "antisense_positions_covered": len(set().union(*(covered.get((sample, contig, "antisense"), set()) for sample in sample_ids))) if sample_ids else 0,
                }
            )
    length_rows: list[dict[str, Any]] = []
    for stats in sample_stats:
        denom = float(stats["denominator_reads"]) or 1.0
        for (length, strand), count in sorted(stats["length_counts"].items()):
            length_rows.append({"sample_id": stats["sample_id"], "group": stats["group"], "replicate": stats["replicate"], "length": length, "strand": strand, "cpm": count * 1_000_000.0 / denom})
    return rows, coverage_by_group, length_rows


def load_dsrna_overlay_coverage(
    samples: list[dict[str, str]],
    reference: Path,
    ref_lengths: dict[str, int],
    outdir: Path,
    threads: int,
    job_id: str | None,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if not samples:
        return {}, []
    if not resolve_executable("minimap2"):
        raise ValueError("minimap2 was not found on PATH; it is required for optional dsRNA-seq overlay plotting.")
    import dsRNA_identification as dsrna_core

    overlay_dir = outdir / "dsrna_overlay"
    overlay_dir.mkdir(parents=True, exist_ok=True)
    coverage_by_group: dict[str, Any] = {}
    stats_rows: list[dict[str, Any]] = []
    for sample in samples:
        sample_id = sample["sample_id"]
        append_job_log(job_id, f"Mapping dsRNA-seq overlay sample {sample_id} ({sample.get('group', 'ungrouped')} replicate {sample.get('replicate', '')}) with minimap2.")
        coverage = dsrna_core.map_sample(
            dsrna_core.SampleInput(
                name=sample_id,
                r1=Path(sample["trimmed_read1"]),
                r2=Path(sample["trimmed_read2"]),
            ),
            reference,
            ref_lengths,
            overlay_dir,
            threads=threads,
            min_mapq=10,
            direction_source="read1",
        )
        stats_rows.append(
            {
                "sample_id": sample_id,
                "group": sample.get("group", sample_id),
                "replicate": sample.get("replicate", ""),
                "library_pairs": coverage.library_pairs,
                "collapsed_pairs": coverage.stats.get("collapsed_pairs", 0),
                "sense_pairs": coverage.stats.get("sense_pairs", 0),
                "antisense_pairs": coverage.stats.get("antisense_pairs", 0),
                "discarded_pairs": coverage.stats.get("discarded_pairs", 0),
            }
        )
        append_job_log(
            job_id,
            f"Done dsRNA-seq overlay sample {sample_id}: {coverage.stats.get('collapsed_pairs', 0):,} collapsed pair(s), "
            f"{coverage.stats.get('sense_pairs', 0):,} sense, {coverage.stats.get('antisense_pairs', 0):,} antisense.",
        )
        scale = 1_000_000.0 / coverage.library_pairs if coverage.library_pairs else 0.0
        for contig, length in ref_lengths.items():
            sense = (dsrna_core.coverage_array(coverage, contig, "sense", length) * scale).tolist()
            antisense = (dsrna_core.coverage_array(coverage, contig, "antisense", length) * scale).tolist()
            coverage_by_group.setdefault(contig, {}).setdefault(sample.get("group", sample_id), {})[sample_id] = {
                "sense": sense,
                "antisense": antisense,
            }
    return coverage_by_group, stats_rows


def plot_srna_outputs(
    outdir: Path,
    summary_rows: list[dict[str, Any]],
    coverage_by_group: dict[str, Any],
    length_rows: list[dict[str, Any]],
    settings: dict[str, Any],
    dsrna_coverage_by_group: dict[str, Any] | None = None,
    sirna_annotations: list[dict[str, Any]] | None = None,
) -> None:
    style_settings = apply_global_plot_settings()
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import numpy as np

    plots_dir = outdir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    dsrna_coverage_by_group = dsrna_coverage_by_group or {}
    sirna_annotations = sirna_annotations or []
    colors = ["#0f766e", "#7c3aed", "#ea580c", "#2563eb", "#be123c", "#16a34a"]
    styles = ["-", "-.", ":"]
    row_by_contig: dict[str, list[dict[str, Any]]] = {}
    for row in summary_rows:
        if float(row.get("total_cpm_mean", 0.0)) >= settings["cpm_threshold"]:
            row_by_contig.setdefault(str(row["contig"]), []).append(row)

    for contig, rows in row_by_contig.items():
        group_data = coverage_by_group.get(contig, {})
        dsrna_group_data = dsrna_coverage_by_group.get(contig, {})
        if not group_data:
            continue
        first_arrays = [
            values
            for sample_values in group_data.values()
            for values in sample_values.values()
            if values.get("sense") and values.get("antisense")
        ]
        if not first_arrays:
            continue
        length = max(len(values["sense"]) for values in first_arrays)
        x = np.arange(1, length + 1)
        fig, ax = plt.subplots(figsize=(style_settings["figure_width"], style_settings["figure_height"]), constrained_layout=True)
        dsrna_ax = ax.twinx() if dsrna_group_data else None
        legend_handles: list[Any] = []
        axis_extrema = {"sRNA": [0.0, 0.0], "dsRNA": [0.0, 0.0]}

        def update_axis_extrema(assay: str, values: np.ndarray) -> None:
            if values.size == 0:
                return
            low, high = axis_extrema[assay]
            axis_extrema[assay] = [min(low, float(np.nanmin(values))), max(high, float(np.nanmax(values)))]

        def apply_symmetric_signed_ylim(target_ax: Any, assay: str) -> None:
            low, high = axis_extrema[assay]
            limit = max(abs(low), abs(high))
            if not math.isfinite(limit) or limit <= 0:
                limit = 1.0
            limit *= 1.08
            target_ax.set_ylim(-limit, limit)

        group_names = sorted(set(group_data) | set(dsrna_group_data))
        for idx, group in enumerate(group_names):
            color = colors[idx % len(colors)]
            style = styles[idx % len(styles)]
            mean_line_width = max(style_settings["line_width"], style_settings["line_width"] * 1.25)
            sem_line_width = max(0.8, style_settings["line_width"] * 0.65)
            for assay_name, sample_values, assay_style, alpha, width_factor in (
                ("sRNA", group_data.get(group, {}), style, 1.0, 1.0),
                ("dsRNA", dsrna_group_data.get(group, {}), (0, (8, 3)), 0.72, 0.82),
            ):
                if not sample_values:
                    continue
                target_ax = dsrna_ax if assay_name == "dsRNA" and dsrna_ax is not None else ax
                sense_arrays = [np.array(values["sense"], dtype=float) for values in sample_values.values()]
                antisense_arrays = [np.array(values["antisense"], dtype=float) for values in sample_values.values()]
                if not sense_arrays:
                    continue
                sense = np.vstack(sense_arrays)
                antisense = np.vstack(antisense_arrays)
                sense_mean = np.mean(sense, axis=0)
                antisense_mean = np.mean(antisense, axis=0)
                sense_sem = np.std(sense, axis=0, ddof=1) / np.sqrt(sense.shape[0]) if sense.shape[0] > 1 else np.zeros_like(sense_mean)
                antisense_sem = np.std(antisense, axis=0, ddof=1) / np.sqrt(antisense.shape[0]) if antisense.shape[0] > 1 else np.zeros_like(antisense_mean)
                antisense_signed = -antisense_mean
                update_axis_extrema(assay_name, sense_mean + sense_sem)
                update_axis_extrema(assay_name, sense_mean - sense_sem)
                update_axis_extrema(assay_name, antisense_signed + antisense_sem)
                update_axis_extrema(assay_name, antisense_signed - antisense_sem)
                target_ax.plot(x, sense_mean, linestyle=assay_style, color=color, linewidth=mean_line_width * width_factor, alpha=alpha)
                target_ax.plot(x, antisense_signed, linestyle=assay_style, color=color, linewidth=mean_line_width * width_factor, alpha=alpha)
                for mean_values, sem_values in ((sense_mean, sense_sem), (antisense_signed, antisense_sem)):
                    target_ax.plot(
                        x,
                        mean_values + sem_values,
                        linestyle=(0, (4, 3)),
                        color=color,
                        linewidth=sem_line_width * width_factor,
                        alpha=0.42 * alpha,
                    )
                    target_ax.plot(
                        x,
                        mean_values - sem_values,
                        linestyle=(0, (4, 3)),
                        color=color,
                        linewidth=sem_line_width * width_factor,
                        alpha=0.42 * alpha,
                    )
            legend_handles.append(Line2D([0], [0], color=color, linestyle=style, linewidth=mean_line_width, label=str(group)))
        apply_symmetric_signed_ylim(ax, "sRNA")
        ax.axhline(0, color="#303030", linewidth=max(0.5, style_settings["line_width"] * 0.45))
        if dsrna_ax is not None:
            apply_symmetric_signed_ylim(dsrna_ax, "dsRNA")
            dsrna_ax.axhline(0, color="#303030", linewidth=max(0.5, style_settings["line_width"] * 0.35), alpha=0.35)
        draw_sirna_annotation_track(
            ax,
            sirna_annotations,
            contig,
            length,
            label_size=max(6.0, float(style_settings["font_size"]) - 1.5),
        )
        ax.set_title("sRNA + dsRNA coverage" if dsrna_ax is not None else "sRNA coverage")
        ax.set_xlabel("Position")
        ax.set_ylabel("sRNA CPM")
        if dsrna_ax is not None:
            dsrna_ax.set_ylabel("dsRNA CPM")
        ax.grid(color="#e5e7eb", linewidth=style_settings["grid_width"], alpha=0.8)
        if legend_handles:
            ax.legend(handles=legend_handles, frameon=False, ncol=min(3, len(legend_handles)))
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", contig)[:120]
        fig.savefig(plots_dir / f"{safe}.positional_coverage.png", dpi=style_settings["dpi"])
        plt.close(fig)

    if length_rows:
        grouped: dict[tuple[str, int, str], list[float]] = {}
        for row in length_rows:
            grouped.setdefault((str(row["group"]), int(row["length"]), str(row["strand"])), []).append(float(row["cpm"]))
        lengths = sorted({int(row["length"]) for row in length_rows})
        groups = sorted({str(row["group"]) for row in length_rows})
        fig, ax = plt.subplots(figsize=(style_settings["figure_width"], style_settings["figure_height"]), constrained_layout=True)
        width = 0.8 / max(1, len(groups))
        base = np.arange(len(lengths), dtype=float)
        for idx, group in enumerate(groups):
            offset = (idx - (len(groups) - 1) / 2) * width
            sense_stats = [group_mean_sem(grouped.get((group, length, "sense"), [])) for length in lengths]
            antisense_stats = [group_mean_sem(grouped.get((group, length, "antisense"), [])) for length in lengths]
            sense = [mean for mean, _sem in sense_stats]
            sense_sem = [sem for _mean, sem in sense_stats]
            antisense = [mean for mean, _sem in antisense_stats]
            antisense_sem = [sem for _mean, sem in antisense_stats]
            color = colors[idx % len(colors)]
            error_style = {
                "ecolor": "#111827",
                "elinewidth": max(1.2, style_settings["line_width"] * 0.85),
                "capsize": 3.5,
                "capthick": max(1.2, style_settings["line_width"] * 0.85),
                "zorder": 5,
            }
            ax.bar(
                base + offset,
                sense,
                width=width,
                color=color,
                alpha=0.58,
                yerr=sense_sem,
                error_kw=error_style,
                label=f"{group} sense",
                zorder=3,
            )
            ax.bar(
                base + offset,
                [-v for v in antisense],
                width=width,
                color=color,
                alpha=0.25,
                yerr=antisense_sem,
                error_kw=error_style,
                label=f"{group} antisense",
                zorder=3,
            )
        ax.axhline(0, color="#303030", linewidth=max(0.5, style_settings["line_width"] * 0.45))
        ax.set_xticks(base)
        ax.set_xticklabels([str(length) for length in lengths])
        ax.set_xlabel("Length (nt)")
        ax.set_ylabel("sRNA CPM")
        ax.set_title("sRNA length distribution")
        ax.grid(axis="y", color="#e5e7eb", linewidth=style_settings["grid_width"])
        ax.legend(frameon=False, ncol=2)
        fig.savefig(plots_dir / "length_distribution.png", dpi=style_settings["dpi"])
        plt.close(fig)


def plot_mapped_srna_length_distribution_rows(outdir: Path, rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    style_settings = apply_global_plot_settings()
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    if not rows:
        return []

    length_dir = outdir / "length_distributions"
    length_dir.mkdir(parents=True, exist_ok=True)
    colors = ["#0f766e", "#7c3aed", "#ea580c", "#2563eb", "#be123c", "#16a34a"]

    def plot_scope(scope_rows: list[dict[str, Any]], label: str, filename: str) -> None:
        grouped: dict[tuple[str, str, int, str], float] = {}
        for row in scope_rows:
            key = (str(row["group"]), str(row["sample_id"]), int(row["length"]), str(row["strand"]))
            grouped[key] = grouped.get(key, 0.0) + float(row["cpm"])
        lengths = sorted({int(row["length"]) for row in scope_rows})
        groups = sorted({str(row["group"]) for row in scope_rows})
        samples_by_group = {
            group: sorted({sample_id for (grp, sample_id, _length, _strand) in grouped if grp == group})
            for group in groups
        }
        fig, ax = plt.subplots(figsize=(style_settings["figure_width"], style_settings["figure_height"]), constrained_layout=True)
        width = 0.78 / max(1, len(groups))
        base = np.arange(len(lengths), dtype=float)
        error_style = {
            "ecolor": "#111827",
            "elinewidth": max(1.1, style_settings["line_width"] * 0.8),
            "capsize": 3.0,
            "capthick": max(1.1, style_settings["line_width"] * 0.8),
            "zorder": 5,
        }
        for idx, group in enumerate(groups):
            offset = (idx - (len(groups) - 1) / 2) * width
            sample_ids = samples_by_group.get(group, [])
            sense_values: list[float] = []
            sense_sem: list[float] = []
            antisense_values: list[float] = []
            antisense_sem: list[float] = []
            for length in lengths:
                sense_reps = [grouped.get((group, sample_id, length, "sense"), 0.0) for sample_id in sample_ids]
                antisense_reps = [grouped.get((group, sample_id, length, "antisense"), 0.0) for sample_id in sample_ids]
                sense_mean, sense_err = group_mean_sem(sense_reps)
                antisense_mean, antisense_err = group_mean_sem(antisense_reps)
                sense_values.append(sense_mean)
                sense_sem.append(sense_err)
                antisense_values.append(antisense_mean)
                antisense_sem.append(antisense_err)
            color = colors[idx % len(colors)]
            ax.bar(
                base + offset,
                sense_values,
                width=width,
                color=color,
                alpha=0.62,
                yerr=sense_sem,
                error_kw=error_style,
                label=f"{group} sense",
                zorder=3,
            )
            ax.bar(
                base + offset,
                [-value for value in antisense_values],
                width=width,
                color=color,
                alpha=0.28,
                yerr=antisense_sem,
                error_kw=error_style,
                label=f"{group} antisense",
                zorder=3,
            )
        ax.axhline(0, color="#303030", linewidth=max(0.5, style_settings["line_width"] * 0.45))
        ax.set_xticks(base)
        ax.set_xticklabels([str(length) for length in lengths])
        ax.set_xlabel("Length (nt)")
        ax.set_ylabel("sRNA CPM")
        ax.set_title(label)
        ax.grid(axis="y", color="#e5e7eb", linewidth=style_settings["grid_width"])
        ax.legend(frameon=False, ncol=2)
        fig.savefig(length_dir / filename, dpi=style_settings["dpi"])
        plt.close(fig)

    plot_scope([row for row in rows if row["scope"] == "overall"], "Mapped sRNA length distribution", "overall.length_distribution.png")
    for contig in sorted({str(row["contig"]) for row in rows if row["scope"] == "contig"}):
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", contig)[:120]
        plot_scope(
            [row for row in rows if row["scope"] == "contig" and row["contig"] == contig],
            f"{contig} mapped sRNA length distribution",
            f"{safe}.length_distribution.png",
        )
    return rows


def plot_mapped_srna_length_distributions(outdir: Path, alignments: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for aln in alignments:
        cpm = float(aln["cpm"])
        base = {
            "sample_id": str(aln["sample_id"]),
            "group": str(aln["group"]),
            "length": int(aln["length"]),
            "strand": str(aln["strand"]),
            "cpm": cpm,
        }
        rows.append({**base, "scope": "overall", "contig": "all_contigs"})
        rows.append({**base, "scope": "contig", "contig": str(aln["contig"])})
    return plot_mapped_srna_length_distribution_rows(outdir, rows)


def format_cpm_token(value: float) -> str:
    if value <= 0:
        return "0"
    if value >= 1:
        return str(max(1, int(round(value))))
    from fractions import Fraction

    fraction = Fraction(value).limit_denominator(100)
    return f"{fraction.numerator}/{fraction.denominator}"


def srna_sample_summary_rows(sample_stats: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for stats in sample_stats:
        raw_mapped = int(stats.get("raw_mapped_alignments", 0))
        simple_filtered = int(stats.get("simple_filtered_alignments", 0))
        retained = int(stats.get("mapped_alignments", 0))
        rows.append(
            {
                "sample_id": stats.get("sample_id", ""),
                "group": stats.get("group", ""),
                "replicate": stats.get("replicate", ""),
                "denominator_reads": int(stats.get("denominator_reads", 0)),
                "raw_mapped_alignments": raw_mapped,
                "simple_filtered_alignments": simple_filtered,
                "focus_length_filtered_alignments": int(stats.get("focus_filtered_alignments", 0)),
                "retained_alignments": retained,
                "simple_filtered_percent_of_mapped": (simple_filtered * 100.0 / raw_mapped) if raw_mapped else 0.0,
            }
        )
    return rows


def srna_export_name(contig: str, start: int, length: int, strand: str, mismatches: int, group: str, mean_cpm: float, replicate_count: int, present_count: int) -> str:
    return f"{contig}_{start}pos_{length}len_{strand}_{mismatches}m_{format_cpm_token(mean_cpm)}CPM"


def unique_srna_table_rows(alignments: list[dict[str, Any]], settings: dict[str, Any], sample_stats: list[dict[str, Any]]) -> list[dict[str, Any]]:
    denominator_by_sample = {str(stats["sample_id"]): int(stats.get("denominator_reads", 0)) for stats in sample_stats}
    group_samples: dict[str, set[str]] = {}
    for stats in sample_stats:
        group_samples.setdefault(str(stats["group"]), set()).add(str(stats["sample_id"]))

    raw_counts: dict[tuple[str, str, str, int, int, int, str, int, str], int] = {}
    group_cpm_values: dict[tuple[str, str, int, int, str, int, str], dict[str, float]] = {}
    total_cpm_by_sample: dict[str, float] = {}
    for aln in alignments:
        if settings["export_length"] and int(aln["length"]) != int(settings["export_length"]):
            continue
        if int(aln["mismatches"]) > int(settings["mismatches"]):
            continue
        sample_id = str(aln["sample_id"])
        group = str(aln["group"])
        contig = str(aln["contig"])
        start = int(aln["start"])
        end = int(aln["end"])
        length = int(aln["length"])
        strand = str(aln["strand"])
        mismatches = int(aln["mismatches"])
        sequence = str(aln["sequence"]).upper().replace("U", "T")
        raw_key = (sample_id, group, contig, start, end, length, strand, mismatches, sequence)
        raw_counts[raw_key] = raw_counts.get(raw_key, 0) + 1
        group_key = (group, contig, start, length, strand, mismatches, sequence)
        sample_values = group_cpm_values.setdefault(group_key, {})
        sample_values[sample_id] = sample_values.get(sample_id, 0.0) + float(aln["cpm"])
        total_cpm_by_sample[sample_id] = total_cpm_by_sample.get(sample_id, 0.0) + float(aln["cpm"])

    group_summary: dict[tuple[str, str, int, int, str, int, str], dict[str, Any]] = {}
    for (group, contig, start, length, strand, mismatches, sequence), sample_values in group_cpm_values.items():
        replicate_count = max(1, len(group_samples.get(group, set(sample_values))))
        mean_cpm = sum(sample_values.values()) / replicate_count
        present_count = len(sample_values)
        group_summary[(group, contig, start, length, strand, mismatches, sequence)] = {
            "name": srna_export_name(contig, start, length, strand, mismatches, group, mean_cpm, replicate_count, present_count),
            "group_mean_CPM": mean_cpm,
            "group_replicates": replicate_count,
            "group_samples_present": present_count,
        }

    rows: list[dict[str, Any]] = []
    for (sample_id, group, contig, start, end, length, strand, mismatches, sequence), raw_count in sorted(raw_counts.items()):
        denominator = denominator_by_sample.get(sample_id, 0)
        cpm = raw_count * 1_000_000.0 / denominator if denominator else 0.0
        total_cpm = total_cpm_by_sample.get(sample_id, 0.0)
        summary = group_summary.get((group, contig, start, length, strand, mismatches, sequence), {})
        rows.append(
            {
                "sample_id": sample_id,
                "group": group,
                "siRNA_name": summary.get("name", ""),
                "sequence": sequence,
                "length": length,
                "mismatches": mismatches,
                "contig": contig,
                "strand": strand,
                "start_1based": start,
                "end_1based": end,
                "position_interval": f"{contig}:{start}-{end}({strand})",
                "raw_count": raw_count,
                "total_reads_for_CPM": denominator,
                "CPM_calculation": f"{raw_count} / {denominator} * 1,000,000",
                "CPM": cpm,
                "total_retained_siRNA_CPM_in_sample": total_cpm,
                "percent_of_total_retained_siRNA_CPM": (cpm * 100.0 / total_cpm) if total_cpm else 0.0,
                "group_mean_CPM": summary.get("group_mean_CPM", 0.0),
                "group_replicates": summary.get("group_replicates", 0),
                "group_samples_present": summary.get("group_samples_present", 0),
            }
        )
    return rows


def srna_export_records(alignments: list[dict[str, Any]], settings: dict[str, Any], sample_stats: list[dict[str, Any]]) -> list[dict[str, Any]]:
    group_samples: dict[str, set[str]] = {}
    for stats in sample_stats:
        group_samples.setdefault(str(stats["group"]), set()).add(str(stats["sample_id"]))
    aggregated: dict[tuple[str, str, int, int, str, int, str], dict[str, float]] = {}
    for aln in alignments:
        if settings["export_length"] and int(aln["length"]) != settings["export_length"]:
            continue
        if int(aln["mismatches"]) > int(settings["mismatches"]):
            continue
        key = (
            str(aln["group"]),
            str(aln["contig"]),
            int(aln["start"]),
            int(aln["length"]),
            str(aln["strand"]),
            int(aln["mismatches"]),
            str(aln["sequence"]).upper().replace("U", "T"),
        )
        sample_values = aggregated.setdefault(key, {})
        sample_id = str(aln["sample_id"])
        sample_values[sample_id] = sample_values.get(sample_id, 0.0) + float(aln["cpm"])

    rows: list[dict[str, Any]] = []
    min_cpm = float(settings["export_min_cpm"])
    for (group, contig, start, length, strand, mismatches, sequence), sample_values in aggregated.items():
        replicate_count = max(1, len(group_samples.get(group, set(sample_values))))
        mean_cpm = sum(sample_values.values()) / replicate_count
        if mean_cpm < min_cpm:
            continue
        present_count = len(sample_values)
        end = start + length - 1
        rows.append(
            {
                "export_rank": 0,
                "name": srna_export_name(contig, start, length, strand, mismatches, group, mean_cpm, replicate_count, present_count),
                "sequence": sequence,
                "length": length,
                "mismatches": mismatches,
                "contig": contig,
                "strand": strand,
                "start_1based": start,
                "end_1based": end,
                "position_interval": f"{contig}:{start}-{end}({strand})",
                "group": group,
                "group_mean_CPM": mean_cpm,
                "group_replicates": replicate_count,
                "group_samples_present": present_count,
                "sample_CPMs": ";".join(f"{sample_id}:{value:.6g}" for sample_id, value in sorted(sample_values.items())),
                "export_min_CPM_filter": min_cpm,
                "export_length_filter": int(settings["export_length"]),
            }
        )

    rows.sort(
        key=lambda row: (
            -float(row["group_mean_CPM"]),
            str(row["group"]),
            str(row["sequence"]),
            str(row["contig"]),
            int(row["start_1based"]),
            str(row["strand"]),
            int(row["mismatches"]),
        )
    )
    top_n = max(1, int(settings.get("export_top_n") or 20))
    selected = rows[:top_n]
    for idx, row in enumerate(selected, start=1):
        row["export_rank"] = idx
        row["export_top_n"] = top_n
    return selected


def write_filtered_srna_fasta(path: Path, export_records: list[dict[str, Any]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in export_records:
            handle.write(f">{row['name']}\n{row['sequence']}\n")
    return path


def start_srna_control_filtering_job(settings: dict[str, Any]) -> str:
    job_id = create_pipeline_job("sRNA control filtering")
    prepare_tool_output_dir("srna-control-filtering")

    def worker() -> None:
        try:
            result = run_srna_control_filtering(settings, job_id)
        except Exception as exc:
            finish_pipeline_job(job_id, "failed", {"message": str(exc)}, f"sRNA control filtering failed: {exc}")
            return
        finish_pipeline_job(job_id, "finished", result, f"sRNA control filtering finished with {result.get('records', 0)} sRNA record(s).")

    threading.Thread(target=worker, name=f"inci-srna-control-{job_id[:8]}", daemon=True).start()
    return job_id


def write_pasted_fasta(text: str, path: Path, default_name: str, allow_plain_sequence: bool = False) -> str:
    cleaned = text.strip()
    if not cleaned:
        return ""
    path.parent.mkdir(parents=True, exist_ok=True)
    if cleaned.startswith(">"):
        path.write_text(cleaned + "\n", encoding="utf-8")
        return str(path)
    if not allow_plain_sequence:
        raise ValueError(f"Paste FASTA records starting with > for {default_name}.")
    sequence = re.sub(r"[^A-Za-z]", "", cleaned).upper().replace("U", "T")
    if not sequence:
        raise ValueError(f"Pasted {default_name} sequence has no nucleotide bases.")
    path.write_text(f">{default_name}\n{sequence}\n", encoding="utf-8")
    return str(path)


def parse_srna_control_payload(payload: dict[str, Any]) -> dict[str, Any]:
    pasted_dir = srna_control_output_dir() / "pasted_inputs"
    unique_fasta_text = str(payload.get("unique_text", "")).strip()
    dsrna_text = str(payload.get("dsrna_reference_text", "")).strip()
    control_genome_text = str(payload.get("control_genome_text", "")).strip()
    control_srna_text = str(payload.get("control_srna_text", "")).strip()
    unique_fasta_value = write_pasted_fasta(unique_fasta_text, pasted_dir / "unique_sRNAs.fasta", "unique_sRNAs") if unique_fasta_text else str(payload.get("unique_fasta", "")).strip()
    dsrna_reference = write_pasted_fasta(dsrna_text, pasted_dir / "dsrna_reference.fasta", "pasted_dsrna_reference", True) if dsrna_text else str(payload.get("dsrna_reference_fasta", "")).strip()
    control_genome = write_pasted_fasta(control_genome_text, pasted_dir / "control_genome.fasta", "control_genome") if control_genome_text else str(payload.get("control_genome_fasta", "")).strip()
    control_srna = write_pasted_fasta(control_srna_text, pasted_dir / "control_srna.fasta", "control_srna") if control_srna_text else str(payload.get("control_srna_fasta", "")).strip()
    if not control_genome and not control_srna and DEFAULT_CONTROL_SRNA_FASTA.exists():
        control_srna = str(DEFAULT_CONTROL_SRNA_FASTA)
    unique_fasta = Path(unique_fasta_value).expanduser()
    if not unique_fasta_value or not unique_fasta.exists():
        raise ValueError(f"Unique sRNA FASTA does not exist: {unique_fasta}")
    if not dsrna_reference:
        raise ValueError("Choose a dsRNA reference FASTA for locus mapping and plots.")
    if not Path(dsrna_reference).expanduser().exists():
        raise ValueError(f"dsRNA reference FASTA does not exist: {dsrna_reference}")
    if not control_genome and not control_srna:
        raise ValueError("Choose at least one control genome or control sRNA FASTA.")
    for label, value in (("Control genome FASTA", control_genome), ("Control sRNA FASTA", control_srna)):
        if value and not Path(value).expanduser().exists():
            raise ValueError(f"{label} does not exist: {value}")
    legacy_mismatches = payload.get("mismatches", None)
    dsrna_mismatches = int(payload.get("dsrna_mismatches", legacy_mismatches if legacy_mismatches is not None else 0))
    control_mismatches = int(payload.get("control_mismatches", legacy_mismatches if legacy_mismatches is not None else 0))
    if dsrna_mismatches not in {0, 1, 2, 3}:
        raise ValueError("dsRNA mapping mismatches must be 0, 1, 2, or 3.")
    if control_mismatches not in {0, 1, 2, 3}:
        raise ValueError("Control search mismatches must be 0, 1, 2, or 3.")
    index_mode = str(payload.get("index_mode", "auto")).strip().lower()
    if index_mode not in {"auto", "fast", "lowmem"}:
        raise ValueError("Index build mode must be automatic, fast, or low memory.")
    length_filter = int(payload.get("length_filter") or 0)
    if length_filter < 0:
        raise ValueError("sRNA length filter must be zero/blank or a positive length.")
    return {
        "unique_fasta": str(unique_fasta),
        "dsrna_reference_fasta": dsrna_reference,
        "control_genome_fasta": control_genome,
        "control_srna_fasta": control_srna,
        "dsrna_mismatches": dsrna_mismatches,
        "control_mismatches": control_mismatches,
        "threads": max(1, int(payload.get("threads") or 4)),
        "index_mode": index_mode,
        "length_filter": length_filter,
        "remove_low_complexity": bool(payload.get("remove_low_complexity", True)),
        "remove_control_mappers": bool(payload.get("remove_control_mappers", True)),
        "collapse_contained": bool(payload.get("collapse_contained", False)),
    }


def read_srna_fasta_records(path: Path) -> list[dict[str, str]]:
    import sRNA_identification as srna_core

    records: list[dict[str, str]] = []
    for record in srna_core.read_fasta(path):
        seq = record.sequence.upper().replace("U", "T")
        if seq:
            records.append({"id": record.name, "sequence": seq})
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def cpm_from_srna_id(name: str) -> float:
    match = re.search(r"_([0-9]+/[0-9]+|[0-9]+(?:(?:p_?|\.)(?:[0-9]+))?)CPM(?:$|_)", name)
    if not match:
        return 0.0
    token = match.group(1)
    if "/" in token:
        numerator, denominator = token.split("/", 1)
        denominator_value = float(denominator)
        return float(numerator) / denominator_value if denominator_value else 0.0
    return float(token.replace("p_", ".").replace("p", "."))


def longest_homopolymer(sequence: str) -> tuple[str, int]:
    if not sequence:
        return "", 0
    best_base = sequence[0]
    best = 1
    current = 1
    for prev, base in zip(sequence, sequence[1:]):
        if base == prev:
            current += 1
        else:
            current = 1
        if current > best:
            best = current
            best_base = base
    return best_base, best


def dominant_kmer(sequence: str, k: int) -> tuple[str, float]:
    if len(sequence) < k:
        return "", 0.0
    counts: dict[str, int] = {}
    total = len(sequence) - k + 1
    for idx in range(total):
        kmer = sequence[idx : idx + k]
        counts[kmer] = counts.get(kmer, 0) + 1
    motif, count = max(counts.items(), key=lambda item: item[1])
    return motif, count / total if total else 0.0


def repeated_unit(sequence: str, max_unit: int = 6) -> str:
    seq = sequence.upper()
    for unit_len in range(1, min(max_unit, len(seq) // 2) + 1):
        unit = seq[:unit_len]
        if unit * (len(seq) // unit_len) + unit[: len(seq) % unit_len] == seq:
            repeats = len(seq) / unit_len
            if repeats >= 3:
                return unit
    return ""


def sequence_entropy(sequence: str) -> float:
    length = len(sequence)
    if not length:
        return 0.0
    entropy = 0.0
    for base in "ACGT":
        count = sequence.count(base)
        if count:
            p = count / length
            entropy -= p * math.log2(p)
    return entropy


def analyze_srna_complexity(sequence: str) -> dict[str, Any]:
    seq = sequence.upper().replace("U", "T")
    length = len(seq)
    counts = {base: seq.count(base) for base in "ACGT"}
    acgt_total = sum(counts.values()) or 1
    gc_fraction = (counts["G"] + counts["C"]) / acgt_total
    at_fraction = (counts["A"] + counts["T"]) / acgt_total
    dominant_base, dominant_base_count = max(counts.items(), key=lambda item: item[1])
    homo_base, homo_len = longest_homopolymer(seq)
    dinuc, dinuc_fraction = dominant_kmer(seq, 2)
    trinuc, trinuc_fraction = dominant_kmer(seq, 3)
    repeat_unit = repeated_unit(seq)
    entropy = sequence_entropy(seq)
    flags: list[str] = []
    if len([count for count in counts.values() if count]) <= 2:
        flags.append("one_or_two_base_alphabet")
    if dominant_base_count / acgt_total >= 0.80:
        flags.append("dominant_single_base")
    if homo_len >= 8 or (length and homo_len / length >= 0.60):
        flags.append("long_homopolymer")
    if dinuc_fraction >= 0.70:
        flags.append("dominant_dinucleotide")
    if trinuc_fraction >= 0.70:
        flags.append("dominant_trinucleotide")
    if repeat_unit:
        flags.append(f"tandem_repeat_{repeat_unit}")
    if entropy < 1.25:
        flags.append("low_entropy")
    if gc_fraction >= 0.85 or gc_fraction <= 0.15:
        flags.append("strong_gc_at_bias")
    if seq.count("N") / length > 0.10 if length else False:
        flags.append("many_N")
    return {
        "length": length,
        "gc_fraction": round(gc_fraction, 4),
        "at_fraction": round(at_fraction, 4),
        "unique_bases": len([base for base, count in counts.items() if count]),
        "entropy_bits": round(entropy, 4),
        "dominant_base": dominant_base,
        "dominant_base_fraction": round(dominant_base_count / acgt_total, 4),
        "longest_homopolymer_base": homo_base,
        "longest_homopolymer_length": homo_len,
        "longest_homopolymer_fraction": round(homo_len / length, 4) if length else 0.0,
        "dominant_dinucleotide": dinuc,
        "dominant_dinucleotide_fraction": round(dinuc_fraction, 4),
        "dominant_trinucleotide": trinuc,
        "dominant_trinucleotide_fraction": round(trinuc_fraction, 4),
        "tandem_repeat_unit": repeat_unit,
        "low_complexity": bool(flags),
        "complexity_flags": ";".join(flags),
    }


def parse_control_mapping_sam(sam_path: Path) -> dict[str, dict[str, Any]]:
    import sRNA_identification as srna_core

    hits: dict[str, dict[str, Any]] = {}
    if not sam_path.exists():
        return hits
    with sam_path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11 or int(fields[1]) & 4:
                continue
            qname = fields[0]
            nm = srna_core.mismatch_count(fields) or 0
            strand = "antisense" if int(fields[1]) & 16 else "sense"
            start = int(fields[3])
            ref_len = srna_core.cigar_reference_length(fields[5], fields[9])
            end = start + ref_len - 1
            record = {
                "contig": fields[2],
                "start": start,
                "end": end,
                "strand": strand,
                "mismatches": nm,
            }
            location = f"{fields[2]}:{start}-{end}:{strand}:{nm}m"
            hit = hits.setdefault(qname, {"best_mismatches": nm, "locations": [], "records": []})
            if nm < int(hit["best_mismatches"]):
                hit["best_mismatches"] = nm
                hit["locations"] = [location]
                hit["records"] = [record]
            elif nm == int(hit["best_mismatches"]) and len(hit["locations"]) < 50:
                hit["locations"].append(location)
                hit["records"].append(record)
    return hits


def bowtie_build_command(reference: Path, index_prefix: Path, settings: dict[str, Any], mode: str | None = None) -> list[str]:
    build_mode = mode or str(settings.get("index_mode", "auto"))
    command = ["bowtie-build", "--threads", str(settings["threads"])]
    if build_mode == "fast":
        command.extend(["--noauto", "--bmaxdivn", "2", "--dcv", "1024"])
    elif build_mode == "lowmem":
        command.append("--packed")
    command.extend([str(reference), str(index_prefix)])
    return command


def map_unique_srna_to_control(query_fasta: Path, reference_fasta: str, outdir: Path, label: str, settings: dict[str, Any], mismatches: int, job_id: str | None) -> dict[str, dict[str, Any]]:
    if not reference_fasta:
        return {}
    reference = Path(reference_fasta).expanduser()
    index_prefix = outdir / "bowtie_index" / label / reference.stem
    index_prefix.parent.mkdir(parents=True, exist_ok=True)
    expected_index = [Path(f"{index_prefix}.{suffix}.ebwt") for suffix in ("1", "2", "3", "4", "rev.1", "rev.2")]
    if not all(path.exists() for path in expected_index):
        mode = str(settings.get("index_mode", "auto"))
        append_job_log(job_id, f"Building {label} Bowtie index with {mode} indexing mode.")
        completed = run_tracked_command(bowtie_build_command(reference, index_prefix, settings, mode), f"Build {label} Bowtie index", job_id)
        if completed.returncode != 0 and mode == "fast":
            append_job_log(job_id, f"Fast Bowtie index build failed for {label}; retrying automatic mode.", "warn")
            completed = run_tracked_command(bowtie_build_command(reference, index_prefix, settings, "auto"), f"Build {label} Bowtie index automatic retry", job_id)
        if completed.returncode != 0:
            raise RuntimeError(completed.stderr or completed.stdout or f"Bowtie index build failed for {label}.")
    bowtie = resolve_executable("bowtie-align-s") or resolve_executable("bowtie")
    if not bowtie:
        raise ValueError("Bowtie1 was not found on PATH.")
    sam_path = outdir / "tables" / f"{label}.mapped.sam"
    completed = run_tracked_command(
        [
            bowtie,
            "-f",
            "-S",
            "--no-unal",
            "-a",
            "-v",
            str(mismatches),
            "-p",
            str(settings["threads"]),
            "-x",
            str(index_prefix),
            str(query_fasta),
            str(sam_path),
        ],
        f"sRNA control Bowtie {label}",
        job_id,
    )
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr or completed.stdout or f"Bowtie failed for {label}.")
    return parse_control_mapping_sam(sam_path)


def hit_contigs(hit: dict[str, Any] | None) -> str:
    if not hit:
        return ""
    return ";".join(sorted({str(record["contig"]) for record in hit.get("records", [])}))


def dsrna_position_summary(hit: dict[str, Any] | None) -> str:
    if not hit:
        return ""
    return ";".join(
        f"{record['contig']}:{record['start']}-{record['end']}:{record['strand']}:{record['mismatches']}m"
        for record in hit.get("records", [])
    )


def truthy_table_value(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def srna_id_strand(srna_id: str, fallback: str) -> str:
    if "_antisense_" in srna_id:
        return "antisense"
    if "_sense_" in srna_id:
        return "sense"
    return fallback


def plot_control_locus_layer(
    outdir: Path,
    dsrna_hits: dict[str, dict[str, Any]],
    control_hits: dict[str, dict[str, Any]],
    ref_lengths: dict[str, int],
    layer_name: str,
    red_label: str,
) -> list[Path]:
    apply_global_plot_settings(load_state())
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plots_dir = outdir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    records_by_contig: dict[str, list[tuple[str, dict[str, Any]]]] = {}
    for srna_id, hit in dsrna_hits.items():
        for record in hit.get("records", []):
            records_by_contig.setdefault(str(record["contig"]), []).append((srna_id, record))

    for contig, items in sorted(records_by_contig.items()):
        length = ref_lengths.get(contig, max(int(record["end"]) for _, record in items))
        fig_height = min(9.0, max(3.2, 1.2 + len(items) * 0.11))
        fig, ax = plt.subplots(figsize=(12, fig_height), constrained_layout=True)
        ax.hlines(0, 1, length, color="#111827", linewidth=1.4)
        for idx, (srna_id, record) in enumerate(items, start=1):
            display_strand = srna_id_strand(srna_id, str(record["strand"]))
            y_base = idx if display_strand == "sense" else -idx
            color = "#dc2626" if srna_id in control_hits else "#2563eb"
            ax.hlines(y_base, int(record["start"]), int(record["end"]), color=color, linewidth=3.0, alpha=0.88)
            ax.vlines([int(record["start"]), int(record["end"])], y_base - 0.18, y_base + 0.18, color=color, linewidth=1.0, alpha=0.75)
        ax.set_xlim(1, length)
        max_y = max(1, len(items))
        ax.set_ylim(-max_y - 1, max_y + 1)
        ax.axhline(0, color="#111827", linewidth=1.0)
        ax.set_xlabel("Position on dsRNA reference (nt)")
        ax.set_ylabel("Mapped sRNAs by strand")
        ax.set_title(f"{contig} sRNAs on dsRNA locus: {red_label}")
        ax.grid(axis="x", color="#e5e7eb", linewidth=0.7, alpha=0.8)
        ax.legend(
            handles=[
                Line2D([0], [0], color="#2563eb", lw=3, label="No control hit"),
                Line2D([0], [0], color="#dc2626", lw=3, label=red_label),
            ],
            frameon=False,
            loc="upper right",
        )
        safe_contig = re.sub(r"[^A-Za-z0-9_.-]+", "_", contig)[:120]
        path = plots_dir / f"{safe_contig}.{layer_name}.png"
        fig.savefig(path, dpi=220)
        plt.close(fig)
        paths.append(path)
    return paths


def plot_readable_control_overview(outdir: Path, summary_rows: list[dict[str, Any]]) -> list[Path]:
    apply_global_plot_settings(load_state())
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plots_dir = outdir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    pattern = re.compile(r"([^:;]+):(\d+)-(\d+):(sense|antisense):(\d+)m")
    rows_by_contig: dict[str, list[dict[str, Any]]] = {}
    for row in summary_rows:
        for match in pattern.finditer(str(row.get("dsrna_reference_positions", ""))):
            contig = match.group(1)
            rows_by_contig.setdefault(contig, []).append(row)

    paths: list[Path] = []
    for contig, rows in sorted(rows_by_contig.items()):
        placements: list[dict[str, Any]] = []
        for row in rows:
            hits = [
                (int(match.group(2)), int(match.group(3)), match.group(4))
                for match in pattern.finditer(str(row.get("dsrna_reference_positions", "")))
                if match.group(1) == contig
            ]
            if not hits:
                continue
            try:
                cpm = float(str(row.get("cpm_from_name", "0")).replace("p", "."))
            except ValueError:
                cpm = 0.0
            weight = (cpm if cpm > 0 else 1.0) / len(hits)
            display_strand = srna_id_strand(str(row.get("srna_id", "")), hits[0][2])
            for start, end, strand in hits:
                placements.append(
                    {
                        "start": start,
                        "end": end,
                        "strand": display_strand,
                        "cpm": weight,
                        "genome": truthy_table_value(row.get("control_genome_maps")),
                        "control_srna": truthy_table_value(row.get("control_srna_maps")),
                    }
                )
        if not placements:
            continue
        length = max(int(item["end"]) for item in placements)
        x = list(range(1, length + 1))

        def coverage(filter_fn):
            sense = [0.0] * length
            antisense = [0.0] * length
            for item in placements:
                if not filter_fn(item):
                    continue
                target = antisense if item["strand"] == "antisense" else sense
                for pos in range(max(1, int(item["start"])), min(length, int(item["end"])) + 1):
                    target[pos - 1] += float(item["cpm"])
            return sense, antisense

        total_sense, total_antisense = coverage(lambda item: True)
        genome_sense, genome_antisense = coverage(lambda item: item["genome"])
        control_sense, control_antisense = coverage(lambda item: item["control_srna"])

        fig, ax = plt.subplots(figsize=(15, 5.8), constrained_layout=True)
        ax.fill_between(x, total_sense, color="#6b7280", alpha=0.22, step="mid", label="All dsRNA-mapped sRNAs")
        ax.fill_between(x, [-value for value in total_antisense], color="#6b7280", alpha=0.22, step="mid")
        ax.plot(x, total_sense, color="#6b7280", alpha=0.35, linewidth=0.8)
        ax.plot(x, [-value for value in total_antisense], color="#6b7280", alpha=0.35, linewidth=0.8)
        if any(genome_sense) or any(genome_antisense):
            ax.fill_between(x, genome_sense, color="#dc2626", alpha=0.72, step="mid", label="Subset also mapping to control genome")
            ax.fill_between(x, [-value for value in genome_antisense], color="#dc2626", alpha=0.72, step="mid")
            ax.plot(x, genome_sense, color="#991b1b", linewidth=1.1)
            ax.plot(x, [-value for value in genome_antisense], color="#991b1b", linewidth=1.1)
        if any(control_sense) or any(control_antisense):
            ax.fill_between(x, control_sense, color="#2563eb", alpha=0.72, step="mid", label="Subset also matching control sRNA FASTA")
            ax.fill_between(x, [-value for value in control_antisense], color="#2563eb", alpha=0.72, step="mid")
            ax.plot(x, control_sense, color="#1e40af", linewidth=1.1)
            ax.plot(x, [-value for value in control_antisense], color="#1e40af", linewidth=1.1)
        else:
            ax.text(0.995, 0.04, "No control-sRNA FASTA matches", transform=ax.transAxes, ha="right", va="bottom", color="#1e40af", fontsize=10)
        ax.axhline(0, color="#111827", linewidth=0.9)
        ax.set_xlim(1, length)
        ax.set_xlabel("Position on dsRNA reference (nt)")
        ax.set_ylabel("CPM coverage\n(+ sense / - antisense)")
        ax.set_title(f"{contig}: sRNA coverage with control-hit overlays")
        ax.grid(axis="x", color="#e5e7eb", linewidth=0.7)
        ax.legend(frameon=False, ncol=2, loc="upper left", fontsize=9)
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", contig)[:120]
        path = plots_dir / f"{safe}.single_overlay_control_hits.png"
        fig.savefig(path, dpi=220)
        plt.close(fig)
        paths.append(path)
    return paths


def collapse_duplicate_and_short_contained_sRNAs(records: list[dict[str, str]]) -> tuple[list[dict[str, str]], list[dict[str, Any]], list[dict[str, Any]]]:
    from collections import Counter, defaultdict

    indexed: list[dict[str, Any]] = []
    for idx, record in enumerate(records):
        sequence = str(record["sequence"]).upper().replace("U", "T")
        indexed.append({"index": idx, "id": record["id"], "sequence": sequence, "length": len(sequence)})

    kept_by_sequence: dict[str, dict[str, Any]] = {}
    duplicate_removed: set[int] = set()
    removal_rows: list[dict[str, Any]] = []
    for record in indexed:
        existing = kept_by_sequence.get(record["sequence"])
        if existing is None:
            kept_by_sequence[record["sequence"]] = record
            continue
        duplicate_removed.add(record["index"])
        removal_rows.append(
            {
                "removed_id": record["id"],
                "removed_sequence": record["sequence"],
                "removed_length": record["length"],
                "removal_reason": "exact_duplicate",
                "parent_id": existing["id"],
                "parent_sequence": existing["sequence"],
                "parent_length": existing["length"],
                "start_in_parent_1based": 1,
                "end_in_parent_1based": record["length"],
                "additional_parent_count_same_or_longer_first_length": "",
            }
        )

    unique_records = [record for record in indexed if record["index"] not in duplicate_removed]
    by_len: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for record in unique_records:
        by_len[int(record["length"])].append(record)
    lengths = sorted(by_len)

    contained_removed: set[int] = set()
    for record in unique_records:
        sequence = record["sequence"]
        parents: list[tuple[dict[str, Any], int]] = []
        for length in lengths:
            if length <= int(record["length"]):
                continue
            for parent in by_len[length]:
                position = str(parent["sequence"]).find(sequence)
                if position != -1:
                    parents.append((parent, position + 1))
            if parents:
                break
        if parents:
            parent, position = parents[0]
            contained_removed.add(record["index"])
            removal_rows.append(
                {
                    "removed_id": record["id"],
                    "removed_sequence": sequence,
                    "removed_length": record["length"],
                    "removal_reason": "shorter_exact_substring",
                    "parent_id": parent["id"],
                    "parent_sequence": parent["sequence"],
                    "parent_length": parent["length"],
                    "start_in_parent_1based": position,
                    "end_in_parent_1based": position + int(record["length"]) - 1,
                    "additional_parent_count_same_or_longer_first_length": len(parents) - 1,
                }
            )

    removed_indexes = duplicate_removed | contained_removed
    kept = [{"id": record["id"], "sequence": record["sequence"]} for record in indexed if record["index"] not in removed_indexes]
    input_by_length = Counter(record["length"] for record in indexed)
    kept_by_length = Counter(len(record["sequence"]) for record in kept)
    removed_duplicate_by_length = Counter(row["removed_length"] for row in removal_rows if row["removal_reason"] == "exact_duplicate")
    removed_contained_by_length = Counter(row["removed_length"] for row in removal_rows if row["removal_reason"] == "shorter_exact_substring")
    summary_rows: list[dict[str, Any]] = [
        {"metric": "input_records_after_other_filters", "value": len(records), "detail": ""},
        {"metric": "kept_records_after_collapse", "value": len(kept), "detail": ""},
        {"metric": "removed_exact_duplicate_records", "value": len(duplicate_removed), "detail": ""},
        {"metric": "removed_shorter_contained_records", "value": len(contained_removed), "detail": ""},
    ]
    for length in sorted(input_by_length):
        summary_rows.append(
            {
                "metric": "length_breakdown",
                "value": length,
                "detail": f"input={input_by_length[length]}; kept={kept_by_length[length]}; duplicate_removed={removed_duplicate_by_length[length]}; shorter_contained_removed={removed_contained_by_length[length]}",
            }
        )
    return kept, removal_rows, summary_rows


def run_srna_control_filtering(settings: dict[str, Any], job_id: str | None = None) -> dict[str, Any]:
    outdir = srna_control_output_dir()
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "tables").mkdir(exist_ok=True)
    plots_dir = outdir / "plots"
    plots_dir.mkdir(exist_ok=True)
    for pattern in ("*.genome_control_hits.png", "*.control_srna_hits.png", "*.readable_control_overview.png", "*.single_overlay_control_hits.png"):
        for stale_plot in plots_dir.glob(pattern):
            stale_plot.unlink()
    unique_fasta = Path(settings["unique_fasta"]).expanduser()
    records = read_srna_fasta_records(unique_fasta)
    append_job_log(job_id, f"Loaded {len(records):,} unique sRNA FASTA records.")
    dsrna_hits = map_unique_srna_to_control(unique_fasta, settings["dsrna_reference_fasta"], outdir, "dsrna_reference", settings, int(settings["dsrna_mismatches"]), job_id)
    genome_hits = map_unique_srna_to_control(unique_fasta, settings["control_genome_fasta"], outdir, "control_genome", settings, int(settings["control_mismatches"]), job_id)
    control_srna_hits = map_unique_srna_to_control(unique_fasta, settings["control_srna_fasta"], outdir, "control_srna", settings, int(settings["control_mismatches"]), job_id)
    plot_paths: list[Path] = []

    summary_rows: list[dict[str, Any]] = []
    complexity_rows: list[dict[str, Any]] = []
    clean_records: list[dict[str, str]] = []
    length_filter = int(settings.get("length_filter") or 0)
    for record in records:
        complexity = analyze_srna_complexity(record["sequence"])
        dsrna = dsrna_hits.get(record["id"])
        genome = genome_hits.get(record["id"])
        control = control_srna_hits.get(record["id"])
        maps_to_control = bool(genome or control)
        low_complexity = bool(complexity["low_complexity"])
        passes_length_filter = not length_filter or len(record["sequence"]) == length_filter
        complexity_row = {"srna_id": record["id"], "sequence": record["sequence"], **complexity}
        complexity_rows.append(complexity_row)
        row = {
            "srna_id": record["id"],
            "sequence": record["sequence"],
            "length": len(record["sequence"]),
            "passes_length_filter": passes_length_filter,
            "requested_length_filter": length_filter or "",
            "cpm_from_name": cpm_from_srna_id(record["id"]),
            "dsrna_reference_maps": bool(dsrna),
            "dsrna_reference_best_mismatches": dsrna.get("best_mismatches", "") if dsrna else "",
            "dsrna_reference_contigs": hit_contigs(dsrna),
            "dsrna_reference_positions": dsrna_position_summary(dsrna),
            "control_genome_maps": bool(genome),
            "control_genome_best_mismatches": genome.get("best_mismatches", "") if genome else "",
            "control_genome_contigs": hit_contigs(genome),
            "control_genome_locations": ";".join(genome.get("locations", [])) if genome else "",
            "control_srna_maps": bool(control),
            "control_srna_best_mismatches": control.get("best_mismatches", "") if control else "",
            "control_srna_contigs": hit_contigs(control),
            "control_srna_locations": ";".join(control.get("locations", [])) if control else "",
            "low_complexity": low_complexity,
            "complexity_flags": complexity["complexity_flags"],
        }
        summary_rows.append(row)
        if not passes_length_filter:
            continue
        if settings["remove_low_complexity"] and low_complexity:
            continue
        if settings["remove_control_mappers"] and maps_to_control:
            continue
        clean_records.append(record)

    plot_rows = [row for row in summary_rows if truthy_table_value(row.get("passes_length_filter", True))]
    plot_paths.extend(plot_readable_control_overview(outdir, plot_rows))
    if settings.get("collapse_contained"):
        pre_collapse_count = len(clean_records)
        clean_records, collapse_removed_rows, collapse_summary_rows = collapse_duplicate_and_short_contained_sRNAs(clean_records)
        append_job_log(job_id, f"Collapsed clean FASTA records from {pre_collapse_count:,} to {len(clean_records):,} by removing exact duplicates and shorter contained reads.")
    else:
        collapse_removed_rows = []
        collapse_summary_rows = [
            {"metric": "collapse_enabled", "value": False, "detail": "Exact duplicate and shorter-contained read collapse was not selected."},
            {"metric": "input_records_after_other_filters", "value": len(clean_records), "detail": ""},
            {"metric": "kept_records_after_collapse", "value": len(clean_records), "detail": ""},
            {"metric": "removed_exact_duplicate_records", "value": 0, "detail": ""},
            {"metric": "removed_shorter_contained_records", "value": 0, "detail": ""},
        ]
    write_tsv(outdir / "tables" / "control_mapping_summary.tsv", summary_rows)
    write_tsv(outdir / "tables" / "complexity_summary.tsv", complexity_rows)
    write_tsv(outdir / "tables" / "short_contained_removed.tsv", collapse_removed_rows)
    write_tsv(outdir / "tables" / "short_contained_collapse_summary.tsv", collapse_summary_rows)
    clean_path = outdir / "filtered_unique_sRNAs.clean.fasta"
    with clean_path.open("w", encoding="utf-8") as handle:
        for record in clean_records:
            handle.write(f">{record['id']}\n{record['sequence']}\n")
    append_job_log(job_id, f"Wrote control filtering outputs to {outdir}; {len(clean_records):,} records passed the selected clean FASTA filters.")
    return {"outdir": str(outdir), "summary_tsv": str(outdir / "tables" / "control_mapping_summary.tsv"), "clean_fasta": str(clean_path), "records": len(records), "clean_records": len(clean_records), "plots": [str(path) for path in plot_paths]}


def degradome_output_dir(state: dict[str, Any] | None = None) -> Path:
    return project_output_root(state) / "degradome-analysis"


def render_degradome_content(state: dict[str, Any]) -> str:
    output_dir = tool_output_root("degradome-analysis", state)
    degradome_path, _ = effective_path(state, "degradome", "degradome-analysis")
    transcriptome_path, _ = effective_path(state, "transcriptome", "degradome-analysis")
    srna_path, _ = effective_path(state, "srna_fasta", "degradome-analysis")
    summary_tsv = output_dir / "tables" / "transcript_summary.tsv"
    transcript_count = table_distinct_value_count(summary_tsv, "transcript")
    plot_links: list[str] = []
    for plot in prioritized_plot_paths((output_dir / "plots").glob("*.png"), single_contig=transcript_count == 1):
        plot_links.append(
            f"""
            <a class="plot-card" href="/reveal?path={url_for(plot)}">
                <strong>{esc(plot.stem)}</strong>
                <img src="{output_file_url(plot)}" alt="{esc(plot.name)}">
            </a>
            """
        )
    plots = "\n".join(plot_links) if plot_links else '<div class="empty-state">Run degradome analysis to generate PNG plots.</div>'
    cleavage_tsv = output_dir / "degradome_cleavage_sites.tsv"
    target_tsv = output_dir / "tables" / "target_predictions.tsv"
    density_tsv = output_dir / "tables" / "degradome_5p_density_by_sample.tsv"
    report_html = output_dir / "degradome_analysis_report.html"
    return f"""
        <section class="hero">
            <div>
                <p class="eyebrow">Analysis tool</p>
                <h2>Degradome Analysis</h2>
                <p class="muted">Extended CleaveLand-style cleavage support with q1-ignored target pairing, CPM-normalized degradome tracks, and replicate-aware PNG plots.</p>
            </div>
            <div class="button-row top-actions">
                <a class="button-link" href="/reveal?path={url_for(output_dir)}">Open Output Folder</a>
            </div>
        </section>
        <div class="notice" id="status">Select degradome libraries above the analysis options. Samples sharing a group are treated as biological replicates.</div>
        {render_module_dataset_selectors(state, 'degradome-analysis')}
        <section class="panel">
            <div class="two-col">
                <div>
                    <p class="eyebrow">sRNAs</p>
                    <h3>Small RNA Queries</h3>
                    <textarea id="degradome-srna-text" placeholder="Paste one sRNA sequence or FASTA records"></textarea>
                    <div class="file-picker">
                        <input id="degradome-srna-fasta" value="{esc(srna_path)}" placeholder="Or select sRNA FASTA">
                        <button type="button" class="secondary" onclick="browseDegradomeInput('degradome-srna-fasta')">Browse</button>
                    </div>
                </div>
                <div>
                    <p class="eyebrow">Transcripts</p>
                    <h3>Transcript Targets</h3>
                    <textarea id="degradome-transcript-text" placeholder="Paste one transcript sequence or FASTA records"></textarea>
                    <div class="file-picker">
                        <input id="degradome-transcript-fasta" value="{esc(transcriptome_path)}" placeholder="Or select transcript FASTA">
                        <button type="button" class="secondary" onclick="browseDegradomeInput('degradome-transcript-fasta')">Browse</button>
                    </div>
                </div>
            </div>
        </section>
        <section class="panel">
            <p class="eyebrow">Analysis options</p>
            <div class="tool-grid">
                <label class="field"><span>MFE ratio cutoff</span><input id="degradome-mfe-cutoff" type="number" min="0" max="1" step="0.01" value="0.70"></label>
                <label class="field"><span>Sort target sites</span><select id="degradome-sort-by"><option value="mfe_ratio" selected>MFE ratio</option><option value="allen">Allen score</option></select></label>
                <label class="field"><span>Slice positions</span><select id="degradome-slice-mode"><option value="10" selected>q10 only</option><option value="9,10,11">q9/q10/q11</option></select></label>
                <label class="field"><span>P-value method</span><select id="degradome-pvalue-method"><option value="transcript_peak_empirical" selected>Transcript-local peak empirical</option><option value="cleaveland_category_rank">CleaveLand category/rank</option></select></label>
                <label class="field"><span>P-value cutoff</span><input id="degradome-pvalue-cutoff" type="number" min="0" max="1" step="0.001" value="0.05"></label>
                <label class="field"><span>Threads</span><input id="degradome-threads" type="number" min="1" value="4"></label>
                <label class="field"><span>Max plots</span><input id="degradome-max-plots" type="number" min="1" value="80"></label>
                <label class="field"><span>Nearby q10 marker window</span><input id="degradome-marker-window" type="number" min="1" max="3" step="1" value="3"></label>
            </div>
            <div class="button-row">
                <label class="inline-check"><input id="degradome-ignore-q1" type="checkbox" checked><span>Ignore first 5' sRNA nucleotide for pairing/MFE</span></label>
                <label class="inline-check"><input id="degradome-compact-markers" type="checkbox" checked><span>Show strongest nearby sRNA marker only</span></label>
            </div>
            {render_analysis_action('degradome-analysis', 'runDegradomeAnalysis()', 'Run Degradome Analysis')}
        </section>
        <section class="panel">
            <div class="panel-heading">
                <div><p class="eyebrow">Outputs</p><h3>Detailed Tables and PNG Plots</h3></div>
                <div class="button-row">
                    <a class="ghost" href="/reveal?path={url_for(cleavage_tsv)}">Cleavages</a>
                    <a class="ghost" href="/reveal?path={url_for(target_tsv)}">Targets</a>
                    <a class="ghost" href="/reveal?path={url_for(summary_tsv)}">Summary</a>
                    <a class="ghost" href="/reveal?path={url_for(density_tsv)}">Density</a>
                    <a class="ghost" href="/reveal?path={url_for(report_html)}">Report</a>
                </div>
            </div>
            <div class="plot-grid">{plots}</div>
        </section>
    """


def parse_degradome_payload(payload: dict[str, Any]) -> dict[str, Any]:
    samples_payload = payload.get("samples", [])
    if not isinstance(samples_payload, list) or not samples_payload:
        raise ValueError("Add at least one degradome sample.")
    samples: list[dict[str, str]] = []
    for index, item in enumerate(samples_payload, start=1):
        if not isinstance(item, dict):
            continue
        path = str(item.get("path", "")).strip()
        if not path:
            continue
        expanded = Path(path).expanduser()
        if not expanded.exists():
            raise ValueError(f"Degradome sample file does not exist: {expanded}")
        sample_id = str(item.get("sample_id", "")).strip() or expanded.stem
        samples.append(
            {
                "sample_id": sample_id,
                "group": str(item.get("group", "")).strip() or sample_id,
                "replicate": str(item.get("replicate", "")).strip() or str(index),
                "path": str(expanded),
            }
        )
    if not samples:
        raise ValueError("Add at least one readable degradome sample.")
    srna_fasta = str(payload.get("srna_fasta", "")).strip()
    transcript_fasta = str(payload.get("transcript_fasta", "")).strip()
    for label, value in (("sRNA FASTA", srna_fasta), ("Transcript FASTA", transcript_fasta)):
        if value and not Path(value).expanduser().exists():
            raise ValueError(f"{label} does not exist: {value}")
    slice_positions = tuple(int(value) for value in str(payload.get("slice_positions", "10")).split(",") if value.strip())
    if not slice_positions or any(value not in {9, 10, 11} for value in slice_positions):
        raise ValueError("Slice positions must be q10 or q9/q10/q11.")
    sort_by = str(payload.get("sort_by", "mfe_ratio"))
    if sort_by not in {"mfe_ratio", "allen"}:
        raise ValueError("Target sort mode must be MFE ratio or Allen score.")
    raw_pvalue_method = str(payload.get("pvalue_method", "transcript_peak_empirical"))
    if raw_pvalue_method == "combined_peak_mfe":
        raise ValueError("Combined peak + MFE p-value has been removed. Choose transcript-local peak empirical or CleaveLand category/rank.")
    if raw_pvalue_method not in {"transcript_peak_empirical", "peak_empirical", "peak_empirical_pvalue", "cleaveland_category_rank", "cleaveland"}:
        raise ValueError("P-value method must be transcript-local peak empirical or CleaveLand category/rank.")
    pvalue_method = normalize_pvalue_method(raw_pvalue_method)
    pvalue_cutoff = float(payload.get("pvalue_cutoff") or 0.05)
    if not 0 < pvalue_cutoff <= 1:
        raise ValueError("P-value cutoff must be >0 and <=1.")
    return {
        "srna_text": str(payload.get("srna_text", "")).strip(),
        "srna_fasta": srna_fasta,
        "transcript_text": str(payload.get("transcript_text", "")).strip(),
        "transcript_fasta": transcript_fasta,
        "samples": samples,
        "ignore_query_pos1": bool(payload.get("ignore_query_pos1", True)),
        "slice_positions": slice_positions,
        "mfe_ratio_cutoff": float(payload.get("mfe_ratio_cutoff") or 0.70),
        "sort_by": sort_by,
        "pvalue_method": pvalue_method,
        "pvalue_cutoff": pvalue_cutoff,
        "threads": max(1, int(payload.get("threads") or 4)),
        "max_transcript_plots": max(1, int(payload.get("max_transcript_plots") or 80)),
        "compact_srna_markers": bool(payload.get("compact_srna_markers", True)),
        "marker_neighborhood_nt": min(3, max(1, int(payload.get("marker_neighborhood_nt") or 3))),
    }


def start_degradome_analysis_job(settings: dict[str, Any]) -> str:
    job_id = create_pipeline_job("degradome analysis")
    prepare_tool_output_dir("degradome-analysis")

    def worker() -> None:
        try:
            apply_global_plot_settings(load_state())
            samples = tuple(
                DegradomeSample(
                    sample_id=row["sample_id"],
                    group=row["group"],
                    replicate=row["replicate"],
                    path=Path(row["path"]),
                )
                for row in settings["samples"]
            )
            config = DegradomeConfig(
                srna_text=settings["srna_text"],
                srna_fasta=Path(settings["srna_fasta"]).expanduser() if settings["srna_fasta"] else None,
                transcript_text=settings["transcript_text"],
                transcript_fasta=Path(settings["transcript_fasta"]).expanduser() if settings["transcript_fasta"] else None,
                samples=samples,
                output_dir=degradome_output_dir(load_state()),
                ignore_query_pos1=settings["ignore_query_pos1"],
                slice_positions=settings["slice_positions"],
                mfe_ratio_cutoff=settings["mfe_ratio_cutoff"],
                sort_by=settings["sort_by"],
                pvalue_method=settings["pvalue_method"],
                pvalue_cutoff=settings["pvalue_cutoff"],
                threads=settings["threads"],
                max_transcript_plots=settings["max_transcript_plots"],
                compact_srna_markers=settings["compact_srna_markers"],
                marker_neighborhood_nt=settings["marker_neighborhood_nt"],
            )

            def log(message: str, level: str = "info") -> None:
                append_job_log(job_id, message, level)

            result = run_degradome_analysis(config, log)
        except Exception as exc:
            finish_pipeline_job(job_id, "failed", {"message": str(exc)}, f"degradome analysis failed: {exc}")
            return
        finish_pipeline_job(job_id, "finished", result, f"Degradome analysis finished with {result.get('cleavage_count', 0)} cleavage-supported site(s).")

    threading.Thread(target=worker, name=f"inci-degradome-{job_id[:8]}", daemon=True).start()
    return job_id


def target_prediction_output_dir(state: dict[str, Any] | None = None) -> Path:
    return project_output_root(state) / "target-prediction"


def render_target_prediction_content(state: dict[str, Any]) -> str:
    output_dir = tool_output_root("target-prediction", state)
    transcriptome_path, _ = effective_path(state, "transcriptome", "target-prediction")
    srna_path, _ = effective_path(state, "srna_fasta", "target-prediction")
    pairs_tsv = output_dir / "target_prediction_pairs.tsv"
    tx_summary = output_dir / "tables" / "transcript_target_summary.tsv"
    srna_summary = output_dir / "tables" / "srna_target_summary.tsv"
    report_html = output_dir / "target_prediction_report.html"
    transcript_count = table_distinct_value_count(tx_summary, "transcript")
    displayed_plots = prioritized_plot_paths((output_dir / "plots").glob("*.png"), single_contig=transcript_count == 1)
    plot_html = "".join(
        f"""
        <a class="plot-card" href="/reveal?path={url_for(plot)}">
            <strong>{esc(plot.stem)}</strong>
            <img src="{output_file_url(plot)}" alt="{esc(plot.name)}">
        </a>
        """
        for plot in displayed_plots
    ) or '<div class="empty-state">Run target prediction to generate plots.</div>'
    return f"""
        <section class="hero">
            <div>
                <p class="eyebrow">Analysis tool</p>
                <h2>Target Prediction</h2>
                <p class="muted">Predict CleaveLand/GSTAr-style sRNA-transcript pairs without degradome-seq evidence.</p>
            </div>
            <div class="button-row top-actions">
                <a class="button-link" href="/reveal?path={url_for(output_dir)}">Open Output Folder</a>
            </div>
        </section>
        <div class="notice" id="status">Paste small sRNA/transcript sets directly, or select FASTA files for larger target-prediction runs.</div>
        <section class="panel">
            <div class="two-col">
                <div>
                    <p class="eyebrow">sRNAs</p>
                    <h3>Small RNA Queries</h3>
                    <textarea id="target-srna-text" placeholder="Paste one sRNA sequence or FASTA records"></textarea>
                    <div class="file-picker">
                        <input id="target-srna-fasta" value="{esc(srna_path)}" placeholder="Or select sRNA FASTA">
                        <button type="button" class="secondary" onclick="browseTargetInput('target-srna-fasta')">Browse</button>
                    </div>
                </div>
                <div>
                    <p class="eyebrow">Transcripts</p>
                    <h3>Transcript Targets</h3>
                    <textarea id="target-transcript-text" placeholder="Paste one transcript sequence or FASTA records"></textarea>
                    <div class="file-picker">
                        <input id="target-transcript-fasta" value="{esc(transcriptome_path)}" placeholder="Or select transcript FASTA">
                        <button type="button" class="secondary" onclick="browseTargetInput('target-transcript-fasta')">Browse</button>
                    </div>
                </div>
            </div>
        </section>
        <section class="panel">
            <p class="eyebrow">Prediction filters</p>
            <div class="tool-grid">
                <label class="field"><span>MFE ratio cutoff</span><input id="target-mfe-cutoff" type="number" min="0" max="1" step="0.01" value="0.70"></label>
                <label class="field"><span>Max Allen score</span><input id="target-max-allen" type="number" min="0" step="0.5" placeholder="No cutoff"></label>
                <label class="field"><span>Max mismatches</span><input id="target-max-mismatches" type="number" min="0" step="1" placeholder="No cutoff"></label>
                <label class="field"><span>Sort target sites</span><select id="target-sort-by"><option value="mfe_ratio" selected>MFE ratio</option><option value="allen">Allen score</option></select></label>
                <label class="field"><span>Transcript plot Y axis</span><select id="target-plot-metric"><option value="mfe_ratio" selected>MFE ratio</option><option value="reversed_allen">Reversed Allen score</option></select></label>
                <label class="field"><span>Maximum transcript plots</span><input id="target-max-transcript-plots" type="number" min="0" step="1" value="100"></label>
                <label class="field"><span>Threads</span><input id="target-threads" type="number" min="1" value="4"></label>
            </div>
            <div class="button-row">
                <label class="inline-check"><input id="target-ignore-q1" type="checkbox" checked><span>Ignore first 5' sRNA nucleotide for pairing/MFE</span></label>
            </div>
            {render_analysis_action('target-prediction', 'runTargetPrediction()', 'Run Target Prediction')}
        </section>
        <section class="panel">
            <div class="panel-heading">
                <div><p class="eyebrow">Outputs</p><h3>Prediction Tables and Priority Plots</h3></div>
                <div class="button-row">
                    <a class="ghost" href="/reveal?path={url_for(pairs_tsv)}">Pairs</a>
                    <a class="ghost" href="/reveal?path={url_for(tx_summary)}">Transcript Summary</a>
                    <a class="ghost" href="/reveal?path={url_for(srna_summary)}">sRNA Summary</a>
                    <a class="ghost" href="/reveal?path={url_for(output_dir / 'plots')}">Transcript Plots</a>
                    <a class="ghost" href="/reveal?path={url_for(report_html)}">Report</a>
                </div>
            </div>
            <div class="plot-grid">{plot_html}</div>
        </section>
    """


def parse_target_prediction_payload(payload: dict[str, Any]) -> dict[str, Any]:
    srna_fasta = str(payload.get("srna_fasta", "")).strip()
    transcript_fasta = str(payload.get("transcript_fasta", "")).strip()
    for label, value in (("sRNA FASTA", srna_fasta), ("Transcript FASTA", transcript_fasta)):
        if value and not Path(value).expanduser().exists():
            raise ValueError(f"{label} does not exist: {value}")
    sort_by = str(payload.get("sort_by", "mfe_ratio"))
    if sort_by not in {"mfe_ratio", "allen"}:
        raise ValueError("Target sort mode must be MFE ratio or Allen score.")
    mfe_ratio_cutoff = float(payload.get("mfe_ratio_cutoff") or 0.70)
    if not 0 < mfe_ratio_cutoff <= 1:
        raise ValueError("MFE ratio cutoff must be >0 and <=1.")
    max_allen_raw = payload.get("max_allen_score", "")
    max_allen_score = None if str(max_allen_raw).strip() == "" else float(max_allen_raw)
    if max_allen_score is not None and max_allen_score < 0:
        raise ValueError("Max Allen score must be non-negative when provided.")
    max_mismatches_raw = payload.get("max_mismatches", "")
    max_mismatches = None if str(max_mismatches_raw).strip() == "" else int(max_mismatches_raw)
    if max_mismatches is not None and max_mismatches < 0:
        raise ValueError("Max mismatches must be non-negative when provided.")
    plot_metric = str(payload.get("plot_metric", "mfe_ratio"))
    if plot_metric not in {"mfe_ratio", "reversed_allen"}:
        raise ValueError("Transcript plot Y axis must be MFE ratio or reversed Allen score.")
    return {
        "srna_text": str(payload.get("srna_text", "")).strip(),
        "srna_fasta": srna_fasta,
        "transcript_text": str(payload.get("transcript_text", "")).strip(),
        "transcript_fasta": transcript_fasta,
        "ignore_query_pos1": bool(payload.get("ignore_query_pos1", True)),
        "mfe_ratio_cutoff": mfe_ratio_cutoff,
        "max_allen_score": max_allen_score,
        "max_mismatches": max_mismatches,
        "sort_by": sort_by,
        "plot_metric": plot_metric,
        "max_transcript_plots": max(0, int(payload.get("max_transcript_plots", 100) or 0)),
        "threads": max(1, int(payload.get("threads") or 4)),
    }


def start_target_prediction_job(settings: dict[str, Any]) -> str:
    job_id = create_pipeline_job("target prediction")
    prepare_tool_output_dir("target-prediction")

    def worker() -> None:
        try:
            apply_global_plot_settings(load_state())
            config = TargetPredictionConfig(
                srna_text=settings["srna_text"],
                srna_fasta=Path(settings["srna_fasta"]).expanduser() if settings["srna_fasta"] else None,
                transcript_text=settings["transcript_text"],
                transcript_fasta=Path(settings["transcript_fasta"]).expanduser() if settings["transcript_fasta"] else None,
                output_dir=target_prediction_output_dir(load_state()),
                ignore_query_pos1=settings["ignore_query_pos1"],
                mfe_ratio_cutoff=settings["mfe_ratio_cutoff"],
                max_allen_score=settings["max_allen_score"],
                max_mismatches=settings["max_mismatches"],
                sort_by=settings["sort_by"],
                plot_metric=settings["plot_metric"],
                max_transcript_plots=settings["max_transcript_plots"],
                threads=settings["threads"],
            )

            def log(message: str, level: str = "info") -> None:
                append_job_log(job_id, message, level)

            result = run_target_prediction(config, log)
        except Exception as exc:
            finish_pipeline_job(job_id, "failed", {"message": str(exc)}, f"target prediction failed: {exc}")
            return
        finish_pipeline_job(job_id, "finished", result, f"Target prediction finished with {result.get('retained_pair_count', 0)} retained pair(s).")

    threading.Thread(target=worker, name=f"inci-target-prediction-{job_id[:8]}", daemon=True).start()
    return job_id


def run_dsrna_plotter(state: dict[str, Any], job_id: str | None = None) -> dict[str, Any]:
    apply_global_plot_settings(state)
    rnaseq_path, rnaseq_source = effective_path(state, "rnaseq", "dsrna-plotter")
    reference_path, reference_source = effective_path(state, "transcriptome", "dsrna-plotter")
    if not rnaseq_path:
        raise ValueError("Add an RNA-seq sample manifest before running dsRNA Plotter.")
    if not reference_path:
        raise ValueError("Add a template FASTA/FNA path before running dsRNA Plotter.")
    if not Path(rnaseq_path).exists():
        raise ValueError(f"RNA-seq sample manifest does not exist: {rnaseq_path}")
    if not Path(reference_path).exists():
        raise ValueError(f"Template FASTA/FNA does not exist: {reference_path}")

    output_dir = project_output_root(state) / "dsrna-plotter"
    output_dir.mkdir(parents=True, exist_ok=True)
    append_job_log(job_id, "Starting grouped paired-end dsRNA directional coverage mapping.")
    command = [
        sys.executable,
        str(APP_DIR / "dsRNA_plotter.py"),
        "--samples-csv",
        rnaseq_path,
        "--reference-fasta",
        reference_path,
        "--outdir",
        str(output_dir),
        "--threads",
        "7",
    ]
    completed = run_tracked_command(command, "dsRNA Plotter", job_id)
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout or "").strip()
        raise RuntimeError(f"dsRNA Plotter failed: {detail}")

    plots_dir = output_dir / "plots"
    plot_count = len(list(plots_dir.glob("*.png"))) if plots_dir.exists() else 0
    return {
        "records": [{}] * plot_count,
        "message": f"dsRNA Plotter finished with {plot_count} contig plot(s).",
        "input_sources": {"rnaseq": rnaseq_source, "transcriptome": reference_source},
        "outputs": {
            "outdir": str(output_dir),
            "plots": str(plots_dir),
            "sample_summary": str(output_dir / "tables" / "sample_summary.csv"),
            "contig_summary": str(output_dir / "tables" / "contig_summary.csv"),
            "manifest": str(output_dir / "run_manifest.json"),
        },
    }


def start_dsrna_plotter_job() -> str:
    """Run the potentially long dsRNA mapping outside the browser request."""
    job_id = create_pipeline_job("dsRNA Plotter")

    def worker() -> None:
        try:
            result = run_dsrna_plotter(load_state(), job_id)
        except Exception as exc:
            append_job_log(job_id, traceback.format_exc(), "error")
            finish_pipeline_job(job_id, "failed", {"message": str(exc)}, f"dsRNA Plotter failed: {exc}")
            return
        finish_pipeline_job(
            job_id,
            "finished",
            result,
            f"dsRNA Plotter finished with {len(result.get('records', []))} contig plot(s).",
        )

    threading.Thread(target=worker, name=f"inci-dsrna-plotter-{job_id[:8]}", daemon=True).start()
    return job_id


def run_dsrna_identification(state: dict[str, Any], params: dict[str, Any] | None = None) -> dict[str, Any]:
    apply_global_plot_settings(state)
    params = params or {}
    reference_path = dsrna_reference_path(state, "dsrna-identification")
    if not reference_path or not Path(reference_path).exists():
        raise ValueError("Paste or select a reference FASTA/FNA before running dsRNA Identification.")

    selected = stored_tool_samples(state, "dsrna-identification", "rnaseq")
    manifest_path = write_tool_sample_manifest(state, "dsrna-identification", "rnaseq") if selected else ""
    if not manifest_path:
        rnaseq_path, _source = effective_path(state, "rnaseq", "dsrna-identification")
        manifest_path = rnaseq_path
    if not manifest_path or not Path(manifest_path).exists():
        raise ValueError("Add paired, preprocessed RNA-seq samples before running dsRNA Identification.")
    annotation_path = optional_sirna_annotation_path(params.get("sirna_annotations_fasta", ""))
    if annotation_path:
        remember_sirna_annotation_path(state, "dsrna-identification", annotation_path)

    def integer(name: str, default: int, minimum: int) -> int:
        try:
            value = int(params.get(name, default))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{name.replace('_', ' ').capitalize()} must be an integer.") from exc
        if value < minimum:
            raise ValueError(f"{name.replace('_', ' ').capitalize()} must be at least {minimum}.")
        return value

    bin_size = integer("bin_size", 250, 25)
    plot_context_bp = integer("plot_context_bp", 500, 0)
    top_n = integer("top_n", 20, 1)
    scoring_group = str(params.get("scoring_group", "")).strip()
    included_groups = {
        str(sample.get("group", ""))
        for sample in selected
        if sample.get("included", True) and sample.get("group")
    }
    if scoring_group not in included_groups:
        scoring_group = ""
    output_dir = project_output_root(state) / "dsrna-identification" / "genome_loci"
    command = [
        sys.executable,
        str(APP_DIR / "scripts" / "run_genome_dsrna_loci.py"),
        "--reference-fasta", reference_path,
        "--samples-csv", manifest_path,
        "--outdir", str(output_dir),
        "--bin-size", str(bin_size),
        "--plot-context-bp", str(plot_context_bp),
        "--top-n", str(top_n),
        "--summary-top-n", "50",
        "--max-reads", "1000000",
        "--threads", "7",
        "--min-mapq", "0",
        "--smooth-span", "25",
    ]
    if scoring_group:
        command.extend(["--scoring-group", scoring_group])
    if annotation_path:
        command.extend(["--sirna-annotations-fasta", annotation_path])
    completed = run_tracked_command(command, "dsRNA locus identification")
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout or "").strip()
        raise RuntimeError(f"dsRNA Identification failed: {detail}")
    all_bins = output_dir / "tables" / "dsrna_loci_all_scored_bins.csv"
    row_count = max(0, sum(1 for _ in all_bins.open(encoding="utf-8")) - 1) if all_bins.exists() else 0
    return {
        "records": [{}] * row_count,
        "message": f"dsRNA Identification finished with {row_count:,} scored bins and {top_n} directional top-hit plots.",
        "outputs": {
            "outdir": str(output_dir),
            "all_bins": str(all_bins),
            "top_hits": str(output_dir / "tables" / "dsrna_loci_top_hits.csv"),
            "summary_plots": str(output_dir / "phase1" / "summary_plots"),
            "plots": str(output_dir / "top_hits" / "plots"),
            "sirna_annotations": str(output_dir / "top_hits" / "tables" / "siRNA_annotation_matches.tsv") if annotation_path else "",
            "manifest": str(output_dir / "run_manifest.json"),
        },
    }


def run_fasta_deduplication(state: dict[str, Any], params: dict[str, Any] | None = None) -> dict[str, Any]:
    params = params or {}
    input_value = str(params.get("input_fasta", "")).strip() or fasta_dedup_input_path(state)
    if not input_value:
        raise ValueError("Choose a nucleotide FASTA file before running FASTA Deduplication.")
    input_fasta = Path(input_value).expanduser()
    if not input_fasta.exists() or not input_fasta.is_file():
        raise ValueError(f"Input FASTA does not exist: {input_fasta}")
    if not resolve_executable("cd-hit-est"):
        raise ValueError("cd-hit-est was not found. Install CD-HIT from Settings before running FASTA Deduplication.")

    try:
        identity_percent = float(params.get("identity_percent", 95))
    except (TypeError, ValueError) as exc:
        raise ValueError("Identity must be a number between 75 and 100.") from exc
    if not 75 <= identity_percent <= 100:
        raise ValueError("Identity must be between 75 and 100.")

    def integer(name: str, default: int, minimum: int) -> int:
        try:
            value = int(params.get(name, default))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{name.replace('_', ' ').capitalize()} must be an integer.") from exc
        if value < minimum:
            raise ValueError(f"{name.replace('_', ' ').capitalize()} must be at least {minimum}.")
        return value

    threads = integer("threads", 7, 1)
    memory_mb = integer("memory_mb", 0, 0)
    word_size = str(params.get("word_size", "auto")).strip().lower()
    if word_size != "auto":
        try:
            parsed_word_size = int(word_size)
        except ValueError as exc:
            raise ValueError("Word size must be Auto or an integer from 4 to 11.") from exc
        if not 4 <= parsed_word_size <= 11:
            raise ValueError("Word size must be between 4 and 11.")
        word_size = str(parsed_word_size)

    strand_mode = str(params.get("strand_mode", "both")).strip().lower()
    if strand_mode not in {"both", "same"}:
        raise ValueError("Strand mode must be both or same.")

    state.setdefault("module_paths", {}).setdefault("fasta-deduplication", {})["transcriptome"] = str(input_fasta)
    save_state(state)

    output_dir = project_output_root(state) / "fasta-deduplication"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_fasta = output_dir / "deduplicated.fasta"
    clusters_csv = output_dir / "clusters.csv"
    log_file = output_dir / "cd-hit-est.log"
    manifest_json = output_dir / "run_manifest.json"
    command = [
        sys.executable,
        str(APP_DIR / "scripts" / "deduplicate_fasta_cdhit_est.py"),
        "--input-fasta",
        str(input_fasta),
        "--output-fasta",
        str(output_fasta),
        "--clusters-csv",
        str(clusters_csv),
        "--log-file",
        str(log_file),
        "--manifest-json",
        str(manifest_json),
        "--identity",
        f"{identity_percent / 100:.5f}",
        "--threads",
        str(threads),
        "--memory-mb",
        str(memory_mb),
        "--force",
    ]
    if word_size != "auto":
        command.extend(["--word-size", word_size])
    if strand_mode == "same":
        command.append("--same-strand-only")

    completed = run_tracked_command(command, "FASTA deduplication with CD-HIT-EST")
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout or "").strip()
        raise RuntimeError(f"FASTA Deduplication failed: {detail}")

    manifest = json.loads(manifest_json.read_text(encoding="utf-8")) if manifest_json.exists() else {}
    representative_count = int(manifest.get("representative_records", 0) or 0)
    removed_count = int(manifest.get("removed_records", 0) or 0)
    return {
        "records": [{}] * representative_count,
        "message": f"FASTA Deduplication finished with {representative_count:,} representative record(s); removed {removed_count:,} redundant record(s).",
        "outputs": {
            "outdir": str(output_dir),
            "deduplicated_fasta": str(output_fasta),
            "clusters_csv": str(clusters_csv),
            "cluster_file": str(Path(str(output_fasta) + ".clstr")),
            "log": str(log_file),
            "manifest": str(manifest_json),
        },
    }


def run_module(module_key: str, state: dict[str, Any], params: dict[str, Any] | None = None) -> dict[str, Any]:
    if module_key == "dsrna-identification":
        return run_dsrna_identification(state, params)
    if module_key == "dsrna-plotter":
        return run_dsrna_plotter(state)
    if module_key == "fasta-deduplication":
        return run_fasta_deduplication(state, params)
    raise ValueError("This pipeline section does not have a runnable action yet.")


def start_module_job(module_key: str, params: dict[str, Any]) -> str:
    state = load_state()
    job_id = create_pipeline_job(module_for_key(module_key).title)

    def worker() -> None:
        try:
            result = run_module(module_key, state, params)
        except BaseException as exc:
            append_job_log(job_id, traceback.format_exc(), 'error')
            finish_pipeline_job(job_id, 'failed', {'message': str(exc)}, f'Analysis failed: {exc}')
            return
        finish_pipeline_job(job_id, 'finished', result, result.get('message', 'Analysis completed.'))

    threading.Thread(target=worker, name=f'inci-module-{job_id[:8]}', daemon=True).start()
    return job_id


def dsrna_reference_path(state: dict[str, Any], module_key: str) -> str:
    module_paths = state.get("module_paths", {}).get(module_key, {})
    if not isinstance(module_paths, dict):
        return ""
    return str(module_paths.get("transcriptome", "")).strip()


def sirna_annotation_path(state: dict[str, Any], module_key: str) -> str:
    module_paths = state.get("module_paths", {}).get(module_key, {})
    if not isinstance(module_paths, dict):
        return ""
    return str(module_paths.get("sirna_annotations", "")).strip()


def optional_sirna_annotation_path(value: Any) -> str:
    selected = str(value or "").strip()
    if not selected:
        return ""
    path = Path(selected).expanduser()
    if not path.exists() or not path.is_file():
        raise ValueError(f"siRNA annotation FASTA does not exist: {path}")
    return str(path)


def remember_sirna_annotation_path(state: dict[str, Any], module_key: str, path: str) -> None:
    state.setdefault("module_paths", {}).setdefault(module_key, {})["sirna_annotations"] = path
    save_state(state)


def save_dsrna_reference(state: dict[str, Any], module_key: str, sequence_text: str, fasta_path: str) -> str:
    if module_key not in {"dsrna-identification", "dsrna-plotter", "srna-dsrna-identification"}:
        raise ValueError("Reference input is only available for dsRNA analysis tools.")
    selected = fasta_path.strip()
    if selected:
        path = Path(selected).expanduser()
        if not path.exists() or not path.is_file():
            raise ValueError(f"Reference FASTA does not exist: {path}")
        reference_path = path
    else:
        text = sequence_text.strip()
        if not text:
            existing = dsrna_reference_path(state, module_key)
            if existing and Path(existing).exists():
                return existing
            raise ValueError("Paste a reference sequence or choose a FASTA file.")
        references_dir = project_output_root(state) / "references"
        references_dir.mkdir(parents=True, exist_ok=True)
        reference_path = references_dir / f"{module_key}_reference.fasta"
        if text.startswith(">"):
            fasta_text = text
        else:
            sequence = re.sub(r"[^A-Za-z]", "", text).upper().replace("U", "T")
            if not sequence or re.search(r"[^ACGTN]", sequence):
                raise ValueError("The pasted reference contains unsupported sequence characters.")
            fasta_text = f">pasted_reference\n{sequence}"
        reference_path.write_text(fasta_text.rstrip() + "\n", encoding="utf-8")
    state.setdefault("module_paths", {}).setdefault(module_key, {})["transcriptome"] = str(reference_path)
    save_state(state)
    return str(reference_path)


def render_dsrna_reference_selector(state: dict[str, Any], module_key: str) -> str:
    current = dsrna_reference_path(state, module_key)
    return f"""
        <section class="panel">
            <div class="panel-heading">
                <div>
                    <p class="eyebrow">Reference</p>
                    <h3>Reference Sequence</h3>
                </div>
                <button type="button" onclick="saveDsrnaReference('{esc(module_key)}')">Save Reference</button>
            </div>
            <div class="two-col">
                <label class="field">
                    <span>Paste sequence or FASTA</span>
                    <textarea id="dsrna-reference-text-{esc(module_key)}" placeholder="Paste a nucleotide sequence or FASTA records"></textarea>
                </label>
                <label class="field">
                    <span>FASTA file</span>
                    <div class="file-picker">
                        <input id="dsrna-reference-file-{esc(module_key)}" value="{esc(current)}" placeholder="Select a FASTA file">
                        <button type="button" class="secondary" onclick="browseDsrnaReference('{esc(module_key)}')">Browse</button>
                    </div>
                </label>
            </div>
        </section>
    """


def render_sirna_annotation_selector(state: dict[str, Any], module_key: str, input_id: str) -> str:
    current = sirna_annotation_path(state, module_key)
    return f"""
        <section class="panel">
            <div class="panel-heading">
                <div>
                    <p class="eyebrow">Optional annotations</p>
                    <h3>siRNA Annotation FASTA</h3>
                    <p class="muted">Exact matches are labelled on positional plots in their correct forward or reverse orientation. This does not change read mapping or scores.</p>
                </div>
            </div>
            <label class="field">
                <span>siRNA FASTA</span>
                <div class="file-picker">
                    <input id="{esc(input_id)}" value="{esc(current)}" placeholder="Optional FASTA containing labelled siRNAs">
                    <button type="button" class="secondary" onclick="browseSirnaAnnotations('{esc(module_key)}', '{esc(input_id)}')">Browse</button>
                </div>
            </label>
        </section>
    """


def fasta_dedup_input_path(state: dict[str, Any]) -> str:
    module_paths = state.get("module_paths", {}).get("fasta-deduplication", {})
    if isinstance(module_paths, dict):
        selected = str(module_paths.get("transcriptome", "")).strip()
        if selected:
            return selected
    path, _source = effective_path(state, "transcriptome", "fasta-deduplication")
    return path


def render_fasta_deduplication_controls(state: dict[str, Any], output_dir: Path) -> str:
    current = fasta_dedup_input_path(state)
    ready = bool(resolve_executable("cd-hit-est"))
    dependency_notice = (
        ""
        if ready
        else '<p class="muted">CD-HIT-EST is not available yet. Install CD-HIT from Settings before running this tool.</p>'
    )
    return f"""
        <section class="panel">
            <div class="panel-heading">
                <div>
                    <p class="eyebrow">Input</p>
                    <h3>Nucleotide FASTA</h3>
                    <p>Select transcriptome, contig, or candidate-region FASTA records to cluster.</p>
                </div>
            </div>
            <label class="field">
                <span>Input FASTA</span>
                <div class="file-picker">
                    <input id="fasta-dedup-input" value="{esc(current)}" placeholder="Select a nucleotide FASTA file">
                    <button type="button" class="secondary" onclick="browseFastaDedupInput()">Browse</button>
                </div>
            </label>
        </section>
        <section class="panel">
            <div class="panel-heading">
                <div>
                    <p class="eyebrow">CD-HIT-EST</p>
                    <h3>Deduplication Settings</h3>
                    <p>Cluster near-identical nucleotide records and write representative sequences plus a cluster membership CSV.</p>
                    {dependency_notice}
                </div>
            </div>
            <div class="settings-grid">
                <label class="field"><span>Identity (%)</span><input id="fasta-dedup-identity" type="number" min="75" max="100" step="0.1" value="95"></label>
                <label class="field"><span>Word size</span><select id="fasta-dedup-word-size"><option value="auto">Auto</option><option value="11">11</option><option value="10">10</option><option value="9">9</option><option value="8">8</option><option value="7">7</option><option value="6">6</option><option value="5">5</option><option value="4">4</option></select></label>
                <label class="field"><span>Threads</span><input id="fasta-dedup-threads" type="number" min="1" step="1" value="7"></label>
                <label class="field"><span>Memory (MB)</span><input id="fasta-dedup-memory" type="number" min="0" step="256" value="0"></label>
                <label class="field"><span>Strand mode</span><select id="fasta-dedup-strand"><option value="both">Both strands</option><option value="same">Same strand only</option></select></label>
            </div>
            {render_analysis_action('fasta-deduplication', "runModule('fasta-deduplication')", 'Deduplicate FASTA', disabled=not ready)}
        </section>
        <section class="panel">
            <div class="panel-heading">
                <div><p class="eyebrow">Outputs</p><h3>CD-HIT Results</h3></div>
                <a class="button-link" href="/reveal?path={url_for(output_dir)}">Open Output Folder</a>
            </div>
            <div class="path-row">
                <span class="badge">{'OK' if (output_dir / 'deduplicated.fasta').exists() else '--'}</span>
                <div><strong>Deduplicated FASTA</strong><small>{esc(output_dir / 'deduplicated.fasta')}</small></div>
                <a href="/reveal?path={url_for(output_dir / 'deduplicated.fasta')}">Reveal</a>
            </div>
            <div class="path-row">
                <span class="badge">{'OK' if (output_dir / 'clusters.csv').exists() else '--'}</span>
                <div><strong>Cluster Summary CSV</strong><small>{esc(output_dir / 'clusters.csv')}</small></div>
                <a href="/reveal?path={url_for(output_dir / 'clusters.csv')}">Reveal</a>
            </div>
            <div class="path-row">
                <span class="badge">{'OK' if (output_dir / 'run_manifest.json').exists() else '--'}</span>
                <div><strong>Run Manifest</strong><small>{esc(output_dir / 'run_manifest.json')}</small></div>
                <a href="/reveal?path={url_for(output_dir / 'run_manifest.json')}">Reveal</a>
            </div>
        </section>
    """


def render_dsrna_identification_results(output_dir: Path) -> str:
    summary_dir = output_dir / "phase1" / "summary_plots"
    table_path = output_dir / "tables" / "dsrna_loci_all_scored_bins.csv"
    annotation_tsv = output_dir / "top_hits" / "tables" / "siRNA_annotation_matches.tsv"
    contig_count = table_distinct_value_count(table_path, "source_contig")
    candidates = [*summary_dir.glob("*.png"), *(output_dir / "top_hits" / "plots").glob("*.png")]
    cards = []
    for plot in prioritized_plot_paths(candidates, single_contig=contig_count == 1):
        cards.append(
            f'<a class="plot-card" href="/reveal?path={url_for(plot)}"><strong>{esc(plot.stem)}</strong><img src="{output_file_url(plot)}" alt="{esc(plot.stem)}"></a>'
        )
    plots = "".join(cards) or '<div class="empty-state">Run the analysis to create the locus summaries.</div>'
    preview = '<div class="empty-state">No scored-bin table yet.</div>'
    if table_path.exists():
        try:
            rows = list(csv.DictReader(table_path.open(encoding="utf-8")))[:20]
            columns = ["rank", "source_contig", "start_1based", "end_1based", "strand_balance_percent", "duplex_depth_cpm", "coverage_score_percent", "total_dsrna_score_percent"]
            header = "".join(f"<th>{esc(column.replace('_', ' '))}</th>" for column in columns)
            body = "".join(
                "<tr>" + "".join(f"<td>{esc(format_preview_value(row.get(column, '')))}</td>" for column in columns) + "</tr>"
                for row in rows
            )
            preview = f'<div class="table-scroll"><table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table></div>'
        except (OSError, csv.Error):
            preview = '<div class="empty-state">The scored-bin table could not be read.</div>'
    return f"""
        <section class="panel">
            <div class="panel-heading">
                <div><p class="eyebrow">First phase</p><h3>Priority Plots</h3></div>
                <a class="button-link" href="/reveal?path={url_for(summary_dir)}">Open Summary Folder</a>
            </div>
            <div class="plot-grid">{plots}</div>
        </section>
        <section class="panel">
            <div class="panel-heading">
                <div><p class="eyebrow">All scored loci</p><h3>dsRNA Locus Table</h3></div>
                <div class="button-row"><a class="button-link" href="/reveal?path={url_for(table_path)}">Open CSV</a><a class="button-link" href="/reveal?path={url_for(annotation_tsv)}">siRNA Annotations</a><a class="button-link" href="/reveal?path={url_for(output_dir / 'top_hits' / 'plots')}">Open Top-Hit Plots</a></div>
            </div>
            {preview}
        </section>
    """


def render_dsrna_identification_content(state: dict[str, Any]) -> str:
    output_dir = tool_output_root("dsrna-identification", state) / "genome_loci"
    selected = stored_tool_samples(state, "dsrna-identification", "rnaseq")
    groups = list(dict.fromkeys(sample["group"] for sample in selected if sample.get("included", True) and sample.get("group")))
    group_options = "".join(f'<option value="{esc(group)}">{esc(group)}</option>' for group in groups)
    no_group = '<option value="">First selected group</option>'
    return f"""
        <section class="hero">
            <div>
                <p class="eyebrow">Analysis tool</p>
                <h2>dsRNA Identification</h2>
                <p class="muted">Rank locally balanced, bidirectionally covered loci from paired-end RNA-seq.</p>
            </div>
            <div class="button-row top-actions">
                <a class="button-link" href="/reveal?path={url_for(output_dir)}">Open Output Folder</a>
            </div>
        </section>
        <div class="notice" id="status">Select paired preprocessed samples, assign biological groups, and choose a genomic or multi-contig FASTA reference.</div>
        {render_module_dataset_selectors(state, 'dsrna-identification')}
        {render_dsrna_reference_selector(state, 'dsrna-identification')}
        {render_sirna_annotation_selector(state, 'dsrna-identification', 'dsrna-sirna-annotations-fasta')}
        <section class="panel">
            <div class="panel-heading">
                <div><p class="eyebrow">Locus screening</p><h3>First-Phase Settings</h3></div>
            </div>
            <div class="settings-grid">
                <label class="field"><span>Scoring group</span><select id="dsrna-score-group">{no_group}{group_options}</select></label>
                <label class="field"><span>Bin length (bp)</span><input id="dsrna-bin-size" type="number" min="25" step="25" value="250"></label>
                <label class="field"><span>Context around hit (bp)</span><input id="dsrna-context-bp" type="number" min="0" step="25" value="500"></label>
                <label class="field"><span>Top hits to plot</span><input id="dsrna-top-n" type="number" min="1" step="1" value="20"></label>
            </div>
            {render_analysis_action('dsrna-identification', "runModule('dsrna-identification')", 'Identify dsRNA Loci')}
        </section>
        {render_dsrna_identification_results(output_dir)}
    """


def render_srna_dsrna_identification_results(output_dir: Path) -> str:
    table_path = output_dir / "tables" / "srna_dsrna_all_scored_bins.tsv"
    top_path = output_dir / "tables" / "srna_dsrna_top_hits.tsv"
    context_fasta = output_dir / "top_hits" / "srna_dsrna_top_hit_contexts.fasta"
    plots_dir = output_dir / "top_hits" / "plots"
    preview = '<div class="empty-state">No scored sRNA-bin table yet.</div>'
    if table_path.exists():
        try:
            with table_path.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))[:20]
            columns = [
                "rank",
                "source_contig",
                "start_1based",
                "end_1based",
                "bidirectional_srna_product_score",
                "sense_depth_cpm",
                "antisense_depth_cpm",
                "strand_balance_percent",
            ]
            header = "".join(f"<th>{esc(column.replace('_', ' '))}</th>" for column in columns)
            body = "".join(
                "<tr>" + "".join(f"<td>{esc(format_preview_value(row.get(column, '')))}</td>" for column in columns) + "</tr>"
                for row in rows
            )
            preview = f'<div class="table-scroll"><table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table></div>'
        except (OSError, csv.Error):
            preview = '<div class="empty-state">The scored-bin table could not be read.</div>'
    plot_cards = []
    if plots_dir.exists():
        contig_count = table_distinct_value_count(table_path, "source_contig")
        for plot in prioritized_plot_paths(plots_dir.glob("*.png"), single_contig=contig_count == 1):
            plot_cards.append(f'<a class="plot-card" href="/reveal?path={url_for(plot)}"><strong>{esc(plot.stem)}</strong><img src="{output_file_url(plot)}" alt="{esc(plot.stem)}"></a>')
    plots = "".join(plot_cards) or '<div class="empty-state">Run the analysis to create top-hit context plots.</div>'
    return f"""
        <section class="panel">
            <div class="panel-heading">
                <div><p class="eyebrow">Top sRNA-supported bins</p><h3>Ranked Bin Table</h3></div>
                <div class="button-row">
                    <a class="button-link" href="/reveal?path={url_for(table_path)}">Open All Bins</a>
                    <a class="button-link" href="/reveal?path={url_for(top_path)}">Open Top Hits</a>
                    <a class="button-link" href="/reveal?path={url_for(context_fasta)}">Open Context FASTA</a>
                </div>
            </div>
            {preview}
        </section>
        <section class="panel">
            <div class="panel-heading">
                <div><p class="eyebrow">Top-hit contexts</p><h3>sRNA Directional Coverage</h3></div>
                <a class="button-link" href="/reveal?path={url_for(plots_dir)}">Open Plot Folder</a>
            </div>
            <div class="plot-grid">{plots}</div>
        </section>
    """


def render_srna_dsrna_identification_content(state: dict[str, Any]) -> str:
    output_dir = tool_output_root("srna-dsrna-identification", state)
    selected = stored_tool_samples(state, "srna-dsrna-identification", "srnaseq")
    groups = list(dict.fromkeys(sample["group"] for sample in selected if sample.get("included", True) and sample.get("group")))
    group_options = "".join(f'<option value="{esc(group)}">{esc(group)}</option>' for group in groups)
    no_group = '<option value="">First selected group</option>'
    return f"""
        <section class="hero">
            <div>
                <p class="eyebrow">Analysis tool</p>
                <h2>sRNA-based dsRNA Identification</h2>
                <p class="muted">Rank reference bins by bidirectional sRNA coverage, then plot top loci with flanking context.</p>
            </div>
            <div class="button-row top-actions">
                <a class="button-link" href="/reveal?path={url_for(output_dir)}">Open Output Folder</a>
            </div>
        </section>
        <div class="notice" id="status">Select preprocessed sRNA-seq samples, assign groups, and choose a genomic or multi-contig FASTA reference.</div>
        {render_module_dataset_selectors(state, 'srna-dsrna-identification')}
        {render_dsrna_reference_selector(state, 'srna-dsrna-identification')}
        <section class="panel">
            <div class="panel-heading">
                <div><p class="eyebrow">sRNA locus screening</p><h3>Bin Ranking Settings</h3></div>
            </div>
            <div class="settings-grid">
                <label class="field"><span>Scoring group</span><select id="srna-dsrna-score-group">{no_group}{group_options}</select></label>
                <label class="field"><span>Bin length (bp)</span><input id="srna-dsrna-bin-size" type="number" min="25" step="25" value="250"></label>
                <label class="field"><span>Context around hit (bp)</span><input id="srna-dsrna-context-bp" type="number" min="0" step="25" value="500"></label>
                <label class="field"><span>Top hits to plot</span><input id="srna-dsrna-top-n" type="number" min="1" step="1" value="20"></label>
                <label class="field"><span>Mismatches</span><select id="srna-dsrna-mismatches"><option selected>0</option><option>1</option><option>2</option><option>3</option></select></label>
                <label class="field"><span>Focus length</span><input id="srna-dsrna-focus-length" type="number" min="0" placeholder="optional, e.g. 21"></label>
                <label class="field"><span>Max multi-mappers</span><input id="srna-dsrna-max-multimappers" type="number" min="1" step="1" value="50"></label>
                <label class="field"><span>Threads</span><input id="srna-dsrna-threads" type="number" min="1" value="4"></label>
            </div>
            <div class="button-row">
                <label class="inline-check"><input id="srna-dsrna-filter-simple" type="checkbox" checked><span>Filter simple mapped reads</span></label>
                <label class="inline-check"><input id="srna-dsrna-report-all" type="checkbox"><span>Report all valid mappings</span></label>
            </div>
            {render_analysis_action('srna-dsrna-identification', 'runSrnaDsrnaIdentification()', 'Identify sRNA-supported dsRNA Loci')}
        </section>
        {render_srna_dsrna_identification_results(output_dir)}
    """


def render_dsrna_plotter_results(output_dir: Path) -> str:
    plots_dir = output_dir / "plots"
    contig_summary = output_dir / "tables" / "contig_summary.csv"
    sample_summary = output_dir / "tables" / "sample_summary.csv"
    contig_count = table_distinct_value_count(contig_summary, "contig")
    plot_cards = []
    for plot in prioritized_plot_paths(plots_dir.glob("*.png"), single_contig=contig_count == 1):
        plot_cards.append(
            f'<a class="plot-card" href="/reveal?path={url_for(plot)}"><strong>{esc(plot.stem)}</strong><img src="{output_file_url(plot)}" alt="{esc(plot.stem)}"></a>'
        )
    plots = "".join(plot_cards) or '<div class="empty-state">Run dsRNA Plotter to create directional coverage plots.</div>'
    return f"""
        <section class="panel">
            <div class="panel-heading">
                <div><p class="eyebrow">Results</p><h3>Directional Coverage Plots</h3></div>
                <a class="button-link" href="/reveal?path={url_for(plots_dir)}">Open Plot Folder</a>
            </div>
            <div class="plot-grid">{plots}</div>
        </section>
        <section class="panel">
            <div class="panel-heading">
                <div><p class="eyebrow">Summaries</p><h3>Mapping Tables</h3></div>
                <a class="button-link" href="/reveal?path={url_for(output_dir / 'tables')}">Open Table Folder</a>
            </div>
            <div class="path-row">
                <span class="badge">{'OK' if sample_summary.exists() else '--'}</span>
                <div><strong>Sample summary</strong><small>{esc(sample_summary)}</small></div>
                <a href="/reveal?path={url_for(sample_summary)}">Reveal</a>
            </div>
            <div class="path-row">
                <span class="badge">{'OK' if contig_summary.exists() else '--'}</span>
                <div><strong>Contig summary</strong><small>{esc(contig_summary)}</small></div>
                <a href="/reveal?path={url_for(contig_summary)}">Reveal</a>
            </div>
        </section>
    """


def render_module_content(module: ModuleSpec, state: dict[str, Any]) -> str:
    if module.key == "preprocessing":
        return render_preprocessing_content(state)
    if module.key == "srna-mapping":
        return render_srna_mapping_content(state)
    if module.key == "srna-control-filtering":
        return render_srna_control_filtering_content(state)
    if module.key == "srna-dsrna-identification":
        return render_srna_dsrna_identification_content(state)
    if module.key == "degradome-analysis":
        return render_degradome_content(state)
    if module.key == "target-prediction":
        return render_target_prediction_content(state)
    if module.key == "dsrna-identification":
        return render_dsrna_identification_content(state)

    output_dir = tool_output_root(module.key, state)

    if module.key == "dsrna-plotter":
        controls = f"""
        <section class="panel">
            <div class="panel-heading">
                <div>
                    <p class="eyebrow">Runnable action</p>
                    <h3>Grouped Directional Coverage Plots</h3>
                    <p>Run paired-end mapping with read-pair collapse, read 1 direction assignment, global CPM normalization, and one sense/antisense group plot per template contig.</p>
                </div>
            </div>
            <p class="muted">Include samples in the analysis sample set above. Samples with the same group are analyzed as biological replicates.</p>
            {render_analysis_action(module.key, f"runModule('{module.key}')", 'Run dsRNA Plotter')}
        </section>
        {render_dsrna_plotter_results(output_dir)}
        """
    elif module.key == "fasta-deduplication":
        controls = render_fasta_deduplication_controls(state, output_dir)
    else:
        raise ValueError("Unknown pipeline section.")

    return f"""
        <section class="hero">
            <div>
                <p class="eyebrow">Pipeline section</p>
                <h2>{esc(module.title)}</h2>
                <p class="muted">{esc(module.subtitle)}</p>
            </div>
            <div class="button-row top-actions">
                <a class="button-link" href="/reveal?path={url_for(output_dir)}">Open Output Folder</a>
            </div>
        </section>
        <div class="notice" id="status">Ready. Review the inputs before starting.</div>
        {render_module_dataset_selectors(state, module.key)}
        {render_dsrna_reference_selector(state, module.key) if module.key in {'dsrna-identification', 'dsrna-plotter'} else ''}
        {controls}
    """


def render_settings_content(state: dict[str, Any]) -> str:
    settings = plot_settings(state)
    tool_rows: list[str] = []
    missing_keys: list[str] = []
    for tool in EXTERNAL_TOOLS:
        status = external_tool_status(tool)
        ready = bool(status["ready"])
        if not ready:
            missing_keys.append(tool.key)
        executable_lines = "".join(
            f"<small><strong>{esc(name)}</strong>: {esc(path or 'not found')}</small>"
            for name, path in status["paths"].items()
        )
        tool_rows.append(
            f"""
            <div class="external-tool-row">
                <span class="badge {'ready' if ready else 'missing'}">{'Ready' if ready else 'Missing'}</span>
                <div>
                    <strong>{esc(tool.title)}</strong>
                    <small>{esc(tool.description)}</small>
                    <small>Used by: {esc(tool.used_by)}</small>
                    <div class="tool-paths">{executable_lines}</div>
                </div>
                <button type="button" class="secondary" onclick="installExternalTools(['{esc(tool.key)}'])" {'disabled' if ready else ''}>{'Installed' if ready else 'Install'}</button>
            </div>
            """
        )
    missing_json = json.dumps(missing_keys)
    return f"""
        <section class="hero">
            <div>
                <p class="eyebrow">Global configuration</p>
                <h2>Settings</h2>
                <p class="muted">Shared plot styling and external dependencies for every project.</p>
            </div>
        </section>
        <div class="notice" id="status">Changes to plot styling apply to plots generated after the settings are saved.</div>
        <section class="panel">
            <div class="panel-heading">
                <div>
                    <p class="eyebrow">Plot defaults</p>
                    <h3>Global Figure Style</h3>
                    <p>These values are used by integrated Matplotlib plots and passed to plotting subprocesses.</p>
                </div>
                <button type="button" onclick="savePlotSettings()">Save Plot Settings</button>
            </div>
            <div class="settings-grid">
                <label class="field"><span>Font family</span><input id="plot-font-family" value="{esc(settings['font_family'])}"></label>
                <label class="field"><span>Base font size</span><input id="plot-font-size" type="number" min="6" max="40" step="0.5" value="{esc(settings['font_size'])}"></label>
                <label class="field"><span>Title size</span><input id="plot-title-size" type="number" min="6" max="48" step="0.5" value="{esc(settings['title_size'])}"></label>
                <label class="field"><span>Line thickness</span><input id="plot-line-width" type="number" min="0.25" max="10" step="0.25" value="{esc(settings['line_width'])}"></label>
                <label class="field"><span>Grid thickness</span><input id="plot-grid-width" type="number" min="0.1" max="5" step="0.1" value="{esc(settings['grid_width'])}"></label>
                <label class="field"><span>Marker size</span><input id="plot-marker-size" type="number" min="1" max="30" step="0.5" value="{esc(settings['marker_size'])}"></label>
                <label class="field"><span>Resolution (DPI)</span><input id="plot-dpi" type="number" min="72" max="600" step="1" value="{esc(settings['dpi'])}"></label>
                <label class="field"><span>Figure width (in)</span><input id="plot-figure-width" type="number" min="4" max="30" step="0.5" value="{esc(settings['figure_width'])}"></label>
                <label class="field"><span>Figure height (in)</span><input id="plot-figure-height" type="number" min="3" max="20" step="0.5" value="{esc(settings['figure_height'])}"></label>
            </div>
        </section>
        <section class="panel">
            <div class="panel-heading">
                <div>
                    <p class="eyebrow">Dependencies</p>
                    <h3>External Tools</h3>
                    <p>Install missing command-line tools into the pipeline's local external-tools environment.</p>
                </div>
                <button type="button" onclick='installExternalTools({missing_json})' {'disabled' if not missing_keys else ''}>Install Missing</button>
            </div>
            <div class="external-tool-list">{''.join(tool_rows)}</div>
        </section>
    """


def render_nav(active: str) -> str:
    links: list[str] = [f'<a href="/" class="nav-main {"active" if active == "home" else ""}">Main Menu</a>']
    for group in MODULE_GROUPS:
        links.append(f'<p class="nav-group">{esc(group.title)}</p>')
        for key in group.keys:
            module = module_for_key(key)
            if module is None:
                continue
            links.append(
                f'<a href="/module/{esc(module.key)}" class="{"active" if active == module.key else ""}">{esc(module.title)}</a>'
            )
    links.append('<p class="nav-group">Configuration</p>')
    links.append(f'<a href="/settings" class="{"active" if active == "settings" else ""}">Settings</a>')
    return "\n".join(links)


def render_scripts() -> str:
    return """
    <script>
        const terminalStorageKey = 'inci-process-terminal';
        const jobLogOffsets = {};
        const jobLastRunning = {};
        const jobLastStatus = {};
        const pageScope = document.body.dataset.project + ':' + window.location.pathname;
        const activeJobKey = 'inci-active-job:' + pageScope;
        const statusKey = 'inci-run-status:' + pageScope;
        let runPending = false;

        function setRunBusy(busy) {
            runPending = busy;
            document.querySelectorAll('.run-primary, .clear-outputs').forEach(button => {
                if (!button.dataset.initialDisabled) button.dataset.initialDisabled = button.disabled ? 'yes' : 'no';
                button.disabled = busy || button.dataset.initialDisabled === 'yes';
            });
            document.querySelectorAll('.analysis-action-group').forEach(group => group.dataset.running = String(busy));
        }
        let lastStatusMessage = '';

        function terminalEntries() {
            try {
                return JSON.parse(localStorage.getItem(terminalStorageKey) || '[]');
            } catch {
                return [];
            }
        }

        function writeTerminal(message, level = 'info') {
            const now = new Date();
            const stamp = now.toLocaleTimeString([], {hour: '2-digit', minute: '2-digit', second: '2-digit'});
            const entries = terminalEntries();
            entries.push({stamp, message, level});
            localStorage.setItem(terminalStorageKey, JSON.stringify(entries.slice(-80)));
            renderTerminal();
        }

        function clearTerminal() {
            localStorage.removeItem(terminalStorageKey);
            renderTerminal();
            writeTerminal('Process console cleared.');
        }

        function toggleTerminal() {
            const terminal = document.getElementById('process-terminal');
            terminal.classList.toggle('collapsed');
            localStorage.setItem('inci-process-terminal-collapsed', terminal.classList.contains('collapsed') ? '1' : '0');
        }

        async function stopPipelineProcesses() {
            setStatus('Stopping running pipeline scripts...', 'warn');
            const response = await fetch('/stop-running', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({})
            });
            const data = await response.json();
            if (!response.ok || !data.ok) {
                setStatus(data.message || 'Could not stop running scripts.', 'error');
                return;
            }
            const count = data.count || 0;
            if (count === 0) {
                setStatus('No running pipeline scripts to stop.', 'warn');
                return;
            }
            setStatus('Stopped ' + count + ' running pipeline script' + (count === 1 ? '.' : 's.'), 'warn');
            if (data.labels && data.labels.length) {
                writeTerminal('Stopped: ' + data.labels.join(', '), 'warn');
            }
        }

        async function pollPipelineJob(jobId, reloadOnFinish = true) {
            if (!jobId) return;
            setRunBusy(true);
            sessionStorage.setItem(activeJobKey, JSON.stringify({
                jobId,
                path: window.location.pathname,
                reloadOnFinish
            }));
            try {
                const response = await fetch('/job-status?job_id=' + encodeURIComponent(jobId), {
                    cache: 'no-store', signal: AbortSignal.timeout(10000)
                });
                const data = await response.json();
                if (!response.ok || !data.ok) {
                    if (response.status === 404) {
                        sessionStorage.removeItem(activeJobKey);
                        setRunBusy(false);
                        setStatus('Run status is unavailable. The server may have restarted. Check the output folder before running again.', 'error');
                        return;
                    }
                    setStatus(data.message || 'Waiting to reconnect to the running analysis...', 'warn');
                    window.setTimeout(() => pollPipelineJob(jobId, reloadOnFinish), 3000);
                    return;
                }
                const job = data.job;
                const logs = job.logs || [];
                const offset = jobLogOffsets[jobId] || 0;
                logs.slice(offset).forEach((entry) => writeTerminal(entry.message, entry.level || 'info'));
                jobLogOffsets[jobId] = logs.length;

                const running = (job.running_processes || []).join(', ');
                if (running && jobLastRunning[jobId] !== running) {
                    jobLastRunning[jobId] = running;
                    writeTerminal('Running now: ' + running);
                }
                if (job.status === 'running') {
                    const current = job.current_message || (running ? 'Running: ' + running : job.label + ' is running...');
                    if (jobLastStatus[jobId] !== current) {
                        jobLastStatus[jobId] = current;
                        setStatus(current, 'success', false);
                    } else {
                        setStatus(current, 'success', false);
                    }
                    window.setTimeout(() => pollPipelineJob(jobId, reloadOnFinish), 2000);
                    return;
                }
                if (job.status === 'finished') {
                    sessionStorage.removeItem(activeJobKey);
                    setRunBusy(false);
                    setStatus(job.current_message || job.label + ' completed.', 'success');
                    if (job.result && job.result.manifest) {
                        writeTerminal('Trimmed manifest: ' + job.result.manifest);
                    }
                    if (job.result && job.label === 'sRNA mapping' && Number(job.result.contigs || 0) === 0) {
                        writeTerminal('sRNA mapping finished, but no contigs passed the current mapping/CPM filters.', 'warn');
                    }
                    if (reloadOnFinish) window.setTimeout(refreshCompletedRun, 500);
                    return;
                }
                sessionStorage.removeItem(activeJobKey);
                setRunBusy(false);
                if (job.status === 'stopped') {
                    setStatus(job.label + ' stopped. Partial outputs may be present in the output folder.', 'error');
                    sessionStorage.removeItem(runFormKey());
                    return;
                }
                setStatus(job.label + ' failed. ' + (job.result?.message || job.current_message || 'Check the process console for details.'), 'error');
                sessionStorage.removeItem(runFormKey());
                if (job.result && job.result.message) writeTerminal(job.result.message, 'error');
            } catch (error) {
                setStatus('Waiting to reconnect to the running analysis...', 'warn');
                window.setTimeout(() => pollPipelineJob(jobId, reloadOnFinish), 3000);
            }
        }

        function setStatus(message, level = 'info', terminal = true) {
            const status = document.getElementById('status');
            if (status) {
                status.textContent = message;
                status.classList.remove('status-info', 'status-success', 'status-warn', 'status-error');
                status.classList.add('status-' + level);
                status.setAttribute('role', level === 'error' ? 'alert' : 'status');
                status.setAttribute('aria-live', level === 'error' ? 'assertive' : 'polite');
            }
            if (level === 'success' || level === 'error') {
                sessionStorage.setItem(statusKey, JSON.stringify({message, level, path: window.location.pathname, timestamp: Date.now()}));
            }
            if (terminal && message !== lastStatusMessage) {
                lastStatusMessage = message;
                writeTerminal(message, level);
            }
        }

        function clearVisibleRunOutputs() {
            document.querySelectorAll('.plot-grid').forEach((container) => {
                container.innerHTML = '<div class="empty-state">New run started. Previous plots were cleared.</div>';
            });
            document.querySelectorAll('.table-scroll').forEach((container) => {
                container.innerHTML = '<div class="empty-state">New run started. Previous table results were cleared.</div>';
            });
            document.querySelectorAll('.sample-output-table').forEach((container) => {
                container.innerHTML = '<div class="empty-state">New run started. Previous output rows were cleared.</div>';
            });
            document.querySelectorAll('.path-row').forEach((row) => {
                const panel = row.closest('.panel');
                const label = panel?.querySelector('.eyebrow')?.textContent || '';
                if (/output|result/i.test(label)) row.hidden = true;
            });
            document.querySelectorAll('.panel').forEach((panel) => {
                if (!panel.querySelector('.plot-grid, .table-scroll, .sample-output-table, .path-row')) return;
                panel.querySelectorAll('.panel-heading a').forEach((link) => link.hidden = true);
            });
        }

        async function requestToolRun(endpoint, payload, failureMessage) {
            if (runPending) return null;
            setRunBusy(true);
            const formSnapshot = captureCurrentFormState();
            try {
                const response = await fetch(endpoint, {
                    method: 'POST',
                    headers: {'Content-Type': 'application/json'},
                    body: JSON.stringify(payload)
                });
                let data = {};
                try {
                    data = await response.json();
                } catch {
                    data = {};
                }
                if (!response.ok || !data.ok) {
                    setRunBusy(false);
                    setStatus(data.message || failureMessage, 'error');
                    if (data.details) writeTerminal(data.details, 'error');
                    return null;
                }
                sessionStorage.setItem(runFormKey(), JSON.stringify(formSnapshot));
                clearVisibleRunOutputs();
                return data;
            } catch (error) {
                setRunBusy(false);
                setStatus(failureMessage + ': ' + (error?.message || error), 'error');
                return null;
            }
        }

        function renderTerminal() {
            const log = document.getElementById('terminal-log');
            if (!log) return;
            const entries = terminalEntries();
            log.innerHTML = entries.map((entry) => (
                `<div class="terminal-line ${escapeHtml(entry.level || 'info')}"><span>${escapeHtml(entry.stamp)}</span>${escapeHtml(entry.message)}</div>`
            )).join('');
            log.scrollTop = log.scrollHeight;
        }

        function escapeHtml(value) {
            return String(value)
                .replaceAll('&', '&amp;')
                .replaceAll('<', '&lt;')
                .replaceAll('>', '&gt;')
                .replaceAll('"', '&quot;')
                .replaceAll("'", '&#039;');
        }

        const clearRestoreKey = () => 'inci-clear-restore:' + pageScope;
        const runFormKey = () => 'inci-run-form:' + pageScope;

        function refreshToolPage() {
            const url = new URL(window.location.href);
            url.searchParams.set('_results', Date.now().toString());
            window.location.replace(url.toString());
        }

        function refreshCompletedRun() {
            const snapshot = sessionStorage.getItem(runFormKey());
            if (snapshot) {
                sessionStorage.setItem(clearRestoreKey(), snapshot);
                sessionStorage.removeItem(runFormKey());
            }
            sessionStorage.removeItem(activeJobKey);
            refreshToolPage();
        }

        function captureCurrentFormState() {
            const fields = {};
            document.querySelectorAll('input[id], select[id], textarea[id]').forEach((element) => {
                if (element.type === 'file') return;
                fields[element.id] = (element.type === 'checkbox' || element.type === 'radio')
                    ? {kind: 'checked', value: element.checked}
                    : {kind: 'value', value: element.value};
            });
            const trimmingRows = Array.from(document.querySelectorAll('[data-trim-row]')).map((row) => ({
                read1: row.querySelector('[data-field="read1"]')?.value || '',
                read2: row.querySelector('[data-field="read2"]')?.value || ''
            }));
            return {fields, trimmingRows};
        }

        function restoreCurrentFormState() {
            const oneShotKey = clearRestoreKey();
            const oneShotRaw = sessionStorage.getItem(oneShotKey);
            const raw = oneShotRaw || sessionStorage.getItem(runFormKey());
            if (!raw) return;
            try {
                const snapshot = JSON.parse(raw);
                Object.entries(snapshot.fields || {}).forEach(([id, saved]) => {
                    const element = document.getElementById(id);
                    if (!element) return;
                    if (saved.kind === 'checked') element.checked = Boolean(saved.value);
                    else element.value = saved.value ?? '';
                });
                const rowsHost = document.getElementById('trim-rows');
                if (rowsHost && Array.isArray(snapshot.trimmingRows) && snapshot.trimmingRows.length) {
                    while (rowsHost.querySelectorAll('[data-trim-row]').length < snapshot.trimmingRows.length) addTrimRow();
                    while (rowsHost.querySelectorAll('[data-trim-row]').length > snapshot.trimmingRows.length) {
                        rowsHost.querySelector('[data-trim-row]:last-child')?.remove();
                    }
                    Array.from(rowsHost.querySelectorAll('[data-trim-row]')).forEach((row, index) => {
                        const saved = snapshot.trimmingRows[index] || {};
                        const read1 = row.querySelector('[data-field="read1"]');
                        const read2 = row.querySelector('[data-field="read2"]');
                        if (read1) read1.value = saved.read1 || '';
                        if (read2) read2.value = saved.read2 || '';
                    });
                }
            } catch {
                // A stale snapshot should never prevent the page from loading.
            } finally {
                if (oneShotRaw) sessionStorage.removeItem(oneShotKey);
            }
        }

        document.addEventListener('DOMContentLoaded', () => {
            const search = document.getElementById('tool-search');
            search?.addEventListener('input', () => {
                const query = search.value.trim().toLowerCase();
                document.querySelectorAll('nav a, .module-card').forEach(link => {
                    link.hidden = !link.textContent.toLowerCase().includes(query);
                });
            });
            document.querySelector('nav a.active')?.setAttribute('aria-current', 'page');
            const sectionNav = document.getElementById('page-sections');
            document.querySelectorAll('main .panel h3').forEach((heading, index) => {
                heading.id = 'section-' + index;
                const link = document.createElement('a');
                link.href = '#' + heading.id;
                link.textContent = heading.textContent;
                sectionNav?.appendChild(link);
            });
            document.querySelectorAll('.compact-field input, .compact-field select').forEach(field => {
                if (!field.labels?.length && !field.getAttribute('aria-label')) {
                    field.setAttribute('aria-label', field.closest('.compact-field').textContent.trim() || 'Sample setting');
                }
            });
            const status = document.getElementById('status');
            const actionStatusSlot = document.querySelector('[data-analysis-status-slot]');
            if (status && actionStatusSlot) {
                status.classList.add('analysis-run-status');
                status.setAttribute('role', 'status');
                status.setAttribute('aria-live', 'polite');
                actionStatusSlot.appendChild(status);
            }
            renderTerminal();
            try {
                const savedStatus = JSON.parse(sessionStorage.getItem(statusKey) || 'null');
                if (savedStatus && savedStatus.path === window.location.pathname) {
                    setStatus(savedStatus.message, savedStatus.level, false);
                    if (new URL(window.location.href).searchParams.has('_results')) {
                        status?.scrollIntoView({block: 'center', behavior: 'instant'});
                    }
                }
            } catch {
                sessionStorage.removeItem(statusKey);
            }
            const terminal = document.getElementById('process-terminal');
            if (localStorage.getItem('inci-process-terminal-collapsed') !== '0') {
                terminal.classList.add('collapsed');
            }
            if (terminalEntries().length === 0) {
                writeTerminal('INCI interface ready.');
            }
            restoreCurrentFormState();
            syncTrimMode();
            syncTrimTool();
            try {
                const activeJob = JSON.parse(sessionStorage.getItem(activeJobKey) || 'null');
                if (activeJob && activeJob.jobId && activeJob.path === window.location.pathname) {
                    pollPipelineJob(activeJob.jobId, activeJob.reloadOnFinish !== false);
                }
            } catch {
                sessionStorage.removeItem(activeJobKey);
            }
        });

        async function savePath(key) {
            const input = document.getElementById(key);
            setStatus('Saving local path for ' + key + '...');
            const response = await fetch('/save', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({key, path: input.value.trim()})
            });
            if (!response.ok) {
                setStatus('Could not save that path.', 'error');
                return;
            }
            writeTerminal('Saved local path for ' + key + '.');
            window.location.reload();
        }

        async function saveProjectName() {
            const input = document.getElementById('project-name');
            const dirInput = document.getElementById('project-dir');
            const name = (input?.value || '').trim();
            const project_dir = (dirInput?.value || '').trim();
            if (!name) {
                setStatus('Enter a project name before saving.', 'warn');
                return;
            }
            setStatus('Saving project name...');
            const response = await fetch('/save-project', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({name, project_dir})
            });
            const data = await response.json();
            if (!response.ok || !data.ok) {
                setStatus(data.message || 'Could not save project name.', 'error');
                return;
            }
            setStatus(data.message || 'Project saved.');
            window.location.reload();
        }

        async function savePlotSettings() {
            const payload = {
                font_family: document.getElementById('plot-font-family')?.value.trim() || 'DejaVu Sans',
                font_size: Number(document.getElementById('plot-font-size')?.value || 11),
                title_size: Number(document.getElementById('plot-title-size')?.value || 14),
                line_width: Number(document.getElementById('plot-line-width')?.value || 2),
                grid_width: Number(document.getElementById('plot-grid-width')?.value || 0.7),
                marker_size: Number(document.getElementById('plot-marker-size')?.value || 6),
                dpi: Number(document.getElementById('plot-dpi')?.value || 220),
                figure_width: Number(document.getElementById('plot-figure-width')?.value || 12),
                figure_height: Number(document.getElementById('plot-figure-height')?.value || 5)
            };
            setStatus('Saving global plot settings...');
            const response = await fetch('/save-settings', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({plot_settings: payload})
            });
            const data = await response.json();
            if (!response.ok || !data.ok) {
                setStatus(data.message || 'Could not save plot settings.', 'error');
                return;
            }
            setStatus(data.message || 'Plot settings saved.');
            writeTerminal('Global plot settings updated.');
        }

        async function installExternalTools(toolKeys) {
            if (!toolKeys || !toolKeys.length) {
                setStatus('All registered external tools are already available.');
                return;
            }
            setStatus('Starting external tool installation...');
            const response = await fetch('/install-tools', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({tool_keys: toolKeys})
            });
            const data = await response.json();
            if (!response.ok || !data.ok) {
                setStatus(data.message || 'Could not start external tool installation.', 'error');
                return;
            }
            setStatus(data.message || 'External tool installation started.');
            if (data.job_id) pollPipelineJob(data.job_id, true);
        }

        async function browseProjectDir() {
            setStatus('Opening local folder chooser for the project directory...');
            const response = await fetch('/browse-project-dir');
            const data = await response.json();
            if (data.path) {
                document.getElementById('project-dir').value = data.path;
                writeTerminal('Selected project directory ' + data.path);
            } else {
                setStatus(data.message || 'No folder selected.', 'warn');
            }
        }

        async function resetCurrentPage(moduleKey) {
            const label = moduleKey === 'preprocessing' ? 'RNA-seq preprocessing' : moduleKey;
            if (!window.confirm('Clear previous ' + label + ' outputs? Current inputs, sample selections, and parameter values will be kept.')) {
                return;
            }
            const restoreKey = clearRestoreKey();
            const snapshot = captureCurrentFormState();
            if (!(await saveVisibleToolSampleSettings())) return;
            sessionStorage.setItem(restoreKey, JSON.stringify(snapshot));
            setStatus('Clearing previous ' + label + ' outputs...', 'warn');
            const response = await fetch('/reset-page', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({module_key: moduleKey})
            });
            const data = await response.json();
            if (!response.ok || !data.ok) {
                sessionStorage.removeItem(restoreKey);
                setStatus(data.message || 'Could not clear previous outputs.', 'error');
                return;
            }
            setStatus(data.message || 'Previous outputs cleared.', 'warn');
            window.location.reload();
        }

        async function browseFile(key) {
            setStatus('Opening local file chooser for ' + key + '...');
            const response = await fetch('/browse?key=' + encodeURIComponent(key));
            const data = await response.json();
            if (data.path) {
                document.getElementById(key).value = data.path;
                writeTerminal('Selected ' + data.path);
                await savePath(key);
            } else {
                setStatus(data.message || 'No file selected.', 'warn');
            }
        }

        async function savePreprocessed(key) {
            const input = document.getElementById('prep-' + key);
            setStatus('Saving preprocessed path for ' + key + '...');
            const response = await fetch('/save-preprocessed', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({key, path: input.value.trim()})
            });
            if (!response.ok) {
                setStatus('Could not save that preprocessed path.', 'error');
                return;
            }
            writeTerminal('Saved preprocessed path for ' + key + '.');
            window.location.reload();
        }

        async function browsePreprocessed(key) {
            setStatus('Opening local file chooser for preprocessed ' + key + '...');
            const response = await fetch('/browse-preprocessed?key=' + encodeURIComponent(key));
            const data = await response.json();
            if (data.path) {
                document.getElementById('prep-' + key).value = data.path;
                writeTerminal('Selected preprocessed file ' + data.path);
                await savePreprocessed(key);
            } else {
                setStatus(data.message || 'No file selected.', 'warn');
            }
        }

        async function browseToolSample(moduleKey, inputKey, field) {
            setStatus('Opening local file chooser...');
            const response = await fetch('/browse-tool-sample?module_key=' + encodeURIComponent(moduleKey) + '&input_key=' + encodeURIComponent(inputKey));
            const data = await response.json();
            if (!response.ok || !data.path) {
                setStatus(data.message || 'No file selected.', 'warn');
                return;
            }
            document.getElementById('tool-' + field + '-' + moduleKey + '-' + inputKey).value = data.path;
            setStatus('File ready to add.');
        }

        async function addToolBatch(moduleKey, inputKey, batchId) {
            setStatus('Adding preprocessing batch...');
            const response = await fetch('/add-tool-batch', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({module_key: moduleKey, input_key: inputKey, batch_id: batchId})
            });
            const data = await response.json();
            if (!response.ok || !data.ok) {
                setStatus(data.message || 'Could not add that batch.', 'error');
                return;
            }
            window.location.reload();
        }

        async function addToolSample(moduleKey, inputKey) {
            const read1 = document.getElementById('tool-read1-' + moduleKey + '-' + inputKey)?.value.trim() || '';
            const read2 = document.getElementById('tool-read2-' + moduleKey + '-' + inputKey)?.value.trim() || '';
            setStatus('Adding sample...');
            const response = await fetch('/add-tool-sample', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({module_key: moduleKey, input_key: inputKey, read1, read2})
            });
            const data = await response.json();
            if (!response.ok || !data.ok) {
                setStatus(data.message || 'Could not add that sample.', 'error');
                return;
            }
            window.location.reload();
        }

        function toolSampleSettings(moduleKey, inputKey) {
            const list = document.getElementById('tool-samples-' + moduleKey + '-' + inputKey);
            const samples = {};
            list?.querySelectorAll('[data-tool-sample]').forEach((row) => {
                samples[row.dataset.sampleId] = {
                    group: row.querySelector('[data-field="group"]')?.value.trim() || '',
                    included: row.querySelector('[data-field="included"]')?.checked || false
                };
            });
            return samples;
        }

        function includedToolSampleIds(moduleKey, inputKey) {
            const list = document.getElementById('tool-samples-' + moduleKey + '-' + inputKey);
            return Array.from(list?.querySelectorAll('[data-tool-sample]') || [])
                .filter((row) => row.querySelector('[data-field="included"]')?.checked)
                .map((row) => row.dataset.sampleName || '')
                .filter(Boolean);
        }

        function refreshToolReplicates(moduleKey, inputKey) {
            const list = document.getElementById('tool-samples-' + moduleKey + '-' + inputKey);
            const counts = {};
            let includedCount = 0;
            list?.querySelectorAll('[data-tool-sample]').forEach((row) => {
                const included = row.querySelector('[data-field="included"]')?.checked || false;
                row.classList.toggle('is-excluded', !included);
                const output = row.querySelector('[data-replicate]');
                if (!included) {
                    if (output) output.textContent = 'Excluded';
                    return;
                }
                includedCount += 1;
                const group = row.querySelector('[data-field="group"]')?.value.trim().toLocaleLowerCase() || row.dataset.sampleId;
                counts[group] = (counts[group] || 0) + 1;
                if (output) output.textContent = counts[group];
            });
            const count = list?.closest('.selected-samples')?.querySelector('[data-included-count]');
            if (count) count.textContent = includedCount;
        }

        async function saveToolSampleSettings(moduleKey, inputKey, reload = true) {
            setStatus('Saving analysis sample set...');
            const response = await fetch('/update-tool-samples', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({module_key: moduleKey, input_key: inputKey, samples: toolSampleSettings(moduleKey, inputKey)})
            });
            const data = await response.json();
            if (!response.ok || !data.ok) {
                setStatus(data.message || 'Could not save the analysis sample set.', 'error');
                return false;
            }
            if (reload) window.location.reload();
            return true;
        }

        async function saveVisibleToolSampleSettings() {
            const lists = Array.from(document.querySelectorAll('[id^="tool-samples-"]'));
            for (const list of lists) {
                const parts = list.id.replace('tool-samples-', '').split('-');
                const inputKey = parts.pop();
                const moduleKey = parts.join('-');
                if (!(await saveToolSampleSettings(moduleKey, inputKey, false))) return false;
            }
            return true;
        }

        async function removeToolSample(moduleKey, inputKey, sampleId) {
            const response = await fetch('/remove-tool-sample', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({module_key: moduleKey, input_key: inputKey, sample_id: sampleId})
            });
            const data = await response.json();
            if (!response.ok || !data.ok) {
                setStatus(data.message || 'Could not remove that sample.', 'error');
                return;
            }
            window.location.reload();
        }

        function dsrnaReferencePayload(moduleKey) {
            return {
                module_key: moduleKey,
                reference_text: document.getElementById('dsrna-reference-text-' + moduleKey)?.value || '',
                reference_fasta: document.getElementById('dsrna-reference-file-' + moduleKey)?.value.trim() || ''
            };
        }

        function moduleRunPayload(moduleKey) {
            if (moduleKey === 'fasta-deduplication') {
                return {
                    input_fasta: document.getElementById('fasta-dedup-input')?.value.trim() || '',
                    identity_percent: Number(document.getElementById('fasta-dedup-identity')?.value || 95),
                    word_size: document.getElementById('fasta-dedup-word-size')?.value || 'auto',
                    threads: Number(document.getElementById('fasta-dedup-threads')?.value || 7),
                    memory_mb: Number(document.getElementById('fasta-dedup-memory')?.value || 0),
                    strand_mode: document.getElementById('fasta-dedup-strand')?.value || 'both'
                };
            }
            if (moduleKey !== 'dsrna-identification') return {};
            return {
                scoring_group: document.getElementById('dsrna-score-group')?.value || '',
                bin_size: Number(document.getElementById('dsrna-bin-size')?.value || 250),
                plot_context_bp: Number(document.getElementById('dsrna-context-bp')?.value || 500),
                top_n: Number(document.getElementById('dsrna-top-n')?.value || 20),
                sirna_annotations_fasta: document.getElementById('dsrna-sirna-annotations-fasta')?.value.trim() || ''
            };
        }

        async function saveDsrnaReference(moduleKey, reload = true) {
            setStatus('Saving reference sequence...');
            const response = await fetch('/save-dsrna-reference', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify(dsrnaReferencePayload(moduleKey))
            });
            const data = await response.json();
            if (!response.ok || !data.ok) {
                setStatus(data.message || 'Could not save the reference.', 'error');
                return false;
            }
            setStatus(data.message || 'Reference saved.');
            if (reload) window.location.reload();
            return true;
        }

        async function browseFastaDedupInput() {
            setStatus('Opening local FASTA chooser...');
            const response = await fetch('/browse-fasta-dedup-input');
            const data = await response.json();
            if (!response.ok || !data.path) {
                setStatus(data.message || 'No FASTA selected.', 'warn');
                return;
            }
            const input = document.getElementById('fasta-dedup-input');
            if (input) input.value = data.path;
            setStatus('FASTA ready for deduplication.');
        }

        async function browseDsrnaReference(moduleKey) {
            setStatus('Opening local FASTA chooser...');
            const response = await fetch('/browse-dsrna-reference?module_key=' + encodeURIComponent(moduleKey));
            const data = await response.json();
            if (!response.ok || !data.path) {
                setStatus(data.message || 'No FASTA selected.', 'warn');
                return;
            }
            window.location.reload();
        }

        async function runModule(moduleKey) {
            if (!(await saveVisibleToolSampleSettings())) return;
            if (document.getElementById('dsrna-reference-file-' + moduleKey) || document.getElementById('dsrna-reference-text-' + moduleKey)) {
                const saved = await saveDsrnaReference(moduleKey, false);
                if (!saved) return;
            }
            setStatus('Running ' + moduleKey + '...');
            const data = await requestToolRun('/run-module', {module_key: moduleKey, params: moduleRunPayload(moduleKey)}, 'Run failed.');
            if (!data) return;
            const outputCount = data.outputs ? Object.keys(data.outputs).length : 0;
            setStatus(data.message || ('Run finished with ' + outputCount + ' output files.'), 'success');
            if (data.job_id) {
                pollPipelineJob(data.job_id, true);
                return;
            }
            if (data.warnings && data.warnings.length) {
                data.warnings.forEach((warning) => writeTerminal(warning, 'warn'));
            }
            window.setTimeout(refreshCompletedRun, 250);
        }

        function syncTrimMode() {
            const mode = document.getElementById('trim-mode')?.value || 'paired';
            document.body.classList.toggle('trim-single', mode === 'single');
            syncTrimTool();
        }

        function syncTrimTool() {
            const tool = document.getElementById('trim-tool')?.value || 'trim_galore';
            document.body.classList.toggle('trim-cutadapt', tool === 'cutadapt');
        }

        function trimPayload() {
            const samples = Array.from(document.querySelectorAll('[data-trim-row]')).map((row) => ({
                read1: row.querySelector('[data-field="read1"]').value.trim(),
                read2: row.querySelector('[data-field="read2"]')?.value.trim() || ''
            })).filter((sample) => sample.read1 || sample.read2);
            return {
                batch_name: document.getElementById('trim-batch-name')?.value.trim() || '',
                data_type: document.getElementById('trim-data-type')?.value || 'rnaseq',
                tool: document.getElementById('trim-tool')?.value || 'trim_galore',
                mode: document.getElementById('trim-mode')?.value || 'paired',
                quality: Number(document.getElementById('trim-quality')?.value || 20),
                length: Number(document.getElementById('trim-length')?.value || 18),
                adapter1: document.getElementById('cutadapt-adapter1')?.value.trim() || '',
                adapter2: document.getElementById('cutadapt-adapter2')?.value.trim() || '',
                cutadapt_error_rate: Number(document.getElementById('cutadapt-error-rate')?.value || 0.1),
                cutadapt_overlap: Number(document.getElementById('cutadapt-overlap')?.value || 3),
                cutadapt_max_n: document.getElementById('cutadapt-max-n')?.value.trim() || '',
                cutadapt_trim_n: document.getElementById('cutadapt-trim-n')?.checked ?? true,
                cutadapt_cores: Number(document.getElementById('cutadapt-cores')?.value || 1),
                cutadapt_pair_filter: document.getElementById('cutadapt-pair-filter')?.value || 'any',
                samples
            };
        }

        function addTrimRow() {
            const rows = document.getElementById('trim-rows');
            if (!rows) return;
            const row = document.createElement('div');
            row.className = 'trim-row';
            row.setAttribute('data-trim-row', '');
            row.innerHTML = `
                <div class="file-picker">
                    <input data-field="read1" placeholder="Read 1 FASTQ">
                    <button type="button" class="secondary" onclick="browseTrimFile(this, 'read1')">Browse</button>
                </div>
                <div class="file-picker paired-only">
                    <input data-field="read2" placeholder="Read 2 FASTQ">
                    <button type="button" class="secondary" onclick="browseTrimFile(this, 'read2')">Browse</button>
                </div>
                <button type="button" class="ghost icon-button" title="Remove row" onclick="removeTrimRow(this)">x</button>
            `;
            rows.appendChild(row);
            syncTrimMode();
            syncTrimTool();
        }

        function removeTrimRow(button) {
            const rows = document.querySelectorAll('[data-trim-row]');
            if (rows.length <= 1) {
                const row = button.closest('[data-trim-row]');
                row.querySelectorAll('input').forEach((input) => input.value = '');
                return;
            }
            button.closest('[data-trim-row]').remove();
        }

        async function browseTrimFile(button, field) {
            setStatus('Opening local file chooser for trimming input...');
            const response = await fetch('/browse-trim-file');
            const data = await response.json();
            if (data.path) {
                const row = button.closest('[data-trim-row]');
                row.querySelector('[data-field="' + field + '"]').value = data.path;
                writeTerminal('Selected trimming input ' + data.path);
            } else {
                setStatus(data.message || 'No file selected.', 'warn');
            }
        }

        async function runTrimming() {
            const payload = trimPayload();
            const toolLabel = payload.tool === 'cutadapt' ? 'Cutadapt' : 'Trim Galore';
            setStatus('Starting ' + toolLabel + ' in the background...');
            const data = await requestToolRun('/run-trimming', payload, toolLabel + ' run failed.');
            if (!data) return;
            setStatus(data.message || (toolLabel + ' batch started.'), 'success');
            if (data.job_id) pollPipelineJob(data.job_id, true);
        }

        async function browseSrnaReference() {
            setStatus('Opening local file chooser for sRNA reference FASTA...');
            const response = await fetch('/browse-srna-reference');
            const data = await response.json();
            if (data.path) {
                document.getElementById('srna-reference-fasta').value = data.path;
                writeTerminal('Selected sRNA reference ' + data.path);
            } else {
                setStatus(data.message || 'No file selected.', 'warn');
            }
        }

        async function browseSirnaAnnotations(moduleKey, inputId) {
            setStatus('Opening local file chooser for siRNA annotations...');
            const response = await fetch('/browse-sirna-annotations?module_key=' + encodeURIComponent(moduleKey));
            const data = await response.json();
            if (data.path) {
                const input = document.getElementById(inputId);
                if (input) input.value = data.path;
                setStatus('siRNA annotation FASTA selected.');
            } else {
                setStatus(data.message || 'No file selected.', 'warn');
            }
        }

        function srnaPayload() {
            const sampleIds = includedToolSampleIds('srna-mapping', 'srnaseq');
            const dsrnaSampleIds = includedToolSampleIds('srna-mapping', 'rnaseq');
            return {
                sample_ids: sampleIds,
                dsrna_sample_ids: dsrnaSampleIds,
                include_dsrna_overlay: dsrnaSampleIds.length > 0,
                reference_text: document.getElementById('srna-reference-text')?.value || '',
                reference_fasta: document.getElementById('srna-reference-fasta')?.value.trim() || '',
                sirna_annotations_fasta: document.getElementById('srna-sirna-annotations-fasta')?.value.trim() || '',
                mismatches: Number(document.getElementById('srna-mismatches')?.value || 0),
                report_all: document.getElementById('srna-report-all')?.checked || false,
                filter_simple: document.getElementById('srna-filter-simple')?.checked !== false,
                focus_length: Number(document.getElementById('srna-focus-length')?.value || 0),
                cpm_threshold: Number(document.getElementById('srna-cpm-threshold')?.value || 0),
                export_min_cpm: Number(document.getElementById('srna-export-min-cpm')?.value || 0),
                export_length: Number(document.getElementById('srna-export-length')?.value || 0),
                export_top_n: Number(document.getElementById('srna-export-top-n')?.value || 20),
                threads: Number(document.getElementById('srna-threads')?.value || 4)
            };
        }

        async function runSrnaMapping() {
            if (!(await saveVisibleToolSampleSettings())) return;
            setStatus('Starting sRNA mapping in the background...');
            const data = await requestToolRun('/run-srna-mapping', srnaPayload(), 'Could not start sRNA mapping.');
            if (!data) return;
            setStatus(data.message || 'sRNA mapping started.', 'success');
            if (data.job_id) pollPipelineJob(data.job_id, true);
        }

        function srnaDsrnaPayload() {
            const sampleIds = includedToolSampleIds('srna-dsrna-identification', 'srnaseq');
            return {
                sample_ids: sampleIds,
                scoring_group: document.getElementById('srna-dsrna-score-group')?.value || '',
                bin_size: Number(document.getElementById('srna-dsrna-bin-size')?.value || 250),
                plot_context_bp: Number(document.getElementById('srna-dsrna-context-bp')?.value || 500),
                top_n: Number(document.getElementById('srna-dsrna-top-n')?.value || 20),
                mismatches: Number(document.getElementById('srna-dsrna-mismatches')?.value || 0),
                focus_length: Number(document.getElementById('srna-dsrna-focus-length')?.value || 0),
                max_multimappers: Number(document.getElementById('srna-dsrna-max-multimappers')?.value || 50),
                filter_simple: document.getElementById('srna-dsrna-filter-simple')?.checked !== false,
                report_all: document.getElementById('srna-dsrna-report-all')?.checked || false,
                threads: Number(document.getElementById('srna-dsrna-threads')?.value || 4)
            };
        }

        async function runSrnaDsrnaIdentification() {
            if (!(await saveVisibleToolSampleSettings())) return;
            const saved = await saveDsrnaReference('srna-dsrna-identification', false);
            if (!saved) return;
            setStatus('Starting sRNA-based dsRNA identification in the background...');
            const data = await requestToolRun('/run-srna-dsrna-identification', srnaDsrnaPayload(), 'Could not start sRNA-based dsRNA identification.');
            if (!data) return;
            setStatus(data.message || 'sRNA-based dsRNA identification started.', 'success');
            if (data.job_id) pollPipelineJob(data.job_id, true);
        }

        async function browseSrnaControlFile(inputId) {
            setStatus('Opening local file chooser...');
            const response = await fetch('/browse-srna-control-file');
            const data = await response.json();
            if (data.path) {
                document.getElementById(inputId).value = data.path;
                writeTerminal('Selected ' + data.path);
            } else {
                setStatus(data.message || 'No file selected.', 'warn');
            }
        }

        function srnaControlPayload() {
            return {
                unique_text: document.getElementById('control-unique-text')?.value || '',
                unique_fasta: document.getElementById('control-unique-fasta')?.value.trim() || '',
                dsrna_reference_text: document.getElementById('control-dsrna-text')?.value || '',
                dsrna_reference_fasta: document.getElementById('control-dsrna-fasta')?.value.trim() || '',
                control_genome_text: document.getElementById('control-genome-text')?.value || '',
                control_genome_fasta: document.getElementById('control-genome-fasta')?.value.trim() || '',
                control_srna_text: document.getElementById('control-srna-text')?.value || '',
                control_srna_fasta: document.getElementById('control-srna-fasta')?.value.trim() || '',
                dsrna_mismatches: Number(document.getElementById('control-dsrna-mismatches')?.value || 0),
                control_mismatches: Number(document.getElementById('control-reference-mismatches')?.value || 0),
                length_filter: Number(document.getElementById('control-length-filter')?.value || 0),
                threads: Number(document.getElementById('control-threads')?.value || 4),
                index_mode: document.getElementById('control-index-mode')?.value || 'auto',
                remove_low_complexity: document.getElementById('control-remove-low-complexity')?.checked || false,
                remove_control_mappers: document.getElementById('control-remove-control-mappers')?.checked || false,
                collapse_contained: document.getElementById('control-collapse-contained')?.checked || false
            };
        }

        async function runSrnaControlFiltering() {
            setStatus('Starting sRNA control mapping and filtering in the background...');
            const data = await requestToolRun('/run-srna-control-filtering', srnaControlPayload(), 'Could not start sRNA control filtering.');
            if (!data) return;
            setStatus(data.message || 'sRNA control filtering started.', 'success');
            if (data.job_id) pollPipelineJob(data.job_id, true);
        }

        async function browseDegradomeInput(inputId) {
            setStatus('Opening local file chooser...');
            const response = await fetch('/browse-degradome-input');
            const data = await response.json();
            if (data.path) {
                document.getElementById(inputId).value = data.path;
                writeTerminal('Selected ' + data.path);
            } else {
                setStatus(data.message || 'No file selected.', 'warn');
            }
        }

        function degradomePayload() {
            const samples = Array.from(document.querySelectorAll('[data-tool-sample]')).map((row) => ({
                sample_id: row.dataset.sampleName || '',
                group: row.querySelector('[data-field="group"]')?.value.trim() || '',
                replicate: row.querySelector('[data-replicate]')?.textContent.trim() || '',
                path: row.dataset.read1 || '',
                included: row.querySelector('[data-field="included"]')?.checked || false
            })).filter((sample) => sample.path && sample.included);
            return {
                srna_text: document.getElementById('degradome-srna-text')?.value || '',
                srna_fasta: document.getElementById('degradome-srna-fasta')?.value.trim() || '',
                transcript_text: document.getElementById('degradome-transcript-text')?.value || '',
                transcript_fasta: document.getElementById('degradome-transcript-fasta')?.value.trim() || '',
                samples,
                ignore_query_pos1: document.getElementById('degradome-ignore-q1')?.checked || false,
                slice_positions: document.getElementById('degradome-slice-mode')?.value || '10',
                mfe_ratio_cutoff: Number(document.getElementById('degradome-mfe-cutoff')?.value || 0.70),
                sort_by: document.getElementById('degradome-sort-by')?.value || 'mfe_ratio',
                pvalue_method: document.getElementById('degradome-pvalue-method')?.value || 'transcript_peak_empirical',
                pvalue_cutoff: Number(document.getElementById('degradome-pvalue-cutoff')?.value || 0.05),
                threads: Number(document.getElementById('degradome-threads')?.value || 4),
                max_transcript_plots: Number(document.getElementById('degradome-max-plots')?.value || 80),
                compact_srna_markers: document.getElementById('degradome-compact-markers')?.checked ?? true,
                marker_neighborhood_nt: Number(document.getElementById('degradome-marker-window')?.value || 3)
            };
        }

        async function runDegradomeAnalysis() {
            if (!(await saveVisibleToolSampleSettings())) return;
            setStatus('Starting degradome analysis in the background...');
            const data = await requestToolRun('/run-degradome-analysis', degradomePayload(), 'Could not start degradome analysis.');
            if (!data) return;
            setStatus(data.message || 'Degradome analysis started.', 'success');
            if (data.job_id) pollPipelineJob(data.job_id, true);
        }

        async function browseTargetInput(inputId) {
            setStatus('Opening local file chooser...');
            const response = await fetch('/browse-degradome-input');
            const data = await response.json();
            if (data.path) {
                document.getElementById(inputId).value = data.path;
                writeTerminal('Selected ' + data.path);
            } else {
                setStatus(data.message || 'No file selected.', 'warn');
            }
        }

        function targetPredictionPayload() {
            return {
                srna_text: document.getElementById('target-srna-text')?.value || '',
                srna_fasta: document.getElementById('target-srna-fasta')?.value.trim() || '',
                transcript_text: document.getElementById('target-transcript-text')?.value || '',
                transcript_fasta: document.getElementById('target-transcript-fasta')?.value.trim() || '',
                ignore_query_pos1: document.getElementById('target-ignore-q1')?.checked || false,
                mfe_ratio_cutoff: Number(document.getElementById('target-mfe-cutoff')?.value || 0.70),
                max_allen_score: document.getElementById('target-max-allen')?.value.trim() || '',
                max_mismatches: document.getElementById('target-max-mismatches')?.value.trim() || '',
                sort_by: document.getElementById('target-sort-by')?.value || 'mfe_ratio',
                plot_metric: document.getElementById('target-plot-metric')?.value || 'mfe_ratio',
                max_transcript_plots: Number(document.getElementById('target-max-transcript-plots')?.value || 0),
                threads: Number(document.getElementById('target-threads')?.value || 4)
            };
        }

        async function runTargetPrediction() {
            setStatus('Starting target prediction in the background...');
            const data = await requestToolRun('/run-target-prediction', targetPredictionPayload(), 'Could not start target prediction.');
            if (!data) return;
            setStatus(data.message || 'Target prediction started.', 'success');
            if (data.job_id) pollPipelineJob(data.job_id, true);
        }
    </script>
    """


def render_shell(title: str, active: str, content: str, port: int) -> str:
    shell_state = load_state()
    shell_project_name = project_name(shell_state)
    shell_project_root = project_output_root(shell_state)
    shell_project_dir = project_custom_dir(shell_state)
    return f"""<!doctype html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <meta name="theme-color" content="#0f5f5a">
    <link rel="icon" href="/favicon.svg" type="image/svg+xml">
    <title>{esc(title)}</title>
    <style>
        :root {{
            --bg: #f4f7f9;
            --ink: #17212b;
            --muted: #667085;
            --line: #d9e3ea;
            --surface: #ffffff;
            --soft: #e8f4f1;
            --accent: #0f766e;
            --accent-dark: #115e59;
            --nav: #102a2a;
        }}
        * {{ box-sizing: border-box; }}
        html {{ scroll-behavior: smooth; }}
        body {{
            margin: 0;
            background: var(--bg);
            color: var(--ink);
            font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
        }}
        .app {{ min-height: 100vh; display: grid; grid-template-columns: 280px minmax(0, 1fr); }}
        aside {{
            position: sticky;
            top: 0;
            height: 100vh;
            background: var(--nav);
            color: white;
            padding: 24px 18px;
            overflow-y: auto;
        }}
        .brand-lockup {{ display: grid; grid-template-columns: 44px minmax(0, 1fr); gap: 11px; align-items: center; text-decoration: none; color: white; }}
        .brand-mark {{ width: 44px; height: 44px; border-radius: 8px; display: grid; place-items: center; background: #d7f3ee; color: #0f5f5a; font-size: 19px; font-weight: 950; }}
        .brand-lockup strong, .brand-lockup small {{ display: block; }}
        .brand-lockup strong {{ font-size: 23px; line-height: 1; }}
        .brand-lockup small {{ margin-top: 4px; color: #a7c7c4; font-size: 11px; line-height: 1.25; }}
        .project-chip {{
            margin: 18px 0 22px;
            padding: 12px;
            border: 1px solid rgba(215, 243, 238, 0.25);
            border-radius: 8px;
            background: rgba(215, 243, 238, 0.08);
        }}
        .project-chip > strong, .project-chip > small {{ display: block; overflow-wrap: anywhere; }}
        .project-chip > small {{ color: #a7c7c4; margin-top: 4px; }}
        .sidebar-project-field {{ display: grid; gap: 5px; margin-top: 11px; }}
        .sidebar-project-field span {{ color: #d7f3ee; font-size: 11px; font-weight: 900; text-transform: uppercase; }}
        .project-chip .sidebar-project-field input {{
            width: 100%;
            min-width: 0;
            border: 1px solid rgba(215, 243, 238, 0.28);
            border-radius: 8px;
            padding: 9px 10px;
            background: rgba(255, 255, 255, 0.1);
            color: white;
            font: inherit;
        }}
        .project-chip .sidebar-project-field input::placeholder {{ color: #8fb5b1; }}
        .sidebar-project-actions {{ display: grid; grid-template-columns: 1fr 1fr; gap: 7px; margin-top: 10px; }}
        .sidebar-project-actions button {{ min-width: 0; padding: 0 8px; font-size: 12px; white-space: normal; }}
        .sidebar-project-open {{ width: 100%; margin-top: 7px; font-size: 12px; }}
        .nav-group {{
            margin: 20px 8px 8px;
            color: #8fd8ce;
            font-size: 11px;
            font-weight: 900;
            letter-spacing: 0;
            text-transform: uppercase;
        }}
        nav a {{
            display: block;
            color: #d7f3ee;
            text-decoration: none;
            padding: 13px 14px;
            border-radius: 8px;
            margin: 6px 0;
            font-weight: 700;
        }}
        nav a:hover, nav a.active {{ background: var(--accent); color: white; }}
        nav a.nav-main {{ margin-bottom: 14px; border: 1px solid rgba(215, 243, 238, 0.22); }}
        main {{ padding: 28px; }}
        .hero {{
            display: flex;
            justify-content: space-between;
            gap: 20px;
            align-items: flex-start;
            margin-bottom: 20px;
        }}
        h2 {{ font-size: 28px; margin: 0 0 8px; letter-spacing: 0; }}
        h3 {{ margin: 0 0 8px; font-size: 17px; letter-spacing: 0; }}
        p {{ line-height: 1.5; }}
        .muted, small {{ color: var(--muted); }}
        .eyebrow {{
            color: var(--accent-dark);
            font-size: 12px;
            font-weight: 800;
            letter-spacing: 0;
            margin: 0 0 8px;
            text-transform: uppercase;
        }}
        .grid {{ display: grid; grid-template-columns: repeat(2, minmax(280px, 1fr)); gap: 14px; }}
        .card, .panel {{
            background: var(--surface);
            border: 1px solid var(--line);
            border-radius: 8px;
            padding: 18px;
            box-shadow: 0 10px 28px rgba(21, 40, 55, 0.05);
        }}
        .input-card {{ display: flex; min-height: 236px; flex-direction: column; justify-content: space-between; }}
        .filename {{
            min-height: 24px;
            margin: 8px 0 10px;
            color: var(--ink);
            font-weight: 700;
            overflow-wrap: anywhere;
        }}
        .prep-status {{
            margin: 0 0 10px;
            color: var(--accent-dark);
            font-weight: 800;
            overflow-wrap: anywhere;
        }}
        input[type="text"], .input-card input, .preprocess-row input {{
            width: 100%;
            border: 1px solid var(--line);
            border-radius: 8px;
            padding: 11px 12px;
            font: inherit;
            color: var(--ink);
            background: #fbfdfe;
        }}
        input[type="number"], select, .field input {{
            width: 100%;
            border: 1px solid var(--line);
            border-radius: 8px;
            padding: 11px 12px;
            font: inherit;
            color: var(--ink);
            background: #fbfdfe;
        }}
        textarea {{
            width: 100%;
            min-height: 170px;
            border: 1px solid var(--line);
            border-radius: 8px;
            padding: 11px 12px;
            font: inherit;
            color: var(--ink);
            background: #fbfdfe;
            resize: vertical;
        }}
        button, .button-link, .ghost {{
            appearance: none;
            border: 0;
            border-radius: 8px;
            background: var(--accent);
            color: white;
            cursor: pointer;
            display: inline-flex;
            align-items: center;
            justify-content: center;
            min-height: 38px;
            padding: 0 13px;
            font: inherit;
            font-weight: 800;
            text-decoration: none;
            white-space: nowrap;
        }}
        button:hover, .button-link:hover {{ background: var(--accent-dark); }}
        button:disabled {{
            cursor: not-allowed;
            opacity: 0.55;
        }}
        button.secondary, .ghost {{
            background: var(--soft);
            color: var(--accent-dark);
        }}
        button.danger {{
            background: #fee2e2;
            color: #991b1b;
        }}
        button.danger:hover {{ background: #fecaca; }}
        .disabled {{ opacity: 0.45; pointer-events: none; }}
        .button-row {{ display: flex; flex-wrap: wrap; gap: 8px; margin-top: 10px; }}
        .analysis-action-group {{ margin-top: 20px; }}
        .analysis-action-bar {{
            display: flex;
            align-items: center;
            justify-content: space-between;
            gap: 20px;
            padding: 15px 16px;
            border: 1px solid #9fcfc7;
            border-left: 4px solid var(--accent);
            border-radius: 8px;
            background: #f2faf8;
        }}
        .analysis-action-copy {{ min-width: 220px; }}
        .analysis-action-copy .eyebrow {{ margin-bottom: 3px; }}
        .analysis-action-copy strong, .analysis-action-copy small {{ display: block; }}
        .analysis-action-copy strong {{ font-size: 16px; color: var(--ink); }}
        .analysis-action-copy small {{ margin-top: 3px; color: var(--muted); }}
        .analysis-action-buttons {{ display: flex; justify-content: flex-end; align-items: center; gap: 8px; flex-wrap: wrap; }}
        button.run-primary {{ min-height: 44px; padding: 0 18px; box-shadow: 0 4px 12px rgba(15, 118, 110, 0.18); }}
        button.clear-outputs {{ background: white; color: #9f1239; border: 1px solid #fecdd3; }}
        button.clear-outputs:hover {{ background: #fff1f2; }}
        .analysis-action-status-slot:empty {{ display: none; }}
        .analysis-action-status-slot .analysis-run-status {{ margin: 8px 0 0; }}
        .ready {{ color: var(--accent-dark); font-weight: 800; margin: 0; }}
        .missing {{ color: #98a2b3; font-weight: 800; margin: 0; }}
        .panel {{ margin-top: 18px; scroll-margin-top: 20px; }}
        .panel-heading {{ display: flex; justify-content: space-between; gap: 20px; align-items: flex-start; margin-bottom: 14px; }}
        .panel-heading p {{ margin: 0; color: var(--muted); }}
        .path-row {{
            display: grid;
            grid-template-columns: 74px minmax(0, 1fr) auto;
            gap: 12px;
            align-items: center;
            padding: 12px 0;
            border-top: 1px solid var(--line);
        }}
        .path-row:first-child {{ border-top: 0; }}
        .path-row strong {{ display: block; overflow-wrap: anywhere; }}
        .path-row small {{ display: block; margin-top: 3px; overflow-wrap: anywhere; }}
        .path-row a {{ color: var(--accent-dark); font-weight: 800; text-decoration: none; }}
        .badge {{
            display: inline-flex;
            justify-content: center;
            border-radius: 999px;
            padding: 6px 8px;
            background: var(--soft);
            color: var(--accent-dark);
            font-weight: 900;
            font-size: 12px;
        }}
        .module-card {{
            display: flex;
            justify-content: space-between;
            gap: 14px;
            min-height: 126px;
            color: inherit;
            text-decoration: none;
            background: var(--surface);
            border: 1px solid var(--line);
            border-radius: 8px;
            padding: 18px;
            box-shadow: 0 10px 28px rgba(21, 40, 55, 0.05);
        }}
        .module-card:hover {{ border-color: var(--accent); }}
        .module-card p {{ margin: 0; color: var(--muted); }}
        .module-card span {{ color: var(--accent-dark); font-weight: 900; white-space: nowrap; }}
        .main-menu-hero h2 {{ font-size: 38px; }}
        .brand-expansion {{ margin: 0; max-width: 680px; color: var(--ink); font-size: 18px; font-weight: 750; }}
        .tool-menu {{ margin-top: 28px; }}
        .section-heading {{ margin-bottom: 18px; }}
        .tool-menu-group {{ padding: 20px 0; border-top: 1px solid var(--line); }}
        .tool-menu-group:first-of-type {{ border-top: 0; }}
        .tool-menu-heading {{ display: flex; justify-content: space-between; align-items: baseline; gap: 16px; margin-bottom: 12px; }}
        .tool-menu-heading h2 {{ font-size: 20px; margin: 0; }}
        .tool-menu-heading > span {{ color: var(--muted); font-size: 12px; font-weight: 800; }}
        .tool-menu-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 12px; }}
        .project-row {{
            display: grid;
            grid-template-columns: minmax(220px, 0.7fr) minmax(320px, 1fr) auto;
            gap: 10px;
            align-items: center;
        }}
        .project-row input {{
            width: 100%;
            border: 1px solid var(--line);
            border-radius: 8px;
            padding: 11px 12px;
            font: inherit;
            background: #fbfdfe;
        }}
        .dataset-library-grid {{
            display: grid;
            grid-template-columns: repeat(3, minmax(0, 1fr));
            gap: 12px;
        }}
        .dataset-library-card {{
            display: grid;
            align-content: start;
            gap: 12px;
            padding: 14px;
            border: 1px solid var(--line);
            border-radius: 8px;
            background: #fbfdfe;
        }}
        .dataset-library-card p {{ margin: 0; }}
        .dataset-options {{
            display: grid;
            grid-template-columns: repeat(2, minmax(0, 1fr));
            gap: 10px;
        }}
        .dataset-option {{
            display: grid;
            grid-template-columns: 20px minmax(0, 1fr);
            gap: 10px;
            align-items: start;
            min-height: 82px;
            padding: 13px;
            border: 1px solid var(--line);
            border-radius: 8px;
            background: #fbfdfe;
            cursor: pointer;
        }}
        .dataset-option.selected {{ border-color: var(--accent); background: var(--soft); }}
        .dataset-option.disabled-option {{ cursor: not-allowed; opacity: 0.62; }}
        .dataset-option input {{ width: 18px; height: 18px; margin: 2px 0 0; accent-color: var(--accent); }}
        .dataset-option strong, .dataset-option small {{ display: block; overflow-wrap: anywhere; }}
        .dataset-option small {{ margin-top: 4px; }}
        .dataset-option.batch-option {{
            grid-template-columns: minmax(0, 1fr) auto;
            align-items: center;
            cursor: default;
        }}
        .dataset-source-grid {{
            display: grid;
            grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);
            gap: 20px;
        }}
        .dataset-source-grid h3 {{ margin-bottom: 10px; }}
        .dataset-options.stacked {{ grid-template-columns: 1fr; }}
        .compact-empty {{ padding: 13px; }}
        .selected-samples {{ margin-top: 20px; border-top: 1px solid var(--line); padding-top: 18px; }}
        .sample-heading {{ margin-bottom: 10px; }}
        .tool-sample-list {{ display: grid; gap: 8px; }}
        .tool-sample-row {{
            display: grid;
            grid-template-columns: 76px minmax(180px, 1fr) minmax(150px, 220px) 118px 38px;
            gap: 12px;
            align-items: end;
            padding: 11px 12px;
            border: 1px solid var(--line);
            border-radius: 8px;
            background: #fbfdfe;
        }}
        .tool-sample-row.is-excluded {{ background: #f4f6f7; color: var(--muted); }}
        .sample-include {{ display: grid; justify-items: start; gap: 5px; align-self: center; font-size: 12px; font-weight: 900; color: var(--ink); }}
        .sample-include input {{ width: 20px; height: 20px; margin: 0; accent-color: var(--accent); }}
        .sample-file strong, .sample-file small, .replicate-value small, .replicate-value strong {{ display: block; overflow-wrap: anywhere; }}
        .sample-file small, .replicate-value small {{ margin-top: 3px; color: var(--muted); }}
        .replicate-value strong {{ margin-top: 7px; }}
        .compact-field {{ gap: 4px; }}
        .settings-grid {{
            display: grid;
            grid-template-columns: repeat(3, minmax(150px, 1fr));
            gap: 12px;
        }}
        .external-tool-list {{ display: grid; gap: 8px; }}
        .external-tool-row {{
            display: grid;
            grid-template-columns: 74px minmax(0, 1fr) auto;
            gap: 12px;
            align-items: center;
            padding: 13px;
            border: 1px solid var(--line);
            border-radius: 8px;
            background: #fbfdfe;
        }}
        .external-tool-row > div > strong, .external-tool-row small {{ display: block; overflow-wrap: anywhere; }}
        .external-tool-row small {{ margin-top: 3px; }}
        .tool-paths {{ margin-top: 7px; }}
        .badge.missing {{ background: #fff1f2; color: #9f1239; }}
        .two-col {{ display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 22px; }}
        .top-actions {{ justify-content: flex-end; margin-top: 0; }}
        .toggle-row {{
            display: flex;
            align-items: center;
            gap: 10px;
            color: var(--ink);
            font-weight: 800;
        }}
        .toggle-row input {{ width: 18px; height: 18px; accent-color: var(--accent); }}
        .preprocess-row {{
            display: grid;
            grid-template-columns: minmax(220px, 0.9fr) minmax(260px, 1fr);
            gap: 12px;
            align-items: center;
            padding: 14px 0;
            border-top: 1px solid var(--line);
        }}
        .preprocess-row:first-of-type {{ border-top: 0; }}
        .preprocess-row strong, .preprocess-row small {{ display: block; margin-top: 4px; overflow-wrap: anywhere; }}
        .preprocess-row .button-row {{ grid-column: 1 / -1; margin-top: 0; }}
        .notice {{
            padding: 12px 14px;
            border-radius: 8px;
            background: #fff8e6;
            color: #8a5a00;
            margin: 12px 0 18px;
            border: 1px solid #f2d996;
        }}
        .notice.status-info {{ background: #eef6ff; color: #174f84; border-color: #b9d7f2; }}
        .notice.status-success {{ background: #ecfdf3; color: #166534; border-color: #86d6a5; }}
        .notice.status-warn {{ background: #fff8e6; color: #8a5a00; border-color: #f2d996; }}
        .notice.status-error {{ background: #fff1f2; color: #a1122f; border-color: #f0a3b2; }}
        .trim-panel code {{
            padding: 2px 5px;
            border-radius: 5px;
            background: #eef5f6;
            color: var(--accent-dark);
        }}
        .trim-controls {{
            display: grid;
            grid-template-columns: minmax(180px, 1fr) repeat(5, minmax(120px, 0.65fr));
            gap: 12px;
            align-items: end;
            margin: 16px 0 20px;
        }}
        .cutadapt-controls {{
            display: none;
            grid-template-columns: repeat(4, minmax(130px, 1fr));
            gap: 12px;
            align-items: end;
            margin: -4px 0 20px;
            padding: 14px;
            border: 1px solid var(--line);
            border-radius: 8px;
            background: #fbfdfe;
        }}
        body.trim-cutadapt .cutadapt-controls {{ display: grid; }}
        .cutadapt-toggle {{
            min-height: 42px;
            align-self: end;
            padding-bottom: 2px;
        }}
        .field {{
            display: grid;
            gap: 6px;
            color: var(--ink);
            font-weight: 800;
        }}
        .field span {{
            font-size: 12px;
            color: var(--muted);
            font-weight: 900;
            text-transform: uppercase;
        }}
        .trim-sheet-head {{
            display: flex;
            justify-content: space-between;
            gap: 18px;
            align-items: flex-start;
            margin-top: 18px;
        }}
        .trim-sheet-head p {{ margin: 0; }}
        .trim-rows {{
            display: grid;
            gap: 10px;
            margin-top: 12px;
        }}
        .trim-row {{
            display: grid;
            grid-template-columns: minmax(220px, 1fr) minmax(220px, 1fr) 38px;
            gap: 10px;
            align-items: center;
            padding: 12px;
            border: 1px solid var(--line);
            border-radius: 8px;
            background: #fbfdfe;
        }}
        .trim-row input {{
            width: 100%;
            border: 1px solid var(--line);
            border-radius: 8px;
            padding: 10px 11px;
            font: inherit;
        }}
        .degradome-rows {{
            display: grid;
            gap: 10px;
        }}
        .degradome-row {{
            display: grid;
            grid-template-columns: minmax(130px, 0.55fr) minmax(130px, 0.55fr) minmax(90px, 0.35fr) minmax(260px, 1fr) 38px;
            gap: 10px;
            align-items: center;
            padding: 12px;
            border: 1px solid var(--line);
            border-radius: 8px;
            background: #fbfdfe;
        }}
        .degradome-row input {{
            width: 100%;
            border: 1px solid var(--line);
            border-radius: 8px;
            padding: 10px 11px;
            font: inherit;
        }}
        .file-picker {{
            display: grid;
            grid-template-columns: minmax(0, 1fr) auto;
            gap: 8px;
        }}
        .icon-button {{
            width: 38px;
            min-height: 38px;
            padding: 0;
        }}
        .inline-check {{
            display: inline-flex;
            align-items: center;
            gap: 8px;
            min-height: 38px;
            font-weight: 800;
            color: var(--ink);
        }}
        .inline-check input {{ width: 18px; height: 18px; accent-color: var(--accent); }}
        .rep-preview {{
            grid-column: 1 / -1;
            margin-top: -4px;
        }}
        .trim-actions {{
            align-items: end;
            margin-top: 16px;
        }}
        .sample-output-table {{
            display: grid;
            gap: 8px;
        }}
        .sample-output-row {{
            display: grid;
            grid-template-columns: 150px minmax(180px, 0.7fr) minmax(260px, 1fr) auto;
            gap: 12px;
            align-items: center;
            padding: 12px;
            border: 1px solid var(--line);
            border-radius: 8px;
            background: #fbfdfe;
        }}
        .sample-output-row strong, .sample-output-row small {{
            display: block;
            overflow-wrap: anywhere;
        }}
        .sample-output-row a {{
            color: var(--accent-dark);
            font-weight: 800;
            text-decoration: none;
        }}
        .tool-grid {{
            display: grid;
            grid-template-columns: repeat(6, minmax(120px, 1fr));
            gap: 12px;
            align-items: end;
            margin-bottom: 12px;
        }}
        .plot-grid {{
            display: grid;
            grid-template-columns: repeat(2, minmax(280px, 1fr));
            gap: 14px;
        }}
        .plot-card {{
            display: grid;
            gap: 8px;
            padding: 10px;
            border: 1px solid var(--line);
            border-radius: 8px;
            color: var(--ink);
            text-decoration: none;
            background: #fbfdfe;
        }}
        .plot-card img {{
            width: 100%;
            border: 1px solid var(--line);
            border-radius: 6px;
            background: white;
        }}
        .empty-state {{
            padding: 18px;
            border: 1px dashed var(--line);
            border-radius: 8px;
            color: var(--muted);
            background: #fbfdfe;
        }}
        .compact-select {{
            width: min(220px, 100%);
        }}
        body.trim-single .paired-only {{
            display: none;
        }}
        .process-terminal {{
            position: fixed;
            right: 22px;
            bottom: 18px;
            z-index: 20;
            width: min(560px, calc(100vw - 44px));
            overflow: hidden;
            border: 1px solid rgba(15, 118, 110, 0.32);
            border-radius: 8px;
            background: rgba(12, 30, 32, 0.94);
            color: #d7f3ee;
            box-shadow: 0 18px 46px rgba(16, 42, 42, 0.24);
            backdrop-filter: blur(10px);
        }}
        .terminal-head {{
            display: flex;
            align-items: center;
            justify-content: space-between;
            gap: 12px;
            min-height: 38px;
            padding: 0 10px 0 12px;
            border-bottom: 1px solid rgba(215, 243, 238, 0.14);
        }}
        .terminal-title {{
            display: flex;
            align-items: center;
            gap: 8px;
            font-size: 12px;
            font-weight: 900;
            text-transform: uppercase;
        }}
        .terminal-dot {{
            width: 8px;
            height: 8px;
            border-radius: 999px;
            background: #22c55e;
            box-shadow: 0 0 0 4px rgba(34, 197, 94, 0.13);
        }}
        .terminal-actions {{ display: flex; gap: 6px; }}
        .terminal-actions button {{
            min-height: 26px;
            padding: 0 9px;
            border-radius: 6px;
            background: rgba(215, 243, 238, 0.12);
            color: #d7f3ee;
            font-size: 12px;
        }}
        .terminal-actions button:hover {{ background: rgba(215, 243, 238, 0.22); }}
        .terminal-actions button.stop-button {{
            background: rgba(220, 38, 38, 0.82);
            color: white;
        }}
        .terminal-actions button.stop-button:hover {{ background: rgba(185, 28, 28, 0.95); }}
        .terminal-log {{
            max-height: 128px;
            overflow: auto;
            padding: 8px 12px 10px;
            font-family: "SFMono-Regular", Consolas, "Liberation Mono", monospace;
            font-size: 12px;
            line-height: 1.45;
        }}
        .terminal-line {{
            display: grid;
            grid-template-columns: 70px minmax(0, 1fr);
            gap: 8px;
            padding: 2px 0;
            overflow-wrap: anywhere;
        }}
        .terminal-line span {{ color: #8fd8ce; }}
        .terminal-line.warn {{ color: #fde68a; }}
        .terminal-line.error {{ color: #fecaca; }}
        .process-terminal.collapsed .terminal-log {{ display: none; }}
        .process-terminal.collapsed {{
            width: min(360px, calc(100vw - 44px));
        }}
        @media (max-width: 980px) {{
            .app {{ grid-template-columns: 1fr; }}
            aside {{ position: static; height: auto; padding: 16px; overflow: hidden; }}
            aside .project-chip {{ margin: 12px 0; }}
            nav {{ display: flex; gap: 6px; overflow-x: auto; padding-bottom: 4px; }}
            nav .nav-group {{ display: none; }}
            nav a, nav a.nav-main {{ flex: 0 0 auto; margin: 0; padding: 10px 12px; }}
            .grid, .two-col {{ grid-template-columns: 1fr; }}
            .hero, .panel-heading {{ flex-direction: column; }}
            .path-row {{ grid-template-columns: 56px minmax(0, 1fr); }}
            .path-row a {{ grid-column: 2; }}
            .preprocess-row {{ grid-template-columns: 1fr; }}
            .preprocess-row .button-row {{ grid-column: auto; }}
            .project-row {{ grid-template-columns: 1fr; }}
            .trim-controls, .cutadapt-controls, .trim-row {{ grid-template-columns: 1fr; }}
            .degradome-row {{ grid-template-columns: 1fr; }}
            .sample-output-row {{ grid-template-columns: 1fr; }}
            .tool-grid, .plot-grid, .dataset-library-grid, .dataset-options, .dataset-source-grid, .settings-grid, .external-tool-row, .tool-menu-grid {{ grid-template-columns: 1fr; }}
            .tool-sample-row {{ grid-template-columns: 76px minmax(0, 1fr) 90px 38px; }}
            .tool-sample-row .sample-file {{ grid-column: 2 / 4; }}
            .tool-sample-row .compact-field {{ grid-column: 1 / 3; }}
            .trim-sheet-head {{ flex-direction: column; }}
            .analysis-action-bar {{ flex-direction: column; align-items: stretch; }}
            .analysis-action-buttons {{ justify-content: stretch; }}
            .analysis-action-buttons > * {{ flex: 1 1 auto; text-align: center; }}
            .process-terminal {{
                right: 12px;
                bottom: 12px;
                width: calc(100vw - 24px);
            }}
        }}
        [hidden] {{ display: none !important; }}
        :focus-visible {{ outline: 3px solid #b45309; outline-offset: 3px; }}
        .skip-link {{ position: fixed; left: 12px; top: -80px; z-index: 100; padding: 12px 18px; background: white; color: #17212b; }}
        .skip-link:focus {{ top: 12px; }}
        body {{ line-height: 1.5; }}
        main {{ min-width: 0; padding-bottom: 100px; }}
        h2, h3, strong, button, a {{ overflow-wrap: anywhere; }}
        .hero h2 {{ font-size: 28px; line-height: 1.2; }}
        .panel {{ border: 0; border-top: 1px solid var(--line); border-radius: 0; box-shadow: none; background: transparent; padding: 24px 0; }}
        .panel h3 {{ font-size: 19px; scroll-margin-top: 24px; }}
        input, select, textarea {{ max-width: 100%; }}
        input:not([type=checkbox]):not([type=radio]), select, button, .button-link, .ghost {{ min-height: 42px; }}
        input[type=checkbox], input[type=radio] {{ width: 18px; height: 18px; accent-color: var(--accent); }}
        button:disabled {{ opacity: .55; cursor: not-allowed; }}
        .run-primary {{ font-weight: 750; border: 2px solid var(--accent-dark); }}
        .analysis-action-bar {{ background: #eaf5ef; border-color: #98bfae; box-shadow: none; }}
        .analysis-action-group[data-running=true] .analysis-action-bar {{ border-color: #0f766e; }}
        .analysis-run-status {{ min-height: 48px; font-weight: 600; overflow-wrap: anywhere; }}
        .status-error {{ border-left: 5px solid #b91c1c; }}
        .status-success {{ border-left: 5px solid #15803d; }}
        .status-warn {{ border-left: 5px solid #b45309; }}
        .tool-search-label {{ display: block; margin: 16px 0 6px; font-size: 13px; color: #e2efed; }}
        #tool-search {{ width: 100%; padding: 10px; border-radius: 6px; border: 1px solid #99b8b1; background: white; color: #17212b; font: inherit; }}
        #page-sections {{ display: flex; gap: 6px 16px; flex-wrap: wrap; margin-bottom: 20px; border-bottom: 1px solid var(--line); padding-bottom: 12px; }}
        #page-sections:empty {{ display: none; }}
        #page-sections a {{ color: #115e59; background: transparent; padding: 4px 0; margin: 0; font-size: 13px; text-decoration: underline; text-underline-offset: 4px; }}
        .plot-card {{ padding: 12px; box-shadow: none; }}
        .plot-card strong {{ font-size: 13px; font-weight: 600; }}
        .table-scroll {{ max-width: 100%; overflow: auto; }}
        .table-scroll th {{ position: sticky; top: 0; background: #edf2f4; }}
        .table-scroll tbody tr:nth-child(even) {{ background: #f0f5f5; }}
        @media (prefers-reduced-motion: reduce) {{ html {{ scroll-behavior: auto; }} * {{ animation: none !important; transition: none !important; }} }}
        @media (max-width: 760px) {{
            main {{ padding: 18px 14px 110px; }}
            .hero h2 {{ font-size: 24px; }}
            .analysis-action-copy {{ min-width: 0; }}
            .terminal-head {{ flex-wrap: wrap; padding: 8px; }}
        }}
    </style>
</head>
<body data-project="{esc(str(shell_project_root))}">
    <a class="skip-link" href="#main-content">Skip to tool</a>
    <div class="app">
        <aside>
            <a class="brand-lockup" href="/" aria-label="INCI main menu">
                <span class="brand-mark">IN</span>
                <span><strong>INCI</strong><small>Insecticidal Cross-Kingdom RNAi</small></span>
            </a>
            <div class="project-chip">
                <strong>Project</strong>
                <small>{esc(shell_project_root)}</small>
                <label class="sidebar-project-field">
                    <span>Project name</span>
                    <input id="project-name" value="{esc(shell_project_name if shell_project_name != 'Untitled project' else '')}" placeholder="Project name">
                </label>
                <label class="sidebar-project-field">
                    <span>Project directory</span>
                    <input id="project-dir" value="{esc(shell_project_dir)}" placeholder="Output directory">
                </label>
                <div class="sidebar-project-actions">
                    <button type="button" class="secondary" onclick="browseProjectDir()">Browse</button>
                    <button type="button" onclick="saveProjectName()">Save Project</button>
                </div>
                <a class="ghost sidebar-project-open" href="/reveal?path={url_for(shell_project_root)}">Open Project Outputs</a>
            </div>
            <label class="tool-search-label" for="tool-search">Find a tool</label>
            <input id="tool-search" type="search" placeholder="Search tools" aria-controls="tool-navigation">
            <nav id="tool-navigation" aria-label="Pipeline tools">
                {render_nav(active)}
            </nav>
        </aside>
        <main id="main-content" tabindex="-1">
            <nav id="page-sections" aria-label="On this page"></nav>
            {content}
        </main>
    </div>
    <section class="process-terminal" id="process-terminal" aria-label="Process console">
        <div class="terminal-head">
            <div class="terminal-title"><span class="terminal-dot"></span><span>Process Console</span></div>
            <div class="terminal-actions">
                <button type="button" class="stop-button" onclick="stopPipelineProcesses()">Stop All Runs</button>
                <button type="button" onclick="clearTerminal()">Clear</button>
                <button type="button" onclick="toggleTerminal()" aria-controls="terminal-log">Show / hide log</button>
            </div>
        </div>
        <div class="terminal-log" id="terminal-log"></div>
    </section>
    {render_scripts()}
</body>
</html>"""


def render_page(state: dict[str, Any], port: int) -> str:
    content = f"""
        <section class="hero main-menu-hero">
            <div>
                <p class="eyebrow">Analysis workspace</p>
                <h2>INCI</h2>
                <p class="brand-expansion">Insecticidal Cross-Kingdom RNA Interference Pipeline</p>
            </div>
            <a class="button-link" href="/reveal?path={url_for(project_output_root(state))}">Open Project Folder</a>
        </section>
        <div class="notice" id="status">Choose a tool below. The active project and its output directory are shared across the pipeline.</div>

        <section class="tool-menu">
            <div class="section-heading"><p class="eyebrow">Tools</p><h2>Pipeline Menu</h2></div>
            {render_modules(state)}
        </section>
    """
    return render_shell("INCI RNAi Pipeline", "home", content, port)


def render_module_page(module: ModuleSpec, state: dict[str, Any], port: int) -> str:
    content = render_module_content(module, state)
    return render_shell(module.title, module.key, content, port)


class InciRequestHandler(BaseHTTPRequestHandler):
    server_version = "INCI/0.1"

    def do_GET(self) -> None:
        parsed = urllib.parse.urlparse(self.path)
        query = urllib.parse.parse_qs(parsed.query)

        if parsed.path == "/":
            state = load_state()
            self.respond_html(render_page(state, self.server.server_port))
            return

        if parsed.path == "/favicon.svg":
            body = FAVICON_SVG.encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "image/svg+xml; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "public, max-age=86400")
            self.end_headers()
            self.wfile.write(body)
            return

        if parsed.path == "/settings":
            state = load_state()
            self.respond_html(render_shell("Settings", "settings", render_settings_content(state), self.server.server_port))
            return

        if parsed.path.startswith("/module/"):
            module_key = parsed.path.removeprefix("/module/").strip("/")
            module = module_for_key(module_key)
            if module is None:
                self.respond_text("Pipeline section not found", status=404)
                return
            self.respond_html(render_module_page(module, load_state(), self.server.server_port))
            return

        if parsed.path == "/browse":
            key = query.get("key", [""])[0]
            if key not in {spec.key for spec in INPUTS}:
                self.respond_json({"path": "", "message": "Unknown input."}, status=400)
                return

            chosen = choose_local_file(f"Choose {input_for_key(key).label} file")
            if chosen:
                state = load_state()
                state["paths"][key] = chosen
                save_state(state)
                self.respond_json({"path": chosen})
            else:
                self.respond_json({"path": "", "message": "No file selected. You can paste a path manually."})
            return

        if parsed.path == "/browse-preprocessed":
            key = query.get("key", [""])[0]
            if key not in {spec.key for spec in INPUTS}:
                self.respond_json({"path": "", "message": "Unknown input."}, status=400)
                return

            chosen = choose_local_file(f"Choose preprocessed {input_for_key(key).label} file")
            if chosen:
                state = load_state()
                state["preprocessed_paths"][key] = chosen
                state["preprocessing_status"][key] = "ready"
                save_state(state)
                self.respond_json({"path": chosen})
            else:
                self.respond_json({"path": "", "message": "No file selected. You can paste a path manually."})
            return

        if parsed.path == "/browse-tool-sample":
            module_key = query.get("module_key", [""])[0]
            input_key = query.get("input_key", [""])[0]
            if not valid_module_dataset(module_key, input_key):
                self.respond_json({"path": "", "message": "Unknown tool input."}, status=400)
                return
            label, _ = sequencing_dataset_spec(input_key)
            chosen = choose_local_file(f"Choose preprocessed {label}")
            if chosen:
                self.respond_json({"path": chosen})
            else:
                self.respond_json({"path": "", "message": "No file selected."})
            return

        if parsed.path == "/browse-dsrna-reference":
            module_key = query.get("module_key", [""])[0]
            if module_key not in {"dsrna-identification", "dsrna-plotter", "srna-dsrna-identification"}:
                self.respond_json({"path": "", "message": "Unknown dsRNA tool."}, status=400)
                return
            chosen = choose_local_file("Choose reference FASTA")
            if chosen:
                state = load_state()
                state.setdefault("module_paths", {}).setdefault(module_key, {})["transcriptome"] = chosen
                save_state(state)
                self.respond_json({"path": chosen})
            else:
                self.respond_json({"path": "", "message": "No FASTA selected."})
            return

        if parsed.path == "/browse-fasta-dedup-input":
            chosen = choose_local_file("Choose nucleotide FASTA for CD-HIT-EST deduplication")
            if chosen:
                state = load_state()
                state.setdefault("module_paths", {}).setdefault("fasta-deduplication", {})["transcriptome"] = chosen
                save_state(state)
                self.respond_json({"path": chosen})
            else:
                self.respond_json({"path": "", "message": "No FASTA selected."})
            return

        if parsed.path == "/browse-trim-file":
            chosen = choose_local_file("Choose FASTQ read file")
            if chosen:
                self.respond_json({"path": chosen})
            else:
                self.respond_json({"path": "", "message": "No file selected. You can paste a path manually."})
            return

        if parsed.path == "/browse-project-dir":
            chosen = choose_local_dir("Choose project output folder")
            if chosen:
                self.respond_json({"path": chosen})
            else:
                self.respond_json({"path": "", "message": "No folder selected. You can paste a path manually."})
            return

        if parsed.path == "/browse-srna-reference":
            chosen = choose_local_file("Choose sRNA reference FASTA")
            if chosen:
                self.respond_json({"path": chosen})
            else:
                self.respond_json({"path": "", "message": "No file selected. You can paste a path manually."})
            return

        if parsed.path == "/browse-sirna-annotations":
            module_key = query.get("module_key", [""])[0]
            if module_key not in {"dsrna-identification", "srna-mapping"}:
                self.respond_json({"path": "", "message": "Unknown siRNA annotation tool."}, status=400)
                return
            chosen = choose_local_file("Choose siRNA annotation FASTA")
            if chosen:
                state = load_state()
                remember_sirna_annotation_path(state, module_key, chosen)
                self.respond_json({"path": chosen})
            else:
                self.respond_json({"path": "", "message": "No FASTA selected. You can paste a path manually."})
            return

        if parsed.path == "/browse-srna-control-file":
            chosen = choose_local_file("Choose sRNA/control FASTA")
            if chosen:
                self.respond_json({"path": chosen})
            else:
                self.respond_json({"path": "", "message": "No file selected. You can paste a path manually."})
            return

        if parsed.path == "/browse-degradome-file":
            chosen = choose_local_file("Choose degradome FASTA/FASTQ")
            if chosen:
                self.respond_json({"path": chosen, "name": Path(chosen).stem})
            else:
                self.respond_json({"path": "", "message": "No file selected. You can paste a path manually."})
            return

        if parsed.path == "/browse-degradome-input":
            chosen = choose_local_file("Choose FASTA input")
            if chosen:
                self.respond_json({"path": chosen})
            else:
                self.respond_json({"path": "", "message": "No file selected. You can paste a path manually."})
            return

        if parsed.path == "/output-file":
            target = Path(query.get("path", [""])[0]).expanduser()
            active_root = project_output_root(load_state())
            try:
                target.relative_to(active_root)
            except ValueError:
                self.respond_text("Forbidden", status=403)
                return
            if not target.exists() or not target.is_file():
                self.respond_text("Not found", status=404)
                return
            content_type = {".png": "image/png", ".svg": "image/svg+xml"}.get(target.suffix.lower(), "application/octet-stream")
            self.respond_bytes(target.read_bytes(), content_type, 200)
            return

        if parsed.path == "/job-status":
            job_id = query.get("job_id", [""])[0]
            snapshot = pipeline_job_snapshot(job_id)
            if snapshot is None:
                self.respond_json({"ok": False, "message": "Job not found."}, status=404)
                return
            snapshot["running_processes"] = running_process_labels()
            self.respond_json({"ok": True, "job": snapshot})
            return

        if parsed.path == "/reveal":
            target = query.get("path", [""])[0]
            reveal_path(target or OUTPUT_DIR)
            back_url = self.headers.get("Referer", "/")
            if not back_url.startswith(f"http://{HOST}:"):
                back_url = "/"
            self.send_response(303)
            self.send_header("Location", back_url)
            self.end_headers()
            return

        self.respond_text("Not found", status=404)

    def do_POST(self) -> None:
        payload = self.read_json()

        if self.path == "/stop-running":
            result = stop_running_processes()
            self.respond_json({"ok": True, **result})
            return

        if self.path == "/save-project":
            name = str(payload.get("name", "")).strip()
            project_dir = str(payload.get("project_dir", "")).strip()
            if not name:
                self.respond_json({"ok": False, "message": "Project name is required."}, status=400)
                return
            state = load_state()
            state["project"] = {"name": name, "dir": project_dir}
            project_output_root(state).mkdir(parents=True, exist_ok=True)
            save_state(state)
            self.respond_json({"ok": True, "message": f"Project saved: {name}", "root": str(project_output_root(state))})
            return

        if self.path == "/save-settings":
            try:
                settings = normalized_plot_settings(payload.get("plot_settings", {}))
            except ValueError as exc:
                self.respond_json({"ok": False, "message": str(exc)}, status=400)
                return
            state = load_state()
            state["plot_settings"] = settings
            save_state(state)
            apply_global_plot_settings(state)
            self.respond_json({"ok": True, "message": "Global plot settings saved.", "plot_settings": settings})
            return

        if self.path == "/install-tools":
            requested = payload.get("tool_keys", [])
            if not isinstance(requested, list):
                self.respond_json({"ok": False, "message": "Choose one or more external tools."}, status=400)
                return
            valid_keys = {tool.key for tool in EXTERNAL_TOOLS}
            tool_keys = [str(key) for key in requested if str(key) in valid_keys]
            try:
                job_id = start_external_tools_install_job(tool_keys)
            except ValueError as exc:
                self.respond_json({"ok": False, "message": str(exc)}, status=400)
                return
            except Exception as exc:
                self.respond_json({"ok": False, "message": f"Could not start installation: {exc}"}, status=500)
                return
            self.respond_json({"ok": True, "message": "External tool installation started.", "job_id": job_id})
            return

        if self.path == "/reset-page":
            module_key = str(payload.get("module_key", "")).strip()
            try:
                result = reset_output_paths_for_page(module_key, load_state())
            except ValueError as exc:
                self.respond_json({"ok": False, "message": str(exc)}, status=400)
                return
            except Exception as exc:
                self.respond_json({"ok": False, "message": f"Reset failed: {exc}"}, status=500)
                return
            count = len(result.get("removed", []))
            self.respond_json({"ok": True, "message": f"Reset {count} output path(s) for this page.", **result})
            return

        if self.path == "/save":
            key = str(payload.get("key", ""))
            path = str(payload.get("path", "")).strip()
            if key not in {spec.key for spec in INPUTS}:
                self.respond_json({"ok": False, "message": "Unknown input."}, status=400)
                return

            state = load_state()
            if path:
                state["paths"][key] = path
            else:
                state["paths"].pop(key, None)
            save_state(state)
            self.respond_json({"ok": True})
            return

        if self.path == "/save-preprocessed":
            key = str(payload.get("key", ""))
            path = str(payload.get("path", "")).strip()
            if key not in {spec.key for spec in INPUTS}:
                self.respond_json({"ok": False, "message": "Unknown input."}, status=400)
                return

            state = load_state()
            if path:
                state["preprocessed_paths"][key] = path
                state["preprocessing_status"][key] = "ready"
            else:
                state["preprocessed_paths"].pop(key, None)
                state["preprocessing_status"].pop(key, None)
            save_state(state)
            self.respond_json({"ok": True})
            return

        if self.path == "/add-tool-batch":
            module_key = str(payload.get("module_key", ""))
            input_key = str(payload.get("input_key", ""))
            try:
                state = load_state()
                added = add_batch_to_tool_samples(state, module_key, input_key, str(payload.get("batch_id", "")))
                save_state(state)
            except ValueError as exc:
                self.respond_json({"ok": False, "message": str(exc)}, status=400)
                return
            message = f"Added {added} sample(s)." if added else "All samples from that batch were already selected."
            self.respond_json({"ok": True, "message": message, "added": added})
            return

        if self.path == "/add-tool-sample":
            module_key = str(payload.get("module_key", ""))
            input_key = str(payload.get("input_key", ""))
            try:
                state = load_state()
                added = add_local_tool_sample(
                    state,
                    module_key,
                    input_key,
                    str(payload.get("read1", "")).strip(),
                    str(payload.get("read2", "")).strip(),
                )
                save_state(state)
            except ValueError as exc:
                self.respond_json({"ok": False, "message": str(exc)}, status=400)
                return
            self.respond_json({"ok": True, "message": "Sample added.", "added": added})
            return

        if self.path == "/update-tool-samples":
            module_key = str(payload.get("module_key", ""))
            input_key = str(payload.get("input_key", ""))
            sample_settings = payload.get("samples", payload.get("groups", {}))
            if not isinstance(sample_settings, dict):
                self.respond_json({"ok": False, "message": "Invalid analysis sample set."}, status=400)
                return
            try:
                state = load_state()
                samples = update_tool_sample_settings(state, module_key, input_key, sample_settings)
                save_state(state)
            except ValueError as exc:
                self.respond_json({"ok": False, "message": str(exc)}, status=400)
                return
            included = sum(1 for sample in samples if sample.get("included", True))
            self.respond_json({"ok": True, "message": "Analysis sample set saved.", "sample_count": len(samples), "included_count": included})
            return

        if self.path == "/remove-tool-sample":
            module_key = str(payload.get("module_key", ""))
            input_key = str(payload.get("input_key", ""))
            try:
                state = load_state()
                removed = remove_tool_sample(state, module_key, input_key, str(payload.get("sample_id", "")))
                save_state(state)
            except ValueError as exc:
                self.respond_json({"ok": False, "message": str(exc)}, status=400)
                return
            if not removed:
                self.respond_json({"ok": False, "message": "Sample was not found."}, status=404)
                return
            self.respond_json({"ok": True, "message": "Sample removed."})
            return

        if self.path == "/save-dsrna-reference":
            module_key = str(payload.get("module_key", ""))
            try:
                path = save_dsrna_reference(
                    load_state(),
                    module_key,
                    str(payload.get("reference_text", "")),
                    str(payload.get("reference_fasta", "")),
                )
            except ValueError as exc:
                self.respond_json({"ok": False, "message": str(exc)}, status=400)
                return
            except Exception as exc:
                self.respond_json({"ok": False, "message": f"Could not save reference: {exc}"}, status=500)
                return
            self.respond_json({"ok": True, "message": "Reference sequence saved.", "path": path})
            return

        if self.path == "/save-trimming":
            try:
                settings = normalize_trimming_payload(payload)
            except ValueError as exc:
                self.respond_json({"ok": False, "message": str(exc)}, status=400)
                return

            state = load_state()
            save_trimming_settings(state, settings)
            self.respond_json({"ok": True, "message": f"Saved {len(settings['samples'])} trimming sample(s)."})
            return

        if self.path == "/run-trimming":
            try:
                settings = normalize_trimming_payload(payload)
                state = load_state()
                save_trimming_settings(state, settings)
                job_id = start_trimming_job(settings)
            except ValueError as exc:
                self.respond_json({"ok": False, "message": str(exc)}, status=400)
                return
            except Exception as exc:
                self.respond_json({"ok": False, "message": f"Could not start trimming: {exc}"}, status=500)
                return

            self.respond_json(
                {
                    "ok": True,
                    "message": f"{trimming_tool_title(settings['tool'])} batch started. Progress will appear in the process console.",
                    "job_id": job_id,
                }
            )
            return

        if self.path == "/install-trimgalore":
            try:
                job_id = start_trim_galore_install_job()
            except ValueError as exc:
                self.respond_json({"ok": False, "message": str(exc)}, status=400)
                return
            except Exception as exc:
                self.respond_json({"ok": False, "message": f"Could not start Trim Galore installation: {exc}"}, status=500)
                return

            self.respond_json(
                {
                    "ok": True,
                    "message": "Trim Galore install/repair started. Progress will appear in the process console.",
                    "job_id": job_id,
                }
            )
            return

        if self.path == "/run-srna-mapping":
            try:
                settings = parse_srna_mapping_payload(payload)
                job_id = start_srna_mapping_job(settings)
            except ValueError as exc:
                self.respond_json({"ok": False, "message": str(exc)}, status=400)
                return
            except Exception as exc:
                self.respond_json({"ok": False, "message": f"Could not start sRNA mapping: {exc}"}, status=500)
                return

            self.respond_json(
                {
                    "ok": True,
                    "message": "sRNA mapping started. Progress will appear in the process console.",
                    "job_id": job_id,
                }
            )
            return

        if self.path == "/run-srna-dsrna-identification":
            try:
                settings = parse_srna_dsrna_payload(payload)
                job_id = start_srna_dsrna_identification_job(settings)
            except ValueError as exc:
                self.respond_json({"ok": False, "message": str(exc)}, status=400)
                return
            except Exception as exc:
                self.respond_json({"ok": False, "message": f"Could not start sRNA-based dsRNA identification: {exc}"}, status=500)
                return

            self.respond_json(
                {
                    "ok": True,
                    "message": "sRNA-based dsRNA identification started. Progress will appear in the process console.",
                    "job_id": job_id,
                }
            )
            return

        if self.path == "/run-srna-control-filtering":
            try:
                settings = parse_srna_control_payload(payload)
                job_id = start_srna_control_filtering_job(settings)
            except ValueError as exc:
                self.respond_json({"ok": False, "message": str(exc)}, status=400)
                return
            except Exception as exc:
                self.respond_json({"ok": False, "message": f"Could not start sRNA control filtering: {exc}"}, status=500)
                return

            self.respond_json(
                {
                    "ok": True,
                    "message": "sRNA control filtering started. Progress will appear in the process console.",
                    "job_id": job_id,
                }
            )
            return

        if self.path == "/run-degradome-analysis":
            try:
                settings = parse_degradome_payload(payload)
                job_id = start_degradome_analysis_job(settings)
            except ValueError as exc:
                self.respond_json({"ok": False, "message": str(exc)}, status=400)
                return
            except Exception as exc:
                self.respond_json({"ok": False, "message": f"Could not start degradome analysis: {exc}"}, status=500)
                return

            self.respond_json(
                {
                    "ok": True,
                    "message": "Degradome analysis started. Progress will appear in the process console.",
                    "job_id": job_id,
                }
            )
            return

        if self.path == "/run-target-prediction":
            try:
                settings = parse_target_prediction_payload(payload)
                job_id = start_target_prediction_job(settings)
            except ValueError as exc:
                self.respond_json({"ok": False, "message": str(exc)}, status=400)
                return
            except Exception as exc:
                self.respond_json({"ok": False, "message": f"Could not start target prediction: {exc}"}, status=500)
                return

            self.respond_json(
                {
                    "ok": True,
                    "message": "Target prediction started. Progress will appear in the process console.",
                    "job_id": job_id,
                }
            )
            return

        if self.path == "/run-module":
            module_key = str(payload.get("module_key", ""))
            if module_for_key(module_key) is None or module_key == "preprocessing":
                self.respond_json({"ok": False, "message": "Unknown pipeline section."}, status=400)
                return

            try:
                params = payload.get("params", {})
                if not isinstance(params, dict):
                    raise ValueError("Invalid module parameters.")
                prepare_tool_output_dir(module_key)
                if module_key == "dsrna-plotter":
                    job_id = start_dsrna_plotter_job()
                    self.respond_json(
                        {
                            "ok": True,
                            "message": "dsRNA Plotter started. Progress will appear in the process console.",
                            "job_id": job_id,
                        }
                    )
                    return
                if module_key not in {'dsrna-identification', 'fasta-deduplication'}:
                    raise ValueError('This tool is not available yet.')
                job_id = start_module_job(module_key, params)
                self.respond_json({'ok': True, 'job_id': job_id, 'message': 'Analysis started.'})
                return
            except ValueError as exc:
                self.respond_json({"ok": False, "message": str(exc)}, status=400)
                return
            except Exception as exc:
                self.respond_json({"ok": False, "message": f"Run failed: {exc}"}, status=500)
                return

        self.respond_text("Not found", status=404)

    def read_json(self) -> dict[str, Any]:
        length = int(self.headers.get("Content-Length", "0"))
        raw = self.rfile.read(length) if length else b"{}"
        try:
            payload = json.loads(raw.decode("utf-8"))
        except json.JSONDecodeError:
            return {}
        return payload if isinstance(payload, dict) else {}

    def respond_html(self, body: str, status: int = 200) -> None:
        self.respond_bytes(body.encode("utf-8"), "text/html; charset=utf-8", status)

    def respond_json(self, payload: dict[str, Any], status: int = 200) -> None:
        self.respond_bytes(json.dumps(payload).encode("utf-8"), "application/json; charset=utf-8", status)

    def respond_text(self, body: str, status: int = 200) -> None:
        self.respond_bytes(body.encode("utf-8"), "text/plain; charset=utf-8", status)

    def respond_bytes(self, body: bytes, content_type: str, status: int) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Cache-Control", "no-store, max-age=0")
        self.send_header("Pragma", "no-cache")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, format: str, *args: object) -> None:
        return


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    apply_global_plot_settings(load_state())
    port = find_free_port()
    server = ThreadingHTTPServer((HOST, port), InciRequestHandler)
    url = f"http://{HOST}:{port}/"
    print(f"INCI pipeline interface is running at {url}")
    print("Keep this PyCharm run process active while using the dashboard.")
    webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nINCI interface stopped.")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
