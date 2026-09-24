#!/usr/bin/env python3
"""Directional dsRNA identification from trimmed paired-end RNA-seq reads.

Reproducibility note for the dsRNA_seq_prelim analysis
=====================================================

Example command used to reproduce Plant_dsNode343 / chloroplast dsRNA hits
with 17_L1 as the reference sample and ShortCut RNase III sample 22_L1 as the
comparison sample:

    ./dsRNA_identification \\
      --pair 17_L1_2_val_2.fq.gz \\
      --pair 22_L1_2_val_2.fq.gz \\
      --library-pairs 17_L1=17440462 \\
      --library-pairs 22_L1=186157194 \\
      --reference-fasta /Volumes/Doga_2TB/dsRNA_seq_prelim/dsRNA_seq_trimmed/fasta/Plant_dsNode343_sequences.fasta \\
      --outdir outputs/dsRNA_identification_17L1_ref_22L1_comp_Plant_dsNode343_200bp \\
      --dsRNA-length 200 \\
      --comparison-sample 22_L1 \\
      --min-reduction 0.20 \\
      --threads 7 \\
      --min-mapq 10 \\
      --direction-source read1

The same run without requiring the 20% ShortCut RNase III reduction filter can
be reproduced by adding:

    --no-require-reduction

Example command used to reproduce the standalone Bnapus_dsNode343 plot and
whole-sequence dsRNA score. Use --mode single when providing this as a FASTA,
because otherwise FASTA input is scanned in fixed windows:

    ./dsRNA_identification \\
      --pair 17_L1 17_L1_2_val_2.fq.gz \\
      --pair 22_L1 22_L1_2_val_2.fq.gz \\
      --library-pairs 17_L1=17440462 \\
      --library-pairs 22_L1=186157194 \\
      --reference-fasta /Volumes/Doga_2TB/dsRNA_seq_prelim/dsRNA_seq_trimmed/fasta/Bnapus_dsNode343.fasta \\
      --mode single \\
      --outdir outputs/dsRNA_identification_17L1_22L1_Bnapus_dsNode343 \\
      --threads 7 \\
      --min-mapq 10 \\
      --direction-source read1

The script maps paired reads to a supplied sequence/reference, collapses each
pair into one fragment footprint, assigns fragment direction from read 1 by
default, and reports/plots all coverage and dsRNA scores in global CPM:

    coverage_cpm = per-base fragment coverage * 1,000,000 / trimmed_read_pairs

Input reads are assumed to be already trimmed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import heapq
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib-codex"))

import matplotlib.pyplot as plt
import numpy as np
from plot_style import apply_matplotlib_style
import pysam


@dataclass(frozen=True)
class SampleInput:
    name: str
    r1: Path
    r2: Path
    library_pairs: int | None = None


@dataclass
class SampleCoverage:
    sample: SampleInput
    library_pairs: int
    sense: dict[str, np.ndarray]
    antisense: dict[str, np.ndarray]
    stats: Counter


@dataclass
class WindowMetrics:
    source_contig: str
    start0: int
    end0: int
    length: int
    mean_sense_cpm: float
    mean_antisense_cpm: float
    mean_duplex_cpm: float
    overlap_balance_percent: float
    bidirectional_fraction_percent: float
    dsRNAness_percent: float
    expression_weighted_dsrna_score: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Identify and plot directional dsRNA candidates from trimmed paired-end reads.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    input_group = parser.add_argument_group("reads")
    input_group.add_argument(
        "--pair",
        nargs=3,
        action="append",
        metavar=("SAMPLE", "R1_FASTQ", "R2_FASTQ"),
        help="Sample name and trimmed paired FASTQ files. Repeat for multiple samples.",
    )
    input_group.add_argument(
        "--read-pair",
        nargs=2,
        action="append",
        metavar=("R1_FASTQ", "R2_FASTQ"),
        help="Trimmed paired FASTQ files; sample name is inferred from R1. Repeat for multiple samples.",
    )
    input_group.add_argument(
        "--library-pairs",
        action="append",
        default=[],
        metavar="SAMPLE=COUNT",
        help="Optional global CPM denominator override. By default records in R1 are counted.",
    )

    reference_group = parser.add_argument_group("reference")
    reference_group.add_argument("--reference-fasta", type=Path, help="FASTA/FNA reference to map against.")
    reference_group.add_argument("--sequence", help="Inline single reference sequence.")
    reference_group.add_argument("--sequence-name", default="query_sequence", help="Name for --sequence.")

    parser.add_argument("--outdir", type=Path, default=Path("dsRNA_identification"))
    parser.add_argument("--dsRNA-length", type=int, default=200, help="Scoring window length for FASTA/FNA scanning.")
    parser.add_argument("--window-step", type=int, default=1, help="Sliding-window step for FASTA/FNA scanning.")
    parser.add_argument("--top-n", type=int, default=100, help="Number of ranked candidates to plot.")
    parser.add_argument(
        "--candidate-buffer",
        type=int,
        default=5000,
        help="Number of high-scoring windows retained before nonredundant top-N selection.",
    )
    parser.add_argument(
        "--max-overlap-fraction",
        type=float,
        default=0.5,
        help="Maximum reciprocal overlap among ranked windows. Use 1 to allow redundant windows.",
    )
    parser.add_argument(
        "--plot-extension-windows",
        type=float,
        default=1.0,
        help="Plot extension upstream/downstream as multiples of dsRNA-length.",
    )
    parser.add_argument("--smooth-span", type=int, default=25, help="Moving-average span for plots.")
    parser.add_argument("--threads", type=int, default=7)
    parser.add_argument("--min-mapq", type=int, default=10)
    parser.add_argument(
        "--direction-source",
        choices=("read1", "read2"),
        default="read1",
        help="Mate whose alignment orientation defines the true fragment direction.",
    )
    parser.add_argument("--comparison-sample", help="Optional second sample used for reduction-adjusted scoring.")
    parser.add_argument(
        "--min-reduction",
        type=float,
        default=0.20,
        help="Minimum fractional reduction in either sense or antisense mean CPM for comparison-adjusted ranking.",
    )
    parser.add_argument(
        "--no-require-reduction",
        action="store_true",
        help="When --comparison-sample is set, keep ranking by reference score even if reduction is below threshold.",
    )
    parser.add_argument(
        "--mode",
        choices=("auto", "single", "scan"),
        default="auto",
        help="single plots/scores whole sequence; scan ranks fixed windows.",
    )
    parser.add_argument("--keep-reference-copy", action="store_true", help="Copy/write the reference FASTA into outdir.")
    return parser.parse_args()


def sanitize_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") or "sample"


def infer_sample_name(r1: Path) -> str:
    name = r1.name
    for suffix in (".fq.gz", ".fastq.gz", ".fq", ".fastq"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    patterns = [
        r"(_R?1(_001)?$)",
        r"(_1_val_1$)",
        r"(_1$)",
    ]
    for pattern in patterns:
        name = re.sub(pattern, "", name)
    return sanitize_name(name)


def parse_library_overrides(items: list[str]) -> dict[str, int]:
    overrides = {}
    for item in items:
        if "=" not in item:
            raise ValueError(f"--library-pairs must be SAMPLE=COUNT, got {item!r}")
        sample, value = item.split("=", 1)
        overrides[sample] = int(value.replace(",", ""))
    return overrides


def load_samples(args: argparse.Namespace) -> list[SampleInput]:
    overrides = parse_library_overrides(args.library_pairs)
    samples: list[SampleInput] = []
    for item in args.pair or []:
        name, r1, r2 = item
        samples.append(SampleInput(sanitize_name(name), Path(r1), Path(r2), overrides.get(sanitize_name(name))))
    for item in args.read_pair or []:
        r1, r2 = map(Path, item)
        name = infer_sample_name(r1)
        samples.append(SampleInput(name, r1, r2, overrides.get(name)))
    if not samples:
        raise ValueError("Provide at least one --pair SAMPLE R1 R2 or --read-pair R1 R2.")
    seen = set()
    for sample in samples:
        if sample.name in seen:
            raise ValueError(f"Duplicate sample name: {sample.name}")
        seen.add(sample.name)
        if not sample.r1.exists() or not sample.r2.exists():
            raise FileNotFoundError(f"Missing FASTQ for sample {sample.name}: {sample.r1}, {sample.r2}")
    return samples


def open_maybe_gzip(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def count_fastq_records(path: Path) -> int:
    lines = 0
    with open_maybe_gzip(path) as handle:
        for lines, _line in enumerate(handle, start=1):
            pass
    if lines % 4:
        raise ValueError(f"FASTQ line count is not divisible by 4: {path}")
    return lines // 4


def write_inline_reference(sequence: str, name: str, outdir: Path) -> Path:
    reference = outdir / "reference" / f"{sanitize_name(name)}.fasta"
    reference.parent.mkdir(parents=True, exist_ok=True)
    seq = re.sub(r"\s+", "", sequence).upper()
    with reference.open("w") as handle:
        handle.write(f">{sanitize_name(name)}\n")
        for idx in range(0, len(seq), 80):
            handle.write(seq[idx : idx + 80] + "\n")
    return reference


def load_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    name: str | None = None
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                sequences[name] = []
            elif name is None:
                raise ValueError(f"FASTA sequence encountered before header in {path}")
            else:
                sequences[name].append(line.upper())
    return {name: "".join(chunks) for name, chunks in sequences.items()}


def copy_reference_if_needed(reference: Path, outdir: Path, requested: bool) -> None:
    if not requested:
        return
    dest = outdir / "reference" / reference.name
    dest.parent.mkdir(parents=True, exist_ok=True)
    if reference.resolve() != dest.resolve():
        shutil.copy2(reference, dest)


def merge_blocks(blocks: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not blocks:
        return []
    blocks = sorted(blocks)
    merged = [blocks[0]]
    for start, end in blocks[1:]:
        old_start, old_end = merged[-1]
        if start <= old_end:
            merged[-1] = (old_start, max(old_end, end))
        else:
            merged.append((start, end))
    return merged


def parse_pair(
    reads: list[pysam.AlignedSegment],
    min_mapq: int,
    direction_source: str,
) -> tuple[str, str, list[tuple[int, int]]] | None:
    usable = [
        read
        for read in reads
        if not read.is_unmapped
        and not read.mate_is_unmapped
        and not read.is_secondary
        and not read.is_supplementary
        and read.mapping_quality >= min_mapq
    ]
    read1s = [read for read in usable if read.is_read1]
    read2s = [read for read in usable if read.is_read2]
    if not read1s or not read2s:
        return None
    direction_reads = read1s if direction_source == "read1" else read2s
    direction_read = max(direction_reads, key=lambda read: read.mapping_quality)
    reference_name = direction_read.reference_name
    if reference_name is None or any(read.reference_name != reference_name for read in usable):
        return None
    blocks = merge_blocks([block for read in usable for block in read.get_blocks()])
    if not blocks:
        return None
    direction = "antisense" if direction_read.is_reverse else "sense"
    return reference_name, direction, blocks


def increment_diff(diff: np.ndarray, start: int, end: int) -> None:
    if start < end:
        diff[start] += 1
        diff[end] -= 1


def map_sample(
    sample: SampleInput,
    reference: Path,
    lengths: dict[str, int],
    outdir: Path,
    threads: int,
    min_mapq: int,
    direction_source: str,
) -> SampleCoverage:
    library_pairs = sample.library_pairs if sample.library_pairs is not None else count_fastq_records(sample.r1)
    logs = outdir / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    log_path = logs / f"{sample.name}.minimap2.log"
    cmd = [
        "minimap2",
        "-ax",
        "sr",
        "-t",
        str(threads),
        "--secondary=no",
        "--sam-hit-only",
        str(reference),
        str(sample.r1),
        str(sample.r2),
    ]
    stats: Counter = Counter()
    diffs: dict[str, dict[str, np.ndarray]] = {}
    pending: dict[str, list[pysam.AlignedSegment]] = defaultdict(list)

    with log_path.open("w") as log_handle:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log_handle)
        assert process.stdout is not None
        with pysam.AlignmentFile(process.stdout, "r") as sam:
            for read in sam:
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                bucket = pending[read.query_name]
                bucket.append(read)
                if not (any(item.is_read1 for item in bucket) and any(item.is_read2 for item in bucket)):
                    continue
                parsed = parse_pair(bucket, min_mapq=min_mapq, direction_source=direction_source)
                del pending[read.query_name]
                if parsed is None:
                    stats["discarded_pairs"] += 1
                    continue
                reference_name, direction, blocks = parsed
                if reference_name not in lengths:
                    continue
                if reference_name not in diffs:
                    diffs[reference_name] = {
                        "sense": np.zeros(lengths[reference_name] + 1, dtype=np.int32),
                        "antisense": np.zeros(lengths[reference_name] + 1, dtype=np.int32),
                    }
                target = diffs[reference_name][direction]
                for start, end in blocks:
                    increment_diff(target, max(0, start), min(end, lengths[reference_name]))
                stats["collapsed_pairs"] += 1
                stats[f"{direction}_pairs"] += 1
        return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(f"minimap2 failed for {sample.name}; see {log_path}")

    sense: dict[str, np.ndarray] = {}
    antisense: dict[str, np.ndarray] = {}
    for contig, arrays in diffs.items():
        sense[contig] = np.cumsum(arrays["sense"][:-1])
        antisense[contig] = np.cumsum(arrays["antisense"][:-1])
    stats["library_pairs"] = library_pairs
    return SampleCoverage(sample, library_pairs, sense, antisense, stats)


def coverage_array(coverage: SampleCoverage, contig: str, direction: str, length: int) -> np.ndarray:
    source = coverage.sense if direction == "sense" else coverage.antisense
    if contig not in source:
        return np.zeros(length, dtype=np.float64)
    return source[contig].astype(np.float64)


def moving_average(values: np.ndarray, span: int) -> np.ndarray:
    if span <= 1:
        return values
    kernel = np.ones(span, dtype=np.float64) / span
    return np.convolve(values, kernel, mode="same")


def score_interval(
    coverage: SampleCoverage,
    contig: str,
    start0: int,
    end0: int,
    lengths: dict[str, int],
) -> WindowMetrics:
    length = end0 - start0
    scale = 1_000_000.0 / coverage.library_pairs if coverage.library_pairs else 0.0
    sense = coverage_array(coverage, contig, "sense", lengths[contig])[start0:end0] * scale
    antisense = coverage_array(coverage, contig, "antisense", lengths[contig])[start0:end0] * scale
    total = sense + antisense
    duplex = np.minimum(sense, antisense)
    total_sum = float(total.sum())
    duplex_sum = float(duplex.sum())
    overlap_balance = 100.0 * 2.0 * duplex_sum / total_sum if total_sum else 0.0
    bidirectional_fraction = 100.0 * np.count_nonzero((sense > 0) & (antisense > 0)) / length if length else 0.0
    dsrna = overlap_balance * math.sqrt(bidirectional_fraction / 100.0) if bidirectional_fraction else 0.0
    mean_duplex = float(np.mean(duplex)) if length else 0.0
    return WindowMetrics(
        source_contig=contig,
        start0=start0,
        end0=end0,
        length=length,
        mean_sense_cpm=float(np.mean(sense)) if length else 0.0,
        mean_antisense_cpm=float(np.mean(antisense)) if length else 0.0,
        mean_duplex_cpm=mean_duplex,
        overlap_balance_percent=overlap_balance,
        bidirectional_fraction_percent=bidirectional_fraction,
        dsRNAness_percent=dsrna,
        expression_weighted_dsrna_score=mean_duplex * dsrna / 100.0,
    )


def metric_row(prefix: str, metrics: WindowMetrics) -> dict[str, object]:
    return {
        f"{prefix}_mean_sense_cpm": metrics.mean_sense_cpm,
        f"{prefix}_mean_antisense_cpm": metrics.mean_antisense_cpm,
        f"{prefix}_mean_duplex_cpm": metrics.mean_duplex_cpm,
        f"{prefix}_overlap_balance_percent": metrics.overlap_balance_percent,
        f"{prefix}_bidirectional_fraction_percent": metrics.bidirectional_fraction_percent,
        f"{prefix}_dsRNAness_percent": metrics.dsRNAness_percent,
        f"{prefix}_expression_weighted_dsrna_score": metrics.expression_weighted_dsrna_score,
    }


def reciprocal_overlap(left: dict[str, object], right: dict[str, object]) -> float:
    if left["source_contig"] != right["source_contig"]:
        return 0.0
    overlap = max(0, min(int(left["end_1based"]), int(right["end_1based"])) - max(int(left["start_1based"]), int(right["start_1based"])) + 1)
    if overlap <= 0:
        return 0.0
    return max(overlap / int(left["length"]), overlap / int(right["length"]))


def select_nonredundant(rows: list[dict[str, object]], top_n: int, max_overlap: float) -> list[dict[str, object]]:
    selected: list[dict[str, object]] = []
    for row in rows:
        if all(reciprocal_overlap(row, old) <= max_overlap for old in selected):
            selected.append(row)
            if len(selected) >= top_n:
                break
    return selected


def iter_reference_windows(lengths: dict[str, int], window_length: int, step: int) -> Iterable[tuple[str, int, int]]:
    for contig, length in lengths.items():
        if length <= 0:
            continue
        if length <= window_length:
            yield contig, 0, length
            continue
        for start0 in range(0, length - window_length + 1, step):
            yield contig, start0, start0 + window_length


def scan_windows(
    reference_coverage: SampleCoverage,
    comparison_coverage: SampleCoverage | None,
    lengths: dict[str, int],
    window_length: int,
    step: int,
    top_n: int,
    candidate_buffer: int,
    max_overlap: float,
    min_reduction: float,
    require_reduction: bool,
) -> list[dict[str, object]]:
    heap: list[tuple[float, int, dict[str, object]]] = []
    counter = 0
    heap_limit = max(top_n, candidate_buffer)
    for contig, start0, end0 in iter_reference_windows(lengths, window_length, step):
        ref = score_interval(reference_coverage, contig, start0, end0, lengths)
        if ref.mean_duplex_cpm <= 0 or ref.dsRNAness_percent <= 0:
            continue
        row: dict[str, object] = {
            "source_contig": contig,
            "start_1based": start0 + 1,
            "end_1based": end0,
            "length": end0 - start0,
            **metric_row("reference", ref),
        }
        rank_score = ref.expression_weighted_dsrna_score
        if comparison_coverage is not None:
            comp = score_interval(comparison_coverage, contig, start0, end0, lengths)
            sense_reduction = 1.0 - comp.mean_sense_cpm / ref.mean_sense_cpm if ref.mean_sense_cpm > 0 else 0.0
            antisense_reduction = 1.0 - comp.mean_antisense_cpm / ref.mean_antisense_cpm if ref.mean_antisense_cpm > 0 else 0.0
            max_reduction = max(sense_reduction, antisense_reduction)
            passes = max_reduction >= min_reduction
            row.update(metric_row("comparison", comp))
            row.update(
                {
                    "sense_reduction_fraction": sense_reduction,
                    "antisense_reduction_fraction": antisense_reduction,
                    "max_reduction_fraction": max_reduction,
                    "passes_reduction_threshold": passes,
                }
            )
            if require_reduction and not passes:
                rank_score = 0.0
            else:
                rank_score *= 1.0 + max(0.0, max_reduction)
        row["rank_score"] = rank_score
        if rank_score <= 0:
            continue
        counter += 1
        item = (rank_score, counter, row)
        if len(heap) < heap_limit:
            heapq.heappush(heap, item)
        elif rank_score > heap[0][0]:
            heapq.heapreplace(heap, item)

    buffered = [item[2] for item in sorted(heap, key=lambda item: item[0], reverse=True)]
    selected = select_nonredundant(buffered, top_n=top_n, max_overlap=max_overlap)
    for rank, row in enumerate(selected, start=1):
        row["rank"] = rank
    return selected


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    if fieldnames is None:
        fieldnames = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_fasta_record(handle, name: str, sequence: str, width: int = 80) -> None:
    handle.write(f">{name}\n")
    for idx in range(0, len(sequence), width):
        handle.write(sequence[idx : idx + width] + "\n")


def write_ranked_fastas(
    ranked: list[dict[str, object]],
    sequences: dict[str, str],
    outdir: Path,
    extension_windows: float,
) -> None:
    fasta_dir = outdir / "fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    scoring_path = fasta_dir / "top_hits.scoring_windows.fasta"
    context_path = fasta_dir / "top_hits.plot_contexts.fasta"
    with scoring_path.open("w") as scoring_handle, context_path.open("w") as context_handle:
        for row in ranked:
            contig = str(row["source_contig"])
            start0 = int(row["start_1based"]) - 1
            end0 = int(row["end_1based"])
            context_start0, context_end0 = context_bounds(start0, end0, len(sequences[contig]), extension_windows)
            rank = int(row["rank"])
            score = float(row["rank_score"])
            dsrna = float(row["reference_dsRNAness_percent"])
            header = (
                f"rank_{rank:03d}|{contig}:{start0 + 1}-{end0}|"
                f"score={score:.6g}|dsRNAness={dsrna:.3f}"
            )
            context_header = (
                f"rank_{rank:03d}|{contig}:{context_start0 + 1}-{context_end0}|"
                f"scored={start0 + 1}-{end0}|score={score:.6g}|dsRNAness={dsrna:.3f}"
            )
            write_fasta_record(scoring_handle, header, sequences[contig][start0:end0])
            write_fasta_record(context_handle, context_header, sequences[contig][context_start0:context_end0])


def plot_interval(
    coverages: list[SampleCoverage],
    lengths: dict[str, int],
    contig: str,
    start0: int,
    end0: int,
    score_start0: int,
    score_end0: int,
    out_png: Path,
    title: str,
    smooth_span: int,
) -> None:
    style = apply_matplotlib_style()
    x = np.arange(start0 + 1, end0 + 1)
    fig, ax = plt.subplots(figsize=(style["figure_width"], style["figure_height"]), constrained_layout=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    local_abs_max = 0.0
    for idx, coverage in enumerate(coverages):
        scale = 1_000_000.0 / coverage.library_pairs if coverage.library_pairs else 0.0
        sense = coverage_array(coverage, contig, "sense", lengths[contig])[start0:end0] * scale
        antisense = coverage_array(coverage, contig, "antisense", lengths[contig])[start0:end0] * scale
        sense = moving_average(sense, smooth_span)
        antisense = moving_average(antisense, smooth_span)
        local_abs_max = max(local_abs_max, float(np.max(np.abs(sense))) if len(sense) else 0.0, float(np.max(np.abs(antisense))) if len(antisense) else 0.0)
        color = colors[idx % len(colors)]
        ax.plot(x, sense, color=color, linewidth=style["line_width"], label=f"{coverage.sample.name} sense")
        ax.plot(x, -antisense, color=color, linewidth=style["line_width"], linestyle="--", label=f"{coverage.sample.name} antisense")
    ax.axhline(0, color="0.25", linewidth=max(0.4, float(style["line_width"]) * 0.4))
    ax.axvspan(score_start0 + 1, score_end0, color="#4C78A8", alpha=0.14)
    if local_abs_max:
        ax.set_ylim(-local_abs_max * 1.25, local_abs_max * 1.25)
    ax.set_title(title, fontsize=style["title_size"])
    ax.set_xlabel(f"Position on {contig}")
    ax.set_ylabel("Global CPM: + sense / - antisense")
    ax.grid(axis="x", color="0.88", linewidth=style["grid_width"])
    ax.legend(ncols=min(4, len(coverages) * 2), frameon=False)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=style["dpi"])
    plt.close(fig)


def context_bounds(start0: int, end0: int, contig_len: int, extension_windows: float) -> tuple[int, int]:
    length = end0 - start0
    extension = int(round(length * extension_windows))
    return max(0, start0 - extension), min(contig_len, end0 + extension)


def plot_ranked_windows(
    coverages: list[SampleCoverage],
    lengths: dict[str, int],
    ranked: list[dict[str, object]],
    outdir: Path,
    extension_windows: float,
    smooth_span: int,
) -> None:
    plot_dir = outdir / "plots" / "top_windows"
    for row in ranked:
        contig = str(row["source_contig"])
        start0 = int(row["start_1based"]) - 1
        end0 = int(row["end_1based"])
        context_start, context_end = context_bounds(start0, end0, lengths[contig], extension_windows)
        rank = int(row["rank"])
        filename = f"rank_{rank:03d}_{sanitize_name(contig)}_{start0 + 1}_{end0}.png"
        reduction_label = ""
        if "sense_reduction_fraction" in row and "antisense_reduction_fraction" in row:
            reduction_label = (
                f" | reduction S {100.0 * float(row['sense_reduction_fraction']):.1f}%, "
                f"AS {100.0 * float(row['antisense_reduction_fraction']):.1f}%"
            )
        title = (
            f"rank {rank}: {contig}:{start0 + 1}-{end0} | "
            f"score {float(row['rank_score']):,.2f} | dsRNAness {float(row['reference_dsRNAness_percent']):.1f}%"
            f"{reduction_label}"
        )
        plot_interval(
            coverages,
            lengths,
            contig,
            context_start,
            context_end,
            start0,
            end0,
            plot_dir / filename,
            title,
            smooth_span,
        )


def write_sample_stats(path: Path, coverages: list[SampleCoverage]) -> None:
    rows = []
    for coverage in coverages:
        for metric in ("library_pairs", "collapsed_pairs", "sense_pairs", "antisense_pairs", "discarded_pairs"):
            rows.append({"sample": coverage.sample.name, "metric": metric, "value": coverage.stats.get(metric, 0)})
    write_csv(path, rows, ["sample", "metric", "value"])


def single_sequence_analysis(
    coverages: list[SampleCoverage],
    lengths: dict[str, int],
    outdir: Path,
    smooth_span: int,
) -> None:
    contig = next(iter(lengths))
    rows = []
    for coverage in coverages:
        metrics = score_interval(coverage, contig, 0, lengths[contig], lengths)
        rows.append({"sample": coverage.sample.name, "source_contig": contig, "start_1based": 1, "end_1based": lengths[contig], **metric_row("whole", metrics)})
    write_csv(outdir / "tables" / "single_sequence_scores.csv", rows)
    plot_interval(
        coverages,
        lengths,
        contig,
        0,
        lengths[contig],
        0,
        lengths[contig],
        outdir / "plots" / f"{sanitize_name(contig)}.whole_sequence.global_cpm.png",
        f"{contig} directional coverage and whole-sequence dsRNA score",
        smooth_span,
    )


def ranked_window_analysis(
    coverages: list[SampleCoverage],
    sequences: dict[str, str],
    lengths: dict[str, int],
    args: argparse.Namespace,
) -> None:
    reference = coverages[0]
    comparison = None
    if args.comparison_sample:
        comparison = next((item for item in coverages if item.sample.name == args.comparison_sample), None)
        if comparison is None:
            raise ValueError(f"--comparison-sample {args.comparison_sample!r} was not among mapped samples.")
    ranked = scan_windows(
        reference,
        comparison,
        lengths,
        window_length=args.dsRNA_length,
        step=args.window_step,
        top_n=args.top_n,
        candidate_buffer=args.candidate_buffer,
        max_overlap=args.max_overlap_fraction,
        min_reduction=args.min_reduction,
        require_reduction=bool(comparison is not None and not args.no_require_reduction),
    )
    write_csv(args.outdir / "tables" / "ranked_windows.csv", ranked)
    write_ranked_fastas(ranked, sequences, args.outdir, args.plot_extension_windows)

    sample_rows = []
    for row in ranked:
        contig = str(row["source_contig"])
        start0 = int(row["start_1based"]) - 1
        end0 = int(row["end_1based"])
        for coverage in coverages:
            metrics = score_interval(coverage, contig, start0, end0, lengths)
            sample_rows.append(
                {
                    "rank": row["rank"],
                    "sample": coverage.sample.name,
                    "source_contig": contig,
                    "start_1based": start0 + 1,
                    "end_1based": end0,
                    **metric_row("sample", metrics),
                }
            )
    write_csv(args.outdir / "tables" / "ranked_windows_by_sample.csv", sample_rows)
    plot_ranked_windows(
        coverages,
        lengths,
        ranked,
        args.outdir,
        extension_windows=args.plot_extension_windows,
        smooth_span=args.smooth_span,
    )


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    if bool(args.reference_fasta) == bool(args.sequence):
        raise ValueError("Provide exactly one of --reference-fasta or --sequence.")
    reference = (
        write_inline_reference(args.sequence, args.sequence_name, args.outdir)
        if args.sequence
        else args.reference_fasta
    )
    assert reference is not None
    copy_reference_if_needed(reference, args.outdir, args.keep_reference_copy)
    sequences = load_fasta(reference)
    lengths = {name: len(seq) for name, seq in sequences.items()}
    if not lengths:
        raise ValueError(f"No sequences found in {reference}")
    if shutil.which("minimap2") is None:
        raise RuntimeError("minimap2 was not found on PATH.")

    samples = load_samples(args)
    coverages = []
    for sample in samples:
        print(f"[dsRNA_identification] mapping {sample.name}", file=sys.stderr)
        coverages.append(
            map_sample(
                sample,
                reference,
                lengths,
                args.outdir,
                threads=args.threads,
                min_mapq=args.min_mapq,
                direction_source=args.direction_source,
            )
        )
    write_sample_stats(args.outdir / "tables" / "sample_stats.csv", coverages)

    mode = args.mode
    if mode == "auto":
        mode = "single" if args.sequence else "scan"
    if mode == "single":
        single_sequence_analysis(coverages, lengths, args.outdir, args.smooth_span)
    else:
        ranked_window_analysis(coverages, sequences, lengths, args)
    print(f"[dsRNA_identification] done: {args.outdir}")


if __name__ == "__main__":
    main()
