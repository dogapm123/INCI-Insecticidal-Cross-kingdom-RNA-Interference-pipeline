#!/usr/bin/env python3
"""Fast whole-genome 250 bp dsRNA candidate screen with paired-end follow-up plots.

Phase one maps R1 reads to the native genome with minimap2, assigns each PAF
alignment to fixed genomic bins, and ranks bins from the mean global-CPM
coverage across ssRNase replicates.  It deliberately keeps only primary PAF
hits for a rapid, conservative prescreen.  Phase two maps both mates to the
small FASTA of top-bin contexts, assigns each pair the strand of R1, and plots
mean directional coverage with SEM for ssRNase and dsRNase groups.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam


os_environ = __import__("os").environ
os_environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib-inci-genome-250bp"))


@dataclass(frozen=True)
class SampleInput:
    name: str
    r1: Path
    r2: Path | None = None


def parse_sample(value: str, paired: bool) -> SampleInput:
    pieces = value.split("=", 1)
    if len(pieces) != 2 or not pieces[0] or not pieces[1]:
        raise argparse.ArgumentTypeError("Sample must be NAME=R1 or NAME=R1,R2")
    name, raw_paths = pieces
    paths = [Path(path) for path in raw_paths.split(",")]
    if paired and len(paths) != 2:
        raise argparse.ArgumentTypeError("Paired sample must be NAME=R1,R2")
    if not paired and len(paths) != 1:
        raise argparse.ArgumentTypeError("Screen sample must be NAME=R1")
    return SampleInput(name=name, r1=paths[0], r2=paths[1] if paired else None)


def open_fastq(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def count_fastq_records(path: Path) -> int:
    with open_fastq(path) as handle:
        lines = sum(1 for _ in handle)
    if lines % 4:
        raise ValueError(f"FASTQ line count is not divisible by four: {path}")
    return lines // 4


def load_fasta_metadata(path: Path) -> tuple[dict[str, int], dict[str, str]]:
    lengths: dict[str, int] = {}
    descriptions: dict[str, str] = {}
    name: str | None = None
    for raw in path.open(errors="replace"):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            header = line[1:]
            name = header.split()[0]
            lengths[name] = 0
            descriptions[name] = header
        elif name is not None:
            lengths[name] += len(line)
    return lengths, descriptions


def ensure_mmi(genome_fasta: Path, index_path: Path, threads: int) -> Path:
    if index_path.exists() and index_path.stat().st_size > 0:
        return index_path
    index_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["minimap2", "-d", str(index_path), str(genome_fasta)]
    print("[screen] building minimap2 index", file=sys.stderr)
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    (index_path.parent / "minimap2_index.log").write_text(result.stderr)
    if result.returncode:
        raise RuntimeError(f"minimap2 index failed; see {index_path.parent / 'minimap2_index.log'}")
    return index_path


def fastq_stream_command(r1: Path, max_reads: int) -> tuple[list[str], bool]:
    if max_reads <= 0:
        return [str(r1)], False
    if shutil.which("seqkit") is None:
        raise RuntimeError("seqkit is required when --max-reads is positive")
    return ["seqkit", "head", "-n", str(max_reads), str(r1)], True


def add_interval(counts: dict[tuple[str, int], np.ndarray], target: str, start0: int, end0: int, bin_size: int, strand_index: int, weight: float = 1.0) -> None:
    """Add one PAF target span to its overlapping fixed bins."""
    if end0 <= start0:
        return
    first_bin = start0 // bin_size
    last_bin = (end0 - 1) // bin_size
    for bin_index in range(first_bin, last_bin + 1):
        bin_start = bin_index * bin_size
        overlap = min(end0, bin_start + bin_size) - max(start0, bin_start)
        if overlap <= 0:
            continue
        values = counts.setdefault((target, bin_index), np.zeros(2, dtype=np.float64))
        values[strand_index] += overlap * weight


def map_r1_to_bins(reference: Path, sample: SampleInput, bin_size: int, threads: int, min_mapq: int, max_reads: int, log_path: Path, include_secondary: bool = False, max_secondary: int = 1000) -> tuple[dict[tuple[str, int], np.ndarray], int, Counter[str]]:
    stream_cmd, streamed = fastq_stream_command(sample.r1, max_reads)
    minimap_cmd = ["minimap2", "-x", "sr", "-t", str(threads)]
    if include_secondary:
        minimap_cmd.extend(["-N", str(max_secondary)])
    else:
        minimap_cmd.append("--secondary=no")
    minimap_cmd.append(str(reference))
    if streamed:
        minimap_cmd.append("-")
    else:
        minimap_cmd.extend(stream_cmd)
    counts: dict[tuple[str, int], np.ndarray] = {}
    run_stats: Counter[str] = Counter()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log_handle:
        source = None
        if streamed:
            source = subprocess.Popen(stream_cmd, stdout=subprocess.PIPE, stderr=log_handle)
            assert source.stdout is not None
            mapper = subprocess.Popen(minimap_cmd, stdin=source.stdout, stdout=subprocess.PIPE, stderr=log_handle, text=True)
            source.stdout.close()
        else:
            mapper = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE, stderr=log_handle, text=True)
        assert mapper.stdout is not None
        current_query: str | None = None
        query_alignments: list[list[str]] = []

        def flush_query() -> None:
            if not query_alignments:
                return
            valid = [fields for fields in query_alignments if int(fields[11]) >= min_mapq]
            run_stats["low_mapq"] += len(query_alignments) - len(valid)
            if not valid:
                return
            weight = 1.0 / len(valid) if include_secondary else 1.0
            for fields in valid:
                target = fields[5]
                start0 = int(fields[7])
                end0 = int(fields[8])
                strand_index = 1 if fields[4] == "-" else 0
                add_interval(counts, target, start0, end0, bin_size, strand_index, weight)
            run_stats["mapped_reads"] += 1
            run_stats["kept_alignments"] += len(valid)
            run_stats["max_valid_mappings_per_read"] = max(run_stats["max_valid_mappings_per_read"], len(valid))
            if len(valid) > 1:
                run_stats["multimapped_reads"] += 1
            if include_secondary and len(valid) >= max_secondary + 1:
                run_stats["reads_at_secondary_cap"] += 1

        for line in mapper.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            query = fields[0]
            if current_query is not None and query != current_query:
                flush_query()
                query_alignments = []
            current_query = query
            query_alignments.append(fields)
        flush_query()
        mapper_rc = mapper.wait()
        source_rc = source.wait() if source is not None else 0
    if mapper_rc or source_rc:
        raise RuntimeError(f"minimap2 failed for {sample.name}; see {log_path}")
    denominator = max_reads if max_reads > 0 else count_fastq_records(sample.r1)
    return counts, denominator, run_stats


def merge_screen_counts(per_sample: list[tuple[dict[tuple[str, int], np.ndarray], int]], contig_lengths: dict[str, int], bin_size: int, descriptions: dict[str, str]) -> pd.DataFrame:
    all_keys: set[tuple[str, int]] = set()
    for counts, _denominator in per_sample:
        all_keys.update(counts)
    rows: list[dict[str, object]] = []
    for target, bin_index in all_keys:
        if target not in contig_lengths:
            continue
        start0 = bin_index * bin_size
        end0 = min(start0 + bin_size, contig_lengths[target])
        length = end0 - start0
        if length <= 0:
            continue
        cpm_areas = []
        for counts, denominator in per_sample:
            raw = counts.get((target, bin_index), np.zeros(2, dtype=np.float64))
            cpm_areas.append(raw * (1_000_000.0 / denominator))
        mean_areas = np.mean(np.vstack(cpm_areas), axis=0)
        sense_area, antisense_area = mean_areas
        duplex_area = min(sense_area, antisense_area)
        total_area = sense_area + antisense_area
        balance = 100.0 * (2.0 * duplex_area / total_area) if total_area else 0.0
        sense_depth = sense_area / length
        antisense_depth = antisense_area / length
        duplex_depth = duplex_area / length
        coverage_score = 100.0 * duplex_depth / (duplex_depth + 100.0) if duplex_depth else 0.0
        candidate_score = balance * coverage_score / 100.0
        target_bin = f"{target}_{start0 + 1}_{end0}"
        rows.append(
            {
                "target_bin": target_bin,
                "source_contig": target,
                "description": descriptions.get(target, target),
                "start_1based": start0 + 1,
                "end_1based": end0,
                "bin_length_nt": length,
                "sense_area_cpm_bp": sense_area,
                "antisense_area_cpm_bp": antisense_area,
                "sense_depth_cpm": sense_depth,
                "antisense_depth_cpm": antisense_depth,
                "duplex_depth_cpm": duplex_depth,
                "strand_balance_percent": balance,
                "coverage_score_percent": coverage_score,
                "candidate_score_percent": candidate_score,
            }
        )
    frame = pd.DataFrame(rows)
    if frame.empty:
        raise ValueError("No phase-one alignments passed the selected mapping threshold.")
    return frame.sort_values(["candidate_score_percent", "duplex_depth_cpm"], ascending=False).reset_index(drop=True)


def node_label_mask(frame: pd.DataFrame) -> pd.Series:
    return frame["source_contig"].eq("Bnapus_dsNode343")


def compact_label(row: pd.Series) -> str:
    if row["source_contig"] == "Bnapus_dsNode343":
        return f"dsNode343 {int(row['start_1based'])}-{int(row['end_1based'])}"
    return f"{row['source_contig']}:{int(row['start_1based'])}-{int(row['end_1based'])}"


def selected_rows_with_node(frame: pd.DataFrame, n: int) -> pd.DataFrame:
    selected = frame.head(n).copy()
    node = frame.loc[node_label_mask(frame)].head(1)
    if not node.empty and node.iloc[0]["target_bin"] not in set(selected["target_bin"]):
        selected = pd.concat([selected, node], ignore_index=True)
    return selected


def plot_summaries(frame: pd.DataFrame, outdir: Path, top_n: int) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    shown = selected_rows_with_node(frame, top_n).iloc[::-1]
    colors = np.where(shown["source_contig"].eq("Bnapus_dsNode343"), "#a51414", "#3e75a8")
    fig, ax = plt.subplots(figsize=(13, max(7, 0.31 * len(shown))), constrained_layout=True)
    ax.barh(np.arange(len(shown)), shown["candidate_score_percent"], color=colors)
    ax.set_yticks(np.arange(len(shown)))
    ax.set_yticklabels([compact_label(row) for _, row in shown.iterrows()], fontsize=7)
    ax.set_xlabel("total dsRNA score (%)")
    ax.set_title("Ranked dsRNA loci")
    ax.grid(axis="x", alpha=0.25)
    fig.savefig(outdir / "genome_250bp_ssRNase_candidate_score_barplot.png", dpi=220)
    plt.close(fig)

    positives = frame.loc[frame["candidate_score_percent"] > 0].copy()
    if len(positives) > 120_000:
        positives = positives.sample(120_000, random_state=17)
    fig, ax = plt.subplots(figsize=(10.5, 8), constrained_layout=True)
    sizes = np.clip(np.sqrt(positives["candidate_score_percent"].to_numpy()) * 3.0, 3, 65)
    ax.scatter(positives["strand_balance_percent"], positives["duplex_depth_cpm"], s=sizes, alpha=0.22, color="#618bb5", linewidths=0)
    top = selected_rows_with_node(frame, 20)
    for _, row in top.iterrows():
        color = "#a51414" if row["source_contig"] == "Bnapus_dsNode343" else "#173e67"
        ax.scatter(row["strand_balance_percent"], row["duplex_depth_cpm"], s=55, color=color, edgecolor="white", linewidth=0.5, zorder=3)
        ax.annotate(compact_label(row), (row["strand_balance_percent"], row["duplex_depth_cpm"]), xytext=(4, 3), textcoords="offset points", fontsize=6.5, color=color)
    ax.set_xlabel("sense/antisense balance (%)")
    ax.set_ylabel("mean duplex depth (global CPM)")
    ax.set_title("Strand balance and bidirectional depth")
    ax.grid(alpha=0.25)
    fig.savefig(outdir / "genome_250bp_ssRNase_balance_vs_duplex_depth.png", dpi=220)
    plt.close(fig)

    heat = selected_rows_with_node(frame, min(30, top_n))
    metric_keys = ["candidate_score_percent", "strand_balance_percent", "coverage_score_percent", "duplex_depth_cpm", "sense_depth_cpm", "antisense_depth_cpm"]
    metric_labels = ["candidate\nscore %", "balance\n%", "coverage\nscore %", "duplex\nCPM", "sense\nCPM", "antisense\nCPM"]
    values = heat[metric_keys].to_numpy(dtype=float)
    scaled = values.copy()
    for column in range(scaled.shape[1]):
        low, high = scaled[:, column].min(), scaled[:, column].max()
        scaled[:, column] = (scaled[:, column] - low) / (high - low) if high > low else 0.0
    fig, ax = plt.subplots(figsize=(10.5, max(7, 0.31 * len(heat))), constrained_layout=True)
    image = ax.imshow(scaled, aspect="auto", cmap="magma")
    ax.set_yticks(np.arange(len(heat)))
    ax.set_yticklabels([compact_label(row) for _, row in heat.iterrows()], fontsize=7)
    ax.set_xticks(np.arange(len(metric_labels)))
    ax.set_xticklabels(metric_labels)
    for row_index in range(len(heat)):
        for col_index, key in enumerate(metric_keys):
            value = values[row_index, col_index]
            color = "white" if scaled[row_index, col_index] < 0.55 else "black"
            ax.text(col_index, row_index, f"{value:.1f}", ha="center", va="center", fontsize=6.5, color=color)
    ax.set_title("Top dsRNA locus score components")
    fig.colorbar(image, ax=ax, label="within-column scaled value")
    fig.savefig(outdir / "genome_250bp_ssRNase_metric_heatmap.png", dpi=220)
    plt.close(fig)


def load_requested_contexts(genome_fasta: Path, node_fasta: Path, selected: pd.DataFrame, context_nt: int) -> dict[str, tuple[str, str, pd.Series]]:
    by_contig: dict[str, list[pd.Series]] = defaultdict(list)
    for _, row in selected.iterrows():
        by_contig[str(row["source_contig"])].append(row)
    records: dict[str, tuple[str, str, pd.Series]] = {}
    current_name: str | None = None
    chunks: list[str] = []

    def finish_record(name: str | None, sequence_chunks: list[str]) -> None:
        if name is None or name not in by_contig:
            return
        sequence = "".join(sequence_chunks).upper()
        for row in by_contig[name]:
            start0 = max(0, int(row["start_1based"]) - 1 - context_nt)
            end0 = min(len(sequence), int(row["end_1based"]) + context_nt)
            target = str(row["target_bin"])
            records[target] = (f"source={name} start={start0 + 1} end={end0} focus={int(row['start_1based'])}-{int(row['end_1based'])}", sequence[start0:end0], row)

    for raw in genome_fasta.open(errors="replace"):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish_record(current_name, chunks)
            current_name = line[1:].split()[0]
            chunks = []
        elif current_name in by_contig:
            chunks.append(line)
    finish_record(current_name, chunks)

    node_lengths, _ = load_fasta_metadata(node_fasta)
    if "Bnapus_dsNode343" in by_contig:
        if "Bnapus_dsNode343" not in node_lengths:
            raise ValueError("The supplied Node343 FASTA does not contain Bnapus_dsNode343")
        seq_chunks: list[str] = []
        found = False
        for raw in node_fasta.open(errors="replace"):
            line = raw.strip()
            if line.startswith(">"):
                found = line[1:].split()[0] == "Bnapus_dsNode343"
            elif found:
                seq_chunks.append(line)
        finish_record("Bnapus_dsNode343", seq_chunks)
    missing = [row["target_bin"] for _, row in selected.iterrows() if row["target_bin"] not in records]
    if missing:
        raise ValueError(f"Could not extract contexts for: {', '.join(map(str, missing[:5]))}")
    return records


def write_context_fasta(records: dict[str, tuple[str, str, pd.Series]], out_fasta: Path) -> None:
    with out_fasta.open("w") as handle:
        for target, (description, sequence, _row) in records.items():
            handle.write(f">{target} {description}\n")
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index:index + 80] + "\n")


def merge_blocks(blocks: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted(blocks):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def parse_paired_reads(reads: list[pysam.AlignedSegment], min_mapq: int) -> tuple[str, bool, list[tuple[int, int]]] | None:
    usable = [read for read in reads if not read.is_unmapped and not read.is_secondary and not read.is_supplementary and read.mapping_quality >= min_mapq]
    read1s = [read for read in usable if read.is_read1]
    read2s = [read for read in usable if read.is_read2]
    if not read1s or not read2s:
        return None
    read1 = max(read1s, key=lambda read: read.mapping_quality)
    target = read1.reference_name
    if target is None or any(read.reference_name != target for read in usable):
        return None
    blocks = merge_blocks([block for read in usable for block in read.get_blocks()])
    return target, read1.is_reverse, blocks if blocks else []


def map_paired_contexts(reference: Path, sample: SampleInput, lengths: dict[str, int], outdir: Path, threads: int, min_mapq: int) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    if sample.r2 is None:
        raise ValueError(f"Missing R2 for {sample.name}")
    sense = {name: np.zeros(length, dtype=float) for name, length in lengths.items()}
    antisense = {name: np.zeros(length, dtype=float) for name, length in lengths.items()}
    pending: dict[str, list[pysam.AlignedSegment]] = defaultdict(list)
    log_path = outdir / "logs" / f"{sample.name}.context.minimap2.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["minimap2", "-ax", "sr", "-t", str(threads), "--secondary=no", "--sam-hit-only", str(reference), str(sample.r1), str(sample.r2)]
    with log_path.open("w") as log_handle:
        mapper = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log_handle)
        assert mapper.stdout is not None
        with pysam.AlignmentFile(mapper.stdout, "r") as sam:
            for read in sam:
                bucket = pending[read.query_name]
                bucket.append(read)
                if not (any(item.is_read1 for item in bucket) and any(item.is_read2 for item in bucket)):
                    continue
                parsed = parse_paired_reads(bucket, min_mapq)
                del pending[read.query_name]
                if parsed is None:
                    continue
                target, reverse, blocks = parsed
                values = antisense[target] if reverse else sense[target]
                for start, end in blocks:
                    values[start:end] += 1.0
        rc = mapper.wait()
    if rc:
        raise RuntimeError(f"Context mapping failed for {sample.name}; see {log_path}")
    denominator = count_fastq_records(sample.r1)
    scale = 1_000_000.0 / denominator
    coverage: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    coverage_dir = outdir / "coverage"
    coverage_dir.mkdir(parents=True, exist_ok=True)
    for target in lengths:
        coverage[target] = (sense[target] * scale, antisense[target] * scale)
        np.savez_compressed(coverage_dir / f"{sample.name}.{target}.npz", sense_cpm=coverage[target][0], antisense_cpm=coverage[target][1], denominator_pairs=np.array([denominator]))
    return coverage


def smooth(values: np.ndarray, span: int) -> np.ndarray:
    if span <= 1:
        return values
    return np.convolve(values, np.ones(span) / span, mode="same")


def mean_sem(rows: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    matrix = np.vstack(rows)
    return matrix.mean(axis=0), matrix.std(axis=0, ddof=1) / np.sqrt(len(rows))


def plot_context_target(target: str, description: str, row: pd.Series, context_nt: int, group_a: list[dict[str, tuple[np.ndarray, np.ndarray]]], group_b: list[dict[str, tuple[np.ndarray, np.ndarray]]], outdir: Path, smooth_span: int) -> Path:
    fig, ax = plt.subplots(figsize=(12, 5.7), constrained_layout=True)
    x = np.arange(len(group_a[0][target][0]))
    settings = [(group_a, "#006400", "ssRNase: RNase If (dsRNA enrichment)"), (group_b, "#8B006B", "dsRNase: ShortCut RNase III (dsRNA depletion)")]
    for group, color, label in settings:
        sense_mean, sense_sem = mean_sem([smooth(item[target][0], smooth_span) for item in group])
        anti_mean, anti_sem = mean_sem([smooth(item[target][1], smooth_span) for item in group])
        ax.plot(x, sense_mean, color=color, linewidth=2.0, label=label)
        ax.plot(x, -anti_mean, color=color, linewidth=2.0, linestyle="--")
        ax.plot(x, sense_mean + sense_sem, color=color, linewidth=0.9, linestyle=":", alpha=0.75)
        ax.plot(x, sense_mean - sense_sem, color=color, linewidth=0.9, linestyle=":", alpha=0.75)
        ax.plot(x, -(anti_mean + anti_sem), color=color, linewidth=0.9, linestyle=":", alpha=0.75)
        ax.plot(x, -(anti_mean - anti_sem), color=color, linewidth=0.9, linestyle=":", alpha=0.75)
    focus_start = context_nt
    focus_end = context_nt + int(row["bin_length_nt"])
    ax.axvspan(focus_start, focus_end, color="#888888", alpha=0.10, zorder=0)
    ax.axvline(focus_start, color="#333333", linewidth=1.1, linestyle="--")
    ax.axvline(focus_end, color="#333333", linewidth=1.1, linestyle="--")
    ax.axhline(0, color="#202020", linewidth=1.2, zorder=0)
    ax.set_xlim(0, len(x) - 1)
    ax.set_xlabel("position in plotted context (nt)")
    ax.set_ylabel("normalized depth (CPM); + sense / - antisense")
    ax.set_title(f"{compact_label(row)} | score {row['candidate_score_percent']:.1f}% | balance {row['strand_balance_percent']:.1f}%")
    ax.legend(loc="lower left", frameon=False)
    ax.grid(axis="x", alpha=0.16)
    outpath = outdir / "plots" / f"{target}.png"
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=220)
    plt.close(fig)
    return outpath


def contact_sheet(paths: list[Path], output: Path, columns: int = 2) -> None:
    from PIL import Image, ImageOps

    images = [Image.open(path).convert("RGB") for path in paths]
    width = max(image.width for image in images)
    height = max(image.height for image in images)
    rows = int(np.ceil(len(images) / columns))
    sheet = Image.new("RGB", (columns * width, rows * height), "white")
    for index, image in enumerate(images):
        canvas = ImageOps.contain(image, (width, height))
        x = (index % columns) * width + (width - canvas.width) // 2
        y = (index // columns) * height + (height - canvas.height) // 2
        sheet.paste(canvas, (x, y))
    sheet.save(output)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--genome-fasta", type=Path, required=True)
    parser.add_argument("--genome-index", type=Path, help="Existing minimap2 genome .mmi index to reuse.")
    parser.add_argument("--node-fasta", type=Path, help="Optional separately mapped comparator FASTA. Omit with --skip-standalone-node for a fair native-reference screen.")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--screen-sample", action="append", required=True, type=lambda value: parse_sample(value, paired=False), help="Repeat NAME=R1 for ssRNase R1 files")
    parser.add_argument("--plot-group-a", action="append", type=lambda value: parse_sample(value, paired=True), help="Repeat NAME=R1,R2 for ssRNase paired files")
    parser.add_argument("--plot-group-b", action="append", type=lambda value: parse_sample(value, paired=True), help="Repeat NAME=R1,R2 for dsRNase paired files")
    parser.add_argument("--bin-size", type=int, default=250)
    parser.add_argument("--max-reads", type=int, default=1_000_000)
    parser.add_argument("--threads", type=int, default=7)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--summary-top-n", type=int, default=50)
    parser.add_argument("--context-nt", type=int, default=250)
    parser.add_argument("--smooth-span", type=int, default=25)
    parser.add_argument("--include-secondary", action="store_true", help="Retain minimap2 secondary R1 placements and split each read's contribution equally among its retained mappings.")
    parser.add_argument("--max-secondary", type=int, default=1000, help="Maximum secondary alignments emitted per read when --include-secondary is used.")
    parser.add_argument("--skip-standalone-node", action="store_true", help="Do not separately map Node343; use this for a fair native-genome-only ranking.")
    parser.add_argument("--phase1-only", action="store_true", help="Write phase-one ranking and summaries, then stop before paired top-hit plotting.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.bin_size != 250:
        print(f"[screen] using requested bin size {args.bin_size} bp", file=sys.stderr)
    args.outdir.mkdir(parents=True, exist_ok=True)
    if not args.skip_standalone_node and args.node_fasta is None:
        raise ValueError("--node-fasta is required unless --skip-standalone-node is used.")
    if not args.phase1_only and (not args.plot_group_a or not args.plot_group_b):
        raise ValueError("--plot-group-a and --plot-group-b are required unless --phase1-only is used.")
    reference_dir = args.outdir / "reference"
    index_path = args.genome_index or ensure_mmi(args.genome_fasta, reference_dir / f"{args.genome_fasta.stem}.mmi", args.threads)
    if not index_path.exists() or not index_path.stat().st_size:
        raise FileNotFoundError(f"Genome index does not exist or is empty: {index_path}")
    genome_lengths, genome_descriptions = load_fasta_metadata(args.genome_fasta)
    node_lengths: dict[str, int] = {}
    node_descriptions: dict[str, str] = {}
    if not args.skip_standalone_node:
        assert args.node_fasta is not None
        node_lengths, node_descriptions = load_fasta_metadata(args.node_fasta)
        if "Bnapus_dsNode343" not in node_lengths:
            raise ValueError("The standalone comparator FASTA must contain a Bnapus_dsNode343 record")

    per_sample: list[tuple[dict[tuple[str, int], np.ndarray], int]] = []
    log_rows: list[dict[str, object]] = []
    for sample in args.screen_sample:
        print(f"[screen] mapping R1: {sample.name}", file=sys.stderr)
        genome_counts, denominator, stats = map_r1_to_bins(index_path, sample, args.bin_size, args.threads, args.min_mapq, args.max_reads, args.outdir / "logs" / f"{sample.name}.genome.minimap2.log", args.include_secondary, args.max_secondary)
        node_counts: dict[tuple[str, int], np.ndarray] = {}
        node_stats: Counter[str] = Counter()
        if not args.skip_standalone_node:
            node_counts, node_denominator, node_stats = map_r1_to_bins(args.node_fasta, sample, args.bin_size, args.threads, args.min_mapq, args.max_reads, args.outdir / "logs" / f"{sample.name}.node343.minimap2.log", args.include_secondary, args.max_secondary)
            if denominator != node_denominator:
                raise RuntimeError("Genome and Node343 denominators differ")
            genome_counts.update(node_counts)
        per_sample.append((genome_counts, denominator))
        log_rows.append({"sample": sample.name, "denominator_reads": denominator, "genome_bins_with_hits": len(genome_counts) - len(node_counts), "node_bins_with_hits": len(node_counts), **{f"genome_{key}": value for key, value in stats.items()}, **{f"node_{key}": value for key, value in node_stats.items()}})

    all_lengths = {**genome_lengths, **({} if args.skip_standalone_node else node_lengths)}
    all_descriptions = {**genome_descriptions, **({} if args.skip_standalone_node else node_descriptions)}
    ranking = merge_screen_counts(per_sample, all_lengths, args.bin_size, all_descriptions)
    ranking.insert(0, "rank", np.arange(1, len(ranking) + 1))
    ranking_path = args.outdir / "whole_genome_250bp_ssRNase_ranked_candidates.csv"
    ranking.to_csv(ranking_path, index=False)
    pd.DataFrame(log_rows).to_csv(args.outdir / "phase1_mapping_summary.csv", index=False)
    plot_summaries(ranking, args.outdir / "summary_plots", args.summary_top_n)

    if args.phase1_only:
        print(f"[screen] ranking: {ranking_path}", file=sys.stderr)
        return

    top = ranking.head(args.top_n).copy()
    contexts = load_requested_contexts(args.genome_fasta, args.node_fasta, top, args.context_nt)
    context_fasta = args.outdir / "top20_250bp_contexts.fasta"
    write_context_fasta(contexts, context_fasta)
    top.to_csv(args.outdir / "top20_250bp_candidates.csv", index=False)
    context_lengths, _ = load_fasta_metadata(context_fasta)

    group_a = [map_paired_contexts(context_fasta, sample, context_lengths, args.outdir, args.threads, args.min_mapq) for sample in args.plot_group_a]
    group_b = [map_paired_contexts(context_fasta, sample, context_lengths, args.outdir, args.threads, args.min_mapq) for sample in args.plot_group_b]
    paths = [plot_context_target(str(row["target_bin"]), contexts[str(row["target_bin"])][0], row, args.context_nt, group_a, group_b, args.outdir, args.smooth_span) for _, row in top.iterrows()]
    contact_sheet(paths, args.outdir / "top20_250bp_directional_coverage_contact_sheet.png")
    print(f"[screen] ranking: {ranking_path}", file=sys.stderr)
    print(f"[screen] plots: {args.outdir / 'plots'}", file=sys.stderr)


if __name__ == "__main__":
    main()
