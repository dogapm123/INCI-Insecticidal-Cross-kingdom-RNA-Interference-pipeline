#!/usr/bin/env python3
"""Map unpaired sRNA reads to template contigs and plot perfect-match hits.

The mapping step uses Bowtie 1 with a configurable mismatch allowance so the
SAM keeps near matches. By default, simple-sequence filtering is applied only
to template-mapped reads to avoid writing a large prefiltered FASTQ. Downstream
scoring and plots use only perfect matches and normalize CPM against every read
in the original FASTQ, before filtering.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import re
import shutil
import subprocess
import sys
import urllib.request
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO

os.environ.setdefault("MPLCONFIGDIR", str(Path(os.environ.get("TMPDIR", "/tmp")) / "matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from plot_style import apply_matplotlib_style


DNA_BASES = set("ACGT")
DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
RUNINFO_URL = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc={accession}"


@dataclass
class FastqRecord:
    name: str
    sequence: str
    plus: str
    quality: str


@dataclass
class Contig:
    name: str
    description: str
    sequence: str


@dataclass
class PerfectAlignment:
    qname: str
    contig: str
    start0: int
    end0: int
    is_reverse: bool


@dataclass
class SamRead:
    qname: str
    sequence: str
    quality: str
    has_mapped_alignment: bool = False
    has_perfect_alignment: bool = False


def open_text(path: Path, mode: str = "rt") -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return path.open(mode)


def iter_fastq(path: Path) -> Iterable[FastqRecord]:
    with open_text(path) as handle:
        while True:
            name = handle.readline()
            if not name:
                return
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not quality:
                raise ValueError(f"Truncated FASTQ record in {path}")
            name = name.rstrip("\n")
            sequence = sequence.rstrip("\n")
            plus = plus.rstrip("\n")
            quality = quality.rstrip("\n")
            if not name.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"Malformed FASTQ record near {name!r} in {path}")
            if len(sequence) != len(quality):
                raise ValueError(
                    f"FASTQ sequence/quality length mismatch for {name} in {path}: "
                    f"{len(sequence)} bases vs {len(quality)} qualities"
                )
            yield FastqRecord(name=name, sequence=sequence, plus=plus, quality=quality)


def read_fasta(path: Path) -> list[Contig]:
    contigs: list[Contig] = []
    name = ""
    description = ""
    chunks: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    contigs.append(Contig(name=name, description=description, sequence="".join(chunks).upper()))
                description = line[1:]
                name = description.split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name:
        contigs.append(Contig(name=name, description=description, sequence="".join(chunks).upper()))
    if not contigs:
        raise ValueError(f"No contigs found in {path}")
    return contigs


def shannon_entropy(sequence: str) -> float:
    counts = Counter(base for base in sequence if base in DNA_BASES)
    total = sum(counts.values())
    if total == 0:
        return 0.0
    return -sum((count / total) * math.log2(count / total) for count in counts.values())


def longest_homopolymer(sequence: str) -> int:
    longest = 0
    current_base = ""
    current = 0
    for base in sequence:
        if base == current_base:
            current += 1
        else:
            current_base = base
            current = 1
        longest = max(longest, current)
    return longest


def dominant_kmer_fraction(sequence: str, k: int) -> float:
    if len(sequence) < k:
        return 0.0
    kmers = [sequence[idx : idx + k] for idx in range(0, len(sequence) - k + 1)]
    kmers = [kmer for kmer in kmers if set(kmer) <= DNA_BASES]
    if not kmers:
        return 0.0
    return max(Counter(kmers).values()) / len(kmers)


def low_complexity_reason(
    sequence: str,
    min_len: int,
    max_len: int,
    max_n_fraction: float,
    max_base_fraction: float,
    min_entropy: float,
    max_run_fraction: float,
    max_run_bases: int,
    max_dinuc_fraction: float,
    max_trinuc_fraction: float,
) -> str | None:
    seq = sequence.upper().replace("U", "T")
    length = len(seq)
    if length < min_len:
        return "too_short"
    if max_len and length > max_len:
        return "too_long"
    if not set(seq) <= (DNA_BASES | {"N"}):
        return "non_acgtn"
    if length and seq.count("N") / length > max_n_fraction:
        return "too_many_n"

    acgt = [base for base in seq if base in DNA_BASES]
    if not acgt:
        return "no_acgt_bases"
    counts = Counter(acgt)
    if len(counts) <= 2:
        return "one_or_two_base_alphabet"
    if max(counts.values()) / len(acgt) >= max_base_fraction:
        return "dominant_single_base"
    if shannon_entropy(seq) < min_entropy:
        return "low_shannon_entropy"

    longest_run = longest_homopolymer(seq)
    if longest_run >= max_run_bases or longest_run / length >= max_run_fraction:
        return "long_homopolymer_run"
    if dominant_kmer_fraction(seq, 2) >= max_dinuc_fraction:
        return "dominant_dinucleotide"
    if dominant_kmer_fraction(seq, 3) >= max_trinuc_fraction:
        return "dominant_trinucleotide"
    return None


def filter_fastq(args: argparse.Namespace, filtered_fastq: Path) -> Counter[str]:
    stats: Counter[str] = Counter()
    filtered_fastq.parent.mkdir(parents=True, exist_ok=True)
    with filtered_fastq.open("wt") as out:
        for record in iter_fastq(args.fastq):
            stats["original_reads"] += 1
            reason = low_complexity_reason(
                record.sequence,
                min_len=args.min_len,
                max_len=args.max_len,
                max_n_fraction=args.max_n_fraction,
                max_base_fraction=args.max_base_fraction,
                min_entropy=args.min_entropy,
                max_run_fraction=args.max_run_fraction,
                max_run_bases=args.max_run_bases,
                max_dinuc_fraction=args.max_dinuc_fraction,
                max_trinuc_fraction=args.max_trinuc_fraction,
            )
            if reason:
                stats[f"filtered_{reason}"] += 1
                continue
            stats["reads_passing_complexity_filter"] += 1
            out.write(f"{record.name}\n{record.sequence.upper().replace('U', 'T')}\n{record.plus}\n{record.quality}\n")
    if stats["original_reads"] == 0:
        raise ValueError(f"No reads found in {args.fastq}")
    return stats


def count_fastq_records(path: Path) -> int:
    lines = 0
    with open_text(path) as handle:
        for lines, _ in enumerate(handle, start=1):
            pass
    if lines % 4:
        raise ValueError(f"FASTQ line count is not divisible by 4 in {path}")
    if lines == 0:
        raise ValueError(f"No reads found in {path}")
    return lines // 4


def write_stats(path: Path, stats: Counter[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"])
        writer.writeheader()
        for metric in sorted(stats):
            writer.writerow({"metric": metric, "value": stats[metric]})


def require_tool(name: str) -> str:
    exe = shutil.which(name)
    if not exe:
        raise SystemExit(f"Required executable not found on PATH: {name}")
    return exe


def build_bowtie_index(fasta: Path, index_prefix: Path, threads: int) -> None:
    require_tool("bowtie-build")
    expected = Path(f"{index_prefix}.1.ebwt")
    if expected.exists():
        return
    index_prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["bowtie-build", "--threads", str(threads), str(fasta), str(index_prefix)]
    subprocess.run(cmd, check=True)


def run_bowtie(
    index_prefix: Path,
    fastq: Path,
    sam_path: Path,
    mismatches: int,
    threads: int,
    report_all: bool = True,
) -> None:
    style = apply_matplotlib_style()
    aligner = shutil.which("bowtie-align-s") or require_tool("bowtie")
    sam_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        aligner,
        "-S",
        "--no-unal",
        "-v",
        str(mismatches),
        "-p",
        str(threads),
        "-x",
        str(index_prefix),
        str(fastq),
        str(sam_path),
    ]
    if report_all:
        cmd.insert(2, "-a")
    subprocess.run(cmd, check=True)


def mapped_query_names(sam_path: Path, perfect_only: bool = False) -> set[str]:
    names: set[str] = set()
    with sam_path.open() as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            flag = int(fields[1])
            if flag & 4:
                continue
            if perfect_only and mismatch_count(fields) != 0:
                continue
            names.add(fields[0])
    return names


def summarize_sam_mismatch_bins(
    sam_path: Path,
    max_mismatches: int,
    prefix: str,
    denominator: int | None = None,
) -> Counter[str]:
    stats: Counter[str] = Counter()

    def flush(qname: str | None, best_nm: int | None) -> None:
        if qname is None:
            return
        stats[f"{prefix}_reads_examined"] += 1
        if best_nm is None:
            stats[f"{prefix}_unmapped_reads"] += 1
            return
        stats[f"{prefix}_mapped_reads"] += 1
        if 0 <= best_nm <= max_mismatches:
            stats[f"{prefix}_best_nm_{best_nm}_reads"] += 1
        else:
            stats[f"{prefix}_best_nm_gt_{max_mismatches}_reads"] += 1

    current_qname: str | None = None
    current_best_nm: int | None = None
    with sam_path.open() as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            qname = fields[0]
            if current_qname is None:
                current_qname = qname
            if qname != current_qname:
                flush(current_qname, current_best_nm)
                current_qname = qname
                current_best_nm = None

            flag = int(fields[1])
            if flag & 4:
                continue
            nm = mismatch_count(fields)
            if nm is None:
                continue
            current_best_nm = nm if current_best_nm is None else min(current_best_nm, nm)
    flush(current_qname, current_best_nm)

    if denominator is None:
        denominator = int(stats[f"{prefix}_reads_examined"])
    stats[f"{prefix}_denominator_reads"] = denominator
    mapped = int(stats[f"{prefix}_mapped_reads"])
    unmapped = int(stats[f"{prefix}_unmapped_reads"])
    stats[f"{prefix}_mapped_percent"] = mapped * 100.0 / denominator if denominator else 0.0
    stats[f"{prefix}_unmapped_percent"] = unmapped * 100.0 / denominator if denominator else 0.0
    cumulative = 0
    for nm in range(max_mismatches + 1):
        count = int(stats[f"{prefix}_best_nm_{nm}_reads"])
        cumulative += count
        stats[f"{prefix}_best_nm_{nm}_percent"] = count * 100.0 / denominator if denominator else 0.0
        stats[f"{prefix}_best_nm_le_{nm}_percent"] = cumulative * 100.0 / denominator if denominator else 0.0
    return stats


def write_mismatch_summary(
    path: Path,
    sample: str,
    stats: Counter[str],
    sections: list[tuple[str, int]],
) -> None:
    rows: list[dict[str, str | int | float]] = []
    for prefix, max_mismatches in sections:
        denominator = int(stats.get(f"{prefix}_denominator_reads", stats.get(f"{prefix}_reads_examined", 0)))
        if denominator == 0 and f"{prefix}_reads_examined" not in stats:
            continue
        mapped = int(stats.get(f"{prefix}_mapped_reads", 0))
        unmapped = int(stats.get(f"{prefix}_unmapped_reads", 0))
        rows.append(
            {
                "sample": sample,
                "target": prefix,
                "category": "mapped_total",
                "reads": mapped,
                "percent_of_denominator": mapped * 100.0 / denominator if denominator else 0.0,
                "denominator_reads": denominator,
            }
        )
        rows.append(
            {
                "sample": sample,
                "target": prefix,
                "category": "unmapped",
                "reads": unmapped,
                "percent_of_denominator": unmapped * 100.0 / denominator if denominator else 0.0,
                "denominator_reads": denominator,
            }
        )
        cumulative = 0
        for nm in range(max_mismatches + 1):
            count = int(stats.get(f"{prefix}_best_nm_{nm}_reads", 0))
            cumulative += count
            rows.append(
                {
                    "sample": sample,
                    "target": prefix,
                    "category": f"best_nm_{nm}",
                    "reads": count,
                    "percent_of_denominator": count * 100.0 / denominator if denominator else 0.0,
                    "denominator_reads": denominator,
                }
            )
            rows.append(
                {
                    "sample": sample,
                    "target": prefix,
                    "category": f"best_nm_le_{nm}",
                    "reads": cumulative,
                    "percent_of_denominator": cumulative * 100.0 / denominator if denominator else 0.0,
                    "denominator_reads": denominator,
                }
            )

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample", "target", "category", "reads", "percent_of_denominator", "denominator_reads"],
        )
        writer.writeheader()
        writer.writerows(rows)


def write_template_mapped_fastq_from_sam(
    sam_path: Path,
    out_fastq: Path,
    args: argparse.Namespace,
    perfect_only: bool = False,
) -> int:
    written = 0
    out_fastq.parent.mkdir(parents=True, exist_ok=True)

    def flush(read: SamRead | None) -> None:
        nonlocal written
        if read is None or not read.has_mapped_alignment:
            return
        if perfect_only and not read.has_perfect_alignment:
            return
        if mapped_read_filter_reason(read, args):
            return
        out.write(f"@{read.qname}\n{read.sequence}\n+\n{read.quality}\n")
        written += 1

    with sam_path.open() as handle, out_fastq.open("wt") as out:
        current_read: SamRead | None = None
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            qname = fields[0]
            if current_read is None:
                sequence, quality = original_read_from_sam_fields(fields)
                current_read = SamRead(qname=qname, sequence=sequence, quality=quality)
            if qname != current_read.qname:
                flush(current_read)
                sequence, quality = original_read_from_sam_fields(fields)
                current_read = SamRead(qname=qname, sequence=sequence, quality=quality)
            flag = int(fields[1])
            if flag & 4:
                continue
            current_read.has_mapped_alignment = True
            if mismatch_count(fields) == 0:
                current_read.has_perfect_alignment = True
        flush(current_read)
    return written


def write_selected_fastq(source_fastq: Path, selected_names: set[str], out_fastq: Path) -> int:
    written = 0
    out_fastq.parent.mkdir(parents=True, exist_ok=True)
    with out_fastq.open("wt") as out:
        for record in iter_fastq(source_fastq):
            qname = record.name[1:].split()[0]
            if qname not in selected_names:
                continue
            out.write(f"{record.name}\n{record.sequence}\n{record.plus}\n{record.quality}\n")
            written += 1
    return written


def parse_tags(fields: list[str]) -> dict[str, str]:
    tags: dict[str, str] = {}
    for field in fields[11:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def cigar_reference_length(cigar: str, query_sequence: str) -> int:
    if cigar == "*":
        return len(query_sequence)
    total = 0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        if op in {"M", "D", "N", "=", "X"}:
            total += int(length)
    return total


def mismatch_count(fields: list[str]) -> int | None:
    tags = parse_tags(fields)
    for tag in ("NM", "XM", "XA"):
        if tag in tags:
            try:
                return int(tags[tag])
            except ValueError:
                return None
    md = tags.get("MD")
    if md is None:
        return None
    mismatches = 0
    idx = 0
    while idx < len(md):
        char = md[idx]
        if char.isdigit():
            idx += 1
        elif char == "^":
            idx += 1
            while idx < len(md) and md[idx].isalpha():
                idx += 1
        else:
            mismatches += 1
            idx += 1
    return mismatches


def original_read_from_sam_fields(fields: list[str]) -> tuple[str, str]:
    sequence = fields[9]
    quality = fields[10]
    flag = int(fields[1])
    if flag & 16:
        sequence = sequence.translate(DNA_COMPLEMENT)[::-1]
        quality = quality[::-1]
    return sequence.upper().replace("U", "T"), quality


def mapped_read_filter_reason(read: SamRead | None, args: argparse.Namespace) -> str | None:
    if read is None or not read.has_mapped_alignment:
        return None
    return low_complexity_reason(
        read.sequence,
        min_len=args.min_len,
        max_len=args.max_len,
        max_n_fraction=args.max_n_fraction,
        max_base_fraction=args.max_base_fraction,
        min_entropy=args.min_entropy,
        max_run_fraction=args.max_run_fraction,
        max_run_bases=args.max_run_bases,
        max_dinuc_fraction=args.max_dinuc_fraction,
        max_trinuc_fraction=args.max_trinuc_fraction,
    )


def flush_alignment_group(
    read: SamRead | None,
    alignments: list[PerfectAlignment],
    mode: str,
    coverage: dict[str, dict[str, np.ndarray]],
    read_counts: dict[str, Counter[str]],
    stats: Counter[str],
    args: argparse.Namespace,
) -> None:
    if read is None or not read.has_mapped_alignment:
        return
    stats["template_mapped_reads_before_complexity_filter"] += 1
    reason = mapped_read_filter_reason(read, args)
    if reason:
        stats[f"mapped_filtered_{reason}"] += 1
        return
    stats["template_mapped_reads_passing_complexity_filter"] += 1
    if not alignments:
        return
    stats["reads_with_perfect_alignment"] += 1
    if mode == "unique" and len(alignments) != 1:
        stats["perfect_multimapper_reads_suppressed"] += 1
        return
    weight = 1.0 / len(alignments) if mode == "fractional" else 1.0
    for aln in alignments:
        direction = "antisense" if aln.is_reverse else "sense"
        coverage[aln.contig][direction][aln.start0] += weight
        coverage[aln.contig][direction][aln.end0] -= weight
        read_counts[aln.contig][direction] += weight
        stats["perfect_alignments_counted"] += weight


def scan_perfect_alignments(
    sam_path: Path,
    contigs: list[Contig],
    multi_hit_mode: str,
    args: argparse.Namespace,
) -> tuple[dict[str, dict[str, np.ndarray]], dict[str, Counter[str]], Counter[str]]:
    contig_lengths = {contig.name: len(contig.sequence) for contig in contigs}
    coverage = {
        contig.name: {
            "sense": np.zeros(len(contig.sequence) + 1, dtype=np.float64),
            "antisense": np.zeros(len(contig.sequence) + 1, dtype=np.float64),
        }
        for contig in contigs
    }
    read_counts: dict[str, Counter[str]] = defaultdict(Counter)
    stats: Counter[str] = Counter()

    current_read: SamRead | None = None
    current_alignments: list[PerfectAlignment] = []
    with sam_path.open() as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            qname = fields[0]
            if current_read is None:
                sequence, quality = original_read_from_sam_fields(fields)
                current_read = SamRead(qname=qname, sequence=sequence, quality=quality)
            if qname != current_read.qname:
                flush_alignment_group(current_read, current_alignments, multi_hit_mode, coverage, read_counts, stats, args)
                sequence, quality = original_read_from_sam_fields(fields)
                current_read = SamRead(qname=qname, sequence=sequence, quality=quality)
                current_alignments = []

            flag = int(fields[1])
            if flag & 4:
                stats["sam_unmapped_records"] += 1
                continue
            current_read.has_mapped_alignment = True
            stats["sam_mapped_records_up_to_mismatch_limit"] += 1
            mismatches = mismatch_count(fields)
            if mismatches != 0:
                continue
            current_read.has_perfect_alignment = True
            contig = fields[2]
            if contig not in contig_lengths:
                continue
            start0 = int(fields[3]) - 1
            end0 = start0 + cigar_reference_length(fields[5], fields[9])
            end0 = min(end0, contig_lengths[contig])
            if start0 < 0 or end0 <= start0:
                continue
            current_alignments.append(
                PerfectAlignment(
                    qname=qname,
                    contig=contig,
                    start0=start0,
                    end0=end0,
                    is_reverse=bool(flag & 16),
                )
            )
    flush_alignment_group(current_read, current_alignments, multi_hit_mode, coverage, read_counts, stats, args)
    return coverage, read_counts, stats


def summarize_hits(
    contigs: list[Contig],
    coverage: dict[str, dict[str, np.ndarray]],
    read_counts: dict[str, Counter[str]],
    original_reads: int,
) -> list[dict[str, float | int | str]]:
    scale = 1_000_000.0 / original_reads
    rows: list[dict[str, float | int | str]] = []
    for contig in contigs:
        sense_depth = np.cumsum(coverage[contig.name]["sense"][:-1])
        antisense_depth = np.cumsum(coverage[contig.name]["antisense"][:-1])
        sense_reads = float(read_counts[contig.name]["sense"])
        antisense_reads = float(read_counts[contig.name]["antisense"])
        total_reads = sense_reads + antisense_reads
        if total_reads == 0:
            continue
        rows.append(
            {
                "contig": contig.name,
                "description": contig.description,
                "length": len(contig.sequence),
                "sense_perfect_reads_weighted": sense_reads,
                "antisense_perfect_reads_weighted": antisense_reads,
                "total_perfect_reads_weighted": total_reads,
                "sense_cpm": sense_reads * scale,
                "antisense_cpm": antisense_reads * scale,
                "total_cpm": total_reads * scale,
                "peak_sense_cpm": float(np.max(sense_depth) * scale),
                "peak_antisense_cpm": float(np.max(antisense_depth) * scale),
                "covered_bases": int(np.count_nonzero((sense_depth + antisense_depth) > 0)),
            }
        )
    rows.sort(key=lambda row: (float(row["total_cpm"]), float(row["peak_sense_cpm"]) + float(row["peak_antisense_cpm"])), reverse=True)
    return rows


def run_control_mapping(args: argparse.Namespace, sample: str, template_sam: Path, stats: Counter[str]) -> None:
    mapped_fastq = args.outdir / "control" / f"{sample}.template_mapped_reads.fq"
    written = write_template_mapped_fastq_from_sam(
        template_sam,
        mapped_fastq,
        args,
        perfect_only=args.control_from_perfect_template_hits,
    )
    stats["template_mapped_reads_for_control"] = written
    stats["template_mapped_reads_fastq_written"] = written
    if not written:
        stats["control_mapped_reads"] = 0
        stats["control_mapped_percent_of_template_mapped"] = 0
        write_stats(args.outdir / "tables" / f"{sample}.control_mapping_stats.csv", stats)
        return

    control_index = args.outdir / "control" / "bowtie_index" / args.control_fasta.stem
    build_bowtie_index(args.control_fasta, control_index, args.threads)
    control_sam = args.outdir / "control" / f"{sample}.template_mapped_vs_control.v{args.control_mismatches}.sam"
    if args.force or not control_sam.exists():
        run_bowtie(control_index, mapped_fastq, control_sam, args.control_mismatches, args.threads, report_all=False)
    control_names = mapped_query_names(control_sam, perfect_only=False)
    control_mismatch_stats = summarize_sam_mismatch_bins(
        control_sam,
        max_mismatches=args.control_mismatches,
        prefix="control",
        denominator=written,
    )
    stats.update(control_mismatch_stats)
    stats["control_mapped_reads"] = len(control_names)
    stats["control_mapped_percent_of_template_mapped"] = (
        len(control_names) * 100.0 / written if written else 0
    )
    write_stats(args.outdir / "tables" / f"{sample}.control_mapping_stats.csv", stats)


def write_summary(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    fields = [
        "contig",
        "description",
        "length",
        "sense_perfect_reads_weighted",
        "antisense_perfect_reads_weighted",
        "total_perfect_reads_weighted",
        "sense_cpm",
        "antisense_cpm",
        "total_cpm",
        "peak_sense_cpm",
        "peak_antisense_cpm",
        "covered_bases",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_top_coverage(
    path: Path,
    top_rows: list[dict[str, float | int | str]],
    coverage: dict[str, dict[str, np.ndarray]],
    original_reads: int,
) -> None:
    scale = 1_000_000.0 / original_reads
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["contig", "position_1based", "sense_cpm", "antisense_cpm", "signed_cpm"],
        )
        writer.writeheader()
        for row in top_rows:
            contig = str(row["contig"])
            sense = np.cumsum(coverage[contig]["sense"][:-1]) * scale
            antisense = np.cumsum(coverage[contig]["antisense"][:-1]) * scale
            for idx, (sense_cpm, antisense_cpm) in enumerate(zip(sense, antisense, strict=True), start=1):
                writer.writerow(
                    {
                        "contig": contig,
                        "position_1based": idx,
                        "sense_cpm": f"{sense_cpm:.6f}",
                        "antisense_cpm": f"{antisense_cpm:.6f}",
                        "signed_cpm": f"{sense_cpm - antisense_cpm:.6f}",
                    }
                )


def plot_hit(
    out: Path,
    row: dict[str, float | int | str],
    coverage: dict[str, dict[str, np.ndarray]],
    original_reads: int,
) -> None:
    scale = 1_000_000.0 / original_reads
    contig = str(row["contig"])
    sense = np.cumsum(coverage[contig]["sense"][:-1]) * scale
    antisense = np.cumsum(coverage[contig]["antisense"][:-1]) * scale
    x = np.arange(1, len(sense) + 1)
    fig, ax = plt.subplots(figsize=(style["figure_width"], style["figure_height"]), constrained_layout=True)
    draw_directional_axis(ax, x, sense, antisense, row)
    fig.savefig(out, dpi=style["dpi"])
    plt.close(fig)


def draw_directional_axis(ax, x, sense, antisense, row) -> None:
    style = apply_matplotlib_style()
    ax.fill_between(x, 0, sense, color="#1F77B4", alpha=0.35, step="mid", label="sense +")
    ax.plot(x, sense, color="#1F77B4", linewidth=style["line_width"])
    ax.fill_between(x, 0, -antisense, color="#D62728", alpha=0.35, step="mid", label="antisense -")
    ax.plot(x, -antisense, color="#D62728", linewidth=style["line_width"])
    ax.axhline(0, color="#303030", linewidth=max(0.4, float(style["line_width"]) * 0.4))
    ax.set_xlim(int(x[0]), int(x[-1]) if len(x) else 1)
    ax.set_ylabel("CPM\n+sense / -antisense")
    ax.set_title(
        f"{row['contig']} | total {float(row['total_cpm']):.2f} CPM "
        f"(+{float(row['sense_cpm']):.2f} / -{float(row['antisense_cpm']):.2f})",
        fontsize=style["title_size"],
    )
    ax.grid(axis="x", color="#E5E5E5", linewidth=style["grid_width"])
    ax.legend(frameon=False, loc="upper right")


def plot_top_hits(
    out: Path,
    top_rows: list[dict[str, float | int | str]],
    coverage: dict[str, dict[str, np.ndarray]],
    original_reads: int,
) -> None:
    if not top_rows:
        return
    style = apply_matplotlib_style()
    scale = 1_000_000.0 / original_reads
    fig_height = max(float(style["figure_height"]), 2.55 * len(top_rows))
    fig, axes = plt.subplots(len(top_rows), 1, figsize=(style["figure_width"], fig_height), constrained_layout=True)
    if len(top_rows) == 1:
        axes = [axes]
    for ax, row in zip(axes, top_rows, strict=True):
        contig = str(row["contig"])
        sense = np.cumsum(coverage[contig]["sense"][:-1]) * scale
        antisense = np.cumsum(coverage[contig]["antisense"][:-1]) * scale
        x = np.arange(1, len(sense) + 1)
        draw_directional_axis(ax, x, sense, antisense, row)
        ax.set_xlabel("position on template contig")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=style["dpi"])
    plt.close(fig)


def run_analysis(args: argparse.Namespace) -> None:
    args.outdir.mkdir(parents=True, exist_ok=True)
    sample = args.sample or args.fastq.name.removesuffix(".gz").removesuffix(".fq").removesuffix(".fastq")
    filter_stats: Counter[str]
    if args.filter_stage == "pre":
        filtered_fastq = args.outdir / "filtered_fastq" / f"{sample}.complexity_filtered.fq"
        filter_stats = filter_fastq(args, filtered_fastq)
        fastq_for_mapping = filtered_fastq
    else:
        filter_stats = Counter({"original_reads": count_fastq_records(args.fastq)})
        fastq_for_mapping = args.fastq
    write_stats(args.outdir / "tables" / f"{sample}.read_filter_stats.csv", filter_stats)

    contigs = read_fasta(args.fasta)
    index_prefix = args.outdir / "bowtie_index" / args.fasta.stem
    build_bowtie_index(args.fasta, index_prefix, args.threads)

    sam_path = args.outdir / "sam" / f"{sample}.v{args.mismatches}.sam"
    if args.force or not sam_path.exists():
        run_bowtie(index_prefix, fastq_for_mapping, sam_path, args.mismatches, args.threads)

    template_mismatch_stats = summarize_sam_mismatch_bins(
        sam_path,
        max_mismatches=args.mismatches,
        prefix="template",
        denominator=int(filter_stats["original_reads"]),
    )
    coverage, read_counts, alignment_stats = scan_perfect_alignments(sam_path, contigs, args.multi_hit_mode, args)
    all_stats = filter_stats + template_mismatch_stats + alignment_stats
    all_stats["template_contigs"] = len(contigs)
    all_stats["complexity_filter_stage"] = args.filter_stage
    mismatch_sections = [("template", args.mismatches)]
    if args.control_fasta:
        run_control_mapping(args, sample, sam_path, all_stats)
        mismatch_sections.append(("control", args.control_mismatches))
    write_mismatch_summary(
        args.outdir / "tables" / f"{sample}.mapping_mismatch_summary.csv",
        sample,
        all_stats,
        mismatch_sections,
    )
    write_stats(args.outdir / "tables" / f"{sample}.analysis_stats.csv", all_stats)

    summary_rows = summarize_hits(contigs, coverage, read_counts, int(filter_stats["original_reads"]))
    top_rows = summary_rows[: args.top_n]
    write_summary(args.outdir / "tables" / f"{sample}.perfect_hit_summary.csv", summary_rows)
    write_summary(args.outdir / "tables" / f"{sample}.top{args.top_n}_perfect_hits.csv", top_rows)
    write_top_coverage(
        args.outdir / "tables" / f"{sample}.top{args.top_n}.perfect_coverage_cpm.csv",
        top_rows,
        coverage,
        int(filter_stats["original_reads"]),
    )

    plots_dir = args.outdir / "plots"
    plot_top_hits(plots_dir / f"{sample}.top{args.top_n}.perfect_hits.directional_cpm.png", top_rows, coverage, int(filter_stats["original_reads"]))
    individual_dir = plots_dir / f"{sample}.top{args.top_n}_individual"
    individual_dir.mkdir(parents=True, exist_ok=True)
    for rank, row in enumerate(top_rows, start=1):
        safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(row["contig"]))[:120]
        plot_hit(individual_dir / f"rank_{rank:03d}_{safe_name}.png", row, coverage, int(filter_stats["original_reads"]))

    print(f"Wrote SAM: {sam_path}")
    print(f"Wrote summaries: {args.outdir / 'tables'}")
    print(f"Wrote plots: {plots_dir}")


def fetch_runinfo(accession: str) -> list[dict[str, str]]:
    url = RUNINFO_URL.format(accession=accession)
    with urllib.request.urlopen(url) as response:
        text = response.read().decode("utf-8")
    reader = csv.DictReader(text.splitlines())
    return list(reader)


def download_sra(args: argparse.Namespace) -> None:
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = fetch_runinfo(args.bioproject)
    if args.organism:
        rows = [row for row in rows if row.get("ScientificName", "").lower() == args.organism.lower()]
    if args.max_runs:
        rows = rows[: args.max_runs]
    if not rows:
        raise SystemExit(f"No SRA runs found for {args.bioproject}")

    metadata_path = args.outdir / f"{args.bioproject}.runinfo.csv"
    with metadata_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    run_list_path = args.outdir / f"{args.bioproject}.runs.txt"
    run_list_path.write_text("\n".join(row["Run"] for row in rows) + "\n")
    print(f"Wrote {len(rows)} run records to {metadata_path}")

    if args.metadata_only:
        return
    require_tool("prefetch")
    require_tool("fasterq-dump")
    sra_dir = args.outdir / "sra"
    fastq_dir = args.outdir / "fastq"
    sra_dir.mkdir(exist_ok=True)
    fastq_dir.mkdir(exist_ok=True)
    for row in rows:
        run = row["Run"]
        print(f"Downloading {run}", flush=True)
        subprocess.run(["prefetch", "--output-directory", str(sra_dir), run], check=True)
        subprocess.run(
            ["fasterq-dump", "--threads", str(args.threads), "--outdir", str(fastq_dir), str(sra_dir / run)],
            check=True,
        )
        if args.gzip:
            for fastq in fastq_dir.glob(f"{run}*.fastq"):
                subprocess.run(["gzip", "-f", str(fastq)], check=True)


def add_run_parser(subparsers) -> None:
    parser = subparsers.add_parser("run", help="Map one unpaired FASTQ and plot top perfect template hits")
    parser.add_argument("--fastq", type=Path, required=True, help="Single-end sRNA FASTQ (.fq/.fastq, optionally .gz)")
    parser.add_argument("--fasta", type=Path, required=True, help="Template contigs FASTA")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--sample", default="", help="Sample name for output files; defaults to FASTQ basename")
    parser.add_argument("--threads", type=int, default=max(1, min(8, os.cpu_count() or 1)))
    parser.add_argument("--mismatches", type=int, default=3, help="Bowtie -v mismatch allowance for SAM generation")
    parser.add_argument("--top-n", type=int, default=100, help="Number of top template hits to table and plot")
    parser.add_argument("--force", action="store_true", help="Rerun Bowtie even if the SAM already exists")
    parser.add_argument(
        "--filter-stage",
        choices=["mapped", "pre"],
        default="mapped",
        help="Apply simple-sequence filtering only to template-mapped reads by default; use pre for legacy whole-FASTQ filtering before Bowtie",
    )
    parser.add_argument("--control-fasta", type=Path, help="Optional control FASTA to remap template-mapped reads against")
    parser.add_argument("--control-mismatches", type=int, default=3, help="Bowtie -v mismatch allowance for the control remap")
    parser.add_argument(
        "--control-from-perfect-template-hits",
        action="store_true",
        help="Use only perfect template-mapping reads for the control remap; default uses all reads mapped to templates within --mismatches",
    )
    parser.add_argument(
        "--multi-hit-mode",
        choices=["fractional", "all", "unique"],
        default="fractional",
        help="How perfect multi-mapping reads contribute to CPM coverage",
    )
    parser.add_argument("--min-len", type=int, default=18)
    parser.add_argument("--max-len", type=int, default=30, help="Use 0 to disable maximum length filtering")
    parser.add_argument("--max-n-fraction", type=float, default=0.10)
    parser.add_argument("--max-base-fraction", type=float, default=0.80)
    parser.add_argument("--min-entropy", type=float, default=1.25)
    parser.add_argument("--max-run-fraction", type=float, default=0.60)
    parser.add_argument("--max-run-bases", type=int, default=8)
    parser.add_argument("--max-dinuc-fraction", type=float, default=0.70)
    parser.add_argument("--max-trinuc-fraction", type=float, default=0.70)
    parser.set_defaults(func=run_analysis)


def add_download_parser(subparsers) -> None:
    parser = subparsers.add_parser("download-sra", help="Fetch SRA run metadata and FASTQ files for a BioProject")
    parser.add_argument("--bioproject", default="PRJNA1400957")
    parser.add_argument("--organism", default="Psylliodes chrysocephalus")
    parser.add_argument("--outdir", type=Path, default=Path("data") / "PRJNA1400957")
    parser.add_argument("--threads", type=int, default=max(1, min(8, os.cpu_count() or 1)))
    parser.add_argument("--metadata-only", action="store_true")
    parser.add_argument("--max-runs", type=int, default=0, help="Limit runs for a quick smoke download")
    parser.add_argument("--gzip", action="store_true", help="gzip FASTQ files after fasterq-dump")
    parser.set_defaults(func=download_sra)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    add_run_parser(subparsers)
    add_download_parser(subparsers)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    args.func(args)


if __name__ == "__main__":
    main()
