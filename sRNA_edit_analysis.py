#!/usr/bin/env python3
"""Analyze sRNA mismatches for substitutions and candidate 3-prime tailing.

The script reads Bowtie1 SAM files produced by sRNA_identification.py and a
template FASTA. It focuses on complexity-passing reads whose best template
alignment has 1..N mismatches, counts substitution spectra, and separately
flags candidate 3-prime tailing events where a 21 nt template-matching anchor
is followed by 1-3 non-template bases at the read's original 3-prime end.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", str(Path(os.environ.get("TMPDIR", "/tmp")) / "matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import sRNA_identification as srna


BASES = "ACGT"
SUBSTITUTIONS = [f"{ref}>{alt}" for ref in BASES for alt in BASES if alt != ref]
SUBSTITUTION_LABELS_RNA = [sub.replace("T", "U") for sub in SUBSTITUTIONS]
RC_MUTATION_FAMILIES = {
    "A>G": "A>G / U>C",
    "T>C": "A>G / U>C",
    "C>T": "C>U / G>A",
    "G>A": "C>U / G>A",
    "A>C": "A>C / U>G",
    "T>G": "A>C / U>G",
    "A>T": "A>U / U>A",
    "T>A": "A>U / U>A",
    "C>G": "C>G / G>C",
    "G>C": "C>G / G>C",
    "G>T": "G>U / C>A",
    "C>A": "G>U / C>A",
}
RC_MUTATION_FAMILY_ORDER = [
    "A>G / U>C",
    "C>U / G>A",
    "A>C / U>G",
    "A>U / U>A",
    "C>G / G>C",
    "G>U / C>A",
]


def substitution_label_rna(substitution: str) -> str:
    return substitution.replace("T", "U")


def complement_base(base: str) -> str:
    return base.translate(srna.DNA_COMPLEMENT).upper()


def oriented_substitution(aln: "SamAlignment", ref_base: str, read_base: str, orientation: str) -> str:
    if orientation == "molecule" and aln.is_reverse:
        return f"{complement_base(ref_base)}>{complement_base(read_base)}"
    return f"{ref_base}>{read_base}"


@dataclass
class SamAlignment:
    qname: str
    flag: int
    contig: str
    start0: int
    cigar: str
    sequence: str
    quality: str
    mismatches: list[tuple[int, str, str, int]]

    @property
    def is_reverse(self) -> bool:
        return bool(self.flag & 16)

    @property
    def end0(self) -> int:
        return self.start0 + len(self.sequence)

    @property
    def nm(self) -> int:
        return len(self.mismatches)


@dataclass
class ReadGroup:
    qname: str
    original_sequence: str
    alignments: list[SamAlignment] = field(default_factory=list)


def simple_m_cigar_length(cigar: str) -> int | None:
    match = re.fullmatch(r"(\d+)M", cigar)
    if not match:
        return None
    return int(match.group(1))


def alignment_from_fields(fields: list[str], contigs: dict[str, srna.Contig]) -> SamAlignment | None:
    flag = int(fields[1])
    if flag & 4:
        return None
    contig = fields[2]
    if contig not in contigs:
        return None
    sequence = fields[9].upper().replace("U", "T")
    cigar_len = simple_m_cigar_length(fields[5])
    if cigar_len is None or cigar_len != len(sequence):
        return None
    start0 = int(fields[3]) - 1
    if start0 < 0 or start0 + len(sequence) > len(contigs[contig].sequence):
        return None
    template = contigs[contig].sequence[start0 : start0 + len(sequence)]
    mismatches = [
        (start0 + offset, ref_base, read_base, offset)
        for offset, (ref_base, read_base) in enumerate(zip(template, sequence, strict=True))
        if ref_base in BASES and read_base in BASES and ref_base != read_base
    ]
    return SamAlignment(
        qname=fields[0],
        flag=flag,
        contig=contig,
        start0=start0,
        cigar=fields[5],
        sequence=sequence,
        quality=fields[10],
        mismatches=mismatches,
    )


def iter_read_groups(sam_path: Path, contigs: dict[str, srna.Contig]) -> Iterable[ReadGroup]:
    current: ReadGroup | None = None
    with sam_path.open() as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            qname = fields[0]
            if current is None:
                original_sequence, _ = srna.original_read_from_sam_fields(fields)
                current = ReadGroup(qname=qname, original_sequence=original_sequence)
            if qname != current.qname:
                yield current
                original_sequence, _ = srna.original_read_from_sam_fields(fields)
                current = ReadGroup(qname=qname, original_sequence=original_sequence)
            alignment = alignment_from_fields(fields, contigs)
            if alignment is not None:
                current.alignments.append(alignment)
    if current is not None:
        yield current


def low_complexity_reason(sequence: str, args: argparse.Namespace) -> str | None:
    return srna.low_complexity_reason(
        sequence,
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


def candidate_tail(aln: SamAlignment, anchor_len: int, max_tail_len: int) -> dict[str, int | str] | None:
    read_len = len(aln.sequence)
    tail_len = read_len - anchor_len
    if tail_len < 1 or tail_len > max_tail_len:
        return None

    mismatch_offsets = {offset for _, _, _, offset in aln.mismatches}
    if aln.is_reverse:
        anchor_offsets = range(tail_len, read_len)
        tail_offsets = range(0, tail_len)
        tail_sequence = aln.sequence[:tail_len].translate(srna.DNA_COMPLEMENT)[::-1].upper()
    else:
        anchor_offsets = range(0, anchor_len)
        tail_offsets = range(anchor_len, read_len)
        tail_sequence = aln.sequence[anchor_len:].upper()

    if any(offset in mismatch_offsets for offset in anchor_offsets):
        return None
    if not all(offset in mismatch_offsets for offset in tail_offsets):
        return None
    if not set(tail_sequence) <= set(BASES):
        return None

    return {
        "tail_length": tail_len,
        "tail_sequence": tail_sequence,
        "tail_is_homopolymer": int(len(set(tail_sequence)) == 1),
        "tail_base": tail_sequence[0] if len(set(tail_sequence)) == 1 else "mixed",
    }


def candidate_one_nt_3prime_tail_from_short_core(aln: SamAlignment, core_len: int) -> dict[str, int | str] | None:
    """Flag a 21 nt alignment that is a perfect 20 nt core plus one 3-prime tail-like base."""
    read_len = len(aln.sequence)
    if read_len != core_len + 1 or aln.nm != 1:
        return None

    mismatch_offsets = {offset for _, _, _, offset in aln.mismatches}
    if aln.is_reverse:
        if mismatch_offsets != {0}:
            return None
        tail_sequence = aln.sequence[0].translate(srna.DNA_COMPLEMENT).upper()
    else:
        if mismatch_offsets != {read_len - 1}:
            return None
        tail_sequence = aln.sequence[-1].upper()

    if not set(tail_sequence) <= set(BASES):
        return None
    return {
        "tail_length": 1,
        "tail_sequence": tail_sequence,
        "tail_is_homopolymer": 1,
        "tail_base": tail_sequence,
    }


def candidate_five_prime_addition(aln: SamAlignment, anchor_len: int, max_addition_len: int) -> dict[str, int | str] | None:
    read_len = len(aln.sequence)
    addition_len = read_len - anchor_len
    if addition_len < 1 or addition_len > max_addition_len:
        return None

    mismatch_offsets = {offset for _, _, _, offset in aln.mismatches}
    if aln.is_reverse:
        anchor_offsets = range(0, anchor_len)
        addition_offsets = range(anchor_len, read_len)
        addition_sequence = aln.sequence[anchor_len:].translate(srna.DNA_COMPLEMENT)[::-1].upper()
    else:
        anchor_offsets = range(addition_len, read_len)
        addition_offsets = range(0, addition_len)
        addition_sequence = aln.sequence[:addition_len].upper()

    if any(offset in mismatch_offsets for offset in anchor_offsets):
        return None
    if not all(offset in mismatch_offsets for offset in addition_offsets):
        return None
    if not set(addition_sequence) <= set(BASES):
        return None

    return {
        "five_prime_length": addition_len,
        "five_prime_sequence": addition_sequence,
        "five_prime_is_homopolymer": int(len(set(addition_sequence)) == 1),
        "five_prime_base": addition_sequence[0] if len(set(addition_sequence)) == 1 else "mixed",
    }


def init_substitution_counters(contig_lengths: dict[str, int]):
    base_exposure: Counter[tuple[str, str]] = Counter()
    substitution_counts: Counter[tuple[str, str]] = Counter()
    pos_coverage = {contig: np.zeros(length, dtype=np.float64) for contig, length in contig_lengths.items()}
    pos_subs = {
        contig: {sub: np.zeros(length, dtype=np.float64) for sub in SUBSTITUTIONS}
        for contig, length in contig_lengths.items()
    }
    return base_exposure, substitution_counts, pos_coverage, pos_subs


def init_strand_position_counters(contig_lengths: dict[str, int]):
    pos_coverage = {
        contig: {strand: np.zeros(length, dtype=np.float64) for strand in ("sense", "antisense")}
        for contig, length in contig_lengths.items()
    }
    pos_subs = {
        contig: {
            strand: {sub: np.zeros(length, dtype=np.float64) for sub in SUBSTITUTIONS}
            for strand in ("sense", "antisense")
        }
        for contig, length in contig_lengths.items()
    }
    return pos_coverage, pos_subs


def add_alignment_to_substitution_counts(
    aln: SamAlignment,
    weight: float,
    base_exposure: Counter[tuple[str, str]],
    substitution_counts: Counter[tuple[str, str]],
    pos_coverage: dict[str, np.ndarray],
    pos_subs: dict[str, dict[str, np.ndarray]],
    contigs: dict[str, srna.Contig],
    orientation: str = "template",
) -> None:
    template = contigs[aln.contig].sequence[aln.start0 : aln.end0]
    for offset, ref_base in enumerate(template):
        if ref_base not in BASES:
            continue
        pos0 = aln.start0 + offset
        exposure_base = complement_base(ref_base) if orientation == "molecule" and aln.is_reverse else ref_base
        base_exposure[(aln.contig, exposure_base)] += weight
        pos_coverage[aln.contig][pos0] += weight
    for pos0, ref_base, read_base, _ in aln.mismatches:
        substitution = oriented_substitution(aln, ref_base, read_base, orientation)
        substitution_counts[(aln.contig, substitution)] += weight
        pos_subs[aln.contig][substitution][pos0] += weight


def add_alignment_to_strand_position_counts(
    aln: SamAlignment,
    weight: float,
    pos_coverage: dict[str, dict[str, np.ndarray]],
    pos_subs: dict[str, dict[str, dict[str, np.ndarray]]],
    contigs: dict[str, srna.Contig],
    orientation: str = "template",
) -> None:
    strand = "antisense" if aln.is_reverse else "sense"
    template = contigs[aln.contig].sequence[aln.start0 : aln.end0]
    for offset, ref_base in enumerate(template):
        if ref_base not in BASES:
            continue
        pos0 = aln.start0 + offset
        pos_coverage[aln.contig][strand][pos0] += weight
    for pos0, ref_base, read_base, _ in aln.mismatches:
        substitution = oriented_substitution(aln, ref_base, read_base, orientation)
        pos_subs[aln.contig][strand][substitution][pos0] += weight


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def substitution_summary_rows(
    sample: str,
    contigs: dict[str, srna.Contig],
    base_exposure: Counter[tuple[str, str]],
    substitution_counts: Counter[tuple[str, str]],
) -> list[dict[str, str | float]]:
    rows = []
    for contig in contigs:
        for substitution in SUBSTITUTIONS:
            ref_base, read_base = substitution.split(">")
            denominator = float(base_exposure[(contig, ref_base)])
            count = float(substitution_counts[(contig, substitution)])
            rows.append(
                {
                    "sample": sample,
                    "contig": contig,
                    "substitution": substitution,
                    "template_base": ref_base,
                    "read_base": read_base,
                    "template_base_appearances_weighted": f"{denominator:.6f}",
                    "substitution_events_weighted": f"{count:.6f}",
                    "percent_of_template_base_appearance": f"{(count * 100.0 / denominator) if denominator else 0.0:.6f}",
                }
            )
    return rows


def substitution_position_rows(
    sample: str,
    contigs: dict[str, srna.Contig],
    pos_coverage: dict[str, np.ndarray],
    pos_subs: dict[str, dict[str, np.ndarray]],
) -> list[dict[str, str | int | float]]:
    rows = []
    for contig_name, contig in contigs.items():
        for pos0, ref_base in enumerate(contig.sequence):
            if ref_base not in BASES:
                continue
            denom = float(pos_coverage[contig_name][pos0])
            for substitution in SUBSTITUTIONS:
                if not substitution.startswith(f"{ref_base}>"):
                    continue
                count = float(pos_subs[contig_name][substitution][pos0])
                rows.append(
                    {
                        "sample": sample,
                        "contig": contig_name,
                        "position_1based": pos0 + 1,
                        "template_base": ref_base,
                        "substitution": substitution,
                        "coverage_weighted": f"{denom:.6f}",
                        "substitution_events_weighted": f"{count:.6f}",
                        "percent_at_position": f"{(count * 100.0 / denom) if denom else 0.0:.6f}",
                    }
                )
    return rows


def substitution_position_strand_rows(
    sample: str,
    contigs: dict[str, srna.Contig],
    pos_coverage: dict[str, dict[str, np.ndarray]],
    pos_subs: dict[str, dict[str, dict[str, np.ndarray]]],
    orientation: str = "template",
) -> list[dict[str, str | int | float]]:
    rows = []
    for contig_name, contig in contigs.items():
        for pos0, ref_base in enumerate(contig.sequence):
            if ref_base not in BASES:
                continue
            for substitution in SUBSTITUTIONS:
                sense_base = ref_base
                antisense_base = complement_base(ref_base) if orientation == "molecule" else ref_base
                sense_allowed = substitution.startswith(f"{sense_base}>")
                antisense_allowed = substitution.startswith(f"{antisense_base}>")
                if not sense_allowed and not antisense_allowed:
                    continue
                sense_cov = float(pos_coverage[contig_name]["sense"][pos0]) if sense_allowed else 0.0
                antisense_cov = float(pos_coverage[contig_name]["antisense"][pos0]) if antisense_allowed else 0.0
                sense_count = float(pos_subs[contig_name]["sense"][substitution][pos0])
                antisense_count = float(pos_subs[contig_name]["antisense"][substitution][pos0])
                rows.append(
                    {
                        "sample": sample,
                        "contig": contig_name,
                        "position_1based": pos0 + 1,
                        "template_base": ref_base,
                        "sense_reference_base": sense_base,
                        "antisense_reference_base": antisense_base,
                        "substitution": substitution,
                        "sense_coverage_weighted": f"{sense_cov:.6f}",
                        "sense_substitution_events_weighted": f"{sense_count:.6f}",
                        "sense_percent_at_position": f"{(sense_count * 100.0 / sense_cov) if sense_cov else 0.0:.6f}",
                        "antisense_coverage_weighted": f"{antisense_cov:.6f}",
                        "antisense_substitution_events_weighted": f"{antisense_count:.6f}",
                        "antisense_percent_at_position": f"{(antisense_count * 100.0 / antisense_cov) if antisense_cov else 0.0:.6f}",
                    }
                )
    return rows


def substitution_summary_strand_rows(
    sample: str,
    contigs: dict[str, srna.Contig],
    pos_coverage: dict[str, dict[str, np.ndarray]],
    pos_subs: dict[str, dict[str, dict[str, np.ndarray]]],
    orientation: str = "template",
) -> list[dict[str, str | float]]:
    rows = []
    for contig_name, contig in contigs.items():
        for strand in ("sense", "antisense"):
            for substitution in SUBSTITUTIONS:
                ref_base, read_base = substitution.split(">")
                denominator = 0.0
                count = 0.0
                for pos0, template_base in enumerate(contig.sequence):
                    if template_base not in BASES:
                        continue
                    oriented_base = complement_base(template_base) if orientation == "molecule" and strand == "antisense" else template_base
                    if oriented_base != ref_base:
                        continue
                    denominator += float(pos_coverage[contig_name][strand][pos0])
                    count += float(pos_subs[contig_name][strand][substitution][pos0])
                rows.append(
                    {
                        "sample": sample,
                        "contig": contig_name,
                        "strand": strand,
                        "substitution": substitution,
                        "reference_base": ref_base,
                        "read_base": read_base,
                        "reference_base_appearances_weighted": f"{denominator:.6f}",
                        "substitution_events_weighted": f"{count:.6f}",
                        "percent_of_reference_base_appearance": f"{(count * 100.0 / denominator) if denominator else 0.0:.6f}",
                    }
                )
    return rows


def aggregate_substitution_summary_strand_rows(rows: list[dict], sample_label: str = "overall") -> list[dict[str, str | float]]:
    counts: Counter[tuple[str, str, str, str, str]] = Counter()
    exposures: Counter[tuple[str, str, str, str, str]] = Counter()
    for row in rows:
        key = (
            str(row["contig"]),
            str(row["strand"]),
            str(row["substitution"]),
            str(row["reference_base"]),
            str(row["read_base"]),
        )
        counts[key] += float(row["substitution_events_weighted"])
        exposures[key] += float(row["reference_base_appearances_weighted"])

    out_rows = []
    for contig, strand, substitution, reference_base, read_base in sorted(counts):
        denominator = float(exposures[(contig, strand, substitution, reference_base, read_base)])
        count = float(counts[(contig, strand, substitution, reference_base, read_base)])
        out_rows.append(
            {
                "sample": sample_label,
                "contig": contig,
                "strand": strand,
                "substitution": substitution,
                "reference_base": reference_base,
                "read_base": read_base,
                "reference_base_appearances_weighted": f"{denominator:.6f}",
                "substitution_events_weighted": f"{count:.6f}",
                "percent_of_reference_base_appearance": f"{(count * 100.0 / denominator) if denominator else 0.0:.6f}",
            }
        )
    return out_rows


def mismatched_read_denominators(read_stat_rows: list[dict]) -> dict[str, float]:
    denominators: dict[str, float] = {}
    for row in read_stat_rows:
        if str(row["metric"]) == "mismatched_reads_analyzed":
            denominators[str(row["sample"])] = float(row["value"])
    return denominators


def substitution_summary_percent_of_mismatched_rows(summary_rows: list[dict], read_stat_rows: list[dict]) -> list[dict[str, str | float]]:
    denominators = mismatched_read_denominators(read_stat_rows)
    rows = []
    for row in summary_rows:
        sample = str(row["sample"])
        denominator = denominators.get(sample, 0.0)
        count = float(row["substitution_events_weighted"])
        out_row = dict(row)
        out_row["substitution_rna"] = substitution_label_rna(str(row["substitution"]))
        out_row["mismatched_sRNAs_denominator"] = f"{denominator:.6f}"
        out_row["percent_of_mismatched_sRNAs"] = f"{(count * 100.0 / denominator) if denominator else 0.0:.6f}"
        rows.append(out_row)
    return rows


def plot_substitution_bar(summary_rows: list[dict], out: Path, title: str) -> None:
    totals: Counter[str] = Counter()
    for row in summary_rows:
        sub = str(row["substitution"])
        totals[sub] += float(row["substitution_events_weighted"])
    by_ref_exposure: Counter[str] = Counter()
    seen_exposures: set[tuple[str, str, str]] = set()
    for row in summary_rows:
        key = (str(row["sample"]), str(row["contig"]), str(row["template_base"]))
        if key in seen_exposures:
            continue
        seen_exposures.add(key)
        by_ref_exposure[str(row["template_base"])] += float(row["template_base_appearances_weighted"])
    values = [totals[sub] * 100.0 / by_ref_exposure[sub[0]] if by_ref_exposure[sub[0]] else 0.0 for sub in SUBSTITUTIONS]

    fig, ax = plt.subplots(figsize=(10.5, 4.8), constrained_layout=True)
    colors = ["#1F77B4" if sub.startswith("A") else "#D62728" if sub.startswith("C") else "#2CA02C" if sub.startswith("G") else "#9467BD" for sub in SUBSTITUTIONS]
    ax.bar(np.arange(len(SUBSTITUTIONS)), values, color=colors)
    ax.set_xticks(np.arange(len(SUBSTITUTIONS)), SUBSTITUTION_LABELS_RNA, rotation=45, ha="right")
    ax.set_ylabel("% of template-base appearances")
    ax.set_title(title)
    ax.grid(axis="y", color="#E5E5E5", linewidth=0.7)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def plot_substitution_percent_of_mismatched(
    summary_rows: list[dict],
    out: Path,
    title: str,
    ylabel: str = "% of sRNAs with >=1 mismatch",
) -> None:
    samples = sorted({str(row["sample"]) for row in summary_rows})
    if not samples:
        return
    values: Counter[tuple[str, str]] = Counter()
    for row in summary_rows:
        values[(str(row["sample"]), str(row["substitution"]))] += float(row["percent_of_mismatched_sRNAs"])

    x = np.arange(len(SUBSTITUTIONS))
    colors = [
        "#1F77B4" if sub.startswith("A") else "#D62728" if sub.startswith("C") else "#2CA02C" if sub.startswith("G") else "#9467BD"
        for sub in SUBSTITUTIONS
    ]
    fig, ax = plt.subplots(figsize=(11, 5.0), constrained_layout=True)
    ymax = 0.0
    for idx, sub in enumerate(SUBSTITUTIONS):
        vals = np.array([values[(sample, sub)] for sample in samples], dtype=np.float64)
        mean_val = float(vals.mean()) if len(vals) else 0.0
        sd_val = float(vals.std(ddof=1)) if len(vals) > 1 else 0.0
        ax.bar(idx, mean_val, color=colors[idx], edgecolor="#303030", linewidth=0.35)
        ax.errorbar([idx], [mean_val], yerr=[sd_val], fmt="none", ecolor="#202020", elinewidth=1.0, capsize=3, zorder=5)
        jitter = np.linspace(-0.16, 0.16, len(samples)) if len(samples) > 1 else np.array([0.0])
        ax.scatter(idx + jitter, vals, color="#202020", s=18, alpha=0.78, zorder=6)
        ymax = max(ymax, mean_val + sd_val, *(vals.tolist() or [0.0]))

    ax.set_xticks(x, SUBSTITUTION_LABELS_RNA, rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(axis="y", color="#E5E5E5", linewidth=0.7)
    if ymax:
        ax.set_ylim(0, ymax * 1.18)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def plot_substitution_template_appearance_sample_bar(summary_rows: list[dict], out: Path, title: str) -> None:
    samples = sorted({str(row["sample"]) for row in summary_rows})
    if not samples:
        return
    values: Counter[tuple[str, str]] = Counter()
    for row in summary_rows:
        values[(str(row["sample"]), str(row["substitution"]))] += float(row["percent_of_template_base_appearance"])

    x = np.arange(len(SUBSTITUTIONS))
    colors = [
        "#1F77B4" if sub.startswith("A") else "#D62728" if sub.startswith("C") else "#2CA02C" if sub.startswith("G") else "#9467BD"
        for sub in SUBSTITUTIONS
    ]
    fig, ax = plt.subplots(figsize=(11, 5.0), constrained_layout=True)
    ymax = 0.0
    for idx, sub in enumerate(SUBSTITUTIONS):
        vals = np.array([values[(sample, sub)] for sample in samples], dtype=np.float64)
        mean_val = float(vals.mean()) if len(vals) else 0.0
        sd_val = float(vals.std(ddof=1)) if len(vals) > 1 else 0.0
        ax.bar(idx, mean_val, color=colors[idx], edgecolor="#303030", linewidth=0.35)
        ax.errorbar([idx], [mean_val], yerr=[sd_val], fmt="none", ecolor="#202020", elinewidth=1.0, capsize=3, zorder=5)
        jitter = np.linspace(-0.16, 0.16, len(samples)) if len(samples) > 1 else np.array([0.0])
        ax.scatter(idx + jitter, vals, color="#202020", s=18, alpha=0.78, zorder=6)
        ymax = max(ymax, mean_val + sd_val, *(vals.tolist() or [0.0]))

    ax.set_xticks(x, SUBSTITUTION_LABELS_RNA, rotation=45, ha="right")
    ax.set_ylabel("% of template nucleotide appearances")
    ax.set_title(title)
    ax.grid(axis="y", color="#E5E5E5", linewidth=0.7)
    if ymax:
        ax.set_ylim(0, ymax * 1.18)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def plot_substitution_percent_of_mismatched_by_strand(
    summary_rows: list[dict],
    out: Path,
    title: str,
    ylabel: str = "% of sRNAs with >=1 mismatch",
) -> None:
    samples = sorted({str(row["sample"]) for row in summary_rows})
    if not samples:
        return
    values: Counter[tuple[str, str, str]] = Counter()
    for row in summary_rows:
        values[(str(row["sample"]), str(row["strand"]), str(row["substitution"]))] += float(row["percent_of_mismatched_sRNAs"])

    x = np.arange(len(SUBSTITUTIONS))
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2), sharey=True, constrained_layout=True)
    colors = ["#0072B2" if sub.startswith(("A", "G")) else "#D55E00" for sub in SUBSTITUTIONS]
    ymax = 0.0
    for ax, strand in zip(axes, ("sense", "antisense"), strict=True):
        for idx, sub in enumerate(SUBSTITUTIONS):
            vals = np.array([values[(sample, strand, sub)] for sample in samples], dtype=np.float64)
            mean_val = float(vals.mean()) if len(vals) else 0.0
            sd_val = float(vals.std(ddof=1)) if len(vals) > 1 else 0.0
            ax.bar(idx, mean_val, color=colors[idx], edgecolor="#303030", linewidth=0.35)
            ax.errorbar([idx], [mean_val], yerr=[sd_val], fmt="none", ecolor="#202020", elinewidth=1.0, capsize=3, zorder=5)
            jitter = np.linspace(-0.16, 0.16, len(samples)) if len(samples) > 1 else np.array([0.0])
            ax.scatter(idx + jitter, vals, color="#202020", s=18, alpha=0.78, zorder=6)
            ymax = max(ymax, mean_val + sd_val, *(vals.tolist() or [0.0]))
        ax.set_title(strand)
        ax.set_xticks(x, SUBSTITUTION_LABELS_RNA, rotation=45, ha="right")
        ax.grid(axis="y", color="#E5E5E5", linewidth=0.7)
    axes[0].set_ylabel(ylabel)
    if ymax:
        axes[0].set_ylim(0, ymax * 1.18)
    fig.suptitle(title, fontsize=14)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def plot_substitution_strand_bar(summary_rows: list[dict], out: Path, title: str) -> None:
    samples = sorted({str(row["sample"]) for row in summary_rows})
    events: Counter[tuple[str, str, str]] = Counter()
    exposures: Counter[tuple[str, str, str]] = Counter()
    seen_exposures: set[tuple[str, str, str, str]] = set()
    for row in summary_rows:
        sample = str(row["sample"])
        strand = str(row["strand"])
        sub = str(row["substitution"])
        ref_base = str(row["reference_base"])
        events[(sample, strand, sub)] += float(row["substitution_events_weighted"])
        key = (sample, str(row["contig"]), strand, ref_base)
        if key in seen_exposures:
            continue
        seen_exposures.add(key)
        exposures[(sample, strand, ref_base)] += float(row["reference_base_appearances_weighted"])

    x = np.arange(len(SUBSTITUTIONS))
    fig, ax = plt.subplots(figsize=(11, 5.2), constrained_layout=True)
    ymax = 0.0
    jitter = np.linspace(-0.14, 0.14, len(samples)) if len(samples) > 1 else np.array([0.0])
    for idx, sub in enumerate(SUBSTITUTIONS):
        ref_base = sub[0]
        sense_vals = np.array(
            [
                events[(sample, "sense", sub)] * 100.0 / exposures[(sample, "sense", ref_base)]
                if exposures[(sample, "sense", ref_base)]
                else 0.0
                for sample in samples
            ],
            dtype=np.float64,
        )
        antisense_vals = np.array(
            [
                events[(sample, "antisense", sub)] * 100.0 / exposures[(sample, "antisense", ref_base)]
                if exposures[(sample, "antisense", ref_base)]
                else 0.0
                for sample in samples
            ],
            dtype=np.float64,
        )
        sense_mean = float(sense_vals.mean()) if len(sense_vals) else 0.0
        sense_sd = float(sense_vals.std(ddof=1)) if len(sense_vals) > 1 else 0.0
        antisense_mean = float(antisense_vals.mean()) if len(antisense_vals) else 0.0
        antisense_sd = float(antisense_vals.std(ddof=1)) if len(antisense_vals) > 1 else 0.0

        ax.bar(idx, sense_mean, width=0.72, color="#0072B2", label="sense" if idx == 0 else None)
        ax.errorbar([idx], [sense_mean], yerr=[sense_sd], fmt="none", ecolor="#202020", elinewidth=1.0, capsize=3, zorder=5)
        ax.scatter(idx + jitter, sense_vals, color="#202020", s=15, alpha=0.72, zorder=6)

        ax.bar(idx, -antisense_mean, width=0.72, color="#D55E00", label="antisense" if idx == 0 else None)
        ax.errorbar([idx], [-antisense_mean], yerr=[antisense_sd], fmt="none", ecolor="#202020", elinewidth=1.0, capsize=3, zorder=5)
        ax.scatter(idx + jitter, -antisense_vals, color="#202020", s=15, alpha=0.72, zorder=6)

        ymax = max(ymax, sense_mean + sense_sd, antisense_mean + antisense_sd, *(sense_vals.tolist() or [0.0]), *(antisense_vals.tolist() or [0.0]))

    ax.axhline(0, color="#303030", linewidth=0.8)
    ax.set_xticks(x, SUBSTITUTION_LABELS_RNA, rotation=45, ha="right")
    ax.set_ylabel("% rate\n+sense / -antisense")
    ax.set_title(title)
    ax.grid(axis="y", color="#E5E5E5", linewidth=0.7)
    ax.legend(frameon=False, ncol=2)
    if ymax:
        ax.set_ylim(-ymax * 1.18, ymax * 1.18)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def collapsed_rc_mutation_family_rows(summary_rows: list[dict]) -> list[dict[str, str | float]]:
    events: Counter[tuple[str, str, str]] = Counter()
    exposures: Counter[tuple[str, str, str]] = Counter()
    for row in summary_rows:
        family = RC_MUTATION_FAMILIES[str(row["substitution"])]
        sample = str(row["sample"])
        strand = str(row["strand"])
        events[(sample, strand, family)] += float(row["substitution_events_weighted"])
        exposures[(sample, strand, family)] += float(row["reference_base_appearances_weighted"])

    samples = sorted({str(row["sample"]) for row in summary_rows})
    out_rows = []
    for sample in samples:
        for strand in ("sense", "antisense"):
            for family in RC_MUTATION_FAMILY_ORDER:
                count = float(events[(sample, strand, family)])
                denominator = float(exposures[(sample, strand, family)])
                out_rows.append(
                    {
                        "sample": sample,
                        "strand": strand,
                        "mutation_family": family,
                        "substitution_events_weighted": f"{count:.6f}",
                        "reference_base_appearances_weighted": f"{denominator:.6f}",
                        "percent_rate": f"{(count * 100.0 / denominator) if denominator else 0.0:.6f}",
                    }
                )
    return out_rows


def plot_collapsed_rc_mutation_family_strand_bar(rows: list[dict], out: Path, title: str) -> None:
    samples = sorted({str(row["sample"]) for row in rows})
    if not samples:
        return
    values: Counter[tuple[str, str, str]] = Counter()
    for row in rows:
        values[(str(row["sample"]), str(row["strand"]), str(row["mutation_family"]))] = float(row["percent_rate"])

    x = np.arange(len(RC_MUTATION_FAMILY_ORDER))
    fig, ax = plt.subplots(figsize=(9.8, 5.0), constrained_layout=True)
    jitter = np.linspace(-0.13, 0.13, len(samples)) if len(samples) > 1 else np.array([0.0])
    ymax = 0.0
    for idx, family in enumerate(RC_MUTATION_FAMILY_ORDER):
        sense_vals = np.array([values[(sample, "sense", family)] for sample in samples], dtype=np.float64)
        antisense_vals = np.array([values[(sample, "antisense", family)] for sample in samples], dtype=np.float64)
        sense_mean = float(sense_vals.mean()) if len(sense_vals) else 0.0
        antisense_mean = float(antisense_vals.mean()) if len(antisense_vals) else 0.0
        sense_sd = float(sense_vals.std(ddof=1)) if len(sense_vals) > 1 else 0.0
        antisense_sd = float(antisense_vals.std(ddof=1)) if len(antisense_vals) > 1 else 0.0
        ax.bar(idx, sense_mean, width=0.70, color="#0072B2", label="sense" if idx == 0 else None)
        ax.errorbar([idx], [sense_mean], yerr=[sense_sd], fmt="none", ecolor="#202020", elinewidth=1.0, capsize=3, zorder=5)
        ax.scatter(idx + jitter, sense_vals, color="#202020", s=18, alpha=0.72, zorder=6)
        ax.bar(idx, -antisense_mean, width=0.70, color="#D55E00", label="antisense" if idx == 0 else None)
        ax.errorbar([idx], [-antisense_mean], yerr=[antisense_sd], fmt="none", ecolor="#202020", elinewidth=1.0, capsize=3, zorder=5)
        ax.scatter(idx + jitter, -antisense_vals, color="#202020", s=18, alpha=0.72, zorder=6)
        ymax = max(
            ymax,
            sense_mean + sense_sd,
            antisense_mean + antisense_sd,
            *(sense_vals.tolist() or [0.0]),
            *(antisense_vals.tolist() or [0.0]),
        )

    ax.axhline(0, color="#303030", linewidth=0.8)
    ax.set_xticks(x, RC_MUTATION_FAMILY_ORDER, rotation=30, ha="right")
    ax.set_ylabel("% rate\n+sense / -antisense")
    ax.set_title(title)
    ax.grid(axis="y", color="#E5E5E5", linewidth=0.7)
    ax.legend(frameon=False, ncol=2)
    if ymax:
        ax.set_ylim(-ymax * 1.2, ymax * 1.2)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def plot_position_heatmap(position_rows: list[dict], contig: srna.Contig, out: Path, title: str) -> None:
    matrix = np.zeros((len(SUBSTITUTIONS), len(contig.sequence)), dtype=np.float64)
    sub_to_idx = {sub: idx for idx, sub in enumerate(SUBSTITUTIONS)}
    for row in position_rows:
        sub = str(row["substitution"])
        pos0 = int(row["position_1based"]) - 1
        matrix[sub_to_idx[sub], pos0] = float(row["percent_at_position"])

    fig, ax = plt.subplots(figsize=(14, 4.8), constrained_layout=True)
    vmax = max(1.0, float(np.percentile(matrix[matrix > 0], 99)) if np.any(matrix > 0) else 1.0)
    im = ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap="magma", vmin=0, vmax=vmax)
    ax.set_yticks(np.arange(len(SUBSTITUTIONS)), SUBSTITUTION_LABELS_RNA)
    ax.set_xlabel(f"position on {contig.name}")
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax, pad=0.01)
    cbar.set_label("% at position")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def plot_position_lollipop_grid(position_rows: list[dict], contig: srna.Contig, out: Path, title: str) -> None:
    grouped: dict[str, dict[int, Counter[str]]] = {
        sub: {pos: Counter() for pos in range(1, len(contig.sequence) + 1)}
        for sub in SUBSTITUTIONS
    }
    for row in position_rows:
        if row["contig"] != contig.name:
            continue
        sub = str(row["substitution"])
        pos = int(row["position_1based"])
        grouped[sub][pos]["sense_events"] += float(row["sense_substitution_events_weighted"])
        grouped[sub][pos]["sense_coverage"] += float(row["sense_coverage_weighted"])
        grouped[sub][pos]["antisense_events"] += float(row["antisense_substitution_events_weighted"])
        grouped[sub][pos]["antisense_coverage"] += float(row["antisense_coverage_weighted"])

    fig, axes = plt.subplots(4, 3, figsize=(18, 12), constrained_layout=True)
    sense_color = "#0072B2"
    antisense_color = "#D55E00"
    for ax, sub in zip(axes.flat, SUBSTITUTIONS, strict=True):
        sense_x, sense_y, anti_x, anti_y = [], [], [], []
        for pos in range(1, len(contig.sequence) + 1):
            counts = grouped[sub][pos]
            if counts["sense_coverage"]:
                pct = counts["sense_events"] * 100.0 / counts["sense_coverage"]
                if pct:
                    sense_x.append(pos)
                    sense_y.append(pct)
            if counts["antisense_coverage"]:
                pct = counts["antisense_events"] * 100.0 / counts["antisense_coverage"]
                if pct:
                    anti_x.append(pos)
                    anti_y.append(-pct)

        if sense_x:
            ax.vlines(sense_x, 0, sense_y, color=sense_color, alpha=0.35, linewidth=0.8)
            ax.scatter(sense_x, sense_y, color=sense_color, s=13, label="sense", zorder=3)
        if anti_x:
            ax.vlines(anti_x, 0, anti_y, color=antisense_color, alpha=0.35, linewidth=0.8)
            ax.scatter(anti_x, anti_y, color=antisense_color, s=13, label="antisense", zorder=3)

        ymax = max([abs(y) for y in sense_y + anti_y] + [0.1])
        ax.set_ylim(-ymax * 1.18, ymax * 1.18)
        ax.set_xlim(1, len(contig.sequence))
        ax.axhline(0, color="#303030", linewidth=0.8)
        ax.set_title(substitution_label_rna(sub), fontsize=11)
        ax.grid(color="#E8E8E8", linewidth=0.6)

    for row_axes in axes:
        row_axes[0].set_ylabel("% rate\n+sense / -antisense")
    for ax in axes[-1, :]:
        ax.set_xlabel("position")
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", label="sense", markerfacecolor=sense_color, markersize=7),
        plt.Line2D([0], [0], marker="o", color="w", label="antisense", markerfacecolor=antisense_color, markersize=7),
    ]
    fig.suptitle(title, fontsize=16)
    fig.legend(handles=handles, loc="upper right", frameon=False, ncol=2)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def plot_tailing_summary(rows: list[dict], out: Path) -> None:
    counts: Counter[tuple[str, str]] = Counter()
    for row in rows:
        counts[(str(row["sample"]), str(row["tail_sequence"]))] += float(row["percent_of_total_sRNAs"])
    samples = sorted({sample for sample, _ in counts})
    top_tails = [tail for tail, _ in Counter({tail: sum(v for (s, t), v in counts.items() if t == tail) for _, tail in counts}).most_common(12)]
    if not samples or not top_tails:
        return
    bottom = np.zeros(len(samples))
    fig, ax = plt.subplots(figsize=(11, 5.2), constrained_layout=True)
    palette = plt.get_cmap("tab20")
    x = np.arange(len(samples))
    for idx, tail in enumerate(top_tails):
        vals = np.array([counts[(sample, tail)] for sample in samples])
        ax.bar(x, vals, bottom=bottom, color=palette(idx % 20), label=tail)
        bottom += vals
    ax.set_xticks(x, samples, rotation=35, ha="right")
    ax.set_ylabel("% of total sRNAs")
    ax.set_title("Candidate 3-prime tailing sequences")
    ax.grid(axis="y", color="#E5E5E5", linewidth=0.7)
    ax.legend(frameon=False, ncol=4, fontsize=8)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def plot_tailing_summary_by_strand(
    rows: list[dict],
    out: Path,
    value_column: str = "percent_of_total_sRNAs",
    ylabel: str = "% of total sRNAs",
    title: str = "Candidate 3-prime tailing sequences by strand",
) -> None:
    counts: Counter[tuple[str, str, str]] = Counter()
    for row in rows:
        counts[(str(row["sample"]), str(row["strand"]), str(row["tail_sequence"]))] += float(row[value_column])
    samples = sorted({sample for sample, _, _ in counts})
    tail_totals: Counter[str] = Counter()
    for (_, _, tail), count in counts.items():
        tail_totals[tail] += count
    tails = [tail for tail, _ in tail_totals.most_common(10)]
    if not samples or not tails:
        return

    labels = [f"{sample}\n{strand}" for sample in samples for strand in ("+", "-")]
    x = np.arange(len(labels))
    bottom = np.zeros(len(labels))
    fig, ax = plt.subplots(figsize=(12, 5.4), constrained_layout=True)
    palette = plt.get_cmap("tab20")
    for idx, tail in enumerate(tails):
        vals = np.array([counts[(sample, strand, tail)] for sample in samples for strand in ("+", "-")])
        ax.bar(x, vals, bottom=bottom, color=palette(idx % 20), label=tail)
        bottom += vals
    ax.set_xticks(x, labels, rotation=35, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(axis="y", color="#E5E5E5", linewidth=0.7)
    ax.legend(frameon=False, ncol=5, fontsize=8)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


TAIL_GROUPS = {
    "U": ["T", "TT", "TTT"],
    "A": ["A", "AA"],
    "C": ["C"],
    "G": ["G"],
}


def grouped_tail_label(tail_sequence: str) -> str | None:
    for label, members in TAIL_GROUPS.items():
        if tail_sequence in members:
            return label
    return None


def rna_tail_label(tail_sequence: str) -> str:
    return tail_sequence.replace("T", "U")


def grouped_terminal_component_summary_rows(
    rows: list[dict],
    sequence_key: str,
    count_key: str,
) -> list[dict[str, str | float]]:
    counts: Counter[tuple[str, str, str, str]] = Counter()
    denominators: dict[tuple[str, str], tuple[float, float]] = {}
    for row in rows:
        sample = str(row["sample"])
        strand = str(row["strand"])
        sequence = str(row[sequence_key])
        group = grouped_tail_label(sequence)
        denominators[(sample, strand)] = (
            float(row["total_sRNAs_denominator"]),
            float(row["mismatched_sRNAs_denominator"]),
        )
        if group is None:
            continue
        counts[(sample, strand, group, sequence)] += float(row[count_key])

    out_rows = []
    for sample, strand in sorted(denominators):
        total_denominator, mismatched_denominator = denominators[(sample, strand)]
        for group, sequences in TAIL_GROUPS.items():
            for sequence in sequences:
                count = float(counts[(sample, strand, group, sequence)])
                out_rows.append(
                    {
                        "sample": sample,
                        "strand": strand,
                        "tail_group": group,
                        "tail_component": rna_tail_label(sequence),
                        "tail_component_dna": sequence,
                        "tail_count_weighted": f"{count:.6f}",
                        "total_sRNAs_denominator": f"{total_denominator:.6f}",
                        "percent_of_total_sRNAs": f"{(count * 100.0 / total_denominator) if total_denominator else 0.0:.9f}",
                        "mismatched_sRNAs_denominator": f"{mismatched_denominator:.6f}",
                        "percent_of_mismatched_sRNAs": f"{(count * 100.0 / mismatched_denominator) if mismatched_denominator else 0.0:.6f}",
                    }
                )
    return out_rows


def grouped_tail_summary_rows(rows: list[dict]) -> list[dict[str, str | float]]:
    counts: Counter[tuple[str, str, str]] = Counter()
    denominators: dict[tuple[str, str], tuple[float, float]] = {}
    for row in rows:
        sample = str(row["sample"])
        strand = str(row["strand"])
        group = grouped_tail_label(str(row["tail_sequence"]))
        denominators[(sample, strand)] = (
            float(row["total_sRNAs_denominator"]),
            float(row["mismatched_sRNAs_denominator"]),
        )
        if group is None:
            continue
        counts[(sample, strand, group)] += float(row["tail_count_weighted"])

    out_rows = []
    for sample, strand in sorted(denominators):
        total_denominator, mismatched_denominator = denominators[(sample, strand)]
        for group in TAIL_GROUPS:
            count = float(counts[(sample, strand, group)])
            out_rows.append(
                {
                    "sample": sample,
                    "strand": strand,
                    "tail_group": group,
                    "grouped_tail_members": ",".join(
                        rna_tail_label(value) for value in sorted(TAIL_GROUPS[group], key=lambda value: (len(value), value))
                    ),
                    "tail_count_weighted": f"{count:.6f}",
                    "total_sRNAs_denominator": f"{total_denominator:.6f}",
                    "percent_of_total_sRNAs": f"{(count * 100.0 / total_denominator) if total_denominator else 0.0:.9f}",
                    "mismatched_sRNAs_denominator": f"{mismatched_denominator:.6f}",
                    "percent_of_mismatched_sRNAs": f"{(count * 100.0 / mismatched_denominator) if mismatched_denominator else 0.0:.6f}",
                }
            )
    return out_rows


def plot_grouped_tail_summary_by_strand(
    rows: list[dict],
    out: Path,
    value_column: str = "percent_of_mismatched_sRNAs",
    ylabel: str = "% of sRNAs with >=1 mismatch",
    title: str = "Grouped candidate 3-prime tailing among mismatched sRNAs",
) -> None:
    samples = sorted({str(row["sample"]) for row in rows})
    if not samples:
        return
    groups = list(TAIL_GROUPS)
    values: Counter[tuple[str, str, str]] = Counter()
    for row in rows:
        values[(str(row["sample"]), str(row["strand"]), str(row["tail_group"]))] = float(row[value_column])

    x = np.arange(len(samples))
    bar_width = 0.10
    offsets = {
        ("+", "U"): -3.5 * bar_width,
        ("+", "A"): -2.5 * bar_width,
        ("+", "C"): -1.5 * bar_width,
        ("+", "G"): -0.5 * bar_width,
        ("-", "U"): 0.5 * bar_width,
        ("-", "A"): 1.5 * bar_width,
        ("-", "C"): 2.5 * bar_width,
        ("-", "G"): 3.5 * bar_width,
    }
    colors = {"U": "#0072B2", "A": "#D55E00", "C": "#009E73", "G": "#CC79A7"}
    hatches = {"+": "", "-": "//"}

    fig, ax = plt.subplots(figsize=(12, 5.5), constrained_layout=True)
    for strand in ("+", "-"):
        for group in groups:
            vals = np.array([values[(sample, strand, group)] for sample in samples])
            ax.bar(
                x + offsets[(strand, group)],
                vals,
                width=bar_width,
                color=colors[group],
                hatch=hatches[strand],
                edgecolor="#303030",
                linewidth=0.3,
                label=f"{strand} {group}",
            )
    ax.set_xticks(x, samples, rotation=25, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(axis="y", color="#E5E5E5", linewidth=0.7)
    ax.legend(frameon=False, ncol=4, fontsize=8)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def plot_terminal_addition_comparison(
    three_prime_rows: list[dict],
    five_prime_rows: list[dict],
    out: Path,
    value_column: str = "percent_of_mismatched_sRNAs",
) -> None:
    datasets = [("3-prime tailing", three_prime_rows), ("5-prime additions", five_prime_rows)]
    samples = sorted({str(row["sample"]) for _, rows in datasets for row in rows})
    if not samples:
        return

    component_colors = {
        "U": "#0072B2",
        "UU": "#56B4E9",
        "UUU": "#A6CEE3",
        "A": "#D55E00",
        "AA": "#E69F00",
        "C": "#009E73",
        "G": "#CC79A7",
    }
    component_order = {group: [rna_tail_label(sequence) for sequence in sequences] for group, sequences in TAIL_GROUPS.items()}
    categories = [(strand, group) for strand in ("+", "-") for group in TAIL_GROUPS]
    x = np.arange(len(categories))

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.0), sharey=True, constrained_layout=True)
    ymax = 0.0
    for ax, (panel_title, rows) in zip(axes, datasets, strict=True):
        values: Counter[tuple[str, str, str, str]] = Counter()
        totals_by_sample: dict[tuple[str, str, str], float] = defaultdict(float)
        for row in rows:
            sample = str(row["sample"])
            strand = str(row["strand"])
            group = str(row["tail_group"])
            component = str(row["tail_component"])
            value = float(row[value_column])
            values[(sample, strand, group, component)] += value
            totals_by_sample[(sample, strand, group)] += value

        for xpos, (strand, group) in enumerate(categories):
            bottom = 0.0
            for component in component_order[group]:
                component_vals = np.array([values[(sample, strand, group, component)] for sample in samples], dtype=np.float64)
                mean_val = float(component_vals.mean()) if len(component_vals) else 0.0
                if mean_val:
                    ax.bar(
                        xpos,
                        mean_val,
                        bottom=bottom,
                        width=0.70,
                        color=component_colors[component],
                        edgecolor="#303030",
                        linewidth=0.35,
                        label=component,
                    )
                bottom += mean_val

            total_vals = np.array([totals_by_sample[(sample, strand, group)] for sample in samples], dtype=np.float64)
            mean_total = float(total_vals.mean()) if len(total_vals) else 0.0
            sd_total = float(total_vals.std(ddof=1)) if len(total_vals) > 1 else 0.0
            ax.errorbar(
                [xpos],
                [mean_total],
                yerr=[sd_total],
                fmt="none",
                ecolor="#202020",
                elinewidth=1.0,
                capsize=3,
                zorder=5,
            )
            jitter = np.linspace(-0.16, 0.16, len(samples)) if len(samples) > 1 else np.array([0.0])
            ax.scatter(
                xpos + jitter,
                total_vals,
                color="#202020",
                s=18,
                alpha=0.78,
                zorder=6,
            )
            ymax = max(ymax, mean_total + sd_total, *(total_vals.tolist() or [0.0]))

        ax.set_title(panel_title)
        ax.set_xticks(x, [group for _, group in categories])
        ax.text(1.5, 1.015, "sense", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=11)
        ax.text(5.5, 1.015, "antisense", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=11)
        ax.grid(axis="y", color="#E5E5E5", linewidth=0.7)
        ax.axvline(3.5, color="#BDBDBD", linewidth=0.8)

    axes[0].set_ylabel("% of sRNAs with >=1 mismatch")
    handles = []
    seen = set()
    for group in TAIL_GROUPS:
        for component in component_order[group]:
            if component in seen:
                continue
            seen.add(component)
            handles.append(plt.Rectangle((0, 0), 1, 1, color=component_colors[component], label=component))
    axes[1].legend(handles=handles, loc="upper right", ncol=2, frameon=False, title="components")
    if ymax:
        axes[0].set_ylim(0, ymax * 1.25)
    fig.suptitle("Terminal additions grouped by base", fontsize=14)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def mismatch_summary_strand_rows(
    sample: str,
    mismatch_counts: Counter[tuple[str, int]],
    strand_denominators: Counter[str],
    max_mismatches: int,
) -> list[dict[str, str | int | float]]:
    rows = []
    for strand in ("sense", "antisense"):
        denominator = float(strand_denominators[strand])
        for nm in range(0, max_mismatches + 1):
            count = float(mismatch_counts[(strand, nm)])
            rows.append(
                {
                    "sample": sample,
                    "strand": strand,
                    "best_mismatches": nm,
                    "reads_weighted": f"{count:.6f}",
                    "strand_mapped_reads_weighted": f"{denominator:.6f}",
                    "percent_within_strand_mapped": f"{(count * 100.0 / denominator) if denominator else 0.0:.6f}",
                }
            )
    return rows


def plot_mismatch_summary_by_strand(rows: list[dict], out: Path) -> None:
    samples = sorted({str(row["sample"]) for row in rows})
    if not samples:
        return
    values: Counter[tuple[str, str, int]] = Counter()
    for row in rows:
        values[(str(row["sample"]), str(row["strand"]), int(row["best_mismatches"]))] = float(row["percent_within_strand_mapped"])

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8), sharey=True, constrained_layout=True)
    colors = {0: "#0072B2", 1: "#009E73", 2: "#E69F00", 3: "#D55E00"}
    ymax = 0.0
    for ax, strand in zip(axes, ("sense", "antisense"), strict=True):
        x = np.arange(4)
        for nm in range(4):
            vals = np.array([values[(sample, strand, nm)] for sample in samples], dtype=np.float64)
            mean_val = float(vals.mean()) if len(vals) else 0.0
            sd_val = float(vals.std(ddof=1)) if len(vals) > 1 else 0.0
            ax.bar(nm, mean_val, width=0.68, color=colors[nm], edgecolor="#303030", linewidth=0.35)
            ax.errorbar([nm], [mean_val], yerr=[sd_val], fmt="none", ecolor="#202020", elinewidth=1.0, capsize=3, zorder=5)
            jitter = np.linspace(-0.14, 0.14, len(samples)) if len(samples) > 1 else np.array([0.0])
            ax.scatter(nm + jitter, vals, color="#202020", s=18, alpha=0.78, zorder=6)
            ymax = max(ymax, mean_val + sd_val, *(vals.tolist() or [0.0]))
        ax.set_title(strand)
        ax.set_xticks(x, ["0", "1", "2", "3"])
        ax.set_xlabel("best mismatches")
        ax.grid(axis="y", color="#E5E5E5", linewidth=0.7)
    axes[0].set_ylabel("% within strand-mapped sRNAs")
    if ymax:
        axes[0].set_ylim(0, min(100.0, ymax * 1.18))
    fig.suptitle("Best template-mismatch ratios by strand", fontsize=14)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def plot_five_prime_summary_by_strand(rows: list[dict], out: Path) -> None:
    counts: Counter[tuple[str, str, str]] = Counter()
    for row in rows:
        counts[(str(row["sample"]), str(row["strand"]), str(row["five_prime_sequence"]))] += float(row["percent_of_total_sRNAs"])
    samples = sorted({sample for sample, _, _ in counts})
    sequence_totals: Counter[str] = Counter()
    for (_, _, sequence), count in counts.items():
        sequence_totals[sequence] += count
    sequences = [sequence for sequence, _ in sequence_totals.most_common(10)]
    if not samples or not sequences:
        return

    labels = [f"{sample}\n{strand}" for sample in samples for strand in ("+", "-")]
    x = np.arange(len(labels))
    bottom = np.zeros(len(labels))
    fig, ax = plt.subplots(figsize=(12, 5.4), constrained_layout=True)
    palette = plt.get_cmap("tab20")
    for idx, sequence in enumerate(sequences):
        vals = np.array([counts[(sample, strand, sequence)] for sample in samples for strand in ("+", "-")])
        ax.bar(x, vals, bottom=bottom, color=palette(idx % 20), label=sequence)
        bottom += vals
    ax.set_xticks(x, labels, rotation=35, ha="right")
    ax.set_ylabel("% of total sRNAs")
    ax.set_title("Candidate 5-prime addition sequences by strand")
    ax.grid(axis="y", color="#E5E5E5", linewidth=0.7)
    ax.legend(frameon=False, ncol=5, fontsize=8)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def analyze_sample(sample: str, sam_path: Path, contigs: dict[str, srna.Contig], args: argparse.Namespace) -> dict[str, list[dict]]:
    contig_lengths = {name: len(contig.sequence) for name, contig in contigs.items()}
    all_exposure, all_counts, all_pos_cov, all_pos_subs = init_substitution_counters(contig_lengths)
    nontail_exposure, nontail_counts, nontail_pos_cov, nontail_pos_subs = init_substitution_counters(contig_lengths)
    all_strand_pos_cov, all_strand_pos_subs = init_strand_position_counters(contig_lengths)
    nontail_strand_pos_cov, nontail_strand_pos_subs = init_strand_position_counters(contig_lengths)
    molecule_all_exposure, molecule_all_counts, molecule_all_pos_cov, molecule_all_pos_subs = init_substitution_counters(contig_lengths)
    molecule_nontail_exposure, molecule_nontail_counts, molecule_nontail_pos_cov, molecule_nontail_pos_subs = init_substitution_counters(contig_lengths)
    molecule_all_strand_pos_cov, molecule_all_strand_pos_subs = init_strand_position_counters(contig_lengths)
    molecule_nontail_strand_pos_cov, molecule_nontail_strand_pos_subs = init_strand_position_counters(contig_lengths)
    tail_rows: list[dict] = []
    five_prime_rows: list[dict] = []
    mismatch_strand_counts: Counter[tuple[str, int]] = Counter()
    mismatch_strand_denominators: Counter[str] = Counter()
    read_stats: Counter[str] = Counter()

    for group in iter_read_groups(sam_path, contigs):
        read_stats["sam_reads_examined"] += 1
        if not group.alignments:
            continue
        read_stats["template_mapped_reads"] += 1
        reason = low_complexity_reason(group.original_sequence, args)
        if reason:
            read_stats[f"mapped_filtered_{reason}"] += 1
            continue
        read_stats["template_mapped_reads_passing_complexity_filter"] += 1

        best_nm = min(aln.nm for aln in group.alignments)
        read_stats[f"best_nm_{best_nm}_reads"] += 1
        best_alignments = [aln for aln in group.alignments if aln.nm == best_nm]
        weight = 1.0 / len(best_alignments)
        if 0 <= best_nm <= args.max_mismatches:
            for aln in best_alignments:
                strand = "antisense" if aln.is_reverse else "sense"
                mismatch_strand_counts[(strand, best_nm)] += weight
                mismatch_strand_denominators[strand] += weight
        if best_nm < args.min_mismatches or best_nm > args.max_mismatches:
            continue
        read_stats["mismatched_reads_analyzed"] += 1

        group_has_tail = False
        tail_infos: list[tuple[SamAlignment, dict[str, int | str]]] = []
        group_has_five_prime = False
        five_prime_infos: list[tuple[SamAlignment, dict[str, int | str]]] = []
        for aln in best_alignments:
            tail_info = candidate_tail(aln, args.tail_anchor_len, args.max_tail_len)
            if tail_info is None and args.short_core_tail_len:
                tail_info = candidate_one_nt_3prime_tail_from_short_core(aln, args.short_core_tail_len)
            if tail_info is not None:
                group_has_tail = True
                tail_infos.append((aln, tail_info))
            five_prime_info = candidate_five_prime_addition(aln, args.tail_anchor_len, args.max_tail_len)
            if five_prime_info is not None:
                group_has_five_prime = True
                five_prime_infos.append((aln, five_prime_info))

        if group_has_tail:
            read_stats["candidate_tail_reads"] += 1
        if group_has_five_prime:
            read_stats["candidate_five_prime_reads"] += 1
        for aln, tail_info in tail_infos:
            tail_rows.append(
                {
                    "sample": sample,
                    "read_id": aln.qname,
                    "contig": aln.contig,
                    "strand": "-" if aln.is_reverse else "+",
                    "template_start_1based": aln.start0 + 1,
                    "template_end_1based": aln.end0,
                    "read_length": len(aln.sequence),
                    "best_nm": best_nm,
                    "alignment_weight": f"{weight:.6f}",
                    **tail_info,
                }
            )
        for aln, five_prime_info in five_prime_infos:
            five_prime_rows.append(
                {
                    "sample": sample,
                    "read_id": aln.qname,
                    "contig": aln.contig,
                    "strand": "-" if aln.is_reverse else "+",
                    "template_start_1based": aln.start0 + 1,
                    "template_end_1based": aln.end0,
                    "read_length": len(aln.sequence),
                    "best_nm": best_nm,
                    "alignment_weight": f"{weight:.6f}",
                    **five_prime_info,
                }
            )

        for aln in best_alignments:
            add_alignment_to_substitution_counts(aln, weight, all_exposure, all_counts, all_pos_cov, all_pos_subs, contigs)
            add_alignment_to_strand_position_counts(aln, weight, all_strand_pos_cov, all_strand_pos_subs, contigs)
            add_alignment_to_substitution_counts(
                aln,
                weight,
                molecule_all_exposure,
                molecule_all_counts,
                molecule_all_pos_cov,
                molecule_all_pos_subs,
                contigs,
                orientation="molecule",
            )
            add_alignment_to_strand_position_counts(
                aln,
                weight,
                molecule_all_strand_pos_cov,
                molecule_all_strand_pos_subs,
                contigs,
                orientation="molecule",
            )
            if not group_has_tail:
                add_alignment_to_substitution_counts(
                    aln,
                    weight,
                    nontail_exposure,
                    nontail_counts,
                    nontail_pos_cov,
                    nontail_pos_subs,
                    contigs,
                )
                add_alignment_to_strand_position_counts(
                    aln,
                    weight,
                    nontail_strand_pos_cov,
                    nontail_strand_pos_subs,
                    contigs,
                )
                add_alignment_to_substitution_counts(
                    aln,
                    weight,
                    molecule_nontail_exposure,
                    molecule_nontail_counts,
                    molecule_nontail_pos_cov,
                    molecule_nontail_pos_subs,
                    contigs,
                    orientation="molecule",
                )
                add_alignment_to_strand_position_counts(
                    aln,
                    weight,
                    molecule_nontail_strand_pos_cov,
                    molecule_nontail_strand_pos_subs,
                    contigs,
                    orientation="molecule",
                )

    read_stat_rows = [{"sample": sample, "metric": metric, "value": value} for metric, value in sorted(read_stats.items())]
    total_srna_denominator = float(read_stats["sam_reads_examined"])
    mismatched_srna_denominator = float(read_stats["mismatched_reads_analyzed"])
    tail_summary_counter: Counter[tuple[int, str, str]] = Counter()
    for row in tail_rows:
        tail_summary_counter[(int(row["tail_length"]), str(row["tail_sequence"]), str(row["tail_base"]))] += float(row["alignment_weight"])
    tail_summary_rows = [
        {
            "sample": sample,
            "tail_length": tail_len,
            "tail_sequence": tail_sequence,
            "tail_base": tail_base,
            "tail_count_weighted": f"{count:.6f}",
            "total_sRNAs_denominator": f"{total_srna_denominator:.6f}",
            "percent_of_total_sRNAs": f"{(count * 100.0 / total_srna_denominator) if total_srna_denominator else 0.0:.9f}",
            "mismatched_sRNAs_denominator": f"{mismatched_srna_denominator:.6f}",
            "percent_of_mismatched_sRNAs": f"{(count * 100.0 / mismatched_srna_denominator) if mismatched_srna_denominator else 0.0:.6f}",
        }
        for (tail_len, tail_sequence, tail_base), count in sorted(tail_summary_counter.items(), key=lambda item: (-item[1], item[0]))
    ]
    tail_summary_strand_counter: Counter[tuple[str, int, str, str]] = Counter()
    for row in tail_rows:
        tail_summary_strand_counter[
            (str(row["strand"]), int(row["tail_length"]), str(row["tail_sequence"]), str(row["tail_base"]))
        ] += float(row["alignment_weight"])
    tail_summary_strand_rows = [
        {
            "sample": sample,
            "strand": strand,
            "tail_length": tail_len,
            "tail_sequence": tail_sequence,
            "tail_base": tail_base,
            "tail_count_weighted": f"{count:.6f}",
            "total_sRNAs_denominator": f"{total_srna_denominator:.6f}",
            "percent_of_total_sRNAs": f"{(count * 100.0 / total_srna_denominator) if total_srna_denominator else 0.0:.9f}",
            "mismatched_sRNAs_denominator": f"{mismatched_srna_denominator:.6f}",
            "percent_of_mismatched_sRNAs": f"{(count * 100.0 / mismatched_srna_denominator) if mismatched_srna_denominator else 0.0:.6f}",
        }
        for (strand, tail_len, tail_sequence, tail_base), count in sorted(
            tail_summary_strand_counter.items(),
            key=lambda item: (item[0][0], -item[1], item[0][1]),
        )
    ]
    five_prime_summary_counter: Counter[tuple[int, str, str]] = Counter()
    for row in five_prime_rows:
        five_prime_summary_counter[
            (int(row["five_prime_length"]), str(row["five_prime_sequence"]), str(row["five_prime_base"]))
        ] += float(row["alignment_weight"])
    five_prime_summary_rows = [
        {
            "sample": sample,
            "five_prime_length": addition_len,
            "five_prime_sequence": sequence,
            "five_prime_base": base,
            "five_prime_count_weighted": f"{count:.6f}",
            "total_sRNAs_denominator": f"{total_srna_denominator:.6f}",
            "percent_of_total_sRNAs": f"{(count * 100.0 / total_srna_denominator) if total_srna_denominator else 0.0:.9f}",
            "mismatched_sRNAs_denominator": f"{mismatched_srna_denominator:.6f}",
            "percent_of_mismatched_sRNAs": f"{(count * 100.0 / mismatched_srna_denominator) if mismatched_srna_denominator else 0.0:.6f}",
        }
        for (addition_len, sequence, base), count in sorted(
            five_prime_summary_counter.items(),
            key=lambda item: (-item[1], item[0]),
        )
    ]
    five_prime_summary_strand_counter: Counter[tuple[str, int, str, str]] = Counter()
    for row in five_prime_rows:
        five_prime_summary_strand_counter[
            (str(row["strand"]), int(row["five_prime_length"]), str(row["five_prime_sequence"]), str(row["five_prime_base"]))
        ] += float(row["alignment_weight"])
    five_prime_summary_strand_rows = [
        {
            "sample": sample,
            "strand": strand,
            "five_prime_length": addition_len,
            "five_prime_sequence": sequence,
            "five_prime_base": base,
            "five_prime_count_weighted": f"{count:.6f}",
            "total_sRNAs_denominator": f"{total_srna_denominator:.6f}",
            "percent_of_total_sRNAs": f"{(count * 100.0 / total_srna_denominator) if total_srna_denominator else 0.0:.9f}",
            "mismatched_sRNAs_denominator": f"{mismatched_srna_denominator:.6f}",
            "percent_of_mismatched_sRNAs": f"{(count * 100.0 / mismatched_srna_denominator) if mismatched_srna_denominator else 0.0:.6f}",
        }
        for (strand, addition_len, sequence, base), count in sorted(
            five_prime_summary_strand_counter.items(),
            key=lambda item: (item[0][0], -item[1], item[0][1]),
        )
    ]

    return {
        "read_stats": read_stat_rows,
        "tail_candidates": tail_rows,
        "tail_summary": tail_summary_rows,
        "tail_summary_strand": tail_summary_strand_rows,
        "five_prime_candidates": five_prime_rows,
        "five_prime_summary": five_prime_summary_rows,
        "five_prime_summary_strand": five_prime_summary_strand_rows,
        "mismatch_summary_strand": mismatch_summary_strand_rows(
            sample,
            mismatch_strand_counts,
            mismatch_strand_denominators,
            args.max_mismatches,
        ),
        "substitution_summary_all": substitution_summary_rows(sample, contigs, all_exposure, all_counts),
        "substitution_summary_non_tail": substitution_summary_rows(sample, contigs, nontail_exposure, nontail_counts),
        "substitution_summary_molecule_all": substitution_summary_rows(sample, contigs, molecule_all_exposure, molecule_all_counts),
        "substitution_summary_molecule_non_tail": substitution_summary_rows(
            sample,
            contigs,
            molecule_nontail_exposure,
            molecule_nontail_counts,
        ),
        "substitution_summary_strand_all": substitution_summary_strand_rows(
            sample,
            contigs,
            all_strand_pos_cov,
            all_strand_pos_subs,
        ),
        "substitution_summary_strand_non_tail": substitution_summary_strand_rows(
            sample,
            contigs,
            nontail_strand_pos_cov,
            nontail_strand_pos_subs,
        ),
        "substitution_summary_strand_molecule_all": substitution_summary_strand_rows(
            sample,
            contigs,
            molecule_all_strand_pos_cov,
            molecule_all_strand_pos_subs,
            orientation="molecule",
        ),
        "substitution_summary_strand_molecule_non_tail": substitution_summary_strand_rows(
            sample,
            contigs,
            molecule_nontail_strand_pos_cov,
            molecule_nontail_strand_pos_subs,
            orientation="molecule",
        ),
        "substitution_positions_all": substitution_position_rows(sample, contigs, all_pos_cov, all_pos_subs),
        "substitution_positions_non_tail": substitution_position_rows(sample, contigs, nontail_pos_cov, nontail_pos_subs),
        "substitution_positions_strand_all": substitution_position_strand_rows(sample, contigs, all_strand_pos_cov, all_strand_pos_subs),
        "substitution_positions_strand_non_tail": substitution_position_strand_rows(sample, contigs, nontail_strand_pos_cov, nontail_strand_pos_subs),
        "substitution_positions_strand_molecule_all": substitution_position_strand_rows(
            sample,
            contigs,
            molecule_all_strand_pos_cov,
            molecule_all_strand_pos_subs,
            orientation="molecule",
        ),
        "substitution_positions_strand_molecule_non_tail": substitution_position_strand_rows(
            sample,
            contigs,
            molecule_nontail_strand_pos_cov,
            molecule_nontail_strand_pos_subs,
            orientation="molecule",
        ),
    }


def sample_to_sam(sample: str, sam_dir: Path, mismatches: int) -> Path:
    direct = sam_dir / f"{sample}.v{mismatches}.sam"
    if direct.exists():
        return direct
    matches = sorted(sam_dir.glob(f"{sample}*.sam"))
    if not matches:
        raise FileNotFoundError(f"No SAM found for sample {sample} in {sam_dir}")
    return matches[0]


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True, help="Template FASTA used for SAM mapping")
    parser.add_argument("--outdir", type=Path, required=True, help="Output directory for edit-analysis tables and plots")
    parser.add_argument("--sam-dir", type=Path, help="Directory containing sample SAM files")
    parser.add_argument("--samples", nargs="*", default=[], help="Sample names to load from --sam-dir")
    parser.add_argument("--sam", action="append", type=Path, default=[], help="Explicit SAM path; may be repeated")
    parser.add_argument("--mismatches", type=int, default=3, help="SAM filename mismatch suffix for --samples lookup")
    parser.add_argument("--min-mismatches", type=int, default=1, help="Minimum best mismatch count to analyze")
    parser.add_argument("--max-mismatches", type=int, default=3, help="Maximum best mismatch count to analyze")
    parser.add_argument("--tail-anchor-len", type=int, default=21, help="Required perfect 5-prime anchor length for tail candidates")
    parser.add_argument("--max-tail-len", type=int, default=3, help="Maximum candidate 3-prime tail length")
    parser.add_argument(
        "--short-core-tail-len",
        type=int,
        default=20,
        help="Also treat one terminal mismatch in a core_len+1 read as a 1 nt 3-prime tail candidate; set 0 to disable",
    )
    parser.add_argument("--min-len", type=int, default=18)
    parser.add_argument("--max-len", type=int, default=30)
    parser.add_argument("--max-n-fraction", type=float, default=0.10)
    parser.add_argument("--max-base-fraction", type=float, default=0.80)
    parser.add_argument("--min-entropy", type=float, default=1.25)
    parser.add_argument("--max-run-fraction", type=float, default=0.60)
    parser.add_argument("--max-run-bases", type=int, default=8)
    parser.add_argument("--max-dinuc-fraction", type=float, default=0.70)
    parser.add_argument("--max-trinuc-fraction", type=float, default=0.70)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "tables").mkdir(exist_ok=True)
    (args.outdir / "plots").mkdir(exist_ok=True)

    contig_list = srna.read_fasta(args.fasta)
    contigs = {contig.name: contig for contig in contig_list}
    sam_paths: list[tuple[str, Path]] = []
    if args.sam_dir and args.samples:
        sam_paths.extend((sample, sample_to_sam(sample, args.sam_dir, args.mismatches)) for sample in args.samples)
    for sam_path in args.sam:
        sample = sam_path.name.split(".")[0]
        sam_paths.append((sample, sam_path))
    if not sam_paths:
        raise SystemExit("Provide --sam-dir with --samples, or one or more --sam paths")

    combined: dict[str, list[dict]] = defaultdict(list)
    for sample, sam_path in sam_paths:
        print(f"Analyzing {sample}: {sam_path}", flush=True)
        result = analyze_sample(sample, sam_path, contigs, args)
        for key, rows in result.items():
            combined[key].extend(rows)
            write_csv(args.outdir / "tables" / f"{sample}.{key}.csv", rows)

    table_fields = {
        "read_stats": ["sample", "metric", "value"],
        "tail_candidates": [
            "sample",
            "read_id",
            "contig",
            "strand",
            "template_start_1based",
            "template_end_1based",
            "read_length",
            "best_nm",
            "alignment_weight",
            "tail_length",
            "tail_sequence",
            "tail_is_homopolymer",
            "tail_base",
        ],
        "five_prime_candidates": [
            "sample",
            "read_id",
            "contig",
            "strand",
            "template_start_1based",
            "template_end_1based",
            "read_length",
            "best_nm",
            "alignment_weight",
            "five_prime_length",
            "five_prime_sequence",
            "five_prime_is_homopolymer",
            "five_prime_base",
        ],
        "tail_summary": [
            "sample",
            "tail_length",
            "tail_sequence",
            "tail_base",
            "tail_count_weighted",
            "total_sRNAs_denominator",
            "percent_of_total_sRNAs",
            "mismatched_sRNAs_denominator",
            "percent_of_mismatched_sRNAs",
        ],
        "tail_summary_strand": [
            "sample",
            "strand",
            "tail_length",
            "tail_sequence",
            "tail_base",
            "tail_count_weighted",
            "total_sRNAs_denominator",
            "percent_of_total_sRNAs",
            "mismatched_sRNAs_denominator",
            "percent_of_mismatched_sRNAs",
        ],
        "tail_summary_grouped_strand": [
            "sample",
            "strand",
            "tail_group",
            "grouped_tail_members",
            "tail_count_weighted",
            "total_sRNAs_denominator",
            "percent_of_total_sRNAs",
            "mismatched_sRNAs_denominator",
            "percent_of_mismatched_sRNAs",
        ],
        "tail_summary_grouped_components_strand": [
            "sample",
            "strand",
            "tail_group",
            "tail_component",
            "tail_component_dna",
            "tail_count_weighted",
            "total_sRNAs_denominator",
            "percent_of_total_sRNAs",
            "mismatched_sRNAs_denominator",
            "percent_of_mismatched_sRNAs",
        ],
        "five_prime_summary": [
            "sample",
            "five_prime_length",
            "five_prime_sequence",
            "five_prime_base",
            "five_prime_count_weighted",
            "total_sRNAs_denominator",
            "percent_of_total_sRNAs",
            "mismatched_sRNAs_denominator",
            "percent_of_mismatched_sRNAs",
        ],
        "five_prime_summary_strand": [
            "sample",
            "strand",
            "five_prime_length",
            "five_prime_sequence",
            "five_prime_base",
            "five_prime_count_weighted",
            "total_sRNAs_denominator",
            "percent_of_total_sRNAs",
            "mismatched_sRNAs_denominator",
            "percent_of_mismatched_sRNAs",
        ],
        "five_prime_summary_grouped_components_strand": [
            "sample",
            "strand",
            "tail_group",
            "tail_component",
            "tail_component_dna",
            "tail_count_weighted",
            "total_sRNAs_denominator",
            "percent_of_total_sRNAs",
            "mismatched_sRNAs_denominator",
            "percent_of_mismatched_sRNAs",
        ],
        "mismatch_summary_strand": [
            "sample",
            "strand",
            "best_mismatches",
            "reads_weighted",
            "strand_mapped_reads_weighted",
            "percent_within_strand_mapped",
        ],
        "substitution_summary_all": [
            "sample",
            "contig",
            "substitution",
            "substitution_rna",
            "template_base",
            "read_base",
            "template_base_appearances_weighted",
            "substitution_events_weighted",
            "percent_of_template_base_appearance",
        ],
        "substitution_summary_non_tail": [
            "sample",
            "contig",
            "substitution",
            "substitution_rna",
            "template_base",
            "read_base",
            "template_base_appearances_weighted",
            "substitution_events_weighted",
            "percent_of_template_base_appearance",
        ],
        "substitution_summary_molecule_all": [
            "sample",
            "contig",
            "substitution",
            "template_base",
            "read_base",
            "template_base_appearances_weighted",
            "substitution_events_weighted",
            "percent_of_template_base_appearance",
        ],
        "substitution_summary_molecule_non_tail": [
            "sample",
            "contig",
            "substitution",
            "template_base",
            "read_base",
            "template_base_appearances_weighted",
            "substitution_events_weighted",
            "percent_of_template_base_appearance",
        ],
        "substitution_summary_molecule_all_percent_of_mismatched": [
            "sample",
            "contig",
            "substitution",
            "substitution_rna",
            "template_base",
            "read_base",
            "template_base_appearances_weighted",
            "substitution_events_weighted",
            "percent_of_template_base_appearance",
            "mismatched_sRNAs_denominator",
            "percent_of_mismatched_sRNAs",
        ],
        "substitution_summary_molecule_non_tail_percent_of_mismatched": [
            "sample",
            "contig",
            "substitution",
            "substitution_rna",
            "template_base",
            "read_base",
            "template_base_appearances_weighted",
            "substitution_events_weighted",
            "percent_of_template_base_appearance",
            "mismatched_sRNAs_denominator",
            "percent_of_mismatched_sRNAs",
        ],
        "substitution_summary_strand_all": [
            "sample",
            "contig",
            "strand",
            "substitution",
            "reference_base",
            "read_base",
            "reference_base_appearances_weighted",
            "substitution_events_weighted",
            "percent_of_reference_base_appearance",
        ],
        "substitution_summary_strand_non_tail": [
            "sample",
            "contig",
            "strand",
            "substitution",
            "reference_base",
            "read_base",
            "reference_base_appearances_weighted",
            "substitution_events_weighted",
            "percent_of_reference_base_appearance",
        ],
        "substitution_summary_strand_molecule_all": [
            "sample",
            "contig",
            "strand",
            "substitution",
            "reference_base",
            "read_base",
            "reference_base_appearances_weighted",
            "substitution_events_weighted",
            "percent_of_reference_base_appearance",
        ],
        "substitution_summary_strand_molecule_non_tail": [
            "sample",
            "contig",
            "strand",
            "substitution",
            "reference_base",
            "read_base",
            "reference_base_appearances_weighted",
            "substitution_events_weighted",
            "percent_of_reference_base_appearance",
        ],
        "substitution_summary_strand_molecule_all_percent_of_mismatched": [
            "sample",
            "contig",
            "strand",
            "substitution",
            "substitution_rna",
            "reference_base",
            "read_base",
            "reference_base_appearances_weighted",
            "substitution_events_weighted",
            "percent_of_reference_base_appearance",
            "mismatched_sRNAs_denominator",
            "percent_of_mismatched_sRNAs",
        ],
        "substitution_summary_strand_molecule_non_tail_percent_of_mismatched": [
            "sample",
            "contig",
            "strand",
            "substitution",
            "substitution_rna",
            "reference_base",
            "read_base",
            "reference_base_appearances_weighted",
            "substitution_events_weighted",
            "percent_of_reference_base_appearance",
            "mismatched_sRNAs_denominator",
            "percent_of_mismatched_sRNAs",
        ],
        "substitution_summary_strand_molecule_non_tail_collapsed_rc_families": [
            "sample",
            "strand",
            "mutation_family",
            "substitution_events_weighted",
            "reference_base_appearances_weighted",
            "percent_rate",
        ],
        "substitution_positions_all": [
            "sample",
            "contig",
            "position_1based",
            "template_base",
            "substitution",
            "coverage_weighted",
            "substitution_events_weighted",
            "percent_at_position",
        ],
        "substitution_positions_non_tail": [
            "sample",
            "contig",
            "position_1based",
            "template_base",
            "substitution",
            "coverage_weighted",
            "substitution_events_weighted",
            "percent_at_position",
        ],
        "substitution_positions_strand_all": [
            "sample",
            "contig",
            "position_1based",
            "template_base",
            "sense_reference_base",
            "antisense_reference_base",
            "substitution",
            "sense_coverage_weighted",
            "sense_substitution_events_weighted",
            "sense_percent_at_position",
            "antisense_coverage_weighted",
            "antisense_substitution_events_weighted",
            "antisense_percent_at_position",
        ],
        "substitution_positions_strand_non_tail": [
            "sample",
            "contig",
            "position_1based",
            "template_base",
            "sense_reference_base",
            "antisense_reference_base",
            "substitution",
            "sense_coverage_weighted",
            "sense_substitution_events_weighted",
            "sense_percent_at_position",
            "antisense_coverage_weighted",
            "antisense_substitution_events_weighted",
            "antisense_percent_at_position",
        ],
        "substitution_positions_strand_molecule_all": [
            "sample",
            "contig",
            "position_1based",
            "template_base",
            "sense_reference_base",
            "antisense_reference_base",
            "substitution",
            "sense_coverage_weighted",
            "sense_substitution_events_weighted",
            "sense_percent_at_position",
            "antisense_coverage_weighted",
            "antisense_substitution_events_weighted",
            "antisense_percent_at_position",
        ],
        "substitution_positions_strand_molecule_non_tail": [
            "sample",
            "contig",
            "position_1based",
            "template_base",
            "sense_reference_base",
            "antisense_reference_base",
            "substitution",
            "sense_coverage_weighted",
            "sense_substitution_events_weighted",
            "sense_percent_at_position",
            "antisense_coverage_weighted",
            "antisense_substitution_events_weighted",
            "antisense_percent_at_position",
        ],
    }
    for key, rows in combined.items():
        write_csv(args.outdir / "tables" / f"all_samples.{key}.csv", rows, table_fields.get(key))

    substitution_percent_keys = {
        "substitution_summary_molecule_all_percent_of_mismatched": "substitution_summary_molecule_all",
        "substitution_summary_molecule_non_tail_percent_of_mismatched": "substitution_summary_molecule_non_tail",
        "substitution_summary_strand_molecule_all_percent_of_mismatched": "substitution_summary_strand_molecule_all",
        "substitution_summary_strand_molecule_non_tail_percent_of_mismatched": "substitution_summary_strand_molecule_non_tail",
    }
    for out_key, source_key in substitution_percent_keys.items():
        if combined[source_key] and combined["read_stats"]:
            rows = substitution_summary_percent_of_mismatched_rows(combined[source_key], combined["read_stats"])
            combined[out_key] = rows
            write_csv(args.outdir / "tables" / f"all_samples.{out_key}.csv", rows, table_fields[out_key])

    if combined["tail_summary_strand"]:
        grouped_tail_rows = grouped_tail_summary_rows(combined["tail_summary_strand"])
        combined["tail_summary_grouped_strand"] = grouped_tail_rows
        write_csv(
            args.outdir / "tables" / "all_samples.tail_summary_grouped_strand.csv",
            grouped_tail_rows,
            table_fields["tail_summary_grouped_strand"],
        )
        grouped_tail_component_rows = grouped_terminal_component_summary_rows(
            combined["tail_summary_strand"],
            sequence_key="tail_sequence",
            count_key="tail_count_weighted",
        )
        combined["tail_summary_grouped_components_strand"] = grouped_tail_component_rows
        write_csv(
            args.outdir / "tables" / "all_samples.tail_summary_grouped_components_strand.csv",
            grouped_tail_component_rows,
            table_fields["tail_summary_grouped_components_strand"],
        )
    if combined["five_prime_summary_strand"]:
        grouped_five_prime_component_rows = grouped_terminal_component_summary_rows(
            combined["five_prime_summary_strand"],
            sequence_key="five_prime_sequence",
            count_key="five_prime_count_weighted",
        )
        combined["five_prime_summary_grouped_components_strand"] = grouped_five_prime_component_rows
        write_csv(
            args.outdir / "tables" / "all_samples.five_prime_summary_grouped_components_strand.csv",
            grouped_five_prime_component_rows,
            table_fields["five_prime_summary_grouped_components_strand"],
        )

    strand_summary_keys = [
        "substitution_summary_strand_all",
        "substitution_summary_strand_non_tail",
        "substitution_summary_strand_molecule_all",
        "substitution_summary_strand_molecule_non_tail",
    ]
    for key in strand_summary_keys:
        if combined[key]:
            overall_rows = aggregate_substitution_summary_strand_rows(combined[key])
            write_csv(args.outdir / "tables" / f"overall.{key}.csv", overall_rows, table_fields.get(key))

    if args.min_mismatches == args.max_mismatches:
        mismatch_phrase = f"sRNAs with exactly {args.min_mismatches} mismatch"
    else:
        mismatch_phrase = f"sRNAs with {args.min_mismatches}-{args.max_mismatches} mismatches"
    mismatch_ylabel = f"% of {mismatch_phrase}"

    if combined["substitution_summary_all"]:
        plot_substitution_bar(
            combined["substitution_summary_all"],
            args.outdir / "plots" / "all_samples.substitution_summary_all.png",
            "All mismatched reads: substitution spectrum",
        )
        plot_substitution_bar(
            combined["substitution_summary_non_tail"],
            args.outdir / "plots" / "all_samples.substitution_summary_non_tail.png",
            "Non-tail-candidate reads: substitution spectrum",
        )
        plot_substitution_bar(
            combined["substitution_summary_molecule_all"],
            args.outdir / "plots" / "all_samples.substitution_summary_molecule_oriented_all.png",
            "All mismatched reads: molecule-oriented substitution spectrum",
        )
        plot_substitution_bar(
            combined["substitution_summary_molecule_non_tail"],
            args.outdir / "plots" / "all_samples.substitution_summary_molecule_oriented_non_tail.png",
            "Non-tail-candidate reads: molecule-oriented substitution spectrum",
        )
        plot_substitution_template_appearance_sample_bar(
            combined["substitution_summary_molecule_all"],
            args.outdir / "plots" / "all_samples.substitution_summary_molecule_oriented_template_appearance_sample_mean_all.png",
            "Molecule-oriented substitutions: % of template nucleotide appearances",
        )
        plot_substitution_template_appearance_sample_bar(
            combined["substitution_summary_molecule_non_tail"],
            args.outdir / "plots" / "all_samples.substitution_summary_molecule_oriented_template_appearance_sample_mean_non_tail.png",
            "Molecule-oriented substitutions, non-tail candidates: % of template nucleotide appearances",
        )
        if combined["substitution_summary_molecule_all_percent_of_mismatched"]:
            plot_substitution_percent_of_mismatched(
                combined["substitution_summary_molecule_all_percent_of_mismatched"],
                args.outdir / "plots" / "all_samples.substitution_summary_molecule_oriented_percent_of_mismatched_sRNAs_all.png",
                f"Molecule-oriented substitutions: {mismatch_ylabel}",
                ylabel=mismatch_ylabel,
            )
        if combined["substitution_summary_molecule_non_tail_percent_of_mismatched"]:
            plot_substitution_percent_of_mismatched(
                combined["substitution_summary_molecule_non_tail_percent_of_mismatched"],
                args.outdir / "plots" / "all_samples.substitution_summary_molecule_oriented_percent_of_mismatched_sRNAs_non_tail.png",
                f"Molecule-oriented substitutions, non-tail candidates: {mismatch_ylabel}",
                ylabel=mismatch_ylabel,
            )
        if combined["substitution_summary_strand_molecule_all_percent_of_mismatched"]:
            plot_substitution_percent_of_mismatched_by_strand(
                combined["substitution_summary_strand_molecule_all_percent_of_mismatched"],
                args.outdir / "plots" / "all_samples.substitution_summary_by_strand_molecule_oriented_percent_of_mismatched_sRNAs_all.png",
                f"Molecule-oriented substitutions by strand: {mismatch_ylabel}",
                ylabel=mismatch_ylabel,
            )
        if combined["substitution_summary_strand_molecule_non_tail_percent_of_mismatched"]:
            plot_substitution_percent_of_mismatched_by_strand(
                combined["substitution_summary_strand_molecule_non_tail_percent_of_mismatched"],
                args.outdir / "plots" / "all_samples.substitution_summary_by_strand_molecule_oriented_percent_of_mismatched_sRNAs_non_tail.png",
                f"Molecule-oriented substitutions by strand, non-tail candidates: {mismatch_ylabel}",
                ylabel=mismatch_ylabel,
            )
        plot_substitution_strand_bar(
            combined["substitution_summary_strand_all"],
            args.outdir / "plots" / "all_samples.substitution_summary_by_strand_all.png",
            "All mismatched reads: substitution rates by strand",
        )
        plot_substitution_strand_bar(
            combined["substitution_summary_strand_non_tail"],
            args.outdir / "plots" / "all_samples.substitution_summary_by_strand_non_tail.png",
            "Non-tail-candidate reads: substitution rates by strand",
        )
        plot_substitution_strand_bar(
            combined["substitution_summary_strand_molecule_all"],
            args.outdir / "plots" / "all_samples.substitution_summary_by_strand_molecule_oriented_all.png",
            "All mismatched reads: molecule-oriented substitution rates by strand",
        )
        plot_substitution_strand_bar(
            combined["substitution_summary_strand_molecule_non_tail"],
            args.outdir / "plots" / "all_samples.substitution_summary_by_strand_molecule_oriented_non_tail.png",
            "Non-tail-candidate reads: molecule-oriented substitution rates by strand",
        )
        collapsed_family_rows = collapsed_rc_mutation_family_rows(combined["substitution_summary_strand_molecule_non_tail"])
        combined["substitution_summary_strand_molecule_non_tail_collapsed_rc_families"] = collapsed_family_rows
        write_csv(
            args.outdir / "tables" / "all_samples.substitution_summary_strand_molecule_non_tail_collapsed_rc_families.csv",
            collapsed_family_rows,
            table_fields["substitution_summary_strand_molecule_non_tail_collapsed_rc_families"],
        )
        plot_collapsed_rc_mutation_family_strand_bar(
            collapsed_family_rows,
            args.outdir / "plots" / "all_samples.substitution_summary_by_strand_molecule_oriented_non_tail_collapsed_rc_families.png",
            "Non-tail reads: collapsed reverse-complement mutation families",
        )
    for contig in contig_list:
        contig_rows_all = [row for row in combined["substitution_positions_all"] if row["contig"] == contig.name]
        contig_rows_non_tail = [row for row in combined["substitution_positions_non_tail"] if row["contig"] == contig.name]
        contig_strand_rows_all = [row for row in combined["substitution_positions_strand_all"] if row["contig"] == contig.name]
        contig_strand_rows_non_tail = [row for row in combined["substitution_positions_strand_non_tail"] if row["contig"] == contig.name]
        contig_molecule_strand_rows_all = [
            row for row in combined["substitution_positions_strand_molecule_all"] if row["contig"] == contig.name
        ]
        contig_molecule_strand_rows_non_tail = [
            row for row in combined["substitution_positions_strand_molecule_non_tail"] if row["contig"] == contig.name
        ]
        if contig_rows_all:
            plot_position_heatmap(
                contig_rows_all,
                contig,
                args.outdir / "plots" / f"{contig.name}.substitution_position_heatmap_all.png",
                f"{contig.name}: substitution percentage by position, all mismatched reads",
            )
            plot_position_heatmap(
                contig_rows_non_tail,
                contig,
                args.outdir / "plots" / f"{contig.name}.substitution_position_heatmap_non_tail.png",
                f"{contig.name}: substitution percentage by position, non-tail candidates",
            )
            plot_position_lollipop_grid(
                contig_strand_rows_all,
                contig,
                args.outdir / "plots" / f"{contig.name}.substitution_position_lollipop_all.png",
                f"{contig.name}: all positional substitution rates",
            )
            plot_position_lollipop_grid(
                contig_strand_rows_non_tail,
                contig,
                args.outdir / "plots" / f"{contig.name}.substitution_position_lollipop_non_tail.png",
                f"{contig.name}: positional substitution rates, non-tail candidates",
            )
            plot_position_lollipop_grid(
                contig_molecule_strand_rows_all,
                contig,
                args.outdir / "plots" / f"{contig.name}.substitution_position_lollipop_molecule_oriented_all.png",
                f"{contig.name}: molecule-oriented positional substitution rates",
            )
            plot_position_lollipop_grid(
                contig_molecule_strand_rows_non_tail,
                contig,
                args.outdir / "plots" / f"{contig.name}.substitution_position_lollipop_molecule_oriented_non_tail.png",
                f"{contig.name}: molecule-oriented positional substitution rates, non-tail candidates",
            )
            for sample in sorted({str(row["sample"]) for row in contig_strand_rows_non_tail}):
                sample_rows = [row for row in contig_strand_rows_non_tail if row["sample"] == sample]
                plot_position_lollipop_grid(
                    sample_rows,
                    contig,
                    args.outdir / "plots" / f"{sample}.{contig.name}.substitution_position_lollipop_non_tail.png",
                    f"{sample} {contig.name}: positional substitution rates, non-tail candidates",
                )
            for sample in sorted({str(row["sample"]) for row in contig_molecule_strand_rows_non_tail}):
                sample_rows = [row for row in contig_molecule_strand_rows_non_tail if row["sample"] == sample]
                plot_position_lollipop_grid(
                    sample_rows,
                    contig,
                    args.outdir / "plots" / f"{sample}.{contig.name}.substitution_position_lollipop_molecule_oriented_non_tail.png",
                    f"{sample} {contig.name}: molecule-oriented positional rates, non-tail candidates",
                )
    if combined["tail_summary"]:
        plot_tailing_summary(combined["tail_summary"], args.outdir / "plots" / "all_samples.tailing_summary.png")
    if combined["mismatch_summary_strand"]:
        plot_mismatch_summary_by_strand(
            combined["mismatch_summary_strand"],
            args.outdir / "plots" / "all_samples.best_mismatch_ratios_by_strand.png",
        )
    if combined["tail_summary_strand"]:
        plot_tailing_summary_by_strand(
            combined["tail_summary_strand"],
            args.outdir / "plots" / "all_samples.tailing_summary_by_strand.png",
        )
        plot_tailing_summary_by_strand(
            combined["tail_summary_strand"],
            args.outdir / "plots" / "all_samples.tailing_summary_by_strand_percent_of_mismatched_sRNAs.png",
            value_column="percent_of_mismatched_sRNAs",
            ylabel="% of sRNAs with >=1 mismatch",
            title="Candidate 3-prime tailing among mismatched sRNAs by strand",
        )
    if combined["tail_summary_grouped_strand"]:
        plot_grouped_tail_summary_by_strand(
            combined["tail_summary_grouped_strand"],
            args.outdir / "plots" / "all_samples.tailing_summary_grouped_by_strand_percent_of_mismatched_sRNAs.png",
        )
    if combined["tail_summary_grouped_components_strand"] and combined["five_prime_summary_grouped_components_strand"]:
        plot_terminal_addition_comparison(
            combined["tail_summary_grouped_components_strand"],
            combined["five_prime_summary_grouped_components_strand"],
            args.outdir / "plots" / "all_samples.terminal_additions_3prime_vs_5prime_grouped_components_percent_of_mismatched_sRNAs.png",
        )
    if combined["five_prime_summary_strand"]:
        plot_five_prime_summary_by_strand(
            combined["five_prime_summary_strand"],
            args.outdir / "plots" / "all_samples.five_prime_summary_by_strand.png",
        )

    print(args.outdir / "tables")
    print(args.outdir / "plots")


if __name__ == "__main__":
    main()
