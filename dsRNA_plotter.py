#!/usr/bin/env python3
"""Grouped dsRNA directional coverage plotter for INCI.

This tool accepts trimmed paired-end RNA-seq samples organized as biological
replicates, maps them to template sequences, collapses each read pair to one
fragment footprint, assigns direction from read 1 by default, normalizes by
global trimmed read pairs, and writes one bidirectional coverage plot per contig.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib-inci-dsrna"))

import matplotlib.pyplot as plt
import numpy as np

from dsRNA_identification import SampleInput, coverage_array, load_fasta, map_sample, sanitize_name
from plot_style import apply_matplotlib_style
from sirna_annotations import draw_sirna_annotation_track, find_sirna_annotations, write_sirna_annotation_table


DEFAULT_COLORS = [
    "#006400",
    "#8B006B",
    "#005F73",
    "#9B2226",
    "#6A4C93",
    "#2F4858",
]


@dataclass(frozen=True)
class GroupedSample:
    name: str
    group: str
    replicate: str
    r1: Path
    r2: Path
    library_pairs: int | None = None


@dataclass(frozen=True)
class Feature:
    contig: str
    start: int
    end: int
    label: str
    kind: str = "span"


def moving_average(values: np.ndarray, span: int) -> np.ndarray:
    if span <= 1:
        return values
    kernel = np.ones(span, dtype=np.float64) / span
    return np.convolve(values, kernel, mode="same")


def mean_and_sem(rows: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    matrix = np.vstack(rows)
    mean = np.mean(matrix, axis=0)
    if matrix.shape[0] == 1:
        return mean, np.zeros(matrix.shape[1], dtype=np.float64)
    sem = np.std(matrix, axis=0, ddof=1) / math.sqrt(matrix.shape[0])
    return mean, sem


def parse_sample_arg(item: list[str]) -> GroupedSample:
    if len(item) != 4:
        raise ValueError("--sample must be exactly GROUP REPLICATE R1 R2")
    group, replicate, r1, r2 = item
    name = sanitize_name(f"{group}_rep_{replicate}")
    return GroupedSample(name, group, replicate, Path(r1), Path(r2))


def read_samples_csv(path: Path) -> list[GroupedSample]:
    samples: list[GroupedSample] = []
    group_counts: dict[str, int] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"Sample table has no header: {path}")
        for row in reader:
            group = (row.get("group") or row.get("condition") or row.get("sample_group") or "").strip()
            sample_id = (row.get("sample_id") or row.get("sample") or row.get("name") or "").strip()
            if not group:
                group = sample_id or "group"
            group_counts[group] = group_counts.get(group, 0) + 1
            replicate = (row.get("replicate") or row.get("rep") or str(group_counts[group])).strip()
            r1 = (row.get("trimmed_read1") or row.get("read1") or row.get("r1") or row.get("R1") or "").strip()
            r2 = (row.get("trimmed_read2") or row.get("read2") or row.get("r2") or row.get("R2") or "").strip()
            if not r1 or not r2:
                continue
            library_text = (row.get("library_pairs") or row.get("trimmed_pairs") or row.get("read_pairs") or "").strip()
            library_pairs = int(library_text.replace(",", "")) if library_text else None
            name = sanitize_name(sample_id or f"{group}_rep_{replicate}")
            samples.append(GroupedSample(name, group, replicate, Path(r1), Path(r2), library_pairs))
    if not samples:
        raise ValueError(f"No paired samples found in {path}")
    return samples


def load_samples(args: argparse.Namespace) -> list[GroupedSample]:
    samples: list[GroupedSample] = []
    for item in args.sample or []:
        samples.append(parse_sample_arg(item))
    if args.samples_csv:
        samples.extend(read_samples_csv(args.samples_csv))
    if not samples:
        raise ValueError("Provide --samples-csv or at least one --sample GROUP REPLICATE R1 R2 [LIBRARY_PAIRS].")
    seen: set[str] = set()
    for sample in samples:
        if sample.name in seen:
            raise ValueError(f"Duplicate sample name after sanitizing: {sample.name}")
        seen.add(sample.name)
        if not sample.r1.exists() or not sample.r2.exists():
            raise FileNotFoundError(f"Missing FASTQ for {sample.name}: {sample.r1}, {sample.r2}")
    return samples


def parse_key_value(items: list[str] | None) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for item in items or []:
        if "=" not in item:
            raise ValueError(f"Expected KEY=VALUE, got {item!r}")
        key, value = item.split("=", 1)
        parsed[key] = value
    return parsed


def parse_feature_text(text: str) -> Feature:
    parts = text.split(":")
    if len(parts) < 4:
        raise ValueError("--feature must be CONTIG:START:END:LABEL[:KIND]")
    contig, start, end, label = parts[:4]
    kind = parts[4] if len(parts) > 4 else "span"
    return Feature(contig, int(start), int(end), label, kind)


def read_features_csv(path: Path) -> list[Feature]:
    features: list[Feature] = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            contig = (row.get("contig") or row.get("sequence") or "*").strip()
            label = (row.get("label") or row.get("name") or "").strip()
            if not label:
                continue
            features.append(
                Feature(
                    contig=contig,
                    start=int(row.get("start", "0")),
                    end=int(row.get("end", "0")),
                    label=label,
                    kind=(row.get("kind") or "span").strip(),
                )
            )
    return features


def load_features(args: argparse.Namespace) -> list[Feature]:
    features = [parse_feature_text(item) for item in args.feature or []]
    if args.features_csv:
        features.extend(read_features_csv(args.features_csv))
    return features


def features_for_contig(features: list[Feature], contig: str) -> list[Feature]:
    return [feature for feature in features if feature.contig in {contig, "*", "all"}]


def safe_contig_filename(contig: str) -> str:
    return sanitize_name(contig)[:180]


def write_per_base_table(
    path: Path,
    contig: str,
    length: int,
    samples: list[GroupedSample],
    coverage_by_sample: dict[str, object],
    smooth_span: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = ["position"]
    sample_arrays: list[tuple[str, np.ndarray, np.ndarray]] = []
    for sample in samples:
        coverage = coverage_by_sample[sample.name]
        scale = 1_000_000.0 / coverage.library_pairs if coverage.library_pairs else 0.0
        sense = moving_average(coverage_array(coverage, contig, "sense", length) * scale, smooth_span)
        antisense = moving_average(coverage_array(coverage, contig, "antisense", length) * scale, smooth_span)
        clean_group = sample.group.replace("\t", " ").strip()
        columns.extend(
            [
                f"{clean_group} sense {sample.replicate}",
                f"{clean_group} antisense {sample.replicate}",
            ]
        )
        sample_arrays.append((sample.name, sense, antisense))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(columns)
        for idx in range(length):
            row: list[object] = [idx]
            for _sample_name, sense, antisense in sample_arrays:
                row.extend([f"{sense[idx]:.6f}", f"{-antisense[idx]:.6f}"])
            writer.writerow(row)


def plot_contig(
    contig: str,
    length: int,
    samples: list[GroupedSample],
    coverage_by_sample: dict[str, object],
    out_png: Path,
    features: list[Feature],
    group_colors: dict[str, str],
    group_labels: dict[str, str],
    smooth_span: int,
    hide_axis_titles: bool,
    xlim: tuple[int, int] | None,
    sirna_annotations: list[dict[str, object]],
) -> dict[str, float | int | str]:
    style = apply_matplotlib_style()
    x = np.arange(length)
    fig, ax = plt.subplots(figsize=(style["figure_width"], style["figure_height"]), constrained_layout=True)

    groups = list(dict.fromkeys(sample.group for sample in samples))
    maxima: list[float] = []
    summary: dict[str, float | int | str] = {"contig": contig, "length": length}

    for group_index, group in enumerate(groups):
        color = group_colors.get(group, DEFAULT_COLORS[group_index % len(DEFAULT_COLORS)])
        label = group_labels.get(group, group)
        group_samples = [sample for sample in samples if sample.group == group]
        sense_rows: list[np.ndarray] = []
        antisense_rows: list[np.ndarray] = []
        for sample in group_samples:
            coverage = coverage_by_sample[sample.name]
            scale = 1_000_000.0 / coverage.library_pairs if coverage.library_pairs else 0.0
            sense_rows.append(moving_average(coverage_array(coverage, contig, "sense", length) * scale, smooth_span))
            antisense_rows.append(moving_average(coverage_array(coverage, contig, "antisense", length) * scale, smooth_span))

        sense_mean, sense_sem = mean_and_sem(sense_rows)
        antisense_mean, antisense_sem = mean_and_sem(antisense_rows)
        signed_antisense_mean = -antisense_mean
        ax.plot(x, sense_mean, color=color, linewidth=style["line_width"], label=label)
        ax.plot(x, signed_antisense_mean, color=color, linewidth=style["line_width"], linestyle="--")
        sem_width = max(0.4, float(style["line_width"]) * 0.43)
        ax.plot(x, sense_mean + sense_sem, color=color, linewidth=sem_width, linestyle=":", alpha=0.9)
        ax.plot(x, sense_mean - sense_sem, color=color, linewidth=sem_width, linestyle=":", alpha=0.9)
        ax.plot(x, signed_antisense_mean + antisense_sem, color=color, linewidth=sem_width, linestyle=":", alpha=0.9)
        ax.plot(x, signed_antisense_mean - antisense_sem, color=color, linewidth=sem_width, linestyle=":", alpha=0.9)

        maxima.extend(
            [
                float(np.nanmax(np.abs(sense_mean + sense_sem))),
                float(np.nanmax(np.abs(sense_mean - sense_sem))),
                float(np.nanmax(np.abs(signed_antisense_mean + antisense_sem))),
                float(np.nanmax(np.abs(signed_antisense_mean - antisense_sem))),
            ]
        )
        duplex = np.minimum(sense_mean, antisense_mean)
        total = sense_mean + antisense_mean
        summary[f"{group}_mean_sense_CPM"] = float(np.mean(sense_mean))
        summary[f"{group}_mean_antisense_CPM"] = float(np.mean(antisense_mean))
        summary[f"{group}_mean_duplex_min_CPM"] = float(np.mean(duplex))
        summary[f"{group}_bidirectional_balance_percent"] = float(100.0 * (2.0 * np.sum(duplex) / np.sum(total))) if np.sum(total) else 0.0

    ax.axhline(0, color="0.10", linewidth=max(0.5, float(style["line_width"]) * 0.75))
    ymax = max(maxima) if maxima else 1.0
    ax.set_ylim(-ymax * 1.18, ymax * 1.18)
    ytop = ax.get_ylim()[1]

    for feature in features_for_contig(features, contig):
        start = max(0, feature.start)
        end = min(length, feature.end)
        if end <= start:
            continue
        kind = feature.kind.lower()
        if kind in {"bounds", "region-bounds", "boundary"}:
            ax.axvline(start, color="0.10", linewidth=style["line_width"], linestyle="--", alpha=0.95)
            ax.axvline(end, color="0.10", linewidth=style["line_width"], linestyle="--", alpha=0.95)
            ax.text(start + max(4, length * 0.003), ytop * 0.93, feature.label, color="0.10", fontsize=style["title_size"], fontweight="bold", ha="left", va="top")
        else:
            alpha = 0.07 if kind in {"span", "dsrna", "dsrna-region"} else 0.045
            ax.axvspan(start, end, color="0.35", alpha=alpha)
            if kind in {"gene", "annotation"}:
                yfactor = 0.92 if start < length * 0.25 else 0.64
            else:
                yfactor = 0.68 if feature.label.upper().startswith("A2") else 0.72
            fontstyle = "italic" if kind in {"gene", "annotation"} else "normal"
            ax.text(start + max(6, (end - start) * 0.08), ytop * yfactor, feature.label, color="0.20", fontsize=style["title_size"], fontstyle=fontstyle, fontweight="bold")

    draw_sirna_annotation_track(
        ax,
        sirna_annotations,
        contig,
        length,
        label_size=max(6.0, float(style["font_size"]) - 1.5),
        # This plot's historic x axis is zero-based, unlike the one-based sRNA mapper.
        position_offset=-1,
    )

    if xlim:
        ax.set_xlim(xlim[0], xlim[1])
    else:
        ax.set_xlim(-max(10, int(length * 0.03)), length)
    if hide_axis_titles:
        ax.set_xlabel("")
        ax.set_ylabel("")
    else:
        ax.set_xlabel("Position (nt)", fontsize=style["font_size"], fontweight="bold")
        ax.set_ylabel("Normalized depth (CPM; + sense / - antisense)", fontsize=style["font_size"], fontweight="bold")
    ax.grid(axis="x", color="0.88", linewidth=style["grid_width"])
    ax.tick_params(axis="both", labelsize=style["font_size"], width=max(0.6, float(style["line_width"]) * 0.6))
    for tick_label in ax.get_xticklabels() + ax.get_yticklabels():
        tick_label.set_fontweight("bold")
    ax.legend(ncols=1, prop={"size": max(6.0, float(style["font_size"]) - 1.0), "weight": "bold"}, frameon=False, loc="lower left")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=style["dpi"])
    plt.close(fig)
    return summary


def write_sample_summary(path: Path, samples: list[GroupedSample], coverage_by_sample: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for sample in samples:
        coverage = coverage_by_sample[sample.name]
        row = {
            "sample": sample.name,
            "group": sample.group,
            "replicate": sample.replicate,
            "library_pairs": coverage.library_pairs,
            "collapsed_pairs": coverage.stats.get("collapsed_pairs", 0),
            "sense_pairs": coverage.stats.get("sense_pairs", 0),
            "antisense_pairs": coverage.stats.get("antisense_pairs", 0),
            "discarded_pairs": coverage.stats.get("discarded_pairs", 0),
        }
        rows.append(row)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_contig_summary(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Map grouped paired-end RNA-seq replicates to template sequences and plot bidirectional dsRNA coverage per contig.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--reference-fasta", type=Path, required=True, help="Template FASTA/FNA containing one or more contigs.")
    parser.add_argument("--samples-csv", type=Path, help="CSV sample sheet. INCI trimming manifests are accepted.")
    parser.add_argument("--sample", nargs=4, action="append", metavar=("GROUP", "REPLICATE", "R1", "R2"), help="One sample row. Repeat as needed. Use --samples-csv to provide library_pairs.")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=7)
    parser.add_argument("--min-mapq", type=int, default=10)
    parser.add_argument("--direction-source", choices=("read1", "read2"), default="read1")
    parser.add_argument("--smooth-span", type=int, default=25)
    parser.add_argument("--group-color", action="append", help="Override group color as GROUP=#RRGGBB. Repeatable.")
    parser.add_argument("--group-label", action="append", help="Override legend label as GROUP=Label. Repeatable.")
    parser.add_argument("--feature", action="append", help="Annotation as CONTIG:START:END:LABEL[:KIND]. CONTIG may be all or *.")
    parser.add_argument("--features-csv", type=Path, help="CSV with contig,start,end,label,kind columns.")
    parser.add_argument("--sirna-annotations-fasta", type=Path, help="Optional siRNA FASTA to label as exact forward or reverse matches on coverage plots.")
    parser.add_argument("--xlim", nargs=2, type=int, metavar=("START", "END"), help="Optional x-axis limits for every contig plot.")
    parser.add_argument("--show-axis-titles", action="store_true", help="Show position and CPM axis titles.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sequences = load_fasta(args.reference_fasta)
    lengths = {name: len(seq) for name, seq in sequences.items()}
    samples = load_samples(args)
    features = load_features(args)
    sirna_annotations = (
        find_sirna_annotations(args.sirna_annotations_fasta, sequences)
        if args.sirna_annotations_fasta
        else []
    )
    group_colors = parse_key_value(args.group_color)
    group_labels = parse_key_value(args.group_label)

    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "plots").mkdir(exist_ok=True)
    (args.outdir / "tables").mkdir(exist_ok=True)
    if args.sirna_annotations_fasta:
        write_sirna_annotation_table(args.outdir / "tables" / "siRNA_annotation_matches.tsv", sirna_annotations)

    coverage_by_sample = {}
    for sample in samples:
        print(f"[dsRNA_plotter] mapping {sample.name}", flush=True)
        coverage_by_sample[sample.name] = map_sample(
            SampleInput(sample.name, sample.r1, sample.r2, sample.library_pairs),
            args.reference_fasta,
            lengths,
            args.outdir,
            threads=args.threads,
            min_mapq=args.min_mapq,
            direction_source=args.direction_source,
        )

    write_sample_summary(args.outdir / "tables" / "sample_summary.csv", samples, coverage_by_sample)

    summary_rows: list[dict[str, float | int | str]] = []
    xlim = tuple(args.xlim) if args.xlim else None
    for contig, length in lengths.items():
        slug = safe_contig_filename(contig)
        write_per_base_table(
            args.outdir / "tables" / f"{slug}.group_replicates.global_CPM.smoothed{args.smooth_span}.tsv",
            contig,
            length,
            samples,
            coverage_by_sample,
            args.smooth_span,
        )
        summary_rows.append(
            plot_contig(
                contig,
                length,
                samples,
                coverage_by_sample,
                args.outdir / "plots" / f"{slug}.dsRNA_directional_group_coverage.png",
                features,
                group_colors,
                group_labels,
                args.smooth_span,
                hide_axis_titles=not args.show_axis_titles,
                xlim=xlim,
                sirna_annotations=sirna_annotations,
            )
        )
    write_contig_summary(args.outdir / "tables" / "contig_summary.csv", summary_rows)
    manifest = {
        "reference_fasta": str(args.reference_fasta),
        "samples": [sample.__dict__ | {"r1": str(sample.r1), "r2": str(sample.r2)} for sample in samples],
        "contigs": lengths,
        "smooth_span": args.smooth_span,
        "direction_source": args.direction_source,
        "global_cpm": "per-base collapsed fragment coverage * 1,000,000 / trimmed read pairs in R1",
        "sirna_annotations_fasta": str(args.sirna_annotations_fasta) if args.sirna_annotations_fasta else "",
        "sirna_annotation_matches": len(sirna_annotations),
    }
    (args.outdir / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[dsRNA_plotter] done: {args.outdir}", flush=True)


if __name__ == "__main__":
    main()
