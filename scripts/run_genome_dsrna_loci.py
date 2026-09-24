#!/usr/bin/env python3
"""INCI first-phase genomic dsRNA locus identification.

The screen uses R1 direction, fractional valid multi-mapping, and global CPM
normalization.  It then maps all selected paired replicates only to the top
ranked genomic contexts for compact directional coverage plots.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent


def load_samples(path: Path) -> list[dict[str, str]]:
    samples: list[dict[str, str]] = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for index, row in enumerate(reader, start=1):
            name = (row.get("sample_id") or row.get("sample") or row.get("name") or f"sample_{index}").strip()
            group = (row.get("group") or row.get("condition") or row.get("sample_group") or name).strip()
            replicate = (row.get("replicate") or row.get("rep") or str(index)).strip()
            r1 = (row.get("trimmed_read1") or row.get("read1") or row.get("r1") or "").strip()
            r2 = (row.get("trimmed_read2") or row.get("read2") or row.get("r2") or "").strip()
            if not r1 or not r2:
                continue
            if not Path(r1).exists() or not Path(r2).exists():
                raise FileNotFoundError(f"Missing paired FASTQ for {name}: {r1}, {r2}")
            samples.append({"name": name, "group": group, "replicate": replicate, "r1": r1, "r2": r2})
    if not samples:
        raise ValueError("The selected sample manifest contains no valid paired-end samples.")
    return samples


def extract_top_contexts(reference_fasta: Path, ranking: pd.DataFrame, context_nt: int, output_fasta: Path, feature_csv: Path) -> None:
    by_contig: dict[str, list[pd.Series]] = defaultdict(list)
    for _, row in ranking.iterrows():
        by_contig[str(row["source_contig"])].append(row)
    records: list[tuple[str, str]] = []
    feature_rows: list[dict[str, object]] = []
    current: str | None = None
    chunks: list[str] = []

    def finish_record(contig: str | None, sequence_chunks: list[str]) -> None:
        if contig not in by_contig:
            return
        sequence = "".join(sequence_chunks).upper()
        for row in by_contig[contig]:
            focus_start0 = int(row["start_1based"]) - 1
            focus_end0 = int(row["end_1based"])
            context_start0 = max(0, focus_start0 - context_nt)
            context_end0 = min(len(sequence), focus_end0 + context_nt)
            target = str(row["target_bin"])
            records.append((target, sequence[context_start0:context_end0]))
            feature_rows.append(
                {
                    "contig": target,
                    "start": focus_start0 - context_start0,
                    "end": focus_end0 - context_start0,
                    "label": f"Rank {int(row['rank'])} locus",
                    "kind": "dsrna-region",
                }
            )

    for raw in reference_fasta.open(errors="replace"):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish_record(current, chunks)
            current = line[1:].split()[0]
            chunks = []
        elif current in by_contig:
            chunks.append(line)
    finish_record(current, chunks)
    if len(records) != len(ranking):
        found = {target for target, _ in records}
        missing = [str(row["target_bin"]) for _, row in ranking.iterrows() if str(row["target_bin"]) not in found]
        raise ValueError(f"Could not extract contexts for {', '.join(missing[:5])}")
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with output_fasta.open("w") as handle:
        for target, sequence in records:
            handle.write(f">{target}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")
    with feature_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["contig", "start", "end", "label", "kind"])
        writer.writeheader()
        writer.writerows(feature_rows)


def run(command: list[str], label: str) -> None:
    print(f"[{label}] {' '.join(command)}", flush=True)
    completed = subprocess.run(command, cwd=SCRIPT_DIR.parent, check=False)
    if completed.returncode:
        raise RuntimeError(f"{label} failed with exit code {completed.returncode}")


def main() -> None:
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--reference-fasta", type=Path, required=True, help="Genome or multi-contig reference FASTA/FNA.")
    parser.add_argument("--samples-csv", type=Path, required=True, help="INCI paired RNA-seq sample manifest.")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--scoring-group", default="", help="Biological group used for first-phase ranking. Defaults to the first selected group.")
    parser.add_argument("--bin-size", type=int, default=250)
    parser.add_argument("--plot-context-bp", type=int, default=500)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--summary-top-n", type=int, default=50)
    parser.add_argument("--max-reads", type=int, default=1_000_000, help="First-phase R1 reads per scoring sample; 0 uses the full library.")
    parser.add_argument("--threads", type=int, default=7)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--smooth-span", type=int, default=25)
    parser.add_argument("--sirna-annotations-fasta", type=Path, help="Optional siRNA FASTA to label as exact matches on top-hit coverage plots.")
    args = parser.parse_args()
    if args.bin_size < 25:
        raise ValueError("Bin size must be at least 25 bp.")
    if args.plot_context_bp < 0 or args.top_n < 1:
        raise ValueError("Plot context must be non-negative and top hits must be at least one.")
    if not args.reference_fasta.exists():
        raise FileNotFoundError(f"Reference FASTA does not exist: {args.reference_fasta}")
    if args.sirna_annotations_fasta and not args.sirna_annotations_fasta.exists():
        raise FileNotFoundError(f"siRNA annotation FASTA does not exist: {args.sirna_annotations_fasta}")

    samples = load_samples(args.samples_csv)
    group_order = list(dict.fromkeys(sample["group"] for sample in samples))
    scoring_group = args.scoring_group.strip() or group_order[0]
    scoring = [sample for sample in samples if sample["group"] == scoring_group]
    if not scoring:
        raise ValueError(f"Scoring group {scoring_group!r} is not in the selected sample set: {', '.join(group_order)}")

    phase_dir = args.outdir / "phase1"
    phase_command = [
        sys.executable,
        str(SCRIPT_DIR / "screen_genome_250bp_dsrna.py"),
        "--genome-fasta", str(args.reference_fasta),
        "--outdir", str(phase_dir),
        "--bin-size", str(args.bin_size),
        "--max-reads", str(args.max_reads),
        "--threads", str(args.threads),
        "--min-mapq", str(args.min_mapq),
        "--summary-top-n", str(args.summary_top_n),
        "--include-secondary",
        "--max-secondary", "1000",
        "--skip-standalone-node",
        "--phase1-only",
    ]
    for sample in scoring:
        phase_command.extend(["--screen-sample", f"{sample['name']}={sample['r1']}"])
    run(phase_command, "first-phase dsRNA screen")

    source_ranking = phase_dir / "whole_genome_250bp_ssRNase_ranked_candidates.csv"
    ranking = pd.read_csv(source_ranking, low_memory=False)
    ranking = ranking.rename(columns={"candidate_score_percent": "total_dsrna_score_percent"})
    ranking.insert(1, "scoring_group", scoring_group)
    ranking_path = args.outdir / "tables" / "dsrna_loci_all_scored_bins.csv"
    ranking_path.parent.mkdir(parents=True, exist_ok=True)
    ranking.to_csv(ranking_path, index=False)
    top = ranking.head(args.top_n).copy()
    top.to_csv(args.outdir / "tables" / "dsrna_loci_top_hits.csv", index=False)

    contexts = args.outdir / "top_hits" / "ranked_dsrna_locus_contexts.fasta"
    features = args.outdir / "top_hits" / "ranked_dsrna_locus_features.csv"
    extract_top_contexts(args.reference_fasta, top, args.plot_context_bp, contexts, features)
    plot_command = [
        sys.executable,
        str(SCRIPT_DIR.parent / "dsRNA_plotter.py"),
        "--reference-fasta", str(contexts),
        "--samples-csv", str(args.samples_csv),
        "--outdir", str(args.outdir / "top_hits"),
        "--threads", str(args.threads),
        "--min-mapq", str(args.min_mapq),
        "--smooth-span", str(args.smooth_span),
        "--features-csv", str(features),
        "--show-axis-titles",
    ]
    if args.sirna_annotations_fasta:
        plot_command.extend(["--sirna-annotations-fasta", str(args.sirna_annotations_fasta)])
    run(plot_command, "top-locus paired coverage plots")

    manifest = {
        "reference_fasta": str(args.reference_fasta),
        "samples_csv": str(args.samples_csv),
        "scoring_group": scoring_group,
        "scoring_samples": [sample["name"] for sample in scoring],
        "all_plot_groups": group_order,
        "bin_size": args.bin_size,
        "plot_context_bp": args.plot_context_bp,
        "top_n": args.top_n,
        "max_reads_per_scoring_sample": args.max_reads,
        "threads": args.threads,
        "multi_mapping": "All valid minimap2 primary/secondary R1 alignments; each read contributes a total fractional weight of 1 across placements.",
        "direction": "Read 1 orientation defines collapsed paired-fragment sense/antisense direction in top-hit plots.",
        "normalization": "Global CPM using trimmed R1 read pairs; phase-one denominator is max_reads when subsampling is enabled.",
        "sirna_annotations_fasta": str(args.sirna_annotations_fasta) if args.sirna_annotations_fasta else "",
    }
    (args.outdir / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[dsRNA-loci] done: {args.outdir}", flush=True)


if __name__ == "__main__":
    main()
