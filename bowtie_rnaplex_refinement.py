"""Refine Bowtie1 candidate pairs for original siRNAs with INCI RNAplex scoring."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

from degradome_analysis import DegradomeConfig, FastaRecord, write_rows
from dsrna_adaptation import (
    _target_rows,
    find_targets_for_candidates,
    make_direct_sirna_guides,
    read_fasta,
)


def candidates_from_bowtie_table(path: Path, transcriptome: str) -> dict[str, set[str]]:
    candidates: dict[str, set[str]] = defaultdict(set)
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["transcriptome"] == transcriptome and row["variant"] == "original":
                candidates[row["siRNA_id"]].add(row["transcript_id"])
    return candidates


def refine(
    sirna_fasta: Path, transcriptome_fasta: Path, bowtie_table: Path, bowtie_label: str,
    output_dir: Path, threads: int, cutoff: float,
) -> tuple[list[dict[str, object]], dict[str, int | float], list[dict[str, object]]]:
    records = read_fasta(sirna_fasta)
    guides = [guide for guide in make_direct_sirna_guides(records, 21, 1, 1) if guide.variant == "original"]
    by_id = {guide.query_id: guide for guide in guides}
    transcriptome = dict(read_fasta(transcriptome_fasta, allow_ambiguous=True))
    candidates = candidates_from_bowtie_table(bowtie_table, bowtie_label)
    query_records = [FastaRecord(guide.query_id, guide.sequence) for guide in guides]
    config = DegradomeConfig(
        srna_text="", srna_fasta=None, transcript_text="", transcript_fasta=None, samples=(), output_dir=output_dir,
        ignore_query_pos1=True, slice_positions=(10,), mfe_ratio_cutoff=cutoff, sort_by="mfe_ratio", threads=threads,
    )
    hits = find_targets_for_candidates(query_records, transcriptome, candidates, config, None)
    rows = _target_rows(hits, by_id, bowtie_label, {})
    hits_by_query: dict[str, list[object]] = defaultdict(list)
    for hit in hits:
        hits_by_query[str(hit.query)].append(hit)
    per_sirna: list[dict[str, object]] = []
    for guide in guides:
        guide_hits = hits_by_query.get(guide.query_id, [])
        ratios = [float(hit.mfe_ratio) for hit in guide_hits]
        site_energies = [float(hit.mfe_site) for hit in guide_hits]
        perfect_energies = [float(hit.mfe_perfect) for hit in guide_hits]
        per_sirna.append({
            "siRNA_id": guide.query_id, "dsrna_id": guide.dsrna_id, "dsrna_position": guide.dsrna_position,
            "strand": guide.strand, "siRNA_sequence": guide.sequence,
            "bowtie_candidate_transcript_count": len(candidates.get(guide.query_id, set())),
            "rnaplex_reported_hit_count": len(guide_hits),
            "mean_mfe_ratio": mean(ratios) if ratios else "", "median_mfe_ratio": median(ratios) if ratios else "",
            "maximum_mfe_ratio": max(ratios) if ratios else "", "mean_mfe_site": mean(site_energies) if site_energies else "",
            "mean_mfe_perfect": mean(perfect_energies) if perfect_energies else "",
        })
    stats = {
        "candidate_pair_count": sum(len(value) for value in candidates.values()),
        "candidate_siRNA_count": sum(bool(candidates.get(guide.query_id)) for guide in guides),
        "retained_mfe_target_hit_count": len(hits),
        "retained_mfe_target_transcript_count": len({str(hit.transcript) for hit in hits}),
        "mean_mfe_ratio": mean(float(hit.mfe_ratio) for hit in hits) if hits else "",
        "median_mfe_ratio": median(float(hit.mfe_ratio) for hit in hits) if hits else "",
        "mean_mfe_site": mean(float(hit.mfe_site) for hit in hits) if hits else "",
    }
    return rows, stats, per_sirna


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sirnas", type=Path, required=True)
    parser.add_argument("--focal-transcriptome", type=Path, required=True)
    parser.add_argument("--tribolium-transcriptome", type=Path, required=True)
    parser.add_argument("--drosophila-transcriptome", type=Path, required=True)
    parser.add_argument("--csfb-bowtie-table", type=Path, required=True)
    parser.add_argument("--tribolium-bowtie-table", type=Path, required=True)
    parser.add_argument("--drosophila-bowtie-table", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mfe-ratio-cutoff", type=float, default=0.70)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    sources = (
        ("CSFB", args.focal_transcriptome, args.csfb_bowtie_table, "focal"),
        ("Tribolium", args.tribolium_transcriptome, args.tribolium_bowtie_table, "control"),
        ("Drosophila", args.drosophila_transcriptome, args.drosophila_bowtie_table, "control"),
    )
    all_rows: list[dict[str, object]] = []
    summary: list[dict[str, object]] = []
    all_per_sirna: list[dict[str, object]] = []
    for name, transcriptome, table, label in sources:
        rows, stats, per_sirna = refine(args.sirnas, transcriptome, table, label, args.output_dir, args.threads, args.mfe_ratio_cutoff)
        for row in rows:
            row["reference_transcriptome"] = name
        all_rows.extend(rows)
        summary.append({"reference_transcriptome": name, **stats})
        for row in per_sirna:
            row["reference_transcriptome"] = name
        all_per_sirna.extend(per_sirna)
    write_rows(args.output_dir / "bowtie_rnaplex_refined_targets.tsv", all_rows)
    write_rows(args.output_dir / "bowtie_rnaplex_refinement_summary.tsv", summary)
    write_rows(args.output_dir / "bowtie_rnaplex_per_sirna_mfe_summary.tsv", all_per_sirna)
    print("RNAplex refinement completed.")


if __name__ == "__main__":
    main()
