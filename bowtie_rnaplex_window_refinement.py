"""RNAplex refinement confined to the exact Bowtie1 alignment windows."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from statistics import mean, median

from degradome_analysis import DegradomeConfig, FastaRecord, target_hits_for_query, write_rows
from dsrna_adaptation import make_direct_sirna_guides, read_fasta


def window_candidates(table: Path, label: str, transcripts: dict[str, str], flank: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with table.open(encoding="utf-8", newline="") as handle:
        for index, row in enumerate(csv.DictReader(handle, delimiter="\t"), start=1):
            if row["transcriptome"] != label or row["variant"] != "original":
                continue
            transcript = str(row["transcript_id"])
            sequence = transcripts.get(transcript)
            if not sequence:
                continue
            bowtie_start = int(row["target_start"])
            bowtie_stop = bowtie_start + len(str(row["aligned_sequence"])) - 1
            window_start = max(1, bowtie_start - flank)
            window_stop = min(len(sequence), bowtie_stop + flank)
            rows.append({
                "window_id": f"bowtie_window_{index}", "siRNA_id": row["siRNA_id"], "transcript_id": transcript,
                "bowtie_start": bowtie_start, "bowtie_stop": bowtie_stop, "window_start": window_start,
                "window_sequence": sequence[window_start - 1:window_stop],
            })
    return rows


def refine_reference(
    name: str, transcriptome_fasta: Path, bowtie_table: Path, bowtie_label: str, guides: dict[str, object],
    cutoff: float, threads: int, flank: int,
) -> tuple[list[dict[str, object]], dict[str, object], list[dict[str, object]]]:
    transcripts = dict(read_fasta(transcriptome_fasta, allow_ambiguous=True))
    windows = window_candidates(bowtie_table, bowtie_label, transcripts, flank)
    config = DegradomeConfig(srna_text="", srna_fasta=None, transcript_text="", transcript_fasta=None, samples=(), output_dir=Path("."), ignore_query_pos1=True, slice_positions=(10,), mfe_ratio_cutoff=cutoff, sort_by="mfe_ratio", threads=1)

    def one(window: dict[str, object]) -> list[dict[str, object]]:
        guide = guides[str(window["siRNA_id"])]
        query = FastaRecord(str(window["window_id"]), guide.sequence)
        hits = target_hits_for_query(query, {str(window["window_id"]): str(window["window_sequence"])}, config)
        rows: list[dict[str, object]] = []
        for hit in hits:
            start = int(window["window_start"]) + int(hit.t_start) - 1
            stop = int(window["window_start"]) + int(hit.t_stop) - 1
            if stop < int(window["bowtie_start"]) or start > int(window["bowtie_stop"]):
                continue
            rows.append({
                "reference_transcriptome": name, "siRNA_id": guide.query_id, "dsrna_id": guide.dsrna_id,
                "dsrna_position": guide.dsrna_position, "strand": guide.strand, "siRNA_sequence": guide.sequence,
                "transcript_id": window["transcript_id"], "bowtie_start": window["bowtie_start"], "bowtie_stop": window["bowtie_stop"],
                "rnaplex_target_start": start, "rnaplex_target_stop": stop, "mfe_ratio": hit.mfe_ratio,
                "mfe_perfect": hit.mfe_perfect, "mfe_site": hit.mfe_site, "allen_score": hit.allen_score,
                "mismatches": hit.mismatch_count, "gu_wobbles": hit.gu_wobble_count, "bulges": hit.bulge_count,
                "match_pattern": hit.match_pattern,
            })
        return rows

    rows: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=min(max(1, threads), len(windows) or 1)) as pool:
        futures = [pool.submit(one, window) for window in windows]
        for future in as_completed(futures):
            rows.extend(future.result())
    unique_rows: dict[tuple[object, ...], dict[str, object]] = {}
    for row in rows:
        key = (row["reference_transcriptome"], row["siRNA_id"], row["transcript_id"], row["bowtie_start"], row["rnaplex_target_start"], row["rnaplex_target_stop"])
        old = unique_rows.get(key)
        if old is None or float(row["mfe_ratio"]) > float(old["mfe_ratio"]):
            unique_rows[key] = row
    rows = list(unique_rows.values())
    by_sirna: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_sirna[str(row["siRNA_id"])].append(row)
    per_sirna: list[dict[str, object]] = []
    for guide in guides.values():
        hits = by_sirna.get(guide.query_id, [])
        ratios = [float(hit["mfe_ratio"]) for hit in hits]
        energies = [float(hit["mfe_site"]) for hit in hits]
        per_sirna.append({
            "reference_transcriptome": name, "siRNA_id": guide.query_id, "dsrna_position": guide.dsrna_position,
            "strand": guide.strand, "siRNA_sequence": guide.sequence,
            "window_overlapping_rnaplex_hit_count": len(hits), "mean_mfe_ratio": mean(ratios) if ratios else "",
            "median_mfe_ratio": median(ratios) if ratios else "", "maximum_mfe_ratio": max(ratios) if ratios else "",
            "mean_mfe_site": mean(energies) if energies else "",
        })
    ratios = [float(row["mfe_ratio"]) for row in rows]
    summary = {
        "reference_transcriptome": name, "bowtie_window_count": len(windows), "overlapping_rnaplex_hit_count": len(rows),
        "unique_target_transcript_count": len({str(row["transcript_id"]) for row in rows}),
        "mean_mfe_ratio": mean(ratios) if ratios else "", "median_mfe_ratio": median(ratios) if ratios else "",
    }
    return rows, summary, per_sirna


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
    parser.add_argument("--flank", type=int, default=10)
    parser.add_argument("--mfe-ratio-floor", type=float, default=0.01)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    records = read_fasta(args.sirnas)
    guides = {guide.query_id: guide for guide in make_direct_sirna_guides(records, 21, 1, 1) if guide.variant == "original"}
    sources = (("CSFB", args.focal_transcriptome, args.csfb_bowtie_table, "focal"), ("Tribolium", args.tribolium_transcriptome, args.tribolium_bowtie_table, "control"), ("Drosophila", args.drosophila_transcriptome, args.drosophila_bowtie_table, "control"))
    all_rows: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []
    per_sirna: list[dict[str, object]] = []
    for source in sources:
        rows, summary, guide_rows = refine_reference(*source, guides, args.mfe_ratio_floor, args.threads, args.flank)
        all_rows.extend(rows); summaries.append(summary); per_sirna.extend(guide_rows)
    write_rows(args.output_dir / "bowtie_window_rnaplex_targets.tsv", all_rows)
    write_rows(args.output_dir / "bowtie_window_rnaplex_summary.tsv", summaries)
    write_rows(args.output_dir / "bowtie_window_rnaplex_per_sirna.tsv", per_sirna)
    print("Window-restricted RNAplex refinement completed.")


if __name__ == "__main__":
    main()
