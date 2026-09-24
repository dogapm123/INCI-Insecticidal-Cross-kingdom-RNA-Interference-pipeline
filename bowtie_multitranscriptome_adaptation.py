"""Rank siRNA adaptation across multiple transcriptomes with exhaustive short-guide mapping.

Every transcriptome is treated symmetrically.  Original siRNAs are compared
with their matched, dinucleotide-preserving shuffle controls within each
transcriptome; densities are normalized by that transcriptome's searchable
guide windows.  Pairwise transcriptome contrasts are then optional outputs,
not inputs that define a focal/control hierarchy.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, pstdev
from typing import Any

from degradome_analysis import write_rows
from dsrna_adaptation import benjamini_hochberg, make_direct_sirna_guides, make_guides, read_fasta, transcriptome_search_space
from bowtie_adaptation import DEFAULT_BOWTIE_INDEX_CACHE, bowtie_map, collapse_wobble_mappings, hamming_scan_map, metrics, wobble_compatible_variants, wobble_hamming_scan_map


CORE_PRETRIMMED_DM_MRNA_FASTA = DEFAULT_BOWTIE_INDEX_CACHE / "poly_a_trimmed" / "Dm_mRNA.polyAtrim10_94101bd99269f494.fasta"
CORE_PRETRIMMED_TC_MRNA_FASTA = DEFAULT_BOWTIE_INDEX_CACHE / "poly_a_trimmed" / "OGS3_mRNA.polyAtrim10_c91caca05aef50a3.fasta"
CORE_PRETRIMMED_CONTROL_FASTAS = frozenset({CORE_PRETRIMMED_DM_MRNA_FASTA, CORE_PRETRIMMED_TC_MRNA_FASTA})


@dataclass(frozen=True)
class MultiTranscriptomeBowtieConfig:
    sirna_fasta: Path
    transcriptomes: tuple[tuple[str, Path], ...]
    output_dir: Path
    input_mode: str = "direct_sirnas"
    guide_length: int = 21
    shuffle_count: int = 3
    mismatches: int = 3
    ignore_query_pos1: bool = True
    threads: int = 8
    seed: int = 1
    global_permutations: int = 100_000
    wobble_max_pairs: int = 0
    trim_terminal_poly_a: bool = False
    poly_a_min_length: int = 10
    mapping_backend: str = "bowtie1"
    index_cache_dir: Path | None = DEFAULT_BOWTIE_INDEX_CACHE


def _empirical_pvalue(observed: float, null: list[float]) -> float:
    return (sum(value >= observed for value in null) + 1) / (len(null) + 1)


def _parse_transcriptome(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("Transcriptomes must be given as NAME=FASTA_PATH.")
    name, raw_path = value.split("=", 1)
    name = name.strip()
    path = Path(raw_path).expanduser()
    if not name:
        raise argparse.ArgumentTypeError("Transcriptome name cannot be empty.")
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"Transcriptome FASTA does not exist: {path}")
    return name, path


def _parse_input_set(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("Input sets must be given as SET_NAME=FASTA_PATH.")
    name, raw_path = value.split("=", 1)
    name = name.strip()
    path = Path(raw_path).expanduser()
    if not name:
        raise argparse.ArgumentTypeError("Input-set name cannot be empty.")
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"Input-set FASTA does not exist: {path}")
    return name, path


def read_direct_sirna_records(path: Path) -> list[tuple[str, str]]:
    """Read strict FASTA, with a narrowly scoped repair for headerless-list exports.

    Some small-RNA export tools write the first FASTA identifier with ``>`` but
    omit it on subsequent identifier lines.  Such a line cannot be DNA and is
    therefore unambiguously a record label; repairing it here keeps source
    files immutable while preserving strict sequence validation.
    """
    try:
        return read_fasta(path)
    except ValueError as error:
        lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
        repaired: list[tuple[str, str]] = []
        index = 0
        while index < len(lines):
            label = lines[index]
            if label.startswith(">"):
                label = label[1:].split()[0]
            elif set(label.upper()) <= {"A", "C", "G", "T", "U", "N"}:
                raise error
            if not label or index + 1 >= len(lines):
                raise error
            sequence = lines[index + 1]
            if sequence.startswith(">") or not set(sequence.upper()) <= {"A", "C", "G", "T", "U", "N"}:
                raise error
            repaired.append((label, sequence.upper().replace("U", "T")))
            index += 2
        if not repaired:
            raise error
        return repaired


def trim_terminal_poly_a(records: list[tuple[str, str]], minimum_length: int) -> tuple[list[tuple[str, str]], int, int]:
    """Remove only a terminal 3′ A run meeting the requested minimum length."""
    trimmed_records: list[tuple[str, str]] = []
    trimmed_transcripts = 0
    removed_bases = 0
    for name, sequence in records:
        trim_start = len(sequence)
        while trim_start and sequence[trim_start - 1] == "A":
            trim_start -= 1
        tail_length = len(sequence) - trim_start
        if tail_length >= minimum_length:
            trimmed_records.append((name, sequence[:trim_start]))
            trimmed_transcripts += 1
            removed_bases += tail_length
        else:
            trimmed_records.append((name, sequence))
    return trimmed_records, trimmed_transcripts, removed_bases


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in records), encoding="utf-8")


def cached_poly_a_trimmed_reference(
    source_fasta: Path, trimmed_records: list[tuple[str, str]], minimum_length: int,
    cache_dir: Path = DEFAULT_BOWTIE_INDEX_CACHE / "poly_a_trimmed",
) -> Path:
    """Persist a deterministic trimmed reference so its Bowtie index is reusable."""
    source = source_fasta.expanduser().resolve()
    stat = source.stat()
    fingerprint = {
        "source": str(source), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns,
        "terminal_poly_a_min_length": minimum_length, "format": "inci-poly-a-trim-v1",
    }
    token = hashlib.sha256(json.dumps(fingerprint, sort_keys=True).encode("utf-8")).hexdigest()[:16]
    stem = "".join(character if character.isalnum() or character in "._-" else "_" for character in source.stem).strip("._") or "reference"
    destination = cache_dir / f"{stem}.polyAtrim{minimum_length}_{token}.fasta"
    manifest = destination.with_suffix(destination.suffix + ".json")
    if destination.is_file() and manifest.is_file():
        try:
            if json.loads(manifest.read_text(encoding="utf-8")) == fingerprint:
                return destination
        except (OSError, json.JSONDecodeError):
            pass
    destination.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(trimmed_records, destination)
    manifest.write_text(json.dumps(fingerprint, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return destination


def is_core_pretrimmed_control(fasta: Path) -> bool:
    return fasta.expanduser().resolve() in {path.resolve() for path in CORE_PRETRIMMED_CONTROL_FASTAS}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_multi_transcriptome_adaptation(config: MultiTranscriptomeBowtieConfig) -> dict[str, Any]:
    if len(config.transcriptomes) < 2:
        raise ValueError("At least two transcriptomes are required for comparative ranking.")
    names = [name for name, _ in config.transcriptomes]
    if len(names) != len(set(names)):
        raise ValueError("Transcriptome names must be unique.")
    if not 0 <= config.mismatches <= 3:
        raise ValueError("Bowtie1 mismatch limit must be between 0 and 3.")
    if config.mapping_backend not in {"hamming_scan", "bowtie1"}:
        raise ValueError("Mapping backend must be hamming_scan or bowtie1.")
    if config.shuffle_count < 1 or config.global_permutations < 1:
        raise ValueError("At least one shuffle and one permutation are required.")
    if config.poly_a_min_length < 1:
        raise ValueError("The terminal poly(A) minimum length must be positive.")

    output = config.output_dir
    (output / "inputs").mkdir(parents=True, exist_ok=True)
    if config.input_mode not in {"direct_sirnas", "dsrna_windows"}:
        raise ValueError("Input mode must be direct_sirnas or dsrna_windows.")
    records = read_direct_sirna_records(config.sirna_fasta) if config.input_mode == "direct_sirnas" else read_fasta(config.sirna_fasta)
    preserve_shuffle_prefix = int(config.ignore_query_pos1)
    guides = (
        make_direct_sirna_guides(records, None, config.shuffle_count, config.seed, preserve_shuffle_prefix)
        if config.input_mode == "direct_sirnas"
        else make_guides(records, config.guide_length, config.shuffle_count, config.seed, preserve_shuffle_prefix)
    )
    originals = [guide for guide in guides if guide.variant == "original"]
    by_id = {guide.query_id: guide for guide in guides}
    effective = {guide.query_id: guide.sequence[1:] if config.ignore_query_pos1 else guide.sequence for guide in guides}
    effective_lengths = {guide.query_id: len(effective[guide.query_id]) for guide in guides}
    effective_length_values = sorted(set(effective_lengths.values()))
    if config.wobble_max_pairs and config.mapping_backend == "bowtie1":
        mapped_queries, variant_metadata = wobble_compatible_variants(effective, config.wobble_max_pairs, preserve_shuffle_prefix)
    else:
        mapped_queries, variant_metadata = effective, None

    search_spaces: dict[str, dict[int, int]] = {}
    transcript_counts: dict[str, int] = {}
    mapped_fastas: dict[str, Path] = {}
    poly_a_reference_modes: dict[str, str] = {}
    trimmed_transcript_counts: dict[str, int] = {}
    trimmed_base_counts: dict[str, int] = {}
    hits_by_transcriptome: dict[str, list[dict[str, Any]]] = {}
    by_query: dict[str, dict[str, list[dict[str, Any]]]] = {}
    hit_rows: list[dict[str, Any]] = []
    for name, fasta in config.transcriptomes:
        reference = read_fasta(fasta, allow_ambiguous=True)
        if is_core_pretrimmed_control(fasta):
            mapped_fastas[name] = fasta
            trimmed_transcript_counts[name] = 0
            trimmed_base_counts[name] = 0
            poly_a_reference_modes[name] = "core_pretrimmed"
        elif config.trim_terminal_poly_a:
            reference, trimmed_transcript_counts[name], trimmed_base_counts[name] = trim_terminal_poly_a(reference, config.poly_a_min_length)
            mapped_fastas[name] = cached_poly_a_trimmed_reference(fasta, reference, config.poly_a_min_length)
            poly_a_reference_modes[name] = "cached_trimmed"
        else:
            mapped_fastas[name] = fasta
            trimmed_transcript_counts[name] = 0
            trimmed_base_counts[name] = 0
            poly_a_reference_modes[name] = "untrimmed"
        search_spaces[name] = {length: transcriptome_search_space(reference, length) for length in effective_length_values}
        transcript_counts[name] = len(reference)
        if config.mapping_backend == "hamming_scan":
            hits = (
                wobble_hamming_scan_map(effective, reference, config.mismatches, config.wobble_max_pairs, preserve_shuffle_prefix)
                if config.wobble_max_pairs
                else hamming_scan_map(effective, reference, config.mismatches)
            )
        else:
            hits = bowtie_map(
                mapped_queries, mapped_fastas[name], output, config.mismatches, config.threads, name, variant_metadata,
                index_cache_dir=config.index_cache_dir,
            )
        if config.wobble_max_pairs and config.mapping_backend == "bowtie1":
            hits = collapse_wobble_mappings(hits)
        for row in hits:
            row["transcriptome"] = name
            row["mapping_backend"] = config.mapping_backend
            guide = by_id[str(row["siRNA_id"])]
            hit_rows.append({**row, "dsrna_id": guide.dsrna_id, "dsrna_position": guide.dsrna_position, "strand": guide.strand, "variant": guide.variant, "siRNA_sequence": guide.sequence,
                             "effective_guide_length": effective_lengths[guide.query_id], "searchable_guide_windows": search_spaces[name][effective_lengths[guide.query_id]]})
        hits_by_transcriptome[name] = hits
        grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in hits:
            grouped[str(row["siRNA_id"])].append(row)
        by_query[name] = grouped

    per_sirna_rows: list[dict[str, Any]] = []
    original_densities: dict[str, float] = {}
    original_counts: dict[str, int] = {}
    original_targets: dict[str, int] = {}
    shuffle_densities: dict[str, list[list[float]]] = {}
    for name in names:
        original_counts[name] = sum(len(by_query[name][guide.query_id]) for guide in originals)
        original_densities[name] = sum(
            float(metrics(by_query[name][guide.query_id], search_spaces[name][effective_lengths[guide.query_id]])["alignment_density_per_million_sites"])
            for guide in originals
        )
        original_targets[name] = len({str(hit["transcript_id"]) for guide in originals for hit in by_query[name][guide.query_id]})
        per_original_shuffle_densities: list[list[float]] = []
        for guide in originals:
            prefix = guide.query_id.rsplit("|", 1)[0]
            space = search_spaces[name][effective_lengths[guide.query_id]]
            original_metric = metrics(by_query[name][guide.query_id], space)
            shuffles = [by_id[f"{prefix}|shuffle{index}"] for index in range(1, config.shuffle_count + 1)]
            shuffle_metrics = [metrics(by_query[name][shuffle.query_id], space) for shuffle in shuffles]
            per_original_shuffle_densities.append([float(item["alignment_density_per_million_sites"]) for item in shuffle_metrics])
            per_sirna_rows.append({
                "transcriptome": name, "mapping_backend": config.mapping_backend,
                "dsrna_id": guide.dsrna_id, "dsrna_position": guide.dsrna_position,
                "strand": guide.strand, "siRNA_id": guide.query_id, "siRNA_sequence": guide.sequence,
                "effective_guide_length": effective_lengths[guide.query_id], "searchable_guide_windows": space,
                "shannon_entropy": guide.entropy,
                "original_alignment_count": original_metric["alignment_count"],
                "original_target_transcript_count": original_metric["target_transcript_count"],
                "original_target_transcripts": original_metric["target_transcripts"],
                "original_density_per_million_sites": original_metric["alignment_density_per_million_sites"],
                "shuffled_alignment_counts": ",".join(str(item["alignment_count"]) for item in shuffle_metrics),
                "shuffled_densities_per_million_sites": ",".join(str(item["alignment_density_per_million_sites"]) for item in shuffle_metrics),
            })
        shuffle_densities[name] = per_original_shuffle_densities

    # Per-siRNA significance has two complementary null models.  The paired
    # model uses only the siRNA's own dinucleotide shuffles (faithful but
    # necessarily discrete with a small shuffle count).  The pooled model
    # pools shuffled siRNAs of the same effective length within a transcriptome
    # to provide a higher-resolution, length-matched empirical null.
    pooled_shuffles_by_transcriptome_and_length: dict[tuple[str, int], list[float]] = defaultdict(list)
    for row in per_sirna_rows:
        key = (str(row["transcriptome"]), int(row["effective_guide_length"]))
        pooled_shuffles_by_transcriptome_and_length[key].extend(
            float(value) for value in str(row["shuffled_densities_per_million_sites"]).split(",") if value
        )
    for row in per_sirna_rows:
        matched_null = [
            float(value) for value in str(row["shuffled_densities_per_million_sites"]).split(",") if value
        ]
        pooled_null = pooled_shuffles_by_transcriptome_and_length[
            (str(row["transcriptome"]), int(row["effective_guide_length"]))
        ]
        observed = float(row["original_density_per_million_sites"])
        matched_mean = mean(matched_null)
        pooled_mean = mean(pooled_null)
        row.update({
            "matched_shuffle_count": len(matched_null),
            "mean_matched_shuffle_density_per_million_sites": matched_mean,
            "paired_adaptation_excess_density_per_million_sites": observed - matched_mean,
            "paired_pvalue_vs_matched_shuffles": _empirical_pvalue(observed, matched_null),
            "pooled_length_matched_shuffle_count": len(pooled_null),
            "mean_pooled_length_matched_shuffle_density_per_million_sites": pooled_mean,
            "pooled_length_matched_adaptation_excess_density_per_million_sites": observed - pooled_mean,
            "pooled_length_matched_pvalue_vs_shuffles": _empirical_pvalue(observed, pooled_null),
        })
    paired_qvalues = benjamini_hochberg([float(row["paired_pvalue_vs_matched_shuffles"]) for row in per_sirna_rows])
    pooled_qvalues = benjamini_hochberg([float(row["pooled_length_matched_pvalue_vs_shuffles"]) for row in per_sirna_rows])
    for row, paired_qvalue, pooled_qvalue in zip(per_sirna_rows, paired_qvalues, pooled_qvalues):
        row["paired_fdr_vs_matched_shuffles"] = paired_qvalue
        row["pooled_length_matched_fdr_vs_shuffles"] = pooled_qvalue

    rng = random.Random(config.seed)
    null_densities = {name: [] for name in names}
    for _ in range(config.global_permutations):
        selected = [rng.randrange(config.shuffle_count) for _ in originals]
        for name in names:
            null_densities[name].append(sum(values[index] for values, index in zip(shuffle_densities[name], selected)))

    global_rows: list[dict[str, Any]] = []
    for name in names:
        null = null_densities[name]
        null_mean = mean(null)
        global_rows.append({
            "transcriptome": name, "input_fasta": str(dict(config.transcriptomes)[name]),
            "mapping_backend": config.mapping_backend,
            "mapped_fasta": str(mapped_fastas[name]), "terminal_poly_a_trimming": config.trim_terminal_poly_a,
            "poly_a_reference_mode": poly_a_reference_modes[name],
            "poly_a_min_length": config.poly_a_min_length if config.trim_terminal_poly_a else 0,
            "poly_a_trimmed_transcript_count": trimmed_transcript_counts[name], "poly_a_bases_removed": trimmed_base_counts[name],
            "transcript_count": transcript_counts[name],
            "searchable_guide_windows": search_spaces[name][effective_length_values[0]] if len(effective_length_values) == 1 else "",
            "searchable_guide_windows_by_effective_length": ",".join(f"{length}:{search_spaces[name][length]}" for length in effective_length_values),
            "normalization_mode": "per-siRNA effective-length search space" if len(effective_length_values) > 1 else "shared effective-length search space",
            "original_alignment_count": original_counts[name], "original_target_transcript_count": original_targets[name],
            "original_density_per_million_sites": original_densities[name],
            "mean_matched_shuffle_density_per_million_sites": null_mean,
            "sd_matched_shuffle_density_per_million_sites": pstdev(null),
            "adaptation_excess_density_per_million_sites": original_densities[name] - null_mean,
            "adaptation_fold_over_matched_shuffle": original_densities[name] / null_mean if null_mean else float("inf"),
            "pvalue_vs_matched_shuffles": _empirical_pvalue(original_densities[name], null),
        })
    qvalues = benjamini_hochberg([float(row["pvalue_vs_matched_shuffles"]) for row in global_rows])
    for row, qvalue in zip(global_rows, qvalues):
        row["fdr_vs_matched_shuffles"] = qvalue
    global_rows.sort(key=lambda row: float(row["adaptation_excess_density_per_million_sites"]), reverse=True)
    for rank, row in enumerate(global_rows, start=1):
        row["adaptation_rank"] = rank

    pairwise_rows: list[dict[str, Any]] = []
    for left_index, left in enumerate(names):
        for right in names[left_index + 1:]:
            observed = original_densities[left] - original_densities[right]
            null = [a - b for a, b in zip(null_densities[left], null_densities[right])]
            pairwise_rows.append({
                "left_transcriptome": left, "right_transcriptome": right,
                "original_density_difference_left_minus_right": observed,
                "mean_shuffle_density_difference_left_minus_right": mean(null),
                "excess_difference_left_minus_right": observed - mean(null),
                "pvalue_left_more_adapted_than_right": _empirical_pvalue(observed, null),
            })
    pairwise_qvalues = benjamini_hochberg([float(row["pvalue_left_more_adapted_than_right"]) for row in pairwise_rows])
    for row, qvalue in zip(pairwise_rows, pairwise_qvalues):
        row["fdr_left_more_adapted_than_right"] = qvalue

    write_rows(output / "tables" / "multi_transcriptome_adaptation_ranking.tsv", global_rows)
    write_rows(output / "tables" / "multi_transcriptome_adaptation_pairwise.tsv", pairwise_rows)
    write_rows(output / "tables" / "multi_transcriptome_adaptation_per_sirna.tsv", per_sirna_rows)
    write_rows(output / "tables" / "multi_transcriptome_adaptation_target_hits.tsv", hit_rows)
    write_rows(output / "tables" / "poly_a_trimming_summary.tsv", [
        {"transcriptome": name, "input_fasta": str(dict(config.transcriptomes)[name]), "mapped_fasta": str(mapped_fastas[name]),
         "mapping_backend": config.mapping_backend,
         "terminal_poly_a_trimming": config.trim_terminal_poly_a, "poly_a_min_length": config.poly_a_min_length if config.trim_terminal_poly_a else 0,
         "poly_a_reference_mode": poly_a_reference_modes[name],
         "trimmed_transcript_count": trimmed_transcript_counts[name], "bases_removed": trimmed_base_counts[name]}
        for name in names
    ])
    return {"outdir": str(output), "transcriptome_count": len(names), "sirna_count": len(originals), "mapping_backend": config.mapping_backend, "ranking": str(output / "tables" / "multi_transcriptome_adaptation_ranking.tsv")}


def run_batch_multi_transcriptome_adaptation(
    input_sets: tuple[tuple[str, Path], ...],
    transcriptomes: tuple[tuple[str, Path], ...],
    output_dir: Path,
    *,
    input_mode: str = "direct_sirnas",
    guide_length: int = 21,
    shuffle_count: int = 3,
    mismatches: int = 3,
    ignore_query_pos1: bool = True,
    threads: int = 8,
    seed: int = 1,
    global_permutations: int = 100_000,
    wobble_max_pairs: int = 0,
    trim_terminal_poly_a: bool = False,
    poly_a_min_length: int = 10,
    mapping_backend: str = "bowtie1",
) -> dict[str, Any]:
    """Run independent, reproducible adaptation analyses for multiple input sets.

    Each input FASTA receives its own matched shuffle null model.  The batch
    tables add an ``input_set`` column so identical siRNA identifiers remain
    distinguishable after results are merged.
    """
    if not input_sets:
        raise ValueError("Provide at least one named siRNA or dsRNA FASTA input set.")
    raw_names = [name for name, _ in input_sets]
    safe_names = ["".join(character if character.isalnum() or character in "._-" else "_" for character in name).strip("._") for name in raw_names]
    if any(not name for name in safe_names) or len(safe_names) != len(set(safe_names)):
        raise ValueError("Input-set names must be unique and contain at least one letter, number, '.', '_' or '-'.")
    for name, path in input_sets:
        if not path.is_file():
            raise ValueError(f"Input set {name!r} does not exist: {path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    tables = output_dir / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    table_names = (
        "multi_transcriptome_adaptation_ranking.tsv",
        "multi_transcriptome_adaptation_pairwise.tsv",
        "multi_transcriptome_adaptation_per_sirna.tsv",
        "multi_transcriptome_adaptation_target_hits.tsv",
    )
    merged: dict[str, list[dict[str, str]]] = {name: [] for name in table_names}
    manifests: list[dict[str, Any]] = []
    set_results: list[dict[str, Any]] = []
    for index, ((set_name, input_fasta), safe_set_name) in enumerate(zip(input_sets, safe_names), start=1):
        set_dir = output_dir / "sets" / safe_set_name
        result = run_multi_transcriptome_adaptation(MultiTranscriptomeBowtieConfig(
            sirna_fasta=input_fasta,
            transcriptomes=transcriptomes,
            output_dir=set_dir,
            input_mode=input_mode,
            guide_length=guide_length,
            shuffle_count=shuffle_count,
            mismatches=mismatches,
            ignore_query_pos1=ignore_query_pos1,
            threads=threads,
            seed=seed + index - 1,
            global_permutations=global_permutations,
            wobble_max_pairs=wobble_max_pairs,
            trim_terminal_poly_a=trim_terminal_poly_a,
            poly_a_min_length=poly_a_min_length,
            mapping_backend=mapping_backend,
        ))
        for table_name in table_names:
            for row in read_rows(set_dir / "tables" / table_name):
                merged[table_name].append({"input_set": set_name, "input_set_fasta": str(input_fasta), **row})
        manifests.append({
            "input_set": set_name,
            "input_fasta": str(input_fasta),
            "input_mode": input_mode,
            "set_output_dir": str(set_dir),
            "original_siRNA_count": result["sirna_count"],
            "transcriptome_count": result["transcriptome_count"],
            "random_seed": seed + index - 1,
            "mapping_backend": mapping_backend,
        })
        set_results.append({"input_set": set_name, "output_dir": str(set_dir), **result})

    for source_name, rows in merged.items():
        destination = tables / f"batch_{source_name}"
        write_rows(destination, rows)
    write_rows(tables / "batch_input_sets.tsv", manifests)
    return {
        "outdir": str(output_dir),
        "input_set_count": len(input_sets),
        "transcriptome_count": len(transcriptomes),
        "sirna_count": sum(int(result["sirna_count"]) for result in set_results),
        "mapping_backend": mapping_backend,
        "set_results": set_results,
        "ranking": str(tables / "batch_multi_transcriptome_adaptation_ranking.tsv"),
        "per_sirna": str(tables / "batch_multi_transcriptome_adaptation_per_sirna.tsv"),
        "target_hits": str(tables / "batch_multi_transcriptome_adaptation_target_hits.tsv"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sirnas", type=Path, help="One siRNA or dsRNA FASTA (single-set mode).")
    parser.add_argument("--sirna-set", action="append", type=_parse_input_set, metavar="SET_NAME=FASTA", help="Repeat to batch multiple siRNA or dsRNA FASTA input sets.")
    parser.add_argument("--transcriptome", required=True, action="append", type=_parse_transcriptome, metavar="NAME=FASTA", help="Repeat for every transcriptome to rank.")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--input-mode", choices=("direct_sirnas", "dsrna_windows"), default="direct_sirnas")
    parser.add_argument("--shuffles", type=int, default=3)
    parser.add_argument("--mismatches", type=int, default=3)
    parser.add_argument("--wobble-max-pairs", type=int, default=0)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--permutations", type=int, default=100_000)
    parser.add_argument("--trim-terminal-poly-a", action="store_true", help="Trim terminal 3′ poly(A) runs from working transcriptome copies.")
    parser.add_argument("--poly-a-min-length", type=int, default=10)
    parser.add_argument("--mapping-backend", choices=("hamming_scan", "bowtie1"), default="bowtie1")
    args = parser.parse_args()
    if args.sirna_set:
        if args.sirnas:
            parser.error("Use either --sirnas or one or more --sirna-set entries, not both.")
        result = run_batch_multi_transcriptome_adaptation(
            tuple(args.sirna_set), tuple(args.transcriptome), args.output_dir,
            input_mode=args.input_mode, shuffle_count=args.shuffles, mismatches=args.mismatches,
            wobble_max_pairs=args.wobble_max_pairs, threads=args.threads, global_permutations=args.permutations,
            trim_terminal_poly_a=args.trim_terminal_poly_a, poly_a_min_length=args.poly_a_min_length,
            mapping_backend=args.mapping_backend,
        )
        print(f"Ranked {result['sirna_count']} original siRNAs from {result['input_set_count']} input set(s) across {result['transcriptome_count']} transcriptomes: {result['ranking']}")
        return
    if args.sirnas is None:
        parser.error("Provide --sirnas for one input set or repeat --sirna-set for batch analysis.")
    result = run_multi_transcriptome_adaptation(MultiTranscriptomeBowtieConfig(
        sirna_fasta=args.sirnas, transcriptomes=tuple(args.transcriptome), output_dir=args.output_dir,
        input_mode=args.input_mode, shuffle_count=args.shuffles, mismatches=args.mismatches,
        wobble_max_pairs=args.wobble_max_pairs, threads=args.threads, global_permutations=args.permutations,
        trim_terminal_poly_a=args.trim_terminal_poly_a, poly_a_min_length=args.poly_a_min_length,
        mapping_backend=args.mapping_backend,
    ))
    print(f"Ranked {result['sirna_count']} original siRNAs across {result['transcriptome_count']} transcriptomes: {result['ranking']}")


if __name__ == "__main__":
    main()
