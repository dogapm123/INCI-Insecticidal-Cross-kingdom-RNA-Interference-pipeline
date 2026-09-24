"""Exhaustive short-guide adaptation statistics for defined siRNAs or dsRNA windows."""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import random
import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from statistics import mean
from typing import Any

from degradome_analysis import build_bowtie_index, write_rows
from dsrna_adaptation import (
    benjamini_hochberg,
    empirical_pvalue,
    make_direct_sirna_guides,
    make_guides,
    read_fasta,
    reverse_complement,
    transcriptome_search_space,
)


_BASE_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
DEFAULT_BOWTIE_INDEX_CACHE = Path(__file__).resolve().parent / "outputs" / "bowtie_reference_indexes"


def _encode_sequence(sequence: str) -> int:
    encoded = 0
    for base in sequence:
        encoded = (encoded << 2) | _BASE_CODE[base]
    return encoded


def _window_encodings(sequence: str, length: int) -> list[int | None]:
    """Encode every fixed-length A/C/G/T window; ambiguous windows are None."""
    if length < 1 or length > len(sequence):
        return []
    mask = (1 << (2 * length)) - 1
    encoded = 0
    valid_run = 0
    windows: list[int | None] = []
    for index, base in enumerate(sequence):
        code = _BASE_CODE.get(base)
        if code is None:
            encoded = 0
            valid_run = 0
        else:
            encoded = ((encoded << 2) | code) & mask
            valid_run += 1
        if index >= length - 1:
            windows.append(encoded if valid_run >= length else None)
    return windows


@lru_cache(maxsize=None)
def _position_mask(length: int) -> int:
    return ((1 << (2 * length)) - 1) // 3


def _hamming_distance_encoded(left: int, right: int, length: int) -> int:
    differences = left ^ right
    return ((differences | (differences >> 1)) & _position_mask(length)).bit_count()


def _encoded_difference_positions(left: int, right: int, length: int) -> list[int]:
    """Return left-to-right base positions that differ between two 2-bit words."""
    differences = left ^ right
    changed = (differences | (differences >> 1)) & _position_mask(length)
    positions: list[int] = []
    while changed:
        lowest = changed & -changed
        positions.append(length - 1 - ((lowest.bit_length() - 1) // 2))
        changed -= lowest
    return positions


@dataclass(frozen=True)
class BowtieAdaptationConfig:
    sirna_fasta: Path
    focal_transcriptome_fasta: Path
    control_transcriptome_fasta: Path
    output_dir: Path
    input_mode: str = "direct_sirnas"
    guide_length: int = 21
    shuffle_count: int = 3
    mismatches: int = 3
    fdr_cutoff: float = 0.05
    ignore_query_pos1: bool = True
    threads: int = 8
    seed: int = 1
    global_permutations: int = 100_000
    wobble_max_pairs: int = 0


def cached_bowtie_index_prefix(reference: Path, cache_dir: Path = DEFAULT_BOWTIE_INDEX_CACHE, threads: int = 8) -> Path:
    """Return a persistent Bowtie1 index, rebuilding only when FASTA changes."""
    source = reference.expanduser().resolve()
    stat = source.stat()
    fingerprint = {
        "source": str(source),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "format": "bowtie1-v1",
    }
    token = hashlib.sha256(json.dumps(fingerprint, sort_keys=True).encode("utf-8")).hexdigest()[:16]
    stem = "".join(character if character.isalnum() or character in "._-" else "_" for character in source.stem).strip("._") or "reference"
    index_dir = cache_dir / f"{stem}_{token}"
    index_prefix = index_dir / "index"
    manifest = index_dir / "manifest.json"
    expected = [Path(f"{index_prefix}.{suffix}.ebwt") for suffix in ("1", "2", "3", "4", "rev.1", "rev.2")]
    if all(path.is_file() for path in expected) and manifest.is_file():
        try:
            if json.loads(manifest.read_text(encoding="utf-8")) == fingerprint:
                return index_prefix
        except (OSError, json.JSONDecodeError):
            pass
    build_bowtie_index(source, index_prefix, threads)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps(fingerprint, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return index_prefix


def bowtie_map(
    query_sequences: dict[str, str], transcriptome_fasta: Path, output_dir: Path,
    mismatches: int, threads: int, label: str, variant_metadata: dict[str, dict[str, object]] | None = None,
    index_cache_dir: Path | None = None,
) -> list[dict[str, Any]]:
    bowtie = shutil.which("bowtie-align-s") or shutil.which("bowtie")
    if not bowtie or not shutil.which("bowtie-build"):
        raise ValueError("Bowtie1 and bowtie-build are required.")
    query_fasta = output_dir / "inputs" / f"{label}_queries.fasta"
    query_fasta.parent.mkdir(parents=True, exist_ok=True)
    query_fasta.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in query_sequences.items()), encoding="utf-8")
    index_prefix = (
        cached_bowtie_index_prefix(transcriptome_fasta, index_cache_dir, threads)
        if index_cache_dir is not None
        else output_dir / "bowtie_index" / transcriptome_fasta.stem
    )
    if index_cache_dir is None:
        build_bowtie_index(transcriptome_fasta, index_prefix, threads)
    run = subprocess.run(
        [bowtie, "-f", "-v", str(mismatches), "-a", "--nofw", "--sam", "--mm", "-p", str(threads), str(index_prefix), str(query_fasta)],
        text=True, capture_output=True, check=False,
    )
    if run.returncode not in {0, 1}:
        raise ValueError(run.stderr.strip() or "Bowtie1 mapping failed.")
    rows: list[dict[str, Any]] = []
    for line in run.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 11 or fields[0] not in query_sequences or fields[2] == "*":
            continue
        tags = {tag.split(":", 2)[0]: tag.split(":", 2)[2] for tag in fields[11:] if tag.count(":") >= 2}
        generic_mismatches = int(tags.get("NM", tags.get("XM", "0")))
        metadata = (variant_metadata or {}).get(fields[0], {"siRNA_id": fields[0], "wobble_count": 0, "wobble_positions": ""})
        rows.append({
            "siRNA_id": metadata["siRNA_id"], "wobble_variant_id": fields[0],
            "wobble_count": metadata["wobble_count"], "wobble_positions": metadata["wobble_positions"],
            "generic_mismatch_count": generic_mismatches,
            "mapping_strand": "-" if int(fields[1]) & 16 else "+", "transcript_id": fields[2],
            "target_start": int(fields[3]), "aligned_sequence": fields[9],
        })
    return rows


def hamming_scan_map(
    query_sequences: dict[str, str], transcriptome_records: list[tuple[str, str]],
    mismatches: int, variant_metadata: dict[str, dict[str, object]] | None = None,
) -> list[dict[str, Any]]:
    """Return every reverse-complement ungapped site within a Hamming limit.

    This is an index-free equivalent of the Bowtie1 arguments used by this
    module: ``-v N -a --nofw``.  The scanner retains every site with at
    most ``mismatches`` substitutions.  A seed is used only to avoid testing
    impossible windows.  Each retained candidate is checked across the full
    query, so the resulting hit set is exhaustive.
    """
    if mismatches < 0:
        raise ValueError("Mismatch limit cannot be negative.")
    valid_bases = frozenset("ACGT")
    target_sequences = {query_id: reverse_complement(sequence) for query_id, sequence in query_sequences.items()}
    query_lengths = {query_id: len(sequence) for query_id, sequence in target_sequences.items()}
    target_encodings = {query_id: _encode_sequence(sequence) for query_id, sequence in target_sequences.items()}
    seed_lookup_by_length: dict[int, dict[str, list[tuple[str, int]]]] = defaultdict(lambda: defaultdict(list))
    unseeded_queries: list[str] = []
    for query_id, target in target_sequences.items():
        if not target or set(target) - valid_bases:
            raise ValueError(f"Query {query_id!r} must contain only A, C, G and T for Hamming scanning.")
        if len(target) <= mismatches:
            unseeded_queries.append(query_id)
            continue
        block_count = min(mismatches + 1, len(target))
        base_length, longer_blocks = divmod(len(target), block_count)
        offset = 0
        for block_index in range(block_count):
            block_length = base_length + int(block_index < longer_blocks)
            seed_lookup_by_length[block_length][target[offset : offset + block_length]].append((query_id, offset))
            offset += block_length

    rows: list[dict[str, Any]] = []
    for transcript_id, reference in transcriptome_records:
        candidate_sites: set[tuple[str, int]] = set()
        windows_by_length = {length: _window_encodings(reference, length) for length in set(query_lengths.values())}
        for seed_length, seed_lookup in seed_lookup_by_length.items():
            for seed_start in range(max(0, len(reference) - seed_length + 1)):
                for query_id, query_offset in seed_lookup.get(reference[seed_start : seed_start + seed_length], ()):
                    target_start = seed_start - query_offset
                    if 0 <= target_start <= len(reference) - query_lengths[query_id]:
                        candidate_sites.add((query_id, target_start))
        for query_id in unseeded_queries:
            target_length = query_lengths[query_id]
            candidate_sites.update((query_id, start) for start in range(max(0, len(reference) - target_length + 1)))
        for query_id, target_start in sorted(candidate_sites):
            target_length = query_lengths[query_id]
            aligned_encoding = windows_by_length[target_length][target_start]
            if aligned_encoding is None:
                continue
            generic_mismatches = _hamming_distance_encoded(target_encodings[query_id], aligned_encoding, target_length)
            if generic_mismatches > mismatches:
                continue
            metadata = (variant_metadata or {}).get(
                query_id, {"siRNA_id": query_id, "wobble_count": 0, "wobble_positions": ""}
            )
            rows.append({
                "siRNA_id": metadata["siRNA_id"], "wobble_variant_id": query_id,
                "wobble_count": metadata["wobble_count"], "wobble_positions": metadata["wobble_positions"],
                "generic_mismatch_count": generic_mismatches, "mapping_strand": "-",
                # Bowtie's SAM POS is 1-based, and with FLAG 16 its SEQ field
                # is reverse-complemented into the transcript orientation.
                "transcript_id": transcript_id, "target_start": target_start + 1,
                "aligned_sequence": target_sequences[query_id],
            })
    return rows


def wobble_hamming_scan_map(
    query_sequences: dict[str, str], transcriptome_records: list[tuple[str, str]],
    mismatches: int, max_pairs: int, original_offset: int,
) -> list[dict[str, Any]]:
    """Map wobble-aware guides without materializing every wobble variant.

    For a retained site, the best variant is determined directly: a permitted
    wobble base is recorded only when it converts that position into a match.
    This is equivalent to expanding the variants and applying
    :func:`collapse_wobble_mappings`, but avoids scanning the same reference
    once for every possible wobble combination.
    """
    if mismatches < 0 or max_pairs < 0:
        raise ValueError("Mismatch and wobble limits cannot be negative.")
    valid_bases = frozenset("ACGT")
    target_sequences = {query_id: reverse_complement(sequence) for query_id, sequence in query_sequences.items()}
    query_lengths = {query_id: len(sequence) for query_id, sequence in target_sequences.items()}
    target_encodings = {query_id: _encode_sequence(sequence) for query_id, sequence in target_sequences.items()}
    wobble_targets: dict[str, dict[int, tuple[str, int]]] = {}
    seed_lookup_by_length: dict[int, dict[str, list[tuple[str, int]]]] = defaultdict(lambda: defaultdict(list))
    unseeded_queries: list[str] = []
    for query_id, query in query_sequences.items():
        target = target_sequences[query_id]
        if not target or set(target) - valid_bases:
            raise ValueError(f"Query {query_id!r} must contain only A, C, G and T for Hamming scanning.")
        alternatives: dict[int, tuple[str, int]] = {}
        for query_position, base in enumerate(query):
            if base == "G":
                alternatives[len(query) - 1 - query_position] = ("T", query_position)
            elif base == "T":
                alternatives[len(query) - 1 - query_position] = ("G", query_position)
        wobble_targets[query_id] = alternatives
        if len(target) <= mismatches:
            unseeded_queries.append(query_id)
            continue
        # Divide into ``mismatches + 1`` blocks.  At least one block has no
        # generic mismatch.  Generating every permitted wobble combination
        # within that one block guarantees an exact seed without expanding
        # every whole-guide wobble variant.
        block_count = min(mismatches + 1, len(target))
        base_length, longer_blocks = divmod(len(target), block_count)
        offset = 0
        for block_index in range(block_count):
            block_length = base_length + int(block_index < longer_blocks)
            seed = target[offset : offset + block_length]
            wobble_options = [
                (target_position - offset, alternative)
                for target_position, (alternative, _) in alternatives.items()
                if offset <= target_position < offset + block_length
            ]
            for wobble_count in range(min(max_pairs, len(wobble_options)) + 1):
                for chosen in itertools.combinations(wobble_options, wobble_count):
                    seed_bases = list(seed)
                    for local_position, alternative in chosen:
                        seed_bases[local_position] = alternative
                    seed_lookup_by_length[block_length]["".join(seed_bases)].append((query_id, offset))
            offset += block_length

    rows: list[dict[str, Any]] = []
    for transcript_id, reference in transcriptome_records:
        candidate_sites: set[tuple[str, int]] = set()
        windows_by_length = {length: _window_encodings(reference, length) for length in set(query_lengths.values())}
        for seed_length, seed_lookup in seed_lookup_by_length.items():
            for seed_start in range(max(0, len(reference) - seed_length + 1)):
                for query_id, query_offset in seed_lookup.get(reference[seed_start : seed_start + seed_length], ()):
                    target_start = seed_start - query_offset
                    if 0 <= target_start <= len(reference) - query_lengths[query_id]:
                        candidate_sites.add((query_id, target_start))
        for query_id in unseeded_queries:
            target_length = query_lengths[query_id]
            candidate_sites.update((query_id, start) for start in range(max(0, len(reference) - target_length + 1)))
        for query_id, target_start in sorted(candidate_sites):
            target = target_sequences[query_id]
            target_length = len(target)
            aligned_encoding = windows_by_length[target_length][target_start]
            if aligned_encoding is None:
                continue
            if _hamming_distance_encoded(target_encodings[query_id], aligned_encoding, target_length) > mismatches + max_pairs:
                continue
            wobble_positions: list[int] = []
            generic_mismatches = 0
            for target_position in _encoded_difference_positions(target_encodings[query_id], aligned_encoding, target_length):
                observed = reference[target_start + target_position]
                alternative = wobble_targets[query_id].get(target_position)
                if alternative and observed == alternative[0]:
                    wobble_positions.append(alternative[1])
                else:
                    generic_mismatches += 1
            if generic_mismatches > mismatches or len(wobble_positions) > max_pairs:
                continue
            wobble_positions.sort()
            labels = [str(position + 1 + original_offset) for position in wobble_positions]
            variant_target = list(target)
            for target_position, (alternative, query_position) in wobble_targets[query_id].items():
                if query_position in wobble_positions:
                    variant_target[target_position] = alternative
            label = "canonical" if not labels else "w" + "-".join(labels)
            rows.append({
                "siRNA_id": query_id, "wobble_variant_id": f"{query_id}|{label}",
                "wobble_count": len(labels), "wobble_positions": ",".join(labels),
                "generic_mismatch_count": generic_mismatches, "mapping_strand": "-",
                "transcript_id": transcript_id, "target_start": target_start + 1,
                "aligned_sequence": "".join(variant_target),
            })
    return rows


def wobble_compatible_variants(query_sequences: dict[str, str], max_pairs: int, original_offset: int) -> tuple[dict[str, str], dict[str, dict[str, object]]]:
    """Expand G:U/U:G pair possibilities for reverse-complement Bowtie mapping.

    In the Bowtie-mapped orientation, guide G:U wobble is encoded by G→A in
    the guide, and U:G wobble by T→C. Exact Bowtie mapping therefore preserves
    canonical pairing while reporting the substituted wobble locations.
    """
    variants: dict[str, str] = {}
    metadata: dict[str, dict[str, object]] = {}
    substitutions = {"G": "A", "T": "C"}
    for query_id, sequence in query_sequences.items():
        positions = [index for index, base in enumerate(sequence) if base in substitutions]
        for count in range(min(max_pairs, len(positions)) + 1):
            for chosen in itertools.combinations(positions, count):
                sequence_list = list(sequence)
                for position in chosen:
                    sequence_list[position] = substitutions[sequence_list[position]]
                label = "canonical" if not chosen else "w" + "-".join(str(position + 1 + original_offset) for position in chosen)
                variant_id = f"{query_id}|{label}"
                variants[variant_id] = "".join(sequence_list)
                metadata[variant_id] = {
                    "siRNA_id": query_id, "wobble_count": len(chosen),
                    "wobble_positions": ",".join(str(position + 1 + original_offset) for position in chosen),
                }
    return variants, metadata


def collapse_wobble_mappings(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Prevent alternative variant encodings from inflating one original site."""
    best: dict[tuple[object, ...], dict[str, Any]] = {}
    for row in rows:
        key = (row["siRNA_id"], row["transcript_id"], row["target_start"], row["mapping_strand"])
        old = best.get(key)
        rank = (int(row["generic_mismatch_count"]), int(row["wobble_count"]), str(row["wobble_variant_id"]))
        old_rank = (int(old["generic_mismatch_count"]), int(old["wobble_count"]), str(old["wobble_variant_id"])) if old else None
        if old is None or rank < old_rank:
            best[key] = row
    return list(best.values())


def metrics(rows: list[dict[str, Any]], search_space: int) -> dict[str, Any]:
    transcripts = sorted({str(row["transcript_id"]) for row in rows})
    return {
        "alignment_count": len(rows), "target_transcript_count": len(transcripts),
        "target_transcripts": ",".join(transcripts),
        "alignment_density_per_million_sites": len(rows) * 1_000_000.0 / max(1, search_space),
    }


def global_permutation_rows(
    summary: list[dict[str, Any]], focal_space: int, control_space: int, seed: int, permutations: int,
) -> list[dict[str, Any]]:
    """Whole-set matched-shuffle tests, using one shuffle per siRNA per draw."""
    if permutations < 1:
        raise ValueError("At least one global permutation is required.")
    original_focal = sum(int(row["focus_bowtie_alignment_count"]) for row in summary)
    original_control = sum(int(row["control_bowtie_alignment_count"]) for row in summary)
    original_focal_density = original_focal * 1_000_000.0 / focal_space
    original_control_density = original_control * 1_000_000.0 / control_space
    original_delta = original_focal_density - original_control_density
    shuffled_focal = [[int(value) for value in row["shuffled_focus_alignment_counts"].split(",")] for row in summary]
    shuffled_control = [[int(value) for value in row["shuffled_control_alignment_counts"].split(",")] for row in summary]
    rng = random.Random(seed)
    focal_at_least = 0
    delta_at_least = 0
    null_focal_values: list[float] = []
    null_delta_values: list[float] = []
    for _ in range(permutations):
        selected = [rng.randrange(len(values)) for values in shuffled_focal]
        focal_count = sum(values[index] for values, index in zip(shuffled_focal, selected))
        control_count = sum(values[index] for values, index in zip(shuffled_control, selected))
        focal_density = focal_count * 1_000_000.0 / focal_space
        delta = focal_density - control_count * 1_000_000.0 / control_space
        null_focal_values.append(focal_density)
        null_delta_values.append(delta)
        focal_at_least += focal_density >= original_focal_density
        delta_at_least += delta >= original_delta
    p_focal = (focal_at_least + 1) / (permutations + 1)
    p_delta = (delta_at_least + 1) / (permutations + 1)
    q_focal, q_delta = benjamini_hochberg([p_focal, p_delta])
    return [{
        "original_siRNA_count": len(summary), "shuffle_count_per_siRNA": len(shuffled_focal[0]) if shuffled_focal else 0,
        "global_permutations": permutations,
        "original_focus_alignment_count": original_focal,
        "mean_matched_shuffle_focus_alignment_count": mean(value * focal_space / 1_000_000.0 for value in null_focal_values),
        "original_control_alignment_count": original_control,
        "mean_matched_shuffle_control_alignment_count": mean((value - delta) * control_space / 1_000_000.0 for value, delta in zip(null_focal_values, null_delta_values)),
        "original_focus_density_per_million_sites": original_focal_density,
        "mean_shuffle_focus_density_per_million_sites": mean(null_focal_values),
        "original_focus_minus_control_density_per_million_sites": original_delta,
        "mean_shuffle_focus_minus_control_density_per_million_sites": mean(null_delta_values),
        "global_pvalue_focus_vs_shuffled": p_focal,
        "global_pvalue_focus_vs_control": p_delta,
        "global_fdr_focus_vs_shuffled": q_focal,
        "global_fdr_focus_vs_control": q_delta,
    }]


def run_bowtie_adaptation(config: BowtieAdaptationConfig) -> dict[str, Any]:
    if config.input_mode not in {"direct_sirnas", "dsrna_windows"}:
        raise ValueError("Input mode must be direct_sirnas or dsrna_windows.")
    if not 0 <= config.mismatches <= 3:
        raise ValueError("Bowtie1 mismatch limit must be between 0 and 3.")
    if config.wobble_max_pairs < 0:
        raise ValueError("Wobble-pair limit cannot be negative.")
    output = config.output_dir
    (output / "inputs").mkdir(parents=True, exist_ok=True)
    records = read_fasta(config.sirna_fasta)
    focal_records = read_fasta(config.focal_transcriptome_fasta, allow_ambiguous=True)
    control_records = read_fasta(config.control_transcriptome_fasta, allow_ambiguous=True)
    preserve_shuffle_prefix = int(config.ignore_query_pos1)
    guides = make_direct_sirna_guides(records, config.guide_length, config.shuffle_count, config.seed, preserve_shuffle_prefix) if config.input_mode == "direct_sirnas" else make_guides(records, config.guide_length, config.shuffle_count, config.seed, preserve_shuffle_prefix)
    originals = [guide for guide in guides if guide.variant == "original"]
    by_id = {guide.query_id: guide for guide in guides}
    effective = {guide.query_id: guide.sequence[1:] if config.ignore_query_pos1 else guide.sequence for guide in guides}
    effective_length = config.guide_length - int(config.ignore_query_pos1)
    focal_space = transcriptome_search_space(focal_records, effective_length)
    control_space = transcriptome_search_space(control_records, effective_length)
    if config.wobble_max_pairs:
        mapped_queries, variant_metadata = wobble_compatible_variants(effective, config.wobble_max_pairs, int(config.ignore_query_pos1))
        mapping_mismatches = config.mismatches
    else:
        mapped_queries, variant_metadata, mapping_mismatches = effective, None, config.mismatches
    focal_hits = bowtie_map(mapped_queries, config.focal_transcriptome_fasta, output, mapping_mismatches, config.threads, "focal", variant_metadata)
    control_hits = bowtie_map(mapped_queries, config.control_transcriptome_fasta, output, mapping_mismatches, config.threads, "control", variant_metadata)
    if config.wobble_max_pairs:
        focal_hits = collapse_wobble_mappings(focal_hits)
        control_hits = collapse_wobble_mappings(control_hits)
    for row in focal_hits:
        row["transcriptome"] = "focal"
    for row in control_hits:
        row["transcriptome"] = "control"
    focal_by_query: dict[str, list[dict[str, Any]]] = defaultdict(list)
    control_by_query: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in focal_hits:
        focal_by_query[str(row["siRNA_id"])].append(row)
    for row in control_hits:
        control_by_query[str(row["siRNA_id"])].append(row)
    focal_metrics = {guide.query_id: metrics(focal_by_query[guide.query_id], focal_space) for guide in guides}
    control_metrics = {guide.query_id: metrics(control_by_query[guide.query_id], control_space) for guide in guides}
    all_shuffle_focal = [float(focal_metrics[g.query_id]["alignment_density_per_million_sites"]) for g in guides if g.variant == "shuffle"]
    all_shuffle_delta = [float(focal_metrics[g.query_id]["alignment_density_per_million_sites"]) - float(control_metrics[g.query_id]["alignment_density_per_million_sites"]) for g in guides if g.variant == "shuffle"]
    summary: list[dict[str, Any]] = []
    for guide in originals:
        prefix = guide.query_id.rsplit("|", 1)[0]
        shuffled = [by_id[f"{prefix}|shuffle{i}"] for i in range(1, config.shuffle_count + 1)]
        focal = focal_metrics[guide.query_id]
        control = control_metrics[guide.query_id]
        fs = float(focal["alignment_density_per_million_sites"])
        cs = float(control["alignment_density_per_million_sites"])
        shuffle_focal = [float(focal_metrics[g.query_id]["alignment_density_per_million_sites"]) for g in shuffled]
        shuffle_control = [float(control_metrics[g.query_id]["alignment_density_per_million_sites"]) for g in shuffled]
        delta = fs - cs
        null_delta = [a - b for a, b in zip(shuffle_focal, shuffle_control)]
        summary.append({
            "dsrna_id": guide.dsrna_id, "dsrna_position": guide.dsrna_position, "strand": guide.strand,
            "siRNA_id": guide.query_id, "siRNA_sequence": guide.sequence, "shannon_entropy": guide.entropy,
            "focus_bowtie_alignment_count": focal["alignment_count"], "focus_target_transcript_count": focal["target_transcript_count"],
            "focus_target_transcripts": focal["target_transcripts"], "focus_alignment_density_per_million_sites": fs,
            "control_bowtie_alignment_count": control["alignment_count"], "control_target_transcript_count": control["target_transcript_count"],
            "control_target_transcripts": control["target_transcripts"], "control_alignment_density_per_million_sites": cs,
            "adaptation_score": delta - mean(null_delta),
            "pooled_pvalue_vs_shuffled": empirical_pvalue(fs, all_shuffle_focal),
            "pooled_pvalue_focus_vs_control": empirical_pvalue(delta, all_shuffle_delta),
            "shuffled_siRNA_sequences": ",".join(g.sequence for g in shuffled),
            "shuffled_focus_alignment_counts": ",".join(str(focal_metrics[g.query_id]["alignment_count"]) for g in shuffled),
            "shuffled_control_alignment_counts": ",".join(str(control_metrics[g.query_id]["alignment_count"]) for g in shuffled),
        })
    summary.sort(key=lambda row: (str(row["dsrna_id"]), int(row["dsrna_position"]), str(row["strand"])))
    q1 = benjamini_hochberg([float(row["pooled_pvalue_vs_shuffled"]) for row in summary])
    q2 = benjamini_hochberg([float(row["pooled_pvalue_focus_vs_control"]) for row in summary])
    for row, a, b in zip(summary, q1, q2):
        row["fdr_vs_shuffled"] = a
        row["fdr_focus_vs_control"] = b
        row["significant_both"] = a <= config.fdr_cutoff and b <= config.fdr_cutoff
    hit_rows = []
    for row in focal_hits + control_hits:
        guide = by_id[str(row["siRNA_id"])]
        hit_rows.append({**row, "dsrna_id": guide.dsrna_id, "dsrna_position": guide.dsrna_position, "strand": guide.strand, "variant": guide.variant, "siRNA_sequence": guide.sequence})
    locus_rows = []
    for locus in sorted({str(row["dsrna_id"]) for row in summary}):
        rows = [row for row in summary if row["dsrna_id"] == locus]
        locus_rows.append({"dsrna_id": locus, "potential_siRNA_count": len(rows), "significant_both_count": sum(bool(row["significant_both"]) for row in rows), "mean_adaptation_score": mean(float(row["adaptation_score"]) for row in rows)})
    global_rows = global_permutation_rows(summary, focal_space, control_space, config.seed, config.global_permutations)
    write_rows(output / "tables" / "bowtie_adaptation_per_sirna.tsv", summary)
    write_rows(output / "tables" / "bowtie_adaptation_target_hits.tsv", hit_rows)
    write_rows(output / "tables" / "bowtie_adaptation_per_locus.tsv", locus_rows)
    write_rows(output / "tables" / "bowtie_adaptation_global_statistics.tsv", global_rows)
    return {"outdir": str(output), "sirna_count": len(summary), "focal_alignment_count": len(focal_hits), "control_alignment_count": len(control_hits), "wobble_max_pairs": config.wobble_max_pairs, "outputs": {"per_sirna": str(output / "tables" / "bowtie_adaptation_per_sirna.tsv"), "target_hits": str(output / "tables" / "bowtie_adaptation_target_hits.tsv"), "per_locus": str(output / "tables" / "bowtie_adaptation_per_locus.tsv"), "global_statistics": str(output / "tables" / "bowtie_adaptation_global_statistics.tsv")}}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sirnas", type=Path, required=True)
    parser.add_argument("--focal-transcriptome", type=Path, required=True)
    parser.add_argument("--control-transcriptome", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--input-mode", choices=("direct_sirnas", "dsrna_windows"), default="direct_sirnas")
    parser.add_argument("--shuffles", type=int, default=3)
    parser.add_argument("--mismatches", type=int, default=3)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--wobble-max-pairs", type=int, default=0, help="Expand up to this many labelled G:U/U:G wobble substitutions per guide before Bowtie mapping")
    args = parser.parse_args()
    result = run_bowtie_adaptation(BowtieAdaptationConfig(args.sirnas, args.focal_transcriptome, args.control_transcriptome, args.output_dir, input_mode=args.input_mode, shuffle_count=args.shuffles, mismatches=args.mismatches, threads=args.threads, wobble_max_pairs=args.wobble_max_pairs))
    print(f"Mapped {result['sirna_count']} original siRNAs: {result['focal_alignment_count']} focal and {result['control_alignment_count']} control Bowtie alignments.")


if __name__ == "__main__":
    main()
