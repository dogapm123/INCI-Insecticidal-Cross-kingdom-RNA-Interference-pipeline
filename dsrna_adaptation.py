"""Per-siRNA dsRNA-adaptation testing using INCI's RNAplex/MFE-ratio engine.

For every one-nucleotide-shifted guide on both dsRNA strands, this module
compares MFE-ratio target density in a focal and control transcriptome.  Each
guide is accompanied by dinucleotide-preserving shuffled controls, allowing
both a composition-controlled test and a focal-versus-control test.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import re
import shutil
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from statistics import mean
from typing import Any, Callable

from degradome_analysis import (
    DegradomeConfig,
    FastaRecord,
    build_bowtie_index,
    find_targets,
    sort_and_dedupe_targets,
    target_hits_for_query,
    write_rows,
)


LogFn = Callable[[str, str], None]
DNA = frozenset("ACGT")
IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")


@dataclass(frozen=True)
class AdaptationConfig:
    dsrna_fasta: Path
    focal_transcriptome_fasta: Path
    control_transcriptome_fasta: Path
    output_dir: Path
    shuffle_count: int = 3
    guide_length: int = 21
    mfe_ratio_cutoff: float = 0.70
    fdr_cutoff: float = 0.05
    ignore_query_pos1: bool = True
    threads: int = 4
    seed: int = 1
    annotation_file: Path | None = None
    mfe_ratio_sweep: tuple[float, ...] = ()
    input_mode: str = "dsrna_windows"
    bowtie_prefilter_mismatches: int | None = None


@dataclass(frozen=True)
class Guide:
    query_id: str
    dsrna_id: str
    dsrna_position: int
    strand: str
    sequence: str
    variant: str
    shuffle_index: int
    entropy: float


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def normalise_sequence(sequence: str, *, label: str, allow_ambiguous: bool = False) -> str:
    cleaned = "".join(sequence.upper().replace("U", "T").split())
    allowed = IUPAC_DNA if allow_ambiguous else DNA
    invalid = sorted(set(cleaned) - allowed)
    if not cleaned:
        raise ValueError(f"{label} is empty")
    if invalid:
        raise ValueError(f"{label} contains unsupported base(s): {''.join(invalid)}")
    return "".join(base if base in DNA else "N" for base in cleaned) if allow_ambiguous else cleaned


def read_fasta(path: Path, *, allow_ambiguous: bool = False) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, normalise_sequence("".join(chunks), label=name, allow_ambiguous=allow_ambiguous)))
                name = line[1:].split()[0]
                if not name:
                    raise ValueError(f"A FASTA header in {path} has no identifier")
                chunks = []
            elif name is None:
                raise ValueError(f"{path} is not FASTA: sequence occurs before a header")
            else:
                chunks.append(line)
    if name is not None:
        records.append((name, normalise_sequence("".join(chunks), label=name, allow_ambiguous=allow_ambiguous)))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def read_gene_annotations(path: Path | None) -> dict[str, str]:
    """Return an optional transcript-ID to gene-ID mapping."""
    if path is None:
        return {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        except csv.Error as exc:
            raise ValueError("Could not determine the annotation-file delimiter.") from exc
        reader = csv.DictReader(handle, dialect=dialect)
        if not reader.fieldnames or "gene_id" not in reader.fieldnames:
            raise ValueError("Annotation file must contain a gene_id column.")
        if "transcript_id" not in reader.fieldnames:
            return {}
        return {
            str(row.get("transcript_id", "")).strip(): str(row.get("gene_id", "")).strip()
            for row in reader
            if str(row.get("transcript_id", "")).strip() and str(row.get("gene_id", "")).strip()
        }


def dinucleotide_shuffle(sequence: str, rng: random.Random, preserve_prefix: int = 0) -> str:
    """Randomize an Eulerian traversal while preserving active-sequence dinucleotides.

    ``preserve_prefix`` leaves leading bases untouched.  This is required when
    a downstream guide-scoring method ignores (for example) the 5' guide base:
    the shuffled part must then be exactly the part that is evaluated.
    """
    if not 0 <= preserve_prefix < len(sequence):
        if preserve_prefix == len(sequence) == 0:
            return sequence
        raise ValueError("The preserved shuffle prefix must be shorter than the sequence.")
    prefix = sequence[:preserve_prefix]
    sequence = sequence[preserve_prefix:]
    if len(sequence) < 3:
        return prefix + sequence
        return sequence
    adjacency: dict[str, list[str]] = defaultdict(list)
    for base, next_base in zip(sequence, sequence[1:]):
        adjacency[base].append(next_base)
    for choices in adjacency.values():
        rng.shuffle(choices)
    stack = [sequence[0]]
    output: list[str] = []
    while stack:
        current = stack[-1]
        if adjacency[current]:
            stack.append(adjacency[current].pop())
        else:
            output.append(stack.pop())
    shuffled = "".join(reversed(output))
    if len(shuffled) != len(sequence):
        raise RuntimeError("Dinucleotide shuffling did not preserve guide length.")
    return prefix + shuffled


def normalized_shannon_entropy(sequence: str) -> float:
    """Base-4 Shannon entropy, where 0 is homopolymeric and 1 is maximal."""
    if not sequence:
        return 0.0
    length = float(len(sequence))
    counts = {base: sequence.count(base) for base in DNA}
    entropy = -sum((count / length) * math.log2(count / length) for count in counts.values() if count)
    return entropy / 2.0


def safe_name(value: str) -> str:
    cleaned = "".join(character if character.isalnum() or character in "._-" else "_" for character in value)
    return cleaned.strip("._")[:140] or "item"


def make_guides(records: list[tuple[str, str]], guide_length: int, shuffle_count: int, seed: int, preserve_shuffle_prefix: int = 0) -> list[Guide]:
    rng = random.Random(seed)
    guides: list[Guide] = []
    for dsrna_id, sequence in records:
        if len(sequence) < guide_length:
            raise ValueError(f"dsRNA record {dsrna_id!r} is shorter than the {guide_length}-nt guide length.")
        for offset in range(len(sequence) - guide_length + 1):
            forward = sequence[offset : offset + guide_length]
            for strand, guide_sequence in (("sense", forward), ("antisense", reverse_complement(forward))):
                prefix = f"{safe_name(dsrna_id)}|p{offset + 1}|{strand}"
                entropy = normalized_shannon_entropy(guide_sequence)
                guides.append(Guide(f"{prefix}|original", dsrna_id, offset + 1, strand, guide_sequence, "original", 0, entropy))
                for shuffle_index in range(1, shuffle_count + 1):
                    guides.append(
                        Guide(
                            f"{prefix}|shuffle{shuffle_index}", dsrna_id, offset + 1, strand,
                            dinucleotide_shuffle(guide_sequence, rng, preserve_shuffle_prefix), "shuffle", shuffle_index, entropy,
                        )
                    )
    return guides


def make_direct_sirna_guides(records: list[tuple[str, str]], guide_length: int | None, shuffle_count: int, seed: int, preserve_shuffle_prefix: int = 0) -> list[Guide]:
    """Use supplied siRNA records directly, retaining encoded position and strand metadata.

    Headers following ``<locus>_<position>pos_<length>len_<strand>_...`` are
    grouped and plotted as a locus.  Other headers remain valid and are placed
    on a sequential direct-siRNA locus.
    """
    rng = random.Random(seed)
    guides: list[Guide] = []
    header_pattern = re.compile(r"^(?P<locus>.+?)_(?P<position>\d+)pos_(?P<length>\d+)len_(?P<strand>sense|antisense)(?:_|$)")
    for ordinal, (record_id, sequence) in enumerate(records, start=1):
        if not 15 <= len(sequence) <= 30:
            raise ValueError(f"Direct siRNA record {record_id!r} has length {len(sequence)}; supported lengths are 15–30 nt.")
        if guide_length is not None and len(sequence) != guide_length:
            raise ValueError(
                f"Direct siRNA record {record_id!r} has length {len(sequence)}; expected the selected {guide_length}-nt guide length."
            )
        match = header_pattern.match(record_id)
        locus_id = match.group("locus") if match else "direct_siRNAs"
        position = int(match.group("position")) if match else ordinal
        strand = match.group("strand") if match else "sense"
        prefix = f"{safe_name(record_id)}|p{position}|{strand}"
        entropy = normalized_shannon_entropy(sequence)
        guides.append(Guide(f"{prefix}|original", locus_id, position, strand, sequence, "original", 0, entropy))
        for shuffle_index in range(1, shuffle_count + 1):
            guides.append(
                Guide(
                    f"{prefix}|shuffle{shuffle_index}", locus_id, position, strand,
                            dinucleotide_shuffle(sequence, rng, preserve_shuffle_prefix), "shuffle", shuffle_index, entropy,
                )
            )
    return guides


def transcriptome_search_space(records: list[tuple[str, str]], effective_guide_length: int) -> int:
    return sum(
        sum(set(sequence[offset : offset + effective_guide_length]) <= DNA for offset in range(max(0, len(sequence) - effective_guide_length + 1)))
        for _, sequence in records
    )


def bowtie_candidate_transcripts(
    query_records: list[FastaRecord],
    transcriptome_fasta: Path,
    transcriptome: dict[str, str],
    output_dir: Path,
    effective_query_sequences: dict[str, str],
    mismatches: int,
    threads: int,
    label: str,
) -> dict[str, set[str]]:
    """Return Bowtie1 candidate transcript IDs for each guide.

    ``--nofw`` retains only the reverse-complement orientation expected for a
    guide binding a transcript. RNAplex remains the authoritative MFE scorer.
    """
    if not shutil.which("bowtie-build") or not (shutil.which("bowtie-align-s") or shutil.which("bowtie")):
        raise ValueError("Bowtie1 and bowtie-build are required for the dsRNA-adaptation prefilter.")
    query_path = output_dir / "inputs" / f"bowtie_prefilter_{safe_name(label)}_queries.fasta"
    query_path.write_text(
        "".join(f">{query.name}\n{effective_query_sequences[query.name]}\n" for query in query_records),
        encoding="utf-8",
    )
    index_prefix = output_dir / "bowtie_index" / safe_name(transcriptome_fasta.stem)
    build_bowtie_index(transcriptome_fasta, index_prefix, threads)
    bowtie = shutil.which("bowtie-align-s") or shutil.which("bowtie")
    completed = subprocess.run(
        [
            bowtie, "-f", "-v", str(mismatches), "-a", "--best", "--nofw",
            "-p", str(max(1, threads)), str(index_prefix), str(query_path),
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode not in {0, 1}:
        raise ValueError(f"Bowtie1 prefilter failed for {label}: {completed.stderr.strip() or completed.stdout.strip()}")
    candidates: dict[str, set[str]] = defaultdict(set)
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) >= 3 and fields[0] in effective_query_sequences and fields[2] in transcriptome:
            candidates[fields[0]].add(fields[2])
    return candidates


def find_targets_for_candidates(
    query_records: list[FastaRecord],
    transcriptome: dict[str, str],
    candidates: dict[str, set[str]],
    config: DegradomeConfig,
    log: LogFn | None,
) -> list[Any]:
    """Use INCI's RNAplex target scorer only on Bowtie1 candidate pairs."""
    if not shutil.which("RNAplex"):
        raise ValueError("RNAplex was not found on PATH.")

    def run_query(query: FastaRecord) -> list[Any]:
        selected = {name: transcriptome[name] for name in candidates.get(query.name, set())}
        return target_hits_for_query(query, selected, config) if selected else []

    worker_count = min(max(1, int(config.threads)), len(query_records))
    if log:
        pair_count = sum(len(candidates.get(query.name, set())) for query in query_records)
        log(f"RNAplex scoring {pair_count:,} Bowtie1 candidate transcript/siRNA pair(s) across {len(query_records):,} guide/control sequence(s).", "info")
    if worker_count == 1:
        hits = [hit for query in query_records for hit in run_query(query)]
    else:
        results: list[list[Any] | None] = [None] * len(query_records)
        with ThreadPoolExecutor(max_workers=worker_count, thread_name_prefix="inci-adaptation-rnaplex") as executor:
            futures = {executor.submit(run_query, query): index for index, query in enumerate(query_records)}
            for future in as_completed(futures):
                results[futures[future]] = future.result()
        hits = [hit for result in results for hit in (result or [])]
    return sort_and_dedupe_targets(hits, config.sort_by)


def benjamini_hochberg(pvalues: list[float]) -> list[float]:
    if not pvalues:
        return []
    size = len(pvalues)
    ordered = sorted(enumerate(pvalues), key=lambda pair: pair[1])
    adjusted = [1.0] * size
    running = 1.0
    for rank, (index, pvalue) in reversed(list(enumerate(ordered, start=1))):
        running = min(running, max(0.0, min(1.0, pvalue * size / rank)))
        adjusted[index] = running
    return adjusted


def empirical_pvalue(observed: float, null_values: list[float]) -> float:
    return (1 + sum(value >= observed for value in null_values)) / (1 + len(null_values))


def target_metrics(
    hits: list[Any],
    search_space: int,
    transcript_to_gene: dict[str, str],
) -> dict[str, Any]:
    ratios = [float(hit.mfe_ratio) for hit in hits]
    transcripts = sorted({str(hit.transcript) for hit in hits})
    genes = sorted({transcript_to_gene.get(transcript, transcript) for transcript in transcripts})
    ratio_sum = sum(ratios)
    scale = 1_000_000.0 / max(1, search_space)
    return {
        "target_hit_count": len(hits),
        "target_transcript_count": len(transcripts),
        "target_gene_count": len(genes),
        "mfe_ratio_sum": ratio_sum,
        "mfe_ratio_density_per_million_sites": ratio_sum * scale,
        "target_hit_density_per_million_sites": len(hits) * scale,
        "target_transcripts": ",".join(transcripts),
        "target_genes": ",".join(genes),
        "mfe_ratios": ",".join(f"{ratio:.5f}" for ratio in ratios),
    }


def locus_cutoff_sweep_rows(
    guides: list[Guide],
    focal_hits: list[Any],
    control_hits: list[Any],
    focal_space: int,
    control_space: int,
    cutoffs: tuple[float, ...],
) -> list[dict[str, Any]]:
    """Summarize each MFE cutoff without repeating RNAplex calculations.

    RNAplex is run once at the lowest requested cutoff. Higher cutoffs simply
    filter those same interactions, making this a quick calibration table.
    """
    rows: list[dict[str, Any]] = []
    dsrna_ids = sorted({guide.dsrna_id for guide in guides})
    for cutoff in sorted(set(cutoffs)):
        def density(hits: list[Any]) -> float:
            return sum(float(hit.mfe_ratio) for hit in hits) * 1_000_000.0

        focal_by_query: dict[str, list[Any]] = defaultdict(list)
        control_by_query: dict[str, list[Any]] = defaultdict(list)
        for hit in focal_hits:
            if float(hit.mfe_ratio) >= cutoff:
                focal_by_query[str(hit.query)].append(hit)
        for hit in control_hits:
            if float(hit.mfe_ratio) >= cutoff:
                control_by_query[str(hit.query)].append(hit)
        for dsrna_id in dsrna_ids:
            original = [guide for guide in guides if guide.dsrna_id == dsrna_id and guide.variant == "original"]
            focus_observed = density([hit for guide in original for hit in focal_by_query.get(guide.query_id, [])]) / max(1, focal_space)
            control_observed = density([hit for guide in original for hit in control_by_query.get(guide.query_id, [])]) / max(1, control_space)
            observed_delta = focus_observed - control_observed
            shuffle_indexes = sorted({guide.shuffle_index for guide in guides if guide.dsrna_id == dsrna_id and guide.variant == "shuffle"})
            null_delta: list[float] = []
            for index in shuffle_indexes:
                shuffled = [guide for guide in guides if guide.dsrna_id == dsrna_id and guide.variant == "shuffle" and guide.shuffle_index == index]
                focus_value = density([hit for guide in shuffled for hit in focal_by_query.get(guide.query_id, [])]) / max(1, focal_space)
                control_value = density([hit for guide in shuffled for hit in control_by_query.get(guide.query_id, [])]) / max(1, control_space)
                null_delta.append(focus_value - control_value)
            rows.append(
                {
                "dsrna_id": dsrna_id,
                "mfe_ratio_cutoff": cutoff,
                "original_focus_density_per_million_sites": focus_observed,
                "original_control_density_per_million_sites": control_observed,
                "original_focus_minus_control": observed_delta,
                "mean_shuffled_focus_minus_control": mean(null_delta) if null_delta else 0.0,
                "whole_locus_adaptation_score": observed_delta - (mean(null_delta) if null_delta else 0.0),
                "whole_locus_empirical_pvalue": empirical_pvalue(observed_delta, null_delta),
                "focus_target_hit_count": sum(len(values) for values in focal_by_query.values()),
                "control_target_hit_count": sum(len(values) for values in control_by_query.values()),
                "shuffle_count": len(null_delta),
                }
            )
    return rows


def _target_rows(hits: list[Any], guides: dict[str, Guide], transcriptome: str, transcript_to_gene: dict[str, str]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for hit in hits:
        guide = guides[str(hit.query)]
        rows.append(
            {
                "dsrna_id": guide.dsrna_id,
                "dsrna_position": guide.dsrna_position,
                "strand": guide.strand,
                "variant": guide.variant,
                "shuffle_index": guide.shuffle_index,
                "siRNA_id": guide.query_id,
                "siRNA_sequence": guide.sequence,
                "transcriptome": transcriptome,
                "transcript_id": hit.transcript,
                "gene_id": transcript_to_gene.get(str(hit.transcript), str(hit.transcript)),
                "target_start": hit.t_start,
                "target_stop": hit.t_stop,
                "mfe_ratio": hit.mfe_ratio,
                "mfe_perfect": hit.mfe_perfect,
                "mfe_site": hit.mfe_site,
                "allen_score": hit.allen_score,
                "mismatches": hit.mismatch_count,
                "gu_wobbles": hit.gu_wobble_count,
                "bulges": hit.bulge_count,
                "match_pattern": hit.match_pattern,
            }
        )
    return rows


def plot_adaptation_landscape(rows: list[dict[str, Any]], output_dir: Path) -> list[Path]:
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    by_dsrna: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_dsrna[str(row["dsrna_id"])].append(row)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        return _write_svg_adaptation_landscapes(by_dsrna, plots_dir)
    colormap = plt.get_cmap("RdYlGn")
    for dsrna_id, locus_rows in by_dsrna.items():
        if not locus_rows:
            continue
        x_values = [int(row["dsrna_position"]) for row in locus_rows]
        y_values = [float(row["plot_adaptation_score"]) for row in locus_rows]
        colors = [colormap(float(row["shannon_entropy"])) for row in locus_rows]
        max_position = max(x_values)
        maximum = max((abs(value) for value in y_values), default=0.0)
        if maximum <= 0:
            maximum = 1.0
        figure, axis = plt.subplots(figsize=(12, 4.8), constrained_layout=True)
        axis.vlines(x_values, 0.0, y_values, colors=colors, linewidth=0.7, alpha=0.78, zorder=2)
        axis.scatter(x_values, y_values, c=colors, s=22, edgecolors="#1f2937", linewidths=0.25, zorder=3)
        axis.axhline(0.0, color="#475569", linewidth=0.8)
        axis.set_xlim(0.5, max_position + 0.5)
        axis.set_ylim(-maximum * 1.1, maximum * 1.1)
        axis.set_xlabel("dsRNA position (5' start of siRNA window, nt)")
        axis.set_ylabel("Strand-signed adaptation score\n(focal-over-shuffle minus control-over-shuffle)")
        axis.set_title(f"{dsrna_id}: per-siRNA adaptation landscape")
        axis.grid(axis="y", color="#e2e8f0", linewidth=0.7)
        colorbar = figure.colorbar(plt.cm.ScalarMappable(cmap=colormap, norm=plt.Normalize(0, 1)), ax=axis, pad=0.015)
        colorbar.set_label("Normalized Shannon entropy (low → red; high → green)")
        path = plots_dir / f"{safe_name(dsrna_id)}.adaptation_landscape.png"
        figure.savefig(path, dpi=300)
        plt.close(figure)
        paths.append(path)
    return paths


def _entropy_color(entropy: float) -> str:
    """Compact red-to-green fallback palette for environments without matplotlib."""
    value = max(0.0, min(1.0, entropy))
    red = int(198 * (1.0 - value) + 38 * value)
    green = int(38 * (1.0 - value) + 140 * value)
    blue = int(38 * (1.0 - value) + 50 * value)
    return f"#{red:02x}{green:02x}{blue:02x}"


def _write_svg_adaptation_landscapes(by_dsrna: dict[str, list[dict[str, Any]]], plots_dir: Path) -> list[Path]:
    """Write a portable vector plot when matplotlib is not installed."""
    paths: list[Path] = []
    width, height, left, right, top, bottom = 1200, 520, 85, 35, 45, 75
    plot_width, plot_height = width - left - right, height - top - bottom
    for dsrna_id, rows in by_dsrna.items():
        max_x = max(int(row["dsrna_position"]) for row in rows)
        limit = max(max(abs(float(row["plot_adaptation_score"])) for row in rows), 1.0)
        baseline = top + plot_height / 2
        points: list[str] = []
        for row in rows:
            x = left + (int(row["dsrna_position"]) - 1) * plot_width / max(1, max_x - 1)
            y = baseline - float(row["plot_adaptation_score"]) / (limit * 1.1) * (plot_height / 2)
            color = _entropy_color(float(row["shannon_entropy"]))
            points.append(f'<line x1="{x:.2f}" y1="{baseline:.2f}" x2="{x:.2f}" y2="{y:.2f}" stroke="{color}" stroke-width="1"/>')
            points.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3" fill="{color}" stroke="#1f2937" stroke-width="0.5"/>')
        svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="white"/>
<text x="{left}" y="25" font-family="Arial" font-size="18" font-weight="bold">{dsrna_id}: per-siRNA adaptation landscape</text>
<line x1="{left}" y1="{baseline:.2f}" x2="{width-right}" y2="{baseline:.2f}" stroke="#475569" stroke-width="1"/>
<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="#334155" stroke-width="1"/>
<line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="#334155" stroke-width="1"/>
<text x="{left}" y="{height-28}" font-family="Arial" font-size="13">dsRNA position (5' start of siRNA window, nt)</text>
<text x="12" y="{baseline-8:.2f}" font-family="Arial" font-size="12">+ sense</text>
<text x="12" y="{baseline+18:.2f}" font-family="Arial" font-size="12">− antisense</text>
{''.join(points)}
<text x="{width-355}" y="{height-28}" font-family="Arial" font-size="12" fill="#c62626">Low entropy</text>
<rect x="{width-265}" y="{height-39}" width="180" height="12" fill="url(#entropy)"/>
<text x="{width-78}" y="{height-28}" font-family="Arial" font-size="12" fill="#268c32">High</text>
<defs><linearGradient id="entropy"><stop offset="0%" stop-color="#c62626"/><stop offset="100%" stop-color="#268c32"/></linearGradient></defs>
</svg>'''
        path = plots_dir / f"{safe_name(dsrna_id)}.adaptation_landscape.svg"
        path.write_text(svg, encoding="utf-8")
        paths.append(path)
    return paths


def run_dsrna_adaptation(config: AdaptationConfig, log: LogFn | None = None) -> dict[str, Any]:
    if config.shuffle_count < 1:
        raise ValueError("At least one dinucleotide-shuffled control per siRNA is required.")
    if config.guide_length < 15 or config.guide_length > 30:
        raise ValueError("Guide length must be between 15 and 30 nt for MFE-ratio target prediction.")
    if not 0 < config.mfe_ratio_cutoff <= 1:
        raise ValueError("MFE-ratio cutoff must be >0 and <=1.")
    if not 0 < config.fdr_cutoff <= 1:
        raise ValueError("FDR cutoff must be >0 and <=1.")
    if config.threads < 1:
        raise ValueError("Threads must be at least 1.")
    if config.input_mode not in {"dsrna_windows", "direct_sirnas"}:
        raise ValueError("Input mode must be 'dsrna_windows' or 'direct_sirnas'.")
    if config.bowtie_prefilter_mismatches is not None and not 0 <= config.bowtie_prefilter_mismatches <= 3:
        raise ValueError("Bowtie1 prefilter mismatches must be between 0 and 3, or disabled.")

    output_dir = config.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "inputs").mkdir(exist_ok=True)
    dsrna_records = read_fasta(config.dsrna_fasta)
    focal_records = read_fasta(config.focal_transcriptome_fasta, allow_ambiguous=True)
    control_records = read_fasta(config.control_transcriptome_fasta, allow_ambiguous=True)
    focal = dict(focal_records)
    control = dict(control_records)
    annotations = read_gene_annotations(config.annotation_file)
    guides = (
        make_guides(dsrna_records, config.guide_length, config.shuffle_count, config.seed, int(config.ignore_query_pos1))
        if config.input_mode == "dsrna_windows"
        else make_direct_sirna_guides(dsrna_records, config.guide_length, config.shuffle_count, config.seed, int(config.ignore_query_pos1))
    )
    guide_by_id = {guide.query_id: guide for guide in guides}
    original_guides = [guide for guide in guides if guide.variant == "original"]
    effective_length = config.guide_length - (1 if config.ignore_query_pos1 else 0)
    focal_space = transcriptome_search_space(focal_records, effective_length)
    control_space = transcriptome_search_space(control_records, effective_length)
    if focal_space <= 0 or control_space <= 0:
        raise ValueError("A transcriptome has no sequence long enough for the selected effective guide length.")
    if log:
        log(f"Generated {len(original_guides):,} original siRNAs and {len(guides) - len(original_guides):,} dinucleotide-shuffled controls.", "info")
        log(f"Focal searchable sites: {focal_space:,}; control searchable sites: {control_space:,}.", "info")

    requested_cutoffs = tuple(sorted(set((config.mfe_ratio_cutoff, *config.mfe_ratio_sweep))))
    if any(not 0 < cutoff <= 1 for cutoff in requested_cutoffs):
        raise ValueError("Every MFE-ratio sweep threshold must be >0 and <=1.")
    engine_config = DegradomeConfig(
        srna_text="", srna_fasta=None, transcript_text="", transcript_fasta=None, samples=(), output_dir=output_dir,
        ignore_query_pos1=config.ignore_query_pos1, slice_positions=(10,), mfe_ratio_cutoff=min(requested_cutoffs),
        sort_by="mfe_ratio", threads=config.threads,
    )
    query_records = [FastaRecord(guide.query_id, guide.sequence) for guide in guides]
    effective_queries = {
        query.name: (query.sequence[1:] if config.ignore_query_pos1 and len(query.sequence) > 1 else query.sequence)
        for query in query_records
    }
    if config.bowtie_prefilter_mismatches is None:
        if log:
            log("Calculating focal-transcriptome siRNA MFE ratios with the INCI target-prediction engine.", "info")
        all_focal_hits = find_targets(query_records, focal, engine_config, log)
        if log:
            log("Calculating control-transcriptome siRNA MFE ratios with the INCI target-prediction engine.", "info")
        all_control_hits = find_targets(query_records, control, engine_config, log)
    else:
        if log:
            log(f"Bowtie1 prefilter: retaining reverse-complement matches with up to {config.bowtie_prefilter_mismatches} mismatch(es).", "info")
        focal_candidates = bowtie_candidate_transcripts(
            query_records, config.focal_transcriptome_fasta, focal, output_dir, effective_queries,
            config.bowtie_prefilter_mismatches, config.threads, "focal",
        )
        control_candidates = bowtie_candidate_transcripts(
            query_records, config.control_transcriptome_fasta, control, output_dir, effective_queries,
            config.bowtie_prefilter_mismatches, config.threads, "control",
        )
        if log:
            log("Calculating focal-transcriptome candidate MFE ratios with the INCI target-prediction engine.", "info")
        all_focal_hits = find_targets_for_candidates(query_records, focal, focal_candidates, engine_config, log)
        if log:
            log("Calculating control-transcriptome candidate MFE ratios with the INCI target-prediction engine.", "info")
        all_control_hits = find_targets_for_candidates(query_records, control, control_candidates, engine_config, log)
    sweep_rows = locus_cutoff_sweep_rows(guides, all_focal_hits, all_control_hits, focal_space, control_space, requested_cutoffs)
    focal_hits = [hit for hit in all_focal_hits if float(hit.mfe_ratio) >= config.mfe_ratio_cutoff]
    control_hits = [hit for hit in all_control_hits if float(hit.mfe_ratio) >= config.mfe_ratio_cutoff]
    hits_by_query: dict[str, dict[str, list[Any]]] = {
        "focal": defaultdict(list), "control": defaultdict(list),
    }
    for hit in focal_hits:
        hits_by_query["focal"][str(hit.query)].append(hit)
    for hit in control_hits:
        hits_by_query["control"][str(hit.query)].append(hit)

    metrics: dict[str, dict[str, dict[str, Any]]] = {"focal": {}, "control": {}}
    for guide in guides:
        metrics["focal"][guide.query_id] = target_metrics(hits_by_query["focal"].get(guide.query_id, []), focal_space, annotations)
        metrics["control"][guide.query_id] = target_metrics(hits_by_query["control"].get(guide.query_id, []), control_space, annotations)

    all_shuffle_focal = [float(metrics["focal"][guide.query_id]["mfe_ratio_density_per_million_sites"]) for guide in guides if guide.variant == "shuffle"]
    all_shuffle_delta = [
        float(metrics["focal"][guide.query_id]["mfe_ratio_density_per_million_sites"])
        - float(metrics["control"][guide.query_id]["mfe_ratio_density_per_million_sites"])
        for guide in guides if guide.variant == "shuffle"
    ]
    summary_rows: list[dict[str, Any]] = []
    for guide in original_guides:
        prefix = guide.query_id.rsplit("|", 1)[0]
        shuffled = [guide_by_id[f"{prefix}|shuffle{index}"] for index in range(1, config.shuffle_count + 1)]
        focal_metric = metrics["focal"][guide.query_id]
        control_metric = metrics["control"][guide.query_id]
        focal_score = float(focal_metric["mfe_ratio_density_per_million_sites"])
        control_score = float(control_metric["mfe_ratio_density_per_million_sites"])
        shuffle_focal = [float(metrics["focal"][item.query_id]["mfe_ratio_density_per_million_sites"]) for item in shuffled]
        shuffle_control = [float(metrics["control"][item.query_id]["mfe_ratio_density_per_million_sites"]) for item in shuffled]
        shuffle_delta = [focus - control_value for focus, control_value in zip(shuffle_focal, shuffle_control)]
        delta = focal_score - control_score
        adaptation_score = (focal_score - mean(shuffle_focal)) - (control_score - mean(shuffle_control))
        row = {
            "dsrna_id": guide.dsrna_id,
            "dsrna_position": guide.dsrna_position,
            "strand": guide.strand,
            "siRNA_id": guide.query_id,
            "siRNA_sequence": guide.sequence,
            "shannon_entropy": guide.entropy,
            "focus_mfe_ratio_density_per_million_sites": focal_score,
            "control_mfe_ratio_density_per_million_sites": control_score,
            "focus_minus_control_density": delta,
            "mean_shuffled_focus_density": mean(shuffle_focal),
            "mean_shuffled_control_density": mean(shuffle_control),
            "adaptation_score": adaptation_score,
            "plot_adaptation_score": (1.0 if guide.strand == "sense" else -1.0) * abs(adaptation_score),
            "within_siRNA_pvalue_vs_shuffled": empirical_pvalue(focal_score, shuffle_focal),
            "within_siRNA_pvalue_focus_vs_control": empirical_pvalue(delta, shuffle_delta),
            "pooled_pvalue_vs_shuffled": empirical_pvalue(focal_score, all_shuffle_focal),
            "pooled_pvalue_focus_vs_control": empirical_pvalue(delta, all_shuffle_delta),
            "focus_target_hit_count": focal_metric["target_hit_count"],
            "focus_target_transcript_count": focal_metric["target_transcript_count"],
            "focus_target_gene_count": focal_metric["target_gene_count"],
            "focus_target_transcripts": focal_metric["target_transcripts"],
            "focus_target_genes": focal_metric["target_genes"],
            "focus_mfe_ratios": focal_metric["mfe_ratios"],
            "control_target_hit_count": control_metric["target_hit_count"],
            "control_target_transcript_count": control_metric["target_transcript_count"],
            "control_target_gene_count": control_metric["target_gene_count"],
            "control_target_transcripts": control_metric["target_transcripts"],
            "control_target_genes": control_metric["target_genes"],
            "control_mfe_ratios": control_metric["mfe_ratios"],
            "shuffled_siRNA_sequences": ",".join(item.sequence for item in shuffled),
            "shuffled_focus_hit_counts": ",".join(str(metrics["focal"][item.query_id]["target_hit_count"]) for item in shuffled),
            "shuffled_control_hit_counts": ",".join(str(metrics["control"][item.query_id]["target_hit_count"]) for item in shuffled),
            "shuffled_focus_target_genes": ",".join(metrics["focal"][item.query_id]["target_genes"] for item in shuffled),
            "shuffled_control_target_genes": ",".join(metrics["control"][item.query_id]["target_genes"] for item in shuffled),
            "shuffled_focus_densities": ",".join(f"{value:.8g}" for value in shuffle_focal),
            "shuffled_control_densities": ",".join(f"{value:.8g}" for value in shuffle_control),
        }
        summary_rows.append(row)
    summary_rows.sort(key=lambda row: (str(row["dsrna_id"]), int(row["dsrna_position"]), str(row["strand"])))
    shuffle_pvalues = [float(row["pooled_pvalue_vs_shuffled"]) for row in summary_rows]
    specificity_pvalues = [float(row["pooled_pvalue_focus_vs_control"]) for row in summary_rows]
    for row, shuffle_fdr, specificity_fdr in zip(summary_rows, benjamini_hochberg(shuffle_pvalues), benjamini_hochberg(specificity_pvalues)):
        row["fdr_vs_shuffled"] = shuffle_fdr
        row["fdr_focus_vs_control"] = specificity_fdr
        row["significant_vs_shuffled"] = shuffle_fdr <= config.fdr_cutoff
        row["significant_focus_vs_control"] = specificity_fdr <= config.fdr_cutoff
        row["significant_both"] = bool(row["significant_vs_shuffled"] and row["significant_focus_vs_control"])

    target_rows = _target_rows(focal_hits, guide_by_id, "focal", annotations) + _target_rows(control_hits, guide_by_id, "control", annotations)
    by_locus: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in summary_rows:
        by_locus[str(row["dsrna_id"])].append(row)
    locus_rows: list[dict[str, Any]] = []
    for dsrna_id, rows in sorted(by_locus.items()):
        focus_genes = {gene for row in rows for gene in str(row["focus_target_genes"]).split(",") if gene}
        control_genes = {gene for row in rows for gene in str(row["control_target_genes"]).split(",") if gene}
        shuffled_gene_counts: list[int] = []
        for shuffle_index in range(1, config.shuffle_count + 1):
            genes = {
                gene
                for guide in guides
                if guide.dsrna_id == dsrna_id and guide.variant == "shuffle" and guide.shuffle_index == shuffle_index
                for gene in str(metrics["focal"][guide.query_id]["target_genes"]).split(",") if gene
            }
            shuffled_gene_counts.append(len(genes))
        locus_rows.append(
            {
                "dsrna_id": dsrna_id,
                "potential_siRNA_count": len(rows),
                "significant_vs_shuffled_count": sum(bool(row["significant_vs_shuffled"]) for row in rows),
                "significant_focus_vs_control_count": sum(bool(row["significant_focus_vs_control"]) for row in rows),
                "significant_both_count": sum(bool(row["significant_both"]) for row in rows),
                "focus_target_gene_count": len(focus_genes),
                "control_target_gene_count": len(control_genes),
                "mean_shuffled_focus_target_gene_count": mean(shuffled_gene_counts),
                "focus_target_genes": ",".join(sorted(focus_genes)),
                "control_target_genes": ",".join(sorted(control_genes)),
                "mean_adaptation_score": mean(float(row["adaptation_score"]) for row in rows),
            }
        )

    input_rows = [
        {"kind": "dsRNA" if config.input_mode == "dsrna_windows" else "direct_siRNA", "id": name, "length": len(sequence), "path": str(config.dsrna_fasta)} for name, sequence in dsrna_records
    ] + [
        {"kind": "focal_transcript", "id": name, "length": len(sequence), "path": str(config.focal_transcriptome_fasta)} for name, sequence in focal_records
    ] + [
        {"kind": "control_transcript", "id": name, "length": len(sequence), "path": str(config.control_transcriptome_fasta)} for name, sequence in control_records
    ]
    write_rows(output_dir / "tables" / "dsrna_adaptation_per_sirna.tsv", summary_rows)
    write_rows(output_dir / "tables" / "dsrna_adaptation_per_locus.tsv", locus_rows)
    write_rows(output_dir / "tables" / "dsrna_adaptation_mfe_target_hits.tsv", target_rows)
    write_rows(output_dir / "tables" / "dsrna_adaptation_mfe_cutoff_sweep.tsv", sweep_rows)
    write_rows(output_dir / "tables" / "dsrna_adaptation_inputs.tsv", input_rows)
    plots = plot_adaptation_landscape(summary_rows, output_dir)
    manifest = {
        "guide_length": config.guide_length, "shuffle_count_per_siRNA": config.shuffle_count,
        "mfe_ratio_cutoff": config.mfe_ratio_cutoff, "fdr_cutoff": config.fdr_cutoff,
        "ignore_query_pos1": config.ignore_query_pos1, "focal_searchable_sites": focal_space,
        "control_searchable_sites": control_space, "original_siRNA_count": len(original_guides),
        "focal_retained_mfe_target_hits": len(focal_hits), "control_retained_mfe_target_hits": len(control_hits),
        "mfe_ratio_sweep": list(requested_cutoffs),
        "input_mode": config.input_mode,
        "bowtie_prefilter_mismatches": config.bowtie_prefilter_mismatches,
    }
    (output_dir / "dsrna_adaptation_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    result = {
        "outdir": str(output_dir), "locus_count": len(locus_rows), "sirna_count": len(summary_rows),
        "significant_both_count": sum(bool(row["significant_both"]) for row in summary_rows),
        "outputs": {
            "per_sirna": str(output_dir / "tables" / "dsrna_adaptation_per_sirna.tsv"),
            "per_locus": str(output_dir / "tables" / "dsrna_adaptation_per_locus.tsv"),
            "mfe_target_hits": str(output_dir / "tables" / "dsrna_adaptation_mfe_target_hits.tsv"),
            "mfe_cutoff_sweep": str(output_dir / "tables" / "dsrna_adaptation_mfe_cutoff_sweep.tsv"),
            "plots": str(output_dir / "plots"), "manifest": str(output_dir / "dsrna_adaptation_manifest.json"),
        },
        **manifest,
    }
    if log:
        log(f"Wrote adaptation tables and {len(plots):,} locus plot(s) to {output_dir}.", "info")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dsrna", required=True, type=Path)
    parser.add_argument("--focal-transcriptome", required=True, type=Path)
    parser.add_argument("--control-transcriptome", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--annotations", type=Path)
    parser.add_argument("--shuffles", type=int, default=3)
    parser.add_argument("--guide-length", type=int, default=21)
    parser.add_argument("--mfe-ratio-cutoff", type=float, default=0.70)
    parser.add_argument("--mfe-ratio-sweep", default="", help="Comma-separated MFE-ratio cutoffs for a quick whole-locus calibration table")
    parser.add_argument("--fdr-cutoff", type=float, default=0.05)
    parser.add_argument("--include-query-pos1", action="store_true")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--input-mode", choices=("dsrna_windows", "direct_sirnas"), default="dsrna_windows")
    parser.add_argument("--bowtie-prefilter-mismatches", type=int, default=None, help="Optional Bowtie1 candidate-screen mismatch limit (0–3)")
    parser.add_argument("--disable-bowtie-prefilter", action="store_true", help="Run RNAplex against every transcript (slow, exhaustive)")
    args = parser.parse_args()
    sweep = tuple(float(value) for value in args.mfe_ratio_sweep.split(",") if value.strip())
    result = run_dsrna_adaptation(
        AdaptationConfig(
            dsrna_fasta=args.dsrna, focal_transcriptome_fasta=args.focal_transcriptome,
            control_transcriptome_fasta=args.control_transcriptome, output_dir=args.output_dir,
            annotation_file=args.annotations, shuffle_count=args.shuffles, guide_length=args.guide_length,
            mfe_ratio_cutoff=args.mfe_ratio_cutoff, fdr_cutoff=args.fdr_cutoff,
            ignore_query_pos1=not args.include_query_pos1, threads=args.threads, seed=args.seed,
            mfe_ratio_sweep=sweep,
            input_mode=args.input_mode,
            bowtie_prefilter_mismatches=None if args.disable_bowtie_prefilter else args.bowtie_prefilter_mismatches,
        ),
        lambda message, _level: print(message),
    )
    print(f"Processed {result['sirna_count']} original siRNAs across {result['locus_count']} dsRNA loci.")


if __name__ == "__main__":
    main()
