#!/usr/bin/env python3
"""Extended CleaveLand-style degradome analysis for INCI.

This module keeps the CleaveLand4 model of degradome 5' evidence at an AGO
slice site, but adds INCI-specific behavior:

* target pairing can ignore the first 5' nucleotide of each sRNA for RNAplex
  MFE and Allen-like scoring while preserving original sRNA coordinates;
* degradome tracks are CPM-normalized per sample;
* multiple samples and biological replicates are summarized as replicate lines
  or mean +/- SEM plots;
* optional exploratory q9/q10/q11 slice-site selection is supported.
"""

from __future__ import annotations

import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
import html
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence

import MFE_ratio as mfe
from plot_style import apply_matplotlib_style


LogFn = Callable[[str, str], None]
os.environ.setdefault("MPLCONFIGDIR", str(Path(os.environ.get("TMPDIR", "/tmp")) / "matplotlib-inci"))


@dataclass(frozen=True)
class FastaRecord:
    name: str
    sequence: str


@dataclass(frozen=True)
class DegradomeSample:
    sample_id: str
    group: str
    replicate: str
    path: Path


@dataclass(frozen=True)
class DegradomeConfig:
    srna_text: str
    srna_fasta: Path | None
    transcript_text: str
    transcript_fasta: Path | None
    samples: tuple[DegradomeSample, ...]
    output_dir: Path
    ignore_query_pos1: bool = True
    slice_positions: tuple[int, ...] = (10,)
    mfe_ratio_cutoff: float = 0.70
    sort_by: str = "mfe_ratio"
    pvalue_method: str = "transcript_peak_empirical"
    pvalue_cutoff: float = 0.05
    max_transcript_plots: int = 80
    compact_srna_markers: bool = True
    marker_neighborhood_nt: int = 3
    threads: int = 4


@dataclass
class TargetHit:
    query: str
    transcript: str
    t_start: int
    t_stop: int
    mfe_perfect: float
    mfe_site: float
    mfe_ratio: float
    allen_score: float
    paired: str
    unpaired: str
    structure: str
    sequence: str
    match_pattern: str
    pair_count: int
    gu_wobble_count: int
    mismatch_count: int
    bulge_count: int
    original_query_length: int
    effective_query_length: int
    ignored_query_pos1: bool
    srna_cpm: float
    slice_sites: dict[int, int]
    rank: int = 0


def normalize_pvalue_method(method: str) -> str:
    if method in {"transcript_peak_empirical", "peak_empirical", "peak_empirical_pvalue"}:
        return "transcript_peak_empirical"
    if method in {"cleaveland_category_rank", "cleaveland"}:
        return "cleaveland_category_rank"
    return "transcript_peak_empirical"


def normalize_seq(sequence: str) -> str:
    return re.sub(r"[^A-Za-z]", "", sequence).upper().replace("T", "U")


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r", encoding="utf-8")


def read_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    name = ""
    chunks: list[str] = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    records.append(FastaRecord(name, normalize_seq("".join(chunks))))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name:
        records.append(FastaRecord(name, normalize_seq("".join(chunks))))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def read_fasta_text(text: str, fallback_name: str) -> list[FastaRecord]:
    text = text.strip()
    if not text:
        return []
    temp = Path(tempfile.mkdtemp(prefix="inci_degradome_text_")) / "input.fasta"
    if text.startswith(">"):
        temp.write_text(text + "\n", encoding="utf-8")
    else:
        temp.write_text(f">{fallback_name}\n{normalize_seq(text)}\n", encoding="utf-8")
    return read_fasta(temp)


def write_fasta(records: Sequence[FastaRecord], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            seq = record.sequence.upper().replace("U", "T")
            handle.write(f">{record.name}\n")
            for idx in range(0, len(seq), 80):
                handle.write(seq[idx : idx + 80] + "\n")


def sequence_records_from_file(path: Path) -> Iterable[tuple[str, str]]:
    first = ""
    with open_text(path) as handle:
        while not first:
            first = handle.readline()
            if not first:
                return
            first = first.strip()
    if first.startswith(">"):
        name = first[1:].split()[0]
        chunks: list[str] = []
        with open_text(path) as handle:
            for raw in handle:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if chunks:
                        yield name, normalize_seq("".join(chunks)).replace("U", "T")
                    name = line[1:].split()[0]
                    chunks = []
                else:
                    chunks.append(line)
            if chunks:
                yield name, normalize_seq("".join(chunks)).replace("U", "T")
        return
    with open_text(path) as handle:
        while True:
            name = handle.readline()
            if not name:
                return
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                raise ValueError(f"Truncated FASTQ record in {path}")
            if not name.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"Expected FASTA or FASTQ records in {path}")
            yield name[1:].strip().split()[0], normalize_seq(seq).replace("U", "T")


def cpm_from_srna_id(name: str) -> float:
    match = re.search(r"_([0-9]+(?:p[0-9]+)?(?:\.[0-9]+)?)CPM(?:$|_)", name)
    if not match:
        return 0.0
    return float(match.group(1).replace("p", "."))


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")[:140] or "item"


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def group_mean_sem(values: Sequence[float]) -> tuple[float, float]:
    if not values:
        return 0.0, 0.0
    mean = sum(values) / len(values)
    if len(values) == 1:
        return mean, 0.0
    variance = sum((value - mean) ** 2 for value in values) / (len(values) - 1)
    return mean, math.sqrt(variance) / math.sqrt(len(values))


def compute_slice_site_at_query_position(
    structure: str, local_pos: Sequence[int], effective_query_position: int
) -> int:
    t_struct, q_struct = structure.split("&")
    t_chars = list(t_struct)
    real_t_pos = local_pos[1] + 1
    real_q_pos = 0
    for qch in q_struct:
        tch = t_chars.pop()
        if tch != "-":
            real_t_pos -= 1
        if qch != "-":
            real_q_pos += 1
        if real_q_pos == effective_query_position:
            return real_t_pos
    return real_t_pos


def allen_score_original_coordinates(structure: str, sequence: str, original_offset: int) -> float:
    score = 0.0
    t_brax = mfe.get_al_array(structure, 0)
    q_brax = mfe.get_al_array(structure, 1)
    t_seq = mfe.get_al_array(sequence, 0)
    q_seq = mfe.get_al_array(sequence, 1)
    effective_q_pos = 0
    for i, qs in enumerate(q_seq):
        if qs != "-":
            effective_q_pos += 1
        original_q_pos = effective_q_pos + original_offset
        qb = q_brax[i]
        ts = t_seq.pop()
        tb = t_brax.pop()
        doubled = 2 <= original_q_pos <= 12
        if tb == "." or qb == ".":
            score += 2 if doubled else 1
        elif tb == "(" and qb == ")" and {ts, qs} == {"G", "U"}:
            score += 1.0 if doubled else 0.5
    return score


def target_hits_for_query(
    query: FastaRecord,
    transcripts: dict[str, str],
    config: DegradomeConfig,
) -> list[TargetHit]:
    original = normalize_seq(query.sequence)
    if config.ignore_query_pos1 and len(original) > 1:
        effective = original[1:]
        original_offset = 1
    else:
        effective = original
        original_offset = 0
    if len(effective) < 14 or len(effective) > 26:
        return []
    perfect_mfe = mfe.get_perfect_mfe(effective)
    option_e = round(config.mfe_ratio_cutoff * perfect_mfe, 2) if perfect_mfe else None
    int_length = len(effective) + 10
    input_text = "".join(
        f">{transcript_name}\n{transcript_seq}\n>query\n{effective}\n"
        for transcript_name, transcript_seq in transcripts.items()
    )
    args = ["RNAplex", "-f", "2", "-z", str(int_length)]
    if option_e is not None:
        args.extend(["-e", f"{option_e:.2f}"])
    output = mfe.run_rnaplex(args, input_text)

    hits: list[TargetHit] = []
    current_transcript = ""
    for line in output.splitlines():
        if line.startswith(">"):
            header = line[1:].strip().split()[0]
            if header != "query":
                current_transcript = header
            continue
        if not current_transcript or not re.match(r"^[.(]+&", line):
            continue
        parsed = mfe.parse_rnaplex_line(line)
        if parsed is None:
            continue
        plex_brax, local_pos, site_mfe = parsed
        mfe_ratio = min(site_mfe / perfect_mfe, 1.0) if perfect_mfe else 0.0
        if mfe_ratio < config.mfe_ratio_cutoff:
            continue
        transcript_seq = transcripts[current_transcript]
        padded, adjusted_pos = mfe.pad_plex_brax(plex_brax, local_pos, effective)
        if adjusted_pos[0] < 1 or adjusted_pos[1] > len(transcript_seq):
            continue
        tx_site_seq = transcript_seq[adjusted_pos[0] - 1 : adjusted_pos[1]]
        if not re.match(r"^[AUCG]+$", tx_site_seq):
            continue
        trimmed_brax, trimmed_seq, adjusted_pos = mfe.no_trailing(
            padded, f"{tx_site_seq}&{effective}", adjusted_pos
        )
        structure, sequence = mfe.gapify(trimmed_brax, trimmed_seq)
        if not mfe.quality_control(structure, sequence):
            continue
        slice_sites: dict[int, int] = {}
        for original_slice_position in sorted(set(config.slice_positions)):
            effective_slice_position = original_slice_position - original_offset
            if 1 <= effective_slice_position <= len(effective):
                slice_sites[original_slice_position] = compute_slice_site_at_query_position(
                    structure, adjusted_pos, effective_slice_position
                )
        if not slice_sites:
            continue
        paired, unpaired = mfe.assess_pairing(structure, adjusted_pos[0], effective)
        metrics = mfe.alignment_metrics(structure, sequence)
        hits.append(
            TargetHit(
                query=query.name,
                transcript=current_transcript,
                t_start=adjusted_pos[0],
                t_stop=adjusted_pos[1],
                mfe_perfect=round(perfect_mfe, 4),
                mfe_site=round(site_mfe, 4),
                mfe_ratio=round(mfe_ratio, 6),
                allen_score=round(allen_score_original_coordinates(structure, sequence, original_offset), 3),
                paired=paired,
                unpaired=unpaired,
                structure=structure,
                sequence=sequence,
                match_pattern=metrics["pattern"],
                pair_count=metrics["pairs"],
                gu_wobble_count=metrics["gu_wobbles"],
                mismatch_count=metrics["mismatches"],
                bulge_count=metrics["bulges"],
                original_query_length=len(original),
                effective_query_length=len(effective),
                ignored_query_pos1=bool(original_offset),
                srna_cpm=cpm_from_srna_id(query.name),
                slice_sites=slice_sites,
            )
        )
    return hits


def sort_and_dedupe_targets(hits: list[TargetHit], sort_by: str) -> list[TargetHit]:
    best: dict[tuple[str, str, tuple[tuple[int, int], ...]], TargetHit] = {}
    for hit in hits:
        key = (hit.query, hit.transcript, tuple(sorted(hit.slice_sites.items())))
        old = best.get(key)
        if old is None:
            best[key] = hit
        elif sort_by == "allen" and (hit.allen_score, -hit.mfe_ratio) < (old.allen_score, -old.mfe_ratio):
            best[key] = hit
        elif sort_by != "allen" and (-hit.mfe_ratio, hit.allen_score) < (-old.mfe_ratio, old.allen_score):
            best[key] = hit
    ordered = sorted(
        best.values(),
        key=(lambda h: (h.allen_score, -h.mfe_ratio)) if sort_by == "allen" else (lambda h: (-h.mfe_ratio, h.allen_score)),
    )
    for index, hit in enumerate(ordered, start=1):
        hit.rank = index
    return ordered


def find_targets(
    srnas: Sequence[FastaRecord],
    transcripts: dict[str, str],
    config: DegradomeConfig,
    log: LogFn | None,
) -> list[TargetHit]:
    if not shutil.which("RNAplex"):
        raise ValueError("RNAplex was not found on PATH.")
    all_hits: list[TargetHit] = []
    worker_count = min(max(1, int(config.threads)), len(srnas))
    if worker_count == 1:
        for index, srna in enumerate(srnas, start=1):
            if log:
                log(f"Target finding for sRNA {index:,}/{len(srnas):,}: {srna.name}", "info")
            all_hits.extend(target_hits_for_query(srna, transcripts, config))
    else:
        if log:
            log(f"Target finding across {len(srnas):,} sRNAs using {worker_count} RNAplex worker(s).", "info")
        with ThreadPoolExecutor(max_workers=worker_count, thread_name_prefix="inci-rnaplex") as executor:
            futures = {
                executor.submit(target_hits_for_query, srna, transcripts, config): (index, srna)
                for index, srna in enumerate(srnas)
            }
            results_by_input_order: list[list[TargetHit] | None] = [None] * len(srnas)
            for index, future in enumerate(as_completed(futures), start=1):
                source_index, srna = futures[future]
                results_by_input_order[source_index] = future.result()
                if log and (index == len(futures) or index % max(1, len(futures) // 20) == 0):
                    log(f"Target finding completed for {index:,}/{len(futures):,} sRNAs (latest: {srna.name}).", "info")
        for result in results_by_input_order:
            all_hits.extend(result or [])
    return sort_and_dedupe_targets(all_hits, config.sort_by)


def build_bowtie_index(reference: Path, index_prefix: Path, threads: int) -> None:
    expected = [Path(f"{index_prefix}.{suffix}.ebwt") for suffix in ("1", "2", "3", "4", "rev.1", "rev.2")]
    if all(path.exists() for path in expected):
        return
    index_prefix.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        ["bowtie-build", "--threads", str(threads), str(reference), str(index_prefix)],
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr or completed.stdout or "bowtie-build failed")


def write_degradome_reads_fasta(source: Path, dest: Path) -> int:
    total = 0
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", encoding="utf-8") as handle:
        for name, seq in sequence_records_from_file(source):
            if not seq:
                continue
            total += 1
            handle.write(f">{safe_name(name) or 'read'}_{total}\n{seq}\n")
    if total == 0:
        raise ValueError(f"No degradome reads found in {source}")
    return total


def parse_sam_density(sam_path: Path, transcript_lengths: dict[str, int]) -> dict[str, dict[int, int]]:
    density: dict[str, dict[int, int]] = {name: {} for name in transcript_lengths}
    with sam_path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("@"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 11 or int(fields[1]) & 4:
                continue
            transcript = fields[2]
            pos = int(fields[3])
            if transcript not in density:
                continue
            density[transcript][pos] = density[transcript].get(pos, 0) + 1
    return density


def density_categories(density: dict[str, dict[int, int]]) -> dict[str, dict[int, int]]:
    categories: dict[str, dict[int, int]] = {}
    for transcript, pos_counts in density.items():
        if not pos_counts:
            categories[transcript] = {}
            continue
        counts = list(pos_counts.values())
        mean_nonzero = sum(counts) / len(counts)
        max_count = max(counts)
        max_positions = sum(1 for count in counts if count == max_count)
        transcript_categories: dict[int, int] = {}
        for pos, count in pos_counts.items():
            if count == 1:
                category = 4
            elif count == max_count and max_positions == 1:
                category = 0
            elif count == max_count:
                category = 1
            elif count > mean_nonzero:
                category = 2
            else:
                category = 3
            transcript_categories[pos] = category
        categories[transcript] = transcript_categories
    return categories


def map_degradome_samples(
    config: DegradomeConfig,
    transcript_fasta: Path,
    transcript_lengths: dict[str, int],
    log: LogFn | None,
) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    if not shutil.which("bowtie") and not shutil.which("bowtie-align-s"):
        raise ValueError("Bowtie1 was not found on PATH.")
    if not shutil.which("bowtie-build"):
        raise ValueError("bowtie-build was not found on PATH.")
    bowtie = shutil.which("bowtie-align-s") or shutil.which("bowtie")
    outdir = config.output_dir
    index_prefix = outdir / "bowtie_index" / transcript_fasta.stem
    build_bowtie_index(transcript_fasta, index_prefix, config.threads)
    sample_data: dict[str, dict[str, Any]] = {}
    peak_rows: list[dict[str, Any]] = []
    for sample in config.samples:
        if log:
            log(f"Mapping degradome sample {sample.sample_id}", "info")
        reads_fasta = outdir / "degradome_reads" / f"{safe_name(sample.sample_id)}.fasta"
        total_reads = write_degradome_reads_fasta(sample.path, reads_fasta)
        sam_path = outdir / "sam" / f"{safe_name(sample.sample_id)}.sam"
        sam_path.parent.mkdir(parents=True, exist_ok=True)
        completed = subprocess.run(
            [
                bowtie,
                "-f",
                "-v",
                "1",
                "--best",
                "-k",
                "1",
                "--norc",
                "--no-unal",
                "-S",
                "-p",
                str(config.threads),
                str(index_prefix),
                str(reads_fasta),
                str(sam_path),
            ],
            text=True,
            capture_output=True,
            check=False,
        )
        if completed.returncode != 0:
            raise RuntimeError(completed.stderr or completed.stdout or f"Bowtie failed for {sample.sample_id}")
        density = parse_sam_density(sam_path, transcript_lengths)
        categories = density_categories(density)
        sample_data[sample.sample_id] = {
            "sample": sample,
            "total_reads": total_reads,
            "density": density,
            "categories": categories,
        }
        for transcript, pos_counts in density.items():
            for pos, count in sorted(pos_counts.items()):
                peak_rows.append(
                    {
                        "sample_id": sample.sample_id,
                        "group": sample.group,
                        "replicate": sample.replicate,
                        "transcript": transcript,
                        "position": pos,
                        "raw_5p_count": count,
                        "degradome_cpm": count * 1_000_000.0 / total_reads,
                        "category": categories.get(transcript, {}).get(pos, ""),
                    }
                )
    return sample_data, peak_rows


def binomial_pvalue_at_least_one(trials: int, probability: float) -> float:
    if trials <= 0 or probability <= 0:
        return 1.0
    if probability >= 1:
        return 1.0
    return 1.0 - ((1.0 - probability) ** trials)


def clamp_pvalue(value: float) -> float:
    if not math.isfinite(value):
        return 1.0
    return min(max(value, 1e-300), 1.0)


def transcript_mean_cpm_track(
    sample_data: dict[str, dict[str, Any]],
    transcript: str,
    length: int,
) -> list[float]:
    values_by_position: list[list[float]] = [[] for _ in range(length)]
    for data in sample_data.values():
        total_reads = float(data["total_reads"]) or 1.0
        density = data["density"].get(transcript, {})
        for pos in range(1, length + 1):
            count = int(density.get(pos, 0))
            values_by_position[pos - 1].append(count * 1_000_000.0 / total_reads)
    means: list[float] = []
    for values in values_by_position:
        means.append(sum(values) / len(values) if values else 0.0)
    return means


def empirical_peak_pvalue(
    sample_data: dict[str, dict[str, Any]],
    transcript_lengths: dict[str, int],
    transcript: str,
    observed_mean_cpm: float,
) -> float:
    length = int(transcript_lengths.get(transcript, 0))
    if length <= 0:
        return 1.0
    track = transcript_mean_cpm_track(sample_data, transcript, length)
    as_extreme = sum(1 for value in track if value >= observed_mean_cpm)
    return clamp_pvalue((1 + as_extreme) / (1 + length))


def build_cleavage_rows(
    hits: Sequence[TargetHit],
    sample_data: dict[str, dict[str, Any]],
    transcript_lengths: dict[str, int],
    config: DegradomeConfig,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    pvalue_method = normalize_pvalue_method(config.pvalue_method)
    effective_size = max(1, sum(transcript_lengths.values()) - len(transcript_lengths) * 20)
    category_totals: dict[int, int] = {}
    for data in sample_data.values():
        for transcript_categories in data["categories"].values():
            for category in transcript_categories.values():
                category_totals[category] = category_totals.get(category, 0) + 1
    for hit in hits:
        qpos_stats: dict[int, dict[str, Any]] = {}
        for qpos, tpos in hit.slice_sites.items():
            sample_cpms: list[float] = []
            sample_counts: list[int] = []
            sample_categories: list[int] = []
            for data in sample_data.values():
                count = int(data["density"].get(hit.transcript, {}).get(tpos, 0))
                sample_counts.append(count)
                sample_cpms.append(count * 1_000_000.0 / float(data["total_reads"]))
                category = data["categories"].get(hit.transcript, {}).get(tpos)
                if category is not None:
                    sample_categories.append(int(category))
            mean_cpm, sem_cpm = group_mean_sem(sample_cpms)
            qpos_stats[qpos] = {
                "tpos": tpos,
                "counts": sample_counts,
                "mean_cpm": mean_cpm,
                "sem_cpm": sem_cpm,
                "total_count": sum(sample_counts),
                "best_category": min(sample_categories) if sample_categories else "",
                "sample_categories": sample_categories,
            }
        selected_qpos, selected = max(
            qpos_stats.items(),
            key=lambda item: (item[1]["mean_cpm"], item[1]["total_count"], -abs(item[0] - 10)),
        )
        if selected["total_count"] <= 0:
            continue
        best_category = selected["best_category"]
        chance = (category_totals.get(int(best_category), 0) / effective_size) if best_category != "" else 0.0
        cleaveland_pvalue = binomial_pvalue_at_least_one(hit.rank, chance)
        peak_pvalue = empirical_peak_pvalue(sample_data, transcript_lengths, hit.transcript, float(selected["mean_cpm"]))
        selected_pvalue = (
            cleaveland_pvalue
            if pvalue_method == "cleaveland_category_rank"
            else peak_pvalue
        )
        rows.append(
            {
                "query": hit.query,
                "transcript": hit.transcript,
                "pvalue_method": pvalue_method,
                "selected_pvalue": selected_pvalue,
                "transcript_peak_empirical_pvalue": peak_pvalue,
                "cleaveland_category_rank_pvalue": cleaveland_pvalue,
                "selected_original_query_position": selected_qpos,
                "selected_slice_site": selected["tpos"],
                "canonical_q10_slice_site": hit.slice_sites.get(10, ""),
                "offset_mode": ",".join(str(v) for v in config.slice_positions),
                "mean_degradome_cpm": selected["mean_cpm"],
                "sem_degradome_cpm": selected["sem_cpm"],
                "total_degradome_5p_count": selected["total_count"],
                "best_category": best_category,
                "mfe_ratio": hit.mfe_ratio,
                "allen_score": hit.allen_score,
                "mfe_perfect": hit.mfe_perfect,
                "mfe_site": hit.mfe_site,
                "srna_cpm": hit.srna_cpm,
                "t_start": hit.t_start,
                "t_stop": hit.t_stop,
                "ignored_query_pos1": hit.ignored_query_pos1,
                "original_query_length": hit.original_query_length,
                "effective_query_length": hit.effective_query_length,
                "paired": hit.paired,
                "unpaired": hit.unpaired,
                "match_pattern": hit.match_pattern,
                "pair_count": hit.pair_count,
                "gu_wobble_count": hit.gu_wobble_count,
                "mismatch_count": hit.mismatch_count,
                "bulge_count": hit.bulge_count,
                "all_candidate_slice_sites": json.dumps(hit.slice_sites, sort_keys=True),
            }
        )
    rows.sort(key=lambda row: (str(row["transcript"]), int(row["selected_slice_site"]), str(row["query"])))
    return rows


def target_rows(hits: Sequence[TargetHit]) -> list[dict[str, Any]]:
    return [
        {
            "rank": hit.rank,
            "query": hit.query,
            "transcript": hit.transcript,
            "t_start": hit.t_start,
            "t_stop": hit.t_stop,
            "slice_sites": json.dumps(hit.slice_sites, sort_keys=True),
            "mfe_ratio": hit.mfe_ratio,
            "allen_score": hit.allen_score,
            "mfe_perfect": hit.mfe_perfect,
            "mfe_site": hit.mfe_site,
            "srna_cpm": hit.srna_cpm,
            "ignored_query_pos1": hit.ignored_query_pos1,
            "original_query_length": hit.original_query_length,
            "effective_query_length": hit.effective_query_length,
            "paired": hit.paired,
            "unpaired": hit.unpaired,
            "structure": hit.structure,
            "sequence": hit.sequence,
            "match_pattern": hit.match_pattern,
            "pair_count": hit.pair_count,
            "gu_wobble_count": hit.gu_wobble_count,
            "mismatch_count": hit.mismatch_count,
            "bulge_count": hit.bulge_count,
        }
        for hit in hits
    ]


def transcript_summary_rows(
    transcripts: dict[str, str],
    hits: Sequence[TargetHit],
    cleavage_rows: Sequence[dict[str, Any]],
    pvalue_cutoff: float,
) -> list[dict[str, Any]]:
    hits_by_tx: dict[str, list[TargetHit]] = {}
    cleavage_by_tx: dict[str, list[dict[str, Any]]] = {}
    for hit in hits:
        hits_by_tx.setdefault(hit.transcript, []).append(hit)
    for row in cleavage_rows:
        cleavage_by_tx.setdefault(str(row["transcript"]), []).append(row)
    rows: list[dict[str, Any]] = []
    for transcript, seq in transcripts.items():
        tx_hits = hits_by_tx.get(transcript, [])
        tx_cleavage = cleavage_by_tx.get(transcript, [])
        significant_cleavage = [
            row
            for row in tx_cleavage
            if float(row.get("selected_pvalue") or 1.0) < pvalue_cutoff
        ]
        selected_positions = [int(row["selected_original_query_position"]) for row in tx_cleavage]
        rows.append(
            {
                "transcript": transcript,
                "length": len(seq),
                "predicted_srna_count": len({hit.query for hit in tx_hits}),
                "predicted_target_site_count": len(tx_hits),
                "cleavage_supported_srna_count": len({str(row["query"]) for row in tx_cleavage}),
                "cleavage_supported_site_count": len(tx_cleavage),
                "significant_srna_count": len({str(row["query"]) for row in significant_cleavage}),
                "significant_cleavage_site_count": len(significant_cleavage),
                "canonical_q10_supported_sites": sum(1 for value in selected_positions if value == 10),
                "offset_q9_supported_sites": sum(1 for value in selected_positions if value == 9),
                "offset_q11_supported_sites": sum(1 for value in selected_positions if value == 11),
                "best_mean_degradome_cpm": max((float(row["mean_degradome_cpm"]) for row in tx_cleavage), default=0.0),
                "best_selected_pvalue": min((float(row["selected_pvalue"]) for row in tx_cleavage), default=1.0),
                "best_mfe_ratio": max((hit.mfe_ratio for hit in tx_hits), default=0.0),
            }
        )
    return rows


def srna_summary_rows(
    srnas: Sequence[FastaRecord],
    hits: Sequence[TargetHit],
    cleavage_rows: Sequence[dict[str, Any]],
    pvalue_cutoff: float,
) -> list[dict[str, Any]]:
    hits_by_srna: dict[str, list[TargetHit]] = {}
    cleavage_by_srna: dict[str, list[dict[str, Any]]] = {}
    for hit in hits:
        hits_by_srna.setdefault(hit.query, []).append(hit)
    for row in cleavage_rows:
        cleavage_by_srna.setdefault(str(row["query"]), []).append(row)

    rows: list[dict[str, Any]] = []
    for srna in srnas:
        srna_hits = hits_by_srna.get(srna.name, [])
        srna_cleavage = cleavage_by_srna.get(srna.name, [])
        significant_cleavage = [
            row
            for row in srna_cleavage
            if float(row.get("selected_pvalue") or 1.0) < pvalue_cutoff
        ]
        rows.append(
            {
                "srna": srna.name,
                "length": len(srna.sequence),
                "predicted_transcript_count": len({hit.transcript for hit in srna_hits}),
                "predicted_target_site_count": len(srna_hits),
                "cleavage_supported_transcript_count": len({str(row["transcript"]) for row in srna_cleavage}),
                "cleavage_supported_site_count": len(srna_cleavage),
                "significant_transcript_count": len({str(row["transcript"]) for row in significant_cleavage}),
                "significant_cleavage_site_count": len(significant_cleavage),
                "best_selected_pvalue": min((float(row["selected_pvalue"]) for row in srna_cleavage), default=1.0),
                "best_mfe_ratio": max((hit.mfe_ratio for hit in srna_hits), default=0.0),
            }
        )
    return rows


def degradome_overall_summary_rows(
    transcript_summary: Sequence[dict[str, Any]],
    srna_summary: Sequence[dict[str, Any]],
    pvalue_method: str,
    pvalue_cutoff: float,
) -> list[dict[str, Any]]:
    predicted_srna_per_tx = [int(row["predicted_srna_count"]) for row in transcript_summary]
    supported_srna_per_tx = [int(row["cleavage_supported_srna_count"]) for row in transcript_summary]
    significant_srna_per_tx = [int(row["significant_srna_count"]) for row in transcript_summary]
    predicted_tx_per_srna = [int(row["predicted_transcript_count"]) for row in srna_summary]
    supported_tx_per_srna = [int(row["cleavage_supported_transcript_count"]) for row in srna_summary]
    significant_tx_per_srna = [int(row["significant_transcript_count"]) for row in srna_summary]

    def mean(values: Sequence[int]) -> float:
        return sum(values) / len(values) if values else 0.0

    return [
        {"section": "settings", "metric": "pvalue_method", "value": pvalue_method},
        {"section": "settings", "metric": "pvalue_cutoff", "value": pvalue_cutoff},
        {"section": "transcript", "metric": "total_transcripts_in_input", "value": len(transcript_summary)},
        {
            "section": "transcript",
            "metric": "transcripts_with_at_least_one_predicted_srna",
            "value": sum(1 for row in transcript_summary if int(row["predicted_srna_count"]) > 0),
        },
        {"section": "transcript", "metric": "mean_sRNAs_per_transcript_predicted", "value": mean(predicted_srna_per_tx)},
        {
            "section": "transcript",
            "metric": "transcripts_with_at_least_one_degradome_supported_srna",
            "value": sum(1 for row in transcript_summary if int(row["cleavage_supported_srna_count"]) > 0),
        },
        {"section": "transcript", "metric": "mean_sRNAs_per_transcript_degradome_supported", "value": mean(supported_srna_per_tx)},
        {
            "section": "transcript",
            "metric": "transcripts_with_at_least_one_significant_srna",
            "value": sum(1 for row in transcript_summary if int(row["significant_srna_count"]) > 0),
        },
        {"section": "transcript", "metric": "mean_sRNAs_per_transcript_significant", "value": mean(significant_srna_per_tx)},
        {"section": "srna", "metric": "total_sRNAs_in_input", "value": len(srna_summary)},
        {
            "section": "srna",
            "metric": "sRNAs_predicted_to_target_at_least_one_transcript",
            "value": sum(1 for row in srna_summary if int(row["predicted_transcript_count"]) > 0),
        },
        {"section": "srna", "metric": "mean_transcripts_per_sRNA_predicted", "value": mean(predicted_tx_per_srna)},
        {
            "section": "srna",
            "metric": "sRNAs_with_degradome_support_on_at_least_one_transcript",
            "value": sum(1 for row in srna_summary if int(row["cleavage_supported_transcript_count"]) > 0),
        },
        {"section": "srna", "metric": "mean_transcripts_per_sRNA_degradome_supported", "value": mean(supported_tx_per_srna)},
        {
            "section": "srna",
            "metric": "sRNAs_significantly_cleaving_at_least_one_transcript",
            "value": sum(1 for row in srna_summary if int(row["significant_transcript_count"]) > 0),
        },
        {"section": "srna", "metric": "mean_transcripts_per_sRNA_significant", "value": mean(significant_tx_per_srna)},
    ]


def sample_track(
    data: dict[str, Any], transcript: str, length: int
) -> list[float]:
    total_reads = float(data["total_reads"]) or 1.0
    arr = [0.0] * length
    for pos, count in data["density"].get(transcript, {}).items():
        if 1 <= pos <= length:
            arr[pos - 1] = count * 1_000_000.0 / total_reads
    return arr


def noncolliding_label_offset(
    x_value: float,
    y_value: float,
    occupied: list[tuple[float, float]],
    x_span: float,
    y_span: float,
) -> tuple[float, float]:
    """Return an annotation offset that avoids nearby existing labels."""

    base_dx = 6.0
    base_dy = 6.0
    if not occupied:
        occupied.append((x_value, y_value))
        return base_dx, base_dy
    x_threshold = max(x_span * 0.018, 1.0)
    y_threshold = max(y_span * 0.055, 1.0)
    lane = 0
    while True:
        candidate_y = y_value + lane * y_threshold
        collides = any(
            abs(x_value - old_x) <= x_threshold and abs(candidate_y - old_y) <= y_threshold
            for old_x, old_y in occupied
        )
        if not collides:
            occupied.append((x_value, candidate_y))
            return base_dx, base_dy + lane * 11.0
        lane += 1
        if lane > 12:
            occupied.append((x_value, candidate_y))
            return base_dx, base_dy + lane * 11.0


def select_transcript_plot_markers(
    rows: Sequence[dict[str, Any]],
    config: DegradomeConfig,
) -> list[dict[str, Any]]:
    """Keep the most informative nearby cleavage markers for a readable plot."""

    if not config.compact_srna_markers:
        return list(rows)

    best_at_exact_site: dict[int, dict[str, Any]] = {}
    for row in rows:
        site = int(row["selected_slice_site"])
        previous = best_at_exact_site.get(site)
        if previous is None or (
            float(row["mfe_ratio"]),
            float(row.get("mean_degradome_cpm") or 0.0),
            float(row.get("srna_cpm") or 0.0),
            str(row["query"]),
        ) > (
            float(previous["mfe_ratio"]),
            float(previous.get("mean_degradome_cpm") or 0.0),
            float(previous.get("srna_cpm") or 0.0),
            str(previous["query"]),
        ):
            best_at_exact_site[site] = row

    candidate_rows = sorted(
        best_at_exact_site.values(),
        key=lambda row: (
            -float(row.get("mean_degradome_cpm") or 0.0),
            -float(row["mfe_ratio"]),
            -float(row.get("total_degradome_5p_count") or 0.0),
            str(row["query"]),
        ),
    )
    selected: list[dict[str, Any]] = []
    neighborhood = max(1, int(config.marker_neighborhood_nt))
    for row in candidate_rows:
        site = int(row["selected_slice_site"])
        if any(abs(site - int(kept["selected_slice_site"])) <= neighborhood for kept in selected):
            continue
        selected.append(row)
    return sorted(selected, key=lambda row: (int(row["selected_slice_site"]), str(row["query"])))


def plot_transcript(
    transcript: str,
    length: int,
    sample_data: dict[str, dict[str, Any]],
    cleavage_rows: Sequence[dict[str, Any]],
    config: DegradomeConfig,
    path: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    style = apply_matplotlib_style()
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors as mpl_colors
    import numpy as np

    x = np.arange(1, length + 1)
    fig, ax = plt.subplots(figsize=(style["figure_width"], style["figure_height"]), constrained_layout=True)
    group_colors = ["#0f766e", "#7c3aed", "#ea580c", "#2563eb", "#be123c", "#16a34a"]
    groups = sorted({data["sample"].group for data in sample_data.values()})
    y_max = 0.0
    plotted_means: list[Any] = []
    for group_index, group in enumerate(groups):
        group_samples = [data for data in sample_data.values() if data["sample"].group == group]
        arrays = np.array([sample_track(data, transcript, length) for data in group_samples], dtype=float)
        if arrays.size == 0:
            continue
        color = group_colors[group_index % len(group_colors)]
        mean = np.mean(arrays, axis=0)
        sem = np.std(arrays, axis=0, ddof=1) / np.sqrt(arrays.shape[0]) if arrays.shape[0] > 1 else np.zeros_like(mean)
        nz = mean > 0
        if not np.any(nz):
            continue
        plotted_means.append(mean)
        upper = mean + sem
        lower = np.maximum(mean - sem, 0)
        y_max = max(y_max, float(np.max(upper[nz])))
        ax.plot(
            x,
            mean,
            color=color,
            linewidth=max(0.34, float(style["line_width"]) * 0.48),
            label=f"{group} mean",
            zorder=3,
        )
        if np.any(sem > 0):
            sem_line_width = max(0.22, float(style["line_width"]) * 0.26)
            ax.plot(
                x,
                upper,
                color=color,
                linewidth=sem_line_width,
                linestyle=(0, (4, 3)),
                alpha=0.72,
                label="_nolegend_",
                zorder=2,
            )
            ax.plot(
                x,
                lower,
                color=color,
                linewidth=sem_line_width,
                linestyle=(0, (4, 3)),
                alpha=0.72,
                label="_nolegend_",
                zorder=2,
            )
    ax.set_xlabel("Transcript position (nt)")
    ax.set_ylabel("Degradome 5' end CPM")
    ax.set_title(f"{transcript} degradome cleavage support")
    ax.grid(axis="y", color="#e5e7eb", linewidth=style["grid_width"], alpha=0.8)

    tx_rows = [row for row in cleavage_rows if row["transcript"] == transcript]
    marker_rows = select_transcript_plot_markers(tx_rows, config)
    srna_cpms = [float(row.get("srna_cpm") or 0.0) for row in marker_rows]
    max_srna_cpm = max(srna_cpms + [1.0])
    cpm_norm = mpl_colors.Normalize(vmin=0.0, vmax=max_srna_cpm)
    cpm_cmap = cm.get_cmap("viridis")
    marker_lift = max(y_max * 0.03, 1.0)
    stacked_marker_lanes: list[float] = []
    label_positions: list[tuple[float, float]] = []
    x_span = max(float(length), 1.0)
    y_span = max(y_max, 1.0)

    def local_peak_marker_position(pos: int, fallback_y: float) -> tuple[float, float]:
        search_radius = max(2, min(5, int(round(x_span * 0.002))))
        left = max(1, pos - search_radius)
        right = min(length, pos + search_radius)
        best_pos = pos
        best_y = fallback_y
        for mean_curve in plotted_means:
            for candidate_pos in range(left, right + 1):
                candidate_y = float(mean_curve[candidate_pos - 1])
                if candidate_y > best_y:
                    best_y = candidate_y
                    best_pos = candidate_pos
        return float(best_pos), best_y

    for idx, row in enumerate(marker_rows):
        pos = int(row["selected_slice_site"])
        srna_cpm = float(row.get("srna_cpm") or 0.0)
        mfe_ratio = float(row["mfe_ratio"])
        peak_y = float(row.get("mean_degradome_cpm") or 0.0)
        x_marker, peak_y = local_peak_marker_position(pos, peak_y)
        nearby_count = sum(1 for old_x in stacked_marker_lanes if abs(old_x - x_marker) <= 1.0)
        stacked_marker_lanes.append(x_marker)
        y_marker = peak_y + marker_lift * (0.12 + nearby_count * 0.34)
        color = cpm_cmap(cpm_norm(srna_cpm))
        y_max = max(y_max, y_marker)
        ax.axvline(pos, color="#8b5cf6", linewidth=max(0.18, style["line_width"] * 0.18), linestyle=":", alpha=0.38)
        ax.scatter([x_marker], [y_marker], s=float(style["marker_size"]) ** 2 * 0.42, color=[color], edgecolor="#202020", linewidth=max(0.18, style["line_width"] * 0.11), zorder=5)
        label = f"{row['query']} {mfe_ratio:.2f}"
        dx, dy = noncolliding_label_offset(float(x_marker), float(y_marker), label_positions, x_span, y_span)
        ax.annotate(
            label,
            (x_marker, y_marker),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=max(6.0, float(style["font_size"]) - 3.0),
            color="#25313d",
            arrowprops={"arrowstyle": "-", "color": "#94a3b8", "linewidth": max(0.3, style["line_width"] * 0.22), "alpha": 0.7} if dy > 18 else None,
        )
    if marker_rows:
        scalar = cm.ScalarMappable(norm=cpm_norm, cmap=cpm_cmap)
        scalar.set_array([])
        cbar = fig.colorbar(scalar, ax=ax, pad=0.012, fraction=0.035)
        cbar.set_label("sRNA CPM")
    if y_max > 0:
        ax.set_ylim(0, y_max * 1.16)
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, frameon=False, loc="upper right")
    fig.savefig(path, dpi=max(int(style["dpi"]), 360))
    plt.close(fig)


def plot_summary(
    transcript_summary: Sequence[dict[str, Any]],
    srna_summary: Sequence[dict[str, Any]],
    pvalue_cutoff: float,
    path: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    style = apply_matplotlib_style()
    import matplotlib.pyplot as plt
    import numpy as np

    if not transcript_summary and not srna_summary:
        return
    transcript_counts = {
        "Input transcripts": len(transcript_summary),
        "With predicted sRNA": sum(1 for row in transcript_summary if int(row["predicted_srna_count"]) > 0),
        "With degradome-supported sRNA": sum(
            1 for row in transcript_summary if int(row["cleavage_supported_srna_count"]) > 0
        ),
        f"With significant sRNA\nP < {pvalue_cutoff:g}": sum(
            1 for row in transcript_summary if int(row["significant_srna_count"]) > 0
        ),
    }
    srna_counts = {
        "Input sRNAs": len(srna_summary),
        "Target >=1 transcript": sum(1 for row in srna_summary if int(row["predicted_transcript_count"]) > 0),
        "With degradome support": sum(
            1 for row in srna_summary if int(row["cleavage_supported_transcript_count"]) > 0
        ),
        f"Significant cleavage\nP < {pvalue_cutoff:g}": sum(
            1 for row in srna_summary if int(row["significant_transcript_count"]) > 0
        ),
    }
    srnas_per_tx_pred = [int(row["predicted_srna_count"]) for row in transcript_summary]
    srnas_per_tx_sig = [int(row["significant_srna_count"]) for row in transcript_summary]
    tx_per_srna_pred = [int(row["predicted_transcript_count"]) for row in srna_summary]
    tx_per_srna_sig = [int(row["significant_transcript_count"]) for row in srna_summary]

    def mean_count(values: Sequence[int]) -> float:
        return sum(values) / len(values) if values else 0.0

    def format_count(value: float) -> str:
        if abs(value - round(value)) < 1e-9:
            return f"{int(round(value)):,}"
        return f"{value:.2f}".rstrip("0").rstrip(".")

    def plot_two_category_bar(
        ax: Any,
        values: Sequence[float],
        ylabel: str,
        title: str,
    ) -> None:
        labels = ["Predicted", f"Significant\nP < {pvalue_cutoff:g}"]
        positions = np.arange(len(labels))
        bars = ax.bar(positions, values, color=[count_color, sig_color], width=0.58)
        ax.set_xticks(positions, labels)
        ax.tick_params(axis="x", length=0)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_ylim(0, max([*values, 1.0]) * 1.24)
        for bar, value in zip(bars, values):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                value,
                format_count(float(value)),
                ha="center",
                va="bottom",
                fontsize=style["font_size"] * 0.95,
            )

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(max(style["figure_width"], 12.0), max(style["figure_height"] * 1.55, 7.2)),
        constrained_layout=True,
    )
    count_color = "#2563eb"
    support_color = "#7c3aed"
    sig_color = "#0f766e"
    total_color = "#64748b"

    labels = list(transcript_counts)
    values = list(transcript_counts.values())
    axes[0, 0].barh(np.arange(len(labels)), values, color=[total_color, count_color, support_color, sig_color])
    axes[0, 0].set_yticks(np.arange(len(labels)), labels)
    axes[0, 0].invert_yaxis()
    axes[0, 0].set_xlabel("Transcript count")
    axes[0, 0].set_title("Transcript-side support")
    axes[0, 0].set_xlim(0, max(values + [1]) * 1.18)
    for idx, value in enumerate(values):
        axes[0, 0].text(value, idx, f" {value:,}", va="center", fontsize=style["font_size"] * 0.9)

    labels = list(srna_counts)
    values = list(srna_counts.values())
    axes[1, 0].barh(np.arange(len(labels)), values, color=[total_color, count_color, support_color, sig_color])
    axes[1, 0].set_yticks(np.arange(len(labels)), labels)
    axes[1, 0].invert_yaxis()
    axes[1, 0].set_xlabel("sRNA count")
    axes[1, 0].set_title("sRNA-side support")
    axes[1, 0].set_xlim(0, max(values + [1]) * 1.18)
    for idx, value in enumerate(values):
        axes[1, 0].text(value, idx, f" {value:,}", va="center", fontsize=style["font_size"] * 0.9)

    plot_two_category_bar(
        axes[0, 1],
        [mean_count(srnas_per_tx_pred), mean_count(srnas_per_tx_sig)],
        "Mean sRNAs per transcript",
        "sRNAs per Transcript",
    )

    plot_two_category_bar(
        axes[1, 1],
        [mean_count(tx_per_srna_pred), mean_count(tx_per_srna_sig)],
        "Mean transcripts per sRNA",
        "Transcripts per sRNA",
    )

    for sub_ax in np.ravel(axes):
        grid_axis = "x" if sub_ax is axes[0, 0] or sub_ax is axes[1, 0] else "y"
        sub_ax.grid(axis=grid_axis, color="#e5e7eb", linewidth=style["grid_width"])
    fig.savefig(path, dpi=max(int(style["dpi"]), 360))
    plt.close(fig)


def write_report(path: Path, result: dict[str, Any]) -> None:
    links = []
    for label, file_path in result["outputs"].items():
        links.append(f"<li><strong>{html.escape(label)}</strong>: {html.escape(str(file_path))}</li>")
    body = f"""<!doctype html>
<html>
<head><meta charset="utf-8"><title>Degradome Analysis Report</title></head>
<body>
<h1>Degradome Analysis Report</h1>
<p>Extended CleaveLand-style run with q1-ignored target pairing enabled where selected.</p>
<ul>
<li>sRNAs: {result['srna_count']}</li>
<li>Transcripts: {result['transcript_count']}</li>
<li>Predicted target sites: {result['target_count']}</li>
<li>Cleavage-supported sites: {result['cleavage_count']}</li>
<li>Significant cleavage-supported sites: {result.get('significant_cleavage_count', 0)}</li>
<li>P-value method: {html.escape(str(result.get('pvalue_method', '')))}</li>
<li>P-value cutoff: {html.escape(str(result.get('pvalue_cutoff', '')))}</li>
<li>Samples: {result['sample_count']}</li>
</ul>
<h2>Outputs</h2>
<ul>{''.join(links)}</ul>
</body>
</html>
"""
    path.write_text(body, encoding="utf-8")


def run_degradome_analysis(config: DegradomeConfig, log: LogFn | None = None) -> dict[str, Any]:
    outdir = config.output_dir
    outdir.mkdir(parents=True, exist_ok=True)
    for folder in ("tables", "plots", "inputs"):
        (outdir / folder).mkdir(exist_ok=True)

    srnas = read_fasta(config.srna_fasta) if config.srna_fasta else read_fasta_text(config.srna_text, "pasted_sRNA")
    transcripts_records = (
        read_fasta(config.transcript_fasta)
        if config.transcript_fasta
        else read_fasta_text(config.transcript_text, "pasted_transcript")
    )
    if not srnas:
        raise ValueError("Provide sRNAs by pasting sequences or selecting an sRNA FASTA.")
    if not transcripts_records:
        raise ValueError("Provide transcripts by pasting a sequence or selecting a transcript FASTA.")
    if not config.samples:
        raise ValueError("Add at least one degradome sample.")

    srna_input = outdir / "inputs" / "srnas.fasta"
    transcript_input = outdir / "inputs" / "transcripts.fasta"
    write_fasta(srnas, srna_input)
    write_fasta(transcripts_records, transcript_input)
    transcripts = {record.name: record.sequence for record in transcripts_records}
    transcript_lengths = {record.name: len(record.sequence) for record in transcripts_records}

    if log:
        log(f"Loaded {len(srnas):,} sRNA(s), {len(transcripts):,} transcript(s), and {len(config.samples):,} degradome sample(s).", "info")
    hits = find_targets(srnas, transcripts, config, log)
    sample_data, peak_rows = map_degradome_samples(config, transcript_input, transcript_lengths, log)
    cleavage = build_cleavage_rows(hits, sample_data, transcript_lengths, config)
    pvalue_method = normalize_pvalue_method(config.pvalue_method)
    summary = transcript_summary_rows(transcripts, hits, cleavage, config.pvalue_cutoff)
    srna_summary = srna_summary_rows(srnas, hits, cleavage, config.pvalue_cutoff)
    overall_summary = degradome_overall_summary_rows(summary, srna_summary, pvalue_method, config.pvalue_cutoff)

    target_path = outdir / "tables" / "target_predictions.tsv"
    cleavage_path = outdir / "degradome_cleavage_sites.tsv"
    peak_path = outdir / "tables" / "degradome_5p_density_by_sample.tsv"
    summary_path = outdir / "tables" / "transcript_summary.tsv"
    srna_summary_path = outdir / "tables" / "srna_summary.tsv"
    overall_summary_path = outdir / "tables" / "degradome_overall_summary.tsv"
    write_rows(target_path, target_rows(hits))
    write_rows(cleavage_path, cleavage)
    write_rows(peak_path, peak_rows)
    write_rows(summary_path, summary)
    write_rows(srna_summary_path, srna_summary)
    write_rows(overall_summary_path, overall_summary)

    plots: list[Path] = []
    summary_plot = outdir / "plots" / "degradome_summary.png"
    plot_summary(summary, srna_summary, config.pvalue_cutoff, summary_plot)
    if summary_plot.exists():
        plots.append(summary_plot)
    if len(transcripts) == 1:
        transcript = next(iter(transcripts))
        plot_path = outdir / "plots" / f"{safe_name(transcript)}.degradome_cpm.png"
        plot_transcript(transcript, len(transcripts[transcript]), sample_data, cleavage, config, plot_path)
        plots.append(plot_path)
    else:
        for row in sorted(summary, key=lambda r: float(r["best_mean_degradome_cpm"]), reverse=True)[: config.max_transcript_plots]:
            if int(row["cleavage_supported_site_count"]) <= 0:
                continue
            transcript = str(row["transcript"])
            plot_path = outdir / "plots" / f"{safe_name(transcript)}.degradome_cpm.png"
            plot_transcript(transcript, len(transcripts[transcript]), sample_data, cleavage, config, plot_path)
            plots.append(plot_path)

    result = {
        "outdir": str(outdir),
        "srna_count": len(srnas),
        "transcript_count": len(transcripts),
        "sample_count": len(config.samples),
        "target_count": len(hits),
        "cleavage_count": len(cleavage),
        "significant_cleavage_count": sum(1 for row in cleavage if float(row.get("selected_pvalue") or 1.0) < config.pvalue_cutoff),
        "pvalue_method": pvalue_method,
        "pvalue_cutoff": config.pvalue_cutoff,
        "plots": [str(path) for path in plots],
        "outputs": {
            "cleavage_sites": str(cleavage_path),
            "target_predictions": str(target_path),
            "degradome_density": str(peak_path),
            "transcript_summary": str(summary_path),
            "srna_summary": str(srna_summary_path),
            "overall_summary": str(overall_summary_path),
            "summary_plot": str(summary_plot),
        },
    }
    report_path = outdir / "degradome_analysis_report.html"
    result["outputs"]["report"] = str(report_path)
    write_report(report_path, result)
    if log:
        log(f"Wrote degradome outputs to {outdir}.", "info")
    return result
