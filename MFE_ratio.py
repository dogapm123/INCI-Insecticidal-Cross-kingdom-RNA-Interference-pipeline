#!/usr/bin/env python3
"""GSTAr-style sRNA/transcript interaction analysis.

This module mirrors the core analysis performed by CleaveLand4's GSTAr:
RNAplex alignment, perfect-match MFE normalization, MFE ratio filtering,
Allen et al. scoring, slice-site estimation, and paired/unpaired region
annotation.  It also adds richer machine-readable output and lightweight
SVG/PNG visualizations of the sRNA-transcript match.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union


GSTAR_COLUMNS = [
    "Query",
    "Transcript",
    "TStart",
    "TStop",
    "TSlice",
    "MFEperfect",
    "MFEsite",
    "MFEratio",
    "AllenScore",
    "Paired",
    "Unpaired",
    "Structure",
    "Sequence",
]

EXTRA_COLUMNS = [
    "MatchPattern",
    "PairCount",
    "GUWobbleCount",
    "MismatchCount",
    "BulgeCount",
    "QueryLength",
    "TranscriptSiteLength",
    "VisualizationSVG",
]


@dataclass
class InteractionHit:
    query: str
    transcript: str
    t_start: int
    t_stop: int
    t_slice: int
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
    query_length: int
    transcript_site_length: int
    alignment_pretty: str
    visualization_svg: str

    def to_gstar_row(self) -> Dict[str, object]:
        return {
            "Query": self.query,
            "Transcript": self.transcript,
            "TStart": self.t_start,
            "TStop": self.t_stop,
            "TSlice": self.t_slice,
            "MFEperfect": self.mfe_perfect,
            "MFEsite": self.mfe_site,
            "MFEratio": self.mfe_ratio,
            "AllenScore": self.allen_score,
            "Paired": self.paired,
            "Unpaired": self.unpaired,
            "Structure": self.structure,
            "Sequence": self.sequence,
        }

    def to_output_row(self, visualization_path: str = "") -> Dict[str, object]:
        row = self.to_gstar_row()
        row.update(
            {
                "MatchPattern": self.match_pattern,
                "PairCount": self.pair_count,
                "GUWobbleCount": self.gu_wobble_count,
                "MismatchCount": self.mismatch_count,
                "BulgeCount": self.bulge_count,
                "QueryLength": self.query_length,
                "TranscriptSiteLength": self.transcript_site_length,
                "VisualizationSVG": visualization_path,
            }
        )
        return row


def compute_mfe_ratio(
    srna_sequence: str,
    transcript_sequence: str,
    query_name: str = "query",
    transcript_name: str = "transcript",
    mfe_ratio_cutoff: Optional[float] = None,
    rnaplex_path: str = "RNAplex",
) -> Dict[str, object]:
    """Analyze one sRNA/transcript pair and return the best GSTAr-style hit.

    The returned dictionary keeps the legacy-friendly ``alignment_pretty`` key
    used by ``main.py`` and includes GSTAr-like columns plus visualization data.
    """

    hits = analyze_interactions(
        {query_name: srna_sequence},
        {transcript_name: transcript_sequence},
        mfe_ratio_cutoff=mfe_ratio_cutoff,
        rnaplex_path=rnaplex_path,
    )
    if not hits:
        perfect = get_perfect_mfe(normalize_rna(srna_sequence), rnaplex_path)
        return {
            "query": query_name,
            "transcript": transcript_name,
            "mfe_perfect": perfect,
            "hit": None,
            "alignment_pretty": "No qualifying RNAplex interaction found.",
        }
    return hit_to_dict(hits[0])


def run_mfe_analysis(
    srna_sequence: Optional[str] = None,
    transcript_sequence: Optional[str] = None,
    srna_fasta: Optional[Union[str, Path]] = None,
    transcript_fasta: Optional[Union[str, Path]] = None,
    query_name: str = "query",
    transcript_name: str = "transcript",
    output_prefix: Optional[Union[str, Path]] = None,
    mfe_ratio_cutoff: Optional[float] = 0.70,
    sort_by: str = "mfe_ratio",
    top_n: int = 0,
    rnaplex_path: str = "RNAplex",
    write_png: bool = False,
    write_visualizations: bool = True,
    visualization_limit: int = 100,
) -> List[Dict[str, object]]:
    """Run GSTAr-style analysis directly from Python.

    Provide either both raw sequences::

        run_mfe_analysis(
            srna_sequence="UGAGGUAGUAGGUUGUAUAGUU",
            transcript_sequence="AAACUAUACAACCUACUACCUCAUUU",
        )

    or both FASTA files::

        run_mfe_analysis(
            srna_fasta="srna.fasta",
            transcript_fasta="transcripts.fasta",
            output_prefix="mfe_ratio_results",
        )

    Returns a list of dictionaries.  If ``output_prefix`` is supplied, the same
    TSV, JSON, SVG, and summary PNG outputs as the command-line interface are
    also written.
    """

    using_sequences = srna_sequence is not None or transcript_sequence is not None
    using_fastas = srna_fasta is not None or transcript_fasta is not None
    if using_sequences and using_fastas:
        raise ValueError("Use either raw sequences or FASTA files, not both.")
    if using_sequences:
        if srna_sequence is None or transcript_sequence is None:
            raise ValueError("Both srna_sequence and transcript_sequence are required.")
        queries = {query_name: srna_sequence}
        transcripts = {transcript_name: transcript_sequence}
    elif using_fastas:
        if srna_fasta is None or transcript_fasta is None:
            raise ValueError("Both srna_fasta and transcript_fasta are required.")
        queries = read_fasta(Path(srna_fasta))
        transcripts = read_fasta(Path(transcript_fasta))
    else:
        raise ValueError(
            "Provide either srna_sequence plus transcript_sequence, or srna_fasta plus transcript_fasta."
        )

    hits = analyze_interactions(
        queries,
        transcripts,
        mfe_ratio_cutoff=mfe_ratio_cutoff,
        sort_by=sort_by,
        rnaplex_path=rnaplex_path,
        write_visualizations=False,
    )
    if top_n:
        hits = hits[:top_n]
    if output_prefix is not None:
        write_outputs(
            hits,
            Path(output_prefix),
            write_png=write_png,
            write_visualizations=write_visualizations,
            visualization_limit=visualization_limit,
        )
    return [hit_to_dict(hit) for hit in hits]


def analyze_interactions(
    queries: Dict[str, str],
    transcripts: Dict[str, str],
    mfe_ratio_cutoff: Optional[float] = 0.70,
    sort_by: str = "mfe_ratio",
    rnaplex_path: str = "RNAplex",
    write_visualizations: bool = False,
) -> List[InteractionHit]:
    """Analyze all query/transcript combinations and return non-redundant hits."""

    if not shutil.which(rnaplex_path):
        raise RuntimeError(
            f"RNAplex was not found at '{rnaplex_path}'. Install ViennaRNA or pass --rnaplex."
        )

    all_hits: List[InteractionHit] = []
    normalized_transcripts = {
        transcript_name: normalize_rna(transcript_seq)
        for transcript_name, transcript_seq in transcripts.items()
    }
    for query_name, query_seq in queries.items():
        query_rna = normalize_rna(query_seq)
        if not 15 <= len(query_rna) <= 26:
            # GSTAr is designed for 15-26 nt queries, but short examples are
            # still useful in this project, so we do not reject them.
            pass
        perfect_mfe = get_perfect_mfe(query_rna, rnaplex_path)
        int_length = len(query_rna) + 10
        all_hits.extend(
            analyze_query_hits(
                query_name=query_name,
                query_seq=query_rna,
                transcripts=normalized_transcripts,
                perfect_mfe=perfect_mfe,
                int_length=int_length,
                mfe_ratio_cutoff=mfe_ratio_cutoff,
                rnaplex_path=rnaplex_path,
                write_visualizations=write_visualizations,
            )
        )

    return sort_and_dedupe_hits(all_hits, sort_by=sort_by)


def analyze_query_hits(
    query_name: str,
    query_seq: str,
    transcripts: Dict[str, str],
    perfect_mfe: float,
    int_length: int,
    mfe_ratio_cutoff: Optional[float],
    rnaplex_path: str,
    write_visualizations: bool = False,
) -> List[InteractionHit]:
    """Analyze one query against all transcripts in one RNAplex process.

    This mirrors GSTAr's fast path: one RNAplex input stream contains every
    transcript paired with the current query.
    """

    option_e = None
    if mfe_ratio_cutoff is not None and perfect_mfe:
        option_e = round(mfe_ratio_cutoff * perfect_mfe, 2)

    input_text = "".join(
        f">{transcript_name}\n{transcript_seq}\n>query\n{query_seq}\n"
        for transcript_name, transcript_seq in transcripts.items()
    )
    args = [rnaplex_path, "-f", "2", "-z", str(int_length)]
    if option_e is not None:
        args.extend(["-e", f"{option_e:.2f}"])
    output = run_rnaplex(args, input_text)

    hits: List[InteractionHit] = []
    current_transcript: Optional[str] = None
    for line in output.splitlines():
        if line.startswith(">"):
            header = line[1:].strip().split()[0]
            if header != "query":
                current_transcript = header
            continue
        if current_transcript is None or not re.match(r"^[.(]+&", line):
            continue
        transcript_seq = transcripts.get(current_transcript)
        if transcript_seq is None:
            continue
        hits.extend(
            interaction_hits_from_rnaplex_line(
                line=line,
                query_name=query_name,
                query_seq=query_seq,
                transcript_name=current_transcript,
                transcript_seq=transcript_seq,
                perfect_mfe=perfect_mfe,
                mfe_ratio_cutoff=mfe_ratio_cutoff,
                write_visualizations=write_visualizations,
            )
        )
    return sort_and_dedupe_hits(hits, sort_by="mfe_ratio")


def analyze_pair(
    query_name: str,
    query_seq: str,
    transcript_name: str,
    transcript_seq: str,
    perfect_mfe: float,
    int_length: int,
    mfe_ratio_cutoff: Optional[float],
    rnaplex_path: str,
) -> Optional[InteractionHit]:
    hits = analyze_pair_hits(
        query_name=query_name,
        query_seq=query_seq,
        transcript_name=transcript_name,
        transcript_seq=transcript_seq,
        perfect_mfe=perfect_mfe,
        int_length=int_length,
        mfe_ratio_cutoff=mfe_ratio_cutoff,
        rnaplex_path=rnaplex_path,
    )
    return hits[0] if hits else None


def analyze_pair_hits(
    query_name: str,
    query_seq: str,
    transcript_name: str,
    transcript_seq: str,
    perfect_mfe: float,
    int_length: int,
    mfe_ratio_cutoff: Optional[float],
    rnaplex_path: str,
    write_visualizations: bool = False,
) -> List[InteractionHit]:
    option_e = None
    if mfe_ratio_cutoff is not None and perfect_mfe:
        option_e = round(mfe_ratio_cutoff * perfect_mfe, 2)

    input_text = f">{transcript_name}\n{transcript_seq}\n>query\n{query_seq}\n"
    args = [rnaplex_path, "-f", "2", "-z", str(int_length)]
    if option_e is not None:
        args.extend(["-e", f"{option_e:.2f}"])
    output = run_rnaplex(args, input_text)

    hits: List[InteractionHit] = []
    for line in output.splitlines():
        if not re.match(r"^[.(]+&", line):
            continue
        hits.extend(
            interaction_hits_from_rnaplex_line(
                line=line,
                query_name=query_name,
                query_seq=query_seq,
                transcript_name=transcript_name,
                transcript_seq=transcript_seq,
                perfect_mfe=perfect_mfe,
                mfe_ratio_cutoff=mfe_ratio_cutoff,
                write_visualizations=write_visualizations,
            )
        )
    return sort_and_dedupe_hits(hits)


def interaction_hits_from_rnaplex_line(
    line: str,
    query_name: str,
    query_seq: str,
    transcript_name: str,
    transcript_seq: str,
    perfect_mfe: float,
    mfe_ratio_cutoff: Optional[float],
    write_visualizations: bool = False,
) -> List[InteractionHit]:
    parsed = parse_rnaplex_line(line)
    if parsed is None or not perfect_mfe:
        return []
    plex_brax, local_pos, site_mfe = parsed
    mfe_ratio = min(site_mfe / perfect_mfe, 1.0)
    if mfe_ratio_cutoff is not None and mfe_ratio < mfe_ratio_cutoff:
        return []

    padded_brax, adjusted_pos = pad_plex_brax(plex_brax, local_pos, query_seq)
    if adjusted_pos[0] < 1 or len(transcript_seq) < adjusted_pos[1]:
        return []

    tx_site_seq = transcript_seq[adjusted_pos[0] - 1 : adjusted_pos[1]]
    if not re.match(r"^[AUCG]+$", tx_site_seq):
        return []

    ungapped = f"{tx_site_seq}&{query_seq}"
    trimmed_brax, trimmed_seq, adjusted_pos = no_trailing(
        padded_brax, ungapped, adjusted_pos
    )
    structure, sequence = gapify(trimmed_brax, trimmed_seq)
    if not quality_control(structure, sequence):
        return []

    slice_site = compute_slice_site(structure, adjusted_pos)
    allen = allen_score(structure, sequence)
    paired, unpaired = assess_pairing(structure, adjusted_pos[0], query_seq)
    metrics = alignment_metrics(structure, sequence)
    pretty = pretty_alignment(
        query_name,
        transcript_name,
        adjusted_pos,
        slice_site,
        perfect_mfe,
        site_mfe,
        mfe_ratio,
        allen,
        paired,
        unpaired,
        structure,
        sequence,
    )
    svg = ""
    if write_visualizations:
        svg = alignment_svg(
            query_name,
            transcript_name,
            adjusted_pos,
            slice_site,
            perfect_mfe,
            site_mfe,
            mfe_ratio,
            allen,
            structure,
            sequence,
        )

    return [
        InteractionHit(
            query=query_name,
            transcript=transcript_name,
            t_start=adjusted_pos[0],
            t_stop=adjusted_pos[1],
            t_slice=slice_site,
            mfe_perfect=round(perfect_mfe, 4),
            mfe_site=round(site_mfe, 4),
            mfe_ratio=round(mfe_ratio, 6),
            allen_score=round(allen, 3),
            paired=paired or "NA",
            unpaired=unpaired or "NA",
            structure=structure,
            sequence=sequence,
            match_pattern=metrics["pattern"],
            pair_count=metrics["pairs"],
            gu_wobble_count=metrics["gu_wobbles"],
            mismatch_count=metrics["mismatches"],
            bulge_count=metrics["bulges"],
            query_length=len(query_seq),
            transcript_site_length=adjusted_pos[1] - adjusted_pos[0] + 1,
            alignment_pretty=pretty,
            visualization_svg=svg,
        )
    ]


def normalize_rna(sequence: str) -> str:
    seq = re.sub(r"\s+", "", sequence.upper()).replace("T", "U")
    if not re.match(r"^[AUCGN]+$", seq):
        raise ValueError(f"Invalid RNA/DNA sequence: {sequence!r}")
    return seq


def reverse_complement_rna(sequence: str) -> str:
    return sequence.translate(str.maketrans("AUCG", "UAGC"))[::-1]


def run_rnaplex(args: Sequence[str], input_text: str) -> str:
    proc = subprocess.run(
        list(args),
        input=input_text,
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"RNAplex failed: {proc.stderr.strip() or proc.stdout.strip()}")
    return proc.stdout


def get_perfect_mfe(query_seq: str, rnaplex_path: str = "RNAplex") -> float:
    perfect = reverse_complement_rna(query_seq)
    output = run_rnaplex(
        [rnaplex_path],
        f">perfect\n{perfect}\n>Query\n{query_seq}\n",
    )
    for line in output.splitlines():
        if re.match(r"^[.(]+&", line):
            parsed = parse_rnaplex_line(line)
            if parsed:
                return parsed[2]
    raise RuntimeError(f"RNAplex did not return a perfect-match MFE for {query_seq}")


def parse_rnaplex_line(line: str) -> Optional[Tuple[str, List[int], float]]:
    mfe_match = re.search(r"\(([\s\-\d.]+)\)", line)
    pos_match = re.search(r"(\d+),(\d+)\s*:\s*(\d+),(\d+)", line)
    brax_match = re.match(r"([().]+&[().]+)", line)
    if not (mfe_match and pos_match and brax_match):
        return None
    return (
        brax_match.group(1),
        [int(pos_match.group(i)) for i in range(1, 5)],
        float(mfe_match.group(1)),
    )


def pad_plex_brax(inbrax: str, pos: List[int], qseq: str) -> Tuple[str, List[int]]:
    pos = pos[:]
    match = re.match(r"^([(.]+)&([).]+)$", inbrax)
    if not match:
        raise ValueError(f"Failed to split RNAplex structure: {inbrax}")
    t_part, q_part = match.groups()
    q_len = len(qseq)

    left_added = 0
    for _ in range(pos[2], 1, -1):
        q_part = "." + q_part
        t_part = t_part + "."
        left_added += 1

    right_added = 0
    for _ in range(pos[3], q_len):
        q_part = q_part + "."
        t_part = "." + t_part
        right_added += 1

    pos[0] -= right_added
    pos[1] += left_added
    pos[2] -= left_added
    pos[3] += right_added

    t_5_dots = len(re.match(r"^\.+", t_part).group(0)) if re.match(r"^\.+", t_part) else 0
    t_3_dots = len(re.search(r"\.+$", t_part).group(0)) if re.search(r"\.+$", t_part) else 0
    q_5_dots = len(re.match(r"^\.+", q_part).group(0)) if re.match(r"^\.+", q_part) else 0
    q_3_dots = len(re.search(r"\.+$", q_part).group(0)) if re.search(r"\.+$", q_part) else 0

    for _ in range(q_3_dots - t_5_dots):
        t_part = "." + t_part
        pos[0] -= 1
    for _ in range(q_5_dots - t_3_dots):
        t_part += "."
        pos[1] += 1
    return f"{t_part}&{q_part}", pos


def no_trailing(brax: str, ungapped: str, pos: List[int]) -> Tuple[str, str, List[int]]:
    pos = pos[:]
    bs = brax.split("&")
    ss = ungapped.split("&")
    n_left_end_dots = len(re.match(r"^\.+", bs[0]).group(0)) if re.match(r"^\.+", bs[0]) else 0
    n_left_mid_dots = len(re.search(r"\.+$", bs[0]).group(0)) if re.search(r"\.+$", bs[0]) else 0
    n_right_end_dots = len(re.search(r"\.+$", bs[1]).group(0)) if re.search(r"\.+$", bs[1]) else 0
    n_right_mid_dots = len(re.match(r"^\.+", bs[1]).group(0)) if re.match(r"^\.+", bs[1]) else 0

    n_end_trim = max(n_left_end_dots - n_right_end_dots, 0)
    n_mid_trim = max(n_left_mid_dots - n_right_mid_dots, 0)

    bleft = list(bs[0])
    sleft = list(ss[0])
    for _ in range(n_end_trim):
        bleft.pop(0)
        sleft.pop(0)
        pos[0] += 1
    for _ in range(n_mid_trim):
        bleft.pop()
        sleft.pop()
        pos[1] -= 1

    return f"{''.join(bleft)}&{bs[1]}", f"{''.join(sleft)}&{ss[1]}", pos


def get_left_right(brax: str) -> Dict[int, int]:
    pairs: Dict[int, int] = {}
    stack: List[int] = []
    for index, char in enumerate(brax, start=1):
        if char == "(":
            stack.append(index)
        elif char == ")" and stack:
            pairs[stack.pop()] = index
    return pairs


def gapify(inbrax: str, inseq: str) -> Tuple[str, str]:
    left_right = get_left_right(inbrax)
    last_left = 0
    last_right = len(inbrax) + 1
    brax_chars = list(inbrax)
    seq_chars = list(inseq)
    out_brax_left = ""
    out_brax_right = ""
    out_seq_left = ""
    out_seq_right = ""

    for i in range(1, len(inbrax) + 1):
        if i in left_right:
            left_delta = i - last_left + 1
            right_delta = last_right - left_right[i] + 1
            delta = left_delta - right_delta

            for x in range(last_left, i - 1):
                out_brax_left += brax_chars[x]
                out_seq_left += seq_chars[x]
            for x in range(last_right - 2, left_right[i] - 1, -1):
                out_brax_right = brax_chars[x] + out_brax_right
                out_seq_right = seq_chars[x] + out_seq_right

            if delta > 0:
                out_brax_right = ("-" * delta) + out_brax_right
                out_seq_right = ("-" * delta) + out_seq_right
            elif delta < 0:
                out_brax_left += "-" * abs(delta)
                out_seq_left += "-" * abs(delta)

            out_brax_left += brax_chars[i - 1]
            out_seq_left += seq_chars[i - 1]
            out_brax_right = brax_chars[left_right[i] - 1] + out_brax_right
            out_seq_right = seq_chars[left_right[i] - 1] + out_seq_right
            last_left = i
            last_right = left_right[i]
        elif brax_chars[i - 1] == "&":
            for x in range(last_left, i - 1):
                out_brax_left += brax_chars[x]
                out_seq_left += seq_chars[x]
            for x in range(last_right - 2, i - 1, -1):
                out_brax_right = brax_chars[x] + out_brax_right
                out_seq_right = seq_chars[x] + out_seq_right
            break

    return f"{out_brax_left}&{out_brax_right}", f"{out_seq_left}&{out_seq_right}"


def get_al_array(string: str, index: int) -> List[str]:
    return list(string.split("&")[index])


def quality_control(structure: str, sequence: str) -> bool:
    t_brax, q_brax = get_al_array(structure, 0), get_al_array(structure, 1)
    t_seq, q_seq = get_al_array(sequence, 0), get_al_array(sequence, 1)
    if not (len(t_brax) == len(q_brax) == len(t_seq) == len(q_seq)):
        return False
    if any(ch not in "AUCG-" for ch in t_seq + q_seq):
        return False
    if any(ch not in "(.-" for ch in t_brax):
        return False
    if any(ch not in ").-" for ch in q_brax):
        return False
    return t_brax.count("(") == q_brax.count(")")


def compute_slice_site(structure: str, local_pos: Sequence[int]) -> int:
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
        if real_q_pos == 10:
            return real_t_pos
    return real_t_pos


def allen_score(structure: str, sequence: str) -> float:
    score = 0.0
    t_brax = get_al_array(structure, 0)
    q_brax = get_al_array(structure, 1)
    t_seq = get_al_array(sequence, 0)
    q_seq = get_al_array(sequence, 1)
    q_pos = 0
    for i, qs in enumerate(q_seq):
        if qs != "-":
            q_pos += 1
        qb = q_brax[i]
        ts = t_seq.pop()
        tb = t_brax.pop()
        doubled = 2 <= q_pos <= 12
        if tb == "." or qb == ".":
            score += 2 if doubled else 1
        elif tb == "(" and qb == ")" and {ts, qs} == {"G", "U"}:
            score += 1.0 if doubled else 0.5
    return score


def assess_pairing(structure: str, tx_start: int, qseq: str) -> Tuple[str, str]:
    paired = check_pairs(structure, tx_start, qseq)
    unpaired = check_unpaired(structure, tx_start, qseq)
    return paired or "NA", annotate_unpaired(unpaired, qseq) if unpaired else "NA"


def check_pairs(structure: str, tx_start: int, qseq: str) -> str:
    t_brax = get_al_array(structure, 0)
    q_brax = list(reversed(get_al_array(structure, 1)))
    t_pos = tx_start - 1
    q_pos = len(qseq) + 1
    blocks: List[Tuple[int, int, int, int]] = []
    current: Optional[List[int]] = None
    for t, q in zip(t_brax, q_brax):
        if t != "-":
            t_pos += 1
        if q != "-":
            q_pos -= 1
        if t == "(" and q == ")":
            if current is None:
                current = [q_pos, q_pos, t_pos, t_pos]
            else:
                current[1] = q_pos
                current[3] = t_pos
        elif current is not None:
            blocks.append(tuple(current))
            current = None
    if current is not None:
        blocks.append(tuple(current))
    return ";".join(f"{q2}-{q1},{t2}-{t1}" for q1, q2, t1, t2 in reversed(blocks))


def check_unpaired(structure: str, tx_start: int, qseq: str) -> str:
    t_brax = get_al_array(structure, 0)
    q_brax = list(reversed(get_al_array(structure, 1)))
    t_pos = tx_start - 1
    q_pos = len(qseq) + 1
    blocks: List[Tuple[object, object, object, object]] = []
    current: Optional[List[object]] = None
    for t, q in zip(t_brax, q_brax):
        if t != "-":
            t_pos += 1
        if q != "-":
            q_pos -= 1
        if t == "(" and q == ")":
            if current is not None:
                blocks.append(tuple(current))
                current = None
            continue
        q_value: object = "x" if q == "-" else q_pos
        t_value: object = "x" if t == "-" else t_pos
        if current is None:
            current = [q_value, q_value, t_value, t_value]
        else:
            current[1] = q_value if current[1] == "x" or q_value == "x" else q_value
            current[3] = t_value if current[3] == "x" or t_value == "x" else t_value
            if current[0] == "x" and q_value != "x":
                current[0] = q_value
            if current[2] == "x" and t_value != "x":
                current[2] = t_value
    if current is not None:
        blocks.append(tuple(current))
    return ";".join(f"{q2}-{q1},{t2}-{t1}" for q1, q2, t1, t2 in reversed(blocks))


def annotate_unpaired(unpaired: str, qseq: str) -> str:
    out: List[str] = []
    max_q = len(qseq)
    for entry in unpaired.split(";"):
        match = re.match(r"^(\S+)-(\S+),(\S+)-(\S+)$", entry)
        if not match:
            out.append(f"{entry}[?]")
            continue
        q_start, q_stop, t_stop, t_start = match.groups()
        code = "?"
        if all(x.isdigit() for x in [q_start, q_stop, t_start, t_stop]):
            q_delta = int(q_stop) - int(q_start) + 1
            t_delta = int(t_stop) - int(t_start) + 1
            if int(q_start) == 1:
                code = "UP5"
            elif int(q_stop) == max_q:
                code = "UP3"
            elif q_delta == t_delta:
                code = "SIL"
            elif q_delta > t_delta:
                code = "AILq"
            else:
                code = "AILt"
        elif q_start == q_stop == "x":
            code = "BULt"
        elif t_start == t_stop == "x":
            code = "BULq"
        out.append(f"{entry}[{code}]")
    return ";".join(out)


def alignment_metrics(structure: str, sequence: str) -> Dict[str, object]:
    t_brax = get_al_array(structure, 0)
    q_brax = list(reversed(get_al_array(structure, 1)))
    t_seq = get_al_array(sequence, 0)
    q_seq = list(reversed(get_al_array(sequence, 1)))
    pattern = []
    pairs = gu = mismatches = bulges = 0
    for tb, qb, ts, qs in zip(t_brax, q_brax, t_seq, q_seq):
        if ts == "-" or qs == "-":
            pattern.append(" ")
            bulges += 1
        elif tb == "(" and qb == ")":
            pairs += 1
            if {ts, qs} == {"G", "U"}:
                pattern.append("o")
                gu += 1
            else:
                pattern.append("|")
        else:
            pattern.append("x")
            mismatches += 1
    return {
        "pattern": "".join(pattern),
        "pairs": pairs,
        "gu_wobbles": gu,
        "mismatches": mismatches,
        "bulges": bulges,
    }


def pretty_alignment(
    query_name: str,
    transcript_name: str,
    local_pos: Sequence[int],
    slice_site: int,
    perfect_mfe: float,
    site_mfe: float,
    mfe_ratio: float,
    allen: float,
    paired: str,
    unpaired: str,
    structure: str,
    sequence: str,
) -> str:
    t_seq = get_al_array(sequence, 0)
    q_seq = list(reversed(get_al_array(sequence, 1)))
    pattern = alignment_metrics(structure, sequence)["pattern"]
    lines = [
        "----------------------------------------------------------------",
        f"5' {''.join(t_seq)} 3'  Transcript: {transcript_name}:{local_pos[0]}-{local_pos[1]}  Slice Site:{slice_site}",
        f"   {pattern}",
        f"3' {''.join(q_seq)} 5'  Query: {query_name}",
        f"MFE of perfect match: {perfect_mfe:.2f}",
        f"MFE of this site: {site_mfe:.2f}",
        f"MFEratio: {mfe_ratio:.4f}",
        f"Allen et al. score: {allen:.3f}",
        f"Paired Regions (query5'-query3',transcript3'-transcript5'): {paired}",
        f"Unpaired Regions (query5'-query3',transcript3'-transcript5'): {unpaired}",
        "----------------------------------------------------------------",
    ]
    return "\n".join(lines)


def alignment_svg(
    query_name: str,
    transcript_name: str,
    local_pos: Sequence[int],
    slice_site: int,
    perfect_mfe: float,
    site_mfe: float,
    mfe_ratio: float,
    allen: float,
    structure: str,
    sequence: str,
) -> str:
    esc = html.escape
    t_brax = get_al_array(structure, 0)
    q_brax = list(reversed(get_al_array(structure, 1)))
    t_seq = get_al_array(sequence, 0)
    q_seq = list(reversed(get_al_array(sequence, 1)))

    annotations: List[Tuple[str, str]] = []
    max_text_len = 1
    for tb, qb, target_base, guide_base in zip(t_brax, q_brax, t_seq, q_seq):
        if target_base == "-" or guide_base == "-":
            label, klass = "gap", "gap"
        elif tb == "(" and qb == ")" and {target_base, guide_base} == {"G", "U"}:
            label, klass = "G-U", "wobble"
        elif tb == "(" and qb == ")":
            label, klass = "|", "match"
        else:
            label, klass = "MM", "mismatch"
        annotations.append((label, klass))
        max_text_len = max(max_text_len, len(label), len(target_base), len(guide_base))

    cell_w = max(34, max_text_len * 11 + 18)
    cell_h = 34
    left = 112
    top = 92
    table_w = cell_w * len(annotations)
    width = max(860, left + table_w + 36)
    height = 270

    target_pos = local_pos[0] - 1
    guide_pos = sum(1 for base in q_seq if base != "-") + 1
    target_pos_labels: List[str] = []
    guide_pos_labels: List[str] = []
    for target_base, guide_base in zip(t_seq, q_seq):
        if target_base != "-":
            target_pos += 1
            target_pos_labels.append(str(target_pos))
        else:
            target_pos_labels.append("")
        if guide_base != "-":
            guide_pos -= 1
            guide_pos_labels.append(str(guide_pos))
        else:
            guide_pos_labels.append("")

    def text_cell(x: int, y: int, value: str, fill: str, color: str = "#20323a", size: int = 16) -> str:
        return (
            f'<rect x="{x}" y="{y}" width="{cell_w}" height="{cell_h}" fill="{fill}" '
            f'stroke="#d8e0de" stroke-width="1"/>'
            f'<text x="{x + cell_w / 2:.1f}" y="{y + 22}" text-anchor="middle" '
            f'font-family="Menlo, Consolas, monospace" font-size="{size}" '
            f'font-weight="700" fill="{color}">{esc(value)}</text>'
        )

    cells: List[str] = []
    for index, (target_base, guide_base) in enumerate(zip(t_seq, q_seq)):
        x = left + index * cell_w
        label, klass = annotations[index]
        ann_fill = {
            "match": "#e8f3ee",
            "wobble": "#fff3d6",
            "mismatch": "#fde2dc",
            "gap": "#eceff3",
        }[klass]
        ann_color = {
            "match": "#1f6f4a",
            "wobble": "#8a5a00",
            "mismatch": "#9f2f22",
            "gap": "#5c6470",
        }[klass]
        cells.append(text_cell(x, top, target_base, "#eef7f8", "#184c5a"))
        cells.append(text_cell(x, top + cell_h, label, ann_fill, ann_color, size=14))
        cells.append(text_cell(x, top + (cell_h * 2), guide_base, "#fff1e7", "#7a3f18"))
        cells.append(
            f'<text x="{x + cell_w / 2:.1f}" y="{top - 8}" text-anchor="middle" '
            f'font-family="Arial, sans-serif" font-size="9" fill="#65747a">{esc(target_pos_labels[index])}</text>'
        )
        cells.append(
            f'<text x="{x + cell_w / 2:.1f}" y="{top + (cell_h * 3) + 15}" text-anchor="middle" '
            f'font-family="Arial, sans-serif" font-size="9" fill="#7a6a5d">{esc(guide_pos_labels[index])}</text>'
        )

    labels = f"""
  <text x="28" y="{top + 22}" font-family="Arial, sans-serif" font-size="13" font-weight="700" fill="#184c5a">Transcript 5' to 3'</text>
  <text x="28" y="{top + cell_h + 22}" font-family="Arial, sans-serif" font-size="13" font-weight="700" fill="#52636a">Annotation</text>
  <text x="28" y="{top + (cell_h * 2) + 22}" font-family="Arial, sans-serif" font-size="13" font-weight="700" fill="#7a3f18">Guide 3' to 5'</text>
  <text x="{left}" y="{top - 24}" font-family="Arial, sans-serif" font-size="10" fill="#65747a">transcript position</text>
  <text x="{left}" y="{top + (cell_h * 3) + 31}" font-family="Arial, sans-serif" font-size="10" fill="#7a6a5d">guide position</text>"""

    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect width="100%" height="100%" fill="#fbfbf8"/>
  <text x="28" y="30" font-family="Arial, sans-serif" font-size="16" font-weight="700" fill="#233">{esc(query_name)} vs {esc(transcript_name)}</text>
  <text x="28" y="54" font-family="Arial, sans-serif" font-size="12" fill="#566">Transcript {local_pos[0]}-{local_pos[1]} | Slice {slice_site} | MFEratio {mfe_ratio:.4f} | Allen {allen:.3f} | MFE {site_mfe:.2f}/{perfect_mfe:.2f}</text>
  {labels}
  {''.join(cells)}
  <text x="28" y="{height - 24}" font-family="Arial, sans-serif" font-size="12" fill="#566">Legend: | canonical pair, G-U wobble, MM mismatch, gap bulge/unpaired column</text>
</svg>
"""


def sort_and_dedupe_hits(hits: List[InteractionHit], sort_by: str = "mfe_ratio") -> List[InteractionHit]:
    best_by_slice: Dict[Tuple[str, str, int], InteractionHit] = {}
    for hit in hits:
        key = (hit.query, hit.transcript, hit.t_slice)
        existing = best_by_slice.get(key)
        if existing is None:
            best_by_slice[key] = hit
        elif sort_by == "allen" and hit.allen_score < existing.allen_score:
            best_by_slice[key] = hit
        elif sort_by != "allen" and hit.mfe_ratio > existing.mfe_ratio:
            best_by_slice[key] = hit

    if sort_by == "allen":
        return sorted(best_by_slice.values(), key=lambda h: (h.allen_score, -h.mfe_ratio))
    return sorted(best_by_slice.values(), key=lambda h: (-h.mfe_ratio, h.allen_score))


def hit_to_dict(hit: InteractionHit) -> Dict[str, object]:
    data = asdict(hit)
    data.update(hit.to_gstar_row())
    data.update(hit.to_output_row(hit.visualization_svg))
    return data


def read_fasta(path: Path) -> Dict[str, str]:
    records: Dict[str, List[str]] = {}
    current: Optional[str] = None
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if not current:
                    raise ValueError(f"Empty FASTA header in {path}")
                records[current] = []
            elif current is None:
                raise ValueError(f"Sequence found before FASTA header in {path}")
            else:
                records[current].append(line)
    return {name: "".join(parts) for name, parts in records.items()}


def write_outputs(
    hits: List[InteractionHit],
    output_prefix: Path,
    write_png: bool = False,
    write_visualizations: bool = False,
    visualization_limit: int = 100,
) -> None:
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    svg_dir = output_prefix.parent / f"{output_prefix.name}_visualizations"
    if write_visualizations:
        svg_dir.mkdir(exist_ok=True)
    rows = []
    json_rows = []
    max_visualizations = max(0, int(visualization_limit))
    for index, hit in enumerate(hits, start=1):
        svg_path = ""
        if write_visualizations and index <= max_visualizations:
            svg = hit.visualization_svg or alignment_svg(
                hit.query,
                hit.transcript,
                [hit.t_start, hit.t_stop],
                hit.t_slice,
                hit.mfe_perfect,
                hit.mfe_site,
                hit.mfe_ratio,
                hit.allen_score,
                hit.structure,
                hit.sequence,
            )
            stem = safe_filename(f"{index:04d}_{hit.query}_{hit.transcript}_{hit.t_slice}")
            svg_file = svg_dir / f"{stem}.svg"
            svg_file.write_text(svg)
            svg_path = str(svg_file)
        rows.append(hit.to_output_row(svg_path))
        item = hit_to_dict(hit)
        item["visualization_svg_path"] = svg_path
        json_rows.append(item)

    tsv_path = output_prefix.with_suffix(".tsv")
    with tsv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=GSTAR_COLUMNS + EXTRA_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    json_path = output_prefix.with_suffix(".json")
    json_path.write_text(json.dumps(json_rows, indent=2))

    if write_png and hits:
        write_summary_plot(hits, output_prefix.with_suffix(".summary.png"))


def write_summary_plot(hits: List[InteractionHit], path: Path) -> None:
    mpl_dir = Path(tempfile.gettempdir()) / "inci_pipeline_matplotlib"
    mpl_dir.mkdir(exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_dir))
    from plot_style import apply_matplotlib_style

    style = apply_matplotlib_style()
    import matplotlib.pyplot as plt

    top_hits = hits[: min(25, len(hits))]
    labels = [f"{h.query}\n{h.transcript}:{h.t_slice}" for h in top_hits]
    ratios = [h.mfe_ratio for h in top_hits]
    allen = [h.allen_score for h in top_hits]

    fig_width = max(float(style["figure_width"]), len(top_hits) * 0.55)
    fig, ax1 = plt.subplots(figsize=(fig_width, style["figure_height"]))
    x = range(len(top_hits))
    ax1.bar(x, ratios, color="#2f7f73", label="MFE ratio")
    ax1.set_ylim(0, 1.05)
    ax1.set_ylabel("MFE ratio")
    ax1.set_xticks(list(x))
    ax1.set_xticklabels(labels, rotation=55, ha="right", fontsize=max(6.0, float(style["font_size"]) - 2.0))
    ax2 = ax1.twinx()
    ax2.plot(list(x), allen, color="#b45f2a", marker="o", label="Allen score")
    ax2.set_ylabel("Allen score")
    ax1.set_title("Top sRNA-transcript matches")
    fig.tight_layout()
    fig.savefig(path, dpi=style["dpi"])
    plt.close(fig)


def safe_filename(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") or "hit"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="GSTAr-style RNAplex analysis for sRNA-transcript matches."
    )
    parser.add_argument("--srna", required=True, type=Path, help="FASTA file of sRNA queries.")
    parser.add_argument(
        "--transcripts", required=True, type=Path, help="FASTA file of transcript sequences."
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("mfe_ratio_results"),
        help="Prefix for TSV, JSON, and visualization outputs.",
    )
    parser.add_argument(
        "--ratio-cutoff",
        type=float,
        default=0.70,
        help="Minimum MFEsite/MFEperfect ratio to keep, matching GSTAr -r behavior.",
    )
    parser.add_argument(
        "--sort-by",
        choices=["mfe_ratio", "allen"],
        default="mfe_ratio",
        help="Sort by descending MFE ratio or ascending Allen score.",
    )
    parser.add_argument("--top-n", type=int, default=0, help="Keep only the top N hits.")
    parser.add_argument("--rnaplex", default="RNAplex", help="Path to RNAplex executable.")
    parser.add_argument("--no-visualizations", action="store_true", help="Skip per-hit SVG visualizations.")
    parser.add_argument(
        "--visualization-limit",
        type=int,
        default=100,
        help="Maximum number of top sorted hits to visualize. Default: 100.",
    )
    parser.add_argument("--write-visualizations", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--write-png", action="store_true", help="Write the summary PNG plot.")
    parser.add_argument("--no-png", action="store_true", help=argparse.SUPPRESS)
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()
    queries = read_fasta(args.srna)
    transcripts = read_fasta(args.transcripts)
    hits = analyze_interactions(
        queries,
        transcripts,
        mfe_ratio_cutoff=args.ratio_cutoff,
        sort_by=args.sort_by,
        rnaplex_path=args.rnaplex,
        write_visualizations=False,
    )
    if args.top_n:
        hits = hits[: args.top_n]
    write_outputs(
        hits,
        args.output_prefix,
        write_png=args.write_png and not args.no_png,
        write_visualizations=not args.no_visualizations,
        visualization_limit=args.visualization_limit,
    )
    print(f"Wrote {len(hits)} hits to {args.output_prefix.with_suffix('.tsv')}")


if __name__ == "__main__":
    main()
