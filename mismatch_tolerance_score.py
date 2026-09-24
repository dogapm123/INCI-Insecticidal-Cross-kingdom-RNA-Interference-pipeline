#!/usr/bin/env python3
"""Score siRNA-guide/transcript mismatches using literature tolerance rules.

This script is intended as a post-processor for ``MFE_ratio.py`` outputs.  It
reads the reported interaction ``Structure`` and ``Sequence`` fields, maps each
aligned column back to a 5' guide-strand coordinate, and marks mismatches with
star scores:

* tolerated nucleotide pair
** poorly tolerated nucleotide pair
additional * when the guide position is a mismatch-sensitive position

The defaults are based mainly on Huang et al. 2009 (NAR 37:7560-7569), with the
broader terminal/central-region pattern from Du et al. 2005, Schwarz et al.
2006, Dahlgren et al. 2008, and Wei et al. 2012.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


CANONICAL_PAIRS = {
    ("A", "U"),
    ("U", "A"),
    ("G", "C"),
    ("C", "G"),
}

# Guide base -> target base -> base-pair star score.
# 1 = comparatively tolerated, 2 = comparatively poorly tolerated.
PAIR_STAR_SCORE: Dict[Tuple[str, str], int] = {
    ("A", "A"): 2,
    ("A", "C"): 1,
    ("A", "G"): 2,
    ("C", "A"): 1,
    ("C", "C"): 2,
    ("C", "U"): 2,
    ("G", "A"): 2,
    ("G", "G"): 2,
    ("G", "U"): 1,
    ("U", "C"): 2,
    ("U", "G"): 1,
    ("U", "U"): 1,
}

PAIR_CLASS: Dict[Tuple[str, str], str] = {
    ("A", "C"): "tolerated_pair:Huang2009_A:C",
    ("C", "A"): "tolerated_pair:Huang2009_C:A",
    ("G", "U"): "tolerated_pair:G:U_wobble",
    ("U", "G"): "tolerated_pair:U:G_wobble",
    ("U", "U"): "tolerated_pair:Huang2009_U:U",
    ("A", "G"): "poorly_tolerated_pair:Huang2009_A:G",
    ("G", "G"): "poorly_tolerated_pair:Huang2009_G:G",
    ("U", "C"): "poorly_tolerated_pair:Huang2009_U:C",
    ("C", "C"): "poorly_tolerated_pair:Huang2009_C:C",
}

# Huang et al. 2009 guide-base-specific sensitive positions.  These are the
# positions where a mismatch is more likely to disrupt silencing and is useful
# for allele discrimination.
GUIDE_BASE_SENSITIVE_POSITIONS: Dict[str, set[int]] = {
    "A": {9, 10, 12, 14, 15},
    "G": {8, 9, 10, 11, 12, 13},
    "U": {9, 11, 12, 13, 15, 16},
    "C": {10, 11, 12, 13, 16},
}

# Literature consensus fallback used when a guide base is ambiguous or when the
# user wants a broad siRNA cleavage-site rule rather than guide-base-specific
# Huang rules.
GENERAL_SENSITIVE_POSITIONS = {9, 10, 11, 12, 13, 16}
TERMINAL_TOLERANT_POSITIONS = {1, 2, 18, 19, 20, 21}


@dataclass
class MismatchCall:
    guide_pos: int
    guide_base: str
    target_base: str
    pair: str
    pair_stars: int
    position_stars: int
    total_stars: int
    reason: str
    target_pos: Optional[int] = None

    @property
    def stars(self) -> str:
        return "*" * self.total_stars

    def compact(self) -> str:
        target = f",t{self.target_pos}" if self.target_pos is not None else ""
        return (
            f"g{self.guide_pos}{target}:{self.pair}:"
            f"{self.stars}:{self.reason}"
        )


@dataclass
class BulgeCall:
    guide_pos: Optional[int]
    target_pos: Optional[int]
    guide_base: str
    target_base: str
    stars: int
    reason: str

    def compact(self) -> str:
        guide = f"g{self.guide_pos}" if self.guide_pos is not None else "g-"
        target = f"t{self.target_pos}" if self.target_pos is not None else "t-"
        return f"{guide},{target}:{self.guide_base}:{self.target_base}:{'*' * self.stars}:{self.reason}"


@dataclass
class PairAnnotation:
    source_row: int
    query: str
    transcript: str
    guide_pos: Optional[int]
    target_pos: Optional[int]
    guide_base: str
    target_base: str
    pair: str
    pair_type: str
    pair_stars: int
    position_stars: int
    total_stars: int
    reason: str

    @property
    def stars(self) -> str:
        return "*" * self.total_stars

    def as_row(self) -> Dict[str, object]:
        return {
            "SourceRow": self.source_row,
            "Query": self.query,
            "Transcript": self.transcript,
            "GuidePosition": self.guide_pos if self.guide_pos is not None else "",
            "TargetPosition": self.target_pos if self.target_pos is not None else "",
            "GuideBase": self.guide_base,
            "TargetBase": self.target_base,
            "Pair": self.pair,
            "PairType": self.pair_type,
            "PairStars": self.pair_stars,
            "PositionStars": self.position_stars,
            "TotalStars": self.total_stars,
            "StarAnnotation": self.stars,
            "Reason": self.reason,
        }


def normalize_rna(sequence: str) -> str:
    return re.sub(r"\s+", "", str(sequence).upper()).replace("T", "U")


def is_sensitive_position(guide_base: str, guide_pos: int, mode: str) -> bool:
    if mode == "huang":
        return guide_pos in GUIDE_BASE_SENSITIVE_POSITIONS.get(
            guide_base, GENERAL_SENSITIVE_POSITIONS
        )
    if mode == "general":
        return guide_pos in GENERAL_SENSITIVE_POSITIONS
    return False


def position_reason(guide_pos: int, guide_base: str, mode: str) -> Tuple[int, str]:
    if is_sensitive_position(guide_base, guide_pos, mode):
        return 1, "sensitive_position"
    if guide_pos in TERMINAL_TOLERANT_POSITIONS:
        return 0, "terminal_or_overhang_tolerant_position"
    return 0, "non_sensitive_position"


def pair_reason(guide_base: str, target_base: str) -> Tuple[int, str]:
    pair = (guide_base, target_base)
    if pair in CANONICAL_PAIRS:
        return 0, "canonical_pair"
    if pair in PAIR_STAR_SCORE:
        return PAIR_STAR_SCORE[pair], PAIR_CLASS.get(pair, "mismatch_pair")
    purines = {"A", "G"}
    pyrimidines = {"C", "U"}
    if guide_base in purines and target_base in purines:
        return 2, "poorly_tolerated_pair:purine:purine"
    if guide_base in pyrimidines and target_base in pyrimidines:
        return 2, "poorly_tolerated_pair:pyrimidine:pyrimidine"
    return 1, "moderately_tolerated_pair:purine:pyrimidine"


def annotate_pair(
    guide_base: str,
    target_base: str,
    guide_pos: Optional[int],
    target_pos: Optional[int],
    position_mode: str,
    source_row: int = 1,
    query: str = "",
    transcript: str = "",
) -> PairAnnotation:
    guide_base = normalize_rna(guide_base)
    target_base = normalize_rna(target_base)
    if guide_base == "-" or target_base == "-":
        return PairAnnotation(
            source_row=source_row,
            query=query,
            transcript=transcript,
            guide_pos=guide_pos,
            target_pos=target_pos,
            guide_base=guide_base,
            target_base=target_base,
            pair=f"{guide_base}:{target_base}",
            pair_type="bulge_or_gap",
            pair_stars=3,
            position_stars=0,
            total_stars=3,
            reason="bulge_or_gap_in_alignment",
        )

    pair_stars, pair_note = pair_reason(guide_base, target_base)
    if pair_stars == 0:
        pair_type = "canonical_match"
    elif pair_stars == 1:
        pair_type = "tolerated_mismatch"
    else:
        pair_type = "less_tolerated_mismatch"

    position_stars = 0
    position_note = "position_not_evaluated"
    if guide_pos is not None:
        position_stars, position_note = position_reason(guide_pos, guide_base, position_mode)

    return PairAnnotation(
        source_row=source_row,
        query=query,
        transcript=transcript,
        guide_pos=guide_pos,
        target_pos=target_pos,
        guide_base=guide_base,
        target_base=target_base,
        pair=f"{guide_base}:{target_base}",
        pair_type=pair_type,
        pair_stars=pair_stars,
        position_stars=position_stars if pair_stars else 0,
        total_stars=pair_stars + (position_stars if pair_stars else 0),
        reason=f"{pair_note}+{position_note}" if pair_stars else pair_note,
    )


def score_pair(
    guide_base: str,
    target_base: str,
    guide_pos: int,
    target_pos: Optional[int],
    position_mode: str,
) -> Optional[MismatchCall]:
    guide_base = normalize_rna(guide_base)
    target_base = normalize_rna(target_base)
    if not guide_base or not target_base or guide_base == "-" or target_base == "-":
        return None
    pair_stars, pair_note = pair_reason(guide_base, target_base)
    if pair_stars == 0:
        return None
    pos_stars, pos_note = position_reason(guide_pos, guide_base, position_mode)
    return MismatchCall(
        guide_pos=guide_pos,
        guide_base=guide_base,
        target_base=target_base,
        pair=f"{guide_base}:{target_base}",
        pair_stars=pair_stars,
        position_stars=pos_stars,
        total_stars=pair_stars + pos_stars,
        reason=f"{pair_note}+{pos_note}",
        target_pos=target_pos,
    )


def split_alignment_field(value: str, field_name: str) -> Tuple[List[str], List[str]]:
    parts = str(value).strip().split("&")
    if len(parts) != 2:
        raise ValueError(f"{field_name} must contain exactly one '&': {value!r}")
    return list(parts[0]), list(parts[1])


def calls_from_mfe_row(
    row: Dict[str, object],
    position_mode: str = "huang",
    include_bulges: bool = True,
    bulge_stars: int = 3,
) -> Tuple[List[MismatchCall], List[BulgeCall]]:
    """Return mismatch and bulge calls from one MFE_ratio.py output row."""
    annotations = pair_annotations_from_mfe_row(
        row,
        source_row=1,
        position_mode=position_mode,
        include_bulges=include_bulges,
        bulge_stars=bulge_stars,
    )
    return calls_from_annotations(annotations, bulge_stars=bulge_stars)


def calls_from_annotations(
    annotations: Sequence[PairAnnotation],
    bulge_stars: int = 3,
) -> Tuple[List[MismatchCall], List[BulgeCall]]:
    mismatches = [
        MismatchCall(
            guide_pos=ann.guide_pos or 0,
            guide_base=ann.guide_base,
            target_base=ann.target_base,
            pair=ann.pair,
            pair_stars=ann.pair_stars,
            position_stars=ann.position_stars,
            total_stars=ann.total_stars,
            reason=ann.reason,
            target_pos=ann.target_pos,
        )
        for ann in annotations
        if ann.pair_type in {"tolerated_mismatch", "less_tolerated_mismatch"}
        and ann.guide_pos is not None
    ]
    bulges = [
        BulgeCall(
            guide_pos=ann.guide_pos,
            target_pos=ann.target_pos,
            guide_base=ann.guide_base,
            target_base=ann.target_base,
            stars=bulge_stars,
            reason=ann.reason,
        )
        for ann in annotations
        if ann.pair_type == "bulge_or_gap"
    ]
    return mismatches, bulges


def pair_annotations_from_mfe_row(
    row: Dict[str, object],
    source_row: int,
    position_mode: str = "huang",
    include_bulges: bool = True,
    bulge_stars: int = 3,
) -> List[PairAnnotation]:
    """Return one annotation row for every aligned guide/target pair."""

    structure = row.get("Structure", "")
    sequence = row.get("Sequence", "")
    if not structure or not sequence:
        raise ValueError("MFE row is missing Structure or Sequence")

    t_struct, q_struct_right = split_alignment_field(str(structure), "Structure")
    t_seq, q_seq_right = split_alignment_field(str(sequence), "Sequence")
    q_struct = list(reversed(q_struct_right))
    q_seq = list(reversed(q_seq_right))

    if not (len(t_struct) == len(q_struct) == len(t_seq) == len(q_seq)):
        raise ValueError(
            "Structure and Sequence alignment lengths disagree after guide reversal"
        )

    query_len = int(row.get("QueryLength") or sum(base != "-" for base in q_seq))
    query = str(row.get("Query", ""))
    transcript = str(row.get("Transcript", ""))
    target_start = safe_int(row.get("TStart"))
    target_pos = target_start - 1 if target_start is not None else None
    guide_pos = query_len + 1
    annotations: List[PairAnnotation] = []

    for tb, qb, target_base, guide_base in zip(t_struct, q_struct, t_seq, q_seq):
        current_target_pos: Optional[int] = None
        current_guide_pos: Optional[int] = None

        if target_base != "-":
            target_pos = target_pos + 1 if target_pos is not None else None
            current_target_pos = target_pos
        if guide_base != "-":
            guide_pos -= 1
            current_guide_pos = guide_pos

        if target_base == "-" or guide_base == "-":
            if include_bulges:
                ann = annotate_pair(
                    guide_base=guide_base,
                    target_base=target_base,
                    guide_pos=current_guide_pos,
                    target_pos=current_target_pos,
                    position_mode=position_mode,
                    source_row=source_row,
                    query=query,
                    transcript=transcript,
                )
                ann.pair_stars = bulge_stars
                ann.total_stars = bulge_stars
                annotations.append(ann)
            continue

        annotations.append(
            annotate_pair(
                guide_base=guide_base,
                target_base=target_base,
                guide_pos=current_guide_pos,
                target_pos=current_target_pos,
                position_mode=position_mode,
                source_row=source_row,
                query=query,
                transcript=transcript,
            )
        )

    return annotations


def pair_annotations_from_direct_sequences(
    guide: str,
    target_site: str,
    position_mode: str = "huang",
    source_row: int = 1,
    query: str = "direct_guide",
    transcript: str = "direct_target",
) -> List[PairAnnotation]:
    """Annotate every guide/target nucleotide pair in direct sequence mode."""

    guide = normalize_rna(guide)
    target_site = normalize_rna(target_site)
    if len(guide) != len(target_site):
        raise ValueError(
            f"Direct mode requires equal lengths, got guide={len(guide)} and target={len(target_site)}"
        )
    annotations: List[PairAnnotation] = []
    for guide_pos, guide_base in enumerate(guide, start=1):
        target_index = len(target_site) - guide_pos
        target_base = target_site[target_index]
        annotations.append(
            annotate_pair(
                guide_base=guide_base,
                target_base=target_base,
                guide_pos=guide_pos,
                target_pos=target_index + 1,
                position_mode=position_mode,
                source_row=source_row,
                query=query,
                transcript=transcript,
            )
        )
    return annotations


def calls_from_direct_sequences(
    guide: str,
    target_site: str,
    position_mode: str = "huang",
) -> List[MismatchCall]:
    """Score a direct guide/target pair.

    ``guide`` should be 5'->3'. ``target_site`` should be transcript 5'->3';
    guide position 1 is compared with the 3' end of the target site.
    """

    annotations = pair_annotations_from_direct_sequences(guide, target_site, position_mode)
    mismatches, _ = calls_from_annotations(annotations)
    return mismatches


def safe_int(value: object) -> Optional[int]:
    try:
        if value is None or value == "":
            return None
        return int(value)
    except (TypeError, ValueError):
        return None


def summarize_calls(
    mismatches: Sequence[MismatchCall],
    bulges: Sequence[BulgeCall] = (),
) -> Dict[str, object]:
    mismatch_stars = sum(call.total_stars for call in mismatches)
    bulge_star_total = sum(call.stars for call in bulges)
    worst = max([call.total_stars for call in mismatches] or [0])
    worst_calls = [call.compact() for call in mismatches if call.total_stars == worst]
    all_details = [call.compact() for call in mismatches]
    bulge_details = [call.compact() for call in bulges]
    return {
        "MismatchToleranceStars": "*" * mismatch_stars,
        "MismatchToleranceScore": mismatch_stars,
        "MismatchCountScored": len(mismatches),
        "WorstMismatchStars": "*" * worst,
        "WorstMismatchScore": worst,
        "WorstMismatches": ";".join(worst_calls) if worst_calls else "NA",
        "MismatchDetails": ";".join(all_details) if all_details else "NA",
        "BulgeToleranceStars": "*" * bulge_star_total,
        "BulgeToleranceScore": bulge_star_total,
        "BulgeDetails": ";".join(bulge_details) if bulge_details else "NA",
        "TotalToleranceScore": mismatch_stars + bulge_star_total,
        "TwoOrMoreMismatchFlag": "yes" if len(mismatches) >= 2 else "no",
    }


def read_tsv(path: Path) -> Tuple[List[Dict[str, object]], List[str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), list(reader.fieldnames or [])


def read_delimited(path: Path) -> Tuple[List[Dict[str, object]], List[str]]:
    sample = path.read_text()[:4096]
    delimiter = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    try:
        delimiter = csv.Sniffer().sniff(sample, delimiters="\t,;").delimiter
    except csv.Error:
        pass
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        return list(reader), list(reader.fieldnames or [])


def read_json_rows(path: Path) -> Tuple[List[Dict[str, object]], List[str]]:
    data = json.loads(path.read_text())
    if isinstance(data, dict):
        rows = data.get("hits") or data.get("rows") or data.get("results") or []
    else:
        rows = data
    if not isinstance(rows, list):
        raise ValueError(f"Could not find a list of rows in {path}")
    normalized = [dict(row) for row in rows]
    fieldnames: List[str] = []
    for row in normalized:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    return normalized, fieldnames


def write_tsv(rows: Sequence[Dict[str, object]], fieldnames: Sequence[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def normalized_key(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", name.lower())


def resolve_column(
    fieldnames: Sequence[str],
    explicit: Optional[str],
    aliases: Sequence[str],
    label: str,
) -> str:
    if explicit:
        if explicit in fieldnames:
            return explicit
        raise ValueError(f"{label} column {explicit!r} was not found.")
    lookup = {normalized_key(name): name for name in fieldnames}
    for alias in aliases:
        found = lookup.get(normalized_key(alias))
        if found:
            return found
    raise ValueError(
        f"Could not find a {label} column. Tried: {', '.join(aliases)}."
    )


def optional_column(
    fieldnames: Sequence[str],
    explicit: Optional[str],
    aliases: Sequence[str],
) -> Optional[str]:
    if explicit:
        if explicit in fieldnames:
            return explicit
        raise ValueError(f"Column {explicit!r} was not found.")
    lookup = {normalized_key(name): name for name in fieldnames}
    for alias in aliases:
        found = lookup.get(normalized_key(alias))
        if found:
            return found
    return None


PAIR_TABLE_COLUMNS = [
    "SourceRow",
    "Query",
    "Transcript",
    "GuidePosition",
    "TargetPosition",
    "GuideBase",
    "TargetBase",
    "Pair",
    "PairType",
    "PairStars",
    "PositionStars",
    "TotalStars",
    "StarAnnotation",
    "Reason",
]


VERTICAL_TABLE_COLUMNS = [
    "SourceRow",
    "Query",
    "Transcript",
    "GuidePosition",
    "TargetPosition",
    "Annotation",
    "siRNA_guide_5to3",
    "transcript_target_3to5",
    "Pair",
    "PairType",
    "TotalStars",
    "Reason",
]


def sorted_annotations(annotations: Sequence[PairAnnotation]) -> List[PairAnnotation]:
    return sorted(
        annotations,
        key=lambda ann: (
            ann.source_row,
            ann.guide_pos if ann.guide_pos is not None else 10_000,
            ann.target_pos if ann.target_pos is not None else 10_000,
        ),
    )


def vertical_rows_from_annotations(
    annotations: Sequence[PairAnnotation],
) -> List[Dict[str, object]]:
    """Build a stable vertical alignment table, one guide position per row."""

    rows: List[Dict[str, object]] = []
    for ann in sorted_annotations(annotations):
        rows.append(
            {
                "SourceRow": ann.source_row,
                "Query": ann.query,
                "Transcript": ann.transcript,
                "GuidePosition": ann.guide_pos if ann.guide_pos is not None else "",
                "TargetPosition": ann.target_pos if ann.target_pos is not None else "",
                "Annotation": ann.stars,
                "siRNA_guide_5to3": ann.guide_base,
                "transcript_target_3to5": ann.target_base,
                "Pair": ann.pair,
                "PairType": ann.pair_type,
                "TotalStars": ann.total_stars,
                "Reason": ann.reason,
            }
        )
    return rows


def text_table(rows: Sequence[Dict[str, object]], columns: Sequence[str]) -> str:
    if not rows:
        return ""
    widths = {
        column: max(
            len(column),
            max(len(str(row.get(column, ""))) for row in rows),
        )
        for column in columns
    }
    header = "  ".join(column.ljust(widths[column]) for column in columns)
    separator = "  ".join("-" * widths[column] for column in columns)
    body = [
        "  ".join(str(row.get(column, "")).ljust(widths[column]) for column in columns)
        for row in rows
    ]
    return "\n".join([header, separator, *body])


def detailed_text_block(
    annotations: Sequence[PairAnnotation],
    title: str,
    bulge_stars: int = 3,
) -> str:
    annotations = sorted_annotations(annotations)
    mismatches, bulges = calls_from_annotations(annotations, bulge_stars=bulge_stars)
    summary = summarize_calls(mismatches, bulges)
    vertical_rows = vertical_rows_from_annotations(annotations)
    table_columns = [
        "GuidePosition",
        "Annotation",
        "siRNA_guide_5to3",
        "transcript_target_3to5",
        "Pair",
        "PairType",
        "Reason",
    ]
    query = annotations[0].query if annotations else ""
    transcript = annotations[0].transcript if annotations else ""
    source_row = annotations[0].source_row if annotations else ""
    lines = [
        "=" * 96,
        title,
        f"SourceRow: {source_row}",
        f"Query: {query}",
        f"Transcript: {transcript}",
        (
            "Summary: "
            f"total_score={summary['TotalToleranceScore']}; "
            f"mismatch_count={summary['MismatchCountScored']}; "
            f"worst={summary['WorstMismatchStars'] or 'none'}; "
            f"two_or_more_mismatches={summary['TwoOrMoreMismatchFlag']}"
        ),
        f"Worst mismatches: {summary['WorstMismatches']}",
        "",
        text_table(vertical_rows, table_columns),
        "",
        f"Mismatch details: {summary['MismatchDetails']}",
        f"Bulge details: {summary['BulgeDetails']}",
    ]
    return "\n".join(lines)


def write_detailed_report(
    annotation_groups: Sequence[Sequence[PairAnnotation]],
    output_path: Path,
    bulge_stars: int = 3,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    blocks = [
        detailed_text_block(group, title=f"Pair {index}", bulge_stars=bulge_stars)
        for index, group in enumerate(annotation_groups, start=1)
    ]
    output_path.write_text("\n\n".join(blocks) + ("\n" if blocks else ""))


def alignment_rows_from_annotations(
    annotations: Sequence[PairAnnotation],
) -> List[Dict[str, object]]:
    """Build three display rows: stars, guide, and paired transcript bases."""

    rows: List[Dict[str, object]] = []
    grouped: Dict[int, List[PairAnnotation]] = {}
    for ann in annotations:
        grouped.setdefault(ann.source_row, []).append(ann)

    for source_row in sorted(grouped):
        group = sorted_annotations(grouped[source_row])
        if not group:
            continue
        max_pos = max((ann.guide_pos or 0) for ann in group)
        columns = [f"Pos{pos:02d}" for pos in range(1, max_pos + 1)]
        query = group[0].query
        transcript = group[0].transcript
        base = {
            "SourceRow": source_row,
            "Query": query,
            "Transcript": transcript,
        }
        by_pos = {ann.guide_pos: ann for ann in group if ann.guide_pos is not None}

        annotation_row = {**base, "RowType": "Annotation"}
        guide_row = {**base, "RowType": "siRNA_guide_5to3"}
        target_row = {**base, "RowType": "transcript_target_3to5"}

        for pos, column in enumerate(columns, start=1):
            ann = by_pos.get(pos)
            annotation_row[column] = ann.stars if ann else ""
            guide_row[column] = ann.guide_base if ann else ""
            target_row[column] = ann.target_base if ann else ""

        rows.extend([annotation_row, guide_row, target_row])
    return rows


def alignment_fieldnames(rows: Sequence[Dict[str, object]]) -> List[str]:
    fields = ["SourceRow", "Query", "Transcript", "RowType"]
    max_pos = 0
    for row in rows:
        for key in row:
            if key.startswith("Pos"):
                max_pos = max(max_pos, int(key[3:]))
    fields.extend(f"Pos{pos:02d}" for pos in range(1, max_pos + 1))
    return fields


def output_fieldnames(input_fieldnames: Sequence[str]) -> List[str]:
    added = [
        "MismatchToleranceStars",
        "MismatchToleranceScore",
        "MismatchCountScored",
        "WorstMismatchStars",
        "WorstMismatchScore",
        "WorstMismatches",
        "MismatchDetails",
        "BulgeToleranceStars",
        "BulgeToleranceScore",
        "BulgeDetails",
        "TotalToleranceScore",
        "TwoOrMoreMismatchFlag",
    ]
    fields = list(input_fieldnames)
    for name in added:
        if name not in fields:
            fields.append(name)
    return fields


def score_file(
    input_path: Path,
    output_path: Path,
    pair_table_output_path: Optional[Path],
    alignment_table_output_path: Optional[Path],
    vertical_table_output_path: Optional[Path],
    position_mode: str,
    include_bulges: bool,
    bulge_stars: int,
) -> int:
    if input_path.suffix.lower() == ".json":
        rows, fields = read_json_rows(input_path)
    else:
        rows, fields = read_tsv(input_path)

    out_rows: List[Dict[str, object]] = []
    pair_rows: List[Dict[str, object]] = []
    all_annotations: List[PairAnnotation] = []
    for row_index, row in enumerate(rows, start=1):
        try:
            annotations = pair_annotations_from_mfe_row(
                row,
                source_row=row_index,
                position_mode=position_mode,
                include_bulges=include_bulges,
                bulge_stars=bulge_stars,
            )
            mismatches, bulges = calls_from_annotations(annotations, bulge_stars=bulge_stars)
            summary = summarize_calls(mismatches, bulges)
            pair_rows.extend(annotation.as_row() for annotation in annotations)
            all_annotations.extend(annotations)
        except Exception as exc:
            summary = {
                "MismatchToleranceStars": "",
                "MismatchToleranceScore": "",
                "MismatchCountScored": "",
                "WorstMismatchStars": "",
                "WorstMismatchScore": "",
                "WorstMismatches": "ERROR",
                "MismatchDetails": f"row_{row_index}: {exc}",
                "BulgeToleranceStars": "",
                "BulgeToleranceScore": "",
                "BulgeDetails": "",
                "TotalToleranceScore": "",
                "TwoOrMoreMismatchFlag": "",
            }
        merged = dict(row)
        merged.update(summary)
        out_rows.append(merged)

    write_tsv(out_rows, output_fieldnames(fields), output_path)
    if pair_table_output_path is not None:
        pair_rows = [annotation.as_row() for annotation in sorted_annotations(all_annotations)]
        write_tsv(pair_rows, PAIR_TABLE_COLUMNS, pair_table_output_path)
    if alignment_table_output_path is not None:
        alignment_rows = alignment_rows_from_annotations(all_annotations)
        write_tsv(alignment_rows, alignment_fieldnames(alignment_rows), alignment_table_output_path)
    if vertical_table_output_path is not None:
        vertical_rows = vertical_rows_from_annotations(all_annotations)
        write_tsv(vertical_rows, VERTICAL_TABLE_COLUMNS, vertical_table_output_path)
    return len(out_rows)


def score_bulk_direct_file(
    bulk_input_path: Path,
    report_output_path: Path,
    position_mode: str,
    guide_column: Optional[str],
    target_column: Optional[str],
    query_column: Optional[str],
    transcript_column: Optional[str],
    bulge_stars: int,
) -> int:
    rows, fieldnames = read_delimited(bulk_input_path)
    guide_col = resolve_column(
        fieldnames,
        guide_column,
        ["guide", "sirna", "sirna_sequence", "guide_sequence", "query_sequence", "srna"],
        "guide",
    )
    target_col = resolve_column(
        fieldnames,
        target_column,
        [
            "target_site",
            "target",
            "target_sequence",
            "transcript_target",
            "transcript_region",
            "transcript_sequence",
        ],
        "target site",
    )
    query_col = optional_column(
        fieldnames,
        query_column,
        ["query", "query_id", "sirna_id", "guide_id", "guide_name", "name"],
    )
    transcript_col = optional_column(
        fieldnames,
        transcript_column,
        ["transcript", "transcript_id", "target_id", "target_name", "transcript_name"],
    )

    annotation_groups: List[List[PairAnnotation]] = []
    for row_index, row in enumerate(rows, start=1):
        guide = str(row.get(guide_col, "")).strip()
        target = str(row.get(target_col, "")).strip()
        if not guide or not target:
            raise ValueError(f"Bulk row {row_index} is missing guide or target sequence.")
        query = str(row.get(query_col, "")).strip() if query_col else f"guide_{row_index}"
        transcript = (
            str(row.get(transcript_col, "")).strip()
            if transcript_col
            else f"target_{row_index}"
        )
        if not query:
            query = f"guide_{row_index}"
        if not transcript:
            transcript = f"target_{row_index}"
        annotations = pair_annotations_from_direct_sequences(
            guide,
            target,
            position_mode=position_mode,
            source_row=row_index,
            query=query,
            transcript=transcript,
        )
        annotation_groups.append(annotations)

    write_detailed_report(annotation_groups, report_output_path, bulge_stars=bulge_stars)
    return len(annotation_groups)


def print_direct_score(
    guide: str,
    target_site: str,
    position_mode: str,
    long_table: bool,
    three_row: bool,
    detail_output: Optional[Path] = None,
) -> None:
    annotations = pair_annotations_from_direct_sequences(
        guide, target_site, position_mode=position_mode
    )
    if detail_output is not None:
        write_detailed_report([annotations], detail_output)
    if not long_table and not three_row:
        writer = csv.DictWriter(
            __import__("sys").stdout,
            fieldnames=VERTICAL_TABLE_COLUMNS,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(vertical_rows_from_annotations(annotations))
        return
    if not long_table:
        rows = alignment_rows_from_annotations(annotations)
        writer = csv.DictWriter(
            __import__("sys").stdout,
            fieldnames=alignment_fieldnames(rows),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)
        return
    writer = csv.DictWriter(
        __import__("sys").stdout,
        fieldnames=PAIR_TABLE_COLUMNS,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(annotation.as_row() for annotation in annotations)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Post-process MFE_ratio.py TSV/JSON outputs and score siRNA guide-target "
            "mismatches by nucleotide-pair tolerance plus position sensitivity."
        )
    )
    parser.add_argument("--input", type=Path, help="MFE_ratio.py TSV or JSON output.")
    parser.add_argument(
        "--bulk-input",
        type=Path,
        help=(
            "TSV/CSV with many direct guide-target pairs. Expected columns include "
            "guide and target_site; column names can be overridden."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output TSV. Defaults to <input>.mismatch_tolerance.tsv.",
    )
    parser.add_argument(
        "--bulk-report-output",
        type=Path,
        help="Combined detailed .txt report for --bulk-input. Defaults to <bulk-input>.detailed.txt.",
    )
    parser.add_argument("--bulk-guide-column", help="Column name for guide siRNA sequence.")
    parser.add_argument("--bulk-target-column", help="Column name for transcript target-site sequence.")
    parser.add_argument("--bulk-query-column", help="Optional column name for guide/query ID.")
    parser.add_argument("--bulk-transcript-column", help="Optional column name for transcript/target ID.")
    parser.add_argument(
        "--pair-table-output",
        type=Path,
        help="Long pair-annotation TSV. Defaults to <input>.pair_annotations.tsv.",
    )
    parser.add_argument(
        "--alignment-table-output",
        type=Path,
        help="Three-row alignment TSV. Defaults to <input>.alignment_annotations.tsv.",
    )
    parser.add_argument(
        "--vertical-table-output",
        type=Path,
        help="Vertical alignment TSV. Defaults to <input>.vertical_annotations.tsv.",
    )
    parser.add_argument(
        "--no-pair-table",
        action="store_true",
        help="Only write the row-level summary TSV.",
    )
    parser.add_argument(
        "--no-alignment-table",
        action="store_true",
        help="Do not write the three-row alignment TSV.",
    )
    parser.add_argument(
        "--no-vertical-table",
        action="store_true",
        help="Do not write the vertical alignment TSV.",
    )
    parser.add_argument("--guide", help="Direct mode: guide siRNA sequence, 5'->3'.")
    parser.add_argument(
        "--target-site",
        help="Direct mode: transcript target site sequence, 5'->3'. Must match guide length.",
    )
    parser.add_argument(
        "--direct-long-table",
        action="store_true",
        help="In direct mode, print the detailed one-row-per-pair table instead of the compact vertical table.",
    )
    parser.add_argument(
        "--direct-three-row",
        action="store_true",
        help="In direct mode, print the older three-row horizontal alignment.",
    )
    parser.add_argument(
        "--detail-output",
        type=Path,
        help="Write a detailed .txt report for the single direct --guide/--target-site pair.",
    )
    parser.add_argument(
        "--position-mode",
        choices=["huang", "general", "none"],
        default="huang",
        help="How to add the extra position star. Default: Huang 2009 guide-base-specific positions.",
    )
    parser.add_argument(
        "--ignore-bulges",
        action="store_true",
        help="Do not add separate bulge/gap annotations.",
    )
    parser.add_argument(
        "--bulge-stars",
        type=int,
        default=3,
        help="Stars assigned to each bulge/gap annotation. Default: 3.",
    )
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    if args.bulk_input:
        if args.guide or args.target_site or args.input:
            parser.error("Use --bulk-input by itself, not with --guide/--target-site or --input.")
        report_output = args.bulk_report_output
        if report_output is None:
            report_output = args.bulk_input.with_suffix(args.bulk_input.suffix + ".detailed.txt")
        count = score_bulk_direct_file(
            bulk_input_path=args.bulk_input,
            report_output_path=report_output,
            position_mode=args.position_mode,
            guide_column=args.bulk_guide_column,
            target_column=args.bulk_target_column,
            query_column=args.bulk_query_column,
            transcript_column=args.bulk_transcript_column,
            bulge_stars=args.bulge_stars,
        )
        print(f"Wrote detailed reports for {count} guide-target pairs to {report_output}")
        return

    if args.guide or args.target_site:
        if not args.guide or not args.target_site:
            parser.error("--guide and --target-site must be supplied together.")
        print_direct_score(
            args.guide,
            args.target_site,
            args.position_mode,
            long_table=args.direct_long_table,
            three_row=args.direct_three_row,
            detail_output=args.detail_output,
        )
        return

    if not args.input:
        parser.error("Provide --input, or use direct mode with --guide and --target-site.")

    output = args.output
    if output is None:
        output = args.input.with_suffix(args.input.suffix + ".mismatch_tolerance.tsv")
    pair_table_output = None
    if not args.no_pair_table:
        pair_table_output = args.pair_table_output
        if pair_table_output is None:
            pair_table_output = args.input.with_suffix(args.input.suffix + ".pair_annotations.tsv")
    alignment_table_output = None
    if not args.no_alignment_table:
        alignment_table_output = args.alignment_table_output
        if alignment_table_output is None:
            alignment_table_output = args.input.with_suffix(
                args.input.suffix + ".alignment_annotations.tsv"
            )
    vertical_table_output = None
    if not args.no_vertical_table:
        vertical_table_output = args.vertical_table_output
        if vertical_table_output is None:
            vertical_table_output = args.input.with_suffix(
                args.input.suffix + ".vertical_annotations.tsv"
            )

    count = score_file(
        input_path=args.input,
        output_path=output,
        pair_table_output_path=pair_table_output,
        alignment_table_output_path=alignment_table_output,
        vertical_table_output_path=vertical_table_output,
        position_mode=args.position_mode,
        include_bulges=not args.ignore_bulges,
        bulge_stars=args.bulge_stars,
    )
    print(f"Wrote mismatch tolerance scores for {count} rows to {output}")
    if pair_table_output is not None:
        print(f"Wrote per-pair annotations to {pair_table_output}")
    if alignment_table_output is not None:
        print(f"Wrote three-row alignment annotations to {alignment_table_output}")
    if vertical_table_output is not None:
        print(f"Wrote vertical alignment annotations to {vertical_table_output}")


if __name__ == "__main__":
    main()
