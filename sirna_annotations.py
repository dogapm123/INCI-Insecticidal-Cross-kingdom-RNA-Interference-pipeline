"""Exact-match siRNA annotation helpers for positional coverage plots.

The annotation FASTA is deliberately independent from read mapping and scoring.
Each record is matched to a plotted reference in both orientations so its given
5-prime to 3-prime sequence direction can be shown accurately on the plot.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Mapping


ANNOTATION_COLORS = (
    "#0F766E",
    "#7C3AED",
    "#C2410C",
    "#BE123C",
    "#0369A1",
    "#4D7C0F",
    "#A21CAF",
)


def reverse_complement(sequence: str) -> str:
    return sequence.upper().replace("U", "T").translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name = ""
    chunks: list[str] = []
    with path.open(encoding="utf-8-sig", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    records.append((name, "".join(chunks)))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name:
        records.append((name, "".join(chunks)))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def find_sirna_annotations(annotation_fasta: Path, references: Mapping[str, str]) -> list[dict[str, Any]]:
    """Find every exact forward or reverse-complement match in each reference."""
    if not annotation_fasta.exists() or not annotation_fasta.is_file():
        raise ValueError(f"siRNA annotation FASTA does not exist: {annotation_fasta}")
    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, str, int, int, str]] = set()
    for sirna_id, raw_sequence in read_fasta(annotation_fasta):
        sequence = "".join(raw_sequence.upper().replace("U", "T").split())
        if not sequence or any(base not in "ACGTN" for base in sequence):
            raise ValueError(f"siRNA {sirna_id!r} contains unsupported sequence characters.")
        orientations = [("forward", sequence)]
        reverse = reverse_complement(sequence)
        if reverse != sequence:
            orientations.append(("reverse", reverse))
        for contig, reference in references.items():
            target = reference.upper().replace("U", "T")
            for direction, query in orientations:
                start0 = target.find(query)
                while start0 >= 0:
                    key = (sirna_id, contig, start0 + 1, start0 + len(query), direction)
                    if key not in seen:
                        seen.add(key)
                        rows.append(
                            {
                                "siRNA_id": sirna_id,
                                "siRNA_sequence": sequence,
                                "contig": contig,
                                "start_1based": start0 + 1,
                                "end_1based": start0 + len(query),
                                "direction_on_reference": direction,
                            }
                        )
                    start0 = target.find(query, start0 + 1)
    return sorted(rows, key=lambda row: (str(row["contig"]), int(row["start_1based"]), str(row["siRNA_id"])))


def write_sirna_annotation_table(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["siRNA_id", "siRNA_sequence", "contig", "start_1based", "end_1based", "direction_on_reference"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def draw_sirna_annotation_track(
    ax: Any,
    annotations: list[dict[str, Any]],
    contig: str,
    length: int,
    *,
    label_size: float = 8.0,
    position_offset: int = 0,
) -> int:
    """Overlay labelled siRNA intervals and orientation arrows on one axes."""
    selected = [row for row in annotations if str(row.get("contig", "")) == contig]
    if not selected:
        return 0
    y_top = float(ax.get_ylim()[1])
    if y_top <= 0:
        return 0
    label_offset = max(18.0, min(length * 0.04, 160.0))
    for index, row in enumerate(selected):
        start = max(1, int(row["start_1based"]))
        end = min(length, int(row["end_1based"]))
        if end < start:
            continue
        plot_start = start + position_offset
        plot_end = end + position_offset
        color = ANNOTATION_COLORS[index % len(ANNOTATION_COLORS)]
        y_fraction = 0.965 - 0.105 * (index % 4)
        arrow_y = y_top * y_fraction
        midpoint = (start + end) / 2.0
        offset = label_offset if index % 2 == 0 else -label_offset
        label_x = min(max(midpoint + position_offset + offset, float(1 + position_offset)), float(length + position_offset))
        ax.axvspan(plot_start - 0.5, plot_end + 0.5, color=color, alpha=0.24, zorder=5)
        if row["direction_on_reference"] == "forward":
            arrow_start, arrow_end = plot_start, plot_end
        else:
            arrow_start, arrow_end = plot_end, plot_start
        ax.annotate(
            "",
            xy=(arrow_end, arrow_y),
            xytext=(arrow_start, arrow_y),
            arrowprops={"arrowstyle": "->", "color": color, "linewidth": 1.2},
            zorder=7,
        )
        ax.annotate(
            str(row["siRNA_id"]),
            xy=(midpoint + position_offset, arrow_y),
            xytext=(label_x, arrow_y + y_top * 0.035),
            color=color,
            fontsize=label_size,
            fontweight="bold",
            ha="center",
            va="bottom",
            arrowprops={"arrowstyle": "-", "color": color, "linewidth": 0.8},
            annotation_clip=True,
            zorder=7,
        )
    return len(selected)
