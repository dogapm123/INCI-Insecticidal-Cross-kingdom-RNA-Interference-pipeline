"""Summarize per-siRNA hits and draw strand-aware multi-transcriptome landscapes."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from xml.sax.saxutils import escape

from degradome_analysis import write_rows
from dsrna_adaptation import normalized_shannon_entropy, safe_name


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def entropy_colour(value: float) -> str:
    low, high = (198, 38, 38), (38, 140, 50)
    channels = [round(a + (b - a) * value) for a, b in zip(low, high)]
    return "#" + "".join(f"{channel:02x}" for channel in channels)


def diverging_colour(value: float, limit: float) -> str:
    if limit <= 0:
        return "#f5f5f5"
    fraction = min(1.0, abs(value) / limit)
    base = (53, 112, 180) if value >= 0 else (196, 75, 75)
    return "#" + "".join(f"{round(245 + (channel - 245) * fraction):02x}" for channel in base)


def landscape_svg(rows: list[dict[str, object]], path: Path, dsrna_id: str) -> None:
    transcriptomes = list(dict.fromkeys(str(row["transcriptome"]) for row in rows))
    positions = [int(row["dsrna_position"]) for row in rows]
    values = [abs(float(row["signed_relative_adaptation_excess_density_per_million_sites"])) for row in rows]
    x_min, x_max = min(positions), max(positions)
    y_max = max(1.0, max(values) * 1.15)
    width, panel_height, left, right, top, bottom = 1100, 235, 145, 35, 36, 52
    height = top + panel_height * len(transcriptomes) + bottom
    plot_width = width - left - right
    def x(value: int) -> float:
        return left + plot_width / 2 if x_min == x_max else left + (value - x_min) / (x_max - x_min) * plot_width
    def y(value: float, panel: int) -> float:
        middle = top + panel * panel_height + (panel_height - 30) / 2
        return middle - value / y_max * ((panel_height - 30) / 2 - 12)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="Strand-aware dsRNA adaptation landscape for {escape(dsrna_id)}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{left}" y="22" font-family="Arial" font-size="18" font-weight="bold">Relative dsRNA adaptation landscape: {escape(dsrna_id)}</text>',
        '<defs><linearGradient id="entropy-gradient"><stop offset="0%" stop-color="#c62626"/><stop offset="100%" stop-color="#268c32"/></linearGradient></defs>',
    ]
    for panel, transcriptome in enumerate(transcriptomes):
        y_zero = y(0.0, panel)
        frame_top = top + panel * panel_height
        parts.extend([
            f'<rect x="{left}" y="{frame_top}" width="{plot_width}" height="{panel_height - 30}" fill="none" stroke="#777" stroke-width="1"/>',
            f'<line x1="{left}" y1="{y_zero:.2f}" x2="{width - right}" y2="{y_zero:.2f}" stroke="#333" stroke-width="1"/>',
            f'<text x="{left + 8}" y="{frame_top + 22}" font-family="Arial" font-size="15" font-weight="bold">{escape(transcriptome)}</text>',
            f'<text x="{left + 8}" y="{frame_top + 42}" font-family="Arial" font-size="11">+ sense / − antisense; relative to other transcriptomes</text>',
            f'<text x="{left - 9}" y="{y( y_max, panel) + 4:.2f}" text-anchor="end" font-family="Arial" font-size="11">+{y_max:.2f}</text>',
            f'<text x="{left - 9}" y="{y(-y_max, panel) + 4:.2f}" text-anchor="end" font-family="Arial" font-size="11">−{y_max:.2f}</text>',
            f'<text x="{left - 9}" y="{y_zero + 4:.2f}" text-anchor="end" font-family="Arial" font-size="11">0</text>',
        ])
        for row in rows:
            if row["transcriptome"] != transcriptome:
                continue
            position = int(row["dsrna_position"])
            score = float(row["signed_relative_adaptation_excess_density_per_million_sites"])
            entropy = float(row["shannon_entropy"])
            xx, yy = x(position), y(score, panel)
            title = f"{row['siRNA_id']} | position {position} | {row['strand']} | relative excess {score:.3f} | entropy {entropy:.3f}"
            parts.extend([
                f'<line x1="{xx:.2f}" y1="{y_zero:.2f}" x2="{xx:.2f}" y2="{yy:.2f}" stroke="#555" stroke-width="1"/>',
                f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="5" fill="{entropy_colour(entropy)}" stroke="#222" stroke-width="0.7"><title>{escape(title)}</title></circle>',
            ])
    x_axis = top + panel_height * len(transcriptomes) - 20
    tick_count = 6
    ticks = [round(x_min + (x_max - x_min) * index / (tick_count - 1)) for index in range(tick_count)]
    for tick in ticks:
        parts.append(f'<text x="{x(tick):.2f}" y="{x_axis + 18}" text-anchor="middle" font-family="Arial" font-size="10">{tick}</text>')
    parts.extend([
        f'<text x="{left + plot_width / 2:.2f}" y="{height - 8}" text-anchor="middle" font-family="Arial" font-size="13">siRNA position in supplied dsRNA locus; y = relative adaptation excess (sites/million windows)</text>',
        f'<text x="{width - 245}" y="22" font-family="Arial" font-size="11">Entropy</text>',
        f'<rect x="{width - 195}" y="12" width="105" height="10" fill="url(#entropy-gradient)"/>',
        f'<text x="{width - 200}" y="34" font-family="Arial" font-size="10">low</text><text x="{width - 85}" y="34" font-family="Arial" font-size="10">high</text>',
        '</svg>',
    ])
    path.write_text("\n".join(parts), encoding="utf-8")


def overall_summary_svg(rows: list[dict[str, object]], ranking: list[dict[str, str]], path: Path) -> list[dict[str, object]]:
    transcriptomes = [str(row["transcriptome"]) for row in ranking]
    sirnas = sorted({str(row["siRNA_id"]) for row in rows}, key=lambda value: next(int(row["dsrna_position"]) for row in rows if row["siRNA_id"] == value))
    by_sirna: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_sirna[str(row["siRNA_id"])].append(row)
    variability_rows: list[dict[str, object]] = []
    for siRNA_id in sirnas:
        values = [float(row["relative_adaptation_excess_density_per_million_sites"]) for row in by_sirna[siRNA_id]]
        reference = by_sirna[siRNA_id][0]
        variability_rows.append({
            "siRNA_id": siRNA_id, "dsrna_position": reference["dsrna_position"], "strand": reference["strand"],
            "siRNA_sequence": reference["siRNA_sequence"], "shannon_entropy": reference["shannon_entropy"],
            "relative_adaptation_range_per_million_sites": max(values) - min(values),
            "most_adapted_transcriptome": max(by_sirna[siRNA_id], key=lambda row: float(row["relative_adaptation_excess_density_per_million_sites"]))["transcriptome"],
        })

    width, height, left, right = 1200, 880, 150, 55
    chart_width = width - left - right
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="Overall siRNA adaptation summaries">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="150" y="28" font-family="Arial" font-size="20" font-weight="bold">Overall siRNA adaptation summaries</text>',
    ]
    # Panel 1: transcriptome-level adaptation excess.
    p1_top, p1_height = 55, 190
    excesses = [float(row["adaptation_excess_density_per_million_sites"]) for row in ranking]
    limit = max(0.5, max(abs(value) for value in excesses) * 1.2)
    zero_x = left + chart_width / 2
    parts.extend([
        f'<text x="{left}" y="{p1_top - 10}" font-family="Arial" font-size="15" font-weight="bold">Transcriptome-level adaptation excess</text>',
        f'<text x="{left}" y="{p1_top + p1_height + 18}" font-family="Arial" font-size="12">Original density minus matched-shuffle density (sites/million windows)</text>',
        f'<line x1="{zero_x:.2f}" y1="{p1_top}" x2="{zero_x:.2f}" y2="{p1_top + p1_height}" stroke="#333" stroke-width="1"/>',
    ])
    for index, row in enumerate(ranking):
        value = float(row["adaptation_excess_density_per_million_sites"])
        y = p1_top + 28 + index * 52
        x_end = zero_x + value / limit * (chart_width / 2 - 55)
        x_rect, bar_width = (zero_x, x_end - zero_x) if value >= 0 else (x_end, zero_x - x_end)
        parts.extend([
            f'<text x="{left - 12}" y="{y + 5}" text-anchor="end" font-family="Arial" font-size="13">{escape(row["transcriptome"])}</text>',
            f'<rect x="{x_rect:.2f}" y="{y - 12}" width="{bar_width:.2f}" height="24" fill="{diverging_colour(value, limit)}"/>',
            f'<text x="{x_end + (8 if value >= 0 else -8):.2f}" y="{y + 5}" text-anchor="{"start" if value >= 0 else "end"}" font-family="Arial" font-size="12">{value:+.3f}; P={float(row["pvalue_vs_matched_shuffles"]):.3g}</text>',
        ])

    # Panel 2: siRNA-by-transcriptome relative-adaptation heatmap.
    p2_top, p2_height, cell_height = 325, 175, 34
    relative_limit = max(0.25, max(abs(float(row["relative_adaptation_excess_density_per_million_sites"])) for row in rows))
    cell_width = chart_width / len(sirnas)
    parts.extend([
        f'<text x="{left}" y="{p2_top - 10}" font-family="Arial" font-size="15" font-weight="bold">Per-siRNA relative adaptation</text>',
        f'<text x="{left}" y="{p2_top + p2_height + 30}" font-family="Arial" font-size="12">siRNAs ordered by supplied position; blue = higher relative adaptation, red = lower</text>',
        f'<rect x="{width - 212}" y="{p2_top - 22}" width="50" height="10" fill="#c44b4b"/><rect x="{width - 162}" y="{p2_top - 22}" width="50" height="10" fill="#f5f5f5"/><rect x="{width - 112}" y="{p2_top - 22}" width="50" height="10" fill="#3570b4"/>',
        f'<text x="{width - 216}" y="{p2_top - 27}" font-family="Arial" font-size="10">−</text><text x="{width - 65}" y="{p2_top - 27}" font-family="Arial" font-size="10">+</text>',
    ])
    lookup = {(str(row["transcriptome"]), str(row["siRNA_id"])): row for row in rows}
    for row_index, transcriptome in enumerate(transcriptomes):
        y = p2_top + row_index * cell_height
        parts.append(f'<text x="{left - 12}" y="{y + 21}" text-anchor="end" font-family="Arial" font-size="12">{escape(transcriptome)}</text>')
        for column, siRNA_id in enumerate(sirnas):
            value = float(lookup[(transcriptome, siRNA_id)]["relative_adaptation_excess_density_per_million_sites"])
            xx = left + column * cell_width
            title = f"{siRNA_id} | {transcriptome} | relative adaptation {value:.3f}"
            parts.append(f'<rect x="{xx:.2f}" y="{y}" width="{cell_width:.2f}" height="{cell_height - 2}" fill="{diverging_colour(value, relative_limit)}"><title>{escape(title)}</title></rect>')
    for index, siRNA_id in enumerate(sirnas):
        if index % 4 == 0 or index == len(sirnas) - 1:
            position = next(row["dsrna_position"] for row in rows if row["siRNA_id"] == siRNA_id)
            xx = left + (index + 0.5) * cell_width
            parts.append(f'<text x="{xx:.2f}" y="{p2_top + p2_height + 15}" text-anchor="middle" font-family="Arial" font-size="10">{position}</text>')

    # Panel 3: sequence complexity versus transcriptome-specific behaviour.
    p3_top, p3_height = 610, 210
    entropies = [float(row["shannon_entropy"]) for row in variability_rows]
    ranges = [float(row["relative_adaptation_range_per_million_sites"]) for row in variability_rows]
    x_min, x_max = min(entropies) - 0.02, min(1.0, max(entropies) + 0.02)
    y_max = max(0.25, max(ranges) * 1.15)
    def sx(value: float) -> float:
        return left + (value - x_min) / (x_max - x_min) * chart_width
    def sy(value: float) -> float:
        return p3_top + p3_height - value / y_max * (p3_height - 35)
    parts.extend([
        f'<text x="{left}" y="{p3_top - 10}" font-family="Arial" font-size="15" font-weight="bold">Sequence complexity versus transcriptome-specificity</text>',
        f'<rect x="{left}" y="{p3_top}" width="{chart_width}" height="{p3_height - 35}" fill="none" stroke="#777" stroke-width="1"/>',
        f'<text x="{left + chart_width / 2}" y="{p3_top + p3_height}" text-anchor="middle" font-family="Arial" font-size="12">Normalized Shannon entropy</text>',
        f'<text x="{left - 14}" y="{p3_top + 12}" text-anchor="end" font-family="Arial" font-size="11">{y_max:.2f}</text><text x="{left - 14}" y="{p3_top + p3_height - 35}" text-anchor="end" font-family="Arial" font-size="11">0</text>',
        f'<text x="24" y="{p3_top + (p3_height - 35) / 2}" transform="rotate(-90 24 {p3_top + (p3_height - 35) / 2})" text-anchor="middle" font-family="Arial" font-size="12">Relative-adaptation range</text>',
    ])
    for tick in (x_min, (x_min + x_max) / 2, x_max):
        parts.append(f'<text x="{sx(tick):.2f}" y="{p3_top + p3_height - 18}" text-anchor="middle" font-family="Arial" font-size="11">{tick:.2f}</text>')
    for row in variability_rows:
        xx, yy = sx(float(row["shannon_entropy"])), sy(float(row["relative_adaptation_range_per_million_sites"]))
        title = f"{row['siRNA_id']} | entropy {float(row['shannon_entropy']):.3f} | relative-adaptation range {float(row['relative_adaptation_range_per_million_sites']):.3f} | most adapted: {row['most_adapted_transcriptome']}"
        shape = '<circle cx="{x:.2f}" cy="{y:.2f}" r="5" fill="{colour}" stroke="#222" stroke-width="0.7">' if row["strand"] == "sense" else '<rect x="{x0:.2f}" y="{y0:.2f}" width="10" height="10" fill="{colour}" stroke="#222" stroke-width="0.7">'
        if row["strand"] == "sense":
            parts.append(shape.format(x=xx, y=yy, colour=entropy_colour(float(row["shannon_entropy"]))) + f'<title>{escape(title)}</title></circle>')
        else:
            parts.append(shape.format(x0=xx - 5, y0=yy - 5, colour=entropy_colour(float(row["shannon_entropy"]))) + f'<title>{escape(title)}</title></rect>')
    parts.extend(['</svg>'])
    path.write_text("\n".join(parts), encoding="utf-8")
    return variability_rows


def adaptability_complexity_svg(rows: list[dict[str, object]], path: Path) -> None:
    """Draw one adaptability-versus-complexity panel per transcriptome."""
    transcriptomes = list(dict.fromkeys(str(row["transcriptome"]) for row in rows))
    width, panel_height, left, right, top, bottom = 1100, 280, 105, 35, 35, 55
    height = top + panel_height * len(transcriptomes) + bottom
    chart_width = width - left - right
    values = [float(row["relative_adaptation_excess_density_per_million_sites"]) for row in rows]
    x_limit = max(0.25, max(abs(value) for value in values) * 1.15)
    entropies = [float(row["shannon_entropy"]) for row in rows]
    y_min, y_max = min(entropies) - 0.025, min(1.0, max(entropies) + 0.015)
    def x(value: float) -> float:
        return left + (value + x_limit) / (2 * x_limit) * chart_width
    def y(value: float, panel: int) -> float:
        return top + panel * panel_height + 20 + (y_max - value) / (y_max - y_min) * (panel_height - 75)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="siRNA adaptability and sequence-complexity plots">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="105" y="22" font-family="Arial" font-size="18" font-weight="bold">siRNA adaptability–complexity plots</text>',
    ]
    for panel, transcriptome in enumerate(transcriptomes):
        frame_top = top + panel * panel_height
        zero_x = x(0.0)
        parts.extend([
            f'<rect x="{left}" y="{frame_top}" width="{chart_width}" height="{panel_height - 50}" fill="none" stroke="#777" stroke-width="1"/>',
            f'<line x1="{zero_x:.2f}" y1="{frame_top}" x2="{zero_x:.2f}" y2="{frame_top + panel_height - 50}" stroke="#333" stroke-width="1"/>',
            f'<text x="{left + 8}" y="{frame_top + 20}" font-family="Arial" font-size="14" font-weight="bold">{escape(transcriptome)}</text>',
            f'<text x="{left - 10}" y="{y(y_max, panel) + 4:.2f}" text-anchor="end" font-family="Arial" font-size="11">{y_max:.2f}</text>',
            f'<text x="{left - 10}" y="{y(y_min, panel) + 4:.2f}" text-anchor="end" font-family="Arial" font-size="11">{y_min:.2f}</text>',
            f'<text x="25" y="{frame_top + (panel_height - 50) / 2}" transform="rotate(-90 25 {frame_top + (panel_height - 50) / 2})" text-anchor="middle" font-family="Arial" font-size="12">Shannon entropy</text>',
            f'<text x="{left}" y="{frame_top + panel_height - 24}" font-family="Arial" font-size="10">− lower relative adaptation</text>',
            f'<text x="{left + chart_width}" y="{frame_top + panel_height - 24}" text-anchor="end" font-family="Arial" font-size="10">higher relative adaptation +</text>',
        ])
        for row in rows:
            if row["transcriptome"] != transcriptome:
                continue
            score = float(row["relative_adaptation_excess_density_per_million_sites"])
            entropy = float(row["shannon_entropy"])
            xx, yy = x(score), y(entropy, panel)
            title = f"{row['siRNA_id']} | relative adaptation {score:.3f} | entropy {entropy:.3f} | {row['strand']}"
            if row["strand"] == "sense":
                parts.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="5" fill="{entropy_colour(entropy)}" stroke="#222" stroke-width="0.7"><title>{escape(title)}</title></circle>')
            else:
                parts.append(f'<rect x="{xx - 5:.2f}" y="{yy - 5:.2f}" width="10" height="10" fill="{entropy_colour(entropy)}" stroke="#222" stroke-width="0.7"><title>{escape(title)}</title></rect>')
    parts.append('</svg>')
    path.write_text("\n".join(parts), encoding="utf-8")


def alignment_quality_summary(
    hit_rows: list[dict[str, str]], per_sirna_rows: list[dict[str, str]], ranking: list[dict[str, str]], path: Path,
) -> list[dict[str, object]]:
    """Summarize mismatch/wobble profiles for originals versus matched shuffles."""
    transcriptomes = [row["transcriptome"] for row in ranking]
    fallback_search_space = {row["transcriptome"]: int(row["searchable_guide_windows"] or 1) for row in ranking}
    original_count = {name: 0 for name in transcriptomes}
    shuffle_count = {name: 0 for name in transcriptomes}
    for row in per_sirna_rows:
        name = row["transcriptome"]
        original_count[name] += 1
        shuffle_count[name] += len(row["shuffled_densities_per_million_sites"].split(","))
    counts: dict[tuple[str, str, int, int], int] = defaultdict(int)
    density_sums: dict[tuple[str, str, int, int], float] = defaultdict(float)
    for row in hit_rows:
        key = (row["transcriptome"], row["variant"], int(row["generic_mismatch_count"]), int(row["wobble_count"]))
        counts[key] += 1
        density_sums[key] += 1_000_000.0 / int(row.get("searchable_guide_windows") or fallback_search_space[row["transcriptome"]])
    summary: list[dict[str, object]] = []
    for transcriptome in transcriptomes:
        for variant, guide_count in (("original", original_count[transcriptome]), ("shuffle", shuffle_count[transcriptome])):
            for generic_mismatches in range(4):
                for wobble_count in sorted({0, 1} | {key[3] for key in counts if key[0] == transcriptome}):
                    hit_count = counts[(transcriptome, variant, generic_mismatches, wobble_count)]
                    summary.append({
                        "transcriptome": transcriptome,
                        "guide_group": "original" if variant == "original" else "mean_matched_shuffle",
                        "generic_mismatch_count": generic_mismatches,
                        "wobble_count": wobble_count,
                        "alignment_count": hit_count,
                        "guide_count": guide_count,
                        "mean_alignments_per_siRNA": hit_count / max(1, guide_count),
                        "mean_alignment_density_per_million_sites_per_siRNA": density_sums[(transcriptome, variant, generic_mismatches, wobble_count)] / max(1, guide_count),
                    })

    categories = [(mismatch, wobble) for mismatch in range(4) for wobble in (0, 1)]
    lookup = {(str(row["transcriptome"]), str(row["guide_group"]), int(row["generic_mismatch_count"]), int(row["wobble_count"])): row for row in summary}
    values = [float(row["mean_alignment_density_per_million_sites_per_siRNA"]) for row in summary]
    limit = max(0.01, max(values) * 1.15)
    width, left, right, top, row_height = 1280, 170, 55, 55, 168
    height = top + row_height * len(transcriptomes) + 55
    chart_width = width - left - right
    bar_step = chart_width / len(categories)
    bar_width = bar_step * 0.34
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="Mismatch and wobble profiles for original siRNAs versus matched shuffles">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="170" y="25" font-family="Arial" font-size="19" font-weight="bold">siRNA mismatch and wobble profiles: originals versus matched shuffles</text>',
        '<rect x="900" y="13" width="14" height="14" fill="#2b6cb0"/><text x="920" y="25" font-family="Arial" font-size="12">Original siRNAs</text>',
        '<rect x="1040" y="13" width="14" height="14" fill="#9aa5b1"/><text x="1060" y="25" font-family="Arial" font-size="12">Mean matched shuffle</text>',
    ]
    for row_index, transcriptome in enumerate(transcriptomes):
        frame_top = top + row_index * row_height
        baseline = frame_top + row_height - 36
        parts.extend([
            f'<rect x="{left}" y="{frame_top}" width="{chart_width}" height="{row_height - 36}" fill="none" stroke="#777" stroke-width="1"/>',
            f'<text x="{left - 12}" y="{frame_top + 20}" text-anchor="end" font-family="Arial" font-size="13" font-weight="bold">{escape(transcriptome)}</text>',
            f'<text x="{left - 12}" y="{frame_top + 38}" text-anchor="end" font-family="Arial" font-size="10">per siRNA / M searchable sites</text>',
            f'<text x="{left - 10}" y="{frame_top + 57}" text-anchor="end" font-family="Arial" font-size="10">{limit:.2f}</text>',
            f'<line x1="{left}" y1="{baseline}" x2="{width - right}" y2="{baseline}" stroke="#333" stroke-width="1"/>',
        ])
        for index, (mismatch, wobble) in enumerate(categories):
            x = left + index * bar_step + bar_step / 2
            for group, colour, offset in (("original", "#2b6cb0", -bar_width / 2), ("mean_matched_shuffle", "#9aa5b1", bar_width / 2)):
                row = lookup[(transcriptome, group, mismatch, wobble)]
                value = float(row["mean_alignment_density_per_million_sites_per_siRNA"])
                bar_height = value / limit * (row_height - 54)
                title = f"{transcriptome} | {group} | {mismatch} generic mismatches | {wobble} wobble pairs | {value:.4f} mean sites/million searchable sites/siRNA"
                parts.append(f'<rect x="{x + offset - bar_width / 2:.2f}" y="{baseline - bar_height:.2f}" width="{bar_width:.2f}" height="{bar_height:.2f}" fill="{colour}"><title>{escape(title)}</title></rect>')
            parts.append(f'<text x="{x:.2f}" y="{baseline + 14}" text-anchor="middle" font-family="Arial" font-size="10">{mismatch}m/{wobble}w</text>')
    parts.extend([
        f'<text x="{left + chart_width / 2:.2f}" y="{height - 8}" text-anchor="middle" font-family="Arial" font-size="12">Alignment category: generic mismatches / explicit wobble pairs. Shuffle bars are averaged per siRNA before comparison.</text>',
        '</svg>',
    ])
    path.write_text("\n".join(parts), encoding="utf-8")
    return summary


def mismatch_wobble_excess_heatmap(summary: list[dict[str, object]], path: Path) -> list[dict[str, object]]:
    """Show directly comparable original-minus-shuffle enrichment by alignment class."""
    transcriptomes = list(dict.fromkeys(str(row["transcriptome"]) for row in summary))
    categories = [(mismatch, wobble) for mismatch in range(4) for wobble in (0, 1)]
    lookup = {
        (str(row["transcriptome"]), str(row["guide_group"]), int(row["generic_mismatch_count"]), int(row["wobble_count"])): row
        for row in summary
    }
    rows: list[dict[str, object]] = []
    for transcriptome in transcriptomes:
        for mismatch, wobble in categories:
            original = float(lookup[(transcriptome, "original", mismatch, wobble)]["mean_alignment_density_per_million_sites_per_siRNA"])
            shuffle = float(lookup[(transcriptome, "mean_matched_shuffle", mismatch, wobble)]["mean_alignment_density_per_million_sites_per_siRNA"])
            rows.append({
                "transcriptome": transcriptome,
                "generic_mismatch_count": mismatch,
                "wobble_count": wobble,
                "original_mean_density_per_million_sites_per_siRNA": original,
                "matched_shuffle_mean_density_per_million_sites_per_siRNA": shuffle,
                "adaptation_excess_density_per_million_sites_per_siRNA": original - shuffle,
            })
    limit = max(0.01, max(abs(float(row["adaptation_excess_density_per_million_sites_per_siRNA"])) for row in rows))
    width, left, right, top, cell_height, cell_width = 1015, 180, 35, 32, 64, 100
    height = top + cell_height * len(transcriptomes) + 18
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="Mismatch and wobble adaptation heat map">',
        '<rect width="100%" height="100%" fill="white"/>',
    ]
    lookup_excess = {(str(row["transcriptome"]), int(row["generic_mismatch_count"]), int(row["wobble_count"])): row for row in rows}
    for column, (mismatch, wobble) in enumerate(categories):
        x = left + column * cell_width
        parts.append(f'<text x="{x + cell_width / 2:.2f}" y="{top - 9}" text-anchor="middle" font-family="Arial" font-size="12">{mismatch}m / {wobble}w</text>')
    for row_index, transcriptome in enumerate(transcriptomes):
        y = top + row_index * cell_height
        parts.append(f'<text x="{left - 14}" y="{y + 34}" text-anchor="end" font-family="Arial" font-size="13" font-weight="bold">{escape(transcriptome)}</text>')
        for column, (mismatch, wobble) in enumerate(categories):
            value = float(lookup_excess[(transcriptome, mismatch, wobble)]["adaptation_excess_density_per_million_sites_per_siRNA"])
            x = left + column * cell_width
            parts.append(f'<rect x="{x}" y="{y}" width="{cell_width - 2}" height="{cell_height - 2}" fill="{diverging_colour(value, limit)}" stroke="#fff" stroke-width="1"/>')
    parts.extend([
        '</svg>',
    ])
    path.write_text("\n".join(parts), encoding="utf-8")
    return rows


def batch_adaptation_summary_svg(ranking: list[dict[str, str]], path: Path) -> None:
    """Draw a compact cross-set transcriptome comparison from batch rankings."""
    input_sets = list(dict.fromkeys(row["input_set"] for row in ranking))
    transcriptomes = list(dict.fromkeys(row["transcriptome"] for row in ranking))
    lookup = {(row["input_set"], row["transcriptome"]): row for row in ranking}
    scores = [float(row["adaptation_excess_density_per_million_sites"]) for row in ranking]
    limit = max(0.25, max(abs(value) for value in scores) * 1.18)
    width, left, right, top = 1280, 105, 55, 62
    panel_height, gap = 330, 68
    chart_width = width - left - right
    height = top + panel_height * 2 + gap + 54
    colours = ["#2563eb", "#f97316", "#16a34a", "#9333ea", "#0891b2", "#dc2626"]
    bar_width = chart_width / max(1, len(transcriptomes)) / max(2, len(input_sets) + 1)

    def bar_y(value: float) -> float:
        return top + (limit - value) / (2 * limit) * (panel_height - 55)

    heatmap_top = top + panel_height + gap
    cell_width = chart_width / max(1, len(transcriptomes))
    cell_height = (panel_height - 35) / max(1, len(input_sets))
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="Batch siRNA set adaptation summary">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{left}" y="25" font-family="Arial" font-size="20" font-weight="bold">siRNA-set adaptation across transcriptomes</text>',
    ]
    for index, input_set in enumerate(input_sets):
        x = left + index * 180
        parts.append(f'<rect x="{x}" y="35" width="13" height="13" fill="{colours[index % len(colours)]}"/><text x="{x + 19}" y="46" font-family="Arial" font-size="12">{escape(input_set)}</text>')
    baseline = bar_y(0.0)
    parts.extend([
        f'<rect x="{left}" y="{top}" width="{chart_width}" height="{panel_height - 55}" fill="none" stroke="#777" stroke-width="1"/>',
        f'<line x1="{left}" y1="{baseline:.2f}" x2="{width - right}" y2="{baseline:.2f}" stroke="#333" stroke-width="1"/>',
        f'<text x="{left - 12}" y="{bar_y(limit) + 4:.2f}" text-anchor="end" font-family="Arial" font-size="11">+{limit:.2f}</text>',
        f'<text x="{left - 12}" y="{bar_y(-limit) + 4:.2f}" text-anchor="end" font-family="Arial" font-size="11">−{limit:.2f}</text>',
        f'<text x="27" y="{top + (panel_height - 55) / 2}" transform="rotate(-90 27 {top + (panel_height - 55) / 2})" text-anchor="middle" font-family="Arial" font-size="12">Adaptation excess per million sites</text>',
    ])
    for transcriptome_index, transcriptome in enumerate(transcriptomes):
        centre = left + chart_width * (transcriptome_index + 0.5) / len(transcriptomes)
        for set_index, input_set in enumerate(input_sets):
            row = lookup.get((input_set, transcriptome))
            if row is None:
                continue
            score = float(row["adaptation_excess_density_per_million_sites"])
            x = centre + (set_index - (len(input_sets) - 1) / 2) * (bar_width + 5)
            y = min(bar_y(score), baseline)
            height_value = abs(bar_y(score) - baseline)
            title = f"{input_set} | {transcriptome} | adaptation excess {score:.4f}; P={float(row['pvalue_vs_matched_shuffles']):.4g}; FDR={float(row['fdr_vs_matched_shuffles']):.4g}"
            parts.append(f'<rect x="{x - bar_width / 2:.2f}" y="{y:.2f}" width="{bar_width:.2f}" height="{height_value:.2f}" fill="{colours[set_index % len(colours)]}"><title>{escape(title)}</title></rect>')
        parts.append(f'<text x="{centre:.2f}" y="{top + panel_height - 34}" text-anchor="middle" font-family="Arial" font-size="12" font-weight="bold">{escape(transcriptome)}</text>')
    parts.extend([
        f'<text x="{left}" y="{heatmap_top - 12}" font-family="Arial" font-size="16" font-weight="bold">Adaptation heat map</text>',
        f'<rect x="{left}" y="{heatmap_top}" width="{chart_width}" height="{panel_height - 35}" fill="none" stroke="#777" stroke-width="1"/>',
    ])
    for column, transcriptome in enumerate(transcriptomes):
        x = left + column * cell_width
        parts.append(f'<text x="{x + cell_width / 2:.2f}" y="{heatmap_top - 22}" text-anchor="middle" font-family="Arial" font-size="12" font-weight="bold">{escape(transcriptome)}</text>')
        for row_index, input_set in enumerate(input_sets):
            row = lookup.get((input_set, transcriptome))
            if row is None:
                continue
            score = float(row["adaptation_excess_density_per_million_sites"])
            y = heatmap_top + row_index * cell_height
            pvalue = float(row["pvalue_vs_matched_shuffles"])
            title = f"{input_set} | {transcriptome} | excess {score:.4f}; P={pvalue:.4g}"
            parts.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{cell_width - 2:.2f}" height="{cell_height - 2:.2f}" fill="{diverging_colour(score, limit)}" stroke="#fff" stroke-width="1"><title>{escape(title)}</title></rect>')
            parts.append(f'<text x="{x + cell_width / 2:.2f}" y="{y + cell_height / 2 + 4:.2f}" text-anchor="middle" font-family="Arial" font-size="11">{score:+.2f}</text>')
    for row_index, input_set in enumerate(input_sets):
        y = heatmap_top + row_index * cell_height + cell_height / 2 + 4
        parts.append(f'<text x="{left - 12}" y="{y:.2f}" text-anchor="end" font-family="Arial" font-size="12" font-weight="bold">{escape(input_set)}</text>')
    parts.append('</svg>')
    path.write_text("\n".join(parts), encoding="utf-8")


def render_batch_adaptation_plots(analysis_dir: Path) -> dict[str, str]:
    tables = analysis_dir / "tables"
    plots = analysis_dir / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    ranking = read_rows(tables / "batch_multi_transcriptome_adaptation_ranking.tsv")
    batch_adaptation_summary_svg(ranking, plots / "batch_adaptation_summary.svg")
    return {"batch_summary": str(plots / "batch_adaptation_summary.svg")}


def run(analysis_dir: Path, input_mode: str = "direct_sirnas") -> dict[str, object]:
    if input_mode not in {"direct_sirnas", "dsrna_windows"}:
        raise ValueError("Input mode must be direct_sirnas or dsrna_windows.")
    tables = analysis_dir / "tables"
    rows = read_rows(tables / "multi_transcriptome_adaptation_per_sirna.tsv")
    ranking = read_rows(tables / "multi_transcriptome_adaptation_ranking.tsv")
    summaries: list[dict[str, object]] = []
    for row in rows:
        shuffled = [float(value) for value in row["shuffled_densities_per_million_sites"].split(",")]
        entropy = float(row.get("shannon_entropy") or normalized_shannon_entropy(row["siRNA_sequence"]))
        excess = float(row["original_density_per_million_sites"]) - sum(shuffled) / len(shuffled)
        signed = excess if row["strand"] == "sense" else -excess
        summaries.append({
            "transcriptome": row["transcriptome"], "dsrna_id": row["dsrna_id"], "dsrna_position": int(row["dsrna_position"]),
            "strand": row["strand"], "siRNA_id": row["siRNA_id"], "siRNA_sequence": row["siRNA_sequence"],
            "shannon_entropy": entropy, "original_alignment_count": int(row["original_alignment_count"]),
            "original_has_hit": int(row["original_alignment_count"]) > 0,
            "original_target_transcript_count": int(row["original_target_transcript_count"]),
            "original_density_per_million_sites": float(row["original_density_per_million_sites"]),
            "mean_matched_shuffle_density_per_million_sites": sum(shuffled) / len(shuffled),
            "adaptation_excess_density_per_million_sites": excess,
            "signed_adaptation_excess_density_per_million_sites": signed,
        })
    summary_rows: list[dict[str, object]] = []
    by_sirna: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in summaries:
        by_sirna[str(row["siRNA_id"])].append(row)
    for group in by_sirna.values():
        for row in group:
            other = [float(item["adaptation_excess_density_per_million_sites"]) for item in group if item is not row]
            relative = float(row["adaptation_excess_density_per_million_sites"]) - sum(other) / len(other)
            row["relative_adaptation_excess_density_per_million_sites"] = relative
            row["signed_relative_adaptation_excess_density_per_million_sites"] = relative if row["strand"] == "sense" else -relative
    for transcriptome in dict.fromkeys(str(row["transcriptome"]) for row in summaries):
        selected = [row for row in summaries if row["transcriptome"] == transcriptome]
        summary_rows.append({
            "transcriptome": transcriptome, "original_siRNA_count": len(selected),
            "original_siRNAs_with_at_least_one_hit": sum(bool(row["original_has_hit"]) for row in selected),
            "percent_original_siRNAs_with_at_least_one_hit": 100 * sum(bool(row["original_has_hit"]) for row in selected) / len(selected),
            "mean_original_siRNA_entropy": sum(float(row["shannon_entropy"]) for row in selected) / len(selected),
        })
    write_rows(tables / "multi_transcriptome_siRNA_summary.tsv", summaries)
    write_rows(tables / "multi_transcriptome_hit_summary.tsv", summary_rows)
    plot_dir = analysis_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    variability_rows = overall_summary_svg(summaries, ranking, plot_dir / "overall_siRNA_summary.svg")
    write_rows(tables / "multi_transcriptome_siRNA_variability.tsv", variability_rows)
    adaptability_complexity_svg(summaries, plot_dir / "siRNA_adaptability_complexity.svg")
    alignment_summary = alignment_quality_summary(
        read_rows(tables / "multi_transcriptome_adaptation_target_hits.tsv"), rows, ranking,
        plot_dir / "siRNA_mismatch_wobble_summary.svg",
    )
    write_rows(tables / "multi_transcriptome_mismatch_wobble_summary.tsv", alignment_summary)
    mismatch_wobble_excess = mismatch_wobble_excess_heatmap(
        alignment_summary, plot_dir / "siRNA_mismatch_wobble_adaptation_heatmap.svg"
    )
    write_rows(tables / "multi_transcriptome_mismatch_wobble_excess.tsv", mismatch_wobble_excess)
    if input_mode == "dsrna_windows":
        for dsrna_id in dict.fromkeys(str(row["dsrna_id"]) for row in summaries):
            landscape_svg([row for row in summaries if row["dsrna_id"] == dsrna_id], plot_dir / f"{safe_name(dsrna_id)}.multi_transcriptome_adaptation.svg", dsrna_id)
    return {"summary": str(tables / "multi_transcriptome_hit_summary.tsv"), "per_sirna": str(tables / "multi_transcriptome_siRNA_summary.tsv"), "mismatch_wobble_summary": str(tables / "multi_transcriptome_mismatch_wobble_summary.tsv"), "mismatch_wobble_excess": str(tables / "multi_transcriptome_mismatch_wobble_excess.tsv"), "plots": str(plot_dir), "adaptability_complexity_plot": str(plot_dir / "siRNA_adaptability_complexity.svg"), "mismatch_wobble_plot": str(plot_dir / "siRNA_mismatch_wobble_summary.svg"), "mismatch_wobble_excess_heatmap": str(plot_dir / "siRNA_mismatch_wobble_adaptation_heatmap.svg")}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", required=True, type=Path)
    parser.add_argument("--input-mode", choices=("direct_sirnas", "dsrna_windows"), default="direct_sirnas")
    args = parser.parse_args()
    print(run(args.analysis_dir, args.input_mode))


if __name__ == "__main__":
    main()
