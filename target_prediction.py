#!/usr/bin/env python3
"""CleaveLand/GSTAr-style sRNA target prediction without degradome evidence."""

from __future__ import annotations

import html
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Sequence

from degradome_analysis import (
    DegradomeConfig,
    FastaRecord,
    find_targets,
    read_fasta,
    read_fasta_text,
    safe_name,
    write_fasta,
    write_rows,
)
from plot_style import apply_matplotlib_style


LogFn = Callable[[str, str], None]
os.environ.setdefault("MPLCONFIGDIR", str(Path(os.environ.get("TMPDIR", "/tmp")) / "matplotlib-inci"))


@dataclass(frozen=True)
class TargetPredictionConfig:
    srna_text: str
    srna_fasta: Path | None
    transcript_text: str
    transcript_fasta: Path | None
    output_dir: Path
    ignore_query_pos1: bool = True
    mfe_ratio_cutoff: float = 0.70
    max_allen_score: float | None = None
    max_mismatches: int | None = None
    sort_by: str = "mfe_ratio"
    plot_metric: str = "mfe_ratio"
    max_transcript_plots: int = 100
    threads: int = 4


def load_prediction_inputs(config: TargetPredictionConfig) -> tuple[list[FastaRecord], dict[str, str]]:
    srnas = read_fasta(config.srna_fasta) if config.srna_fasta else read_fasta_text(config.srna_text, "pasted_sRNA")
    transcript_records = (
        read_fasta(config.transcript_fasta)
        if config.transcript_fasta
        else read_fasta_text(config.transcript_text, "pasted_transcript")
    )
    if not srnas:
        raise ValueError("Provide sRNAs by pasting sequences or selecting an sRNA FASTA.")
    if not transcript_records:
        raise ValueError("Provide transcripts by pasting sequences or selecting a transcript FASTA.")
    return srnas, {record.name: record.sequence for record in transcript_records}


def retained_pair_rows(hits: Sequence[Any], config: TargetPredictionConfig) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for hit in hits:
        if float(hit.mfe_ratio) < config.mfe_ratio_cutoff:
            continue
        if config.max_allen_score is not None and float(hit.allen_score) > config.max_allen_score:
            continue
        if config.max_mismatches is not None and int(hit.mismatch_count) > config.max_mismatches:
            continue
        rows.append(
            {
                "rank": hit.rank,
                "srna": hit.query,
                "transcript": hit.transcript,
                "t_start": hit.t_start,
                "t_stop": hit.t_stop,
                "canonical_q10_slice_site": hit.slice_sites.get(10, ""),
                "mfe_ratio": hit.mfe_ratio,
                "mfe_perfect": hit.mfe_perfect,
                "mfe_site": hit.mfe_site,
                "allen_score": hit.allen_score,
                "mismatches": hit.mismatch_count,
                "gu_wobbles": hit.gu_wobble_count,
                "bulges": hit.bulge_count,
                "paired_regions": hit.paired,
                "unpaired_regions": hit.unpaired,
                "dot_bracket_structure": hit.structure,
                "aligned_sequence": hit.sequence,
                "match_pattern": hit.match_pattern,
                "pair_count": hit.pair_count,
                "srna_cpm_from_name": hit.srna_cpm,
                "ignored_query_pos1": hit.ignored_query_pos1,
                "original_query_length": hit.original_query_length,
                "effective_query_length": hit.effective_query_length,
            }
        )
    rows.sort(key=lambda row: (-float(row["mfe_ratio"]), float(row["allen_score"]), int(row["mismatches"])))
    return rows


def transcript_summary_rows(transcripts: dict[str, str], rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    by_tx: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        by_tx.setdefault(str(row["transcript"]), []).append(row)
    summary: list[dict[str, Any]] = []
    for transcript, sequence in transcripts.items():
        tx_rows = by_tx.get(transcript, [])
        summary.append(
            {
                "transcript": transcript,
                "length": len(sequence),
                "targeting_srna_count": len({str(row["srna"]) for row in tx_rows}),
                "retained_pair_count": len(tx_rows),
                "best_mfe_ratio": max((float(row["mfe_ratio"]) for row in tx_rows), default=0.0),
                "best_allen_score": min((float(row["allen_score"]) for row in tx_rows), default=""),
                "best_mismatch_count": min((int(row["mismatches"]) for row in tx_rows), default=""),
            }
        )
    summary.sort(key=lambda row: (int(row["retained_pair_count"]), int(row["targeting_srna_count"])), reverse=True)
    return summary


def srna_summary_rows(srnas: Sequence[FastaRecord], rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    by_srna: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        by_srna.setdefault(str(row["srna"]), []).append(row)
    summary: list[dict[str, Any]] = []
    for srna in srnas:
        srna_rows = by_srna.get(srna.name, [])
        summary.append(
            {
                "srna": srna.name,
                "length": len(srna.sequence),
                "targeted_transcript_count": len({str(row["transcript"]) for row in srna_rows}),
                "retained_pair_count": len(srna_rows),
                "best_mfe_ratio": max((float(row["mfe_ratio"]) for row in srna_rows), default=0.0),
                "best_allen_score": min((float(row["allen_score"]) for row in srna_rows), default=""),
                "best_mismatch_count": min((int(row["mismatches"]) for row in srna_rows), default=""),
            }
        )
    summary.sort(key=lambda row: (int(row["targeted_transcript_count"]), int(row["retained_pair_count"])), reverse=True)
    return summary


def plot_target_prediction_summary(
    transcript_summary: Sequence[dict[str, Any]],
    srna_summary: Sequence[dict[str, Any]],
    path: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    style = apply_matplotlib_style()
    import matplotlib.pyplot as plt

    targeted = sum(1 for row in transcript_summary if int(row["targeting_srna_count"]) > 0)
    untargeted = max(0, len(transcript_summary) - targeted)
    top_tx = [row for row in transcript_summary if int(row["retained_pair_count"]) > 0][:10]
    top_srna = [row for row in srna_summary if int(row["targeted_transcript_count"]) > 0][:10]

    def short_label(value: Any, limit: int = 32) -> str:
        text = str(value)
        return text if len(text) <= limit else text[: limit - 3] + "..."

    row_count = max(len(top_tx), len(top_srna), 4)
    dynamic_height = min(14.0, max(float(style["figure_height"]), row_count * 0.42 + 2.4))
    fig, axes = plt.subplots(1, 3, figsize=(style["figure_width"], dynamic_height), constrained_layout=True)
    axes[0].bar(["Targeted", "No hit"], [targeted, untargeted], color=["#0f766e", "#94a3b8"])
    axes[0].set_ylabel("Transcript count")
    axes[0].set_title("Transcripts with >=1 predicted sRNA")

    tx_labels = [short_label(row["transcript"]) for row in reversed(top_tx)]
    tx_values = [int(row["retained_pair_count"]) for row in reversed(top_tx)]
    axes[1].barh(tx_labels, tx_values, color="#2563eb")
    axes[1].set_xlabel("Retained sRNA-transcript pairs")
    axes[1].set_title("Top transcripts")

    srna_labels = [short_label(row["srna"]) for row in reversed(top_srna)]
    srna_values = [int(row["targeted_transcript_count"]) for row in reversed(top_srna)]
    axes[2].barh(srna_labels, srna_values, color="#ea580c")
    axes[2].set_xlabel("Targeted transcripts")
    axes[2].set_title("Top sRNAs")

    for ax in axes:
        ax.grid(axis="x" if ax is not axes[0] else "y", color="#e5e7eb", linewidth=style["grid_width"])
        ax.tick_params(axis="y", labelsize=max(6.0, float(style["font_size"]) - 2.0))
    fig.savefig(path, dpi=style["dpi"])
    plt.close(fig)


def plot_transcript_target_predictions(
    transcript: str,
    length: int,
    rows: Sequence[dict[str, Any]],
    plot_metric: str,
    path: Path,
) -> None:
    """Plot retained target predictions at their canonical q10 positions."""

    import matplotlib

    matplotlib.use("Agg")
    style = apply_matplotlib_style()
    import matplotlib.pyplot as plt
    import numpy as np

    points: list[tuple[int, float]] = []
    for row in rows:
        q10 = row.get("canonical_q10_slice_site")
        if q10 in (None, ""):
            continue
        x_value = int(q10)
        if not 1 <= x_value <= length:
            continue
        y_value = float(row["mfe_ratio"]) if plot_metric == "mfe_ratio" else -float(row["allen_score"])
        points.append((x_value, y_value))
    if not points:
        return

    points.sort(key=lambda point: (point[0], point[1]))
    x_values = np.array([point[0] for point in points], dtype=float)
    y_values = np.array([point[1] for point in points], dtype=float)
    color = "#2563eb" if plot_metric == "mfe_ratio" else "#ea580c"
    ylabel = "MFE ratio" if plot_metric == "mfe_ratio" else "Reversed Allen score (-penalty)"

    fig, ax = plt.subplots(figsize=(style["figure_width"], style["figure_height"]), constrained_layout=True)
    ax.vlines(
        x_values,
        0.0,
        y_values,
        color=color,
        linewidth=max(0.3, float(style["line_width"]) * 0.42),
        alpha=0.72,
        zorder=2,
    )
    ax.scatter(
        x_values,
        y_values,
        s=float(style["marker_size"]) ** 2 * 0.42,
        color=color,
        edgecolor="#202020",
        linewidth=max(0.18, float(style["line_width"]) * 0.1),
        alpha=0.84,
        zorder=3,
    )
    ax.axhline(0.0, color="#64748b", linewidth=max(0.24, float(style["line_width"]) * 0.2), zorder=1)
    ax.set_xlim(0.5, length + 0.5)
    if plot_metric == "mfe_ratio":
        ax.set_ylim(0.0, max(1.0, float(np.max(y_values)) * 1.08))
    else:
        lower = min(float(np.min(y_values)) * 1.12, -0.25)
        ax.set_ylim(lower, 0.12)
    ax.set_xlabel("Transcript position (nt); marker/stem at canonical q10 site")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{transcript}: {len(points):,} retained sRNA target pair(s)")
    ax.grid(axis="y", color="#e5e7eb", linewidth=style["grid_width"], alpha=0.85)
    fig.savefig(path, dpi=max(int(style["dpi"]), 300))
    plt.close(fig)


def write_report(path: Path, result: dict[str, Any]) -> None:
    links = "".join(
        f"<li><strong>{html.escape(label)}</strong>: {html.escape(str(file_path))}</li>"
        for label, file_path in result["outputs"].items()
    )
    body = f"""<!doctype html>
<html>
<head><meta charset="utf-8"><title>Target Prediction Report</title></head>
<body>
<h1>Target Prediction Report</h1>
<p>CleaveLand/GSTAr-style sRNA-transcript target prediction without degradome evidence.</p>
<ul>
<li>sRNAs: {result['srna_count']}</li>
<li>Transcripts: {result['transcript_count']}</li>
<li>RNAplex target sites passing the MFE-ratio cutoff: {result['raw_target_count']}</li>
<li>Retained sRNA-transcript pairs after optional Allen/mismatch filtering: {result['retained_pair_count']}</li>
<li>Transcript prediction plots: {result['transcript_plot_count']}</li>
</ul>
<h2>Outputs</h2>
<ul>{links}</ul>
</body>
</html>
"""
    path.write_text(body, encoding="utf-8")


def run_target_prediction(config: TargetPredictionConfig, log: LogFn | None = None) -> dict[str, Any]:
    outdir = config.output_dir
    outdir.mkdir(parents=True, exist_ok=True)
    for folder in ("inputs", "plots", "tables"):
        (outdir / folder).mkdir(exist_ok=True)

    srnas, transcripts = load_prediction_inputs(config)
    write_fasta(srnas, outdir / "inputs" / "srnas.fasta")
    write_fasta([FastaRecord(name, sequence) for name, sequence in transcripts.items()], outdir / "inputs" / "transcripts.fasta")
    if log:
        log(f"Loaded {len(srnas):,} sRNA(s) and {len(transcripts):,} transcript(s).", "info")

    engine_config = DegradomeConfig(
        srna_text="",
        srna_fasta=None,
        transcript_text="",
        transcript_fasta=None,
        samples=(),
        output_dir=outdir,
        ignore_query_pos1=config.ignore_query_pos1,
        slice_positions=(10,),
        mfe_ratio_cutoff=config.mfe_ratio_cutoff,
        sort_by=config.sort_by,
        threads=config.threads,
    )
    hits = find_targets(srnas, transcripts, engine_config, log)
    rows = retained_pair_rows(hits, config)
    tx_summary = transcript_summary_rows(transcripts, rows)
    srna_summary = srna_summary_rows(srnas, rows)

    pairs_path = outdir / "target_prediction_pairs.tsv"
    tx_summary_path = outdir / "tables" / "transcript_target_summary.tsv"
    srna_summary_path = outdir / "tables" / "srna_target_summary.tsv"
    plot_path = outdir / "plots" / "target_prediction_summary.png"
    report_path = outdir / "target_prediction_report.html"
    write_rows(pairs_path, rows)
    write_rows(tx_summary_path, tx_summary)
    write_rows(srna_summary_path, srna_summary)
    plot_target_prediction_summary(tx_summary, srna_summary, plot_path)

    rows_by_transcript: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        rows_by_transcript.setdefault(str(row["transcript"]), []).append(row)
    targeted_transcripts = [
        str(row["transcript"])
        for row in tx_summary
        if int(row["retained_pair_count"]) > 0
    ]
    plot_limit = int(config.max_transcript_plots)
    if plot_limit > 0:
        targeted_transcripts = targeted_transcripts[:plot_limit]
    transcript_plot_paths: list[Path] = []
    for index, transcript in enumerate(targeted_transcripts, start=1):
        transcript_plot_path = outdir / "plots" / f"{safe_name(transcript)}.target_predictions.png"
        plot_transcript_target_predictions(
            transcript,
            len(transcripts[transcript]),
            rows_by_transcript[transcript],
            config.plot_metric,
            transcript_plot_path,
        )
        if transcript_plot_path.exists():
            transcript_plot_paths.append(transcript_plot_path)
        if log and (index == len(targeted_transcripts) or index % 25 == 0):
            log(f"Wrote prediction plot {index:,}/{len(targeted_transcripts):,}: {transcript}", "info")

    result = {
        "outdir": str(outdir),
        "srna_count": len(srnas),
        "transcript_count": len(transcripts),
        "raw_target_count": len(hits),
        "retained_pair_count": len(rows),
        "transcript_plot_count": len(transcript_plot_paths),
        "outputs": {
            "retained_pairs": str(pairs_path),
            "transcript_summary": str(tx_summary_path),
            "srna_summary": str(srna_summary_path),
            "summary_plot": str(plot_path),
            "transcript_plots": str(outdir / "plots"),
            "report": str(report_path),
        },
    }
    write_report(report_path, result)
    if log:
        log(f"Wrote target prediction outputs to {outdir}.", "info")
    return result
