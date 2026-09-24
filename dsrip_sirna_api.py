"""API adapter for dsRIP siRNA/dsRNA efficiency prediction.

This module keeps the INCI pipeline decoupled from the dsRIP Flask app.  The
original dsRIP implementation expects to run from ``dsRIP_web/main_site`` and
mostly communicates through files; this adapter turns that behavior into a
callable function that returns structured results and INCI-friendly outputs.
"""

from __future__ import annotations

import html
import importlib
import os
import sys
import json
import argparse
import zipfile
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any


APP_DIR = Path(__file__).resolve().parent
DSRIP_MAIN_SITE = APP_DIR / "dsRIP" / "dsRIP_web" / "main_site"


class DsRipApiError(RuntimeError):
    """Raised when dsRIP prediction cannot be run."""


@dataclass(frozen=True)
class SiRnaPredictionConfig:
    """Configuration for dsRIP siRNA/dsRNA efficiency prediction."""

    fasta_path: Path
    output_dir: Path
    off_target_species: str = ""
    min_length: int = 280
    max_length: int = 310
    sirna_length: int = 21
    buffer_size: int = 10
    max_genes: int = 5
    only_orf: bool = True
    orf_correction: bool = False
    safety: bool = False
    safety_essential_coeff: int = 20
    efficacy_coeff: int = 50
    safety_coeff: int = 50
    user_mismatch: int = 0
    new_off_target_on: bool = False
    new_off_target_species_input_name: str = "user_transcriptome"
    sirna_parameter_overrides: dict[str, float] | None = None


def _install_numpy_compatibility_aliases() -> None:
    """Provide deprecated NumPy aliases still used by older dsRIP dependencies."""

    try:
        import numpy as np
    except ModuleNotFoundError:
        return

    if not hasattr(np, "float"):
        np.float = float  # type: ignore[attr-defined]


@contextmanager
def _dsrip_runtime() -> Any:
    """Temporarily provide dsRIP's expected import path and working directory."""

    if not DSRIP_MAIN_SITE.exists():
        raise DsRipApiError(f"dsRIP main_site folder was not found: {DSRIP_MAIN_SITE}")

    previous_cwd = Path.cwd()
    added_path = False
    main_site_str = str(DSRIP_MAIN_SITE)
    if main_site_str not in sys.path:
        sys.path.insert(0, main_site_str)
        added_path = True

    try:
        os.chdir(DSRIP_MAIN_SITE)
        yield
    finally:
        os.chdir(previous_cwd)
        if added_path:
            try:
                sys.path.remove(main_site_str)
            except ValueError:
                pass


def _load_dsrip_efficiency() -> Any:
    try:
        _install_numpy_compatibility_aliases()
        with _dsrip_runtime():
            return importlib.import_module("dsRNA_efficiency")
    except ModuleNotFoundError as exc:
        missing = exc.name or "unknown dependency"
        raise DsRipApiError(
            "dsRIP siRNA prediction is missing a Python dependency: "
            f"{missing}. Run INCI with an environment that has the dsRIP "
            "requirements installed."
        ) from exc
    except Exception as exc:
        raise DsRipApiError(f"Could not load dsRIP siRNA prediction: {exc}") from exc


def _clean_number(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.3f}".rstrip("0").rstrip(".")
    return str(value)


def _write_tsv(path: Path, records: list[dict[str, Any]]) -> None:
    columns = [
        "gene_name",
        "input_length",
        "orf",
        "start_siRNA",
        "end_siRNA",
        "window_size",
        "selected_dsRNA_total_score",
        "average_gene_total_score",
        "selected_dsRNA_efficiency_score",
        "selected_dsRNA_safety_score",
        "optimized_dsRNA_length",
        "optimized_dsRNA",
    ]
    lines = ["\t".join(columns)]
    for record in records:
        lines.append("\t".join(_clean_number(record.get(column, "")) for column in columns))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_fasta(path: Path, records: list[dict[str, Any]]) -> None:
    chunks = []
    for record in records:
        sequence = record.get("optimized_dsRNA") or ""
        if sequence:
            chunks.append(f">{record['gene_name']}\n{sequence}")
    path.write_text("\n".join(chunks) + ("\n" if chunks else ""), encoding="utf-8")


def _write_off_target_summary(path: Path, safety: bool, records: list[dict[str, Any]]) -> None:
    lines = ["gene_name\tsafety_enabled\tselected_dsRNA_safety_score\tnote"]
    note = "Off-target scoring was enabled." if safety else "Off-target scoring was not enabled for this INCI run."
    for record in records:
        lines.append(
            "\t".join(
                [
                    str(record.get("gene_name", "")),
                    str(bool(safety)),
                    _clean_number(record.get("selected_dsRNA_safety_score", "")),
                    note,
                ]
            )
        )
    if not records:
        lines.append(f"\t{bool(safety)}\t\t{note}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_report(path: Path, records: list[dict[str, Any]], errors: list[str], outputs: dict[str, Path]) -> None:
    rows = []
    for record in records:
        rows.append(
            "<tr>"
            f"<td>{html.escape(str(record.get('gene_name', '')))}</td>"
            f"<td>{html.escape(_clean_number(record.get('optimized_dsRNA_length', '')))}</td>"
            f"<td>{html.escape(_clean_number(record.get('selected_dsRNA_total_score', '')))}</td>"
            f"<td>{html.escape(str(record.get('start_siRNA', '')))}</td>"
            f"<td>{html.escape(str(record.get('end_siRNA', '')))}</td>"
            "</tr>"
        )

    error_items = "".join(f"<li>{html.escape(error)}</li>" for error in errors)
    output_items = "".join(
        f"<li><strong>{html.escape(label)}</strong>: {html.escape(str(path_value))}</li>"
        for label, path_value in outputs.items()
    )
    path.write_text(
        f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>dsRNA Enhancement Report</title>
  <style>
    body {{ font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 32px; color: #17212b; }}
    table {{ border-collapse: collapse; width: 100%; margin-top: 16px; }}
    th, td {{ border-bottom: 1px solid #d9e3ea; padding: 9px; text-align: left; vertical-align: top; }}
    th {{ color: #115e59; }}
    code {{ overflow-wrap: anywhere; }}
  </style>
</head>
<body>
  <h1>dsRNA Enhancement Report</h1>
  <p>Generated from the dsRIP siRNA efficiency prediction module.</p>
  <h2>Predicted Regions</h2>
  <table>
    <thead><tr><th>Gene</th><th>Length</th><th>Total score</th><th>Start siRNA</th><th>End siRNA</th></tr></thead>
    <tbody>{''.join(rows) if rows else '<tr><td colspan="5">No regions were predicted.</td></tr>'}</tbody>
  </table>
  <h2>Outputs</h2>
  <ul>{output_items}</ul>
  <h2>Warnings</h2>
  <ul>{error_items if error_items else '<li>No warnings were reported.</li>'}</ul>
</body>
</html>
""",
        encoding="utf-8",
    )


def _zip_outputs(zip_path: Path, output_dir: Path, outputs: dict[str, Path]) -> None:
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in outputs.values():
            if path.exists() and path.is_file():
                archive.write(path, path.relative_to(output_dir))

        efficiency_dir = output_dir / "efficiency"
        if efficiency_dir.exists():
            for path in efficiency_dir.rglob("*"):
                if path.is_file():
                    archive.write(path, path.relative_to(output_dir))


def predict_sirna_regions(config: SiRnaPredictionConfig) -> dict[str, Any]:
    """Run dsRIP siRNA efficiency prediction and return structured outputs."""

    fasta_path = config.fasta_path.expanduser().resolve()
    if not fasta_path.exists():
        raise DsRipApiError(f"Input FASTA was not found: {fasta_path}")

    output_dir = config.output_dir.expanduser().resolve()
    efficiency_dir = output_dir / "efficiency"
    off_target_dir = efficiency_dir / "off_target"
    efficiency_dir.mkdir(parents=True, exist_ok=True)
    off_target_dir.mkdir(parents=True, exist_ok=True)

    dsrip = _load_dsrip_efficiency()
    errors: list[str] = []
    records: list[dict[str, Any]] = []
    fasta_dict: dict[str, dict[str, Any]] = {}

    with _dsrip_runtime():
        if config.new_off_target_on:
            try:
                dsrip.run_cd_hit_and_indexer(
                    str(efficiency_dir), c_param=0.95, min_length=200,
                    user_species_name=config.new_off_target_species_input_name
                )
            except Exception as exc:
                errors.append(f"Off-target index preparation failed: {exc}")

        sequence_data = dsrip.read_fasta_file(str(fasta_path))
        if not sequence_data:
            raise DsRipApiError(f"No FASTA records could be read from {fasta_path}")

        min_window = int(config.min_length) + 1 - (int(config.sirna_length) - 1)
        max_window = int(config.max_length) + 1 - (int(config.sirna_length) - 1)

        for raw_name, raw_sequence in sequence_data[: config.max_genes]:
            name = dsrip.sanitize_sheet_title(str(raw_name))
            try:
                sequence = str(raw_sequence).upper().replace("N", "")
                if not dsrip.check_valid_dna(sequence):
                    errors.append(f"{name}: skipped because the sequence contains non-DNA/RNA bases.")
                    continue

                if config.orf_correction:
                    corrected_sequence = str(dsrip.correct_sequence_based_on_ORF(sequence, name))
                else:
                    corrected_sequence = sequence

                gene_orf = str(dsrip.ORF_finder_module.get_ORF(corrected_sequence, name)[4])
                rna_sequence = corrected_sequence.upper().replace("T", "U")
                orf_info = dsrip.ORF_siRNA(rna_sequence)
                sirnas, sirnas_dna_off_target = dsrip.generate_siRNA(rna_sequence, siRNA_length=config.sirna_length)

                off_target_report = []
                if config.safety:
                    try:
                        _, _, _, off_target_report, _ = dsrip.off_target_siRNA_all(
                            sirnas_dna_off_target,
                            len(corrected_sequence),
                            name,
                            config.off_target_species,
                            str(output_dir),
                            siRNA_length=config.sirna_length,
                            user_mismatch=config.user_mismatch,
                            new_off_target_species_input_name=config.new_off_target_species_input_name,
                        )
                    except Exception as exc:
                        errors.append(f"{name}: off-target scoring failed, continuing without it ({exc}).")

                access = dsrip.mRNA_accessibility(rna_sequence, siRNA_length=config.sirna_length)
                features = dsrip.siRNA_feature_prediction(
                    sirnas,
                    rna_sequence,
                    name,
                    config.sirna_length,
                    access,
                    orf_info,
                    off_target_report,
                )
                gene_dict = {
                    "gene_name": name,
                    "input_gene_seq": sequence,
                    "corrected_gene_seq": corrected_sequence,
                    "gene_ORF": gene_orf,
                    "gene_RNA_seq": rna_sequence,
                    "siRNA_features": features,
                }
                scored_features = dsrip.score_and_update_siRNA_features(
                    gene_dict,
                    safety=config.safety,
                    safety_essential_coeff=config.safety_essential_coeff,
                    parameter_overrides=config.sirna_parameter_overrides,
                )
                best_region = dsrip.find_best_siRNA_region(
                    scored_features,
                    min_window=min_window,
                    max_window=max_window,
                    only_ORF=config.only_orf,
                    safety=config.safety,
                    safety_priority=config.safety_coeff,
                    efficiency_priority=config.efficacy_coeff,
                    user_dir=str(output_dir),
                )
                optimized = dsrip.extract_subsequence(
                    corrected_sequence,
                    best_region,
                    bufsi=config.buffer_size,
                )

                record = {
                    "gene_name": name,
                    "input_length": len(sequence),
                    "orf": gene_orf,
                    "optimized_dsRNA": optimized or "",
                    "optimized_dsRNA_length": len(optimized or ""),
                    **best_region,
                }
                records.append(record)
                fasta_dict[name] = {**gene_dict, "siRNA_features": scored_features["siRNA_features"]}
            except Exception as exc:
                errors.append(f"{name}: prediction failed ({exc}).")

        if fasta_dict:
            try:
                dsrip.write_fasta_dict_to_excel(fasta_dict, str(output_dir))
            except Exception as exc:
                errors.append(f"Could not write siRNA prediction workbook: {exc}")

    outputs = {
        "enhanced_designs": output_dir / "enhanced_dsrna_designs.tsv",
        "best_regions_fasta": efficiency_dir / "best_dsRNA_regions.fasta",
        "off_target_summary": output_dir / "off_target_summary.tsv",
        "report": output_dir / "dsrna_enhancement_report.html",
        "prediction_workbook": efficiency_dir / "siRNA_predictions.xlsx",
        "result_manifest": output_dir / "dsrip_efficiency_outputs.json",
    }

    _write_tsv(outputs["enhanced_designs"], records)
    _write_fasta(outputs["best_regions_fasta"], records)
    _write_off_target_summary(outputs["off_target_summary"], config.safety, records)
    _write_report(outputs["report"], records, errors, outputs)
    outputs["result_manifest"].write_text(
        json.dumps(
            {
                "records": records,
                "errors": errors,
                "outputs": {key: str(value) for key, value in outputs.items()},
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    return {
        "records": records,
        "errors": errors,
        "outputs": {key: str(value) for key, value in outputs.items()},
    }


def run_efficiency_from_fasta(
    fasta_path: str | Path,
    output_path: str | Path,
    *,
    min_length: int = 280,
    max_length: int = 310,
    sirna_length: int = 21,
    buffer_size: int = 10,
    max_genes: int = 5,
    only_orf: bool = True,
    orf_correction: bool = False,
    safety: bool = False,
    off_target_species: str = "",
    output_zip: str | Path | None = None,
    sirna_parameter_overrides: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Run the dsRIP web-style dsRNA efficiency workflow on a FASTA file.

    ``output_path`` is the output folder. The main dsRIP-style outputs are:
    ``efficiency/best_dsRNA_regions.fasta`` and
    ``efficiency/siRNA_predictions.xlsx``.
    """

    result = predict_sirna_regions(
        SiRnaPredictionConfig(
            fasta_path=Path(fasta_path),
            output_dir=Path(output_path),
            off_target_species=off_target_species,
            min_length=min_length,
            max_length=max_length,
            sirna_length=sirna_length,
            buffer_size=buffer_size,
            max_genes=max_genes,
            only_orf=only_orf,
            orf_correction=orf_correction,
            safety=safety,
            sirna_parameter_overrides=sirna_parameter_overrides,
        )
    )

    if output_zip:
        zip_path = Path(output_zip).expanduser().resolve()
        output_dir = Path(output_path).expanduser().resolve()
        _zip_outputs(zip_path, output_dir, {key: Path(value) for key, value in result["outputs"].items()})
        result["outputs"]["zip"] = str(zip_path)

    return result


def _parse_parameter_override(values: list[str]) -> dict[str, float]:
    overrides: dict[str, float] = {}
    for value in values:
        if "=" not in value:
            raise argparse.ArgumentTypeError(f"Parameter override must look like key=value: {value}")
        key, raw_number = value.split("=", 1)
        overrides[key.strip()] = float(raw_number.strip())
    return overrides


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run dsRIP web-style dsRNA efficiency prediction on a FASTA file."
    )
    parser.add_argument("fasta", help="Input FASTA file.")
    parser.add_argument("output", help="Output folder for dsRIP efficiency files.")
    parser.add_argument("--zip", dest="output_zip", help="Optional ZIP file containing the generated outputs.")
    parser.add_argument("--min-length", type=int, default=280, help="Minimum dsRNA length. Default: 280.")
    parser.add_argument("--max-length", type=int, default=310, help="Maximum dsRNA length. Default: 310.")
    parser.add_argument("--sirna-length", type=int, default=21, help="siRNA length. Default: 21.")
    parser.add_argument("--buffer-size", type=int, default=10, help="Extra nt around selected siRNA window. Default: 10.")
    parser.add_argument("--max-genes", type=int, default=5, help="Maximum FASTA records to process. Default: 5.")
    parser.add_argument("--include-utr", action="store_true", help="Allow selected regions outside the ORF.")
    parser.add_argument("--orf-correction", action="store_true", help="Apply dsRIP ORF correction before prediction.")
    parser.add_argument("--safety", action="store_true", help="Enable off-target safety scoring.")
    parser.add_argument("--off-target-species", default="", help="Comma-separated dsRIP off-target species names.")
    parser.add_argument(
        "--set-param",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="Override one built-in siRNA scoring parameter. Can be used more than once.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    try:
        result = run_efficiency_from_fasta(
            args.fasta,
            args.output,
            min_length=args.min_length,
            max_length=args.max_length,
            sirna_length=args.sirna_length,
            buffer_size=args.buffer_size,
            max_genes=args.max_genes,
            only_orf=not args.include_utr,
            orf_correction=args.orf_correction,
            safety=args.safety,
            off_target_species=args.off_target_species,
            output_zip=args.output_zip,
            sirna_parameter_overrides=_parse_parameter_override(args.set_param),
        )
    except (DsRipApiError, ValueError, argparse.ArgumentTypeError) as exc:
        parser.exit(1, f"Error: {exc}\n")

    print(json.dumps(result["outputs"], indent=2))
    if result["errors"]:
        print("\nWarnings:")
        for error in result["errors"]:
            print(f"- {error}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
