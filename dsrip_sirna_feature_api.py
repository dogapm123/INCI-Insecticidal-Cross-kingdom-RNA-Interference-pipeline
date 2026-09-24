"""Single-siRNA feature API for dsRIP calculations.

Use this module when you already have one antisense siRNA sequence and want the
same per-siRNA feature dictionary produced by dsRIP's ``siRNA_feature_prediction``.
"""

from __future__ import annotations

import argparse
import sqlite3
import json
from pathlib import Path
from typing import Any

from dsrip_sirna_api import DSRIP_MAIN_SITE, DsRipApiError


LOOKUP_DB = DSRIP_MAIN_SITE / "constant_files" / "lookup.db"
SIRNA_PARAMETERS_FILE = DSRIP_MAIN_SITE / "constant_files" / "siRNA_parameters.txt"

DEFAULT_SIRNA_PARAMETERS = {
    "asymScore_weight": 5.0,
    "self_fold_energy_weight": 3.0,
    "mRNA_accessibility_weight": 0.0,
    "edge_asyym_bonus": 10.0,
    "anti_GC_lower_bound": 0.0,
    "anti_GC_upper_bound": 65.0,
    "anti_GC_bonus": 0.0,
    "ORF_info_bonus": 0.0,
    "anti_10th_A_bonus": 10.0,
    "anti_9_14_GC_bonus": 5.0,
    "anti_9_14_GC_lower_bound": 35.0,
    "min_window": 250.0,
    "max_window": 350.0,
    "only_ORF": 0.0,
    "GGGG_CCCC_penalty": 0.0,
}


def normalize_antisense_sequence(sequence: str) -> str:
    """Return an uppercase RNA sequence after validating siRNA bases."""

    normalized = "".join(str(sequence).split()).upper().replace("T", "U")
    if not normalized:
        raise ValueError("Antisense siRNA sequence is empty.")

    invalid = sorted(set(normalized) - set("AUCG"))
    if invalid:
        raise ValueError(f"Antisense siRNA contains invalid base(s): {', '.join(invalid)}")

    if len(normalized) < 19:
        raise ValueError("Antisense siRNA must be at least 19 nt for dsRIP asymmetry calculations.")

    return normalized


def complement_rna(sequence: str) -> str:
    complement = str.maketrans({"A": "U", "U": "A", "C": "G", "G": "C"})
    return sequence.translate(complement)


def count_nucleotides(sequence: str) -> dict[str, int]:
    return {
        "A": sequence.count("A"),
        "U": sequence.count("U"),
        "G": sequence.count("G"),
        "C": sequence.count("C"),
    }


def find_gc_repeat(sequence: str) -> str | None:
    for substring in ("GGGG", "CCCC"):
        if substring in sequence:
            return substring
    return None


def lookup_thermo_value(key: str, lookup_db: Path = LOOKUP_DB) -> float | None:
    if not lookup_db.exists():
        return None
    with sqlite3.connect(lookup_db) as conn:
        row = conn.execute("SELECT value FROM lookup WHERE key=?", (key,)).fetchone()
    return float(row[0]) if row else None


def read_sirna_parameters(
    file_path: Path = SIRNA_PARAMETERS_FILE,
    overrides: dict[str, float] | None = None,
) -> dict[str, float]:
    """Return dsRIP siRNA scoring parameters with optional overrides."""

    params = dict(DEFAULT_SIRNA_PARAMETERS)
    if file_path.exists():
        with file_path.open("r", encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line or line.startswith("#") or "=" not in line:
                    continue
                key, value = line.split("=", 1)
                params[key.strip()] = float(value.strip())

    if overrides:
        for key, value in overrides.items():
            params[key] = float(value)
    return params


def normalize_score(
    value: float,
    min_val: float,
    max_val: float,
    target_min: float,
    target_max: float,
) -> float:
    return (value - min_val) / (max_val - min_val) * (target_max - target_min) + target_min


def add_total_score(
    features: dict[str, Any],
    *,
    safety: bool = False,
    safety_essential_coeff: int = 20,
    parameter_overrides: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Add dsRIP-style normalized feature scores and ``total_score`` in-place."""

    params = read_sirna_parameters(overrides=parameter_overrides)
    asym_weight = params["asymScore_weight"]
    self_fold_weight = params["self_fold_energy_weight"]
    combined_weight = asym_weight + self_fold_weight
    if combined_weight == 0:
        raise ValueError("asymScore_weight and self_fold_energy_weight cannot both be zero.")

    features["asymScore_normalized"] = normalize_score(
        float(features["new_asymScore"]), 5, -5, 100, 0
    )
    features["self_fold_energy_normalized"] = normalize_score(
        float(features["self_fold_energy"]), 0, -22, 100, 0
    )

    accessibility = float(features.get("accessibility", 0) or 0)
    if accessibility > 0:
        import math

        log_access = math.log10(accessibility)
        features["accessibility_normalized"] = normalize_score(
            log_access, math.log10(1e-10), math.log10(1), 0, 100
        )
    else:
        features["accessibility_normalized"] = 0

    base_score = (
        features["asymScore_normalized"] * asym_weight
        + features["self_fold_energy_normalized"] * self_fold_weight
    ) / combined_weight

    additional_score = 0.0
    gc_repeat = features.get("GGGG_CCCC") or features.get("G_C_repeat")
    if gc_repeat in {"GGGG", "CCCC"}:
        additional_score -= params["GGGG_CCCC_penalty"]
    if features.get("edge_asyym"):
        additional_score += params["edge_asyym_bonus"]
    if float(features.get("anti_GC", 0) or 0) <= params["anti_GC_upper_bound"]:
        additional_score += params["anti_GC_bonus"]
    if params["anti_9_14_GC_lower_bound"] <= float(features.get("anti_9_14_GC", 0) or 0):
        additional_score += params["anti_9_14_GC_bonus"]
    if features.get("ORF_info") == "ORF":
        additional_score += params["ORF_info_bonus"]
    if int(features.get("anti_10th_A", 0) or 0) == 1:
        additional_score += params["anti_10th_A_bonus"]

    features["total_score"] = round(base_score + additional_score, 1)

    if safety:
        all_off_targets = int(features.get("all_off_targets", 0) or 0)
        lethal_off_targets = int(features.get("lethal_off_targets", 0) or 0)
        weighted_off_targets = (all_off_targets - lethal_off_targets) + (
            lethal_off_targets * int(safety_essential_coeff)
        )
        features["off_target_weighted_count"] = weighted_off_targets

    return features


def thermo_asymmetry_score(antisense_5_3: str, sense_3_5: str, lookup_db: Path = LOOKUP_DB) -> float | None:
    antisense_value = lookup_thermo_value(antisense_5_3[:4], lookup_db)
    sense_value = lookup_thermo_value((sense_3_5[-4:])[::-1], lookup_db)
    if antisense_value is None or sense_value is None:
        return None
    return round(antisense_value - sense_value, 3)


def new_asymmetry_calc(antisense_5_3: str) -> float:
    positions = [0, 1, 2, 3, 4, 14, 15, 16, 17, 18]
    au_values = [1 if antisense_5_3[pos] in {"A", "U"} else 0 for pos in positions]
    diffs = [
        au_values[0] - au_values[9],
        au_values[1] - au_values[8],
        au_values[2] - au_values[7],
        au_values[3] - au_values[6],
        au_values[4] - au_values[5],
    ]
    weights = [1.6, 1.3, 0.3, 0.5, 0.6]
    return round(sum(difference * weight for difference, weight in zip(diffs, weights)), 5)


def fold_antisense(antisense_5_3: str, temperature: float) -> tuple[str, float]:
    try:
        import RNA
    except ModuleNotFoundError as exc:
        raise DsRipApiError("ViennaRNA Python bindings are required for self-fold calculation.") from exc

    md = RNA.md()
    md.temperature = temperature
    fc = RNA.fold_compound(antisense_5_3, md)
    structure, energy = fc.mfe()
    return structure, round(float(energy), 4)


def predict_sirna_features(
    antisense_sequence: str,
    *,
    accessibility: float = 1.0,
    orf_status: str = "ORF",
    sense_overhang: str = "TT",
    name: str = "siRNA",
    temperature: float = 25,
    all_off_targets: int | None = None,
    lethal_off_targets: int | None = None,
    calculate_total_score: bool = True,
    safety: bool = False,
    safety_essential_coeff: int = 20,
    parameter_overrides: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Return dsRIP per-siRNA feature calculations for one antisense sequence.

    The original dsRIP function needs gene-position context. For direct
    antisense input, ``accessibility`` and ``orf_status`` are supplied as
    explicit values. The sense strand is inferred as the RNA complement of the
    antisense sequence, matching dsRIP's ``sense_3_5_`` convention. Because the
    sense strand is stored 3' to 5', the sense 3' overhang is prefixed.
    """

    if orf_status not in {"ORF", "5_UTR", "3_UTR", "partial_ORF"}:
        raise ValueError("orf_status must be one of: ORF, 5_UTR, 3_UTR, partial_ORF.")

    antisense = normalize_antisense_sequence(antisense_sequence)
    overhang = "".join(str(sense_overhang).split()).upper().replace("U", "T")
    invalid_overhang = sorted(set(overhang) - set("ATCG"))
    if invalid_overhang:
        raise ValueError(f"Sense overhang contains invalid DNA base(s): {', '.join(invalid_overhang)}")

    sense_3_5 = overhang + complement_rna(antisense)
    antisense_9_14 = antisense[8:14]
    antisense_9_14_counts = count_nucleotides(antisense_9_14)
    antisense_counts = count_nucleotides(antisense)
    fold_structure, fold_energy = fold_antisense(antisense, temperature)

    result: dict[str, Any] = {
        "antisense_5_3": antisense,
        "sense_3_5_": sense_3_5,
        "antisense_5_3_DNA": antisense.replace("U", "T"),
        "anti_GC": round(((antisense_counts["C"] + antisense_counts["G"]) / len(antisense) * 100), 1),
        "anti_9_14_GC": round(((antisense_9_14_counts["C"] + antisense_9_14_counts["G"]) / 6 * 100), 1),
        "asymScore": thermo_asymmetry_score(antisense, sense_3_5),
        "new_asymScore": new_asymmetry_calc(antisense),
        "self_fold_energy": fold_energy,
        "accessibility": float(accessibility),
        "anti_10th_A": "1" if antisense[9] == "A" else "0",
        "ORF_info": orf_status,
        "G_C_repeat": find_gc_repeat(antisense),
        "GGGG_CCCC": find_gc_repeat(antisense),
        "edge_asyym": 1 if (antisense[0] in ["A", "U"]) and (antisense[18] == "G" or antisense[0] == "C") else 0,
        "self_fold_structure": fold_structure,
    }

    if all_off_targets is not None or lethal_off_targets is not None:
        result["all_off_targets"] = int(all_off_targets or 0)
        result["lethal_off_targets"] = int(lethal_off_targets or 0)

    if calculate_total_score:
        add_total_score(
            result,
            safety=safety,
            safety_essential_coeff=safety_essential_coeff,
            parameter_overrides=parameter_overrides,
        )

    return result


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Calculate dsRIP features for one antisense siRNA sequence.")
    parser.add_argument("antisense", help="Antisense siRNA sequence, 5' to 3'. DNA T is accepted and converted to U.")
    parser.add_argument("--accessibility", type=float, default=1.0, help="mRNA accessibility value. Default: 1.0.")
    parser.add_argument(
        "--orf-status",
        default="ORF",
        choices=["ORF", "5_UTR", "3_UTR", "partial_ORF"],
        help="Context label to report in ORF_info. Default: ORF.",
    )
    parser.add_argument("--sense-overhang", default="TT", help="Sense 3' overhang. Default: TT.")
    parser.add_argument("--name", default="siRNA", help="Name used internally for the one-siRNA record.")
    parser.add_argument("--temperature", type=float, default=25, help="RNAfold temperature. Default: 25.")
    parser.add_argument("--all-off-targets", type=int, help="Optional all-off-target count to include.")
    parser.add_argument("--lethal-off-targets", type=int, help="Optional lethal-off-target count to include.")
    parser.add_argument("--no-total-score", action="store_true", help="Skip dsRIP-style total_score calculation.")
    parser.add_argument("--safety", action="store_true", help="Include weighted off-target count if off-target columns are supplied.")
    parser.add_argument("--safety-essential-coeff", type=int, default=20, help="Lethal off-target coefficient. Default: 20.")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    try:
        result = predict_sirna_features(
            args.antisense,
            accessibility=args.accessibility,
            orf_status=args.orf_status,
            sense_overhang=args.sense_overhang,
            name=args.name,
            temperature=args.temperature,
            all_off_targets=args.all_off_targets,
            lethal_off_targets=args.lethal_off_targets,
            calculate_total_score=not args.no_total_score,
            safety=args.safety,
            safety_essential_coeff=args.safety_essential_coeff,
        )
    except (DsRipApiError, ValueError) as exc:
        parser.exit(1, f"Error: {exc}\n")

    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
