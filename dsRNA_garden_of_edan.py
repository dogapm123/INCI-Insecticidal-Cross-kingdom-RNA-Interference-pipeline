"""Generate and process every siRNA duplex represented by a dsRNA.

This is a dependency-free rewrite of the duplex-generation part of dsRIP's
``dsRNA_efficiency.generate_siRNA``.  It keeps the original 2-nt 3' overhang
geometry while making strand orientation and 1-based coordinates explicit.
The input sequence is the sense strand of the dsRNA, written 5' to 3'.  DNA or
RNA input is accepted and normalized to RNA.  For a conventional 21-nt siRNA,
each duplex has a 19-bp paired region and a 2-nt 3' overhang on each strand.

Example
-------
    python dsRNA_garden_of_edan.py \
        --name dsRNA_1 \
        --sequence AUGCAUGCAUGCAUGCAUGCAUGCAUGC

FASTA input is also supported:

    python dsRNA_garden_of_edan.py --fasta dsRNAs.fasta --output garden.json
"""

from __future__ import annotations

import argparse
import json
import re
from collections.abc import Callable, Iterable, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, TypeVar


RNA_BASES = frozenset("ACGU")
DEFAULT_SIRNA_LENGTH = 21
DEFAULT_OVERHANG_LENGTH = 2
BOUNDARY_BASE = "X"


class DsRNAInputError(ValueError):
    """Raised when a dsRNA name or sequence cannot be used."""


@dataclass(frozen=True)
class DsRNARecord:
    """One named dsRNA sense strand, normalized to RNA and oriented 5' to 3'."""

    name: str
    sequence_5_to_3: str

    def __post_init__(self) -> None:
        clean_name = str(self.name).strip()
        if not clean_name:
            raise DsRNAInputError("The dsRNA name must not be empty.")
        object.__setattr__(self, "name", clean_name)
        object.__setattr__(self, "sequence_5_to_3", normalize_rna(self.sequence_5_to_3))

    @classmethod
    def from_sequence(cls, name: str, sequence: str) -> "DsRNARecord":
        return cls(name, sequence)


@dataclass(frozen=True)
class SiRNADuplex:
    """A dsRIP-style siRNA duplex with explicit strand directions.

    ``position_start`` and ``position_end`` are 1-based, inclusive positions of
    the antisense guide's target window on the input dsRNA sense strand.  The
    passenger/sense strand is shifted downstream by ``overhang_length`` bases.
    At the right boundary, dsRIP uses ``X`` placeholders for unavailable bases;
    ``is_complete`` makes those two terminal cases easy to filter if desired.
    """

    sirna_id: str
    number: int
    dsrna_name: str
    position_start: int
    position_end: int
    sense_position_start: int
    sense_position_end: int
    antisense_5_to_3: str
    sense_5_to_3: str
    sense_3_to_5: str
    overhang_length: int
    paired_length: int
    is_complete: bool

    @property
    def antisense_3prime_overhang(self) -> str:
        return self.antisense_5_to_3[-self.overhang_length :]

    @property
    def sense_3prime_overhang(self) -> str:
        return self.sense_3_to_5[: self.overhang_length]

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["antisense_3prime_overhang"] = self.antisense_3prime_overhang
        result["sense_3prime_overhang"] = self.sense_3prime_overhang
        return result


@dataclass(frozen=True)
class DsRNAGardenResult:
    """The dsRNA, all generated duplexes, and results from their processor."""

    dsrna: DsRNARecord
    sirna_length: int
    overhang_length: int
    duplexes: tuple[SiRNADuplex, ...]
    processed_duplexes: tuple[Any, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "dsrna": asdict(self.dsrna),
            "sirna_length": self.sirna_length,
            "overhang_length": self.overhang_length,
            "duplex_count": len(self.duplexes),
            "sirna_duplexes": [duplex.to_dict() for duplex in self.duplexes],
            "processed_duplexes": list(self.processed_duplexes),
        }


ProcessorResult = TypeVar("ProcessorResult")
DuplexProcessor = Callable[[SiRNADuplex], ProcessorResult]


def normalize_rna(sequence: str) -> str:
    """Return uppercase RNA after removing whitespace and validating bases."""

    if sequence is None:
        raise DsRNAInputError("The dsRNA sequence must not be None.")

    normalized = re.sub(r"\s+", "", str(sequence)).upper().replace("T", "U")
    if not normalized:
        raise DsRNAInputError("The dsRNA sequence must not be empty.")

    invalid = sorted(set(normalized) - RNA_BASES)
    if invalid:
        raise DsRNAInputError(
            "The dsRNA sequence contains invalid base(s): " + ", ".join(invalid)
        )
    return normalized


def complement_rna(sequence: str) -> str:
    """Return the RNA complement without changing strand direction."""

    return sequence.translate(str.maketrans("ACGU", "UGCA"))


def reverse_complement_rna(sequence: str) -> str:
    """Return the RNA reverse complement."""

    return complement_rna(sequence)[::-1]


def _identifier_name(name: str) -> str:
    """Make a stable identifier component while retaining readable names."""

    identifier = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip()).strip("_")
    return identifier or "dsRNA"


def generate_sirna_duplexes(
    dsrna: DsRNARecord,
    *,
    sirna_length: int = DEFAULT_SIRNA_LENGTH,
    overhang_length: int = DEFAULT_OVERHANG_LENGTH,
) -> tuple[SiRNADuplex, ...]:
    """Generate every dsRIP-style siRNA duplex in 5'-to-3' input order.

    This is algebraically equivalent to dsRIP's padded 23-nt sliding windows:
    the antisense guide targets ``start..end`` and the sense/passenger strand
    comes from ``start + overhang..end + overhang``.  Missing bases beyond the
    dsRNA's right boundary are represented by ``X``, as in the original code.
    """

    if isinstance(sirna_length, bool) or not isinstance(sirna_length, int):
        raise TypeError("sirna_length must be an integer.")
    if isinstance(overhang_length, bool) or not isinstance(overhang_length, int):
        raise TypeError("overhang_length must be an integer.")
    if sirna_length < 1:
        raise ValueError("sirna_length must be at least 1.")
    if not 0 < overhang_length < sirna_length:
        raise ValueError("overhang_length must be between 1 and sirna_length - 1.")
    if len(dsrna.sequence_5_to_3) < sirna_length:
        raise DsRNAInputError(
            f"dsRNA '{dsrna.name}' is {len(dsrna.sequence_5_to_3)} nt long; "
            f"at least {sirna_length} nt are required."
        )

    sequence = dsrna.sequence_5_to_3
    identifier_name = _identifier_name(dsrna.name)
    duplexes: list[SiRNADuplex] = []

    for zero_based_start in range(len(sequence) - sirna_length + 1):
        number = zero_based_start + 1
        position_start = number
        position_end = zero_based_start + sirna_length

        target_window = sequence[zero_based_start:position_end]
        antisense = reverse_complement_rna(target_window)

        sense_start = zero_based_start + overhang_length
        sense_end = sense_start + sirna_length
        available_sense = sequence[sense_start:sense_end]
        missing_bases = sirna_length - len(available_sense)
        sense_5_to_3 = available_sense + (BOUNDARY_BASE * missing_bases)
        sense_3_to_5 = sense_5_to_3[::-1]

        duplexes.append(
            SiRNADuplex(
                sirna_id=(
                    f"{identifier_name}_si_pos_{position_start}-{position_end}"
                ),
                number=number,
                dsrna_name=dsrna.name,
                position_start=position_start,
                position_end=position_end,
                sense_position_start=position_start + overhang_length,
                sense_position_end=position_end + overhang_length,
                antisense_5_to_3=antisense,
                sense_5_to_3=sense_5_to_3,
                sense_3_to_5=sense_3_to_5,
                overhang_length=overhang_length,
                paired_length=sirna_length - overhang_length,
                is_complete=missing_bases == 0,
            )
        )

    return tuple(duplexes)


def dummy_process_sirna_duplex(duplex: SiRNADuplex) -> dict[str, Any]:
    """Example processing hook; replace its body with a real duplex analysis."""

    gc_count = duplex.antisense_5_to_3.count("G") + duplex.antisense_5_to_3.count("C")
    return {
        "sirna_id": duplex.sirna_id,
        "status": "dummy_processed",
        "antisense_gc_percent": round(100 * gc_count / len(duplex.antisense_5_to_3), 2),
        "complete_duplex": duplex.is_complete,
    }


def process_all_sirna_duplexes(
    duplexes: Iterable[SiRNADuplex],
    processor: DuplexProcessor[ProcessorResult] = dummy_process_sirna_duplex,
) -> tuple[ProcessorResult, ...]:
    """Run one processor function over every duplex, preserving input order."""

    if not callable(processor):
        raise TypeError("processor must be callable.")
    return tuple(processor(duplex) for duplex in duplexes)


def garden_of_edan(
    dsrna: DsRNARecord,
    *,
    sirna_length: int = DEFAULT_SIRNA_LENGTH,
    overhang_length: int = DEFAULT_OVERHANG_LENGTH,
    processor: DuplexProcessor[Any] = dummy_process_sirna_duplex,
) -> DsRNAGardenResult:
    """Generate every duplex for one dsRNA and process each with ``processor``."""

    duplexes = generate_sirna_duplexes(
        dsrna,
        sirna_length=sirna_length,
        overhang_length=overhang_length,
    )
    processed = process_all_sirna_duplexes(duplexes, processor)
    return DsRNAGardenResult(
        dsrna=dsrna,
        sirna_length=sirna_length,
        overhang_length=overhang_length,
        duplexes=duplexes,
        processed_duplexes=processed,
    )


def read_fasta(path: str | Path) -> tuple[DsRNARecord, ...]:
    """Read one or more dsRNAs from FASTA without external dependencies."""

    fasta_path = Path(path)
    records: list[DsRNARecord] = []
    name: str | None = None
    sequence_parts: list[str] = []
    seen_names: set[str] = set()

    def finish_record() -> None:
        nonlocal name, sequence_parts
        if name is None:
            return
        if name in seen_names:
            raise DsRNAInputError(f"Duplicate FASTA record name: {name}")
        records.append(DsRNARecord.from_sequence(name, "".join(sequence_parts)))
        seen_names.add(name)

    try:
        with fasta_path.open("r", encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    finish_record()
                    header = line[1:].strip()
                    if not header:
                        raise DsRNAInputError(
                            f"Empty FASTA header at line {line_number} in {fasta_path}."
                        )
                    name = header.split()[0]
                    sequence_parts = []
                elif name is None:
                    raise DsRNAInputError(
                        f"Sequence found before the first FASTA header at line {line_number}."
                    )
                else:
                    sequence_parts.append(line)
    except OSError as exc:
        raise DsRNAInputError(f"Could not read FASTA file '{fasta_path}': {exc}") from exc

    finish_record()
    if not records:
        raise DsRNAInputError(f"No FASTA records found in '{fasta_path}'.")
    return tuple(records)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate all dsRIP-style siRNA duplexes for one or more dsRNAs."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--sequence", help="One dsRNA sense-strand sequence, 5' to 3'.")
    source.add_argument("--fasta", type=Path, help="FASTA file containing one or more dsRNAs.")
    parser.add_argument("--name", default="dsRNA", help="Name used with --sequence.")
    parser.add_argument("--sirna-length", type=int, default=DEFAULT_SIRNA_LENGTH)
    parser.add_argument("--overhang-length", type=int, default=DEFAULT_OVERHANG_LENGTH)
    parser.add_argument("--output", type=Path, help="JSON output file; default is standard output.")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    try:
        records = (
            read_fasta(args.fasta)
            if args.fasta is not None
            else (DsRNARecord.from_sequence(args.name, args.sequence),)
        )
        gardens = [
            garden_of_edan(
                record,
                sirna_length=args.sirna_length,
                overhang_length=args.overhang_length,
            ).to_dict()
            for record in records
        ]
    except (DsRNAInputError, TypeError, ValueError) as exc:
        parser.exit(2, f"Error: {exc}\n")

    payload: dict[str, Any] = {
        "garden_count": len(gardens),
        "gardens": gardens,
    }
    rendered = json.dumps(payload, indent=2) + "\n"

    if args.output is None:
        print(rendered, end="")
    else:
        try:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(rendered, encoding="utf-8")
        except OSError as exc:
            parser.exit(1, f"Error: could not write '{args.output}': {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
