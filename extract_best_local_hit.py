#!/usr/bin/env python3
"""
Extract the best local target sequence hit for a short DNA query.

The script maps a short query sequence, usually around 2-3 kb, against a
larger FASTA file with minimap2. It then extracts only the target interval
covered by the best mapping, rather than the whole contig, using bedtools.
"""

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


DNA_RE = re.compile(r"^[ACGTRYSWKMBDHVNacgtryswkmbdhvn\s]+$")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Map a short DNA sequence to a larger FASTA and extract the best "
            "local matching target stretch."
        )
    )
    parser.add_argument(
        "query",
        help=(
            "Query FASTA path, or a raw DNA sequence. Raw sequences are useful "
            "for short 2-3 kb inputs copied directly into the command."
        ),
    )
    parser.add_argument("target_fasta", help="Larger FASTA file to search against.")
    parser.add_argument(
        "-o",
        "--output",
        default="best_local_hit.fasta",
        help="Output FASTA path for the extracted local target sequence.",
    )
    parser.add_argument(
        "--paf",
        default="best_local_hit.minimap2.paf",
        help="Path to write all minimap2 PAF alignments.",
    )
    parser.add_argument(
        "--bed",
        default="best_local_hit.bed",
        help="Path to write the selected BED interval.",
    )
    parser.add_argument(
        "--flank",
        type=int,
        default=0,
        help="Extra bases to include on each side of the mapped interval.",
    )
    parser.add_argument(
        "--preset",
        default="asm5",
        help="minimap2 preset to use. For similar DNA sequences, asm5 is a good default.",
    )
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help="Minimum mapping quality required for the selected hit.",
    )
    parser.add_argument(
        "--target-orientation",
        action="store_true",
        help=(
            "Keep the extracted sequence in target FASTA orientation. By default, "
            "bedtools -s is used so minus-strand hits are reverse-complemented."
        ),
    )
    parser.add_argument(
        "--minimap2",
        default="minimap2",
        help="Path to minimap2 executable.",
    )
    parser.add_argument(
        "--bedtools",
        default="bedtools",
        help="Path to bedtools executable.",
    )
    return parser.parse_args()


def require_executable(executable):
    if shutil.which(executable) or Path(executable).exists():
        return
    raise FileNotFoundError(f"Required executable not found: {executable}")


def existing_path(value):
    try:
        return Path(value) if Path(value).exists() else None
    except OSError:
        return None


def write_query_fasta(query_arg, work_dir):
    query_path = existing_path(query_arg)
    if query_path is not None:
        return query_path

    if not DNA_RE.match(query_arg):
        raise ValueError(
            "Query is neither an existing FASTA path nor a valid raw DNA sequence."
        )

    sequence = re.sub(r"\s+", "", query_arg).upper()
    if not sequence:
        raise ValueError("Query sequence is empty.")

    temp_query = Path(work_dir) / "query.fasta"
    with temp_query.open("w") as handle:
        handle.write(">query\n")
        for i in range(0, len(sequence), 80):
            handle.write(sequence[i : i + 80] + "\n")
    return temp_query


def run_command(command):
    completed = subprocess.run(
        command,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if completed.returncode != 0:
        sys.stderr.write(completed.stderr)
        raise subprocess.CalledProcessError(
            completed.returncode,
            command,
            output=completed.stdout,
            stderr=completed.stderr,
        )
    return completed


def run_minimap2(args, query_fasta, target_fasta):
    command = [
        args.minimap2,
        "-x",
        args.preset,
        str(target_fasta),
        str(query_fasta),
    ]
    result = run_command(command)
    paf_path = Path(args.paf)
    paf_path.write_text(result.stdout)
    return paf_path


def paf_fields(line):
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 12:
        raise ValueError(f"Invalid PAF line with fewer than 12 fields: {line!r}")
    return {
        "query_name": fields[0],
        "query_length": int(fields[1]),
        "query_start": int(fields[2]),
        "query_end": int(fields[3]),
        "strand": fields[4],
        "target_name": fields[5],
        "target_length": int(fields[6]),
        "target_start": int(fields[7]),
        "target_end": int(fields[8]),
        "matching_bases": int(fields[9]),
        "alignment_block_length": int(fields[10]),
        "mapping_quality": int(fields[11]),
        "raw": line.rstrip("\n"),
    }


def choose_best_hit(paf_path, min_mapq):
    hits = []
    with paf_path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            hit = paf_fields(line)
            if hit["mapping_quality"] >= min_mapq:
                hits.append(hit)

    if not hits:
        raise RuntimeError("No minimap2 hits passed the mapping-quality filter.")

    return sorted(
        hits,
        key=lambda hit: (
            hit["mapping_quality"],
            hit["matching_bases"],
            hit["alignment_block_length"],
            hit["query_end"] - hit["query_start"],
        ),
        reverse=True,
    )[0]


def write_bed(best_hit, bed_path, flank):
    start = max(0, best_hit["target_start"] - flank)
    end = min(best_hit["target_length"], best_hit["target_end"] + flank)
    name = (
        f"{best_hit['query_name']}|{best_hit['target_name']}:"
        f"{start + 1}-{end}({best_hit['strand']})"
    )
    score = best_hit["mapping_quality"]
    bed_line = "\t".join(
        [
            best_hit["target_name"],
            str(start),
            str(end),
            name,
            str(score),
            best_hit["strand"],
        ]
    )
    Path(bed_path).write_text(bed_line + "\n")


def run_bedtools(args, target_fasta):
    command = [
        args.bedtools,
        "getfasta",
        "-fi",
        str(target_fasta),
        "-bed",
        str(args.bed),
        "-fo",
        str(args.output),
        "-name",
    ]
    if not args.target_orientation:
        command.append("-s")
    run_command(command)


def main():
    args = parse_args()
    target_fasta = Path(args.target_fasta)

    if not target_fasta.exists():
        raise FileNotFoundError(f"Target FASTA not found: {target_fasta}")
    if args.flank < 0:
        raise ValueError("--flank must be 0 or greater.")

    require_executable(args.minimap2)
    require_executable(args.bedtools)

    with tempfile.TemporaryDirectory() as work_dir:
        query_fasta = write_query_fasta(args.query, work_dir)
        paf_path = run_minimap2(args, query_fasta, target_fasta)
        best_hit = choose_best_hit(paf_path, args.min_mapq)

    write_bed(best_hit, args.bed, args.flank)
    run_bedtools(args, target_fasta)

    print("Best local hit extracted.")
    print(f"Target: {best_hit['target_name']}")
    print(
        "Target interval: "
        f"{best_hit['target_start'] + 1}-{best_hit['target_end']} "
        f"({best_hit['strand']} strand)"
    )
    print(
        "Query interval: "
        f"{best_hit['query_start'] + 1}-{best_hit['query_end']} "
        f"of {best_hit['query_length']} bp"
    )
    print(f"Mapping quality: {best_hit['mapping_quality']}")
    print(f"PAF alignments: {args.paf}")
    print(f"BED interval: {args.bed}")
    print(f"Extracted FASTA: {args.output}")


if __name__ == "__main__":
    main()

