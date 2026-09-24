#!/usr/bin/env python3
"""Deduplicate nucleotide FASTA records with CD-HIT-EST."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
from pathlib import Path


def count_fasta_records(path: Path) -> int:
    count = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def auto_word_size(identity: float) -> int:
    if identity >= 0.95:
        return 10
    if identity >= 0.90:
        return 8
    if identity >= 0.88:
        return 7
    if identity >= 0.85:
        return 6
    if identity >= 0.80:
        return 5
    return 4


def parse_cluster_line(line: str) -> tuple[str, int, bool]:
    body = line.split("\t", 1)[1]
    if "nt, " in body:
        size_text, rest = body.split("nt, ", 1)
    else:
        size_text, rest = body.split("aa, ", 1)
    name = rest.split("...", 1)[0].lstrip(">")
    return name, int(size_text), line.rstrip().endswith("*")


def parse_cdhit_clusters(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    cluster_id = ""
    members: list[tuple[str, int, bool]] = []

    def flush() -> None:
        if not cluster_id or not members:
            return
        representative = next((name for name, _size, is_rep in members if is_rep), members[0][0])
        member_names = [name for name, _size, _is_rep in members]
        for name, size, is_rep in members:
            rows.append(
                {
                    "cluster_id": cluster_id,
                    "representative_id": representative,
                    "member_id": name,
                    "member_length": size,
                    "is_representative": "yes" if is_rep else "no",
                    "cluster_size": len(members),
                    "cluster_members": ";".join(member_names),
                }
            )

    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">Cluster "):
                flush()
                cluster_id = f"cluster_{line.rsplit(' ', 1)[-1]}"
                members = []
            else:
                members.append(parse_cluster_line(raw_line.rstrip("\n")))
    flush()
    return rows


def run_cdhit_est(args: argparse.Namespace, word_size: int) -> None:
    executable = args.cd_hit_est or shutil.which("cd-hit-est")
    if not executable:
        raise RuntimeError("cd-hit-est was not found on PATH. Install CD-HIT first.")

    command = [
        executable,
        "-i",
        str(args.input_fasta),
        "-o",
        str(args.output_fasta),
        "-c",
        f"{args.identity:.5f}".rstrip("0").rstrip("."),
        "-n",
        str(word_size),
        "-T",
        str(args.threads),
        "-M",
        str(args.memory_mb),
        "-d",
        str(args.description_length),
        "-r",
        "1" if args.both_strands else "0",
    ]
    if args.extra_args:
        command.extend(args.extra_args)

    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    args.log_file.parent.mkdir(parents=True, exist_ok=True)
    args.log_file.write_text((completed.stdout or "") + (completed.stderr or ""), encoding="utf-8")
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout or "").strip()
        raise RuntimeError(f"cd-hit-est failed with exit code {completed.returncode}: {detail}")


def positive_int(value: str) -> int:
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("value must be at least 1")
    return number


def identity_value(value: str) -> float:
    number = float(value)
    if number > 1.0:
        number = number / 100.0
    if not 0.75 <= number <= 1.0:
        raise argparse.ArgumentTypeError("identity must be between 0.75 and 1.0, or 75 and 100")
    return number


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-fasta", type=Path, required=True)
    parser.add_argument("--output-fasta", type=Path, required=True)
    parser.add_argument("--clusters-csv", type=Path, required=True)
    parser.add_argument("--log-file", type=Path, required=True)
    parser.add_argument("--manifest-json", type=Path)
    parser.add_argument("--identity", type=identity_value, default=0.95)
    parser.add_argument("--word-size", type=int, choices=range(4, 12))
    parser.add_argument("--threads", type=positive_int, default=4)
    parser.add_argument("--memory-mb", type=int, default=0, help="Memory limit passed to -M. Use 0 for CD-HIT default/unlimited.")
    parser.add_argument("--description-length", type=int, default=0, help="Description length passed to -d. Use 0 to keep full headers.")
    parser.add_argument("--same-strand-only", action="store_true", help="Only compare sequences on the same strand.")
    parser.add_argument("--cd-hit-est", help="Path to cd-hit-est executable.")
    parser.add_argument("--force", action="store_true", help="Overwrite existing output files.")
    parser.add_argument("extra_args", nargs=argparse.REMAINDER, help="Additional arguments passed after -- to cd-hit-est.")
    args = parser.parse_args()

    args.input_fasta = args.input_fasta.expanduser()
    args.output_fasta = args.output_fasta.expanduser()
    args.clusters_csv = args.clusters_csv.expanduser()
    args.log_file = args.log_file.expanduser()
    if args.manifest_json:
        args.manifest_json = args.manifest_json.expanduser()
    args.both_strands = not args.same_strand_only
    if args.extra_args and args.extra_args[0] == "--":
        args.extra_args = args.extra_args[1:]

    if not args.input_fasta.exists():
        raise FileNotFoundError(f"Input FASTA does not exist: {args.input_fasta}")
    if args.memory_mb < 0:
        raise ValueError("memory-mb must be 0 or greater.")
    if args.description_length < 0:
        raise ValueError("description-length must be 0 or greater.")

    clstr_path = Path(str(args.output_fasta) + ".clstr")
    protected_outputs = [args.output_fasta, clstr_path, args.clusters_csv, args.log_file]
    if args.manifest_json:
        protected_outputs.append(args.manifest_json)
    existing = [path for path in protected_outputs if path.exists()]
    if existing and not args.force:
        raise FileExistsError(f"Output already exists; use --force to overwrite: {existing[0]}")
    for path in existing:
        path.unlink()

    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.clusters_csv.parent.mkdir(parents=True, exist_ok=True)
    word_size = args.word_size or auto_word_size(args.identity)

    input_count = count_fasta_records(args.input_fasta)
    run_cdhit_est(args, word_size)
    if not clstr_path.exists():
        raise RuntimeError(f"CD-HIT cluster file was not created: {clstr_path}")

    rows = parse_cdhit_clusters(clstr_path)
    with args.clusters_csv.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "cluster_id",
            "representative_id",
            "member_id",
            "member_length",
            "is_representative",
            "cluster_size",
            "cluster_members",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    representative_count = count_fasta_records(args.output_fasta)
    manifest = {
        "input_fasta": str(args.input_fasta),
        "output_fasta": str(args.output_fasta),
        "cluster_file": str(clstr_path),
        "clusters_csv": str(args.clusters_csv),
        "log_file": str(args.log_file),
        "input_records": input_count,
        "representative_records": representative_count,
        "removed_records": max(0, input_count - representative_count),
        "identity": args.identity,
        "word_size": word_size,
        "threads": args.threads,
        "memory_mb": args.memory_mb,
        "both_strands": args.both_strands,
    }
    if args.manifest_json:
        args.manifest_json.parent.mkdir(parents=True, exist_ok=True)
        args.manifest_json.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")

    print(f"Input records: {input_count}")
    print(f"Representative records: {representative_count}")
    print(f"Removed records: {manifest['removed_records']}")
    print(f"Wrote {args.output_fasta}")
    print(f"Wrote {args.clusters_csv}")


if __name__ == "__main__":
    main()
