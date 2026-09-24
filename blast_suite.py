import subprocess
from tempfile import NamedTemporaryFile
from pathlib import Path
import csv
import os

# -----------------------------
# Inputs
# -----------------------------
sequence = "GCAAACGCCGACGAAAAGATCGCTTGGATTTAGTAAAAGGTCGCTTAGATTTATTCATAATTTTTTATTCCTTACCTATAAAATCCTATTTTTGATCGGATCAAAATAAAAACGATTGCTATTTCTATATAGGATAGAGTTTCTATATAGAATTAAGAATATAGGATTTGATTTCTTTCTATTATTTTATTTGTTTATATTTATATATATGTAATAATACGTAGATACTCTCATTATTCCTGTTCTTAAAATAAAATTTGACATACAAGCACTCAATTTTATCTATTTCTTAATTTCCCCATGAATTGTATGTTTATAACAATAGGGACAGAATTTTCTCAATTCCAATCGACTAGGAGTGTTATGCCGATTCTTTTGAGTAATATATCTGGAAATTCCAGCCGATTCTTTCTTAATATCATTTCGAACACAACTGGTACATTCCAAAATAATTGTTACTCGAACATCTTTACCTTTGGCCATGAAACCTCCTTTGGATTTATGATTCACCCCCAATCTTCTATTTTAATTTTAGGGTCCATAAAGAACAAGAACAAAAAAAATAGAAGCTGTTAATTAACTCAAAAAAATGTCTTAAATATTTCGTTGCAATGAATATTAATTAGTAATTAAATAAAAAGAAAATAATAAGTAATATAATGATTTTTTGAAAACTATATTTCAAATCTTTGTACAGTATTTTTTTTAAGTATCTCGTAGGATACTCTTTCCCTAGCAACACTTTTTTTTTAATTGGATCCTGAACCCTCGTTTTTATGACCTTTCGGCCCCCCTTTTTTTTTTCTAATTCTTTTTTTAGAATTAGAAAAAA"  # put your query sequence here
database = Path("/Users/dogacedden/Desktop/Postdoc/natural_RNAi_postdoc/transcriptomes/Pchr_corr.fasta")

makeblastdb_path = Path("external_tools/ncbi-blast-2.17.0+/bin/makeblastdb")
blastn_path = Path("external_tools/ncbi-blast-2.17.0+/bin/blastn")


db_name = "mydb_2"
results_file = Path("results.tsv")
summary_file = Path("blast_summary.txt")

# -----------------------------
# Validate input
# -----------------------------
def load_fasta_sequences(fasta_path):
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)

                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())

        if current_id is not None:
            sequences[current_id] = "".join(current_seq)

    return sequences


# -----------------------------
# Validate input
# -----------------------------
sequence = sequence.replace(" ", "").replace("\n", "").upper()

if not sequence:
    raise ValueError("Query sequence is empty.")

if not database.exists():
    raise FileNotFoundError(f"Database FASTA not found: {database}")

if not makeblastdb_path.exists():
    raise FileNotFoundError(f"makeblastdb not found: {makeblastdb_path}")

if not blastn_path.exists():
    raise FileNotFoundError(f"blastn not found: {blastn_path}")


# -----------------------------
# Create temporary query FASTA
# -----------------------------
with NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as tmp:
    tmp.write(">query\n")
    tmp.write(sequence + "\n")
    query_path = tmp.name


try:
    # -----------------------------
    # Create BLAST database
    # -----------------------------
    subprocess.run([
        str(makeblastdb_path),
        "-in", str(database),
        "-dbtype", "nucl",
        "-out", db_name
    ], check=True)

    # -----------------------------
    # Run BLAST
    # -----------------------------
    outfmt_fields = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]

    subprocess.run([
        str(blastn_path),
        "-query", query_path,
        "-db", db_name,
        "-out", str(results_file),
        "-outfmt", "6 " + " ".join(outfmt_fields)
    ], check=True)

    # -----------------------------
    # Add headers to results.tsv
    # -----------------------------
    with open(results_file, "r") as original:
        data = original.read()

    with open(results_file, "w") as modified:
        modified.write("\t".join(outfmt_fields) + "\n")
        modified.write(data)

    # -----------------------------
    # Parse BLAST results
    # -----------------------------
    hits = []

    with open(results_file, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            row["pident"] = float(row["pident"])
            row["length"] = int(row["length"])
            row["mismatch"] = int(row["mismatch"])
            row["gapopen"] = int(row["gapopen"])
            row["qstart"] = int(row["qstart"])
            row["qend"] = int(row["qend"])
            row["sstart"] = int(row["sstart"])
            row["send"] = int(row["send"])
            row["evalue"] = float(row["evalue"])
            row["bitscore"] = float(row["bitscore"])
            hits.append(row)

    # -----------------------------
    # Load full database sequences
    # -----------------------------
    database_sequences = load_fasta_sequences(database)

    # -----------------------------
    # Print top 3 BLAST results
    # -----------------------------
    top_hits = sorted(hits, key=lambda x: (x["evalue"], -x["bitscore"]))[:3]

    if top_hits:
        print("\nTop 3 BLAST results:")

        for i, hit in enumerate(top_hits, start=1):
            subject_id = hit["sseqid"]
            full_subject_sequence = database_sequences.get(subject_id, "Sequence not found")

            print(f"\nHit {i}:")
            print(f"  Subject ID: {subject_id}")
            print(f"  Percent identity: {hit['pident']}%")
            print(f"  Alignment length: {hit['length']} bp")
            print(f"  E-value: {hit['evalue']}")
            print(f"  Bit score: {hit['bitscore']}")
            print(f"  Query range: {hit['qstart']}-{hit['qend']}")
            print(f"  Subject range: {hit['sstart']}-{hit['send']}")
            print(f"  Full subject sequence: {full_subject_sequence}")
    else:
        print("No BLAST hits found.")

    # -----------------------------
    # Summarize results
    # -----------------------------
    with open(summary_file, "w") as out:
        out.write("BLAST Summary\n")
        out.write("=============\n\n")
        out.write(f"Database FASTA: {database}\n")
        out.write(f"Query length: {len(sequence)} bp\n")
        out.write(f"Results file: {results_file}\n\n")

        if not hits:
            out.write("No BLAST hits found.\n")
        else:
            best_hit = top_hits[0]
            best_subject_id = best_hit["sseqid"]
            best_full_sequence = database_sequences.get(best_subject_id, "Sequence not found")

            out.write(f"Number of hits: {len(hits)}\n\n")
            out.write("Best hit:\n")
            out.write(f"  Subject ID: {best_subject_id}\n")
            out.write(f"  Percent identity: {best_hit['pident']}%\n")
            out.write(f"  Alignment length: {best_hit['length']} bp\n")
            out.write(f"  E-value: {best_hit['evalue']}\n")
            out.write(f"  Bit score: {best_hit['bitscore']}\n")
            out.write(f"  Query range: {best_hit['qstart']}-{best_hit['qend']}\n")
            out.write(f"  Subject range: {best_hit['sstart']}-{best_hit['send']}\n")
            out.write(f"  Full subject sequence: {best_full_sequence}\n")

            print("\nBLAST completed.")
            print(f"Number of hits: {len(hits)}")
            print(f"Best hit: {best_subject_id}")
            print(f"Identity: {best_hit['pident']}%")
            print(f"E-value: {best_hit['evalue']}")
            print(f"Full results written to: {results_file}")
            print(f"Summary written to: {summary_file}")

finally:
    if os.path.exists(query_path):
        os.remove(query_path)