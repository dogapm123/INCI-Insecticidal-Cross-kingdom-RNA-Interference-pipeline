from __future__ import annotations

import glob
import os
import shutil
import subprocess
import tempfile


def dinucleotide_shuffling(sequence: str, number: int = 10, seed: int = 1) -> str:
    """
    Return dinucleotide-shuffled sequences as FASTA text.

    The input can be either a raw sequence:
        "AACGGAGGTATTCTATAGTTATA"

    or a FASTA record:
        ">OSR_dsNode343\\nAACGGAGGTATTCTATAGTTATA"
    """
    lines = [line.strip() for line in sequence.splitlines() if line.strip()]
    if not lines:
        raise ValueError("sequence must not be empty")

    if lines[0].startswith(">"):
        sequence_name = lines[0][1:].split()[0] or "sequence"
        sequence = "".join(lines[1:]).upper()
    else:
        sequence_name = "sequence"
        sequence = "".join(sequence.split()).upper()

    if not sequence:
        raise ValueError("sequence must not be empty")
    if number < 1:
        raise ValueError("number must be at least 1")

    command = (
        shutil.which("fasta-dinucleotide-shuffle-py3")
        or shutil.which("fasta-dinucleotide-shuffle")
    )
    if command is None and os.environ.get("MEME_BIN"):
        candidate = os.path.join(os.environ["MEME_BIN"], "fasta-dinucleotide-shuffle")
        if os.path.exists(candidate):
            command = candidate
    if command is None:
        matches = sorted(glob.glob("/opt/local/libexec/meme-*/fasta-dinucleotide-shuffle"))
        if matches:
            command = matches[-1]
    if command is None:
        raise FileNotFoundError("MEME Suite fasta-dinucleotide-shuffle was not found")

    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as fasta_file:
        fasta_file.write(f">{sequence_name}\n{sequence}\n")
        fasta_path = fasta_file.name

    try:
        result = subprocess.run(
            [command, "-f", fasta_path, "-c", str(number), "-s", str(seed)],
            check=True,
            capture_output=True,
            text=True,
        )
    finally:
        os.unlink(fasta_path)

    return result.stdout


fasta_output = dinucleotide_shuffling("AACGGAGGTATTCTATAGTTATAACCGATTCATACACAATCTGTACAAGAAGCAATTACTTCTTAATCGAAAAATACTTGCACAAATAGCTCTATTAAATAGGAGTTGTCTTTATACGATTTCAAATGAGATCAAAAAATGAGGGGATTGGAAGAAATCCACTAAAATATTTGAAATAGAGTTCTTAAGTAGAATAAGCTCGGGGAGGGTAGAGTAGAAATTAGTATTATAAAAAAAGTCGTCAAAACAAAACGAGGGTATATAGTCGCTAAAAAAAGACTTATTCTTTTTCTTTTCGACGATTCTTTTCTTTTTTTTTTTTTTTTTTTATGACAGAAACAGACGTAACAATCAAATTTTGATTCTTGATTGGATTTGTCGAACAAAATATTAACCTCTTTTTTAATTCCTTCAGATTAGTTATTCAATTGAAGAATAAGTCTATTTTTTTCTGGTTCTAAGGCTAGTAGTTCTAGGGGTCGACTCACTTCTTTCAAATTGTTTCTGATTATTAAGAAAAGGTAACAAAGATAAAATACGAGCTTGTTTTATAGCAATAGTAATTAATCGTTGTTGTTTTAAAGTTACTCTATTCACCCGTCTAGATAATATTTTTCCTTGTTCACTAATAAATCGACTAATTAAACTCATGTTTCTATAATCAATTCGATCCCCCGATTGGATCGGGGGCAAACGCCGACGAAAAGATCGCTTGGATTTAGTAAAAGGTCGCTTAGATTTATTCATAATTTTTTATTCCTTACCTATAAAATCCTATTTTTGATCGGATCAAAATAAAAACGATTGCTATTTCTATATAGGATAGAGTTTCTATATAGAATTAAGAATATAGGATTTGATTTCTTTCTATTATTTTATTTGTTTATATTTATATATATGTAATAATACGTAGATACTCTCATTATTCCTGTTCTTAAAATAAAATTTGACATACAAGCACTCAATTTTATCTATTTCTTAATTTCCCCATGAATTGTATGTTTATAACAATAGGGACAGAATTTTCTCAATTCCAATCGACTAGGAGTGTTATGCCGATTCTTTTGAGTAATATATCTGGAAATTCCAGCCGATTCTTTCTTAATATCATTTCGAACACAACTGGTACATTCCAAAATAATTGTTACTCGAACATCTTTACCTTTGGCCATGAAACCTCCTTTGGATTTATGATTCACCCCCAATCTTCTATTTTAATTTTAGGGTCCATAAAGAACAAGAACAAAAAAAATAGAAGCTGTTAATTAACTCAAAAAAATGTCTTAAATATTTCGTTGCAATGAATATTAATTAGTAATTAAATAAAAAGAAAATAATAAGTAATATAATGATTTTTTGAAAACTATATTTCAAATCTTTGTACAGTATTTTTTTTAAGTATCTCGTAGGATACTCTTTCCCTAGCAACACTTTTTTTTTAATTGGATCCTGAACCCTCGTTTTTATGACCTTTCGGCCCCCCTTTTTTTTTTCTAATTCTTTTTTTAGAATTAGAAAAAA",
                                      number=90)
print(fasta_output)