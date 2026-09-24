from collections import Counter
from pathlib import Path
from random import Random
from types import SimpleNamespace

import dsrna_adaptation as adaptation


def test_dinucleotide_shuffle_preserves_adjacent_pairs():
    sequence = "ACGTTGCAACGTTGCAAC"
    shuffled = adaptation.dinucleotide_shuffle(sequence, Random(8))
    assert len(shuffled) == len(sequence)
    assert Counter(zip(shuffled, shuffled[1:])) == Counter(zip(sequence, sequence[1:]))


def test_dinucleotide_shuffle_preserves_only_the_scored_guide_when_pos1_is_excluded():
    sequence = "AACGTTGCAACGTTGCAACGT"
    shuffled = adaptation.dinucleotide_shuffle(sequence, Random(8), preserve_prefix=1)
    assert shuffled[0] == sequence[0]
    assert Counter(zip(shuffled[1:], shuffled[2:])) == Counter(zip(sequence[1:], sequence[2:]))


def test_make_guides_generates_both_strands_and_controls():
    guides = adaptation.make_guides([("locus", "A" * 21)], guide_length=21, shuffle_count=3, seed=1)
    assert len(guides) == 8  # one position × sense/antisense × original + three controls
    assert {guide.strand for guide in guides} == {"sense", "antisense"}
    assert sum(guide.variant == "original" for guide in guides) == 2
    assert adaptation.normalized_shannon_entropy("A" * 21) == 0
    assert adaptation.normalized_shannon_entropy("ACGT" * 5 + "A") > 0.95


def test_mfe_adaptation_outputs_per_sirna_statistics(tmp_path, monkeypatch):
    dsrna = tmp_path / "dsrna.fa"
    focal = tmp_path / "focal.fa"
    control = tmp_path / "control.fa"
    dsrna.write_text(">locus_1\nACGTACGTACGTACGTACGTA\n", encoding="utf-8")
    focal.write_text(">focus_gene\nACGTACGTACGTACGTACGTACGTACGTACGT\n", encoding="utf-8")
    control.write_text(">control_gene\nACGTACGTACGTACGTACGTACGTACGTACGT\n", encoding="utf-8")

    def fake_find_targets(queries, transcripts, _config, _log):
        if "control_gene" in transcripts:
            return []
        hits = []
        for query in queries:
            if query.name.endswith("|original"):
                hits.append(
                    SimpleNamespace(
                        query=query.name, transcript="focus_gene", t_start=1, t_stop=21,
                        mfe_ratio=0.9, mfe_perfect=-30.0, mfe_site=-27.0, allen_score=0.0,
                        mismatch_count=0, gu_wobble_count=0, bulge_count=0, match_pattern="|" * 21,
                    )
                )
        return hits

    monkeypatch.setattr(adaptation, "find_targets", fake_find_targets)
    result = adaptation.run_dsrna_adaptation(
        adaptation.AdaptationConfig(
            dsrna_fasta=dsrna, focal_transcriptome_fasta=focal, control_transcriptome_fasta=control,
            output_dir=tmp_path / "out", guide_length=21, shuffle_count=10, fdr_cutoff=0.05,
        )
    )

    assert result["sirna_count"] == 2
    per_sirna = (tmp_path / "out" / "tables" / "dsrna_adaptation_per_sirna.tsv").read_text(encoding="utf-8")
    assert "pooled_pvalue_vs_shuffled" in per_sirna
    assert "fdr_focus_vs_control" in per_sirna
    assert "focus_gene" in per_sirna
    plot_dir = tmp_path / "out" / "plots"
    assert (plot_dir / "locus_1.adaptation_landscape.png").exists() or (plot_dir / "locus_1.adaptation_landscape.svg").exists()
