from __future__ import annotations

from pathlib import Path

import pytest

from pybarrnap import Barrnap
from pybarrnap.utils import load_example_fasta_file
from tests.marker import skipif_cmscan_not_installed


def test_bacteria_run(tmp_path: Path):
    """Test pybarrnap run for bacteria"""
    fasta_file = load_example_fasta_file("bacteria.fna")
    barrnap = Barrnap(fasta_file)
    result = barrnap.run()

    expected_rrna_count = 22
    assert len(result.get_rrna_seq_records()) == expected_rrna_count

    gff_file = tmp_path / "rrna.gff"
    result.write_gff(gff_file, incseq=False)
    assert gff_file.exists()

    gff_incseq_file = tmp_path / "rrna_incseq.gff"
    result.write_gff(gff_incseq_file, incseq=True)
    assert gff_incseq_file.exists()

    rrna_fasta_file = tmp_path / "rrna.fna"
    result.write_fasta(rrna_fasta_file)
    assert rrna_fasta_file.exists()


def test_archaea_run():
    """Test pybarrnap run for archaea"""
    fasta_file = load_example_fasta_file("archaea.fna")
    barrnap = Barrnap(fasta_file, kingdom="arc")
    result = barrnap.run()

    expected_rrna_count = 6
    assert len(result.get_rrna_seq_records()) == expected_rrna_count


def test_eukaryote_run():
    """Test pybarrnap run for eukaryote (fungus)"""
    fasta_file = load_example_fasta_file("fungus.fna")
    barrnap = Barrnap(fasta_file, kingdom="euk")
    result = barrnap.run()

    expected_rrna_count = 13
    assert len(result.get_rrna_seq_records()) == expected_rrna_count


def test_gzip_fasta_run():
    """Test barrnap run for gzip fasta file"""
    gzip_fasta_file = load_example_fasta_file("minimum.fna.gz")
    barrnap = Barrnap(gzip_fasta_file)
    result = barrnap.run()

    expected_rrna_count = 3
    assert len(result.get_rrna_seq_records()) == expected_rrna_count


def test_null_fasta_run_failed():
    """Test pybarrnap run for null fasta file (failed)"""
    null_fasta_file = load_example_fasta_file("null.fna")
    with pytest.raises(ValueError):
        Barrnap(null_fasta_file).run()


def test_dupid_fasta_run_failed():
    """Test pybarrnap run for duplication id fasta file (failed)"""
    dupid_fasta_file = load_example_fasta_file("dupid.fna")
    with pytest.raises(ValueError):
        Barrnap(dupid_fasta_file).run()


def test_protein_fasta_run_failed():
    """Test pybarrnap run for protein fasta file (failed)"""
    protein_fasta_file = load_example_fasta_file("protein.fna")
    with pytest.raises(ValueError):
        Barrnap(protein_fasta_file).run()


def test_empty_fasta_run():
    """Test pybarrnap run for empty (0 length) fasta file"""
    empty_fasta_file = load_example_fasta_file("empty.fna")
    barrnap = Barrnap(empty_fasta_file)
    result = barrnap.run()

    expected_rrna_count = 0
    assert len(result.get_rrna_seq_records()) == expected_rrna_count


def test_nohits_fasta_run():
    """Test pybarrnap run for nohits fasta file"""
    nohits_fasta_file = load_example_fasta_file("nohits.fna")
    barrnap = Barrnap(nohits_fasta_file)
    result = barrnap.run()

    expected_rrna_count = 0
    assert len(result.get_rrna_seq_records()) == expected_rrna_count


def test_accurate_with_no_cmscan_installed_run_failed(monkeypatch):
    """Test accurate option with no cmscan installed (failed)"""
    monkeypatch.setattr("pybarrnap.utils.is_cmscan_installed", lambda: False)
    fasta_file = load_example_fasta_file("minimum.fna")
    with pytest.raises(RuntimeError):
        Barrnap(fasta_file, accurate=True).run()


def test_kingdom_all_with_accurate_false_run_failed():
    """Test kingdom=`all` with accurate=False (failed)"""
    fasta_file = load_example_fasta_file("minimum.fna")
    with pytest.raises(ValueError):
        Barrnap(fasta_file, kingdom="all", accurate=False).run()


@skipif_cmscan_not_installed
def test_bacteria_accurate_run(tmp_path: Path):
    """Test pybarrnap accurate run for bacteria"""
    fasta_file = load_example_fasta_file("bacteria.fna")
    barrnap = Barrnap(fasta_file, accurate=True)
    result = barrnap.run()

    expected_rrna_count = 22
    assert len(result.get_rrna_seq_records()) == expected_rrna_count

    gff_file = tmp_path / "rrna.gff"
    result.write_gff(gff_file, incseq=False)
    assert gff_file.exists()

    gff_incseq_file = tmp_path / "rrna_incseq.gff"
    result.write_gff(gff_incseq_file, incseq=True)
    assert gff_incseq_file.exists()

    rrna_fasta_file = tmp_path / "rrna.fna"
    result.write_fasta(rrna_fasta_file)
    assert rrna_fasta_file.exists()


@skipif_cmscan_not_installed
def test_archaea_accurate_run():
    """Test pybarrnap accurate run for archaea"""
    fasta_file = load_example_fasta_file("archaea.fna")
    barrnap = Barrnap(fasta_file, kingdom="arc", accurate=True)
    result = barrnap.run()

    expected_rrna_count = 6
    assert len(result.get_rrna_seq_records()) == expected_rrna_count


@skipif_cmscan_not_installed
def test_accurate_all_run():
    """Test pybarrnap accurate run for all kingdom"""
    fasta_file = load_example_fasta_file("minimum.fna")
    barrnap = Barrnap(fasta_file, kingdom="all", accurate=True)
    result = barrnap.run()

    expected_rrna_count = 3
    assert len(result.get_rrna_seq_records()) == expected_rrna_count
