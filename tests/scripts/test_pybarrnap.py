from __future__ import annotations

import subprocess as sp
from pathlib import Path

import pytest

import pybarrnap
from pybarrnap.utils import load_example_fasta_file


def test_cli_with_min_option():
    """Test cli default run"""
    fasta_file = load_example_fasta_file("bacteria.fna")
    cmd = f"pybarrnap {fasta_file}"
    result = sp.run(cmd, shell=True, capture_output=True)

    if result.returncode != 0:
        pytest.fail(str(result))


def test_cli_with_all_option(tmp_path: Path):
    """Test cli with all option"""
    fasta_file = load_example_fasta_file("mitochondria.fna")
    rrna_outfile = tmp_path / "rrna.fna"

    cmd = f"pybarrnap {fasta_file} "
    cmd += f"-e 1e-6 -l 0.8 -r 0.25 -t 1 -k mito -o {rrna_outfile} -i -q"
    result = sp.run(cmd, shell=True, capture_output=True)

    if result.returncode != 0:
        pytest.fail(str(result))

    assert rrna_outfile.exists()


def test_cli_with_invalid_option_value():
    """Test cli with invalid option value"""
    fasta_file = load_example_fasta_file("mitochondria.fna")
    # Invalid evalue (expected: 'v > 0')
    cmd = f"pybarrnap {fasta_file} --evalue 0"
    result = sp.run(cmd, shell=True)
    assert result.returncode == 2
    # Invalid lencutoff (expected: '0 <= v <= 1')
    cmd = f"pybarrnap {fasta_file} --lencutoff 10"
    result = sp.run(cmd, shell=True)
    assert result.returncode == 2
    # Invalid reject (expected '0 <= v <= 1')
    cmd = f"pybarrnap {fasta_file} --reject 10"
    result = sp.run(cmd, shell=True)
    assert result.returncode == 2
    # Invalid threads (expected '1 <= v <= max_threads')
    cmd = f"pybarrnap {fasta_file} --threads -1"
    result = sp.run(cmd, shell=True)
    assert result.returncode == 2


def test_cli_fasta_input_cases(tmp_path: Path):
    """Test cli fasta input cases"""
    fasta_file = load_example_fasta_file("mitochondria.fna")
    # Case1: pipe
    cmd = f"cat {fasta_file} | pybarrnap"
    result = sp.run(cmd, shell=True)
    assert result.returncode == 0
    # Case2: pipe with '-'
    cmd = f"cat {fasta_file} | pybarrnap -"
    result = sp.run(cmd, shell=True)
    assert result.returncode == 0
    # Case3: stdin with '<'
    cmd = f"pybarrnap < {fasta_file} > {tmp_path / 'rrna.gff'}"
    result = sp.run(cmd, shell=True)
    assert result.returncode == 0


def test_cli_version():
    """Test cli -v option"""
    result = sp.run("pybarrnap -v", shell=True, capture_output=True, text=True)
    assert result.stdout.rstrip() == f"v{pybarrnap.__version__}"
    assert result.returncode == 0
