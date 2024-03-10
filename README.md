# pybarrnap: Python implementation of barrnap

![Python3](https://img.shields.io/badge/Language-Python3-steelblue)
![OS](https://img.shields.io/badge/OS-_Mac_|_Linux-steelblue)
![License](https://img.shields.io/badge/license-GPLv3-blue)
[![Latest PyPI version](https://img.shields.io/pypi/v/pybarrnap.svg)](https://pypi.python.org/pypi/pybarrnap)
[![CI](https://github.com/moshi4/pybarrnap/actions/workflows/ci.yml/badge.svg)](https://github.com/moshi4/pybarrnap/actions/workflows/ci.yml)

## Table of contents

- [Overview](#overview)
- [Installation](#installation)
- [CLI Usage](#cli-usage)
- [API Usage](#api-usage)
- [LICENSE](#license)

## Overview

pybarrnap is a python implementation of [barrnap](https://github.com/tseemann/barrnap) (Bacterial ribosomal RNA predictor).
pybarrnap provides a CLI compatible with barrnap and also provides a python API for running rRNA prediction and retrieving predicted rRNA.
pybarrnap default mode depends only on the python library and not on the external command-line tools nhmmer and bedtools.
As an additional feature from barrnap, accurate mode is available by installing the external command line tool cmscan([infernal](http://eddylab.org/infernal/)).

> [!NOTE]
> Barrnap v0.9 uses the HMM profile database created from older releases of Rfam(11.0) and SILVA(128).
> On the other hand, pybarrnap default mode uses the HMM profile database created from the Rfam(14.10).
> Therefore, there will be some differences in results between Barrnap v0.9 and pybarrnap default mode.

## Installation

`Python 3.8 or later` is required for installation.
pybarrnap depends on [pyhmmer](https://github.com/althonos/pyhmmer) and [biopython](https://github.com/biopython/biopython) python library.
If accurate mode is required, please install [infernal](http://eddylab.org/infernal/) additionally.

**Install PyPI package:**

    pip install pybarrnap

**Install bioconda package:**

    conda install -c conda-forge -c bioconda pybarrnap

**Use Docker ([Image Registry](https://github.com/moshi4/pybarrnap/pkgs/container/pybarrnap)):**

    docker run -it --rm ghcr.io/moshi4/pybarrnap:latest pybarrnap -h

## CLI Usage

### Basic Command

    pybarrnap genome.fna > genome_rrna.gff

### Options

    $ pybarrnap --help
    usage: pybarrnap [options] genome.fna[.gz] > genome_rrna.gff

    Python implementation of barrnap (Bacterial ribosomal RNA predictor)

    positional arguments:
      fasta              Input fasta file (or stdin)

    optional arguments:
      -e , --evalue      E-value cutoff (default: 1e-06)
      -l , --lencutoff   Proportional length threshold to label as partial (default: 0.8)
      -r , --reject      Proportional length threshold to reject prediction (default: 0.25)
      -t , --threads     Number of threads (default: 1)
      -k , --kingdom     Target kingdom [bac|arc|euk] (default: 'bac')
      -o , --outseq      Output rRNA hit seqs as fasta file (default: None)
      -i, --incseq       Include FASTA input sequences in GFF output (default: OFF)
      -a, --accurate     Use cmscan instead of pyhmmer.nhmmer (default: OFF)
      -q, --quiet        No print log on screen (default: OFF)
      -v, --version      Print version information
      -h, --help         Show this help message and exit

> [!TIP]
> If `--accurate` option is set, cmscan(infernal) is used for rRNA search instead of pyhmmer.nhmmer.
> Although cmscan is slower than pyhmmer.nhmmer, it is expected to give more accurate results because it performs rRNA searches using RNA secondary structure profiles.

### CLI Example

Click [here](https://github.com/moshi4/pybarrnap/raw/main/examples/examples.zip) to download examples dataset.

#### CLI Example 1

Print rRNA prediction result on screen

    pybarrnap examples/bacteria.fna

#### CLI Example 2

Output rRNA predition result to file

    pybarrnap examples/archaea.fna -k arc --outseq rrna.fna --incseq > rrna_incseq.gff

#### CLI Example 3

With pipe stdin

    cat examples/fungus.fna | pybarrnap -q -k euk | grep 28S

## API Usage

pybarrnap provides simple API for running rRNA prediction and retrieving predicted rRNA.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/moshi4/pybarrnap/blob/main/notebooks/pybarrnap.ipynb)

```python
from pybarrnap import Barrnap
from pybarrnap.utils import load_example_fasta_file

# Get example fasta file path
fasta_file = load_example_fasta_file("bacteria.fna")

# Run pybarrnap rRNA prediction
barrnap = Barrnap(
    fasta_file,
    evalue=1e-6,
    lencutoff=0.8,
    reject=0.25,
    threads=1,
    kingdom="bac",
    accurate=False,
    quiet=False,
)
result = barrnap.run()

# Output rRNA GFF file
result.write_gff("bacteria_rrna.gff")
# Output rRNA GFF file (Include input fasta sequence)
result.write_gff("bacteria_rrna_incseq.gff", incseq=True)
# Output rRNA fasta file
result.write_fasta("bacteria_rrna.fna")

# Get rRNA GFF text and print
print("\n========== Print rRNA GFF ==========")
print(result.get_gff_text())

# Get rRNA features and print
print("\n========== Print rRNA features ==========")
for rec in result.seq_records:
    for feature in rec.features:
        print(feature.id, feature.type, feature.location, feature.qualifiers)

# Get rRNA sequences and print
print("\n========== Print rRNA sequences ==========")
for rec in result.get_rrna_seq_records():
    print(f">{rec.id}\n{rec.seq}")
```

## LICENSE

pybarrnap was reimplemented in python based on the perl implementation of Barrnap v0.9.
HMM(Hidden Marcov Model) and CM(Covariance Model) profile database for pybarrnap was created from Rfam(14.10).

- pybarrnap: [GPLv3](https://github.com/moshi4/pybarrnap/blob/main/LICENSE)  
- Barrnap([v0.9](https://github.com/tseemann/barrnap/tree/0.9)): [GPLv3](https://github.com/moshi4/pybarrnap/blob/main/src/pybarrnap/db/LICENSE.Barrnap)
- Rfam([14.10](https://ftp.ebi.ac.uk/pub/databases/Rfam/14.10/)): [CC0](https://github.com/moshi4/pybarrnap/blob/main/src/pybarrnap/db/LICENSE.Rfam)
