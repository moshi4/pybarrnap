from __future__ import annotations

from pathlib import Path


def load_example_fasta_file(filename: str) -> Path:
    """Load example fasta file from local package data

    List of example fasta file name

    - `bacteria.fna`
    - `archaea.fna`
    - `fungus.fna`
    - `mitochondria.fna`
    - `mitochondria.fna.gz`

    Parameters
    ----------
    filename : str
        Fasta file name

    Returns
    -------
    fasta_file_path : Path
        Fasta file path
    """
    fasta_dir = Path(__file__).parent / "example_data"
    fasta_files = list(fasta_dir.glob("*.fna")) + list(fasta_dir.glob("*.fna.gz"))
    fasta_filenames = [f.name for f in fasta_files]

    if filename.lower() in fasta_filenames:
        return fasta_dir / filename.lower()
    else:
        err_msg = f"{filename=} is not found.\n"
        err_msg += f"Available filenames = {fasta_filenames}"
        raise FileNotFoundError(err_msg)
