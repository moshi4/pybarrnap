from __future__ import annotations

import math
import re
import shutil
import subprocess as sp
from pathlib import Path
from typing import TypeVar


def load_example_fasta_file(filename: str) -> Path:
    """Load example fasta file from local package data

    List of example fasta file name

    - `bacteria.fna`
    - `archaea.fna`
    - `fungus.fna`

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


T = TypeVar("T")


def array_split(
    ary: list[T],
    n: int,
) -> list[list[T]]:
    """Split array into multiple sub arrays

    Parameters
    ----------
    ary : list[T]
        Array
    n : int
        Split num

    Returns
    -------
    split_ary : list[list[T]]
        Split array
    """
    size = math.ceil(len(ary) / n)
    return [ary[idx : idx + size] for idx in range(0, len(ary), size)]


def is_cmscan_installed() -> bool:
    """Check cmscan is installed or not"""
    return True if shutil.which("cmscan") else False


def get_cmscan_version() -> str:
    """Get cmscan version (vX.X.X)"""
    if not is_cmscan_installed():
        raise RuntimeError("cmscan is not installed!!")
    version_unknown = "X.X.X"
    try:
        cmd_res = sp.run(["cmscan", "-h"], capture_output=True, text=True)
        if cmd_res.returncode == 0:
            pattern = r"# INFERNAL\s+(\d+\.\d+\.\d+)"
            match = re.search(pattern, cmd_res.stdout, flags=re.MULTILINE)
            version = str(match.group(1))
            return version
        else:
            return version_unknown
    except Exception:
        return version_unknown
