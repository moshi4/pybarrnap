from __future__ import annotations

import argparse
import gzip
import os
import signal
import sys
import time
from io import BufferedReader, TextIOWrapper
from pathlib import Path
from typing import Any

from Bio.SeqRecord import SeqRecord

import pybarrnap
from pybarrnap import Barrnap
from pybarrnap.config import KINGDOMS
from pybarrnap.logger import get_logger


def main():
    """Main function called from CLI"""
    args = get_args()
    run(**args.__dict__)


def run(
    fasta: str | Path | TextIOWrapper | BufferedReader | SeqRecord | list[SeqRecord],
    *,
    evalue: float = 1e-6,
    lencutoff: float = 0.8,
    reject: float = 0.25,
    threads: int = 1,
    kingdom: str = "bac",
    outseq: str | Path | None = None,
    incseq: bool = False,
    quiet: bool = False,
) -> None:
    """
    Parameters
    ----------
    fasta : str | Path | TextIOWrapper | BufferedReader | SeqRecord | list[SeqRecord]
        Fasta file (handle) or SeqRecord or list[SeqRecord]
    evalue : float, optional
        E-value cutoff
    lencutoff : float, optional
        Proportional length threshold to label as partial
    reject : float, optional
        Proportional length threshold to reject prediction
    threads : int, optional
        Number of threads
    kingdom : str, optional
        Target kingdom (`bac`|`arc`|`euk`)
    outseq : str | Path | None, optional
        Output rRNA hit seqs as fasta file
    incseq : bool, optional
        Include fasta input sequences in GFF output
    quiet : bool, optional
        If True, print log on screen
    """
    start_time = time.time()
    logger = get_logger(__name__, quiet=quiet)

    opts: dict[str, Any] = dict(
        evalue=evalue,
        lencutoff=lencutoff,
        reject=reject,
        threads=threads,
        kingdom=kingdom,
        quiet=quiet,
    )

    try:
        if isinstance(fasta, BufferedReader):
            GZIP_MAGIC_NUM = b"\x1f\x8b"
            if fasta.peek().startswith(GZIP_MAGIC_NUM):
                with gzip.open(fasta, "rb") as f:
                    result = Barrnap(TextIOWrapper(f), **opts).run()
            else:
                result = Barrnap(TextIOWrapper(fasta), **opts).run()
        else:
            result = Barrnap(fasta, **opts).run()
    except KeyboardInterrupt:
        logger.error("Interrupted")
        sys.exit(-signal.SIGINT)
    except Exception as err:
        logger.error(f"Error: {err}")
        sys.exit(getattr(err, "errno", 1))

    # Write rRNA fasta
    if outseq:
        try:
            result.write_fasta(outseq)
            logger.info(f"Write rRNA fasta file '{outseq}'")
        except OSError as err:
            logger.error(f"Failed to write rRNA fasta file '{outseq}'")
            sys.exit(err.errno)

    # Print rRNA GFF on screen
    logger.info("Sorting features and outputting rRNA GFF...")
    if incseq:
        print(result.get_gff_genome_fasta_text(), end="")
    else:
        print(result.get_gff_text(), end="")

    elapsed_time = time.time() - start_time
    logger.info(f"Done (elapsed time: {elapsed_time:.2f}[s])")


def get_args() -> argparse.Namespace:
    """Get arguments

    Returns
    -------
    args : argparse.Namespace
        Argument parameters
    """
    description = "Python implementation of barrnap (Bacterial ribosomal RNA predictor)"
    parser = argparse.ArgumentParser(
        description=description,
        usage="pybarrnap [options] genome.fna[.gz] > genome_rrna.gff",
        add_help=False,
        allow_abbrev=False,
    )

    parser.add_argument(
        "fasta",
        nargs="?",
        type=argparse.FileType("rb"),
        help="Input fasta file (or stdin)",
        default=sys.stdin,
    )
    default_evalue = 1e-6
    parser.add_argument(
        "-e",
        "--evalue",
        type=float,
        help=f"E-value cutoff (default: {default_evalue})",
        default=default_evalue,
        metavar="",
    )
    default_lencutoff = 0.8
    parser.add_argument(
        "-l",
        "--lencutoff",
        type=float,
        help="Proportional length threshold to label as partial "
        f"(default: {default_lencutoff})",
        default=default_lencutoff,
        metavar="",
    )
    default_reject = 0.25
    parser.add_argument(
        "-r",
        "--reject",
        type=float,
        help="Proportional length threshold to reject prediction "
        f"(default: {default_reject})",
        default=default_reject,
        metavar="",
    )
    default_threads = 1
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        help=f"Number of threads (default: {default_threads})",
        default=default_threads,
        metavar="",
    )
    default_kingdom = "bac"
    parser.add_argument(
        "-k",
        "--kingdom",
        type=str,
        help=f"Target kingdom [bac|arc|euk] (default: '{default_kingdom}')",
        default=default_kingdom,
        choices=KINGDOMS,
        metavar="",
    )
    default_outseq = None
    parser.add_argument(
        "-o",
        "--outseq",
        type=Path,
        help=f"Output rRNA hit seqs as fasta file (default: {default_outseq})",
        default=default_outseq,
        metavar="",
    )
    parser.add_argument(
        "-i",
        "--incseq",
        help="Include FASTA input sequences in GFF output (default: OFF)",
        action="store_true",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        help="No print log on screen (default: OFF)",
        action="store_true",
    )
    parser.add_argument(
        "-v",
        "--version",
        version=f"v{pybarrnap.__version__}",
        help="Print version information",
        action="version",
    )
    parser.add_argument(
        "-h",
        "--help",
        help="Show this help message and exit",
        action="help",
    )

    args = parser.parse_args()

    # Check evalue argument
    evalue = args.evalue
    if not evalue > 0:
        parser.error(f"--evalue must be 'value > 0' ({evalue=}).")
    # Check lencutoff argument
    lencutoff = args.lencutoff
    if not 0 <= lencutoff <= 1.0:
        parser.error(f"--lencutoff must be '0 <= value <= 1.0' ({lencutoff=}).")
    # Check reject argument
    reject = args.reject
    if not 0 <= reject <= 1.0:
        parser.error(f"--reject must be '0 <= value <= 1.0' ({reject=}).")
    # Check threads argument
    threads = args.threads
    max_threads = os.cpu_count()
    if max_threads is None:
        max_threads = 1
    if not 1 <= threads <= max_threads:
        parser.error(f"--threads must be '1 <= value <= {max_threads}' ({threads=})")

    # If input fasta is not seekable (not file) and waiting for stdin,
    # print help and exit
    fasta: BufferedReader = args.fasta
    if not fasta.seekable() and sys.stdin.isatty():
        parser.print_help()
        parser.exit(1)

    return args


if __name__ == "__main__":
    main()
