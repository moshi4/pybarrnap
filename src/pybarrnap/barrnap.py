from __future__ import annotations

import gzip
import io
import platform
import sys
from copy import deepcopy
from pathlib import Path

import Bio
import pyhmmer
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pyhmmer import nhmmer
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence
from pyhmmer.plan7 import HMMFile

import pybarrnap
from pybarrnap.config import KINGDOM2HMM_FILE, KINGDOMS, MAXLEN, SEQTYPE2LEN
from pybarrnap.logger import get_logger
from pybarrnap.record import HmmRecord
from pybarrnap.result import BarrnapResult


class Barrnap:
    """Barrnap rRNA Prediction Class"""

    def __init__(
        self,
        fasta: str | Path | io.TextIOWrapper | SeqRecord | list[SeqRecord],
        *,
        evalue: float = 1e-6,
        lencutoff: float = 0.8,
        reject: float = 0.25,
        threads: int = 1,
        kingdom: str = "bac",
        quiet: bool = False,
    ) -> None:
        """
        Parameters
        ----------
        fasta : str | Path | TextIOWrapper | SeqRecord | list[SeqRecord]
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
            Target kingdom (`bac`|`arc`|`euk`|`mito`)
        quiet : bool, optional
            If True, print log on screen
        """
        # Load fasta as SeqRecord list
        if isinstance(fasta, (str, Path)):
            if Path(fasta).suffix == ".gz":
                with gzip.open(fasta, "rt") as f:
                    seq_records: list[SeqRecord] = list(SeqIO.parse(f, format="fasta"))
            else:
                seq_records: list[SeqRecord] = list(SeqIO.parse(fasta, format="fasta"))
        elif isinstance(fasta, io.TextIOWrapper):
            seq_records: list[SeqRecord] = list(SeqIO.parse(fasta, format="fasta"))
        elif isinstance(fasta, SeqRecord):
            seq_records = [fasta]
        else:
            seq_records = fasta
        if len(seq_records) == 0:
            raise ValueError("No sequence found in input records")
        fasta_name_list = [str(rec.name) for rec in seq_records]
        if len(set(fasta_name_list)) != len(fasta_name_list):
            raise ValueError("Duplicate names in input records")
        self._seq_records = seq_records

        # Convert SeqRecord to DigitalSequenceBlock for pyhmmer.nhmmer execution
        try:
            alphabet = Alphabet.rna()
            self._seqs = DigitalSequenceBlock(alphabet)
            for rec in seq_records:
                name, description = rec.name.encode(), rec.description.encode()
                seq = TextSequence(name, description, sequence=str(rec.seq))
                self._seqs.append(seq.digitize(alphabet))
        except ValueError as e:
            raise ValueError(
                "Failed to convert nucleotide sequences, maybe input contains proteins?"
            ) from e

        # Set parameters
        self._evalue = evalue
        self._lencutoff = lencutoff
        self._reject = reject
        self._threads = threads
        self._quiet = quiet

        if kingdom not in KINGDOMS:
            raise ValueError(f"{kingdom=} is invalid ({KINGDOMS}).")
        self._kingdom = kingdom
        self._hmm_file = KINGDOM2HMM_FILE[kingdom]

        with HMMFile(self._hmm_file) as hf:
            self._hmms = list(hf)

    def run(self) -> BarrnapResult:
        """Run rRNA prediction

        Returns
        -------
        result : BarrnapResult
            Barrnap result
        """
        logger = get_logger(__name__, quiet=self._quiet)
        logger.info(f"Run pybarrnap v{pybarrnap.__version__}")
        logger.info(f"Operating System: {sys.platform}")
        logger.info(f"Python Version: v{platform.python_version()}")
        logger.info(f"Check Dependencies: pyhmmer v{pyhmmer.__version__} is installed")
        logger.info(f"Check Dependencies: biopython v{Bio.__version__} is installed")
        logger.info(f"Set Option: evalue={self._evalue}")
        logger.info(f"Set Option: lencutoff={self._lencutoff}")
        logger.info(f"Set Option: reject={self._reject}")
        logger.info(f"Set Option: threads={self._threads}")
        logger.info(f"Set Option: kingdom='{self._kingdom}'")
        logger.info(f"Number of Target Sequence = {len(self._seq_records)}")
        for idx, rec in enumerate(self._seq_records, 1):
            name, length, description = rec.name, len(str(rec.seq)), rec.description
            logger.info(f"Seq{idx}. {name=}, {length=:,}, {description=}")
        logger.info(f"Use HMM DB: {self._hmm_file}")

        logger.info("Run pyhmmer.nhmmer")
        all_hmm_records: list[HmmRecord] = []
        opts = dict(cpus=self._threads, E=self._evalue, window_length=MAXLEN)
        for hits in nhmmer(self._hmms, self._seqs, **opts):  # type: ignore
            # Extract results and filter by length threshold
            for hit in hits.reported:
                rec = HmmRecord.from_hit(hit)
                if rec.length < int(SEQTYPE2LEN[rec.query_name] * self._reject):
                    logger.info(f"Reject: {rec}")
                    continue
                logger.info(f"Found: {rec}")
                all_hmm_records.append(rec)
        logger.info(f"Found {len(all_hmm_records)} ribosomal RNA features")

        return BarrnapResult(
            all_hmm_records,
            deepcopy(self._seq_records),
            self._kingdom,
            self._evalue,
            self._lencutoff,
            self._reject,
        )
