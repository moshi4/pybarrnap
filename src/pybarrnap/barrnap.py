from __future__ import annotations

import gzip
import io
import logging
import platform
import shlex
import subprocess as sp
import sys
import textwrap
from pathlib import Path
from tempfile import TemporaryDirectory

import Bio
import pyhmmer
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pyhmmer import nhmmer
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence
from pyhmmer.plan7 import HMMFile

import pybarrnap
from pybarrnap.config import (
    KINGDOM2CM_DB,
    KINGDOM2HMM_DB,
    KINGDOMS,
    MAXLEN,
    SEQTYPE2LEN,
)
from pybarrnap.logger import get_logger
from pybarrnap.record import ModelRecord
from pybarrnap.result import BarrnapResult
from pybarrnap.utils import get_cmscan_version, is_cmscan_installed


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
        accurate: bool = False,
        quiet: bool = True,
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
            Target kingdom (`bac`|`arc`|`euk`|`all`)
            kingdom=`all` is available only when set with `accurate=True`
        accurate : bool, optional
            If True, use cmscan(infernal) instead of pyhmmer.nhmmer.
            cmscan installation is required to enable this option.
        quiet : bool, optional
            If True, no print log on screen
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

        # Set parameters
        self._evalue = evalue
        self._lencutoff = lencutoff
        self._reject = reject
        self._threads = threads
        self._accurate = accurate
        self._quiet = quiet

        if kingdom not in KINGDOMS:
            raise ValueError(f"{kingdom=} is invalid ({KINGDOMS}).")
        if kingdom == "all" and not accurate:
            raise ValueError("kingdom='all' must be set with accurate=True")
        self._kingdom = kingdom

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
        if self._accurate:
            if is_cmscan_installed():
                version = f"v{get_cmscan_version()}"
                logger.info(f"Check Dependencies: cmscan {version} is installed")
            else:
                err_msg = textwrap.dedent(
                    """
                    Check Dependencies: cmscan is not installed!!

                    Please install cmscan(infernal) to enable pybarrnap accurate mode.

                    # Install bioconda package
                    $ conda install -c bioconda infernal

                    # Install ubuntu(debian) package
                    $ sudo apt-get install -y infernal

                    # Install homebrew science package
                    $ brew tap brewsci/bio
                    $ brew install infernal

                    """
                )
                for err_msg_line in err_msg.splitlines():
                    logger.error(err_msg_line)
                raise RuntimeError("Failed to run pybarrnap accurate mode.")
        logger.info(f"Set Option: evalue={self._evalue}")
        logger.info(f"Set Option: lencutoff={self._lencutoff}")
        logger.info(f"Set Option: reject={self._reject}")
        logger.info(f"Set Option: threads={self._threads}")
        logger.info(f"Set Option: kingdom='{self._kingdom}'")
        logger.info(f"Set Option: accurate={self._accurate}")
        logger.info(f"Number of Target Sequence = {len(self._seq_records)}")
        for idx, rec in enumerate(self._seq_records, 1):
            name, length, description = rec.name, len(str(rec.seq)), rec.description
            logger.info(f"Seq{idx}. {name=}, {length=:,}, {description=}")

        if self._accurate:
            mdl_records = self._run_cmscan(logger)
        else:
            mdl_records = self._run_nhmmer(logger)

        filtered_mdl_records = []
        for rec in mdl_records:
            if rec.length / int(SEQTYPE2LEN[rec.query_name]) < self._reject:
                logger.info(f"Reject: {rec}")
            else:
                logger.info(f"Found: {rec}")
                filtered_mdl_records.append(rec)
        logger.info(f"Found {len(filtered_mdl_records)} ribosomal RNA features")

        return BarrnapResult(
            filtered_mdl_records,
            self._seq_records,
            self._kingdom,
            self._evalue,
            self._lencutoff,
            self._reject,
        )

    def _run_nhmmer(self, logger: logging.Logger) -> list[ModelRecord]:
        """Run pyhmmer.nhmmer (faster, lower accuracy for rRNA prediction)"""
        # Setup HMM database
        hmm_db = KINGDOM2HMM_DB[self._kingdom]
        with HMMFile(hmm_db) as hf:
            hmms = list(hf)
        logger.info(f"Use HMM DB: {hmm_db}")

        try:
            # Convert SeqRecord to DigitalSequenceBlock
            alphabet = Alphabet.rna()
            seq_block = DigitalSequenceBlock(alphabet)
            for rec in self._seq_records:
                name, description = rec.name.encode(), rec.description.encode()
                seq = TextSequence(name, description, sequence=str(rec.seq))
                seq_block.append(seq.digitize(alphabet))
        except ValueError as e:
            raise ValueError(
                "Failed to convert nucleotide sequences, maybe input contains proteins?"
            ) from e

        # Run pyhmmer.nhmmer
        logger.info("Running pyhmmer.nhmmer...")
        mdl_records: list[ModelRecord] = []
        opts = dict(cpus=self._threads, E=self._evalue, window_length=MAXLEN)
        for hits in nhmmer(hmms, seq_block, **opts):  # type: ignore
            for hit in hits.reported:
                mdl_records.append(ModelRecord.from_hit(hit))
        return mdl_records

    def _run_cmscan(self, logger: logging.Logger) -> list[ModelRecord]:
        """Run cmscan (slower, higher accuracy for rRNA prediction)"""
        # Setup CM database
        cm_db = KINGDOM2CM_DB[self._kingdom]
        logger.info(f"Use CM DB: {cm_db}")

        # Run cmscan
        logger.info("Running cmscan...")
        with TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            # Write tmp seq file
            seq_fasta_file = tmpdir / "seq.fna"
            SeqIO.write(self._seq_records, seq_fasta_file, format="fasta")
            # Run cmscan
            result_file = tmpdir / "result.tblout"
            total_seq_len = sum([len(rec.seq) for rec in self._seq_records])
            Z = 2 * total_seq_len / 1000000
            cmd = f"cmscan --rfam --nohmmonly --noali --cut_ga --oskip --fmt 2 --cpu {self._threads} -Z {Z} --tblout {result_file} {cm_db} {seq_fasta_file}"  # noqa: E501
            logger.info(f"$ {cmd}")
            cmd_args = shlex.split(cmd)
            cmd_res = sp.run(cmd_args, capture_output=True, text=True)

            if cmd_res.returncode == 0:
                return ModelRecord.parse_from_cmscan_table(result_file, self._evalue)
            else:
                logger.error("Failed to run command below!!")
                logger.error(f"$ {cmd}")
                stdout_lines = cmd_res.stdout.splitlines()
                if len(stdout_lines) > 0:
                    logger.error("STDOUT:")
                    for line in stdout_lines:
                        logger.error(line)
                stderr_lines = cmd_res.stderr.splitlines()
                if len(stderr_lines) > 0:
                    logger.error("STDERR:")
                    for line in stderr_lines:
                        logger.error(line)
                raise RuntimeError("Failed to run cmscan.")
