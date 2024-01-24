from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pybarrnap.record import HmmRecord


@dataclass
class BarrnapResult:
    """Barrnap Result Class"""

    hmm_records: list[HmmRecord]
    seq_records: list[SeqRecord]
    kingdom: str
    evalue: float
    lencutoff: float
    reject: float

    def __post_init__(self):
        # Sort hmm records (1. fasta record order, 2. rRNA feature location order)
        name2hmm_records: dict[str, list[HmmRecord]] = defaultdict(list)
        for hmm_rec in self.hmm_records:
            name2hmm_records[hmm_rec.target_name].append(hmm_rec)
        sorted_all_hmm_records: list[HmmRecord] = []
        for seq_rec in self.seq_records:
            hmm_records = name2hmm_records[seq_rec.name]
            sorted_hmm_records = sorted(hmm_records, key=lambda rec: rec.start)
            name2hmm_records[seq_rec.name] = sorted_hmm_records
            sorted_all_hmm_records.extend(sorted_hmm_records)
        self.hmm_records = sorted_all_hmm_records
        # Add features to SeqRecord
        for seq_rec in self.seq_records:
            hmm_records = name2hmm_records[seq_rec.name]
            for hmm_rec in hmm_records:
                seq_rec.features.append(hmm_rec.to_feature(self.lencutoff))

    def get_gff_text(self) -> str:
        """Get rRNA GFF text"""
        text = "##gff-version 3\n"
        for hmm_rec in self.hmm_records:
            text += hmm_rec.to_gff_line(self.lencutoff) + "\n"
        return text

    def get_gff_genome_fasta_text(self) -> str:
        """Get rRNA GFF + genome FASTA text"""
        text = self.get_gff_text()
        text += "##FASTA\n"
        for rec in self.seq_records:
            seq = str(rec.seq)
            wrap_seq = "\n".join([seq[x : x + 70] for x in range(0, len(seq), 70)])
            text += f">{rec.description}\n{wrap_seq}\n"
        return text

    def get_rrna_seq_records(self) -> list[SeqRecord]:
        """Get rRNA SeqRecord list"""
        rrna_seq_records = []
        for seq_rec in self.seq_records:
            for feature in seq_rec.features:
                start = int(feature.location.start)  # type: ignore
                end = int(feature.location.end)  # type: ignore
                strand = "-" if feature.location.strand == -1 else "+"
                seq = str(feature.extract(str(seq_rec.seq)))
                name = str(feature.qualifiers.get("Name", [None])[0])
                desc = f"{name}::{seq_rec.name}:{start}-{end}({strand})"
                rrna_seq_records.append(SeqRecord(Seq(seq), id=desc, description=desc))
        return rrna_seq_records

    def write_gff(self, outfile: str | Path, *, incseq: bool = False) -> None:
        """Write rRNA GFF file

        Parameters
        ----------
        outfile : str | Path
            Output file path
        incseq : bool, optional
            Include fasta input sequences in GFF output
        """
        with open(outfile, "w") as f:
            if incseq:
                f.write(self.get_gff_genome_fasta_text())
            else:
                f.write(self.get_gff_text())

    def write_fasta(self, outfile: str | Path) -> None:
        """Write rRNA fasta file

        Parameters
        ----------
        outfile : str | Path
            Output file path
        """
        rrna_seq_records = self.get_rrna_seq_records()
        SeqIO.write(rrna_seq_records, handle=outfile, format="fasta-2line")
