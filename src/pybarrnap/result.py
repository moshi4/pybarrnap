from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pybarrnap.record import ModelRecord


@dataclass
class BarrnapResult:
    """Barrnap Result Class"""

    mdl_records: list[ModelRecord]
    seq_records: list[SeqRecord]
    kingdom: str
    evalue: float
    lencutoff: float
    reject: float

    def __post_init__(self):
        # Sort model records (1. fasta record order, 2. rRNA feature location order)
        name2mdl_records: dict[str, list[ModelRecord]] = defaultdict(list)
        for mdl_rec in self.mdl_records:
            name2mdl_records[mdl_rec.target_name].append(mdl_rec)
        sorted_all_mdl_records: list[ModelRecord] = []
        for seq_rec in self.seq_records:
            mdl_records = name2mdl_records[seq_rec.name]
            sorted_mdl_records = sorted(mdl_records, key=lambda rec: rec.start)
            name2mdl_records[seq_rec.name] = sorted_mdl_records
            sorted_all_mdl_records.extend(sorted_mdl_records)
        self.mdl_records = sorted_all_mdl_records
        # Add features to SeqRecord
        for seq_rec in self.seq_records:
            mdl_records = name2mdl_records[seq_rec.name]
            for mdl_rec in mdl_records:
                seq_rec.features.append(mdl_rec.to_feature(self.lencutoff))

    def get_gff_text(self) -> str:
        """Get rRNA GFF text"""
        text = "##gff-version 3\n"
        for mdl_rec in self.mdl_records:
            text += mdl_rec.to_gff_line(self.lencutoff) + "\n"
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
