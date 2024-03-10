from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio.SeqFeature import SeqFeature, SimpleLocation
from pyhmmer.plan7 import Hit

import pybarrnap
from pybarrnap.config import SEQTYPE2LEN


@dataclass
class ModelRecord:
    """Model Record Class (HMM or CM)

    HMM: Hidden Marcov Model, CM: Covariance Model
    """

    target_name: str
    target_acc: str
    query_name: str
    query_acc: str
    mdl_from: int
    mdl_to: int
    ali_from: int
    ali_to: int
    strand: str
    evalue: float
    score: float
    bias: float
    description: str

    @property
    def start(self) -> int:
        """Start position"""
        return self.ali_from if self.strand == "+" else self.ali_to

    @property
    def end(self) -> int:
        """End position"""
        return self.ali_to if self.strand == "+" else self.ali_from

    @property
    def length(self) -> int:
        """Length"""
        return self.end - self.start + 1

    @property
    def product(self) -> str:
        """product"""
        return self.query_name.replace("_r", " ribosomal ").replace("5_8", "5.8")

    @staticmethod
    def from_hit(hit: Hit) -> ModelRecord:
        """Create a new record from a PyHMMER ``Hit``"""
        query_name = hit.hits.query_name.decode()
        query_acc = (
            "-"
            if hit.hits.query_accession is None
            else hit.hits.query_accession.decode()
        )
        dom = hit.best_domain
        ali = dom.alignment
        target_name = hit.name.decode()
        target_acc = "-" if hit.accession is None else hit.accession.decode()
        desc = "-" if hit.description is None else hit.description.decode()
        return ModelRecord(
            target_name=target_name,
            target_acc=target_acc,
            query_name=query_name,
            query_acc=query_acc,
            mdl_from=ali.hmm_from,
            mdl_to=ali.hmm_to,
            ali_from=ali.target_from,
            ali_to=ali.target_to,
            strand=dom.strand,  # type: ignore
            evalue=hit.evalue,
            score=hit.score,
            bias=dom.bias,
            description=desc,
        )

    @staticmethod
    def parse_from_cmscan_table(tbl_file: str | Path) -> list[ModelRecord]:
        """Parse from cmscan result table (format=2)"""
        mdl_records = []
        with open(tbl_file) as f:
            for line in f.read().splitlines():
                if line.startswith("#"):
                    continue
                # Order of target and query is reversed in cmscan and nhmmer
                split_line = line.split()
                mdl_record = ModelRecord(
                    target_name=split_line[3],
                    target_acc=split_line[4],
                    query_name=split_line[1],
                    query_acc=split_line[2],
                    mdl_from=int(split_line[7]),
                    mdl_to=int(split_line[8]),
                    ali_from=int(split_line[9]),
                    ali_to=int(split_line[10]),
                    strand=split_line[11],
                    evalue=float(split_line[17]),
                    score=float(split_line[16]),
                    bias=float(split_line[15]),
                    description=" ".join(split_line[26:]),
                )
                mdl_records.append(mdl_record)
        return mdl_records

    def is_partial(self, lencutoff: float = 0.8) -> bool:
        """Check partial or not"""
        return self.length / SEQTYPE2LEN[self.query_name] < lencutoff

    def to_gff_line(self, lencutoff: float = 0.8) -> str:
        """Convert to gff line

        Parameters
        ----------
        lencutoff : float, optional
            Proportional length threshold to label as partial

        Returns
        -------
        gff_line : str
            GFF line
        """
        if self.is_partial(lencutoff):
            tags = f"Name={self.query_name};product={self.product} (partial)"
            perc = self.length / SEQTYPE2LEN[self.query_name] * 100
            tags += f";note=aligned only {perc:.2f} percent of the {self.product}"
        else:
            tags = f"Name={self.query_name};product={self.product}"

        return "\t".join(
            (
                self.target_name,
                f"pybarrnap:{pybarrnap.__version__}",
                "rRNA",
                str(self.start),
                str(self.end),
                "0" if self.evalue == 0 else f"{self.evalue:.1e}",
                self.strand,
                ".",
                tags,
            )
        )

    def to_feature(self, lencutoff: float = 0.8) -> SeqFeature:
        """Convert to BioPython's SeqFeature

        Parameters
        ----------
        lencutoff : float, optional
            Proportional length threshold to label as partial

        Returns
        -------
        feature : SeqFeature
            SeqFeature object
        """
        strand = -1 if self.strand == "-" else 1
        if self.is_partial(lencutoff):
            perc = self.length / SEQTYPE2LEN[self.query_name] * 100
            qualifiers = dict(
                Name=[self.query_name],
                product=[f"{self.product} (partial)"],
                note=[f"aligned only {perc:.2f} percent of the {self.product}"],
            )
        else:
            qualifiers = dict(Name=[self.query_name], product=[self.product])

        # 1-based start is converted to 0-based
        return SeqFeature(
            location=SimpleLocation(self.start - 1, self.end, strand),
            type="rRNA",
            id=self.target_name,
            qualifiers=qualifiers,
        )

    def __repr__(self):
        return str(self)

    def __str__(self):
        perc = self.length / SEQTYPE2LEN[self.query_name] * 100
        result = ""
        result += f"{self.query_name} {self.target_name} "
        result += f"{self.start}..{self.end}({self.strand}) "
        result += f"L={self.length}/{SEQTYPE2LEN[self.query_name]}({perc:.2f}%)"
        return result
