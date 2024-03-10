from __future__ import annotations

from pathlib import Path

# HMM profile database config
_db_path = Path(__file__).parent / "db"

KINGDOMS = ["arc", "bac", "euk"]

SEQTYPE2LEN: dict[str, int] = {
    "5S_rRNA": 119,
    "16S_rRNA": 1585,
    "23S_rRNA": 3232,
    "5_8S_rRNA": 156,
    "18S_rRNA": 1869,
    "28S_rRNA": 2912,
}

MAXLEN = int(max(SEQTYPE2LEN.values()) * 1.2)
