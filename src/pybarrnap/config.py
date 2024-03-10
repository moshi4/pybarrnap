from __future__ import annotations

from pathlib import Path

# HMM/CM profile database config
_db_path = Path(__file__).parent / "db"
_hmm_db_path = _db_path / "hmm"
_cm_db_path = _db_path / "cm"

KINGDOMS = ["arc", "bac", "euk"]
KINGDOM2HMM_DB = {k: _hmm_db_path / f"{k}.hmm" for k in KINGDOMS}
KINGDOM2CM_DB = {k: _cm_db_path / f"{k}.cm" for k in KINGDOMS}

SEQTYPE2LEN: dict[str, int] = {
    "5S_rRNA": 119,
    "16S_rRNA": 1585,
    "23S_rRNA": 3232,
    "5_8S_rRNA": 156,
    "18S_rRNA": 1869,
    "28S_rRNA": 2912,
}

MAXLEN = int(max(SEQTYPE2LEN.values()) * 1.2)
