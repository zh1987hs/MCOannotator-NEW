from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional

from Bio import SeqIO


@dataclass
class SequenceRecord:
    seq_id: str
    sequence: str
    description: str = ""
    label: Optional[int] = None  # 1 positive, 0 strong negative, None unlabeled



def read_fasta(path: str | Path, label: Optional[int] = None) -> List[SequenceRecord]:
    records: List[SequenceRecord] = []
    for rec in SeqIO.parse(str(path), "fasta"):
        records.append(
            SequenceRecord(
                seq_id=rec.id,
                sequence=str(rec.seq).upper(),
                description=rec.description,
                label=label,
            )
        )
    return records



def write_fasta(records: Iterable[SequenceRecord], path: str | Path) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for r in records:
            f.write(f">{r.seq_id} {r.description}\n{r.sequence}\n")



def load_dataset(pos_fasta: str, unl_fasta: str, neg_fasta: Optional[str] = None):
    positives = read_fasta(pos_fasta, label=1)
    unlabeled = read_fasta(unl_fasta, label=None)
    negatives = read_fasta(neg_fasta, label=0) if neg_fasta else []
    return positives, negatives, unlabeled
