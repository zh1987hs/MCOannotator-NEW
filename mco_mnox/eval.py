from __future__ import annotations

from typing import Dict, List

import numpy as np
from sklearn.metrics import average_precision_score
from sklearn.model_selection import GroupKFold

from .data import SequenceRecord


def kmer_set(seq: str, k: int = 3):
    return {seq[i : i + k] for i in range(max(0, len(seq) - k + 1))}


def cluster_by_kmer_jaccard(records: List[SequenceRecord], threshold: float = 0.7, k: int = 3) -> np.ndarray:
    sets = [kmer_set(r.sequence, k=k) for r in records]
    groups = -np.ones(len(records), dtype=int)
    gid = 0
    for i in range(len(records)):
        if groups[i] >= 0:
            continue
        groups[i] = gid
        for j in range(i + 1, len(records)):
            if groups[j] >= 0:
                continue
            inter = len(sets[i] & sets[j])
            union = max(1, len(sets[i] | sets[j]))
            jac = inter / union
            if jac >= threshold:
                groups[j] = gid
        gid += 1
    return groups


def evaluate_pu_ranking(y_true_holdout: np.ndarray, scores: np.ndarray, top_fracs=(0.01, 0.05, 0.1)) -> Dict[str, float]:
    out: Dict[str, float] = {}
    if len(np.unique(y_true_holdout)) > 1:
        out["pr_auc"] = float(average_precision_score(y_true_holdout, scores))
    else:
        out["pr_auc"] = float("nan")
    order = np.argsort(scores)[::-1]
    for frac in top_fracs:
        k = max(1, int(len(scores) * frac))
        out[f"top_{int(frac*100)}pct_hit_rate"] = float(y_true_holdout[order[:k]].mean())
    return out


def group_kfold_indices(records: List[SequenceRecord], n_splits: int = 3, threshold: float = 0.7):
    labels = np.array([1 if r.label == 1 else 0 for r in records])
    groups = cluster_by_kmer_jaccard(records, threshold=threshold)
    gkf = GroupKFold(n_splits=n_splits)
    return list(gkf.split(np.zeros(len(records)), labels, groups=groups))
