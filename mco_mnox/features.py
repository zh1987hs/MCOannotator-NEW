from __future__ import annotations

import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP

from .data import SequenceRecord
from .utils import seq_hash

HYDRO = set("AILMFWVYCGTP")
AA = "ACDEFGHIKLMNPQRSTVWY"

MOTIF_PATTERNS = {
    "motif1": r"H.A?HG",  # HXHG tolerant
    "motif2": r"H.A?H",  # HXH tolerant
    "motif3": r"H..H.H",  # HXXHXH
    "motif4": r"HC.H..?.?H",  # HCHXXXH tolerant
}


def _find_first(pattern: str, seq: str):
    m = re.search(pattern, seq)
    if m:
        return m.start(), m.group(0)
    return -1, ""


def motif_features(seq: str) -> Dict[str, float]:
    out: Dict[str, float] = {}
    positions = {}
    for name, pat in MOTIF_PATTERNS.items():
        pos, hit = _find_first(pat, seq)
        positions[name] = pos
        out[f"{name}_hit"] = float(pos >= 0)
        out[f"{name}_pos"] = float(pos)
        out[f"{name}_len"] = float(len(hit))
    ordered = [positions[f"motif{i}"] for i in range(1, 5)]
    for i in range(3):
        a, b = ordered[i], ordered[i + 1]
        out[f"motif_gap_{i+1}_{i+2}"] = float(b - a) if a >= 0 and b >= 0 else -1.0
    out["mco_motif_completeness"] = float(sum(out[f"motif{i}_hit"] for i in range(1, 5)) / 4.0)
    return out


def _window(seq: str, center: int, size: int = 10) -> str:
    if center < 0:
        return ""
    s = max(0, center - size)
    e = min(len(seq), center + size)
    return seq[s:e]


def t1_proxy_features(seq: str, motif_dict: Dict[str, float]) -> Dict[str, float]:
    out: Dict[str, float] = {}
    m1_pos = int(motif_dict.get("motif1_pos", -1))
    win = _window(seq, m1_pos + 2, size=12)
    if not win:
        out.update({"t1_hydrophobicity": 0.0, "t1_aromatic_ratio": 0.0, "t1_mlf_ratio": 0.0})
        return out
    out["t1_hydrophobicity"] = sum(aa in HYDRO for aa in win) / len(win)
    out["t1_aromatic_ratio"] = sum(aa in set("FWY") for aa in win) / len(win)
    out["t1_mlf_ratio"] = sum(aa in set("MLF") for aa in win) / len(win)
    return out


def signal_peptide_features(seq: str) -> Dict[str, float]:
    n = seq[:35]
    hydro_ratio = sum(aa in HYDRO for aa in n) / max(1, len(n))
    axa = 1.0 if re.search(r"A.A", n[18:35]) else 0.0
    rr = 1.0 if re.search(r"RR..FLK", n) else 0.0
    return {
        "nterm_hydro_ratio": hydro_ratio,
        "signal_peptide_flag": float(hydro_ratio > 0.42 and axa > 0),
        "tat_rr_flag": rr,
    }


def tm_helix_flag(seq: str, win: int = 19) -> float:
    for i in range(max(0, len(seq) - win + 1)):
        chunk = seq[i : i + win]
        if sum(aa in HYDRO for aa in chunk) / win > 0.68:
            return 1.0
    return 0.0


def acidic_features(seq: str) -> Dict[str, float]:
    cnt = Counter(seq)
    acidic = (cnt.get("D", 0) + cnt.get("E", 0)) / max(1, len(seq))
    acidic_cluster = 1.0 if re.search(r"[DE]{3,}", seq) else 0.0
    try:
        pI = float(IP(seq).pi())
    except Exception:
        pI = 7.0
    return {"acidic_ratio": acidic, "acidic_cluster_flag": acidic_cluster, "estimated_pI": pI}


def composition_features(seq: str) -> Dict[str, float]:
    c = Counter(seq)
    out = {
        "length": float(len(seq)),
        "his_ratio": c.get("H", 0) / max(1, len(seq)),
        "cys_ratio": c.get("C", 0) / max(1, len(seq)),
        "cys_count": float(c.get("C", 0)),
        "low_complexity_ratio": _low_complexity_ratio(seq),
        "tm_helix_flag": tm_helix_flag(seq),
    }
    return out


def _low_complexity_ratio(seq: str, k: int = 3) -> float:
    if len(seq) < k:
        return 0.0
    kmers = [seq[i : i + k] for i in range(len(seq) - k + 1)]
    uniq = len(set(kmers))
    return 1.0 - uniq / max(1, len(kmers))


def kmer_embed(seq: str, k: int = 3) -> np.ndarray:
    # deterministic fallback embedding; keeps CPU-only offline usability
    dim = len(AA) ** 2
    vec = np.zeros(dim, dtype=np.float32)
    if len(seq) < k:
        return vec
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if any(ch not in AA for ch in kmer):
            continue
        idx = (AA.index(kmer[0]) * len(AA) + AA.index(kmer[1])) % dim
        vec[idx] += 1
    vec /= max(1.0, vec.sum())
    return vec


class EmbeddingExtractor:
    def __init__(self, embedder: str = "none", cache_path: str | Path = ".cache/embeddings.json"):
        self.embedder = embedder
        self.cache_path = Path(cache_path)
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        self.cache: Dict[str, List[float]] = {}
        if self.cache_path.exists():
            self.cache = json.loads(self.cache_path.read_text(encoding="utf-8"))

    def _save(self):
        self.cache_path.write_text(json.dumps(self.cache), encoding="utf-8")

    def _embed_single(self, seq: str) -> np.ndarray:
        # Pluggable names kept for CLI compatibility; currently uses robust fallback.
        # You may extend with ESM/ProtT5 loading when weights available.
        return kmer_embed(seq, k=3)

    def embed(self, records: List[SequenceRecord]) -> np.ndarray:
        embs = []
        changed = False
        for r in records:
            h = seq_hash(r.sequence)
            if h not in self.cache:
                self.cache[h] = self._embed_single(r.sequence).tolist()
                changed = True
            embs.append(np.array(self.cache[h], dtype=np.float32))
        if changed:
            self._save()
        return np.vstack(embs)


def handcrafted_features(records: List[SequenceRecord]) -> Tuple[np.ndarray, List[Dict[str, float]], List[str]]:
    rows: List[Dict[str, float]] = []
    for r in records:
        m = motif_features(r.sequence)
        feats = {}
        feats.update(m)
        feats.update(t1_proxy_features(r.sequence, m))
        feats.update(signal_peptide_features(r.sequence))
        feats.update(acidic_features(r.sequence))
        feats.update(composition_features(r.sequence))
        rows.append(feats)
    feat_names = sorted(rows[0].keys()) if rows else []
    mat = np.array([[row[n] for n in feat_names] for row in rows], dtype=np.float32)
    return mat, rows, feat_names


def build_features(records: List[SequenceRecord], embedder: str, cache_path: str | Path):
    emb = EmbeddingExtractor(embedder=embedder, cache_path=cache_path).embed(records)
    hand, hand_rows, hand_names = handcrafted_features(records)
    x = np.hstack([emb, hand])
    feat_names = [f"emb_{i}" for i in range(emb.shape[1])] + hand_names
    return x, feat_names, hand_rows, emb
