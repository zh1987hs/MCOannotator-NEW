from __future__ import annotations

import csv
import json
import re
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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

ESM_MODEL_MAP = {
    "esm2_t12_35M": "facebook/esm2_t12_35M_UR50D",
    "esm2_t33_650M": "facebook/esm2_t33_650M_UR50D",
}


def _is_probable_local_path(path_str: str) -> bool:
    if re.match(r"^[A-Za-z]:[\/]", path_str):
        return True
    if path_str.startswith("\\"):
        return True
    return path_str.startswith(("./", ".\\", "../", "..\\", "/", "~"))


def _resolve_local_model_dir(local_dir: str) -> Path:
    p = Path(local_dir).expanduser()
    if not p.is_absolute():
        p = p.resolve()
    return p


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
    return {
        "length": float(len(seq)),
        "his_ratio": c.get("H", 0) / max(1, len(seq)),
        "cys_ratio": c.get("C", 0) / max(1, len(seq)),
        "cys_count": float(c.get("C", 0)),
        "low_complexity_ratio": _low_complexity_ratio(seq),
        "tm_helix_flag": tm_helix_flag(seq),
    }


def _low_complexity_ratio(seq: str, k: int = 3) -> float:
    if len(seq) < k:
        return 0.0
    kmers = [seq[i : i + k] for i in range(len(seq) - k + 1)]
    uniq = len(set(kmers))
    return 1.0 - uniq / max(1, len(kmers))


def kmer_embed(seq: str, k: int = 3) -> np.ndarray:
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


def _sequence_chunks(seq: str, chunk_size: int, overlap: int) -> List[str]:
    if len(seq) <= chunk_size:
        return [seq]
    step = max(1, chunk_size - max(0, overlap))
    chunks = []
    for i in range(0, len(seq), step):
        s = seq[i : i + chunk_size]
        if not s:
            break
        chunks.append(s)
        if i + chunk_size >= len(seq):
            break
    return chunks


def load_structure_features(
    records: List[SequenceRecord], structure_features_path: Optional[str]
) -> Tuple[np.ndarray, List[str], List[Dict[str, float]]]:
    if not structure_features_path:
        return np.zeros((len(records), 0), dtype=np.float32), [], [{} for _ in records]

    p = Path(structure_features_path)
    if not p.exists():
        raise FileNotFoundError(f"Structure feature TSV not found: {structure_features_path}")

    mapping: Dict[str, Dict[str, float]] = {}
    all_cols: List[str] = []
    with open(p, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if "seq_id" not in (reader.fieldnames or []):
            raise ValueError("Structure feature TSV must include 'seq_id' column")
        all_cols = [c for c in reader.fieldnames if c != "seq_id"]
        for row in reader:
            sid = row["seq_id"]
            mapping[sid] = {}
            for c in all_cols:
                try:
                    mapping[sid][c] = float(row[c])
                except Exception:
                    mapping[sid][c] = 0.0

    rows = []
    for r in records:
        d = dict(mapping.get(r.seq_id, {}))
        d["structure_feature_missing"] = float(r.seq_id not in mapping)
        rows.append(d)

    cols = sorted({k for d in rows for k in d.keys()})
    mat = np.array([[d.get(c, 0.0) for c in cols] for d in rows], dtype=np.float32)
    return mat, cols, rows


class EmbeddingExtractor:
    def __init__(
        self,
        embedder: str = "none",
        cache_path: str | Path = ".cache/embeddings.json",
        embedder_kwargs: Optional[Dict] = None,
    ):
        self.embedder = embedder
        self.cache_path = Path(cache_path)
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        self.cache: Dict[str, List[float]] = {}
        if self.cache_path.exists():
            self.cache = json.loads(self.cache_path.read_text(encoding="utf-8"))
        self.embedder_kwargs = embedder_kwargs or {}
        self._esm_tokenizer = None
        self._esm_model = None
        self._esm_device = None

    def _save(self):
        self.cache_path.write_text(json.dumps(self.cache), encoding="utf-8")

    def _init_esm(self) -> None:
        if self._esm_model is not None:
            return
        try:
            import torch
            from transformers import AutoModel, AutoTokenizer
        except Exception as e:
            raise ImportError(
                "ESM2 embedding requires torch + transformers. "
                "Install them and retry, or use --embedder none."
            ) from e

        local_dir = self.embedder_kwargs.get("esm_local_dir")
        force_local = bool(self.embedder_kwargs.get("esm_force_local", False))

        target: str
        local_files_only: bool
        if local_dir and str(local_dir).strip():
            local_path_raw = str(local_dir).strip()
            local_path = _resolve_local_model_dir(local_path_raw)
            if not local_path.exists() or not local_path.is_dir():
                raise FileNotFoundError(
                    "ESM local directory not found: "
                    f"{local_path_raw}. Please download model files first, or clear esm_local_dir."
                )
            target = str(local_path)
            local_files_only = True
        else:
            target = ESM_MODEL_MAP[self.embedder]
            local_files_only = force_local

        try:
            self._esm_tokenizer = AutoTokenizer.from_pretrained(target, local_files_only=local_files_only)
            self._esm_model = AutoModel.from_pretrained(target, local_files_only=local_files_only)
        except OSError as e:
            if _is_probable_local_path(str(target)):
                raise OSError(
                    f"Failed to load local ESM model from '{target}'. "
                    "Ensure config.json/tokenizer files/model weights are present in that folder."
                ) from e
            raise

        device = self.embedder_kwargs.get("esm_device", "cpu")
        if device == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"
        self._esm_device = torch.device(device)
        self._esm_model.to(self._esm_device)
        self._esm_model.eval()

    def _embed_batch_esm(self, seqs: List[str]) -> np.ndarray:
        import torch

        self._init_esm()
        max_len = int(self.embedder_kwargs.get("esm_max_length", 1024))
        encoded = self._esm_tokenizer(
            seqs,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=max_len,
        )
        encoded = {k: v.to(self._esm_device) for k, v in encoded.items()}
        with torch.no_grad():
            out = self._esm_model(**encoded)
            hidden = out.last_hidden_state
            mask = encoded["attention_mask"].unsqueeze(-1)
            pooled = (hidden * mask).sum(dim=1) / mask.sum(dim=1).clamp(min=1)
        return pooled.detach().cpu().numpy().astype(np.float32)

    def _embed_single_esm(self, seq: str) -> np.ndarray:
        max_len = int(self.embedder_kwargs.get("esm_max_length", 1024))
        strategy = str(self.embedder_kwargs.get("esm_long_strategy", "chunk_mean"))
        chunk_size = int(self.embedder_kwargs.get("esm_chunk_size", max_len))
        chunk_overlap = int(self.embedder_kwargs.get("esm_chunk_overlap", 128))

        if len(seq) <= max_len or strategy == "truncate":
            return self._embed_batch_esm([seq])[0]

        chunks = _sequence_chunks(seq, chunk_size=chunk_size, overlap=chunk_overlap)
        chunk_vecs = self._embed_batch_esm(chunks)
        weights = np.array([len(c) for c in chunks], dtype=np.float32)
        weights /= max(weights.sum(), 1.0)
        return (chunk_vecs * weights[:, None]).sum(axis=0).astype(np.float32)

    def _embed_single(self, seq: str) -> np.ndarray:
        if self.embedder == "none":
            return kmer_embed(seq, k=3)
        if self.embedder in ESM_MODEL_MAP:
            return self._embed_single_esm(seq)
        return kmer_embed(seq, k=3)

    def embed(self, records: List[SequenceRecord]) -> np.ndarray:
        embs = []
        changed = False

        missing_indices: List[int] = []
        missing_hashes: List[str] = []
        missing_seqs: List[str] = []

        for i, r in enumerate(records):
            h = seq_hash(r.sequence)
            if h in self.cache:
                embs.append(np.array(self.cache[h], dtype=np.float32))
                continue
            embs.append(None)
            missing_indices.append(i)
            missing_hashes.append(h)
            missing_seqs.append(r.sequence)

        if missing_indices:
            if self.embedder in ESM_MODEL_MAP:
                max_len = int(self.embedder_kwargs.get("esm_max_length", 1024))
                strategy = str(self.embedder_kwargs.get("esm_long_strategy", "chunk_mean"))
                short_pack = [(i, s) for i, s in enumerate(missing_seqs) if len(s) <= max_len or strategy == "truncate"]
                long_pack = [(i, s) for i, s in enumerate(missing_seqs) if len(s) > max_len and strategy != "truncate"]

                gen = [None] * len(missing_seqs)
                if short_pack:
                    bs = int(self.embedder_kwargs.get("esm_batch_size", 4))
                    idxs = [i for i, _ in short_pack]
                    seqs = [s for _, s in short_pack]
                    cursor = 0
                    for j in range(0, len(seqs), bs):
                        batch = seqs[j : j + bs]
                        vecs = self._embed_batch_esm(batch)
                        for k in range(len(batch)):
                            gen[idxs[cursor]] = vecs[k]
                            cursor += 1
                for i, s in long_pack:
                    gen[i] = self._embed_single_esm(s)
                gen_emb = np.vstack(gen)
            else:
                gen_emb = np.vstack([self._embed_single(s) for s in missing_seqs])

            for idx, h, vec in zip(missing_indices, missing_hashes, gen_emb):
                self.cache[h] = vec.tolist()
                embs[idx] = np.array(vec, dtype=np.float32)
                changed = True

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


def build_features(records: List[SequenceRecord], embedder: str, cache_path: str | Path, embedder_kwargs: Optional[Dict] = None):
    embedder_kwargs = embedder_kwargs or {}
    emb = EmbeddingExtractor(embedder=embedder, cache_path=cache_path, embedder_kwargs=embedder_kwargs).embed(records)
    hand, hand_rows, hand_names = handcrafted_features(records)

    struct_path = embedder_kwargs.get("structure_features_path", "")
    struct_mat, struct_names, struct_rows = load_structure_features(records, struct_path)
    for i in range(len(hand_rows)):
        hand_rows[i].update(struct_rows[i])

    mats = [emb, hand]
    feat_names = [f"emb_{i}" for i in range(emb.shape[1])] + hand_names
    if struct_mat.shape[1] > 0:
        mats.append(struct_mat)
        feat_names += [f"struct_{n}" for n in struct_names]

    x = np.hstack(mats)
    return x, feat_names, hand_rows, emb
