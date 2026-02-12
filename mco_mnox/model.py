from __future__ import annotations

import pickle
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity

from .data import SequenceRecord, load_dataset
from .features import build_features
from .pu import BaggingPU, ReliableNegativePU, explain_linear
from .utils import save_json


@dataclass
class TrainedModel:
    config: Dict
    strategy: str
    feature_names: List[str]
    positives: List[SequenceRecord]
    pos_embeddings: np.ndarray
    handcrafted_rows: List[Dict[str, float]]
    bagging_models: Optional[List] = None
    rn_model: Optional[object] = None
    scaler: Optional[object] = None



def train_model(pos_fasta: str, unl_fasta: str, outdir: str, config: Dict, neg_fasta: Optional[str] = None) -> TrainedModel:
    outp = Path(outdir)
    outp.mkdir(parents=True, exist_ok=True)
    positives, negatives, unlabeled = load_dataset(pos_fasta, unl_fasta, neg_fasta)
    train_records = positives + negatives + unlabeled

    x_all, feat_names, hand_rows, emb = build_features(
        train_records,
        config["features"]["embedder"],
        outp / "embedding_cache.json",
        embedder_kwargs=config.get("features", {}),
    )

    n_p, n_n = len(positives), len(negatives)
    x_pos = x_all[:n_p]
    x_neg = x_all[n_p : n_p + n_n] if n_n > 0 else np.zeros((0, x_all.shape[1]))
    x_unl = x_all[n_p + n_n :]

    strategy = config["pu"]["strategy"]
    if strategy == "bagging":
        pu = BaggingPU(
            n_estimators=config["pu"]["bagging"]["n_estimators"],
            unlabeled_sample_ratio=config["pu"]["bagging"]["unlabeled_sample_ratio"],
            seed=config["seed"],
            c=config["model"]["logreg_c"],
        )
        if len(x_neg) > 0:
            x_unl = np.vstack([x_unl, x_neg])
        res = pu.fit_predict(x_pos, x_unl)
        model = TrainedModel(
            config=config,
            strategy=strategy,
            feature_names=feat_names,
            positives=positives,
            pos_embeddings=emb[:n_p],
            handcrafted_rows=hand_rows[:n_p],
            bagging_models=res.models,
            scaler=res.scaler,
        )
    elif strategy == "rn":
        rn = ReliableNegativePU(rn_quantile=config["pu"]["rn"]["rn_quantile"], c=config["model"]["logreg_c"])
        if len(x_neg) > 0:
            x_unl = np.vstack([x_unl, x_neg])
        m, s, _ = rn.fit_predict(x_pos, x_unl)
        model = TrainedModel(
            config=config,
            strategy=strategy,
            feature_names=feat_names,
            positives=positives,
            pos_embeddings=emb[:n_p],
            handcrafted_rows=hand_rows[:n_p],
            rn_model=m,
            scaler=s,
        )
    else:
        raise ValueError(f"Unknown PU strategy: {strategy}")

    with open(outp / "model.pkl", "wb") as f:
        pickle.dump(model, f)
    save_json(config, outp / "config.used.json")
    return model



def load_model(path: str) -> TrainedModel:
    with open(path, "rb") as f:
        return pickle.load(f)



def _predict_scores(model: TrainedModel, x: np.ndarray):
    x_s = model.scaler.transform(x)
    if model.strategy == "bagging":
        pred = np.vstack([m.predict_proba(x_s)[:, 1] for m in model.bagging_models])
        return pred.mean(axis=0), pred.std(axis=0), model.bagging_models[0]
    p = model.rn_model.predict_proba(x_s)[:, 1]
    return p, np.zeros_like(p), model.rn_model



def predict(model: TrainedModel, fasta: str, out_tsv: str):
    from .data import read_fasta
    records = read_fasta(fasta, label=None)
    x, feat_names, hand_rows, emb = build_features(
        records,
        model.config["features"]["embedder"],
        Path(out_tsv).parent / "embedding_cache.predict.json",
        embedder_kwargs=model.config.get("features", {}),
    )
    scores, stds, expl_model = _predict_scores(model, x)
    pos_sim = cosine_similarity(emb, model.pos_embeddings)
    near_idx = np.argmax(pos_sim, axis=1)

    explanations = explain_linear(expl_model, feat_names, topn=model.config["predict"]["top_explain_features"])

    df = pd.DataFrame(
        {
            "seq_id": [r.seq_id for r in records],
            "length": [len(r.sequence) for r in records],
            "score_mean": scores,
            "score_std": stds,
            "nearest_positive": [model.positives[i].seq_id for i in near_idx],
            "nearest_positive_cosine": [float(pos_sim[j, i]) for j, i in enumerate(near_idx)],
        }
    )
    for k in [
        "signal_peptide_flag",
        "tat_rr_flag",
        "tm_helix_flag",
        "acidic_ratio",
        "estimated_pI",
        "mco_motif_completeness",
        "motif1_hit",
        "motif2_hit",
        "motif3_hit",
        "motif4_hit",
    ]:
        df[k] = [row.get(k, np.nan) for row in hand_rows]
    df["noncanonical_mco_flag"] = (df["mco_motif_completeness"] < model.config["predict"]["min_motif_completeness"]).astype(int)
    df["score_adjusted"] = df["score_mean"] * np.where(df["noncanonical_mco_flag"] == 1, 0.7, 1.0)
    df = df.sort_values("score_adjusted", ascending=False).reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)
    df["explanation"] = str(explanations)
    df.to_csv(out_tsv, sep="\t", index=False)

    out_json = {
        "high_score_low_uncertainty": df[(df.score_mean > 0.7) & (df.score_std < 0.1)]["seq_id"].tolist(),
        "high_score_high_uncertainty": df[(df.score_mean > 0.7) & (df.score_std >= 0.1)]["seq_id"].tolist(),
        "boundary_with_key_motifs": df[
            (df.score_mean.between(0.4, 0.7)) & (df.mco_motif_completeness >= 0.75)
        ]["seq_id"].tolist(),
    }
    save_json(out_json, Path(out_tsv).with_suffix(".json"))
    return df, out_json
