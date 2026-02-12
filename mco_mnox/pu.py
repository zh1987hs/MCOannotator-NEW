from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler


@dataclass
class PUResult:
    scores_mean: np.ndarray
    scores_std: np.ndarray
    models: List[LogisticRegression]
    scaler: StandardScaler


class BaggingPU:
    def __init__(self, n_estimators: int = 30, unlabeled_sample_ratio: float = 1.0, seed: int = 42, c: float = 1.0):
        self.n_estimators = n_estimators
        self.unlabeled_sample_ratio = unlabeled_sample_ratio
        self.seed = seed
        self.c = c

    def fit_predict(self, x_pos: np.ndarray, x_unl: np.ndarray) -> PUResult:
        rng = np.random.default_rng(self.seed)
        models = []
        preds = []
        scaler = StandardScaler()
        x_all = np.vstack([x_pos, x_unl])
        scaler.fit(x_all)
        x_pos_s = scaler.transform(x_pos)
        x_unl_s = scaler.transform(x_unl)

        n_u = len(x_unl)
        take = max(1, int(n_u * self.unlabeled_sample_ratio * min(1.0, len(x_pos) / max(1, n_u))))
        for _ in range(self.n_estimators):
            idx = rng.choice(n_u, size=take, replace=(take > n_u))
            x_neg_tmp = x_unl_s[idx]
            x_train = np.vstack([x_pos_s, x_neg_tmp])
            y_train = np.array([1] * len(x_pos_s) + [0] * len(x_neg_tmp))
            clf = LogisticRegression(C=self.c, max_iter=2000, class_weight="balanced", solver="liblinear")
            clf.fit(x_train, y_train)
            models.append(clf)
            preds.append(clf.predict_proba(x_unl_s)[:, 1])
        pred_mat = np.vstack(preds)
        return PUResult(pred_mat.mean(axis=0), pred_mat.std(axis=0), models, scaler)


class ReliableNegativePU:
    def __init__(self, rn_quantile: float = 0.15, c: float = 1.0):
        self.rn_quantile = rn_quantile
        self.c = c

    def fit_predict(self, x_pos: np.ndarray, x_unl: np.ndarray) -> Tuple[LogisticRegression, StandardScaler, np.ndarray]:
        scaler = StandardScaler()
        x_all = np.vstack([x_pos, x_unl])
        scaler.fit(x_all)
        x_pos_s = scaler.transform(x_pos)
        x_unl_s = scaler.transform(x_unl)

        weak = LogisticRegression(C=self.c, max_iter=2000, class_weight="balanced", solver="liblinear")
        x_w = np.vstack([x_pos_s, x_unl_s])
        y_w = np.array([1] * len(x_pos_s) + [0] * len(x_unl_s))
        weak.fit(x_w, y_w)
        s = weak.predict_proba(x_unl_s)[:, 1]

        thr = float(np.quantile(s, self.rn_quantile))
        rn_idx = np.where(s <= thr)[0]
        if len(rn_idx) < 3:
            rn_idx = np.argsort(s)[: min(3, len(s))]

        x_rn = x_unl_s[rn_idx]
        x_final = np.vstack([x_pos_s, x_rn])
        y_final = np.array([1] * len(x_pos_s) + [0] * len(x_rn))
        final = LogisticRegression(C=self.c, max_iter=2000, class_weight="balanced", solver="liblinear")
        final.fit(x_final, y_final)
        pred = final.predict_proba(x_unl_s)[:, 1]
        return final, scaler, pred


def explain_linear(model: LogisticRegression, feature_names: List[str], topn: int = 8) -> Dict[str, float]:
    coefs = model.coef_.ravel()
    order = np.argsort(np.abs(coefs))[::-1][:topn]
    return {feature_names[i]: float(coefs[i]) for i in order}
