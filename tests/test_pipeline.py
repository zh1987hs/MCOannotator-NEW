from pathlib import Path

from mco_mnox.model import load_model, predict, train_model
from mco_mnox.utils import load_config


def test_train_predict_toy(tmp_path: Path):
    cfg = load_config("configs/default.yaml")
    cfg["pu"]["bagging"]["n_estimators"] = 5
    outdir = tmp_path / "run"
    train_model(
        "data/toy/positives.fasta",
        "data/toy/unlabeled.fasta",
        str(outdir),
        cfg,
        "data/toy/negatives.fasta",
    )
    model = load_model(str(outdir / "model.pkl"))
    out_tsv = tmp_path / "results.tsv"
    df, summary = predict(model, "data/toy/unlabeled.fasta", str(out_tsv))
    assert out_tsv.exists()
    assert len(df) == 5
    assert "high_score_low_uncertainty" in summary
    assert "score_mean" in df.columns
