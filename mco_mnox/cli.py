from __future__ import annotations

import argparse
from pathlib import Path

from .model import load_model, predict, train_model
from .utils import load_config, set_seed, setup_logging


def main():
    parser = argparse.ArgumentParser("mco_mnox")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_train = sub.add_parser("train")
    p_train.add_argument("--pos", required=True)
    p_train.add_argument("--unl", required=True)
    p_train.add_argument("--neg", default=None)
    p_train.add_argument("--outdir", required=True)
    p_train.add_argument("--config", required=True)

    p_pred = sub.add_parser("predict")
    p_pred.add_argument("--model", required=True)
    p_pred.add_argument("--fasta", required=True)
    p_pred.add_argument("--out", required=True)

    args = parser.parse_args()

    if args.cmd == "train":
        config = load_config(args.config)
        set_seed(config.get("seed", 42))
        setup_logging(args.outdir)
        train_model(args.pos, args.unl, args.outdir, config, args.neg)
    elif args.cmd == "predict":
        m = load_model(args.model)
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        predict(m, args.fasta, args.out)


if __name__ == "__main__":
    main()
