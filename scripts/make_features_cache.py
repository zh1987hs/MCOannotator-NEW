import argparse
from pathlib import Path

from mco_mnox.data import read_fasta
from mco_mnox.features import build_features


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--embedder", default="none")
    ap.add_argument("--cache", default=".cache/embeddings.json")
    args = ap.parse_args()

    records = read_fasta(args.fasta, label=None)
    build_features(records, embedder=args.embedder, cache_path=Path(args.cache))
    print(f"Cached features for {len(records)} sequences to {args.cache}")


if __name__ == "__main__":
    main()
