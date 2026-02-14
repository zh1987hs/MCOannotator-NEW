"""Download ESM model files with huggingface_hub API (Windows-friendly)."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Download ESM2 model snapshot to local dir.")
    parser.add_argument("--repo-id", required=True, help="e.g. facebook/esm2_t33_650M_UR50D")
    parser.add_argument("--local-dir", required=True, help="output directory")
    parser.add_argument("--token", default=None, help="HF token if needed")
    parser.add_argument(
        "--allow-pattern",
        action="append",
        default=None,
        help="optional allow pattern; can pass multiple times",
    )
    args = parser.parse_args()

    try:
        from huggingface_hub import snapshot_download
    except Exception as e:
        raise SystemExit(
            "huggingface_hub is required. Install with: pip install huggingface_hub"
        ) from e

    local_dir = Path(args.local_dir)
    local_dir.mkdir(parents=True, exist_ok=True)

    snapshot_download(
        repo_id=args.repo_id,
        local_dir=str(local_dir),
        token=args.token,
        resume_download=True,
        local_dir_use_symlinks=False,
        allow_patterns=args.allow_pattern,
    )
    print(f"Downloaded {args.repo_id} to {local_dir}")


if __name__ == "__main__":
    main()
