"""Extract simple structure-derived features from PDB files into TSV.

Expected input layout:
- --pdb-dir contains files named <seq_id>.pdb (or .cif not handled here)

Output columns:
seq_id, structure_mean_bfactor, structure_rg, structure_ca_count
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


def _calc_features_for_pdb(path: Path):
    from Bio.PDB import PDBParser

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(path.stem, str(path))

    cas = []
    b_factors = []
    for atom in structure.get_atoms():
        if atom.get_id() == "CA":
            coord = atom.get_coord()
            cas.append(coord)
            b_factors.append(float(atom.get_bfactor()))

    if not cas:
        return 0.0, 0.0, 0

    n = len(cas)
    cx = sum(v[0] for v in cas) / n
    cy = sum(v[1] for v in cas) / n
    cz = sum(v[2] for v in cas) / n
    rg = math.sqrt(sum((v[0] - cx) ** 2 + (v[1] - cy) ** 2 + (v[2] - cz) ** 2 for v in cas) / n)
    mean_b = sum(b_factors) / max(1, len(b_factors))
    return mean_b, rg, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    pdb_dir = Path(args.pdb_dir)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for p in sorted(pdb_dir.glob("*.pdb")):
        mean_b, rg, n_ca = _calc_features_for_pdb(p)
        rows.append((p.stem, mean_b, rg, n_ca))

    with open(out, "w", encoding="utf-8") as f:
        f.write("seq_id\tstructure_mean_bfactor\tstructure_rg\tstructure_ca_count\n")
        for sid, b, rg, n in rows:
            f.write(f"{sid}\t{b:.6f}\t{rg:.6f}\t{n}\n")

    print(f"Wrote {len(rows)} rows to {out}")


if __name__ == "__main__":
    main()
