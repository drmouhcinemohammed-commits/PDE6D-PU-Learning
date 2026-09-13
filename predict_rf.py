#!/usr/bin/env python3
# ============================================================
# predict_rf.py
#
# Command-line tool to predict the probability that a molecule
# is active against PDE6D, using the trained Random Forest model
# (RandomForest_PDE6D.pkl).
#
# Usage
# -----
# Predict for a single SMILES:
#     python predict_rf.py -i "CCOc1ccc(cc1)C(=O)O" -m RandomForest_PDE6D.pkl
#
# Predict for a CSV file of SMILES:
#     python predict_rf.py -f molecules.csv -o predictions.csv -m RandomForest_PDE6D.pkl
#
# The input CSV must contain a column with SMILES strings
# (default column name: "SMILES", override with --smiles-col).
#
# The output CSV contains the original columns plus:
#     - Canonical_SMILES : RDKit-canonicalized SMILES
#     - P_active         : predicted probability of PDE6D activity
#     - valid            : False if the SMILES could not be parsed
# ============================================================

import argparse
import sys
import numpy as np
import pandas as pd
import joblib

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

# Silence RDKit's warnings about invalid SMILES; we handle them ourselves.
RDLogger.DisableLog("rdApp.*")

RADIUS = 2
N_BITS = 2048


def compute_fingerprint(smiles: str):
    """Return a (2048,) Morgan fingerprint array, or None if the SMILES is invalid."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    canonical = Chem.MolToSmiles(mol)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, RADIUS, nBits=N_BITS)
    arr = np.zeros((N_BITS,), dtype=np.int8)
    Chem.DataStructs.ConvertToNumpyArray(fp, arr)
    return canonical, arr


def load_model(model_path: str):
    """Load the trained Random Forest model (joblib/pickle)."""
    return joblib.load(model_path)


def predict_dataframe(df: pd.DataFrame, smiles_col: str, model) -> pd.DataFrame:
    canonical_list = []
    valid_list = []
    fps = []

    for smi in df[smiles_col].astype(str):
        canonical, fp = compute_fingerprint(smi)
        if fp is None:
            canonical_list.append(None)
            valid_list.append(False)
            fps.append(np.zeros((N_BITS,), dtype=np.int8))
        else:
            canonical_list.append(canonical)
            valid_list.append(True)
            fps.append(fp)

    X = np.vstack(fps)
    probabilities = model.predict_proba(X)[:, 1]

    out = df.copy()
    out["Canonical_SMILES"] = canonical_list
    out["P_active"] = probabilities
    out["valid"] = valid_list

    # Invalid SMILES get no meaningful prediction.
    out.loc[~out["valid"], "P_active"] = np.nan

    n_invalid = (~out["valid"]).sum()
    if n_invalid:
        print(f"Warning: {n_invalid} SMILES could not be parsed and were skipped "
              f"(P_active = NaN).", file=sys.stderr)

    return out


def main():
    parser = argparse.ArgumentParser(
        description="Predict PDE6D activity probability from SMILES using the "
                    "trained Random Forest model (RandomForest_PDE6D.pkl)."
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "-f", "--file",
        help="Input CSV file containing a column of SMILES strings."
    )
    input_group.add_argument(
        "-i", "--smiles",
        help="A single SMILES string to predict."
    )
    parser.add_argument(
        "-o", "--output",
        default="predictions_rf.csv",
        help="Output CSV file (used only with -f). Default: predictions_rf.csv"
    )
    parser.add_argument(
        "-m", "--model",
        default="RandomForest_PDE6D.pkl",
        help="Path to the trained Random Forest model (.pkl). "
             "Default: RandomForest_PDE6D.pkl"
    )
    parser.add_argument(
        "--smiles-col",
        default="SMILES",
        help="Name of the SMILES column in the input CSV. Default: SMILES"
    )
    parser.add_argument(
        "--sep",
        default=",",
        help="CSV column separator for the input file. Default: ','"
    )
    args = parser.parse_args()

    print(f"Loading model: {args.model}")
    model = load_model(args.model)

    if args.smiles:
        canonical, fp = compute_fingerprint(args.smiles)
        if fp is None:
            print(f"Error: could not parse SMILES: {args.smiles}", file=sys.stderr)
            sys.exit(1)
        proba = model.predict_proba(fp.reshape(1, -1))[:, 1][0]
        print(f"SMILES            : {args.smiles}")
        print(f"Canonical_SMILES  : {canonical}")
        print(f"P_active          : {proba:.6g}")
        return

    print(f"Reading input file: {args.file}")
    df = pd.read_csv(args.file, sep=args.sep)

    if args.smiles_col not in df.columns:
        print(f"Error: column '{args.smiles_col}' not found in {args.file}. "
              f"Available columns: {list(df.columns)}", file=sys.stderr)
        sys.exit(1)

    print(f"Predicting {len(df)} molecules...")
    result = predict_dataframe(df, args.smiles_col, model)

    result = result.sort_values(by="P_active", ascending=False)
    result.to_csv(args.output, index=False)

    n_valid = result["valid"].sum()
    print(f"Done. {n_valid}/{len(result)} molecules successfully predicted.")
    print(f"Results saved to: {args.output}")


if __name__ == "__main__":
    main()
