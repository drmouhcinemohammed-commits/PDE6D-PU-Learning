# PDE6D PU-Learning Predictors

Command-line tools to predict the probability that a molecule is active
against **PDE6D**, using the two models trained under the PU-Learning
framework:

- `predict_pu.py` — MLP model (`PU_model.keras`)
- `predict_rf.py` — Random Forest model (`RandomForest_PDE6D.pkl`)

Both scripts share the same interface and the same featurization: a
Morgan fingerprint (radius = 2, 2048 bits), computed with RDKit —
identical to the one used to train the models — so outputs from the two
scripts can be merged directly (e.g. to build a consensus score).

## Installation

```bash
pip install -r requirements.txt
```

You will also need the corresponding trained model file
(`PU_model.keras` and/or `RandomForest_PDE6D.pkl`), not included in this
repository.

## Usage

### Predict a single SMILES

```bash
python predict_pu.py -i "CCOc1ccc(cc1)C(=O)O" -m PU_model.keras
python predict_rf.py -i "CCOc1ccc(cc1)C(=O)O" -m RandomForest_PDE6D.pkl
```

Output:

```
SMILES            : CCOc1ccc(cc1)C(=O)O
Canonical_SMILES  : CCOc1ccc(C(=O)O)cc1
P_active          : 0.452084
```

### Predict a CSV file of SMILES

```bash
python predict_pu.py -f molecules.csv -o predictions_mlp.csv -m PU_model.keras
python predict_rf.py -f molecules.csv -o predictions_rf.csv -m RandomForest_PDE6D.pkl
```

The input CSV must contain a column of SMILES strings (default column
name: `SMILES`; override with `--smiles-col` if your file uses a
different name, e.g. `--smiles-col canonical_smiles`).

The output CSV contains the original columns plus:

| Column            | Description                                              |
|-------------------|-----------------------------------------------------------|
| `Canonical_SMILES`| RDKit-canonicalized SMILES                                 |
| `P_active`        | Predicted probability of PDE6D activity (`NaN` if invalid) |
| `valid`           | `False` if the SMILES could not be parsed by RDKit         |

Results are sorted by `P_active` in descending order.

## Options

Both `predict_pu.py` and `predict_rf.py` accept the same flags (only the
default model filename and output filename differ):

| Flag            | Description                                              | Default (predict_pu.py) | Default (predict_rf.py)    |
|-----------------|-----------------------------------------------------------|---------------------------|------------------------------|
| `-f, --file`    | Input CSV file (mutually exclusive with `-i`)              | —                         | —                            |
| `-i, --smiles`  | A single SMILES string (mutually exclusive with `-f`)       | —                         | —                            |
| `-o, --output`  | Output CSV path (used only with `-f`)                       | `predictions.csv`        | `predictions_rf.csv`        |
| `-m, --model`   | Path to the trained model                                    | `PU_model.keras`         | `RandomForest_PDE6D.pkl`    |
| `--smiles-col`  | Name of the SMILES column in the input CSV                  | `SMILES`                  | `SMILES`                     |
| `--sep`         | CSV column separator for the input file                     | `,`                        | `,`                           |

## Notes

- Invalid or unparsable SMILES are skipped and reported with a warning;
  they receive `P_active = NaN` rather than crashing the run.
- Fingerprints are recomputed from scratch (radius = 2, 2048 bits) to
  guarantee consistency with the model's training featurization —
  do not reuse pre-computed fingerprints from another source (e.g.
  ZINC/Deep Docking's 1024-bit fingerprints), as they are not compatible.
