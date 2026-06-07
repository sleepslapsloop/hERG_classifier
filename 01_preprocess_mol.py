"""
This script processes the hERG ChEMBL dataset stored in Data/hERG_ChEMBL.csv.
Relation Schema:
    R1(Chembl attr) -----This Script----> R2(Smiles: str, Standard Type: str, Standard Relation: str, Standard Value: float, Standard Units: str, mol: rdkit.Chem.rdchem.Mol)
    (Only Smiles, pIC50, mol are saved)
Standard Type refers to the type of standard used for the assay (e.g., IC50, Ki, Kd), here it is IC50
Standard Relation refers to the relation between the standard value and the target (e.g., '='), here '=' has been selected
Standard Units refers to the units of the standard value (e.g., nM, μM, mM), here it is nM

smiles_to_mol() converts a SMILES string to an rdkit.Chem.rdchem.Mol object
It has been implemented in ProcessMol.py
the final dataframe is stored in Data/Mol/hERG_ChEMBL_mol.pkl as a pickle file

Here are the specifications of the env for reproducibility:
    Python 3.11.15
    Pandas 3.0.3
    Numpy 2.4.6
    Rdkit 2026.03.3

In case the pickle file cannot be loaded, Data/Mol/hERG_ChEMBL_Mol.csv can be used instead.
While it does not contain the mol column, it contains pIC50 values and can be used to load the SMILES strings and convert them to mol objects using smiles_to_mol().
"""

from typing import cast

import pandas as pd
from numpy import log10

from ProcessMol import smiles_to_mol

# load the csv into a pandas DataFrame
df = pd.read_csv("Data/hERG_ChEMBL.csv", sep=";", low_memory=False)

# keep only necessary columns
df_col: pd.DataFrame = df.loc[
    :,
    [
        "Smiles",
        "Standard Type",
        "Standard Relation",
        "Standard Value",
        "Standard Units",
    ],
].copy()

df_filtered: pd.DataFrame = cast(
    pd.DataFrame,
    df_col[
        (df_col["Standard Type"] == "IC50")
        & (df_col["Standard Relation"] == "'='")
        & (df_col["Standard Units"] == "nM")
    ].copy(),
)


len1 = df_filtered.shape[0]

df_filtered["pIC50"] = 9 - log10(df_filtered["Standard Value"])
df_filtered["mol"] = df_filtered["Smiles"].apply(smiles_to_mol)
df_mol = df_filtered.dropna(subset=["mol"]).copy()

len2 = df_mol.shape[0]

df_mol = df_mol.groupby("Smiles", as_index=False).agg({"pIC50": "mean", "mol": "first"})

len3 = df_mol.shape[0]

print(f"Number of rows before mol generation = {len1}")
print(f"Number of rows after dropping failures = {len2}")
print(f"Duplicates squashed = {len2 - len3}")
print(f"Final unique molecules for model = {len3}")

# save the mol column to a csv and pickle file iff all rows have been converted
# this has been commented out because the files have already been saved
# if len1 == len2:
#     df_mol.to_pickle("Data/Mol/hERG_ChEMBL_Mol.pkl")
#     df_mol.drop(columns=["mol"]).to_csv("Data/Mol/hERG_ChEMBL_Mol.csv", index=False)
# else:
#     print(f"Warning: {len1 - len2} rows could not be converted to mol objects")
#     print("Enter y to save anyway, or n to exit")
#     choice = input()
#     if choice == "y":
#         df_mol.to_pickle("Data/Mol/hERG_ChEMBL_Mol.pkl")
#         df_mol.drop(columns=["mol"]).to_csv("Data/Mol/hERG_ChEMBL_Mol.csv", index=False)
#     else:
#         exit()
