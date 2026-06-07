"""
This script preprocesses the preveously generated hERG_ChEMBL_Mol.pkl dataset by embedding 3D molecular structures.
Relation Schema:
    R2(Smiles: str, pIC50: float, Mol: rdkit.Chem.rdchem.Mol ) ----This Script----> R3(Smiles: str, pIC50: float, Mol3D: rdkit.Chem.rdchem.Mol)

embed_3d() converts an rdkit.Chem.rdchem.Mol (mol2d) object to an rdkit.Chem.rdchem.Mol object with 3D coordinates (mol3d)
It has been implemented in ProcessMol.py
the final dataframe is stored in Data/Mol3D/hERG_ChEMBL_Mol3D.pkl as a pickle file

Here are the specifications of the env for reproducibility:
    Python 3.11.15
    Pandas 3.0.3
    Numpy 2.4.6
    Rdkit 2026.03.3
"""

from typing import cast

import pandas as pd

from ProcessMol import embed_3d

df = cast(pd.DataFrame, pd.read_pickle("Data/Mol/hERG_ChEMBL_Mol.pkl"))

print(f"Starting with {len(df)} 2D molecules...")

df["Mol3D"] = df["mol"].apply(embed_3d)

df_3d = df.dropna(subset=["Mol3D"]).copy()
df_final = df_3d[["Smiles", "pIC50", "Mol3D"]].copy()

print("\nFinal Dataset Preview:")
print(df_final.head())

df_final.to_pickle("Data/Mol3D/hERG_ChEMBL_Mol3D.pkl")

print(f"\nSuccessfully saved {len(df_final)} 3D molecules!")
