from typing import cast

import pandas as pd
from numpy import log10

from ProcessMol import smiles_to_mol

# load the csv into a pandas DataFrame
df = pd.read_csv("hERG_ChEMBL.csv", sep=";", low_memory=False)

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
df_filtered["mol"] = df_filtered["Smiles"].apply(smiles_to_mol).dropna()

len2 = df_filtered.shape[0]
print(f"len1={len1}, len2={len2}, diff={len1 - len2}")
