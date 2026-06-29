import numpy as np
import pandas as pd
from ProcessMol import compute_descriptors

df = pd.read_pickle("../Data/Mol3D/hERG_ChEMBL_Mol3D.pkl")

descriptors_series = df["Mol3D"].apply(compute_descriptors)

feature_names = [
    "MolWt",
    "LogP",
    "TPSA",
    "HBD",
    "HBA",
    "RotB",
    "FormalCharge",
    "RadiusGyr",
    "Asphericity",
    "Spherocity",
    "PMI1",
    "PMI2",
    "PMI3",
]

bit_names = [f"Morgan_Bit_{i}" for i in range(2048)]
all_column_names = feature_names + bit_names

valid_mask = descriptors_series.notnull()
valid_series = descriptors_series[valid_mask]

desc_df = pd.DataFrame(
    valid_series.tolist(),  # .tolist() is the fastest way to expand arrays in pandas
    columns=all_column_names,
    index=valid_series.index,  # Preserve the original index to merge accurately
)

df_final = df[valid_mask].join(desc_df)

print(f"Original shape: {df.shape}")
print(f"Final shape with descriptors: {df_final.shape}")

df_final.to_pickle("../Data/Mol3D/hERG_ChEMBL_Mol3D_descriptors.pkl")

print(df_final.head())
