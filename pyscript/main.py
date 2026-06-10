import pandas as pd

df = pd.read_pickle("Data/Docked/hERG_ChEMBL_Mol3D_docked.pkl")

print(df, len(df.dropna()))
