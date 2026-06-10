from typing import cast

import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

from bashscript.master import (
    cleanup,
    dock,
    get_top_energies,
    hERG,
    ligPrep,
    proteinPrep,
)

# Prepare the protein for docking
proteinPrep("8ZYP")


def largest_fragment(smiles: str) -> str:
    """
    Strips counterions and salts, returning the SMILES of the largest fragment.
    E.g. "molecule.[Br-]"  ->  "molecule"
         "molecule.Cl.Cl"  ->  "molecule"
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles  # let ligPrep fail and be caught downstream
    chooser = rdMolStandardize.LargestFragmentChooser()
    largest = chooser.choose(mol)
    return Chem.MolToSmiles(largest)


def compute_energies(smiles: str) -> tuple[float, ...] | None:
    """
    Returns a 3-tuple of binding energies, or None if any step fails/times out.
    Always runs cleanup regardless of success or failure.
    """
    clean_smiles = largest_fragment(smiles)
    try:
        ligPrep(clean_smiles)
        dock(hERG.center_x, hERG.center_y, hERG.center_z, 25, 25, 25)
        energies = get_top_energies()
        return tuple(energies) if energies else None
    except Exception as e:
        print(f"  ⚠️  Skipping — {e}")
        return None
    finally:
        cleanup()  # always runs


df = cast(pd.DataFrame, pd.read_pickle("Data/Mol3D/hERG_ChEMBL_Mol3D.pkl"))

print(df.iloc[79]["Smiles"])

OUTPUT_PATH = "Data/Docked/hERG_ChEMBL_Mol3D_docked.pkl"

MODE_COLS = ["energy_mode1", "energy_mode2", "energy_mode3"]
try:
    df_out = pd.read_pickle(OUTPUT_PATH)
    print(
        f"▶️  Resuming — {df_out[MODE_COLS[0]].notna().sum()} / {len(df_out)} rows already done."
    )
except FileNotFoundError:
    df_out = df.copy()
    for col in MODE_COLS:
        df_out[col] = None


pending = df_out[df_out[MODE_COLS[0]].isna()]
total = len(pending)
print(f"📋  {total} molecules to dock.\n")

for count, (idx, row) in enumerate(pending.iterrows(), 1):
    print(f"[{count}/{total}] idx={idx}  {row['Smiles'][:60]}...")
    result = compute_energies(row["Smiles"])  # type: ignore

    if result and len(result) == 3:
        df_out.loc[idx, MODE_COLS] = list(result)
    else:
        # Store NaN explicitly so resume skips it next time
        df_out.loc[idx, MODE_COLS] = [float("nan")] * 3

    # Checkpoint every 50 molecules so a crash doesn't lose all progress
    if count % 50 == 0:
        df_out.to_pickle(OUTPUT_PATH)
        print(f"  💾  Checkpoint saved ({count} processed).")

# Final derived columns and save
df_out[MODE_COLS] = df_out[MODE_COLS].astype(float)
df_out["energy_range"] = df_out[MODE_COLS].max(axis=1) - df_out[MODE_COLS].min(axis=1)
df_out["energy_mean"] = df_out[MODE_COLS].mean(axis=1)
df_out["energy_std"] = df_out[MODE_COLS].std(axis=1)

df_out.to_pickle(OUTPUT_PATH)
print(f"\n✅  Done. Saved to {OUTPUT_PATH}")
print(df_out[MODE_COLS + ["energy_mean", "energy_std"]].describe())


# mode_cols = ["energy_mode1", "energy_mode2", "energy_mode3"]
# df[mode_cols] = pd.DataFrame(
#     df["Smiles"].apply(compute_energies).tolist(), index=df.index
# )

# df["energy_range"] = df[mode_cols].max(axis=1) - df[mode_cols].min(axis=1)
# df["energy_mean"] = df[mode_cols].mean(axis=1)
# df["energy_std"] = df[mode_cols].std(axis=1)

# df.to_pickle("Data/Docked/hERG_ChEMBL_Mol3D_docked.pkl")

# print(df)
