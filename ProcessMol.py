from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize


def smiles_to_mol(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        return mol
    except Exception:
        return None


def embed_3d(mol):
    if mol is None:
        return None

    try:
        # 1. Strip Salts (Keep the largest fragment)
        chooser = rdMolStandardize.LargestFragmentChooser()
        clean_mol = chooser.choose(mol)

        # 2. Neutralize the naked charges!
        uncharger = rdMolStandardize.Uncharger()
        neutral_mol = uncharger.uncharge(clean_mol)

        # 3. Add explicit hydrogens
        mol_h = Chem.AddHs(neutral_mol)

        # 4. Generate Initial 3D Coordinates (ETKDGv3)
        params = AllChem.ETKDGv3()  # type: ignore
        params.randomSeed = 42
        if AllChem.EmbedMolecule(mol_h, params) == -1:  # type: ignore
            return None  # If it can't even guess a 3D shape, drop it

        # 5. The Fallback Optimization Pipeline
        # Try MMFF94 first (Best for drugs)
        if AllChem.MMFFOptimizeMolecule(mol_h, maxIters=5000) != -1:  # type: ignore
            return mol_h

        # Fallback to UFF (More generic periodic table coverage)
        elif AllChem.UFFOptimizeMolecule(mol_h, maxIters=5000) != -1:  # type: ignore
            return mol_h

        # If both fail (exotic sulfur/metals), return the ETKDGv3 shape anyway!
        else:
            return mol_h

    except Exception:
        return None
