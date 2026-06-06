from rdkit import Chem
from rdkit.Chem import AllChem


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
    try:
        params = AllChem.ETKDGv3()  # type: ignore
        params.randomSeed = 42
        AllChem.EmbedMolecule(mol, params)  # type: ignore
        AllChem.UFFOptimizeMolecule(mol)  # type: ignore
        return mol
    except Exception:
        return None
