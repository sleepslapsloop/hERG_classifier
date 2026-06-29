import numpy as np
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Descriptors3D,
    Lipinski,
    rdFingerprintGenerator,
)
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
        chooser = rdMolStandardize.LargestFragmentChooser()
        clean_mol = chooser.choose(mol)

        uncharger = rdMolStandardize.Uncharger()
        neutral_mol = uncharger.uncharge(clean_mol)

        mol_h = Chem.AddHs(neutral_mol)

        params = AllChem.ETKDGv3()  # type: ignore
        params.randomSeed = 42
        if AllChem.EmbedMolecule(mol_h, params) == -1:  # type: ignore
            return None
        if AllChem.MMFFOptimizeMolecule(mol_h, maxIters=5000) != -1:  # type: ignore
            return mol_h
        elif AllChem.UFFOptimizeMolecule(mol_h, maxIters=5000) != -1:  # type: ignore
            return mol_h
        else:
            return mol_h

    except Exception:
        return None


def get_morgan_fingerprint(mol, radius=2):
    if mol is None:
        return [0] * 2048

    mfp1gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius)
    fp = mfp1gen.GetFingerprint(mol)
    return list(fp)


def compute_descriptors(mol):
    if mol is None:
        return None
    if mol.GetNumConformers() == 0:
        return None

    features = {
        "MolWt": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "RotB": Descriptors.NumRotatableBonds(mol),
        "FormalCharge": Chem.GetFormalCharge(mol),
        "RadiusGyr": Descriptors3D.RadiusOfGyration(mol),
        "Asphericity": Descriptors3D.Asphericity(mol),
        "Spherocity": Descriptors3D.SpherocityIndex(mol),
        "PMI1": Descriptors3D.PMI1(mol),
        "PMI2": Descriptors3D.PMI2(mol),
        "PMI3": Descriptors3D.PMI3(mol),
    }

    descriptor_values = list(features.values())
    fp_values = get_morgan_fingerprint(mol, radius=2)

    return np.array(descriptor_values + fp_values)  # type: ignore
