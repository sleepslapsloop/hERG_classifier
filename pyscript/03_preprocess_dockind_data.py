from bashscript.master import (
    cleanup,
    dock,
    get_top_energies,
    hERG,
    ligPrep,
    proteinPrep,
)

# Prepare the protein for docking
# Can be run multiple times, but running it once is enough
proteinPrep("8ZYP")


def compute_energies(smiles: str):
    ligPrep(smiles)
    dock(
        hERG.center_x, hERG.center_y, hERG.center_z, 25, 25, 25
    )  # box size of 25x25x25 = 15,625 Angstroms^3
    energies = get_top_energies()
    cleanup()
