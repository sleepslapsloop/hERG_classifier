import subprocess
from dataclasses import dataclass
from pathlib import Path

# Gets the absolute path of the directory containing master.py
BASE_DIR = Path(__file__).parent.resolve()


@dataclass(frozen=True)
class hERG:
    center_of_mass: tuple[float, float, float] = (143.219, 153.455, 132.858)
    center_x: float = center_of_mass[0]
    center_y: float = center_of_mass[1]
    center_z: float = center_of_mass[2]


def _run_script(script_name: str, args: list):
    """Internal helper function to execute bash scripts securely."""
    script_path = BASE_DIR / script_name
    command = ["bash", str(script_path)] + [str(arg) for arg in args]

    try:
        # cwd=BASE_DIR ensures the relative paths in the bash scripts (like ../Data/)
        # don't break when master.py is imported and run from a completely different folder.
        subprocess.run(
            command, capture_output=True, text=True, check=True, cwd=BASE_DIR
        )
    except subprocess.CalledProcessError as e:
        print(f"❌ Error in {script_name}:\n{e.stderr}")
        raise


def proteinPrep(inputName: str):
    """Prepares the protein for docking."""
    print(f"Preparing protein: {inputName}...")
    _run_script("proteinPrep.sh", [inputName])


def ligPrep(smiles: str):
    """Prepares the ligand for docking from a SMILES string."""
    print("Preparing ligand...")
    _run_script("ligPrep.sh", [smiles])


def dock(
    center_x: float,
    center_y: float,
    center_z: float,
    size_x: float,
    size_y: float,
    size_z: float,
):
    """Runs Vina docking with the specified grid box coordinates and dimensions."""
    print("Running Vina docking...")
    args = [center_x, center_y, center_z, size_x, size_y, size_z]
    _run_script("vinaDock.sh", args)


def get_top_energies() -> list:
    """Reads 1.txt/2.txt/3.txt and returns a list of the top 3 binding energies."""
    dock_dir = (BASE_DIR / "../Data/Docking/dock").resolve()
    top_energies = []

    for i in range(1, 4):
        file_path = dock_dir / f"{i}.txt"
        try:
            with open(file_path, "r") as file:
                line = file.read()
                if line:
                    line = line.strip()
                    energy = float(line[8:14].strip())
                    top_energies.append(energy)
        except FileNotFoundError:
            continue

    return top_energies


def cleanup():
    """
    Deletes all temporary and output files generated during a docking run:
      - log.txt, 1.txt, 2.txt, 3.txt  (from ../Data/Docking/dock/)
      - All ligand files               (from ../Data/Docking/ligand/)
      - Vina *_out.pdbqt output files  (from the dock directory)

    The receptor (../Data/Docking/receptor/) is intentionally left untouched.
    """
    dock_dir = (BASE_DIR / "../Data/Docking/dock").resolve()
    ligand_dir = (BASE_DIR / "../Data/Docking/ligand").resolve()

    deleted = []
    missing = []

    # --- Dock directory: log + mode result files ---
    for filename in ["log.txt", "1.txt", "2.txt", "3.txt"]:
        f = dock_dir / filename
        if f.exists():
            f.unlink()
            deleted.append(str(f))
        else:
            missing.append(str(f))

    # --- Dock directory: Vina *_out.pdbqt files ---
    for f in dock_dir.glob("*_out.pdbqt"):
        f.unlink()
        deleted.append(str(f))

    # --- Ligand directory: all files ---
    if ligand_dir.exists():
        for f in ligand_dir.iterdir():
            if f.is_file():
                f.unlink()
                deleted.append(str(f))

    # if deleted:
    #     print(f"🧹 Cleaned up {len(deleted)} file(s).")
    # if missing:
    #     print(f"⚠️  {len(missing)} expected file(s) were already absent: {missing}")


if __name__ == "__main__":
    proteinPrep("8ZYP")
    ligPrep("Fc1ccc(cc1)Cn2c5ccccc5nc2NC4CCN(CCc3ccc(OC)cc3)CC4")
    dock(
        143.219, 153.455, 132.858, 20, 20, 20
    )  # Center of Mass of bound ligand in 8ZYP.cif: [ 143.219, 153.455, 132.858]

    energies = get_top_energies()
    print(f"Top 3 binding modes: {energies}")
    cleanup()
