import subprocess
from pathlib import Path

# Gets the absolute path of the directory containing master.py
BASE_DIR = Path(__file__).parent.resolve()


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
    """Reads log.txt and returns a list of the top 3 binding energies."""
    # Based on vinaDock.sh, the log is saved in the same directory as execution
    log_path = BASE_DIR / "log.txt"
    top_energies = []

    try:
        with open(log_path, "r") as file:
            lines = file.readlines()

        # Vina's log outputs a table. We need to find where the separator line is.
        table_start_index = -1
        for i, line in enumerate(lines):
            if line.startswith("-----+------------+----------+----------"):
                table_start_index = i + 1
                break

        if table_start_index != -1:
            # Grab up to 3 lines immediately following the table header
            for i in range(table_start_index, min(table_start_index + 3, len(lines))):
                columns = lines[i].split()
                if len(columns) >= 2:
                    try:
                        # The second column is the affinity (energy)
                        energy = float(columns[1])
                        top_energies.append(energy)
                    except ValueError:
                        continue

        return top_energies

    except FileNotFoundError:
        print(f"❌ Could not find Vina log file at {log_path}")
        return []


def get_top_energies_and_cleanup() -> list:
    """Reads log.txt, extracts top 3 energies, and cleans up docking files."""
    log_path = BASE_DIR / "log.txt"
    top_energies = []

    # 1. Extract the energies
    try:
        with open(log_path, "r") as file:
            lines = file.readlines()

        table_start_index = -1
        for i, line in enumerate(lines):
            if line.startswith("-----+------------+----------+----------"):
                table_start_index = i + 1
                break

        if table_start_index != -1:
            for i in range(table_start_index, min(table_start_index + 3, len(lines))):
                columns = lines[i].split()
                if len(columns) >= 2:
                    try:
                        top_energies.append(float(columns[1]))
                    except ValueError:
                        continue
    except FileNotFoundError:
        print(f"❌ Could not find Vina log file at {log_path}")

    # 2. Safely Clean Up the Workspace
    print("🧹 Cleaning up temporary and output files...")

    # Delete the log.txt file
    if log_path.exists():
        log_path.unlink()

    # Delete all generated ligand files (ligand.pdbqt, ligand_temp.mol2)
    ligand_dir = (BASE_DIR / "../Data/Docking/ligand").resolve()
    if ligand_dir.exists():
        for file in ligand_dir.glob("*"):
            if file.is_file():
                file.unlink()

    # Delete Vina PDBQT output files (Vina usually creates an *_out.pdbqt file)
    # This checks the directory where master.py was run from
    for out_file in BASE_DIR.glob("*_out.pdbqt"):
        out_file.unlink()

    # We do NOT touch ../Data/Docking/receptor, keeping receptor.pdbqt safe!

    return top_energies


# Optional: You can keep a test block here.
# Code inside here won't run when you import this file from another script.
if __name__ == "__main__":
    proteinPrep("my_target")
    ligPrep("CC(=O)OC1=CC=CC=C1C(=O)O")
    dock(10.5, 12.1, -5.4, 20, 20, 20)

    energies = get_top_energies()
    print(f"Top 3 binding modes: {energies}")
