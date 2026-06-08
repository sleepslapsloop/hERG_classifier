###################################################
# This script prepares the ligand for docking.    #
# Lig.pdb ----This Script----> Ligand.pdbqt       #
# Ligand.pdbqt is ready for docking               #
###################################################

dir_curr=$(pwd)
cd ../Data/Docking/ligand

obabel -:"$1" -O ligand_temp.mol2 --gen3d --best -p 7.4
obabel ligand_temp.mol2 -O ligand.pdbqt --minimize --ff MMFF94 --steps 2000 --sd

cd $dir_curr
