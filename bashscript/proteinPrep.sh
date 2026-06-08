###################################################
# This script prepares the protein for docking.   #
# Prot.pdb ----This Script----> Receptor.pdbqt    #
# Receptor.pdbqt is ready for docking             #
###################################################

dir_curr=$(pwd)
cd ../Data/Docking/receptor

grep "^ATOM" $1.pdb > protein_clean.pdb
obabel protein_clean.pdb -O receptor.pdbqt -xr -p 7.4

rm -f *.pdb
cd $dir_curr
