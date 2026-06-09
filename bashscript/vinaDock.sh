###################################################
# This script docks the ligand using Vina.        #
# It uses vinardo scoring                         #
#                                                 #
###################################################

curr=$(pwd)
cd ../Data/Docking/dock

vina --receptor ../receptor/receptor.pdbqt \
     --ligand ../ligand/ligand.pdbqt \
     --center_x $1 --center_y $2 --center_z $3 \
     --size_x $4 --size_y $5 --size_z $6 \
     --scoring vinardo \
     --exhaustiveness 10 \
     --num_modes 9 > log.txt

grep "   1      " log.txt > 1.txt
grep "   2      " log.txt > 2.txt
grep "   3      " log.txt > 3.txt

cd $curr
