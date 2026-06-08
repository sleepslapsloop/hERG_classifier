###################################################
# This script docks the ligand using Vina.        #
# It uses vinardo scoring                         #
#                                                 #
###################################################

vina --receptor ../receptor/receptor.pdbqt \
     --ligand ../ligand/ligand.pdbqt \
     --center_x $1 --center_y $2 --center_z $3 \
     --size_x $4 --size_y $5 --size_z $6 \
     --scoring vinardo \
     --exhaustiveness 64 \
     --num_modes 9 \
     --energy_range 5 \
     --log log.txt
