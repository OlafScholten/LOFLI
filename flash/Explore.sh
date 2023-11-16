#!/bin/bash
# 
# 
source  /home/olaf/LOFLI/ShortCuts.sh
Prog="LOFAR-Imag"

# start Exploration run
source ${LL_scripts}RunProgram.sh "Explore.in"
exit
# read -rsp $'Press enter to continue...\n'
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata
