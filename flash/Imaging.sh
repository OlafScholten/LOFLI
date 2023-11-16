#!/bin/bash
# 
# 
source  /home/olaf/LOFLI/ShortCuts.sh
Prog="LOFAR-Imag"

# start Impulsive-Imager run
source ${LL_scripts}RunProgram.sh "Imaging.in"
exit
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata


