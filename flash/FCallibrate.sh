#!/bin/bash
# 
# 
source  /home/olaf/LOFLI/ShortCuts.sh
Prog="LOFAR-Imag"

# start Field Calibration run
source ${LL_scripts}RunProgram.sh "FCallibrate.in"
exit
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata


