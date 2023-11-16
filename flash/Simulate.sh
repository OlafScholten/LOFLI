#!/bin/bash
# 
# 
source  /home/olaf/LOFLI/ShortCuts.sh
Prog="SimulateData"

source ${LL_scripts}RunProgram.sh "Simulate.in"

## start Interferometry run
Prog="LOFAR-Imag"
source ${LL_scripts}RunProgram.sh "SimIntf.in"
#source ${LL_scripts}RunProgram.sh "SimData.in"

exit
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata


