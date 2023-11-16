#!/bin/bash
# 
# 
source  ${LL_BaseDir}/ShortCuts.sh
Prog="LOFAR-Imag"

# start TRI-D Interferometry run
source ${LL_scripts}RunProgram.sh "Interferometry.in"
exit
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata


