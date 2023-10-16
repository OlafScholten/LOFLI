#!/bin/bash
# 
# 
source ../ShortCuts.sh
Prog="LOFAR-Imag"
source ${UtilDir}compileHDF5program.sh

# start Exploration run
${ProgramDir}$Prog  ${FlashFolder}  <Explore.in
exit
# read -rsp $'Press enter to continue...\n'
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata
