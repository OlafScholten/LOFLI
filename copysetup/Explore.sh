#!/bin/bash
# 
# 
source ../ShortCuts.sh
#Prog="LOFAR-Imag"${Vrsn}
Prog="LOFAR-Imag"
source ${UtilDir}compile.sh

# start Exploration run
${ProgramDir}$Prog  <Explore.in
exit
# read -rsp $'Press enter to continue...\n'
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata
