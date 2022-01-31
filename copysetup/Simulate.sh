#!/bin/bash
# 
# 
source ../ShortCuts.sh
Prog="SimulateData"
source ${UtilDir}compile.sh

# start Simulate run
${ProgramDir}$Prog  <Simulate.in
#
#exit
# start Interferometry run
Prog="LOFAR-Imag"
source ${UtilDir}/compile.sh
${ProgramDir}/$Prog  ${FlashFolder} <SimIntf.in

exit
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata


