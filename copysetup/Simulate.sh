#!/bin/bash
# 
# 
source ../ShortCuts.sh

## start Simulate run
Prog="SimulateData"
source ${UtilDir}compile.sh
${ProgramDir}$Prog  <Simulate.in
#
#exit

## start Interferometry run
Prog="LOFAR-Imag"
source ${UtilDir}compileHDF5program.sh
#${ProgramDir}/$Prog  ${FlashFolder} <SimData.in
${ProgramDir}$Prog  ${FlashFolder} <SimIntf.in

exit
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata


