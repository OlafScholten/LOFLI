#!/bin/bash
# 
# 

# make sure 
#  source  ${LL_Base}/scripts/RunProgram.sh LOFAR-Imag  SimData.in
#was run previourly

#rm ${LL_Base}/bin/SimulateData
source  ${LL_Base}/scripts/RunProgram.sh SimulateData  Simulate.in


rm ${LL_Base}/bin/LOFAR-Imag

source  ${LL_Base}/scripts/RunProgram.sh LOFAR-Imag  SimIntf.in

exit
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata


