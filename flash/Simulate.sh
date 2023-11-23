#!/bin/bash
# 
# 

source  ${LL_Base}/scripts/RunProgram.sh SimulateData  Simulate.in


source  ${LL_Base}/scripts/RunProgram.sh LOFAR-Imag  SimIntf.in

exit
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata


