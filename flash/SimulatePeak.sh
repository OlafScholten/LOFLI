#!/bin/bash
# 
# 

# make sure 
#  source  ${LL_Base}/scripts/RunProgram.sh LOFAR-Imag  SimData.in
#was run previourly

rm ${LL_Base}/bin/SimulateData

source  ${LL_Base}/scripts/RunProgram.sh SimulateData  Simulate.in

rm ${LL_Base}/bin/LOFAR-Imag

source  ${LL_Base}/scripts/RunProgram.sh LOFAR-Imag  PkIntfSim.in
source  ${LL_Base}/scripts/RunProgram.sh LOFAR-Imag  PkIntfSim2.in

rm ${LL_Base}/bin/AirplaneSources

source  ${LL_Base}/scripts/RunProgram.sh AirplaneSources  AirplaneSim.in

exit
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata


