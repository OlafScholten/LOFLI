#!/bin/bash
# 
# 
# start TRI-D Interferometry run

rm ${LL_Base}/bin/AirplaneSources

source  ${LL_Base}/scripts/RunProgram.sh AirplaneSources  AirplaneSources.in
exit
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata


