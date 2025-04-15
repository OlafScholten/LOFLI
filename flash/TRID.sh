#!/bin/bash
# 
# 
# start TRI-D Interferometry run

#rm ${LL_Base}/bin/LOFAR-Imag

source  ${LL_Base}/scripts/RunProgram.sh LOFAR-Imag  TRID.in
exit
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata


