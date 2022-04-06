#!/bin/bash
# 
# 
source ../ShortCuts.sh
Prog="LOFAR-Imag"
source ${UtilDir}/compile.sh

# start Imaging run
${ProgramDir}/$Prog  ${FlashFolder} <PkIntf.in
#./GLE-plots.sh
exit
#gle -d pdf plot-corr_01.gle
#gle -d pdf plot-corr_02.gle
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata


