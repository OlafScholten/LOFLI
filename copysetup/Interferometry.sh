#!/bin/bash
# 
# 
source ../ShortCuts.sh
Prog="LOFAR-Imag"
#source ${UtilDir}compile.sh
source ${UtilDir}compileHDF5program.sh

# start Interferometry run
${ProgramDir}$Prog  ${FlashFolder} <Interferometry.in
#gle -d pdf -o InterfContour.pdf ${UtilDir}/InterfContour.gle ${FlashFolder}/
#source GLE-plots.sh
exit
#gle -d pdf plot-corr_01.gle
#gle -d pdf plot-corr_02.gle
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata


