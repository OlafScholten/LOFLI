#!/bin/bash
# 
# 
source ../ShortCuts.sh
Prog="LOFAR-Imag"${Vrsn}
source ${UtilDir}compile.sh

# start Imaging run
${ProgramDir}$Prog  ${FlashFolder} <Imaging.in
exit
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata


