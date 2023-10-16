#!/bin/bash
# 
# 
source ../ShortCuts.sh
#Prog="LOFAR-Imag"${Vrsn}
Prog="LOFAR-Imag"
#source ${UtilDir}compile.sh
source ${UtilDir}compileHDF5program.sh

# start Imaging run
${ProgramDir}$Prog  ${FlashFolder} <Imaging.in
exit
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata


