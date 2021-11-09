#!/bin/bash
# 
# 
source ../ShortCuts.sh
#Vrsn="-v21"  # version number
Prog="Track"${Vrsn}
source ${UtilDir}/compile.sh
#
#cd ${ProgramDir}
#gfortran -o TrackExe.exe Track${Vrsn}.f90 ${LIBRARY}
#cd ${FlashFolder}
#
${ProgramDir}/$Prog  <FlashImage.in
exit
