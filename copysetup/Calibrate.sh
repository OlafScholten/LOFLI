#!/bin/bash
# 
# 
source ../ShortCuts.sh
#Prog="LOFAR-Imag"${Vrsn}
Prog="LOFAR-Imag"
source ${UtilDir}compile.sh

# start Calibration run
${ProgramDir}$Prog  ${FlashFolder} <Calibrate.in
exit
# read -rsp $'Press enter to continue...\n'
