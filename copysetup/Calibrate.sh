#!/bin/bash
# 
# 
source ../ShortCuts.sh
#Prog="LOFAR-Imag"${Vrsn}
Prog="LOFAR-Imag"
source ${UtilDir}compileHDF5program.sh

# start Calibration run
${ProgramDir}$Prog  ${FlashFolder} <Calibrate.in
exit
# read -rsp $'Pausing. Press enter to continue...\n'
