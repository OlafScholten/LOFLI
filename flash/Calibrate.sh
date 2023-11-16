#!/bin/bash
# 
# 
source  ${LL_BaseDir}/ShortCuts.sh
Prog="LOFAR-Imag"

# start Calibration run
source ${LL_scripts}RunProgram.sh "Calibrate.in"
exit
#exit
# read -rsp $'Pausing. Press enter to continue...\n'
