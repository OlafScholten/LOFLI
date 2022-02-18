#!/bin/bash
# 
# 
source ../ShortCuts.sh

## start Compare Calibrations run
Prog="CompareCalibrations."
source ${UtilDir}compile.sh
${ProgramDir}$Prog  <CompareCalibr.in
#
exit
