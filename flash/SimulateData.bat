@echo off
Set LL_BaseDir="C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\LOFLI"
call %LL_BaseDir%\ShortCuts.bat
set Prog=SimulateData

echo "SimulateData" functions only on linux systems
pause
exit

echo run: call %LL_scripts%\RunProgram.bat "Simulate.in"
call %LL_scripts%\RunProgram.bat "Simulate.in"

echo run: call %LL_scripts%\RunProgram.bat Simulate.in
call %LL_scripts%\RunProgram.bat Simulate.in
pause
rem @%LL_scripts%\RunSimulateData.bat Simulate.in

echo finished
pause

exit


#!/bin/bash
#
#
source ../ShortCuts.sh
Prog="SimulateData"
#source ${UtilDir}compile.sh
source ${UtilDir}compileHDF5program.sh

# start Simulate run
${ProgramDir}$Prog  <Simulate.in
#
#exit
# start Interferometry run
Prog="LOFAR-Imag"
#source ${UtilDir}/compile.sh
source ${UtilDir}compileHDF5program.sh
${ProgramDir}/$Prog  ${FlashFolder} <SimIntf.in

exit
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata
