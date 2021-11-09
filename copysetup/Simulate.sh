#!/bin/bash
# 
# 
source ../ShortCuts.sh
Prog="SimulateData"
source ${UtilDir}compile.sh

# start Simulate run
${ProgramDir}$Prog  <Simulate.in
exit
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata


