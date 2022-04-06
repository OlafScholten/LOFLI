#!/bin/bash
# 
# 
source ../ShortCuts.sh
Prog="PowerTrace"
source ${UtilDir}/compile.sh

# start Exploration run
${ProgramDir}/$Prog  <PwT.in
#gle -d pdf -o Map_Explore.pdf ${UtilDir}/LeaderTrack.gle ${FlashFolder}/Explore
exit
# read -rsp $'Press enter to continue...\n'
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata
#GLE_USRLIB=/home/olaf/NumLib/gle
#gle -d pdf -o Map_Explore.pdf ${GLE_USRLIB}/LeaderTrack.gle ${PWD}/Explore
