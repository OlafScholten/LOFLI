#!/bin/bash
# 
# 
source  /home/olaf/LOFLI/ShortCuts.sh
Prog="RFI_Mitigation"
echo "Flash Folder:" ${FlashFolder}
echo "start RFI mitigation"

source ${LL_scripts}RunProgram.sh "2"

echo "Start exploration run"
nohup ./Explore.sh  >Explore.log 2>&1  &   
echo "Explore was submitted, now wait."

exit
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata
 olaf@APP01:~/kapdatalink/lightning_data/2018


