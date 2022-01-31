#!/bin/bash
# 
# 
source ../ShortCuts.sh
#Prog="RFI_Mitigation"${Vrsn}
Prog="RFI_Mitigation"
source ${UtilDir}compile.sh
echo "Flash Folder:" ${FlashFolder}
echo "start RFI mitigation"
${ProgramDir}$Prog ${AntennaFieldsDir} 2
#source GLE-plots.sh
echo "Start exploration run"
nohup ./Explore.sh  >Explore.log 2>&1  &   
echo "Explore was submitted, now wait."

exit
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata
 olaf@APP01:~/kapdatalink/lightning_data/2018


