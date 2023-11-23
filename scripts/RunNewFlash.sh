#!/bin/bash
# 
# 
 source  ${LL_Base}/ShortCuts.sh
# Create a new flash-bub folder inside the present directory

ParentDir=${LL_Base}/flash/
#FlashFolder='18F-1'
#FlashFolder='18D-1'
#FlashFolder='20B-10'
#FlashFolder='21C-1'
#FlashFolder='21C-9'
#FlashFolder='21C-7'
#FlashFolder='21D-1'
#FlashFolder='21D-3'
#FlashFolder='21D-4'
#FlashFolder='21E-1'
NewFlashFolder=$1


#Year="2020"
#D20200814T143241.768Z
#FlashFolder="D20200814T143241"
#Extension="000Z"
#Year="2017"
#D20190424T210306.154Z
#FlashFolder="D20190424T210306"
#ArchiveDir="/home/olaf/kaptdata/lightning_data/2019/"${FlashFolder}".154Z"
#FlashFolder="D20190424T213055"
#ArchiveDir="/home/olaf/kaptdata/lightning_data/2019/"${FlashFolder}".182Z   !!!!/2018/D20180813T153001.413Z
#ArchiveDir="/home/olaf/kaptdata/lightning_data/"${Year}"/"${FlashFolder}"."${Extension}


MainDir=$(pwd)
echo Present directory= ${MainDir}
echo Copy ${ParentDir} to folder ${MainDir}/${NewFlashFolder}
source ${LL_scripts}/FlashID.sh

mkdir ${NewFlashFolder}
cd ${NewFlashFolder}
mkdir files
mkdir Book

cp ${ParentDir}/Book/{StationCalibrations.dat,Calibrations*.dat,*.cal} ${MainDir}/${NewFlashFolder}/Book
cp ${ParentDir}/{*.sh,*.in,*.bat} ${MainDir}/${NewFlashFolder}
 cd ${NewFlashFolder}
chmod 764 *.sh
ls -l
# copying is finished now, time to become productive
ls ${ArchiveDir}/*.h5  >directory.out
#ls ~/kaptdata/lightning_data/2019/D20190424T210306.154Z/*.h5  >directory.out

nohup ./RFI_Mitigation.sh  >RFI_Mitigation.log 2>&1  &   
echo RFI_Mitigation was submitted, now wait.
 exit
#read -rsp $'Press enter to continue...\n'
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata
# olaf@APP01:~/kapdatalink/lightning_data/2018


