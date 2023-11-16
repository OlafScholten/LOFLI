#!/bin/bash
# 
# 
#  It is essential not to have spaces around the = sign in the previous commands
#echo "Flash Folder:" ${FlashFolder}
#
#FlashFolder='19A-4'
FlaID=${NewFlashFolder}
FlaLine="$(grep ${FlaID} ${LL_scripts}list.ssv)"
FlaArray=($FlaLine)
FlaUTC=${FlaArray[2]}
FlaYear=${FlaArray[1]}
export ArchiveDir=${ArchiveBase}${FlaYear}/${FlaUTC}
#
echo "Flash: "${FlaID}":"
echo " Year: "$FlaYear
echo " @UTC: "$FlaUTC":"
echo " Arch: "$ArchiveDir":"
#exit
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata
#sentence="this is a story"
#stringarray=($sentence)
#echo ${stringarray[0]}
#
#grep "19A-4" list.ssv >lineD

