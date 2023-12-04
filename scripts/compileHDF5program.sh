#!/bin/bash
# 
# 
# Define some importants paths
#source  /home/olaf/LOFLI/ShortCuts.sh

# If necessary recompile program
echo "compile code:" ${Prog}
cd ${LL_src}



# make sure library is updated, rebuild library when modules have been changed.
#  f90split is used (a home brew version) and un-comment the following line to 
#  make the code if not present already
#cd ${LL_mod}
make -f ${LL_scripts}/f90split.make
#  If "GenerateDepd.exe" does not exist, configure it by un-commenting the following line.
make -f ${LL_scripts}/GenerateDepd.make

# it is important to have a space after each filename in "OrderedLibraryFiles.lst"
export LibrarySources=$(<${LL_scripts}/OrderedLibraryFiles.lst)
#echo "content of OrderedLibraryFiles: $LibrarySources"
make -f ${LL_scripts}/UpdateLibraryLOFLI.make

export MainProg=${Prog}
#  "MainProgDepd" is set to a strange value to make sure it is re-generated later this script.
export MainProgDepd="???????.f90"
# the dependencies (including extensions) of the code will be written in "MainProgDepd".
#  The following code-call searches the include statements in the file "${Prog}.f90" and 
#  makes the list of dependencies (in file *_depd.sh) that will be checked for updates in 
#  the subsequent compilation by the "make" command.
GenerateDepd ${Prog}.f90
# This generated the dependencies list, included by the following shell
source ${Prog}_depd.sh
echo "Dependencies for" ${Prog} "are" ${MainProgDepd}
make -f ${LL_scripts}/HDF5-Compile.make

#Clean
rm *.o
rm *.mod
# End compilation, return to the FlashFolder
cd ${FlashFolder}
