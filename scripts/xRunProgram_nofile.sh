#!/bin/bash
# 
# 
source ../ShortCuts.sh
echo "program run is:" ${Prog} with argument "$1" 

if ! test -f ${LL_bin}/${Prog}; then
   echo ${Prog} "does not exist and will be generated."
   source ${LL_scripts}compileHDF5program.sh
fi

echo "command: ${LL_bin}/${Prog}  $1 "
# pause
${Prog}  $1
