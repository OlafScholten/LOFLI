#!/bin/bash
# 
# 
echo "program run is:" ${Prog} with argument "$1" 

if ! test -f ${LL_bin}/${Prog}; then
   echo ${Prog} "does not exist and will be generated."
   source ${LL_scripts}compileHDF5program.sh
fi

if test -f $1; then
   echo "File $1 exists"
   echo "executing command: ${LL_bin}/${Prog}  <$1 "
   ${Prog}  <$1
else
   echo "arguments $1 is not an existing file"
   echo "executing command: ${LL_bin}/${Prog}  $1 "
   ${Prog}  $1
fi
