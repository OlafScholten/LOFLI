#!/bin/bash
# 
# 
source  ${LL_Base}/ShortCuts.sh
export Prog=$1
export input=$2

echo "program run is: ${Prog} with argument ${input}" 

if ! test -f ${LL_bin}/${Prog}; then
   echo ${Prog} "does not exist and will be generated."
   source ${LL_scripts}/compileHDF5program.sh
fi

if test -f ${input}; then
   echo "File ${input} exists"
   echo "executing command: ${LL_bin}/${Prog}  <${input} "
   ${Prog}  <${input}
else
   echo "arguments ${input} is not an existing file"
   echo "executing command: ${LL_bin}/${Prog}  ${input} "
   ${Prog}  ${input}
fi
