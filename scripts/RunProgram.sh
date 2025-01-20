#!/bin/bash 
# 
# 
source  ${LL_Base}/ShortCuts.sh
export Prog=$1

echo "program run is: ${Prog} with argument $2" 

if ! test -f ${LL_bin}/${Prog}; then
   echo ${Prog} "does not exist and will be generated."
   source ${LL_scripts}/compileHDF5program.sh
fi

if test -f $2; then
   echo "File $2 exists"
   export input=' <'$2   
   echo "executing command: ${LL_bin}/${Prog}  ${input} "
   ${LL_bin}/${Prog}  <$2   # This is apparently very different from    ${LL_bin}/${Prog}  ${input}
else
   echo "arguments $2 is not an existing file"
   export input=$2
   echo "executing command: ${LL_bin}/${Prog}  ${input} "
   ${LL_bin}/${Prog}  $2
fi

#echo "executing command: ${LL_bin}/${Prog}  ${input} "
#${LL_bin}/${Prog}  ${input}
