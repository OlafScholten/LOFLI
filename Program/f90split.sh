#!/bin/bash
#  usage ./f90split.sh
#
gfortran -c f90split_OS.f90
if [ $? -ne 0 ]; then
  echo "Errors compiling f90split_OS.f90"
  exit
fi
#
gfortran -o f90split f90split_OS.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading f90split.o"
  exit
fi
rm f90split_OS.o
#
chmod ugo+x f90split
echo "Program installed as   f90split"  

#mv f90split ../Utilities/f90split
#echo "Program installed as ../Utilities/f90split"  
