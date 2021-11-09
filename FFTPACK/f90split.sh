#!/bin/bash
#  usage ./f90split.sh
#
gfortran -c f90split.f90
if [ $? -ne 0 ]; then
  echo "Errors compiling f90split.f90"
  exit
fi
#
gfortran -o f90split f90split.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading f90split.o"
  exit
fi
rm f90split.o
#
chmod ugo+x f90split
mv f90split ~/NumLib/bin/f90split
#
echo "Program installed as ~/NumLib/bin/f90split"  
