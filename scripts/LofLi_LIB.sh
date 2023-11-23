#!/bin/bash
#
#source = fftpack5.1d
mkdir temp
cd temp
~/NumLib/bin/f90split ../fftpack5.1d.f90
#
for FILE in `ls -1 *.f90`;
do
  gfortran -c -O3 $FILE -fPIC
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
done
rm *.f90
#
ar cr libfftpack5.1d.a *.o
rm *.o
#
mv libfftpack5.1d.a ~/NumLib/bin
cd ..
rmdir temp
#
echo "Library installed as ~/NumLib/bin/$ARCH/libfftpack5.1.a."
