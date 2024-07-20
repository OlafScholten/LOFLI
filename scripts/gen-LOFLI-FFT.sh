#!/bin/bash
#
# source = fftpack5.1d should be in  ${FFTPackBase}/FFTPACK
# Make sure the folder  ${FFTPackBase}/bin/  exists

cd ${LL_src}
make -f ${LL_scripts}/f90split.make

cd ${FFTPackBase}/FFTPACK
mkdir temp
cd temp
${LL_bin}/f90split ../fftpack5.1d.f90
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
mv libfftpack5.1d.a ${FFTPackBase}/bin
cd ..
rmdir temp
#
echo "Library installed as ${FFTPackBase}/bin/libfftpack5.1.a."
