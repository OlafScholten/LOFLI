#!/bin/bash
# 
# 

cd modules
echo A-Me-Hoela [=NoProblem]: $1

for FileA in $1; do
	echo $'\n' "******" ${FileA} $'\n'

	f90split ../${FileA}
	#
	for FILE in `ls -1 *.f90`; do
	  gfortran ${HDF5Compile} -c -O ${FILE} ${FCFLAGS}
	  if [ $? -ne 0 ]; then
		echo "Errors compiling " ${FILE}
		exit
	  fi
	done
	rm *.f90
	ar rcsv ${LL_bin}/libLOFLI.a *.o
	rm *.o

done
cd ..
#rmdir temp
