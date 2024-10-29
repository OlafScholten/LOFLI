#!/bin/bash 
# 
# 

cd modules
echo A-Me-Hoela [=NoProblem]: $1

ErrorFlag=0

for FileA in $1; do
	echo $'\n' "******" ${FileA} $'\n'

	f90split ../${FileA}
	#
	for FILE in `ls -1 *.f90`; do
	  gfortran ${HDF5Compile} -c -O ${FILE} ${FCFLAGS}
	  if [ $? -ne 0 ]; then
		echo "Errors compiling " ${FILE}
      ErrorFlag=1
		#break
	  fi
	done
	rm *.f90

done

#echo "ErrorFlag=" ${ErrorFlag}

if [ ${ErrorFlag} -eq 0 ]; then
   ar rcsv ${LL_bin}/libLOFLI.a *.o
   #break
else
   echo "Library not updated due to compiling errors " 
fi

rm *.o

cd ..
#rmdir temp
