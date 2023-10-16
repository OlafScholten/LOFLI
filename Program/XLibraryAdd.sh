#!/bin/bash
# 
# 
#mkdir temp
cd temp

echo Update library with $1
#echo nrvars $#
for FileA in $1; do
  #echo filea $FileA
  #echo -e '\n fileA $FileA'
#for FileA in $UpdateFilesf90 ;
#do
	echo $'\n' "working on " $FileA $'\n'

	# FileA := fftpack5.1d.f90
	# LibraryUpdate := ~/NumLib/bin/libfftpack5.1d.a
#	~/NumLib/bin/f90split ../${FileA}
	../../Utilities/f90split ../${FileA}
	#
	for FILE in `ls -1 *.f90`;
	do
	  gfortran -c -O $FILE -fPIC
	  if [ $? -ne 0 ]; then
		echo "Errors compiling " $FILE
		exit
	  fi
	done
#	rm *.f90
	#
	ar rcsv ../libLOFLI.a *.o
#	rm *.o
#
done
cd ..
#rmdir temp
