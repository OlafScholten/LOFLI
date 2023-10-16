# to run:    make -f GenerateDepd.make
FC = gfortran
#FCFLAGS = -fcheck=all -fno-automatic -finit-local-zero
#FCFLAGS =


GenerateDepd.exe : GenerateDepd.o
	$(FC) -o GenerateDepd.exe GenerateDepd.f90

GenerateDepd.o : GenerateDepd.f90
	$(FC)  -c GenerateDepd.f90

