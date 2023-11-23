# to run:    make -f GenerateDepd.make
# run in any directory
FC = gfortran
#FCFLAGS = -fcheck=all -fno-automatic -finit-local-zero
#FCFLAGS =


${LL_bin}/GenerateDepd : GenerateDepd.o
	$(FC) -o ${LL_bin}/GenerateDepd GenerateDepd.o

GenerateDepd.o : ${LL_src}/GenerateDepd.f90
	$(FC)  -c ${LL_src}/GenerateDepd.f90

