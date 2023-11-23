# to run:    make -f f90split.make
# run in any directory
FC = gfortran
#FCFLAGS = -fcheck=all -fno-automatic -finit-local-zero -static-libgfortran
#FCFLAGS =


${LL_bin}/f90split : f90split_OS.o
	$(FC)  -o ${LL_bin}/f90split f90split_OS.o  

f90split_OS.o : ${LL_src}/f90split_OS.f90
	$(FC)  -c ${LL_src}/f90split_OS.f90

