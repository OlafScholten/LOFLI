# to run:    make -f f90split.make
FC = gfortran
#FCFLAGS = -fcheck=all -fno-automatic -finit-local-zero
#FCFLAGS =

#all: ${MainProg}

${MainProg} : ${MainProg}.o libLOFLI.a
	$(FC) ${HDF5Compile} ${MainProg}.o -o ${MainProg} ${LOFLIlib} ${FCFLAGS} ${HDF5Link}

${MainProg}.o : ${MainProg}.f90 ${MainProgDepd}
	$(FC)  ${HDF5Compile} ${LOFLIinc} -c ${MainProg}.f90 ${FCFLAGS}

#gfortran ${HDF5Compile} ${LOFLIinc} -c ${SourceCode} ${FCFLAGS}
#gfortran  ${HDF5Compile} ${Prog}.o -o ${Prog} ${LOFLIlib} ${FCFLAGS} ${HDF5Link}
