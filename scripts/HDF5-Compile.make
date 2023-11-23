# to run in fortran-source directory    
FC = gfortran
#FCFLAGS = -fcheck=all -fno-automatic -finit-local-zero
#FCFLAGS =

#all: ${MainProg}

${LL_bin}/${MainProg} : ${MainProg}.o ${LL_bin}/libLOFLI.a
	$(FC) ${HDF5Compile} ${MainProg}.o -o ${LL_bin}/${MainProg} ${LOFLIlib} ${FCFLAGS} ${HDF5Link}

${MainProg}.o : ${MainProg}.f90 ${MainProgDepd}
	$(FC)  ${HDF5Compile} ${LOFLIinc} -c ${MainProg}.f90 ${FCFLAGS}

#gfortran ${HDF5Compile} ${LOFLIinc} -c ${SourceCode} ${FCFLAGS}
#gfortran  ${HDF5Compile} ${Prog}.o -o ${Prog} ${LOFLIlib} ${FCFLAGS} ${HDF5Link}
