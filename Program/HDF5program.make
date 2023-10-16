# to run:    make -f f90split.make
FC = gfortran
#FCFLAGS = -fcheck=all -fno-automatic -finit-local-zero
#FCFLAGS =
Prog = LOFAR-Imag
Prog_Include = ECallibrOption.f90
$(info 'hdf5:' ${HDF5Compile})

all: ${Prog}.exe

${Prog}.exe : ${Prog}.o
	$(FC) ${HDF5Compile} ${Prog}.o -o ${Prog} ${LOFLIlib} ${FCFLAGS} ${HDF5Link}

${Prog}.o : ${Prog}.f90 ${Prog_Include}
	$(FC)  ${HDF5Compile} ${LOFLIinc} -c ${SourceCode} ${FCFLAGS}

#gfortran ${HDF5Compile} ${LOFLIinc} -c ${SourceCode} ${FCFLAGS}
#gfortran  ${HDF5Compile} ${Prog}.o -o ${Prog} ${LOFLIlib} ${FCFLAGS} ${HDF5Link}
