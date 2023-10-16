# Define some importants paths
# If necessary recompile program
echo "compile code:" ${Prog}

#ProgramDir="../test/"
export LOFLIinc="-I./modules"
export LOFLIlib="-lm ./libLOFLI.a ${FFTLIB}"
export FCFLAGS="-fcheck=bounds"

cd ${ProgramDir}
# make sure library is updates, rebuild library when modules have been changed
#make -f f90split.make
make -f UpdateLibraryLOFLI.make
#
export MainProg=${Prog}
export MainProgDepd="???????.f90"
# the dependencies (including extensions)of the code should be written in:
#source ${Prog}_depd.sh
#make -f GenerateDepd.make
./GenerateDepd.exe ${Prog}.f90
make -f HDF5-Compile.make

## for  HF5Compilation with automatic search for dependencies
#prefix="/usr"
#PARALLEL=no             # Am I in parallel mode?
#AR="ar"
#RANLIB="ranlib"
#H5TOOL_BIN="${prefix}/bin/h5fc"   # The path of the tool binary
####LIBRARY="-lm /home/olaf/NumLib/bin/libfftpack5.1d.a"
# ${H5TOOL_BIN} -o $Prog $SourceCode ${LOFLILIB} ${FCFLAGS}
#${H5TOOL_BIN} -o $Prog $SourceCode ${LOFLILIB}  ${FCFLAGS} -echo  # gives something like
#gfortran ${HDF5Compile} -I./modules -c ${SourceCode} ${FCFLAGS}
#gfortran  ${HDF5Compile} -I./modules ${Prog}.o -o ${Prog} -lm ./libLOFLI.a ${LIBRARY} ${FCFLAGS} ${HDF5Link}

# End compilation, return to the FlashFolder
cd ${FlashFolder}
