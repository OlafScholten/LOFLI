# Define some importants paths
# If necessary recompile program
echo "compile code:" ${Prog}
prefix="/usr/local/hdf5"
PARALLEL=no             # Am I in parallel mode?
AR="ar"
RANLIB="ranlib"
H5TOOL_BIN="${prefix}/bin/h5fc"   # The path of the tool binary
#LIBRARY="-lm ~/NumLib/bin/libfftpack5.1d.a"
SourceCode="${Prog}.f90"
cd ${ProgramDir}
${H5TOOL_BIN} -o $Prog $SourceCode ${LIBRARY}
# End compilation, return to the FlashFolder
cd ${FlashFolder}
