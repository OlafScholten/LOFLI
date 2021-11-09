echo "compile:" ${Prog}
#ProgramDir="/home/olaf/LMA-fit/LMA2019/program"
#ProgramDir="../program"
#HDF5prefix="/usr/local/hdf5"
PARALLEL=no             # Am I in parallel mode?
AR="ar"
RANLIB="ranlib"
H5TOOL_BIN="${HDF5prefix}/bin/h5fc"   # The path of the tool binary
#LIBRARY="-lm ~/NumLib/bin/libfftpack5.1d.a"
SourceCode="${Prog}.f90"
cd ${ProgramDir}
${H5TOOL_BIN} -o $Prog $SourceCode ${LIBRARY}
