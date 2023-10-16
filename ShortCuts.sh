# Define some importants shortcuts
ProgramDir="../Program/"
UtilDir="../Utilities/"
# point to directory where HDF5 utilities are stored
export HDF5lib="/usr/lib/x86_64-linux-gnu/hdf5/serial"
export HDF5Compile="-g -O2 -fdebug-prefix-map=/build/hdf5-X9JKIg/hdf5-1.10.0-patch1+docs=. -fstack-protector-strong -I/usr/include/hdf5/serial"
export HDF5Link="-L${HDF5lib} ${HDF5lib}/libhdf5hl_fortran.a ${HDF5lib}/libhdf5_hl.a ${HDF5lib}/libhdf5_fortran.a ${HDF5lib}/libhdf5.a -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,${HDF5lib}"

export FFTLIB="-lm /home/olaf/NumLib/bin/libfftpack5.1d.a"  # FFT library, double precision

export AntennaFun="../AntenFunct/v2-"
#export AntennaFun="../AntenFunct/"
export AntennaFieldsDir="../AntennaFields/"
FlashFolder=$(pwd)
Vrsn=""  # version number
#  It is essential not to have spaces around the = sign in the previous commands
echo "FlashFolder:" ${FlashFolder}
echo "Antenna:" ${AntennaFieldsDir}
