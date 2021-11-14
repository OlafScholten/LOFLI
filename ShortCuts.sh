# Define some importants shortcuts
ProgramDir="../Program/"
UtilDir="../Utilities/"
# point to directory where HDF5 utilities are stored
HDF5prefix="/usr/local/hdf5"   
LIBRARY="-lm ~/NumLib/bin/libfftpack5.1d.a"
export AntennaFun="../AntenFunct/v2-"
#export AntennaFun="../AntenFunct/"
export AntennaFieldsDir="../AntennaFields/"
FlashFolder=$(pwd)
Vrsn=""  # version number
#  It is essential not to have spaces around the = sign in the previous commands
echo "FlashFolder:" ${FlashFolder}
echo "Antenna:" ${AntennaFieldsDir}
