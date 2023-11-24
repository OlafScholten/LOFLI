#  Define some importants shortcuts for running the LOFLI codes (LL=Lofar Lightning)
#@echo off is in LINUX:  scriptname > /dev/null

# Test if the LOFLI system was installed already
if [[ ! -z "${LL_mod}" ]]; then
   echo LOFLI system variables have already been defined
else
   # or using a short-hand version    [[ -z "${DEPLOY_ENV}" ]] && MyVar='default' || MyVar="${DEPLOY_ENV}"
   echo Define LOFLI system variables and PATH extensions
  
   # Main directory of the LOFLI installation
   export LL_BaseDir="LL_Base"   
   
   ###### User installation dependent settings ###############################
   # Preambule for the files "LBA_Vout_theta.txt" and "LBA_Vout_phi.txt" containing the antenna function as a table
   export AntennaFun="${LL_Base}/AntenFunct/v2-"
   
   # Directory with the coordinates of the LOFAR antennas (needed by the RFI-Mitigation program)
   export AntennaFieldsDir="${LL_Base}/AntennaFields"
   
   # Base directory from where the raw data are retrieved, needed when running "NewFlash" 
   export ArchiveBase="/home/olaf/kaptdata/lightning_data"
   
   # should point to the "fftpack5.1d" library, containing FFT routines
   FFTLIB="-lm /home/olaf/NumLib/bin/libfftpack5.1d.a"  # FFT library, double precision
   
   # Link to "lapack" linear algebra library; to find, use command     find / -xdev -name *lapack*
   LAPACKlib="/usr/lib/x86_64-linux-gnu/lapack/liblapack.a"
   
   # Link to "blas" library for basic algebra functions
   BLASlib="/usr/lib/x86_64-linux-gnu/blas/libblas.a"
   
   # point to directory where HDF5 utilities are stored
   export HDF5lib="/usr/lib/x86_64-linux-gnu/hdf5/serial"

   #############################################################################
   # Rest are LOFLI system definitions and should not be touched
   export  LL_WrkDir=$(pwd)

   export LL_bin=${LL_Base}/bin
   export LL_src=${LL_Base}/FORTRANsrc
   export LL_mod=${LL_src}/modules
   export LL_scripts=${LL_Base}/scripts
   export LL_utilities=${LL_Base}/GLEsrc
   export PATH=${PATH}:${LL_bin}  
   #echo "path:" ${PATH}
   
   #   Windows:
   # set LL_generalLib=C:\OlafsUtil\NumLib\bin
   # set LAPACKlib=%LL_generalLib%\liblapack.a
   # set BLASlib=%LL_generalLib%\libblas.a
   # set FFTLIB=%LL_generalLib%\libFFTPack-d.a
   #   LINUX:
   export LOFLIlib="-lm ${LL_bin}/libLOFLI.a ${FFTLIB} ${LAPACKlib} ${BLASlib}"
   export LOFLIinc="-I${LL_mod}"
   
   export FCFLAGS="-fcheck=bounds"
   export HDF5Compile="-g -O2 -fdebug-prefix-map=/build/hdf5-X9JKIg/hdf5-1.10.0-patch1+docs=. -fstack-protector-strong -I/usr/include/hdf5/serial"
   export HDF5Link="-L${HDF5lib} ${HDF5lib}/libhdf5hl_fortran.a ${HDF5lib}/libhdf5_hl.a ${HDF5lib}/libhdf5_fortran.a ${HDF5lib}/libhdf5.a -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,${HDF5lib}"
   
   
   #REM  Old definitions
   export ProgramDir="${LL_bin}/"
   export UtilDir="${LL_utilities}/"
   #:: echo %Path%
   Vrsn=""  # version number
   #  It is essential not to have spaces around the = sign in the previous commands
   echo "Antenna:" ${AntennaFieldsDir}

fi
export FlashFolder=$(pwd)
echo "FlashFolder:" ${FlashFolder}


