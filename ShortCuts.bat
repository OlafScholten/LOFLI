rem Define some importants shortcuts for running the LOFLI codes (LL=Lofar Lightning)
@echo off
:: setlocal enableDelayedExpansion

REM Test if the LOFLI system was installed already
If NOT [%LL_mod%]==[]  (
   Echo LOFLI system variables have already been defined
   GOTO :EOF
)

echo Define LOFLI system variables and PATH extensions
:: ##### User installation dependent settings #################################
::	lines starting with :: are just comment lines
:: Main directory of the LOFLI installation
:: Enter in the windows 'environement variables' the following:
::  variable           LL_Base
::  pointing to        "C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\LOFLI"
::  The quotes are needed is the base path contains spaces
:: Set LL_BaseDir="C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\LOFLI"\
Set LL_BaseDir=%LL_Base%\

:: Directory containing the antenna function as a table
:: Apparently in windows  \"  unfolds to only "  while "\ stays
set AntenFun="C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\Imaging\LMA\LMA2019\AntenFunct"

:: Directory with the coordinates of the LOFAR antennas (needed by the RFI-Mitigation program)
set AntennaFieldsDir=%LL_Base%\AntennaFields

:: Base directory where the raw data are stored, needed when running "NewFlash"
set ArchiveBase=%LL_Base%\lightning_data

:: should point to the "fftpack5.1d" library, contains FFT routines, double precision
set LL_generalLib=C:\OlafsUtil\NumLib\bin
set FFTLIB=%LL_generalLib%\libFFTPack-d.a

:: Link to "lapack" linear algebra library; to find, use command     find / -xdev -name *lapack*
set LAPACKlib=%LL_generalLib%\liblapack.a

:: Link to "blas" library for basic algebra functions
set BLASlib=%LL_generalLib%\libblas.a

:: point to directory where HDF5 utilities are stored (not used in implemented windows version)
set HDF5lib="C:\Program Files\HDF_Group\HDF5\1.10.5\lib"

:: ############################################################################
:: Rest are LOFLI system definitions and should not be touched

set LL_bin=%LL_Base%\bin
set LL_src=%LL_Base%\FORTRANsrc
set LL_mod=%LL_src%\modules
set LL_scripts=%LL_Base%\scripts
set LL_utilities=%LL_Base%\GLEsrc
set PATH=%PATH%;%LL_bin%

REM  Not implemented in windows version
set HDF5Compile="-g -O2 -fdebug-prefix-map=/build/hdf5-X9JKIg/hdf5-1.10.0-patch1+docs=. -fstack-protector-strong -I/usr/include/hdf5/serial"
set HDF5Link="-L${HDF5lib} ${HDF5lib}/libhdf5hl_fortran.a ${HDF5lib}/libhdf5_hl.a ${HDF5lib}/libhdf5_fortran.a ${HDF5lib}/libhdf5.a -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,${HDF5lib}"

set LOFLIinc=-I%LL_mod%
set LOFLIlib=-lm %LL_bin%\libLOFLI.a %FFTLIB% %LAPACKlib% %BLASlib%
set FCFLAGS="-fcheck=bounds"

REM  Old definitions
set ProgramDir=%LL_src%\
set UtilDir=%LL_Utilities%\
set "Vrsn=-v21"
:: echo %Path%

set FlashFolder=%cd%
set LL_WrkDir=%cd%
Echo "FlashFolder:" %FlashFolder%

:EOF
