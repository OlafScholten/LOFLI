Set LL_BaseDir="C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\LOFLI"
call %LL_BaseDir%\ShortCuts.bat
::@echo on
set input=%1
echo program run is: %Prog% with argument %input%

if NOT exist %LL_bin%libLOFLI.a (
  echo Create %LL_bin%libLOFLI.a
  call %LL_scripts%LOFLI_Lib_Create.bat
  cd %LL_bin%
  del %Prog%.exe
)


cd %LL_bin%
set R_src=..\FORTRANsrc
set Rinc=-I%R_src%\modules
set Rlib=-lm libLOFLI.a  %FFTLIB%  %LAPACKlib%  %BLASlib%

if NOT exist %prog%.exe (
  echo Create %prog%.exe
::  gfortran  %FCFLAGS% -o %Prog%.exe %R_src%\%Prog%.f90 %Rinc%  %Rlib%
  gfortran  %FCFLAGS% -o %Prog%.exe %LL_src%%Prog%.f90 %LOFLIinc%  %LOFLIlib%
)

cd %LL_WrkDir%
echo executing %LL_bin%%Prog%.exe  "<"%input%
%LL_bin%%Prog%.exe  <%input%
