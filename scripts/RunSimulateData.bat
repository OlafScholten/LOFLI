@echo off
Set LL_BaseDir="C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\LOFLI"
call %LL_BaseDir%\ShortCuts.bat
set prog=SimulateData

cd %LL_bin%
if NOT exist %prog%.exe (
  echo Create %prog%.exe
  gfortran  %FCFLAGS% -o %prog%.exe %LL_src%\%prog%.f90 %LOFLIinc%  %LOFLIlib%
)

cd %LL_WrkDir%
echo working directory: %cd%
echo Starting: %prog% with input file %1
@echo on
 pause
%prog%.exe  <%1
