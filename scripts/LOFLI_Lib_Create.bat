@echo off
setlocal enableDelayedExpansion

:: Set LL_BaseDir="C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\LOFLI\"
:: call %LL_BaseDir%ShortCuts.bat

echo fortransrc %LL_src%
cd %LL_BaseDir%
if exist bin\f90split.exe (
  echo Yes, %LL_bin%\f90split.exe exists
) else (
   cd %LL_bin%
   gfortran -c %LL_src%\f90split_OS.f90
   gfortran -o f90split f90split_OS.o
  echo Just created %LL_bin%\f90split.exe
)

cd %LL_src%
if exist modules\ (
  echo Yes, %LL_mod% exists
) else (
  mkdir modules
  echo Just created directory %LL_mod%
)

cd modules
REM  The order matters in this list because of loading the modules at the right time.
for /F %%x in (..\..\scripts\OrderedLibraryFiles.lst) do (
  set FILENAME=..\%%x
  REM echo ===========================  Splitting !FILENAME! ===========================
  REM set BareFILENAME=%%~nx
  REM echo Bare file name is: !BareFILENAME!
  %LL_bin%\f90split.exe "!FILENAME!"
  REM pause
   for /F %%y in ('dir /B/D *.f90') do (
      REM echo compiling file %%y as extracted from !FILENAME!
      gfortran %FCFLAGS% -c %%y
      REM pause
   )
   ar rcs libLOFLI.a *.o
   del *.o
   del *.f90
)
copy libLOFLI.a %LL_bin%\libLOFLI.a
::pause
