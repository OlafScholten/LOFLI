:: = rem
call ..\ShortCuts.bat

  cd %ProgramDir%
  gfortran -o Simulate.exe SimulateData.f90 %LIBRARY%
  cd %FlashFolder%

%ProgramDir%Simulate.exe  <Simulate.in
:: call RunGLEplots.bat
 pause
exit
echo off
set EHTO=C:/Users/Olaf/Documents/AstroPhys/GeoMagn/REP/xxxx.dat
set main=C:\Users\Olaf\Documents\AstroPhys\GeoMagn
echo on "%FlashFolder%"
