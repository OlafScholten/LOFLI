:: = rem
call ..\ShortCuts.bat

  cd %ProgramDir%
::  gfortran -o ShowCal.exe ShowCalibrations.f90  %LIBRARY%
  gfortran -o ShowCal.exe ShowCalibrations.f90
  cd %FlashFolder%

%ProgramDir%ShowCal.exe  <ShowCal.in
:: pause
call RunGLEplots.bat
:: pause
exit
echo off
set EHTO=C:/Users/Olaf/Documents/AstroPhys/GeoMagn/REP/xxxx.dat
set main=C:\Users\Olaf\Documents\AstroPhys\GeoMagn
echo on "%FlashFolder%"
