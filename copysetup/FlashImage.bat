:: = rem
call ..\ShortCuts.bat

  cd %ProgramDir%
  gfortran -o TrackExe.exe Track%Vrsn%.f90 %LIBRARY%
  cd %FlashFolder%

%ProgramDir%TrackExe  <FlashImage.in
call RunGLEplots.bat
:: pause
exit
echo off
set EHTO=C:/Users/Olaf/Documents/AstroPhys/GeoMagn/REP/xxxx.dat
set main=C:\Users\Olaf\Documents\AstroPhys\GeoMagn
echo on "%FlashFolder%"
