:: = rem
::cd ..\
gfortran -o AntFuncExe.exe AntFunc-v3.f90
::cd 2018
:: pause
:: "C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\Lght_papers\Brian\leaves\Pictures\rewriteEvent"  <rewriteEvent.in >rewriteEvent.out
::..\AntFuncExe  <AntFunc.in >AntFunc.out
AntFuncExe
 pause
exit
 RunGLEplots.bat
 pause
exit
echo off
set EHTO=C:/Users/Olaf/Documents/AstroPhys/GeoMagn/REP/xxxx.dat
set main=C:\Users\Olaf\Documents\AstroPhys\GeoMagn
echo on
