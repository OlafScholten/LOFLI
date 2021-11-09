:: = rem
::cd ..\
gfortran -o AntTst.exe AntennaTest.f90 -lm C:\OlafsUtil\NumLib\bin\libFFTPack-d.a
::cd 2018
:: pause
:: "C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\Lght_papers\Brian\leaves\Pictures\rewriteEvent"  <rewriteEvent.in >rewriteEvent.out
::..\AntFuncExe  <AntFunc.in >AntFunc.out
AntTst
 pause
exit
 RunGLEplots.bat
 pause
exit
echo off
set EHTO=C:/Users/Olaf/Documents/AstroPhys/GeoMagn/REP/xxxx.dat
set main=C:\Users\Olaf\Documents\AstroPhys\GeoMagn
echo on
