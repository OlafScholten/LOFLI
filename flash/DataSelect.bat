@echo off
Set LL_BaseDir="C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\LOFLI"
call %LL_BaseDir%\ShortCuts.bat
set Prog=DataSelect

:: del  %LL_bin%%Prog%.exe

call %LL_scripts%\RunProgram.bat DataSelect.in
:: pause
exit
