@echo off

:: del  %LL_Base%\bin\DataSelect.exe

call %LL_Base%\scripts\RunProgram.bat DataSelect DataSelect.in
:: pause
exit
