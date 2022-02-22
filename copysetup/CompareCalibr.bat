:: = rem   -lm C:\OlafsUtil\NumLib\bin\libFFTPack-d.a
call ..\ShortCuts.bat

  cd %ProgramDir%
  gfortran -o ComCal.exe CompareCalibrations.f90
  cd %FlashFolder%

%ProgramDir%ComCal  <CompareCalibr.in
:: pause
:: pause
exit
