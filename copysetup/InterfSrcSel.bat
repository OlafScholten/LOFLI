:: = rem   -lm C:\OlafsUtil\NumLib\bin\libFFTPack-d.a
call ..\ShortCuts.bat

  cd %ProgramDir%
  gfortran -o InterfSrcSelect.exe InterfSrcSelect%Vrsn%.f90
  cd %FlashFolder%

%ProgramDir%InterfSrcSelect  <InterfSrcSel.in
:: pause
call RunGLEplots.bat
:: pause
exit
