call %LL_Base%\ShortCuts.bat
::@echo on
set Prog=%1
set input=%2
echo program run is: %Prog% with argument %input%

if NOT exist %LL_bin%\libLOFLI.a (
  echo Create %LL_bin%\libLOFLI.a
  call %LL_scripts%\LOFLI_Lib_Create.bat
  cd %LL_bin%
  del %Prog%.exe
)


cd %LL_bin%

if NOT exist %prog%.exe (
  echo Create %prog%.exe
  gfortran  %FCFLAGS% -o %Prog%.exe %LL_src%\%Prog%.f90 %LOFLIinc%  %LOFLIlib%
  del *.mod
)

cd %LL_WrkDir%
if exist %input% (
   echo "File %input% exists"
   echo executing %LL_bin%\%Prog%.exe  "<"%input%
   %LL_bin%\%Prog%.exe  <%input%
) else (
   echo "%input% is not an existing file, assumed to be direct input"
   echo executing %LL_bin%\%Prog%.exe  %input%
   %LL_bin%\%Prog%.exe  %input%
)
