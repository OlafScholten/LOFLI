Set LL_BaseDir="C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\LOFLI"
call %LL_BaseDir%\ShortCuts.bat
::@echo on
set input=%1
echo program run is: %Prog% with argument %input%

if NOT exist %LL_bin%libLOFLI.a (
  echo Create %LL_bin%libLOFLI.a
  call %LL_scripts%LOFLI_Lib_Create.bat
  cd %LL_bin%
  del %Prog%.exe
)


cd %LL_bin%

if NOT exist %prog%.exe (
  echo Create %prog%.exe  -cpp  -Dvariable=value
  gfortran  %FCFLAGS%  -o %Prog%.exe %LL_src%%Prog%.F90 %LOFLIinc%  %LOFLIlib%
  del *.mod
)

cd %LL_WrkDir%
echo executing %LL_bin%%Prog%.exe  "<"%input%
%LL_bin%%Prog%.exe  <%input%
