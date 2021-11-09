rem my FFTPACK library creation

cd SRCd
C:\OlafsUtil\NumLib\F-split\f90split   ..\fftpack5.1d.f90
gfortran -c *.f90
ar rcs libFFTPack.a *.o
copy libFFTPack.a ..\..\bin\libFFTPack-d.a

pause
exit
