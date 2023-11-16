Set LL_BaseDir="C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\LOFLI"
call %LL_BaseDir%\ShortCuts.bat
Set LL_BaseDir="C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\LOFLI"
set LL_Utilities=%LL_BaseDir%\GLEsrc\
set UtilDir=%LL_Utilities%

@echo on
set glefile="C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\LOFLI\GLEsrc\TrackScatt.gle"
set gleinpt="C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\\LOFLI\flash\files\flash-all"
echo LL_Utilities: %LL_Utilities%
echo UtilDir: %UtilDir%
echo FlashFolder: %FlashFolder%
call GLE -d jpg -r 300 -o TrSc_flash-all.jpg %glefile% %gleinpt%
pause
call GLE -d jpg -r 300 -o TrSc_flash-all.jpg %UtilDir%\TrackScatt.gle "%FlashFolder%\files/flash-all"
pause
call GLE -d jpg -r 300 -o flash-allAmplFit.jpg %UtilDir%Intensity.gle "%FlashFolder%\files/AmplFitflash-allMx_d"
call GLE -d jpg -r 300 -o Imp_flash-all.jpg %UtilDir%SourcesPlot.gle "%FlashFolder%\files/flash-all"
pause
