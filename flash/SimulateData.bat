@echo off

call %LL_Base%\scripts\RunProgram.bat SimulateData Simulate.in

call %LL_Base%\scripts\RunProgram.bat LOFAR-Imag SimIntf.in

exit
#sshfs scholten@kapteyn.astro.rug.nl:/net/dataserver3/data/users/hare/ ~/kaptdata
