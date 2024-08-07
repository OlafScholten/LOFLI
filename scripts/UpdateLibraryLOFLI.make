# Update Library -- make file
# Should check if a source-file is more recent than the library file and update the library if this is the case.
# To update the library with new code, the sourcefile is split and stored in a temporary folder and each is component compiled.
#   - the .o are replaced in the archive.
#   - the .mod are placed at ...???
# command to run:   make -f UpdateLibraryLOFLI.make
# run in Fortran-source directory

#LibrarySources = ConstantsModules.f90 
#LibrarySources = ConstantsModules.f90 FFT_routines.f90 ParamModules.f90 AntFuncCnst.f90 AntFunct.f90 InterferomPars.f90 \
# PredefTrack.f90 MappingUtilitiesModules.f90 CalibrationsMod.f90 LOFLI_InputHandling.f90 HDF5_LOFAR_Read.f90 MappingUtilities.f90 GLEplotUtil.f90 \
# Ant-Read.f90 CrossCorr.f90 FindCallibr.f90 MGMR3D_spline.f90 FitParams.f90 nl2sol.f90 Fitter_CorrMax.f90 FindSources.f90 \
# SourceTryal.f90 ExplorationOption.f90 CurtainPlotOption.f90 \
# ImpulsImagOption.f90 InterferometryOptSbRtns.f90 InterferometryOption.f90 EIOption.f90 EICallRtns.f90 EIFitter.f90 \
# ECallibrOption.f90 System_Utilities.f90 PeakInterferoOption.f90
#

$(info Updating LOFLI Library)

#all: libLOFLI.a

${LL_bin}/libLOFLI.a: ${LibrarySources}
	$(info 'Complete list:' $(LibrarySources))
	$(info 'Updates used:' $?)
	-Messages=$(shell ${LL_scripts}/LibraryAdd.sh '$?')   # this is where the actual work is done
	$(info Messages: $(Messages)) # still creates error messages

.PHONY: clean
clean:
	rm -f *.o
