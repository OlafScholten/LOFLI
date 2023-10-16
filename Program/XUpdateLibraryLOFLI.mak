# Update Library -- make file
# Should check if a source-file is more recent than the library file and update the library if this is the case.
# To update the library with new code, the sourcefile is split and stored in a temporary folder and each is component compiled.
#   - the .o are replaced in the archive.
#   - the .mod are placed at ...???

#LibrarySources = ConstantsModules.f90 
LibrarySources = ConstantsModules.f90 FFT_routines.f90 ParamModules.f90 AntFuncCnst.f90 AntFunct.f90 InterferomPars.f90 \
 MappingUtilities.f90 
#
# FileA := fftpack5.1d.f90
# LibraryUpdate := ~/NumLib/bin/libfftpack5.1d.a
#UpdateLib = LibraryAdd.sh

$(info 'I am here' now)
#echo 'I am here'

all: libLOFLI.a

libLOFLI.a: ${LibrarySources}
#	UpdateFilesf90 := $(subst .f90,., $?)
#	$(info 'param1' $@)
	$(info 'Updates used:' $?)
#	$(info 'param3' $(UpdateFilesf90))
	$(info 'Complete list:' $(LibrarySources))
#	Ux := $(subst r,x, $(LibrarySources))
#	Ux := $(subst r,x,random terrr $(LibrarySources))
#	$(info 'param5' $(Ux))
#	@echo 'param5- $(Ux)'
#	$(UpdateLib) -f
	$(shell ./LibraryAdd.sh '$?')
#	$(shell ./LibraryAdd.sh $(UpdateFilesf90))   # does not work

