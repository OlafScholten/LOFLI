FC = gfortran
#FCFLAGS = -ggdb -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid 
FCFLAGS =

EXECUTE = ttret
PROGRAM = $(EXECUTE)-v30f
DEPENDENCIES1 = FitFoot.f90 Fit_RadioFoot.f90 TTret_shower-v31f.f90 TTret_analyse-v30f.f90 SetParams.f90 RFootPars.f90
DEPENDENCIES2 = ../TTret_subr-v30.f90   ~/NumLib/LSQ/nl2sol.f90 
LIBRARY = ~/NumLib/bin/libfftpack5.1d.a

all: $(EXECUTE)
OBJECTS = $(PROGRAM).o
SOURCES = $(PROGRAM).f90

$(EXECUTE): $(OBJECTS)
	$(FC) $(FCFLAGS) -o $(EXECUTE) $(OBJECTS) -lm $(LIBRARY)
$(OBJECTS): $(SOURCES) $(DEPENDENCIES1) $(DEPENDENCIES2)
	$(FC) $(FCFLAGS) -c $(SOURCES)
clean:
	rm  .mod

# @gfortran ttret-v30f.f90 -o FIT_T_Tret.exe -lm C:\OlafsUtil\NumLib\bin\libFFTPack-d.a
# rem @gfortran -ggdb -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid ttret-v18.f90 -o T_Tret.exe -lm C:\OlafsUtil\NumLib\bin\libFFTPack-d.a
 
