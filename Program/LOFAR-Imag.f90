 !----------------
 !   Main program
 !-------------------------------------------
 !-------------------------------------------
 !------- Source file:LOFAR-Imag.f90 ----
 !   Include 'ConstantsModules.f90'
 !   I   Module constants
 !   I   Module DataConstants
 !   I   M   Character(len=100), save :: ProgramFolder='' !  path to the directory containing the fortran programs
 !   Include 'FFT_routines.f90'
 !   I   Module FFT
 !   I   M   Subroutine RFTransform_su(N_T)
 !   I   M   Subroutine DAssignFFT()
 !   I   M   Subroutine RFTransform_CF(A,Cnu)
 !   I   M   Subroutine RFTransform_CF_Filt(A,F,t_shft,Cnu)
 !   I   M   Subroutine RFTransform_CF2CT(Cnu,CD)
 !   I   M   Subroutine R2C(R,C)
 !   I   M   Subroutine  C2R(R,C)
 !   I   M   Subroutine DownSamlple(E_t,E,padding,tTrace_dim,FF_dim, filt, E_nu_dwn, inui,inum)
 !   Include 'ParamModules.f90'  ! v18d: Cnu storage changed  for InterfEngineB
 !   I   Module Chunk_AntInfo
 !   I   M   Subroutine Alloc_Chunk_AntInfo
 !   I   Module ThisSource
 !   I   M   Subroutine Alloc_ThisSource
 !   I   Module FitParams
 !   I   Module Explore_Pars
 !   I   M   Subroutine Alloc_Explore_Pars
 !   I   Module Interferom_Pars
 !   I   M   Subroutine Alloc_SInterferom_Pars
 !   I   M   Subroutine Alloc_EInterferom_Pars
 !   I   M   Subroutine DeAlloc_SInterferom_Pars
 !   I   M   Subroutine DeAlloc_EInterferom_Pars
 !   Include 'MappingUtilities.f90'
 !   I   Include 'HDF5_LOFAR_Read.f90'
 !   I   I   Module HDF5_LOFAR_Read
 !   I   I   M   Subroutine GetFileName(filenr, nxx)
 !   I   I   M   Subroutine ListGroups
 !   I   I   M   Subroutine CloseFile
 !   I   I   M   Subroutine ListGroupStructure(GroupName)
 !   I   I   M   Subroutine CloseGroup
 !   I   I   M   Subroutine ListDataAtt(GroupName,DSetName, prnt)
 !   I   I   M   Subroutine CloseDSet
 !   I   I   M   Subroutine GetData(Chunk, DSet_offset, DSet_dim)
 !   I   I   M   Subroutine AttrRead(file_id,Obj_id,Obj_name,prn)
 !   I   I   M   Subroutine DataRead(dset_id, Chunk, DSet_offset, DSet_dim)
 !   I   I   M   Subroutine GetDataChunk(GroupName,DSetName, Chunk, DSet_offset, DSet_dim, prnt, DataReadErr)
 !   I   I   M   Subroutine CloseDataFiles()
 !   I   Module StationMnemonics
 !   I   M   Subroutine Station_ID2Mnem(STATION_ID,Station_Mnem)
 !   I   M   Character(len=5) Function Statn_ID2Mnem(STATION_ID)
 !   I   M   Subroutine Station_Mnem2ID(Station_Mnem,STATION_ID)
 !   I   M   Integer Function Statn_Mnem2ID(Station_Mnem)
 !   I   Subroutine Station_ID2Calib(STATION_ID,Ant_ID,StatAnt_Calib)
 !   I   Subroutine ITRF2LOFARConstruct()
 !   I   Subroutine ITRF2LOFAR(ITRF,LOFAR)
 !   I   Subroutine RelDist(Source,LFRAnt,RDist)
 !   I   Real(kind=8) Function SubRelDist(SrcPos,i_ant,i_chunk)
 !   I   Module unque
 !   I   M   Subroutine unique(vec,vec_unique)
 !   I   M   Subroutine Selection_sort(a)
 !   I   M   Subroutine Double_sort(a)
 !   I   M   Subroutine Double_IR_sort(N,a,R)
 !   I   M   Subroutine Double_RI_sort(N,a,R)
 !   I   M   Subroutine sort(n, a)
 !   I   Subroutine GetNonZeroLine(lineTXT)
 !   I   Module ansi_colors
 !   I   M   Function color(str, code) result(out)
 !   I   Subroutine Convert2m(CenLoc)
 !   Include 'GLEplotUtil.f90'
 !   I   Module GLEplots
 !   I   M   use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows
 !   I   M   Subroutine GLEplotControl(PlotType, PlotName, PlotDataFile, SpecialCmnd, Submit)
 !   Include 'AntFunct.f90'
 !   I   Module AntFunCconst
 !   I   Subroutine Cheb_odd(Phi,Cheb)
 !   I   Subroutine Legendre_odd(x,Lgdr)
 !   I   Subroutine AntFun(thet_d,phi_d)
 !   I   Subroutine AntFun_Inv(thet_d,phi_d)
 !   I   Subroutine AntFieParGen()
 !   I   WRITE (2,*) "Using environment variable: AntennaFun; Antenna Function read from files:",' "',TRIM(AntFunFile)//'LBA_Vout_*.txt"
 !   Include 'Ant-Read.f90'
 !   I   Subroutine AntennaRead(i_chunk,SourceGuess)
 !   I   Subroutine PlotSpectra(i_chunk)
 !   I   Subroutine GLE_script(unt, file, Ant_nr, i_chunk)
 !   I   Subroutine SimulationRead(SourceGuess)
 !   Include 'CrossCorr.f90'
 !   I   Subroutine BuildCC(StatMax,DistMax)
 !   I   Subroutine GetCorrSingAnt( i_ant, J_Corr, i_eo, i_chunk)
 !   I   Real(kind = 8) Function CCorr_der(t_max,j_corr,i_Peak)
 !   I   Subroutine CrossCorrPol(CnuRt, i_ant, i_chunk, StLoc, T2_dim, i_peak, t_shft, CrCort, CrCorp)
 !   I   Subroutine CrossCorr(CnuR, i_ant, i_chunk, StLoc, T2_dim, t_shft, CrCor)
 !   I   Subroutine GetAntPol(Thet_d, Phi_d, i_ant, i_chunk, StLoc, T_dim, T2_dim, PolBasis, alpha, W,CnuPol)
 !   I   Subroutine GetAnt(i_ant, i_chunk, StLoc, T_dim, T2_dim, Cnu)
 !   I   Subroutine GetAntSourceAng(i_ant, i_chunk, SourceLoc, Thet_d, Phi_d, Dist)
 !   I   Subroutine CrossCorr_Max(i_ant,i_chunk, i_Peak, CrCor, TrueCC, RtMax, CCval, CCPhase, SearchRange, Error)
 !   I   Subroutine ReImAtMax(i_ant,i_chunk, i_Peak, CrCor, ACCorr, RCCorr, RtMax, Aval, Rval, Ival, SearchRange, Error)
 !   I   Subroutine GetRefAnt
 !   Include 'FindCallibr.f90'   ! Station Callibration
 !   I   Subroutine FindCallibr(SourceGuess)
 !   I   Subroutine ReadPeakInfo(ReadErr)
 !   I   Subroutine FitCycle(FitFirst,StatMax,DistMax,FitNoSources)
 !   I   Subroutine GetLargePeaks(i_ant, i_chunk, Peakposs)
 !   I   Subroutine MergeFine
 !   Include 'FitParams.f90'
 !   I   Include 'MGMR3D_spline.f90'  ! ../../NumLib
 !   I   Subroutine SetFitParamsLMA(X,first,FitPos)
 !   I   Subroutine PrntFitPars(X)
 !   I   Subroutine PrntNewSources()
 !   I   Subroutine PrntCompactSource(i_Peak,i_eo)
 !   I   Subroutine X2Source(X)
 !   I   Subroutine Find_unique_StatAnt()
 !   I   Integer Function XIndx(i,i_peak)
 !   Include 'Fitter_CorrMax.f90'     ! uses chi-square mini
 !   I   Include 'nl2sol.f90'
 !   I   I   write ( pu, '(a)' ) 'Function evaluation limit.'
 !   I   Subroutine FitCCorr(X)
 !   I   Subroutine JacobianCorrTime ( meqn, nvar, X_p, nf, Jacobian, uiparm, urparm, ufparm )
 !   I   Subroutine CompareCorrTime ( meqn, nvar, X_p, nf, R, uiparm, Jacobian, ufparm )
 !   I   Subroutine ufparm ( meqn, nvar, x )
 !   Include 'FindSources.f90'
 !   I   Subroutine SourceFind(TimeFrame,SourceGuess,units)
 !   I   Subroutine SourceFitCycle(StatMax,DistMax)
 !   I   Subroutine CleanPeak(i_eo_in,PeakPos,SourcePos,Wl,Wu)
 !   I   Subroutine DualPeakFind(PeakS_dim, i_chunk, PeakD_nr, PeakSP, PeakSWl, PeakSWu, PeakSAmp)
 !   Include 'SourceTryal.f90'  ! v16.f90' = v17.f90'  ; v17a.f90' uses grid search
 !   Include 'ExplorationOption.f90'
 !   I   Subroutine ExplorationRun
 !   I   Subroutine GetSpectStats(i_chunk, i_expl)
 !   I   Subroutine AnalSpectStats
 !   Include 'CurtainPlotOption.f90'
 !   I   Subroutine PlotAllCurtainSpectra(CurtainWidth)
 !   I   Subroutine GLEscript_CurtainPlot(unt, file, CurtainWidth, i_Peak)
 !   I   Subroutine GLE_Corr()
 !   Include 'InterferometryOption.f90'  ! d: Cnu storage changed for InterfEngineB
 !   I   Include 'InterferometryOptSbRtns.f90'
 !   I   I   Subroutine Pol2Carth(LocPol,LocCth)
 !   I   I   Subroutine Carth2Pol(LocCth,LocPol)
 !   I   I   Subroutine PixBoundingBox()
 !   I   I   Subroutine OutputIntfPowrTotal(RefAntSAI, i_eo)
 !   I   I   Subroutine OutputIntfPowrMxPos(i_eo)
 !   I   I   Real*8 Function Sq(x)
 !   I   I   Subroutine FindInterpolMx(i,d_gr,SMPow,Qualty)
 !   I   I   Function Paraboloid(y0,Ay,By,Ry,d_gr) result(SMPow)
 !   I   I   Subroutine FindBarycenter(i,Baryd,SMPow,Qualty)
 !   I   I   Subroutine OutputIntfSlices(i_eo)
 !   I   Subroutine InterferometerRun
 !   I   'cp ${ProgramDir}LOFAR-Imag-v20 ./Imag-'//TRIM(OutFileLabel)//'.exe',  &  ! Just to display a more intelligent name when
 !   Include 'EIOption.f90'
 !   I   Subroutine EI_Run
 !   I   Subroutine Inverse33(A,B, eig, PolBasis)
 !   I   Subroutine EigvecS33(A,E,ev)
 !   I   Subroutine EISelectAntennas()
 !   I   Subroutine EIPixBoundingBox()
 !   I   Subroutine EISetupSpec(Nr_IntFer, IntfNuDim, CMCnu)
 !   I   Subroutine EIEngine(Nr_IntFer, IntfNuDim, CMCnu, CMTime_pix)
 !   I   Subroutine EIAnalyzePixelTTrace(i_N, i_E, i_h, SumWindw, IntfNuDim, CMTime_pix)
 !   Include 'EIStokes-testW.f90'
 !   I   Subroutine EI_PolarizW(Nr_IntFer, IntfNuDim, i_slice)
 !   I   Subroutine matinv3(A,B)
 !   Include 'SIOption.f90'
 !   I   Subroutine SI_Run
 !   I   Subroutine InterfEngineB(i_eo, i_N, i_E, i_h, Nr_IntFer)
 !   I   Subroutine SISelectAntennas(i_eo)
 !   I   Subroutine SISetupSpec(i_eo)
 !   I   Subroutine SIAnalyzePixelTTrace(i_eo, i_N, i_E, i_h)
 !   Include 'System_Utilities.f90'
 !   I   Subroutine System_Initiation(Utility, release, ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows)
 !   I   Character(len=*), intent(out) :: ProgramFolder, UtilitiesFolder, FlashFolder, FlashName
 !   I   ProgramFolder='%ProgramDir%'
 !   I   ProgramFolder='${ProgramDir}'
 !   I   Subroutine System_MemUsage(valueRSS)
 !   I   end Subroutine System_MemUsage
 !   Program LOFAR_Imaging
 !   use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, DataFolder, FlashName, Windows
 !   Call System_Initiation(Utility, release, ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows)  ! set Flash name & f
 !   Call GLEplotControl(SpecialCmnd='cd '//TRIM(ProgramFolder))
 !   Call GLEplotControl(SpecialCmnd=TRIM(ProgramFolder)//'TrackExe.exe  <'//TRIM(lname))
 !   Subroutine AntennaSanity()
 !-------------------------------------------
 !------------------------------------------
    Include 'ConstantsModules.f90'
    Include 'FFT_routines.f90'
    Include 'ParamModules.f90'  ! v18d: Cnu storage changed  for InterfEngineB
    Include 'HDF5_LOFAR_Read.f90'
    Include 'MappingUtilities.f90'
    Include 'GLEplotUtil.f90'
    Include 'AntFunct.f90'
    Include 'Ant-Read.f90'
    Include 'CrossCorr.f90'
    Include 'FindCallibr.f90'   ! Station Callibration
    Include 'FitParams.f90'
    Include 'Fitter_CorrMax.f90'     ! uses chi-square mini
    !Include 'StationPolar-v10.f90'
    Include 'FindSources.f90'
    Include 'SourceTryal.f90'  ! v16.f90' = v17.f90'  ; v17a.f90' uses grid search
    Include 'ExplorationOption.f90'
    Include 'CurtainPlotOption.f90'
    Include 'InterferometryOption.f90'  ! d: Cnu storage changed for InterfEngineB
    Include 'EIOption.f90'
    !Include 'EIStokes-v21-test.f90'
    !Include 'EIStokes-v21.f90'
    Include 'EIStokes-testW.f90'
    Include 'SIOption.f90'
    Include 'System_Utilities.f90'
!-----------------------------------
Program LOFAR_Imaging
!
!   v4: allow for removing bad antennas and polarization flips
!   v5: allow for fitting peaks in different chunks of data
!       Loose compatability with image-quality fitting, allow for correlation-time fitting only
!   v6: (aug19) account for noise-error in the pulse of the reference antenna by a constant offset
!   v7: (15aug19) allow for fixing equal source position for even and odd antennas
!   v7: (15aug19) allow for fitting delays of individual antennas
!   v8: Include antenna function
!   v8: Average signals antennas of one station
!   v8: Move FFT-setup to higher up in the program
!   v8: PeakPos is the position of a peak in a virtual antenna at the center of CS002
!   v9: Separate read-in of spectra to separate subroutine
!   v9: Start with real imaging
!
!  v13: Include Kalman filter and taken out again. Several refinements in 'SourceFind'
!  v13: Relative distance of candidate pulses tuned
!  v13: Double sources are ordered in magnitude
!  v14: Make sure doble=dual works with zero distance between the two ref. antennas
!  v15: Increase value chunk-length
!  v15: Move explore to separate subroutine
!  v15: Include more sophisticated norm to background in Explore
!  v15: Better antenna position in RFI-suppression program
!  v15: Automatic Chunk_dim
!  v15: Modify CurtainPlot, make several
!  v15: Write pulse strength to file (width and amplitude)
!To do
!  Take flash label from antenna files
!  write to .sh   gle -d pdf GLE_file.gle
!  v16: include interferometry
!  v17: work on locating algorithm versions 17 and 17a
!  v17a: NL2SOL minimal stepsize improved
!  v18: Optimize locating algorithm and esstablish best of 17 and 17a
!  File Units in use: (2=output),10,11,(12=RFI filter),13,(14=H5 data),
!     (15,16,17=source locations),19,20,21,(25,26,27=peakfind stat),29,30,31
!  v18: Peakfinding optimized further
!  v18: Zero-data criterium as well as normalization modified (in 'AntRead')
!  v18: Produce plots of interferometric maxima for automatic-slices; units= 28, 29 used locally
!  v19: Plotting control unified
!
    use constants, only : dp,pi,ci,sample,Refrac
    use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, DataFolder, FlashName, Windows
    use DataConstants, only : Time_dim, Cnu_dim, Production, RunMode
    use DataConstants, only : PeakNr_dim, ChunkNr_dim, Diagnostics, EdgeOffset, Calibrations, OutFileLabel
    use DataConstants, only : Polariz
    use ThisSource, only : Alloc_ThisSource, Dual, RealCorrelation, t_ccorr, Safety, CurtainPlot
    use ThisSource, only : CCShapeCut_lim, ChiSq_lim, EffAntNr_lim, PeakPos
    use FitParams, only : FitIncremental, Fit_AntOffset, WriteCalib, FullAntFitPrn, ImagingRun, AntennaRange
    use FitParams, only : PeakS_dim, MaxFitAntDistcs, MaxFitAntD_nr, Explore, Sigma_AntT, SearchRangeFallOff
    use FitParams, only : FullSourceSearch, SigmaGuess
    use Chunk_AntInfo, only : ExcludedStatID, Start_time, BadAnt_nr, BadAnt_SAI, DataReadError, TimeFrame
    use Chunk_AntInfo, only : NoiseLevel, PeaksPerChunk, TimeBase, Simulation, WriteSimulation
    use Chunk_AntInfo, only : ExcludedStat_max, SgnFlp_nr, PolFlp_nr, SignFlp_SAI, PolFlp_SAI, BadAnt_SAI, AntennaNrError
    use Chunk_AntInfo, only : Alloc_Chunk_AntInfo, ExcludedStat   !  subroutine
    use FFT, only : RFTransform_su,DAssignFFT
    use StationMnemonics, only : Statn_ID2Mnem, Statn_Mnem2ID
    use Explore_Pars, only : NMin, NMax, Emin, EMax
    Use Interferom_Pars, only : IntfPhaseCheck, IntfSmoothWin, ChainRun, PixPowOpt
    use GLEplots, only : GLEplotControl
    Implicit none
    Character(len=20) :: Utility, release
    INTEGER :: DATE_T(8),i
    !CHARACTER*12 :: REAL_C(3)
    !Logical :: CurtainPlot=.false.
    !Character(len=5) :: ExcludedStat(ExcludedStat_max)
    character(len=2) :: chj
    !
    Integer :: j,i_chunk, ChunkNr_start, ChunkNr_stop, units(0:2), FitRange_Samples, CurtainWidth
    Real*8 :: StartTime_ms, StartingTime, StoppingTime, D
    Real*8 :: SourceGuess(3,10) ! = (/ 8280.01,  -15120.48,    2618.37 /)     ! 1=North, 2=East, 3=vertical(plumbline)
    CHARACTER(LEN=6) :: txt,Version
    CHARACTER(LEN=10) :: Sources
    Character(LEN=180) :: Sources1, Sources2
    Character(LEN=180) :: lname
    Character(LEN=250) :: TxtIdentifier, TxtImagingPars
    Character(LEN=250) :: Txt20Identifier, Txt20ImagingPars
    Logical :: XcorelationPlot, prnt=.false. ! CrossCorrelate=.false. ,
    Logical :: Interferometry=.false.
    Logical :: E_FieldsCalc=.true.
    Integer :: i_dist, i_guess, nxx, valueRSS
    NAMELIST /Parameters/ Explore, ImagingRun, Interferometry, FullSourceSearch, CurtainPlot, XcorelationPlot &
         , IntfPhaseCheck, IntfSmoothWin, TimeBase  &
         , Diagnostics, Dual, FitIncremental, Fit_AntOffset, RealCorrelation &
         , Simulation, WriteSimulation, SignFlp_SAI, PolFlp_SAI, BadAnt_SAI, PeakNr_dim, Calibrations, WriteCalib &
         , ExcludedStat, FitRange_Samples, FullAntFitPrn, AntennaRange, E_FieldsCalc, PixPowOpt, OutFileLabel, ChainRun &
         , CCShapeCut_lim, ChiSq_lim, EffAntNr_lim, Sigma_AntT, SearchRangeFallOff, NoiseLevel, PeaksPerChunk    !  ChunkNr_dim,
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Version='v21'
   release='v21 (Aug, 2021)'
   Utility='LOFAR_Lightn.Imaging'
   CALL get_environment_variable("LIBRARY", lname)
   WRITE (*,*) "LIBRARY=",TRIM(lname)
   CALL get_environment_variable("HOME", lname)
   WRITE (*,*) "HOME=",TRIM(lname)
   !
   ExcludedStat='  '
   ExcludedStatID=0
   SignFlp_SAI=-1
   BadAnt_SAI=-1
   PolFlp_SAI=-2
   XcorelationPlot=.false.
   FitRange_Samples=Safety
   CCShapeCut_lim=0.6 ; ChiSq_lim=80. ; EffAntNr_lim=0.8
   NoiseLevel=70.
   Simulation=""
   !reshape((/ -1.D-10,0.d0,0.d0,0.d0,-1.D-10,0.d0,0.d0,0.d0,-1.D-10 /), shape(array))
   !AntennaRange_km=AntennaRange
   read(*,NML = Parameters)
   !If(CurtainPlot) Then
   !   ChunkNr_dim=1
   !   PeakNr_dim=2
   !   Explore=.false.
   !   ImagingRun=.false.
   !   OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='CurtainPlots'//TRIM(OutFileLabel)//'.out')
   !ElseIf(Explore) Then
   If(Explore) Then
      Write(*,*) 'Explore, OutFileLabel=',TRIM(OutFileLabel)
      RunMode=1
      ImagingRun=.false.
      Production=.true.
      Interferometry=.false.
      !Production=.false.
      FitRange_Samples=70
      AntennaRange=2.49999999
      Dual=.false.
      ChunkNr_dim=1
      PeakNr_dim=10
      NoiseLevel=20.
      If(Simulation.ne."") then
         write(2,*) 'Explore is not compatible with running with a Simulation=', TRIM(Simulation)
         stop 'Incompatible options'
      EndIf
      OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='Explore'//TRIM(OutFileLabel)//'.out')
   ElseIf(ImagingRun) Then
      Write(*,*) 'ImagingRun, OutFileLabel=',TRIM(OutFileLabel)
      RunMode=3
      ChunkNr_dim=1
      PeakNr_dim=2
      FitIncremental=.true.
      Interferometry=.false.
      !FitRange_Samples=170
      FitRange_Samples=86 ! 90  !  some cases 70 is better, sometimes 90, 80 seems to be worse i.e. highly non-linear!
      OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='Imaging'//TRIM(OutFileLabel)//'.out')
   ElseIf(Interferometry) Then
      Write(*,*) 'Interferometry, OutFileLabel=',TRIM(OutFileLabel)
      RunMode=4
      FitRange_Samples=7
      Polariz=Dual  ! Affects reading in LOFAR Data; even/odd pairs only & equal delay for the pair
      ChunkNr_dim=1
      PeakNr_dim=2
      OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='Interferometer'//TRIM(OutFileLabel)//'.out')
   Else     ! Calibration run
      Write(*,*) 'Calibration, OutFileLabel=',TRIM(OutFileLabel)
      RunMode=2
      ChunkNr_dim=0
      TxtIdentifier="(i3,2i2,I8,3(F10.2,1x),F8.2)"
      Polariz=(Dual .and. RealCorrelation) .and. E_FieldsCalc
      If(Simulation.ne."") then
         write(2,*) 'Calibration is not compatible with running with a Simulation=', TRIM(Simulation)
         stop 'Incompatible options'
      EndIf
      Do i=1,10
         Call GetNonZeroLine(lname)
         Read(lname,*,iostat=nxx) StartTime_ms, SourceGuess(:,i)  ! just dummy arguments
         Call Convert2m(SourceGuess(:,i))
         If(nxx.ne.0) exit
         Read(lname,TxtIdentifier,iostat=nxx)  DATE_T(1:4),StartTime_ms, SourceGuess(:,i) ! just dummy arguments
         If(nxx.eq.0) exit
         ChunkNr_dim=i   ! this was a genuine chunk card
      EndDo
      If(ChunkNr_dim.ge.10) stop 'Chunk number too large'
      Rewind(unit=5) ! Standard input
      read(*,NML = Parameters) ! to reposition correctly
      OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='Calibrate'//TRIM(OutFileLabel)//'.out')
   Endif
   If(FitRange_Samples.lt.5) FitRange_Samples=5
   !AntennaRange=AntennaRange_km
   !
   !write(2,"(3x,5(1H-),1x,'LOFAR_Imaging version ',A20,25(1H-))") release
   !CALL DATE_AND_TIME (Values=DATE_T)
   !WRITE(2,"(3X,5(1H-),1x,'run on ',I2,'/',I2,'/',I4,' , started at ',&
   !    I2,':',I2,':',I2,'.',I3,1X,25(1H-))") &
   !    DATE_T(3),DATE_T(2),DATE_T(1),(DATE_T(i),i=5,8)
   !
   Call System_Initiation(Utility, release, ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows)  ! set Flash name & folder names
   !
   CALL DATE_AND_TIME (Values=DATE_T)
   Write(TxtIdentifier,"('Release ',A,', run on ',I2,'/',I2,'/',I4,' , at ',I2,':',I2,':',I2,'.',I3,', Flash: ',A )") &
       TRIM(release),DATE_T(3),DATE_T(2),DATE_T(1),(DATE_T(i),i=5,8),TRIM(FlashName)
    !
    write(2,NML = Parameters)
    !
    Safety=FitRange_Samples
    Call Alloc_Chunk_AntInfo
    Call Alloc_ThisSource
    Do i=-Safety,Safety
        t_ccorr(i)=i ! time axis needed for spline interpolation of CrossCorrelation
    Enddo
    !Do j=1,SgnFlp_nr_max
    !  SignFlp_SAI(j).lt.0) exit
    !Enddo
    write(2,*)'refractivity=',Refrac-1.
    SgnFlp_nr=count(SignFlp_SAI.gt.0)
    PolFlp_nr=count(PolFlp_SAI.gt.0)  !first pole-flip, then bad antenna
    BadAnt_nr=count(BadAnt_SAI.gt.0)
    If(any(mod(PolFlp_SAI,2) .eq. 1,1) ) stop 'PolFlp-Ant_nr should be even'
    !write(2,*) 'signflip=',SgnFlp_nr
    write(*,*) 'PolarizationFlips=',PolFlp_nr,'; BadAntennas=',BadAnt_nr ! , &
    !'Color test'//achar(27)//'[95m pink '//achar(27)//'[0m.'
!
    Do j=1,ExcludedStat_max
      If(ExcludedStat(j).eq.'     ') exit
      ExcludedStatID(j)=Statn_Mnem2ID(ExcludedStat(j))
      write(2,*) j,ExcludedStatID(j),ExcludedStat(j)
    enddo
    !
    !Open(unit=12,STATUS='unknown',ACTION='read',FORM ="unformatted", FILE = 'Book/RFI_Filters-v18.uft')
    !Open(unit=14,STATUS='old',ACTION='read', FILE = 'Book/LOFAR_H5files_Structure-v18.dat')
    !
    Do i_dist=1,MaxFitAntD_nr  ! set the maximal distance to be fitted to AntennaRange
      If(MaxFitAntDistcs(i_dist) .gt. AntennaRange) then
         MaxFitAntD_nr=i_dist
         MaxFitAntDistcs(i_dist)=AntennaRange
         exit
      endif
    EndDo ! i_dist=1,MaxFitAntD_nr
    Flush(unit=6)  ! standard out
    ! =====================
    If(Explore) Then
    !  need to op
      Open(unit=12,STATUS='unknown',ACTION='read',FORM ="unformatted", FILE = 'Book/RFI_Filters-v18.uft')
      Open(unit=14,STATUS='old',ACTION='read', FILE = 'Book/LOFAR_H5files_Structure-v18.dat')
      Call ExplorationRun
      stop 'Normal exploration end'
    ElseIf(Interferometry) Then
      Call InterferometerRun
      !Call EI_Run
    ElseIf(ImagingRun) Then
    ! =================================================
      Call GetNonZeroLine(lname)
      Read(lname,*) StartTime_ms, SourceGuess(:,1), StartingTime, StoppingTime  ! Start time offset = 1150[ms] for 2017 event
      Call Convert2m(SourceGuess(:,1))
      write(2,*) 'Imaging input line-1:',lname
      Write(TxtImagingPars,"(A,I3,A,F6.1,A,F4.1,A,F4.1,A,F4.0,A,F5.0,A,F5.1,A,F6.1,A,3F5.1,A,i5)") &
         'FitRange=',FitRange_Samples,'[samples],',AntennaRange,'[km], Sigma_AntT=',Sigma_AntT, &
         '[samples], SearchRangeFallOff-factor=',SearchRangeFallOff,', CCShapeCut_lim=',CCShapeCut_lim*100., &
         '%, ChiSq_lim=',ChiSq_lim,'[ns^2], EffAntNr_lim=',EffAntNr_lim*100., '%, Noise=',NoiseLevel, &
         ', Guess(N,E,h)=(', SourceGuess(:,1)/1000.,')[km], PeaksPerChunk=',PeaksPerChunk
      If(.not.FullSourceSearch) then
         Call GetNonZeroLine(lname)
         write(2,*) 'Imaging input line-2:',lname
         read(lname,*,iostat=nxx) SigmaGuess !, PeakPos(1), PeakPos(2)  ! Start time offset = 1150[ms] for 2017 event
         If(nxx.ne.0) SigmaGuess=(/ 3000.d0,3000.d0,2000.d0 /)
         Write(2,"('Preferred source location:',3F8.2,'[km]')") SourceGuess(:,1)/1000.
         Write(2,"('Initial error on location:',3F8.3,'[km]')") SigmaGuess(:)/1000.
      Endif
      Start_time(1)=(StartTime_ms/1000.)/sample  ! in sample's
      write(*,*) 'Start Source Finding'
      lname=' [s] ; ID, (North, East, vertical at core) [m] , t [s] , chi^2, sigma(N,E,v),    N_EffAnt, Max_EffAnt ;'
      lname=lname//' height, l&r widths'
      write(Txt20ImagingPars,"(A,1x,F8.2,1x,F8.2,1x, F8.2)") TRIM(FlashName), StartTime_ms, StartingTime, StoppingTime
      write(Txt20Identifier,"('Chnk#  , t[ms] , max,found, <5^2 ,',2x,A)") TRIM(release)
      Sources='Srcs'//release(2:3)
      If(Dual) then
         Sources1=TRIM(Sources)//'-dbl'//TRIM(OutFileLabel)
         Open(unit=17,STATUS='unknown',ACTION='write', FILE = TRIM(Sources1)//'.csv')
         Write(17,"('! ',A,A,A)") TRIM(TxtIdentifier),'; ',TRIM(Calibrations)
         Write(17,"('! ',A)") TRIM(TxtImagingPars)
         Write(17,"(f12.9,2A)") StartTime_ms/1000.,TRIM(lname),' even&odd antenna numbers'
         Open(unit=27,STATUS='unknown',ACTION='write', FILE = TRIM(DataFolder)//TRIM(Sources1)//'_stat.csv')
         write(27,*) TRIM(Txt20ImagingPars)
         write(27,"('! ',A)") TRIM(Txt20Identifier)
         units(2)=17
      Else
         Sources1=TRIM(Sources)//'-even'//TRIM(OutFileLabel)
         Open(unit=15,STATUS='unknown',ACTION='write', FILE = TRIM(Sources1)//'.csv')
         Write(15,"('! ',A,A,A)") TRIM(TxtIdentifier),'; ',TRIM(Calibrations)
         Write(15,"('! ',A)") TRIM(TxtImagingPars)
         Write(15,"(f12.9,2A)") StartTime_ms/1000.,TRIM(lname),' even antenna numbers'
         Open(unit=25,STATUS='unknown',ACTION='write', FILE = TRIM(DataFolder)//TRIM(Sources1)//'_stat.csv')
         write(25,*) TRIM(Txt20ImagingPars)
         write(25,"('! ',A)") TRIM(Txt20Identifier)
         units(0)=15
         !
         Sources2=TRIM(Sources)//'-odd'//TRIM(OutFileLabel)
         Open(unit=16,STATUS='unknown',ACTION='write', FILE = TRIM(Sources2)//'.csv')
         Write(16,"('! ',A,A,A)") TRIM(TxtIdentifier),'; ',TRIM(Calibrations)
         Write(16,"('! ',A)") TRIM(TxtImagingPars)
         Write(16,"(f12.9,2A)") StartTime_ms/1000.,TRIM(lname),' odd antenna numbers'
         Open(unit=26,STATUS='unknown',ACTION='write', FILE = TRIM(DataFolder)//TRIM(Sources2)//'_stat.csv')
         write(26,*) TRIM(Txt20ImagingPars)
         write(26,"('! ',A)") TRIM(Txt20Identifier)
         units(1)=16
      Endif
      !
      ChunkNr_start=StartingTime/(1000.*sample*(Time_dim-2*EdgeOffset))+1
      ChunkNr_stop=StoppingTime/(1000.*sample*(Time_dim-2*EdgeOffset))+1
      !ChunkNr_start=1501 ; ChunkNr_stop=1900 ! No Time Frame == ChunkNr
      write(*,"(A,i5,A)") achar(27)//'[45m # of data-blocks read in=',(ChunkNr_stop-ChunkNr_start),achar(27)//'[0m'
      If((ChunkNr_stop-ChunkNr_start).gt.3) Production=.true.
      Do TimeFrame=ChunkNr_start, ChunkNr_stop
         AntennaNrError=0
         write(*,"(A,i6,A)", ADVANCE='NO') achar(27)//'[31m Processing block# ',TimeFrame,achar(27)//'[0m'  ! [1000D    !  //achar(27)//'[0m.'
         i_chunk=1
         Start_Time(1)=(StartTime_ms/1000.)/sample + (TimeFrame-1)*(Time_dim-2*EdgeOffset)  ! in sample's
         !Write(2,*) 'StartTime_ms=',StartTime_ms,', TimeFrame=',TimeFrame,', Start_Time=',Start_Time(1)
         Write(2,"(A,I5,A,F8.3,A,F12.6,A,I11,A,3F10.1,A)") 'Start time for slice',TimeFrame &
            ,' (dt=',Start_time(1)*1000.d0*sample-StartTime_ms,' [ms]) is',Start_time(1)*1000.d0*sample &
            ,' [ms]=',Start_time(1),' [Samples]; SourceGuess=',SourceGuess(:,i_chunk),';  1=North, 2=East, 3=vertical(plumbline)'
         Call RFTransform_su(Time_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
           !Source_Crdnts= (/ 10000 , 16000 , 4000 /)    ! 1=North, 2=East, 3=vertical(plumbline)
         Call AntennaRead(i_chunk,SourceGuess(:,i_chunk))
         Call DAssignFFT()
         If(DataReadError.ne.0) cycle  !  set in AntennaRead
         If(AntennaNrError.gt.0) cycle  !  set in AntennaRead
         !
         Call SourceFind(TimeFrame,SourceGuess,units)
         !
         !
      EndDo !  i_chunk
      !
      If(Dual) Then
         Close(Unit=17)
         Close(Unit=27)
      Else
         Close(Unit=15)
         Close(Unit=16)
         Close(Unit=25)
         Close(Unit=26)
      Endif
      close(unit=14)
      close(unit=12)
      !
      write(*,*) NMin, NMax, Emin, EMax
      If((NMax-Nmin).gt.(Emax-Emin)) Then
         NMin=NMin -.5
         NMax=NMax +.5
         D=(EMax+Emin)/2.
         EMax=D + (NMax-Nmin)/2.
         EMin=D - (NMax-Nmin)/2.
      Else
         EMin=EMin -.5
         EMax=EMax +.5
         D=(NMax+Nmin)/2.
         NMax=D + (EMax-Emin)/2.
         NMin=D - (EMax-Emin)/2.
      Endif
      lname='AFlashImage.in'  ! temporary name holder
      Open(UNIT=10,STATUS='unknown',ACTION='WRITE',FILE=TRIM(lname))
      If(Dual) then
         write(10,"(' ',A,' 6. 2. 5.0 55  ',A,'-all ',4f5.0,' 0.0 12. ',2f7.1,' NoBox')") &
               '"'//TRIM(Sources1)//'"', TRIM(FlashName)//'d'//TRIM(OutFileLabel), &
               Emin, EMax, NMin, NMax, StartingTime, StoppingTime !  plotting command
         write(10,*) ' 0.1 0.5 0.01 " " 1.  1.   0   0.  ! no plots, MaxTrackDist[km], Wtr[0.5], TimeWin[ms^2]'
         write(10,"('======================================================')")
      Else
         write(10,"(' ',A,' 6. 2. 5.0 55  ',A,' ',4f5.0,' 0.0 12. ',2f7.1,' NoBox')") &
               '"'//TRIM(Sources1)//'"', TRIM(FlashName)//'e'//TRIM(OutFileLabel), &
               Emin, EMax, NMin, NMax, StartingTime, StoppingTime !  plotting command
         write(10,*) ' 0.1 0.5 0.01 " " 1.  1.   0   0. !     no plots, MaxTrackDist[km], Wtr[0.5], TimeWin[ms^2]'
         write(10,"(' ',A,' 6. 2. 5.0 55  ',A,' ',4f5.0,' 0.0 12. ',2f7.1,' NoBox')") &
               '"'//TRIM(Sources2)//'"', TRIM(FlashName)//'o'//TRIM(OutFileLabel), &
               Emin, EMax, NMin, NMax, StartingTime, StoppingTime !  plotting command
         write(10,*) ' 0.1 0.5 0.01 " " 1.  1.   0   0. !   no plots, MaxTrackDist[km], Wtr[0.5], TimeWin[ms^2]'
         write(10,"('======================================================')")
      EndIf
      Close(unit=10)
      !
      If(Dual) then
         Call GLEplotControl(PlotType='PeakStats', PlotName='PeakStatsD'//TRIM(OutFileLabel), &
            PlotDataFile=TRIM(DataFolder)//TRIM(Sources1) )
      Else
         Call GLEplotControl(PlotType='PeakStats', PlotName='PeakStatsE'//TRIM(OutFileLabel), &
            PlotDataFile=TRIM(DataFolder)//TRIM(Sources1) )
         Call GLEplotControl(PlotType='PeakStats', PlotName='PeakStatsO'//TRIM(OutFileLabel), &
            PlotDataFile=TRIM(DataFolder)//TRIM(Sources2) )
      EndIf
      !
      Call GLEplotControl(SpecialCmnd='cd '//TRIM(ProgramFolder))
      If(Windows) Then
         Call GLEplotControl(SpecialCmnd='gfortran -o TrackExe.exe Track-'//TRIM(Version)//'.f90 %LIBRARY%')
      Else
         Call GLEplotControl(SpecialCmnd='gfortran -o TrackExe.exe Track-'//TRIM(Version)//'.f90 ${LIBRARY}')
      EndIf
      Call GLEplotControl(SpecialCmnd='cd '//TRIM(FlashFolder))
      Call GLEplotControl(SpecialCmnd=TRIM(ProgramFolder)//'TrackExe.exe  <'//TRIM(lname))
      !
      If(Dual) then
         Call GLEplotControl(PlotType='SourcesPlot', PlotName='Map_'//TRIM(FlashName)//'d'//TRIM(OutFileLabel), &
            PlotDataFile=TRIM(DataFolder)//TRIM(FlashName)//'d'//TRIM(OutFileLabel), Submit=.true.)
      Else
         Call GLEplotControl(PlotType='SourcesPlot', PlotName='Map_'//TRIM(FlashName)//'e'//TRIM(OutFileLabel), &
            PlotDataFile=TRIM(DataFolder)//TRIM(FlashName)//'e'//TRIM(OutFileLabel))
         Call GLEplotControl(PlotType='SourcesPlot', PlotName='Map_'//TRIM(FlashName)//'o'//TRIM(OutFileLabel), &
            PlotDataFile=TRIM(DataFolder)//TRIM(FlashName)//'o'//TRIM(OutFileLabel), Submit=.true.)
      EndIf
      !
    ! ============================================================
    Else   !  Calibrate & plot
       Call RFTransform_su(Time_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
       write(*,"(A,i3,A)") achar(27)//'[45m # of data-blocks read in=',ChunkNr_dim,achar(27)//'[0m'
       Do i_chunk=1, ChunkNr_dim
           !StartTime_ms=1200 ! i_chunk*100.
           !
           Call GetNonZeroLine(lname)
           ! PeakPos is used only for the plot option
           read(lname,*) StartTime_ms, SourceGuess(:,i_chunk) !, PeakPos(1), PeakPos(2)  ! Start time offset = 1150[ms] for 2017 event
           Call Convert2m(SourceGuess(:,i_chunk))
           write(2,*) 'StartTime_ms, SourceGuess', StartTime_ms, SourceGuess(:,i_chunk)
           If(NINT(SourceGuess(2,i_chunk)).eq.1) write(*,"(A,i3,A)") achar(27)//'[45m # Check ChunkNr_dim',i_chunk,achar(27)//'[0m'
           Start_time(i_chunk)=(StartTime_ms/1000.)/sample  ! in sample's
           !Source_Crdnts= (/ 10000 , 16000 , 4000 /)    ! 1=North, 2=East, 3=vertical(plumbline)
           Call AntennaRead(i_chunk,SourceGuess(:,i_chunk))
           !
       EndDo !  i_chunk
       !
       close(unit=14)
       close(unit=12)
       call DAssignFFT()
       !Ant_ID=MaxLoc(AntMnem, mask=AntMnem.eq.DSet_Names(3))
       !
      write(*,*) 'start fitting'
      Call FindCallibr(SourceGuess) ! Find Station Callibrations
      !
      i_chunk=0 ! scratch
      If(XcorelationPlot) Call GLE_Corr()  ! opens unit 10
      If(CurtainPlot) then
         read(lname,*,iostat=nxx) StartTime_ms, SourceGuess(:,i_chunk), CurtainWidth !, PeakPos(1), PeakPos(2)  ! Start time offset = 1150[ms] for 2017 event
         Call Convert2m(SourceGuess(:,i_chunk))
         write(2,*) 'curtainplot input: ',lname
         If(nxx.ne.0 .or. CurtainWidth.le.10) CurtainWidth=300
         Call PlotAllCurtainSpectra(CurtainWidth) ! closes unit 10
      EndIf
      Call GLEplotControl( Submit=.true.)
      !
   Endif ! ImagingRun/Explore/Calibrate

    stop
end program LOFAR_Imaging
!=================================

!=================================
Subroutine AntennaSanity()
    use Chunk_AntInfo, only : BadAnt_nr, BadAnt_nr_max, BadAnt_SAI
    use Chunk_AntInfo, only : PolFlp_nr, PolFlp_nr_max, PolFlp_SAI
    use DataConstants, only : DataFolder
    Implicit none
    Integer :: StAntID
    integer :: nxx
    !
    BadAnt_nr=0
    Open(unit=9,STATUS='old',ACTION='read', FILE = 'Book/BadAntennas.dat', IOSTAT=nxx)
    If(nxx.ne.0) then
        write(2,*) '******problems with file=','Book/BadAntennas.dat'
    else
        Do
            read(9,*,IOSTAT=nxx) StAntID
            if(nxx.ne.0) then
                exit
            else
                If(BadAnt_nr.eq.BadAnt_nr_max) stop 'BadAnt_nr max exceeded'
                BadAnt_nr=BadAnt_nr+1
                BadAnt_SAI(BadAnt_nr) = StAntID  ! StAntID
            endif
        Enddo
    endif
    write(2,*) 'BadAnt_nr=',BadAnt_nr
    write(2,*) 'BadAnt_SAI',BadAnt_SAI(1:BadAnt_nr)
    close(unit=9)
    !
    PolFlp_nr=0
    Open(unit=9,STATUS='old',ACTION='read', FILE = 'Book/PolarFlipAntennas.dat', IOSTAT=nxx)
    If(nxx.ne.0) then
        write(2,*) '******problems with file=','Book/PolarFlipAntennas.dat'
    else
        Do
            read(9,*,IOSTAT=nxx) StAntID
            if(nxx.ne.0) then
                exit
            else
                If(PolFlp_nr.eq.PolFlp_nr_max) stop 'PolFlp_nr max exceeded'
                If(mod(StAntID,2) .eq. 1) stop 'PolFlp-Ant_nr should be even'
                PolFlp_nr=PolFlp_nr+1
                PolFlp_SAI(PolFlp_nr) = StAntID
            endif
        Enddo
    endif
    write(2,*) 'PolFlp_nr=',PolFlp_nr
    write(2,*) 'PolFlp_SAI',PolFlp_SAI(1:PolFlp_nr)
    close(unit=9)
    !
    Return
End Subroutine AntennaSanity
