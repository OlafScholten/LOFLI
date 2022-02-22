 !----------------
 !   Main program
 !-------------------------------------------
 !-----------------16/02/2022@14:04:57.575--------------------------
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
 !   I   M   Subroutine RFTransform_CF2RT(Cnu,RD)
 !   I   M   Subroutine RFTransform_CF2CT(Cnu,CD)
 !   I   M   Pure Subroutine R2C(R,C)
 !   I   M   Pure Subroutine  C2R(R,C)
 !   I   M   Subroutine DownSamlple(E_t,E,padding,tTrace_dim,FF_dim, filt, E_nu_dwn, inui,inum)
 !   Include 'ParamModules.f90'  ! v18d: Cnu storage changed  for InterfEngineB
 !   I   Module Chunk_AntInfo
 !   I   M   Subroutine Alloc_Chunk_AntInfo
 !   I   Module ThisSource
 !   I   M   Subroutine Alloc_ThisSource
 !   I   Module FitParams
 !   I   Module Explore_Pars
 !   I   M   Subroutine Alloc_Explore_Pars
 !   Include 'AntFunct.f90'
 !   I   Module AntFunCconst
 !   I   Subroutine Cheb_odd(Phi,Cheb)
 !   I   Subroutine Legendre_odd(x,Lgdr)
 !   I   Subroutine AntFun(thet_d,phi_d)
 !   I   Subroutine AntFun_Inv(thet_d,phi_d)
 !   I   Subroutine AntFieParGen()
 !   I   WRITE (2,*) "Using environment variable: AntennaFun; Antenna Function from files:",' "',TRIM(AntFunFile)//'LBA_Vout_*.txt"'
 !   Include 'InterferomPars.f90'
 !   I   Module Interferom_Pars
 !   I   M   Subroutine Alloc_EInterfAll_Pars
 !   I   M   Subroutine Alloc_EInterfCalib_Pars
 !   I   M   Subroutine Alloc_EInterfImag_Pars
 !   I   M   Subroutine DeAlloc_EInterfImag_Pars
 !   Include 'LOFLI_InputHandling.f90'
 !   I   Subroutine PrintIntArray(Var, Var_name, Label)
 !   I   Subroutine PrintIntVal(Var, Var_name, Label)
 !   I   Subroutine PrintRealVal(Var, Var_name, Label)
 !   I   Subroutine PrintLogiVal(Var, Var_name, Label)
 !   I   Subroutine PrintChArray(Var, Var_name, Label)
 !   I   Subroutine PrintChVal(Var, Var_name, Label)
 !   I   Subroutine PrintParIntro(Label, RunOption, OutFileLabel)
 !   I   Subroutine ReadSourceTimeLoc(StartTime_ms, CenLoc)
 !   Include 'HDF5_LOFAR_Read.f90'
 !   I   Module HDF5_LOFAR_Read
 !   I   M   Subroutine GetFileName(filenr, nxx)
 !   I   M   Subroutine ListGroups
 !   I   M   Subroutine CloseFile
 !   I   M   Subroutine ListGroupStructure(GroupName)
 !   I   M   Subroutine CloseGroup
 !   I   M   Subroutine ListDataAtt(GroupName,DSetName, prnt)
 !   I   M   Subroutine CloseDSet
 !   I   M   Subroutine GetData(Chunk, DSet_offset, DSet_dim)
 !   I   M   Subroutine AttrRead(file_id,Obj_id,Obj_name,prn)
 !   I   M   Subroutine DataRead(dset_id, Chunk, DSet_offset, DSet_dim)
 !   I   M   Subroutine GetDataChunk(GroupName,DSetName, Chunk, DSet_offset, DSet_dim, prnt, DataReadErr)
 !   I   M   Subroutine CloseDataFiles()
 !   Include 'MappingUtilities.f90'
 !   I   Module unque
 !   I   M   Subroutine unique(vec,vec_unique)
 !   I   M   Subroutine Selection_sort(a)
 !   I   M   Subroutine Double_sort(a)
 !   I   M   Pure Subroutine Double_IR_sort(N,a,R)
 !   I   M   Pure Subroutine Double_RI_sort(N,a,R)
 !   I   M   Subroutine sort(n, a)
 !   I   Module StationMnemonics
 !   I   M   Subroutine Station_ID2Mnem(STATION_ID,Station_Mnem)
 !   I   M   Character(len=5) Function Statn_ID2Mnem(STATION_ID)
 !   I   M   Subroutine Station_Mnem2ID(Station_Mnem,STATION_ID)
 !   I   M   Integer Function Statn_Mnem2ID(Station_Mnem)
 !   I   Module Calibration
 !   I   M   Subroutine ReadCalib()
 !   I   M   Subroutine Station_ID2Calib(STATION_ID,Ant_ID,StatAnt_Calib, Calibrated)
 !   I   M   Subroutine WriteCalibration ! MergeFine
 !   I   Subroutine ITRF2LOFARConstruct()
 !   I   Subroutine ITRF2LOFAR(ITRF,LOFAR)
 !   I   Subroutine RelDist(Source,LFRAnt,RDist)
 !   I   Real(kind=8) Function SubRelDist(SrcPos,i_ant,i_chunk)
 !   I   Subroutine GetNonZeroLine(lineTXT)
 !   I   Subroutine GetMarkedLine(Mark,lineTXT)
 !   I   Module ansi_colors
 !   I   M   Function color(str, code) result(out)
 !   I   Subroutine Convert2m(CenLoc)
 !   I   Pure Subroutine SetSmooth(N_Smth, Smooth)
 !   I   Function random_stdnormal() Result(x)
 !   I   end Function random_stdnormal
 !   I   Subroutine random_stdnormal3D(x)
 !   I   end Subroutine random_stdnormal3D
 !   Include 'GLEplotUtil.f90'
 !   I   Module GLEplots
 !   I   M   use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows
 !   I   M   Subroutine GLEplotControl(PlotType, PlotName, PlotDataFile, SpecialCmnd, Submit)
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
 !   I   Subroutine GetLargePeaks(i_ant, i_chunk, Peakposs)  ! Obsolete, not used anymore
 !   I   Subroutine GetStationFitOption(FP_s, FitNoSources)
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
 !   I   Subroutine Inv5StarPk(HEnvel, PeakPos, Windw, AvePos, Ampl, StDev)
 !   Include 'SourceTryal.f90'  ! v16.f90' = v17.f90'  ; v17a.f90' uses grid search
 !   Include 'ExplorationOption.f90'
 !   I   Subroutine ExplorationRun
 !   I   Subroutine GetSpectStats(i_chunk, i_expl)
 !   I   Subroutine AnalSpectStats
 !   Include 'CurtainPlotOption.f90'
 !   I   Subroutine PlotAllCurtainSpectra(CurtainWidth)
 !   I   Subroutine GLEscript_CurtainPlot(unt, file, CurtainWidth, i_Peak, UsePeakNr)
 !   I   Subroutine GLE_Corr()
 !   I   Subroutine GLEscript_Curtains(unt, file, WWidth, i_chunk, FileA, Label, dChi_ap, dChi_at, Power_p, Power_t, Chi2pDF)
 !   Include 'ImpulsImagOption.f90'
 !   I   Subroutine ImpulsImagRun
 !   I   use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, DataFolder, FlashName, Windows
 !   I   Call GLEplotControl(SpecialCmnd='cd '//TRIM(ProgramFolder))
 !   I   Call GLEplotControl(SpecialCmnd=TRIM(ProgramFolder)//'TrackExe.exe  <'//TRIM(lname))
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
 !   I   I   Subroutine EI_PolarizSlice(i_slice)
 !   I   I   Subroutine matinv3(A,B)
 !   I   Subroutine InterferometerRun
 !   I   'cp ${ProgramDir}LOFAR-Imag-v20 ./Imag-'//TRIM(OutFileLabel)//'.exe',  &  ! Just to display a more intelligent name when
 !   Include 'EIOption.f90'
 !   I   Subroutine EI_Run
 !   I   Subroutine Inverse33(A,B, eig, PolBasis)
 !   I   Subroutine EigvecS33(A,E,ev)
 !   I   Subroutine EISelectAntennas(i_chunk)
 !   I   Subroutine EIPixBoundingBox()
 !   I   Subroutine EISetupSpec(Nr_IntFer, IntfNuDim, CMCnu)
 !   I   Subroutine EIEngine(Nr_IntFer, IntfNuDim, CMCnu, CMTime_pix)
 !   I   Subroutine EIAnalyzePixelTTrace(i_N, i_E, i_h, SumWindw, IntfNuDim, CMTime_pix)
 !   Include 'ECallibrOption.f90'
 !   I   Include 'EICallRtns.f90'
 !   I   I   Subroutine EI_PrntFitPars(X)
 !   I   I   Subroutine EIPrntNewSources()
 !   I   I   Subroutine EIPrntCompactSource(i_Peak)
 !   I   I   Subroutine EIX2Source(X)
 !   I   I   Subroutine EI_PolarizPeak(i_Peak)
 !   I   I   Subroutine EISource2X_Stoch(X,Time_width, Space_Spread)
 !   I   Include 'EIFitter.f90'
 !   I   I   Subroutine EI_Fitter(X, Time_width, Space_Spread)
 !   I   I   Subroutine CompareEI( meqn, nvar, X, nf, R, uiparm, urparm, ufparm )
 !   I   I   Subroutine EI_PolGridDel(Nr_IntFer, FitDelay, i_sample, i_chunk , VoxLoc, AntPeak_OffSt, &
 !   I   I   Subroutine EI_PolSetUp(Nr_IntFer, IntfBase, i_chunk, VoxLoc, AntPeak_OffSt, Cnu0, Cnu1, W_ap, W_at)
 !   I   I   Subroutine TimeTracePlot(j_IntFer, IntfBase, i_chunk, VoxLoc,Windw, Label)
 !   I   I   Subroutine  EI_Weights(Nr_IntFer, i_sample, i_chunk, PixLoc, AntPeak_OffSt, W_ap, W_at)
 !   I   I   Pure Subroutine GetInterfFitDelay(i_chunk, FitDelay)
 !   I   I   Subroutine WriteDelChiPeak(i_chunk, DelChi,PartChiSq,PartChi2Int)
 !   I   Subroutine E_Callibr()
 !   I   Subroutine EIReadPeakInfo()
 !   Include 'System_Utilities.f90'
 !   I   Subroutine System_Initiation(Utility, release, ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows)
 !   I   Character(len=*), intent(out) :: ProgramFolder, UtilitiesFolder, FlashFolder, FlashName
 !   I   ProgramFolder='%ProgramDir%'
 !   I   ProgramFolder='${ProgramDir}'
 !   I   Subroutine System_MemUsage(valueRSS)
 !   I   end Subroutine System_MemUsage
 !   I   Subroutine CreateNewFolder(FileName)
 !   Program LOFAR_Imaging
 !   use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, DataFolder, FlashName, Windows
 !   Call System_Initiation(Utility, release, ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows)  ! set Flash name & f
 !   Subroutine AntennaSanity()
 !-----------------16/02/2022@14:04:57.575--------------------------
 !------------------------------------------
    Include 'ConstantsModules.f90'
    Include 'FFT_routines.f90'
    Include 'ParamModules.f90'  ! v18d: Cnu storage changed  for InterfEngineB
    Include 'AntFunct.f90'
    Include 'InterferomPars.f90'
    Include 'LOFLI_InputHandling.f90'
    Include 'HDF5_LOFAR_Read.f90'
    Include 'MappingUtilities.f90'
    Include 'GLEplotUtil.f90'
    Include 'Ant-Read.f90'
    Include 'CrossCorr.f90'
    Include 'FindCallibr.f90'   ! Station Callibration
    Include 'FitParams.f90'
    Include 'Fitter_CorrMax.f90'     ! uses chi-square mini
    Include 'FindSources.f90'
    Include 'SourceTryal.f90'  ! v16.f90' = v17.f90'  ; v17a.f90' uses grid search
    Include 'ExplorationOption.f90'
    Include 'CurtainPlotOption.f90'
    Include 'ImpulsImagOption.f90'
    Include 'InterferometryOption.f90'  ! d: Cnu storage changed for InterfEngineB
    Include 'EIOption.f90'
    !Include 'EIStokes-testW.f90'
    !Include 'SIOption.f90'
    Include 'ECallibrOption.f90'
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
! v22.01: chi^2 fitting of antenna delays in TRI-D mode for calibration
    use constants, only : dp,pi,ci,sample,Refrac
    use LOFLI_Input
    use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, DataFolder, FlashName, Windows
    use DataConstants, only : Time_dim, Cnu_dim, Production, RunMode, Utility, release
    use DataConstants, only : PeakNr_dim, ChunkNr_dim, Diagnostics, EdgeOffset, Calibrations, OutFileLabel
    use DataConstants, only : Polariz
    use ThisSource, only : Alloc_ThisSource, Dual, RealCorrelation, t_ccorr, Safety, CurtainHalfWidth
    use ThisSource, only : CCShapeCut_lim, ChiSq_lim, EffAntNr_lim, NrP, PeakNrTotal ! , PeakPos
    use FitParams, only : FitIncremental, WriteCalib, FullAntFitPrn,  AntennaRange
    use FitParams, only : MaxFitAntDistcs, MaxFitAntD_nr, Sigma_AntT, SearchRangeFallOff  ! ,PeakS_dim, Explore
    use FitParams, only : FullSourceSearch ! , SigmaGuess
    use Chunk_AntInfo, only : ExcludedStatID, Start_time, BadAnt_nr, BadAnt_SAI, DataReadError, TimeFrame
    use Chunk_AntInfo, only : NoiseLevel, PeaksPerChunk, TimeBase, Simulation, WriteSimulation, CalibratedOnly
    use Chunk_AntInfo, only : ExcludedStat_max, SgnFlp_nr, PolFlp_nr, SignFlp_SAI, PolFlp_SAI, BadAnt_SAI, AntennaNrError
    use Chunk_AntInfo, only : Alloc_Chunk_AntInfo, ExcludedStat, SaturatedSamplesMax, N_Chunk_max
    use FFT, only : RFTransform_su,DAssignFFT
    use StationMnemonics, only : Statn_ID2Mnem, Statn_Mnem2ID
    Use Interferom_Pars, only : IntfPhaseCheck, N_smth, ChainRun, PixPowOpt, N_fit
    use GLEplots, only : GLEplotControl
    Implicit none
    !
    Integer :: i,j,i_chunk, i_Peak, units(0:2), FitRange_Samples, IntfSmoothWin !, CurtainHalfWidth
    Real*8 :: StartTime_ms, StartingTime, StoppingTime, D
    Real*8 :: SourceGuess(3,N_Chunk_max) ! = (/ 8280.01,  -15120.48,    2618.37 /)     ! 1=North, 2=East, 3=vertical(plumbline)
    Integer, parameter :: lnameLen=180
    CHARACTER(LEN=1) :: Mark
    CHARACTER(LEN=6) :: txt,Version
    CHARACTER(LEN=10) :: Sources
    Character(LEN=30) :: RunOption, FMTSrces
    Character(LEN=180) :: Sources1, Sources2
    Character(LEN=lnameLen) :: lname
    Character(LEN=250) :: TxtIdentifier, TxtImagingPars
    Character(LEN=250) :: Txt20Identifier, Txt20ImagingPars
    Logical :: XcorelationPlot, prnt=.false. ! CrossCorrelate=.false. ,
    Logical :: Interferometry=.false.
    Logical :: E_FieldsCalc=.true.
    Integer :: i_dist, i_guess, nxx, valueRSS
    NAMELIST /Parameters/ RunOption &
         !,  Explore, ImagingRun, Interferometry !, RealCorrelation, E_FieldsCalc &
         , FullSourceSearch, CurtainHalfWidth, XcorelationPlot &
         , IntfPhaseCheck, IntfSmoothWin, TimeBase  &
         , Diagnostics, Dual, FitIncremental &
         , Simulation, SignFlp_SAI, PolFlp_SAI, BadAnt_SAI, SaturatedSamplesMax, Calibrations, WriteCalib, CalibratedOnly &
         , ExcludedStat, FitRange_Samples, FullAntFitPrn, AntennaRange, PixPowOpt, OutFileLabel, ChainRun &
         , CCShapeCut_lim, ChiSq_lim, EffAntNr_lim, Sigma_AntT, SearchRangeFallOff, NoiseLevel, PeaksPerChunk    !  ChunkNr_dim,
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Version='v22.01'
   release='22.01 (Jan, 2022)'
   Utility='LOFLi-Imaging'
   CALL get_environment_variable("LIBRARY", lname)
   !WRITE (*,*) "LIBRARY=",TRIM(lname)
   CALL get_environment_variable("HOME", lname)
   WRITE (*,*) "HOME=",TRIM(lname),', ',Trim(Utility), ', release=', TRIM(release)
   !
   ExcludedStat='  '
   ExcludedStatID=0
   SignFlp_SAI=-1
   BadAnt_SAI=-1
   PolFlp_SAI=-2
   XcorelationPlot=.false.
   FitRange_Samples=Safety
   RunOption="None"
   CCShapeCut_lim=0.6 ; ChiSq_lim=80. ; EffAntNr_lim=0.8
   NoiseLevel=70.
   IntfSmoothWin=20
   CurtainHalfWidth=-1
   SaturatedSamplesMax=5
   Simulation=""
   CalibratedOnly=.true.
   !reshape((/ -1.D-10,0.d0,0.d0,0.d0,-1.D-10,0.d0,0.d0,0.d0,-1.D-10 /), shape(array))
   !AntennaRange_km=AntennaRange
   read(*,NML = Parameters,iostat=nxx)
   j = iachar(RunOption(1:1))  ! convert to upper case if not already
   if (j>= iachar("a") .and. j<=iachar("z") ) then
      RunOption(1:1) = achar(j-32)
   end if
      SELECT CASE (RunOption(1:1))
   ! Explore
   ! Calibrate
   ! SelectData
   ! XY-Interferometry
   ! PolInterferometry==TRI-D imager
   ! ImpulsiveImager
      CASE("S")  ! SelectData
         Write(*,*) 'SelectData run, OutFileLabel=',TRIM(OutFileLabel)
         OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='SelectData'//TRIM(OutFileLabel)//'.out')
         RunMode=0 ! GLE utilities are not called.
         Sources1='SelectData'
      CASE("E")  ! Explore
         Write(*,*) 'Explore run, OutFileLabel=',TRIM(OutFileLabel)
         OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='Explore'//TRIM(OutFileLabel)//'.out')
         RunMode=1
         Sources1='Explore'
      CASE("C")  ! Calibrate
         Write(*,*) 'Calibration run, OutFileLabel=',TRIM(OutFileLabel)
         OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='Calibrate'//TRIM(OutFileLabel)//'.out')
         RunMode=2
         Sources1='Calibration'
      CASE("I")  ! ImpulsiveImager
         Write(*,*) 'Impulsive Imager run, OutFileLabel=',TRIM(OutFileLabel)
         OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='Imaging'//TRIM(OutFileLabel)//'.out')
         RunMode=3
         Sources1='Impulsive Imager'
      CASE("X")  ! XY-Interferometry
         Write(*,*) 'Interferometry, OutFileLabel=',TRIM(OutFileLabel)
         OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='Interferometer'//TRIM(OutFileLabel)//'.out')
         RunMode=14
         Write(2,*) 'This option has been discontinued, use "T" instead.'
         Stop 'Discontinued option'
         Sources1='XY- TRI-D Imager (Interferometry)'
      CASE("T")  ! TRI-D imager==PolInterferometry
         Write(*,*) 'Interferometry, OutFileLabel=',TRIM(OutFileLabel)
         OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='Interferometer'//TRIM(OutFileLabel)//'.out')
         RunMode=24
         Sources1='Polarimetric TRI-D Imager (Interferometry)'
      CASE("F")  ! FieldCalibrate for TRI-D imager
         Write(*,*) 'Interferometric Calibration for TRI-D mode, OutFileLabel=',TRIM(OutFileLabel)
         OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='FC'//TRIM(OutFileLabel)//'.out')
         RunMode=7
         Sources1='Field Calibration for TRI-D Imager (Interferometry)'
      CASE DEFAULT  ! Help
         Write(*,*) 'specified RunOption: "',Trim(RunOption),'", however the possibilities are:'
         Write(*,*) '- "Explore" for first exploration of this flash to get some idea of the layout and timing'
         Write(*,*) '- "Calibrate" for performing time calibration using the Hilbert envelopes of the cross correlations'
         Write(*,*) '- "ImpulsiveImager" for the impulsive Imager'
         Write(*,*) '- "FieldCalibrate" Field Calibration for the TRI-D interferometric imager'
         Write(*,*) '- "TRI-D" for the TRI-D imager with polarization observables, accounting for antenna function'
         Write(*,*) '- "SelectData" to select real data, possibly for setup of simulation runs using program "SimulateData"'
         Stop 'No valid RunOption specified in namelist input'
   End SELECT
   !
   Call System_Initiation(Utility, release, ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows)  ! set Flash name & folder names
   Call PrintParIntro(TRIM(Sources1), RunOption,OutFileLabel)
   !
   If(nxx.ne.0) Then
      write(2,*) 'Solely use parameters in the namelist from the following list'
      write(2,NML = Parameters)
      Stop 'Namelist reading problem'
   EndIf
   !
   N_fit=0
   ChunkNr_dim=1
   PeakNr_dim=1
   FMTSrces="(1x,3i2,I8,3(F10.2,1x))"
   SELECT CASE (RunOption(1:1))
   ! option          runmode
   ! Explore            1
   ! Calibrate          2
   ! SelectData         0, same as for RFI-Mitigate
   ! FieldCalibrate     7
   ! TRI-D imager       4
   ! ImpulsiveImager    3
   ! Other, TrackScatt  5
   ! Other, InterfSelct 6
   ! Other, RFI-Mitigat 0
      CASE("S")  ! SelectData                          RunMode=0
         If(Simulation.eq."") Simulation="Sims/Sims"
         If(WriteSimulation(1).le.0) WriteSimulation(1)=200  ! a reasonable save buffer
         If(WriteSimulation(2).le.0) WriteSimulation(2)=3
         Call PrintValues(Simulation, 'Simulation', 'Naming for the "Simulation" files with the antenna info and data.' ) !
      CASE("E")  ! Explore                          RunMode=1
         !ImagingRun=.false.
         Production=.true.
         Interferometry=.false.
         !Production=.false.
         FitRange_Samples=70
         AntennaRange=2.49999999
         Dual=.false.
         PeakNr_dim=10
         NoiseLevel=20.
         WriteSimulation(2)=-1
         Simulation=""
      CASE("C")  ! Calibrate                          RunMode=2
         ChunkNr_dim=0
         Do i=1,N_Chunk_max ! perform some pre-scanning of the input to now the number of chunks that will be used
            Call GetNonZeroLine(lname)
            Read(lname(2:lnameLen),*,iostat=nxx) StartTime_ms, SourceGuess(:,1)  ! just dummy arguments
            Call Convert2m(SourceGuess(:,1))
            If(nxx.ne.0) exit
            Read(lname,FMTSrces,iostat=nxx) &
                i_dist, i_guess,j ,i_chunk,  SourceGuess(:,1) ! just dummy arguments
            !write(2,*) i, nxx, lname
            If(nxx.eq.0) exit
            ChunkNr_dim=i   ! this was a genuine chunk card
         EndDo
         If(ChunkNr_dim.eq.0) Then
            Write(2,*) lname
            write(2,*) 'ChunkNr_dim:',ChunkNr_dim, 'last line:', StartTime_ms,i_dist, i_guess,j,i_chunk
            stop 'Chunk number too small'
         EndIf
         If(ChunkNr_dim.gt.N_Chunk_max) Then
            Write(2,*) lname
            write(2,*) 'ChunkNr_dim:',ChunkNr_dim, 'last line:', StartTime_ms,i_dist, i_guess,j,i_chunk
            stop 'Chunk number too large'
         EndIf
         !  Determine number of peaks/sorces that are included in the calibration search
         write(2,*) 'number of calibration chunks:',ChunkNr_dim
         i_peak=0
         Do
            Read(lname,FMTSrces,iostat=nxx) &
               i_dist, i_guess,j ,i_chunk, SourceGuess(:,1) ! just dummy arguments
            If(nxx.ne.0) exit
            i_peak=i_peak+1
            Call GetNonZeroLine(lname)
            read(lname,*)  txt  ! check for possible 'exclude' line following this
            If(trim(txt).eq.'exclud') Then
               Call GetNonZeroLine(lname)
            EndIf
         Enddo
         If(i_peak.eq.0) Then
            PeakNr_dim=2*NrP*ChunkNr_dim
         Else
            PeakNr_dim=i_peak
         EndIf
         !
         write(2,*) 'number of calibration sources:',PeakNr_dim
         Rewind(unit=5) ! Standard input
         read(*,NML = Parameters) ! to reposition correctly
         !
         WriteSimulation(2)=-1
         CalibratedOnly=.False.
         E_FieldsCalc=.false.
         RealCorrelation=.False.
         Polariz= E_FieldsCalc
         Call PrintValues(CurtainHalfWidth,'CurtainHalfWidth', 'Produce a "Curtain" plot when positive.')  ! width of plot?
         Call PrintValues(XcorelationPlot,'XcorelationPlot', &
         'Produce a plot of the cross-correlation functions (real or absolute).')
         Call PrintValues(Diagnostics,'Diagnostics', 'Print diagnostics information, creates much output.')
         Call PrintValues(FullAntFitPrn,'FullAntFitPrn', &
         'Print detailed information for the time deviations per antenna per pulse.')
         Call PrintValues(Simulation,'Simulation', 'Run on simulated data from such files.' ) !
         Call PrintValues(Dual,'Dual', 'Make a combined analysis of the even (Y-) and odd (X-) numbered dipoles.')
         Call PrintValues(FullSourceSearch,'FullSourceSearch', &
         'Perform without any preferred direction, otherwise take the sourceguess as a preference.')
         Call PrintValues(FitIncremental,'FitIncremental', &
         'Slowly increasing the range of the antennas to find the source position.')
         Call PrintValues(AntennaRange,'AntennaRange', &
            'Maximum distance (from the core, in [km]) for  antennas to be included.')
         Call PrintValues(FitRange_Samples,'FitRange_Samples', &
         'Maximum deviation (in samples) for a pulse in an antenna to be included in source position fitting.')
         !Call PrintValues(Fit_AntOffset,'Fit_AntOffset', 'Off-sets per antenna are searched and fitted.')
         Call PrintValues(WriteCalib,'WriteCalib', 'Write out an updated calibration-data file.')
         !Call PrintValues(PeakNr_dim,'PeakNr_dim', &
         !'Maximum number of sources (counting even and odd dipoles separately) included in the input.'//&
         !'Note: this is a bit of an annoying variable meant to allow to fit more offsets simultaneously.')
         !
         ! Pre-process inputdata for sources to be used in calibration
         !
      CASE("I")  ! ImpulsiveImager                          RunMode=3
         !ImagingRun=.true.
         ChunkNr_dim=1
         PeakNr_dim=2
         Production=.true.
         CurtainHalfWidth=-1
         Diagnostics=.false.
         FullAntFitPrn=.false.
         XcorelationPlot=.false.
         FitIncremental=.true.
         Interferometry=.false.
         Polariz=.false.  ! ????
         RealCorrelation=.false.
         !FitRange_Samples=170
         FitRange_Samples=86 ! 90  !  some cases 70 is better, sometimes 90, 80 seems to be worse i.e. highly non-linear!
         Call PrintValues(Dual,'Dual', 'Fix pulses in the even (Y-) and odd (X-) numbered dipoles at same source position.')
         Call PrintValues(PeaksPerChunk,'PeaksPerChunk', &
         'Maximum number of sources searched for per chunk (of 0.3 ms).' ) !
         Call PrintValues(NoiseLevel,'NoiseLevel', 'Any weaker sources will not be imaged.' ) !
         Call PrintValues(CalibratedOnly,'CalibratedOnly', &
            'Use only antennas that have been calibrated.')
         Call PrintValues(AntennaRange,'AntennaRange', &
            'Maximum distance (from the core, in [km]) for  antennas to be included.')
         Call PrintValues(FullSourceSearch,'FullSourceSearch', &
         'Perform without any preferred direction, otherwise take the sourceguess as a preference.')
         Call PrintValues(Simulation,'Simulation', 'Run on simulated data from such files.')  !
         Call PrintValues(EffAntNr_lim,'EffAntNr_lim', &
            'minimal fraction of the total number of antennas for which the pulse is located.' ) !
         Call PrintValues(CCShapeCut_lim,'CCShapeCut_lim', &
            'Maximum ratio of the width of the cross correlation function by that of the self correlation.')
         Call PrintValues(ChiSq_lim,'ChiSq_lim', 'Do not save any sources with a worse chi-square.' ) !
         Call PrintValues(SearchRangeFallOff,'SearchRangeFallOff', &
            'Multiplier for the (parabolic) width of the pulse-search window.')
         Call PrintValues(Sigma_AntT,'Sigma_AntT', 'Constant added to width of search window for pulses.' ) !
      CASE("F")  ! Field Calibration; Interferometric               RunMode=7
         ! Pre-process inputdata for sources to be used in calibration
         ChunkNr_dim=0
         Do i=1,N_Chunk_max ! perform some pre-scanning of the input to now the number of chunks that will be used
            Call GetNonZeroLine(lname)
            Read(lname(2:lnameLen),*,iostat=nxx) StartTime_ms, SourceGuess(:,i)  ! just dummy arguments
            Call Convert2m(SourceGuess(:,i))
            !Write(2,*) nxx, 'A:', lname
            If(nxx.ne.0) exit
            Read(lname,FMTSrces,iostat=nxx) &
               i_dist, i_guess, i_chunk,j, SourceGuess(:,i) ! just dummy arguments
            !Write(2,*) nxx, 'B:', lname
            If(nxx.eq.0) exit
            ChunkNr_dim=i   ! this was a genuine chunk card
         EndDo
         If(ChunkNr_dim.eq.0) Then
            write(2,*) 'ChunkNr_dim:',ChunkNr_dim, 'last line:', StartTime_ms,i_dist, i_guess,j,i_chunk
            stop 'Chunk number too small for FC'
         EndIf
         If(ChunkNr_dim.gt.N_Chunk_max) Then
            write(2,*) 'ChunkNr_dim:',ChunkNr_dim, 'last line:', StartTime_ms,i_dist, i_guess,j,i_chunk
            stop 'Chunk number too large for FC'
         EndIf
         !  Determine number of peaks/sorces that are included in the calibration search
         write(2,*) 'number of calibration chunks:',ChunkNr_dim
         i_peak=0
         Do
            Read(lname,FMTSrces,iostat=nxx) &
               i_dist, i_guess, i_chunk,j,  SourceGuess(:,1) ! just dummy arguments
            !write(2,*) 'nxx:',nxx,trim(lname)
            If(nxx.ne.0) exit
            If(I_guess.eq.0) i_peak=i_peak+1
            Call GetNonZeroLine(lname)
            read(lname,*)  txt  ! check for possible 'exclude' line following this
            If(trim(txt).eq.'exclud') Then
               Call GetNonZeroLine(lname)
            EndIf
         Enddo
         PeakNr_dim=i_peak
         PeakNrTotal=i_peak  ! one of these is obsolete now
         !
         write(2,*) 'number of calibration sources:',PeakNr_dim
         Rewind(unit=5) ! Standard input
         read(*,NML = Parameters) ! to reposition correctly
         !
         Interferometry=.true.
         FitRange_Samples=7
         Dual=.false.
         CurtainHalfWidth=-1
         Polariz=Dual  !
         CalibratedOnly=.False.
         If(IntfSmoothWin.lt.3) IntfSmoothWin=3
         N_smth=IntfSmoothWin
         !Call PrintValues(CurtainHalfWidth,'CurtainHalfWidth', 'Produce a "Curtain" plot when positive.')  ! width of plot?
         Call PrintValues(XcorelationPlot,'XcorelationPlot', &
            'Produce a plot of the cross-correlation functions (real or absolute).')
         Call PrintValues(Diagnostics,'Diagnostics', 'Print diagnostics information, creates much output.')
         Call PrintValues(AntennaRange,'AntennaRange', &
            'Maximum distance (from the core, in [km]) for  antennas to be included.')
         Call PrintValues(WriteCalib,'WriteCalib', 'Write out an updated calibration-data file.')
         Call PrintValues(IntfSmoothWin,'IntfSmoothWin', 'Width (in samples) of the slices for TRI-D imaging.')
         !
      CASE("T")  ! PolInterferometry==TRI-D imager         RunMode=24
         Interferometry=.true.
         FitRange_Samples=7
         CurtainHalfWidth=-1
         Dual=.true.
         Polariz=Dual  ! Affects reading in LOFAR Data; even/odd pairs only & equal delay for the pair
         PeakNr_dim=2
         If(IntfSmoothWin.lt.3) IntfSmoothWin=3
         N_smth=IntfSmoothWin
         Call PrintValues(AntennaRange,'AntennaRange', &
            'Maximum distance (from the core, in [km]) for  antennas to be included.')
         Call PrintValues(TimeBase,'TimeBase', &
            'Time-offset from the start of the data, best if kept the same for all analyses for this flash')
         !Call PrintValues(CurtainHalfWidth,'CurtainHalfWidth', 'Produce a "Curtain" plot when positive.' ) !
         Call PrintValues(Simulation,'Simulation', 'Run on simulated data from such files.' )
         Call PrintValues(ChainRun,'ChainRun', &
            'Automatically start jobs (for - previous or + following timeslots, made to follow a negative leader.' )
         Call PrintValues(IntfSmoothWin,'IntfSmoothWin', 'Width (in samples) of the slices for TRI-D imaging.' )
         Call PrintValues(PixPowOpt,'PixPowOpt', &
            '=0=default: Intensity=sum two transverse polarizations only; '//&
            '=1: Intensity=sum all polarizations weighted with alpha == intensity of F vector; '//&
            '=2: Intensity=sum all three polarizations, including longitudinal with the full weight')
         Call PrintValues(CalibratedOnly,'CalibratedOnly', &
            'Use only antennas that have been calibrated.')
         Call PrintValues(NoiseLevel,'NoiseLevel', 'Any weaker sources will not be imaged.' ) !
      CASE DEFAULT  ! Help
         Write(2,*) 'Should never reach here!'
   End SELECT
   !
   Call PrintValues(Calibrations,'Calibrations', &
   'The antenna time calibration file. Not used when running on simulated data!' )
   Call PrintValues(SignFlp_SAI,'SignFlp_SAI', &
   'Station-Antenna Identifiers for those where the sign of the signal should be reversed.')
   Call PrintValues(PolFlp_SAI,'PolFlp_SAI', &
   'Station-Antenna Identifiers for those where the even-odd signals should be interchanged.')
   Call PrintValues(BadAnt_SAI,'BadAnt_SAI', &
   'Station-Antenna Identifiers for those that are malfunctioning.' )
   Call PrintValues(ExcludedStat,'ExcludedStat', &
   'Mnemonics of the stations that should be excluded.' )
   Call PrintValues(SaturatedSamplesMax,'SaturatedSamplesMax', &
      'Maximum number of saturates time-samples per chunk of data')
   write(2,*) '&end'
   write(2,"(20(1x,'='))")
   !     Some remaining setup
   If(FitRange_Samples.lt.5) FitRange_Samples=5
   Safety=FitRange_Samples
   Call Alloc_Chunk_AntInfo !  uses ChunkNr_dim
   Call Alloc_ThisSource    !  uses ChunkNr_dim
   Do i=-Safety,Safety
     t_ccorr(i)=i ! time axis needed for spline interpolation of CrossCorrelation
   Enddo
   !write(2,*)'refractivity=',Refrac-1.
   SgnFlp_nr=count(SignFlp_SAI.gt.0)
   PolFlp_nr=count(PolFlp_SAI.gt.0)  !first pole-flip, then bad antenna
   BadAnt_nr=count(BadAnt_SAI.gt.0)
   If(any(mod(PolFlp_SAI,2) .eq. 1,1) ) Then
      write(2,*) 'PolFlp-Ant_nr should be even'
      write(*,*) 'PolFlp-Ant_nr should be even'
   EndIf
   write(2,*) '# PolarizationFlips=',PolFlp_nr,'; # BadAntennas=',BadAnt_nr,CHAR(13)//CHAR(10) ! , &
    !'Color test'//achar(27)//'[95m pink '//achar(27)//'[0m.'
!
   Do j=1,ExcludedStat_max
      If(ExcludedStat(j).eq.'     ') exit
      ExcludedStatID(j)=Statn_Mnem2ID(ExcludedStat(j))
      write(2,*) j,ExcludedStatID(j),ExcludedStat(j)
   enddo
   !
   Do i_dist=1,MaxFitAntD_nr  ! set the maximal distance to be fitted to AntennaRange
      If(MaxFitAntDistcs(i_dist) .gt. AntennaRange) then
         MaxFitAntD_nr=i_dist
         MaxFitAntDistcs(i_dist)=AntennaRange
         exit
      endif
   EndDo ! i_dist=1,MaxFitAntD_nr
   flush(unit=2)
   !
   ! Start the more serious stuff
   SELECT CASE (MOD(RunMode,10))
      ! ========0 0 0 0 0 0 0 0 0 0 0 0 0
      CASE(0)  ! SelectData                          RunMode=0
         i_chunk=1
         Call ReadSourceTimeLoc(StartTime_ms, SourceGuess(:,i_chunk))
         !Read(lname,*) StartTime_ms, SourceGuess(:,1), StartingTime, StoppingTime  ! Start time offset = 1150[ms] for 2017 event
         !Call Convert2m(SourceGuess(:,1))
         !write(2,*) 'Simulation setup input line-1:',lname
         !Call GetNonZeroLine(lname)
         Call GetMarkedLine(Mark,lname)
         !Read(lname(2:lnameLen),*,iostat=nxx) WriteSimulation(1:2)
         !write(2,"(A)") 'Input line-2: "'//lname(1:1)//'|'//TRIM(lname(2:lnameLen))// &
         !   '" !  First/Median| sample number & window length of data written to "Simulation" files'
         write(2,"(A)") 'Input line-2: "'//Mark//'|'//TRIM(lname)// &
            '" !  First/Median| sample number & window length of data written to "Simulation" files'
         Read(lname,*,iostat=nxx) WriteSimulation(1:2)
         If(  WriteSimulation(1) .le. 0. .or.  WriteSimulation(2).le. 0) Then
            write(2,*) 'Both window definition numbers should be positive'
            stop 'Window definition error'
         EndIf
         If(nxx.ne.0) then
            write(2,*) 'Parsing error in 2nd line:'
            stop 'Interferometry; reading error'
         Endif
         If(Mark.eq.'M') Then
            WriteSimulation(1)=WriteSimulation(1)-WriteSimulation(2)/2
         EndIf
         Write(2,*) 'Window start at sample#=',WriteSimulation(1), ', median at ',WriteSimulation(1)+WriteSimulation(2)/2 &
            ,', window=',WriteSimulation(2),' [samples]'
         Start_Time(1)=(StartTime_ms/1000.)/sample
         write(2,"(20(1x,'='),/)")
         Call RFTransform_su(Time_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         Call AntennaRead(i_chunk,SourceGuess(:,i_chunk))
         Call DAssignFFT()
         Stop 'SelectData'
      ! ========1 1 1 1 1 1 1 1 1 1 1
      CASE(1)  ! Explore                          RunMode=1
         Open(unit=12,STATUS='unknown',ACTION='read',FORM ="unformatted", FILE = 'Book/RFI_Filters-v18.uft')
         Open(unit=14,STATUS='old',ACTION='read', FILE = 'Book/LOFAR_H5files_Structure-v18.dat')
         Call ExplorationRun
         stop 'Normal exploration end'
    ! ========4 4 4 4 4 4 4 4 4 4 4
      CASE(4)  ! PTRI-D imager         RunMode=14 & 24
         Call InterferometerRun
      ! =============3 3 3 3 3 3 3 3 3
      CASE(3)  ! ImpulsiveImager                          RunMode=3
         Call ImpulsImagRun
       ! =====2 2 2 2 2 2 2 2
      CASE(2)  ! Calibrate                          RunMode=2
          Call RFTransform_su(Time_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          write(*,"(A,i3,A)") achar(27)//'[45m # of data-blocks read in=',ChunkNr_dim,achar(27)//'[0m'
          Do i_chunk=1, ChunkNr_dim
              !StartTime_ms=1200 ! i_chunk*100.
              !
              !Call GetNonZeroLine(lname)
              !read(lname,*) StartTime_ms, SourceGuess(:,i_chunk) !
              !Call Convert2m(SourceGuess(:,i_chunk))
              Call ReadSourceTimeLoc(StartTime_ms, SourceGuess(:,i_chunk))
              !write(2,*) 'StartTime_ms, SourceGuess', StartTime_ms, SourceGuess(:,i_chunk)
              If(NINT(SourceGuess(2,i_chunk)).eq.1) &
                     write(*,"(A,i3,A)") achar(27)//'[45m # Check ChunkNr_dim',i_chunk,achar(27)//'[0m'
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
         If(CurtainHalfWidth.gt.0) then
            If(CurtainHalfWidth.le.10) CurtainHalfWidth=300
            write(2,*) 'curtainplot half width: ',CurtainHalfWidth
            Call PlotAllCurtainSpectra(CurtainHalfWidth) ! closes unit 10
         EndIf
         Call GLEplotControl( Submit=.true.)
         !
      ! =====7 7 7 7 7 7 7 7 7
      CASE(7)  ! Field Calibrate                          RunMode=7
         Call E_Callibr()
         !
         !
         i_chunk=0 ! scratch
         If(XcorelationPlot) Call GLE_Corr()  ! opens unit 10
         Call GLEplotControl( Submit=.true.)
         !
  End SELECT

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
!
