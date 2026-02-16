!======================
Module Chunk_AntInfo
    use constants, only : dp
    use DataConstants, only : Time_dim !, Station_nrMax, HStation_nrMax  ! Parameters, fixed, Station_nrMax applies to LBA only
    use DataConstants, only : ChunkNr_dim, Ant_nrMax, Used_StationNr  ! applies to all
    use DataConstants, only : LAnt_nrMax, Used_LStationNr, HAnt_nrMax, Used_HStationNr  ! applies to LBA or HBA only
    Integer, parameter :: LOFARNrMax=39  ! smaller or equal to Station_nrMax
    Character(len=5), save :: LOFAR_name(1:LOFARNrMax)=(/&
        "CS001","CS002","CS003","CS004","CS005","CS006","CS007","CS011","CS013","CS017",&
        "CS021","CS024","CS026","CS030","CS031","CS032","CS101","CS103","CS201","CS301",&
        "CS401","CS501","RS106","RS205","RS208","RS210","RS305","RS306","RS307","RS310",&
        "RS406","RS407","RS409","RS503","RS508","CS302","CS028","RS509","XS   "  /)
    Integer, save :: LOFAR_number(1:LOFARNrMax)=&
        (/   1 ,     2 ,     3 ,     4 ,     5 ,     6 ,     7 ,    11 ,    13 ,    17 ,&
            21 ,    24 ,    26 ,    30 ,    31 ,    32 ,   101 ,   103 ,   121 ,   141 ,&
           161 ,   181 ,   106 ,   125 ,   128 ,   130 ,   145 ,   146 ,   147 ,   150 ,&
           166 ,   167 ,   169 ,   183 ,   188 ,   142 ,    28 ,   189 ,     0    /)
    Character(len=80), save :: Simulation
    Integer, save :: WriteSimulation(2)=(/ -1,-1 /)  ! start sample, window size
    ! Convention for i_ant
    ! i_ant runs from 1 till  Ant_nr(i_chunk), a sequential number for all antennas for which there are data in this chunk
    !    where i_ant=[1, LBA_nr(i_chunk)] are the LBA antennas
    !    and i_ant=[LBA_nr(i_chunk)+1, Ant_nr(i_chunk)]  are the HBA antennas
    !
    ! Ant_IDs(i_ant, i_chunk) gives the antenna ID [1,94]
    ! Ant_Stations(i_ant, i_chunk)=Station_ID  is ID number for the stations =10*LOFAR_ID + (0 for LBA .or. 1 for HBA)
    !    This will distinguish HBA and LBA
    !    MODULO(Unique_StatID,10)  Gives the LOFAR-ID of the combined LBA-HBA station
    !    i_type=MOD(Unique_StatID,10)  Gives the antenna type; 0 for LBA & 1 for HBA
    ! AntType(i_ant, i_chunk) gives 0 for LBA and 1 for HBA; may not be a necessary parameter
    !
    ! Be careful with i_stat since
    ! i_stCh is the sequential station number for the stations in a particular chunk
    ! i_stUq corresponds to the sequential station number in Unique_StatID
    !
    ! ======= used in Calibration fitting =====================================
    ! Tot_UniqueAnt, Unique_SAI are updated in "Find_unique_StatAnt()" in FitParams.f90  --or-- in "EISelectAntennas" in EIOption.f90, nowhere else seemingly
    ! Unique_StatID(i_stat)
    !       Function Statn_ID2Mnem(Unique_StatID(i_stat)) returns something like "CS002L"
    !
    ! Unique_SAI(i_ant)=100*Unique_StatID(i_stat)+Ant_IDs(i_ant, i_chunk)  ! constructed in "Find_unique_StatAnt" in FitParams.f90
    ! Tot_UniqueAnt(i_stat)    ! highest ranknr in Unique_SAI of unique_antenna for station Unique_StatID(i_stat)
    !
    ! Ant_nr(i_chunk)   ! Total number of LBA + HBA antennas for this chunk
    ! LBA_nr(i_chunk)   ! Total number of LBA antennas
    ! HBA_nr(i_chunk)   ! Total number of HBA antennas
    ! RefAnt(i_chunk, i_eo, i_type)      ! RefAnt(1:ChunkNr_dim,0:1, 0:1), constructed in Subroutine GetRefAnt(ChunkNr) in CrossCorr.f90
    ! HRefAnt(i_chunk, i_eo)
    ! ChunkStTime_ms(i_chunk) ! Chunk Start Time
    ! CTime_spectr(1:Time_dim, i_ant, i_chunk) ! stores RFI-cleaned time-traces for LBA
    !       where i_ant=[1, LBA_nr(i_chunk)] are the LBA antennas
    ! CTime_Hspectr(1:2*Time_dim, i_HAnt, i_chunk) ! stores RFI-cleaned time-traces for HBA
    !       where i_HAnt=i_ant-LBA_nr(i_chunk)=[1, HBA_nr(i_chunk)]  are the HBA antennas
!
    Integer, save, allocatable :: Ant_Stations(:,:)
    Integer, save, allocatable :: Ant_IDs(:,:)
    Integer, save, allocatable :: Ant_nr(:)
    Integer, save, allocatable :: LBA_nr(:)
    Integer, save, allocatable :: HBA_nr(:)
    Integer, save, allocatable :: RefAnt(:,:,:)
    Integer, save, allocatable :: StartT_sam(:)  ! ChunkNr(:)  !
    Real(dp), save, allocatable  :: ChunkStTime_ms(:), ChunkFocus(:,:)
    Real(dp), save, allocatable :: Ant_RawSourceDist(:,:) ! distances that have been taken into account already when filling out this chunk
    Real(dp), save, allocatable :: Ant_pos(:,:,:)
    Complex(dp), save, allocatable, target :: CTime_spectr(:,:,:)
    Integer, save, allocatable :: HRefAnt(:,:)
    Complex(dp), save, allocatable, target :: CTime_Hspectr(:,:,:)
    Integer, parameter :: ExcludedStat_max=25
    integer :: N_Chunk_max=1    !  Make this dynamic, a pretty large number
    Integer, save :: Nr_UniqueStat=0
    Integer, save :: Nr_UniqueAnt=0
    Integer, save :: Unique_StatID(0:2*LOFARNrMax), Tot_UniqueAnt(0:2*LOFARNrMax)
    Integer, save, Allocatable :: Unique_SAI(:)      ! highest ranknr in Unique_SAI of unique_antenna for station Unique_StatID(i_stat)
    Integer, save :: Nr_UniqueHStat=0      ! constructed in "Find_unique_StatAnt" in FitParams.f90
    Integer, save :: Nr_UniqueHAnt=0
    Integer, save :: BadAnt_nr, PolFlp_nr
    Integer, parameter :: BadAnt_nr_max=100, PolFlp_nr_max=20, SgnFlp_nr_max=40 ! BadAnt_nr_max=100 on March 19, 2021
    Integer, save :: BadAnt_SAI(BadAnt_nr_max)  ! SAI = Station+Antenna_ID = station_ID*1000 + Ant_nr
    Integer, save :: PolFlp_SAI(PolFlp_nr_max)
    Integer, save :: SignFlp_SAI(SgnFlp_nr_max),SgnFlp_nr
    Integer, save :: ExcludedStatID(ExcludedStat_max)
    Character(len=6), save :: ExcludedStat(ExcludedStat_max)
    Integer, save :: SaturatedSamplesMax  ! Max nr of saturated samples per chunk
    Logical, save :: CalibratedOnly
    Integer, save :: RefStat_ID
    Integer, save :: DataReadError
    Integer, save :: AntennaNrError=0
    Integer, save :: TimeFrame
    Integer, save :: PeaksPerChunk=500  ! actual max number of peaks that will be used (set on input)
    Integer, save, allocatable :: Station_OutOfRange(:) ! number of time cleanPeak finds the trace for this station is insufficient to hold the pulse
    Integer, parameter :: MaxPeaksPerChunk=500             ! Maximal number of peaks per chunck
    !Integer, save :: PeakD_nr                      ! obtained number of peaks while in 'Dual' option
    !Integer, save, allocatable :: PeakSP(:,:)     ! peak Sample nr
    !Integer, save, allocatable :: PeakSWl(:,:)    ! Peak lower width
    !Integer, save, allocatable :: PeakSWu(:,:)    ! Peak upper width
    !Integer, save, allocatable :: PeakSAmp(:,:)   ! Peak Amplitude in ref station
    !Integer, allocatable :: PeakSPw(:,:)          ! local array in 'DualPeakFind'
    !Integer, allocatable :: SPeak(:,:,:)        ! local array in 'DualPeakFind'
    Real(dp), save :: NoiseLevel=80.
    Real*8 :: TimeBase=0.0  ! [ms]; To be added to the time-offset read-in on the first line. Mostly introduced for interferometry.
    !Logical, save :: IntfPhaseCheck=.true. ! Print phase info when performing interference imaging
    Real(dp),save :: NormOdd, NormEven, NormHOdd, NormHEven !Powr_eo(0:1)
    !Integer,save :: NAnt_eo(0:1)
    contains
! --------------
    Subroutine Alloc_Chunk_AntInfo
    Implicit none
    write(2,*) 'Assign array space with (Time_dim, Ant_nrMax, ChunkNr_dim)=', Time_dim, Ant_nrMax, ChunkNr_dim
    Flush(unit=2)
    Allocate( Station_OutOfRange(1:Used_StationNr) )
    Allocate( Ant_Stations(1:Ant_nrMax,1:ChunkNr_dim) )
    Allocate( Ant_IDs(1:Ant_nrMax,1:ChunkNr_dim) )
    Allocate( Ant_RawSourceDist(1:Ant_nrMax ,1:ChunkNr_dim) ) ! distances that have been taken into account already when filling out this chunk
    Allocate( Ant_pos(3,1:Ant_nrMax ,1:ChunkNr_dim) )
    Allocate( CTime_spectr(1:Time_dim,1:LAnt_nrMax,1:ChunkNr_dim) )
    !Allocate( Unique_StatID(1:Used_StationNr), Tot_UniqueAnt(0:Used_StationNr) )
    Allocate( Unique_SAI(1:Ant_nrMax))
    If(HAnt_nrMax.gt.0) Then
      !Allocate( Unique_HStatID(1:HStation_nrMax), Unique_HSAI(1:HAnt_nrMax), Tot_UniqueHAnt(0:HStation_nrMax) )
      Allocate( CTime_Hspectr(1:2*Time_dim,1:HAnt_nrMax,1:ChunkNr_dim) )
      Allocate( HRefAnt(1:ChunkNr_dim,0:1) )
    EndIf
    Allocate( Ant_nr(1:ChunkNr_dim),  LBA_nr(1:ChunkNr_dim),  HBA_nr(1:ChunkNr_dim) )
    Ant_nr(1:ChunkNr_dim)=0
    LBA_nr(1:ChunkNr_dim)=0
    HBA_nr(1:ChunkNr_dim)=0
    Allocate( RefAnt(1:ChunkNr_dim,0:1, 0:1) )
    Allocate( StartT_sam(1:ChunkNr_dim) )
    !Allocate( PeakSP(PeaksPerChunk,0:2) )
    !Allocate( PeakSWl(PeaksPerChunk,0:2) )
    !Allocate( PeakSWu(PeaksPerChunk,0:2) )
    !Allocate( PeakSAmp(PeaksPerChunk,0:2) )
    !Allocate( PeakSPw(PeaksPerChunk,0:1) )
    !Allocate( SPeak(PeaksPerChunk,1:4,0:1) )
    End Subroutine Alloc_Chunk_AntInfo
   !=====================================
   Subroutine LHBASet(i_ant,i_chunk,CTspec_a_c, i_type, Sample_ms, toSample)
      use Constants, only : dp, LBAs_ms, HBAs_ms
!      use Chunk_AntInfo, only : CTime_spectr, CTime_Hspectr, Ant_Stations
      Implicit none
      Integer, intent(in) ::  i_ant, i_chunk
      Complex(dp), Pointer, intent(out) :: CTspec_a_c(:)
      Integer, intent(out) :: i_type, toSample
      Real(dp), intent(out) :: Sample_ms
      !
      i_type=MOD(Ant_Stations(i_ant,i_chunk),10)
      If(i_type.eq. 1) Then
         CTspec_a_c(1:2*Time_dim) => CTime_Hspectr(1:2*Time_dim,i_ant, i_chunk)      ! HBA
         Sample_ms=HBAs_ms
         toSample=2
      Else
         CTspec_a_c(1:Time_dim) => CTime_spectr(1:Time_dim,i_ant, i_chunk)       ! LBA
         Sample_ms=LBAs_ms
         toSample=1
      EndIf
      Return
   End Subroutine LHBASet
    !
End Module Chunk_AntInfo
!=================================
Module ThisSource
    use DataConstants, only : PeakNr_dim
    use DataConstants, only : ChunkNr_dim, Ant_nrMax, Used_StationNr  ! applies to all
    use DataConstants, only : LAnt_nrMax, Used_LStationNr, HAnt_nrMax, Used_HStationNr  ! applies to LBA or HBA only
    use constants, only : dp
    !
    ! Length of 'dead end', [5 ns], at either side to calculate the correlation time trace
    Integer, save :: Safety=20 ! [5 ns].
    Integer, save :: SafeHL(0:1) ! [LBA/HBA samples].
    ! Place holder, set on input by 'FitRange_Samples'
    !
    ! length of reference pulse [5 ns], including Hann window of 2x5
    Integer, parameter :: Tref_dim=20 !  [5 ns]
    Integer :: PulseDimHL(0:1)=(/Tref_dim,2*Tref_dim/)  ! for LBA & HBA
    ! Not (yet) dynamic.
    !
    ! Used for calculating the correlation functions:
    ! The total length of the time trace, sum of two 'dead ends' and pulse length.
    Integer, save :: T2_dim=128 ! in [samples] for H/L, Place holder
    !   T2_dim = Tref_dim + Safety*2  !  Length of trace for calculating cross correlations,
    !             better be multiple of 2, and easies when same length for LBA and HBA (saves FFT-initiations)
    !             Used for LBA: T2_dim = 20 + 2* 20=64, thus for HBA double = 128
    ! The actual value is set in 'Alloc_ThisSource' to the max value, i.e. value for 'i_type=1': HBA
    !
    !Integer, parameter :: NrP=4  ! Max nr of peaks per spectrum during calibration run
    Integer :: NrP=4  ! Max nr of peaks per spectrum during calibration run
    integer, parameter :: UpSamplFact=1
    Integer, save, allocatable :: CorrAntNrs(:,:,:)
    Integer, save, allocatable :: Nr_Corr(:,:)
    Integer, save, allocatable ::  PeakNr(:,:)       ! number of peaks for this (i_eo,i_chunk)
    Integer, save, allocatable ::  TotPeakNr(:,:)    ! last peak# for this (i_eo,i_chunk)
    !
    Integer, save, allocatable  ::  Peak_eo(:), Peak_ty(:), ChunkNr(:) ! Peak_start(:) ! , Peak_Offst(1:PeakNr_dim)
    !Integer, save, allocatable  ::  Peak_RefSAI(:), Peak_IRef(:)
    Real(dp), save, allocatable  :: Peak_Offst(:) ! Set in 'GetCorrSingAnt',
    !          pulse arrival time shift (w.r.t. core CS002) of pulse in reference antenna depending on source position
    Integer, save, allocatable  ::  ExclStatNr(:,:)
    integer, save, allocatable  :: PeakPos(:)  ! peak position in the reference antenna
    integer, save, allocatable  :: PeakCoret(:)  ! peak position in time trace of a virtual antenna at the core
    ! Discussion:
    !    Time at core allows for more straightforward calculation of expected peak-times in other stations, however:
    !    this may no linger correspond to the correct time in the reference station where the peak is found once the source position starts moving in fitting procedure.
    Real(dp), save, allocatable  :: SourcePos(:,:)
    Real(dp), save, allocatable  :: RefAntErr(:)
    Real(dp), save, allocatable  :: PeakRMS(:)
    Real(dp), save, allocatable  :: PeakChiSQ(:)
    Real(dp), save, allocatable  :: CCorr_L(:,:,:)  ! used for plotting only
    Real(dp), save, allocatable  :: CCorr_H(:,:,:)  ! used for plotting only
!    Real(dp), save, allocatable  :: CCorr_pp(:,:,:)
!    Real(dp), save, allocatable  :: CCNorm(:,:)  ! Value of the max of the correlation function
    Real(dp), save, allocatable  :: CCorr_val(:,:)  ! Position of the maximum of the correlation function
    Real(dp), save, allocatable  :: CCorr_tSh(:,:)  ! Position of the maximum of the correlation function
    Real(dp), save :: Error_norm=0.2 ! default minimal timing-error for an antenna, in [LBA samples]
    Logical, save :: ImpImWeights=.false.  ! Implement Weights factors in error estimates for individual peaks and antennas
    Real(dp), save, allocatable  :: CCorr_Err(:,:)  ! Error (inverse weight) used in fitting, min= Error_norm, in [LBA samples]
    Real(dp), save, allocatable :: MeanErr(:) ! effective number of antennas, where each antenna is weighted with the errorbar^(-2)
    Real(dp), save, allocatable  :: T_Offset(:,:) ! Set in 'GetCorrSingAnt',
    Real(dp), save, allocatable  :: T_Offset_ms(:,:) ! Set in 'GetCorrSingAnt',
    !Real(dp), save, allocatable  :: PolBasis(:,:,:)
    !          pulse arrival time difference w.r.t. reference antenna and depends on source position
!    Real(dp), save, allocatable  :: CCPhase_ave(:,:)
    Integer, save, allocatable :: Dropped(:,:)
    Real(dp), save :: StStdDevMax_ns=100.d0  ! Maximum allowed standard deviation in antenna-arraval times for a pulse in one station while calibrating
    Complex(dp), save, allocatable :: CnuRt(:,:,:)
    Real(dp), save, allocatable :: t_CCorr(:)
    !integer, save, allocatable :: Unique_pos(:)  ! Used in Xindx to store the peakpositions to find double peaks
    Real(dp), save, allocatable :: CC_WidRef(:)  ! Estimate of the width of the cross correlation for the ref station, needed for 'GetCorrSingAnt'
!    Real(dp), save, allocatable :: CC_Qual(:,:)  ! Quality of the cross correlation
    Real(dp), save :: CCShapeCut=0.1            ! Shape Quality cut parameter use to cut correlation spectra that differ from the self-correlation
    Real(dp), save :: CCShapeCut_lim=0.6        ! limit set on CCShapeCut
    Real(dp), save :: ChiSq_lim=80.              ! limit set on the chi-square (FitQual)
    Real(dp), save :: EffAntNr_lim=0.8        ! limit set on the ratio of effective number of used antennas v.s. available number
    Real(dp), save :: XFrameEi, XFrameEf, XFrameNi, XFrameNf, XFrameh       ! Frame used in Explore
    !
    Integer, save ::  PeakNrTotal
    Logical, save ::  RealCorrelation=.false.
    Logical, save ::  PrntCCNrm=.false.
    ! PlotCCPhase is disabled for HBA, hardly used for LBA
    !Logical, save ::  PlotCCPhase=.false. ! produce the data for plotting a map of the CrossCorr phase at max
    Integer, save :: CurtainHalfWidth=-1  ! Make curtain plot when positive
    Logical, save ::  Dual=.true.!.false. ! .true.: set source locations identical for even and odd antennas when pulse-position is same in ref. antennas
    Real(dp), save, allocatable  :: Stat_pos(:,:)
!    Real(dp), save :: CC_Wid, CC_Max
    !Integer, save :: FineOffset_StatID(1:Station_nrMax)
    contains
! --------------
    Subroutine Alloc_ThisSource
       Implicit none
       Integer :: i,j, i_type=1
       Allocate( CorrAntNrs(1:Ant_nrMax,0:1,1:ChunkNr_dim))
       Allocate( Nr_Corr(0:1,1:ChunkNr_dim))
       !
       SafeHL(0:1)=(/Safety,2*Safety/)  ! for LBA & HBA
       T2_dim = PulseDimHL(i_type) + SafeHL(i_type)*2  !  Length of trace for calculating cross correlations, better be multiple of 2
       !write(2,*) '! Alloc_ThisSource1:',Safety, T2_dim, Tref_dim
       Do i=1,11 ! construct as pure multiple of 2
         j=ISHFT(T2_dim, -i)
         !write(2,*) i,j,t2_dim
         If(j.eq.0) Then
            j=1
            T2_dim=ISHFT(j, i)
            exit
         EndIf
       EndDo
       !write(2,*) '! Alloc_ThisSource2:',Safety, T2_dim, Tref_dim
       Allocate( t_CCorr(-SafeHL(i_type):SafeHL(i_type)) )
   Do i=-SafeHL(i_type),SafeHL(i_type)
     t_ccorr(i)=i ! time axis needed for spline interpolation of CrossCorrelation
   Enddo
       If(.not. Allocated(PeakNr) ) Then
          Allocate( PeakNr(0:1,1:ChunkNr_dim) )      ! number of peaks for this (i_eo,i_chunk)
          Allocate( TotPeakNr(0:1,1:ChunkNr_dim))    ! last peak# for this (i_eo,i_chunk)
          Allocate( Peak_eo(1:PeakNr_dim), Peak_ty(1:PeakNr_dim), ChunkNr(1:PeakNr_dim)) ! , Peak_Offst(1:PeakNr_dim)
          Allocate( PeakPos(1:PeakNr_dim))  ! peak time [LBA-samples in this chunk] in the reference antenna
          ! could also be Peak_TiRC = Time in Referenceantenna Chunk [LBA samples]
          Allocate( SourcePos(3,1:PeakNr_dim))
          Allocate( RefAntErr(1:PeakNr_dim))
       EndIf
       If(.not. Allocated(ExclStatNr) ) Then
          Allocate( ExclStatNr(1:Used_StationNr,1:PeakNr_dim))
          ExclStatNr(:,:)=0
       EndIf
       !
       Allocate( Stat_pos(1:3,1:Used_StationNr)  )
       !Allocate( Peak_RefSAI(1:PeakNr_dim))
       !Allocate( Peak_IRef(1:PeakNr_dim))
       Allocate( Peak_Offst(1:PeakNr_dim)) ! [LBA-samples]; Set in 'GetCorrSingAnt',
!       Allocate( Peak_Offst_ms(1:PeakNr_dim)) ! [ms]; Set in 'GetCorrSingAnt',
       !          pulse arrival time shift (w.r.t. core CS002) of pulse in reference antenna depending on source position
       Allocate( PeakRMS(1:PeakNr_dim))
       Allocate( PeakChiSQ(1:PeakNr_dim))
       Allocate( MeanErr(1:PeakNr_dim))
       Allocate( CnuRt(0:T2_dim/2,1:3,1:PeakNr_dim) ) ! used in GetCorrSingAnt, but needs set centrally because of 'save'
       Allocate( CCorr_L(-SafeHL(0):SafeHL(0), LAnt_nrMax, 1:PeakNr_dim))
       Allocate( CCorr_H(-SafeHL(1):SafeHL(1), HAnt_nrMax, 1:PeakNr_dim))
!       Allocate( CCorr_pp(-Safety:Safety,Ant_nrMax,1:PeakNr_dim))
!       Allocate( CCNorm(Ant_nrMax,1:PeakNr_dim))  ! Value of the max of the correlation function
       Allocate( CCorr_val(Ant_nrMax,1:PeakNr_dim))  ! Peak value of cross-correlation, normalized by the auto correlations
       Allocate( CCorr_tSh(Ant_nrMax,1:PeakNr_dim))  ! [LBA-samples]; time-shift of the maximum of the correlation function
       Allocate( CCorr_Err(Ant_nrMax,1:PeakNr_dim))  ! [LBA-samples]; Error (inverse weight) used in fitting,
       Allocate( T_Offset(1:Ant_nrMax,1:PeakNr_dim)) ! [LBA-samples]; Set in 'GetCorrSingAnt',
       Allocate( T_Offset_ms(1:Ant_nrMax,1:PeakNr_dim)) ! [ms]; Set in 'GetCorrSingAnt',
       !          pulse arrival time difference w.r.t. reference antenna and depends on source position
!       Allocate( CCPhase_ave!(1:PeakNr_dim,1:Used_StationNr))
       Allocate( Dropped(1:Used_StationNr,1:PeakNr_dim) ) ! Number of antennas effectively dropped from the fitting
       Dropped(:,:)=0
       !Allocate( Unique_pos(1:PeakNr_dim) ) ! Used in Xindx to store the peakpositions to find double peaks
       Allocate( CC_WidRef(1:PeakNr_dim) ) ! [LBA-samples]; Estimate of the width of the cross correlation
       !Allocate( CC_Qual(1:Ant_nrMax,1:PeakNr_dim) ) ! Quality of the cross correlation
       !write(2,"(A,o10,i7)") 'from Alloc_ThisSource, T2_dim:', T2_dim, T2_dim
    End Subroutine Alloc_ThisSource
End Module ThisSource
!=================================
Module FitParams
    use constants, only : dp
    use DataConstants, only : Ant_nrMax, Used_StationNr
    use Chunk_AntInfo, only : LOFARNrMax
    integer, parameter :: N_FitPar_max=1400
    Logical, save :: Fit_AntOffset=.false.  ! .false. ! when .true. fit the off-sets for each antenna in the stations to be fitted
    Logical, save :: FitIncremental=.true.
    Logical, save :: CalcHessian=.false.
    Logical, save :: PulsPosCore=.true.
    Logical, save :: WriteCalib=.false.
    Logical, save :: FullAntFitPrn=.false.
    !Logical, save :: ImagingRun=.False.  ! Also used in RFI-suppression program
    !Logical, save :: Explore=.false.
    Logical, save :: Kalman=.false.
    Logical, save :: FullSourceSearch=.true.
    Logical, save :: MeanCircPol=.false.
    Complex(dp), save, allocatable  :: Cnu01(:,:)
    Integer :: i_SAI
    Integer, save :: Fit_PeakNrTotal  ! total number of sources that are fitted (some may be the same for even and odd antennas)
    Integer, save :: FitParam(N_FitPar_max), N_FitPar, N_FitStatTim, Nr_TimeOffset !, X_Offset(0:Station_nrMax)
    Integer, save ::  X_Offset(0:N_FitPar_max)
    Integer, save :: PeakS_dim=80  ! Number of pulses searched for in the even/odd numbered reference antenna
    Integer, save :: N_EffAnt, Max_EffAnt
    Real(dp), save :: ParamScaleFac(N_FitPar_max)
    Real(dp), save :: SigmaGuess(3)
    Real(dp), save :: Fit_TimeOffsetStat(0:2*LOFARNrMax)
    Real(dp), save, allocatable :: Fit_TimeOffsetAnt(:)
    Real(dp), save :: FitQual    ! chi^2 value which is optimized
    Real(dp), save :: Sigma(4)    ! SQRT(Diagonal) of the covariance matrix
    Real(dp), save :: SpaceCov(1:3,1:3)=reshape((/ -1.D-10,0.d0,0.d0,0.d0,-1.D-10,0.d0,0.d0,0.d0,-1.D-10 /), (/3,3/))    ! The space components of the covariance matrix [m^2]
    Real(dp), save :: Sigma_AntT=3.  ! Intrinsic timing accuracy for pulse-finding in a single antenna [sample]
    Real(dp), save :: SearchRangeFallOff=4.
    Real(dp), save :: AntennaRange=100.  ! Maximum distance (in km) measured from the core for including stations in the fit
    !Integer, save :: MaxFitAntD_nr=11
    !Real(dp), save :: MaxFitAntDistcs(12)=(/0.05, 0.30, 0.6, 1.2, 2.5, 5., 10., 20., 30., 50., 100., -1. /)
    Integer, save :: MaxFitAntD_nr=12
    Real(dp), save :: MaxFitAntDistcs(12)=(/0.05, 0.30, 0.6, 1.0, 2.0, 5., 10., 20., 30., 60., 80., 100. /)
    ! Maximum distance from the core for the stations to be included in the fit for source location
    contains
! --------------
    Subroutine Alloc_FitParams
       Implicit none
       integer, save :: Previous_Used_StationNr=0, Previous_Ant_nrMax=0
      !If(Previous_Used_StationNr .lt. Used_StationNr) Then
      !   If(Allocated(Fit_TimeOffsetStat) ) DeAllocate( Fit_TimeOffsetStat )
      !   Previous_Used_StationNr = Used_StationNr
      !   Allocate( Fit_TimeOffsetStat(1:Used_StationNr) )
      !EndIf
      If(Previous_Ant_nrMax .lt. Ant_nrMax) Then
         If(Allocated(Fit_TimeOffsetAnt) ) DeAllocate( Fit_TimeOffsetAnt )
         Previous_Ant_nrMax = Ant_nrMax
         Allocate( Fit_TimeOffsetAnt(1:Ant_nrMax) )
      EndIf
       Fit_TimeOffsetStat(:)=0.
       Fit_TimeOffsetAnt(:)=0.
    End Subroutine Alloc_FitParams
End Module FitParams
!=================================
!======================
Module Explore_Pars
    use constants, only : dp
    use DataConstants, only : Ant_nrMax
    Integer, save, allocatable :: StatsStore_SAI(:,:)
    Real(dp), save, allocatable :: StatsStore_Ave(:,:)
    Real(dp), save, allocatable :: StatsStore_RMS(:,:)
    Real(dp), save, allocatable :: StatsStore_Peak(:,:,:)
    Integer, parameter :: PeakS=4
    Integer, save :: N_ExplTimes
    Real(dp), save :: NMin=+99., NMax=-99., Emin=+99., EMax=-99. ! Extreme Northing and Easting [km]
    contains
! --------------
    Subroutine Alloc_Explore_Pars
       Implicit none
       Allocate( StatsStore_SAI(1:N_ExplTimes,1:Ant_nrMax) )
       Allocate( StatsStore_Ave(1:N_ExplTimes,1:Ant_nrMax) )
       Allocate( StatsStore_RMS(1:N_ExplTimes,1:Ant_nrMax) )
       Allocate( StatsStore_Peak(1:N_ExplTimes,1:Ant_nrMax,1:PeakS) )
    End Subroutine Alloc_Explore_Pars
    !
End Module Explore_Pars
!=================================
