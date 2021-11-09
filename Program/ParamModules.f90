!======================
Module Chunk_AntInfo
    use constants, only : dp
    use DataConstants, only : Ant_nrMax, Time_dim,Station_nrMax, ChunkNr_dim
    Character(len=5), save :: Station_name(1:Station_nrMax)=(/&
        "CS001","CS002","CS003","CS004","CS005","CS006","CS007","CS011","CS013","CS017",&
        "CS021","CS024","CS026","CS030","CS031","CS032","CS101","CS103","CS201","CS301",&
        "CS401","CS501","RS106","RS205","RS208","RS210","RS305","RS306","RS307","RS310",&
        "RS406","RS407","RS409","RS503","RS508","CS302","CS028","RS509","RS   ","XS   "  /)
    Integer, save :: Station_number(1:Station_nrMax)=&
        (/   1 ,     2 ,     3 ,     4 ,     5 ,     6 ,     7 ,    11 ,    13 ,    17 ,&
            21 ,    24 ,    26 ,    30 ,    31 ,    32 ,   101 ,   103 ,   121 ,   141 ,&
           161 ,   181 ,   106 ,   125 ,   128 ,   130 ,   145 ,   146 ,   147 ,   150 ,&
           166 ,   167 ,   169 ,   183 ,   188 ,   142 ,    28 ,   189 ,     0 ,     0   /)
    Character(len=80), save :: Simulation
    Integer, save :: WriteSimulation(2)=(/ -1,-1 /)  ! start sample, window size
    Integer, save, allocatable :: Ant_Stations(:,:)
    Integer, save, allocatable :: Ant_IDs(:,:)
    Integer, save, allocatable :: Ant_nr(:)
    Integer, save, allocatable :: RefAnt(:,:)
    Integer, save, allocatable :: Start_time(:)
    Real(dp), save, allocatable :: Ant_RawSourceDist(:,:) ! distances that have been taken into account already when filling out this chunk
    Real(dp), save, allocatable :: Ant_pos(:,:,:)
    Complex(dp), save, allocatable :: CTime_spectr(:,:,:)
    Integer, parameter :: ExcludedStat_max=35
    Integer, save :: Unique_StatID(1:Station_nrMax), Nr_UniqueStat=0
    Integer, save :: Unique_SAI(1:Ant_nrMax), Nr_UniqueAnt=0
    Integer, save :: Tot_UniqueAnt(0:Station_nrMax)      ! highest ranknr in Unique_SAI of unique_antenna for station Unique_StatID(i_stat)
    Integer, save :: BadAnt_nr, PolFlp_nr
    Integer, parameter :: BadAnt_nr_max=100, PolFlp_nr_max=20, SgnFlp_nr_max=40 ! BadAnt_nr_max=100 on March 19, 2021
    Integer, save :: BadAnt_SAI(BadAnt_nr_max)  ! SAI = Station+Antenna_ID = station_ID*1000 + Ant_nr
    Integer, save :: PolFlp_SAI(PolFlp_nr_max)
    Integer, save :: SignFlp_SAI(SgnFlp_nr_max),SgnFlp_nr
    Integer, save :: ExcludedStatID(ExcludedStat_max)
    Character(len=5), save :: ExcludedStat(ExcludedStat_max)
    Integer, save :: RefStat_ID
    Integer, save :: DataReadError
    Integer, save :: AntennaNrError=-1
    Integer, save :: TimeFrame
    Integer, save :: PeaksPerChunk=500  ! actual max number of peaks that will be used (set on input)
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
    Real(dp),save :: Powr_eo(0:1)
    Integer,save :: NAnt_eo(0:1)
    contains
! --------------
    Subroutine Alloc_Chunk_AntInfo
    Implicit none
    Allocate( Ant_Stations(1:Ant_nrMax,1:ChunkNr_dim) )
    Allocate( Ant_IDs(1:Ant_nrMax,1:ChunkNr_dim) )
    Allocate( Ant_nr(1:ChunkNr_dim) )
    Allocate( RefAnt(1:ChunkNr_dim,0:1) )
    Allocate( Start_time(1:ChunkNr_dim) )
    Allocate( Ant_RawSourceDist(1:Ant_nrMax,1:ChunkNr_dim) ) ! distances that have been taken into account already when filling out this chunk
    Allocate( Ant_pos(3,1:Ant_nrMax,1:ChunkNr_dim) )
    Allocate( CTime_spectr(1:Time_dim,1:Ant_nrMax,1:ChunkNr_dim) )
    !Allocate( PeakSP(PeaksPerChunk,0:2) )
    !Allocate( PeakSWl(PeaksPerChunk,0:2) )
    !Allocate( PeakSWu(PeaksPerChunk,0:2) )
    !Allocate( PeakSAmp(PeaksPerChunk,0:2) )
    !Allocate( PeakSPw(PeaksPerChunk,0:1) )
    !Allocate( SPeak(PeaksPerChunk,1:4,0:1) )
    End Subroutine Alloc_Chunk_AntInfo
    !
End Module Chunk_AntInfo
!=================================
Module ThisSource
    use DataConstants, only : Ant_nrMax, Station_nrMax,PeakNr_dim, ChunkNr_dim
    use constants, only : dp
    ! to do: make Safety & Tref_dim dynamic
    ! to do: make PeakNr_dim & ChunkNr_dim dynamic
    Integer, save :: Safety=20 ! 20 40     ! length of time spectrum used for correlation = 2xSafety + Tref_dim
    Integer, parameter :: Tref_dim=20 ! 40  ! length of reference pulse, including Hann window of 2x5
    Integer, parameter :: NrP=4  ! Max nr of peaks per spectrum during calibration run
    integer, parameter :: UpSamplFact=1
    integer, save :: T2_dim ! = Tref_dim + Safety*2
    Integer, save, allocatable :: CorrAntNrs(:,:,:)
    Integer, save, allocatable :: Nr_Corr(:,:)
    Integer, save, allocatable ::  PeakNr(:,:)       ! number of peaks for this (i_eo,i_chunk)
    Integer, save, allocatable ::  TotPeakNr(:,:)    ! last peak# for this (i_eo,i_chunk)
    !
    Integer, save, allocatable  ::  Peak_eo(:), Peak_start(:) ! , Peak_Offst(1:PeakNr_dim)
    Integer, save, allocatable  ::  Peak_RefSAI(:), Peak_IRef(:)
    Real(dp), save, allocatable  :: Peak_Offst(:) ! Set in 'GetCorrSingAnt',
    !          pulse arrival time shift (w.r.t. core CS002) of pulse in reference antenna depending on source position
    Integer, save, allocatable  ::  ExclStatNr(:,:)
    integer, save, allocatable  :: PeakPos(:)  ! peak position in the reference antenna
    Real(dp), save, allocatable  :: SourcePos(:,:)
    Real(dp), save, allocatable  :: RefAntErr(:)
    Real(dp), save, allocatable  :: PeakRMS(:)
    Real(dp), save, allocatable  :: PeakChiSQ(:)
    Real(dp), save, allocatable  :: CCorr(:,:,:)
    Real(dp), save, allocatable  :: CCorr_pp(:,:,:)
    Real(dp), save, allocatable  :: CCNorm(:,:)  ! Value of the max of the correlation function
    Real(dp), save, allocatable  :: CCorr_max(:,:)  ! Position of the maximum of the correlation function
    Real(dp), save, allocatable  :: CCorr_Err(:,:)  ! Error (inverse weight) used in fitting, default=1 ns
    Real(dp), save, allocatable  :: T_Offset(:,:) ! Set in 'GetCorrSingAnt',
    Real(dp), save, allocatable  :: PolBasis(:,:,:)
    !          pulse arrival time difference w.r.t. reference antenna and depends on source position
    Real(dp), save, allocatable  :: CCPhase_ave(:,:)
    Integer, save, allocatable :: Dropped(:,:)
    Complex(dp), save, allocatable :: CnuRt(:,:,:)
    Real(dp), save, allocatable :: t_CCorr(:)
    integer, save, allocatable :: Unique_pos(:)  ! Used in Xindx to store the peakpositions to find double peaks
    Real(dp), save, allocatable :: CC_WidRef(:)  ! Estimate of the width of the cross correlation
    Real(dp), save, allocatable :: CC_Qual(:,:)  ! Quality of the cross correlation
    Real(dp), save :: CCShapeCut=0.1            ! Shape Quality cut parameter use to cut correlation spectra that differ from the self-correlation
    Real(dp), save :: CCShapeCut_lim=0.6        ! limit set on CCShapeCut
    Real(dp), save :: ChiSq_lim=80.              ! limit set on the chi-square (FitQual)
    Real(dp), save :: EffAntNr_lim=0.8        ! limit set on the ratio of effective number of used antennas v.s. available number
    Real(dp), save :: XFrameEi, XFrameEf, XFrameNi, XFrameNf, XFrameh       ! Frame used in Explore
    !
    Integer, save ::  PeakNrTotal
    Logical, save ::  RealCorrelation=.false.
    Logical, save ::  PrntCCNrm=.false.
    Logical, save ::  PlotCCPhase=.false. ! produce the data for plotting a map of the CrossCorr phase at max
    Logical, save :: CurtainPlot=.false.  ! Make curtain plot
    Logical, save ::  Dual=.true.!.false. ! .true.: set source locations identical for even and odd antennas when pulse-position is same in ref. antennas
    !Real(dp), save :: t_CCorr(-Safety:Safety)
    Real(dp), save :: Stat_pos(1:3,1:Station_nrMax)
    Real(dp), save :: CC_Wid, CC_Max, CC_Int
    !Integer, save :: FineOffset_StatID(1:Station_nrMax)
    contains
! --------------
    Subroutine Alloc_ThisSource
       Implicit none
       T2_dim = Tref_dim + Safety*2
       Allocate( CorrAntNrs(1:Ant_nrMax,0:1,1:ChunkNr_dim))
       Allocate( Nr_Corr(0:1,1:ChunkNr_dim))
       Allocate( PeakNr(0:1,1:ChunkNr_dim) )      ! number of peaks for this (i_eo,i_chunk)
       Allocate( TotPeakNr(0:1,1:ChunkNr_dim))    ! last peak# for this (i_eo,i_chunk)
       !
       Allocate( Peak_eo(1:PeakNr_dim), Peak_start(1:PeakNr_dim)) ! , Peak_Offst(1:PeakNr_dim)
       Allocate( Peak_RefSAI(1:PeakNr_dim))
       Allocate( Peak_IRef(1:PeakNr_dim))
       Allocate( Peak_Offst(1:PeakNr_dim)) ! Set in 'GetCorrSingAnt',
       !          pulse arrival time shift (w.r.t. core CS002) of pulse in reference antenna depending on source position
       Allocate( ExclStatNr(1:10,1:PeakNr_dim))
       Allocate( PeakPos(1:PeakNr_dim))  ! peak position in the reference antenna
       Allocate( SourcePos(3,1:PeakNr_dim))
       Allocate( RefAntErr(1:PeakNr_dim))
       Allocate( PeakRMS(1:PeakNr_dim))
       Allocate( PeakChiSQ(1:PeakNr_dim))
       Allocate( CnuRt(0:T2_dim/2,1:3,1:PeakNr_dim) )
       Allocate( CCorr(-Safety:Safety,Ant_nrMax,1:PeakNr_dim))
       Allocate( CCorr_pp(-Safety:Safety,Ant_nrMax,1:PeakNr_dim))
       Allocate( CCNorm(Ant_nrMax,1:PeakNr_dim))  ! Value of the max of the correlation function
       Allocate( CCorr_max(Ant_nrMax,1:PeakNr_dim))  ! Position of the maximum of the correlation function
       Allocate( CCorr_Err(Ant_nrMax,1:PeakNr_dim))  ! Error (inverse weight) used in fitting, default=1 ns
       Allocate( T_Offset(1:Ant_nrMax,1:PeakNr_dim)) ! Set in 'GetCorrSingAnt',
       !          pulse arrival time difference w.r.t. reference antenna and depends on source position
       Allocate( CCPhase_ave(1:PeakNr_dim,1:Station_nrMax))
       Allocate( Dropped(1:Station_nrMax,1:PeakNr_dim) ) ! Number of antennas effectively dropped from the fitting
       Allocate( t_CCorr(-Safety:Safety) )
       Allocate( Unique_pos(1:PeakNr_dim) ) ! Used in Xindx to store the peakpositions to find double peaks
       Allocate( CC_WidRef(1:PeakNr_dim) ) ! Estimate of the width of the cross correlation
       Allocate( CC_Qual(1:Ant_nrMax,1:PeakNr_dim) ) ! Quality of the cross correlation
       ExclStatNr(:,:)=0
       write(2,"(A,o10,i7)") 'from Alloc_ThisSource, T2_dim:', T2_dim, T2_dim
    End Subroutine Alloc_ThisSource
End Module ThisSource
!=================================
Module FitParams
    use constants, only : dp
    use DataConstants, only : Ant_nrMax, Station_nrMax
    integer, parameter :: N_FitPar_max=200
    Logical, save :: Fit_AntOffset=.false.  ! .false. ! when .true. fit the off-sets for each antenna in the stations to be fitted
    Logical, save :: FitIncremental=.true.
    Logical, save :: CalcHessian=.false.
    Logical, save :: PulsPosCore=.true.
    Logical, save :: WriteCalib=.false.
    Logical, save :: FullAntFitPrn=.false.
    Logical, save :: ImagingRun=.False.  ! Also used in RFI-suppression program
    Logical, save :: Explore=.false.
    Logical, save :: Kalman=.false.
    Logical, save :: FullSourceSearch=.true.
    Logical, save :: MeanCircPol=.false.
    Complex(dp), save, allocatable  :: Cnu01(:,:)
    Integer :: i_SAI
    Integer, save :: Fit_PeakNrTotal  ! total number of sources that are fitted (some may be the same for even and odd antennas)
    Integer, save :: FitParam(N_FitPar_max), N_FitPar, N_FitStatTim, Nr_TimeOffset, X_Offset(Station_nrMax)
    Integer, save :: PeakS_dim=80  ! Number of pulses searched for in the even/odd numbered reference antenna
    Integer, save :: N_EffAnt, Max_EffAnt
    Real(dp), save :: ParamScaleFac(N_FitPar_max)
    Real(dp), save :: SigmaGuess(3)
    Real(dp), save :: Fit_TimeOffsetStat(1:Station_nrMax), Fit_TimeOffsetAnt(1:Ant_nrMax)
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
Module Interferom_Pars
   use constants, only : dp
   use DataConstants, only : Time_Dim  !, Cnu_dim, DataFolder, OutFileLabel
   use DataConstants, only : Ant_nrMax
   Implicit none
   Complex(dp) :: CTime_sum(1:Time_dim)  ! should indeed be the large dimension
   integer, parameter :: Nmax_IntFer=250  ! Max number of interference antennas
   integer, save :: Nr_IntFer
   Integer, Save ::  IntfSmoothWin=20
   Integer :: IntFer_ant(Nmax_IntFer)
   Integer :: i_chunk
   Real(dp), save :: dtAnt_tp(0:Ant_nrMax/2)  ! Step size for pixels
   Integer, save :: NSrc_tp(1:Ant_nrMax/2)  ! Step size for pixels
   Real(dp), allocatable :: IntfNorm(:)
   Real(dp), save :: d_loc(1:3)  ! Step size for pixels
   Real(dp), save :: diff(3) ! max abs. time shift per pixel in direction i
   integer, save :: N_pix(3,2) ! Min & maxvalue of range
   Real(dp), save :: StartTime_ms
   Real(dp), save :: PixLoc(1:3)
   Real(dp), save :: CenLoc(1:3), CenLocPol(1:3)
   Complex(dp), allocatable :: Cnu(:,:), Cnu_pix(:), CTime_pix(:)
   Real(dp), allocatable :: RTime(:)
   !Integer, parameter :: NS_max=20 ! max number of slices in the SumWindw
   Real(dp), allocatable :: AveInten(:,:), AveIntenE(:,:), AveIntenN(:,:), SlcInten(:), RimInten(:,:)
   Real(dp), allocatable :: MaxSlcInten(:), MaxSlcIntenLoc(:,:)
   Real(dp), save :: t_shft, t_offsetPow, PowerScale=1.d0
   Integer, save :: N_sum, SumStrt, SumWindw, NrSlices, SliceLen
   Integer, save :: IntfDim, IntfNuDim, IntfLead, IntfBase
   Real(dp), save :: MaxIntfInten, MaxIntfIntenLoc(1:3)
   Real(dp), save :: NewCenLoc(1:3)=(/0.d0,0.d0,0.d0/)
   Integer,save :: ChainRun=0
   Logical, save :: Polar=.false. ! .true.
   !
   Real(dp), allocatable :: RefSmPowTr(:)
   !   Real(dp), allocatable :: PixelPower(:), MaxSmPow(:), MaxPowPix(:,:), PixSmPowTr(:,:,:,:)  ! MaxPowPix(:,i)==PixLoc(:)
   Real(dp), allocatable :: PixelPower(:), MaxSmPow(:)  !
   Real(dp), allocatable :: MaxSmPowQ(:), MaxSmPowU(:), MaxSmPowV(:)
   Real(dp), allocatable ::  MaxSmPowI3(:), MaxSmPowU1(:), MaxSmPowV1(:), MaxSmPowU2(:), MaxSmPowV2(:)
   Real, allocatable :: PixSmPowTr(:,:,:,:)  ! Array can be huge, and single precision should be sufficient
   Integer :: N_smth=40  ! 40 gives much better localization than 20, less radial scatter (checked 5 Febr 2021)
   Real(dp), allocatable, save :: smooth(:)
   Logical, save :: IntPowSpec=.true.
   Integer, save :: NrPixSmPowTr
   Integer, allocatable :: MaxSmPowGrd(:,:)
   Real(dp), save :: xMin, xMax, yMin, yMax, zMin, zMax, tMin, tMax, AmpltPlot
   Logical, save :: IntfPhaseCheck=.true. ! Print phase info when performing interference imaging
   Real(dp), allocatable :: RimSmPow(:)
   Real(dp), save :: RatMax=1./1.2 ! parameters used in barycentric analysis
   Real(dp), save :: alpha(3), PolBasis(3,3)  ! polbasis(:,k) belongs to eigenvalue alpha(k)
   Integer, save :: PixPowOpt=0   ! used in "EIAnalyzePixelTTrace"
   Real(dp),save :: StI, StI12, StQ, StU, StV, StI3, StU1, StV1, StU2, StV2
   Real(dp),save :: dStI, dStI12, dStQ, dStU, dStV, dStI3, dStU1, dStV1, dStU2, dStV2, Chi2pDF
   !Real(dp) :: EleAng,AziAng  elevation & azimuth
   !
   contains
! --------------
   Subroutine Alloc_SInterferom_Pars
      use FFT, only : RFTransform_su
      write(2,*) 'Alloc_Interferom_Pars:',IntfNuDim,Nr_IntFer,IntfDim,NrSlices
      Allocate( Cnu(1:Nr_IntFer,0:IntfNuDim) )   !  indices interchanged !!!
      Allocate( IntfNorm(1:Nr_IntFer) )
      Allocate( Cnu_pix(0:IntfNuDim) )
      Allocate( CTime_pix(1:IntfDim) )
      Allocate( Rtime(1:IntfDim) )
      Allocate( AveInten(1:NrSlices,N_pix(3,1):N_pix(3,2)) )
      Allocate( AveIntenE(1:NrSlices,N_pix(3,1):N_pix(3,2)) )
      Allocate( AveIntenN(1:NrSlices,N_pix(3,1):N_pix(3,2)))
      Allocate( RimInten(1:NrSlices,N_pix(3,1):N_pix(3,2)) )
      Allocate( SlcInten(1:NrSlices), MaxSlcInten(1:NrSlices), MaxSlcIntenLoc(1:NrSlices,3))
      Allocate( PixelPower(1:SumWindw) )
      NrPixSmPowTr=(SumWindw-1)/N_smth-1
      If(NrPixSmPowTr.lt.1) then
         NrPixSmPowTr=0
         IntPowSpec=.false.
      Endif
      Allocate( MaxSmPow(0:NrPixSmPowTr) ) !, MaxSmPowPix(3,1:NrPixSmPowTr) )
      Allocate( RefSmPowTr(0:NrPixSmPowTr) )
      Allocate( PixSmPowTr(0:NrPixSmPowTr,N_pix(1,1):N_pix(1,2),N_pix(2,1):N_pix(2,2),N_pix(3,1):N_pix(3,2)) )
      ! since Real(4), single precision, requires 4 bytes.
      ! Integer calculation may cause overflow to the sign bit, or worse, thus convert to real early on
      write(2,"('Storing pixel traces takes ',F10.6,' Gbytes')")  &
            NrPixSmPowTr* (N_pix(1,2)-N_pix(1,1)+1.)* (N_pix(2,2)-N_pix(2,1)+1.)* (N_pix(3,2)-N_pix(3,1)+1.)*4./1.E9
      Allocate( MaxSmPowGrd(1:3,0:NrPixSmPowTr) )
      allocate(  RimSmPow(N_pix(3,1):N_pix(3,2)) )
      Call RFTransform_su(IntfDim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   End Subroutine Alloc_SInterferom_Pars
! --------------
   Subroutine Alloc_EInterferom_Pars
      use FFT, only : RFTransform_su
      write(2,*) 'Alloc_Interferom_Pars:',IntfNuDim,Nr_IntFer,IntfDim,NrSlices
      Allocate( AveInten(1:NrSlices,N_pix(3,1):N_pix(3,2)) )
      Allocate( AveIntenE(1:NrSlices,N_pix(3,1):N_pix(3,2)) )
      Allocate( AveIntenN(1:NrSlices,N_pix(3,1):N_pix(3,2)))
      Allocate( RimInten(1:NrSlices,N_pix(3,1):N_pix(3,2)) )
      Allocate( SlcInten(1:NrSlices), MaxSlcInten(1:NrSlices), MaxSlcIntenLoc(1:NrSlices,3))
      NrPixSmPowTr=(SumWindw-1)/N_smth-1
      If(NrPixSmPowTr.lt.1) then
         NrPixSmPowTr=0
         IntPowSpec=.false.
      Endif
      Allocate( MaxSmPow(0:NrPixSmPowTr) ) !, MaxSmPowPix(3,1:NrPixSmPowTr) )
      Allocate( MaxSmPowQ(0:NrPixSmPowTr) )
      Allocate( MaxSmPowU(0:NrPixSmPowTr) )
      Allocate( MaxSmPowV(0:NrPixSmPowTr) )
      Allocate( MaxSmPowI3(0:NrPixSmPowTr) )
      Allocate( MaxSmPowU1(0:NrPixSmPowTr) )
      Allocate( MaxSmPowV1(0:NrPixSmPowTr) )
      Allocate( MaxSmPowU2(0:NrPixSmPowTr) )
      Allocate( MaxSmPowV2(0:NrPixSmPowTr) )
      Allocate( RefSmPowTr(0:NrPixSmPowTr) )
      Allocate( PixSmPowTr(0:NrPixSmPowTr,N_pix(1,1):N_pix(1,2),N_pix(2,1):N_pix(2,2),N_pix(3,1):N_pix(3,2)) )
      ! since Real(4), single precision, requires 4 bytes.
      ! Integer calculation may cause overflow to the sign bit, or worse, thus convert to real early on
      write(2,"('Storing pixel traces takes ',F10.6,' Gbytes')")  &
            NrPixSmPowTr* (N_pix(1,2)-N_pix(1,1)+1.)* (N_pix(2,2)-N_pix(2,1)+1.)* (N_pix(3,2)-N_pix(3,1)+1.)*4./1.E9
      Allocate( MaxSmPowGrd(1:3,0:NrPixSmPowTr) )
      allocate(  RimSmPow(N_pix(3,1):N_pix(3,2)) )
      Call RFTransform_su(IntfDim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   End Subroutine Alloc_EInterferom_Pars
! -----------------------------
   Subroutine DeAlloc_SInterferom_Pars
      use FFT, only : DAssignFFT
      Call DAssignFFT()                   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DeAllocate( Cnu )
      DeAllocate( IntfNorm )
      DeAllocate( Cnu_pix )
      DeAllocate( CTime_pix )
      DeAllocate( Rtime )
      DeAllocate( AveInten, AveIntenE, AveIntenN)
      DeAllocate( RimInten, SlcInten, MaxSlcInten, MaxSlcIntenLoc)
      DeAllocate( PixelPower )
      DeAllocate( MaxSmPow, RefSmPowTr, PixSmPowTr, MaxSmPowGrd ) ! , MaxSmPowPix
      DeAllocate( RimSmPow )
   End Subroutine DeAlloc_SInterferom_Pars
! -----------------------------
   Subroutine DeAlloc_EInterferom_Pars
      use FFT, only : DAssignFFT
      Call DAssignFFT()                   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DeAllocate( AveInten, AveIntenE, AveIntenN)
      DeAllocate( RimInten, SlcInten, MaxSlcInten, MaxSlcIntenLoc)
      !DeAllocate( PixelPower )
      DeAllocate( MaxSmPow, RefSmPowTr, PixSmPowTr, MaxSmPowGrd ) ! , MaxSmPowPix
      DeAllocate( RimSmPow )
   End Subroutine DeAlloc_EInterferom_Pars
End Module Interferom_Pars
!=================================
