!=================================
Module Interferom_Pars
   use constants, only : dp
   use DataConstants, only : Time_Dim  !, Cnu_dim, DataFolder, OutFileLabel
   use DataConstants, only : Ant_nrMax, PeakNr_dim, ChunkNr_dim
   Implicit none
   Complex(dp) :: CTime_sum(1:Time_dim)  ! should indeed be the large dimension
   integer, parameter :: Nmax_IntFer=250  ! Max number of interference antennas
   integer, allocatable, save :: Nr_IntFerCh(:)   ! nr of antenna pairs included per chunk
   integer, save :: Nr_IntFerMx=0 ! maximum over all chunks of nr of antenna pairs included, set in EISelectAntennas
   integer, save :: Nr_IntFer=0 ! set equal to Nr_IntFerMx, to become obsolete
   Integer, allocatable, save :: IntFer_ant(:,:)
   !Integer, Save ::  IntfSmoothWin=20   ! set equal to N_Smth
   Integer :: i_chunk
   Real(dp), save :: dtAnt_tp(0:Ant_nrMax/2)  ! Step size for pixels; obsolete but used in CrossCorr
   Integer, save :: NSrc_tp(0:Ant_nrMax/2)  !! difference in timing for theta and phi pol; Obsolete
   !Real(dp), allocatable :: IntfNorm(:)
   Real(dp), save :: d_loc(1:3)  ! Step size for pixels
   Real(dp), save :: diff(3) ! max abs. time shift per pixel in direction i
   integer, save :: N_pix(3,2) ! Min & maxvalue of range
   !Real(dp), save :: StartTime_ms ! replaced by  StartTime_ms=StartT_sam(i_chunk)*sample*1000.d0
   Real(dp), save :: PixLoc(1:3)
   Real(dp), save :: CenLoc(1:3), CenLocPol(1:3)
   !Complex(dp), allocatable :: Cnu(:,:), Cnu_pix(:), CTime_pix(:)
   !Real(dp), allocatable :: RTime(:)
   !Integer, parameter :: NS_max=20 ! max number of slices in the SumWindw
   Real(dp), allocatable :: AveInten(:,:), AveIntenE(:,:), AveIntenN(:,:), SlcInten(:), RimInten(:,:)
   Real(dp), allocatable :: MaxSlcInten(:), MaxSlcIntenLoc(:,:)
   Real(dp), save :: t_shft, t_offsetPow !, PowerScale=1.d0
   Integer, save :: N_sum, SumStrt, SumWindw, NrSlices, SliceLen
   !     In many cases  NrSlices=1 and  SliceLen=SumWindw  use NrPixSmPowTr ans N_smth instead
   Integer, save :: IntfDim, IntfNuDim, IntfLead, IntfBase
   Real(dp), save :: MaxIntfInten, MaxIntfIntenLoc(1:3)
   Real(dp), save :: NewCenLoc(1:3)=(/0.d0,0.d0,0.d0/)
   Integer,save :: ChainRun=0
   Logical, save :: Polar=.false. ! .true.
   !
   Real(dp), allocatable :: RefSmPowTr(:)
   !   Real(dp), allocatable :: PixelPower(:), MaxSmPow(:), MaxPowPix(:,:), PixSmPowTr(:,:,:,:)  ! MaxPowPix(:,i)==PixLoc(:)
   Real(dp), allocatable :: PixelPower(:), MaxSmPow(:), MaxSmPowLoc(:,:)  !
   Real(dp), allocatable :: MaxSmPowI12(:), MaxSmPowQ(:), MaxSmPowU(:), MaxSmPowV(:)
   Real(dp), allocatable ::  MaxSmPowI3(:), MaxSmPowU1(:), MaxSmPowV1(:), MaxSmPowU2(:), MaxSmPowV2(:)
   Real, allocatable :: PixSmPowTr(:,:,:,:)  ! Array can be huge, and single precision should be sufficient
   Real(dp), save :: RMSGridInterpol(1:3)=(/0.d0,0.d0,0.d0/)  ! to keep statistics on shift from gridpoint.
   Integer, save :: NGridInterpol=0  ! the number of sources used for calculating RMSGridInterpol
   Integer, save :: NGridExtrapol=0  ! Number pf times the grid interpolation goes beyond the set limit of 0.75 (in 3D grid spacing) and thus looks like extrapolate
   Integer :: N_fit ! sets window used for calibration fitting
   Integer :: N_smth=40  ! 40 gives much better localization than 20, less radial scatter (checked 5 Febr 2021)
   Real(dp), allocatable, save :: smooth(:)
   logical, save :: ParabolicSmooth=.true.
   Logical, save :: IntPowSpec=.true.
   Integer, save :: NrPixSmPowTr
   Integer, allocatable :: MaxSmPowGrd(:,:)
   Real(dp), save :: xMin, xMax, yMin, yMax, zMin, zMax, tMin, tMax, AmpltPlot
   Logical, save :: IntfPhaseCheck=.true. ! Print phase info when performing interference imaging
   Real(dp), allocatable :: RimSmPow(:)
   Real(dp), save :: RatMax=1./1.2 ! parameters used in barycentric analysis
   Real(dp), save :: alpha(3), PolBasis(3,3)  ! polbasis(:,k) belongs to eigenvalue alpha(k)
   Integer, save :: PixPowOpt=0   ! used in "EIAnalyzePixelTTrace"
   Integer, save :: N_best=6
   Integer, allocatable :: Grd_best(:,:,:)
   Real(dp), allocatable :: Int_best(:,:)
   Real(dp),save :: StI, StI12, StQ, StU, StV, StI3, StU1, StV1, StU2, StV2, P_un, P_lin, P_circ
   Real(dp),save :: dStI, dStI12, dStQ, dStU, dStV, dStI3, dStU1, dStV1, dStU2, dStV2, Chi2pDF
   Complex, save :: Stk_NEh(1:3,1:3)  ! Not double; Stokes parameters in Cartesian frame
   Real, save :: PolZen(3), PolAzi(3), PolMag(3), PoldOm(3) ! analyzed polarization direction, from "PolTestCath"
   Logical, save :: RefinePolarizObs ! Calculate polarization observables in cartesian coordinates at interpolated location & write pol. obs.
   Real(dp),save ::  BoundingBox(1:2,1:4)
   !Real(dp) :: EleAng,AziAng  elevation & azimuth
   !
   Real(dp), Allocatable :: D_a(:),Ph_a(:)  ! weighted p and t directions in polarization basis (_PB)
   Real(dp), Allocatable :: wap_PB(:,:),wat_PB(:,:)  ! weighted p and t directions in polarization basis (_PB)
   Complex(dp), Allocatable :: wEnu_p(:), wEnu_t(:)
   Complex(dp), Allocatable :: wETime_p(:), wETime_t(:)
   Complex(dp), Allocatable :: wEtime_ap(:,:), wEtime_at(:,:) ! weighted E fields, time dependent
   !Complex(dp), Allocatable :: dChi_ap(:), dChi_at(:)
   Real(dp), Allocatable :: dChi_ap(:), dChi_at(:)
   Complex(dp), Allocatable, save :: Cnu_p0(:,:,:), Cnu_t0(:,:,:), Cnu_p1(:,:,:), Cnu_t1(:,:,:)
   Real(dp), Allocatable :: TIntens_pol(:,:)   !  Intensity trace for different pol directions as seen at the core, For use in "PeakInterferoOption"
   !       array TIntens_pol is build in subroutine EI_PolGridDel
   Real(dp), Allocatable :: TIntens000(:)   !  Intensity trace for central pixel, For use in "PeakInterferoOption", build in EIAnalyzePixelTTrace

   !Real(dp), Allocatable, save :: W_ap(:,:), W_at(:,:) !created in EI_Weights for all chunks
   !Real(dp), Allocatable, save :: Noise_p(:), Noise_t(:) !needed in EI_Weights for present chunk for generating weights
   !Complex(dp), Allocatable, save :: CTime_p(:,:), CTime_t(:,:) !needed in EI_Weights for present chunk for generating weights
   Real(dp), Allocatable, save :: AntPeak_OffSt(:,:)  ! Calculated in EI_PolSetUp
   Real(dp), save :: dnu
   Integer, save :: inu1, inu2
   Real(dp),save :: NormEven, NormOdd
   Logical ,save :: FirstTimeInterf=.true.
   !
   contains
! --------------
! --------------
   Subroutine Alloc_EInterfAll_Pars
      !  Initialize many arrays, call should be followed at some stage by call to "DeAlloc_EInterfAll_Pars"
      use AntFunCconst, only : Freq_min, Freq_max, Gain
      use FFT, only : RFTransform_su
      Logical :: GainFactor=.true.
      Real(dp) :: Nrm
      Logical, save :: AntennaFunctionDone=.false.
      Integer, save :: N_smth_previous=-1
      If(.not. AntennaFunctionDone) Then
         Call AntFieParGen()
         If(.not. GainFactor) Then
            Nrm=SUM(Gain(:))/(Freq_max - Freq_min)
           write(2,*) 'E-field traces do not contain an additional frequency-dependent Gain factor, average=',NRM
            Gain(:)=Nrm
         EndIf
         AntennaFunctionDone=.true.
      EndIf
      IntfDim=2*IntfNuDim  !  should be about =512 and may be different on subsequent calls
      dnu=100./IntfNuDim   ! [MHz] Jones matrix is stored on 1MHz grid
      inu1=Int(Freq_min/dnu)+1
      inu2=Int(Freq_max/dnu)-1
      Nr_IntFer=Nr_IntFerMx
      !Allocate(  CTime_p(1:2*IntfNuDim,1:Nr_IntFerMx), CTime_t(1:2*IntfNuDim,1:Nr_IntFerMx) )
      !Allocate(  Noise_p(1:Nr_IntFerMx), Noise_t(1:Nr_IntFerMx) )
      If(.not. Allocated(AntPeak_OffSt) ) Allocate(  AntPeak_OffSt(1:Ant_nrMax,1:PeakNr_dim) )  ! keep for possible later use
      !Allocate( NEh2PB(1:3,1:3,1:PeakNr_dim))
      !write(2,*) '!N_smth_previous, N_smth:', N_smth_previous, N_smth
      If(N_smth_previous .ne. N_smth) Then
         If(Allocated(Smooth) ) DeAllocate( Smooth )
         N_smth_previous = N_smth
         If(Allocated(TIntens_pol) ) DeAllocate( TIntens_pol )
         If(Allocated(TIntens000) ) DeAllocate( TIntens000 )
         Allocate( Smooth(-N_smth:N_smth) )
         Allocate( TIntens000(1:SumWindw) )  !  Intensity trace for different pol directions as seen at the core, For use in "PeakInterferoOption"
         Allocate( TIntens_pol(1:3,-N_smth:N_smth) )  !  Intensity trace for different pol directions as seen at the core, For use in "PeakInterferoOption"
         Call SetSmooth(N_Smth, ParabolicSmooth, Smooth)
         Write(2,"(1x,A,I3,A,I3,I3,A,3F6.3)") 'Width smoothing fie=', N_smth, &
            '; relative values in range (',N_smth/2-1,N_smth/2+1,') = ',Smooth(N_smth/2-1:N_smth/2+1)/Smooth(0)
      EndIf
      Call RFTransform_su(IntfDim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   End Subroutine Alloc_EInterfAll_Pars
! --------------
   Subroutine Alloc_EInterfCalib_Pars
      !use FFT, only : RFTransform_su
      IntfNuDim=512  !
      Allocate( Cnu_p0(0:IntfNuDim,1:Nr_IntFerMx,1:PeakNr_dim) )   !  indices interchanged !!!
      Allocate( Cnu_p1(0:IntfNuDim,1:Nr_IntFerMx,1:PeakNr_dim) )   !  indices interchanged !!!
      Allocate( Cnu_t0(0:IntfNuDim,1:Nr_IntFerMx,1:PeakNr_dim) )   !  indices interchanged !!!
      Allocate( Cnu_t1(0:IntfNuDim,1:Nr_IntFerMx,1:PeakNr_dim) )   !  indices interchanged !!!
      Call Alloc_EInterfAll_Pars
      !write(2,*) 'Allocate_InterferometricFieldCalibration_Parameters:',IntfNuDim,Nr_IntFerMx,IntfDim,NrSlices
      !
   End Subroutine Alloc_EInterfCalib_Pars
! --------------
   Subroutine Alloc_EInterfImag_Pars
      !  Initialize many arrays, call should be followed at some stage by call to "DeAlloc_EInterfImag_Pars"
      Integer, save :: NrPixSmPowTr_Previous=-1
      Allocate( Cnu_p0(0:IntfNuDim,1:Nr_IntFerMx,1) )   !  indices interchanged !!!
      Allocate( Cnu_p1(0:IntfNuDim,1:Nr_IntFerMx,1) )   !  indices interchanged !!!
      Allocate( Cnu_t0(0:IntfNuDim,1:Nr_IntFerMx,1) )   !  indices interchanged !!!
      Allocate( Cnu_t1(0:IntfNuDim,1:Nr_IntFerMx,1) )   !  indices interchanged !!!
      Call Alloc_EInterfAll_Pars
      !write(2,*) 'Alloc_EInterfImag_Pars:', SumWindw, N_smth, NrSlices,N_pix(3,1), N_pix(3,2)
      !Flush(unit=2)
      NrPixSmPowTr=(SumWindw-1)/N_smth-1  ! smallest window has length=(2*N_smth +1)
      N_fit=-1
      If(NrPixSmPowTr.lt.1) then
         NrPixSmPowTr=0
         IntPowSpec=.false.
      Endif
      !write(2,*) 'Alloc_Interferom_Pars:',IntfNuDim,Nr_IntFerMx,IntfDim,NrSlices, NrPixSmPowTr, N_smth
      !     In many cases  NrSlices=1 and  SliceLen=SumWindw  use NrPixSmPowTr and N_smth instead
      Allocate( AveInten(1:NrSlices,N_pix(3,1):N_pix(3,2)) )
      Allocate( AveIntenE(1:NrSlices,N_pix(3,1):N_pix(3,2)) )
      Allocate( AveIntenN(1:NrSlices,N_pix(3,1):N_pix(3,2)))
      Allocate( RimInten(1:NrSlices,N_pix(3,1):N_pix(3,2)) )
      Allocate( SlcInten(1:NrSlices), MaxSlcInten(1:NrSlices), MaxSlcIntenLoc(1:NrSlices,3))
      !
      Allocate( MaxSmPowI12(0:NrPixSmPowTr) ) !, MaxSmPowPix(3,1:NrPixSmPowTr) )
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
      If(FirstTimeInterf) write(2,"('Storing pixel traces takes ',F10.6,' Gbytes')")  &
            NrPixSmPowTr* (N_pix(1,2)-N_pix(1,1)+1.)* (N_pix(2,2)-N_pix(2,1)+1.)* (N_pix(3,2)-N_pix(3,1)+1.)*4./1.E9
      allocate(  RimSmPow(N_pix(3,1):N_pix(3,2)) )
      If(NrPixSmPowTr_Previous .ne. NrPixSmPowTr) Then
         write(2,*) '!Reallocate arrays MaxSmPow..', NrPixSmPowTr_Previous, NrPixSmPowTr
         If( Allocated( MaxSmPowGrd)) DeAllocate( MaxSmPowGrd )
         Allocate( MaxSmPowGrd(1:3,0:NrPixSmPowTr) )  ! keep for possible later use
         If( Allocated( MaxSmPow)) DeAllocate( MaxSmPow)
         Allocate( MaxSmPow(0:NrPixSmPowTr) )  ! keep for possible later use
         If( Allocated( MaxSmPowLoc)) DeAllocate( MaxSmPowLoc)
         Allocate( MaxSmPowLoc(1:3,0:NrPixSmPowTr) )  ! keep for possible later use
         If( Allocated( Grd_best)) DeAllocate( Grd_best)
         Allocate( Grd_best(1:3,1:N_best,0:NrPixSmPowTr) )  ! keep for possible later use
         If( Allocated( Int_best)) DeAllocate( Int_best)
         Allocate( Int_best(1:N_best,0:NrPixSmPowTr) )  ! keep for possible later use
      EndIf
   End Subroutine Alloc_EInterfImag_Pars
! -----------------------------
   Subroutine DeAlloc_EInterfImag_Pars
      Call DeAlloc_EInterfAll_Pars
      DeAllocate( AveInten, AveIntenE, AveIntenN)
      DeAllocate( RimInten, SlcInten, MaxSlcInten, MaxSlcIntenLoc)
      DeAllocate( MaxSmPowI12, MaxSmPowQ, MaxSmPowU, MaxSmPowV, MaxSmPowI3, MaxSmPowU1, MaxSmPowV1, MaxSmPowU2, MaxSmPowV2 )
      DeAllocate( RefSmPowTr, PixSmPowTr, MaxSmPowGrd ) ! ,MaxSmPow,  MaxSmPowPix
      DeAllocate( RimSmPow )
   End Subroutine DeAlloc_EInterfImag_Pars
! ------------------------------------------------
   Subroutine DeAlloc_EInterfAll_Pars
      use FFT, only : DAssignFFT
      Call DAssignFFT()                   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DeAllocate( Cnu_p0, Cnu_t0, Cnu_p1, Cnu_t1 )
   End Subroutine DeAlloc_EInterfAll_Pars
End Module Interferom_Pars
