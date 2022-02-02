!=================================
Subroutine EI_Fitter(X)
!   Use method of steepest descent by calling routine PRAXIS
    !use constants, only : dp
!  v13: Introduce a dynamic search window based on the covariance matrix 'SpaceCov'
!        - Window is implemented as a parabolic multiplier normalized to unity and is zero
!           at a certain multiple of the search window. This factor is named 'SearchRangeFallOff',
!           implemented in 'ReImAtMax'.
!  v13: 'SpaceCov' is updated after each fit, either by the new covariance matrix, or by reducing the old one.
    ! v17 : offset in Nr_TimeOffset reduced by 1 in nvar
!
    use constants, only : dp
    use FitParams, only : FitParam, N_FitPar, N_FitStatTim, N_FitStatTim, Fit_PeakNrTotal, Nr_TimeOffset, N_FitPar_max
    use FitParams, only : X_Offset
    use Interferom_Pars, only :  N_Smth, N_fit, Nr_IntFerMx, Nr_IntferCh ! the latter gives # per chunk
    !use DataConstants, only : ChunkNr_dim, Production
    use Interferom_Pars, only : W_ap, W_at, Cnu0, Cnu1, IntfNuDim, CTime_p, CTime_t, Noise_p, Noise_t, AntPeak_OffSt
    use ThisSource, only :  PeakNr, PeakNrTotal, ChunkNr, PeakPos, sourcepos, PeakChiSQ
    Implicit none
    real ( kind = 8 ), intent(inout) :: X(N_FitPar_max)
    integer ( kind = 4 ) :: i, j, k
    integer :: v_dim, i_chunk, error
    integer ( kind = 4 ) :: nvar    ! number of parameters
    integer ( kind = 4 ) :: meqn  ! Number of data points to fit
    integer ( kind = 4 ) iv(60+N_FitPar_max)
    external CompareEI  ! subroutine that compares FitFunction to Data
    external ufparm  ! Dummy external routine
    integer ( kind = 4 ) :: uiparm(1)    ! Not really used
    real ( kind = 8 ) :: urparm(1)    ! Not really used
    real ( kind = 8 ),allocatable :: v(:) !  in NL@SOL: real ( kind = 8 ) v(93 + n*p + 3*n + (p*(3*p+33))/2)
    integer ( kind = 4 ) :: prin, NF
    Integer :: i_peak, IntfBase
    !
    !
!stop
    !  Nr_TimeOffset = the number of antenna/station timings that are to be fitted
    !  PeakNrTotal   = number of pulses to be fitted
    !  N_FitPar      =  ??
    !write(2,*) 'N_FitPar', N_FitPar,X_Offset(N_FitPar)
    nvar=X_Offset(N_FitPar)-1  ! total number of parameters that are fitted
    !Call RFTransform_su(2*IntfNuDim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !
    Meqn=0          ! total number of data-points that enter in the chi^2 search
   N_fit=N_Smth/2+1   ! number of smaples on each side of pulse to be included in fitting
   !
    Do i_Peak=1,PeakNrTotal
      i_chunk=ChunkNr(i_Peak)
      IntfBase= Peakpos(i_Peak) - IntfNuDim
      Meqn = Meqn + Nr_IntFerCh(i_chunk)*2*(1+2*N_fit)     ! number of equations
      Call EI_PolSetUp(Nr_IntFerCh(i_chunk), IntfBase, i_chunk, SourcePos(:,i_Peak), &
         AntPeak_OffSt(1,i_Peak), Cnu0(0,1,i_Peak), Cnu1(0,1,i_Peak), W_ap(1,i_Peak), W_at(1,i_Peak))
      Call EI_Weights(Nr_IntFerCh(i_chunk), IntfNuDim, i_chunk, SourcePos(:,i_Peak), &
         AntPeak_OffSt(1,i_Peak), W_ap(1,i_Peak), W_at(1,i_Peak))
    EndDo
    !
    v_dim=93 + Meqn*nvar + 3*Meqn + nvar*(3*nvar+33)/2
    !write(*,*) 'Meqn=',Meqn,nvar,Nr_Corr(0,1),Nr_Corr(1,1),i,k
    !write(2,*) 'Nr_TimeOffset:',Nr_TimeOffset,N_FitPar,N_FitStatTim,Fit_PeakNrTotal,k
    write(2,*) 'Fitting statistics: nvar=',nvar,', Meqn=',Meqn,', buffer size=',v_dim
    allocate( v(v_dim) )
    NF=-1
    call CompareEI ( meqn, nvar, x, NF, v, uiparm, urparm, ufparm )
    !
    call dfault( iv, v)
    iv(1) = 12 ! 12= do not call dfault again
    iv(14) = 0 ! 1: means print a covariance matrix at the solution.
    iv(15) = 0 ! if = 1 or 2, then a finite difference hessian approximation h is obtained.
               ! if =0 nothing is done
               ! if positive: with step sizes determined using v(delta0=44), a multiplicative factor)
               ! If negative: then only function values are used with step sizes determined using v(dltfdc=40)
    iv(19) = 0 ! controls the number and length of iteration summary lines printed
    iv(21) = 0 !2 ! is the output unit number on which all printing is done.
    iv(22) = 0!1 ! print out the value of x returned (as well as the corresponding gradient and scale vector d).
    iv(23) = 1 ! print summary statistics upon returning.
    iv(24) = 0 ! print the initial x and scale vector d
    v(32) =1.d-4 ! is the relative function convergence tolerance
    v(40) =1.40d-3 ! the step size used when computing the covariance matrix when iv(covreq=15) = -1 or -2, step size = v(dltfdc=40) * max(abs(x(i)), 1/d(i))
   ! v(44) =1.d-3 ! the factor used in choosing the finite difference step size used in computing the covariance matrix when
                !    iv(covreq=15) = 1 or 2, step size = v(delta0=44) * max(abs(x(i)), 1/d(i)) * sign(x(i))
    !
    ! iv(nfcall)... iv(6) is the number of calls so far made on calcr (i.e.,
    !             function evaluations, including those used in computing
    !             the covariance).
    ! iv(mxfcal)... iv(17) gives the maximum number of function evaluations
    !             (calls on calcr, excluding those used to compute the co-
    !             variance matrix) allowed.  if this number does not suf-
    !             fice, then nl2sol returns with iv(1) = 9.  default = 200.
    ! iv(mxiter)... iv(18) gives the maximum number of iterations allowed.
    !             it also indirectly limits the number of gradient evalua-
    !             tions (calls on calcj, excluding those used to compute
    !             the covariance matrix) to iv(mxiter) + 1.  if iv(mxiter)
    !             iterations do not suffice, then nl2sol returns with
    !             iv(1) = 10.  default = 150.
    !    iv(18)=1
   !  iv(26) if (iv(covmat) is positive, then the lower triangle of the covariance matrix is stored rowwise in v starting at
   !             v(iv(covmat)).  if iv(covmat) = 0, then no attempt was made.
    ! iv(nfcov).... iv(40) is the number of calls made on calcr when
    !             trying to compute covariance matrices.
    ! iv(ngcall)... iv(30) is the number of gradient evaluations (calls on
    !             calcj) so far done (including those used for computing
    !             the covariance).
    ! iv(ngcov).... iv(41) is the number of calls made on calcj when
    !             trying to compute covariance matrices.
    ! iv(niter).... iv(31) is the number of iterations performed.
    !
   if((nvar .gt. 0) .and. (meqn .gt. nvar)) then   ! otherwise the system is underdetermined
         !NF=-1
         !call CompareCorrTime ( meqn, nvar, x, NF, v, uiparm, Jacobian, ufparm )
         ! subroutine nl2sno ( n, p, x, calcr, iv, v, uiparm, urparm, ufparm, error )
         !    CALCR must be declared external in the calling program.
         !    It is invoked by
         !      call calcr ( n, p, x, nf, r, uiparm, urparm, ufparm )
         !    When CALCR is called, NF is the invocation count for CALCR.  It is
         ! uiparm... (input) user integer ( kind = 4 ) parameter array passed without
         ! change to calcr and calcj.
         !
         ! urparm... (input) user floating-point parameter array passed without
         !                  change to calcr and calcj.
         !
         ! ufparm... (input) user external subroutine or function passed without
         !                  change to calcr and calcj.
         error=0
         call nl2sno ( meqn, nvar, X, CompareEI, iv, v, uiparm, urparm, ufparm, error )
         !
         If(error.ge.10) then
            write(2,*) 'parameter in "nl2itr" is NaN, error=', error,X(1)
            error=10
            goto 9
         endif
         If(error.gt.0) then
            write(2,*) 'endless loop in NL2SOL truncated'
            goto 9
         endif
         !
         !Write(*,*) 'end nl2sol'
          write(2,"('Result, chi^2/ndf=',F9.2)", ADVANCE='NO') 2*v(10)/(meqn-nvar)
          write(*,"('chi^2/ndf=',F9.2)") 2*v(10)/(meqn-nvar)
          write(2,"(', # of fie calls=',i3,' and # of iterations=',i3)") iv(6),iv(31)
   endif
   NF=-1
   call CompareEI ( meqn, nvar, x, NF, v, uiparm, urparm, ufparm )
   Call EI_PrntFitPars(X)
   !N_fit=-1
   Do i_Peak=1,PeakNrTotal
      Call EI_PolarizPeak(i_Peak)
   Enddo
9  Continue
    !
    Deallocate( v )
    N_fit=-1
    !Call DAssignFFT
    !
    !Call EI_PrntFitPars(X)
    Return
    !
End Subroutine EI_Fitter
!=================================================
!=================================================
Subroutine CompareEI( meqn, nvar, X, nf, R, uiparm, urparm, ufparm )
    use constants, only : dp,sample,c_mps,Refrac
    use DataConstants, only : Station_nrMax, Ant_nrMax
!    use ThisSource, only : Nr_Corr, CCorr_max, CCorr_Err
    use ThisSource, only :  SourcePos, PeakNrTotal, PeakPos, Peak_eo, ChunkNr, PeakRMS, PeakChiSQ
    use Interferom_Pars, only : W_ap, W_at, Cnu0, Cnu1, N_fit, AntPeak_OffSt, Chi2pDF
    Use Interferom_Pars, only : IntfNuDim, IntFer_ant, Nr_IntFerMx, Nr_IntferCh ! the latter gives # per chunk
    use Chunk_AntInfo, only : Ant_Stations, Ant_pos
    use Chunk_AntInfo, only : Unique_StatID, Nr_UniqueStat, Ant_IDs, Unique_SAI, Tot_UniqueAnt
    use FitParams, only : N_FitPar, N_FitStatTim, FitParam, X_Offset
    use FitParams, only : Fit_TimeOffsetStat, Fit_TimeOffsetAnt, Fit_AntOffset, FullAntFitPrn
    use StationMnemonics, only : Station_ID2Mnem
!*****************************************************************************80
!
!  Discussion:
!    Given the value of the vector X, this routine computes the value of the residues R=(t_measured-t_X) used for
!     chi^2 fitting of source positions and antenna timing calibration for LOFAR Lightning Imaging
!  Modified:
!    August 2019
!  Author:
!    Olaf Scholten
!
!  Parameters:
!    Input, integer ( kind = 4 ) MEQN, the number of functions.
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!    Input, real ( kind = 8 ) X(NVAR), the current value of the variables.
!    Input, integer ( kind = 4 ) NF, the number of times the residual routine
!       has been called so far.
!    Output, real ( kind = 8 ) R(MEQN), the residual vector, that is, the
!       value of the functions for the given input value of the variables.
!    Input, integer ( kind = 4 ) UIPARM(*), a user array.
!    Input, real ( kind = 8 ) URPARM(*), a user array.
!    Input, external UFPARM, an external reference to a user subroutine
!       or function.
!   v6: account for noise-error in the pulse of the reference antenna by a constant offset
!   v7: account for timing off-set to be fitted per antenna
!
  implicit none
    real ( kind = 8 ), intent(in) :: X(1:nvar)
  integer ( kind = 4 ), intent(in) :: meqn
  integer ( kind = 4 ), intent(in) :: nvar
  integer ( kind = 4 ), intent(in) :: nf
  real ( kind = 8 ), intent(out) :: R(meqn)
  external ufparm
  integer ( kind = 4 ), intent(in) :: uiparm(*)
  real ( kind = 8 ), intent(out) :: urparm(*)
  !real ( kind = 8 ) :: D1,D2,Rns
  Real(dp) :: FitDelay(1:Ant_nrMax)
  Real(dp) :: DelChi(-N_fit:+N_fit,1:2*Nr_IntFerMx)
    !
    integer :: j, i_SAI, i_eqn, i_stat, i_ant, i_Peak, i_chunk, Station_ID, j_IntFer, Antenna_SAI, Outpt
    Character(len=5) :: Station_Mnem
    Character(len=9) :: Label
    !
    Outpt=0
    if(NF.lt.0) Outpt=1
    !
    Call EIX2Source(X)
    !write(2,*) 'X', X(1:5)
    !
    i_eqn=0
    !if(nvar.eq.1) write(*,*) 'X_p',X_p
    Do i_Peak=1,PeakNrTotal
        !i_eo=Peak_eo(i_peak)
        i_chunk=ChunkNr(i_peak)
        !
         Do j_IntFer=1,Nr_IntferCh(i_chunk)
            i_ant=IntFer_ant(j_IntFer,i_chunk)
            Station_ID = Ant_Stations(i_ant,i_chunk)
            Antenna_SAI= 1000*Ant_Stations(i_ant,i_chunk) + Ant_IDs(i_ant,i_chunk)
            !write(*,*) 'j_corr=',j_corr,Antenna_SAI
            Do i_stat=1, Nr_UniqueStat      ! Get station number from the Unique_StatID list
                If(Unique_StatID(i_stat).eq. Station_ID) exit
            enddo
            !if(nvar.eq.1) write(*,*) j_corr,i_ant, Station_ID, Antenna_SAI
            !if(prn) write(2,*) 'Unique_StatID(i_stat):',Unique_StatID(i_stat),Station_ID,i_stat
            FitDelay(i_ant)=Fit_TimeOffsetStat(i_stat)
            FitDelay(i_ant+1)=Fit_TimeOffsetStat(i_stat)
            If(Fit_AntOffset) then
               Do i_SAI=Tot_UniqueAnt(i_stat-1)+1,Tot_UniqueAnt(i_stat)      ! Get antenna number from the Unique_Antennas list
                   If(Unique_SAI(i_SAI).eq. Antenna_SAI) exit
               enddo
               !If(FullAntFitPrn) write(2,*) 'StatFineOff',i_SAI,Antenna_SAI,StatFineOff,Fit_TimeOffsetAnt(i_SAI)
               FitDelay(i_ant)=FitDelay(i_ant) + Fit_TimeOffsetAnt(i_SAI)
               FitDelay(i_ant+1)=FitDelay(i_ant+1) + Fit_TimeOffsetAnt(i_SAI+1)
            endif
         Enddo
         write(Label,"(' Peak',i3)") i_Peak
         Call EI_PolGridDel(Nr_IntFerCh(i_chunk), FitDelay, IntfNuDim, i_chunk, SourcePos(1,i_peak), AntPeak_OffSt(1,i_Peak), &
               W_ap(1,i_peak), W_at(1,i_peak), Cnu0(0,1,i_peak), Cnu1(0,1,i_peak), Outpt, DelChi, Label)
         PeakChiSQ(i_peak) = Chi2pDF
           !
         Do j=-N_fit,N_fit
            Do j_IntFer=1,Nr_IntFerCh(i_Chunk)   ! Loop over selected antennas
               i_eqn=i_eqn+2
               R(i_eqn-1)=DelChi(j,2*j_IntFer-1)
               R(i_eqn)=DelChi(j,2*j_IntFer)
            Enddo
         Enddo
   EndDo
    !val=sum(ChiSq)
    !If(prn) Write(2,*) ' total chisq ',sum(ChiSq),(sum(ChiSq(:,i_Peak)),i_Peak=1,PeakNrTotal)
    !stop  'CompareCorrTime'
    !
    ! If(N_FitStatTim .gt. 0) stop
    !OPEN(UNIT=4,STATUS='unknown',FILE=trim(FileFitResult)//'.dat' ) !'plot/FitResult.dat')
    !write(4,"('!CoreDist[m], phi[rad]',3x,'I',11x,'I_calc',8x,'sigma_I',7x,&
    !    'Q',9x,'Q_calc',4x,'sigma_Q',3x,'U',9x,'U_calc',4x,'sigma_U',3x,'V',9x,'V_calc',4x,'sigma_V')")
    !close(unit=4)
    return
End Subroutine CompareEI
!==================================================
!-----------------------------------------------
Subroutine EI_PolGridDel(Nr_IntFer, FitDelay, i_sample, i_chunk , VoxLoc, AntPeak_OffSt, &
            W_ap, W_at, Cnu0, Cnu1, Outpt, DelChi, Label)
   ! Needs: Ant_RawSourceDist as used when reading in data (when fitted time-shifts = 0)
   !alculates: 1) trace in frequency space for even and odd antennas symmetrically around
   !  calculated pulse time (based on suggested source time and position) for each antenna pair
   !           2) noise background level for each antenna, needed for chi^2 calculation
   !Purpose: to be used in interferometry calculation for a voxel in the vicinity of the
   !  guessed source and with fitted antenna timings
   ! i_s=i_sample=IntfLead + 1+i_slice*N_smth
   ! i_sample is the central sample for wrapping the smooth window, counten for the time-subtrace of length=2*IntfNuDim.
   !--------------------------------------------
   use constants, only : dp, pi, ci
   use DataConstants, only : DataFolder, OutFileLabel
   use Chunk_AntInfo, only : Ant_Stations, Ant_IDs, Ant_nr, Ant_pos
   use Interferom_Pars, only :  Nr_IntFerMx, Nr_IntferCh ! the latter gives # per chunk
   Use Interferom_Pars, only : IntFer_ant, N_fit, N_smth, smooth
   Use Interferom_Pars, only : StI, StI12, StQ, StU, StV, StI3, StU1, StV1, StU2, StV2, P_un, P_lin, P_circ
   Use Interferom_Pars, only : dStI, dStI12, dStQ, dStU, dStV, dStI3, dStU1, dStV1, dStU2, dStV2, Chi2pDF
   use AntFunCconst, only : Freq_min, Freq_max,Ji_p0,Ji_t0,Ji_p1,Ji_t1, Gain  !J_0p,J_0t,J_1p,J_1t,
   use GLEplots, only : GLEplotControl
   use Interferom_Pars, only :IntfNuDim, dnu, inu1, inu2
   use FFT, only : RFTransform_CF, RFTransform_CF2CT
   use StationMnemonics, only : Statn_ID2Mnem, Station_ID2Mnem
   Implicit none
   Integer, intent(in) :: Nr_IntFer, i_sample,i_chunk
   Real(dp), intent(in) :: VoxLoc(1:3), W_ap(1:Nr_IntFer), W_at(1:Nr_IntFer)
   Real(dp), intent(in) :: FitDelay(*), AntPeak_OffSt(*)
   Complex(dp), intent(in) :: Cnu0(0:IntfNuDim,Nr_IntFer), Cnu1(0:IntfNuDim,Nr_IntFer)
   Integer, intent(in) :: Outpt ! <=0: minimal ; =1: 1line ; =2: chi^2 info to files
   Character(len=8), intent(in) :: Label
   Real(dp), intent(out) :: DelChi(-N_fit:N_fit,*)
   Real(dp) :: Power_p(1:Nr_IntFer), Power_t(1:Nr_IntFer)
   complex(dp), parameter :: ipi=ci*pi
   integer :: i_ant, j_IntFer, i, j, i_freq, i_nu, m, n, i_s, n_s
   Real(dp) :: Vec_p(1:3), Vec_t(1:3), Ras(1:3),  p_PB(1:3), t_PB(1:3)
   Real(dp) :: NEh2PB(1:3,1:3) ! NEh to Polarization Basis; NEh2VHR(NEh,VHR); formerly known as VHRBasis(1:3,1:3)
   Real(dp) :: HorDist, Thet_r ,Phi_r, dfreq, D, W, nu
   Real(dp) :: thet_d, Phi_d  ! AntSourceD(1:Nr_IntFer),
   Complex(dp) :: Sp, St
   Complex(dp) :: FTime_PB(1:2*IntfNuDim,1:3), Fnu_PB(0:IntfNuDim,1:3)
   Complex(dp) :: Phase_0, dPhase_0, Phase_1, dPhase_1, nu_p(0:IntfNuDim), nu_t(0:IntfNuDim)
   Real(dp) :: RDist, dt_AntPix_0, dt_AntPix_1, sw_p, sw_t, W_p, W_t
   Real(dp) :: AMat(1:3,1:3), Ai(1:3,1:3)
   Complex :: Stk(3,3), AiF(1:3)
   Real(dp) :: Esq_ak, FdotI
   Real(dp) :: St8
   logical :: TestCh2
   Real(dp) :: D_a(1:Nr_IntFer), Ph_a(1:Nr_IntFer)   ! just for printing & plotting
   !Real(dp) :: TestW_p(1:Nr_IntFer), TestW_t(1:Nr_IntFer)   ! just for printing & plotting
   Real(dp) :: wap_PB(1:3,1:Nr_IntFer), wat_PB(1:3,1:Nr_IntFer)  ! weighted p and t directions in polarization basis (_PB)
   Complex(dp) :: wEnu_p(0:IntfNuDim), wEnu_t(0:IntfNuDim)
   Complex(dp) :: wETime_p(1:2*IntfNuDim), wETime_t(1:2*IntfNuDim)
   Complex(dp) :: wEtime_ap(-N_smth:N_smth,1:Nr_IntFer), wEtime_at(-N_smth:N_smth,1:Nr_IntFer) ! weighted E fields, time dependent
   Real(dp) :: dChi_ap(1:Nr_IntFer), dChi_at(1:Nr_IntFer)
   Real(dp) :: AveAmp_ap(1:Nr_IntFer), AveAmp_at(1:Nr_IntFer)
   !Complex(dp) :: SumDiff, SumSq, Del_p, Del_t
   Real(dp) :: SumDiff, SumSq, Del_p, Del_t
   character(len=8) :: txt
   character(len=20) :: FMT
   CHARACTER(LEN=60) :: GLE_file
   !
   TestCh2=.false.
   If(Outpt.ge.2) TestCh2=.true.
   ! Basis for Stokes parameters
!   D=SUM(PixLoc(:)*PixLoc(:))  ; HVRBasis(:,3)=PixLoc(:)/sqrt(D)  ! radial, out
!   HVRBasis(1,1)=-PixLoc(2)  ; HVRBasis(2,1)=PixLoc(1) ; HVRBasis(3,1)=0.  ! horizontal
!   HorDist=HVRBasis(1,1)*HVRBasis(1,1)+ HVRBasis(2,1)*HVRBasis(2,1)  ;  HVRBasis(:,1)=HVRBasis(:,1)/sqrt(HorDist)
!   HVRBasis(1,2)=HVRBasis(2,3)*HVRBasis(3,1)- HVRBasis(3,3)*HVRBasis(2,1) ! vertical
!   HVRBasis(2,2)=HVRBasis(3,3)*HVRBasis(1,1)- HVRBasis(1,3)*HVRBasis(3,1)
!   HVRBasis(3,2)=HVRBasis(1,3)*HVRBasis(2,1)- HVRBasis(2,3)*HVRBasis(1,1)
!   NEh2VHR: VHRBasis(NEh,VHR)
   D=SUM(VoxLoc(:)*VoxLoc(:))  ; NEh2PB(:,3)=VoxLoc(:)/sqrt(D)  ! radial, out (North,E=0,h)=> VHR(3)~(N,0,h)
   NEh2PB(1,2)=-VoxLoc(2)  ; NEh2PB(2,2)=VoxLoc(1) ; NEh2PB(3,2)=0.  ! horizontal  VHR(2)~(0,N,0) to the right
   HorDist=NEh2PB(1,2)*NEh2PB(1,2)+ NEh2PB(2,2)*NEh2PB(2,2)  ;  NEh2PB(:,2)=NEh2PB(:,2)/sqrt(HorDist)
   NEh2PB(1,1)=NEh2PB(2,3)*NEh2PB(3,2)- NEh2PB(3,3)*NEh2PB(2,2) ! vertical tilted toward
   NEh2PB(2,1)=NEh2PB(3,3)*NEh2PB(1,2)- NEh2PB(1,3)*NEh2PB(3,2)  ! VHR(1)~(-h,0,N)
   NEh2PB(3,1)=NEh2PB(1,3)*NEh2PB(2,2)- NEh2PB(2,3)*NEh2PB(1,2)
   !
   !write(2,*) VoxLoc(:), NEh2PB
   !write(2,*) 'w_ap',w_ap
   !write(2,*) 'w_at',w_at   !
   wEnu_p(:)=0.
   wEnu_t(:)=0.
   !i_s=IntfLead  +  1+i_slice*N_smth  ! approximately correct, upto rounding errors for dt
   i_s=i_sample
   !write(2,*) 'EI_PolGridDel:',i_s,VoxLoc(:)
   !flush(unit=2)
   Fnu_PB(:,:)=0.
   AMat(:,:)=0.
   SumSq=0.
   SumDiff=0.
   Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
      !  Get interferometry phase shifts
      i_ant=IntFer_ant(j_IntFer,i_chunk)
      Call RelDist(VoxLoc(1),Ant_pos(1,i_ant,i_chunk),RDist)
      dt_AntPix_0 =Rdist - AntPeak_OffSt(i_ant) + FitDelay(i_ant)
      dt_AntPix_1 =Rdist - AntPeak_OffSt(i_ant+1) + FitDelay(i_ant+1)
      dphase_0 = exp(-ipi*dt_AntPix_0 /IntfNuDim)
      dphase_1 = exp(-ipi*dt_AntPix_1 /IntfNuDim)
      Phase_0 =exp(-ipi*dt_AntPix_0 *inu1/IntfNuDim)
      Phase_1 =exp(-ipi*dt_AntPix_1 *inu1/IntfNuDim)
      !
      !If(j_intFer.eq.110) write(2,*) 'dt_AntPix_0:',i_ant, j_IntFer, dt_AntPix_0, Rdist, AntPeak_OffSt(i_ant), dt_AntPix_1
      !If(j_intFer.eq.110) write(2,*) 'dt_AntPix_1:',i_ant+1, j_IntFer, dt_AntPix_1, Rdist, AntPeak_OffSt(i_ant+1)
      Ras(1)=(VoxLoc(1)-Ant_pos(1,i_ant,i_chunk))/1000.  ! \vec{R}_{antenna to source}
      Ras(2)=(VoxLoc(2)-Ant_pos(2,i_ant,i_chunk))/1000.
      Ras(3)=(VoxLoc(3)-Ant_pos(3,i_ant,i_chunk))/1000.
      HorDist= Ras(1)*Ras(1) + Ras(2)*Ras(2)  ! =HYPOT(X,Y)
      D=sqrt(HorDist + Ras(3)*Ras(3))
      !AntSourceD =D
      HorDist=sqrt( HorDist ) ! =HYPOT(X,Y)
      Thet_r=atan(HorDist,Ras(3))  ! Zenith angle
      Phi_r=atan2( Ras(2) , Ras(1) )  ! \phi=0 = north
      Thet_d =Thet_r*180/pi
      Phi_d =Phi_r*180/pi
      D_a(j_IntFer)=D      ! Just for plotting
      Ph_a(j_IntFer)=Phi_d
      !
      Call AntFun_Inv(thet_d ,Phi_d ) ! sets ,Ji_p0,Ji_t0,Ji_p1,Ji_t1; Inverse Jones; gives _p & _t polarized fields when contracted with (0.1) voltages
      Vec_p(1)=sin(Phi_r)              ; Vec_p(2)=-cos(Phi_r)        ; Vec_p(3)=0.
      Vec_t(1)=-cos(Thet_r)*Vec_p(2)  ; Vec_t(2)=cos(Thet_r)*Vec_p(1) ; Vec_t(3)=-sin(Thet_r)
      !
      Do i=1,3  ! convert t&p orientations from (NEh) to polar base
         p_PB(i)=SUM( NEh2PB(:,i)*Vec_p(:) )
         t_PB(i)=SUM( NEh2PB(:,i)*Vec_t(:) )
      Enddo
      !write(2,"(4(A,3F7.3),3F7.3)") 'p_PB=',p_PB(:),', t_PB=',t_PB(:)! &
      !write(2,*) Vec_p(:), Ras(:)
      !flush(unit=2)
      !
      w_p=W_ap(j_IntFer)
      w_t=W_at(j_IntFer)  ! do not bother about scale factor, divides out anyhow.
      sw_p=sqrt(W_ap(j_IntFer))
      sw_t=sqrt(W_at(j_IntFer))
      wap_PB(:,j_IntFer)=p_PB(:)*sw_p/D
      wat_PB(:,j_IntFer)=t_PB(:)*sw_t/D
      !
      Do i_nu=inu1,inu2   ! Increment frequency spectrum with the spectrum from this antenna
         nu=i_nu*dnu
         i_freq=Int(nu)
         dfreq=nu-i_freq
         Sp=((1.-dfreq)*Ji_p0(i_freq) + dfreq*Ji_p0(i_freq+1)) *Cnu0(i_nu,j_IntFer)*phase_0 + &
            ((1.-dfreq)*Ji_p1(i_freq) + dfreq*Ji_p1(i_freq+1)) *Cnu1(i_nu,j_IntFer)*phase_1
         St=((1.-dfreq)*Ji_t0(i_freq) + dfreq*Ji_t0(i_freq+1)) *Cnu0(i_nu,j_IntFer)*phase_0 + &
            ((1.-dfreq)*Ji_t1(i_freq) + dfreq*Ji_t1(i_freq+1)) *Cnu1(i_nu,j_IntFer)*phase_1
         Fnu_PB(i_nu,:)=Fnu_PB(i_nu,:) + ( p_PB(:)* Sp*w_p + t_PB(:)* St*w_t ) *Gain(i_freq)/D
         wEnu_p(i_nu)=  Sp*sw_p *Gain(i_freq)
         wEnu_t(i_nu)=  St*sw_t *Gain(i_freq)
         phase_0 =phase_0 *dphase_0
         phase_1 =phase_1 *dphase_1
      Enddo
      !
      !write(2,*) 'j_IntFer:',j_IntFer, i_ant, wEnu_p(inu1+1), abs(SP), phase_0, dphase_0
      !flush(unit=2)
      ! Calculate contribution to AMat       == A matrix
      Do i=1,3
         Do j=1,3  ! a single factor distance since the other is already included in W (as in W_p & W_t)
            AMat(i,j)=AMat(i,j)+ ( p_PB(i)*p_PB(j)*w_p + t_PB(i)*t_PB(j)*w_t)/(D*D)
         Enddo
      Enddo
      !
      Call RFTransform_CF2CT(wEnu_p(0),wETime_p(1) ) !time trace measured E-field for antenna pair
      Call RFTransform_CF2CT(wEnu_t(0),wETime_t(1) )
      Do j=-N_smth,N_smth ! copy for later use
         wEtime_ap(j,j_IntFer)=wETime_P(i_s+j)  !time trace measured E-field for antenna pair
         wEtime_at(j,j_IntFer)=wETime_t(i_s+j)
      Enddo
      !
   Enddo     ! j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
   !write(2,*) 'Sp',i_ant, i_s, abs(Sp), abs(St), phase_0, phase_1
   !flush(unit=2)
   !
   Call matinv3(AMat, Ai)
   Do m=1,3  ! convert to time
      Call RFTransform_CF2CT(Fnu_PB(0,m),FTime_PB(1,m) )
   Enddo
   !
   !write(2,*) 'EI_PolGridDel;Ai',Ai
   !flush(unit=2)
   !write(2,*) 'FTime_PB(1:10,1):',FTime_PB(1:10,1)
   !
   Stk(:,:)=0.
   FdotI=0.
   SumSq=0.
   If(TestCh2) then
      !Write(2,*) 'wap:',W_ap(1:10)
      !Write(2,*) 'gain:',Gain(Freq_min:Freq_min+10)
      !Write(2,*) 'FTime_PB:',FTime_PB(100:110,1)
      !write(2,*) 'EI_PolGridDel;FTime_PB(1:10,1):',FTime_PB(1:10,1)
      !flush(unit=2)
      Do i=1,8
         If( Label(i:i).eq.' ') then
            txt(i:i)='-'
         Else
            txt(i:i)=Label(i:i)
         EndIf
      EndDo
      !i=188 ! 111   !            <<<<<<<<<<<<<<<<<<<<<<<<
      !OPEN(UNIT=30,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'DelChi2_'//TRIM(txt)//'.dat')
      !write(30,"(3I4,I7,' 0')") i_sample, Nr_IntFer, i, i_s
      !Write(30,"(A,9x,A,13x,';   ', 9(A,'/I [%]   ;   '),A )") &
       !  '! i_ant, j; Ep-data;  model , Et-data, model, d_chi_p','d_chi_t'
      !
      write(FMT,"(A,I4,A)") '(I4,',Nr_IntFer,'G13.3)'
      OPEN(UNIT=31,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'PhiDat_'//TRIM(txt)//'.dat') ! for curtain plots
      OPEN(UNIT=32,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'ThDat_'//TRIM(txt)//'.dat')
      OPEN(UNIT=33,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'PhiMod_'//TRIM(txt)//'.dat')
      OPEN(UNIT=34,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'ThMod_'//TRIM(txt)//'.dat')
      Power_p(:)=0.
      Power_t(:)=0.
      AveAmp_ap(:)=0.
      AveAmp_at(:)=0.
   EndIf
   If(N_fit.gt.0) then
      N_s=N_fit
   Else
      N_s=N_smth
   EndIf
   dChi_ap(:)=0.  ;     dChi_at(:)=0.
   Do j=-N_s,N_s ! fold time-dependence with smoothing function
      Do m=1,3  ! convert to time
         AiF(m)=SUM( Ai(m,:)*FTime_PB(i_s+j,:) )
      Enddo
      !write(2,*) 'FTime_PB(i_s+j,:)',i_s,j,FTime_PB(i_s+j,:)
      !flush(unit=2)
      FdotI=FdotI + smooth(j)*Real( SUM( FTime_PB(i_s+j,:)*Conjg(AiF) ) )  ! Imag is zero
      Do m=1,3
         Do n=1,3
            Stk(m,n)= Stk(m,n) + smooth(j)*AiF(m)*Conjg(AiF(n))
         Enddo
      Enddo
      W=smooth(j)*N_smth  ! to compensate the averaging in the smoothing fie and to get the value per degr-of-freedm
      Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
      !   Del_p = ABS( wEtime_ap(j,j_IntFer) -SUM( wap_PB(:,j_IntFer)*AiF(:) ) ) ! true difference (measured-calculated)for this antenna
      !   Del_t = ABS( wEtime_at(j,j_IntFer) -SUM( wat_PB(:,j_IntFer)*AiF(:) ) )
      ! NB should absolute difference be fitted or the real part?
         Del_p = Real( wEtime_ap(j,j_IntFer) -SUM( wap_PB(:,j_IntFer)*AiF(:) ) ) ! true difference (measured-calculated)for this antenna
         Del_t = Real( wEtime_at(j,j_IntFer) -SUM( wat_PB(:,j_IntFer)*AiF(:) ) )
         dChi_ap(j_IntFer)=dChi_ap(j_IntFer) + W* Del_p**2 ! wEtime_ap(j,j_IntFer) * Conjg( Del_p  ) ! Imag is zero ??
         dChi_at(j_IntFer)=dChi_at(j_IntFer) + W* Del_t**2 ! wEtime_at(j,j_IntFer) * Conjg( Del_t  ) ! Imag is zero ??
         ! SumDiff=SumDiff+ sqrt(W)*(Del_p + Del_t)
         If(N_fit.le.0) cycle
         DelChi(j,2*j_IntFer-1)=sqrt(W)*Del_p
         DelChi(j,2*j_IntFer)  =sqrt(W)*Del_t
         If(TestCh2) then
            Power_p(j_IntFer)=Power_p(j_IntFer) + W*ABS(wEtime_ap(j,j_IntFer))**2
            Power_t(j_IntFer)=Power_t(j_IntFer) + W*ABS(wEtime_at(j,j_IntFer))**2
            AveAmp_ap(j_IntFer)=AveAmp_ap(j_IntFer) + sqrt(W)*ABS(wEtime_ap(j,j_IntFer))
            AveAmp_at(j_IntFer)=AveAmp_at(j_IntFer) + sqrt(W)*ABS(wEtime_at(j,j_IntFer))
         EndIf
      Enddo
      !write(2,*) 'j',j,smooth(j), Del_p, Del_t,AiF(1), FdotI
      !flush(unit=2)
      If(TestCh2) then
         !j_IntFer=i
         !write(2,*) 'EI_PolGridDel;Del_p:', j, Del_p, wEtime_ap(j,j_IntFer), SUM( wap_PB(:,j_IntFer)*AiF(:) )
         !flush(unit=2)
         !write(30,"(I4,I5,10G13.3)")  j_IntFer, j, &
         !   wEtime_ap(j,j_IntFer), SUM( wap_PB(:,j_IntFer)*AiF(:) ), &
         !   wEtime_at(j,j_IntFer), SUM( wat_PB(:,j_IntFer)*AiF(:) ), &
          !  dChi_ap(j_IntFer), dChi_at(j_IntFer)
         write(31,FMT)  j, REAL( wEtime_ap(j,1:Nr_IntFer) )  ! Data
         write(32,FMT)  j, REAL( wEtime_at(j,1:Nr_IntFer) )
         write(33,FMT)  j,( REAL( SUM( wap_PB(:,j_IntFer)*AiF(:) ) ), j_IntFer=1,Nr_IntFer)  ! model
         write(34,FMT)  j,( REAL( SUM( wat_PB(:,j_IntFer)*AiF(:) ) ), j_IntFer=1,Nr_IntFer)
      EndIf
   Enddo
   dChi_ap(1:Nr_IntFer)=dChi_ap(1:Nr_IntFer)/(2*N_s+1.)
   dChi_at(1:Nr_IntFer)=dChi_at(1:Nr_IntFer)/(2*N_s+1.)
   SumSq=SUM(dChi_ap(1:Nr_IntFer)) +SUM(dChi_at(1:Nr_IntFer)) ! smooth(j)* (Del_p**2 + Del_t**2) ! True contribution to chi^2 from this antenna
   Chi2pDF=SumSq/(2.*Nr_IntFer)
   If(TestCh2) then
      Power_p(:)=Power_p(:)/W_ap(:)
      Power_t(:)=Power_t(:)/W_at(:)
      AveAmp_ap(:)=AveAmp_ap(:)/sqrt(W_ap(:))
      AveAmp_at(:)=AveAmp_at(:)/sqrt(W_at(:))
      !Close(UNIT=30)
      Close(UNIT=31)
      Close(UNIT=32)
      Close(UNIT=33)
      Close(UNIT=34)
      Write(GLE_file,"('CuP-',A,A)") TRIM(OutFileLabel),trim(txt)
      Call GLEscript_Curtains(30, GLE_file, N_s, i_chunk, trim(DataFolder)//TRIM(OutFileLabel), TRIM(txt), &
               dChi_ap, dChi_at, Power_p, Power_t, Chi2pDF)
      Write(2,"(16x,A,5x,A,5x,A,6x,A,2x,A,4x,A)") 'D', 'Power','weight','delta_chi^2'!,'Power','d_chi/sqrt(P)'
      !Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
      !   i_ant=IntFer_ant(j_IntFer,i_chunk)
      !   Write(2,"(I3,1x,A5,I4, F6.1,1x,A,F10.1,F10.3,F9.2, F10.2,F9.2)") &
      !                  j_IntFer, Statn_ID2Mnem(Ant_Stations(i_ant,i_chunk)),Ant_IDs(i_ant,i_chunk),D_a(j_IntFer), &
      !                  'p:', Power_p(j_IntFer),W_ap(j_IntFer),dChi_ap(j_IntFer), &
      !                  AveAmp_ap(j_IntFer)**2/Power_p(j_IntFer)
      !   Write(2,"(20x,A,F10.1,F10.3,F9.2, F10.2,F9.2)") 't:', Power_t(j_IntFer),W_at(j_IntFer),dChi_at(j_IntFer), &
      !                  AveAmp_at(j_IntFer)**2/Power_t(j_IntFer)
      !Enddo
   EndIf
   !
   !Chi2pDF=(Esq_ak-FdotI)/(2.*Nr_IntFer)    !  Reasonable approximation to chi^2
   StI = Real(Stk(1,1) + Stk(2,2) + Stk(3,3))  ! =\Lambda_0  from  PHYSICAL REVIEW E 66, 016615 ~2002!
   StI12=Real(Stk(1,1) + Stk(2,2))             ! =
   StI3 =  Real(Stk(3,3))                      ! =>\Lambda_8 = StI12 *sq(3)/2 - StI3 *sq(3)
   StQ = Real(Stk(1,1)-Stk(2,2))               ! =\Lambda_3 *2/3
   StU = 2*Real(Stk(1,2))                      ! =\Lambda_1 *2/3
   StV =  2*Imag(Stk(1,2))                     ! =\Lambda_2 *2/3
   StU2 = 2*Real(Stk(3,1))                     ! =\Lambda_4 *2/3
   StV2 =  2*Imag(Stk(3,1))                    ! =\Lambda_5 *2/3
   StU1 = 2*Real(Stk(2,3))                     ! =\Lambda_6 *2/3
   StV1 =  2*Imag(Stk(2,3))                    ! =\Lambda_7 *2/3
   St8 = (StI12  - 2* StI3)/sqrt(3.)           ! =\Lambda_8 *2/3
   P_un = 1.- (3./4.)*(StQ**2+StU**2+StV**2+StU1**2+StV1**2+StU2**2+StV2**2+St8**2)/(StI**2)
   P_lin = (3./4.)*(StQ**2+StU**2+StU1**2+StU2**2+St8**2)/(StI**2)
   P_circ= (3./4.)*(StV**2+StV1**2+StV2**2)/(StI**2)
   !
   Do m=1,3
      Do n=1,3
         AMat(m,n)=SUM( Ai(m,:)*Ai(:,n) )*Chi2pDF**2  /4. ! this should correspond to the square of the error
         ! arbitrary factor 100 included to stretch error bars and make them visible
      EndDo
    !  write(2,*) 'Correlation matrix(m,:)=',AMat(m,:)
   EndDo
   dStI = SQRT(SUM(AMat(:,:) ))
   dStI12=SQRT(AMat(1,1) + AMat(2,2))
   dStI3 =  SQRT( AMat(3,3) )
   dStQ = dStI12
   dStU = 2*SQRT( ABS(AMat(1,2)) )
   dStV =  dStU
   dStU2 = 2*SQRT( ABS(AMat(3,1)) )
   dStV2 =  dStU2
   dStU1 = 2*SQRT( ABS(AMat(2,3)) )
   dStV1 =  dStU1
   !
   If(Outpt.ge.1) Then
      W=ATAN2(StU,StQ)/2.*180./pi
      If(W.lt.-90.) W=W+180.
      write(2,"(A,2(A,G12.4),10(A,F7.2))") Label,': I123=',StI,', I12=',StI12, ', Q/I=',StQ/StI, &
            ', U/I=',StU/StI, ', V/I=',StV/StI, ', I3/I=',StI3/StI,', angle=',W,', chi^2=',Chi2pDF &
            ,', P_unpol=',P_un, ', P_lin=',P_lin, ', P_circ=',P_circ
      If(TestCh2 .and. Outpt.lt.1) then  !  i.e. always skip this part
         SumSq=SumSq/(2.*Nr_IntFer)
         !txt(1:2)=Label(7:8) !  write(txt,"(I2.2)") i_sample
         OPEN(UNIT=30,STATUS='unknown',ACTION='WRITE', &
               FILE=trim(DataFolder)//TRIM(OutFileLabel)//'IntfSpecChi2_'//TRIM(txt)//'.dat')
         !write(2,*) 'Chi^2 plot info written to file: "', &
         !   trim(DataFolder)//TRIM(OutFileLabel)//'IntfSpecChi2_'//TRIM(txt)//'.dat','"'
         !Flush(unit=2)
         write(30,"(2I4,F9.2,2G12.3,I7,' 0')") i_sample, Nr_IntFer, SumSq, StI12, StI, i_s
         Write(30,"(A,9x,A,13x,';   ', 9(A,'/I [%]   ;   '),A )") &
            '! nr, i_ant, Ant#, St#;  D[km] , phi[deg], d_chi_ph, d_chi_th','StI','StI12'
         w=(W_ap(1)+W_at(1))/2.
         Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
            i_ant=IntFer_ant(j_IntFer,i_chunk)
            write(30,"(I4,I5,I6,I6,F9.3,F8.1,5G12.3)")  j_IntFer, i_ant, &
               Ant_IDs(I_ant,i_chunk), Ant_Stations(I_ant,i_chunk), D_a(j_IntFer), Ph_a(j_IntFer), &
               dChi_ap(j_IntFer), dChi_at(j_IntFer), &
               W_ap(j_IntFer)/W, W_at(j_IntFer)/W
         Enddo
         Close(Unit=30)
         !Call GLEplotControl(PlotType='EI_delta-chisq', PlotName='EI_delchisq_'//TRIM(txt)//TRIM(OutFileLabel), &
         !      PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel)//' "'//TRIM(txt)//'"' ) ! does NOT work in windows
         !Call GLEplotControl(SpecialCmnd='rm '//TRIM(DataFolder)//TRIM(OutFileLabel)//'DelChi2_'//TRIM(txt)//'.dat') !  Command to delete files
         !Plot of chi^2 value per antenna-pair:
         Call GLEplotControl(PlotType='EI_chisq', PlotName='EI_chisq_'//TRIM(txt)//TRIM(OutFileLabel), &
               PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel)//' "'//TRIM(txt)//'"' ) ! does NOT work in windows
         Call GLEplotControl(SpecialCmnd='rm '//TRIM(DataFolder)//TRIM(OutFileLabel)//'IntfSpecChi2_'//TRIM(txt)//'.dat') !  Command to delete files
      EndIf
   EndIf
   !
   Return
End Subroutine EI_PolGridDel
!-----------------------------------------------

Subroutine EI_PolSetUp(Nr_IntFer, IntfBase, i_chunk, VoxLoc, AntPeak_OffSt, Cnu0, Cnu1, W_ap, W_at)
   !--------------------------------------------
   use constants, only : dp, pi, ci
   use DataConstants, only : Time_dim, DataFolder, OutFileLabel, Ant_nrMax
   use Chunk_AntInfo, only : CTime_spectr, Ant_Stations, Ant_IDs, Ant_nr, Ant_pos, Ant_RawSourceDist
   use Chunk_AntInfo, only : NormOdd, NormEven !Powr_eo !, NAnt_eo
   Use Interferom_Pars, only :  SumWindw, N_smth, smooth  !, NrPixSmPowTr
   Use Interferom_Pars, only : IntFer_ant, Nr_IntFerMx, Nr_IntFerCh
   Use Interferom_Pars, only : CTime_p, CTime_t, Noise_p, Noise_t  !PixLoc, CenLoc,
   use AntFunCconst, only : Freq_min, Freq_max,Ji_p0,Ji_t0,Ji_p1,Ji_t1, Gain  !J_0p,J_0t,J_1p,J_1t,
   use Interferom_Pars, only :IntfNuDim, dnu, inu1, inu2
   use GLEplots, only : GLEplotControl
   use FFT, only : RFTransform_CF, RFTransform_CF2CT
   Implicit none
   Integer, intent(in) :: Nr_IntFer, IntfBase, i_chunk  ! IntfLead
   Real(dp), intent(in) :: VoxLoc(1:3)
   Real(dp), intent(out) :: AntPeak_OffSt(1:Ant_nrMax)
   Complex(dp), intent(out) :: Cnu0(0:IntfNuDim,Nr_IntFerMx), Cnu1(0:IntfNuDim,Nr_IntFerMx)
   Real(dp), intent(out) :: W_ap(1:Nr_IntFer), W_at(1:Nr_IntFer)
   complex(dp), parameter :: ipi=ci*pi
   integer :: i_ant, j_IntFer, i, j, i_freq, i_nu, IntfDim, i_s, m, n, SamplOff_0, SamplOff_1
   Real(dp) :: Ras(1:3), HorDist, Thet_r ,Phi_r, dfreq, D, W,  nu
   Real(dp) :: thet_d, Phi_d  ! AntSourceD(1:Nr_IntFer),
   Real(dp) :: Rtime(1:2*IntfNuDim)
   Complex(dp) :: Sp, St
   Complex(dp) :: Phase, dPhase, nu_p(0:IntfNuDim), nu_t(0:IntfNuDim)
   Real(dp) :: RDist, dt_AntPix
   Real(dp) :: W_p, W_t  !, SumDiff, SumSq, Del_p, Del_t, sw_p, sw_t
   !
   IntfDim=2*IntfNuDim  !  should be about =512
   !
   !PowerScale=1.
   !NormEven=sqrt(PowerScale*2.*Powr_eo(0)/(Powr_eo(0)+Powr_eo(1)) )/100.
   !NormOdd=sqrt(PowerScale*2.*Powr_eo(1)/(Powr_eo(0)+Powr_eo(1)) )/100.
   Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas to detremine t- and p- polarized trace for central pixel; needed for noise estimate
      i_ant=IntFer_ant(j_IntFer,i_chunk)
      !write(2,*) i_ant, j_IntFer, Ant_pos(:,i_ant,i_chunk)
      Ras(1)=(VoxLoc(1)-Ant_pos(1,i_ant,i_chunk))/1000.
      Ras(2)=(VoxLoc(2)-Ant_pos(2,i_ant,i_chunk))/1000.
      Ras(3)=(VoxLoc(3)-Ant_pos(3,i_ant,i_chunk))/1000.
      HorDist=sqrt(  Ras(1)*Ras(1) + Ras(2)*Ras(2)  )
      Thet_r=atan(HorDist,Ras(3))  ! Zenith angle
      Phi_r=atan2( Ras(2) , Ras(1) )  ! \phi=0 = north
      Thet_d =Thet_r*180/pi
      Phi_d =Phi_r*180/pi
      !If(j_IntFer.eq.1) Then
      !   I_Scale= (HorDist*HorDist + Ras(3)*Ras(3))/Nr_IntFer ! scales value of chi-square, not intensity, and value of F
      !   !write(2,*) 'AmplitudeW NormEven,odd',NormEven,NormOdd,', ratio=',NormEven/NormOdd, Noise
      !EndIf
      !Calculation of noise-power level
      Call AntFun_Inv(thet_d ,Phi_d ) ! sets ,Ji_p0,Ji_t0,Ji_p1,Ji_t1; Inverse Jones; gives _p & _t polarized fields when contracted with (0.1) voltages
      w_p=SUM( ((NormEven*ABS(Ji_p0(Freq_min:Freq_max)))**2+(NormOdd*ABS(Ji_p1(Freq_min:Freq_max)))**2)* &
         Gain(Freq_min:Freq_max)**4)
      w_t=SUM( ((NormEven*ABS(Ji_t0(Freq_min:Freq_max)))**2+(NormOdd*ABS(Ji_t1(Freq_min:Freq_max)))**2)* &
         Gain(Freq_min:Freq_max)**4) ! This to account for the noise level in the theta direction when pol is unfolded in the proper direction
      !   noise=(100.)^2*W; Factor 100 to copy the factor that was introduces in Antenna_Read
      Noise_p(j_IntFer)=W_p*5.d6 ! Additional /10 (approx) because effective bandwidth / full, thus degrees of freedom are over estimated
      W_ap(j_IntFer)=2.d-7/W_p!=1./W_p*1.0d3 *200.*25 !Noise_p(j_IntFer) ! =1./(Noise_p(j_IntFer)*sqrt(Pow_p/Noise_p(j_IntFer)+1.))
      Noise_t(j_IntFer)=W_t*5.d6 ! factor 200 obtained by comparing with noise trace         !
      W_at(j_IntFer)=2.d-7/W_t
      !
      Call RelDist(VoxLoc(1),Ant_pos(1,i_ant,i_chunk),RDist) ! same for even and odd antenna
      !dt_AntPix_0 =Rdist - Ant_RawSourceDist(i_ant,i_chunk)
      SamplOff_0=Int(Rdist - Ant_RawSourceDist(i_ant,i_chunk))
      !If(j_IntFer.eq.1) write(2,*) 'SamplOff_0:',SamplOff_0,Rdist, Ant_RawSourceDist(i_ant,i_chunk)
      !write(2,*) 'SamplOff_0:',SamplOff_0,j_IntFer,IntfBase,Ant_RawSourceDist(i_ant,i_chunk), i_ant,i_chunk
      !Flush(unit=2)
      AntPeak_OffSt(i_ant)= SamplOff_0 + Ant_RawSourceDist(i_ant,i_chunk)
      !dt_AntPix_0 =Rdist - Ant_OffSt(i_ant)
      SamplOff_0=IntfBase + SamplOff_0
      SamplOff_1=Int(Rdist - Ant_RawSourceDist(i_ant+1,i_chunk))
      AntPeak_OffSt(i_ant+1)= SamplOff_1 + Ant_RawSourceDist(i_ant+1,i_chunk)
      SamplOff_1=IntfBase + SamplOff_1
      !If(j_intFer.eq.110) write(2,*) 'setup:',i_ant, j_IntFer, Rdist, Ant_RawSourceDist(i_ant,i_chunk), &
      !      AntPeak_OffSt(i_ant), SamplOff_0-IntfBase
      !If(j_intFer.eq.110) write(2,*) 'setup:',i_ant+1, j_IntFer, Rdist, Ant_RawSourceDist(i_ant+1,i_chunk), &
      !      AntPeak_OffSt(i_ant+1), SamplOff_1-IntfBase
      !
      If(SamplOff_0.lt.0 .or. SamplOff_1.lt.0) Then
         Write(2,*) 'Position in reference antenna too close to beginning of time-trace for source @',VoxLoc,' in chunk#',i_chunk
         Stop 'EI_PolSetUp: too small leading buffer'
      EndIf
      If(SamplOff_0.gt.(Time_dim-IntfDim) .or. SamplOff_1.gt.(Time_dim-IntfDim)) Then
         Write(2,*) 'Position in reference antenna too close to end of time-trace for source @',VoxLoc,' in chunk#',i_chunk
         Stop 'EI_PolSetUp: too small trailing buffer'
      EndIf
      RTime(:)=REAL(CTime_spectr(SamplOff_0+1:SamplOff_0+IntfDim,i_ant,i_chunk))*NormEven  ! IntfBase=SumStrt - IntfLead
      Call RFTransform_CF(RTime,Cnu0(0,j_IntFer))
      RTime(:)=REAL(CTime_spectr(SamplOff_1+1:SamplOff_1+IntfDim,i_ant+1,i_chunk))*NormOdd
      Call RFTransform_CF(RTime,Cnu1(0,j_IntFer))
      !
      nu_p(:)=0.
      nu_t(:)=0.
      Do i_nu=inu1,inu2   ! Increment frequency spectrum with this antenna
         nu=i_nu*dnu
         i_freq=Int(nu)
         dfreq=nu-i_freq ! phase-shifts are zero for centran pixel, no timeshift!!
         nu_p(i_nu) =(  ((1.-dfreq)*Ji_p0(i_freq) + dfreq*Ji_p0(i_freq+1)) *Cnu0(i_nu,j_IntFer) + &
                        ((1.-dfreq)*Ji_p1(i_freq) + dfreq*Ji_p1(i_freq+1)) *Cnu1(i_nu,j_IntFer) ) *Gain(i_freq)
         nu_t(i_nu) =(  ((1.-dfreq)*Ji_t0(i_freq) + dfreq*Ji_t0(i_freq+1)) *Cnu0(i_nu,j_IntFer) + &
                        ((1.-dfreq)*Ji_t1(i_freq) + dfreq*Ji_t1(i_freq+1)) *Cnu1(i_nu,j_IntFer) ) *Gain(i_freq)
         ! Gain(i_freq)=sqrt( SUM(J^2) )/(Freq_max-Freq_min) ; to neutralize effect Ji on frequency, calculated in antenna_function
      Enddo
      ! convert to time
      Call RFTransform_CF2CT(nu_p(0),CTime_p(1,j_IntFer) )  ! needed for the peak/background estimate
      Call RFTransform_CF2CT(nu_t(0),CTime_t(1,j_IntFer) )
   Enddo    !  j_IntFer=1,Nr_IntFer
   !
   !Write(2,*) 'EI_PolSetUp:CTime_p(1:10_IntFer)=',CTime_p(1:10,1)
   !Write(2,*) 'EI_PolSetUp:NormOdd', IntfBase, NormOdd
   !Write(2,*) 'EI_PolSetUp:Noise_p(1:10)', Noise_p(1:10)
   !
   Return
End Subroutine EI_PolSetUp
!--------------------
!--------------------
Subroutine TimeTracePlot(j_IntFer, IntfBase, i_chunk, VoxLoc,Windw, Label)
   use constants, only : dp
   use DataConstants, only : DataFolder, OutFileLabel
   use Chunk_AntInfo, only : CTime_spectr, Ant_pos, Ant_RawSourceDist
   use Chunk_AntInfo, only : NormOdd, NormEven  ! Powr_eo !, NAnt_eo
   Use Interferom_Pars, only :  IntFer_ant, IntfNuDim, CTime_p, CTime_t
   !Use Interferom_Pars, only :  PowerScale
   use GLEplots, only : GLEplotControl
   Implicit none
   Integer, intent(in) :: j_IntFer, i_chunk, IntfBase, Windw
   Real(dp), intent(in) :: VoxLoc(1:3)
   Character(LEN=8), intent(in) :: Label
   Real(dp) :: RDist
   Integer :: i_ant, i_s, i, base
   ! Just for checking causality on the first pass
   !j_IntFer=1
   !PowerScale=1.
   !NormEven=sqrt(PowerScale*2.*Powr_eo(0)/(Powr_eo(0)+Powr_eo(1)))/100.  ! to undo the factor 100 (and more) that was introduced in antenna-read
   !NormOdd=sqrt(PowerScale*2.*Powr_eo(1)/(Powr_eo(0)+Powr_eo(1)))/100.  ! to undo the factor that was introduced in antenna-read
   i_ant=IntFer_ant(j_IntFer,i_chunk)
   Call RelDist(VoxLoc(1),Ant_pos(1,i_ant,i_chunk),RDist)
   Rdist =Rdist - Ant_RawSourceDist(i_ant,i_chunk)
   i_s=NINT(Rdist)+IntfNuDim-Windw/2
   base=IntfBase + i_s
   !write(txt,"(I2.2)") j_IntFer
   !write(2,*) 'E-field hilbert envelope for the central pixel for the antenna pair ',TRIM(Label),' written to:',&
   !      trim(DataFolder)//TRIM(OutFileLabel)//'Trace_'//TRIM(Label)//'.dat',' shift=',NINT(Rdist ),' samples'
   OPEN(UNIT=30,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'Trace_'//TRIM(Label)//'.dat')
   Write(30,*) 1,Windw
      Write(2,*) 'In TimeTracePlot, data written to file:', trim(DataFolder)//TRIM(OutFileLabel)//'Trace_'//TRIM(Label)//'.dat'
      write(2,*) 'Header: base,base+Windw, data range=',1, Windw
      Write(2,*) 1,Windw
   Do i=1, Windw
      write(30,"(i5,40G12.4)") i,ABS(CTime_spectr(Base+i,i_ant,i_chunk))*NormEven &
         ,ABS(CTime_spectr(Base+i,i_ant+1,i_chunk))*NormOdd &
         ,Abs(CTime_p(i_s+i,j_IntFer)), Abs(CTime_t(i_s+i,j_IntFer))
   EndDo
   Close(Unit=30)
   Call GLEplotControl(PlotType='PlotTimeTraces', PlotName='TimeTraces_'//TRIM(Label)//TRIM(OutFileLabel), &
         PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel)//' "'//TRIM(Label)//'"' ) ! does NOT work in windows
   !
   Return
End Subroutine TimeTracePlot
!--------------------
Subroutine  EI_Weights(Nr_IntFer, i_sample, i_chunk, PixLoc, AntPeak_OffSt, W_ap, W_at)
   use constants, only : dp
   use Chunk_AntInfo, only :  Ant_pos, Ant_RawSourceDist
   Use Interferom_Pars, only :  N_smth, smooth, IntFer_ant, IntfNuDim, CTime_p, CTime_t, Noise_p, Noise_t
   Use Interferom_Pars, only :  Nr_IntFerMx, Nr_IntFerCh
   Implicit none
   Integer, intent(in) :: Nr_IntFer, i_sample, i_chunk  ! IntfLead
   Real(dp), intent(in) :: PixLoc(1:3), AntPeak_OffSt(*)
   Real(dp), intent(out) :: W_ap(1:Nr_IntFer), W_at(1:Nr_IntFer)
   !Complex(dp), intent(in) :: CTime_p(1:2*IntfNuDim,1:Nr_IntFer), CTime_t(1:2*IntfNuDim,1:Nr_IntFer)
   !Real(dp),intent(in) :: Noise_p(1:Nr_IntFer), Noise_t(1:Nr_IntFer)
   Real(dp) :: RDist, pow_p, pow_t
   Integer :: j_IntFer, i_ant, i_s, j
   !
   !Nr_IntFer=Nr_IntFerCh(i_chunk)
   j_IntFer=1
   i_ant=IntFer_ant(j_IntFer,i_chunk)
   Call RelDist(PixLoc(1),Ant_pos(1,i_ant,i_chunk),RDist)
   Rdist =Rdist - AntPeak_OffSt(i_ant)
   i_s=i_sample + NINT(Rdist)  ! approximately correct, upto rounding errors for dt
   !Write(2,*) 'CTime_p(i_s+j,j_IntFer):',i_s,CTime_p(i_s:i_s+10,1)
  ! flush(unit=2)
   !Write(2,*) 'CTime_t(i_s+j,j_IntFer):',i_s,CTime_t(i_s,:)
   !
   Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
      !  Get interferometry phase shifts
      i_ant=IntFer_ant(j_IntFer,i_chunk)
      Call RelDist(PixLoc(1),Ant_pos(1,i_ant,i_chunk),RDist)
      !Rdist - Ant_RawSourceDist(i_ant,i_chunk) + FitDelay(i_ant)
      !Note that for calculational speed the fitdelays are assumed to be of minor importance for determining the
      ! signal-over-noise values, and thus the wieghting factors.
      Rdist =Rdist - AntPeak_OffSt(i_ant) ! Used only for weighting, fitted shifts not included
      i_s=i_sample + NINT(Rdist)  ! approximately correct, upto rounding errors for dt
   !write(2,*) 'i_s',i_s,j_intFer, i_ant, rdist, PixLoc(:)
   !flush(unit=2)
      Pow_p=0.
      Pow_t=0.
      Do j=-N_smth,N_smth ! pulse power, fold with smoothing function
         Pow_p= Pow_p + smooth(j)*Abs(CTime_p(i_s+j,j_IntFer))**2
         Pow_t= Pow_t + smooth(j)*Abs(CTime_t(i_s+j,j_IntFer))**2
      Enddo
      !The amplitudes are fit, so the sqrt(W) should be equal to the inverse of the error on the amplitude, which is the amplitude of the noise
      W_ap(j_IntFer)=1./(Noise_p(j_IntFer)+0.01*Pow_p) ! =1./(Noise_p(j_IntFer)*sqrt(Pow_p/Noise_p(j_IntFer)+1.))
      W_at(j_IntFer)=1./(Noise_t(j_IntFer)+0.01*Pow_t) ! =1./(Noise_t(j_IntFer)*sqrt(Pow_t/Noise_t(j_IntFer)+1.))
      !w_p=I_Scale*W_ap(j_IntFer)  ! to have alpha of order unity
      !w_t=I_Scale*W_at(j_IntFer)
      ! For chi-square calculation, sum_ak[ W_ak * E_ak^2] :
      !Esq_ak=Esq_ak + (Pow_p*W_ap(j_IntFer) + Pow_t*W_at(j_IntFer))  ! Used in approximation to chi^2
      !write(2,*) 'Noise_p(',Noise_p(j_IntFer),Pow_p
   Enddo
   !write(2,*) 'Noise_p(',Noise_p(1:10)
   !write(2,*) 'W_ap(',W_ap(1:10)
   !Stop
   !
   Return
End Subroutine EI_Weights
! ------------------------
Pure Subroutine GetInterfFitDelay(i_chunk, FitDelay)
    use constants, only : dp
    use DataConstants, only : Station_nrMax, Ant_nrMax
    use ThisSource, only :  PeakNrTotal, PeakPos, Peak_eo, ChunkNr
    Use Interferom_Pars, only : IntfNuDim, IntFer_ant, Nr_IntFerMx, Nr_IntferCh ! the latter gives # per chunk
    use Chunk_AntInfo, only : Ant_Stations, Ant_pos
    use Chunk_AntInfo, only : Unique_StatID, Nr_UniqueStat, Ant_IDs, Unique_SAI, Tot_UniqueAnt
    !use FitParams, only : N_FitPar, N_FitStatTim, FitParam, X_Offset
    use FitParams, only : Fit_TimeOffsetStat, Fit_TimeOffsetAnt, Fit_AntOffset
    use StationMnemonics, only : Station_ID2Mnem
  implicit none
  real(dp), intent(out) :: FitDelay(*)
  integer, intent(in) :: i_chunk
    !
    integer :: i_SAI, i_stat, i_ant, Station_ID,j_IntFer, Antenna_SAI
    !
   Do j_IntFer=1,Nr_IntferCh(i_chunk)
      i_ant=IntFer_ant(j_IntFer,i_chunk)
      Station_ID = Ant_Stations(i_ant,i_chunk)
      Antenna_SAI= 1000*Ant_Stations(i_ant,i_chunk) + Ant_IDs(i_ant,i_chunk)
      Do i_stat=1, Nr_UniqueStat      ! Get station number from the Unique_StatID list
          If(Unique_StatID(i_stat).eq. Station_ID) exit
      enddo
      FitDelay(i_ant)=Fit_TimeOffsetStat(i_stat)
      FitDelay(i_ant+1)=Fit_TimeOffsetStat(i_stat)
      If(Fit_AntOffset) then
         Do i_SAI=Tot_UniqueAnt(i_stat-1)+1,Tot_UniqueAnt(i_stat)      ! Get antenna number from the Unique_Antennas list
             If(Unique_SAI(i_SAI).eq. Antenna_SAI) exit
         enddo
         FitDelay(i_ant)=FitDelay(i_ant) + Fit_TimeOffsetAnt(i_SAI)
         FitDelay(i_ant+1)=FitDelay(i_ant+1) + Fit_TimeOffsetAnt(i_SAI+1)
      endif
   Enddo
End Subroutine GetInterfFitDelay
! ------------------------
Subroutine WriteDelChiPeak(i_chunk, DelChi,PartChiSq,PartChi2Int)
   use constants, only : dp
   use DataConstants, only : Station_nrMax, Ant_nrMax
   use ThisSource, only :  PeakNrTotal, PeakPos, Peak_eo, ChunkNr
   Use Interferom_Pars, only : IntfNuDim, IntFer_ant, Nr_IntFerMx, Nr_IntferCh, N_fit
   use Chunk_AntInfo, only : Ant_Stations, Ant_pos
   use Chunk_AntInfo, only : Unique_StatID, Nr_UniqueStat, Ant_IDs, Unique_SAI, Tot_UniqueAnt
   !use FitParams, only : N_FitPar, N_FitStatTim, FitParam, X_Offset
   !use FitParams, only : Fit_TimeOffsetStat, Fit_TimeOffsetAnt, Fit_AntOffset
   use StationMnemonics, only : Statn_ID2Mnem, Station_ID2Mnem
   use unque,only : Double_RI_sort
   implicit none
   real(dp), intent(in) :: DelChi(-N_fit:+N_fit,1:2*Nr_IntFerMx)
   integer, intent(in) :: i_chunk
   real(dp), intent(out) :: PartChiSq(1:Nr_IntFerMx)
   integer, intent(out) :: PartChi2Int(1:Nr_IntFerMx)
   Real(dp) :: ChiSq
   real(dp) :: PartChiSq_p(1:Nr_IntFerMx), PartChiSq_t(1:Nr_IntFerMx)
   !
   integer :: i, i_SAI, i_stat, i_ant, Station_ID,j_IntFer, Antenna_SAI
   !
   ChiSq=0.
   PartChiSq(:)=0
   PartChiSq_p(:)=0. ; PartChiSq_t(:)=0.
   Do j_IntFer=1,Nr_IntferCh(i_chunk)
      Do i=-N_fit, +N_fit
         PartChiSq_p(j_IntFer)=PartChiSq_p(j_IntFer) + DelChi(i,2*j_IntFer-1)*DelChi(i,2*j_IntFer-1)
         PartChiSq_t(j_IntFer)=PartChiSq_t(j_IntFer) + DelChi(i,2*j_IntFer)*DelChi(i,2*j_IntFer)
      EndDo
      PartChiSq_p(j_IntFer)=PartChiSq_p(j_IntFer)/(2*N_fit+1.)
      PartChiSq_t(j_IntFer)=PartChiSq_t(j_IntFer)/(2*N_fit+1.)
      PartChiSq(j_IntFer)=(PartChiSq_p(j_IntFer)+PartChiSq_t(j_IntFer))/2.
      PartChi2Int(j_IntFer)=j_IntFer
      ChiSq=ChiSq + PartChiSq(j_IntFer)
   Enddo
   ChiSq=ChiSq/Nr_IntferCh(i_chunk)
      !
   Do j_IntFer=1,Nr_IntferCh(i_chunk) ! do the printing
      i_ant=IntFer_ant(j_IntFer,i_chunk)
      Station_ID = Ant_Stations(i_ant,i_chunk)
      Antenna_SAI= 1000*Ant_Stations(i_ant,i_chunk) + Ant_IDs(i_ant,i_chunk)
      !Do i_stat=1, Nr_UniqueStat      ! Get station number from the Unique_StatID list
      !    If(Unique_StatID(i_stat).eq. Station_ID) exit
      !enddo
      If(PartChiSq(j_IntFer).gt. (2*ChiSq)) Then
         Write(2,*) j_IntFer, Statn_ID2Mnem(Station_ID), Antenna_SAI, PartChiSq(j_IntFer),&
               PartChiSq_p(j_IntFer), PartChiSq_t(j_IntFer)
      EndIf
   Enddo
   Call Double_RI_sort(Nr_IntferCh(i_chunk),PartChiSq,PartChi2Int)
   write(2,*) 'chi^2/ndf=',ChiSq
   Return
End Subroutine WriteDelChiPeak
