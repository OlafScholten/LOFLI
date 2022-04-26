    !include 'C:/OlafsUtil/Minimizer_ASA/asa047.f90'
    !Include '../../NumLib/LSQ/nl2sol.f90'
    Include 'nl2sol.f90'
!=================================
Subroutine FitCCorr(X)
!   Use method of steepest descent by calling routine PRAXIS
    !use constants, only : dp
!  v13: Introduce a dynamic search window based on the covariance matrix 'SpaceCov'
!        - Window is implemented as a parabolic multiplier normalized to unity and is zero
!           at a certain multiple of the search window. This factor is named 'SearchRangeFallOff',
!           implemented in 'ReImAtMax'.
!  v13: 'SpaceCov' is updated after each fit, either by the new covariance matrix, or by reducing the old one.
    ! v17 : offset in Nr_TimeOffset reduced by 1 in nvar
!
    use FitParams, only : FitParam, N_FitPar, N_FitStatTim, N_FitStatTim, Fit_PeakNrTotal, Nr_TimeOffset, N_FitPar_max
    use FitParams, only : CalcHessian, Sigma, SpaceCov, FitQual, ParamScaleFac, N_EffAnt, Max_EffAnt !, ImagingRun
    use DataConstants, only : ChunkNr_dim, Production
    use ThisSource, only : Nr_Corr, PeakNr, PeakNrTotal, CCorr_Err, Peak_eo, Safety
    Implicit none
    real ( kind = 8 ), intent(inout) :: X(N_FitPar_max)
    integer ( kind = 4 ) :: i, j, k
    !    real ( kind = 8 ), external :: ImageQuality  ! function that compares FitFunction to Data
    !  ImageQuality  calculates the goodness of fit for a given set of parameters (source pos & stat delays)
    real ( kind = 8 ) :: X0(N_FitPar_max) ! parameters that are optimized
    !
    !integer ( kind = 4 ), parameter :: v_dim = 20000 ! dim=93 + n*p + 3*n + p*(3*p+33)/2 ,n-meqn, p=nvar
    integer :: v_dim, i_chunk, error
    integer ( kind = 4 ) :: nvar    ! number of parameters
    integer ( kind = 4 ) :: meqn  ! Number of data points to fit
    integer ( kind = 4 ) iv(60+N_FitPar_max)
    external CompareCorrTime  ! subroutine that compares FitFunction to Data
    External JacobianCorrTime
    external ufparm  ! Dummy external routine
    integer ( kind = 4 ) uiparm(1)    ! Not really used
    real ( kind = 8 ),allocatable :: Jacobian(:,:)
    !real ( kind = 8 ) v(v_dim) ! dim=93 + n*p + 3*n + p*(3*p+33)/2 ,n-meqn, p=nvar
    real ( kind = 8 ),allocatable :: v(:) !  in NL@SOL: real ( kind = 8 ) v(93 + n*p + 3*n + (p*(3*p+33))/2)
    integer ( kind = 4 ) :: prin, NF
    Real(kind=8) :: XStore(1:3),Dist
    !
    !If(.not. Production) Call PrntFitPars(X)
    !
!stop
    i=N_FitPar-N_FitStatTim  ! number of parameters for a single source
    k=0
    If(i.gt.1) then
      If(any(FitParam(N_FitStatTim+1:N_FitStatTim+i)==4)) k=PeakNrTotal-Fit_PeakNrTotal ! for each pulse the timings need be fitted
    endif
    nvar=Nr_TimeOffset-1 + (N_FitPar-N_FitStatTim)*Fit_PeakNrTotal+k  ! total number of parameters that are fitted
    If(nvar.gt.N_FitPar_max) then
      write(*,*) 'Fitting stat: nvar=',nvar,', Meqn=',Meqn,', buffer size=',v_dim
      write(2,*) '*****Too many fitparameters=',nvar,'should not exceed ',N_FitPar_max
      write(2,*) 'breakdown, Nr_TimeOffset=',Nr_TimeOffset-1,'peak-related=',nvar-Nr_TimeOffset+1
      stop 'FitCCorr'
    endif
    !write(2,*) 'Nr_TimeOffset + (N_FitPar-N_FitStatTim)*Fit_PeakNrTotal+k',Nr_TimeOffset, (N_FitPar-N_FitStatTim),Fit_PeakNrTotal,k
    Meqn=0          ! total number of data-points that enter in the chi^2 search
    Do i_chunk=1, ChunkNr_dim
      i=MAX(Nr_Corr(0,i_chunk)-1,0)
      k=MAX(Nr_Corr(1,i_chunk)-1,0)
      If(Nr_Corr(0,i_chunk).eq.0 .and. PeakNr(0,i_chunk).ne.0) then
         write(2,*) '*****FitCCorr even',i_chunk,Nr_Corr(0,i_chunk),PeakNr(0,i_chunk),Nr_Corr(1,i_chunk),PeakNr(1,i_chunk)
         return
      endif
      If(Nr_Corr(1,i_chunk).eq.0 .and. PeakNr(1,i_chunk).ne.0) then
         write(2,*) '*****FitCCorr odd',i_chunk,Nr_Corr(1,i_chunk),PeakNr(1,i_chunk)
         return
      Endif
      Meqn = Meqn + i*PeakNr(0,i_chunk) + k*PeakNr(1,i_chunk)      ! number of equations
    enddo
    v_dim=93 + Meqn*nvar + 3*Meqn + nvar*(3*nvar+33)/2
    !write(*,*) 'Meqn=',Meqn,nvar,Nr_Corr(0,1),Nr_Corr(1,1),i,k
    !write(2,*) 'Nr_TimeOffset:',Nr_TimeOffset,N_FitPar,N_FitStatTim,Fit_PeakNrTotal,k
    IF(.not. Production) write(2,*) 'Fitting stat: nvar=',nvar,', Meqn=',Meqn,', buffer size=',v_dim
    !write(*,*) 'Meqn,nvar',Meqn,nvar
    allocate( Jacobian(Meqn,nvar) )
    allocate( v(v_dim) )
    !
    ParamScaleFac(:) = 1.
    !Do i= 1, N_FitStatTim
    !    ParamScaleFac(i) = 0.01  !  Bring all fit variables to similar order of magnitude
    !Enddo
    !Do i=1,N_FitPar-N_FitStatTim
    !    If(FitParam(N_FitStatTim+i).eq. 3) then
    !        k=N_FitStatTim+(i-1)*2*PeakNr
    !        ParamScaleFac(k+1:k+2*PeakNr) = 10
    !    Endif
    !enddo
    ! max step-length in Station-offset = 1[sample]; in x & y = 100[meter]; in z = 1000[meter]
    !
    X(1:nvar) = X(1:nvar) / ParamScaleFac(1:nvar) !  Bring all fit variables to similar order of magnitude
    !
    !write(2,*) 'X(1):', X(1)
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
    If(CalcHessian) then
      iv(14) = 1
      iv(15) = 2
      v(32) =1.d-8 ! is the relative function convergence tolerance
    endif
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
    FitQual=-1
   if((nvar .gt. 0) .and. (meqn .gt. nvar)) then   ! otherwise the system is underdetermined
         !NF=-1
         !call CompareCorrTime ( meqn, nvar, x, NF, v, uiparm, Jacobian, ufparm )
         !call nl2sno ( meqn, nvar, x, CompareCorrTime, iv, v, uiparm, Jacobian, ufparm )
         Call nl2sol ( meqn, nvar, x, CompareCorrTime, JacobianCorrTime, iv, v, uiparm, Jacobian, ufparm, error )
         If(error.ge.10) then
            write(2,*) 'parameter is NaN'
            error=10
            goto 9
         endif
         If(error.gt.0) then
            write(2,*) 'endless loop in NL2SOL truncated'
            goto 9
         endif
         !
         !Write(*,*) 'end nl2sol'
         FitQual=2*v(10)/(meqn-nvar)
         IF(.not. Production) write(2,"('Result, chi^2/ndf=',F9.2)", ADVANCE='NO') 2*v(10)/(meqn-nvar)
         IF(.not. Production) write(*,"('chi^2/ndf=',F9.2)") 2*v(10)/(meqn-nvar)
         IF(.not. Production) write(2,"(', # of fie calls=',i3,' and # of iterations=',i3,' Covariance at',i5)") iv(6),iv(31)
         sigma(:)=999.999
         N_EffAnt=-nvar
         Do k=1,PeakNrTotal
            N_EffAnt=N_EffAnt+COUNT(CCorr_Err(2:Nr_Corr(Peak_eo(k),1),k).lt.1.)
         Enddo
         Max_EffAnt= meqn-nvar
         k=NVar
         If(NVar.gt.3) k=3
         !If(NVar.gt.3) SpaceCov=SpaceCov/4.  ! Reduce the default errors by a factor 2, in case Hessian was not calculated
         !write(2,*) 'Fitter, SpaceCov(1:k,1:k)',SpaceCov(1:k,1:k)
         !SpaceCov(1:k,1:k)=SpaceCov(1:k,1:k)/9.  ! Reduce the default errors by a factor 3, in case Hessian was not calculated
         SpaceCov(1:k,1:k)=SpaceCov(1:k,1:k)/2.  ! Reduce the default errors by a factor sqrt(2),(=roughly ratio of distMax) in case Hessian was not calculated
         !Write(2,*) 'if iv(covmat): ',iv(26)
!           if iv(covmat) = -1, then the finite difference hessian was indefinite
!              if iv(covmat) = -2, then a successful finite differencing step could not be found for some component of x
         If( (iv(26) .gt.0) .and. (Fit_PeakNrTotal.eq.1) .and. (NVar.ge.3) ) then
            IF(.not. Production) write(2,*) 'Covariance matrix = chi^2/ndf * Hessian^{-1}'
            k=NVar
            If(NVar.gt.4) k=4
            Do i=1,k
               !write(2,"(10f10.4)") (V(iv(26)+k+i*(i-1)/2-1),k=1,i)
               IF(.not. Production) write(2,"(10g12.3)") V(iv(26)+i*(i-1)/2:iv(26)+i*(i+1)/2-1)
               Sigma(i)=sqrt(V(iv(26)+i*(i+1)/2-1))
               !If(Sigma(i) .gt. 1000.) Sigma(i)=999.999
            Enddo
            !IF(.not. Production .and. ImagingRun) then
            !   XStore(1:3)=X(1:3)
            !   Dist= sqrt(X(1)*X(1) + X(2)*X(2) + X(3)*X(3))
            !   NF=-1
            !   write(2,*) 'Source distance moved to ',Dist+2.,'[m]'
            !   X(1:3)=XStore(1:3)*(Dist+2.)/Dist
            !   !write(2,*) 'Source distance moved to height',XStore(3)+2.5,'[m]'
            !   !X(3)=XStore(3)+2.5
            !   call CompareCorrTime ( meqn, nvar, x, NF, v, uiparm, Jacobian, ufparm )
            !   write(2,*) 'Source distance moved to ',Dist-2.,'[m]'
            !   X(1:3)=XStore(1:3)*(Dist-2.)/Dist
            !   !write(2,*) 'Source distance moved to height',XStore(3)-2.5,'[m]'
            !   !X(3)=XStore(3)-2.5
            !   call CompareCorrTime ( meqn, nvar, x, NF, v, uiparm, Jacobian, ufparm )
            !   X(1:3)=XStore(1:3)
            !   write(2,*) 'Source distance moved back to ',Dist,'[m]'
            !Endif
            Do i=1,3
               Do j=1,i
               SpaceCov(i,j)=V(iv(26)+i*(i-1)/2+j-1)
               SpaceCov(j,i)=SpaceCov(i,j)
               Enddo
            Enddo
            !If(k.le.3) then  !  Assume that RefTimeErr is  fitted
            !   Do i=k,3
            !      SpaceCov(i,i)=1.D7 ! about 3 km
            !   Enddo
            !Endif
         Endif
         !If(Sigma(3).lt.99.)
         IF(.not. Production) then
            Write(2,"(A,2I6,A,F8.2,A)", ADVANCE='NO') 'N_EffAnt=',N_EffAnt,Max_EffAnt &
               ,', chi^2/ndf=',2*v(10)/(meqn-nvar), ', stations marked when timing StDev > 1.5'
               !,', Sigma=',Sigma(1:3)
            If(SpaceCov(1,1) .gt. 0.) then
               Write(2,"(A,G11.3)", ADVANCE='NO') ', SQRT{SpaceCov(i,i)}[km]:',sqrt(SpaceCov(1,1))/1000.
               If(SpaceCov(2,2) .gt. 0.) then
                  Write(2,"(G11.3)", ADVANCE='NO') sqrt(SpaceCov(2,2))/1000.
                  If(SpaceCov(3,3) .gt. 0.) Write(2,"(G11.3)", ADVANCE='NO') sqrt(SpaceCov(3,3))/1000.
               Endif
            Endif
            write(2,*) ' '
         EndIf
         !Write(2,*) 'SQRT{SpaceCov(i,i)}[km]:',(sqrt(SpaceCov(i,i))/1000.,i=1,3)
         !
         NF=-1
         call CompareCorrTime ( meqn, nvar, x, NF, v, uiparm, Jacobian, ufparm )
         X(1:nvar) = X(1:nvar) * ParamScaleFac(1:nvar)  !  Undo factor in COMPASS
         If(.not. Production) Call PrntFitPars(X)
   endif
9  Continue
    !
    Deallocate( Jacobian)
    Deallocate( v )
    !
    !Call PrntFitPars(X)
    Return
    !
End Subroutine FitCCorr
!=================================================
Subroutine JacobianCorrTime ( meqn, nvar, X_p, nf, Jacobian, uiparm, urparm, ufparm )
  implicit none
  real ( kind = 8 ), intent(in) :: X_p(1:nvar)
  integer ( kind = 4 ), intent(in) :: meqn
  integer ( kind = 4 ), intent(in) :: nvar
  integer ( kind = 4 ), intent(in) :: nf
  external ufparm
  integer ( kind = 4 ), intent(in) :: uiparm(*)
  real ( kind = 8 ), intent(in) :: urparm(*)
  real ( kind = 8 ), intent(out) :: Jacobian(meqn,nvar) ! in calling this routine this space has been reserved as part of v
  real ( kind = 8 ) :: R(meqn), Jacob(meqn,nvar)
  !
  !write(2,*) 'X_p',X_p(1:3)
  !write(*,*) 'calc Jac'
  Call CompareCorrTime ( meqn, nvar, X_p, nf, R, uiparm, Jacob, ufparm )
  !write(2,*) 'Jacob2',Jacob(2,:)
  !write(2,*) 'Jacob3',Jacob(3,:)
  !write(*,*) 'JacobN',meqn,Jacob(meqn,:)
  Jacobian(:,:)=Jacob(:,:)
  Return
  End Subroutine JacobianCorrTime
!=================================================
Subroutine CompareCorrTime ( meqn, nvar, X_p, nf, R, uiparm, Jacobian, ufparm )
    use constants, only : dp,sample,c_mps
    use DataConstants, only : Station_nrMax, Ant_nrMax, Production
    use DataConstants, only : Polariz
    use ThisSource, only : Nr_Corr, CCorr_max, CCorr_Err, PeakNrTotal, PeakPos, Peak_eo, ChunkNr, PeakRMS, PeakChiSQ
    use ThisSource, only : CorrAntNrs, T_Offset, SourcePos, RefAntErr, Peak_Offst, Dropped, StStdDevMax_ns
    use Chunk_AntInfo, only : Ant_Stations, Ant_pos, Ant_RawSourceDist
    use Chunk_AntInfo, only : Unique_StatID, Nr_UniqueStat, Ant_IDs, Unique_SAI, Tot_UniqueAnt
    use FitParams, only : ParamScaleFac, N_FitPar, N_FitStatTim, FitParam, X_Offset
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
    real ( kind = 8 ), intent(in) :: X_p(1:nvar)
  integer ( kind = 4 ), intent(in) :: meqn
  integer ( kind = 4 ), intent(in) :: nvar
  integer ( kind = 4 ), intent(in) :: nf
  real ( kind = 8 ), intent(out) :: R(meqn)
  external ufparm
  integer ( kind = 4 ), intent(in) :: uiparm(*)
  real ( kind = 8 ), intent(out) :: Jacobian(meqn,nvar)
  real ( kind = 8 ) :: X(nvar),D1,D2,Rns

   real ( kind = 8 ) ::   StatFineOff, T_shft, val
   Real(dp) :: ChiSq(1:Station_nrMax,1:PeakNrTotal), Ave(1:Station_nrMax,1:PeakNrTotal), RMS(1:Station_nrMax,1:PeakNrTotal)
   Real(dp) :: SumJac(1:Station_nrMax,1:PeakNrTotal), ChiSqWeight(1:Station_nrMax,1:PeakNrTotal)
   integer :: Cnt(1:Station_nrMax), CPE_eqn(Ant_nrMax),C_A(Ant_nrMax)  ! , Dropped(1:Station_nrMax)
   !
   integer :: i,j,k, i_SAI, j_corr, i_eqn, i_stat, i_ant, i_Peak, i_eo, i_chunk, Station_ID, i_ant1, i_xyz, i_StOff, Antenna_SAI
   real(dp) :: RDist, T_shft1
   integer, external :: XIndx
   logical :: prn, StCal
   Character(len=5) :: Station_Mnem
   Real(dp) :: IndxRefrac
   Real(dp), external :: RefracIndex
    !
    prn=.false.
    StCal=.false.  ! Calculate fit statistics per source
    if(NF.lt.0) prn=.true.
    if(NF.lt.0) StCal=.true.
    IF( Production) prn=.false.
    !if(NF.lt.0) FullAntFitPrn=.true.
    !write(2,*) 'NF=',NF
    !If(N_FitStatTim .gt. 0) prn=.true.
    !If(N_FitStatTim .gt. 0) FullAntFitPrn=.true.
    !write(2,*) 'NF',NF, X_p(1:nvar)
    !prn=.true.
    !
    !if(nvar.eq.1) write(*,*) 'CompareCorrTime',meqn, nvar, nf
    !write(2,*) 'X---Corr', nf, nvar, X_p(1:10),';', ParamScaleFac(1:10)
    X(1:nvar) = X_p(1:nvar) * ParamScaleFac(1:nvar)  !
    !write(2,*) 'X(1:nvar)',X(1:nvar)
    Call X2Source(X)
    !
    !If(nf.eq.1)         Call PrntFitPars(X)
    !write(2,"(A,10F11.2)") 'RefAntErr_CCT(i_Peak):',(RefAntErr(i_Peak),i_Peak=1,PeakNrTotal)
    !write(2,"(A,10F11.2)") 'X_CCT=',(X(i),i=1,10)
    !write(2,"(A,10F11.2)") 'X_CCT=',(X(i),i=11,nvar)
    !write(2,200) X
    !    write(*,"('Fitting-try:',15(G12.6,','))") X
    !
    i_eqn=0
    ChiSq=0.
    Ave=0.
    RMS=0.
    SumJac=0.
    ChiSqWeight(:,:)=0.
    Dropped = 0
    !if(nvar.eq.1) write(*,*) 'X_p',X_p
    Do i_Peak=1,PeakNrTotal
        Cnt=0
        i_eo=Peak_eo(i_peak)
        i_chunk=ChunkNr(i_peak)
        !
        Do j_corr=1,Nr_Corr(i_eo,i_chunk)
            i_ant=CorrAntNrs(j_corr,i_eo,i_chunk)
            Station_ID = Ant_Stations(i_ant,i_chunk)
            Antenna_SAI= 1000*Ant_Stations(i_ant,i_chunk) + Ant_IDs(i_ant,i_chunk)
            !write(*,*) 'j_corr=',j_corr,Antenna_SAI
            Do i_stat=1, Nr_UniqueStat      ! Get station number from the Unique_StatID list
                If(Unique_StatID(i_stat).eq. Station_ID) exit
            enddo
            !if(nvar.eq.1) write(*,*) j_corr,i_ant, Station_ID, Antenna_SAI
            !if(prn) write(2,*) 'Unique_StatID(i_stat):',Unique_StatID(i_stat),Station_ID,i_stat
            StatFineOff=Fit_TimeOffsetStat(i_stat)
            If(Fit_AntOffset) then
               Do i_SAI=Tot_UniqueAnt(i_stat-1)+1,Tot_UniqueAnt(i_stat)      ! Get antenna number from the Unique_Antennas list
                   If(Unique_SAI(i_SAI).eq. Antenna_SAI) exit
               enddo
               !If(FullAntFitPrn) write(2,*) 'StatFineOff',i_SAI,Antenna_SAI,StatFineOff,Fit_TimeOffsetAnt(i_SAI)
               StatFineOff=StatFineOff + Fit_TimeOffsetAnt(i_SAI)
            endif
            !
            i_StOff=0
            If(N_FitStatTim .gt. 0) then
                Do i=1,N_FitStatTim     ! When fitted, change fine-offset to fit-value i_StOff
                    If(FitParam(i) .eq. i_stat) then
                        If(Fit_AntOffset) then
                           i_StOff=X_Offset(i) + i_SAI-Tot_UniqueAnt(i_stat-1)-1 ! position of this parameter in array X
                           If(polariz) then  ! This needs to agree with the settings in X2Source
                              i_StOff=X_Offset(i) + (i_SAI-Tot_UniqueAnt(i_stat-1)-1)/2 ! position of this parameter in array X, count even only
                           EndIf
                           !write(2,*) 'i_StOff:',i_StOff,i,i_stat,Antenna_SAI, X_Offset(i), i_SAI, Tot_UniqueAnt(i_stat-1)
                        else
                           i_StOff=X_Offset(i)       ! Station offset for fitting station timing
                        endif
                        ! i_StOff is needed for calculating  Jacobian(i_eqn,i_StOff)=[d(R(i_ant,i_Peak)0/d(X(i_StOff))]
                        ! where 'X(i_StOff)' refers to a station or antenna timing off-set
                        ! This is finite only when 'i_StOff' corresponds to the station or antenna for which 'R' is calculated
                        !if(prn) write(2,*) 'i_StOff:',i_StOff,i,FitParam(i)
                        exit
                    endif
                Enddo
            Endif
            !
            Call RelDist(SourcePos(1,i_Peak),Ant_pos(1,i_ant,i_chunk),RDist)
            Rdist=Rdist - Ant_RawSourceDist(i_ant,i_chunk) + StatFineOff ! - INT(Peak_Offst(i_Peak))
            If(j_corr .eq. 1) then
               !T_shft1=Rdist + RefAntErr(i_Peak)  ! Before April 2022, works better than when correcting for the off-set
                T_shft1=Rdist - CCorr_max(1,i_Peak)  + RefAntErr(i_Peak)  ! April 2022
               ! One would expect that the self-correlation, i.e. CCorr_max(1,i_Peak), peaks at zero, however,
               !  this seems not the case for asymmetric pulses. For these asymmetric pulses the real part peaks at
               !  (close to) zero, as expected, but the imaginary part is not anti-symmetric w.r.t. zero. This appears
               !  to pull the max of the Hilbert transform away from zero, much to my surprise. To compensate for this
               !  the off-set is corrected for this asymmetry that should persist for all cross correlations. Nasty!!
               !If(FullAntFitPrn) write(2,*) 'i_Peak,j_corr,i_ant,Antenna_SAI, CCorr_max(j_corr,i_Peak):' &
               !   ,i_Peak,j_corr,i_ant,Antenna_SAI, CCorr_max(1,i_Peak), CCorr_max(2,i_Peak)
               i_ant1=i_ant
               D1=sum((SourcePos(:,i_Peak)-Ant_pos(:,i_ant,i_chunk))*(SourcePos(:,i_Peak)-Ant_pos(:,i_ant,i_chunk)))
               D1=sqrt(D1)     ! True distance from source to reference antena
            !
               cycle
            EndIf
            T_shft=Rdist - T_shft1 - T_Offset(j_corr,i_Peak) ! shift in cross-correlation peak due to changed source position while fitting
            ! where 'T_Offset' is the shift between the this and the reference antenna-timings which is
            !   already accounted for when calculating the cross correlations,
            !   which, of course, is based on the source positions before the coming fitting round. (set in GetCorrSingAnt)
            !
            i_eqn=i_eqn+1
            !if(nvar.eq.1) write(*,*) i_eqn
            R(i_eqn)= (CCorr_max(j_corr,i_Peak) - T_shft)/CCorr_Err(j_corr,i_Peak)
            If(CCorr_Err(j_corr,i_Peak).gt.1000) R(i_eqn)=0.  ! For excluded stations no peak-position is determined
            !
            !If(FullAntFitPrn) write(2,*) 'i_Peak,j_corr,i_ant,Antenna_SAI, CCorr_max(j_corr,i_Peak):' &
            !   ,i_Peak,j_corr,i_ant,Antenna_SAI, CCorr_max(j_corr,i_Peak)- T_shft
            !
            ChiSq(i_stat,i_Peak) = ChiSq(i_stat,i_Peak) + R(i_eqn)*R(i_eqn)
            Cnt(i_stat) = Cnt(i_stat) + 1
            !
            Jacobian(i_eqn,:) = 0.
            If(i_StOff.gt.0) Jacobian(i_eqn,i_StOff) = -1./CCorr_Err(j_corr,i_Peak)
            D2=sum((SourcePos(:,i_Peak)-Ant_pos(:,i_ant,i_chunk))*(SourcePos(:,i_Peak)-Ant_pos(:,i_ant,i_chunk)))
            D2=sqrt(D2)     ! True distance from source to present antena
            !
            val=0.
            !if(nvar.eq.1) write(*,*) Cnt(i_stat), i_stat
            IndxRefrac = RefracIndex( SourcePos(3,i_Peak) )
            Do i=1 , N_FitPar-N_FitStatTim !Jacobian(i_eqn,i) = d R(i_eqn) / d X(i)
                i_xyz=FitParam(N_FitStatTim+i)
                k=XIndx(i,i_Peak)
                If(i_xyz .eq. 4) then
                    Jacobian(i_eqn,k)=+1./CCorr_Err(j_corr,i_Peak) ! Due to RefAntErr(i_Peak)
                else
                    Jacobian(i_eqn,k)=(SourcePos(i_xyz,i_Peak)-Ant_pos(i_xyz,i_ant,i_chunk))/D2 &
                        - (SourcePos(i_xyz,i_Peak)-Ant_pos(i_xyz,i_ant1,i_chunk))/d1
                    Jacobian(i_eqn,k)=-Jacobian(i_eqn,k)*IndxRefrac /(c_mps*sample) ! Convert to units of [samples]
                    val=val+Jacobian(i_eqn,k)*Jacobian(i_eqn,k)*25.  ! 25 to convert from [samples^2] to [ns^2]
                    Jacobian(i_eqn,k)=Jacobian(i_eqn,k)/CCorr_Err(j_corr,i_Peak)
                endif
            enddo
            !if(nvar.eq.1) write(*,*) k
            !
            If(StCal) then
                If(CCorr_Err(j_corr,i_Peak) .gt. 10.) then  ! modified May 16, 2021 to allow double peaks in real coorelation with reduced weight
                    Dropped(i_stat,i_Peak)= Dropped(i_stat,i_Peak) + 1
                    !write(2,*) 'CCorr_Err(j_corr,i_Peak)=',CCorr_Err(j_corr,i_Peak),j_corr,i_Peak,i_stat,Dropped(i_stat,i_Peak)
                endif
                CPE_eqn(j_corr)=i_eqn
                C_A(j_corr)=i_ant
                Rns = R(i_eqn)*5.*CCorr_Err(j_corr,i_Peak) ! put this in [ns]
                Ave(i_stat,i_Peak)=Ave(i_stat,i_Peak) + Rns ! put this in [ns]
                RMS(i_stat,i_Peak)=RMS(i_stat,i_Peak) + Rns*Rns ! put this in [ns]
                SumJac(i_stat,i_Peak)=SumJac(i_stat,i_Peak)+ val
                ChiSqWeight(i_stat,i_Peak)=ChiSqWeight(i_stat,i_Peak) +0.04/(CCorr_Err(j_corr,i_Peak)*CCorr_Err(j_corr,i_Peak))
            endif
            !
            !Call spline_cubic_val ( 2*Safety+1, t_CCorr(-Safety), &
            !    CCorr(-Safety,j_corr,i_Peak,i_eo), CCorr_pp(-Safety,j_corr,i_Peak,i_eo), &
            !    T_shft, val)
            !If(prn) write(2,"(A,I3,A,I2,I2,A,f7.2,A,F6.4)") 'j_corr=',j_corr,', i_Peak=',i_Peak,i_eo,&
            !    ', T_shft=', T_shft,', correlation=', val
        enddo ! j_corr=1,Nr_Corr(i_eo)
        If(StCal) then
            !PeakChiSQ(i_Peak)= sum(ChiSq(:,i_Peak))/(sum(Cnt(:))-sum(Dropped(:,i_Peak)))
            PeakChiSQ(i_Peak)= sum(ChiSq(:,i_Peak))/sum(ChiSqWeight(:,i_Peak))
            PeakRMS(i_Peak)= sqrt(sum(RMS(:,i_Peak))/sum(Cnt(:)))
        Endif
        If(prn) then
            write(2,"(A,I3,I2,A,I7,A, 2f7.2,A,3(F11.2,','),A,F7.3)", ADVANCE='NO') 'i_Peak=',i_Peak,i_eo, &
                ', PeakPos=',PeakPos(i_Peak),', Chi^2/DegrF=', PeakChiSQ(i_Peak), &
                sum(ChiSq(:,i_Peak))/(sum(Cnt(:))-sum(Dropped(:,i_Peak))), &
                ', source position:',SourcePos(:,i_Peak),' RefAntTimeErr:',RefAntErr(i_Peak)
            If(FullAntFitPrn) then
                Station_ID = -1 ! Ant_Stations(i_ant)
                Do j_corr=2,Nr_Corr(i_eo,i_chunk)
                    i_ant=C_A(j_corr)
                    i=CPE_eqn(j_corr) ! should not be named i_eqn
                    If(Station_ID .ne. Ant_Stations(i_ant,i_chunk)) then  ! continue printing
                        Station_ID = Ant_Stations(i_ant,i_chunk)
                        Call Station_ID2Mnem(Station_ID,Station_Mnem)
                        write(2,"(/,1x,A5,';')", ADVANCE='NO')   Station_Mnem
                    endif
                    val=R(i)*5.*CCorr_Err(j_corr,i_Peak)  !  [ns]
                    If(abs(val) .gt. 10.) then
                        write(2,"(i2.2,'*',F5.0,'*;')", ADVANCE='NO') Ant_IDs(i_ant,i_chunk),val
                    else
                        write(2,"(i2.2,F7.2,';')", ADVANCE='NO') Ant_IDs(i_ant,i_chunk),val
                    endif
                enddo
            endif
            write(2,*)
            write(2,"(A)", ADVANCE='NO') 'Stat Nr ='
            Do i_stat=1, Station_nrMax
                If(Cnt(i_stat) .ne. 0) then
                    Call Station_ID2Mnem(Unique_StatID(i_stat),Station_Mnem)
                    write(2,"(1x,A5,';')", ADVANCE='NO') &
                        Station_Mnem
                endif
            enddo
            write(2,*) '-----'
            write(2,"(A)", ADVANCE='NO') 'Dropped#='
            Do i_stat=1, Station_nrMax
                If(Cnt(i_stat) .ne. 0) then
                    write(2,"(I3,'/',i2,';')", ADVANCE='NO') &
                        Dropped(i_stat,i_Peak),Cnt(i_stat)
                endif
            enddo
            write(2,*) '-----'
            !write(2,"(A)", ADVANCE='NO') 'Chi^2/df='
            !Do i_stat=1, Station_nrMax
            !    If(Cnt(i_stat) .ne. 0) then
            !        write(2,"(F6.1,';')", ADVANCE='NO') &
            !            ChiSq(i_stat,i_Peak)/ChiSqWeight(i_stat,i_Peak)  ! Cnt(i_stat)
            !    endif
            !enddo
            !write(2,*) '-----'
            !write(2,"(A)", ADVANCE='NO') 'RMS [ns]='
            !Do i_stat=1, Station_nrMax
            !    If(Cnt(i_stat) .ne. 0) then
            !        write(2,"(F6.1,';')", ADVANCE='NO') &
            !            sqrt(RMS(i_stat,i_Peak)/Cnt(i_stat))
            !    endif
            !enddo
            !write(2,*) '-----'
            write(2,"(A)", ADVANCE='NO') 'Avrg[ns]='
            Do i_stat=1, Station_nrMax
                If(Cnt(i_stat) .ne. 0) then
                    write(2,"(F6.1,';')", ADVANCE='NO') &
                        (Ave(i_stat,i_Peak))/Cnt(i_stat)
                endif
            enddo
            write(2,*) '-----'
            !write(2,"(A)", ADVANCE='NO') 'DAve[ns]='
            !Do i_stat=1, Station_nrMax
            !    If(Cnt(i_stat) .ne. 0) then
            !        write(2,"(F6.1,';')", ADVANCE='NO') &
            !           -Ave(i_stat,i_Peak)/Cnt(i_stat)+SIGN(sqrt(RMS(i_stat,i_Peak)/Cnt(i_stat)),Ave(i_stat,i_Peak))
            !    endif
            !enddo
            !write(2,*) '-----'
            write(2,"(A)", ADVANCE='NO') 'StD [ns]='
            Do i_stat=1, Station_nrMax
                If(Cnt(i_stat) .ne. 0) then
                    val=SQRT(RMS(i_stat,i_Peak)/Cnt(i_stat) - (Ave(i_stat,i_Peak)/Cnt(i_stat))**2)
                    If(val.gt.StStdDevMax_ns) Then
                     Dropped(i_stat,i_Peak)= Dropped(i_stat,i_Peak) + 1
                    EndIf
                    If(val.lt.1.5) Then
                        write(2,"(F6.1,';')", ADVANCE='NO') val
                    Else
                        write(2,"('**',F4.1,';')", ADVANCE='NO') val
                    EndIf
                endif
            enddo
            !write(2,*) '-----'
            !write(2,"(A)", ADVANCE='NO') 'RMS(Jac)='
            !Do i_stat=1, Station_nrMax
            !    If(Cnt(i_stat) .ne. 0) then
            !        write(2,"(F6.2,';')", ADVANCE='NO') &
            !            sqrt(SumJac(i_stat,i_Peak)/Cnt(i_stat))
            !    endif
            !enddo
            write(2,"(A,/)") '-----'
        endif ! prn
        !if(nvar.eq.1) write(*,*) 'SourcePos(1,i_Peak)',SourcePos(1,i_Peak),i_peak
    Enddo  ! i_Peak=1,PeakNr
    !    if(nvar.eq.1) write(*,*) 'R',R
    !
    !val=sum(ChiSq)
    If(prn) Write(2,*) ' total chisq ',sum(ChiSq),(sum(ChiSq(:,i_Peak)),i_Peak=1,PeakNrTotal)
    !stop  'CompareCorrTime'
    !
    ! If(N_FitStatTim .gt. 0) stop
    !OPEN(UNIT=4,STATUS='unknown',FILE=trim(FileFitResult)//'.dat' ) !'plot/FitResult.dat')
    !write(4,"('!CoreDist[m], phi[rad]',3x,'I',11x,'I_calc',8x,'sigma_I',7x,&
    !    'Q',9x,'Q_calc',4x,'sigma_Q',3x,'U',9x,'U_calc',4x,'sigma_U',3x,'V',9x,'V_calc',4x,'sigma_V')")
    !close(unit=4)
    return
End Subroutine CompareCorrTime
!==================================================
Subroutine ufparm ( meqn, nvar, x )
!*****************************************************************************80
!! UFPARM is a user-supplied external routine.
!
!  Discussion:
!    The name of the routine, the argument list, and even whether
!       it is a function or subroutine, are left to the user.
!    NL2SOL simply passes the external reference from the calling
!       program through to the residual and jacobian routines.
!    If the user has no need for this facility, then a dummy
!       routine like this one may be used.
!
!  Modified:
!    07 February 2003
!
!  Parameters:
!    Input, integer ( kind = 4 ) MEQN, the number of functions.
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!    Input, real ( kind = 8 ) X(NVAR), the current value of the variables.
!
  implicit none
  integer ( kind = 4 ) meqn
  integer ( kind = 4 ) nvar
  real ( kind = 8 ) x(nvar)
  return
End Subroutine ufparm
