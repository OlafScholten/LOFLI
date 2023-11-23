!=================================
Subroutine MDD_Fitter(FitQual,error)
!   Use method of steepest descent to tit the spectra in all antennas as the superposition of a number of dipole-delta sources
!     A dipole-delta (DD) sources in this context is a point-impulsive source (zero extent in space and time) emitting like an ideal dipole.
!
!  Input:
!     - suggested DD parameters (number DD, each [position, time, direction})
!
!  Output:
!     - best fitting DD parameters (number DD, each [position, time, direction})
!
!  Procedural steps:
!   1)
!
!
   use constants, only : dp
   Use MDD_Pars, only : NSrc, dCenT, dCenLoc, T_Range, DDChiSQ, MDDLabel
   Use MDD_Pars, only : WORK, ReWORK, LWORK, G_IR, Xscale, RealDip   ! Used in CompareMDD
   Use Interferom_Pars, only : IntFer_ant, Nr_IntFerCh, Nr_IntFerMx
   use DataConstants, only : RunMode
    Implicit none
    Integer, intent(out) :: error
    Real(dp), Intent(OUT) :: FitQual
    integer ( kind = 4 ) :: i, j, k
    integer ( kind = 4 ) :: nvar    ! number of parameters
    integer ( kind = 4 ) :: meqn  ! Number of data points to fit
    external ufparm  ! Dummy external routine
    integer ( kind = 4 ) :: uiparm(1)    ! Not really used
    real ( kind = 8 ) :: urparm(1)    ! Not really used
    integer ( kind = 4 ) :: NF
    real ( kind = 8 ),allocatable :: v(:) !  in NL@SOL: real ( kind = 8 ) v(93 + n*p + 3*n + (p*(3*p+33))/2)
    integer ( kind = 4 ), allocatable :: iv(:)
    real ( kind = 8 ), allocatable :: X(:)
    external CompareMDD  ! subroutine that compares FitFunction to Data
    Integer :: m, IntfBase, v_dim
    !Real(dp) :: Chi2_start !, Chi2StopRatio=1.0
    !
    Meqn=T_Range*Nr_IntFerCh(1)*2*2   ! total number of data-points that enter in the chi^2 search one *2 for 2 pol directions, one *2 for real/imag parts
    nvar=4*NSrc
    v_dim=93 + Meqn*nvar + 3*Meqn + nvar*(3*nvar+33)/2
    write(2,*) 'Fitting statistics: nvar=',nvar,', Meqn=',Meqn,', buffer size=',v_dim
    allocate( v(v_dim) )
    allocate( Iv(60+nvar) )
    allocate( X(nvar) )
    Allocate( G_IR(-T_Range/2:T_Range/2,1:NSrc,1:Nr_IntFerCh(1)) )
    !
    Do m=1,Nsrc
      X(4*m-3)=dCenT(m)/Xscale
      X(4*m-2:4*m)=dCenloc(1:3,m)/Xscale
    EndDo
    !
    MDDLabel(5:5)=achar(iachar("a"))
    NF=-2
    call CompareMDD ( meqn, nvar, x, NF, v, uiparm, urparm, ufparm )  ! prepare allocation workspace
    !
    !Chi2_start=SUM(DDChiSQ(1:NSrc,1))/NSrc
    write(2,*) 'Initial values before time & position search'
    !write(2,*) nf,' mean chi^2=',SUM(DDChiSQ(1:Nr_IntFerCh(1),1))/Nr_IntFerCh(1), &
    !  SUM(DDChiSQ(1:Nr_IntFerCh(1),2))/Nr_IntFerCh(1)
    Call MDD_write()
    !Time_width=2.  ! [samples]
    !Space_Spread(1:3)=0. ! Time_width*(/ 1.d0,1.d0,3.d0 /)     ![m]
    !========================================================================
    !write(2,*) 'stop with MDD_Fitter ==============================================='
    !Goto 8
    !========================================================================
1   Continue
    !
    call dfault( iv, v)
    write(2,*) v(41),v(38),v(40), v(37), v(nvar+87), v(39)
    iv(1) = 12 ! 12= do not call dfault again
    iv(14) = 0 ! 1: means print a covariance matrix at the solution.
    iv(15) = 0 ! if = 1 or 2, then a finite difference hessian approximation h is obtained.
               ! if =0 nothing is done
               ! if positive: with step sizes determined using v(delta0=44), a multiplicative factor)
               ! If negative: then only function values are used with step sizes determined using v(dltfdc=40)
    iv(19) = 0 ! controls the number and length of iteration summary lines printed
    iv(21) = 0 !2 ! is the output unit number on which all printing is done.
    iv(22) = 0 !1 ! print out the value of x returned (as well as the corresponding gradient and scale vector d).
    iv(23) = 1 ! print summary statistics upon returning.
    iv(24) = 0 ! print the initial x and scale vector d
    v(40) =1.40d-1 ! the step size used when computing the covariance matrix when iv(covreq=15) = -1 or -2, step size = v(dltfdc=40) * max(abs(x(i)), 1/d(i))
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
   iv(18)=50
   !  iv(26) if (iv(covmat) is positive, then the lower triangle of the covariance matrix is stored rowwise in v starting at
   !             v(iv(covmat)).  if iv(covmat) = 0, then no attempt was made.
    ! iv(nfcov).... iv(40) is the number of calls made on calcr when trying to compute covariance matrices.
    ! iv(ngcall)... iv(30) is the number of gradient evaluations (calls on
    !             calcj) so far done (including those used for computing
    !             the covariance).
    ! iv(ngcov).... iv(41) is the number of calls made on calcj when
    !             trying to compute covariance matrices.
    ! iv(niter).... iv(31) is the number of iterations performed.
    !
    v(32) =1.d-3 ! v(rfctol)... v(32) is the relative function convergence tolerance.
!             if the current model predicts a maximum possible function
!             reduction (see v(nreduc)) of at most v(rfctol)*abs(f0) at
!             the start of the current iteration, where  f0  is the
!             then current function value, and if the last step attempt-
!             ed achieved no more than twice the predicted function
!             decrease, then nl2sol returns with iv(1) = 4 (or 5).
!             default = max(10**-10, machep**(2/3)), where machep is
!             the unit roundoff.
!   write(2,*) 'default v(33):',v(33)
!    v(33) =0.01/Xscale ! v(xctol).... v(33) is the x-convergence tolerance.  if a Newton step
!             (see v(nreduc)) is tried that has v(reldx) <= v(xctol)
!             and if this step yields at most twice the predicted func-
!             tion decrease, then nl2sol returns with iv(1) = 3 (or 5).
!             (see the description of v(reldx) below.) default = machep**0.5, where machep is the unit roundoff.
!   write(2,*) 'set to v(33):',v(33)
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
         !v(39)=0.1
         !v(87:87+nvar*nvar)=1.
         !write(2,*) v(41),v(38),v(40), v(37), v(nvar+87), v(39)
         error=0
         call nl2sno ( meqn, nvar, X, CompareMDD, iv, v, uiparm, urparm, ufparm, error )
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
         FitQual=2*v(10)/(meqn-nvar)
          write(2,"('Result, chi^2/ndf=',F12.5)", ADVANCE='NO') FitQual
          write(*,"('chi^2/ndf=',F12.5)") FitQual
          write(2,"(', # of fie calls=',i3,' and # of iterations=',i3)") iv(6),iv(31)
   endif
   NF=-1
    MDDLabel(5:5)=achar(iachar("a")+1)
   call CompareMDD ( meqn, nvar, x, NF, v, uiparm, urparm, ufparm )
   Call MDD_write()
9  Continue
    !
    If(error.gt.0) Then
      goto 8
    EndIf
8  Continue
    If(RealDip) Then
      Deallocate( v, Iv, X, ReWork, G_IR )
    Else
      Deallocate( v, Iv, X, Work, G_IR )
    EndIf
    LWORK=-1
    !
    Return
    !
End Subroutine MDD_Fitter
!=================================================
!=================================================
Subroutine CompareMDD( meqn, nvar, X, nf, R, uiparm, urparm, ufparm )
!*****************************************************************************80
!
!  Discussion:
!    This routine is required for runninf n2sol chi-2 fittingantennas
!    Given the value of the vector X (=the positions parameters for all DD source), this routine computes the value of
!     the residues R used for chi^2 fitting the spectra in all antennas.
!     The residues are calculated for gain_multiplied E-fields at the antenna, R[s,ant,k]=(E_measured[s,nt,k]-A_DDmodel[s,nt,k]),
!     where the E-fields from the data are calculated for the central location only, multiplied by the gain to dissolve
!     possible end-of-spectrum divergencies, for all samples [s], antennas [ant], and polarizations [k].
!    Much of the nomenclature ad the structure follows that of EI-Interferometry.
!    All calculations in this routine are made in the local time, where E-fields are retarded from antenna to
!     local time, t_l, where t_l=0 corresponds to the central source time at the center of the source region.
!    All antenna functions are unfolded from the data w.r.t. the central source location implementing an overal gain (=sensitivity).
!  Author:
!    Olaf Scholten, Janary 2023
!
!  Constans (for this routine):
!     - Weight factors for antennas and may depend on polarisation, W[ant] == Weight(1:Nr_IntFer),
!        calculated in first part of "EISetupSpec(,,,)" in 'EIOption.f90'
!     - E[t_c,ant,k] where [k]=[theta,phi] polarization directions and time shifted to the center of the source region,
!        theta[ant,i]  & phi[ant,i] i=Cartesian coordinates
!     - G(t) - upsampled pulse-response of system due to the frequency dependence of the gain.
!        Note that we compare gain-multiplied E-fields.
!
!  Intut:
!     -X(1:nvar)=SourceStructure(0:3,1:M)  ! time and position w.r.t. (C_Tms,CenLoc), units=samples,m

!  Output:
!     - R[s,ant,k]=R(meqn)
!
!  Procedure:
!     - First steps (0....4) are to construct the dipole vectors for the M sources
!        0) G[m,ant,t_c]=G(t_c-t_m  +  \hat{r}_{as}\cdot \vec{x}_m/c) / R_{as}
!           taking care of time shift and time-grid
!        1) A_{im,jn}= \sum_{a} W(a) [ [\sum_{k} \hat{r}_{ak,i} \hat{r}_{ak,j}] [\sum_{t_c} G[n,ant,t_c]  G[m,ant,t_c]] ]
!        In the process A^-1 is calculated where A = [3*M,3*M] matrix and may have zero eigenvalues
!         USe DSYEVR from the LAPACK library; in the linker add : -llapack -lblas
!           check using     find / -xdev -name *lapack*
!           gives:  /usr/lib/x86_64-linux-gnu/lapack/liblapack.a
!         see:  https://netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen_gaeed8a131adf56eaa2a9e5b1e0cce5718.html#gaeed8a131adf56eaa2a9e5b1e0cce5718
!         sample application from https://www.intel.com/content/www/us/en/develop/documentation/onemkl-lapack-examples/top/least-squares-and-eigenvalue-problems/symmetric-eigenproblems/syevr-function/dsyevr-example/dsyevr-example-fortran.html
!        Alternatively A x=B can be solved with
!           dsysv (UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO)
!         or with
!           dsysvx (FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, LWORK, IWORK, INFO)
!         to provide error estimates.
!        2) What to do for (near zero) eigenvalues?
!        3) F[i,m]= \sum_{ant}  w[ant] \left[\sum_{t_c} \sum_{k} \hat{r}_{ak} E[t_c,ant,k] \, G[m,ant,t_c] \right]_i
!        4) I[i,m] = A^-1  F
!     - Calculate R[s=t_c,ant,k] = W(ant) (E[t_c,ant,k]- \sum_m  \left( \hat{r}_{ak}\cdot \vec{I}_m )\right) G[m,ant,t_c]
!
!
!        above taken from https://www.intel.com/content/www/us/en/develop/documentation/onemkl-lapack-examples/top/least-squares-and-eigenvalue-problems/symmetric-eigenproblems/syevr-function/dsyevr-example/dsyevr-example-fortran.html
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
!---------------------------------------------------------------------------------------
   use constants, only : dp,sample,c_mps, c_l  ! speed of light in Vacuum [m/nsec]
   use DataConstants, only : DataFolder, OutFileLabel
    use DataConstants, only : Ant_nrMax
!    use ThisSource, only : Nr_Corr, CCorr_max, CCorr_Err
!    use ThisSource, only :  SourcePos, PeakNrTotal, PeakPos, ChunkNr, PeakChiSQ, ExclStatNr
 !   use Interferom_Pars, only : Cnu_p0, Cnu_t0, Cnu_p1, Cnu_t1,  N_fit, AntPeak_OffSt, Chi2pDF
   Use Interferom_Pars, only : Nr_IntferCh ! the latter gives # per chunk
    !use Chunk_AntInfo, only : Ant_Stations, Ant_pos
    !use Chunk_AntInfo, only : Unique_StatID, Nr_UniqueStat, Ant_IDs, Unique_SAI, Tot_UniqueAnt
   ! use FitParams, only : N_FitPar, N_FitStatTim, FitParam, X_Offset
   ! use FitParams, only : Fit_TimeOffsetStat, Fit_TimeOffsetAnt, Fit_AntOffset, FullAntFitPrn
    !use StationMnemonics, only : Station_ID2Mnem
   Use MDD_Pars, only : Weight, SQWght, CTime_p, CTime_t, Vec_p, Vec_t, Vec_l, AntSourceD  ! Input
   Use MDD_Pars, only : G_dim, G_ns !, IRFW_s
   Use MDD_Pars, only : NSrc, dCenT, dCenloc, T_Range, MDDLabel, MDDnuDim, Bias_inv, DDXcorr_depl, RealDip  !  Input
   Use MDD_Pars, only : WORK, ReWORK, LWORK,G_IR, Xscale  ! Local workspace
   Use MDD_Pars, only : DDChiSQ  ! Output
   Use MDD_Pars, only : DelDip  !  Output
    use GLEplots, only : GLEplotControl
   implicit none
   real ( kind = 8 ), intent(in) :: X(1:nvar)
   integer ( kind = 4 ), intent(in) :: meqn
   integer ( kind = 4 ), intent(in) :: nvar
   integer ( kind = 4 ), intent(in) :: nf
   real ( kind = 8 ), intent(out) :: R(meqn)
   external ufparm
   integer ( kind = 4 ), intent(in) :: uiparm(*)
   real ( kind = 8 ), intent(out) :: urparm(*)
   !
   Integer :: Nt_l  ! T_Range should be an odd number
   Integer ::  i, j, m, n, j_IntFer, idt, it_G, it_l, i_chunk=1, i_s, Nr_IntFer
   Real(dp) :: Power_p(1:Ant_nrMax), Power_t(1:Ant_nrMax)
   Real(dp) :: dt, ddt, D, W
   Real(dp), save ::  X_old(1:100), Ch1, CH1_old=-1.d0, Ch2, CH2_old=-1.d0
   Integer, save :: NF_old=0
   Complex(dp) ::  G, CE_MDD_p(1:Ant_nrMax), CE_MDD_t(1:Ant_nrMax), E_p, E_t, Gp, Gt
   Real(dp) :: St_MDD(1:4), St_dat(1:4), St_chi2(1:4)
   Integer :: INFO, NRhs, IPIV(1:3*Nvar/4)
   Character*1 :: UPLO
   Complex(dp) :: A(1:3*Nvar/4,1:3*Nvar/4), F(1:3*Nvar/4,1)
   Real(dp) :: ReA(1:3*Nvar/4,1:3*Nvar/4), ReF(1:3*Nvar/4,1)
   Complex(dp) :: Ck(1:3*Nvar/4), B(1:3*Nvar/4)
   Logical :: MakeCurtain=.false.
   Character(len=20) :: Txt
   character(len=20) :: FMT
   CHARACTER(LEN=60) :: GLE_file
   !Real(dp), External :: CalcStokes
   !
   MakeCurtain=.false.
   if(NF.eq.-1) MakeCurtain=.true.
   if(NF.eq.-2) MakeCurtain=.true.
   !
   NSrc=nvar/4
   Nr_IntFer=Nr_IntFerCh(1)
   Nt_l=T_Range/2  ! T_Range should be an odd number
   i=1
   Do m=1,Nsrc
      dCenT(m)=X(i)*Xscale
      dCenloc(1:3,m)=X(i+1:i+3)*Xscale
      i=i+4
   Enddo
   If(NF.ge.0) Then
      !write(2,"(A,2I5,100F7.2)") 'CompareMDD:',nvar, nf !,X(1:nvar)
      !write(2,"(100F9.6)") (X_old(1:nvar)-X(1:nvar))*Xscale
      !write(2,"(A,10G13.6)") 'change in pars:',X_old(1:8)-X(1:8)
      X_old(1:nvar)=X(1:nvar)
   EndIf
   !===============================================
   !     - First steps (0....4) are to construct the dipole vectors for the M sources
   !        0) G[m,ant,t_c]=G(t_c-t_m  +  \hat{r}_{as}\cdot \vec{x}_m/c) / R_{as}
   !           taking care of time shift and time-grid
   !=======================================================
   Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
      !i_ant=IntFer_ant(j_IntFer,i_chunk)  ! Refers to even antenna
      D=AntSourceD(j_IntFer)
      Do m=1,Nsrc
         dt= -(dCenT(m)*5. +SUM(Vec_l(1:3,j_IntFer)*dCenloc(1:3,m))/c_l)  ! in [ns]  ! Should be +SUM as Ras is (source-ant); further source=later time
         idt=FLOOR(dt)
         ddt=dt - idt  ! always betwee (0,1)
         !If(j_IntFer.eq.1) write(2,*) m,'dt',dt,idt,ddt, G_dim/2.
         Do it_l=-Nt_l, Nt_l ! local in [samples]
            it_G=G_dim/2.+ it_l*5 + idt  ! time in [ns]
            If(it_G.ge.G_dim .or. it_G.le.1) then
               G=0.
            Else
               G=((1.-ddt)*G_ns(it_G) + ddt*G_ns(it_G+1))/D
            EndIf
            !If(j_IntFer.eq.1) write(2,*) m,'G',G_ns(it_G),it_G,G
            G_IR(it_l,m,j_IntFer)=G
         EndDo
      Enddo
   EndDo
   !write(2,*) 'G_ns:', G_IR(1:15,1,1)
   !write(2,*) 'ABS(G_ns):', ABS( G_IR(1:15,1,1) )
   !write(2,*) 'G_IR(it_l,m,j_IntFer):', G_IR(-Nt_l:15-Nt_l,1,1)
   !write(2,*) 'ABS(G_IR(it_l,m,j_IntFer)):', ABS( G_IR(-Nt_l:15-Nt_l,1,1) )*AntSourceD(1)
   !===============================================
   !     - First steps (0....4) are to construct the dipole vectors for the M sources
   !        1a) A_{im,jn}= \sum_{a} W(a) [ [\sum_{k} \hat{r}_{ak,i} \hat{r}_{ak,j}] [\sum_{t_c} G[n,ant,t_c]  G[m,ant,t_c]] ]
   !=======================================================
   A(:,:)=0.
   Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
      !i_ant=IntFer_ant(j_IntFer,i_chunk)  ! Refers to even antenna
      W=Weight(j_IntFer)
      Do m=1,Nsrc
         Do n=1,Nsrc
            G=SUM(G_IR(:,n,j_IntFer)*CONJG(G_IR(:,m,j_IntFer)) )
            If(n.ne.m .and. (DDXcorr_depl.gt.0.d0)) G=G*(1.d0-DDXcorr_depl)
            !If(j_IntFer.eq.1) write(2,*)'A;m,n,G', m,n,G
            Do i=1,3
               Do j=1,3
                  A(3*m+i-3,3*n+j-3)=A(3*m+i-3,3*n+j-3) &
                     +G*W*(Vec_p(i,j_IntFer) * Vec_p(j,j_IntFer) + Vec_t(i,j_IntFer)* Vec_t(j,j_IntFer))
               EndDo
            EndDo ! i=1,3
         EndDo !  n=1,Nsrc
      EndDo !  m=1,Nsrc
      !write(2,*) 'vec_p,t:',j_IntFer, Vec_p(1:3,j_IntFer), Vec_t(1:3,j_IntFer)
   EndDo ! j_IntFer=1,Nr_IntFer
   Do i=1,3*Nsrc
      A(i,i)=A(i,i)+Bias_inv
      !write(2,*) 'A(i,1:3)',i,A(i,1:12)
   EndDo
   !===============================================
   !     - First steps (0....4) are to construct the dipole vectors for the M sources
   !        1b)=3) F(im)=\sum_{a} w_{a} [\sum_{t_c} \sum_{k} \hat{r}_{ak} E[t_c,ant,k] G[m,ant,t_c] ]
   !=======================================================
   F(:,1)=0.
   Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
      !i_ant=IntFer_ant(j_IntFer,i_chunk)  ! Refers to even antenna
      W=SQWght(j_IntFer)  ! Ctime_ are also weighted by SQWght
      Do m=1,Nsrc
         Gp=SUM( CONJG( G_IR(-Nt_l:Nt_l,m,j_IntFer) ) * ( CTime_p(MDDnuDim-Nt_l:MDDnuDim+Nt_l,j_IntFer) ) )
         Gt=SUM( CONJG( G_IR(-Nt_l:Nt_l,m,j_IntFer) ) * ( CTime_t(MDDnuDim-Nt_l:MDDnuDim+Nt_l,j_IntFer) ) )
         !If(j_IntFer.eq.1) write(2,*)'F;m,G', m,Gp,Gt
         Do i=1,3
            F(3*m+i-3,1)=F(3*m+i-3,1) + W*(Vec_p(i,j_IntFer) * Gp + Vec_t(i,j_IntFer)* Gt )
         EndDo ! i=1,3
      EndDo !  m=1,Nsrc
      !write(2,*) 'F(1:3):', j_IntFer, F(1:3,1), Gp, Gt
   EndDo ! j_IntFer=1,Nr_IntFer
   !write(2,*) 'F(1:3):', j_IntFer, F(1:3,1),W
   B(:)=F(:,1)
!
   !===============================================
   !     - First steps (0....4) are to construct the dipole vectors for the M sources
   !        1c)=4) Solve A X = F
   !         ZHESV computes the solution to a complex system of linear equations
   !           A * X = B,
   !         where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS matrices.
   ! https://netlib.org/lapack/explore-html/d3/d9d/group__complex16_h_esolve_gac2def3254215d3a7c56bc162483193d4.html#gac2def3254215d3a7c56bc162483193d4
   !=======================================================
   !FACT = 'N'
   UPLO =  'L'  ! 'U'  (Upper or lower half of A is used (L implies first index is .le. second)
   N = 3*Nsrc
   NRHS = 1  !  The number of columns in matrix B
   ! Real(dp) :: Work(  !  WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   ! Integer :: LWork(N*N)  ! If LWORK = -1, then a workspace query is assumed
   ! INFO = 0: successful exit
   !       < 0: if INFO = -i, the i-th argument had an illegal value
   !       > 0: if INFO = i, D(i,i) is exactly zero.  The factorization has been completed, but the block diagonal
   !            matrix D is exactly singular, so the solution could not be computed.
   If(RealDip) Then
      ReA(:,:)=RealPart(A(:,:))
      ReF(:,:)=RealPart(F(:,:))
      If(NF.eq.-2) Then
      LWORK=-1
         Call dsysv (UPLO, N, NRHS, ReA, N, IPIV, ReF, N, D, LWORK, INFO)
         LWORK=nint(D)
         Allocate (ReWork(LWORK))
      EndIf
      Call dsysv (UPLO, N, NRHS, ReA, N, IPIV, ReF, N, ReWORK, LWORK, INFO)
      If(INFO.gt.0) Then
         write(2,*) 'Convergence criterium:',INFO
         Do i=1,N
            Write(2,*) i, IPIV(i), ReA(i,1:i)
         EndDo
      EndIf
      F(:,1)=ReF(:,1)
   Else
      If(NF.le.-2) Then
      LWORK=-1
         Call ZHESV (UPLO, N, NRHS, A, N, IPIV, F, N, Gp, LWORK, INFO)
         LWORK=nint(RealPart(Gp))
         !write(2,*) 'ZHESV: workspace=', LWORK
         Allocate (Work(LWORK))
      EndIf
      !Do i=1,N
      !   Write(2,*) 'A:',i, A(i,1:i)
      !EndDo
      !Write(2,*) 'F:',F(1:N,1)
      Call ZHESV (UPLO, N, NRHS, A, N, IPIV, F, N, WORK, LWORK, INFO)
      !DeAllocate(Work(LWORK))
      ! on return: A contails factorization,  F the solution
      !Write(2,*) 'I:',F(1:N,1)
      !
      !===============================================
      !     - First steps (0....4) are to construct the dipole vectors for the M sources
      !        2) What to do for (near zero) eigenvalues?
      !=======================================================
      If(INFO.gt.0) Then
         write(2,*) 'Convergence criterium:',INFO
      !  If IPIV(k) > 0, then rows and columns k and IPIV(k) were interchanged, and D(k,k) is a 1-by-1 diagonal block.
      !  If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and columns k-1 and -IPIV(k) were interchanged and
      !     D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
      !  If UPLO = 'L' and IPIV(k) = IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k)
      !   were interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
         Do i=1,N
            Write(2,*) i, IPIV(i), A(i,1:i)
         EndDo
      EndIf
      !
   EndIf ! If(RealDip)
   !===============================================
   !     - Calculate R[s=t_c,ant,k] = W(ant) (E[t_c,ant,k]- \sum_m  \left( \hat{r}_{ak}\cdot \vec{I}_m )\right) G[m,ant,t_c]
   !
   !===============================================
   ! note that the solution is residing in F
   DelDip(:,:)=0.
   Do m=1,Nsrc
      DelDip(1:3,m)=F(3*m-2:3*m,1)
   EndDo
   !write(2,*) 'DelDip(1:3,m):',DelDip(1:3,1)
   !stop
   !
   ! Check result
   goto 1
   Ck(:)=0.
   Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
      W=Weight(j_IntFer)
      Do m=1,Nsrc
         E_p=0.
         E_t=0.
         Do n=1,Nsrc
            E_p=E_p+ SUM( Vec_p(:,j_IntFer)*DelDip(:,n) )* &
                  SUM( CONJG( G_IR(-Nt_l:Nt_l,m,j_IntFer) ) * ( G_IR(-Nt_l:Nt_l,n,j_IntFer) ) )
            E_t=E_t+ SUM( Vec_t(:,j_IntFer)*DelDip(:,n) )* &
                  SUM( CONJG( G_IR(-Nt_l:Nt_l,m,j_IntFer) ) * ( G_IR(-Nt_l:Nt_l,n,j_IntFer) ) )
         EndDo !  n=1,Nsrc
         Do i=1,3
            Ck(3*m+i-3)=Ck(3*m+i-3) + W*(Vec_p(i,j_IntFer) * E_p + Vec_t(i,j_IntFer)* E_t )
         EndDo ! i=1,3
      EndDo !  m=1,Nsrc
      !write(2,*) 'F(3):', j_IntFer, F(3,1), Gp, Gt
   EndDo ! j_IntFer=1,Nr_IntFer
   write(2,*) 'Ck:',Ck(:)
   write(2,*) ' B:', B(:)
1  Continue
   !
   If(MakeCurtain) then
      txt=TRIM(OutFileLabel)//TRIM(MDDLabel)
      write(FMT,"(A,I4,A)") '(I4,',Nr_IntFer,'G13.3)'
      OPEN(UNIT=31,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//'PhiDat_'//TRIM(txt)//'.dat') ! for curtain plots
      OPEN(UNIT=32,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//'ThDat_'//TRIM(txt)//'.dat')
      OPEN(UNIT=33,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//'PhiMod_'//TRIM(txt)//'.dat')
      OPEN(UNIT=34,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//'ThMod_'//TRIM(txt)//'.dat')
      Power_p(:)=0.
      Power_t(:)=0.
   EndIf    ! If(MakeCurtain)
   !
   DDChiSQ(:,:)=0.
   St_chi2(1:4)=0.
   Do I_s=-Nt_l,Nt_l
      Do j_IntFer=1,Nr_IntFer !=Nr_IntFer   ! Loop over selected antennas
         E_p=0.
         E_t=0.
         Do m=1,Nsrc
            E_p=E_p + SUM( Vec_p(:,j_IntFer)*DelDip(:,m) ) * G_IR(I_s,m,j_IntFer)
            E_t=E_t + SUM( Vec_t(:,j_IntFer)*DelDip(:,m) ) * G_IR(I_s,m,j_IntFer)
         EndDo !  m=1,Nsrc
         CE_MDD_p(j_IntFer)=SQWght(j_IntFer)*E_p
         CE_MDD_t(j_IntFer)=SQWght(j_IntFer)*E_t
         j= Nt_l+I_s+ (j_IntFer-1)*(2*T_Range)
         !R(2*j+1)           =RealPart( (CTime_p(MDDnuDim+I_s,j_IntFer) - SQWght(j_IntFer)*E_p) )
         !R(2*(j+ T_Range)+1)=RealPart( (CTime_t(MDDnuDim+I_s,j_IntFer) - SQWght(j_IntFer)*E_t) )
         !R(2*j+2)           =ImagPart( (CTime_p(MDDnuDim+I_s,j_IntFer) - SQWght(j_IntFer)*E_p) )
         !R(2*(j+ T_Range)+2)=ImagPart( (CTime_t(MDDnuDim+I_s,j_IntFer) - SQWght(j_IntFer)*E_t) )
         Call CalcStokes(CE_MDD_p(j_IntFer),CE_MDD_t(j_IntFer), St_MDD)
         Call CalcStokes(CTime_p(MDDnuDim+I_s,j_IntFer),CTime_t(MDDnuDim+I_s,j_IntFer), St_dat)
         !
         R(2*j+1)           =SQWght(j_IntFer)*( St_dat(1)-St_MDD(1) )
         R(2*(j+ T_Range)+1)=SQWght(j_IntFer)*( St_dat(2)-St_MDD(2) )
         R(2*j+2)           =SQWght(j_IntFer)*( St_dat(3)-St_MDD(3) )
         R(2*(j+ T_Range)+2)=SQWght(j_IntFer)*( St_dat(4)-St_MDD(4) )
      If(J_intFer.eq.-1 .and. I_s.gt.-20 .and. I_s.lt.-10) Then
            write(2,*) 'I_s,St_MDD(1:4)', I_s,St_MDD(1),St_dat(1),St_MDD(2:4)/St_MDD(1),St_dat(2:4)/St_dat(1)
            Write(2,*) '  E_p (mod..dat)', SQWght(j_IntFer)*E_p,CTime_p(MDDnuDim+I_s,j_IntFer)
            Write(2,*) '  E_t (mod..dat)', SQWght(j_IntFer)*E_t,CTime_t(MDDnuDim+I_s,j_IntFer)
      EndIf
         DDChiSQ(j_IntFer,1)=DDChiSQ(j_IntFer,1) + R(2*j+1)**2 + R(2*j+2)**2
         DDChiSQ(j_IntFer,2)=DDChiSQ(j_IntFer,2) + R(2*(j+ T_Range)+1)**2  + R(2*(j+ T_Range)+2)**2
         St_chi2(1:4)=St_chi2(1:4)+ (SQWght(j_IntFer)*( St_dat(1:4)-St_MDD(1:4) ))**2
      EndDo ! j_IntFer=1,Nr_IntFer=Nr_IntFerCh(1)
      If(MakeCurtain) then
         Power_p(1:Nr_IntFer)=Power_p(1:Nr_IntFer) + ABS( CTime_p(MDDnuDim+I_s,1:Nr_IntFer) )**2
         Power_t(1:Nr_IntFer)=Power_t(1:Nr_IntFer) + ABS( CTime_t(MDDnuDim+I_s,1:Nr_IntFer) )**2
         write(31,FMT)  I_s, RealPart( CTime_p(MDDnuDim+I_s,1:Nr_IntFer) ) ! Data
         write(32,FMT)  I_s, RealPart( CTime_t(MDDnuDim+I_s,1:Nr_IntFer) )
         write(33,FMT)  I_s, RealPart( CE_MDD_p(1:Nr_IntFer) ) ! model
         write(34,FMT)  I_s, RealPart( CE_MDD_t(1:Nr_IntFer) )
      EndIf
   EndDo !  I_s=-Nt_l,Nt_l
   !stop
   If(NF.ge. (NF_old+20)) Then
      W=SUM(Weight(1:Nr_IntFer))  ! *T_Range
      CH1=SUM(DDChiSQ(1:Nr_IntFer,1))/W
      CH2=SUM(DDChiSQ(1:Nr_IntFer,2))/W
      !write(2,*) nf,' mean chi^2=',CH1, CH2, ' previous:', CH1_old, CH2_old
      write(2,*) 'Steepest descent iteration#',nf,', Bias_inv=', Bias_inv
      !write(2,*) 'Weight(1:10)', Weight(1:10)
      !write(2,*) 'DelDip(:,1:3)', DelDip(:,1:3)
      !   Do i=1,N
      !      Write(2,*) i, IPIV(i), ReA(i,1:i)
       !  EndDo
      Call MDD_write()
      NF_old=NF
      CH1_old=CH1
      CH2_old=CH2
   EndIf
   !
   !Write(2,*) 'St_chi2(1:4):', St_chi2(1:4)
   If(MakeCurtain) then
      !Power_p(:)=Power_p(:)/(2*N_s+1.)
      !Power_t(:)=Power_t(:)/(2*N_s+1.)
      Close(UNIT=31)
      Close(UNIT=32)
      Close(UNIT=33)
      Close(UNIT=34)
      CH1=SUM(DDChiSQ(1:Nr_IntFer,1:2))/(2*SUM(Weight(1:Nr_IntFer)))
      Write(GLE_file,"('CuP-',A)") trim(txt)
!      Call GLEscript_Curt_MDD(30, GLE_file, Nt_l, i_chunk, trim(DataFolder)//TRIM(OutFileLabel), TRIM(txt), &
!               DDChiSQ(1,1), DDChiSQ(1,2), Power_p, Power_t, Ch1, dCenloc(:,1) )
      Call GLEscript_Curtains(30, GLE_file, Nt_l, i_chunk, trim(DataFolder), &
               TRIM(txt), DDChiSQ(1,1), DDChiSQ(1,2), Power_p, Power_t, Ch1, dCenloc(:,1) )
      write(2,*) 'Make CurtainPlot ',trim(txt)
      NF_old=0
      CH1_old=-1.d0
      CH2_old=-1.d0
      Call MDD_SrcsFile(txt)
   EndIf
   !
   !Call GLEplotControl( Submit=.true.)  ;
   !stop 'CompareMDD'
   Return
End Subroutine CompareMDD
!-----------------------------
!INTERFACE
pure Subroutine CalcStokes(E_p,E_t,St)
   use constants, only : dp
   Implicit none
   Complex(dp), intent(in) :: E_p,E_t
   real(dp), intent(out)  :: St(4)
   St(1)=E_p*CONJG(E_p)+E_t*CONJG(E_t)   ! StI
   St(2)=E_p*CONJG(E_p)-E_t*CONJG(E_t)   ! StQ
   St(3)=2*REALPART( E_p*CONJG(E_t) )   ! StU
   St(4)=2*IMAGPART( E_p*CONJG(E_t) )   ! StV
   Return
end Subroutine CalcStokes
!End INTERFACE
!==================================================
Subroutine MDD_QuickFit(FitQual,error)
!
!
   use constants, only : dp
   Use MDD_Pars, only : NSrc, dCenT, dCenLoc, T_Range, DDChiSQ
   Use MDD_Pars, only : WORK, ReWORK, LWORK, G_IR, Xscale, RealDip   ! Used in CompareMDD
   Use Interferom_Pars, only : IntFer_ant, Nr_IntFerCh, Nr_IntFerMx
   use DataConstants, only : RunMode
   Implicit none
   Integer, intent(out) :: error
   Real(dp), Intent(OUT) :: FitQual
   integer ( kind = 4 ) :: i, j, k
   integer ( kind = 4 ) :: nvar    ! number of parameters
   integer ( kind = 4 ) :: meqn  ! Number of data points to fit
   integer ( kind = 4 ) :: uiparm(1)    ! Not really used
   real ( kind = 8 ) :: urparm(1)    ! Not really used
   external ufparm  ! Dummy external routine
   integer ( kind = 4 ) :: NF
   real ( kind = 8 ),allocatable :: R(:) !
   real ( kind = 8 ), allocatable :: X(:)
   external CompareMDD  ! subroutine that compares FitFunction to Data
   Integer :: m
   !Real(dp) :: Chi2_start !, Chi2StopRatio=1.0
   !
   Meqn=T_Range*Nr_IntFerCh(1)*2*2   ! total number of data-points that enter in the chi^2 search one *2 for 2 pol directions, one *2 for real/imag parts
   nvar=4*NSrc
   !write(2,*) 'Fitting statistics: nvar=',nvar,', Meqn=',Meqn
   flush(unit=2)
   allocate( R(Meqn) )
   allocate( X(nvar) )
   Allocate( G_IR(-T_Range/2:T_Range/2,1:NSrc,1:Nr_IntFerCh(1)) )
   !
   Do m=1,Nsrc
      X(4*m-3)=dCenT(m)/Xscale
      X(4*m-2:4*m)=dCenloc(1:3,m)/Xscale
   EndDo
   !
   NF=-2 ! make plot
   call CompareMDD ( meqn, nvar, x, NF, R, uiparm, urparm, ufparm )  ! prepare allocation workspace
   !
   If(RealDip) Then
      Deallocate( R, X, ReWork, G_IR )
   Else
      Deallocate( R, X, Work, G_IR )
   EndIf
   LWORK=-1
   !
   Return
End Subroutine MDD_QuickFit
!-----------------------------------------------
