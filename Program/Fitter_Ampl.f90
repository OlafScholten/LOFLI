    !include 'C:/OlafsUtil/Minimizer_ASA/asa047.f90'
    !Include '../../NumLib/LSQ/nl2sol.f90'
    Include 'nl2sol.f90'
!===============================================
Subroutine FitAmpl(Meqn, nvar, X,ChiSQDF, BestFit)
   Use Tracks,  only : i_AmplTh, Nmax_Ampl, AmplScale, Ampl_Hist, FitFunc
   ! i_AmplTh=first entry that will be included in the fit
    Implicit none
    Integer( kind = 4 ), intent(inout) :: Meqn
    integer ( kind = 4 ),intent(in) :: nvar    ! number of parameters
    real ( kind = 8 ), intent(inout) :: X(4)
    real ( kind = 8 ), intent(out) :: ChiSQDF
    integer ( kind = 4 ) :: i, j, k
    integer :: v_dim, error
    !integer ( kind = 4 ) :: meqn  ! Number of data points to fit
    integer ( kind = 4 ) iv(60+nvar)
    external CompareAmpl  ! subroutine that compares FitFunction to Data
    External JacobianAmpl
    external ufparm  ! Dummy external routine
    integer ( kind = 4 ) uiparm(1)    ! Not really used
    real ( kind = 8 ) :: Jacobian(Meqn,nvar), BestFit(Meqn)
    !real ( kind = 8 ) v(v_dim) ! dim=93 + n*p + 3*n + p*(3*p+33)/2 ,n-meqn, p=nvar
    real ( kind = 8 ),allocatable :: v(:) !  in NL@SOL: real ( kind = 8 ) v(93 + n*p + 3*n + (p*(3*p+33))/2)
    integer ( kind = 4 ) :: NF
    !
    !     write(2,*) 'entering fitting', nvar, Meqn
    !      flush(unit=2)
    If(nvar.ge.Meqn) then
      Write(2,*) '****** too few values to fit, ',Meqn, FitFunc
      Meqn=-1
      !    flush(unit=2)
      Return
    EndIf
    v_dim=93 + Meqn*nvar + 3*Meqn + nvar*(3*nvar+33)/2
    !allocate( Jacobian(Meqn,nvar) )
    allocate( v(v_dim) )
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
    !If(CalcHessian) then
    !  iv(14) = 1
    !  iv(15) = 2
    !  v(32) =1.d-8 ! is the relative function convergence tolerance
    !endif
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
   !call nl2sno ( meqn, nvar, x, CompareAmpl, iv, v, uiparm, Jacobian, ufparm )
   flush(unit=2)
   Call nl2sol ( meqn, nvar, x, CompareAmpl, JacobianAmpl, iv, v, uiparm, Jacobian, ufparm, error )
   If(error.gt.0) then
      write(2,*) 'endless loop in NL2SOL truncated'
      goto 9
   endif
   !
   ChiSQDF=2*v(10)/(meqn-nvar)
   NF=-1  ! Store fitted function in R
   call CompareAmpl ( meqn, nvar, X, NF, BestFit, uiparm, Jacobian, ufparm )
9  Continue
    !
    !Deallocate( Jacobian)
    Deallocate( v )
    !
    ! write(2,"(A,I2,A,I2,A,F7.2,A,4G12.5)") 'Final result fit function=', FitFunc, &
    !  ', First amplitude=', i_AmplTh,', Chi^2=', ChiSQDF,', Parameters=', X(1:nvar)
    Return
    !
End Subroutine FitAmpl
!=================================================
Subroutine JacobianAmpl ( meqn, nvar, X_p, nf, Jacobian, uiparm, urparm, ufparm )
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
  Call CompareAmpl ( meqn, nvar, X_p, nf, R, uiparm, Jacob, ufparm )
  !write(*,*) 'JacobN',meqn,Jacob(meqn,:)
  Jacobian(:,:)=Jacob(:,:)
  Return
  End Subroutine JacobianAmpl
!=================================================
Subroutine CompareAmpl ( meqn, nvar, X_p, nf, R, uiparm, Jacobian, ufparm )
!*****************************************************************************80
!
!  Discussion:
!    Given the value of the vector X, this routine computes the value of the residues R=(t_measured-t_X) used for chi^2 fitting
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
!    Input, real ( kind = 8 ) URPARM(*)==Jacobian, a user array.
!    Input, external UFPARM, an external reference to a user subroutine
!       or function.
!
   Use Tracks,  only : i_AmplTh, Ampl, Ampl_Hist, Ampl_Hist_W, FitFunc
!
  implicit none
    real ( kind = 8 ), intent(in) :: X_p(1:nvar)
  integer ( kind = 4 ), intent(in) :: meqn
  integer ( kind = 4 ), intent(in) :: nvar
  integer ( kind = 4 ), intent(in) :: nf
  real ( kind = 8 ), intent(out) :: R(meqn)
  external ufparm
  integer ( kind = 4 ), intent(in) :: uiparm(*)
  real ( kind = 8 ), intent(out) :: Jacobian(meqn,nvar)  ! derivative of R w.r.t. parameters
  real ( kind = 8 ) :: X(nvar), A, CalcN, Weight, chisq
  real ( kind = 8 ), parameter :: ZeroErr=0.5
   !
   integer :: i,j,k, i_eqn
   !
   ChiSq=0.
   Select Case(FitFunc)
      Case (5)   ! F(A)=b A^(-a)*exp(-c/A) = exp(-a ln(A) +b' -c/A)
         Do i_eqn=1,meqn
            k=i_AmplTh-1+i_eqn
            Weight = Ampl_Hist_W(k)
            A= Ampl(k)
            CalcN = A**(-X_p(1)) * exp(-X_p(3)/A)
            R(i_eqn)= (Ampl_Hist(k) -X_p(2)* CalcN)*Weight
            ChiSq = ChiSq + R(i_eqn)*R(i_eqn)
            !
            Jacobian(i_eqn,1) =  +X_p(2)*CalcN *Weight*log(A)
            Jacobian(i_eqn,2) =  - CalcN *Weight
            Jacobian(i_eqn,3) =  + X_p(2)*CalcN *Weight/A
            !If(nf.lt.0) write(2,*) i_eqn,k, A, Ampl_Hist(k), X_p(2)* CalcN, R(i_eqn)*Weight, Weight
            If(nf.lt.0) R(i_eqn)=X_p(2)* CalcN
         Enddo  ! i_eqn=1,meqn
      Case (4)   ! F(A)=b (A+1)^(-a)*exp(-c/A^2)
         Do i_eqn=1,meqn
            k = i_AmplTh-1+i_eqn
            Weight = Ampl_Hist_W(k)
            A = Ampl(k)
            CalcN = A**(-X_p(1)) * exp(-X_p(3)/(A*A))
            R(i_eqn)= (Ampl_Hist(k) -X_p(2)* CalcN)*Weight
            ChiSq = ChiSq + R(i_eqn)*R(i_eqn)
            !
            Jacobian(i_eqn,1) =  +X_p(2)*CalcN *Weight*log(A)
            Jacobian(i_eqn,2) =  - CalcN *Weight
            Jacobian(i_eqn,3) =  + X_p(2)*CalcN *Weight/(A*A)
            If(nf.lt.0) R(i_eqn)=X_p(2)* CalcN
         Enddo  ! i_eqn=1,meqn
      Case (3)     !  exp(-x(1)*A+b-x(3)/A^1.5)
         Do i_eqn=1,meqn
            k = i_AmplTh-1+i_eqn
            Weight = Ampl_Hist_W(k)
            A = Ampl(k)
            CalcN = exp(-X_p(1)*A+X_p(2)-X_p(3)/(A*A))
            R(i_eqn)= (Ampl_Hist(k) -CalcN)*Weight
            !ChiSq = ChiSq + R(i_eqn)*R(i_eqn)
            !
            Jacobian(i_eqn,1) =  +A*CalcN *Weight  ! =dR/da
            Jacobian(i_eqn,2) =  - CalcN *Weight      ! =dR/db
            Jacobian(i_eqn,3) =  + CalcN *Weight/(A*A)      ! =dR/dc
            !If(nf.lt.0) write(2,*) i_eqn,k, A, Ampl_Hist(k), CalcN, R(i_eqn)*Weight, Weight
            If(nf.lt.0) R(i_eqn)= CalcN
         Enddo  ! i_eqn=1,meqn
      Case (2)
         Do i_eqn=1,meqn !  exp(-x(1)*A+b-x(3)/A)
            k = i_AmplTh-1+i_eqn
            Weight = Ampl_Hist_W(k)
            A = Ampl(k)
            CalcN = exp(-X_p(1)*A+X_p(2)-X_p(3)/A)
            R(i_eqn)= (Ampl_Hist(k) -CalcN)*Weight
            !
            Jacobian(i_eqn,1) =  +A*CalcN *Weight  ! =dR/da
            Jacobian(i_eqn,2) =  - CalcN *Weight      ! =dR/db
            Jacobian(i_eqn,3) =  + CalcN *Weight/A      ! =dR/dc
            If(nf.lt.0) R(i_eqn)= CalcN
         Enddo  ! i_eqn=1,meqn
      Case (1)  ! powerlaw
         Do i_eqn=1,meqn
            k = i_AmplTh-1+i_eqn
            Weight = Ampl_Hist_W(k)
            A = Ampl(k)
            CalcN = A**(-X_p(1))
            R(i_eqn)= (Ampl_Hist(k) -X_p(2)* CalcN)*Weight
            ChiSq = ChiSq + R(i_eqn)*R(i_eqn)
            !
            Jacobian(i_eqn,1) =  +X_p(2)*CalcN *Weight*log(A)
            Jacobian(i_eqn,2) =  - CalcN *Weight
            If(nf.lt.0) R(i_eqn)=X_p(2)* CalcN
         Enddo  ! i_eqn=1,meqn
      Case Default          !  exp(-x(1)*A+b)
         Do i_eqn=1,meqn !  exp
            k = i_AmplTh-1+i_eqn
            Weight = Ampl_Hist_W(k)
            A = Ampl(k)
            CalcN = exp(-X_p(1)*A+X_p(2))
            R(i_eqn)= (Ampl_Hist(k) -CalcN)*Weight
            !ChiSq(i_stat,i_Peak) = ChiSq(i_stat,i_Peak) + R(i_eqn)*R(i_eqn)
            !
            Jacobian(i_eqn,1) =  +A*CalcN *Weight
            Jacobian(i_eqn,2) =  - CalcN *Weight
            !If(nf.lt.0) write(2,*) i_eqn,k, A, Ampl_Hist(k), CalcN, R(i_eqn)*Weight, Weight
            If(nf.lt.0) R(i_eqn)= CalcN
         Enddo  ! i_eqn=1,meqn
   End Select
   !write(2,*) 'ChiSq:',ChiSq/(meqn-nvar),nf,X_P(1:nvar)
   !stop
   return
End Subroutine CompareAmpl
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
