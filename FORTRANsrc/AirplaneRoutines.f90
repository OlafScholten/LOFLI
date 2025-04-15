!=================================
Subroutine PolPCACathCon(Stk_NEh_Con, PolZen, PolAzi, PolMag, PoldOm, prin)
   !Perform PolPCACath on the stokes matrix given in condensed form
   !NOTE: The Stokes matrix is single precision complex
   use constants, only : dp,pi
   Implicit none
   Complex, intent(in) :: Stk_NEh_Con(1:6)
   Real(dp), intent(out) :: PolZen(1:3), PolAzi(1:3), PolMag(1:3), PoldOm(1:3)
   Logical, intent(in) :: prin
   complex  :: Stk_NEh(1:3,1:3) !  ,TrackLenMax)
   !
   Stk_NEh(1,1:3)=Stk_NEh_Con(1:3)
   Stk_NEh(2,2:3)=Stk_NEh_Con(4:5)
   Stk_NEh(3,3)=Stk_NEh_Con(6)
   Stk_NEh(2:3,1)=conjg(Stk_NEh(1,2:3))
   Stk_NEh(3,2)=conjg(Stk_NEh(2,3))
   Call PolPCACath(Stk_NEh, PolZen, PolAzi, PolMag, PoldOm, prin)
   Return
End Subroutine PolPCACathCon
!-----------------------------------------------
Subroutine PolPCACath(Stk_NEh, PolZen, PolAzi, PolMag, PoldOm, prin)
   ! Here we use the full complex Stokes matrix, however one expects that for
   !  averaging over long slices, or over many slices (the same thing) the
   !  Real part of the Stokes matrix will dominate because it is assumed that
   !  the imaginary part will average to zero for a large enough ensemble
   !  i.e. assumes there is no preferred circular polarization direction.
   ! In Pricipal Component Analysis one maximizes
   !  sum_i( p_i . P)^2  for (th) and (ph) where
   !  P=[cos(ph) cos(yh) , sin(ph) cos(th) , sin(th)]
   !  and where Stk_NEh(x,y)= sum_i( p^x_i p^y_i)
   ! Note that performing a  Pricipal Component Analysis on the polarization vector is
   !  identical to finding the eigenvector for the largest eigenvalue of the Stokes matrix.
   use constants, only : dp,pi
   Implicit none
   Complex, intent(in) :: Stk_NEh(1:3,1:3)
   Real(dp), intent(out) :: PolZen(1:3), PolAzi(1:3), PolMag(1:3), PoldOm(1:3)
   Logical, intent(in) :: prin
!   Real(dp), intent(in) :: VoxLoc(1:3)
   Real(dp) :: A, theta, phi, N, E, h, dOmg, cPolVec(1:3), I_circ, I_lin
   Complex :: C, EigVec(3,3)
   Integer :: i ,j,k,m
   Character(len=1) :: JOBZ= 'V'   !  Compute eigenvalues and eigenvectors.,
   Character(len=1) :: UPLO= 'L'   !  Lower triangle of A is stored., a(1,1),(2,1),(3,1)
   Integer :: Adim=3   ! = 'L':  Lower triangle of A is stored.,
   Complex :: Amat(3,3), Rmat(3,3)
   Integer :: LDA=3
   Real :: WW(3) ! is REAL array, dimension (N),          If INFO = 0, the eigenvalues in ascending order.,
   Complex :: WORK(99)  ! is COMPLEX array, dimension (MAX(1,LWORK))  On exit, if INFO = 0, WORK(1) returns the optimal
   Integer :: LWORK=99 !The length of the array WORK.  LWORK >= max(1,2*N-1).
   !       For optimal efficiency, LWORK >= (NB+1)*N,   where NB is the blocksize for CHETRD returned by ILAENV.
   !       If LWORK = -1, then a workspace query is assumed;.,
   Real :: RWORK(7) ! is REAL array, dimension (max(1, 3*N-2)),
   Integer :: INFO
   !
   !  CHEEV computes all eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix A.
   Amat(:,:)=Stk_NEh(:,:) !  A is COMPLEX array, dimension (LDA, N)
   !       On entry, the Hermitian matrix A.  If UPLO = 'U', the          leading N-by-N upper triangular part of A contains the
   !       upper triangular part of the matrix A.  If UPLO = 'L',          the leading N-by-N lower triangular part of A contains
   !       the lower triangular part of the matrix A.
   !       On exit, if JOBZ = 'V', then if INFO = 0, A contains the          orthonormal eigenvectors of the matrix A.
   Call cheev (JOBZ, UPLO, Adim, Amat, LDA, WW, WORK, LWORK, RWORK, INFO) ! from LAPACK library
   !write(2,*) 'Work(1)',Work(1), WW
   !Amat(:,:)=NEh2PB(:,:)
   !Call RotTensor(Stk_NEh, Amat, Stk_ev) ! just a test
   I_circ=0.
   I_lin=0.
   Do i=1,3
      !j=MAXLOC(ABS(Amat(:,i)),1)
      !C=CONJG(Amat(j,i))/ABS(Amat(j,i))  ! make the largest component real
      !EigVec(:,i)=Amat(:,i)*C
      !! get angle mean polarization direction
      !N=Real(EigVec(1,i))  ; E=Real(EigVec(2,i))  ; h=Real(EigVec(3,i))
      !A=N*N + E*E
      !Theta=atan2(sqrt(A),h)*180./pi  ; Phi=atan2(E,N)*180./pi
      !dOmg=acos(sqrt(A+h*h))*180./pi
      !write(2,*) 'Eig Vec(i)', WW(i), ' ; ', EigVec(:,i), ' ; th,ph,dO=', theta, Phi, dOmg
      !
      ! Different approach, maximize the real component of EigVec(:,i) by zeroing Im(EigVec(:,i)*EigVec(:,i)*e^i\phi)
      ! Find the phase that maximizes the length of the real part
      C=sqrt(SUM(Amat(:,i)*Amat(:,i)))  ! not conjugated here!!!!
      C=CONJG(C)/ABS(C)
      EigVec(:,i)=Amat(:,i)*C
      ! get angle mean polarization direction
      N=Real(EigVec(1,i))  ; E=Real(EigVec(2,i))  ; h=Real(EigVec(3,i))
      Call Carth2Angle(N,E,h,Theta, Phi, dOmg)
      PolZen(4-i)=Theta  ; PolAzi(4-i)=Phi  ; PolMag(4-i)=WW(i)  ; PoldOm(4-i)=dOmg
      !
      If(prin) Then
         A=cos(2*dOmg*pi/180.)
         I_circ=I_circ+(WW(i)*sin(2*dOmg*pi/180.))**2
         I_lin=I_lin + (WW(i)*A)**2
         write(2,"(A,F10.2,A,F10.2, A, 3F7.1)") 'Vec: I_lin=', abs(WW(i)*A), ' ; I_circ=', &
            abs(WW(i)*sin(2*dOmg*pi/180.)), ' ; th(3),ph(1),dO=', theta, Phi, dOmg
      EndIf
   EndDo !i=1,3
   If(prin) Then


!   StI = Real(Stk(1,1) + Stk(2,2) + Stk(3,3))  ! =\Lambda_0  from  PHYSICAL REVIEW E 66, 016615 ~2002!
!   StI12=Real(Stk(1,1) + Stk(2,2))             ! =
!   StI3 =  Real(Stk(3,3))                      ! =>\Lambda_8 = StI12 *sq(3)/2 - StI3 *sq(3)
!   StQ = Real(Stk(1,1)-Stk(2,2))               ! =\Lambda_3 *2/3
!   StU = 2*Real(Stk(1,2))                      ! =\Lambda_1 *2/3
!   StV =  2*Imag(Stk(1,2))                     ! =\Lambda_2 *2/3
!   StU2 = 2*Real(Stk(3,1))                     ! =\Lambda_4 *2/3
!   StV2 =  2*Imag(Stk(3,1))                    ! =\Lambda_5 *2/3
!   StU1 = 2*Real(Stk(2,3))                     ! =\Lambda_6 *2/3
!   StV1 =  2*Imag(Stk(2,3))                    ! =\Lambda_7 *2/3
!   St8 = (StI12  - 2* StI3)/sqrt(3.)           ! =\Lambda_8 *2/3
!   P_un = 1.- (3./4.)*(StQ**2+StU**2+StV**2+StU1**2+StV1**2+StU2**2+StV2**2+St8**2)/(StI**2)
!   P_lin = (3./4.)*(StQ**2+StU**2+StU1**2+StU2**2+St8**2)/(StI**2)
!      write(2,*) Stk_NEh(1,1), Stk_NEh(1,2), Stk_NEh(1,3), Stk_NEh(2,2), Stk_NEh(2,3), Stk_NEh(3,3)


      I_Circ=3.*( (IMAG(Stk_NEh(1,2)))**2 +(IMAG(Stk_NEh(1,3)))**2 + (IMAG(Stk_NEh(2,3)))**2 )
      I_Lin =3.*( (REAL(Stk_NEh(1,2)))**2 +(REAL(Stk_NEh(1,3)))**2 + (REAL(Stk_NEh(2,3)))**2 ) &
         + (3.*(Real(Stk_NEh(1,1)-Stk_NEh(2,2)) )**2 + (Real(Stk_NEh(1,1) + Stk_NEh(2,2) -2*Stk_NEh(3,3)))**2)/4.
      A=( Real( Stk_NEh(1,1) + Stk_NEh(2,2) + Stk_NEh(3,3) ) )**2
      !SourceUn(i_Peak)=1.-SourceLin(i_Peak)-SourceCirc(i_Peak)
      write(2,"(A,G12.4,2(A,F6.2))") 'total intensity_PCA=',sqrt(A), ', % lin=', 100.*I_lin/A, ', % circ=', 100.*I_circ/A
   EndIf
   !
   Return
End Subroutine PolPCACath
!--------------------------------
Subroutine Carth2Angle(N,E,h,Theta, Phi, dOmg)
   use constants, only : dp,pi
   Implicit none
   Real(dp), intent(in) :: N, E, h
   Real(dp), intent(out) :: Theta, Phi, dOmg
   Real(dp) :: A
   A=N*N + E*E
   If(h.lt.0) Then
      Theta=atan2(sqrt(A),-h)*180./pi  ;
      Phi=atan2(-E,-N)*180./pi ! The real part corresponds to lin polarized
   Else
      Theta=atan2(sqrt(A),h)*180./pi  ;
      Phi=atan2(E,N)*180./pi ! The real part corresponds to lin polarized
   EndIf
   dOmg=acos(sqrt(A+h*h))*180./pi  ! the angle between the real and the complex vector
End Subroutine Carth2Angle
! ----------------------------------
! ----------------------------------
Subroutine RRotCTensor(Stk_NEh, NEh2PB, Stk_PB)
   ! Stk_NEh=NEh2PB^T Stk_PB NEh2PB
   !First index in NEh2PB is the old (NEh) basis, second the new PB one
   use constants, only : dp
   Implicit none
   Complex, intent(in) :: Stk_PB(3,3)
   Real(dp), intent(in) :: NEh2PB(3,3)
   Complex, intent(out) :: Stk_NEh(3,3)
   Complex :: C
   Integer :: i,j,k !,m
   Stk_NEh(:,:)=0
   Do i=1,3
      Do j=1,3
         Do k=1,3
            Stk_NEh(i,j)=Stk_NEh(i,j)+ Sum(NEh2PB(i,:)*Stk_PB(:,k))*NEh2PB(j,k)
         EndDo !k=1,3
      EndDo !j=1,3
   EndDo !i=1,3
End Subroutine RRotCTensor
!-----------------------------------------------
!===============================================
    Subroutine GetNonZeroLine(lineTXT)
    implicit none
    character(LEN=*), intent(inout) :: lineTXT
    Character(len=6) :: FMT
    integer :: eof,N
      N=LEN(lineTXT)
      write(FMT,"(A,I3.3,A)") '(A',N,')'
      !write(2,*) 'GetNonZeroLine2:', FMT, N
      lineTXT=''
      Do while (trim(lineTXT).eq.'')
         read(*,FMT,iostat=eof) lineTXT
         if(eof.lt.0) exit
         If(lineTXT(1:1).eq. '!') Then
            write(2,"(A)") lineTXT
            lineTXT=''
         EndIf
      enddo
    End Subroutine GetNonZeroLine
!===============================================
    Subroutine GetMarkedLine(Mark,lineTXT)
    implicit none
    character(LEN=*), intent(inout) :: lineTXT
    character(LEN=*), intent(out) :: Mark
    Character(len=6) :: FMT
    integer :: eof,N,j
      N=LEN(lineTXT)
      write(FMT,"(A,I3.3,A)") '(A',N,')'
      !write(2,*) 'GetNonZeroLine2:', FMT, N
   1  continue
      lineTXT=''
      Do while (trim(lineTXT).eq.'')
         read(*,FMT,iostat=eof) lineTXT
         if(eof.lt.0) Return
      enddo
      If(lineTXT(1:1).eq.'!') goto 1   ! skip comment lines
      j = iachar(lineTXT(1:1))
      if (j>= iachar("a") .and. j<=iachar("z") ) then
         j = (j-32)   ! convert to upper case if not already
      end if
      Mark=achar(j)
      lineTXT(1:N-1)=lineTXT(2:N)
    End Subroutine GetMarkedLine
!==========================================
!==========================================
Subroutine Convert2m(CenLoc)
   Implicit none
   real*8, intent(inout) :: CenLoc(3)
   If((abs(CenLoc(1)) .lt. 180.) .and. (abs(CenLoc(2)) .lt. 180.) .and. (abs(CenLoc(3)) .lt. 20.) ) Then
      CenLoc(:)=CenLoc(:)*1000.  ! convert from [km] to [m]
   Endif
   Return
End Subroutine Convert2m
! ==========================
Subroutine Inverse33(A,B, eig, PolBasis)
   ! B= A^-1, from https://en.wikipedia.org/wiki/Eigenvalue_algorithm
   ! derivation, see https://www.geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf
   use constants, only : dp, pi
   Implicit none
   Real(dp), intent(in) :: A(3,3)
   Real(dp), intent(out) :: B(3,3),eig(3), PolBasis(3,3)
   integer, parameter :: c(3)=(/ 2,3,1 /) , r(3)=(/ 3,1,2 /)
   Real(dp) :: Det
   integer :: i,j,k
   Real(dp) :: p,q,s, phi
   s = A(1,2)**2 + A(1,3)**2 + A(2,3)**2 !Eigen values real symmetric 3x3 matrix
   if (s .eq. 0) then !      % A is diagonal.
      eig(1) = A(1,1)
      eig(2) = A(2,2)
      eig(3) = A(3,3)
   else
      q = (A(1,1)+A(2,2)+A(3,3))/3      !         % trace(A) is the sum of all diagonal values
      p = sqrt(((A(1,1) - q)**2 + (A(2,2) - q)**2 + (A(3,3) - q)**2 + 2 * s)/6.)
      s=(A(1,1)-q)*(A(2,2)-q)*(A(3,3)-q) + 2.*A(1,2)*A(2,3)*A(1,3) - &
         (A(1,1)-q)*A(2,3)*A(2,3) - (A(2,2)-q)*A(1,3)*A(1,3) - (A(3,3)-q)*A(1,2)*A(1,2)
      s=s/(2.*p*p*p)
      if (s .le. -1) then
         phi = pi / 3
      elseif (s .ge. 1) then
         phi = 0
      else
         phi = acos(s) / 3
      endif
      !write(2,*) 'phi',phi,s
      !% the eigenvalues satisfy eig3 <= eig2 <= eig1
      eig(1) = q + 2 * p * cos(phi)
      eig(3) = q + 2 * p * cos(phi + (2*pi/3))
      eig(2) = 3 * q - eig(1) - eig(3)   !  % since trace(A) = eig1 + eig2 + eig3
   endif
   !Write(2,*) 'eigenvalues', eig(:)
   !
   Call EigvecS33(A,eig,PolBasis)
   !
   Det=0.
   B(:,:)=0.
   Do i=1,3 ! get the Adjoint of A (transpose of subdeterminant matrices)
   Do j=1,3
   B(j,i)=A(c(i),c(j))*A(r(i),r(j))-A(c(i),r(j))*A(r(i),c(j))
   Enddo
   Enddo
   Do i=1,3
   Det=Det + A(i,1)*B(i,1)
   Enddo
   B(:,:)=B(:,:)/Det
   ! check
   !write(2,*) 'check inverse:'
   !Do i=1,3
   !Do j=1,3
   !   s=0
   !   Do k=1,3
   !      s=s+B(i,k)*A(k,j)
   !      !write(2,*) k,s,B(i,k),A(k,j)
   !   enddo
   !   write(2,*) i,j,B(i,j), A(i,j),s
   !Enddo
   !Enddo
   !
   Return
End Subroutine Inverse33
Subroutine EigvecS33(A,E,ev)
!  construct the eigen vectors ev(:,k) given the eigenvalues E(k) of the real symmetric matrix A.
!  Explicitly orthogonnormalize to the eigenvector for the smallest eigenvalue to minimize num. errors,
!  which could be necessary since EigVal could be like (1,1,0), i.e. the first two almost degenerate.
!  eigenvalues are assumed to be ordered, smallest last.
   use constants, only : dp
   Implicit none
   Real(dp), intent(in) :: A(3,3),E(3)
   Real(dp), intent(out) :: ev(3,3)
   Real(dp) :: s
   integer :: k
  ! check
   Do k=1,3
      ev(1,k)=1.
      ev(2,k)=(E(k)*A(2,3)-A(1,1)*A(2,3)+A(1,2)*A(1,3)) / (A(1,2)*A(2,3)+(E(k)-A(2,2))*A(1,3))
      ev(3,k)=(E(k)*A(2,3)-A(1,1)*A(2,3)+A(1,2)*A(1,3)) / (A(1,3)*A(2,3)+(E(k)-A(3,3))*A(1,2))
      s=SUM(ev(:,k)*ev(:,k))
      ev(:,k)=ev(:,k)/sqrt(s)
   Enddo
   s=SUM(ev(:,1)*ev(:,3))
   ev(:,1)=(ev(:,1)-s*ev(:,3))/sqrt(1-s*s)  ! ev 1 is orthonormalized to 3
   ev(1,2)=ev(2,1)*ev(3,3)-ev(3,1)*ev(2,3)  !  2= 1 x 3
   ev(2,2)=ev(3,1)*ev(1,3)-ev(1,1)*ev(3,3)  !  2= 1 x 3
   ev(3,2)=ev(1,1)*ev(2,3)-ev(2,1)*ev(1,3)  !  2= 1 x 3
   !write(2,*) 'check eigenvexcors (=0)', (SUM(ev(:,k)*ev(:,k))-1, k=1,3)
   !
   Return
End Subroutine EigvecS33
!--------------------
