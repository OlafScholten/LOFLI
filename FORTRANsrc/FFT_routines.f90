      Module FFT
    use constants, only : dp,pi,ci
    implicit none
    integer ( kind = 4 ), save :: lensav
    integer ( kind = 4 ), save :: lenwrk
    integer ( kind = 4 ), save :: N_Time
    integer ( kind = 4 ), save :: N_nu
    integer, save :: HW_size=32
    !  Allocate the work arrays.
    real ( kind = 8 ), allocatable, dimension ( : ) :: work
    real ( kind = 8 ), allocatable, dimension ( : ) :: wsave
    real (dp), allocatable, dimension ( : ) :: Hann
      contains
!      ------------------------------------
    Subroutine RFTransform_su(N_T)
! Initialize real-FFT transform and allocate workspace
! N_Time= Length of input array in time domain.
  implicit none
  integer , intent (in) :: N_T
  integer ( kind = 4 ) :: ier, i
  N_Time=N_T
  lensav = N_Time + int ( log ( real ( N_Time, kind = 8 ) ) / log ( 2.0D+00 ) ) + 4
  lenwrk = N_Time ! Agrees with example given
  allocate ( work(1:lenwrk) )
  allocate ( wsave(1:lensav) )
  allocate ( Hann(1:N_Time) )
!    write(2,*) 'FFTransform_su',lenwrk,lensav,N_t
  call Rfft1i ( N_Time, wsave, lensav, ier )
  N_nu=N_Time/2
  HW_size=min(N_Time/10,32)
  !write(2,*) 'HW_size=',HW_size,pi
    Hann(:)=1.
    Do i=0,HW_size
        Hann(i+1) = sin(0.5*i*pi/HW_size)**2
        Hann(N_Time-i)=Hann(i+1)
    EndDo
  return
  End Subroutine RFTransform_su
!      ------------------------------------
!      ------------------------------------
    Subroutine DAssignFFT()
       deallocate ( work )
       deallocate ( wsave )
       deallocate ( Hann )
    End Subroutine DAssignFFT
!      ------------------------------------
  Subroutine RFTransform_CF(A,Cnu)
!  Subroutine FFTransform_FF(N_t,A,N_f,f,t_shft,C, wsave, lensav, work, lenwrk)
! Fourier transform real timetrace to complex-frequency
! A= input array in time domain, length=N_t
! Cnu=output array (complex) in frequency, length=N_f
!    use constants, only : dp,pi,ci
  implicit none
  complex(dp), intent(out) :: Cnu(0:N_nu)
  REAL(dp), intent(in) :: A(1:N_Time)
  !
  integer :: i
  REAL(dp) :: D(1:N_Time)
  complex(dp) :: C(0:N_nu)
  integer ( kind = 4 ) :: ier, inr=1, lenR
  !ipi=ci*pi  ! This is the real constant i*pi
  !DO I=inui,inum
  !  ph_shft=exp(ipi*t_shft*i/nuTrace_dim)
  !  Cnu(I)=C(I)*F(I) *ph_shft
  !ENDDO
!
!  Compute the FFT coefficients.
!
  !  write(*,*) 'FFTransform_FF',lenwrk,lensav
  D(:)=A(:) !  *Hann(:)
  call Rfft1f ( N_Time, inR, D, N_Time, wsave, lensav, work, lenwrk, ier )
  CALL R2C(D,Cnu)
!
  return
  End Subroutine RFTransform_CF
!
!      ------------------------------------
  Subroutine RFTransform_CF_Filt(A,F,t_shft,Cnu)
!  Subroutine FFTransform_FF(N_t,A,N_f,f,t_shft,C, wsave, lensav, work, lenwrk)
! Fourier transform to frequency, apply filter and timeshift
! A= input array in time domain, length=N_t
! F=input real-valued filter function  in frequency
! t_shft= shift of pulse in time in units of sample-time
! Cnu=output array (complex) in frequency, length=N_f; after filter
!    use constants, only : dp,pi,ci
  implicit none
  complex(dp), intent(out) :: Cnu(0:N_nu)
  REAL(dp), intent(in) :: A(1:N_Time),t_shft,F(0:N_nu)
  integer :: i
  REAL(dp) :: D(1:N_Time)
  complex(dp) :: C(0:N_nu),ph_shft,ipi
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inR
  integer ( kind = 4 ) lenR
  ipi=ci*pi  ! This is the complex constant i*pi
    !
  D(:)=A(:) ! *Hann(:)
  inR = 1
  lenR = N_Time
  call Rfft1f ( N_Time, inR, D, N_Time, wsave, lensav, work, lenwrk, ier )
  CALL R2C(D,C)
!  write(2,*) '!RFTransform_CF_Filt', t_shft, N_nu
!
!    write(2,*) '0',abs(B(0)),ph(B(0)), f(0)
!    write(2,*) '1',abs(B(1)),ph(B(1)), f(1)
!    write(2,*) 'n/2',abs(B(N_f/2)),ph(B(N_f/2)) , f(n_f/2)
!    write(2,*) 'n-1',abs(B(N_f-1)),ph(B(N_f-1))
!    write(2,*) 'n',abs(B(N_f)),ph(B(N_f))
  DO I=0,N_nu
    ph_shft=exp(-ipi*t_shft*i/N_nu)
    Cnu(I)=C(I)*F(I) *ph_shft
  ENDDO
  return
  End Subroutine RFTransform_CF_Filt
!
!      ------------------------------------
!      ------------------------------------
  Subroutine RFTransform_CF2RT(Cnu,RD)
! Fourier transform back to time domain
! Cnu=input array (complex) in frequency, length=N_f
! RD= Real output array in time domain, length=N_t
!  use constants, only : dp,ci
  implicit none
  complex(dp), intent(in) :: Cnu(0:N_nu)
  real(dp), intent(out) :: RD(1:N_Time)
  !
  integer ( kind = 4 ) :: ier, inr=1, lenR
  !
  lenR = N_Time
!
!  Compute inverse FFT of coefficients.
  CALL C2R(RD,Cnu)
  call Rfft1b ( N_Time, inR, RD, lenR, wsave, lensav, work, lenwrk, ier )
!
  return
  End Subroutine RFTransform_CF2RT
!
!      ------------------------------------
  Subroutine RFTransform_CF2CT(Cnu,CD)
! Fourier transform back to time domain
! Cnu=input array (complex) in frequency, length=N_f
! CD= complex output array in time domain, length=N_t
!  use constants, only : dp,ci
  implicit none
  complex(dp), intent(in) :: Cnu(0:N_nu)
  complex(dp), intent(out) :: CD(1:N_Time)
  !
  INTEGER*4 i
  REAL*8 D(1:N_Time),ID(1:N_Time)
  complex*16 C(0:N_nu),IC(0:N_nu)
  integer ( kind = 4 ) :: ier, inr=1, lenR
  !
  lenR = N_Time
  IC(:)=ci*Cnu(:)
  c(:)=Cnu(:)
  !DO I=0,N_nu ! inui,inum
  !  C(i)=Cnu(I)
  !  IC(I)= ci*Cnu(i)
  !ENDDO
!
!  Compute inverse FFT of coefficients.
  CALL C2R(D,C)
  call Rfft1b ( N_Time, inR, D, lenR, wsave, lensav, work, lenwrk, ier )
  CALL C2R(ID,IC)
  call Rfft1b ( N_Time, inR, ID, lenR, wsave, lensav, work, lenwrk, ier )
!
  DO I=1,N_Time
    CD(i)=D(i) + ci*ID(i)
  ENDDO
  return
  End Subroutine RFTransform_CF2CT
!
!------------------------------------------------
Pure Subroutine R2C(R,C)
! Nc=Nr/2 note that zero and max frequency components are always real
    IMPLICIT NONE
    INTEGER I
    REAL*8, intent(in) :: R(1:N_Time)
    COMPLEX*16, intent(out) :: C(0:N_nu)
  C(0)=R(1)
  DO I=1,N_nu-1
    C(I)=CMPLX(R(2*I),R(2*I+1))
  ENDDO
  C(N_nu)=R(2*N_nu)
  RETURN
END
Pure Subroutine  C2R(R,C)
! Nc=Nr/2 note that zero and max frequency components are always real
    IMPLICIT NONE
    INTEGER I
    REAL*8, intent(out) :: R(1:N_Time)
    COMPLEX*16, intent(in) :: C(0:N_nu)
  R(1)=DREAL(C(0))
  DO I=1,N_nu-1
    R(2*I)=DREAL(C(I))
    R(2*I+1)=DIMAG(C(I))
  ENDDO
  R(2*N_nu)=DREAL(C(N_nu))
  RETURN
END
!      ------------------------------------
Pure function ph(zzz)
!    use constants, only : pi
	implicit none
	complex*16, intent(in) :: zzz
	real*8 :: re,im,ph
!	common / constants / CI,pi
!	  COMPLEX*16 CI
!	  REAL*8 pi
	re=REALpart(zzz)
	im=IMAGpart(zzz)
	if(im.eq.0) then
		ph=0.
	else
		ph=datan(im/re)
	endif
	if(re.lt.0.) ph=ph+pi
	if(ph.gt.pi) ph=ph-2*pi
	ph=ph/pi
	return
	end
!      ------------------------------------
!      ------------------------------------
  Subroutine DownSamlple(E_t,E,padding,tTrace_dim,FF_dim, filt, E_nu_dwn, inui,inum)

!  Subroutine FFTransform_FF(N_t,A,N_f,f,t_shft,C, wsave, lensav, work, lenwrk)
! Fourier transform to frequency, apply filter and timeshift
! A= input array in time domain, length=N_t
! F=input filter function (complex) in frequency
! C_nu=output array (complex) in frequency, length=N_f; after filter
!    use constants, only : dp,pi,ci
  implicit none
  integer, intent(in) :: inui,inum,padding,tTrace_dim, FF_dim
  complex(dp), intent(out) :: E_nu_dwn(inui:inum)
  REAL(dp), intent(in) :: E_t(1:tTrace_dim)
  REAL(dp), intent(inout) :: E(1:FF_dim)
  complex(dp), intent(in) :: filt(inui:inum)
  integer :: i,nu_dim
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inR
  integer ( kind = 4 ) lenR
!
!  Compute the FFT coefficients.
!
  inR = 1
  lenR = FF_dim
  nu_dim=ff_dim/2
  E(1:padding)=0.0d0 ;  E(padding+1:tTrace_dim+padding)=E_t(1:tTrace_dim) ;  E(tTrace_dim+padding+1:FF_dim)=0.0d0
  call Rfft1f (FF_dim, inR, E, FF_dim, wsave, lensav, work, lenwrk, ier )
  DO I=inui,inum
    E_nu_dwn(I)=filt(I)*CMPLX(E(2*I),E(2*I+1))
  ENDDO
  return
  End Subroutine DownSamlple
!
      End Module FFT
! ==========================================================
