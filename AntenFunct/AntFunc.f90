    Include '../ConstantsModules-v7.f90'
    !   Module constants
    !   Module DataConstants
    Include 'SpecialFie.f90'
Program AntennaFunction
!     Phi=0 is along one of the arms; theta=0 virtical
!
!     Parametrixe Antenna function in therms of spherical harmonics, d\Omega =\sin(\theta) d\theta d\phi
!
!     Due to symmetry the odd Chebyshev polynominals (in sin(\phi)) are expected for \phi-dependence for \phi polarized incoming and
!        odd Chebyshev polynominals (in cos(\phi)) for \theta-polarized incoming waves, and
!        odd Legendre polynominals for x=cos(\theta) dependence.
!     Note 0<\theta<90 so integral runs only over half of the cos(\theta) interval and orthoganality of even and odd poly is broken.
!
!     Legendre polynominals obey \int_{-1}^{+1} P_n(x) P_m(x) dx= \delta(n,m) 2/(2n+1)
!     where x=\cos(\theta) and dx=\sin(\theta) d\theta
!     Chebyshev polynominals obey \int_{-1}^{+1} T_n(y) T_m(y) d\phi= \delta(n,m) \pi/2 (=\pi when n=0=m)
!     where y=\cos(\phi) and dx=\sin(\phi) d\phi or d\phi=(1-x^2)^{-1/2} dx
!     T_n(sin\phi)=sin(m\phi) which gives the relation to spherical harmonics
!
!     T_0(y)=1
!     T_1(y)=y
!     T_2(y)= 2y^2 - 1
!     T_3(y)= 4y^3 - 3y
!     T_4(y)= 8y^4 - 8y^2 + 1
!     T_5(y)=16y^5 -20y^3 + 5y
!     T_6(y)=32y^6 -48y^4 +18y^2 - 1
!
!     P_0(x) =  1
!     P_1(x) =  x
!     P_2(x) =( 3x^2 - 1  )/2
!     P_3(x) =( 5x^3 - 3x )/2
!     P_4(x) =(35x^4 -30x^2 + 3  )/8
!     P_5(x) =(63x^5 -70x^3 +15x )/8
!     P_6(x) =(231x^6 - 315x^4 +105x^2 - 5 )/16
!     P_7(x) =(429x^7 -693x^5 +315x^3 -35x )/16
!
   use constants, only : dp,pi,ci
   Implicit none
   Real(dp) :: Ch0,Ch1,Ch2,Ch3,Ch4,Ch5,Ch6
   Real(dp) :: Freq, Thet, Phi, Vp_r, Vp_i
   Real(dp) :: Freq0, Thet0, Phi0, Vt_r, Vt_i
   Real(dp) :: Int0r, Int1r, Int2r, Int3r, Int4r, Int5r, Int6r
   Real(dp) :: Int0i, Int1i, Int2i, Int3i, Int4i, Int5i, Int6i
   Real(dp) :: T1r(20), T3r(20), T1i(20), T3i(20)
   Real(dp) :: th,x, dx, Lgdr0,Lgdr1,Lgdr2,Lgdr3,Lgdr4,Lgdr5,Lgdr6,Lgdr7
   Real(dp) :: Y11r,Y31r,Y51r, Y71r, Y13r,Y33r,Y53r, Y73r, check, dthet
   Real(dp) :: Y11i,Y31i,Y51i,Y13i,Y33i,Y53i
   Integer :: nxx, it, nthet
   Character(len=20) :: release
   INTEGER :: DATE_T(8),i
   CHARACTER*12 :: REAL_C(3)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='AntFunc.out')
   release='version 7 (15Aug-19)'
   write(2,"(3x,5(1H-),1x,'Antenna Fun release of ',A22,25(1H-))") release
   CALL DATE_AND_TIME (REAL_C(1),REAL_C(2),REAL_C(3),DATE_T)
   WRITE(2,"(3X,5(1H-),1x,'run on ',I2,'/',I2,'/',I4,' , started at ',&
       I2,':',I2,':',I2,'.',I3,1X,25(1H-))") &
       DATE_T(3),DATE_T(2),DATE_T(1),(DATE_T(i),i=5,8)
   !
   Open(Unit=10,Status='old',Action='read',File='LBA_Vout_phi.txt')
   read(10,*)
   !Thet0=-11
   !Phi0=-11
   !Freq0=-11
   Int0r=0  ; Int0i=0
   Int1r=0  ; Int1i=0
   Int2r=0  ; Int2i=0
   Int3r=0  ; Int3i=0
   Int4r=0  ; Int4i=0
   Int5r=0  ; Int5i=0
   Int6r=0  ; Int6i=0
   read(10,*,iostat=nxx) Freq, Thet, Phi, Vp_r, Vp_i
   Do while(nxx.eq.0)
      Vp_r=Vp_r/18.      ! to normalize the \phi integrals
      Vp_i=Vp_i/18.      ! to normalize the \phi integrals
      Call Cheb(Phi*pi/180.,Ch0,Ch1,Ch2,Ch3,Ch4,Ch5,Ch6)
      !Vp_r=Ch3/18.
      !Vp_i=Ch4/18.
      !write(2,*) Freq, Thet, Phi, Vp_r, Vp_i
      If(Phi .lt. 359.) then ! integrate over \phi since \phi=0 is same as 360
         Int0r=Int0r+Ch0*Vp_r  ; Int0i=Int0i+Ch0*Vp_i
         Int1r=Int1r+Ch1*Vp_r  ; Int1i=Int1i+Ch1*Vp_i
         Int2r=Int2r+Ch2*Vp_r  ; Int2i=Int2i+Ch2*Vp_i
         Int3r=Int3r+Ch3*Vp_r  ; Int3i=Int3i+Ch3*Vp_i
         Int4r=Int4r+Ch4*Vp_r  ; Int4i=Int4i+Ch4*Vp_i
         Int5r=Int5r+Ch5*Vp_r  ; Int5i=Int5i+Ch5*Vp_i
         Int6r=Int6r+Ch6*Vp_r  ; Int6i=Int6i+Ch6*Vp_i
      else ! Phi=360, print and initialize integral
            !Write(2,*) Freq, Thet, Int0r, Int1r, Int2r, Int3r, Int4r, Int5r, Int6r, &
            !   ' ; ',Int0i, Int1i, Int2i, Int3i, Int4i, Int5i, Int6i
         Write(2,"(F5.1,F7.1,3F10.4,A,4F10.4)") Freq, Thet,  Int1r,  Int3r,  Int5r,  &
            ' ; ', Int1i,  Int3i, Int5i
         it=nint(1 + thet/5.)
         T1r(it)=Int1r
         T3r(it)=Int3r
         T1i(it)=Int1i
         T3i(it)=Int3i
         If(Thet.gt.89) then ! do \theta integral
            !Call Legendre(1.d0,Lgdr0,Lgdr1,Lgdr2,Lgdr3,Lgdr4,Lgdr5,Lgdr6,Lgdr7) ! x=cos(theta)
            Y11r=0. ! 5*T1r(1)*Lgdr1*3/18.  ! norm: *(2n+1)/18; 0.5 since interval at th=0 is only half width
            Y31r=0. ! 5*T1r(1)*Lgdr3*7/18.
            Y51r=0. ! 5*T1r(1)*Lgdr5*11/18.
            Y71r=0. ! 5*T1r(1)*Lgdr7*15/18.
            Y13r=0. ! 5*T3r(1)*Lgdr1*3/18.  ! norm: *(2n+1)/18; 0.5 since interval at th=0 is only half width
            Y33r=0. ! 5*T3r(1)*Lgdr3*7/18.
            Y53r=0. ! 5*T3r(1)*Lgdr5*11/18.
            Y73r=0. ! 5*T3r(1)*Lgdr7*15/18.
            check=0.
            dthet=5.
            nthet=nint(90./dthet)
            Do it=2,nthet ! should integrate over d[x]=d[cos(th)]=sin(th) d[th]
            ! because of the 5 deg spacing the accuracy of this integral is about 1% of the largest number
               th=(it-1)*dthet*pi/180.
               x=cos(th)
               dx=sin(th)*pi*dthet/180.
               Call Legendre(x,Lgdr0,Lgdr1,Lgdr2,Lgdr3,Lgdr4,Lgdr5,Lgdr6,Lgdr7) ! x=cos(theta)
               !T1r(it)=Lgdr1  !  will have to add half of the 90^deg bin
               !T3r(it)=Lgdr5  !  will have to add half of the 90^deg bin
               !check=check+dx
               Y11r=Y11r + T1r(it)*Lgdr1*3*dx
               Y31r=Y31r + T1r(it)*Lgdr3*7*dx
               Y51r=Y51r + T1r(it)*Lgdr5*11*dx
               Y71r=Y71r + T1r(it)*Lgdr7*15*dx
               Y13r=Y13r + T3r(it)*Lgdr1*3*dx
               Y33r=Y33r + T3r(it)*Lgdr3*7*dx
               Y53r=Y53r + T3r(it)*Lgdr5*11*dx
               Y73r=Y73r + T3r(it)*Lgdr7*15*dx
            Enddo
            !check=check+0.5*pi/36.
            write(2,"(F5.1,A,4F10.4,A,4F10.4)") freq,' real:=',Y11r,Y31r,Y51r,Y71r,' , ',Y13r,Y33r,Y53r,Y73r
            !write(2,*) 'check',check,nxx
            !write(2,*) freq,' Imag:=',Y11r,Y31r,Y51r,' , ',Y13r,Y33r,Y53r
         Endif
         !Endif
         Int0r=0  ; Int0i=0
         Int1r=0  ; Int1i=0
         Int2r=0  ; Int2i=0
         Int3r=0  ; Int3i=0
         Int4r=0  ; Int4i=0
         Int5r=0  ; Int5i=0
         Int6r=0  ; Int6i=0
      Endif
      read(10,*,iostat=nxx) Freq, Thet, Phi, Vp_r, Vp_i
   Enddo
   !
   !Call AntFun_Inv(thet,phi,Ji_p0,Ji_t0,Ji_p1,Ji_t1)
   stop
End Program AntennaFunction
! =================================
Subroutine AntFun(thet,phi,J_0p,J_0t,J_1p,J_1t)
! Jones; gives _0 & _1 voltages when contracted with (t,p) polarized field
!  Angle \phi is calculated from the orientation of the 0-LOFAR dipole,
!        \theta is from the local vertical.
!  The elements of the Jones matrix are complex and functions of frequency (in 1MHz steps)
!  Uses the parametrized antenna function J_i,[t,p]=\sum_jm Ant_[t,p] (freq,j,m) Cheb(m,phi) Lgdr(j,thet)
!  where depending on i=[0,1] and [t,p] phi is shifted by 90; [t,p] denote polarization incoming wave.
   use constants, only : dp,pi,ci
   Implicit none
   Integer :: j,m,i_freq
   real(dp) :: Cheb(3),Cheb_s(3),Lgdr(5)
   Call Cheb_odd(Phi*pi/180.,Cheb)
   Call Cheb_odd((Phi-90.)*pi/180.,Cheb_s)
   x=cos(thet*pi/180.)
   Call Legendre_odd(x,Lgdr) ! x=cos(theta)
   Do i_freq=20,90
      J_0p(i_freq)=0.
      Do j=1,4 ! corresponding to 1,3,5,7, Legendre polynominals
         Do m=1,2 ! corresponding to 1,3,5 Chebyshev polynominals
            J_0p(i_freq)=J_0p(i_freq) + Ant_p(i_freq,j,m)*Lgdr(j)*Cheb(m)
         Enddo ! m=1,2
      Enddo ! j=1,4
   Enddo ! i_freq
End Subroutine AntFun
