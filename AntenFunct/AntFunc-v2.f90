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
   use AntFunCconst
   Implicit none
   Real(dp) :: Cheb(m_dim), Cheb_s(m_dim), Lgdr(j_dim)
   Real(dp) :: Freq0, Thet0, Phi0, Vp_r, Vp_i
   Real(dp) :: Freq, Thet, Phi, Vt_r, Vt_i
   Complex(dp) :: J_0p(Freq_min:Freq_max),J_0t(Freq_min:Freq_max),J_1p(Freq_min:Freq_max),J_1t(Freq_min:Freq_max)
   Complex(dp) :: PhIn_p(m_dim,20), PhIn_t(m_dim,20)
   Complex(dp) :: PhIn_0p(m_dim), PhIn_0t(m_dim)
   Real(dp) :: th, x, dx, dthet, S1p25, c1p25
   Complex(dp) :: Vp, Vt
   Integer :: nxx, it, nthet, i_freq, j
   Character(len=20) :: release
   INTEGER :: DATE_T(8),i
   CHARACTER*12 :: REAL_C(3)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='AntFunc-v2.out')
   release='version 7 (15Aug-19)'
   write(2,"(3x,5(1H-),1x,'Antenna Fun release of ',A22,25(1H-))") release
   CALL DATE_AND_TIME (REAL_C(1),REAL_C(2),REAL_C(3),DATE_T)
   WRITE(2,"(3X,5(1H-),1x,'run on ',I2,'/',I2,'/',I4,' , started at ',&
       I2,':',I2,':',I2,'.',I3,1X,25(1H-))") &
       DATE_T(3),DATE_T(2),DATE_T(1),(DATE_T(i),i=5,8)
   !
   S1p25=sin(1.25*pi/180.) ! needed for average angle first and last half bin
   c1p25=cos(1.25*pi/180.)
   Open(Unit=10,Status='old',Action='read',File='LBA_Vout_phi.txt')
   Open(Unit=11,Status='old',Action='read',File='LBA_Vout_theta.txt')
   read(10,*)
   read(11,*)
   PhIn_p(:,:)=0.
   PhIn_t(:,:)=0.
   Freq=0.
   nxx=0
   Do while(nxx.eq.0)
      Vp=cmplx(Vp_r,Vp_i)
      Vt=cmplx(Vt_r,Vt_i)
      it=nint(1. + thet/5.)
      i_freq=nint(Freq)
      If(i_freq .lt. Freq_min) goto 1
      If(i_freq .gt. Freq_max) exit
      !write(2,*) Freq, Thet, Phi, Vp_r, Vp_i
      Call Cheb_odd(Phi*pi/180.,Cheb)
      Call Cheb_odd((Phi+90.)*pi/180.,Cheb_s)
      If(Phi .lt. 359.) then ! integrate over \phi since \phi=0 is same as 360
         PhIn_p(:,it)=PhIn_p(:,it)+Cheb(:)*Vp/18.  ! 1/18=2*10/360
         PhIn_t(:,it)=PhIn_t(:,it)+Cheb_s(:)*Vt/18.  ! 1/18=2*10/360
      else ! Phi=360, print and initialize integral
         Write(2,"(F5.1,F7.1,4F10.4,4F10.4)") Freq, Thet,  PhIn_p(:,it), PhIn_t(:,it)
         If(Thet.gt.89) then ! do \theta integral
            dthet=5.
            nthet=nint(90./dthet)
            dx=0.5*S1p25*pi*dthet/180.  ! get the average for the first bin (1/2 size)
            Call Legendre_odd(C1p25,Lgdr) ! x=cos(theta)
            Do j=1,j_dim
               Ant_p(i_freq,j,1:m_dim)= PhIn_p(1:m_dim,1)*Lgdr(j)*(4.*j-1.)*dx ! norm: (2n+1)/18 where n=2j-1; 0.5 since interval at th=0 is only half width
               Ant_t(i_freq,j,1:m_dim)= PhIn_t(1:m_dim,1)*Lgdr(j)*(4.*j-1.)*dx ! norm: (2n+1)/18 where n=2j-1; 0.5 since interval at th=0 is only half width
            Enddo
            Do it=1,nthet ! should integrate over d[x]=d[cos(th)]=sin(th) d[th]
            ! because of the 5 deg spacing the accuracy of this integral is about 1% of the largest number
               If(it.lt.nthet) then
                  th=it*dthet*pi/180.
                  x =cos(th)
                  dx=sin(th)*pi*dthet/180.
               else ! take care of last half bin
                  x =S1p25
                  dx=0.5*C1p25*pi*dthet/180.
               endif
               Call Legendre_odd(x,Lgdr) ! x=cos(theta)
               Do j=1,j_dim
                  Ant_p(i_freq,j,1:m_dim)=Ant_p(i_freq,j,1:m_dim) + PhIn_p(1:m_dim,it)*Lgdr(j)*(4.*j-1.)*dx ! norm: (2n+1)/18 where n=2j-1; 0.5 since interval at th=0 is only half width
                  Ant_t(i_freq,j,1:m_dim)=Ant_t(i_freq,j,1:m_dim) + PhIn_t(1:m_dim,it)*Lgdr(j)*(4.*j-1.)*dx ! norm: (2n+1)/18 where n=2j-1; 0.5 since interval at th=0 is only half width
               Enddo
            Enddo
            write(2,"(F5.1,A,4(2F10.4,2x))") freq,' Phi(re,im):=',Ant_p(i_freq,:,1)
            write(2,"(F5.1,A,4(2F10.4,2x))") freq,' Thet(re,im):=',Ant_t(i_freq,:,1)
            PhIn_p(:,:)=0.
            PhIn_t(:,:)=0.
            !stop
         Endif
      Endif
   1  continue
      read(10,*,iostat=nxx) Freq0, Thet0, Phi0, Vp_r, Vp_i
      read(11,*,iostat=nxx) Freq, Thet, Phi, Vt_r, Vt_i
      If( nint(Freq0*1000.+ Thet0*100.+ Phi0) .ne. nint(Freq*1000.+ Thet*100.+ Phi)) then
         write(2,*) 'mismatch',nint(Freq0*1000.+ Thet0*100.+ Phi0),nint(Freq*1000.+ Thet*100.+ Phi)
         stop
      endif
   Enddo
   !stop
   !  Checking
   rewind(unit=10)
   read(10,*)
   rewind(unit=11)
   read(11,*)
   Freq=0.
   nxx=0
   PhIn_p(:,:)=0.
   PhIn_t(:,:)=0.
   PhIn_0p(:)=0.
   PhIn_0t(:)=0.
   Write(2,"(6x,A,7x,A,11x,A)") 'Freq, Thet, Phi,','V(data)', 'J_0(i_freq)'
   Do while(nxx.eq.0)
      i_freq=nint(Freq)
      it=1
      If(i_freq .lt. Freq_min) goto 2
      Call AntFun(thet,phi-45,J_0p,J_0t,J_1p,J_1t) ! Jones; gives _0 & _1 voltages when contracted with (t,p) polarized field
      Call Cheb_odd(Phi*pi/180.,Cheb)
      Call Cheb_odd((Phi+90.)*pi/180.,Cheb_s)
      If(Phi .lt. 359.) then ! integrate over \phi since \phi=0 is same as 360
         PhIn_p(:,it)=PhIn_p(:,it)+Cheb(:)*Vp/18.  ! 1/18=2*10/360
         PhIn_t(:,it)=PhIn_t(:,it)+Cheb_s(:)*Vt/18.  ! 1/18=2*10/360
         PhIn_0p(:)=PhIn_0p(:)+Cheb(:)*J_0p(i_freq)/18.  ! 1/18=2*10/360
         PhIn_0t(:)=PhIn_0t(:)+Cheb_s(:)*J_0t(i_freq)/18.  ! 1/18=2*10/360
      else ! Phi=360, print and initialize integral
         Write(2,"(F5.1,F7.1,4(2F10.4,2x))") Freq, Thet,  PhIn_p(:,it), PhIn_t(:,it)
         Write(2,"(5x,A,4(2F10.4,2x))") 'model:', PhIn_0p(:), PhIn_0t(:)
         PhIn_p(:,:)=0.
         PhIn_t(:,:)=0.
         PhIn_0p(:)=0.
         PhIn_0t(:)=0.
      Endif
      !
      If(Phi.gt.90) goto 2
      j=nint(Phi)
      !If(j.ne.0 .and. j.ne.90) goto 2
      If(i_freq .gt. Freq_max) exit
      Call AntFun(thet,phi-45,J_0p,J_0t,J_1p,J_1t) ! Jones; gives _0 & _1 voltages when contracted with (t,p) polarized field
      If(abs(Vp-J_0p(i_freq)) .gt. abs(Vp+J_0p(i_freq))/30.) then
         If(abs(J_0p(i_freq)) .gt. 0.001) &
         write(2,"(A,3F6.1,2F9.4,2x,2F9.4)") ' phi', Freq, Thet, Phi, Vp, J_0p(i_freq)
      endif
      If(abs(Vt-J_0t(i_freq)) .gt. abs(Vt+J_0t(i_freq))/30.) then
         If(abs(J_0t(i_freq)) .gt. 0.001) &
         write(2,"(A,3F6.1,2F9.4,2x,2F9.4)") 'thet', Freq, Thet, Phi, Vt, J_0t(i_freq)
      endif
   2  continue
      read(10,*,iostat=nxx) Freq, Thet, Phi, Vp_r, Vp_i
      read(11,*,iostat=nxx) Freq, Thet, Phi, Vt_r, Vt_i
      Vp=cmplx(Vp_r, Vp_i)
      Vt=cmplx(Vt_r, Vt_i)
      if(i_freq .eq. (Freq_min+1)) exit
   Enddo
   !
   Close(unit=10)
   Close(unit=11)
   !
   i_Freq=20
   Write(2,"(A,i2,4x,A,20x,A)") '@frequency=',i_Freq,'J_0p, J_0t @ phi=0','J_0p, J_0t'
   Do j=0,90,5
      thet=j
      phi=0
      Call AntFun(thet,phi-45,J_0p,J_0t,J_1p,J_1t) ! Jones; gives _0 & _1 voltages when contracted with (t,p) polarized field
      write(2,"(A,F6.1,2F9.4,2x,2F9.4)", ADVANCE='NO') 'thet=',thet,J_0p(i_Freq),J_0t(i_Freq)
      phi=90
      Call AntFun(thet,phi-45,J_0p,J_0t,J_1p,J_1t) ! Jones; gives _0 & _1 voltages when contracted with (t,p) polarized field
      write(2,"(A,2F9.4,2x,2F9.4)") ', phi=90:',J_0p(i_Freq),J_0t(i_Freq)
   Enddo
   stop
End Program AntennaFunction
