! ============================================
!=======================================
Subroutine Cheb_odd(Phi,Cheb)
!     Chebyshev polynominals obey \int_{0}^{+\pi} T_n(y) T_m(y) d\phi= \delta(n,m) \pi/2 (=\pi when n=0=m)
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
   use constants, only : dp,pi,ci
   use AntFuncConst
   Implicit none
   real(dp), intent(in) :: Phi  ! in radian
   real(dp), intent(out) :: Cheb(m_dim)  ! corresponding to 1,3,5 Chebyshev polynominals
   Real(dp) :: y,y2,y3,y5
   y=sin(Phi)
   y2=y*y
   y3=y*y2  ; y5=y3*y2
   !Ch0=1.
   Cheb(1)=y
   !Ch2=2*y2 - 1.
   Cheb(2)=4*y3 -3*y
   !!Ch4=8*y4 - 8*y2 + 1
   !Cheb(3)=16*y5 - 20*y3 + 5*y
   !!Ch6=32*y6 - 48*y4 + 18*y2 -1
   Return
End Subroutine Cheb_odd
! ===============
Subroutine Legendre_odd(x,Lgdr)
!     d\Omega =\sin(\theta) d\theta d\phi
!     Legendre polynominals obey \int_{-1}^{+1} P_n(x) P_m(x) dx= \delta(n,m) 2/(2n+1)
!     where x=\cos(\theta) and dx=\sin(\theta) d\theta
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
   real(dp), intent(in) :: x
   real(dp), intent(out) :: Lgdr(j_dim)  ! corresponding to 1,3,5,7, Legendre polynominals
   real(dp) :: x2,x3,x4,x5,x6,x7
   x2=x*x
   x3=x*x2  ; x5=x3*x2  ; x7=x5*x2
   !x4=x2*x2 ; x6=x4*x2
   !Lgdr0=1.
   Lgdr(1)=x
   !Lgdr2=(  3*x2 -  1.  )/2.
   Lgdr(2)=(  5*x3 -  3*x )/2.
   !Lgdr4=( 35*x4 - 30*x2 +  3  )/8.
   Lgdr(3)=( 63*x5 - 70*x3 + 15*x )/8.
   !Lgdr6=(231*x6 -315*x4 +105*x2 - 5 )/16.
   Lgdr(4)=(429*x7 -693*x5 +315*x3 -35*x )/16.
   Return
End Subroutine Legendre_odd
! =================================
Subroutine AntFun(thet_d,phi_d)
! Jones; gives _0 & _1 voltages when contracted with (t,p) polarized field
!  Assumed that the 1-dipole (=X) is rotated over +90 w.r.t. 0-dipole (=Y)
!  See Figure 3.1 in:
! https://books.google.nl/books?id=C3FyDwAAQBAJ&pg=PA38&lpg=PA38&dq=lofar+x+dipole+and+y+dipole&source=bl&ots=D1BmD_X3px&sig=ACfU3U0wSFRcijLcXbDgvVTXV06g_2VO0Q&hl=en&sa=X&ved=2ahUKEwjpos_7nfPkAhXDIlAKHYrnC3cQ6AEwAHoECAkQAQ#v=onepage&q=lofar%20x%20dipole%20and%20y%20dipole&f=false
!     Which says that the X-dipole (odd) is NE(-)--SW(+), and Y (even) is NW(-)--SE(+)
!  Brian:  I see from the documentation that odd is the X dipole and even is Y dipole for LBA-outer
!  Angle \phi is calculated from the orientation of the 0-LOFAR dipole,
!        \theta is from the local vertical.
!  The elements of the Jones matrix are complex and functions of frequency (in 1MHz steps)
!  Uses the parametrized antenna function J_i,[t,p]=\sum_jm Ant_[t,p] (freq,j,m) Cheb(m,phi) Lgdr(j,thet)
!  where depending on i=[0,1] and [t,p] phi is shifted by 90; [t,p] denote polarization incoming wave.
   !Call AntFun_Inv(thet,phi,Ji_p0,Ji_t0,Ji_p1,Ji_t1)
   use constants, only : dp,pi,ci
   use AntFunCconst
   Implicit none
   Real(dp), intent(in) :: thet_d, phi_d  ! in degree phi=0=North, phi=90=East, thet=0=vertical
   !Complex(dp), intent(out) :: J_0p(Freq_min:Freq_max),J_0t(Freq_min:Freq_max),J_1p(Freq_min:Freq_max),J_1t(Freq_min:Freq_max)
   Integer :: j,m,i_freq
   real(dp) :: Cheb(m_dim), Cheb_s(m_dim), Lgdr(j_dim), Phi_r
   Real(dp) :: x
   ! Phi_r in radian, where Phi_r=0=along Y-dipole
   Phi_r=(phi_d+45)*pi/180.
   Call Cheb_odd(Phi_r,Cheb)
   Call Cheb_odd((Phi_r+pi/2),Cheb_s)
   x=cos(thet_d*pi/180.)
   Call Legendre_odd(x,Lgdr) ! x=cos(theta)
   Do i_freq=Freq_min,Freq_max
      J_0p(i_freq)=0.
      J_1p(i_freq)=0.
      J_0t(i_freq)=0.
      J_1t(i_freq)=0.
      Do j=1,j_dim ! corresponding to 1,3,5,7, Legendre polynominals
         Do m=1,m_dim ! corresponding to 1,3,5 Chebyshev polynominals
            J_0p(i_freq)=J_0p(i_freq) + Ant_p(i_freq,j,m)*Lgdr(j)*Cheb(m)
            J_0t(i_freq)=J_0t(i_freq) + Ant_t(i_freq,j,m)*Lgdr(j)*Cheb_s(m)
            J_1p(i_freq)=J_1p(i_freq) - Ant_p(i_freq,j,m)*Lgdr(j)*Cheb_s(m) ! sign change because of rotation over 180 deg
            J_1t(i_freq)=J_1t(i_freq) + Ant_t(i_freq,j,m)*Lgdr(j)*Cheb(m)
         Enddo ! m=1,2
      Enddo ! j=1,4
   Enddo ! i_freq
   Return
End Subroutine AntFun
! =================================
Subroutine AntFun_Inv(thet_d,phi_d)
! Inverse Jones; gives _p & _t polarized fields when contracted with (0.1) voltages
!  Assumed that the 1-dipole is rotated over +90 w.r.t. 0-dipole
!  Angle \phi is calculated from the orientation of the 0-LOFAR dipole,
!        \theta is from the local vertical.
!  The elements of the Jones matrix are complex and functions of frequency (in 1MHz steps)
!  Uses the parametrized antenna function J_i,[t,p]=\sum_jm Ant_[t,p] (freq,j,m) Cheb(m,phi) Lgdr(j,thet)
!  where depending on i=[0,1] and [t,p] phi is shifted by 90; [t,p] denote polarization incoming wave.
   !Call AntFun_Inv(thet,phi,Ji_p0,Ji_t0,Ji_p1,Ji_t1)
   use constants, only : dp
   use AntFunCconst
   Implicit none
   Real(dp), intent(in) :: thet_d, phi_d  ! in degree
   !Complex(dp), intent(out) :: Ji_p0(Freq_min:Freq_max),Ji_t0(Freq_min:Freq_max),Ji_p1(Freq_min:Freq_max),Ji_t1(Freq_min:Freq_max)
   !Complex(dp) :: J_0p(Freq_min:Freq_max),J_0t(Freq_min:Freq_max),J_1p(Freq_min:Freq_max),J_1t(Freq_min:Freq_max)
   Complex(dp) :: Det(Freq_min:Freq_max)
   Call AntFun(thet_d,phi_d) ! Jones; gives _0 & _1 voltages when contracted with (t,p) polarized field
   Det(:)=J_0p(:)*J_1t(:) - J_1p(:)*J_0t(:)
   Ji_t1(:)=  J_0p(:)/Det(:)
   Ji_t0(:)=-J_1p(:)/Det(:)
   Ji_p1(:)=-J_0t(:)/Det(:)
   Ji_p0(:)=  J_1t(:)/Det(:)
   Return
End Subroutine AntFun_Inv
! ===========================
Subroutine AntFieParGen()
   use constants, only : dp,pi,ci
   use AntFunCconst
   Implicit none
   Real(dp) :: Cheb(m_dim), Cheb_s(m_dim), Lgdr(j_dim)
   Real(dp) :: Freq0, Thet0, Phi0, Vp_r, Vp_i
   Real(dp) :: Freq, Thet, Phi, Vt_r, Vt_i
   !Complex(dp) :: J_0p(Freq_min:Freq_max),J_0t(Freq_min:Freq_max),J_1p(Freq_min:Freq_max),J_1t(Freq_min:Freq_max)
   Complex(dp) :: PhIn_p(m_dim,20), PhIn_t(m_dim,20)
   Real(dp) :: th, x, dx, dthet, S1p25, c1p25
   Complex(dp) :: Vp, Vt
   Integer :: nxx, it, nthet, i_freq, j
   Character(len=150) :: AntFunFile
   !Complex(dp) :: J_0p(Freq_min:Freq_max),J_0t(Freq_min:Freq_max),J_1p(Freq_min:Freq_max),J_1t(Freq_min:Freq_max) ! only for testing
   !Complex(dp) :: Ji_p0(Freq_min:Freq_max),Ji_t0(Freq_min:Freq_max),Ji_p1(Freq_min:Freq_max),Ji_t1(Freq_min:Freq_max)
   !
   !WRITE (*,*) "LIBRARY=",TRIM(ldata)
   CALL get_environment_variable("AntennaFun", AntFunFile)
   WRITE (2,*) "Using environment variable: AntennaFun; Antenna Function from files:",' "',TRIM(AntFunFile)//'LBA_Vout_*.txt"'
   If(AntTst) Then
      open(Unit=80,Status='unknown',Action='write',File='V_t-nu.dat') !fix th=0, phi=0, nu-dependence
      open(Unit=81,Status='unknown',Action='write',File='V_p-nu.dat') !fix th=0, phi=90, nu-dependence
      open(Unit=82,Status='unknown',Action='write',File='V_t60-Th.dat') !th-dep, phi=0, fix nu=60
      open(Unit=83,Status='unknown',Action='write',File='V_p60-Th.dat') !th-dep, phi=90, fix nu=60
   Endif
   S1p25=sin(1.25*pi/180.) ! needed for average angle first and last half bin
   c1p25=cos(1.25*pi/180.)
   Open(Unit=10,Status='old',Action='read',File=TRIM(AntFunFile)//'LBA_Vout_phi.txt')
   Open(Unit=11,Status='old',Action='read',File=TRIM(AntFunFile)//'LBA_Vout_theta.txt')
   read(10,*)
   read(11,*)  ! skip the first comment line
   PhIn_p(:,:)=0.
   PhIn_t(:,:)=0.
   Freq=0.
   nxx=0
   write(2,*) 'causal ant fie'
   Do while(nxx.eq.0)
      i_freq=nint(Freq)
      If(i_freq .lt. Freq_min) goto 1  ! this is the case on the first pass
      If(i_freq .gt. Freq_max) exit
      !write(2,*) Freq, Thet, Phi, Vp_r, Vp_i
      it=nint(1. + thet/5.)
      Call Cheb_odd(Phi*pi/180.,Cheb)
      Call Cheb_odd((Phi+90.)*pi/180.,Cheb_s)
      If(Phi .lt. 359.) then ! integrate over \phi since \phi=0 is same as 360
         PhIn_p(:,it)=PhIn_p(:,it)+Cheb(:)*Vp/18.  ! 1/18=2*10/360
         PhIn_t(:,it)=PhIn_t(:,it)+Cheb_s(:)*Vt/18.  ! 1/18=2*10/360
      else ! Phi=360, print and initialize integral
         ! Write(2,"(F5.1,F7.1,4F10.4,4F10.4)") Freq, Thet,  PhIn_p(:,it), PhIn_t(:,it)
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
                  Ant_p(i_freq,j,1:m_dim)=Ant_p(i_freq,j,1:m_dim) + PhIn_p(1:m_dim,it+1)*Lgdr(j)*(4.*j-1.)*dx ! norm: (2n+1)/18 where n=2j-1; 0.5 since interval at th=0 is only half width
                  Ant_t(i_freq,j,1:m_dim)=Ant_t(i_freq,j,1:m_dim) + PhIn_t(1:m_dim,it+1)*Lgdr(j)*(4.*j-1.)*dx ! norm: (2n+1)/18 where n=2j-1; 0.5 since interval at th=0 is only half width
               Enddo
            Enddo
            ! write(2,"(F5.1,A,4F10.4,4F10.4)") freq,' Phi(re,im):=',Ant_p(i_freq,:,1)
            ! write(2,"(F5.1,A,4F10.4,4F10.4)") freq,' Thet(re,im):=',Ant_t(i_freq,:,1)
            PhIn_p(:,:)=0.
            PhIn_t(:,:)=0.
            !stop
         Endif ! (Thet.gt.89)
      Endif ! (Phi .lt. 359.)
   1  continue
      read(10,*,iostat=nxx) Freq0, Thet0, Phi0, Vp_r, Vp_i
      read(11,*,iostat=nxx) Freq, Thet, Phi, Vt_r, Vt_i
      Vp=cmplx(Vp_r,-Vp_i) ! minus is because of a different convention between Python and FFTPack
      Vt=cmplx(Vt_r,-Vt_i)
      If(AntTst) Then
         !open(Unit=80,Status='unknown',Action='write',File='V_t-nu.dat') !fix th=0, phi=0, nu-dependence
         If((NINT(Thet) .eq. 0) .and. (NINT(Phi) .eq. 0)) then
            write(80,"(F6.1, 2(2f10.4,f10.2,f8.1))") Freq, Vt_r, Vt_i, Vt_r*Vt_r + Vt_i*Vt_i, atan2(Vt_r, Vt_i)*180/pi
            If(NINT(Freq).eq.30) write(2,*) '@30,vt=',vt
         Endif
         !open(Unit=81,Status='unknown',Action='write',File='V_p-nu.dat') !fix th=0, phi=90, nu-dependence
         If((NINT(Thet) .eq. 0) .and. (NINT(Phi) .eq. 90)) then
            write(81,"(F6.1, 2(2f10.4,f10.2,f8.1))") Freq, Vp_r, Vp_i, Vp_r*Vp_r + Vp_i*Vp_i, atan2(Vp_r, Vp_i)*180/pi
         Endif
         !open(Unit=82,Status='unknown',Action='write',File='V_t60-Th.dat') !th-dep, phi=0, fix nu=60
         If((NINT(freq) .eq. 60) .and. (NINT(Phi) .eq. 0)) then
            write(82,"(F6.1, 2(2f10.4,f10.2,f8.1))") Thet, Vt_r, Vt_i, Vt_r*Vt_r + Vt_i*Vt_i, atan2(Vt_r, Vt_i)*180/pi
         Endif
         !open(Unit=83,Status='unknown',Action='write',File='V_p60-Th.dat') !th-dep, phi=90, fix nu=60
         If(NINT(freq) .eq. 60 .and. NINT(Phi) .eq. 90) then
            write(83,"(F6.1, 2(2f10.4,f10.2,f8.1))") Thet, Vp_r, Vp_i, Vp_r*Vp_r + Vp_i*Vp_i, atan2(Vp_r, Vp_i)*180/pi
         Endif
      EndIf
      If( nint(Freq0*1000.+ Thet0*100.+ Phi0) .ne. nint(Freq*1000.+ Thet*100.+ Phi)) then
         write(2,*) 'mismatch',nint(Freq0*1000.+ Thet0*100.+ Phi0),nint(Freq*1000.+ Thet*100.+ Phi)
         stop 'AntFieParGen'
      endif
   Enddo
   !
   Close(unit=10)
   Close(unit=11)
   !Ant_p=conjg(Ant_p)
   !Ant_t=conjg(Ant_t)
   !   !
   !   !
   !   Open(Unit=10,Status='old',Action='read',File=TRIM(AntFunFile)//'LBA_Vout_phi.txt')
   !   Open(Unit=11,Status='old',Action='read',File=TRIM(AntFunFile)//'LBA_Vout_theta.txt')
   !   read(10,*)
   !   read(11,*)
   !   Freq=0.
   !   nxx=0
   !   Do while(nxx.eq.0)
   !      !read(10,*,iostat=nxx) Freq0, Thet0, Phi0, Vp_r, Vp_i
   !      read(11,*,iostat=nxx) Freq, Thet, Phi, Vt_r, Vt_i
   !      i_freq=nint(Freq)
   !      If(i_freq .lt. Freq_min) cycle
   !      If(i_freq .gt. Freq_max) exit
   !      Call AntFun(thet,phi-45.) ! sets ,J_0p,J_0t,J_1p,J_1t ; Jones; gives _0 & _1 voltages when contracted with (t,p) polarized field
   !      write(2,*) Freq, Thet, Phi, Vt_r, Vt_i, J_0t(i_freq) !, J_1t(i_freq)
   !   Enddo
   !   stop
   !   !
   !!
   !Write(2,*) 'Ant_p(60,1,1)',Ant_p(60,1,1) !,' ********** Modified by a factor 1/2, not sure why !!!!!!!!!!!!!'
   !Write(2,*) 'Ant_t(60,1,1)',Ant_t(60,1,1)
   Write(2,*) 'Ant_t(30,1,1)',Ant_t(30,1,1)
   If(AntTst) Then
      Close(Unit=80)
      Close(Unit=81)
      Close(Unit=82)
      Close(Unit=83)
      open(Unit=80,Status='unknown',Action='write',File='aV_t-nu.dat') !fix th=0, phi=0, nu-dependence
      !open(Unit=81,Status='unknown',Action='write',File='aV_p-nu.dat') !fix th=0, phi=90, nu-dependence
      open(Unit=82,Status='unknown',Action='write',File='aV_t60-Th.dat') !th-dep, phi=0, fix nu=60
      !open(Unit=83,Status='unknown',Action='write',File='aV_p60-Th.dat') !th-dep, phi=90, fix nu=60
      phi=-45.
      thet=0. ! in degrees
      Call AntFun(thet,phi) ! sets ,J_0p,J_0t,J_1p,J_1t ; Jones; gives _0 & _1 voltages when contracted with (t,p) polarized field
      Call AntFun_Inv(thet,Phi) ! sets ,Ji_p0,Ji_t0,Ji_p1,Ji_t1; Inverse Jones; gives _p & _t polarized fields when contracted with (0.1) voltages
      write(2,*) 'J0t@30',J_0t(30)
      Do i_freq=Freq_min,Freq_max
         Vp_r=-Real(J_1p(i_freq)) ; Vp_i=-Imag(J_1p(i_freq))
         Vt_r=Real(J_0t(i_freq)) ; Vt_i=Imag(J_0t(i_freq))
         write(80,"(I4, 2(2f10.4,f10.2,f8.1))") i_freq, Vt_r, Vt_i, Vt_r*Vt_r + Vt_i*Vt_i, atan2(Vt_r, Vt_i)*180/pi &
            ,Vp_r, Vp_i, Vp_r*Vp_r + Vp_i*Vp_i, atan2(Vp_r, Vp_i)*180/pi
         vp=J_0p(i_freq)*Ji_p0(i_freq)+J_0t(i_freq)*Ji_t0(i_freq)
         vt=J_0p(i_freq)*Ji_p0(i_freq)+J_1p(i_freq)*Ji_p1(i_freq)
         write(2,*) 'orthonormality', vp,vt,J_0p(i_freq)*Ji_p1(i_freq)+J_0t(i_freq)*Ji_t1(i_freq)
      Enddo
      i_freq=60
      Do it=0,89,1
         thet=it
         Call AntFun(thet,phi) ! sets ,J_0p,J_0t,J_1p,J_1t ; Jones; gives _0 & _1 voltages when contracted with (t,p) polarized field
         Vp_r=-Real(J_1p(i_freq)) ; Vp_i=-Imag(J_1p(i_freq))
         Vt_r=Real(J_0t(i_freq)) ; Vt_i=Imag(J_0t(i_freq))
         write(82,"(I4, 2(2f10.4,f10.2,f8.1) )") it, Vt_r, Vt_i, Vt_r*Vt_r + Vt_i*Vt_i, atan2(Vt_r, Vt_i)*180/pi &
            ,Vp_r, Vp_i, Vp_r*Vp_r + Vp_i*Vp_i, atan2(Vp_r, Vp_i)*180/pi
      Enddo
      Close(Unit=80)
      Close(Unit=82)
   !!
   Endif
   !
   Call AntFun(0.d0,0.d0) ! sets ,Ji_p0,Ji_t0,Ji_p1,Ji_t1; Inverse Jones; gives _p & _t polarized fields when contracted with (0.1) voltages
   Do i_freq=Freq_min,Freq_max   ! Increment frequency spectrum with this antenna
      Gain(i_freq)=abs(J_0p(i_freq))**2+ abs(J_0t(i_freq))**2+ abs(J_1p(i_freq))**2+ abs(J_1t(i_freq))**2
      Gain(i_freq)=sqrt(Gain(i_freq))/(Freq_max-Freq_min)
      !write(2,*) i_freq, Gain(i_freq)
   Enddo
   Return
End Subroutine AntFieParGen
!=====================
