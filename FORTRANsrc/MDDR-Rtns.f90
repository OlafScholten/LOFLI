!-----------------------------------------------
Subroutine MDDSetup()
!
!  Input:
!     - IntfNuDim giving IntfDim
!     C_Tms, CenLoc(1:3)
!     - Central location and time of composite source
!  Output:
!     Weight(j_IntFer)
!     - Weight factors for antennas and may depend on polarisation, W[ant] == Weight(1:Nr_IntFer),
!        same as calculated in first part of "EISetupSpec(,,,)" in 'EIOption.f90'
!     CTime_p(1:IntfDim,j_IntFer), CTime_t(1:IntfDim,j_IntFer)
!     - E[t_c,ant,k] where [k]=[theta,phi] polarization directions and time shifted to the center of the source region,
!     Vec_p(j_IntFer,i), Vec_t(j_IntFer,i)
!     - theta[ant,i]  & phi[ant,i] i=Cartesian coordinates
!     AntSourceD(j_IntFer)
!     - Distance to center source region
!
!--------------------------------------------
   use constants, only : dp, pi, ci, Sample, c_mps
   use DataConstants, only : Time_Dim
   use Chunk_AntInfo, only : CTime_spectr, Ant_Stations, Ant_IDs, Ant_nr, Ant_pos, Ant_RawSourceDist
   use Chunk_AntInfo, only : NormOdd, NormEven
   Use Interferom_Pars, only : i_chunk, IntFer_ant, Nr_IntFerCh
   Use Interferom_Pars, only : StartTime_ms
   Use MDD_Pars, only : MDDtDim, MDDnuDim, Weight, RTime_p, RTime_t, Vec_p, Vec_t, Vec_l, AntSourceD  ! Output
   Use MDD_Pars, only : DDChiSQ
   Use MDD_Pars, only : CR_Tms, C_Tms, C_Loc  ! Input
   use AntFunCconst, only : Freq_min, Freq_max,Ji_p0,Ji_t0,Ji_p1,Ji_t1, Gain ! J_0p,J_0t,J_1p,J_1t,
   use FFT, only :  RFTransform_su, RFTransform_CF, RFTransform_CF2CT, RFTransform_CF2RT
   Implicit none
   Integer :: Nr_IntFer
   integer :: i_ant, j_IntFer, i, j, i_freq, i_nu, IntfDim, inu1, inu2, N_Offset
   Real(dp) :: Ras(1:3), WNrm
   Real(dp) :: HorDist, Thet_r ,Phi_r, dfreq, D, W, dnu, nu
   Real(dp) :: thet_d, Phi_d, Rdist, SamplOff, dt, AntennaNorm=1.! e-4
   Real(dp), allocatable :: Rtime(:)
   Complex(dp), allocatable :: Cnu0(:), Cnu1(:), Cnu_p(:), Cnu_t(:)
   Complex(dp) :: Q, phase, dphase
   complex(dp), parameter :: ipi=ci*pi
   !
   Nr_IntFer=Nr_IntFerCh(i_chunk)
   dnu=100./MDDnuDim    ! Top frequncy is 100MHz, ! [MHz] Jones matrix is stored on 1MHz grid
   inu1=Int(Freq_min/dnu)+1
   inu2=Int(Freq_max/dnu)-1
   !---- global:
   allocate (RTime_p(1:MDDtDim,1:Nr_IntFer))
   allocate (RTime_t(1:MDDtDim,1:Nr_IntFer))
   allocate (Vec_p(1:3,1:Nr_IntFer), Vec_t(1:3,1:Nr_IntFer), Vec_l(1:3,1:Nr_IntFer))
   allocate (Weight(1:Nr_IntFer),  AntSourceD(1:Nr_IntFer))
   allocate (DDChiSQ(1:Nr_IntFer,1:2) )
   !----- local:
   allocate (RTime(1:MDDtDim) )
   allocate ( Cnu0(0:MDDnuDim), Cnu1(0:MDDnuDim), Cnu_p(0:MDDnuDim), Cnu_t(0:MDDnuDim) )
   !
   Call RFTransform_su(MDDtDim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ! Call DAssignFFT()                   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !
   write(2,*) 'Nr_IntFer=',Nr_IntFer
   Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
      i_ant=IntFer_ant(j_IntFer,i_chunk)  ! Refers to even antenna
      !write(2,*) i_ant, j_IntFer, Ant_pos(:,i_ant,i_chunk)
      Ras(1)=(C_Loc(1)-Ant_pos(1,i_ant,i_chunk))/1000.
      Ras(2)=(C_Loc(2)-Ant_pos(2,i_ant,i_chunk))/1000.
      Ras(3)=(C_Loc(3)-Ant_pos(3,i_ant,i_chunk))/1000.
      HorDist= Ras(1)*Ras(1) + Ras(2)*Ras(2)  ! =HYPOT(X,Y)
      AntSourceD(j_IntFer)=sqrt(HorDist + Ras(3)*Ras(3))
      HorDist=sqrt( HorDist ) ! =HYPOT(X,Y)
      Thet_r=atan(HorDist,Ras(3))  ! Zenith angle
      Phi_r=atan2( Ras(2) , Ras(1) )  ! \phi=0 = north
      Thet_d=Thet_r*180/pi
      Phi_d=Phi_r*180/pi
      If(j_IntFer.eq.1) Then
         WNrm = AntSourceD(1)*AntSourceD(1)/(Ras(3)*Ras(3))  ! to normalize at unity=(1/Nr_IntFer) for the reference antenna and A=(1,1,0)
         !WNrm = AntSourceD(1)*AntSourceD(1)/(Ras(3)*Ras(3)*Nr_IntFer)  ! to normalize at unity=(1/Nr_IntFer) for the reference antenna and A=(1,1,0)
      EndIf
      Weight(j_IntFer)= Ras(3)*Ras(3)* WNrm /(AntSourceD(j_IntFer)*AntSourceD(j_IntFer)) ! changed Nov 2022 to make it more similar to 1/noise power
      If((AntSourceD(j_IntFer)/AntSourceD(1)).lt.1. ) Then  ! The factor between E-field and signal is approx square of this
         Weight(j_IntFer) = Weight(1)
      EndIf  ! this gives a much smoother and narrower interference max
      D=AntSourceD(j_IntFer)*AntSourceD(j_IntFer)
      !
      Vec_p(1,j_IntFer)=sin(Phi_r)             ; Vec_p(2,j_IntFer)=-cos(Phi_r)          ; Vec_p(3,j_IntFer)=0.
      Vec_t(1,j_IntFer)=-cos(Thet_r)*Vec_p(2,j_IntFer)
      Vec_t(2,j_IntFer)=cos(Thet_r)*Vec_p(1,j_IntFer)
      Vec_t(3,j_IntFer)=-sin(Thet_r)
      Vec_l(1:3,j_IntFer)=Ras(1:3)/AntSourceD(j_IntFer)
      !
      i_ant=IntFer_ant(j_IntFer,i_chunk)
      Call RelDist(C_Loc(1),Ant_pos(1,i_ant,i_chunk),RDist)
      !
      SamplOff=((CR_Tms-StartTime_ms)/1000.d0)/sample + Rdist - Ant_RawSourceDist(i_ant,i_chunk) - MDDnuDim
      N_Offset=NINT(SamplOff)
      dt=SamplOff-N_Offset
      ! check dt_AntPix in range
      !write(2,*) 'SamplOff:', (CR_Tms-StartTime_ms), Rdist - Ant_RawSourceDist(i_ant,i_chunk), MDDnuDim
      If(N_Offset.lt.2) then
         write(2,*) 'Warning, time-shift=',SamplOff,'[samples] below range, i_ant=',i_ant,(CR_Tms-StartTime_ms),'[ms]',&
            (Rdist - Ant_RawSourceDist(i_ant,i_chunk)),'[s]'
         Return
      Endif
      If((N_Offset+2*MDDnuDim).gt.(Time_dim-2)) then
         write(2,*) 'Warning, time-shift=',SamplOff,'[samples] above range, i_ant=',i_ant,(CR_Tms-StartTime_ms),'[ms]',&
            (Rdist - Ant_RawSourceDist(i_ant,i_chunk)),'[s]',N_Offset+2*MDDnuDim
         Return
      Endif
      !
      W=AntennaNorm*Weight(j_IntFer)*NormEven
      RTime(:)=REALPART(CTime_spectr(N_Offset+1:N_Offset+2*MDDnuDim,i_ant,i_chunk))*W
      Call RFTransform_CF(RTime,Cnu0(0))
      W=AntennaNorm*Weight(j_IntFer)*NormOdd
      RTime(:)=REALPART(CTime_spectr(N_Offset+1:N_Offset+2*MDDnuDim,i_ant+1,i_chunk))*W
      Call RFTransform_CF(RTime,Cnu1(0))
      !
      Call AntFun_Inv(thet_d,Phi_d) ! sets ,Ji_p0,Ji_t0,Ji_p1,Ji_t1; Inverse Jones; gives _p & _t polarized fields when contracted with (0.1) voltages
      dphase = exp(-ipi*dt/MDDnuDim)
      Phase =exp(-ipi*dt*inu1/MDDnuDim)
      Cnu_p(:)=0.d0
      Cnu_t(:)=0.d0
      Do i_nu=inu1,inu2   ! Increment frequency spectrum with this antenna
         nu=i_nu*dnu
         i_freq=Int(nu)  ! [MHz] Jones matrix is stored on 1MHz grid
         dfreq=nu-i_freq
         Q = phase * ((1.-dfreq)*Gain(i_freq) + dfreq*Gain(i_freq+1)) ! Gain = amplitude gain
         Cnu0(i_nu)=Cnu0(i_nu) * Q
         Cnu1(i_nu)=Cnu1(i_nu) * Q
         Cnu_p(i_nu)=((1.-dfreq)*Ji_p0(i_freq) + dfreq*Ji_p0(i_freq+1)) *Cnu0(i_nu) + &
               ((1.-dfreq)*Ji_p1(i_freq) + dfreq*Ji_p1(i_freq+1)) *Cnu1(i_nu)
         Cnu_t(i_nu)=((1.-dfreq)*Ji_t0(i_freq) + dfreq*Ji_t0(i_freq+1)) *Cnu0(i_nu) + &
               ((1.-dfreq)*Ji_t1(i_freq) + dfreq*Ji_t1(i_freq+1)) *Cnu1(i_nu)
         phase=phase*dphase
      EndDo
      !
      ! convert to time
      Call RFTransform_CF2RT(Cnu_p(0),RTime_p(1,j_IntFer) )
      Call RFTransform_CF2RT(Cnu_t(0),RTime_t(1,j_IntFer) )
      !
   Enddo         ! loop antennas for Interferometry
   !
   !write(2,*) 'CMCnu(Nr_IntFer/2,1:2,MDDnuDim/2)', MDDnuDim/2,  IntFer_ant(Nr_IntFer/2),  CMCnu(Nr_IntFer/2,1,MDDnuDim/2) &
   !   , IntFer_ant(Nr_IntFer/2+1),CMCnu(Nr_IntFer/2+1,1,MDDnuDim/2)
   !Flush(unit=2)
   !----- local:
   Deallocate (RTime)
   Deallocate ( Cnu0, Cnu1, Cnu_p, Cnu_t )
   !
   Return
End Subroutine MDDSetup
!-----------------------------------------------
Subroutine MDDFitSU()
   use constants, only : dp, sample, c_mps
   Use MDD_Pars, only : NSrc_max, IRFW_s, T_Range  ! Input
   Use MDD_Pars, only : NSrc, dCenT, dCenloc  ! Output
   Implicit none
   Real(dp) :: dt, dd, D
   Integer :: i
   !
   dt=IRFW_s  ! [samples]
   dd=IRFW_s*c_mps*sample  ! [samples] * [m/s] * [s/sample] =[m]
   dCenT(:)=0.
   dCenloc(1:3,:)=0.
   dCenT(2)=0. ;  dCenloc(1,2)=+dd
   dCenT(3)=0. ;  dCenloc(1,3)=-dd
   dCenT(4)=0. ;  dCenloc(2,4)=+dd
   dCenT(5)=0. ;  dCenloc(2,5)=-dd
   dCenT(6)=0. ;  dCenloc(3,6)=+dd
   dCenT(7)=0. ;  dCenloc(3,7)=-dd
   dCenT(8)=+dt
   dCenT(9)=-dt
   If(T_Range/2. .lt. sqrt(2.)*dt) Then
      Nsrc=9
      Return
   EndIf
   dCenT(10)=+dt ;  dCenloc(1,10)=+dd
   dCenT(11)=+dt ;  dCenloc(1,11)=-dd
   dCenT(12)=-dt ;  dCenloc(1,12)=+dd
   dCenT(13)=-dt ;  dCenloc(1,13)=-dd
   dCenT(14)=+dt ;  dCenloc(2,14)=dd
   dCenT(15)=+dt ;  dCenloc(2,15)=-dd
   dCenT(16)=-dt ;  dCenloc(2,16)=dd
   dCenT(17)=-dt ;  dCenloc(2,17)=-dd
   Nsrc=17
   NSrc=2
   Return
End Subroutine MDDFitSU
!-----------------------------------------------
Subroutine ImpResp_ns()
!
!     - G(t) - upsampled pulse-response of system due to the frequency dependence of the gain.
!        Note that we compare gain-multiplied E-fields.
!        t=0ns corresponds to element G(G_dim/2)
!-----------------------------------------------
   use constants, only : dp
   Use MDD_Pars, only : G_dim, G_ns, IRFW_s  ! Output
   use AntFunCconst, only : Freq_min, Freq_max, Gain ! J_0p,J_0t,J_1p,J_1t,
   use FFT, only : RFTransform_su, DAssignFFT, RFTransform_CF, RFTransform_CF2CT
   Implicit none
   !Integer, intent(IN) :: G_dim
   !Real(dp), intent(OUT) :: G_ns(G_dim)
   Complex(dp) :: Gnu(0:G_dim/2),CG_ns(1:G_dim)
   Integer :: GnuDim, i_nu, inu1, inu2, i_freq,i
   Real(dp) :: dnu, Q, nu, dfreq, ScaleFactor
   !
   Call RFTransform_su(G_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   G_ns(:)=0.d0
   G_ns(G_dim/2)=1.d3
   Call RFTransform_CF(G_ns,Gnu(0))
   GnuDim=G_Dim/2
   dnu=500./GnuDim  ! Top frequncy is taken here at 500MHz, i.e. Nyquist sampled at 1ns
   inu1=Int(Freq_min/dnu)+1  ! [MHz] Jones matrix is stored on 1MHz grid
   inu2=Int(Freq_max/dnu)-1
   Gnu(0:inu1)=0.   ! Increment frequency spectrum with this antenna
   Do i_nu=inu1,inu2  ! Increment frequency spectrum with this antenna
      nu=i_nu*dnu
      i_freq=Int(nu)  ! [MHz] Jones matrix is stored on 1MHz grid
      dfreq=nu-i_freq
      Q = ((1.-dfreq)*Gain(i_freq) + dfreq*Gain(i_freq+1)) ! Gain = amplitude gain
      Gnu(i_nu)=Gnu(i_nu) * Q
   EndDo
   Gnu(inu2+1:GnuDim )=0.   ! Increment frequency spectrum with this antenna
   Call RFTransform_CF2CT(Gnu(0),CG_ns(1) )
   Q=Abs(CG_ns(G_dim/2))
   Do i=G_dim/2,G_dim
      nu=Abs(CG_ns(i))
      If(nu/Q .lt. 0.5d0) exit
   EndDo
   IRFW_s=(2*i-G_dim)/5.
   write(2,*) 'Impulse Response FWHM(amplitude)=',IRFW_s,'[samples], peak=', Q
   G_ns(:)=Real(CG_ns(:))
   !write(2,"(10F10.5)") G_ns(:)
   !
   call DAssignFFT()                   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   Return
End Subroutine ImpResp_ns
!-----------------------------------------------
Subroutine MDD_write()
!
!-----------------------------------------------
   use constants, only : dp
   Use Interferom_Pars, only : IntFer_ant, i_chunk, Nr_IntFerCh
   Use MDD_Pars, only : IRFW_s
   Use MDD_Pars, only : NSrc, dCenT, dCenloc, T_Range
   Use MDD_Pars, only : DDChiSQ
   Use MDD_Pars, only : DelDip
   Use MDD_Pars, only : C_Tms, C_Loc
   Implicit none
   Integer :: i, m
   write(2,"(A,F10.5,A,2(F7.3,','),F6.3,A,I4,A,I3,A )") 'C-time=',C_Tms,'[ms], C-loc (N,E,h)=(',C_Loc(1:)/1000., &
      ') [km], Fit window=',T_Range,' samples with',NSrc, ' Delta Dipole sources'
   Write(2,"('  #  dt[smpl]',T15,'dDeltDip (N,E,h) [m]', T43, 'dipole moments (N,E,h) [?]')")
   Do m=1,NSrc
      Write(2,"(I3,F8.2, 3F8.2, 3F12.4)") m,dCenT(m), dCenloc(1:3,m), DelDip(1:3,m)
   EndDo
   Return
End Subroutine MDD_write
!-------------------------------------
Subroutine MDDClose()
   Use MDD_Pars, only : Weight, RTime_p, RTime_t, Vec_p, Vec_t, Vec_l, AntSourceD
   Use MDD_Pars, only : DDChiSQ
   use FFT, only : DAssignFFT
   !
   call DAssignFFT()                   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   Deallocate (RTime_p, RTime_t, Vec_p, Vec_t, Vec_l, Weight,  AntSourceD, DDChiSQ )
   !
   Return
End Subroutine MDDClose
