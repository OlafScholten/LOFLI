!-----------------------------------------------
Subroutine EffAntSr()
!  Determine effective opening angle of a LOFAR antenna for galactic background
   use constants, only : dp, pi, ci, Sample
   use AntFunCconst, only : Freq_min, Freq_max,J_0p,J_0t,J_1p,J_1t !, Gain ! J_0p,J_0t,J_1p,J_1t,
   Implicit none
   Integer :: i_ct, i_ph
   Integer, parameter :: N_ct=21, N_ph=51, Freq=60
   Real(dp) :: Sr_f, ct, d_ct, th_d, d_ph, ph_d, Gain_0, Gain_sr
   Sr_f=0.
   th_d=0.
   ph_d=0.
   Call AntFun(th_d,Ph_d) !
   Gain_0=abs(J_0p(Freq))**2 + abs(J_0t(Freq))**2
   gain_sr=0.
   d_ct=1./N_ct
   d_ph=2*pi/N_ph
   Do i_ct=1, N_ct
      Ct=d_ct*(i_ct-0.5)
      th_d=acos(ct)*180/pi
      Do i_ph=1, N_ph
         ph_d=d_ph*(i_ph-0.5)
         Sr_f=Sr_f+ d_ct*d_ph
         Call AntFun(th_d,Ph_d) !
         gain_sr=gain_sr + (abs(J_0p(Freq))**2 + abs(J_0t(Freq))**2)*d_ct*d_ph
      EndDo
   EndDo
   write(2,*) 'EffAntSr:',Sr_f/pi, Freq, Gain_0, gain_sr/gain_0
End Subroutine EffAntSr
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
   use constants, only : dp, pi, ci, Sample
   use DataConstants, only : Time_Dim
   use Chunk_AntInfo, only : CTime_spectr, Ant_Stations, Ant_IDs, Ant_nr, Ant_pos, Ant_RawSourceDist
   use Chunk_AntInfo, only : NormOdd, NormEven
   Use Interferom_Pars, only : i_chunk, IntFer_ant, Nr_IntFerCh
   Use Interferom_Pars, only : StartTime_ms
   Use MDD_Pars, only : MDDtDim, MDDnuDim, Weight, SQWght, CTime_p, CTime_t, Vec_p, Vec_t, Vec_l, AntSourceD  ! Output
   Use MDD_Pars, only : DDChiSQ  ! allocated only
   Use MDD_Pars, only : CR_Tms, C_Tms, C_Loc, DepleteCoreWeigth, DepleteDistFrac  ! Input
   use AntFunCconst, only : Freq_min, Freq_max,Ji_p0,Ji_t0,Ji_p1,Ji_t1, Gain ! J_0p,J_0t,J_1p,J_1t,
   use FFT, only :  RFTransform_su, RFTransform_CF, RFTransform_CF2CT
   use StationMnemonics, only : Station_ID2Mnem ! Station_Mnem2ID,
   Implicit none
   Integer :: Nr_IntFer
   integer :: i_ant, j_IntFer, i, j, i_freq, i_nu, IntfDim, inu1, inu2, N_Offset
   Real(dp) :: Ras(1:3), WNrm
   Real(dp) :: HorDist, Thet_r ,Phi_r, dfreq, W, dnu, nu ! , D
   Real(dp) :: thet_d, Phi_d, Rdist, SamplOff, dt, AntennaNorm=1.! e-4
   Real(dp), allocatable :: Rtime(:)
   Complex(dp), allocatable :: Cnu0(:), Cnu1(:), Cnu_p(:), Cnu_t(:)
   Complex(dp) :: Q, phase, dphase
   complex(dp), parameter :: ipi=ci*pi
   Character(len=5) :: Station_MN
   !
   Nr_IntFer=Nr_IntFerCh(i_chunk)
   dnu=100./MDDnuDim    ! Top frequncy is 100MHz, ! [MHz] Jones matrix is stored on 1MHz grid
   inu1=Int(Freq_min/dnu)+1
   inu2=Int(Freq_max/dnu)-1
   !---- global:
   allocate (CTime_p(1:MDDtDim,1:Nr_IntFer))
   allocate (CTime_t(1:MDDtDim,1:Nr_IntFer))
   allocate (Vec_p(1:3,1:Nr_IntFer), Vec_t(1:3,1:Nr_IntFer), Vec_l(1:3,1:Nr_IntFer))
   allocate (Weight(1:Nr_IntFer), SQWght(1:Nr_IntFer),  AntSourceD(1:Nr_IntFer))
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
      If( AntSourceD(j_IntFer) .lt. DepleteDistFrac*AntSourceD(1)) Then  ! Keep weights constant for antennas that are closer to the source than this fraction of the distance to the core
         Weight(j_IntFer) = Weight(1)
      EndIf  ! this gives a much smoother and narrower interference max
      !
      ! Decrease weight for core stations
      If(DepleteCoreWeigth.lt.1.) Then
         Call Station_ID2Mnem(Ant_Stations(i_ant,1),Station_MN)
         If(Station_MN(1:2).eq.'CS') Then
            Weight(j_IntFer)=Weight(j_IntFer)*DepleteCoreWeigth
         EndIf ! (Station_MN(1:2).eq.'CS')
      EndIf ! (DepleteCoreWeigth.lt.1.)
      !D=AntSourceD(j_IntFer)*AntSourceD(j_IntFer)
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
      Call RFTransform_CF2CT(Cnu_p(0),CTime_p(1,j_IntFer) )
      Call RFTransform_CF2CT(Cnu_t(0),CTime_t(1,j_IntFer) )
      !
   Enddo         ! loop antennas for Interferometry
   SQWght(:)=sqrt( Weight(:) )
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
Subroutine MDD_IniManySrcs()
   use constants, only : dp, sample, c_mps
   Use Interferom_Pars, only : i_chunk, Nr_IntFerCh
   Use MDD_Pars, only : NSrc_max, IRFW_s, T_Range, DelDip, Weight, SQWght, Vec_l, Trace_DeadEnds, GridScale  ! Input
   Use MDD_Pars, only : NSrc, dCenT, dCenloc, SpaceGrid  ! Output
   Use MDD_Pars, only : Bias_inv_base, DDXcorr_depl, Trace_DeadEnds, RealDip_logic  ! input
   Implicit none
   Real(dp) :: D, del_t, Ave_l(1:3)
   Real(dp) :: Gr(0:3,0:3), Sig(1:3), B, Length(0:3)
   Integer :: j_IntFer, m, Nr_IntFer, i
   Integer :: i_gr(0:3), Ngrd(0:3)
   Integer, parameter :: Ngrd_max=20
   !Integer, parameter :: Tailend=25
   Real(dp), parameter :: dD_max=300.d0 ! [m]  ! max distance from center
   !
   ! Determine relative size space grid from
   !   del_t= -dCenT(m) +SUM(Vec_l(1:3,j_IntFer)*dCenloc(1:3,m))/(c_mps*sample)   ! in [samples]
   Nr_IntFer=Nr_IntFerCh(i_chunk)
   Gr(0:3,0:3)=0.
   Ave_l(:)=0.
   Do j_IntFer=1,Nr_IntFer
      Ave_l(:)= Ave_l(:) + Vec_l(:,j_IntFer)/Nr_IntFer
   EndDo ! j_IntFer=1,Nr_IntFer
   D=sum(Ave_l(:)*Ave_l(:)) ;  Gr(1,1:3)=Ave_l(:)/sqrt(D)  !  will probably be pointing in the radial direction similar change as dt
   write(2,*) 'sqrtD', sqrt(D)
   Gr(2,1)=Ave_l(2) ; Gr(2,2)=-Ave_l(1) ; Gr(2,3)=0.    ! Ave_p . Ave_l =0.
   D=sum(Gr(2,1:3)*Gr(2,1:3)) ;  Gr(2,1:3)=Gr(2,1:3)/sqrt(D)
   Gr(3,1)=Gr(2,2)*Ave_l(3) ; Gr(3,2)=-Gr(2,1)*Ave_l(3) ; Gr(3,3)=Gr(2,1)*Ave_l(2)-Gr(2,2)*Ave_l(1)  ! Ave_t =Ave_p x Ave_l
   D=sum(Gr(3,1:3)*Gr(3,1:3)) ;  Gr(3,1:3)=Gr(3,1:3)/sqrt(D)
   Sig(1:3)=0.
   B=0.
   Do j_IntFer=1,Nr_IntFer  ! determine sigma for the 3 directions
      Do i=1,3
         Sig(i)= Sig(i) + (SUM(Vec_l(1:3,j_IntFer)*Gr(i,1:3)))**2
      EndDo
      B=B+SUM(Vec_l(1:3,j_IntFer)*Gr(1,1:3))/Nr_IntFer
   EndDo ! j_IntFer=1,Nr_IntFer
   Sig(1)=sqrt(Sig(1)/Nr_IntFer-B*B)  ! subtract average^2
   Sig(2)=sqrt(Sig(2)/Nr_IntFer)
   Sig(3)=sqrt(Sig(3)/Nr_IntFer)
   !   Delta_t= -dCenT(m) +SUM(Vec_l(1:3,j_IntFer)*dCenloc(1:3,m))/(c_mps*sample)   ! in [samples]
   !dd=dt*c_mps*sample  ! [samples] * [m/s] * [s/sample] =[m]
   Gr(0,0)= 1./(1-B*B)    ;  Gr(0,1:3)=Gr(1,1:3)*B*c_mps*sample/(1-B*B)
   Gr(1,0)= B/Sig(1)      ;  Gr(1,1:3)=Gr(1,1:3)*c_mps*sample/Sig(1)
   Gr(2,1:3)=Gr(2,1:3)*c_mps*sample/Sig(2)
   Gr(3,1:3)=Gr(3,1:3)*c_mps*sample/Sig(3)
   !
   ! Adepted grid???
   Gr(0,0)= 1.    ;  Gr(0,1:3)=0.
   !
   Gr(0,0:3)=Gr(0,0:3)*IRFW_s *GridScale ! *3.
   Gr(1:3,0:3)=Gr(1:3,0:3)*IRFW_s *GridScale  ! *2.
   !
   write(2,*) ' MDD_IniManySrcs, B, sig',B, sig(1:3)
   Length(0)=Gr(0,0)**2 - SUM(Gr(0,1:3)*Gr(0,1:3))/((c_mps*sample)**2)
   Do i=1,3
      Length(i)=SUM(Gr(i,1:3)*Gr(i,1:3))  ! take only space part
   EndDo
   Length(0)=sqrt(Length(0))
   Length(1:3)=sqrt(Length(1:3))
   SpaceGrid(1:3)=Length(1:3)   ! need for curtain plots
   !
   !SpaceGridScale(:)=sqrt(D/SpaceGridScale(:))  !  sqrt creates large dependencies (=over fitting)
   Write(2,"(A,F5.1,A,3F6.1,3F6.1)")  'Impulse response width=',IRFW_s, &
      '[samples], Grid spacings (t,N,E,h) in [samples,m,m,m]; Length'
   Write(2,"(2x,F6.1,3F6.1,2x,F6.1,f9.3)")   Gr(0,0:3), Length(0), ( Gr(0,0)-SUM(Ave_l(1:3)*Gr(0,1:3))/(c_mps*sample))
   Write(2,"(2x,F6.1,3F6.1,2x,F6.1,f9.3)")   Gr(1,0:3), Length(1), ( Gr(1,0)-SUM(Ave_l(1:3)*Gr(1,1:3))/(c_mps*sample))
   Write(2,"(2x,F6.1,3F6.1,2x,F6.1)")   Gr(2,0:3), Length(2)
   Write(2,"(2x,F6.1,3F6.1,2x,F6.1)")   Gr(3,0:3), Length(3)
   !Write(2,*) 'Bias_inv_base, DDXcorr_depl, Trace_DeadEnds:',Bias_inv_base, DDXcorr_depl, Trace_DeadEnds, RealDip_logic
   !Write(2,"(I2,F6.1,3F6.1)")  ( i,Gr(i,0:3), i=0:3)
   !
   !
   Ngrd(0)=(T_Range-2*Trace_DeadEnds)/Length(0)  !
   If(Ngrd(0).lt.1) Ngrd(0)=1
   If(Ngrd(0).gt.9) Ngrd(0)=9
   !
   ! Determine steps in space directions from average dot product longitudinal vector
   NSrc_max=Ngrd(0)
   Do i=1,3
      Ngrd(i) = NINT(2*dD_max/Length(i))  ;   If(Ngrd(i).lt.1) Ngrd(i)=1
      Ngrd(i)=2*(Ngrd(i)/2)+1  !  make it an odd number to have at least one grid point at the center
      If(Ngrd(i).gt.9) Ngrd(i)=9
      NSrc_max=NSrc_max*Ngrd(i)
   EndDo
   !
   Allocate(  dCenT(1:NSrc_max) ) ! time difference (in [samples]) of each DD from central time  C_Tms
   Allocate( dCenloc(1:3,1:NSrc_max) ) ! location of each DD w.r.t. C_Loc,  in [m]
   Allocate( DelDip(1:3,1:NSrc_max) ) ! Dipole moment of each DD
   !
   m=0
   i_gr(0)=0
   Do while (i_gr(0).lt.Ngrd(0))
      i_gr(0)=i_gr(0)+1
      i_gr(1)=0
      Do while (i_gr(1).lt.Ngrd(1))
         i_gr(1)=i_gr(1)+1
         i_gr(2)=0
         Do while (i_gr(2).lt.Ngrd(2))
            i_gr(2)=i_gr(2)+1
            i_gr(3)=0
            Do while (i_gr(3).lt.Ngrd(3))
               i_gr(3)=i_gr(3)+1
               m=m+1
               dCenT(m)=SUM( Gr(0:3,0)*(i_gr(0:3)-0.5*(Ngrd(0:3)+1.)) ) ! makes evenly spread around 0
               dCenloc(1,m)=SUM( Gr(0:3,1)*(i_gr(0:3)-0.5*(Ngrd(0:3)+1.)) )
               dCenloc(2,m)=SUM( Gr(0:3,2)*(i_gr(0:3)-0.5*(Ngrd(0:3)+1.)) )
               dCenloc(3,m)=SUM( Gr(0:3,3)*(i_gr(0:3)-0.5*(Ngrd(0:3)+1.)) )
               Do j_IntFer=1,Nr_IntFer
                  del_t= -dCenT(m) +SUM(Vec_l(1:3,j_IntFer)*dCenloc(1:3,m))/(c_mps*sample)   ! in [samples]
                  If(abs(Del_t) .gt. (T_Range-2.*Trace_DeadEnds)/2.) Then
                     !write(2,*) 'cancelled:',m,j_IntFer,dCenT(m), dCenloc(:,m), del_t
                     m=m-1  ! Throw away this dd-location as the pulse will be too close to the edge of a spectrum
                     exit
                  EndIf
               EndDo ! j_IntFer=1,Nr_IntFer
               !If(j_IntFer .eq. 1+Nr_IntFer) Then  ! not ued previous exit statement
               !   write(2,"(A,i5,4F5.1,2x,4(F7.0))") 'Kept gridpoint:',m,i_gr(0:3)-0.5*(Ngrd(0:3)+1.),dCenT(m),dCenloc(1:3,m)
               !EndIf
            EndDo ! i_3=0,Ngrd
         EndDo ! i_2=0,Ngrd
      EndDo ! i_1=0,Ngrd
   EndDo ! i_t=0,Ngrd
   Nsrc=m
   Write(2,"(A,I5,A,I8,A,4I3)") 'MDD_IniManySrcs: N_after=', Nsrc, ', N_before space cuts=',NSrc_max, &
      ', N_grid(t,N,E,h)=', Ngrd(0:3)
   Return
End Subroutine MDD_IniManySrcs
!-----------------------------------------------
Subroutine MDD_EliminateSrcs()
!  Keep the NSrc_max brightest sources
   use constants, only : dp, sample, c_mps
   Use MDD_Pars, only : NSrc_max, IRFW_s, T_Range, Bias_inv_base, NSrc_min, IntensityFrac_keep  ! Input
   Use MDD_Pars, only : NSrc, dCenT, dCenloc, Bias_inv, MDDLabel  ! Output
   Use MDD_Pars, only : DDChiSQ, Weight
   Use MDD_Pars, only : DelDip, RealDip, RealDip_logic
   Use unque, only : Double_RI_sort  ! (N,Real(dp),integer)
   Implicit none
   !Real(dp) :: dt, dd, D
   Real(dp), allocatable :: Temp(:), dtsort(:), dlsort(:,:)
   Complex(dp), allocatable :: ddsort(:,:)
   Integer, allocatable :: Permut(:)
   Integer :: m,n, NSrc_new, ReductionRound
   Real(dp) :: Q_max
   Real(dp) :: FitQual
   Integer :: Error
   !
   NSrc_max=NSrc_min
   RealDip=.false.
   If((RealDip_logic.eq.1) .or. (RealDip_logic.eq.3)) RealDip=.true.
   !
   Do ReductionRound=1,9
      If(NSrc .le. NSrc_min) goto 9  ! no readj needs to be done
      !write(MDDLabel(5:6),"('-',i1)") ReductionRound
      MDDLabel(5:5)=achar(iachar("a")+ReductionRound-1)
      Bias_inv=NSrc*NSrc*Bias_inv_base
      write(2,*) 'ReductionRound:',MDDLabel(4:5), Bias_inv
      Call MDD_QuickFit(FitQual,error)
      Call MDD_write()
      !
      Allocate( Temp(1:NSrc), Permut(1:NSrc) )
      Do m=1,NSrc
         Permut(m)=m
         Temp(m)=SUM( ABS( DelDip(:,m) )**2 )
      EndDo
      Call Double_RI_sort(NSrc,Temp,Permut)
      !
      Do m=1,NSrc  ! find Frac from max
         If(Temp(m) .gt. IntensityFrac_keep*Temp(NSrc)) exit
      EndDo ! m=1,NSrc
      !NSrc_new=NSrc-m
      NSrc_new=NSrc/2
      If((NSrc-m).lt.NSrc_new) NSrc_new=NSrc-m
      If(NSrc_new.lt.NSrc_min) NSrc_new=NSrc_min
      If(ReductionRound.eq.9) NSrc_new=NSrc_min
      Allocate( dtsort(1:NSrc_new), dlsort(1:3,1:NSrc_new), ddsort(1:3,1:NSrc_new) )
      !
      Write(2,"(A,A2, A, 1pG10.2,A,1pG10.2,A, I4,1pG10.2)") 'ReductionRound:',MDDLabel(4:5), &
         ', range from I=', Temp(NSrc),' till', Temp(1), ', kept ', NSrc_new, Temp(NSrc-NSrc_new+1)
      Write(*,"(A,1pG10.2,A,1pG10.2)") 'Rastered kept I=', Temp(NSrc),' till', Temp(NSrc-NSrc_new+1)
      Do m=1,NSrc_new
         n=Permut(NSrc-m+1)
         dtsort(m)=dCenT(n)
         dlsort(:,m)=dCenloc(:,n)
         ddsort(:,m)=DelDip(:,n)
      EndDo ! m=1,NSrc_max
      !
      ! Get organized for the next round
      DeAllocate( dCenT, dCenloc, DelDip )
      Allocate(  dCenT(1:NSrc_new) ) ! time difference (in [samples]) of each DD from central time  C_Tms
      Allocate( dCenloc(1:3,1:NSrc_new) ) ! location of each DD w.r.t. C_Loc,  in [m]
      Allocate( DelDip(1:3,1:NSrc_new) ) ! Dipole moment of each DD
      dCenT(:)=     dtsort(:)
      dCenloc(:,:)= dlsort(:,:)
      DelDip(:,:)=ddsort(:,:)
      DeAllocate( Temp, Permut, dtsort, dlsort, ddsort )
      Nsrc=NSrc_new
      DDChiSQ(:,:)=0.
   EndDo ! ReductionRound=1,10
   !Call MDD_write()
9  Continue
   Bias_inv=Bias_inv_base
   MDDLabel(5:6)='  '
   RealDip=.false.
   If((RealDip_logic.eq.2) .or. (RealDip_logic.eq.3)) RealDip=.true.
   Return
End Subroutine MDD_EliminateSrcs
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
   Integer :: GnuDim, i_nu, inu1, inu2, i_freq,i,j
   Real(dp) :: dnu, Q, nu, dfreq, RG_ns(1:G_dim)
   !
   Call RFTransform_su(G_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   RG_ns(:)=0.d0
   RG_ns(G_dim/2)=1.d3  ! such have a decent normalisation
   Call RFTransform_CF(RG_ns,Gnu(0))
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
   Call RFTransform_CF2CT(Gnu(0),G_ns(1) )
   Q=Abs(G_ns(G_dim/2))
   G_ns(:)=G_ns(:)
   Do i=G_dim/2,G_dim
      nu=Abs(G_ns(i))/Q
      If(nu .lt. 0.5d0) exit
   EndDo
   IRFW_s=(2*i-G_dim)/5.
   write(2,"(A,20F5.2)") 'Abs(G_ns):',Abs(G_ns(G_dim/2:G_dim/2+19))*10./Q
   write(2,"(A,20F5.2)") ' --40:',Abs(G_ns(G_dim/2+20:G_dim/2+39))*10./Q
   write(2,"(A,20F5.2)") ' --60:',Abs(G_ns(G_dim/2+40:G_dim/2+59))*10./Q
   Do i=G_dim/2,G_dim
      nu=RealPart(G_ns(i))
      If(nu .lt. 0.) exit
   EndDo
   write(2,"(A,50('(',2F5.2,')',1x))") 'G_ns(max:i[ns]):',G_ns(G_dim/2:i)
   write(2,*) 'Impulse Response FWHM(amplitude)=',IRFW_s,'[samples], peak=', Q, ', Real Zeros@',(2*i-G_dim)/5.
   !G_ns(:)=Real(CG_ns(:))
   !write(2,"(100F5.2)") ABS(G_ns(G_dim/2:G_dim/2+40))
   !
   call DAssignFFT()                   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   Return
End Subroutine ImpResp_ns
!-----------------------------------------------
Subroutine ReadMDD_C(lname)
! Read some key steering parameters for MDD mode
   use constants, only : dp, pi
   Use MDD_Pars, only : Bias_inv_base, DDXcorr_depl, Trace_DeadEnds, RealDip_logic, NSrc_min, IntensityFrac_keep  ! output
   Use MDD_Pars, only : DepleteCoreWeigth, DepleteDistFrac, GridScale  ! output
   Use Interferom_Pars, only : IntFer_ant, i_chunk, Nr_IntFerCh
   use Chunk_AntInfo, only : Ant_Stations
   use StationMnemonics, only : Station_Mnem2ID, Station_ID2Mnem
   Implicit none
   Character(LEN=*), Intent(INout) :: lname
   Integer :: TxtLength, nxx, i
   NAMELIST /MDD_cntrl/ Bias_inv_base, DDXcorr_depl, Trace_DeadEnds, RealDip_logic, NSrc_min, IntensityFrac_keep, &
      DepleteCoreWeigth, DepleteDistFrac, GridScale
   !
   TxtLength=LEN_TRIM(lname)
   Do i=TxtLength,1,-1
      lname(i+1:i+1)=lname(i:i)
   EndDo
   lname(1:1)='&'
   !
   Read(lname, NML = MDD_cntrl, iostat=nxx)
   If(nxx.ne.0) Then
      write(2,*) 'Error in reading &MDD_cntrl: ',lname
      write(2,*) 'Expected one of the following parameters'
      write(2,NML = MDD_cntrl)
      Flush(unit=2)
   EndIf
   If(Bias_inv_base.lt.0.d0) Bias_inv_base=tiny(Bias_inv_base) ! smallest number just larger than zero
   If(DDXcorr_depl.lt.0.d0) DDXcorr_depl=EPSILON(DDXcorr_depl) ! smallest number to add to 1 and make a difference
   If(DDXcorr_depl.gt.1.d0) DDXcorr_depl=1.d0-EPSILON(DDXcorr_depl) ! smallest number to add to 1 and make a difference
   If((RealDip_logic.lt.0) .or. (RealDip_logic.gt.3)) RealDip_logic=0
   If(Trace_DeadEnds.lt.1) Trace_DeadEnds=10
   If(NSrc_min.lt.1) NSrc_min=10
   If( (IntensityFrac_keep.lt.0.) .or. (IntensityFrac_keep.gt.1.) ) IntensityFrac_keep=0.3
   If( (DepleteCoreWeigth.lt.0) .or. (DepleteCoreWeigth.gt. 1.) ) DepleteCoreWeigth=1.d0
   If( (DepleteDistFrac.lt.0) .or. (DepleteDistFrac.gt.1.d0) ) DepleteDistFrac=1.d0
   !
   Return
End Subroutine ReadMDD_C
!-----------------------------------------------
Subroutine WriteMDD_C()
! Read some key steering parameters for MDD mode
   !use constants, only : dp, pi
   Use MDD_Pars, only : Bias_inv_base, DDXcorr_depl, Trace_DeadEnds, RealDip_logic, NSrc_min, IntensityFrac_keep  ! input
   Use MDD_Pars, only : DepleteCoreWeigth, DepleteDistFrac, GridScale  ! input
   Implicit none
   !Character(LEN=*), Intent(INout) :: lname
   !Integer :: TxtLength, nxx, i
   !
   write(2,"(A,1pG12.3)") ' - Impulse-response grid spacing multiplied by GridScale=',GridScale
   write(2,"(A,1pG12.3)") ' - Bias added to the diagonal of A, Bias_inv_base=',Bias_inv_base
   write(2,"(A,1pG12.3)") ' - (1-DDXcorr_depl) is factor to diminish off-diagonal correlations',DDXcorr_depl
   write(2,"(A,I4)") ' - Trailing parts of traces devoid of sources, Trace_DeadEnds=',Trace_DeadEnds
   write(2,"(A,I4,A,A)") ' - RealDip_logic=',RealDip_logic,' 0: never set RealDip=.t.; 1: RealDip=.t. in gridsearch;',&
      ' 2: RealDip=.t. in final;  3: all RealDip=.t.;'
   write(2,"(A,1pG12.3)") ' - Kept raction of max source intensity IntensityFrac_keep=', IntensityFrac_keep
   write(2,"(A,I2)") ' - Number of sources kept for fine location search, NSrc_min=', NSrc_min
   write(2,"(A,1pG12.3)") ' - Weights of the CS stations is decreased by factor=',DepleteCoreWeigth
   write(2,"(A,1pG12.3)") ' - Weights kept constant for antennas closer than the core with distance fraction=',DepleteDistFrac
   !
   Return
End Subroutine WriteMDD_C
!-----------------------------------------------
Subroutine ReadTeststations(lname)
   use constants, only : dp, pi
   Use MDD_Pars, only : Teststations_MN, TeststAntenna_ID, NTestSt_max, NTestAnt  ! output
   Use Interferom_Pars, only : IntFer_ant, i_chunk, Nr_IntFerCh
   use Chunk_AntInfo, only : Ant_Stations
   use StationMnemonics, only : Station_Mnem2ID, Station_ID2Mnem
   Implicit none
   Character(LEN=*), Intent(IN) :: lname
   Integer :: i,k, j_IntFer, i_ant, nxx
   !
   NTestAnt=1
   TeststAntenna_ID(NTestAnt)=1
   Call Station_ID2Mnem(Ant_Stations(TeststAntenna_ID(NTestAnt),1),Teststations_MN(NTestAnt))
   Teststations_MN(2:)='     '
   read(lname,*,iostat=nxx) Teststations_MN(2:) ! option=abut,only
   Do i=2,NTestSt_max
      !write(2,*) 'FP_Mnem(i):',i,FP_Mnem(i)
      If(Teststations_MN(i).eq.'     ') exit
      If(TRIM(Teststations_MN(i)).eq.'!') exit
      Call Station_Mnem2ID(Teststations_MN(i),k)
      !write(2,*) 'Teststations:',i, Teststations_MN(i),k
      If(k.eq.0) cycle ! exit
      !write(2,*) i,FP_Mnem(i),FP_s(i)
      Do j_IntFer=1,Nr_IntFerCh(1)   ! Find antennas in selected stations
         i_ant=IntFer_ant(j_IntFer,1)  ! Refers to even antenna
         If(k.eq.Ant_Stations(i_ant,1)) Then
            NTestAnt=NTestAnt+1
            TeststAntenna_ID(NTestAnt)=j_IntFer
            Teststations_MN(NTestAnt)=Teststations_MN(i)  ! keep only the names of existing stations
            exit
         EndIf
      EndDo !  j_IntFer=1,Nr_IntFerCh(1)
   Enddo !  i=1,NTestSt_max
   write(2,"(11(A,1x))") 'Test Stations:', Teststations_MN(1:NTestAnt)
   !
End Subroutine ReadTeststations
!-----------------------------------------------
Subroutine MDD_write()
!
!-----------------------------------------------
   use constants, only : dp, pi, c_mps, sample
   Use Interferom_Pars, only : IntFer_ant, i_chunk, Nr_IntFerCh
   Use MDD_Pars, only : IRFW_s
   Use MDD_Pars, only : NSrc, dCenT, dCenloc, T_Range, Vec_l
   Use MDD_Pars, only : DDChiSQ, Weight
   Use MDD_Pars, only : DelDip
   Use MDD_Pars, only : C_Tms, C_Loc
   Use MDD_Pars, only : Teststations_MN, TeststAntenna_ID, NTestAnt, NTestSt_max  ! input
   Use unque, only : Double_RI_sort  ! (N,Real(dp),integer)
   Implicit none
   Integer :: i, n, m, Nr_IntFer, j_IntFer, i_ant
   Real(dp) :: Q, dt_TestStation(1:NTestSt_max)
   Real(dp), allocatable :: Temp(:)
   Integer, allocatable :: Permut(:)
   Character(len=100) :: FMT1, FMT2, FMT3
   !
   Allocate( Temp(1:NSrc), Permut(1:NSrc) )
   j_IntFer=1
   Do i=1,NSrc
      Permut(i)=i
      !Temp(i)=dCenT(i)
      Temp(i)=(dCenT(i) - SUM(Vec_l(1:3,j_IntFer)*dCenloc(1:3,i))/(c_mps*sample))  ! arrival time in the core
   EndDo
   Call Double_RI_sort(NSrc,Temp,Permut)
   !
   !Write(2,"('  #  d_smpl',T15,'d_(N,E,h) [m]', T43, 'dipole moments (N,E,h) [?]')")
   write(FMT1,"(A,I0,A)") "('  #  d_smpl',T15,'d_(N,E,h) [m]', T33,'t in:',", &
      NTestAnt,"A7,8x, 'dipole moments (N,E,h) [?]')"
   !write(2,*) 'FMT1:',FMT1
   !Write(2,"(I3,F8.2, 3F8.2,2x,F9.1,A,'[',3('(',F5.2,','F5.2,')',1x),']' )") &
   write(FMT2,"(A,I0,A)") '(I3,F8.2, 3F8.2,2x,', &
      NTestAnt,"F7.1,2x,F9.1,A,'[',3('(',F5.2,','F5.2,')',1x),']')"
   !write(2,*) 'FMT2:',FMT2
   !Flush(unit=2)
   !
   write(2,"(A,F10.5,A,2(F7.3,','),F6.3,A,I4,A,I3,A )") 'C-time=',C_Tms,'[ms], C-loc (N,E,h)=(',C_Loc(1:)/1000., &
      ') [km], Fit window=',T_Range,' samples with',NSrc, ' Delta Dipole sources'
   Write(2,FMT1) Teststations_MN(1:NTestAnt)
   Flush(unit=2)
   Do n=1,NSrc
      m=Permut(n)
      !Write(2,"(I3,F8.2, 3F8.2, 2x,3('(',F9.2,SP,F9.2,'i)',1x) )") m,dCenT(m), dCenloc(1:3,m), DelDip(1:3,m)
      Q=sqrt(SUM( ABS( DelDip(:,m) )**2 ) )
      !Write(2,"(I3,F8.2, 3F8.2 &
      !!   ,2x,F9.1,1x,3(F9.2,'(',F5.0,')',1x)  &
      !   ,2x,F9.1,A,'[',3('(',F5.2,','F5.2,')',1x),']'   &
      !    )") &
      Do i=1,NTestAnt
         j_IntFer=TeststAntenna_ID(i)
         dt_TestStation(i) = dCenT(m) - SUM(Vec_l(1:3,j_IntFer)*dCenloc(1:3,m))/(c_mps*sample)  ! in [samples]
      EndDo ! i=1,NTestAnt
      Write(2,FMT2) n,dCenT(m), dCenloc(1:3,m), dt_TestStation(1:NTestAnt) &
      !   , Q, (ABS( DelDip(i,m) ) , Atan2(ImagPart(DelDip(i,m)),RealPart(DelDip(i,m)))*180./pi,i=1,3)  &
      !   , Q,' x ', (RealPart(DelDip(i,m))/Q,ImagPart(DelDip(i,m))/Q,i=1,3)
         , Q*Q, ' ', (RealPart(DelDip(i,m))/Q,ImagPart(DelDip(i,m))/Q,i=1,3)
      !If(NSrc.gt.100) then
      If(N.gt.225) then
         write(2,*) 'More than 100 DD sources, N=',NSrc
         exit
      EndIf
   EndDo
   Nr_IntFer=Nr_IntFerCh(1)
   Q=SUM(Weight(1:Nr_IntFer))  ! *T_Range
   write(2,"(A,2F10.4,A)") 'Mean chi^2 _phi & _th:',SUM(DDChiSQ(1:Nr_IntFer,1))/Q, &
      SUM(DDChiSQ(1:Nr_IntFer,2))/Q, ' [amplitude^2/sample]'
   flush(unit=2)
   DeAllocate( Temp, Permut )
   Return
End Subroutine MDD_write
!-----------------------------------------------
Subroutine MDD_GLEplot(peaknr)
!
!-----------------------------------------------
   use constants, only : dp, pi, c_mps, sample
   !Use Interferom_Pars, only : IntFer_ant, i_chunk, Nr_IntFerCh
   use DataConstants, only : DataFolder, OutFileLabel
   !Use MDD_Pars, only : IRFW_s
   Use MDD_Pars, only : NSrc, dCenT, dCenloc, T_Range, SpaceGrid
   !Use MDD_Pars, only : DDChiSQ, Weight
   Use MDD_Pars, only : DelDip
   Use MDD_Pars, only : C_Tms, C_Loc
   use GLEplots, only : GLEplotControl
   Implicit none
   Integer, Intent(IN) :: peaknr
   Integer :: i, n, m !, Nr_IntFer, Permut(100)
   Real(dp) :: dd, dt, Q, T_unts !, Temp(100)
   character(len=1) :: txt
   !
   write(txt,"(I1)") peaknr
   dt=T_Range*sample/2  ! [samples] * [s/sample] =[s]
   dd=dt*c_mps  ! [s] * [m/s]  =[m]
   T_unts = 1.d6  !  *T_unts  converts [s] to [{\mu}s]
   dt=dt*T_unts ! [ms]
   Write(2,*) 'Box-size for plot:',dt,' [\mu s],',dd,' [m]x',SpaceGrid
   OPEN(UNIT=29,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//txt//'MDDPow.dat')
   Q=SpaceGrid(1)
   If(SpaceGrid(2).gt.Q) Q=SpaceGrid(2)
   If(Q.gt. 500.) Q=500.
   write(29,"(6F8.3,2F9.3,A,F12.6,i3,' 0')") &
      (C_Loc(2)-Q)/1000., (C_Loc(2)+Q)/1000., (C_Loc(1)-Q)/1000., (C_Loc(1)+Q)/1000., &
      (C_Loc(3)-SpaceGrid(3))/1000., (C_Loc(3)+SpaceGrid(3))/1000., -dt, +dt, ' NoBox ', C_Tms, 0 ! gleRead:  xMin xMax yMin yMax zMin zMax tMin tMax ZmBx$ t_offst
   Write(29,*) '0 ',NSrc,' 0 0 ',TRIM(OutFileLabel)//txt, 10.,' 0 0 0 0 0 0 0 ',0.1 ,  '1.0 !' ! gleRead:  NTracks EventNr Q t_start label$ AmplitudePlot a1 b1 c1 a2 b2 c2 d2
   !
   write(29,"('! MDD#  time[mu s]',T22,'x=East',T35,'y=North',T51,'z=h',T63,'Ampl',T77,'(Re,Im) pol (N,E,h)')")
   Do m=1,NSrc
     Q=sqrt(SUM( ABS( DelDip(:,m) )**2 ) )
      write(29,"(i6,' ',f9.3, 3(f12.5,' '),F11.2,' ',6f9.2)") &
         m, dCenT(m)*sample*T_unts, (C_Loc(2)+dCenloc(2,m))/1000., (C_Loc(1)+dCenloc(1,m))/1000., (C_Loc(3)+dCenloc(3,m))/1000. &
         , Q, DelDip(1:3,m)
   EndDo
   Close(UNIT=29)
   !
   !
   Call GLEplotControl(PlotType='SourcesPlotBckgr', PlotName='MDDPk'//TRIM(OutFileLabel)//txt, &  !  SourcesPlotBckgr
      PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel)//txt//'MDDPow' )
   Return
End Subroutine MDD_GLEplot
!-------------------------------------
Subroutine MDDClose()
   Use MDD_Pars, only : Weight, SQWght, CTime_p, CTime_t, Vec_p, Vec_t, Vec_l, AntSourceD
   Use MDD_Pars, only : DDChiSQ
   Use MDD_Pars, only : dCenT, dCenloc, DelDip
   use FFT, only : DAssignFFT
   !
   call DAssignFFT()                   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   Deallocate (CTime_p, CTime_t, Vec_p, Vec_t, Vec_l, Weight, SQWght, AntSourceD, DDChiSQ )
   DeAllocate( dCenT, dCenloc, DelDip )
   !
   Return
End Subroutine MDDClose
