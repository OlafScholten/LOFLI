    !
!-----------------------------------------------
Subroutine EI_PolarizW(Nr_IntFer, IntfNuDim, i_slice)
   ! Needs:
   !   ?
   !--------------------------------------------
   use constants, only : dp, pi, ci
   use DataConstants, only : Time_Dim, Cnu_dim
   use DataConstants, only : DataFolder, OutFileLabel
   use Chunk_AntInfo, only : CTime_spectr, Ant_Stations, Ant_IDs, Ant_nr, Ant_pos, Ant_RawSourceDist
   use Chunk_AntInfo, only : Start_time, TimeBase
   use Chunk_AntInfo, only : Powr_eo,NAnt_eo
   use FFT, only : RFTransform_CF, RFTransform_CF2CT
   Use Interferom_Pars, only : IntfBase, IntfLead, N_smth, smooth, PowerScale
   Use Interferom_Pars, only : IntFer_ant
   Use Interferom_Pars, only : i_chunk, t_shft, PixLoc, CenLoc
   Use Interferom_Pars, only : StI, StI12, StQ, StU, StV, StI3, StU1, StV1, StU2, StV2
   Use Interferom_Pars, only : dStI, dStI12, dStQ, dStU, dStV, dStI3, dStU1, dStV1, dStU2, dStV2, Chi2pDF
   use AntFunCconst, only : Freq_min, Freq_max,Ji_p0,Ji_t0,Ji_p1,Ji_t1, Gain  !J_0p,J_0t,J_1p,J_1t,
   use GLEplots, only : GLEplotControl
   Implicit none
   Integer, intent(in) :: Nr_IntFer, IntfNuDim, i_slice
   complex(dp), parameter :: ipi=ci*pi
   integer :: i_ant, j_IntFer, i, j, i_freq, i_nu, IntfDim, inu1, inu2, i_s, m, n
   Real(dp) :: Vec_p(1:3), Vec_t(1:3), Ras(1:3),  p_PB(1:3), t_PB(1:3)
   Real(dp) :: NEh2PB(1:3,1:3) ! NEh to Polarization Basis; NEh2VHR(NEh,VHR); formerly known as VHRBasis(1:3,1:3)
   Real(dp) :: HorDist, Thet_r ,Phi_r, dfreq, D, W, dnu, nu
   Real(dp) :: thet_d, Phi_d  ! AntSourceD(1:Nr_IntFer),
   Real(dp) :: Rtime(1:2*IntfNuDim)
   Complex(dp) :: Sp, St
   Complex(dp) :: FTime_PB(1:2*IntfNuDim,1:3), Fnu_PB(0:IntfNuDim,1:3)
   Complex(dp) :: Phase, dPhase, nu_p(0:IntfNuDim), nu_t(0:IntfNuDim)
   Real(dp) :: NormEven, NormOdd, RDist, dt_AntPix
   Complex(dp), Allocatable, save :: CTime_p(:,:), CTime_t(:,:)
   Complex(dp), Allocatable, save :: Cnu0(:,:), Cnu1(:,:)
   Real(dp), Allocatable, save :: Noise_p(:), Noise_t(:), Power_p(:), Power_t(:)
   Real(dp) :: Pow_p, Pow_t, W_p, W_t
   Real(dp), save :: I_scale
   Real(dp) :: Cur2E(1:3,1:3), Ai(1:3,1:3)
   Complex :: Stk(3,3), AiF(1:3)
   logical :: First=.true.
   Real(dp) :: Esq_ak, FdotI
   !
   !logical :: TestCh2=.true.
   logical :: TestCh2=.false.
   Real(dp) :: sw_p, sw_t
   Real(dp), Allocatable :: D_a(:),Ph_a(:)  ! weighted p and t directions in polarization basis (_PB)
   Real(dp), Allocatable :: wap_PB(:,:),wat_PB(:,:)  ! weighted p and t directions in polarization basis (_PB)
   Complex(dp), Allocatable :: wEnu_p(:), wEnu_t(:)
   Complex(dp), Allocatable :: wETime_p(:), wETime_t(:)
   Complex(dp), Allocatable :: wEtime_ap(:,:), wEtime_at(:,:) ! weighted E fields, time dependent
   !Complex(dp), Allocatable :: dChi_ap(:), dChi_at(:)
   Real(dp), Allocatable :: dChi_ap(:), dChi_at(:)
   Real(dp), Allocatable :: W_ap(:), W_at(:)
   !Complex(dp) :: SumDiff, SumSq, Del_p, Del_t
   Real(dp) :: SumDiff, SumSq, Del_p, Del_t
   character(len=2) :: txt
   !
   IntfDim=2*IntfNuDim
   dnu=100./IntfNuDim   ! [MHz] Jones matrix is stored on 1MHz grid
   inu1=Int(Freq_min/dnu)
   inu2=Int(Freq_max/dnu-0.5)+1  ! set integration regime to frequencies within filter band
   If(First) Then  ! preparation
      First=.false.
      Allocate( CTime_p(1:2*IntfNuDim,Nr_IntFer), CTime_t(1:2*IntfNuDim,Nr_IntFer) )
      Allocate( Cnu0(0:IntfNuDim,Nr_IntFer), Cnu1(0:IntfNuDim,Nr_IntFer) )
      Allocate( Noise_p(1:Nr_IntFer), Noise_t(1:Nr_IntFer) )
      Allocate( Power_p(0:Nr_IntFer), Power_t(0:Nr_IntFer) )
      Power_p(:)=0.
      Power_t(:)=0.
      !
      NormEven=sqrt(2.*Powr_eo(0)/(PowerScale*(Powr_eo(0)+Powr_eo(1))) )/100.
      NormOdd=sqrt(2.*Powr_eo(1)/(PowerScale*(Powr_eo(0)+Powr_eo(1))) )/100.
      Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas to detremine t- and p- polarized trace for central pixel; needed for noise estimate
         i_ant=IntFer_ant(j_IntFer)
         !write(2,*) i_ant, j_IntFer, Ant_pos(:,i_ant,i_chunk)
         Ras(1)=(CenLoc(1)-Ant_pos(1,i_ant,i_chunk))/1000.
         Ras(2)=(CenLoc(2)-Ant_pos(2,i_ant,i_chunk))/1000.
         Ras(3)=(CenLoc(3)-Ant_pos(3,i_ant,i_chunk))/1000.
         HorDist=sqrt(  Ras(1)*Ras(1) + Ras(2)*Ras(2)  )
         Thet_r=atan(HorDist,Ras(3))  ! Zenith angle
         Phi_r=atan2( Ras(2) , Ras(1) )  ! \phi=0 = north
         Thet_d =Thet_r*180/pi
         Phi_d =Phi_r*180/pi
         If(j_IntFer.eq.1) Then
            I_Scale= (HorDist*HorDist + Ras(3)*Ras(3))/Nr_IntFer ! scales value of chi-square, not intensity, and value of F
            !write(2,*) 'AmplitudeW NormEven,odd',NormEven,NormOdd,', ratio=',NormEven/NormOdd, Noise
         EndIf
         !
         !Calculation of noise-power level
         Call AntFun_Inv(thet_d ,Phi_d ) ! sets ,Ji_p0,Ji_t0,Ji_p1,Ji_t1; Inverse Jones; gives _p & _t polarized fields when contracted with (0.1) voltages
         w_p=SUM( ((NormEven*ABS(Ji_p0(Freq_min:Freq_max)))**2+(NormOdd*ABS(Ji_p1(Freq_min:Freq_max)))**2)* &
            Gain(Freq_min:Freq_max)**4)
         w_t=SUM( ((NormEven*ABS(Ji_t0(Freq_min:Freq_max)))**2+(NormOdd*ABS(Ji_t1(Freq_min:Freq_max)))**2)* &
            Gain(Freq_min:Freq_max)**4)
         !   noise=100./W  ! in 100./W Factor 100 to copy the factor that was introduces in Antenna_Read
         Noise_p(j_IntFer)=W_p*1.4d5  ! *10.
         Noise_t(j_IntFer)=W_t*1.4d5  ! *10. ! such a scaling factor changes the chi-square value but not the error-bars         !
         RTime(:)=REAL(CTime_spectr(IntfBase:IntfBase+IntfDim,i_ant,i_chunk))*NormEven
         Call RFTransform_CF(RTime,Cnu0(0,j_IntFer))
         RTime(:)=REAL(CTime_spectr(IntfBase:IntfBase+IntfDim,i_ant+1,i_chunk))*NormOdd
         Call RFTransform_CF(RTime,Cnu1(0,j_IntFer))
         !
         nu_p(:)=0.
         nu_t(:)=0.
         Do i_nu=inu1,inu2   ! Increment frequency spectrum with this antenna
            nu=i_nu*dnu
            i_freq=Int(nu)
            dfreq=nu-i_freq ! phase-shifts are zero for centran pixel, no timeshift!!
            nu_p(i_nu) =(  ((1.-dfreq)*Ji_p0(i_freq) + dfreq*Ji_p0(i_freq+1)) *Cnu0(i_nu,j_IntFer) + &
                           ((1.-dfreq)*Ji_p1(i_freq) + dfreq*Ji_p1(i_freq+1)) *Cnu1(i_nu,j_IntFer) ) *Gain(i_freq)
            nu_t(i_nu) =(  ((1.-dfreq)*Ji_t0(i_freq) + dfreq*Ji_t0(i_freq+1)) *Cnu0(i_nu,j_IntFer) + &
                           ((1.-dfreq)*Ji_t1(i_freq) + dfreq*Ji_t1(i_freq+1)) *Cnu1(i_nu,j_IntFer) ) *Gain(i_freq)
            ! Gain(i_freq)=sqrt( SUM(J^2) )/(Freq_max-Freq_min) ; to neutralize effect Ji on frequency, calculated in antenna_function
         Enddo
         ! convert to time
         Call RFTransform_CF2CT(nu_p(0),CTime_p(1,j_IntFer) )
         Call RFTransform_CF2CT(nu_t(0),CTime_t(1,j_IntFer) )
      Enddo    !  j_IntFer=1,Nr_IntFer
      !
   EndIf  ! End preparation
   !
   ! Just for checking causality
   write(2,*) 'central pixel from slice center for first antenna pair:'
   j_IntFer=1
   i_ant=IntFer_ant(j_IntFer)
   Call RelDist(PixLoc(1),Ant_pos(1,i_ant,i_chunk),RDist)
   dt_AntPix =Rdist - Ant_RawSourceDist(i_ant,i_chunk)
   i_s=1+i_slice*N_smth+NINT(dt_AntPix )  ! approximately correct, upto rounding errors for dt
   !write(2,*) 'p_up',  (Abs(CTime_p(IntfLead+i_s+j,j_IntFer))**2, j=0,N_smth)
   !write(2,*) 'p_dwn', (Abs(CTime_p(IntfLead+i_s-j,j_IntFer))**2, j=0,N_smth)
   !write(2,*) 't_up',  (Abs(CTime_t(IntfLead+i_s+j,j_IntFer))**2, j=0,N_smth)
   !write(2,*) 't_dwn', (Abs(CTime_t(IntfLead+i_s-j,j_IntFer))**2, j=0,N_smth)
   write(2,*) 'p/t_up',  (CTime_p(IntfLead+i_s+j,j_IntFer)/CTime_t(IntfLead+i_s+j,j_IntFer), j=0,N_smth)
   write(2,*) 'p/t_dwn', (CTime_p(IntfLead+i_s-j,j_IntFer)/CTime_t(IntfLead+i_s-j,j_IntFer), j=0,N_smth)
   !
   If(i_slice.ge.100) TestCh2=.false.
   If(TestCh2) then
      Allocate (D_a(1:Nr_IntFer), Ph_a(1:Nr_IntFer) )
      Allocate (wap_PB(1:3,1:Nr_IntFer), wat_PB(1:3,1:Nr_IntFer) )
      Allocate (wEnu_p(0:IntfNuDim), wEnu_t(0:IntfNuDim) )
      Allocate (wETime_p(1:2*IntfNuDim), wETime_t(1:2*IntfNuDim) )
      Allocate (wEtime_ap(-N_smth:N_smth,1:Nr_IntFer), wEtime_at(-N_smth:N_smth,1:Nr_IntFer) )
      Allocate (dChi_ap(1:Nr_IntFer), dChi_at(1:Nr_IntFer) )
      Allocate (W_ap(1:Nr_IntFer), W_at(1:Nr_IntFer) )
      dChi_ap(:)=0.  ;     dChi_at(:)=0.
      wEnu_p(:)=0.   ;     wEnu_t(:)=0.
      SumDiff=0      ;     SumSq=0.
   EndIf
   !
   ! Basis for Stokes parameters
!   D=SUM(PixLoc(:)*PixLoc(:))  ; HVRBasis(:,3)=PixLoc(:)/sqrt(D)  ! radial, out
!   HVRBasis(1,1)=-PixLoc(2)  ; HVRBasis(2,1)=PixLoc(1) ; HVRBasis(3,1)=0.  ! horizontal
!   HorDist=HVRBasis(1,1)*HVRBasis(1,1)+ HVRBasis(2,1)*HVRBasis(2,1)  ;  HVRBasis(:,1)=HVRBasis(:,1)/sqrt(HorDist)
!   HVRBasis(1,2)=HVRBasis(2,3)*HVRBasis(3,1)- HVRBasis(3,3)*HVRBasis(2,1) ! vertical
!   HVRBasis(2,2)=HVRBasis(3,3)*HVRBasis(1,1)- HVRBasis(1,3)*HVRBasis(3,1)
!   HVRBasis(3,2)=HVRBasis(1,3)*HVRBasis(2,1)- HVRBasis(2,3)*HVRBasis(1,1)
!   NEh2VHR: VHRBasis(NEh,VHR)
   D=SUM(PixLoc(:)*PixLoc(:))  ; NEh2PB(:,3)=PixLoc(:)/sqrt(D)  ! radial, out (North,E=0,h)=> VHR(3)~(N,0,h)
   NEh2PB(1,2)=-PixLoc(2)  ; NEh2PB(2,2)=PixLoc(1) ; NEh2PB(3,2)=0.  ! horizontal  VHR(2)~(0,N,0) to the right
   HorDist=NEh2PB(1,2)*NEh2PB(1,2)+ NEh2PB(2,2)*NEh2PB(2,2)  ;  NEh2PB(:,2)=NEh2PB(:,2)/sqrt(HorDist)
   NEh2PB(1,1)=NEh2PB(2,3)*NEh2PB(3,2)- NEh2PB(3,3)*NEh2PB(2,2) ! vertical tilted toward
   NEh2PB(2,1)=NEh2PB(3,3)*NEh2PB(1,2)- NEh2PB(1,3)*NEh2PB(3,2)  ! VHR(1)~(-h,0,N)
   NEh2PB(3,1)=NEh2PB(1,3)*NEh2PB(2,2)- NEh2PB(2,3)*NEh2PB(1,2)
   !
   Fnu_PB(:,:)=0.
   Cur2E(:,:)=0.
   Esq_ak=0.
   Power_p(0)=Power_p(0)+1.
   Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
      !  Get interferometry phase shifts
      i_ant=IntFer_ant(j_IntFer)
      Call RelDist(PixLoc(1),Ant_pos(1,i_ant,i_chunk),RDist)
      dt_AntPix =Rdist - Ant_RawSourceDist(i_ant,i_chunk)
      dphase = exp(-ipi*dt_AntPix /IntfNuDim)
      Phase =exp(-ipi*dt_AntPix *inu1/IntfNuDim)
      !
      !write(2,*) i_ant, j_IntFer, Ant_pos(:,i_ant,i_chunk)
      Ras(1)=(PixLoc(1)-Ant_pos(1,i_ant,i_chunk))/1000.  ! \vec{R}_{antenna to source}
      Ras(2)=(PixLoc(2)-Ant_pos(2,i_ant,i_chunk))/1000.
      Ras(3)=(PixLoc(3)-Ant_pos(3,i_ant,i_chunk))/1000.
      HorDist= Ras(1)*Ras(1) + Ras(2)*Ras(2)  ! =HYPOT(X,Y)
      D=sqrt(HorDist + Ras(3)*Ras(3))
      !AntSourceD =D
      HorDist=sqrt( HorDist ) ! =HYPOT(X,Y)
      Thet_r=atan(HorDist,Ras(3))  ! Zenith angle
      Phi_r=atan2( Ras(2) , Ras(1) )  ! \phi=0 = north
      Thet_d =Thet_r*180/pi
      Phi_d =Phi_r*180/pi
      !
      Call AntFun_Inv(thet_d ,Phi_d ) ! sets ,Ji_p0,Ji_t0,Ji_p1,Ji_t1; Inverse Jones; gives _p & _t polarized fields when contracted with (0.1) voltages
      Vec_p(1)=sin(Phi_r)              ; Vec_p(2)=-cos(Phi_r)        ; Vec_p(3)=0.
      Vec_t(1)=-cos(Thet_r)*Vec_p(2)  ; Vec_t(2)=cos(Thet_r)*Vec_p(1) ; Vec_t(3)=-sin(Thet_r)
      !
      Do i=1,3  ! convert t&p orientations from (NEh) to polar base
         p_PB(i)=SUM( NEh2PB(:,i)*Vec_p(:) )
         t_PB(i)=SUM( NEh2PB(:,i)*Vec_t(:) )
      Enddo
      !write(2,"(4(A,3F7.3),3F7.3)") 'p_PB=',p_PB(:),', t_PB=',t_PB(:)! &
      !
      ! determine signal strength in slice
      i_s=1+i_slice*N_smth+NINT(dt_AntPix )  ! approximately correct, upto rounding errors for dt
      Pow_p=0.
      Pow_t=0.
      Do j=-N_smth,N_smth ! fold with smoothing function
         Pow_p= Pow_p + smooth(j)*Abs(CTime_p(IntfLead+i_s+j,j_IntFer))**2
         Pow_t= Pow_t + smooth(j)*Abs(CTime_t(IntfLead+i_s+j,j_IntFer))**2
      Enddo
      Power_p(j_IntFer)=Power_p(j_IntFer)+Pow_p
      Power_t(j_IntFer)=Power_t(j_IntFer)+Pow_t
      !write(2,*) 'pow=',i_slice,j_IntFer,Power_p(j_IntFer)/Power_p(0), Noise_p(j_IntFer), &
      !            Power_t(j_IntFer)/Power_p(0), Noise_t(j_IntFer)
      w=I_Scale ! to have alpha of order unity
      w_p=W/(Noise_p(j_IntFer)*sqrt(Pow_p/Noise_p(j_IntFer)+1.)) ! the denominator should be due to 'experimental accuracy' in power
      w_t=W/(Noise_t(j_IntFer)*sqrt(Pow_t/Noise_t(j_IntFer)+1.))
      ! For chi-square calculation, sum_ak[ W_ak * E_ak^2]
      Esq_ak=Esq_ak + (Pow_p*w_p + Pow_t*w_t)  !
      If(TestCh2) Then
         D_a(j_IntFer)=D
         Ph_a(j_IntFer)=Phi_d
         sw_p=sqrt(w_p/I_Scale)
         sw_t=sqrt(w_t/I_Scale)
         wap_PB(:,j_IntFer)=p_PB(:)*sw_p/D
         W_ap(j_IntFer)=w_p/I_Scale
         wat_PB(:,j_IntFer)=t_PB(:)*sw_t/D
         W_at(j_IntFer)=w_t/I_Scale
      EndIf
      !write(2,*) 'Esq_ak=',Esq_ak
      !
      Do i_nu=inu1,inu2   ! Increment frequency spectrum with the spectrum from this antenna
         nu=i_nu*dnu
         i_freq=Int(nu)
         dfreq=nu-i_freq
         Sp=((1.-dfreq)*Ji_p0(i_freq) + dfreq*Ji_p0(i_freq+1)) *Cnu0(i_nu,j_IntFer) + &
            ((1.-dfreq)*Ji_p1(i_freq) + dfreq*Ji_p1(i_freq+1)) *Cnu1(i_nu,j_IntFer)
         St=((1.-dfreq)*Ji_t0(i_freq) + dfreq*Ji_t0(i_freq+1)) *Cnu0(i_nu,j_IntFer) + &
            ((1.-dfreq)*Ji_t1(i_freq) + dfreq*Ji_t1(i_freq+1)) *Cnu1(i_nu,j_IntFer)
         Fnu_PB(i_nu,:)=Fnu_PB(i_nu,:) + phase *( p_PB(:)* Sp*w_p + t_PB(:)* St*w_t )*Gain(i_freq)/D
         If(TestCh2) Then
            wEnu_p(i_nu)= wEnu_p(i_nu)+ phase * Sp*sw_p *Gain(i_freq)
            wEnu_t(i_nu)= wEnu_t(i_nu)+ phase * St*sw_t *Gain(i_freq)
         EndIf
         phase =phase *dphase
      Enddo
      !
      ! Calculate contribution to Cur2E == A matrix
      Do i=1,3
         Do j=1,3  ! a single factor distance since the other is already included in W (as in W_p & W_t)
            Cur2E(i,j)=Cur2E(i,j)+ ( p_PB(i)*p_PB(j)*w_p + t_PB(i)*t_PB(j)*w_t)/(D*D)
         Enddo
      Enddo
      !
      If(TestCh2) Then
         Call RFTransform_CF2CT(wEnu_p(0),wETime_p(1) )
         Call RFTransform_CF2CT(wEnu_t(0),wETime_t(1) )
         wEnu_p(:)=0.
         wEnu_t(:)=0.
         Do j=-N_smth,N_smth ! copy
            wEtime_ap(j,j_IntFer)=wETime_P(IntfLead+i_s+j)
            wEtime_at(j,j_IntFer)=wETime_t(IntfLead+i_s+j)
         Enddo
      EndIf
      !
   Enddo     ! j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
   !
   Call matinv3(Cur2E, Ai)
   Do m=1,3  ! convert to time
      Call RFTransform_CF2CT(Fnu_PB(0,m),FTime_PB(1,m) )
   Enddo
   !
   i_s=1+i_slice*N_smth
   Stk(:,:)=0.
   FdotI=0.
   If(TestCh2) then
      i=100   !            <<<<<<<<<<<<<<<<<<<<<<<<
      write(txt,"(I2.2)") i_slice
      OPEN(UNIT=30,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'DelChi2_'//TRIM(txt)//'.dat')
      write(30,"(3I4,I7,' 0')") i_slice, Nr_IntFer, i, IntfBase+IntfLead+i_s
      Write(30,"(A,9x,A,13x,';   ', 9(A,'/I [%]   ;   '),A )") &
         '! i_ant, j; Ep-data;  model , Et-data, model, d_chi_p','d_chi_t'
   EndIf
   Do j=-N_smth,N_smth ! fold time-dependence with smoothing function
      Do m=1,3  ! convert to time
         AiF(m)=SUM( Ai(m,:)*FTime_PB(IntfLead+i_s+j,:) )
      Enddo
      FdotI=FdotI + smooth(j)*Real( SUM( FTime_PB(IntfLead+i_s+j,:)*Conjg(AiF) ) )  ! Imag is zero
      Do m=1,3
         Do n=1,3
            Stk(m,n)= Stk(m,n) + smooth(j)*AiF(m)*Conjg(AiF(n))
         Enddo
      Enddo
      If(TestCh2) then
         Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
            Del_p = ABS( wEtime_ap(j,j_IntFer) -SUM( wap_PB(:,j_IntFer)*AiF(:) ) )
            Del_t = ABS( wEtime_at(j,j_IntFer) -SUM( wat_PB(:,j_IntFer)*AiF(:) ) )
            dChi_ap(j_IntFer)=dChi_ap(j_IntFer) + smooth(j)* Del_p**2 ! wEtime_ap(j,j_IntFer) * Conjg( Del_p  ) ! Imag is zero ??
            dChi_at(j_IntFer)=dChi_at(j_IntFer) + smooth(j)* Del_t**2 ! wEtime_at(j,j_IntFer) * Conjg( Del_t  ) ! Imag is zero ??
            ! SumDiff=SumDiff+ sqrt(smooth(j))*(Del_p + Del_t)
            SumSq=SumSq + smooth(j)* (Del_p**2 + Del_t**2)
         Enddo
         j_IntFer=i
         write(30,"(I4,I5,10G13.3)")  j_IntFer, j, &
            wEtime_ap(j,j_IntFer), SUM( wap_PB(:,j_IntFer)*AiF(:) ), &
            wEtime_at(j,j_IntFer), SUM( wat_PB(:,j_IntFer)*AiF(:) ), &
            dChi_ap(j_IntFer), dChi_at(j_IntFer)
      EndIf
   Enddo
   If(TestCh2) then
      Close(UNIT=30)
   EndIf
   !
   Chi2pDF=(Esq_ak-FdotI)/(2.*Nr_IntFer)/ I_Scale   ! undo arbitrary scaling factor
   !Write(2,*) Esq_ak,FdotI,'chi-sq=', Esq_ak-FdotI, ', /DegFree=', (Esq_ak-FdotI)/(2.*Nr_IntFer), I_Scale
   !
   StI = Real(Stk(1,1) + Stk(2,2) + Stk(3,3))
   StI12=Real(Stk(1,1) + Stk(2,2))
   StI3 =  Real(Stk(3,3))
   StQ = Real(Stk(1,1)-Stk(2,2))
   StU = 2*Real(Stk(1,2))
   StV =  2*Imag(Stk(1,2))
   StU2 = 2*Real(Stk(3,1))
   StV2 =  2*Imag(Stk(3,1))
   StU1 = 2*Real(Stk(2,3))
   StV1 =  2*Imag(Stk(2,3))
   W=ATAN2(StU,StQ)/2.*180./pi
   If(W.lt.-90.) W=W+180.
   write(2,"(I4,A,G12.4,7(A,F7.2))") i_slice,': Wgt: I123=',StI,', I12/I=',StI12/StI, ', Q/I=',StQ/StI, &
         ', U/I=',StU/StI, ', V/I=',StV/StI, ', I3/I=',StI3/StI,', angle=',W,', chi^2=',Chi2pDF
   Do m=1,3
      Do n=1,3
         Cur2E(m,n)=SUM( Ai(m,:)*Ai(:,n) )*Chi2pDF**2  /4. ! this should correspond to the square of the error
         ! arbitrary factor 100 included to stretch error bars and make them visible
      EndDo
    !  write(2,*) 'Correlation matrix(m,:)=',Cur2E(m,:)
   EndDo
   dStI = SQRT(SUM(Cur2E(:,:) ))
   dStI12=SQRT(Cur2E(1,1) + Cur2E(2,2))
   dStI3 =  SQRT( Cur2E(3,3) )
   dStQ = dStI12
   dStU = 2*SQRT( ABS(Cur2E(1,2)) )
   dStV =  dStU
   dStU2 = 2*SQRT( ABS(Cur2E(3,1)) )
   dStV2 =  dStU2
   dStU1 = 2*SQRT( ABS(Cur2E(2,3)) )
   dStV1 =  dStU1
   !Write(2,*) Esq_ak,FdotI,'chi-sq=', Esq_ak-FdotI, ', /DegFree=', (Esq_ak-FdotI)/(2.*Nr_IntFer), I_Scale
   !write(2,"(A, 3G12.2,A,2F7.2,I7)") 'Variance(Stokes), [A(m,:)^-1 x chi^2/DoF]^2=',Cur2E(1,1), Cur2E(2,2), Cur2E(3,3),&
   !         '; Chi^2/DoF=',Chi2pDF, SumSq/(2.*Nr_IntFer), IntfBase+IntfLead+i_s
   !
   If(TestCh2) then
      SumSq=SumSq/(2.*Nr_IntFer)
      write(txt,"(I2.2)") i_slice
      OPEN(UNIT=30,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'IntfSpecChi2_'//TRIM(txt)//'.dat')
      write(30,"(2I4,F9.2,2G12.3,I7,' 0')") i_slice, Nr_IntFer, SumSq, StI12, StI, IntfBase+IntfLead+i_s
      Write(30,"(A,9x,A,13x,';   ', 9(A,'/I [%]   ;   '),A )") &
         '! nr, i_ant, Ant#, St#;  D[km] , phi[deg], d_chi_ph, d_chi_th','StI','StI12'
      w=(W_ap(1)+W_at(1))/2.
      Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
         i_ant=IntFer_ant(j_IntFer)
         write(30,"(I4,I5,I6,I6,F9.3,F8.1,5G12.3)")  j_IntFer, i_ant, &
            Ant_IDs(I_ant,i_chunk), Ant_Stations(I_ant,i_chunk), D_a(j_IntFer), Ph_a(j_IntFer), &
            dChi_ap(j_IntFer), dChi_at(j_IntFer), &
            W_ap(j_IntFer)/W, W_at(j_IntFer)/W
      Enddo
      Close(Unit=30)
      !write(2,*) 'dChi_ap:', dChi_ap(:)
      !write(2,*) 'dChi_at:', dChi_at(:)
      DeAllocate (wap_PB, wat_PB, wEnu_p, wEnu_t, wETime_p, wETime_t, wEtime_ap, wEtime_at )
      DeAllocate (dChi_ap, dChi_at )
!call GLE -d pdf -o NLMx_dChi2_01.pdf %UtilDir%EI_chisq.gle "../19A-5/files/NLtst" "01"
      Call GLEplotControl(PlotType='EI_delta-chisq', PlotName='EI_delchisq_'//TRIM(txt)//TRIM(OutFileLabel), &
            PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel)//' "'//TRIM(txt)//'"' ) ! does NOT work in windows
      Call GLEplotControl(PlotType='EI_chisq', PlotName='EI_chisq_'//TRIM(txt)//TRIM(OutFileLabel), &
            PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel)//' "'//TRIM(txt)//'"' ) ! does NOT work in windows
      Call GLEplotControl(SpecialCmnd='rm '//TRIM(DataFolder)//TRIM(OutFileLabel)//'IntfSpecChi2_'//TRIM(txt)//'.dat') !  Command to delete files
      Call GLEplotControl(SpecialCmnd='rm '//TRIM(DataFolder)//TRIM(OutFileLabel)//'DelChi2_'//TRIM(txt)//'.dat') !  Command to delete files
   EndIf
   Return
End Subroutine EI_PolarizW
!-----------------------------------------------

 Subroutine matinv3(A,B)
   use constants, only : dp
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    !complex(wp), intent(in) :: A(3,3)   !! Matrix
    !complex(wp)             :: B(3,3)   !! Inverse matrix
    !complex(wp)             :: detinv
    Real(dp), intent(in) :: A(3,3)   !! Matrix
    Real(dp)             :: B(3,3)   !! Inverse matrix
    Real(dp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    Return
End Subroutine matinv3
