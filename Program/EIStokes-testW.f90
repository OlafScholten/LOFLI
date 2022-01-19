    !
!-----------------------------------------------
Subroutine EI_PolarizSlice(i_slice)
   use constants, only : dp
   use DataConstants, only : Ant_nrMax
   use Interferom_Pars, only :  Nr_IntFerMx, Nr_IntferCh ! the latter gives # per chunk
   use Interferom_Pars, only : W_ap, W_at, Cnu0, Cnu1, IntfNuDim, CTime_p, CTime_t, Noise_p, Noise_t, AntPeak_OffSt
   Use Interferom_Pars, only : IntfBase, IntfLead, SumWindw, N_Smth
   Use Interferom_Pars, only : i_chunk, PixLoc, CenLoc
   Implicit none
   Integer, intent(in) :: i_slice
   Integer :: i_sample, i_1, i_2, Output, j_IntFer
   Real(dp) :: ChiSq, DelChi(-N_Smth:+N_Smth,1:Nr_IntFerMx)
   Real(dp) :: FitDelay(1:Ant_nrMax), VoxLoc(1:3), del_1, del_2
   Character(len=8) :: Label
   logical :: First=.true.
   If(First) Then
      First=.false.
      Write(2,*) ' IntfBase, IntfNuDim:', IntfBase, IntfNuDim, Nr_IntFerMx
      Call EI_PolSetUp(Nr_IntFerCh(i_chunk), IntfBase, i_chunk, CenLoc(:), &
         AntPeak_OffSt(1,1), Cnu0(0,1,1), Cnu1(0,1,1))
      j_IntFer=1
      write(Label,"(i2.2)") j_IntFer
      Call TimeTracePlot(j_IntFer, IntfBase, i_chunk, CenLoc, SumWindw, Label) ! window width around central sample
   EndIf
   !
   Output=2
   If(SumWindw/N_Smth.ge.5) Output=1
   i_sample=IntfLead + 1+i_slice*N_smth
   FitDelay(:)=0.
   write(Label,"('Slc ',i4.2)") i_slice
   write(2,"(1x,A,i4,A,I5,A,2(F9.4,','),F9.4,A)") 'Slice',i_slice,', Ref. ant. sample=',IntfBase+i_sample, &
      ', Max Intensity @ (N,E,h)=(',PixLoc(:)/1000.,') [km]'
   Call EI_Weights(Nr_IntFerCh(i_chunk), i_sample, i_chunk, PixLoc(:), AntPeak_OffSt(1,1), W_ap, W_at)
   !write(2,*) ' i_sample, i_chunk, :', i_sample, i_chunk, Nr_IntFerCh(:),'fitdelay:',FitDelay
   Call EI_PolGridDel(Nr_IntFerCh(i_chunk), FitDelay, i_sample, i_chunk, PixLoc(:), &
      AntPeak_OffSt(1,1), W_ap(1,1), W_at(1,1), Cnu0(0,1,1), Cnu1(0,1,1), Output, DelChi, Label)
   !
   Return
   If(SumWindw/N_Smth.ge.5) Return
   !i_sample=IntfLead + i_slice*N_smth  ! approximately correct, upto rounding errors for dt
   Del_1=2  ! in meter
   Del_2=2
   Output=1
      write(2,*) i_slice, i_sample,'Del_1=',Del_1,', Del_2=',Del_2, ' [m]'
   Do i_1=-2,2
      Do i_2=-2,2
         VoxLoc(1)=PixLoc(1) + i_1*del_1
         VoxLoc(2)=PixLoc(2) + i_2*del_2
         VoxLoc(3)=PixLoc(3)
         write(Label,"('Gr ',I2,',',i2)") i_1,i_2
         Call EI_Weights(Nr_IntFerCh(i_chunk), i_sample, i_chunk, VoxLoc(:), AntPeak_OffSt(1,1), W_ap, W_at)
         Call EI_PolGridDel(Nr_IntFerCh(i_chunk), FitDelay, i_sample, i_chunk, VoxLoc(:), &
            AntPeak_OffSt(1,1), W_ap(1,1), W_at(1,1), Cnu0(0,1,1), Cnu1(0,1,1), Output, DelChi, Label)
      Enddo
   Enddo
   !
   Return
End Subroutine EI_PolarizSlice
!-----------------------------------------------
Subroutine EI_PolarizW(Nr_IntFer, IntfNuDim, i_slice)
   ! Needs: Ant_RawSourceDist as used when reading in data (when fitted time-shifts = 0)
   !alculates: 1) trace in frequency space for even and odd antennas symmetrically around
   !  calculated pulse time (based on suggested source time and position) for each antenna pair
   !           2) noise background level for each antenna, needed for chi^2 calculation
   !Purpose: to be used in interferometry calculation for a voxel in the vicinity of the
   !  guessed source and with fitted antenna timings
   !--------------------------------------------
   use constants, only : dp, pi, ci
   use DataConstants, only : Time_Dim, Cnu_dim
   use DataConstants, only : DataFolder, OutFileLabel
   use Chunk_AntInfo, only : CTime_spectr, Ant_Stations, Ant_IDs, Ant_nr, Ant_pos, Ant_RawSourceDist
   use Chunk_AntInfo, only : Start_time! , TimeBase
   use Chunk_AntInfo, only : NormOdd, NormEven ! Powr_eo,NAnt_eo
   use FFT, only : RFTransform_CF, RFTransform_CF2CT
   Use Interferom_Pars, only : IntfBase, IntfLead, SumWindw, N_smth, smooth, NrPixSmPowTr
   Use Interferom_Pars, only : IntFer_ant
   Use Interferom_Pars, only : i_chunk, PixLoc, CenLoc
   Use Interferom_Pars, only : StI, StI12, StQ, StU, StV, StI3, StU1, StV1, StU2, StV2, P_un, P_lin, P_circ
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
   Real(dp) :: RDist, dt_AntPix
   Complex(dp), Allocatable, save :: CTime_p(:,:), CTime_t(:,:)
   Complex(dp), Allocatable, save :: Cnu0(:,:), Cnu1(:,:)
   !Complex(dp), save :: Cnu0(0:IntfNuDim,Nr_IntFer), Cnu1(0:IntfNuDim,Nr_IntFer)
   Real(dp), Allocatable, save :: Noise_p(:), Noise_t(:)
   Real(dp) :: Pow_p, Pow_t
   Real(dp), save :: I_scale
   Real(dp) :: AMat(1:3,1:3), Ai(1:3,1:3)
   Complex :: Stk(3,3), AiF(1:3)
   logical :: First=.true.
   Real(dp) :: Esq_ak, FdotI
   Real(dp) :: St8
   !
   logical :: TestCh2
   Real(dp), Allocatable :: D_a(:),Ph_a(:)  ! weighted p and t directions in polarization basis (_PB)
   Real(dp), Allocatable :: wap_PB(:,:),wat_PB(:,:)  ! weighted p and t directions in polarization basis (_PB)
   Complex(dp), Allocatable :: wEnu_p(:), wEnu_t(:)
   Complex(dp), Allocatable :: wETime_p(:), wETime_t(:)
   Complex(dp), Allocatable :: wEtime_ap(:,:), wEtime_at(:,:) ! weighted E fields, time dependent
   !Complex(dp), Allocatable :: dChi_ap(:), dChi_at(:)
   Real(dp), Allocatable :: dChi_ap(:), dChi_at(:)
   Real(dp) :: W_ap(1:Nr_IntFer), W_at(1:Nr_IntFer)
   Real(dp) :: SumDiff, SumSq, Del_p, Del_t, W_p, W_t, sw_p, sw_t
   Integer :: i_1, i_2
   real(dp) :: Del_1, Del_2, VoxLoc(1:3)
   character(len=2) :: txt
   !
   TestCh2=.true.
   IntfDim=2*IntfNuDim
   If(First) Then  ! preparation
      First=.false.
      If(NrPixSmPowTr.ge.10) TestCh2=.false.
      Allocate( CTime_p(1:2*IntfNuDim,Nr_IntFer), CTime_t(1:2*IntfNuDim,Nr_IntFer) )
      Allocate( Cnu0(0:IntfNuDim,Nr_IntFer), Cnu1(0:IntfNuDim,Nr_IntFer) )
      Allocate( Noise_p(1:Nr_IntFer), Noise_t(1:Nr_IntFer) )
      !Allocate( W_p(0:Nr_IntFer), W_t(0:Nr_IntFer) )
      !
      dnu=100./IntfNuDim   ! [MHz] Jones matrix is stored on 1MHz grid
      inu1=Int(Freq_min/dnu)+1
      inu2=Int(Freq_max/dnu)-1
      !PowerScale=1.
      !NormEven=sqrt(PowerScale*2.*Powr_eo(0)/(Powr_eo(0)+Powr_eo(1)) )/100.
      !NormOdd=sqrt(PowerScale*2.*Powr_eo(1)/(Powr_eo(0)+Powr_eo(1)) )/100.
      Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas to detremine t- and p- polarized trace for central pixel; needed for noise estimate
         i_ant=IntFer_ant(j_IntFer,i_chunk)
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
         RTime(:)=REAL(CTime_spectr(IntfBase+1:IntfBase+IntfDim,i_ant,i_chunk))*NormEven  ! IntfBase=SumStrt - IntfLead
         Call RFTransform_CF(RTime,Cnu0(0,j_IntFer))
         RTime(:)=REAL(CTime_spectr(IntfBase+1:IntfBase+IntfDim,i_ant+1,i_chunk))*NormOdd
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
         Call RFTransform_CF2CT(nu_p(0),CTime_p(1,j_IntFer) )  ! needed for the peak/background estimate
         Call RFTransform_CF2CT(nu_t(0),CTime_t(1,j_IntFer) )
      Enddo    !  j_IntFer=1,Nr_IntFer
      !
    !write(2,*) 'Noise_p:',Noise_p
    !write(2,*) 'Noise_t:',Noise_t
     ! Just for checking causality on the first pass
      j_IntFer=1
      i_ant=IntFer_ant(j_IntFer,i_chunk)
      Call RelDist(PixLoc(1),Ant_pos(1,i_ant,i_chunk),RDist)
      dt_AntPix =Rdist - Ant_RawSourceDist(i_ant,i_chunk)
      i_s=1+NINT(dt_AntPix )  ! approximately correct, upto rounding errors for dt
      write(txt,"(I2.2)") j_IntFer
      write(2,*) 'E-field hilbert envelope for the central pixel for the antenna pair ',TRIM(txt),' written to:',&
            trim(DataFolder)//TRIM(OutFileLabel)//'Trace_'//TRIM(txt)//'.dat',' shift=',NINT(dt_AntPix ),' samples'
      OPEN(UNIT=30,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'Trace_'//TRIM(txt)//'.dat')
      Write(30,*) i_s,i_s+SumWindw
      Write(2,*) 'In EI_PolarizW, data written to file:', trim(DataFolder)//TRIM(OutFileLabel)//'Trace_'//TRIM(txt)//'.dat'
      write(2,*) 'Header: i_s,i_s+SumWindw; Data range=',i_s,i_s+SumWindw
      Write(2,*) i_s,i_s+SumWindw
      Do i=i_s,i_s+SumWindw
         write(30,"(i5,40G12.4)") i,ABS(CTime_spectr(IntfBase+IntfLead+i,i_ant,i_chunk))*NormEven &
            ,ABS(CTime_spectr(IntfBase+IntfLead+i,i_ant+1,i_chunk))*NormOdd &
            ,Abs(CTime_p(IntfLead+i,j_IntFer)), Abs(CTime_t(IntfLead+i,j_IntFer))
      EndDo
      Close(Unit=30)
      Call GLEplotControl(PlotType='PlotTimeTraces', PlotName='TimeTraces_'//TRIM(txt)//TRIM(OutFileLabel), &
            PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel)//' "'//TRIM(txt)//'"' ) ! does NOT work in windows
      !i_s=1+i_slice*N_smth+NINT(dt_AntPix )  ! approximately correct, upto rounding errors for dt
      !write(2,"(A,40G12.4)") 'p_dwn', (Abs(CTime_p(IntfLead+i_s-j,j_IntFer))**2, j=0,N_smth/2+1)
      !write(2,"(A,40G12.4)") 'p_up',  (Abs(CTime_p(IntfLead+i_s+j,j_IntFer))**2, j=0,N_smth/2+1)
      !write(2,"(A,40G12.4)") 't_dwn', (Abs(CTime_t(IntfLead+i_s-j,j_IntFer))**2, j=0,N_smth/2+1)
      !write(2,"(A,40G12.4)") 't_up',  (Abs(CTime_t(IntfLead+i_s+j,j_IntFer))**2, j=0,N_smth/2+1)
      !write(2,"(A,40G12.4)") 'p/t_dwn', (CTime_p(IntfLead+i_s-j,j_IntFer)/CTime_t(IntfLead+i_s-j,j_IntFer), j=0,N_smth/2+1)
      !write(2,"(A,40G12.4)") 'p/t_up',  (CTime_p(IntfLead+i_s+j,j_IntFer)/CTime_t(IntfLead+i_s+j,j_IntFer), j=0,N_smth/2+1)
      !
   EndIf  ! End preparation
   !
   If(i_slice.ge.10) TestCh2=.false.
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
   i_s=1+i_slice*N_smth
   !Write(2,*) 'CTime_p(i_s+j,j_IntFer):',IntfLead+i_s,CTime_p(IntfLead+i_s,:)
   !Write(2,*) 'CTime_t(i_s+j,j_IntFer):',i_s,CTime_t(IntfLead+i_s,:)
   Esq_ak=0.
   !write(2,*) PixLoc(:), NEh2PB
   Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
      !  Get interferometry phase shifts
      i_ant=IntFer_ant(j_IntFer,i_chunk)
      Call RelDist(PixLoc(1),Ant_pos(1,i_ant,i_chunk),RDist)
      dt_AntPix =Rdist - Ant_RawSourceDist(i_ant,i_chunk)
      i_s=1+i_slice*N_smth+NINT(dt_AntPix )  ! approximately correct, upto rounding errors for dt
      Pow_p=0.
      Pow_t=0.
      Do j=-N_smth,N_smth ! pulse power, fold with smoothing function
         Pow_p= Pow_p + smooth(j)*Abs(CTime_p(IntfLead+i_s+j,j_IntFer))**2
         Pow_t= Pow_t + smooth(j)*Abs(CTime_t(IntfLead+i_s+j,j_IntFer))**2
      Enddo
      W_ap(j_IntFer)=1./(Noise_p(j_IntFer)*sqrt(Pow_p/Noise_p(j_IntFer)+1.)) ! the denominator should be due to 'experimental accuracy' in power
      W_at(j_IntFer)=1./(Noise_t(j_IntFer)*sqrt(Pow_t/Noise_t(j_IntFer)+1.))
      !w_p=I_Scale*W_ap(j_IntFer)  ! to have alpha of order unity
      !w_t=I_Scale*W_at(j_IntFer)
      ! For chi-square calculation, sum_ak[ W_ak * E_ak^2] :
      Esq_ak=Esq_ak + (Pow_p*W_ap(j_IntFer) + Pow_t*W_at(j_IntFer))  ! Used in approximation to chi^2
   Enddo
   !write(2,*) 'w_ap',w_ap
   !write(2,*) 'w_at',w_at   !
   Call EI_PolGrid(Nr_IntFer, IntfNuDim, i_slice, PixLoc, W_ap, W_at, Cnu0, Cnu1, NEh2PB, TestCh2)

   !stop


   Return
   !
   If(i_slice.ge.10) Return
   i_s=1+i_slice*N_smth  ! approximately correct, upto rounding errors for dt
   Del_1=2
   Del_2=2
   TestCh2=.false.
   Do i_1=-2,2
      write(2,*) i_slice, i_s,'i_1=',i_1,', i_2=-2,+2'
      Do i_2=-2,2
         VoxLoc(1)=PixLoc(1) + i_1*del_1
         VoxLoc(2)=PixLoc(2) + i_2*del_2
         VoxLoc(3)=PixLoc(3)
         Call EI_PolGrid(Nr_IntFer, IntfNuDim, i_slice, VoxLoc, W_ap, W_at, Cnu0, Cnu1, NEh2PB, TestCh2)
      Enddo
   Enddo
   !
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
!----------------------------------------------
!-----------------------------------------------
Subroutine EI_PolGrid(Nr_IntFer, IntfNuDim, i_slice, VoxLoc, W_ap, W_at, Cnu0, Cnu1, NEh2PB, TestCh2)
   ! Needs: Ant_RawSourceDist as used when reading in data (when fitted time-shifts = 0)
   !alculates: 1) trace in frequency space for even and odd antennas symmetrically around
   !  calculated pulse time (based on suggested source time and position) for each antenna pair
   !           2) noise background level for each antenna, needed for chi^2 calculation
   !Purpose: to be used in interferometry calculation for a voxel in the vicinity of the
   !  guessed source and with fitted antenna timings
   !--------------------------------------------
   use constants, only : dp, pi, ci
   use DataConstants, only : Time_Dim, Cnu_dim
   use DataConstants, only : DataFolder, OutFileLabel
   use Chunk_AntInfo, only : CTime_spectr, Ant_Stations, Ant_IDs, Ant_nr, Ant_pos, Ant_RawSourceDist
   use Chunk_AntInfo, only : Start_time! , TimeBase
   !use Chunk_AntInfo, only : Powr_eo,NAnt_eo
   use FFT, only : RFTransform_CF, RFTransform_CF2CT
   Use Interferom_Pars, only : IntfBase, IntfLead, SumWindw, N_smth, smooth, NrPixSmPowTr
   Use Interferom_Pars, only : IntFer_ant
   Use Interferom_Pars, only : i_chunk
   Use Interferom_Pars, only : StI, StI12, StQ, StU, StV, StI3, StU1, StV1, StU2, StV2, P_un, P_lin, P_circ
   Use Interferom_Pars, only : dStI, dStI12, dStQ, dStU, dStV, dStI3, dStU1, dStV1, dStU2, dStV2, Chi2pDF
   use AntFunCconst, only : Freq_min, Freq_max,Ji_p0,Ji_t0,Ji_p1,Ji_t1, Gain  !J_0p,J_0t,J_1p,J_1t,
   use GLEplots, only : GLEplotControl
   Implicit none
   Integer, intent(in) :: Nr_IntFer, IntfNuDim, i_slice
   Real(dp), intent(in) :: VoxLoc(1:3), W_ap(1:Nr_IntFer), W_at(1:Nr_IntFer), NEh2PB(1:3,1:3)
   Complex(dp), intent(in) :: Cnu0(0:IntfNuDim,Nr_IntFer), Cnu1(0:IntfNuDim,Nr_IntFer)
   logical, intent(in) :: TestCh2
   complex(dp), parameter :: ipi=ci*pi
   integer :: i_ant, j_IntFer, i, j, i_freq, i_nu, IntfDim, inu1, inu2, m, n, i_s
   Real(dp) :: Vec_p(1:3), Vec_t(1:3), Ras(1:3),  p_PB(1:3), t_PB(1:3)
   !Real(dp) :: NEh2PB(1:3,1:3) ! NEh to Polarization Basis; NEh2VHR(NEh,VHR); formerly known as VHRBasis(1:3,1:3)
   Real(dp) :: HorDist, Thet_r ,Phi_r, dfreq, D, W, dnu, nu
   Real(dp) :: thet_d, Phi_d  ! AntSourceD(1:Nr_IntFer),
   Complex(dp) :: Sp, St
   Complex(dp) :: FTime_PB(1:2*IntfNuDim,1:3), Fnu_PB(0:IntfNuDim,1:3)
   Complex(dp) :: Phase, dPhase, nu_p(0:IntfNuDim), nu_t(0:IntfNuDim)
   Real(dp) :: RDist, dt_AntPix
   Real(dp) :: Pow_p, Pow_t, W_p, W_t
   Real(dp) :: AMat(1:3,1:3), Ai(1:3,1:3)
   Complex :: Stk(3,3), AiF(1:3)
   Real(dp) :: Esq_ak, FdotI
   Real(dp) :: St8
   !
   Real(dp) :: sw_p, sw_t, D_a(1:Nr_IntFer), Ph_a(1:Nr_IntFer)   ! just for printing & plotting
   Real(dp) :: wap_PB(1:3,1:Nr_IntFer), wat_PB(1:3,1:Nr_IntFer)  ! weighted p and t directions in polarization basis (_PB)
   Complex(dp) :: wEnu_p(0:IntfNuDim), wEnu_t(0:IntfNuDim)
   Complex(dp) :: wETime_p(1:2*IntfNuDim), wETime_t(1:2*IntfNuDim)
   Complex(dp) :: wEtime_ap(-N_smth:N_smth,1:Nr_IntFer), wEtime_at(-N_smth:N_smth,1:Nr_IntFer) ! weighted E fields, time dependent
   !Complex(dp), Allocatable :: dChi_ap(:), dChi_at(:)
   Real(dp) :: dChi_ap(1:Nr_IntFer), dChi_at(1:Nr_IntFer)
   !Complex(dp) :: SumDiff, SumSq, Del_p, Del_t
   Real(dp) :: SumDiff, SumSq, Del_p, Del_t
   character(len=2) :: txt
   !
   IntfDim=2*IntfNuDim
   dnu=100./IntfNuDim   ! [MHz] Jones matrix is stored on 1MHz grid
   inu1=Int(Freq_min/dnu)+1
   !inu2=Int(Freq_max/dnu-0.5)+1  ! set integration regime to frequencies within filter band
   inu2=Int(Freq_max/dnu)-1
   !
   wEnu_p(:)=0.
   wEnu_t(:)=0.
   i_s=1+i_slice*N_smth
   !write(2,*) i_s
   !flush(unit=2)
   Fnu_PB(:,:)=0.
   AMat(:,:)=0.
   SumSq=0.
   SumDiff=0.
   dChi_ap(:)=0.  ;     dChi_at(:)=0.
   Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
      !  Get interferometry phase shifts
      i_ant=IntFer_ant(j_IntFer,i_chunk)
      Call RelDist(VoxLoc(1),Ant_pos(1,i_ant,i_chunk),RDist)
      dt_AntPix =Rdist - Ant_RawSourceDist(i_ant,i_chunk)
      dphase = exp(-ipi*dt_AntPix /IntfNuDim)
      Phase =exp(-ipi*dt_AntPix *inu1/IntfNuDim)
      !If(j_intFer.eq.110) write(2,*) 'dt_AntPix:',i_ant, j_IntFer, dt_AntPix, Rdist, Ant_RawSourceDist(i_ant,i_chunk)
      !
      !write(2,*) i_ant, j_IntFer, Ant_pos(:,i_ant,i_chunk)
      Ras(1)=(VoxLoc(1)-Ant_pos(1,i_ant,i_chunk))/1000.  ! \vec{R}_{antenna to source}
      Ras(2)=(VoxLoc(2)-Ant_pos(2,i_ant,i_chunk))/1000.
      Ras(3)=(VoxLoc(3)-Ant_pos(3,i_ant,i_chunk))/1000.
      HorDist= Ras(1)*Ras(1) + Ras(2)*Ras(2)  ! =HYPOT(X,Y)
      D=sqrt(HorDist + Ras(3)*Ras(3))
      !AntSourceD =D
      HorDist=sqrt( HorDist ) ! =HYPOT(X,Y)
      Thet_r=atan(HorDist,Ras(3))  ! Zenith angle
      Phi_r=atan2( Ras(2) , Ras(1) )  ! \phi=0 = north
      Thet_d =Thet_r*180/pi
      Phi_d =Phi_r*180/pi
      D_a(j_IntFer)=D      ! Just for plotting
      Ph_a(j_IntFer)=Phi_d
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
      w_p=W_ap(j_IntFer)
      w_t=W_at(j_IntFer)  ! do not bother about scale factor, divides out anyhow.
      sw_p=sqrt(W_ap(j_IntFer))
      sw_t=sqrt(W_at(j_IntFer))
      wap_PB(:,j_IntFer)=p_PB(:)*sw_p/D
      wat_PB(:,j_IntFer)=t_PB(:)*sw_t/D
      !
      Do i_nu=inu1,inu2   ! Increment frequency spectrum with the spectrum from this antenna
         nu=i_nu*dnu
         i_freq=Int(nu)
         dfreq=nu-i_freq
         Sp=((1.-dfreq)*Ji_p0(i_freq) + dfreq*Ji_p0(i_freq+1)) *Cnu0(i_nu,j_IntFer) + &
            ((1.-dfreq)*Ji_p1(i_freq) + dfreq*Ji_p1(i_freq+1)) *Cnu1(i_nu,j_IntFer)
         St=((1.-dfreq)*Ji_t0(i_freq) + dfreq*Ji_t0(i_freq+1)) *Cnu0(i_nu,j_IntFer) + &
            ((1.-dfreq)*Ji_t1(i_freq) + dfreq*Ji_t1(i_freq+1)) *Cnu1(i_nu,j_IntFer)
         Fnu_PB(i_nu,:)=Fnu_PB(i_nu,:) + phase *( p_PB(:)* Sp*w_p + t_PB(:)* St*w_t ) *Gain(i_freq)/D
         wEnu_p(i_nu)= phase * Sp*sw_p *Gain(i_freq)
         wEnu_t(i_nu)= phase * St*sw_t *Gain(i_freq)
         phase =phase *dphase
      Enddo
      !
      !write(2,*) 'j_IntFer:',j_IntFer, wEnu_p(inu1+1), abs(SP), phase, dphase
      ! Calculate contribution to AMat       == A matrix
      Do i=1,3
         Do j=1,3  ! a single factor distance since the other is already included in W (as in W_p & W_t)
            AMat(i,j)=AMat(i,j)+ ( p_PB(i)*p_PB(j)*w_p + t_PB(i)*t_PB(j)*w_t)/(D*D)
         Enddo
      Enddo
      !
      Call RFTransform_CF2CT(wEnu_p(0),wETime_p(1) ) !time trace measured E-field for antenna pair
      Call RFTransform_CF2CT(wEnu_t(0),wETime_t(1) )
      Do j=-N_smth,N_smth ! copy for later use
         wEtime_ap(j,j_IntFer)=wETime_P(IntfLead+i_s+j)  !time trace measured E-field for antenna pair
         wEtime_at(j,j_IntFer)=wETime_t(IntfLead+i_s+j)
      Enddo
      !
   Enddo     ! j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
   !write(2,*) 'Sp',i_ant, IntfLead+i_s, abs(Sp), abs(St), phase
   !
   Call matinv3(AMat, Ai)
   Do m=1,3  ! convert to time
      Call RFTransform_CF2CT(Fnu_PB(0,m),FTime_PB(1,m) )
   Enddo
   !
   !write(2,*) 'FTime_PB(1:10,1):',FTime_PB(1:10,1)
   !
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
   !write(2,*) 'Ai',Ai
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
      Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
         Del_p = ABS( wEtime_ap(j,j_IntFer) -SUM( wap_PB(:,j_IntFer)*AiF(:) ) ) ! true difference (measured-calculated)for this antenna
         Del_t = ABS( wEtime_at(j,j_IntFer) -SUM( wat_PB(:,j_IntFer)*AiF(:) ) )
         dChi_ap(j_IntFer)=dChi_ap(j_IntFer) + smooth(j)* Del_p**2 ! wEtime_ap(j,j_IntFer) * Conjg( Del_p  ) ! Imag is zero ??
         dChi_at(j_IntFer)=dChi_at(j_IntFer) + smooth(j)* Del_t**2 ! wEtime_at(j,j_IntFer) * Conjg( Del_t  ) ! Imag is zero ??
         ! SumDiff=SumDiff+ sqrt(smooth(j))*(Del_p + Del_t)
         SumSq=SumSq + smooth(j)* (Del_p**2 + Del_t**2) ! True contribution to chi^2 from this antenna
      Enddo
      !write(2,*) 'j',j,smooth(j), Del_p, Del_t,AiF(1), FdotI
      If(TestCh2) then
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
   !Chi2pDF=(Esq_ak-FdotI)/(2.*Nr_IntFer)    !  Reasonable approximation to chi^2
   Chi2pDF=SumSq/(2.*Nr_IntFer)
   StI = Real(Stk(1,1) + Stk(2,2) + Stk(3,3))  ! =\Lambda_0  from  PHYSICAL REVIEW E 66, 016615 ~2002!
   StI12=Real(Stk(1,1) + Stk(2,2))             ! =
   StI3 =  Real(Stk(3,3))                      ! =>\Lambda_8 = StI12 *sq(3)/2 - StI3 *sq(3)
   StQ = Real(Stk(1,1)-Stk(2,2))               ! =\Lambda_3 *2/3
   StU = 2*Real(Stk(1,2))                      ! =\Lambda_1 *2/3
   StV =  2*Imag(Stk(1,2))                     ! =\Lambda_2 *2/3
   StU2 = 2*Real(Stk(3,1))                     ! =\Lambda_4 *2/3
   StV2 =  2*Imag(Stk(3,1))                    ! =\Lambda_5 *2/3
   StU1 = 2*Real(Stk(2,3))                     ! =\Lambda_6 *2/3
   StV1 =  2*Imag(Stk(2,3))                    ! =\Lambda_7 *2/3
   St8 = (StI12  - 2* StI3)/sqrt(3.)           ! =\Lambda_8 *2/3
   P_un = 1.- (3./4.)*(StQ**2+StU**2+StV**2+StU1**2+StV1**2+StU2**2+StV2**2+St8**2)/(StI**2)
   P_lin = (3./4.)*(StQ**2+StU**2+StU1**2+StU2**2+St8**2)/(StI**2)
   P_circ= (3./4.)*(StV**2+StV1**2+StV2**2)/(StI**2)
   W=ATAN2(StU,StQ)/2.*180./pi
   If(W.lt.-90.) W=W+180.
   write(2,"(I4,2(A,G12.4),10(A,F7.2))") i_slice,': Wgt: I123=',StI,', I12=',StI12, ', Q/I=',StQ/StI, &
         ', U/I=',StU/StI, ', V/I=',StV/StI, ', I3/I=',StI3/StI,', angle=',W,', chi^2=',Chi2pDF &
         ,', P_unpol=',P_un, ', P_lin=',P_lin, ', P_circ=',P_circ
   !
   !
   !
   Do m=1,3
      Do n=1,3
         AMat(m,n)=SUM( Ai(m,:)*Ai(:,n) )*Chi2pDF**2  /4. ! this should correspond to the square of the error
         ! arbitrary factor 100 included to stretch error bars and make them visible
      EndDo
    !  write(2,*) 'Correlation matrix(m,:)=',AMat(m,:)
   EndDo
   dStI = SQRT(SUM(AMat(:,:) ))
   dStI12=SQRT(AMat(1,1) + AMat(2,2))
   dStI3 =  SQRT( AMat(3,3) )
   dStQ = dStI12
   dStU = 2*SQRT( ABS(AMat(1,2)) )
   dStV =  dStU
   dStU2 = 2*SQRT( ABS(AMat(3,1)) )
   dStV2 =  dStU2
   dStU1 = 2*SQRT( ABS(AMat(2,3)) )
   dStV1 =  dStU1
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
         i_ant=IntFer_ant(j_IntFer,i_chunk)
         write(30,"(I4,I5,I6,I6,F9.3,F8.1,5G12.3)")  j_IntFer, i_ant, &
            Ant_IDs(I_ant,i_chunk), Ant_Stations(I_ant,i_chunk), D_a(j_IntFer), Ph_a(j_IntFer), &
            dChi_ap(j_IntFer), dChi_at(j_IntFer), &
            W_ap(j_IntFer)/W, W_at(j_IntFer)/W
      Enddo
      Close(Unit=30)
      !call GLE -d pdf -o NLMx_dChi2_01.pdf %UtilDir%EI_chisq.gle "../19A-5/files/NLtst" "01"
      Call GLEplotControl(PlotType='EI_delta-chisq', PlotName='EI_delchisq_'//TRIM(txt)//TRIM(OutFileLabel), &
            PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel)//' "'//TRIM(txt)//'"' ) ! does NOT work in windows
      Call GLEplotControl(PlotType='EI_chisq', PlotName='EI_chisq_'//TRIM(txt)//TRIM(OutFileLabel), &
            PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel)//' "'//TRIM(txt)//'"' ) ! does NOT work in windows
      Call GLEplotControl(SpecialCmnd='rm '//TRIM(DataFolder)//TRIM(OutFileLabel)//'IntfSpecChi2_'//TRIM(txt)//'.dat') !  Command to delete files
      Call GLEplotControl(SpecialCmnd='rm '//TRIM(DataFolder)//TRIM(OutFileLabel)//'DelChi2_'//TRIM(txt)//'.dat') !  Command to delete files
   EndIf
   !
   Return
End Subroutine EI_PolGrid
!-----------------------------------------------
