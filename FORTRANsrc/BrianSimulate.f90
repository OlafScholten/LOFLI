    Include 'ConstantsModules.f90'
    Include 'FFT_routines.f90'
    Include 'ParamModules.f90'  ! v18d: Cnu storage changed  for InterfEngineB
    Include 'MappingUtilities.f90'
    !Include 'GLEplotUtil.f90'
    Include 'AntFunct.f90'
    !Include 'Ant-Read.f90'
    !Include 'CrossCorr.f90'
    !Include 'CurtainPlotOption.f90'
    Include 'System_Utilities.f90'
!-----------------------------------
Function random_stdnormal() Result(x)
!  https://masuday.github.io/fortran_tutorial/random.html
! General interest: https://en.wikibooks.org/wiki/Fortran/Fortran_procedures_and_functions
   implicit none
   real :: x
   real,parameter :: pi=3.14159265
   real :: u1,u2
   ! call random_number(r) gives 0=< r < 1 i.e. including 0, excluding 1
   call random_number(u1)
   call random_number(u2)
   x = sqrt(-2*log(1-u1))*cos(2*pi*u2)
   Return
end Function random_stdnormal
!-----------------------------------
Program Simulate_Data
   use constants, only : dp,sample,pi, Refrac, c_mps ! ,Refrac,pi,ci
   use DataConstants, only : Station_nrMax
   use AntFunCconst, only : Freq_min, Freq_max
   use AntFunCconst, only : J_0p, J_0t, J_1p, J_1t,Gain
   use AntFunCconst, only : Ji_p0,Ji_t0,Ji_p1,Ji_t1
   use FFT, only : RFTransform_CF, RFTransform_CF2CT, RFTransform_CF_Filt, RFTransform_su, DAssignFFT, RFTransform_CF2RT
!
   Implicit none
   !
   Integer, parameter :: AntpStat=20  ! Max Nr antennas per station
   Real(dp) :: Ant_pos(1:3,1:AntpStat)
   !
   Integer :: SrcNrMax=100
   Real(dp), allocatable :: SourcesListLocNEh(:,:), SourcesListTms(:), SourcesListAmpNEh(:,:)
   Real(dp), allocatable :: SourcesListTS(:)  !  Souce time in samples
   !
   Integer :: NtSamples= 2**10 ! =1024 !Time_dim !
   Integer :: NnuSamples
   Real(dp), allocatable :: RTime_s(:), RTime_0(:), RTime_1(:)
   Complex(dp), allocatable :: FTime_0(:,:), FTime_1(:,:)
   Complex(dp), allocatable :: CNu_s(:), CNu_0(:), CNu_1(:), CNu_p(:), CNu_t(:)
   Real(dp), allocatable :: nu_Fltr(:) !, Av, Bv, FiltFact
   !
   Real(dp) :: Ras(1:3), Vec_p(1:3), Vec_t(1:3)
   Real(dp) :: D, HorDist, A_p,A_t, FracGalacNoisePow, GN, IstN, TotalGain, NormPulseAmp0Bckg
   Real(dp) :: PulseRespWidth=20. !samples
   Real(dp) :: dFreq, nu, dnu, Phi_d, Phi_r, Thet_d, Thet_r, SourceGuess(3), Start_time,t_shft
   Character(len=12) :: Station_name, Lab1, Lab2, Lab3, Lab4
   Integer :: i_file, i_ant, j_ant, Ant_nr, i_sample, k, Ant_IDs(1:AntpStat)
   Integer,save :: EvenOdd=-1
   Integer :: inu1, inu2, i_eo, STATION_ID, Ant_ID, Sample_Offset, nxx
   Real(dp) :: StatStartTime, SubSample_Offset, LFRAnt_crdnts(3), RDist, T_Offset, Powr !,StatAnt_Calib !,StartTime_ms
   Integer :: i_freq, i_nu, i_src, SrcNr, i_sampl
   Character(LEN=180) :: lname
   Character(len=20) :: Utility, release, Antennas, Simulation, OutFileLabel, Folder
   INTEGER :: DATE_T(8),i
   Real :: random_stdnormal
   !
   Integer( kind = 4 ) :: trace_length !! will be length of output trace
   Integer( kind = 4 ) :: IER !! error code for FFT
   Real( kind = 8 ), DIMENSION(:), allocatable :: WSAVE(:), WORK(:)
   Integer( kind = 4 ) :: LENSAV, LENWRK
   Complex( kind = 8 ), DIMENSION(:),  allocatable :: waveform(:)   !! the waveform in frequeny or time space we are working with. Not specifying precision so is compatible with FFT
   !
   NAMELIST /Parameters/ Antennas, Simulation, FracGalacNoisePow, OutFileLabel, SrcNrMax, NtSamples !, &
   !   SourcesListLocNEh, SourcesListTms, SourcesListAmpNEh
   !
   Open(unit=2,STATUS='unknown',ACTION='write', FILE ="Simulate.out")
   !
   SrcNrMax=10
   allocate(SourcesListLocNEh(1:3,SrcNrMax), SourcesListTms(SrcNrMax), SourcesListAmpNEh(3,SrcNrMax))
   allocate(SourcesListTS(SrcNrMax) ) !  Souce time in samples
   !
   NtSamples=1024
   NtSamples=64
   NnuSamples=NtSamples/2
   Write(2,*) 'Number of samples in the time trace = NtSamples =',NtSamples, ' =2^',i_sampl
   allocate( RTime_s(1:NtSamples), RTime_0(1:NtSamples), RTime_1(1:NtSamples))
   allocate( FTime_0(1:NtSamples,1:AntpStat), FTime_1(1:NtSamples,1:AntpStat))
   allocate( CNu_s(0:NnuSamples), CNu_0(0:NnuSamples), CNu_1(0:NnuSamples), CNu_p(0:NnuSamples), CNu_t(0:NnuSamples))
   allocate( nu_Fltr(0:NnuSamples) )!, Av, Bv, FiltFact
      allocate( waveform(1:NtSamples) )
   !
   !
   write(2,*) 'calling antennapars'
   Call AntFieParGen()
   call random_seed()
   Call RFTransform_su(NtSamples)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   SrcNr=1
   SourcesListLocNEh(:,SrcNr)=(/-20000., +40120.5,   7713./)
   SourcesListAmpNEh(:,SrcNr)=(/40000.0, 20000., 0000.0/)
   !
   dnu=100./NnuSamples   ! [MHz] Jones matrix is stored on 1MHz grid
   inu1=Int(Freq_min/dnu)
   inu2=Int(Freq_max/dnu-0.5)+1  ! set integration regime to frequencies within filter band
   nu_Fltr(:)=0.
   nu_Fltr(inu1:inu2)=1.
   !
      j_ant=1
      Ant_IDs(j_ant)=10*j_ant
      Ant_pos(:,j_ant)=LFRAnt_crdnts(:)
      Ant_nr=j_ant
      !
      ! Start constructing traces
      !
      Do j_ant=1,Ant_nr  ! write spectra
         !stop
         Cnu_0(:)=0.
         Cnu_1(:)=0.
         Do i_src=1,SrcNr
            !write(2,*) 'Galactic background=',SUM(ABS(Cnu_p(:))**2),SUM(ABS(Cnu_t(:))**2)
            T_Offset=200
            Sample_Offset = INT(T_Offset) ! in units of sample size
            SubSample_Offset = T_Offset - Sample_Offset ! in units of sample size
            RTime_s(:)=0.
            RTime_s(Sample_Offset)=1.
            Call RFTransform_CF_Filt(RTime_s,nu_fltr,-SubSample_Offset,Cnu_s(0)) ! delta peak at right time
            ! get t and p oriented signals
            !
            Ras(:)=(SourcesListLocNEh(:,i_src)-Ant_pos(:,j_ant))/1000.  ! \vec{R}_{antenna to source}
            HorDist= Ras(1)*Ras(1) + Ras(2)*Ras(2)  ! =HYPOT(X,Y)
            D=sqrt(HorDist + Ras(3)*Ras(3))
            HorDist=sqrt( HorDist ) ! =HYPOT(X,Y)
            Thet_r=atan(HorDist,Ras(3))  ! Zenith angle
            Phi_r=atan2( Ras(2) , Ras(1) )  ! \phi=0 = north
            Thet_d =Thet_r*180/pi
            Phi_d =Phi_r*180/pi
            !
            write(2,*) thet_d,phi_d
            Call AntFun_Inv(thet_d,phi_d) ! sets J_0p,J_0t,J_1p,J_1t as well as their inverse;  Jones; gives _p & _t polarized fields when contracted with (0.1) voltages
            Vec_p(1)=sin(Phi_r)              ; Vec_p(2)=-cos(Phi_r)        ; Vec_p(3)=0.
            Vec_t(1)=-cos(Thet_r)*Vec_p(2)  ; Vec_t(2)=cos(Thet_r)*Vec_p(1) ; Vec_t(3)=-sin(Thet_r)
            ! Amplitude in theta & phi directions
            A_p=SUM( SourcesListAmpNEh(:,i_src)*Vec_p(:) )/D
            A_t=SUM( SourcesListAmpNEh(:,i_src)*Vec_t(:) )/D
            write(2,*) 'A_p,A_t:',A_p,A_t
            Cnu_s(:)=1.
            Do i_nu=inu1,inu2   ! Increment frequency spectrum of the antenna with the signal from this source
               nu=i_nu*dnu
               i_freq=Int(nu)
               dfreq=nu-i_freq
               Cnu_0(i_nu)=((1.-dfreq)*J_0p(i_freq) + dfreq*J_0p(i_freq+1)) * (A_p*Cnu_s(i_nu)) + &
                           ((1.-dfreq)*J_0t(i_freq) + dfreq*J_0t(i_freq+1)) * (A_t*Cnu_s(i_nu))
               Cnu_1(i_nu)=((1.-dfreq)*J_1p(i_freq) + dfreq*J_1p(i_freq+1)) * (A_p*Cnu_s(i_nu)) + &
                           ((1.-dfreq)*J_1t(i_freq) + dfreq*J_1t(i_freq+1)) * (A_t*Cnu_s(i_nu))
               !waveform(NnuSamples+i_nu)=((1.-dfreq)*J_1t(i_freq) + dfreq*J_1t(i_freq+1))
            Enddo
            !write(2,*) 'Src=',i_src,SUM(ABS(Cnu_0(:))**2)*NtSamples/PulseRespWidth, &
            !      SUM(ABS(Cnu_1(:))**2)*NtSamples/PulseRespWidth, A_p, A_t, Sample_Offset
            !stop
         EndDo ! i_src=1,SrcNr
         !CNu(:)=
         Call RFTransform_CF2CT(Cnu_0(0),FTime_0(1,j_ant) )
         write(2,*) 'H', Abs(FTime_0(:,j_ant))
         !RTime_s(:)=Real(FTime_0(:,j_ant))
         !SubSample_Offset=0.
         !Call RFTransform_CF_Filt(RTime_s,nu_fltr,-SubSample_Offset,Cnu_s(0)) !
         !Call RFTransform_CF2CT(Cnu_s(0),FTime_1(1,j_ant) )
         !write(2,*) '0', FTime_0(Sample_Offset-10:Sample_Offset+10,j_ant)
         !write(2,*) 'r',FTime_0(Sample_Offset-10:Sample_Offset+10,j_ant)/FTime_1(Sample_Offset-10:Sample_Offset+10,j_ant)
         !write(2,*) 'a0', Abs(FTime_0(Sample_Offset-10:Sample_Offset+10,j_ant))
         !write(2,*) 'ar',Abs(FTime_0(Sample_Offset-10:Sample_Offset+10,j_ant)/FTime_1(Sample_Offset-10:Sample_Offset+10,j_ant))
         !
         trace_length=NtSamples
         LENSAV = 2*trace_length + int ( log ( real ( trace_length, kind = 8 ) ) / log ( 2.0D+00 ) ) + 4
         LENWRK = 2*trace_length
         allocate( WSAVE(LENSAV) )
         allocate( WORK(LENWRK) )
         !
         IER = 0
         call CFFT1I (trace_length, WSAVE, LENSAV, IER) !! initialize FFT
         if ( IER .ne. 0  ) THEN
             stop 'cannot initialize FFT'
         END IF
         !
            waveform(:)=0.
            Do i_nu=inu1,inu2   ! Increment frequency spectrum of the antenna with the signal from this source
               waveform(NtSamples-i_nu)=Cnu_0(i_nu)
            Enddo
         call CFFT1B (trace_length, 1, waveform, trace_length,   WSAVE,  LENSAV,  WORK, LENWRK, IER) !! DO IFFT
         if ( IER .ne. 0 ) THEN
             stop 'cannot do iFFT'
         END IF
         Write(2,*) 'neg freq',abs(waveform(:))
            waveform(:)=0.
            Do i_nu=inu1,inu2   ! Increment frequency spectrum of the antenna with the signal from this source
               waveform(i_nu)=Cnu_0(i_nu)
            Enddo
         call CFFT1B (trace_length, 1, waveform, trace_length,   WSAVE,  LENSAV,  WORK, LENWRK, IER) !! DO IFFT
         if ( IER .ne. 0 ) THEN
             stop 'cannot do iFFT'
         END IF
         Write(2,*) 'pos freq',abs(waveform(:))
         stop
         !
         !Call RFTransform_CF2CT(Cnu_0(0),FTime_0(1,j_ant) )  ! RFTransform_CF2RT(Cnu,RD)
         Call RFTransform_CF2RT(Cnu_0(0),RTime_s(1) )  ! RFTransform_CF2RT(Cnu,RD)
         write(2,*) 'r',Real(FTime_0(Sample_Offset-10:Sample_Offset+10,j_ant))/RTime_s(Sample_Offset-10:Sample_Offset+10)
         stop 'simulate test'
         Call RFTransform_CF2CT(Cnu_1(0),FTime_1(1,j_ant) )
      Enddo ! j_ant=1,NrAnt
      !
   Call DAssignFFT()
   !
   Stop
End Program Simulate_Data
!=====================================
