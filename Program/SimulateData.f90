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
   !use DataConstants, only : Time_dim, Cnu_dim, Station_nrMax
   use DataConstants, only : Station_nrMax
   !use Chunk_AntInfo, only : Start_time
   !use Chunk_AntInfo, only : Ant_Stations, Ant_IDs, Ant_nr,
   !use Chunk_AntInfo, only :  Ant_RawSourceDist
   !use Chunk_AntInfo, only : Powr_eo,NAnt_eo
   use AntFunCconst, only : Freq_min, Freq_max
   use AntFunCconst, only : J_0p, J_0t, J_1p, J_1t,Gain
   use AntFunCconst, only : Ji_p0,Ji_t0,Ji_p1,Ji_t1
   !use StationMnemonics, only : Statn_ID2Mnem, Statn_Mnem2ID
   use FFT, only : RFTransform_CF, RFTransform_CF2CT, RFTransform_CF_Filt, RFTransform_su, DAssignFFT
!
   Implicit none
   !
   Integer, parameter :: NtSamples= 2**10 ! =1024 !Time_dim !
   Integer, parameter :: NnuSamples=NtSamples/2
   Integer, parameter :: SrcNrMax=100
   Integer, parameter :: AntpStat=20  ! Max Nr antennas per station
   Real(dp) :: RTime_s(1:NtSamples), RTime_0(1:NtSamples), RTime_1(1:NtSamples)
   Complex(dp) :: FTime_0(1:NtSamples,1:AntpStat), FTime_1(1:NtSamples,1:AntpStat)
   Complex(dp) :: CNu_s(0:NnuSamples), CNu_0(0:NnuSamples), CNu_1(0:NnuSamples), CNu_p(0:NnuSamples), CNu_t(0:NnuSamples)
   Real(dp) :: nu_Fltr(0:NnuSamples) !, Av, Bv, FiltFact
   Real(dp) :: Ant_pos(1:3,1:AntpStat)
   Character(len=12) :: Station_name, Lab1, Lab2, Lab3, Lab4
   Integer :: i_file, i_ant, j_ant, Ant_nr, i_sample, k, Ant_IDs(1:AntpStat)
   Integer,save :: EvenOdd=-1
   Integer :: inu1, inu2, i_eo, STATION_ID, Ant_ID, Sample_Offset, nxx
   Real(dp) :: StatStartTime, SubSample_Offset, LFRAnt_crdnts(3), RDist, T_Offset, Powr !,StatAnt_Calib !,StartTime_ms
   Real(dp) :: SourcesListLocNEh(3,SrcNrMax), SourcesListTms(SrcNrMax), SourcesListAmpNEh(3,SrcNrMax)
   Real(dp) :: SourcesListTS(SrcNrMax)  !  Souce time in samples
   Real(dp) :: Time_Interval, Sources_TS, Sources_LocNEh(1:3), Sources_AmpNEh(1:3)
   Integer :: Time_nr
   !Character(len=4),save :: Signature='----'
   Real(dp) :: Ras(1:3), Vec_p(1:3), Vec_t(1:3)
   Real(dp) :: D, HorDist, A_p,A_t, FracGalacNoisePow, GN, IstN, TotalGain, NormPulseAmp0Bckg
   Real(dp) :: PulseRespWidth=20. !samples
   Real(dp) :: dFreq, nu, dnu, Phi_d, Phi_r, Thet_d, Thet_r, SourceGuess(3), Start_time,t_shft
   Integer :: i_freq, i_nu, i_src, SrcNr, i_sampl
   Character(LEN=180) :: lname
   Character(len=20) :: Utility, release, Antennas, Simulation, OutFileLabel
   INTEGER :: DATE_T(8),i
   Real :: random_stdnormal
   NAMELIST /Parameters/ Antennas, Simulation, FracGalacNoisePow, OutFileLabel !, &
   !   SourcesListLocNEh, SourcesListTms, SourcesListAmpNEh
   !
   Open(unit=2,STATUS='unknown',ACTION='write', FILE ="Simulate.out")
   FracGalacNoisePow=0.5
   OutFileLabel=""
   Read(*,NML = Parameters)
   If(FracGalacNoisePow.gt.1. .or. FracGalacNoisePow.lt. 0.) Then
      FracGalacNoisePow=0.5
      write(2,*) '***** FracGalacNoisePow changed to',FracGalacNoisePow
   Endif
   write(2,NML = Parameters)
   Call AntFieParGen()
   TotalGain=0.
   Do i_freq=Freq_min,Freq_max
      TotalGain=TotalGain + Gain(i_freq)*Gain(i_freq)
   EndDo
   TotalGain=sqrt(TotalGain)*(Freq_max-Freq_min) ! take out bandwidth to get total power gain
   If( TRIM(Simulation).eq.TRIM(Antennas) ) then
      write(2,*) 'Simulation and Antennas should be different: ', TRIM(Simulation),'=', TRIM(Antennas)
      Stop
   endif
   Open(unit=14,STATUS='old',ACTION='read', FILE = 'files/'//TRIM(Antennas)//'_Structure.dat')
   call random_seed()
   Open(unit=24,STATUS='unknown',ACTION='write', FILE = 'files/'//TRIM(Simulation)//'_Structure.dat')
   !
   SrcNr=0
   Do !i_src=1,SrcNrMax
      If(SrcNr.eq.SrcNrMax) exit
      Call GetNonZeroLine(lname)
      Read(lname,*,IOSTAT=nxx) Lab1
      If(nxx.ne.0) exit
      If((Lab1(1:3) .eq. 'Rep') .or. (Lab1(1:3) .eq. 'rep')) Then
         Read(lname,*,IOSTAT=nxx) Lab1, Time_Interval, Time_nr, Sources_TS, Sources_LocNEh(:), Sources_AmpNEh(:)
         If(nxx.ne.0) then
            Write(2,*) 'Expected: "Lab1, Time_Interval, Time_nr, Sources_TS, Sources_LocNEh(1:3), Sources_AmpNEh(1:3)" but got '&
               ,TRIM(lname)
            stop
         Endif
         If(Time_nr.gt.(SrcNrMax-SrcNr)) Time_nr=(SrcNrMax-SrcNr)
         Do i_src=0,Time_nr-1
            SrcNr=SrcNr+1
            SourcesListTS(SrcNr)=Sources_TS + i_src*Time_Interval
            SourcesListLocNEh(:,SrcNr)=Sources_LocNEh(:)
            SourcesListAmpNEh(:,SrcNr)=Sources_AmpNEh(:)
         EndDo
         cycle
      EndIf
      Read(lname,*,IOSTAT=nxx) Sources_TS, Sources_LocNEh(:), Sources_AmpNEh(:)
      If(nxx.ne.0) Then
         exit  ! terminate reading sequence
      EndIf
      SrcNr=SrcNr+1
      SourcesListTS(SrcNr)=Sources_TS
      SourcesListLocNEh(:,SrcNr)=Sources_LocNEh(:)
      SourcesListAmpNEh(:,SrcNr)=Sources_AmpNEh(:)
      ! constants tunes as to have background (sum square freq. spectrum) =1.
      ! source at distance D [km] with intensity=D gives power per 20 samples = 1 = backgrnd
   EndDo
   If(SrcNr.le.0) Then
      write(2,*) 'No valid sources have been entered:',lname
   EndIf
   Do i_src=1,SrcNr
      t_shft=sqrt(SUM(SourcesListLocNEh(:,i_src)*SourcesListLocNEh(:,i_src)))*1000.*Refrac/c_mps ! in [ms] due to signal travel distance
      SourcesListTms(i_src)=SourcesListTS(i_src)*(Sample*1000.)+t_shft  ! source time at core
      write(2,*) i_src,(SourcesListTms(i_src)-SourcesListTms(1))/(Sample*1000.), &
         SourcesListTS(i_src), SourcesListLocNEh(:,i_src), SourcesListAmpNEh(:,i_src)
   Enddo
   write(2,*) 'Nr of sources=',SrcNr, Station_NrMax, ', time first source in ref antenna=',SourcesListTms(1)
   Flush(unit=2)
   !
   Call RFTransform_su(NtSamples)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   dnu=100./NnuSamples   ! [MHz] Jones matrix is stored on 1MHz grid
   inu1=Int(Freq_min/dnu)
   inu2=Int(Freq_max/dnu-0.5)+1  ! set integration regime to frequencies within filter band
   nu_Fltr(:)=0.
   nu_Fltr(inu1:inu2)=1.
   !
   Station_name='CXnnn'
   SourceGuess(:)=SourcesListLocNEh(:,1)
   Do i_file=1,Station_NrMax ! 10 !80
      read(14,*,IOSTAT=nxx) STATION_ID, Station_name
      If(nxx.ne.0) exit
      write(24,*,IOSTAT=nxx) STATION_ID, Station_name
      !Station_number(i_file)=STATION_ID
      !Flush(unit=2)
      Open(unit=12,STATUS='old',ACTION='read', FILE = 'files/'//TRIM(Antennas)//'_'//TRIM(Station_name)//'.dat')
      If(EvenOdd.eq.-1) Then  ! Check is station numbers were specified for the first time reading a file like this
         Read(12,*) Lab1
         Read(12,*) Lab1
         read(Lab1,*,IOSTAT=nxx) Ant_ID  ! is an integer
         If(nxx.eq.0) Then ! there was a number at this place, all antennas are pairs infact
            EvenOdd=1
            !Signature='Dual' ! It is assumed that antanna-pairs are meant, not single antennas. Number is ignored.
         Else
            EvenOdd=0
            Write(2,*) 'antennas should be paired', i_file,Lab1
            stop 'antennas should be paired'
         EndIf
         rewind(unit=12)
      EndIf
      Open(unit=22,STATUS='unknown',ACTION='write', FILE = 'files/'//TRIM(Simulation)//'_'//TRIM(Station_name)//'.dat')
      !write(2,*) 'EvenOdd',EvenOdd,Signature
      !flush(Unit=2)
      Read(12,*) Lab1! , StatStartTime, Lab2, i_src, lab3, NrAnt, lab4, powr  ![ms],,&  =noise power
      !write(2,*) Lab1,STATION_ID, Station_name
      !write(2,*) Lab1, StatStartTime, Lab2, NrSamples, lab3, NrAnt, lab4, powr
      !if(powr.lt.1.d-36) Powr=1.d-36  ! =noise power
      j_ant=0
      Do
         !If(NrAnt.eq.j_ant-1) exit
         Read (12,*,IOSTAT=nxx) Ant_ID, lab1, LFRAnt_crdnts(1:3)
         !write(2,*) 'input',Ant_ID, lab1, LFRAnt_crdnts(1:3)
         ! value of Ant_ID is ignored. all antennas are taken as pairs.
         !Read (12,*,IOSTAT=nxx) Signature, Ant_ID,  LFRAnt_crdnts
         If(nxx.ne.0) exit ! end of antenna list signaled.
         j_ant=j_ant+1
         !Call RelDist(SourceGuess,LFRAnt_crdnts,RDist)
         Ant_IDs(j_ant)=10*j_ant
         Ant_pos(:,j_ant)=LFRAnt_crdnts(:)
         !Write(2,*) j_ant,Ant_pos(:,j_ant)
         If(j_ant.ge.AntpStat) then
            write(2,*) 'nr of antennas too large in SimulationRead:',i_file,STATION_ID,j_ant
            exit
         EndIf
         !
         !flush(Unit=2)
      EndDo
      Close(unit=12)
      Ant_nr=j_ant
      !
      ! Start constructing traces
      !
      Call RelDist(SourceGuess,Ant_pos(:,1),RDist)
      StatStartTime=SourcesListTms(1)/(Sample*1000.)+RDist-200 ! place first source at sample 50 in trace
      write(22,*) 'StartTime[ms]=', StatStartTime*(Sample*1000.), 'N_samples=', NtSamples, &
            'N_ant=', Ant_nr, 'P_noise= 1.0'  ![ms],,&  =noise power
      ! reconstruction should be done best with Start_time(i_chunk) (=time in reference antenna)=StatStartTime for reference antenna
      write(2,"(A,I3,I9,1x,A5,A,I2,A,A,F10.5,A,F9.6)") 'Station=',i_file,STATION_ID, Station_name, &
         ' uses',j_ant,' antenna pairs.', ' Start time=',StatStartTime*(Sample*1000.),'[ms]'
      Flush(unit=2)
      Do j_ant=1,Ant_nr  ! write spectra
         write(22,*) Ant_IDs(j_ant), 'NEh=', Ant_pos(:,j_ant)
         !i_ant=AntNr_lw+j_ant
         ! Instrumental noise:
         Do i_sampl=1,NtSamples
            RTime_0(i_sampl)=random_stdnormal()
            ! sqrt(2.) appears to correct for the difference in effective bandwidth (or peak-spreading width) due to the antenna function
            RTime_1(i_sampl)=random_stdnormal()
         EndDo
         !write(2,*) 'Instrumental t-background=',SUM(RTime_0(:)),SUM(RTime_0(:)**2)/NtSamples, &
         !   SUM(RTime_1(:)) ,SUM(RTime_1(:)**2)/NtSamples !,TotalGain
         Call RFTransform_CF(RTime_0,Cnu_0(0))
         Call RFTransform_CF(RTime_1,Cnu_1(0)) ! pure noise, NOT modelled by antenna function
         !write(2,*) 'Instrumental background=',SUM(ABS(Cnu_0(:))**2),SUM(ABS(Cnu_1(:))**2) !,TotalGain
         !
         ! Model galactic noise, coming from zenith and modified by antenna function
         Do i_sampl=1,NtSamples
            RTime_0(i_sampl)=random_stdnormal()
            RTime_1(i_sampl)=random_stdnormal()
         EndDo
         Call RFTransform_CF(RTime_0,Cnu_p(0))
         Call RFTransform_CF(RTime_1,Cnu_t(0)) ! pure noise, flat frequency spectrum
         ! Get Zenithal antenna function
         !Cnu_0(:)=0.
         !Cnu_1(:)=0.
         Thet_d =0.
         Phi_d =0.
         Call AntFun(thet_d ,Phi_d ) ! sets J_0p,J_0t,J_1p,J_1t;  Jones; gives _p & _t polarized fields when contracted with (0.1) voltages
         J_0p(:)=CONJG(J_0p(:))
         J_1p(:)=CONJG(J_1p(:))
         J_0t(:)=CONJG(J_0t(:))
         J_1t(:)=CONJG(J_1t(:))
         IstN=sqrt((1.-FracGalacNoisePow)*0.33)
         GN=sqrt(FracGalacNoisePow*12.) ! phenomenological factor to have similar power in GalNoise as in InstNoise (when GalacNoise=1.)
         GN= GN*(inu2+inu1)*dnu/TotalGain
         Do i_nu=inu1,inu2   ! Pull through antenna function with 1/nu frequency spectrum
            nu=i_nu*dnu
            i_freq=Int(nu)
            !write(2,*) i_freq, Ji_p0(i_freq)*J_0p(i_freq)+Ji_p1(i_freq)*J_1p(i_freq), &
            !      Ji_t0(i_freq)*J_0p(i_freq)+Ji_t1(i_freq)*J_1p(i_freq)
            dfreq=nu-i_freq
            Cnu_0(i_nu)=IstN*Cnu_0(i_nu) + &
                        ((1.-dfreq)*J_0p(i_freq) + dfreq*J_0p(i_freq+1)) * Cnu_p(i_nu)*GN/nu + &
                        ((1.-dfreq)*J_0t(i_freq) + dfreq*J_0t(i_freq+1)) * Cnu_t(i_nu)*GN/nu
            Cnu_1(i_nu)=IstN*Cnu_1(i_nu) + &
                        ((1.-dfreq)*J_1p(i_freq) + dfreq*J_1p(i_freq+1)) * Cnu_p(i_nu)*GN/nu + &
                        ((1.-dfreq)*J_1t(i_freq) + dfreq*J_1t(i_freq+1)) * Cnu_t(i_nu)*GN/nu
         Enddo
         !write(2,*) 'Galactic background added=',SUM(ABS(Cnu_0(:))**2),SUM(ABS(Cnu_1(:))**2)
         !
         !Cnu_0(:)=0.
         !Cnu_1(:)=0.
         NormPulseAmp0Bckg=PulseRespWidth*sqrt(14.)/TotalGain
         Do i_src=1,SrcNr
            !write(2,*) 'Galactic background=',SUM(ABS(Cnu_p(:))**2),SUM(ABS(Cnu_t(:))**2)
            Call RelDist(SourcesListLocNEh(:,i_src),Ant_pos(:,j_ant),RDist)
            T_Offset=RDist - StatStartTime +SourcesListTms(i_src)/(Sample*1000.)  ! in units of samples
            Sample_Offset = INT(T_Offset) ! in units of sample size
            If(Sample_Offset.le.0 .or. Sample_Offset.ge.NtSamples) Then
               write(2,*) 'Sample_Offset out of range for src=', i_src, Sample_Offset
               Stop 'Sample_Offset out of range'
            EndIf
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
            Call AntFun(thet_d ,Phi_d ) ! sets J_0p,J_0t,J_1p,J_1t;  Jones; gives _p & _t polarized fields when contracted with (0.1) voltages
            Vec_p(1)=sin(Phi_r)              ; Vec_p(2)=-cos(Phi_r)        ; Vec_p(3)=0.
            Vec_t(1)=-cos(Thet_r)*Vec_p(2)  ; Vec_t(2)=cos(Thet_r)*Vec_p(1) ; Vec_t(3)=-sin(Thet_r)
            J_0p(:)=CONJG(J_0p(:))
            J_1p(:)=CONJG(J_1p(:))
            J_0t(:)=CONJG(J_0t(:))
            J_1t(:)=CONJG(J_1t(:))
            !Call AntFun_Inv(thet_d, Phi_d ) ! sets ,Ji_p0,Ji_t0,Ji_p1,Ji_t1; Inverse Jones; gives _p & _t polarized fields when contracted with (0.1) voltages
            ! Amplitude in theta & phi directions
            A_p=SUM( SourcesListAmpNEh(:,i_src)*Vec_p(:) )*NormPulseAmp0Bckg/D
            A_t=SUM( SourcesListAmpNEh(:,i_src)*Vec_t(:) )*NormPulseAmp0Bckg/D
            !If(i_file.eq.1) then
              ! Write(2,*) Ras(:), SUM(Ras(:)*SourcesListLocNEh(:,i_src)), D*D
               !Write(2,*) SUM(SourcesListAmpNEh(:,i_src)*SourcesListLocNEh(:,i_src)),Ras(:),SourcesListLocNEh(:,i_src)
               !Write(2,*) j_ant,Ant_pos(:,j_ant)
               !write(2,*) 'PolDirection',i_src,Vec_p(:), A_p, Vec_t(:), A_t, SUM( SourcesListAmpNEh(:,i_src)*Ras(:) )
               !write(2,*) 'PolDirection',i_src, A_p, A_t, SUM( SourcesListAmpNEh(:,i_src)*Ras(:) )
            !EndIf
            Do i_nu=inu1,inu2   ! Increment frequency spectrum of the antenna with the signal from this source
               nu=i_nu*dnu
               i_freq=Int(nu)
               !write(2,*) i_freq, Ji_p0(i_freq)*J_0p(i_freq)+Ji_p1(i_freq)*J_1p(i_freq), &
               !      Ji_t0(i_freq)*J_0p(i_freq)+Ji_t1(i_freq)*J_1p(i_freq)
               dfreq=nu-i_freq
               Cnu_0(i_nu)=Cnu_0(i_nu) + &
                           ((1.-dfreq)*J_0p(i_freq) + dfreq*J_0p(i_freq+1)) * (A_p*Cnu_s(i_nu)) + &
                           ((1.-dfreq)*J_0t(i_freq) + dfreq*J_0t(i_freq+1)) * (A_t*Cnu_s(i_nu))
               Cnu_1(i_nu)=Cnu_1(i_nu) + &
                           ((1.-dfreq)*J_1p(i_freq) + dfreq*J_1p(i_freq+1)) * (A_p*Cnu_s(i_nu)) + &
                           ((1.-dfreq)*J_1t(i_freq) + dfreq*J_1t(i_freq+1)) * (A_t*Cnu_s(i_nu))
            Enddo
            !write(2,*) 'Src=',i_src,SUM(ABS(Cnu_0(:))**2)*NtSamples/PulseRespWidth, &
            !      SUM(ABS(Cnu_1(:))**2)*NtSamples/PulseRespWidth, A_p, A_t, Sample_Offset
            !stop
         EndDo
         Call RFTransform_CF2CT(Cnu_0(0),FTime_0(1,j_ant) )
         Call RFTransform_CF2CT(Cnu_1(0),FTime_1(1,j_ant) )
         !
      Enddo ! j_ant=1,NrAnt
      !
      Do i_sampl=1,NtSamples
         write(22,*) (REAL(FTime_0(i_sampl,j_ant)), REAL(FTime_1(i_sampl,j_ant)), j_ant=1,Ant_nr)
         !write(2,*) i_sample, Sample_Offset
      Enddo
      !
   Enddo !  i_file=1,StationNrMax
   Close(unit=14)
   Call DAssignFFT()
   !
   Stop
End Program Simulate_Data
