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
!-----------------------------------
Program Simulate_Data
   use constants, only : dp,sample,pi, Refrac, c_mps ! ,Refrac,pi,ci
   use DataConstants, only : Station_nrMax
   use AntFunCconst, only : Freq_min, Freq_max
   use AntFunCconst, only : J_0p, J_0t, J_1p, J_1t,Gain
   use AntFunCconst, only : Ji_p0,Ji_t0,Ji_p1,Ji_t1
   use FFT, only : RFTransform_CF, RFTransform_CF2CT, RFTransform_CF_Filt, RFTransform_su, DAssignFFT, RFTransform_CF2RT
! from: https://gcc.gnu.org/onlinedocs/gcc-4.9.4/gfortran/SELECTED_005fCHAR_005fKIND.html
         !Description:
         !    SELECTED_CHAR_KIND(NAME) returns the kind value for the character set named NAME, if a character set with such a name is supported, or -1 otherwise. Currently, supported character sets include “ASCII” and “DEFAULT”, which are equivalent, and “ISO_10646” (Universal Character Set, UCS-4) which is commonly known as Unicode.
         !
         !            use iso_fortran_env
         !            !    use, intrinsic :: iso_fortran_env, only: output_unit
         !            implicit none
         !            integer, parameter :: ascii = selected_char_kind ("ascii")
         !            integer, parameter :: ucs4  = selected_char_kind ('ISO_10646')
         !            character(kind=ascii, len=26) :: alphabet
         !            character(kind=ucs4,  len=30) :: hello_world
         !            alphabet = ascii_"abcdefghijklmnopqrstuvwxyz"
         !            hello_world = ucs4_'Hello World and Ni Hao -- ' &
         !                          // char (int (z'4F60'), ucs4)     &
         !                          // char (int (z'597D'), ucs4)
         !            write (*,*) alphabet
         !            open (output_unit, encoding='UTF-8')  !  appears not to do what I want using winedt for the output
         !            write (*,*) trim (hello_world)
         !    use, intrinsic :: iso_fortran_env, only: output_unit
         !    integer, parameter :: ucs2 = selected_char_kind('ISO_10646')
         !    character(kind=ucs2, len=:), allocatable :: str
         !    str = ucs2_'Unicode character: \u263B'
         !  open (output_unit, encoding='UTF-8')
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
   Complex(dp), allocatable :: CTime_0(:), CTime_1(:)
   Real(dp), allocatable :: FTime_0(:,:), FTime_1(:,:)
   Complex(dp), allocatable :: CNu_s(:), CNu_0(:), CNu_1(:), CNu_p(:), CNu_t(:)
   Real(dp), allocatable :: nu_Fltr(:) !, Av, Bv, FiltFact
   !
   Real(dp) :: Ras(1:3), Vec_p(1:3), Vec_t(1:3)
   Real(dp) :: D, HorDist, A_p,A_t, FracGalacNoisePow, GN, IstN, TotalGain, NormPulseAmp0Bckg, TimingErr_ns
   Real(dp) :: PulseRespWidth=20. !samples
   Real(dp) :: dFreq, nu, dnu, Phi_d, Phi_r, Thet_d, Thet_r, SourceGuess(3), StartT_sam,t_shft
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
   NAMELIST /Parameters/ Antennas, Simulation, FracGalacNoisePow, OutFileLabel, SrcNrMax, NtSamples, TimingErr_ns !, &
   !   SourcesListLocNEh, SourcesListTms, SourcesListAmpNEh
   !
   Open(unit=2,STATUS='unknown',ACTION='write', FILE ="Simulate.out")
   FracGalacNoisePow=0.5
   OutFileLabel=""
   TimingErr_ns = 1 ! [ns]
   !
   Read(*,NML = Parameters)
   !
   If(FracGalacNoisePow.gt.1. .or. FracGalacNoisePow.lt. 0.) Then
      FracGalacNoisePow=0.5
      write(2,*) '***** FracGalacNoisePow changed to',FracGalacNoisePow
   Endif
   !
   If(SrcNrMax.lt.1) Then
      SrcNrMax=100
   EndIf
   Write(2,*) 'Max number of sources = SrcNrMax =', SrcNrMax
   allocate(SourcesListLocNEh(1:3,SrcNrMax), SourcesListTms(SrcNrMax), SourcesListAmpNEh(3,SrcNrMax))
   allocate(SourcesListTS(SrcNrMax) ) !  Souce time in samples
   !
   If(NtSamples.lt.400) Then
      NtSamples=500
   Else
   NtSamples=NtSamples-1
   EndIf
   Do i_sampl=1,14  ! Time_dim in the imaging progrm is set to 2^16, this better be less
      NtSamples=shiftr( NtSamples, 1 )  ! basically dividing by 2
      If(NtSamples.le.0) exit
   Enddo
   NtSamples=shiftl( 1, i_sampl )  ! Make sure it equals some multiple of 2
   NnuSamples=NtSamples/2
   Write(2,*) 'Number of samples in the time trace = NtSamples =',NtSamples, ' =2^',i_sampl
   !write(2,"(A,o10,o10,i7)") 'i:',i_sampl,NtSamples
   allocate( RTime_s(1:NtSamples), RTime_0(1:NtSamples), RTime_1(1:NtSamples))
   allocate( CTime_0(1:NtSamples), CTime_1(1:NtSamples) )
   allocate( FTime_0(1:NtSamples,1:AntpStat), FTime_1(1:NtSamples,1:AntpStat))
   allocate( CNu_s(0:NnuSamples), CNu_0(0:NnuSamples), CNu_1(0:NnuSamples), CNu_p(0:NnuSamples), CNu_t(0:NnuSamples))
   allocate( nu_Fltr(0:NnuSamples) )!, Av, Bv, FiltFact
   !
   If( TRIM(Simulation).eq.TRIM(Antennas) ) then
      write(2,*) 'Simulation and Antennas should be different: ', TRIM(Simulation),'=', TRIM(Antennas)
      Stop
   endif
   !
   write(2,NML = Parameters)
   Call AntFieParGen()
   !
   TotalGain=0.
   Do i_freq=Freq_min,Freq_max
      TotalGain=TotalGain + Gain(i_freq)*Gain(i_freq)
   EndDo
   TotalGain=sqrt(TotalGain)*(Freq_max-Freq_min) ! take out bandwidth to get total power gain
   !
   Open(unit=14,STATUS='old',ACTION='read', FILE = 'files/'//TRIM(Antennas)//'_Structure.dat')
   !
   Call CreateNewFolder(Simulation) ! create new folder when needed
   Open(unit=24,STATUS='unknown',ACTION='write', FILE = 'files/'//TRIM(Simulation)//'_Structure.dat')
   !
   call random_seed()
   Call GetSources(SrcNrMax,SourcesListLocNEh, SourcesListTms, SourcesListAmpNEh, SourcesListTS,SrcNr)
   !
   Call RFTransform_su(NtSamples)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   dnu=100./NnuSamples   ! [MHz] Jones matrix is stored on 1MHz grid
   inu1=Int(Freq_min/dnu)+1
   inu2=Int(Freq_max/dnu)-1
   !inu2=Int(Freq_max/dnu-0.5)+1  ! set integration regime to frequencies within filter band
   nu_Fltr(:)=0.
   nu_Fltr(inu1:inu2)=1.
   !
   Station_name='CXnnn'
   SourceGuess(:)=SourcesListLocNEh(:,1)
   Do i_file=1,Station_NrMax ! 10 !80
      read(14,*,IOSTAT=nxx) STATION_ID, Station_name
      If(nxx.ne.0) exit
      write(24,*,IOSTAT=nxx) STATION_ID, Station_name
      Open(unit=12,STATUS='old',ACTION='read', FILE = 'files/'//TRIM(Antennas)//'_'//TRIM(Station_name)//'.dat')
      !If(EvenOdd.eq.-1) Then  ! Check is station numbers were specified for the first time reading a file like this
      !   Read(12,*) Lab1
      !   Read(12,*) Lab1
      !   read(Lab1,*,IOSTAT=nxx) Ant_ID  ! is an integer
      !   If(nxx.eq.0) Then ! there was a number at this place, all antennas are pairs infact
      !      EvenOdd=1
      !      !Signature='Dual' ! It is assumed that antanna-pairs are meant, not single antennas. Number is ignored.
      !   Else
      !      EvenOdd=0
      !      Write(2,*) 'antennas should be paired', i_file,Lab1
      !      stop 'antennas should be paired'
      !   EndIf
      !   rewind(unit=12)
      !EndIf
      Open(unit=22,STATUS='unknown',ACTION='write', FILE = 'files/'//TRIM(Simulation)//'_'//TRIM(Station_name)//'.dat')
      Read(12,*) Lab1! , StatStartTime, Lab2, i_src, lab3, NrAnt, lab4, powr  ![ms],,&  =noise power
      j_ant=0
      Do    ! count number of antenna pairs
         Read(12,"(A180)",IOSTAT=nxx) lname
         If(nxx.ne.0) exit ! end of antenna list signaled.
         Read (lname,*,IOSTAT=nxx) lab1, lab2, LFRAnt_crdnts(1:3)
         If(nxx.ne.0) exit ! end of antenna list signaled.
         Read(lab2,*,IOSTAT=nxx) D ! when possible then started reading timespectrum, thus quit
         If(nxx.eq.0) exit ! end of antenna list signaled.
         If(trim(lab1).eq.'Even' .or. trim(lab1).eq.'Odd') cycle
         write(2,*) 'Valid antenna pair:',trim(lname)
         !write(2,*) 'input',Ant_ID, lab1, LFRAnt_crdnts(1:3)
         ! value of Ant_ID is ignored. all antennas are taken as pairs.
         j_ant=j_ant+1
         Ant_IDs(j_ant)=10*j_ant
         Ant_pos(:,j_ant)=LFRAnt_crdnts(:)
         If(j_ant.ge.AntpStat) then
            write(2,*) 'nr of antennas too large in SimulationRead:',i_file,STATION_ID,j_ant
            exit
         EndIf
         !flush(Unit=2)
      EndDo
      Close(unit=12)
      Ant_nr=j_ant
      !
      ! Start constructing traces
      !
      Call RelDist(SourceGuess,Ant_pos(:,1),RDist)
      StatStartTime=SourcesListTms(1)/(Sample*1000.)+RDist-200 ! place first source at sample 50 in trace
      StatStartTime=StatStartTime + TimingErr_ns * random_stdnormal()/5. ! to convert timing error to samples
      write(22,*) 'StartTime[ms]=', StatStartTime*(Sample*1000.), 'N_samples=', NtSamples, &
            'N_ant=', Ant_nr, 'P_noise= 1.'  ![ms],,&  =noise power
      write(2,"(A,I3,I9,1x,A5,A,I2,A,A,F10.5,A,F9.6)") 'Station=',i_file,STATION_ID, Station_name, &
         ' uses',j_ant,' antenna pairs.', ' Start time=',StatStartTime*(Sample*1000.),'[ms]'
      Flush(unit=2)
      Do j_ant=1,Ant_nr  ! Calculate and write spectra
         write(22,*) Ant_IDs(j_ant), 'NEh=', Ant_pos(:,j_ant)
         ! Instrumental noise:
         Do i_sampl=1,NtSamples
            RTime_0(i_sampl)=random_stdnormal()
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
         !J_0p(:)=CONJG(J_0p(:))
         !J_1p(:)=CONJG(J_1p(:))
         !J_0t(:)=CONJG(J_0t(:))
         !J_1t(:)=CONJG(J_1t(:))
         IstN=sqrt((1.-FracGalacNoisePow)*0.33)
         GN=sqrt(FracGalacNoisePow*12.) ! phenomenological factor to have similar power in GalNoise as in InstNoise (when GalacNoise=1.)
         GN= GN*(inu2+inu1)*dnu/TotalGain
         Do i_nu=inu1,inu2   ! Pull through antenna function with 1/nu frequency spectrum
            nu=i_nu*dnu
            i_freq=Int(nu)
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
         !stop
         !Cnu_0(:)=0.
         !Cnu_1(:)=0.
         NormPulseAmp0Bckg=PulseRespWidth*sqrt(14.)/TotalGain
         Do i_src=1,SrcNr
            !write(2,*) 'Galactic background=',SUM(ABS(Cnu_p(:))**2),SUM(ABS(Cnu_t(:))**2)
            Call RelDist(SourcesListLocNEh(:,i_src),Ant_pos(:,j_ant),RDist)
            T_Offset=RDist - StatStartTime +SourcesListTms(i_src)/(Sample*1000.)  ! in units of samples
            Sample_Offset = INT(T_Offset) ! in units of sample size
            If(Sample_Offset.le.0 .or. Sample_Offset.ge.NtSamples) Then
               !write(2,*) 'Sample_Offset out of range for src=', i_src, Sample_Offset
               cycle
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
            Call AntFun_Inv(thet_d,phi_d) ! sets J_0p,J_0t,J_1p,J_1t as well as their inverse;  Jones; gives _p & _t polarized fields when contracted with (0.1) voltages
            !Call AntFun(thet_d ,Phi_d ) ! sets J_0p,J_0t,J_1p,J_1t;  Jones; gives _p & _t polarized fields when contracted with (0.1) voltages
            Vec_p(1)=sin(Phi_r)              ; Vec_p(2)=-cos(Phi_r)        ; Vec_p(3)=0.
            Vec_t(1)=-cos(Thet_r)*Vec_p(2)  ; Vec_t(2)=cos(Thet_r)*Vec_p(1) ; Vec_t(3)=-sin(Thet_r)
            !J_0p(:)=CONJG(J_0p(:))
            !J_1p(:)=CONJG(J_1p(:))
            !J_0t(:)=CONJG(J_0t(:))
            !J_1t(:)=CONJG(J_1t(:))
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
            If((j_ant.eq.3) .and.  (i_src.eq. 1)) Then
               write(2,*) 'thet_d,phi_d', thet_d,phi_d, Ant_pos(:,j_ant), A_p/A_t
               i_freq=(Freq_min + Freq_max)/2
               write(2,*) i_freq, 'p', Ji_p0(i_freq)*J_0p(i_freq)+Ji_p1(i_freq)*J_1p(i_freq), &
                  Ji_t0(i_freq)*J_0p(i_freq)+Ji_t1(i_freq)*J_1p(i_freq),  &
               't', Ji_t0(i_freq)*J_0t(i_freq)+Ji_t1(i_freq)*J_1t(i_freq), Ji_p0(i_freq)*J_0t(i_freq)+Ji_p1(i_freq)*J_1t(i_freq)
            EndIf
            Do i_nu=inu1,inu2   ! Increment frequency spectrum of the antenna with the signal from this source
               nu=i_nu*dnu
               i_freq=Int(nu)
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
         EndDo ! i_src=1,SrcNr
         !Call RFTransform_CF2CT(Cnu_0(0),FTime_0(1,j_ant) )
         Call RFTransform_CF2RT(Cnu_0(0),FTime_0(1,j_ant) )
         !RTime_s(:)=Real(FTime_0(:,j_ant))
         !SubSample_Offset=0.
         !Call RFTransform_CF_Filt(RTime_s,nu_fltr,-SubSample_Offset,Cnu_s(0)) !
         !Call RFTransform_CF2CT(Cnu_s(0),FTime_1(1,j_ant) )
         !write(2,*) '0', FTime_0(Sample_Offset-10:Sample_Offset+10,j_ant)
         !write(2,*) 'r',FTime_0(Sample_Offset-10:Sample_Offset+10,j_ant)/FTime_1(Sample_Offset-10:Sample_Offset+10,j_ant)
         !write(2,*) 'a0', Abs(FTime_0(Sample_Offset-10:Sample_Offset+10,j_ant))
         !write(2,*) 'ar',Abs(FTime_0(Sample_Offset-10:Sample_Offset+10,j_ant)/FTime_1(Sample_Offset-10:Sample_Offset+10,j_ant))
         !stop
         !
         !Call RFTransform_CF2CT(Cnu_0(0),FTime_0(1,j_ant) )  ! RFTransform_CF2RT(Cnu,RD)
         !Call RFTransform_CF2RT(Cnu_0(0),RTime_s(1) )  ! RFTransform_CF2RT(Cnu,RD)
         !write(2,*) 'r',Real(FTime_0(Sample_Offset-10:Sample_Offset+10,j_ant))/RTime_s(Sample_Offset-10:Sample_Offset+10)
         !stop
         Call RFTransform_CF2RT(Cnu_1(0),FTime_1(1,j_ant) )
         If(i_file.eq.1 .and. j_ant.eq.1) Then
         Call RFTransform_CF2CT(Cnu_0(0),CTime_0(1) )  ! RFTransform_CF2RT(Cnu,RD)
         Call RFTransform_CF2CT(Cnu_1(0),CTime_1(1) )
            Open(unit=23,STATUS='unknown',ACTION='write', FILE = 'files/'//TRIM(Simulation)//'_RefAntTrace.dat')
            write(2,*) 'spectrum of ant=',j_ant,' written to:','files/'//TRIM(Simulation)//'_RefAntTrace.dat'
            Do i_sampl=1,NtSamples
               write(23,*) i_sampl,REAL(CTime_0(i_sampl)), ABS(CTime_0(i_sampl)), REAL(CTime_1(i_sampl)), ABS(CTime_1(i_sampl))
            Enddo
         EndIf
      Enddo ! j_ant=1,NrAnt
      !
      Do i_sampl=1,NtSamples
         write(22,*) (FTime_0(i_sampl,j_ant), FTime_1(i_sampl,j_ant), j_ant=1,Ant_nr)
         !write(2,*) i_sample, Sample_Offset
      Enddo
      Close(unit=22)
      Close(unit=12)
      !
   Enddo !  i_file=1,StationNrMax
   Close(unit=14)
   Close(unit=24)
   Call DAssignFFT()
   !
   Stop
End Program Simulate_Data
!=====================================
Subroutine GetSources(SrcNrMax,SourcesListLocNEh, SourcesListTms, SourcesListAmpNEh, SourcesListTS,SrcNr)
   use constants, only : dp,sample,pi, Refrac, c_mps ! ,Refrac,pi,ci
   !Use RandomFunc, only : random_stdnormal, random_stdnormal3D
   Implicit none
   !
   Integer, intent(in) :: SrcNrMax
   Real(dp), intent(out) :: SourcesListLocNEh(3,SrcNrMax), SourcesListTms(SrcNrMax), SourcesListAmpNEh(3,SrcNrMax)
   Real(dp), intent(out) :: SourcesListTS(SrcNrMax)  !  Souce time in samples
   Integer, intent(out) :: SrcNr
   !
   Character(len=12) :: Lab1  !, Lab2, Lab3, Lab4  Station_name,
   Real(dp) :: StatStartTime, SubSample_Offset, LFRAnt_crdnts(3), RDist, T_Offset, Powr !,StatAnt_Calib !,StartTime_ms
   Real(dp) :: Time_Interval, Sources_TS, Sources_LocNEh(1:3), Sources_AmpNEh(1:3),t_shft
   Integer :: NConstruct, nxx, i_src
   Real(dp) :: D, Phi, Thet, Time_width, Space_width, PolSpread, R,Epsln=epsilon(Epsln)
   Character(LEN=180) :: lname
   Real :: x,y,z, RGV(1:3) ! Random Gaussian distributed Vector
   Real :: random_stdnormal
   !
   SrcNr=0
   Do !i_src=1,SrcNrMax
      If(SrcNr.eq.SrcNrMax) exit
      Call GetNonZeroLine(lname)
      Read(lname,*,IOSTAT=nxx) Lab1
      If(nxx.ne.0) exit
      ! Check for repetitive source in time
      If((Lab1(1:3) .eq. 'Rep') .or. (Lab1(1:3) .eq. 'rep')) Then  ! Repetitive
         Read(lname,*,IOSTAT=nxx) Lab1, Time_Interval, NConstruct
         If(nxx.ne.0) then
            Write(2,*) 'Expected: "Lab1, Time_Interval, NConstruct" but got ',TRIM(lname)
            stop
         Endif
         Call GetNonZeroLine(lname)
         Read(lname,*,IOSTAT=nxx) Sources_TS, Sources_LocNEh(:), Sources_AmpNEh(:)
         If(nxx.ne.0) then
            Write(2,*) 'Expected: "Sources_TS[samples], Sources_LocNEh(1:3)[m], Sources_AmpNEh(1:3)" but got ' ,TRIM(lname)
            stop
         Endif
         Call Convert2m(Sources_LocNEh)
         If(NConstruct.gt.(SrcNrMax-SrcNr)) NConstruct=(SrcNrMax-SrcNr)
         Do i_src=0,NConstruct-1
            SrcNr=SrcNr+1
            SourcesListTS(SrcNr)=Sources_TS + i_src*Time_Interval
            SourcesListLocNEh(:,SrcNr)=Sources_LocNEh(:)
            SourcesListAmpNEh(:,SrcNr)=Sources_AmpNEh(:)
         EndDo
         cycle
      EndIf
      ! check for source cloud
      If((Lab1(1:3) .eq. 'Clo') .or. (Lab1(1:3) .eq. 'clo')) Then ! Cloud
         Read(lname,*,IOSTAT=nxx) Lab1, Time_width, Space_width, PolSpread, NConstruct
         If(nxx.ne.0) then
            Write(2,*) 'Expected: "Lab1, Time_width[samples], Space_width, PolSpread, NConstruct" but got ',TRIM(lname)
            stop
         Endif
         Call GetNonZeroLine(lname)
         Read(lname,*,IOSTAT=nxx) Sources_TS, Sources_LocNEh(:), Sources_AmpNEh(:)
         If(nxx.ne.0) then
            Write(2,*) 'Expected: "Sources_TS[samples], Sources_LocNEh(1:3)[m], Sources_AmpNEh(1:3)" but got ',TRIM(lname)
            stop
         Endif
         Call Convert2m(Sources_LocNEh)
         write(2,*) 'Form a cloud around:', Sources_TS, Sources_LocNEh(:)/1000., Sources_AmpNEh(:)
         write(2,*) 'Spread in time, pos, pol=',Time_width, Space_width, PolSpread
         If(NConstruct.gt.(SrcNrMax-SrcNr)) NConstruct=(SrcNrMax-SrcNr)
         Do i_src=0,NConstruct-1
            SrcNr=SrcNr+1
            !call random_number(x)  ! may include zero
            !call random_number(y)
            !call random_number(z)
            !R=sqrt((0.5-x)**2+(0.5-y)**2+(0.5-z)**2+Epsln) ! to prevent zero
            !D=random_stdnormal()! can be negative
            !D=((abs(d))**(1/3.))*Space_width  ! to have distances distributed like [d^2 x gaussian(d)]
            !!D=sqrt(abs(d)) * Space_width  ! to have distances distributed like [d x gaussian(d)]
            SourcesListTS(SrcNr)=Sources_TS + random_stdnormal()*Time_width
            !SourcesListLocNEh(1,SrcNr)=Sources_LocNEh(1)+ D*(0.5-x)/R
            !SourcesListLocNEh(2,SrcNr)=Sources_LocNEh(2)+ D*(0.5-y)/R
            SourcesListLocNEh(3,SrcNr)=Sources_LocNEh(3)+ D*(0.5-z)/R
            Call random_stdnormal3D(RGV)
            SourcesListLocNEh(:,SrcNr)=Sources_LocNEh(:)+ Space_width*RGV(:)
            Call random_stdnormal3D(RGV)
            SourcesListAmpNEh(:,SrcNr)=Sources_AmpNEh(:) + PolSpread*RGV(:)
            ! write(3,*) i_src,D,SourcesListLocNEh(:,SrcNr) ! just for checking the cloud structure using  Sources_Plots.gle
         EndDo
         cycle
      EndIf
      ! Check for single source, otherwise quit
      Read(lname,*,IOSTAT=nxx) Sources_TS, Sources_LocNEh(:), Sources_AmpNEh(:)
      If(nxx.ne.0) Then
         exit  ! terminate reading sequence
      EndIf
      Call Convert2m(Sources_LocNEh)
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
      !write(2,*) i_src,t_shft
      SourcesListTms(i_src)=SourcesListTS(i_src)*(Sample*1000.)+t_shft  ! source time at core
      write(2,*) i_src,(SourcesListTms(i_src)-SourcesListTms(1))/(Sample*1000.), &
         SourcesListTS(i_src), SourcesListLocNEh(:,i_src)/1000., SourcesListAmpNEh(:,i_src)
   Enddo
   write(2,*) 'Nr of sources=',SrcNr, ', time of first source pulse at the core=',SourcesListTms(1)
   Flush(unit=2)
End Subroutine GetSources
