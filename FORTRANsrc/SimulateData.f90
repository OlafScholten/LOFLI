    Include 'H_ConstantsModules.f90'
    Include 'AntFuncCnst.f90'
    Include 'FFT_routines.f90'
    Include 'H_ParamModules.f90'  ! v18d: Cnu storage changed  for InterfEngineB
    Include 'AntFunct.f90'
    Include 'H_InterferomPars.f90'
    Include 'PredefTrack.f90'
    Include 'H_MappingUtilitiesModules.f90'
    Include 'H_GLEplotUtil.f90'
    Include 'H_LOFLI_InputHandling.f90'
    Include 'H_MappingUtilities.f90'
    Include 'System_Utilities.f90'
    Include 'H_InterferometryOptSbRtns.f90'  ! d: Cnu storage changed for InterfEngineB
    Include 'H_EIOption.f90'
    Include 'H_FitParams.f90'
    Include 'HDF5_LOFAR_Read.f90'
    Include 'H_CalibrationsMod.f90'
    Include 'H_Ant-Read.f90'
    Include 'H_EICallRtns.f90'
    Include 'H_EIFitter.f90'
    Include 'H_InterferometryOption.f90'  ! d: Cnu storage changed for InterfEngineB
    Include 'H_PeakInterferoOption.f90'
    Include 'H_CurtainPlotOption.f90'
    Include 'H_CrossCorr.f90'
!    Include 'H_FindSources.f90'
    Include 'H_FindCallibr.f90'   ! Station Callibration
    Include 'MGMR3D_spline.f90'   ! Station Callibration
    Include 'H_SourceTryal.f90'  ! v16.f90' = v17.f90'  ; v17a.f90' uses grid search
    Include 'H_Fitter_CorrMax.f90'     ! uses chi-square mini
    Include 'H_FindSources.f90'
    Include 'nl2sol.f90'
!-----------------------------------
!-----------------------------------
Program Simulate_Data
   use constants, only : dp,sample,pi, ci, c_mps ! ,Refrac,pi,ci
   use DataConstants, only : Station_nrMax
   use AntFunCconst, only : Freq_min, Freq_max
   use AntFunCconst, only : J_0p, J_0t, J_1p, J_1t,Gain
   use AntFunCconst, only : Ji_p0,Ji_t0,Ji_p1,Ji_t1
   use FFT, only : RFTransform_CF, RFTransform_CF2CT, RFTransform_su, DAssignFFT, RFTransform_CF2RT  ! , RFTransform_CF_Filt
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
   !
   Integer :: SrcNrMax=100
   Real(dp), allocatable :: SourcesListLocNEh(:,:), SourcesListTms(:), SourcesListAmpNEh(:,:)
   Real(dp), allocatable :: SourcesListTS(:)  !  Souce time in samples
   !
   !Integer :: NtSamples= 2**10 ! =1024 !Time_dim !
   Integer :: NtSamples= 2**12 ! =4k ! equivalent to Time_dim !
   Integer :: NnuSamples
   Real(dp), allocatable :: RTime_s(:), RTime_0(:), RTime_1(:)
   Complex(dp), allocatable :: CTime_0(:), CTime_1(:)
   Real(dp), allocatable :: FTime_0(:,:), FTime_1(:,:)
   Complex(dp), allocatable :: CNu_s(:), CNu_0(:), CNu_1(:), CNu_p(:), CNu_t(:)
   Complex(dp) :: phShift, dphShift, CX
   Complex(dp), parameter ::   ipi=ci*pi  ! This is the complex constant i*pi
   !Real(dp), allocatable :: nu_Fltr(:) !, Av, Bv, FiltFact
   !
   Real(dp) :: Ras(1:3), Vec_p(1:3), Vec_t(1:3)
   Real(dp) :: D, HorDist, A_p,A_t, FracGalacNoisePow, GN, IstN, TotalGain, NormPulseAmp0Bckg, TimingErr_ns
   Real(dp) :: PulseRespWidth=20. !samples
   Real(dp) :: dFreq, nu, dnu, Phi_d, Phi_r, Thet_d, Thet_r, SourceGuess(3), StartT_sam,t_shft
   Character(len=12) :: Station_name, Lab1, Lab2, Lab3, Lab4
   Integer :: i_file, i_ant, j_ant, Ant_nr, i_sample, k
   Real(dp) :: Ant_pos(1:3,1:2*AntpStat), Ant_RawSourceDist(1:2*AntpStat)
   Integer :: Ant_Stations(1:2*AntpStat), Ant_IDs(1:2*AntpStat), NrAnt
   Integer,save :: EvenOdd=-1
   Integer :: inu1, inu2, i_eo, STATION_ID, Ant_ID, Sample_Offset, nxx, N_even, N_odd, N_StAnt, Unt
   Real(dp) :: StatStartTime, SubSample_Offset, LFRAnt_crdnts(3), RDist, T_Offset, Powr, N_TRID, NrmStI
   Real(dp) :: BckgrPwr_0, BckgrPwr_1, NormEvenOdd=1.d2 !,StatAnt_Calib !,StartTime_ms
   Integer :: i_freq, i_nu, i_src, SrcNr, i_sampl, IntfSmoothWin
   Character(LEN=180) :: lname
   Character(len=20) :: Antennas, Simulation, OutFileLabel ! , Folder ! Utility, release,
   INTEGER :: DATE_T(8),i
   Logical :: DualPosOnly=.true.
   Real :: random_stdnormal
   NAMELIST /Parameters/ Antennas, Simulation, FracGalacNoisePow, OutFileLabel, SrcNrMax, NtSamples, TimingErr_ns, IntfSmoothWin !, &
   !   SourcesListLocNEh, SourcesListTms, SourcesListAmpNEh
   !
   FracGalacNoisePow=0.5
   OutFileLabel=""
   TimingErr_ns = 1 ! [ns]
   IntfSmoothWin=20
   N_even=0
   N_odd=0
   BckgrPwr_0=0.
   BckgrPwr_1=0.
   !EmpNoiseNorm=sqrt(sqrt(4/pi)/2.) ! empirical fudge factor to set background=1., do not really understand, possibly related to adding random numbers for background
   ! factor 1/2 to
   !NormEvenOdd=1.d2   ! to undo the action of NormEven & NormOdd in AntennaRead
   !
   Read(*,NML = Parameters)
   Open(unit=2,STATUS='unknown',ACTION='write', FILE ="Simulate"//TRIM(OutFileLabel)//".out")
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
   !allocate( nu_Fltr(0:NnuSamples) )!, Av, Bv, FiltFact
   !
   If( TRIM(Simulation).eq.TRIM(Antennas) ) then
      write(2,*) 'Simulation and Antennas should be different: ', TRIM(Simulation),'=', TRIM(Antennas)
      Stop
   endif
   !
   write(2,NML = Parameters)
   Call AntFieParGen()
   Call RFTransform_su(NtSamples)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   dnu=100./NnuSamples   ! [MHz] Jones matrix is stored on 1MHz grid
   inu1=Int(Freq_min/dnu)+1
   inu2=Int(Freq_max/dnu)-1
   !inu2=Int(Freq_max/dnu-0.5)+1  ! set integration regime to frequencies within filter band
   !nu_Fltr(:)=0.
   !nu_Fltr(inu1:inu2)=1.
   !
   TotalGain=0.
   Do i_freq=Freq_min,Freq_max
      TotalGain=TotalGain + Gain(i_freq)*Gain(i_freq)
   EndDo
   TotalGain=sqrt(TotalGain)*(Freq_max-Freq_min) ! take out bandwidth to get total power gain
   NormPulseAmp0Bckg=PulseRespWidth*sqrt(14.)/TotalGain  !  14 probably effective frequency band LOFAR
   !
   ! Calculation of norm to agree with the TRI-D imager
   Sample_Offset = NtSamples/2 ! in units of sample size
   SubSample_Offset = 0.
   RTime_s(:)=0.
   RTime_s(Sample_Offset)=1.
   Call RFTransform_CF(RTime_s,Cnu_s(0)) ! delta peak at right time
   !Call RFTransform_CF2CT(Cnu_s(0),CTime_1(1) )
   !write(2,*) sqrt(SUM(abs(CTime_1(:))**2)),CTime_1(Sample_Offset-1),CTime_1(Sample_Offset),CTime_1(Sample_Offset+1)
   Cnu_s(0:inu1-1)=0.
   Cnu_s(inu2+1:NnuSamples)=0.
   ! previous is same as:
   !nu_Fltr(:)=1.
   !Call RFTransform_CF_Filt(RTime_s,nu_fltr,-SubSample_Offset,Cnu_s(0)) ! delta peak at right time
   !Cnu_s(0:inu1-1)=0.
   !Cnu_s(inu2+1:NnuSamples)=0.
   Do i_nu=inu1,inu2   ! Increment frequency spectrum with this antenna
      nu=i_nu*dnu
      i_freq=Int(nu)
      dfreq=nu-i_freq ! phase-shifts are zero for centran pixel, no timeshift!!
      Cnu_s(i_nu) = ((1.-dfreq)*Gain(i_freq) + dfreq*Gain(i_freq+1))* Cnu_s(i_nu)
   Enddo
   Call RFTransform_CF2CT(Cnu_s(0),CTime_1(1) )
   N_TRID=0.
   Do i_sampl=-IntfSmoothWin/2, IntfSmoothWin/2
      N_TRID=N_TRID + (ABS(CTime_1(Sample_Offset+i_sampl)))**2  ! smooth(j)*
      !write(2,*) i_sampl, ABS(CTime_1(Sample_Offset+i_sampl))
   EndDo
   N_TRID=(N_TRID/IntfSmoothWin)

   NrmStI=N_TRID*(NormPulseAmp0Bckg**2)
   i_sampl=IntfSmoothWin/2
   write(2,"(A,F7.4,A,F7.4,A,F9.3,A,I3)") 'TRI-D re-norm', sqrt(NrmStI),', Intensity at window edge=', 100* &
      (0.5*(ABS(CTime_1(Sample_Offset-i_sampl))+ABS(CTime_1(Sample_Offset+i_sampl)))/ABS(CTime_1(Sample_Offset)))**2 &
      ,'% of peak, TotalGain=', TotalGain, ', slicing window IntfSmoothWin=',IntfSmoothWin
   !
   Open(unit=14,STATUS='old',ACTION='read', FILE = 'files/'//TRIM(Antennas)//'_Structure_LBA.dat',IOSTAT=nxx)
   If(nxx.eq.0) Then
      write(2,*) 'opened file:', 'files/'//TRIM(Antennas)//'_Structure.dat'
   Else
      write(2,*) 'tried to open file:','files/'//TRIM(Antennas)//'_Structure.dat'
      write(2,*) 'error nr',nxx,', most probably the folder or the file needed did not exist'
      write(2,*) 'Advice: run "SimData" first.'
      goto 9
   EndIf
   !
   Call CreateNewFolder(Simulation) ! create new folder when needed
   !
   Open(unit=24,STATUS='unknown',ACTION='write', FILE = 'files/'//TRIM(Simulation)//'_Structure_LBA.dat',IOSTAT=nxx)
   If(nxx.eq.0) Then
      write(2,*) 'Writing to file:', 'files/'//TRIM(Simulation)//'_Structure.dat'
   Else
      write(2,*) 'tried to open file:','files/'//TRIM(Simulation)//'_Structure.dat'
      write(2,*) 'error nr',nxx,', most probably the folder did not exist'
      goto 9
   EndIf
   !
   call random_seed()
   Call GetSources(SrcNrMax,SourcesListLocNEh, SourcesListTms, SourcesListAmpNEh, SourcesListTS,SrcNr,NrmStI)
   !
   Station_name='CXnnn'
   SourceGuess(:)=SourcesListLocNEh(:,1)
   !AntNr_lw=1
   Do i_file=1,Station_NrMax ! 10 !80
      read(14,*,IOSTAT=nxx) STATION_ID, Station_name, N_StAnt
      If(nxx.ne.0) exit
      If(N_StAnt.gt.2*AntpStat) Then
         write(2,*) '**** mismatch antenna nr',N_StAnt,AntpStat
         N_StAnt=2*AntpStat
      EndIf
      ! write(24,*,IOSTAT=nxx) STATION_ID, Station_name
      Unt=12
      Open(unit=Unt,STATUS='old',ACTION='read', Form='unformatted', &
         FILE = 'files/'//TRIM(Antennas)//'_'//TRIM(Station_name)//'.udt')
      Call ReadSimulStatHead(Unt, SourceGuess, STATION_ID, N_StAnt,  StatStartTime, &
         Ant_Stations, Ant_IDs, Ant_pos, Ant_RawSourceDist, NrAnt, DualPosOnly)
      Close(unit=12)
      !
      Ant_nr= Ant_nr + NrAnt
      write(24,*,IOSTAT=nxx) STATION_ID, Station_name, NrAnt,  NtSamples
      !
      ! Start constructing traces
      !
      Call RelDist(SourceGuess,Ant_pos(:,1),RDist)
      StatStartTime=SourcesListTms(1)/(Sample*1000.)+RDist-200 ! place first source at sample 200 in trace
      StatStartTime=StatStartTime + TimingErr_ns * random_stdnormal()/5. ! to convert timing error to samples
      !
      !Open(unit=22,STATUS='unknown',ACTION='write', FILE = 'files/'//TRIM(Simulation)//'_'//TRIM(Station_name)//'.dat')
      OPEN( Unit=22, STATUS='unknown',ACTION='write', Form='unformatted', &
         FILE = 'files/'//TRIM(Simulation)//'_'//TRIM(Station_name)//'.udt')
      ! Structure of this unformatted file
      ! rec 1: StatStartTime [ms],  NtSamples, NrAnt [=2xnumber of different dipole positions]
      !  followed by 'j_ant=1:ant_nr/2' data for each dipole
      ! recs 2: 'Dual       ', Ant_IDs(j_ant), Ant_pos(:,j_ant) ! info for dipole 'j_ant'
      !  followed by 'i_ant=1:ant_nr' data for each antenna
      ! recs 3: time trace even antenna (only if 'dual' or 'Even', real*4) for antenna 'i_ant'
      ! at present all dipoles are necessarily 'Dual'
      !
      !write(22,*) 'StartTime[ms]=', StatStartTime*(Sample*1000.), 'N_samples=', NtSamples, &
      !      'N_ant=', Ant_nr, 'P_noise= 1.'  ![ms],,&  =noise power
      write(22)  StatStartTime*(Sample*1000.),  NtSamples, NrAnt
      write(2,"(A,I3,I9,1x,A5,A,I2,A,A,F10.5,A,F9.6)") 'Station=',i_file,STATION_ID, Station_name, &
         ' uses',NrAnt/2,' antenna pairs.', ' Start time=',StatStartTime*(Sample*1000.),'[ms]'
      Flush(unit=2)
      Do j_ant=1,NrAnt/2  ! Calculate and write spectra
         i_ant=2*j_ant - 1
         !write(22,*) Ant_IDs(j_ant), 'NEh=', Ant_pos(:,j_ant)
         write(22) 'Dual', Ant_IDs(i_ant), Ant_pos(:,i_ant)
         !write(2,*) i_ant, Ant_IDs(i_ant), Ant_pos(:,i_ant)
      Enddo
      Do j_ant=1,NrAnt/2  ! Calculate and write spectra
         i_ant=2*j_ant - 1
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
         !
         ! Get Zenithal antenna function
         Thet_d =0.
         Phi_d =0.
         Call AntFun(thet_d ,Phi_d ) ! sets J_0p,J_0t,J_1p,J_1t;  Jones; gives _p & _t polarized fields when contracted with (0.1) voltages
         !
         IstN=sqrt((1.-FracGalacNoisePow)*NnuSamples/(inu2-inu1))
         GN=(sqrt((FracGalacNoisePow)*NnuSamples/(inu2-inu1)))/sqrt(0.0273) ! *4.62 ! phenomenological factor to have similar power in GalNoise as in InstNoise (when GalacNoise=1.)
         !write(2,*) 'FracGalacNoisePow', FracGalacNoisePow, IstN, GN, (inu2-inu1)*1./NnuSamples
         GN= GN*(inu2+inu1)*dnu/TotalGain
         !write(2,*) 'GalacNoisePow', (inu2+inu1),dnu,TotalGain, (inu2+inu1)*dnu/TotalGain, (inu2-inu1)*1./NnuSamples
         !  GalacNoisePow        2252   4.8828125000000000E-002   17.869174845587548        6.1536662129170869       0.699218750
         !
         Cnu_0(0:inu1-1)=0.
         Cnu_0(inu2+1:NnuSamples)=0.
         Cnu_1(0:inu1-1)=0.
         Cnu_1(inu2+1:NnuSamples)=0.
         Do i_nu=inu1,inu2   ! Pull through antenna function with 1/nu frequency spectrum
            nu=i_nu*dnu
            i_freq=Int(nu)
            dfreq=nu-i_freq
            Cnu_0(i_nu)=(IstN*Cnu_0(i_nu) + &
                        ((1.-dfreq)*J_0p(i_freq) + dfreq*J_0p(i_freq+1)) * Cnu_p(i_nu)*GN/nu + &
                        ((1.-dfreq)*J_0t(i_freq) + dfreq*J_0t(i_freq+1)) * Cnu_t(i_nu)*GN/nu)
            Cnu_1(i_nu)=(IstN*Cnu_1(i_nu) + &
                        ((1.-dfreq)*J_1p(i_freq) + dfreq*J_1p(i_freq+1)) * Cnu_p(i_nu)*GN/nu + &
                        ((1.-dfreq)*J_1t(i_freq) + dfreq*J_1t(i_freq+1)) * Cnu_t(i_nu)*GN/nu)
         Enddo
         BckgrPwr_0=BckgrPwr_0 + SUM(ABS(Cnu_0(:))**2)/2.
         BckgrPwr_1=BckgrPwr_1 + SUM(ABS(Cnu_1(:))**2)/2.  !=mean power per time sample
         N_even=N_even +1
         N_odd=N_odd +1
         !write(2,*) 'Ave background power=', BckgrPwr_0/N_even, BckgrPwr_1/N_odd !  (should be about 1. to be consistent with Ant-read.f90)
         !
         !stop
         !Cnu_0(:)=0.
         !Cnu_1(:)=0.
         Do i_src=1,SrcNr
            !write(2,*) 'Galactic background=',SUM(ABS(Cnu_p(:))**2),SUM(ABS(Cnu_t(:))**2)
            Call RelDist(SourcesListLocNEh(:,i_src),Ant_pos(:,i_ant),RDist)
            T_Offset=RDist - StatStartTime +SourcesListTms(i_src)/(Sample*1000.)  ! in units of samples
            Sample_Offset = INT(T_Offset) ! in units of sample size
            !write(2,*) 'Sample_Offset for src=', i_src, Sample_Offset, &
            !      ' ,antenna=', Ant_IDs(j_ant), 'NEh=', Ant_pos(:,j_ant)
            If(Sample_Offset.le.0 .or. Sample_Offset.ge.NtSamples) Then
               write(2,*) '*****Sample_Offset out of range for src=', i_src, Sample_Offset, &
                     ' ,antenna=', Ant_IDs(i_ant), 'NEh=', Ant_pos(:,i_ant)
               cycle
            EndIf
            SubSample_Offset = T_Offset - Sample_Offset ! in units of sample size
            RTime_s(:)=0.
            RTime_s(Sample_Offset)=1.
            Call RFTransform_CF(RTime_s,Cnu_s(0)) ! delta peak at abou right time, filtering done while selecting proper frequencies.
            !
            ! get t and p oriented signals
            !
            Ras(:)=(SourcesListLocNEh(:,i_src)-Ant_pos(:,i_ant))/1000.  ! \vec{R}_{antenna to source}
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
            If((j_ant.eq.3) .and.  (i_src.eq. -1)) Then
               write(2,*) 'thet_d,phi_d', thet_d,phi_d, Ant_pos(:,j_ant), A_p/A_t
               i_freq=(Freq_min + Freq_max)/2
               write(2,*) i_freq, 'p', Ji_p0(i_freq)*J_0p(i_freq)+Ji_p1(i_freq)*J_1p(i_freq), &
                  Ji_t0(i_freq)*J_0p(i_freq)+Ji_t1(i_freq)*J_1p(i_freq),  &
               't', Ji_t0(i_freq)*J_0t(i_freq)+Ji_t1(i_freq)*J_1t(i_freq), Ji_p0(i_freq)*J_0t(i_freq)+Ji_p1(i_freq)*J_1t(i_freq)
            EndIf
            !  for sub-sample shift of peak position:
            dphShift=exp(ipi*SubSample_Offset/NnuSamples )
            phShift=exp(ipi*inu1*SubSample_Offset/NnuSamples )
            Do i_nu=inu1,inu2   ! Increment frequency spectrum of the antenna with the signal from this source
               nu=i_nu*dnu
               i_freq=Int(nu)
               dfreq=nu-i_freq
               CX=phShift*Cnu_s(i_nu)
               Cnu_0(i_nu)=Cnu_0(i_nu) + &
                           ((1.-dfreq)*J_0p(i_freq) + dfreq*J_0p(i_freq+1)) * A_p*CX + &
                           ((1.-dfreq)*J_0t(i_freq) + dfreq*J_0t(i_freq+1)) * A_t*CX
               Cnu_1(i_nu)=Cnu_1(i_nu) + &
                           ((1.-dfreq)*J_1p(i_freq) + dfreq*J_1p(i_freq+1)) * A_p*CX + &
                           ((1.-dfreq)*J_1t(i_freq) + dfreq*J_1t(i_freq+1)) * A_t*CX
               phShift=phShift*dphShift
            Enddo
            !write(2,*) 'Src=',i_src,SUM(ABS(Cnu_0(:))**2)*NtSamples/PulseRespWidth, &
            !      SUM(ABS(Cnu_1(:))**2)*NtSamples/PulseRespWidth, A_p, A_t, Sample_Offset
            !stop
         EndDo ! i_src=1,SrcNr
         !
         Call RFTransform_CF2RT(Cnu_0(0),FTime_0(1,j_ant) )
         Call RFTransform_CF2RT(Cnu_1(0),FTime_1(1,j_ant) )
         !write(2,*) '!writing:', NtSamples, j_ant
         !write(2,*) '!before time trace:', FTELL(22)
         write(22) REAL(NormEvenOdd *FTime_0(1:NtSamples,j_ant),4)
         !write(2,*) '!between time trace:', FTELL(22)
         write(22) REAL(NormEvenOdd *FTime_1(1:NtSamples,j_ant),4)
         !write(2,*) '!after time trace:', FTELL(22)
         If(i_file.eq.1 .and. j_ant.eq.1) Then
            Call RFTransform_CF2CT(Cnu_0(0),CTime_0(1) )  ! RFTransform_CF2RT(Cnu,RD)
            Call RFTransform_CF2CT(Cnu_1(0),CTime_1(1) )
            write(2,*) '* background power/sample=', BckgrPwr_0, BckgrPwr_1 !  (should be about 1. to be consistent with Ant-read.f90)
            write(2,*) '* total power in background=', NtSamples*BckgrPwr_0, NtSamples*BckgrPwr_1
            write(2,*) '* power in sources only:', &
               NtSamples*(SUM(ABS(Cnu_0(:))**2)/2.-BckgrPwr_0), NtSamples*(SUM(ABS(Cnu_1(:))**2)/2.-BckgrPwr_1)
               ! this is the first cycle so N_odd=N_even=1
            write(2,*) '* total power, backgr + sources=',SUM(FTime_0(:,j_ant)**2),SUM(FTime_1(:,j_ant)**2)
            Open(unit=23,STATUS='unknown',ACTION='write', FILE = 'files/'//TRIM(Simulation)//'_RefAntTrace.dat')
            write(2,*) 'spectrum of ant=',j_ant,' written to:','files/'//TRIM(Simulation)//'_RefAntTrace.dat'
            Do i_sampl=1,NtSamples
               write(23,*) i_sampl,REAL(CTime_0(i_sampl)), ABS(CTime_0(i_sampl)), REAL(CTime_1(i_sampl)), ABS(CTime_1(i_sampl))
            Enddo
            close(unit=23)
         EndIf
      Enddo ! j_ant=1,NrAnt
      !
      !NormEvenOdd   ! to undo the action of NormEven & NormOdd in AntennaRead
      !Do i_sampl=1,NtSamples
      !   !write(22,*) (NormEvenOdd *FTime_0(i_sampl,j_ant), NormEvenOdd *FTime_1(i_sampl,j_ant), j_ant=1,Ant_nr)
      !   write(22) (NormEvenOdd *FTime_0(i_sampl,j_ant), NormEvenOdd *FTime_1(i_sampl,j_ant), j_ant=1,Ant_nr)
      !   !write(2,*) i_sample, Sample_Offset
      !Enddo
      Close(unit=22)
      !
      !AntNr_lw=   AntNr_lw+NrAnt
   Enddo !  i_file=1,StationNrMax
   write(2,*) 'Ave background power (all antennas)=', BckgrPwr_0/N_even, BckgrPwr_1/N_odd, N_even, N_odd, NtSamples !  (should be about 1. to be consistent with Ant-read.f90)
   Close(unit=14)
   Close(unit=24)
   Call DAssignFFT()
   !
9  Continue
   Stop "normal ending of 'Simulate_Data'"
End Program Simulate_Data
!=====================================
Subroutine GetSources(SrcNrMax,SourcesListLocNEh, SourcesListTms, SourcesListAmpNEh, SourcesListTS,SrcNr, NrmStI)
   use constants, only : dp,sample_ms,pi, c_mps ! ,Refrac,pi,ci
   !Use RandomFunc, only : random_stdnormal, random_stdnormal3D
   Implicit none
   !
   Integer, intent(in) :: SrcNrMax
   Real(dp), intent(out) :: SourcesListLocNEh(3,SrcNrMax), SourcesListTms(SrcNrMax), SourcesListAmpNEh(3,SrcNrMax)
   Real(dp), intent(out) :: SourcesListTS(SrcNrMax)  !  Souce time in samples
   Integer, intent(out) :: SrcNr
   Real(dp), intent(in) :: NrmStI
   !
   Character(len=12) :: Lab1  !, Lab2, Lab3, Lab4  Station_name,
   !Real(dp) :: StatStartTime, SubSample_Offset, LFRAnt_crdnts(3), RDist, T_Offset, Powr !,StatAnt_Calib !,StartTime_ms
   Real(dp) :: Time_Interval, Sources_TS, Sources_LocNEh(1:3), Sources_AmpNEh(1:3),t_shft
   Integer :: NConstruct, NCloud, nxx, i_src, i_rep
   Real(dp) :: D, Phi, Thet, Time_width, Space_width, PolSpread, R,Epsln=epsilon(Epsln)
   Character(LEN=180) :: lname
   Real :: x,y,z, RGV(1:3) ! Random Gaussian distributed Vector
   Real :: random_stdnormal
   Real(dp) :: IndxRefrac, I12, I123
   Real(dp), external :: RefracIndex
   !
   SrcNr=0
   Do !i_src=1,SrcNrMax
      If(SrcNr.eq.SrcNrMax) exit
      Call GetNonZeroLine(lname)
      write(2,"(A,i4,A,A,A)") 'input line for source',SrcNr+1,':"',TRIM(lname),'"'
      Read(lname,*,IOSTAT=nxx) Lab1
      If(nxx.ne.0) exit
      ! Check for repetitive source in time
      If((Lab1(1:3) .eq. 'Rep') .or. (Lab1(1:3) .eq. 'rep')) Then  ! Repetitive
         Read(lname,*,IOSTAT=nxx) Lab1, Time_Interval, NConstruct
         If(nxx.ne.0) then
            Write(2,*) 'Expected: "Rep..., Time_Interval, NConstruct" but got ',TRIM(lname)
            stop
         Endif
         !
         Call GetNonZeroLine(lname)
         write(2,"(A,i4,A,A,A)") 'input line for Rep-source',SrcNr+1,':"',TRIM(lname),'"'
         Read(lname,*,IOSTAT=nxx) Lab1
         If((Lab1(1:3) .eq. 'Clo') .or. (Lab1(1:3) .eq. 'clo')) Then ! Cloud
            Read(lname,*,IOSTAT=nxx) Lab1, Time_width, Space_width, PolSpread, NCloud
            If(nxx.ne.0) then
               Write(2,*) 'Expected: "Clo..., Time_width[samples], Space_width, PolSpread, NConstruct" but got ',TRIM(lname)
               stop
            Endif
            Call GetNonZeroLine(lname)
            write(2,"(A,i4,A,A,A)") 'input line for Rep-Cloud-source',SrcNr+1,':"',TRIM(lname),'"'
            Read(lname,*,IOSTAT=nxx) Sources_TS, Sources_LocNEh(:), Sources_AmpNEh(:)
            If(nxx.ne.0) then
               Write(2,*) 'Expected: "Sources_TS[samples], Sources_LocNEh(1:3)[m], Sources_AmpNEh(1:3)" but got ',TRIM(lname)
               stop
            Endif
            Call Convert2m(Sources_LocNEh)
            If(NCloud*NConstruct.gt.(SrcNrMax-SrcNr)) NCloud=(SrcNrMax-SrcNr)/NConstruct
            write(2,"(A,I5,A, 3F9.3,A, 3F9.3)") 'Form a cloud of',NCloud,' sources around:', Sources_LocNEh(:)/1000., &
               '[km] with strength=', Sources_AmpNEh(:)
            write(2,*) 'Spread in time, pos, pol=',Time_width, Space_width, PolSpread
            Do i_rep=0,NConstruct-1
               Do i_src=0,NCloud-1
                  SrcNr=SrcNr+1
                  write(2,"(A,F9.1,A)") '   @ central time=',Sources_TS,' samples'
                  SourcesListTS(SrcNr)=Sources_TS + random_stdnormal()*Time_width
                  Call random_stdnormal3D(RGV)
                  SourcesListLocNEh(:,SrcNr)=Sources_LocNEh(:)+ Space_width*RGV(:)
                  Call random_stdnormal3D(RGV)
                  SourcesListAmpNEh(:,SrcNr)=Sources_AmpNEh(:) + PolSpread*RGV(:)
                  ! write(3,*) i_src,D,SourcesListLocNEh(:,SrcNr) ! just for checking the cloud structure using  Sources_Plots.gle
               EndDo
               Sources_TS=Sources_TS + Time_Interval
            EndDo
         Else
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
      ElseIf((Lab1(1:3) .eq. 'Clo') .or. (Lab1(1:3) .eq. 'clo')) Then ! Cloud
         Read(lname,*,IOSTAT=nxx) Lab1, Time_width, Space_width, PolSpread, NConstruct
         If(nxx.ne.0) then
            Write(2,*) 'Expected: "Clo..., Time_width[samples], Space_width, PolSpread, NConstruct" but got ',TRIM(lname)
            stop
         Endif
         Call GetNonZeroLine(lname)
         write(2,"(A,i4,A,A,A)") 'input line for Cloud-source',SrcNr+1,':"',TRIM(lname),'"'
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
            !SourcesListLocNEh(3,SrcNr)=Sources_LocNEh(3)+ D*(0.5-z)/R
            Call random_stdnormal3D(RGV)
            SourcesListLocNEh(:,SrcNr)=Sources_LocNEh(:)+ Space_width*RGV(:)
            Call random_stdnormal3D(RGV)
            SourcesListAmpNEh(:,SrcNr)=Sources_AmpNEh(:) + PolSpread*RGV(:)
            ! write(3,*) i_src,D,SourcesListLocNEh(:,SrcNr) ! just for checking the cloud structure using  Sources_Plots.gle
         EndDo
         cycle
      Else
         ! Check for single source, otherwise quit
         Read(lname,*,IOSTAT=nxx) Sources_TS, Sources_LocNEh(:), Sources_AmpNEh(:)
         !write(2,*) '!nxx=',nxx,Sources_TS, Sources_LocNEh(:), Sources_AmpNEh(:)
         If(nxx.ne.0) Then
            exit  ! terminate reading sequence
         EndIf
         Call Convert2m(Sources_LocNEh)
         SrcNr=SrcNr+1
         SourcesListTS(SrcNr)=Sources_TS
         SourcesListLocNEh(:,SrcNr)=Sources_LocNEh(:)
         SourcesListAmpNEh(:,SrcNr)=Sources_AmpNEh(:)
      EndIf
      ! constants tunes as to have background (sum square freq. spectrum) =1.
      ! source at distance D [km] with intensity=D gives power per 20 samples = 1 = backgrnd
   EndDo
   If(SrcNr.le.0) Then
      write(2,*) 'No valid sources have been entered:',lname
   EndIf !
   write(2,"(T4,'#',T7,'t_r[ms]',T18,'t[ms]',T30,'(N,E,h) [km]', T61,'Ampl(N,E,h)',T87,'I123_TRI-D')")
   Do i_src=1,SrcNr
      !t_shft=sqrt(SUM(SourcesListLocNEh(:,i_src)*SourcesListLocNEh(:,i_src)))*1000.*Refrac/c_mps ! in [ms] due to signal travel distance
      t_shft=sqrt(SUM(SourcesListLocNEh(:,i_src)*SourcesListLocNEh(:,i_src)))*1000.*RefracIndex(SourcesListLocNEh(3,i_src))/c_mps ! in [ms] due to signal travel distance
      !write(2,*) i_src,t_shft
      x=sqrt(SUM(SourcesListLocNEh(:,i_src)*SourcesListLocNEh(:,i_src)))
      z=sqrt(SUM(SourcesListAmpNEh(:,i_src)*SourcesListAmpNEh(:,i_src)))
      y=SUM(SourcesListLocNEh(:,i_src)*SourcesListAmpNEh(:,i_src))/(x*z)
      thet=ACOS(y)*180./pi
      SourcesListTms(i_src)=SourcesListTS(i_src)*Sample_ms+t_shft  ! source time at core
      I123=SUM(SourcesListAmpNEh(:,i_src)*SourcesListAmpNEh(:,i_src))*NrmStI
      write(2,"(I4,F11.5,F11.5,2x, 3F10.4,2x,3F8.2,2x,F12.2, 2x, F8.2)") i_src,(SourcesListTms(i_src)-SourcesListTms(1)), &
         SourcesListTS(i_src)*Sample_ms, SourcesListLocNEh(:,i_src)/1000., SourcesListAmpNEh(:,i_src), I123, thet
   Enddo
   write(2,*) 'Nr of sources=',SrcNr, ', time of first source pulse at the core=',SourcesListTms(1)
   Flush(unit=2)
End Subroutine GetSources
