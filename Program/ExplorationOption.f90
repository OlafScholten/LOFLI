! ---------------------------------
Subroutine ExplorationRun
   use constants, only : dp,sample
   use DataConstants, only : Time_dim, DataFolder, FlashName
   use ThisSource, only : XFrameEi, XFrameEf, XFrameNi, XFrameNf, XFrameh
   use Chunk_AntInfo, only : StartT_sam, TimeFrame
   use FitParams, only : SpaceCov
   use FitParams, only : PeakS_dim, MaxFitAntDistcs, MaxFitAntD_nr !,
   use Explore_Pars, only : PeakS, StatsStore_SAI, StatsStore_Ave, StatsStore_RMS, StatsStore_Peak, N_ExplTimes
   use Explore_Pars, only : Alloc_Explore_Pars
   use FFT, only : RFTransform_su,DAssignFFT
   use GLEplots, only : GLEplotControl
   !use StationMnemonics, only : Statn_ID2Mnem, Statn_Mnem2ID
   Implicit none
   !
   Integer :: j,i_chunk, ChunkNr_start, ChunkNr_stop, units(0:2), FitRange_Samples, i_expl
   Real*8 :: StartTime_ms, StartingTime, StoppingTime, dt_ms
   Real*8 :: SourceGuess(3,1) ! = (/ 8280.01,  -15120.48,    2618.37 /)     ! 1=North, 2=East, 3=vertical(plumbline); dimension needed for "Call FindCallibr(SourceGuess)"
   Integer :: ir_file,ir_grp,ir_dst, DATA_LENGTH, SAMPLE_NUMBER_first, i_dist, i_guess, TimeFr
   Real*8 :: LFRAnt_crdnts(3), powr, T_Offset !,StatAnt_Calib !,StartTime_ms
!
   i_chunk=1
   SourceGuess(:,1)=(/40000.,40000.,5000./)
   units(:)=-1
   PeakS_dim=10
   TimeFr=0
   Read (12) ir_file,ir_grp,ir_dst, LFRAnt_crdnts, powr, DATA_LENGTH, SAMPLE_NUMBER_first!, Absolute_TIME, DIPOLE_CALIBRATION_DELAY,nu_fltr !nu_fltr
   StartTime_ms= 1.d0*NINT(SAMPLE_NUMBER_first *sample*1000.d0)
   StoppingTime=StartTime_ms+DATA_LENGTH*sample*1000.d0
   write(2,*) 'start & stop times=',StartTime_ms,StoppingTime
   !StartTime_ms=784.
   XFrameEi=-80  ; XFrameEf=+80
   XFrameNi=-80  ; XFrameNf=+80
   XFrameh=15
   OPEN(UNIT=29,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//'Explore.dat')
   write(29,"(6(1x,f8.3),2(1x,f8.3),1x,A20)") XFrameEi, XFrameEf, XFrameNi, XFrameNf, &
      0.,XFrameh,StartTime_ms,StoppingTime,' NoBox  0 ! Explore '
   Write(29,*) ' 0 1 0 0 Explore-', TRIM(FlashName), ' 0 0 0 0 0 0 0 0 0.1 !'
        !If(Ant_nr(1).eq.0) write(2,*) ', TimeOffset', T_Offset, Dset_Offset,DATA_LENGTH, Sample_Offset - SAMPLE_NUMBER_first
   dt_ms=30.
   TimeFr=0
   StartTime_ms=StartTime_ms + dt_ms*TimeFr
   !StoppingTime=2000
   N_ExplTimes=  INT((StoppingTime - StartTime_ms)/dt_ms)
   StartTime_ms=StartTime_ms - dt_ms*2/3.
   Call Alloc_Explore_Pars
   Do i_expl=1,N_ExplTimes
      write(2,"(50('= '))")
      StartTime_ms=StartTime_ms + dt_ms
      StartT_sam(i_chunk)=(StartTime_ms/1000.d0)/sample  ! in sample's
      If(StartTime_ms.gt.StoppingTime) exit
      Do i_guess=0,3
         TimeFrame=TimeFr+ 10*i_expl + i_guess    ! just a label
         j=MOD(i_guess,2)
         SourceGuess(1,i_chunk)= (2*j-1)* 40000.  ! [km]  ,  toggle between +/- 40 km
         SourceGuess(2,i_chunk)= (i_guess-j-1)* 40000.  ! [km]  ,  toggle between +/- 40 km
         Call RFTransform_su(Time_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         Call AntennaRead(i_chunk,SourceGuess(:,i_chunk))
         Call DAssignFFT()
         !Call Find_unique_StatAnt()
         Call FindCallibr(SourceGuess)
         !Write(2,*) 'number of close peaks in even&odd:',PeakD_nr
         !Call SourceFind(TimeFrame,SourceGuess,units)
      Enddo !  i_guess=1,4
      Call GetSpectStats(i_chunk, i_expl)
      !stop 'Explore test'
   enddo
   Close(unit=29)
   !
   Call GLEplotControl(PlotType='SourcesPlot', PlotName='Img_Explore', &
      PlotDataFile=TRIM(DataFolder)//'Explore', Submit=.true.)
   Return
   !Call AnalSpectStats  ! Obsolete now, a more thorough job is done in RFI_Mitigation-v18
End Subroutine ExplorationRun
!=====================================
Subroutine GetSpectStats(i_chunk, i_expl)
   use DataConstants, only : Station_nrMax,PeakNr_dim,Time_dim, Production
   use Chunk_AntInfo, only : CTime_spectr, Ant_Stations, Ant_IDs, Ant_nr !
   use Chunk_AntInfo, only : Unique_SAI, Unique_StatID, Nr_UniqueStat
   use constants, only : dp,pi,ci,sample
   use StationMnemonics, only : Station_ID2Mnem, Statn_ID2Mnem
   use Explore_Pars, only : PeakS, StatsStore_SAI, StatsStore_Ave, StatsStore_RMS, StatsStore_Peak
   Implicit none
   Integer, intent(in) ::  i_chunk, i_expl
   !
   Integer ::  Edge, i_ant, PSP, Wu, Wl
   Integer :: j, i_loc(1), SAI
   Real(dp) :: HEnvel(Time_dim), PeakSAmp(0:PeakS), Ave, RMS
   Real(dp), parameter :: safe=1.d-8 ! Almost =0, used for bookkeeping of signals found

   integer :: i_SAI, StLoc, t_Max,i_stat, Antenna_SAI
   Integer, save :: PeakNr1, i_stat_old=0, N_ant
   character(len=5) :: Station_Mnem
   logical :: prnt= .false.
   !
   !
   Edge=100
   Do i_ant=1,Ant_nr(i_chunk)
      HEnvel(:)=abs(CTime_spectr(:,i_ant, i_chunk))
      ! write(label,"(I2.2,I3.3)") TimeFrame,Ant_IDs(i_ant,i_chunk)
      Do j=1,PeakS
         i_loc=MaxLoc( HEnvel(Edge:Time_dim-Edge) )
         PSP=i_loc(1) + Edge - 1
         !write(*,*) 'j=',j,psp
         !Write(2,*) 'PSP',j,i_peak,PSP
         !Write(2,"(11F12.1)")  HEnvel(PSP- 10:PSP + 10)
         Wu=20  ! determine the position of the min in the spectrum above the peak
         Wl=20 ! determine the position of the min in the spectrum below the peak
         StatsStore_Peak(i_expl,i_ant,j)=HEnvel(PSP)
         HEnvel(PSP - Wl:PSP + Wu )=-safe  ! zero the test spectrum around the peak just found over an interval that matches that of "CleanPeak"
      Enddo   !  j=1,2*PeakS_dim
      !Write(2,"(A,i3,A)", ADVANCE='NO') 'Peak(',i_peak,'):'
      !Write(2,"(20i6)") PeakSP(1:i_peak,i_eo)
      !Write(2,"(A)", ADVANCE='NO') 'Peakval = '
      !Write(2,"(20f6.0)") HEnvel(PeakSP(1:i_peak,i_eo))
      !?? should   nu_fltr  be accounted for?
      HEnvel(:)=abs(CTime_spectr(:,i_ant, i_chunk))
      Ave=SUM( HEnvel(Edge:Time_dim-Edge) ) /(Time_dim-2*Edge)
      RMS=SUM( HEnvel(Edge:Time_dim-Edge)*HEnvel(Edge:Time_dim-Edge) ) /(Time_dim-2*Edge)
      SAI= 1000*Ant_Stations(i_ant,i_chunk) + Ant_IDs(i_ant,i_chunk)
      StatsStore_SAI(i_expl,i_ant)=SAI
      StatsStore_Ave(i_expl,i_ant)=Ave
      StatsStore_RMS(i_expl,i_ant)=SQRT(RMS)
      !StatsStore_Peak(i_expl,i_ant,1:PeakS)=PeakSAmp(1:PeakS)
      !
   EndDo ! i_ant=1,Ant_nr(i_chunk)
   ! norm =2int x e^{x^2/s^2} = 2 s sqrt(2/pi) = N_samples  ! note 2 parameters = real & imag
   ! ave= (1/norm) 2 int x^2 e^{x^2/s^2} = 2 s^2/(2 s sqrt(2/pi))= s sqrt (pi/2)
   ! RMS^2= (1/norm) 2 int x^3 e^{x^2/s^2} = 2 s^3 2 sqrt(2/pi)/(2 s sqrt(2/pi))= 2 s^2
   ! RMS = sqrt(2) s
   ! ave/RMS = sqrt(pi/4)
   !Chance amplitude > 4 sigma = 1/15787.2 !	
   !with 65536-2*Edge samples: numerically: ~ 4 times amplitude gt 3*RMS=4.24 s
   !write(2,*) 'stats:',i_expl,StatsStore_SAI(i_expl,1),StatsStore_Ave(i_expl,1),StatsStore_RMS(i_expl,1) &
   !   ,StatsStore_Peak(i_expl,1,1:PeakS)
   Return
End Subroutine GetSpectStats
! =============================
Subroutine AnalSpectStats
   use constants, only : dp,pi
   use Explore_Pars, only : PeakS, StatsStore_SAI, StatsStore_Ave, StatsStore_RMS, StatsStore_Peak, N_ExplTimes, Ant_nrMax
   use Chunk_AntInfo, only : Unique_SAI, Nr_UniqueAnt, Ant_nrMax
   Implicit none
   Integer :: i_expl, i_ant
   Integer :: i_SAI, i_loc(1), SAI, N_bkgr, i_eo, N_eo(0:1)
   Real(dp) :: Ave, RMS, sigma, RMS_bkgr, RMS_eo(0:1), RMS_SAI(1:Ant_nrMax)
   !
   write(2,*) "i_expl,  i_ant,  sigma, Ave/(sigma*SQRT(pi/2)) ,  StatsStore_Peak(i_expl,i_ant,1:PeakS)/sigma",&
      PeakS,"=PeakS"
   RMS_eo(:)=0.
   N_eo(:)=0.
   Do i_SAI=1,Nr_UniqueAnt
      SAI=Unique_SAI(i_SAI)
      N_bkgr=0
      RMS_bkgr=0.
      !Write(2,*) 'summary stats:',N_ExplTimes,SAI
      Do i_expl=1,N_ExplTimes
         i_loc=MAXLOC(StatsStore_SAI(i_expl,:), MASK = StatsStore_SAI(i_expl,:) .eq. SAI)
         If(i_loc(1).gt.0) then
            i_ant=i_loc(1)
            Ave=StatsStore_Ave(i_expl,i_ant)
            RMS=StatsStore_RMS(i_expl,i_ant)
            sigma=RMS/sqrt(2.)
            ! ave/RMS = sqrt(pi/4) in case of real noise, in addition (StatsStore_Peak(i_expl,i_ant,1)/sigma .lt. 5)
            If(StatsStore_Peak(i_expl,i_ant,1)/sigma .lt. 5.) Then
               N_bkgr=N_bkgr+1
               RMS_bkgr=RMS_bkgr+RMS
               !write(2,*) i_expl,i_ant,sigma, Ave/(sigma*SQRT(pi/2)) ,StatsStore_Peak(i_expl,i_ant,1:PeakS)/sigma
            EndIf
         EndIf ! i_ant=SAI
      EndDo ! i_expl=1,N_ExplTimes
      If(N_bkgr .lt. 3) Then
         Write(2,*) 'For',SAI, 'too few backgrouds=', N_bkgr
         RMS_SAI(i_SAI)=-1.
      Else
         !Write(2,*)'RMS_bkgr/N_bkgr',RMS_bkgr/N_bkgr,N_bkgr
         RMS_SAI(i_SAI)=RMS_bkgr/N_bkgr
         i_eo=MOD(SAI,2)
         RMS_eo(i_eo)=RMS_eo(i_eo) +RMS_bkgr/N_bkgr
         N_eo(i_eo)=N_eo(i_eo) +1
      EndIf
   EndDo ! j=1,Nr_UniqueAnt
   Do i_eo=0,1
      RMS_eo(i_eo)=RMS_eo(i_eo)/N_eo(i_eo)
   Enddo
   Write(2,*) 'Background',RMS_eo(:),N_eo(:)
   Do i_SAI=1,Nr_UniqueAnt
      SAI=Unique_SAI(i_SAI)
      i_eo=MOD(SAI,2)
      If( (RMS_SAI(i_SAI).gt.RMS_eo(i_eo)*1.1) .or. (RMS_SAI(i_SAI).lt.RMS_eo(i_eo)/1.1) ) Then
      Write(2,*) 'For SAI=',SAI, 'more than 10% difference in RMS',RMS_SAI(i_SAI),RMS_eo(i_eo)
      Endif
      RMS_SAI(i_SAI)=RMS_SAI(i_SAI)/RMS_eo(i_eo)
   EndDo
   !
   Open(unit=12,STATUS='unknown',ACTION='write',FORM ="unformatted", FILE = 'Book/NormFactors.uft') ! trim(DataFolder)//
   Write(12) Nr_UniqueAnt,RMS_eo(:)
   Write(12) Unique_SAI(1:Nr_UniqueAnt)
   Write(12) RMS_SAI(1:Nr_UniqueAnt)
   Close(Unit=12)
   Return
End Subroutine AnalSpectStats
