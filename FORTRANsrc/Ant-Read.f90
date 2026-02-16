Subroutine AntennaRead(i_chunk,SourceGuess)
!  v18 : Normalization is changed
    ! Convention for i_ant
    ! i_ant runs from 1 till  Ant_nr(i_chunk), a sequential number for all antennas for which there are data in this chunk
    !
    ! Ant_IDs(i_ant, i_chunk) gives the antenna ID [1,94]
    ! Ant_Stations(i_ant, i_chunk)=Station_ID  is ID number for the stations =10*LOFAR-ID + (0 for LBA .or. 1 for HBA)
    !    This will distinguish HBA and LBA
    !    MODULO(Unique_StatID,10)  Gives the LOFAR-ID of the combined LBA-HBA station
    !    MOD(Unique_StatID,10)  Gives the antenna type; 0 for LBA & 1 for HBA
    ! AntType(i_ant, i_chunk) gives 0 for LBA and 1 for HBA; may not be a necessary parameter
    !
   use constants, only : dp,sample, sample_ms, HeightCorrectIndxRef, ci, pi
   use DataConstants, only : Time_dim, Diagnostics, Production, OutFileLabel, RunMode , ChunkNr_dim
   use DataConstants, only : Ant_nrMax, HAnt_nrMax, LAnt_nrMax, Used_HStationNr
   use Chunk_AntInfo, only : ExcludedStatID, StartT_sam, BadAnt_nr, BadAnt_SAI, DataReadError, AntennaNrError
   use Chunk_AntInfo, only : ExcludedStat_max, SgnFlp_nr, PolFlp_nr, SignFlp_SAI, PolFlp_SAI, BadAnt_SAI, SaturatedSamplesMax
   use Chunk_AntInfo, only : Ant_Stations, Ant_IDs, Ant_nr, LBA_nr, HBA_nr, Ant_NrMax, Ant_pos, CalibratedOnly
   use Chunk_AntInfo, only : CTime_spectr, CTime_Hspectr, Ant_RawSourceDist, Simulation, WriteSimulation
   use Chunk_AntInfo, only : NormOdd, NormEven, NormHOdd, NormHEven ! Powr_eo,NAnt_eo
   !use StationMnemonics, only : Statn_ID2Mnem, Statn_Mnem2ID
   use HDF5_LOFAR_Read, only : filename, Group_names, Group_nr ! , Group_max
   use HDF5_LOFAR_Read, only : DSet_names, DSet_nr, Ant_ID, STATION_ID !, DSet_max, ANTENNA_POSITION
   use HDF5_LOFAR_Read, only : DATA_LENGTH, SAMPLE_NUMBER_first, Absolute_TIME, DIPOLE_CALIBRATION_DELAY
   !use HDF5_LOFAR_Read, only : GetFileName, ListGroups, ListDataAtt, ListGroupStructure
   use HDF5_LOFAR_Read, only : GetDataChunk
   Use StationMnemonics, only : Station_ID2Mnem, Statn_ID2Mnem
   use FFT, only : RFTransform_CF, RFTransform_su, DAssignFFT, RFTransform_CF2CT, RFTransform_CF_Filt,Hann
   Use Calibration, only : Station_ID2Calib
   Use FitParams, only : Alloc_FitParams
   use ThisSource, only : Alloc_ThisSource
   use GLEplots, only : GLEplotControl
   use CPU_timeUsage, only : CPU_usage
   Implicit none
   Integer, intent(in) :: i_chunk
   Real*8, intent(in) :: SourceGuess(3) ! = (/ 8280.01,  -15120.48,    2618.37 /)     ! 1=North, 2=East, 3=vertical(plumbline)
   !
   logical :: file14open=.false., file15open=.false.
   Integer :: U_bk, U_flt
   Integer*2 :: Chunk(1:Time_dim)
   Real(dp) :: RTime_s(1:Time_dim)
   Complex(dp) :: CNu_s(0:Time_dim/2), CTime_s(1:2*Time_dim)
   Complex(dp) :: HCNu_s(0:Time_dim), ph_shft
   integer :: Cnu_dim, Dset_offset, Sample_Offset, DataReadErr, NAnt_eo(0:1, 0:1)
   Real*8 :: Powr, Powr_eo(0:1,0:1), NormEvenOdd=1.d2
   Real(dp) :: nu_Fltr(0:Time_dim/2) !, Av, Bv, FiltFact
   !Integer, save :: WallCount
   !Real,save :: CPUstartTime, CPUstopTime=-1., WallstartTime=0, WallstopTime=-1.
   !Real,save :: TotCPURead=0., TotCPUFit=0., TotWallRead=0., TotWallFit=0.
   !
   Integer :: j,NZero, N_one, ChMax, ChMin,i_file,i_grp,i_dst, ir_file,ir_grp,ir_dst, i_type !, it_min, it_max
   Integer :: i_ant, i_time, i_eo, unt, nxx, StAntID, sgn,  AntNr_lw, AntNr_up
   Real*8 :: SubSample_Offset, LFRAnt_crdnts(3), RDist, T_Offset,StatAnt_Calib !,StartTime_ms
   Logical :: Dubbel,test, Calibrated
   Logical, save :: FirstPass=.true., HBA, LBA
   character*8 :: Ant_name
   !Logical ::    QuickCurtain=.true.
   Logical ::    QuickCurtain=.false.
   Integer ::   i_sample, i, nu_i, nu_f
   Real(dp) ::   Location(1:3), X, Av,Bv
   Character(len=6) :: Station_Mnem
   !
   Cnu_dim=Time_dim/2
   ! time recording
   !call cpu_time(CPUstartTime)
   !CALL SYSTEM_CLOCK(WallCount, powr)  !  Wallcount_rate)
   Call CPU_usage(Message='start reading' )
   !
   !Diagnostics=.true.
   !
   !WallstartTime=WallCount/powr
   !If(WallstopTime.lt.0) WallstopTime=WallstartTime
   !If(CPUstopTime.lt.0) CPUstopTime=CPUstartTime
   !WRITE(*,"(A,F9.3,A)") 'Wall clock time:',WallstartTime, '[s]'  ! count_rate, count_max
   !TotCPUFit=TotCPUFit - CPUstopTime + CPUstartTime
   !TotWallFit=TotWallFit - WallstopTime +WallstartTime
   !WRITE(*,"(A,F9.3,F12.6,A)") 'Totals Fit, Wall & CPU:',TotCPUFit, TotWallFit, '[s]'  ! count_rate, count_max
   !
   !Source_Crdnts= (/ 10000 , 16000 , 4000 /)    ! 1=North, 2=East, 3=vertical(plumbline)
   If(FirstPass) Write(2,"(A,F12.6,A,I11,A,3F10.1,A)") 'Start time for this chunk is set at ',StartT_sam(i_chunk)*sample_ms &
      ,' [ms] =',StartT_sam(i_chunk),'Samples, SourceGuess=',SourceGuess,';  1=North, 2=East, 3=vertical(plumbline@core)'
   !Write(2,*) 'SourceGuess=',SourceGuess,';  1=North, 2=East, 3=vertical(plumbline)'
   write(*,"(A,i3,A,f8.2,A)") achar(27)//'[33m',i_chunk,' @',StartT_sam(i_chunk)*sample_ms,'[ms]'//achar(27)//'[0m'  ! [1000D    !  //achar(27)//'[0m.'
   ! In main program, after option selection:
   Inquire(unit=14, opened=file14open)
   Inquire(unit=15, opened=file15open)
   !write(2,*) '!Ant-Read, Ant_nrMax:', LAnt_nrMax, HAnt_nrMax, Ant_nrMax
   HBA=.false.
   If(HAnt_nrMax.gt.0) HBA=.true.
   LBA=.false.
   If(LAnt_nrMax.gt.0) LBA=.true.
   If((Simulation.eq."")) WriteSimulation(2)=-1
   !write(2,*) 'Simulation,WriteSimulation(2):', Simulation,';',WriteSimulation(:), 'file14open:',file14open
   If((Simulation.eq."") .or. (WriteSimulation(2).gt.0)) Then  ! imaging on real data (="") or 'SelectData' (.gt.0) options
      If(LBA) Then
         If(file14open) then
            Rewind(unit=14)
            Rewind(unit=12)
         Else
            Open(unit=12,STATUS='unknown',ACTION='read',FORM ="unformatted", FILE = 'Book/RFI_Filters-v18.uft')
            Open(unit=14,STATUS='old',ACTION='read', FILE = 'Book/LOFAR_H5files_Structure-v18.dat')
         EndIf
      EndIf
      If(HBA) Then
         If(file15open) then
            Rewind(unit=15)
            Rewind(unit=13)
         Else
            Open(unit=13,STATUS='unknown',ACTION='read',FORM ="unformatted", FILE = 'Book/HRFI_Filters-v18.uft')
            Open(unit=15,STATUS='old',ACTION='read', FILE = 'Book/HLOFAR_H5files_Structure-v18.dat')
         EndIf
         !write(2,*) '!File 14 opened:', 'Book/LOFAR_H5files_Structure-v18.dat'
         !flush(unit=2)
      EndIf
   Else  ! simulation ne '' and width <0; run on simulation data
      If(file14open) then
         write(2,*) '****Error in SimulationRead, file 14 already open'
         Stop ' SimulationRead'
      EndIf
      Open(unit=14,STATUS='old',ACTION='read', FILE = 'files/'//TRIM(Simulation)//'_Structure_LBA.dat')
      If(HBA) Then
         If(file15open) then
            write(2,*) '****Error in SimulationRead, file 15 already open'
            Stop ' SimulationRead'
         EndIf
         Open(unit=15,STATUS='old',ACTION='read', FILE = 'files/'//TRIM(Simulation)//'_Structure_HBA.dat')
      EndIf
      Call SimulationRead(SourceGuess)
      Powr_eo(:,:)=1.d0  ! only matters that Powr_eo(0) and Powr_eo(1) are equal
      goto 9
   Endif
   If((Simulation.ne."") .and. (WriteSimulation(2).gt.0)) Then  ! 'SelectData' option
      Call CreateNewFolder(Simulation) ! create new folder when needed
      Open(Unit=30,STATUS='unknown',ACTION='write', FILE = 'files/'//TRIM(Simulation)//'_Structure_LBA.dat')
      write(2,*) 'created simulation file:','files/'//TRIM(Simulation)//'_Structure_LBA.dat'
      write(*,*) 'created simulation file:','files/'//TRIM(Simulation)//'_Structure_LBA.dat'
      If(HBA) Then
         Open(Unit=32,STATUS='unknown',ACTION='write', FILE = 'files/'//TRIM(Simulation)//'_Structure_HBA.dat')
         write(2,*) 'created simulation file:','files/'//TRIM(Simulation)//'_Structure_HBA.dat'
         write(*,*) 'created simulation file:','files/'//TRIM(Simulation)//'_Structure_HBA.dat'
      EndIf
   EndIf
   Ant_nr(i_chunk)=0
   LBA_nr(i_chunk)=0
   HBA_nr(i_chunk)=0
   AntNr_lw=1
   Powr_eo(:,:)=0.  ; NAnt_eo(:,:)=0.
   test=.false.
   Do i_type=0,1
      ir_file=0 ; ir_grp=0 ; ir_dst=0
      If(i_type.eq.0) Then
         If(.not. LBA) cycle
         U_bk=14
         U_flt=12
      ElseIf(i_type.eq.1) Then
         If(.not. HBA) cycle
         U_bk=15
         U_flt=13
      Endif
      Do i_file=1,130 ! 10 !80
      !i_file=1 ; rewind(unit=14)
         read(U_bk,*,IOSTAT=nxx)  Group_nr, filename
         If(nxx.ne.0) exit
         If(FirstPass .and. (i_file.eq.1))  write(2,*) 'First data file name=',filename
         !If(Production) write(2,*) 'file #=',i_file,', name=',filename
         !write(2,*) 'file #=',i_file,', name=',filename
         If(.not. Production) write(*,"(A,i3)", ADVANCE='NO') achar(27)//'[100D'//'file#=',i_file  ! [1000D    !  //achar(27)//'[0m.'
         Do i_grp=1,Group_nr
             read(U_bk,*) DSet_nr, Group_Names(i_grp)
             If(Diagnostics) Write(2,"(A,i0,A,A)") 'Group#',i_grp,' : ',trim(Group_names(i_grp))
             Do i_dst=1,DSet_nr
                 read(U_bk,*) DSet_Names(i_dst), STATION_ID, Ant_ID
                 !write(*,*) DSet_nr, Group_Names(i_grp), DSet_Names(i_dst), STATION_ID, Ant_ID
                 !If(STATION_ID .eq. 150) cycle  !  throuw out RS310  ****************************
                 If( COUNT(ExcludedStatID .eq. STATION_ID) .gt. 0) cycle
                 STATION_ID=STATION_ID*10 + i_type
                 StAntID= STATION_ID * 100 + Ant_ID
                 !write(2,*) DSet_nr, Group_Names(i_grp), DSet_Names(i_dst), STATION_ID, Ant_ID
                 i_eo=Mod(Ant_ID,2)
                 If( COUNT((PolFlp_SAI(1:PolFlp_nr)+i_eo) .eq. StAntID) .gt. 0) then
                     If(Diagnostics) &
                        write(2,*) '****Pole-Flip antenna:',STATION_ID,Ant_ID
                     ! cycle
                     Ant_ID=Ant_ID+1-2*i_eo
                     StAntID= STATION_ID * 100 + Ant_ID
                 endif
                 !
                 If( COUNT(BadAnt_SAI(1:BadAnt_nr) .eq. StAntID) .gt. 0) then
                     If(Diagnostics) write(2,*) '*****Bad antenna:',STATION_ID,Ant_ID, StAntID
                     cycle
                 endif
                 !
                 Do I_ant=1,Ant_nr(i_chunk)
                     If( (Ant_IDs(I_ant,i_chunk) .eq. Ant_ID) .and. (Ant_Stations(I_ant,i_chunk) .eq. STATION_ID) ) then
                         If(Diagnostics) write(2,*) '*****Double Data',STATION_ID,Ant_ID
                         goto 1  ! get next DSet
                     endif
                 enddo
                 !
                 Do while((i_file .gt. ir_file) .or. (i_grp .gt. ir_grp) .or. (i_dst .gt. ir_dst) )
                     !Write(2,*) 'reading Filter data',i_file,i_grp,i_dst,ir_file,ir_grp,ir_dst
                     Read(U_flt) ir_file,ir_grp,ir_dst, LFRAnt_crdnts, powr,&
                     DATA_LENGTH, SAMPLE_NUMBER_first, Absolute_TIME, DIPOLE_CALIBRATION_DELAY,nu_fltr
                 enddo
                 !write(2,*) '!', StAntID,'LFRAnt_crdnts, powr',LFRAnt_crdnts, powr
                 if(powr.lt.1.d-6) cycle
                 If((i_file .ne. ir_file) .or. (i_grp .ne. ir_grp) .or. (i_dst .ne. ir_dst) ) then
                     Write(2,*) '****Filter data not found',i_file,i_grp,i_dst, ir_file,ir_grp,ir_dst
                     cycle   ! get next DSet
                 Endif
                 !
                 Call RelDist(SourceGuess,LFRAnt_crdnts,RDist)
                 Ant_RawSourceDist(Ant_nr(i_chunk)+1,i_chunk)=RDist         ! units of samples
                 If(WriteSimulation(2).gt.0) Then
                     !write(2,*) 'RDist',RDist,AntNr_lw,Ant_nr(i_chunk)+1
                     RDist=Ant_RawSourceDist(AntNr_lw,i_chunk)
                 EndIf
                 !RDist=0.
                 !Get LOFAR Station calibrations = negative of delays
                 Call Station_ID2Calib(STATION_ID,Ant_ID,StatAnt_Calib, Calibrated) ! StatAnt_Calib in units of samples
                 If(CalibratedOnly .and. (.not. Calibrated) .and. Diagnostics) Then
                  write(2,*) 'STATION_ID,Ant_ID was not calibrated',STATION_ID,Ant_ID,StatAnt_Calib
                  cycle
                 EndIf
                 !Write(2,*) '!',Calibrated,'calibration:', STATION_ID,Ant_ID,StatAnt_Calib, DIPOLE_CALIBRATION_DELAY/Sample
                 !
                 Powr_eo(i_eo,i_type)=Powr_eo(i_eo,i_type)+Powr
                 NAnt_eo(i_eo,i_type)=NAnt_eo(i_eo,i_type)+1
                 !
                 !Absolute_TIME  ! should be the same for all
                 !SAMPLE_NUMBER_first  ! Number of samples, after 'Absolute_TIME' for the first recording
                 !IF(Ant_ID.ge.90 .and. Ant_ID.le.93)
                 !write(2,*) '!StatAnt_Calibratio',StAntID, StAntID, &
                 !          DIPOLE_CALIBRATION_DELAY/Sample, StatAnt_Calib, DIPOLE_CALIBRATION_DELAY/Sample + StatAnt_Calib
                 !If(i_chunk.eq.1) Write(2,*) 'AntennaRead, timecorrection',StAntID, &
                 !    DIPOLE_CALIBRATION_DELAY/Sample, StatAnt_Calib, DIPOLE_CALIBRATION_DELAY/Sample+StatAnt_Calib
                 StatAnt_Calib=DIPOLE_CALIBRATION_DELAY/Sample + StatAnt_Calib  ! in units of samples
                 T_Offset=RDist + StatAnt_Calib  ! in units of samples
                 Sample_Offset = INT(T_Offset) ! in units of sample size
                 SubSample_Offset = T_Offset - Sample_Offset ! in units of sample size
                 Dset_offset=StartT_sam(i_chunk) + Sample_Offset - SAMPLE_NUMBER_first
                 !write(2,*) 'StartT_sam(i_chunk),:',StartT_sam(i_chunk), Sample_Offset, SAMPLE_NUMBER_first
                 !If(Ant_nr(1).eq.0) write(2,*) ', TimeOffset', T_Offset, Dset_Offset,DATA_LENGTH, Sample_Offset - SAMPLE_NUMBER_first
                 If(Dset_offset .gt.DATA_LENGTH) then
                     If(Diagnostics) write(2,*) '****DATA_LENGTH=',Dset_offset, ' greater than ',DATA_LENGTH
                     cycle       ! get next DSet
                 endif
                 If(Dset_offset .lt.0) then
                     If(Diagnostics) write(2,*) '****Dset_offset=',Dset_offset, ' less than zero'
                     cycle       ! get next DSet
                 endif
                 !
                 ! Retrieve a chunk of data
                 DataReadError=0
                 Call GetDataChunk(Group_Names(i_grp),DSet_Names(i_dst), Chunk, Dset_offset, Time_dim, DataReadErr=DataReadErr)
                 If(DataReadErr.ne.0) then
                  write(2,*) '************ error in data-read',DataReadErr,', file #=',i_file,', name=',filename
                  DataReadError=-1
                  write(*,*) 'error in data-read'
                  cycle
                  !return
                 endif
                 !Call ListDataAtt(Group_Names(i_grp),DSet_Names(i_dst), Chunk, Dset_offset, Time_dim)
                 ! write(2,*) '!Ant-Read, start zeros:', Ant_nr(i_chunk)
                 !
                  NZero=COUNT(Chunk.eq.0)
                  !Ntwo=COUNT(Chunk.eq.2)
                  N_one=COUNT(Chunk.eq.1)
                  If(N_one*1.0/NZero .lt. .5) then
                     If(Diagnostics) Write(2,"(I6,A,F7.4,A,2i8)") StAntID,' ****Ratio of twos/zeros in this chunk=',N_one*1.0/NZero
                     ! This Dset has too many zeros.
                     cycle       ! get next DSet
                 endif
                 !
                 If( RunMode.eq.7 .or. RunMode.eq.24 .or. Runmode.eq.2) Then
                     N_one=COUNT( Chunk(:).gt.2045 .or. Chunk(:).lt.-2045)  !! saturates at 2047
                     !NZero=COUNT( Chunk(:).lt.-2020)  !! saturates at 2047
                     If(N_one.gt.SaturatedSamplesMax) Then
                        write(2,"(A,I4,', chunk#',I4,A,I4,', @')",ADVANCE='NO') &
                              Statn_ID2Mnem(STATION_ID), Ant_ID, i_chunk, ' excluded, #Saturated Samples=',N_one
                        If(N_one.lt.10) then
                           Do i_time=1,Time_dim
                              If(Chunk(i_time).gt.2020 .or. Chunk(i_time).lt.-2020) Then
                                 Write(2,"(1x,I6,':',I5)",ADVANCE='NO') i_time, Chunk(i_time)
                              EndIf
                           EndDo
                           Write(2,*) ' '
                        Else
                           Write(2,*) ' too many to list!'
                        EndIf
                        cycle
                     Else
                        If(N_one .ge. SaturatedSamplesMax/2) Then
                           write(2,"(A,I4,', chunk#',I4,A,I3,', @')",ADVANCE='NO') &
                                 Statn_ID2Mnem(STATION_ID), Ant_ID, i_chunk, ' used with #Saturated Samples=',N_one
                           If(N_one.lt.10) then
                              Do i_time=1,Time_dim
                                 If(Chunk(i_time).gt.2020 .or. Chunk(i_time).lt.-2020) Then
                                    Write(2,"(1x,I6,':',I5)",ADVANCE='NO') i_time, Chunk(i_time)
                                 EndIf
                              EndDo
                              Write(2,*) ' '
                           Else
                              Write(2,*) ' too many to list!'
                           EndIf
                        EndIf
                     EndIf
                 EndIf
                 If(Diagnostics) Then
                    ChMax=MaxVal(Chunk)
                    ChMin=MinVal(Chunk)
                    write(2,*) StAntID,' min & max=', ChMin, ChMax,', RDist=',RDist, &
                        ', Calibr=', - DIPOLE_CALIBRATION_DELAY/Sample + StatAnt_Calib, &
                        ', SampNrFrst=',SAMPLE_NUMBER_first,', norm=',sqrt(Powr)
                 EndIf
                 !
                 Ant_nr(i_chunk)=Ant_nr(i_chunk)+1
                 If(Ant_nr(i_chunk) .gt. (Ant_nrMax+HAnt_nrMax) ) then
                     Write(2,*) '****Max nr of antennas exceeded@ file,grp,dst',i_file,i_grp,i_dst,'Ant_nr',Ant_nr(i_chunk)
                     stop 'Antenna numbers exceeded'
                 Endif
                 !
                 ! write(2,*) '!Ant-Read, Ant_nrMax:', Ant_nr(i_chunk)
                 Ant_IDs(Ant_nr(i_chunk),i_chunk)=Ant_ID
                 Ant_Stations(Ant_nr(i_chunk),i_chunk)=STATION_ID
                 Ant_pos(:,Ant_nr(i_chunk),i_chunk)=LFRAnt_crdnts(:)
                 !write(2,*) Ant_nr(i_chunk), ' ,position: ', LFRAnt_crdnts(:)
                 !
                 sgn=+1
                 If(any(SignFlp_SAI(1:SgnFlp_nr) .eq. StAntID,1) ) sgn=-1 ! Take care of sign flips in I*2
                 If(WriteSimulation(2).le.0) Then
                    If(RunMode.eq.1) then
                     RTime_s(:)=Chunk(:)*Hann(:)*sgn
                    Else
                     RTime_s(:)=Chunk(:)*Hann(:)*sgn*NormEvenOdd/sqrt(Powr) !  sqrt(power) level is normalized to 100.
      If(STATION_ID.eq.-61) Then
         X=0.
         Do i=1,Time_dim
            X=X+ 1.0*Chunk(i)*Chunk(i)! SUM(Chunk(:)*Chunk(:)), overflows, and thus nonsense
         EndDo
         X=2*X*(sgn*NormEvenOdd/sqrt(Powr))**2/100.
         Av=MaxVal(Chunk) ! RTime_s)
         Bv=Av-MinVal(Chunk)
         Av=Av-Bv/2.
         nu_i=27*(Time_dim/2)/100
         nu_f=91*(Time_dim/2)/100
         Call RFTransform_CF(RTime_s,Cnu_s)
         !write(2,"(/A,2I6,G14.5,10(A,G14.5))") '!AntennaRead, power:',STATION_ID,Ant_ID,100.*X, &
         write(2,"(A,2I4,G11.3,5(A,F6.1),i2,5(A,F6.1))") '!AntennaRead, power:',STATION_ID,Ant_ID,100.*X, & ! ', *Hann=',2*SUM(RTime_s(:)*RTime_s(:))/X,
            ', no-filt(nu):', SUM(abs(Cnu_s(0:Time_dim/2))**2)*Time_dim/X, &
            '% Freq_window:',SUM(abs(Cnu_s(nu_i:nu_f))**2)*Time_dim/X, '% Power(nu=0)=', abs(Cnu_s(0))**2*Time_dim/X, &
            '% bias(t)=', Av, ', range(t)=', Bv, i_eo,', norm=', sqrt(Powr) ! , ', Ampl(nu=0)=', Real(Cnu_s(0))
      EndIf
                    Endif
                 Else
                    RTime_s(WriteSimulation(1):WriteSimulation(1)+WriteSimulation(2))= &
                        Chunk(WriteSimulation(1):WriteSimulation(1)+WriteSimulation(2))*sgn*NormEvenOdd/sqrt(Powr) ! /Time_dim)
                 EndIf
                 !
                 !SubSample_Offset=0.
                 !nu_fltr(:)=1.
                 Call RFTransform_CF_Filt(RTime_s,nu_fltr,SubSample_Offset,Cnu_s)
                 !
      If(STATION_ID.eq.-61) Then
         write(2,"(A,8G14.5)") '!AntennaRead, filt(nu):', SUM(abs(Cnu_s(0:Time_dim/2))**2)*Time_dim/X
         !write(2,"(A,8G14.5)") '!AntennaRead, nfilt(nu):', SUM(abs(Cnu_s(nu_i:nu_f))**2)*Time_dim
      EndIf
                 !
                 If(i_type.eq.0) Then
                    LBA_nr(i_chunk)=LBA_nr(i_chunk)+1
                    Call RFTransform_CF2CT(Cnu_s,CTime_spectr(1,LBA_nr(i_chunk),i_chunk) )
                 !    Av=MaxVal(RTime_s)
                 !    Bv=MinVal(RTime_s)
                 !    write(2,*) 'Normalized: min & max=', Bv, Av
                 ! write(2,*) 'STATION_ID,Ant_ID',STATION_ID,Ant_ID,StatAnt_Calib
                 ElseIf(i_type.eq.1) Then
                    HBA_nr(i_chunk)=HBA_nr(i_chunk)+1
                    !SubSample_Offset=0.5
                    !SubSample_Offset=0.  !  SubSample_Offset may not work properly  here for HBA
        !             Call RFTransform_CF_Filt(RTime_s,nu_fltr,SubSample_Offset,Cnu_s)
                    ! spectrum is shifted in the negative time direction by SubSample_Offset
                    ! =============================== HBA circus
                    ! cnu_dim from 0 -- timedim/2;  Hcnu_dim from 0 -- time_dim
                    !  test
                    !SubSample_Offset=1.
                    !test=.true.
               !     Call RFTransform_CF(RTime_s,Cnu_s)
                !    Cnu_s(0:cnu_dim)=Cnu_s(0:cnu_dim)*nu_fltr(0:cnu_dim)
                    !------end test
                    !
                    HCnu_s(0:cnu_dim)=0.d0  ! zero all frequencies below 100MHz
                    Do i_time=0,cnu_dim
                        HCnu_s(cnu_dim+i_time)=conjg(Cnu_s(cnu_dim-i_time))
                    EndDo
                    !
                    Call DAssignFFT()                      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    Call RFTransform_su(2*Time_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                    Call RFTransform_CF2CT(HCnu_s,CTime_Hspectr(1,HBA_nr(i_chunk),i_chunk) )
                    Call DAssignFFT()                      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    Call RFTransform_su(Time_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!      If(STATION_ID.eq.1211) write(2,"(A,2I6,3G14.5)") '!AntennaRead, shifted(t):',STATION_ID,Ant_ID,  &
!            SUM(abs(CTime_Hspectr(:,HBA_nr(i_chunk),i_chunk))**2)/2.
                    !
                    If(test) then  !   Just for testing: --------------------------
                        write(Ant_name,"(I6.6)") StAntID
                        Open(unit=20, STATUS='unknown',ACTION='write', FILE = 'files/Testspectr'//TRIM(Ant_name)//'.dat')
                        Call RFTransform_CF2RT(Cnu_s,RTime_s )
                        write(20,*) SubSample_Offset, Cnu_s(0), Cnu_s(cnu_dim)
                        Do i_time=1,Time_dim
                           write(20,*) i_time, RTime_s(i_time), Real(CTime_Hspectr(2*i_time-1,HBA_nr(i_chunk),i_chunk)), &
                                 Real(CTime_Hspectr(2*i_time,HBA_nr(i_chunk),i_chunk))
                        EndDo
                        Close(unit=20)
                    ! shows that RTime_s(i_time) and  Real(CTime_Hspectr(2*i_time-1,HAnt_nr(i_chunk),i_chunk)) are equal
                    EndIf  ! ---------------------------------------------
                    ! =================================================
                 EndIf
                 !
             1   continue
             EndDo  ! i_dst=1,DSet_nr
             ! re-order neighboring antennas according to increasing antenna ID (assumed in interferometric imaging)
             ! re-order : CTime_spectr(1,i_ant,i_chunk) ; Ant_IDs(i_ant,i_chunk) ; Ant_pos(:,i_ant,i_chunk), only for this STATION_ID
             Do i_ant=2,Ant_nr(i_chunk)
               If(Ant_Stations(i_ant-1,i_chunk).eq.Ant_Stations(i_ant,i_chunk)) Then ! apply to antennas from a single station only
                  If(Ant_IDs(i_ant-1,i_chunk).gt.Ant_IDs(i_ant,i_chunk)) Then  ! flip order
                     !write(2,*) '!change order ',i_ant-1,Ant_IDs(i_ant-1,i_chunk),' and ',i_ant,Ant_IDs(i_ant,i_chunk), &
                     !   ' for station',Ant_Stations(i_ant,i_chunk)
                     If(i_type.eq.0) Then
                        CTime_s(1:Time_dim)=CTime_spectr(1:Time_dim,i_ant,i_chunk)
                        CTime_spectr(1:Time_dim,i_ant,i_chunk)=CTime_spectr(1:Time_dim,i_ant-1,i_chunk)
                        CTime_spectr(1:Time_dim,i_ant-1,i_chunk)=CTime_s(1:Time_dim)
                     ElseIf(i_type.eq.1) Then
                        j=i_ant-LBA_nr(i_chunk)
                        CTime_s(1:2*Time_dim)=CTime_Hspectr(1:2*Time_dim, j,i_chunk)
                        CTime_Hspectr(1:2*Time_dim, j,i_chunk)=CTime_Hspectr(1:2*Time_dim, j-1,i_chunk)
                        CTime_Hspectr(1:2*Time_dim, j-1,i_chunk)=CTime_s(1:2*Time_dim)
                     EndIf
                     Ant_ID=Ant_IDs(i_ant,i_chunk)
                     LFRAnt_crdnts(:)=Ant_pos(:,i_ant,i_chunk)
                     RDist=Ant_RawSourceDist(i_ant,i_chunk)
                     !
                     Ant_IDs(i_ant,i_chunk)=Ant_IDs(i_ant-1,i_chunk)
                     Ant_pos(:,i_ant,i_chunk)=Ant_pos(:,i_ant-1,i_chunk)
                     Ant_RawSourceDist(i_ant,i_chunk)=Ant_RawSourceDist(i_ant-1,i_chunk)
                     !
                     Ant_IDs(i_ant-1,i_chunk)=Ant_ID
                     Ant_pos(:,i_ant-1,i_chunk)=LFRAnt_crdnts(:)
                     Ant_RawSourceDist(i_ant-1,i_chunk)=RDist
                  EndIf
               EndIf
             Enddo
             If((Simulation.ne."") .and. (WriteSimulation(2).gt.0)) Then !  'SelectData' option
               AntNr_up=Ant_nr(i_chunk)
               If((AntNr_up-AntNr_lw+1) .gt.0) Then
                  write(30,*) STATION_ID,TRIM(Statn_ID2Mnem(STATION_ID)),AntNr_up-AntNr_lw+1, WriteSimulation(2)
                  !Call CreateNewFolder(Simulation) ! create new folder when needed
                  OPEN( Unit=31, STATUS='unknown',ACTION='write', Form='unformatted', &
                     FILE = 'files/'//TRIM(Simulation)//'_'//TRIM(Statn_ID2Mnem(STATION_ID))//'.udt')
                  ! Structure of this unformatted file
                  ! rec 1: StatStartTime [ms],  NtSamples, Ant_nr [=number of different antennas (2 for each dipole)]
                  !  followed by 'j_ant=1:ant_nr' data for each antenna
                  ! rec 2..N: ('Dual' or 'Even' or 'Odd'), Ant_IDs(j_ant), Ant_pos(:,j_ant) ! info for dipole 'j_ant'
                  ! rec N+j_ant: time trace antenna 'j_ant' (real*4)
                  ! these records are repeated.
                  ! at present all dipoles are necessarily 'Dual'
                  !
                  !Open(Unit=31,STATUS='unknown',ACTION='write', &
                  !      FILE = 'files/'//TRIM(Simulation)//'_'//TRIM(Statn_ID2Mnem(STATION_ID))//'.dat')
                  !Write(31,*) 'StartTime_ms= ',(StartT_sam(i_chunk)+WriteSimulation(1)+RDist)*sample*1000.d0, &
                  !      ' N_samp= ', WriteSimulation(2), ' noise_power= 1.'
                  Write(31) (StartT_sam(i_chunk)+WriteSimulation(1)+RDist)*sample*1000.d0, &
                        WriteSimulation(2), AntNr_up-AntNr_lw+1
                  Dubbel=.false.
                  Do i_ant=AntNr_lw,AntNr_up
                     i_eo=mod(Ant_IDs(i_ant,i_chunk),2)
                     If((i_ant.eq.AntNr_up) .and. (i_eo.eq.0)) then
                        Write(31) 'Even ', 'NEh= ',Ant_pos(:,i_ant,i_chunk)  !  Ant_IDs(i_ant,i_chunk),
                        cycle
                     ElseIf(i_eo.eq.0) then
                        If(Ant_IDs(i_ant+1,i_chunk).eq.(Ant_IDs(i_ant,i_chunk)+1)) Then
                           Dubbel=.true.
                        Else
                           Write(31) 'Even', Ant_IDs(i_ant,i_chunk),Ant_pos(:,i_ant,i_chunk)  !  Ant_IDs(i_ant,i_chunk),
                           cycle
                        EndIf
                     Else
                        If(Dubbel) Then
                           Write(31) 'Dual', Ant_IDs(i_ant,i_chunk)-1, Ant_pos(:,i_ant,i_chunk)
                           Dubbel=.false.
                        Else
                           Write(31) 'Odd ', Ant_IDs(i_ant,i_chunk), Ant_pos(:,i_ant,i_chunk)
                        EndIf
                     EndIf
                  EndDo
                  !Write(31,*) '===   0 0. 0. 0.'
                  !write(2,*) WriteSimulation(1),(WriteSimulation(1)+WriteSimulation(2)),AntNr_lw,AntNr_up,i_chunk
                  !Do i_time=WriteSimulation(1),(WriteSimulation(1)+WriteSimulation(2))
                  !   write(31,*) Real( CTime_spectr(i_time,AntNr_lw:AntNr_up,i_chunk),4 )
                  Do i_ant=AntNr_lw,AntNr_up
                     write(31) Real( CTime_spectr(WriteSimulation(1):(WriteSimulation(1)+WriteSimulation(2)),i_ant,i_chunk),4 )
                  EndDo
                  Close(Unit=31)
                  AntNr_lw=AntNr_up+1
               EndIf
             EndIf
             !
         Enddo ! i_grp=1,Group_nr
   !   !
   !   Close(unit=22)
   !   !
         !
      Enddo ! i_file=1,....
      !write(2,*) 'End i_type=', i_type
   Enddo ! i_type=1,....
!   Close(unit=24)
   !
   If((Simulation.ne."") .and. (WriteSimulation(2).gt.0)) Then
      Close(Unit=30)
      write(2,*) 'true data written to files to be used in sumulation runs'
      Stop ' Normal end Simulation or SelectData'
   EndIf
   !
   If(Ant_nr(i_chunk).le.6) then
      write(2,*) Ant_nr(i_chunk),' are too few active LBA-antennas for chunk starting at ',StartT_sam(i_chunk)
      write(*,*) Ant_nr(i_chunk),' are too few active LBA-antennas for chunk starting at ',StartT_sam(i_chunk)
      AntennaNrError=AntennaNrError+1
      If(AntennaNrError.gt.10) then
         write(2,*) 'AntennaNrError:',AntennaNrError
         stop 'Too few antennas fo a while'
      Endif
   Else
      AntennaNrError=0
   Endif
   !
   write(*,*) ' '
   Powr_eo(:,:)=Powr_eo(:,:)/NAnt_eo(:,:)
9  Continue  ! entry after simulation read
   NormEven=sqrt(2.*Powr_eo(0,0)/(Powr_eo(0,0)+Powr_eo(1,0)))/NormEvenOdd  ! to undo the factor 100 (and more) that was introduced in antenna-read
   NormOdd=sqrt(2.*Powr_eo(1,0)/(Powr_eo(0,0)+Powr_eo(1,0)))/NormEvenOdd  ! to undo the factor that was introduced in antenna-read
   If(HBA) Then
      NormHEven=sqrt(2.*Powr_eo(0,1)/(Powr_eo(0,1)+Powr_eo(1,1)))/NormEvenOdd  ! to undo the factor 100 (and more) that was introduced in antenna-read
      NormHOdd=sqrt(2.*Powr_eo(1,1)/(Powr_eo(0,1)+Powr_eo(1,1)))/NormEvenOdd  ! to undo the factor that was introduced in antenna-read
   EndIf
   If(FirstPass) Then
      write(2,*) 'ratio of average background amplitude for LBA even/odd antennas',NormEven/NormOdd,&
         (Powr_eo(0,0)+Powr_eo(1,0))/2.
      If(HBA) write(2,*) 'ratio of average background amplitude for HBA even/odd antennas',NormHEven/NormHOdd, &
         (Powr_eo(0,1)+Powr_eo(1,1))/2
      !write(2,*) sum(ABS(CTime_spectr(:,10,1))**2),Time_dim,Powr
      write(2,*) 'Height-Corrected Index of Refraction=', HeightCorrectIndxRef
      !write(2,*) 'NAnt_eo:',NAnt_eo(:)
      FirstPass=.false.
   EndIf
   !
   ! time recording
   !call cpu_time(CPUstopTime)
   !CALL SYSTEM_CLOCK(WallCount, powr)  !  Wallcount_rate)
   Call CPU_usage(Message='finish reading' )
   !WRITE(*,*) 'Wall clock time:',WallCount, ', countrate:', powr  ! count_rate, count_max
   !WallstopTime=WallCount/powr
   !WRITE(*,"(A,F9.3,A)") 'Wall clock time:',WallstartTime, '[s]'  ! count_rate, count_max
   !TotCPURead=TotCPURead + CPUstopTime-CPUstartTime
   !TotWallRead=TotWallRead + WallstopTime-WallstartTime
   !WRITE(*,"(A,F9.3,F12.6,A)") 'Totals Read, Wall & CPU:',TotCPURead, TotWallRead, '[s]'  ! count_rate, count_max
   !WRITE(*,"(A,F9.3,F12.6,A)") 'Totals Fit, Wall & CPU:',TotCPUFit, TotWallFit, '[s]'  ! count_rate, count_max
   !WRITE(*,"(A,F10.3,F13.6,A,F10.3,F13.6,A)") 'TotWall&CPU; Read:', TotWallRead,TotCPURead,'[s], Fit:', TotWallFit,TotCPUFit,'[s]'
   !
   !write(2,*) '!QuickCurtain:', QuickCurtain
   If(QuickCurtain) Then
      i_sample=Time_dim/2
      Location(1:3)=SourceGuess(1:3)
      write(2,*) 'QuickCurtain@loc', Location(1:3)
      Flush(unit=2)
      !Call Alloc_ThisSource
      Call Alloc_FitParams
      write(2,*) 'QuickCurtain@i_sample', i_sample
      Flush(unit=2)
      Call Find_unique_StatAnt
      Flush(unit=2)
      Call GetRefAnt(ChunkNr_dim)
      Flush(unit=2)
      Call PlotsingleCurtainSpectrum(i_sample, Location)
      Call GLEplotControl(Submit=.true.)
      stop 'QuickCurtain'
   EndIf
   !
   Return
   !
End Subroutine AntennaRead
! ========================
Subroutine AntennaCount(N_LAnt, N_LStat, N_HAnt, N_HStat)
   use Chunk_AntInfo, only : ExcludedStatID, BadAnt_nr, BadAnt_SAI
   use Chunk_AntInfo, only : Simulation, WriteSimulation
   Implicit none
   Integer, intent(out) :: N_LAnt, N_LStat, N_HAnt, N_HStat
   Integer :: i_file,i_grp,i_dst, nxx
   Integer :: Group_nr, DSet_nr
   Integer :: Ant_ID, STATION_ID, Nant
   Integer :: Used_STATION_ID(0:50)
   Character(len=50) ::  txt
   Logical :: H_exists=.false.
   Logical :: L_exists
   Logical :: RdSimul=.false.
   Real :: X
   If((Simulation.eq."") .or. (WriteSimulation(2).gt.0)) Then  ! imaging on real data (="") or 'SelectData' (.gt.0) options
      INQUIRE(FILE = 'Book/LOFAR_H5files_Structure-v18.dat', exist=L_exists)  ! in the main Fsash directory
      If(L_exists) Open(unit=14,STATUS='old',ACTION='read', FILE = 'Book/LOFAR_H5files_Structure-v18.dat')
      INQUIRE(FILE = 'Book/HLOFAR_H5files_Structure-v18.dat', exist=H_exists)  ! in the main Fsash directory
      If(H_exists) Open(unit=15,STATUS='old',ACTION='read', FILE = 'Book/HLOFAR_H5files_Structure-v18.dat')
   Else  ! simulation ne '' and width <0; run on simulation data
      RdSimul=.true.
      INQUIRE(FILE = 'files/'//TRIM(Simulation)//'_Structure_LBA.dat', exist=L_exists)  ! in the main Fsash directory
      If(L_exists) Open(unit=14,STATUS='old',ACTION='read', FILE = 'files/'//TRIM(Simulation)//'_Structure_LBA.dat')
      !INQUIRE(FILE = 'files/'//TRIM(Simulation)//'_Structure_HBA.dat', exist=H_exists)  ! in the main Fsash directory
      !If(H_exists) Open(unit=15,STATUS='old',ACTION='read', FILE = 'files/'//TRIM(Simulation)//'_Structure_HBA.dat')
   Endif
   N_LAnt=0
   N_LStat=0
   Used_STATION_ID(0)=0
   If(RdSimul) Then
   ! Simulation is implemented for LBA only!!!
      Do i_file=1,130
         read(14,*,IOSTAT=nxx) STATION_ID, txt,  Nant ! , LOFAR_name(i_file)
         If(nxx.ne.0) exit
         !LOFAR_number(i_file)=STATION_ID
         !inquire(FILE = 'files/'//TRIM(Simulation)//'_'//LOFAR_name(i_file)//'.udt', exist=UnformattedSimul)
         If( COUNT(ExcludedStatID .eq. STATION_ID) .gt. 0) cycle
         N_LAnt=N_LAnt + Nant
         If( COUNT(Used_STATION_ID(0:N_LStat) .eq. STATION_ID) .gt. 0) cycle
         N_LStat=N_LStat+1
         Used_STATION_ID(N_LStat) = STATION_ID
      EndDo
      N_HAnt=0
      N_HStat=0
      Close(Unit=14)
   Else
      If(L_exists) Then
         Do i_file=1,130 ! 10 !80
            read(14,*,IOSTAT=nxx)  Group_nr, txt
            If(nxx.ne.0) exit
            Do i_grp=1,Group_nr
                read(14,*,IOSTAT=nxx) DSet_nr, txt
                If(nxx.ne.0) exit
                Do i_dst=1,DSet_nr
                    read(14,*) txt, STATION_ID, Ant_ID
                    !write(*,*) DSet_nr, Group_Names(i_grp), DSet_Names(i_dst), STATION_ID, Ant_ID
                    !If(STATION_ID .eq. 150) cycle  !  throuw out RS310  ****************************
                    If( COUNT(ExcludedStatID .eq. STATION_ID) .gt. 0) cycle
                    N_LAnt=N_LAnt+1
                    If( COUNT(Used_STATION_ID(0:N_LStat) .eq. STATION_ID) .gt. 0) cycle
                    N_LStat=N_LStat+1
                    Used_STATION_ID(N_LStat) = STATION_ID
                EndDo
            EndDo
         EndDo
         Close(Unit=14)
      EndIf
      !
      N_HAnt=0
      N_HStat=0
      If(H_exists) Then
         Used_STATION_ID(0)=0
         Do i_file=1,130 ! 10 !80
            read(15,*,IOSTAT=nxx)  Group_nr, txt
            If(nxx.ne.0) exit
            Do i_grp=1,Group_nr
                read(15,*) DSet_nr, txt
                Do i_dst=1,DSet_nr
                    read(15,*) txt, STATION_ID, Ant_ID
                    !write(2,*) DSet_nr, i_grp, i_dst, STATION_ID, Ant_ID, N_HStat, N_HAnt
                    !If(STATION_ID .eq. 150) cycle  !  throuw out RS310  ****************************
                    If( COUNT(ExcludedStatID .eq. STATION_ID) .gt. 0) cycle
                    N_HAnt=N_HAnt+1
                    If( COUNT(Used_STATION_ID(0:N_HStat) .eq. STATION_ID) .gt. 0) cycle
                    N_HStat=N_HStat+1
                    Used_STATION_ID(N_HStat) = STATION_ID
                EndDo
            EndDo
         EndDo
         Close(Unit=15)
      EndIf
   EndIf
   Return
End Subroutine AntennaCount
! ========================
! ========================
Subroutine PlotSpectra(i_chunk)
   use constants, only : dp,sample
   use DataConstants, only : Time_dim, DataFolder, Diagnostics ! , ChunkNr_dim
   use Chunk_AntInfo, only : Ant_nr, CTime_spectr
   Implicit none
   Integer, intent(in) :: i_chunk
   Integer :: unt
   !
   Real(dp) :: B ! ,Freq, Freq_s(0:Cnu_dim),
   !Integer :: j, ChMax, ChMin,i_file,i_grp,i_dst, ir_file,ir_grp,ir_dst
   !character(len=2) :: chj
   !Real(dp) :: Freq, Freq_s(0:Cnu_dim), B
   !
   Integer :: i_ant, i_time !, i_eo, unt, nxx, StAntID !, sgn
   !Real*8 :: SubSample_Offset, LFRAnt_crdnts(3), RDist, T_Offset,StatAnt_Calib !,StartTime_ms
   CHARACTER(LEN=40) :: GLE_file
   CHARACTER(LEN=6) :: txt
      write(*,*) 'MakePlots'
      write(2,*) 'number of antennas=',Ant_nr(i_chunk)
      Do I_ant=1,Ant_nr(i_chunk)
          write(txt,"(i3.3,'-'i2.2)") I_ant,i_chunk
          !write(2,*) txt, Ant_IDs(I_ant,i_chunk), Ant_Stations(I_ant,i_chunk)
          !
          B=maxval(real(CTime_spectr(:,i_Ant,i_chunk)))
          Open(Unit=9,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//'LOFAR_Time'//txt//'.dat')
          Do i_time=1,Time_dim
              Write(9,*) i_time, Real(CTime_spectr(i_time,i_Ant,i_chunk))/B, Imag(CTime_spectr(i_time,i_Ant,i_chunk))/B
          Enddo
          close(unit=9)
      Enddo
      !unt=9 ; GLE_file='GLE_file'  ! CurtainPlot
      unt=9 ; GLE_file='CurtainPlot'  !
      Call GLE_script(unt, GLE_file, Ant_nr(i_chunk),i_chunk)
      Open(UNIT=10,STATUS='unknown',ACTION='WRITE',FILE='GLE-plots.sh')
      !write(10,"('gle -d pdf GLE_file.gle')")
      Write(10,"(A)") 'gle -d pdf '//trim(GLE_file)//'.gle'
      !write(10,"('gle -d pdf CurtainPlot.gle')") !  Command to produce curtain plots
      Close(unit=10)
      !
      Return
End Subroutine PlotSpectra
! ===========================================================================
Subroutine GLE_script(unt, file, Ant_nr, i_chunk)
    use DataConstants, only : Time_dim, DataFolder
    use Chunk_AntInfo, only : CTime_spectr, Ant_Stations, Ant_IDs, Ant_pos, Ant_RawSourceDist
    use DataConstants, only : Used_StationNr
    use StationMnemonics, only : Station_ID2Mnem
    use ThisSource, only : PeakPos
    Implicit none
    Integer, intent(in) :: unt, Ant_nr, i_chunk
    Character(Len=30), intent(in) :: file
    integer :: I_ant, i_stat, Stat_nr, Stations(1:Used_StationNr), i_time, i_c
    integer :: A(1), Plot_scale,ch1,ch2,VLine
    Character(len=6) :: txt
    Character(len=6) :: Station_Mnem
    real*8 :: plot_offset,b
    !
    !    B=maxval(real(CTime_spectr(:,1))) ! only normalized spectra are written
    Plot_scale= 1. ! int(B) + 2.
    VLine=Time_dim/2
    ch1=100
    ch2= Time_dim-100
    If( (PeakPos(1).gt. ch1) .and. (PeakPos(1).lt. ch2) ) Then
      VLine=PeakPos(1)
      ch1=Vline-500
      ch2=VLine+500
      If( (PeakPos(2).lt. Time_dim/2) ) Then
         ch1=Vline-PeakPos(2)
         ch2=VLine+PeakPos(2)
      EndIf
    EndIf
    !write(2,*) 'B', B!,Ant_Stations
    !
    Open(UNIT=20,STATUS='unknown',ACTION='WRITE',FILE=trim(file)//'-peak_even.dat')
    Open(UNIT=21,STATUS='unknown',ACTION='WRITE',FILE=trim(file)//'-peak_odd.dat')
    Open(UNIT=UNT,STATUS='unknown',ACTION='WRITE',FILE=trim(file)//'.gle')
    Write(unt,"(A)") '! COMMAND:  gle -d pdf '//trim(file)//'.gle'
    Write(unt,"(A,3(/A),3(/A,I0))") 'size 63 82','set font pstr fontlwidth 0.08 hei 1.2 just CC',&
        'include "../DiscreteColor.gle"',&
        'set lwidth 0.1','t_min = ',ch1,'t_max = ',ch2,'scl = ',Plot_scale
!        'set lwidth 0.1','t_min = 13250 !12500 !0','t_max = 13350 !15000 !',Time_dim,'scl = ',Plot_scale
    Write(unt,"(A,12(/A))")  'amove 4 4','begin graph','  size 55 75','  vscale 1','  hscale 1',&
        '   title  "Time Spectra"','   xtitle "time [samples]"','   ytitle "Amplitude"',&
        '   xaxis min t_min max t_max ! nticks 10','   yaxis  min 0 max 2 dticks 5 dsubticks 1', &
        '   x2labels on',' end graph'
    Write(unt,"(A,i5.5,A,1(/A))")  'amove xg(',Vline,') 4','rline 0 75'
    !
    Stat_nr=0
    Stations(:)=0
    write(2,*) 'maximum between channels',ch1,ch2
    Do I_ant=1,Ant_nr
        B=maxval(abs(CTime_spectr(ch1:ch2,I_ant,i_chunk)))
        A=maxloc(abs(CTime_spectr(ch1:ch2,I_ant,i_chunk)))+CH1-1
        !If( (mod(Ant_IDs(I_ant,i_chunk),2).eq.0) ) then
        !    write(2,"(A,i7,F9.1,9x,I4,I4)") 'maxvl-e',A,B,Ant_Stations(I_ant,i_chunk),Ant_IDs(I_ant,i_chunk)
        !else
        !    write(2,"(A,i7,9x,F9.1,I4,I4)") 'maxvl-o',A,B,Ant_Stations(I_ant,i_chunk),Ant_IDs(I_ant,i_chunk)
        !endif
        !
        If(mod(Ant_IDs(I_ant,i_chunk),2) .eq. 1) then
            write(21,"(I3,I3,I7,4(1x,F8.1))") Ant_Stations(I_ant,i_chunk),Ant_IDs(I_ant,i_chunk), A, &
                Ant_RawSourceDist(I_ant,i_chunk), Ant_pos(:,I_ant,i_chunk)
        else
            write(20,"(I3,I3,I7,4(1x,F8.1))") Ant_Stations(I_ant,i_chunk),Ant_IDs(I_ant,i_chunk), A, &
                Ant_RawSourceDist(I_ant,i_chunk), Ant_pos(:,I_ant,i_chunk)
        endif
        A=MinLoc(Ant_Stations(:,i_chunk), mask=Ant_Stations(:,i_chunk).eq.Ant_Stations(I_ant,i_chunk))
        !write(2,*) I_ant,Ant_Stations(I_ant),A(1)
        If(A(1).eq.I_ant) then
            Stat_nr=Stat_nr+1
            Stations(Stat_nr)=Ant_Stations(I_ant,i_chunk)
            i_stat=Stat_nr
            !write(2,*) 'stat_nr=',Stat_nr
        else
            A=MaxLoc(Stations, mask=Stations.eq.Ant_Stations(I_ant,i_chunk))
            i_stat=A(1)
            !write(2,*) 'i_stat=',i_Stat,Stations(1:Stat_nr)
        endif
        write(txt,"(i3.3,'-'i2.2)") I_ant,i_chunk
        plot_offset=2*i_stat + 2 + mod(Ant_IDs(I_ant,i_chunk),2)
        i_c = I_Ant-11*int(I_Ant/11)
        Write(Unt,901) plot_offset,trim(DataFolder), txt,I_Ant,I_Ant,i_c
901     Format('amove 4 ',F5.2,/'begin graph',/'  size 55 2',/'  vscale 1',&
            /'  hscale 1',/'   NoBox',&
            /'   xaxis min t_min max t_max',/'   x1axis off',&
            /'   yaxis  min -scl max scl',/'   y1axis off',&
            /'     data "',A,'LOFAR_Time',A6,'.dat" d',I0,'=c1,c3',&
            /'     d',I0,' line lwidth 0.04 color MyCol',i0,'$',&
            /' end graph')
        !     let d9=0
        !     d9 line lstyle 1 lwidth 0.01 color black
        !set font pstr fontlwidth 0.06 hei 0.8 just CC
        !set lwidth 0.1
        !begin key
        !nobox
        !position tr
        !text "r=0-1 m" line color red lstyle 1 lwidth 0.1
        !end key
    Enddo
    Do i_stat=1,Stat_nr
        Call Station_ID2Mnem(Stations(I_stat),Station_Mnem)
        plot_offset=2*i_stat  + 3
        Write(unt,902) plot_offset,plot_offset,Station_Mnem
902 Format('amove 4 ',F5.2,/'aline 59 ',F5.2,/'rmove 2 0 ',/'write "',A5,'"')
    Enddo
    !
    Close(unit=unt)
    Close(unit=20)
    Close(unit=21)
    return
End Subroutine GLE_script
!=================================
Subroutine SimulationRead(SourceGuess)
!  simulation input, structured like
!   file: 'SOMELABEL_Structure.dat'  containing several lines with "STATION_ID, Station_name" like
!   001  CS001
!   002  CS002
!
!   For each station a file is expected named like: 'SOMELABEL_CS001.dat' with on
!   line 1:
!   StatStartTime, NrSamples, NrAnt, powr
!   with StatStartTime in milisec, NrSamples=length of time trace, NrAnt=nr of antennas for this station, powr=level of noise power
!   line 2 -- (NrAnt+1), each:
!   Ant_ID, LFRAnt_crdnts=antenna position in (N, E, h)
!   The following should give the timetraces for the antennas where each column corresponds to an antenna
!   line (NrAnt+2) -- (NrAnt+2+NrSamples) each having NrAnt (real) numbers
!  The antenna positions are the same for the odd and even antennas, which is why they are not doubled.
!  There are now (for each station) 4 columns giving the time traces,
!  col1: even antenna at location 1
!  col2: odd antenna at location 1
!  col3: even antenna at location 2
!  and so on.
!
   use constants, only : dp,sample
   use DataConstants, only : Time_dim, Used_StationNr
   Use Chunk_AntInfo, only : LOFAR_name, LOFARNrMax, LOFAR_number, Simulation
   use Chunk_AntInfo, only : StartT_sam
   use Chunk_AntInfo, only : Ant_Stations, Ant_IDs, Ant_nr, Ant_pos
   use Chunk_AntInfo, only : CTime_spectr, Ant_RawSourceDist
   use AntFunCconst, only : Freq_min, Freq_max
   use FFT, only : RFTransform_CF, RFTransform_CF2CT, RFTransform_CF_Filt!,Hann
   Implicit none
   Real*8, intent(in) :: SourceGuess(3) ! = (/ 8280.01,  -15120.48,    2618.37 /)     ! 1=North, 2=East, 3=vertical(plumbline)
   !
   Real(dp) :: RTime_s(1:Time_dim)
   Real*4, allocatable :: RTime_read(:)
   Complex(dp) :: CNu_s(0:Time_dim/2)  !,  CTime_s(1:Time_dim),
   Real(dp) :: nu_Fltr(0:Time_dim/2) !, Av, Bv, FiltFact
   Character(len=12) :: Lab1, Lab2, Lab3, Lab4
   Character(len=6) :: Station_name
   Character(len=5) :: LOFAR_St
   Integer :: i_file, i_ant, i_LOFAR, j_ant, i_sample, k,  AntNr_lw, AntNr_up
   Integer,save :: i_chunk=1, EvenOdd=-1
   Integer :: Cnu_dim, inu1, inu2, i_eo, STATION_ID, Ant_ID, Sample_Offset, NrAnt, NrSamples, nxx, NrAntLine
   Real*8 :: dnu, StatStartTime, SubSample_Offset, LFRAnt_crdnts(3), RDist, T_Offset, Powr !,StatAnt_Calib !,StartTime_ms
   Real(dp) :: Spec(64)
   Character(len=4),save :: Signature='----'
   Logical :: DualPosOnly=.false.
   !INTEGER, PARAMETER :: SEEK_SET = 0, SEEK_CUR = 1, SEEK_END = 2
   INTEGER ::  Unt, Start_read, final_read, NrSamples_prev, StatAntNr! , Read_offset, ierr
   !
   Open(unit=14,STATUS='old',ACTION='read', FILE = 'files/'//TRIM(Simulation)//'_Structure_LBA.dat')
   Cnu_dim=Time_dim/2
   dnu=100./Cnu_dim   ! [MHz] Jones matrix is stored on 1MHz grid
   inu1=Int(Freq_min/dnu)+1
   !inu2=Int(Freq_max/dnu-0.5)+1  ! set integration regime to frequencies within filter band
   inu2=Int(Freq_max/dnu)-1
   nu_Fltr(:)=0.
   nu_Fltr(inu1:inu2)=1.
   !
   Ant_nr(i_chunk)=0
   LOFAR_name(:)='CXnnn'
   LOFAR_number(:)=0
   i_LOFAR=0
   LOFAR_name(:)='X'
   LOFAR_number(:)=0
   Do i_file=1,2*LOFARNrMax! 10 !80
      read(14,*,IOSTAT=nxx) STATION_ID, Station_name, StatAntNr, NrSamples !LOFAR_name(i_file)
      If(nxx.ne.0) exit
      If(StatAntNr.eq.0) cycle
      ! To allow for immaginative station names, however still ending with L or H
      LOFAR_St(1:5)=Station_name(1:5)
      If( COUNT(LOFAR_name(1:i_LOFAR) .eq. LOFAR_St) .eq. 0) Then  ! add this station to the list
         i_LOFAR=i_LOFAR+1
         LOFAR_name(i_LOFAR)=LOFAR_St
         LOFAR_number(i_LOFAR)=STATION_ID/10
      EndIf
      !Open(unit=12,STATUS='old',ACTION='read', Form='unformatted', &
      !   FILE = 'files/'//TRIM(Simulation)//'_'//Station_name//'.udt')
      If(i_file.eq.1) Then
         write(2,*) 'Reading simulated data from files like= "', &
               'files/'//TRIM(Simulation)//'_'//Station_name//'.udt','"'
      EndIf
      !
      EvenOdd=1
      Signature='Dual'
      !Read(12)  StatStartTime,  NrSamples,  NrAnt    ![ms],,&  =noise power
      powr=1.
      !write(2,*) '!SimulationRead-1:', STATION_ID, Station_name, StatAntNr, NrSamples
      !Flush(unit=2)
      !
      If( i_file .eq. 1 ) Then
         Allocate( RTime_read(1:NrSamples))
         NrSamples_prev=NrSamples
      ElseIf(NrSamples_prev .ne. NrSamples) Then
         write(2,*) '********** NrSamples should be equal for all traces:', NrSamples, NrSamples_prev
         stop 'NrSamples unequal in Simulation read'
      EndIf
      !
      Unt=12
      OPEN( Unit=Unt, STATUS='unknown',ACTION='read', Form='unformatted', &
         FILE = 'files/'//TRIM(Simulation)//'_'//Station_name//'.udt')
      AntNr_lw=   Ant_nr(i_chunk)+1
      Call ReadSimulStatHead(Unt, SourceGuess, STATION_ID, StatAntNr, StatStartTime, &
         Ant_Stations(AntNr_lw,i_chunk), Ant_IDs(AntNr_lw,i_chunk), Ant_pos(1,AntNr_lw,i_chunk), &
         Ant_RawSourceDist(AntNr_lw,i_chunk), NrAnt, DualPosOnly)
      !
      !write(2,*) '!SimulationRead-2,', AntNr_lw, StatAntNr, NrAnt, i_chunk, STATION_ID, NrSamples
      !write(2,*) '!SimulationRead-3,', Ant_Stations(AntNr_lw:AntNr_lw+StatAntNr-1,i_chunk)
      !write(2,*) '!SimulationRead-4,', Ant_nr(i_chunk), StatStartTime
      !Flush(unit=2)
      !
      Do j_ant=1,StatAntNr  ! read spectra
      i_ant=Ant_nr(i_chunk)+j_ant
         RDist=Ant_RawSourceDist(i_ant,i_chunk)
         T_Offset=-RDist + StatStartTime/(Sample*1000.d0) -StartT_sam(i_chunk)  ! in units of samples
         !write(2,*) 'T_Offset', T_Offset,RDist, StatStartTime/(Sample*1000.d0)-StartT_sam(i_chunk) , i_ant
         !flush(Unit=2)
         Sample_Offset = INT(T_Offset) ! in units of sample size
         SubSample_Offset = T_Offset - Sample_Offset ! in units of sample size
         !write(2,*) 'Sample_Offset',Sample_Offset,j_ant,NrSamples,Time_dim
         If( (1+Sample_Offset).lt.1) Then
            write(2,*) 'For ',TRIM(Station_name),' first sample ',(1+Sample_Offset), &
               ' too low, reading only from 1 till',NrSamples+Sample_Offset
            If(NrSamples.gt.Time_dim/2) write(2,*) 'warning messages are created (but may be ignored)',&
               ' since the simulated trace is of the order of the internal buffer size',NrSamples,Time_dim
            If((NrSamples+Sample_Offset).lt.100) Then
               write(2,*) 'skip this Station completely (less than 100 readable samples)'
               write(*,*) 'Station ',TRIM(Station_name),' dropped'
               NrAnt=0
               exit
            EndIf
         EndIf
         If( (NrSamples+Sample_Offset).gt.Time_dim) Then
            write(2,*)  'For ',TRIM(Station_name),' last sample ',(NrSamples+Sample_Offset), &
               ' too high, reading only from',Sample_Offset,' till',Time_dim
            If(NrSamples.gt.Time_dim/2) write(2,*) 'warning messages are created (but may be ignored)',&
               ' since the simulated trace is of the order of the internal buffer size',NrSamples,Time_dim
            If((Sample_Offset+100).gt.Time_dim) Then
               write(2,*) 'skip this Station completely (less than 100 readable samples)'
               write(*,*) 'Station ',TRIM(Station_name),' dropped'
               NrAnt=0
               exit
            EndIf
         EndIf
         !
         RTime_s(:)=0.
         Start_read = 1-Sample_Offset
         If(Start_read.lt.1) Start_read=1
         final_read = Time_dim-Sample_Offset
         If(final_read .gt. NrSamples) final_read=NrSamples
         If(Start_read .gt. final_read-50) Then
            write(2,*) 'Too few samples to read; Start, final:',Start_read, final_read, Sample_Offset, &
                  j_ant, STATION_ID, Station_name
            stop 'Too few samples in simulation read'
         EndIf
         read(Unt,IOSTAT=nxx) RTime_read(1:NrSamples)
         RTime_s(Start_read+Sample_Offset: final_read+Sample_Offset)=RTime_read(Start_read: final_read)/sqrt(Powr)  ! Since power level is normalized to 100. in rest of code
         !write(2,*) 'total power:', RTime_read(Start_read:Start_read+4)
         !write(2,*) 'total power:', RTime_read(Start_read), RTime_read(final_read), &
         !   sum(RTime_read(Start_read: final_read)**2)/(final_read-Start_read), Powr
         Call RFTransform_CF_Filt(RTime_s,nu_fltr,-SubSample_Offset,Cnu_s)
         Call RFTransform_CF2CT(Cnu_s,CTime_spectr(1,i_ant,i_chunk) )
         !write(2,*) 't-trace:', i_ant,i_chunk, CTime_spectr(Start_read,i_ant,i_chunk), CTime_spectr(final_read,i_ant,i_chunk)
         !
      Enddo ! j_ant=1,NrAnt
      Ant_nr(i_chunk)=Ant_nr(i_chunk) + StatAntNr
      !
   9  Continue
      Close(Unit=Unt)

   Enddo !  i_file=1,StationNrMax
   Close(unit=14)
   DeAllocate( RTime_read )
   If(Ant_nr(i_chunk).lt.8) Then
      write(2,*) 'number of antennas=',Ant_nr(i_chunk),' less than lower limit.'
      stop 'too few antennas'
   EndIf
   !
   ! stop
   !
   Return
   !
End Subroutine SimulationRead
!++++++++++++++++++++++++++++++++++++=
   Subroutine ReadSimulStatHead(Unt, SourceGuess, STATION_ID, StatAntNr,  StatStartTime, &
      Ant_Stations, Ant_IDs, Ant_pos, Ant_RawSourceDist, NrAnt, DualPosOnly)
   use constants, only : dp,sample
   !use DataConstants, only : Time_dim, Used_StationNr
   !Use Chunk_AntInfo, only : LOFAR_name, LOFARNrMax, LOFAR_number, Simulation
   !use Chunk_AntInfo, only : StartT_sam
   Implicit none
   Integer, intent(in) :: Unt, StatAntNr,  STATION_ID
   Real(dp), intent(in) :: SourceGuess(1:3)
   Real(dp), intent(out) :: StatStartTime, Ant_pos(1:3,1:StatAntNr), Ant_RawSourceDist(1:StatAntNr)
   Integer, intent(out) :: Ant_Stations(1:StatAntNr), Ant_IDs(1:StatAntNr), NrAnt
   Logical, intent(in) :: DualPosOnly
   Integer :: NrtSamples, Ant_ID, i_ant, i_eo, nxx
   Character(len=4) :: EvOd
   Real(dp) :: LFRAnt_crdnts(1:3), RDist
   !
   ! Structure of this unformatted file
   ! rec 1: StatStartTime [ms],  NtSamples, Ant_nr [=number of different antennas (2 for each dipole)]
   !  followed by 'j_ant=1:ant_nr' data for each antenna
   ! rec 2..N: ('Dual' or 'Even' or 'Odd'), Ant_IDs(j_ant), Ant_pos(:,j_ant) ! info for dipole 'j_ant'
   ! rec N+j_ant: time trace antenna 'j_ant' (real*4)
   ! these records are repeated.
   ! at present all dipoles are necessarily 'Dual'
   !
   Read(Unt)  StatStartTime,  NrtSamples,  NrAnt    ![ms],,&  =noise power
   !IF(NrAnt.ne.StatAntNr) Then
   !   write(2,*) 'Mismatch in number of antennas',NrAnt,StatAntNr,' for this station: ',STATION_ID
   !   Stop 'ReadSimulStatData'
   !EndIf
   i_ant=0
   Do
      Read (Unt,IOSTAT=nxx) EvOd, Ant_ID, LFRAnt_crdnts
      !write(2,*) i_ant, STATION_ID, EvOd, Ant_ID, LFRAnt_crdnts
      If(DualPosOnly) Then
         If(nxx.ne.0) exit
         If((EvOd.eq."Even") .or. (EvOd.eq."Odd ")) cycle
         If(EvOd.ne."Dual") exit
      EndIf
      Call RelDist(SourceGuess,LFRAnt_crdnts,RDist)
      i_ant=i_ant+1
      Ant_Stations(i_ant)=STATION_ID
      Ant_IDs(i_ant )=Ant_ID
      Ant_pos(:, i_ant )=LFRAnt_crdnts(:)
      Ant_RawSourceDist(i_ant)=RDist         ! units of samples
      i_eo=modulo(Ant_ID,2)
      If( (i_eo.eq.1) .and. (EvOd .ne. 'Odd ')) Then
         Write(2,*) '!ReadSimulStatData: parity problem with antenna', EvOd, Ant_ID
         Stop 'ReadSimulStatData'
      EndIf
      !
      If(EvOd .eq. 'Dual') Then
         i_ant=i_ant+1
         Ant_Stations(i_ant )=STATION_ID
         Ant_IDs(i_ant )=Ant_ID+1
         Ant_pos(:, i_ant )=LFRAnt_crdnts(:)
         Ant_RawSourceDist(i_ant )=RDist         ! units of samples
      EndIf
      NrAnt=i_ant
      If(i_ant.ge.StatAntNr) exit
   EndDo
   return
End Subroutine ReadSimulStatHead
!+++++++++++++++++=======================================
!However, if we only know the size of each item but not the size of the records on a file we can use recl=1 on the OPEN statement to have the I/O list itself determine how many items to read:
!
!Example: Direct-access read, variable-length records, recl=1:
!
!demo% cat Direct3.f
!	integer u, v, w, x, y, z
!	open( 1, access='DIRECT', recl=1 )
!	read( 1, rec=1 ) u, v, w
!	read( 1, rec=13 ) x, y, z
!	write(*,*) u, v, w, x, y, z
!	end
!demo% f77 -silent Direct3.f
!demo% a.out
!------------------
!    open(unit=1, file='bytes.bin', access='stream', status='replace', &
!       & action='write', iostat=status)
!
!    write(1, iostat=status) i8
!    inquire(1, POS=my_pos)
!    print *, "Second byte will be written at the position", my_pos
!    write(1, iostat=status) i8
!    write(1, iostat=status) i8
!    write(1, iostat=status) i8
!
!    inquire(1, POS=final_pos)
!    print *, "Temp final position (unwritten) ", final_pos

!INTEGER :: file_pos
!INQUIRE(UNIT=iu, POS=file_pos)   !  stream access
!READ (iu,"()", ADVANCE='NO', POS=file_pos)
!
! ========================
