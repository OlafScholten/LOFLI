Subroutine AntennaRead(i_chunk,SourceGuess)
!  v18 : Normalization is changed
   use constants, only : dp,sample ! ,Refrac,pi,ci
   use DataConstants, only : Time_dim, Cnu_dim, Diagnostics, Production, OutFileLabel, RunMode  ! , ChunkNr_dim
   !use DataConstants, only : Polariz
   !use FitParams, only : Explore
    use Chunk_AntInfo, only : ExcludedStatID, Start_time, BadAnt_nr, BadAnt_SAI, DataReadError, AntennaNrError
    use Chunk_AntInfo, only : ExcludedStat_max, SgnFlp_nr, PolFlp_nr, SignFlp_SAI, PolFlp_SAI, BadAnt_SAI, SaturatedSamplesMax
    use Chunk_AntInfo, only : Ant_Stations, Ant_IDs, Ant_nr, Ant_NrMax, Ant_pos, CalibratedOnly
    use Chunk_AntInfo, only : CTime_spectr, Ant_RawSourceDist, Simulation, WriteSimulation
    use Chunk_AntInfo, only : NormOdd, NormEven ! Powr_eo,NAnt_eo
   !use StationMnemonics, only : Statn_ID2Mnem, Statn_Mnem2ID
   use HDF5_LOFAR_Read, only : filename, Group_names, Group_nr ! , Group_max
   use HDF5_LOFAR_Read, only : DSet_names, DSet_nr, Ant_ID, STATION_ID !, DSet_max, ANTENNA_POSITION
   use HDF5_LOFAR_Read, only : DATA_LENGTH, SAMPLE_NUMBER_first, Absolute_TIME, DIPOLE_CALIBRATION_DELAY
   !use HDF5_LOFAR_Read, only : GetFileName, ListGroups, ListDataAtt, ListGroupStructure
   use HDF5_LOFAR_Read, only : GetDataChunk
   Use StationMnemonics, only : Statn_ID2Mnem
   use FFT, only : RFTransform_CF, RFTransform_CF2CT, RFTransform_CF_Filt,Hann
   Use Calibration, only : Station_ID2Calib
   Implicit none
   Integer, intent(in) :: i_chunk
   logical :: file14open=.false.
   Real*8, intent(in) :: SourceGuess(3) ! = (/ 8280.01,  -15120.48,    2618.37 /)     ! 1=North, 2=East, 3=vertical(plumbline)
   !
   Integer*2 :: Chunk(1:Time_dim)
   Real(dp) :: RTime_s(1:Time_dim)
   Complex(dp) :: CNu_s(0:Cnu_dim), CTime_s(1:Time_dim)
   integer :: Dset_offset, Sample_Offset, DataReadErr, NAnt_eo(0:1)
   Real*8 :: Powr, Powr_eo(0:1)
   Real(dp) :: nu_Fltr(0:Cnu_dim) !, Av, Bv, FiltFact
   Integer, save :: WallCount
   Real,save :: CPUstartTime, CPUstopTime=-1., WallstartTime=0, WallstopTime=-1.
   Real,save :: TotCPURead=0., TotCPUFit=0., TotWallRead=0., TotWallFit=0.
   !
   Integer :: j,NZero, N_one, ChMax, ChMin,i_file,i_grp,i_dst, ir_file,ir_grp,ir_dst
   Integer :: i_ant, i_time, i_eo, unt, nxx, StAntID, sgn, Calibrated,  AntNr_lw, AntNr_up
   Real*8 :: SubSample_Offset, LFRAnt_crdnts(3), RDist, T_Offset,StatAnt_Calib !,StartTime_ms
   Logical :: Dubbel
   Logical, save :: FirstPass=.true.
   !character*80 :: lname
   !
   ! time recording
   call cpu_time(CPUstartTime)
   CALL SYSTEM_CLOCK(WallCount, powr)  !  Wallcount_rate)
   WallstartTime=WallCount/powr
   If(WallstopTime.lt.0) WallstopTime=WallstartTime
   If(CPUstopTime.lt.0) CPUstopTime=CPUstartTime
   !WRITE(*,"(A,F9.3,A)") 'Wall clock time:',WallstartTime, '[s]'  ! count_rate, count_max
   TotCPUFit=TotCPUFit - CPUstopTime + CPUstartTime
   TotWallFit=TotWallFit - WallstopTime +WallstartTime
   !WRITE(*,"(A,F9.3,F12.6,A)") 'Totals Fit, Wall & CPU:',TotCPUFit, TotWallFit, '[s]'  ! count_rate, count_max
   !
   !
   !Source_Crdnts= (/ 10000 , 16000 , 4000 /)    ! 1=North, 2=East, 3=vertical(plumbline)
   Write(2,"(A,F12.6,A,I11,A,3F10.1,A)") 'Start time for this chunk is set at ',Start_time(i_chunk)*1000.d0*sample &
      ,' [ms] =',Start_time(i_chunk),'Samples, SourceGuess=',SourceGuess,';  1=North, 2=East, 3=vertical(plumbline)'
   !Write(2,*) 'SourceGuess=',SourceGuess,';  1=North, 2=East, 3=vertical(plumbline)'
   write(*,"(A,f8.2,A)") achar(27)//'[33m @',Start_time(i_chunk)*1000.*sample,'[ms]'//achar(27)//'[0m'  ! [1000D    !  //achar(27)//'[0m.'
   ! In main program, after option selection:
   Inquire(unit=14, opened=file14open)
   If((Simulation.eq."")) WriteSimulation(2)=-1
   !write(2,*) 'Simulation,WriteSimulation(2):', Simulation,WriteSimulation(:), 'file14open:',file14open
   If((Simulation.eq."") .or. (WriteSimulation(2).gt.0)) Then
      If(file14open) then
         Rewind(unit=14)
         Rewind(unit=12)
      Else
         Open(unit=12,STATUS='unknown',ACTION='read',FORM ="unformatted", FILE = 'Book/RFI_Filters-v18.uft')
         Open(unit=14,STATUS='old',ACTION='read', FILE = 'Book/LOFAR_H5files_Structure-v18.dat')
      EndIf
   Else  ! simulation ne '' and width <0; run on simulation data
      If(file14open) then
         write(2,*) '****Error in SimulationRead, file 14 already open'
         Stop ' SimulationRead'
      EndIf
      Open(unit=14,STATUS='old',ACTION='read', FILE = 'files/'//TRIM(Simulation)//'_Structure.dat')
      Call SimulationRead(SourceGuess)
      goto 9
   Endif
   If((Simulation.ne."") .and. (WriteSimulation(2).gt.0)) Then
      Call CreateNewFolder(Simulation) ! create new folder when needed
      Open(Unit=30,STATUS='unknown',ACTION='write', FILE = 'files/'//TRIM(Simulation)//'_Structure.dat')
      write(2,*) 'created simulation file:','files/'//TRIM(Simulation)//'_Structure.dat'
      !write(*,*) 'created simulation file:','files/'//TRIM(Simulation)//'_Structure.dat'
   EndIf
   Ant_nr(i_chunk)=0
   AntNr_lw=1
   ir_file=0 ; ir_grp=0 ; ir_dst=0
   Powr_eo(:)=0.  ; NAnt_eo(:)=0
   Do i_file=1,130 ! 10 !80
   !i_file=1 ; rewind(unit=14)
      read(14,*,IOSTAT=nxx)  Group_nr, filename
      If(nxx.ne.0) exit
      If(FirstPass) Then
         FirstPass=.false.
         write(2,*) 'First data file #=',i_file,', name=',filename
      EndIf
      !If(Production) write(2,*) 'file #=',i_file,', name=',filename
      If(.not. Production) write(*,"(A,i3)", ADVANCE='NO') achar(27)//'[100D'//'file#=',i_file  ! [1000D    !  //achar(27)//'[0m.'
      Do i_grp=1,Group_nr
          read(14,*) DSet_nr, Group_Names(i_grp)
          If(Diagnostics) Write(2,"(A,i0,A,A)") 'Group#',i_grp,' : ',trim(Group_names(i_grp))
          Do i_dst=1,DSet_nr
              read(14,*) DSet_Names(i_dst), STATION_ID, Ant_ID
              !write(*,*) DSet_nr, Group_Names(i_grp), DSet_Names(i_dst), STATION_ID, Ant_ID
              !If(STATION_ID .eq. 150) cycle  !  throuw out RS310  ****************************
              If( COUNT(ExcludedStatID .eq. STATION_ID) .gt. 0) cycle
              StAntID= STATION_ID * 1000 + Ant_ID
              i_eo=Mod(Ant_ID,2)
              If( COUNT((PolFlp_SAI(1:PolFlp_nr)+i_eo) .eq. StAntID) .gt. 0) then
                  If(Diagnostics) &
                     write(2,*) '****Pole-Flip antenna:',STATION_ID,Ant_ID
                  ! cycle
                  Ant_ID=Ant_ID+1-2*i_eo
                  StAntID= STATION_ID * 1000 + Ant_ID
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
              !    Write(2,*) 'reading Filter data',i_file,i_grp,i_dst,ir_file,ir_grp,ir_dst
                  Read (12) ir_file,ir_grp,ir_dst, LFRAnt_crdnts, powr,&
                  DATA_LENGTH, SAMPLE_NUMBER_first, Absolute_TIME, DIPOLE_CALIBRATION_DELAY,nu_fltr !nu_fltr
              enddo
              if(powr.lt.1.d-6) cycle
              If((i_file .ne. ir_file) .or. (i_grp .ne. ir_grp) .or. (i_dst .ne. ir_dst) ) then
                  Write(2,*) '****Filter data not found',i_file,i_grp,i_dst,ir_file,ir_grp,ir_dst
                  cycle   ! get next DSet
              Endif
              !
              Powr_eo(i_eo)=Powr_eo(i_eo)+Powr
              NAnt_eo(i_eo)=NAnt_eo(i_eo)+1
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
              If(CalibratedOnly .and. Calibrated.eq.0) Then
               !write(2,*) 'STATION_ID,Ant_ID was not calibrated',STATION_ID,Ant_ID,StatAnt_Calib
               cycle
              EndIf
              !
              !Absolute_TIME  ! should be the same for all
              !SAMPLE_NUMBER_first  ! Number of samples, after 'Absolute_TIME' for the first recording
              !IF(Ant_ID.ge.90 .and. Ant_ID.le.93) write(2,*) 'StatAnt_Calib',StAntID, &
              !          DIPOLE_CALIBRATION_DELAY/Sample, StatAnt_Calib, DIPOLE_CALIBRATION_DELAY/Sample + StatAnt_Calib
              !If(i_chunk.eq.1) Write(2,*) 'AntennaRead, timecorrection',StAntID, &
              !    DIPOLE_CALIBRATION_DELAY/Sample, StatAnt_Calib, DIPOLE_CALIBRATION_DELAY/Sample+StatAnt_Calib
              StatAnt_Calib=DIPOLE_CALIBRATION_DELAY/Sample + StatAnt_Calib  ! in units of samples
              T_Offset=RDist + StatAnt_Calib  ! in units of samples
              Sample_Offset = INT(T_Offset) ! in units of sample size
              SubSample_Offset = T_Offset - Sample_Offset ! in units of sample size
              Dset_offset=Start_time(i_chunk) + Sample_Offset - SAMPLE_NUMBER_first
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
               write(2,*) '************ error in data-read!, file #=',i_file,', name=',filename
               DataReadError=-1
               write(*,*) 'error in data-read'
               return
              endif
              !Call ListDataAtt(Group_Names(i_grp),DSet_Names(i_dst), Chunk, Dset_offset, Time_dim)
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
                     write(2,"(A,I4,', chunk#',I4,A,I3,', @')",ADVANCE='NO') &
                           Statn_ID2Mnem(STATION_ID), Ant_ID, i_chunk, ', #Saturated Samples=',N_one
                     Do i_time=1,Time_dim
                        If(Chunk(i_time).gt.2020 .or. Chunk(i_time).lt.-2020) Then
                           Write(2,"(1x,I6,':',I5)",ADVANCE='NO') i_time, Chunk(i_time)
                        EndIf
                     EndDo
                     Write(2,*) ' '
                     cycle
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
              If(Ant_nr(i_chunk) .gt. Ant_nrMax) then
                  Write(2,*) '****Max nr of antennas exceeded@ file,grp,dst',i_file,i_grp,i_dst,'Ant_nr',Ant_nr(i_chunk)
                  stop 'Antenna numbers exceeded'
              Endif
              !
              Ant_IDs(Ant_nr(i_chunk),i_chunk)=Ant_ID
              Ant_Stations(Ant_nr(i_chunk),i_chunk)=STATION_ID
              Ant_pos(:,Ant_nr(i_chunk),i_chunk)=LFRAnt_crdnts(:)
              !
              sgn=+1
              If(any(SignFlp_SAI(1:SgnFlp_nr) .eq. StAntID,1) ) sgn=-1 ! Take care of sign flips
              If(WriteSimulation(2).le.0) Then
                 If(RunMode.eq.1) then
                  RTime_s(:)=Chunk(:)*Hann(:)*sgn
                 Else
                  RTime_s(:)=Chunk(:)*Hann(:)*sgn*100./sqrt(Powr) !  sqrt(power) level is normalized to 100.
                 Endif
              Else
                 RTime_s(WriteSimulation(1):WriteSimulation(1)+WriteSimulation(2))= &
                     Chunk(WriteSimulation(1):WriteSimulation(1)+WriteSimulation(2))*sgn*100./sqrt(Powr) ! /Time_dim)
              EndIf
              !SubSample_Offset=0.
              Call RFTransform_CF_Filt(RTime_s,nu_fltr,SubSample_Offset,Cnu_s)
              Call RFTransform_CF2CT(Cnu_s,CTime_spectr(1,Ant_nr(i_chunk),i_chunk) )
              !    Av=MaxVal(RTime_s)
              !    Bv=MinVal(RTime_s)
              !    write(2,*) 'Normalized: min & max=', Bv, Av
              ! write(2,*) 'STATION_ID,Ant_ID',STATION_ID,Ant_ID,StatAnt_Calib
              !
          1   continue
          EndDo  ! i_dst=1,DSet_nr
          ! re-order neighboring antennas according to increasing antenna ID (assumed in interferometric imaging)
          ! re-order : CTime_spectr(1,i_ant,i_chunk) ; Ant_IDs(i_ant,i_chunk) ; Ant_pos(:,i_ant,i_chunk), only for this STATION_ID
          Do i_ant=2,Ant_nr(i_chunk)
            If(Ant_Stations(i_ant-1,i_chunk).eq.Ant_Stations(i_ant,i_chunk)) Then ! apply to antennas from a single station only
               If(Ant_IDs(i_ant-1,i_chunk).gt.Ant_IDs(i_ant,i_chunk)) Then  ! flip order
                  !write(2,*) 'change order ',i_ant-1,Ant_IDs(i_ant-1,i_chunk),' and ',i_ant,Ant_IDs(i_ant,i_chunk), &
                  !   ' for station',Ant_Stations(i_ant,i_chunk)
                  CTime_s(1:Time_dim)=CTime_spectr(1:Time_dim,i_ant,i_chunk)
                  Ant_ID=Ant_IDs(i_ant,i_chunk)
                  LFRAnt_crdnts(:)=Ant_pos(:,i_ant,i_chunk)
                  CTime_spectr(1:Time_dim,i_ant,i_chunk)=CTime_spectr(1:Time_dim,i_ant-1,i_chunk)
                  Ant_IDs(i_ant,i_chunk)=Ant_IDs(i_ant-1,i_chunk)
                  Ant_pos(:,i_ant,i_chunk)=Ant_pos(:,i_ant-1,i_chunk)
                  CTime_spectr(1:Time_dim,i_ant-1,i_chunk)=CTime_s(1:Time_dim)
                  Ant_IDs(i_ant-1,i_chunk)=Ant_ID
                  Ant_pos(:,i_ant-1,i_chunk)=LFRAnt_crdnts(:)
               EndIf
            EndIf
          Enddo
          If((Simulation.ne."") .and. (WriteSimulation(2).gt.0)) Then
               write(30,*) STATION_ID,TRIM(Statn_ID2Mnem(STATION_ID))
               AntNr_up=Ant_nr(i_chunk)
               !Call CreateNewFolder(Simulation) ! create new folder when needed
               Open(Unit=31,STATUS='unknown',ACTION='write', &
                     FILE = 'files/'//TRIM(Simulation)//'_'//TRIM(Statn_ID2Mnem(STATION_ID))//'.dat')
               Write(31,*) 'StartTime_ms= ',(Start_time(i_chunk)+WriteSimulation(1)+RDist)*sample*1000., &
                     ' N_samp= ', WriteSimulation(2), ' noise_power= 100.'
               Do i_ant=AntNr_lw,AntNr_up
                  i_eo=mod(Ant_IDs(i_ant,i_chunk),2)
                  If(i_eo.eq.0) then
                     If(Ant_IDs(i_ant+1,i_chunk).eq.(Ant_IDs(i_ant,i_chunk)+1)) Then
                        Dubbel=.true.
                     Else
                        Write(31,*) 'Even ', 'NEh= ',Ant_pos(:,i_ant,i_chunk)  !  Ant_IDs(i_ant,i_chunk),
                        cycle
                     EndIf
                  Endif
                  If(i_eo.eq.1) Then
                     If(Dubbel) Then
                        Write(31,*) 'Dual ', 'NEh= ', Ant_pos(:,i_ant,i_chunk)
                        Dubbel=.false.
                     Else
                        Write(31,*) 'Odd  ', 'NEh= ', Ant_pos(:,i_ant,i_chunk)
                     EndIf
                  EndIf
               EndDo
               !Write(31,*) '===   0 0. 0. 0.'
               !write(2,*) WriteSimulation(1),(WriteSimulation(1)+WriteSimulation(2)),AntNr_lw,AntNr_up,i_chunk
               Do i_time=WriteSimulation(1),(WriteSimulation(1)+WriteSimulation(2))
                  write(31,*) Real( CTime_spectr(i_time,AntNr_lw:AntNr_up,i_chunk) )
               EndDo
               AntNr_lw=AntNr_up+1
               Close(Unit=31)
          EndIf
          !
      Enddo ! i_grp=1,Group_nr
!   !
!   Close(unit=22)
!   !
      !
   Enddo ! i_file=1,....
!   Close(unit=24)
   !
   If((Simulation.ne."") .and. (WriteSimulation(2).gt.0)) Then
      Close(Unit=30)
      write(2,*) 'true data written to files to be used in sumulation runs'
      Stop 'simulation'
   EndIf
   !
   If(Ant_nr(i_chunk).le.6) then
      write(2,*) Ant_nr(i_chunk),' are too few active antennas for chunck starting at ',Start_time(i_chunk)
      If(AntennaNrError.lt.0) then
         stop 'AntennaRead'
      Else
         AntennaNrError=1
      Endif
   Endif
   !
   write(*,*) ' '
   Powr_eo(:)=Powr_eo(:)/NAnt_eo(:)
   NormEven=sqrt(2.*Powr_eo(0)/(Powr_eo(0)+Powr_eo(1)))/100.  ! to undo the factor 100 (and more) that was introduced in antenna-read
   NormOdd=sqrt(2.*Powr_eo(1)/(Powr_eo(0)+Powr_eo(1)))/100.  ! to undo the factor that was introduced in antenna-read
9  Continue
   write(2,*) 'ratio of average background amplitude for all even/odd antennas',NormEven/NormOdd
   !write(2,*) 'NAnt_eo:',NAnt_eo(:)
   !
   ! time recording
   call cpu_time(CPUstopTime)
   CALL SYSTEM_CLOCK(WallCount, powr)  !  Wallcount_rate)
   !WRITE(*,*) 'Wall clock time:',WallCount, ', countrate:', powr  ! count_rate, count_max
   WallstopTime=WallCount/powr
   !WRITE(*,"(A,F9.3,A)") 'Wall clock time:',WallstartTime, '[s]'  ! count_rate, count_max
   TotCPURead=TotCPURead + CPUstopTime-CPUstartTime
   TotWallRead=TotWallRead + WallstopTime-WallstartTime
   !WRITE(*,"(A,F9.3,F12.6,A)") 'Totals Read, Wall & CPU:',TotCPURead, TotWallRead, '[s]'  ! count_rate, count_max
   !WRITE(*,"(A,F9.3,F12.6,A)") 'Totals Fit, Wall & CPU:',TotCPUFit, TotWallFit, '[s]'  ! count_rate, count_max
   WRITE(*,"(A,F10.3,F13.6,A,F10.3,F13.6,A)") 'TotWall&CPU; Read:', TotWallRead,TotCPURead,'[s], Fit:', TotWallFit,TotCPUFit,'[s]'
   !
   !
   ! stop
   !
   Return
   !
End Subroutine AntennaRead
! ========================
Subroutine PlotSpectra(i_chunk)
   use constants, only : dp,sample ! ,Refrac,pi,ci
   use DataConstants, only : Time_dim, Cnu_dim, DataFolder, Diagnostics ! , ChunkNr_dim
   use Chunk_AntInfo, only : Ant_nr, CTime_spectr
   !use StationMnemonics, only : Statn_ID2Mnem, Statn_Mnem2ID
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
    use DataConstants, only : Station_nrMax
    use StationMnemonics, only : Station_ID2Mnem
    use ThisSource, only : PeakPos
    Implicit none
    Integer, intent(in) :: unt, Ant_nr, i_chunk
    Character(Len=30), intent(in) :: file
    integer :: I_ant, i_stat, Stat_nr, Stations(1:Station_nrMax), i_time, i_c
    integer :: A(1), Plot_scale,ch1,ch2,VLine
    Character(len=6) :: txt
    Character(len=5) :: Station_Mnem
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
   use constants, only : dp,sample ! ,Refrac,pi,ci
   use DataConstants, only : Time_dim, Cnu_dim, Station_nrMax
   Use Chunk_AntInfo, only : Station_name, Station_number, Simulation
   use Chunk_AntInfo, only : Start_time
   use Chunk_AntInfo, only : Ant_Stations, Ant_IDs, Ant_nr, Ant_NrMax, Ant_pos
   use Chunk_AntInfo, only : CTime_spectr, Ant_RawSourceDist
   use Chunk_AntInfo, only : NormOdd, NormEven !Powr_eo,NAnt_eo
   use AntFunCconst, only : Freq_min, Freq_max
   !use StationMnemonics, only : Statn_ID2Mnem, Statn_Mnem2ID
   use FFT, only : RFTransform_CF, RFTransform_CF2CT, RFTransform_CF_Filt!,Hann
   Implicit none
   Real*8, intent(in) :: SourceGuess(3) ! = (/ 8280.01,  -15120.48,    2618.37 /)     ! 1=North, 2=East, 3=vertical(plumbline)
   !
   Real(dp) :: RTime_s(1:Time_dim)
   Complex(dp) :: CNu_s(0:Cnu_dim)  ! CTime_s(1:Time_dim),
   Real(dp) :: nu_Fltr(0:Cnu_dim) !, Av, Bv, FiltFact
   Character(len=12) :: Lab1, Lab2, Lab3, Lab4
   Integer :: i_file, i_ant, j_ant, i_sample, k,  AntNr_lw, AntNr_up
   Integer,save :: i_chunk=1, EvenOdd=-1
   Integer :: inu1, inu2, i_eo, STATION_ID, Ant_ID, Sample_Offset, NrAnt, NrSamples, nxx, NrAntLine
   Real*8 :: dnu, StatStartTime, SubSample_Offset, LFRAnt_crdnts(3), RDist, T_Offset, Powr !,StatAnt_Calib !,StartTime_ms
   Real(dp) :: Spec(64)
   Character(len=4),save :: Signature='----'
   !
   !Open(unit=14,STATUS='old',ACTION='read', FILE = TRIM(Simulation)//'_Structure.dat')
   !
   dnu=100./Cnu_dim   ! [MHz] Jones matrix is stored on 1MHz grid
   inu1=Int(Freq_min/dnu)+1
   !inu2=Int(Freq_max/dnu-0.5)+1  ! set integration regime to frequencies within filter band
   inu2=Int(Freq_max/dnu)-1
   nu_Fltr(:)=0.
   nu_Fltr(inu1:inu2)=1.
   Ant_nr(i_chunk)=0
   Station_name(:)='CXnnn'
   Station_number(:)=0
   !Powr_eo(:)=0.  ; NAnt_eo(:)=0
   AntNr_lw=0
   AntNr_up=0
   Do i_file=1,Station_NrMax ! 10 !80
      read(14,*,IOSTAT=nxx) STATION_ID, Station_name(i_file)
      If(nxx.ne.0) exit
      Station_number(i_file)=STATION_ID
      If(i_file.eq.1) Then
         write(2,*) 'Reading simulated data from files like= "', &
               'files/'//TRIM(Simulation)//'_'//Station_name(i_file)//'.dat','"'
      EndIf
      Open(unit=12,STATUS='old',ACTION='read', FILE = 'files/'//TRIM(Simulation)//'_'//Station_name(i_file)//'.dat')
      If(EvenOdd.eq.-1) Then
         Read(12,*) Lab1
         Read(12,*) Lab1
         read(Lab1,*,IOSTAT=nxx) Ant_ID
         If(nxx.eq.0) Then ! there was a number at this place
            EvenOdd=1
            Signature='Dual'
         Else
            EvenOdd=0
         EndIf
         rewind(unit=12)
      EndIf
      !write(2,*) 'EvenOdd',EvenOdd,Signature
      !flush(Unit=2)
      If(EvenOdd.eq.0) then
         Read(12,*) Lab1, StatStartTime, Lab2, NrSamples, lab4, powr  ![ms],,&  =noise power
         !write(2,*) Lab1, StatStartTime, Lab2, NrSamples, lab4, powr
      else
         Read(12,*) Lab1, StatStartTime, Lab2, NrSamples, lab3, NrAnt, lab4, powr  ![ms],,&  =noise power
         !write(2,*) Lab1, StatStartTime, Lab2, NrSamples, lab3, NrAnt, lab4, powr
      endIf
      if(powr.lt.1.d-36) Powr=1.d-36  ! =noise power
      NrAntLine=0
      Do
         NrAntLine=NrAntLine+1
         If(EvenOdd.eq.0) then
            Read (12,*,IOSTAT=nxx) Signature, lab1, LFRAnt_crdnts
         else
            If(NrAnt.eq.NrAntLine-1) exit
            Read (12,*,IOSTAT=nxx) Ant_ID, lab1, LFRAnt_crdnts
         endIf
         !Read (12,*,IOSTAT=nxx) Signature, Ant_ID,  LFRAnt_crdnts
         If(nxx.ne.0) exit
         Call RelDist(SourceGuess,LFRAnt_crdnts,RDist)
         AntNr_up=AntNr_up+1
         Ant_Stations(AntNr_up ,i_chunk)=STATION_ID
         Ant_IDs(AntNr_up ,i_chunk)=10*NrAntLine
         Ant_pos(:,AntNr_up ,i_chunk)=LFRAnt_crdnts(:)
         Ant_RawSourceDist(AntNr_up ,i_chunk)=RDist         ! units of samples
         If(Signature.eq.'Odd ') Then
            Ant_IDs(AntNr_up ,i_chunk)=10*NrAntLine+1
         ElseIf(Signature.eq.'Dual') Then
            AntNr_up=AntNr_up+1
            Ant_Stations(AntNr_up ,i_chunk)=STATION_ID
            Ant_IDs(AntNr_up ,i_chunk)=10*NrAntLine+1
            Ant_pos(:,AntNr_up ,i_chunk)=LFRAnt_crdnts(:)
            Ant_RawSourceDist(AntNr_up ,i_chunk)=RDist         ! units of samples
         ElseIf(Signature.ne.'Even') Then
            AntNr_up=AntNr_up-1
            exit
         EndIf
         If(AntNr_up-AntNr_lw.gt.64) then
            write(2,*) 'nr of antennas too large in SimulationRead:',AntNr_up-AntNr_lw
            stop 'SimulationRead'
         EndIf
         !write(2,*) 'AntNr_up,AntNr_lw',AntNr_up,AntNr_lw, Signature, Ant_IDs(AntNr_up ,i_chunk), NrAntLine
         !
         !flush(Unit=2)
      EndDo
      !Powr_eo(0)=Powr_eo(0)+ (AntNr_up-AntNr_lw)*Powr
      !NAnt_eo(0)=NAnt_eo(0)+ AntNr_up-AntNr_lw
      Ant_nr(i_chunk)=AntNr_up
      !write(2,*) 'Ant_IDs',Ant_IDens(1:Ant_nr(i_chunk),i_chunk)
      RTime_s(:)=0.
      Do j_ant=1,AntNr_up-AntNr_lw  ! read spectra
         Rewind(unit=12) ! reposition file
         Do k=1,NrAntLine
            read(12,*) Lab1
         EndDo
         i_ant=AntNr_lw+j_ant
         RDist=Ant_RawSourceDist(i_ant,i_chunk)
         !StatAnt_Calib=StatStartTime/Sample ! StatAnt_Calib in units of samples
         !T_Offset=-RDist - StatStartTime/(Sample*1000.) +Start_time(i_chunk)  ! in units of samples
         T_Offset=-RDist + StatStartTime/(Sample*1000.) -Start_time(i_chunk)  ! in units of samples
         !write(2,*) 'T_Offset', T_Offset,RDist, StatStartTime/(Sample*1000.)-Start_time(i_chunk) , i_ant
         !flush(Unit=2)
         Sample_Offset = INT(T_Offset) ! in units of sample size
         SubSample_Offset = T_Offset - Sample_Offset ! in units of sample size
         !write(2,*) 'Sample_Offset',Sample_Offset,j_ant
         If( (1+Sample_Offset).lt.1) write(2,*) 'too low',(1+Sample_Offset)
         If( (NrSamples+Sample_Offset).gt.Time_dim) write(2,*) 'too high',(NrSamples+Sample_Offset)
         Do i_sample=1,NrSamples
            read(12,*) spec(1:AntNr_up-AntNr_lw)
            !write(2,*) i_sample, Sample_Offset
            If( (i_sample+Sample_Offset).lt.1) cycle
            If( (i_sample+Sample_Offset).gt.Time_dim) exit
            RTime_s(i_sample+Sample_Offset)=spec(j_ant)*100./sqrt(Powr)  ! Since power level is normalized to 100. in rest of code
         Enddo
         Call RFTransform_CF_Filt(RTime_s,nu_fltr,-SubSample_Offset,Cnu_s)
         Call RFTransform_CF2CT(Cnu_s,CTime_spectr(1,i_ant,i_chunk) )
         !
      Enddo ! j_ant=1,NrAnt
      Close(unit=12)
      AntNr_lw=AntNr_up
      !
   Enddo !  i_file=1,StationNrMax
   Close(unit=14)
   If(AntNr_up.lt.8) Then
      write(2,*) 'number of antennas=',AntNr_up,' less than lower limit.'
      stop 'too few antennas'
   EndIf
   !Powr_eo(1)=Powr_eo(0)
   !NAnt_eo(1)=NAnt_eo(0)
   NormOdd=1. ; NormEven=1.
   !
   ! stop
   !
   Return
   !
End Subroutine SimulationRead
! ========================
