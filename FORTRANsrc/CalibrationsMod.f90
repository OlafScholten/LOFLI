Module Calibration
    use constants, only : dp
    Integer, parameter :: Station_nrMax=80
    Integer, parameter :: MaximumAntennas=1500
    ! Station mnemonics of stations that are read-in from existing calibration file:
    Character(len=6), save :: Fine_STMnm(1:Station_nrMax)  ! station mnemonics for which Ever Calibrations Have Been Generated (ECHBG)
    Real(dp), save :: Fine_STDelay(1:Station_nrMax)! time calibration data for stations
    Real(dp), save :: Fine_AntDelay(1:MaximumAntennas)! time calibration data for antennas
    Integer, save :: SAI_AntDelay(1:MaximumAntennas)  ! SAI of antennas for which ECHBG
    Logical, save :: AntInCal(1:MaximumAntennas)  ! =1 if this antenna was used when writing cal-file
    Integer, save :: Nr_AntDelay, Nr_StatDelay  ! Nr of antennas or stations for which ECHBG
    Logical, save :: StationInCal(Station_nrMax)  ! =1 when this station was used in calculation when writing cal-file
Contains
! -------------------------------------------------
Subroutine ReadCalib()
    use constants, only : dp,sample, HeightCorrectIndxRef
    use DataConstants, only : DataFolder, Calibrations, RunMode ! , Diagnostics
    !Use Chunk_AntInfo, only : Fine_STMnm, Fine_STDelay, Fine_AntDelay, SAI_AntDelay, Nr_AntDelay, Nr_StatDelay, AntInCal ! all these are output from this routine
    Implicit none
    Character(len=5) :: txt !, Station_Mnem
    Character(len=6) :: tHL !, Station_Mnem
    Real(dp) :: Delay
    Character(len=5), save :: LOFAR_STMnm(Station_nrMax) ! Archaic
    Real(dp), save :: LOFAR_STDelay(Station_nrMax) ! Archaic
    integer :: nxx, i, SAI
    Integer :: i_LOFAR, i_fine, n, CalVerion
    Logical :: Old, FMT2022, IndxRef
    Character(len=80) :: tst !, Station_Mnem
    Character(len=3) :: Mrk
    Character(len=25) :: date,TxtIndx
    !
   i_LOFAR=0
   Open(unit=9,STATUS='old',ACTION='read', FILE = 'Book/StationCalibrations.dat', IOSTAT=nxx) ! trim(DataFolder)//
   ! distilled from LOFAR data-base with "CalibrationTablesGenerator.f90"
   If(nxx.ne.0) then
      write(2,*) 'system station calibrations might be missing! expected file =', &
         'StationCalibrations.dat'
      !stop 'ReadCalib'
   else
      Do
          read(9,*,IOSTAT=nxx) txt, Delay
          if(nxx.ne.0) exit
              i_LOFAR=i_LOFAR+1
              LOFAR_STMnm(i_LOFAR) = txt
              LOFAR_STDelay(i_LOFAR)= Delay/sample  ! Convert to samples
!              write(2,*) '! ReadCalib:',i_LOFAR,LOFAR_STMnm(i_LOFAR), LOFAR_STDelay(i_LOFAR)
      Enddo
   endif
   close(unit=9)
   !
   Nr_AntDelay=0
   Old=(trim(Calibrations).eq.'')
   If(Old) then
      Open(unit=9,STATUS='old',ACTION='read', FILE = 'Book/AntCalibrations.dat', IOSTAT=nxx) ! trim(DataFolder)//
      If(nxx.ne.0) then
         write(2,*) 'problems with file=','Book/AntCalibrations.dat', ', specify Calibrations parameter'
         stop 'ReadCalib:old'
      EndIf
      CalVerion=0
   Else
      Open(unit=9,STATUS='old',ACTION='read', FILE = 'Book/'//trim(Calibrations), IOSTAT=nxx)
      If(nxx.ne.0) then
         write(2,*) 'problems with file=','Book/'//trim(Calibrations),', probably does not exist'
         stop 'ReadCalib'
      EndIf
      CalVerion=1
   Endif
   !
   read(9,"(A60)") tst  ! Check for a new-style calibration file
   read(tst,*,IOSTAT=nxx) Mrk, date, TxtIndx, IndxRef  ! Check for Height correction Index Ref
   !write(2,*) 'tst:',tst
   !write(2,*) 'From ReadCalib:nxx=',nxx,Mrk,i,' ', date,' ', TxtIndx,' ', IndxRef
   If(nxx.eq.0) Then
      FMT2022=.true.
      i_LOFAR=0  ! the data of  'Book/StationCalibrations.dat' (Archaic) have already been included in the station delays
      If( Mrk(2:3) .eq. '04' ) Then  ! Account for presence of HBA antennas
         CalVerion=4
         write(2,*) 'From ReadCalib, use LBA & HBA calibrations:',Mrk,' ', date,' ', TxtIndx,' ', IndxRef
      Else
         CalVerion=3
         write(2,*) 'From ReadCalib, assume only LBA calibrations::',Mrk,' ', date,' ', TxtIndx,' ', IndxRef
      Endif
   Else
      IndxRef=.false.  ! In the old days no correction for the height dependence of the index of refraction was used
      read(tst,*,IOSTAT=nxx) SAI, Delay,n  ! Check for a new-style calibration file
      If(nxx.eq.0) Then
         CalVerion=2
         FMT2022=.true.
         i_LOFAR=0  ! the data of  'Book/StationCalibrations.dat' (Archaic) have already been included in the station delays
      Else
         CalVerion=1 ! was already set as such
         FMT2022=.false.
         If(i_LOFAR.eq.0) Then
            write(2,*) 'In "ReadCalib": something went wrong reading the file StationCalibrations.dat'
            Stop 'Error in ReadCalib'
         EndIf
      EndIf
      rewind(unit=9)
   EndIf  ! file 9 is now positioned to first line with antenna data
   ! CalVerion=0, same as
   !  Old=.true. !  REALLY, REALLY old
   !
   ! CalVerion=1, same as
   !  i_LOFAR > 0    ! 'Book/StationCalibrations.dat' is needed, Name calfile can be specified
   !
   ! CalVerion=2, same as
   !  i_LOFAR=0  ! the data of  'Book/StationCalibrations.dat' (Archaic) have already been included in the station delays
   !
   ! CalVerion=3, same as
   !  FMT2022=.true. ! the 2022 version where track is kept of calibrated stations, uses i_LOFAR=0
   !
   ! CalVerion=4, same as
   !  FMT2025=.true. ! superceeds 2022 version where HBA & LBA are included
   !
   If(runmode.ne.2 .and. runmode.ne.7) Then  ! do not overwrite parameter when set on input
      HeightCorrectIndxRef= IndxRef
      write(2,*) 'runmode, HeightCorrectIndxRef:',runmode, HeightCorrectIndxRef
   EndIf
   !
 ! ReadCalib:          35 RS305    1384.3654367279998
 !ReadCalib           1          35 F F
   write(2,*) '!ReadCalib',CalVerion, i_lofar, FMT2022, old
   n=-1
   AntInCal(:)=.false.
   Do       ! Read single antenna calibrations
      If(CalVerion.ge.3) Then
         read(9,*,IOSTAT=nxx) SAI, Delay,n
      Else
         read(9,*,IOSTAT=nxx) SAI, Delay
      EndIf
      if(nxx.ne.0) exit
      Nr_AntDelay=Nr_AntDelay+1
      SAI_AntDelay(Nr_AntDelay) = SAI
      Fine_AntDelay(Nr_AntDelay)= Delay ! already in [samples]
      If(n.gt.0) AntInCal(Nr_AntDelay)=.true.
   enddo
   !
   nxx=0
   If(CalVerion.eq.0) Then
      Close(unit=9)
      Open(unit=9,STATUS='old',ACTION='read', FILE = 'Book/FineCalibrations.dat', IOSTAT=nxx) ! trim(DataFolder)//
      If(nxx.ne.0) then
         write(2,*) 'problems with file=','FineCalibrations.dat'
         stop 'ReadCalib, Reading problem old style'
      EndIf
   EndIf
   !
   Nr_StatDelay=0
   StationInCal(:)=.false.
   Do  ! get station calibrations
      read(9,"(A60)",IOSTAT=nxx) tst  ! Check for a new-style calibration file
      if(nxx.ne.0) exit
      If(CalVerion.eq.4) Then
         read(tst,*,IOSTAT=nxx) tHL, Delay, n
         If(tHL(2:2).ne.'S') exit
         Nr_StatDelay=Nr_StatDelay+1
         Fine_STMnm(Nr_StatDelay) = tHL
         Fine_STDelay(Nr_StatDelay)= Delay ! already in samples
         StationInCal(Nr_StatDelay)=.true.
         tHL=''
      Else
         read(tst,*,IOSTAT=nxx) txt, Delay, n
         !write(2,*) nxx,'"',txt,'"', Delay, n,'"',tst
         if(nxx.ne.0) n=0
         If(txt(2:2).ne.'S') exit
         Nr_StatDelay=Nr_StatDelay+1
         Fine_STMnm(Nr_StatDelay) = txt//'L'
         Fine_STDelay(Nr_StatDelay)= Delay ! already in samples
         StationInCal(Nr_StatDelay)=.true.
         Nr_StatDelay=Nr_StatDelay+1
         Fine_STMnm(Nr_StatDelay) = txt//'H'
         Fine_STDelay(Nr_StatDelay)= Delay ! already in samples
         StationInCal(Nr_StatDelay)=.false.
         txt=''  !  to catch empty in-read line
      EndIf
   enddo
   close(unit=9)
   !write(2,*) 'Nr_StatDelay:', Nr_StatDelay, Fine_STMnm(Nr_StatDelay-3:Nr_StatDelay)
   !
   If(CalVerion.le.3) Then  ! need to do something for HBA
      If(CalVerion.le.2) Then  ! Archaic delays need to be merged
         Do i_fine=1,Nr_StatDelay    ! Write stationdelays
            Do i=1,i_LOFAR  ! include old station delays from a previous age
               If(LOFAR_STMnm(i).eq.Fine_STMnm(i_fine)(1:5)) then
                  Fine_STDelay(i_fine)  =Fine_STDelay(i_fine) + LOFAR_STDelay(i)
                  exit
               EndIf
            EndDo
         EndDo
      EndIf
   EndIf
   !write(2,*) 'Nr_AntDelay=',Nr_AntDelay
   Return
End Subroutine ReadCalib
!=================================
Subroutine Station_ID2Calib(STATION_ID,Ant_ID,StatAnt_Calib, Calibrated)
    use constants, only : dp
    use StationMnemonics, only : Station_ID2Mnem
    !Use Chunk_AntInfo, only : Fine_STMnm, Fine_STDelay, Fine_AntDelay, SAI_AntDelay, Nr_AntDelay, Nr_StatDelay, AntInCal ! all these are generated in 'ReadCalib'
    Implicit none
    Integer, intent(in) :: STATION_ID
    Integer, intent(in) :: Ant_ID
    Logical, intent(out) :: Calibrated
    Real(dp), intent(out) :: StatAnt_Calib
    Character(len=6) :: Station_Mnem
    integer :: SAI
    Integer, save ::  i_fine
    Logical, save :: First=.true.
    integer :: i_Ant
     !
    If(first) then
      Call ReadCalib
      first=.false.
    endif   ! (first)
    !
    StatAnt_Calib=0
    Call Station_ID2Mnem(STATION_ID,Station_Mnem)
    Do i_Fine=1,Nr_StatDelay
        If(Fine_STMnm(i_Fine).eq.Station_Mnem) then
            StatAnt_Calib= Fine_STDelay(i_Fine)
            exit
        endif
    Enddo
    !
    If(Nr_AntDelay.eq.0) Stop 'There are no antenna calibrations specified !!!!!!!!!!!!!!'
    SAI=100*station_ID + Ant_ID
    i_Ant=MAXLOC(SAI_AntDelay(1:Nr_AntDelay), MASK = SAI_AntDelay(1:Nr_AntDelay) .eq. SAI, DIM=1)
    Calibrated=.false.
    !write(2,*) '!Station_ID2Calib', i_ant, SAI, i_Fine, Station_Mnem
    If(i_Ant.gt.0) Then  ! antenna was found
      !If(SAI_AntDelay(Ant(1)) .eq. SAI) Then
         StatAnt_Calib= StatAnt_Calib + Fine_AntDelay(i_Ant)
         Calibrated= AntInCal(i_Ant)
      !EndIf
    EndIf
    !write(2,*) '!Station_ID2Calib:',STATION_ID,Ant_ID,StatAnt_Calib, Calibrated, i_Ant,Station_Mnem, i_fine
    Return
End Subroutine Station_ID2Calib
!============================
Subroutine WriteCalibration(ZeroCalOffset) ! MergeFine
! Merge values of FineOffset with input from 'FineCalibrations.dat'
! Nr_AntDelay is increased with the antennas used in this calculation that were not on the read-in cal file
! MaximumAntennas is increased with the antennas used in this calculation that were not on the read-in cal file
   use constants, only : dp, sample, HeightCorrectIndxRef
   use DataConstants, only : Calibrations, RunMode  ! , DataFolder
   use DataConstants, only : FlashName, OutFileLabel
   use Chunk_AntInfo, only : Unique_StatID, Nr_UniqueStat, Unique_SAI, Tot_UniqueAnt, Nr_UniqueAnt, RefAnt, Ant_Stations
   use Chunk_AntInfo, only : ExcludedStat_max, SgnFlp_nr, PolFlp_nr, BadAnt_nr
   use Chunk_AntInfo, only : ExcludedStat, SaturatedSamplesMax, SignFlp_SAI, PolFlp_SAI, BadAnt_SAI
   !Use Chunk_AntInfo, only : Fine_STMnm, Fine_STDelay, Fine_AntDelay, SAI_AntDelay, Nr_AntDelay, Nr_StatDelay, AntInCal ! all these are generated in 'ReadCalib'
   Use Interferom_Pars, only : IntFer_ant
   use FitParams, only : Fit_AntOffset, Fit_TimeOffsetStat, Fit_TimeOffsetAnt, FitQual
   use unque,only : Double_IR_sort
   use StationMnemonics, only : Statn_ID2Mnem,Station_ID2Mnem
   Implicit none
   Real(dp), intent(out) :: ZeroCalOffset
   Real(dp) :: mean(1:Station_nrMax), UpDate_StDelay(1:Station_nrMax)=0.
   character(len=6) :: txt, Station_Mnem!, Fine_STMnm(Station_nrMax), LOFAR_STMnm(Station_nrMax)
   integer :: i, k, nxx, i_fine, SAI, n, i_ant, Stati(1:MaximumAntennas), NrMean(Station_nrMax)
   INTEGER :: DATE_T(8), i_stat, i_SAI, Nr_WriteStat, i_Ref, OutUnit, i_type, StatID
   Character(len=12) :: Date_mn
   Character(len=7) :: cmnt
   Character(len=70) :: CalibrationFileName, FMT
   Logical ::  AntInRun(1:MaximumAntennas)  !  =Antenna used in this calcularion
   !
   CALL DATE_AND_TIME (Values=DATE_T)
   WRITE(Date_mn,"(I4,I2.2,I2.2, I2.2,I2.2)") &
       DATE_T(1),DATE_T(2),DATE_T(3),(DATE_T(i),i=5,6)
   !
   AntInRun(:)=.false.
   Do i_stat=1,Nr_UniqueStat  ! set fine-offset to values obtained from present fit
     !Stat_ID=Unique_StatID(i_stat)
     !write(2,*) '!WriteCalibration, i-asai',i_stat,Unique_StatID(i_stat), Tot_UniqueAnt(i_stat-1)+1,Tot_UniqueAnt(i_stat)
     Do i_SAI=Tot_UniqueAnt(i_stat-1)+1,Tot_UniqueAnt(i_stat)      ! Get antenna number from the Unique_Antennas list
         SAI=Unique_SAI(i_SAI)  !     SAI=1000*station_ID + Ant_ID
         !  Unique_SAI(1:Nr_UniqueAnt) : the SAI for all antennas in this run
         !  SAI_AntDelay(i)   : the SAI for all antennas were calibration data were read, now supplemented with those in run
         i_Ant=MAXLOC(SAI_AntDelay(1:Nr_AntDelay), MASK = SAI_AntDelay(1:Nr_AntDelay) .eq. SAI, DIM=1) ! obtain position in the calibration data array
         !write(2,*) 'i_SAI=',i_SAI,SAI,SAI_AntDelay(Ant(1)),Ant(1)
         If(i_Ant .eq. 0) then  ! append this antenna to the end of the list
            Nr_AntDelay=Nr_AntDelay+1
            If(Nr_AntDelay .gt. MaximumAntennas) Then
               write(2,*) 'extent length antenna-delay list, ',MaximumAntennas,' is too small'
               stop 'WriteCalibration:extent'
            Endif
            i_Ant=Nr_AntDelay
            SAI_AntDelay(i_Ant) = SAI
            Fine_AntDelay(i_Ant)= 0.
         endif
         AntInRun(i_Ant)=.true.
         If(Fit_AntOffset) Fine_AntDelay(i_Ant)=Fine_AntDelay(i_Ant)+ Fit_TimeOffsetAnt(i_SAI)  ! Fit_TimeOffsetAnt is non-zero only when this option is active
         !write(2,*) 'UpDate_AntDelay',i_SAI,SAI,Nr_AntDelay, UpDate_AntDelay(Nr_AntDelay)
     Enddo
   Enddo
   !
   !Do further processing only for those antennas that were used in the present run
   i_ant=0
   Do i=1,Nr_AntDelay
      If(AntInRun(i)) Then
         i_ant=i_ant+1
         SAI_AntDelay(i_Ant)=SAI_AntDelay(i)
         Fine_AntDelay(i_Ant)=Fine_AntDelay(i)
      EndIf
   Enddo
   Nr_AntDelay=i_Ant
   ! Sort antennas
   Call Double_IR_sort(Nr_AntDelay,SAI_AntDelay,Fine_AntDelay)
   !
   ! set mean offset per station to zero
   !write(2,*) '!WriteCalibration, Nr_UniqueStat',Nr_UniqueStat, Nr_AntDelay
   !write(2,*) '!WriteCalibration, UpDate_AntDelay(1:30)',UpDate_AntDelay(1:30)
   Mean(:)=0.
   i_stat=-1 !NINT(SAI_AntDelay(1)/1000.)
   Nr_WriteStat=Nr_UniqueStat

   StatID=-1
   Mean(:)=0.
   NrMean(:)=0
   Do i=1,Nr_AntDelay  !  loop over all antenna delays, old & new
      k=SAI_AntDelay(i)/100
      If(k.ne.StatID) then  ! get i_stat for this station
         StatID=k          ! Store ID new station
         i_stat=MAXLOC(Unique_StatID(1:Nr_UniqueStat), MASK = Unique_StatID(1:Nr_WriteStat) .eq. StatID, DIM=1) ! FINDLOC(Unique_StatID(1:Nr_UniqueStat), i_stat) !
      EndIf
      Stati(i)=i_stat
      Mean(i_stat)=Mean(i_stat)+ Fine_AntDelay(i)
      NrMean(i_stat)=NrMean(i_stat) + 1
      !write(2,*) '!WriteCalibration, i',i, SAI_AntDelay(i), StatID, i_stat, Nr_WriteStat
   Enddo  !   i=1,Nr_AntDelay
   Do i_stat=1,Nr_UniqueStat
      If(NrMean(i_stat).gt.0) Mean(i_stat)=Mean(i_stat)/NrMean(i_stat)
   EndDo
   Do i=1,Nr_AntDelay  !  loop over all antenna delays, old & new
      If(AntInRun(i)) Fine_AntDelay(i)=Fine_AntDelay(i)-Mean(Stati(i))
   EndDo
   !
   If(RunMode.eq.2) Then
      i_type=0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      If(RefAnt(1,0,i_type).le.0) i_type=1
      StatID=Ant_Stations(RefAnt(1,0,i_type),1)
   Else
      StatID=Ant_Stations(IntFer_ant(1,1),1) ! If(allocated(IntFer_ant))
   EndIf
   i_fine=Nr_StatDelay  ! all stations read from calibration file; will be increased with the newly found ones
   !write(2,*) 'Nr_StatDelay:', Nr_StatDelay, Nr_WriteStat
   !Flush(unit=2)
   Do i_stat=1,Nr_UniqueStat    ! Write stationdelays
     Call Station_ID2Mnem(Unique_StatID(i_stat),Station_Mnem)
     If( Unique_StatID(i_stat).eq. StatID ) i_Ref=i_stat
     !write(2,*) 'mnem',k,Unique_StatID(k),Station_Mnem
     !flush(unit=2)
     !core=((RunMode.eq.1) .and. (Station_Mnem(1:2) .eq. 'CS'))  ! zero the calibration timings for the core stations
     nxx=0
     Do i=1,i_fine
         If(Station_Mnem .eq. Fine_STMnm(i)) then
             UpDate_StDelay(i)=Fit_TimeOffsetStat(i_stat) + Mean(i_stat)
             Fine_StDelay(i)  =Fine_StDelay(i) + UpDate_StDelay(i)
             !If(core) Fine_STDelay(i) = 0
             nxx=i
             exit
         endif
     enddo
     If(nxx.eq.0) then       ! New station that was not in the original file 'FineCalibrations.dat'
         i_fine=i_fine+1
         Fine_STMnm(i_fine) = Station_Mnem
         UpDate_StDelay(i_fine)=Fit_TimeOffsetStat(i_stat)+ Mean(i_stat)
         Fine_STDelay(i_fine)  = UpDate_StDelay(i_fine)
         !If(core) Fine_STDelay(i_fine) = 0
         nxx=i_fine
         !write(2,*) 'i_fine',i_fine,UpDate_StDelay(i_fine), Station_Mnem, k
     endif
   Enddo
   !
   ! set shift for a reference station = 0
   ZeroCalOffset=Fine_STDelay(i_ref)
   ZeroCalOffset=0.0  ! dangerous to change this to non-zero since it will shift effective peak positions
   Nr_StatDelay=i_fine  ! all stations, ECHBG, increased with the newly found ones
   !write(2,*) 'Calibration constant zeroed for ref-station:', Fine_STMnm(k_ref)

   !  ============== start writing
   ! get unique filename
   cmnt="Hilbert"
   CalibrationFileName="Hil"
   If(RunMode.eq.7) Then
      cmnt="FldCal"
      CalibrationFileName="Fld"
   EndIf
   CalibrationFileName=TRIM(CalibrationFileName)//TRIM(FlashName)//TRIM(OutFileLabel)//'-'//Date_mn//'.cal'
   Open(unit=19,STATUS='unknown',ACTION='write', FILE = 'Book/'//TRIM(CalibrationFileName), IOSTAT=nxx)  !  trim(DataFolder)//
   write(2,"(A,1x,L,F6.2,'  ')") '  Calibrations="'//TRIM(CalibrationFileName)//'" ! '//cmnt,HeightCorrectIndxRef, FitQual  ! trim(DataFolder)//
   Write(19,"(A4,1x,A,A,L)") '!04 ',Date_mn,' HeightCorrectIndxRef= ',HeightCorrectIndxRef
   n=1
   Do i=1,Nr_AntDelay
      ! Check wether this antenna is part of the present analysis, i.e. in Unique_SAI(1:Nr_UniqueAnt)
      !n=COUNT(Unique_SAI(1:Nr_UniqueAnt).eq. SAI_AntDelay(i))
      If(abs(Fine_AntDelay(i)).lt. 99.) Then
         write(19,"(I8, F10.4, I3 )") SAI_AntDelay(i),Fine_AntDelay(i), n
      Else
         write(19,*) SAI_AntDelay(i),Fine_AntDelay(i), n
         write(2,*) '******** Ridiculosly large antenna timing calibration:',SAI_AntDelay(i),Fine_AntDelay(i)
      EndIf
   Enddo
   !
   write(19,*) '============ =========== =============== ======= all in units of LBA-samples'
   !
   !write(2,*) 'Nr_StatDelay_2:', Nr_StatDelay, Nr_WriteStat
   write(2,"(' station',1x,'NewCalibrations; Updated with [samples]; overall calibration shift=',F8.1)") ZeroCalOffset
   Do i=1,i_fine
      n=0
      Do k=1,Nr_UniqueStat
         If(Statn_ID2Mnem(Unique_StatID(k)) .ne. Fine_STMnm(i) ) cycle
         n=1
         exit
      Enddo
      !n=COUNT(Unique_SAI(1:Nr_UniqueAnt).eq. Fine_STMnm(i))
      ! The offset in the station delays will shift the peak location in the time-traces
      write(2,"(1x,A6,2F13.3,i5)") Fine_STMnm(i), Fine_STDelay(i)-ZeroCalOffset, UpDate_StDelay(i),n
      write(19,"(1x,A6,F14.4,i5)") Fine_STMnm(i), Fine_STDelay(i)-ZeroCalOffset, n
   enddo
   !
   write(19,*) '!========================================'
   n=0
   Do i=1,ExcludedStat_max
      If(TRIM(ExcludedStat(i)).eq.'') exit
      n=n+1
   enddo
   !If(n.gt.0) Then
   !   Date_mn=15H('"',A,'"')
   !   write(FMT,"(A,I2,A,')'   )") n,Date_mn
   !   write(2,*) 'FMT=',FMT
   !   Flush(unit=2)
   !   write(19, FMT) 'ExcludedStat=', ExcludedStat(1:n)
   !EndIf
   write(19,"(  50A)") 'ExcludedStat=', ('"', ExcludedStat(i), '",', i=1,n)
   write(19,"( A,40(I7,',') )") 'BadAnt_SAI= ', BadAnt_SAI(1:BadAnt_nr)
   write(19,"( A,40(I7,',') )") 'SignFlp_SAI=', SignFlp_SAI(1:SgnFlp_nr)
   write(19,"( A,40(I7,',') )") 'PolFlp_SAI= ', PolFlp_SAI(1:PolFlp_nr)
   OutUnit=19
   Call PrntNewSources(ZeroCalOffset, OutUnit)
   Close(unit=19)
   !
   Return
End Subroutine WriteCalibration
!----------------------------------
End Module Calibration
!=================================
