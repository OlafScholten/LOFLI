Module Calibration
    use constants, only : dp
    use DataConstants, only : Ant_nrMax, Station_nrMax
    Character(len=5), save :: Fine_STMnm(1:Station_nrMax)          ! station mnemonics for which Ever Calibrations Have Been Generated (ECHBG)
    Real(dp), save :: Fine_STDelay(1:Station_nrMax)! time calibration data for stations
    Real(dp), save :: Fine_AntDelay(1:Ant_nrMax)! time calibration data for antennas
    Integer, save :: SAI_AntDelay(1:Ant_nrMax)  ! SAI of antennas for which ECHBG
    Integer, save :: CalibrDelay(1:Ant_nrMax)  ! =1 if this antenna was used when writing cal-file
    Integer, save :: Nr_AntDelay, Nr_StatDelay  ! Nr of antennas or stations for which ECHBG
    Integer, save :: StationInCal(Station_nrMax)  ! =1 if this station was used when writing cal-file
Contains
! -------------------------------------------------
Subroutine ReadCalib()
    use constants, only : dp,sample, HeightCorrectIndxRef
    use DataConstants, only : Station_nrMax, DataFolder, Ant_nrMax, Calibrations, RunMode ! , Diagnostics
    !Use Chunk_AntInfo, only : Fine_STMnm, Fine_STDelay, Fine_AntDelay, SAI_AntDelay, Nr_AntDelay, Nr_StatDelay, CalibrDelay ! all these are output from this routine
    Implicit none
    Character(len=5) :: txt !, Station_Mnem
    Real(dp) :: Delay
    Character(len=5), save :: LOFAR_STMnm(Station_nrMax) ! Archaic
    Real(dp), save :: LOFAR_STDelay(Station_nrMax) ! Archaic
    integer :: nxx, i, SAI
    Integer :: i_LOFAR, i_fine, n
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
      Enddo
   endif
   close(unit=9)
   !
   Nr_AntDelay=0
   Old=(trim(Calibrations).eq.'')
   If(Old) then
     Open(unit=9,STATUS='old',ACTION='read', FILE = 'Book/AntCalibrations.dat', IOSTAT=nxx) ! trim(DataFolder)//
     If(nxx.ne.0) then
         write(2,*) 'problems with file=','AntCalibrations.dat'
         stop 'ReadCalib:old'
     EndIf
   Else
     Open(unit=9,STATUS='old',ACTION='read', FILE = 'Book/'//trim(Calibrations), IOSTAT=nxx)
     If(nxx.ne.0) then
         write(2,*) 'problems with file=','Book/'//trim(Calibrations),', probably does not exist'
         stop 'ReadCalib'
     EndIf
   Endif
   !
   read(9,"(A60)") tst  ! Check for a new-style calibration file
   read(tst,*,IOSTAT=nxx) Mrk, date, TxtIndx, IndxRef  ! Check for Height correction Index Ref
   !write(2,*) 'tst:',tst
   !write(2,*) 'From ReadCalib:nxx=',nxx,Mrk,i,' ', date,' ', TxtIndx,' ', IndxRef
   If(nxx.eq.0) Then
      FMT2022=.true.
      i_LOFAR=0  ! the data of  'Book/StationCalibrations.dat' (Archaic) have already been included in the station delays
      write(2,*) 'From ReadCalib:',Mrk,' ', date,' ', TxtIndx,' ', IndxRef
   Else
      IndxRef=.false.  ! In the old days no correction for the height dependence of the index of refraction was used
      read(tst,*,IOSTAT=nxx) SAI, Delay,n  ! Check for a new-style calibration file
      If(nxx.eq.0) Then
         FMT2022=.true.
         i_LOFAR=0  ! the data of  'Book/StationCalibrations.dat' (Archaic) have already been included in the station delays
      Else
         FMT2022=.false.
         If(i_LOFAR.eq.0) Then
            write(2,*) 'In "ReadCalib": something went wrong reading the file StationCalibrations.dat'
            Stop 'Error in ReadCalib'
         EndIf
      EndIf
      rewind(unit=9)
   EndIf  ! file 9 is now positioned to first line with antenna data
   If(runmode.ne.2 .and. runmode.ne.7) Then  ! do not overwrite parameter when set on input
      HeightCorrectIndxRef= IndxRef
      write(2,*) 'runmode, HeightCorrectIndxRef:',runmode, HeightCorrectIndxRef
   EndIf
   !
   n=-1
   Do       ! Read single antenna calibrations
      If(FMT2022) Then
         read(9,*,IOSTAT=nxx) SAI, Delay,n
      Else
         read(9,*,IOSTAT=nxx) SAI, Delay
      EndIf
      if(nxx.ne.0) exit
      Nr_AntDelay=Nr_AntDelay+1
      SAI_AntDelay(Nr_AntDelay) = SAI
      Fine_AntDelay(Nr_AntDelay)= Delay ! already in [samples]
      CalibrDelay(Nr_AntDelay)=n
   enddo
   !
   nxx=0
   If(Old) Then
      Close(unit=9)
      Open(unit=9,STATUS='old',ACTION='read', FILE = 'Book/FineCalibrations.dat', IOSTAT=nxx) ! trim(DataFolder)//
      If(nxx.ne.0) then
         write(2,*) 'problems with file=','FineCalibrations.dat'
         stop 'ReadCalib, Reading problem old style'
      EndIf
   EndIf
   !
   Nr_StatDelay=0
   StationInCal(:)=0
   Do  ! get station calibrations
      read(9,"(A60)",IOSTAT=nxx) tst  ! Check for a new-style calibration file
       if(nxx.ne.0) exit
       read(tst,*,IOSTAT=nxx) txt, Delay, n
       !write(2,*) nxx,txt, Delay, n,tst
       if(nxx.ne.0) n=0
           Nr_StatDelay=Nr_StatDelay+1
           Fine_STMnm(Nr_StatDelay) = txt
           Fine_STDelay(Nr_StatDelay)= Delay ! already in samples
           StationInCal(Nr_StatDelay)=n
   enddo
   close(unit=9)
   !
   If(.not.FMT2022) Then  ! Archaic delays need to be merged
      Do i_fine=1,Nr_StatDelay    ! Write stationdelays
         Do i=1,i_LOFAR  ! include old station delays from a previous age
            If(LOFAR_STMnm(i).eq.Fine_STMnm(i_fine)) then
               Fine_STDelay(i_fine)  =Fine_STDelay(i_fine) + LOFAR_STDelay(i)
               exit
            EndIf
         EndDo
      EndDo
   EndIf
   write(2,*) 'Nr_AntDelay=',Nr_AntDelay
   Return
End Subroutine ReadCalib
!=================================
Subroutine Station_ID2Calib(STATION_ID,Ant_ID,StatAnt_Calib, Calibrated)
    use constants, only : dp
    use StationMnemonics, only : Station_ID2Mnem
    !Use Chunk_AntInfo, only : Fine_STMnm, Fine_STDelay, Fine_AntDelay, SAI_AntDelay, Nr_AntDelay, Nr_StatDelay, CalibrDelay ! all these are generated in 'ReadCalib'
    Implicit none
    Integer, intent(in) :: STATION_ID
    Integer, intent(in) :: Ant_ID
    Integer, intent(out) :: Calibrated
    Real(dp), intent(out) :: StatAnt_Calib
    Character(len=5) :: Station_Mnem
    integer :: SAI
    Integer, save ::  i_fine
    Logical, save :: First=.true.
    integer :: Ant(1)
     !
    If(first) then
      Call ReadCalib
      first=.false.
    endif   ! (first)
    !
    StatAnt_Calib=0.
    Call Station_ID2Mnem(STATION_ID,Station_Mnem)
    Do i_Fine=1,Nr_StatDelay
        If(Fine_STMnm(i_Fine).eq.Station_Mnem) then
            StatAnt_Calib= Fine_STDelay(i_Fine)
            exit
        endif
    Enddo
    !
    If(Nr_AntDelay.eq.0) Stop 'There are no antenna calibrations specified !!!!!!!!!!!!!!'
    SAI=1000*station_ID + Ant_ID
    Ant=MAXLOC(SAI_AntDelay(1:Nr_AntDelay), MASK = SAI_AntDelay(1:Nr_AntDelay) .eq. SAI)
    If(Ant(1).gt.0) Then  ! antenna was found
      !If(SAI_AntDelay(Ant(1)) .eq. SAI) Then
         StatAnt_Calib= StatAnt_Calib + Fine_AntDelay(Ant(1))
         Calibrated= CalibrDelay(Ant(1))
      !EndIf
    EndIf
    Return
End Subroutine Station_ID2Calib
!============================
Subroutine WriteCalibration(ZeroCalOffset) ! MergeFine
! Merge values of FineOffset with input from 'FineCalibrations.dat'
   use constants, only : dp, sample, HeightCorrectIndxRef
   use DataConstants, only : Station_nrMax, Ant_nrMax, Calibrations, RunMode  ! , DataFolder
   use DataConstants, only : FlashName, OutFileLabel
   use Chunk_AntInfo, only : Unique_StatID, Nr_UniqueStat, Unique_SAI, Tot_UniqueAnt, Nr_UniqueAnt, RefAnt, Ant_Stations
   Use Interferom_Pars, only : IntFer_ant
   !Use Chunk_AntInfo, only : Fine_STMnm, Fine_STDelay, Fine_AntDelay, SAI_AntDelay, Nr_AntDelay, Nr_StatDelay, CalibrDelay ! all these are generated in 'ReadCalib'
   use FitParams, only : Fit_AntOffset, Fit_TimeOffsetStat, Fit_TimeOffsetAnt, FitQual
   use unque,only : Double_IR_sort
   use StationMnemonics, only : Statn_ID2Mnem,Station_ID2Mnem
   Implicit none
   Real(dp), intent(out) :: ZeroCalOffset
   Real(dp) :: Delay, mean(Station_nrMax), UpDate_STDelay(1:Ant_nrMax)=0.
   character(len=5) :: txt, Station_Mnem!, Fine_STMnm(Station_nrMax), LOFAR_STMnm(Station_nrMax)
   integer :: i, k, nxx, i_fine, SAI, n, i_unq
   INTEGER :: DATE_T(8), Ant(1), i_stat, i_SAI, Nr_WriteStat, k_Ref
   Character(len=12) :: Date_mn
   Character(len=7) :: cmnt
   Character(len=70) :: CalibrationFileName
   !
   CALL DATE_AND_TIME (Values=DATE_T)
   WRITE(Date_mn,"(I4,I2.2,I2.2, I2.2,I2.2)") &
       DATE_T(1),DATE_T(2),DATE_T(3),(DATE_T(i),i=5,6)
   !
   i_fine=Nr_AntDelay ! all antennas, ECHBG
   UpDate_STDelay(:)=0.  ! merge with station/antenna delays obtained from fit
   If(Fit_AntOffset) then
      Do i_stat=1,Nr_UniqueStat  ! set fine-offset to values obtained from present fit
        !Stat_ID=Unique_StatID(i_stat)
        !write(*,*) 'i-asai',Tot_UniqueAnt(i_stat-1)+1,Tot_UniqueAnt(i_stat)
        Do i_SAI=Tot_UniqueAnt(i_stat-1)+1,Tot_UniqueAnt(i_stat)      ! Get antenna number from the Unique_Antennas list
            SAI=Unique_SAI(i_SAI)  !     SAI=1000*station_ID + Ant_ID
            Ant=MAXLOC(SAI_AntDelay(1:Nr_AntDelay), MASK = SAI_AntDelay(1:Nr_AntDelay) .eq. SAI) ! obtain position in the calibration data array
            !write(2,*) 'i_SAI=',i_SAI,SAI,SAI_AntDelay(Ant(1)),Ant(1)
            If(Ant(1) .ne. 0) then  ! this antenna was already in the list
               UpDate_STDelay(Ant(1))= Fit_TimeOffsetAnt(i_SAI)
            else  ! append this antenna to the end of the list
               Nr_AntDelay=Nr_AntDelay+1
               If(Nr_AntDelay .gt. Ant_nrMax) Then
                  write(2,*) 'extent length antenna-delay list, ',Ant_nrMax,' is too small'
                  stop 'WriteCalibration:extent'
               Endif
               SAI_AntDelay(Nr_AntDelay) = SAI
               UpDate_STDelay(Nr_AntDelay)= Fit_TimeOffsetAnt(i_SAI)
               Fine_AntDelay(Nr_AntDelay)= 0.
            endif
            !write(2,*) 'UpDate_STDelay',i_SAI,SAI,Nr_AntDelay, UpDate_STDelay(Nr_AntDelay)
        Enddo
      Enddo
   Endif    ! (Fit_AntOffset)
   ! Sort antennas and update the calibration
   Fine_AntDelay(1:Nr_AntDelay)=Fine_AntDelay(1:Nr_AntDelay) + UpDate_STDelay(1:Nr_AntDelay)
   Call Double_IR_sort(Nr_AntDelay,SAI_AntDelay,Fine_AntDelay)
   !
   ! set mean offset per station to zero
   Mean(:)=0.
   i_stat=-1 !NINT(SAI_AntDelay(1)/1000.)
   i_unq=1
   n=1
   Nr_WriteStat=Nr_UniqueStat
   Do i=1,Nr_AntDelay  !  loop over all antenna delays, old & new
      k=NINT(SAI_AntDelay(i)/1000.)
      If(k.eq.i_stat) then  ! update running sum for this station
         Mean(i_unq)=Mean(i_unq)+ Fine_AntDelay(i)
         n=n+1
      Else
         Mean(i_unq)=Mean(i_unq)/n ! Calculate mean from running sum for previous station
         i_stat=k          ! Store ID new station
         If(i-n.ge.1) then
            Do k=i-n,i-1      ! Set mean antenna delay to zero for previous stattion
             Fine_AntDelay(k)=Fine_AntDelay(k)-Mean(i_unq)
            Enddo
         endif
         Ant=MAXLOC(Unique_StatID(1:Nr_WriteStat), MASK = Unique_StatID(1:Nr_WriteStat) .eq. i_stat) ! FINDLOC(Unique_StatID(1:Nr_UniqueStat), i_stat) !
         i_unq=Ant(1)
         If(i_unq.eq.0) then  ! If station not in the list, add it
            Nr_WriteStat=Nr_WriteStat+1
            If(Nr_WriteStat.gt.Station_nrMax) then
               write(2,*) 'nr of unique stations exceeded for',i_stat
               stop 'WriteCalibration:exceed'
            Endif
            Unique_StatID(Nr_WriteStat)= i_stat ! add to the unique station list
            Fit_TimeOffsetStat(Nr_WriteStat)=0.
            i_unq=Nr_WriteStat
         Endif
         n=1
         Mean(i_unq)=Fine_AntDelay(i)  ! Start running sum for this station
      Endif
   Enddo  !   i=1,Nr_AntDelay
   Mean(i_unq)=Mean(i_unq)/n  ! for the last one
   Do k=Nr_AntDelay-n,Nr_AntDelay
    Fine_AntDelay(k)=Fine_AntDelay(k)-Mean(i_unq)
   Enddo
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
   Write(19,"(A4,1x,A,A,L)") '!01 ',Date_mn,' HeightCorrectIndxRef= ',HeightCorrectIndxRef
   Do i=1,Nr_AntDelay
      !k=NINT(SAI_AntDelay(i)/1000.)
      !Call Station_ID2Mnem(k,Station_Mnem)
      ! Check wether this antenna is part of the present analysis, i.e. in Unique_SAI(1:Nr_UniqueAnt)
      n=COUNT(Unique_SAI(1:Nr_UniqueAnt).eq. SAI_AntDelay(i))
      !Fine_AntDelay(i)= Fine_AntDelay(i) +UpDate_STDelay(i)
      write(19,*) SAI_AntDelay(i),Fine_AntDelay(i), n
   Enddo
   write(19,*) '============ =========== =============== ======= all in units of samples'
   If(RunMode.eq.2) Then
      n=Ant_Stations(RefAnt(1,0),1)
   Else
      n=Ant_Stations(IntFer_ant(1,1),1) ! If(allocated(IntFer_ant))
   EndIf
   !
   i_fine=Nr_StatDelay  ! all stations, ECHBG; will be increased with the newly found ones
   Do k=1,Nr_WriteStat    ! Write stationdelays
     If(Unique_StatID(k).le. 0) exit  ! should not happen; The first are the stations
     Call Station_ID2Mnem(Unique_StatID(k),Station_Mnem)
     !write(2,*) 'mnem',k,Unique_StatID(k),Station_Mnem
     !core=((RunMode.eq.1) .and. (Station_Mnem(1:2) .eq. 'CS'))  ! zero the calibration timings for the core stations
     If( Unique_StatID(k).eq. n ) k_Ref=k
     nxx=0
     Do i=1,i_fine
         If(Station_Mnem .eq. Fine_STMnm(i)) then
             UpDate_STDelay(i)=Fit_TimeOffsetStat(k) + Mean(k)
             Fine_STDelay(i)  =Fine_STDelay(i) + UpDate_STDelay(i)
             !If(core) Fine_STDelay(i) = 0
             nxx=i
             exit
         endif
     enddo
     If(nxx.eq.0) then       ! New station that was not in the original file 'FineCalibrations.dat'
         i_fine=i_fine+1
         Fine_STMnm(i_fine) = Station_Mnem
         UpDate_STDelay(i_fine)=Fit_TimeOffsetStat(k)+ Mean(k)
         Fine_STDelay(i_fine)  = UpDate_STDelay(i_fine)
         !If(core) Fine_STDelay(i_fine) = 0
         nxx=i_fine
         !write(2,*) 'i_fine',i_fine,UpDate_STDelay(i_fine)
     endif
   Enddo
   ZeroCalOffset=Fine_STDelay(k_ref)
   ZeroCalOffset=0.0  ! dangerous to change this to non-zero since it will shift effective peak positions
   Nr_StatDelay=i_fine  ! all stations, ECHBG, increased with the newly found ones
   !write(2,*) 'Calibration constant zeroed for ref-station:', Fine_STMnm(k_ref)
   !flush(unit=2)
   ! RefAnt(i_chunk,i_eo)
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
      write(2,"(1x,A5,2F13.3,i5)") Fine_STMnm(i), Fine_STDelay(i)-ZeroCalOffset, UpDate_STDelay(i),n
      write(19,"(1x,A5,F14.4,i5)") Fine_STMnm(i), Fine_STDelay(i)-ZeroCalOffset, n
   enddo
   Close(unit=19)
   !
   Return
End Subroutine WriteCalibration
!----------------------------------
End Module Calibration
!=================================
