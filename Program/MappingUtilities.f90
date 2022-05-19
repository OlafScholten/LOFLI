!=================================
Module unque
contains
Subroutine unique(vec,vec_unique)
! copied from http://degenerateconic.com/unique/
! Return only the unique values from vec.

implicit none

integer,dimension(:),intent(in) :: vec
integer,dimension(:),allocatable,intent(out) :: vec_unique

integer :: i,num
logical,dimension(size(vec)) :: mask

mask = .false.

do i=1,size(vec)

    !count the number of occurrences of this element:
    num = count( vec(i)==vec )

    if (num==1) then
        !there is only one, flag it:
        mask(i) = .true.
    else
        !flag this value only if it hasn't already been flagged:
        if (.not. any(vec(i)==vec .and. mask) ) mask(i) = .true.
    End if

End do

!return only flagged elements:
allocate( vec_unique(count(mask)) )
vec_unique = pack( vec, mask )

!if you also need it sorted, then do so.
call Selection_sort(vec_unique)

End Subroutine unique
!-----------------------
  Subroutine Selection_sort(a)
  ! From:  https://rosettacode.org/wiki/Sorting_algorithms/Selection_sort#Fortran
    INTEGER, INTENT(IN OUT) :: a(:)
    INTEGER :: i, minIndex, temp

    DO i = 1, SIZE(a)-1
       minIndex = MINLOC(a(i:), 1) + i - 1
       IF (a(i) > a(minIndex)) THEN
          temp = a(i)
          a(i) = a(minIndex)
          a(minIndex) = temp
       End IF
    End DO
  End Subroutine Selection_sort
!-----------------------
Subroutine Double_sort(a)
! Adepted From:  https://rosettacode.org/wiki/Sorting_algorithms/Selection_sort#Fortran
!  Double sorts a double array a according to the first element
    INTEGER, INTENT(IN OUT) :: a(:,:)
    INTEGER :: i, minIndex, temp(SIZE(a,2))
    !
    DO i = 1, SIZE(a,1)-1
       minIndex = MINLOC(a(i:,1), 1) + i - 1
       IF (a(i,1) > a(minIndex,1)) THEN ! sort on first column
          temp(:) = a(i,:)
          a(i,:)= a(minIndex,:)
          a(minIndex,:) = temp(:)
       End IF
    End DO
End Subroutine Double_sort
!-----------------------
Pure Subroutine Double_IR_sort(N,a,R)
! Adepted From:  https://rosettacode.org/wiki/Sorting_algorithms/Selection_sort#Fortran
!  sorts integer array a according to the first dimension
! sort also the additional real array
    use constants, only : dp
    implicit none
    INTEGER, INTENT(IN) :: N  ! length of arrays a and R that should be sorted, actual size may be larger
    INTEGER, INTENT(IN OUT) :: a(:)
    Real(dp), INTENT(IN OUT) :: R(:)
    INTEGER :: i, minIndex, temp, Rtemp

    DO i = 1, N-1
       minIndex = MINLOC(a(i:N), 1) + i - 1
       IF (a(i) > a(minIndex)) THEN
          temp = a(i)
          a(i) = a(minIndex)
          a(minIndex) = temp
          Rtemp = R(i)
          R(i) = R(minIndex)
          R(minIndex) = Rtemp
       End IF
    End DO
End Subroutine Double_IR_sort
!-----------------------
Pure Subroutine Double_RI_sort(N,a,R)
! Adepted From:  https://rosettacode.org/wiki/Sorting_algorithms/Selection_sort#Fortran
!  sorts integer array a (min first) and reorders the additional real array in the same way
    use constants, only : dp
    implicit none
    INTEGER, INTENT(IN) :: N  ! length of arrays a and R that should be sorted, actual size may be larger
    Real(dp), INTENT(IN OUT) :: a(:)
    Integer, INTENT(IN OUT) :: R(:)
    INTEGER :: i, minIndex
    Real(dp) :: temp
    Integer :: Rtemp

    DO i = 1, N-1
       minIndex = MINLOC(a(i:N), 1) + i - 1
       IF (a(i) > a(minIndex)) THEN
          temp = a(i)
          a(i) = a(minIndex)
          a(minIndex) = temp
          Rtemp = R(i)
          R(i) = R(minIndex)
          R(minIndex) = Rtemp
       End IF
    End DO
End Subroutine Double_RI_sort
!-------------------------------------
Subroutine sort(n, a)
! From  https://rosettacode.org/wiki/Sorting_algorithms/Insertion_sort#Fortran
    implicit none
    integer :: n, i, j
    real :: a(n), x

    do i = 2, n
        x = a(i)
        j = i - 1
        do while (j >= 1)
            if (a(j) <= x) exit
            a(j + 1) = a(j)
            j = j - 1
        End do
        a(j + 1) = x
    End do
End Subroutine sort
! -----------------------------------
End Module unque
!=================================
Module StationMnemonics
    Use Chunk_AntInfo, only : Station_name, Station_number
    Integer, parameter :: StationNrMax=39  ! smaller or equal to Station_nrMax
    !Character(len=5), parameter :: Station_name(1:StationNrMax)=(/&
    !    "CS001","CS002","CS003","CS004","CS005","CS006","CS007","CS011","CS013","CS017",&
    !    "CS021","CS024","CS026","CS030","CS031","CS032","CS101","CS103","CS201","CS301",&
    !    "CS401","CS501","RS106","RS205","RS208","RS210","RS305","RS306","RS307","RS310",&
    !    "RS406","RS407","RS409","RS503","RS508","CS302","CS028","RS509","RS   "  /)
    !Integer, parameter :: Station_number(1:StationNrMax)=&
    !    (/   1 ,     2 ,     3 ,     4 ,     5 ,     6 ,     7 ,    11 ,    13 ,    17 ,&
    !        21 ,    24 ,    26 ,    30 ,    31 ,    32 ,   101 ,   103 ,   121 ,   141 ,&
    !       161 ,   181 ,   106 ,   125 ,   128 ,   130 ,   145 ,   146 ,   147 ,   150 ,&
    !       166 ,   167 ,   169 ,   183 ,   188 ,   142 ,    28 ,   189 ,     0      /)
    !  name number
    !   00x  0n
    !   01x  1n
    !   02x  2n
    !   03x  3n
    !   10x 10n
    !   20x 12n
    !   21x 13n
    !   30x 14n
    !   31x 15n
    !   40x 16n
    !   50x 18n
Contains
Subroutine Station_ID2Mnem(STATION_ID,Station_Mnem)
    Implicit none
    Integer, intent(in) :: STATION_ID
    Character(len=5), intent(out) :: Station_Mnem
    Integer :: i
    Station_Mnem="XSNNN"
    Do i=1,StationNrMax
        If(STATION_ID.eq.Station_number(i)) then
            Station_Mnem=Station_name(i)
            exit
        endif
    enddo
    !Write(2,*) 'Station_Mnem',Station_Mnem,i,STATION_ID
    return
End Subroutine Station_ID2Mnem
!
Character(len=5) Function Statn_ID2Mnem(STATION_ID)
    Implicit none
    Integer, intent(in) :: STATION_ID
    Character(len=5) :: Station_Mnem
    Call Station_ID2Mnem(STATION_ID,Station_Mnem)
    Statn_ID2Mnem=Station_Mnem
End Function Statn_ID2Mnem
!
Subroutine Station_Mnem2ID(Station_Mnem,STATION_ID)
    Implicit none
    Integer, intent(out) :: STATION_ID
    Character(len=5), intent(in) :: Station_Mnem
    Integer :: i
    Station_ID=0
    Do i=1,StationNrMax
        If(Station_Mnem(3:5) .eq. Station_name(i)(3:5)) then
            STATION_ID=Station_number(i)
            exit
        endif
    enddo
    If(Station_ID.eq.0) then
        Write(2,*) '*******Station_Mnem not found',Station_Mnem
    Endif
    return
End Subroutine Station_Mnem2ID
!
Integer Function Statn_Mnem2ID(Station_Mnem)
    Implicit none
    Character(len=5), intent(in) :: Station_Mnem
    Integer :: STATION_ID
    Call Station_Mnem2ID(Station_Mnem,STATION_ID)
    Statn_Mnem2ID=Station_ID
End Function Statn_Mnem2ID
!
End Module StationMnemonics
!=================================
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
   EndIf
   If(runmode.ne.2 .and. runmode.ne.7) Then
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
Subroutine WriteCalibration ! MergeFine
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
   Nr_StatDelay=i_fine  ! all stations, ECHBG, increased with the newly found ones
   !write(2,*) 'Calibration constant zeroed for ref-station:', Fine_STMnm(k_ref)
   !flush(unit=2)
   ! RefAnt(i_chunk,i_eo)
   write(2,"(' station',1x,'NewCalibrations; Updated with [samples]')")
   Do i=1,i_fine
      n=0
      Do k=1,Nr_UniqueStat
         If(Statn_ID2Mnem(Unique_StatID(k)) .ne. Fine_STMnm(i) ) cycle
         n=1
         exit
      Enddo
      !n=COUNT(Unique_SAI(1:Nr_UniqueAnt).eq. Fine_STMnm(i))
      write(2,"(1x,A5,2F13.3,i5)") Fine_STMnm(i), Fine_STDelay(i)-Fine_STDelay(k_ref), UpDate_STDelay(i),n
      write(19,"(1x,A5,F14.4,i5)") Fine_STMnm(i), Fine_STDelay(i)-Fine_STDelay(k_ref), n
   enddo
   Close(unit=19)
   !
   Return
End Subroutine WriteCalibration
!----------------------------------
End Module Calibration
!=================================

Subroutine ITRF2LOFARConstruct()
    use constants, only : dp,pi,RotMat,center_CS002
    implicit none
!    real(dp), intent(out) :: RotMat(3,3) = reshape( (/ 0.8056597 ,  0.        ,  0.59237863 , &
!                              -0.05134149, -0.99623707,  0.06982658 , &
!                              -0.59014956,  0.08667007,  0.80262806 /), (/ 3,3/))
!'''
! 52.99006128232752 deg North and 6.350944049999953 deg East
!rotation_matrix = np.array([[-0., -0.8056597, 0.59237863],
!       [0.99623707, 0.05134149, 0.06982658],
!       [-0.08667007, 0.59014956, 0.80262806]])
!
!# Center of station CS002 in ITRF coordinates, substract this before rotating
!    real(dp) :: center_CS002(3)= (/ 3826577.5, 461021.3125, 5064893.0 /)
    real(dp) :: x,a,b,Rotmatphi(3,3),Rotmatth(3,3),LOFAR(3)
    real(dp) :: theta=(52.915141-90.d0)*pi/180., phi=6.8698134123520900*pi/180. ! corrects for omegaxomegaxr-->0.19878 deg
!    real(dp) :: theta=(52.729828630863530-90.d0)*pi/180., phi=6.8698134123520900*pi/180. ! CS002 in Z-direction

    integer :: i,j,k
!    RotMat = (/ 0.8056597 ,  0. , 0.59237863 , -0.05134149, -0.99623707,  0.06982658 ,  -0.59014956,  0.08667007,  0.80262806 /)
    Rotmatphi=0.
    RotMatphi(1,1)=cos(phi)  ; RotMatphi(1,2)=sin(phi)
    RotMatphi(2,1)=-sin(phi) ; RotMatphi(2,2)=cos(phi)
    RotMatphi(3,3)=1
    RotMatTh=0.
    RotMatTh(2,2)=1.
    RotMatTh(1,1)=-cos(theta)  ; RotMatTh(1,3)=-sin(theta)
    RotMatTh(3,1)=-sin(theta) ; RotMatTh(3,3)=cos(theta)
    Do i=1,3
        Do j=1,3
            x=0.
            Do k=1,3
                x=x+RotMatTh(i,k)*RotMatPhi(k,j)
            enddo
            RotMat(i,j)=x
        enddo
    enddo

    a=0.
    Do i=1,3
        a=a + center_CS002(i)**2
        x=0.
        Do j=1,3
            x=x + RotMat(i,j) * center_CS002(j)
        enddo
        LOFAR(i)=x
    enddo
!    write(2,*) 'cs002=',LOFAR
    b=center_CS002(1)**2 + center_CS002(2)**2
    a=sqrt(a)
    b=sqrt(b)
!    write(2,*) 'phi=',b,(center_CS002(2)/b),asin(center_CS002(2)/b) *180./pi
!    write(2,*) 'theta=',a,(center_CS002(3)/a),asin(center_CS002(3)/a) *180./pi,asin(b/a)*180./pi
!    write(2,*) 'RotMat=',RotMat
    return
End Subroutine ITRF2LOFARConstruct
!-------------------------------------
!-----------------------------------------
Subroutine ITRF2LOFAR(ITRF,LOFAR)
    use constants, only : dp,pi,RotMat,center_CS002
    implicit none
    real(dp), intent(in) :: ITRF(3)
    real(dp), intent(out) :: LOFAR(3)
!    real(dp) :: RotMat(3,3) = reshape( (/ 0.8056597 ,  0.        ,  0.59237863 , &
!                              -0.05134149, -0.99623707,  0.06982658 , &
!                              -0.59014956,  0.08667007,  0.80262806 /), (/ 3,3/))
!'''
! 52.99006128232752 deg North and 6.350944049999953 deg East
!rotation_matrix = np.array([[-0., -0.8056597, 0.59237863],
!       [0.99623707, 0.05134149, 0.06982658],
!       [-0.08667007, 0.59014956, 0.80262806]])
!
!# Center of station CS002 in ITRF coordinates, substract this before rotating
!    real(dp) :: center_CS002(3)= (/ 3826577.5, 461021.3125, 5064893.0 /)
    real(dp) :: x
    integer :: i,j
    Do i=1,3
        x=0.
        Do j=1,3
            x=x + RotMat(i,j) * (ITRF(j)-center_CS002(j))
        enddo
        LOFAR(i)=x  ! 1=North, 2=East, 3=vertical(plumbline)
    enddo
    return
End Subroutine ITRF2LOFAR
!================================================
Subroutine RelDist(Source,LFRAnt,RDist)
!   RDist is distance relative to center of CS002, in units of time-samples
!  RDist= dt_(core-source) - dt_(ant-source)=dt_{cs} - dt_{as}
! t_core=t_s + dt_cs ; t_ant=t_s + dt_as  ; t_s=t_ant-dt_as
! thus: t_core=t_ant + Rdist
    use constants, only : dp,sample,c_mps
    implicit none
    real(dp), intent(in) :: Source(3)
    real(dp), intent(in) :: LFRAnt(3)
    real(dp), intent(out) :: RDist
    Real(dp) :: D
   Real(dp) :: IndxRefrac
   Real(dp), external :: RefracIndex
    !write(2,*) 'source',source
    !write(2,*) 'LFRAnt',LFRAnt
    D=sum((Source(:)-LFRAnt(:))*(Source(:)-LFRAnt(:)))
    RDist = sqrt(D)
    D=sum(Source(:)*Source(:))  ! core=center station CS002 is at zero
    RDist = RDist - sqrt(D)         ! units of [meter]
    RDist = Rdist*RefracIndex( Source(3) )/(c_mps*sample) ! Convert to units of [samples]
    Return
End Subroutine RelDist
!================================================
Real(kind=8) Function tShift_ms(Source)  !
    use constants, only : dp,sample,c_mps
    implicit none
    real(dp), intent(in) :: Source(1:3)
    Real(dp), external :: RefracIndex
    !real(dp), intent(out) :: tShift
    tShift_ms = RefracIndex(Source(3))*sqrt(Source(1)*Source(1)+Source(2)*Source(2)+Source(3)*Source(3))*1000.d0/c_mps
    Return
End Function tShift_ms
!================================================    Real(dp), external ::
Real(kind=8) Function tShift_smpl(Source)  !
    use constants, only : dp,sample,c_mps
    implicit none
    real(dp), intent(in) :: Source(1:3)
    Real(dp), external :: RefracIndex
    !real(dp), intent(out) :: tShift
    tShift_smpl = RefracIndex(Source(3))*sqrt(Source(1)*Source(1)+Source(2)*Source(2)+Source(3)*Source(3))/(c_mps*sample)
    Return
End Function tShift_smpl
!================================================
Real(kind=8) Function RefracIndex(h)
!        Parameters of atmospheric density as used in Corsika
    use constants, only : dp,HeightCorrectIndxRef
    implicit none
    integer, parameter :: AtmHei_Dim=1000
    Real(dp), parameter :: AtmHei_step=10.d0 ! [m]
    real(dp), intent(in) :: h
    Real(dp), save :: Xi(0:AtmHei_Dim)
    Real(dp) :: hi
    Logical, save :: First=.true.
    integer :: i
    !
    If (First) Then ! Calculate atmosphere (from MGMR3D)
      First=.false.
      Call AveIndxRefr(AtmHei_dim,AtmHei_step, xi)
    EndIf
    If(HeightCorrectIndxRef) Then
       hi=h/AtmHei_step
       i=INT(hi)
       If(i.ge.AtmHei_Dim) then
         RefracIndex=xi(AtmHei_Dim)
       ElseIf(i.gt.0) then
         RefracIndex=xi(i)*(hi-i) + xi(i+1)*(i+1.-hi)
       Else
         RefracIndex=Xi(0)
       EndIf
    Else
       RefracIndex=Xi(0)
    EndIf
    !Write(2,*) 'Function RefracIndex:',h,RefracIndex
End Function RefracIndex
!=========================================
Subroutine AveIndxRefr(AtmHei_dim,AtmHei_step, xi)
   use constants, only : dp
   implicit none
   integer, intent(in) :: AtmHei_Dim
   Real(dp), intent(in) :: AtmHei_step ! [m]
   Real(dp), intent(out) :: Xi(0:AtmHei_Dim)
   !REAL(dp), parameter :: Refrac=1.000267394d0 ! Index of refraction at sea level (20^0, 100% hum) From https://emtoolbox.nist.gov/Wavelength/Documentation.asp
   REAL(dp), parameter :: Refrac= 1.0003     ! Index of refraction at sea level
   real(dp) :: height
   real(dp):: mass,b,c,RefOrho
   integer :: i
   !
   xi(0)=Refrac
   mass=0.
   RefOrho=(xi(0)-1.d0)*9941.8638d0/1222.6562d0
   do i = 1, AtmHei_dim
      height=i*AtmHei_step  ! distance to ground along shower axis
       if (height.ge.10d3) then
            b = 1305.5948; c = 6361.4304
       elseif (height.ge.4d3) then
            b = 1144.9069d0; c = 8781.5355d0
       else
            b = 1222.6562d0; c = 9941.8638d0
       endif
      mass = mass + (b/c) * exp(-height/c) * AtmHei_step ! assume density is about constant
      !write(2,*) height, a, PenDepth(0)-X_rh-RPenDepth(i)
      !
      ! calculate  averaged refractivity per meter
      xi(i) = 1.d0+ RefOrho* mass/height   ! mean index of refraction for rays from height to 0
    end do
    !write(2,*) 'Xi:',Xi(0),Xi(1),Xi(10),Xi(100),Xi(1000)
    Return
End Subroutine AveIndxRefr
!=========================================
Real(kind=8) Function SubRelDist(SrcPos,i_ant,i_chunk)
   use Chunk_AntInfo, only : Ant_pos, Ant_RawSourceDist
   use constants, only : dp
   Implicit none
   Real(dp), intent(in) ::  SrcPos(3)
   Integer, intent(in) :: i_ant,i_chunk
   Real(dp) :: RDist
      Call RelDist(SrcPos(1),Ant_pos(1,i_ant,i_chunk),RDist)  ! source position may have changed compared to previous
      SubRelDist=Rdist - Ant_RawSourceDist(i_ant,i_chunk)
End Function SubRelDist
!===============================================
    Subroutine GetNonZeroLine(lineTXT)
    implicit none
    character(LEN=*), intent(inout) :: lineTXT
    Character(len=6) :: FMT
    integer :: eof,N
      N=LEN(lineTXT)
      write(FMT,"(A,I3.3,A)") '(A',N,')'
      !write(2,*) 'GetNonZeroLine2:', FMT, N
      lineTXT=''
      Do while (trim(lineTXT).eq.'')
         read(*,FMT,iostat=eof) lineTXT
         if(eof.lt.0) exit
      enddo
    End Subroutine GetNonZeroLine
!===============================================
    Subroutine GetMarkedLine(Mark,lineTXT)
    implicit none
    character(LEN=*), intent(inout) :: lineTXT
    character(LEN=*), intent(out) :: Mark
    Character(len=6) :: FMT
    integer :: eof,N,j
      N=LEN(lineTXT)
      write(FMT,"(A,I3.3,A)") '(A',N,')'
      !write(2,*) 'GetNonZeroLine2:', FMT, N
   1  continue
      lineTXT=''
      Do while (trim(lineTXT).eq.'')
         read(*,FMT,iostat=eof) lineTXT
         if(eof.lt.0) Return
      enddo
      If(lineTXT(1:1).eq.'!') goto 1   ! skip comment lines
      j = iachar(lineTXT(1:1))
      if (j>= iachar("a") .and. j<=iachar("z") ) then
         j = (j-32)   ! convert to upper case if not already
      end if
      Mark=achar(j)
      lineTXT(1:N-1)=lineTXT(2:N)
    End Subroutine GetMarkedLine
!==========================================
!==========================================
Module ansi_colors
  implicit none

  character(len=1), parameter :: c_esc = achar(27)  !  =  \u001b
  character(len=2), parameter :: c_start = c_esc // '['
  character(len=1), parameter :: c_end = 'm'
  character(len=*), parameter :: c_black = '30'
  character(len=*), parameter :: c_red = '31'
  character(len=*), parameter :: c_green = '32'
  character(len=*), parameter :: c_yellow = '33'
  character(len=*), parameter :: c_blue = '34'
  character(len=*), parameter :: c_magenta = '35'
  character(len=*), parameter :: c_cyan = '36'
  character(len=*), parameter :: c_white = '37'
  character(len=*), parameter :: c_clear = c_start // '0' // c_end

! ANSI color code:  Black: \u001b[30m
!    [90m=dark grey           [30m=black
![91m=peach               [31m=red
![92m=light green         [32m=green
![93m=light yellow        [33m=yellow
![94m=light blue          [34m=blue
![95m=pink                [35m=purple
![96m=light aqua          [36m=aqua
![97m=pearl white
!For example, the 8 background colors correspond to the codes:
!
!    Background Black: \u001b[40m
!    Background Red: \u001b[41m
!    Background Green: \u001b[42m
!    Background Yellow: \u001b[43m
!    Background Blue: \u001b[44m
!    Background Magenta: \u001b[45m
!    Background Cyan: \u001b[46m
!    Background White: \u001b[47m
!With the bright versions being (however this seems not to work for me):
!    Background Bright Black: \u001b[40;1m
!    Background Bright Red: \u001b[41;1m
!    Background Bright Green: \u001b[42;1m
!    Background Bright Yellow: \u001b[43;1m
!    Background Bright Blue: \u001b[44;1m
!    Background Bright Magenta: \u001b[45;1m
!    Background Bright Cyan: \u001b[46;1m
!    Background Bright White: \u001b[47;1m
!    Bold: \u001b[1m        see also https://en.wikipedia.org/wiki/ANSI_escape_code
!    Underline: \u001b[4m
!    Reversed: \u001b[7m
!Or together print u"\u001b[1m\u001b[4m\u001b[7m BOLD Underline Reversed \u001b[0m"
! cursor: (see: http://www.lihaoyi.com/post/BuildyourownCommandLinewithANSIescapecodes.html)
!    Up: \u001b[{n}A
!    Down: \u001b[{n}B
!    Right: \u001b[{n}C
!    Left: \u001b[{n}D
!as used in  "\u001b[1000D" + str(i + 1) + "%"
!And reset is the same:
!    Reset: \u001b[0m
contains

  Function color(str, code) result(out)
    character(len=*), intent(in) :: str
    character(len=*), intent(in) :: code
    character(len=:), allocatable :: out
    out = c_start // code // c_end // str // c_clear
  End Function color

End Module ansi_colors
!==========================================
Subroutine Convert2m(CenLoc)
   Implicit none
   real*8, intent(inout) :: CenLoc(3)
   If((abs(CenLoc(1)) .lt. 180.) .and. (abs(CenLoc(2)) .lt. 180.) .and. (abs(CenLoc(3)) .lt. 20.) ) Then
      CenLoc(:)=CenLoc(:)*1000.  ! convert from [km] to [m]
   Endif
   Return
End Subroutine Convert2m
!--------------------
Pure Subroutine SetSmooth(N_Smth, Smooth)
   use constants, only : dp
   Implicit none
   Integer, intent(in) :: N_Smth
   Real(dp), intent(out) :: Smooth(-N_smth:N_smth)
   integer :: i
   !Allocate( Smooth(-N_smth:N_smth) )
   Smooth(:)=0.
   Smooth(0)=1.
   Smooth(N_smth/2)=0.5  ! needed for block profile for even values of N_smth
   Do i=1,N_smth
      !Smooth(i)=(1.-(Real(i)/N_smth)**2)**2.5  ! Parabola to power to make it =0.5 at N_smth/2
      If(i.lt.(N_smth+1)/2) Smooth(i)=1.      ! Block gives indistinguishable results from parabola; parabola has a tiny bit smaller volume around peak.
      Smooth(-i)=Smooth(i)
   EndDo
   !SmPow=SUM(Smooth(:))
   Smooth(:)=Smooth(:)/SUM(Smooth(:))
   Return
End Subroutine SetSmooth
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
Subroutine random_stdnormal3D(x)
!  https://masuday.github.io/fortran_tutorial/random.html
! General interest: https://en.wikibooks.org/wiki/Fortran/Fortran_procedures_and_functions
   implicit none
   real, intent(out) :: x(1:3)
   real :: R,D, Epsln=epsilon(Epsln)
   real :: u1,u2,u3
   Real :: random_stdnormal
   ! call random_number(r) gives 0=< r < 1 i.e. including 0, excluding 1
   call random_number(u1)  ! may include zero
   call random_number(u2)
   call random_number(u3)
   R=sqrt((0.5-u1)**2+(0.5-u2)**2+(0.5-u3)**2+Epsln) ! to prevent zero
   D=random_stdnormal()! can be zero
   D=((abs(d))**(1/3.))  ! to have distances distributed like [d^2 x gaussian(d)]
   !D=sqrt(abs(d)) * Space_width  ! to have distances distributed like [d x gaussian(d)]
   x(1)= D*(0.5-u1)/R
   x(2)= D*(0.5-u2)/R
   x(3)= D*(0.5-u3)/R
   Return
end Subroutine random_stdnormal3D
