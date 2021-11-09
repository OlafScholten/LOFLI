!    Include 'HDF5_LOFAR_Read.f90'
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
Subroutine Station_ID2Calib(STATION_ID,Ant_ID,StatAnt_Calib)
    use constants, only : dp,sample
    use DataConstants, only : Station_nrMax, DataFolder, Diagnostics, Ant_nrMax, Calibrations
    use DataConstants, only : Polariz
    use StationMnemonics, only : Station_ID2Mnem
    use FitParams, only : Explore
    Implicit none
    Integer, intent(in) :: STATION_ID
    Integer, intent(in) :: Ant_ID
    Real(dp), intent(out) :: StatAnt_Calib
    Character(len=5) :: txt, Station_Mnem
    Real(dp) :: Delay
    Character(len=5), save :: LOFAR_STMnm(Station_nrMax), Fine_STMnm(Station_nrMax)
    Real(dp), save :: LOFAR_STDelay(Station_nrMax), Fine_STDelay(Station_nrMax)
    Real(dp), save :: Fine_AntDelay(1:Ant_nrMax)
    Integer, save :: SAI_AntDelay(1:Ant_nrMax),Nr_AntDelay
    integer :: nxx, i, SAI
    Integer, save :: i_LOFAR, i_fine
    Logical, save :: First=.true.
    integer :: Ant(1),SAI_partner
    Logical :: Old,core
    !
    If(first) then
        Old=.false.
        If(trim(Calibrations).eq.'') Old=.true.
        i_LOFAR=0
        Open(unit=9,STATUS='old',ACTION='read', FILE = 'Book/StationCalibrations.dat', IOSTAT=nxx) ! trim(DataFolder)//
        ! distilled from LOFAR data-base with "CalibrationTablesGenerator.f90"
        If(nxx.ne.0) then
            write(2,*) 'system station calibrations missing! expected file =', &
               'StationCalibrations.dat'
            stop 'Station_ID2Calib'
        else
            Do
                read(9,*,IOSTAT=nxx) txt, Delay
                if(nxx.ne.0) then
                    exit
                else
                    i_LOFAR=i_LOFAR+1
                    LOFAR_STMnm(i_LOFAR) = txt
                    LOFAR_STDelay(i_LOFAR)= Delay/sample  ! Convert to samples
                endif
            Enddo
        endif
        close(unit=9)
        !
        Nr_AntDelay=0
        If(Old) then
            Open(unit=9,STATUS='old',ACTION='read', FILE = 'Book/AntCalibrations.dat', IOSTAT=nxx) ! trim(DataFolder)//
        Else
            Open(unit=9,STATUS='old',ACTION='read', FILE = 'Book/'//trim(Calibrations), IOSTAT=nxx)
        Endif
        ! Generated by fitting discreate sources
        If(nxx.ne.0) then
            write(2,*) 'problems with file=','AntCalibrations.dat'
            Flush(unit=2)
        else
            Do
                read(9,*,IOSTAT=nxx) SAI, Delay
                if(nxx.ne.0) then
                    exit
                else
                    Nr_AntDelay=Nr_AntDelay+1
                    SAI_AntDelay(Nr_AntDelay) = SAI
                    !If(Polariz) then  ! set equal delays for even-odd pairs, independent of reading order
                    !    if(mod(SAI,2) .eq.1) then
                    !       SAI_partner=SAI-1    !compare with SAI-1
                    !    Else
                    !       SAI_partner=SAI+1    !compare with SAI-1
                    !    EndIf
                    !    Ant=MAXLOC(SAI_AntDelay(1:Nr_AntDelay), MASK = SAI_AntDelay .eq. SAI_partner)
                    !    If(SAI_AntDelay(Ant(1)) .eq. SAI_partner) Then
                    !       If(ABS(Delay - Fine_AntDelay(Ant(1))).gt.0.2) write(2,*) 'Station_ID2Calib',SAI, SAI_partner,&
                    !             (Delay - Fine_AntDelay(Ant(1))),(Delay + Fine_AntDelay(Ant(1)))/2.
                    !       Delay=(Delay + Fine_AntDelay(Ant(1)))/2.
                    !       Fine_AntDelay(Ant(1))= Delay
                    !    EndIf
                    !EndIf   ! Polariz
                    Fine_AntDelay(Nr_AntDelay)= Delay ! already in [samples]
                endif
            enddo
        endif
        If(Old) close(unit=9)
        !
        nxx=0
        i_fine=0
        If(Old) Open(unit=9,STATUS='old',ACTION='read', FILE = 'Book/FineCalibrations.dat', IOSTAT=nxx) ! trim(DataFolder)//
        ! Generated by fitting discreate sources
        If(nxx.ne.0) then
            write(2,*) 'problems with file=','FineCalibrations.dat'
        else
            Do
                read(9,*,IOSTAT=nxx) txt, Delay
                if(nxx.ne.0) then
                    exit
                else
                    i_fine=i_fine+1
                    Fine_STMnm(i_fine) = txt
                    Fine_STDelay(i_fine)= Delay ! already in samples
                endif
            enddo
        endif
        close(unit=9)
        !
        write(2,*) 'Nr_AntDelay=',Nr_AntDelay
        first=.false.
    endif   ! (first)
    !
    nxx=1
    Call Station_ID2Mnem(STATION_ID,Station_Mnem)
    core=.false.    !  (Explore .and. (Station_Mnem(1:2) .eq. 'CS'))
    ! zero additional time-calibration offsets for the core in the very first run
    Do i=1,i_LOFAR
        If(LOFAR_STMnm(i).eq.Station_Mnem) then
            StatAnt_Calib=LOFAR_STDelay(i)
            nxx=0
            exit
        endif
    Enddo
    if(nxx.ne.0) then
        StatAnt_Calib=0.
        If(Diagnostics) Write(2,*) '**** No LOFAR-calibration timing found for',STATION_ID,Station_Mnem
    endif
    !
    If(core) then  !  Skip adding fine_delays
      Return
    Endif
    !write(2,*) 'i_Fine=',i_Fine
    Do i=1,i_Fine
        !write(2,*) '-',Fine_STMnm(i),'-',Station_Mnem,'-'
        If(Fine_STMnm(i).eq.Station_Mnem) then
            !write(2,*) 'station',i,Station_Mnem,StatAnt_Calib, Fine_STDelay(i)
            StatAnt_Calib= StatAnt_Calib + Fine_STDelay(i)
            exit
        endif
    Enddo
    !
    If(Nr_AntDelay.eq.0) then
      write(2,*) 'There are no antenna delays specified !!!!!!!!!!!!!!'
      StatAnt_Calib= StatAnt_Calib
      Return
    Endif
    SAI=1000*station_ID + Ant_ID
    Ant=MAXLOC(SAI_AntDelay(1:Nr_AntDelay), MASK = SAI_AntDelay .eq. SAI)
    If(SAI_AntDelay(Ant(1)) .eq. SAI) StatAnt_Calib= StatAnt_Calib + Fine_AntDelay(Ant(1))
    Return
End Subroutine Station_ID2Calib
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
    use constants, only : dp,sample,c_mps,Refrac
    implicit none
    real(dp), intent(in) :: Source(3)
    real(dp), intent(in) :: LFRAnt(3)
    real(dp), intent(out) :: RDist
    Real(dp) :: D
    !write(2,*) 'source',source
    !write(2,*) 'LFRAnt',LFRAnt
    D=sum((Source(:)-LFRAnt(:))*(Source(:)-LFRAnt(:)))
    RDist = sqrt(D)
    D=sum(Source(:)*Source(:))  ! core=center station CS002 is at zero
    RDist = RDist - sqrt(D)         ! units of [meter]
    RDist = Rdist*Refrac/(c_mps*sample) ! Convert to units of [samples]
    Return
End Subroutine RelDist
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
Subroutine Double_IR_sort(N,a,R)
! Adepted From:  https://rosettacode.org/wiki/Sorting_algorithms/Selection_sort#Fortran
!  sorts integer array a according to the first element
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
Subroutine Double_RI_sort(N,a,R)
! Adepted From:  https://rosettacode.org/wiki/Sorting_algorithms/Selection_sort#Fortran
!  sorts integer array a according to the first element
! sort also the additional real array
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
   If((abs(CenLoc(1)) .lt. 100.) .and. (abs(CenLoc(1)) .lt. 100.) .and. (abs(CenLoc(1)) .lt. 100.) ) Then
      CenLoc(:)=CenLoc(:)*1000.  ! convert from [km] to [m]
   Endif
   Return
End Subroutine Convert2m
