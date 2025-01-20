!=================================
Module unque

   INTERFACE GenSort ! sorts according to ascending values
      MODULE PROCEDURE Selection_sort  ! [1d integer array]
                           ! sort(n, a)
      MODULE PROCEDURE Double_sort     ! the first column (*,1) of a [2d integer array] and rearrange the others accordingly
      MODULE PROCEDURE Double_IR_sort  ! [N] elements of [1d integer] and also a [1d real(dp)] array
                           ! Double_IR_sort(N,a,R)
                           !  sorts integer array a according to the first dimension & reorders the additional real array
      MODULE PROCEDURE Double_RI_sort  ! N elements of 1d real(dp) and also a 1d integer array
                           !  sorts Real array a (min first) and reorders the additional Integer array in the same way
      MODULE PROCEDURE sort            ! N elements of 1r real
      MODULE PROCEDURE Double_dp_sort  ! the first column of a [2d real(dp) array] and rearrange the others accordingly
      MODULE PROCEDURE HPSORT_mult_RI  ! HPSORT_mult_RI(RA,IA,N); RA(*,1:N) [real] and IA(*,1:N) [integer]	 rearranged like  RA(1,1:N)
      MODULE PROCEDURE SortPerm        ! SortPerm(RA(*,1:N),N,IA)  ;  keeps track of the permutations such that RA(IA(j)) =original RA(j)

   END INTERFACE GenSort


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
!  Double sorts a 2d integer array a according to the first element
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
Subroutine Double_dp_sort(a)
! Adepted From:  https://rosettacode.org/wiki/Sorting_algorithms/Selection_sort#Fortran
!  Double sorts a double array [a] according to the first column, i.e. A(*,1)
    use constants, only : dp
    Real(dp), INTENT(IN OUT) :: a(:,:)
    INTEGER :: i, minIndex
    Real(dp) :: temp(SIZE(a,2))
    !
    DO i = 1, SIZE(a,1)-1
       minIndex = MINLOC(a(i:,1), 1) + i - 1
       IF (a(i,1) > a(minIndex,1)) THEN ! sort on first column
          temp(:) = a(i,:)
          a(i,:)= a(minIndex,:)
          a(minIndex,:) = temp(:)
       End IF
    End DO
End Subroutine Double_dp_sort
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
!  sorts Real array a (min first) and reorders the additional Integer array in the same way
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
SUBROUTINE HPSORT_mult_RI(RA,IA,N)  ! should move elsewhere at some point
!*****************************************************
!*  Sorts an array RA of length N in ascending order *
!*                by the Heapsort method             *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*	    N	  size of table RA                                                                  *
!*          RA(1,1:N)	  table to be sorted, real values                                    *
!*          RA(*,1:N) [real] and IA(*,1:N) [integer]	 rearranged like  RA(1,1:N)            *
!* OUTPUT:                                                                                   *
!*	    RA    table sorted in ascending order                                                 *
!*	    IA    table rearranged following order RA                                             *
!*                                                   *
!* NOTE: The Heapsort method is a N Log2 N routine,  *
!*       and can be used for very large arrays.      *
!*****************************************************
  IMPLICIT none
  integer, intent(in) :: N
  real*8 RA(:,:)
  Integer IA(:,:)
  real*8 RRA(SIZE(RA,1))
  Integer RIA(SIZE(IA,1))
  integer :: d1,d2,L,IR,I,J
  d1=SIZE(RA,1)
  d2=SIZE(IA,1)
  L=N/2+1
  IR=N
  !The index L will be decremented from its initial value during the
  !"hiring" (heap creation) phase. Once it reaches 1, the index IR
  !will be decremented from its initial value down to 1 during the
  !"retirement-and-promotion" (heap selection) phase.
10 continue
  if(L > 1)then
    L=L-1
    RRA(1:d1)=RA(1:d1,L)
    RIA(1:d2)=IA(1:d2,L)
  else
    RRA(1:d1)  =RA(1:d1,IR)
    RA(1:d1,IR)=RA(1:d1,1)
    RIA(1:d2)  =IA(1:d2,IR)
    IA(1:d2,IR)=IA(1:d2,1)
    IR=IR-1
    if(IR.eq.1)then
      RA(1:d1,1)=RRA(1:d1)
      IA(:,1)=RIA(:)
      return
    end if
  end if
  I=L
  J=L+L
20 if(J.le.IR)then
  if(J < IR)then
    if(RA(1,J) < RA(1,J+1))  J=J+1
  end if
  if(RRA(1) < RA(1,J))then
    RA(1:d1,I)=RA(1:d1,J)
    IA(:,I)=IA(:,J)
    I=J; J=J+J
  else
    J=IR+1
  end if
  goto 20
  end if
  RA(1:d1,I)=RRA(1:d1)
  IA(:,I)=RIA(:)
  goto 10
END SUBROUTINE HPSORT_mult_RI
!-------------------------------------
Subroutine SortPerm(RA,N,IA)
!*****************************************************
!*  Sorts an array RA of length N in ascending order *
!*                by the Heapsort method             *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*	    N	  size of table RA                                                                  *
!*          RA(1:N)	  table to be sorted, real values                                    *
!*          IA(1:N) [real] and IA(*,1:N) [integer]	 rearranged like  RA(1,1:N)            *
!* OUTPUT:                                                                                   *
!*	    RA    table sorted in ascending order                                                 *
!*	    IA    table rearranged following order RA                                             *
!*                                                   *
!* NOTE: The Heapsort method is a N Log2 N routine,  *
!*       and can be used for very large arrays.      *
!*****************************************************
   use constants, only : dp
   implicit none
   INTEGER, INTENT(IN) :: N  ! length of arrays a and R that should be sorted, actual size may be larger
   Real(dp), INTENT(IN OUT) :: RA(:,:)
   Integer, INTENT(OUT) :: IA(:)  ! keeps track of the permutations such that RA(IA(j)) =original RA(j)
   real(dp) :: RRA(SIZE(RA,1))
   integer :: d1,L,IR,I,J
   Integer SI
   d1=SIZE(RA,1)
   Do i=1,N
      IA(i)=i
   EndDo
   L=N/2+1
   IR=N
   !The index L will be decremented from its initial value during the
   !"hiring" (heap creation) phase. Once it reaches 1, the index IR
   !will be decremented from its initial value down to 1 during the
   !"retirement-and-promotion" (heap selection) phase.
10 continue
   if(L > 1)then
      L=L-1
      RRA(1:d1)=RA(1:d1,L)
      SI=IA(L)
   else
      RRA(1:d1)  =RA(1:d1,IR)
      RA(1:d1,IR)=RA(1:d1,1)
      SI  =IA(IR)
      IA(IR)=IA(1)
      IR=IR-1
      if(IR.eq.1)then
         RA(1:d1,1)=RRA(1:d1)
         IA(1)=SI
         return
      end if
   end if
   I=L
   J=L+L
20 if(J.le.IR)then
      if(J < IR)then
         if(RA(1,J) < RA(1,J+1))  J=J+1
      end if
      if(RRA(1) < RA(1,J))then
         RA(1:d1,I)=RA(1:d1,J)
         IA(I)=IA(J)
         I=J; J=J+J
      else
         J=IR+1
      end if
      goto 20
   end if
   RA(1:d1,I)=RRA(1:d1)
   IA(I)=SI
   goto 10
!
End Subroutine SortPerm
!=================================
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
