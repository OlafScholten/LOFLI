MODULE LOFLI_Input

!   USE ISO_FORTRAN_ENV
    use constants, only : dp

   IMPLICIT NONE

   INTERFACE PrintValues
      MODULE PROCEDURE PrintIntArray
      MODULE PROCEDURE PrintIntVal
      MODULE PROCEDURE PrintRealVal
      MODULE PROCEDURE PrintLogiVal
      MODULE PROCEDURE PrintChArray
      MODULE PROCEDURE PrintChVal
   END INTERFACE PrintValues

   CONTAINS

   Subroutine PrintIntArray(Var, Var_name, Label)
      Integer, dimension(:), Intent(in) :: Var
      Character(len=*), Intent(in) :: Var_name
      Character(len=*), Intent(in), optional :: Label
      Integer :: N,i
      Character(len=20) :: VarCh
      N=size(Var)
      write(2,"(1x,A,A)", ADVANCE='NO') TRIM(Var_Name),'= '
      Do i=1,N
         write(VarCh,*) Var(i)
         If(Var(i).eq.Var(N)) exit
         write(2,"(A,A)", ADVANCE='NO') Trim(ADJUSTL(VarCh)),', '
      Enddo
      If(present(Label)) Then
         write(2,"(I3,A,A,A,5x,'! ',A)") N-i+1,'*( ',Trim(ADJUSTL(VarCh)),' )',Label
      Else
         write(2,"(I3,A,A,A)") N-i+1,'*( ',Trim(ADJUSTL(VarCh)),' )'
      EndIf
   End Subroutine PrintIntArray
!
   Subroutine PrintIntVal(Var, Var_name, Label)
      Integer, Intent(in) :: Var
      Character(len=*), Intent(in) :: Var_name
      Character(len=*), Intent(in), optional :: Label
      Character(len=20) :: VarCh
      write(VarCh,*) Var
      If(present(Label)) Then
         write(2,"(1x,3A,5x,'! ',A)") TRIM(Var_Name),'= ', Trim(ADJUSTL(VarCh)),Label
      Else
         write(2,"(1x,3A)") TRIM(Var_Name),'= ', Trim(ADJUSTL(VarCh))
      EndIf
      !write(2,"(1x,A,A,A)") TRIM(Var_Name),'= ', Trim(ADJUSTL(VarCh))
   End Subroutine PrintIntVal
!
   Subroutine PrintRealVal(Var, Var_name, Label)
      Real(dp), Intent(in) :: Var
      Character(len=*), Intent(in) :: Var_name
      Character(len=*), Intent(in), optional :: Label
      Character(len=50) :: VarCh
      write(VarCh,*) Var
      If(present(Label)) Then
         write(2,"(1x,3A,5x,'! ',A)") TRIM(Var_Name),'= ', Trim(ADJUSTL(VarCh)),Label
      Else
         write(2,"(1x,3A)") TRIM(Var_Name),'= ', Trim(ADJUSTL(VarCh))
      EndIf
      !write(2,"(1x,A,A,A)") TRIM(Var_Name),'= ', Trim(ADJUSTL(VarCh))
   End Subroutine PrintRealVal
!
   Subroutine PrintLogiVal(Var, Var_name, Label)
      Logical, Intent(in) :: Var
      Character(len=*), Intent(in) :: Var_name
      Character(len=*), Intent(in), optional :: Label
      Character(len=20) :: VarCh
      write(VarCh,*) Var
      If(present(Label)) Then
         write(2,"(1x,3A,5x,'! ',A)") TRIM(Var_Name),'= ', Trim(ADJUSTL(VarCh)),Label
      Else
         write(2,"(1x,3A)") TRIM(Var_Name),'= ', Trim(ADJUSTL(VarCh))
      EndIf
      !write(2,"(1x,A,A,A)") TRIM(Var_Name),'= ', Trim(ADJUSTL(VarCh))
   End Subroutine PrintLogiVal
!
   Subroutine PrintChArray(Var, Var_name, Label)
      Character(len=*), dimension(:), Intent(in) :: Var
      Character(len=*), Intent(in) :: Var_name
      Character(len=*), Intent(in), optional :: Label
      Integer :: N,i
      Character(len=20) :: VarCh
      N=size(Var)
      write(2,"(1x,A,A)", ADVANCE='NO') TRIM(Var_Name),'= '
      Do i=1,N
         write(VarCh,"(A)") Var(i)
         If(Var(i).eq.Var(N)) exit
         write(2,"(A,A,A)", ADVANCE='NO') '"',Trim(VarCh),'", '
      Enddo
      If(present(Label)) Then
         write(2,"(I3,A,A,A,5x,'! ',A)") N-i+1,'*( "',Trim(VarCh),'" )',Label
      Else
         write(2,"(I3,A,A,A)") N-i+1,'*( "',Trim(VarCh),'" )'
      EndIf
   End Subroutine PrintChArray
!
   Subroutine PrintChVal(Var, Var_name, Label)
      Character(len=*), Intent(in) :: Var
      Character(len=*), Intent(in) :: Var_name
      Character(len=*), Intent(in), optional :: Label
      Character(len=200) :: VarCh
      write(VarCh,"(A)") Var
      If(present(Label)) Then
         write(2,"(1x,4A,5x,'! ',A)") TRIM(Var_Name),'= "',Trim(VarCh),'"',Label
      Else
         write(2,"(1x,4A)") TRIM(Var_Name),'= "',Trim(VarCh),'"'
      EndIf
      !write(2,"(1x,4A)") TRIM(Var_Name),'= "',Trim(VarCh),'"'
   End Subroutine PrintChVal
!
   Subroutine PrintParIntro(Label, RunOption, OutFileLabel)
      Character(len=*), Intent(in) :: Label
      Character(len=*), Intent(in) :: RunOption
      Character(len=*), Intent(in) :: OutFileLabel
      write(2,"(20(1x,'='))")
      write(2,*) Label,' run with following useful parameters in namelist: " &Parameters":'
      write(2,"(1x,A)", ADVANCE='NO') '&Parameters '
      Call PrintValues(RunOption,'RunOption')
      Call PrintValues(OutFileLabel,'OutFileLabel')
   End Subroutine PrintParIntro
!
   Subroutine ReadSourceTimeLoc(StartTime_ms, CenLoc)
      use constants, only : dp, Sample, c_mps
      use Chunk_AntInfo, only : TimeBase
      Implicit none
      Real(dp), intent(OUT) :: StartTime_ms, CenLoc(1:3)
      Integer, parameter ::lnameLen=180
      Character(LEN=lnameLen) :: lname
      Real(dp) :: t_shft
      Integer :: j
      Real(dp), external :: tShift_ms
      !
      Call GetNonZeroLine(lname)
      Read(lname(2:lnameLen),*) StartTime_ms, CenLoc  ! Start time
      !Call PrintValues(lname,'input line-1', 'time & position' )
      write(2,"(A)") 'Input line-1: "'//lname(1:1)//'|'//TRIM(lname(2:lnameLen))// &
         '" !  Core/Source-| time, & position'
      Call Convert2m(CenLoc)
      StartTime_ms=StartTime_ms+TimeBase
      t_shft=tShift_ms(CenLoc(:)) ! sqrt(SUM(CenLoc(:)*CenLoc(:)))*1000.*Refrac/c_mps ! in mili seconds due to signal travel distance
      !
      j = iachar(lname(1:1))  ! convert to upper case if not already
      if (j>= iachar("a") .and. j<=iachar("z") ) then
         lname(1:1) = achar(j-32)
      end if
      SELECT CASE (lname(1:1))
         CASE("S")  ! time at Source (central voxel) is given
            StartTime_ms=StartTime_ms+t_shft
         CASE DEFAULT  ! time at reference antenna is given
      End SELECT
      write(2,"(A,F12.6,A,F12.6,A,A,2(F9.4,','),F9.4,A)") &
         ' True start time trace, adding base, at core (at source)=', StartTime_ms, ' (',StartTime_ms-t_shft,') [ms]',&
         ' for source @(N,E,h)=(', CenLoc(1:3)/1000., ' ) km'
      Return
   End Subroutine ReadSourceTimeLoc
   !
END MODULE LOFLI_Input
