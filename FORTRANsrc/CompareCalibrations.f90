    Include 'ConstantsModules.f90'
    Include 'ParamModules.f90'
    Include 'MappingUtilities.f90'
!-----------------------------------
!=================================
Module Interferom_Pars  ! Just placeholder
   Integer, allocatable, save :: IntFer_ant(:,:)
End Module Interferom_Pars
!-----------------------------------
Program CompareCalibr
   use constants, only : dp
   use DataConstants, only : Ant_nrMax, Station_nrMax
   use DataConstants, only : Station_nrMax, DataFolder, Calibrations! , Diagnostics
   use Calibration, only : ReadCalib, Fine_STMnm, Fine_STDelay, StationInCal
   Implicit none
   !
   Real(dp) :: Fine_STDelayA(1:Station_nrMax), Fine_STDelayB(1:Station_nrMax)! time calibration data
   Real(dp) :: offset
   Character(len=85) :: CalibrationsA, CalibrationsB
   Integer :: i,StationInCalA(Station_nrMax),StationInCalB(Station_nrMax)
   !
   Open(unit=2,STATUS='unknown',ACTION='write', FILE ="CalComp.out")
   Read(*,*)  CalibrationsA
   Read(*,*)  CalibrationsB
   Fine_STMnm(:)=''
   !
   !Note: it is assumed that the order of the station is the same for the two calibration files
   Calibrations=CalibrationsA
   Call ReadCalib()
   Fine_STDelayA(:)=Fine_STDelay(:)
   StationInCalA(:)=StationInCal(:)
   !
   Calibrations=CalibrationsB
   Call ReadCalib()
   Fine_STDelayB(:)=Fine_STDelayA(:) - Fine_STDelay(:)
   StationInCalB(:)=StationInCal(:)
   !
   Do i=1,Station_nrMax
      If(StationInCalA(i).eq.1 .and. StationInCalB(i).eq.1) Then
         offset=Fine_STDelayB(i)
         Fine_STDelayB(:)=Fine_STDelayB(:)-offset
         exit
      EndIf
   Enddo
   !
   write(2,*) TRIM(CalibrationsA),' minus ', TRIM(CalibrationsB), ' in [ns], offset=',offset*5
   !write(2,"(40A7)") Fine_STMnm(:)
   !write(2,"(40F7.2)") 5*Fine_STDelayB(:)
   Do i=1,Station_nrMax
      If(Fine_STMnm(i)(1:1).eq.'C') Then
         write(2,"(A6)", ADVANCE='no') Fine_STMnm(i)
      EndIf
   Enddo
   write(2,*) '  '
   !
   Do i=1,Station_nrMax
      If(Fine_STMnm(i)(1:1).eq.'C') Then
         write(2,"(F6.2)", ADVANCE='no') 5*Fine_STDelayB(i)
      EndIf
   Enddo
   write(2,*) '  '
   !
   Do i=1,Station_nrMax
      If(Fine_STMnm(i)(1:1).eq.'C') Then
         If(StationInCalA(i).eq.1 .and. StationInCalB(i).eq.1) Then
            write(2,"(3x,'1*1')", ADVANCE='no')
         Else
            write(2,"(2x,i2,i2)", ADVANCE='no') StationInCalA(i),StationInCalB(i)
         EndIf
      EndIf
   Enddo
   write(2,*) '  '
   !
   Do i=1,Station_nrMax
      If(Fine_STMnm(i)(1:1).eq.'R') Then
         write(2,"(A6)", ADVANCE='no') Fine_STMnm(i)
      EndIf
   Enddo
   write(2,*) '  '
   !
   Do i=1,Station_nrMax
      If(Fine_STMnm(i)(1:1).eq.'R') Then
         write(2,"(F6.2)", ADVANCE='no') 5*Fine_STDelayB(i)
      EndIf
   Enddo
   write(2,*) '  '
   !
   Do i=1,Station_nrMax
      If(Fine_STMnm(i)(1:1).eq.'R') Then
         If(StationInCalA(i).eq.1 .and. StationInCalB(i).eq.1) Then
            write(2,"(3x,'1*1')", ADVANCE='no')
         Else
            write(2,"(2x,i2,i2)", ADVANCE='no') StationInCalA(i),StationInCalB(i)
         EndIf
      EndIf
   Enddo
   write(2,*) '  '
   !
   Stop
   !
End Program CompareCalibr
