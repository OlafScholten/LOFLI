    Include 'ConstantsModules.f90'
    Include 'ParamModules.f90'
    Include 'MappingUtilities.f90'
!-----------------------------------
!-----------------------------------
Program CompareCalibr
   use constants, only : dp
   use DataConstants, only : Ant_nrMax, Station_nrMax
   use DataConstants, only : Station_nrMax, DataFolder, Calibrations! , Diagnostics
   use Calibration, only : ReadCalib, Fine_STMnm, Fine_STDelay
   Implicit none
   !
   Real(dp) :: Fine_STDelayA(1:Station_nrMax), Fine_STDelayB(1:Station_nrMax)! time calibration data
   Character(len=85) :: CalibrationsA, CalibrationsB
   Integer :: i
   !
   Open(unit=2,STATUS='unknown',ACTION='write', FILE ="CalComp.out")
   Read(*,*)  CalibrationsA
   Read(*,*)  CalibrationsB
   Fine_STMnm(:)=''
   !
   Calibrations=CalibrationsA
   Call ReadCalib()
   Fine_STDelayA(:)=Fine_STDelay(:)
   !
   Calibrations=CalibrationsB
   Call ReadCalib()
   Fine_STDelayB(:)=Fine_STDelay(:)-Fine_STDelayA(:)
   !
   write(2,*) TRIM(CalibrationsB),' minus ',TRIM(CalibrationsA), ' in [ns]'
   write(2,"(40A7)") Fine_STMnm(:)
   write(2,"(40F7.2)") 5*Fine_STDelayB(:)
   Do i=1,Station_nrMax
      If(Fine_STMnm(i)(1:1).eq.'C') Then
         write(2,"(A6)", ADVANCE='no') Fine_STMnm(i)
      EndIf
   Enddo
   write(2,*) '  '
   Do i=1,Station_nrMax
      If(Fine_STMnm(i)(1:1).eq.'C') Then
         write(2,"(F6.2)", ADVANCE='no') 5*Fine_STDelayB(i)
      EndIf
   Enddo
   write(2,*) '  '
   Do i=1,Station_nrMax
      If(Fine_STMnm(i)(1:1).eq.'R') Then
         write(2,"(A6)", ADVANCE='no') Fine_STMnm(i)
      EndIf
   Enddo
   write(2,*) '  '
   Do i=1,Station_nrMax
      If(Fine_STMnm(i)(1:1).eq.'R') Then
         write(2,"(F6.2)", ADVANCE='no') 5*Fine_STDelayB(i)
      EndIf
   Enddo
   write(2,*) '  '
   !
   Stop
   !
End Program CompareCalibr
