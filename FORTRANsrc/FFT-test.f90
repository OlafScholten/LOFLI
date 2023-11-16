!
    Include '../Program/ConstantsModules-v15.f90'
    Include '../Program/ParamModules-v15.f90'
    Include '../Program/FFT_routines-v2.f90'
PROGRAM FFT_test
! copied from C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\Imaging\LOFAR-MAP
   use constants, only : dp,pi,ci, sample
   use DataConstants, only : DataFolder, Time_dim, Cnu_dim
   use FFT, only : RFTransform_su,DAssignFFT, RFTransform_CF, RFTransform_CF2CT, Hann
   IMPLICIT none
   !
   Real(dp) :: RTime_s(1:Time_dim)
   Complex(dp) :: CTime_s(1:Time_dim), CNu_s(0:Cnu_dim)
   integer :: Dset_offset
   Real(dp) :: nu_Fltr(0:Cnu_dim), Av, Bv, FiltFact
   Integer :: nu, nu_i, nu_f, dnu, i
   !
   call RFTransform_su(Time_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !
   OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='FFT-test.out')
   Do i=1,50
      RTime_s(i)=sin(i*0.5)
   Enddo
   Call RFTransform_CF(RTime_s,Cnu_s)
   Call RFTransform_CF2CT(Cnu_s,CTime_s)
   Do i=1,20
      Write(2,*) Cnu_s(i),RTime_s(i),'?=?',CTime_s(i)
   Enddo
   !
   Stop
!--------------------------------------------------
END PROGRAM FFT_test
!===================================================================
!=========================================
!==========================================
! ----------------------------------------------------------------------------------------
!     shellin = 'gle /d pdf '//trim(dummy3(i))//'.gle'
!     CALL system(shellin)
!     shellin = 'epstopdf '//trim(dummy3(i))//'.eps'
!     call system(shellin)
