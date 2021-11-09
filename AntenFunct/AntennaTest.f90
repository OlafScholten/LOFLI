    Include '../Program/ConstantsModules-v18.f90'
    Include '../Program/FFT_routines-v2.f90'
    Include '../Program/ParamModules-v20.f90'  ! v18d: Cnu storage changed  for InterfEngineB
    Include '../Program/AntFunct-v20.f90'
    !
!-----------------------------------------------
Program AntennaTest
   ! CMCnu is the frequency dependent Current Moment Contribution of each antenna pair
   ! Needs:
   !   ?
   !--------------------------------------------
   use constants, only : dp, pi, Sample, Refrac, c_mps
   use AntFunCconst, only : Freq_min, Freq_max
   Use AntFunCconst, only : AntTst
   Implicit none
   integer :: i_ant, j_IntFer, MLoc, Mloc_all, i, j, i_s, i_freq, i_nu, IntfDim
   !Real(dp) :: angle, angle_all, HorAng, AziAng, D, RD,  AbsAmp, NRM, SmPow  ! elevation
   !Complex(dp) :: CNu_s(0:Cnu_dim)  ! CTime_s(1:Time_dim),
   !
   !
   ! Get antenna function parameters
   open(unit=2, status='unknown', file='AntennaTest.out')
   AntTst=.true.
   Call AntFieParGen()
   !
   Return
End Program AntennaTest
!-----------------------------------------------
! =========================
