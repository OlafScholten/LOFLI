!=================================
Module MDD_Pars
   use constants, only : dp
   !use DataConstants, only : Time_Dim  !, Cnu_dim, DataFolder, OutFileLabel
   !use DataConstants, only : Ant_nrMax
   Implicit none
   Integer, save :: MDDtDim, MDDnuDim ! time and frequency dimensions of part of trace used
   Integer, save :: T_Range  ! Length of time window [samples] involved in fitting, should be an odd number
   !
   Integer, parameter :: G_dim=2000
   Real(dp), save :: G_ns(G_dim) ! Impulse response function at [ns] timing to make for simple interpolation
   Real(dp), save :: IRFW_s  ! Full width at half max of amplitude of the Impulse Response in [samples]
   !
   Character(len=9), save :: MDDLabel
   Real(dp), save, Allocatable :: RTime_p(:,:), RTime_t(:,:)
   Real(dp), save, Allocatable :: Vec_p(:,:), Vec_t(:,:), Vec_l(:,:)
   Real(dp), save, Allocatable :: Weight(:),  AntSourceD(:)
   Integer, parameter :: NSrc_max=20  ! Maximum number of DD in composite source
   Real(dp), save ::  CR_Tms, C_Tms, dCenT(1:NSrc_max) ! time difference (in [samples]) of each DD from central time  C_Tms (in [ms]) or time @ref ant (=CR_Tms)
   Real(dp), save :: C_Loc(1:3), dCenloc(1:3,1:NSrc_max) ! location of each DD w.r.t. C_Loc, both in [m]
   Integer, save :: NSrc ! Number of DD sources used
   Real(dp), save :: DelDip(1:3,1:NSrc_max) ! dipole vectors for each of the delta sources
   !
   Real(dp), save, allocatable :: DDChiSQ(:,:)  ! Chi^2 as calculated in CompareMDD per antenna, per polarization
   Real(dp), allocatable :: WORK(:)  ! needed in CompareMDD
   Integer, save :: LWORk ! Length of array WORK
   Real(dp), allocatable :: G_IR(:,:,:)  ! greens functions for antennas and DD sources
   !Real(dp), allocatable :: RE_MDD_p(:), RE_MDD_p(:)
!   contains
! --------------
! --------------
End Module MDD_Pars
