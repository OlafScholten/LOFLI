!=================================
Module MDD_Pars
   use constants, only : dp
   !use DataConstants, only : Time_Dim  !, Cnu_dim, DataFolder, OutFileLabel
   !use DataConstants, only : Ant_nrMax
   Implicit none
   Integer, save :: MDDtDim, MDDnuDim ! time and frequency dimensions of part of trace used
   Integer, save :: T_Range  ! Length of time window [samples] involved in fitting, should be an odd number
   Integer, save :: Trace_DeadEnds=25 ! Begin and Ends of T_range that should be void of DD sources
   Integer, save :: RealDip_logic=0 ! 0: never set RealDip=.t.; 1: RealDip=.t. in gridsearch;
                                    ! 2: RealDip=.t. in final;  3: all RealDip=.t.;
   Integer, save :: NSrc_min=15  !  The number of sources retained in the final iterations
   Real(dp), save :: GridScale=1.d0 ! Scale factor for grid
   Real(dp), save :: IntensityFrac_keep=0.3d0 ! the fraction of the max intensity below which the sources are thrown out
   Real(dp), save :: DepleteCoreWeigth=1.0d0 ! Factor with which the weights of the CS stations is decreased
   Real(dp), save :: DepleteDistFrac=1.0d0  ! Keep weights constant for antennas that are closer to the source than DepleteDistFrac*(distance source-core)
   !
   Integer, parameter :: G_dim=2000
   Complex(dp), save :: G_ns(G_dim) ! Impulse response function at [ns] timing to make for simple interpolation
   Real(dp), save :: IRFW_s  ! Full width at half max of amplitude of the Impulse Response in [samples]
   !
   Character(len=9), save :: MDDLabel
   Complex(dp), save, Allocatable :: CTime_p(:,:), CTime_t(:,:)
   Real(dp), save, Allocatable :: Vec_p(:,:), Vec_t(:,:), Vec_l(:,:)
   Real(dp), save, Allocatable :: Weight(:), SQWght(:),  AntSourceD(:)
   Real(dp), save :: WeightNorm
   Real(dp), save :: SpaceGrid(1:3)
   Integer :: NSrc_max=20  ! Maximum number of DD in composite source
   Real(dp), save ::  CR_Tms, C_Tms !  central time  C_Tms (in [ms]) or time @ref ant (=CR_Tms)
   Real(dp), save :: C_Loc(1:3)! Central location C_Loc,  in [m]
   Real(dp), save, allocatable ::  dCenT(:) ! time difference (in [samples]) of each DD from central time  C_Tms
   Real(dp), save, allocatable :: dCenloc(:,:) ! location of each DD w.r.t. C_Loc,  in [m]
   Integer, save :: NSrc ! Number of DD sources used
   Complex(dp), save, allocatable :: DelDip(:,:) ! dipole vectors for each of the delta sources
   Logical, save :: RealDip=.false.
   !
   ! In the DD fitter:
   Real(dp), save, allocatable :: DDChiSQ(:,:)  ! Chi^2 as calculated in CompareMDD per antenna, per polarization
   Complex(dp), allocatable :: WORK(:)  ! needed in CompareMDD for Call ZHESV
   Real(dp), allocatable :: ReWORK(:)  ! needed in CompareMDD for Call dsysv
   Integer, save :: LWORK ! Length of array WORK
   Real(dp), parameter :: Xscale=1.d7  ! Factor to compensate automatic setting of first step
   Real(dp), save :: Bias_inv=1.d-8 , Bias_inv_base=1.d-8  ! Bias added to the diagonal matrix elements of A in "CompareMDD" before inversion
   Real(dp), save :: DDXcorr_depl=0.0d0 ! Depletion factor for off-diagonal cross correlations, works only when <1.
   Complex(dp), allocatable :: G_IR(:,:,:)  ! greens functions for antennas and DD sources
   !
   Integer, Parameter :: NTestSt_max=10
   Character(len=5), save :: Teststations_MN(1:NTestSt_max)
   Integer, save :: TeststAntenna_ID(1:NTestSt_max), NTestAnt
!   contains
! --------------
! --------------
End Module MDD_Pars
