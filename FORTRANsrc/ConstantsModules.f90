Module constants
    integer, parameter:: dp=kind(0.d0)
    COMPLEX(dp), parameter :: CI=CMPLX(0._dp,1._dp)
    REAL(dp), parameter :: pi=2.d0*asin(1.d0)
    REAL(dp), parameter :: c_l=0.299792458d0     ! speed of light in Vacuum [m/nsec]
    REAL(dp), parameter :: c_mps=c_l * 1.d9 ! =0.299792458d9     ! speed of light in Vacuum [m/sec]
    !REAL(dp), parameter :: Refrac=1.0003     ! set in routine 'AveIndxRefr'
    real(dp), save :: RotMat(3,3) ! Rotation matrix from ITRF to local Superterp coordinates
    real(dp), parameter :: center_CS002(3)= (/ 3826577.5, 461021.3125, 5064893.0 /)
    real(dp), parameter :: sample=5.d-9   ! Sample length in seconds [s/sample]
    real(dp), parameter :: sample_ms=5.d-6   ! Sample length in miliseconds
    !Logical, save :: HeightCorrectIndxRef=.false.
    Logical, save :: HeightCorrectIndxRef=.true.
End Module constants
!-----------------------------------
Module DataConstants
    Integer, parameter :: Time_dim=  65536 ! 65536=2^16 ; 32768=2^15 ! 2048=2^11
    Integer, parameter :: EdgeOffset=11000 ! about 15km distance margin
    !Integer, parameter :: Time_dim=32768 ! 32768=2^15 ! 2048=2^11
    !nteger, parameter :: EdgeOffset=7000 !=5km/c /sample=5 10^3 /(3 10^8 5 10^-9)=5 10^3/1.5=7000
    Integer, parameter :: Cnu_dim=Time_dim/2
    Integer, parameter :: Ant_nrMax=1500
    Integer, parameter :: Station_nrMax=40
    Integer, save :: PeakNr_dim=25 ! number of peaks while calibrating (set on input)
    Integer, save :: ChunkNr_dim=1 ! Number of starting times for blocks of data that are used
    !Integer, parameter :: ChunkNr_dim=1 ! Number of starting times for blocks of data that are used
    Character(len=25), save :: DataFolder='files/'  ! points to the subfolder storing produced data (plotting etc)
    Character(len=100), save :: UtilitiesFolder='' !  path to the directory containing utilities (like the .gle files)
    Character(len=100), save :: ProgramFolder='' !  path to the directory containing the fortran programs
    Character(len=300), save :: FlashFolder='' ! Absolute path to the working (=flash) directory
    Character(len=40), save :: FlashName='' ! Short name for this flash, like 20A-1, set equal to the name of the folder
    Character(len=25), save :: OutFileLabel=''
    Character(len=85), save :: Calibrations=''  ! Name of the file containing the calibration data
    Character(len=20), save :: Utility, release
    Logical, save :: Diagnostics=.false.  !  Print diagnostic info on LOFAR-data
    Logical, save :: Production=.false.
    Logical, save :: Windows=.false.  ! Check if run on a windows machine
    Integer, save :: RunMode=0 ! 1=Explore; 2=Calibrate; 3=Impulsive imager; 4=Interferometry
    !Logical, save :: Polariz=.false.
    Logical, save :: Interferometry=.false.
End Module DataConstants
! --------------
