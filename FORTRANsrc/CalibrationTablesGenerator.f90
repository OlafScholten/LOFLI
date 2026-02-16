    module constants
        integer, parameter:: dp=kind(0.d0)
        COMPLEX(dp), parameter :: CI=CMPLX(0._dp,1._dp)
        REAL(dp), parameter :: pi=2.d0*asin(1.d0)
        REAL(dp), save :: c_l=0.299792458d0     ! speed of light in Vacuum [m/nsec]
      real(dp), save :: RotMat(3,3) ! Rotation matrix from ITRF to local Superterp coordinates
      real(dp), parameter :: center_CS002(3)= (/ 3826577.5, 461021.3125, 5064893.0 /)
      real(dp), parameter :: sample=5.d-9   ! Sample length in seconds
    end module constants
!-----------------------------------------
program CTG
!   CalibrationTablesGenerator.f90
    use constants, only : dp,pi,ci
    Implicit none
    Character(len=20) :: release
    INTEGER :: DATE_T(8),i
    CHARACTER*12 :: REAL_C(3)
    !
    Integer :: nxx
    character(len=5) :: txt, Station(45)
    Real(dp) :: Delay, DelayX(45), DelayY(45), DelayX002, DelayY002
    CHARACTER(LEN=140) :: line
    Character(len=31) :: Key
    character(len=2) :: core
    Integer :: Nr_cs,Nr_tot
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='CalibrationTablesGenerator.out')
    release='Most Recent'
    write(2,"(3x,5(1H-),1x,'CalibrationTablesGenerator release of ',A22,25(1H-))") release
    CALL DATE_AND_TIME (REAL_C(1),REAL_C(2),REAL_C(3),DATE_T)
    WRITE(2,"(3X,5(1H-),1x,'run on ',I2,'/',I2,'/',I4,' , started at ',&
          I2,':',I2,':',I2,'.',I3,1X,25(1H-))") &
          DATE_T(3),DATE_T(2),DATE_T(1),(DATE_T(i),i=5,8)
    !
    Open(UNIT=9,STATUS='old',ACTION='READ',FILE='StationCalibration.parset') ! LBA.LBA_OUTER.LBA_10_90.phase0
    Open(UNIT=10,STATUS='old',ACTION='READ',FILE='station_delays_2.txt')
    Open(unit=13,STATUS='unknown',ACTION='WRITE', FILE = 'StationCalibrations.dat')
    i=0
    nxx=1
    Do !while(nxx.ne.0)
        read(9,*,IOSTAT=nxx) line, txt, Delay
        !write(2,*) delay,txt,line
        if(nxx.ne.0) exit
! PIC.Core.CS006LBA.LBA_OUTER.LBA_10_90.phase0.X               = 0.000000e+00
        Read(line,"(9x,A5,A32)") txt,Key
        !write(2,*) txt,' , "',key,'"',nxx
        If(Key.eq.'LBA.LBA_OUTER.LBA_10_90.delay.X') then
            write(*,*) txt,i,delay
            i=i+1
            DelayX(i)=Delay
            Station(i)=txt
            read(9,*) line, txt, Delay
            DelayY(i)=Delay
        endif
    enddo
    Nr_cs=i
    delay=DelayX(2)
    Do i=1,Nr_cs
        DelayX(i)=Delay-DelayX(i)
        DelayY(i)=Delay-DelayY(i)
    Enddo
    i=Nr_cs
    Do
        read(10,*,IOSTAT=nxx) txt,delay
        if(nxx.ne.0) exit
        read(txt,*) core
        If(core.eq.'RS') then
            i=i+1
            DelayX(i)=Delay
            Station(i)=txt
            DelayY(i)=Delay
        endif
    Enddo
    Nr_tot=i
    write(2,*) Nr_cs,Nr_tot
    Do i=1,Nr_tot
        write(2,*) i,Station(i),DelayX(i),DelayY(i)
        write(13,*) Station(i),DelayX(i),DelayY(i)
    Enddo
    Close(unit=9)
    Close(unit=10)
    Close(unit=13)
    stop

end program CTG
!=================================
