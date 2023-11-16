!=================================
!================================
!=================================
Subroutine MultipleDDSourcesOption()
!   Chi-square interferometric imaging of peaks found by the impulsive imager
!     Fit the spectra in all antennas as the superposition of a number of dipole-delta sources
!     A dipole-delta (DD) sources in this context is a point-impulsive source (zero extent in space and time) emitting like an ideal dipole.
!
!  Input: ??
!     -reads a list of sorces?
!
!  Output:
!     - list of dipole-delta parameters for each source
!
!   Procedural steps:
!     1) read source parameters
!     2) determine the range in time and space over whic DD sorces should be fit
!
!Copy here from  Subroutine InterferometerRun  (in InterferometryOption.f90)
   use constants, only : dp, ci, pi, Sample, sample_ms
   Use MDD_Pars, only : CR_Tms, C_Tms, C_Loc, T_Range, MDDtDim, MDDnuDim, MDDLabel  ! output
   Use MDD_Pars, only : Trace_DeadEnds ! input
   Use Interferom_Pars, only : StartTime_ms, CenLoc, IntFer_ant, i_chunk, Nr_IntFerCh
   use DataConstants, only : OutFileLabel, DataFolder, Time_Dim !, Cnu_dim,
   use Chunk_AntInfo, only : StartT_sam, RefAnt, Simulation, TimeBase
   use FFT, only : RFTransform_su, DAssignFFT !, RFTransform_CF, RFTransform_CF2CT
   use HDF5_LOFAR_Read, only : CloseDataFiles
   use GLEplots, only : GLEplotControl
   !use HDF5_LOFAR_Read, only : Group_names, DSet_names, Group_nr ! , Group_max
   !use HDF5_LOFAR_Read, only : DSet_nr, Ant_ID, STATION_ID !, DSet_max, ANTENNA_POSITION
   !use HDF5_LOFAR_Read, only : DATA_LENGTH, SAMPLE_NUMBER_first, Absolute_TIME, DIPOLE_CALIBRATION_DELAY
   use HDF5_LOFAR_Read, only : Absolute_TIME
   !use HDF5_LOFAR_Read, only : GetFileName, ListGroups, ListDataAtt, ListGroupStructure
   Implicit none
   integer :: i, error, peaknr, nxx  !   j, k, i_s, i_loc(1), i_ant, j_corr, i_eo, Station, j_IntFer, i_nu,  nxx
   Real(dp) :: FitQual, t_shft, ERA
   Integer :: Range, Nr_IntFer
   Integer, parameter ::lnameLen=180
   !character(len=5) :: txt
   character(Len=1) :: Mark
   Character(LEN=lnameLen) :: lname
   Integer :: Date_T(8)
   !Character*7 :: PreAmble
   Real(dp), external :: tShift_ms
   Logical :: Beamforming=.true.
   Complex(dp), allocatable :: CMCnu(:,:,:)
   !    .
   ! Get center voxel
   Call ReadSourceTimeLoc(StartTime_ms, CenLoc)  ! general area
   !
   !  Read the time traces
   CALL DATE_AND_TIME (Values=DATE_T)
   WRITE(2,"(1X,I2,':',I2.2,':',I2.2,'.',I3.3,A)") (DATE_T(i),i=5,8), 'initializing'
   i_chunk=1
   StartT_sam(i_chunk)=(StartTime_ms/1000.d0)/sample  ! in sample's
   Call RFTransform_su(Time_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   Call AntennaRead(i_chunk, CenLoc)
   !
   ! Get time of the recording
   Date_t(1)=Absolute_TIME/(365*24*60*60)  ! year
   Date_t(1)=(Absolute_TIME-(1451606400))/(24*60*60)   ! days since GMT Jan 1  2016
   !write(2,*) 'days since 1 jan 2016=', Date_t(1)
   Date_t(2)=Date_t(1)/365  ! year since 2016
   Date_t(3)=Date_t(1)-365*Date_t(2) - (date_t(2)+3)/4  ! days since jan 1 current year, corrected for leap years
   Date_t(4)=(Absolute_TIME-(1451606400))-(24*60*60)*Date_t(1) ! seconds since midnight
   Date_t(6)=Date_t(3) ! Day of the month
   If(Date_t(3).gt. (58+(date_t(2)+4)/4) ) Date_t(6)=Date_t(6) + 3-(date_t(2)+4)/4 ! Correct for february
   If(Date_t(6).lt.154) Then
      Date_t(6)=Date_t(6)- 31*(Date_t(6)/31) + Date_t(6)/61  ! Day of the month
      Date_t(5)=(Date_t(3)-Date_t(6))/29  ! month of the year
   Else
      Date_t(6)=Date_t(6)-154  ! Days after July 31
      i=Date_t(6)- 31*(Date_t(6)/31) + Date_t(6)/61  ! Day of the month
      Date_t(5)=7+(Date_t(6)-i)/30  ! month of the year
      Date_t(6)=i
   EndIf
   !
   write(2,*) 'TIME of recording: Year=',2016+Date_t(2),', day=',Date_t(3)+1, &
      ', Month=',Date_t(5)+1, ', Day=', Date_t(6)+1,', hour=', Date_t(4)/(60*60.)
   !
   ERA=0.7790572732640d0 + 1.00273781191135448d0 * (Absolute_TIME-946728000.d0)/(24*60*60)
   i=FLOOR(ERA +19.2 -0.779)
   write(2,"(A, F5.3,A,F5.1,A)") 'ERA-fraction (=fraction of siderial day)', ERA-FLOOR(ERA) &
      , ', LST=', 24*(ERA +19.2 -0.779-i),' [h]'
   ! $t_U$ is time from 12:00 (midday) Terrestrial Time of January 1, 2000, Unix Timestamp = 946728000.
   !
   If(Simulation.eq."") Then
      CALL CloseDataFiles()    ! no more reading of antenna data, only close when reading real data
   EndIf
   call DAssignFFT()                   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !
   Call EISelectAntennas(i_chunk)  ! select antennas for which there is an even and an odd one.
   Nr_IntFer=Nr_IntFerCh(i_chunk)
   Call AntFieParGen()  !  Prepare antenna function and gain
   !Call EffAntSr()  ;   stop
   Call ImpResp_ns()    !  Construct Impulse response function
   !
   peaknr=0
   Do
      Do
         Call GetMarkedLine(Mark,lname)
         SELECT CASE (Mark)
            CASE("T")  ! Read Test Stations for timing info
               Call ReadTeststations(lname)
            CASE("&")  ! Read some key steering parameters for MDD mode
               Call ReadMDD_C(lname)
            CASE DEFAULT  ! Continue calculations
               Exit
         End SELECT
      EndDo
      read(lname,*,IOSTAT=nxx) C_Tms, C_Loc(1:3), Range !  full range in [s] that should be considered
      If(nxx.ne.0) exit
      peaknr=peaknr+1
      Call Convert2m(C_Loc(:))
      C_Tms=C_Tms+TimeBase
      t_shft=tShift_ms(C_Loc(:)) ! sqrt(SUM(CenLoc(:)*CenLoc(:)))*1000.*Refrac/c_mps ! in mili seconds due to signal travel distance
      CR_Tms=C_Tms+t_shft  ! time in reference antenna
      Write(2,*) '============= CR_Tms:',CR_Tms, C_Tms, peaknr,'===================='
      write(MDDLabel,"(A,I1)") 'MDD',peaknr
      !
      T_Range=(Range/2)*2 +1  ! make sure it is an odd number
      If(Trace_DeadEnds*2 .gt. T_Range) Then
         write(2,*) '******  Trace_DeadEnds*2 exceeds T_Range', Trace_DeadEnds*2, T_Range
         cycle
      EndIf
      Call WriteMDD_C()
      !
      MDDtDim=50+2*T_Range ! lead+tail =window lendth length + 50 for safety,
      MDDnuDim=MDDtDim/2  ! Length of frequency trace used for interference
      i = shiftr( MDDnuDim, 9 )+1
      MDDnuDim=shiftl( i, 9 )  ! Make sure it approches some multiple of 2
      !write(2,"(A,o10,o10,i7)") 'i:',i,IntfNuDim,IntfNuDim
      MDDtDim=2*MDDnuDim ! Length of time trace used for interference
      !write(*,"(A,o10,b20)") '(lead,dim):', IntfLead, IntfDim
      !
      Allocate( CMCnu(Nr_IntFerCh(i_chunk),0:MDDnuDim,1:3)  ) ! MDDnuDim=IntfNuDim
      Call MDDSetup(Nr_IntFer, CMCnu)
      !
      !Beamforming=.false.
      If(Beamforming) Then
         Call EI_scan(Nr_IntFer, MDDnuDim, CMCnu)
         DeAllocate( CMCnu  )       !
      Else
         DeAllocate( CMCnu  )       !
         Call MDD_IniManySrcs()
         !
         !Call MDD_QuickFit(FitQual,error)
         !write(2,*) 'Values without time & position search:'
         !Call MDD_write()
         !Call MDD_ResortSrcs()
         !
         Call MDD_EliminateSrcs()
      EndIF
      !
      Call MDD_Fitter(FitQual,error)
      !Call MDD_write()
      Call WriteMDD_C()
      Call MDD_GLEplot(peaknr)
      !
      Call MDDClose()
   EndDo
   Return
End Subroutine MultipleDDSourcesOption
!=================================
!=================================

!============================
