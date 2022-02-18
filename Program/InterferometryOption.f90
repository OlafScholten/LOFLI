    Include 'InterferometryOptSbRtns.f90'
    !
Subroutine InterferometerRun
    !+++++++++ What is used:
    !- Loop over antennas
    ! - make frequency spectra
    !- End antennaLoop
    !- Loop over all pixels
    !  - calculate phase shift Each antenna
    !   - Loop frequency
    !    - Sum(antenna) spectra
    ! - End pixelloop
    !+++++++++ Alternative:
    !- Loop over antennas
    ! - make frequency spectra
    ! - Loop over all pixels
    !    - calculate phase shift
    !    - Sum frequency: update spectra with this antenna
    ! - End pixelLoop
    !- End antennaLoop
    !++++++++++++++++++++++
   Use Interferom_Pars, only : Polar, NewCenLoc, ChainRun, StartTime_ms, AmpltPlot, i_chunk
   Use Interferom_Pars, only : N_pix, d_loc, CenLocPol, CenLoc, PixLoc, IntFer_ant
   Use Interferom_Pars, only : SumStrt, SumWindw, NrSlices, SliceLen
   Use Interferom_Pars, only : IntfPhaseCheck, N_smth, IntPowSpec, NrPixSmPowTr
   use DataConstants, only : OutFileLabel, DataFolder, Calibrations, Time_Dim !, Cnu_dim,
   !use Chunk_AntInfo, only : CTime_spectr, Ant_Stations, Ant_IDs, Ant_nr, Ant_pos, Ant_RawSourceDist
   use DataConstants, only : Polariz
   use Chunk_AntInfo, only : Start_time, TimeFrame, TimeBase, ExcludedStat, NoiseLevel, RefAnt, Simulation
   use ThisSource, only : Dual, CurtainHalfWidth, PeakNrTotal, Peak_eo, ChunkNr, Peakpos, SourcePos
   use FitParams, only : AntennaRange
   use constants, only : dp, ci, pi, Sample, Refrac, c_mps
   use FFT, only : RFTransform_su, DAssignFFT !, RFTransform_CF, RFTransform_CF2CT
   !use StationMnemonics, only : Statn_ID2Mnem
   use HDF5_LOFAR_Read, only : CloseDataFiles
   use Chunk_AntInfo, only : SignFlp_SAI, PolFlp_SAI, BadAnt_SAI
   use GLEplots, only : GLEplotControl
    use LOFLI_Input, only : ReadSourceTimeLoc
   !use mod_test
   Implicit none
   integer :: i,j, k, i_s, i_loc(1), i_ant, j_corr, i_eo, Station, j_IntFer, i_nu, IntfSmoothWin, nxx
   Integer, parameter ::lnameLen=180
   character(len=5) :: txt
   character(Len=1) :: Mark
   Character(LEN=lnameLen) :: lname
   Real(dp) :: SmPow, t_shft
   Integer :: Date_T(8)
   character*100 :: shellin, IntfRun_file
   Character*7 :: PreAmble
   Logical :: Interferometry=.true., E_FieldsCalc
   NAMELIST /Parameters/ Interferometry, Dual, E_FieldsCalc  &  ! just those that are of interest for interferometry
         , IntfPhaseCheck, IntfSmoothWin, TimeBase  &
         , SignFlp_SAI, PolFlp_SAI, BadAnt_SAI, Calibrations &
         , ExcludedStat, OutFileLabel, ChainRun, NoiseLevel, AntennaRange    !  ChunkNr_dim,
   !    .
   E_FieldsCalc=Polariz
   IntfSmoothWin=N_smth
   !
   ! Get center voxel
   Call ReadSourceTimeLoc(StartTime_ms, CenLoc)
   !Call GetNonZeroLine(lname)
   !Read(lname(2:lnameLen),*) StartTime_ms, CenLoc  ! Start time
   !write(2,"(A)") 'Interferometry input line-1: "'//lname(1:1)//'|'//TRIM(lname(2:lnameLen))// &
   !   '" !  Source/Reference-| time, & position'
   !flush(unit=2)
   !Call Convert2m(CenLoc)
   !StartTime_ms=StartTime_ms+TimeBase
   !t_shft=sqrt(SUM(CenLoc(:)*CenLoc(:)))*1000.*Refrac/c_mps ! in mili seconds due to signal travel distance
   !j = iachar(lname(1:1))  ! convert to upper case if not already
   !if (j>= iachar("a") .and. j<=iachar("z") ) then
   !   lname(1:1) = achar(j-32)
   !end if
   !SELECT CASE (lname(1:1))
   !   CASE("S")  ! time at Source (central voxel) is given
   !      StartTime_ms=StartTime_ms+t_shft
   !   CASE DEFAULT  ! time at reference antenna is given
   !End SELECT
   !write(2,"(A,F12.6,A,F12.6,A)") &
   !   ' Ttrue start time trace, adding base, in ref antenna (at source)=', StartTime_ms, ' (',StartTime_ms-t_shft,') [ms]'
   !
   ! Get grid:
   Call GetMarkedLine(Mark,lname)
   write(2,"(A)") 'Input line-2: "'//Mark//'|'//TRIM(lname)// &
      '" !  Polar(Phi,Th,R)/Carthesian(N,E,h) | 3x(grid spacing, #gridpoints)'
   Read(lname,*,iostat=nxx) N_pix(1,2), d_loc(1), N_pix(2,2), d_loc(2), N_pix(3,2), d_loc(3)
   If(nxx.ne.0) then
      write(2,*) 'reading error in 2nd line:'
      stop 'Interferometry; reading error'
   Endif
   SELECT CASE (Mark(1:1))
      CASE("P")  ! polar option is selected
         Polar=.true.
      CASE("C")   ! Cartesian option
         Polar=.false.
      CASE DEFAULT
         Polar=.true.
         If(d_loc(1).gt.0.5) Polar=.false.
   End SELECT
   !d_loc(1)=d_N  ;  d_loc(2)=d_E   ;  d_loc(3)=d_h
   N_pix(:,1)=-N_pix(:,2)
   If(polar) then
   write(2,*) 'Polar coordinates used for grid!!'
   d_loc(1)=d_loc(1)*pi/180.   ! convert to radian
   d_loc(2)=d_loc(2)*pi/180.   ! convert to radian
   Endif
   !
   ! Start with interferometry setup
   !
   ! Summing regions
   Call GetMarkedLine(Mark,lname)
   !Read(lname,*,iostat=nxx) SumStrt, SumWindw, NrSlices, AmpltPlot
   Read(lname,*,iostat=nxx) SumStrt, SumWindw, AmpltPlot ! rough slicing by NrSlices, is depreciated, Dec 2021
   write(2,"(A)") 'Input line-3: "'//Mark//'|'//TRIM(lname)//'"  !  First/Median| SumStrt, SumWindw, AmpltPlot'
   If(nxx.ne.0) then
      write(2,*) 'reading error in 3rd line:'
      stop 'Interferometry; reading error'
   Endif
   NrSlices=1  !  depreciated, Dec 2021; option may be removed completely someday
   !If(NrSlices.gt.NS_max) NrSlices=NS_max ;  NrPixSmPowTr=(SumWindw-1)/N_smth-1
   !If(NrSlices.lt.1) NrSlices=1  ! This is for the rough slicing, still implemented but depreciated
   If(SumWindw .lt. 3*N_smth) Then
      SumWindw=3*N_smth
      write(2,*) 'SumWindw was smaller than 3 times slicing window of',N_smth
   EndIf
   NrPixSmPowTr=SumWindw/N_smth-1  ! The number of fine slices for TRI-D imaging
   SumWindw=N_smth*(NrPixSmPowTr+1) ! set to an integer multiple of N_smth and allow for initial and final trailing
   If(Mark.eq.'M') Then
      If(MOD(NrPixSmPowTr,2).eq.0) NrPixSmPowTr=NrPixSmPowTr+1 ! Make total number of slices=NrPixSmPowTr odd
      SumWindw=N_smth*(NrPixSmPowTr+1) ! set to an integer multiple of N_smth
      SumStrt=SumStrt-SumWindw/2-1 ! make sure there is a slice centered around this sample
   EndIf
   SliceLen=SumWindw/NrSlices ! obsolete
   If(SliceLen.lt.1) SliceLen=1 ! obsolete
   NrSlices=SumWindw/SliceLen ! obsolete
   write(2,"(A,I5,A,I5,A,F6.3)") 'Adjusted TRI-D window, SumStrt=',SumStrt, ', SumWindw=', SumWindw, ', AmpltPlot=', AmpltPlot
   write(2,"(20(1x,'='))")
   !
   !  Read the time traces
   CALL DATE_AND_TIME (Values=DATE_T)
   WRITE(2,"(1X,I2,':',I2.2,':',I2.2,'.',I3.3,A)") (DATE_T(i),i=5,8), 'initializing'
   i_chunk=1
   Start_time(i_chunk)=(StartTime_ms/1000.)/sample  ! in sample's
   TimeFrame=i_chunk
   Call RFTransform_su(Time_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   Call AntennaRead(i_chunk, CenLoc)
   If(Simulation.eq."") Then
      CALL CloseDataFiles()    ! no more reading of antenna data
   EndIf
   call DAssignFFT()                   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !
   If(NoiseLevel.lt.0) NoiseLevel=0.2
   !
   ! ----------------------------------------------------
   ! Extract polar coordinates;  1=N=phi; 2=E=theta; 3=h=Radial distance
   ! Check limits on ranges      N_hlow, N_Eup, N_Elow
   write(2,*) 'Central Carthesian coordinates (N,E,h)=',CenLoc(:)
   If(polar) then
      Call Carth2Pol(CenLoc,CenLocPol)
      !CenLocPol(3)=sqrt(SUM(CenLoc(:)*CenLoc(:)))  ! distance [m]    ;h
      !CenLocPol(2)=asin(CenLoc(3)/CenLocPol(3))    ! elevation angle=theta [radian]  ;E range: 0^o< th <80^o
      !CenLocPol(1)=atan2(CenLoc(2),CenLoc(1))      ! phi [radian]    ;N
      write(2,*) 'Central polar coordinates (R,th,ph)=',CenLocPol(3),CenLocPol(2)*180./pi,CenLocPol(1)*180/pi
      ! check ranges:
      If(CenLocPol(3)+N_pix(3,1)*d_loc(3) .lt. 0) N_pix(3,1)= 1-CenLocPol(3)/d_loc(3) ! negative number generally ! 1-: not to have it ridiculously close
      If(CenLocPol(2)+N_pix(2,1)*d_loc(2) .lt. 0) N_pix(2,1)= -CenLocPol(2)/d_loc(2) ! negative number generally
      If(CenLocPol(2)+N_pix(2,2)*d_loc(2) .gt. pi*4/9.) N_pix(2,2)= (pi*4/9.-CenLocPol(2))/d_loc(2) ! pos number generally
      write(2,*) 'distance range:',CenLocPol(3)+N_pix(3,1)*d_loc(3),CenLocPol(3)+N_pix(3,2)*d_loc(3)
      write(2,*) 'Elevation angle range:', (CenLocPol(2)+N_pix(2,1)*d_loc(2))*180./pi , (CenLocPol(2)+N_pix(2,2)*d_loc(2))*180./pi
      write(2,*) 'Azimuth angle range:', (CenLocPol(1)+N_pix(1,1)*d_loc(1))*180./pi, (CenLocPol(1)+N_pix(1,2)*d_loc(1))*180./pi
   Else
      If(CenLoc(3)+N_pix(3,1)*d_loc(3) .lt. 0) N_pix(3,1)= -CenLoc(3)/d_loc(3) ! negative number generally
   Endif
   !N_pix(2,1)=N_Elow   ;  N_pix(3,1)=N_hlow  ;  N_pix(2,2)=N_Eup
   write(2,*) 'N_hlow, N_Eup, N_Elow',N_pix(:,1), N_pix(:,2)
   !
   !trim(DataFolder)//TRIM(OutFileLabel)//'IntfSpecPowMx'//TRIM(txt)//'.dat'
   !------------------------------------
   ! main loop over direction antennas
!   If(Dual) Then
      Write(2,*) 'performing E-field Interferometry'
      Call EI_run
      PreAmble='EI'
!   Else
!      Write(2,*) 'performing even/odd antenna-Signal Interferometry'
!      Call SI_run
!      PreAmble='Interf'
!   EndIf


   !DeAllocate( Smooth )
   !
   !------------------------------------------
   !
   !trim(DataFolder)//TRIM(OutFileLabel)//'IntfSpecPowMx'//TRIM(txt)//'.dat'
   !Open(UNIT=10,STATUS='unknown',ACTION='WRITE',FILE='GLE-plots.sh')
   !write(2,*) 'PlotAllCurtainSpectra=',CurtainPlot,SumStrt,SumWindw
   If(CurtainHalfWidth.gt.0) then
      !read(lname,*,iostat=nxx) StartTime_ms, SourceGuess(:,i_chunk), CurtainWidth !, PeakPos(1), PeakPos(2)  ! Start time offset = 1150[ms] for 2017 event
      write(2,*) 'PlotAllCurtainSpectra:',CurtainHalfWidth
      Flush(unit=2)
      PeakNrTotal=2
      Peak_eo(1)=0   ; Peak_eo(2)=1  ! only the i_eo=0 peaks are used for curtainplotting
      ChunkNr(1)=1   ; ChunkNr(2)=1 ! chunk nr
      RefAnt(1,0)=IntFer_ant(1,1)
      !RefAnt(i_chunk,i_eo)
      !Unique_StatID(i_stat)  ! stations that will be used ???
      PeakPos(1)=SumStrt+SumWindw/2  ; PeakPos(2)=SumStrt+SumWindw/2  ! take a single peak at the center of the window
      !SourcePos(1,i_Peak)
      Call PlotAllCurtainSpectra(CurtainHalfWidth)
   EndIf
   !
   If(NrSlices.gt.1 .and. Polar) Call GLEplotControl(PlotType=TRIM(PreAmble)//'Track', &
            PlotName=TRIM(PreAmble)//'Track'//TRIM(OutFileLabel), PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel) )
   Call GLEplotControl(PlotType=TRIM(PreAmble)//'Contour', PlotName=TRIM(PreAmble)//'Contour'//TRIM(OutFileLabel), &
            PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel), Submit=.true.)
   !
   write(2,*) 'NewCenLoc', NewCenLoc
   If((ChainRun.ne.0).and.(NewCenLoc(1).ne.0.d0)) then
      If(ChainRun.gt.0) Then
         ChainRun=ChainRun-1
         StartTime_ms=StartTime_ms+0.3
         write(txt,"(':+',I2.2)") ChainRun
      Else
         ChainRun=ChainRun+1
         StartTime_ms=StartTime_ms-0.3
         write(txt,"(':-',I2.2)") -ChainRun
      EndIf
      !write(2,*) 'txt1',txt,ChainRun
      j=LEN(TRIM(OutFileLabel))
      k=j
      Do i=1,j
         If(OutFileLabel(i:i) .eq. ':') Then
            k=i-1
            exit
         Endif
      Enddo
      If(i_nu.eq.0) k=1
      OutFileLabel(k+1:k+4)=txt(1:4)
      IntfRun_file='A-Intf_'//TRIM(OutFileLabel)
      write(2,*) 'jk',j, OutFileLabel, i_nu, txt, ChainRun, IntfRun_file
      Open(UNIT=10,STATUS='unknown',ACTION='WRITE',FILE=TRIM(IntfRun_file)//'.in' )
      write(10,NML = Parameters)
      !
      Write(10,"(F9.3, 3F10.4, F6.1)") StartTime_ms-TimeBase, NewCenLoc/1000. ! Start time
!  239.7 , -25.42  37.28   5.54,  50 !-i- PID=25216;  StartTime_ms, (N,E,h)CenLoc, DistMax[km]
      lname(1:1)="C"
      If(Polar) then
         lname(1:1)="P"
         d_loc(1)=d_loc(1)*180./pi   ! convert to degree
         d_loc(2)=d_loc(2)*180./pi   ! convert to degree
      EndIf
      write(10,"(A,3(I5,F9.5))") lname(1:1),N_pix(1,2), d_loc(1), N_pix(2,2), d_loc(2), N_pix(3,2), d_loc(3)
!P  40 .003, 15 .01 , 20 10. ! N_phi, d_N, N_theta, d_E, N_R, d_h[m]  ; VeryFine resolution extended grid(vfe)
      Write(10,"(2I7,I6,F7.3,F9.4)") SumStrt, SumWindw, NrSlices, AmpltPlot
!   2000  60040 10 0.01 0.1 ! Ini [samples] starting loc & window for summing power & slice number & powerthresh & PlotAmplitude
      Close(Unit=10)
      Open(UNIT=10,STATUS='unknown',ACTION='WRITE',FILE=TRIM(IntfRun_file)//'.sh' )
      Write(10,"(9(A,/) )") &
         '#!/bin/bash -v', &
         '#','#', &
         'source ../ShortCuts.sh' ,   &  ! defines ${ProgramDir} and ${FlashFolder}
!         'source ${UtilDir}/compile.sh',  &
         'cp ${ProgramDir}LOFAR-Imag-v20 ./Imag-'//TRIM(OutFileLabel)//'.exe',  &  ! Just to display a more intelligent name when runnung 'ps' or 'top'
         './Imag-'//TRIM(OutFileLabel)//'.exe  ${FlashFolder} <'//TRIM(IntfRun_file)//'.in',  &
         'rm Imag-'//TRIM(OutFileLabel)//'.exe'
      Close(unit=10)
      shellin = 'chmod 755 '//TRIM(IntfRun_file)//'.sh'
      CALL system(shellin)
      shellin = 'nohup ./'//TRIM(IntfRun_file)//'.sh  >'//TRIM(IntfRun_file)//'.log 2>&1  & '
      write(2,*) shellin
      write(*,*) TRIM(IntfRun_file)//' submitted'
      CALL system(shellin)
   EndIf
   Return
    !
End Subroutine InterferometerRun
