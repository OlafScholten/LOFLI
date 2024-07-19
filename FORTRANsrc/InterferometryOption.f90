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
   Use Interferom_Pars, only : Polar, NewCenLoc, ChainRun,  StartTime_ms, AmpltPlot, i_chunk
   Use Interferom_Pars, only : N_pix, d_loc, CenLocPol, CenLoc, PixLoc, IntFer_ant
   Use Interferom_Pars, only : SumStrt, SumWindw, NrSlices, SliceLen
   Use Interferom_Pars, only : IntfPhaseCheck, N_smth, IntPowSpec, NrPixSmPowTr
   !Use Interferom_Pars, only : PixPowOpt
   use DataConstants, only : OutFileLabel, DataFolder, Time_Dim !, Cnu_dim,
   !use Chunk_AntInfo, only : CTime_spectr, Ant_Stations, Ant_IDs, Ant_nr, Ant_pos, Ant_RawSourceDist
   !use DataConstants, only : Polariz
   use Chunk_AntInfo, only : StartT_sam, TimeFrame,   NoiseLevel, RefAnt, Simulation
   use ThisSource, only : CurtainHalfWidth, PeakNrTotal, Peak_eo, ChunkNr, Peakpos, SourcePos
   !use FitParams, only : AntennaRange
   use constants, only : dp, ci, pi, Sample, c_mps, sample_ms
   use FFT, only : RFTransform_su, DAssignFFT !, RFTransform_CF, RFTransform_CF2CT
   !use StationMnemonics, only : Statn_ID2Mnem
   use HDF5_LOFAR_Read, only : CloseDataFiles
   use GLEplots, only : GLEplotControl
   use PredefinedTracks, only : PreDefTrackFile
   !use LOFLI_Input, only : ReadSourceTimeLoc
   !use mod_test
   Implicit none
   integer :: i,j, k, i_s, i_loc(1), i_ant, j_corr, i_eo, Station, j_IntFer, i_nu,  nxx
   Integer, parameter ::lnameLen=180
   character(len=5) :: txt
   character(Len=1) :: Mark
   Character(LEN=lnameLen) :: lname
   !Real(dp) :: SmPow,  MidTime_ms
   Integer :: Date_T(8)
   Character*7 :: PreAmble
   Real(dp), external :: tShift_ms
   !    .
   ! Get center voxel
   Call ReadSourceTimeLoc(StartTime_ms, CenLoc)  ! routine resides in "LOFLI_InputHandling.f90
   !  A predefined track name can be entered.
   !  This may be a .dat file having the same format at image-source files, i.e. label, t, x y z
   !        or a .trc file with format   t, N, E, h
   !
   ! Get grid:
   Call GetMarkedLine(Mark,lname)
   write(2,"(A)") 'Input line-2: "'//Mark//'|'//TRIM(lname)// &
      '" !  Polar(Phi,Th,R)/Carthesian(N,E,h) | 3x(#gridpoints, grid spacing)'
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
   If(SumStrt.lt.1000) Then
      write(2,*) 'Window starting sample should be 1000'
      Stop 'Window starting sample should be sufficiently large'
   EndIf

   If(PreDefTrackFile.ne.'') Then  ! time was mid-time on track
      StartTime_ms=StartTime_ms-sample_ms*(SumStrt+SumWindw/2)  !Middle of segmet in local time
   EndIf

   write(2,"(20(1x,'='))")
   !
   !  Read the time traces
   CALL DATE_AND_TIME (Values=DATE_T)
   WRITE(2,"(1X,I2,':',I2.2,':',I2.2,'.',I3.3,A)") (DATE_T(i),i=5,8), 'initializing'
   i_chunk=1
   StartT_sam(i_chunk)=(StartTime_ms/1000.d0)/sample  ! in sample's
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
   write(2,*) 'NewCenLoc', NewCenLoc, ' , ChainRun=', ChainRun
   If(ChainRun.ne.0) Call ChainRuns(NewCenLoc)
   Return
    !
End Subroutine InterferometerRun
!-------------------------------------------------------------
Subroutine ChainRuns(NewCenLoc)
! Run several TRI-D calculations in a series
! If no TrackFile is specified, the new center location is taken as the bary-center of the present run (=ewCenLoc)
! With a TrackFile there are two options,
!     t_track is not given (=negative):
!        Inch along the track to see the sources on the track with the TRI-D imager
!     t_track is given (=positive):
!        Inch along the track to see possible sources on the track with the TRI-D imager that is
!            run at a fixed time as given by StartTime_ms.
   Use Interferom_Pars, only : Polar, ChainRun,  StartTime_ms, AmpltPlot
   Use Interferom_Pars, only : N_pix, d_loc, CenLocPol, CenLoc
   Use Interferom_Pars, only : SumStrt, SumWindw, N_smth
   !Use Interferom_Pars, only : IntfPhaseCheck, N_smth, IntPowSpec, NrPixSmPowTr
   Use Interferom_Pars, only : PixPowOpt
   use DataConstants, only : OutFileLabel, DataFolder, Calibrations !, Cnu_dim,
   !use Chunk_AntInfo, only : CTime_spectr, Ant_Stations, Ant_IDs, Ant_nr, Ant_pos, Ant_RawSourceDist
   !use DataConstants, only : Polariz
   use Chunk_AntInfo, only : TimeBase, ExcludedStat, NoiseLevel,  Simulation, SaturatedSamplesMax
   !use ThisSource, only : CurtainHalfWidth,  Peakpos, SourcePos
   use FitParams, only : AntennaRange
   use constants, only : dp, ci, pi, Sample, c_mps, sample_ms
   !use FFT, only : RFTransform_su, DAssignFFT !, RFTransform_CF, RFTransform_CF2CT
   !use StationMnemonics, only : Statn_ID2Mnem
   !use HDF5_LOFAR_Read, only : CloseDataFiles
   use Chunk_AntInfo, only : SignFlp_SAI, PolFlp_SAI, BadAnt_SAI
   use PredefinedTracks, only : GetTrPos, GetTrTime, PreDefTrackFile, t_track
   Implicit none
   Real(dp), intent(inout) :: NewCenLoc(1:3)
   integer :: i,j, k, i_s, i_loc(1), i_ant, j_corr, i_eo, Station, j_IntFer, i_nu, IntfSmoothWin, nxx
   Integer, parameter ::lnameLen=180
   character(len=5) :: txt
   character(Len=1) :: Mark
   Character(LEN=lnameLen) :: lname
   Real(dp) :: t_shft, MidTime_ms, t_new , StartTime_local
   character*100 :: shellin, IntfRun_file
   Character*5 :: RunOption='TRI-D' !, E_FieldsCalc
   Real(dp), external :: tShift_ms
   Logical :: FollowCenterIntensity=.false.
   NAMELIST /Parameters/ RunOption &!, E_FieldsCalc  &  ! just those that are of interest for interferometry
         , IntfSmoothWin, TimeBase, PixPowOpt  &
         , SignFlp_SAI, PolFlp_SAI, BadAnt_SAI, Calibrations, SaturatedSamplesMax &
         , ExcludedStat, OutFileLabel, ChainRun, NoiseLevel, AntennaRange    !  ChunkNr_dim,
   !
   If(ChainRun.eq.0) Return
   write(2,*) StartTime_ms, t_track
   t_shft=tShift_ms(CenLoc(:)) ! sqrt(SUM(CenLoc(:)*CenLoc(:)))*1000.*Refrac/c_mps ! in mili seconds due to signal travel distance
   StartTime_ms=StartTime_ms-t_shft-TimeBase
   StartTime_ms=StartTime_ms+sample_ms*(SumStrt+SumWindw/2)  !Middle of segmet in local time
   If(t_track.le.0.) Then  ! Move along a track when present
      If(ChainRun.gt.0) Then
         ChainRun=ChainRun-1
         StartTime_ms=StartTime_ms+sample_ms*SumWindw
         !write(txt,"('#+',I2.2)") ChainRun
         txt='#+01'
      Else
         ChainRun=ChainRun+1
         StartTime_ms=StartTime_ms-sample_ms*SumWindw
         !write(txt,"('#-',I2.2)") -ChainRun
         txt='#-01'
      EndIf
      If(PreDefTrackFile.ne.'') Then  ! With a track defined StartTime_ms is actually midtime
         Call  GetTrPos(StartTime_ms, NewCenLoc)
         write(2,*) 'NewCenLoc from track', StartTime_ms, NewCenLoc
      Else If(.not.  FollowCenterIntensity) Then
         NewCenLoc(:)=CenLoc(:) ! Stick to the old position
         write(2,*) 'old position for image cube kept @',CenLoc(:)
      EndIf ! otherwise take the NewCenLoc as has been calculated (=centroid of strength in present run)
      t_new=-1.
   Else ! move image window along track but keep time fixed (used for looking for positive leaders)
      If(ChainRun.le.0) Then
         write(2,*) 'Chainrun should be positive for the "Follow a track at fixed time" (FTFT) option '
         Stop 'Negative Chainrun with FTFT'
      Endif
      ChainRun=ChainRun-1
      Call GetTrTime( t_track, t_new, NewCenLoc)
      write(2,*) 't_track(old)=',t_track,', t_track(new)=',t_new,', NewCenLoc from track =',  NewCenLoc
      txt='#+01'
   EndIf
   !
   !
   If(NewCenLoc(1).ne.0.d0) then
      !write(2,*) 'txt1',txt,ChainRun
      j=LEN(TRIM(OutFileLabel))
      k=j
      Do i=1,j
         If(OutFileLabel(i:i) .eq. '#') Then
            k=i-1
            read(OutFileLabel(k+3:k+4),*) i_s
            i_s=i_s+1
            write(OutFileLabel(k+3:k+4),"(I2.2)") i_s
            exit
         Endif
      Enddo
      If(k.eq.j) OutFileLabel(k+1:k+4)=txt(1:4)
      IntfRun_file='A-Intf_'//TRIM(OutFileLabel)
      write(2,*) 'OutFileLabel: ',TRIM(OutFileLabel),', ChainRun=', ChainRun, ', IntfRun_file:', IntfRun_file
      Open(UNIT=10,STATUS='unknown',ACTION='WRITE',FILE=TRIM(IntfRun_file)//'.in' )
      !
      IntfSmoothWin=N_smth
      write(10,NML = Parameters)
      !
      StartTime_ms=StartTime_ms + t_shft - tShift_ms(NewCenLoc(:))  ! correct for local time shift, keeping same time in the core
      If(PreDefTrackFile.ne.'') Then
         Write(10,"('S',F12.6,1x, A, F12.6)") StartTime_ms, TRIM(PreDefTrackFile), t_new ! Start time at the source location
      Else
         Write(10,"('S',F12.6, 3F10.4, F6.1)") StartTime_ms, NewCenLoc/1000. ! Start time at the source location
         !  239.7 , -25.42  37.28   5.54,  50 !-i- PID=25216;  StartTime_ms, (N,E,h)CenLoc, DistMax[km]
      EndIf
      lname(1:1)="C"
      If(Polar) then
         lname(1:1)="P"
         d_loc(1)=d_loc(1)*180./pi   ! convert to degree
         d_loc(2)=d_loc(2)*180./pi   ! convert to degree
      EndIf
      write(10,"(A,3(I5,F9.5))") lname(1:1),N_pix(1,2), d_loc(1), N_pix(2,2), d_loc(2), N_pix(3,2), d_loc(3)
!P  40 .003, 15 .01 , 20 10. ! N_phi, d_N, N_theta, d_E, N_R, d_h[m]  ; VeryFine resolution extended grid(vfe)
      Write(10,"('F',2I7,F7.3,F9.4)")  SumStrt, SumWindw,  AmpltPlot
!   2000  60040 10 0.01 0.1 ! Ini [samples] starting loc & window for summing power & slice number & powerthresh & PlotAmplitude
      Close(Unit=10)
      Open(UNIT=10,STATUS='unknown',ACTION='WRITE',FILE=TRIM(IntfRun_file)//'.sh' )
      Write(10,"(9(A,/) )") &
         '#!/bin/bash', &  !   -v', &
         '#','#', &
         'source  ${LL_Base}/ShortCuts.sh' ,   &  ! defines ${ProgramDir} and ${FlashFolder}
!         'source ${UtilDir}/compile.sh',  &
         'cp ${LL_bin}/LOFAR-Imag ./Imag-'//TRIM(OutFileLabel)//'.exe',  &  ! Just to display a more intelligent name when runnung 'ps' or 'top'
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
End Subroutine ChainRuns
!=========================================
