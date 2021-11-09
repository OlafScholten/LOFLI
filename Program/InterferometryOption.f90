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
   Use Interferom_Pars, only : N_pix, d_loc, CenLocPol, CenLoc, PixLoc, IntFer_ant  !, Nr_IntFer
   Use Interferom_Pars, only : Smooth, SumStrt, SumWindw, NrSlices, SliceLen
   Use Interferom_Pars, only : IntfPhaseCheck, IntfSmoothWin, N_smth, IntPowSpec
   use DataConstants, only : OutFileLabel, DataFolder, Calibrations, Time_Dim !, Cnu_dim,
   !use Chunk_AntInfo, only : CTime_spectr, Ant_Stations, Ant_IDs, Ant_nr, Ant_pos, Ant_RawSourceDist
   use DataConstants, only : Polariz
   use Chunk_AntInfo, only : Start_time, TimeFrame, TimeBase, ExcludedStat, NoiseLevel, RefAnt, Simulation
   use ThisSource, only : Dual, CurtainPlot, PeakNrTotal, Peak_eo, Peak_start, Peakpos, SourcePos
   use FitParams, only : AntennaRange
   use constants, only : dp, ci, pi, Sample, Refrac, c_mps
   use FFT, only : RFTransform_su, DAssignFFT !, RFTransform_CF, RFTransform_CF2CT
   !use StationMnemonics, only : Statn_ID2Mnem
   use HDF5_LOFAR_Read, only : CloseDataFiles
   use Chunk_AntInfo, only : SignFlp_SAI, PolFlp_SAI, BadAnt_SAI
   use GLEplots, only : GLEplotControl
   !use mod_test
   Implicit none
   !real ( kind = 8 ) :: Dist
   integer :: i,j, k, i_s, i_loc(1), i_ant, j_corr, i_eo, Station, j_IntFer, i_nu, nxx
   character(len=5) :: txt
   Character(LEN=180) :: lname
   !Real(dp) :: RDist, angle, angle_all, dt_AntPix
   Real(dp) :: SmPow
   !Real(dp) :: RMax, Rmin, Epsln=epsilon(Epsln)
   !Complex(dp) :: dPhase, Phase, A
   Integer :: Date_T(8)

   !integer :: i_N, i_E, i_h !, N_h,  N_hlow, N_Eup, N_Elow
   character*100 :: shellin, IntfRun_file
   Character*7 :: PreAmble
   Logical :: Interferometry=.true., E_FieldsCalc
   NAMELIST /Parameters/ Interferometry, Dual, E_FieldsCalc, CurtainPlot &  ! just those that are of interest for interferometry
         , IntfPhaseCheck, IntfSmoothWin, TimeBase  &
         , SignFlp_SAI, PolFlp_SAI, BadAnt_SAI, Calibrations &
         , ExcludedStat, OutFileLabel, ChainRun, NoiseLevel, AntennaRange    !  ChunkNr_dim,
   !    .
   E_FieldsCalc=Polariz
   N_smth=IntfSmoothWin
   If(N_smth.lt.2) N_smth=20
   Smooth(:)=0.
   Smooth(0)=1.
   Smooth(N_smth/2)=0.5  ! needed for block profile for even values of N_smth
   Do i=1,N_smth
      !Smooth(i)=(1.-(Real(i)/N_smth)**2)**2.5  ! Parabola to power to make it =0.5 at N_smth/2
      If(i.lt.(N_smth+1)/2) Smooth(i)=1.      ! Block gives indistinguishable results from parabola; parabola has a tiny bit smaller volume around peak.
      Smooth(-i)=Smooth(i)
   EndDo
   Write(2,"(1x,A,I3,A,I3,I3,A,3F6.3)") 'Width smoothing fie=', N_smth, &
      '; relative values in range (',N_smth/2-1,N_smth/2+1,') = ',Smooth(N_smth/2-1:N_smth/2+1)/Smooth(0)
   SmPow=SUM(Smooth(:))
   Smooth(:)=Smooth(:)/SmPow
   !
   Call GetNonZeroLine(lname)
   Read(lname,*) StartTime_ms, CenLoc  ! Start time
   write(2,"(A,'"',A,'"')") 'Interferometry input line-1; time & position:',,TRIM(lname)
   StartTime_ms=StartTime_ms+TimeBase
   Call Convert2m(CenLoc)
   !If(abs(CenLoc(1)) .lt. 100. ) CenLoc(:)=CenLoc(:)*1000.  ! convert from [km] to [m]
   write(2,*) 'true start time (adding base)=', StartTime_ms, '[ms]'
   Call GetNonZeroLine(lname)
   Read(lname(2:180),*,iostat=nxx) N_pix(1,2), d_loc(1), N_pix(2,2), d_loc(2), N_pix(3,2), d_loc(3)
   SELECT CASE (lname(1:1))
      CASE("P")  ! polar option is selected
         Polar=.true.
      CASE("C")   ! Cartesian option
         Polar=.false.
      CASE DEFAULT
         Polar=.true.
         If(d_loc(1).gt.0.5) Polar=.false.
   End SELECT
   write(2,"(A,'"',A,'"')") 'Interferometry input line-2; grid:',TRIM(lname)
   !d_loc(1)=d_N  ;  d_loc(2)=d_E   ;  d_loc(3)=d_h
   N_pix(:,1)=-N_pix(:,2)
   If(polar) then
   write(2,*) 'Polar coordinates used for grid!!'
   d_loc(1)=d_loc(1)*pi/180.   ! convert to radian
   d_loc(2)=d_loc(2)*pi/180.   ! convert to radian
   Endif
   !
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
   ! Start with interferometry setup
   !
   ! Summing regions
   Call GetNonZeroLine(lname)
   Read(lname,*,iostat=nxx) SumStrt, SumWindw, NrSlices, AmpltPlot
   write(2,"(A,'"',A,'"')") 'Interferometry input line-3: SumStrt, SumWindw, NrSlices, AmpltPlot=:',TRIM(lname)
   If(nxx.ne.0) then
      write(2,*) 'reading error in 3rd line:',lname
      stop 'Interferometry; reading error'
   Endif
   If(NoiseLevel.lt.0) NoiseLevel=0.2
   !If(NrSlices.gt.NS_max) NrSlices=NS_max
   If(NrSlices.lt.1) NrSlices=1
   SliceLen=NINT(SumWindw*1./NrSlices)
   If(SliceLen.lt.1) SliceLen=1
   SumWindw=NrSlices*SliceLen
   write(2,*) 'Adapted window slices',SumWindw, NrSlices, SliceLen, AmpltPlot
   If(SumWindw .lt. 5*N_smth) IntPowSpec=.false.
   write(2,*) 'SumWindw:',SumWindw,IntPowSpec
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
   If(Dual) Then
      Write(2,*) 'performing E-field Interferometry'
      Call EI_run
      PreAmble='EI'
   Else
      Write(2,*) 'performing even/odd antenna-Signal Interferometry'
      Call SI_run
      PreAmble='Interf'
   EndIf


   DeAllocate( Smooth )
   !
   !------------------------------------------
   !
   !trim(DataFolder)//TRIM(OutFileLabel)//'IntfSpecPowMx'//TRIM(txt)//'.dat'
   !Open(UNIT=10,STATUS='unknown',ACTION='WRITE',FILE='GLE-plots.sh')
   !write(2,*) 'PlotAllCurtainSpectra=',CurtainPlot,SumStrt,SumWindw
   If(CurtainPlot) then
      !read(lname,*,iostat=nxx) StartTime_ms, SourceGuess(:,i_chunk), CurtainWidth !, PeakPos(1), PeakPos(2)  ! Start time offset = 1150[ms] for 2017 event
      PeakNrTotal=2
      Peak_eo(1)=0   ; Peak_eo(2)=1  ! only the i_eo=0 peaks are used for curtainplotting
      Peak_start(1)=1   ; Peak_start(2)=1 ! chunk nr
      RefAnt(1,0)=IntFer_ant(1)
      !RefAnt(i_chunk,i_eo)
      !Unique_StatID(i_stat)  ! stations that will be used ???
      PeakPos(1)=SumStrt+SumWindw/2  ; PeakPos(2)=SumStrt+SumWindw/2  ! take a single peak at the center of the window
      !SourcePos(1,i_Peak)
      Call PlotAllCurtainSpectra(SumWindw/2)
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
         '#!/bin/bash', &
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
