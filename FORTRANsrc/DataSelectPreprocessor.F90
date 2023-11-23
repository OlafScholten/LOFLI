! Update/revision of
!  C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\Imaging\LMA\LMA2019\Imaging\Track-New.f90
!  Selecting sources usung quality indicators
!  Arranging sources in tracks
!  Allow for pre-defined tracks
!  Perform track analysis
!     - lateral spread
!     - velocity distributions
!     - Source-time distributions
!  Output:
!     - RunGLEplots.bat (Windows) contains commands to
!        - to run "Utilities/SourcesPlot.gle"
!        - to run "Utilities/TrackScatt.gle"
!
!
!    Include 'TrackFind.f90'
    Include 'DS_SbRtns.f90'
    Include 'DS_Correlator.f90'
!===========================================================
PROGRAM DataSelect
! short DS
! copied from C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\Imaging\LOFAR-MAP
   Use constants, only : dp, CI, pi, c_l
   Use TrackConstruct, only : PreDefTrackNr, PreDefTrackFile
   Use TrackConstruct, only : TrackNr, LongTrackNr, LongTrack_Min, TrackNrMax, NLongTracksMax, TrackLenMax
   Use TrackConstruct, only : Wtr, MaxTrackDist, TimeWin, HeightFact, dt_MTL, Aweight, TrackENr, TrackE, TrackNrLim, Aweight
   Use DS_Select, only : Image, BckgrFile, datafile, ZoomBox, RMS_ns
   Use DS_Select, only : RA, maxd, SourcTotNr, Label, t_offset, tCutl, tCutu, xyztBB
   Use DS_Select, only : DS_Mode, MaxAmplFitPercent, SMPowCut
   Use DS_Select, only : Corr_dD, Corr_Dnr, Corr_dtau, QualPlot, AmplScale,  AmplitudePlot
   Use DS_Select, only : Nrm, a1,b1,c1,d1,ChiSq1,Nrm1, a2,b2,c2,d2,ChiSq2,Nrm2, a3,b3,c3,d3,ChiSq3 ! parameters for describing the amplitude distribution
   use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, DataFolder, FlashName, Windows, RunMode
   use GLEplots, only : GLEplotControl
   Use DS_Select, only : DS_ReadCntrl
   Use unque, only : HPSORT_mult_RI !,  WriDir
   Use AmpFit, only : AmplitudeFit,AmplitudeFitTracks
   Use TrackConstruct, only : Assign2Tracks, ConstructTracks, BinTracks, AnalyzeBurst
   IMPLICIT none
   character*150 :: SystemCommand, SelFileName, PlotFile
   Character*20 :: Utility, release
   Character(len=180) :: lineTXT
   integer :: i,j,k, nxx, d_Ampl, Fini, SuccessImages  ! unit,
   integer :: kk, jk, wrunt, N_EffAnt !, Max_EffAnt
   Logical :: EndInput
   ! For preprocessor options, see
   !   http://ahamodel.uib.no/intel/GUID-F6619F6F-7D70-4B06-A266-6F39EF6D51B7.html
   !  https://cyber.dabamos.de/programming/modernfortran/preprocessor.html
   !  https://fortranwiki.org/fortran/show/Preprocessors
   ! The FPP preprocessor:
   !  https://www.cita.utoronto.ca/~merz/intel_f10b/main_for/mergedProjects/bldaps_for/common/bldaps_use_fpp.htm
   !  https://www.smcm.iqfr.csic.es/docs/intel/compiler_f/main_for/bldaps_for/common/bldaps_use_ffpdir.htm
   ! Info on make and cmake:
   !  https://cyber.dabamos.de/programming/modernfortran/build-automation.html
   ! Suggestions for source documentation:
   !  https://cyber.dabamos.de/programming/modernfortran/source-code-documentation.html
#if defined(WINDOWS)
character(20), parameter :: ccOS = "cOS_WIN"
#elif   ! defined(LINUX)
character(20), parameter :: ccOS = "cOS_LINUX"
#endif
   Real :: Qual
  !
  OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='DataSelect.out')
  write(2,*) 'Program DataSelectPreprocessor', ccOS
  !    ,', limited version with source statistics disabled because of FFT calls needing library.'
  Utility='Track&Qual.Control'
  release='v23'
  Call System_Initiation(Utility, release, ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows)  ! set Flash name & folder names
  RunMode=5  !  used in the gle-plotting routine
  SuccessImages=0  ! Number of succesful Images (more than 1 point)
  !
  !
   Do
      ! Read in selection criteria from stdin and apply to located sources
      Call DS_ReadCntrl(Fini)
      If(Fini.ne.0) exit
      ! Read in selection criteria from stdin and apply to located sources
      Call DS_ReadSelData
      !
      write(*,"(A,A,A,I6)") '# of sources ',trim(Image),' =',SourcTotNr
      write(2,"(A,A,A,I6)") '# of sources ',trim(Image),' =',SourcTotNr
      If(SourcTotNr.le.1) then
        write(2,*)'SourcTotNr=',SourcTotNr
        write(2,*) '======================================= next figure'
        cycle  !  Read data for next plot
      endif
      !
      SuccessImages=SuccessImages+1
      !Sorting data
      CALL HPSORT_mult_RI(RA,Label,SourcTotNr)  ! RA(1:4,*) are the accepted source locatopns
      ! write(2,"(A,2f12.5,A,f12.3)") 'Time span data=', RA(1,1), RA(1,SourcTotNr), ', SMPowCut=', SMPowCut
      If(tCutu.gt.tCutl) Then
         write(2,*) 't-notch filter:', tCutl, tCutu
      EndIf
      ! write(*,*) 'Sorting Done'
      !
      !-------------- Start putting sources on tracks ----------------------
      If(NLongTracksMax.Gt.0) Then
         If(TimeWin .lt. 0.) TimeWin = ((RA(1,SourcTotNr)-RA(1,1))*10./SourcTotNr) ! *((RA(1,SourcTotNr)-RA(1,1))*10./SourcTotNr)
         ! 2017:  0.03 0.3 0.01  ! MaxTrackDist[km], Wtr[0.5], TimeWin[ms^2]
         ! 2018:  0.12 0.5 0.05  ! MaxTrackDist[km], Wtr[0.5], TimeWin[ms^2]
         If(NLongTracksMax.lt.0) NLongTracksMax=0
         If(NLongTracksMax.gt.9) NLongTracksMax=0
         TrackNrLim=NLongTracksMax
         write(2,"(A,F7.4,A,F4.1,A,F5.2,A,F11.8,A,i2)") 'MaxTrackDist=',MaxTrackDist,'[km], HeightFact=',HeightFact, &
            '; Weight new event=', Wtr, ', TimeWin=', TimeWin,'[ms], Max Nr Long Tracks=', NLongTracksMax
         !call flush(2)
         !
         LongTrackNr=0
         If(NLongTracksMax.gt.0) then
            Call Assign2Tracks(RA, Label, SourcTotNr)
         EndIf
         ! TrackE(k,i)=j :  Source TrackE(k,i) is the k^th source assigned to track i
         ! TrackENr(i) : Total number of sources assigned to track i
         If(LongTrackNr.gt.0) then
            PlotFile=TRIM(DataFolder)//trim(Image)
            Write(2,"(A,A,A,I6)") 'PlotFile: ',trim(PlotFile)
            Flush(unit=2)
            Call ConstructTracks(RA, Label, SourcTotNr, PlotFile)
            !
            Call AmplitudeFitTracks(TrackNr,TrackLenMax,TrackNrMax,TrackENR,TrackE)
            write(2,*) 'AmplitudeFitTracks(AmplitudeFit) was called!!!!!!'
            !
            If (dt_MTL.gt.0.) then
               Call BinTracks(dt_MTL, RA, SourcTotNr, PlotFile)  ! t_Mean Track Location
            EndIf
               !
            If (dt_MTL.gt.0.) then   ! dt_MTL may be set to negative in "BinTracks" when too few bins
               Call AnalyzeBurst(RA, Label, SourcTotNr, PlotFile) ! Analyze how many sources occur within a certain time-resolution
               Call GLEplotControl(PlotType='TrackScatt', PlotName=TRIM(Image)//'TrSc', &
                  PlotDataFile=TRIM(PlotFile), Submit=.false.)
               !SystemCommand="call GLE -d pdf -o TrSc_"//trim(datfile)//".pdf ../Utilities/TrackScatt.gle "//trim(pars)
               !      CALL SYSTEM(SystemCommand,stat)
            EndIf
         EndIf
      EndIf  ! (NLongTracksMax.Gt.0)
      !
      !
      If(MaxAmplFitPercent.gt.0.) Then
         SelFileName=TRIM(DataFolder)//'AmplFit'//trim(Image)
         Call AmplitudeFit(SourcTotNr, SelFileName=SelFileName)
         !  SelFileName set to '' when there are too few amplitudes to meke a reasonable fit
         !write(2,"(1x,A, f6.2, A, 2g10.4, g11.3, A,2g11.3)") '$b *exp(-a*A-c/A^2); \chi^2=$',ChiSq1, &
         !                        ',with  a,b,c= ',a1,b1/Nrm,c1,'; nrm=',Nrm1, Nrm
         !write(2,"(1x,A, f6.2, A, 2g10.4, g11.3, A,2g11.3)") '$b *A^-a *exp(-c/A); \chi^2=$',ChiSq2, &
         !                        ',with  a,b,c= ',a2,b2/Nrm,c2,'; nrm=',Nrm2, Nrm
         !
         If(SelFileName.NE. '') Then
            Call GLEplotControl(PlotType='Intensity', PlotName=trim(Image)//'AmplFit', &
                  PlotDataFile=trim(SelFileName), Submit=.false.)
         EndIf
      EndIf ! (MaxAmplFitPercent.gt.0.)
      !
      !  Produce Images of sources
      If(DS_Mode .eq. 1) Then  ! Impulsive imager data
         SelFileName=TRIM(DataFolder)//trim(Image)
         PlotFile=TRIM(Image)
         Qual=RMS_ns
      Else  ! TRI-D Imager
         SelFileName=TRIM(DataFolder)//trim(Image)
         PlotFile=trim(Image)
         Qual=SMPowCut
         ! like  H2f#+Mx_dHIntfSpecSel.pdf  -->  H2f#+HIntfSpecSel.pdf  -->  H2f#+H_ISS.pdf
      EndIf
      write(2,"(A,2f12.5,A,f12.3)") 'Time span data=', RA(1,1), RA(1,SourcTotNr), ', QualityCut=', Qual
      OPEN(unit=28,FILE=TRIM(SelFileName)//'.dat', FORM='FORMATTED',STATUS='unknown',ACTION='write') ! space separated values
      write(28,"(6F9.3,2F9.3,1x,A,1x,F7.1,' 0 ')") xyztBB(:), trim(ZoomBox), t_offset         ! gleRead:  xMin xMax yMin yMax zMin zMax tMin tMax ZmBx$ t_offst
      write(28,"(1x,I2,I8,F11.2, ' 0 ',2x,A,1x,F6.2, 6(1x,g10.4),f7.3,1x,A)")  &
         LongTrackNr, SourcTotNr, Qual, trim(Image), AmplitudePlot,  &
         a1,b1,c1,a2,b2,c2,d2, ' ! by DataSelect' ! gleRead:  NTracks SourcTotNr Q t_start label$ AmplitudePlot a1 b1 c1 a2 b2 c2 d2
      write(28,"('! SourceLab',T15,'time',T29,'x=East',T46,'y=North',T61,'z=h',T80,'Ampl',T87,'Wl',T92,'Wr')")
      !
      Do j=1,SourcTotNr  !  finally write all accepted source to file
         write(28,"(1x,i8,4(2x,g14.8),1x,F12.2,2x,I3,2x,I3)")  Label(1,j),RA(1:4,j), Label(2,j)/AmplScale
      enddo
      close(unit=28)
      write(2,*) 'testing trim(ZoomBox):', trim(ZoomBox)
      write(2,*) 'testing LongTrackNr, SourcTotNr, Qual,AmplitudePlot:', LongTrackNr, SourcTotNr, Qual,AmplitudePlot
      write(2,*) SourcTotNr,' sources in plots:',trim(Image),' ; to file: ',TRIM(SelFileName)//'.dat'
      !
      If(BckgrFile.ne.'') then
         Call GLEplotControl(PlotType='SourcesPlotBckgr', PlotName=trim(PlotFile), &
               PlotDataFile=trim(SelFileName), Submit=.false., Bckgr=TRIM(DataFolder)//trim(BckgrFile) )
      Else
         Call GLEplotControl(PlotType='SourcesPlot', PlotName=trim(PlotFile), &
               PlotDataFile=trim(SelFileName), Submit=.false.)
      EndIf
!Used as  $gle -d pdf -o IntfMx_dH2t3.pdf ${UtilDir}SourcesPlot.gle ${FlashFolder}/files/H2t3IntfSpecPowMx_d
      !  Corr (=correlate sources following Frankfurt AI paper)
      If((Corr_dD.gt.0) .and. (Corr_Dnr .gt. 0)) Then
            Call ApplyCorrelator(RA, Label, SourcTotNr, Aweight, Corr_dD, Corr_Dnr, Corr_dtau)
      EndIf
      !  Qual (=analysis of quality factors, make 2D scatter plots)
      If(QualPlot) Then
            Call ApplyQualityAna
      endif
      !
998   write(2,*) '======================================='
   enddo !  end of the main reading loop
   !
999 continue
   If(SuccessImages.gt.0) Then
      Flush(Unit=2)
      Call GLEplotControl(Submit=.true.)
   EndIf
End PROGRAM DataSelect
!--------------------------------------------------
Subroutine ApplyQualityAna()  ! du -sh  directory size
   Use constants, only : dp
   use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, DataFolder, FlashName, Windows, RunMode
   Use DS_Select, only : Image !,  WriDir
   use GLEplots, only : GLEplotControl
   Use DS_Select, only : SourcTotNr, Label, QualIndic
   IMPLICIT none
   Integer :: i, i_src
   !
         Write(2,*) 'option: ','Qual'
   !
   OPEN(unit=29,FILE=TRIM(DataFolder)//trim(Image)//'QuAna.dat',FORM='FORMATTED',STATUS='unknown')  ! will contain all accepted sources
   write(29,"(2x,A,1x,I6,1x,F6.3,i5,' !')") trim(Image), SourcTotNr
   Write(2,*) 'Quality Indicator analysis, writing ',trim(Image), SourcTotNr
   Do i_src=1,SourcTotNr
      write(29,*) Label(1:4,i_src), QualIndic(1:4,i_src)
   EndDo ! i_src=1,SourcTotNr-1
   close(unit=29)
   Write(2,*) 'ApplyQualityAna; File:',TRIM(DataFolder)//trim(Image)//'QuAna.dat',' with ',SourcTotNr, ' data lines'
   !
   Call GLEplotControl(PlotType='QualContrPlot', PlotName=TRIM(Image)//'QC', &
         PlotDataFile=TRIM(DataFolder)//trim(Image), Submit=.false.)
End Subroutine ApplyQualityAna
!-----------------------------------------------------
! ----------------------------------------------------------------------------------------
!     shellin = 'gle /d pdf '//trim(dummy3(i))//'.gle'
!     CALL system(shellin)
!     shellin = 'epstopdf '//trim(dummy3(i))//'.eps'
!     call system(shellin)
