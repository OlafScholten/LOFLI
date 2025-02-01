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
!    Include 'DTrackFind.f90'
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
   Use DS_Select, only : Image, BckgrFile, datafile, ZoomBox, RMS_ns, IntensSpread
   Use DS_Select, only : RA, maxd, SourcTotNr, Label, QualIndic, t_offset, tCutl, tCutu, xyztBB, NEhtBB
   Use DS_Select, only : SrcWidth, Iperm, SrcChi2, SrcI20, SrcI3, SrcUn, SrcLin, SrcCirc, SrcPZen, SrcPAzi, SrcISpr
   Use DS_Select, only : SrcI20_r
   Use DS_Select, only : Stk_NEh, PolarAna, IPerm
   Use DS_Select, only : DS_Mode, MaxAmplFitPercent, SMPowCut
   Use DS_Select, only : Corr_dD, Corr_Dnr, Corr_dtau, QualPlot,  AmplitudePlot
   Use DS_Select, only : Nrm, a1,b1,c1,d1,ChiSq1,Nrm1, a2,b2,c2,d2,ChiSq2,Nrm2, a3,b3,c3,d3,ChiSq3 ! parameters for describing the amplitude distribution
   use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, DataFolder, FlashName, Windows, RunMode
   use GLEplots, only : GLEplotControl
   Use DS_Select, only : DS_ReadCntrl
   Use unque, only : SortPerm !, HPSORT_mult_RI, sort
   Use AmpFit, only : AmplitudeFit,AmplitudeFitTracks
   Use TrackConstruct, only : Assign2Tracks, ConstructTracks, BinTracks, AnalyzeBurst
   IMPLICIT none
   character*150 :: SystemCommand, SelFileName, PlotFile
   Character*20 :: Utility, release
   Character*10 :: extension
   Character(len=180) :: lineTXT
   integer :: i,j,k, nxx, d_Ampl, Fini, SuccessImages  ! unit,
   integer :: kk, jk, wrunt, N_EffAnt !, Max_EffAnt
   Logical :: EndInput
   Character*45 :: Qual
   Real(dp) :: Lin1, PolN, PolE, Polh
   !
   OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='DataSelect.out')
   write(2,*) 'Program DataSelect'
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
         DeAllocate(SrcI20)
        goto 998  !  Read data for next plot
      endif
      !
      SuccessImages=SuccessImages+1
      !Sorting data
 !     CALL HPSORT_mult_RI(RA,Label,SourcTotNr)  ! RA(1:4,*) are the accepted source locatopns
      Call SortPerm(RA,SourcTotNr,IPerm)
      ! Sorted has index _s; before sorting index _u
      ! Now:   RA_s(i) = RA_u(Iperm(i))
      !
      Allocate( SrcI20_r(1:SourcTotNr) )
      Do i=1,SourcTotNr
         SrcI20_r(i)=SrcI20(IPerm(i))
      EndDo
      DeAllocate(SrcI20)
      ! write(2,"(A,2f12.5,A,f12.3)") 'Time span data=', RA(1,1), RA(1,SourcTotNr), ', SMPowCut=', SMPowCut
      If(tCutu.gt.tCutl) Then
         write(2,*) 't-notch filter:', tCutl, tCutu
      EndIf
      ! write(*,*) 'Sorting Done'
      !
      !-------------- Start putting sources on tracks ----------------------
      LongTrackNr=0
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
         If(NLongTracksMax.gt.0) then
            Call Assign2Tracks(RA, SrcI20_r, SourcTotNr)
         EndIf
         ! TrackE(k,i)=j :  Source TrackE(k,i) is the k^th source assigned to track i
         ! TrackENr(i) : Total number of sources assigned to track i
         If(LongTrackNr.gt.0) then
            PlotFile=TRIM(DataFolder)//trim(Image)
            Write(2,"(A,A,A,I6)") 'PlotFile: ',trim(PlotFile)
            Flush(unit=2)
            Call ConstructTracks(RA, SrcI20_r, IPerm, SourcTotNr, PlotFile, PolarAna, Stk_NEh)
            !
            Call AmplitudeFitTracks(TrackNr,TrackLenMax,TrackNrMax,TrackENR,TrackE)
            write(2,*) 'AmplitudeFitTracks(AmplitudeFit) was called!!!!!!'
            !
            If (dt_MTL.gt.0.) then
               Call BinTracks(dt_MTL, RA, SourcTotNr, PlotFile)  ! t_Mean Track Location
            EndIf
               !
            If (dt_MTL.gt.0.) then   ! dt_MTL may be set to negative in "\emph{BinTracks}" when too few bins
               Call AnalyzeBurst(RA, SrcI20_r, SourcTotNr, PlotFile) ! Analyze how many sources occur within a certain time-resolution
               Call GLEplotControl(PlotType='SrcsTrScatt', PlotName=TRIM(Image)//'TrSc', &
                  PlotDataFile=TRIM(PlotFile), Submit=.false.)
               If(PolarAna) Call GLEplotControl(PlotType='SrcsTrAngles', PlotName=TRIM(Image)//'Angl', &
                  PlotDataFile=TRIM(PlotFile), Submit=.false.)
               !SystemCommand="call GLE -d pdf -o TrSc_"//trim(datfile)//".pdf ../Utilities/TrackScatt.gle "//trim(pars)
               !      CALL SYSTEM(SystemCommand,stat)
            EndIf
         EndIf
      EndIf  ! (NLongTracksMax.Gt.0)
      !
      If(PolarAna) Call PolarizationAna()
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
         SelFileName=TRIM(DataFolder)//trim(Image)
         PlotFile=TRIM(Image)
      If(DS_Mode .eq. 1) Then  ! Impulsive imager data
         write(Qual,"(A,F5.2,A)") '"RMS<',RMS_ns,' ns"'
      ElseIf(DS_Mode .eq. 2) Then  ! TRI-D Imager
         write(Qual,"(A,F8.1,A)") '"I_{12}>',SMPowCut,',"'
         ! like  H2f#+Mx_dHIntfSpecSel.pdf  -->  H2f#+HIntfSpecSel.pdf  -->  H2f#+H_ISS.pdf
      Else  ! ATRI-D Imager
         write(Qual,"(A,F5.2,'ns & I>',F8.1,' & ISpr<',F8.1,A)") '"Q<',RMS_ns, SMPowCut, IntensSpread,',"'
         ! like  H2f#+Mx_dHIntfSpecSel.pdf  -->  H2f#+HIntfSpecSel.pdf  -->  H2f#+H_ISS.pdf
      EndIf
      write(2,"(A,2f12.5,A,A)") 'Time span data=', RA(1,1), RA(1,SourcTotNr), ', QualityCut:', trim(Qual)
      !
      If(DS_Mode .eq. 3) Then  ! ATRID
         OPEN(unit=28,FILE=TRIM(SelFileName)//'.plt', FORM='FORMATTED',STATUS='unknown',ACTION='write') ! space separated values
         write(28,"(6F9.3,2F9.3,1x,A,1x,F7.1,' 0 ')") NEhtBB(:), trim(ZoomBox), t_offset         ! gleRead:  xMin xMax yMin yMax zMin zMax tMin tMax ZmBx$ t_offst
         write(28,"(1x,I2,I8,1x,A, ' 0 ',2x,A,1x,F6.2, 6(1x,g10.4),f7.3,1x,A)")  &
            LongTrackNr, SourcTotNr, Trim(Qual), trim(Image), AmplitudePlot,  &
            a1,b1,c1,a2,b2,c2,d2, ' ! by DataSelect' ! gleRead:  NTracks SourcTotNr Q t_start label$ AmplitudePlot a1 b1 c1 a2 b2 c2 d2
         write(28,"(2A)") '! Src_lab  t_src[ms]          Position (N,E,h)[km]           Int  Width  Chi^2  Sprd', &
            ' %Un- %lin- %circpol  I_lin   Pol_direction(NEh)'
         !
         Do j=1,SourcTotNr  !  finally write all accepted source to file
            Lin1= SrcI20_r(j)* SrcLin(Iperm(j))/100.
            polN=sin(SrcPZen(Iperm(j))*pi/180.) * cos(SrcPAzi(Iperm(j))*pi/180.)
            PolE=sin(SrcPZen(Iperm(j))*pi/180.) * sin(SrcPAzi(Iperm(j))*pi/180.)
            Polh=cos(SrcPZen(Iperm(j))*pi/180.)
            write(28,"(I4,F14.6, 3F12.5, F12.1, I5, 2F7.2,  3F6.1, F10.1, 3f7.3)") &
!               i_peak, SourceTime_ms(i_Peak), SourcePos(1:3,i_Peak)/1000.,  ABS(SourceIntensity(i_Peak)), &
               Label(1,Iperm(j)),RA(1:4,j), SrcI20_r(j), &
               SrcWidth(Iperm(j)), SrcChi2(Iperm(j)), SrcISpr(Iperm(j)), &
               SrcUn(Iperm(j)), SrcLin(Iperm(j)), SrcCirc(Iperm(j)), Lin1, PolN, PolE, Polh
         enddo
         close(unit=28)
         Call GLEplotControl(PlotType='SrcsPltPol', PlotName=trim(PlotFile)//'IPol', &
               PlotDataFile=trim(SelFileName), Submit=.false.)
         Call GLEplotControl(PlotType='SrcsPltLoc', PlotName=trim(PlotFile), &
               PlotDataFile=trim(SelFileName), Submit=.false.)
      Else
         OPEN(unit=28,FILE=TRIM(SelFileName)//'.plt', FORM='FORMATTED',STATUS='unknown',ACTION='write') ! space separated values
         write(28,"(6F9.3,2F9.3,1x,A,1x,F7.1,' 0 ')") NEhtBB(:), trim(ZoomBox), t_offset         ! gleRead:  xMin xMax yMin yMax zMin zMax tMin tMax ZmBx$ t_offst
         write(28,"(1x,I2,I8,1x,A, ' 0 ',2x,A,1x,F6.2, 6(1x,g10.4),f7.3,1x,A)")  &
            LongTrackNr, SourcTotNr, Trim(Qual), trim(Image), AmplitudePlot,  &
            a1,b1,c1,a2,b2,c2,d2, ' ! by DataSelect' ! gleRead:  NTracks SourcTotNr Q t_start label$ AmplitudePlot a1 b1 c1 a2 b2 c2 d2
         write(28,"('! SourceLab',T15,'time',T29,'North=y',T46,'East=x',T61,'h=z',T80,'Ampl',T87,'Wl',T92,'Wr')")
         !
         Do j=1,SourcTotNr  !  finally write all accepted source to file
            write(28,"(1x,i8,4(1x,f13.6),1x,F12.2)")  Label(1,Iperm(j)), RA(1:4,j), SrcI20_r(j)
   !         write(28,"(1x,i8,4(1x,f13.6),1x,F12.2,2x,I3,2x,I3)")  Label(1,Iperm(j)),RA(1:4,j), SrcI20_r(j)
   !         write(29,"(i6,',',4(f12.6,','),g13.6,',',3f6.2,',',5g13.3)") i_slice, t_ms-t_shft, &
         enddo
         close(unit=28)
         If(BckgrFile.ne.'') then
            Call GLEplotControl(PlotType='SrcsPltLocBckgr', PlotName=trim(PlotFile), &
                  PlotDataFile=trim(SelFileName), Submit=.false., Bckgr=TRIM(DataFolder)//trim(BckgrFile) )
         Else
            Call GLEplotControl(PlotType='SrcsPltLoc', PlotName=trim(PlotFile), &
                  PlotDataFile=trim(SelFileName), Submit=.false.)
         EndIf
      EndIf
      !
      !
      write(2,*) 'testing trim(ZoomBox):', trim(ZoomBox)
      write(2,*) 'testing LongTrackNr, SourcTotNr, Qual,AmplitudePlot:', LongTrackNr, SourcTotNr, trim(Qual),AmplitudePlot
      write(2,*) SourcTotNr,' sources in plots:',trim(Image),' ; to file: ',TRIM(SelFileName)//'.plt'
      !
!Used as  $gle -d pdf -o IntfMx_dH2t3.pdf ${UtilDir}SourcesPlot.gle ${FlashFolder}/files/H2t3IntfSpecPowMx_d
      !  Corr (=correlate sources following Frankfurt AI paper)
      If((Corr_dD.gt.0) .and. (Corr_Dnr .gt. 0)) Then
            Call ApplyCorrelator(RA, SrcI20_r, SourcTotNr, Aweight, Corr_dD, Corr_Dnr, Corr_dtau)
      EndIf
      !  Qual (=analysis of quality factors, make 2D scatter plots)
      If(QualPlot) Then
            Call ApplyQualityAna
      endif
      !
      ! Clean some temporary files
      j=0
      do i=1,TrackNr
         If(TrackENr(i).lt. LongTrack_Min) cycle
         j=j+1
         if(j.gt.LongTrackNr) exit
         write(extension,"(i1,A4)") j,'.plt' !,&
         Call GLEplotControl(CleanPlotFile=TRIM(Image)//trim(extension))
         If(PolarAna) Call GLEplotControl(CleanPlotFile=TRIM(Image)//'Angls'//trim(extension))
      EndDo
      !
      DeAllocate( SrcI20_r )
998   write(2,*) '======================================='
      DeAllocate( RA, Label, Iperm )
      If(DS_Mode .eq. 1) Then  ! Impulsive imager data
         DeAllocate( QualIndic )
      ElseIf(DS_Mode .eq. 2) Then ! TRI-D Imager
         DeAllocate( Stk_NEh )
      ElseIf(DS_Mode .eq. 3) Then ! ATRID
         DeAllocate( SrcWidth, SrcChi2, SrcI3, SrcUn, SrcLin, SrcISpr, SrcCirc, SrcPZen, SrcPAzi, Stk_NEh )
      EndIf
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
   Use DS_Select, only : SourcTotNr, Label, QualIndic, Iperm
   IMPLICIT none
   Integer :: i, i_src
   !
         Write(2,*) 'option: ','Qual'
   !
   OPEN(unit=29,FILE=TRIM(DataFolder)//trim(Image)//'QuAna.dat',FORM='FORMATTED',STATUS='unknown')  ! will contain all accepted sources
   write(29,"(2x,A,1x,I6,1x,F6.3,i5,' !')") trim(Image), SourcTotNr
   Write(2,*) 'Quality Indicator analysis, writing ',trim(Image), SourcTotNr
   Do i_src=1,SourcTotNr
      write(29,*) Label(1:4,i_src), QualIndic(1:4,Iperm(i_src))  ! 'not sure the two are correctly sorted',
   EndDo ! i_src=1,SourcTotNr-1
   close(unit=29)
   Write(2,*) 'ApplyQualityAna; File:',TRIM(DataFolder)//trim(Image)//'QuAna.dat',' with ',SourcTotNr, ' data lines'
   !
   Call GLEplotControl(PlotType='QualContrPlot', PlotName=TRIM(Image)//'QC', &
         PlotDataFile=TRIM(DataFolder)//trim(Image), Submit=.false.)
End Subroutine ApplyQualityAna
!--------------------------------------------------
Subroutine PolarizationAna()  ! du -sh  directory size
   Use constants, only : dp
   !use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, DataFolder, FlashName, Windows, RunMode
   !Use DS_Select, only : Image !,  WriDir
   !use GLEplots, only : GLEplotControl
   Use DS_Select, only : Stk_NEh,SourcTotNr, Label, Iperm
   IMPLICIT none
   Integer :: i, i_src
   Complex :: AveStks(1:3,1:3)
   Real :: PolZen(1:3), PolAzi(1:3), PolMag(1:3), PoldOm(1:3)
   logical :: prin=.true.
   !
   !
   !OPEN(unit=29,FILE=TRIM(DataFolder)//trim(Image)//'QuAna.dat',FORM='FORMATTED',STATUS='unknown')  ! will contain all accepted sources
  ! write(29,"(2x,A,1x,I6,1x,F6.3,i5,' !')") trim(Image), SourcTotNr
   Write(2,*) 'Polarization analysis,  ', SourcTotNr
   AveStks(:,:)=0
   Do i_src=1,SourcTotNr ! Add all stokes matrices
      AveStks(1,1:3)=AveStks(1,1:3)+ Stk_NEh(1:3,Iperm(i_src))
      AveStks(2,2:3)=AveStks(2,2:3)+ Stk_NEh(4:5,Iperm(i_src))
      AveStks(3,3)=AveStks(3,3)+ Stk_NEh(6,Iperm(i_src))
   EndDo ! i_src=1,SourcTotNr
   AveStks(2:3,1)=AveStks(1,2:3)
   AveStks(3,2)=AveStks(2,3)
   !
   Call PolPCACath(AveStks, PolZen, PolAzi, PolMag, PoldOm, prin)  ! in EIFitter
   !Write(2,*) 'PolarizationAna:'
   !
   !Call GLEplotControl(PlotType='QualContrPlot', PlotName=TRIM(Image)//'QC', &
   !      PlotDataFile=TRIM(DataFolder)//trim(Image), Submit=.false.)
End Subroutine PolarizationAna
!--------------------------------------------------
!-----------------------------------------------------
! ----------------------------------------------------------------------------------------
!     shellin = 'gle /d pdf '//trim(dummy3(i))//'.gle'
!     CALL system(shellin)
!     shellin = 'epstopdf '//trim(dummy3(i))//'.eps'
!     call system(shellin)
