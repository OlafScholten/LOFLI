Module DS_Select
!   Use Tracks, only : PreDefTrack(0:3,tNr,1:PreDefTrackNr), tNrS_max(i_Pre), PreDefTrackNr
   Use constants, only : dp, CI, pi, c_l
   use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, DataFolder, FlashName, Windows, RunMode
   !
   !Character*130 :: datfile
   Character*100, save :: WriDir
   !
   integer, parameter :: maxd=123000
   real*8 :: RA(4,maxd)  ! 1=t [ms]  2-4= E,N,h in [km]
   Integer :: Label(4,maxd)  ! Unique source label and amplitude and widths
   real :: QualIndic(4,maxd)  ! Quality indicators
   Integer :: SourcTotNr  ! number of sources stored in RA, passing the selection criteria
   !
   Integer, parameter :: Nmax_Ampl=5000
   Real(dp) :: MaxAmplFitPercent=100.
   Integer :: FitFunc=0  ! Fitfunftion=0 : exp ; =1 powerlaw
   Real(dp) :: Ampl_Hist(Nmax_Ampl)  ! Amplitude histogram filled out in "AmplitudeFitTracks"
   Real(dp) :: Ampl_Hist_W(Nmax_Ampl)  ! Amplitude histogram filled out in "AmplitudeFitTracks"
   Real(dp) :: Ampl(Nmax_Ampl)  ! Amplitude histogram filled out in "AmplitudeFitTracks"
   Integer :: i_AmplTh  ! Calculated in "AmplitudeFitTracks" to the index of the first amplitude for which the probability spectrum is fitted
   Real(dp) :: Nrm, a1,b1,c1,d1,ChiSq1,Nrm1, a2,b2,c2,d2,ChiSq2,Nrm2, a3,b3,c3,d3,ChiSq3 ! parameters for describing the amplitude distribution
   Real*8 :: AmplScale=1.d3 ! Scaling factor for amplitude to convert from integer to unity-normalized real
   !Real*8 :: d_AmplScale
   !
   Real(dp), save :: xMin,xMax,yMin,yMax,zMin,zMax,tMin,tMax, SMPowCut ! RatMax=ratio (Max rim)/peak
   Real(dp), save :: xyztBB(8), NEhtBB(8) ! RatMax=ratio (Max rim)/peak
   Real(dp), save :: Tiny=epsilon(Tiny), Small
   Integer, parameter :: NDataFil=50
   character*20 :: datafile(NDataFil)
   character*100 :: PlotName, ZoomBox, Image
   Character(len=180) :: BckgrFile
   Real(dp), save :: TimeBase, t_offset, AmplitudePlot
   Real(dp), save :: tCutl, tCutu, CutSigmaH, LinCutH, RMS_ns
   Real(dp), save ::  Corr_dD, Corr_dtau
   Integer :: Corr_Dnr
   Integer, save :: FileStatN, DelNEff
   Integer, save :: DS_Mode
   Logical, save :: QualPlot
   Logical, save :: StatNCut
contains
!===========================================================
Subroutine DS_ReadCntrl(Fini)
   Use constants, only : dp, CI, pi, c_l
   Use TrackConstruct, only : PreDefTrackNr, PreDefTrackFile, SrcDensTimeResol
   Use TrackConstruct, only : TrackNr, LongTrackNr, LongTrack_Min, TrackNrMax, NLongTracksMax, TrackLenMax
   Use TrackConstruct, only : Wtr, Aweight, MaxTrackDist, TimeWin, HeightFact, dt_MTL, TrackENr, TrackE, TrackNrLim
   use LOFLI_Input
   IMPLICIT none
   Integer, intent(out) :: Fini
   Integer :: j, nxx
   Character(len=20) :: RunOption !, shellin
   Character(len=120) :: PlotFile
   NAMELIST /Parameters/ datafile, BckgrFile, PlotName, TimeBase, &
      SMPowCut, DelNEff, MaxAmplFitPercent, FileStatN, StatNCut, SrcDensTimeResol, &
      tCutl, tCutu, RunOption , xyztBB, NEhtBB, CutSigmaH, LinCutH, RMS_ns, DelNEff, ZoomBox &
   , MaxTrackDist, Wtr, Aweight, TimeWin, PreDefTrackFile, dt_MTL, HeightFact, NLongTracksMax, AmplitudePlot &
   , Corr_dD, Corr_Dnr, Corr_dtau, QualPlot
   !
   OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='DataSelect.out')
   !
      !  Read additional options
      !  Corr (=correlate sources following Frankfurt AI paper)
      !  Qual (=analysis of quality factors, make 2D scatter plots)
      !  Bckg (=make plot including data from indicated file in grey, in folder 'files')
   AmplitudePlot=Tiny ! set it to some strange value, just to determine if it has been set on input
   xMin=Tiny; xMax=Tiny; yMin=Tiny; yMax=Tiny; zMin=Tiny; zMax=Tiny; tMin=Tiny; tMax=Tiny
   xyztBB(:)=Tiny
   NEhtBB(:)=Tiny
   TimeBase=tiny  ! to initialize it to a value that will not be put on input, to be able to check if set at all
   tCutl=0.; tCutu=0.
   datafile(:)=''
   RunOption=''
   MaxAmplFitPercent=0.1
   SrcDensTimeResol=0.1 ! [ms]
   BckgrFile=''
   datafile(:)=''
   PlotName=''
   SMPowCut=0.
   CutSigmaH=-1.  !  No cut of sigma(h)
   LinCutH = 0.
   RMS_ns = 10.
   StatNCut=.false.
   MaxTrackDist= 0.1 ! [km]
   Wtr= 0.5      ! Weight for next point on track
   dt_MTL=0.5
   TimeWin=0.1
   HeightFact=3.
   NLongTracksMax=0
   Aweight=tiny
   !
   Fini=0
   read(*,NML = Parameters,IOSTAT=nxx)
   If(nxx.ne.0) Then
      Backspace(Unit=5)  !   (unit = 5, file = 'stdin')
      Write(2,NML = Parameters)
      Read(*,"(A100)") BckgrFile
      Write(2,*) TRIM(BckgrFile)
      Write(2,*) '======================  End input or namelist "Parameters" error'
      !write(2,NML = Parameters,IOSTAT=nxx)
      Fini=1
      return
   EndIf
   j = iachar(RunOption(1:1))  ! convert to upper case if not already
   if (j>= iachar("a") .and. j<=iachar("z") ) then
      RunOption(1:1) = achar(j-32)
   end if
   If(SrcDensTimeResol.eq.tiny) SrcDensTimeResol=TimeWin
   write(2,*) 'RunOption:', RunOption
      SELECT CASE (RunOption(1:1))
      CASE("I")  ! Impulsive Imager
         Write(2,*) 'Image data from Impulsive Imager, PlotName=',TRIM(PlotName)
         DS_Mode=1
         If(AmplitudePlot.eq.tiny) AmplitudePlot=10.
         write(2,*) 'Source selection criteria:'
         Call PrintValues(datafile,'datafile', 'Raw Imager data-file label (Impulsive).')  ! width of plot?
         If(CutSigmaH.gt. 0.) Then
            Call PrintValues(CutSigmaH,'CutSigmaH', 'Condition: [Sigma(h) < CutSigmaH], unit [m].')  ! width of plot?
            Call PrintValues(LinCutH,'LinCutH', 'Condition: [Sigma(h) * h/LinCutH < CutSigmaH] when [h < LinCutH], unit [km].')  ! width of plot?
         Else
            write(2,*) 'Cut of sigma(h) requires positive value for "CutSigmaH"'
         EndIf
         Call PrintValues(RMS_ns,'RMS_ns', 'Condition: [sqrt(Chi^2) < RMS_ns], unit [ns].')  ! width of plot?
         Call PrintValues(DelNEff,'DelNEff', 'Condition: [(Max_EffAnt - N_EffAnt) <= DelNEff], hard when [DelNEff > 0]'//&
            'otherwise soft: [(sqrt(Chi^2)-RMS_ns*(Max_EffAnt-N_EffAnt)/DelNEff) < RMS_ns] ')
         Call PrintValues(tCutl,'tCutl', 'Lower end of block filter in time.')  ! width of plot?
         Call PrintValues(tCutu,'tCutu', 'Upper end of block filter in time.')  ! width of plot?
         write(2,*) 'Image control:'
         Call PrintValues(PlotName,'PlotName', 'Added to flash name to label plots')  ! width of plot?
         Call PrintValues(BckgrFile,'BckgrFile', 'Name of file that is displayes as grey background plot.')  ! width of plot?
         Call PrintValues(ZoomBox,'ZoomBox', 'Name of zoom-in plot.')  ! width of plot?
         Call PrintValues(AmplitudePlot,'AmplitudePlot', 'Base-size for source-dots in plots.')  ! width of plot?
         Call PrintValues(NLongTracksMax,'NLongTracksMax', 'Maximum number of long tracks to be plotted.')  ! width of plot?
         If(NLongTracksMax.gt.0) Then
            Call PrintValues(MaxTrackDist,'MaxTrackDist', 'Max. distance between sources to include on a track.')  ! width of plot?
            Call PrintValues(Wtr,'Wtr', 'Weight of newest source for track centroid.')  ! width of plot?
            Call PrintValues(Aweight,'Aweight', 'Importance of intensity in weight of newest source for track centroid.')  ! width of plot?
            Call PrintValues(TimeWin,'TimeWin', '[ms]. Width of gaussian in time to weigh the sources for mean track position.')  ! width of plot?
            Call PrintValues(PreDefTrackFile,'PreDefTrackFile', 'Label of file that contains a track to be included.')  ! width of plot?
            Call PrintValues(dt_MTL,'dt_MTL', '[ms]. Time-step for constructing tracks.')  ! width of plot?
            Call PrintValues(HeightFact,'HeightFact', 'Relative height scale for calculating distances.')  ! width of plot?
         EndIf!(NLongTracksMax.gt.0)
         If((Corr_dD.gt.0) .and. (Corr_Dnr .gt. 0)) Then
            Call PrintValues(Corr_dD,'Corr_dD', 'Distance step size for time-distance correlator.')  ! width of plot?
            Call PrintValues(Corr_Dnr,'Corr_Dnr', 'max. number of distance bins for time-distance correlator.')  ! width of plot?
            Call PrintValues(Corr_dtau,'Corr_dtau', 'Time step size for time-distance correlator.')  ! width of plot?
         Else
            Write(2,*) 'TD correlators require positive values for "Corr_dD" and "Corr_Dnr"'
         EndIf
         Call PrintValues(QualPlot,'QualPlot', 'make scattter plot of different quality indicators.')  ! width of plot?
      CASE("T")  ! TRI-D Imager
         Write(2,*) 'Image data from TRI-D Imager, PlotName extension=',TRIM(PlotName)
         DS_Mode=2
         write(2,*) 'Source selection criteria:'
         Call PrintValues(datafile,'datafile', 'Raw Imager data-file label (TRI_D).')  ! width of plot?
         Call PrintValues(TimeBase,'TimeBase', 'Time offset.')  ! width of plot?
         Call PrintValues(SMPowCut,'SMPowCut', 'Threshol Intensity.')  ! width of plot?
         !Call PrintValues(RatMaxCut,'RatMaxCut', 'ratio (Max rim)/peak.')  ! width of plot?
         Call PrintValues(StatNCut,'StatNCut', 'Keep only the "FileStatN" strongest sources per TRI-D image.')  ! width of plot?
         Call PrintValues(tCutl,'tCutl', 'Lower end of block filter in time.')  ! width of plot?
         Call PrintValues(tCutu,'tCutu', 'Upper end of block filter in time.')  ! width of plot?
         write(2,*) 'Image control:'
         Call PrintValues(PlotName,'PlotName', 'Added to data name to label plots')  ! width of plot?
         Call PrintValues(BckgrFile,'BckgrFile', 'Name of file that is displayes as grey background plot.')  ! width of plot?
         Call PrintValues(ZoomBox,'ZoomBox', 'Name of zoom-in plot.')  ! width of plot?
         Call PrintValues(AmplitudePlot,'AmplitudePlot', 'Base-size for source-dots in plots.')  ! width of plot?
         Call PrintValues(MaxAmplFitPercent,'MaxAmplFitPercent', 'Max amplitude fitted with modifies power law.')  ! width of plot?
         Call PrintValues(FileStatN,'FileStatN', 'Print the strength of the N^th strongest source.')  ! width of plot?
         If(NLongTracksMax.gt.0) Then
            Call PrintValues(MaxTrackDist,'MaxTrackDist', 'Max. distance between sources to include on a track.')  ! width of plot?
            Call PrintValues(Wtr,'Wtr', 'Weight of newest source for track centroid.')  ! width of plot?
            Call PrintValues(Aweight,'Aweight', 'Importance of amplitude in weight of newest source for track centroid.')  ! width of plot?
            Call PrintValues(TimeWin,'TimeWin', '[ms]. Width of gaussian in time to weigh the sources for mean track position.')  ! width of plot?
            Call PrintValues(PreDefTrackFile,'PreDefTrackFile', 'Label of file that contains a track to be included.')  ! width of plot?
            Call PrintValues(dt_MTL,'dt_MTL', '[ms]. Time-step for constructing tracks.')  ! width of plot?
            Call PrintValues(HeightFact,'HeightFact', 'Relative height scale for calculating distances.')  ! width of plot?
            Call PrintValues(SrcDensTimeResol,'SrcDensTimeResol', '[ms]. Source Density Time Resolution.')  ! width of plot?
         EndIf!(NLongTracksMax.gt.0)
         If((Corr_dD.gt.0) .and. (Corr_Dnr .gt. 0)) Then
            Call PrintValues(Corr_dD,'Corr_dD', '[km]. Distance step size for time-distance correlator.')  ! width of plot?
            Call PrintValues(Corr_Dnr,'Corr_Dnr', 'max. number of distance bins for time-distance correlator.')  ! width of plot?
            Call PrintValues(Corr_dtau,'Corr_dtau', '[ms]. Time step size for time-distance correlator.')  ! width of plot?
         Else
            Write(2,*) 'TD correlators require positive values for "Corr_dD" and "Corr_Dnr"'
         EndIf
      CASE DEFAULT  ! Help
         Write(*,*) 'specified RunOption: "',Trim(RunOption),'", however the possibilities are:'
         Write(*,*) '- "Impulsive" for the impulsive Imager'
         Write(*,*) '- "TRI-D" for the TRI-D imager with polarization observables, accounting for antenna function'
         Stop 'No valid RunOption specified in namelist input'
   End SELECT
   If(NEhtBB(1).ne.Tiny) Then
      xyztBB(1:2)=NEhtBB(3:4)
      xyztBB(3:4)=NEhtBB(1:2)
      xyztBB(5:8)=NEhtBB(5:8)
      write(2,*) 'equivalent xyztBB=', xyztBB(:)
   Else
      NEhtBB(3:4)=xyztBB(1:2)
      NEhtBB(1:2)=xyztBB(3:4)
      NEhtBB(5:8)=xyztBB(5:8)
   EndIf
   write(2,"(A,3(2F9.3,2x),A,2F10.4,A)") 'used NEhtBB=', xyztBB(1:6),'[km]', xyztBB(7:8),' [ms]'
   !
   If(ZoomBox.ne.'NoBox')  ZoomBox=TRIM(DataFolder)//ZoomBox
   !Write(2,"(A,F7.2,A,F7.2,A,F7.2)") 'quality cuts: sigma(z)=', ICVal,'[m] @', Height_1,'[km]'
   !Write(2,"(A,I2,A)") 'Effective antenna cut=',DelNEff,' less than the maximum'
   Do j=1,4
      If(xyztBB(2*j-1).gt.xyztBB(2*j)) Then
         write(2,*) 'Min & Max given in wrong order:',xyztBB(2*j-1),xyztBB(2*j)
         stop "wrong limits"
      EndIf
   EndDo
   !
   Return
End Subroutine DS_ReadCntrl
!===========================================================
End Module DS_Select
!===========================================================
!===============================================
Subroutine SelectBB(xMinR, xMaxR, yMinR, yMaxR, zMinR, zMaxR, tMinR, tMaxR)
!     Possibly overwrite the BB (xyztRBB) read from file
   Use constants, only : dp
   Use DS_Select, only : xyztBB, tiny
   implicit none
   Real(dp), intent(in) :: xMinR, xMaxR, yMinR, yMaxR, zMinR, zMaxR, tMinR, tMaxR
   !write(2,*) tiny,'xyztBB before:', xyztBB(:)
   !write(2,*) 'BB:', xMinR, xMaxR, yMinR, yMaxR, zMinR, zMaxR, tMinR, tMaxR
   If(xyztBB(1).eq.Tiny) xyztBB(1)=xMinR
   If(xyztBB(2).eq.Tiny) xyztBB(2)=xMaxR
   If(xyztBB(3).eq.Tiny) xyztBB(3)=yMinR
   If(xyztBB(4).eq.Tiny) xyztBB(4)=yMaxR
   If(xyztBB(5).eq.Tiny) xyztBB(5)=zMinR
   If(xyztBB(6).eq.Tiny) xyztBB(6)=zMaxR
   If(xyztBB(7).eq.Tiny) xyztBB(7)=tMinR
   If(xyztBB(8).eq.Tiny) Then
      xyztBB(8)=tMaxR
      write(2,*) 'new xyztBB=', xyztBB(:)
   EndIf
End Subroutine SelectBB
!===========================================================
Subroutine DS_ReadSelData
!  Read Selected Data
   ! Apply quality criteria to select sources
   ! On return:   Selected sources in RA
   Use DS_Select, only : DS_Mode, datafile, PlotName, Image
   Use constants, only : dp
   use DataConstants, only : FlashFolder, FlashName
   IMPLICIT none
   If(DS_Mode .eq. 1) Then  ! Impulsive imager data
      Call DS_ReadSelData_Imp
      Image='I-'//TRIM(FlashName)//TRIM(PlotName)
   Else  ! TRI-D Imager
      Call DS_ReadSelData_TRID
      Image='T-'//TRIM(datafile(1))//TRIM(PlotName)
   EndIf
   !write(2,*) 'basis Images naming:',TRIM(Image)
   Return
   !
End Subroutine DS_ReadSelData
!===========================================================
Subroutine DS_ReadSelData_Imp
!  Read Selected Data
   Use constants, only : dp
   use DataConstants, only : DataFolder !, ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows, RunMode
   Use DS_Select, only : tiny, ZoomBox, tCutl, tCutu, CutSigmaH, LinCutH, RMS_ns,  DelNEff
   Use DS_Select, only : QualIndic, RA, Maxd, Label, SourcTotNr, xyztBB, t_offset, datafile
   Use DS_Select, only : Nrm, a1,b1,c1,d1,ChiSq1,Nrm1, a2,b2,c2,d2,ChiSq2,Nrm2, a3,b3,c3,d3,ChiSq3 ! parameters for describing the amplitude distribution
   IMPLICIT none
   !Real(dp) ::
   !Integer ::
   Character(len=180) :: lineTXT
   Character(len=100) :: txt
   Integer :: i,j, Nxx
   real*8 :: x,y,z,t, val, sigma(1:3), Q, Chi, t_start
   integer :: N_EffAnt, Max_EffAnt, Ampl, Wl, Wu
   !
   SourcTotNr=0
      !If(.not. Windows) call system_mem_usage(valueRSS)
   OPEN(unit=28,FILE=trim(datafile(1))//'.csv', FORM='FORMATTED',STATUS='OLD',ACTION='READ', IOSTAT=nxx) ! space separated values
   If(nxx.ne.0) then   ! check weather this particular file can be opened
      Write(2,*) 'Problems opening file:"',trim(datafile(1))//'.csv','"'
      Close(Unit=28)
      return
   endif
   !
   Do
      Read(28,"(A180)",IOSTAT=nxx) lineTXT
      If(nxx.ne.0) stop 'DS_ReadSelData_Imp-readingC'
      If(lineTXT(1:1).ne.'!') exit
      !Write(2,"(A)") TRIM(lineTXT)
   Enddo
   !write(2,*) 'from data file, lineTXT:', lineTXT
   read(lineTXT,*,end=998) t_start,txt     ! from the first line that does not start with !
   write(2,*) 'time offset=',t_start,'[s]'
   !
   ! Start actual reading
   j=0
   do    ! loop over all source points in this file
      NXX=0
      read(28,"(A180)",end=998)  lineTXT
      read(lineTXT,*,iostat=nxx)  i,y,x,z,t,val,sigma, N_EffAnt, Max_EffAnt, Ampl, Wl, Wu
      If(nxx.ne.0) then
         !write(2,*) lineTXT
         !flush(unit=2)
         read(lineTXT,*,iostat=nxx)  i,y,x,z,t,val,sigma, N_EffAnt, Max_EffAnt
         !write(2,*) Nxx,lineTXT
         If(nxx.ne.0) Then
            read(28,*,end=998,iostat=nxx)  i,y,x,z,t,val,sigma
         EndIf
         Ampl=100+mod(i,10) ; Wl=1 ; Wu=1
         !TOrder=+1  ! older version
         !write(2,*) 'Read error:',i,y,x,z,t,val,sigma, N_EffAnt, Max_EffAnt, Ampl, Wl, Wu
      Endif
      x=x/1000.d0                   ! convert to [km]
      y=y/1000.d0                   ! convert to [km]
      z=z/1000.d0                   ! convert to [km]
      t = (t-t_start)*1000.  ! convert to [ms]
      z=abs(z)
      if(x.le.xyztBB(1) .or. x.ge.xyztBB(2)) cycle
      if(y.le.xyztBB(3) .or. y.ge.xyztBB(4)) cycle   ! [km]
      if(z.le.xyztBB(5) .or. z.ge.xyztBB(6)) cycle
      if(t.le.xyztBB(7) .or. t.ge.xyztBB(8)) cycle   ! [ms]
      if(t.gt.tCutl .and. t.lt.tCutu) cycle  ![ms]
      !
!     CutSigma_h  <-- ICVal
!     LinCutH <-- Height_1
!     RMS_ns <-- Chi_max
         If(CutSigmaH.gt.0.) Then
            If(z.gt. LinCutH) then
               Q=sigma(3)
            Else
               Q=sigma(3)*z/LinCutH
            Endif
            If(Q.gt. CutSigmaH .and. val.gt.CutSigmaH) cycle
         EndIf
         Q=sqrt(val)
         If(DelNEff .lt. 0) then
            Q=Q - RMS_ns*(Max_EffAnt-N_EffAnt)/DelNEff
         else
            If( (Max_EffAnt - N_EffAnt) .gt. DelNEff) cycle
         Endif
         If(Q.gt. RMS_ns) cycle
         !
         j=j+1
         RA(1,j)=t ;  RA(2,j)=x ;  RA(3,j)=y ;  RA(4,j)=z   ! units [ms], [km], [km], [km]
         Label(1,j)=i ;Label(2,j)=Ampl ; Label(3,j)=Wl ; Label(4,j)=Wu
         QualIndic(1,j)=Chi ; QualIndic(2:4,j)=sigma(1:3)
         if(j.eq.maxd) then
            write(*,*) 'Max. dimension reached of', maxd,' at',i
            write(2,*) 'Max. dimension reached of', maxd,' at',i
            exit
         EndIf
   enddo     ! loop over all source points in this file
998 continue
   Close(Unit=28)  ! following is mainly important for multiple files-reads for TRI-D
   !For Impulsive:
    write(2,*) 'Last event read:',i,' @ t=',t,'[s]'
    !Qual=RMS_ns
    SourcTotNr=j
    t_offset=t_start*1000.
    Return
End Subroutine DS_ReadSelData_Imp
!===========================================================
Subroutine DS_ReadSelData_TRID
!  Read Selected Data
   Use constants, only : dp, CI, pi, c_l
   use DataConstants, only : DataFolder !, ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows, RunMode
   Use DS_Select, only : tiny
   Use DS_Select, only : NDataFil, datafile,  TimeBase, t_offset, AmplitudePlot, SMPowCut
   Use DS_Select, only : tCutl, tCutu, SourcTotNr, FileStatN, StatNCut
   Use DS_Select, only : QualIndic, RA, Maxd, Label, AmplScale, xyztBB
   Use DS_Select, only : Nrm, a1,b1,c1,d1,ChiSq1,Nrm1, a2,b2,c2,d2,ChiSq2,Nrm2, a3,b3,c3,d3,ChiSq3 ! parameters for describing the amplitude distribution
   Use unque, only : Selection_sort
   IMPLICIT none
   Character*20 :: BoxOption !, Utility, release
   character*100 ::  FileMain, extension, datfile !, OutFileLabel, SelFileName, PlotName, shellin
   Character(len=180) :: lineTXT !, OS
   integer :: i,j, nxx, i_datfil, PrevSourcTotNr
   Real*8 :: t_start
   Real(dp) :: t,x,y,z,SMPow
   Real(dp) :: Q, AmpltPlotRead
   Real(dp) :: xMinR,xMaxR,yMinR,yMaxR,zMinR,zMaxR,tMinR,tMaxR
   integer :: i_sequence, LastChar  ! i_eo,i_BM,i_datfil, valueRSS, i_sequence, LastChar
!   Logical :: ZoomClip=.false.
   Logical :: First28=.true.
   Integer :: SourceNrIncr, ScalAmpl(1:4000), AmplMax, AmplStatN, FileStatNp, AmplStatNp, FileStatNm, AmplStatNm
   !
   extension='IntfSpecPowMx_d.dat' ! For TRI-D (beamformed) data file
   First28=.true.  ! Read timebase from sources file
   SourcTotNr=0
   PrevSourcTotNr=0
            !If(.not. Windows) call system_mem_usage(valueRSS)
   Do i_datfil=1,NDataFil  ! select one particular file (may be sequence)
      If((i_datfil.gt.1).and. (datafile(i_datfil).eq."")) exit
      nxx=0
      i_sequence=0
      LastChar=LEN_TRIM(datafile(i_datfil))
      Do While (i_sequence.lt.100)  ! run over the full sequnce of files for this name
         !write(2,*) trim(datafile(i_datfil)), LastChar, datafile(i_datfil)(LastChar:LastChar)
         datfile=datafile(i_datfil)
         If(LastChar.ge.2) Then
            If(datafile(i_datfil)(LastChar-1:LastChar-1).eq.'#') then
               If(i_sequence.eq.0) Then
                  datfile=datafile(i_datfil)(1:LastChar-2)
               Else
                  write(datfile,"(A,i2.2)") datafile(i_datfil)(1:LastChar), i_sequence
               EndIf
               i_sequence=i_sequence+1 ! exit when =0
               !write(2,*) datfile
            EndIf
         EndIf
         FileMain=datfile  ! Base for output file
         datfile=TRIM(DataFolder)//trim(datfile)//TRIM(extension)
         OPEN(unit=28,FILE=trim(datfile), FORM='FORMATTED',STATUS='OLD',ACTION='READ', IOSTAT=nxx) ! space separated values
         If(nxx.ne.0) then   ! check weather this particular file can be opened
            Write(2,*) 'Probelems opening file:"',trim(datfile),'"'
            Close(Unit=28)
            exit
         endif
         !
         Do
            Read(28,"(A180)",IOSTAT=nxx) lineTXT
            If(nxx.ne.0) stop 'InterfSrcSel-readingC'
            If(lineTXT(1:1).ne.'!') exit
            !Write(2,"(A)") TRIM(lineTXT)
         Enddo
         !write(2,*) 'from data file, lineTXT:', lineTXT
         read(lineTXT,*,IOSTAT=nxx) xMinR, xMaxR, yMinR, yMaxR, zMinR, zMaxR, tMinR, tMaxR, BoxOption, t_start     ! from the first line that does not start with !
         Call SelectBB(xMinR, xMaxR, yMinR, yMaxR, zMinR, zMaxR, tMinR, tMaxR)
         If(First28) Then
            If(TimeBase.eq.tiny) then
               t_offset=t_start
            Else
               t_offset=TimeBase
            EndIf
            First28=.false.
         EndIf
         Read(28,*) i, j, Q, SMPow, lineTXT, AmpltPlotRead, a1, b1, c1, a2, b2, c2, d2
         If(AmplitudePlot.eq.Tiny) AmplitudePlot=AmpltPlotRead
         !Write(2,*)  Q, SMPow, TRIM(FileLabel), AmplitudePlot, a1, b1, c1, a2, b2, c2, d2
         !write(2,*) TRIM(BoxOption)
         If(TRIM(BoxOption).eq.'IntfBox') then
            !write(2,*) TRIM(BoxOption)
            Do i=1,8 ! skip 8 records
               Read(28,"(A180)",IOSTAT=nxx) lineTXT
            Enddo
         EndIf
         !
         t_start = t_start-t_offset
         !
         ! Start actual reading
         j=0
         do    ! loop over all source points in this file
            NXX=0
!         write(29,"(i6,',',4(f11.5,','),g13.6,',',3f6.2,',',5g13.3)") i_slice, t_ms-t_shft, &
!         PixLoc(2)/1000., PixLoc(1)/1000., PixLoc(3)/1000., SMPowMx, &   ! MaxPowPix(:,i)
!                                                             MaxSmPowQ(i_slice), MaxSmPowU(i_slice), MaxSmPowV(i_slice), MaxSmPowI3(i_slice) ,s,X
!         s=sqrt(MaxSmPowQ(i_slice)*MaxSmPowQ(i_slice)+MaxSmPowU(i_slice)*MaxSmPowU(i_slice))
!         X=(ATAN2(MaxSmPowU(i_slice),MaxSmPowQ(i_slice))/2.+AngOff)*180./pi
!         MaxSmPowI3  = I3/SMPowMx  where  SMPowMx depends on  PixPowOpt
            read(28,*,iostat=nxx)  i,t,x,y,z,SMPow  ! already in proper units for plotting
            t=t+t_start  ! Same t-scale now as for first file
            If(nxx.gt.0) then
               write(2,*) 'Read error:',i,t,x,y,z,SMPow
            ElseIf (nxx.lt.0) Then
               !Write(2,*) 'EOF reached'
               exit
            Endif
            if(x.le.xyztBB(1) .or. x.ge.xyztBB(2)) cycle
            if(y.le.xyztBB(3) .or. y.ge.xyztBB(4)) cycle   ! [km]
            if(z.le.xyztBB(5) .or. z.ge.xyztBB(6)) cycle
            if(t.le.xyztBB(7) .or. t.ge.xyztBB(8)) cycle   ! [ms]
            if(t.gt.tCutl .and. t.lt.tCutu) cycle  ![ms]
            !
            !write(2,*)  'DS_ReadSelData_TRID', i,t,x,y,z,SMPow, SMPowCut  ! already in proper units for plotting
            If( SMPow .lt. SMPowCut) cycle
            SourcTotNr=SourcTotNr+1
            RA(1,SourcTotNr)=t ;  RA(2,SourcTotNr)=x ;  RA(3,SourcTotNr)=y ;  RA(4,SourcTotNr)=z   ! units [ms], [km], [km], [km]
            Label(1,SourcTotNr)=i
            Label(2,SourcTotNr)=SMPow*AmplScale
            if(SourcTotNr.eq.maxd) then
               write(*,*) 'Max. dimension reached of', maxd,' at',i
               write(2,*) 'Max. dimension reached of', maxd,' at',i,t
               exit
            EndIf
         enddo     ! loop over all source points in this file
         Close(Unit=28)  ! following is mainly important for multiple files-reads for TRI-D
         !For Impulsive:
         !For TRI-D:
         SourceNrIncr=SourcTotNr-PrevSourcTotNr
         AmplMax=-1
         AmplStatNm=-AmplScale
         AmplStatNp=-AmplScale
         AmplStatN =-AmplScale
         If((FileStatN.gt.1) .and. (SourceNrIncr.ge.1)) then ! determine max amplitude and cuts
            SourceNrIncr=SourcTotNr-PrevSourcTotNr
            ScalAmpl(1:SourceNrIncr)=Label(2,PrevSourcTotNr+1:SourcTotNr)
            Call Selection_sort(ScalAmpl(1:SourceNrIncr))
            AmplMax=ScalAmpl(SourceNrIncr)
            If(FileStatNm.lt.SourceNrIncr) AmplStatNm=ScalAmpl(SourceNrIncr-FileStatNm)
            If(FileStatNp.lt.SourceNrIncr) AmplStatNp=ScalAmpl(SourceNrIncr-FileStatNp)
            If(FileStatN.lt.SourceNrIncr)  AmplStatN =ScalAmpl(SourceNrIncr-FileStatN)
            If(StatNCut) Then
               j=PrevSourcTotNr
               Do i=PrevSourcTotNr+1,SourcTotNr
                  If(Label(2,i) .gt. AmplStatN) Then
                  j=j+1
                  RA(1:4,j)=RA(1:4,i)
                  Label(1:2,j)=Label(1:2,i)
                  EndIf
               EndDo
               SourcTotNr=j
            EndIf
         EndIf
         Write(2,"(A, A15,A,I6,A,I6, 4(1pg11.3))") 'After file: ',trim(FileMain), ', SourcTotNr=',SourcTotNr, &
               ', increment=',SourceNrIncr &
               , AmplStatNp/AmplScale, AmplStatN/AmplScale, AmplStatNm/AmplScale, AmplMax/AmplScale
         PrevSourcTotNr=SourcTotNr
         If(i_sequence.eq.0) exit ! to exit sequence loop
      EndDo ! i_sequence loop
   Enddo
   write(2,*) 'Last event read:',i,' @ t=',t,'[s]'
   !
   !
    Return
End Subroutine DS_ReadSelData_TRID
!=========================================
Module AmpFit
contains
!===============================================
Subroutine AmplitudeFitTracks(TrackNr,TrackLenMax,TrackNrMax,TrackENR,TrackE)
   ! Make Source amplitude histogram and fit this with exp or power law
   !      write(29,"(1x,i8,4(2x,g14.8),3x,F7.2,2x,I3,2x,I3)")  Label(1,j),RA(1:4,j), Label(2,j)/100., Label(3:4,j)
   IMPLICIT none
   Integer, intent(in) :: TrackNr,TrackLenMax,TrackNrMax
   Integer, intent(in) :: TrackENr(TrackLenMax), TrackE(TrackLenMax,TrackNrMax)
   Integer :: i ! ,k, i_Ampl, Count
   !Character*8 :: extension
   !   write(2,*) 'long track',i_LongTr,' has ',t_MTL,'bins of ',dt_MTL,'[ms]'
   do i=1,TrackNr
      !write(extension,"(A2,i1,A4)") '_s',nxx,'.dat' !,&
      !write(extension,"(i1,A4)") i,'.dat' !,&
      !OPEN(unit=27,FILE='files/'//trim(Image)//'_??_'//trim(extension),FORM='FORMATTED',STATUS='unknown')
      !write(27,"('! nr', 4(1x,A14),3x,3(1x,A12),2x,A10)") 't [ms]','E [km]','N','h','v_E','v_N','v_h','v [10^6m/s]'
      !
      write(2,"(A,I2)") 'fitting source strengths along track#',TrackNr
      Call AmplitudeFit(TrackENr(i), SourceNrs=TrackE(1,i))
      !write(2,*) 'AmplitudeFitTracks: AmplitudeFit(',TrackENr(i), TrackE(1,i),i
      !Call Flush(2)
      !close(unit=27)
   enddo
    !
   Return
End Subroutine AmplitudeFitTracks
!===============================================
Subroutine AmplitudeFit(NrSources, SourceNrs, SelFileName)
   ! Make Source amplitude histogram and fit this with exp or power law
   ! d_Ampl is stepsize in amplitude, integer
   !      write(29,"(1x,i8,4(2x,g14.8),3x,F7.2,2x,I3,2x,I3)")  Label(1,j),RA(1:4,j), Label(2,j)/100., Label(3:4,j)
   ! Writes files  trim(datfile)//??//'.dat'
   ! When SourceNrs(i) present it contains the list of source numbers that are to be considered
   ! NrSources: number to be considered
   !Use Tracks, only :  Label
   ! contains amplitudes as integers
   !Use Tracks, only : Nmax_Ampl, AmplScale, MaxAmplFitPercent
   ! Nmax_Ampl: Maximum number of bins in histogram, determines max amplitude that can be fitted.
   ! AmplScale: factor (>1) used to convert amplitudes to integers as needed for Label(2,:)
   ! d_AmplScale: bin-size for histogramming the scaled amplitudes before fitting
   ! MaxAmplFitPercent: Maximum percentile of sources in overflow bin
   !Use Tracks, only : a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3
   ! extracted parameters of designated fit functions
   ! Presently: 1== (modified) exponential; 2== (modified)powerlaw; 3=not used
   !
   ! internal AmplitudeFit parameters =======================
   !Use Tracks, only : Selection_sort, i_AmplTh, Ampl_Hist, FitFunc
   ! i_AmplTh: Index of first bin in scaled amplitude histogram that is fitted.
   ! Ampl_Hist: actual histogram of amplitudes
   ! FitFunc: Function that is used
   Use DS_Select, only : Label, MaxAmplFitPercent, Nmax_Ampl, Ampl_Hist, Ampl_Hist_W, AmplScale, FitFunc
   Use DS_Select, only : Nrm, a1,b1,c1,d1,ChiSq1,Nrm1, a2,b2,c2,d2,ChiSq2,Nrm2, a3,b3,c3,d3,ChiSq3 ! parameters for describing the amplitude distribution
   Use DS_Select, only : i_AmplTh, Ampl
   Use unque, only : Selection_sort
   IMPLICIT none
   Integer, intent(in) :: NrSources
   Integer, intent(in), optional :: SourceNrs(1:NrSources)
   character(len=*), intent(in), optional :: SelFileName
   Integer :: SourceAmp(1:NrSources)  !  Amplitudes of sources in the i^th track
   !
   Integer :: i,k, i_Ampl, i_AmplMax, i_cut, i_MaxFit,N_bins
   Integer ( kind = 4 ) :: meqn, nvar  ! Number of data points to fit & number of parameters
   Character*8 :: extension
   Real*8 :: Max_Ampl, Ampl10, X(4),ChiSQDF, AmplBin, Counts2BinValue, ZeroErr=1., HistNrm ! 0.25
   Real*8 :: Fit_exp(Nmax_Ampl), Fit_pow(Nmax_Ampl),a,b,c, d_AmplScale
   Real*8, save :: Tiny=epsilon(Tiny)
   !
   !write(extension,"(A2,i1,A4)") '_s',nxx,'.dat' !,&
   !write(extension,"(i1,A4)") i,'.dat' !,&
   !OPEN(unit=27,FILE='files/'//trim(datfile)//'_??_'//trim(extension),FORM='FORMATTED',STATUS='unknown')
   !write(27,"('! nr', 4(1x,A14),3x,3(1x,A12),2x,A10)") 't [ms]','E [km]','N','h','v_E','v_N','v_h','v [10^6m/s]'
   !
   !If(present(SelFileName)) write(2,*) 'A; SelFileName:"',SelFileName,'"'
   If(present(SourceNrs)) Then
      Do i=1,NrSources
         SourceAmp(i)=Label(2,SourceNrs(i))
      Enddo
   Else
      SourceAmp(1:NrSources)=Label(2,1:NrSources)
      write(2,"(A,I6,A)") 'Fit amplitude spectrum of all',NrSources,' sources'
   EndIf
   CALL Selection_sort(SourceAmp)  ! Re-order  sources according to pulse amplitude
   !
   Fit_exp(:)=0.0
   Fit_pow(:)=0.0
   N_bins=-1
   !
   If(MaxAmplFitPercent.lt.0.) MaxAmplFitPercent=0.1
   Max_Ampl=SourceAmp(NrSources)/AmplScale
   Ampl10=SourceAmp(NrSources*9/10)/AmplScale  ! 10 percentile amplitude
   write(2,"(A,I8,A,F9.1,A,F10.2,A,F9.2,A,I4)") 'Amplitudes (#=',NrSources,'), max@', Max_Ampl, &
      ', 5-pctile@', SourceAmp(NrSources*95/100)/AmplScale, ', 10-pctile@', Ampl10
   i_cut=NrSources*(1.-MaxAmplFitPercent/100)
   If(i_cut.ge.0) then
      If(i_cut.eq.0) i_cut=1
      write(2,"(A,F7.3,A,F10.2,A,I4,F8.2)") 'Amplitude with ', MaxAmplFitPercent, &
         'pctile as included in fit @', SourceAmp(i_cut)/AmplScale,', i_cut=', i_cut, SourceAmp(1)/AmplScale
   Else
      i_cut=NrSources  !  No sources in fit
   Endif
   !
   d_AmplScale=1.1  ! factor for log bins
   ! Automatic scaling
   If(NrSources.gt.50) Then
      N_bins=Nmax_Ampl-2  ! -2 for a stable fit with analytic function
      If(NrSources.lt.N_bins) N_bins=NrSources/10
      d_AmplScale=log(1.d0*SourceAmp(NrSources)/SourceAmp(1))/N_bins
      d_AmplScale=exp(d_AmplScale)
      !write(2,*) 'Automatic amplitude histogram scaling;', d_AmplScale, N_bins, &
      !      SourceAmp(1)/AmplScale, SourceAmp(NrSources)/AmplScale
      !write(2,*) ';', d_AmplScale**N_bins, 1.* SourceAmp(NrSources)/SourceAmp(1)
   EndIf
   !
   AmplBin=SourceAmp(1)/d_AmplScale  ! scaled up by factor 1/AmplScale, always start in the first bin
   i_Ampl=1
   i_MaxFit=0
   Ampl_Hist(:)=0
   HistNrm=0.  ! =\int I N(I) dI
   !Ampl10=0.
   !write(2,*) 'AmplitudeFit: SourceAmp(1)', SourceAmp(1), AmplBin, i_cut
   Do k=1,NrSources
      Do While( SourceAmp(k) .gt. AmplBin)  ! start from the smallest amplitude and fill up bins
         Counts2BinValue=(AmplBin*(d_AmplScale-1.)/AmplScale)
         Ampl_Hist_W(i_Ampl)= Counts2BinValue/sqrt(Ampl_Hist(i_Ampl)+ZeroErr)   ! Calculate fitting weight for previous bin
         Ampl(i_Ampl)=AmplBin*(1.+d_AmplScale)/(AmplScale*2.)  ! In true units
         Ampl_Hist(i_Ampl)=Ampl_Hist(i_Ampl)/Counts2BinValue
         HistNrm=HistNrm + Ampl_Hist(i_Ampl)*Ampl(i_Ampl)*AmplBin*(d_AmplScale-1.)/AmplScale
         !Ampl10=Ampl10 + Ampl_Hist(i_Ampl)/AmplScale*AmplBin*(d_AmplScale-1.)  ! equals tot nr of sources, as it should
         ! chi^2 = sum_i (Ampl_Hist(i) - model) * Ampl_Hist_W(i)
         !write(2,*) k, i_Ampl, AmplBin,Ampl(i_Ampl), Ampl_Hist(i_Ampl)*Counts2BinValue, Ampl_Hist_W(i_Ampl)/Counts2BinValue
         AmplBin=AmplBin*d_AmplScale
         !write(2,*) 'AmplitudeFit: SourceAmp(k)',k, SourceAmp(k),i_Ampl, AmplBin
         If(i_cut .gt. k) i_MaxFit=i_Ampl+1   !  Still include this bin in fitting
         i_Ampl=i_Ampl+1  ! start to fill next bin
         If(i_Ampl .ge. Nmax_Ampl) goto 9
      Enddo
      Ampl_Hist(i_Ampl)=Ampl_Hist(i_Ampl)+1.
   EndDo
9  Continue
   !write(2,*) NrSources, SourceAmp(NrSources) , AmplBin, Ampl_Hist(i_Ampl), i_Ampl, Ampl10, HistNrm
   Counts2BinValue=(AmplBin*(d_AmplScale-1.)/AmplScale)
   Ampl_Hist_W(i_Ampl)= Counts2BinValue/sqrt(Ampl_Hist(i_Ampl)+ZeroErr)   ! Calculate fitting weight for previous bin
   Ampl(i_Ampl)=AmplBin*(1.+d_AmplScale)/(2.*AmplScale)  ! In true units
   Ampl_Hist(i_Ampl)=Ampl_Hist(i_Ampl)/Counts2BinValue
   !write(2,*) N_bins, i_Ampl, Ampl(i_Ampl), Ampl_Hist(i_Ampl), Ampl_Hist_W(i_Ampl)
   i_AmplMax=i_Ampl
   If(i_Ampl .ge. Nmax_Ampl) Then
      Ampl_Hist(i_Ampl)=Ampl_Hist(i_Ampl)+NrSources-k  !  Fill out the highest bin with remainder, when necessary
      i_AmplMax=Nmax_Ampl-1
      write(2,*) 'i_Ampl .ge. Nmax_Ampl', i_Ampl,Nmax_Ampl, AmplBin,Ampl(i_Ampl)
   EndIf
   If(i_MaxFit.eq.0) i_MaxFit=i_AmplMax
   !write(2,*) i_AmplMax, NrSources, ', i_MaxFit:',i_MaxFit
   !write(2,*) Ampl(i_AmplMax), SourceAmp(NrSources), AmplBin
   !
   !If(present(SelFileName)) write(2,*) '(B,SelFileName:"',SelFileName,'"'
   !flush(unit=2)
   !
   !write(2,*) 'Ampl_Hist',Ampl_Hist(1:20)
   !write(2,*) 'bins of max etc:',Max_Ampl, i_AmplZero, i_AmplMax, d_Ampl,Ampl_Hist(11)
   !write(2,*) d_Ampl,AmplScale,d_AmplScale
   a1=0.; b1=0.; c1=0.
   !goto 5
   !
   i_AmplTh=1 ! first index for fitting in histogram
   Fit_exp(1:i_AmplTh)=0.
   Meqn = i_MaxFit-i_AmplTh+1      ! number of equations
   !FitFunc=0  ! b*exp(-a*Ampl)  a=x(1)/d_Ampl
   !FitFunc=2 ; nvar= 3 ! b*exp(-a*A-c/A)  a=x(1)/d_Ampl
   FitFunc=3 ; nvar= 3 ! b*exp(-a*A-c/A^2)  a=x(1)/d_Ampl
   x(1)=1/20.   ! CalcN = exp(-X_p(1)*Ampl+X_p(2))
   X(2)=log(Ampl_Hist(i_AmplTh))+ X(1)*ampl(i_AmplTh)
   X(3)=0.
   Call FitAmpl(Meqn, nvar, X,ChiSq1, Fit_exp(i_AmplTh) )
   a1=X(1)
   b1=(x(2))
   c1=x(3)
   d1=0.
   If(Meqn.lt.0) Then ! no fitting possible
      write(2,*) 'no fitting possible of a1--d1'
      !flush(unit=2)
      goto 5
   EndIf
   b1=exp(x(2))
   !goto 8
   If(c1.lt.0.) then
      i_AmplTh=4 ! first index for fitting in histogram
      Fit_exp(1:i_AmplTh)=0.
      Meqn = i_MaxFit-i_AmplTh+1      ! number of equations
      FitFunc=0 ; nvar= 2 ! exp(b-a*A)  a=x(1)/d_Ampl
      x(1)=1/20.   ! CalcN = exp(-X_p(1)*Ampl+X_p(2))
      X(2)=log(Ampl_Hist(i_AmplTh))+ X(1)*ampl(i_AmplTh)
      X(3)=0.
      Call FitAmpl(Meqn, nvar, X,ChiSq1, Fit_exp(i_AmplTh) )
      a1=X(1)
      b1=exp(x(2))
      c1=0
   Endif
   !Goto 1 ! --------------------
   !
5   FitFunc=5 ; nvar= 3  ! F(A)=b (A)^(-a)*exp(-c/A)
   !In GLE: d8 = d_bin*b2*(x)^(-a2)*exp(-c2/x-d2/x^2) with x is the amplitude
   i_AmplTh=1 ! first index for fitting in histogram
   Fit_pow(1:i_AmplTh)=0.
   Meqn = i_MaxFit-i_AmplTh+1      ! number of equations
   X(1)=2.
   X(2)=Ampl_Hist(i_AmplTh) * (ampl(i_AmplTh)**X(1))
   X(3)=0.
   Call FitAmpl(Meqn, nvar, X,ChiSq2, Fit_pow(i_AmplTh) )
   a2=X(1)
   b2=x(2)!*d_Ampl**(X(1)-1)
   c2=x(3)!*d_Ampl
   d2=0.
   If(Meqn.lt.0) Then ! no fitting possible
      write(2,*) 'no fitting possible of a2--d2'
      !flush(unit=2)
      goto 8
   EndIf
   If(c2.lt.0.) then
1     FitFunc=1 ; nvar=2 ! F(A)=b (A)^(-a) with   where A is true amplitude
      !      CalcN = -X(2) * (i_Ampl)**(-X(1))
      !In GLE: d8  = d_bin*b2 * (x)^(-a2)        with x is the amplitude
      i_AmplTh=4  ! first index for fitting in histogram
      Fit_pow(1:i_AmplTh)=0.
      Meqn = i_MaxFit-i_AmplTh+1      ! number of equations
      X(1)=4.
      X(2)=Ampl_Hist(i_AmplTh) * (ampl(i_AmplTh)**X(1))
      x(3)=0.
      !write(2,*) ' initial guess: ',X(1),X(2), Meqn
      Call FitAmpl(Meqn, nvar, X,ChiSq2, Fit_pow(i_AmplTh) )
      a2=X(1)
      b2=x(2)
      c2=0.
      d2=0.
      !write(2,*) X(1),X(2), a2, b2
   EndIf
   !
8  continue
   !   write(2,*) 'entering SelFileName', present(SelFileName)
   !   flush(unit=2)
   If(present(SelFileName) .and. (N_bins.gt.1) .and. (SelFileName.ne."")) then
      ! write obtained results to file for plotting
      !write(2,*) (d_AmplScale-1.),Ampl(1)
      !write(2,*) 'entering pre-logscaling', N_bins
      !flush(unit=2)
      a=Ampl(1) ; b=Ampl(N_bins)  ! set scale for plots to cover at least factor 10
      !write(2,*) 'AmplFit: entering logscaling', a,b, N_bins
      !flush(unit=2)
      If(b/a.lt.10.) Then
         c=sqrt(10.*a/b)
         a=a/c ; b=b*c
      EndIf
      c=log(a)/log(10.)
      !write(2,*) 'c=log(a)/log(10.):',c,a
      If(c.lt.0) Then
         i=int(c)-1 ;
      Else
         i=int(c)
      EndIf
      c=10.**i  ! next lower power of 10
      If(a.gt.5.*c) Then
         a=5.*c
      ElseIf(a.gt.3.*c) Then
         a=3.*c
      ElseIf(a.gt.2.*c) Then
         a=2.*c
      Else
         a=c
      EndIf
      !write(2,*) i,c,a
      If(a.lt.0.01) a=0.01
      c=log(b)/log(10.)
      !write(2,*) 'c=log(b)/log(10.):',c,b
      If(c.lt.0) Then
         i=int(c)-1 ;
      Else
         i=int(c)
      EndIf
      c=10.**i  ! next lower power of 10
      If(b.lt.2.*c) Then
         b=2.*c
      ElseIf(b.lt.3.*c) Then
         b=3.*c
      ElseIf(b.lt.5.*c) Then
         b=5.*c
      Else
         b=10.*c
      EndIf
      !write(2,*) 'AmplFit: end logscaling', a,b, c, i, 'SelFileName="',TRIM(SelFileName), '"'
      !  End log amplitude scale
      Nrm1=0.  ; Nrm2=0.
      !write(2,*) 'N_bins:',N_bins
      OPEN(unit=28,FILE=TRIM(SelFileName)//'.dat', FORM='FORMATTED',STATUS='unknown',ACTION='write') ! space separated values
      write(28,"(1x,F8.2,2F12.2, 6(1x,g10.4),f7.3,1x,A,A)") a, Ampl(i_MaxFit), b, &
         a1,b1,c1,a2,b2,c2,d2, TRIM(SelFileName),' ! by AmplitudeFit' ! gleRead:  NTracks SourcTotNr Q t_start label$ AmplitudePlot a1 b1 c1 a2 b2 c2 d2
      Do i_Ampl=1,i_AmplMax !N_bins ! i_MaxFit ! i_AmplMax  !  finally write all accepted source to file
         !Ampl10=(Ampl(1)/NrSources*(d_AmplScale-1.))*Ampl(i_Ampl)*Ampl(i_Ampl)
         Ampl10= Ampl(i_Ampl)*Ampl(i_Ampl)/HistNrm  ! norm changed 18 Lan 2023
         !Ampl10=1./(NrSources*Ampl(1)) * Ampl(i_Ampl)*Ampl(i_Ampl)
         !write(2,"(1x,I5,7(1x,g14.6))") i_Ampl, Ampl(i_Ampl), Ampl10*Ampl_Hist(i_Ampl), Ampl10, Ampl_Hist(i_Ampl)
         write(28,"(1x,I5,7(1x,g14.6))") i_Ampl, Ampl(i_Ampl), Ampl10*Ampl_Hist(i_Ampl),Ampl10/Ampl_Hist_W(i_Ampl), &
               Ampl10*Fit_exp(i_Ampl)+Tiny, Ampl10*Fit_pow(i_Ampl)+Tiny, Fit_pow(i_Ampl)
         Nrm1=Nrm1+Fit_exp(i_Ampl)*(Ampl(i_Ampl+1)-Ampl(i_Ampl))
         Nrm2=Nrm2+Fit_pow(i_Ampl)*(Ampl(i_Ampl+1)-Ampl(i_Ampl))
      enddo
      close(unit=28)
      Write(2,*) 'File written:',TRIM(SelFileName)//'.dat',' with ',i_AmplMax, ' lines'
   EndIf
   Nrm=NrSources*Ampl(1)
   write(2,"(1x,A, f6.2, A, f7.4, g10.3, g11.3, A,2g11.3)") '$b *exp(-a*A-c/A^2); \chi^2=$',ChiSq1, &
                                    ',with  a,b,c= ',a1,b1/Nrm,c1,'; nrm=', Nrm ! , HistNrm
   write(2,"(1x,A, f6.2, A, f7.3, g10.3, g11.3, A,2g11.3)") '$b *A^-a *exp(-c/A); \chi^2=$',ChiSq2, &
                                    ',with  a,b,c= ',a2,b2/Nrm,c2,'; nrm=', Nrm
   Return
   !
4   FitFunc=4 ; nvar= 3  ! F(A)=b A^(-a)*exp(-c/A^2) with i_Ampl=A/d_Ampl  where A is plotted
   i_AmplTh=3 !SourceAmp(NrSources*1/10)*d_AmplScale  ! Threshold for fitting in histogram index
   Meqn = Max_Ampl-i_AmplTh-1      ! number of equations
   X(1)=4.
   X(2)=Ampl_Hist(i_AmplTh)*(i_AmplTh/10.)**X(1)
   !Call FitAmpl(Meqn, nvar, X,ChiSQDF)
   a2=X(1)
   b2=x(2)
   d2=x(3)
   c2=0.
   Return
   !
3   FitFunc=3 ; nvar= 3 ! b*exp(-a*A-c/A^1.5)  a=x(1)/d_Ampl
   !Call FitAmpl(Meqn, nvar, X,ChiSQDF)
   a2=X(1)
   b2=exp(x(2))
   c2=x(3)
   d2=0.
   !
!   Return
   !write(2,*) MeanTrackLoc(1,1:t_MTL)
   !Call Flush(2)
   !close(unit=27)
    !
   Return
End Subroutine AmplitudeFit
!===============================================
Subroutine FitAmpl(Meqn, nvar, X,ChiSQDF, BestFit)
   !Use AmpFit, only :  jacobianampl
   Use DS_Select,  only : i_AmplTh, Nmax_Ampl, AmplScale, Ampl_Hist, FitFunc
   ! i_AmplTh=first entry that will be included in the fit
    Implicit none
    Integer( kind = 4 ), intent(inout) :: Meqn
    integer ( kind = 4 ),intent(in) :: nvar    ! number of parameters
    real ( kind = 8 ), intent(inout) :: X(4)
    real ( kind = 8 ), intent(out) :: ChiSQDF
    integer ( kind = 4 ) :: i, j, k
    integer :: v_dim, error
    !integer ( kind = 4 ) :: meqn  ! Number of data points to fit
    integer ( kind = 4 ) iv(60+nvar)
    !External CompareAmpl  ! subroutine that compares FitFunction to Data
    !External JacobianAmpl
    external ufparm  ! Dummy external routine
    integer ( kind = 4 ) uiparm(1)    ! Not really used
    real ( kind = 8 ) :: Jacobian(Meqn,nvar), BestFit(Meqn)
    !real ( kind = 8 ) v(v_dim) ! dim=93 + n*p + 3*n + p*(3*p+33)/2 ,n-meqn, p=nvar
    real ( kind = 8 ),allocatable :: v(:) !  in NL@SOL: real ( kind = 8 ) v(93 + n*p + 3*n + (p*(3*p+33))/2)
    integer ( kind = 4 ) :: NF
    !
    !     write(2,*) 'entering fitting', nvar, Meqn
    !      flush(unit=2)
    If(nvar.ge.Meqn) then
      Write(2,*) '****** too few values to fit, ',Meqn, FitFunc
      Meqn=-1
      !    flush(unit=2)
      Return
    EndIf
    v_dim=93 + Meqn*nvar + 3*Meqn + nvar*(3*nvar+33)/2
    !allocate( Jacobian(Meqn,nvar) )
    allocate( v(v_dim) )
    !
    call dfault( iv, v)
    iv(1) = 12 ! 12= do not call dfault again
    iv(14) = 0 ! 1: means print a covariance matrix at the solution.
    iv(15) = 0 ! if = 1 or 2, then a finite difference hessian approximation h is obtained.
               ! if =0 nothing is done
               ! if positive: with step sizes determined using v(delta0=44), a multiplicative factor)
               ! If negative: then only function values are used with step sizes determined using v(dltfdc=40)
    iv(19) = 0 ! controls the number and length of iteration summary lines printed
    iv(21) = 0 !2 ! is the output unit number on which all printing is done.
    iv(22) = 0!1 ! print out the value of x returned (as well as the corresponding gradient and scale vector d).
    iv(23) = 1 ! print summary statistics upon returning.
    iv(24) = 0 ! print the initial x and scale vector d
    v(32) =1.d-4 ! is the relative function convergence tolerance
    v(40) =1.40d-3 ! the step size used when computing the covariance matrix when iv(covreq=15) = -1 or -2, step size = v(dltfdc=40) * max(abs(x(i)), 1/d(i))
    ! v(44) =1.d-3 ! the factor used in choosing the finite difference step size used in computing the covariance matrix when
                !    iv(covreq=15) = 1 or 2, step size = v(delta0=44) * max(abs(x(i)), 1/d(i)) * sign(x(i))
    !If(CalcHessian) then
    !  iv(14) = 1
    !  iv(15) = 2
    !  v(32) =1.d-8 ! is the relative function convergence tolerance
    !endif
    !
    ! iv(nfcall)... iv(6) is the number of calls so far made on calcr (i.e.,
    !             function evaluations, including those used in computing
    !             the covariance).
    ! iv(mxfcal)... iv(17) gives the maximum number of function evaluations
    !             (calls on calcr, excluding those used to compute the co-
    !             variance matrix) allowed.  if this number does not suf-
    !             fice, then nl2sol returns with iv(1) = 9.  default = 200.
    ! iv(mxiter)... iv(18) gives the maximum number of iterations allowed.
    !             it also indirectly limits the number of gradient evalua-
    !             tions (calls on calcj, excluding those used to compute
    !             the covariance matrix) to iv(mxiter) + 1.  if iv(mxiter)
    !             iterations do not suffice, then nl2sol returns with
    !             iv(1) = 10.  default = 150.
    !    iv(18)=1
    !  iv(26) if (iv(covmat) is positive, then the lower triangle of the covariance matrix is stored rowwise in v starting at
    !             v(iv(covmat)).  if iv(covmat) = 0, then no attempt was made.
    ! iv(nfcov).... iv(40) is the number of calls made on calcr when
    !             trying to compute covariance matrices.
    ! iv(ngcall)... iv(30) is the number of gradient evaluations (calls on
    !             calcj) so far done (including those used for computing
    !             the covariance).
    ! iv(ngcov).... iv(41) is the number of calls made on calcj when
    !             trying to compute covariance matrices.
    ! iv(niter).... iv(31) is the number of iterations performed.
    !
   !call nl2sno ( meqn, nvar, x, CompareAmpl, iv, v, uiparm, Jacobian, ufparm )
   flush(unit=2)
   Call nl2sol ( meqn, nvar, x, CompareAmpl, JacobianAmpl, iv, v, uiparm, Jacobian, ufparm, error )
   If(error.gt.0) then
      write(2,*) 'endless loop in NL2SOL truncated'
      goto 9
   endif
   !
   ChiSQDF=2*v(10)/(meqn-nvar)
   NF=-1  ! Store fitted function in R
   call CompareAmpl ( meqn, nvar, X, NF, BestFit, uiparm, Jacobian, ufparm )
9  Continue
    !
    !Deallocate( Jacobian)
    Deallocate( v )
    !
    ! write(2,"(A,I2,A,I2,A,F7.2,A,4G12.5)") 'Final result fit function=', FitFunc, &
    !  ', First amplitude=', i_AmplTh,', Chi^2=', ChiSQDF,', Parameters=', X(1:nvar)
    Return
    !
End Subroutine FitAmpl
!=================================================
Subroutine JacobianAmpl ( meqn, nvar, X_p, nf, Jacobian, uiparm, urparm, ufparm )
  implicit none
  real ( kind = 8 ), intent(in) :: X_p(1:nvar)
  integer ( kind = 4 ), intent(in) :: meqn
  integer ( kind = 4 ), intent(in) :: nvar
  integer ( kind = 4 ), intent(in) :: nf
  external ufparm
  integer ( kind = 4 ), intent(in) :: uiparm(*)
  real ( kind = 8 ), intent(in) :: urparm(*)
  real ( kind = 8 ), intent(out) :: Jacobian(meqn,nvar) ! in calling this routine this space has been reserved as part of v
  real ( kind = 8 ) :: R(meqn), Jacob(meqn,nvar)
  !
  !write(2,*) 'X_p',X_p(1:3)
  !write(*,*) 'calc Jac'
  Call CompareAmpl ( meqn, nvar, X_p, nf, R, uiparm, Jacob, ufparm )
  !write(*,*) 'JacobN',meqn,Jacob(meqn,:)
  Jacobian(:,:)=Jacob(:,:)
  Return
  End Subroutine JacobianAmpl
!=================================================
Subroutine CompareAmpl ( meqn, nvar, X_p, nf, R, uiparm, Jacobian, ufparm )
!*****************************************************************************80
!
!  Discussion:
!    Given the value of the vector X, this routine computes the value of the residues R=(t_measured-t_X) used for chi^2 fitting
!  Author:
!    Olaf Scholten
!
!  Parameters:
!    Input, integer ( kind = 4 ) MEQN, the number of functions.
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!    Input, real ( kind = 8 ) X(NVAR), the current value of the variables.
!    Input, integer ( kind = 4 ) NF, the number of times the residual routine
!       has been called so far.
!    Output, real ( kind = 8 ) R(MEQN), the residual vector, that is, the
!       value of the functions for the given input value of the variables.
!    Input, integer ( kind = 4 ) UIPARM(*), a user array.
!    Input, real ( kind = 8 ) URPARM(*)==Jacobian, a user array.
!    Input, external UFPARM, an external reference to a user subroutine
!       or function.
!
   Use DS_Select,  only : i_AmplTh, Ampl, Ampl_Hist, Ampl_Hist_W, FitFunc
!
  implicit none
    real ( kind = 8 ), intent(in) :: X_p(1:nvar)
  integer ( kind = 4 ), intent(in) :: meqn
  integer ( kind = 4 ), intent(in) :: nvar
  integer ( kind = 4 ), intent(in) :: nf
  real ( kind = 8 ), intent(out) :: R(meqn)
  external ufparm
  integer ( kind = 4 ), intent(in) :: uiparm(*)
  real ( kind = 8 ), intent(out) :: Jacobian(meqn,nvar)  ! derivative of R w.r.t. parameters
  real ( kind = 8 ) :: X(nvar), A, CalcN, Weight, chisq
  real ( kind = 8 ), parameter :: ZeroErr=0.5
   !
   integer :: i,j,k, i_eqn
   !
   ChiSq=0.
   Select Case(FitFunc)
      Case (5)   ! F(A)=b A^(-a)*exp(-c/A) = exp(-a ln(A) +b' -c/A)
         Do i_eqn=1,meqn
            k=i_AmplTh-1+i_eqn
            Weight = Ampl_Hist_W(k)
            A= Ampl(k)
            CalcN = A**(-X_p(1)) * exp(-X_p(3)/A)
            R(i_eqn)= (Ampl_Hist(k) -X_p(2)* CalcN)*Weight
            ChiSq = ChiSq + R(i_eqn)*R(i_eqn)
            !
            Jacobian(i_eqn,1) =  +X_p(2)*CalcN *Weight*log(A)
            Jacobian(i_eqn,2) =  - CalcN *Weight
            Jacobian(i_eqn,3) =  + X_p(2)*CalcN *Weight/A
            !If(nf.lt.0) write(2,*) i_eqn,k, A, Ampl_Hist(k), X_p(2)* CalcN, R(i_eqn)*Weight, Weight
            If(nf.lt.0) R(i_eqn)=X_p(2)* CalcN
         Enddo  ! i_eqn=1,meqn
      Case (4)   ! F(A)=b (A+1)^(-a)*exp(-c/A^2)
         Do i_eqn=1,meqn
            k = i_AmplTh-1+i_eqn
            Weight = Ampl_Hist_W(k)
            A = Ampl(k)
            CalcN = A**(-X_p(1)) * exp(-X_p(3)/(A*A))
            R(i_eqn)= (Ampl_Hist(k) -X_p(2)* CalcN)*Weight
            ChiSq = ChiSq + R(i_eqn)*R(i_eqn)
            !
            Jacobian(i_eqn,1) =  +X_p(2)*CalcN *Weight*log(A)
            Jacobian(i_eqn,2) =  - CalcN *Weight
            Jacobian(i_eqn,3) =  + X_p(2)*CalcN *Weight/(A*A)
            If(nf.lt.0) R(i_eqn)=X_p(2)* CalcN
         Enddo  ! i_eqn=1,meqn
      Case (3)     !  exp(-x(1)*A+b-x(3)/A^1.5)
         Do i_eqn=1,meqn
            k = i_AmplTh-1+i_eqn
            Weight = Ampl_Hist_W(k)
            A = Ampl(k)
            CalcN = exp(-X_p(1)*A+X_p(2)-X_p(3)/(A*A))
            R(i_eqn)= (Ampl_Hist(k) -CalcN)*Weight
            !ChiSq = ChiSq + R(i_eqn)*R(i_eqn)
            !
            Jacobian(i_eqn,1) =  +A*CalcN *Weight  ! =dR/da
            Jacobian(i_eqn,2) =  - CalcN *Weight      ! =dR/db
            Jacobian(i_eqn,3) =  + CalcN *Weight/(A*A)      ! =dR/dc
            !If(nf.lt.0) write(2,*) i_eqn,k, A, Ampl_Hist(k), CalcN, R(i_eqn)*Weight, Weight
            If(nf.lt.0) R(i_eqn)= CalcN
         Enddo  ! i_eqn=1,meqn
      Case (2)
         Do i_eqn=1,meqn !  exp(-x(1)*A+b-x(3)/A)
            k = i_AmplTh-1+i_eqn
            Weight = Ampl_Hist_W(k)
            A = Ampl(k)
            CalcN = exp(-X_p(1)*A+X_p(2)-X_p(3)/A)
            R(i_eqn)= (Ampl_Hist(k) -CalcN)*Weight
            !
            Jacobian(i_eqn,1) =  +A*CalcN *Weight  ! =dR/da
            Jacobian(i_eqn,2) =  - CalcN *Weight      ! =dR/db
            Jacobian(i_eqn,3) =  + CalcN *Weight/A      ! =dR/dc
            If(nf.lt.0) R(i_eqn)= CalcN
         Enddo  ! i_eqn=1,meqn
      Case (1)  ! powerlaw
         Do i_eqn=1,meqn
            k = i_AmplTh-1+i_eqn
            Weight = Ampl_Hist_W(k)
            A = Ampl(k)
            CalcN = A**(-X_p(1))
            R(i_eqn)= (Ampl_Hist(k) -X_p(2)* CalcN)*Weight
            ChiSq = ChiSq + R(i_eqn)*R(i_eqn)
            !
            Jacobian(i_eqn,1) =  +X_p(2)*CalcN *Weight*log(A)
            Jacobian(i_eqn,2) =  - CalcN *Weight
            If(nf.lt.0) R(i_eqn)=X_p(2)* CalcN
         Enddo  ! i_eqn=1,meqn
      Case Default          !  exp(-x(1)*A+b)
         Do i_eqn=1,meqn !  exp
            k = i_AmplTh-1+i_eqn
            Weight = Ampl_Hist_W(k)
            A = Ampl(k)
            CalcN = exp(-X_p(1)*A+X_p(2))
            R(i_eqn)= (Ampl_Hist(k) -CalcN)*Weight
            !ChiSq(i_stat,i_Peak) = ChiSq(i_stat,i_Peak) + R(i_eqn)*R(i_eqn)
            !
            Jacobian(i_eqn,1) =  +A*CalcN *Weight
            Jacobian(i_eqn,2) =  - CalcN *Weight
            !If(nf.lt.0) write(2,*) i_eqn,k, A, Ampl_Hist(k), CalcN, R(i_eqn)*Weight, Weight
            If(nf.lt.0) R(i_eqn)= CalcN
         Enddo  ! i_eqn=1,meqn
   End Select
   !write(2,*) 'ChiSq:',ChiSq/(meqn-nvar),nf,X_P(1:nvar)
   !stop
   return
End Subroutine CompareAmpl
!==================================================
Subroutine ufparm ( meqn, nvar, x )
!*****************************************************************************80
!! UFPARM is a user-supplied external routine.
!
!  Discussion:
!    The name of the routine, the argument list, and even whether
!       it is a function or subroutine, are left to the user.
!    NL2SOL simply passes the external reference from the calling
!       program through to the residual and jacobian routines.
!    If the user has no need for this facility, then a dummy
!       routine like this one may be used.
!
!  Modified:
!    07 February 2003
!
!  Parameters:
!    Input, integer ( kind = 4 ) MEQN, the number of functions.
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!    Input, real ( kind = 8 ) X(NVAR), the current value of the variables.
!
  implicit none
  integer ( kind = 4 ) meqn
  integer ( kind = 4 ) nvar
  real ( kind = 8 ) x(nvar)
  return
End Subroutine ufparm
End Module AmpFit
