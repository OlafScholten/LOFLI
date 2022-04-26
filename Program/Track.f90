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
    Include 'TrackModule.f90'
    Include 'Fitter_Ampl.f90'
    Include 'FFT_routines.f90'
    Include 'System_Utilities.f90'
    Include 'GLEplotUtil.f90'
    !
PROGRAM rewriteEvent
! copied from C:\Users\Olaf Scholten\Documents\AstroPhys\Lightning\Imaging\LOFAR-MAP
   Use Tracks, only : PreDefTrackNr, PreDefTrackFile
   Use Tracks, only : TrackNr, LongTrackNr, LongTrack_Min, TrackNrMax, NLongTracksMax, TrackLenMax, TrackENr, TrackE, TrackNrLim
   Use Tracks, only : Wtr, MaxTrackDist, TimeWin, HeightFact, datfile !,  WriDir
   Use Tracks, only : RA, maxd, EventNr, Label, HPSORT_mult_RI
   Use Tracks, only : AmplScale, d_AmplScale, AmplitudeFit
   Use Tracks, only : a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3
   use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, DataFolder, FlashName, Windows, RunMode
   use GLEplots, only : GLEplotControl
   IMPLICIT none
   character*150 :: SystemCommand
   Character*8 :: extension
   character*100 :: pars
   Character*20 :: Utility, release
  Character(len=180) :: lineTXT
  integer :: i,j,k, nxx, d_Ampl  ! unit,
  Real*8 :: Qual, t_start, dt_MTL, AmplitudePlot
  integer :: kk, jk, wrunt, DelNEff, N_EffAnt !, Max_EffAnt
  Logical :: EndInput
  !
  !unit=12
  d_Ampl=1  ! bin size for amplitude histogram
   If(d_Ampl.le.1) d_Ampl=1
   d_AmplScale=AmplScale/d_Ampl
  !OPEN(UNIT=unit,STATUS='REPLACE',FILE='RunGLEplots.bat')
  OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='FlashImage.out')
  !write(2,*) 'Program RewriteEvent alias Track-v19 alias FlashImage'
  Utility='Track&Qual.Control'
  release='v21'
  Call System_Initiation(Utility, release, ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows)  ! set Flash name & folder names
  RunMode=5  !  used in the gle-plotting routine
  !
  !WriDir=TRIM(DataFolder)
  !
  Do
      !
      ! Read in selection criteria from stdin and apply to located sources
      ! Unit=29 is opened
      Call SelSources(Qual, t_start, EndInput)
      If(EndInput) exit
      If(EventNr.le.1) then
        write(2,*)'EventNr=',EventNr
        Read(*,*) MaxTrackDist, Wtr, TimeWin
        close(unit=29)
        write(2,*) '======================================= quick'
        cycle  !  Read data for next plot
      endif
      !
      !
      CALL HPSORT_mult_RI(RA,Label,EventNr)  ! RA(1:4,*) are the accepted source locatopns
      ! write(*,*) 'Sorting Done'
      !
      !-------------- Start putting sources on tracks ----------------------
      !
      MaxTrackDist= 0.1 ! [km]
      Wtr= 0.5      ! Weight for next point on track
      Read(*,*,iostat=nxx) MaxTrackDist, Wtr, TimeWin, PreDefTrackFile, dt_MTL, HeightFact, NLongTracksMax, AmplitudePlot !  dt_MTL should be about 2*TimeWin to generate comparable results
      If(nxx.ne.0) then
         dt_MTL=0.5
         HeightFact=3.
         NLongTracksMax=0
         AmplitudePlot=-1. ! do not plot amplitude info
         write(2,*) '** No value given to time for calculating Mean Track Location or HeightFact ****'
      Endif
      If(TimeWin .lt. 0.) TimeWin = ((RA(1,EventNr)-RA(1,1))*10./EventNr) ! *((RA(1,EventNr)-RA(1,1))*10./EventNr)
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
         Call Assign2Tracks
      EndIf
      ! TrackE(k,i)=j :  Source TrackE(k,i) is the k^th source assigned to track i
      ! TrackENr(i) : Total number of sources assigned to track i
      !
      !  unit=29 is opened in    Call SelSources(Qual, t_start)
      !OPEN(unit=29,FILE=trim(datfile)//'.dat',FORM='FORMATTED',STATUS='unknown')  ! will contain all accepted sources
      !write(29,"(6(1x,f8.3),2(1x,f8.3),1x,A,f9.3,A)") xMin, xMax, yMin,yMax,zMin,zMax,tMin,tMax,trim(ZoomBox),tp*1000.,'  ! '
      !
      Call AmplitudeFit(EventNr)
      !write(2,*) 'a1,b1,c1,a2,b2,c2',a1,b1,c1,a2,b2,c2
      Flush(unit=2)
      !
      write(29,"(1x,I3,I8,F7.2,F8.4,2x,A,1x,F6.2, 2(f8.3,f10.0,f8.3),f8.3,1x,A)") LongTrackNr, EventNr, Qual, t_start, &
          TRIM(FlashName)//'/'//TRIM(datfile), AmplitudePlot, a1,b1,c1,a2,b2,c2,d2, ' !'
      Do j=1,EventNr  !  finally write all accepted source to file
         write(29,"(1x,i8,4(2x,g14.8),3x,F7.2,2x,I3,2x,I3)")  Label(1,j),RA(1:4,j), Label(2,j)*AmplScale, Label(3:4,j)
      enddo
      close(unit=29)
      !
      If(LongTrackNr.gt.0) then
         Call ConstructTracks
         !
         Call AmplitudeFitTracks
         !
         If (dt_MTL.gt.0.) then
            Call BinTracks(dt_MTL)  ! t_Mean Track Location
         EndIf
         !
         Call AnalyzeBurst ! Analyze how many sources occur within a certain time-resolution
      EndIf
      !
      !    write(pars,"(1x,A20)") FlashFolder//'/'//trim(datfile)//'_s'
      !    write(pars,"(1x,A20)") trim(datfile)//'_s'
      !write(pars,"(A)") TRIM(DataFolder)//trim(datfile)
      !SystemCommand="call GLE -d pdf -o Map_"//trim(datfile)//".pdf ../Utilities/LeaderTrack.gle "//trim(pars)
      !SystemCommand="call GLE -d pdf -o Imp_"//trim(datfile)//".pdf ../Utilities/SourcesPlot.gle "//trim(pars)
      Call GLEplotControl(PlotType='SourcesPlot', PlotName='Imp_'//TRIM(datfile), &
            PlotDataFile=TRIM(DataFolder)//trim(datfile), Submit=.false.)
      ! call flush(unit_number)
      !      CALL SYSTEM(SystemCommand,stat)
      !write(2,*) SystemCommand
      !write(unit,*) SystemCommand
      !
      !Call GLEplotControl(SpecialCmnd=TRIM(ProgramFolder)//'TrackExe.exe  <'//TRIM(lname))

      If(LongTrackNr.gt.0) then
         Call GLEplotControl(PlotType='TrackScatt', PlotName='TrSc_'//TRIM(datfile), &
            PlotDataFile=TRIM(DataFolder)//trim(datfile), Submit=.false.)
         !SystemCommand="call GLE -d pdf -o TrSc_"//trim(datfile)//".pdf ../Utilities/TrackScatt.gle "//trim(pars)
         !      CALL SYSTEM(SystemCommand,stat)
         !write(2,*) SystemCommand
         !write(unit,*) SystemCommand
      EndIf
      !
998   write(2,*) '======================================='
   enddo !  end of the main reading loop
   !
999 continue
   Call GLEplotControl(Submit=.true.)
   stop
!--------------------------------------------------
  Contains
!===============================================
Subroutine SelSources(Qual,t_start, EndInput)
   ! Apply quality criteria to select sources
   ! On return:   Selected sources in RA
   !              Opened file 29 with header written
   Use Tracks, only : TrackNr, TOrder
   Use Tracks, only : datfile
   Use Tracks, only : RA, maxd, EventNr, Label
   use DataConstants, only : FlashFolder  !, DataFolder, FlashName, Windows
   IMPLICIT none
   Real*8, intent(out) :: Qual,t_start
   Logical, intent(out) :: EndInput
   Character*100, save :: datafile
   Logical :: Old
   Character(len=180) :: lineTXT
   Character*130 :: ZoomBox
   real*8 :: ICVal, Height_1, Chi_max, R_ex, xMin,xMax,yMin,yMax,zMin,zMax,tMin,tMax
   Integer :: i,j, Nxx, DelNEff
   real*8 :: x,y,z,t, val, sigma(1:3), Q
   character*100 :: txt
   integer :: N_EffAnt, Max_EffAnt, Ampl, Wl, Wu
   !
   !
   Chi_max=100.
   Height_1=20.
   DelNEff=40
   Old=.true.
   EndInput=.false.
1  Continue
   Call GetNonZeroLine(lineTXT)
   read(lineTXT,*,iostat=nxx) datafile, ICVal, datfile, xMin,xMax,yMin,yMax,zMin,zMax,tMin,tMax, ZoomBox
   if(nxx.ne.0) then
      Old=.false.
      R_ex=-1.
      read(lineTXT,*,iostat=nxx) datafile, ICVal, Height_1, Chi_max, DelNEff, &
         datfile, xMin,xMax,yMin,yMax,zMin,zMax,tMin,tMax, ZoomBox
      If(nxx.ne.0) then
         read(lineTXT,*,iostat=nxx) datafile, ICVal, Height_1, Chi_max, R_ex, &
            datfile, xMin,xMax,yMin,yMax,zMin,zMax,tMin,tMax, ZoomBox
         If(nxx.eq.0) Then
            Write(2,"(A,F5.2,A,F5.1,A)") 'Quality factor =',Chi_max,' [ns] used with R_ex=',R_ex
            DelNEff=-1
         EndIf
      Else
         Write(2,"(A,F5.2,A,I2,A)") 'RMS =',Chi_max,' [ns] used with Effective antenna cut=',DelNEff,' less than the maximum'
      Endif
      If(nxx.ne.0) Then
         If(datafile(1:1).eq."!") then
            write(2,*) 'presumed commentline ignored on input: ',TRIM(lineTXT)
            goto 1
         Else
            EndInput=.true.
            Return
         EndIf
      EndIf
   Endif
   If(ZoomBox.ne.'NoBox')  ZoomBox=TRIM(DataFolder)//ZoomBox
   Write(2,"(A,F7.2,A,F7.2,A,F7.2)") 'quality cuts: sigma(z)=', ICVal,'[m] @', Height_1,'[km]'
   !Write(2,"(A,I2,A)") 'Effective antenna cut=',DelNEff,' less than the maximum'
   write(2,*) trim(datafile),' ', TRIM(datfile), xMin,xMax,yMin,yMax,zMin,zMax,tMin,tMax
   !
   OPEN(unit=28,FILE=trim(datafile)//'.csv',FORM='FORMATTED',STATUS='OLD',ACTION='READ') ! space separated values
   Do
      Read(28,"(A180)",end=998) lineTXT
      If(lineTXT(1:1).ne.'!') exit
      Write(2,"(A)") TRIM(lineTXT)
   Enddo
   read(lineTXT,*,end=998) t_start,txt     ! from the first line that does not start with !
   write(2,*) 'time offset=',t_start,'[s]'
   !
   j=0
   do
      NXX=0
      If(Old) then
         read(28,*,end=998,iostat=nxx)  i,y,x,z,t,val,sigma
      Else
         read(28,"(A180)",end=998)  lineTXT
         read(lineTXT,*,iostat=nxx)  i,y,x,z,t,val,sigma, N_EffAnt, Max_EffAnt, Ampl, Wl, Wu
         If(nxx.ne.0) then
            !write(2,*) lineTXT
            !flush(unit=2)
            read(lineTXT,*,iostat=nxx)  i,y,x,z,t,val,sigma, N_EffAnt, Max_EffAnt
            !write(2,*) Nxx,lineTXT
            If(nxx.ne.0) exit
            Ampl=100+mod(i,10) ; Wl=1 ; Wu=1
            TOrder=+1  ! older version
            !write(2,*) 'Read error:',i,y,x,z,t,val,sigma, N_EffAnt, Max_EffAnt, Ampl, Wl, Wu
         Endif
      EndIf
      !
      t = (t-t_start)*1000.  ! convert to [ms]
      z=abs(z)
      if(x.le.xMin*1000. .or. x.ge.xMax*1000.) cycle
      if(y.le.yMin*1000. .or. y.ge.yMax*1000.) cycle
      if(z.le.zMin*1000. .or. z.ge.zMax*1000.) cycle
      if(t.le.tMin .or. t.ge.tMax) cycle
      !if(ICVal.gt.val) cycle
      !If(sigma(3).gt. ICVal) cycle
      z=z/1000.d0                   ! convert to [km]
      If(z.gt. Height_1) then
         Q=sigma(3)*Height_1
      Else
         Q=sigma(3)*z
      Endif
      If(Q.gt. ICVal .and. val.gt.ICVal) cycle
      If(.not. old) then
         Q=sqrt(val)
         If(DelNEff .lt. 0) then
            Q=Q + R_ex*(Max_EffAnt-N_EffAnt)/Max_EffAnt
         else
            If( (Max_EffAnt - N_EffAnt) .gt. DelNEff) cycle
         Endif
         If(Q.gt. Chi_max) cycle
      endif
      j=j+1
      RA(1,j)=t ;  RA(2,j)=x/1000.d0 ;  RA(3,j)=y/1000.d0 ;  RA(4,j)=z   ! units [ms], [km], [km], [km]
      Label(1,j)=i ;Label(2,j)=Ampl ; Label(3,j)=Wl ; Label(4,j)=Wu
      if(j.eq.maxd) then
         write(*,*) 'Max. dimension reached of', maxd,' at',i
         write(2,*) 'Max. dimension reached of', maxd,' at',i
         exit
      EndIf
      !Q=sqrt(val) + 35.*(Max_EffAnt-N_EffAnt)/Max_EffAnt
      !write(30,"(1x,i7,4(2x,g14.8),3x,f8.2,i4)") i,sigma,Q,sqrt(val), Max_EffAnt-N_EffAnt
    enddo
    !Close(Unit=30)
    write(2,*) 'Last event read:',i,' @ t=',t,'[s]'
998 continue
    close(unit=28)
    !
    Qual=Chi_max
    EventNr=j
    write(*,"(A,A,A,I6)") '# of sources ',trim(datfile),' =',EventNr
    write(2,"(A,A,A,I6)") '# of sources ',trim(datfile),' =',EventNr
    OPEN(unit=29,FILE=TRIM(DataFolder)//trim(datfile)//'.dat',FORM='FORMATTED',STATUS='unknown',ERR=9)  ! will contain all accepted sources
    write(29,"(6(1x,f8.3),2(1x,f8.3),1x,A,f9.3,A)") xMin, xMax, yMin,yMax,zMin,zMax,tMin,tMax,trim(ZoomBox),t_start*1000.,'  ! '
    Return
9  Continue
   write(2,*) 'file "',TRIM(DataFolder)//trim(datfile)//'.dat','" could not be opened'
   stop 'file open problem'
End Subroutine SelSources
!===============================================
Subroutine Assign2Tracks
! version 17a, use reverse time order
   !  On return:
   ! TrackE(k,i)=j :  Source TrackE(k,i) is the k^th source assigned to track i
   ! TrackENr(i) : Total number of sources assigned to track i
   ! Used:
   ! MaxTrackDist = maximal distance for a source to be removed from the track-head to be assigned to a track
   ! wtr = weigth for new found source position to determine position update of track-head
   Use Tracks, only : PreDefTrackNr, PreDefTrackFile, TOrder
   Use Tracks, only : TrackNr, LongTrackNr, LongTrack_Min, TrackNrMax, TrackLenMax, TrackENr, TrackE, TrackTimeLim, TrackNrLim
   Use Tracks, only : Wtr, MaxTrackDist, HeightFact, NLongTracksMax
   Use Tracks, only : RA, maxd,  EventNr!, Label
   IMPLICIT none
   Real*8 :: TrackPos(1:3,TrackLenMax)  !  TrackPos(1:3,i): position of the head of the i^th track
   Integer :: i_cls(1), n_cls, Cls_track(3)
   Real*8 :: Cls_dist(1:3), dist
   integer :: i,j,k, TotLongTrackEve, j_i, j_f
   !
   !-- Get predefined tracks (if any)
   If(TrackNrLim.gt.TrackNrMax) TrackNrLim=TrackNrMax
   Call ReadPreDefTr(PreDefTrackFile)
    ! Reconstruct TrackPosition at time of source i
    !
   TrackNr=PreDefTrackNr
   TrackENr(:)=0
   !tNrS(:)=1
   Write(2,*) 'PreDefTrackNr=',PreDefTrackNr
   !TrackE(TrackENr(TrackNr),TrackNr)=j
   ! TrackPos(1:3,TrackNr)=RA(2:4,j)
   !TOrder=-1
   !TrackTimeLim=2.
   If(TOrder.eq.1) then
      j_i=1
      j_f=EventNr
   Else
      j_f=1
      j_i=EventNr
      TOrder=-1
   EndIf
   Do j=j_i,j_f,TOrder ! loop over all sources
      Call GetPreDefTrPos(RA(1,j),TrackPos) ! fills positions of all predifined tracks at time RA(1,j)
      n_cls=0  ! number of close-lying tracks
      !write(2,*) j,RA(1,j),TrackPos(1:3,1)
      Do i=1,TrackNr
         dist=0.         ! distance to the next point (in time)
         dist=((TrackPos(3,i)-RA(4,j))/HeightFact)**2  ! weigh vertical distance by factor "HeightFact" less
         !dist=dist/25.      ! used to weigh vertical distance by factor 5 less
         dist=dist +(TrackPos(1,i)-RA(2,j))**2  ! changed to factor 1 less weight for z
         dist=dist +(TrackPos(2,i)-RA(3,j))**2
         Dist=sqrt(dist)       ! distance to the next point (in time)
         !If(j.lt.5) Write(2,*) j,Dist,i,TrackPos(3,i)
         If(dist .gt. MaxTrackDist) cycle
         ! write(*,*) TrackNr,j,i,dist,TrackPos(1:3,i)
         ! Add to existing Track
         n_cls=n_cls+1
         Cls_track(n_cls)=i
         Cls_dist(n_cls)=Dist ! store distances to closest track
         If(n_cls.eq.3) exit
         !goto 10
      enddo
      If(n_cls.eq.0) then  ! no close-lying track, start a new one
         ! Open new track
         !write(2,*) 'n_cls=0',RA(2:4,j),dist
         If(TrackNr.ge. TrackNrLim) cycle
         TrackNr=TrackNr+1           ! found a new track
         TrackENr(TrackNr)=1        ! number of sources on this track
         TrackE(TrackENr(TrackNr),TrackNr)=j    ! source number
         TrackPos(1:3,TrackNr)=RA(2:4,j)
      Else
         !write(2,*) 'n_cls=',n_cls,Cls_track(1:n_cls),Cls_dist(1:n_cls)
         i_cls=minloc(Cls_dist(1:n_cls)) ! get closest track
         ! i_cls=1  !                        get first track that is close
         i=Cls_track(i_cls(1))           ! Track-number of closest track is i
         If(TrackENr(i) .lt. TrackLenMax) then   ! Check if maximal track-length will be exceeded
             TrackENr(i)=TrackENr(i)+1   ! last index for this track
             TrackE(TrackENr(i),i)=j     ! store event-number j in top-position of track i
         endif
         !If(i .gt. PreDefTrackNr) Then  ! do this only for non-predefined tracks
            TrackPos(1:3,i)=(TrackPos(1:3,i)+Wtr*RA(2:4,j))/(1.d0+Wtr) ! update track position
         !EndIf
      Endif
      ! Check for long non-active tracks
      Do i = 1, TrackNr
         k=TrackE(TrackENr(i),i)  ! event number of last member
         If(TOrder*(RA(1,j)-RA(1,k)).lt. TrackTimeLim) cycle ! track should continue; take care of proper time-order
         If(TrackENr(i).lt. LongTrack_Min) then ! delete this track
            If(i.lt.TrackNr) then
               Do k=i+1,TrackNr  ! move tracks
                  TrackENr(k-1)=TrackENr(k)
                  TrackE(1:TrackENr(k-1),k-1)=TrackE(1:TrackENr(k),k)
                  TrackPos(1:3,k-1)=TrackPos(1:3,k)
               Enddo
            EndIf
            TrackNr=TrackNr-1 ! delete last one
            exit
         endif
      Enddo
      !
      !If(RA(1,j).gt.25) Write(2,*) j,RA(1,j),n_cls,i,TrackNr
   enddo ! loop over events
   !
   ! Clean-up short tracks
   LongTrackNr=TrackNr
   TotLongTrackEve=0
   Do i = 1, TrackNr
      If(TrackENr(i).lt. LongTrack_Min) then ! delete this track
         LongTrackNr=LongTrackNr-1
         If(TrackENr(i).le.0) cycle
      endif
      j=TrackE(1,i)  ! source number for first source on this track
      k=TrackE(TrackENr(i),i)  ! Source number of last source on this track
      write(2,"(A,i2,A,I4,A,2F7.1,2I5)") 'TrackNr=',i,', has # pixels=',TrackENr(i),'; t_i, t_f =',RA(1,j),RA(1,k)!,j,k
      TotLongTrackEve=TotLongTrackEve + TrackENr(i)
    enddo
    !
    write(2,"(A,i5,A,I5)") 'Nr of events=',EventNr,' , total # in tracks=',TotLongTrackEve
    write(2,"(A,i2,A,I2)") 'Nr of tracks=',TrackNr,'; # of long tracks=', LongTrackNr
    if(LongTrackNr.gt.NLongTracksMax) then
        LongTrackNr=NLongTracksMax
        write(*,*) 'LongTrackNr reduced to 9 for plotting reasons'
    endif
   !Call Flush(2)
   !
   Return
End Subroutine Assign2Tracks
!===============================================
Subroutine ConstructTracks
   ! Construct mean, smooth, tracks from assigned source locations on return in  LeaderPos
   ! Calculate lateral deviations of source locations from mean track
   ! Writes files  trim(datfile)//i_LongTr//'.dat'
   ! Used:
   ! TimeWin = width of gaussian in time to weigh the sources to obtain mean position of track
   ! Discussion:
   ! The parameter "TimeWin" is essential to smooth out the plotted track. When too small the curve moves back and
   !  forth, following the individual sources, when too big many corners may be cut-off
   Use Tracks, only : TrackNr, LongTrackNr, LongTrack_Min, TrackENr, TrackE, LeaderPos
   Use Tracks, only : Wtr, MaxTrackDist, TimeWin, datfile, TOrder
   Use Tracks, only : RA, maxd, EventNr!, Label
   use DataConstants, only : DataFolder  !, FlashFolder, FlashName, Windows
   IMPLICIT none
   Real*8 :: TrackPos(1:3,TrackLenMax)  !  TrackPos(1:3,i): position of the head of the i^th track
   !
   Integer :: i_LongTr, i,j,k,kk,jk
   Character*8 :: extension
   Real*8 :: Xsq, Ysq, Zsq, tw, w, Pos(1:3), dt
   Real*8 :: RadDev, HorDev, HDist, Velocity, TW2
   !Integer, Parameter :: bin_max=30
   !Integer :: i_bin
   !Real*8 :: Hist1(0:bin_max), Hist2(0:bin_max)
   !Real*8 :: Bin1size=0.0001, Bin2size=0.020
   !
    Write(2,"(A,f7.4,A)") 'Time window for track smoothing=',TimeWin,'[ms]'
   Tw2=TimeWin*TimeWin
    i_LongTr=0
   !Call Flush(2)
   !    Hist1=0.000001 ; Hist2=0.000001
    do i=1,TrackNr
        If(TrackENr(i).lt. LongTrack_Min) cycle
        i_LongTr=i_LongTr+1
        if(i_LongTr.ge.10) exit
        !write(extension,"(A2,i1,A4)") '_s',nxx,'.dat' !,&
        write(extension,"(i1,A4)") i_LongTr,'.dat' !,&
        OPEN(unit=29,FILE=TRIM(DataFolder)//trim(datfile)//trim(extension),FORM='FORMATTED',STATUS='unknown')
        !write(29,*) 'TrackNr=',i
        !MeanTrackLoc(1:3,:)=0.   !, N_MTL, Nmax_MTL=500
        Xsq=0.
        Ysq=0.
        Zsq=0.
        Do k=1,TrackENr(i)
            j=TrackE(k,i) ! j= rank-number of k-th source on track i
            ! Calculate Estimated Leader position
            tw=0.
            Pos(1:3)=0.
            Do kk=1,TrackENr(i) ! average positions of sources along this track with gaussian-in-time weighting
                jk=TrackE(kk,i) ! jk= rank-number of other sources on track i
                dt=(RA(1,j)-RA(1,jk))**2        ! dt has units of t^2 !!!
                If(dt.gt.8.*Tw2) cycle
                w=exp(-dt/Tw2)
                tw=tw+w
                Pos(1:3)=Pos(1:3)+w*RA(2:4,jk)
                ! write(*,*) k,kk,dt,w
            enddo
            LeaderPos(1,k,i)=RA(1,j)         ! Time of the leader head = the time of the last source
            LeaderPos(2:4,k,i)=Pos(1:3)/tw  ! position of the leader head
            !
            RadDev=(RA(4,j)-LeaderPos(4,k,i))*RA(4,j)  ! deviation along the axis from souce to reference antenna
            HorDev=0.d0
            HDist=0.d0
            Do kk=2,3
               RadDev=RadDev+(RA(kk,j)-LeaderPos(kk,k,i))*RA(kk,j)
               HorDev=HorDev+(RA(kk,j)-LeaderPos(kk,k,i))*RA(5-kk,j)*(-1)**kk  ! cross product in horizontal plane
               HDist=HDist+ RA(kk,j)*RA(kk,j)            ! length of vector, distance^2 to reference
            EndDo
            RadDev=RadDev/SQRT(HDist+ RA(4,j)*RA(4,j))
            HorDev=HorDev/SQRT(HDist)
            tw=sqrt(SUM( (LeaderPos(2:4,k,i)-RA(2:4,j))**2 ))  ! place holder for distance from tip
            Velocity=TOrder*SQRT(SUM((LeaderPos(2:4,k,i)-LeaderPos(2:4,k-1,i))**2))/(LeaderPos(1,k,i)-LeaderPos(1,k-1,i))  ! units: 10^6 m/s
            write(29,"(1x,4(2x,g14.8),3x,3(2x,g14.8),3x,g14.8,3(2x,g14.8),3x,g12.3,1x,F8.2)") &
                 LeaderPos(1:4,k,i), RA(2:4,j), LeaderPos(1,k,i)-LeaderPos(1,k-1,i) &
                 ,RadDev,HorDev,RA(4,j)-LeaderPos(4,k,i), Velocity, tw*1000.
            !write(29,"(1x,3i4,4(2x,g14.8),3x,f8.6)")  i,k,j, RA(1:4,j)
            !i_bin=(LeaderPos(1,k,i)-LeaderPos(1,k-1,i))/Bin1size
            !if(i_bin.gt.bin_max) i_bin=bin_max
            !Hist1(i_bin)=Hist1(i_bin) + 1.   ! Histogram for time between subsequent steps in a leader
            !i_bin=(LeaderPos(1,k,i)-LeaderPos(1,k-1,i))/Bin2size
            !if(i_bin.gt.bin_max) i_bin=bin_max
            !Hist2(i_bin)=Hist2(i_bin) + 1.   ! Histogram for time between subsequent steps in a leader
            Xsq=Xsq+(LeaderPos(2,k,i)-RA(2,j))*(LeaderPos(2,k,i)-RA(2,j))
            Ysq=Ysq+(LeaderPos(3,k,i)-RA(3,j))*(LeaderPos(3,k,i)-RA(3,j))
            Zsq=Zsq+(LeaderPos(4,k,i)-RA(4,j))*(LeaderPos(4,k,i)-RA(4,j))
            !If(SQRT(SUM((LeaderPos(2:4,k,i)-RA(2:4,j))*(LeaderPos(2:4,k,i)-RA(2:4,j)))) .gt. 0.06) then
            !   write(2,*) k,j,i,(LeaderPos(2:4,k,i)-RA(2:4,j)), RadDev, HorDev
            !Endif
        enddo  ! k=1,TrackENr(i) loop over sources in track i
        Write(2,"(A,i2,I5,A,3F7.1)") 'Track#',i,TrackENr(i),', RMS(E,N,H)[m]:', &
         sqrt(Xsq/TrackENr(i))*1000., sqrt(Ysq/TrackENr(i))*1000., sqrt(Zsq/TrackENr(i))*1000.
        close(unit=29)
    enddo ! i=1,TrackNr
    !
   Return
End Subroutine ConstructTracks
!===============================================
Subroutine AmplitudeFitTracks
   ! Make Source amplitude histogram and fit this with exp or power law
   !      write(29,"(1x,i8,4(2x,g14.8),3x,F7.2,2x,I3,2x,I3)")  Label(1,j),RA(1:4,j), Label(2,j)/100., Label(3:4,j)
   ! Writes files  trim(datfile)//??//'.dat'
   Use Tracks, only : TrackNr,  TrackENr !, TrackE,  datfile
   Use Tracks, only :  Label
   Use Tracks, only : AmplScale, d_AmplScale, AmplitudeFit
   Use Tracks, only : Selection_sort
   IMPLICIT none
   Real*8 :: TrackAmp(1:1,TrackLenMax)  !  Amplitudes of sources in the i^th track
   Integer :: TrackSourceAmp(1:1,TrackLenMax)  !  rank number of sources on the i^th track
   !
   Real*8 :: Ampl, Max_Ampl   !, dt_MTL=0.2  ! [ms]
   Integer :: i,k, i_Ampl, Count
   Integer ( kind = 4 ) :: meqn  ! Number of data points to fit
   Character*8 :: extension
   Real*8 :: Ampl10, AmplThr, X(2),ChiSQDF, CalcN
   !
   !   write(2,*) 'long track',i_LongTr,' has ',t_MTL,'bins of ',dt_MTL,'[ms]'
   do i=1,TrackNr
      !write(extension,"(A2,i1,A4)") '_s',nxx,'.dat' !,&
      write(extension,"(i1,A4)") i,'.dat' !,&
      !OPEN(unit=27,FILE='files/'//trim(datfile)//'_??_'//trim(extension),FORM='FORMATTED',STATUS='unknown')
      !write(27,"('! nr', 4(1x,A14),3x,3(1x,A12),2x,A10)") 't [ms]','E [km]','N','h','v_E','v_N','v_h','v [10^6m/s]'
      !
      Call AmplitudeFit(TrackENr(i), TrackE(1,i))
      !write(2,*) MeanTrackLoc(1,1:t_MTL)
      !Call Flush(2)
      !close(unit=27)
   enddo
    !
   Return
End Subroutine AmplitudeFitTracks
!===============================================
Subroutine BinTracks(dt_MTL)
   ! Construct mean tracks from assigned source locations on return in  LeaderPos
   ! Calculate lateral deviations of source locations from mean track
   ! Writes files  trim(datfile)//i_LongTr//'.dat'
   Use Tracks, only : TrackNr, LongTrackNr, LongTrack_Min, TrackENr, TrackE, LeaderPos, datfile
   Use Tracks, only : Wtr, MaxTrackDist, TimeWin, datfile, TOrder
   Use Tracks, only : RA, maxd,  EventNr!, Label
   use DataConstants, only : DataFolder  !, FlashFolder, FlashName, Windows
   IMPLICIT none
   Real*8, intent(in) :: dt_MTL
   Real*8 :: TrackPos(1:3,TrackLenMax)  !  TrackPos(1:3,i): position of the head of the i^th track
   !
   Integer, parameter :: Nmax_MTL=500
   Real*8 :: MeanTrackLoc(1:4,Nmax_MTL)   !, dt_MTL=0.2  ! [ms]
   Integer :: N_MTL, i_MTL, t_MTL, t_old, i,k, i_LongTr
   Character*8 :: extension
   Real*8 :: dt,dist,t0
   !
   i_LongTr=0
   !   write(2,*) 'long track',i_LongTr,' has ',t_MTL,'bins of ',dt_MTL,'[ms]'
   do i=1,TrackNr
      If(TrackENr(i).lt. LongTrack_Min) cycle
      i_LongTr=i_LongTr+1
      if(i_LongTr.ge.10) exit
      !write(extension,"(A2,i1,A4)") '_s',nxx,'.dat' !,&
      write(extension,"(i1,A4)") i_LongTr,'.dat' !,&
      OPEN(unit=27,FILE=TRIM(DataFolder)//trim(datfile)//'_bin_'//trim(extension),FORM='FORMATTED',STATUS='unknown')
      write(27,"('! nr', 4(1x,A14),3x,3(1x,A12),2x,A10)") 't [ms]','E [km]','N','h','v_E','v_N','v_h','v [10^6m/s]'
      t_old=-99
      i_MTL=0  ! just counter for filled bins
      MeanTrackLoc(1:4,:)=0.
      !write(2,*) 'long track',i_LongTr,' has ',t_MTL,'bins of ',dt_MTL,'[ms]'
      !Call Flush(2)
      t0=RA(1,TrackE(1,i)) ! start time of this track
      Do k=1,TrackENr(i)
         j=TrackE(k,i) ! j= rank-number of k-th source on track i
         t_MTL=TOrder*(RA(1,j)-t0)/dt_MTL
         !write(2,*) i_old,t_MTL,i,j,MeanTrackLoc(1,1:t_MTL)
         !Call Flush(2)
         If(t_MTL.gt.t_old) then
            i_MTL=i_MTL+1
            !If(i_old.gt.0) MeanTrackLoc(1:4,i_old)=MeanTrackLoc(1:4,i_old)/N_MTL
            if(i_MTL.gt.Nmax_MTL) Then
               i_MTL=i_MTL-1
               exit
            Endif
            MeanTrackLoc(1:4,i_MTL)=RA(1:4,j)
            t_old=t_MTL
            N_MTL=1
         Else
            MeanTrackLoc(1:4,i_MTL)=(N_MTL*MeanTrackLoc(1:4,i_MTL)+RA(1:4,j))/(N_MTL+1)  ! running average
            N_MTL=N_MTL+1
         Endif
      Enddo
      write(2,*) 'long track',i_LongTr,' has ',i_MTL,'bins of ',dt_MTL,'[ms] filled'
      !write(2,*) MeanTrackLoc(1,1:t_MTL)
      !Call Flush(2)
      !i_old=0
      Do k=2,i_MTL
         !If(MeanTrackLoc(1,k).gt. 0.) then
            !If(t_old.gt.0) then
               dt=TOrder*(MeanTrackLoc(1,k)-MeanTrackLoc(1,k-1))
               dist=sqrt(SUM((MeanTrackLoc(2:4,k)-MeanTrackLoc(2:4,k-1))**2))
               write(27,"(I4,1x,4(1x,g14.8),3x,3(1x,g12.4),2x,g10.4,2x,g10.4)") k, &
                  MeanTrackLoc(1:4,k),(MeanTrackLoc(2:4,k)-MeanTrackLoc(2:4,k-1))/dt, dist/dt, &
                  (MeanTrackLoc(1,k)+MeanTrackLoc(1,k-1))/2.
            !Endif
            !t_old=k
         !Endif
      Enddo
      close(unit=27)
   enddo
    !
   Return
End Subroutine BinTracks
!===============================================
Subroutine AnalyzeBurst()
    ! Make a histogram of the source times in a track where each pulse is smeared with a gaussian
    ! with half width of t_resol, where this resolution changes. The histogram is stored in TBurst_trace(T)
    ! Then check how often >3.5 (N1) or >2.5 (N2) sources are within the t_resol
    Use Tracks, only : TOrder, datfile
    Use constants, only : dp
   use DataConstants, only : DataFolder  !, FlashFolder, FlashName, Windows
    use FFT, only : RFTransform_su, DAssignFFT, RFTransform_CF
    IMPLICIT none
    integer, parameter :: TB_max=131072, N_nu=TB_max/2  ! 65536 ! 65536=2^16 ; 32768=2^15 ! 2048=2^11 ! 262144=2^18 ! 131072=2^17
    real*8 :: t_resol,T_max, T_min, TimeDur, A, B
    Real*8 :: TBurst_trace(TB_max), d_nu, F_Th, nu, BinSize
    Complex(dp) :: Cnu(0:N_nu)
    integer :: i,k,N_k,NTSampl, T_C, T, T1, T2, N1, N2, NM2, NM3, tr
    Character*8 :: extension
    !
    i=1
    write(2,*) TOrder,RA(1,TrackE(1,i)),TrackE(1,i),RA(1,TrackE(TrackENr(i),i)),TrackE(TrackENr(i),i)
    Do i=1,1   ! analyze only these tracks
      If(TrackENr(i).lt. LongTrack_Min) cycle ! check if this track qualifies as a long one
      If(TOrder.gt.0) then
         T_min=RA(1,TrackE(1,i))                 ! time of first event on track i
         T_max=RA(1,TrackE(TrackENr(i),i))       ! time of last event on track i
      Else
         T_max=RA(1,TrackE(1,i))                 ! time of first event on track i
         T_min=RA(1,TrackE(TrackENr(i),i))       ! time of last event on track i
      EndIf
      TimeDur= (T_max - T_min)
      t_resol=1  ! mili second
      Do tr=1,9
         t_resol=t_resol/2
         NTSampl=int(TimeDur*2/t_resol)
         if(NTSampl .gt. TB_max) exit
         TBurst_trace=0.
         Do k=1,TrackENr(i) ! loop over all events in track i
             j=TrackE(k,i)   ! get event number of k-th track member
             T_C=Int((RA(1,j)-T_min)*2/t_resol)
             T1=T_C-8
             if(T1.lt.1) T1=1
             T2=T_C+8
             if(T2.gt.NTSampl) T2=NTSampl
             Do T=T1,T2
                 TBurst_trace(T)=TBurst_trace(T)+ exp(-(T/2.-(RA(1,j)-T_min)/t_resol)**2)
             Enddo
         enddo
         A=MaxVal(TBurst_trace(:))
         ! B=MinVal(TBurst_trace(:))
         N1=COUNT(TBurst_trace .gt. 3.5)
         N2=COUNT(TBurst_trace .gt. 2.5)
         NM2=COUNT(TBurst_trace .gt. A/2.)
         NM3=COUNT(TBurst_trace .gt. A/3.)
         Write(2,"(A,i2,A,f7.4,A,F7.1,A,i5,A,i5,A,I6,3I6)") 'track #=',i,', sample=',t_resol,'[ms], maximum=',A, &
             ', > 3.5:',N1,', > 2.5:',N2,', #samples=',NTSampl,NM2,NM3
         If(i.eq.1) then
            write(extension,"(i2.2,A4)") tr,'.dat' !,& SourceTimeDistribution (STD)
            OPEN(unit=27,FILE=TRIM(DataFolder)//trim(datfile)//'_STD-'//trim(extension),FORM='FORMATTED',STATUS='unknown')
            write(27,"('! ', 4(1x,A14),3x,3(1x,A12),2x,A10)") 't [ms]','Amplitude','N','h','v_E','v_N','v_h','v [10^6m/s]'
            Do T=1,NTSampl
               write(27,*) T_min+T*t_resol/2.,TBurst_trace(T)
            enddo
            Close(unit=27)
         endif
      enddo  ! tr
      !
      Call RFTransform_su(TB_max)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !
      Call RFTransform_CF(TBurst_trace,Cnu)
      OPEN(unit=27,FILE=TRIM(DataFolder)//trim(datfile)//'_STD-nu.dat',FORM='FORMATTED',STATUS='unknown')
      write(27,"('! ', 4(1x,A14),3x,3(1x,A12),2x,A10)") 'nu [kHz]','Amplitude','N','h','v_E','v_N','v_h','v [10^6m/s]'
      d_nu=1./(t_resol*TB_max/2)
      F_Th=d_nu/2
      N_k=0
      A=0.
      B=0.
      BinSize=d_nu/128
      BinSize=BinSize/128
      Do T=0,TB_max/2
         nu=T*d_nu
         If(nu.gt.F_th) then
            write(27,*) B/N_k,A/N_k, nu, BinSize, N_k
            F_th=1.05*nu
            N_k=1
            A=ABS(Cnu(T))
            B=nu
         Else
            A=A+ABS(Cnu(T))
            B=B+nu
            N_k=N_k+1
         Endif
         !write(27,*) T*d_nu,ABS(Cnu(T))
      enddo
      Close(unit=27)
      Call DAssignFFT()         !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !
    enddo ! i, track number
end Subroutine AnalyzeBurst
END PROGRAM rewriteEvent
!===================================================================
Subroutine ReadPreDefTr(PreDefTrackFile)
   Use Tracks, only : PreDefTrack, tNrS_max, PreDefTrackNr, PreDefTrackNr_Max, NrPreDefTrackPoints_Max, tNrS
   IMPLICIT none
   Character(len=100), Intent(IN) :: PreDefTrackFile
   Integer :: nxx, i_track, i_t
   Character(len=5) :: Extension
   !
   PreDefTrackNr=0
   tNrS(1:PreDefTrackNr_Max)=1  ! needed in GetPreDefTrPos
   Do i_track=1,PreDefTrackNr_Max
      write(Extension,"(i1,A4)") i_track,'.dat' !,&
      OPEN(UNIT=22,STATUS='old',ACTION='Read',FILE=trim(PreDefTrackFile)//trim(extension),IOSTAT=nxx)
      If(nxx .ne. 0) then
         write(2,*) 'No predefined track file:',trim(PreDefTrackFile)//trim(extension)
         exit
      Endif
      !write(2,*) 'ReadPreDefTr:',trim(PreDefTrackFile)//trim(extension)
      Do i_t=1,NrPreDefTrackPoints_Max
         read(22,*,iostat=nxx) PreDefTrack(0:3, i_t, i_track)
         !write(2,*) PreDefTrack(0:3, i_t, i_track), i_t, i_track
         If(nxx.ne.0) exit
      Enddo
      Close(UNIT=22)
      write(2,*) 'ReadPreDefTr:',trim(PreDefTrackFile)//trim(extension), i_t-1
      tNrS_max(i_track)=i_t-1
      PreDefTrackNr=i_track
   Enddo
   !
   Return
End Subroutine ReadPreDefTr

!=========================================
Subroutine GetPreDefTrPos(t_source,TrackPos)
   Use Tracks, only : TOrder, PreDefTrack, tNrS_max, PreDefTrackNr, PreDefTrackNr_Max,tNrS
   IMPLICIT none
   Real*8,intent(IN) :: t_source
   Real*8,intent(OUT) :: TrackPos(1:3,*)
   !
   Integer :: tNr, i_Pre
   Real*8 ::DtP, DtM
   !
   !If(TOrder) less than 0, use reverse t-order
   If(PreDefTrackNr .lt. 1 ) Return
   !write(2,*) 't_source',t_source,PreDefTrackNr
   Do i_Pre=1,PreDefTrackNr
      !write(2,*) 'tNrS_max(i_Pre), tNrS(i_Pre)',i_pre,tNrS_max(i_Pre), tNrS(i_Pre)
   1  Continue
      tNr=tNrS(i_Pre)  ! Index of last used point on this track
      If(tNr .gt. tNrS_max(i_Pre)) cycle ! Do not update position of this track
      DtM= TOrder*(t_source - PreDefTrack(0,tNr,i_Pre))
      !write(2,*) 'DtM:',DtM,tNr,t_source, PreDefTrack(0,tNr,i_Pre)
      If(DtM .lt. 0) then ! track not yet at this position
         cycle
      Endif
      DtP=TOrder*(PreDefTrack(0,tNr+1,i_Pre)- t_source)
      !write(2,*) 'tNr',tNr,DtM,DtP
      If(DtP .lt. 0 ) then
         tNrS(i_Pre)=tNrS(i_Pre)+1
         goto 1
      Endif
      TrackPos(1:3,i_Pre)= &
      (PreDefTrack(1:3,tNr,i_Pre)*DtP + PreDefTrack(1:3,tNr+1,i_Pre)*DtM)/(DtM+DtP) ! average between earlier and later positions
   EndDo
   Return
End Subroutine GetPreDefTrPos
! ----------------------------------------------------------------------------------------
!     shellin = 'gle /d pdf '//trim(dummy3(i))//'.gle'
!     CALL system(shellin)
!     shellin = 'epstopdf '//trim(dummy3(i))//'.eps'
!     call system(shellin)
