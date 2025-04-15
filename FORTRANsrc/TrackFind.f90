Module TrackConstruct
   Use constants, only : dp
   Character(len=100) :: PreDefTrackFile
   Integer, parameter :: PreDefTrackNr_Max=10  !  Maximum number of Pre-defined Tracks
   Integer, parameter :: NrPreDefTrackPoints_Max=900   ! Maximum number of defining times per Pre-defined Track
   Real(dp), save :: PreDefTrack(0:3,1:NrPreDefTrackPoints_Max,1:PreDefTrackNr_Max)
   Integer, save :: tNrS_max(1:NrPreDefTrackPoints_Max), PreDefTrackNr
   Integer, save :: tNrS(1:PreDefTrackNr_Max)=1
   Integer :: TrackNr, LongTrackNr, LongTrack_Min=20 !12
   Integer, parameter :: TrackNrMax=19, TrackLenMax=5000
   Integer :: TrackENr(TrackLenMax), TrackE(TrackLenMax,TrackNrMax)
   Real*8 :: LeaderPos(1:4,0:TrackLenMax,TrackNrMax)
   Real*8 :: Wtr  ! weighting of new source location for track definition
   Real*8 :: Aweight ! Determines importance of amplitude in weighting new point for track (=0 is not important)
   Real*8 :: MaxTrackDist ! maximal distance from trackhead
   Real*8 :: TimeWin   !  Time window for track smoothing
   Real*8 :: HeightFact   !  Factor to increase "MaxTrackDist" for height, used for track smoothing
   Real*8 :: TrackTimeLim=5. ! maximal distance between sources to mark short tracks for deletion
   Real*8 :: SrcDensTimeResol
   Integer :: TrackNrLim=1
   Integer :: TOrder=-1  ! set to -1 for running tracks backward
   Integer :: NLongTracksMax
   Real*8 :: MeanTrackLoc(1:3,500), dt_MTL
   Integer :: Nmax_MTL=500
contains
!End Module TrackData
!===============================================
!  Contains
!===============================================
Subroutine Assign2Tracks(RA, SrcI20_r, SourcTotNr)
! version 17a, use reverse time order
   !  On return:
   ! TrackE(k,i)=j :  Source TrackE(k,i) is the k^th source assigned to track i
   ! TrackENr(i) : Total number of sources assigned to track i
   ! Used:
   ! MaxTrackDist = maximal distance for a source to be removed from the track-head to be assigned to a track
   ! wtr = weigth for new found source position to determine position update of track-head
    Use constants, only : dp
   !Use PredefinedTracks, only : ReadPreDefTr
   !Use TrackConstruct, only : TrackNr, LongTrackNr, LongTrack_Min, TrackNrMax, TrackLenMax
   !Use TrackConstruct, only : Wtr, MaxTrackDist, HeightFact, NLongTracksMax, TrackENr, TrackE, TrackTimeLim, TrackNrLim
   IMPLICIT none
   Real(dp), Intent(in) :: RA(4,*)
   Real, intent(in) :: SrcI20_r(*)
   Integer, intent(in) :: SourcTotNr ! , Label(4,*)
   Real*8 :: TrackPos(1:3,TrackNrLim)  !  TrackPos(1:3,i): position of the head of the i^th track
   Real*8 :: TrackWeight(TrackNrLim), AmplWeight
   Integer :: i_cls(1), n_cls, Cls_track(3)
   Real*8 :: Cls_dist(1:3), dist, PWtr
   integer :: i,j,k, TotLongTrackEve, j_i, j_f
   Character(len=120) :: PlotFile
   !
   !-- Get predefined tracks (if any)
   If(TrackNrLim.gt.TrackNrMax) TrackNrLim=TrackNrMax
   Call ReadPreDefTrFile(PreDefTrackFile)
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
      j_f=SourcTotNr
   Else
      j_f=1
      j_i=SourcTotNr
      TOrder=-1
   EndIf
   !write(2,*) 'Assign2Tracks: SourcTotNr=',SourcTotNr, j_i,j_f
   Do j=j_i,j_f,TOrder ! loop over all sources
      Call GetPreDefTrPos(RA(1,j),TrackPos) ! fills positions of all predifined tracks at time RA(1,j)
      n_cls=0  ! number of close-lying tracks
      !write(2,*) j,RA(1,j),TrackPos(1:3,1)
      Do i=1,TrackNr
         !dist=0.         ! distance to the next point (in time)
         dist=((TrackPos(3,i)-RA(4,j))/HeightFact)**2  ! weigh vertical distance by factor "HeightFact" less
         !dist=dist/25.      ! used to weigh vertical distance by factor 5 less
         dist=dist +(TrackPos(1,i)-RA(2,j))**2  ! changed to factor 1 less weight for z
         dist=dist +(TrackPos(2,i)-RA(3,j))**2
         Dist=sqrt(dist)       ! distance of track-head(i) to the next (in time) source(j)
         !If(i.lt.2) Write(2,*) j,RA(1,j),Dist,i,TrackPos(3,i)
         If(dist .gt. MaxTrackDist) cycle
         ! write(2,*) TrackNr,j,i,dist,TrackPos(1:3,i)
         ! Add to existing Track
         n_cls=n_cls+1
         Cls_track(n_cls)=i
         Cls_dist(n_cls)=Dist ! store distances to closest track
         If(n_cls.eq.3) exit
         !goto 10
      enddo
      AmplWeight=(SrcI20_r(j)*Aweight+1.)  ! Amplitude determined weight of the new source
      !write(2,*) 'Assign2Tracks: Source=',j,n_cls,dist,TrackNr, TrackNrLim
      !
      If(n_cls.eq.0) then  ! no close-lying track, start a new one
         ! Open new track
         !write(2,*) 'n_cls=0',RA(2:4,j),dist
         If(TrackNr.ge. TrackNrLim) cycle
         TrackNr=TrackNr+1           ! found a new track
         TrackENr(TrackNr)=1        ! number of sources on this track
         TrackE(TrackENr(TrackNr),TrackNr)=j    ! source number
         TrackPos(1:3,TrackNr)=RA(2:4,j)
         TrackWeight(TrackNr)=AmplWeight  !  Intensity
         !write(2,*) 'Assign2Tracks: source,newtrack=',j,TrackNr
      Else
         !write(2,*) 'n_cls=',n_cls,Cls_track(1:n_cls),Cls_dist(1:n_cls)
         i_cls=minloc(Cls_dist(1:n_cls)) ! get closest track
         ! i_cls=1  !                        get first track that is close
         i=Cls_track(i_cls(1))           ! Track-number of closest track is i
         If(TrackENr(i) .lt. TrackLenMax) then   ! Check if maximal track-length will be exceeded
             TrackENr(i)=TrackENr(i)+1   ! last index for this track
             TrackE(TrackENr(i),i)=j     ! store event-number j in top-position of track i
         endif
         PWtr=Wtr*AmplWeight/TrackWeight(i)
         TrackWeight(i)=(TrackWeight(i)+PWtr*AmplWeight)/(1.d0+PWtr)  !  update Intensity trackhead
         TrackPos(1:3,i)=(TrackPos(1:3,i)+PWtr*RA(2:4,j))/(1.d0+PWtr) ! update track position
         !
         ! or ?
 !        TrackPos(1:3,i)=(TrackPos(1:3,i)+Wtr*RA(2:4,j))/(1.d0+Wtr) ! update track position
         !write(2,*) 'Assign2Tracks: source,track=',j,i, sqrt(sum((TrackPos(1:3,i)-RA(2:4,j))**2 ))
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
    write(2,"(A,i5,A,I5)") 'Nr of events=',SourcTotNr,' , total # in tracks=',TotLongTrackEve
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
Subroutine ConstructTracks(DS_Mode, RA, SrcI20_r, SrcPZen, SrcPAzi, SrcWidth, IPerm, SourcTotNr, PlotFile, PolarAna, Stk_NEh)
   ! Construct mean, smooth, tracks from assigned source locations on return in  LeaderPos
   ! Calculate lateral deviations of source locations from mean track
   ! Used:
   ! TimeWin = width of gaussian in time to weigh the sources to obtain mean position of track
   ! Discussion:
   ! The parameter "TimeWin" is essential to smooth out the plotted track. When too small the curve moves back and
   !  forth, following the individual sources, when too big many corners may be cut-off
   !Use TrackConstruct, only : TrackNr, LongTrackNr, LongTrack_Min, TrackENr, TrackE, LeaderPos
   !Use TrackConstruct, only : Wtr, MaxTrackDist, TimeWin, TOrder, TrackLenMax
   Use constants, only : dp, pi
   IMPLICIT none
   Real(dp), Intent(in) :: RA(4,*)
   Real, intent(in) :: SrcI20_r(*), SrcPZen(*), SrcPAzi(*)
   Integer, intent(in) :: DS_Mode, SourcTotNr, SrcWidth(*), IPerm(*) ! , Label(4,*)
   Logical, intent(in) :: PolarAna
   Complex, intent(in) :: Stk_NEh(1:6,*)
   Character(len=*) :: PlotFile
   !Real*8 :: TrackPos(1:3,TrackLenMax)  !  TrackPos(1:3,i): position of the head of the i^th track
   !
   Integer :: i_LongTr, i,j,k,kk,jk
   Character*8 :: extension
   Real*8 :: Xsq, Ysq, Zsq, tw, w, Pos(1:3), dt
   Real*8 :: RadDev, HorDev, HDist, Velocity, TW2, V_hor, V_th, V_ph
   complex  :: LeaderStks(1:6), AveStks(1:3,1:3) !  ,TrackLenMax)
   Real(dp) :: PolZen(1:3), PolAzi(1:3), PolMag(1:3), PoldOm(1:3), Thet1, Phi1
   Real(dp) :: Polv(1:3), Trv(1:3), distv(1:3), vv(1:3), PdothatV, PprphatV, PdothatTr, PdotTr, PdotD, dist, Trd,PdotR
   Real(dp) :: Vh(1:3), Qv(1:3), Qd, PdotQ, Pdoth, wp,twp
   logical :: prin
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
   prin=.false.   !  Do not print polarization observables
   !prin=.true.   !  Do not print polarization observables
   PolZen(1:3)=0. ; PolAzi(1:3)=0. ; PolMag(1:3)=0. ; PoldOm(1:3)=0.
   do i=1,TrackNr
      !write(2,*) 'ConstructTracks',i,TrackENr(i), LongTrack_Min, PolarAna
      If(TrackENr(i).lt. LongTrack_Min) cycle
      i_LongTr=i_LongTr+1
      if(i_LongTr.ge.10) exit
      !write(extension,"(A2,i1,A4)") '_s',nxx,'.dat' !,&
      write(extension,"(i1,A4)") i_LongTr,'.plt' !,&
      OPEN(unit=29,FILE=TRIM(PlotFile)//trim(extension),FORM='FORMATTED',STATUS='unknown')
      write(29,"('!',T22,'Running Average Leader Head position',T73,'Source position',T121,'Derived quantities')")
      write(29,"('!',T7,'time',T22,'North=y',T38,'East=x',T53,'height=z',T73,'North=y',T89,'East=x',T104,'z=height', &
         T121,'Delta_t[ms]',T137,'Delta_Rad[km]',T153,'Delta_Rxz',T169,'Delta_h',T188,'v[km/ms]',T202,'Dist.toTip[m]')")
      If(PolarAna) Then
         OPEN(unit=30,FILE=TRIM(PlotFile)//'Angls'//trim(extension),FORM='FORMATTED',STATUS='unknown')
         write(30,"('!',T10,'Running Average Leader Head position',T47,'Velocity angle',T77,&
            'Polarization amplitudes & angles')")
         write(30,"('!',T5,'time',T17,'North=y',T28,'East=x',T39,'z=height',T50,'Th(Zen),Azim(N)',T72,'Ampl', &
            T82,'Th(Zen), Azim(N), d(omga)',T116,'Repeated twice more', T185,'F(P.V), F(P.R), Dist, P.Dist,   TrDist,   P.Tr')")
      EndIf
      !MeanTrackLoc(1:3,:)=0.   !, N_MTL, Nmax_MTL=500
      Xsq=0.
      Ysq=0.
      Zsq=0.
      vv(1:3)=0.
      Do k=1,TrackENr(i)
            j=TrackE(k,i) ! j= rank-number of k-th source on track i
            ! Calculate Estimated Leader position
            tw=0.
            twp=0.
            Pos(1:3)=0.
            LeaderStks(1:6)=0.
            Do kk=1,TrackENr(i) ! average positions of sources along this track with gaussian-in-time weighting
                jk=TrackE(kk,i) ! jk= rank-number of other sources on track i
                dt=(RA(1,j)-RA(1,jk))**2        ! dt has units of t^2 !!!
                If(dt.gt.8.*Tw2) cycle
                w=exp(-dt/Tw2)*(SrcI20_r(j)*Aweight+1.)  ! Include intensity in weight factor
                tw=tw+w
                Pos(1:3)=Pos(1:3)+w*RA(2:4,jk)
                If(PolarAna) Then
                  If(DS_Mode .eq. 3) Then
                     wp=w*SrcWidth(Iperm(jk))
                  Else
                     wp=w
                  EndIf
                  twp=twp+wp
                  LeaderStks(1:6)=LeaderStks(1:6)+wp*Stk_NEh(1:6,Iperm(jk))
              EndIf
            enddo
            LeaderPos(1,k,i)=RA(1,j)         ! Time of the leader head = the time of the last source
            LeaderPos(2:4,k,i)=Pos(1:3)/tw  ! position of the leader head
            If(PolarAna) Then
               LeaderStks(1:6)=LeaderStks(1:6)/twp
               Call PolPCACathCon(LeaderStks, PolZen, PolAzi, PolMag, PoldOm, prin)
            EndIf
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
            dist=1000.*sqrt(SUM( (LeaderPos(2:4,k,i)-RA(2:4,j))**2 ))  ! distance from tip; in [m]
            If(k.ge.2) then  ! otherwise LeaderPos(1,k-1,i) is not defined
               If(abs(LeaderPos(1,k,i)-LeaderPos(1,k-1,i)).lt. 1.D-7) Then
                  If((k.le.3) ) Then
                     Velocity=0.  ! units: 10^6 m/s
                  Else
                     If((LeaderPos(1,k,i)-LeaderPos(1,k-2,i)).ne. 0.) Then
                        Velocity=TOrder*SQRT(SUM((LeaderPos(2:4,k,i)-LeaderPos(2:4,k-2,i))**2))/ &
                              (LeaderPos(1,k,i)-LeaderPos(1,k-2,i))  ! units: 10^6 m/s
                     EndIf
                     !write(2,*) 'Mean Velocity',k, Velocity, SQRT(SUM((LeaderPos(2:4,k,i)-LeaderPos(2:4,k-2,i))**2)), &
                     !   (LeaderPos(1,k,i)-LeaderPos(1,k-2,i))
                  EndIf
               Else
                  Velocity=TOrder*SQRT(SUM((LeaderPos(2:4,k,i)-LeaderPos(2:4,k-1,i))**2))/(LeaderPos(1,k,i)-LeaderPos(1,k-1,i))  ! units: 10^6 m/s
                  Vv(1:3)=(LeaderPos(2:4,k,i)-LeaderPos(2:4,k-1,i))/(LeaderPos(1,k,i)-LeaderPos(1,k-1,i))
                  V_hor=sqrt(Vv(1)*Vv(1) + Vv(2)*Vv(2))
                  V_th=atan2(V_hor,Vv(3))*180./pi  ; V_ph=atan2(Vv(2),Vv(1))*180./pi
                  !write(2,*) 'Velocity',k, Velocity, SQRT(SUM((LeaderPos(2:4,k,i)-LeaderPos(2:4,k-1,i))**2)), &
                  !   (LeaderPos(1,k,i)-LeaderPos(1,k-1,i))
               EndIf
               If(PolarAna .and. (DS_Mode .eq. 3)) Then  ! ATRID
                  Vh(1:3)=Vv(1:3)/Velocity     ! normalized to unity
                  polv(1)=sin(SrcPZen(Iperm(j))*pi/180.) * cos(SrcPAzi(Iperm(j))*pi/180.)
                  Polv(2)=sin(SrcPZen(Iperm(j))*pi/180.) * sin(SrcPAzi(Iperm(j))*pi/180.)
                  Polv(3)=cos(SrcPZen(Iperm(j))*pi/180.)
                  PdothatV=abs(sum(Polv(1:3)*Vh(1:3)) )
                  distv(1:3)=1000.*(LeaderPos(2:4,k,i)-RA(2:4,j))  ! in [m]
                  PdotD=abs(sum(Polv(1:3)*distv(1:3))/(dist+0.0001))    ! Add 0.1 mm to avoid /0.
                  Trv(1:3)=distv(1:3)-Vh(1:3) * sum(distv(1:3)*Vh(1:3))   ! transverse distance vector from leader; in [m]
                  Trd=sqrt(sum(Trv(1:3)*Trv(1:3)) ) !+ 0.000001)   ! Add 0.001 mm^2 to avoid /0
                  Trd=TrD+ 0.000001   ! Add 0.001 mm to avoid /0
                  Qv(1) = Trv(2)*Vh(3)-Trv(3)*Vh(2)  ! Q= (Tr x V)
                  Qv(2) = Trv(3)*Vh(1)-Trv(1)*Vh(3)
                  Qv(3) = Trv(1)*Vh(2)-Trv(2)*Vh(1)
                  Qd=sqrt(sum(Qv(1:3)*Qv(1:3)))
                  Qd=Qd+ 0.00001   ! Add 0.001 mm to avoid /0; should equal Trd
                  PdotTr=abs(sum(Polv(1:3)*Trv(1:3))/Trd )
                  PdotQ=abs(sum(Polv(1:3)*Qv(1:3))/Qd )
                  PdotR=abs(sum(Polv(1:3)*RA(2:4,j))/sqrt(sum(RA(2:4,j)*RA(2:4,j))) )
                  Pdoth=abs(Polv(3))
                  If(TrD.lt.0.2) Then
                     PdotQ=0.
                     PdotTr=0.
                  EndIf
               EndIf
               If(k.eq.2) Then
                  Thet1=V_th ;  Phi1=V_ph
               EndIf
               If( (abs(PolAzi(1)-Phi1) .gt. 90.) .and. (abs(PolAzi(1)-Phi1) .lt. 270.) ) Then
                  PolAzi(1)=PolAzi(1)+180.
                  If(PolAzi(1) .gt. 180.) PolAzi(1)=PolAzi(1)-360.
                  PolZen(1)=180.-PolZen(1)
               EndIf
               Phi1=PolAzi(1)
               write(29,"(1x,4(2x,g14.8),3x,3(2x,g14.8),3x,g14.8,3(2x,g14.8),3x,g12.3,1x,F8.2)") &
                    LeaderPos(1:4,k,i), RA(2:4,j), LeaderPos(1,k,i)-LeaderPos(1,k-1,i) &
                    ,RadDev,HorDev,RA(4,j)-LeaderPos(4,k,i), Velocity, dist
               If(PolarAna) write(30,"(1x,F11.6,3(F11.5), 2x, 2(',',f7.2) ,2x, 3(' , ',g12.4, 3(',',f7.2)), &
                     2x,2F6.3,4(F10.3,F9.4) )") &
                    LeaderPos(1:4,k,i), V_th, V_ph, (PolMag(j), PolZen(j), PolAzi(j), PoldOm(j), j=1,3), &
                    PdothatV, PdotR, dist, PdotD, Trd, PdotTr, Qd, PdotQ, Pdoth
               !If(PolarAna) write(2,"(1x,F11.6,3(F11.5), 2x, 2(',',f7.2) ,2x, 3(' , ',g12.4, 3(',',f7.2)) )") &
               !     LeaderPos(1:4,k,i), V_th, V_ph, (PolMag(j), PolZen(j), PolAzi(j), PoldOm(j), j=1,3)  !V_hor, V_N, V_E, V_h!
               !write(2,*) 'ConstructTracks', PolMag(:), PolZen(:)
            EndIf
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
        !
        If(PolarAna) Then
            LeaderStks(1:6)=0
            twp=0
            Do k=1,TrackENr(i)
               j=TrackE(k,i) ! j= rank-number of k-th source on track i
               LeaderStks(1:6)=LeaderStks(1:6)+Stk_NEh(1:6,Iperm(j))*SrcWidth(Iperm(j))
               twp=twp+SrcWidth(Iperm(j))
               !write(2,"('!',4I4,10F8.3)") k,j,Iperm(j),SrcWidth(Iperm(j)),Real(LeaderStks(1)+LeaderStks(4)+LeaderStks(6))/twp, &
               !   Imag(LeaderStks(2))/twp,Imag(LeaderStks(3))/twp,Imag(LeaderStks(5))/twp
            EndDo
            LeaderStks(1:6)=LeaderStks(1:6)/twp
            prin=.true.   !  Print polarization observables
            write(2,*) 'Track averages integrated Stokes [= Width x (Stokes/sample)]:'
            Call PolPCACathCon(LeaderStks, PolZen, PolAzi, PolMag, PoldOm, prin)
            prin=.false.   !  Do not print polarization observables
         EndIf
        !
        Write(2,"(A,i2,I5,A,3F7.1)") 'Track#',i,TrackENr(i),', RMS(E,N,H)[m]:', &
         sqrt(Xsq/TrackENr(i))*1000., sqrt(Ysq/TrackENr(i))*1000., sqrt(Zsq/TrackENr(i))*1000.
        Close(unit=29)
        If(PolarAna) Close(unit=30)
        Write(2,*) 'ConstructTracks, File written:',TRIM(PlotFile)//trim(extension),' with ',TrackENr(i), ' data lines'
    enddo ! i=1,TrackNr
    !
   Return
End Subroutine ConstructTracks
!===============================================
Subroutine BinTracks(dt_MTL, RA, SourcTotNr, PlotFile)
   ! Construct mean tracks from assigned source locations on return in  LeaderPos
   ! Calculate lateral deviations of source locations from mean track
   Use constants, only : dp
   IMPLICIT none
   Real*8, intent(inout) :: dt_MTL
   Real(dp), Intent(in) :: RA(4,*)
   Integer, intent(in) :: SourcTotNr
   Character(len=*) :: PlotFile
   !Real*8 :: TrackPos(1:3,TrackLenMax)  !  TrackPos(1:3,i): position of the head of the i^th track
   !
   Integer, parameter :: Nmax_MTL=500
   Real*8 :: MeanTrackLoc(1:4,Nmax_MTL)   !, dt_MTL=0.2  ! [ms]
   Integer :: N_MTL, i_MTL, t_MTL, t_old, i,j, k, i_LongTr, count
   Character*8 :: extension
   Real*8 :: dt,dist,t0
   !
   count=0
   i_LongTr=0
   do i=1,TrackNr
      If(TrackENr(i).lt. LongTrack_Min) cycle
      i_LongTr=i_LongTr+1
      if(i_LongTr.ge.10) exit
      !write(extension,"(A2,i1,A4)") '_s',nxx,'.dat' !,&
      write(extension,"(i1,A4)") i_LongTr,'.plt' !,&
      OPEN(unit=27,FILE=TRIM(PlotFile)//'_bin_'//trim(extension),FORM='FORMATTED',STATUS='unknown')
      write(27,"('!',F5.3,'Binned, t [ms]',T24,'N [km]',T39,'E',T53,'h',T73,'v_N',T86,'v_E',T99,'v_h',T110,'v[km/ms]', &
         T124,'Ave_t[ms]')")   dt_MTL
1     Continue
      t_old=-99
      i_MTL=0  ! just counter for filled bins
      MeanTrackLoc(1:4,:)=0.
      !Call Flush(2)
      t0=RA(1,TrackE(1,i)) ! start time of this track
      Do k=1,TrackENr(i)  ! average location over bin with length dt_MTL
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
      If((i_LongTr.eq.1) .and. (i_MTL.lt.3) .and. (count.lt.5)) Then
         dt_MTL=dt_MTL/2.
         count=count+1
         goto 1
      ElseIf((i_LongTr.eq.1) .and. (i_MTL.lt.3)) Then
         Write(2,*) 'Too few bins at dt_MTL=',dt_MTL,'[ms] for the first long track'
         dt_MTL=-dt_MTL
         close(unit=27)
         return
      EndIf
      !write(2,*) MeanTrackLoc(1,1:t_MTL)
      !Call Flush(2)
      !i_old=0
      Do k=2,i_MTL
         !If(MeanTrackLoc(1,k).gt. 0.) then
            !If(t_old.gt.0) then
               dt=TOrder*(MeanTrackLoc(1,k)-MeanTrackLoc(1,k-1))
               dist=sqrt(SUM((MeanTrackLoc(2:4,k)-MeanTrackLoc(2:4,k-1))**2))
               write(27,"(I4,1x,4(1x,g14.8),3x,3(1x,g12.4),2x,g10.4,2x,g13.7)") k, &
                  MeanTrackLoc(1:4,k),(MeanTrackLoc(2:4,k)-MeanTrackLoc(2:4,k-1))/dt, dist/dt, &
                  (MeanTrackLoc(1,k)+MeanTrackLoc(1,k-1))/2.
            !Endif
            !t_old=k
         !Endif
      Enddo
      close(unit=27)
      Write(2,*) 'BinTracks; File:',TRIM(PlotFile)//'_bin_'//trim(extension),' with ',i_MTL, ' data lines'
      if(i_LongTr.eq.1) Then  ! write track for the first, may be used for interferometric tracing
         OPEN(unit=27,FILE=TRIM(PlotFile)//'.trc',FORM='FORMATTED',STATUS='unknown')
         Do k=1,i_MTL
            If(TOrder.gt.0) Then  ! order in increasing time in tNEh notation
               write(27,"(4(1x,g14.8),3x,3(1x,g12.4),2x,g10.4,2x,g13.7)") MeanTrackLoc(1,k), &
                  MeanTrackLoc(2:4,k)
            Else
               write(27,"(4(1x,g14.8),3x,3(1x,g12.4),2x,g10.4,2x,g13.7)") MeanTrackLoc(1,i_MTL-k+1), &
                  MeanTrackLoc(2:4,i_MTL-k+1)
            EndIf
         Enddo
         close(unit=27)
      Write(2,*) 'BinTracks; File:',TRIM(PlotFile)//'.trc',' with ',i_MTL, ' data lines'
      EndIf
   enddo
    !
   Return
End Subroutine BinTracks
!===============================================
Subroutine AnalyzeBurst(RA, SrcI20_r, SourcTotNr, PlotFile)
   ! Make a histogram of the source times in a track where each pulse is smeared with a gaussian
   ! with half width of t_resol, where this resolution changes. The histogram is stored in TBurst_trace(T)
   ! Then check how often >3.5 (N1) or >2.5 (N2) sources are within the t_resol
   Use constants, only : dp
   !Use TrackConstruct, only : TOrder, LongTrack_Min, TrackE, TrackENr
   use FFT, only : RFTransform_su, DAssignFFT, RFTransform_CF
   IMPLICIT none
   Real(dp), Intent(in) :: RA(4,*)
   Integer, intent(in) :: SourcTotNr !, Label(4,*)  !  SrcI20_r(Label(2,j)*Aweight+1.)  ! Include intensity in weight factor
   real, intent(in) :: SrcI20_r(*)  !  SrcI20_r = (Label(2,j)*Aweight+1.)  ! Include intensity in weight factor
   Character(len=*) :: PlotFile
   integer, parameter :: TB_max=131072, N_nu=TB_max/2  ! 65536 ! 65536=2^16 ; 32768=2^15 ! 2048=2^11 ! 262144=2^18 ! 131072=2^17
   real*8 :: t_resol,T_max, T_min, TimeDur, A, B, Norm
   Real*8 :: TBurst_trace(0:TB_max), d_nu, F_Th, nu, BinSize, ResScalFact
   Real*8, allocatable :: TBurst_trace1(:)
   real*8 :: dt1, A1, RT1, d1
   Complex(dp) :: Cnu(0:N_nu)
   integer :: i,j, k,N_k,NTSampl, T_C, T, T1, T2, N1, N2, NM2, NM3, tr, iT1
   Character*8 :: extension
    !
    i=1
    write(2,*) TOrder,RA(1,TrackE(1,i)),TrackE(1,i),RA(1,TrackE(TrackENr(i),i)),TrackE(TrackENr(i),i)
    extension=''
    ResScalFact=5.
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
      t_resol= SrcDensTimeResol  ! mili second
      If(t_resol .gt. TimeDur) t_resol=TimeDur
      write(2,*) 'AnalyzeBurst: t_resol ', t_resol, T_min, T_max
      NTSampl=0
      !
      Do tr=1,5
         t_resol=t_resol/ResScalFact
         NTSampl=int(TimeDur/t_resol)
         if(NTSampl .gt. TB_max) Then
            NTSampl=TB_max
         EndIf
         TBurst_trace=0.
         Norm=0.
         Do k=1,TrackENr(i) ! loop over all events in track i
             j=TrackE(k,i)   ! get event number of k-th track member
             T_C=Int((RA(1,j)-T_min)/t_resol)  ! central time for this averaging period
             norm=norm + (SrcI20_r(j)*Aweight+1.) ! Use intensity weighting
             T1=T_C-8
             if(T1.lt.0) T1=0
             T2=T_C+8
             if(T2.gt.NTSampl) T2=NTSampl
             If(T1.ge.T2) exit
             Do T=T1,T2 ! update the +/- 8 time bins around the central
                 TBurst_trace(T)=TBurst_trace(T)+ exp(-(T-(RA(1,j)-T_min)/t_resol)**2)*(SrcI20_r(j)*Aweight+1.)
             Enddo
         enddo
         TBurst_trace(:)=TBurst_trace(:)*TrackENr(i)/(norm)  !  /(t_resol*norm)
         A=MaxVal(TBurst_trace(:))
         If(tr.eq.1) Then
            dt1=t_resol
            A1=A
            TBurst_trace(:)= TBurst_trace(:)/A
            allocate( TBurst_trace1(0:NTSampl+2) )
            TBurst_trace1(0:NTSampl+2)=TBurst_trace(0:NTSampl+2)
         EndIf
         ! B=MinVal(TBurst_trace(:))
         N1=COUNT(TBurst_trace .gt. 3.5)
         N2=COUNT(TBurst_trace .gt. 2.5)
         NM2=COUNT(TBurst_trace .gt. A/2.)
         NM3=COUNT(TBurst_trace .gt. A/3.)
         !
         write(extension,"(i2.2,A4)") tr,'.plt' !,& SourceTimeDistribution (STD)
         OPEN(unit=27,FILE=TRIM(PlotFile)//'_STD-'//trim(extension),FORM='FORMATTED',STATUS='unknown')
         write(27,*) t_resol*1000., A1, A
         write(27,"('! ', T8,'time[ms]',T33,'PulseDensity',T49,'resol[ms]=',F9.6)") t_resol
         If(NTSampl.le.1) Then
            NTSampl=1
            TBurst_trace(1)=0.
         EndIf
         Do T=0,NTSampl
            RT1=T*t_resol/dt1
            it1=int(RT1)
            d1=RT1-iT1
            norm=((1.-d1)*TBurst_trace1(it1) + d1*TBurst_trace1(it1+1))*A1
            write(27,*) T_min+T*t_resol, TBurst_trace(T), TBurst_trace(T)/norm
            !TBurst_trace(T)=TBurst_trace(T)/norm
         enddo
         Close(unit=27)
         Write(2,"(A,i2,A,f7.4,A,F7.1,A,i5,A,i5,A,3I6,2A)") 'track #=',i,', sample=',t_resol,'[ms], maximum=',A, &
             ', > 3.5:',N1,', > 2.5:',N2,', #samples=',NTSampl,NM2,NM3,'; File:',TRIM(PlotFile)//'_STD-'//trim(extension)
      enddo  ! tr
      deallocate( TBurst_trace1)
      !
      Call RFTransform_su(NTSampl+1)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !
      Call RFTransform_CF(TBurst_trace,Cnu)
      OPEN(unit=27,FILE=TRIM(PlotFile)//'_STD-nu.plt',FORM='FORMATTED',STATUS='unknown')
      write(27,"('! ',T6,'freq[1/ms]',T31,'Ampl',T57,'Not rebinned')")
      d_nu=2./(t_resol*(NTSampl+1))
      F_Th=d_nu/2
      N_k=0
      A=0.
      B=0.
      BinSize=d_nu/128
      BinSize=BinSize/128
      write(2,*) 'AnalyzeBurst: frequency ', d_nu, d_nu*TB_max/2
      Do T=0,(NTSampl+1)/2
         nu=T*d_nu
         If(nu.gt.F_th) then
            write(27,*) B/N_k,A/N_k, nu, BinSize, N_k
            F_th=1.05*nu ! to re-bin results
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
      Write(2,*) 'AnalyzeBurst; File:',TRIM(PlotFile)//'_STD-nu.plt',' with <',(NTSampl+1)/2, ' data lines'
      Call DAssignFFT()         !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !
    enddo ! i, track number
end Subroutine AnalyzeBurst
!===================================================================
Subroutine ReadPreDefTrFile(PreDefTrackFile)
   !Use TrackConstruct, only : PreDefTrack, tNrS_max, PreDefTrackNr, PreDefTrackNr_Max, NrPreDefTrackPoints_Max, tNrS
   IMPLICIT none
   Character(len=100), Intent(IN) :: PreDefTrackFile
   Integer :: nxx, i_track, i_t, lab,i
   Character(len=5) :: Extension
   character(len=80) :: line
   !
   PreDefTrackNr=0
   tNrS(1:PreDefTrackNr_Max)=1  ! needed in GetPreDefTrPos
   Do i_track=1,PreDefTrackNr_Max
      write(Extension,"(i1,A4)") i_track,'.plt' !,&
      OPEN(UNIT=22,STATUS='old',ACTION='Read',FILE=trim(PreDefTrackFile)//trim(extension),IOSTAT=nxx)
      If(nxx .ne. 0) then
         write(2,*) 'No predefined track file:',trim(PreDefTrackFile)//trim(extension)
         exit
      Endif
      !write(2,*) 'ReadPreDefTr:',trim(PreDefTrackFile)//trim(extension)
      i_t=0
      Do i=1,NrPreDefTrackPoints_Max
         read(22,"(A80)",iostat=nxx) line
         If(nxx.ne.0) exit
         if(line(1:1).eq.'!') cycle
         i_t=i_t+1
         read(line,*,iostat=nxx) lab, PreDefTrack(0:3, i_t, i_track)
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
End Subroutine ReadPreDefTrFile

!=========================================
Subroutine GetPreDefTrPos(t_source,TrackPos)
   !Use TrackConstruct, only : TOrder, PreDefTrack, tNrS_max, PreDefTrackNr, PreDefTrackNr_Max,tNrS
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
End Module TrackConstruct
