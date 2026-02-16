Module PredefinedTracks
!   Use Tracks, only : PreDefTrack(0:3,tNr,1:PreDefTrackNr), tNrS_max(i_Pre), PreDefTrackNr
   Use constants, only : dp, CI, pi, c_l
   Character(len=100), save :: PreDefTrackFile=' '
   Integer, parameter :: NrPreDefTrackPoints_Max=900
   Real(dp), save :: PreDefTrack(NrPreDefTrackPoints_Max,0:3)
   Real(dp), save :: t_track
   Integer, save :: tNrS  ! number of points on the track
contains
!-------------------------------------
Subroutine GetTrPos(t_source, TrackPos)
   use constants, only : dp
   IMPLICIT none
   !Character(len=100), intent(in) :: TrackFile
   Real(dp),intent(IN) :: t_source
   Real(dp),intent(OUT) :: TrackPos(1:3)
   Real(dp) ::DtP, DtM
   !Integer :: i_track=0
   Integer :: LowPos
   Call ReadPreDefTr()
   !Call GenSort( PreDefTrack(1:tNrS,0:3) )
   !write(2,*) PreDefTrack(1:tNrS,0)
   If(tNrS.lt.2) Then
      write(2,*) 'Predefined track is too short, length=',tNrS,TRIM(PreDefTrackFile)
      stop 'too short predefined track'
   EndIf
   Call OrderedSearch(PreDefTrack(1,0), tNrS, t_source, LowPos)
   DtM= (t_source - PreDefTrack(LowPos,0))
   DtP= (PreDefTrack(LowPos+1,0)-t_source)
   !write(2,*) 'DtM:',DtM,tNr,t_source, PreDefTrack(0,tNr,i_Pre)
   TrackPos(1:3)= &
      (PreDefTrack(LowPos,1:3)*DtP + PreDefTrack(LowPos+1,1:3)*DtM)/(DtM+DtP) ! average between earlier and later positions
   Return
End Subroutine GetTrPos
!======================================
Subroutine ReadPreDefTr( i_track)
!  i_track is a possible label, to keep consistency with "Track.f90"
!Subroutine ReadPreDefTr(TrackFile, PreDefTrack, tNrS_max, NrPreDefTrackPoints_Max, i_track)
!      NrPreDefTrackPoints_Max=900   ! Maximum number of defining times per Pre-defined Track
!        tNrS(1:PreDefTrackNr_Max)=1  ! needed in GetPreDefTrPos, not sure why
!   Use Tracks, only : PreDefTrack, tNrS_max, PreDefTrackNr, PreDefTrackNr_Max, NrPreDefTrackPoints_Max, tNrS
   use constants, only : dp
   IMPLICIT none
   !Character(len=100), intent(in) :: TrackFile
   !Real(dp), intent(out) :: PreDefTrack(1:NrPreDefTrackPoints_Max,0:3)
   !Integer, intent(out) :: tNrS_max
   Integer,optional, intent(in) :: i_track ! , NrPreDefTrackPoints_Max
   Integer :: nxx, indx, i_t, Len
   Character(len=5) :: Extension
   !Character(len=100) :: TrackFile
   Real(dp) :: t_ms, NEh(1:3)
   Logical :: ENh
   !
   !TrackFile=PreDefTrackFile
   Len=LEN_TRIM(PreDefTrackFile)
   Extension=''
   ENh=.false.
   if((PreDefTrackFile(Len-3:Len) .eq. '.dat')) Then  ! use a plot file
      ENh=.true.
      OPEN(UNIT=22,STATUS='old',ACTION='Read',FILE=trim(PreDefTrackFile),IOSTAT=nxx)
      read(22,*,iostat=nxx) t_ms
      read(22,*,iostat=nxx) t_ms
   Elseif((PreDefTrackFile(Len-3:Len) .eq. '.trc')) Then ! Use a special track file
      OPEN(UNIT=22,STATUS='old',ACTION='Read',FILE=trim(PreDefTrackFile),IOSTAT=nxx)
   Elseif((PreDefTrackFile(Len-3:Len) .eq. '.plt')) Then ! Use a more modern plot file
      OPEN(UNIT=22,STATUS='old',ACTION='Read',FILE=trim(PreDefTrackFile),IOSTAT=nxx)
   ElseIf(present(i_track)) Then
      write(Extension,"(i1,A4)") i_track,'.plt' !,&
      OPEN(UNIT=22,STATUS='old',ACTION='Read',FILE=trim(PreDefTrackFile)//trim(extension),IOSTAT=nxx)
   Else
      OPEN(UNIT=22,STATUS='old',ACTION='Read',FILE=trim(PreDefTrackFile),IOSTAT=nxx)
   EndIf
   If(nxx .ne. 0) then
      write(2,*) 'No predefined track file:',trim(PreDefTrackFile)//trim(extension)
      stop
   Endif
   !
   Do i_t=1,NrPreDefTrackPoints_Max
      If(ENh) Then
         read(22,*,iostat=nxx) indx, t_ms, NEh(:)
         PreDefTrack(i_t, 0)=t_ms
         PreDefTrack(i_t, 1)=NEh(2)*1000.d0
         PreDefTrack(i_t, 2)=NEh(1)*1000.d0
         PreDefTrack(i_t, 3)=NEh(3)*1000.d0
      Else
         read(22,*,iostat=nxx) t_ms, NEh(:)
         Call Convert2m(NEh(:))
         PreDefTrack(i_t, 0)=t_ms
         PreDefTrack(i_t, 1:3)=NEh(1:3)
      EndIf
      !write(2,*) PreDefTrack(0:3, i_t, i_track), i_t, i_track
      If(nxx.ne.0) exit
   Enddo
   Close(UNIT=22)
   write(2,*) 'ReadPreDefTr: ',trim(PreDefTrackFile)//trim(extension), i_t-1
   tNrS=i_t-1
   !Enddo
   !
   Return
End Subroutine ReadPreDefTr
!=========================================
Subroutine OrderedSearch(SearchArr, High, SearchEl,LowPos)
! It is not assumed that values are in increasing order of magnitude, but the values should be ordered.
   use constants, only : dp
   IMPLICIT none
   Real(dp),intent(IN) :: SearchArr(*), SearchEl
   Integer, intent(in) :: High
   Integer, intent(out) :: LowPos
   Integer :: Low, Upp, middle, TOrder, cnt
   TOrder=+1
   cnt=1
   If(SearchArr(1) .gt. SearchArr(High)) TOrder=-1  ! check for reverse ordering
   IF (TOrder*SearchEl .le. TOrder*SearchArr(1)) Then
      write(2,*) SearchEl,' before track start @', SearchArr(1)
      stop 'Track-start problem'
   EndIf
   IF (TOrder*SearchEl .ge. TOrder*SearchArr(High)) Then
      write(2,*) SearchEl,' after track end @', SearchArr(High)
      stop 'Track-end problem'
   EndIf
   Low=1
   Upp=High
   DO WHILE(low .lt. Upp )
      middle = (low + Upp)/2
      cnt=cnt+1
      ! WRITE(2,*)"Now searching element: ", middle, low, Upp, TOrder*SearchArr(middle),cnt, SearchEl
      IF (TOrder*SearchEl .lt. TOrder*SearchArr(middle)) THEN
         Upp = middle
      ELSEIF (TOrder*SearchEl .lt. TOrder*SearchArr(middle+1)) THEN !SearchEl between SearchArr(middle) and SearchArr(middle+1)
         exit
      Else
         low = middle
      ENDIf
      If (cnt.gt.100) exit
   END DO
   LowPos=middle
   WRITE(2,*)"time ", SearchEl,' between times-on-track of', SearchArr(middle), SearchArr(middle+1), &
      ', t-order=', TOrder, ', cnt=', cnt
   If((TOrder*SearchEl .lt. TOrder*SearchArr(middle)) .or. (TOrder*SearchEl .gt. TOrder*SearchArr(middle+1))) Then
      write(2,*) 'interpolation problem', cnt
      stop 'Sorting problem'
   EndIf
   Return
End Subroutine OrderedSearch
!=========================================
Subroutine GetTrTime(t_old, t_new, TrackPos)
!  Get the time on track corresponding to a track-source that is at a certain distance from CenLoc
   use constants, only : dp
   Use Interferom_Pars, only : N_pix, d_loc, CenLoc
   IMPLICIT none
   Real(dp),intent(IN) :: t_old
   Real(dp),intent(OUT) :: t_new, TrackPos(1:3)
   !Integer, parameter :: NrPreDefTrackPoints_Max=900
   !Real(dp) :: PreDefTrack(NrPreDefTrackPoints_Max,0:3)
   Real(dp) :: DtP, DelMax(1:3), Del(1:3), d,  DtM
   Integer :: i_track, j, k
   Integer :: LowPos, Niter
   Call ReadPreDefTr()
! assume Cathesian coordinate system
!  maximum deviation from old position
   Do j=1,3
      DelMax(j)=(N_pix(j,2)-N_pix(j,1)-3)*d_loc(j) ! -4 to have some overlap of boxes
   EndDo
   !
   Call OrderedSearch(PreDefTrack(1,0), tNrS, t_old, LowPos)
   !DtM= (t_guess - PreDefTrack(LowPos,0))
   !DtP= (PreDefTrack(LowPos+1,0)-t_guess)
   !Del(1:3)= ABS(CenLoc(1:3) -&
   !   (PreDefTrack(LowPos,1:3)*DtP + PreDefTrack(LowPos+1,1:3)*DtM)/(DtM+DtP) )! average between earlier and later positions
   DtP=0.d0
   Niter=20
   !write(2,*) 'cenloc', LowPos, CenLoc(1:3), DelMax(1:3)
   Do i_track=LowPos,tNrS-1  ! step along the track till distance becomes too large
      ! assume Cathesian coordinate system:
      Del(1:3)= ABS(CenLoc(1:3) - PreDefTrack(i_track+1,1:3)) !( track position for DtP=0
      j=MaxLOC(Del(1:3)-DelMax(1:3),1)    ! find largest step (compared to max allowed) (voxels should have overlap)
      !  largest step should not exceed max allowed
      !         write(2,*) j,(Del(1:3)-DelMax(1:3))
      !         write(2,*) i_track, (PreDefTrack(i_track,1:3)-CenLoc(1:3)), (PreDefTrack(i_track+1,1:3)-CenLoc(1:3))
      If((Del(j)-DelMax(j)) .gt. 0.) Then ! largest step should not exceeds max allowed ==> finetune step
         If(i_track.eq.LowPos) Then  ! step between track-points is larger than grid-size
            DtP=(t_old-PreDefTrack(i_track+1,0))/(PreDefTrack(i_track,0) - PreDefTrack(i_track+1,0))
            d=1.d0-DtP
            write(2,*) 'same track point', DtP,d
         Else
            d=1.d0
         EndIf
         Do k=1,Niter
            ! assume Cathesian coordinate system:
            Del(1:3)= ABS(CenLoc(1:3) -  &
              (PreDefTrack(i_track,1:3)*DtP + PreDefTrack(i_track+1,1:3)*(1.-DtP)) )! average between earlier and later positions
            j=MaxLOC(Del(1:3)-DelMax(1:3),1)    ! find largest step (compared to max allowed) (voxels should have overlap)
            d=d/2.d0
            If((Del(j)-DelMax(j)) .gt. 0.) Then ! largest step should not exceeds max allowed ==> decrease step
               DtP=DtP+d  !  Move closer to earlier track position
            ElseIf((Del(j)-DelMax(j)) .lt. -d_loc(j)) Then ! largest step could be made larger ==> increase step
               DtP=DtP-d  !  Move closer to later track position
            Else !Converged
               write(2,*) 'converged', j,(Del(1:3)-DelMax(1:3)),d_loc(1:3),k
               write(2,*) DtP, (PreDefTrack(i_track,1:3)-CenLoc(1:3)), (PreDefTrack(i_track+1,1:3)-CenLoc(1:3))
               exit
            EndIf
         EndDo
         DtM=1.d0-DtP
         !converged
         If(k.ge.Niter) Then
            write(2,*) 'Not really converged; exit', j,(Del(1:3)-DelMax(1:3)),d_loc(1:3)
            write(2,*) k,Niter,d,DtP
            write(2,*) DtP, (PreDefTrack(i_track,1:3)-CenLoc(1:3)), (PreDefTrack(i_track+1,1:3)-CenLoc(1:3))
         EndIf
         exit
      EndIf
      !If(ANY(Del(1:3).gt.DelMax(1:3))) Then !  shorten step
   EndDo
   t_new=PreDefTrack(i_track,0)*DtP + PreDefTrack(i_track+1,0)*(1.d0-DtP)
   TrackPos(1:3)= PreDefTrack(i_track,1:3)*DtP + PreDefTrack(i_track+1,1:3)*(1.d0-DtP) ! average between earlier and later positions
   Return
End Subroutine GetTrTime
End Module PredefinedTracks
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=    !
