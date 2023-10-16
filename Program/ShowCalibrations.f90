    Include 'ConstantsModules.f90'
    !Include 'ParamModules.f90'
    !Include 'MappingUtilitiesModules.f90'
    !Include 'FFT_routines.f90'
    Include 'System_Utilities.f90'
    Include 'GLEplotUtil.f90'
!-----------------------------------
!=================================
!Module Utilities
!contains
!    Include 'MappingUtilities.f90'
!End Module Utilities
!-----------------------------------
Program ShowCalibr
   use constants, only : dp,sample!,Refrac,c_mps
   use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, DataFolder, FlashName, Windows, RunMode
   use GLEplots, only : GLEplotControl
   !Use Utilities, only : GetNonZeroLine
   Implicit none
   !
   Integer, parameter :: ChunkNr_dim=50, PeakNr_dim=500, lnameLen=180
   Real(dp) :: T1, NEh(1:3), xMin, xMax, yMin,yMax,zMin,zMax,tMin,tMax, t
   integer :: i_eo, i_chunk, i_dist, i, j, k, i_c,i_ca,i_eoa, i_d, n_shft, i_s, i_pos, i_peak
   Real(dp) ::  SourcePos(1:3,PeakNr_dim)
   Integer :: ChunkNr(PeakNr_dim),  ChunkStTime_ms(ChunkNr_dim), PeakNrTotal, n_chunk, nxx
   character*150 :: SystemCommand
   Character*8 :: extension,ZoomBox
   character*100 :: pars, datfile
   Character*20 :: Utility, release
   Character(len=180) :: lname
   character*10 :: txt
   character*1 :: Label
   character*35 :: Frmt
   !
   Open(unit=2,STATUS='unknown',ACTION='write', FILE ="ShowCal.out")
   Utility='ShowCalibration'
   release='v22'
   Call System_Initiation(Utility, release, ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows)  ! set Flash name & folder names
   RunMode=5  !  used in the gle-plotting routine
   ZoomBox='NoBox'
   !t_start=0.
   !
   Call GetNonZeroLine(lname)
   read(lname,*,iostat=nxx)  datfile
   datfile=TRIM(FlashName)//datfile
   write(2,*) 'datfile:',datfile
   flush(unit=2)
   i_chunk=0
   n_chunk=0
   i_ca=-1
   i_eoa=-1
   i_peak=0
   t=0.
   !
   Frmt="(A1,3i2,I8,3(F10.2,1x),F8.2)"
   Call GetNonZeroLine(lname)
   Do
      Read(lname(2:lnameLen),*,iostat=nxx) T1,NEh(:)  ! just dummy arguments
      !write(2,*)  nxx, lname
      !flush(unit=2)
      If(nxx.ne.0) exit  ! Neither a new-chunk line or a pulse-info line
      read(lname,Frmt,iostat=nxx)  Label,i_s,i_eo,i_c, i_pos,  NEh(:) !,t
      !write(2,*) 'pulse read;', nxx, Label,i_s,i_eo,i_c, i_pos,  NEh(:)
      If(nxx.eq.0 .and. (i_s.ge.0) .and. (i_eo.ge.0 .and. i_eo.le.2) .and. (i_c.gt.0) .and. (i_pos.gt.1000) ) Then ! check for pulse-info
         !  Determine number of peaks/sorces that are included in the calibration search
         i_peak=i_peak+1 ! always assume both polarities for FieldCalibrate
         If(i_peak .gt. PeakNr_dim) then
            write(2,*) 'PeakNr_dim will be exceeded for #',i_peak
            Write(2,*) 'Expected:   Core/Reference antenna |i |i_eo |i_c |k |N |E |h | dt " in format ',FrMT
            Write(2,*) 'Culprit: "',trim(lname),'"'
            stop 'ReadPeakInfo problem'
         endif
         If(i_ca.ne.i_c) then
             i_chunk=i_chunk+1
         ElseIf(i_eoa .ne. i_eo) then
             If(i_eo.eq. 0) i_chunk=i_chunk+1
         Endif
         i_eoa=i_eo
         i_ca=i_c
         If(i_eo.eq.2) then
            i_eo=0 ! same peak info for even&odd
         EndIf
         If(i_chunk .gt. n_chunk)  Then
            write(2,*) 'EIReadPeakInfo: Chunk-info line for chunk',i_chunk,' should preceed pulse-info', n_chunk,lname
            stop 'Chunk-info line should preceed pulse-info lines'
         EndIf
         SourcePos(:,i_Peak)=NEh(:)
         ChunkNr(i_Peak)=i_chunk
         PeakNrTotal=i_peak
         !write(2,*) 'i_peak,i_eo,i_c',i_peak,i_eo,i_c,i_eoa,i_ca,';',PeakNrTotal
         Call GetNonZeroLine(lname)
         read(lname,* ,iostat=nxx)  txt
         If(trim(txt).eq.'exclude') Then
            Call GetNonZeroLine(lname)
         EndIf
      Else
         n_chunk=n_chunk+1   ! this was a genuine chunk card
         If(n_chunk .gt. ChunkNr_dim) Then
            write(2,*) 'ChunkNr_dim=',ChunkNr_dim,' will be exceeded'
            stop 'ChunkNr_dim'
         EndIf
         !write(2,*) 'chunk',n_chunk,NEh(:),T1
         Call Convert2m(NEh(:))
         ChunkStTime_ms(n_chunk)=T1
         Call GetNonZeroLine(lname)
         i_ca=-1
      EndIf
   EndDo

   xMin=MINVAL(SourcePos(2,1:PeakNrTotal))/1000.-1.
   xMax=MAXVAL(SourcePos(2,1:PeakNrTotal))/1000.+1.
   yMin=MINVAL(SourcePos(1,1:PeakNrTotal))/1000.-1.
   yMax=MAXVAL(SourcePos(1,1:PeakNrTotal))/1000.+1.
   !write(2,"(6(1x,f8.3),2(1x,f8.3),1x,A,f9.3,A)") xMin, xMax, yMin,yMax
   T1=MAX((xMax-xMin), (yMax-yMin), 100.d0)/2.
   t=(yMax+yMin)/2.
   yMax=t+T1
   yMin=t-T1
   t=(xMax+xMin)/2.
   xMax=t+T1
   xMin=t-T1
   !write(2,"(6(1x,f8.3),2(1x,f8.3),1x,A,f9.3,A)") xMin, xMax, yMin,yMax,T1
   zMin=MINVAL(SourcePos(3,1:PeakNrTotal))/1000.-1.
   zMax=MAXVAL(SourcePos(3,1:PeakNrTotal))/1000.+1.
   tMin=MINVAL(ChunkStTime_ms(1:n_chunk))-1.
   tMax=MAXVAL(ChunkStTime_ms(1:n_chunk))+1.


   OPEN(unit=29,FILE=TRIM(DataFolder)//trim(datfile)//'.dat',FORM='FORMATTED',STATUS='unknown')  !
   write(29,"(6(1x,f8.3),2(1x,f8.3),1x,A,f9.3,A)") xMin, xMax, yMin,yMax,zMin,zMax,tMin,tMax,trim(ZoomBox),0.,'  ! '
   write(29,"(1x,I3,I8,F7.2,F8.4,2x,A,1x,A)") 0, PeakNrTotal, 0.0, 0.0, &
       TRIM(FlashName)//'/'//TRIM(datfile),  '0 0 0 0 00 0  0 0 0 0 0 0 0 0 0 0  !'
   Do i_Peak=1,PeakNrTotal  !  finally write all accepted source to file
      write(29,"(1x,i8,4(2x,g14.8),3x,F7.2)")  ChunkNr(i_Peak),ChunkStTime_ms(ChunkNr(i_Peak)), &
         SourcePos(2,i_Peak)/1000.,SourcePos(1,i_Peak)/1000.,SourcePos(3,i_Peak)/1000., 1.0
   enddo
   close(unit=29)
   !
   Call GLEplotControl(PlotType='SourcesPlot', PlotName='CalSrcs_'//TRIM(datfile), &
         PlotDataFile=TRIM(DataFolder)//trim(datfile), Submit=.true.)
   !
   Stop
   !
End Program ShowCalibr
!===============================================
!===============================================
    Subroutine GetNonZeroLine(lineTXT)
    implicit none
    character(LEN=*), intent(inout) :: lineTXT
    Character(len=6) :: FMT
    integer :: eof,N
      N=LEN(lineTXT)
      write(FMT,"(A,I3.3,A)") '(A',N,')'
      !write(2,*) 'GetNonZeroLine2:', FMT, N
      lineTXT=''
      Do while (trim(lineTXT).eq.'')
         read(*,FMT,iostat=eof) lineTXT
         if(eof.lt.0) exit
      enddo
    End Subroutine GetNonZeroLine
!==========================================
Subroutine Convert2m(CenLoc)
   Implicit none
   real*8, intent(inout) :: CenLoc(3)
   If((abs(CenLoc(1)) .lt. 180.) .and. (abs(CenLoc(2)) .lt. 180.) .and. (abs(CenLoc(3)) .lt. 20.) ) Then
      CenLoc(:)=CenLoc(:)*1000.  ! convert from [km] to [m]
   Endif
   Return
End Subroutine Convert2m
!===============================================
