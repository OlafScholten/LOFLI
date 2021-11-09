! Interferometry Sources Selction
!
   Include 'TrackModule.f90'
   Include 'Fitter_Ampl.f90'
   Include 'System_Utilities.f90'
   Include 'GLEplotUtil.f90'
!
Module IntfSrcSelMod
   Use constants, only :  dp
   Real(dp), save :: xMin,xMax,yMin,yMax,zMin,zMax,tMin,tMax, SMPowCut=0., RatMaxCut=1.0 ! RatMax=ratio (Max rim)/peak
   Real(dp), save :: Tiny=epsilon(Tiny), Small
!contains
End Module IntfSrcSelMod
    !
PROGRAM InterfSrcSel
! Make an additional selection to the sources found by the interferometric method
   Use constants, only :  dp
   Use Tracks, only : maxd, RA, Label, EventNr
   Use Tracks, only : Nmax_Ampl, AmplScale, d_AmplScale, MaxAmplFitPercent, AmplitudeFit, HPSORT_mult_RI
   ! Nmax_Ampl: Maximum number of bins in histogram, determines max amplitude that can be fitted.
   ! AmplScale: factor (>1) used to convert amplitudes to integers as needed for Label(2,:)
   ! d_AmplScale: bin-size for histogramming the scaled amplitudes before fitting
   ! MaxAmplFitPercent: Maximum percentile of sources in overflow bin
   Use Tracks, only : Nrm, a1,b1,c1,d1,ChiSq1,Nrm1, a2,b2,c2,d2,ChiSq2,Nrm2, a3,b3,c3,d3
   ! extracted parameters of designated fit functions
   ! Presently: 1== (modified) exponential; 2== (modified)powerlaw; 3=not used
   !
   use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, DataFolder, FlashName, Windows, RunMode
   use GLEplots, only : GLEplotControl
   Use IntfSrcSelMod
   IMPLICIT none
   !character*150 :: SystemCommand
   Character*8 :: extension
   Character*20 :: BoxOption, Utility, release
   Character*3 :: BarMx(1:2)=(/ 'Bar', 'Mx ' /)
   Integer, parameter :: NDataFil=50
   character*100 :: datafile(NDataFil),datfile, OutFileLabel, SelFileName, PlotName, FileLabel, shellin
   Character(len=180) :: lineTXT, OS
   integer :: i,j, nxx, SourcTotNr
   Real*8 :: t_start
   Real(dp) :: t,x,y,z,SMPow, RatMax ! , Tiny=epsilon(Tiny)
   Real(dp) :: Q, IntScale=1000., AmpltPlot, AmpltPlotRead
   Real(dp) :: xMinR,xMaxR,yMinR,yMaxR,zMinR,zMaxR,tMinR,tMaxR, t_offset
   integer :: i_eo,i_BM,i_datfil, valueRSS
   Logical :: ZoomClip=.false.
   Logical :: First28=.true.
   NAMELIST /Parameters/ OutFileLabel, datafile, &
      xMin,xMax,yMin,yMax,zMin,zMax,tMin,tMax, ZoomClip, SMPowCut, RatMaxCut, AmpltPlot, MaxAmplFitPercent
   !
   OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='IntfSrcSel.out')
   !write(2,*) 'Program InterfSrcSel-v19'
   Utility='Interf.Select'
   release='V-19'
   Call System_Initiation(Utility, release, ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows)
   Small=sqrt(Tiny)
   RunMode=6  !  used in the gle-plotting routine
   !
   Do
      AmpltPlot=Tiny ! set it to some strange value, just to determine if it has been set on input
      xMin=Tiny; xMax=Tiny; yMin=Tiny; yMax=Tiny; zMin=Tiny; zMax=Tiny; tMin=Tiny; tMax=Tiny
      datafile(:)=''
      MaxAmplFitPercent=0.1
      AmplScale=1./IntScale
      read(*,NML = Parameters,IOSTAT=nxx)
      If(nxx.ne.0) Then
         Write(2,*) '======================  End input'
         exit
      EndIf
      If(SMPowCut.le.0.01) SMPowCut=0.01
      If(RatMaxCut.lt.1.) RatMaxCut=1.
      !write(2,*) 'time offset=',t_start,'[ms]', RatMaxCut
      RatMaxCut=1./RatMaxCut
      !
      Do i_BM=1,2
         Do i_eo=0,2
            Write(extension,"(A,'_',i1)") TRIM(BarMx(i_BM)),i_eo
            If(i_eo .eq. 2) extension=TRIM(BarMx(i_BM))//'_d'
            SourcTotNr=0
            write(2,*) '======================  ',extension
            If(.not. Windows) call system_mem_usage(valueRSS)
            First28=.true.
            Do i_datfil=1,NDataFil
               If(datafile(i_datfil).eq."") exit
               datfile=TRIM(DataFolder)//trim(datafile(i_datfil))//'IntfSpecPow'//TRIM(extension)//'.dat'  ! IntfSpecPowMx
               OPEN(unit=28,FILE=trim(datfile), FORM='FORMATTED',STATUS='OLD',ACTION='READ', IOSTAT=nxx) ! space separated values
               If(nxx.ne.0) then
                  Write(2,*) 'Problems opening file:"',trim(datfile),'"'
                  Close(Unit=28)
                  cycle
               endif
               Read(28,"(A180)",IOSTAT=nxx) lineTXT
               If(nxx.ne.0) stop 'InterfSrcSel-readingC'
               read(lineTXT,*,IOSTAT=nxx) xMinR, xMaxR, yMinR, yMaxR, zMinR, zMaxR, tMinR, tMaxR, BoxOption, t_start     ! from the first line that does not start with !
               Call SelectBB(xMinR, xMaxR, yMinR, yMaxR, zMinR, zMaxR, tMinR, tMaxR)
               If(First28) Then
                  t_offset=t_start
                  First28=.false.
               EndIf
               Read(28,*) i, j, Q, SMPow, FileLabel, AmpltPlotRead, a1, b1, c1, a2, b2, c2, d2
               If(AmpltPlot.eq.Tiny) AmpltPlot=AmpltPlotRead
               !Write(2,*)  Q, SMPow, TRIM(FileLabel), AmpltPlot, a1, b1, c1, a2, b2, c2, d2
               !write(2,*) TRIM(BoxOption)
               If(TRIM(BoxOption).eq.'IntfBox') then
                  !write(2,*) TRIM(BoxOption)
                  Do i=1,8 ! skip 8 records
                     Read(28,"(A180)",IOSTAT=nxx) lineTXT
                  Enddo
               EndIf
               !
               t_start = t_start-t_offset
               do
                  NXX=0
                  read(28,*,iostat=nxx)  i,t,x,y,z,SMPow, RatMax  ! already in proper units for plotting
                  t=t+t_start  ! Same t-scale now as for first file
                  If(nxx.gt.0) then
                     write(2,*) 'Read error:',i,t,x,y,z,SMPow
                  ElseIf (nxx.lt.0) Then
                     !Write(2,*) 'EOF reached'
                     exit
                  Endif
                  If(ZoomClip) then
                     if(x.le.xMinR .or. x.ge.xMaxR) cycle
                     if(y.le.yMinR .or. y.ge.yMaxR) cycle
                     if(z.le.zMinR .or. z.ge.zMaxR) cycle
                     if(t.le.tMinR .or. t.ge.tMaxR) cycle
                  EndIf
                  If( SMPow .lt. SMPowCut) cycle
                  !If(RatMax.gt.RatMaxCut) cycle ! allow only those sources where (background/peak < 1/RatMaxCut)
                  !write(2,*) i,x,xMin, xMax
                  SourcTotNr=SourcTotNr+1
                  RA(1,SourcTotNr)=t ;  RA(2,SourcTotNr)=x ;  RA(3,SourcTotNr)=y ;  RA(4,SourcTotNr)=z   ! units [ms], [km], [km], [km]
                  Label(1,SourcTotNr)=i
                  Label(2,SourcTotNr)=SMPow*IntScale
                  if(SourcTotNr.eq.maxd) then
                     write(*,*) 'Max. dimension reached of', maxd,' at',i
                     write(2,*) 'Max. dimension reached of', maxd,' at',i,t
                     exit
                  EndIf
               enddo
               Close(Unit=28)
               Write(2,*) 'After file:',trim(datfile), ', SourcTotNr=',SourcTotNr
            Enddo
            If(SourcTotNr.le.0) cycle
            !
            !
            CALL HPSORT_mult_RI(RA,Label,SourcTotNr)  ! RA(1:4,*) are the accepted source locatopns
            write(2,*) 'Time span data=', RA(1,1), RA(1,SourcTotNr)
            ! write(*,*) 'Sorting Done'
            !
            Flush(unit=2)
            !
            If(SMPowCut.gt.1.) then
               d_AmplScale=AmplScale !(If too small there are rounding-off problems since there is a conversion to interegs involved)
            Else
               d_AmplScale=AmplScale/SMPowCut  !  Needed for fitting amplitude spectrum
            EndIf
            SelFileName=TRIM(DataFolder)//'AmplFit'//trim(datafile(1))//TRIM(OutFileLabel)//TRIM(extension)
            PlotName=trim(datafile(1))//TRIM(extension)//TRIM(OutFileLabel)//'AmplFit'
            Call AmplitudeFit(SourcTotNr, SelFileName=SelFileName)
            write(2,"(1x,A, f6.2, A, 2g10.4, g11.3, A,2g11.3)") '$b *exp(-a*A-c/A^2); \chi^2=$',ChiSq1, &
                                    ',with  a,b,c= ',a1,b1/Nrm,c1,'; nrm=',Nrm1, Nrm
            write(2,"(1x,A, f6.2, A, 2g10.4, g11.3, A,2g11.3)") '$b *A^-a *exp(-c/A); \chi^2=$',ChiSq2, &
                                    ',with  a,b,c= ',a2,b2/Nrm,c2,'; nrm=',Nrm2, Nrm
            Call GLEplotControl(PlotType='Intensity', PlotName=trim(PlotName), &
                  PlotDataFile=trim(SelFileName), Submit=.false.)
            Flush(unit=2)
            !
            SelFileName=TRIM(DataFolder)//'IntfSpecSel'//trim(datafile(1))//TRIM(OutFileLabel)//TRIM(extension)
            OPEN(unit=28,FILE=TRIM(SelFileName)//'.dat', FORM='FORMATTED',STATUS='unknown',ACTION='write') ! space separated values
            write(28,"(6F9.3,2F9.3,A,F7.1,' 0 ')") xMinR, xMaxR, yMinR, yMaxR, zMinR, zMaxR, &
                                    tMinR, tMaxR, ' NoBox ', t_offset         ! gleRead:  xMin xMax yMin yMax zMin zMax tMin tMax ZmBx$ t_offst
            write(28,"(1x,'0 ',I8,F9.2, ' 0 ',2x,A,1x,F6.2, 6(1x,g10.4),f7.3,1x,A)")  &
               SourcTotNr, SMPowCut, trim(datafile(1))//TRIM(OutFileLabel)//TRIM(extension), AmpltPlot,  &
               a1,b1,c1,a2,b2,c2,d2, ' ! by InterfSrcSel' ! gleRead:  NTracks EventNr Q t_start label$ AmplitudePlot a1 b1 c1 a2 b2 c2 d2
            !write(29,"(1x,I3,I8,F7.2,F8.4,2x,A,1x,F6.2, 2(f7.3,f8.0,f7.3),f7.3,1x,A)") LongTrackNr, EventNr, Qual, t_start, &
            !    TRIM(Directory)//'/'//TRIM(datfile), AmplitudePlot, a1,b1,c1,a2,b2,c2,d2, ' !'
            Do j=1,SourcTotNr  !  finally write all accepted source to file
               write(28,"(1x,i8,4(2x,g14.8),1x,F9.2,2x,I3,2x,I3)")  Label(1,j),RA(1:4,j), Label(2,j)/IntScale
            enddo
            close(unit=28)
            !
            PlotName=trim(datafile(1))//TRIM(extension)//TRIM(OutFileLabel)//'IntfSpecSel'
            write(2,*) 'number of sources in plots:',SourcTotNr,SMPowCut,trim(PlotName)
            !
            If(SourcTotNr.le.0) then
               write(*,*) 'no sources for ',trim(SelFileName)
               cycle
            End If
            Call GLEplotControl(PlotType='SourcesPlot', PlotName=trim(PlotName), &
                  PlotDataFile=trim(SelFileName), Submit=.false.)
         EndDo ! i_eo
      Enddo ! i_BM
   Enddo
   !write(10,*) 'pause'

999 continue
   Call GLEplotControl( Submit=.true.)
    stop
!--------------------------------------------------
End Program InterfSrcSel
!===============================================
Subroutine SelectBB(xMinR, xMaxR, yMinR, yMaxR, zMinR, zMaxR, tMinR, tMaxR)
   Use IntfSrcSelMod, only : dp
   Use IntfSrcSelMod, only : xMin,xMax,yMin,yMax,zMin,zMax,tMin,tMax, Tiny, Small
   implicit none
   Real(dp), intent(inout) :: xMinR, xMaxR, yMinR, yMaxR, zMinR, zMaxR, tMinR, tMaxR
   If(xMin.ne.Tiny) xMinR=xMin
   If(xMax.ne.Tiny) xMaxR=xMax
   If(yMin.ne.Tiny) yMinR=yMin
   If(yMax.ne.Tiny) yMaxR=yMax
   If(zMin.ne.Tiny) zMinR=zMin
   If(zMax.ne.Tiny) zMaxR=zMax
   If(tMin.ne.Tiny) tMinR=tMin
   If(tMax.ne.Tiny) tMaxR=tMax
   !write(2,*) 'Box=', xMinR, xMaxR, yMinR, yMaxR, zMinR, zMaxR, tMinR, tMaxR
   !write(2,*) 'test', Tiny,xmin, Small
End Subroutine SelectBB
! ----------------------------------------------------------------------------------------
!     shellin = 'gle /d pdf '//trim(dummy3(i))//'.gle'
!     CALL system(shellin)
!     shellin = 'epstopdf '//trim(dummy3(i))//'.eps'
!     call system(shellin)
!===============================================
Subroutine system_mem_usage(valueRSS)
!          * VmPeak: Peak virtual memory size.
!          * VmSize: Virtual memory size.   total program size
!  VMsize is the "address space" that the process has in use: the number of available adresses. These addresses do not have to have any physical memory attached to them. ( Attached physical memory is the RSS figure)
!  VmSize includes RSS, plus things like shared libraries and memory mapped files (which don't actually use RAM), as well as allocated, but unused, memory.
!          * VmLck: Locked memory size (see mlock(3)).
!          * VmHWM: Peak resident set size ("high water mark").
!          * VmRSS: Resident set size.
!  VmRSS is the measure of how much RAM the process is actually using.
!          * VmData, VmStk, VmExe: Size of data, stack, and text segments.
   implicit none
   integer, intent(out) :: valueRSS
   character(len=200):: filename=' '
   character(len=80) :: line
   character(len=8)  :: pid_char=' '
   integer :: pid
   logical :: ifxst

   valueRSS=-1    ! return negative number if not found

   !--- get process ID

   pid=getpid()
   write(pid_char,'(I8)') pid
   filename='/proc/'//trim(adjustl(pid_char))//'/status'

   !--- read system file

   inquire (file=filename,exist=ifxst)
   if (.not.ifxst) then
     write (*,*) 'system file does not exist'
     return
   endif

   open(unit=100, file=filename, action='read')
   do
     read (100,'(a)',end=120) line
     if (line(1:6).eq.'VmRSS:') then
        read (line(7:),*) valueRSS
        exit
     endif
   enddo
   120 continue
   close(100)

   Write(2,*) 'pid=', pid,', mem usage=',valueRSS,' kB'

   return
end subroutine system_mem_usage
