!=========================
!Routines for the ATRID imager
Module AReadPKI
Contains
Subroutine ReadFlashImageDat(PeakStartChunk_ms, PeakWidth, StartTrace_ms, BoxCenter, WindowTime_ms)
   !
   use constants, only : dp,sample, c_mps, sample_ms
   use DataConstants, only : PeakNr_dim, ChunkNr_dim, DataFolder, EdgeOffset, Time_dim
   use Chunk_AntInfo, only : Station_nrMax, StartT_sam, TimeBase
   use Chunk_AntInfo, only : ChunkStTime_ms, ChunkFocus
   Use Interferom_Pars, only : N_Smth, BoundingBox
   use ThisSource, only : SourcePos,  TotPeakNr !NrP, t_ccorr,
   use ThisSource, only : PeakNrTotal !,PeakNr,  PlotCCPhase, Safety, Nr_corr
   use ThisSource, only : PeakPos, ChunkNr
   use ThisSource, only : Alloc_ThisSource
   use PredefinedTracks, only : PreDefTrackFile, t_track
   use PredefinedTracks, only : GetTrPos, PreDefTrackFile, t_track
   !use StationMnemonics, only : Statn_ID2Mnem, Station_Mnem2ID
   !#ifdef f2003
   !use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
   !                                          stdout=>output_unit, &
   !                                          stderr=>error_unit
   !#else
   !#define stdin  5
   !#define stdout 6
   !#define stderr 0
   !#endif   Implicit none
   !
   Implicit none
   Real(dp), allocatable, intent(out) :: PeakStartChunk_ms(:)
   Integer, allocatable, intent(out) :: PeakWidth(:)
   Real(dp), intent(out) :: StartTrace_ms, BoxCenter(1:3), WindowTime_ms
   integer :: i_chunk, i, k, i_nr=1 !, i_c ,i_ca  !,i_eoa
   Character(len=5) :: Station_Mnem, ExclStMnem(1:Station_nrMax)
   integer :: i_Peak, nxx, InUnit, PkW
   Integer, parameter :: stdin=5
   Integer, parameter :: SourceRead_N=4000
   Real(dp) :: Source_time(1:SourceRead_N), Source_pos(1:3,1:SourceRead_N)
   Integer :: Source_Width(1:SourceRead_N)
   Logical :: RS, RefSource(1:SourceRead_N)
   integer, parameter :: InLineLength=380
   character(len=InLineLength) :: lname
   character*50 :: FlashFile
   character*10 :: txt
   character*35 :: Frmt
   real(dp) :: x1,x2,x3, t_s, t_r, A, TimeShft
   Real(dp), external :: tShift_smpl, tShift_ms   !
   Integer :: FileLengthChar
   Logical :: OldDatFile, PlotFile, csvFile, ATRIDinp
   !
   i_chunk=1
   i_peak=0
   ChunkNr_dim=1
   Allocate(ChunkStTime_ms(1:ChunkNr_dim))  !  Not really used ???
   Allocate(ChunkFocus(1:3,1:ChunkNr_dim))  !  Not really used ???
   InUnit=17
   WindowTime_ms=-1.0
   Do while (InUnit .eq. 17)
      Call GetNonZeroLine(lname)
      write(2,"('Input line #',I2,3A)") i_nr,': "',trim(lname),'"'
      i_nr=i_nr+1
      Read(lname,*,iostat=nxx) FlashFile
      If(nxx.ne.0) exit    ! end of input reached
   !      inquire(FILE = 'files/'//TRIM(Simulation)//'_'//Station_name(i_file)//'.udt', exist=UnformattedSimul)
      Open(Unit=InUnit, STATUS='old',ACTION='read', FILE=TRIM(datafolder)//TRIM(FlashFile),iostat=nxx)
      If(nxx.ne.0) Then  ! failed to open the file, it does not exist.
         InUnit=stdin
         OldDatFile=.false.
         PlotFile=.false.
         ATRIDinp=.false.
         !  Check for Auto run of ATRID on first pass
         If(i_peak.eq.0) Then
            read(lname(2:InLineLength), *,iostat=nxx) k,t_r,x1,x2,x3  ! [--,tsource_ms,N_km,E_km,h_km]
            If(nxx.ne.0) Then    ! Try Auto Atrid input format
               Read(lname(2:InLineLength),*,iostat=nxx) StartTrace_ms, BoxCenter(1:3), WindowTime_ms  ! Start time



               If(nxx.ne.0) Then  ! no CenLoc was given, but probably a track-name
                  Read(lname(2:InLineLength),*,iostat=nxx) StartTrace_ms, PreDefTrackFile, t_track, WindowTime_ms
                  write(2,"(A)") 'Input line-1: "'//lname(1:1)//'|'//TRIM(lname(2:InLineLength))// &
                     '" !  Core/Source-| time, & position --or--   -| Midtrack_time, Track_file, track_time', &
                     'followed by WindowTime_ms'
                  write(2,*) 'Try track definition at input: '
                  Write(2,*) '*** StartTrace_ms, Track, t_track, WindowTime_ms but got a non-fitting format'
                  Write(2,*) '***   for t_track>0: TRI-D @ position=[track @ t=t_track, increment], t=[StartTrace, fixed]'
                  Write(2,*) '***   for t_track<0: TRI-D @ position=[track @ t=StartTrace], t=[StartTrace, increment]'
                  If(nxx.ne.0) Then
                     Write(2,*) '*** error, expected: StartTrace_ms, [BoxCenter(1:3) or track], WindowTime_ms '
                     Stop 'ATRID reading error'
                  EndIf
                  If( t_track.gt.0.0 ) Then  !  t_track was specified
                     Call  GetTrPos( t_track, BoxCenter(1:3))
                  Else
                     Call  GetTrPos( StartTrace_ms, BoxCenter(1:3))
                  EndIf
                  lname(1:1)='S'
               EndIf



               Call Convert2m(BoxCenter(1:3))
               Return
               !StartTrace_ms ! in source time
               !TraceLength ! in samples
               !BoxCenter(1:3) ! in [m]
            EndIf
         EndIf
         !
         Write(2,*) 'Not found!, no Sources-Data File: "',TRIM(datafolder)//TRIM(FlashFile),'", nxx=',nxx, &
            '. Assume lines from standard input where TimeBase is defined in Namelist'
         read(lname(2:InLineLength), *,iostat=nxx) k,t_r,x1,x2,x3, PkW ! [--,tsource_ms,N_km,E_km,h_km]
         If(nxx.eq.0) ATRIDinp=.true.  !  assume widths are given
         If((PkW.gt.-7) .and. (PkW.lt.0)) write(2,*) '*************** Ridiculously small fixed value for width:', PkW
         !Stop 'SourcesDataFile not found'
      Else           ! Read from file
         Write(2,*) 'Reading Sources-Data File: "',TRIM(datafolder)//TRIM(FlashFile),'"'
         FileLengthChar=Len_Trim(FlashFile)
         OldDatFile=.false.
         PlotFile=.false.
         If(FlashFile(FileLengthChar-3:FileLengthChar).eq. '.dat') Then
             OldDatFile=.true.
         ElseIf(FlashFile(FileLengthChar-3:FileLengthChar).eq. '.plt') Then
            PlotFile=.true.
         ElseIf(FlashFile(FileLengthChar-3:FileLengthChar).eq. '.csv') Then
            csvFile=.true.
            If(FlashFile(FileLengthChar-8:FileLengthChar-4).eq. 'ATRID') Then
               ATRIDinp=.true.
            Else
               ATRIDinp=.false.
            EndIf
         Else
            Write(2,*) 'Not a valid file type'
            exit
         EndIf
         write(2,*) FlashFile(FileLengthChar-8:FileLengthChar),', OldDatFile=',OldDatFile,', PlotFile=',PlotFile, csvFile
         ! file like  "T-Be.dat"  created by the TRI-D imager
         ! file like  "BdATRID.csv"  created by the ATRID imager in an earlier run
         Read(InUnit,"(A180)",iostat=nxx) lname
         If(nxx.ne.0) Then
            write(2,*) 'reading from an empty file? :"',trim(lname),'"'
            Stop 'Sources-Data File reading problem'
         EndIf
         Read(lname,*,iostat=nxx) BoundingBox(1:2,1:4),txt,TimeBase ! StartTime_ms
         If(nxx.eq.0)  Then      ! Reading from a plot-file
            write(2,*) 'Assume reading a .dat or .plt file, as generated by TRI-D or DataSelect.'
            Read(InUnit,*) txt
            Read(InUnit,*) txt
            Do i=1,9
               Read(InUnit,"(A180)",iostat=nxx) lname  !
               If(nxx.ne.0) exit
               read(lname, *,iostat=nxx) x2,x1,x3,k,PkW, i_chunk
               !write(2,*) 'line:',nxx, lname
               If(nxx.ne.0) exit   ! part of the frame for a 'SpecPowMx_d.dat'file
               !write(2,*) 'line:',nxx, lname
            EndDo
         Else
            If(csvFile) Then
               If(ATRIDinp) Then
                  Read(lname,*,iostat=nxx) TimeBase ! StartTime_ms
                  write(2,*) 'Assume reading a .csv Sources-Data File, generated by ATRID.'
                  Read(InUnit,*) txt
                  Read(InUnit,"(A180)",iostat=nxx) lname
               Else
                  Read(lname,*,iostat=nxx) x1,x2, TimeBase ! StartTime_ms
                  write(2,*) 'Assume reading a .csv Sources-Data File, generated by TRI-D.'
                  Read(InUnit,*) txt
                  Read(InUnit,"(A180)",iostat=nxx) lname
               EndIf
            Else
               write(2,*) 'Sources-Data File "',TRIM(datafolder)//TRIM(FlashFile),&
                  '", should have been a .csv (from ATRID or TRI-D) file.'
               stop 'Sources-Data File not compatible'
            EndIf
         endIf
      EndIf
      !
      !StartTime_ms=TimeBase
      !write(2,*) '!TimeBase=', TimeBase, ', Max nr of sources=', SourceRead_N, N_Smth, InUnitv
      !write(2,"(A, 8F11.3)") 'Bounding Box:', BoundingBox(1:2,1:4) ! 274.63
      !
      flush(unit=2)
      !
      !  'lname' needs to be loaded before this loop
      Do while (i_peak.lt.SourceRead_N)
         write(2,*) 'input line:"', trim(lname),'"'
         If(lname.eq.'') Then
            goto 1
         ElseIf(lname(1:1) .eq. '!') Then
            !write(2,*) '::', lname
            goto 1
         ElseIf(OldDatFile) Then
            RS =.false.  ! read-in time is really t_s, not. t_r
            read(lname, *,iostat=nxx) k,t_r,x2,x1,x3 ! [--,tsource_ms,E_km,N_km,h_km, reading a plot file]
            PkW = N_Smth
         ElseIf(PlotFile) Then
            RS =.false.  ! read-in time is really t_s, not. t_r
            read(lname, *,iostat=nxx) k,t_r,x1,x2,x3 ! [--,tsource_ms,E_km,N_km,h_km, reading a plot file]
            PkW = N_Smth
         ElseIf(lname(1:1) .eq. 'R') Then
            RS =.true.
            read(lname, *,iostat=nxx) txt, k,t_r,x1,x2,x3, PkW ! [--,tref_ms,N_km,E_km,h_km]
         Else
            RS =.false.  ! read-in time is really t_s, not. t_r
            If(TRIM(lname(2:InLineLength)).eq. '') goto 1
            If(ATRIDinp) Then
               read(lname(2:InLineLength), *,iostat=nxx) k,t_r,x1,x2,x3, PkW ! [--,tsource_ms,N_km,E_km,h_km]
            Else
               read(lname(2:InLineLength), *,iostat=nxx) k,t_r,x1,x2,x3 ! [--,tsource_ms,N_km,E_km,h_km]
               PkW = N_Smth
            EndIf
         EndIf
         !write(2,*) '!ReadFlashImageDat; txt, "', lname(1:1),'", "',TRIM(lname(2:InLineLength)),'"' ,RS, nxx, i_peak, k
         !flush(unit=2)
         If(nxx.ne.0) Then
            exit
         EndIf
         !
         !write(2,*) 'logicals:', OldDatFile, PlotFile, csvFile, ATRIDinp, PkW
         i_peak=i_peak+1
         If(PkW .eq. 0 ) PkW=N_Smth
         Source_Width(i_Peak)  = PkW
         Source_time(i_Peak) = t_r
         RefSource(i_Peak) = RS
         Source_pos(1,i_Peak)=x1
         Source_pos(2,i_Peak)=x2
         Source_pos(3,i_Peak)=x3
         !write(2,*) '!ReadFlashImageDat; txt, k,t_r,x1,x2,x3, PkW=',i_Peak, k,t_r,x1,x2,x3, PkW
         !flush(unit=2)
         Call Convert2m(Source_pos(1:3,i_Peak))
      1  continue
         lname=''
         Read(InUnit,"(A180)",iostat=nxx) lname  !
         If(nxx.ne.0) exit  ! last line
         !write(2,*) 'line:',nxx, lname
      EndDo
      If(InUnit .eq. 17) Close(Unit=17)  ! finished reading
   EndDo   ! Look for more input sources
   PeakNrTotal=i_peak
   PeakNr_dim=i_peak
   !write(2,*) '!ReadFlashImageDat; PeakNr_dim=', PeakNr_dim
   write(2,*) 'TimeBase=', TimeBase, ', total of sources=', PeakNr_dim,', N_smooth=', N_Smth
   flush(unit=2)
   !
   Allocate( PeakStartChunk_ms(1:PeakNr_dim), PeakWidth(1:PeakNr_dim) )
   Call Alloc_ThisSource    !  uses ChunkNr_dim & PeakNr_dim
   !
   PeakNrTotal = PeakNr_dim
   SourcePos(1:3,1:PeakNr_dim)=Source_Pos(1:3,1:PeakNr_dim)
   PeakWidth(1:PeakNr_dim) = Source_Width(1:PeakNr_dim)
   Do i_peak=1, PeakNrTotal
      TimeShft=tShift_ms(SourcePos(:,i_Peak))
      If(RefSource(i_Peak)) Then
         t_r=Source_time(i_Peak)
         t_s=t_r - TimeShft
      Else
         t_s=Source_time(i_Peak)
         t_r=t_s + TimeShft
      EndIf
      PeakStartChunk_ms(i_peak)=TimeBase + t_r - Time_dim*sample_ms/2.  ! in ms, the start of the chunk for this peak
      PeakPos(i_Peak) = Time_dim/2 ! NINT((TimeBase+t)/(Sample*1000.d0) - StartT_sam(ChunkNr(i_Peak))+ TimeShft) ! better be positive, <Time_dim
      !write(2,"(I4,A,I6, 2(A,F11.6),A,3F9.4, A, I4)") i_peak, ', PeakPos', PeakPos(i_Peak) &
      !   ,', t_source=', t_s, '[ms], t_ref=', t_r, '[ms], position=', &
      !   SourcePos(:,i_Peak)/1000., ', Window=', PeakWidth(i_Peak)
      !
   EndDo
   TotPeakNr(0,:)=0
   ChunkNr(1:PeakNr_dim)=1
   !
   !
   Return
End Subroutine ReadFlashImageDat
End Module AReadPKI
!====================================================================
!
!      Call AtridAuto(Chi2_Lim, Spread_Lim, BoxFineness_coarse, BoxSize_coarse, StartTrace_ms, TraceLength, BoxCenter, &
Module ATRIDA
Contains
Subroutine AtridAuto(Chi2_Lim, IVol_Lim, BoxFineness_coarse, BoxSize_coarse, StartTrace_ms, WindowTime_ms, BoxCenter, &
         PeakStartChunk_ms, SourceTime_ms,  &
         OrignlSourcePos, OrignlPeakpos, PeakWidth, Chi2, SourceIntensity, SourceIntSpread, SourceIntVol, I3, SourceStI, &
         SourceUn, SourceLin, SourceCirc, SourcePolMag, SourcePolZen, SourcePolAzi, SourcePoldOm, SourcesStk_NEh)
   use constants, only : dp, sample_ms
   use DataConstants, only : PeakNr_dim, ChunkNr_dim, RunMode, Time_Dim
   use DataConstants, only : DataFolder, OutFileLabel
   use Chunk_AntInfo, only : Unique_SAI, StartT_sam, Tot_UniqueAnt, Ant_Stations, TimeBase
   use Chunk_AntInfo, only : Unique_StatID,  Nr_UniqueStat, Nr_UniqueAnt, N_Chunk_max
   use Chunk_AntInfo, only : StartT_sam, TimeFrame, NoiseLevel, RefAnt, Simulation
   use DataConstants, only : Ant_nrMax
   use ThisSource, only : ChunkNr, SourcePos, PeakPos, Dual
   use ThisSource, only : PeakNr, PeakNrTotal, PeakChiSQ, TotPeakNr
   Use Interferom_Pars, only :  TIntens000
   Use Interferom_Pars, only : ChainRun
   Use Interferom_Pars, only :  N_Smth, N_fit, smooth, i_chunk, BoundingBox, TRIDFile
   Use Interferom_Pars, only : SumStrt, SumWindw, NrSlices, SliceLen
   Use Interferom_Pars, only :  SlicePnts, PulsePos, PulseWidth, MaxNrSlPnts
   use Interferom_Pars, only : Alloc_EInterfCalib_Pars
   use Interferom_Pars, only :  MaxSmPow, MaxSmPowLoc, NrPixSmPowTr
   Use Interferom_Pars, only : N_pix, d_loc, CenLocPol, CenLoc, PixLoc, IntFer_ant
   Use Interferom_Pars, only : N_best, Grd_best, Int_best
   Use Interferom_Pars, only : StI, StI12, StQ, StU, StV, StI3, StU1, StV1, StU2, StV2, P_un, P_lin, P_circ, Chi2pDF
   Use Interferom_Pars, only : PolZen, PolAzi, PolMag, PoldOm, Stk_NEh  !  calculated in "PolTestCath", called from EI_PolGridDel when (Outpt.ge.1); (PolMag(k), PolZen(k), PolAzi(k), PoldOm(k), k=1,3)
   use FFT, only : RFTransform_su,DAssignFFT, RFTransform_CF2CT
   use GLEplots, only : GLEplotControl
   use ThisSource, only : Alloc_ThisSource
   use HDF5_LOFAR_Read, only : CloseDataFiles
   use CPU_timeUsage, only : CPU_usage
   !
   Implicit none
   Real(dp), intent(inout) :: Chi2_Lim, IVol_Lim
   Real(dp), intent(in) :: BoxFineness_coarse, BoxSize_coarse
   Real(dp), intent(in) :: StartTrace_ms, WindowTime_ms ! in source time
   Integer :: TraceLength ! in samples
   Real(dp), intent(in) :: BoxCenter(1:3) ! in [m]
   real(dp), allocatable, intent(out) :: OrignlSourcePos(:,:), SourceIntensity(:), PeakStartChunk_ms(:)
   real(dp), allocatable, intent(out) :: SourceTime_ms(:), SourceIntSpread(:), SourceIntVol(:)
   real(dp), allocatable, intent(out) :: Chi2(:), I3(:), SourceStI(:), SourceUn(:), SourceLin(:), SourceCirc(:)
   real(dp), allocatable, intent(out) :: SourcePolZen(:,:), SourcePolAzi(:,:), SourcePolMag(:,:), SourcePoldOm(:,:)
   Integer, allocatable, intent(out) :: PeakWidth(:)
   Complex, allocatable, intent(out) :: SourcesStk_NEh(:,:,:)
   Integer, allocatable, intent(out) :: OrignlPeakpos(:)
   !
   Real(dp) :: SourceGuess(3,N_Chunk_max)
   integer :: j, k, nxx, DeadEnd=50, i_s1, i_s2 ! i,
   integer :: i_Peak, Sampl_shiftA, N_Smth_new !, N_slice_default
   Integer, parameter :: Nmin=5000
   Integer :: Npix(1:3,1:2), DataUnit,  N_minima       !  N_Int,
   Real(dp) :: dpix(1:3), Spread_Grd, Spread_Dist, Spread_Aspect, GridVolume(1:3), StartTref_ms, TimeShft, ST_ms  ! Spread_SclLong,
   Logical :: ContourPlot, SourceLabel_C(1:Nmin), verbose   ! , ListTable
   !Integer ::  SourcePos_C(1:Nmin), SourceWidth_C(1:Nmin)  !, SlicingPoints(0,3000)
   !character(len=2) :: txt
   Real(dp) :: BoxSize, BoxFineness !, Spread_IntDens
   Real(dp), external :: tShift_ms   !
   !
   !StartTrace_ms ! in source time
   !TraceLength ! in samples
   !BoxCenter(1:3) ! in [m]
   !
   write(2,"(A)") '===== Running Atrid in Auto-search mode with parameters:'
   write(2,"(A,F11.5)") 'middle time [ms] of the trace:',StartTrace_ms
   write(2,"(A,F8.6,A)") 'Length [ms] of the trace to be searched:',WindowTime_ms,', should not exeed 0.3 ms'
   write(2,"(A,3F9.4)") 'Center of the box (NEh, [km]) where search is performed:',BoxCenter(1:3)/1000.
   !
   TraceLength=WindowTime_ms/sample_ms
   If( (Time_dim-TraceLength-2*DeadEnd) .lt. 1500) Then
      write(2,*) 'too long window', (Time_dim-TraceLength)
      TraceLength=Time_dim-1500-2*DeadEnd
      write(2,"(A,I5,A,F8.6,A)") 'Search tracelength shortened to', TraceLength, ' samples=', TraceLength*sample_ms,' ms'
   EndIf
   TimeShft=tShift_ms(BoxCenter(1:3))
   StartTref_ms=TimeBase + StartTrace_ms + TimeShft - (Time_dim-TraceLength)*sample_ms/2.  ! in ms, the start of the chunk for this peak
   !
   NrSlices=1  !  depreciated, Dec 2021; used in making summed traces in 'EIAnalyzePixelTTrace' and in 'OutputIntfSlices'
   i_chunk=1
   StartT_sam(i_chunk)=StartTref_ms/sample_ms
   SumWindw=(Tracelength/2+DeadEnd)*2
   SumStrt=(Time_dim-SumWindw)/2  !offset from base of chunk
   CenLoc(1:3) = BoxCenter(1:3)
   Call RFTransform_su(Time_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   Call AntennaRead(i_chunk, BoxCenter)
   If(Simulation.eq."") CALL CloseDataFiles()    ! no more reading of antenna data
   call DAssignFFT()                   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !
   Call AntFieParGen()
   !
   If(NoiseLevel.lt.0) NoiseLevel=0.1
   !N_Smth=2  ! because summing starts from N_Smth/2 which should be >0
   ! Note: in 'Alloc_EInterfImag_Pars':    If(.not. SlicePnts)  NrPixSmPowTr=(SumWindw-1)/N_smth-1  ! smallest window has length=(2*N_smth +1)
   SlicePnts=.true.
   MaxNrSlPnts=Nmin
   If(Nmin.gt. MaxNrSlPnts) stop 'MaxNrSlPnts too small'
   Allocate( PulsePos(0:MaxNrSlPnts), PulseWidth(0:MaxNrSlPnts))  ! index 0 should cover the full range of all pulses
   NrPixSmPowTr=1
   PulsePos(1)=SumWindw/2
   PulseWidth(1)=SumWindw
   PulsePos(0)=PulsePos(1)
   PulseWidth(0)=PulseWidth(1)
   ContourPlot=.false.
   BoxSize=0.2
   BoxFineness=1.
   GridVolume(1:3)=BoxSize
   Call EI_run(GridVolume, BoxFineness)
!      write(2,*) 'dimensions for interferometry out of range: ', IntfBase, IntfBase+IntfDim, Time_dim
!      write(2,*) 'input for interferometry: ', SumStrt, SumWindw
   !
   !write(2,*) 'search window:', SumStrt,SumWindw, DeadEnd
   If(N_Smth.lt.7) Then
      Call WindowCutting(SumWindw, TIntens000(1:SumWindw), DeadEnd, Nmin, PulsePos, PulseWidth, SourceLabel_C, N_minima, &
         NoiseLevel )
      PeakNr_dim=N_minima
   Else
      Call WindowCutting_N_Smth(N_Smth, SumWindw, TIntens000(1:SumWindw), DeadEnd, Nmin, PulsePos, PulseWidth, PeakNr_dim, &
         NoiseLevel )
!Subroutine WindowCutting_N_Smth(N_Smth, N_Window, Intens, DeadEnd, Nmin, SourcePos_C, SourceWidth_C, N_Pks)
      SourceLabel_C(:)=.true.
   EndIf
   !
   PeakNrTotal=PeakNr_dim
   Call Alloc_ThisSource    !  uses ChunkNr_dim & PeakNr_dim
   !
   TotPeakNr(0,:)=0
   ChunkNr(1:PeakNr_dim)=1
   N_best=1000
   !N_slice_default=N_Smth
   !
   Call GLEplotControl(SpecialCmnd='echo on')  ! to generate unit=10 with the proper name
   !
   SlicePnts=.true.
   BoxSize=BoxSize_coarse
   BoxFineness=BoxFineness_coarse
   GridVolume(1:3)=BoxSize
   GridVolume(3)=2*BoxSize
   NrPixSmPowTr=PeakNr_dim
   Call CPU_usage
   flush(unit=2)
   Call EI_run(GridVolume, BoxFineness)
   !
   Write(2,"(A,3F6.2,A,3I3,A,3F6.1,A,F5.2,10(' -'))") 'Coarse-grid round finished on grd-D:' &
         ,d_loc(1:3),', grd-#:',N_pix(1:3,2), ', Box size: p/m',N_pix(1:3,2)*d_loc(1:3),' m, Fineness=',BoxFineness
   Allocate( PeakWidth(1:PeakNr_dim) )
   Call CPU_usage
   !
   i_peak=0
   write(2,"(A,F7.3,A)") 'Selecting peaks with max I_12 larger than ', NoiseLevel, ' @ coarse imagecube'
   Flush(unit=2)
   Do k=1,NrPixSmPowTr ! Now perform fine search and determine chi^2
      ST_ms= StartTref_ms + (SumStrt+PulsePos(k))*sample_ms - tShift_ms(MaxSmPowLoc(:,k)) - TimeBase
      If(MaxSmPow(k) .lt. 0.) Then ! otherwise border case or below noise
         write(2,"(A,I5,A,F11.6, A,F10.2)") ' pulse#',k,' @ t_source=', ST_ms, &
            '[ms] is borderline case or below noise level, I_12=',MaxSmPow(k)
      Else
      i_peak=i_peak+1
      SourcePos(1:3,i_Peak)=MaxSmPowLoc(1:3,k)
      PeakWidth(i_peak)=PulseWidth(k)
      Peakpos(i_Peak) = SumStrt+PulsePos(k)
      write(2,"(I4,A,I5,A,I6, A,F11.6, A,3F9.4, A,I4, A,F10.2)") i_peak, ' = pulse#',k,', sample=', PeakPos(i_Peak), &
         ', t_source=', ST_ms, '[ms], NEh=', SourcePos(:,i_Peak)/1000., ', Window=', PeakWidth(i_Peak), ', I_12=',MaxSmPow(k)
      EndIf
      Flush(unit=2)
   EndDo
   PeakNrTotal=i_peak
   !
   write(2,*) 'TimeBase=', TimeBase, ', total of ', PeakNrTotal,' sources kept out of', NrPixSmPowTr
   Call CPU_usage
   flush(unit=2)
   Allocate( PeakStartChunk_ms(1:PeakNrTotal) )
   Allocate( OrignlSourcePos(1:3,1:PeakNrTotal), SourceIntensity(1:PeakNrTotal) )
   Allocate( SourceTime_ms(1:PeakNrTotal), Chi2(1:PeakNrTotal), I3(1:PeakNrTotal) )
   Allocate( SourcePolZen(1:3,1:PeakNrTotal), SourcePolAzi(1:3,1:PeakNrTotal), SourcePolMag(1:3,1:PeakNrTotal), &
            SourcePoldOm(1:3,1:PeakNrTotal), SourceUn(1:PeakNrTotal), SourceLin(1:PeakNrTotal), SourceCirc(1:PeakNrTotal) )
   Allocate( OrignlPeakpos(1:PeakNrTotal), SourceStI(1:PeakNrTotal) ) !, SrcQual(1:10,1:PeakNrTotal) )
   Allocate( SourcesStk_NEh(1:3,1:3,1:PeakNrTotal), SourceIntSpread(1:PeakNrTotal), SourceIntVol(1:PeakNrTotal) )
   !
   Do i_Peak=1,PeakNrTotal
      OrignlSourcePos(1:3,i_Peak) = BoxCenter(1:3)
      OrignlPeakpos(i_peak) = Peakpos(i_peak)
      PeakStartChunk_ms(i_peak)= StartTref_ms
      SlicePnts=.false.
      N_Smth=abs(PeakWidth(i_peak))
      !write(txt,"(I2.2)") i_Peak
      write(*,*) 'Peak#', i_Peak
      !
      write(2,"(1x,A,i4,A,F11.6,A,2(F9.4,','),F9.4,A)") '++ ATRI-D, auto mode: Peak',i_peak,', t_Ref.ant=',&
         PeakStartChunk_ms(i_peak)+PeakPos(i_Peak)*sample_ms-TimeBase,'[ms], at (N,E,h)=(',SourcePos(:,i_Peak)/1000.,') [km]'
      !
      Call RefinementStep(BoxSize_coarse, BoxFineness_coarse, NoiseLevel, PeakNrTotal, PeakStartChunk_ms, i_Peak, Peakpos, &
         SourcePos, SourceIntensity, SourceTime_ms, SourceIntSpread, SourceIntVol, Chi2, I3, &
         SourceStI, SourceUn, SourceLin, SourceCirc, SourcePolMag, SourcePolZen, SourcePolAzi, SourcePoldOm, SourcesStk_NEh)
      !
      Write(2,*) 'Peak#',i_peak,  'finalized ================================================'
      Flush(unit=2)
      !Call GLEplotControl(SpecialCmnd='ls') !  ${FlashFolder}/
      !
   EndDo

End Subroutine AtridAuto
End Module ATRIDA
!====================================================================
Subroutine RefinementStep(BoxSize_coarse, BoxFineness_coarse, NoiseLevel, PeakNrTotal, PeakStartChunk_ms, i_Peak, Peakpos, &
   SourcePos, SourceIntensity, SourceTime_ms, SourceIntSpread, SourceIntVol, Chi2, I3, &
   SourceStI, SourceUn, SourceLin, SourceCirc, SourcePolMag, SourcePolZen, SourcePolAzi, SourcePoldOm, SourcesStk_NEh)
   use constants, only : dp, sample_ms
   !use DataConstants, only : DataFolder, OutFileLabel, Time_dim
   !use GLEplots, only : GLEplotControl
   use Chunk_AntInfo, only : TimeBase
   use Interferom_Pars, only :  MaxSmPow, MaxSmPowLoc
   Use Interferom_Pars, only : N_best, Grd_best, Int_best
   Use Interferom_Pars, only : StI, StI12, StQ, StU, StV, StI3, StU1, StV1, StU2, StV2, P_un, P_lin, P_circ, Chi2pDF
   Use Interferom_Pars, only : PolZen, PolAzi, PolMag, PoldOm, Stk_NEh  !  calculated in "PolTestCath", called from EI_PolGridDel when (Outpt.ge.1); (PolMag(k), PolZen(k), PolAzi(k), PoldOm(k), k=1,3)
   Use Interferom_Pars, only : ChainRun
   Implicit none
   Integer, intent(in) :: PeakNrTotal,  Peakpos(1:PeakNrTotal)
   Real(dp), intent(in) :: BoxSize_coarse, BoxFineness_coarse, NoiseLevel
   Real(dp), intent(in) :: PeakStartChunk_ms(1:PeakNrTotal)
   Real(dp), intent(out) :: SourcePos(1:3,1:PeakNrTotal), Chi2(1:PeakNrTotal), SourceIntensity(1:PeakNrTotal), &
      SourceIntSpread(1:PeakNrTotal), SourceIntVol(1:PeakNrTotal), I3(1:PeakNrTotal), SourceStI(1:PeakNrTotal), &
      SourceTime_ms(1:PeakNrTotal), SourceUn(1:PeakNrTotal), SourceLin(1:PeakNrTotal), SourceCirc(1:PeakNrTotal), &
      SourcePolMag(1:3,1:PeakNrTotal), SourcePolZen(1:3,1:PeakNrTotal), SourcePolAzi(1:3,1:PeakNrTotal), &
      SourcePoldOm(1:3,1:PeakNrTotal)
   Complex, intent(out) :: SourcesStk_NEh(1:3,1:3,1:PeakNrTotal)
   Integer :: i_peak, Npix(1:3,1:2), N_Int
   Real(dp) :: dpix(1:3), Spread_Grd, Spread_Dist, Spread_SclLong, Spread_Aspect, Spread_IntVol
   !Integer ::    !, k, Sam_peak, DataUnit, DataUnitp, Width_max, keep
   Logical :: ContourPlot, verbose   ! , ListTable
   Real(dp) :: Boxsize, BoxFineness   !
    Real(dp), external :: tShift_ms   !
  !
      BoxSize=(10.+BoxSize_coarse/3.)/2.
      BoxFineness=BoxFineness_coarse/3.
      If(PeakNrTotal.lt.15 .and. ChainRun.eq.0) Then
         ContourPlot=.true.
         !verbose=.true.
         verbose=.false.
      Else
         ContourPlot=.false.
         verbose=.false.
      EndIf
      Call PeakInterferometerRun(PeakStartChunk_ms, i_Peak, BoxSize, BoxFineness, dpix, Npix, ContourPlot )
      write(2,"('Peak#',I4, ', Intensity=',F11.5, ', (N,E,h)=',3F11.5,'[km], diffs:', 3F8.3,'[m]')") i_peak, &
         MaxSmPow(1), MaxSmPowLoc(1:3,1)/1000., MaxSmPowLoc(1:3,1)-SourcePos(1:3,i_Peak)
      !
      !      i_s1=1+PeakWidth(i_peak)/2   ! if width=odd cut is tampered by smooth
      !      i_s2=1+3*PeakWidth(i_peak)/2
      !write(2,*) 'search window:', SumStrt,SumWindw,PeakWidth(i_peak), i_s1, i_s2
      !write(2,*) 'lo:', TIntens000(i_s1-2:i_s1+2)
      !write(2,*) 'hi:',TIntens000(i_s2-2:i_s2+2)
      If(MaxSmPowLoc(3,1) .gt. NoiseLevel) Then ! otherwise border case
         SourcePos(1:3,i_Peak)=MaxSmPowLoc(1:3,1)
      Else
         write(2,*) '********************** no change in source position ********************************'
      EndIf
      SourceIntensity(i_Peak)=MaxSmPow(1)
      SourceTime_ms(i_Peak)= PeakStartChunk_ms(i_peak)+PeakPos(i_Peak)*sample_ms - tShift_ms(SourcePos(:,i_Peak)) - TimeBase
      Write(2,"(A,I4,A,F11.6,A,3F6.2,A,3I3,A,3F6.1,A,F5.2,10(' -'))") 'Peak#',i_peak, ', @',SourceTime_ms(i_Peak), &
            ' [ms], Fine-grid round finished on grd-D:' &
            ,dpix(1:3),', grd-N:',Npix(1:3,2), ', Box size: p/m',Npix(1:3,2)*dpix(1:3),' m, Fineness=',BoxFineness
      Call IntensitySpread(Grd_best(1,1,1), Int_best(1,1), dpix, N_best, N_Int, Npix, BoxFineness, &
         Spread_Grd, Spread_Dist, Spread_SclLong, Spread_Aspect, Spread_IntVol, verbose)
      !   Spread_Grd, Spread_Dist, Spread_SclLong, Spread_Aspect, Spread_IntVol, verbose)
      !write(2,*) '!Source Spread:', N_Int, Spread_Grd, Spread_Dist, Int_best(N_Int,1)/Int_best(1,1)
      !Spread_SclLong=Spread_SclLong*BoxFineness  ! to make it independent of the gridding coarseness
      SourceIntSpread(i_Peak)=Spread_SclLong
      SourceIntVol(i_Peak)=Spread_IntVol  ! Spread_Aspect*100. ! in fact inverse intensity density
      !write(2,*) 'Long spreading axis, normalized for Fineness and Intensity cut, =',Spread_SclLong, &
      !
      ! prepare & produce curtain plot
      ! Calculate polarization observables at interpolated position
      Call PeakCurtain(i_Peak)
!      SrcQual(6,i_Peak)=NINT(Chi2pDF*10.)
      Chi2(i_Peak) = Chi2pDF
      I3(i_Peak) = StI3/StI
      SourceStI(i_Peak) = StI
      SourceUn(i_Peak) = P_un
      SourceLin(i_Peak) = P_lin
      SourceCirc(i_Peak) = P_circ
      SourcePolMag(1:3,i_Peak) = PolMag(1:3)
      SourcePolZen(1:3,i_Peak) = PolZen(1:3)
      SourcePolAzi(1:3,i_Peak) = PolAzi(1:3)
      SourcePoldOm(1:3,i_Peak) = PoldOm(1:3)      !
      SourcesStk_NEh(1:3,1:3,i_Peak) = Stk_NEh(1:3,1:3)
End Subroutine RefinementStep
!=================================
Subroutine PeakInterferoOption(Chi2_Lim, IVol_Lim, BoxFineness_coarse, BoxSize_coarse) !CurtainHalfWidth
!   Adaptive TRI-D (ATRI-D) imaging of peaks found by the impulsive (or the TRI-D) imager
!
!   Procedural steps:
!   1) Read sources selected by FlashImage
!   2) Use TRID at the central location to Adept the window from the coherent time trace
!   3) Use TRI-D twice with different Fineness factors.
!
!     -
   use constants, only : dp, sample_ms
   use DataConstants, only : PeakNr_dim, ChunkNr_dim, RunMode, Time_Dim
   use DataConstants, only : DataFolder, OutFileLabel
   use Chunk_AntInfo, only : Unique_SAI, StartT_sam, Tot_UniqueAnt, Ant_Stations, TimeBase
   use Chunk_AntInfo, only : Unique_StatID, Nr_UniqueStat, Nr_UniqueAnt, N_Chunk_max, NoiseLevel
   use DataConstants, only : Ant_nrMax
   use ThisSource, only : ChunkNr, SourcePos, PeakPos, Dual
   use ThisSource, only : PeakNr, PeakNrTotal, PeakChiSQ
   Use Interferom_Pars, only :  TIntens000
   Use Interferom_Pars, only : ChainRun, CenLoc
   Use Interferom_Pars, only :  N_Smth, N_fit, smooth, i_chunk, BoundingBox, TRIDFile
   Use Interferom_Pars, only : IntFer_ant,  Nr_IntferCh ! the latter gives # per chunk
   use Interferom_Pars, only : Alloc_EInterfCalib_Pars
   use Interferom_Pars, only :  MaxSmPow, MaxSmPowLoc
   Use Interferom_Pars, only : N_best, Grd_best, Int_best
   Use Interferom_Pars, only : StI, StI12, StQ, StU, StV, StI3, StU1, StV1, StU2, StV2, P_un, P_lin, P_circ, Chi2pDF
   Use Interferom_Pars, only : PolZen, PolAzi, PolMag, PoldOm, Stk_NEh  !  calculated in "PolTestCath", called from EI_PolGridDel when (Outpt.ge.1); (PolMag(k), PolZen(k), PolAzi(k), PoldOm(k), k=1,3)
   !use StationMnemonics, only : Statn_ID2Mnem !, Station_Mnem2ID
   use FFT, only : RFTransform_su,DAssignFFT, RFTransform_CF2CT
   use GLEplots, only : GLEplotControl
   use AReadPKI, only : ReadFlashImageDat
   use ATRIDA, only : AtridAuto
   use CPU_timeUsage, only : CPU_usage
   !
   Implicit none
   Real(dp), intent(inout) :: Chi2_Lim, IVol_Lim
   Real(dp), intent(in) :: BoxFineness_coarse, BoxSize_coarse
   !
   Real(dp) :: SourceGuess(3,N_Chunk_max), StartTrace_ms, BoxCenter(1:3), WindowTime_ms
   integer :: j, k, nxx ! i,
   integer :: i_Peak, Sampl_shiftA, N_Smth_new !, N_slice_default
   Integer :: Npix(1:3,1:2), N_Int, DataUnit, N_minima
   Real(dp) :: MinSrcDist, dpix(1:3)
   real(dp), allocatable :: OrignlSourcePos(:,:), SourceIntensity(:), PeakStartChunk_ms(:)
   real(dp), allocatable :: SourceTime_ms(:), SourceIntSpread(:), SourceIntVol(:) !, MovingSys(:,:)
   real(dp), allocatable :: Chi2(:), I3(:), SourceStI(:), SourceUn(:), SourceLin(:), SourceCirc(:)
   real(dp), allocatable :: SourcePolZen(:,:), SourcePolAzi(:,:), SourcePolMag(:,:), SourcePoldOm(:,:)
   Integer, allocatable :: PeakWidth(:)
   Integer, allocatable :: OrignlPeakpos(:) !, SrcQual(:,:)
   Complex, allocatable :: SourcesStk_NEh(:,:,:)
   Logical :: ContourPlot, SourceLabel_C(0:20), verbose   ! , ListTable
   Integer ::  SourcePos_C(0:20), SourceWidth_C(0:20), Nmin
   Character(len=180):: Message
   Real(dp) :: time, x_plot, y_plot, z_plot, Cnst, BoxSize, BoxFineness, Chi2Mean, Chi2Var
   Real(dp) :: Spread_Grd, Spread_Dist, Spread_SclLong, Spread_Aspect, Spread_IntVol
   Real(dp), external :: tShift_ms   !
   !
   !Call GetNonZeroLine(Message)
   !Read(Message,*,iostat=nxx) Chi2_Lim, Spread_Lim, BoxFineness_coarse, BoxSize_coarse !, Track ! StartTime_ms
   !If(nxx.ne.0) Then
   !   Write(2,*) 'Read: "', Trim(Message),'" but expected:'
   !   write(2,*) 'Extra input parameters, Chi2_Lim=', Chi2_Lim, ', Spread_Lim=', Spread_Lim, &
   !      ', BoxFineness_coarse=', BoxFineness_coarse, ', BoxSize_coarse=',BoxSize_coarse,'[m], cubic.'
   !   stop 'reading error'
   !EndIf
   write(2,"(A,F5.2, A,F6.1, A, F5.2, A,f5.1,A)") 'Extra input parameters, Chi2_Lim=', Chi2_Lim, &
         ', IVol_Lim=', IVol_Lim, &
         ', BoxFineness_coarse=', BoxFineness_coarse, ', BoxSize_coarse=',BoxSize_coarse,'[m], cubic.'
   !
   Dual=.true.
   TRIDFile=.false.  !  do not produce files when running the trid imager
   Call ReadFlashImageDat( PeakStartChunk_ms, PeakWidth, StartTrace_ms, BoxCenter, WindowTime_ms)
   If(WindowTime_ms.gt. 0.) Then ! run in Auto-Atrid mode
      Call AtridAuto(Chi2_Lim, IVol_Lim, BoxFineness_coarse, BoxSize_coarse, StartTrace_ms, WindowTime_ms, BoxCenter, &
         PeakStartChunk_ms, SourceTime_ms,  &
         OrignlSourcePos, OrignlPeakpos, PeakWidth, Chi2, SourceIntensity, SourceIntSpread, SourceIntVol, I3, SourceStI, &
         SourceUn, SourceLin, SourceCirc, SourcePolMag, SourcePolZen, SourcePolAzi, SourcePoldOm, SourcesStk_NEh)
      goto 9
   EndIf
   flush(unit=2)
   Allocate( OrignlSourcePos(1:3,1:PeakNr_dim), SourceIntensity(1:PeakNr_dim) )
   Allocate( SourceTime_ms(1:PeakNr_dim), Chi2(1:PeakNr_dim), I3(1:PeakNr_dim) )
   Allocate( SourcePolZen(1:3,1:PeakNr_dim), SourcePolAzi(1:3,1:PeakNr_dim), SourcePolMag(1:3,1:PeakNr_dim), &
            SourcePoldOm(1:3,1:PeakNr_dim), SourceUn(1:PeakNr_dim), SourceLin(1:PeakNr_dim), SourceCirc(1:PeakNr_dim) )
   Allocate( OrignlPeakpos(1:PeakNr_dim), SourceStI(1:PeakNr_dim) ) !, SrcQual(1:10,1:PeakNr_dim) )
   Allocate( SourcesStk_NEh(1:3,1:3,1:PeakNrTotal), SourceIntSpread(1:PeakNr_dim), SourceIntVol(1:PeakNr_dim) )
   !
   Call AntFieParGen()
   !
   N_best=3000
   !N_slice_default=N_Smth
   !
   Call GLEplotControl(SpecialCmnd='echo on')  ! to generate unit=10 with the proper name
   !DataUnit=28
   !OPEN(UNIT=DataUnit,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'IntfPkSources.dat')
   !Write(DataUnit,*) PeakNrTotal,TRIM(OutFileLabel), 10.,' 0 0 0 0 0 0 0 ',0.1 ,  '1.0 !' ! gleRead:  NTracks EventNr Q t_start label$ AmplitudePlot a1 b1 c1 a2 b2 c2 d2
   write(2,*) 'Starting source analysis for',PeakNrTotal,' sources'
   Call CPU_usage
   flush(unit=2)
   Do i_Peak=1,PeakNrTotal
!      SrcQual(10,i_Peak)=0
!      SrcQual(5,i_Peak) =0
      OrignlSourcePos(1:3,i_Peak) = SourcePos(1:3,i_Peak)
      OrignlPeakpos(i_peak) = Peakpos(i_peak)
      N_Smth=abs(PeakWidth(i_peak))
      !write(txt,"(I2.2)") i_Peak
      write(*,*) 'Peak#', i_Peak
      !
      write(2,"(1x,A,i4,A,F11.6,A,2(F9.4,','),F9.4,A)") '++ ATRID: Peak',i_peak,', t_Ref.ant=',&
         PeakStartChunk_ms(i_peak)+PeakPos(i_Peak)*sample_ms-TimeBase,'[ms], at (N,E,h)=(',SourcePos(:,i_Peak)/1000.,') [km]'
      !
      !dpix(1:3) = (/ 4.d0, 3.d0, 6.d0 /)
      !Npix(1:3,2)=(/ 0, 0, 0 /)
      !Npix(1:3,1)=-Npix(1:3,2)
      ContourPlot=.false.
      BoxSize=0.2
      BoxFineness=1.
      Call PeakInterferometerRun(PeakStartChunk_ms, i_Peak, BoxSize, BoxFineness, dpix, Npix, ContourPlot )
      !
      !
      !write(2,*) 'TIntens000(1:2*N_smth+1)', TIntens000(1:2*N_smth+1)
      !Call WindowShift(N_Smth, smooth, TIntens000(1:2*N_smth+1), Sampl_shiftA, txt)
      Call WindowOptimal(N_Smth, smooth, TIntens000(1:2*N_smth+1), N_Smth_new, Sampl_shiftA, PeakWidth(i_peak), NoiseLevel)
      Write(2,"(A,I4,A,I3,A,I4,A,I4,A)") 'Peak#',i_peak, ', time shift:',Sampl_shiftA,', window-length changed from ', &
         N_Smth,' to ',N_Smth_new, ' ---------------------------------'
      N_Smth=N_Smth_new
      Peakpos(i_Peak) = Peakpos(i_Peak)+ Sampl_shiftA
      If(PeakWidth(i_peak) .gt. 0) PeakWidth(i_Peak) = N_Smth  ! don't change in case value was set to a negative value [window=abs()]
!      If(abs(Sampl_shiftA) .eq. N_Smth )  SrcQual(5,i_Peak)=1
      !
      BoxSize=BoxSize_coarse
      BoxFineness=BoxFineness_coarse
      If(PeakNrTotal.lt.15 ) Then
         ContourPlot=.true.
         !verbose=.true.
         verbose=.false.
      Else
         ContourPlot=.false.
         verbose=.false.
      EndIf
      Call PeakInterferometerRun(PeakStartChunk_ms, i_Peak, BoxSize, BoxFineness, dpix, Npix, ContourPlot )
      write(2,"('Peak#',I4, ', Intensity=',F11.5, ', (N,E,h)=',3F11.5,'[km], diffs:', 3F8.3,'[m]')") i_peak, &
         MaxSmPow(1), MaxSmPowLoc(1:3,1)/1000., MaxSmPowLoc(1:3,1)-SourcePos(1:3,i_Peak)
      Call IntensitySpread(Grd_best(1,1,1), Int_best(1,1), dpix, N_best, N_Int, Npix, BoxFineness, &
         Spread_Grd, Spread_Dist, Spread_SclLong, Spread_Aspect, Spread_IntVol, verbose)
         !
      SourcePos(1:3,i_Peak)=MaxSmPowLoc(1:3,1)
      If(MaxSmPow(1) .lt. 0.0) Then ! otherwise border case
         write(2,*) '********************** Borderline case, change in source position ********************************'
      EndIf
      Write(2,"(A,I4,A,3F6.2,A,3I3,A,3F6.1,A,F5.2,10(' -'))") 'Peak#',i_peak, ', Intensity round 1 finished on grd-D:' &
            ,dpix(1:3),'m , grd-N:',Npix(1:3,2), ', Box size: p/m',Npix(1:3,2)*dpix(1:3),' m, Fineness=',BoxFineness
      !write(2,*) MaxSmPow(1),' ; ', MaxSmPowLoc(1:3,1)
      !
      Call RefinementStep(BoxSize_coarse, BoxFineness_coarse, NoiseLevel, PeakNrTotal, PeakStartChunk_ms, i_Peak, Peakpos, &
         SourcePos, SourceIntensity, SourceTime_ms, SourceIntSpread, SourceIntVol, Chi2, I3, &
         SourceStI, SourceUn, SourceLin, SourceCirc, SourcePolMag, SourcePolZen, SourcePolAzi, SourcePoldOm, SourcesStk_NEh)
      !
      Write(2,*) 'Peak#',i_peak,  'finalized ================================================'
      Flush(unit=2)
      !Call GLEplotControl(SpecialCmnd='ls') !  ${FlashFolder}/
      !
   EndDo
   !Close(Unit=DataUnit)
9  continue
   !
   Call CPU_usage
   MinSrcDist=BoxSize_coarse
   Call OutputPkInterfer(PeakNrTotal, Chi2_Lim, IVol_Lim, MinSrcDist, NoiseLevel, TimeBase, PeakStartChunk_ms, SourceTime_ms, &
         Peakpos, SourcePos,OrignlSourcePos, OrignlPeakpos, PeakWidth, Chi2, SourceIntensity, SourceIntSpread, SourceIntVol, I3,&
         SourceStI, SourceUn, SourceLin, SourceCirc, SourcePolMag, SourcePolZen, SourcePolAzi, SourcePoldOm,  SourcesStk_NEh )
   !
   Call GLEplotControl(CleanPlotFile=TRIM(OutFileLabel)//'Interferometer*.z')
   Call GLEplotControl(CleanPlotFile=TRIM(OutFileLabel)//'IntfSpecPowMx_d.dat')
   Call GLEplotControl(CleanPlotFile=TRIM(OutFileLabel)//'_EISpec*.dat')
   Call GLEplotControl(CleanPlotFile=TRIM(OutFileLabel)//'IntfSpecWin_d.csv')
   !
   Call GLEplotControl( Submit=.true.)
   !
   If(ChainRun.ne.0) Call ChainRuns(CenLoc, WindowTime_ms)
   Return
   !
   !====================================================================
   !
End Subroutine PeakInterferoOption
!-----------------------------------------
Subroutine OutputPkInterfer(PeakNrTotal, Chi2_set, IVol_set, MinSrcDist, NoiseLevel, TimeBase, &
         PeakStartChunk_ms, SourceTime_ms, Peakpos, SourcePos, &
         OrignlSourcePos, OrignlPeakpos, PeakWidth, Chi2, SourceIntensity, SourceIntSpread, SourceIntVol, I3, SourceStI, &
         SourceUn, SourceLin, SourceCirc, SourcePolMag, SourcePolZen, SourcePolAzi, SourcePoldOm,  Stk_NEh )
   use constants, only : dp, pi, sample_ms
   use DataConstants, only : DataFolder, OutFileLabel, FlashName, Time_dim
   use GLEplots, only : GLEplotControl
   Implicit none
   Integer, intent(in) :: PeakNrTotal,  PeakWidth(1:PeakNrTotal), Peakpos(1:PeakNrTotal), OrignlpeakPos(1:PeakNrTotal)
   Real(dp), intent(in) :: Chi2_set, IVol_set
   Real(dp), intent(in) :: MinSrcDist, NoiseLevel, TimeBase, PeakStartChunk_ms(1:PeakNrTotal), SourceTime_ms(1:PeakNrTotal), &
      SourcePos(1:3,1:PeakNrTotal), Chi2(1:PeakNrTotal), SourceIntensity(1:PeakNrTotal), SourceIntSpread(1:PeakNrTotal), &
      SourceIntVol(1:PeakNrTotal), I3(1:PeakNrTotal), SourceStI(1:PeakNrTotal), &
      SourceUn(1:PeakNrTotal), SourceLin(1:PeakNrTotal), SourceCirc(1:PeakNrTotal), &
      SourcePolMag(1:3,1:PeakNrTotal), SourcePolZen(1:3,1:PeakNrTotal), SourcePolAzi(1:3,1:PeakNrTotal), &
      SourcePoldOm(1:3,1:PeakNrTotal), OrignlSourcePos(1:3,1:PeakNrTotal)
   Complex, intent(in) :: Stk_NEh(1:3,1:3,1:PeakNrTotal)
   Integer :: i_peak, k, Sam_peak, DataUnit, DataUnitp, Width_max, keep
   Real(dp) :: I20, StI, Chi2Mean, Chi2Var, SpreadMean, SpreadVar, Lin1, PolN, PolE, Polh, PolPlotScale, I3_lim  ! , Aspct_lim
   Real(dp) :: IVolMean, IVolVar, Spread_lim, Chi2_Lim, IVol_Lim  ! , Aspct_lim
   Real(dp) :: tMin, tMax, BBx_min(1:3), BBx_max(1:3)
   character(len=2) :: txt
   character(len=4) :: txt2
   character(len=8) :: txt3
   Character*30 :: Qual
   Character(len=180):: Message
   Logical :: DeSelect(1:PeakNrTotal)
   !
   DataUnit=27  ! for storing more complete info on analyzed sources
   DataUnit=28  ! for storing information for plots
   write(2,*) 'starting Output for the ATRID Imager for',PeakNrTotal, ' sources'
   If(PeakNrTotal.eq.0) return
   !
   Chi2Mean=0.
   Chi2Var=0.
   SpreadMean=0.
   SpreadVar=0.
   IVolMean=0.
   IVolVar=0.
   Width_max=0.
   Do i_Peak=1,PeakNrTotal
      write(2,"('Peak#',I4, A,F11.6,A,F11.6, A,F6.1, ', (N,E,h)=',3F11.5,'[km], diffs:', I3, 3F8.3,'[m]')") &
         i_peak, ', t_src=', SourceTime_ms(i_Peak), ', t_ref=', PeakStartChunk_ms(i_peak)+PeakPos(i_Peak)*sample_ms -TimeBase, &
         '[ms],  Intensity=',SourceIntensity(i_Peak), SourcePos(1:3,i_Peak)/1000., Peakpos(i_Peak)-OrignlPeakpos(i_Peak), &
         SourcePos(1:3,i_Peak)-OrignlSourcePos(1:3,i_Peak)
      Chi2Mean=Chi2Mean + Chi2(i_Peak)
      Chi2Var=Chi2Var + Chi2(i_Peak)*Chi2(i_Peak)
      SpreadMean=SpreadMean + SourceIntSpread(i_Peak)
      SpreadVar=SpreadVar + SourceIntSpread(i_Peak)*SourceIntSpread(i_Peak)
      IVolMean=IVolMean + SourceIntVol(i_Peak)
      IVolVar=IVolVar + SourceIntVol(i_Peak)*SourceIntVol(i_Peak)
      If(PeakWidth(i_Peak).gt. Width_max) Width_max= PeakWidth(i_Peak)
   EndDo
   Chi2Mean=Chi2Mean/PeakNrTotal
   Chi2Var=sqrt((Chi2Var - Chi2Mean*Chi2Mean*PeakNrTotal)/(PeakNrTotal-1))
   If(Chi2_set.lt.0) Then
      !write(2,*) 'The input negative Chi2_Lim is set equal to Chi2_Lim=Mean+Sigma=',Chi2Mean+Chi2Var
      Chi2_Lim=Chi2Mean+Chi2Var
   Else
      Chi2_Lim=Chi2_set  ! Not to spoil limit setting when using ChainRun
   EndIf
   SpreadMean=SpreadMean/PeakNrTotal
   SpreadVar=sqrt((SpreadVar - SpreadMean*SpreadMean*PeakNrTotal)/(PeakNrTotal-1))
   IVolMean=IVolMean/PeakNrTotal
   IVolVar=sqrt((IVolVar - IVolMean*IVolMean*PeakNrTotal)/(PeakNrTotal-1))
!   IVol_lim=10.
   If(IVol_set.lt.0) Then
      !write(2,*) 'The input negative Chi2_Lim is set equal to Chi2_Lim=Mean+Sigma=',Chi2Mean+Chi2Var
      IVol_Lim=IVolMean+IVolVar
   Else
      IVol_Lim=IVol_set
   EndIf
   If(IVol_Lim.gt.888) Then
      !write(2,*) 'The input negative Chi2_Lim is set equal to Chi2_Lim=Mean+Sigma=',Chi2Mean+Chi2Var
      write(2,*) 'IVol limit tuned down from ', IVol_Lim,' to 888. which is still ridiculously large.'
      IVol_Lim=888 ! when larger than 999.5 problems with the format of the .plt file
   EndIf
   !If(Spread_Lim.lt.0) Then
   !   !write(2,*) 'The input negative Chi2_Lim is set equal to Chi2_Lim=Mean+Sigma=',Chi2Mean+Chi2Var
   !   Spread_Lim=SpreadMean+SpreadVar
   !EndIf
   write(2,"(A, 2F7.2,A,F7.2)") 'chi^2 mean and st.dev=', Chi2Mean, Chi2Var,'; Chi2_Lim =',Chi2_Lim
   write(2,"(A, 2F7.2,A,F7.1)") 'Long-Spread Intensity: mean and st.dev=', SpreadMean, SpreadVar,'; Spread_Lim =',Spread_Lim
   write(2,"(A, 2F7.1,A,F7.1)") 'Intensity Volume: mean and st.dev=', IVolMean, IVolVar,'; IVol_Lim =',IVol_Lim
   write(2,"(A, I4)") 'Maximum Width=', Width_max
   !
   !     Disqualify sources that do not obey the goodness criteria ===================================
   DeSelect(:)=.false.
   I3_lim=0.9
   !Aspct_lim=5.
   Do i_Peak=1,PeakNrTotal  ! Hard deselect criteria:
      If(Chi2(i_Peak) .gt. Chi2_lim) DeSelect(i_Peak)=.true.
      !If(SourceIntSpread(i_Peak) .gt. Spread_Lim)DeSelect(i_Peak)=.true.
      If(SourceIntensity(i_Peak) .lt. NoiseLevel) DeSelect(i_Peak)=.true.
      !If(I3(i_Peak) .gt. I3_lim) DeSelect(i_Peak)=.true.
      If(SourceIntVol(i_Peak) .gt. IVol_lim) DeSelect(i_Peak)=.true.
   EndDo
   !
   ! Determine bounding box for all sources
   BBx_min(1:3)=+1.e7
   BBx_max(1:3)=-1.e7
   tmin=SourceTime_ms(1)
   tmax=SourceTime_ms(1)
   keep=0
   Do i_Peak=1,PeakNrTotal
      If(DeSelect(i_Peak)) cycle
      keep=keep+1
      If(SourceTime_ms(i_Peak) .gt. tmax) tmax=SourceTime_ms(i_Peak)
      If(SourceTime_ms(i_Peak) .lt. tmin) tmin=SourceTime_ms(i_Peak)
      Do k=1,3
         If(SourcePos(k,i_Peak).gt.BBx_max(k)) BBx_max(k)=SourcePos(k,i_Peak)
         If(SourcePos(k,i_Peak).lt.BBx_min(k)) BBx_min(k)=SourcePos(k,i_Peak)
      EndDo
   EndDo
   tmax= ( NINT(tmax*1000. + 0.5) )/1000. ! integer mu second units [ms]
   tmin= ( NINT(tmin*1000. - 0.5) )/1000.
   Do k=1,3
      !write(2,*) 'max,min,k=',k, BBx_max(k), BBx_min(k)
      BBx_max(k)=NINT(BBx_max(k) + 1.)/1000.  ! 1 m units and convert to [km]; 6F9.3
      BBx_min(k)=NINT(BBx_min(k) - 1.)/1000.  ! could have been +/- 0.5
      !write(2,*) 'max,min,k=',k, BBx_max(k), BBx_min(k)
   Enddo
   !
   OPEN(UNIT=DataUnit,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'ATRID.csv') ! comma separated value
   Write(DataUnit,*) TimeBase
!  #, Source t[ms] ;    North      East        h [km]   Width  Chi^2    St_I/20     I_3%   P_un, P_lin, P_circ%  Zen , azi    ;  Stk_NEh (N,N) (N,E) (N,h) (E,E) (E,h) (h,h)
   Write(DataUnit,"(A,T25,A,T36,A,T48,A,T57,A, T64,A,T72,A,T81,A,T93,A,  T100,A , T122,A,   T134,A,A)") '!  #, Source t[ms] ;', &
      'North','East', 'h [km]','Width','Chi^2','IntVol','St_I/20',  'I_3%','P_un, P_lin, P_circ%','Zen , azi ', &
      ' ;  Stk_NEh (N,N) (N,E) (N,h) (E,E) (E,h) (h,h)'
   !
   write(Qual,"(A,F5.2,'ns & I_V<',F8.1,A)") '"Q<',Chi2_Lim, IVol_lim,',"'
   OPEN(UNIT=DataUnitp,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'ATRID.plt')
   write(DataUnitp,"(6F9.3,2F9.3,A,F7.1,i3,' 0')") BBx_min(1),BBx_max(1),BBx_min(2),BBx_max(2),BBx_min(3),BBx_max(3), tMin, tMax,&
         ' NoBox ', TimeBase, 0 ! gleRead:  xMin xMax yMin yMax zMin zMax tMin tMax ZmBx$ t_offst
   Write(DataUnitp,"(A,I5,1x,A,F6.1,A,F7.1,F7.0,A, 2F7.1, A)") &
      ' 0 ',keep, Trim(Qual)//' 0 '//TRIM(FlashName)//':'//TRIM(OutFileLabel),10.,' -1. ', &
      IVolMean, IVolVar, ' -1. ', Chi2Mean, Chi2Var, ' -1. 0 0 -1. 0 0 0 1.  1.0 !' ! gleRead:  NTracks EventNr Q t_start label$ AmplitudePlot a1 b1 c1 a2 b2 c2 d2
   write(DataUnitp,"(2A)") '!  #   t_src[ms]                   (NEh)[km]                Im_Int Width  Chi^2  Sprd', &
      ' PUn- Plin- Pcircpol  I_lin   Pol_direction(NEh)'
   !
   !OPEN(UNIT=29,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'ATRIDSrc.dat')  ! for plotting
   !write(29,"(6F8.2,2F9.3,A,F7.1,i3,' 0')") BBx_min(2),BBx_max(2),BBx_min(1),BBx_max(1),BBx_min(3),BBx_max(3), tMin, tMax,  &
   !      ' NoBox ', TimeBase, 0 ! gleRead:  xMin xMax yMin yMax zMin zMax tMin tMax ZmBx$ t_offst
   !Write(29,"(A,I5, A,F6.1,A, F7.1, F7.0, A, 2F7.1, A)") ' 0 ',PeakNrTotal, Trim(Qual)//' 0 '//TRIM(OutFileLabel),10.,' -1. ', &
   !   Chi2_Lim, Spread_Lim, ' -1. ', Chi2Mean, Chi2Var, ' -1. 0 0 -1. 0 0 0 1.  1.0 !' ! gleRead:  NTracks EventNr Q t_start label$ AmplitudePlot a1 b1 c1 a2 b2 c2 d2
   !Write(29,"(A,A)") '! sourceNr,  time,     x=east,    y=north,    z=height, Intensity,  ', &
   !   '    spread,       1/rho,        chi^2,      I3/I '
   !
   write(2,"('!',I5,A,A)") PeakNrTotal,', t_src [ms] ,    north=y ,   east=x , height[km], Width,   I_12 ISprd  IVol ',&
      ' Chi^2,   I3%  P_un,P_lin,P_circ; Diffs (samples, [m]) ; St_I/20,       Zen, Azi  , t_ref'
   keep=0
   PolPlotScale=1.
   Do i_Peak=1,PeakNrTotal
      txt='S'
      Message=''
      If(DeSelect(i_Peak)) Then
         Write(2,"(2(A,F7.2), 2(A,F8.1), 2(A,F6.1),2(A,F8.1),2(A,F6.3), A)")   &
            '! Chi^2=', Chi2(i_Peak), ' exceeds max allowed of ',Chi2_lim, &  !  2(A,F7.2),
            !', and/or normalized Intensity spread=',SourceIntSpread(i_Peak), ' exceeds ', Spread_Lim, & !  2(A,F6.3),
            ', and/or I_12=',SourceIntensity(i_Peak), ' below noise level of ', NoiseLevel, & !  2(A,F8.1),
            !', and/or % radial polarization=',I3(i_Peak)*100., ' exceeds ', I3_lim*100., & !  2(A,F6.3),
            ', and/or % Intensty Vol',SourceIntVol(i_Peak), ' exceeds ', IVol_lim, &  !  2(A,F6.1),
            ' ----- '
         txt='!'
      Else
         Sam_peak=PeakStartChunk_ms(i_peak)/sample_ms + PeakPos(i_Peak)  !=ref station time, different from SourceTime_ms
         Do k=1,PeakNrTotal
            If(k.eq.i_Peak) cycle
            If(( (abs(PeakStartChunk_ms(k)/sample_ms + PeakPos(k) -Sam_peak) .lt. (PeakWidth(i_Peak)+PeakWidth(k))/4)) .and. &
            !If(( (abs(SourceTime_ms(k) -SourceTime_ms(i_Peak)) .lt. (PeakWidth(i_Peak)+PeakWidth(k))/4)) .and. &
                  (.not. DeSelect(k)) ) Then ! check for two similar sources
               Write(txt2,"(I4)") k
               !write(2,*) '!Message:"',TRIM(message),'"', k
               !flush(unit=2)
               If(i_Peak.gt.k) Then
                  If(Message .eq. '') Then
                     Message='Similar position as source #'//txt2
                  Else
                     Message=TRIM(message)//','//txt2
                  EndIf
               EndIf
               If( ABS(PeakWidth(i_Peak)-PeakWidth(k)).lt. 20) Then ! Similar position and width; check chi^2
                  StI=sqrt( sum((SourcePos(1:3,i_Peak)-SourcePos(1:3,k))**2 ) ) ! just a dummy
                  If(Chi2(i_Peak) .gt. Chi2(k) .and. StI.lt.MinSrcDist ) Then  ! Disqualify this source
                     Write(txt3,"(F8.2)") StI
                     Message='Disqualify; Similar position & width as source #'//txt2//', rel.dist='//txt3//'m, but worse chi^2'
                     txt='!'
                     DeSelect(i_Peak)=.true.
                     exit
                  EndIf
               EndIf
            EndIf
         Enddo
         If(Message .ne. '') write(2,"(A,A)") '! Note: ', trim(Message)
      EndIf
      !write(2,*) 'OutputPkInterfer: txt,Message',i_Peak, txt,Message
      !Flush(unit=2)
      StI = Real(Stk_NEh(1,1,i_Peak) + Stk_NEh(2,2,i_Peak) + Stk_NEh(3,3,i_Peak))  ! =\Lambda_0  from  PHYSICAL REVIEW E 66, 016615 ~2002!
      ! Note that the stokes parameters are averaged over the window, summed and divided by the window length.
      ! Stokes thus gives the intensity per sample, averaged over the complete window.
      I20=StI*PeakWidth(i_Peak)/20.
      ! Since StI is per sample, I20 is proportional to the total (summed) intensity in the window.
      If(txt.ne. '!') keep=keep+1
      !
      write(2,"(A,I4, ', ',F11.6,', ',3F11.5,', ',I4, ', ',F8.1,2F6.1,F6.2,', ',4F6.1, I3, 3F8.3 ,1pg12.4, 0p2f7.2, F11.6)") txt, &
         i_peak, SourceTime_ms(i_Peak), SourcePos(1:3,i_Peak)/1000.,  PeakWidth(i_Peak), &
         SourceIntensity(i_Peak), SourceIntSpread(i_Peak), SourceIntVol(i_Peak), Chi2(i_Peak), &
         I3(i_Peak)*100., SourceUn(i_Peak)*100., SourceLin(i_Peak)*100., SourceCirc(i_Peak)*100., &
         Peakpos(i_Peak)-OrignlPeakpos(i_Peak), SourcePos(1:3,i_Peak)-OrignlSourcePos(1:3,i_Peak), &
         I20, SourcePolZen(1,i_Peak), SourcePolAzi(1,i_Peak), &
         PeakStartChunk_ms(i_peak)+PeakPos(i_Peak)*sample_ms -TimeBase
      !
      !  For storing data for later processing in whatever way:
      write(DataUnit,"(A,I4, ',',F11.6,3(',',F11.5),', ',I4, ',', F7.2,',', F7.1, ',',1pg12.4, ',',0pF6.1, 3(',',F6.1), &
          2(',',f7.2),1p, 6(' , (',g12.4','g12.4,')' ))") txt, &
         i_peak, SourceTime_ms(i_Peak), SourcePos(1:3,i_Peak)/1000.,  PeakWidth(i_Peak), Chi2(i_Peak), SourceIntVol(i_Peak), &
         I20, I3(i_Peak)*100., SourceUn(i_Peak)*100., SourceLin(i_Peak)*100., SourceCirc(i_Peak)*100., &
         SourcePolZen(1,i_Peak), SourcePolAzi(1,i_Peak),  Stk_NEh(1,1:3,i_Peak), Stk_NEh(2,2:3,i_Peak), Stk_NEh(3,3,i_Peak)
      !
      ! For making a plot of the just analyzed sources:
      !write(29,"(i6,' ',4(f11.5,' '),g13.6, 2x ,4(1pg13.4), I5)") i_Peak, SourceTime_ms(i_Peak), SourcePos(2,i_Peak)/1000., &
      !   SourcePos(1,i_Peak)/1000., SourcePos(3,i_Peak)/1000., ABS(SourceIntensity(i_Peak)), &
      !   SourceIntSpread(i_Peak), SourceIntVol(i_Peak), Chi2(i_Peak), I3(i_Peak)*100.,  PeakWidth(i_Peak)
      !
      If(.not. DeSelect(i_Peak)) Then
         Lin1= I20* SourceLin(i_Peak)  ! /100.
         polN=PolPlotScale*sin(SourcePolZen(1,i_Peak)*pi/180.) * cos(SourcePolAzi(1,i_Peak)*pi/180.)
         PolE=PolPlotScale*sin(SourcePolZen(1,i_Peak)*pi/180.) * sin(SourcePolAzi(1,i_Peak)*pi/180.)
         Polh=PolPlotScale*cos(SourcePolZen(1,i_Peak)*pi/180.)
         write(DataUnitp,"(I4,F14.6, 3F12.5, F12.1, I5, 2F7.2,  3F6.3, F10.1, 3f7.3)") &
            i_peak, SourceTime_ms(i_Peak), SourcePos(1:3,i_Peak)/1000.,  ABS(SourceIntensity(i_Peak)), &
            PeakWidth(i_Peak), Chi2(i_Peak), SourceIntVol(i_Peak), &
            SourceUn(i_Peak), SourceLin(i_Peak), SourceCirc(i_Peak), Lin1, PolN, PolE, Polh
      EndIf
      !write(2,*) 'OutputPkInterfer: DeSelect', i_Peak, DeSelect(i_Peak)
      Flush(unit=2)
   EndDo
   Write(2,*) 'Figure-Data File written: "', trim(DataFolder)//TRIM(OutFileLabel)//'ATRID.plt','"'
   Write(2,*) 'Complete Sources-Data File written: "', trim(DataFolder)//TRIM(OutFileLabel)//'ATRID.csv','"'
   Write(2,*) keep,' sources kept out of ',PeakNrTotal
   !
   !Close(Unit=29)
   Close(Unit=DataUnitp)
   Close(Unit=DataUnit)
   !
   If(PeakNrTotal.gt.1) then
      Call GLEplotControl(PlotType='SrcsPltPol', PlotName='ATRID'//TRIM(OutFileLabel)//'IPol', &
         PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel)//'ATRID' )
      Call GLEplotControl(PlotType='SrcsPltPol', PlotName='ATRID'//TRIM(OutFileLabel)//'IPolEA', &
         PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel)//'ATRID  1' )
      Call GLEplotControl(PlotType='SrcsPltChi2', PlotName='ATRID'//TRIM(OutFileLabel)//'Chi2', &
         PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel) )
      Call GLEplotControl(PlotType='SrcsPltISpr', PlotName='ATRID'//TRIM(OutFileLabel)//'ISpr', &
         PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel) )
      Call GLEplotControl(PlotType='SrcsPltLoc', PlotName='ATRID'//TRIM(OutFileLabel), &
         PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel)//'ATRID' )
      ! write(10,"('gle -d pdf -o ',A,'-InfImaMx_',I1,'.pdf ${UtilDir}/SourcesPlot.gle ${FlashFolder}/files/',A)") &!
      !   TRIM(OutFileLabel), i_eo, TRIM(OutFileLabel)//'IntfSpecPowMx'//TRIM(txt)
      !Call GLEplotControl(PlotType='EIPolariz', PlotName='IntfPol'//TRIM(txt)//TRIM(OutFileLabel), &
      !   PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel) )
   EndIf

End Subroutine OutputPkInterfer
!-----------------------------------------
Subroutine PeakCurtain(i_peak)
   use Constants, only : dp !,sample,c_mps, pi, sample_ms
   use DataConstants, only : Ant_nrMax
   use Interferom_Pars, only :  Nr_IntFerMx, Nr_IntferCh ! the latter gives # per chunk
   use Interferom_Pars, only : Cnu_p0, Cnu_t0, Cnu_p1, Cnu_t1, IntfNuDim, AntPeak_OffSt
   Use Interferom_Pars, only : N_smth, N_fit, i_chunk
   Use Interferom_Pars, only : dStI, dStI12, dStQ, dStU, dStV, dStI3, dStU1, dStV1, dStU2, dStV2, Chi2pDF
   Use Interferom_Pars, only : PolZen, PolAzi, PolMag, PoldOm, Stk_NEh  !  calculated in "PolTestCath", called from EI_PolGridDel when (Outpt.ge.1); (PolMag(k), PolZen(k), PolAzi(k), PoldOm(k), k=1,3)
   use ThisSource, only : SourcePos, PeakPos, ExclStatNr, PeakNrTotal
   Use Interferom_Pars, only : Alloc_EInterfImag_Pars, DeAlloc_EInterfImag_Pars
   Implicit none
   Integer, intent(in) :: i_Peak
   Real(dp) :: FitDelay(1:Ant_nrMax)
   Real(dp) :: ChiSq, DelChi(-N_smth:+N_smth,1:2*Nr_IntFerMx)
   Integer :: IntfBase, Outpt
   Character(len=8) :: Label  ! should not be longer
   !
   !Outpt=2  ! will generate also a curtain plot
   Outpt=1  ! will not a curtain plot, but calculate carthesial pol. observables
   If(PeakNrTotal.le.15) Outpt=2
   !Outpt=0  ! minimal output
   N_fit=N_smth
   write(Label,"('Peak ',i3.3)") i_Peak
   IntfBase= Peakpos(i_Peak) - IntfNuDim
   Call Alloc_EInterfImag_Pars
   Call EISelectAntennas(i_chunk)  ! select antennas for which there is an even and an odd one.
   Call GetInterfFitDelay(i_chunk, FitDelay)
   Call EI_PolSetUp(Nr_IntFerCh(i_chunk), IntfBase, i_chunk, SourcePos(1,i_Peak), &
      AntPeak_OffSt(1,1), Cnu_p0(0,1,1), Cnu_t0(0,1,1), Cnu_p1(0,1,1), Cnu_t1(0,1,1))  ! created variables
   Call EI_PolGridDel(Nr_IntFerCh(i_chunk), FitDelay, IntfNuDim, i_chunk, SourcePos(1,i_Peak), &
      AntPeak_OffSt(1,1), Cnu_p0(0,1,1), Cnu_t0(0,1,1), Cnu_p1(0,1,1), Cnu_t1(0,1,1), &
      Outpt, DelChi, Label, ExclStatNr )
   Call DeAlloc_EInterfImag_Pars
   !
   Return
   !
End Subroutine PeakCurtain
!-----------------------------------------
!==========================================
!-----------------------------------------------
Subroutine BestOnes(I_1,I_2,I_3, A, I_best, A_best, N_best)
   use constants, only : dp
   Implicit none
   Integer, intent(in) :: I_1,I_2,I_3, N_best
   Real(dp), intent(in) :: A
   Integer, intent(inout) :: I_best(1:3,1:N_best)
   Real(dp), intent(inout) :: A_best(1:N_best)
   Integer :: k,m
   !
   Do k=1,N_best
      If(A .gt. A_best(k) ) Then
         Do m=N_best-1,k,-1
            A_best(m+1)=A_best(m)
            I_best(1:3,m+1)=I_best(1:3,m)
         EndDo
         A_best(k)=A
         I_best(1,k)=I_1
         I_best(2,k)=I_2
         I_best(3,k)=I_3
         Exit
      EndIf
   EndDo
   Return
End Subroutine BestOnes
!=====================================================
Subroutine WindowShift(N_Smth, smooth, Intens, Sampl_shift)
   use constants, only : dp
   use DataConstants, only : DataFolder, OutFileLabel
   Implicit none
   Integer, intent(in) :: N_Smth
   Real(dp), intent(in) :: smooth(-N_Smth:N_Smth), Intens(-N_Smth:N_Smth)
   Integer, intent(out) :: Sampl_shift
   !Character(len=*), intent(in) :: txt
   Integer :: i, i_d, i_u, move
   Real(dp) :: A_s, A_d, A_u, Min
   !Integer, parameter :: PeakCrit =1 ! Maximize intensity in window
   !Integer, parameter :: PeakCrit =2 ! Minimize intensity at window borders, incremental
   Integer, parameter :: PeakCrit =3 ! Minimize intensity at window borders, over the full smoothing window
   !Find window edges
   !Intens(-N_Smth:N_Smth)=TIntens_pol(1,-N_Smth:N_Smth) + TIntens_pol(2,-N_Smth:N_Smth)
   !
   Do i=-N_Smth,0
      If(smooth(i) .gt. 0.45/N_smth) exit
   EndDo
   i_d=i ! down side of window
   Do i=N_Smth,0,-1
      If(smooth(i) .gt. 0.45/N_smth) exit
   EndDo
   i_u=i ! upper side of window
   !
   i=0
   If(PeakCrit.eq.1) then
      A_d=Intens(i_d-1) - Intens(i_u) ! Gain in intesity in window when it moves down
      A_u=Intens(i_u+1) - Intens(i_d) ! Gain in intesity in window when it moves up
      Min=(Intens(i_u)+Intens(i_d))/40. ! Minimum intensity-change to make moving worthwhile
      If(A_u.gt. A_d) Then ! move window up,
         Do i=0,N_smth ! Move=+1
            If((i_u+i) .eq. N_smth) exit
            If((Intens(i_u+i+1) - Intens(i_d +i)) .lt. Min) exit
         EndDo
      Else  ! Move down
         Do i=0,-N_smth,-1 ! Move down; i is negative
            If((i_d+i) .eq. -N_smth) exit
            If((Intens(i_d +i-1)-Intens(i_u+i))  .lt. Min) exit
         EndDo
      EndIf
      Sampl_shift=i
   ElseIf(PeakCrit.eq.2) Then
      A_s=Intens(i_d) + Intens(i_u) ! Sum border intesities of window
      A_d=Intens(i_d-1) + Intens(i_u-1) ! Sum border intesities of window when it moves down
      A_u=Intens(i_d+1) + Intens(i_u+1) ! Sum border intesities of window when it moves up
      If(A_u.lt. A_d) Then ! move window up,
         Do i=1,N_smth ! Move=+1
            If((i_u+i) .eq. N_smth) exit
            If(A_u .gt. A_s) exit
            A_s=A_u
            A_u=Intens(i_u+i+1) + Intens(i_d +i+1)
         EndDo
         Sampl_shift=i
      Else  ! Move down
         Do i=-1,-N_smth,-1 ! Move down; i is negative
            If((i_d+i) .eq. -N_smth) exit
            If(A_d .gt. A_s) exit
            A_s=A_d
            A_d=Intens(i_d+i-1) + Intens(i_u +i-1)
         EndDo
         Sampl_shift=i
      EndIf
   ElseIf(PeakCrit.eq.3) Then
      A_s=(Intens(i_d) + Intens(i_u))
      Sampl_shift=0
      Do i=(-N_smth-i_d), (N_smth-i_u)
         If(A_s .gt. (Intens(i_d+i) + Intens(i_u+i)) ) Then
            Sampl_shift=i
            A_s = (Intens(i_d+i) + Intens(i_u+i) )
            !write(2,*) '!t-shift option3:', i, A_s
         EndIf
      EndDo
   EndIf
   !write(2,*) '!Shift window', Sampl_shift, i_d, i_u, A_s,A_u, A_d
   !Flush(unit=2)
   !Write(2,*) 'Window time shift=', Sampl_shift, i_d, i_u, Sum(Intens(i_d+Sampl_shift:i_u+Sampl_shift)), &
   !      Intens(i_d+Sampl_shift), Intens(i_d+Sampl_shift+1), Intens(i_u+Sampl_shift-1), Intens(i_u+Sampl_shift)
   !write(2,"(A, I3, 80F10.3)") 'Smooth',N_Smth, smooth(-N_Smth:N_Smth)
   !write(2,"(A,  80F10.3)") 'TIntens_pol1', Intens(-N_Smth:N_Smth)
   !write(2,"(A,  80F10.3)") 'TIntens_pol2', TIntens_pol(2,-N_Smth:N_Smth)
   !write(2,"(A,  80F10.3)") 'TIntens_pol3', TIntens_pol(3,-N_Smth:N_Smth)
   !
   !OPEN(UNIT=30,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'PkIntf-'//txt//'.csv')
   !write(30,"(4(I3,','), 3(A,',') )")  N_Smth, Sampl_shift, i_d, i_u, ' ! '
   !write(30,"(A)")  ' ! i,    windowing,  Intensity,  Intensity for polarization directio = 1, 2, 3 '
   !Do i=-N_Smth,N_Smth
   !   write(30,"(I3,',',5(g12.4,','))") i,smooth(i), Intens(i)
   !EndDo
   !Close(UNIT=30)
   !
   Return
End Subroutine WindowShift
!---------------------------------------
!=====================================================
Subroutine WindowOptimal(N_Smth, smooth, Intens, N_Smth_new, Sampl_shift, PkW, A_Bckgr)
   ! Optimize window by
   ! a)  Find sample with max intensity
   ! b)  For up and down sides (within N_Smth):
   !  1) Find min intensity
   !  2) Check for other min between max_sample and min sample, at least 3 away from each
   !  3) When new min less than 1.0=bckgr, adopt this and repeat step 2 & 3
   !     This procedure is the same as taking the first min of less than 1.0=bckgr, or the best one case this is not there.
   ! c)  Take N_smooth as difference between up and down minima
   ! d)  Put center at half
   use constants, only : dp
   use DataConstants, only : DataFolder, OutFileLabel
   Implicit none
   Integer, intent(in) :: N_Smth, PkW
   Real(dp), intent(in) :: smooth(-N_Smth:N_Smth), Intens(-N_Smth:N_Smth), A_Bckgr
   Integer, intent(out) :: N_Smth_new, Sampl_shift
   !Character(len=*), intent(in) :: txt
   Integer :: i, j, i_d, i_u, i_max, MinLength, i_min
   Real(dp) :: A_max, A_d, A_u, A, A_0
   Integer, parameter :: N_Smth_mh=4, N_Smth_min =2*N_Smth_mh ! Minimize size of window, should be even
   !Real(dp), parameter :: A_Bckgr=1.d0
   !
   N_Smth_new=N_Smth
   Sampl_shift=0
   !write(2,*) '!WindowOptimal on entry:', N_Smth, PkW, N_Smth_min
   If(PkW.lt. N_Smth_min) Then
      write(2,*) '*** Do not optimize window width, keep',N_Smth
      Call WindowShift(N_Smth, smooth, Intens, Sampl_shift)
      Return
   EndIf
   !
   A_max=0  !  Find position max within bounds
   Do i=-N_Smth/2, N_Smth/2
      A=Intens(i)
      If(A.gt. A_max) Then
         i_max=i
         A_max=A
      EndIf
   EndDo
   !write(2,*) '!WindowOptimal; max', i_max,A_max, A_Bckgr
   !
   A_d=A_max ! find down side of window
   i_min=-N_Smth_mh
   If(i_max.lt.0) i_min=i_max-1
!   If(i_min.gt.0) i_min=0
   i=i_min+1
   i_d=i ! down side of window
   A_0=Intens(i)
   A_d=A_0
   Do i=i_min-1, -N_Smth, -1
      A=Intens(i)
      !write(2,*) '!WindowOptimal; d', i,A, A_d, ' i_d=', i_d
      If(A.lt.A_Bckgr) Then
         Do j=i, -N_Smth, -1  ! 06/05/2025 changed: find first minimum in this stretch, othewise take j=0
            A=Intens(j)
            If(A.gt.A_0) cycle  ! running uphill
            If(A.gt.A_d) exit  ! running uphill, check time to stop
            A_d=A
            i_d=j ! down side of window
            !write(2,*) '!WindowOptimal; dloop', i,A, A_d, i_d
         EndDo
         !write(2,*) '!WindowOptimal; d1', j,A, A_d, ' i_d=', i_d
         exit
      ElseIf(A.le. A_d) Then
         i_d=i ! down side of window
         A_d=A
         !write(2,*) '!WindowOptimal; d2', A_d, ' i_d=', i_d
      ElseIf(A.gt. A_max) Then
         exit
      EndIf
   EndDo
   !
   A_u=A_max ! find up side of window
   i_u=-N_Smth
   i_min=i_d+N_Smth_min
   If(i_min.lt.i_max) i_min=i_max+1 ! make sure to start beyond max, and minimal distance from low end
   If(i_min.lt.1) i_min=1
   i=i_min-1
   i_u=i ! down side of window
   A_0=Intens(i)
   A_u=A_0
   Do i=i_min+1, N_Smth
      A=Intens(i)
      !write(2,*) '!WindowOptimal; u', i,A, A_u, ' i_u=', i_u
      If(A.lt.A_Bckgr) Then
         Do j=i, N_Smth  ! 06/05/2025 changed: find first minimum in this stretch, othewise take j=0
            A=Intens(j)
            If(A.gt.A_0) cycle  ! running uphill
            If(A.gt.A_u) exit  ! running uphill again, time to stop
            A_u=A
            i_u=j ! up side of window
         EndDo
         !write(2,*) '!WindowOptimal; u1', j,A, A_u, ' i_u=', i_u
         exit
      ElseIf(A.le. A_u) Then
         i_u=i ! up side of window
         A_u=A
         !write(2,*) '!WindowOptimal; u2', A_u, ' i_u=', i_u
      ElseIf(A.gt. A_max) Then
         If(i_u.le.-N_Smth) Then
            i_u=N_Smth_mh
            Do while(Intens(i_u).gt.Intens(i_u+1))
               i_u=i_u+1
            Enddo
            i_u=i_u-1
            If(i_d .gt. i_u-N_Smth_min) i_d=i_u-N_Smth_min
         EndIf
         !write(2,*) '!WindowOptimal; u',  A_max, ' i_u, i_d=', i_u, i_d
         exit
      EndIf
   EndDo
   !
   Sampl_shift=(i_u+i_d)/2
   N_Smth_new=i_u-i_d+1
   !write(2,*) '!Optimal window', Sampl_shift, N_Smth_new, i_d, i_max, i_u, A_d , A_max,A_u
   !write(2,*) '!Pulse-low:', i_d, i_d+1, Intens(i_d:i_d+1)
   !write(2,*) '!Pulse-hi:', i_u-1, i_u, Intens(i_u-1:i_u)
   !Flush(unit=2)
! on entry WindowOptimal          13          13
 !Optimal window          -1          11          -6          -3           4   4.5137780861166404E-002   5.9222742791737136E-002   4.8830912471730346E-002
 !Pulse-low   7.0790495640559892E-002   4.8984191882427690E-002   4.5137780861166404E-002   5.2131777129244980E-002   5.7873606233532236E-002
 !Pulse-hi   1.9314045380176840E-002   3.2668149596622488E-002   4.8830912471730346E-002   5.7516430333449442E-002   5.2371389649467304E-002
! 5,   106  119; new:   106  116   11

! on entry WindowOptimal          13          13
! !Optimal window           4          11          -1           4           9   3.3274522110728011E-002   5.6717755444165900E-002   6.4469943025060664E-002
! !Pulse-low   3.5349032670011216E-002   3.7887499617648576E-002   3.3274522110728011E-002   2.7534792273900720E-002   2.8214192711843764E-002
! !Pulse-hi   2.8934081923755951E-002   4.2179271183759151E-002   6.4469943025060664E-002   7.5166462330850642E-002   6.1842622164272873E-002
!11,   184  197; new:   189  199   11
   !
   Return
End Subroutine WindowOptimal
!---------------------------------------
Subroutine WindowCutting_N_Smth(N_Smth, N_Window, Intens, DeadEnd, Nmin, SourcePos_C, SourceWidth_C, N_Pks, NoiseLevel)
!  This cutting routine selects the windows close in length to N_Smth
   use constants, only : dp
   !use DataConstants, only : DataFolder, OutFileLabel
   Implicit none
   Integer, intent(in) :: N_Smth, N_Window, DeadEnd, Nmin
   Real(dp), intent(in) :: Intens(1:N_Window), NoiseLevel
   Integer, intent(out) :: SourcePos_C(0:Nmin), SourceWidth_C(0:Nmin)
   Integer, intent(out) :: N_Pks
   Integer :: i,k, N_Smth_new, NSmth, Sampl_shift
   Integer, parameter :: MinLength=7
   Real(dp) :: smooth(-N_Smth:N_Smth)  ! a dummy
   NSmth=2*(N_Smth/2)  ! make sure N_Smth is even
   k=0
   If(NSmth.lt.MinLength) stop "N_Smth too small"
   If(NSmth.gt.2*DeadEnd) stop "Dead end"
   i=DeadEnd+NSmth/2
   !write(2,*) i, DeadEnd, NSmth
   Do while ( i.lt.(N_Window-2*DeadEnd))
      !write(2,*) '!WindowOptimal; Minimab',i,k
      !      Flush(unit=2)
      If(k .ge. Nmin) exit
      k=k+1
      ! Do not skip any part of the trace searching for the max
      !write(2,*) 'Cutting_N_',k,' search window:',i-NSmth/2, i+NSmth+NSmth/2+1
      Call WindowOptimal(NSmth, smooth, Intens(i-NSmth/2:i+NSmth+NSmth/2+1), N_Smth_new, Sampl_shift, NSmth, NoiseLevel)
      SourcePos_C(k)=i + NSmth/2 +Sampl_shift
      SourceWidth_C(k)=N_Smth_new
      !write(2,*) 'Cutting_N_',k,' found cutting:',SourcePos_C(k)-N_Smth_new/2, SourcePos_C(k)+N_Smth_new/2
      i=SourcePos_C(k)+ N_Smth_new/2
   EndDo
   N_Pks=k
End Subroutine WindowCutting_N_Smth
!---------------------------------------
Subroutine WindowCutting(N_Window, Intens, DeadEnd, Nmin, SourcePos_C, SourceWidth_C, SourceLabel_C, N_minima, NoiseLevel)
!  This cutting routine selects the windows with minimal length to cover a pulse
   use constants, only : dp
   !use DataConstants, only : DataFolder, OutFileLabel
   Implicit none
   Integer, intent(in) :: N_Window, DeadEnd, Nmin
   Real(dp), intent(in) :: Intens(1:N_Window), NoiseLevel
   Integer, intent(out) :: SourcePos_C(0:Nmin), SourceWidth_C(0:Nmin)
   Logical, intent(out) :: SourceLabel_C(0:Nmin)
   Integer, intent(out) :: N_minima
   Integer :: Minima(0:Nmin)
   !
   Integer :: N_Smth_new, Sampl_shift
   Real(dp), allocatable :: smooth_T(:), Intens_T(:)
   Integer :: i, k, i_max, MinLength, N_Smth_T, i_s1, i_s2
   Real(dp) :: A_max, A_d, A_u, A
   !Integer, parameter :: N_Smth_mh=5, N_Smth_min =2*N_Smth_mh ! Minimize size of window, should be even
   !Real(dp), parameter :: A_Bckgr=1.d0
   !
   !  Find all possible cutting (=minima) places in new time trace
   MinLength=10
   N_minima=0
   Minima(N_minima)=DeadEnd
   i=Minima(N_minima)+MinLength
!   write(2,*) '!WindowOptimal; Minimaa',i,Sampl_shift, N_Smth_new,N_Smth
!   Flush(unit=2)
   Do while ( i.lt.(N_Window-2*DeadEnd))
!      write(2,*) '!WindowOptimal; Minimab',i,Sampl_shift, N_Smth_new,N_minima
!      Flush(unit=2)
      A_d=Intens(i-1)
      A=Intens(i)
      A_u=Intens(i+1)
      If((A.lt. A_u) .and. (A.lt. A_d)) Then  ! found a minimum
         N_minima=N_minima+1
         Minima(N_minima)=i
         i=i+MinLength
         If(N_minima+1 .ge. Nmin-1) exit
      Else
         i=i+1
      EndIf
   EndDo
   !Minima(1:N_minima)= Minima(1:N_minima)
   !  Add the end point:
   If(Minima(N_minima) .lt. (N_Window-MinLength-DeadEnd)) N_minima=N_minima+1
   Minima(N_minima)=N_Window-DeadEnd
   !write(2,"(A,I3,A,20I4 )") '!Found',N_minima, ' minima in time trace@', Minima(1:N_minima)
   !write(2,"(A,25x,20I4 )") '!Slicing lengths:', Minima(2:N_minima)-Minima(1:N_minima-1)
   !
   ! Testing:
   If(N_minima.gt.2) Then
      !write(2,*) 'suggestions for fine cutting of this source:'
      Do k=1,N_minima
         i_max=(Minima(k) + Minima(k-1))/2
         N_Smth_T= Minima(k) - Minima(k-1)
         Allocate( smooth_T(-N_Smth_T:N_Smth_T), Intens_T(-N_Smth_T:N_Smth_T) )
         smooth_T(:)=0.
         smooth_T(-N_Smth_T/2:N_Smth_T/2)=1./N_Smth_T
         Intens_T(-N_Smth_T:N_Smth_T)=Intens(i_max-N_Smth_T:i_max+N_Smth_T)
         Call WindowOptimal(N_Smth_T, smooth_T, Intens_T, N_Smth_new, Sampl_shift, N_Smth_T, NoiseLevel)
         DeAllocate( smooth_T, Intens_T )
         SourcePos_C(k)=i_max+Sampl_shift
         SourceWidth_C(k)=N_Smth_new
         i_s1=SourcePos_C(k)-(SourceWidth_C(k)-1)/2
         i_s2=SourcePos_C(k)+  (SourceWidth_C(k))/2
         !write(2,"(I5, ', ',2I6,'; new: ',4I6, I4 )") k, Minima(k-1), Minima(k), i_s1, i_s2, SourceWidth_C(k)
         !write(2,*) '!orig-lo', Intens(i_s1-2:i_s1+2)
         !write(2,*) '!orig-hi', Intens(i_s2-2:i_s2+2)
         !
      EndDo
      SourceLabel_C(1)=.true.
      Do i=2,N_minima-1
         SourceLabel_C(i)=.true.
         Do k=1,i-1
            If(( abs(SourcePos_C(k)-SourcePos_C(i)) .lt. (SourceWidth_C(i)+SourceWidth_C(k))/4) .and. &
               ( ABS(SourceWidth_C(i)-SourceWidth_C(k)).lt. MinLength) ) Then ! Similar position and width
               SourceLabel_C(i)=.false.
               write(2,*) '! sub-peak #',i, 'not passable because too similar to #',k
               exit
            EndIf
         Enddo
      EndDo
   EndIf
   !SourcePos_C(1:N_minima)=SourcePos_C(1:N_minima)+N_Smth/2
   !stop
   !
   Return
End Subroutine WindowCutting
!---------------------------------------
! ======================================================================
! ======================================================================
Subroutine IntensitySpread(Grd_best, Int_best, Grd, N_best, N_Int, Npix, BoxFineness, &
         Spread_Grd, Spread_Dist, Spread_SclLong, Spread_Aspect, Spread_IntVol, verbose)
!  Considerations (March 2025):
!     N_best; Should be large enough to (hardly) hit the intensity limit. The intended correction factor for hitting
!        this limit, "Int_scale", is not doing a good job. The number of voxels above a certain threshold appears to
!        grow exponentially with intensity, i.e. much faster than the supposed ^3.
!     Frac_Int=0.8d0; This should be small enough to test with a reasonable statistics the intensity contours around
!        the max. However, for a too small value the number of higher-intensity sources will over saturate and, more
!        importantly, the calculated spread becomes too sensitive to what happens at the edges of the image cube.
!     Npix; The image cube should be large enough that it can contain a typical point=spread function. This implies
!        at least +/- 10m in the horizontal plane and +/- 20m vertically.
!  Dependencies for an impulsive source on a dart-leader background:
!     Changing Frac_Int results in change of I-voxels Prop to about (1-Frac_Int)^3 or a bit faster
!     Changing (@Frac_Int=0.8) fine grid size from ~20 to ~10 results in 1 to 19 edge cases and 458 to 402 I_voxels
!     Changing (@Frac_Int=0.8) Coarse grid size from ~90 to ~40 results in 0 to 0 edge cases and 23 to 19 I_voxels
!     For 0.75 the numbers are worse.
!     One may conclude that 1 edge case are equivalent to 3-4 I-Voxels
    use constants, only : dp
   Implicit none
   Integer, intent(in) :: Grd_best(1:3,1:N_best), N_best, Npix(1:3,2)
   Real(dp), intent(in) :: BoxFineness, Int_best(1:N_best), Grd(1:3)
   Integer, intent(out) :: N_Int
   Real(dp), intent(out) :: Spread_Grd, Spread_Dist, Spread_SclLong, Spread_Aspect, Spread_IntVol
   logical, intent(in) :: verbose
   Integer :: i,j, N_side
   Real(dp) :: Frac_Int=0.8d0, Max_grd, Max_dist, DGrd, DGrd_i, DGrd_j, Dist, curv_dist, curv_Grd, curv, Int_scale
   Real(dp) :: WGrd_best(1:N_best), Int_side, Imax_Side
   !
   Spread_Grd=0.
   Spread_Dist=0.
   N_Int=1
   Int_scale=1.
   Do i=2, N_best
      If( Int_best(i) .lt. Frac_int*Int_best(1)) exit
      N_Int=i
      Do j=1,i-1
         Spread_Grd=Spread_Grd + SUM( (Grd_best(1:3,i)-Grd_best(1:3,j))**2 )
         Spread_Dist=Spread_Dist + SUM( ((Grd_best(1:3,i)-Grd_best(1:3,j))*Grd(1:3))**2 )
      EndDo
   EndDo
   If(   N_Int .gt. 1) Then
      Spread_Grd=sqrt( Spread_Grd/(0.5*N_Int*(N_Int-1)) )
      Spread_Dist=sqrt( Spread_Dist/(0.5*N_Int*(N_Int-1)) )
      !Int_scale=sqrt((1.-Int_best(N_Int)/Int_best(1)))  ! Intensity falls with square distance
      Int_scale=(1.-Int_best(N_Int)/Int_best(1))  ! On the basis of phenomenology (March 2025)
   EndIf
   !
   !  Explore the best ones that lie on the sides of intensity contour cube
   N_side=0
   Int_side=0.
   Imax_Side=0.
   Do i=1,3
      Do j=1,N_Int
         If((Grd_best(i,j).eq. Npix(i,1)) .or. (Grd_best(i,j).eq. Npix(i,2)) ) Then
            N_Side=N_side+1
            If(Int_best(j).gt.Imax_Side) Imax_Side=Int_best(j)
            Int_Side=Int_Side+Int_best(j)
         EndIf
      EndDo
   EndDo
   Int_Side=Int_Side/Int_best(1)
   Imax_Side=Imax_Side/Int_best(1)
   If(N_Side .gt. 0) Then
      Write(2,"(I5,A,F4.1,A,F4.1, A,F4.1)") N_Side,' side-voxels have an intensity above the limit of ',Frac_int*100,&
         '% with max of ', Imax_Side*100., '% and average of ', 100.*Int_Side/N_side
   Else
      Write(2,"(A,F4.1,A)") 'Fractional intensity at edges of image cube below limit of ',Frac_int*100.,'%'
   EndIf
   !
   !
   Max_grd=0.
   Max_dist=0.
   curv_grd=0.
   curv_dist=0.
   curv=0.
   Do i=2, N_Int
      Do j=1,i-1
         DGrd= SUM( (Grd_best(1:3,i)-Grd_best(1:3,j))**2 )
         If(DGrd.gt.Max_grd) Max_grd=DGrd
         Dist=SUM( ((Grd_best(1:3,i)-Grd_best(1:3,j))*Grd(1:3))**2 )
         If(Dist.gt.Max_Dist) Max_Dist=Dist
         If(j.eq.1) Then
            DGrd_i= DGrd
            Curv_Grd=Curv_Grd + (1.-Int_best(i)/Int_best(1))/DGrd_i
            Curv_Dist=Curv_Dist + (1.-Int_best(i)/Int_best(1))/Dist
         Else
            DGrd_j= SUM( (Grd_best(1:3,j)-Grd_best(1:3,1))**2 )
            If(DGrd_i .ne. Dgrd_j) Then
               Curv=Curv + (Int_best(j)-Int_best(i))/(DGrd_i-Dgrd_j)
               !write(2,*) 'non-equal distance:', DGrd_i-Dgrd_j,Int_best(j)-Int_best(i),Curv
            !Else
            !   write(2,*) 'equal distance:', i,j,Int_best(j)-Int_best(i)
            EndIf
         EndIf
      EndDo
   EndDo
   !
   !write(2,"(A, 100F10.3)")    '!Intensity: ',Int_best(:)
   !write(2,"(A, 100(3I3,';'))") '!IntensityGrd:',Grd_best(:,:)
   Spread_IntVol=(N_Int*BoxFineness**3+ 13*N_Side*BoxFineness**2)/(10.*(10.*Int_scale)**4)  ! 13, ^4 are phenomenologic
   Write(2,"(A,2F5.2, A,I4,F11.2,F9.2, A,I4,F11.2,F9.2, A,F6.1)") 'Frac Int=',Int_best(N_Int)/Int_best(1), Imax_Side, &
      ', N_side=',N_Side,N_Side*BoxFineness**2, (13*N_Side*BoxFineness**2)/(10.*(10.*Int_scale)**4), &
      ', N_vol=', N_Int, N_Int*BoxFineness**3, (N_Int*BoxFineness**3)/(10.*(10.*Int_scale)**4), &
      ', Result I_vol=', Spread_IntVol
   If(verbose) Then
      write(2,"(A,F5.1,A,F5.1,A,F5.1,A,I3,A,F5.3)") 'Spread [grd], noWeight=', Spread_Grd,', IntWeight=', &
            1/sqrt(Curv_Grd/(N_Int-1)), ', maximal dist=',sqrt(Max_grd), ', #=',N_Int, ', Int rat=',Int_best(N_Int)/Int_best(1)
      Write(2,*) 'Intensity scale factor=', Int_scale, ', scaled grid distances', Spread_Grd/Int_scale, sqrt(Max_grd)/Int_scale
      write(2,"(A,F5.1,A,F5.1,A,F5.1,A,I3,A,F5.3)") '  Spread [m], noWeight=', Spread_Dist,', IntWeight=', &
            1/sqrt(Curv_Dist/(N_Int-1)), ', maximal dist=',sqrt(Max_Dist), ', #=',N_Int, ', Int rat=',Int_best(N_Int)/Int_best(1)
      !
      WGrd_best(1:N_best)=(Int_best(1:N_best)-Int_best(N_best))*(Int_best(1:N_best)-Int_best(N_best))
      write(2,*) Int_best(1), Int_best(2:5)/Int_best(1)
      write(2,"(20(3x,3I4,11x))") Grd_best(1:3,1:6)
      Call WeightedIPCA(N_Int, Grd_best, WGrd_best, Spread_SclLong, Spread_Aspect, verbose)  ! long axis of density spread
   !
   Else
      write(2,"(A,F5.1,A,F5.1,A,F5.1,A,I4,A,F5.3,A,F10.3)") 'Spread [grd], Int. scaled=', Spread_Grd/Int_scale,&
            ', Int. Weight=', 1/sqrt(Curv_Grd/(N_Int-1)), ', maximal dist (Scaled)=',sqrt(Max_grd)/Int_scale, ', #=',N_Int, &
            ', Int. ratio=',Int_best(N_Int)/Int_best(1), ', Int. diffuseness=',Spread_IntVol
   EndIf
!   write(2,*) 'not-weighted Spread:', Spread_Grd,'[grd],', Spread_Dist,'[m]', N_Int, Int_best(N_Int)/Int_best(1)
!   write(2,*) 'maximal distances=',sqrt(Max_grd),'[grd],', sqrt(Max_Dist),'[m],', N_Int
!   !write(2,*) 'curvatures=', Curv/(0.5*(N_Int-2)*(N_Int-1)), Curv_Grd/(N_Int-1), Curv_Dist/(N_Int-1)
!   write(2,*) 'Int-Weighted Spread=', & ! 1/sqrt(Curv/(0.5*(N_Int-2)*(N_Int-1))),'[grd],',
!      1/sqrt(Curv_Grd/(N_Int-1)),'[grd],', 1/sqrt(Curv_Dist/(N_Int-1)),'[m]'
   Call LocationPCA(N_Int, Grd_best, Spread_SclLong, Spread_Aspect, verbose)  ! long axis of density spread
   Spread_SclLong=BoxFineness*Spread_SclLong/sqrt(Int_scale)
   !
   Return
End Subroutine IntensitySpread
!-----------------------------------------
!-----------------------------------------
Subroutine WeightedIPCA(PeakNrTotal, Locations, Weight, Spread_SclLong, Spread_Aspect, verbose)
   use Constants, only : dp, pi !,sample,c_mps, pi, sample_ms
   !use DataConstants, only : DataFolder, OutFileLabel
   Implicit none
   Integer, intent(in) :: PeakNrTotal
   Integer, intent(in) :: Locations(1:3,1:PeakNrTotal)
   Real(dp), intent(in) :: Weight(1:PeakNrTotal)
   Real(dp), intent(out) :: Spread_SclLong, Spread_Aspect
   Logical, intent(in) :: verbose
   Integer :: i_source, i,j
   Real(dp) :: LocatCovariance(1:3,1:3), LocateMean(1:3), LocCor(1:3), W,W2
   Real(dp) :: B(3,3), eigval(1:3), eigvec(1:3,1:3)
   !
   If(verbose) write(2,*) 'perform Weighted PCA-analysis on location distribution'
   LocateMean(:)=0.
   W=0.
   Do i_source=1,PeakNrTotal
      LocateMean(:)=LocateMean(:) + Locations(:,i_source)*Weight(i_source)
      W=W + Weight(i_source)
   EndDo
   LocateMean(:)=LocateMean(:)/W
   If(verbose) write(2,*) 'Average position of sources=',LocateMean(:),', #=',PeakNrTotal
   !
   LocatCovariance(:,:)=0.
   W2=0.
   Do i_source=1,PeakNrTotal
      LocCor(:)=(Locations(:,i_source)-LocateMean(:))*Weight(i_source)
      W2=W2 + Weight(i_source)*Weight(i_source)
      Do i=1,3
         LocatCovariance(:,i)=LocatCovariance(:,i) + LocCor(:)*LocCor(i)
      EndDo
   EndDo
   LocatCovariance(:,:)=LocatCovariance(:,:)/W2
   !
   Call Inverse33(LocatCovariance(1,1), B, eigval, eigvec) ! vec(1:3,i) with val(i)
   If(verbose) Then
      Do i=1,3 ! loop over all eigenvalues
         write(2,"(A,I1, A,F8.2, A,3F6.3, A,3F6.3)") 'Axis:',i, &
            ', sqrt(Cov)=',sqrt(eigval(i)),', orientation:',eigvec(1:3,i)
      EndDo
      write(2,*) 'ratio short/long=',eigval(3)/eigval(1)
   EndIf
   Spread_SclLong=sqrt(eigval(1))
   Spread_Aspect=eigval(3)/eigval(1)
   !
   Return
End Subroutine WeightedIPCA
!-----------------------------------------
Subroutine LocationPCA(PeakNrTotal, Locations, Spread_SclLong, Spread_Aspect, verbose)
   use Constants, only : dp, pi !,sample,c_mps, pi, sample_ms
   !use DataConstants, only : DataFolder, OutFileLabel
   Implicit none
   Integer, intent(in) :: PeakNrTotal
   !Real(dp), intent(in) :: Locations(1:3,1:PeakNrTotal)
   Integer, intent(in) :: Locations(1:3,1:PeakNrTotal)
   Real(dp), intent(out) :: Spread_SclLong, Spread_Aspect
   Logical, intent(in) :: verbose
   Integer :: i_source, i,j
   Real(dp) :: LocatCovariance(1:3,1:3), LocateMean(1:3), LocCor(1:3)
   Real(dp) :: B(3,3), eigval(1:3), eigvec(1:3,1:3)
   !
   If(verbose) write(2,*) 'perform PCA-analysis on location distribution'
   LocateMean(:)=0.
   Do i_source=1,PeakNrTotal
      LocateMean(:)=LocateMean(:) + Locations(:,i_source)
   EndDo
   LocateMean(:)=LocateMean(:)/PeakNrTotal
   If(verbose) write(2,*) 'Average position of sources=',LocateMean(:),', #=',PeakNrTotal
   !
   LocatCovariance(:,:)=0.
   Do i_source=1,PeakNrTotal
      LocCor(:)=Locations(:,i_source)-LocateMean(:)
      Do i=1,3
         LocatCovariance(:,i)=LocatCovariance(:,i) + LocCor(:)*LocCor(i)/PeakNrTotal
      EndDo
   EndDo
   !
   Call Inverse33(LocatCovariance(1,1), B, eigval, eigvec) ! vec(1:3,i) with val(i)
   If(verbose) Then
      Do i=1,3 ! loop over all eigenvalues
         write(2,"(A,I1, A,F8.2, A,3F6.3, A,3F6.3)") 'Axis:',i, &
            ', sqrt(Cov)=',sqrt(eigval(i)),', orientation:',eigvec(1:3,i)
      EndDo
      write(2,*) 'ratio short/long=',eigval(3)/eigval(1)
   EndIf
   Spread_SclLong=sqrt(eigval(1))
   Spread_Aspect=eigval(3)/eigval(1)
   !
   Return
End Subroutine LocationPCA
!-----------------------------------------
!=================
Subroutine PeakInterferometerRun(PeakStartTref_ms, i_Peak, BoxSize, BoxFineness, dpix, Npix, ContourPlot )
   use constants, only : dp, ci, pi, c_mps, sample_ms
   Use Interferom_Pars, only : Polar, NewCenLoc,  AmpltPlot
   Use Interferom_Pars, only : N_pix, d_loc, CenLocPol, CenLoc, PixLoc, IntFer_ant
   Use Interferom_Pars, only : SumStrt, SumWindw, NrSlices, SliceLen
   Use Interferom_Pars, only : N_smth, NrPixSmPowTr, i_chunk
   use ThisSource, only : ChunkNr, SourcePos, PeakPos
   use DataConstants, only : OutFileLabel, DataFolder, Time_Dim !, Cnu_dim,
   use Chunk_AntInfo, only : StartT_sam, TimeFrame,   NoiseLevel, RefAnt, Simulation
   use FFT, only : RFTransform_su, DAssignFFT !, RFTransform_CF, RFTransform_CF2CT
   use HDF5_LOFAR_Read, only : CloseDataFiles
   use GLEplots, only : GLEplotControl
   Implicit none
   Integer, intent(in) :: i_Peak
   Real(dp), intent(in) :: PeakStartTref_ms(*), BoxSize, BoxFineness
   Real(dp), intent(out) :: dpix(1:3)
   Integer, intent(out) :: Npix(1:3,1:2)
   logical, intent(in) :: ContourPlot
   integer :: i  !,j, k, i_s, i_loc(1), i_ant, j_corr, i_eo, Station, j_IntFer, i_nu,  nxx
   character(len=3) :: txt
   character(len=100) :: OutFileLabel_original
   Real(dp), save :: CenLoc_old(1:3), GridVolume(1:3)
   Integer, save :: Start_old
   Integer :: Date_T(8)
   Real(dp), external :: tShift_ms
   !    .
   i_chunk=1
   StartT_sam(1)=PeakStartTref_ms(i_peak)/sample_ms
   ! Get grid:
   Polar=.false.
   !
   SumWindw=2*N_smth+1  ! should be even to have window centered at peak of interest
   NrPixSmPowTr=(SumWindw-1)/N_smth-1  ! The number of fine slices for TRI-D imaging
   SumStrt=PeakPos(i_Peak) - SumWindw/2-1
   CenLoc(1:3) = SourcePos(1:3,i_Peak)
   AmpltPlot = 10.
   !write(2,"(A,I5,A,I5,A,F6.3)") 'TRI-D window: SumStrt=',SumStrt, ', SumWindw=', SumWindw, ', AmpltPlot=', AmpltPlot
   If(SumStrt.lt.1000) Then
      write(2,*) 'Window starting sample should be 1000'
      Stop 'Window starting sample should be sufficiently large'
   EndIf
   !
   NrSlices=1  !  depreciated, Dec 2021; used in making summed traces in 'EIAnalyzePixelTTrace' and in 'OutputIntfSlices'
   SliceLen=SumWindw/NrSlices ! same as NrSlices
   !
   CALL DATE_AND_TIME (Values=DATE_T)
   !WRITE(2,"(1X,I2,':',I2.2,':',I2.2,'.',I3.3,A)") (DATE_T(i),i=5,8), 'initializing'
   !i_chunk=ChunkNr(i_Peak)
   TimeFrame=i_chunk  ! not sure this is really used
   !
   !write(2,*) '!Start_old:', Start_old, StartT_sam(1)
   If((CenLoc_old(1) .ne. CenLoc(1)) .or. (CenLoc_old(2) .ne. CenLoc(2)) .or. &
      (CenLoc_old(3) .ne. CenLoc(3)) .or. (Start_old .ne. StartT_sam(1))) Then
      !  Read the time traces
      !write(2,*) 'CenLoc_old(1:3)', i_chunk_old, CenLoc_old(1:3)
      Call RFTransform_su(Time_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      Call AntennaRead(i_chunk, CenLoc)
      If(Simulation.eq."") CALL CloseDataFiles()    ! no more reading of antenna data
      call DAssignFFT()                   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CenLoc_old(1:3) = CenLoc(1:3)
      Start_old = StartT_sam(1)
   EndIf
   !
   If(NoiseLevel.lt.0) NoiseLevel=0.2
   !
   !   Allocate( Smooth(-N_smth:N_smth) )
   !   Allocate( TIntens_pol(1:3,-N_smth:N_smth) )  !  Intensity trace for different pol directions as seen at the core, For use in "PeakInterferoOption"
   !   Call SetSmooth(N_Smth, ParabolicSmooth, Smooth)
   ! ----------------------------------------------------
   If(ContourPlot) Then
      If(BoxFineness.gt.2) Then
         write(txt,"('+',I2.2)") i_peak
      Else
         write(txt,"('-',I2.2)") i_peak
      EndIf
      OutFileLabel_original=OutFileLabel
      OutFileLabel=TRIM(OutFileLabel)//txt
   EndIf
   ! main loop over direction antennas
   GridVolume(1:3)=BoxSize
   GridVolume(3)=2*BoxSize
   Call EI_run(GridVolume, BoxFineness)
   !
   dpix(1:3)=d_loc(1:3)
   Npix(1:3,1:2)=N_pix(1:3,1:2)
   !------------------------------------------
   If(ContourPlot) Then
      !write(2,*) 'prepare plot:', TRIM(OutFileLabel), '--------------'
      !Call GLEplotControl(PlotType='EIContour-QUV', PlotName='EIContourQUV'//TRIM(OutFileLabel), &
      !         PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel))
      Call GLEplotControl(PlotType='EIContour', PlotName='EIContour'//TRIM(OutFileLabel), &
               PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel))
      Call GLEplotControl(CleanPlotFile=TRIM(OutFileLabel)//'_EISpec*.dat')
      Call GLEplotControl(CleanPlotFile=TRIM(OutFileLabel)//'Interferometer*.z')
      Call GLEplotControl(CleanPlotFile=TRIM(OutFileLabel)//'IntfSpecPowMx_d.dat')
      Call GLEplotControl(CleanPlotFile=TRIM(OutFileLabel)//'IntfSpecWin_d.csv')
      OutFileLabel=OutFileLabel_original
   EndIf
   !
   Return
    !
End Subroutine PeakInterferometerRun
!-------------------------------------------------------------
!=========================================

!=========================================
