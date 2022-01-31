!=================================
    Include 'EICallRtns.f90'
    Include 'EIFitter.f90'
!================================
!=================================
Subroutine E_Callibr()
!   Fine-callibrate the timing of the LOFAR stations by fitting a few lightning sources.
!   This should be done for even as well as odd numbered antennas.
!
!   Procedural steps:
!   1) The first (odd or even numbered) antenna is taken as the reference antenna (=RefAnt).
!   2) A few peaks on the RefAnt are selected, preferrably in the same chunk.
!   3) The direction for each source is found by fitting to the superterp, keeping station delays fixed to the LOFAR values.
!   4) The source positions are found by fitting to all core stations, keeping LOFAR station delays.
!   5) Fit core station delays, keeping positions fixed?
!   6) Fit all station delays
!
!   For fitting the Abs(corr) or Real(corr) can be mixed.
!
!  v20: when Real(corr)=.true. and Dual=.true. the polarizations mode will be switched on, "Polariz=.true.".
!        It makes sense to run this mode only after rather accurate source positions have been found in earlier runs.
!        In Polariz mode the cross correlations are calculated for the t (i_tp=1=Theta) and p (i_tp=0=phi) modes.
!      1) This requires treating the even and odd antennas at the same level and it will be assumed that the
!        i_eo=0 (even=Y-dipoles) peak positions hold for the two antenna orientations.
!      2) This also requires that the Y- (i_ant=even) and X- (i_ant=even+1) dipoles are both active at the same
!        antenna location (in the same chunck).
   use constants, only : dp,sample
   use DataConstants, only : PeakNr_dim, ChunkNr_dim, RunMode, Time_Dim
   use Chunk_AntInfo, only : Unique_SAI, Start_time, Tot_UniqueAnt, Ant_Stations
   use Chunk_AntInfo, only : Unique_StatID,  Nr_UniqueStat, Nr_UniqueAnt
   use DataConstants, only : Ant_nrMax
   use ThisSource, only : SourcePos
   use ThisSource, only : PeakNr, PeakNrTotal
   use ThisSource, only : ChunkNr
   use FitParams, only : N_FitPar, N_FitPar_max, Nr_TimeOffset, WriteCalib
   use FitParams, only : Fit_TimeOffsetStat, Fit_TimeOffsetAnt, Fit_AntOffset, N_FitStatTim, FitParam, X_Offset
   Use Interferom_Pars, only : IntFer_ant, StartTime_ms, Nr_IntferCh ! the latter gives # per chunk
   use Interferom_Pars, only : Alloc_EInterfCalib_Pars
   use StationMnemonics, only : Statn_ID2Mnem !, Station_Mnem2ID
    use LOFLI_Input, only : ReadSourceTimeLoc
   use FFT, only : RFTransform_su,DAssignFFT, RFTransform_CF2CT
   Use Calibration, only : WriteCalibration ! was MergeFine
   Implicit none
   Real(dp) :: SourceGuess(3,10)
   integer, save ::  FP_s(0:N_FitPar_max)!, FP(0:4)
   real ( kind = 8 ) :: X(N_FitPar_max)
   !
   integer :: i, j, k, i_ant, i_eo, i_chunk, i_stat, StatID
   integer :: i_loc(1), i_Peak
   Logical :: FitNoSources
   character*10 :: txt
   !character*35 :: Frmt
   real(dp) :: x1,x2,x3,t
   !real(dp) :: dt,B, MTC, dtI
   Integer :: i_SAI
   !
   !       Read appropriate data chunks
   Call RFTransform_su(Time_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   write(*,"(A,i3,A)") achar(27)//'[45m # of data-blocks read in=',ChunkNr_dim,achar(27)//'[0m'
   Do i_chunk=1, ChunkNr_dim
      Call ReadSourceTimeLoc(StartTime_ms, SourceGuess(:,i_chunk))
      If(NINT(SourceGuess(2,i_chunk)).eq.1) &
            write(*,"(A,i3,A)") achar(27)//'[45m # Check ChunkNr_dim',i_chunk,achar(27)//'[0m'
      Start_time(i_chunk)=(StartTime_ms/1000.)/sample  ! in sample's
      Call AntennaRead(i_chunk,SourceGuess(:,i_chunk))
      !
      Call EISelectAntennas(i_chunk)  ! select antennas for which there is an even and an odd one.
      !
   EndDo !  i_chunk
   close(unit=14)
   close(unit=12)
   call DAssignFFT()  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !
   Write(2,"(A,I3,2x,50(1x,I5))") 'Nr_UniqueStat=',Nr_UniqueStat, Unique_StatID(1:Nr_UniqueStat)
   write(2,"(A,I4,2x,50(1x,A5))") 'Nr_UniqueAnt=',Nr_UniqueAnt, (Statn_ID2Mnem(Unique_StatID(k)),k=1,Nr_UniqueStat)
    !
   Call Alloc_EInterfCalib_Pars()
   Call EIReadPeakInfo() ! get peak positions as they would be seen in a trace at the center of the core
   !Get fitoption:
   Call GetStationFitOption(FP_s, FitNoSources)
   write(2,*) 'reading for calibration finalized, FitNoSources=', FitNoSources
   !
   Fit_TimeOffsetStat(:)=0.
   Fit_TimeOffsetAnt(:)=0.
   Tot_UniqueAnt(0)=0
   StatID=Unique_SAI(1)/1000
   i_stat=1
   Do k=1,Nr_UniqueAnt
        If(Unique_SAI(k)/1000 .ne. StatID) then
            StatID=Unique_SAI(k)/1000
            i_stat=i_stat+1
        endif
        Tot_UniqueAnt(i_stat)=k       ! highest ranknr of unique_antenna for station Unique_StatID(i_stat)
   Enddo
   N_FitStatTim=0
   Nr_TimeOffset=1
   X_Offset(N_FitStatTim)=Nr_TimeOffset
   if(FP_s(1).GT.0) then  ! fit timeoffsets for this station
      Do i=1,N_FitPar_max
         if(FP_s(i) .le. 0) exit
         !If(.not. any(FP_s(i)==Ant_Stations(IntFer_ant(:,1),1)) ) cycle  ! check if this station is included in present fit for first chunk Nr_IntferCh
         !k=FINDLOC(Unique_StatID,FP_s(i),1)
         k=MAXLOC(Unique_StatID(:), MASK = Unique_StatID .eq. FP_s(i),DIM=1)
         If(k.le.0) cycle  !  ! Data from this station is not used in the present calculation
         If(Unique_StatID(k).ne.FP_s(i)) Then  ! Data from this station is not used in the present calculation
            write(2,*) 'i,k:',i,k,FP_s(i),Unique_StatID(:)
            stop 'problem in E_Callibr'
         Endif
         N_FitStatTim= N_FitStatTim + 1 ! position in 'FitParam'
         FitParam(N_FitStatTim)= k ! the timing of station "Unique_StatID(k)" will be fitted.
         X_Offset(N_FitStatTim)=Nr_TimeOffset   ! position in 'X' for this station, or the first antanna for this station
         If(Fit_AntOffset) then
              Do j= 1,Tot_UniqueAnt(k)-Tot_UniqueAnt(k-1)
                 X(Nr_TimeOffset)=Fit_TimeOffsetAnt(Tot_UniqueAnt(k-1)+j)
                 Nr_TimeOffset=Nr_TimeOffset+1
              Enddo
         Else
              X(Nr_TimeOffset)=Fit_TimeOffsetStat(k)
              Nr_TimeOffset=Nr_TimeOffset+1
         Endif
         X_Offset(N_FitStatTim)=Nr_TimeOffset  ! start location in X of the next batch of parameters
         !write(2,*) 'X:',i,FP_s(i),N_FitStatTim, X_Offset(N_FitStatTim),X(X_Offset(N_FitStatTim-1):X_Offset(N_FitStatTim)-1)
      Enddo
   endif
   !
   N_FitPar=N_FitStatTim
   X_Offset(N_FitStatTim)=Nr_TimeOffset  ! start location in X of the next batch of parameters
   !write(*,*) 'FP',FP(:),N_FitPar
   !write(2,"(A,40('(',i2,i3,')'))") '(FitParam,X_Offset)=',(FitParam(i),X_Offset(i),i=1,N_FitStatTim+1)
   If(.not.FitNoSources) Then
      Do i_Peak=1,PeakNrTotal
         k=X_Offset(N_FitPar)
            !write(2,*) 'from XIndx:', k,i,i_peak
          If(k .gt. N_FitPar_max) Then
            write(2,*) 'Too many fit parameters, ',k, &
               'exceeds max of', N_FitPar_max,' for i,i_peak=',i,i_Peak
            stop 'SetFitParamsLMA'
          Endif
          X( k:k+2 ) = SourcePos(1:3,i_Peak)
          N_FitPar= N_FitPar +1
          X_Offset(N_FitPar)=k+3
      enddo
   EndIf
    !
    write(2,*) '# of sources fitted ',N_FitPar-N_FitStatTim , &
       ', # of source parameters fitted=',X_Offset(N_FitPar)-X_Offset(N_FitStatTim)-1, &
       ', # of stations fitted',N_FitStatTim,' , # of station parameters fitted=',X_Offset(N_FitStatTim)-1
    flush(unit=2)
    !
    !stop 'FitCycle'
!    write(2,"(A,I4,A,I3,A,i4,A,F7.2,A)") 'N_FitPar=',N_FitPar ,', N_FitStatTim=', N_FitStatTim, &
!        ', max station nr in fit=',StatMax, ', max distance to ref. station=',DistMax,'[km]'
   Call EI_PrntFitPars(X)
    Call EI_Fitter(X)  ! fit source
    flush(unit=2)
    Call EIX2Source(X)
    !
   Call EI_PrntFitPars(X)
    !
   Write(2,"(//,A)") ' ==== Summary of new parameters ===='
   Call EIPrntNewSources
   ! Merge values for "Fit_TimeOffsetStat" with input from 'FineCalibrations.dat'
   If(WriteCalib) then
      Call WriteCalibration ! was MergeFine
   Endif
   !
   Return
    !
End Subroutine E_Callibr
!=================================
Subroutine EIReadPeakInfo()
   !
    use constants, only : dp,sample
    use DataConstants, only : PeakNr_dim, ChunkNr_dim
    use Chunk_AntInfo, only : Station_nrMax
    use ThisSource, only : SourcePos,  ExclStatNr, TotPeakNr !NrP, t_ccorr,
    use ThisSource, only : PeakNrTotal !,PeakNr,  PlotCCPhase, Safety, Nr_corr
    use ThisSource, only : PeakPos, ChunkNr
!    use FitParams
    use StationMnemonics, only : Statn_ID2Mnem, Station_Mnem2ID
    Implicit none
    !
    integer :: i_eo, i_chunk, i, k, i_c !,i_ca,i_eoa
    !logical :: Fitfirst !,TraceFirst
    Character(len=5) :: Station_Mnem, ExclStMnem(1:Station_nrMax)
    !Integer, parameter :: PeakS_dim=NrP
    integer :: i_Peak, nxx
    character*80 :: lname  ! Should be 80 to be compatible with "GetNonZeroLine"
    character*10 :: txt
    character*35 :: Frmt
    real(dp) :: x1,x2,x3,t
    !
   Call GetNonZeroLine(lname)
!      Write(2,"(i3,2i2,I8)", ADVANCE='NO') i_Peak,Peak_eo(i_Peak),ChunkNr(i_Peak),PeakPos(i_Peak)
!      Write(2,"(3(F10.2,','))", ADVANCE='NO') SourcePos(:,i_Peak)
!      Write(2,"(F8.2)", ADVANCE='NO') RefAntErr(i_Peak)
   Frmt="(i3,2x,i2,I8,3(F10.2,1x))"
   i_peak=0
   TotPeakNr(0,:)=0
   !PeakNr(:)=0
   i_chunk=1
   PeakNrTotal=PeakNr_dim
   Do i_Peak=1,PeakNr_dim       ! Read source positions from input. There should be at least one un-readable line in the input.
      !write(2,*) 'lname="',lname,'"'
      read(lname,Frmt,iostat=nxx)  i, i_c, k, x1,x2,x3
      !write(2,*) k, x1,x2,x3,t
      If(nxx.ne.0) Then
         Write(2,*) 'error when reading source info for #',i_Peak
         Write(2,*) 'Culprit: "',trim(lname),'"'
         Stop 'EI-sources read'
      EndIf
      SourcePos(1,i_Peak)=x1
      SourcePos(2,i_Peak)=x2
      SourcePos(3,i_Peak)=x3
      Peakpos(i_Peak)=k
      ChunkNr(i_Peak)=i_c
      if(i_chunk.gt.ChunkNr_dim .or. i_c.lt.i_chunk) Then
         Write(2,*) 'Chunk nr error when reading source info for #',i_Peak
         Write(2,*) 'Culprit: "',trim(lname),'"'
         Stop 'EI-sources chunk nr error'
      EndIf
      i_chunk=i_c
     TotPeakNr(0,i_chunk)=i_peak ! last peak# for this (i_eo,i_chunk); not really used and obsolete
     !PeakNr(i_chunk)=i_peak-PeakNr1        ! number of peaks for this (i_chunk)
     !
     ExclStatNr(:,i_peak)=0
     Call GetNonZeroLine(lname)
     !write(2,*) 'lname="',lname,'"'
     !read(lname,"(A7,i3,10i5)",iostat=nxx)  txt,k,ExclStatNr(1:k,i_peak)
     ExclStMnem='     '
     read(lname,*,iostat=nxx)  txt,ExclStMnem
     !write(2,*) 'excl',txt,';',ExclStMnem
     !If(nxx.ne.0 ) write(2,*) 'nxx=',nxx  ! always = -1
     If(trim(txt).eq.'exclude') Then
        Do k=1,Station_nrMax
            If(ExclStMnem(k).eq.'     ') exit
            Call Station_Mnem2ID(ExclStMnem(k),ExclStatNr(k,i_peak))
            If(ExclStatNr(k,i_peak).eq.0) exit
        Enddo
        k=k-1
        !read(lname,*,iostat=nxx)  txt,k,ExclStatNr(1:k,i_peak)  Statn_Mnem2ID
        !If(nxx.ne.0 .or. txt.ne.'exclude') cycle
        write(2,"(A,I3,A,I2,20(I4,A,A,',  '))") 'excluded for source',i_peak,' #',k, &
               (ExclStatNr(i,i_peak),'=',Statn_ID2Mnem(ExclStatNr(i,i_peak)), i=1,k)
        Call GetNonZeroLine(lname)
     EndIf
   Enddo
    write(2,*) 'PeakNrTotal=',PeakNrTotal  ! total # of pulses
    !write(2,*) 'TotPeakNr',TotPeakNr(0,:)         ! array giving number of peaks for this (i_eo,i_chunk)
    !
    !write(2,*) 'SourcePos=',SourcePos
    !
    !
    Return
End Subroutine EIReadPeakInfo
!================================
!=====================================
!============================
