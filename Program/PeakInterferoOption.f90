!=================================
!================================
!=================================
Subroutine PeakInterferoOption()
!   Chi-square interferometric imaging of peaks found by the impulsive imager
!
!   Procedural steps:
!   1) Read sources selected by FlashImage
!   2) Use a similar fit procedure as used for FCalibrate
!
   use constants, only : dp,sample
   use DataConstants, only : PeakNr_dim, ChunkNr_dim, RunMode, Time_Dim
   use Chunk_AntInfo, only : Unique_SAI, StartT_sam, Tot_UniqueAnt, Ant_Stations
   use Chunk_AntInfo, only : Unique_StatID,  Nr_UniqueStat, Nr_UniqueAnt, N_Chunk_max
   use DataConstants, only : Ant_nrMax
   use ThisSource, only : SourcePos
   use ThisSource, only : PeakNr, PeakNrTotal
   use ThisSource, only : ChunkNr
   use FitParams, only : N_FitPar, N_FitPar_max, Nr_TimeOffset, WriteCalib
   use FitParams, only : Fit_TimeOffsetStat, Fit_TimeOffsetAnt, Fit_AntOffset, N_FitStatTim, FitParam, X_Offset
   Use Interferom_Pars, only : IntFer_ant,  Nr_IntferCh ! the latter gives # per chunk
   use Interferom_Pars, only : Alloc_EInterfCalib_Pars
   use StationMnemonics, only : Statn_ID2Mnem !, Station_Mnem2ID
    use LOFLI_Input, only : ReadSourceTimeLoc
   use FFT, only : RFTransform_su,DAssignFFT, RFTransform_CF2CT
   Use Calibration, only : WriteCalibration ! was MergeFine
   Implicit none
   Real(dp) :: SourceGuess(3,N_Chunk_max)
   integer, save ::  FP_s(0:N_FitPar_max)!, FP(0:4)
   real ( kind = 8 ) :: X(N_FitPar_max)
   !
   integer :: i, j, k, i_ant, i_eo, i_chunk, i_stat, StatID, nxx
   integer :: i_loc(1), i_Peak
   Logical :: FitNoSources
   character*10 :: StochOption
   character*80 :: lname
   real(dp) :: x1,x2,x3,t, Time_width, Space_Spread(1:3)
   !real(dp) :: dt,B, MTC, dtI
   Integer :: i_SAI
   !
   Call ReadFlashImageDat(SourceGuess)
   !       Read appropriate data chunks
   Call RFTransform_su(Time_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   write(*,"(A,i3,A)") achar(27)//'[45m # of data-blocks read in=',ChunkNr_dim,achar(27)//'[0m'
   Do i_chunk=1, ChunkNr_dim
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
   !Call EIReadPeakInfo() ! get peak positions as they would be seen in a trace at the center of the core
   !Get fitoption:
   !Call GetStationFitOption(FP_s, FitNoSources)
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
   N_FitPar=N_FitStatTim
   X_Offset(N_FitStatTim)=Nr_TimeOffset  ! start location in X of the next batch of parameters
   !write(*,*) 'FP',FP(:),N_FitPar
   !write(2,"(A,40('(',i2,i3,')'))") '(FitParam,X_Offset)=',(FitParam(i),X_Offset(i),i=1,N_FitStatTim+1)
   Do i_Peak=1,PeakNrTotal
      k=X_Offset(N_FitPar)
         !write(2,*) 'from XIndx:', k,i,i_peak
       If(k .gt. N_FitPar_max) Then
         write(2,*) 'Too many fit parameters, ',k, &
            'exceeds max of', N_FitPar_max,' for i,i_peak=',i,i_Peak
         stop 'ECallibrate: too many parameters'
       Endif
       X( k:k+2 ) = SourcePos(1:3,i_Peak)
       N_FitPar= N_FitPar +1
       X_Offset(N_FitPar)=k+3
   enddo
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
    Time_width=-1.
    Call EI_Fitter(X, Time_width, Space_Spread)  ! fit source
    flush(unit=2)
    Call EIX2Source(X)
    !
   Call EI_PrntFitPars(X)
    !
   Write(2,"(//,A)") ' ==== Summary of new parameters ===='
   Call PrntNewSources
   !
   !
   Return
    !
End Subroutine PeakInterferoOption
!=================================
Subroutine WriteInterfRslts(i_Peak)
   use Constants, only : dp,sample,c_mps, pi
   use DataConstants, only : DataFolder, OutFileLabel
   Use Interferom_Pars, only : StI, StI12, StQ, StU, StV, StI3, StU1, StV1, StU2, StV2, P_un, P_lin, P_circ
   Use Interferom_Pars, only : dStI, dStI12, dStQ, dStU, dStV, dStI3, dStU1, dStV1, dStU2, dStV2, Chi2pDF, BoundingBox, StartTime_ms
   use ThisSource, only : PeakNrTotal, ChunkNr, PeakChiSQ, PeakRMS, PeakPos, SourcePos
   use Chunk_AntInfo, only : StartT_sam
   use GLEplots, only : GLEplotControl
   Implicit none
   Integer, intent(in) :: i_Peak  ! =29, 28
   Real(dp) :: time, x, y, z
   Real(dp), external :: tShift_ms   !
   !
   !Time=SQRT(sum(SourcePos(:,i_Peak)*SourcePos(:,i_Peak)))
   !Time = Time*Refrac/(c_mps*sample) ! Convert to units of [samples]
   !Time=( StartT_sam(ChunkNr(i_Peak)) + PeakPos(i_Peak) - Time )*Sample*1000. - StartTime_ms
   Time=( StartT_sam(ChunkNr(i_Peak)) + PeakPos(i_Peak) )*Sample*1000.d0 - tShift_ms(SourcePos(:,i_Peak)) - StartTime_ms
   x=SourcePos(2,i_Peak)/1000.
   y=SourcePos(1,i_Peak)/1000.
   z=SourcePos(3,i_Peak)/1000.
   !
   If(i_peak.eq.1) then
      OPEN(UNIT=28,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'IntfPkPolrz.dat')
      write(28,"(2F11.3,4F10.2,1x,A,I3,F7.1,I3,' 0')") Time-.0005, Time+1., &
         0.0 ,  SourcePos(:,i_Peak) ,TRIM(OutFileLabel), 0 , 0.0, PeakNrTotal  ! gleRead:  xMin xMax yMin yMax zMin zMax tMin tMax ZmBx$ t_offst
      Write(28,"(A,9x,A,13x,';   ', 9(A,'/I [%]   ;   '),A )") &
         '!   nr,  time[ms] ;','StI','StI12', '  StQ','  StU','  StV',' StI3',' StU1',' StV1',' StU2',' StV2', 'Chi^2/DoF'
      OPEN(UNIT=29,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'IntfPkPow.dat')
      write(29,"(6F8.2,2F9.3,A,F7.1,i3,' 0')") BoundingBox(1:2,1:4), ' NoBox ', StartTime_ms, 0 ! gleRead:  xMin xMax yMin yMax zMin zMax tMin tMax ZmBx$ t_offst
      Write(29,*) '0 ',PeakNrTotal,' 0 0 ',TRIM(OutFileLabel), 10.,' 0 0 0 0 0 0 0 ',0.1 ,  '1.0 !' ! gleRead:  NTracks EventNr Q t_start label$ AmplitudePlot a1 b1 c1 a2 b2 c2 d2
   EndIf
   !
   write(29,"(i6,' ',4(f11.5,' '),g13.6,' ',3f6.2,' ',5g13.3)") i_Peak, Time, &
      x, y, z, StI12, StQ/StI, StU/StI, StV/StI, StI3 ,P_lin,( ATAN2(StU,StQ)/2. )*180./pi
   !
   Write(28,"(i6,' ',f11.5,' ', 2(g12.4,' '), 22(f8.3,' ') )") i_Peak, Time, &
      StI, dStI, 100*StI12/StI, 100*dStI12/StI, 100*StQ/StI, 100*dStQ/StI, &
      100*StU/StI, 100*dStU/StI, 100*StV/StI, 100*dStV/StI, &
      100*StI3/StI, 100*dStI3/StI, 100*StU1/StI, 100*dStU1/StI, 100*StV1/StI, 100*dStV1/StI, &
      100*StU2/StI, 100*dStU2/StI, 100*StV2/StI , 100*dStV2/StI, &
      Chi2pDF, 100.*P_un, 100.*P_lin, 100.*P_circ
   !
   If(i_Peak.eq.PeakNrTotal) Then
      Close(Unit=28)
      Close(Unit=29)
      Call GLEplotControl(PlotType='SourcesPlot', PlotName='IntfPk'//TRIM(OutFileLabel), &
         PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel)//'IntfPkPow' )
      Call GLEplotControl(PlotType='EIPolariz', PlotName='IntfPkPol'//TRIM(OutFileLabel), &
         PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel)//'IntfPkPolrz' )
   EndIf
   !
   Return
End Subroutine WriteInterfRslts
!=========================
Subroutine ReadFlashImageDat(SourceGuess)
   !
    use constants, only : dp,sample, c_mps
    use DataConstants, only : PeakNr_dim, ChunkNr_dim, DataFolder, EdgeOffset, Time_dim
    use Chunk_AntInfo, only : Station_nrMax, StartT_sam
    Use Interferom_Pars, only : BoundingBox, StartTime_ms
    use ThisSource, only : SourcePos,  TotPeakNr !NrP, t_ccorr,
    use ThisSource, only : PeakNrTotal !,PeakNr,  PlotCCPhase, Safety, Nr_corr
    use ThisSource, only : PeakPos, ChunkNr
!    use FitParams
    use StationMnemonics, only : Statn_ID2Mnem, Station_Mnem2ID
    Implicit none
    !
    Real(dp), intent(out) :: SourceGuess(1:3,*)
    integer :: i_eo, i_chunk, i, k!, i_c ,i_ca  !,i_eoa
    Character(len=5) :: Station_Mnem, ExclStMnem(1:Station_nrMax)
    integer :: i_Peak, nxx, TimeFrame, TimeFrame0
    character*80 :: lname
    character*50 :: FlashFile
    character*10 :: txt
    character*35 :: Frmt
    real(dp) :: x1,x2,x3,t!, TimeShft
    Real(dp), external :: tShift_smpl   !
    !
   Call GetNonZeroLine(lname)
   Read(lname,*) FlashFile
   Open(Unit=17, STATUS='old',ACTION='read', FILE=TRIM(datafolder)//TRIM(FlashFile), IOSTAT=nxx)
   If(nxx.ne.0) Then
      Write(2,*) 'Not found!, source information file: "',TRIM(datafolder)//TRIM(FlashFile),'"'
      Stop 'SourcesDataFile not found'
   Else
      Write(2,*) 'Source information read from file: "',TRIM(datafolder)//TRIM(FlashFile),'"'
   EndIf
!  -20.830  -20.770  -16.280  -16.220   10.700   10.900   10.500   35.000 NoBox  300.000  !
!   0      26   3.50  0.3000  21C-1/21C1e-E2aZ3   0.00   0.556      672.  84.021   6.633181725864.  39.148   0.000  !
!    53003   11.492543      -20.807222      -16.259872       10.738346          3.34    6   10

   TimeFrame0=-9999
   i_chunk=0
   i_peak=0
   TotPeakNr(0,:)=0
   Read(17,*,iostat=nxx) BoundingBox(1:2,1:4),txt,StartTime_ms
   write(2,*) 'StartTime_ms', StartTime_ms
   Read(17,*) txt
   Do while (i_peak.lt.PeakNr_dim)
      Read(17,*,iostat=nxx) k,t,x2,x1,x3
      If(nxx.ne.0) Then
         exit
      EndIf
      TimeFrame=k/1000
      If(TimeFrame.ne.TimeFrame0) then
         i_chunk=i_chunk+1
         if(i_chunk.gt.ChunkNr_dim) Then
            Write(2,*) 'Chunk nr error when reading source info for #',i_Peak
            Write(2,*) 'Culprit: "',trim(lname),'"'
            Stop 'EI-sources chunk nr error'
         EndIf
         TimeFrame0=TimeFrame
         StartT_sam(i_chunk)=(StartTime_ms/1000.d0)/sample + (TimeFrame-1)*(Time_dim-2*EdgeOffset)  ! in sample's
         SourceGuess(1,i_chunk)=x1*1000.
         SourceGuess(2,i_chunk)=x2*1000.
         SourceGuess(3,i_chunk)=x3*1000.
      EndIf
      i_peak=i_peak+1
      SourcePos(1,i_Peak)=x1*1000.
      SourcePos(2,i_Peak)=x2*1000.
      SourcePos(3,i_Peak)=x3*1000.
      ChunkNr(i_Peak)=i_chunk
      TotPeakNr(0,i_chunk)=i_peak ! last peak# for this (i_eo,i_chunk); not really used and obsolete
      !TimeShft=SQRT(sum(SourcePos(:,i_Peak)*SourcePos(:,i_Peak)))
      !TimeShft = TimeShft*Refrac/(c_mps*sample) ! Convert to units of [samples]
      !PeakPos(i_Peak) = (StartTime_ms+t)/(Sample*1000.) - StartT_sam(ChunkNr(i_Peak))+ TimeShft
      PeakPos(i_Peak) = (StartTime_ms+t)/(Sample*1000.d0) - StartT_sam(ChunkNr(i_Peak))+ tShift_smpl(SourcePos(:,i_Peak))
      write(2,*) k,i_peak,i_chunk, StartT_sam(i_chunk), x1, x2, x3, PeakPos(i_Peak)
      !write(2,*) t/(Sample*1000.), StartT_sam(ChunkNr(i_Peak)), TimeShft

   EndDo
   PeakNrTotal=i_peak
   ChunkNr_dim=i_chunk
   Close(Unit=17)
   !
   !
   Return
End Subroutine ReadFlashImageDat
!================================
!=====================================
!============================
