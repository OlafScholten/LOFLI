!================================
!=================================
Subroutine FindCallibr(SourceGuess, XcorelationPlot)
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
   use Chunk_AntInfo, only : Unique_SAI, StartT_sam, TimeFrame, tot_uniqueant
   use DataConstants, only : Ant_nrMax
   use ThisSource, only : Dual, RealCorrelation
   use ThisSource, only : NrP, t_ccorr, SourcePos, RefAntErr, ExclStatNr, Peak_eo, ChunkNr, TotPeakNr
   use ThisSource, only : Nr_corr, PeakNr, PeakNrTotal, PlotCCPhase, Safety, PolBasis, Tref_dim, T2_Dim
   use ThisSource, only : PeakPos, Peak_eo, PeakRMS, PeakChiSQ, StStdDevMax_ns
   use ThisSource, only : XFrameEi, XFrameEf, XFrameNi, XFrameNf, XFrameh
   use FitParams, only : MaxFitAntDistcs, MaxFitAntD_nr, FitIncremental, PulsPosCore, WriteCalib, N_EffAnt, Max_EffAnt
   use FitParams, only : FullSourceSearch
   use constants, only : dp,sample
   use DataConstants, only : PeakNr_dim, ChunkNr_dim, RunMode
   use FitParams, only : MeanCircPol, Cnu01, Fit_TimeOffsetAnt, Fit_AntOffset, N_FitStatTim, FitParam, X_Offset
   use FFT, only : RFTransform_su, RFTransform_CF2CT
   !use StationMnemonics, only : Statn_ID2Mnem, Station_Mnem2ID
   Use Calibration, only : WriteCalibration ! was MergeFine
   Implicit none
   Real(dp), intent(in) :: SourceGuess(3,*)
   Logical, intent(in) :: XcorelationPlot
   !
   integer :: i, j, k, i_ant, j_corr, i_eo, i_chunk, i_dist, MinFitAntD_nr
   logical, save :: Fitfirst=.true.
   !Character(len=5) :: Station_Mnem, ExclStMnem(1:Station_nrMax)
   integer :: i_loc(1), i_Peak, StLoc, StatMax, ReadErr, PeakNr1
   Integer :: PeakSP(2*NrP,0:2), PeakSWl(2*NrP,0:2), PeakSWu(2*NrP,0:2), PeakSAmp(2*NrP,0:2), PeakD_nr
   Real(dp) :: DistMax
   Integer, parameter :: w_low=5 ! Min slices between peaks in odd&even to be merged; should be about twice the value set in DualPeakFind
   !character*80 :: lname  ! Should be 80 to be compatible with "GetNonZeroLine"
   character*10 :: txt
   character*35 :: Frmt
   real(dp) :: x1,x2,x3,t, ZeroCalOffset
   logical :: FitNoSources=.true.
   logical :: NewPeak
    real(dp) :: dt,B, MTC, dtI
    Complex(dp), Allocatable :: ST01nu(:)
    Complex(dp), Allocatable :: ST01T(:)
    Complex(dp) :: ST01_max
    Integer :: UpSF= 16, i_SAI, OutUnit ! UpSampleFactor
   !
   !       Initialize
   Call Find_unique_StatAnt()
   !
   Do i=-Safety,Safety
     t_ccorr(i)=i ! time axis needed for spline interpolation of CrossCorrelation
   Enddo
 !   Do i_Peak=1,PeakNr_dim      ! initialize
 !       SourcePos(:,i_Peak) = SourceGuess(:)
 !       RefAntErr(i_Peak) = 0.
 !   Enddo
!   write(2,*) 'sourceguess', NrP, SourceGuess
!   flush(unit=2)
   Nr_Corr=0
 !  ExclStatNr(:,:)=0
   Call GetRefAnt(ChunkNr_dim)
   !
   i_peak=0
   PeakNr1=0
  ! If(RunMode.eq.1) then ! To prevent reading from input in Explore
  !    ReadErr=-1
  ! Else
  !    If(TotPeakNr(1,ChunkNr_dim).eq.0)
  ! Endif
   !write(2,*) 'ReadErrh: ',ReadErr
   !flush(unit=2)
   write(2,*) 'FindCallibr:',ChunkNr_dim,TotPeakNr(1,ChunkNr_dim)
   If((RunMode.eq.1) .or. (TotPeakNr(1,ChunkNr_dim).eq.0)) then  ! Find peak positions instead of reading them from input
      Do i_chunk=1, ChunkNr_dim     ! get the NrP strongest peaks in this time block
         Call DualPeakFind(2*NrP, i_chunk, PeakD_nr, PeakSP, PeakSWl, PeakSWu, PeakSAmp) ! used for imaging production
         If(Dual) Then
            If((PeakNr1 + NrP) .gt. PeakNr_dim) then
               write(2,*) 'PeakNr_dim will be exceeded',PeakNr1, NrP
               stop 'FindCallibr-1'
            endif
            !Peakpos(PeakNr1+1:PeakNr1+NrP)=PeakSP(1:NrP,i_eo)
            Do j=1,NrP/2
                i_peak=PeakNr1+j
                Peakpos(i_peak)=PeakSP(j,0)
                Peak_eo(i_peak)=0
                ChunkNr(i_peak)=i_chunk
                SourcePos(:,i_Peak) = SourceGuess(:,i_chunk)
                RefAntErr(i_Peak) = 0.
                !write(2,*) 'Peakpos(i_peak)',i_peak,Peakpos(i_peak)
            Enddo
            Do i=1, NrP  ! add the positions from odd antennas
               !write(2,*) 'i',i
               NewPeak=.true.
               Do j=PeakNr1+1,PeakNr1 + NrP/2
                  If(abs(Peakpos(j)-PeakSP(i,1)) .lt. w_low) Then
                     Peakpos(j)=(Peakpos(j)+PeakSP(i,1))/2
                     NewPeak=.false.
                     !write(2,*) 'Peakpos(j)',j,Peakpos(j)
                     exit
                  EndIf
               EndDo
               If(NewPeak .and. (i_peak-PeakNr1).lt.NrP) Then
                   i_peak=i_peak+1
                   Peakpos(i_peak)=PeakSP(i,1)
                   Peak_eo(i_peak)=0
                   ChunkNr(i_peak)=i_chunk
                   SourcePos(:,i_Peak) = SourceGuess(:,i_chunk)
                   RefAntErr(i_Peak) = 0.
                   !write(2,*) 'Peakpos(i_peak);1',i_peak,Peakpos(i_peak)
               EndIf
            EndDo
            !write(2,*) 'PeakNr1=',PeakNr1,i_eo,i_chunk
            TotPeakNr(0,i_chunk)=i_peak ! last peak# for this (i_eo,i_chunk)
            PeakNr(0,i_chunk)=i_peak-PeakNr1        ! number of peaks for this (i_eo,i_chunk)
            Call DualReadPeakInfo(2,i_chunk)  !  Copy to the odd antennas
            PeakNr1=TotPeakNr(1,i_chunk)
         Else
            Do i_eo=0,1
               If((PeakNr1 + NrP) .gt. PeakNr_dim) then
                  write(2,*) 'PeakNr_dim will be exceeded',PeakNr1, NrP
                  stop 'FindCallibr-1'
               endif
               Peakpos(PeakNr1+1:PeakNr1+NrP)=PeakSP(1:NrP,i_eo)
               Do i_peak=PeakNr1+1,PeakNr1 + NrP
                   Peakpos(i_peak)=PeakSP(i_peak-PeakNr1,i_eo)
                   Peak_eo(i_peak)=i_eo
                   ChunkNr(i_peak)=i_chunk
                  SourcePos(:,i_Peak) = SourceGuess(:,i_chunk)
                  RefAntErr(i_Peak) = 0.
               Enddo
               PeakNr1=PeakNr1 + NrP
               !write(2,*) 'PeakNr1=',PeakNr1,i_eo,i_chunk
               TotPeakNr(i_eo,i_chunk)=PeakNr1! last peak# for this (i_eo,i_chunk)
               PeakNr(i_eo,i_chunk)=NrP        ! number of peaks for this (i_eo,i_chunk)
            enddo
         EndIf
      EndDo  ! i_chunk=1, ChunkNr_dim
      PeakNrTotal=PeakNr1
   endif
   write(2,*) 'PeakNrTotal:',PeakNrTotal
   write(2,*) 'Peakpos:',Peakpos
   write(2,*) 'Peak_eo:',Peak_eo
   write(2,*) 'ChunkNr:',ChunkNr
   write(2,"(1x,A,100F12.3)") 'SourcePos1:',SourcePos(1,1:PeakNrTotal)
   write(2,"(1x,A,100F12.3)") 'SourcePos2:',SourcePos(2,1:PeakNrTotal)
   write(2,"(1x,A,100F12.3)") 'SourcePos3:',SourcePos(3,1:PeakNrTotal)
   write(2,*) TotPeakNr
   write(2,*) PeakNr
   flush(unit=2)
   !
   !write(2,*) 'FullSourceSearch: ', FullSourceSearch
   !flush(unit=2)
   FitNoSources=.not.FullSourceSearch
   PulsPosCore=.false.
   If(FitIncremental) then
      MinFitAntD_nr=1
   Else
      MinFitAntD_nr=MaxFitAntD_nr
   Endif
   StatMax=1050 ! Ant_Stations(1,1)
   !FitNoSources=.false.
   PeakRMS(1)=0.
   Do i_dist=MinFitAntD_nr,MaxFitAntD_nr
      DistMax=MaxFitAntDistcs(i_dist) ! in kilometers
      If(DistMax.eq.-1) exit
    !write(2,*) 'RefAntErr---X', RefAntErr(1)
      Call FitCycle(FitFirst,StatMax,DistMax,FitNoSources)
    !write(2,*) 'RefAntErr---Y', RefAntErr(1)
      !Call PrntNewSources()
   !write(2,*) 'i_dist: ', i_dist
   !flush(unit=2)
      !Write(2,*) 'Average RMS=',SUM(PeakRMS(1:PeakNrTotal))/PeakNrTotal
      If((RunMode.eq.1) .and. (MINVAL(PeakRMS(1:PeakNrTotal)).gt.30.)) then
         write(2,"(A,F7.2,A,10F6.0)") 'chi-square too poor for max antenna distance=',DistMax &
            ,'[m], PeakRMS(1:PeakNrTotal):',PeakRMS(1:PeakNrTotal)
         Return
      Endif
   EndDo ! i_dist=1,MaxFitAntD_nr
   !flush(unit=2)
    !
    If(RunMode.eq.1) then  ! =Explore
      !write(2,*) 'PeakNrTotal=',PeakNrTotal
      Do i_peak=1,PeakNrTotal
         If(PeakRMS(i_Peak).gt. 25.) cycle
         write(2,"(i2,i7,', N,E,h=', 3F9.0,', RMS=',F7.1)") i_peak,Peakpos(i_peak),SourcePos(:,i_Peak),PeakRMS(i_Peak)
         If((SourcePos(2,i_Peak)/1000. .lt. XFrameEi) .or. (SourcePos(2,i_Peak)/1000. .gt. XFrameEf)) cycle
         If((SourcePos(1,i_Peak)/1000. .lt. XFrameNi) .or. (SourcePos(1,i_Peak)/1000. .gt. XFrameNf)) cycle
         If(SourcePos(1,i_Peak) .ne. SourcePos(1,i_Peak)) cycle  ! check for NaN
         If(abs(SourcePos(3,i_Peak))/1000. .gt. XFrameh) cycle
         write(29,"(1x,i5,4(2x,g14.8),3x,f9.6,2I5)")  &
            TimeFrame*10+i_Peak,StartT_sam(1)*1000.d0*sample, &
            SourcePos(2,i_Peak)/1000.,SourcePos(1,i_Peak)/1000.,abs(SourcePos(3,i_Peak)/1000.),PeakRMS(i_Peak), &
            N_EffAnt, Max_EffAnt
      Enddo
      Return
    endif
!   flush(unit=2)
    StatMax=1050
    DistMax= MaxFitAntDistcs(MaxFitAntD_nr) ! in kilometers
    !FitNoSources=.false.
    !FitNoSources=.true.
    Call FitCycle(FitFirst,StatMax,DistMax,FitNoSources)
    !
1   continue
    !
    !If(Polariz) Then
    !  Allocate( Cnu01(0:T2_Dim/2,0:Ant_nrMax/2) )
    !  Cnu01(:,:)=0.d0
    !  MeanCircPol=.true.
    !EndIf
    MeanCircPol=.false.
    PlotCCPhase=XcorelationPlot
    Call BuildCC(StatMax,DistMax)
   !
   If(MeanCircPol) Then  ! Obsolete option
      call RFTransform_su(UpSF*T2_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      Allocate( ST01T(1:UpSF*T2_dim) )
      Allocate( ST01nu(0:UpSF*T2_dim/2) )
      !    Write(2,*) 'Mean t-p dt:',SUM(dtAnt_tp(:))/SUM(NSrc_tp(:)),'Polariz',Polariz, T2_dim, Tref_dim
      Do i=1, N_FitStatTim
         k = FitParam(i)   ! =Station number
      !   If(Fit_AntOffset .and. Polariz) then ! The other cases are dealt with in X2Source
         If(Fit_AntOffset) then ! The other cases are dealt with in X2Source
         ! Here we want to take care of the timing difference between odd and even antennas only
            Do j= 1,(Tot_UniqueAnt(k)-Tot_UniqueAnt(k-1))/2  ! Keep antenna delay difference for even and odd
               i_SAI=Tot_UniqueAnt(k-1)+2*j-1
               If(MOD(Unique_SAI(i_SAI),2).eq.0) Then
                  ST01nu(0:T2_dim/2)=Cnu01(0:T2_dim/2,i_SAI/2)  ! as generated in " Call BuildCC(StatMax,DistMax)"
                  ST01nu(T2_dim/2+1:UpSF*T2_dim/2)=0.
                  Call RFTransform_CF2CT(ST01nu,ST01T )
                  i_loc=MaxLoc(REAL(ST01T(1+UpSF*Tref_dim/2:UpSF*T2_dim-UpSF*Tref_dim/2) )) +UpSF*Tref_dim/2
                  B=(REAL(ST01T(i_loc(1)+1))-REAL(ST01T(i_loc(1)-1)))/2.   ! fit with Y(x)=Y(0)+B x - MTC/2 x^2
                  MTC=2*REAL(ST01T(i_loc(1))) - REAL(ST01T(i_loc(1)+1)) - REAL(ST01T(i_loc(1)-1))
                  dt=B/MTC  ! additional odd-even timeing difference from circ pol
                  ST01_max= ST01T(i_loc(1))*(1-dt*dt) +  &
                     dt*( ST01T(i_loc(1)+1)*(dt+1) + ST01T(i_loc(1)-1)*(dt-1) )/2.
                  B=(IMAG(ST01T(i_loc(1)+1))-IMAG(ST01T(i_loc(1)-1)))/2.   ! fit with Y(x)=Y(0)+B x - MTC/2 x^2
                  MTC=2*IMAG(ST01T(i_loc(1))) - IMAG(ST01T(i_loc(1)+1)) - IMAG(ST01T(i_loc(1)-1))
                  dtI=sqrt(B*B+2*IMAG(ST01T(i_loc(1)))*MTC)
                  DTI=(B-sign(DTI,B))/MTC
                  !Write(2,"(A,F6.2,A,F6.2,3(';',2F8.3))") 'Real-Max@',dt,', Imag=0@',dtI,ST01T(i_loc(1)-1:i_loc(1)+1) ! fit with Y(x)=Y(0)+B x - MTC/2 x^2
                  dt=(i_loc(1)+dt)/UpSF -T2_dim/2.
                  write(2,"('antennas=',I7,',',I7,' at dt(Max(U))=',f6.2,'[samples] has Stokes U+iV=(',2g10.4,')')", &
                        ADVANCE='NO') Unique_SAI(i_SAI),Unique_SAI(i_SAI+1),dt,ST01_max
                  !i_loc=MinLoc(REAL(ST01T(1+UpSF*Tref_dim/2:UpSF*T2_dim-UpSF*Tref_dim/2)  ))+UpSF*Tref_dim/2
                  dtI=(i_loc(1)+dtI)/UpSF -T2_dim/2.
                  write(2,"(', and at dt=',f6.2,' has Stokes V=0')") dtI
                  !i_loc=MaxLoc(Abs(ST01T(1+UpSF*Tref_dim/2:UpSF*T2_dim-UpSF*Tref_dim/2)  ))+UpSF*Tref_dim/2
                  !write(2,"(', and at dt(Max|U+iV|)=',f6.2,' has U+iV Stokes',2g10.4)") &
                  !                                 1.*i_loc(1)/UpSF-T2_dim/2.,ST01T(i_loc(1))
                  !write(2,*) 'dt=' , Fit_TimeOffsetAnt(i_SAI),Fit_TimeOffsetAnt(i_SAI+1),dt !  = 0 !!!!!!!!!!!!!!!!               !
                  !dt=Fit_TimeOffsetAnt(i_SAI)-Fit_TimeOffsetAnt(i_SAI+1) - dt
                  Fit_TimeOffsetAnt(i_SAI)  = Fit_TimeOffsetAnt(i_SAI)  - dt/2.
                  Fit_TimeOffsetAnt(i_SAI+1)= Fit_TimeOffsetAnt(i_SAI+1) +dt/2.
                  !write(2,*) 'FindCallibr', i,j,Unique_SAI(i_SAI),Fit_TimeOffsetAnt(i_SAI)&
                  !      ,Unique_SAI(i_SAI+1), Fit_TimeOffsetAnt(i_SAI+1)
               Else
                  write(2,*) 'This SAI is not allowed for this case', Unique_SAI(i_SAI),'should be even'
               EndIf
            Enddo
         EndIf
         !Stop
      Enddo
      DeAllocate( ST01T )
      DeAllocate( ST01nu )
      DeAllocate( PolBasis)
      DeAllocate( Cnu01 )
      MeanCircPol=.false.
   Endif
   Write(2,"(//,A,F5.1)") ' ==== Summary of new parameters, stations will be dropped when timing-StDev exceeds ', StStdDevMax_ns
   ! Merge values for "Fit_TimeOffsetStat" with input from 'FineCalibrations.dat'
   ZeroCalOffset=0.
   If(WriteCalib) then
      Call WriteCalibration(ZeroCalOffset) ! was MergeFine
   Endif
   OutUnit=2
   Call PrntNewSources(ZeroCalOffset, OutUnit)
   !
   PlotCCPhase=.false.
   Return
    !
End Subroutine FindCallibr
!=================================
!================================
Subroutine FitCycle(FitFirst,StatMax,DistMax,FitNoSources)
    !use DataConstants, only : ChunkNr_dim
    !use Chunk_AntInfo, only : Ant_Stations, Ant_IDs, Ant_nr, Ant_pos
    use ThisSource, only : Nr_Corr, PeakRMS, PeakChiSQ, PeakNrTotal
    use FitParams, only : N_FitPar, N_FitPar_max, N_FitStatTim, CalcHessian, FitQual, SpaceCov
    use constants, only : dp
    Implicit none
    logical, intent(inout) :: FitFirst
    logical, intent(in) :: FitNoSources
    integer, intent(in) :: StatMax
    Real(dp), intent(in) :: DistMax
    real ( kind = 8 ) :: X(N_FitPar_max) !, Dist
    integer :: FitPos(4) !, i, j, k, i_ant, j_corr, i_eo, i_chunk, Station
    !logical :: FitNoSources=.false.
    !
    !
    !write(2,*) '=== New Round ============================================================'
    flush(unit=2)
    CalcHessian=.false.
    SpaceCov(1,1)=-.1 ; SpaceCov(2,2)=-.1 ; SpaceCov(3,3)=-.1 ;
    Call BuildCC(StatMax,DistMax)
    !
    FitPos(1)=1 ; FitPos(2)=2 ; FitPos(3)=3 ; FitPos(4)=-4  !(4=timing error, 1-3 is N,E,h)
    If(DistMax.le. 0.30) then
        FitPos(1)=1 ; FitPos(2)=2 ; FitPos(3)=-1 ; FitPos(4)=-1 ! single station, fit direction only
    ElseIf(DistMax.le. 1.0) then
        FitPos(1)=1 ; FitPos(2)=2 ; FitPos(3)=3 ; FitPos(4)=-1 ! Superterp, fit direction and curvature front
    ElseIf(DistMax.le. 5) then
        FitPos(1)=1 ; FitPos(2)=2 ; FitPos(3)=3 ; FitPos(4)=-4
    EndIf
    If(FitNoSources) then
        FitPos(1)=-1 ; FitPos(2)=-2 ; FitPos(3)=-3 ; FitPos(4)=-4
    endif
    Call SetFitParamsLMA(X,Fitfirst,FitPos) ! Important to call 'SetFitParamsLMA' after 'GetCorrSingAnt'
    !write(2,*) 'N_FitPar1 , N_FitStatTim',N_FitPar , N_FitStatTim
    !flush(unit=2)
    !
    !stop 'FitCycle'
    !write(2,*) 'Fit with Nr_Corr=',Nr_Corr,'============================================================'
    !write(2,"(A,I4,A,I3,A,i4,A,F7.2,A)") 'N_FitPar=',N_FitPar ,', N_FitStatTim=', N_FitStatTim, &
    !    ', max station nr in fit=',StatMax, ', max distance to ref. station=',DistMax,'[km]'
   If(StatMax.lt.50) then
      write(2,"(A,2I4,A,I4,A,i4,A,F7.2,A)") '==== Fit with Nr_Corr=',Nr_Corr(0:1,1),', N_FitPar=',N_FitPar , &
         ', max station nr in fit=',StatMax, ', max distance to ref. station=',DistMax,'[km]'
   else
      write(2,"(A,2I4,A,I4,A,F7.2,A)") '==== Fit with Nr_Corr=',Nr_Corr(0:1,1),', N_FitPar=',N_FitPar, &
         ', max distance to ref. station=',DistMax,'[km]'
   endif
!    write(2,"(A,I4,A,I3,A,i4,A,F7.2,A)") 'N_FitPar=',N_FitPar ,', N_FitStatTim=', N_FitStatTim, &
!        ', max station nr in fit=',StatMax, ', max distance to ref. station=',DistMax,'[km]'
    Call FitCCorr(X)  ! fit source
    write(2,*) 'chi-square=',FitQual, ', Average RMS=',SUM(PeakRMS(1:PeakNrTotal))/PeakNrTotal,&
      ', per source:', PeakRMS(1:PeakNrTotal)
    flush(unit=2)
    Call X2Source(X)
    !
End Subroutine FitCycle
!=====================================
Subroutine GetLargePeaks(i_ant, i_chunk, Peakposs)  ! Obsolete, not used anymore
    use DataConstants, only : Time_dim
    use Chunk_AntInfo, only : CTime_spectr, Ant_IDs
    use ThisSource, only : NrP  ! PeakNr
    use constants, only : dp
    Implicit none
    Integer, intent(in) :: i_ant, i_chunk
    Integer, intent(out) :: Peakposs(1:NrP)
    !
    real(dp) :: HEnvel(Time_dim), peakStr(1:9), AvePos, TotAmpl, StDev
    integer :: i, j, k, PeakN,Offst
    integer :: i_loc(1), PeakPos(0:9), NPeakPos(0:9)
    !
    HEnvel(:)=abs(CTime_spectr(:,i_ant, i_chunk))
    PeakN=1
    PeakPos(0)=0
    PeakPos(PeakN)=Time_dim
    Offst=103
    Do i=1,3 ! get 3 largest peaks
        k=0
        NPeakPos(0)= 0
        Do j=1,PeakN
            If( (PeakPos(j)-PeakPos(j-1)) .gt. 300) then
                i_loc=MaxLoc( HEnvel(PeakPos(j-1)+Offst:PeakPos(j)-Offst) )
                !A= Maxval( HEnvel(PeakPos(j-1)+Offst:PeakPos(j)-Offst) )
                k=k+1
                NPeakPos(k) = PeakPos(j-1)+Offst-1 + i_loc(1)
                !write(2,*) 'i_loc',i_loc(1),HEnvel(NPeakPos(k)),A
                !write(2,*) HEnvel(PeakPos(j-1)+Offst+i_loc(1)-5:PeakPos(j-1)+Offst+i_loc(1)+5)
                !write(2,*) 'NPeakPos(k)',NPeakPos(0:k)
            endif
            k=k+1
            NPeakPos(k) = PeakPos(j)
        enddo
        PeakPos(:)=NPeakPos(:)
        PeakN=k
        write(2,*) i,k,'PeakPos(:)',PeakPos(0:k)
    enddo ! i
    !
    peakStr=0
    Do j=1,PeakN-1
        Call Inv5StarPk(HEnvel, PeakPos(j), Offst, AvePos, TotAmpl, StDev)
        write(2,*) j,'position offset=',PeakPos(j)-AvePos,', StDev=',StDev, &
         ', power=',StDev*peakStr(j), ' or', TotAmpl
    enddo
    !Write(2,*) 'PeakPos',peakpos
    !Write(2,*) 'PeakStr',peakstr
    !
    Do i=1,NrP  ! get the 4 strongest peaks
        i_loc=MaxLoc( peakStr(1:PeakN-1) )
        PeakPoss(i) = PeakPos( i_loc(1) ) ! The peak position for this as well as a vistual antenna in the center of CS002
        ! The latter is true because of the way the RawSourceDistance has been calculated and corrected for by shifting the spectrum
        peakStr(i_loc(1))=0.
    Enddo
    !
    Write(2,"('data block#',i2,', Peak position=', 5I9)") i_chunk, PeakPoss(1:NrP)
    Write(2,"('   Ant-ID=',i3,', Value Envelope=', 5F9.1)") Ant_IDs(i_ant,i_chunk),HEnvel(PeakPoss(1:NrP))
    !
    Return
End Subroutine GetLargePeaks
!============================
!----------------------------------
Subroutine GetStationFitOption(FP_s, FitNoSources)
   use DataConstants, only : Station_nrMax
   use FitParams, only : Fit_AntOffset, N_FitPar_max
   use Chunk_AntInfo, only : Nr_UniqueStat, Unique_StatID
   use StationMnemonics, only : Station_Mnem2ID
   Implicit none
   integer, Intent(OUT)::  FP_s(0:N_FitPar_max)
   Logical, Intent(out) ::  FitNoSources
   Character(len=5) :: FP_MNem(0:N_FitPar_max)
   logical,dimension(Nr_UniqueStat) :: mask
   character*300 :: lname
   character*4 :: option
   Integer :: i, nxx, k, N_FitPar
   !
   FP_s(:)=0
   FP_MNem(:)=' '
   FitNoSources=.false.
   Fit_AntOffset=.false.
   N_FitPar=0
   Call GetNonZeroLine(lname)
   write(2,*) 'timing offset-fit option="',lname,'"'
   write(2,*) 'possible options are: "only" & "abut" & "antenna", followed by station names and/or "NoSrc", ending with "!" '
   read(lname,*,iostat=nxx) option, FP_MNem(1:N_FitPar_max) ! option=abut,only
   Do i=1,N_FitPar_max
      !write(2,*) 'FP_Mnem(i):',i,FP_Mnem(i)
      If(FP_Mnem(i).eq.'     ') exit
      If(TRIM(FP_Mnem(i)).eq.'!') exit
      If(FP_Mnem(i).eq.'NoSrc') Then
          FitNoSources=.true.
          cycle
      EndIf
      Call Station_Mnem2ID(FP_Mnem(i),k)
      !write(2,*) 'FP_Mnem(i)2:',i, FP_Mnem(i),k
      If(k.eq.0) cycle ! exit
      N_FitPar=N_FitPar+1
      FP_s(N_FitPar)=k
      !write(2,*) i,FP_Mnem(i),FP_s(i)
   Enddo
   !read(lname,*,iostat=nxx) option, FP_s(1:N_FitPar_max) ! option=abut,only
   If(option .eq. 'abut') then
      mask=.true.
      Do i=1,Nr_UniqueStat
          if ( any(Unique_StatID(i)==FP_s) ) mask(i) = .false.
          !write(2,*) i,Unique_StatID(i),mask(i)
      enddo
      i=count(mask)
      !write(2,*) N_FitPar_max, Nr_UniqueStat, size(Unique_StatID), Station_nrMax,i
      !Flush(unit=2)
      !write(2,*) 'mask=',mask(1:size(Station_IDs))
      !write(2,*) 'Station_IDs=',Station_IDs(1:size(Station_IDs))
      FP_s(1:i)=pack(Unique_StatID(1:Nr_UniqueStat), mask(1:Nr_UniqueStat) )
      !write(2,*) 'FP_s',FP_s(0:Nr_UniqueStat)
      !Flush(unit=2)
      N_FitPar=i
   ElseIf(option.eq. 'ante') Then
      Fit_AntOffset=.true.
      write(2,*) 'fit calibration for antennas in stations separately'
   endif
   !
   Write(2,"(A,I3,A)") "Fitting timings for",N_FitPar," stations"
   If(FitNoSources) Write(2,*) 'No source parameters are being fitted'
   Return
   !
End Subroutine GetStationFitOption
