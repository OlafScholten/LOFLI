!=================================
Subroutine SourceFind(TimeFrame,SourceGuess,units)
!   Find source locations for peaks in a time-chunk
!   This should be done for even as well as odd numbered antennas.
!
!   Procedural steps:
!  1) find peaks in reference antenna
!  2) Follow peaks in other antennas at increasing distances
!
!   For fitting the Abs(corr) or Real(corr) can be mixed.
!  V13: Introduce Kalman filtering, however this does not give the expected improvements
!  v13: Introduce a dynamic search window based on the covariance matrix
!        - Window is implemented as a parabolic multiplier normalized to unity and is zero
!           at a certain multiple of the search window. This factor is named 'SearchRangeFallOff',
!           implemented in 'ReImAtMax'.
!  v13: Initial choise for source location is based on a double grid search, see 'SourceTryal'.
!        This routine also produces a first guess for the diagonal covariance matrix. It requires an
!        maximal window set by 'Safety' to be at least 50.
!  v13: Introduced an estimator for the shape of the cross correlation kept in array 'CC_Qual', as
!        calculated in subroutine 'GetCorrSingAnt'. The shape is calculated in 'ReImAtMax' as the width
!        of the cross correlation function. The quality is calculated as the ratio with the width of the
!        reference antenna. If 'CC_Qual' is worse than 'CCShapeCut' (unity is perfect) the antenna is
!        cut in 'GetCorrSingAnt'.
!  v16: include noise level cut-off
!  v18: Peak location optimized, in particular confirmed that 'i_distmin=2' works well
!        and that 'FitRange_Samples=??' is optimal
!  v18: Peak-location statistics file written
!  v18: 'PeaksPerChunk' introduced as a parameter
   use constants, only : dp,sample,c_mps
   use DataConstants, only : Production, Time_dim, EdgeOffset
   use Chunk_AntInfo, only : Ant_Stations, StartT_sam, NoiseLevel
   use Chunk_AntInfo, only : MaxPeaksPerChunk, PeaksPerChunk!, PeakSP, PeakSWl, PeakSWu, PeakSAmp, PeakD_nr
   use Chunk_AntInfo, only : Station_OutOfRange, Station_nrMax, Nr_UniqueStat
   use ThisSource, only : Nr_Corr, PeakNr, PeakNrTotal, TotPeakNr, PeakPos, Peak_eo, ChunkNr
   use ThisSource, only : RefAntErr, SourcePos, Dual, Peak_Offst, CCShapeCut
   use ThisSource, only : CCShapeCut_lim, ChiSq_lim, EffAntNr_lim
   use FitParams, only : MaxFitAntDistcs, FitIncremental, Fit_PeakNrTotal, FitQual, Sigma, MaxFitAntD_nr
   use FitParams, only : N_FitStatTim, Nr_TimeOffset, PulsPosCore, CalcHessian, N_EffAnt, Max_EffAnt  ! N_FitPar_max,
   use FitParams, only : FullSourceSearch, SigmaGuess, SpaceCov
   use StationMnemonics, only : Statn_ID2Mnem, Station_Mnem2ID
   use Explore_Pars, only : NMin, NMax, Emin, EMax
   Implicit none
   Real(dp), intent(in) :: SourceGuess(3)
   Integer, intent(in) :: TimeFrame  !  used to generate unique pulse-number
   Integer, intent(in) :: units(0:2)   ! used for writing fitresults
   !Integer, parameter :: Separation=Safety ! determines part of the ref spectrum that is zeroed after finding a peak; Safety from 'ThisSource'
   !real ( kind = 8 ) :: X(N_FitPar_max)
   Real(dp), parameter :: Dist_Lim=2.5d5
   !Integer, parameter :: PeakS_dim=350
   !
   integer :: i, j, k, i_ant, i_eo, i_chunk=1, i_eo_s=0, i_eo_f=1
   integer :: i_Peak, StLoc, StatMax, nxx, i_peakS, i_dist, i_distmin
   Integer :: OutOfRange(1:Station_nrMax)
   Real(dp) :: DistMax, KalSource(0:3), KalCoVariance(0:3,0:3)
   Integer :: Date_T(8)
   Integer :: PeakSP(MaxPeaksPerChunk,0:2), PeakSWl(MaxPeaksPerChunk,0:2), PeakSWu(MaxPeaksPerChunk,0:2)
   Integer :: PeakSAmp(MaxPeaksPerChunk,0:2), PeakD_nr
   Integer :: Wu,Wl, CalSource(10),i_cal, PeakNrSearched, PeakNrFound, PeakNrGood, i_FvStr
   Integer :: Peakpos_0 ! pulse position for an virtual antenna at the core of CS002
   integer :: NoSourceFound
   Real(dp), external :: RefracIndex
   !integer, external :: XIndx
   !       Initialize
   PulsPosCore=.false.  ! Pulse positions are on the ref. antenna
   !
   Call Find_unique_StatAnt()
   Call GetRefAnt(i_chunk)
   !
   Call DualPeakFind(PeaksPerChunk, i_chunk, PeakD_nr, PeakSP, PeakSWl, PeakSWu, PeakSAmp)
   !
   N_FitStatTim= 0   ! Do not fit any station/antenna timing
   Nr_TimeOffset=1
   Fit_PeakNrTotal=1     ! Fit only a single source at a time
   i_cal=0
   !
   ChunkNr(1) =1
   ChunkNr(2) =1
   Peak_eo(2)=1
   If(Dual) then
      i_eo_s=2
      i_eo_f=2
   Endif
   Do i_eo=i_eo_s,i_eo_f  !  0,1
!      Dual=.false.
      PeakNr(0:1,1)=1
      TotPeakNr(0:1,1)=1
      PeakNrTotal=1
      If(i_eo .eq. 0) then
         PeakNr(1,1)=0
         Peak_eo(1)=i_eo
      ElseIf(i_eo .eq. 1) then
         PeakNr(0,1)=0
         TotPeakNr(0,1)=0
         Peak_eo(1)=i_eo
      Else     ! i_eo = 2
!         Dual=.true.
         TotPeakNr(1,1)=2
         Peak_eo(1)=0
         PeakNrTotal=2
      Endif
      !Do i=1,4
      !   write(2,*) 'i=',i,XIndx(i,1)
      !Enddo
      !write(2,*) '   Fit_PeakNrTotal=1',    Fit_PeakNrTotal
      !
      PeakNrSearched=0  ; PeakNrFound=0  ; PeakNrGood=0
      Do i_peakS=1,PeaksPerChunk ! ,10
         !FitNoSources=.false.
         If(PeakSAmp(i_peakS,i_eo).lt.NoiseLevel) cycle
         SourcePos(:,1)=SourceGuess(:)
         Peakpos(1)=PeakSP(i_peakS,i_eo)
         If(peakPos(1).le.0) exit
         SourcePos(:,2)=SourceGuess(:)
         Peakpos(2)=PeakSP(i_peakS,i_eo)
         RefAntErr(1)=0.
         RefAntErr(2)=0.
         If(.not.production) write(*,"(A,i2,i4)", ADVANCE='NO') &
            achar(27)//'[92m'//achar(27)//'[100D'//'peak#=', i_eo, i_peakS !, achar(27)//'[0m.'
         !Write(*,"(A,I4)") 'Peak#',i_peakS
         !write(2,"(A,i2,I4,A,I6,3F11.2,8('-'))") 'Par, Peak#', i_eo,i_peakS,', @:',Peakpos(1), SourcePos(:,1)
         !If(Dual) write(2,"(20x,I6,3F11.2,8('-'))") Peakpos(2), SourcePos(:,2)
         !StatMax=Ant_Stations(1,1)
         StatMax=1050
         CalcHessian=.false.
         If(.not.production) write(2,"(A,i3,A,i6,A,i6,A)") 'FindSource; i_peak:',i_peakS,'@',Peakpos(1),&
            ', ampl=',PeakSAmp(i_peakS,i_eo),'=================='
         !production=.false.
         !
         CCShapeCut=0.1  ! do not cut on shape of Cross Correlation function
         If(FullSourceSearch) then
            i_distmin=2
            DistMax=MaxFitAntDistcs(i_distmin) ! Distance will effectively be limited in 'SourTryal' by the value of 'FitRange'
      !      Call SourceTryal(DistMax,i_peakS, NoSourceFound)
      !      Call SourceTryal(DistMax,i_peakS, NoSourceFound)
      !      Call SourceTryal(DistMax,i_peakS, NoSourceFound)
            Call SourceTryal_v2(DistMax,i_peakS, NoSourceFound)
            !write(2,*) 'NoSourceFound:',NoSourceFound
            If(NoSourceFound.ne.0) Then
               If(.not.production) write(2,*) '================',i_peakS,' LMA search failed ==============='
               cycle
            EndIf
            !Call SourceTryal(MaxFitAntDistcs(i_distmin))
            i_distmin=2
         Else
            i_distmin=2
            Do i_Peak=1,PeakNrTotal
               SourcePos(:,i_Peak)=SourceGuess
            Enddo
            Do i=1,3
               SpaceCov(i,i)=SigmaGuess(i)*SigmaGuess(i)
            Enddo
         EndIf
         !Write(2,*) 'i_peakS:',i_peakS
         !
         Do i_dist=i_distmin,MaxFitAntD_nr
            DistMax=MaxFitAntDistcs(i_dist) ! in kilometers
            If(DistMax.eq.-1) exit
            !If(i_dist.eq.MaxFitAntD_nr) CalcHessian=.true.
            CalcHessian=.true.
            !If(i_dist .gt. i_distmin) production=.true.
            Call SourceFitCycle(StatMax,DistMax)
            ! Real(dp), save :: CCShapeCut_lim=0.8        ! limit set on CCShapeCut
            CCShapeCut=CCShapeCut_lim ! Should be lower for Dual option
            !stop 'FindStatCall-1'
            !if(i_peakS.eq.2) stop 'FindStatCall-1'
            !If(DistMax.gt.25.) write(2,*) '(sqrt(SpaceCov(i,i)),i=1,3):',DistMax,(sqrt(SpaceCov(i,i)),i=1,3),Sigma
         EndDo ! i_dist=1,MaxFitAntD_nr
         !When the Hesian was calculated the covariance matrix is stored, as well as Sigma, otherwise sigma=large
         ! as this throws out many good sources, we repair this here (April 2022)
         Do i=1,3
            sigma(i)=sqrt(SpaceCov(i,i))
         Enddo
         !
         If(units(i_eo).lt.0) then
            Write(2,"(I8,', ',3(f12.4,', '),f9.3)") & ! ID (x,y,z) t chi^2 Diag(covariance)
               (TimeFrame*1000+i_peakS),SourcePos(:,1), FitQual
            cycle
         Endif
   !   STOP
         !
         !SourceSPos(:,i_eo,i_peakS)=SourcePos(:,1)
         !Peakpos_0=Peakpos(1) + NINT(Peak_Offst(1)) ! latter is calculated in GetCorrSingAnt
         Peakpos_0=Peakpos(1) - NINT(Peak_Offst(1)) ! latter is calculated in GetCorrSingAnt; sign changed to - March 7, 2022
         DistMax=sqrt(sum(SourcePos(:,1)*SourcePos(:,1)))
         PeakNrSearched= PeakNrSearched +1
         !N_EffAnt, Max_EffAnt
         ! ChiSq_lim Limit should be higher for Dual option;
         ! Real(dp), save :: ChiSq_lim=10              ! limit set on the chi-square (FitQual)
         ! Real(dp), save :: EffAntNr_lim=0.8        ! limit set on the ratio of effective number of used antennas v.s. available number
         If(Max_EffAnt.lt.10) exit
         If(i_peakS.eq.1 .and. i_eo.ne.1) Then
            i_FvStr=0
            Write(18,"(' ', F13.6, 3F11.2,i9)")  &
                        1000.d0*StartT_sam(1)*sample, SourcePos(:,1)/1000., (TimeFrame*1000+i_peakS)
            Flush(Unit=18)
         EndIf
         If(.not.production) write(2,*) 'selectioncriterium', DistMax, Dist_Lim, (N_EffAnt*1./Max_EffAnt), EffAntNr_lim, &
            Sigma(1), Sigma(2), Sigma(3) , 10.*max(99.,1000./(abs(SourcePos(3,1))+0.001)) , Sigma(3), &
            FitQual, ChiSq_lim
         If((DistMax .lt. Dist_Lim) .and. ( (N_EffAnt*1./Max_EffAnt) .gt. EffAntNr_lim) &
            .and. ( Sigma(1) .lt. 999. ) .and. ( Sigma(2) .lt. 999.) &
            .and. ( Sigma(3) .lt. 10.*max(99.,1000./(abs(SourcePos(3,1))+0.001)) ) .and. ( Sigma(3) .gt. 0.) &
            .and. (FitQual .lt. ChiSq_lim)) then
               If( Sigma(3) .gt. 999.99) Sigma(3)=999.990
               Wl=PeakSWl(i_peakS,i_eo)
               Wu=PeakSWu(i_peakS,i_eo)
               !write(2,*) 'NINT(Peak_Offst(1)):',NINT(Peak_Offst(1))
               Write(units(i_eo),"(I8,', ',3(f12.3,', '),f12.9,', ',f9.3,',',3(f8.4,','),2(I4,','),I7,',',2(I3,','))") & ! ID (N,E,z) t chi^2 Diag(covariance)
                  (TimeFrame*1000+i_peakS),SourcePos(:,1), &
                  (StartT_sam(1)+Peakpos_0)*sample-DistMax*RefracIndex(SourcePos(3,1))/c_mps, FitQual &
                  , Sigma(1:3), N_EffAnt, Max_EffAnt, PeakSAmp(i_peakS,i_eo), Wl, Wu
               !write(*,*) 'Wu,Wl',Peakpos_0,Wl,Wu
               Call CleanPeak(i_eo,Peakpos_0,SourcePos(1,1),Wl,Wu, OutOfRange) ! clean peak only when decent source location was obtained
               PeakNrFound=PeakNrFound+1
               Station_OutOfRange(1:Nr_UniqueStat)=Station_OutOfRange(1:Nr_UniqueStat)+ OutOfRange(1:Nr_UniqueStat)
               If(FitQual.lt.25.) then
                  NMin = Min(NMin,SourcePos(1,1)/1000.) ! Extremes in [km]
                  EMin = Min(EMin,SourcePos(2,1)/1000.)
                  NMax = Max(NMax,SourcePos(1,1)/1000.)
                  EMax = Max(EMax,SourcePos(2,1)/1000.)
                  PeakNrGood=PeakNrGood+1
               Endif
               !
               !searching for good peaks to be used for improving calibration
               !If(i_peakS.lt.10 .and. wl.le.10 .and. Wu.le.10) &
               If(i_peakS.lt.10 .and. (FitQual .lt. ChiSq_lim/2.) ) Then
                  i_FvStr=i_FvStr+1
                  Write(18,"('C',i2,' 2 1',I8,3(F10.2,','),F12.5,';',f9.3,',',2(I4,','),I7,',',2(I3,','),2I3)") &
                     i_FvStr, Peakpos_0, SourcePos(:,1), &
                     (StartT_sam(1)+Peakpos_0)*sample-DistMax*RefracIndex(SourcePos(3,1))/c_mps , FitQual &
                     ,  N_EffAnt, Max_EffAnt, PeakSAmp(i_peakS,i_eo), Wl, Wu, i_eo, i_peakS
               EndIf
               If(wl.le.5 .and. Wu.le.5 .and. FitQual.lt. 20. .and. i_cal .lt.10) Then ! see for a possible 5-star source
                  If(PeakSAmp(i_peakS,i_eo).gt.PeakSAmp(1,i_eo)/10.) then
               !StartT_sam(1)=(StartTime_ms/1000.)/sample + (TimeFrame-1)*(Time_dim-2*EdgeOffset)  ! in sample's
                  If(i_eo.eq.0) Then
                     i_cal=i_cal+1
                     if(i_cal.le.10) CalSource(i_cal)=Peakpos_0
                  Else
                     If(i_cal.gt.10) i_cal=10
                     Do i=1,i_cal
                        If(abs(Peakpos_0-CalSource(i_cal)).lt.6) Then
                           write(2,*) 'got a good one!!!!!!!!!!!!!!!'
                           write(18,"(A, F13.6, F12.6, I7, 3F9.3, f5.1, f6.1, I7, 3I3)") 'Start time[ms]:', &
                              1000.d0*StartT_sam(1)*sample, (TimeFrame-1)*(Time_dim-2*EdgeOffset)*1000.d0*sample, Peakpos_0, &
                              SourcePos(:,1)/1000.d0, FitQual, 100.*N_EffAnt/Max_EffAnt, PeakSAmp(i_peakS,i_eo), Wl, Wu,i_eo
                           Write(18,"(A,I8,':', F13.5, 3F11.2,i3)") '***** label=', (TimeFrame*1000+i_peakS), &
                              1000.d0*StartT_sam(1)*sample, SourcePos(:,1), i_cal
                           Write(18,"(A,I8,3(F10.2,','),F12.5,'; 0.0')") 'R 1 2 1',  (Peakpos_0+CalSource(i_cal))/2, &
                              SourcePos(:,1), (StartT_sam(1)+Peakpos_0)*sample-DistMax*RefracIndex(SourcePos(3,1))/c_mps
                           Flush(unit=18)
                           exit
                        EndIf
                     Enddo
                  EndIf
               EndIf
            EndIf
            !
            If(.not.production) write(2,*) '================',i_peakS,' is located ===============', Sigma(1:3), &
               N_EffAnt, Max_EffAnt, PeakSAmp(i_peakS,i_eo), Wl, Wu
         Else
            If(.not.production) write(2,*) '================',i_peakS,' out of bounds ===============', Sigma(1:3), &
               N_EffAnt, Max_EffAnt, PeakSAmp(i_peakS,i_eo), Wl, Wu
         Endif
         !If(i_peakS.gt.6) stop 'SourceFind'
         !If(i_peakS.gt.20) stop 'SourceFind'
      EndDo !  i_peakS=1,PeaksPerChunk
      write(units(i_eo)+10,"(I8,',',F8.2,',',I4,',',I4,',',I4)") TimeFrame, StartT_sam(i_chunk)*1000.d0*sample, &
         PeakNrSearched, PeakNrFound, PeakNrGood
      write(2,*) 'located pulses', PeakNrSearched, PeakNrFound, PeakNrGood
      write(*,*) 'located pulses', PeakNrSearched, PeakNrFound, PeakNrGood
      !
      CALL DATE_AND_TIME (Values=DATE_T)
      WRITE(*,"(1X,I2,':',I2.2,':',I2.2,'.',I3.3,A)") (DATE_T(i),i=5,8), achar(27)//'[0m'
      !write(*,*) ' ', achar(27)//'[0m'
   Enddo ! i_eo=0,1
   flush(unit=18)
    !
    !
    ! Merge values for "Fit_TimeOffsetStat" with input from 'FineCalibrations.dat'
!1   continue
    !
    !Call CleanPeak(i_eo)
    !
    !
    !Call AntFieParGen()
    !Call AllStatStokes()
    !
    Return
End Subroutine SourceFind
!=================================
Subroutine SourceFitCycle(StatMax,DistMax)
   use DataConstants, only : Production
    use ThisSource, only : Nr_Corr, RefAntErr, SourcePos, PeakNrTotal !, Peakpos, PlotCCPhase
    use FitParams, only : N_FitPar, N_FitPar_max, N_FitStatTim, FitParam, SpaceCov, FitQual, N_EffAnt, Max_EffAnt
    use constants, only : dp
    Implicit none
    !logical, intent(inout) :: FitFirst
    !logical, intent(in) :: FitNoSources
    integer, intent(in) :: StatMax
    Real(dp), intent(in) :: DistMax
    !Real(dp), intent(out) :: FitQual
    Real ( kind = 8 ) :: X(N_FitPar_max)
    integer :: FitPos(4), i, i_Peak ! , j, k, i_ant, j_corr, i_eo, i_chunk, Station
    integer, external :: XIndx
    !Real(dp) :: SourcePos(3,1)
    !logical :: FitNoSources=.false.
    !
    !
    If(.not. Production) write(2,*) '=== New Round ============================================================'
    !write(2,*) 'Call BuildCC(StatMax,DistMax):',StatMax,DistMax
    Call BuildCC(StatMax,DistMax)
    !
    FitPos(1)=1 ; FitPos(2)=2 ; FitPos(3)=3 ; FitPos(4)=4
    !CalcHessian=.false.
    If(DistMax.le. 0.05) then
        FitPos(4)=-1 ; FitPos(1)=1 ; FitPos(2)=-1 ; FitPos(3)=-1 ! single station, fit direction only
    ElseIf(DistMax.le. 0.3) then
        FitPos(4)=-1 ; FitPos(1)=1 ; FitPos(2)=2 ; FitPos(3)=-1 ! Superterp, fit direction and curvature front
    ElseIf(DistMax.le. 1.6) then
        FitPos(4)=-1 ; FitPos(1)=1 ; FitPos(2)=2 ; FitPos(3)=3
    Endif
    !
   N_FitPar=0
   Do i=1,4        ! aug19: changed from 3 to 4 where 4 is the error in the reference antenna due to noise for this peak
      if((FitPos(i) .le. 0) .or. (FitPos(i) .gt. 4) ) exit
      N_FitPar = N_FitPar +1
      FitParam(N_FitPar)=FitPos(i)
      Do i_Peak=1,PeakNrTotal
         If(FitPos(i).eq.4) then
           X(  XIndx(i,i_Peak) ) = RefAntErr(i_Peak)
         elseIf(FitPos(i).eq.3) then
           X(  XIndx(i,i_Peak) ) = abs(SourcePos(FitPos(i),i_Peak))
         Else
           X(  XIndx(i,i_Peak) ) = SourcePos(FitPos(i),i_Peak)
         EndIf
      Enddo
   enddo
   If(.not. Production) write(2,*) 'X',X( 1:N_FitPar )
    !write(2,*) 'N_FitPar1 , N_FitStatTim',N_FitPar , N_FitStatTim
    !
    !stop 'FindStatCall'
   If(.not. Production)  &
      write(2,*) 'New fit starting, Nr_Corr=',Nr_Corr,'============================================================'
   If(.not. Production)  then
      If(StatMax.lt.50) then
         write(2,"(A,I4,A,i4,A,F7.2,A)") 'N_FitPar=',N_FitPar , &
            ', max station nr in fit=',StatMax, ', max distance to ref. station=',DistMax,'[km]'
      else
         write(2,"(A,I4,A,F7.2,A)") 'N_FitPar=',N_FitPar,', max distance to ref. station=',DistMax,'[km]'
      endif
   endif
   Call FitCCorr(X)  ! fit source
   Call X2Source(X)
   !If(N_FitPar.eq.4) Then
   !   If(Production .and. DistMax.gt.20.)  write(2,"(A,F7.2,A,4F9.1, A,F8.2, A,2I4, A,3G11.3)") &
   !      'Max distance=',DistMax,'[km], fitted:', X( 1:N_FitPar ),', chi^2/ndf=', FitQual &
   !      ,' N_Ant=',N_EffAnt,Max_EffAnt, ', sigma[m]:',(sqrt(SpaceCov(i,i)),i=1,3)
   !Else
   !   If(Production .and. DistMax.gt.20.)  write(2,"(A,F7.2,A,3F9.1, A,F8.2, A,2I4, A,3G11.3)") &
   !      'Max distance=',DistMax,'[km], fitted:', X( 1:3 ),', chi^2/ndf=', FitQual &
   !      ,' N_Ant=',N_EffAnt,Max_EffAnt, ', sigma[m]:',(sqrt(SpaceCov(i,i)),i=1,3)
   !EndIf
   !
End Subroutine SourceFitCycle
!=====================================
Subroutine CleanPeak(i_eo_in,PeakPos,SourcePos,Wl,Wu, OutOfRange)
!  Zero spectrum at the good peak inorder not to find it twice
!  v13:  Zeroing window is made peak dependent
   use ThisSource, only : Nr_Corr, CorrAntNrs, Tref_dim, RefAntErr, Peak_Offst
   use Chunk_AntInfo, only : CTime_spectr, Ant_pos, Ant_RawSourceDist, Ant_Stations, Station_nrMax
   use DataConstants, only : Time_dim
   use constants, only : dp,pi
   Implicit none
   integer, intent(in) :: i_eo_in, PeakPos, Wl, Wu
   Real(dp), intent(in) :: SourcePos(3)
   integer, intent(out) :: OutOfRange(1:Station_nrMax)
   Integer :: j_corr, i_chunk=1, i_ant, i_Peak=1, StLoc, i, HW_size=5, Zero_dim, i_eo, i_eo_s,i_eo_f, Wi, i_stat
   Real(dp) :: Rdist,Hann(1:2*Tref_dim)

   !write(*,*) i_eo,PeakPos,SourcePos,Tref_dim
   Hann=1.
   Zero_dim=Wl+Wu+2
   If(Zero_dim .gt. 2*Tref_dim) then
      Zero_dim=2*Tref_dim
      Wi=Tref_dim
   else
      Wi=Wl+1
   endif
   !write(*,*) 'PeakPos',PeakPos
   !write(*,*) 'Zero_dim=',Zero_dim
   Do i=0,HW_size
        Hann(i+1) = sin(0.5*i*pi/HW_size)**2
        Hann(Zero_dim-i)=Hann(i+1)
    EndDo

   !write(*,*) Hann
   If(i_eo_in.eq.2) then
      i_eo_s=0 ; i_eo_f=1
   Else
      i_eo_s=i_eo_in ; i_eo_f=i_eo_in
   Endif
   OutOfRange(:)=0
   Do i_eo=i_eo_s,i_eo_f
      Do j_corr=1,Nr_Corr(i_eo,i_chunk) ! assumes PulsPosCore=.false.
         i_ant=CorrAntNrs(j_corr,i_eo,i_chunk)
         Call RelDist(SourcePos,Ant_pos(1,i_ant,i_chunk),RDist)
         Rdist=Rdist - Ant_RawSourceDist(i_ant,i_chunk) +  RefAntErr(i_Peak)! - INT(Peak_Offst(i_Peak))
         StLoc=PeakPos + INT(Rdist) - Peak_Offst(i_Peak) - Wi
         !write(*,*) j_corr,StLoc
         If((StLoc.lt.1).or.(StLoc.gt.Time_dim-Zero_dim)) then
            write(2,*) 'ERROR, CleanPeak out of bounds',StLoc,Time_dim,', for:', &
                  j_corr, Ant_Stations(i_Ant,i_chunk),RDist,'[samples]'
            Call Conv_Station_ID2i(Ant_Stations(i_Ant,i_chunk), i_stat)
            OutOfRange(i_stat)=1
            cycle
         Endif
         !If(PeakPos.eq.10657) write(*,*) i_eo,J_corr,i_ant,StLoc
         Do i=1,Zero_dim
            !write(*,*) 'i',i
            CTime_spectr(StLoc+i,i_ant,i_chunk) = CTime_spectr(StLoc+i,i_ant,i_chunk)*(1.-Hann(i))
         Enddo
      Enddo
   Enddo
   Return
End Subroutine CleanPeak
!=====================================
Subroutine DualPeakFind(PeakS_dim, i_chunk, PeakD_nr, PeakSP, PeakSWl, PeakSWu, PeakSAmp)
! To optimize:
!  - reduce "Safety" to 10 ??? v Taken care of in v14
!  - "Time_dim=20" seems pretty realistic for the total width of the pulse, maybe reduce to 15?
!  - Modify "Separation" to Safety+Time_dim  v Taken care of in v14
!  v13: Excluded region around peak adjusted
!
   use Chunk_AntInfo, only : PeaksPerChunk !, PeakSP, PeakSWl, PeakSWu, PeakSAmp, PeakD_nr
   !use Chunk_AntInfo, only : PeakSPw, SPeak
   use Chunk_AntInfo, only : CTime_spectr, Ant_IDs, Ant_Pos, Ant_nrMax, RefAnt, TimeFrame !, NoiseLevel
   use ThisSource, only : Safety, Dual, Tref_dim
   use unque, only : Double_sort
   use constants, only : dp
   use DataConstants, only : EdgeOffset, Time_dim, RunMode  ! , ChunkNr_dim
   use DataConstants, only : Production
   Implicit none
   Integer, intent(in) :: PeakS_dim, i_chunk
   Integer, intent(out) :: PeakSWu(PeakS_dim,0:2), PeakSWl(PeakS_dim,0:2), PeakSP(PeakS_dim,0:2)  ! Peaks positions
   Integer, intent(out) :: PeakSAmp(PeakS_dim,0:2)  ! Peaks amplitude in ref station
   Integer, intent(out) :: PeakD_nr    ! Nr of peaks found in both even & odd antennas
   Integer :: Separation ! =Safety ! determines part of the ref spectrum that is zeroed after finding a peak; Safety from 'ThisSource'
   real(dp), parameter :: safe=1.d-8 ! Almost =0, used for bookkeeping of signals found
   !Integer, parameter :: W_low= 6 ! half the smallest distance between peaks
   Integer, parameter :: W_low= 3 ! Min value for Wl & Wu; (Wl,Wu) + Tref_dim/2 = half the smallest distance between peaks
   !
   Integer :: PeakSPw(PeakS_dim,0:1) , SPeak(PeakS_dim,1:4,0:1)

   integer :: i, j, k, i_ant, i_eo
   integer :: i_loc(1), i_Peak
   real(dp) :: HEnvel(Time_dim), PeakTail
   Integer :: PSP,Wu,Wl,w, peak0,peak1,Peak2, pos0, pos1, WindW
   real(dp) :: AvePos,TotAmpl, StDev
   !
   Separation=Safety
   !
   PeakSP(:,:)=0
   Do i_eo=0,1
      i_ant= RefAnt(i_chunk,i_eo)
      If(.not.production) &
         write(2,*) i_eo, ', Reference antenna=',i_ant,Ant_IDs(i_ant,i_chunk),', @position ',Ant_pos(:,i_ant,i_chunk)
       HEnvel(:)=abs(CTime_spectr(:,i_ant, i_chunk))
      ! write(label,"(I2.2,I3.3)") TimeFrame,Ant_IDs(i_ant,i_chunk)
      ! OPEN(unit=30,FILE='TimeTrace'//label//'.dat',FORM='FORMATTED',STATUS='unknown')
      ! Do j=1,Time_dim
      !    write(30,*) j,Real(CTime_spectr(j,i_ant, i_chunk)),Imag(CTime_spectr(j,i_ant, i_chunk)),HEnvel(j)
      ! Enddo
      ! Close(Unit=30)
       i_peak=1
       Do j=1,2*PeaksPerChunk
         i_loc=MaxLoc( HEnvel(EdgeOffset:Time_dim-EdgeOffset) )
         PSP=i_loc(1) + EdgeOffset - 1
         !write(2,*) 'j=',j,psp,HEnvel(PSP)
         PeakTail=HEnvel(PSP)/2.  ! pulse should be reduced to less than PeakTail to determine width
         PeakSAmp(i_peak,i_eo)=NINT(HEnvel(PSP))
         !Write(2,*) 'PSP',j,i_peak,PSP
         !Write(2,"(11F12.1)")  HEnvel(PSP- 10:PSP + 10)
         Wu=20  ! determine the position of the min in the spectrum above the peak
         Do i=W_low+1,20
            If(HEnvel(PSP+i) .gt. PeakTail ) cycle ! only proceed when less that (/2) peak
            If(HEnvel(PSP+i) .le.0. ) then   ! reached close to previous pulse
               Wu=W_low  ! i-1
               Exit
            Elseif(HEnvel(PSP+i) .lt. (HEnvel(PSP+i+1)+safe) ) then ! check if local min is reached.
               Wu=i
               Exit
            Endif
         Enddo
         Wl=20 ! determine the position of the min in the spectrum below the peak
         Do i=W_low+1,20
            If(HEnvel(PSP-i) .gt. PeakTail ) cycle
            If(HEnvel(PSP-i) .le.0. ) then
               Wl=W_low  ! i-1
               Exit
            Elseif(HEnvel(PSP-i) .lt. (HEnvel(PSP-i-1)+safe) ) then
               Wl=i
               Exit
            Endif
         Enddo
         If((Wu.eq.W_low) .or. (Wl.eq.W_low)) then  ! Peak too close to a previously located one
            HEnvel(PSP - Wl:PSP + Wu)=safe
         Else
            !PeakSW(i_peak)=Wu+Wl+1
            PeakSP(i_peak,i_eo)=PSP
            PeakSWu(i_peak,i_eo)=Wu
            PeakSWl(i_peak,i_eo)=Wl
            PeakSPw(i_peak,i_eo)=PSP + Wu-Wl
            !HEnvel(PSP - Separation:PSP + Separation)=-safe  ! zero the test spectrum around the peak just found
            !HEnvel(PSP - Wl - Tref_dim/2:PSP + Wu + Tref_dim/2)=-safe  ! zero the test spectrum around the peak just found
            HEnvel(PSP - Wl-W_low:PSP + Wu +W_low )=-safe  ! zero the test spectrum around the peak just found over an interval that matches that of "CleanPeak"
            i_peak=i_peak+1
         Endif
         if(i_peak.eq.(PeakS_dim+1)) exit  ! PeaksPerChunk+1)) exit
      Enddo   !  j=1,2*PeaksPerChunk
      i_peak=i_peak-1
      HEnvel(:)=abs(CTime_spectr(:,i_ant, i_chunk))  ! to get rid of all zeroed parts
      If(.not.production) Write(2,"(A,i3,A)", ADVANCE='NO') 'Peak(',i_peak,'):'
      If(.not.production) Write(2,"(20i6)") PeakSP(1:i_peak,i_eo)
      If(.not.production) Write(2,"(A)", ADVANCE='NO') 'Peakval = '
      If(.not.production) Write(2,"(20i6)") PeakSAmp(1:i_peak,i_eo)
      !Write(2,"(A)", ADVANCE='NO') 'PeakSWl='
      !Write(2,"(20i6)") PeakSWl(:,i_eo)
      !Write(2,"(20f6.0)") HEnvel(PeakSP(:,i_eo)-PeakSWl(:,i_eo))
      !Write(2,"(A)", ADVANCE='NO') 'PeakSWu='
      !Write(2,"(20i6)") PeakSWu(:,i_eo)
      !Write(2,"(20f6.0)") HEnvel(PeakSP(:,i_eo)+PeakSWu(:,i_eo))
      !Write(2,"(A)", ADVANCE='NO') 'PeakSd='
      !Write(2,"(20i6)") PeakSP(i_eo,:)-PeakSPw(i_eo,:)
      SPeak(:,1,i_eo)=0
      SPeak(1:i_peak,1,i_eo)=PeakSP(1:i_peak,i_eo)
      SPeak(1:i_peak,2,i_eo)=PeakSWu(1:i_peak,i_eo)
      SPeak(1:i_peak,3,i_eo)=PeakSWl(1:i_peak,i_eo)
      !SPeak(:,4,i_eo)=NINT(HEnvel(PeakSP(1:i_peak,i_eo)))
      SPeak(1:i_peak,4,i_eo)=PeakSAmp(1:i_peak,i_eo)
      Call Double_sort(SPeak(1:i_peak,:,i_eo))  ! re-order on peak position; only for internal use
      If(RunMode.eq. 2) Then   ! test for calibration quality
         Windw=100
         i=0
         j=0
         Do j=1,i_peak
            If(.not.production) write(2,"(I2)", Advance='no') j
            Call Inv5StarPk(HEnvel, PeakSP(j,i_eo), Windw, AvePos, TotAmpl, StDev)
            !write(2,*) '***** ',j,PeakSP(j,i_eo),'position offset=',AvePos,', StDev=',StDev, &
            !', power=',StDev*PeakSAmp(j,i_eo), ' or', TotAmpl
            If(StDev.lt.25) Then
               i=i+1
               PSP=PeakSP(i,i_eo)
               Wu=PeakSWu(i,i_eo)
               Wl=PeakSWl(i,i_eo)
               w=PeakSAmp(i,i_eo)
               !
               PeakSP(i,i_eo)=PeakSP(j,i_eo)
               PeakSWu(i,i_eo)=PeakSWu(j,i_eo)
               PeakSWl(i,i_eo)=PeakSWl(j,i_eo)
               PeakSAmp(i,i_eo)=PeakSAmp(j,i_eo)
               !
               PeakSP(j,i_eo)=PSP
               PeakSWu(j,i_eo)=Wu
               PeakSWl(j,i_eo)=Wl
               PeakSAmp(j,i_eo)=w
            EndIf
         EndDo
         If(.not.production) write(2,*)
      EndIf
      !write(*,*) 'done:',i_eo
   enddo  ! i_eo=0,1
   !            PeakSWl(i_eo,i_peak)=Wl
   If(.not. Dual) then
      PeakD_nr=0
      Return
   endif
   !
   !stop
   !Flush(unit=2)
   !write(*,*) 'end polarity'
   If(.not.production) write(2,"(A)", ADVANCE='NO') 'even sort:'
   If(.not.production) Write(2,"(20i6)") SPeak(:,1,0)
   !SPeak(:,1,1)=PeakSP(:,1)
   !SPeak(:,2,1)=PeakSWu(:,1)
   !SPeak(:,3,1)=PeakSWl(:,1)
   !Call Double_sort(SPeak(:,:,1))
   If(.not.production) write(2,"(A)", ADVANCE='NO') ' odd sort:'
   If(.not.production) Write(2,"(20i6)") SPeak(:,1,1)
   peak0=1
   peak1=1
   peak2=0
   !PeakDP(:)=0
   !write(2,*) 'PeaksPerChunk',PeaksPerChunk, PeakS_dim
   Do While(max(peak0,peak1,peak2).le.PeakS_dim)
   !write(*,*) peak0,peak1,Peak2
   !write(2,*) peak0,peak1,Peak2
   !Flush(unit=2)
      pos0=SPeak(peak0,1,0)
      if(pos0.lt.EdgeOffset) exit
      pos1=SPeak(peak1,1,1)
      if(pos1.lt.EdgeOffset) exit
      If(pos0.lt.pos1) then ! peak 0 earlier than peak 1
         w=MAX(SPeak(peak0,2,0),SPeak(peak1,3,1))  ! require that the tail of one overlaps with the other
         !w=min(SPeak(peak0,2,0),SPeak(peak1,3,1)) ! require that upper tail of the earlier overlaps with the later peak and v.v.
         If((pos1-pos0) .lt. w) then
            Peak2=Peak2+1
            SPeak(Peak2,1,0)=-(SPeak(Peak0,4,0)+ SPeak(Peak1,4,1)) ! Peak strength; first and negative for sorting purposes
            SPeak(Peak2,4,0)=(pos0+pos1)/2 ! position is the average of the early and late one
            SPeak(Peak2,2,0)=MAX( (pos0+SPeak(Peak0,2,0)),(pos1+SPeak(Peak1,2,1)) ) - (pos0+pos1)/2   ! Wu
            SPeak(Peak2,3,0)=(pos0+pos1)/2 - MIN( (pos0-SPeak(Peak0,3,0)),(pos1-SPeak(Peak1,3,1)) )   ! Wl
            Peak0=Peak0+1
            Peak1=Peak1+1
         Else ! Peaks are too far apart to count as common peak
            Peak0=Peak0+1  ! advance the earlier one
         Endif
      else  ! peak 1 earlier than peak 0 (pos1 < pos0)
         w=MAX(SPeak(peak1,2,1),SPeak(peak0,3,0))
         !w=min(SPeak(peak1,2,1),SPeak(peak0,3,0))
         If((pos0-pos1) .lt. w) then  ! peak in both polarizations
            Peak2=Peak2+1
            SPeak(Peak2,1,0)=-(SPeak(Peak0,4,0)+ SPeak(Peak1,4,1)) ! Peak strength, Peak2 is necessarily smaller than Peak0 or Peak1
            SPeak(Peak2,4,0)=(pos0+pos1)/2
            SPeak(Peak2,2,0)=MAX( (pos0+SPeak(Peak0,2,0)),(pos1+SPeak(Peak1,2,1)) ) - (pos0+pos1)/2   ! Wu
            SPeak(Peak2,3,0)=(pos0+pos1)/2 - MIN( (pos0-SPeak(Peak0,3,0)),(pos1-SPeak(Peak1,3,1)) )   ! Wl
            Peak0=Peak0+1
            Peak1=Peak1+1
         Else ! Peaks are too far apart to count as common peak
            Peak1=Peak1+1
         Endif
      endif
   Enddo
   !write(*,*) 'PeakD_nr=',Peak2
   !PeakSP(:,2)=PeakDP(:)
   PeakD_nr=Peak2
   SPeak(PeakD_nr+1:PeakS_dim,1,0)=0     ! set-up for sorting on pulse-strength
   Call Double_sort(SPeak(:,:,0))
   Do Peak2=1,PeakD_nr
      PeakSP(Peak2,2)=SPeak(Peak2,4,0)
      PeakSWu(Peak2,2)=SPeak(Peak2,2,0)
      PeakSWl(Peak2,2)=SPeak(Peak2,3,0)
      PeakSAmp(Peak2,2)=-SPeak(Peak2,1,0)/2
   EndDo
   !
   If(.not.production) Write(2,"(A)", ADVANCE='NO') 'PeakSP(double)='
   If(.not.production) Write(2,"(20i6)") PeakSP(1:PeakD_nr,2)
   If(.not.production) Write(2,"(A)", ADVANCE='NO') 'Peakval(Dual)='
   If(.not.production) Write(2,"(20I6)") -SPeak(1:PeakD_nr,1,0)/2
   If(.not.production) Write(2,"(A,i3)", ADVANCE='NO') 'PeakSWl(Dual)='
   If(.not.production) Write(2,"(20i6)") PeakSWl(1:PeakD_nr,2)
   If(.not.production) Write(2,"(A,i3)", ADVANCE='NO') 'PeakSWu(Dual)='
   If(.not.production) Write(2,"(20i6)") PeakSWu(1:PeakD_nr,2)
   !Write(2,"(A,i3)", ADVANCE='NO') 'ordered pol0='
   !Write(2,"(20i6)") SPeak(:,1,0)
   !Write(2,"(A,i3)", ADVANCE='NO') 'ordered pol1='
   !Write(2,"(20i6)") SPeak(:,1,1)
   !write(*,*) 'end DualPeakFind',PeakD_nr
   If(.not.production) Write(2,*) ' PeakD_nr',PeakD_nr
   Return
End Subroutine DualPeakFind
!==================
Subroutine Inv5StarPk(HEnvel, PeakPos, Windw, AvePos, Ampl, StDev)
!  Unvestigate 5-star peak structure, i.e. good for calibration
   use constants, only : dp
   Implicit none
   Real(dp), intent(in) :: HEnvel(*)
   Integer, intent(in) :: PeakPos, Windw  ! Peaks positions
   Real(dp), intent(out) :: AvePos, Ampl, StDev
   integer :: i
   real(dp) :: AvePosSq,W, TotAmpl
   AvePos=0. ; AvePosSq=0.  ; TotAmpl=0.
   Do i=-Windw,Windw
      W=HEnvel(i+PeakPos)**2
      AvePos=AvePos + i*W
      AvePosSq=AvePosSq + i*i*W
      TotAmpl=TotAmpl + W
   Enddo
   AvePos=AvePos/TotAmpl
   AvePosSq=AvePosSq/TotAmpl
   StDev= sqrt(AvePosSq-AvePos**2)
   Ampl=sqrt(TotAmpl/(2*Windw+1))
   If(StDev.lt.25) Then
      write(2,"(A,I6,A,F5.1,A,F5.1,A,F7.0,A,F7.0,A,F5.2,A,I4)") '*****',PeakPos,' position offset=',AvePos,', StDev=',StDev, &
         ', Peak Ampl=',sqrt(TotAmpl/StDev), ' or',HEnvel(PeakPos),' ratio:',HEnvel(PeakPos)/sqrt(TotAmpl/StDev), &
         ', halfwindow=',Windw
   endif
End Subroutine Inv5StarPk
!=============================
Subroutine Conv_Station_ID2i(Station_ID,i_stat)
   use Chunk_AntInfo, only : Unique_StatID, Nr_UniqueStat
   Implicit none
   integer, intent(in) :: Station_ID
   integer, intent(out) :: i_stat
   integer :: i
   i_stat=1
   Do i=1, Nr_UniqueStat      ! Get station number from the Unique_StatID list
       If(Unique_StatID(i).eq. Station_ID) Then
          i_stat=i
          exit
       EndIf
   enddo
   Return
End Subroutine Conv_Station_ID2i
