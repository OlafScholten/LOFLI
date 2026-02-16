Subroutine BuildCC(DistMax)
! Build all the Cross-correlation traces
!  v16:  True real correlation included
!  v16: self correlation included in the fit and its peak in time related to RefAntErr
!  v20: Include Polariz mode where cross correlations are calculated for the t (i_tp=1=Theta) and p (i_tp=0=phi) modes
!        This requires treating the even and odd antennas at the same level and it will be assumed that the
!        i_eo=0 (even=Y-dipoles) peak positions hold for the two antenna orientations.
    use DataConstants, only : ChunkNr_dim, DataFolder, RunMode  ! , Production
    !use DataConstants, only : Polariz
    use Chunk_AntInfo, only : Ant_Stations, Ant_IDs, Ant_nr, Ant_pos, Nr_UniqueStat, RefAnt, Nr_UniqueAnt, Unique_SAI
    use ThisSource, only : Nr_Corr, Peakpos, Stat_pos, TotPeakNr  ! , PlotCCPhase, CCPhase_ave
    use ThisSource, only : T2_dim, PeakNrTotal  !
    use constants, only : dp
    use FitParams, only : Fit_TimeOffsetAnt ! ImagingRun,
    use FFT, only : RFTransform_su,DAssignFFT, RFTransform_CF2CT
    !Use Interferom_Pars, only : dtAnt_tp, NSrc_tp
    Implicit none
    !integer, intent(in) :: StatMax
    Real(dp), intent(in) :: DistMax
    real ( kind = 8 ) :: Dist
    integer :: i_loc(1), i_ant, j_corr, i_eo, i_chunk, Station, i_Peak, ReferenceAnt, i_eo1, i_SAI, i_type
    character(len=2) :: txt
    logical, save :: first=.true., RefAntCase
    !
   !
   !
   i_eo1=1
   !If(Polariz) Then
   !   i_eo1=0
   !EndIf
   !dtAnt_tp(:)=0.
   !NSrc_tp(:)=0
   If(first) then
      Call AntFieParGen()
      first=.false.
   endif
    !
    call RFTransform_su(T2_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !write(2,*) '!entering BuildCC:', T2_dim
    !flush(unit=2)
    Nr_Corr(:,:)=0
    Do i_chunk=1, ChunkNr_dim
        Do i_eo=0,i_eo1
            Station=0
            j_corr=0
            i_type = -1 !MOD(Ant_Stations(1,i_chunk),10)
            !write(2,*) '!BuildCC:', i_type
            Do i_ant=1,Ant_nr(i_chunk)
               ! N.B. It loops over all antenna's not knowing the polarity of the peak
               If(i_type .ne. MOD(Ant_Stations(i_ant,i_chunk),10) ) Then
                  i_type = MOD(Ant_Stations(i_ant,i_chunk),10)
                  ReferenceAnt=RefAnt(i_chunk,i_eo, i_type)
                  !write(2,*) '!BuildCC, refant-set', i_ant, ReferenceAnt, i_type, J_Corr
                  RefAntCase=.true.
                  Call GetCorrSingAnt(ReferenceAnt, J_Corr, i_eo, i_chunk, RefAntCase) ! first call with ref ant
  !                J_RefCorr(i_type)=J_Corr
                  RefAntCase=.false.
               EndIf
               If(i_ant.eq.ReferenceAnt) cycle
               Dist=sqrt(sum((Ant_pos(:,i_ant,i_chunk)-Ant_pos(:,ReferenceAnt,i_chunk))**2))  ! should be ReferenceAnt and not 1 ???
               If(Dist .gt. DistMax*1000.) cycle  ! Limit to antennas near the superterp (dist in [m], DistMax in [km])
               if(mod(Ant_IDs(i_ant,i_chunk),2) .ne. i_eo) cycle       ! keep antenna orientation
               !
               Call GetCorrSingAnt( i_ant, J_Corr, i_eo, i_chunk, RefAntCase) ! will infact be for a coupl when "Polariz=.true."
               !
            Enddo
            ! write(2,*) 'Nr_Corr(i_eo,i_chunk):',i_eo,i_chunk,j_corr
            ! flush(unit=2)
            Nr_Corr(i_eo,i_chunk)=j_corr
            !If(PlotCCPhase) then
            !   Close(UNIT=30)
            !endif
        Enddo
    EndDo  !  i_chunk=1, ChunkNr_dim
    !
    call DAssignFFT()                   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !stop 'BuildCC'
    !
   !
   Return
    !stop 'BuildCC'
End Subroutine BuildCC
!=====================================
Subroutine GetCorrSingAnt( i_ant, J_Corr, i_eo, i_chunk, RefAntCase)
!   Get the Correlation for a single antenna and all peaks in this chunk
! Peak_Offst(i_Peak)= The offset in timing [LBA-samples] between LOFAR-core and reference antenna (depends on source position)
! PeakPos(i_Peak) = position [LBA-samples] of this source in data-chunk for the reference station (not dependent on source position)
   use DataConstants, only : PeakNr_dim,Time_dim, Production
   use Chunk_AntInfo, only : CTime_spectr, CTime_Hspectr, Ant_Pos, Ant_RawSourceDist, Ant_Stations, Ant_IDs
   use Chunk_AntInfo, only : Unique_SAI, Unique_StatID, Nr_UniqueStat, Tot_UniqueAnt
   use Chunk_AntInfo, only : Ant_nr ! , RefAnt
   use ThisSource, only : CorrAntNrs, SourcePos, Stat_pos, PulseDimHL  ! , Peak_IRef, Peak_RefSAI
   use ThisSource, only : Peak_Offst, T_Offset, T_Offset_ms
   use ThisSource, only : CCorr_Err, CCorr_tSh, CnuRt, CCorr_L, CCorr_H, ImpImWeights !, Error_norm, CCNorm, CCPhase_ave, RealCorrelation, PlotCCPhase
   use ThisSource, only : SafeHL, T2_Dim, prntCCNrm, TotPeakNr, PeakPos, RefAntErr !
   use ThisSource, only : CCShapeCut, CCorr_Val, CC_WidRef
   use FitParams, only : MeanCircPol, i_SAI
   use FitParams, only : Fit_AntOffset, Fit_TimeOffsetAnt, Fit_TimeOffsetStat, PulsPosCore  ! N_FitPar_max,
   use constants, only : dp,pi,ci,sample
   use StationMnemonics, only : Station_ID2Mnem, Statn_ID2Mnem
   Use Interferom_Pars, only : dtAnt_tp, NSrc_tp
   use Chunk_AntInfo, only : LHBASet
   Implicit none
   Integer, intent(inout) :: J_Corr
   Integer, intent(in) ::  i_ant, i_eo, i_chunk
   Logical, intent(in) :: RefAntCase
   !
   Real(dp) :: RDist, SubSample_Offset, StatFineOff, sample_ms, CC_Wid, CC_val, RtMax, SearchRange, Error, Norm
   Real(dp) :: ACCorr(-SafeHL(1):SafeHL(1)) !, CCPhase(1:PeakNr_dim)  ! SafeHL(1) to fit HBA also
   !Real(dp), save :: CC_WidRef(1:PeakNr_dim)
   Integer ::  i_Peak, Sample_Offset, i, i_type, i_RefAnt, SafeRange, toSampl
   Complex(dp) ::  CrCor(T2_dim)
   !Complex(dp), save :: CnuRt(0:T2_dim/2,1:3,1:PeakNr_dim) !  Fourier transform for reference antenna
   integer :: i_loc(1), StLoc, t_Max,i_stat, Antenna_SAI  !i_SAI,
   Integer, save :: PeakNr1, i_stat_old=0, N_ant, RefAntenna
   Integer, save :: j_corr_L, j_corr_H ! highest j_corr value for this antenna type
   character(len=6) :: Station_Mnem
   !    Real(dp), external :: CCorr_der
   logical :: prnt= .false.
   !logical :: prnt= .true.
   Complex(dp), Pointer :: CTspec_a_c(:)
   !
   !=============================
   !Safety is used also in [SourceTryal.f90] & [GLE plotting of correlations]
   !   128=T2_dim = Tref_dim + SafeHL(1)*2
   !Real(dp) :: ( CCNorm(Ant_nrMax,1:PeakNr_dim)),
   !Real(dp), save, allocatable  ::  CCorr(-Safety:Safety,Ant_nrMax,1:PeakNr_dim))
   !============================================
   !
   If(j_corr.eq.0) Then
      j_corr_L=0
      j_corr_H=0
   EndIf
   j_corr=j_corr+1
   CorrAntNrs(j_corr,i_eo, i_chunk)=i_ant
   if(j_corr.eq. 1) then
     If((i_eo.eq.0) .and. (i_chunk.eq.1)) then
         PeakNr1=1
     elseif(i_eo.eq.0) then
         PeakNr1=TotPeakNr(1,i_chunk-1) + 1
     Else
         PeakNr1=TotPeakNr(0,i_chunk) + 1
     Endif
   Endif
   If(PeakNr1.gt.TotPeakNr(i_eo,i_chunk)) return
   !
   Do i_stat=1,Nr_UniqueStat  ! Get station number from unique list to retrieve timing-offset
     If(Unique_StatID(i_stat).eq. Ant_Stations(i_ant,i_chunk)) exit
   Enddo
   StatFineOff=Fit_TimeOffsetStat(i_stat)
   Antenna_SAI=100*Ant_Stations(i_ant,i_chunk) + Ant_IDs(i_ant,i_chunk)
   Do i_SAI=Tot_UniqueAnt(i_stat-1)+1,Tot_UniqueAnt(i_stat)      ! Get antenna number from the Unique_Antennas list
       If(Unique_SAI(i_SAI).eq. Antenna_SAI) exit
   enddo
   If(Fit_AntOffset) then
      Do i_SAI=Tot_UniqueAnt(i_stat-1)+1,Tot_UniqueAnt(i_stat)      ! Get antenna number from the Unique_Antennas list
          If(Unique_SAI(i_SAI).eq. Antenna_SAI) goto 1
      enddo
      write(*,*) 'unknown antenna in GetCorrSingAnt: ',Antenna_SAI
      write(2,*) 'unknown antenna in GetCorrSingAnt: ',Antenna_SAI, Unique_StatID(i_stat), Ant_Stations(i_ant,i_chunk)
      write(2,*) i_stat,i_ant, Tot_UniqueAnt(i_stat-1)+1,Tot_UniqueAnt(i_stat)
      write(2,*) Unique_SAI(Tot_UniqueAnt(i_stat-1)+1:Tot_UniqueAnt(i_stat))
      stop 'GetCorrSingAnt'
    1 Continue
      StatFineOff=Fit_TimeOffsetStat(i_stat) + Fit_TimeOffsetAnt(i_SAI)
   endif
   Call LHBASet(i_ant, i_chunk, CTspec_a_c, i_type, Sample_ms, toSampl)
   SafeRange=SafeHL(i_type)
   If(i_type.eq.0) Then
      j_corr_L=j_corr_L+1
   Else
      j_corr_H=j_corr_H+1
   EndIf
   !
   !write(2,*) 'in GetCorrSingAnt, PeakNr1=',PeakNr1,i_eo,i_chunk,j_corr, i_ant, RefAntCase,PeakPos(PeakNr1)
   !flush(unit=2)
   !Production=.false.
   Do i_Peak=PeakNr1,TotPeakNr(i_eo,i_chunk) ! This assumes that the loop over antenna number goes fastest (constant [i_eo,i_chunk])
      Call RelDist(SourcePos(1,i_Peak),Ant_pos(1,i_ant,i_chunk),RDist)
      Rdist=Rdist - Ant_RawSourceDist(i_ant,i_chunk) + StatFineOff ! - INT(Peak_Offst(i_Peak))
      If(j_corr.eq.1) Then
         Peak_Offst(i_Peak)=(Rdist + RefAntErr(i_Peak))  !=T_shft1 in 'CompareCorrTime' !  In units of [LBA-samples]
      end If
      ! Correction is needed because pulse position is found in reference antenna, while all is usually calculated w.r.t. the core
      If(RefAntCase) then ! all is calculated w.r.t. first antenna
         If(.not. Production) then
            Write(2,"(i3,2i2,A,I8,A,3F11.1,A,F6.2,A)") i_Peak,i_eo,i_chunk,', PeakPos=',PeakPos(i_Peak), &
                 ', source-pos=(',SourcePos(1:3,i_Peak),')[m] , RefAntErr=', RefAntErr(i_Peak),'[samples]'
         Endif ! (.not. Production)
         StLoc=toSampl*PeakPos(i_Peak) - PulseDimHL(i_type)/2  ! integer sample number for LBA
         If(PulsPosCore) StLoc=StLoc+ NINT(toSampl*Peak_Offst(i_Peak)) ! Peak_Offset corrects for ref-antanna not at the center of CS002;
         !write(2,*) '!GetCorrSingAnt, refant:',i_ant,RDist,Ant_RawSourceDist(i_ant,i_chunk),(StLoc+ PulseDimHL(i_type)/2)/toSampl
         !  Correction should be negligible for ref antennas in CS002 since the position shift is crudely corrected for by RawSourceDist
         !write(2,*) 'StLoc=',StLoc
         !flush(unit=2)
         RefAntenna=i_ant
         Call GetAnt(CTspec_a_c, i_ant, i_chunk, StLoc, PulseDimHL(i_type), T2_dim, CnuRt(0,1,i_peak), Norm) ! get nu-spectrum for reference antenna
         CnuRt(:,1,i_peak)=CnuRt(:,1,i_peak)/Norm
         !write(2,*) 'ref:',i_ant,i_peak, StLoc,  T2_dim, SubSample_Offset
      EndIf  ! (RefAntCase)
      !write(2,*) 'start corr',i_ant, J_Corr, i_eo, i_chunk
      !    flush(unit=2)
      !
      RDist=RDist - Peak_Offst(i_Peak) !  In units of [LBA-samples]
!      T_Offset_ms(j_corr,i_Peak)=RDist*LBAs_ms  !  All station delays are presumably already included when extracting the time-trace
      T_Offset(j_corr,i_Peak)=RDist   !  In units of [LBA-samples]
      ! where 'T_Offset' is the shift between the this and the reference antenna-timings which is
      !   already accounted for when calculating the cross correlations,
      !   which, of course, is based on the source positions before the coming fitting round.
      !  T_Offset  is used in fitting routine to determine time-shifts due to source-position shifts (and other)
      ! RDist is measured in LBA samples
      RDist=toSampl*RDist  ! convert to HBA samples when needed
      Sample_Offset = INT(RDist) ! in units of sample size for LBA/HBA
      SubSample_Offset = RDist - Sample_Offset ! in units of sample size
      !
      !write(2,*) 'Sample_Offset',J_corr, i_peak, Sample_Offset, SubSample_Offset
      StLoc=toSampl*PeakPos(i_Peak) + Sample_Offset - T2_dim/2  ! Check - sign 'Sample_Offset'
      If(PulsPosCore) StLoc=StLoc+ NINT(Peak_Offst(i_Peak)) ! Peak_Offset corrects for ref-antanna not at the center of CS002;
      !
      !If(j_corr.gt. 140)        write(2,*) 'I_ant=71936:',j_corr, i_ant, PulsPosCore, StLoc, Time_dim, T2_dim
      If( (StLoc.lt.100) .or. (StLoc.gt.(toSampl*Time_dim-T2_dim-100) )) then
         If(.not. Production) &
            write(2,*) '!****StLoc=',StLoc,' reaches [below lower/ above upper] limit for peak', &
               i_Peak,' in ',Statn_ID2Mnem(Unique_StatID(i_stat)), Rdist, Ant_RawSourceDist(i_ant,i_chunk), StatFineOff
         Error=2000. ! Huge and thus do not include in fitting
         CCorr_Val(j_corr,i_Peak) = 0.
         CCorr_tSh(j_corr,i_Peak) = 0. ! in units [LBA-samples]
         ACCorr(:)=0.
      Else
         ! write(2,*) '! Before  CrossCorr', SearchRange, i_ant, j_corr, i_Peak, Peakpos(i_Peak), SourcePos(1,i_Peak)
         Call CrossCorr(CTspec_a_c, CnuRt(0,1,i_peak), i_ant, i_chunk, StLoc, T2_dim, SubSample_Offset, CrCor)
         !If(i_ant.lt.5) write(2,*) i_ant,i_peak,j_corr, StLoc, T2_dim, SubSample_Offset,CrCor(T2_dim),CrCor(1:2)
         ! Checked: for the self correlation the max of the real part is (almost) at 1 (delta-t=0, where you would expect it,
         !  however the max of the hilbert envelope may be moved from zero. Seems that for asymmetric pulses the max&min of the
         !  imaginary part is not symmetrix, which may be the reason that the Hilbert transform peaks elsewhere than delta-t=0
         ! For a true self correlation this cannot happen, however the reference-antenna spectrum lenght 'PulseDimHL(i_type)' is padded
         !  with zeros to reach the length 'T2_dim' of the spectra for which the cross correlations are calculated. This
         !  difference in length is the culprit.
         !
         Call SearchWin(RefAntenna, i_ant, i_chunk, SourcePos(1,i_Peak), SearchRange)
         SearchRange=toSampl*SearchRange
         ! write(2,*) '! Before  CrossCorr_Max', SearchRange, i_ant, j_corr, i_Peak, SourcePos(1,i_Peak)
         !    flush(unit=2)
         If(RefAntCase) SearchRange=SafeRange !  for the reference antenna in units of LBA samples
         Call CrossCorr_Max(i_ant,i_chunk, i_Peak, CrCor, ACCorr, RtMax, CC_Wid, CC_val, SearchRange, Error, SafeRange)
   !      CCorr(:,j_corr,i_Peak)=ACCorr(:)/CC_val  !  /maxval(CCorr(:,j_corr,i_Peak)); used for plotting only
         !       stop
         error=error
         RtMax = RtMax/toSampl  ! convert to units [LBA-samples]
         CC_Wid= CC_Wid/toSampl  ! convert to units [LBA-samples]
         If(j_corr .eq. 1) then
            RtMax=RtMax-RefAntErr(i_Peak)  ! in units [LBA-samples]
            CC_WidRef(i_Peak)=CC_Wid  ! set in call to CC_Max
            !write(2,"(2I3,2F8.2,A,11(i2,':(',G12.3,'), '))") i_ant, i_Peak, RtMax, RDist,' ABScrcor:', &
            !   ( i-1,ABS(crcor(T2_dim+i)),i=-2,0,1), (i-1,ABS(crcor(i)), i=1,4)
         Endif
         CC_Wid=CC_Wid/CC_WidRef(i_Peak) ! relative width compared to that of the auto correlation
         CCorr_Val(j_corr,i_Peak) = CC_val
         CCorr_tSh(j_corr,i_Peak) = RtMax ! in units [LBA-samples]
         !
         If(CC_wid .lt. 1) CC_wid=1./CC_wid  ! too small is also not good
         If(ImpImWeights) error=error*0.5*(CC_Wid/CC_val)  ! Est. Error in time-shift in units [LBA samples]; squaring seems not to make a significant difference.
         !If(ImpImWeights) error=error*0.5*(CC_Wid/CC_val)**2  ! square gives huge values; Est. Error in time-shift in units [LBA samples]; squaring seems not to make a significant difference.
         !If((CC_Wid.lt. CCShapeCut) .or. (CC_Wid.gt. 1./CCShapeCut)) then
         !   Error=Error+30
         !   If(.not. Production) write(2,*) 'Quality:',Error, CC_Wid, i_ant, &
         !      Statn_ID2Mnem(Ant_Stations(i_ant,i_chunk)), Ant_IDs(i_ant,i_chunk)
         !   write(2,"(A,I4,I3,A7,I3,', delta_t=',F7.2,', RelWidth=',F6.2,', norm=',F7.3, 2F9.3)") '! **** Quality:',j_corr, &
         !      i_Peak,Statn_ID2Mnem(Ant_Stations(i_ant,i_chunk)), Ant_IDs(i_ant,i_chunk), RtMax, CC_Wid, CC_val,Error,CCShapeCut
         !Endif
         If(((CC_Wid .gt. 1.6) .or. (CC_val.lt. 0.6)  ) .and. prnt ) Then
            write(2,"(A,I4,I3,A7,I3,', delta_t=',F7.2,', RelWidth=',F6.2,', norm=',F7.3, 2F9.3)") '! Poor Quality:',j_corr, &
               i_Peak,Statn_ID2Mnem(Ant_Stations(i_ant,i_chunk)), Ant_IDs(i_ant,i_chunk), RtMax, CC_Wid, CC_val, Error
         EndIf
      EndIf
      CCorr_Err(j_corr,i_Peak) = Error  ! in units [LBA-samples]
      If(i_type.eq.0) Then
         CCorr_L(-SafeHL(0):SafeHL(0),j_corr_L,i_Peak)=ACCorr(-SafeHL(0):SafeHL(0))/CC_val
      Else
         CCorr_H(:,j_corr_H,i_Peak)=ACCorr(:)/CC_val
      EndIf
      !
      If( (i_Peak.eq.1) .and. prnt ) then
            Call Station_ID2Mnem(Ant_Stations(i_ant,i_chunk),Station_Mnem)
            write(2,"(A,I3,A,i3,2x,A5,A,i2.2,A,F8.2,A,I2,A,F7.2,A,I7,f12.5)") '!GetCorrSingAnt, j_corr=',J_corr, &
                ', i_ant=',i_ant, Station_Mnem,'-',Ant_IDs(i_ant,i_chunk), ', offset=',RDist,', i_eo=',i_eo, &
                ', RtMax=',RtMax,', startLoc=',StLoc + T2_dim/2,Ant_RawSourceDist(i_ant,i_chunk)
            !write(2,*) Sample_Offset,StLoc+T2_dim/2
        endif
        !write(2,"(10g10.3,' ; ',g10.3,' ; ',10g10.3,' --')") RTTrace2(T2_dim/2-10:T2_dim/2+10) ! Self correlation, correct for upsampling
        !write(2,"(10g10.3,' ; ',g10.3,' ; ',10g10.3,' --')") &
        !    abs(CTime_spectr(StLoc + T2_dim/2-10:StLoc + T2_dim/2+10,i_ant)) ! Self correlation, correct for upsampling
        !If(prnt) write(2,"(i2,10f7.3,' ; ',f7.3,' ; ',10f7.3,' --')") i_Peak,CCorr(-10:10,j_corr,i_Peak) ! Self correlation, correct for upsampling
        ! If(Error.gt. 1. ) Write(2,*) 'GetCorrSingAnt:',j_corr,Antenna_SAI,i_Peak,Error
        !If(j_corr.eq.10) write(2,"(A,I3,3F6.2, 10(F5.1,';'))") '!GetCorrSingAnt, error',i_Peak, CC_Wid, CC_val, CC_Wid/CC_val, &
        !    CCorr_Err(1:10,i_Peak)
    enddo  ! ,i_Peak
    !
    NULLIFY( CTspec_a_c )
    !
    Return
End Subroutine GetCorrSingAnt
!=================================
! ==========================! ==========================
! ==========================! ==========================
Subroutine CrossCorr(CTspec_a_c, CnuR, i_ant, i_chunk, StLoc, T2_dim, t_shft, CrCor)
! all timings are in units of samples for this particular antenna type
    use FFT, only : RFTransform_CF2CT
    use constants, only : dp,pi,ci
    Implicit none
    integer, parameter :: UpSamplFact=1
    integer, parameter :: HW_size=5
    !integer, parameter :: HW_size=3  ! changed on Aug 30, 2020
    Complex(dp), intent(in) :: CTspec_a_c(*)
    Integer, Intent(in) ::  i_ant, i_chunk, StLoc, T2_dim  ! T2_dim is small window around peak
    Real(dp), intent(in) :: t_shft
    Complex(dp), intent(in) :: CnuR(0:T2_dim/2)  ! Frequence spectrum reference peak
    Complex(dp), intent(out) :: CrCor(T2_dim)
    integer :: i
    Complex(dp) :: Cnu(0:T2_dim/2)
    complex(dp), parameter :: tipi=2*ci*pi
    Real(dp) :: Norm
    !
    Call GetAnt(CTspec_a_c, i_ant, i_chunk, StLoc, T2_dim, T2_dim, Cnu, Norm)
    !
!    Cnu_s(0:T2_dim/2)=Cnu2_s(0:T2_dim/2) * conjg(Cnu1_s(0:T2_dim/2))
    DO i=0,T2_dim/2
        Cnu(i)=Cnu(i) * conjg(CnuR(i)) * exp(-tipi*t_shft*i/T2_dim) /Norm
    ENDDO
    !
    Call RFTransform_CF2CT(Cnu,CrCor )
    !call DAssignFFT()                   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !
    Return
End Subroutine CrossCorr
! ==========================! ==========================
! ==========================! ==========================
Subroutine GetAnt(CTspec_a_c, i_ant, i_chunk, StLoc, T_dim, T2_dim, Cnu, Norm)
! all timings are in units of samples for this particular antenna type
!  T_dim : dimension oof window around pulse (from StLoc till StLoc+t_dim), will be centered in middle of T2_dim,
!  T2_dim : dimension larger array used for cross correlations
   use FFT, only : RFTransform_CF
   !use Chunk_AntInfo, only : CTime_spectr
   use constants, only : dp,pi,ci
   Implicit none
    Complex(dp), intent(in) :: CTspec_a_c(*)  ! = CTime_spectr(*,i_ant,  i_chunk) from "use Chunk_AntInfo, only : CTime_spectr"
   Integer, intent(in) :: i_ant, i_chunk, StLoc, T_dim, T2_dim
   ! Real(dp), intent(in) :: RTTrace(T_dim)
    Complex(dp), intent(out) :: Cnu(0:T2_dim/2)
    Real(dp), intent(out) :: Norm
    integer, parameter :: HW_size=5
    Real(dp) :: Hi, RTime(T2_dim)
    integer :: i, strt
   !
   strt=(T2_dim-T_dim)/2
   If(strt.lt.0) then
     write(2,*) '*******strt1=',strt
     stop 'GetAnt'
   endif
   RTime(1:strt+1)=0.  ! padding with zeros
   Rtime(strt+HW_size:strt+T_dim-HW_size+1) = Real(CTspec_a_c(StLoc+HW_size:StLoc+T_dim+1-HW_size))
   !write(2,*) 'Rtime(strt+HW_size:strt+T_dim-HW_size+1)',strt,Rtime(strt+HW_size:strt+T_dim-HW_size+1)
   RTime(strt+T_dim:T2_dim)=0.
   Do i=1,HW_size-1
     Hi=sin(0.5*(i-0.5)*pi/HW_size)**2 ! Calculate Hann window
     RTime(strt+i) = Hi*Real(CTspec_a_c(StLoc+i))
     RTime(strt+T_dim+1-i)= Hi*Real(CTspec_a_c(StLoc+T_dim+1-i))
   EndDo
   RTime(strt+T_dim+1:T2_dim)=0.
   !
   Call RFTransform_CF(RTime(1),Cnu(0))
   Norm= sqrt(SUM( Cnu(:)*conjg(Cnu(:)) ))
   Return
End Subroutine GetAnt
! ==========================! ==========================
! ==========================! ==========================
!=================================
Subroutine CrossCorr_Max(i_ant,i_chunk, i_Peak, CrCor, TrueCC, RtMax, CC_Wid, CC_val, &
      SearchRange, Error, Safety)
! all timings are in units of samples for this particular antenna type
!   Get the value of real and imaginary parts of the X-Correlation at the position of the Max in abs. value
!   Use a dynamic window to search for pulses
   !use Chunk_AntInfo
   use ThisSource, only : T2_dim, SafeHL, ExclStatNr, t_CCorr, Error_norm
   use ThisSource, only : RealCorrelation
   !use ThisSource, only : PlotCCPhase
   use Chunk_AntInfo, only : Ant_IDs, Ant_Stations
   use DataConstants, only : Production
   use FitParams, only : SearchRangeFallOff
   use constants, only : dp ,pi !,ci,sample
   use StationMnemonics, only : Statn_ID2Mnem
   Implicit none
   !Integer, intent(in) ::  i_ant, i_chunk, j_corr,i_Peak
   Complex(dp), intent(in) ::  CrCor(1:T2_dim)
   Integer, intent(in) :: i_ant,i_chunk,i_Peak, Safety
   Real(dp), intent(in) :: SearchRange
   Real(dp), intent(out) :: TrueCC(-SafeHL(1):SafeHL(1)) ! as declared in calling
   Real(dp), intent(out) :: Error, RtMax, CC_Wid, CC_Val !, CCPhase
   !
   integer :: i_loc(1), t_Max, i,k, Range, j
   !
   Real(dp) :: TrueCC_pp(-Safety:Safety), RCCorr_pp(-Safety:Safety), ICCorr_pp(-Safety:Safety)
   Real(dp) :: RCCorr(-Safety:Safety), ICCorr(-Safety:Safety), Max_Time(1:2*Safety), Max_Val(1:2*Safety)
   Real(dp) :: W, Ws, Wx, Wxx
   Real(dp), save :: tA,B,Yp, dt, Rval, Ival
   !
   !write(2,*) '! Entering CrossCorr_Max', i_ant,i_chunk, i_Peak, SearchRange, Safety !, ', CrCor=', CrCor
   !flush(unit=2)
   If(RealCorrelation) then
      TrueCC(0:Safety) = Real(CrCor(1:Safety+1))
      TrueCC(-Safety:-1) = Real(CrCor(T2_dim-Safety+1:T2_dim))  ! correct for upsampling
   else
      TrueCC(0:Safety) = abs(CrCor(1:Safety+1))
      TrueCC(-Safety:-1) = abs(CrCor(T2_dim-Safety+1:T2_dim))  ! correct for upsampling
   endif
   !
   Range=NINT(SearchRange*SearchRangeFallOff)
   if(Range.gt.Safety) Range=Safety
   CC_Val = MAXVAL( TrueCC(-Range:Range) )
   !
   Do i=1,Safety ! (flat part upto error, lin decreasing beyond) is changed to parabola
      B=(1. - i*i/(SearchRange*SearchRange*SearchRangeFallOff*SearchRangeFallOff))
      TrueCC(i) = TrueCC(i)*B
      TrueCC(-i) = TrueCC(-i)*B
   Enddo
   !
   !CCPhase = 0.
   Error=0.
   If(count((ExclStatNr(:,i_peak)-Ant_Stations(i_ant,i_chunk)).eq.0,1) .ge. 1) then
      !write(2,*) 'excluded station=',Ant_Stations(i_ant,i_chunk),i_ant, ', for i_peak=',i_peak
      RtMax=0.
      Error=2000
      Return
   EndIf
   !
!   Call spline_cubic_set( 2*Safety+1, t_ccorr(-Safety), TrueCC(-Safety), TrueCC_pp(-Safety) )
   If(Range.le.0) Then
      write(2,*) 'CrossCorr_Max; pline_cubic_set, Range problems', Range, SearchRange
      Flush(unit=2)
   EndIf
   !
   Call spline_cubic_set( 2*Range+1, t_ccorr(-Range), TrueCC(-Range), TrueCC_pp(-Range) )
   !
   If(RealCorrelation) then ! first find all local maxima, then get the global one
      CC_Wid = sqrt(SUM( TrueCC(-Range:Range)*TrueCC(-Range:Range) )) / CC_val
      j=1
      Max_Time(j)=-Range
      Max_Val(j) = TrueCC(-Range)
      Do i=-Range,Range-1
         Yp=TrueCC(i+1)-TrueCC(i) - TrueCC_pp(i+1)/6. -TrueCC_pp(i)/3.    !=y'(left=i)
         if(Yp.lt.0.) then
             cycle
         endif
         tA=TrueCC_pp(i+1) - TrueCC_pp(i)   ! 2 * A
         B=TrueCC_pp(i)
         CC_Val=B*B - 2.*tA * Yp
         If(CC_Val.le.0.) cycle
         CC_Val=sqrt(CC_Val)
         dt= ( -B - CC_Val )/tA
         !If(abs(tA).lt.1.e-20) dt=-yp/B
         !write(2,*) i,'yp:',yp,( -B - sqrt(B*B - 2.*tA * Yp) )/tA,( -B + sqrt(B*B - 2.*tA * Yp) )/tA
         If(dt.lt.0. .or. dt.gt.1) then
            dt= ( -B + CC_Val )/tA
            !If(dt.lt.0. .or. dt.gt.1) then
            If(dt.gt.0. .and. dt.lt.1) then
               write(2,*) 'dt',dt,yp,B,tA,CC_Val,'found THE exception'
               cycle
            Endif
         Endif
         CC_Val= TrueCC(i) &
             + dt*( Yp + dt*( 0.5D0*TrueCC_pp(i) + dt*( (TrueCC_pp(i+1)-TrueCC_pp(i))/6.D0 ) ) )
         !If(CC_Val.lt.TrueCC(i+1)) write(2,*) 'CC_Val?;',CC_Val, TrueCC(i+1)
         If(CC_Val.lt.0.) cycle
         j=j+1
         Max_Time(j)=i+dt
         Max_Val(j) = CC_Val
      Enddo
      j=j+1
      Max_Time(j)= Range
      Max_Val(j) = TrueCC(Range)
      i_loc=MaxLoc(Max_Val(1:j) )
      RtMax=Max_Time(i_loc(1))
      t_max=NINT(RtMax)
      CC_Val=Max_Val(i_loc(1))
      i=COUNT(Max_Val(1:j).gt. 0.95*CC_Val)
      If(i.gt.1) then
         !If(.not. Production) write(2,"(5('*'),A,i3,A,i2,A,i4,' ',A)") 'For peak#', i_Peak,' there are ',i, &
         !   ' large correlation maxima found in antenna',Ant_IDs(i_ant,i_chunk),Statn_ID2Mnem(Ant_Stations(i_ant,i_chunk))
         Error=(i-1)*20.    !  Peak becomes ambguous as there are more equivalent ones
      EndIf
      !If(i_ant.le.20) then
      !   write(2,"(5x,A,3i3,A,F7.3,g11.3)") 'i_Ant & Peak:',i_ant, i_Peak,j,' max:',RtMax,CC_Val
      !   write(2,"(A,20F9.3)") 'max_time:',Max_Time(1:j)
      !   write(2,"(A,20F9.1)") 'max_val: ',Max_Val(1:j)
      !EndIf
   Else !  For the case of the Hilbert envelope of the cross correlation function
      t_max=MaxLoc(TrueCC(-Range:Range), dim=1 ) - Range -1 ! works only for a very smooth function, not for real correlation
      If((t_max.le. -Range) .or. (t_max .ge. Range)) then
         If(.not. Production) write(2,"(A,I5,I3,I4,A,A5,A,I2)") 'maximum in correlation function at bound',t_max, &
             Ant_IDs(i_ant,i_chunk),Ant_Stations(i_ant,i_chunk),'=',Statn_ID2Mnem(Ant_Stations(i_ant,i_chunk)),', peak#=',i_Peak
         RtMax=t_max
         Error=Error+300  ! time samples
         CC_Wid=range
      Else
         Yp=TrueCC(t_max+1)-TrueCC(t_max) - TrueCC_pp(t_max+1)/6. -TrueCC_pp(t_max)/3.    !CCorr_der(t_max,j_corr,i_Peak)
         if(Yp.lt.0.) then
             t_max=t_max-1  !  shift to the point left of the real maximum
             Yp=TrueCC(t_max+1)-TrueCC(t_max) - TrueCC_pp(t_max+1)/6. -TrueCC_pp(t_max)/3. !CCorr_der(t_max,j_corr,i_Peak)
         endif
         tA=TrueCC_pp(t_max+1) - TrueCC_pp(t_max)   ! 2 * A
         B=TrueCC_pp(t_max)
         dt= ( -B - sqrt(B*B - 2.*tA * Yp) )/tA
         RtMax=t_max+dt  ! position of the maximum in the correlation function
         ! for a precise determination of value of cross correlation:
         !Call spline_cubic_val( 2*Range+1, t_ccorr(-Range), TrueCC(-Range), TrueCC_pp(-Range), RtMax, CC_Val)
         ! Hhowever, precision for CC_Val is not important.
         !
         !write(2,*) '!CrossCorr_Max:RtMax=',RtMax, range
         !flush(unit=2)
      endif
      !
      ! determine width of the cross-correlation
      Ws=0.
      Wx=0.
      Wxx=0.
      j=max(t_max-15,-Range)
      k=min(t_max+15,Range)
      Do i=j,k    ! if range is large then sensitive to background
         W=TrueCC(i)*TrueCC(i)
         Ws=Ws+ W
         Wx=Wx+ i*W
         Wxx=Wxx+i*i*W
      EndDo
      Wx=Wx/Ws
      Wxx=Wxx/Ws
      CC_Wid=Wxx-Wx*Wx
      !write(2,*) '!CrossCorr_Max, wdth:', CC_Wid,Wx, Wxx, range
      !
      !If(i_ant.lt.5) write(2,*) i_ant,i_peak,t_max,'TrueCC(t_max-1:t_max+1)',TrueCC(t_max-1:t_max+1)
   EndIf
   !
   If(abs(RtMax).gt.2.*SearchRange)  Then
      Error = Error+30  ! time samples
      If(.not. Production) write(2,*) 'Large deviation:',SearchRange,RtMax, Statn_ID2Mnem(Ant_Stations(i_ant,i_chunk)),i_ant
   Else
      Error= Error+Error_norm  ! in units of LBA samples
   EndIf
   !write(2,*) 'Exit  CrossCorr_Max', RtMax, CC_Val, CCPhase, Error, '; TrueCC=', TrueCC
   !flush(unit=2)
   !
   Return
End Subroutine CrossCorr_Max
!=====================================
!=================================
!=====================================
Subroutine GetRefAnt(ChunkNr)
   ! with the dual option the relative distance between the two polarity antennas should be small
   use DataConstants, only : ChunkNr_dim
!   use DataConstants, only : Polariz
   use Chunk_AntInfo, only : Ant_IDs, Ant_pos, RefAnt, LBA_nr, HBA_nr, Ant_Stations
   use ThisSource, only : Dual
   use constants, only : dp
   use unque, only : Double_sort
    use StationMnemonics, only : Statn_ID2Mnem
   Implicit none
   Integer, intent(in) :: ChunkNr
   Integer, parameter :: MaxRefAnt=20, i_e=0, i_o=1, MA=100  ! max nr of antenna pairs tested for close proximity
   Integer :: Ant(1:Ma,0:2), j, i_ante, i_anto, i_chunk, N_chunk, i_type
   Integer :: Antnr(0:1,1:2)
   Real(dp) :: Distnr(0:1,1:2), Dist
   If((ChunkNr .le. 0) .or. (ChunkNr .gt. ChunkNr_dim)) Then
      N_chunk=ChunkNr_dim
   Else
      N_chunk=ChunkNr
   EndIf
   RefAnt(:,:,:)=0
   Do i_chunk=1, N_chunk
      Antnr(0,1)=1  ! Antnr(i_type, min/max)
      Antnr(0,2)=min(LBA_nr(i_chunk),MaxRefAnt)
      Antnr(1,1)=LBA_nr(i_chunk)+1
      Antnr(1,2)=LBA_nr(i_chunk)+min(HBA_nr(i_chunk),MaxRefAnt)
  !    If(Dual) then
         Do i_type=0,1 ! MOD(Ant_Stations(1, i_chunk),10)
            j=0   ! counts the number of relative antenna distances
            If( (Antnr(i_type,2)-Antnr(i_type,1)).lt.4) cycle
            Do i_ante=Antnr(i_type,1),Antnr(i_type,2) ! build table of relative distances
               if(mod(Ant_IDs(i_ante,i_chunk),2) .ne. i_e) cycle       ! check this is an even antenna
               Do i_anto=1,MaxRefAnt
                  If(j.ge. MA ) exit
                  if(mod(Ant_IDs(i_anto,i_chunk),2) .ne. i_o) cycle       ! check this is an odd antenna
                  Dist=sqrt(sum((Ant_pos(:,i_ante,i_chunk)-Ant_pos(:,i_anto,i_chunk))**2))
                  !write(2,*) 'even:', Ant_pos(:,i_ante,i_chunk)
                  !write(2,*) 'odd: ', Ant_pos(:,i_anto,i_chunk)
                  !write(2,*) dist, j
                  !If(Dist.lt.0.001) Dist=0.
                  j=j+1
                  Ant(j,0)=NINT(10.*Dist) ! store distance between antenna pair
                  Ant(j,2)=i_ante
                  Ant(j,1)=i_anto
                  If(Dist.lt.0.5) goto 1  ! the two antennas are close enough that delays between then are less than 1 sample
               Enddo
            Enddo
            If(j.eq.0) Then
               Write(2,*) 'No reference antenna found, Dual for type', i_type
               cycle
            End If
            Call Double_sort(Ant(1:j,0:2))
            write(2,*) 'relative reference antenna distance [0.1m] and numbers:',Ant(1,:)
            j=1
         1  continue
            RefAnt(i_chunk,0,i_type)=Ant(j,2)
            RefAnt(i_chunk,1,i_type)=Ant(j,1)
            !write(2,*) 'even&odd reference antennas in station', Statn_ID2Mnem(Ant_Stations(Ant(j,2),i_chunk)), ' & ',&
            !   Statn_ID2Mnem(Ant_Stations(Ant(j,1),i_chunk))
         EndDo
   !   Else
   !      Do i_type=0,1 ! MOD(Ant_Stations(1, i_chunk),10)
   !         Do i_ante=Antnr(i_type,1),Antnr(i_type,2) ! do NOT build table of relative distances
   !            if(mod(Ant_IDs(i_ante,i_chunk),2) .ne. i_e) cycle       ! check this is an even antenna
   !            RefAnt(i_chunk,0,i_type)=i_ante
   !         EndDo
   !         Do i_anto=Antnr(i_type,1),Antnr(i_type,2)
   !            if(mod(Ant_IDs(i_anto,i_chunk),2) .ne. i_o) cycle       ! check this is an odd antenna
   !            RefAnt(i_chunk,1,i_type)=i_anto
   !         Enddo
   !      Enddo
   !      !stop 'No reference antenna found, Single'
   !   Endif  ! Dual
      !If(Dual) write(2,*) j,'reference antenna* distancebetween',Dist,'[m], and numbers:',RefAnt(i_chunk,0,i_type),RefAnt(i_chunk,1,i_type)
      !write(2,*) '!GetRefAnt2@',i_chunk,' RefAnt(i_chunk,i_eo=0:1,i_type=0:1):',RefAnt(i_chunk,0:1,0:1)
      !flush(unit=2)
   Enddo  ! i_chunk
   Return
End Subroutine GetRefAnt
!=================================
