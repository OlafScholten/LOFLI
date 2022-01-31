Subroutine BuildCC(StatMax,DistMax)
! Build all the Cross-correlation traces
!  v16:  True real correlation included
!  v16: self correlation included in the fit and its peak in time related to RefAntErr
!  v20: Include Polariz mode where cross correlations are calculated for the t (i_tp=1=Theta) and p (i_tp=0=phi) modes
!        This requires treating the even and odd antennas at the same level and it will be assumed that the
!        i_eo=0 (even=Y-dipoles) peak positions hold for the two antenna orientations.
    use DataConstants, only : ChunkNr_dim, Production, DataFolder, RunMode
    use DataConstants, only : Polariz
    use Chunk_AntInfo, only : Ant_Stations, Ant_IDs, Ant_nr, Ant_pos, Nr_UniqueStat, RefAnt, Nr_UniqueAnt, Unique_SAI
    use ThisSource, only : Nr_Corr, Peakpos, PlotCCPhase, CCPhase_ave, Stat_pos, TotPeakNr
    use ThisSource, only : T2_dim, Tref_dim, PeakNrTotal
    use constants, only : dp
    use FitParams, only : Fit_TimeOffsetAnt ! ImagingRun,
    use FFT, only : RFTransform_su,DAssignFFT, RFTransform_CF2CT
    Use Interferom_Pars, only : dtAnt_tp, NSrc_tp
    Implicit none
    integer, intent(in) :: StatMax
    Real(dp), intent(in) :: DistMax
    real ( kind = 8 ) :: Dist
    integer :: i_loc(1), i_ant, j_corr, i_eo, i_chunk, Station, i_Peak, ReferenceAnt, i_eo1, i_SAI
    character(len=2) :: txt
    logical, save :: first=.true.
    !
   !
   !Real,save :: TotCPUOut=0., TotCPUIn=0., CPUstartTime, CPUstopTime=-1.
   !
   ! beginning time recording
   !call cpu_time(CPUstartTime)
   !If(CPUstopTime.lt.0) CPUstopTime=CPUstartTime
   !TotCPUOut=TotCPUOut - CPUstopTime + CPUstartTime
   !WRITE(*,"(A,2F13.3,A)") 'Total CPU, (in/out)-side "BuildCC":',TotCPUIn, TotCPUOut, '[s]'
   !
   !
   !
   i_eo1=1
   If(Polariz) Then
      i_eo1=0
   EndIf
   dtAnt_tp(:)=0.
   NSrc_tp(:)=0
   If(first) then
      Call AntFieParGen()
      first=.false.
   endif
    If(PlotCCPhase) then
      Open(UNIT=31,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//'CCPeakPhaseAve.dat')
    endif
    !
    call RFTransform_su(T2_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !write(2,*) 'BuildCC:', T2_dim
    !flush(unit=2)
    Nr_Corr(:,:)=0
    Do i_chunk=1, ChunkNr_dim
        Do i_eo=0,i_eo1
            Station=0
            j_corr=0
            ReferenceAnt=RefAnt(i_chunk,i_eo)
            If(PlotCCPhase) then
               write(txt,"(I1,I1)") i_chunk,i_eo
               Open(UNIT=30,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//'CCPeakPhase-'//txt//'.dat')
            endif
            Call GetCorrSingAnt(ReferenceAnt, J_Corr, i_eo, i_chunk) ! first call with ref ant
            Do i_ant=1,Ant_nr(i_chunk)
               ! N.B. It loops over all antenna's not knowing the polarity of the peak
               If(i_ant.eq.ReferenceAnt) cycle
                If(Ant_Stations(i_ant,i_chunk) .ne. Station) then
                    Dist=sqrt(sum((Ant_pos(:,i_ant,i_chunk)-Ant_pos(:,1,i_chunk))**2))
                    Station=Ant_Stations(i_ant,i_chunk)
                endif
                ! write(2,*) 'i_ant, Station, Dist',i_ant, Station, Dist
                ! flush(unit=2)
                !
                If(Ant_Stations(i_ant,i_chunk) .gt. StatMax) cycle  ! Limit to antennas near the superterp
                If(Dist .gt. DistMax*1000.) cycle  ! Limit to antennas near the superterp
                if(mod(Ant_IDs(i_ant,i_chunk),2) .ne. i_eo) cycle       ! limit to odd antennas
                If(Polariz) Then
                  if((Ant_IDs(i_ant+1,i_chunk)-Ant_IDs(i_ant,i_chunk)) .ne. 1) cycle       ! limit to even-odd antenna pairss
                EndIf
                ! write(2,*) 'GetCorrSingAnt:',i_ant, J_Corr, i_eo, i_chunk
                ! flush(unit=2)
                !
                Call GetCorrSingAnt( i_ant, J_Corr, i_eo, i_chunk) ! will infact be for a coupl when "Polariz=.true."
                !
            Enddo
            ! write(2,*) 'Nr_Corr(i_eo,i_chunk):',i_eo,i_chunk,j_corr
            ! flush(unit=2)
            Nr_Corr(i_eo,i_chunk)=j_corr
            If(PlotCCPhase) then
               Close(UNIT=30)
            endif
        Enddo
    EndDo  !  i_chunk=1, ChunkNr_dim
    !
    call DAssignFFT()                   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !stop 'BuildCC'
    !
    If(PlotCCPhase) then
      Write(31,*) '!  Average station phases'
      Do station=1,Nr_UniqueStat
            Write(31,"(3F11.4,26f8.1)") Stat_pos(:,station)/1000., &
         (CCPhase_ave(i_Peak,station),i_Peak=1,TotPeakNr(1,ChunkNr_dim))
      Enddo
      Close(UNIT=31)
    endif
   !
   ! ending time recording
   !call cpu_time(CPUstopTime)
   !TotCPUIn=TotCPUIn + CPUstopTime-CPUstartTime
   !WRITE(2,"(A,2F13.3,A)") 'Total CPU, (in/out)-side "BuildCC":',TotCPUIn, TotCPUOut, '[s]'
   !
   Return
    !stop 'BuildCC'
End Subroutine BuildCC
!=====================================
Subroutine GetCorrSingAnt( i_ant, J_Corr, i_eo, i_chunk)
!   Get the Correlation for a single antenna and all peaks in this chunk
    use DataConstants, only : Station_nrMax,PeakNr_dim,Time_dim, Production
    use DataConstants, only : Polariz, Ant_nrMax
    use Chunk_AntInfo, only : CTime_spectr, Ant_Pos, Ant_RawSourceDist, Ant_Stations, Ant_IDs
    use Chunk_AntInfo, only : Unique_SAI, Unique_StatID, Nr_UniqueStat, Tot_UniqueAnt
    use Chunk_AntInfo, only : Ant_nr
    use ThisSource, only : Peak_IRef, Peak_RefSAI, CorrAntNrs, SourcePos, Stat_pos
    use ThisSource, only : Peak_Offst,  T_Offset, CnuRt, PolBasis
    use ThisSource, only : CCPhase_ave, CCorr, CCorr_Err, CCorr_max, CCNorm, PlotCCPhase !, RealCorrelation
    use ThisSource, only : Safety, T2_Dim, Tref_dim, prntCCNrm, TotPeakNr, PeakPos, RefAntErr
    use ThisSource, only : CC_Wid, CC_Max, CC_Int, CC_WidRef, CC_Qual, CCShapeCut
    use FitParams, only : MeanCircPol, i_SAI
    use FitParams, only : N_FitPar_max, Fit_AntOffset, Fit_TimeOffsetAnt, Fit_TimeOffsetStat, PulsPosCore
    use constants, only : dp,pi,ci,sample
    use StationMnemonics, only : Station_ID2Mnem, Statn_ID2Mnem
    Use Interferom_Pars, only : dtAnt_tp, NSrc_tp
    Implicit none
    Integer, intent(inout) :: J_Corr
    Integer, intent(in) ::  i_ant, i_eo, i_chunk
    !
    !Real(dp) :: RTTraceRef(Tref_dim), RTTraceRef1(Tref_dim), RTTrace2(T2_dim), RTTrace21(T2_dim)
    Real(dp) :: RDist, SubSample_Offset, StatFineOff
    Real(dp) :: ACCorr(-Safety:Safety), CCPhase(1:PeakNr_dim)
    Real(dp) :: alpha(1:3)= (/ 1.d0, 1.d0, 1.d0/),   W=1.d0, Thet_d, Phi_d, Dist
    Real(dp), save :: Error, Errorp, CC_val, CC_Phase, RtMax, RtMaxp, SearchRange
    Integer ::  i_Peak, Sample_Offset, i
    Complex(dp) ::  CrCor(T2_dim), CrCorp(T2_dim)
    integer :: i_loc(1), StLoc, t_Max,i_stat, Antenna_SAI  !i_SAI,
    Integer, save :: PeakNr1, i_stat_old=0, N_ant
    character(len=5) :: Station_Mnem
    Real(dp), external :: CCorr_der
    logical :: prnt= .false.
    !
   !
   ! beginning time recording Ant_nr(i_chunk)
   !call cpu_time(CPUstartTime)
   !If(CPUstopTime.lt.0) CPUstopTime=CPUstartTime
   !TotCPUOut=TotCPUOut - CPUstopTime + CPUstartTime
   !If(i_ant.ge.(Ant_nr(i_chunk)-2)) then
   !   WRITE(2,"(A,F13.3,A,A,3F13.3)")  &
   !      'Total CPU, (in/out)-side "GetCorrSingAnt":',TotCPUIn, '[s]',', inside CrossCorr:',TotCC,TotSW, TotCCM
   !endif
   !
   !
    !prnt = .true.
    !prntCCNrm = .true.
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
    Antenna_SAI=1000*Ant_Stations(i_ant,i_chunk) + Ant_IDs(i_ant,i_chunk)
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
    If(PlotCCPhase .and. (j_corr.eq.1)) then
      Write(30,"(A,2I4)") '!  ',PeakNr1,TotPeakNr(i_eo,i_chunk)
    Endif
    !
    !write(2,*) 'in GetCorrSingAnt, PeakNr1=',PeakNr1,i_eo,i_chunk,j_corr, i_ant
    !flush(unit=2)
    !Production=.false.
    Do i_Peak=PeakNr1,TotPeakNr(i_eo,i_chunk) ! This assumes that the loop over antenna number goes fastest (constant [i_eo,i_chunk])
        Call RelDist(SourcePos(1,i_Peak),Ant_pos(1,i_ant,i_chunk),RDist)
        Rdist=Rdist - Ant_RawSourceDist(i_ant,i_chunk) + StatFineOff ! - INT(Peak_Offst(i_Peak))
        If(j_corr .eq. 1) then ! all is calculated w.r.t. first antenna
            If(.not. Production) then
               If(PlotCCPhase ) then
                  Write(30,"(A,i3,3F11.4)") '!  ',i_Peak,SourcePos(1:3,i_Peak)/1000.
               Else
                  Write(2,"(i3,2i2,A,I8,A,3F11.1,A,F6.2,A)") i_Peak,i_eo,i_chunk,', PeakPos=',PeakPos(i_Peak), &
                       ', source-pos=(',SourcePos(1:3,i_Peak),')[m] , RefAntErr=', RefAntErr(i_Peak),'[samples]'
               Endif
            Endif
            Peak_Offst(i_Peak)=Rdist + RefAntErr(i_Peak)  !=RDist1  !  In units of samples
            Peak_RefSAI(i_Peak)=i_SAI
            Peak_IRef(i_Peak)=i_ant
            StLoc=PeakPos(i_Peak) - Tref_dim/2
            If(PulsPosCore) StLoc=StLoc+ NINT(Peak_Offst(i_Peak)) ! Peak_Offset corrects for ref-antanna not at the center of CS002;
            !  Correction should be negligible for ref antennas in CS002 since the position shift is crudely corrected for by RawSourceDist
            !write(2,*) 'StLoc=',StLoc
            !flush(unit=2)
            !RTTraceRef(1:Tref_dim) = Real(CTime_spectr(StLoc+1:StLoc+Tref_dim,i_ant,i_chunk))
            If(Polariz) Then
               !RTTraceRef1(1:Tref_dim) = Real(CTime_spectr(StLoc+1:StLoc+Tref_dim,i_ant+1,i_chunk))
               Call GetAntSourceAng(i_ant, i_chunk, SourcePos(1,i_Peak), Thet_d, Phi_d, Dist)
               PolBasis(:,2,i_Peak)= (/sin(Phi_d*pi/180.) , -cos(Phi_d*pi/180.) , 0.d0 /)
               PolBasis(:,1,i_Peak)= (/-cos(Thet_d*pi/180.)*PolBasis(2,2,i_Peak) ,cos(Thet_d*pi/180.)*PolBasis(1,2,i_Peak) ,&
                        sin(Thet_d*pi/180.) /)
               PolBasis(:,3,i_Peak)= (/-sin(Thet_d*pi/180.)*PolBasis(2,2,i_Peak) ,sin(Thet_d*pi/180.)*PolBasis(1,2,i_Peak) ,&
                        -cos(Thet_d*pi/180.) /)
               Call GetAntPol(Thet_d, Phi_d, i_ant, i_chunk, StLoc, Tref_dim, &
                     T2_dim, PolBasis(1,1,i_Peak), alpha, W, CnuRt(0,1,i_peak)) ! get FFT spectrum around peak for (t,p) directions
            Else
               Call GetAnt(i_ant, i_chunk, StLoc, Tref_dim, T2_dim, CnuRt(0,1,i_peak))
            EndIf
            If(.not. Production) then
               CC_val=abs(CTime_spectr(StLoc+Tref_dim/2,i_ant,i_chunk))/2.
               CC_Phase=abs(CTime_spectr(StLoc+1,i_ant,i_chunk))
               RtMax=abs(CTime_spectr(StLoc+Tref_dim,i_ant,i_chunk))
               If(CC_Phase.gt.CC_val .or. CC_val.gt.CC_val) then
                  write(2,*) 'Peak',i_Peak,' NotWellCentered:(begin, center, end)=', CC_Phase, CC_val*2., RtMax
                  CC_val=abs(CTime_spectr(StLoc+Tref_dim/2+1,i_ant,i_chunk))/2.
                  CC_Phase=abs(CTime_spectr(StLoc+2,i_ant,i_chunk))
                  RtMax=abs(CTime_spectr(StLoc+Tref_dim+1,i_ant,i_chunk))
                  write(2,*) 'Moved Up:(begin+1, center+1, end+1)=', CC_Phase, CC_val*2., RtMax
              endif
           endif
        EndIf  ! (j_corr .eq. 1)
        !write(2,*) 'start corr',i_ant, J_Corr, i_eo, i_chunk
        !    flush(unit=2)
         !
        RDist=RDist - Peak_Offst(i_Peak)
        T_Offset(j_corr,i_Peak)=RDist  !  All station delays are presumably already included when extracting the time-trace
        ! where 'T_Offset' is the shift between the this and the reference antenna-timings which is
        !   already accounted for when calculating the cross correlations,
        !   which, of course, is based on the source positions before the coming fitting round.
        !  T_Offset  is used in fitting routine to determine time-shifts due to source-position shifts (and other)
        Sample_Offset = INT(RDist) ! in units of sample size
        SubSample_Offset = RDist - Sample_Offset ! in units of sample size
        !
        !write(2,*) 'Sample_Offset',J_corr, i_peak, Sample_Offset, SubSample_Offset
        StLoc=PeakPos(i_Peak) + Sample_Offset - T2_dim/2  ! Check - sign 'Sample_Offset'
        If(PulsPosCore) StLoc=StLoc+ NINT(Peak_Offst(i_Peak)) ! Peak_Offset corrects for ref-antanna not at the center of CS002;
        RtMax=0.
        Error=2000
        CC_val =0.
        CC_Phase =0.
        !CCorr(:,j_corr,i_Peak)=0.
        If(StLoc.lt.100) then
            If(.not. Production) write(2,*) '****StLoc=',StLoc,' reaches below 100, ', &
               Statn_ID2Mnem(Unique_StatID(i_stat)),i_Peak
            CCorr(:,j_corr,i_Peak)=0
            goto 9
        endif
        If(StLoc.gt.(Time_dim-T2_dim-100) ) then
            If(.not. Production) write(2,*) '****StLoc=',StLoc, ' reaches above safe upper-limit by', &
               StLoc-(Time_dim-T2_dim-100),'; ', Statn_ID2Mnem(Unique_StatID(i_stat)),i_Peak
            CCorr(:,j_corr,i_Peak)=0
            goto 9
        endif
        !RTTrace2(1:T2_dim)= Real(CTime_spectr(StLoc + 1:StLoc + T2_dim,i_ant,i_chunk))
         If(Polariz) Then
            !RTTrace21(1:T2_dim)= Real(CTime_spectr(StLoc + 1:StLoc + T2_dim,i_ant+1,i_chunk))
            Call CrossCorrPol(CnuRt(0,1,i_peak), i_ant, i_chunk, StLoc, T2_dim, &
               i_Peak, SubSample_Offset, CrCor, CrCorp)
         Else
            Call CrossCorr(CnuRt(0,1,i_peak), i_ant, i_chunk, StLoc, T2_dim, SubSample_Offset, CrCor)
         EndIf
        !
        Call SearchWin(Peak_IRef(i_Peak), i_ant, i_chunk, SourcePos(1,i_Peak), SearchRange)
        ! write(2,*) 'Before  CrossCorr_Max', SearchRange, Polariz, Ant_nrMax, j_corr,Safety
        !    flush(unit=2)
   !call cpu_time(CPUCC)
   !TotSW=TotSW+CPUCC
   !TotCCM=TotCCM-CPUCC
        If(j_corr .eq. 1) SearchRange=Safety !  i_Peak, CrCor, ACCorr, RCCorr, RtMax, Aval, Rval, Ival, SearchRange, Error)
        If(Polariz) Then
            Call CrossCorr_Max(i_ant,i_chunk, i_Peak, CrCorp, ACCorr, RtMaxp, CC_val, CC_Phase, SearchRange, Errorp)
            !write(2,*) 'j_corr+Ant_nrMax/2=', j_corr+Ant_nrMax/2, i_Peak, CC_val, Ant_nrMax, j_corr
            !flush(unit=2)
            CCorr(:,j_corr+Ant_nrMax/2,i_Peak)=ACCorr(:)/CC_val  !  /maxval(CCorr(:,j_corr,i_Peak)); used for plotting
            !write(2,*) 'inbetween  CrossCorr_Max', SearchRange, Polariz
            !flush(unit=2)
            Call CrossCorr_Max(i_ant,i_chunk, i_Peak, CrCor, ACCorr, RtMax, CC_val, CC_Phase, SearchRange, Error)
            CCorr(:,j_corr,i_Peak)=ACCorr(:)/CC_val  !  /maxval(CCorr(:,j_corr,i_Peak)); used for plotting
            !If(Abs(RtMaxp-RtMax) .gt. .1) then
            !   write(2,*) 'timing difference (t - p) polarizations', i_Peak, Antenna_SAI, RtMax-RtMaxp, 'samples'&
            !         ,Errorp-Error
            !EndIf
            If(Abs(RtMaxp-RtMax) .lt. 1.) then
               dtAnt_tp(i_SAI/2)=dtAnt_tp(i_SAI/2)+ RtMax-RtMaxp  ! difference in timing for theta and phi pol; Obsolete
               NSrc_tp(i_SAI/2)=NSrc_tp(i_SAI/2)+  1
            endif
        Else
            Call CrossCorr_Max(i_ant,i_chunk, i_Peak, CrCor, ACCorr, RtMax, CC_val, CC_Phase, SearchRange, Error)
            CCorr(:,j_corr,i_Peak)=ACCorr(:)/CC_val  !  /maxval(CCorr(:,j_corr,i_Peak)); used for plotting
        EndIf
   !call cpu_time(CPUCC)
   !TotCCM=TotCCM+CPUCC
        !If(searchrange.gt.50) write(2,*) j_corr,'SearchRange:',SearchRange, RtMax
        !If(j_corr .eq. 1) &
        !If(i_ant .le. 5 .and. i_chunk.le.1) &
        !    write(2,"(2I3,2F8.2,A,11(i2,':(',2G12.3,'), '))") i_ant, i_Peak, RtMax, RDist,' crcor:', &
        !       ( i-1,crcor(T2_dim+i),i=-2,0,1), (i-1,crcor(i), i=1,4)
        !       stop
        If(j_corr .eq. 1) then
            CC_WidRef(i_Peak)=CC_Wid  ! set in call to CC_Max
            !write(2,"(2I3,2F8.2,A,11(i2,':(',G12.3,'), '))") i_ant, i_Peak, RtMax, RDist,' ABScrcor:', &
            !   ( i-1,ABS(crcor(T2_dim+i)),i=-2,0,1), (i-1,ABS(crcor(i)), i=1,4)
        Endif
        CC_Qual(j_corr,i_Peak)=CC_Wid/CC_WidRef(i_Peak)
        If((CC_Qual(j_corr,i_Peak).lt. CCShapeCut) .or. (CC_Qual(j_corr,i_Peak).gt. 1./CCShapeCut)) then
            Error=Error+30
            If(.not. Production) write(2,*) 'Quality:',Error, CC_Qual(j_corr,i_Peak), i_ant, &
               Statn_ID2Mnem(Ant_Stations(i_ant,i_chunk)), Ant_IDs(i_ant,i_chunk)
        Endif
        !
        !If(j_corr.eq.2 .and. abs(RtMax).gt.5.) write(2,*) 'GetCorrSingAnt; CCorr(:,j_corr,i_Peak)'
        !If(j_corr.eq.2 .and. abs(RtMax).gt.5.) write(2,"(10F7.2)")  CCorr(:,j_corr,i_Peak)
    9    continue
        CCorr_max(j_corr,i_Peak) = RtMax
        CCorr_Err(j_corr,i_Peak) = Error
        CCNorm(j_corr,i_Peak)  =  CC_val
        CCPhase(i_Peak)  = CC_Phase
        If( (i_Peak.eq.1) .and. prnt ) then
            Call Station_ID2Mnem(Ant_Stations(i_ant,i_chunk),Station_Mnem)
            write(2,"(A,I3,A,i3,2x,A5,A,i2.2,A,F8.2,A,I2,A,F7.2,A,I7,f12.5)") 'j_corr=',J_corr,', i_ant=',i_ant, &
                Station_Mnem,'-',Ant_IDs(i_ant,i_chunk), ', offset=',RDist,', i_eo=',i_eo, &
                ', RtMax=',RtMax,', startLoc=',StLoc + T2_dim/2,Ant_RawSourceDist(i_ant,i_chunk)
            !write(2,*) Sample_Offset,StLoc+T2_dim/2
        endif
        !write(2,"(10g10.3,' ; ',g10.3,' ; ',10g10.3,' --')") RTTrace2(T2_dim/2-10:T2_dim/2+10) ! Self correlation, correct for upsampling
        !write(2,"(10g10.3,' ; ',g10.3,' ; ',10g10.3,' --')") &
        !    abs(CTime_spectr(StLoc + T2_dim/2-10:StLoc + T2_dim/2+10,i_ant)) ! Self correlation, correct for upsampling
        If(prnt) write(2,"(i2,10f7.3,' ; ',f7.3,' ; ',10f7.3,' --')") i_Peak,CCorr(-10:10,j_corr,i_Peak) ! Self correlation, correct for upsampling
        ! If(Error.gt. 1. ) Write(2,*) 'GetCorrSingAnt:',j_corr,Antenna_SAI,i_Peak,Error
    enddo  ! ,i_Peak
    !
    If(prntCCNrm) write(2,*) 'CCNorm ',Ant_Stations(i_ant,i_chunk), &
      (CCNorm(j_corr,i_Peak)/CCNorm(1,i_Peak),i_Peak=PeakNr1,TotPeakNr(i_eo,i_chunk))
    If(PlotCCPhase ) then
      Write(30,"(3F11.4,6f8.1)") Ant_pos(:,i_ant,i_chunk)/1000., &
         (CCPhase(i_Peak),i_Peak=PeakNr1,TotPeakNr(i_eo,i_chunk))
      If(i_stat_old .eq. i_stat) then
         N_ant=N_ant+1
         Do i_Peak=PeakNr1,TotPeakNr(i_eo,i_chunk)
            CCPhase_ave(i_Peak,i_stat)=(CCPhase_ave(i_Peak,i_stat)*(N_ant-1.) + CCPhase(i_Peak))/N_ant
         Enddo
         Stat_pos(:,i_stat)= (Stat_pos(:,i_stat)*(N_ant-1.) + Ant_pos(:,i_ant,i_chunk))/N_ant
      Else
         N_ant=1
         i_stat_old=i_stat
         Do i_Peak=PeakNr1,TotPeakNr(i_eo,i_chunk)
            CCPhase_ave(i_Peak,i_stat)=CCPhase(i_Peak)
         Enddo
         Stat_pos(:,i_stat)= Ant_pos(:,i_ant,i_chunk)
      Endif
    Endif
    !
    !stop 'GetCorrSingAnt'
    !If(j_corr.eq.1) Write(2,*) 'ACCorr',ACCorr(-2:+2)
    !Write(2,*) 'GetCorrSingAnt:',j_corr,i_Peak,
   !
   ! ending time recording
   !call cpu_time(CPUstopTime)
   !TotCPUIn=TotCPUIn + CPUstopTime-CPUstartTime
   !If(i_ant.ge.(Ant_nr(i_chunk)-1)) then
   !      WRITE(2,"(A,2F13.3,A,A,F13.3)")  &
   !      'Total CPU, (in/out)-side "GetCorrSingAnt":',TotCPUIn, TotCPUOut, '[s]',', inside CrossCorr:',TotCC
   !endif
    !
    Return
End Subroutine GetCorrSingAnt
!=================================
Real(kind = 8) Function CCorr_der(t_max,j_corr,i_Peak)
    use ThisSource, only : CCorr,CCorr_pp
    implicit none
    Integer, intent(in) :: t_max,j_corr,i_Peak
    CCorr_der=CCorr(t_max+1,j_corr,i_Peak)-CCorr(t_max,j_corr,i_Peak) &
        -CCorr_pp(t_max+1,j_corr,i_Peak)/6. -CCorr_pp(t_max,j_corr,i_Peak)/3.
    !  ypval = y(right) - y(left) - ( ypp(right) / 6. + ypp(left) / 3.0 )  &
    !       + dt * ypp(left) + dt*dt * ( ypp(right) - ypp(left) ) / 2.      !   [for h=1]
End Function CCorr_der
! ==========================! ==========================
Subroutine CrossCorrPol(CnuRt, i_ant, i_chunk, StLoc, T2_dim, i_peak, t_shft, CrCort, CrCorp)
    use FFT, only : RFTransform_CF2CT
    use ThisSource, only : SourcePos, PolBasis
    use constants, only : dp,pi,ci
    Implicit none
    !integer, parameter :: HW_size=3  ! changed on Aug 30, 2020
    Integer, Intent(in) :: i_ant, i_chunk, StLoc, T2_dim, i_peak
    Real(dp), intent(in) ::  t_shft
    Complex(dp), intent(in) :: CnuRt(0:T2_dim/2,1:3)
    Complex(dp), intent(out) :: CrCort(T2_dim), CrCorp(T2_dim)
    integer :: i
    Complex(dp) :: Cnut(0:T2_dim/2,1:3)
    complex(dp), parameter :: tipi=2*ci*pi
    Real(dp) :: alpha(1:3)= (/ 1.d0, 1.d0, 1.d0/),   W=1.d0, Thet_d, Phi_d, Dist
    !

   Call GetAntSourceAng(i_ant, i_chunk, SourcePos(1,i_Peak), Thet_d, Phi_d, Dist)
    !  Subroutine GetAntPol(Thet_d, Phi_d, RTTraceE, RTTraceO, T_dim, T2_dim, PolBasis, alpha, W,CnuPol)
   Call GetAntPol(Thet_d, Phi_d, i_ant, i_chunk, StLoc,T2_dim, T2_dim, PolBasis(1,1,i_peak), alpha, W, Cnut)

    !
!    Cnu_s(0:T2_dim/2)=Cnu2_s(0:T2_dim/2) * conjg(Cnu1_s(0:T2_dim/2))
    DO i=0,T2_dim/2
        Cnut(i,1)=Cnut(i,1) * conjg(CnuRt(i,1)) * exp(-tipi*t_shft*i/T2_dim)
        Cnut(i,2)=Cnut(i,2) * conjg(CnuRt(i,2)) * exp(-tipi*t_shft*i/T2_dim)
    ENDDO
    !
    Call RFTransform_CF2CT(Cnut(0,1),CrCort )
    Call RFTransform_CF2CT(Cnut(0,2),CrCorp )
    !call DAssignFFT()                   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !
    Return
End Subroutine CrossCorrPol
! ==========================! ==========================
Subroutine CrossCorr(CnuR, i_ant, i_chunk, StLoc, T2_dim, t_shft, CrCor)
    use FFT, only : RFTransform_CF2CT
    use constants, only : dp,pi,ci
    Implicit none
    integer, parameter :: UpSamplFact=1
    integer, parameter :: HW_size=5
    !integer, parameter :: HW_size=3  ! changed on Aug 30, 2020
    Integer, Intent(in) ::  i_ant, i_chunk, StLoc, T2_dim
    Real(dp), intent(in) :: t_shft
    Complex(dp), intent(in) :: CnuR(0:T2_dim/2)
    Complex(dp), intent(out) :: CrCor(T2_dim)
    integer :: i
    Complex(dp) :: Cnu(0:T2_dim/2)
    complex(dp), parameter :: tipi=2*ci*pi
    !

   Call GetAnt(i_ant, i_chunk, StLoc, T2_dim, T2_dim, Cnu)

    !
!    Cnu_s(0:T2_dim/2)=Cnu2_s(0:T2_dim/2) * conjg(Cnu1_s(0:T2_dim/2))
    DO i=0,T2_dim/2
        Cnu(i)=Cnu(i) * conjg(CnuR(i)) * exp(-tipi*t_shft*i/T2_dim)
    ENDDO
    !
    Call RFTransform_CF2CT(Cnu,CrCor )
    !call DAssignFFT()                   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !
    Return
End Subroutine CrossCorr
! ==========================! ==========================
Subroutine GetAntPol(Thet_d, Phi_d, i_ant, i_chunk, StLoc, T_dim, T2_dim, PolBasis, alpha, W,CnuPol)
! Get incoming e-field along 3 pre-defined orthonormal basis vectors in frequency
!   input:  Thet_d, Phi_d
!   input:  PolBasis(:,i), alpha(i) ;  i=1,3        w=AntennaNorm*Weight(j_IntFer)/AntSourceD(j_IntFer)
!   input: REAL(CTime_spectr(IntfBase:IntfBase+IntfDim,i_ant,i_chunk)) for even and odd , RTTraceE, RTTraceO, , T_dim, T2_dim,  ; NuDim=T2_dim/2
!            RTTraceRef(1:Tref_dim) = Real(CTime_spectr(StLoc+1:StLoc+Tref_dim,i_ant,i_chunk))
!  T_dim : dimension this array
!  T2_dim : dimension larger array
   use FFT, only : RFTransform_CF,RFTransform_CF2CT
   use constants, only : dp,pi,ci
    use Chunk_AntInfo, only : Ant_Stations, Ant_IDs
    use FitParams, only : MeanCircPol, Cnu01, i_SAI
   use Chunk_AntInfo, only : CTime_spectr
   use Chunk_AntInfo, only : NormOdd, NormEven !Powr_eo,NAnt_eo
   use AntFunCconst, only : Freq_min, Freq_max,Ji_p0,Ji_t0,Ji_p1,Ji_t1, Gain !, J_0p,J_0t,J_1p,J_1t
   Implicit none
   Real(dp), intent(in) :: Thet_d, Phi_d
   Integer, intent(in) :: i_ant, i_chunk, StLoc, T_dim, T2_dim
   Real(dp), intent(in) :: PolBasis(1:3,1:3), alpha(1:3), W  ! PolBasis(:,i), alpha(i) ;  i=1,3
   Complex(dp), intent(out) :: CnuPol(0:T2_dim/2,1:3)
   integer, parameter :: HW_size=5
   Real(dp) :: Hi, RTime(T2_dim,0:1), dNu !, NormEven, NormOdd
   Real(dp) :: Thet_r, Phi_r, dfreq
   Real(dp) :: Vec_p(1:3), Vec_t(1:3), Aip(1:3), Ait(1:3), Pp, Pt
   integer :: i, strt, NuDim, i_freq, i_nu, inu1, inu2, inuPh
   Complex(dp) :: Cnu0(0:T2_dim/2), Cnu1(0:T2_dim/2), Sp, St !,CTst(1:T2_dim)
   ! Get antenna function parameters
   ! Call AntFieParGen()
   !PowerScale=1.
   !NormEven=sqrt(PowerScale*2.*Powr_eo(0)/(Powr_eo(0)+Powr_eo(1)))/100.  ! to undo the factor 100 (and more) that was introduced in antenna-read
   !NormOdd=sqrt(PowerScale*2.*Powr_eo(1)/(Powr_eo(0)+Powr_eo(1)))/100.  ! to undo the factor that was introduced in antenna-read
   !NormEven=sqrt(2.*Powr_eo(0)/(Powr_eo(0)+Powr_eo(1)))/1000.  ! to undo the factor 100 (and more) that was introduced in antenna-read
   !NormOdd=sqrt(2.*Powr_eo(1)/(Powr_eo(0)+Powr_eo(1)))/1000.  ! to undo the factor that was introduced in antenna-read
   Call AntFun_Inv(thet_d,Phi_d) ! sets ,Ji_p0,Ji_t0,Ji_p1,Ji_t1; Inverse Jones; gives _p & _t polarized fields when contracted with (0.1) voltages
   Thet_r = Thet_d*pi/180  ! Zenith angle
   Phi_r  = Phi_d*pi/180   ! \phi=0 = north
   Vec_p(1)=sin(Phi_r)              ; Vec_p(2)=-cos(Phi_r)        ; Vec_p(3)=0.
   Vec_t(1)=-cos(Thet_r)*Vec_p(2)  ; Vec_t(2)=cos(Thet_r)*Vec_p(1) ; Vec_t(3)=sin(Thet_r)
   !Write(2,*) 'Vec_p(:)',Vec_p(:),Vec_t(:)
   Do i=1,3  ! current moment is calculated in basis given by PolBasis
      Aip(i)=SUM(PolBasis(:,i)*Vec_p(:))/alpha(i)
      Ait(i)=SUM(PolBasis(:,i)*Vec_t(:))/alpha(i)
   Enddo
   !If(i_ant.lt.10) write(2,*) 'GetAntPol', i_ant,ait(:),aip(:)
   !
   strt=(T2_dim-T_dim)/2
   If(strt.lt.0) then
     write(2,*) '*******strt=',strt
     stop 'GetAtPol'
   endif
   !write(2,*) 'GetAtPol',strt,HW_size, T_dim, T2_dim
   !Flush(unit=2)
   RTime(1:strt+1,:)=0.  ! padding with zeros
   RTime(strt+T_dim:T2_dim,:)=0.
   Rtime(strt+HW_size:strt+T_dim-HW_size+1,0) = NormEven*Real(CTime_spectr(StLoc+HW_size:StLoc+T_dim+1-HW_size,i_ant,  i_chunk))
   Rtime(strt+HW_size:strt+T_dim-HW_size+1,1) = NormOdd* Real(CTime_spectr(StLoc+HW_size:StLoc+T_dim+1-HW_size,i_ant+1,i_chunk))
   Do i=1,HW_size-1
     Hi=sin(0.5*(i-0.5)*pi/HW_size)**2 ! Calculate Hann window
     RTime(strt+i,0) = Hi*NormEven*Real(CTime_spectr(StLoc+i,i_ant,  i_chunk))
     RTime(strt+i,1) = Hi*NormOdd *Real(CTime_spectr(StLoc+i,i_ant+1,i_chunk))
     RTime(strt+T_dim+1-i,0)= Hi*NormEven*Real(CTime_spectr(StLoc+T_dim+1-i,i_ant,  i_chunk))
     RTime(strt+T_dim+1-i,1)= Hi*NormOdd *Real(CTime_spectr(StLoc+T_dim+1-i,i_ant+1,i_chunk))
   EndDo
   RTime(strt+T_dim+1:T2_dim,:)=0.
   !
   Call RFTransform_CF(RTime(1,0),Cnu0(0))
   Call RFTransform_CF(RTime(1,1),Cnu1(0))
   !
   CnuPol(:,:)=0.
   NuDim=T2_dim/2
   dnu=100./NuDim   ! [MHz] Jones matrix is stored on 1MHz grid
   inu1=Int(Freq_min/dnu)+1
   inu2=Int(Freq_max/dnu)-1  ! set integration regime to frequencies within filter band
   Pp=0.  ; Pt=0.
   Do i_nu=inu1,inu2
      i_freq=Int(i_nu*dnu)
      dfreq=i_nu*dnu-i_freq
      Sp=(((1.-dfreq)*Ji_p0(i_freq) + dfreq*Ji_p0(i_freq+1)) *Cnu0(i_nu) + &
         ((1.-dfreq)*Ji_p1(i_freq) + dfreq*Ji_p1(i_freq+1)) *Cnu1(i_nu))*W*Gain(i_freq)
      St=(((1.-dfreq)*Ji_t0(i_freq) + dfreq*Ji_t0(i_freq+1)) *Cnu0(i_nu) + &
         ((1.-dfreq)*Ji_t1(i_freq) + dfreq*Ji_t1(i_freq+1)) *Cnu1(i_nu))*W*Gain(i_freq)
      CnuPol(i_nu,:)=( Aip(:)* Sp + Ait(:)* St )
      Cnu0(i_nu)=Sp*conjg(St)
      Pp=Pp+ Sp*conjg(Sp)
      Pt=Pt+ St*conjg(St)
      !Pp=Pp+ Cnu0(i_nu)*conjg(Cnu0(i_nu)) ! Appears not to give a reasonable result
      !Pt=Pt+ Cnu1(i_nu)*conjg(Cnu1(i_nu))
   Enddo
   !
   !For odd-even timing check:
   If(.not. MeanCircPol) Return
   If(i_SAI.le.0) Return
   !Do i_nu=inu1,inu2
   !   Cnu0(i_nu)=Cnu0(i_nu)*conjg(Cnu1(i_nu)) ! Appears not to give a reasonable result
   !Enddo
   Cnu0(0:inu1-1)=0.
   Cnu0(inu2+1:T2_dim/2)=0.
   Sp=SUM(Cnu0(inu1:inu2)) ! The value of the stokes parameters U+iV at delta_t=0
   inuPh=-1
   If(REAL(Sp).gt.0.) inuPh=+1
   dnu=inuPh/sqrt(Pp+Pt)
   !If(i_SAI.le. 1) then
   !   write(2,*) 'inuPh',i_sai,Sp, inuPh, Pp+Pt, Pp-Pt/(Pp+Pt), ABS(Sp)*ABS(Sp)/(Pp+Pt)
   !   !Call RFTransform_CF2CT(Cnu0,CTst )
   !   !Write(2,*) Real(CTst(:))
   !EndIf
   Do i_nu=0,inu2
      Cnu01(i_nu,i_SAI/2)=Cnu01(i_nu,i_SAI/2)+Cnu0(i_nu)*dnu
      dnu=-dnu  ! to move delta_t=0 to middle of timerange
   Enddo
   !
   Return
End Subroutine GetAntPol
! ==========================! ==========================
Subroutine GetAnt(i_ant, i_chunk, StLoc, T_dim, T2_dim, Cnu)
!  T_dim : dimension this array
!  T2_dim : dimension larger array
   use FFT, only : RFTransform_CF
   use Chunk_AntInfo, only : CTime_spectr
   use constants, only : dp,pi,ci
   Implicit none
   Integer, intent(in) :: i_ant, i_chunk, StLoc, T_dim, T2_dim
   ! Real(dp), intent(in) :: RTTrace(T_dim)
    Complex(dp), intent(out) :: Cnu(0:T2_dim/2)
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
   Rtime(strt+HW_size:strt+T_dim-HW_size+1) = Real(CTime_spectr(StLoc+HW_size:StLoc+T_dim+1-HW_size,i_ant,  i_chunk))
   !write(2,*) 'Rtime(strt+HW_size:strt+T_dim-HW_size+1)',strt,Rtime(strt+HW_size:strt+T_dim-HW_size+1)
   RTime(strt+T_dim:T2_dim)=0.
   Do i=1,HW_size-1
     Hi=sin(0.5*(i-0.5)*pi/HW_size)**2 ! Calculate Hann window
     RTime(strt+i) = Hi*Real(CTime_spectr(StLoc+i,i_ant,  i_chunk))
     RTime(strt+T_dim+1-i)= Hi*Real(CTime_spectr(StLoc+T_dim+1-i,i_ant,  i_chunk))
   EndDo
   RTime(strt+T_dim+1:T2_dim)=0.
   !
   Call RFTransform_CF(RTime(1),Cnu(0))
   Return
End Subroutine GetAnt
! ==========================! ==========================
Subroutine GetAntSourceAng(i_ant, i_chunk, SourceLoc, Thet_d, Phi_d, Dist)
   use constants, only : dp, pi
   use Chunk_AntInfo, only : Ant_pos
   Implicit none
   Integer, intent(in) :: i_ant, i_chunk
   Real(dp), intent(in) :: SourceLoc(1:3)
   Real(dp), intent(out) :: Thet_d, Phi_d, Dist
   Real(dp) :: Ras(1:3), HorDist
   !
      Ras(1)=(SourceLoc(1)-Ant_pos(1,i_ant,i_chunk))/1000. !
      Ras(2)=(SourceLoc(2)-Ant_pos(2,i_ant,i_chunk))/1000.
      Ras(3)=(SourceLoc(3)-Ant_pos(3,i_ant,i_chunk))/1000.
      HorDist= Ras(1)*Ras(1) + Ras(2)*Ras(2)  ! =HYPOT(X,Y)
      Dist=sqrt(HorDist + Ras(3)*Ras(3))
      HorDist=sqrt( HorDist ) ! =HYPOT(X,Y)
      Thet_d=atan(HorDist,Ras(3))*180./pi  ! Zenith angle
      Phi_d=atan2( Ras(2) , Ras(1) )*180./pi  ! \phi=0 = north
   Return
End Subroutine GetAntSourceAng
! ==========================! ==========================
!=================================
Subroutine CrossCorr_Max(i_ant,i_chunk, i_Peak, CrCor, TrueCC, RtMax, CCval, CCPhase, SearchRange, Error)
!   Get the value of real and imaginary parts of the X-Correlation at the position of the Max in abs. value
!   Use a dynamic window to search for pulses
   !use Chunk_AntInfo
   use ThisSource, only : T2_dim, Safety, t_CCorr, ExclStatNr
   use ThisSource, only : CC_Wid, CC_Max, CC_Int, RealCorrelation
   use ThisSource, only : PlotCCPhase
   use Chunk_AntInfo, only : Ant_IDs, Ant_Stations
   use DataConstants, only : Production
   use FitParams, only : SearchRangeFallOff
   use constants, only : dp ,pi !,ci,sample
   use StationMnemonics, only : Station_ID2Mnem, Statn_ID2Mnem
   Implicit none
   !Integer, intent(in) ::  i_ant, i_chunk, j_corr,i_Peak
   Complex(dp), intent(in) ::  CrCor(1:T2_dim)
   Integer, intent(in) :: i_ant,i_chunk,i_Peak
   Real(dp), intent(in) :: SearchRange
   Real(dp), intent(out) :: TrueCC(-Safety:Safety)
   Real(dp), intent(out) :: Error, RtMax, CCval, CCPhase
   !
   integer :: i_loc(1), t_Max, i,k, Range, j
   Real(dp) :: TrueCC_pp(-Safety:Safety), RCCorr_pp(-Safety:Safety), ICCorr_pp(-Safety:Safety)
   Real(dp) :: RCCorr(-Safety:Safety), ICCorr(-Safety:Safety), Max_Time(1:2*Safety), Max_Val(1:2*Safety)
   !Real(dp) :: SearchRangeFallOff=5.D0 ! as multiple of SearchRange
   Real(dp), save :: tA,B,Yp, dt, Rval, Ival
   !
   !write(2,*) 'Entering CrossCorr_Max', i_ant,i_chunk, i_Peak, SearchRange !, ', CrCor=', CrCor
   !flush(unit=2)
   If(RealCorrelation) then
      TrueCC(0:Safety) = Real(CrCor(1:Safety+1))
      TrueCC(-Safety:-1) = Real(CrCor(T2_dim-Safety+1:T2_dim))  ! correct for upsampling
   else
      TrueCC(0:Safety) = abs(CrCor(1:Safety+1))
      TrueCC(-Safety:-1) = abs(CrCor(T2_dim-Safety+1:T2_dim))  ! correct for upsampling
   endif
   !
    i=Safety/3
    CC_Int = sqrt(SUM( TrueCC(-i:i)*TrueCC(-i:i) ))
    CC_Max = MAXVAL( TrueCC(-i:i) )
    CC_Wid = CC_Max / CC_Int
   !
   Do i=1,Safety ! (flat part upto error, lin decreasing beyond) is changed to parabola
      B=(1. - i*i/(SearchRange*SearchRange*SearchRangeFallOff*SearchRangeFallOff))
      TrueCC(i) = TrueCC(i)*B
      TrueCC(-i) = TrueCC(-i)*B
   Enddo
   Range=NINT(SearchRange*SearchRangeFallOff)
   if(Range.gt.Safety) Range=Safety
   !if(i_ant.eq.1) write(2,*) 'TrueCC', TrueCC(-Range:Range)
   !if(i_ant.eq.1) write(2,*) 'CrCor', crcor(T2_dim), crcor(1), crcor(2)
   CCPhase = 0.
   Error=0.
   If(count((ExclStatNr(:,i_peak)-Ant_Stations(i_ant,i_chunk)).eq.0,1) .ge. 1) then
      !write(2,*) 'excluded station=',Ant_Stations(i_ant,i_chunk),i_ant, ', for i_peak=',i_peak
      CCval=maxval(TrueCC(-Range:Range))
      RtMax=0.
      Error=2000
      Return
   EndIf
   !
!   Call spline_cubic_set( 2*Safety+1, t_ccorr(-Safety), TrueCC(-Safety), TrueCC_pp(-Safety) )
   Call spline_cubic_set( 2*Range+1, t_ccorr(-Range), TrueCC(-Range), TrueCC_pp(-Range) )
   !
   !
   If(RealCorrelation) then ! first find all local maxima, then get the global one
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
         CCval=B*B - 2.*tA * Yp
         If(CCval.le.0.) cycle
         CCval=sqrt(CCval)
         dt= ( -B - CCval )/tA
         !If(abs(tA).lt.1.e-20) dt=-yp/B
         !write(2,*) i,'yp:',yp,( -B - sqrt(B*B - 2.*tA * Yp) )/tA,( -B + sqrt(B*B - 2.*tA * Yp) )/tA
         If(dt.lt.0. .or. dt.gt.1) then
            dt= ( -B + CCval )/tA
            !If(dt.lt.0. .or. dt.gt.1) then
            If(dt.gt.0. .and. dt.lt.1) then
               write(2,*) 'dt',dt,yp,B,tA,CCval,'found THE exception'
               cycle
            Endif
         Endif
         CCval= TrueCC(i) &
             + dt*( Yp + dt*( 0.5D0*TrueCC_pp(i) + dt*( (TrueCC_pp(i+1)-TrueCC_pp(i))/6.D0 ) ) )
         !If(CCval.lt.TrueCC(i+1)) write(2,*) 'CCval?;',CCval, TrueCC(i+1)
         If(CCval.lt.0.) cycle
         j=j+1
         Max_Time(j)=i+dt
         Max_Val(j) = CCval
      Enddo
      j=j+1
      Max_Time(j)= Range
      Max_Val(j) = TrueCC(Range)
      i_loc=MaxLoc(Max_Val(1:j) )
      RtMax=Max_Time(i_loc(1))
      t_max=NINT(RtMax)
      CCval=Max_Val(i_loc(1))
      i=COUNT(Max_Val(1:j).gt. 0.95*CCval)
      If(i.gt.1) then
         !If(.not. Production) write(2,"(5('*'),A,i3,A,i2,A,i4,' ',A)") 'For peak#', i_Peak,' there are ',i, &
         !   ' large correlation maxima found in antenna',Ant_IDs(i_ant,i_chunk),Statn_ID2Mnem(Ant_Stations(i_ant,i_chunk))
         Error=(i-1)*5.    !  Peak becomes ambguous as there are more equivalent ones
      EndIf
      !If(i_ant.le.20) then
      !   write(2,"(5x,A,3i3,A,F7.3,g11.3)") 'i_Ant & Peak:',i_ant, i_Peak,j,' max:',RtMax,CCval
      !   write(2,"(A,20F9.3)") 'max_time:',Max_Time(1:j)
      !   write(2,"(A,20F9.1)") 'max_val: ',Max_Val(1:j)
      !EndIf
   Else
      i_loc=MaxLoc(TrueCC(-Range:Range) ) ! works only for a very smooth function, not for real correlation
      t_max=i_loc(1) - Range -1
      CCval=-1.
   EndIf
   !If(i_peak.eq.1) write(2,*) i_ant,TrueCC(:)
   !if(i_ant.eq.5) stop
   !Write(2,*) 't_max:',t_max,TrueCC(t_max),CC_Max
   !         flush(unit=2)
   !If(Unique_StatID(i_stat).eq.147) Write(2,*) '147,i,j',i_ant,j_corr,'Corrected Rdist=',Rdist,t_max
   !
   CCPhase = 0.
   If((t_max.le. -Range) .or. (t_max .ge. Range)) then
      If(.not. Production) write(2,"(A,I5,I3,I4,A,A5,A,I2)") 'maximum in correlation function at bound',t_max, &
          Ant_IDs(i_ant,i_chunk),Ant_Stations(i_ant,i_chunk),'=',Statn_ID2Mnem(Ant_Stations(i_ant,i_chunk)),', peak#=',i_Peak
      CCval=maxval(TrueCC(-Range:Range))
      RtMax=t_max
      Error=Error+200  ! time samples
   Else
      If(CCval.lt.0.) then ! it was not set before
         Yp=TrueCC(t_max+1)-TrueCC(t_max) - TrueCC_pp(t_max+1)/6. -TrueCC_pp(t_max)/3.    !CCorr_der(t_max,j_corr,i_Peak)
         if(Yp.lt.0.) then
             t_max=t_max-1  !  shift to the point left of the real maximum
             Yp=TrueCC(t_max+1)-TrueCC(t_max) - TrueCC_pp(t_max+1)/6. -TrueCC_pp(t_max)/3. !CCorr_der(t_max,j_corr,i_Peak)
         endif
         tA=TrueCC_pp(t_max+1) - TrueCC_pp(t_max)   ! 2 * A
         B=TrueCC_pp(t_max)
         dt= ( -B - sqrt(B*B - 2.*tA * Yp) )/tA
         RtMax=t_max+dt  ! position of the maximum in the correlation function
!         Call spline_cubic_val( 2*Safety+1, t_ccorr(-Safety), TrueCC(-Safety), TrueCC_pp(-Safety), RtMax, CCval)
         Call spline_cubic_val( 2*Range+1, t_ccorr(-Range), TrueCC(-Range), TrueCC_pp(-Range), RtMax, CCval)
      EndIf
      If(abs(RtMax).gt.2.*SearchRange)  Then
         Error = Error+20  ! time samples
         If(.not. Production) write(2,*) 'Large deviation:',SearchRange,RtMax, Statn_ID2Mnem(Ant_Stations(i_ant,i_chunk)),i_ant
      Else
         Error= Error+0.2  ! time samples
      EndIf
      !
      !calculate the phase of the cross correlation at the point of interest
      If( PlotCCPhase ) then
         ICCorr(0:Safety) = Imag(CrCor(1:Safety+1))
         ICCorr(-Safety:-1) = Imag(CrCor(T2_dim-Safety+1:T2_dim))  ! correct for upsampling
         Call spline_cubic_set( 2*Safety+1, t_ccorr(-Safety), ICCorr(-Safety), ICCorr_pp(-Safety) )
         Call spline_cubic_val( 2*Safety+1, t_ccorr(-Safety), ICCorr(-Safety), ICCorr_pp(-Safety), RtMax, Ival) !, ypval, yppval )
         If(RealCorrelation) then
            Rval=CCval  ! CCval is the value of the real cc at max
         else
            RCCorr(0:Safety) = Real(CrCor(1:Safety+1))
            RCCorr(-Safety:-1) = Real(CrCor(T2_dim-Safety+1:T2_dim))  ! correct for upsampling
            Call spline_cubic_set( 2*Safety+1, t_ccorr(-Safety), RCCorr(-Safety), RCCorr_pp(-Safety) )
            Call spline_cubic_val( 2*Safety+1, t_ccorr(-Safety), RCCorr(-Safety), RCCorr_pp(-Safety), RtMax, Rval) !, ypval, yppval )
         endif
         CCPhase  = atan2(Ival, Rval)*180./pi
      Endif
      !If(i_ant.lt.10 .and. i_peak.lt.3) write(2,*) i_ant, i_Peak,'cc-corr:',Rval,Ival,RtMax
      !write(2,"(A,i3,I5,A,A5,A,I2,F7.1,F11.1,f8.1)") ' ant#=',Ant_IDs(i_ant,i_chunk),Ant_Stations(i_ant,i_chunk),'=', &
      !      Statn_ID2Mnem(Ant_Stations(i_ant,i_chunk)),&
      !      ', peak#=',i_Peak, RtMax, CCval, atan2(Ival, Rval)*180./pi
      !flush(unit=2)
   endif
   !write(2,*) 'Exit  CrossCorr_Max', RtMax, CCval, CCPhase, Error, '; TrueCC=', TrueCC
   !flush(unit=2)
   !
   Return
End Subroutine CrossCorr_Max
!=====================================
!=================================
Subroutine ReImAtMax(i_ant,i_chunk, i_Peak, CrCor, ACCorr, RCCorr, RtMax, Aval, Rval, Ival, SearchRange, Error)
!   Get the value of real and imaginary parts of the X-Correlation at the position of the Max in abs. value
!   Use a dynamic window to search for pulses
   !use Chunk_AntInfo
   use ThisSource, only : T2_dim, Safety, t_CCorr, ExclStatNr
   use ThisSource, only : CC_Wid, CC_Max, CC_Int
   use Chunk_AntInfo, only : Ant_IDs, Ant_Stations
   use DataConstants, only : Production
   use FitParams, only : SearchRangeFallOff
   use constants, only : dp ,pi !,ci,sample
   use StationMnemonics, only : Station_ID2Mnem, Statn_ID2Mnem
   Implicit none
   !Integer, intent(in) ::  i_ant, i_chunk, j_corr,i_Peak
   Complex(dp), intent(in) ::  CrCor(1:T2_dim)
   Integer, intent(in) :: i_ant,i_chunk,i_Peak
   Real(dp), intent(out) :: ACCorr(-Safety:Safety), RCCorr(-Safety:Safety)
   Real(dp), intent(out) :: RtMax, Aval, Rval, Ival
   Real(dp), intent(inout) :: Error, SearchRange
   !
   integer :: i_loc(1), t_Max, i,k
   Real(dp) :: ACCorr_pp(-Safety:Safety), RCCorr_pp(-Safety:Safety), ICCorr_pp(-Safety:Safety), ICCorr(-Safety:Safety)
   !Real(dp) :: SearchRangeFallOff=5.D0 ! as multiple of SearchRange
   Real(dp), save :: tA,B,Yp, dt
   !
   ACCorr(0:Safety) = abs(CrCor(1:Safety+1))
   ACCorr(-Safety:-1) = abs(CrCor(T2_dim-Safety+1:T2_dim))  ! correct for upsampling
   ICCorr(0:Safety) = Imag(CrCor(1:Safety+1))
   ICCorr(-Safety:-1) = Imag(CrCor(T2_dim-Safety+1:T2_dim))  ! correct for upsampling
   RCCorr(0:Safety) = Real(CrCor(1:Safety+1))
   RCCorr(-Safety:-1) = Real(CrCor(T2_dim-Safety+1:T2_dim))  ! correct for upsampling
   !
    i=Safety/3
    CC_Int = SUM( ACCorr(-i:i) )
    CC_Max = MAXVAL( ACCorr(-i:i) )
    CC_Wid = CC_Max / CC_Int
   !
   Do i=1,Safety ! (flat part upto error, lin decreasing beyond) is changed to parabola
      B=(1. - i*i/(SearchRange*SearchRange*SearchRangeFallOff*SearchRangeFallOff))
      ACCorr(i) = ACCorr(i)*B
      ACCorr(-i) = ACCorr(-i)*B
   Enddo
   !
   Call spline_cubic_set( 2*Safety+1, t_ccorr(-Safety), ACCorr(-Safety), ACCorr_pp(-Safety) )
   Call spline_cubic_set( 2*Safety+1, t_ccorr(-Safety), ICCorr(-Safety), ICCorr_pp(-Safety) )
   Call spline_cubic_set( 2*Safety+1, t_ccorr(-Safety), RCCorr(-Safety), RCCorr_pp(-Safety) )
   !
   i_loc=MaxLoc(ACCorr(:) )
   t_max=i_loc(1) - Safety -1
   !Write(2,*) 't_max:',t_max
   !If(Unique_StatID(i_stat).eq.147) Write(2,*) '147,i,j',i_ant,j_corr,'Corrected Rdist=',Rdist,t_max
   If(count((ExclStatNr(:,i_peak)-Ant_Stations(i_ant,i_chunk)).eq.0,1) .ge. 1) then
      !write(2,*) 'excluded station=',Ant_Stations(i_ant,i_chunk),i_ant, ', for i_peak=',i_peak
      Aval=maxval(ACCorr(:))
      RtMax=0.
      Error=2000
   ElseIf((t_max.le. -Safety) .or. (t_max .ge. Safety)) then
      If(.not. Production) write(2,"(A,I5,I3,I4,A,A5,A,I2)") 'maximum in correlation function at bound',t_max, &
          Ant_IDs(i_ant,i_chunk),Ant_Stations(i_ant,i_chunk),'=',Statn_ID2Mnem(Ant_Stations(i_ant,i_chunk)),', peak#=',i_Peak
      Aval=maxval(ACCorr(:))
      RtMax=t_max
      Error=200  ! time samples
   Else
      Yp=ACCorr(t_max+1)-ACCorr(t_max) - ACCorr_pp(t_max+1)/6. -ACCorr_pp(t_max)/3.    !CCorr_der(t_max,j_corr,i_Peak)
      if(Yp.lt.0.) then
          t_max=t_max-1  !  shift to the point left of the real maximum
          Yp=ACCorr(t_max+1)-ACCorr(t_max) - ACCorr_pp(t_max+1)/6. -ACCorr_pp(t_max)/3. !CCorr_der(t_max,j_corr,i_Peak)
      endif
      tA=ACCorr_pp(t_max+1) - ACCorr_pp(t_max)   ! 2 * A
      B=ACCorr_pp(t_max)
      dt= ( -B - sqrt(B*B - 2.*tA * Yp) )/tA
      RtMax=t_max+dt  ! position of the maximum in the correlation function
      !write(2,*) 'ReImAtMax:',t_max,dt,B,tA,Yp
      If(abs(RtMax).gt.2.*SearchRange)  Then
         Error = 20  ! time samples
         If(.not. Production) write(2,*) 'Large deviation:',SearchRange,RtMax, Statn_ID2Mnem(Ant_Stations(i_ant,i_chunk)),i_ant
      Else
         Error= 0.2  ! time samples
      EndIf
      !If(abs(RtMax).gt. Safety/4.)  Error=20.  ! October 2019; error adjustment
      Call spline_cubic_val( 2*Safety+1, t_ccorr(-Safety), ACCorr(-Safety), ACCorr_pp(-Safety), RtMax, Aval)
      Call spline_cubic_val( 2*Safety+1, t_ccorr(-Safety), ICCorr(-Safety), ICCorr_pp(-Safety), RtMax, Ival) !, ypval, yppval )
      Call spline_cubic_val( 2*Safety+1, t_ccorr(-Safety), RCCorr(-Safety), RCCorr_pp(-Safety), RtMax, Rval) !, ypval, yppval )
      !write(2,"(A,i3,I5,A,A5,A,I2,F7.1,F11.1,f8.1)") ' ant#=',Ant_IDs(i_ant,i_chunk),Ant_Stations(i_ant,i_chunk),'=', &
      !      Statn_ID2Mnem(Ant_Stations(i_ant,i_chunk)),&
      !      ', peak#=',i_Peak, RtMax, Aval, atan2(Ival, Rval)*180./pi
   endif
   !
   Return
End Subroutine ReImAtMax
!=====================================
Subroutine GetRefAnt
   ! with the dual option the relative distance between the two polarity antennas should be small
   use DataConstants, only : ChunkNr_dim
   use DataConstants, only : Polariz
   use Chunk_AntInfo, only : Ant_IDs, Ant_pos, RefAnt
   use ThisSource, only : Dual
   use constants, only : dp
   use unque, only : Double_sort
   Implicit none
   Integer, parameter :: MaxRefAnt=20, i_e=0, i_o=1, MA=100  ! max nr of antenna pairs tested for close proximity
   Integer :: Ant(1:Ma,0:2), j, i_ante, i_anto, i_chunk
   Real(dp) :: Dist
   Do i_chunk=1, ChunkNr_dim
      If(Dual) then
         j=0
         Do i_ante=1,MaxRefAnt ! build table of relative distances
            if(mod(Ant_IDs(i_ante,i_chunk),2) .ne. i_e) cycle       ! check this is an even antenna
            If(Polariz) Then  ! the reference should be at the same position
               If((Ant_IDs(i_ante,i_chunk)+1) .eq. Ant_IDs(i_ante+1,i_chunk)) Then ! Reference antenna found
                  RefAnt(i_chunk,0)=i_ante
                  RefAnt(i_chunk,1)=i_ante+1
                  Goto 9
               Else
                  cycle
               EndIf
            EndIf
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
         If(j.eq.0) stop 'No reference antenna found, Dual'
         Call Double_sort(Ant(1:j,0:2))
         write(2,*) 'relative reference antenna distance [0.1m] and numbers:',Ant(1,:)
         j=1
      Else
         Do i_ante=1,MaxRefAnt ! do NOT build table of relative distances
            if(mod(Ant_IDs(i_ante,i_chunk),2) .ne. i_e) cycle       ! check this is an even antenna
            Do i_anto=1,MaxRefAnt
               if(mod(Ant_IDs(i_anto,i_chunk),2) .ne. i_o) cycle       ! check this is an odd antenna
               RefAnt(i_chunk,0)=i_ante
               RefAnt(i_chunk,1)=i_anto
               goto 9
            Enddo
         Enddo
         If(j.eq.0) stop 'No reference antenna found, Single'
      Endif
   1  continue
!      write(2,*) 'GetRefAnt2',j,Ant(j,:)
!      flush(unit=2)
      RefAnt(i_chunk,0)=Ant(j,2)
      RefAnt(i_chunk,1)=Ant(j,1)
   9  Continue
      !If(Dual) write(2,*) j,'reference antenna* distancebetween',Dist,'[m], and numbers:',RefAnt(i_chunk,0),RefAnt(i_chunk,1)
   Enddo
   Return
End Subroutine GetRefAnt
!=================================
