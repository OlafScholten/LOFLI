!=================================
!================================
!=================================
Subroutine FindCallibr(SourceGuess)
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
   use Chunk_AntInfo, only : Unique_SAI, Start_time, TimeFrame, tot_uniqueant
   use DataConstants, only : Polariz, Ant_nrMax
   use ThisSource, only : Dual, RealCorrelation
   use ThisSource, only : NrP, t_ccorr, SourcePos, RefAntErr, ExclStatNr, Peak_eo, ChunkNr, TotPeakNr
   use ThisSource, only : Nr_corr, PeakNr, PeakNrTotal, PlotCCPhase, Safety, PolBasis, Tref_dim, T2_Dim
   use ThisSource, only : PeakPos, Peak_eo, PeakRMS, PeakChiSQ
   use ThisSource, only : XFrameEi, XFrameEf, XFrameNi, XFrameNf, XFrameh
   use FitParams, only : MaxFitAntDistcs, MaxFitAntD_nr, FitIncremental, PulsPosCore, WriteCalib, N_EffAnt, Max_EffAnt
   use FitParams, only : FullSourceSearch
   use constants, only : dp,sample
   use DataConstants, only : PeakNr_dim, ChunkNr_dim, RunMode
   use FitParams, only : MeanCircPol, Cnu01, Fit_TimeOffsetAnt, Fit_AntOffset, N_FitStatTim, FitParam, X_Offset
   use FFT, only : RFTransform_su,DAssignFFT, RFTransform_CF2CT
   !use StationMnemonics, only : Statn_ID2Mnem, Station_Mnem2ID
   Implicit none
   Real(dp), intent(in) :: SourceGuess(3,10)
   !
   integer :: i, j, k, i_ant, j_corr, i_eo, i_chunk, i_dist, MinFitAntD_nr
   logical, save :: Fitfirst=.true.
   !Character(len=5) :: Station_Mnem, ExclStMnem(1:Station_nrMax)
   integer :: i_loc(1), i_Peak, StLoc, StatMax, ReadErr, PeakNr1
   Integer :: PeakSP(2*NrP,0:2), PeakSWl(2*NrP,0:2), PeakSWu(2*NrP,0:2), PeakSAmp(2*NrP,0:2), PeakD_nr
   Real(dp) :: DistMax
   !character*80 :: lname  ! Should be 80 to be compatible with "GetNonZeroLine"
   character*10 :: txt
   character*35 :: Frmt
   real(dp) :: x1,x2,x3,t
   logical :: FitNoSources=.true.
    real(dp) :: dt,B, MTC, dtI
    Complex(dp), Allocatable :: ST01nu(:)
    Complex(dp), Allocatable :: ST01T(:)
    Complex(dp) :: ST01_max
    Integer :: UpSF= 16, i_SAI ! UpSampleFactor
   !
   !       Initialize
   !Polariz=(Dual .and. RealCorrelation)
   Call Find_unique_StatAnt()  ! works in polariz mode
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
   ExclStatNr(:,:)=0
   Call GetRefAnt      ! works in polariz mode
   !
   i_peak=0
   PeakNr1=0
   If(RunMode.eq.1) then ! To prevent reading from input in Explore
      ReadErr=-1
   Else
      Call ReadPeakInfo(ReadErr)
   Endif
   !write(2,*) 'ReadErrh: ',ReadErr
   !flush(unit=2)
   If(ReadErr.ne.0) then  ! Find peak positions instead of reading them from input
      Do i_chunk=1, ChunkNr_dim     ! get the NrP strongest peaks in this time block
        Call DualPeakFind(2*NrP, i_chunk, PeakD_nr, PeakSP, PeakSWl, PeakSWu, PeakSAmp) ! used for imaging production
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
      enddo ! i_chunk=1, ChunkNr_dim
      PeakNrTotal=PeakNr1
   endif
   !write(2,*) 'PeakNrTotal:',PeakNrTotal
   !write(2,*) 'SourcePos:',SourcePos
   !
   If(Polariz) then
      write(2,*) 'Use theta polarization for antenna correlations'
      Allocate( PolBasis(1:3,1:3,1:PeakNrTotal) )
   Endif
   write(2,*) 'FullSourceSearch: ', FullSourceSearch
   flush(unit=2)
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
   write(2,*) 'i_dist: ', i_dist
   flush(unit=2)
      !Write(2,*) 'Average RMS=',SUM(PeakRMS(1:PeakNrTotal))/PeakNrTotal
      If((RunMode.eq.1) .and. (MINVAL(PeakRMS(1:PeakNrTotal)).gt.30.)) then
         write(2,"(A,F7.2,A,10F6.0)") 'chi-square too poor for max antenna distance=',DistMax &
            ,'[m], PeakRMS(1:PeakNrTotal):',PeakRMS(1:PeakNrTotal)
         Return
      Endif
   EndDo ! i_dist=1,MaxFitAntD_nr
   flush(unit=2)
    !
    If(RunMode.eq.1) then  ! =Explore
      !write(2,*) 'PeakNrTotal=',PeakNrTotal
      Do i_peak=1,PeakNrTotal
         If(PeakRMS(i_Peak).gt. 25.) cycle
         write(2,"(i2,i7,', x,y,z=', 3F9.0,', RMS=',F7.1)") i_peak,Peakpos(i_peak),SourcePos(:,i_Peak),PeakRMS(i_Peak)
         If((SourcePos(2,i_Peak)/1000. .lt. XFrameEi) .or. (SourcePos(2,i_Peak)/1000. .gt. XFrameEf)) cycle
         If((SourcePos(1,i_Peak)/1000. .lt. XFrameNi) .or. (SourcePos(1,i_Peak)/1000. .gt. XFrameNf)) cycle
         If(SourcePos(1,i_Peak) .ne. SourcePos(1,i_Peak)) cycle  ! check for NaN
         If(abs(SourcePos(3,i_Peak))/1000. .gt. XFrameh) cycle
         write(29,"(1x,i5,4(2x,g14.8),3x,f9.6,2I5)")  &
            TimeFrame*10+i_Peak,Start_time(1)*1000.*sample, &
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
    If(Polariz) Then
      Allocate( Cnu01(0:T2_Dim/2,0:Ant_nrMax/2) )
      Cnu01(:,:)=0.d0
      MeanCircPol=.true.
    EndIf
    PlotCCPhase=.true.
    Call BuildCC(StatMax,DistMax)
   !
   If(MeanCircPol) Then
      call RFTransform_su(UpSF*T2_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      Allocate( ST01T(1:UpSF*T2_dim) )
      Allocate( ST01nu(0:UpSF*T2_dim/2) )
      !    Write(2,*) 'Mean t-p dt:',SUM(dtAnt_tp(:))/SUM(NSrc_tp(:)),'Polariz',Polariz, T2_dim, Tref_dim
      Do i=1, N_FitStatTim
         k = FitParam(i)   ! =Station number
         If(Fit_AntOffset .and. Polariz) then ! The other cases are dealt with in X2Source
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
   Write(2,"(//,A)") ' ==== Summary of new parameters ===='
   Call PrntNewSources
   ! Merge values for "Fit_TimeOffsetStat" with input from 'FineCalibrations.dat'
   If(WriteCalib) then
      Call MergeFine
   Endif
   !
   !Call GLE_Corr()
   PlotCCPhase=.false.
   Return
    !
End Subroutine FindCallibr
!=================================
Subroutine ReadPeakInfo(ReadErr)
!   Read-in the peakpos and source guesses for individual pulses.
!
!
    use Chunk_AntInfo, only : Station_nrMax
    use ThisSource, only : SourcePos, RefAntErr, ExclStatNr, Peak_eo, TotPeakNr, ChunkNr !NrP, t_ccorr,
    use ThisSource, only : Dual, PeakNr, PeakNrTotal !, PlotCCPhase, Safety, Nr_corr
    use ThisSource, only : PeakPos, Peak_eo
!    use FitParams
    use constants, only : dp,sample
    use DataConstants, only : PeakNr_dim, ChunkNr_dim
    use DataConstants, only : Polariz
    use StationMnemonics, only : Statn_ID2Mnem, Station_Mnem2ID
    Implicit none
    !Character(len=80), intent(out) :: lname
    Integer, intent(out) :: ReadErr
    !logical, intent(in) :: WriteCalib
    !real ( kind = 8 ) :: X(N_FitPar_max)
    !
    integer :: i_eo, i_chunk, i_dist, i, k, i_c,i_ca,i_eoa
    !logical :: Fitfirst !,TraceFirst
    Character(len=5) :: Station_Mnem, ExclStMnem(1:Station_nrMax)
    !Integer, parameter :: PeakS_dim=NrP
    integer :: i_Peak, nxx, PeakNr1
    character*80 :: lname  ! Should be 80 to be compatible with "GetNonZeroLine"
    character*10 :: txt
    character*35 :: Frmt
    real(dp) :: x1,x2,x3,t
    !
   ReadErr=0
   Call GetNonZeroLine(lname)
!      Write(2,"(i3,2i2,I8)", ADVANCE='NO') i_Peak,Peak_eo(i_Peak),ChunkNr(i_Peak),PeakPos(i_Peak)
!      Write(2,"(3(F10.2,','))", ADVANCE='NO') SourcePos(:,i_Peak)
!      Write(2,"(F8.2)", ADVANCE='NO') RefAntErr(i_Peak)
   Frmt="(i3,2i2,I8,3(F10.2,1x),F8.2)"
   read(lname,Frmt,iostat=nxx)  i,i_eoa,i_ca, k, x1,x2,x3,t
   If(nxx.ne.0) then
      write(2,*) 'in:"',lname,'"'
      ReadErr=-1
      Return
   Endif
    !
    !
   Do i_Peak=1,PeakNr_dim      ! initialize
     !SourcePos(:,i_Peak) = SourceGuess(:)
     RefAntErr(i_Peak) = 0.
   Enddo
   i_peak=0
   TotPeakNr=0
   PeakNr=0
   i_chunk=1
   PeakNr1=0
   Do      ! Read source positions from input. There should be at least one un-readable line in the input.
      !write(2,*) 'lname="',lname,'"'
      read(lname,Frmt,iostat=nxx)  i,i_eo,i_c, k, x1,x2,x3,t
      !write(2,*) k, x1,x2,x3,t
      If(nxx.ne.0) exit
      if((i_eo.lt.0) .or. (i_eo.gt.1) ) exit
      i_peak=i_peak+1
      If(i_peak .gt. PeakNr_dim) then
         write(2,*) 'PeakNr_dim will be exceeded',i_peak
         stop 'FindCallibr:ReadPeakInfo'
      endif
      SourcePos(1,i_Peak)=x1
      SourcePos(2,i_Peak)=x2
      SourcePos(3,i_Peak)=x3
      RefAntErr(i_Peak) = t
      Peakpos(i_Peak)=k
      Peak_eo(i_peak)=i_eo
      If(i_eoa .ne. i_eo) then
          PeakNr1= i_peak-1
          If(i_eo.eq. 0) i_chunk=i_chunk+1
      Elseif(i_ca.ne.i_c) then ! i_eo did not jump
           TotPeakNr(1,i_chunk)=i_peak-1 ! last peak# for this (i_eo,i_chunk)
           PeakNr(1,i_chunk)=0        ! number of peaks for this (i_eo,i_chunk)
         i_chunk=i_chunk+1
      Endif
      if(i_chunk.gt.ChunkNr_dim) exit
      If(Dual .and. (i_eo.eq.1)) then
      Do i=0,PeakNr(0,i_chunk)-1
         !write(2,*) 'peakpos-i',i,TotPeakNr(0,i_chunk)-i,Peakpos(TotPeakNr(0,i_chunk)-i),k
         If(Peakpos(TotPeakNr(0,i_chunk)-i).eq. k) then
            SourcePos(:,i_Peak)=SourcePos(:,TotPeakNr(0,i_chunk)-i)
            exit
         Endif
      Enddo
     endif
     !write(2,*) 'Source ',i_Peak, ' searched near ',SourcePos(:,i_Peak)
     ChunkNr(i_peak)=i_chunk
     TotPeakNr(i_eo,i_chunk)=i_peak ! last peak# for this (i_eo,i_chunk)
     PeakNr(i_eo,i_chunk)=i_peak-PeakNr1        ! number of peaks for this (i_eo,i_chunk)
     PeakNrTotal=i_peak
     !write(2,*) i_eo,i_chunk,'PeakNr(i_eo,i_chunk)=',PeakNr(i_eo,i_chunk)
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
     If(Polariz .AND. (i_eo.eq.1) ) then !  take i_eo=0 peak info only
        i_peak=TotPeakNr(0,i_chunk)
        TotPeakNr(1,i_chunk)=i_peak ! last peak# for this (i_eo,i_chunk)
        PeakNr(1,i_chunk)=0        ! number of peaks for this (i_eo,i_chunk)
        PeakNrTotal=i_peak
        PeakNr1=i_peak
     EndIf
      i_eoa=i_eo
      i_ca=i_c
   Enddo
    write(2,*) 'PeakNrTotal=',PeakNrTotal
    write(2,*) 'TotPeakNr',TotPeakNr
    !
    !write(2,*) 'SourcePos=',SourcePos
    !
    !
    Return
End Subroutine ReadPeakInfo
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
    write(2,*) '=== New Round ============================================================'
    flush(unit=2)
    CalcHessian=.false.
    SpaceCov(1,1)=-.1 ; SpaceCov(2,2)=-.1 ; SpaceCov(3,3)=-.1 ;
    Call BuildCC(StatMax,DistMax)
    !
    FitPos(1)=1 ; FitPos(2)=2 ; FitPos(3)=3 ; FitPos(4)=4  !(4=timing error, 1-3 is N,E,h)
    If(DistMax.le. 0.30) then
        FitPos(1)=1 ; FitPos(2)=2 ; FitPos(3)=-1 ; FitPos(4)=-1 ! single station, fit direction only
    ElseIf(DistMax.le. 1.0) then
        FitPos(1)=1 ; FitPos(2)=2 ; FitPos(3)=3 ; FitPos(4)=-1 ! Superterp, fit direction and curvature front
    ElseIf(DistMax.le. 5) then
        FitPos(1)=1 ; FitPos(2)=2 ; FitPos(3)=3 ; FitPos(4)=4
    EndIf
    If(FitNoSources) then
        FitPos(1)=-1 ; FitPos(2)=-2 ; FitPos(3)=-3 ; FitPos(4)=-4
    endif
    Call SetFitParamsLMA(X,Fitfirst,FitPos) ! Important to call 'SetFitParamsLMA' after 'GetCorrSingAnt'
    write(2,*) 'N_FitPar1 , N_FitStatTim',N_FitPar , N_FitStatTim
    flush(unit=2)
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
    write(2,*) 'chi-square=',FitQual, ', Average RMS=',SUM(PeakRMS(1:PeakNrTotal))/PeakNrTotal, PeakRMS(1:PeakNrTotal)
    flush(unit=2)
    Call X2Source(X)
    !
End Subroutine FitCycle
!=====================================
Subroutine GetLargePeaks(i_ant, i_chunk, Peakposs)
    use DataConstants, only : Time_dim
    use Chunk_AntInfo, only : CTime_spectr, Ant_IDs
    use ThisSource, only : NrP  ! PeakNr
    use constants, only : dp
    Implicit none
    Integer, intent(in) :: i_ant, i_chunk
    Integer, intent(out) :: Peakposs(1:NrP)
    !
    real(dp) :: HEnvel(Time_dim), peakStr(1:9),A
    integer :: i, j, k, PeakN,Offst
    integer :: i_loc(1), PeakPos(0:9), NPeakPos(0:9)
    !
    HEnvel(:)=abs(CTime_spectr(:,i_ant, i_chunk))
    PeakN=1
    PeakPos(0)=0
    PeakPos(PeakN)=Time_dim
    Offst=103
    Do i=1,3
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
        !write(2,*) 'PeakPos(:)',PeakPos(0:k)
    enddo ! i
    !
    peakStr=0
    Do j=1,PeakN-1
        peakStr(j) = HEnvel(PeakPos(j))
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
Subroutine MergeFine
! Merge values of FineOffset with input from 'FineCalibrations.dat'
    use DataConstants, only : Station_nrMax, Ant_nrMax, Calibrations, RunMode  ! , DataFolder
    use Chunk_AntInfo, only : Unique_StatID, Nr_UniqueStat, Unique_SAI, Tot_UniqueAnt
    use FitParams, only : Fit_AntOffset, Fit_TimeOffsetStat, Fit_TimeOffsetAnt
    use unque,only : Double_IR_sort
    use constants, only : dp
    use StationMnemonics, only : Station_ID2Mnem
    Implicit none
    Real(dp) :: Delay, Fine_STDelay(Station_nrMax), UpDate_STDelay(1:Ant_nrMax)=0.
    character(len=5) :: txt, Station_Mnem, Fine_STMnm(Station_nrMax)
    integer :: i, k, nxx, i_fine, SAI, n, i_unq
    Real(dp) :: Fine_AntDelay(1:Ant_nrMax), mean(Station_nrMax)
    Integer :: SAI_AntDelay(1:Ant_nrMax),Nr_AntDelay, Ant(1), i_stat, i_SAI
    INTEGER :: DATE_T(8)
    Character(len=12) :: Date_mn
    Logical :: Old, Core
    !
    CALL DATE_AND_TIME (Values=DATE_T)
    WRITE(Date_mn,"(I4,I2.2,I2.2, I2.2,I2.2)") &
          DATE_T(1),DATE_T(2),DATE_T(3),(DATE_T(i),i=5,6)
    !
   i_fine=0
   Fine_AntDelay=0.
   SAI_AntDelay=0
   Old=.false.
   If(trim(Calibrations).eq.'') Old=.true.
   If(Old) then
      Open(unit=9,STATUS='old',ACTION='read', FILE = 'Book/AntCalibrations.dat', IOSTAT=nxx) ! trim(DataFolder)//
   Else
      write(2,*) 'Calibration input from: ',trim(Calibrations)
      Open(unit=9,STATUS='old',ACTION='read', FILE = 'Book/'//trim(Calibrations), IOSTAT=nxx)
   Endif
   If(nxx.ne.0) then
     write(2,*) 'problems with file=','Book/AntCalibrations.dat', ' or ' , trim(Calibrations)
     write(*,*) 'problems with file=','AntCalibrations.dat'
   else
     Do
         read(9,*,IOSTAT=nxx) SAI, Delay
         if(nxx.ne.0) then
             If(Old) close(unit=9)
             exit
         else
             i_fine=i_fine+1
             !write(2,*) i_fine, SAI, Delay
             SAI_AntDelay(i_fine) = SAI
             Fine_AntDelay(i_fine)= Delay
         endif
     enddo
   endif
   !
   UpDate_STDelay=0.
   Nr_AntDelay=i_fine
   If(Fit_AntOffset) then
      Do i_stat=1,Nr_UniqueStat  ! set fine-offset to obtained value from previous fit
        !Stat_ID=Unique_StatID(i_stat)
        !write(*,*) 'i-asai',Tot_UniqueAnt(i_stat-1)+1,Tot_UniqueAnt(i_stat)
        Do i_SAI=Tot_UniqueAnt(i_stat-1)+1,Tot_UniqueAnt(i_stat)      ! Get antenna number from the Unique_Antennas list
            SAI=Unique_SAI(i_SAI)  !     SAI=1000*station_ID + Ant_ID
            Ant=MAXLOC(SAI_AntDelay(1:i_fine), MASK = SAI_AntDelay(1:i_fine) .eq. SAI)
            !write(2,*) 'i_SAI=',i_SAI,SAI,SAI_AntDelay(Ant(1)),Ant(1)
            If(SAI_AntDelay(Ant(1)) .eq. SAI) then  ! this antenna was already in the list
               UpDate_STDelay(Ant(1))= Fit_TimeOffsetAnt(i_SAI)
            else  ! append this antenna to the end of the list
               Nr_AntDelay=Nr_AntDelay+1
               If(Nr_AntDelay .gt. Ant_nrMax) Then
                  write(2,*) 'extent length antenna-delay list, ',Ant_nrMax,' is too small'
                  stop 'MergeFine'
               Endif
               SAI_AntDelay(Nr_AntDelay) = SAI
               UpDate_STDelay(Nr_AntDelay)= Fit_TimeOffsetAnt(i_SAI)
            endif
            !write(2,*) 'UpDate_STDelay',i_SAI,SAI,Nr_AntDelay, UpDate_STDelay(Nr_AntDelay)
        Enddo
      Enddo
   Endif    ! (Fit_AntOffset)
   ! Sort antennas
   Fine_AntDelay(1:Nr_AntDelay)=Fine_AntDelay(1:Nr_AntDelay) + UpDate_STDelay(1:Nr_AntDelay)
   Call Double_IR_sort(Nr_AntDelay,SAI_AntDelay,Fine_AntDelay)
   !
   ! set mean offset per station to zero
   Mean=0.
   i_stat=-1 !NINT(SAI_AntDelay(1)/1000.)
   i_unq=1
   n=1
   Do i=1,Nr_AntDelay
      k=NINT(SAI_AntDelay(i)/1000.)
      If(k.eq.i_stat) then  ! update running sum for this station
         Mean(i_unq)=Mean(i_unq)+ Fine_AntDelay(i)
         n=n+1
      Else
         Mean(i_unq)=Mean(i_unq)/n ! Calculate mean from running sum for previous station
         i_stat=k          ! Store ID new station
         If(i-n.ge.1) then
            Do k=i-n,i-1      ! Set neam antenna delay to zero for previous stattion
             Fine_AntDelay(k)=Fine_AntDelay(k)-Mean(i_unq)
            Enddo
         endif
         Ant=MAXLOC(Unique_StatID(1:Nr_UniqueStat), MASK = Unique_StatID(1:Nr_UniqueStat) .eq. i_stat) ! FINDLOC(Unique_StatID(1:Nr_UniqueStat), i_stat) !
         i_unq=Ant(1)
         If(i_unq.eq.0) then  ! If station not in the list, add it
            Nr_UniqueStat=Nr_UniqueStat+1
            If(Nr_UniqueStat.gt.Station_nrMax) then
               write(2,*) 'nr of unique stations exceeded for',i_stat
               stop 'MergeFine'
            Endif
            Unique_StatID(Nr_UniqueStat)= i_stat ! add to the unique station list
            Fit_TimeOffsetStat(Nr_UniqueStat)=0.
            i_unq=Nr_UniqueStat
         Endif
         n=1
         Mean(i_unq)=Fine_AntDelay(i)  ! Start running sum for this station
      Endif
   Enddo  !   i=1,Nr_AntDelay
   Mean(i_unq)=Mean(i_unq)/n  ! for the last one
   Do k=Nr_AntDelay-n,Nr_AntDelay
    Fine_AntDelay(k)=Fine_AntDelay(k)-Mean(i_unq)
   Enddo
   ! get unique filename
   Open(unit=19,STATUS='unknown',ACTION='write', FILE = 'Book/Calibrations'//Date_mn//'.dat', IOSTAT=nxx)  !  trim(DataFolder)//
   write(2,*) 'New: Calibrations="','Calibrations'//Date_mn//'.dat','"'  ! trim(DataFolder)//
   Do i=1,Nr_AntDelay
      k=NINT(SAI_AntDelay(i)/1000.)
      Call Station_ID2Mnem(k,Station_Mnem)
      If((RunMode.eq.1) .and. (Station_Mnem(1:2) .eq. 'CS')) Then
         Fine_AntDelay(i)=0.
      Endif
      !Fine_AntDelay(i)= Fine_AntDelay(i) +UpDate_STDelay(i)
      write(19,*) SAI_AntDelay(i),Fine_AntDelay(i)
   Enddo
   write(19,*) '============ =========== =============== ======='
   !Do i=1,Nr_UniqueStat
   !   write(2,*) '',i,Unique_StatID(i),Mean(i)
   !enddo
   !close(unit=19)
    !
    i_fine=0
    nxx=0
    If(Old) Open(unit=9,STATUS='old',ACTION='read', FILE = 'Book/FineCalibrations.dat', IOSTAT=nxx) ! trim(DataFolder)//
    If(nxx.ne.0) then
        write(2,*) 'problems with file=','Book/FineCalibrations.dat'
    else
        Do
            read(9,*,IOSTAT=nxx) txt, Delay
            if(nxx.ne.0) then
                close(unit=9)
                exit
            else
                ! write(2,*) txt, Delay
                i_fine=i_fine+1
                Fine_STMnm(i_fine) = txt
                Fine_STDelay(i_fine)= Delay
            endif
        enddo
    endif
    !
    Do k=1,Nr_UniqueStat
        If(Unique_StatID(k).le. 0) exit
        Call Station_ID2Mnem(Unique_StatID(k),Station_Mnem)
        !write(2,*) 'mnem',k,Unique_StatID(k),Station_Mnem
        core=((RunMode.eq.1) .and. (Station_Mnem(1:2) .eq. 'CS'))  ! zero the calibration timings for the core stations
        nxx=1
        Do i=1,i_fine
            If(Station_Mnem .eq. Fine_STMnm(i)) then
                UpDate_STDelay(i)=Fit_TimeOffsetStat(k) + Mean(k)
                Fine_STDelay(i)  =Fine_STDelay(i) + UpDate_STDelay(i)
                If(core) Fine_STDelay(i) = 0
               !write(2,*) 'i',i,UpDate_STDelay(i)
                nxx=0
                exit
            endif
        enddo
        If(nxx.ne.0) then       ! New station that was not in the original file 'FineCalibrations.dat'
            i_fine=i_fine+1
            Fine_STMnm(i_fine) = Station_Mnem
            UpDate_STDelay(i_fine)=Fit_TimeOffsetStat(k)+ Mean(k)
            Fine_STDelay(i_fine)  = UpDate_STDelay(i_fine)
            If(core) Fine_STDelay(i_fine) = 0
            !write(2,*) 'i_fine',i_fine,UpDate_STDelay(i_fine)
        endif
    enddo
    !
    write(2,"(' station',1x,'NewFineCalibrations',5x,'Updated with')")
    Do i=1,i_fine
        write(2,"(1x,A5,3F13.3)") Fine_STMnm(i), Fine_STDelay(i), UpDate_STDelay(i)
         write(19,"(1x,A5,F14.4)") Fine_STMnm(i), Fine_STDelay(i)
    enddo
    Close(unit=19)
    !
    Return
End Subroutine MergeFine
!----------------------------------
Subroutine GetStationFitOption(FP_s, FitNoSources)
   use DataConstants, only : Station_nrMax
   use FitParams, only : N_FitPar_max
   use Chunk_AntInfo, only : Nr_UniqueStat, Unique_StatID
   use StationMnemonics, only : Station_Mnem2ID
   Implicit none
   integer, Intent(OUT)::  FP_s(0:N_FitPar_max)
   Logical, Intent(out) ::  FitNoSources
   Character(len=5) :: FP_MNem(0:N_FitPar_max)
   logical,dimension(Station_nrMax) :: mask
   character*180 :: lname
   character*4 :: option
   Integer :: i, nxx, k, N_FitPar
   !
   FP_s(:)=0
   FP_MNem(:)=' '
   FitNoSources=.false.
   N_FitPar=0
   Call GetNonZeroLine(lname)
   write(2,*) 'timing offset-fit option="',lname,'"'
   read(lname,*,iostat=nxx) option, FP_MNem(1:N_FitPar_max) ! option=abut,only
   Do i=1,N_FitPar_max
      If(FP_Mnem(i).eq.'     ') exit
      If(TRIM(FP_Mnem(i)).eq.'!') exit
      If(FP_Mnem(i).eq.'NoSrc') Then
          FitNoSources=.true.
          cycle
      EndIf
      Call Station_Mnem2ID(FP_Mnem(i),k)
      If(k.eq.0) exit
      N_FitPar=N_FitPar+1
      FP_s(N_FitPar)=k
      !write(2,*) i,FP_Mnem(i),FP_s(i)
   Enddo
   !read(lname,*,iostat=nxx) option, FP_s(1:N_FitPar_max) ! option=abut,only
   If(option .eq. 'abut') then
      mask=.true.
      Do i=1,Nr_UniqueStat
          if ( any(Unique_StatID(i)==FP_s) ) mask(i) = .false.
      enddo
      !write(2,*) 'mask=',mask(1:size(Station_IDs))
      !write(2,*) 'Station_IDs=',Station_IDs(1:size(Station_IDs))
      FP_s(1:Nr_UniqueStat)=pack(Unique_StatID, mask)
      !write(2,*) 'FP_s',FP_s(0:size(Station_IDs))
   endif
   Return
   !
End Subroutine GetStationFitOption
