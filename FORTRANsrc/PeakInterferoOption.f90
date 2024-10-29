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
!  ToDo:
!     - Center windows in time
!     - Apply intensity optimization using fine voxels
!     - Apply intensity optimization using chi-square optimization for cross correlation with beamed result
!     - Generalize read-in to wotk for TRI-D imager results (do not use peak ID)
!     -
   use constants, only : dp, sample_ms
   use DataConstants, only : PeakNr_dim, ChunkNr_dim, RunMode, Time_Dim
   use DataConstants, only : DataFolder, OutFileLabel
   use Chunk_AntInfo, only : Unique_SAI, StartT_sam, Tot_UniqueAnt, Ant_Stations, TimeBase
   use Chunk_AntInfo, only : Unique_StatID,  Nr_UniqueStat, Nr_UniqueAnt, N_Chunk_max
   use DataConstants, only : Ant_nrMax
   use ThisSource, only : ChunkNr, SourcePos, PeakPos, Dual
   use ThisSource, only : PeakNr, PeakNrTotal, PeakChiSQ
   use FitParams, only : N_FitPar, N_FitPar_max, Nr_TimeOffset, WriteCalib
   use FitParams, only : Fit_TimeOffsetStat, Fit_TimeOffsetAnt, Fit_AntOffset, N_FitStatTim, FitParam, X_Offset
   Use Interferom_Pars, only :  TIntens000
   Use Interferom_Pars, only :  N_Smth, N_fit, smooth, i_chunk, BoundingBox
   Use Interferom_Pars, only : IntFer_ant,  Nr_IntferCh ! the latter gives # per chunk
   use Interferom_Pars, only : Alloc_EInterfCalib_Pars
   use Interferom_Pars, only :  MaxSmPow, MaxSmPowLoc
   Use Interferom_Pars, only : N_best, Grd_best, Int_best
   use StationMnemonics, only : Statn_ID2Mnem !, Station_Mnem2ID
   ! use LOFLI_Input, only : ReadSourceTimeLoc
   use FFT, only : RFTransform_su,DAssignFFT, RFTransform_CF2CT
   use GLEplots, only : GLEplotControl
   !Use Calibration, only : WriteCalibration ! was MergeFine
   Implicit none
   Real(dp) :: SourceGuess(3,N_Chunk_max)
   integer, save ::  FP_s(0:N_FitPar_max)!, FP(0:4)
   real ( kind = 8 ) :: X(N_FitPar_max)
   !
   integer :: j, k, i_ant, i_eo, i_stat, StatID, nxx ! i,
   integer :: i_loc(1), i_Peak, Sampl_shiftA
   Integer :: Npix(1:3,1:2), N_Int
   Real(dp) :: dpix(1:3), Spread_Grd, Spread_Dist
   Logical :: FitNoSources
   character*10 :: StochOption
   character*80 :: lname
   real(dp) :: x1,x2,x3,t, Time_width, Space_Spread(1:3), ZeroCalOffset, Del_I(1:3)
   real(dp) :: OrignlSourcePos(1:3,1:PeakNr_dim), SourceIntensity(1:PeakNr_dim), PeakStartChunk_ms(1:PeakNr_dim)
   Integer :: OrignlPeakpos(1:PeakNr_dim), SrcQual(1:10,1:PeakNr_dim)
   !real(dp) :: dt,B, MTC, dtI
   Integer :: i_SAI, OutUnit
   Logical :: Chi2Fit=.false. , ContourPlot
   character(len=2) :: txt
   Real(dp) :: time, x_plot, y_plot, z_plot
   Real(dp), external :: tShift_ms   !
   !
   !
   Dual=.true.
   Call ReadFlashImageDat(PeakStartChunk_ms)
   Call AntFieParGen()
   N_best=40
   !
   Call GLEplotControl(SpecialCmnd='echo on')  ! to generate unit=10 with the proper name
   Do i_Peak=1,PeakNrTotal
      SrcQual(10,i_Peak)=0
      SrcQual(5,i_Peak) =0
      OrignlSourcePos(1:3,i_Peak) = SourcePos(1:3,i_Peak)
      OrignlPeakpos(i_peak) = Peakpos(i_peak)
      write(txt,"(I2.2)") i_Peak
      write(*,*) 'Peak#', i_Peak
      !
      write(2,"(1x,A,i4,A,F11.6,A,2(F9.4,','),F9.4,A)") '++ EI_GridPeakFineSearch: Peak',i_peak,', t_Ref.ant=',&
         PeakStartChunk_ms(i_peak)+PeakPos(i_Peak)*sample_ms-TimeBase,'[ms], at (N,E,h)=(',SourcePos(:,i_Peak)/1000.,') [km]'
      !
      dpix(1:3) = (/ 1.d0, 2.d0, 4.d0 /)
      Npix(1:3,2)=(/ 0, 0, 0 /)
      Npix(1:3,1)=-Npix(1:3,2)
      ContourPlot=.false.
      Call PeakInterferometerRun(PeakStartChunk_ms, i_Peak, dpix, Npix, ContourPlot )
      !
      !
      !write(2,*) 'TIntens000(1:2*N_smth+1)', TIntens000(1:2*N_smth+1)
      Call WindowShift(N_Smth, smooth, TIntens000(1:2*N_smth+1), Sampl_shiftA, txt)
      Write(2,*) 'Peak#',i_peak, 'time shift:',Sampl_shiftA,' ---------------------------------------'
      Peakpos(i_Peak)=Peakpos(i_Peak)+ Sampl_shiftA
      If(abs(Sampl_shiftA) .eq. N_Smth )  SrcQual(5,i_Peak)=1
      !
      !Call PeakCurtain(i_Peak)
      !
      Write(2,*) 'Peak#',i_peak,  'Intensity distribution, round 1: ---------------------------------------'
      dpix(1:3) = (/ 2.d0, 3.d0, 6.d0 /)
      dpix(1:3) = dpix(1:3)/2.
      Npix(1:3,2)=(/ 30, 20, 10 /)
      Npix(1:3,1)=-Npix(1:3,2)
      !i_chunk=ChunkNr(i_Peak)
      ContourPlot=.false.
      Call PeakInterferometerRun(PeakStartChunk_ms, i_Peak, dpix, Npix, ContourPlot )
      write(2,"('Peak#',I2, ', Intensity=',F11.5, ', (N,E,h)=',3F11.5,'[km], diffs:', 3F8.3,'[m]')") i_peak, &
         MaxSmPow(1), MaxSmPowLoc(1:3,1)/1000., MaxSmPowLoc(1:3,1)-SourcePos(1:3,i_Peak)
      SourcePos(1:3,i_Peak)=MaxSmPowLoc(1:3,1)
      If(MaxSmPow(1) .gt. 0.0) Then ! otherwise border case
         SrcQual(10,i_Peak)=SrcQual(10,i_Peak) + 1
         !Call WindowShift(N_Smth, smooth, TIntens000(1:2*N_smth+1), Sampl_shiftA, txt)
         !Peakpos(i_Peak)=Peakpos(i_Peak)+ Sampl_shiftA
      Else
         SrcQual(10,i_Peak)=SrcQual(10,i_Peak) - 1
         write(2,*) '********************** Borderline case, change in source position ********************************'
      EndIf
      !write(2,*) MaxSmPow(1),' ; ', MaxSmPowLoc(1:3,1)
      !
      Write(2,*) 'Peak#',i_peak,  'Intensity distribution, round 2: ---------------------------------------'
      If(SrcQual(10,i_Peak).gt.0) dpix(1:3) = dpix(1:3)/2.
      Npix(1:3,2)=(/ 30, 20, 10 /)
      Npix(1:3,1)=-Npix(1:3,2)
      !i_chunk=ChunkNr(i_Peak)
      ContourPlot=.true.
      Call PeakInterferometerRun(PeakStartChunk_ms, i_Peak, dpix, Npix, ContourPlot )
      write(2,"('Peak#',I2, ', Intensity=',F11.5, ', (N,E,h)=',3F11.5,'[km], diffs:', 3F8.3,'[m]')") i_peak, &
         MaxSmPow(1), MaxSmPowLoc(1:3,1)/1000., MaxSmPowLoc(1:3,1)-SourcePos(1:3,i_Peak)
      !
      If(MaxSmPowLoc(3,1) .gt. 1.0) Then ! otherwise border case
         SourcePos(1:3,i_Peak)=MaxSmPowLoc(1:3,1)
         SrcQual(10,i_Peak)=SrcQual(10,i_Peak) + 1
         !Call WindowShift(N_Smth, smooth, TIntens000(1:2*N_smth+1), Sampl_shiftA, txt)
         !Peakpos(i_Peak)=Peakpos(i_Peak)+ Sampl_shiftA
      Else
         write(2,*) '********************** no change in source position ********************************'
         SrcQual(10,i_Peak)=SrcQual(10,i_Peak) - 1
      EndIf
      Call IntensitySpread(Grd_best(1,1,1), Int_best(1,1), dpix, N_best, N_Int, Spread_Grd, Spread_Dist)
      !write(2,*) '!Source Spread:', N_Int, Spread_Grd, Spread_Dist, Int_best(N_Int,1)/Int_best(1,1)
      SrcQual(1,i_Peak)=NINT(Spread_Grd)
      SrcQual(2,i_Peak)=NINT(Spread_Dist)
      SrcQual(3,i_Peak)=NINT(Int_best(1,1)/(Int_best(1,1)-Int_best(N_Int,1)))
      SrcQual(4,i_Peak)=NINT(Spread_Grd*Int_best(1,1)/(Int_best(1,1)-Int_best(N_Int,1)))
      SourceIntensity(i_Peak)=MaxSmPow(1)
      !
      ! prepare & produce curtain plot
      Call PeakCurtain(i_Peak)
      !
      Write(2,*) 'Peak#',i_peak,  'finalized ================================================'
      !Call GLEplotControl(SpecialCmnd='ls') !  ${FlashFolder}/
      !
   EndDo
   !
   Do i_Peak=1,PeakNrTotal
      write(2,"('Peak#',I2, ', Q=',6I5,A,F11.6, A,F6.1, ', (N,E,h)=',3F11.5,'[km], diffs:', I3, 3F8.3,'[m]')") &
         i_peak, SrcQual(1:5,i_Peak), SrcQual(10,i_Peak), &
         ', rel t ref ant=', PeakStartChunk_ms(i_peak)+PeakPos(i_Peak)*sample_ms -TimeBase, &
         '[ms],  Intensity=',SourceIntensity(i_Peak), SourcePos(1:3,i_Peak)/1000., Peakpos(i_Peak)-OrignlPeakpos(i_Peak), &
         SourcePos(1:3,i_Peak)-OrignlSourcePos(1:3,i_Peak)
      !
   EndDo
   !
   OPEN(UNIT=29,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'IntfPkQPow.dat')
   write(29,"(6F8.2,2F9.3,A,F7.1,i3,' 0')") BoundingBox(1:2,1:4), ' NoBox ', TimeBase, 0 ! gleRead:  xMin xMax yMin yMax zMin zMax tMin tMax ZmBx$ t_offst
   Write(29,*) '0 ',PeakNrTotal,' 0 0 ',TRIM(OutFileLabel), 10.,' 0 0 0 0 0 0 0 ',0.1 ,  '1.0 !' ! gleRead:  NTracks EventNr Q t_start label$ AmplitudePlot a1 b1 c1 a2 b2 c2 d2
   Do i_Peak=1,PeakNrTotal
      write(2,"('Peak#',I2, ', Q=',I2,A,F11.6, A,F6.1, ', (N,E,h)=',3F11.5,'[km], diffs:', I3, 3F8.3,'[m]')") &
         i_peak, SrcQual(10,i_Peak),', abs t ref ant=', PeakStartChunk_ms(i_peak)+PeakPos(i_Peak)*sample_ms, &
         '[ms],  Intensity=',SourceIntensity(i_Peak), SourcePos(1:3,i_Peak)/1000., Peakpos(i_Peak)-OrignlPeakpos(i_Peak), &
         SourcePos(1:3,i_Peak)-OrignlSourcePos(1:3,i_Peak)
      !
      !
      Time=PeakStartChunk_ms(i_peak)+PeakPos(i_Peak)*sample_ms - tShift_ms(SourcePos(:,i_Peak)) - TimeBase ! StartTime_ms
      x_plot=SourcePos(2,i_Peak)/1000.
      y_plot=SourcePos(1,i_Peak)/1000.
      z_plot=SourcePos(3,i_Peak)/1000.
      write(29,"(i6,' ',4(f11.5,' '),g13.6,' ',3f6.2,' ',5g13.3)") i_Peak, Time, &
         x_plot, y_plot, z_plot, ABS(SourceIntensity(i_Peak))
      !
   EndDo
   Close(Unit=29)
   !
   Call GLEplotControl(PlotType='SourcesPlot', PlotName='IntfPkQ'//TRIM(OutFileLabel), &
      PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel)//'IntfPkQPow' )
   !
   Call GLEplotControl( Submit=.true.)
   !
   write(2,*) '#   Qual,  Int , t_src [ms]  , x=east  ,  y=north  ,  z=height  [km]'
   Do i_Peak=1,PeakNrTotal
      write(2,"(I2, ', ',F11.6,', ',3F11.5,', ',I2, ', ', F6.1)") &
         i_peak, PeakStartChunk_ms(i_peak)+PeakPos(i_Peak)*sample_ms - tShift_ms(SourcePos(:,i_Peak)) - TimeBase,  &
         SourcePos(2,i_Peak)/1000., SourcePos(1,i_Peak)/1000., SourcePos(3,i_Peak)/1000. , &
         SrcQual(10,i_Peak), SourceIntensity(i_Peak)
      !
   EndDo
   !
   Return
   !
   !====================================================================
   !
End Subroutine PeakInterferoOption
!-----------------------------------------
Subroutine PeakCurtain(i_peak)
   use Constants, only : dp !,sample,c_mps, pi, sample_ms
   use DataConstants, only : Ant_nrMax
   use Interferom_Pars, only :  Nr_IntFerMx, Nr_IntferCh ! the latter gives # per chunk
   use Interferom_Pars, only : Cnu_p0, Cnu_t0, Cnu_p1, Cnu_t1, IntfNuDim, AntPeak_OffSt
   Use Interferom_Pars, only : N_smth, N_fit, i_chunk
   use ThisSource, only : SourcePos, PeakPos, ExclStatNr
   Use Interferom_Pars, only : Alloc_EInterfImag_Pars, DeAlloc_EInterfImag_Pars
   Implicit none
   Integer, intent(in) :: i_Peak
   Real(dp) :: FitDelay(1:Ant_nrMax)
   Real(dp) :: ChiSq, DelChi(-N_smth:+N_smth,1:2*Nr_IntFerMx)
   Integer :: IntfBase, Outpt
   Character(len=8) :: Label  ! should not be longer
   !
   Outpt=2  ! will generate also a curtain plot
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
!=================================
Subroutine WriteInterfRslts(i_Peak, SrcQual) ! original called from EIFitter for mode=8 only
   use Constants, only : dp,sample,c_mps, pi, sample_ms
   use DataConstants, only : DataFolder, OutFileLabel
   Use Interferom_Pars, only : StI, StI12, StQ, StU, StV, StI3, StU1, StV1, StU2, StV2, P_un, P_lin, P_circ
   Use Interferom_Pars, only : dStI, dStI12, dStQ, dStU, dStV, dStI3, dStU1, dStV1, dStU2, dStV2, Chi2pDF, BoundingBox
   use ThisSource, only : PeakNrTotal, ChunkNr, PeakPos, SourcePos
   use Chunk_AntInfo, only : StartT_sam, TimeBase
   use GLEplots, only : GLEplotControl
   Implicit none
   Integer, intent(in) :: i_Peak, SrcQual(1:10)  ! =29, 28
   Real(dp) :: time, x, y, z
   Real(dp), external :: tShift_ms   !
   !
   !Time=SQRT(sum(SourcePos(:,i_Peak)*SourcePos(:,i_Peak)))
   !Time = Time*Refrac/(c_mps*sample) ! Convert to units of [samples]
   !Time=( StartT_sam(ChunkNr(i_Peak)) + PeakPos(i_Peak) - Time )*Sample*1000. - StartTime_ms
   Time=( StartT_sam(ChunkNr(i_Peak)) + PeakPos(i_Peak) )*sample_ms - tShift_ms(SourcePos(:,i_Peak)) - TimeBase ! StartTime_ms
   x=SourcePos(2,i_Peak)/1000.
   y=SourcePos(1,i_Peak)/1000.
   z=SourcePos(3,i_Peak)/1000.
   !
   If(i_peak.eq.1) then
      OPEN(UNIT=29,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'IntfPkQPow.dat')
      write(29,"(6F8.2,2F9.3,A,F7.1,i3,' 0')") BoundingBox(1:2,1:4), ' NoBox ', TimeBase, 0 ! gleRead:  xMin xMax yMin yMax zMin zMax tMin tMax ZmBx$ t_offst
      Write(29,*) '0 ',PeakNrTotal,' 0 0 ',TRIM(OutFileLabel), 10.,' 0 0 0 0 0 0 0 ',0.1 ,  '1.0 !' ! gleRead:  NTracks EventNr Q t_start label$ AmplitudePlot a1 b1 c1 a2 b2 c2 d2
   EndIf
   !
   If(SrcQual(10).gt.0)  Then
      write(29,"(i6,' ',4(f11.5,' '),g13.6,' ',3f6.2,' ',5g13.3)") i_Peak, Time, &
         x, y, z, StI12, StQ/StI, StU/StI, StV/StI, StI3 ,P_lin,( ATAN2(StU,StQ)/2. )*180./pi
   EndIf
   !
   !
   If(i_Peak.eq.PeakNrTotal) Then
      !Close(Unit=28)
      Close(Unit=29)
      Call GLEplotControl(PlotType='SourcesPlot', PlotName='IntfPkQ'//TRIM(OutFileLabel), &
         PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel)//'IntfPkQPow' )
      !Call GLEplotControl(PlotType='EIPolariz', PlotName='IntfPkQPol'//TRIM(OutFileLabel), &
      !   PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel)//'IntfPkQPolrz' )
   EndIf
   !
   Return
End Subroutine WriteInterfRslts
!=========================
Subroutine ReadFlashImageDat(PeakStartChunk_ms)
   !
   use constants, only : dp,sample, c_mps, sample_ms
   use DataConstants, only : PeakNr_dim, ChunkNr_dim, DataFolder, EdgeOffset, Time_dim
   use Chunk_AntInfo, only : Station_nrMax, StartT_sam, TimeBase
   use Chunk_AntInfo, only : ChunkStTime_ms, ChunkFocus
   Use Interferom_Pars, only : BoundingBox
   use ThisSource, only : SourcePos,  TotPeakNr !NrP, t_ccorr,
   use ThisSource, only : PeakNrTotal !,PeakNr,  PlotCCPhase, Safety, Nr_corr
   use ThisSource, only : PeakPos, ChunkNr
   use StationMnemonics, only : Statn_ID2Mnem, Station_Mnem2ID
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
   Real(dp), intent(out) :: PeakStartChunk_ms(1:PeakNr_dim)
   integer :: i_eo, i_chunk, i, k!, i_c ,i_ca  !,i_eoa
   Character(len=5) :: Station_Mnem, ExclStMnem(1:Station_nrMax)
   integer :: i_Peak, nxx, InUnit
   Integer, parameter :: stdin=5
   character*180 :: lname
   character*50 :: FlashFile
   character*10 :: txt
   character*35 :: Frmt
   real(dp) :: x1,x2,x3,t, A, TimeShft
   Real(dp), external :: tShift_smpl, tShift_ms   !
   !
   Allocate(ChunkStTime_ms(1:ChunkNr_dim))
   Allocate(ChunkFocus(1:3,1:ChunkNr_dim))
   i_chunk=1
   i_peak=0
   TotPeakNr(0,:)=0
   !
   Call GetNonZeroLine(lname)
   Read(lname,*,iostat=nxx) BoundingBox(1:2,1:4),txt ! StartTime_ms
   If(nxx.ne. 0) Then
      Read(lname,*) FlashFile
      InUnit=17
      Open(Unit=InUnit, STATUS='old',ACTION='read', FILE=TRIM(datafolder)//TRIM(FlashFile), IOSTAT=nxx)
      If(nxx.ne.0) Then
         Write(2,*) 'Not found!, source information file: "',TRIM(datafolder)//TRIM(FlashFile),'"'
         Stop 'SourcesDataFile not found'
      Else
         Write(2,*) 'Source information read from file: "',TRIM(datafolder)//TRIM(FlashFile),'"'
         Read(InUnit,*,iostat=nxx) BoundingBox(1:2,1:4),txt,TimeBase ! StartTime_ms
         Read(InUnit,*) txt
         Read(InUnit,*) txt
      EndIf
   Else
      InUnit=stdin
   EndIf
   !StartTime_ms=TimeBase
   write(2,*) 'TimeBase', TimeBase
   !
   Do while (i_peak.lt.PeakNr_dim)
      Read(InUnit,*,iostat=nxx) k,t,x2,x1,x3 ! [--,tsource_ms,N_km,E_km,h_km]
      !write(2,*) '! ReadFlashImageDat1:',k,t,x2,x1,x3
      If(nxx.ne.0) Then
         exit
      EndIf
      !
      i_peak=i_peak+1
      SourcePos(1,i_Peak)=x1*1000.
      SourcePos(2,i_Peak)=x2*1000.
      SourcePos(3,i_Peak)=x3*1000.
      TimeShft=tShift_ms(SourcePos(:,i_Peak))
      PeakStartChunk_ms(i_peak)=TimeBase+t + TimeShft - Time_dim*sample_ms/2.  ! in ms, the start of the chunk for this peak
      ChunkNr(i_Peak)=1
      PeakPos(i_Peak) = Time_dim/2 ! NINT((TimeBase+t)/(Sample*1000.d0) - StartT_sam(ChunkNr(i_Peak))+ TimeShft) ! better be positive, <Time_dim
      write(2,"(I2,A,I6, A,F11.6,A,3F9.4, A,F12.6, A,F9.6)") i_peak, ', PeakPos', PeakPos(i_Peak) &
         ,' @ t_ref=', PeakStartChunk_ms(i_peak)+PeakPos(i_Peak)*sample_ms-TimeBase,'[ms], position=', &
         x1, x2, x3, ', t_source=', t, '[ms], t-shift from ref=', TimeShft*sample_ms
      !
      !write(2,*) '!',k,i_peak,i_chunk, StartT_sam(i_chunk), x1, x2, x3, PeakPos(i_Peak)
      !write(2,*) '! ReadFlashImageDat[sampl]:', PeakPos(i_Peak), (StartTime_ms+t)/(Sample*1000.d0), StartT_sam(ChunkNr(i_Peak))
      !write(2,*) '! ReadFlashImageDat[ms]:',  (StartTime_ms+t)+ TimeShft*(Sample*1000.d0),  &
      !      PeakPos(i_Peak)*(Sample*1000.d0)
      !write(2,*) t/(Sample*1000.), StartT_sam(ChunkNr(i_Peak)), TimeShft

   EndDo
   PeakNrTotal=i_peak
   ChunkNr_dim=1
   Close(Unit=17)
   !
   Return
End Subroutine ReadFlashImageDat
!================================
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
Subroutine WindowShift(N_Smth, smooth, Intens, Sampl_shift, txt)
   use constants, only : dp
   use DataConstants, only : DataFolder, OutFileLabel
   Implicit none
   Integer, intent(in) :: N_Smth
   Real(dp), intent(in) :: smooth(-N_Smth:N_Smth), Intens(-N_Smth:N_Smth)
   Integer, intent(out) :: Sampl_shift
   Character(len=*), intent(in) :: txt
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
            write(2,*) '!t-shift option3:', i, A_s
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
   OPEN(UNIT=30,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'PkIntf-'//txt//'.csv')
   write(30,"(4(I3,','), 3(A,',') )")  N_Smth, Sampl_shift, i_d, i_u, ' ! '
   write(30,"(A)")  ' ! i,    windowing,  Intensity,  Intensity for polarization directio = 1, 2, 3 '
   Do i=-N_Smth,N_Smth
      write(30,"(I3,',',5(g12.4,','))") i,smooth(i), Intens(i)
   EndDo
   Close(UNIT=30)
   !
   Return
End Subroutine WindowShift
!---------------------------------------
! ======================================================================
! ======================================================================
Subroutine IntensitySpread(Grd_best, Int_best, Grd, N_best, N_Int, Spread_Grd, Spread_Dist)
   use constants, only : dp
   Implicit none
   Integer, intent(in) :: Grd_best(1:3,1:N_best), N_best
   Real(dp), intent(in) :: Int_best(1:N_best), Grd(1:3)
   Integer, intent(out) :: N_Int
   Real(dp), intent(out) :: Spread_Grd, Spread_Dist
   Integer :: i,j
   Real(dp) :: Frac_int=0.9, Max_grd, Max_dist, DGrd, DGrd_i, DGrd_j, Dist, curv_dist, curv_Grd, curv
   Spread_Grd=0.
   Spread_Dist=0.
   N_Int=1
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
   EndIf
   !
   write(2,"(A, 100F10.3)")    '!Intensity: ',Int_best(:)
   write(2,"(A, 100(3I3,';'))") '!IntensityGrd:',Grd_best(:,:)
   write(2,*) 'not-weighted Spread:', Spread_Grd,'[grd],', Spread_Dist,'[m]', N_Int, Int_best(N_Int)/Int_best(1)
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
   write(2,*) 'maximal distances=',sqrt(Max_grd),'[grd],', sqrt(Max_Dist),'[m],', N_Int
   !write(2,*) 'curvatures=', Curv/(0.5*(N_Int-2)*(N_Int-1)), Curv_Grd/(N_Int-1), Curv_Dist/(N_Int-1)
   write(2,*) 'Int-Weighted Spread=', & ! 1/sqrt(Curv/(0.5*(N_Int-2)*(N_Int-1))),'[grd],',
      1/sqrt(Curv_Grd/(N_Int-1)),'[grd],', 1/sqrt(Curv_Dist/(N_Int-1)),'[m]'
   !
   Return
End Subroutine IntensitySpread
!=================
Subroutine PeakInterferometerRun(PeakStartTref_ms, i_Peak, dpix, Npix, ContourPlot )
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
   Real(dp), intent(in) :: dpix(1:3), PeakStartTref_ms(*)
   Integer, intent(in) :: Npix(1:3,1:2)
   logical, intent(in) :: ContourPlot
   integer :: i  !,j, k, i_s, i_loc(1), i_ant, j_corr, i_eo, Station, j_IntFer, i_nu,  nxx
   character(len=3) :: txt
   character(len=100) :: OutFileLabel_original
   Real(dp), save :: CenLoc_old(1:3)
   Integer, save :: Start_old
   Integer :: Date_T(8)
   Real(dp), external :: tShift_ms
   !    .
   i_chunk=1
   StartT_sam(1)=PeakStartTref_ms(i_peak)/sample_ms
   ! Get grid:
   Polar=.false.
   d_loc(1:3)=dpix(1:3)
   N_pix(1:3,1:2)=Npix(1:3,1:2)
   !
   SumWindw=2*N_smth+1  ! should be even to have window centered at peak of interest
   NrPixSmPowTr=(SumWindw-1)/N_smth-1  ! The number of fine slices for TRI-D imaging
   SumStrt=PeakPos(i_Peak) - SumWindw/2
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
      write(txt,"('-',I2.2)") i_peak
      OutFileLabel_original=OutFileLabel
      OutFileLabel=TRIM(OutFileLabel)//txt
   EndIf
   ! main loop over direction antennas
   Call EI_run
   !
   !------------------------------------------
   If(ContourPlot) Then
      write(2,*) 'prepare plot:', TRIM(OutFileLabel), '--------------'
      Call GLEplotControl(PlotType='EIContour', PlotName='EIContour'//TRIM(OutFileLabel), &
               PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel))
      Call GLEplotControl(CleanPlotFile=TRIM(OutFileLabel)//'_EISpec'//'*.csv')
      Call GLEplotControl(CleanPlotFile=TRIM(OutFileLabel)//'Interferometer'//'*.z')
      OutFileLabel=OutFileLabel_original
   EndIf
   !
   Return
    !
End Subroutine PeakInterferometerRun
!-------------------------------------------------------------
!=========================================
!==========================
