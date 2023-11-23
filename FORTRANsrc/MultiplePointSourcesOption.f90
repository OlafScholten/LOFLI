!=================================
!================================
!=================================
Subroutine MultipleDDSourcesOption()
!   Chi-square interferometric imaging of peaks found by the impulsive imager
!     Fit the spectra in all antennas as the superposition of a number of dipole-delta sources
!     A dipole-delta (DD) sources in this context is a point-impulsive source (zero extent in space and time) emitting like an ideal dipole.
!
!  Input: ??
!     -reads a list of sorces?
!
!  Output:
!     - list of dipole-delta parameters for each source
!
!   Procedural steps:
!     1) read source parameters
!     2) determine the range in time and space over whic DD sorces should be fit
!
Copy here from  Subroutine InterferometerRun  (in InterferometryOption.f90)
   use FFT, only : RFTransform_su, DAssignFFT !, RFTransform_CF, RFTransform_CF2CT
   use HDF5_LOFAR_Read, only : CloseDataFiles
   !    .
   ! Get center voxel
   Call ReadSourceTimeLoc(StartTime_ms, CenLoc)  ! general area
   !
   !  Read the time traces
   CALL DATE_AND_TIME (Values=DATE_T)
   WRITE(2,"(1X,I2,':',I2.2,':',I2.2,'.',I3.3,A)") (DATE_T(i),i=5,8), 'initializing'
   i_chunk=1
   StartT_sam(i_chunk)=(StartTime_ms/1000.d0)/sample  ! in sample's
   Call RFTransform_su(Time_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   Call AntennaRead(i_chunk, CenLoc)
   If(Simulation.eq."") Then
      CALL CloseDataFiles()    ! no more reading of antenna data, only close when reading real data
   EndIf
   call DAssignFFT()                   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !
   Call AntFieParGen()
   Call ImpResp_ns()
   !
   !
   !
   Call GetMarkedLine(Mark,lname)
   read(lname,*) CenT_ms, CenLoc(1:3)
   !
   !  Get source region in space and time
   !  setup arrays
   IntfLead=10+(diff(1)*N_pix(1,2)+diff(2)*N_pix(2,2)+diff(3)*N_pix(3,2)) ! Estimate lead and trail buffer length,
   IntfNuDim=SumWindw/2 + IntfLead +1 ! Length of frequency trace used for interference
   i = shiftr( IntfNuDim, 9 )+1
   IntfNuDim=shiftl( i, 9 )  ! Make sure it approches some multiple of 2
   !write(2,"(A,o10,o10,i7)") 'i:',i,IntfNuDim,IntfNuDim
   IntfDim=2*IntfNuDim ! Length of time trace used for interference
   IntfLead=(IntfDim-SumWindw)/2
   IntfBase=SumStrt - IntfLead
   !write(*,"(A,o10,b20)") '(lead,dim):', IntfLead, IntfDim
   If((IntfBase+IntfDim .gt. Time_dim) .or. (IntfBase .lt. 1)) then
      write(2,*) 'dimensions for interferometry out of range: ', IntfBase, IntfBase+IntfDim, Time_dim
      write(2,*) 'input for interferometry: ', SumStrt, SumWindw
      stop 'Intferom dimension problem'
   Endif
   !
   Call MDDSetup(Nr_IntFer, IntfNuDim, CenT_ms, CenLoc)
   !
this should call MDD_run
!=================================
Subroutine MDD_run()
similar to  Subroutine EI_Run  (in EIOption)
!   Chi-square interferometric imaging of peaks found by the impulsive imager
!     Fit the spectra in all antennas as the superposition of a number of dipole-delta sources
!     A dipole-delta (DD) sources in this context is a point-impulsive source (zero extent in space and time) emitting like an ideal dipole.
!
!  Input:
!     CenT_ms, CenLoc(1:3)
!     - Central location and time
!
!   Procedural steps:
!   1) determine the range in time and space over whic DD sorces should be fit
!        How: interactively?  or some standard range (~150 samples?)
!   2) put 3x7 sources on c+(+0-dx,+0-dy,+0-dz,+0-dt) excluding in space ++ and +-
!   3) start fitting round, fitting DD source parameters, excluding non-critical ones??
!        Call MDD_Fitter
!



   Complex(dp), allocatable :: CMCnu(:,:,:), CMTime_Pix(:,:)
   !
   CALL DATE_AND_TIME (Values=DATE_T)
   WRITE(2,"(1X,I2,':',I2.2,':',I2.2,'.',I3.3,A,i1)") (DATE_T(i),i=5,8),&
      ' start Multiple Delta Dipoles (MDD) composite sources fit'
   WRITE(*,"(1X,I2,':',I2.2,':',I2.2,'.',I3.3,A)") (DATE_T(i),i=5,8), achar(27)//'[0m'
   !
   i_chunk=1
   Call EISelectAntennas(i_chunk)  ! select antennas for which there is an even and an odd one.
   Nr_IntFer=Nr_IntFerCh(i_chunk)
   Call Alloc_EInterfImag_Pars
   Allocate( CMCnu(Nr_IntFer,0:IntfNuDim,1:3),  CMTime_pix(1:2*IntfNuDim,1:3) )
   !
   Call MDDSetupSpec(Nr_IntFer, IntfNuDim, CMCnu)
   !Nr_IntFerMx=Nr_IntFerCh(i_chunk) since there is only a single chunk
   !
   !------------------------------------



End Subroutine MDD_run
!=================================
Subroutine WriteInterfRslts(i_Peak)
   use Constants, only : dp,sample,c_mps, pi
   use DataConstants, only : DataFolder, OutFileLabel
   Use Interferom_Pars, only : StI, StI12, StQ, StU, StV, StI3, StU1, StV1, StU2, StV2, P_un, P_lin, P_circ
   Use Interferom_Pars, only : dStI, dStI12, dStQ, dStU, dStV, dStI3, dStU1, dStV1, dStU2, dStV2, Chi2pDF, BoundingBox, StartTime_ms
   use ThisSource, only : PeakNrTotal, ChunkNr, PeakChiSQ,  PeakPos, SourcePos
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
