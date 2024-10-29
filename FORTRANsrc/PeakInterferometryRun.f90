Subroutine PeakInterferometerRun( PeakPos, SourcePos )
!   Call PeakInterferometerRun( PeakPos(i_Peak), SourcePos(1:3,i_Peak),  )
   use constants, only : dp, ci, pi, Sample, c_mps, sample_ms
   Use Interferom_Pars, only : Polar, NewCenLoc, ChainRun,  AmpltPlot, i_chunk
   Use Interferom_Pars, only : N_pix, d_loc, CenLocPol, CenLoc, PixLoc, IntFer_ant
   Use Interferom_Pars, only : SumStrt, SumWindw, NrSlices, SliceLen
   Use Interferom_Pars, only : N_smth, NrPixSmPowTr
   use DataConstants, only : OutFileLabel, DataFolder, Time_Dim !, Cnu_dim,
   use Chunk_AntInfo, only : StartT_sam, TimeFrame,   NoiseLevel, RefAnt, Simulation
   use ThisSource, only : CurtainHalfWidth, PeakNrTotal, Peak_eo, ChunkNr, Peakpos, SourcePos
   use FFT, only : RFTransform_su, DAssignFFT !, RFTransform_CF, RFTransform_CF2CT
   use HDF5_LOFAR_Read, only : CloseDataFiles
   use GLEplots, only : GLEplotControl
   !use PredefinedTracks, only : PreDefTrackFile
   Implicit none
   Integer, intent(in) :: PeakPos
   !integer :: i,j, k, i_s, i_loc(1), i_ant, j_corr, i_eo, Station, j_IntFer, i_nu,  nxx
   !Integer, parameter ::lnameLen=180
   !character(len=5) :: txt
   !character(Len=1) :: Mark
   !Character(LEN=lnameLen) :: lname
   Integer :: Date_T(8)
   Real(dp), external :: tShift_ms
   !    .
   ! Get grid:
   Polar=.false.
   !StartTime_ms=StartT_sam(i_chunk)*sample*1000.d0  ! in ms &  use Chunk_AntInfo, only : StartT_sam
   d_loc(1:3)=(/ 1.d0, 2.d0, 4.d0 /)
   N_pix(1:3,2)=(/ 30, 15, 7 /)
   N_pix(1:3,1)=-N_pix(1:3,2)
   !
   SumWindw=2*N_smth  ! should be even to have window centered at peak of interest
   NrPixSmPowTr=SumWindw/N_smth-1  ! The number of fine slices for TRI-D imaging
   SumStrt=PeakPos - SumWindw/2
   CenLoc(1:3) = SourcePos(1:3)
   AmpltPlot = 10.
   write(2,"(A,I5,A,I5,A,F6.3)") 'TRI-D window: SumStrt=',SumStrt, ', SumWindw=', SumWindw, ', AmpltPlot=', AmpltPlot
   If(SumStrt.lt.1000) Then
      write(2,*) 'Window starting sample should be 1000'
      Stop 'Window starting sample should be sufficiently large'
   EndIf
   !
   NrSlices=1  !  depreciated, Dec 2021; used in making summed traces in 'EIAnalyzePixelTTrace' and in 'OutputIntfSlices'
   SliceLen=SumWindw/NrSlices ! same as NrSlices
   !
   CALL DATE_AND_TIME (Values=DATE_T)
   WRITE(2,"(1X,I2,':',I2.2,':',I2.2,'.',I3.3,A)") (DATE_T(i),i=5,8), 'initializing'
   i_chunk=1
   TimeFrame=i_chunk  ! not sure this is really used
   !
   !  Read the time traces
   Call RFTransform_su(Time_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   Call AntennaRead(i_chunk, CenLoc)
   CALL CloseDataFiles()    ! no more reading of antenna data
   CALL CloseDataFiles()    ! no more reading of antenna data
   call DAssignFFT()                   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !
   If(NoiseLevel.lt.0) NoiseLevel=0.2
   !
   !   Allocate( Smooth(-N_smth:N_smth) )
   !   Allocate( TIntens_pol(1:3,-N_smth:N_smth) )  !  Intensity trace for different pol directions as seen at the core, For use in "PeakInterferoOption"
   !   Call SetSmooth(N_Smth, Smooth)
   ! ----------------------------------------------------
   ! main loop over direction antennas
   Call EI_run
   !DeAllocate( Smooth )
   !DeAllocate( TIntens_pol )  !  Intensity trace for different pol directions as seen at the core, For use in "PeakInterferoOption"
   !
   !------------------------------------------
   Call GLEplotControl(PlotType='EIContour', PlotName='EIContour'//TRIM(OutFileLabel), &
            PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel), Submit=.true.)
   Call GLEplotControl(CleanPlotFile=trim(DataFolder)//TRIM(OutFileLabel)//'_EISpec'//'*.csv')
   Call GLEplotControl(CleanPlotFile=trim(DataFolder)//TRIM(OutFileLabel)//'Interferometer'//'*.z')
   !
   Return
    !
End Subroutine PeakInterferometerRun
!-------------------------------------------------------------
!=========================================
