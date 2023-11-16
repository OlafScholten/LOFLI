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
!Copy here from  Subroutine InterferometerRun  (in InterferometryOption.f90)
   use constants, only : dp, ci, pi, Sample, c_mps, sample_ms
   Use MDD_Pars, only : CR_Tms, C_Tms, C_Loc, T_Range, MDDtDim, MDDnuDim, MDDLabel  ! output
   Use Interferom_Pars, only : StartTime_ms, CenLoc, IntFer_ant, i_chunk, Nr_IntFerCh
   use DataConstants, only : OutFileLabel, DataFolder, Time_Dim !, Cnu_dim,
   use Chunk_AntInfo, only : StartT_sam, RefAnt, Simulation, TimeBase
   use FFT, only : RFTransform_su, DAssignFFT !, RFTransform_CF, RFTransform_CF2CT
   use HDF5_LOFAR_Read, only : CloseDataFiles
   use GLEplots, only : GLEplotControl
   Implicit none
   integer :: i, error  !   j, k, i_s, i_loc(1), i_ant, j_corr, i_eo, Station, j_IntFer, i_nu,  nxx
   Real(dp) :: FitQual, t_shft
   Integer :: Range, Nr_IntFer
   Integer, parameter ::lnameLen=180
   !character(len=5) :: txt
   character(Len=1) :: Mark
   Character(LEN=lnameLen) :: lname
   Integer :: Date_T(8)
   !Character*7 :: PreAmble
   Real(dp), external :: tShift_ms
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
   Call EISelectAntennas(i_chunk)  ! select antennas for which there is an even and an odd one.
   Nr_IntFer=Nr_IntFerCh(i_chunk)
   Call AntFieParGen()  !  Prepare antenna function and gain
   Call ImpResp_ns()    !  Construct Impulse response function
   !
   Call GetMarkedLine(Mark,lname)
   read(lname,*) C_Tms, C_Loc(1:3), Range !  full range in [s] that should be considered
   Call Convert2m(C_Loc(:))
   C_Tms=C_Tms+TimeBase
   t_shft=tShift_ms(C_Loc(:)) ! sqrt(SUM(CenLoc(:)*CenLoc(:)))*1000.*Refrac/c_mps ! in mili seconds due to signal travel distance
   CR_Tms=C_Tms+t_shft  ! time in reference antenna
   Write(2,*) 'CR_Tms:',CR_Tms, C_Tms

   T_Range=(Range/2)*2 +1  ! make sure it is an odd number
   !
   MDDtDim=50+2*T_Range ! lead+tail =window lendth length + 50 for safety,
   MDDnuDim=MDDtDim/2  ! Length of frequency trace used for interference
   i = shiftr( MDDnuDim, 9 )+1
   MDDnuDim=shiftl( i, 9 )  ! Make sure it approches some multiple of 2
   !write(2,"(A,o10,o10,i7)") 'i:',i,IntfNuDim,IntfNuDim
   MDDtDim=2*MDDnuDim ! Length of time trace used for interference
   !write(*,"(A,o10,b20)") '(lead,dim):', IntfLead, IntfDim
   !
   Call MDDSetup()
   MDDLabel='MDD'
   !
   Call MDDFitSU()
   Call MDD_Fitter(FitQual,error)
   Call MDD_write()
   !
   Call MDDClose()
   Return
End Subroutine MultipleDDSourcesOption
!=================================
!=================================

!============================
