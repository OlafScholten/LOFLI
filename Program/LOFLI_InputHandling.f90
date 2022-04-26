MODULE LOFLI_Input

!   USE ISO_FORTRAN_ENV
    use constants, only : dp

   IMPLICIT NONE

   INTERFACE PrintValues
      MODULE PROCEDURE PrintIntArray
      MODULE PROCEDURE PrintIntVal
      MODULE PROCEDURE PrintRealVal
      MODULE PROCEDURE PrintLogiVal
      MODULE PROCEDURE PrintChArray
      MODULE PROCEDURE PrintChVal
   END INTERFACE PrintValues

   CONTAINS

   Subroutine PrintIntArray(Var, Var_name, Label)
      Integer, dimension(:), Intent(in) :: Var
      Character(len=*), Intent(in) :: Var_name
      Character(len=*), Intent(in), optional :: Label
      Integer :: N,i
      Character(len=20) :: VarCh
      N=size(Var)
      write(2,"(1x,A,A)", ADVANCE='NO') TRIM(Var_Name),'= '
      Do i=1,N
         write(VarCh,*) Var(i)
         If(Var(i).eq.Var(N)) exit
         write(2,"(A,A)", ADVANCE='NO') Trim(ADJUSTL(VarCh)),', '
      Enddo
      If(present(Label)) Then
         write(2,"(I3,A,A,A,5x,'! ',A)") N-i+1,'*( ',Trim(ADJUSTL(VarCh)),' )',Label
      Else
         write(2,"(I3,A,A,A)") N-i+1,'*( ',Trim(ADJUSTL(VarCh)),' )'
      EndIf
   End Subroutine PrintIntArray
!
   Subroutine PrintIntVal(Var, Var_name, Label)
      Integer, Intent(in) :: Var
      Character(len=*), Intent(in) :: Var_name
      Character(len=*), Intent(in), optional :: Label
      Character(len=20) :: VarCh
      write(VarCh,*) Var
      If(present(Label)) Then
         write(2,"(1x,3A,5x,'! ',A)") TRIM(Var_Name),'= ', Trim(ADJUSTL(VarCh)),Label
      Else
         write(2,"(1x,3A)") TRIM(Var_Name),'= ', Trim(ADJUSTL(VarCh))
      EndIf
      !write(2,"(1x,A,A,A)") TRIM(Var_Name),'= ', Trim(ADJUSTL(VarCh))
   End Subroutine PrintIntVal
!
   Subroutine PrintRealVal(Var, Var_name, Label)
      Real(dp), Intent(in) :: Var
      Character(len=*), Intent(in) :: Var_name
      Character(len=*), Intent(in), optional :: Label
      Character(len=50) :: VarCh
      write(VarCh,*) Var
      If(present(Label)) Then
         write(2,"(1x,3A,5x,'! ',A)") TRIM(Var_Name),'= ', Trim(ADJUSTL(VarCh)),Label
      Else
         write(2,"(1x,3A)") TRIM(Var_Name),'= ', Trim(ADJUSTL(VarCh))
      EndIf
      !write(2,"(1x,A,A,A)") TRIM(Var_Name),'= ', Trim(ADJUSTL(VarCh))
   End Subroutine PrintRealVal
!
   Subroutine PrintLogiVal(Var, Var_name, Label)
      Logical, Intent(in) :: Var
      Character(len=*), Intent(in) :: Var_name
      Character(len=*), Intent(in), optional :: Label
      Character(len=20) :: VarCh
      write(VarCh,*) Var
      If(present(Label)) Then
         write(2,"(1x,3A,5x,'! ',A)") TRIM(Var_Name),'= ', Trim(ADJUSTL(VarCh)),Label
      Else
         write(2,"(1x,3A)") TRIM(Var_Name),'= ', Trim(ADJUSTL(VarCh))
      EndIf
      !write(2,"(1x,A,A,A)") TRIM(Var_Name),'= ', Trim(ADJUSTL(VarCh))
   End Subroutine PrintLogiVal
!
   Subroutine PrintChArray(Var, Var_name, Label)
      Character(len=*), dimension(:), Intent(in) :: Var
      Character(len=*), Intent(in) :: Var_name
      Character(len=*), Intent(in), optional :: Label
      Integer :: N,i
      Character(len=20) :: VarCh
      N=size(Var)
      write(2,"(1x,A,A)", ADVANCE='NO') TRIM(Var_Name),'= '
      Do i=1,N
         write(VarCh,"(A)") Var(i)
         If(Var(i).eq.Var(N)) exit
         write(2,"(A,A,A)", ADVANCE='NO') '"',Trim(VarCh),'", '
      Enddo
      If(present(Label)) Then
         write(2,"(I3,A,A,A,5x,'! ',A)") N-i+1,'*( "',Trim(VarCh),'" )',Label
      Else
         write(2,"(I3,A,A,A)") N-i+1,'*( "',Trim(VarCh),'" )'
      EndIf
   End Subroutine PrintChArray
!
   Subroutine PrintChVal(Var, Var_name, Label)
      Character(len=*), Intent(in) :: Var
      Character(len=*), Intent(in) :: Var_name
      Character(len=*), Intent(in), optional :: Label
      Character(len=200) :: VarCh
      write(VarCh,"(A)") Var
      If(present(Label)) Then
         write(2,"(1x,4A,5x,'! ',A)") TRIM(Var_Name),'= "',Trim(VarCh),'"',Label
      Else
         write(2,"(1x,4A)") TRIM(Var_Name),'= "',Trim(VarCh),'"'
      EndIf
      !write(2,"(1x,4A)") TRIM(Var_Name),'= "',Trim(VarCh),'"'
   End Subroutine PrintChVal
!
   Subroutine PrintParIntro(Label, RunOption, OutFileLabel)
      Character(len=*), Intent(in) :: Label
      Character(len=*), Intent(in) :: RunOption
      Character(len=*), Intent(in) :: OutFileLabel
      write(2,"(20(1x,'='))")
      write(2,*) Label,' run with following useful parameters in namelist: " &Parameters":'
      write(2,"(1x,A)", ADVANCE='NO') '&Parameters '
      Call PrintValues(RunOption,'RunOption')
      Call PrintValues(OutFileLabel,'OutFileLabel')
   End Subroutine PrintParIntro
!
   Subroutine ReadSourceTimeLoc(StartTime_ms, CenLoc)
      use constants, only : dp, Sample, c_mps
      use Chunk_AntInfo, only : TimeBase
      Implicit none
      Real(dp), intent(OUT) :: StartTime_ms, CenLoc(1:3)
      Integer, parameter ::lnameLen=180
      Character(LEN=lnameLen) :: lname
      Real(dp) :: t_shft
      Integer :: j
      Real(dp), external :: tShift_ms
      !
      Call GetNonZeroLine(lname)
      Read(lname(2:lnameLen),*) StartTime_ms, CenLoc  ! Start time
      !Call PrintValues(lname,'input line-1', 'time & position' )
      write(2,"(A)") 'Input line-1: "'//lname(1:1)//'|'//TRIM(lname(2:lnameLen))// &
         '" !  Core/Source-| time, & position'
      Call Convert2m(CenLoc)
      StartTime_ms=StartTime_ms+TimeBase
      t_shft=tShift_ms(CenLoc(:)) ! sqrt(SUM(CenLoc(:)*CenLoc(:)))*1000.*Refrac/c_mps ! in mili seconds due to signal travel distance
      !
      j = iachar(lname(1:1))  ! convert to upper case if not already
      if (j>= iachar("a") .and. j<=iachar("z") ) then
         lname(1:1) = achar(j-32)
      end if
      SELECT CASE (lname(1:1))
         CASE("S")  ! time at Source (central voxel) is given
            StartTime_ms=StartTime_ms+t_shft
         CASE DEFAULT  ! time at reference antenna is given
      End SELECT
      write(2,"(A,F12.6,A,F12.6,A,A,2(F9.4,','),F9.4,A)") &
         ' True start time trace, adding base, at core (at source)=', StartTime_ms, ' (',StartTime_ms-t_shft,') [ms]',&
         ' for source @(N,E,h)=(', CenLoc(1:3)/1000., ' ) km'
      Return
   End Subroutine ReadSourceTimeLoc
   !
END MODULE LOFLI_Input
! ============================================================
Subroutine PreReadPeakFitInfo(ChunkNr_dim,i_peak)
   use constants, only : dp
   use DataConstants, only : RunMode
   use ThisSource, only : PeakNrTotal ! , PeakPos
   use Chunk_AntInfo, only : N_Chunk_max
 !  Use Interferom_Pars, only :  N_fit
   Implicit none
   Integer, intent(out) :: ChunkNr_dim,i_peak
   Integer :: nxx
   Integer :: i_s,i_eo,i_c,i_pos
   Real(dp) :: T1, NEh(1:3)
   CHARACTER(LEN=6) :: txt
   Character(LEN=30) :: FMTSrces
   Integer, parameter :: lnameLen=180
   Character(LEN=lnameLen) ::  lname

   ChunkNr_dim=0
   i_peak=0
   FMTSrces="(1x,3i2,I8,3(F10.2,1x))"
   Call GetNonZeroLine(lname)
   Do  ! perform some pre-scanning of the input to now the number of chunks that will be used
      Read(lname(2:lnameLen),*,iostat=nxx) T1,NEh(:)  ! just dummy arguments
      !write(2,*) nxx, lname
      If(nxx.ne.0) exit  ! Neither a new-chunk line or a pulse-info line
      Read(lname,FMTSrces,iostat=nxx) &
          i_s,i_eo,i_c,i_pos,  NEh(:) ! just dummy arguments
      If(nxx.eq.0 .and. (i_s.ge.0) .and. (i_eo.ge.0 .and. i_eo.lt.3) .and. (i_c.gt.0) .and. (i_pos.gt.1000) ) Then ! check for pulse-info
         !  Determine number of peaks/sorces that are included in the calibration search
         If(.not.(RunMode.eq.7 .and. i_eo.eq.1)) i_peak=i_peak+1
         If((RunMode.eq.2) .and. i_eo.eq.2) i_peak=i_peak+1 ! assume both polarities for Dual
         If(ChunkNr_dim.eq.0) Then
            write(2,*) 'PreReadPeakFitInfo: Chunk-info line should preceed pulse-info lines ',lname
            stop 'Chunk-info line should preceed pulse-info lines'
         EndIf
         Call GetNonZeroLine(lname)
         read(lname,*)  txt  ! check for possible 'exclude' line following this
         If(trim(txt).eq.'exclud') Then
            Call GetNonZeroLine(lname)
         EndIf
      Else
         ChunkNr_dim=ChunkNr_dim+1   ! this was a genuine chunk card
         Call GetNonZeroLine(lname)
      EndIf
      !write(2,*) ChunkNr_dim, i_peak, i_s, i_eo, i_c, i_pos
   EndDo
         !write(2,*) 'test, ChunkNr_dim:',ChunkNr_dim, 'last line:', StartTime_ms,i_dist, i_guess,j,i_chunk, SourceGuess(:,1)
   If(ChunkNr_dim.eq.0) Then
      Write(2,*) lname
      write(2,*) 'ChunkNr_dim:',ChunkNr_dim, 'too small, last line:', lname
      stop 'Chunk number too small'
   EndIf
   If(ChunkNr_dim.gt.N_Chunk_max) Then
      Write(2,*) lname
      write(2,*) 'ChunkNr_dim:',ChunkNr_dim, 'too large, last line:', lname
      stop 'Chunk number too large'
   EndIf
   !
   write(2,*) 'number of calibration chunks:',ChunkNr_dim
End Subroutine PreReadPeakFitInfo
!=================================
Subroutine ReadPeakFitInfo()
!   Read-in the peakpos and source guesses for individual pulses.
!
!
   use constants, only : dp,sample!,Refrac,c_mps
   use Chunk_AntInfo, only : Station_nrMax, RefAnt, Ant_pos, Ant_RawSourceDist
   use Chunk_AntInfo, only : Ant_Stations, Ant_nr, Ant_IDs, StartT_sam, ChunkStTime_ms, ChunkFocus
   use ThisSource, only : SourcePos, RefAntErr, ExclStatNr, Peak_eo, PeakPos, TotPeakNr, ChunkNr !NrP, t_ccorr,
   use ThisSource, only : Dual, PeakNr, PeakNrTotal !, PlotCCPhase, Safety, Nr_corr
   use DataConstants, only : Time_dim, PeakNr_dim, ChunkNr_dim
   use DataConstants, only : RunMode
   use StationMnemonics, only : Statn_ID2Mnem, Station_Mnem2ID
   use FFT, only : RFTransform_su,DAssignFFT
   Implicit none
   integer :: i_eo, i_chunk, i_dist, i, j, k, i_c,i_ca,i_eoa, i_d
   Character(len=5) :: Station_Mnem, ExclStMnem(1:Station_nrMax)
   integer :: i_Peak, nxx
   Integer, parameter :: lnameLen=210
   Character(LEN=lnameLen) ::  lname
   character*10 :: txt
   character*1 :: Label
   character*35 :: Frmt
   Integer :: i_s,i_pos, n_chunk
   Real(dp) :: T1, NEh(1:3)
   Real(dp) :: D,t
    !

   Allocate(ChunkStTime_ms(1:ChunkNr_dim))
   Allocate(ChunkFocus(1:3,1:ChunkNr_dim))
   If(.not. Allocated(ExclStatNr) ) Then
      Allocate( ExclStatNr(1:30,1:PeakNr_dim))
      ExclStatNr(:,:)=0
   EndIf

   !
   i_d=0
   If(Dual) i_d=2
   Frmt="(A1,3i2,I8,3(F10.2,1x),F8.2)"
   Call RFTransform_su(Time_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !
   ExclStatNr(:,:)=0
   TotPeakNr(0,:)=0
   RefAntErr(1:PeakNr_dim) = 0.
   PeakNr(:,:)=0
   i_chunk=0
   n_chunk=0
   i_ca=-1
   i_eoa=-1
   i_peak=0
   t=0.
   !
   Call GetNonZeroLine(lname)
   Do  ! perform some pre-scanning of the input to now the number of chunks that will be used
      Read(lname(2:lnameLen),*,iostat=nxx) T1,NEh(:)  ! just dummy arguments
      !write(2,*)  nxx, lname
      If(nxx.ne.0) exit  ! Neither a new-chunk line or a pulse-info line
      j = iachar(lname(1:1))  ! convert to upper case if not already
      if (j>= iachar("a") .and. j<=iachar("z") ) then
         lname(1:1) = achar(j-32)
      end if
      read(lname,Frmt,iostat=nxx)  Label,i_s,i_eo,i_c, i_pos,  NEh(:) !,t
      !write(2,*) nxx, Label,i_s,i_eo,i_c, i_pos,  NEh(:) ,t
      If(nxx.eq.0 .and. (i_s.ge.0) .and. (i_eo.ge.0 .and. i_eo.le.2) .and. (i_c.gt.0) .and. (i_pos.gt.1000) ) Then ! check for pulse-info
         !  Determine number of peaks/sorces that are included in the calibration search
         If(.not.(RunMode.eq.7 .and. i_eo.eq.1)) Then
            i_peak=i_peak+1 ! always assume both polarities for FieldCalibrate
            If(i_peak .gt. PeakNr_dim) then
               write(2,*) 'PeakNr_dim will be exceeded for #',i_peak
               Write(2,*) 'Expected:   Core/Reference antenna |i |i_eo |i_c |k |N |E |h | dt " in format ',FrMT
               Write(2,*) 'Culprit: "',trim(lname),'"'
               stop 'FindCallibr:ReadPeakInfo problem'
            endif
            If(i_ca.ne.i_c) then
                i_chunk=i_chunk+1
                If(i_eoa.eq.2 .and. (i_chunk.gt.1) .and. (RunMode.eq.2)) Then
                   Call DualReadPeakInfo(i_eoa,i_chunk-1) ! i_eo=2; double the information of the even antennas to the odd ones when Dual=.t. (i_d=2)
                   i_peak=i_peak +PeakNr(0,i_chunk-1)
                EndIf
            ElseIf(i_eoa .ne. i_eo) then
                If(i_eo.eq. 0) i_chunk=i_chunk+1
            Endif
            i_eoa=i_eo
            i_ca=i_c
            If(i_eo.eq.2) then
               i_eo=0 ! same peak info for even&odd
               t=0.  ! assume input came from FCal where t has a different meaning.
            EndIf
            If(i_chunk .gt. n_chunk)  Then
               write(2,*) 'EIReadPeakInfo: Chunk-info line for chunk',i_chunk,' should preceed pulse-info', n_chunk,lname
               stop 'Chunk-info line should preceed pulse-info lines'
            EndIf
            SourcePos(:,i_Peak)=NEh(:)
            !flush(unit=2)
            !write(2,*) 'i_chunk, i_eo,i_peak:',i_chunk, i_eo,i_peak,SourcePos(:,i_Peak)
            !flush(unit=2)
            !write(2,*) 'RefAnt:',RefAnt(i_chunk,i_eo)
            If(Label.eq.'C' .and. RunMode.eq.2) Then
               Call RelDist(SourcePos(:,i_Peak),Ant_pos(:,RefAnt(i_chunk,i_eo),i_chunk),D) !distance to ant - distance to core in samples
               i_pos=i_pos+NINT(D-Ant_RawSourceDist(RefAnt(i_chunk,i_eo),i_chunk)) ! in samples due to signal travel distance to reference
               !write(2,*) i_peak,RefAnt(i_chunk,i_eo),k,x3, Ant_RawSourceDist(RefAnt(i_chunk,i_eo),i_chunk),  &
               !      ' ref station=',Ant_Stations(RefAnt(i_chunk,i_eo),i_chunk), Ant_IDs(RefAnt(i_chunk,i_eo),i_chunk)
               !write(2,*) Label,i,i_eo,i_c,i_ca,i_chunk
            EndIf
            If(Dual .and. (i_eo.eq.1)) then
               Do i=0,PeakNr(0,i_chunk)-1
                  !write(2,*) 'peakpos-i',i,TotPeakNr(0,i_chunk)-i,Peakpos(TotPeakNr(0,i_chunk)-i),k
                  If(Peakpos(TotPeakNr(0,i_chunk)-i).eq. k) then
                     SourcePos(:,i_Peak)=SourcePos(:,TotPeakNr(0,i_chunk)-i)
                     exit
                  Endif
               Enddo
            endif
            Peakpos(i_Peak)=i_pos
            ChunkNr(i_Peak)=i_chunk
            RefAntErr(i_Peak) = 0.
            Peak_eo(i_peak)=i_eo
            PeakNr(i_eo,i_chunk)=PeakNr(i_eo,i_chunk)+1        ! number of peaks for this (i_eo,i_chunk)
            TotPeakNr(i_eo,i_chunk)=i_peak ! last peak# for this (i_eo,i_chunk)
            PeakNrTotal=i_peak
            !write(2,*) 'i_peak,i_eo,i_c',i_peak,i_eo,i_c,i_eoa,i_ca,';',PeakNr
         EndIf
         Call GetNonZeroLine(lname)
         ExclStMnem='     '
         read(lname,* ,iostat=nxx)  txt,ExclStMnem
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
      Else
         n_chunk=n_chunk+1   ! this was a genuine chunk card
         Call Convert2m(NEh(:))
         ChunkStTime_ms(n_chunk)=T1
         ChunkFocus(1:3,n_chunk)=NEh(:)
         write(*,"(A,F7.2,A)") achar(27)//'[45m reading t=',T1,' ms'//achar(27)//'[0m'
         StartT_sam(n_chunk)=(T1/1000.d0)/sample  ! in sample's
         Call AntennaRead(n_chunk,NEh(:))
          !
         If(RunMode.eq.2) Then
            Call Find_unique_StatAnt()
            Call GetRefAnt(n_chunk)
         ElseIf(RunMode.eq.7) Then
          !Ant_ID=MaxLoc(AntMnem, mask=AntMnem.eq.DSet_Names(3))
            Call EISelectAntennas(n_chunk)  ! select antennas for which there is an even and an odd one.
         EndIf

         Call GetNonZeroLine(lname)
         i_ca=-1
      EndIf
   EndDo
   !
   If(i_eoa.eq.0)    Call DualReadPeakInfo(i_d,i_chunk)
   If((RunMode.eq.2) .and. i_eoa.eq.2)    Call DualReadPeakInfo(i_eoa,i_chunk)
   !write(2,*) 'ChunkNr_dim:',ChunkNr_dim
   !write(2,*) 'PeakNr',PeakNr
   !write(2,*) 'PeakNrTotal=',PeakNrTotal
   !write(2,*) 'TotPeakNr',TotPeakNr
   !
   !
   close(unit=14)
   close(unit=12)
   call DAssignFFT()
   !
   !
   Return
End Subroutine ReadPeakFitInfo
!!================================
!================================
Subroutine DualReadPeakInfo(i_eo,i_chunk)
    use constants, only : dp,sample!,Refrac,c_mps
    use ThisSource, only : SourcePos, RefAntErr, ExclStatNr, Peak_eo, TotPeakNr, ChunkNr !NrP, t_ccorr,
    use ThisSource, only : Dual, PeakNr, PeakNrTotal !, PlotCCPhase, Safety, Nr_corr
    use ThisSource, only : PeakPos, Peak_eo
    use DataConstants, only : PeakNr_dim, ChunkNr_dim
    use DataConstants, only : Polariz
    use StationMnemonics, only : Statn_ID2Mnem, Station_Mnem2ID
    Implicit none
    Integer, intent(in) :: i_eo,i_chunk
    Integer :: i
    !write(2,*) 'i_eo:',i_eo,i_chunk,TotPeakNr(0,i_chunk),PeakNr(0,i_chunk)
    !flush(unit=2)
    If(i_eo.eq.2) Then ! same for even&odd
       PeakNr(1,i_chunk)=PeakNr(0,i_chunk)                        ! number of peaks for this (i_eo,i_chunk)
       TotPeakNr(1,i_chunk)=TotPeakNr(0,i_chunk)+PeakNr(1,i_chunk) ! last peak# for this (i_eo,i_chunk)
       Do i=0,PeakNr(0,i_chunk)-1 ! copy
         Peakpos(TotPeakNr(1,i_chunk)-i) = Peakpos(TotPeakNr(0,i_chunk)-i)
         SourcePos(:,TotPeakNr(1,i_chunk)-i) = SourcePos(:,TotPeakNr(0,i_chunk)-i)
         RefAntErr(TotPeakNr(1,i_chunk)-i) = RefAntErr(TotPeakNr(0,i_chunk)-i)
         Peak_eo(TotPeakNr(1,i_chunk)-i)=1
         ChunkNr(TotPeakNr(1,i_chunk)-i)=i_chunk
         PeakNrTotal=TotPeakNr(1,i_chunk)
       EndDo
 !   ElseIf(i_eo.eq.0) Then
 !      PeakNr(1,i_chunk)=0        ! number of peaks for this (i_eo,i_chunk)
  !     TotPeakNr(1,i_chunk)=TotPeakNr(0,i_chunk)+PeakNr(1,i_chunk) ! last peak# for this (i_eo,i_chunk)
    EndIf
    Return
End Subroutine DualReadPeakInfo
!==========================================
Subroutine PrntNewSources()
   use ThisSource, only : PeakNrTotal, ChunkNr, Peak_eo, ChunkNr, Dual, Dropped
   use ThisSource, only : Nr_Corr, CorrAntNrs, ExclStatNr
   use Chunk_AntInfo, only : Ant_Stations, Unique_StatID, Nr_UniqueStat, StartT_sam, ChunkStTime_ms, ChunkFocus
   use Constants, only : Sample
   use DataConstants, only : Polariz, RunMode, Station_nrMax
   use StationMnemonics, only : Statn_ID2Mnem
   Implicit none
   integer ( kind = 4 ) :: i,j, i_Peak, i_chunk, i_eo, i_stat, i_ant, j_corr, count(0:1,1:Station_nrMax), Station_ID
   integer ( kind = 4 ) :: Peaks(0:1)
   integer, external :: XIndx
   Character(len=5) :: Station_Mnem
   Character(len=1) :: FitParam_Mnem(4)=(/'N','E','h','t'/)
   !
   If(RunMode.eq.2 ) Then
      ! Count the number of pulses that are included in the calibration for a station
      Count(:,:)=0
      Peaks(0:1)=0
      Do i_Peak=1,PeakNrTotal
         i_eo=Peak_eo(i_peak)
         i_chunk=ChunkNr(i_peak)
         Peaks(i_eo)=Peaks(i_eo)+1
         !
         Station_ID=0
         Do j_corr=1,Nr_Corr(i_eo,i_chunk)
            j= Station_ID
            Station_ID = Ant_Stations(CorrAntNrs(j_corr,i_eo,i_chunk),i_chunk)
            If(Station_ID .eq. j ) cycle
            Do i_stat=1, Nr_UniqueStat      ! Get station number from the Unique_StatID list
                If(Unique_StatID(i_stat).eq. Station_ID) exit
            enddo
            !write(2,*) 'Station_ID:', Station_ID, i_peak, i_stat,Dropped(i_stat,i_peak)
            If(Dropped(i_stat,i_peak).eq.0) count(i_eo,i_stat)=count(i_eo,i_stat)+1
         EndDo
      EndDo
      write(2,"('Station(',i2,'): ')", ADVANCE='NO') Nr_UniqueStat
      Do i=1,Nr_UniqueStat
         Write(2,"(1x,A5)", ADVANCE='NO') Statn_ID2Mnem(Unique_StatID(i))
      Enddo
      Do i_eo=0,1
         write(2,"(/,'#pulss',i1,'(',I3,') ')", ADVANCE='NO') i_eo,Peaks(i_eo)
         Do i_stat=1,Nr_UniqueStat
            Write(2,"(1x,I5)", ADVANCE='NO') count(i_eo,i_stat)
         Enddo
      Enddo
      write(2,*) ' '
   Else
      Do i_Peak=1,PeakNrTotal
         Do i=1,Station_nrMax
            If(ExclStatNr(i,i_peak).eq.0) exit
            Do i_stat=1, Nr_UniqueStat      ! Get station number from the Unique_StatID list
                If(Unique_StatID(i_stat).eq. ExclStatNr(i,i_peak)) Then
                  Dropped(i_stat,i_peak)=1
                  exit
                EndIf
            Enddo
         Enddo
      EndDo
   EndIf
   !
   !
   Write(2,*) 'Nr,eo,Blk,PPos,(Northing,   Easting,    height, -----); RMS[ns], sqrt(chi^2/df), Excluded: '
      If(Polariz) Then  ! Obsolete?
         write(2,*) 'PrntNewSources: Was not sure if this options is used anywhere, but obviously it is'
         i_chunk=1
         j=0
         Do i_Peak = 1,PeakNrTotal
            If(i_chunk.eq.ChunkNr(i_Peak)) Then
               j=j+1
            Else
               i_eo=1
               Do i=i_Peak-j, i_Peak-1
                  Call PrntCompactSource(i)
               Enddo
               j=1
            EndIf
            i_chunk=ChunkNr(i_Peak)
            i_eo=0
            Call PrntCompactSource(i_Peak)
         Enddo  !  i_Peak = 1,PeakNrTotal
         i_eo=1
         Do i=PeakNrTotal-j+1, PeakNrTotal
            Call PrntCompactSource(i)
         Enddo
      Else
         i_chunk=0
         Do i_Peak = 1,PeakNrTotal
            If(i_chunk.ne.ChunkNr(i_Peak)) then
               i_chunk=ChunkNr(i_Peak)
               Write(2,"(5x,F13.6, 3F10.3,i9)")  ChunkStTime_ms(i_chunk), ChunkFocus(:,i_chunk)/1000, i_chunk
             !  write(2,"(A,i2,A,I11,A,f10.3,A)") 'block=',i_chunk,', starting time=',StartT_sam(i_chunk),&
             !     '[samples]=',StartT_sam(i_chunk)*Sample*1000.d0,'[ms]'
            Endif
            Call PrntCompactSource(i_Peak)
         Enddo  !  i_Peak = 1,PeakNrTotal
      EndIf
   !EndIf
End Subroutine PrntNewSources
!==========================================
Subroutine PrntCompactSource(i_Peak)
   use constants, only : dp,Sample_ms
   use FitParams, only : station_nrMax
   use DataConstants, only : RunMode
   use ThisSource, only : PeakPos, Peak_eo, RefAntErr, PeakNrTotal, ChunkNr, PeakChiSQ, PeakRMS, ExclStatNr, Dropped, SourcePos
   use Chunk_AntInfo, only : Unique_StatID, Ant_pos, Ant_RawSourceDist, RefAnt, StartT_sam
   use StationMnemonics, only : Statn_ID2Mnem
   Implicit none
   Integer, intent(in) :: i_Peak
   integer ( kind = 4 ) :: i,k, i_chunk,i_eo
   Character(len=5) :: Station_Mnem
   Real(dp) :: RDist,Time
   Real(dp), external :: tShift_smpl
   !Character(len=1) :: FitParam_Mnem(4)=(/'N','E','h','t'/)
   !
   i_chunk=ChunkNr(i_Peak)
   If(RunMode.le.3) Then
      i_eo=Peak_eo(i_Peak)
      ! translate peakposition in the reference antenna to one for a (virtual) antenna at the core
      Call RelDist(SourcePos(:,i_Peak),Ant_pos(:,RefAnt(i_chunk,i_eo),i_chunk),RDist) !distance to ant - distance to core in samples
      k=PeakPos(i_Peak)-NINT(RDist-Ant_RawSourceDist(RefAnt(i_chunk,i_eo),i_chunk)) ! in samples due to signal travel distance to reference
   ElseIf(RunMode.ge.7) Then
      i_eo=2
      k=PeakPos(i_Peak)
   EndIf
   Time=( StartT_sam(ChunkNr(i_Peak)) + PeakPos(i_Peak) - tShift_smpl(SourcePos(:,i_Peak)) )*Sample_ms
   !
   If(i_Peak.lt.100) Then
      Write(2,"('C',3i2,I8)", ADVANCE='NO') i_Peak,i_eo,ChunkNr(i_Peak),k
   Else
      Write(2,"('C',i2.2,2i2,I8)", ADVANCE='NO') i_Peak-100,i_eo,ChunkNr(i_Peak),k
   EndIf
   Write(2,"(3(F10.2,','))", ADVANCE='NO') SourcePos(:,i_Peak)
   If(RunMode.le.3) Then
      !Write(2,"(F8.2)", ADVANCE='NO') RefAntErr(i_Peak)
      write(2,"(F12.5';',F7.2,',',F7.2)", ADVANCE='NO') Time, PeakRMS(i_Peak),sqrt(PeakChiSQ(i_Peak))
   ElseIf(RunMode.ge.7) Then
      write(2,"(F12.5';',F7.2)", ADVANCE='NO') Time, PeakChiSQ(i_Peak)
   EndIf
   If(Station_nrMax.gt.0 .and. ((SUM(Dropped(:,i_Peak)).gt.0).or.(ExclStatNr(1,i_peak) .ne.0))) then
      Do i=1,Station_nrMax
         If(ExclStatNr(i,i_peak).eq.0) exit
         Write(2,"(1x,A5)", ADVANCE='NO') Statn_ID2Mnem(ExclStatNr(i,i_peak))
      Enddo
      write(2,"(/,'exclude  ')", ADVANCE='NO')
      Do i=1,Station_nrMax
         If(Dropped(i,i_peak).eq.0) cycle
         Write(2,"(1x,A5)", ADVANCE='NO') Statn_ID2Mnem(Unique_StatID(i))
      Enddo
   Endif
   write(2,*) ' '
   Return
End Subroutine PrntCompactSource
!===========================
