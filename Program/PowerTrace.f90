!-----------------------------------------
 !-------------------------------------------
 !------------------------------------------
    Include 'ConstantsModules.f90'
    Include 'FFT_routines.f90'
    Include 'ParamModules.f90'
    Include 'HDF5_LOFAR_Read.f90'
    Include 'MappingUtilities.f90' ! h5-read
!-----------------------------------
Program PowerTrace
!
!
    use constants, only : dp,sample
    use DataConstants, only : Calibrations
    use DataConstants, only : Time_dim
    use FFT, only : RFTransform_su,DAssignFFT
    use StationMnemonics, only : Station_Mnem2ID, Station_ID2Mnem
    Implicit none
    Character(len=20) :: Utility, release
    INTEGER :: DATE_T(8),i
    !CHARACTER*12 :: REAL_C(3)
    !
    Integer :: j, nxx, ReqStat_ID, ReqStat_eo
    Real*8 :: SourceCoor(3) = (/ 8280.01,  -15120.48,    2618.37 /)     ! 1=North, 2=East, 3=vertical(plumbline)
    Character(len=1) :: EveOdd(0:1)=(/'e','o'/)
    Character(len=5) :: Station_Mnem
    CHARACTER(LEN=6) :: txt
    Character(LEN=100) :: OutDataFileName
    Character(LEN=180) :: lname
    Integer :: UpS
    real(dp) :: StartTime_ms,dTime_ms
    !Character(LEN=180) :: TxtIdentifier, TxtImagingPars
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   release='v14 (April 10, 2020)'
   Utility='PowerTrace'
   !
   OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='PowerTrace.out')
    write(2,"(3x,5(1H-),1x,'PowerTrace version ',A20,25(1H-))") release
    CALL DATE_AND_TIME (Values=DATE_T)
    WRITE(2,"(3X,5(1H-),1x,'run on ',I2,'/',I2,'/',I4,' , started at ',&
          I2,':',I2,':',I2,'.',I3,1X,25(1H-))") &
          DATE_T(3),DATE_T(2),DATE_T(1),(DATE_T(i),i=5,8)
    !
    ! =====================
   StartTime_ms= 1032
   dTime_ms=10.
   UpS=64   ! should be power of 2, not larger 128
   Call GetNonZeroLine(lname)
   Read(lname,*) Calibrations
   write(2,*) 'Calibration file: ',trim(Calibrations)
   Open(unit=12,STATUS='unknown',ACTION='read',FORM ="unformatted", FILE = 'Book/RFI_Filters-v18.uft')
   Open(unit=14,STATUS='old',ACTION='read', FILE = 'Book/LOFAR_H5files_Structure-v18.dat')
1  Continue
   Call GetNonZeroLine(lname)
   Read(lname,*,IOSTAT=nxx) Station_Mnem, StartTime_ms, dTime_ms, UpS, OutDataFileName, ReqStat_eo
   If(nxx.ne.0) goto 9
!Nr_UniqueStat=   1     2     3     4     5     6     7    11    13    17    21    24    26    28    30    31    32   101   103   106
!   125   128   130   141   142   145   146   147   150   161   166   167   169   181   183   188   189
!Nr_UniqueAnt=CS001 CS002 CS003 CS004 CS005 CS006 CS007 CS011 CS013 CS017 CS021 CS024 CS026 CS028 CS030 CS031 CS032 CS101 CS103 RS106
! RS205 RS208 RS210 CS301 CS302 RS305 RS306 RS307 RS310 CS401 RS406 RS407 RS409 CS501 RS503 RS508 RS509
!  RS20n = 12n; RS21n=13n; RS/CS30n=14n; RS31n=15n; RS/CS40n=16n; RS/CS50n=18n
   Call Station_Mnem2ID(Station_Mnem,ReqStat_ID)
   Call Station_ID2Mnem(ReqStat_ID,Station_Mnem)
   Write(2,*) 'ReqStat_ID number=',ReqStat_ID,'=',Station_Mnem,', i_eo=',ReqStat_eo
   If((ReqStat_eo.ne.0) .and. (ReqStat_eo.ne.1)) ReqStat_eo=0
   write(2,*) 'start & stop times=',StartTime_ms,StartTime_ms+dTime_ms, UpS,'[samples] '!,ReqStat_eo
   Call GetNonZeroLine(lname)
   Read(lname,*) SourceCoor
   OutDataFileName=trim(OutDataFileName)//Station_Mnem//EveOdd(ReqStat_eo)//'.dat'
   write(*,"(A,A)") 'file=',trim(OutDataFileName)
   write(2,*) ' SourceCoor=', SourceCoor, trim(OutDataFileName)
   !StartTime_ms=784.
   OPEN(UNIT=29,STATUS='unknown',ACTION='WRITE',FILE=trim(OutDataFileName))
   Call RFTransform_su(Time_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   Call MkPwrTrace( SourceCoor,StartTime_ms,dTime_ms,UpS,ReqStat_ID, ReqStat_eo)
   Call DAssignFFT()
   write(2,*) 'out data file: ',trim(OutDataFileName)
   Close(unit=29)
   goto 1
9  Continue
   write(2,*) 'Last input line=',lname
   Close(unit=12)
   Close(unit=14)
   Stop 'Normal end'
end program PowerTrace
!=================================
Subroutine MkPwrTrace(SourceCoor,StartTime_ms,dTime_ms,UpS,ReqStat_ID, ReqStat_eo)
   use constants, only : dp,sample, Refrac, c_mps  !pi,ci
   use DataConstants, only : Time_dim, Cnu_dim  ! , ChunkNr_dim
   use HDF5_LOFAR_Read, only : filename, Group_names, Group_nr ! , Group_max
   use HDF5_LOFAR_Read, only : DSet_names, DSet_nr, Ant_ID, STATION_ID !, DSet_max, ANTENNA_POSITION
   use HDF5_LOFAR_Read, only : DATA_LENGTH, SAMPLE_NUMBER_first, Absolute_TIME, DIPOLE_CALIBRATION_DELAY
   use HDF5_LOFAR_Read, only : GetDataChunk
   use FFT, only : RFTransform_CF, RFTransform_CF2CT, RFTransform_CF_Filt,Hann
   Use Calibration, only : Station_ID2Calib
   Implicit none
   Integer, intent(inout) :: UpS ! Upsample number
   Integer, intent(inout) :: ReqStat_ID, ReqStat_eo ! Station ID for which power spectrum should be made
   Real(dp), intent(in) :: SourceCoor(3),StartTime_ms,dTime_ms ! = (/ 8280.01,  -15120.48,    2618.37 /)     ! 1=North, 2=East, 3=vertical(plumbline)
   !
   Integer :: i_chunk, OverLap, i_chunkMax, Calibrated
   Integer*2 :: Chunk(1:Time_dim)
   Real(dp) :: RTime_s(1:Time_dim)
   Complex(dp) :: CNu_s(0:Cnu_dim), CTime_s(1:Time_dim)
   integer :: Dset_offset, Sample_Offset, DataReadErr
   Real(dp) :: nu_Fltr(0:Cnu_dim) !, Av, Bv, FiltFact
   Real*8 :: Powr, p, t, x, DistSRC
   !
   Integer :: i,i_file,i_grp,i_dst, ir_file,ir_grp,ir_dst, i_sampl, i_samplM
   Integer :: i_ant, i_time, i_eo, unt, nxx, StAntID, sgn
   Real*8 :: SubSample_Offset, LFRAnt_crdnts(3), RDist, T_Offset,StatAnt_Calib
   !
   Rewind(unit=14)
   Rewind(unit=12)
   ir_file=0 ; ir_grp=0 ; ir_dst=0
   !           call flush(2)
   Do i_file=1,130 ! 10 !80
      read(14,*,IOSTAT=nxx)  Group_nr, filename
      If(nxx.ne.0) exit
      write(*,"(A,i3)", ADVANCE='NO') achar(27)//'[100D'//'file#=',i_file  ! [1000D    !  //achar(27)//'[0m.'
      Do i_grp=1,Group_nr
         read(14,*) DSet_nr, Group_Names(i_grp)
         Do i_dst=1,DSet_nr
              read(14,*) DSet_Names(i_dst), STATION_ID, Ant_ID
              If(STATION_ID .ne. ReqStat_ID) cycle
              i_eo=Mod(Ant_ID,2)
              If(i_eo.ne.ReqStat_eo) cycle
              Write(2,*) DSet_nr, Group_Names(i_grp), DSet_Names(i_dst), STATION_ID, Ant_ID
             Do while((i_file .gt. ir_file) .or. (i_grp .gt. ir_grp) .or. (i_dst .gt. ir_dst) )
                  !Write(2,*) 'reading Filter data',i_file,i_grp,i_dst,ir_file,ir_grp,ir_dst
                  Read (12) ir_file,ir_grp,ir_dst, LFRAnt_crdnts, powr,&
                  DATA_LENGTH, SAMPLE_NUMBER_first, Absolute_TIME, DIPOLE_CALIBRATION_DELAY,nu_fltr !nu_fltr
              enddo
              If (Powr.le.1.d-6) then
                  Write(2,*) 'reading Filter data',i_file,i_grp,i_dst,ir_file,ir_grp,ir_dst
                  write(2,*) 'powr:', powr,&
                     DATA_LENGTH, SAMPLE_NUMBER_first, Absolute_TIME, DIPOLE_CALIBRATION_DELAY!,nu_fltr
                  cycle
              EndIf
              If((i_file .ne. ir_file) .or. (i_grp .ne. ir_grp) .or. (i_dst .ne. ir_dst) ) then
                  Write(2,*) '****Filter data not found',i_file,i_grp,i_dst,ir_file,ir_grp,ir_dst
                  cycle   ! get next DSet
              Endif
              !
              Call Station_ID2Calib(STATION_ID,Ant_ID,StatAnt_Calib, Calibrated) ! StatAnt_Calib in units of samples
              !
              Write(2,"(A,I4,I3,A,3F9.1)") 'Antenna coordinates of ', STATION_ID, Ant_ID,' (N,E,h)',LFRAnt_crdnts
              DistSRC=sqrt(sum(SourceCoor(:)*SourceCoor(:)))
              RDist=sqrt(sum((SourceCoor(:)-LFRAnt_crdnts(:))**2))
              write(2,*) 'distances [m]',DistSRC,RDist
              Call RelDist(SourceCoor,LFRAnt_crdnts,RDist)
      !         -DistMax*Refrac/c_mps ! convert to s
              !Call Station_ID2Calib(STATION_ID,Ant_ID,StatAnt_Calib) ! StatAnt_Calib in units of samples
              !write(2,*) 'STATION_ID,Ant_ID,StatAnt_Calib',STATION_ID,Ant_ID,StatAnt_Calib
              !
      ! antRead: Dset_offset=Start_time(i_chunk) + RDist + DIPOLE_CALIBRATION_DELAY/Sample + StatAnt_Calib - SAMPLE_NUMBER_first
              Sample_Offset= INT((StartTime_ms/1000. + DIPOLE_CALIBRATION_DELAY)/Sample)+ StatAnt_Calib- SAMPLE_NUMBER_first   ! in units of samples
              Sample_Offset= Sample_Offset+Rdist + DistSRC*Refrac/(c_mps*Sample)  ! Changed to -DistSRC March 10,2022
              write(2,*) 'Sample_Offset', Sample_Offset, Rdist,  DistSRC*Refrac/(c_mps*Sample)
              SubSample_Offset = 0
              OverLap=128
              If(Mod(UpS,2).ne.0) UpS=UpS+1
              i_samplM=INT((Time_dim-2*OverLap)/UpS)  ! number of sections in a block of data
              If(i_samplM.le.0) stop 'i_samplM'
              Overlap=(Time_dim-i_samplM*UpS)/2 ! Recalculate the overlap to have it match the re-sample size
              i_samplM = i_samplM -1  ! to later start counting from zero
              write(2,*) 'i_samplM=',i_samplM,UpS, OverLap, i_samplM+1- (Time_dim-2*OverLap)/UpS ! last should be =0
              i_chunkMax=INT(dTime_ms/(1000.*sample*(Time_dim-2*OverLap)))
              write(2,*) 'i_chunkMax=',i_chunkMax, Sample_Offset, Time_dim, OverLap
              call flush(2)
              Do i_chunk=0,i_chunkMax
                 Dset_offset= Sample_Offset +i_chunk*(Time_dim-2*OverLap+1)
                 !If(Ant_nr(1).eq.0) write(2,*) ', TimeOffset', T_Offset, Dset_Offset,DATA_LENGTH, Sample_Offset - SAMPLE_NUMBER_first
                 If(Dset_offset .gt.DATA_LENGTH) then
                     Write(2,*) 'run out of data'
                     stop 'run out of data'
                 endif
                 !
                 ! Retrieve a chunk of data
                 DataReadErr=0
                 Call GetDataChunk(Group_Names(i_grp),DSet_Names(i_dst), Chunk, Dset_offset, Time_dim, DataReadErr=DataReadErr)
                 If(DataReadErr.ne.0) then
                  write(2,*) '************ error in data-read!, file #=',i_file,', for chunk#',i_chunk,', name=',filename
                  write(*,*) 'error at chunk#',i_chunk
                  stop 'DataReadError'
                 endif
                 !
                 !write(2,*) 'powr',Powr
                 !write(2,*) 'Chunk(1:10):',Chunk(1:10)
                 !write(2,*) 'Hann(:):',Hann(1:10)
                 RTime_s(:)=Chunk(:)*Hann(:)/sqrt(Powr/7.)  !  7=2 pi ???
                 Call RFTransform_CF_Filt(RTime_s,nu_fltr,SubSample_Offset,Cnu_s)
                 !write(2,*) 'RTime_s(1:10):',RTime_s(1:10)
                 Call RFTransform_CF2CT(Cnu_s,CTime_s )
                 !Write(2,*) 'i_samplM=', Dset_offset+OverLap, Dset_offset+i_samplM*UpS+UpS+OverLap,p
                 Do i_sampl=0,i_samplM
                     p=0.
                     Do i=1,UpS
                        !If( (i_sampl*UpS+i+OverLap).gt. Time_dim) write(*,*) (i_sampl*UpS+i+OverLap),i_samplM
                        x=ABS(CTime_s(i_sampl*UpS+i+OverLap))
                        p=p+x*x
                        !write(2,*) 'px',p,x,i_sampl*UpS+i+OverLap,CTime_s(i_sampl*UpS+i+OverLap)
                     Enddo
                     p=p/UpS
                     t=StartTime_ms + (i_chunk*(Time_dim-2*OverLap) + OverLap + (i_sampl+0.5)*UpS)*sample*1000.
                     write(29,"(f11.6,g14.3)") t,p
                     !write(2,*) 't,p,UpS',t,p,UpS
                     !stop
                 Enddo
              Enddo  ! i_chunk
              Write(2,*) 'power data written to file'
              goto 9
              !
          1   continue
          EndDo  ! i_dst=1,DSet_nr
          !
      Enddo ! i_grp=1,Group_nr
      !
   Enddo ! i_file=1,....
   !
9  continue
   write(*,*) ' '
   !
   Return
   !
End Subroutine MkPwrTrace
! ========================
