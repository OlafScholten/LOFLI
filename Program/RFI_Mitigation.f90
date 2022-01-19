!-------------------------------------
    Include 'ConstantsModules.f90'
    Include 'FFT_routines.f90'
    Include 'ParamModules.f90'
    Include 'HDF5_LOFAR_Read.f90'
    Include 'MappingUtilities.f90'
    Include 'System_Utilities.f90'
    Include 'GLEplotUtil.f90'
!    Include 'ConstantsModules-v15.f90'
!    Include 'ParamModules-v15.f90'
!    Include 'MappingUtilities-v15.f90'
!-----------------------------------------
Module RFI_MitPars
    use constants, only : dp
    use DataConstants, only : DataFolder, Time_dim, Cnu_dim, Ant_nrMax
    Implicit none
    integer :: NBackgr ! Number of chunck used for frequency filter determination
    Integer :: Filtring ! number of filtered frequencies
    Integer :: FiltringRef=-1 ! set to the number for the first antenna, to put semi-edicated limits on the range
    Integer ::  RefOdd=0, RefEve=0 ! Antenna-ID of reference antennas
    Integer, allocatable :: MaxAmpChunk(:)
    Integer :: i_chunkMax
    Integer :: MinAmp ! Amplitude for determining background
    Integer :: N_zeroChunk  ! number of chunks with too many zeros
    Real(dp) :: Powr
    Real(dp) :: nu_Fltr(0:Cnu_dim) ! RFI filter
    Real(dp) :: Freq_s(0:Cnu_dim) ! Average frequency spectrum for background
    Logical :: MakePlots
End Module RFI_MitPars
! ---------------------------------------------
program RFI_mitigation
! v18: include search for part of spectrum with no lightning activity
    use constants, only : dp,pi,ci, sample
    use HDF5_LOFAR_Read, only : filename,Group_names, Group_nr, Group_max
    use HDF5_LOFAR_Read, only : DSet_names, DSet_nr, DSet_max, Ant_ID, STATION_ID
    use HDF5_LOFAR_Read, only : ANTENNA_POSITION, DATA_LENGTH, SAMPLE_NUMBER_first, Absolute_TIME, DIPOLE_CALIBRATION_DELAY
    use HDF5_LOFAR_Read, only : GetFileName, ListGroups, ListDataAtt, ListGroupStructure, CloseFile, CloseGroup, CloseDSet!, GetData
    use FFT, only : RFTransform_su,DAssignFFT !, RFTransform_CF, RFTransform_CF2CT, Hann, RFTransform_CF_Filt
    use Chunk_AntInfo, only : BadAnt_SAI, BadAnt_nr_max, Unique_SAI
    use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, DataFolder, FlashName, Windows
   use DataConstants, only : RunMode ! 1=Explore; 2=Calibrate; 3=Impulsive imager; 4=Interferometry
    Use RFI_MitPars
    use GLEplots, only : GLEplotControl
    Implicit none
    Character(len=20) :: Utility, release !, Char_AbsTimeOffset
    Character(len=5) :: Version='v18  '
    INTEGER :: i, VeryFirstSampl, MostLatestSampl=0
    Character(len=5) :: txt
    !
    Integer :: nxx, NFail=0 ! number of antennas for which RFI-determination failed
    Integer :: i_even=0, i_odd=0, MinAmp_Ant(Ant_nrMax), Powr_Ant(Ant_nrMax)
    Real(dp) :: LFRAnt_crdnts(3), MinAmpAv, MinAmpSq, sigma, powrAv, powrSq, powrSigm, power_Even, power_Odd
    !
    Integer :: j,i_fil,i_filR,i_file,i_grp,i_dst, i_chunk, i_ant !
    Character(LEN=100) :: AntennaFieldsDir='../AntennaFields/'   ! Directory where info on antenna positions
      !           can be found if not present in the HDF5 files of the dat
    !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   i_filR= 2  ! The file containing the reference antenna
   j=IARGC()
   If(j.gt. 0) then
      CALL getarg(1, AntennaFieldsDir)
      write(*,*) 'argument 1:"',TRIM(AntennaFieldsDir),'"'
      !Read(SystemCommand,*) AntennaFieldsDir
   Endif
   If(j.gt. 1) then
      CALL getarg(2, txt)
      write(*,*) 'argument 2:"',TRIM(txt),'"'
      Read(txt,*) i_filR
   Endif
   !CALL get_environment_variable("AntennaFieldsDir", AntennaFieldsDir)
   !write(*,*) 'AntennaFieldsDir:"',TRIM(AntennaFieldsDir),'"'
   !
   OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='RFI_mit.out')
   Utility='RFI-Mitigation'
   release='v19, April 2021'   ! TRIM(Version)//'; Nov, 2020'
   RunMode=0
   !write(2,"(3x,5(1H-),1x,'RFI_mitigation release of ',A22,25(1H-))") release
   !CALL DATE_AND_TIME (REAL_C(1),REAL_C(2),REAL_C(3),DATE_T)
   !WRITE(2,"(3X,5(1H-),1x,'run on ',I2,'/',I2,'/',I4,' , started at ',&
    !   I2,':',I2,':',I2,'.',I3,1X,25(1H-))") &
   !    DATE_T(3),DATE_T(2),DATE_T(1),(DATE_T(i),i=5,8)
   Call System_Initiation(Utility, release, ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows)  ! set Flash name & folder names
   !
   VeryFirstSampl=5./sample ! well beyond possible upper limit of 5 seconds
   i_ant=0
   MinAmpAv=0.
   MinAmpSq=0.
   powrAv=0.
   powrSq=0.
   power_Even=0.
   power_Odd=0.
   !
   Write(2,*) 'Antenna information from files in folder: ',TRIM(AntennaFieldsDir)
   !Call GETCWD(FlashFolder)  ! Get current working directory
   write(2,*) 'Reference antenna from file#',i_filR,' in the list at: /directory.out' ! trim(DataFolder)
   flush(unit=2)
   !
   Open(unit=12,STATUS='unknown',ACTION='WRITE',FORM ="unformatted", &
      FILE = 'Book/RFI_Filters-'//TRIM(Version)//'.uft')
   Open(unit=14,STATUS='unknown',ACTION='WRITE', &
      FILE = 'Book/LOFAR_H5files_Structure-'//TRIM(Version)//'.dat')
   Open(UNIT=16,STATUS='unknown',ACTION='WRITE',  FILE=trim(DataFolder)//'RFI_EvenStat.dat')
   Open(UNIT=17,STATUS='unknown',ACTION='WRITE',  FILE=trim(DataFolder)//'RFI_OddStat.dat')
   Open(UNIT=18,STATUS='unknown',ACTION='WRITE',  FILE=trim(DataFolder)//'RFI_EvenZero.dat')
   Open(UNIT=19,STATUS='unknown',ACTION='WRITE',  FILE=trim(DataFolder)//'RFI_OddZero.dat')
   !
    call RFTransform_su(Time_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Do i_fil=1, 140
      If(i_fil.eq.1) Then ! swap numbers of the first and the reference antenna
         i_file=i_filR
      Elseif(i_fil.eq.i_filR) Then
         i_file=1
      Else
         i_file=i_fil
      Endif
      Call GetFileName(i_file, nxx)
      If(nxx .ne. 0) exit
      write(*,"(A,I0,': ',A)") 'i_file=',i_file,trim(filename)
      !write(14,"(I3,A)") i_file, trim(filename)
      Call ListGroups ! >>>>>>>>>>>>>>>>>>>> opens HDF5 file & extracts group-structure
      write(14,"(i3,2x,3A)")  Group_nr, '"',trim(filename),'"'
      Do i_grp=1,Group_nr
         Write(2,"(A,i0,A,A)") 'Group#',i_grp,' :',trim(Group_names(i_grp))
         Call ListGroupStructure(Group_Names(i_grp)) ! >>>>>>>>>>>>>>>>>>>> opens HDF5 Group in the file
         write(14,"(i3,1x,A)") DSet_nr, trim(Group_Names(i_grp))
         Do i_dst=1,DSet_nr
            Call ListDataAtt(Group_Names(i_grp),DSet_Names(i_dst)) ! >>>>>>>>>>>>>>>>>>>>
            write(14,"(1x,A,I6,I6)")  trim(DSet_Names(i_dst)), STATION_ID, Ant_ID
            !
            If(VeryFirstSampl.gt.SAMPLE_NUMBER_first) VeryFirstSampl=SAMPLE_NUMBER_first
            !
            !Scan file for non-lightning sections
            ! Note that the noise is NOT gaussian since RFI has not been removed yet !
            !DATA_LENGTH=1./sample  ! to shorten time to 1 second for testing
            Powr=-1.
            nu_fltr=0.
            i_chunkMax=DATA_LENGTH/Time_dim -1
            If(i_chunkMax.lt.100) goto 999 ! not worthwhile
            If(((RefOdd.eq.0) .and. (mod(Ant_ID,2).eq.1)) .or. ((RefEve.eq.0) .and. (mod(Ant_ID,2).eq.0))) then
               MakePlots=.true.
            Else
               MakePlots=.false.
            EndIf
            !MakePlots=.true.
            Allocate( MaxAmpChunk(0:i_chunkMax))
            Call RoughStatistics(i_dst)
            !
            Call AccumulateBackgrFreq()
            DeAllocate (MaxAmpChunk)
            !
            ! Determine averaged frequency spectrum for low-MaxAmpli chuncks
            Call BuildRFIFilter()
            !
            !
            i_ant=i_ant+1
            If(i_ant.le. Ant_nrMax) then
               MinAmp_Ant(i_ant)=MinAmp
               Powr_Ant(i_ant)=NINT(powr)
               Unique_SAI(i_ant)=STATION_ID*1000+ Ant_ID
            EndIf
            If(powr.le.0.5) then
               Nfail=Nfail+1
               If(Nfail.lt.BadAnt_nr_max) BadAnt_SAI(Nfail)=STATION_ID*1000+ Ant_ID
               goto 999
            Else
               MinAmpAv=MinAmpAv + MinAmp
               MinAmpSq=MinAmpSq + MinAmp*MinAmp
               powrAv=powrAv + powr
               powrSq=powrSq + powr *powr
            Endif
            !
            If(VeryFirstSampl.gt.SAMPLE_NUMBER_first) VeryFirstSampl=SAMPLE_NUMBER_first
            If(MostLatestSampl.lt.(SAMPLE_NUMBER_first+DATA_LENGTH)) MostLatestSampl=(SAMPLE_NUMBER_first+DATA_LENGTH)
            !
            ! Test effect filter
            If(i_file.eq.i_filR) then  ! The reference station
               If((RefOdd.eq.0) .and. (mod(Ant_ID,2).eq.1)) then
                  RefOdd=Ant_ID
                  write(2,*) 'odd reference antenna',trim(DSet_Names(i_dst))
                  Call PlotFreqFilt(i_dst)  ! Frequency spectra
                  Call TestFilter(i_dst)     ! Filtered time spectra
               Else If((RefEve.eq.0) .and. (mod(Ant_ID,2).eq.0)) then
                  RefEve=Ant_ID
                   write(2,*) 'even reference antenna',trim(DSet_Names(i_dst))
                  Call PlotFreqFilt(i_dst)  ! Frequency spectra
                  Call TestFilter(i_dst)     ! filtered spectra
               EndIf
            EndIf
            !
999         continue
            !
            !DeAllocate (MaxAmpChunk)
            If(mod(Ant_ID,2).eq.0) then
               i_even=i_even+1
               If(powr.gt.0) then
                  Write(16,"(i4,1x,A10,1x,f7.2,1x,I5,1x,I4,1x,F6.2)") i_even, trim(DSet_Names(i_dst)),&
                              sqrt(powr), MinAmp, Filtring, 100.*N_zeroChunk/i_chunkMax
                  power_Even=power_Even +powr
               Endif
                  Write(18,"(i4,1x,A10,1x,F6.2,1x,F6.2)") i_even, trim(DSet_Names(i_dst)), 100.*N_zeroChunk/i_chunkMax, powr
                  Write(2,"(A,i4,A,I3,A,A10,A,I6,F6.2,A,F6.2,A,I4)") &
                     'i_ant=',i_ant,', i_even=', i_even,', antenna:', trim(DSet_Names(i_dst))&
                     ,'; # of zero-chunks=', N_zeroChunk,100.*N_zeroChunk/i_chunkMax,'%; power norm=', powr, &
                     '; # filtered frequencies:', Filtring
            Else
               i_odd=i_odd+1
               If(powr.gt.0) then
                  Write(17,"(i4,1x,A10,1x,f7.2,1x,I5,1x,I4,1x,F6.2)") i_odd, trim(DSet_Names(i_dst)), &
                              sqrt(powr), MinAmp, Filtring, 100.*N_zeroChunk/i_chunkMax
                  power_Odd= power_Odd+powr
               Endif
                  Write(19,"(i4,1x,A10,1x,F6.2,1x,F6.2)") i_odd, trim(DSet_Names(i_dst)), 100.*N_zeroChunk/i_chunkMax, powr
                  Write(2,"(A,i4,A,I3,A,A10,A,I6,F6.2,A,F6.2,A,I4)") &
                     'i_ant=',i_ant,', i_odd=', i_odd,', antenna:', trim(DSet_Names(i_dst))&
                     ,'; # of zero-chunks=', N_zeroChunk,100.*N_zeroChunk/i_chunkMax,'%; power norm=', powr, &
                     '; # filtered frequencies:', Filtring
            EndIf
            !
            Call GetAntLofarAntPos(LFRAnt_crdnts, AntennaFieldsDir)
            Write(2,*) 'LOFAR coordinates=',LFRAnt_crdnts
            Write(12) i_fil,i_grp,i_dst,LFRAnt_crdnts, powr, &
              DATA_LENGTH, SAMPLE_NUMBER_first, Absolute_TIME, DIPOLE_CALIBRATION_DELAY,nu_fltr
            Call CloseDSet !<<<<<<<<<<<<<<<<<<<<
         Enddo ! i_dst=1,DSet_nr
         Call CloseGroup !<<<<<<<<<<<<<<<<<<<<
        EndDo ! ListGroups
        Call CloseFile !<<<<<<<<<<<<<<<<<<<<
   Enddo
998 Continue
   Close(unit=14)
   Close(unit=15)
   Close(unit=16)
   Close(unit=17)
   Close(unit=18)
   Close(unit=19)
   Close(unit=12)
   call DAssignFFT()
   write(2,*) Power_Even/i_even, power_Odd/i_odd, Power_Even/power_Odd
   !
   Call GLEplotControl(PlotType='RFI_Plots', PlotName='RFI_Plots'//TRIM(FlashName), &
            PlotDataFile=TRIM(DataFolder), Submit=.true.)
   !
   !Write(2,*) 'total nr antennas and stations'
   !Write(2,*) 'Check for bad antennas occurring twice; indicate for bad antennas why they are bad'
   write(2,*) 'time range for which there are data [ms]',VeryFirstSampl*sample*1000., MostLatestSampl*sample*1000.
   Write(2,*) '!!!!!!!!! For ',NFail,' Antennas RFI-mitigation failed, out of a total of',i_ant,' !!!!!!!!!!!!'
   !write(2,*) MinAmp_Ant(1:i_ant)
   MinAmpAv=MinAmpAv/(i_ant-Nfail)
   MinAmpSq=MinAmpSq/(i_ant-Nfail)
   sigma=sqrt(MinAmpSq-MinAmpAv*MinAmpAv)
   powrAv=powrAv/(i_ant-Nfail)
   powrSq=powrSq/(i_ant-Nfail)
   powrSigm=sqrt(powrSq-powrAv*powrAv)
   If(Nfail.gt.BadAnt_nr_max) Then
      write(2,*) '!!! Warning, not all bad antennas are listed:'
      Nfail=BadAnt_nr_max
   EndIf
   !write(2,*) BadAnt_SAI(1:Nfail)
   write(2,"(A,I7,9(',',I7),10(/12x,I7,9(',',I7)) )") ' BadAnt_SAI=', BadAnt_SAI(1:Nfail)
   write(2,*) 'MinAmp statistics:', MinAmpAv,', sigma=', sigma!, sigma, i_ant-Nfail
   write(2,*) 'statistics powr=sq(raw/Norm):', powrAv,', sigma=', powrSigm!, powrsigm, i_ant-Nfail
   If(i_ant.gt. Ant_nrMax) then
      Write(2,*) '!!!!!!!!!! Nr of antennas exceeds max:',I_ant
      i_ant=Ant_nrMax
   EndIf
   Nfail=0
   Do i=1,i_ant
      !write(2,*) i,Unique_SAI(i),MinAmp_Ant(i),powr_Ant(i)
      If(powr_Ant(i).le. 0) cycle
!      write(2,*) i,Unique_SAI(i),MinAmp_Ant(i),powr_Ant(i)
      If( abs(MinAmp_Ant(i)-MinAmpAv) .gt. 2*sigma) then
         write(2,*) '!!!!!!!! Bad antenna',Unique_SAI(i), &
            ' based on Raw Amplitude deviation of', (MinAmp_Ant(i)-MinAmpAv),' exceedig > 2sigma', &
            ', antenna powr dev=', (powr_Ant(i)-powrAv)
         If(Nfail.lt.BadAnt_nr_max) then
            Nfail=Nfail+1
            BadAnt_SAI(Nfail)=Unique_SAI(i)
         Endif
      Else If( abs(powr_Ant(i)-powrAv) .gt. 2*powrsigm) then
         write(2,*) '!!!!!!!! Bad antenna',Unique_SAI(i), &
            ' based on powr deviation of',(powr_Ant(i)-powrAv)
         If(Nfail.lt.BadAnt_nr_max) then
            Nfail=Nfail+1
            BadAnt_SAI(Nfail)=Unique_SAI(i)
         Endif
      EndIf
   Enddo
   write(2,"(A,I7,9(',',I7),10(/12x,I7,9(',',I7)) )") ' BadAnt_SAI=', BadAnt_SAI(1:Nfail)
   !
   !CALL SYSTEM(SystemCommand,nxx)
   !write(2,*) nxx,'=zero for success '
   stop
   !
end program RFI_mitigation
!=================================
Subroutine GetAntLofarAntPos(LFRAnt_crdnts, AntennaFieldsDir)
! Base files from https://github.com/lofar-astron/lofar-antenna-positions
   use constants, only : dp, pi
   use HDF5_LOFAR_Read, only : Ant_ID, STATION_ID
   use StationMnemonics, only : Station_ID2Mnem
   Implicit none
   Real(dp), intent(out) :: LFRAnt_crdnts(1:3)
   Character(Len=100), intent(in) :: AntennaFieldsDir
   Character(len=5) :: Station_Mnem
   Character(len=5) :: txt
   Character(len=3) :: LBA
   Logical, save :: First=.true.
   Real(dp), save :: ETRS_RotMat(3,3),Center_CS002(3)
   real(dp) :: Rotmatphi(3,3),Rotmatth(3,3)
   real(dp) :: x,theta, phi
   Real(dp) :: ETRS(1:3)
   INTEGER :: i,j,k,nxx
   Real(dp) :: central(3),RelPos(1:3,0:1)
   !
   Open(unit=11, STATUS='old',ACTION='read', FILE = TRIM(AntennaFieldsDir)//'etrs-antenna-positions.csv')
   !
   If(First) Then
      !    real(dp), intent(out) :: RotMat(3,3) = reshape( (/ 0.8056597 ,  0.        ,  0.59237863 , &
      !                              -0.05134149, -0.99623707,  0.06982658 , &
      !                              -0.59014956,  0.08667007,  0.80262806 /), (/ 3,3/))
      !# Center of station CS002 in ITRF coordinates, substract this before rotating
      !    real(dp) :: center_CS002(3)= (/ 3826577.5, 461021.3125, 5064893.0 /)
      theta=(52.915141-90.d0)*pi/180.
      phi=6.8698134123520900*pi/180. ! corrects for omegaxomegaxr-->0.19878 deg
      Rotmatphi=0.
      RotMatphi(1,1)=cos(phi)  ; RotMatphi(1,2)=sin(phi)
      RotMatphi(2,1)=-sin(phi) ; RotMatphi(2,2)=cos(phi)
      RotMatphi(3,3)=1
      RotMatTh=0.
      RotMatTh(2,2)=1.
      RotMatTh(1,1)=-cos(theta)  ; RotMatTh(1,3)=-sin(theta)
      RotMatTh(3,1)=-sin(theta) ; RotMatTh(3,3)=cos(theta)
      Do i=1,3
        Do j=1,3
            x=0.
            Do k=1,3
                x=x+RotMatTh(i,k)*RotMatPhi(k,j)
            enddo
            ETRS_RotMat(i,j)=x
        enddo
      enddo
      !    write(2,*) 'phi=',b,(center_CS002(2)/b),asin(center_CS002(2)/b) *180./pi
      !    write(2,*) 'theta=',a,(center_CS002(3)/a),asin(center_CS002(3)/a) *180./pi,asin(b/a)*180./pi
      !    write(2,*) 'RotMat=',RotMat
      txt=' '
      Read(11,*,IOSTAT=nxx) txt
      Read(11,*,IOSTAT=nxx) txt
      Do while(txt.ne.'CS002')
        Read(11,*,IOSTAT=nxx) txt,LBA,j,center_CS002
        If(nxx.ne.0) stop 'GetAntLofarAntPos'
      enddo
      Rewind(unit=11)
      First=.false.
   EndIf  ! First
   !
   Call Station_ID2Mnem(STATION_ID,Station_Mnem)
   txt=' '
   Do while(txt.ne.Station_Mnem)
     Read(11,*,IOSTAT=nxx) txt
     If(nxx.ne.0) stop 'GetAntLofarAntPos'
   enddo
   BackSpace(unit=11)
   !
   k=49+Ant_ID/2
   Do i =1,k
     Read(11,*) txt,LBA,j,ETRS
   enddo
   Close(unit=11)
   !Write(2,*) txt,Ant_ID,k,LBA,j,ETRS
   !
   Do i=1,3
     x=0.
     Do j=1,3
         x=x + ETRS_RotMat(i,j) * (ETRS(j)-center_CS002(j))
     enddo
     LFRAnt_crdnts(i)=x  ! 1=North, 2=East, 3=vertical(plumbline)
   enddo
   !
   Return
End Subroutine GetAntLofarAntPos
!------------------
Subroutine RoughStatistics(i_dst)
!Scan file for non-lightning sections
! Note that the noise is NOT gaussian since RFI has not been removed yet !
   use constants, only : dp,pi,ci, sample
   use DataConstants, only : DataFolder, Time_dim
   use HDF5_LOFAR_Read, only : DSet_names!, DATA_LENGTH
   use HDF5_LOFAR_Read, only : SAMPLE_NUMBER_first, DIPOLE_CALIBRATION_DELAY
   use HDF5_LOFAR_Read, only : GetData
   Use RFI_MitPars, only : i_chunkMax, MaxAmpChunk, MakePlots, N_zeroChunk
   Implicit none
   Integer, intent(in) :: i_dst
   Integer*2 :: Chunk(1:Time_dim)
   Real(dp) :: time, p, pw
   integer :: Dset_offset
   Integer :: j,NZero, N_one, Ntwo, ChMax, ChMin, i_chunk
   !
   If(MakePlots) Open(UNIT=10,STATUS='unknown',ACTION='WRITE', &
   !   FILE=trim(DataFolder)//'Stat'//trim(DSet_Names(i_dst))//'.dat')
        FILE=trim(DataFolder)//'RFI_AntR'//trim(DSet_Names(i_dst))//'.dat')
   !flush(unit=2)
   N_zeroChunk=0
   Do i_chunk=0,i_chunkMax
      Dset_offset= i_chunk*Time_dim
      Call GetData(Chunk, Dset_offset, Time_dim) ! ListDataAtt should have been called first
      NZero=COUNT(Chunk.eq.0)
      !Ntwo=COUNT(Chunk.eq.2)
      N_one=COUNT(Chunk.eq.1)
      ChMax=MaxVal(Chunk)
      ChMin=MinVal(Chunk)
      If(N_one*1.0/NZero .lt. .5) then
         !Write(2,"(A,i11,F7.4,A,2i8)") 'Ratio of twos/zeros in chunk=',i_chunk, &
         !  Ntwo*1.0/NZero,', min & max=', ChMin, ChMax
         N_zeroChunk = N_zeroChunk +1
         MaxAmpChunk(i_chunk)=5500  ! greater than saturation
         cycle
      EndIf
      !p=SUM(Chunk(:)*Chunk(:)*1.)/Time_dim
      pw=0.
      Do j=1,Time_dim
         pw=pw+Chunk(j)*Chunk(j)
      Enddo
      pw=sqrt(pw/Time_dim)
      time=1000.*((i_chunk*Time_dim+SAMPLE_NUMBER_first)* Sample -DIPOLE_CALIBRATION_DELAY) ! in [ms]
      If(MakePlots) write(10,*) i_chunk, time, pw, (ChMax-ChMin), N_one, Ntwo
      MaxAmpChunk(i_chunk)=ChMax-ChMin
   EndDo
   !write(2,*) 'number of zero-chunks:',N_zeroChunk,' out of',i_chunkMax
   If(MakePlots) Close(Unit=10)
   !stop
   Return
End Subroutine RoughStatistics
! ----------------------------------------
Subroutine AccumulateBackgrFreq()
! Determine average frequency spectrum using only chuncks for which the raw amplitude is small.
! Since only chuncks for which the (max amplitude) was determined (below above saturation limit),
!      no checks for missing data needs to be done.
   use constants, only : dp,pi,ci, sample
   use DataConstants, only : Time_dim, Cnu_dim
   use HDF5_LOFAR_Read, only : GetData
   use FFT, only : RFTransform_CF, Hann
   Use RFI_MitPars, only : i_chunkMax, MaxAmpChunk, NBackgr, Freq_s, MinAmp
   Implicit none
   Integer, parameter :: NBackgrMax=50 !200  ! Important for background
   Integer*2 :: Chunk(1:Time_dim)
   Real(dp) :: RTime_s(1:Time_dim), q
   Complex(dp) :: CNu_s(0:Cnu_dim)
   Integer :: j, i_chunk, Dset_offset
   Integer :: nu_i, nu_f
   !
   Freq_s(:)=0.
   NBackgr=0
   nu_i=25*Cnu_dim/100
   nu_f=79*Cnu_dim/100
   MinAmp=MinVal(MaxAmpChunk)
   If(MinAmp.gt.3000.) then
      write(2,*) '!!!! Noisy natennas, MinAmp=',MinAmp
      Return  ! Antenna is overloaded or bad
   EndIf
   q=1.2
   j=0
   Do while (j.lt.NBackgrMax*2) ! build in safety margin
      q=q+.05
      j=COUNT(MaxAmpChunk.lt.MinAmp*q)
      !write(2,*) 'count:',j,q
      If(q.gt.2.) then
         Write(2,*) '!!!! need unreasonably large range of amplitudes for background; MinAmp, range:', j, q,MinAmp
         return
      EndIf
   Enddo
   !write(2,*) 'i_chunkMax, j, q, MinAmp:',i_chunkMax, j, q,MinAmp
   MinAmp=MinAmp*q
   !
   Do i_chunk=1,i_chunkMax-1
      If(MaxAmpChunk(i_chunk) .gt. MinAmp) cycle
      If(MaxAmpChunk(i_chunk-1) .gt. MinAmp*q) cycle ! select non-lightning area
      If(MaxAmpChunk(i_chunk+1) .gt. MinAmp*q) cycle
      !write(2,*) 'MaxAmpChunk(i_chunk)', i_chunk, MaxAmpChunk(i_chunk), MinAmp
      Dset_offset=i_chunk*Time_dim
      Call GetData(Chunk, Dset_offset, Time_dim) ! ListDataAtt should be called first
      RTime_s(:)=Chunk(:)*Hann(:)
      Call RFTransform_CF(RTime_s,Cnu_s)
      NBackgr=NBackgr+1
      Freq_s(:)=Freq_s(:) + abs(Cnu_s(:))
      !If(NBackgr.ge.NBackgrMax) exit
   EndDo
   If(NBackgr.gt.1) Then
      Freq_s(:)=Freq_s(:)/NBackgr
   Endif
   write(2,*) 'NBackgr, MinAmp, q:',NBackgr, MinAmp, q !  ,'average of squares, which is larger than square of average'
   !
   Return
End Subroutine AccumulateBackgrFreq
!--------------------------------------------
Subroutine BuildRFIFilter()
! Determine RFI-frequency filter
! Determine power normalization; power in filtered spectrum should be less.
!     This has been checked for some well behaving cases (100 filtered frequencies) and the
!     difference is less than 1%, hardly noticable.
! A more important effect is that the power of a single chunk is given by the sum(over frequency) of squares.
!     The normalization 'powr' is calculated as a sum of squares of average amplitudes. Since an average (or sum) of squares
!     is larger than the square of average amplitudes, the norm, 'powr', is thus too small. Phenomenologically it should be
!     increased by a factor of about 1.12
!
   use constants, only : dp,pi ! ,ci, sample
   use DataConstants, only : Time_dim, Cnu_dim
   Use RFI_MitPars, only : NBackgr, Freq_s, nu_Fltr, powr, Filtring, FiltringRef
   Implicit none
   Real(dp) :: Av, Bv, FiltFact !, FiltPwr
   Integer :: i, nu, nu_i, nu_f, dnu
   !
   If(NBackgr.lt.9) then ! Too low statistics
      write(2,*) '****Insufficient nr of non-zero background blocks****'
      Return
   endif
   nu_i=25*Cnu_dim/100
   nu_f=79*Cnu_dim/100
   dnu=2*Cnu_dim/100
   nu_Fltr(0:nu_i)=0.
   nu_Fltr(nu_i:nu_f)=1.
   nu_Fltr(nu_f:Cnu_dim)=0.
   FiltFact=1.6   ! important for workings of filtering
   Filtring=0
   Powr=0.
   !FiltPwr=0
   !write(2,*) 'power before filter',SUM(Freq_s(nu_i:nu_f)*Freq_s(nu_i:nu_f))
   Do nu=nu_i-dnu,nu_f!+dnu Av+Bv*(
      Av=sum(Freq_s(nu-dnu:nu))/(dnu)
      Bv=0.
      Do i=nu-dnu,nu
         Bv=Bv+Freq_s(i)*(i-nu+dnu/2.)
      Enddo
      !write(2,*) nu,Freq_s(nu), Av, Bv*6/(dnu*dnu), Bv*12/(dnu*dnu*dnu)
      Av=Av+Bv*6/(dnu*dnu)
      !Bv=sum(Freq_s(nu-dnu:nu))/(dnu)
      if(Freq_s(nu).gt.FiltFact*Av) then
         !FiltPwr=FiltPwr + Freq_s(nu)*Freq_s(nu)
         nu_fltr(nu)=0.
         !Freq_s(nu)=Av
         Filtring=Filtring + 1
      Endif
         Powr=Powr + Freq_s(nu)*Freq_s(nu)*nu_fltr(nu)
   Enddo
   ! apply (sum of square) in stead of (square of sum) factor (4/pi)
   Powr=Powr*4./pi
   Freq_s(:)=Freq_s(:)/sqrt(Powr)  ! Normalize the average frequency spectrum
   !Write(2,*) '# filtered frequencies=',Filtring,', RMS=',sqrt(Powr) !,FiltPwr
   !Write(2,*) '# filtered frequencies=',Filtring,', powerN=',sqrt(Powr),dnu
   !If(FiltringRef.lt.0) FiltringRef=Filtring
   !If(Filtring.lt.FiltringRef/5) then
   !   powr=-2.
   !   Write(2,*) '!!!! Unreasonable small number of RFI-lines',Filtring,' below limit of', FiltringRef/5
   !EndIf
   Return
End Subroutine BuildRFIFilter
! -----------------------------------------------
Subroutine PlotFreqFilt(i_dst)
   use constants, only : dp !,pi,ci, sample
   use DataConstants, only : DataFolder, Time_dim, Cnu_dim
   use HDF5_LOFAR_Read, only : DSet_names
   Use RFI_MitPars, only : Freq_s, nu_Fltr
   Implicit none
   Integer, intent(in) :: i_dst
   Integer :: nu, nu_i, nu_f, dnu, nxx
   !Integer :: Filtring
   !Integer, allocatable :: MaxAmpChunk(:)
   !Real*8 :: Powr, LFRAnt_crdnts(3)
   !Logical :: MakePlots=.true.
   Integer :: i,j, Bin_Size
   Real(dp) :: Freq,  Bin
   !CHARACTER(LEN=140) :: SystemCommand
   !
   Open(UNIT=9,STATUS='unknown',ACTION='WRITE', &
      FILE=trim(DataFolder)//'RFI_AntF'//trim(DSet_Names(i_dst))//'.dat')
   !write(2,*) 'output file:',trim(DataFolder)//'RFI_AntF'//trim(DSet_Names(i_dst))//'.dat'
   Bin_size=10
   Do i=0,Cnu_dim/Bin_size-1
      freq=Bin_size*(i+0.5)*100./Cnu_dim
      Bin=0.
      Do j=1,Bin_size
          Bin=Bin+Freq_s(i*Bin_size+j)
      Enddo
      write(9,*) Freq,100.*Bin/Bin_size, minval(nu_fltr(i*Bin_size+1:(i+1)*Bin_size))+.000001, Bin*(nu_fltr(i)+.000001)
   EndDo
   Close(unit=9)
   Return
End Subroutine PlotFreqFilt
!  ---------------------------------
Subroutine TestFilter(i_dst)
   use constants, only : dp,pi, sample
   use DataConstants, only : DataFolder, Time_dim, Cnu_dim
   use HDF5_LOFAR_Read, only : DSet_names
   use HDF5_LOFAR_Read, only : SAMPLE_NUMBER_first, DIPOLE_CALIBRATION_DELAY
   use HDF5_LOFAR_Read, only : GetData
   Use RFI_MitPars, only : powr, nu_Fltr, i_chunkMax    !,    MaxAmpChunk, MinAmp, NBackgr
   use FFT, only : RFTransform_CF2CT, Hann, RFTransform_CF_Filt
   Implicit none
   Integer, intent(in) :: i_dst
   !
   Integer*2 :: Chunk(1:Time_dim)
   Real(dp) :: time, p, pw, RTime_s(1:Time_dim), SubSample_Offset, AmpNrm, ChMax, ChMin
   Complex(dp) :: Cnu_s(0:Cnu_dim), CTime_s(1:Time_dim)
   integer :: Dset_offset, i_chunk
   !Integer :: j,NZero, Ntwo
   !Integer :: nc
   !Real(dp) :: Freq_sum(0:Cnu_dim), p1,p2, AvAmpl
   !
   SubSample_Offset=0.
   AmpNrm=sqrt(1./powr) ! To normalize the complex power per sample to unity
   Open(UNIT=10,STATUS='unknown',ACTION='WRITE', &
      FILE=trim(DataFolder)//'RFI_AntT'//trim(DSet_Names(i_dst))//'.dat')
   !Freq_sum(:)=0 ;   nc=0 ;   p1=0. ;   p2=0.  ;  AvAmpl=0.
   Do i_chunk=0,i_chunkMax
      Dset_offset= i_chunk*Time_dim
      Call GetData(Chunk, Dset_offset, Time_dim) ! ListDataAtt should have been called first
      !
      !NZero=COUNT(Chunk.eq.0)
      !Ntwo=COUNT(Chunk.eq.2)
      !If(Ntwo*1.0/NZero .lt. .66) then
      !   cycle
      !EndIf
      !
      RTime_s(:)=Chunk(:)*Hann(:)*AmpNrm
      Call RFTransform_CF_Filt(RTime_s,nu_fltr,SubSample_Offset,Cnu_s)
      Call RFTransform_CF2CT(Cnu_s, CTime_s )
      RTime_s(:)=REAL(CTime_s(:))
      ChMax=MaxVal(RTime_s)
      ChMin=MinVal(RTime_s)
      p=SUM(RTime_s(:)*RTime_s(:))*2./Time_dim ! factor 2 to account for the imaginary part
      !If(MaxAmpChunk(i_chunk) .lt. MinAmp) then
      !   nc=nc+1
      !   if(nc.le.NBackgr) then
      !      Freq_sum(:)=Freq_sum(:) + abs(Cnu_s(:))
      !      p1=p1+sum(abs(Cnu_s(:))*abs(Cnu_s(:)))
      !      AvAmpl=AvAmpl + (ChMax-ChMin)
      !   EndIf
      !Endif
      time=1000.*((i_chunk*Time_dim+SAMPLE_NUMBER_first)* Sample -DIPOLE_CALIBRATION_DELAY) ! in [ms]
      write(10,*) i_chunk, time, sqrt(p), (ChMax-ChMin) !,pw/p
   EndDo
   Close(Unit=10)
   !Freq_sum(:)=Freq_sum(:)/NBackgr
   !p2=Sum(Freq_sum(:)*Freq_sum(:))*4./pi ! to account for different (sum)^2 and sum(^2)
   !p1=p1/NBackgr
   !AvAmpl=AvAmpl/NBackgr
   !write(2,*) 'AvMaxAmpl,p1,p2:', AvAmpl, p1,p2,', ratio of powers:',p1/p2
   flush(unit=2)
   Return
End Subroutine TestFilter
