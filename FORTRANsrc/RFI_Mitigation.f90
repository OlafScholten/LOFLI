!-------------------------------------
!    Include 'H_ConstantsModules.f90'
!    Include 'FFT_routines.f90'
!    Include 'H_ParamModules.f90'
!    Include 'AntFuncCnst.f90'
!    Include 'AntFunct.f90'
!!    Include 'H_InterferomPars.f90'
!    Include 'HDF5_LOFAR_Read.f90'
!    Include 'H_MappingUtilitiesModules.f90'
!    Include 'H_MappingUtilities.f90'
!    Include 'System_Utilities.f90'
!    Include 'H_GLEplotUtil.f90'
!-----------------------------------------
Module RFI_MitPars
    use constants, only : dp
    use DataConstants, only : DataFolder, Time_dim  !!
    Implicit none
    integer :: NBackgr ! Number of chunck used for frequency filter determination
    Integer :: Filtring ! number of filtered frequencies
    Integer :: FiltringRef=-1 ! set to the number for the first antenna, to put semi-edicated limits on the range
    Integer ::  RefOdd=0, RefEve=0 ! Antenna-ID of reference antennas
    Logical, allocatable :: Background(:)
    Integer, allocatable :: MaxAmpChunk(:)
    Real(dp), allocatable :: powerChunk(:)
    Real(dp) :: MinPow, qp
    Integer :: i_chunkMax
    Integer :: MinAmp ! Amplitude for determining background
    Integer :: N_zeroChunk  ! number of chunks with too many zeros
    Real(dp) :: Powr, q    ! q is a renorm factor to determine when a chunck is background
    Real(dp) :: Freq_lo=10, Freq_hi=90  ! bandwidth [MHz]
    Real(dp) :: FreqNtch_lo=90, FreqNtch_hi=10  ! bandwidth [MHz] for notch filter
    Real(dp) :: nu_Fltr(0:Time_dim/2) ! RFI filter
    Real(dp) :: Freq_s(0:Time_dim/2) ! Average frequency spectrum for background
    Real(dp) :: HBAtile(1:7,0:47,1:3)
    Logical :: MakePlots
    Integer, parameter :: MaximumAntennas=1500
End Module RFI_MitPars
! ---------------------------------------------
program RFI_mitigation
! v18: include search for part of spectrum with no lightning activity
   use constants, only : dp,pi,ci, sample
   use HDF5_LOFAR_Read, only : filename,Group_names, Group_nr, Group_max
   use HDF5_LOFAR_Read, only : DSet_names, DSet_nr, DSet_max, Ant_ID !, STATION_ID
   use HDF5_LOFAR_Read, only : ANTENNA_POSITION, DATA_LENGTH, SAMPLE_NUMBER_first, Absolute_TIME, DIPOLE_CALIBRATION_DELAY
   use HDF5_LOFAR_Read, only : GetFileName, ListGroups, ListDataAtt, ListGroupStructure, CloseFile, CloseGroup, CloseDSet!, GetData
   use FFT, only : RFTransform_su,DAssignFFT !, RFTransform_CF, RFTransform_CF2CT, Hann, RFTransform_CF_Filt
   use Chunk_AntInfo, only : BadAnt_nr_max  ! , Unique_SAI  ! BadAnt_SAI,
   use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, DataFolder, FlashName, Windows
   use DataConstants, only : RunMode ! 1=Explore; 2=Calibrate; 3=Impulsive imager; 4=Interferometry
   Use RFI_MitPars
   use GLEplots, only : GLEplotControl
   Implicit none
   Character(len=20) :: Utility, release !, Char_AbsTimeOffset
   Character(len=5) :: Version='v18  '
   INTEGER :: i, VeryFirstSampl, MostLatestSampl=0
   Character(len=5) :: txt
   Logical :: AllBckgrSpec=.false.  ! produce data files for plotting the backgroung frequency spectra for the reference antennas only
   !Logical :: AllBckgrSpec=.true. ! produce data files for plotting the backgroung frequency spectra for all antennas.
   Logical :: exists, PlentyData, AllCoreAnt=.false.
   !
   Integer :: AntType, i_filR
   Integer :: i_even=0, i_odd=0, MinAmp_Ant(MaximumAntennas), Powr_Ant(MaximumAntennas)
   Real(dp) :: LFRAnt_crdnts(3), MinAmpAv, MinAmpSq, sigma, powrAv, powrSq, powrSigm, power_Even, power_Odd
   !
   Integer :: j,i_fil,i_file,i_grp,i_dst, i_chunk, i_ant !
   Character(LEN=100) :: AntennaFieldsDir='../AntennaFields'   ! Directory where info on antenna positions
   !           can be found if not present in the HDF5 files of the dat
   !
   !AllCoreAnt=.true.
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   i_filR= 2  ! The file containing the reference antenna
   j=IARGC()
   CALL get_environment_variable("AntennaFieldsDir", AntennaFieldsDir)
   write(*,*) 'AntennaFieldsDir:"',TRIM(AntennaFieldsDir),'"'
   !If(j.gt. 0) then
   !   CALL getarg(1, AntennaFieldsDir)
   !   write(*,*) 'argument 1:"',TRIM(AntennaFieldsDir),'"'
   !   !Read(SystemCommand,*) AntennaFieldsDir
   !Endif
   If(j.gt. 1) then
      CALL getarg(2, txt)
      write(*,*) 'argument 2:"',TRIM(txt),'"'
      Read(txt,*) i_filR
   Endif
   !
   OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='RFI_mit.out')
   Utility='RFI-Mitigation'
   release='v25.0, Sept 2025'   ! TRIM(Version)//'; Nov, 2020'
   RunMode=0
   !write(2,"(3x,5(1H-),1x,'RFI_mitigation release of ',A22,25(1H-))") release
   !CALL DATE_AND_TIME (REAL_C(1),REAL_C(2),REAL_C(3),DATE_T)
   !WRITE(2,"(3X,5(1H-),1x,'run on ',I2,'/',I2,'/',I4,' , started at ',&
    !   I2,':',I2,':',I2,'.',I3,1X,25(1H-))") &
   !    DATE_T(3),DATE_T(2),DATE_T(1),(DATE_T(i),i=5,8)
   Call System_Initiation(Utility, release, ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows)  ! set Flash name & folder names
   Write(2,*) 'Antenna information from files in folder: ',TRIM(AntennaFieldsDir)
   !
   If(AllCoreAnt) Then
      Call GetCoordCoreAnt(AntennaFieldsDir)
   EndIf
   AntType=0      ! LBA-outer
   Freq_lo=30
   Freq_hi=80
   Call OneKindRFI_mitigation(AntennaFieldsDir, AntType, i_filR)
   !
   AntType=1      ! HBA
   Freq_lo=200.-170. ! 27
   Freq_hi=200.-109.  ! real freq=200-Freq
   FreqNtch_lo=200.-134.7
   FreqNtch_hi=200.-134.  ! center at nu=134.375
   i_filR=1
   Call OneKindRFI_mitigation(AntennaFieldsDir, AntType, i_filR)
   !
   Call GLEplotControl( CleanPlotFile='*RFI_Ant*.plt')
   Call GLEplotControl( CleanPlotFile='*RFI_Even*.plt')
   Call GLEplotControl( CleanPlotFile='*RFI_Odd*.plt')
   Call GLEplotControl( Submit=.true.)
   !
end program RFI_mitigation
!=================================
Subroutine GetCoordCoreAnt(AntennaFieldsDir)
   use constants, only : dp,pi,ci, sample
   Use RFI_MitPars, only : HBAtile
   Implicit none
   INTEGER :: i, i_stat, i_Type, i_Ant, Ant_ID, STATION_ID  ! , STATIONtype_ID
   Real(dp) :: LFRAnt_crdnts(3)
   Character(LEN=100), intent(in) :: AntennaFieldsDir ! ='../AntennaFields'   ! Directory where info on antenna positions
   i_type=1
   Do i_stat=1,7
      STATION_ID=i_stat*10+i_type
      Do i_Ant=0,47
         Ant_ID=2*i_Ant
         Call GetAntLofarAntPos(LFRAnt_crdnts, AntennaFieldsDir, i_Type, Ant_ID, STATION_ID)
         write(2,"(I4,3F8.2,3x,3F8.2)") i_stat*100+Ant_ID, HBAtile(i_stat,i_Ant,1:3), LFRAnt_crdnts(1:3)
      EndDo
   EndDo
   Stop 'GetCoordCoreAnt'
End Subroutine GetCoordCoreAnt
! ---------------------------------------------
Subroutine OneKindRFI_mitigation(AntennaFieldsDir, AntType, i_filR)
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
   use CPU_timeUsage, only : CPU_usage
   Implicit none
   Character(LEN=100), intent(in) :: AntennaFieldsDir  ! ='../AntennaFields'   ! Directory where info on antenna positions
   integer, intent(in) :: AntType, i_filR
   !Character(len=20) :: Utility, release !, Char_AbsTimeOffset
   Character(len=5) :: Version='v18  '
   INTEGER :: i, VeryFirstSampl, MostLatestSampl  ! , STATIONtype_ID
   Character(len=5) :: txt
   Character(len=1) :: LH
   Character(len=20) :: FileListing
   Character(len=8) :: Label
   Character(len=9) :: LabelE
   Character(len=9) :: LabelO
   Logical :: AllBckgrSpec=.false.  ! produce data files for plotting the backgroung frequency spectra for the reference antennas only
   !Logical :: AllBckgrSpec=.true. ! produce data files for plotting the backgroung frequency spectra for all antennas.
   Logical :: exists, PlentyData
   !
   Integer :: nxx, NFail ! number of antennas for which RFI-determination failed
   Integer :: i_even, i_odd, MinAmp_Ant(MaximumAntennas), Powr_Ant(MaximumAntennas), i_OddEven
   Real(dp) :: LFRAnt_crdnts(3), MinAmpAv, MinAmpSq, sigma, powrAv, powrSq, powrSigm, power_Even, power_Odd
   Real(dp) :: renorm
   !
   Integer :: j,i_fil, i_file,i_grp,i_dst, i_chunk, i_ant, Nr_ants !
   !           can be found if not present in the HDF5 files of the dat
   !
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   VeryFirstSampl=5./sample ! well beyond possible upper limit of 5 seconds
   NFail=0
   i_even=0
   i_odd=0
   i_ant=0
   MostLatestSampl=0
   MinAmpAv=0.
   MinAmpSq=0.
   powrAv=0.
   powrSq=0.
   power_Even=0.
   power_Odd=0.
   !
   If(AntType.eq.0) Then
      LH=' '
      FileListing = 'Book/directory.dat'  ! in the main Fsash directory
      INQUIRE(FILE = TRIM(FileListing), exist=exists)  ! in the main Fsash directory
      If(.not. exists) Then
         Write(2,*) 'No LBA antennas, non existing file "', TRIM(FileListing),'"'
         Return
      EndIf
      write(2,*) '============================= HBA ====================================='
      AllBckgrSpec=.true. ! produce data files for plotting the backgroung frequency spectra for all antennas.
   ElseIf(AntType.eq.1) Then
      LH='H'
      FileListing = 'Book/H_directory.dat'  ! in the main Fsash directory
      INQUIRE(FILE = TRIM(FileListing), exist=exists)  ! in the main Fsash directory
      If(.not. exists) Then
         Write(2,*) 'No HBA antennas, non existing file "', TRIM(FileListing),'"'
         Return
      EndIf
      write(2,*) '============================= HBA ====================================='
      AllBckgrSpec=.true. ! produce data files for plotting the backgroung frequency spectra for all antennas.
   Else
      Write(2,*) '**** Invalid Antenna type:',AntType
      return
   EndIf
   !
   ! Allocate( Unique_StatID(1:Station_nrMax), Unique_SAI(1:Ant_nrMax), Tot_UniqueAnt(0:Station_nrMax) )
   Allocate(  Unique_SAI(1:MaximumAntennas) )
   Unique_SAI(:)=0
   !
   !Call GETCWD(FlashFolder)  ! Get current working directory
   write(2,*) 'Reference antenna from file#',i_filR,' in the list at: /'//FileListing ! trim(DataFolder)
   flush(unit=2)
   !
   If(AntType.eq.0) Then
      Open(unit=12,STATUS='unknown',ACTION='WRITE',FORM ="unformatted", &
         FILE = 'Book/RFI_Filters-'//TRIM(Version)//'.uft')
      Open(unit=14,STATUS='unknown',ACTION='WRITE', &
         FILE = 'Book/LOFAR_H5files_Structure-'//TRIM(Version)//'.dat')
   ElseIf(AntType.eq.1) Then
      Open(unit=12,STATUS='unknown',ACTION='WRITE',FORM ="unformatted", &
         FILE = 'Book/HRFI_Filters-'//TRIM(Version)//'.uft')
      Open(unit=14,STATUS='unknown',ACTION='WRITE', &
         FILE = 'Book/HLOFAR_H5files_Structure-'//TRIM(Version)//'.dat')
   EndIf
   Open(UNIT=16,STATUS='unknown',ACTION='WRITE',  FILE=trim(DataFolder)//TRIM(LH)//'RFI_EvenStat.plt')
   Open(UNIT=17,STATUS='unknown',ACTION='WRITE',  FILE=trim(DataFolder)//TRIM(LH)//'RFI_OddStat.plt')
   Open(UNIT=18,STATUS='unknown',ACTION='WRITE',  FILE=trim(DataFolder)//TRIM(LH)//'RFI_EvenZero.plt')
   Open(UNIT=19,STATUS='unknown',ACTION='WRITE',  FILE=trim(DataFolder)//TRIM(LH)//'RFI_OddZero.plt')
   Open(Unit=20,STATUS='old',ACTION='READ',FILE=TRIM(FileListing),FORM ="formatted")  ! main change is name ************************88
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
   !   Call GetFileName(i_file, nxx)
      Rewind(unit=20)
      Do i=1,i_file
         read(20,"(A)",IOSTAT=nxx) filename
         !Write(2,*) '!file=',filename, nxx,i,i_fil
         if(nxx.ne.0) exit  ! stop('EOF reached in directory')
      Enddo
      if(nxx.ne.0) exit  ! stop('EOF reached in directory')
      WRITE(2,'(A,i3,A,A)') 'Filenr: ',i_file,': ',trim(filename)
      INQUIRE(FILE = trim(filename), exist=exists)  ! in the main Fsash directory
      !If(nxx .ne. 0) exit
      If(exists) Then
         write(2,*) 'Start working on i_file=',i_file,' ,"',trim(filename), '"'
      Else
         write(2,*) 'Not existing: i_file=',i_file,' ,"',trim(filename), '"'
         exit
      EndIf
      write(*,"(A,I0,': ',A)") 'i_file=',i_file,trim(filename)
      !write(14,"(I3,A)") i_file, trim(filename)
      Call ListGroups ! >>>>>>>>>>>>>>>>>>>> opens HDF5 file & extracts group-structure
      write(14,"(i3,2x,3A)")  Group_nr, '"',trim(filename),'"'
      Do i_grp=1,Group_nr
         Write(2,"(A,i0,A,A)") 'Group#',i_grp,' :',trim(Group_names(i_grp))
         Call ListGroupStructure(Group_Names(i_grp)) ! >>>>>>>>>>>>>>>>>>>> opens HDF5 Group in the file
         write(14,"(i3,1x,A)") DSet_nr, trim(Group_Names(i_grp))
         Do i_dst=1, DSet_nr
            Call ListDataAtt(Group_Names(i_grp),DSet_Names(i_dst)) ! >>>>>>>>>>>>>>>>>>>>
            write(14,"(1x,A,I6,I6)")  trim(DSet_Names(i_dst)), STATION_ID, Ant_ID
            STATION_ID=STATION_ID*10 + AntType
            !
            Call CPU_usage(Message='starting '//trim(DSet_Names(i_dst)) )
            If(VeryFirstSampl.gt.SAMPLE_NUMBER_first) VeryFirstSampl=SAMPLE_NUMBER_first
            !
            !Scan file for non-lightning sections
            ! Note that the noise is NOT gaussian since RFI has not been removed yet !
            !DATA_LENGTH=1./sample  ! to shorten time to 1 second for testing
            Powr=-1.
            nu_fltr=0.
            i_chunkMax=DATA_LENGTH/Time_dim -1
            PlentyData=.true.
            If(i_chunkMax.lt.100) Then
               Write(2,"(A,I4)") 'Only',i_chunkMax,' Chunks, data set skipped, all data counted as zeros.'
               N_zeroChunk=i_chunkMax
               powr=0.
               MinAmp=0
               PlentyData=.false. ! not worthwhile
            EndIf
            !
            If(i_ant.lt. MaximumAntennas) then
               i_ant=i_ant+1
               Unique_SAI(i_ant)=STATION_ID*100+ Ant_ID
            Else
               Write(2,*) '***** number of antennas exceeds',MaximumAntennas
               exit
            EndIf
            !
            write(Label,"(I7.7)") Unique_SAI(i_ant)
            !
            If( PlentyData) Then
               If(((RefOdd.eq.0) .and. (mod(Ant_ID,2).eq.1)) .or. ((RefEve.eq.0) .and. (mod(Ant_ID,2).eq.0))) then
                  MakePlots=.true.  !  write statistics to file for the reference station only
               Else
                  MakePlots=.false.
               EndIf
               If(AllBckgrSpec) MakePlots=.true.
               Allocate( MaxAmpChunk(0:i_chunkMax), powerChunk(0:i_chunkMax), Background(0:i_chunkMax))
               Call RoughStatistics(LH,Label)
               !     N_zeroChunk= (number of chucks that have an excess in zero valued samples) is determined for this antenna
               !     MaxAmpChunk(i_chunk) is determined for this antenna
               !
               Call AccumulateBackgrFreq()  !  MinAmp is determined here
            !   DeAllocate (MaxAmpChunk)
               !
               ! Determine averaged frequency spectrum for low-MaxAmpli chuncks
               Call BuildRFIFilter()
               !  Powr is calculated, used for re-scaling amplitude.
               !  This works well for the re-scaling constant when there are no stron RFI-lines in raw spectrum, i.e.
               !     the selected chunks have indeed low background and not only low intensity
               !     in the RFI-lines.
            EndIf
            !
            !Call CPU_usage
            Call GetNorm(renorm)
            write(2,*) 'Renormalisation power changed from',powr,' to',renorm
            powr=renorm
            Call CPU_usage(Message='power normalized')
            MinAmp_Ant(i_ant)=MinAmp
            Powr_Ant(i_ant)=NINT(powr)
            !
            If((powr.le.0.5) .or. .not.PlentyData) then
               Nfail=Nfail+1
               If(Nfail.lt.BadAnt_nr_max) BadAnt_SAI(Nfail)=Unique_SAI(i_ant)
               Write(2,*) 'Quit further annalysis of this antenna:', BadAnt_SAI(Nfail), powr
            Else
               MinAmpAv=MinAmpAv + MinAmp
               MinAmpSq=MinAmpSq + MinAmp*MinAmp
               powrAv=powrAv + powr
               powrSq=powrSq + powr *powr
               !
               write(2,"(A, I3,A,I7)") ' Unique_SAI(',i_ant,')=', Unique_SAI(i_ant)
               !
               If(VeryFirstSampl.gt.SAMPLE_NUMBER_first) VeryFirstSampl=SAMPLE_NUMBER_first
               If(MostLatestSampl.lt.(SAMPLE_NUMBER_first+DATA_LENGTH)) MostLatestSampl=(SAMPLE_NUMBER_first+DATA_LENGTH)
               !
               ! Test effect filter
               If(i_file.eq.i_filR .and. ((RefOdd.eq.0) .and. (mod(Ant_ID,2).eq.1)) ) then  ! The reference station
                  RefOdd=i_ant
                  write(2,*) 'odd reference antenna',trim(DSet_Names(i_dst))
    !              Call PlotFreqFilt(LH,Label)  ! Frequency spectra
                  i_OddEven=i_odd+1
                  Call TestFilter(LH,Label, i_OddEven)     ! Filtered time spectra
               Else If(i_file.eq.i_filR .and. ((RefEve.eq.0) .and. (mod(Ant_ID,2).eq.0)) ) then
                  RefEve=i_ant
                  write(2,*) 'even reference antenna',trim(DSet_Names(i_dst))
      !            Call PlotFreqFilt(LH,Label)  ! Frequency spectra
                  Call TestFilter(LH,Label, i_OddEven)     ! filtered spectra
                  i_OddEven=i_even+1
               Else If(AllBckgrSpec) Then
                  !write(2,*) 'Working on frequency spectrum for ',trim(DSet_Names(i_dst))
       !           Call PlotFreqFilt(LH,Label)  ! Frequency spectra
                  If(mod(Ant_ID,2).eq.1) Then
                     i_OddEven=i_odd+1
                  Else
                     i_OddEven=i_even+1
                  EndIf
                  Call TestFilter(LH,Label, i_OddEven)    ! Filtered time spectra
               EndIf
            Endif
            !
            If( PlentyData ) DeAllocate (MaxAmpChunk, powerChunk, Background)
            !Call CPU_usage
            !
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
            Call GetAntLofarAntPos(LFRAnt_crdnts, AntennaFieldsDir, AntType, Ant_ID, STATION_ID)
            Write(2,*) 'LOFAR coordinates=',LFRAnt_crdnts
            Write(12) i_fil,i_grp,i_dst,LFRAnt_crdnts, powr, &
              DATA_LENGTH, SAMPLE_NUMBER_first, Absolute_TIME, DIPOLE_CALIBRATION_DELAY,nu_fltr
            Call CloseDSet !<<<<<<<<<<<<<<<<<<<<
            Call CPU_usage(Message='Antenna finalized' )
         Enddo ! i_dst=1,DSet_nr
         Call CloseGroup !<<<<<<<<<<<<<<<<<<<<
         Flush(unit=14)
        EndDo ! ListGroups
        Call CloseFile !<<<<<<<<<<<<<<<<<<<<
   Enddo
! 998 Continue
   Close(unit=20)
   Close(unit=14)
   !Close(unit=15)
   Close(unit=16)
   Close(unit=17)
   Close(unit=18)
   Close(unit=19)
   Close(unit=12)
   call DAssignFFT()
   write(2,*) 'Ave. power in even & odd antennas:',Power_Even/i_even, power_Odd/i_odd, Power_Even/power_Odd
   Flush(unit=2)
   !
   If(AllBckgrSpec) Then
      Do i_ant=1,i_odd+i_even-1
         If(mod(Unique_SAI(i_ant),2).eq.1) cycle
         write(2,*) i_ant, Unique_SAI(i_ant), Unique_SAI(i_ant+1)
         If(Unique_SAI(i_ant)+1 .ne. (Unique_SAI(i_ant+1)) ) cycle
         If(Powr_Ant(i_ant).lt.1) cycle
         If(Powr_Ant(i_ant+1).lt.1) cycle
         write(LabelO,"(1x,I7.7)") Unique_SAI(i_ant+1)
         write(LabelE,"(1x,I7.7)") Unique_SAI(i_ant)
         write(Label,"(I3.3,'_',I2.2)") Unique_SAI(i_ant)/1000, MODULO(Unique_SAI(i_ant),100)
         write(2,*) 'plot antenna',i_ant, LabelE, LabelO,' , PlotName=','RFI_ant-'//TRIM(FlashName)//TRIM(LH)//TRIM(Label)
         Call GLEplotControl(PlotType='H_RFI_Plots', PlotName='RFI_ant-'//TRIM(FlashName)//TRIM(LH)//TRIM(Label), &
                  PlotDataFile=TRIM(DataFolder)//LH//LabelE//LabelO)
      EndDo
   Else
      write(LabelE,"(1x,I7.7)") Unique_SAI(RefEve)
      write(LabelO,"(1x,I7.7)") Unique_SAI(RefOdd) ! RefOdd=0, RefEve=0
      Call GLEplotControl(PlotType='H_RFI_Plots', PlotName='RFI_Plots'//TRIM(FlashName)//TRIM(LH), &
               PlotDataFile=TRIM(DataFolder)//LH//LabelE//LabelO)
      write(2,*) '===== config written for:', 'RFI_Plots'//TRIM(FlashName)//TRIM(LH)
   EndIf
   !
   !Write(2,*) 'total nr antennas and stations'
   !Write(2,*) 'Check for bad antennas occurring twice; indicate for bad antennas why they are bad'
   write(2,*) 'time range for which there are data [ms]',VeryFirstSampl*sample*1000, MostLatestSampl*sample*1000.
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
   If(i_ant.gt. MaximumAntennas) then
      Write(2,*) '!!!!!!!!!! Nr of antennas exceeds max:',I_ant
      i_ant=MaximumAntennas
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
   DeAllocate(  Unique_SAI )
   !CALL SYSTEM(SystemCommand,nxx)
   !write(2,*) nxx,'=zero for success '
   RefOdd=0
   RefEve=0
   Return
   !
End Subroutine OneKindRFI_mitigation
!=================================
Subroutine GetAntLofarAntPos(LFRAnt_crdnts, AntennaFieldsDir, AntType, Ant_ID, STATION_ID)
! Base files from https://github.com/lofar-astron/lofar-antenna-positions
   use constants, only : dp, pi
   Use Chunk_AntInfo, only : LOFAR_number, LOFARNrMax
   use StationMnemonics, only : Station_ID2Mnem
   Implicit none
   Real(dp), intent(out) :: LFRAnt_crdnts(1:3)
   Integer, intent(in) :: AntType, Ant_ID, STATION_ID
   Character(Len=100), intent(in) :: AntennaFieldsDir
   !Character(len=6) :: Station_Mnem
   !Character(len=5) :: txt, LOFAR_Mnem
   !Character(len=3) :: LBA
   Logical, save :: First=.true.
   !Real(dp), save :: ETRS_RotMat(3,3),Center_CS002(3)
   Real(dp), save :: LBA_crdnts(1:38,0:47,1:3),HBA_crdnts(1:38,0:47,1:3)
   !real(dp) :: Rotmatphi(3,3),Rotmatth(3,3)
   !real(dp) :: x,theta, phi
   !Real(dp) :: ETRS(1:3)
   !INTEGER :: i,j,k,nxx,X_ant, Y_ant, i_LOFAR
   INTEGER :: i_LOFAR
   !Real(dp) :: central(3), LFRAnt_crdnts_new(1:3)
   !
   If(First) Then
      Call GetAllLofarAntPos(LBA_crdnts,HBA_crdnts, AntennaFieldsDir)
      first=.false.
   EndIf
   !
   Do i_LOFAR=1,LOFARNrMax
      If(STATION_ID/10 .eq. LOFAR_number(i_LOFAR)) exit
   enddo
   !
   If(AntType.eq.0) Then
      LFRAnt_crdnts(1:3)=LBA_crdnts(i_LOFAR,Ant_ID/2,1:3) ! LBA-outer only
   Else
      LFRAnt_crdnts(1:3)=HBA_crdnts(i_LOFAR,Ant_ID/2,1:3) ! HBA-dipole
   EndIf
   !write(2,*) '!GetAntLofarAntPos, new:', i_LOFAR,Ant_ID,AntType, LFRAnt_crdnts_new(1:3)
   !
   !
   !write(2,*) '!GetAntLofarAntPos:',STATION_ID,Station_Mnem,' ', LOFAR_Mnem,LFRAnt_crdnts-LFRAnt_crdnts_new
   !
   Return
End Subroutine GetAntLofarAntPos
!------------------
!=================================
Subroutine GetAllLofarAntPos(LBA_crdnts,HBA_crdnts, AntennaFieldsDir)
! Base files from https://github.com/lofar-astron/lofar-antenna-positions
   use constants, only : dp, pi
   use HDF5_LOFAR_Read, only : Ant_ID, STATION_ID
   use StationMnemonics, only : Station_ID2Mnem, LOFAR_Mnem2ID
   Use Chunk_AntInfo, only : LOFAR_name, LOFAR_number, LOFARNrMax
   Use RFI_MitPars, only : HBAtile
   Implicit none
   Real(dp), intent(out) :: LBA_crdnts(1:38,0:47,1:3),HBA_crdnts(1:38,0:47,1:3)
 !  Integer, intent(in) :: AntType
   Character(Len=100), intent(in) :: AntennaFieldsDir
 !  Character(len=6) :: Station_Mnem
   Character(len=5) :: txt, LOFAR_Mnem
   Character(len=3) :: LBA
  ! Logical, save :: First=.true.
   Real(dp) :: ETRS_RotMat(3,3),Center(3)
   real(dp) :: Rotmatphi(3,3),Rotmatth(3,3)
   real(dp) :: x,theta, phi
   Real(dp) :: ETRS(1:3),LBA_ETRS(1:38,0:95,1:3),HBA_ETRS(1:38,0:47,1:3)
   INTEGER :: i,j,k,nxx, i_ant, i_LOFAR, Idip(0:3,0:3)
   Real(dp) :: central(3),HDipxy(0:47,0:1), DipoleSpacing(0:1,0:1,1:3)
   !
   Open(unit=11, STATUS='old',ACTION='read', FILE = TRIM(AntennaFieldsDir)//'/etrs-antenna-positions.csv')
   !
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
   LOFAR_Mnem=' '
   Do
      Read(11,*,IOSTAT=nxx) txt,LBA,i_ant,ETRS
      If(nxx.ne.0) exit
      If(txt.ne.LOFAR_Mnem) Then
          LOFAR_Mnem=txt
          Do i_LOFAR=1,LOFARNrMax
              If(LOFAR_Mnem(3:5) .eq. LOFAR_name(i_LOFAR)(3:5)) then
                  exit
              endif
          enddo
      EndIf
      If(LBA.eq.'LBA') Then
         LBA_ETRS(i_LOFAR,i_ant,1:3)=ETRS(1:3)
      ElseIf(LBA.eq.'HBA') Then
         HBA_ETRS(i_LOFAR,i_ant,1:3)=ETRS(1:3)
      Else
         Write(2,*) '!GetAllLofarAntPos, error??', LBA
      EndIf
   enddo
   Close(unit=11)
   !
   Open(unit=11, STATUS='old',ACTION='read', FILE = TRIM(AntennaFieldsDir)//'/HBAs.txt')
   Do i_ant=0,47
      Read(11,*) Idip !
      Do i=0,3
         Do j=0,3
            If(Idip(i,j).eq.1) goto 1
         EndDo
      EndDo
   1  continue
      HDipxy(i_ant,0)=i-1.5
      HDipxy(i_ant,1)=j-1.5
      Read(11,*) Idip ! skip the even records
   EndDo
   Close(unit=11)
   !
   !
   center(1:3)=LBA_ETRS(2,0,1:3)  ! center CS002-LBA dipole 0&1
   Do i_LOFAR=1,38
      Do i_ant=0,47
         Do i=1,3
           x=0.
           Do j=1,3
               x=x + ETRS_RotMat(i,j) * (LBA_ETRS(i_LOFAR,i_ant+48,j)-center(j))
           enddo
           LBA_crdnts(i_LOFAR,i_ant,i)=x  ! 1=North, 2=East, 3=vertical(plumbline)
         enddo
         !If(i_LOFAR.eq.2) write(2,*) '!GetAllLofarAntPos, LBA:', i_LOFAR,i_ant, LBA_crdnts(i_LOFAR,i_ant,1:3)
      EndDo
      !
      ! Get HBA-tile centers
      Do i_ant=0,47
         Do i=1,3
           x=0.
           Do j=1,3
               x=x + ETRS_RotMat(i,j) * (HBA_ETRS(i_LOFAR,i_ant,j)-center(j))
           enddo
           HBA_crdnts(i_LOFAR,i_ant,i)=x  ! 1=North, 2=East, 3=vertical(plumbline)
         enddo
      EndDo
      !
      If(i_LOFAR.le.7) HBAtile(i_LOFAR,:,:)=HBA_crdnts(i_LOFAR,:,:)
      !
      !  Get single HBA-dipole positions
      If(i_LOFAR.lt.22) Then ! core stations have two ears
         DipoleSpacing(0,0,1:3)=(HBA_crdnts(i_LOFAR,11,1:3)-HBA_crdnts(i_LOFAR,6,1:3))/20.      ! Ear number 0, x
         DipoleSpacing(0,1,1:3)=(HBA_crdnts(i_LOFAR,22,1:3)-HBA_crdnts(i_LOFAR,0,1:3))/20.      ! Ear number 0, y
         DipoleSpacing(1,0,1:3)=(HBA_crdnts(i_LOFAR,35,1:3)-HBA_crdnts(i_LOFAR,30,1:3))/20.      ! Ear number 1, x
         DipoleSpacing(1,1,1:3)=(HBA_crdnts(i_LOFAR,46,1:3)-HBA_crdnts(i_LOFAR,24,1:3))/20.      ! Ear number 1, y
      Else ! Remote stations have a single ear
         DipoleSpacing(0,0,1:3)=(HBA_crdnts(i_LOFAR,23,1:3)-HBA_crdnts(i_LOFAR,16,1:3))/28.      ! Ear number 0, x
         DipoleSpacing(0,1,1:3)=(HBA_crdnts(i_LOFAR,45,1:3)-HBA_crdnts(i_LOFAR,1,1:3))/28.      ! Ear number 0, y
         DipoleSpacing(1,0,1:3)=DipoleSpacing(0,0,1:3)      ! Ear number 1, x
         DipoleSpacing(1,1,1:3)=DipoleSpacing(0,1,1:3)      ! Ear number 1, y
      Endif
      !
      ! Place the HBA-dipoles at the right position
      Do i_ant=0,47
         k=i_ant/24
         !If(i_LOFAR.eq.2) write(2,*) '!GetAllLofarAntPos, HBA-tile:', i_ant, HBA_crdnts(i_LOFAR,i_ant,1:3)
         HBA_crdnts(i_LOFAR,i_ant,1:3)=HBA_crdnts(i_LOFAR,i_ant,1:3) + &
            HDipxy(i_ant,0)*DipoleSpacing(k,0,1:3) + HDipxy(i_ant,1)*DipoleSpacing(k,1,1:3)
         !If(i_LOFAR.eq.2) write(2,*) '!GetAllLofarAntPos, HBA-dipl:', HDipxy(i_ant,0:1), &
         !   HDipxy(i_ant,0)*DipoleSpacing(k,0,1:3) + HDipxy(i_ant,1)*DipoleSpacing(k,1,1:3)
      EndDo
   EndDo
   !
   Return
End Subroutine GetAllLofarAntPos
!------------------
Subroutine RoughStatistics(LH, Label)
!Scan file for non-lightning sections
! Note that the noise is NOT gaussian since RFI has not been removed yet !
!     N_zeroChunk= (number of chucks that have an excess in zero valued samples) is determined for this antenna
!     MaxAmpChunk(i_chunk) is determined for this antenna
   use constants, only : dp,pi,ci, sample
   use DataConstants, only : DataFolder, Time_dim
   use HDF5_LOFAR_Read, only : DSet_names!, DATA_LENGTH
   use HDF5_LOFAR_Read, only : SAMPLE_NUMBER_first, DIPOLE_CALIBRATION_DELAY
   use HDF5_LOFAR_Read, only : GetData
   use FFT, only : RFTransform_CF, Hann
   Use RFI_MitPars, only : i_chunkMax, MaxAmpChunk, MakePlots, N_zeroChunk
   Use RFI_MitPars, only : powerChunk, MinPow, Background
   Use RFI_MitPars, only : Freq_lo, Freq_hi! , FreqNtch_lo, FreqNtch_hi !, nu_Fltr ! , powr, Filtring, FiltringRef
   Implicit none
   Character(len=*), intent(in) :: LH, Label
   Integer*2 :: Chunk(1:Time_dim)
   Real(dp) :: RTime_s(1:Time_dim)
   Complex(dp) :: CNu_s(0:Time_dim/2)
   Real(dp) :: time, p, pw, SubSample_Offset
   integer :: Dset_offset
   Integer :: j,NZero, N_one, Ntwo, ChMax, ChMin, i_chunk
   Integer ::  nu_i, nu_f !, nuc_i, nuc_f
   !
   If(MakePlots) Open(UNIT=11,STATUS='unknown',ACTION='WRITE', &
   !   FILE=trim(DataFolder)//'Stat'//trim(DSet_Names(i_dst))//'.dat')
        FILE=trim(DataFolder)//TRIM(LH)//'RFI_AntR'//trim(Label)//'.plt')
   !If(MakePlots) write(2,*) '! writing plotdata to file:',trim(DataFolder)//TRIM(LH)//'RFI_AntR'//trim(Label)//'.plt', MakePlots
   !flush(unit=2)
   nu_i=Freq_lo*(Time_dim/2)/100
   nu_f=Freq_hi*(Time_dim/2)/100
!   nuc_i=FreqNtch_lo*(Time_dim/2)/100
!   nuc_f=FreqNtch_hi*(Time_dim/2)/100
   ! include notch filter!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   nu_Fltr(0:nu_i)=0.
!   nu_Fltr(nu_i:nu_f)=1.
!   nu_Fltr(nu_f:Time_dim/2)=0.
   !
   N_zeroChunk=0
   SubSample_Offset=0.
   Do i_chunk=0,i_chunkMax
      Dset_offset= i_chunk*Time_dim
      Call GetData(Chunk, Dset_offset, Time_dim) ! ListDataAtt should have been called first
      NZero=COUNT(Chunk.eq.0)
      !Ntwo=COUNT(Chunk.eq.2)
      N_one=COUNT(Chunk.eq.1)
      If(N_one*1.0/NZero .lt. .5) then  !  too many zero samples in this trace
         !Write(2,"(A,i11,F7.4,A,2i8)") 'Ratio of twos/zeros in chunk=',i_chunk, &
         !  Ntwo*1.0/NZero,', min & max=', ChMin, ChMax
         N_zeroChunk = N_zeroChunk +1
         MaxAmpChunk(i_chunk)=5500  ! greater than saturation
         powerChunk(i_chunk)=5500
         cycle
      EndIf
      ChMax=MaxVal(Chunk)
      ChMin=MinVal(Chunk)
      !
      RTime_s(:)=Chunk(:)*Hann(:)
      Call RFTransform_CF(RTime_s,Cnu_s)
      !Call RFTransform_CF2CT(Cnu_s, CTime_s )
      powerChunk(i_chunk)=SUM( Cnu_s(nu_i:nu_f)*conjg(Cnu_s(nu_i:nu_f)) )! *2./Time_dim ! factor 2 to account for the imaginary part
      !
      pw=0.
      Do j=1,Time_dim
         pw=pw+Chunk(j)*Chunk(j)
      Enddo
      pw=sqrt(pw/Time_dim)
      time=1000*((i_chunk*Time_dim+SAMPLE_NUMBER_first)* Sample -DIPOLE_CALIBRATION_DELAY) ! in [ms]
      If(MakePlots) write(11,*) i_chunk, time, pw, (ChMax-ChMin), N_one, Ntwo, powerChunk(i_chunk)
      MaxAmpChunk(i_chunk)=ChMax-ChMin
   EndDo
   !write(2,*) 'number of zero-chunks:',N_zeroChunk,' out of',i_chunkMax
   If(MakePlots) Close(Unit=11)
   MinPow=MinVal(powerChunk(0:i_chunkMax) )
   !stop
   Return
End Subroutine RoughStatistics
! ----------------------------------------
Subroutine AccumulateBackgrFreq()
! Determine average frequency spectrum using only chuncks for which the raw amplitude is small.
! Since only chuncks for which the (max amplitude) was determined (below above saturation limit),
!      no checks for missing data needs to be done.
   use constants, only : dp,pi,ci, sample
   use DataConstants, only : Time_dim
   use HDF5_LOFAR_Read, only : GetData
   use FFT, only : RFTransform_CF, Hann
   Use RFI_MitPars, only : i_chunkMax, MaxAmpChunk, NBackgr, Freq_s, MinAmp, q
   Use RFI_MitPars, only : powerChunk, MinPow, qp,    Background
   Implicit none
   Integer, parameter :: NBackgrMax=50 !200  ! Important for background
   Integer, parameter :: NBackgrMin=1000  ! =300 ms !Important for background
   Integer*2 :: Chunk(1:Time_dim)
   Real(dp) :: RTime_s(1:Time_dim), BkgrThreh
   Complex(dp) :: CNu_s(0:Time_dim/2)
   Integer :: j, i_chunk, Dset_offset
   !
   Freq_s(:)=0.
   NBackgr=0
!   MinAmp=MinVal(MaxAmpChunk)
!   If(MinAmp.gt.3000.) then
!      write(2,*) '!!!! Noisy natennas, MinAmp=',MinAmp
!      Return  ! Antenna is overloaded or bad
!   EndIf
   q=1.9 ! taking this too low has the property of missing some strong RFI lines that are there for only part of the time
   ! even steeing this to 5 does not resolve sufficiently the problem with the strong & imtermittet RFI line
   j=0
   Do while (j.lt.NBackgrMin) ! build in safety margin and set the amplitude range over which the noise is normalized
!   Do while (j.lt.NBackgrMax*2) ! build in safety margin and set the amplitude range over which the noise is normalized
      q=q+.1
!      j=COUNT(MaxAmpChunk.lt.MinAmp*q)
      j=COUNT(powerChunk(:).lt.MinPow*q)
      !write(2,*) 'count:',j,q
      If(q.gt.10.) then
         Write(2,*) '!!!! need unreasonably large range of amplitudes for background; MinAmp, range:', j, q,MinPow,MinAmp
   !      return
         exit
      EndIf
   Enddo
   !write(2,*) '! Determine range of maximum amplitudes per chunk for which background power is determined:',i_chunkMax, q,MinAmp,&
   !   '(Obsolete)'
!   BkgrThreh=MinAmp*q
   BkgrThreh=MinPow*q
   !
   !qp=1.2
   !BkgrThreh=MinPow+5.
   !q=BkgrThreh/MinPow
   !write(2,*) '!Filtering on power'
   Background(:)=.false.
   Do i_chunk=1,i_chunkMax-1
   !   If(MaxAmpChunk(i_chunk) .gt. MinAmp) cycle
   !   If(MaxAmpChunk(i_chunk-1) .gt. MinAmp*q) cycle ! select non-lightning area
   !   If(MaxAmpChunk(i_chunk+1) .gt. MinAmp*q) cycle
      If(powerChunk(i_chunk) .gt. BkgrThreh) cycle
      If(powerChunk(i_chunk-1) .gt. BkgrThreh) cycle ! select non-lightning area
      If(powerChunk(i_chunk+1) .gt. BkgrThreh) cycle
      Background(i_chunk)=.true.
      !write(2,*) 'MaxAmpChunk(i_chunk)', i_chunk, MaxAmpChunk(i_chunk), MinAmp
      Dset_offset=i_chunk*Time_dim
      Call GetData(Chunk, Dset_offset, Time_dim) ! ListDataAtt should be called first
      RTime_s(:)=Chunk(:)*Hann(:)
      Call RFTransform_CF(RTime_s,Cnu_s)
      NBackgr=NBackgr+1
      Freq_s(:)=Freq_s(:) + abs(Cnu_s(:))
      !If(NBackgr.ge.NBackgrMax) exit
      !write(2,*) '!AccumulateBackgrFreq, i_chunk=', i_chunk
   EndDo
   If(NBackgr.gt.1) Then
      Freq_s(:)=Freq_s(:)/NBackgr
   Endif
!   write(2,"(A,I5,A,I6, A,F6.2, A,F5.0,A)") 'Background spectrum calculated using', NBackgr, ' out of', i_chunkMax, &
!      ' data chunks for which the (unfiltered max. amlitude/chunk) is below a factor ',q, &
!      ' greater than minimum of', MinAmp/q,' [ADC counts].'
   write(2,"(A,I5,A,I6, A,F6.2, A,F8.3,A)") 'Background spectrum calculated using', NBackgr, ' out of', i_chunkMax, &
      ' data chunks for which the (unfiltered power/chunk) is below a factor ',q, &
      ' greater than minimum of', MinPow
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
!     increased by a factor of about 1.12 = sqrt(4/pi) since Freq_s(:)=Freq_s(:) + abs(Cnu_s(:))
!
   use constants, only : dp,pi ! ,ci, sample
   use DataConstants, only : Time_dim
   Use RFI_MitPars, only : NBackgr, Freq_lo, Freq_hi, FreqNtch_lo, FreqNtch_hi, Freq_s, nu_Fltr, powr, Filtring, FiltringRef
   Implicit none
   Real(dp) :: Av, Bv, FiltFact !, FiltPwr
   Integer :: i, nu, nu_i, nu_f, dnu, nuc_i, nuc_f
   !
   If(NBackgr.lt.9) then ! Too low statistics
      write(2,*) '****Insufficient nr of non-zero background blocks****'
      Return
   endif
   !
   nuc_i=FreqNtch_lo*(Time_dim/2)/100
   nuc_f=FreqNtch_hi*(Time_dim/2)/100
   nu_i=Freq_lo*(Time_dim/2)/100
   nu_f=Freq_hi*(Time_dim/2)/100
   dnu=2*Time_dim/2/100
   nu_Fltr(0:nu_i)=0.
   nu_Fltr(nu_i:nu_f)=1.
   nu_Fltr(nu_f:Time_dim/2)=0.
   FiltFact=1.6   ! important for workings of filtering
   Filtring=0
   Powr=0.
   !write(2,*) '!BuildRFIFilter, lo,hi:',200.-nu_i*200./Time_dim,200.-nu_f*200./Time_dim
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
      If(nu.gt.nuc_i .and. nu .lt. nuc_f) Av=0 ! 21504 ! takes out RFI-line that causes jump in the power spectrum
      if(Freq_s(nu).gt.FiltFact*Av) then
         !FiltPwr=FiltPwr + Freq_s(nu)*Freq_s(nu)
         nu_fltr(nu)=0.
         !Freq_s(nu)=Av
         Filtring=Filtring + 1
         !write(2,*) '!BuildRFIFilter, nu=',Filtring, nu, 200.-nu*200./Time_dim
         !write(2,"(A,I5,F4.1,I6,F9.3)") '!BuildRFIFilter, nu=',Filtring, nu_fltr(nu), nu, 200.-nu*200./Time_dim
      Endif
      !If(Filtring.lt. 300.) nu_fltr(nu)=1
      Powr=Powr + Freq_s(nu)*Freq_s(nu)*nu_fltr(nu)
   Enddo  ! nu=nu_i-dnu,nu_f!
   ! apply (sum of square) in stead of (square of sum) factor (4/pi) since Freq_s(:)=Freq_s(:) + abs(Cnu_s(:))  and later /NBackgr
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
!-----------------------------------------
Subroutine GetNorm(renorm)
   use constants, only : dp !,pi, sample
   use DataConstants, only : Time_dim
   use HDF5_LOFAR_Read, only : GetData
   Use RFI_MitPars, only : powr, nu_Fltr, i_chunkMax!, MinAmp, q, MaxAmpChunk
   Use RFI_MitPars, only : Background!, NBackgr
   use FFT, only : RFTransform_CF, Hann
   Implicit none
   !
   Integer*2 :: Chunk(1:Time_dim)
   Real(dp) :: pwr, renorm, BckgrLev
   !Real(dp) :: time, p, pw, RTime_s(1:Time_dim), SubSample_Offset, AmpNrm, ChMax, ChMin, BckgrLev, BckgrPow
   Real(dp) :: RTime_s(1:Time_dim), SubSample_Offset
   Complex(dp) :: Cnu_s(0:Time_dim/2)
   integer :: Dset_offset, i_chunk, NB
   SubSample_Offset=0.
   BckgrLev=1.1
   NB=0
   renorm=0.
   Do i_chunk=0,i_chunkMax
      If(.not. Background(i_chunk) ) cycle
      Dset_offset= i_chunk*Time_dim
      Call GetData(Chunk, Dset_offset, Time_dim) ! ListDataAtt should have been called first
      !
      RTime_s(:)=Chunk(:)*Hann(:)
      Call RFTransform_CF(RTime_s,Cnu_s)
      Cnu_s(:)=Cnu_s(:)*nu_Fltr(:)
      pwr=SUM( Cnu_s(:)*conjg(Cnu_s(:)) )! *2./Time_dim ! factor 2 to account for the imaginary part
      !
      If(pwr.gt.powr*BckgrLev) cycle
      NB=NB+1
      renorm=renorm+Pwr
   EndDo
   If(NB.lt.1) Then
      Renorm=-1
      write(2,*) '!renormalization not possible'
   Else
      renorm=renorm/NB
      !write(2,*) '!renorm:',renorm, powr, NB
   EndIf
   return
End Subroutine GetNorm
! -----------------------------------------------
Subroutine PlotFreqFilt(LH, Label, FreqBackgr)
   use constants, only : dp !,pi,ci, sample
   use DataConstants, only : DataFolder, Time_dim
   !use HDF5_LOFAR_Read, only : DSet_names
   Use RFI_MitPars, only : Freq_s, nu_Fltr
   Implicit none
   Character(len=*), intent(in) :: LH, Label
   Real(dp), intent(in) :: FreqBackgr(*)
   !Integer :: nu, nu_i, nu_f, dnu, nxx
   !Integer :: Filtring
   !Integer, allocatable :: MaxAmpChunk(:)
   !Real*8 :: Powr, LFRAnt_crdnts(3)
   !Logical :: MakePlots=.true.
   Integer :: i,j, Bin_Size
   Real(dp) :: Freq,  Bin, BinB
   !CHARACTER(LEN=140) :: SystemCommand
   !
   Open(UNIT=9,STATUS='unknown',ACTION='WRITE', &
      FILE=trim(DataFolder)//TRIM(LH)//'RFI_AntF'//trim(Label)//'.plt')
   !write(2,*) 'plot data to file:',trim(DataFolder)//TRIM(LH)//'RFI_AntF'//trim(Label)//'.plt'
   Bin_size=10
   Do i=0,(Time_dim/2)/Bin_size-1
      If(LH.eq.'H') Then
         freq=200-Bin_size*(i+0.5)*100./(Time_dim/2)
      Else
         freq=Bin_size*(i+0.5)*100./(Time_dim/2)
      EndIf
      Bin=0.
      BinB=0.
      Do j=1,Bin_size
          Bin=Bin+Freq_s(i*Bin_size+j)
          BinB=BinB+FreqBackgr(i*Bin_size+j)
      Enddo
      write(9,*) Freq,100.*Bin/Bin_size, minval(nu_fltr(i*Bin_size+1:(i+1)*Bin_size))+.000001, Bin*(nu_fltr(i)+.000001), &
             100.*BinB/Bin_size
   EndDo
   Close(unit=9)
   Return
End Subroutine PlotFreqFilt
!  ---------------------------------
Subroutine TestFilter(LH, Label,i_OddEven)
   use constants, only : dp,pi, sample
   use DataConstants, only : DataFolder, Time_dim
   !use HDF5_LOFAR_Read, only : DSet_names
   use HDF5_LOFAR_Read, only : SAMPLE_NUMBER_first, DIPOLE_CALIBRATION_DELAY
   use HDF5_LOFAR_Read, only : GetData
   use HDF5_LOFAR_Read, only : Ant_ID, STATION_ID
   use StationMnemonics, only : Station_ID2Mnem, Statn_ID2Mnem
   Use RFI_MitPars, only : powr, nu_Fltr, i_chunkMax, MinAmp, q, MaxAmpChunk
   Use RFI_MitPars, only : powerChunk, MinPow, qp,    Background, NBackgr, Filtring
   use FFT, only : RFTransform_CF2CT, Hann, RFTransform_CF_Filt
   Implicit none
   Character(len=*), intent(in) :: LH, Label
   Integer, intent(in) :: i_OddEven
   !
   Integer*2 :: Chunk(1:Time_dim)
   Real(dp) :: time, p, pw, RTime_s(1:Time_dim), SubSample_Offset, AmpNrm, ChMax, ChMin, BckgrLev, BckgrPow
   Complex(dp) :: Cnu_s(0:Time_dim/2), CTime_s(1:Time_dim)
   integer :: Dset_offset, i_chunk, NBackgrInc, NB
   !
   Real(dp) :: FreqBackgr(0:Time_dim/2) ! Average frequency spectrum for background
   Integer :: i,j, Bin_Size
   Real(dp) :: Freq,  Bin
   !
   SubSample_Offset=0.
   AmpNrm=sqrt(1./powr) ! To normalize the complex power per sample to unity
   Open(UNIT=11,STATUS='unknown',ACTION='WRITE', &
      FILE=trim(DataFolder)//TRIM(LH)//'RFI_AntT'//trim(Label)//'.plt')
   !write(2,*) '! writing plotdata to file:',trim(DataFolder)//TRIM(LH)//'RFI_AntT'//trim(Label)//'.plt'
   !Freq_sum(:)=0 ;   nc=0 ;   p1=0. ;   p2=0.  ;  AvAmpl=0.
   Write(11,"(2x,I5,2x,A6,2x,I7,1x,F7.2,1x,I4,1x,A,I5)") MinAmp, Statn_ID2Mnem(STATION_ID), Ant_ID, sqrt(powr), &
         Filtring, LH, i_OddEven+1
   write(11,"(A)") '! i_chunk,  time,   sqrt(p),   (ChMax-ChMin)'
   !
   BckgrLev=MinAmp
   BckgrPow=MinPow*qp
   NBackgrInc=0
   NB=0
   FreqBackgr(:)=0.
   Do i_chunk=0,i_chunkMax
      If( MaxAmpChunk(i_chunk) .gt. 4500) cycle  !  gt than saturation due to too many zeroes.
      Dset_offset= i_chunk*Time_dim
      Call GetData(Chunk, Dset_offset, Time_dim) ! ListDataAtt should have been called first
      !
      RTime_s(:)=Chunk(:)*Hann(:)*AmpNrm
      Call RFTransform_CF_Filt(RTime_s,nu_fltr,SubSample_Offset,Cnu_s)
      Call RFTransform_CF2CT(Cnu_s, CTime_s )
      RTime_s(:)=REAL(CTime_s(:))
      ChMax=MaxVal(RTime_s)
      ChMin=MinVal(RTime_s)
      p=SUM(RTime_s(:)*RTime_s(:))*2./Time_dim ! factor 2 to account for the imaginary part
      !
      time=1000.d0*((i_chunk*Time_dim+SAMPLE_NUMBER_first)* Sample -DIPOLE_CALIBRATION_DELAY) ! in [ms]
      write(11,*) i_chunk, time, sqrt(p), (ChMax-ChMin) !,pw/p
   !  Diagnostics for amplitude-norm factor calculation
      If( (sqrt(p) .lt. 1.2) ) Then
         If(BckgrLev.lt. MaxAmpChunk(i_chunk)) Then
            BckgrLev=MaxAmpChunk(i_chunk)
            !write(2,*) '!amp',i_chunk, time, sqrt(p), Background(i_chunk), MaxAmpChunk(i_chunk), powerChunk(i_chunk)
         EndIf
         If(BckgrPow.lt. powerChunk(i_chunk)) Then
            BckgrPow=powerChunk(i_chunk)
            !write(2,*) '!pwr',i_chunk, time, sqrt(p), Background(i_chunk), MaxAmpChunk(i_chunk), powerChunk(i_chunk)
         EndIf
         NB=NB+1
         FreqBackgr(:)=FreqBackgr(:) + abs(Cnu_s(:))
         If(Background(i_chunk) ) NBackgrInc=NBackgrInc+1
      EndIf
   EndDo
   write(2,*) 'Max adc counts for background signals=',BckgrLev, ' q=',q*BckgrLev/MinAmp, BckgrPow, MinPow
   write(2,*) 'Tot # Chunks with sqrt(p)<1.2 =',NB, ' of these, # included in norm calc=',NBackgrInc, &
      ' = ',NBackgrInc*100./NBackgr,' %'
   Close(Unit=11)
   !
   If(NB.gt.1) Then
      FreqBackgr(:)=FreqBackgr(:)/NB
   Endif
   Call PlotFreqFilt(LH, Label, FreqBackgr)
   !
   flush(unit=2)
   Return
End Subroutine TestFilter
