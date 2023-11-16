Module GLEplots
   use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows
   use DataConstants, only : OutFileLabel, DataFolder
   use DataConstants, only : RunMode ! 1=Explore; 2=Calibrate; 3=Impulsive imager; 4=Interferometry
   implicit none
Contains
Subroutine GLEplotControl(PlotType, PlotName, PlotDataFile, SpecialCmnd, Submit, Bckgr)
   character(LEN=*), intent(in), optional :: PlotType, PlotName, PlotDataFile, SpecialCmnd, Bckgr
   Logical, intent(in), optional :: Submit
   Character(len=40), save :: BatchFile=''
   character*100 :: shellin
   !Character(len=19), parameter :: CallGLE='call GLE -d pdf -o '  ! space at end is important
   !Character(len=5), parameter :: PJ='.pdf ' ! space at end is important
   Character(len=26), parameter :: CallGLE='call GLE -d jpg -r 300 -o '  ! space at end is important
   Character(len=5), parameter :: PJ='.jpg ' ! space at end is important
   !
   !inquire(unit=10, opened=itsopen)  ; if ( itsopen ) then okay  ; logical itsopen
   If(BatchFile.eq. '' .and. (present(PlotType) .or. present(SpecialCmnd))) Then  ! Open unit 10
      SELECT CASE (MOD(RunMode,10))  ! 1=Explore; 2=Calibrate; 3=Impulsive imager; 4=Interferometry
         CASE(0)  ! RFI-Mitigation
            BatchFile='Afig-RFIM'//trim(OutFileLabel)  ! Automatic Run of GLE -
         CASE(1)  ! Explore
            BatchFile='Afig-Expl'//trim(OutFileLabel)  ! Automatic Run of GLE -
         CASE(2)  ! Calibrate
            BatchFile='Afig-Cal'//trim(OutFileLabel)
         CASE(3)  ! Impulsive Imager
            BatchFile='Afig-ImpIm'//trim(OutFileLabel)
         CASE(4)  ! Interf Imager TRID
            BatchFile='Afig-Intf'//TRIM(OutFileLabel)
         CASE(5)  ! TrackScatter; Impulsive imager
            BatchFile='AImp'//trim(OutFileLabel)
         CASE(6)  ! InterfSelect; TRID imager
            BatchFile='ATRID'//trim(OutFileLabel)
         CASE(7)  ! InterfSelect; TRID imager
            BatchFile='Afig-FldCal'//trim(OutFileLabel)
         CASE(8)  ! PeakInterferometry
            BatchFile='Afig-PkInt'//trim(OutFileLabel)
         CASE(9)  ! PeakInterferometry
            BatchFile='Afig-MDD'//trim(OutFileLabel)
         CASE DEFAULT
            write(2,*) 'not a foreseen plotting RunMode:',RunMode
            write(*,*) 'not a foreseen plotting RunMode:',RunMode
            stop
      End SELECT
      If(windows) then
         OPEN(UNIT=10,STATUS='REPLACE',ACTION='WRITE',FILE='RunGLEplots.bat')
         !Write(10,"(7(A,/) )") 'call ../ShortCuts.bat'  ! define the proper shortcuts; not needed here
      Else
         Open(UNIT=10,STATUS='unknown',ACTION='WRITE',FILE=TRIM(BatchFile)//'.sh')   ! It is assumed that this scipt is launched in the FlashFolder
         Write(10,"(7(A,/) )") &
            '#!/bin/bash -v', &
            '#','#' , 'source  ${LL_BaseDir}/ShortCuts.sh'     ! ,  'FlashFolder=$(pwd)',  'ProgramDir="/home/olaf/LMA-fit/LMA2019/program"'
      EndIf
   EndIf
   If(BatchFile.eq. '') Return
   !
   If(present(PlotType)) Then  ! PlotName, PlotDataFile
            !PlotDataFile='files/AmplFit'//trim(datafile(1))//TRIM(OutFileLabel)//TRIM(extension)
            !PlotName=trim(datafile(1))//TRIM(extension)//TRIM(OutFileLabel)//'AmplFit'
            If(Windows) Then
               If(present(Bckgr)) Then
                  Write(10,"(A)") &
                  CallGLE//trim(PlotName)//PJ//TRIM(UtilitiesFolder)//trim(PlotType)//".gle "&
                        //'"'//TRIM(FlashFolder)//trim(PlotDataFile)//'"'//' "'//TRIM(FlashFolder)//trim(Bckgr)//'"'
               Else
                  Write(10,"(A)") &
                  CallGLE//trim(PlotName)//PJ//TRIM(UtilitiesFolder)//trim(PlotType)//".gle "&
                        //'"'//TRIM(FlashFolder)//trim(PlotDataFile)//'"'
                  !"call GLE -d pdf -o "//trim(PlotName)//".pdf "//TRIM(UtilitiesFolder)//trim(PlotType)//".gle "&
                  !      //'"'//TRIM(FlashFolder)//trim(PlotDataFile)//'"'
               EndIf
            Else
               If(present(Bckgr)) Then
                  Write(10,"(A)") &
                  "gle -d pdf -o "//trim(PlotName)//".pdf "//TRIM(UtilitiesFolder)//trim(PlotType)//".gle "&
                        //TRIM(FlashFolder)//TRIM(PlotDataFile)//' '//TRIM(FlashFolder)//trim(Bckgr) !//'"'
               Else
                  Write(10,"(A)") &
                  "gle -d pdf -o "//trim(PlotName)//".pdf "//TRIM(UtilitiesFolder)//trim(PlotType)//".gle "&
                        //TRIM(FlashFolder)//TRIM(PlotDataFile)
               EndIf
               If((PlotType.eq.'InterfContour') .or. (PlotType.eq.'EIContour')) then
                  write(10,"('rm ',A,'{*.z,*-cvalues.dat,*-clabels.dat,*-cdata.dat}')") TRIM(DataFolder)
                  write(10,"('rm ',A,'{*SpecWin_*.csv,*IntfTrack_*.csv,*Interferometer*.csv,*_EISpec*.csv}')") &
                              TRIM(DataFolder)//trim(OutFileLabel)
               EndIf
            EndIf
   EndIf
   !
   If(present(SpecialCmnd)) Then  ! PlotName, PlotDataFile
            If(.not. Windows) Then
               write(10,"(A)") trim(SpecialCmnd)
            EndIf
   EndIf
   !
   If(present(Submit)) Then
      If(Submit) Then
         If(windows) then
            Close(unit=10)
            shellin = 'RunGLEplots.bat'
            CALL system(shellin)
            !write(2,*) 'system command submitted'
            Write(2,*) 'GLE should be running now'
         Else
            If(RunMode.eq.2) Then
              write(10,"('rm ',A,'{LOFAR_Corr*.dat,CCPeakPhase*.dat}')") TRIM(DataFolder)
      !!         write(10,"('rm ',A,'{LOFAR_Time*.dat,LOFAR_Freq*.dat}')") TRIM(DataFolder)
            EndIf
            Close(unit=10)
            shellin = 'chmod 755 '//TRIM(BatchFile)//'.sh'
            CALL system(shellin)
            shellin = 'nohup ./'//TRIM(BatchFile)//'.sh  >'//TRIM(BatchFile)//'.log 2>&1  & '
            CALL system(shellin)
            Write(2,*) TRIM(BatchFile)//'.sh  is submitted'
            Write(*,*) TRIM(BatchFile)//'.sh  is submitted'
         EndIf
         BatchFile=''
      EndIf
   EndIf
   Return
End Subroutine GLEplotControl

End Module GLEplots
