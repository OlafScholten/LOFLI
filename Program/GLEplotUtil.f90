Module GLEplots
   use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows
   use DataConstants, only : OutFileLabel, DataFolder
   use DataConstants, only : RunMode ! 1=Explore; 2=Calibrate; 3=Impulsive imager; 4=Interferometry
   implicit none
Contains
Subroutine GLEplotControl(PlotType, PlotName, PlotDataFile, SpecialCmnd, Submit)
   character(LEN=*), intent(in), optional :: PlotType, PlotName, PlotDataFile, SpecialCmnd
   Logical, intent(in), optional :: Submit
   Character(len=40), save :: BatchFile=''
   character*100 :: shellin
   !
   !inquire(unit=10, opened=itsopen)  ; if ( itsopen ) then okay  ; logical itsopen
   If(BatchFile.eq. '' .and. (present(PlotType) .or. present(SpecialCmnd))) Then  ! Open unit 10
      SELECT CASE (RunMode)  ! 1=Explore; 2=Calibrate; 3=Impulsive imager; 4=Interferometry
         CASE(0)  ! RFI-Mitigation
            BatchFile='Afig-RFIM'  ! Automatic Run of GLE -
         CASE(1)  ! Explore
            BatchFile='Afig-Expl'  ! Automatic Run of GLE -
         CASE(2)  ! Calibrate
            BatchFile='Afig-Cal'
         CASE(3)  ! Explore
            BatchFile='Afig-ImpIm'
         CASE(4)  ! Interf Imager TRID
            BatchFile='Afig-Intf'//TRIM(OutFileLabel)
         CASE(5)  ! TrackScatter; Impulsive imager
            BatchFile='AImp'
         CASE(6)  ! Select; TRID imager
            BatchFile='ATRID'
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
            '#','#' , 'source ../ShortCuts.sh'     ! ,  'FlashFolder=$(pwd)',  'ProgramDir="/home/olaf/LMA-fit/LMA2019/program"'
      EndIf
   EndIf
   If(BatchFile.eq. '') Return
   !
   If(present(PlotType)) Then  ! PlotName, PlotDataFile
            !PlotDataFile='files/AmplFit'//trim(datafile(1))//TRIM(OutFileLabel)//TRIM(extension)
            !PlotName=trim(datafile(1))//TRIM(extension)//TRIM(OutFileLabel)//'AmplFit'
            If(Windows) Then
               Write(10,"(A)") &
                  "call GLE -d pdf -o "//trim(PlotName)//".pdf "//TRIM(UtilitiesFolder)//trim(PlotType)//".gle "&
                        //'"'//TRIM(FlashFolder)//trim(PlotDataFile)//'"'
            Else
               write(10,"(A)") &
                  "gle -d pdf -o "//trim(PlotName)//".pdf "//TRIM(UtilitiesFolder)//trim(PlotType)//".gle "&
                        //TRIM(FlashFolder)//TRIM(PlotDataFile)
               If((PlotType.eq.'InterfContour') .or. (PlotType.eq.'EIContour')) then
                  write(10,"('rm ',A,'{*.z,*-cvalues.dat,*-clabels.dat,*-cdata.dat}')") TRIM(DataFolder)
                  write(10,"('rm ',A,'{*SpecWin_*.csv,*IntfTrack_*.csv,*Interferometer*.csv,*_EISpec*.csv}')") TRIM(DataFolder)
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
            Write(2,*) 'GLE should be running now'
         Else
            If(RunMode.eq.2) Then
               write(10,"('rm ',A,'{LOFAR_Corr*.dat,CCPeakPhase*.dat}')") TRIM(DataFolder)
               write(10,"('rm ',A,'{LOFAR_Time*.dat,LOFAR_Freq*.dat}')") TRIM(DataFolder)
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
