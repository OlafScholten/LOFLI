Subroutine System_Initiation(Utility, release, ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows)
   !use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows
   Implicit none
   Character(len=*), intent(in) :: Utility, release
   Character(len=*), intent(out) :: ProgramFolder, UtilitiesFolder, FlashFolder, FlashName
   Logical, intent(out) :: Windows
   Integer :: i,j,i_guess
   character(len=1) :: FChar='/'
   INTEGER :: DATE_T(8)
   !
   write(2,"(3x,5(1H-),1x,'Utility: ',A,', release ',A20,25(1H-))") TRIM(Utility), release
   CALL DATE_AND_TIME (Values=DATE_T)
   WRITE(2,"(3X,5(1H-),1x,'run on ',I2,'/',I2,'/',I4,' , started at ',&
       I2,':',I2,':',I2,'.',I3,1X,25(1H-))") &
       DATE_T(3),DATE_T(2),DATE_T(1),(DATE_T(i),i=5,8)
   !
   ! EXECUTE_COMMAND_LINE:	  	Execute a shell command
   ! GETCWD:	  	Get current working directory
   ! GETENV:	  	Get an environmental variable
   !CALL GET_ENVIRONMENT_VARIABLE(NAME[, VALUE, LENGTH, STATUS, TRIM_NAME)
   CALL GET_ENVIRONMENT_VARIABLE("OS", FlashName)
   If(FlashName(1:7).eq.'Windows') Windows=.true.
   If(Windows) FChar='\'
   !
   Call GETCWD(FlashFolder)  ! Get current working directory
   j=LEN_TRIM(FlashFolder)
   write(2,*) 'FlashFolder="',TRIM(FlashFolder),'" with',j, ' characters'
   flush(unit=2)
   i_guess=0
   Do i=1,j-1
      If(FlashFolder(i:i) .eq. FChar) Then
         i_guess=i ! search for the bottom-most folder
         cycle
      Endif
   Enddo
   If(j.eq.0) then
      FlashName='Anonymous'   !  set equal to the name of the folder
   Else
      FlashName=Trim(FlashFolder(i_guess+1:j))   !  set equal to the name of the folder
   Endif
   !
   If(Windows) Then
      ProgramFolder='%ProgramDir%'
      UtilitiesFolder='%UtilDir%'
      FlashFolder='%FlashFolder%'//TRIM(FChar)  ! should be equivalent in usage to the original FlashName (with ending /)
   Else
      call System_MemUsage(j)  !  j is just dummy here
      ProgramFolder='${ProgramDir}'
      UtilitiesFolder='${UtilDir}'
      FlashFolder='${FlashFolder}'//'/'  ! '../'//TRIM(FlashName)//'/' is shorter, but should be equivalent in usage to the original FlashName (with ending /)
   EndIf
   Write(2,*) 'FlashName: "',TRIM(FlashName) ! ,'", FlashFolder: "',TRIM(FlashFolder),'"'
   !
   Return
End Subroutine System_Initiation
!===============================================
Subroutine System_MemUsage(valueRSS)
!          * VmPeak: Peak virtual memory size.
!          * VmSize: Virtual memory size.   total program size
!  VMsize is the "address space" that the process has in use: the number of available adresses. These addresses do not have to have any physical memory attached to them. ( Attached physical memory is the RSS figure)
!  VmSize includes RSS, plus things like shared libraries and memory mapped files (which don't actually use RAM), as well as allocated, but unused, memory.
!          * VmLck: Locked memory size (see mlock(3)).
!          * VmHWM: Peak resident set size ("high water mark").
!          * VmRSS: Resident set size.
!  VmRSS is the measure of how much RAM the process is actually using.
!          * VmData, VmStk, VmExe: Size of data, stack, and text segments.
   implicit none
   integer, intent(out) :: valueRSS
   character(len=30):: filename=' '
   character(len=30) :: line
   character(len=8)  :: pid_char=' '
   integer :: pid
   logical :: ifxst
   !
   valueRSS=-1    ! return negative number if not found
   !
   !--- get process ID
   pid=getpid()
   write(pid_char,'(I8)') pid
   filename='/proc/'//trim(adjustl(pid_char))//'/status'
   !
   !--- read system file
   inquire (file=filename,exist=ifxst)
   if (.not.ifxst) then
     write (*,*) 'system file does not exist'
     return
   endif
   open(unit=100, file=filename, action='read')
   do
     read (100,'(a)',end=120) line
     if (line(1:6).eq.'VmRSS:') then
        read (line(7:),*) valueRSS
        exit
     endif
   enddo
120 continue
   close(100)
   !
   Write(2,*) 'pid=', pid,', mem usage=',valueRSS,' kB'
   !
   return
end Subroutine System_MemUsage
!======================
