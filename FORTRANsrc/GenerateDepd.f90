      Module LFS
    implicit none
   Character(len=130) :: line
   Character(len=40) :: IncFile
   Character(len=7) :: StIncl='Include'
   Character(len=1) :: StExcl='!'
   Character(len=3) :: StEnd='End'
   Character(len=10) :: StSubr='Subroutine'
   Character(len=8) :: StFunc='Function'
   Character(len=7) :: StProg='Program'
   Character(len=6) :: StModu='Module'
   Character(len=16) :: Lead='!'
   Integer, Parameter :: zeroOff=2
   Integer :: nxx=0, Indnt=0
   Integer :: PosExcl, PosEnd, PosIncl, Pos1, Pos2, pos
   Logical :: Prnt
      contains
! +++++++++++++++++
Recursive Subroutine SubFun(Source, Unt)
   implicit none
   Character(len=30), intent(in) :: source
   Integer, intent(in) :: Unt
   Integer :: Unt2
   Logical :: Found
   !
   nxx=0
   Open(unit=Unt,STATUS='unknown',ACTION='read', FILE = trim(Source), IOStat=nxx)
   !write(12,*) 'unt',Unt,nxx,trim(Source)
   If(nxx.ne.0) return
   Do
      line='!'
      read(Unt,"(A130)",IOSTAT=nxx)  line
      !If(line.eq.' ') cycle
      line = '   '//ADJUSTL(line)
      !write(*,*) line
      If(nxx.gt.0) then
         exit
      Endif
      PosExcl = Scan(line, StExcl )
      If(PosExcl.eq.0) PosExcl=140
      PosIncl = INDEX(line, StIncl )
      If(PosIncl.eq.0) PosIncl=140
      If(PosIncl.ge.PosExcl) then
         PosIncl=140
      Else
         Write(2,*) PosIncl,PosExcl
         Pos1 = INDEX(line, "'" )+1
         Pos2 = INDEX(line, "'" ,.true.)-1
         !Write(*,*) 'i',Pos1,Pos2
         IncFile=line(Pos1:Pos2)
         Write(2,*) '--',trim(IncFile),'---'
         Write(12,*) trim(Lead),trim(line)
         Write(13,"(A,1x)", ADVANCE='NO') trim(IncFile)  ! List of dependencies
         Indnt=Indnt+1
         !write(12,*) '+zeroOff-I',zeroOff+(Indnt-1)*4
         Lead(zeroOff+(Indnt-1)*4:zeroOff+Indnt*4)="   I"
         Unt2=Unt+3
         Call SubFun(IncFile, Unt2)
         !write(12,*) '-zeroOff-I',zeroOff+(Indnt-1)*4
         Lead(zeroOff+(Indnt-1)*4:zeroOff+Indnt*4)="    "
         Indnt=Indnt-1
         cycle
      Endif
      Call Recognize(StProg,Found)
      If(Found) Cycle
      !Pos = INDEX(line, StProg )
      !If(Pos.gt.0 .and. PosExcl.gt.Pos) then
      !   Write(*,*) 'Program:',pos
      !   PosEnd = INDEX(line, StEnd )
      !   If(PosEnd.eq.0) PosEnd=140
      !   If(PosEnd.gt.Pos) Write(12,*) trim(Lead),line
      !   Cycle
      !Endif
      Call Recognize(StSubr,Found)
      If(Found) Cycle
      Call Recognize(StFunc,Found)
      If(Found) Cycle
      Pos = INDEX(line, StModu )
      If(Pos.gt.0 .and. PosExcl.gt.Pos) then
         Write(2,*) 'Module:',pos
         PosEnd = INDEX(line, StEnd )
         If(PosEnd.eq.0) PosEnd=140
         If(PosEnd.gt.Pos) then
            Write(12,*) trim(Lead),line
            Indnt=Indnt+1
            !write(12,*) '+zeroOff-M',zeroOff+(Indnt-1)*4
            Lead(zeroOff+(Indnt-1)*4:zeroOff+Indnt*4)="   M"
         Else
            !write(12,*) '-zeroOff-M',zeroOff+(Indnt-1)*4,line
            Lead(zeroOff+(Indnt-1)*4:zeroOff+Indnt*4)="    "
            Indnt=Indnt-1
         Endif
      Endif
   Enddo
   Close(Unit=Unt)
   Return
   End Subroutine SubFun
!==================
Subroutine Recognize(St,Found) ! search for string 'ST' in 'line' and return position if found
! In case found and before an '!' (position given by PosExcl) and there is no 'StEnd', then write line to unit=12 and set pos
   Implicit none
   Character(len=*), intent(in) :: St
   Logical, intent(out) :: Found
   Found=.false.
   Pos = INDEX(line, St )
   If(Pos.gt.0 .and. PosExcl.gt.Pos) then
      !Write(*,*) 'Program:',pos
      PosEnd = INDEX(line, StEnd )
      If(PosEnd.eq.0) PosEnd=140
      If(PosEnd.gt.Pos) Write(12,*) trim(Lead),line
      Found=.true.
   EndIf
   Return
End Subroutine Recognize
!=================
End  Module LFS
!=================
Program ListFortStructure
!• ADJUSTL:	  	Left adjust a string
!• ADJUSTR:	  	Right adjust a string
! • INDEX:	  	Position of a substring within a string
!• LNBLNK:	  	Index of the last non-blank character in a string
!• REPEAT:	  	Repeated string concatenation
!• SCAN:	  	Scan a string for the presence of a set of characters
!• VERIFY:	  	Scan a string for the absence of a set of characters
   use LFS, only : SubFun
   Implicit none
   Character(len=30) :: source,Date,NoExt
   Integer :: Unt=11,i
   INTEGER :: DATE_T(8)
   !
   i= IARGC()
   Source = 'LOFAR-Imag-v8.f90'
   If(i.ge.1) Then
      CALL getarg(1, Source)
   endif
   OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='LFS.out')
   write(2,*) 'source code =',trim(Source)
   NoExt='?'
   Do i=1,30
      If(Source(i:i).eq.'.') exit
      NoExt(i:i)=Source(i:i)
   EndDo
   CALL DATE_AND_TIME (Values=DATE_T)
   !   ----- run on 16/ 2/2022 , started at 12:17:22.719 -------------------------
   WRITE(Date,"(I2,'/',I2.2,'/',I4,'@',I2,':',I2.2,':',I2.2,'.',I3.3)") &
         DATE_T(3),DATE_T(2),DATE_T(1),(DATE_T(i),i=5,8)
   OPEN(UNIT=13,STATUS='unknown',ACTION='WRITE',FILE=trim(NoExt)//'_depd.sh')
   Write(13,"(A)", ADVANCE='NO') 'export MainProgDepd="'
   OPEN(UNIT=12,STATUS='unknown',ACTION='WRITE',FILE=trim(NoExt)//'.LFS')
   !Open(unit=11,STATUS='unknown',ACTION='read', FILE = trim(Source))
   write(12,*) '!-----------------',Trim(Date),'--------------------------'
   Write(12,*) '!------- Source file:',trim(Source),' ---- '
   Call SubFun(Source,Unt)
   Write(13,"(A)") '"'
   write(12,*) '!-----------------',Trim(Date),'--------------------------'
   !
   Stop
End Program ListFortStructure
!===============
