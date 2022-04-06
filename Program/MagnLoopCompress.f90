!
PROGRAM Compress
  IMPLICIT none
  character*100 :: datafile
  Character*40 :: directory
  Character*8 :: extension
  real*8 :: t,ts,v,vs,dt,t_start,t_end,vr,vi
  integer :: i,Ns,nxx,j
  character(len=80) :: FileIn,FileOut
  Character(len=11) :: txt
  !
  Read(*,*) FileIn
  Read(*,*) FileOut
  !OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='Compress.out') !
    !OPEN(unit=28,FILE='loop_out.txt',FORM='FORMATTED',STATUS='OLD',ACTION='READ') ! space separated values
   dt=5d-6 ! [ms] ! convert from samples to milliseconds
 !  goto 1 ! LOFAR
    Ns=4 ! for 20 n sec average
   ! Ns=40 ! for .2 mu sec
   t=0.d0 ! 1215.138  ! [ms]
   j=0
   t_start=40 ; t_end=50
   t_start=0 ; t_end=168
   t_start=0 ; t_end=168 ! ML_213056_cleaned_4_i
   Read(*,*) Ns, t_start, t_end
   write(txt,"(I3.3,'-',I3.3,'-',I3.3)") Ns,NINT(t_start),NINT(t_end)
!  Magnetic Loop Broad-band antenna -------------------------------
    OPEN(unit=28,FILE='MagnLoop/'//TRIM(FileIn),FORM='FORMATTED',STATUS='OLD',ACTION='READ') ! space separated values
    !OPEN(unit=28,FILE='event_20190424_210306_cleaned.txt',FORM='FORMATTED',STATUS='OLD',ACTION='READ') ! space separated values
    !OPEN(unit=30,FILE='brian_20190424_213056_2_0-10.ssv',FORM='FORMATTED',STATUS='unknown') ! space separated values
    !OPEN(unit=30,FILE='brian_20190424_213056_2_35-40.ssv',FORM='FORMATTED',STATUS='unknown') ! space separated values
    !OPEN(unit=30,FILE='ML_213056_cleaned_40_all.ssv',FORM='FORMATTED',STATUS='unknown') ! space separated values
    OPEN(unit=30,FILE='MagnLoop/Compressed'//TRIM(FileOut)//txt//'.ssv',FORM='FORMATTED',STATUS='unknown') ! space separated values
    !OPEN(unit=30,FILE='ML_210306_cleaned_40_all.ssv',FORM='FORMATTED',STATUS='unknown') ! space separated values
    !OPEN(unit=30,FILE='ML_210306_cleaned_4_50.ssv',FORM='FORMATTED',STATUS='unknown') ! space separated values
    !OPEN(unit=30,FILE='ML_210306_cleaned_4_180.ssv',FORM='FORMATTED',STATUS='unknown') ! space separated values
     !
   !t_start=40 ; t_end=80  ! ML_210306_cleaned_4_50.ssv
   !t_start=120 ; t_end=168  ! ML_210306_cleaned_4_180.ssv
   do
       vs=0.
      Do i=1,Ns
         !read(28,*,iostat=nxx) v
         read(28,*,iostat=nxx) ts,v
         If(nxx.ne.0) goto 999
         !ts=ts+t
         vs=vs+v
      enddo
      t=dt*(j+0.5d0)*Ns  ! for ML to convert from sec to msec
      v=vs/Ns
      j=j+1
      if(t.lt.t_start) cycle
      write(30,"(1x,f11.6,3x,g13.6)") t,v ! for ML to print in ms
      If(t.gt.t_end) goto 999
    enddo
    goto 999
1  continue
    OPEN(unit=28,FILE='LOFAR_Time001-01.dat',FORM='FORMATTED',STATUS='OLD',ACTION='READ') ! space separated values
    !OPEN(unit=29,FILE='ML2019_1mus.dat',FORM='FORMATTED',STATUS='unknown')
    OPEN(unit=30,FILE='VHFpower331.ssv',FORM='FORMATTED',STATUS='unknown') ! space separated values
   Ns=100 ! for 500 n sec average
   do
      vs=0.
      ts=0.
      Do i=1,Ns
         read(28,*,iostat=nxx) j,vr,vi
         If(nxx.ne.0) stop
         !t=dt*i  ! for ML to convert from sec to msec
         ts=ts+dt*j
         vs=vs+vr*vr+vi*vi
      enddo
      v=vs/Ns
      t=ts/Ns
      write(30,"(1x,f11.6,3x,g13.6)") t,v ! for ML to print in ms
      !write(30,"(1x,f11.2,3x,F8.1)") t,v ! for VHF to print in mus
    enddo
999 continue
    write(*,*) i,', number of lines written=',j,', last time written=',t,' [ms]',', last record read',ts
    stop
!--------------------------------------------------
  END PROGRAM Compress
