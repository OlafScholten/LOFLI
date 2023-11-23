    Include 'ConstantsModules.f90'
    !Include 'ParamModules.f90'
    !Include 'MappingUtilities.f90'
!-----------------------------------
!=================================
!-----------------------------------
Program AnalyzeLeader
! Read-in .dat plot files from sequences of TRI-D images along the path of a recoil leader
! compare the source densities at different times for multiple sections of this track
! Slice in Northing coordinate since that applies to the case of interest
!
   use constants, only : dp
   Implicit none
   !
   Integer, parameter :: N_data_max=999
   Integer, parameter :: N_slices_max=50
   Integer, parameter :: N_files_max=9
   Real(dp) :: time(N_data_max),East(N_data_max),North(N_data_max),heigh(N_data_max)
   Integer :: Numb(N_slices_max,N_files_max)
   Real(dp) :: Step, E_i,E_f,N_i,N_f
   Integer :: Lab, i,j,m, Ndata, N_slices, nxx, N_files
   Character(len=85) :: datFile
   !
   Open(unit=2,STATUS='unknown',ACTION='write', FILE ="AnaLead.out")
   !
   !Note: it is assumed that the order of the station is the same for the two calibration files

   Read(*,*) Step
   N_files=N_files_max
   Do m=1,N_files_max
      Read(*,*,iostat=nxx) datFile
      If(nxx.ne.0) then
         N_files=m
         exit
      EndIf
      write(2,"('! data',I2,' =',A)") m,TRIM(datFile)
      Open(unit=10,STATUS='old',ACTION='read', FILE ="files/IntfSpecSel"//TRIM(datFile)//"Mx_d.dat")
      Read(10,*) E_i,E_f,N_i,N_f
      Read(10,*) lab
      Do i=1,N_data_max
         Read(10,*,iostat=nxx) lab,time(i),East(i),North(i),heigh(i)
         If(nxx.ne.0) exit
         Ndata=i
      Enddo
      !
      Numb(:,m)=0
      N_slices=1+(N_f-N_i)/Step
      Do i=1,Ndata
         j=1+(North(i)-N_i)/Step
         Numb(j,m)=Numb(j,m)+1
      EndDo
   EndDo
   !
   write(2,"('!',1x,'[km], N=',i4,15i6)") ((i),i=1,N_files)
   Do j=1,N_slices
      write(2,"(1x,f6.2, 1x ,15i6)") N_i+(j-1)*step,Numb(j,1:N_files)
   EndDo
    !
!         write(2,"(A6)", ADVANCE='no') Fine_STMnm(i)
   write(2,*) '  '
   !
   Stop
   !
End Program AnalyzeLeader
