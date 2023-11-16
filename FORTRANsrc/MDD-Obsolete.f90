!-----------------------------------------------
Subroutine MDD_IniSrcs()         !  Obsolete
   use constants, only : dp, sample, c_mps
   Use MDD_Pars, only : NSrc_max, IRFW_s, T_Range  ! Input
   Use MDD_Pars, only : NSrc, dCenT, dCenloc  ! Output
   Implicit none
   Real(dp) :: dt, dd, D
   Integer :: i
   !
   NSrc_max=17
   Allocate(  dCenT(1:NSrc_max) ) ! time difference (in [samples]) of each DD from central time  C_Tms
   Allocate( dCenloc(1:3,1:NSrc_max) ) ! location of each DD w.r.t. C_Loc,  in [m]
   dt=IRFW_s*4  ! [samples]
   dd=dt*c_mps*sample  ! [samples] * [m/s] * [s/sample] =[m]
   dCenT(:)=0.
   dCenloc(1:3,:)=0.
   dCenT(2)=0. ;  dCenloc(1,2)=+dd
   dCenT(3)=0. ;  dCenloc(1,3)=-dd
   dCenT(4)=0. ;  dCenloc(2,4)=+dd
   dCenT(5)=0. ;  dCenloc(2,5)=-dd
   dCenT(6)=0. ;  dCenloc(3,6)=+dd
   dCenT(7)=0. ;  dCenloc(3,7)=-dd
   dCenT(8)=+dt
   dCenT(9)=-dt
   If(T_Range/2. .lt. sqrt(2.)*dt) Then
      Nsrc=9
      Return
   EndIf
   dCenT(10)=+dt ;  dCenloc(1,10)=+dd
   dCenT(11)=+dt ;  dCenloc(1,11)=-dd
   dCenT(12)=-dt ;  dCenloc(1,12)=+dd
   dCenT(13)=-dt ;  dCenloc(1,13)=-dd
   dCenT(14)=+dt ;  dCenloc(2,14)=dd
   dCenT(15)=+dt ;  dCenloc(2,15)=-dd
   dCenT(16)=-dt ;  dCenloc(2,16)=dd
   dCenT(17)=-dt ;  dCenloc(2,17)=-dd
   Nsrc=17
   Return
End Subroutine MDD_IniSrcs
!-----------------------------------------------
Subroutine MDD_ResortSrcs()
!  Keep the NSrc_max brightest sources
   use constants, only : dp, sample, c_mps
   Use MDD_Pars, only : NSrc_max, IRFW_s, T_Range, Bias_inv_base, NSrc_min  ! Input
   Use MDD_Pars, only : NSrc, dCenT, dCenloc, Bias_inv  ! Output
   Use MDD_Pars, only : DDChiSQ, Weight
   Use MDD_Pars, only : DelDip
   Use unque, only : Double_RI_sort  ! (N,Real(dp),integer)
   Implicit none
   !Real(dp) :: dt, dd, D
   Real(dp), allocatable :: Temp(:), dtsort(:), dlsort(:,:)
   Complex(dp), allocatable :: ddsort(:,:)
   Integer, allocatable :: Permut(:)
   Integer :: m,n !i,
   !
   NSrc_max=NSrc_min
   Bias_inv=Bias_inv_base
   If(NSrc .le. NSrc_max) goto 9  ! no readj needs to be done
   Allocate( Temp(1:NSrc), Permut(1:NSrc) )
   Allocate( dtsort(1:NSrc_max), dlsort(1:3,1:NSrc_max), ddsort(1:3,1:NSrc_max) )
   Do m=1,NSrc
      Permut(m)=m
      Temp(m)=SUM( ABS( DelDip(:,m) )**2 )
   EndDo
   Call Double_RI_sort(NSrc,Temp,Permut)
!
   Write(2,"(A,1pG10.2,A,1pG10.2)") 'Rastered sources, selected range from I=', Temp(NSrc),' till', Temp(NSrc-NSrc_max+1)
   Write(*,"(A,1pG10.2,A,1pG10.2)") 'Rastered selected I=', Temp(NSrc),' till', Temp(NSrc-NSrc_max+1)
   Do m=1,NSrc_max
      n=Permut(NSrc-m+1)
      dtsort(m)=dCenT(n)
      dlsort(:,m)=dCenloc(:,n)
      ddsort(:,m)=DelDip(:,n)
   EndDo ! m=1,NSrc_max
   DeAllocate( dCenT, dCenloc, DelDip )
   Allocate(  dCenT(1:NSrc_max) ) ! time difference (in [samples]) of each DD from central time  C_Tms
   Allocate( dCenloc(1:3,1:NSrc_max) ) ! location of each DD w.r.t. C_Loc,  in [m]
   Allocate( DelDip(1:3,1:NSrc_max) ) ! Dipole moment of each DD
   dCenT(:)=     dtsort(:)
   dCenloc(:,:)= dlsort(:,:)
   DelDip(:,:)=ddsort(:,:)
   DeAllocate( Temp, Permut, dtsort, dlsort, ddsort )
   !
   Nsrc=NSrc_max
   DDChiSQ(:,:)=0.
   !Call MDD_write()
9  Return
End Subroutine MDD_ResortSrcs
