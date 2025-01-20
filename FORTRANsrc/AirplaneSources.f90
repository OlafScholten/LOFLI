!=================================
Module ReadPlane
Contains
Subroutine ReadAirplane(OutFileLabel, BackgrndData, PeakNrTotal, SourceTime_ms, SourcePos, SourceIntensity, PeakWidth, &
         Chi2, SrcQual, Stk_NEh, N_backgr)
   !
   use constants, only : dp
   use DataConstants, only : DataFolder
   !
   Implicit none
   Integer, intent(out) :: PeakNrTotal, N_backgr
   Real(dp), allocatable, intent(out) :: SourceTime_ms(:), SourcePos(:,:), SourceIntensity(:), Chi2(:)
   Complex, allocatable, intent(out) :: Stk_NEh(:,:,:)
   Integer, allocatable, intent(out) :: PeakWidth(:), SrcQual(:,:)
   Character(len=*), intent(in) :: OutFileLabel, BackgrndData
   Character(len=1) :: Marker
   integer :: i_Peak, nxx, DataUnit, i, ip, j
   Complex :: Stk(1:6)
   Logical :: exists
   real :: time, I3,Un,Lin,Circ,Zen,Azi
   !
   DataUnit=28
   INQUIRE(FILE = trim(DataFolder)//TRIM(OutFileLabel)//'IntfPkSources.dat', exist=exists)  ! in the main Fsash directory
   If(exists) then
      OPEN(UNIT=DataUnit,STATUS='old',ACTION='read',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'IntfPkSources.dat', IOSTAT=nxx)
      If(nxx.ne.0) Then
         write(2,*) 'file not found: "',trim(DataFolder)//TRIM(OutFileLabel)//'IntfPkSources.dat','"'
         stop 'no datafile'
      EndIf
      read(DataUnit,*) PeakNrTotal
      Allocate(SourceTime_ms(-7:PeakNrTotal), SourcePos(1:3,-7:PeakNrTotal), SourceIntensity(-7:PeakNrTotal), &
            PeakWidth(-7:PeakNrTotal), Chi2(-7:PeakNrTotal), SrcQual(1:2,-7:PeakNrTotal), Stk_NEh(1:3,1:3,-7:PeakNrTotal) )
      !
      ip=PeakNrTotal
      Do i_Peak=1,ip
         Read(DataUnit,*,iostat=nxx) i, SourceTime_ms(i_Peak), SourcePos(1:3,i_Peak), SourceIntensity(i_Peak), PeakWidth(i_peak), &
            Chi2(i_Peak), SrcQual(1:2,i_Peak), Stk_NEh(1,1:3,i_Peak), Stk_NEh(2,2:3,i_Peak), Stk_NEh(3,3,i_Peak)
         !write(2,*) i_Peak, SourceTime_ms(i_Peak), SourcePos(1:3,i_Peak)
         Stk_NEh(2,1,i_Peak)=conjg(Stk_NEh(1,2,i_Peak))
         Stk_NEh(3,1,i_Peak)=conjg(Stk_NEh(1,3,i_Peak))
         Stk_NEh(3,2,i_Peak)=conjg(Stk_NEh(2,3,i_Peak))
         If(nxx.ne.0) then
            PeakNrTotal=i_peak-1
            exit
         EndIf
      EndDo
   Else
      OPEN(unit=DataUnit,FILE=trim(DataFolder)//TRIM(OutFileLabel)//'PkInt.csv', STATUS='OLD',ACTION='READ', IOSTAT=nxx) ! space separated values
      If(nxx.ne.0) then   ! check weather this particular file can be opened
         Write(2,*) 'Probelems opening file:"',trim(DataFolder)//TRIM(OutFileLabel)//'PkInt.csv','"'
         stop 'no datafile'
      endif
      read(DataUnit,*,iostat=nxx) Marker, PeakNrTotal ! skip first comment line
      Allocate(SourceTime_ms(-7:PeakNrTotal), SourcePos(1:3,-7:PeakNrTotal), SourceIntensity(-7:PeakNrTotal), &
            PeakWidth(-7:PeakNrTotal), Chi2(-7:PeakNrTotal), SrcQual(1:2,-7:PeakNrTotal), Stk_NEh(1:3,1:3,-7:PeakNrTotal) )
      ip=PeakNrTotal
      i_peak=1
      Do j=1,ip
         nxx=0
         read(DataUnit,*,iostat=nxx) Marker, i,SourceTime_ms(i_Peak), SourcePos(1:3,i_Peak), PeakWidth(i_peak), Chi2(i_Peak), &
             SourceIntensity(i_Peak), I3,Un,Lin,Circ,Zen,Azi, Stk_NEh(1,1:3,i_Peak), Stk_NEh(2,2:3,i_Peak), Stk_NEh(3,3,i_Peak)
         If(nxx.gt.0) then
            write(2,*) 'Read error:',i
         ElseIf (nxx.lt.0) Then
            Write(2,*) 'EOF reached'
            exit
         Endif
         If(marker.eq.'!') Then
            cycle
         ElseIf(marker.ne.'S') Then
            write(2,*) 'problem: ',Marker, i
            cycle
         EndIf
         !
         Stk_NEh(2,1,i_Peak)=conjg(Stk_NEh(1,2,i_Peak))
         Stk_NEh(3,1,i_Peak)=conjg(Stk_NEh(1,3,i_Peak))
         Stk_NEh(3,2,i_Peak)=conjg(Stk_NEh(2,3,i_Peak))
         SourcePos(1:3,i_Peak)=SourcePos(1:3,i_Peak)*1000.
         PeakNrTotal=i_peak
         i_peak=i_peak+1
      enddo     ! loop over all source points in this file
   EndIf
   Close(Unit=DataUnit)
   Write(2,*) 'analysis done for', PeakNrTotal
   Write(*,*) 'analysis done for', PeakNrTotal
   !
   Stk_NEh(1:3,1:3,0)=0.
   N_backgr=0
   If(BackgrndData.ne."") Then
      OPEN(UNIT=DataUnit,STATUS='old',ACTION='read',FILE=trim(DataFolder)//TRIM(BackgrndData), IOSTAT=nxx)
      If(nxx.ne.0) Then
         write(2,*) 'file not found: "',trim(DataFolder)//TRIM(BackgrndData),'"'
         write(2,*) 'no background study included'
         return
      EndIf
      read(DataUnit,*) marker
      Do
         Read(DataUnit,*,iostat=nxx) i, ip, time, Stk(1:6)
         If(nxx.ne.0) exit
         Stk_NEh(1,1:3,0)=Stk_NEh(1,1:3,0)+ Stk(1:3)
         Stk_NEh(2,2:3,0)=Stk_NEh(2,2:3,0)+ Stk(4:5)
         Stk_NEh(3,3,0)=  Stk_NEh(3,3,0)  + Stk(6)
         N_backgr=N_backgr +1
      Enddo
      Stk_NEh(2,1,0)=conjg(Stk_NEh(1,2,0))
      Stk_NEh(3,1,0)=conjg(Stk_NEh(1,3,0))
      Stk_NEh(3,2,0)=conjg(Stk_NEh(2,3,0))
   EndIf
   !
   Return
End Subroutine ReadAirplane
End Module ReadPlane
!================================
Program AirplaneSources
!   Chi-square interferometric imaging of peaks found by the impulsive imager
!
!   Procedural steps:
!   1) Read sources selected by FlashImage
!   2) Use a similar fit procedure as used for FCalibrate
!
!  ToDo:
!     - Center windows in time
!     - Apply intensity optimization using fine voxels
!     - Apply intensity optimization using chi-square optimization for cross correlation with beamed result
!     - Generalize read-in to wotk for TRI-D imager results (do not use peak ID)
!     -
   use constants, only : dp, pi
   use DataConstants, only : RunMode  ! needed for GLE-control
   use DataConstants, only : DataFolder, OutFileLabel, ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows
   use GLEplots, only : GLEplotControl
   use ReadPlane, only : ReadAirplane
   Implicit none
   !
   Real(dp), allocatable :: SourceTime_ms(:), SourcePos(:,:), SourceIntensity(:), Chi2(:)
   Complex, allocatable :: Stk_NEh(:,:,:)
   Integer, allocatable :: PeakWidth(:), SrcQual(:,:)
   !Character(len=80) :: OutFileLabel
   character*180 :: lname
   Real(dp), allocatable :: MovingSys(:,:), SourceUn(:), SourceLin(:), SourceCirc(:), SourceI3(:), SourceStI(:)
   real(dp), allocatable :: SourcePolZen(:,:), SourcePolAzi(:,:), SourcePolMag(:,:), SourcePoldOm(:,:)
   Integer, allocatable :: Cath(:)
   Complex :: Stk_ltu(1:3,1:3)
   Complex :: Stk_ltu_Cath(1:3,1:3,-3:+3)
   !
   Character(len=80):: Track, BackgrndData
   Real(dp), parameter :: FlightAngle=73.429
   Real(dp), parameter :: BankAngle=14.0
   Real(dp) :: PlanePos(1:3), PlaneVel(1:3)
   Real(dp) :: ltu2NEh(1:3,1:3)
   integer :: PeakNrTotal, i_Peak, j, k, nxx, DataUnit, DataUnitp, DataUnits, N_Cat(-3:3), N_backgr ! i,
   !
   real(dp) :: Cnst, StI, PolZen(1:3), PolAzi(1:3), PolMag(1:3), PoldOm(1:3)
   real(dp) :: x_plot,y_plot, z_plot, time, A, B, lon, tra, up, PolPlotScale
   logical :: prin=.true.
   character(len=2) :: txt
   Character(len=20), parameter :: Utility='LOFLi-Airplane', release='24.11 (Nov, 2024)'
   !
   Call GetNonZeroLine(lname)
   Read(lname,*,iostat=nxx) OutFileLabel, Track, BackgrndData
   OPEN(UNIT=2,STATUS='unknown',ACTION='WRITE',FILE='Airplane'//trim(OutFileLabel)//'.out')
   Call System_Initiation(Utility, release, ProgramFolder, UtilitiesFolder, FlashFolder, FlashName, Windows)  ! set Flash name & folder names
   RunMode=8
   !
   write(2,*) 'OutFileLabel:"',TRIM(OutFileLabel),'"'
   If(nxx.ne.0) Then
      write(2,*) 'Input line:',lname
      stop 'read problem'
   EndIf
   Call ReadAirplane(OutFileLabel, BackgrndData, PeakNrTotal, SourceTime_ms, SourcePos, SourceIntensity, PeakWidth, &
         Chi2, SrcQual, Stk_NEh, N_backgr)
   OutFileLabel=TRIM(OutFileLabel)//'Analyze'
   !
   Allocate( MovingSys(1:3,-7:PeakNrTotal), Cath(-7:PeakNrTotal) )
   Allocate( SourcePolZen(1:3,-7:PeakNrTotal), SourcePolAzi(1:3,-7:PeakNrTotal), SourcePolMag(1:3,-7:PeakNrTotal), &
            SourcePoldOm(1:3,-7:PeakNrTotal), SourceStI(-7:PeakNrTotal) )
   Allocate( SourceUn(-7:PeakNrTotal), SourceLin(-7:PeakNrTotal), SourceCirc(-7:PeakNrTotal), SourceI3(-7:PeakNrTotal) )
   !
   Call LinearCorrectSetup(track, FlightAngle, BankAngle, PlanePos, PlaneVel, ltu2NEh)
   !Call GLEplotControl(SpecialCmnd='echo on')  ! to generate unit=10 with the proper name
   write(2,*) 'plane position @ t=0. [m]:',PlanePos(1:3), PlaneVel(1:3)
   write(2,*) 'plane velocity vector [m]:',PlanePos(1:3), PlaneVel(1:3)
   Stk_ltu_Cath(1:3,1:3,-3:+3)=0.
   !
   N_Cat(-3:3)=0
   Do i_Peak=1,PeakNrTotal
      write(txt,"(I2.2)") i_Peak
      !write(*,*) 'Peak#', i_Peak
      !
      write(2,*) 'sourcepos:',SourcePos(1:3,i_Peak),', planepos:', PlanePos
      Call LinearCorrect(i_peak, PlanePos, PlaneVel, ltu2NEh, SourceTime_ms(i_Peak), SourcePos(1:3,i_Peak), &
         MovingSys(1:3,i_Peak), Cath(i_Peak))
      StI = Real(Stk_NEh(1,1,i_Peak) + Stk_NEh(2,2,i_Peak) + Stk_NEh(3,3,i_Peak))  ! =\Lambda_0  from  PHYSICAL REVIEW E 66, 016615 ~2002!
      SourceI3(i_Peak) =  Real(Stk_NEh(3,3,i_Peak)) / StI                     ! =>\Lambda_8 = StI12 *sq(3)/2 - StI3 *sq(3)
      Sourcelin(i_Peak) = (3./4.)*( (Real(Stk_NEh(1,1,i_Peak)-Stk_NEh(2,2,i_Peak)))**2+ 4.*(Real(Stk_NEh(1,2,i_Peak)))**2+ &
         4.*(Real(Stk_NEh(1,3,i_Peak)))**2 + 4.*(Real(Stk_NEh(2,3,i_Peak)))**2 + &
         ((Real(Stk_NEh(1,1,i_Peak)+ Stk_NEh(2,2,i_Peak)-2*Stk_NEh(3,3,i_Peak)))**2)/3.)/(StI**2)
      Sourcecirc(i_Peak)= (3.)*( (Imag(Stk_NEh(2,3,i_Peak)))**2+(Imag(Stk_NEh(1,3,i_Peak)))**2 + &
         (Imag(Stk_NEh(2,1,i_Peak)))**2)/(StI**2)
      SourceUn(i_Peak) = 1.- Sourcelin(i_Peak) -Sourcecirc(i_Peak)
      SourceStI(i_Peak) = StI
      !
      !write(2,*) SourceTime_ms(i_Peak), SourcePos(1:3,i_Peak), MovingSys(1:3,i_Peak)
      write(2,"(A,g13.1,3F6.1, f9.3)") 'Polarization:',SourceStI(i_Peak), 100.*SourceUn(i_Peak), 100.*Sourcelin(i_Peak) &
         , 100.*Sourcecirc(i_Peak), StI*PeakWidth(i_Peak)
      Call RRotCTensor(Stk_ltu, ltu2NEh, Stk_NEh(1,1,i_Peak))
      Call PolPCACath(Stk_ltu, SourcePolZen(1,i_Peak) , SourcePolAzi(1,i_Peak), SourcePolMag(1,i_Peak), &
         SourcePoldOm(1,i_Peak), prin)
      !
      If(Cath(i_Peak).le.3 .and. Cath(i_Peak).ge.-3) Then
         If(Cath(i_Peak).eq.3) Cath(i_Peak)=-3
         Stk_ltu_Cath(1:3,1:3,Cath(i_Peak))=Stk_ltu_Cath(1:3,1:3,Cath(i_Peak)) + Stk_ltu(1:3,1:3)*PeakWidth(i_Peak)
         N_Cat(Cath(i_Peak)) = N_Cat(Cath(i_Peak)) + PeakWidth(i_Peak)
      EndIf
      Write(2,*) 'Peak#',i_peak,  'finalized ================================================'
      Flush(unit=2)
      !Call GLEplotControl(SpecialCmnd='ls') !  ${FlashFolder}/
      !
   EndDo
   !
   !  Do a PCA position analysis
   Call LocationCatPCA(PeakNrTotal, MovingSys, Cath, ltu2NEh)
   !
   !  Determine centers of the plots
   x_plot=sum(MovingSys(1,1:PeakNrTotal))/PeakNrTotal !-40.
   If(abs(x_plot) .gt. 50.) x_plot=0.
   y_plot=sum(MovingSys(2,1:PeakNrTotal))/PeakNrTotal
   If(abs(y_plot) .gt. 30.) y_plot=0.
   z_plot=sum(MovingSys(3,1:PeakNrTotal))/PeakNrTotal
   If(abs(z_plot) .gt. 20.) z_plot=0.
   !  Determine size plots
   Cnst=0
   Do j=1,PeakNrTotal
      If(Cath(j) .eq. -3) cycle
      time=(MovingSys(1,j)-x_plot)**2 + (MovingSys(2,j)-y_plot)**2 + (MovingSys(3,j)-z_plot)**2
      If(time.gt.Cnst) Cnst=time
   endDo
   cnst=sqrt(cnst) !/sqrt(2.)
   cnst=nint(2*cnst/5.)*5.+5.  ! plot size
   !
   ! work on summarized pol observables per category
   SourceIntensity(-7:0)=0
   Chi2(-7:0)=0.
   PeakWidth(-7:0)=0
   SourceStI(-7:0)=0.
   SourceI3(-7:0)=0
   Call RRotCTensor(Stk_ltu_Cath(1,1,3), ltu2NEh, Stk_NEh(1,1,0) )  ! background results to cath=3; already per sample
   N_Cat(3)= N_backgr
   Do k=-3,3
      i_Peak=k-3
      SourceStI(i_Peak)=1.
      StI = Real(Stk_ltu_Cath(1,1,k) + Stk_ltu_Cath(2,2,k) + Stk_ltu_Cath(3,3,k))  ! =\Lambda_0  from  PHYSICAL REVIEW E 66, 016615 ~2002!
      If(StI.lt.0.1) cycle
      SourceStI(i_Peak)=StI
      SourceLin(i_Peak)=3.*((REAL(Stk_ltu_Cath(1,2,k)))**2 +(REAL(Stk_ltu_Cath(1,3,k)))**2 + &
                              (REAL(Stk_ltu_Cath(2,3,k)))**2 )/(StI**2)
      SourceLin(i_Peak)=SourceLin(i_Peak)+ (3.*(Real(Stk_ltu_Cath(1,1,k)-Stk_ltu_Cath(2,2,k)) )**2 + &
                  (Real(Stk_ltu_Cath(1,1,k) + Stk_ltu_Cath(2,2,k)) -2*Real(Stk_ltu_Cath(3,3,k)))**2)/(4.*StI**2)
      SourceCirc(i_Peak)=3.*((IMAG(Stk_ltu_Cath(1,2,k)))**2 +(IMAG(Stk_ltu_Cath(1,3,k)))**2 + &
                              (IMAG(Stk_ltu_Cath(2,3,k)))**2 )/(StI**2)
      SourceUn(i_Peak)=1.-SourceLin(i_Peak)-SourceCirc(i_Peak)
      Cath(i_Peak)=k
      Chi2(i_Peak)=1.
      SourceTime_ms(i_Peak)=0
      MovingSys(1:3,i_Peak)=k*10.
      Call PolPCACath(Stk_ltu_Cath(1,1,k), SourcePolZen(1,i_Peak) , SourcePolAzi(1,i_Peak), SourcePolMag(1,i_Peak), &
         SourcePoldOm(1,i_Peak), prin)
      write(2,"(I4, F12.1,  3(' , ',f6.3, 3(',',f7.2)))") &
         k, StI,(SourcePolMag(j,i_Peak)/STI, SourcePolZen(j,i_Peak), SourcePolAzi(j,i_Peak), SourcePoldOm(j,i_Peak), j=1,3)
      write(2,"(2I4, 3(A,f7.2), A, F9.2)") k,N_Cat(k),', P_unpol=',SourceUn(i_Peak)*100., '%, P_lin=', SourceLin(i_Peak)*100., &
         '%, P_circ=', SourceCirc(i_Peak)*100., '%, I/sample=', StI/N_Cat(k)
   EndDo
   !
   DataUnitp=28
   DataUnits=29
   OPEN(UNIT=DataUnitp,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'CorrectedSources.dat')
   write(DataUnitp,"(I4,4F8.0, A)") PeakNrTotal, x_plot, y_plot, z_plot, Cnst, '   !'
   write(DataUnitp,"(2A)") '!  #   t_src[ms]  Length[m]  side[m]    up[m]       StI    chi^2 Categ  PrimeLinPol', &
      '    Long    Trans   Upward ;  %Unpol %linpol %circpol '
   OPEN(UNIT=DataUnits,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'CorrectedSummary.dat')
   write(DataUnits,"(I4,4F8.0, A)") PeakNrTotal, x_plot, y_plot, z_plot, Cnst, '   !'
   write(DataUnits,"(2A)") '!  #   t_src[ms]  Length[m]  side[m]    up[m]       StI    chi^2 Categ  PrimeLinPol', &
      '    Long    Trans   Upward ;  %Unpol %linpol %circpol '
   write(2,"(I4,4F8.0, A)") PeakNrTotal, x_plot, y_plot, z_plot, Cnst, '   !'
   write(2,"(I4,2A)") PeakNrTotal,'   t_src[ms]  Length[m]  side[m]   up[m]  Intens   chi^2 Cath Width  unpol[%]', &
      '  lin[%]  circ[%]    I3[%]     Primary{I[%], th_1, ph_1, d_omg};   Secondary{ }                ; Tertiary{ }'
   PolPlotScale=5.
   DataUnit=DataUnits
   Do i_Peak=-5,PeakNrTotal
      If(i_peak.gt.0) DataUnit=DataUnitp
      If(SourceStI(i_Peak).eq.0) cycle
      !Cnst=100./sqrt(sum(SourcePolMag(1:3,i_Peak)*SourcePolMag(1:3,i_Peak)))
      Cnst=100./(sum( SourcePolMag(1:3,i_Peak) ))
      A= SourcePolMag(1,i_Peak)* cos(2.*SourcePoldOm(1,i_Peak)*pi/180.)
      lon=PolPlotScale*sin(SourcePolZen(1,i_Peak)*pi/180.) * cos(SourcePolAzi(1,i_Peak)*pi/180.)
      tra=PolPlotScale*sin(SourcePolZen(1,i_Peak)*pi/180.) * sin(SourcePolAzi(1,i_Peak)*pi/180.)
      up=PolPlotScale*cos(SourcePolZen(1,i_Peak)*pi/180.)
      write(2,"(I4,F14.6, 3F9.2, F8.1, F8.2, 2I5, 4(' , ',f6.1), 3(' , ',f8.1, 3(',',f7.2)))") &
         i_peak, SourceTime_ms(i_Peak), MovingSys(1:3,i_Peak), SourceIntensity(i_Peak), Chi2(i_Peak), Cath(i_Peak) &
         , PeakWidth(i_Peak), SourceUn(i_Peak)*100., SourceLin(i_Peak)*100., SourceCirc(i_Peak)*100., SourceI3(i_Peak)*100.  &
         , (SourcePolMag(j,i_Peak)*Cnst, SourcePolZen(j,i_Peak), SourcePolAzi(j,i_Peak), SourcePoldOm(j,i_Peak), j=1,3)
      If(Cath(i_Peak).le.-3) cycle
      write(DataUnit,"(I4,F14.6, 3F9.2, F12.1, F8.2, I5, F12.1, 3(' , ',f6.2), 3(' , ',f6.1))") &
         i_peak, SourceTime_ms(i_Peak), MovingSys(1:3,i_Peak), SourceStI(i_Peak), Chi2(i_Peak), Cath(i_Peak) &
         , A, lon, tra, up, SourceUn(i_Peak)*100., SourceLin(i_Peak)*100., SourceCirc(i_Peak)*100.
   EndDo
   Close(unit=DataUnit)
   Call GLEplotControl(PlotType='AirPlanePolariz', PlotName='AirplanePol'//TRIM(OutFileLabel), &
      PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel) )
   Call GLEplotControl(PlotType='AirplaneEvents', PlotName='Airplane'//TRIM(OutFileLabel), &
      PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel) )
   Call GLEplotControl(PlotType='AirplaneEventsPol', PlotName='AirplanePol'//TRIM(OutFileLabel), &
      PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel) )
   !Call GLEplotControl(CleanPlotFile=TRIM(OutFileLabel)//'Interferometer*.z')
   !Call GLEplotControl(CleanPlotFile=TRIM(OutFileLabel)//'IntfSpecPowMx_d.dat')
   !
   SourceUn(-3:3)=0.0001 ! just to make it non-zero
   SourceLin(-3:3)=0.
   SourceCirc(-3:3)=0.
   Do i_Peak=1,PeakNrTotal
      If(SourceStI(i_Peak).eq.0) cycle
      k=Cath(i_Peak)
      SourceUn(k)=SourceUn(k)+1
      SourceLin(k)=SourceLin(k) + SourceLin(i_Peak)
      SourceCirc(k)=SourceCirc(k) + SourceCirc(i_Peak)
   EndDo
   Do k=-3,3
      write(2,"(A,I2,A,F6.1,A,F6.1, A, I5)") 'Average polarization % for Cath=',k, ', lin=',100.*SourceLin(k)/SourceUn(k) &
         , ', circ=',100.*Sourcecirc(k)/SourceUn(k), ', #=',NINT(SourceUn(k))
   EndDo
   !
   Call GLEplotControl( Submit=.true.)
   !
   Return
   !
   !====================================================================
   !
End Program AirplaneSources
!-----------------------------------------
Subroutine LocationCatPCA(PeakNrTotal, MovingSys, Cath, ltu2NEh)
   use Constants, only : dp, pi !,sample,c_mps, pi, sample_ms
   !use DataConstants, only : DataFolder, OutFileLabel
   Implicit none
   Integer, intent(in) :: PeakNrTotal
   Real(dp), intent(in) :: MovingSys(1:3,-7:PeakNrTotal)
   Integer, intent(in) :: Cath(-7:PeakNrTotal)
   Real(dp), intent(in) :: ltu2NEh(1:3,1:3)
   Integer :: i_source, i,j,k
   Real(dp) :: LocatCovariance(1:3,1:3,-3:3), LocateMean(1:3,-3:3), LocCor(1:3)
   Real(dp) :: B(3,3), eigval(1:3), eigvec(1:3,1:3), LocNEh(1:3)
   Integer :: NCat(-3:3)
   !
   write(2,*) 'perform PCA-analysis on location distribution'
   LocateMean(:,:)=0.
   NCat(:)=0
   Do i_source=1,PeakNrTotal
      k=Cath(i_source)
      If(k.lt.-3 .or. k.gt.3) cycle
      NCat(k)=NCat(k)+1
      LocateMean(:,k)=LocateMean(:,k) + MovingSys(:,i_source)
   EndDo
   Do k=-3,3
      If(NCat(k).le.0) cycle
      LocateMean(:,k)=LocateMean(:,k)/NCat(k)
      write(2,*) 'Average position of sources of Category=',k,LocateMean(:,k),NCat(k)
   Enddo
   !
   LocatCovariance(:,:,:)=0.
   Do i_source=1,PeakNrTotal
      k=Cath(i_source)
      If(k.lt.-3 .or. k.gt.3) cycle
      If(NCat(k).le.3) cycle
      LocCor(:)=MovingSys(:,i_source)-LocateMean(:,k)
      Do i=1,3
         LocatCovariance(:,i,k)=LocatCovariance(:,i,k) + LocCor(:)*LocCor(i)/NCat(k)
      EndDo
   EndDo
   !
   Do k=-3,3
      write(2,*) 'Categ=',k, LocatCovariance(1,1,k), NCat(k)
      If(NCat(k).le.5) cycle
      Call Inverse33(LocatCovariance(1,1,k), B, eigval, eigvec) ! vec(1:3,i) with val(i)
      Do i=1,3 ! loop over all eigenvalues
         Do j=1,3 ! rotate to Earth-oriented frame NEh
            LocNEh(j)=sum(ltu2NEh(1:3,j)*eigvec(1:3,i))
         EndDo
         write(2,"(A, I2,A,I1, A,F8.2, A,3F6.3, A,3F6.3)") 'category=',k,', axis:',i, &
            ', sqrt(Cov)=',sqrt(eigval(i)),'[m], orientation_Airplane:',eigvec(1:3,i),', orientation_NEh:',LocNEh(1:3)
      EndDo
   EndDo
   !
   Return
End Subroutine LocationCatPCA
!-----------------------------------------
Subroutine LinearCorrectSetup(track, FlightAngle, BankAngle, d_0, d_t, ltu2NEh)
   use Constants, only : dp, pi !,sample,c_mps, pi, sample_ms
   use DataConstants, only : DataFolder, OutFileLabel
   Implicit none
   Character(len=80), intent(in) :: Track  ! should not be longer
   real(dp), intent(in) :: FlightAngle, BankAngle
   Real(dp), intent(out) :: ltu2NEh(1:3,1:3)  ! rotation matrix to body-aligned coordinates
   Real(dp), intent(out) :: d_0(1:3), d_t(1:3)
   Real(dp) :: S_t, S_tt, S_d(1:3), S_dt(1:3), t_trc, d_trc(1:3), RotMat(1:3,1:3)
   Real(dp) :: VelPlane, Bank_rad
   Integer :: N_trc, nxx, i
   !
   d_0(1:3)=(/  23.22366d0,  -41.86376d0,  8.05660d0 /)
   d_0(1:3)=d_0(1:3)*1000.
   d_t(1:3)=0.
   ltu2NEh(1:3,1:3)=0.
   Do i=1,3
      ltu2NEh(i,i)=1.
   EndDo
   !
   If(track.eq.'') Return
   !
   OPEN(UNIT=29,STATUS='unknown',ACTION='read',FILE=trim(DataFolder)//TRIM(Track))
   N_trc=0
   nxx=0
   S_t =0.
   S_tt=0.
   S_d(1:3) =0.
   S_dt(1:3)=0.
   write(2,*) 'reading airplane track from:',trim(DataFolder)//TRIM(Track)
   Do    ! fit d_i=d_0 + d_t*t_i
      read(29,*,iostat=nxx) i, t_trc, d_trc(1:3)
      If(nxx.ne.0) exit
      N_trc=N_trc+1
      !t_test(N_trc)= t_trc
      !write(2,*)  t_trc, d_trc(1:3)
      Call Convert2m( d_trc(1:3) )
      S_t =S_t  + t_trc
      S_tt=S_tt + t_trc*t_trc
      S_d(1:3) =S_d(1:3)  + d_trc(1:3)
      S_dt(1:3)=S_dt(1:3) + d_trc(1:3)*t_trc
   EndDo
   Close(Unit=29)
   !
   d_0(1:3) = (S_dt(1:3)*S_t  - S_d(1:3) * S_tt)/(S_t*S_t - N_trc*S_tt)
   d_t(1:3) = (S_d(1:3) * S_t - S_dt(1:3)*N_trc)/(S_t*S_t - N_trc*S_tt)
   write(2,*) 'track:', d_0, d_t
   write(2,*) 'speed:', 1000.*sqrt(sum(d_t(1:3)*d_t(1:3))), ' [m/s]', atan2(d_t(2),d_t(1))*180./pi
   d_t(1)=d_t(2)/tan(FlightAngle*pi/180.)
   d_t(1:3) = d_t(1:3)*1.015  ! re-adjust the velocity
   VelPlane=sqrt(sum(d_t(1:3)*d_t(1:3)))
   write(2,*) 'Re-adjust angle:', 1000.*VelPlane, ' [m/s]', atan2(d_t(2),d_t(1))*180./pi
   !
   ! Construct rotation matrix RotMat
   ! vec_(plane-aligned)(i) = RotMat(i,j) vec_{BdyCentr}(j) ! ltu2NEh first is ltu second NEh
   RotMat(1,1:3)=d_t(1:3)/VelPlane
   RotMat(2,1)=d_t(2)
   RotMat(2,2)=-d_t(1)
   RotMat(2,3)=0
   RotMat(2,1:3)=RotMat(2,1:3)/sqrt(sum(RotMat(2,1:3)*RotMat(2,1:3)))
   RotMat(3,1)= -d_t(1)*d_t(3)                  ! -d_trc(2)*d_t(3) + d_trc(3)*d_t(2) !
   RotMat(3,2)= -d_t(2)*d_t(3)                 ! -d_trc(3)*d_t(1) + d_trc(1)*d_t(3)
   RotMat(3,3)= +d_t(2)*d_t(2) + d_t(1)*d_t(1)  ! -d_trc(1)*d_t(2) + d_trc(2)*d_t(1)
   RotMat(3,1:3)= RotMat(3,1:3)/sqrt(sum(RotMat(3,1:3)*RotMat(3,1:3)))
   !
   ! correct for banking angle
   Bank_rad=BankAngle*pi/180.
   ! l_c= l  ; t_c= t cos(b) - u * cin(b)  ; u_c= t sin(b) + u * cos(b)
   !(l,t,u)_c= sum(ltu2NEh(ltu_c,1:3)*shft(1:3))
   ltu2NEh(1,1:3) = RotMat(1,1:3) ! gives new longitudinal
   ltu2NEh(2,1:3) = cos(Bank_rad)*RotMat(2,1:3)-sin(Bank_rad)*RotMat(3,1:3) ! gives new transverse
   ltu2NEh(3,1:3) = cos(Bank_rad)*RotMat(3,1:3)+sin(Bank_rad)*RotMat(2,1:3) ! gives new up
   !S_t=-1000.
   !write(2,*) 'Airplane @:', S_t, ' [ms]', d_0(1:3)+ d_t(1:3)*S_t
   !S_t=-00.
   !write(2,*) 'Airplane @:', S_t, ' [ms]', d_0(1:3)+ d_t(1:3)*S_t
   !S_t=+1000.
   !write(2,*) 'Airplane @:', S_t, ' [ms]', d_0(1:3)+ d_t(1:3)*S_t
   !Do i=1,N_trc
   !   S_t= t_test(i)
   !   write(2,"(A, F13.6, A, 3F10.4)") 'Airplane @:', S_t, ' [ms]', (d_0(1:3)+ d_t(1:3)*S_t)/1000.
   !
   Return
   !
End Subroutine LinearCorrectSetup
!-----------------------------------------
Subroutine LinearCorrect(i_peak, d_0, d_t, ltu2NEh, SourceTime, SourcePos, Rot, Cath)  !SourceIntensity, Chi2, &
   !      PeakWidth, SourcePolMag, SourcePolZen, SourcePolAzi, SourcePoldOm, SourceUn, SourceLin, SourceCirc, I3, DataUnit )
   use Constants, only : dp, pi !,sample,c_mps, pi, sample_ms
   !use DataConstants, only : DataFolder, OutFileLabel
   Implicit none
   Integer, intent(in) :: i_Peak
   real(dp), intent(in) :: d_0(1:3), d_t(1:3), SourceTime, SourcePos(1:3)  !, SourceIntensity, Chi2
   Real(dp), intent(in) :: ltu2NEh(1:3,1:3)  ! rotation matrix to body-aligned coordinates
   real(dp), intent(out) :: Rot(1:3)  ! rotated coordinates
   Integer, intent(out) :: Cath
   Real(dp) :: A, B, F
   Real(dp), parameter :: SourceScatt(1:3)=(/ 0.489 , -0.323 , -0.810 /) ! from PCA analysis left engine
   Real(dp):: shft(1:3)
   Integer :: N_trc, nxx, i
   !
   shft(1:3) = (SourcePos(1:3)- d_0(1:3)- d_t(1:3)*SourceTime)  ! still in NEh system, corrected for velocity
   Do i=1,3 ! ! ltu2NEh first is ltu second NEh !NEH2AP First is NEh orientation, second is airplane-oriented frame,
      Rot(i)=sum(ltu2NEh(i,1:3)*shft(1:3))
   EndDo
   !
   !  Determine source
   !
   A=rot(1) + rot(2)  ! perpendicular to stripes
   B=rot(1) - rot(2) - 2*rot(3)  ! along stripes
   Cath=0
   If(rot(2) .lt. -20) then
      Cath=-3
   ElseIf(rot(2) .gt. 20) then
      Cath=-3
   ElseIf(A.gt. 55) Then
      Cath=2
   ElseIf(A.gt. 50) Then
      Cath=1
   ElseIf(A.gt. 34) Then
      Cath=-2
   ElseIf(A.gt. 20) Then
      Cath=-1
   EndIf
   !
   Write(2,"(A, I3, 3F7.2, A, 3F7.2, A, I3,2F7.2)") 'rotated along airplane:', i_peak, rot(1:3), & !sqrt(sum(rot(1:3)*rot(1:3)))
         ' [m], (N,E,h)[m]=', shft(1:3),', Category:',Cath, A, B
   !
   F=rot(3)/SourceScatt(3)  ! distance from plane along scatter
   A=rot(1)-SourceScatt(1)*F ! Along distance when projecting on up=0 plane
   B=rot(2)-SourceScatt(2)*F ! Transverse distance when projecting on up=0 plane
   If(ABS(F) .gt. 20) Then
      If(Cath.ne.-3) Write(2,*) '!!!!!!!!!!!!!!!!!!Recategorize', Cath, F, A, B, ' to', -3
      Cath=-3
   ElseIf(A .lt. 20) Then ! tail
      If(ABS(B) .gt. 10) then
         If(Cath.ne.-3) Write(2,*) '!!!!!!!!!!!!!!!!!!Recategorize', Cath, F, A, B, ' to', -3
         Cath=-3
      Else
         If(Cath.ne.0) Write(2,*) '!!!!!!!!!!!!!!!!!!Recategorize', Cath, F, A, B, ' to', 0
         Cath=0
      EndIf
   ElseIf(B .lt. 0) Then  ! right side
      If(A .lt. 41) Then   ! right wing
         If(Cath.ne.-1) Write(2,*) '!!!!!!!!!!!!!!!!!!Recategorize', Cath, F, A, B, ' to', -1
         Cath=-1
      Else                 ! right engine
         If(Cath.ne.-2) Write(2,*) '!!!!!!!!!!!!!!!!!!Recategorize', Cath, F, A, B, ' to',-2
         Cath=-2
      EndIf
   Else
      If(A .lt. 41) Then   ! left wing
         If(Cath.ne.+1) Write(2,*) '!!!!!!!!!!!!!!!!!!Recategorize', Cath, F, A, B, ' to', 1
         Cath=1
      Else                 ! left engine
         If(Cath.ne.+2) Write(2,*) '!!!!!!!!!!!!!!!!!!Recategorize', Cath, F, A, B, ' to', 2
         Cath=2
      EndIf
   EndIf
   Return
   !
End Subroutine LinearCorrect
!=================================
!-----------------------------------------
Subroutine PolarCorrect(i_peak, FlightAngle, SourceTime, SourcePos, Rot, Cath)  !SourceIntensity, Chi2, &
   !      PeakWidth, SourcePolMag, SourcePolZen, SourcePolAzi, SourcePoldOm, SourceUn, SourceLin, SourceCirc, I3, DataUnit )
   use Constants, only : dp, pi !,sample,c_mps, pi, sample_ms
   use DataConstants, only : DataFolder, OutFileLabel
   !use Interferom_Pars, only :  Nr_IntFerMx, Nr_IntferCh ! the latter gives # per chunk
   !use Interferom_Pars, only : Cnu_p0, Cnu_t0, Cnu_p1, Cnu_t1, IntfNuDim, AntPeak_OffSt
   !Use Interferom_Pars, only : N_smth, N_fit, i_chunk
   !Use Interferom_Pars, only : dStI, dStI12, dStQ, dStU, dStV, dStI3, dStU1, dStV1, dStU2, dStV2, Chi2pDF
   !Use Interferom_Pars, only : Alloc_EInterfImag_Pars, DeAlloc_EInterfImag_Pars
   Implicit none
   Integer, intent(in) :: i_Peak
   real(dp), intent(out) :: Rot(1:3)  ! rotated coordinates
   real(dp), intent(in) :: SourcePos(1:3), SourceTime, FlightAngle
   Integer, intent(out) :: Cath
   !real(dp), intent(in) :: SourcePolMag(1:3), SourcePolZen(1:3), SourcePolAzi(1:3), SourcePoldOm(1:3)
   !real(dp), intent(in) :: SourceUn, SourceLin, SourceCirc, I3
   !Integer, intent(in) :: PeakWidth, DataUnit
   !logical, intent(in) :: List
   Real(dp) :: S_t, S_tt, S_d(1:3), S_dt(1:3), t_trc, d_trc(1:3)
   Real(dp), save :: d_0(1:3), d_t(1:3), shft(1:3)  , t_test(1:100)
   Integer :: N_trc, nxx, i
   !
   rot(1)= sum(shft(1:3)*d_t(1:3))/sqrt(sum(d_t(1:3)*d_t(1:3)))
   d_trc(1)=  +d_t(2)
   d_trc(2)=  -d_t(1)
   d_trc(3)=0
   rot(2)= sum(shft(1:3)*d_trc(1:3))/sqrt(sum(d_trc(1:3)*d_trc(1:3)))
   d_trc(1)= -d_t(1)*d_t(3)                  ! -d_trc(2)*d_t(3) + d_trc(3)*d_t(2) !
   d_trc(2)= -d_t(2)*d_t(3)                 ! -d_trc(3)*d_t(1) + d_trc(1)*d_t(3)
   d_trc(3)= +d_t(2)*d_t(2) + d_t(1)*d_t(1)  ! -d_trc(1)*d_t(2) + d_trc(2)*d_t(1)
   rot(3)= sum(shft(1:3)*d_trc(1:3))/sqrt(sum(d_trc(1:3)*d_trc(1:3)))
   !If(DataUnit.gt.2) Then
   !   S_t=100./sqrt(sum(SourcePolMag(1:3)*SourcePolMag(1:3)))
   !   write(DataUnit,"(I4,F14.6, 3F9.2, F8.1, F8.2, I5, 4(' , ',f6.1), 3(' , ',f8.1, 3(',',f7.2)))") &
   !      i_peak, SourceTime, rot(1:3), SourceIntensity &
   !      , Chi2, PeakWidth, SourceUn*100., SourceLin*100., SourceCirc*100., I3*100.  &
   !      , (SourcePolMag(i)*S_t, SourcePolZen(i), SourcePolAzi(i), SourcePoldOm(i), i=1,3)
   !   write(2,"(I4,F14.6, 3F9.2, F8.1, F8.2, I5, 4(' , ',f6.1), 3(' , ',f8.1, 3(',',f7.2)))") &
   !      i_peak, SourceTime, rot(1:3), SourceIntensity &
   !      , Chi2, PeakWidth, SourceUn*100., SourceLin*100., SourceCirc*100., I3*100.  &
   !      , (SourcePolMag(i)*S_t, SourcePolZen(i), SourcePolAzi(i), SourcePoldOm(i), i=1,3)
   !Else
   !
   !  Determine source
   !
   S_t=rot(1) + rot(2)  ! perpendicular to stripes
   S_tt=rot(1) - rot(2) - 2*rot(3)  ! along stripes
   Cath=0
   If(rot(2) .lt. -20) then
      Cath=-3
   ElseIf(rot(2) .gt. 20) then
      Cath=3
   ElseIf(S_t.gt. 50) Then
      Cath=2
   ElseIf(S_t.gt. 40) Then
      Cath=1
   ElseIf(S_t.gt. 35) Then
      Cath=-2
   ElseIf(S_t.gt. 410) Then
      Cath=-1
   EndIf
   !
      Write(2,"(A, I3, 3F7.2, A, 3F7.2, A, I3,2F7.2)") 'rotated along airplane:', i_peak, rot(1:3), & !sqrt(sum(rot(1:3)*rot(1:3)))
         ' [m], (N,E,h)[m]=', shft(1:3),', Category:',Cath, S_t, S_tt
   !EndIf
   !
   Return
   !
End Subroutine PolarCorrect
!==========================================

!================================
!-----------------------------------------------
!=====================================================
!-------------------------------------------------------------
!=========================================
!==========================
