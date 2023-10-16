Subroutine SI_Run
    !+++++++++ What is used:
    !- Loop over antennas
    ! - make frequency spectra
    !- End antennaLoop
    !- Loop over all pixels
    !  - calculate phase shift Each antenna
    !   - Loop frequency
    !    - Sum(antenna) spectra
    ! - End pixelloop
    !+++++++++ Alternative:
    !- Loop over antennas
    ! - make frequency spectra
    ! - Loop over all pixels
    !    - calculate phase shift
    !    - Sum frequency: update spectra with this antenna
    ! - End pixelLoop
    !- End antennaLoop
    !++++++++++++++++++++++
    !CTime_sum, SumStrt, SumWindw, IntFer_ant, IntfLead,
   Use Interferom_Pars, only : polar, N_pix, d_loc, CenLocPol, CenLoc, PixLoc, Nr_IntFer, IntFer_ant
   Use Interferom_Pars, only : AveInten, AveIntenE, AveIntenN, RimInten, MaxIntfInten !, MaxIntfIntenLoc
   Use Interferom_Pars, only : MaxSlcInten, MaxSmPow, IntPowSpec
   use constants, only : dp, ci, pi !, Sample, Refrac, c_mps
   use StationMnemonics, only : Statn_ID2Mnem
   use Chunk_AntInfo, only : SignFlp_SAI, PolFlp_SAI, BadAnt_SAI
   use GLEplots, only : GLEplotControl
   Use Interferom_Pars, only : Alloc_SInterferom_Pars,DeAlloc_SInterferom_Pars
   use mod_test
   Implicit none
   integer :: i_eo, i !,j, i_s, i_loc(1), i_ant, j_corr, Station, j_IntFer, i_nu, nxx
   Real(dp) :: PixLocPol(1:3) !, SmPow
   Integer :: Date_T(8)
   integer :: i_N, i_E, i_h !, N_h,  N_hlow, N_Eup, N_Elow
   !------------------------------------
   ! main loop over direction antennas
   Do i_eo=0, 1
      CALL DATE_AND_TIME (Values=DATE_T)
      WRITE(2,"(1X,I2,':',I2.2,':',I2.2,'.',I3.3,A,i1)") (DATE_T(i),i=5,8), ' start Signal Interferometry for antennas=',i_eo
      WRITE(*,"(1X,I2,':',I2.2,':',I2.2,'.',I3.3,A)") (DATE_T(i),i=5,8), achar(27)//'[0m'
      !
      Call SISelectAntennas(i_eo)
      Call PixBoundingBox()      !  allocate array space
      Call Alloc_SInterferom_Pars
      Call SISetupSpec(i_eo)
      !
      ! -----------------------------------------------
      ! main loop over Voxel location
      ! do interferometry
      MaxIntfInten=0.
      MaxSlcInten(:)=0.
      AveInten(:,:)=0.
      AveIntenE(:,:)=0.
      AveIntenN(:,:)=0.
      RimInten(:,:)=0.
      MaxSmPow(:)=0.
      Do i_h= N_pix(3,1), N_pix(3,2) ! or distance
         CALL DATE_AND_TIME (Values=DATE_T)
         WRITE(2,"(1X,I2,':',I2.2,':',I2.2,'.',I3.3,A,i3)") (DATE_T(i),i=5,8), ' height',i_h
         WRITE(*,"(A,i3,1X,I2,':',I2.2,':',I2.2,'.',I3.3,A)") ' height',i_h,(DATE_T(i),i=5,8), achar(27)//'[0m'
         !
         Do i_N=N_pix(1,1),N_pix(1,2)   ! or Phi  ! Loop over pixels
         Do i_E=N_pix(2,1),N_pix(2,2)   ! or elevation angle
            If(polar) then
               PixLocPol(1)=CenLocPol(1)+i_N*d_loc(1) ! azimuth
               PixLocPol(2)=CenLocPol(2)+i_E*d_loc(2) ! elevation angle
               PixLocPol(3)=CenLocPol(3)+i_h*d_loc(3)
               If(PixLocPol(2).gt.pi/2 .or. PixLocPol(2).lt.0.) stop 'Interferometry; theta problem'
               Call Pol2Carth(PixLocPol,PixLoc)
            else
               PixLoc(1)=CenLoc(1)+i_N*d_loc(1)
               PixLoc(2)=CenLoc(2)+i_E*d_loc(2)
               PixLoc(3)=CenLoc(3)+i_h*d_loc(3)
            Endif
            If(PixLoc(3).lt.0.) stop 'Interferometry; height problem'
            !
            ! sum frequency spectra
            !If(i_N.eq.0 .and. i_E.eq.0) then
               !CALL DATE_AND_TIME (Values=DATE_T)
               !WRITE(2,"(1X,I2,':',I2.2,':',I2.2,'.',I3.3,A,i2)") (DATE_T(i),i=5,8), ' A'
               !WRITE(*,"(A,1X,I2,':',I2.2,':',I2.2,'.',I3.3,A)") ' A:',(DATE_T(i),i=5,8), achar(27)//'[0m'
            !endif
            !
            Call InterfEngineB(i_eo, i_N, i_E, i_h, Nr_IntFer)
            !InterfEngine
         EndDo ! i_E
         EndDo ! i_N
      EndDo   ! Loop over distances or heights
      !
      !----------------------------------------
      Call OutputIntfPowrTotal(IntFer_ant(1,1), i_eo)
      If(IntPowSpec) Call OutputIntfPowrMxPos(i_eo)
      Call OutputIntfSlices(i_eo)
      !
      CALL DATE_AND_TIME (Values=DATE_T)
      WRITE(2,"(1X,I2,':',I2.2,':',I2.2,'.',I3.3,A,i1)") (DATE_T(i),i=5,8), ' end Interferometry',i_eo
      Call DeAlloc_SInterferom_Pars
   Enddo !  i_eo=0,1
   Return
End Subroutine SI_Run
!------------------------------------------
! =========================
Subroutine InterfEngineB(i_eo, i_N, i_E, i_h, Nr_IntFer)
!
   use constants, only : dp, ci, pi, Sample, Refrac, c_mps
   use DataConstants, only : Time_Dim, Cnu_dim, DataFolder, OutFileLabel, Calibrations
   use Chunk_AntInfo, only : Ant_pos, Ant_RawSourceDist
   Use Interferom_Pars, only : i_chunk, IntFer_ant, IntfDim, IntfNuDim, IntfLead, PixLoc  ! ,Nr_IntFer
   Use Interferom_Pars, only : Cnu, Cnu_pix, CTime_pix
   use FFT, only : RFTransform_CF2CT
   Implicit none
   Integer, intent(in) :: i_eo, i_N, i_E, i_h, Nr_IntFer
   Integer :: i_ant, i_nu, j_IntFer
   Complex(dp) :: Phase(Nr_IntFer), dPhase(Nr_IntFer)
   complex(dp), parameter :: tipi=2*ci*pi
   Real(dp) :: RDist, dt_AntPix
   !InterfEngine
   Cnu_pix(:)=0.
   Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
      i_ant=IntFer_ant(j_IntFer,1)
      Call RelDist(PixLoc(1),Ant_pos(1,i_ant,i_chunk),RDist)
      dt_AntPix=Rdist - Ant_RawSourceDist(i_ant,i_chunk)
      ! check dt_AntPix in range
      If(ABS(dt_AntPix).gt.IntfLead) then
         write(2,*) 'Warning, time-shift=',dt_AntPix,'[samples] out of range for pixel',i_N, i_E, i_h
         Return
      Endif
      dphase(j_IntFer)= exp(-tipi*dt_AntPix/IntfDim)
   EndDo ! j_IntFer
   Phase(:)=1.d0
   ! Cnu produced in Subroutine SetupIntfSpec, should be modified for version B
   Do i_nu=0,IntfNuDim   ! Increment frequency spectrum with this antenna
      Cnu_pix(i_nu)=Cnu_pix(i_nu) + SUM(phase(:)*Cnu(:,i_nu))  ! Note: order indices changed in Cnu  needs v18d (for otions and module)
      phase(:)=phase(:)*dphase(:)
   EndDo
   !
   Call RFTransform_CF2CT(Cnu_pix,CTime_pix )
   !
   Call SIAnalyzePixelTTrace(i_eo, i_N, i_E, i_h)
   !
   Return
End Subroutine InterfEngineB
!-----------------------------------------------
Subroutine SISelectAntennas(i_eo)
   ! Select antennas within a range of 'AntennaRange' from the core
   !--------------------------------------------
   use constants, only : dp !, pi, Sample, Refrac, c_mps
   use Chunk_AntInfo, only : Ant_Stations, Ant_IDs, Ant_nr, Ant_pos !, Ant_RawSourceDist
   use Chunk_AntInfo, only : Unique_StatID,  Nr_UniqueStat, RefAnt  ! for curtainplot
   Use Interferom_Pars, only : Nmax_IntFer
   Use Interferom_Pars, only : i_chunk, IntFer_ant, Nr_IntFer
   use FitParams, only : AntennaRange
   use StationMnemonics, only : Statn_ID2Mnem
   Implicit none
   Integer, intent(in) :: i_eo
   real ( kind = 8 ) :: Dist
   integer :: i_ant, Station, j_IntFer ! , MLoc, Mloc_all
   !
   Station=0
   j_IntFer=0
   Nr_UniqueStat=0
   Do i_ant=1,Ant_nr(i_chunk)         ! Select the antennas for Interferometry
      If(Ant_Stations(i_ant,i_chunk) .ne. Station) then
        Dist=sqrt(sum((Ant_pos(:,i_ant,i_chunk)-Ant_pos(:,1,i_chunk))**2))
        Station=Ant_Stations(i_ant,i_chunk)
        Nr_UniqueStat=Nr_UniqueStat+1
        Unique_StatID(Nr_UniqueStat)=Station  ! for Curtainplot
        !write(2,*) Statn_ID2Mnem(Station), Ant_pos(:,i_ant,i_chunk)/1000.
      endif
      !write(2,*) 'i_ant, Station, Dist',i_ant, Station, Dist
      If(Dist .gt. AntennaRange*1000.) cycle  ! Limit to antennas near the superterp
      if(mod(Ant_IDs(i_ant,i_chunk),2) .ne. i_eo) cycle       ! limit to odd antennas
      j_IntFer=j_IntFer+1
      IntFer_ant(j_IntFer,i_chunk)=i_ant
      If(j_IntFer .ge. Nmax_IntFer) then
         Write(2,*) 'max nr of interferometry-antennas reached, station=',Station, Ant_IDs(i_ant,i_chunk),', at distance ', Dist
         exit
      Endif
   EndDo
   !stop
   Nr_Intfer=j_IntFer
!   SlcNrm=1./(SliceLen*Nr_IntFer*Nr_IntFer)
   RefAnt(i_chunk,i_eo)=IntFer_ant(1,i_chunk)  ! for curtainplot
   Return
End Subroutine SISelectAntennas
!-----------------------------------------------
Subroutine SISetupSpec(i_eo)
   ! Get Summed spectrum for the central pixel
   ! Print some pahes diagnostics for the central pixel
   ! Needs:
   !   Call Alloc_Interferom_Pars  first
   !--------------------------------------------
   use constants, only : dp, pi, Sample, Refrac, c_mps
   use DataConstants, only : Time_Dim, Cnu_dim
   use Chunk_AntInfo, only : CTime_spectr, Ant_Stations, Ant_IDs, Ant_nr, Ant_pos !, Ant_RawSourceDist
   use Chunk_AntInfo, only : Start_time, TimeBase
   use FFT, only : RFTransform_CF !, RFTransform_CF2CT
   use StationMnemonics, only : Statn_ID2Mnem
   Use Interferom_Pars, only : IntfBase, IntfDim, IntfPhaseCheck, SumStrt, IntfLead, SumWindw
   Use Interferom_Pars, only : CTime_sum, Cnu, RTime, IntFer_ant, Nr_IntFer, IntfNorm, IntfNuDim
   Use Interferom_Pars, only : CenLoc, t_shft, t_offsetPow, tMin, tMax
   Use Interferom_Pars, only : i_chunk, IntPowSpec, NrPixSmPowTr, smooth, N_smth, RefSmPowTr
   Implicit none
   Integer, intent(in) :: i_eo
   integer :: i_loc(1), i_ant, j_IntFer, MLoc, Mloc_all, i, j, i_s
   Real(dp) :: angle, angle_all, HorAng, AziAng, D, RD, North, East, Heigh, AbsAmp, NRM, SmPow  ! elevation
   Real(dp) :: IntfNormMax=1.2
   Character (len=5) :: quadr
   Complex(dp) :: CNu_s(0:Cnu_dim)  ! CTime_s(1:Time_dim),
   !
   CTime_sum(:)=0.
   Do j_IntFer=1,Nr_Intfer
      i_ant=IntFer_ant(j_IntFer,1)
      !write(2,*) i_ant, j_IntFer, Ant_pos(:,i_ant,i_chunk)
      CTime_spectr(:,i_ant,i_chunk)=CTime_spectr(:,i_ant,i_chunk)/100.  ! to undo the factor that was introduced in antenna-read
      RTime(:)=REAL(CTime_spectr(IntfBase:IntfBase+IntfDim,i_ant,i_chunk))
      North=CenLoc(1)-Ant_pos(1,i_ant,i_chunk)
      East=CenLoc(2)-Ant_pos(2,i_ant,i_chunk)
      Heigh=CenLoc(3)-Ant_pos(3,i_ant,i_chunk)
      RD=North*North +East*East
      D=sqrt(RD+Heigh*Heigh)
      RD=sqrt(RD)
      AziAng=atan2(East,North)*180./pi
      HorAng=atan2(Heigh, RD)
      If(AziAng .lt. -45.) AziAng=AziAng+360.
      !write(2,*) AziAng, (AziAng+45.)/90.,INT((AziAng+45.)/90.),NINT(AziAng/90.)
      Select Case(INT((AziAng+45.)/90.))
         Case (1)
         Quadr='East'
         Case(2)
         Quadr='South'
         Case(3)
         Quadr='West'
         Case(0)
         Quadr='North'
         Case Default
         write(2,*) (AziAng+45.)/90.,INT((AziAng+45.)/90.),NINT(AziAng/90.)
         Quadr='???'
      End Select
      If(j_IntFer.eq.1) NRM=sin(HorAng)/D
      IntfNorm(j_IntFer)=NRM*D/sin(HorAng)
      If(IntfNorm(j_IntFer) .gt. IntfNormMax) IntfNorm(j_IntFer)=IntfNormMax
      HorAng=HorAng*180./pi
      !IntfNorm(j_IntFer)=1.  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !Call RFTransform_CF(RTime,Cnu(0,j_IntFer))
      Call RFTransform_CF(RTime,Cnu_s(0))
      Cnu(j_IntFer,0:IntfNuDim)=Cnu_s(0:IntfNuDim)*IntfNorm(j_IntFer)
      CTime_sum(:)=CTime_sum(:) + CTime_spectr(:,i_ant,i_chunk)*IntfNorm(j_IntFer)
      If(IntfPhaseCheck) Then
         If(j_IntFer.eq.1) Then
            i_loc=MaxLoc( ABS( CTime_sum(SumStrt:SumStrt+SumWindw) ))
            MLoc=i_loc(1)+SumStrt-1
            !write(2,*) 'mloc:', MLoc, i_loc(1), IntfBase, IntfDim, IntfLead, SumStrt, SumWindw
            i_loc=MaxLoc( ABS( CTime_sum(:) ))
            Mloc_all=i_loc(1)
         Endif
         AbsAmp=ABS(CTime_spectr(MLoc,i_ant,i_chunk))
         angle=atan2(IMAG(CTime_spectr(MLoc,i_ant,i_chunk)), REAL(CTime_spectr(MLoc,i_ant,i_chunk)))*180./pi
         angle_all=atan2(IMAG(CTime_spectr(MLoc_all,i_ant,i_chunk)), REAL(CTime_spectr(MLoc_all,i_ant,i_chunk)))*180./pi
         write(2,"(3x,A,2x,I3.3,I3.3,' Ampl=(',2F9.2,')=(A,phi;nrmA)=(',F9.1,',',F6.1, ';',F7.1, &
            A ,F6.1,'[km], (th,phi)=(',F5.1,',',F6.1,')',1x,A5)") &
            Statn_ID2Mnem(Ant_Stations(i_ant,i_chunk)),Ant_Stations(i_ant,i_chunk), Ant_IDs(i_ant,i_chunk), &
            CTime_spectr(MLoc,i_ant,i_chunk), AbsAmp, angle, AbsAmp*IntfNorm(j_IntFer), &
            '), Antenna @', D/1000.,HorAng,AziAng,Quadr
      EndIf
       !
   Enddo         ! Select the antennas for Interferometry
   !write(2,*) 'Cnu(Nr_IntFer/2,IntfNuDim/2)', IntfNuDim/2, IntFer_ant(Nr_IntFer/2), Cnu(Nr_IntFer/2,IntfNuDim/2) &
   !   ,  IntFer_ant(Nr_IntFer/2+1), Cnu(Nr_IntFer/2+1,IntfNuDim/2)
   CTime_sum(:)=CTime_sum(:)/Nr_IntFer
   !
   i_ant=IntFer_ant(1,1)
   t_shft=sqrt(SUM(CenLoc(:)*CenLoc(:)))*Refrac/c_mps ! in seconds due to signal travel distance
   write(2,*) 'Nr_IntFer:',Nr_IntFer,', ref. ant.= ', Ant_IDs(i_ant,i_chunk), Statn_ID2Mnem(Ant_Stations(i_ant,i_chunk)),&
      ', time difference with central pixel=',t_shft*1000.,'[ms]'
   If(IntfPhaseCheck) Then
      write(2,*) 'Source time highest peak in selected spectrum [ms]:',((Start_time(i_chunk)+MLoc)*sample-t_shft)*1000., &
         '[ms], @sample', MLoc
      write(2,*) 'Source time highest peak in whole chunk [ms]:',((Start_time(i_chunk)+MLoc_all)*sample-t_shft)*1000., &
         '[ms], @sample', MLoc_all
   EndIf
   !
   If(IntPowSpec) Then  ! Build summed power spectrum for ref antenna to check coherence
      Do i=1,NrPixSmPowTr ! N_smth+1,SumWindw-N_smth,N_smth
         i_s=1+i*N_smth
         SmPow=0.
         Do j=-N_smth,N_smth
            SmPow=SmPow+smooth(j)*(ABS(CTime_spectr(SumStrt+i_s+j,IntFer_ant(1,i_chunk),i_chunk)))**2
         Enddo
         RefSmPowTr(i)=SmPow
      Enddo
   EndIf
   !
   tMin=((Start_time(i_chunk)+SumStrt)*sample-t_shft)*1000.
   tMax=((Start_time(i_chunk)+SumStrt+SumWindw)*sample-t_shft)*1000.
   write(2,*) 'Source time window start@', tMin, '[ms], end@', tMax
   t_offsetPow=((Start_time(i_chunk)+SumStrt)*sample)*1000.-TimeBase
   !
   Return
End Subroutine SISetupSpec
!-----------------------------------------------
Subroutine SIAnalyzePixelTTrace(i_eo, i_N, i_E, i_h)
   ! Analyze the trace of this pixel
   !--------------------------------------------
   use constants, only : dp !, pi, Sample, Refrac, c_mps
   Use Interferom_Pars, only : CTime_pix, Nr_IntFer, SumWindw, IntfLead, PixLoc, polar, N_pix
   Use Interferom_Pars, only : IntPowSpec, PixelPower, MaxSmPow, N_smth, smooth
   Use Interferom_Pars, only : AveInten, AveIntenE, AveIntenN, RimInten, MaxIntfInten, MaxIntfIntenLoc
   Use Interferom_Pars, only : SlcInten, NrSlices, SliceLen, MaxSlcInten, MaxSlcIntenLoc
   Use Interferom_Pars, only : NrPixSmPowTr, MaxSmPowGrd, PixSmPowTr, RefSmPowTr
   Implicit none
   Integer, intent(in) :: i_eo, i_N, i_E, i_h
   integer :: i, j, i_s  !_loc(1), i_ant, j_IntFer, MLoc, Mloc_all
   Real(dp) :: IntfInten, SmPow
   !
   PixelPower(1:SumWindw)=ABS(CTime_pix(IntfLead+1:IntfLead+SumWindw))**2/Nr_IntFer**2
   !
   ! Fold intensity with smooth function and determine sequential positions at the peaks of the smoothing function
   If(IntPowSpec) Then
      Do i=1,NrPixSmPowTr ! N_smth+1,SumWindw-N_smth,N_smth
         i_s=1+i*N_smth
         SmPow=0.
         Do j=-N_smth,N_smth ! fold with smoothing function
            SmPow=SmPow+smooth(j)*PixelPower(i_s+j)
         Enddo
         !If( SmPow.gt. 0.6*RefSmPowTr(i) ) SmPow=0.
         PixSmPowTr(i,i_N, i_E, i_h)=SmPow
         If( SmPow.gt. MaxSmPow(i) ) Then
            MaxSmPow(i)=SmPow
            MaxSmPowGrd(1,i)=i_N
            MaxSmPowGrd(2,i)=i_E
            MaxSmPowGrd(3,i)=i_h
         EndIf
      Enddo
   EndIf
   ! Extract interesting numbers such as pixel intensity
   IntfInten=0.
   i=0
   SlcInten(:)=0.
   Do i_s=1,NrSlices
      Do j=1,SliceLen
         i=i+1
         SlcInten(i_s)=SlcInten(i_s) + PixelPower(i)
      Enddo
      SlcInten(i_s)=SlcInten(i_s)/SliceLen
      IntfInten = IntfInten + SlcInten(i_s)
      If(SlcInten(i_s) .gt. MaxSlcInten(i_s)) then  ! find max intensity location
         MaxSlcInten(i_s)=SlcInten(i_s)
         MaxSlcIntenLoc(i_s,:)=PixLoc(:)
      Endif
      SlcInten(i_s)=SlcInten(i_s)*SlcInten(i_s) ! /SliceLen
      AveInten(i_s,i_h)=AveInten(i_s,i_h) + SlcInten(i_s)
      AveIntenE(i_s,i_h)=AveIntenE(i_s,i_h) + i_E*SlcInten(i_s)
      AveIntenN(i_s,i_h)=AveIntenN(i_s,i_h) + i_N*SlcInten(i_s)
      If((i_E.eq.N_pix(2,1)) .or. (i_E.eq.N_pix(2,2)) .or. (i_N.eq.N_pix(1,1)) .or. (i_N.eq.N_pix(1,2))) Then ! rim of picture
         RimInten(i_s,i_h)=RimInten(i_s,i_h) + SlcInten(i_s)
      Endif
   Enddo
   IntfInten = IntfInten/NrSlices
!   If(IntfInten .gt. MaxIntfInten) then  ! find max intensity location
!      MaxIntfInten=IntfInten
!      MaxIntfIntenLoc(:)=PixLoc(:)
!   Endif
   !
   i=0  ! Store totals information
   PixSmPowTr(i,i_N, i_E, i_h)=IntfInten
   If( IntfInten .gt. MaxSmPow(i) ) Then
      MaxSmPow(i)=IntfInten
      MaxSmPowGrd(1,i)=i_N
      MaxSmPowGrd(2,i)=i_E
      MaxSmPowGrd(3,i)=i_h
      MaxIntfInten=IntfInten
      MaxIntfIntenLoc(:)=PixLoc(:)
   EndIf
   !
   Return
End Subroutine SIAnalyzePixelTTrace
