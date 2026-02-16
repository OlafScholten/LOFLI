module mod_test
contains
  function ilog2_b(val ) result( res )
    integer, intent(in) :: val
    integer             :: res
    integer             :: tmp

    res = -1
    ! Negativ values not allowed
    if ( val < 1 ) return

    tmp = val
    do while (tmp > 0)
      res = res + 1
      tmp = shiftr( tmp, 1 )
    enddo
  end function
subroutine test
 ! use mod_test
  print *,'Bitshifting: ', ilog2_b(12345)
  print *,'Formula:     ', floor( log(real(12345) ) / log(2.) )
end subroutine test
end module
!
!
Subroutine Pol2Carth(LocPol,LocCth)
   use constants, only : dp
   Implicit none
   Real(dp), Intent(in) :: LocPol(3)
   Real(dp), Intent(out) :: LocCth(3)
   LocCth(1)=LocPol(3)*cos(LocPol(2))*cos(LocPol(1))
   LocCth(2)=LocPol(3)*cos(LocPol(2))*sin(LocPol(1))
   LocCth(3)=LocPol(3)*sin(LocPol(2))
   Return
End Subroutine Pol2Carth
!
Subroutine Carth2Pol(LocCth,LocPol)
   use constants, only : dp
   Implicit none
   Real(dp), Intent(in) :: LocCth(3)
   Real(dp), Intent(out) :: LocPol(3)
   LocPol(3)=sqrt(SUM(LocCth(:)*LocCth(:)))  ! distance [m]    ;h
   LocPol(2)=asin(LocCth(3)/LocPol(3))    ! elevation angle=theta [radian]  ;E range: 0^o< th <80^o
   LocPol(1)=atan2(LocCth(2),LocCth(1))      ! phi [radian]    ;N
   Return
End Subroutine Carth2Pol
!------------------------
Subroutine GridPhaseChange(RMSW, IO_control)
! Prepatory work for PixBoundingBox
   ! Check max time shift for neighboring pixels; should ideally differ by 1 sample for neighboring pixels
   ! Calculate approximate bounding box of pixel image
   !  Needs
   !     call SelectIntfAntennas  first
   ! ------------------------------------
   use constants, only : dp
   use DataConstants, only : Time_Dim
   use Chunk_AntInfo, only : Ant_pos, Ant_RawSourceDist  ! CTime_spectr, Ant_Stations, Ant_IDs, Ant_nr,
   Use Interferom_Pars, only : i_chunk, IntFer_ant, Nr_IntFerMx, Nr_IntFerCh
   Use Interferom_Pars, only : CenLoc, CenLocPol, Polar, d_loc, Diff, N_pix
   Use Interferom_Pars, only : xMin, xMax, yMin, yMax, zMin, zMax
   !Use Interferom_Pars, only : SumStrt, SumWindw, IntfDim, IntfNuDim, IntfLead, IntfBase
   Implicit none
   Real(dp), intent(out) :: RMSW(1:3)
   Logical, intent(in) :: IO_control
   Real(dp) :: PixLoc(1:3), PixLocPol(1:3), t_shft, RDist
   Real(dp) :: t_n, t_p, mean, RMS, meanW, Weight, D, W, AW2
   Integer :: i, i_ant, j_IntFer, Nr_IntFer
   !
   PixLoc(:)=CenLoc(:)
   xMin=+99999 ; xMax=-99999; yMin=+99999; yMax=-99999; zMin=+99999; zMax=-99999 !; tMin, tMax,
   !Err(:)=0.
   Nr_IntFer=Nr_IntFerCh(i_chunk)
   Do i=1,3
      If(polar) then
         PixLocPol(:)=CenLocPol(:)
         PixLocPol(i)=CenLocPol(i)+d_loc(i)
         Call Pol2Carth(PixLocPol,PixLoc)
      else
         PixLoc(:)=CenLoc(:)
         PixLoc(i)=CenLoc(i)+d_loc(i)
      Endif
      !write(2,*) 'PixLoc-carthesian',PixLoc(:)
      t_n=0.   ; t_p=0
      Mean=0. ; RMS=0.  ; AW2=0.
      MeanW=0. ; RMSW(i)=0. ; Weight=0
      Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
         i_ant=IntFer_ant(j_IntFer,i_chunk)
         Call RelDist(PixLoc(1),Ant_pos(1,i_ant,i_chunk),RDist)
         t_shft=Rdist - Ant_RawSourceDist(i_ant,i_chunk)    ! units of samples
         if(t_shft.lt.t_n) t_n=t_shft
         if(t_shft.gt.t_p) t_p=t_shft
         D=sqrt(sum( (PixLoc(1:3)-Ant_pos(1:3,i_ant,i_chunk))**2 ))
         !Weight=Weight + 1./D
         !MeanW=MeanW+t_shft/D
         !RMSW(i)=RMSW(i) + t_shft*t_shft/(D*D)
         ! 1 Jan 2026: more correctly
         W=1/D       ! may also be 1/D^2
         Weight=Weight + W
         MeanW=MeanW + t_shft*W
         RMSW(i)=RMSW(i) + t_shft*t_shft*W
         W=1.  !  1./(D*D)       ! may also be 1/D^2
         AW2=AW2+ W
         Mean=Mean+t_shft*W
         RMS=RMS + t_shft*t_shft*W
      EndDo ! j_IntFer
      Mean=Mean/AW2
      RMS=RMS/AW2
      RMS=sqrt(RMS-Mean*Mean)
      !MeanW=MeanW/Nr_IntFer
      !Weight=Weight/Nr_IntFer
      !RMSW(i)=RMSW(i)/Nr_IntFer
      !RMSW(i)=5.*sqrt(RMSW(i)-MeanW*MeanW)/Weight  ! Factor is bandwidth*sample_time=10 10^6 [1/s] * 5 10^-9 [s] *100%=5 %
      MeanW=MeanW/Weight
      RMSW(i)=RMSW(i)/Weight
      RMSW(i)=5.*sqrt(RMSW(i)-MeanW*MeanW)  ! Factor is bandwidth*sample_time=10 10^6 [1/s] * 5 10^-9 [s] *100%=5 %
      Weight=Weight/Nr_IntFer
      If(IO_control) write(2,"(A,i2,2f7.2,A,f6.2,f7.3,A,f9.5,A,3f7.2,A)") 'min & max time shift [samples] per pixel',i, &
         t_n, t_p,', RMS=', RMS, RMSW(i), '%, for d(i)=',d_loc(i),', d(N,E,h)=',Pixloc(:)-CenLoc(:),'[m]'
      diff(i) =t_p      ! needed to calculate lead-time for complete grid
      !Err(:)=Err(:) + (Pixloc(:)-CenLoc(:))**2/(RMS*RMS)
      xMin=min(xmin,CenLoc(2)+N_pix(i,1)*abs(CenLoc(2)-PixLoc(2))); xMax=max(xMax,CenLoc(2)+N_pix(i,2)*abs(CenLoc(2)-PixLoc(2)))
      yMin=min(ymin,CenLoc(1)+N_pix(i,1)*abs(CenLoc(1)-PixLoc(1))); yMax=max(yMax,CenLoc(1)+N_pix(i,2)*abs(CenLoc(1)-PixLoc(1)))
      zMin=min(zmin,CenLoc(3)+N_pix(i,1)*abs(CenLoc(3)-PixLoc(3))); zMax=max(zMax,CenLoc(3)+N_pix(i,2)*abs(CenLoc(3)-PixLoc(3)))
      If(-t_n.gt.t_p) diff(i) =-t_n
      PixLoc(i)=CenLoc(i)  ! set back to center for Cartesian coordinates
   Enddo ! i
   If(zMin.lt.0.) zMin=0.
End Subroutine GridPhaseChange
!------------------------
Subroutine PixBoundingBox(GridVolume, BoxFineness)
   ! Check max time shift for neighboring pixels; should ideally differ by 1 sample for neighboring pixels
   ! Calculate approximate bounding box of pixel image
   ! Setup arrays
   !  Needs
   !     call SelectIntfAntennas  first
   ! ------------------------------------
   use constants, only : dp,pi
   use DataConstants, only : Time_Dim
   use Chunk_AntInfo, only : Ant_pos, Ant_RawSourceDist  ! CTime_spectr, Ant_Stations, Ant_IDs, Ant_nr,
   Use Interferom_Pars, only : i_chunk, IntFer_ant, Nr_IntFerMx, Nr_IntFerCh
   Use Interferom_Pars, only : CenLoc, CenLocPol, Polar, d_loc, Diff, N_pix
   Use Interferom_Pars, only : xMin, xMax, yMin, yMax, zMin, zMax, FirstTimeInterf
   Use Interferom_Pars, only : SumStrt, SumWindw, IntfDim, IntfNuDim, IntfLead, IntfBase
   Implicit none
   !Integer, intent(in) :: i_chunk
   Real(dp), intent(in) :: GridVolume(1:3), BoxFineness
   Real(dp) :: PixLoc(1:3), PixLocPol(1:3), t_shft, RDist
   Real(dp) ::  RMSW(1:3) !  t_n, t_p,  RMS, Err(3),  Weight, D
   Integer :: i, i_ant, j_IntFer, Nr_IntFer
   logical :: IO_control
   !
   !write(2,*) 'Error ellips (N,E,h)=',sqrt(Err(:)),' [m]'
   If((GridVolume(1)*GridVolume(2)*GridVolume(3) .gt. 0.1) .and. (BoxFineness .gt. 0.01)) Then
      If(FirstTimeInterf) Write(2,*) 'Fineness=',BoxFineness,', for a grid volume of:',GridVolume(1:3)
      IO_control=.false.
      !IO_control=.true.
      d_loc(1:3)=1.
      Call GridPhaseChange(RMSW, IO_control)
      d_loc(1:3)=BoxFineness/RMSW(1:3)
      N_pix(1:3,2)=NINT(GridVolume(1:3)/d_loc(1:3))
      !write(2,*) '!PixBoundingBox:', BoxFineness, GridVolume(:), d_loc(1:3), N_pix(1:3,1)
      !Flush(unit=2)
   EndIf
   !
   !d_loc(1)=d_N  ;  d_loc(2)=d_E   ;  d_loc(3)=d_h
   N_pix(:,1)=-N_pix(:,2)
   If(polar) then
      write(2,*) 'Polar coordinates used for grid!!'
      d_loc(1)=d_loc(1)*pi/180.   ! convert to radian
      d_loc(2)=d_loc(2)*pi/180.   ! convert to radian
   Endif
   ! ----------------------------------------------------
   ! Extract polar coordinates;  1=N=phi; 2=E=theta; 3=h=Radial distance
   ! Check limits on ranges      N_hlow, N_Eup, N_Elow
   If(FirstTimeInterf) write(2,*) 'Central Carthesian coordinates (N,E,h)=',CenLoc(:)
   If(polar) then
      Call Carth2Pol(CenLoc,CenLocPol)
      !CenLocPol(3)=sqrt(SUM(CenLoc(:)*CenLoc(:)))  ! distance [m]    ;h
      !CenLocPol(2)=asin(CenLoc(3)/CenLocPol(3))    ! elevation angle=theta [radian]  ;E range: 0^o< th <80^o
      !CenLocPol(1)=atan2(CenLoc(2),CenLoc(1))      ! phi [radian]    ;N
      If(FirstTimeInterf) write(2,*) 'Central polar coordinates (R,th,ph)=',CenLocPol(3),CenLocPol(2)*180./pi,CenLocPol(1)*180/pi
      ! check ranges:
      If(CenLocPol(3)+N_pix(3,1)*d_loc(3) .lt. 0) N_pix(3,1)= 1-CenLocPol(3)/d_loc(3) ! negative number generally ! 1-: not to have it ridiculously close
      If(CenLocPol(2)+N_pix(2,1)*d_loc(2) .lt. 0) N_pix(2,1)= -CenLocPol(2)/d_loc(2) ! negative number generally
      If(CenLocPol(2)+N_pix(2,2)*d_loc(2) .gt. pi*4/9.) N_pix(2,2)= (pi*4/9.-CenLocPol(2))/d_loc(2) ! pos number generally
      If(FirstTimeInterf) Then
         write(2,*) 'distance range:',CenLocPol(3)+N_pix(3,1)*d_loc(3),CenLocPol(3)+N_pix(3,2)*d_loc(3)
         write(2,*) 'Elevation angle range:',(CenLocPol(2)+N_pix(2,1)*d_loc(2))*180./pi,(CenLocPol(2)+N_pix(2,2)*d_loc(2))*180./pi
         write(2,*) 'Azimuth angle range:', (CenLocPol(1)+N_pix(1,1)*d_loc(1))*180./pi, (CenLocPol(1)+N_pix(1,2)*d_loc(1))*180./pi
      EndIf
   Else
      If(CenLoc(3)+N_pix(3,1)*d_loc(3) .lt. 0) N_pix(3,1)= -CenLoc(3)/d_loc(3) ! negative number generally
   Endif
   !N_pix(2,1)=N_Elow   ;  N_pix(3,1)=N_hlow  ;  N_pix(2,2)=N_Eup
   If(FirstTimeInterf) write(2,*) 'N_hlow, N_Eup, N_Elow',N_pix(:,1), N_pix(:,2)
   !
   !------------------------------------
   !
   Call GridPhaseChange(RMSW, FirstTimeInterf)
   !
   !-----------------------------------------------
   !  setup arrays
   IntfLead=10+(diff(1)*N_pix(1,2)+diff(2)*N_pix(2,2)+diff(3)*N_pix(3,2)) ! Estimate lead and trail buffer length,
   IntfNuDim=SumWindw/2 + IntfLead +1 ! Length of frequency trace used for interference
   i = shiftr( IntfNuDim, 9 )+1
   IntfNuDim=shiftl( i, 9 )  ! Make sure it approches some multiple of 2
   !write(2,"(A,o10,o10,i7)") 'i:',i,IntfNuDim,IntfNuDim
   IntfDim=2*IntfNuDim ! Length of time trace used for interference
   IntfLead=(IntfDim-SumWindw)/2
   IntfBase=SumStrt - IntfLead
   !write(*,"(A,o10,b20)") '(lead,dim):', IntfLead, IntfDim
   If((IntfBase+IntfDim .gt. Time_dim) .or. (IntfBase .lt. 1)) then
      write(2,*) 'dimensions for interferometry out of range: ', IntfBase, IntfBase+IntfDim, Time_dim
      write(2,*) 'input for interferometry: ', SumStrt, SumWindw
      Flush(unit=2)
      stop 'Intferom dimension problem'
   Endif
   !
   !write(2,*) diff(:)
   !
   Return
End Subroutine PixBoundingBox
!-----------------------------------------------
Subroutine OutputIntfPowrTotal(RefAntSAI, i_eo)
   ! Output power info for complete trace as function of pixel
   !--------------------------------------------
   use constants, only : dp, pi, Sample
   use DataConstants, only : DataFolder, OutFileLabel, Time_Dim
   !use Chunk_AntInfo, only : TimeBase
   use Chunk_AntInfo, only : CTime_spectr, StartT_sam
   Use Interferom_Pars, only : CTime_sum, SumStrt, SumWindw, polar, N_pix, d_loc, IntfLead, CenLocPol, CenLoc, RimSmPow
   Use Interferom_Pars, only : PixelPower, N_smth, FirstTimeInterf ! , MaxSmPow
   Use FitParams, only : AntennaRange
   Use ThisSource, only : Dual
   !Use Interferom_Pars, only : xMin, xMax, yMin, yMax, zMin, zMax, tMin, tMax
   Use Interferom_Pars, only : i_chunk, MaxSmPowGrd, PixSmPowTr, t_shft, NewCenLoc
   Implicit none
   Integer, intent(in) :: RefAntSAI, i_eo
   integer :: i, j, i_N, i_E, i_h, SSm=5 ! Resum Hilbert envelope
   character(len=4) :: txt
   character(len=2) :: txt1
   Real(dp) :: SMPowMx, d_Mx(1:3), QualMx(10), SMPowBar, d_Bar(1:3), QualBar(10)
   Real(dp) :: Pref, Psum, half
   Real(dp) :: PixLocPol(1:3), PixLoc(1:3), StartTime_ms
   Logical :: Obsolete=.false.
   !
   !  write intensity distribution for pixel-planes with i_h=-1,0,+1 needed for contour plots
   If(Dual) then
      txt1='d_'
   Else
      Write(txt1,"(i1,'_')") i_eo
   Endif
   If(Obsolete) Then
      i_h=0
      !write(2,*) 'started:',txt
      !flush(unit=2)
      txt=txt1//'NE'
      OPEN(UNIT=30,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'Interferometer'//txt//'.csv')
      Do i_N=N_pix(1,1),N_pix(1,2)   ! or Phi  ! Loop over pixels
         Do i_E=N_pix(2,1),N_pix(2,2)   ! or elevation angle
            If(polar) then
               write(30,"(i5,',',i5,',',g12.4)") i_N,i_E, PixSmPowTr(0,i_N,i_E,i_h) ! write the sum of the square of the Hilbert transform for this pixel
            Else
               write(30,"(i5,',',i5,',',g12.4)") i_N,i_E, PixSmPowTr(0,i_N,i_E,i_h) ! write the sum of the square of the Hilbert transform
            Endif
         Enddo
      Enddo
      Close(UNIT=30)
      i_N=0
      txt=txt1//'Eh'
      OPEN(UNIT=30,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'Interferometer'//txt//'.csv')
      Do i_E=N_pix(2,1),N_pix(2,2)   ! or Phi  ! Loop over pixels
         Do i_h=N_pix(3,1),N_pix(3,2)   ! or elevation angle
            If(polar) then
               write(30,"(i5,',',i5,',',g12.4)") i_E,i_h, PixSmPowTr(0,i_N,i_E,i_h) ! write the sum of the square of the Hilbert transform for this pixel
            Else
               write(30,"(i5,',',i5,',',g12.4)") i_E,i_h, PixSmPowTr(0,i_N,i_E,i_h) ! write the sum of the square of the Hilbert transform
            Endif
         Enddo
      Enddo
      Close(UNIT=30)
      i_E=0
      txt=txt1//'Nh'
      OPEN(UNIT=30,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'Interferometer'//txt//'.csv')
      Do i_N=N_pix(1,1),N_pix(1,2)   ! or Phi  ! Loop over pixels
         Do i_h=N_pix(3,1),N_pix(3,2)   ! or elevation angle
            If(polar) then
               write(30,"(i5,',',i5,',',g12.4)") i_N,i_h, PixSmPowTr(0,i_N,i_E,i_h) ! write the sum of the square of the Hilbert transform for this pixel
            Else
               write(30,"(i5,',',i5,',',g12.4)") i_N,i_h, PixSmPowTr(0,i_N,i_E,i_h) ! write the sum of the square of the Hilbert transform
            Endif
         Enddo
      Enddo
      Close(UNIT=30)
   Else
      i=0
      Call WriteBeamingIntensities(i, txt1)
   EndIf
   !write(2,*) 'Done:',txt
   !flush(unit=2)
   !
   !---------------
   ! Do maximum pixel analysis for complete trace
   i=0
   Call FindInterpolMx(i,d_Mx, SMPowMx,QualMx)
   !
   If(FirstTimeInterf) Then
      write(2,*) 'on first pass calculate bary center complete trace', FirstTimeInterf
      If(SMPowMx.gt.0.) then  ! check max not at rim
         If(polar) then
            PixLocPol(:)=CenLocPol(:)+d_Mx(:)*d_loc(:)
            Call Pol2Carth(PixLocPol,PixLoc)
         else
            PixLoc(:)=CenLoc(:)+d_Mx(:)*d_loc(:)
            Call Carth2Pol(PixLoc,PixLocPol)
         Endif
         NewCenLoc(:)=PixLoc(:)
         write(2,"('(N,E,h)=',3(f11.5,','),' Mx=',g13.6,', Q=',3g13.6,', grid=',3F7.2)")   &
            PixLoc(:)/1000., SMPowMx, QualMx(1:3), d_Mx(1:3) !, BarySet
         !
         ! Use the barycenter to interpolate in distance and use only those pixels that have a large intensity
         Call FindBarycenter(i,d_Bar,SMPowBar,QualBar)
         If(SMPowBar.gt.0.) then
            If(polar) then
               PixLocPol(:)=CenLocPol(:)+d_Bar(:)*d_loc(:)
               Call Pol2Carth(PixLocPol,PixLoc)
            else
               PixLoc(:)=CenLoc(:)+d_Bar(:)*d_loc(:)
               Call Carth2Pol(PixLoc,PixLocPol)
            Endif
            !
            write(2,"('(N,E,h)=',3(f11.5,','),' Bar=',g13.6,', Q=',2g13.6,12x,', grid=',3F7.2)")   &
               PixLoc(:)/1000., SMPowBar, QualBar(1), QualBar(2), d_Bar(1:3)
         EndIf
      Else
         Write(2,*) 'Borderline case, Max @',MaxSmPowGrd(:,i)
      EndIf
   EndIf
   ! =======================================================================
   !--------------------------------------------
   ! output powerspectrum, resummed
   !SSm=5  ! Resum Hilbert envelope
   i=0   ;  If(polar) i=1 ! 0=false; 1=true
   half=(SSm+1.)/2.
   StartTime_ms=StartT_sam(1)*sample*1000.d0  ! in ms
   If(Dual) then
      OPEN(UNIT=30,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'_EISpeceo.dat')
      Write(30,"(f9.3,3(',',f8.3),',',f5.1,',',i2,',',i2,F8.3,g12.4,I4, ' 0')") StartTime_ms, CenLoc/1000., AntennaRange, &
         SSm, i, t_shft*1000., ABS(SMPowMx), N_smth
      Do i=0,Time_Dim/SSm-1
         Pref=0.
         Psum=0.
         Do j=1,SSm
            Pref=Pref + (ABS(CTime_spectr(j+i*SSm,RefAntSAI,i_chunk)))**2  !!/100.  ! to undo the factor that was introduced in antenna-read
            Psum=Psum + (ABS(CTime_spectr(j+i*SSm,RefAntSAI+1,i_chunk)))**2
         Enddo
         write(30,"(f8.3,',',g12.4,',',g12.4)") (i*SSm+half)/200., Pref*1e-4/SSm, Psum*1e-4/SSm
      Enddo
      Close(unit=30)
   Else
      Write(txt,"('_',i1)") i_eo
      OPEN(UNIT=30,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'IntfSpec'//TRIM(txt)//'.csv')
      Write(30,"(f9.3,3(',',f8.3),',',f5.1,',',i2,',',i2,F8.3,F9.2,' 0')") StartTime_ms, CenLoc/1000., AntennaRange, SSm, i, &
         t_shft*1000., SMPowMx
      Do i=0,Time_Dim/SSm-1
         Pref=0.
         Psum=0.
         Do j=1,SSm
            Pref=Pref + (ABS(CTime_spectr(j+i*SSm,RefAntSAI,i_chunk)))**2
            Psum=Psum + ABS(CTime_sum(j+i*SSm))**2
         Enddo
         write(30,"(f8.3,',',g12.4,',',g12.4)") (i*SSm+half)/200., Pref/SSm, Psum/SSm
      Enddo
      Close(unit=30)
   Endif
   !
   !--------------------------------------------------------------
   !
   Return
End Subroutine OutputIntfPowrTotal
!-----------------------------------------------
Subroutine OutputIntfPowrMxPos(i_eo)
   ! Write to file the positions on the maximum intensity in the interferometry plots
   ! Needs
   !     File 30 opened
   !--------------------------------------------
   use constants, only : dp, pi, Sample,  c_mps
   use DataConstants, only : DataFolder, OutFileLabel, FlashName
   use Chunk_AntInfo, only : TimeBase
   Use Interferom_Pars, only : SumStrt, SumWindw, polar, N_pix, d_loc, IntfLead, CenLocPol, CenLoc, RimSmPow
   Use Interferom_Pars, only : PixLoc, PixelPower, MaxSmPow, MaxSmPowLoc, N_smth, Nr_IntFerMx, Nr_IntFerCh, i_chunk, IntfNuDim
   Use Interferom_Pars, only : MaxSmPowI12, MaxSmPowQ, MaxSmPowU, MaxSmPowV
   Use Interferom_Pars, only : MaxSmPowI3, MaxSmPowU1, MaxSmPowV1, MaxSmPowU2, MaxSmPowV2
   Use Interferom_Pars, only : PolBasis, NGridInterpol, RMSGridInterpol, NGridExtrapol, RefinePolarizObs, TRIDFile
   Use Interferom_Pars, only : xMin, xMax, yMin, yMax, zMin, zMax, tMin, tMax, AmpltPlot, t_offsetPow
   Use Interferom_Pars, only : NrPixSmPowTr, MaxSmPowGrd, PixSmPowTr, RatMax, PixPowOpt, FirstTimeInterf
   Use Interferom_Pars, only : StI, StI12, StQ, StU, StV, StI3, StU1, StV1, StU2, StV2, P_un, P_lin, P_circ
   Use Interferom_Pars, only : dStI, dStI12, dStQ, dStU, dStV, dStI3, dStU1, dStV1, dStU2, dStV2, Chi2pDF
   Use Interferom_Pars, only : PolZen, PolAzi, PolMag, PoldOm, Stk_NEh  !  calculated in "PolTestCath"
   use Chunk_AntInfo, only : NoiseLevel
   use ThisSource, only : Dual
   use GLEplots, only : GLEplotControl
   Implicit none
   Integer, intent(in) :: i_eo
   integer :: i_slice, i_1, i_2, i_3, i_gr(1:3), k
   character(len=5) :: txt
   Real(dp) :: PixLocPol(1:3), t_ms, t_shft, y0
   Integer :: NBar, NMx
   Real(dp) :: SMPowMx, d_Mx(1:3), QualMx(10), SMPowBar, d_Bar(1:3), QualBar(10)
   Real(dp) :: s, X, Y, AngOff
   Real(dp) :: HorVec(1:3), VerVec(1:3)
   !Complex(dp), allocatable :: CMTime_Pix(:,:)
   Real(dp), external :: tShift_ms   !
   Character(len=30) :: Qual
   !
   Write(txt,"('_',i1)") i_eo  ! should not occur, obsolete option
   If(Dual) Then
      txt='x'
      !Allocate( CMTime_pix(1:2*IntfNuDim,1:3) )
   EndIf
   !  Calculate boundingBox
   If((xMax-xMin) .gt. (yMax-yMin)) Then
      y0=(yMax+yMin)/2.
      yMax=y0+(xMax-xMin)/2.
      yMin=y0-(xMax-xMin)/2.
   Else
      y0=(xMax+xMin)/2.
      xMax=y0+(yMax-yMin)/2.
      xMin=y0-(yMax-yMin)/2.
   EndIf
   !write(2,*) xMin/1000.-.05, xMax/1000.+.05, yMin/1000.-.05, yMax/1000.+.05, zMin/1000.-.05, zMax/1000.+.05
   !
   !If(Dual .and. RefinePolarizObs) Then
   !
   If(TRIDFile) Then
   !  Storing data for later processing in DataSelect:
      OPEN(UNIT=28,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'TRID'//TRIM(txt)//'.csv')
      write(28,"(2F11.3,4F10.2,1x,A,I3,L2,I3,I5,F7.2' 0')") tMin-.0005-TimeBase, tMax+.0005-TimeBase, &
         TimeBase,  CenLoc(:), TRIM(FlashName)//':'//TRIM(OutFileLabel), PixPowOpt, RefinePolarizObs, &
         N_smth, NrPixSmPowTr,NoiseLevel  ! gleRead:  xMin xMax yMin yMax zMin zMax tMin tMax ZmBx$ t_offst
      Write(28,"(A,20x,A )") '!   nr, Src_time[ms] ;       Location (N,E,h) [km]            Intens  ; chi^2    StI ', &
         'Stk_NEh(1,1:3), Stk_NEh(2,2:3), Stk_NEh(3,3)'
 !   nr, Src_time[ms] ;       Location (N,E,h) [km]            Intens  ; chi^2    StI  Stk_NEh(1,1:3), Stk_NEh(2,2:3), Stk_NEh(3,3)
!    1,  690.710076,   -7.66095,  -32.94734,    7.82216,         21.6,   1.40,   37.98     , (   4.904    ,  1.1848E-08) , (  -2.853    ,   4.850    ) , (   1.146    ,  0.6805    ) , (   24.08    , -1.7859E-08) , (   6.393    ,  -8.329    ) , (   9.000    , -5.3593E-08)
!               i_slice, t_ms-t_shft, PixLoc(1:3)/1000., SMPowMx, Chi2pDF, StI,   Stk_NEh(1,1:3), Stk_NEh(2,2:3), Stk_NEh(3,3)
      Write(2,*) 'Source info written to file:', &
            trim(DataFolder)//TRIM(OutFileLabel)//'TRID'//TRIM(txt)//'.csv'
      write(2,*) 'Header:tMin-.0005-TimeBase, tMax+.0005-TimeBase,  TimeBase,  CenLoc(:), ', &
            'TRIM(OutFileLabel)//TRIM(txt), PixPowOpt, 0.0, N_smth; data range=', 1, (1+NrPixSmPowTr)*N_smth
      write(2,"(2F11.3,4F10.2,1x,A,I3,F7.1,I3,' 0')") tMin-.0005-TimeBase, tMax+.0005-TimeBase, &
         TimeBase,  CenLoc(:),TRIM(OutFileLabel)//TRIM(txt), PixPowOpt, 0.0, N_smth  ! gleRead:  xMin xMax yMin yMax zMin zMax tMin tMax ZmBx$ t_offst
     ! OPEN(UNIT=27,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'IntfSpecPolAng'//TRIM(txt)//'.dat')
     ! Write(27,"(A,T19,A,T30,A,T38,A,T48,A, T60,A,T132,A)") '! slice#, samp# ;', &
     !    'Intensty','Zenith', 'Azimuth','domg','Twice repeated','time[ms] ;, Stk_NEh (N,N) (N,E) (N,h) (E,E) (E,h) (h,h)'
     ! OPEN(UNIT=26,STATUS='unknown',ACTION='WRITE', &
     !    FILE=trim(DataFolder)//TRIM(OutFileLabel)//'IntfSpecCartStokes'//TRIM(txt)//'.dat')
     ! Write(26,"(A,T19,A)") '! slice#, samp# ;', 'time[ms] ;, Stk_NEh (N,N) (N,E) (N,h) (E,E) (E,h) (h,h)'
  ! Else IF(.not. Dual) Then
  !    OPEN(UNIT=28,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'IntfSpecPowBar'//TRIM(txt)//'.dat')
  !    write(28,"(6F8.2,2F9.3,A,F7.1,' 0')") xMin/1000.-.005, xMax/1000.+.005, yMin/1000.-.005, yMax/1000.+.005, &
  !       zMin/1000.-.05, zMax/1000.+.05, tMin-.0005-TimeBase, tMax+.0005-TimeBase, ' IntfBox ', TimeBase         ! gleRead:  xMin xMax yMin yMax zMin zMax tMin tMax ZmBx$ t_offst
  !    Write(28,"(A,I7,A,A,F8.3,A,F6.3)") &
  !       '0 ',NrPixSmPowTr,' 0 0 ',TRIM(OutFileLabel)//TRIM(txt),AmpltPlot,' 0 0 0 0 0 0 0 ',NoiseLevel ! gleRead:  NTracks EventNr Q t_start label$ AmplitudePlot a1 b1 c1 a2 b2 c2 d2
  ! EndIf
  ! OPEN(UNIT=29,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'IntfSpecPowMx'//TRIM(txt)//'.dat')
   ! Data for results plots of this imager run
      write(Qual,"(A,F6.2,A)") '"I>',NoiseLevel,'"'
      OPEN(UNIT=29,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'TRID'//TRIM(txt)//'.plt')
      write(29,"(6F8.2,2F9.3,A,F7.1,i3,' 0')") yMin/1000.-.005, yMax/1000.+.005, xMin/1000.-.005, xMax/1000.+.005, &
         zMin/1000.-.005, zMax/1000.+.005, tMin-.0005-TimeBase, tMax+.0005-TimeBase, ' IntfBox ', TimeBase, PixPowOpt ! gleRead:  xMin xMax yMin yMax zMin zMax tMin tMax ZmBx$ t_offst
      Write(29,*) '0 ',NrPixSmPowTr,TRIM(Qual),N_smth, TRIM(FlashName)//':'//TRIM(OutFileLabel), AmpltPlot, &
         ' 0 0 0 0 0 0 0 ',NoiseLevel,  '1.0 !' ! gleRead:  NTracks EventNr Q t_start label$ AmplitudePlot a1 b1 c1 a2 b2 c2 d2
      ! True bounding Box
      Do i_1=1,2   ! or azimuth Phi
         Do i_2=1,2   ! or elevation angle theta
            Do i_3=1,2   ! or Distance R
               i_gr(1)=N_pix(1,i_1) ; i_gr(2)=N_pix(2,i_2) ; i_gr(3)=N_pix(3,i_3)
               If(polar) then
                  PixLocPol(:)=CenLocPol(:)+i_gr(:)*d_loc(:)
                  Call Pol2Carth(PixLocPol,PixLoc)
               Else
                  PixLoc(:)=CenLoc(:)+i_gr(:)*d_loc(:)
               Endif
               Write(29,"(1x,3F9.3,2x,3I3,' ;')") PixLoc(1:3)/1000.,i_1,i_2,i_3 ! (N,E,h)=(x,y,z)
            Enddo
         Enddo
      Enddo
     !    Write(29,"(A,9x,A,13x,';   ', 9(A,'/I [%]   ;   '),A )") &
     !       '!   nr,  time[ms] ;','StI','StI12', '  StQ','  StU','  StV',' StI3',' StU1',' StV1',' StU2',' StV2', &
     !       'Chi^2/DoF, P_un, P_lin, P_circ [%]'
      Write(29,"(A,7x,A )") '!   nr, Src_time[ms] ;       Location (N,E,h) [km]           Intens; chi^2 sampl    StI', &
         'I12[%], I3    dI3, P_Un P_Lin P_Circ    (Magn(k),    Zen(k), Azi(k), dOm(k),    k=1,3)'
      !
   EndIf
   !
   ! Analyze Smoothed-power positions
   NMx=0
   NBar=0
   !
   ! Setup for angle linear polarization
   ! Vec is system with (N,E,h)
   !  VerVec points up as seen from reference antenna
   !  HorVec points to the right,
   !  VerVec x HorVec point in radial direction away from RefAnt
   !   Use Interferom_Pars, only : alpha, PolBasis, PowerScale
   ! First check to what extent PolBasis(:,3) is indeed in the radial direction from RefAnt
   s=SUM(CenLoc(:)*CenLoc(:))  ; X=SUM(PolBasis(:,3)*CenLoc(:))/sqrt(s)
   If(FirstTimeInterf) Write(2,"(A,F8.3,A)") 'overlap radial distance and longitudinal pol=',X*100.,'%'
   HorVec(1)=-CenLoc(2)  ; HorVec(2)=CenLoc(1) ; HorVec(3)=0.
   s=HorVec(1)*HorVec(1)+ HorVec(2)*HorVec(2)  ;  HorVec(:)=HorVec(:)/sqrt(s)
   VerVec(1)=CenLoc(2)*HorVec(3)- CenLoc(3)*HorVec(2) !
   VerVec(2)=CenLoc(3)*HorVec(1)- CenLoc(1)*HorVec(3)
   VerVec(3)=CenLoc(1)*HorVec(2)- CenLoc(2)*HorVec(1)
   s=VerVec(1)*VerVec(1)+ VerVec(2)*VerVec(2) + VerVec(3)*VerVec(3)  ;  VerVec(:)=VerVec(:)/sqrt(s)
!      s=SUM(ev(:,k)*ev(:,k))
!      ev(:,k)=ev(:,k)/sqrt(s)
   ! get angle of PolBasis(:,1) with VerVec(:)
   X=SUM(PolBasis(:,1)*VerVec(:)) ; Y=SUM(PolBasis(:,1)*HorVec(:)) ; AngOff=ATAN2(Y,X)
   ! N.B. PolBasis is by construction lefthanded in the real world
   !
   !write(2,"(8x,A,4x,A)") "i, AveSmPow/y0, RimMx/y0,  Mx,  volume;      Barymetric pixel ;    Maximal pixel", &
   !   ';  Closest edge'
   If(FirstTimeInterf) write(2,"(6x,A, 4x,A,9x,';', 4x,A, 4x,A,I5)") "i, AveSmPow/y0, RimMx/y0,  Mx;", &
      "t[ms], Position=(N,E,h)[km]","Q/I, U/I, V/I [%] ;" , "I3, Linpol*100., angle ; nr slices=", NrPixSmPowTr
   !If(FirstTimeInterf) write(2,*) '!OutputIntfPowrMxPos;NrPixSmPowTr:',NrPixSmPowTr
   MaxSmPowLoc(:,:) = 0.0 ! finite value signal that it was not a border or below noise case
   Do i_slice=1,NrPixSmPowTr ! N_smth+1,SumWindw-N_smth,N_smth
      t_ms=(1+i_slice*N_smth)*sample*1000.d0+t_offsetPow  ! time stamp
!!!!!!!!!   t_shft=sqrt(SUM(CenLoc(:)*CenLoc(:)))*Refrac/c_mps ! in seconds due to signal travel distance
!!!!!!!!!!!!!   t_offsetPow=((StartT_sam(i_chunk)+SumStrt)*sample-t_shft)*1000.-TimeBase
      !
      ! Do a 3D parabolic fit to interpolate around the maximal pixel
      Call FindInterpolMx(i_slice,d_Mx, SMPowMx,QualMx)
      !
      If(polar) then
         PixLocPol(:)=CenLocPol(:)+d_Mx(:)*d_loc(:)
         Call Pol2Carth(PixLocPol,PixLoc)
      else
         PixLoc(:)=CenLoc(:)+d_Mx(:)*d_loc(:)
         Call Carth2Pol(PixLoc,PixLocPol)
      Endif
      MaxSmPowLoc(1:3,i_slice)=PixLoc(:) ! SMPowMx=-1 When at border
      t_shft=tShift_ms(PixLoc(:))  ! sqrt(SUM(PixLoc(:)*PixLoc(:)))*1000.*Refrac/c_mps ! in seconds due to signal travel distance
      !write(29,"(i6,',',4(f11.5,','),4g13.6,',',3g13.5,2x,L1)") i, t_ms-t_shft, &
      !write(2,*) '!OutputIntfPowrMxPos;SMPowMx:',i_slice, SMPowMx, MaxSmPow(i_slice), NoiseLevel
      If(SMPowMx.lt.NoiseLevel) Then
         MaxSmPow(i_slice)=-MaxSmPow(i_slice)  ! for borderline cases SMPowMx is even negative!
         cycle  ! also filters out the border cases
      EndIf
      !write(2,*) 'd_Mx', i_slice, d_Mx(:), SMPowMx,QualMx(1:3)
      NMx=NMx+1
      If(.not. Dual) Then
         write(2,*) '!!!!!!!!!!!!!! Obsolete option, dual=.false. !!!!!!!!!!!!!'
      EndIf
         s=sqrt(MaxSmPowQ(i_slice)*MaxSmPowQ(i_slice)+MaxSmPowU(i_slice)*MaxSmPowU(i_slice))
         X=(ATAN2(MaxSmPowU(i_slice),MaxSmPowQ(i_slice))/2.+AngOff)*180./pi
         Y=t_ms-t_shft
      If(TRIDFile) Then
         If(RefinePolarizObs) Then
            !
            MaxSmPow(i_slice)=SMPowMx
            Call EI_PolarizSlice(i_slice)   ! Refine polarization analysis by calculating observables at interpolated position
            !
            !  For storing data for later processing in DataSelect:
            write(28,"(I5, ','F12.6,3(','F11.5),', ',F12.1, ',', F7.2,',', 1pg12.4, 6(' ,  ',g12.4','g12.4,' ' ))")  &
               i_slice, t_ms-t_shft, PixLoc(1:3)/1000., SMPowMx, Chi2pDF, StI,   Stk_NEh(1,1:3), Stk_NEh(2,2:3), Stk_NEh(3,3)
            !
            ! For making a plot of the just analyzed sources:
            write(29,"(I5,F13.6, 3F12.5, F13.1, F7.2, I6, g13.4, 3F6.1, ',' 3F6.1, 3(' , ',g12.4, 3(',',f7.2))  )") &
               i_slice, t_ms-t_shft, PixLoc(1:3)/1000.,  SMPowMx, &
               Chi2pDF, (1+i_slice*N_smth), StI, 100.*StI12/StI, 100*StI3/StI, dStI3/StI, &
               P_Un*100., P_Lin*100., P_Circ*100., (PolMag(k), PolZen(k), PolAzi(k), PoldOm(k), k=1,3)
        Else
            !  For storing data for later processing in whatever way:
            write(28,"(I5, ','F12.6,3(','F11.5),', ',F12.1, ',',  1pg12.4, 4(','g12.4) )")  &
               i_slice, t_ms-t_shft, PixLoc(1:3)/1000.,  SMPowMx, MaxSmPowI12(i_slice), &
               MaxSmPowQ(i_slice), MaxSmPowU(i_slice), MaxSmPowV(i_slice), MaxSmPowI3(i_slice)
            !
            ! For making a plot of the just analyzed sources:
            write(29,"(I7,F13.6, 3F12.5, F13.1  )") &
               (1+i_slice*N_smth), t_ms-t_shft, PixLoc(1:3)/1000.,  SMPowMx
         EndIf
         write(2,"(1x,A,i4,A,f12.6,A,2(F9.4,','),F9.4,A ,g13.6,A,3f6.2,A,5g13.3)") 'Slice#',i_slice,  &
            ', time=', t_ms-t_shft, '[ms], Max Intensity @ (N,E,h)=(',PixLoc(:)/1000.,') [km], I12_center=', &
            MaxSmPowI12(i_slice),', Q,U,V=',MaxSmPowQ(i_slice)*100., MaxSmPowU(i_slice)*100., MaxSmPowV(i_slice)*100., &
            '%, I3=', MaxSmPowI3(i_slice)
      EndIf
      If(.not. Dual) Then     ! obsolete
         write(2,"(A,I6, 2F7.3, F10.2, A, 4F10.5, ';',3F7.1, '%;',3F7.1)") &
                           'MxPow', i_slice, QualBar(3), QualBar(1), SMPowMx, ';', &
                           Y, PixLoc(1)/1000., PixLoc(2)/1000., PixLoc(3)/1000.
      EndIf
      flush(unit=2)
   EndDo
   If(TRIDFile) Then
      Close(Unit=29)
      If(Dual .and. RefinePolarizObs) Then
         !DeAllocate( CMTime_pix )
         Close(Unit=28)
      !   Close(Unit=27)
      !   Close(Unit=26)
      Else IF(.not. Dual) Then
         Close(Unit=28)
      EndIf
   EndIf
   If ( (NGridInterpol.gt.0) .and. FirstTimeInterf) Then
      write(2,"(A,I5,A,3F5.2,A,I4,A)") 'GridInterpol:', NGridInterpol, ', RMS of interpolation [grid spacing]=' &
         ,sqrt(RMSGridInterpol(:)),'; ', NGridExtrapol &
         ,' interpolations beyond 3-D limit of 0.75; ideal interpolation distance=sqrt(1/12)=0.29'
   EndIf
   If(FirstTimeInterf) Then
      Write(2,*) 'Source info written to file:', trim(DataFolder)//TRIM(OutFileLabel)//'TRID'//TRIM(txt)//'.csv'
      Write(2,*) 'Info for plots written to file:', trim(DataFolder)//TRIM(OutFileLabel)//'TRID'//TRIM(txt)//'.plt'
      write(2,"(A,2I6,2F6.3,A,3F7.3,A,I7,A,F9.1)") 'number of sources in plots:',NMx,NBar, NoiseLevel, RatMax
   EndIf
   If(NMx.gt.1 .and. TRIDFile) then
      Call GLEplotControl(PlotType='SrcsPltLoc', PlotName='TRID'//TRIM(OutFileLabel), &
         PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel)//'TRID'//TRIM(txt) )
      ! write(10,"('gle -d pdf -o ',A,'-InfImaMx_',I1,'.pdf ${UtilDir}/SourcesPlot.gle ${FlashFolder}/files/',A)") &!
      !   TRIM(OutFileLabel), i_eo, TRIM(OutFileLabel)//'IntfSpecPowMx'//TRIM(txt)
      If(RefinePolarizObs) Then
         Call GLEplotControl(PlotType='TRIDPol', PlotName='TRIDPol'//TRIM(OutFileLabel), &
            PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel) )
      !Else
      !   If(NBar.gt.1) &
      !      Call GLEplotControl(PlotType='SourcesPlot', PlotName='IntfBar'//TRIM(txt)//TRIM(OutFileLabel), &
      !         PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel)//'IntfSpecPowBar'//TRIM(txt) )
         ! write(10,"('gle -d pdf -o ',A,'-InfImgBar_',I1,'.pdf ${UtilDir}/SourcesPlot.gle ${FlashFolder}/files/',A)") &!
         !   TRIM(OutFileLabel), i_eo, TRIM(OutFileLabel)//'IntfSpecPowBar'//TRIM(txt)
         !--------------------------------------------------------------
      EndIf
   EndIf
   !
   Return
End Subroutine OutputIntfPowrMxPos
!-----------------------------------------------
Real*8 Function Sq(x)
Real*8 x
Sq=x
End Function Sq
!-----------------------------------------------
Subroutine FindInterpolMx(i,d_gr,SMPow,Qualty)
   ! Do a 3D parabolic fit to interpolate around the maximal pixel
   !--------------------------------------------
   use constants, only : dp
   Use Interferom_Pars, only : SumStrt, SumWindw, polar, N_pix, FirstTimeInterf, TRIDFile !, IntFer_ant
   Use Interferom_Pars, only : PixelPower, MaxSmPow, MaxSmPowGrd, PixSmPowTr, RMSGridInterpol, NGridInterpol, NGridExtrapol
   Implicit none
   Integer, intent(in) :: i
   Real(dp), intent(out) :: d_gr(1:3), SMPow, Qualty(1:3)
   integer :: j, i_gr(1:3)
   Integer, parameter :: jp(1:3)=(/2,3,1/), jm(1:3)=(/3,1,2/)
   Real(dp) :: y0, Ay(1:3), By(1:3), yp, ym, Det3, Ry(1:3), c_gr(1:3), Parb(-2:+2), Valu(-2:+2)
   Integer :: i_h, i_E, i_N, Nxx
   Real(dp) :: A, B
   Real(dp), parameter :: D0=0.75 ! maximal distance (in grid spacings) for the interpolated point from the max-intensity voxel
   Real(dp) :: Paraboloid ! , Sq
   !Real(dp), save :: RMSGridInterpol(1:3)=(/0.d0,0.d0,0.d0/)
   !Integer, save :: NGridInterpol=0
   Logical :: Check=.false.
   !Logical :: Check=.true.  !  .false.
   !
   !write(2,*) '!FindInterpolMx',i, MaxSmPowGrd(:,i)
   !flush(unit=2)
   i_gr(:)=MaxSmPowGrd(:,i)
   !write(2,*) '!FindInterpolMx; i_gr(:)',i, i_gr(:), MaxSmPow(0:i)
   flush(unit=2)
   y0=(PixSmPowTr(i,i_gr(1),i_gr(2),i_gr(3)))
   SMPow=-y0
   !write(2,*) 'y0',i,y0,i_gr(:)
   !flush(unit=2)
   !
   Do j=1,3  ! fit y=Ax^2/2+Bx+C; A=y"/d^2 & B=y'/2d & x_max= -B/A= -d y'/(2y")
      ! 3D case: y=Sum[A(i)X(i)^2/2 + R'(i,j)X(i)X(j)+B(i)x(i)+C]
      If((i_gr(j).eq.N_pix(j,1)) .or. (i_gr(j).eq.N_pix(j,2))) Then
         If((N_pix(1,2)*N_pix(2,2)*N_pix(3,2)).gt.1) Then  ! skip this for case of central pixel calculation for ATRID
            If(i.eq.0) Then
               Write(2,"(A,F9.2,A,3(I5,','))") &
                  ' Bordercase full window, Intensty=',y0,', Max @',MaxSmPowGrd(:,i)
               !write(2,*) "! BordercaseTest",j,i_gr(j),N_pix(j,1),N_pix(j,2)
            ElseIf(TRIDFile) then
               Write(2,"(A,i5,A,F9.2,A,3(I5,','))") &
                  ' Bordercase slice',i,', Intensty=',y0,', Max @',MaxSmPowGrd(:,i)
            EndIf
         EndIf
         d_gr(:)=i_gr(:) ! get rough position of the maximum on the border
         Return ! max pixel should not be at the edge of the hypercibe
      EndIf
      i_gr(j)=i_gr(j)+1
      yp=(PixSmPowTr(i,i_gr(1),i_gr(2),i_gr(3)) )
      i_gr(j)=i_gr(j)-2
      ym=(PixSmPowTr(i,i_gr(1),i_gr(2),i_gr(3)) )
      i_gr(j)=i_gr(j)+1    ! back to the initial value
      Ay(j)=(yp+ym-2*y0)  ! in units of 1/d^2(2*d_loc(j)*d_loc(j))  ! negative
      By(j)=(yp-ym)/2. ! in units of 1/d
      !d_gr(j)=-By/Ay ! in units of d
      !SMPow= SMPow -By*By/(2*Ay)
      !write(2,*) 'yp,ym',j,yp/y0, ym/y0, Ay(j)/y0,By(j)/y0
   Enddo
   !
   !B=(PixSmPowTr(i,i_gr(1)+1,i_gr(2),i_gr(3)+1)) - (PixSmPowTr(i,i_gr(1)-1,i_gr(2),i_gr(3)+1)) &
   !   +(PixSmPowTr(i,i_gr(1)+1,i_gr(2),i_gr(3)-1)) - (PixSmPowTr(i,i_gr(1)-1,i_gr(2),i_gr(3)-1)) &
   !   +(PixSmPowTr(i,i_gr(1)+1,i_gr(2)+1,i_gr(3))) - (PixSmPowTr(i,i_gr(1)-1,i_gr(2)+1,i_gr(3))) &
   !   +(PixSmPowTr(i,i_gr(1)+1,i_gr(2)-1,i_gr(3))) - (PixSmPowTr(i,i_gr(1)-1,i_gr(2)-1,i_gr(3)))
   !write(2,*) 'B(1)',B/8./y0,By(1)/y0
   !By(1)=(By(1)+B/4)/3.
   !B=(PixSmPowTr(i,i_gr(1),i_gr(2)+1,i_gr(3)+1)) - (PixSmPowTr(i,i_gr(1),i_gr(2)-1,i_gr(3)+1)) &
   !   +(PixSmPowTr(i,i_gr(1),i_gr(2)+1,i_gr(3)-1)) - (PixSmPowTr(i,i_gr(1),i_gr(2)-1,i_gr(3)-1)) &
   !   +(PixSmPowTr(i,i_gr(1)+1,i_gr(2)+1,i_gr(3))) - (PixSmPowTr(i,i_gr(1)+1,i_gr(2)-1,i_gr(3))) &
   !   +(PixSmPowTr(i,i_gr(1)-1,i_gr(2)+1,i_gr(3))) - (PixSmPowTr(i,i_gr(1)-1,i_gr(2)-1,i_gr(3)))
   !write(2,*) 'B(2)',B/8./y0,By(2)/y0
   !By(2)=(By(2)+B/4)/3.
   !B=(PixSmPowTr(i,i_gr(1),i_gr(2)+1,i_gr(3)+1)) - (PixSmPowTr(i,i_gr(1),i_gr(2)+1,i_gr(3)-1)) &
   !   +(PixSmPowTr(i,i_gr(1),i_gr(2)-1,i_gr(3)+1)) - (PixSmPowTr(i,i_gr(1),i_gr(2)-1,i_gr(3)-1)) &
   !   +(PixSmPowTr(i,i_gr(1)+1,i_gr(2),i_gr(3)+1)) - (PixSmPowTr(i,i_gr(1)+1,i_gr(2),i_gr(3)-1)) &
   !   +(PixSmPowTr(i,i_gr(1)-1,i_gr(2),i_gr(3)+1)) - (PixSmPowTr(i,i_gr(1)-1,i_gr(2),i_gr(3)-1))
   !write(2,*) 'B(3)',B/8./y0,By(3)/y0
   !By(3)=(By(3)+B/4)/3.
   !
   Ry(1)=((PixSmPowTr(i,i_gr(1),i_gr(2)+1,i_gr(3)+1)) + (PixSmPowTr(i,i_gr(1),i_gr(2)-1,i_gr(3)-1)) &
         -(PixSmPowTr(i,i_gr(1),i_gr(2)+1,i_gr(3)-1)) - (PixSmPowTr(i,i_gr(1),i_gr(2)-1,i_gr(3)+1)))/4.
   Ry(2)=((PixSmPowTr(i,i_gr(1)+1,i_gr(2),i_gr(3)+1)) + (PixSmPowTr(i,i_gr(1)-1,i_gr(2),i_gr(3)-1)) &
         -(PixSmPowTr(i,i_gr(1)+1,i_gr(2),i_gr(3)-1)) - (PixSmPowTr(i,i_gr(1)-1,i_gr(2),i_gr(3)+1)))/4.
   Ry(3)=((PixSmPowTr(i,i_gr(1)+1,i_gr(2)+1,i_gr(3))) + (PixSmPowTr(i,i_gr(1)-1,i_gr(2)-1,i_gr(3))) &
         -(PixSmPowTr(i,i_gr(1)+1,i_gr(2)-1,i_gr(3))) - (PixSmPowTr(i,i_gr(1)-1,i_gr(2)+1,i_gr(3))))/4.
   !
   If(Check) Then
      !  Checking the area around the max
      nxx=0
!      j=3
      Do j=1,3  ! fit y=Ax^2/2+Bx+C; A=y"/d^2 & B=y'/2d & x_max= -B/A= -d y'/(2y")
         If((i_gr(j).le.N_pix(j,1)+2) .or. (i_gr(j).ge.N_pix(j,2)-2)) nxx=nxx+1 ! max pixel should not be at the edge of the hypercibe
      Enddo
      !write(2,*) 'nxx:',nxx, i_gr(:), N_pix(:,1)+2,N_pix(:,2)-2
      If(nxx.eq. 0) then
         !write(2,*) 'max pixel=',i_gr(:),y0
         Do i_h=-2,2
            B=i_h
            c_gr(3)=B ! SIGN(SQRT(ABS(B)),B)
            Do i_N= -2,2
               !i_N=0
               B=i_N
               c_gr(1)=B !SIGN(SQRT(ABS(B)),B)
               Do i_E= -2,2
                  !i_E=0
                  B=i_E
                  c_gr(2)=B !SIGN(SQRT(ABS(B)),B)
                  !Valu(i_E)=1.-(PixSmPowTr(i,i_gr(1)+i_N,i_gr(2)+i_E,i_gr(3)+i_h))/y0
                  !Parb(i_E)=1.-Paraboloid(y0,Ay,By,Ry,c_gr)/y0
                  Valu(i_h)=1.-(PixSmPowTr(i,i_gr(1)+i_N,i_gr(2)+i_E,i_gr(3)+i_h))/y0
                  Parb(i_h)=1.-Paraboloid(y0,Ay,By,Ry,c_gr)/y0
               Enddo
               !write(2,"(A, 2i3,5F6.2,'; ',5F6.2)") 'area i_h,i_n',i_h,i_N, Valu(:), Parb(:)
               !flush(unit=2)
            Enddo
         Enddo
         !      write(2,"(/,A, 2i3,5F6.2,'; ',5F6.2,'; ',F9.1)") 'area i_h=-2:2,i_n=0',i_h,i_N, Valu(:), Parb(:),y0
      Endif
   EndIf
   !
   ! Solve for location maximum
   Det3=Ay(1)*Ay(2)*Ay(3)+2*Ry(1)*Ry(2)*Ry(3)-Ay(1)*Ry(1)*Ry(1)-Ay(2)*Ry(2)*Ry(2)-Ay(3)*Ry(3)*Ry(3) !-Ay()*R()*R()
   If(Det3.ge.0.) then  ! all eigenvalues should be negative for a real maximum (should always be fulfilled)
      B=Ay(1)*Ay(2)*Ay(3)
      !Write(2,*) i,'Det3', Det3, Ay(:), (4*Ry(j)*Ry(j)*Ay(j)/B ,j=1,3)
      !flush(unit=2)
      !write(2,*) 'B:',B, Ry(:),i_gr(:)
      !write(2,*) 'surroundings:',  PixSmPowTr(i,i_gr(1),i_gr(2)+1,i_gr(3)+1)-y0, PixSmPowTr(i,i_gr(1),i_gr(2)-1,i_gr(3)-1)-y0 &
      !   , PixSmPowTr(i,i_gr(1),i_gr(2)+1,i_gr(3)-1)-y0, PixSmPowTr(i,i_gr(1),i_gr(2)-1,i_gr(3)+1)-y0 &
      !   , PixSmPowTr(i,i_gr(1)+1,i_gr(2),i_gr(3)+1)-y0, PixSmPowTr(i,i_gr(1)-1,i_gr(2),i_gr(3)-1)-y0 &
      !   , PixSmPowTr(i,i_gr(1)+1,i_gr(2),i_gr(3)-1)-y0, PixSmPowTr(i,i_gr(1)-1,i_gr(2),i_gr(3)+1)-y0 &
      !   , PixSmPowTr(i,i_gr(1)+1,i_gr(2)+1,i_gr(3))-y0, PixSmPowTr(i,i_gr(1)-1,i_gr(2)-1,i_gr(3))-y0 &
      !   , PixSmPowTr(i,i_gr(1)+1,i_gr(2)-1,i_gr(3))-y0, PixSmPowTr(i,i_gr(1)-1,i_gr(2)+1,i_gr(3))-y0
      !
      !Return
   EndIf
   Do j=1,3 ! get relative position of the maximum
      d_gr(j)=-(By(j)*(Ay(jp(j))*Ay(jm(j))-Ry(j)*Ry(j)) - By(jp(j))*(Ay(jm(j))*Ry(jm(j))-Ry(jp(j))*Ry(j)) &
         - By(jm(j))*(Ay(jp(j))*Ry(jp(j))-Ry(jm(j))*Ry(j)) )/Det3
   Enddo
   !Do i_h=-1,1
   !   Do i_N=-1,1
   !      Do i_E=-1,1
   !         c_gr(1)=d_gr(1)+0.1*i_N; c_gr(2)=d_gr(2)+0.1*i_E; c_gr(3)=d_gr(3)+0.1*i_h
   !         If(Paraboloid(y0,Ay,By,Ry,c_gr).ge.SMPow) write(2,*) 'real max@',i_N,i_E,i_h
   !      Enddo
   !   Enddo
   !Enddo
   !write(2,"(3f7.3,' ; ',3F6.2)") Ry(:)/y0,d_gr(:)
   !If(i.eq.i_s) then
   !Write(2,"(A,F9.1,4(3F6.2,';'),F6.3,f9.4)") 'out:y0,Ay,By,Ry',y0,Ay/y0,By/y0,Ry/y0,d_gr,SMPow/y0,Det3/(y0**3)
   !EndIf
   !
   !write(2,*) 'y0:',y0,det3
   !write(2,*) 'FindInterpolMx:', i_gr(:), d_gr(:)
   RMSGridInterpol(:)=(RMSGridInterpol(:)*NGridInterpol+d_gr(:)*d_gr(:))/(NGridInterpol+1)
   NGridInterpol=NGridInterpol+1
   !write(2,*) 'GridInterpol:', NGridInterpol, RMSGridInterpol(:)
   B=sqrt(SUM(d_gr(:)*d_gr(:)))  ! Limit the max distance from the most intense voxel
   If(B.gt.D0) Then
      write(2,"(A,F5.2,A,3F5.2)") 'Interpolation distance shortend:',B, ', stepsizes were (in grid spacings)', d_gr(:)
      d_gr(:)=d_gr(:)*D0/B
      NGridExtrapol=NGridExtrapol+1
   EndIf
   SMPow=Paraboloid(y0,Ay,By,Ry,d_gr)  ! construct intensity at intepolated position of max
   !
   d_gr(:)=i_gr(:)+d_gr(:) ! get true position of the maximum
   Qualty(1:3)=Ay(1:3)/y0
   !
   Return
End Subroutine FindInterpolMx
!------------------------------
Function Paraboloid(y0,Ay,By,Ry,d_gr) result(SMPow)
! SMPow=Sum[A(i)X(i)^2/2 + R'(i,j)X(i)X(j)+B(i)x(i)+C]
   use constants, only : dp
   Implicit none
   Real(dp), intent(in) :: d_gr(1:3), y0, Ay(1:3), By(1:3), Ry(1:3)
   Real(dp) :: SMPow
   integer :: j
         SMPow=y0 + Ry(1)*d_gr(2)*d_gr(3) + Ry(2)*d_gr(3)*d_gr(1) + Ry(3)*d_gr(1)*d_gr(2)
         Do j=1,3
            SMPow=SMPow + Ay(j)*d_gr(j)*d_gr(j)/2. + By(j)*d_gr(j)
         Enddo
         !
End Function Paraboloid
!-----------------------------------------------
Subroutine FindBarycenter(i,Baryd,SMPow,Qualty)
   ! Use the barycenter to interpolate in distance and use only those pixels that have a large intensity
   ! Needs
   !     File 30 opened
   !--------------------------------------------
   use constants, only : dp, pi, Sample
   Use Interferom_Pars, only : SumStrt, SumWindw, polar, N_pix, d_loc, IntfLead, CenLocPol, CenLoc, RimSmPow
   Use Interferom_Pars, only : PixelPower, MaxSmPow !, N_smth
   Use Interferom_Pars, only : xMin, xMax, yMin, yMax, zMin, zMax, tMin, tMax
   Use Interferom_Pars, only : MaxSmPowGrd, PixSmPowTr, RatMax
   use Chunk_AntInfo, only : NoiseLevel
   Implicit none
   Integer, intent(in) :: i
   Real(dp), intent(out) :: Baryd(1:3), SMPow, Qualty(1:3)
   integer :: j, i_gr(1:3)
   Real(dp) :: y0, Ay(1:3), By(1:3), c_gr(1:3)
   Integer :: i_h, i_E, i_N, Volume, RimMaxPix(1:3)
   Real(dp) :: AveSmPow, A, B, dh, BaryH, BackGr
   Real(dp) :: RimMax, BackgrThresh, Epsln=epsilon(Epsln)  ! Only those pixels with power=y0*RatMax are used for barycentric analysis
   !
   i_gr(:)=MaxSmPowGrd(:,i)
   y0=PixSmPowTr(i,i_gr(1),i_gr(2),i_gr(3))
   SMPow=-1.
   !
   Qualty(1:3)=-1.
   If(y0.lt.NoiseLevel) Return
   BackgrThresh=y0*RatMax
9  continue
   AveSmPow=0.
   RimSmPow(:)=0.
   Volume=0
   Baryd(:)=0.
   RimMax= - BackgrThresh
   RimMaxPix(:)=0
   Do i_h=N_pix(3,1),N_pix(3,2)
      Do i_N=N_pix(1,1),N_pix(1,2)
         Do i_E=N_pix(2,1),N_pix(2,2)
            A = PixSmPowTr(i,i_N,i_E,i_h) - BackgrThresh
            If(A .gt. 0.) then ! consider pixels with large intensity only
               AveSmPow  = AveSmPow  + A  ! take care of only the part sticking out above 'background'
               Baryd(1) = Baryd(1) + i_N*A
               Baryd(2) = Baryd(2) + i_E*A
               Baryd(3) = Baryd(3) + i_h*A
               Volume=Volume+1
            EndIf
            If((i_h.eq.N_pix(3,1)) .or. (i_h.eq.N_pix(3,2)) .or. (i_E.eq.N_pix(2,1)) .or. (i_E.eq.N_pix(2,2)) .or. &
                   (i_N.eq.N_pix(1,1)) .or. (i_N.eq.N_pix(1,2))) Then ! rim of picture
               RimSmPow(i_h)=RimSmPow(i_h) + A
               If(A.gt.RimMax) then
                  RimMax=A
                  RimMaxPix(1)=i_N;  RimMaxPix(2)=i_E;  RimMaxPix(3)=i_h
               EndIf
            Endif
         Enddo ! i_E
      Enddo ! i_N
   EndDo  ! i_h
   If(RimMax.gt.0.) then  ! "backgound" sticking out to the rim and should be subtracted (in a way, background level is raised)
      !Write(2,*) 'RimMax:',i, RimMax,y0
      BackgrThresh=(BackgrThresh + RimMax)*(1.+ Epsln)
      Qualty(1) = BackgrThresh/y0
      If(BackgrThresh/y0.gt.0.95) Return
      goto 9
   EndIf
   Baryd(:)=Baryd(:)/AveSmPow  ! now independent of Volume barycenter of max in pixel units
   !RimSmPow(:)=RimSmPow(:)*0.5/(N_pix(1,2)-N_pix(1,1) +N_pix(2,2)-N_pix(2,1)) ! edges accounted for by +1's in sum-range
   SMPow = y0
   !BackGr=BackgrThresh + MINVAL(RimSmPow(N_pix(3,1):N_pix(3,2)))
   !write(2,*) BackgrThresh/y0, RimMax/y0
8  Continue
   Qualty(1) = (BackgrThresh + RimMax)/y0
   Qualty(2) = Volume*1.
   Qualty(3) = (BackgrThresh + AveSmPow/Volume)/y0
   !write(2,"(A,I6, 2F7.3, F8.2,I7,A,10F7.2)") 'AveSmPow', i_s, AveSmPow/y0, BackGr/y0, y0, volume, ';', RimSmPow(:)/y0
   !write(2,"(A,2x,3F7.1,';',3F7.1)") 'Baryd(:)',
   !
   Return
End Subroutine FindBarycenter
!-----------------------------------------------
Subroutine PlotMultContours()
   use constants, only : dp
   use DataConstants, only : DataFolder, OutFileLabel
   Use Interferom_Pars, only : NrPixSmPowTr !, SumWindw  !, polar, N_pix, d_loc, IntFer_ant, IntfLead, CenLocPol, CenLoc
   !Use Interferom_Pars, only : AveInten, AveIntenE, AveIntenN, RimInten, MaxIntfInten, MaxIntfIntenLoc, FirstTimeInterf
   !Use Interferom_Pars, only : i_chunk, SlcInten, NrSlices, SliceLen, MaxSlcInten, MaxSlcIntenLoc   !IntfBase, IntfDim, IntfPhaseCheck, SumStrt, SumWindw
   use GLEplots, only : GLEplotControl
   Implicit none
   character(len=2) :: txt1
   Integer :: i_Slice
   Do i_slice=1,NrPixSmPowTr ! N_smth+1,SumWindw-N_smth,N_smth
      write(txt1,"(I2.2)") i_Slice
      Call  WriteBeamingIntensities(i_Slice, txt1)
   EndDo
   Call GLEplotControl(PlotType='MultContour', PlotName='MultCont'//TRIM(OutFileLabel), &
            PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel))
   !Call GLEplotControl(CleanPlotFile=TRIM(OutFileLabel)//'IntfSpecWin'//'*.csv')
   !Call GLEplotControl(CleanPlotFile=TRIM(OutFileLabel)//'_EISpec'//'*.dat')
   !Call GLEplotControl(CleanPlotFile=TRIM(OutFileLabel)//'Interferometer'//'*.z', Submit=.true.)
End Subroutine PlotMultContours
!-----------------------------------------------
Subroutine WriteBeamingIntensities(i_Slice, txt1)
   use constants, only : dp
   use DataConstants, only : DataFolder, OutFileLabel
   Use Interferom_Pars, only : N_pix, PixSmPowTr, MaxSmPow
   Implicit none
   Integer, intent(in) :: i_Slice
   character(len=2), intent(in) :: txt1
   integer :: i_N, i_E, i_h
   character(len=4) :: txt
   Character(len=20) :: FMT
   Real :: SlMax
   !
   txt=txt1//'EN'
   SlMax=1.
   If(i_Slice.ne.0) Then
      SlMax=abs(MaxSmPow(i_slice))/100.
   EndIf
   i_h=0
   write(FMT,"(A,I3,A)") '(',N_pix(2,2)-N_pix(2,1)+1,'(g12.4,1x) )'
   OPEN(UNIT=30,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'Interferometer'//txt//'.z')
   Write(30,"(10(A,I4))") '! nx ', N_pix(2,2)-N_pix(2,1)+1, ' ny ', N_pix(1,2)-N_pix(1,1)+1,  &
      ' xmin ', N_pix(2,1), ' xmax ', N_pix(2,2) , ' ymin ', N_pix(1,1) , ' ymax ', N_pix(1,2)
   !If(FirstTimeInterf) write(2,*) 'OutputIntfPowrTotal:', trim(DataFolder)//TRIM(OutFileLabel)//'Interferometer'//txt//'.z'
   Do i_N=N_pix(1,1),N_pix(1,2)   ! or Phi  ! Loop over pixels
      write(30,FMT) PixSmPowTr(i_Slice, i_N, N_pix(2,1):N_pix(2,2), i_h)/SlMax
   Enddo
   Close(UNIT=30)
   i_N=0
   txt=txt1//'Eh'
   !write(FMT,"(A,I3,A)") '(',N_pix(2,2)-N_pix(2,1)+1,'(g12.4,1x) )'
   OPEN(UNIT=30,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'Interferometer'//txt//'.z')
   Write(30,"(10(A,I4))") '! nx ', N_pix(2,2)-N_pix(2,1)+1, ' ny ', N_pix(3,2)-N_pix(3,1)+1,  &
      ' xmin ', N_pix(2,1), ' xmax ', N_pix(2,2) , ' ymin ', N_pix(3,1) , ' ymax ', N_pix(3,2)
   Do i_h=N_pix(3,1),N_pix(3,2)   ! or elevation angle
      write(30,FMT) PixSmPowTr(i_Slice, i_N, N_pix(2,1):N_pix(2,2), i_h)/SlMax
   Enddo
   Close(UNIT=30)
   i_E=0
   txt=txt1//'hN'
   write(FMT,"(A,I3,A)") '(',N_pix(3,2)-N_pix(3,1)+1,'(g12.4,1x) )'
   OPEN(UNIT=30,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'Interferometer'//txt//'.z')
   Write(30,"(10(A,I4))") '! nx ', N_pix(3,2)-N_pix(3,1)+1, ' ny ', N_pix(1,2)-N_pix(1,1)+1,  &
      ' xmin ', N_pix(3,1), ' xmax ', N_pix(3,2) , ' ymin ', N_pix(1,1) , ' ymax ', N_pix(1,2)
   Do i_N=N_pix(1,1),N_pix(1,2)   ! or Phi  ! Loop over pixels
      write(30,FMT) PixSmPowTr(i_Slice, i_N, i_E, N_pix(3,1):N_pix(3,2))/SlMax
   Enddo
   Close(UNIT=30)
End Subroutine WriteBeamingIntensities
!-----------------------------------------------
Subroutine OutputIntfSlices(i_eo)
   ! Analyze the Slices
   ! Needs
   !--------------------------------------------
   use constants, only : dp, pi, Sample
   use DataConstants, only : DataFolder, OutFileLabel
   use Chunk_AntInfo, only : CTime_spectr
   use ThisSource, only : CurtainHalfWidth, Dual, SourcePos
   Use Interferom_Pars, only : CTime_sum, SumStrt, SumWindw, polar, N_pix, d_loc, IntFer_ant, IntfLead, CenLocPol, CenLoc
   Use Interferom_Pars, only : AveInten, AveIntenE, AveIntenN, RimInten, MaxIntfInten, MaxIntfIntenLoc, FirstTimeInterf
   Use Interferom_Pars, only : i_chunk, SlcInten, NrSlices, SliceLen, MaxSlcInten, MaxSlcIntenLoc   !IntfBase, IntfDim, IntfPhaseCheck, SumStrt, SumWindw
   Implicit none
   Integer, intent(in) :: i_eo
   integer :: i, j, i_s  !_loc(1), i_ant, j_IntFer, MLoc, Mloc_all
   character(len=5) :: txt
   Real(dp) :: Rdist, PixLocPol(1:3), PixLoc(1:3)
   !
   ! write windowed spectra to file
   Write(txt,"('_',i1)") i_eo
   If(Dual) txt='_d'
!   OPEN(UNIT=29,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'IntfSpecWin'//TRIM(txt)//'.csv')
!   If(polar) Then
!      write(29,"(g12.4,5(',',i4),',',i8,',',i5,3(',',f9.4))") MaxIntfInten, N_pix(1,2), N_pix(2,1), N_pix(2,2), &
!         N_pix(3,1), N_pix(3,2), SumStrt,SumWindw, d_loc(1:2)*180./pi ,d_loc(3)
!   Else
!      write(29,"(g12.4,5(',',i4),',',i8,',',i5,3(',',f9.4))") MaxIntfInten, N_pix(1,2), N_pix(2,1), N_pix(2,2), &
!         N_pix(3,1), N_pix(3,2), SumStrt,SumWindw, d_loc(1:3)
!   EndIf
!   Do i=0, SumWindw
!      j=i+SumStrt
!      write(29,"(i6,',',g12.4,',',g12.4)") i, (ABS(CTime_spectr(j,IntFer_ant(1,i_chunk),i_chunk)) )**2, ABS(CTime_sum(j))**2
!   Enddo
!   Close(Unit=29)
   !
   !--------------------------------------------------------------
   ! Write general info for this picture
   If(FirstTimeInterf) write(2,"(A,G11.3,A,3f9.4,A)", ADVANCE='NO') 'Maximum ',MaxIntfInten, &
         ' @ (N,E,h)=(',MaxIntfIntenLoc(:)/1000.,') [km]'
   Call Carth2Pol(MaxIntfIntenLoc(:),PixLocPol)
   If(CurtainHalfWidth.gt.0) SourcePos(:,i_eo+1)=MaxIntfIntenLoc(:)  ! for curtainplot but interferes with other usages of SourcesPos
   If(FirstTimeInterf) write(2,"(A,f9.4,f8.2,f7.2)") ' = (ph,th,R)=',PixLocPol(1)*180/pi,PixLocPol(2)*180./pi,PixLocPol(3)/1000.
   ! Info per slice
   AveIntenN(:,:)=AveIntenN(:,:)*d_loc(1)/AveInten(:,:)
   AveIntenE(:,:)=AveIntenE(:,:)*d_loc(2)/AveInten(:,:)
   RimInten(:,:)=RimInten(:,:)*0.5/(N_pix(1,2)-N_pix(1,1) +N_pix(2,2)-N_pix(2,1)) ! edges accounted for by +1's in sum-range
   AveInten(:,:)=AveInten(:,:)/((N_pix(1,2)-N_pix(1,1)+1)*(N_pix(2,2)-N_pix(2,1)+1))
   If(NrSlices.gt.1) Then  !  this is usually =1; obsolete thus
      OPEN(UNIT=29,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'IntfTrack'//TRIM(txt)//'.csv')
      Write(29,"(2i8,3f9.4,f9.4,f7.2,f8.2,2f8.4)") SliceLen,IntfLead,MaxIntfIntenLoc(:)/1000., &
         PixLocPol(3)/1000.,PixLocPol(2)*180./pi,PixLocPol(1)*180/pi, d_loc(1:2)*180./pi
      Write(29,"(A,8x,A,11x,A,21x,A)") '! Layer; MaxIntensity;','(N,E,h)|Max;','(R,th,ph)|Max;','(R,th,ph)|Ave'
   EndIf
   Do i_s=1,NrSlices
      write(2,"(A,i4,A,G11.3,A,3f9.4,A)", ADVANCE='NO') 'slice number=',i_s,', max intensity=',MaxSlcInten(i_s),&
         ' @ (N,E,h)=',MaxSlcIntenLoc(i_s,:)/1000.,'[km]'
      If(polar) then
         Call Carth2Pol(MaxSlcIntenLoc(i_s,:),PixLocPol)
         write(2,"(A,f9.4,f8.2,f7.2)", ADVANCE='NO') ' = (ph,th,R)=',PixLocPol(1)*180./pi,PixLocPol(2)*180/pi,PixLocPol(3)/1000.
         If(NrSlices.gt.1) write(29,"(i5,',',G11.3,',',3(f9.4,','),4x,f9.4,',',f9.4,',',f10.4)", ADVANCE='NO') i_s, &
            MaxSlcInten(i_s), MaxSlcIntenLoc(i_s,:)/1000., PixLocPol(3)/1000., PixLocPol(2)*180./pi, PixLocPol(1)*180/pi
         If((AveInten(i_s,0)/RimInten(i_s,0)).gt.2.) then ! rather arbitrary criterion to determine if a peak is significant
            PixLocPol(1)=CenLocPol(1) + AveIntenN(i_s,0)*AveInten(i_s,0)/(AveInten(i_s,0)-RimInten(i_s,0))
            PixLocPol(2)=CenLocPol(2) + AveIntenE(i_s,0)*AveInten(i_s,0)/(AveInten(i_s,0)-RimInten(i_s,0))
            write(2,"(A,f7.2,f8.2)") ', Barycenter @ (th,ph)=', PixLocPol(2)*180./pi,PixLocPol(1)*180/pi
         Else
            PixLocPol(1)=CenLocPol(1) !+ AveIntenN(i_s,0)*AveInten(i_s,0)/(AveInten(i_s,0)-RimInten(i_s,0))
            PixLocPol(2)=CenLocPol(2) !+ AveIntenE(i_s,0)*AveInten(i_s,0)/(AveInten(i_s,0)-RimInten(i_s,0))
            write(2,*) ', Barycenter not calculated, (AveInten(i_s,0)/RimInten(i_s,0))=', (AveInten(i_s,0)/RimInten(i_s,0))
         EndIf
         PixLocPol(3)=CenLocPol(3)
         If(NrSlices.gt.1) write(29,"( ',',4x,f9.4,',',f9.4,',',f10.4)") &
            PixLocPol(3)/1000., PixLocPol(2)*180./pi, PixLocPol(1)*180/pi
      else
         !write(2,*) ''
         PixLoc(1)=CenLoc(1)+AveIntenN(i_s,0)
         PixLoc(2)=CenLoc(2)+AveIntenE(i_s,0)
         PixLoc(3)=CenLoc(3)
         write(2,"(A,3f9.4)") 'Barycenter intensity @ (N,E,h)=',PixLoc(1),PixLoc(2)
         If(NrSlices.gt.1) write(29,"( i5,',',G11.3,',',3(f9.4,',') )") i_s,MaxSlcInten(i_s),MaxSlcIntenLoc(i_s,:)/1000.
      Endif
      !Do i_h=-N_h,N_h
      !i_h=0
         !write(2,*) 'height nr:',i_h, sqrt(RimInten(i_s,i_h)), sqrt(AveInten(i_s,i_h)) !,AveIntenE(i_s,i_h),AveIntenN(i_s,i_h)
         !write(2,*) 'height nr:',i_h, RimInten(i_s,i_h)/AveInten(i_s,i_h), SQRT(AveInten(i_s,i_h)) !,AveIntenE(i_s,i_h),AveIntenN(i_s,i_h)
      !Enddo
   Enddo
   If(NrSlices.gt.1) Close(Unit=29)
   !
   Return
End Subroutine OutputIntfSlices
!-----------------------------------------------
Subroutine EI_PolarizSlice(i_slice)
   use constants, only : dp
   use DataConstants, only : Ant_nrMax
   use Interferom_Pars, only :  Nr_IntFerMx, Nr_IntferCh ! the latter gives # per chunk
   use Interferom_Pars, only : Cnu_p0, Cnu_t0, Cnu_p1, Cnu_t1, IntfNuDim, AntPeak_OffSt
   Use Interferom_Pars, only : IntfBase, IntfLead, SumWindw, N_Smth, FirstTimeInterf
   Use Interferom_Pars, only : i_chunk, PixLoc, CenLoc
   Implicit none
   Integer, intent(in) :: i_slice
   Integer :: i_sample, i_1, i_2, Output, j_IntFer, ExclStat(1:30)= (/ (0.d0, I_1 = 1, 30) /) !
   Real(dp) :: ChiSq, DelChi(-N_Smth:+N_Smth,1:Nr_IntFerMx)
   Real(dp) :: FitDelay(1:Ant_nrMax), VoxLoc(1:3), del_1, del_2
   Character(len=8) :: Label
   logical :: First=.true.
   If(First) Then
      First=.false.
      Write(2,*) ' IntfBase, IntfNuDim:', IntfBase, IntfNuDim, Nr_IntFerMx
      Call EI_PolSetUp(Nr_IntFerCh(i_chunk), IntfBase, i_chunk, CenLoc(:), &
         AntPeak_OffSt(1,1), Cnu_p0(0,1,1), Cnu_t0(0,1,1), Cnu_p1(0,1,1), Cnu_t1(0,1,1))
      j_IntFer=1
      write(Label,"(i2.2)") j_IntFer
   EndIf
   !
   Output=2
   If(SumWindw/N_Smth.ge.11) Output=1
   i_sample=IntfLead + 1+i_slice*N_smth
   FitDelay(:)=0.
   write(Label,"('Slc ',i4.2)") i_slice
   write(2,"(1x,A,i4,A,I5,A,2(F9.4,','),F9.4,A)") 'Slice',i_slice,', Ref. ant. sample=',IntfBase+i_sample, &
      ', Max Intensity @ (N,E,h)=(',PixLoc(:)/1000.,') [km]'
   !write(2,*) 'EI_PolarizSlice: i_sample, i_chunk:', i_sample, i_chunk, Output, Nr_IntFerCh(i_chunk) ! ,'fitdelay:',FitDelay
   Call EI_PolGridDel(Nr_IntFerCh(i_chunk), FitDelay, i_sample, i_chunk, PixLoc(:), AntPeak_OffSt(1,1), &
      Cnu_p0(0,1,1), Cnu_t0(0,1,1), Cnu_p1(0,1,1), Cnu_t1(0,1,1), Output, DelChi, Label, ExclStat)
   !
   Return
   If(SumWindw/N_Smth.ge.5) Return
   !i_sample=IntfLead + i_slice*N_smth  ! approximately correct, upto rounding errors for dt
   Del_1=2  ! in meter
   Del_2=2
   Output=1
   write(2,*) i_slice, i_sample,'Del_1=',Del_1,', Del_2=',Del_2, ' [m]'
   Do i_1=-2,2
      Do i_2=-2,2
         VoxLoc(1)=PixLoc(1) + i_1*del_1
         VoxLoc(2)=PixLoc(2) + i_2*del_2
         VoxLoc(3)=PixLoc(3)
         write(Label,"('Gr ',I2,',',i2)") i_1,i_2
         Call EI_PolGridDel(Nr_IntFerCh(i_chunk), FitDelay, i_sample, i_chunk, VoxLoc(:), AntPeak_OffSt(1,1), &
            Cnu_p0(0,1,1), Cnu_t0(0,1,1), Cnu_p1(0,1,1), Cnu_t1(0,1,1), Output, DelChi, Label, ExclStat)
      Enddo
   Enddo
   !
   Return
End Subroutine EI_PolarizSlice
!-----------------------------------------------
!-----------------------------------------------

 Subroutine matinv3(A,B)
   use constants, only : dp
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    !complex(wp), intent(in) :: A(3,3)   !! Matrix
    !complex(wp)             :: B(3,3)   !! Inverse matrix
    !complex(wp)             :: detinv
    Real(dp), intent(in) :: A(3,3)   !! Matrix
    Real(dp)             :: B(3,3)   !! Inverse matrix
    Real(dp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    Return
End Subroutine matinv3
!----------------------------------------------
