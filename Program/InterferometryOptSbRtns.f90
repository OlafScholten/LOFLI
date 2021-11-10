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
Subroutine PixBoundingBox()
   ! Check max time shift for neighboring pixels; should ideally differ by 1 sample for neighboring pixels
   ! Calculate approximate bounding box of pixel image
   ! Setup arrays
   !  Needs
   !     call SelectIntfAntennas  first
   ! ------------------------------------
   use constants, only : dp !, pi, Sample, Refrac, c_mps
   use DataConstants, only : Time_Dim
   use Chunk_AntInfo, only : Ant_pos, Ant_RawSourceDist  ! CTime_spectr, Ant_Stations, Ant_IDs, Ant_nr,
   Use Interferom_Pars, only : i_chunk, IntFer_ant, Nr_IntFer
   Use Interferom_Pars, only : CenLoc, CenLocPol, Polar, d_loc, Diff, N_pix
   Use Interferom_Pars, only : xMin, xMax, yMin, yMax, zMin, zMax
   Use Interferom_Pars, only : SumStrt, SumWindw, IntfDim, IntfNuDim, IntfLead, IntfBase
   Implicit none
   !Integer, intent(in) :: i_chunk
   Real(dp) :: PixLoc(1:3), PixLocPol(1:3), t_shft, RDist
   Real(dp) :: t_n, t_p, mean, RMS, Err(3)
   Integer :: i, i_ant, j_IntFer
   !
   PixLoc(:)=CenLoc(:)
   xMin=+99999 ; xMax=-99999; yMin=+99999; yMax=-99999; zMin=+99999; zMax=-99999 !; tMin, tMax,
   Err(:)=0.
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
      Mean=0. ; RMS=0.
      Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
         i_ant=IntFer_ant(j_IntFer)
         Call RelDist(PixLoc(1),Ant_pos(1,i_ant,i_chunk),RDist)
         t_shft=Rdist - Ant_RawSourceDist(i_ant,i_chunk)
         Mean=Mean+t_shft
         RMS=RMS + t_shft*t_shft
         if(t_shft.lt.t_n) t_n=t_shft
         if(t_shft.gt.t_p) t_p=t_shft
         !write(2,*) i_ant,'t_shft',t_shft, Rdist, Ant_RawSourceDist(i_ant,i_chunk)
      EndDo ! j_IntFer
      Mean=Mean/Nr_IntFer
      RMS=RMS/Nr_IntFer
      RMS=sqrt(RMS-Mean*Mean)
      write(2,"(A,i2,2f7.2,A,f6.2,A,f9.5,A,3f7.2,A)") 'min & max time shift [samples] per pixel',i, t_n, t_p,', RMS=', RMS, &
         ', for d(i)=',d_loc(i),', d(N,E,h)=',Pixloc(:)-CenLoc(:),'[m]'
      diff(i) =t_p
      Err(:)=Err(:) + (Pixloc(:)-CenLoc(:))**2/(RMS*RMS)
      xMin=min(xmin,CenLoc(2)+N_pix(i,1)*abs(CenLoc(2)-PixLoc(2))); xMax=max(xMax,CenLoc(2)+N_pix(i,2)*abs(CenLoc(2)-PixLoc(2)))
      yMin=min(ymin,CenLoc(1)+N_pix(i,1)*abs(CenLoc(1)-PixLoc(1))); yMax=max(yMax,CenLoc(1)+N_pix(i,2)*abs(CenLoc(1)-PixLoc(1)))
      zMin=min(zmin,CenLoc(3)+N_pix(i,1)*abs(CenLoc(3)-PixLoc(3))); zMax=max(zMax,CenLoc(3)+N_pix(i,2)*abs(CenLoc(3)-PixLoc(3)))
      If(-t_n.gt.t_p) diff(i) =-t_n
      PixLoc(i)=CenLoc(i)  ! set back to center for Cartesian coordinates
   Enddo ! i
   !write(2,*) 'Error ellips (N,E,h)=',sqrt(Err(:)),' [m]'
   If(zMin.lt.0.) zMin=0.
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
   !write(2,*) 'diagnostics;(lead,dim):', IntfLead, IntfDim
   !write(*,"(A,o10,b20)") '(lead,dim):', IntfLead, IntfDim
   If((IntfBase+IntfDim .gt. Time_dim) .or. (IntfBase .lt. 1)) then
      write(2,*) 'dimensions for interferometry out of range: ', IntfBase, IntfBase+IntfDim, Time_dim
      write(2,*) 'input for interferometry: ', SumStrt, SumWindw
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
   use constants, only : dp, pi, Sample  ! , Refrac, c_mps
   use DataConstants, only : DataFolder, OutFileLabel, Time_Dim
   use Chunk_AntInfo, only : TimeBase
   use Chunk_AntInfo, only : CTime_spectr
   Use Interferom_Pars, only : CTime_sum, SumStrt, SumWindw, polar, N_pix, d_loc, IntFer_ant, IntfLead, CenLocPol, CenLoc, RimSmPow
   Use Interferom_Pars, only : PixelPower, MaxSmPow, N_smth, StartTime_ms
   Use FitParams, only : AntennaRange
   Use ThisSource, only : Dual
   !Use Interferom_Pars, only : xMin, xMax, yMin, yMax, zMin, zMax, tMin, tMax
   Use Interferom_Pars, only : i_chunk, NrPixSmPowTr, MaxSmPowGrd, PixSmPowTr, t_shft, NewCenLoc
   !Use Interferom_Pars, only : SlcInten, NrSlices, SliceLen, MaxSlcInten, MaxSlcIntenLoc   !IntfBase, IntfDim, IntfPhaseCheck, SumStrt, SumWindw
   Implicit none
   Integer, intent(in) :: RefAntSAI, i_eo
   integer :: i, j, i_N, i_E, i_h, SSm=5 ! Resum Hilbert envelope
   character(len=4) :: txt
   character(len=1) :: txt1
   Real(dp) :: SMPowMx, d_Mx(1:3), QualMx(10), SMPowBar, d_Bar(1:3), QualBar(10)
   Real(dp) :: Pref, Psum, half
   Real(dp) :: PixLocPol(1:3), PixLoc(1:3)
   !
   !  write intensity distribution for pixel-planes with i_h=-1,0,+1 needed for contour plots
   If(Dual) then
      txt1='d'
   Else
      Write(txt1,"(i1)") i_eo
   Endif
   i_h=0
   !write(2,*) 'started:',txt
   !flush(unit=2)
   txt=txt1//'_NE'
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
   !write(2,*) 'Done:',txt
   !flush(unit=2)
   i_N=0
   txt=txt1//'_Eh'
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
   !write(2,*) 'Done:',txt
   !flush(unit=2)
   i_E=0
   txt=txt1//'_Nh'
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
   !write(2,*) 'Done:',txt
   !flush(unit=2)
   !
   !---------------
   ! Do maximum pixel analysis for complete trace
   i=0
   Call FindInterpolMx(i,d_Mx, SMPowMx,QualMx)
   !
   If(SMPowMx.gt.0.) then  ! max not at rim
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
   ! =======================================================================
   !--------------------------------------------
   ! output powerspectrum, resummed
   !SSm=5  ! Resum Hilbert envelope
   i=0   ;  If(polar) i=1 ! 0=false; 1=true
   half=(SSm+1.)/2.
   If(Dual) then
      OPEN(UNIT=30,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'_EISpeceo.csv')
      Write(30,"(f9.3,3(',',f8.3),',',f5.1,',',i2,',',i2,F8.3,F9.2,' 0')") StartTime_ms, CenLoc/1000., AntennaRange, SSm, i, &
         t_shft*1000., SMPowMx
      Do i=0,Time_Dim/SSm
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
      Do i=0,Time_Dim/SSm
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
   use constants, only : dp, pi, Sample, Refrac, c_mps
   use DataConstants, only : DataFolder, OutFileLabel
   use Chunk_AntInfo, only : TimeBase
   Use Interferom_Pars, only : SumStrt, SumWindw, polar, N_pix, d_loc, IntFer_ant, IntfLead, CenLocPol, CenLoc, RimSmPow
   Use Interferom_Pars, only : PixLoc, PixelPower, MaxSmPow, N_smth, Nr_IntFer, IntfNuDim
   Use Interferom_Pars, only : MaxSmPowQ, MaxSmPowU, MaxSmPowV, MaxSmPowI3, MaxSmPowU1, MaxSmPowV1, MaxSmPowU2, MaxSmPowV2
   Use Interferom_Pars, only : alpha, PolBasis, PowerScale
   Use Interferom_Pars, only : xMin, xMax, yMin, yMax, zMin, zMax, tMin, tMax, AmpltPlot, t_offsetPow
   Use Interferom_Pars, only : NrPixSmPowTr, MaxSmPowGrd, PixSmPowTr, RatMax, PixPowOpt
   Use Interferom_Pars, only : StI, StI12, StQ, StU, StV, StI3, StU1, StV1, StU2, StV2
   Use Interferom_Pars, only : dStI, dStI12, dStQ, dStU, dStV, dStI3, dStU1, dStV1, dStU2, dStV2, Chi2pDF
   use Chunk_AntInfo, only : NoiseLevel
   use ThisSource, only : Dual
   use GLEplots, only : GLEplotControl
   !Use Interferom_Pars, only : SlcInten, NrSlices, SliceLen, MaxSlcInten, MaxSlcIntenLoc   !IntfBase, IntfDim, IntfPhaseCheck, SumStrt, SumWindw
   Implicit none
   Integer, intent(in) :: i_eo
   integer :: i_slice, i_1, i_2, i_3, i_gr(1:3)
   character(len=5) :: txt
   Real(dp) :: PixLocPol(1:3), t_ms, t_shft, y0
   Integer :: NBar, NMx
   Real(dp) :: SMPowMx, d_Mx(1:3), QualMx(10), SMPowBar, d_Bar(1:3), QualBar(10)
   Real(dp) :: s, X, Y, AngOff
   Real(dp) :: HorVec(1:3), VerVec(1:3)
   Complex(dp), allocatable :: CMTime_Pix(:,:)
   !
   Write(txt,"('_',i1)") i_eo
   If(Dual) Then
      txt='_d'
      Allocate( CMTime_pix(1:2*IntfNuDim,1:3) )
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
   write(2,*) xMin/1000.-.05, xMax/1000.+.05, yMin/1000.-.05, yMax/1000.+.05, zMin/1000.-.05, zMax/1000.+.05
   !
   If(Dual) Then
      OPEN(UNIT=28,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'IntfSpecPolrz'//TRIM(txt)//'.dat')
      write(28,"(2F11.3,4F10.2,1x,A,I3,F7.1,' 0')") tMin-.0005-TimeBase, tMax+.0005-TimeBase, &
         TimeBase,  CenLoc(:),TRIM(OutFileLabel)//TRIM(txt), PixPowOpt , 0.0        ! gleRead:  xMin xMax yMin yMax zMin zMax tMin tMax ZmBx$ t_offst
      Write(28,"(A,9x,A,13x,';   ', 9(A,'/I [%]   ;   '),A )") &
         '!   nr,  time[ms] ;','StI','StI12', '  StQ','  StU','  StV',' StI3',' StU1',' StV1',' StU2',' StV2', 'Chi^2/DoF'
   Else
      OPEN(UNIT=28,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'IntfSpecPowBar'//TRIM(txt)//'.dat')
      write(28,"(6F8.2,2F9.3,A,F7.1,' 0')") xMin/1000.-.005, xMax/1000.+.005, yMin/1000.-.005, yMax/1000.+.005, &
         zMin/1000.-.05, zMax/1000.+.05, tMin-.0005-TimeBase, tMax+.0005-TimeBase, ' IntfBox ', TimeBase         ! gleRead:  xMin xMax yMin yMax zMin zMax tMin tMax ZmBx$ t_offst
      Write(28,"(A,I7,A,A,F8.3,A,F6.3)") &
         '0 ',NrPixSmPowTr,' 0 0 ',TRIM(OutFileLabel)//TRIM(txt),AmpltPlot,' 0 0 0 0 0 0 0 ',NoiseLevel ! gleRead:  NTracks EventNr Q t_start label$ AmplitudePlot a1 b1 c1 a2 b2 c2 d2
   EndIf
   OPEN(UNIT=29,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'IntfSpecPowMx'//TRIM(txt)//'.dat')
   write(29,"(6F8.2,2F9.3,A,F7.1,i3,' 0')") xMin/1000.-.005, xMax/1000.+.005, yMin/1000.-.005, yMax/1000.+.005, &
      zMin/1000.-.05, zMax/1000.+.05, tMin-.0005-TimeBase, tMax+.0005-TimeBase, ' IntfBox ', TimeBase, PixPowOpt ! gleRead:  xMin xMax yMin yMax zMin zMax tMin tMax ZmBx$ t_offst
   Write(29,*) '0 ',NrPixSmPowTr,' 0 0 ',TRIM(OutFileLabel)//TRIM(txt),AmpltPlot,' 0 0 0 0 0 0 0 ',NoiseLevel, PowerScale, ' !' ! gleRead:  NTracks EventNr Q t_start label$ AmplitudePlot a1 b1 c1 a2 b2 c2 d2
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
            If(.not.Dual) Then  ! write frame
               Write(28,"(1x,3F9.3,2x,3I3,' ;')") PixLoc(2)/1000.,PixLoc(1)/1000.,PixLoc(3)/1000.,i_1,i_2,i_3 ! (E,N,h)=(x,y,z)
            EndIf
            Write(29,"(1x,3F9.3,2x,3I3,' ;')") PixLoc(2)/1000.,PixLoc(1)/1000.,PixLoc(3)/1000.,i_1,i_2,i_3 ! (E,N,h)=(x,y,z)
         Enddo
      Enddo
   Enddo
   !
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
   Write(2,"(A,F8.3,A)") 'overlap radial distance and longitudinal pol=',X*100.,'%'
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
   Write(2,*) 'CenLoc(:)',CenLoc(:)
   Write(2,*) 'HorVec(:)',HorVec(:)
   Write(2,*) 'VerVec(:)',VerVec(:), AngOff*180./pi
   !y0Max=0.
   !
   !write(2,"(8x,A,4x,A)") "i, AveSmPow/y0, RimMx/y0,  Mx,  volume;      Barymetric pixel ;    Maximal pixel", &
   !   ';  Closest edge'
   write(2,"(6x,A, 4x,A,9x,';', 4x,A, 4x,A)") "i, AveSmPow/y0, RimMx/y0,  Mx;","t[ms], Position=(N,E,h)[km]",&
      "Q/I, U/I, V/I [%] ;" , "I3, Linpol*100., angle "
   Do i_slice=1,NrPixSmPowTr ! N_smth+1,SumWindw-N_smth,N_smth
      t_ms=(1+i_slice*N_smth)*sample*1000.+t_offsetPow  ! time stamp
!!!!!!!!!   t_shft=sqrt(SUM(CenLoc(:)*CenLoc(:)))*Refrac/c_mps ! in seconds due to signal travel distance
!!!!!!!!!!!!!   t_offsetPow=((Start_time(i_chunk)+SumStrt)*sample-t_shft)*1000.-TimeBase
      !
      ! Do a 3D parabolic fit to interpolate around the maximal pixel
      Call FindInterpolMx(i_slice,d_Mx, SMPowMx,QualMx)
      !
      If(SMPowMx.lt.NoiseLevel) cycle
      !write(2,*) 'd_Mx', i_slice, d_Mx(:), SMPowMx,QualMx(1:3)
      NMx=NMx+1
      If(polar) then
         PixLocPol(:)=CenLocPol(:)+d_Mx(:)*d_loc(:)
         Call Pol2Carth(PixLocPol,PixLoc)
      else
         PixLoc(:)=CenLoc(:)+d_Mx(:)*d_loc(:)
         Call Carth2Pol(PixLoc,PixLocPol)
      Endif
      t_shft=sqrt(SUM(PixLoc(:)*PixLoc(:)))*1000.*Refrac/c_mps ! in seconds due to signal travel distance
      !write(29,"(i6,',',4(f11.5,','),4g13.6,',',3g13.5,2x,L1)") i, t_ms-t_shft, &
      If(Dual) Then
         s=sqrt(MaxSmPowQ(i_slice)*MaxSmPowQ(i_slice)+MaxSmPowU(i_slice)*MaxSmPowU(i_slice))
         X=(ATAN2(MaxSmPowU(i_slice),MaxSmPowQ(i_slice))/2.+AngOff)*180./pi
         Y=t_ms-t_shft
         write(29,"(i6,',',4(f11.5,','),g13.6,',',3f6.2,',',5g13.3)") i_slice, t_ms-t_shft, &
         PixLoc(2)/1000., PixLoc(1)/1000., PixLoc(3)/1000., SMPowMx, &   ! MaxPowPix(:,i)
         MaxSmPowQ(i_slice), MaxSmPowU(i_slice), MaxSmPowV(i_slice), MaxSmPowI3(i_slice) ,s,X
      !   PixLocPol(3)/1000.,PixLocPol(2)*180./pi,PixLocPol(1)*180/pi, d_Mx(:) !, BarySet
         ! calculate polarization for maximum intensity location
         !Call EI_Polariz(Nr_IntFer, IntfNuDim, i_slice)
         Call EI_PolarizW(Nr_IntFer, IntfNuDim, i_slice)
         Write(28,"(i6,',',f11.5,',', 2(g12.4,','), 19(f8.3,',') )") i_slice, t_ms-t_shft, &
            StI, dStI, 100*StI12/StI, 100*dStI12/StI, 100*StQ/StI, 100*dStQ/StI, &
            100*StU/StI, 100*dStU/StI, 100*StV/StI, 100*dStV/StI, &
            100*StI3/StI, 100*dStI3/StI, 100*StU1/StI, 100*dStU1/StI, 100*StV1/StI, 100*dStV1/StI, &
            100*StU2/StI, 100*dStU2/StI, 100*StV2/StI , 100*dStV2/StI, &
            Chi2pDF
         !
      Else
         ! =======================================================================
         write(29,"(i6,',',4(f11.5,','),g13.6,',',3f6.2,',',5g13.3)") i_slice, t_ms-t_shft, &
         PixLoc(2)/1000., PixLoc(1)/1000., PixLoc(3)/1000., SMPowMx
         ! Use the barycenter to interpolate in distance and use only those pixels that have a large intensity
         Call FindBarycenter(i_slice,d_Bar,SMPowBar,QualBar)
         If(SMPowBar.gt.0.) then
            NBar=NBar+1
            If(polar) then
               PixLocPol(:)=CenLocPol(:)+d_Bar(:)*d_loc(:)
               Call Pol2Carth(PixLocPol,PixLoc)
            else
               PixLoc(:)=CenLoc(:)+d_Bar(:)*d_loc(:)
               Call Carth2Pol(PixLoc,PixLocPol)
            Endif
            !
            t_shft=sqrt(SUM(PixLoc(:)*PixLoc(:)))*1000.*Refrac/c_mps ! in seconds due to signal travel distance
            write(28,"(i6,',',4(f11.5,','),4g13.6,',',3g13.5)") i_slice, t_ms-t_shft, &
               PixLoc(2)/1000., PixLoc(1)/1000., PixLoc(3)/1000., SMPowBar, QualBar(1)
            !   PixLocPol(3)/1000.,PixLocPol(2)*180./pi,PixLocPol(1)*180/pi, d_gr(:)
         Else
            d_Bar(:)=9999.
         EndIf
      EndIf
      !
      !i_gr(:)=MaxSmPowGrd(:,i)
      !y0=PixSmPowTr(i,i_gr(1),i_gr(2),i_gr(3))
      !write(2,"(A,I6, 2F7.3, F10.2,f8.0,A,2x,3F7.1,';',3F7.1,';',3I4)") &
      !                  'AveSmPow', i, QualBar(3), QualBar(1), SMPowMx, QualBar(2), ';', &
      !                  d_Bar(:), d_Mx(:)!, RimMaxPix(1:3)
      If(Dual) Then
         write(2,"(A,I6, 2F7.3, F10.2, A, 4F10.5, ';',3F7.1, '%;',3F7.1)") &
                           'MxPow', i_slice, QualBar(3), QualBar(1), SMPowMx, ';', &
                           Y, PixLoc(1)/1000., PixLoc(2)/1000., PixLoc(3)/1000., &
                           100.*MaxSmPowQ(i_slice), 100.*MaxSmPowU(i_slice), 100.*MaxSmPowV(i_slice), &
                           MaxSmPowI3(i_slice) ,100.*s,X
      Else
         write(2,"(A,I6, 2F7.3, F10.2, A, 4F10.5, ';',3F7.1, '%;',3F7.1)") &
                           'MxPow', i_slice, QualBar(3), QualBar(1), SMPowMx, ';', &
                           Y, PixLoc(1)/1000., PixLoc(2)/1000., PixLoc(3)/1000.
      EndIf
      flush(unit=2)
   EndDo
   Close(Unit=28)
   Close(Unit=29)
   If(Dual) Then
      DeAllocate( CMTime_pix )
   EndIf
   write(2,"(A,2I6,2F6.3,A,3F7.3,A,I7,A,F9.1)") 'number of sources in plots:',NMx,NBar, NoiseLevel, RatMax
   If(NMx.gt.1) then
      Call GLEplotControl(PlotType='SourcesPlot', PlotName='IntfMx'//TRIM(txt)//TRIM(OutFileLabel), &
         PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel)//'IntfSpecPowMx'//TRIM(txt) )
      ! write(10,"('gle -d pdf -o ',A,'-InfImaMx_',I1,'.pdf ${UtilDir}/SourcesPlot.gle ${FlashFolder}/files/',A)") &!
      !   TRIM(OutFileLabel), i_eo, TRIM(OutFileLabel)//'IntfSpecPowMx'//TRIM(txt)
      If(Dual) Then
         Call GLEplotControl(PlotType='EIPolariz', PlotName='IntfPol'//TRIM(txt)//TRIM(OutFileLabel), &
            PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel) )
      Else
         If(NBar.gt.1) &
            Call GLEplotControl(PlotType='SourcesPlot', PlotName='IntfBar'//TRIM(txt)//TRIM(OutFileLabel), &
               PlotDataFile=TRIM(DataFolder)//TRIM(OutFileLabel)//'IntfSpecPowBar'//TRIM(txt) )
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
Subroutine FindInterpolMx(i,d_gr,SMPow,Qualty)
   ! Do a 3D parabolic fit to interpolate around the maximal pixel
   !--------------------------------------------
   use constants, only : dp
   Use Interferom_Pars, only : SumStrt, SumWindw, polar, N_pix, IntFer_ant
   Use Interferom_Pars, only : PixelPower, MaxSmPow
   !Use Interferom_Pars, only : PowrThresh
   Use Interferom_Pars, only : NrPixSmPowTr, MaxSmPowGrd, PixSmPowTr
   !Use Interferom_Pars, only : SlcInten, NrSlices, SliceLen, MaxSlcInten, MaxSlcIntenLoc   !IntfBase, IntfDim, IntfPhaseCheck, SumStrt, SumWindw
   Implicit none
   Integer, intent(in) :: i
   Real(dp), intent(out) :: d_gr(1:3), SMPow, Qualty(1:3)
   integer :: j, i_gr(1:3)
   Integer, parameter :: jp(1:3)=(/2,3,1/), jm(1:3)=(/3,1,2/)
   Real(dp) :: y0, Ay(1:3), By(1:3), yp, ym, Det3, Ry(1:3), c_gr(1:3), Parb(-2:+2), Valu(-2:+2)
   Integer :: i_h, i_E, i_N, Nxx
   Real(dp) :: A, B
   Real(dp) :: Paraboloid ! , Sq
   Logical :: Check=.false.
   !Logical :: Check=.true.  !  .false.
   !
   !write(2,*) 'FindInterpolMx',i
   !flush(unit=2)
   i_gr(:)=MaxSmPowGrd(:,i)
   !write(2,*) 'i_gr(:)',i_gr(:), MaxSmPow(0:2)
   !flush(unit=2)
   y0=(PixSmPowTr(i,i_gr(1),i_gr(2),i_gr(3)))
   SMPow=-1.
   !write(2,*) 'y0',i,y0,i_gr(:)
   !flush(unit=2)
   !
   ! Testing for i=1
   !If(i.eq.i_s) Then  ! just for testing
   !   i_gr(1)=1; i_gr(2)=0; i_gr(3)=-1
   !   y0=12345
   !   Do j=1,3
   !      Ay(j)=-j*y0/5.
   !      By(j)=(j-1)*y0/10.
   !      Ry(j)= j*y0/20.
   !   Enddo
   !   Do i_h=-1,1
   !      Do i_N=-1,1
   !         Do i_E=-1,1
   !            d_gr(1)=i_N; d_gr(2)=i_E; d_gr(3)=i_h
   !            PixSmPowTr(i,i_gr(1)+i_N,i_gr(2)+i_E,i_gr(3)+i_h)=Paraboloid(y0,Ay,By,Ry,d_gr)
   !         Enddo
   !      Enddo
   !   Enddo
   !   Write(2,"(A,F9.1,3(3F6.2,','))") 'in:y0,Ay,By,Ry',y0,Ay/y0,By/y0,Ry/y0
   !EndIf
   !write(2,*) 'FindInterpolMx @i_slice=:',i,y0,i_gr(:)
   Do j=1,3  ! fit y=Ax^2/2+Bx+C; A=y"/d^2 & B=y'/2d & x_max= -B/A= -d y'/(2y")
      ! 3D case: y=Sum[A(i)X(i)^2/2 + R'(i,j)X(i)X(j)+B(i)x(i)+C]
      If((i_gr(j).eq.N_pix(j,1)) .or. (i_gr(j).eq.N_pix(j,2))) Then
         Write(2,"(A,i5,A,F9.2,A,3(I5,','))") &
            ' Bordercase',i,', Intensty=',y0,', Max @',MaxSmPowGrd(:,i)
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
               write(2,"(A, 2i3,5F6.2,'; ',5F6.2)") 'area i_h,i_n',i_h,i_N, Valu(:), Parb(:)
               flush(unit=2)
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
      Write(2,*) i,'Det3', Det3, Ay(:), (4*Ry(j)*Ry(j)*Ay(j)/B ,j=1,3)
      Return
   EndIf
   Do j=1,3 ! get relative position of the maximum
      d_gr(j)=-(By(j)*(Ay(jp(j))*Ay(jm(j))-Ry(j)*Ry(j)) - By(jp(j))*(Ay(jm(j))*Ry(jm(j))-Ry(jp(j))*Ry(j)) &
         - By(jm(j))*(Ay(jp(j))*Ry(jp(j))-Ry(jm(j))*Ry(j)) )/Det3
   Enddo
   SMPow=Paraboloid(y0,Ay,By,Ry,d_gr)  ! construct intensity at intepolated position of max
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
   B=SUM(d_gr(:)*d_gr(:))
   If(B.gt.0.5) Then
      write(2,*) 'Interpolation distance shortend:',B
      B=2.*B
      d_gr(:)=d_gr(:)/B
   EndIf
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
   use constants, only : dp, pi, Sample  ! , Refrac, c_mps
   Use Interferom_Pars, only : SumStrt, SumWindw, polar, N_pix, d_loc, IntFer_ant, IntfLead, CenLocPol, CenLoc, RimSmPow
   Use Interferom_Pars, only : PixelPower, MaxSmPow, N_smth
   Use Interferom_Pars, only : xMin, xMax, yMin, yMax, zMin, zMax, tMin, tMax
   Use Interferom_Pars, only : NrPixSmPowTr, MaxSmPowGrd, PixSmPowTr, RatMax
   use Chunk_AntInfo, only : NoiseLevel
   !Use Interferom_Pars, only : SlcInten, NrSlices, SliceLen, MaxSlcInten, MaxSlcIntenLoc   !IntfBase, IntfDim, IntfPhaseCheck, SumStrt, SumWindw
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
Subroutine OutputIntfSlices(i_eo)
   ! Analyze the Slices
   ! Needs
   !--------------------------------------------
   use constants, only : dp, pi, Sample  ! , Refrac, c_mps
   use DataConstants, only : DataFolder, OutFileLabel
   use Chunk_AntInfo, only : CTime_spectr
   use ThisSource, only : Dual, SourcePos
   Use Interferom_Pars, only : CTime_sum, SumStrt, SumWindw, polar, N_pix, d_loc, IntFer_ant, IntfLead, CenLocPol, CenLoc
   Use Interferom_Pars, only : AveInten, AveIntenE, AveIntenN, RimInten, MaxIntfInten, MaxIntfIntenLoc
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
   OPEN(UNIT=29,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'IntfSpecWin'//TRIM(txt)//'.csv')
   If(polar) Then
      write(29,"(g12.4,5(',',i4),',',i8,',',i5,3(',',f9.4))") MaxIntfInten, N_pix(1,2), N_pix(2,1), N_pix(2,2), &
         N_pix(3,1), N_pix(3,2), SumStrt,SumWindw, d_loc(1:2)*180./pi ,d_loc(3)
   Else
      write(29,"(g12.4,5(',',i4),',',i8,',',i5,3(',',f9.4))") MaxIntfInten, N_pix(1,2), N_pix(2,1), N_pix(2,2), &
         N_pix(3,1), N_pix(3,2), SumStrt,SumWindw, d_loc(1:3)
   EndIf
   Do i=0, SumWindw
      j=i+SumStrt
      write(29,"(i6,',',g12.4,',',g12.4)") i, (ABS(CTime_spectr(j,IntFer_ant(1),i_chunk)) )**2, ABS(CTime_sum(j))**2
   Enddo
   Close(Unit=29)
   !
   !--------------------------------------------------------------
   ! Write general info for this picture
   write(2,"(A,G11.3,A,3f9.4,A)", ADVANCE='NO') 'Maximum ',MaxIntfInten,' @ (N,E,h)=(',MaxIntfIntenLoc(:)/1000.,') [km]'
   Call Carth2Pol(MaxIntfIntenLoc(:),PixLocPol)
   SourcePos(:,i_eo+1)=MaxIntfIntenLoc(:)  ! for curtainplot
   write(2,"(A,f9.4,f8.2,f7.2)") ' = (ph,th,R)=',PixLocPol(1)*180/pi,PixLocPol(2)*180./pi,PixLocPol(3)/1000.
   ! Info per slice
   AveIntenN(:,:)=AveIntenN(:,:)*d_loc(1)/AveInten(:,:)
   AveIntenE(:,:)=AveIntenE(:,:)*d_loc(2)/AveInten(:,:)
   RimInten(:,:)=RimInten(:,:)*0.5/(N_pix(1,2)-N_pix(1,1) +N_pix(2,2)-N_pix(2,1)) ! edges accounted for by +1's in sum-range
   AveInten(:,:)=AveInten(:,:)/((N_pix(1,2)-N_pix(1,1)+1)*(N_pix(2,2)-N_pix(2,1)+1))
   If(NrSlices.gt.1) Then
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
