!=================================
!Fine Tuning of source locations, previously found by TRI-D, using the TRI-D techniques.
!- 1) chi^2 optimization of position after finding an improved source time
!- 2) Intensity optimization, using a grid ?or? mimimization procedure using cross corr for $\delta I$.
!==========================================
!-----------------------------------------------
Subroutine EI_GridPeak(i_Peak) ! based on "Subroutine EI_PolarizPeak(i_Peak)" in "EICallRtns.f90"
   use constants, only : dp
   use DataConstants, only : Ant_nrMax
   use Interferom_Pars, only :  Nr_IntFerMx, Nr_IntferCh ! the latter gives # per chunk
   use Interferom_Pars, only : Cnu_p0, Cnu_t0, Cnu_p1, Cnu_t1, IntfNuDim, AntPeak_OffSt
   Use Interferom_Pars, only :  N_Smth, N_fit, Chi2pDF
   use ThisSource, only : ChunkNr, PeakPos, sourcepos, PeakChiSQ, ExclStatNr
   Implicit none
   Integer, intent(in) :: i_Peak
   Integer :: IntfBase, i_chunk, i_1, i_2, i_3, Outpt, j_IntFer, Windw
   Real(dp) :: ChiSq, DelChi(-N_fit:+N_fit,1:2*Nr_IntFerMx)
   Real(dp) :: PartChiSq(1:Nr_IntFerMx) !
   Integer :: PartChiSqInt(1:Nr_IntFerMx)
   Real(dp) :: FitDelay(1:Ant_nrMax), VoxLoc(1:3), del_1, del_2, del_3
   Character(len=8) :: Label
   !
   Outpt=2
   i_chunk=ChunkNr(i_Peak)
   !
   write(2,"(1x,A,i4,A,I5,A,2(F9.4,','),F9.4,A)") 'Peak',i_peak,', Ref. ant. sample=',Peakpos(i_Peak), &
      ', at (N,E,h)=(',SourcePos(:,i_Peak)/1000.,') [km]'
   !
   IntfBase= Peakpos(i_Peak) - IntfNuDim
   Windw=3*N_Smth
   Call GetInterfFitDelay(i_chunk, FitDelay)
   Call EI_PolSetUp(Nr_IntFerCh(i_chunk), IntfBase, i_chunk, SourcePos(1,i_Peak), AntPeak_OffSt(1,i_Peak), &
      Cnu_p0(0,1,i_peak), Cnu_t0(0,1,i_peak), Cnu_p1(0,1,i_peak), Cnu_t1(0,1,i_peak))
   Windw=3*N_Smth
   !write(2,*) 'EI_PolarizPeak: Outpt=',Outpt
   write(Label,"('Pk ',i4.2)") i_Peak
   Call EI_PolGridDel(Nr_IntFerCh(i_chunk), FitDelay, IntfNuDim, i_chunk, SourcePos(1,i_Peak), AntPeak_OffSt(1,i_Peak), &
      Cnu_p0(0,1,i_peak), Cnu_t0(0,1,i_peak), Cnu_p1(0,1,i_peak), Cnu_t1(0,1,i_peak), &
      Outpt, DelChi, Label, ExclStatNr(:,i_peak) )
   PeakChiSQ(i_Peak)=Chi2pDF
   Call WriteDelChiPeak(i_chunk, DelChi,PartChiSq,PartChiSqInt)
   !j_IntFer=PartChiSqInt(Nr_IntFerCh(i_chunk))
   !write(Label,"(i2.2,i3.3)") i_Peak,j_IntFer
   !Call TimeTracePlot(j_IntFer, IntfBase, i_chunk, SourcePos(1,i_Peak), Windw, Label)
   !
   Del_1=1  ! in meter
   Del_2=1
   Del_3=2
   Outpt=1
   write(2,*) i_Peak, 'Del_1=',Del_1,', Del_2=',Del_2, ' [m]'
   Do i_1=-2,2
      Do i_2=-2,2
      Do i_3=-2,2
         VoxLoc(1)=SourcePos(1,i_Peak) + i_1*del_1
         VoxLoc(2)=SourcePos(2,i_Peak) + i_2*del_2
         VoxLoc(3)=SourcePos(3,i_Peak) + i_3*del_3
         write(Label,"('Gr ',I2,',',i2)") i_1,i_2
         Call EI_PolGridDel(Nr_IntFerCh(i_chunk), FitDelay, IntfNuDim, i_chunk, VoxLoc(:), AntPeak_OffSt(1,i_Peak), &
            Cnu_p0(0,1,i_peak), Cnu_t0(0,1,i_peak), Cnu_p1(0,1,i_peak), Cnu_t1(0,1,i_peak),  &
            Outpt, DelChi, Label, ExclStatNr(:,i_peak) )
      Enddo
      Enddo
   Enddo
   !
   Return
End Subroutine EI_GridPeak
!==========================================


!===================================
