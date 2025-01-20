Subroutine LMA_Layout(DistLMA,LMAAnt, NoLMAantennas) !
! Constructs layout of antennas to use in an LMA-like search for the initial guess of the peak positions.
! Needs to be called once before call to LMA_Interfer where actual guess is made
   use ThisSource, only : PeakNrTotal, Peak_eo  ! Nr_Corr, PeakNrTotal, PeakPos, Peak_eo
   !use ThisSource, only : SourcePos
   use Chunk_AntInfo, only : Ant_IDs, Ant_nr, Ant_pos, RefAnt
   use constants, only : dp !,pi,c_mps, sample
   use DataConstants, only : Production
   Implicit none
   Real(dp), Intent(out) :: DistLMA
   Integer, intent(out) :: LMAAnt(1:3,0:1)
   Logical, intent(out) :: NoLMAantennas
   integer :: i_chunk, i_peak, i_eo, i_ant
   integer :: Peak_nr
   !Real(dp) :: RDist, D, D0=6.0d4, H0=1.0d4
   Real(dp) :: LMA_Pos(1:3,0:1)
   Integer, parameter :: i_sel_Max=10
   Real(dp), parameter :: cos_LMA_Max=0.9
   Integer :: SelectAnt(1:i_sel_Max), i_sel
   Real(dp) :: cos_LMA(0:1), D_LMA(0:1), Dist
   Real(dp) :: D_ab, D_ac, D_bc, x_ab, x_ac, x_bc, y_ab, y_ac, y_bc, cos_a, cos_b, cos_c, cos_max
   Integer :: i_a, i_b, i_c, i_LMA
   !real(dp) :: phi(1:2), thet(1:2), cos_th, dt_ab, dt_ac, dx, dy
   !
   ! Basic assumptions:
   ! - really only a single source position is searched for, i.e. the same for all peaks
   !Production=.false.  ! produce printout & stop
   i_chunk=1
   NoLMAantennas=.false.
   ! find three antennas that have largest sum distance, supposedly a best triangle
   ! other possibility: find largest spanned area
   !
   !write(2,*) 'Production:',Production
   DistLMA=50.  ! Keep antennas within LMA distance [m] only from reference antenna
   goto 8
 9 continue ! increase search distance for LMA antennas
   If(.not.production) Write(2,*) i_eo,', cos_LMA(i_eo):', cos_LMA(i_eo),D_LMA(i_eo),', antennas:',LMAAnt(1:3,i_eo)
   If(i_sel.eq.i_sel_Max) Then
      NoLMAantennas=.true.
      Return
   EndIf
   DistLMA=DistLMA+10.  ! Keep antennas within LMA distance [m] only from reference antenna
   If(.not.production) Write(2,*) 'nr antennas=',i_sel,i_eo,', increase LMA-search distance to:',DistLMA
   flush(unit=2)
 8 continue
   Peak_nr=PeakNrTotal
   If(Peak_nr.gt.2) Then
      write(2,*) 'PeakNrTotal:',PeakNrTotal,' exceeds max expected of 2 in LMA_Interfer'
      stop 'LMA_Interfer'
   EndIf
   Do i_Peak=1,Peak_nr ! may be even or odd or both
      i_eo=Peak_eo(i_peak)
      !ReferenceAnt=RefAnt(i_chunk,i_eo)
      i_sel=1
      SelectAnt(i_sel)=RefAnt(i_chunk,i_eo)
      Do i_ant=2,Ant_nr(i_chunk)
         if(mod(Ant_IDs(i_ant,i_chunk),2) .ne. i_eo) cycle       ! keep antenna orientation
         Dist=sqrt(sum( ( Ant_pos(:,i_ant,i_chunk)-Ant_pos(:,RefAnt(i_chunk,i_eo),i_chunk) )**2 ))  ! [m]
         If(Dist .gt. DistLMA) cycle
         i_sel=i_sel+1
         SelectAnt(i_sel)=i_ant
         If(i_sel.eq.i_sel_Max) exit
      EndDo
      If(i_sel.lt.3) goto 9  ! increase DistLMA
      If(.not.production) write(2,*) 'selected:', SelectAnt(1:i_sel)
      ! On almost equilateral triangle ?
      ! and assume planar, i.e. Ant_pos(3,i_ant,i_chunk)=0
      cos_LMA(i_eo)=2.  ! larger than 1.
      i_a=1  ! this should always be the reference antenna
      !Do i_a=1,i_sel-2
         Do i_b=i_a+1,i_sel-1
            x_ab=Ant_pos(1,SelectAnt(i_a),i_chunk)-Ant_pos(1,SelectAnt(i_b),i_chunk)
            y_ab=Ant_pos(2,SelectAnt(i_a),i_chunk)-Ant_pos(2,SelectAnt(i_b),i_chunk)
            D_ab=sqrt(x_ab**2+y_ab**2)  ! [m]
            Do i_c=i_b+1,i_sel
               x_ac=Ant_pos(1,SelectAnt(i_a),i_chunk)-Ant_pos(1,SelectAnt(i_c),i_chunk)
               y_ac=Ant_pos(2,SelectAnt(i_a),i_chunk)-Ant_pos(2,SelectAnt(i_c),i_chunk)
               D_ac=sqrt(x_ac**2+y_ac**2)  ! [m]
               x_bc=Ant_pos(1,SelectAnt(i_b),i_chunk)-Ant_pos(1,SelectAnt(i_c),i_chunk)
               y_bc=Ant_pos(2,SelectAnt(i_b),i_chunk)-Ant_pos(2,SelectAnt(i_c),i_chunk)
               D_bc=sqrt(x_bc**2+y_bc**2)  ! [m]
               cos_c=abs(x_ac*x_bc+y_ac*y_bc)/(D_ac*D_bc)
               cos_a=abs(x_ac*x_ab+y_ac*y_ab)/(D_ac*D_ab)
               cos_b=abs(x_ab*x_bc+y_ab*y_bc)/(D_ab*D_bc)
               cos_max=MAX(cos_a, cos_b, cos_c)  ! should be minimal
               If(cos_max.lt.cos_LMA(i_eo) ) then
                  cos_LMA(i_eo)=cos_max
                  LMAAnt(1:3,i_eo)=(/SelectAnt(i_a), SelectAnt(i_b), SelectAnt(i_c)/)
                  D_LMA(i_eo)=MAX(D_ab, D_ac, D_bc)  ! should be minimal
               EndIF
            EndDo ! i_c=i_b+1,i_sel
         EndDo ! i_b=i_a+1,i_sel-1
      !EndDo ! i_a=1,i_sel-2
      If(cos_LMA(i_eo).gt.cos_LMA_Max) goto 9
      If(.not.production) Write(2,*) i_Peak,i_eo,', cos_LMA(i_eo):', cos_LMA(i_eo),D_LMA(i_eo),', antennas:',LMAAnt(1:3,i_eo)
   EndDo ! i_Peak=1,Peak_nr
   flush(unit=2)
   ! result: LMAAnt(1:3, i_eo)
   Return
End Subroutine LMA_Layout
!=====================================
Subroutine LMA_Interfer(LMAAnt,Thet_LMA,Phi_LMA, NoSourceFound) ! (LMA_pos)
! called from 'SourceFind' with peak nr set to 1 or 2 as well as parities
!  Find an initial guess for the source position through a grid search.
!  Works only in imaging mode since covariance matrix is stored in this mode only. For this
!     the link " X(  XIndx(i,i_Peak) ) = SourcePos(FitPos(i),i_Peak)" is not implemented.
   !use Chunk_AntInfo, only : CTime_spectr, Ant_pos, Time_dim, Ant_Stations, Ant_IDs, Ant_nr
   use ThisSource, only : T2_dim
   use ThisSource, only : Nr_Corr, PeakNrTotal, PeakPos, Peak_eo  ! , Peak_eo, ChunkNr
   use ThisSource, only : CCorr_max, SourcePos
   !use ThisSource, only : Safety
   use FitParams, only : SpaceCov, Sigma_AntT, SearchRangeFallOff
   use Chunk_AntInfo, only : Ant_IDs, Ant_nr, Ant_pos, RefAnt
   use constants, only : dp,pi,c_mps, sample
   use DataConstants, only : Production
    use FFT, only : RFTransform_su,DAssignFFT
   Implicit none
   !Real(dp), intent(out) ::
   Real(dp), Intent(out) :: Thet_LMA,Phi_LMA
   Integer, Intent(in) :: LMAAnt(1:3,0:1)
   Integer, Intent(out) :: NoSourceFound
   integer :: i_chunk, i_peak, i_eo, i_ant, j_corr,i_iter   ! , FitRange_Samples
   !integer :: Peak_nr  !i,j,k,i_loc(1),
   Real(dp) :: LMA_Pos(1:3,0:1), H0=1.0d4  ! RDist, D, D0=6.0d4,
   !Real(dp) :: LMA_Pos(1:3,0:1)
   !Integer, parameter :: i_sel_Max=10
   !Real(dp), parameter :: cos_LMA_Max=0.9
   !Integer :: SelectAnt(1:i_sel_Max), i_sel, LMAAnt(1:3,0:1)
   !Real(dp) :: cos_LMA(0:1), D_LMA(0:1), Dist
   Real(dp) :: D_ab, D_ac, D_bc, x_ab, x_ac, x_bc, y_ab, y_ac, y_bc, cos_a, cos_b, cos_c, cos_max
   Integer :: i_a, i_b, i_c, i_LMA
   real(dp) :: phi(1:2), thet(1:2), cos_th, dt_ab, dt_ac, dx, dy, scl, Num, Denom
   !
   ! Basic assumptions:
   ! - really only a single source position is searched for, i.e. the same for all peaks
   !Production=.false.  ! produce printout & stop
   i_chunk=1
   !write(2,*) 'Safety:',Safety, Sigma_AntT
   !FitRange_Samples=safety
   !safety=20
   ! Find source guess using an LMA-Time of Arrival Difference method, assume a planar array
   NoSourceFound=0
   SpaceCov(:,:)=0.0 ; SpaceCov(1,1)=0.1 ; SpaceCov(2,2)=0.1 ; SpaceCov(3,3)=0.1 ;
   Call RFTransform_su(T2_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   Do i_Peak=1,PeakNrTotal ! may be even or odd or both
      SourcePos(1:3,i_Peak)=(/0.d0,0.d0,H0/)  ! prepare already for later LMA search
      i_eo=Peak_eo(i_peak)
      !ReferenceAnt=LMAAnt(1,i_eo)
      j_corr=0  ! not really used here
      !     free after BuildCC:
      i_a=LMAAnt(1,i_eo)
      i_b=LMAAnt(2,i_eo)
      i_c=LMAAnt(3,i_eo)
      x_ab=Ant_pos(1,i_a,i_chunk)-Ant_pos(1,i_b,i_chunk)
      y_ab=Ant_pos(2,i_a,i_chunk)-Ant_pos(2,i_b,i_chunk)
         D_ab=sqrt(x_ab**2+y_ab**2)  ! [m]
      x_ac=Ant_pos(1,i_a,i_chunk)-Ant_pos(1,i_c,i_chunk)
      y_ac=Ant_pos(2,i_a,i_chunk)-Ant_pos(2,i_c,i_chunk)
         D_ac=sqrt(x_ac**2+y_ac**2)  ! [m]
      j_corr=0  ! not really used here
      Call GetCorrSingAnt(LMAAnt(1,i_eo), J_Corr, i_eo, i_chunk) ! first call with ref ant
      NoSourceFound=0
      Do i_iter=1,1 ! was 1,2 however, not more acceptabl sources found in second run
         j_corr=1  ! not really used here
         scl=0.5*i_iter
         Sigma_AntT=D_ab/(c_mps*sample*SearchRangeFallOff*scl) !*i_iter) ! shorten range for second tryal
         Do i_LMA=2,3
            i_ant=LMAAnt(i_LMA,i_eo)
            Call GetCorrSingAnt( i_ant, J_Corr, i_eo, i_chunk) ! will infact be for a coupl when "Polariz=.true."
            Sigma_AntT=D_ac/(c_mps*sample*SearchRangeFallOff*scl) !*i_iter)
            ! Result in CCorr_max(j_corr,i_Peak) = RtMax
         EndDo !  i_LMA=2,3
         Sigma_AntT=3
         dt_ab=CCorr_max(2,i_Peak)
         dt_ac=CCorr_max(3,i_Peak)
         dx=dt_ac*x_ab-dt_ab*x_ac
         dy=dt_ab*y_ac-dt_ac*y_ab
         Phi(i_Peak)=atan2(dy,dx)
         Num=c_mps*sample*dt_ab
         Denom=x_ab*sin(phi(i_Peak))+y_ab*cos(phi(i_Peak))
         cos_th=Num/Denom
         If(.not.production) write(2,*) 'LMA: phi, cos_th:',phi(i_Peak)*180./pi, cos_th &
            , ', direction:',H0*cos_th*sin(phi(i_Peak)), H0*cos_th*cos(phi(i_Peak)), H0*sqrt(1.d0-cos_th**2)
         !   write(2,*) 'dx,dy',dx,dy,dt_ac*x_ab, dt_ab*x_ac, dt_ab*y_ac, dt_ac*y_ab, x_ab*sin(phi(i_Peak)),y_ab*cos(phi(i_Peak))
         !write(2,*) 'LMA_Interfer:',cos_th,Num,Denom
         If((Abs(Denom) .lt. tiny(Denom)) .or. (cos_th .lt. -0.99) .or. (cos_th .gt. 1.1)) Then
            !D_ab=sqrt(x_ab**2+y_ab**2)  ! [m]
            !write(2,*) 'ab:',c_mps*sample*dt_ab,D_ab,x_ab*cos_th*sin(phi(i_Peak))+y_ab*cos_th*cos(phi(i_Peak))
            !D_ac=sqrt(x_ac**2+y_ac**2)  ! [m]
            !write(2,*) 'ac:',c_mps*sample*dt_ac,D_ac,x_ac*cos_th*sin(phi(i_Peak))+y_ac*cos_th*cos(phi(i_Peak))
            !write(2,*) 'cos_th ac',c_mps*sample*dt_ac/(x_ac*sin(phi(i_Peak))+y_ac*cos(phi(i_Peak)))
            cos_th=0.99
            NoSourceFound=1
         Else If(cos_th .gt.  0.999) Then
            cos_th=0.999
         Else
            exit
         EndIf
      EndDo
      Thet(i_Peak)=acos(cos_th)
      If(.not. Production) Then
         LMA_pos(1,i_eo)=H0*cos_th*sin(phi(i_Peak))
         LMA_pos(2,i_eo)=H0*cos_th*cos(phi(i_Peak))
         LMA_pos(3,i_eo)=H0*sqrt(1.d0-cos_th**2)
         write(2,*) i_Peak,i_eo, 'LMA_estimate:',LMA_pos(:,i_eo)
         D_ab=sqrt(x_ab**2+y_ab**2)  ! [m]
         write(2,*) 'ab:',c_mps*sample*dt_ab,D_ab,x_ab*cos_th*sin(phi(i_Peak))+y_ab*cos_th*cos(phi(i_Peak))
         D_ac=sqrt(x_ac**2+y_ac**2)  ! [m]
         write(2,*) 'ac:',c_mps*sample*dt_ac,D_ac,x_ac*cos_th*sin(phi(i_Peak))+y_ac*cos_th*cos(phi(i_Peak))
      EndIf
   EndDo ! i_Peak=1,Peak_nr
   Call DAssignFFT()                   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   If(PeakNrTotal.eq.2) Then
      Phi_LMA=(Phi(1)+Phi(2))/2.
      Thet_LMA=(Thet(1)+Thet(2))/2.
      NoSourceFound=0
   Else
      Phi_LMA=Phi(1)
      Thet_LMA=Thet(1)
   EndIf
   !safety=FitRange_Samples
   !stop
   Return
End Subroutine LMA_Interfer
!=====================================
!=====================================
!Subroutine SourceTryal(DistMax,i_peakS, NoSourceFound)
Subroutine SourceTryal_v2(DistMax,i_peakS, NoSourceFound)
!  Find an initial guess for the source position through a grid search.
!  Works only in imaging mode since covariance matrix is stored in this mode only. For this
!     the link " X(  XIndx(i,i_Peak) ) = SourcePos(FitPos(i),i_Peak)" is not implemented.
   !use Chunk_AntInfo, only : CTime_spectr, Ant_pos, Time_dim, Ant_Stations, Ant_IDs, Ant_nr
   use DataConstants, only : Production
   use ThisSource, only : Nr_Corr, PeakNr, PeakNrTotal , PeakPos, Peak_eo  ! , Peak_eo, ChunkNr
   use ThisSource, only : CorrAntNrs, CCorr_max, SourcePos, Ant_nrMax, CCorr_Err
   use FitParams, only : Sigma, SpaceCov, CalcHessian
   use DataConstants, only : Production
   use constants, only : dp,pi
   Implicit none
   Real(dp), intent(in) :: DistMax
   Integer, intent(in) :: i_peakS
   Integer, Intent(out) :: NoSourceFound
   Integer, parameter :: NAng= 1
   integer :: i_loc(1),i_chunk, i_peak, i_eo, i_ant, j_corr, StatMax=2000
   integer :: i,j,k, N_ActAnt, NTry, CorrAntNrs_store(1:Ant_nrMax,0:1)
   Integer, save :: LMAAnt(1:3,0:1)
   !logical :: WildGuess=.false.
   !logical :: WildGuess=.true.
   Real(dp) :: SourceTrPos(3,1:(4*NAng+1)*(2*NAng+1)**2), D, D0=8.0d4, Dv(1:3)
   Real(dp) :: CS(1:3*NAng), ChiSq, T_shft,ph, th, dph, dth, ChiSq_min, ChiSq_thr
   Real(dp) :: KalSource(0:3),KalCoVariance(0:3,0:3), ErrMax
   Real(dp) :: DistLMA,Thet_LMA,Phi_LMA,cth,sth,cph,sph,dD
   Integer, save :: cyc=0
   Logical, save :: NoLMAantennas
   !Real(dp), external :: SubRelDist
   !        X(  XIndx(i,i_Peak) ) = SourcePos(FitPos(i),i_Peak)
   !Production=.false.  ! produce printout & stop
   !
   NoSourceFound=1
   Do i_Peak=1,PeakNrTotal
      SourceTrPos(1:3,i_peak)=SourcePos(1:3,i_Peak)
   Enddo
   If(i_peakS.eq.1) Call LMA_Layout(DistLMA,LMAAnt, NoLMAantennas)
   If(NoLMAantennas) return
   Call LMA_Interfer(LMAAnt,Thet_LMA,Phi_LMA, NoSourceFound)  ! changes SourcePos
   Do i_Peak=1,PeakNrTotal
      SourcePos(1:3,i_Peak)=SourceTrPos(1:3,i_peak)
   Enddo
   If(NoSourceFound.ne.0) return
   !Production=.true.  ! less printout
   ! Basic assumptions:
   ! - PeakPos is the same for all peaks
   ! - really only a single source position is searched for, i.e. the same for all peaks
   i_chunk=1
   !
   dph=2.*pi/180.  ! 2 degrees
   dth=2.*pi/180.  ! 2 degrees
   dD=D0*(2.**(-4*Nang))
   NTry=0
   Do i=-NAng,NAng
      th=Thet_LMA+i*dth
      sth=sin(th)
      if(sth.lt.0) cycle
      cth=cos(th)
      Do j=-NAng,NAng
         ph=Phi_LMA+j*dph
         sph=sin(ph)
         cph=cos(ph)
         D=dD
         Do k=1,4*NAng+1
            NTry=NTry+1
            SourceTrPos(1:3,NTry)=(/D*cth*sph,D*cth*cph,D*sth/)
            If(.not. Production) write(2,*) 'ntry:',ntry,D,ph,th
            D=D*2
         enddo
      EndDo
   EndDo
   !
   !Call CheckAntInRange(SourceTrPos,3*NAng,N_ActAnt)
   !
   Dv(1)= NAng*dph/2.
   Dv(2)=D/2.
   Dv(3)=D*sin(Thet_LMA)
   Call CalcCovarianc(Phi_LMA,Dv(2),Dv(3),Dv)
!   Sigma(1)=ABS(SourceTrPos(1,1+NAng)-SourceTrPos(1,2*NAng))/2. !  ! Estimated error
!   Sigma(2)=ABS(SourceTrPos(2,1+NAng)-SourceTrPos(2,2*NAng))/2. !  ! Estimated error
!   Sigma(3)=H/2.  ! Estimated error
!   Do i=1,3
!      SpaceCov(i,i)=Sigma(i)*Sigma(i)
!   Enddo
   If(.not. Production) write(2,*) 'SpaceCov(i,i)a',(SpaceCov(i,i),i=1,3)
   Call BuildCC(StatMax,DistMax)  ! needs guess for covariance matrix
   Call CheckAntInRange(SourceTrPos,NTry,N_ActAnt)
   If(.not. Production) write(2,*) 'N_ActAnt=',N_ActAnt
   Call OptSrcPos(SourceTrPos,NTry,ChiSq_min)
   !
   ph=atan2(SourcePos(1,1),SourcePos(2,1))
   D=sqrt(SourcePos(1,1)*SourcePos(1,1)+SourcePos(2,1)*SourcePos(2,1))
   Dv(1)=dph/2.
   Dv(2)=D/2.
   Dv(3)=SourcePos(3,1)
   Call CalcCovarianc(ph,D,SourcePos(3,1),Dv)
   If(.not. Production) write(2,*) 'SpaceCov(i,i)b',(SpaceCov(k,k),k=1,3)
   !If(ChiSq_min.gt.ChiSq_thr) then
      !Call KalmanFilt(KalSource,KalCoVariance)
      !sourcepos(1:3,1)=KalSource(1:3)*4000./KalSource(3)
      !SpaceCov(1:3,1:3)=KalCoVariance(1:3,1:3)   ! 3.*
      !Write(2,*) 'Kalman guess',sourcepos(:,1),ChiSq_Thr,ChiSq_min
      !DistMax=0.1
      !CalcHessian=.true.
      !Call SourceFitCycle(StatMax,DistMax)
   !EndIf
   !Stop 'SourceTryal2-end'
   !Production=.true.  ! produce printout
   If(.not. Production) Then
      If( cyc .gt. 10) stop "SourceTryal"
      cyc=cyc+1
   EndIf
   Return
!End Subroutine SourceTryal
End Subroutine SourceTryal_v2
!=====================================
!=====================================
Subroutine SourceTryal_v1(DistMax,i_peakS, NoSourceFound)
!  Obsolete ??
!Subroutine SourceTryal(DistMax,i_peakS, NoSourceFound)
!   Integer, intent(in) :: i_peakS
!   Integer, Intent(out) :: NoSourceFound
!  Find an initial guess for the source position through a grid search.
!  Works only in imaging mode since covariance matrix is stored in this mode only. For this
!     the link " X(  XIndx(i,i_Peak) ) = SourcePos(FitPos(i),i_Peak)" is not implemented.
   !use Chunk_AntInfo, only : CTime_spectr, Ant_pos, Time_dim, Ant_Stations, Ant_IDs, Ant_nr
   use DataConstants, only : Production
   use ThisSource, only : Nr_Corr, PeakNr, PeakNrTotal , PeakPos, Peak_eo  ! , Peak_eo, ChunkNr
   use ThisSource, only : CorrAntNrs, CCorr_max, SourcePos, Ant_nrMax, CCorr_Err
   use FitParams, only : Sigma, SpaceCov, CalcHessian
   use DataConstants, only : Production
   use constants, only : dp,pi
   Implicit none
   Real(dp), intent(in) :: DistMax
   Integer, intent(in) :: i_peakS
   Integer, Intent(out) :: NoSourceFound
   !Integer, parameter :: NAng= 16
   Integer, parameter :: NAng= 32
   integer :: i_loc(1),i_chunk, i_peak,i_ref, i_eo, i_ant, j_corr, StatMax=2000
   integer :: i,j,k, N_ActAnt, NTry, CorrAntNrs_store(1:Ant_nrMax,0:1)
   logical :: WildGuess=.false.
   !logical :: WildGuess=.true.
   Real(dp) :: SourceTrPos(3,1:6*NAng), RDist, D, D0=6.0d4, h, H0=4.0d3, Dv(1:3)
   Real(dp) :: CS(1:3*NAng), ChiSq, T_shft,ph, ph0, dph, Fact, ChiSq_min, ChiSq_thr
   Real(dp) :: KalSource(0:3),KalCoVariance(0:3,0:3), ErrMax
   Real(dp) :: LMA_pos(1:3,0:1)
   !Real(dp), external :: SubRelDist
   !        X(  XIndx(i,i_Peak) ) = SourcePos(FitPos(i),i_Peak)
   !integer, external :: XIndx
   !
   !Do i_Peak=1,PeakNrTotal
   !   SourceTrPos(1:3,i_peak)=SourcePos(1:3,i_Peak)
   !Enddo
   !Call LMA_Interfer(LMA_pos)  ! changes SourcePos
   !Do i_Peak=1,PeakNrTotal
   !   SourcePos(1:3,i_Peak)=SourceTrPos(1:3,i_peak)
   !Enddo
   ! Basic assumptions:
   ! - PeakPos is the same for all peaks
   ! - really only a single source position is searched for, i.e. the same for all peaks
   ChiSq_thr=50.
   i_chunk=1
   ! DistMax=0.3
   !Production=.false.  ! produce printout & stop
   ph0=atan2(SourcePos(1,1),SourcePos(2,1))
   If(WildGuess) then
      h=H0
      Do i_Peak=1,PeakNrTotal
         SourcePos(1,i_Peak)=0
         SourcePos(2,i_Peak)=0
         SourcePos(3,i_Peak)=H
      Enddo
      D=D0
      Sigma(1)=D  ! Estimated error
      Sigma(2)=D  ! Estimated error
      Sigma(3)=H  ! Estimated error
      SpaceCov=0.0
      Do i=1,3
         SpaceCov(i,i)=Sigma(i)*Sigma(i)
      Enddo
      dph=2.*pi/NAng
   Else  ! Use source pos as central guess and search over 180 deg
      D=sqrt(SourcePos(1,1)*SourcePos(1,1)+SourcePos(2,1)*SourcePos(2,1))
      h=SourcePos(3,1)
      Dv(1)=pi/2.
      Dv(2)=D ! /2.
      Dv(3)=h ! /2.
      Call CalcCovarianc(ph0,D,h,Dv)
      ! dph=pi/NAng !original   ;
      dph=2.*pi/NAng  ! modified 02/05/2023, appears to give more imaged sources than the old choice
   Endif
   !
   Do i_Peak=1,PeakNrTotal  ! to start following search in the prescribed direction
      SourcePos(2,i_Peak)=-SourcePos(2,i_Peak) ! changes signeach time
   Enddo
   Do k=1,3 ! Search for a decent source location at three distances
      Do i=1,4 ! Search for a decent source location in four quadrants
         j=2*mod(i,2)-1
         Do i_Peak=1,PeakNrTotal
            SourcePos(1,i_Peak)=j*SourcePos(1,i_Peak) ! changes sign for even values of i
            SourcePos(2,i_Peak)=-SourcePos(2,i_Peak) ! changes signeach time
         Enddo
         Call BuildCC(StatMax,DistMax/2.)
         If(.not. Production) write(2,*) 'SourceTryal(DistMax), BuildCC called:', DistMax,PeakNrTotal
         flush(unit=2)
         ErrMax=0.
         Do i_Peak=1,PeakNrTotal
            i_eo=Peak_eo(i_peak)
            CorrAntNrs_store(1:Nr_Corr(i_eo,i_chunk),i_eo) =CorrAntNrs(1:Nr_Corr(i_eo,i_chunk),i_eo, i_chunk) ! needed since modified in CheckAntInRange
            If(.not. Production) write(2,*) 'CC-errors', Nr_Corr(i_eo,i_chunk), CCorr_Err(1:Nr_Corr(i_eo,i_chunk),i_Peak)
            ErrMax=Max(ErrMax,maxval(CCorr_Err(1:Nr_Corr(i_eo,i_chunk),i_Peak)))
         EndDo
         If(.not. Production) write(2,*) 'ErrMax:', ErrMax
         flush(unit=2)
         If(ErrMax.lt.1) goto 3 ! this is likely a decent source and exit the search
      enddo  ! loop over azimuth directions
      Do i_Peak=1,PeakNrTotal
         H=H*3
         SourcePos(3,i_Peak)=H
      EndDo
   EndDo ! loop over sourceheights
   !flush(unit=2)
   !
 3 NTry=0
 1 Continue
   NTry=Ntry+1
   ph=ph0-0.5*NAng*dph
   If(.not. Production) write(2,*) 'D,H:',D,H,DistMax/2., ph,ph+NAng*dph, NAng
   Do k=1,Nang
   SourceTrPos(1,k)=D*sin(ph+k*dph)
   SourceTrPos(2,k)=D*cos(ph+k*dph)
   SourceTrPos(3,k)=H
   enddo
   !
   If(.not. Production) write(2,*) 'SourceTryal,NTry:', NTry, ph, D
   Call CheckAntInRange(SourceTrPos,NAng,N_ActAnt)
   If(.not. Production) write(2,*) 'N_ActAnt=',N_ActAnt, Nr_Corr(i_eo,i_chunk)
   If(N_ActAnt.lt.4) then
      If(.not. Production) write(2,*) 'too few active antennas,',N_ActAnt,', for (D,ph)=',D,ph*180./pi
      If(NTry.gt. 4)  Then
         write(2,*) 'Too many tryals in SourceTryal, formerly: stop SourceTryal2'
         Return
      EndIf    !stop 'SourceTryal2'
      If(WildGuess) then
         D=D/2.
      Else
         !dph=dph/2.
         D=D/3.
      Endif
      !write(2,*) 'restore CorrAntNrs'
      flush(unit=2)
      Do i_Peak=1,PeakNrTotal
         i_eo=Peak_eo(i_peak)
         CorrAntNrs(1:Nr_Corr(i_eo,i_chunk),i_eo, i_chunk)= CorrAntNrs_store(1:Nr_Corr(i_eo,i_chunk),i_eo) ! needed since modified in CheckAntInRange
      EndDo
      !write(2,*) 'restored CorrAntNrs'
      goto 1
   Endif
   !ChiSq_min=ChiSq_thr
   !
   Call OptSrcPos(SourceTrPos,NAng,ChiSq_min)  ! replaces SourcePos
   !
   !If(ChiSq_min.gt.ChiSq_thr) then
   !   R=R/3.
   !   If(R.lt.1000.) stop 'SourceTryal'
   !   goto 1
   !EndIf
   ! Next round
   If(.not. Production) write(2,*) 'source:',SourcePos(1,1),SourcePos(2,1)
   ph0=atan2(SourcePos(1,1),SourcePos(2,1))
   If(WildGuess) then
      dph=1.8*dph/NAng
      D=D0/3.
   Else
      dph=1.8*dph/NAng
      D=D/2.
   Endif
   If(.not. Production) write(2,*) 'direction=',Ph0,dph,D,WildGuess
   !   Call KalmanFilt(KalSource,KalCoVariance)
   !   sourcepos(1:3,1)=KalSource(1:3)*4000./KalSource(3)
   !   SpaceCov(1:3,1:3)=KalCoVariance(1:3,1:3)   ! 3.*
   !   Write(2,*) 'Kalman guess',sourcepos(:,1),ChiSq_Thr,ChiSq_min
 2 Continue
   ph=ph0 - 0.5*NAng*dph
   Do k=1,NAng
   SourceTrPos(1,k)=D*sin(ph+dph*k)
   SourceTrPos(2,k)=D*cos(ph+dph*k)
   SourceTrPos(3,k)=H
   SourceTrPos(1,k+NAng)=2.*D*sin(ph +dph*k)
   SourceTrPos(2,k+NAng)=2.*D*cos(ph +dph*k)
   SourceTrPos(3,k+NAng)=H
   SourceTrPos(1,k+2*NAng)=3.*D*sin(ph +dph*k)
   SourceTrPos(2,k+2*NAng)=3.*D*cos(ph +dph*k)
   SourceTrPos(3,k+2*NAng)=H
   enddo
   SourceTrPos(:,3*NAng+1:6*NAng)=  2.* SourceTrPos(:,1:3*NAng)
   !
   !Call CheckAntInRange(SourceTrPos,3*NAng,N_ActAnt)
   !
   Dv(1)= NAng*dph/2.
   Dv(2)=D
   Dv(3)=h
   Call CalcCovarianc(ph0,D,h,Dv)
!   Sigma(1)=ABS(SourceTrPos(1,1+NAng)-SourceTrPos(1,2*NAng))/2. !  ! Estimated error
!   Sigma(2)=ABS(SourceTrPos(2,1+NAng)-SourceTrPos(2,2*NAng))/2. !  ! Estimated error
!   Sigma(3)=H/2.  ! Estimated error
!   Do i=1,3
!      SpaceCov(i,i)=Sigma(i)*Sigma(i)
!   Enddo
   If(.not. Production) write(2,*) 'ph0',ph0,'SpaceCov(i,i)a',(SpaceCov(i,i),i=1,3)
   Call BuildCC(StatMax,DistMax)  ! needs guess for covariance matrix
   Call CheckAntInRange(SourceTrPos,6*NAng,N_ActAnt)
   If(.not. Production) write(2,*) 'N_ActAnt=',N_ActAnt
   Call OptSrcPos(SourceTrPos,6*NAng,ChiSq_min)
   !
   ph0=atan2(SourcePos(1,1),SourcePos(2,1))
   D=sqrt(SourcePos(1,1)*SourcePos(1,1)+SourcePos(2,1)*SourcePos(2,1))
   Dv(1)=dph/2.
   Dv(2)=D/2.
   Dv(3)=h/2.
   Call CalcCovarianc(ph0,D,SourcePos(3,1),Dv)
   If(.not. Production) write(2,*) 'SpaceCov(i,i)b',(SpaceCov(k,k),k=1,3)
   !If(ChiSq_min.gt.ChiSq_thr) then
      !Call KalmanFilt(KalSource,KalCoVariance)
      !sourcepos(1:3,1)=KalSource(1:3)*4000./KalSource(3)
      !SpaceCov(1:3,1:3)=KalCoVariance(1:3,1:3)   ! 3.*
      !Write(2,*) 'Kalman guess',sourcepos(:,1),ChiSq_Thr,ChiSq_min
      !DistMax=0.1
      !CalcHessian=.true.
      !Call SourceFitCycle(StatMax,DistMax)
   !EndIf
   !Stop 'SourceTryal2-end'
   !Production=.true.  ! produce printout
   If(.not. Production) stop "SourceTryal"
   Return
!End Subroutine SourceTryal
End Subroutine SourceTryal_v1
!=====================================
!=====================================
Subroutine SourceTryal_v0(DistMax,i_peakS, NoSourceFound)
! Obsolete ??
!   Integer, intent(in) :: i_peakS
!   Integer, Intent(out) :: NoSourceFound
!  Find an initial guess for the source position through a grid search.
!  Works only in imaging mode since covariance matrix is stored in this mode only. For this
!     the link " X(  XIndx(i,i_Peak) ) = SourcePos(FitPos(i),i_Peak)" is not implemented.
   !use Chunk_AntInfo, only : CTime_spectr, Ant_pos, Time_dim, Ant_Stations, Ant_IDs, Ant_nr
   use DataConstants, only : Production
   use ThisSource, only : Nr_Corr, PeakNr, PeakNrTotal , PeakPos, Peak_eo  ! , Peak_eo, ChunkNr
   use ThisSource, only : CorrAntNrs, CCorr_max, SourcePos, Ant_nrMax
   use FitParams, only : Sigma, SpaceCov, CalcHessian
   use DataConstants, only : Production
   use constants, only : dp,pi
   Implicit none
   Real(dp), intent(in) :: DistMax
   Integer, intent(in) :: i_peakS
   Integer, Intent(out) :: NoSourceFound
   Integer, parameter :: NAng= 16
   integer :: i_loc(1),i_chunk, i_peak,i_ref, i_eo, i_ant, j_corr, StatMax=2000
   integer :: i,j,k, N_ActAnt, NTry, CorrAntNrs_store(1:Ant_nrMax,0:1)
   logical :: WildGuess=.false.
   !logical :: WildGuess=.true.
   Real(dp) :: SourceTrPos(3,1:6*NAng), RDist, D, D0=6.0d4, h, H0=4.0d3, Dv(1:3)
   Real(dp) :: CS(1:3*NAng), ChiSq, T_shft,ph, ph0, dph, Fact, ChiSq_min, ChiSq_thr
   Real(dp) :: KalSource(0:3),KalCoVariance(0:3,0:3)
   !Real(dp), external :: SubRelDist
   !        X(  XIndx(i,i_Peak) ) = SourcePos(FitPos(i),i_Peak)
   !integer, external :: XIndx
   !
   ! Basic assumptions:
   ! - PeakPos is the same for all peaks
   ! - really only a single source position is searched for, i.e. the same for all peaks
   ChiSq_thr=50.
   i_chunk=1
   ! DistMax=0.3
   !Production=.false.  ! produce printout
   ph0=atan2(SourcePos(1,1),SourcePos(2,1))
   If(WildGuess) then
      h=H0
      Do i_Peak=1,PeakNrTotal
         SourcePos(1,i_Peak)=0
         SourcePos(2,i_Peak)=0
         SourcePos(3,i_Peak)=H
      Enddo
      D=D0
      Sigma(1)=D  ! Estimated error
      Sigma(2)=D  ! Estimated error
      Sigma(3)=H  ! Estimated error
      SpaceCov=0.0
      Do i=1,3
         SpaceCov(i,i)=Sigma(i)*Sigma(i)
      Enddo
      dph=2.*pi/NAng
   Else  ! Use source pos as central guess and search over 180 deg
      D=sqrt(SourcePos(1,1)*SourcePos(1,1)+SourcePos(2,1)*SourcePos(2,1))
      h=SourcePos(3,1)
      Dv(1)=pi/2.
      Dv(2)=D ! /2.
      Dv(3)=h ! /2.
      Call CalcCovarianc(ph0,D,h,Dv)
      dph=pi/NAng
   Endif
   Call BuildCC(StatMax,DistMax/2.)
   !write(2,*) 'SourceTryal(DistMax), BuildCC called:', DistMax,PeakNrTotal
   !flush(unit=2)
   !
   NTry=0
   Do i_Peak=1,PeakNrTotal
      i_eo=Peak_eo(i_peak)
      CorrAntNrs_store(1:Nr_Corr(i_eo,i_chunk),i_eo) =CorrAntNrs(1:Nr_Corr(i_eo,i_chunk),i_eo, i_chunk) ! needed since modified in CheckAntInRange
   EndDo
 1 Continue
   NTry=Ntry+1
   ph=ph0-0.5*NAng*dph
   !write(2,*) 'D,H:',D,H,DistMax/2.
   Do k=1,Nang
   SourceTrPos(1,k)=D*sin(ph+k*dph)
   SourceTrPos(2,k)=D*cos(ph+k*dph)
   SourceTrPos(3,k)=H
   enddo
   !write(2,*) ((ph+k*dph),k=1,Nang)
   !
   !write(2,*) 'SourceTryal,NTry:', NTry, ph, D
   Call CheckAntInRange(SourceTrPos,NAng,N_ActAnt)
   If(.not. Production) write(2,*) 'N_ActAnt=',N_ActAnt
   If(N_ActAnt.lt.4) then
      write(2,*) 'too few active antennas,',N_ActAnt,', for (D,ph)=',D,ph*180./pi
      If(NTry.gt. 4)  Then
         write(2,*) 'Too many tryals in SourceTryal, formerly: stop SourceTryal2'
         Return
      EndIf    !stop 'SourceTryal2'
      If(WildGuess) then
         D=D/2.
      Else
         dph=dph/2.
      Endif
      !write(2,*) 'restore CorrAntNrs'
      flush(unit=2)
      Do i_Peak=1,PeakNrTotal
         i_eo=Peak_eo(i_peak)
         CorrAntNrs(1:Nr_Corr(i_eo,i_chunk),i_eo, i_chunk)= CorrAntNrs_store(1:Nr_Corr(i_eo,i_chunk),i_eo) ! needed since modified in CheckAntInRange
      EndDo
      !write(2,*) 'restored CorrAntNrs'
      goto 1
   Endif
   !ChiSq_min=ChiSq_thr
   !
   Call OptSrcPos(SourceTrPos,NAng,ChiSq_min)  ! replaces SourcePos
   !
   !If(ChiSq_min.gt.ChiSq_thr) then
   !   R=R/3.
   !   If(R.lt.1000.) stop 'SourceTryal'
   !   goto 1
   !EndIf
   ! Next round
   !write(2,*) 'source:',SourcePos(1,1),SourcePos(2,1)
   ph0=atan2(SourcePos(1,1),SourcePos(2,1))
   If(WildGuess) then
      dph=1.8*dph/NAng
      D=D0/3.
   Else
      dph=1.8*dph/NAng
      D=D/2.
   Endif
   !write(2,*) 'direction=',th,NAng*th/(2.*pi),fact
   !   Call KalmanFilt(KalSource,KalCoVariance)
   !   sourcepos(1:3,1)=KalSource(1:3)*4000./KalSource(3)
   !   SpaceCov(1:3,1:3)=KalCoVariance(1:3,1:3)   ! 3.*
   !   Write(2,*) 'Kalman guess',sourcepos(:,1),ChiSq_Thr,ChiSq_min
 2 Continue
   ph=ph0 - 0.5*NAng*dph
   Do k=1,NAng
   SourceTrPos(1,k)=D*sin(ph+dph*k)
   SourceTrPos(2,k)=D*cos(ph+dph*k)
   SourceTrPos(3,k)=H
   SourceTrPos(1,k+NAng)=2.*D*sin(ph +dph*k)
   SourceTrPos(2,k+NAng)=2.*D*cos(ph +dph*k)
   SourceTrPos(3,k+NAng)=H
   SourceTrPos(1,k+2*NAng)=3.*D*sin(ph +dph*k)
   SourceTrPos(2,k+2*NAng)=3.*D*cos(ph +dph*k)
   SourceTrPos(3,k+2*NAng)=H
   enddo
   SourceTrPos(:,3*NAng+1:6*NAng)=  2.* SourceTrPos(:,1:3*NAng)
   !
   !Call CheckAntInRange(SourceTrPos,3*NAng,N_ActAnt)
   !
   Dv(1)= NAng*dph/2.
   Dv(2)=D
   Dv(3)=h
   Call CalcCovarianc(ph0,D,h,Dv)
!   Sigma(1)=ABS(SourceTrPos(1,1+NAng)-SourceTrPos(1,2*NAng))/2. !  ! Estimated error
!   Sigma(2)=ABS(SourceTrPos(2,1+NAng)-SourceTrPos(2,2*NAng))/2. !  ! Estimated error
!   Sigma(3)=H/2.  ! Estimated error
!   Do i=1,3
!      SpaceCov(i,i)=Sigma(i)*Sigma(i)
!   Enddo
   !write(2,*) 'ph0',ph0,'SpaceCov(i,i)a',(SpaceCov(i,i),i=1,3)
   Call BuildCC(StatMax,DistMax)  ! needs guess for covariance matrix
   Call CheckAntInRange(SourceTrPos,6*NAng,N_ActAnt)
   If(.not. Production) write(2,*) 'N_ActAnt=',N_ActAnt
   Call OptSrcPos(SourceTrPos,6*NAng,ChiSq_min)
   !
   ph0=atan2(SourcePos(1,1),SourcePos(2,1))
   D=sqrt(SourcePos(1,1)*SourcePos(1,1)+SourcePos(2,1)*SourcePos(2,1))
   Dv(1)=dph/2.
   Dv(2)=D/2.
   Dv(3)=h/2.
   Call CalcCovarianc(ph0,D,SourcePos(3,1),Dv)
   !write(2,*) 'SpaceCov(i,i)b',(SpaceCov(k,k),k=1,3)
   !If(ChiSq_min.gt.ChiSq_thr) then
      !Call KalmanFilt(KalSource,KalCoVariance)
      !sourcepos(1:3,1)=KalSource(1:3)*4000./KalSource(3)
      !SpaceCov(1:3,1:3)=KalCoVariance(1:3,1:3)   ! 3.*
      !Write(2,*) 'Kalman guess',sourcepos(:,1),ChiSq_Thr,ChiSq_min
      !DistMax=0.1
      !CalcHessian=.true.
      !Call SourceFitCycle(StatMax,DistMax)
   !EndIf
   !Stop 'SourceTryal2-end'
   !Production=.true.  ! produce printout
   Return
End Subroutine SourceTryal_v0
!=====================================
!=====================================
Subroutine CheckAntInRange(SourceTrPos,N_try,N_ActAnt)
! Checks which antennat in the list "CorrAntNrs(2:j_corr,i_eo, i_chunk)" are able to see pulses from
!     all source locations given in "SourceTrPos(3,1:N_try)".
   !use Chunk_AntInfo, only : CTime_spectr, Ant_pos, Time_dim, Ant_Stations, Ant_IDs, Ant_nr
   use ThisSource, only : Nr_Corr, PeakNr, PeakNrTotal, PeakPos, Peak_eo, ChunkNr
   use ThisSource, only : CorrAntNrs, CCorr_max, SourcePos, Safety
   use constants, only : dp
   Implicit none
   Integer, intent(IN) :: N_try
   Real(dp), intent(IN) :: SourceTrPos(3,1:N_try)
   Integer, intent(OUT) :: N_ActAnt
   integer :: i_peak, i_ref, i_Ant, i_chunk, j_corr, i_eo, k
   Real(dp) :: TS_max,TS_min
   !
   !write(2,*) 'N_try=',N_try
   N_ActAnt=0
   Do i_Peak=1,PeakNrTotal
      i_eo=Peak_eo(i_peak)
      i_chunk=ChunkNr(i_peak)
      j_corr=1
      i_ref=CorrAntNrs(j_corr,i_eo, i_chunk)
      !write(2,*) 'CheckAntInRange:', Safety, i_Peak, PeakNrTotal
      !write(2,*) 'CheckAntInRange:',Nr_Corr(i_eo,i_chunk),';', CorrAntNrs(2:Nr_Corr(i_eo,i_chunk),i_eo, i_chunk)
      Do j_corr=2,Nr_Corr(i_eo,i_chunk)
         !i_ant=CorrAntNrs(j_corr,i_eo, i_chunk)
         Call AntInRangeChk(i_peak,j_corr,i_ref,SourceTrPos,N_try,TS_max,TS_min)
         !write(2,*) 'TS',j_corr,i_ref,TS_max,TS_min
         If((TS_max.gt.Safety) .or. (TS_min.lt.-Safety)) then
            CorrAntNrs(j_corr,i_eo, i_chunk)=-1
         Else
            N_ActAnt=N_ActAnt+1
         Endif
      Enddo
      !write(2,*) 'CorrAntNrs:',i_eo,Nr_Corr(i_eo,i_chunk), (CorrAntNrs(j_corr,i_eo, i_chunk),j_corr=2,Nr_Corr(i_eo,i_chunk))
   EndDo
   Return
End Subroutine CheckAntInRange
!=====================================
Subroutine AntInRangeChk(i_peak,j_corr,i_ref,SourceTrPos,N_try,TS_max,TS_min)
!  Calculates max and min of pulse position shift for a given antenna (w.r.t. reference) and
!     list of source positions
   !use Chunk_AntInfo, only : CTime_spectr, Ant_pos, Time_dim, Ant_Stations, Ant_IDs, Ant_nr
   !use ThisSource, only : Nr_Corr, PeakNr, PeakNrTotal, PeakPos, Peak_eo, ChunkNr
   !use ThisSource, only : CorrAntNrs, CCorr_max, SourcePos, Safety
   use ThisSource, only : Nr_Corr, PeakNr, PeakNrTotal, PeakPos, Peak_eo, ChunkNr, T_Offset, CorrAntNrs
   use constants, only : dp
   Implicit none
   integer, intent(IN) :: i_peak,j_corr,i_ref, N_try
   Real(dp), intent(IN) :: SourceTrPos(3,1:N_try)
   Real(dp), intent(OUT) :: TS_max,TS_min
   integer :: k, i_chunk, i_ant
   Real(dp) :: TS(1:N_try)
   Real(dp), external :: SubRelDist
   !
   !i_ref=CorrAntNrs(j_corr,i_eo, i_chunk)
   i_chunk=ChunkNr(i_peak)
   i_ant=CorrAntNrs(j_corr,Peak_eo(i_peak), i_chunk)
   !write(2,*) 'AntInRangeChk:', i_ant, j_corr, i_peak, Peak_eo(i_peak), i_chunk,N_try, T_Offset(j_corr,i_Peak)
   !flush (unit=2)
   Do k=1,N_try
      TS(k)=SubRelDist(SourceTrPos(1,k),i_ant,i_chunk) - SubRelDist(SourceTrPos(1,k),i_ref,i_chunk) - T_Offset(j_corr,i_Peak)
      !k=k+1
   Enddo
   TS_max=MAXVAL(TS)
   TS_min=MINVAL(TS)
   !write(2,*) TS_max,TS_min, 'TS:',TS
   !   i_loc=MinLoc(CS(0:3*NAng) ) - 1 ! to correct for 0 as first element
End Subroutine AntInRangeChk
!=====================================
Subroutine OptSrcPos(SourceTrPos,N_try,ChiSq_min)
   !use Chunk_AntInfo, only : CTime_spectr, Ant_pos, Time_dim, Ant_Stations, Ant_IDs, Ant_nr
   use DataConstants, only : Production
   use ThisSource, only : Nr_Corr, PeakNrTotal, PeakPos, Peak_eo, ChunkNr ! , PeakNr
   use ThisSource, only : CorrAntNrs, CCorr_max, SourcePos, T_Offset ! , Safety
   use constants, only : dp
   Implicit none
   integer, intent(IN) :: N_try
   Real(dp), intent(IN) :: SourceTrPos(3,1:N_try)
   Real(dp), intent(InOUT) :: ChiSq_min
   integer :: k, i_chunk, i_eo, j_corr, i_peak,i_ref, i_Ant
   Integer :: i_loc(1)
   Real(dp) :: CS(1:N_try), RDist, T_shft  , R(50)
   Real(dp), external :: SubRelDist
   !
   CS(:)=0.
   Do i_Peak=1,PeakNrTotal
      i_eo=Peak_eo(i_peak)
      i_chunk=ChunkNr(i_peak)
      j_corr=1
      i_ref=CorrAntNrs(j_corr,i_eo, i_chunk)
      If(.not. Production) write(2,*) 'OptSrcPos:', i_ant, j_corr, i_peak, Peak_eo(i_peak), i_chunk
      !flush (unit=2)
      Do k=1,N_try  ! calculate chi^2 for the different source tryal positions
         Rdist=SubRelDist(SourceTrPos(1,k),i_ref,i_chunk)
         Do j_corr=2,Nr_Corr(i_eo,i_chunk)
            i_ant=CorrAntNrs(j_corr,i_eo, i_chunk)
            If(i_ant.lt.0) cycle
            !Call RelDist(SourceTrPos(1,k),Ant_pos(1,i_ant,i_chunk),RDist)
            T_shft=SubRelDist(SourceTrPos(1,k),i_ant,i_chunk) - Rdist - T_Offset(j_corr,i_Peak) ! subtract the value for the reference antenna
            !if(k.eq.1) write(2,*) k,j_corr,T_shft, Rdist, T_Offset(j_corr,i_Peak),CCorr_max(j_corr,i_Peak)
            !If(J_corr.lt.10) R(j_corr)=CCorr_max(j_corr,i_Peak) - T_shft
            CS(k)=CS(k) + (CCorr_max(j_corr,i_Peak) - T_shft)*(CCorr_max(j_corr,i_Peak) - T_shft)
            !R(j_corr)=(CCorr_max(j_corr,i_Peak)-T_shft)*5
         Enddo
         !Write(2,*) k,'residuals:',CS(k),R(2:20)
         !write(2,*) k,SourceTrPos(:,k),CS(k)
      Enddo
   EndDo  ! i_Peak=1,PeakNrTotal
   i_loc=MinLoc(CS(1:N_try) )
   !write(2,*) 'T_shft=',CCorr_max(2:Nr_Corr(i_eo,i_chunk),i_peak)
   If(.not. Production) then
      write(2,'(A)') ' chi^2:'
      write(2,'(16F9.1)') CS(1:N_try)
   endif
   !write(2,*) 'chi^2=',CS(NAng:2*NAng)
   !write(2,*) 'chi^2=',CS(2*NAng:3*NAng)
   If(.not. Production) &
         write(2,"(A,g12.4,A,i3,A,3f7.1,A,2i3)") 'chi^2_min=',CS(i_loc(1)),', grid=', i_loc(1), &
         ', source@',SourceTrPos(:,i_loc(1))/1000., '[km], N_antenna=',Nr_Corr(0,i_chunk),Nr_Corr(1,i_chunk)
   !If(CS(i_loc(1)).lt. ChiSq_min) then
   Do i_Peak=1,PeakNrTotal
      SourcePos(:,i_Peak) = SourceTrPos(:,i_loc(1)) ! source position for first iteration
   Enddo
   !EndIf
   ChiSq_min=CS(i_loc(1))
   !
End Subroutine OptSrcPos
! ======================================
Subroutine SearchWin(i_ref, i_ant, i_chunk, SrcPos, EEst)
   use Chunk_AntInfo, only : Ant_pos !, Time_dim, Ant_Stations, Ant_IDs, Ant_nr
   use DataConstants, only : Production
   use constants, only : dp,pi,sample,c_mps
   use ThisSource, only : Safety
   use FitParams, only : SpaceCov, Sigma_AntT
   Implicit none
   Integer, Intent(in) :: i_ref, i_ant, i_chunk
   Real(dp), Intent(in) :: SrcPos(1:3)
   Real(dp), Intent(out) :: EEst
   !integer :: i_chunk, i_peak,i_ref, i_eo, i_ant, j_corr
   integer :: i,j
   Real(dp) :: Jac(1:3), Da, Dr, Drel, Er1, MaxTimeWin, IndxRefrac, small=1.d-9, Large=1.d11![samples]
   Real(dp), external :: RefracIndex
   !  small=epsilon(small)
   !
   Dr=0. ; Da=0. ; Drel=0.
   Er1=Sigma_AntT   ! default accuracy of pulse-peak
   Do i=1,3
      If(SpaceCov(i,i).lt.0.) Then
         Er1=Er1+Safety
      EndIf
      Dr=Dr+ (SrcPos(i)-Ant_pos(i,i_ref,i_chunk))*(SrcPos(i)-Ant_pos(i,i_ref,i_chunk))    ! True distance from source to reference antena
      Da=Da+ (SrcPos(i)-Ant_pos(i,i_ant,i_chunk))*(SrcPos(i)-Ant_pos(i,i_ant,i_chunk))  ! True distance from source to measurement antena
      Drel=Drel+ (Ant_pos(i,i_ref,i_chunk)-Ant_pos(i,i_ant,i_chunk))*(Ant_pos(i,i_ref,i_chunk)-Ant_pos(i,i_ant,i_chunk))  ! relative distance to reference antena
   Enddo
   Da=sqrt(Da)
   Dr=sqrt(Dr)
   IndxRefrac = RefracIndex(SrcPos(3))
   MaxTimeWin= sqrt(Drel)*IndxRefrac/(c_mps*sample) ! In [samples]
   Do i=1,3
      Jac(i)=((SrcPos(i)-Ant_pos(i,i_ant,i_chunk))/Da -(SrcPos(i)-Ant_pos(i,i_ref,i_chunk))/Dr) &
      *IndxRefrac/(c_mps*sample)
   EndDo ! the space components of the Jacobian; F in the notes as given in \eqref{RealKalmanJacobian}
   !write(2,*) 'Jac:',Jac, SigmaSpace
   EEst=Er1*Er1  !  Estimated error, sum of error in estimate (\eqref{RealKalmanEstimateError} in notes) and intrinsic measurement error
   Do i=1,3
   Do j=1,3
      EEst=EEst + Jac(i)*SpaceCov(i,j)*Jac(j)
   Enddo ; enddo  ! calculated the denominator of \eqref{RealKalmanWeight} of the notes
   !If(.not. Production)
   !If(EEst.lt. 1.) then
   !   !write(2,*) i_ref,i_ant,'SearchWin:',sqrt(EEst),MaxTimeWin,Da,Dr,Er1
   !   write(2,*) 'SearchWin:',SpaceCov(1,1),SpaceCov(2,2),SpaceCov(3,3),Jac
   !   write(2,*) i_ref,i_ant,'SearchWin:',sqrt(EEst),MaxTimeWin,Da,Dr,Er1,Safety, Sigma_AntT
   !EndIf
   !EEst=sqrt(EEst)  ! Estimate of the error in [samples]
   !write(2,*) 'SearchWin:', EEst,MaxTimeWin, da, dr, Er1, Jac(1:3)
   !Flush(unit=2)
   EEst=MIN(sqrt(EEst),MaxTimeWin)  ! Reduction of MaxTimeWin due to the fact that the actual window is larger than SearchRange
End Subroutine SearchWin
!========================================
Subroutine CalcCovarianc(ph,D,h,Dv)
   !Cov(i,j)=sum_k[(dx_i/dv_k)*(dx_j/dv_k) Delta(v_k)^] when errors in v_k are independent
   use constants, only : dp
   use FitParams, only :  SpaceCov
   Implicit none
   Real(dp), Intent(in) :: ph,D,h
   Real(dp), Intent(in) :: Dv(1:3)
   !Real(dp), Intent(out) :: Cov(1:3,1:3)
   integer :: i,j,k
   Real(dp) :: dXdv(1:3,1:3), cph,sph
   !  small=epsilon(small)
   !
   ! Polar
   !Loc(1)=R*cos(th)*cos(ph)  ! North
   !Loc(2)=R*cos(th)*sin(ph)  !  East
   !Loc(3)=R*sin(th)     ! Height
   ! v1=ph=azimuth;  v2=th=horizon angle; v3=R=distance (following interferometry notation)
   !cth=cos(th)  ; sth=sin(th)
   !cph=cos(ph)  ; sph=sin(ph)
   !dXdv(1,1)=-R*cth*sph
   !dXdv(1,2)=-R*sth*cph
   !dXdv(1,3)=  cth*cph
   !dXdv(2,1)=+R*cth*cph
   !dXdv(2,2)=-R*sth*sph
   !dXdv(2,3)=  cth*sph
   !dXdv(3,1)=0.
   !dXdv(3,2)= R*cth
   !dXdv(3,3)=  sth
   !
   ! Planar
   !Loc(1)=D*cos(ph)  ! North
   !Loc(2)=D*sin(ph)  !  East
   !Loc(3)=h     ! Height
   ! v1=ph=azimuth;  v2=D=plane distance;  v3=h=height
   !write(2,*) 'ph',ph,D,H,Dv(:)
   cph=cos(ph)  ; sph=sin(ph)
   dXdv(1,1)= D*cph
   dXdv(1,2)= sph
   dXdv(1,3)=  0.
   dXdv(2,1)=-D*sph
   dXdv(2,2)= cph
   dXdv(2,3)=  0.
   dXdv(3,1)=0.
   dXdv(3,2)= 0.
   dXdv(3,3)=  1.
   !
   Do i=1,3
      dXdv(i,:)=dXdv(i,:)*Dv(:)
   Enddo
   SpaceCov(:,:)=0.
   Do i=1,3
   Do j=1,3
      Do k=1,3
         SpaceCov(i,j)=SpaceCov(i,j) + dXdv(i,k)*dXdv(j,k)
      Enddo
   Enddo
   Enddo
   Return
End Subroutine CalcCovarianc
