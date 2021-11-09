!=====================================
Subroutine SourceTryal(DistMax)
!  Find an initial guess for the source position through a grid search.
!  Works only in imaging mode since covariance matrix is stored in this mode only. For this
!     the link " X(  XIndx(i,i_Peak) ) = SourcePos(FitPos(i),i_Peak)" is not implemented.
   !use Chunk_AntInfo, only : CTime_spectr, Ant_pos, Time_dim, Ant_Stations, Ant_IDs, Ant_nr
   use DataConstants, only : Production
   use ThisSource, only : Nr_Corr, PeakNr, PeakNrTotal , PeakPos  ! , Peak_eo, Peak_start
   use ThisSource, only : CorrAntNrs, CCorr_max, SourcePos
   use FitParams, only : Sigma, SpaceCov, CalcHessian
   use DataConstants, only : Production
   !use FitParams, only : N_FitPar_max, N_FitStatTim, Nr_TimeOffset, PulsPosCore, CalcHessian, Kalman
   use constants, only : dp,pi
   Implicit none
   Real(dp), intent(in) :: DistMax
   Integer, parameter :: NAng= 16
   integer :: i_loc(1),i_chunk, i_peak,i_ref, i_eo, i_ant, j_corr, StatMax=2000
   integer :: i,j,k, N_ActAnt, NTry
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
   !
   NTry=0
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
   Call CheckAntInRange(SourceTrPos,NAng,N_ActAnt)
   If(.not. Production) write(2,*) 'N_ActAnt=',N_ActAnt
   If(N_ActAnt.lt.4) then
      write(2,*) 'too few active antennas',N_ActAnt,' for',D
      If(NTry.gt. 4)  Then
         write(2,*) 'Too many tryals in SourceTryal, formerly: stop SourceTryal2'
         Return
      EndIf    !stop 'SourceTryal2'
      If(WildGuess) then
         D=D/2.
      Else
         dph=dph/2.
      Endif
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
End Subroutine SourceTryal
!=====================================
Subroutine CheckAntInRange(SourceTrPos,N_try,N_ActAnt)
! Checks which antennat in the list "CorrAntNrs(2:j_corr,i_eo, i_chunk)" are able to see pulses from
!     all source locations given in "SourceTrPos(3,1:N_try)".
   !use Chunk_AntInfo, only : CTime_spectr, Ant_pos, Time_dim, Ant_Stations, Ant_IDs, Ant_nr
   use ThisSource, only : Nr_Corr, PeakNr, PeakNrTotal, PeakPos, Peak_eo, Peak_start
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
      i_chunk=Peak_start(i_peak)
      j_corr=1
      i_ref=CorrAntNrs(j_corr,i_eo, i_chunk)
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
   !use ThisSource, only : Nr_Corr, PeakNr, PeakNrTotal, PeakPos, Peak_eo, Peak_start
   !use ThisSource, only : CorrAntNrs, CCorr_max, SourcePos, Safety
   use ThisSource, only : Nr_Corr, PeakNr, PeakNrTotal, PeakPos, Peak_eo, Peak_start, T_Offset, CorrAntNrs
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
   i_chunk=Peak_start(i_peak)
   i_ant=CorrAntNrs(j_corr,Peak_eo(i_peak), i_chunk)
   Do k=1,N_try
      TS(k)=SubRelDist(SourceTrPos(1,k),i_ant,i_chunk) - SubRelDist(SourceTrPos(1,k),i_ref,i_chunk) - T_Offset(j_corr,i_Peak)
      !k=k+1
   Enddo
   !write(2,*) 'TS:',TS
   TS_max=MAXVAL(TS)
   TS_min=MINVAL(TS)
   !   i_loc=MinLoc(CS(0:3*NAng) ) - 1 ! to correct for 0 as first element
End Subroutine AntInRangeChk
!=====================================
Subroutine OptSrcPos(SourceTrPos,N_try,ChiSq_min)
   !use Chunk_AntInfo, only : CTime_spectr, Ant_pos, Time_dim, Ant_Stations, Ant_IDs, Ant_nr
   use DataConstants, only : Production
   use ThisSource, only : Nr_Corr, PeakNrTotal, PeakPos, Peak_eo, Peak_start ! , PeakNr
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
      i_chunk=Peak_start(i_peak)
      j_corr=1
      i_ref=CorrAntNrs(j_corr,i_eo, i_chunk)
    !  R(:)=0.
      !write(2,*) 'Nr_Corr(i_eo,i_chunk):',Nr_Corr(i_eo,i_chunk),i_eo,SourceTrPos(:,1)
      Do k=1,N_try
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
   !If(.not. Production)
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
   use constants, only : dp,pi,sample,c_mps,Refrac
   use ThisSource, only : Safety
   use FitParams, only : SpaceCov, Sigma_AntT
   Implicit none
   Integer, Intent(in) :: i_ref, i_ant, i_chunk
   Real(dp), Intent(in) :: SrcPos(1:3)
   Real(dp), Intent(out) :: EEst
   !integer :: i_chunk, i_peak,i_ref, i_eo, i_ant, j_corr
   integer :: i,j
   Real(dp) :: Jac(1:3), Da, Dr, Drel, Er1, MaxTimeWin, small=1.d-9, Large=1.d11![samples]
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
   MaxTimeWin=sqrt(Drel)*Refrac/(c_mps*sample) ! In [samples]
   Do i=1,3
      Jac(i)=((SrcPos(i)-Ant_pos(i,i_ant,i_chunk))/Da -(SrcPos(i)-Ant_pos(i,i_ref,i_chunk))/Dr) &
      *Refrac/(c_mps*sample)
   EndDo ! the space components of the Jacobian; F in the notes as give in \eqref{RealKalmanJacobian}
   !write(2,*) 'Jac:',Jac, SigmaSpace
   EEst=Er1*Er1  !  Estimated error, sum of error in estimate (\eqref{RealKalmanEstimateError} in notes) and intrinsic measurement error
   Do i=1,3
   Do j=1,3
      EEst=EEst + Jac(i)*SpaceCov(i,j)*Jac(j)
   Enddo ; enddo  ! calculated the denominator of \eqref{RealKalmanWeight} of the notes
   !write(2,*) i_ref,i_ant,'SearchWin:',sqrt(EEst),MaxTimeWin,Da,Dr,Er1
   !EEst=sqrt(EEst)  ! Estimate of the error in [samples]
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
