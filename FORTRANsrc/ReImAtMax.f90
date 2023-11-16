Subroutine ReImAtMax(CrCor, ACCorr, RCCorr, RtMax, Aval, Rval, Ival, Error)
!   Get the value of real and imaginary parts of the X-Correlation at the position of the Max in abs. value
   use Chunk_AntInfo
   use ThisSource
   !use FitParams, only : N_FitPar_max, Fit_AntOffset, Fit_TimeOffsetAnt, Fit_TimeOffsetStat
   use constants, only : dp !,pi,ci,sample
   !use StationMnemonics, only : Station_ID2Mnem, Statn_ID2Mnem, Statn_ID2Mnem
   Implicit none
   Integer, intent(in) ::  i_ant, i_start, j_corr,i_Peak
   Complex(dp), intent(in) ::  CrCor(1:T2_dim)
   Real(dp), intent(out) :: ACCorr(-Safety:Safety), RCCorr(-Safety:Safety)
   Real(dp), intent(out) :: RtMax, Aval, Rval, Ival, Error
   !
   integer :: i_loc(1), t_Max
   Real(dp) :: ACCorr_pp(-Safety:Safety), RCCorr_pp(-Safety:Safety), ICCorr_pp(-Safety:Safety), ICCorr(-Safety:Safety)
   Real(dp), save :: tA,B,Yp, dt
   !
   ACCorr(0:Safety) = abs(CrCor(1:Safety+1))
   ACCorr(-Safety:-1) = abs(CrCor(T2_dim-Safety+1:T2_dim))  ! correct for upsampling
   ICCorr(0:Safety) = Imag(CrCor(1:Safety+1))
   ICCorr(-Safety:-1) = Imag(CrCor(T2_dim-Safety+1:T2_dim))  ! correct for upsampling
   RCCorr(0:Safety) = Real(CrCor(1:Safety+1))
   RCCorr(-Safety:-1) = Real(CrCor(T2_dim-Safety+1:T2_dim))  ! correct for upsampling
   Call spline_cubic_set( 2*Safety+1, t_ccorr(-Safety), ACCorr(-Safety), ACCorr_pp(-Safety) )
   Call spline_cubic_set( 2*Safety+1, t_ccorr(-Safety), ICCorr(-Safety), ICCorr_pp(-Safety) )
   Call spline_cubic_set( 2*Safety+1, t_ccorr(-Safety), RCCorr(-Safety), RCCorr_pp(-Safety) )
   !
   i_loc=MaxLoc(ACCorr(:) )
   t_max=i_loc(1) - Safety -1
   !If(Unique_StatID(i_stat).eq.147) Write(2,*) '147,i,j',i_ant,j_corr,'Corrected Rdist=',Rdist,t_max
   If(count((ExclStatNr(:,i_peak)-Ant_Stations(i_ant,i_start)).eq.0,1) .ge. 1) then
      !write(2,*) 'excluded station=',Ant_Stations(i_ant,i_start),i_ant, ', for i_peak=',i_peak
      RtMax=0.
      Error=2000
   ElseIf((t_max.le. -Safety) .or. (t_max .ge. Safety)) then
      write(2,"(A,I5,I7,A,A5,A,I2)") 'maximum in correlation function at bound',t_max, &
          Ant_Stations(i_ant,i_start),'=',Statn_ID2Mnem(Ant_Stations(i_ant,i_start)),', peak#=',i_Peak
      RtMax=t_max
      Error=200  ! time samples
   else
      Yp=ACCorr(t_max+1)-ACCorr(t_max) - ACCorr_pp(t_max+1)/6. -ACCorr_pp(t_max)/3.    !CCorr_der(t_max,j_corr,i_Peak)
      if(Yp.lt.0.) then
          t_max=t_max-1  !  shift to the point left of the real maximum
          Yp=CCorr_der(t_max,j_corr,i_Peak)
      endif
      tA=CCorr_pp(t_max+1,j_corr,i_Peak) -CCorr_pp(t_max,j_corr,i_Peak)   ! 2 * A
      B=CCorr_pp(t_max,j_corr,i_Peak)
      dt= ( -B - sqrt(B*B - 2.*tA * Yp) )/tA
      RtMax=t_max+dt  ! position of the maximum in the correlation function
      Error= 0.2  ! time samples
   endif
   Call spline_cubic_val( 2*Safety+1, t_ccorr(-Safety), ACCorr(-Safety), ACCorr_pp(-Safety), RtMax, Aval) )
   !Call spline_cubic_set( 2*Safety+1, t_ccorr(-Safety), ICCorr(-Safety), ICCorr_pp(-Safety) )
   Call spline_cubic_val( 2*Safety+1, t_ccorr(-Safety), ICCorr(-Safety), ICCorr_pp(-Safety), RtMax, Ival) !, ypval, yppval )
   Call spline_cubic_val( 2*Safety+1, t_ccorr(-Safety), RCCorr(-Safety), RCCorr_pp(-Safety), RtMax, Rval) !, ypval, yppval )
   !
   Return
End Subroutine ReImAtMax
