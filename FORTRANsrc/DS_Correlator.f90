! ----------------------------------------------------------------------------------------
Subroutine ApplyCorrelator(RA, Label, SourcTotNr, Aweight, dD, Dnr, tauMax)
!  In all distributions  an source-pair (i,j) carries a weight AmplWeight=(A_i*Aweight+1.)*(A_j*Aweight+1.)
!     where A denotes the amplitude
!  T_trace(i_t,i_d) is the sum of weights of source pairs at fine-binned relative time and coarse-binned distance
!     separation of (i_t*dT , i_d*TD_dd ) where dT=tauMax/(N_it=200) and TD_dd=dD*Dnr/(4.*N_id=4*5=20).
!
   Use constants, only : dp, CI, pi, c_l
   use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, DataFolder, FlashName, Windows, RunMode
   Use DS_Select, only : Image !,  WriDir
   use GLEplots, only : GLEplotControl
   IMPLICIT none
   Integer, parameter :: N_it=200, N_id=5  ! for the TD-trace plots
   real*8, intent(IN) :: RA(4,*)  ! 1=t [ms]  2-4= E,N,h in [km]
   Integer, intent(in) :: Label(4,*)  ! (2,*) contains intensity
   Integer, intent(IN) :: SourcTotNr  ! number of sources stored in RA, passing the selection criteria
   real*8, intent(IN) :: dD  ! bin width for fine-distance
   real*8, intent(IN) :: tauMax  ! Max range for t-grid for calculating \zeta (coarse t-grid)
   Integer, intent(IN) :: Dnr  ! number of fine-distance points
   !Integer, intent(IN) :: Tnr  ! number of fine-time points range will be full range; not used, set equal to N_it
   Integer :: i, i_src, j_src, m_src, n_src
   real*8 :: TD_corr(0:Dnr), T2D_corr(0:Dnr), T4D_corr(0:Dnr), nD_dens(0:Dnr), TD_var(0:Dnr), T2D_var(0:Dnr), nrm, TD_dd
   real*8 :: T_span, D_span, dT, t_ij, t_mn, d_ij, d_mn
   Integer :: nD_corr(0:Dnr), T_trace(0:N_it,0:N_id), i_t, i_d,i_m
   Real*8 :: TDCC(0:N_it,0:N_id), t_trace_CC(0:N_it,0:N_id), W_d, Omg_d, Omg_t, W_t
   Integer :: n_corr_t(0:N_it), i_Ct
   Real(dp) :: D_corr_t(0:N_it), n_dens_t(0:N_it), t_range_Ct, dt_Ct, GrandNorm
   !
   t_range_Ct=5.  ! [ms]
   dt_Ct=t_range_Ct/N_it  ! [ms]
   nD_corr(0:Dnr)=0
   TD_corr(0:Dnr)=0.
   T2D_corr(0:Dnr)=0.
   T4D_corr(0:Dnr)=0.
   D_corr_t(:)=0.
   n_corr_t(:)=0
   T_span=RA(1,SourcTotNr)-RA(1,1)+1.D-9  ! for rough t-binning
   If(tauMax.gt. 0.) Then
      T_span=tauMax
   EndIf
   dT=T_span/N_it  ! for  t-binning
   D_span=dD*Dnr  ! for rough d-binning
   write(2,*) T_span, dT, D_span
   T_trace(:,:)=0
   Omg_t=dT
   Omg_d=dD
   TDCC(:,:)=0.
   TD_dd = D_span/(4.*N_id)  ! d_bin size for t-trace plots; factor 4 to get better resolution for the smaller distances
   !
   GrandNorm=0.
   Do i_src=1,SourcTotNr-1
      Do j_src=i_src+1,SourcTotNr
         AmplWeight=(Label(2,i_src)*Aweight+1.)*(Label(2,j_src)*Aweight+1.)  ! Amplitude determined weight of the new source
         GrandNorm=GrandNorm + AmplWeight   ! independent of later grid-cuts
         !
         ! Calculate <\tau^k>(d) for narrowly binned distances where \tau=t_ij
         ! dist: bin size=dD with Dnr bins (subroutine input); fine binning
         d_ij=sqrt( SUM((RA(2:4,i_src)-RA(2:4,j_src))**2) )
         i=FLOOR(d_ij/dD)
         If(i.gt.Dnr) i=Dnr
         t_ij=ABS(RA(1,j_src)-RA(1,i_src))  ! Is already positive because of the presorting in RA
         TD_corr(i)=TD_corr(i) + AmplWeight*t_ij
         T2D_corr(i)=T2D_corr(i) + AmplWeight*t_ij**2
         T4D_corr(i)=T4D_corr(i) + AmplWeight*t_ij**4
         nD_corr(i)=nD_corr(i) + AmplWeight
         !
         !  Calculate T_trace(i_t,i_d) equal to zeta function notes (=product of delta functions in distance and time, binned)
         !  time: bin size=dT=tauMax/(N_it=200) with (N_it=200) bins; fine binning
         !  dist: bin size=TD_dd=dD*Dnr/(4.*N_id=20) with (N_id=5) bins;  coarse binning
         i_t=FLOOR(t_ij/dT)
         If(i_t.gt.N_it) i_t=N_it
         i_d=Floor( d_ij/TD_dd )
         If(i_d.gt.N_id) Then
            If(i_t.gt.N_it) cycle
            i_d=N_id
         EndIf
         T_trace(i_t,i_d)=T_trace(i_t,i_d)+ AmplWeight
         !
         ! Calculate <d_ij>(t) for narrowly binned time differences
         !  time: bin size=dt_Ct=(t_range_Ct=5)/(N_it=200) with (N_it=200) bins; fine binning
         i_Ct=FLOOR(t_ij/dt_Ct)
         If(i_Ct .gt. N_it) i_Ct=N_it
         n_corr_t(i_Ct) = n_corr_t(i_Ct) + AmplWeight
         D_corr_t(i_Ct) = D_corr_t(i_Ct) + AmplWeight*d_ij
         !
         If(SourcTotNr.lt.-200) Then
         !     Disabled    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !  Cross correlation between the zeta distribution and itself using gaussian smoothing (width) in distance (Omg_d) and time (Omg_t)
         !  time: bin size=dT=tauMax/(N_it=200) with (N_it=200) bins; fine binning.
         !  dist: bin size=2.*(D_span=dD*Dnr)/(N_id=5) with (N_id=5) bins;  coarse binning.
         !        Seems not very useful, Array-limit checking must be done.
         !        Later repeated, replacing gaussian distance weighting with coarse distance binning.
            Do m_src=1,SourcTotNr-1
               Do n_src=i_src+1,SourcTotNr
                  d_mn=sqrt( SUM((RA(2:4,m_src)-RA(2:4,n_src))**2) )
                  If(abs(d_ij-d_mn) .gt. 3.*Omg_d) cycle
                  t_mn=ABS(RA(1,m_src)-RA(1,n_src))
                  i_d=Floor( N_id*(d_ij+d_mn)/(2.*D_span) )
                  w_d=exp(-((d_ij-d_mn)/Omg_d)**2)
                  Call AddTDCC(t_ij, t_mn, Omg_t, dt, w_d, TDCC(0,i_d), N_it)
                  Call AddTDCC(t_ij,-t_mn, Omg_t, dt, w_d, TDCC(0,i_d), N_it)
               EndDo ! n_src=i_src+1,SourcTotNr
            EndDo ! m_src=1,SourcTotNr-1
         EndIf
      EndDo ! j_src=i_src+1,SourcTotNr
      !If(SourcTotNr.lt.-200) write(2,*) 'ApplyCorrelator: i_src=', i_src, TDCC(:,1)
      !flush(unit=2)
   EndDo ! i_src=1,SourcTotNr-1
   !
   write(2,*) 'tD_correlations, distance step=',dD
   !
   !Normalize
   !nT_corr(0:N_it)=sum(T_trace(0:N_it,:))
   D_corr_t(:)=D_corr_t(:)/(n_corr_t(:)+0.00001)
   n_dens_t(:)=n_corr_t(:)/(dt_Ct*GrandNorm)
   TD_corr(:)=TD_corr(:)/(nD_corr(:)+0.001)
   T2D_corr(:)=T2D_corr(:)/(nD_corr(:)+0.001)
   T4D_corr(:)=T4D_corr(:)/(nD_corr(:)+0.001)
   TD_var(:)=sqrt(T2D_corr(:)-TD_corr(:)**2)
   T2D_var(:)=sqrt(T4D_corr(:)-T2D_corr(:)**2)
   T2D_corr(:)=sqrt(T2D_corr(:))
   nD_dens(:)=nD_corr(:)/(dD*GrandNorm)
   T_trace(:,:)=T_trace(:,:) /(dt*dD*GrandNorm)  ! same as the binned \zeta density from the notes
   !
   OPEN(unit=29,FILE=TRIM(DataFolder)//trim(Image)//'TDcorr.dat',FORM='FORMATTED',STATUS='unknown')  ! will contain all accepted sources
   write(29,"(2x,A,1x,I6,1x,F6.3,i5,' !')") trim(Image), SourcTotNr, dD, Dnr
   Do i=0,Dnr
      nrm=nD_corr(i)+0.001
      !write(2,*) (i+0.5)*dD,nD_corr(i),TD_corr(i), T2D_corr(i), TD_var(i), sqrt(T2D_var(i))
      write(29,*) (i+0.5)*dD,nD_dens(i),TD_corr(i), T2D_corr(i), TD_var(i), sqrt(T2D_var(i))
   EndDo ! i=0,Dnr
   close(unit=29)
   Write(2,*) 'File written:',TRIM(DataFolder)//trim(Image)//'TDcorr.dat',' with ',Dnr, ' lines'
   !
   OPEN(unit=29,FILE=TRIM(DataFolder)//trim(Image)//'TDtrace.dat',FORM='FORMATTED',STATUS='unknown')  ! will contain all accepted sources
   write(29,"(2x,A,1x,F8.3,1x,F6.3,i5,' !')") trim(Image), T_span, TD_dd, N_id
   Do i_t=0,N_it-1
      !Write(2,*) i_t,T_trace(i_t,:)
      Write(29,*) (i_t+0.5)*dt_Ct, n_dens_t(i_t), D_corr_t(i_t), (i_t+0.5)*dT, T_trace(i_t,:)/(dt*dD*(SourcTotNr-1)*SourcTotNr)
   EndDo ! i_t=0,N_it
   close(unit=29)
   Write(2,*) 'File written:',TRIM(DataFolder)//trim(Image)//'TDtrace.dat',' with ',N_it, ' lines'
   !
   !  Time cross-correlation between the zeta distribution (=T_trace(i_t,i_d)) and itself using gaussian smoothing in time (width=Omg_t)
   !  time: bin size=dT=tauMax/(N_it=200) with (N_it=200) bins; fine binning
   !  dist: bin size=TD_dd=dD*Dnr/(4.*N_id=20) with (N_id=5) bins;  coarse binning
   Omg_t=dT  ! smoothing window should be order of the bin-width for the trace
   t_trace_CC(:,:)=0
   Do i_d=0, N_id
      W_t=0.
      !W_t=SUM(T_trace(:,i_d))
      !write(2,*) 'i_d:', i_d, W_t*W_t/((dt*dD*(SourcTotNr-1)*SourcTotNr)**2)
      !W_t=W_t/(N_it+1.)
      Do i_t= 0, N_it
         Do i_m= 0, N_it
            w_d=(T_trace(i_t,i_d)-W_t)*(T_trace(i_m,i_d)-W_t)
            Call AddTDCC(i_t*dt, i_m*dt, Omg_t, dt, w_d, t_trace_CC(0,i_d), N_it)
            Call AddTDCC(i_t*dt, -i_m*dt, Omg_t, dt, w_d, t_trace_CC(0,i_d), N_it)
            Call AddTDCC(-i_t*dt, i_m*dt, Omg_t, dt, w_d, t_trace_CC(0,i_d), N_it)
            Call AddTDCC(-i_t*dt, -i_m*dt, Omg_t, dt, w_d, t_trace_CC(0,i_d), N_it)
         EndDo
      EndDo
      !W_t=SUM(T_trace(:,i_d))
      !write(2,*) 'i_d:', i_d, W_t*W_t/((dt*dD*(SourcTotNr-1)*SourcTotNr)**2)
      !t_trace_CC(:,i_d)=(t_trace_CC(:,i_d)-2*W_t*W_t/(N_it+1.))/((dt*dD*(SourcTotNr-1)*SourcTotNr)**2)
      t_trace_CC(:,i_d)=t_trace_CC(:,i_d)*dD
      W_t = SUM(t_trace_CC(:,i_d))
      write(2,*) 'SUM(t_trace_CC(:,i_d)):', W_t
      !t_trace_CC(:,i_d)=t_trace_CC(:,i_d) - w_d
   EndDo
   !
   OPEN(unit=29,FILE=TRIM(DataFolder)//trim(Image)//'TtrCC.dat',FORM='FORMATTED',STATUS='unknown')  ! will contain all accepted sources
   Do i_t=0,N_it
      !If(SourcTotNr.lt.-200) Write(2,*) i_t,TDCC(i_t,:)
      !Write(2,*) i_t,t_trace_CC(i_t,:)
      Write(29,*) i_t*dT, T_trace_CC(i_t,:)
   EndDo ! i_t=0,N_it
   close(unit=29)
   Write(2,*) 'File written:',TRIM(DataFolder)//trim(Image)//'TtrCC.dat',' with ',N_it+1, ' lines'
   !
   Call GLEplotControl(PlotType='TD_corrPlot', PlotName=TRIM(Image)//'TDcorr', &
         PlotDataFile=TRIM(DataFolder)//trim(Image), Submit=.false.)
End Subroutine ApplyCorrelator
!-----------------------------------------------------
Subroutine AddTDCC(t_ij,t_mn, Omg_t, dt, w_d, TDCC, tCC_dim)
   Use constants, only : dp, CI, pi, c_l
   IMPLICIT none
   real(dp), intent(IN) :: t_ij, t_mn, Omg_t, dt, w_d
   Integer, intent(IN) :: tCC_dim
   Real(dp), intent(INOUT) :: TDCC(0:tCC_dim)
   Real(dp) :: Del_t, W_t, t_cc
   Integer :: it_low, it_upp, i_t,i
   Real(dp), save :: OmgT_old=-1.
   Real(dp), save :: Wt_i(0:100), d_Wt
   If(OmgT_old.ne.Omg_t) Then ! make template for gaussian folding, store in wt_i
      d_Wt=3./99.
      Do i_t=0,100
         t_cc=i_t*d_Wt ! to cover the full range for 3 sigma in close to 100 intervals
         Wt_i(i_t)=exp(-(t_cc)**2)
      EndDo ! i_t=it_low, it_upp
      OmgT_old = Omg_t
      W_t=(Wt_i(0)+2.*SUM( Wt_i(1:100) ))*d_Wt
      Wt_i(:)=Wt_i(:)/W_t
      write(2,*) 'Wt_i', Wt_i(:)
   EndIf
   Del_t=(t_ij-t_mn)  !   t_ij >0   exp(-(T_cc-(t_ij-t_mn))**2/omg**2)
   it_low=CEILING((-3*Omg_t+Del_t)/dt)
   If(it_low.lt.0) it_low=0
   If(it_low.gt.tCC_dim) return
   it_upp=FLOOR((3*Omg_t+Del_t)/dt)
   If(it_upp.lt.0) return
   If(it_upp.gt.tCC_dim) it_upp=tCC_dim
   Do i_t=it_low, it_upp
      t_cc=dt*i_t
      i=NINT( ABS((t_cc-Del_t)/(d_Wt*Omg_t)) )
      !W_t=exp(-((t_cc-Del_t)/Omg_t)**2)
      W_t=Wt_i(i)*dt/Omg_t
      TDCC(i_t)=TDCC(i_t) + W_d*W_t
   EndDo ! i_t=it_low, it_upp
   !write(2,*) 'it_low, it_upp', t_ij,t_mn, it_low, it_upp, TDCC(it_low: it_upp)
   Return
End Subroutine AddTDCC
! ----------------------------------------------------------------------------------------
