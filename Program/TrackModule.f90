    Include 'ConstantsModules.f90'
!-----------------------------------
Module Tracks
!   Use Tracks, only : PreDefTrack(0:3,tNr,1:PreDefTrackNr), tNrS_max(i_Pre), PreDefTrackNr
   Use constants, only : dp, CI, pi, c_l
   !
   Character(len=100) :: PreDefTrackFile
   Character*130 :: datfile
   Character*100, save :: WriDir
   Integer, parameter :: PreDefTrackNr_Max=10  !  Maximum number of Pre-defined Tracks
   Integer, parameter :: NrPreDefTrackPoints_Max=900   ! Maximum number of defining times per Pre-defined Track
   Real(dp), save :: PreDefTrack(0:3,1:NrPreDefTrackPoints_Max,1:PreDefTrackNr_Max)
   Integer, save :: tNrS_max(1:NrPreDefTrackPoints_Max), PreDefTrackNr
   Integer, save :: tNrS(1:PreDefTrackNr_Max)=1
   Integer :: TrackNr, LongTrackNr, LongTrack_Min=20 !12
   Integer, parameter :: TrackNrMax=19, TrackLenMax=5000
   Integer :: TrackENr(TrackLenMax), TrackE(TrackLenMax,TrackNrMax)
   Real*8 :: LeaderPos(1:4,0:TrackLenMax,TrackNrMax)
   Real*8 :: Wtr  ! weighting of new source location for track definition
   Real*8 :: MaxTrackDist ! maximal distance from trackhead
   Real*8 :: TimeWin   !  Time window for track smoothing
   Real*8 :: HeightFact   !  Factor to increase "MaxTrackDist" for height, used for track smoothing
   Real*8 :: TrackTimeLim=5. ! maximal distance between sources to mark short tracks for deletion
   Integer :: TrackNrLim=1
   Integer :: TOrder=-1  ! set to -1 for running tracks backward
   Integer :: NLongTracksMax
   Real*8 :: MeanTrackLoc(1:3,500), dt_MTL
   Integer :: Nmax_MTL=500
   !
   integer, parameter :: maxd=123000
   real*8 :: RA(4,maxd)  ! 1=t [ms]  2-4= E,N,h in [km]
   Integer :: Label(4,maxd)  ! Unique source label and amplitude and widths
   real :: QualIndic(4,maxd)  ! Quality indicators
   Integer :: EventNr  ! number of sources stored in RA, passing the selection criteria
   !
   Integer, parameter :: Nmax_Ampl=5000
   Real(dp) :: MaxAmplFitPercent=100.
   Integer :: FitFunc=0  ! Fitfunftion=0 : exp ; =1 powerlaw
   Real(dp) :: Ampl_Hist(Nmax_Ampl)  ! Amplitude histogram filled out in "AmplitudeFitTracks"
   Real(dp) :: Ampl_Hist_W(Nmax_Ampl)  ! Amplitude histogram filled out in "AmplitudeFitTracks"
   Real(dp) :: Ampl(Nmax_Ampl)  ! Amplitude histogram filled out in "AmplitudeFitTracks"
   Integer :: i_AmplTh  ! Calculated in "AmplitudeFitTracks" to the index of the first amplitude for which the probability spectrum is fitted
   Real(dp) :: Nrm, a1,b1,c1,d1,ChiSq1,Nrm1, a2,b2,c2,d2,ChiSq2,Nrm2, a3,b3,c3,d3,ChiSq3 ! parameters for describing the amplitude distribution
   Real*8 :: AmplScale=1./100., d_AmplScale ! Scaling factor for amplitude to convert from integer to unity-normalized real
contains
SUBROUTINE HPSORT_mult_RI(RA,IA,N)
!*****************************************************
!*  Sorts an array RA of length N in ascending order *
!*                by the Heapsort method             *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*	    N	  size of table RA                                                                  *
!*          RA(1,1:N)	  table to be sorted, real values                                    *
!*          RA(*,1:N) [real] and IA(*,1:N) [integer]	 rearranged like  RA(1,1:N)            *
!* OUTPUT:                                                                                   *
!*	    RA    table sorted in ascending order                                                 *
!*	    IA    table rearranged following order RA                                             *
!*                                                   *
!* NOTE: The Heapsort method is a N Log2 N routine,  *
!*       and can be used for very large arrays.      *
!*****************************************************
  IMPLICIT none
  integer, intent(in) :: N
  real*8 RA(:,:)
  Integer IA(:,:)
  real*8 RRA(SIZE(RA,1))
  Integer RIA(SIZE(IA,1))
  integer :: d1,d2,L,IR,I,J
  d1=SIZE(RA,1)
  d2=SIZE(IA,1)
  L=N/2+1
  IR=N
  !The index L will be decremented from its initial value during the
  !"hiring" (heap creation) phase. Once it reaches 1, the index IR
  !will be decremented from its initial value down to 1 during the
  !"retirement-and-promotion" (heap selection) phase.
10 continue
  if(L > 1)then
    L=L-1
    RRA(1:d1)=RA(1:d1,L)
    RIA(1:d2)=IA(1:d2,L)
  else
    RRA(1:d1)  =RA(1:d1,IR)
    RA(1:d1,IR)=RA(1:d1,1)
    RIA(1:d2)  =IA(1:d2,IR)
    IA(1:d2,IR)=IA(1:d2,1)
    IR=IR-1
    if(IR.eq.1)then
      RA(1:d1,1)=RRA(1:d1)
      IA(:,1)=RIA(:)
      return
    end if
  end if
  I=L
  J=L+L
20 if(J.le.IR)then
  if(J < IR)then
    if(RA(1,J) < RA(1,J+1))  J=J+1
  end if
  if(RRA(1) < RA(1,J))then
    RA(1:d1,I)=RA(1:d1,J)
    IA(:,I)=IA(:,J)
    I=J; J=J+J
  else
    J=IR+1
  end if
  goto 20
  end if
  RA(1:d1,I)=RRA(1:d1)
  IA(:,I)=RIA(:)
  goto 10
END SUBROUTINE HPSORT_mult_RI
!-----------------------
  Subroutine Selection_sort(a)  ! same as in MappingUtilities.f90
  ! Smallest value first
  ! From:  https://rosettacode.org/wiki/Sorting_algorithms/Selection_sort#Fortran
    INTEGER, INTENT(IN OUT) :: a(:)
    INTEGER :: i, minIndex, temp

    DO i = 1, SIZE(a)-1
       minIndex = MINLOC(a(i:), 1) + i - 1
       IF (a(i) > a(minIndex)) THEN
          temp = a(i)
          a(i) = a(minIndex)
          a(minIndex) = temp
       End IF
    End DO
  End Subroutine Selection_sort
!=========================================
Subroutine AmplitudeFit(NrSources, SourceNrs, SelFileName)
   ! Make Source amplitude histogram and fit this with exp or power law
   ! d_Ampl is stepsize in amplitude, integer
   !      write(29,"(1x,i8,4(2x,g14.8),3x,F7.2,2x,I3,2x,I3)")  Label(1,j),RA(1:4,j), Label(2,j)/100., Label(3:4,j)
   ! Writes files  trim(datfile)//??//'.dat'
   ! When SourceNrs(i) present it contains the list of source numbers that are to be considered
   ! NrSources: number to be considered
   !Use Tracks, only :  Label
   ! contains amplitudes as integers
   !Use Tracks, only : Nmax_Ampl, AmplScale, d_AmplScale, MaxAmplFitPercent
   ! Nmax_Ampl: Maximum number of bins in histogram, determines max amplitude that can be fitted.
   ! AmplScale: factor (>1) used to convert amplitudes to integers as needed for Label(2,:)
   ! d_AmplScale: bin-size for histogramming the scaled amplitudes before fitting
   ! MaxAmplFitPercent: Maximum percentile of sources in overflow bin
   !Use Tracks, only : a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3
   ! extracted parameters of designated fit functions
   ! Presently: 1== (modified) exponential; 2== (modified)powerlaw; 3=not used
   !
   ! internal AmplitudeFit parameters =======================
   !Use Tracks, only : Selection_sort, i_AmplTh, Ampl_Hist, FitFunc
   ! i_AmplTh: Index of first bin in scaled amplitude histogram that is fitted.
   ! Ampl_Hist: actual histogram of amplitudes
   ! FitFunc: Function that is used
   IMPLICIT none
   Integer, intent(in) :: NrSources
   Integer, intent(in), optional :: SourceNrs(1:NrSources)
   character*100, intent(in), optional :: SelFileName
   Integer :: SourceAmp(1:NrSources)  !  Amplitudes of sources in the i^th track
   !
   Integer :: i,k, i_Ampl, i_AmplMax, i_cut, i_MaxFit,N_bins
   Integer ( kind = 4 ) :: meqn, nvar  ! Number of data points to fit & number of parameters
   Character*8 :: extension
   Real*8 :: Max_Ampl, Ampl10, X(4),ChiSQDF, AmplBin, Counts2BinValue, ZeroErr=1., HistNrm ! 0.25
   Real*8 :: Fit_exp(Nmax_Ampl), Fit_pow(Nmax_Ampl),a,b,c
   Real(dp), save :: Tiny=epsilon(Tiny)
   !
   !write(extension,"(A2,i1,A4)") '_s',nxx,'.dat' !,&
   !write(extension,"(i1,A4)") i,'.dat' !,&
   !OPEN(unit=27,FILE='files/'//trim(datfile)//'_??_'//trim(extension),FORM='FORMATTED',STATUS='unknown')
   !write(27,"('! nr', 4(1x,A14),3x,3(1x,A12),2x,A10)") 't [ms]','E [km]','N','h','v_E','v_N','v_h','v [10^6m/s]'
   !
   If(present(SourceNrs)) Then
      Do i=1,NrSources
         SourceAmp(i)=Label(2,SourceNrs(i))
      Enddo
   Else
      SourceAmp(1:NrSources)=Label(2,1:NrSources)
   EndIf
   CALL Selection_sort(SourceAmp)  ! Re-order  sources according to pulse amplitude
   !
   Fit_exp(:)=0.0
   Fit_pow(:)=0.0
   N_bins=-1
   !
   If(MaxAmplFitPercent.lt.0.) MaxAmplFitPercent=0.1
   Max_Ampl=SourceAmp(NrSources)*AmplScale
   Ampl10=SourceAmp(NrSources*9/10)*AmplScale  ! 10 percentile amplitude
   write(2,"(A,I8,A,F9.1,A,F8.2,A,F8.2,A,I4)") 'Amplitudes (#=',NrSources,'), max@', Max_Ampl, &
      ', 5-pctile@', SourceAmp(NrSources*95/100)*AmplScale, ', 10-pctile@', Ampl10
   i_cut=NrSources*(1.-MaxAmplFitPercent/100)
   If(i_cut.ge.0) then
      If(i_cut.eq.0) i_cut=1
      write(2,"(A,F7.3,A,F8.2,A,I4,F8.2)") 'Amplitude with ', MaxAmplFitPercent, &
         'pctile as included in fit @', SourceAmp(i_cut)*AmplScale,', i_cut=', i_cut, SourceAmp(1)*AmplScale
   Else
      i_cut=NrSources  !  No sources in fit
   Endif
   !
   d_AmplScale=1.1  ! factor for log bins
   ! Automatic scaling
   If(NrSources.gt.50) Then
      N_bins=Nmax_Ampl-2  ! -2 for a stable fit with analytic function
      If(NrSources.lt.N_bins) N_bins=NrSources/10
      d_AmplScale=log(1.d0*SourceAmp(NrSources)/SourceAmp(1))/N_bins
      d_AmplScale=exp(d_AmplScale)
      write(2,*) 'Automatic amplitude histogram scaling;', d_AmplScale, N_bins, &
            SourceAmp(1)*AmplScale, SourceAmp(NrSources)*AmplScale
      !write(2,*) ';', d_AmplScale**N_bins, 1.* SourceAmp(NrSources)/SourceAmp(1)
   EndIf
   !
   AmplBin=SourceAmp(1)*d_AmplScale  ! scaled up by factor 1/AmplScale, always start in the first bin
   i_Ampl=1
   i_MaxFit=0
   Ampl_Hist(:)=0
   HistNrm=0.  ! =\int I N(I) dI
   !Ampl10=0.
   Do k=1,NrSources
      Do While( SourceAmp(k) .gt. AmplBin)  ! start from the smallest amplitude and fill up bins
         Counts2BinValue=(AmplScale*AmplBin*(d_AmplScale-1.))
         Ampl_Hist_W(i_Ampl)= Counts2BinValue/sqrt(Ampl_Hist(i_Ampl)+ZeroErr)   ! Calculate fitting weight for previous bin
         Ampl(i_Ampl)=AmplScale*AmplBin*(1.+d_AmplScale)/2.  ! In true units
         Ampl_Hist(i_Ampl)=Ampl_Hist(i_Ampl)/Counts2BinValue
         HistNrm=HistNrm + Ampl_Hist(i_Ampl)*Ampl(i_Ampl)*AmplScale*AmplBin*(d_AmplScale-1.)
         !Ampl10=Ampl10 + Ampl_Hist(i_Ampl)*AmplScale*AmplBin*(d_AmplScale-1.)  ! equals tot nr of sources, as it should
         ! chi^2 = sum_i (Ampl_Hist(i) - model) * Ampl_Hist_W(i)
         !write(2,*) k, i_Ampl, AmplBin,Ampl(i_Ampl), Ampl_Hist(i_Ampl)*Counts2BinValue, Ampl_Hist_W(i_Ampl)/Counts2BinValue
         AmplBin=AmplBin*d_AmplScale
         If(i_cut .gt. k) i_MaxFit=i_Ampl+1   !  Still include this bin in fitting
         i_Ampl=i_Ampl+1  ! start to fill next bin
         If(i_Ampl .ge. Nmax_Ampl) goto 9
      Enddo
      Ampl_Hist(i_Ampl)=Ampl_Hist(i_Ampl)+1.
   EndDo
9  Continue
   !write(2,*) NrSources, SourceAmp(NrSources) , AmplBin, Ampl_Hist(i_Ampl), i_Ampl, Ampl10, HistNrm
   Counts2BinValue=(AmplScale*AmplBin*(d_AmplScale-1.))
   Ampl_Hist_W(i_Ampl)= Counts2BinValue/sqrt(Ampl_Hist(i_Ampl)+ZeroErr)   ! Calculate fitting weight for previous bin
   Ampl(i_Ampl)=AmplScale*AmplBin*(1.+d_AmplScale)/2.  ! In true units
   Ampl_Hist(i_Ampl)=Ampl_Hist(i_Ampl)/Counts2BinValue
   !write(2,*) N_bins, i_Ampl, Ampl(i_Ampl), Ampl_Hist(i_Ampl), Ampl_Hist_W(i_Ampl)
   i_AmplMax=i_Ampl
   If(i_Ampl .ge. Nmax_Ampl) Then
      Ampl_Hist(i_Ampl)=Ampl_Hist(i_Ampl)+NrSources-k  !  Fill out the highest bin with remainder, when necessary
      i_AmplMax=Nmax_Ampl-1
      write(2,*) 'i_Ampl .ge. Nmax_Ampl', i_Ampl,Nmax_Ampl, AmplBin,Ampl(i_Ampl)
   EndIf
   If(i_MaxFit.eq.0) i_MaxFit=i_AmplMax
   !write(2,*) i_AmplMax, NrSources, ', i_MaxFit:',i_MaxFit
   !write(2,*) Ampl(i_AmplMax), SourceAmp(NrSources), AmplBin
   !
   flush(unit=2)
   !
   !write(2,*) 'Ampl_Hist',Ampl_Hist(1:20)
   !write(2,*) 'bins of max etc:',Max_Ampl, i_AmplZero, i_AmplMax, d_Ampl,Ampl_Hist(11)
   !write(2,*) d_Ampl,AmplScale,d_AmplScale
   a1=0.; b1=0.; c1=0.
   !goto 5
   !
   i_AmplTh=1 ! first index for fitting in histogram
   Fit_exp(1:i_AmplTh)=0.
   Meqn = i_MaxFit-i_AmplTh+1      ! number of equations
   !FitFunc=0  ! b*exp(-a*Ampl)  a=x(1)/d_Ampl
   !FitFunc=2 ; nvar= 3 ! b*exp(-a*A-c/A)  a=x(1)/d_Ampl
   FitFunc=3 ; nvar= 3 ! b*exp(-a*A-c/A^2)  a=x(1)/d_Ampl
   x(1)=1/20.   ! CalcN = exp(-X_p(1)*Ampl+X_p(2))
   X(2)=log(Ampl_Hist(i_AmplTh))+ X(1)*ampl(i_AmplTh)
   X(3)=0.
   Call FitAmpl(Meqn, nvar, X,ChiSq1, Fit_exp(i_AmplTh) )
   a1=X(1)
   b1=(x(2))
   c1=x(3)
   d1=0.
   If(Meqn.lt.0) Then ! no fitting possible
      write(2,*) 'no fitting possible of a1--d1'
      !flush(unit=2)
      goto 5
   EndIf
   b1=exp(x(2))
   !goto 8
   If(c1.lt.0.) then
      i_AmplTh=4 ! first index for fitting in histogram
      Fit_exp(1:i_AmplTh)=0.
      Meqn = i_MaxFit-i_AmplTh+1      ! number of equations
      FitFunc=0 ; nvar= 2 ! exp(b-a*A)  a=x(1)/d_Ampl
      x(1)=1/20.   ! CalcN = exp(-X_p(1)*Ampl+X_p(2))
      X(2)=log(Ampl_Hist(i_AmplTh))+ X(1)*ampl(i_AmplTh)
      X(3)=0.
      Call FitAmpl(Meqn, nvar, X,ChiSq1, Fit_exp(i_AmplTh) )
      a1=X(1)
      b1=exp(x(2))
      c1=0
   Endif
   !Goto 1 ! --------------------
   !
5   FitFunc=5 ; nvar= 3  ! F(A)=b (A)^(-a)*exp(-c/A)
   !In GLE: d8 = d_bin*b2*(x)^(-a2)*exp(-c2/x-d2/x^2) with x is the amplitude
   i_AmplTh=1 ! first index for fitting in histogram
   Fit_pow(1:i_AmplTh)=0.
   Meqn = i_MaxFit-i_AmplTh+1      ! number of equations
   X(1)=2.
   X(2)=Ampl_Hist(i_AmplTh) * (ampl(i_AmplTh)**X(1))
   X(3)=0.
   Call FitAmpl(Meqn, nvar, X,ChiSq2, Fit_pow(i_AmplTh) )
   a2=X(1)
   b2=x(2)!*d_Ampl**(X(1)-1)
   c2=x(3)!*d_Ampl
   d2=0.
   If(Meqn.lt.0) Then ! no fitting possible
      write(2,*) 'no fitting possible of a2--d2'
      !flush(unit=2)
      goto 8
   EndIf
   If(c2.lt.0.) then
1     FitFunc=1 ; nvar=2 ! F(A)=b (A)^(-a) with   where A is true amplitude
      !      CalcN = -X(2) * (i_Ampl)**(-X(1))
      !In GLE: d8  = d_bin*b2 * (x)^(-a2)        with x is the amplitude
      i_AmplTh=4  ! first index for fitting in histogram
      Fit_pow(1:i_AmplTh)=0.
      Meqn = i_MaxFit-i_AmplTh+1      ! number of equations
      X(1)=4.
      X(2)=Ampl_Hist(i_AmplTh) * (ampl(i_AmplTh)**X(1))
      x(3)=0.
      !write(2,*) ' initial guess: ',X(1),X(2), Meqn
      Call FitAmpl(Meqn, nvar, X,ChiSq2, Fit_pow(i_AmplTh) )
      a2=X(1)
      b2=x(2)
      c2=0.
      d2=0.
      !write(2,*) X(1),X(2), a2, b2
   EndIf
   !
8  continue
   !   write(2,*) 'entering SelFileName', present(SelFileName)
   !   flush(unit=2)
   If(present(SelFileName) .and. N_bins.gt.1) then
      !write(2,*) (d_AmplScale-1.),Ampl(1)
      !write(2,*) 'entering pre-logscaling', N_bins
      !flush(unit=2)
      a=Ampl(1) ; b=Ampl(N_bins)  ! set scale for plots to cover at least factor 10
      !write(2,*) 'entering logscaling', a,b
      !flush(unit=2)
      If(b/a.lt.10.) Then
         c=sqrt(10.*a/b)
         a=a/c ; b=b*c
      EndIf
      c=log(a)/log(10.)
      !write(2,*) 'c=log(a)/log(10.):',c,a
      If(c.lt.0) Then
         i=int(c)-1 ;
      Else
         i=int(c)
      EndIf
      c=10.**i  ! next lower power of 10
      If(a.gt.5.*c) Then
         a=5.*c
      ElseIf(a.gt.3.*c) Then
         a=3.*c
      ElseIf(a.gt.2.*c) Then
         a=2.*c
      Else
         a=c
      EndIf
      !write(2,*) i,c,a
      If(a.lt.0.01) a=0.01
      c=log(b)/log(10.)
      !write(2,*) 'c=log(b)/log(10.):',c,b
      If(c.lt.0) Then
         i=int(c)-1 ;
      Else
         i=int(c)
      EndIf
      c=10.**i  ! next lower power of 10
      If(b.lt.2.*c) Then
         b=2.*c
      ElseIf(b.lt.3.*c) Then
         b=3.*c
      ElseIf(b.lt.5.*c) Then
         b=5.*c
      Else
         b=10.*c
      EndIf
      !write(2,*) i,c,b
      !  End log amplitude scale
      Nrm1=0.  ; Nrm2=0.
      !write(2,*) 'N_bins:',N_bins
      OPEN(unit=28,FILE=TRIM(SelFileName)//'.dat', FORM='FORMATTED',STATUS='unknown',ACTION='write') ! space separated values
      write(28,"(1x,F8.2,2F10.2, 6(1x,g10.4),f7.3,1x,A,A)") a, Ampl(i_MaxFit), b, &
         a1,b1,c1,a2,b2,c2,d2, TRIM(SelFileName),' ! by AmplitudeFit' ! gleRead:  NTracks EventNr Q t_start label$ AmplitudePlot a1 b1 c1 a2 b2 c2 d2
      Do i_Ampl=1,i_AmplMax !N_bins ! i_MaxFit ! i_AmplMax  !  finally write all accepted source to file
         !Ampl10=(Ampl(1)/NrSources*(d_AmplScale-1.))*Ampl(i_Ampl)*Ampl(i_Ampl)
         Ampl10= Ampl(i_Ampl)*Ampl(i_Ampl)/HistNrm  ! norm changed 18 Lan 2023
         !Ampl10=1./(NrSources*Ampl(1)) * Ampl(i_Ampl)*Ampl(i_Ampl)
         !write(2,"(1x,I5,7(1x,g14.6))") i_Ampl, Ampl(i_Ampl), Ampl10*Ampl_Hist(i_Ampl), Ampl10, Ampl_Hist(i_Ampl)
         write(28,"(1x,I5,7(1x,g14.6))") i_Ampl, Ampl(i_Ampl), Ampl10*Ampl_Hist(i_Ampl),Ampl10/Ampl_Hist_W(i_Ampl), &
               Ampl10*Fit_exp(i_Ampl)+Tiny, Ampl10*Fit_pow(i_Ampl)+Tiny, Fit_pow(i_Ampl)
         Nrm1=Nrm1+Fit_exp(i_Ampl)*(Ampl(i_Ampl+1)-Ampl(i_Ampl))
         Nrm2=Nrm2+Fit_pow(i_Ampl)*(Ampl(i_Ampl+1)-Ampl(i_Ampl))
      enddo
      close(unit=28)
   EndIf
   Nrm=NrSources*Ampl(1)
   write(2,"(1x,A, f6.2, A, f7.4, g10.3, g11.3, A,2g11.3)") '$b *exp(-a*A-c/A^2); \chi^2=$',ChiSq1, &
                                    ',with  a,b,c= ',a1,b1/Nrm,c1,'; nrm=', Nrm, HistNrm
   write(2,"(1x,A, f6.2, A, f7.3, g10.3, g11.3, A,2g11.3)") '$b *A^-a *exp(-c/A); \chi^2=$',ChiSq2, &
                                    ',with  a,b,c= ',a2,b2/Nrm,c2,'; nrm=', Nrm
   Return
   !
4   FitFunc=4 ; nvar= 3  ! F(A)=b A^(-a)*exp(-c/A^2) with i_Ampl=A/d_Ampl  where A is plotted
   i_AmplTh=3 !SourceAmp(NrSources*1/10)*d_AmplScale  ! Threshold for fitting in histogram index
   Meqn = Max_Ampl-i_AmplTh-1      ! number of equations
   X(1)=4.
   X(2)=Ampl_Hist(i_AmplTh)*(i_AmplTh/10.)**X(1)
   !Call FitAmpl(Meqn, nvar, X,ChiSQDF)
   a2=X(1)
   b2=x(2)
   d2=x(3)
   c2=0.
   Return
   !
3   FitFunc=3 ; nvar= 3 ! b*exp(-a*A-c/A^1.5)  a=x(1)/d_Ampl
   !Call FitAmpl(Meqn, nvar, X,ChiSQDF)
   a2=X(1)
   b2=exp(x(2))
   c2=x(3)
   d2=0.
   !
!   Return
   !write(2,*) MeanTrackLoc(1,1:t_MTL)
   !Call Flush(2)
   !close(unit=27)
    !
   Return
End Subroutine AmplitudeFit
!=================================
End Module Tracks
!===============================================
!===============================================
    Subroutine GetNonZeroLine(lineTXT)
    implicit none
    character(LEN=*), intent(inout) :: lineTXT
    Character(len=6) :: FMT
    integer :: eof,N
      N=LEN(lineTXT)
      write(FMT,"(A,I3.3,A)") '(A',N,')'
      !write(2,*) 'GetNonZeroLine2:', FMT, N
      lineTXT=''
      Do while (trim(lineTXT).eq.'')
         read(*,FMT,iostat=eof) lineTXT
         if(eof.lt.0) exit
      enddo
    End Subroutine GetNonZeroLine
!==========================================
!
