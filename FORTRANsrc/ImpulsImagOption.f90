Subroutine ImpulsImagRun
   use constants, only : dp,pi,ci,sample,c_mps
   use DataConstants, only : ProgramFolder, UtilitiesFolder, FlashFolder, DataFolder, FlashName, Windows
   use DataConstants, only : Time_dim, Production, release  ! , Cnu_dim, RunMode, Utility
   use DataConstants, only : Calibrations, OutFileLabel, EdgeOffset
   use ThisSource, only : CCShapeCut_lim, ChiSq_lim, EffAntNr_lim, Dual, Safety !, PeakPos
   use FitParams, only : FullSourceSearch, SigmaGuess, AntennaRange, SearchRangeFallOff, Sigma_AntT
   use Chunk_AntInfo, only : StartT_sam, AntennaNrError, DataReadError, TimeFrame !, BadAnt_nr, BadAnt_SAI ExcludedStatID,
   use Chunk_AntInfo, only : NoiseLevel, PeaksPerChunk !TimeBase, Simulation, WriteSimulation
   use Chunk_AntInfo, only : Station_OutOfRange, Unique_StatID, Nr_UniqueStat
   use FFT, only : RFTransform_su,DAssignFFT
   use Explore_Pars, only : NMin, NMax, Emin, EMax
   use StationMnemonics, only : Statn_ID2Mnem
   use GLEplots, only : GLEplotControl
   Implicit none
   INTEGER :: DATE_T(8),i
   Integer :: j,i_chunk, ChunkNr_start, ChunkNr_stop, units(0:2), FitRange_Samples!, CurtainWidth
   Integer :: i_dist, i_guess, nxx, valueRSS
   Real*8 :: StartTime_ms, StartingTime, StoppingTime, D
   Real*8 :: SourceGuess(3,1) ! = (/ 8280.01,  -15120.48,    2618.37 /)     ! 1=North, 2=East, 3=vertical(plumbline)
   CHARACTER(LEN=1) :: Mark
   CHARACTER(LEN=10) :: Sources
   Character(LEN=180) :: Sources1, Sources2
   Character(LEN=180) :: lname
   Character(LEN=250) :: TxtIdentifier, TxtImagingPars
   Character(LEN=250) :: Txt20Identifier, Txt20ImagingPars
   Real(dp), external :: tShift_ms
   !
   CALL DATE_AND_TIME (Values=DATE_T)
   Write(TxtIdentifier,"('Release ',A,', run on ',I2,'/',I2,'/',I4,' , at ',I2,':',I2,':',I2,'.',I3,', Flash: ',A )") &
       TRIM(release),DATE_T(3),DATE_T(2),DATE_T(1),(DATE_T(i),i=5,8),TRIM(FlashName)
   !
   !Call GetNonZeroLine(lname)
   Call GetMarkedLine(Mark,lname)
   write(2,"(A)") 'Input line-1: "'//Mark//'|'//TRIM(lname)// &
         '" !  Reference/Source-| t_0 & position, ImagingStartingTime[ms] (after t_0), ImagingStoppingTime[ms] '
   Read(lname,*,iostat=nxx) StartTime_ms, SourceGuess(:,1), StartingTime, StoppingTime  ! Start time offset = 1150[ms] for 2017 event
   Write(2,*) 'Better points to the center of the region to be imaged, to minimize out-of-range antennas.'
   If(nxx.ne.0) then
      write(2,*) 'Parsing error in input line'
      stop 'ImpulsiveImager; reading error'
   Endif
   Call Convert2m(SourceGuess(:,1))
   If(Mark.eq."S") Then
      !StartTime_ms=StartTime_ms + sqrt(SUM(SourceGuess(:,1)*SourceGuess(:,1)))*1000.*RefracIndex(SourceGuess(3,1))/c_mps ! convert to reference antenna time
      StartTime_ms=StartTime_ms + tShift_ms(SourceGuess(:,1)) ! convert to reference antenna time
   EndIf
   Write(TxtImagingPars,"(A,I3,A,F6.1,A,F4.1,A,F4.1,A,F4.0,A,F5.0,A,F5.1,A,F6.1,A,3F5.1,A,i5)") &
      'FitRange=',Safety,'[samples],',AntennaRange,'[km], Sigma_AntT=',Sigma_AntT, &
      '[samples], SearchRangeFallOff-factor=',SearchRangeFallOff,', CCShapeCut_lim=',CCShapeCut_lim*100., &
      '%, ChiSq_lim=',ChiSq_lim,'[ns^2], EffAntNr_lim=',EffAntNr_lim*100., '%, Noise=',NoiseLevel, &
      ', Guess(N,E,h)=(', SourceGuess(:,1)/1000.,')[km], PeaksPerChunk=',PeaksPerChunk
   If(.not.FullSourceSearch) then
      Call GetNonZeroLine(lname)
      write(2,*) 'Imaging input line-2:',lname
      read(lname,*,iostat=nxx) SigmaGuess !, PeakPos(1), PeakPos(2)  ! Start time offset = 1150[ms] for 2017 event
      If(nxx.ne.0) SigmaGuess=(/ 3000.d0,3000.d0,2000.d0 /)
      Write(2,"('Preferred source location:',3F8.2,'[km]')") SourceGuess(:,1)/1000.
      Write(2,"('Initial error on location:',3F8.3,'[km]')") SigmaGuess(:)/1000.
   Endif
   StartT_sam(1)=(StartTime_ms/1000.d0)/sample  ! in sample's
   write(*,*) 'Start Source Finding'
   lname=' [s] ; ID, (North, East, vertical at core) [m] , t [s] , chi^2, sigma(N,E,v),    N_EffAnt, Max_EffAnt ;'
   lname=lname//' height, l&r widths'
   write(Txt20ImagingPars,"(A,1x,F8.2,1x,F8.2,1x, F8.2)") TRIM(FlashName), StartTime_ms, StartingTime, StoppingTime
   write(Txt20Identifier,"('Chnk#  , t[ms] , max,found, <5^2 ,',2x,A)") TRIM(release)
   Sources='Srcs'//release(1:2)
   If(Dual) then
      Sources1=TRIM(Sources)//'-dbl'//TRIM(OutFileLabel)
      Open(unit=17,STATUS='unknown',ACTION='write', FILE = TRIM(Sources1)//'.csv')
      Write(17,"('! ',A,A,A)") TRIM(TxtIdentifier),'; ',TRIM(Calibrations)
      Write(17,"('! ',A)") TRIM(TxtImagingPars)
      Write(17,"(f12.9,2A)") StartTime_ms/1000.d0,TRIM(lname),' even&odd antenna numbers'
      Open(unit=27,STATUS='unknown',ACTION='write', FILE = TRIM(DataFolder)//TRIM(Sources1)//'_stat.csv')
      write(27,*) TRIM(Txt20ImagingPars)
      write(27,"('! ',A)") TRIM(Txt20Identifier)
      units(2)=17
   Else
      Sources1=TRIM(Sources)//'-even'//TRIM(OutFileLabel)
      Open(unit=15,STATUS='unknown',ACTION='write', FILE = TRIM(Sources1)//'.csv')
      Write(15,"('! ',A,A,A)") TRIM(TxtIdentifier),'; ',TRIM(Calibrations)
      Write(15,"('! ',A)") TRIM(TxtImagingPars)
      Write(15,"(f12.9,2A)") StartTime_ms/1000.d0,TRIM(lname),' even antenna numbers'
      Open(unit=25,STATUS='unknown',ACTION='write', FILE = TRIM(DataFolder)//TRIM(Sources1)//'_stat.csv')
      write(25,*) TRIM(Txt20ImagingPars)
      write(25,"('! ',A)") TRIM(Txt20Identifier)
      units(0)=15
      !
      Sources2=TRIM(Sources)//'-odd'//TRIM(OutFileLabel)
      Open(unit=16,STATUS='unknown',ACTION='write', FILE = TRIM(Sources2)//'.csv')
      Write(16,"('! ',A,A,A)") TRIM(TxtIdentifier),'; ',TRIM(Calibrations)
      Write(16,"('! ',A)") TRIM(TxtImagingPars)
      Write(16,"(f12.9,2A)") StartTime_ms/1000.d0,TRIM(lname),' odd antenna numbers'
      Open(unit=26,STATUS='unknown',ACTION='write', FILE = TRIM(DataFolder)//TRIM(Sources2)//'_stat.csv')
      write(26,*) TRIM(Txt20ImagingPars)
      write(26,"('! ',A)") TRIM(Txt20Identifier)
      units(1)=16
   Endif
   !
   ChunkNr_start=StartingTime/(1000.d0*sample*(Time_dim-2*EdgeOffset))+1
   ChunkNr_stop=StoppingTime/(1000.d0*sample*(Time_dim-2*EdgeOffset))+1
   Station_OutOfRange(:)=0
   !ChunkNr_start=1501 ; ChunkNr_stop=1900 ! No Time Frame == ChunkNr
   write(*,"(A,i5,A)") achar(27)//'[45m # of data-blocks read in=',(ChunkNr_stop-ChunkNr_start),achar(27)//'[0m'
   If((ChunkNr_stop-ChunkNr_start).gt.3) Production=.true.
   Open(unit=18,STATUS='unknown',ACTION='write', FILE = TRIM(Sources)//'-5star-'//TRIM(OutFileLabel)//'.dat')
   write(18,*) '! StartT_sam[ms]:', '1000.d0*StartT_sam(1)*sample, (TimeFrame-1)*(Time_dim-2*EdgeOffset)*1000.*sample, ', &
      'Peakpos_0, SourcePos(:,1)/1000., FitQual, 100.*N_EffAnt/Max_EffAnt, PeakSAmp(i_peakS,i_eo), Wl, Wu,i_eo'
   Write(18,"(' ', A13, 3A11,A9)")  'StartT','N[km]','E[km]','h[km]','label'
   Write(18,"('C # 2 1',A8,3(A10,','),A12,';',A9,',',3(A8,','),2(A4,','),A7,',',2(A3,','),A3)") &
      'Peakpos','N[m]','E[m]','h[m]','t[ms]','Chi^2','N_a','N_mx','Ampl','W_l','W_u','eo'
!                  Write(18,"('C',i2,' 2 1',I8,3(F10.2,','),F12.5,';',f9.3,',',2(I4,','),I7,',',2(I3,','),2I3)") &
!                     i_FvStr, Peakpos_0, SourcePos(:,1), &
!                     (StartT_sam(1)+Peakpos_0)*sample-DistMax*RefracIndex(SourcePos(3,1))/c_mps , FitQual &
!                     ,  N_EffAnt, Max_EffAnt, PeakSAmp(i_peakS,i_eo), Wl, Wu, i_eo, i_peakS
   Do TimeFrame=ChunkNr_start, ChunkNr_stop
      write(*,"(A,i6,A)", ADVANCE='NO') achar(27)//'[31m Processing block# ',TimeFrame,achar(27)//'[0m'  ! [1000D    !  //achar(27)//'[0m.'
      i_chunk=1
      StartT_sam(1)=(StartTime_ms/1000.d0)/sample + (TimeFrame-1)*(Time_dim-2*EdgeOffset)  ! in sample's
      !Write(2,*) 'StartTime_ms=',StartTime_ms,', TimeFrame=',TimeFrame,', StartT_sam=',StartT_sam(1)
      Write(2,"(A,I5,A,F8.3,A,F12.6,A,I11,A,3F10.1,A)") 'Start time for slice',TimeFrame &
         ,' (dt=',StartT_sam(1)*1000.d0*sample-StartTime_ms,' [ms]) is',StartT_sam(1)*1000.d0*sample &
         ,' [ms]=',StartT_sam(1),' [Samples]; SourceGuess=',SourceGuess(:,i_chunk),';  1=North, 2=East, 3=vertical(plumbline)'
      Call RFTransform_su(Time_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        !Source_Crdnts= (/ 10000 , 16000 , 4000 /)    ! 1=North, 2=East, 3=vertical(plumbline)
      Call AntennaRead(i_chunk,SourceGuess(:,i_chunk))
      Call DAssignFFT()
      If(DataReadError.ne.0) cycle  !  set in AntennaRead
      If(AntennaNrError.gt.0) cycle  !  set in AntennaRead
      !
      !write(18,*) '5star: processing chunk nr',TimeFrame
      Call SourceFind(TimeFrame,SourceGuess,units)
      !
      flush(unit=2)
      !
   EndDo !  i_chunk
   !
   If(Dual) Then
      Close(Unit=17)
      Close(Unit=27)
   Else
      Close(Unit=15)
      Close(Unit=16)
      Close(Unit=25)
      Close(Unit=26)
   Endif
   close(unit=14)
   close(unit=12)
   !
   Do j=1,Nr_UniqueStat
      write(2,*) 'Station_included:',Statn_ID2Mnem(Unique_StatID(j)),Unique_StatID(j) ! , ' #:',Station_OutOfRange(j)
   EndDo
   !write(*,*) 'Station_OutOfRange:',Station_OutOfRange(1:Nr_UniqueStat)
   write(*,*) 'BoundingBox (N,E):',NMin, NMax, Emin, EMax
   If((NMax-Nmin).gt.(Emax-Emin)) Then
      NMin=NMin -.5
      NMax=NMax +.5
      D=(EMax+Emin)/2.
      EMax=D + (NMax-Nmin)/2.
      EMin=D - (NMax-Nmin)/2.
   Else
      EMin=EMin -.5
      EMax=EMax +.5
      D=(NMax+Nmin)/2.
      NMax=D + (EMax-Emin)/2.
      NMin=D - (EMax-Emin)/2.
   Endif
   lname='AFlashImage.in'  ! temporary name holder
   Open(UNIT=10,STATUS='unknown',ACTION='WRITE',FILE=TRIM(lname))
   If(Dual) then
      write(10,"(A,A)") '&Parameters  datafile="'//TRIM(Sources1)//'", PlotName="'//TRIM(OutFileLabel)//'d", ', &
         'CutSigmaH=17., LinCutH=2., RMS_ns= 3.5 , DelNEff=25,  AmplitudePlot=10. '
      write(10,"(A,4f6.1,' 0.0 12. ',2f7.1,' &End')") ' NEhtBB= ', &
            NMin, NMax, Emin, EMax, StartingTime, StoppingTime !  plotting command
      write(10,"('======================================================')")
   Else
      write(10,"(A,A)") '&Parameters  datafile="'//TRIM(Sources1)//'", PlotName="'//TRIM(OutFileLabel)//'e", ', &
         'CutSigmaH=17., LinCutH=2., RMS_ns= 3.5 , DelNEff=25,  AmplitudePlot=10. '
      write(10,"(A,4f6.1,' 0.0 12. ',2f7.1,' &End')") ' NEhtBB= ', &
            NMin, NMax, Emin, EMax, StartingTime, StoppingTime !  plotting command
      write(10,"(A,A)") '&Parameters  datafile="'//TRIM(Sources2)//'", PlotName="'//TRIM(OutFileLabel)//'o", ', &
         'CutSigmaH=17., LinCutH=2., RMS_ns= 3.5 , DelNEff=25,  AmplitudePlot=10. '
      write(10,"(A,4f6.1,' 0.0 12. ',2f7.1,' &End')") ' NEhtBB= ', &
            NMin, NMax, Emin, EMax, StartingTime, StoppingTime !  plotting command
      write(10,"('======================================================')")
   EndIf
   Close(unit=10)
   !
   If(Dual) then
      Call GLEplotControl(PlotType='PeakStats', PlotName='PeakStatsD'//TRIM(OutFileLabel), &
         PlotDataFile=TRIM(DataFolder)//TRIM(Sources1) )
   Else
      Call GLEplotControl(PlotType='PeakStats', PlotName='PeakStatsE'//TRIM(OutFileLabel), &
         PlotDataFile=TRIM(DataFolder)//TRIM(Sources1) )
      Call GLEplotControl(PlotType='PeakStats', PlotName='PeakStatsO'//TRIM(OutFileLabel), &
         PlotDataFile=TRIM(DataFolder)//TRIM(Sources2) )
   EndIf
   !
   If(Windows) Then
      Call GLEplotControl(SpecialCmnd='call %LL_Base%\scripts\RunProgram.bat DataSelect '//TRIM(lname))
      !call %LL_Base%\scripts\RunProgram.bat DataSelect DataSelect.in
   Else
      Call GLEplotControl(SpecialCmnd='source  ${LL_Base}/scripts/RunProgram.sh DataSelect '//TRIM(lname))
      !source  ${LL_Base}/scripts/RunProgram.sh DataSelect  DataSelect.in
   EndIf
   !
   !If(Dual) then
   !   Call GLEplotControl(PlotType='SourcesPlot', PlotName='Map_'//TRIM(FlashName)//'d'//TRIM(OutFileLabel), &
   !      PlotDataFile=TRIM(DataFolder)//TRIM(FlashName)//'d'//TRIM(OutFileLabel), Submit=.true.)
   !Else
   !   Call GLEplotControl(PlotType='SourcesPlot', PlotName='Map_'//TRIM(FlashName)//'e'//TRIM(OutFileLabel), &
   !      PlotDataFile=TRIM(DataFolder)//TRIM(FlashName)//'e'//TRIM(OutFileLabel))
   !   Call GLEplotControl(PlotType='SourcesPlot', PlotName='Map_'//TRIM(FlashName)//'o'//TRIM(OutFileLabel), &
   !      PlotDataFile=TRIM(DataFolder)//TRIM(FlashName)//'o'//TRIM(OutFileLabel), Submit=.true.)
   !EndIf
   !
   Return
End Subroutine ImpulsImagRun
