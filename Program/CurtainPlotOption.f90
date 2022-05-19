! ========================
Subroutine PlotAllCurtainSpectra(CurtainWidth)
!  Make Curtain Plots
!  Option 'Dual' is assumed and thus only the i_eo=0 peaks are taken
   ! peaknr is included in the naming convention when plots for multple peaks are needed (=calibration runs)
   use constants, only : dp,sample
   use DataConstants, only : Station_nrMax, Time_dim, DataFolder, OutFileLabel  ! , ChunkNr_dim
   use ThisSource, only : PeakNrTotal, Peak_eo, ChunkNr, Peakpos, SourcePos
   use Chunk_AntInfo, only : CTime_spectr, Ant_Pos, Ant_nr, Ant_RawSourceDist, Unique_StatID,  Nr_UniqueStat, Ant_Stations, RefAnt
   use StationMnemonics, only : Statn_ID2Mnem, Statn_Mnem2ID
   use FitParams, only : Fit_TimeOffsetStat ! , N_FitPar_max, Fit_AntOffset, Fit_TimeOffsetAnt, PulsPosCore
   use GLEplots, only : GLEplotControl
   Implicit none
   !
   integer, intent(in) :: CurtainWidth
   Real(dp) :: B, RDist, OSt
   Integer :: j, Lower, Upper, Ch1, Ch2, OffSet
   Integer :: unt
   Integer :: i_ant, i_time, i_eo, i_chunk, i_stat, i_peak
   Logical :: UsePeakNr
   CHARACTER(LEN=60) :: GLE_file
   ! Character(len=25), save :: OutFileLabel=''
   CHARACTER(LEN=32) :: txt      ! 7 more than length of  OutFileLabel
   !
   write(*,*) 'MakeAllCurtainPlots for "'//TRIM(OutFileLabel)//'"'
   write(2,*) 'curtain width=',CurtainWidth
   !Open(UNIT=10,STATUS='unknown',ACTION='WRITE',FILE='GLE-plots.sh')
      !write(2,*) 'number of antennas=',Ant_nr(i_chunk)
   ! check if curtain plots for multiple peaks are to be made
   UsePeakNr=.true.
   If(PeakNrTotal.le.2) Then
      If(PeakNrTotal.eq.1) Then
         UsePeakNr=.false.
      Else
         If(Peak_eo(PeakNrTotal).eq.1) UsePeakNr=.false.
      EndIf
   EndIf
   Do i_peak=1,PeakNrTotal
      i_eo=Peak_eo(i_peak)
      !write(2,*) 'CurtainPlot:I_peak-A',i_peak,PeakNrTotal,I_ant,Nr_UniqueStat,i_chunk,i_eo
      If(i_eo .ne. 0) cycle !  Option 'Dual' is assumed and thus only the i_eo=0 peaks are taken
      i_chunk=ChunkNr(i_peak)
      I_ant=RefAnt(i_chunk,i_eo)
      !write(2,*) 'CurtainPlot:I_peak',i_peak,PeakNrTotal,I_ant,Nr_UniqueStat,i_chunk,i_eo
      !flush(unit=2)
      Do i_stat=1,Nr_UniqueStat  ! Get station number from unique list to retrieve timing-offset
         If(Unique_StatID(i_stat).eq. Ant_Stations(i_ant,i_chunk)) exit
      Enddo
      Call RelDist(SourcePos(1,i_Peak),Ant_pos(1,i_ant,i_chunk),RDist)
      OSt=(Rdist - Ant_RawSourceDist(i_ant,i_chunk) + Fit_TimeOffsetStat(i_stat)) ! to correct for off-center position reference antenna
      Do I_ant=1,Ant_nr(i_chunk)
         !write(2,*) 'CurtainPlot:I_ant',I_ant,Ant_nr(i_chunk)
         !flush(unit=2)
         Do i_stat=1,Nr_UniqueStat  ! Get station number from unique list to retrieve timing-offset
            If(Unique_StatID(i_stat).eq. Ant_Stations(i_ant,i_chunk)) exit
         Enddo
         !   A=MaxLoc(Stations, mask=Stations.eq.Ant_Stations(I_ant,i_chunk))
         !   i_stat=A(1)
         !
         Call RelDist(SourcePos(1,i_Peak),Ant_pos(1,i_ant,i_chunk),RDist)
         OffSet=NINT(Rdist - Ant_RawSourceDist(i_ant,i_chunk) + Fit_TimeOffsetStat(i_stat) -OSt)
         !write(2,*) 'curtainoffset',OffSet, Rdist - Ant_RawSourceDist(i_ant,i_chunk)
         Ch2=PeakPos(i_Peak) + OffSet
         Lower=MIN(CurtainWidth, Ch2-1)
         Upper=MIN(CurtainWidth,(Time_dim - Ch2))
         Ch1=Ch2-Lower
         Ch2=Ch2+Upper
         !write(2,*) I_ant,i_peak, Ch1, Ch2, Lower, Upper,Rdist, Ant_RawSourceDist(i_ant,i_chunk), Fit_TimeOffsetStat(i_stat), OSt
         If(Upper.le.0) then
            Ch1=Time_dim/2
            Ch2=Ch1+2
         Endif
         If(Ch2.lt.Ch1) Ch2=Ch1+2
         If(UsePeakNr) Then  !  "'//TRIM(OutFileLabel)//'"
            write(txt,"(A,i3.3,'-',i3.3)") TRIM(OutFileLabel),I_ant,i_peak
         Else
            write(txt,"(A,i3.3)") TRIM(OutFileLabel),I_ant
         EndIf
         !write(2,*) 'working on curtainplot: "','LOFAR_Time'//TRIM(txt),'"'
         !
         B=maxval(real(CTime_spectr(ch1:ch2,i_Ant,i_chunk)))
         Open(Unit=9,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//'LOFAR_Time'//TRIM(txt)//'.dat')
         Do i_time=ch1,ch2
           Write(9,*) i_time-Offset, Real(CTime_spectr(i_time,i_Ant,i_chunk))/B ! written time corresponds to time at core for this source
         Enddo
         close(unit=9)
      Enddo
      unt=9
      write(2,*) 'SourcePos(:,i_Peak)',SourcePos(:,i_Peak)
      If(UsePeakNr) Then  !  "'//TRIM(OutFileLabel)//'"
         Write(GLE_file,"('CuP-',A,i3.3)") TRIM(OutFileLabel),i_peak
      Else
         Write(GLE_file,"('CuP-',A)") TRIM(OutFileLabel)
      EndIf
      Write(2,*) 'GLE-file name:',TRIM(GLE_file)
      Call GLEscript_CurtainPlot(unt, GLE_file, CurtainWidth, i_Peak, UsePeakNr)
      !Write(10,"(A)") 'gle -d pdf '//trim(GLE_file)//'.gle'
   Enddo ! i_peak
   !
   Call GLEplotControl(SpecialCmnd='rm '//TRIM(DataFolder)//'LOFAR_Time'//TRIM(OutFileLabel)//'*.dat') !  Command to delete curtain files
   !   Call GLEplotControl(SpecialCmnd='rm '//TRIM(DataFolder)//TRIM(GLE_file)//'.gle') !  Command to delete curtain files
   Return
End Subroutine PlotAllCurtainSpectra
! ===========================================================================
Subroutine GLEscript_CurtainPlot(unt, file, CurtainWidth, i_Peak, UsePeakNr)
   use constants, only : dp
    use DataConstants, only : Time_dim, DataFolder, OutFileLabel
    use Chunk_AntInfo, only : CTime_spectr, Ant_Stations, Ant_nr, Ant_IDs, Ant_pos, Ant_RawSourceDist, Nr_UniqueStat
    use DataConstants, only : Station_nrMax
    use StationMnemonics, only : Station_ID2Mnem
    use ThisSource, only : PeakPos, ChunkNr, ExclStatNr
    use GLEplots, only : GLEplotControl
    Implicit none
    Integer, intent(in) :: unt, CurtainWidth, i_Peak
    Logical, intent(in) :: UsePeakNr
    Character(Len=*), intent(in) :: file
    integer :: I_ant, i_stat, Stat_nr, Stations(1:Station_nrMax), i_time, i_c, lstyle
    integer :: A(1), Plot_scale, ch1, ch2, VLine, i_chunk, k
    Character(len=42) :: txt
    Character(len=5) :: Station_Mnem
    real*8 :: plot_offset,b
    !
    i_chunk=ChunkNr(i_peak)
    Plot_scale= 2. ! int(B) + 2.
    VLine=PeakPos(i_Peak)
    ch1=PeakPos(i_Peak)-CurtainWidth
    ch2=PeakPos(i_Peak)+CurtainWidth
    !vsize=Nr_UniqueStat*2 34
    !
    !write(2,*) 'GLEscript_CurtainPlot:Nr_UniqueStat=', Nr_UniqueStat
    Open(UNIT=UNT,STATUS='unknown',ACTION='WRITE',FILE=trim(file)//'.gle')
    Write(unt,"(A)") '! COMMAND:  gle -d pdf '//trim(file)//'.gle'
    Write(unt,"(A,I2,3(/A),3(/A,I0))") 'size 63 ',2*Nr_UniqueStat+8,'set font pstr fontlwidth 0.08 hei 1.2 just CC',&
        'include "../Utilities/DiscreteColor.gle"',&
        'set lwidth 0.1','t_min = ',ch1,'t_max = ',ch2,'scl = ',Plot_scale
!        'set lwidth 0.1','t_min = 13250 !12500 !0','t_max = 13350 !15000 !',Time_dim,'scl = ',Plot_scale
     If(UsePeakNr) Then  !  "'//TRIM(OutFileLabel)//'"
         write(txt,"(A,i1.1,A,i3.3)") 'Chunk=',i_chunk,', peak=',i_peak
     Else
         write(txt,"(A)") TRIM(OutFileLabel)
     EndIf
    Write(unt,"(A,2(/A),i2,3(/A),7(/A))")  'amove 4 4','begin graph', &
        '   size 55 ',2*Nr_UniqueStat+1,'   vscale 1','  hscale 1',&
        '   title  "Time Spectra, '//TRIM(txt)//'"', &
        '   xtitle "time [samples]"','   ytitle "Amplitude"',&
        '   xaxis min t_min max t_max ! nticks 10','   yaxis  min 0 max 2 dticks 5 dsubticks 1', &
        '   x2labels on',' end graph'
    Write(unt,"(A,6(/A))" ) 'begin key','position tc','nobox','   offset 20. -1.5', &
         '   text "Excluded" lstyle 4 color black lwidth 0.1', &
         'end key','set just cc'
    Write(unt,"(A,i5.5,A,1(/A),I2)")  'amove xg(',Vline,') 4','rline 0 ',2*Nr_UniqueStat+1
    !
    Stat_nr=0
    Stations(:)=0
    Do I_ant=1,Ant_nr(i_chunk)
        A=MinLoc(Ant_Stations(:,i_chunk), mask=Ant_Stations(:,i_chunk).eq.Ant_Stations(I_ant,i_chunk))
        !write(2,*) I_ant,Ant_Stations(I_ant),A(1)
        If(A(1).eq.I_ant) then
            Stat_nr=Stat_nr+1
            Stations(Stat_nr)=Ant_Stations(I_ant,i_chunk)
            i_stat=Stat_nr
            !write(2,*) 'stat_nr=',Stat_nr
        else
            A=MaxLoc(Stations, mask=Stations.eq.Ant_Stations(I_ant,i_chunk))
            i_stat=A(1)
            !write(2,*) 'i_stat=',i_Stat,Stations(1:Stat_nr)
        endif
        If(UsePeakNr) Then  !  "'//TRIM(OutFileLabel)//'"
            write(txt,"(A,i3.3,'-',i3.3)") TRIM(OutFileLabel),I_ant,i_peak
        Else
            write(txt,"(A,i3.3)") TRIM(OutFileLabel),I_ant
        EndIf
        !write(txt,"(i3.3,'-'i2.2)") I_ant,i_peak
        plot_offset=2*i_stat + 2 + mod(Ant_IDs(I_ant,i_chunk),2)
        i_c = I_Ant-11*int(I_Ant/11)
        !
        lstyle=0
        Do k=1,Station_nrMax
            If(ExclStatNr(k,i_peak).eq.0) exit
            If(Stations(I_stat).eq.ExclStatNr(k,i_peak)) Then
               lstyle=4  !  ! Marks excluded station, see "ReImAtMax"
               exit
            EndIf
        Enddo
        !
        Write(Unt,901) plot_offset,trim(DataFolder), TRIM(txt),I_Ant,I_Ant,lstyle,i_c
901     Format('amove 4 ',F5.2,/'begin graph',/'  size 55 2',/'  vscale 1',&
            /'  hscale 1',/'   NoBox',&
            /'   xaxis min t_min max t_max',/'   x1axis off',&
            /'   yaxis  min -scl max scl',/'   y1axis off',&
            /'     data "',A,'LOFAR_Time',A,'.dat" d',I0,'=c1,c2',&
            /'     d',I0,' line lwidth 0.04 lstyle ',i1,' color MyCol',i0,'$',&
            /' end graph')
        !     let d9=0
        !     d9 line lstyle 1 lwidth 0.01 color black
        !set font pstr fontlwidth 0.06 hei 0.8 just CC
        !set lwidth 0.1
        !begin key
        !nobox
        !position tr
        !text "r=0-1 m" line color red lstyle 1 lwidth 0.1
        !end key
    Enddo
    Do i_stat=1,Stat_nr
        Call Station_ID2Mnem(Stations(I_stat),Station_Mnem)
        plot_offset=2*i_stat  + 3
        Write(unt,902) plot_offset,plot_offset,Station_Mnem
902 Format('amove 4 ',F5.2,/'aline 59 ',F5.2,/'rmove 2 0 ',/'write "',A5,'"')
    Enddo
    !
    Close(unit=unt)
    !
    Call GLEplotControl(SpecialCmnd='gle -d pdf '//trim(file)//'.gle') !  Command to produce curtain plots
    Call GLEplotControl(SpecialCmnd='rm '//TRIM(file)//'.gle') !  Command to delete curtain files
    !
    return
End Subroutine GLEscript_CurtainPlot
!====================================
Subroutine GLE_Corr()
    use DataConstants, only : Station_nrMax, DataFolder
    use DataConstants, only :  Ant_nrMax
    use Chunk_AntInfo
    use ThisSource, only : CCorr, PeakNrTotal, CorrAntNrs, Nr_Corr, Safety, Peak_eo, ChunkNr, Peakpos, CCorr_Err, RefAntErr
    use FitParams, only : N_FitPar_max
    use constants, only : dp
    use StationMnemonics, only : Station_ID2Mnem
    use GLEplots, only : GLEplotControl
    Implicit none
    Integer :: J_Corr
    Integer ::  i_ant, i_tp, i_tp2, pbase
    !
    Integer ::  i_Peak, i_eo, Station, i_c, i_chunk, lstyle
    Integer :: i, i_max, k, J_corr_st(0:Station_nrMax),Height,lr
    Logical :: nra,nrb
    character(len=5) :: Station_Mnem,Lab_a,Lab_b
    character(len=10) :: txt
    character(len=13) :: XCorrPlot
    Character(len=2) :: EveOdd(0:1)=(/'Ev','Od'/), tp(0:1)=(/'Th','Ph'/), ext
    Real(dp) :: plot_offset, lwidth
    !
    !Open(UNIT=10,STATUS='unknown',ACTION='WRITE',FILE='GLE-plots.sh')
    i_tp2=0
    !If(polariz) i_tp2=1
    Do i_peak=1,PeakNrTotal
      Do i_tp=0,i_tp2
        i_eo=Peak_eo(i_peak)
        i_chunk=ChunkNr(i_peak)
        Station=Ant_Stations(CorrAntNrs(1,i_eo,i_chunk),i_chunk)
        pBase=0
        !If(polariz) then
        !    ext=tp(i_tp)
        !    If(i_tp.eq.1) pBase=Ant_nrMax/2
        !    write(XCorrPlot,"('XCP_',I3.3,A2,'.gle')") i_peak,ext
        !Else
            ext=EveOdd(i_eo)
            write(XCorrPlot,"('XCP_',I3.3,'.gle')") i_peak
        !EndIf
        i=0
        J_corr_st(i)=1
        Do j_corr=1,Nr_Corr(i_eo,i_chunk)
            i_ant=CorrAntNrs(j_corr,i_eo,i_chunk)
            If(Ant_Stations(i_ant,i_chunk) .ne. station) then
                i=i+1
                J_corr_st(i)=j_corr ! Number of spectra for this station
                Station=Ant_Stations(i_ant,i_chunk)
            endif
        enddo
        i_max=i+1   ! total number of stations for which spectra that will be plotted
        J_corr_st(i_max)=Nr_Corr(i_eo,i_chunk)+1
        Height=4 + (i_max+Nr_Corr(i_eo,i_chunk))/3
        Open(UNIT=9,STATUS='unknown',ACTION='WRITE',FILE=trim(XCorrPlot))
        Write(9,"(A,i3,3(/A),3(/A,I0))") 'size 63 ',Height+8,'set font pstr fontlwidth 0.08 hei 1.2 just CC',&
            'include "../Utilities/DiscreteColor.gle"',&
            'set lwidth 0.1','t_min = ',-Safety,'t_max = ',Safety,'Height = ',Height
    !        'set lwidth 0.1','t_min = 13250 !12500 !0','t_max = 13350 !15000 !',Time_dim,'scl = ',Plot_scale
        Write(9,"(A,5(/A),I3,'=',i3,':',i5.5,'(',A,')',A,6(/A))") &
            'amove 4 4','begin graph','  size 55 Height','  vscale 1','  hscale 1',&
            '   title  "Time Spectra, Peak# ',i_peak,i_chunk,Peakpos(i_peak),TRIM(ext),'"', &
            '   xtitle "time [samples]"', &! '   ytitle "Abs cross corr"',&
            '   xaxis min t_min max t_max ! nticks 10','   yaxis  min 0 max 2 dticks 5 dsubticks 1', &
            '   x2labels on',' end graph'
        Write(9,"(A,5(/A))" ) 'begin key','position tc','nobox','   offset 5. 0.0', &
            '   text "Excluded" lstyle 4 color black lwidth 0.1', &
            'end key'
        Write(9,"(A,4(/A),I2.2,A,/A)" ) 'begin key','position tc','nobox','   offset 14. 0.0', &
            '   text ">',Safety,' off" lstyle 2 color black lwidth 0.1', &
            'end key'
        Write(9,"(A,4(/A),I2.2,A,/A)" ) 'begin key','position tc','nobox','   offset 23. 0.0', &
            '   text ">',Safety/4,' off" lstyle 9 color black lwidth 0.1', &
            'end key'
        !Write(9,"(/A)" ) 'begin key','position tr','nobox', & ! '   offset -0.2 -0.5', &
        !    '   text ">Safe/4" lstyle 9 color black lwidth 0.1', &
        !    'end key'
        Write(9,"('amove xg(',f4.1,') 4',/'aline xg(',f4.1,') Height+4',/'set hei 0.7')") RefAntErr(i_Peak), RefAntErr(i_Peak) !vertical line at zero
        !
        plot_offset=4.+1./3.
        nra=.true.
        Do i=1,i_max  ! Loop over stations
            write(txt,"(i4.4,A,I3.3,A2)") Ant_Stations(CorrAntNrs(J_corr_st(i-1),i_eo,i_chunk),i_chunk),'_',i_peak,ext
            !write(2,*) i,txt,J_corr_st(i)-1-J_corr_st(i-1)
            Open(Unit=11,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//'LOFAR_Corr-'//txt//'.dat')
            write(11,*) '!',Ant_Stations(J_corr_st(i-1),i_chunk), &
                (Ant_IDs(CorrAntNrs(k,i_eo,i_chunk),i_chunk),k=J_corr_st(i-1), J_corr_st(i)-1)
            Do k=-Safety,+Safety
                write(11,*) k,(CCorr(k,pBase+J_corr_st(i-1):pBase+J_corr_st(i)-1,i_Peak))
            enddo
            close(Unit=11)
            nrb=(J_corr_st(i)-1-J_corr_st(i-1)).ge.1  !  is #antennas sufficient for label printing?
            If(i.eq.i_max) nrb=.true.
            If(nra) lr=0   ! previous #antennas sufficient for label printing
            lr=-lr
            If(.not.nrb) then   ! this #antennas is too few for station label
               If(lr.eq.0) then
                  lr=-1
               endif
            Endif
            ! write(2,*) nra,nrb,lr,(J_corr_st(i)-1-J_corr_st(i-1))
            nra=nrb
            !
            Call Station_ID2Mnem(Ant_Stations(CorrAntNrs(J_corr_st(i-1),i_eo,i_chunk),i_chunk),Station_Mnem)
            Lab_a=' '  ; If(lr.ge.0) Lab_a=Station_Mnem
            Lab_b=' '  ; If(lr.le.0) Lab_b=Station_Mnem
            Write(9,902) plot_offset,plot_offset,Lab_a
902 Format('amove 50 ',F6.2,/'aline 59 ',F6.2,/'rmove 0.5 0 ',/'write "',A5,'"')
            Write(9,903) plot_offset,Lab_b,plot_offset
903 Format('amove 2.0 ',F6.2,/'write "',A5,'"',/'amove 4 ',F6.2,/'rline 5 0')
            lwidth=0.1
            Do J_corr=J_corr_st(i-1),J_corr_st(i)-1
                !plot_offset=4.+(i+J_corr-1)/3.
                !i_c = J_corr-11*int(J_corr/11) ! Color
                k=J_corr-J_corr_st(i-1)+1
                i_c = k - 11*int(k/11) ! Color
                lstyle=0
                If(CCorr_Err(j_corr,i_Peak).ge.1900) then  ! Marks excluded station, see "ReImAtMax"
                  lstyle=4  ! 3
                ElseIf(CCorr_Err(j_corr,i_Peak).ge.190) then  ! Marks beyond Safety, see "ReImAtMax"
                   lstyle=2
                ElseIf(CCorr_Err(j_corr,i_Peak).ge.19) then  ! Marks beyond Safety/4, see "ReImAtMax"
                   lstyle=9
                ElseIf(CCorr_Err(j_corr,i_Peak).ge.10) then  ! Marks multiple peaks in correlation function, see "CrossCorr_Max"
                   lstyle=5
                Endif
                Write(9,901) plot_offset,trim(DataFolder), TRIM(txt),k,k+1,k,lwidth,lstyle,i_c
        901     Format('amove 4 ',F6.2,/'begin graph',/'  size 55 2',/'  vscale 1',&
                    /'  hscale 1',/'   NoBox',&
                    /'   xaxis min t_min max t_max',/'   x1axis off',&
                    /'   yaxis  min 0 max 1',/'   y1axis off',&
                    /'     data "',A,'LOFAR_Corr-',A,'.dat" d',I0,'=c1,c',I0,&
                    /'     d',I0,' line lwidth ',f4.2,' lstyle ',i1,' color MyCol',i0,'$',&
                    /' end graph')
                lwidth=0.05  ! Reset to normal width
                plot_offset=plot_offset+1/3.
            Enddo
                plot_offset=plot_offset+1/3.  ! creates an extra offset separating stations
        enddo   ! i=1,i_max
        !
        Close(unit=9)
         Call GLEplotControl(SpecialCmnd='gle -d pdf '//trim(XCorrPlot)) !  Command to produce curtain plots
      Enddo  ! i_tp=0,i_tp2
    enddo  !  i_peak=1,PeakNrTotal
    Return
End Subroutine GLE_Corr
!=================================
! ===========================================================================
Subroutine GLEscript_Curtains(unt, file, WWidth, i_chunk, FileA, Label, dChi_ap, dChi_at, Power_p, Power_t, Chi2pDF, VoxLoc)
!  To make curtain plots for the polarized fields for the antennas
   use constants, only : dp,pi
   !use DataConstants, only : DataFolder!, OutFileLabel
   !use Chunk_AntInfo, only : Ant_Stations, Ant_nr, Ant_IDs, Ant_pos, Ant_RawSourceDist, Nr_UniqueStat
   use Chunk_AntInfo, only : Ant_Stations, Nr_UniqueStat, Unique_StatID,  Ant_IDs, Ant_pos
   use DataConstants, only : Station_nrMax
   Use Interferom_Pars, only : IntFer_ant,  Nr_IntferCh ! the latter gives # per chunk
   use ThisSource, only : ChunkNr
   use GLEplots, only : GLEplotControl
   use StationMnemonics, only : Statn_ID2Mnem
   Implicit none
   Integer, intent(in) :: unt, WWidth, i_chunk
   Real(dp), intent(in) :: dChi_ap(*), dChi_at(*), Chi2pDF, Power_p(*), Power_t(*), VoxLoc(1:3)
   Character(Len=*), intent(in) :: FileA, Label
   Character(Len=*), intent(in) :: file
   integer :: I_ant, i_stat, Stat_ID, counter
   integer :: j_IntFer, dn
   !Character(len=41) :: txt
   Character(len=5) :: Station_Mnem
   real*8 :: dChi_sp(1:Nr_UniqueStat), dChi_st(1:Nr_UniqueStat)
   Real*8 :: StPowr_p(1:Nr_UniqueStat), StPowr_t(1:Nr_UniqueStat), phi(1:Nr_UniqueStat)
   real*8 :: plot_offset, HOffSt_p, HOffSt_t, PlotW, Sep, whiteness
   !
   !vsize=Nr_UniqueStat*2 34
   PlotW=27.
   Sep=1.
   HOffSt_p=4
   HOffSt_t=HOffSt_p+PlotW+Sep
   !
   Write(2,*) 'GLEscript_Curtains:Nr_UniqueStat=',Nr_UniqueStat
   Open(UNIT=UNT,STATUS='unknown',ACTION='WRITE',FILE=trim(file)//'.gle')
   Write(unt,"(A)") '! COMMAND:  gle -d pdf '//trim(file)//'.gle'
   Write(unt,"(A,I2,3(/A),2(/A,I0))") 'size 63 ',Nr_UniqueStat+9,'set font pstr fontlwidth 0.08 hei 1.2 just CC',&
      'include "../Utilities/DiscreteColor.gle"',&
      'set lwidth 0.1','t_min = ', -WWidth,'t_max = ', WWidth
   Write(unt,"(A,F5.1,A,2(/A),i2,3(/A),7(/A))")  'amove',HOffSt_p,' 4','begin graph', &
      '   size 27 ',Nr_UniqueStat+1,'   vscale 1','  hscale 1',&
      '   title  "E_\phi time traces, '//TRIM(Label)//'"', &
      '   xtitle "time [samples]"', '   ytitle "Amplitude"',&
      '   xaxis min t_min max t_max ! nticks 10','   yaxis  min 0 max 2 dticks 5 dsubticks 1', &
      '   x2labels on',' end graph'
   Write(unt,"(A,6(/A))" ) 'begin key','position tl','nobox', '   text \phi', &
      !         '   text "Excluded" lstyle 4 color black lwidth 0.1', &
      'end key','set just cc'
   Write(unt,"(A,/A,I2)")  'amove xg(0) 4','rline 0 ',Nr_UniqueStat+1
   Write(unt,"(A,F5.1,A,2(/A),i2,3(/A),7(/A))")  'amove',HOffSt_t,' 4','begin graph', &
      '   size 27 ',Nr_UniqueStat+1,'   vscale 1','  hscale 1',&
      '   title  "E_\theta time traces, '//TRIM(Label)//'"', &
      '   xtitle "time [samples]"',  & !'   ytitle "Amplitude"',&
      '   xaxis min t_min max t_max ! nticks 10','   yaxis  min 0 max 2 dticks 5 dsubticks 1', &
      '   x2labels on',' end graph'
   Write(unt,"(A,6(/A))" ) 'begin key','position tr','nobox', & !'   offset 20. -1.5', &
      '   text \theta', &
      !         '   text "Excluded" lstyle 4 color black lwidth 0.1', &
      'end key','set just cc'
   Write(unt,"(A,/A,I2)")  'amove xg(0) 4','rline 0 ',Nr_UniqueStat+1
   !
   !write(2,*) 'Nr_UniqueStat',Nr_UniqueStat, i_chunk, Label, Nr_IntFerCh(i_chunk)
   dChi_sp(:)=0.
   dChi_st(:)=0.
   StPowr_p(:)=0.
   StPowr_t(:)=0.
   !   write(2,*) 'dChi_ap',dChi_ap(:)
   Do i_stat=1,Nr_UniqueStat
      Stat_ID=Unique_StatID(i_stat)
      !write(2,*) 'i_stat:',i_stat, Statn_ID2Mnem(Unique_StatID(i_stat)), Unique_StatID(i_stat)
      plot_offset=i_stat + 2
      Write(Unt,901) HOffSt_p, plot_offset!,trim(FileA), TRIM(Label),j_IntFer, trim(FileA), TRIM(Label),j_IntFer, i_c, i_c ! subgraph for this antenna pair ! , ADVANCE='NO'
      dn=0
      counter=0
      Do j_IntFer=1,Nr_IntFerCh(i_chunk)
         i_ant=IntFer_ant(j_IntFer,i_chunk)
         !write(2,*) 'j_IntFer',j_IntFer,i_ant, Ant_Stations(i_ant,i_chunk), i_stat, Ant_IDs(i_ant,i_chunk)
         If(Ant_Stations(i_ant,i_chunk).ne.Stat_ID) cycle
         !write(2,*) 'j_IntFer-used',j_IntFer,i_ant, Ant_Stations(i_ant,i_chunk), Stat_ID
         dn=dn+1
         Write(Unt,902) trim(FileA)//'PhiDat_'//TRIM(Label), dn, j_IntFer+1, dn, '0', MODULO(i_ant,11)
         dn=dn+1
         Write(Unt,902) trim(FileA)//'PhiMod_'//TRIM(Label), dn, j_IntFer+1, dn, '4', MODULO(i_ant,11)
         counter=counter+1
         dChi_sp(i_stat)=dChi_sp(i_stat)*(counter-1)/counter + dChi_ap(j_IntFer)/counter  ! keep running mean
         StPowr_p(i_stat)=StPowr_p(i_stat)*(counter-1)/counter + Power_p(j_IntFer)/counter  ! keep running mean
      EndDo
      !write(2,*) 'dn:',dn,j_IntFer
      Write(Unt,"(' end graph')")
      plot_offset=i_stat + 2
      Write(Unt,901) HOffSt_t, plot_offset!,trim(FileA), TRIM(Label),j_IntFer, trim(FileA), TRIM(Label),j_IntFer, i_c, i_c ! subgraph for this antenna pair ! , ADVANCE='NO'
      dn=0
      counter=0
      Do j_IntFer=1,Nr_IntFerCh(i_chunk)
         i_ant=IntFer_ant(j_IntFer,i_chunk)
         If(Ant_Stations(I_ant,i_chunk).ne.Stat_ID) cycle
         dn=dn+1
         Write(Unt,902) trim(FileA)//'ThDat_'//TRIM(Label), dn, j_IntFer+1, dn, '0', MODULO(i_ant,11)
         dn=dn+1
         Write(Unt,902) trim(FileA)//'ThMod_'//TRIM(Label), dn, j_IntFer+1, dn, '4', MODULO(i_ant,11)
         counter=counter+1
         dChi_st(i_stat)=dChi_st(i_stat)*(counter-1)/counter + dChi_at(j_IntFer)/counter  ! keep running mean
         StPowr_t(i_stat)=StPowr_t(i_stat)*(counter-1)/counter + Power_t(j_IntFer)/counter  ! keep running mean
         Phi(i_stat)=180.+atan2( VoxLoc(2)-Ant_pos(2,i_ant,i_chunk) , VoxLoc(1)-Ant_pos(1,i_ant,i_chunk) ) *180./pi ! \phi=0 = south
      EndDo
      Write(Unt,"(' end graph')")


      !Ras(1)=(VoxLoc(1)-Ant_pos(1,i_ant,i_chunk))/1000.  ! \vec{R}_{antenna to source}
      !Ras(2)=(VoxLoc(2)-Ant_pos(2,i_ant,i_chunk))/1000.
      !Ras(3)=(VoxLoc(3)-Ant_pos(3,i_ant,i_chunk))/1000.
      !HorDist= Ras(1)*Ras(1) + Ras(2)*Ras(2)  ! =HYPOT(X,Y)
      !D=sqrt(HorDist + Ras(3)*Ras(3))
      !HorDist=sqrt( HorDist ) ! =HYPOT(X,Y)
      !Thet_r=atan(HorDist,Ras(3))  ! Zenith angle
      !Thet_d =Thet_r*180/pi
      !write(2,*) i_stat,Phi(i_stat),MOD(phi(i_stat),90.d0),i_ant,i_chunk,Ant_pos(2,i_ant,i_chunk)



   EndDo
901   Format('amove ',2F5.1,/'begin graph',/'  size 27 2',/'  vscale 1',&
         /'  hscale 1',/'   NoBox',&
         /'   xaxis min t_min max t_max',/'   x1axis off',&
         /'   yaxis ',/'   y1axis off')
!         /'   yaxis  min -scl max scl',/'   y1axis off')
902   Format('     data "',A,'.dat" d',i0,'=c1,c',I0,&
         /'     d',i0,' line lwidth 0.04 lstyle ',A,' color MyCol',i0,'$')
    !
    Write(unt,"(A,2F5.1,/A,F6.2,A)") 'amove ',HOffSt_t,Nr_UniqueStat+7.5, 'write "Chi^2/dof=',Chi2pDF,'"'
    Write(unt,"('set lstyle 0',/'set lwidth 0.01',/'set hei 0.7')")
    !write(2,*) 'Chi2pDF',Chi2pDF
    !write(2,*) 'dChi_sp(i_stat):',dChi_sp(:)
    Do i_stat=1,Nr_UniqueStat  ! put labels
!   HOffSt_t=HOffSt_p+PlotW+Sep
        plot_offset=i_stat  + 3
        Write(unt,"('amove ',2F5.1,/'rline ',F5.1' 0')") HOffSt_p, plot_offset,    PlotW  ! hline
        Write(unt,"('amove ',2F5.1,/A,A5,A)") HOffSt_p-1.5, plot_offset,'write "',Statn_ID2Mnem(Unique_StatID(i_stat)),'"'  ! text on left
        whiteness=1/(1+2*dChi_sp(i_stat)/Chi2pDF)
        Write(unt,904) HOffSt_p-4.0+PlotW, plot_offset+.5, 1.-whiteness, whiteness, dChi_sp(i_stat)            ! text in right side
        Write(unt,905) HOffSt_p-1.5+PlotW, plot_offset+.5,  sqrt(StPowr_p(i_stat))            ! text in right side
        Write(unt,"('amove ',2F5.1,/'rline ',F5.1' 0')") HOffSt_t, plot_offset,    PlotW
        Write(unt,"('amove ',2F5.1,/A,A5,A)") HOffSt_t+PlotW+1.5, plot_offset,'write "',Statn_ID2Mnem(Unique_StatID(i_stat)),'"'
        whiteness=1/(1+2*dChi_st(i_stat)/Chi2pDF)
        Write(unt,904) HOffSt_t+1.5, plot_offset+.5, 1.-whiteness, whiteness, dChi_st(i_stat)
        Write(unt,905) HOffSt_t+3.5, plot_offset+.5,  sqrt(StPowr_t(i_stat))            ! text in right side
        Write(unt,905) HOffSt_t+5.5, plot_offset+.5,  MOD(phi(i_stat),90.d0)-45.           ! text in right side
    Enddo
904   Format('amove ',2F5.1,/'set color rgb(',F3.1,',0.0,',F3.1,') ',/'write "',F5.1,'"',/'Set color black')
905   Format('amove ',2F5.1,/'write "',F5.0,'"')
!        Write(unt,903) plot_offset,plot_offset,Statn_ID2Mnem(Unique_StatID(i_stat)), &
!            plot_offset,Statn_ID2Mnem(Unique_StatID(i_stat))
!903 Format('amove 4 ',F5.2,/'aline 59',F5.1,/'rmove 1.5 0 ',/'write "',A5,'"',/'amove 2.5 ',F5.2,/'write "',A5,'"')
   !
   Close(unit=unt)
   !stop
   !
   Call GLEplotControl(SpecialCmnd='gle -d pdf '//trim(file)//'.gle') !  Command to produce curtain plots
   write(10,"('rm ',A,'.dat')") trim(FileA)//'PhiDat_'//TRIM(Label) !    remove files; clean-up
   write(10,"('rm ',A,'.dat')") trim(FileA)//'PhiMod_'//TRIM(Label)
   write(10,"('rm ',A,'.dat')") trim(FileA)//'ThDat_'//TRIM(Label)
   write(10,"('rm ',A,'.dat')") trim(FileA)//'ThMod_'//TRIM(Label)
   Call GLEplotControl(SpecialCmnd='rm '//TRIM(file)//'.gle') !  Command to delete curtain files
 !  !write(10,"('rm ',A,'{*SpecWin_*.csv,*IntfTrack_*.csv,*Interferometer*.csv,*_EISpec*.csv}')") TRIM(DataFolder) !    remove files; clean-up
   !
   return
End Subroutine GLEscript_Curtains
!====================================
