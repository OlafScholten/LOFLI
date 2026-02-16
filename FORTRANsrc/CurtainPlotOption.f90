! ========================
Subroutine PlotAllCurtainSpectra(CurtainWidth)
!  Make Curtain Plots
!  Option 'Dual' is assumed and thus only the i_eo=0 peaks are taken
   ! peaknr is included in the naming convention when plots for multple peaks are needed (=calibration runs)
   use constants, only : dp,sample
   use DataConstants, only : Time_dim, DataFolder, OutFileLabel  ! , ChunkNr_dim
   use ThisSource, only : PeakNrTotal, Peak_eo, ChunkNr, Peakpos, SourcePos
   use Chunk_AntInfo, only : CTime_spectr, Ant_Pos, Ant_nr, Ant_RawSourceDist, Unique_StatID,  Nr_UniqueStat, Ant_Stations, RefAnt
   use Chunk_AntInfo, only : CTime_Hspectr, LBA_nr, HBA_nr, HRefAnt
   use Chunk_AntInfo, only : Tot_UniqueAnt, Unique_SAI ! constructed in "Find_unique_StatAnt" in FitParams.f90
   use Chunk_AntInfo, only : Ant_IDs
   use StationMnemonics, only : Statn_ID2Mnem
   use FitParams, only : Fit_TimeOffsetStat, Fit_TimeOffsetAnt ! , N_FitPar_max, Fit_AntOffset, Fit_TimeOffsetAnt, PulsPosCore
   use GLEplots, only : GLEplotControl
   use CPU_timeUsage, only : CPU_usage
   Implicit none
   !
   integer, intent(in) :: CurtainWidth
   Real(dp) :: B, B_a(1000), RDist, OSt, Ampl_a(1000), OffSet, Ipk, Ibckg
   Real(dp) :: P_pulse(1000), P_Bckg(1000)
   Integer :: j, Lower, Upper, Ch1, Ch2, Ch1_a(1000), Ch2_a(1000), OffSet_a(1000)
   Integer :: Sfact, unt
   Integer :: i_ant, i_time, i_eo, i_chunk, i_stat, i_peak, i_type, i_HAnt, j_type
   Integer :: SAI, i_SAI, k
   Logical :: UsePeakNr
   CHARACTER(LEN=60) :: GLE_file
   Character(len=25) :: FMT=''
   CHARACTER(LEN=7) :: txt
   !
   write(*,*) 'MakeAllCurtainPlots for "'//TRIM(OutFileLabel)//'"'
   P_pulse(:)=0.
   P_Bckg(:)=0.
   Call PlotSources
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
   Call CPU_usage
   j_type=0
   Do i_peak=1,PeakNrTotal
      i_eo=Peak_eo(i_peak)
      !write(2,*) 'CurtainPlot:I_peak-A',i_peak,PeakNrTotal,I_ant,Nr_UniqueStat,i_chunk,i_eo
      If(i_eo .ne. 0) cycle !  Option 'Dual' is assumed and thus only the i_eo=0 peaks are taken
      i_chunk=ChunkNr(i_peak)
      i_ant= RefAnt(i_chunk,i_eo, j_type)
      If(i_ant.le.0) i_ant= RefAnt(i_chunk,i_eo, 1-j_type)
      !write(2,*) 'CurtainPlot:I_peak',i_peak,PeakNrTotal,I_ant,Nr_UniqueStat,i_chunk,i_eo, j_type
      !flush(unit=2)
      Do i_stat=1,Nr_UniqueStat  ! Get station number from unique list to retrieve timing-offset
         If(Unique_StatID(i_stat).eq. Ant_Stations(i_ant,i_chunk)) exit
      Enddo
      Call RelDist(SourcePos(1,i_Peak),Ant_pos(1,i_ant,i_chunk),RDist)
      OSt=(Rdist - Ant_RawSourceDist(i_ant,i_chunk) + Fit_TimeOffsetStat(i_stat)) ! to correct for off-center position reference antenna
      !setup for the different antennas
      If(Ant_nr(i_chunk).gt.1000) Then
         Write(2,*) '!!!!! Can make curtain plot only for fewer than 1000 antennas **************', i_peak,Ant_nr(i_chunk)
         Return
      EndIf
      Do I_ant=1,Ant_nr(i_chunk)
         !write(2,*) 'CurtainPlot:I_ant',I_ant,Ant_nr(i_chunk)
         !flush(unit=2)
         Do i_stat=1,Nr_UniqueStat  ! Get station number from unique list to retrieve timing-offset
            If(Unique_StatID(i_stat).eq. Ant_Stations(i_ant,i_chunk)) exit
         Enddo
         !   A=MaxLoc(Stations, mask=Stations.eq.Ant_Stations(I_ant,i_chunk))
         !   i_stat=A(1)
         !
         !write(2,*) 'SAI:',I_ant,i_stat, Tot_UniqueAnt(i_stat-1)+1,Tot_UniqueAnt(i_stat) &
         !      ,Fit_TimeOffsetAnt(i_ant), Ant_IDs(i_ant,i_chunk), Unique_StatID(i_stat)
         i_type=MOD(Unique_StatID(i_stat),10)
         SAI=100*Ant_Stations(i_ant,i_chunk) + Ant_IDs(i_ant,i_chunk)
         Do i_SAI=Tot_UniqueAnt(i_stat-1)+1,Tot_UniqueAnt(i_stat)      ! Get antenna number from the Unique_Antennas list
            If(SAI .eq. Unique_SAI(i_SAI)) Then !     found the i_SAI for this antenna, needed for the delay
               exit
            EndIf
         EndDo
         !write(2,*) 'SAI_AntDelay(:', k, SAI, i_SAI, Unique_SAI(i_SAI), Fit_TimeOffsetAnt(i_SAI)
         !
         Call RelDist(SourcePos(1,i_Peak),Ant_pos(1,i_ant,i_chunk),RDist)
         Sfact=1+i_type
         If(i_type.eq.1) Then
            i_HAnt=i_ant-LBA_nr(i_chunk)
         Endif
         OffSet=(Rdist - Ant_RawSourceDist(i_ant,i_chunk) + Fit_TimeOffsetStat(i_stat) +Fit_TimeOffsetAnt(i_SAI) -OSt)
         !write(2,*) Ant_Stations(i_ant,i_chunk), 'curtainoffset:',OffSet, Rdist, Ant_RawSourceDist(i_ant,i_chunk), &
         !      Fit_TimeOffsetStat(i_stat), Fit_TimeOffsetAnt(i_SAI), OSt
         Ch2=NINT(Sfact*(PeakPos(i_Peak) + OffSet))
         Lower=MIN(Sfact*CurtainWidth, Ch2-1)
         Upper=MIN(Sfact*CurtainWidth,(Sfact*Time_dim - Ch2))
         Ch1=Ch2-Lower
         Ch2=Ch2+Upper
         !write(2,*) I_ant,i_peak, Ch1, Ch2, Lower, Upper,Rdist, Ant_RawSourceDist(i_ant,i_chunk), Fit_TimeOffsetStat(i_stat), OSt
         If(Upper.le.0) then
            Ch1=Time_dim/2
            Ch2=Ch1+2
         Endif
         If(Ch2.lt.Ch1) Ch2=Ch1+2
         !write(2,*) 'working on curtainplot: "','LOFAR_Time'//TRIM(txt),'"'
         !
         If(i_type.eq.0) Then
            B=maxval(real(CTime_spectr(ch1:ch2,i_Ant,i_chunk)))
            i_time=MaxLoc(abs(CTime_spectr(ch1:ch2,i_Ant,i_chunk)), DIM=1)+ch1-1
            Lower=MAX(ch1,i_time-10)
            Upper=MIN(ch2,i_time+10)
            j=Upper-Lower
            Ipk=SUM(abs(CTime_spectr(Lower:Upper,i_Ant,i_chunk))**2)
            If((ch2-ch1) .le. 22) Then
               Ibckg=Ipk/j
            Else
               Ibckg=(SUM(abs(CTime_spectr(ch1:ch2,i_Ant,i_chunk))**2)-Ipk)/(ch2-ch1-j)
            EndIf
            Ipk=Ipk/j
         ElseIf(i_type.eq.1) Then
            B=maxval(abs(CTime_Hspectr(ch1:ch2,i_HAnt,i_chunk)))
            i_time=MaxLoc(abs(CTime_Hspectr(ch1:ch2,i_HAnt,i_chunk)), DIM=1)+ch1-1
            Lower=MAX(ch1,i_time-20)
            Upper=MIN(ch2,i_time+20)
            j=Upper-Lower
            Ipk=SUM(abs(CTime_Hspectr(Lower:Upper,i_HAnt,i_chunk))**2)
            If((ch2-ch1) .le. 44) Then
               Ibckg=Ipk/j
            Else
               Ibckg=(SUM(abs(CTime_Hspectr(ch1:ch2,i_HAnt,i_chunk))**2)-Ipk)/(ch2-ch1-j)
            EndIf
            Ipk=Ipk/j
            !write(2,*) ch1, ch2, i_time, Lower, Upper, Ipk, Ibckg
         EndIf
         P_pulse(i_Ant)=sqrt(Ipk)
         P_Bckg(i_Ant)=sqrt(Ibckg)
         ! Need to store:
         Ch1_a(i_Ant)=Ch1
         Ch2_a(i_Ant)=Ch2
         Offset_a(i_Ant)=NINT(Sfact*Offset)
         B_a(i_Ant)=B
      EndDo ! I_ant=1,Ant_nr(i_chunk)
      !stop
      If(UsePeakNr) Then  !  "'//TRIM(OutFileLabel)//'"
         write(txt,"('-',i3.3)") i_peak
      Else
         txt='X'
      EndIf
      !
      !write(2,*) '!PlotAllCurtainSpectra:',i_chunk, LBA_nr(i_chunk), HBA_nr(i_chunk)
      IF(LBA_nr(i_chunk).gt.0) Then
         Open(Unit=9,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//'LOFAR_Time'//TRIM(OutFileLabel)//TRIM(txt)//'L.dat')
         write(FMT,"(A,I3,A)") '(I6,1x,',Ant_nr(i_chunk),'(F6.3,1x) )'
         Do j=PeakPos(i_Peak)-CurtainWidth,PeakPos(i_Peak)+CurtainWidth
            Do I_ant=1,LBA_nr(i_chunk)
               i_time=j+Offset_a(I_ant)
               If((i_time.ge.Ch1_a(i_Ant)) .and. (i_time.le.Ch2_a(i_Ant)) ) Then
                  Ampl_a(I_ant)= Real(CTime_spectr(i_time,i_Ant,i_chunk))/B_a(i_Ant)
               Else
                  Ampl_a(I_ant)=0.
               EndIf
            EndDo
            Write(9,FMT) j, Ampl_a(1:LBA_nr(i_chunk)) ! written time corresponds to time at core for this source
         Enddo
         close(unit=9)
      EndIf
      !
      IF(HBA_nr(i_chunk).gt.0) Then
         Open(Unit=9,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//'LOFAR_Time'//TRIM(OutFileLabel)//TRIM(txt)//'H.dat')
         write(FMT,"(A,I3,A)") '(I6,1x,',Ant_nr(i_chunk),'(F6.3,1x) )'
         Sfact=2
         Do j=Sfact*PeakPos(i_Peak)-Sfact*CurtainWidth,Sfact*PeakPos(i_Peak)+Sfact*CurtainWidth
            Do i_HAnt=1,HBA_nr(i_chunk)
               i_ant=LBA_nr(i_chunk)+i_HAnt
               i_time=j+Offset_a(I_ant)
               If((i_time.ge.Ch1_a(i_Ant)) .and. (i_time.le.Ch2_a(i_Ant)) ) Then
                  !Ampl_a(I_ant)= Real(CTime_Hspectr(i_time,i_HAnt,i_chunk))/B_a(i_Ant)
                  Ampl_a(I_ant)= ABS(CTime_Hspectr(i_time,i_HAnt,i_chunk))/B_a(i_Ant) ! since real part oscillates too much for HBA
               Else
                  Ampl_a(I_ant)=0.
               EndIf
            EndDo
            Write(9,FMT) j, Ampl_a(LBA_nr(i_chunk)+1:Ant_nr(i_chunk)) ! written time corresponds to time at core for this source
         Enddo
         close(unit=9)
      EndIf
      !
      unt=9
      write(2,*) 'SourcePos(:,i_Peak)',SourcePos(:,i_Peak)
      If(UsePeakNr) Then  !  "'//TRIM(OutFileLabel)//'"
         Write(GLE_file,"('CuP-',A,i3.3)") TRIM(OutFileLabel),i_peak
      Else
         Write(GLE_file,"('CuP-',A)") TRIM(OutFileLabel)
      EndIf
      Write(2,*) 'GLE-file name:',TRIM(GLE_file)
      Call GLEscript_CurtainPlot(unt, GLE_file, CurtainWidth, i_Peak, UsePeakNr)
      !
      Call GLEscript_AntCurtainPlot(unt, GLE_file, CurtainWidth, i_Peak, UsePeakNr, B_a, P_pulse, P_Bckg)
      !
   Enddo ! i_peak
   !
   !stop '!PlotAllCurtainSpectra'
   !Call GLEplotControl(SpecialCmnd='rm '//TRIM(DataFolder)//'LOFAR_Time'//TRIM(OutFileLabel)//'*.dat') !  Command to delete curtain files
   Return
End Subroutine PlotAllCurtainSpectra
! ===========================================================================
Subroutine PlotsingleCurtainSpectrum(i_sample, Location)
! Needs first a call to  Find_unique_StatAnt()
   use constants, only : dp
   use ThisSource, only : PeakNrTotal, Peak_eo, ChunkNr, Peakpos, SourcePos
   use ThisSource, only : CurtainHalfWidth
   Implicit none
   Integer, intent(in) :: i_sample
   Real(dp), intent(in) :: Location(1:3)
   PeakNrTotal=1
   Peak_eo(1)=0     ! only the i_eo=0 peaks are used for curtainplotting
   ChunkNr(1)=1   ! chunk nr
   !RefAnt(1,0,0:1)=1
   !RefAnt(i_chunk,i_eo, i_type)
   PeakPos(1)=i_sample   ! take a single peak at the center of the window
   SourcePos(1:3,1)=Location(1:3)
   Call PlotAllCurtainSpectra(CurtainHalfWidth)  ! uses TotPeakNr(0,i_chunk), PeakNr(0,i_chunk)
   Return
End Subroutine PlotsingleCurtainSpectrum
! ===========================================================================
Subroutine GLEscript_CurtainPlot(unt, file, CurtainWidth, i_Peak, UsePeakNr)
   use constants, only : dp
    use DataConstants, only : Time_dim, DataFolder, OutFileLabel
    use Chunk_AntInfo, only : Ant_Stations, Ant_nr, Ant_IDs, Ant_pos, Ant_RawSourceDist, Nr_UniqueStat, Unique_StatID
    use DataConstants, only : Used_StationNr
    use StationMnemonics, only : Station_ID2Mnem
    use ThisSource, only : PeakPos, ChunkNr, ExclStatNr
   use ThisSource, only : PeakNrTotal, Dropped, PeakChiSQ, PeakRMS
    use GLEplots, only : GLEplotControl
    Implicit none
    Integer, intent(in) :: unt, CurtainWidth, i_Peak
    Logical, intent(in) :: UsePeakNr
    Character(Len=*), intent(in) :: file
    integer :: I_ant, i_stat, Stat_nr, Stations(1:Used_StationNr), i_time, i_c, lstyle
    integer :: A(1), Plot_scale, ch1, ch2, VLine, i_chunk, k
   Integer ::  i_src, i_Peak_o, i_eo, k_peak, i_type
    Character(len=45) :: txt
    Character(len=6) :: Station_Mnem
    real*8 :: plot_offset,b
   Integer, external :: SourceNr  ! source code in  LOFLI_InputHandling.f90
    !
    i_chunk=ChunkNr(i_peak)
    Plot_scale= 2. ! int(B) + 2.
    VLine=PeakPos(i_Peak)
    ch1=PeakPos(i_Peak)-CurtainWidth
    ch2=PeakPos(i_Peak)+CurtainWidth
    !vsize=Nr_UniqueStat*2 34
   !
   ! Check if there is another peak (with other i_eo) is associated with this source
   i_src=SourceNr(i_Peak)
   i_Peak_o=i_Peak
   If(i_src.gt.0) then
      Do k=i_Peak+1,PeakNrTotal
         If(i_src .eq. SourceNr(k) ) Then
            i_Peak_o=k
            Exit
         EndIf
      EndDo
   end if
   !write(2,*) PeakNrTotal, i_src, i_peak, i_Peak_o
   !write(2,"(40I3)") Dropped(1:Nr_UniqueStat,i_peak)
   !write(2,"(40I3)") Dropped(1:Nr_UniqueStat,i_peak_o)
    !
    !write(2,*) 'GLEscript_CurtainPlot:Nr_UniqueStat=', Nr_UniqueStat
    Open(UNIT=UNT,STATUS='unknown',ACTION='WRITE',FILE=trim(file)//'.gle')
    Write(unt,"(A)") '! COMMAND:  gle -d pdf '//trim(file)//'.gle'
    ! make sure:  export GLE_USRLIB=/Users/users/scholten/LOFLI/GLEsrc
    Write(unt,"(A,I2,3(/A),3(/A,I0))") 'size 64 ',2*Nr_UniqueStat+9,'set font pstr fontlwidth 0.08 hei 1.2 just CC',&
        'include "DiscreteColor.gle"',&
        'set lwidth 0.1','t_min = ',ch1,'t_max = ',ch2,'scl = ',Plot_scale
!        'set lwidth 0.1','t_min = 13250 !12500 !0','t_max = 13350 !15000 !',Time_dim,'scl = ',Plot_scale
     If(UsePeakNr) Then  !  "'//TRIM(OutFileLabel)//'"
         write(txt,"(A,i2,A,i3.3,A,2F5.1)") 'Chunk=',modulo(i_chunk,100),', peak=',i_peak,', RMS=',PeakRMS(i_Peak),PeakRMS(i_Peak_o)
     Else
         write(txt,"(A)") TRIM(OutFileLabel)
     EndIf
    Write(unt,"(A,2(/A),i2,3(/A),7(/A))")  'amove 4 4','begin graph', &
        '   size 55 ',2*Nr_UniqueStat+1,'   vscale 1','  hscale 1',&
        '   title  "Time Traces, '//TRIM(txt)//'"', &
        '   xtitle "time [5 ns] with shifts rounded off"','   ytitle "Amplitude"',&
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
            i_c=-1   ! initialize coloring scheme for every station
            !write(2,*) 'stat_nr=',Stat_nr
        else
            A=MaxLoc(Stations, mask=Stations.eq.Ant_Stations(I_ant,i_chunk))
            i_stat=A(1)
            !write(2,*) 'i_stat=',i_Stat,Stations(1:Stat_nr)
        endif
        !write(2,"(A, 40I3)") 'i_stat=',i_Stat,Stat_nr, Ant_Stations(I_ant,i_chunk), Stations(1:Stat_nr)
        If(UsePeakNr) Then  !  "'//TRIM(OutFileLabel)//'"
            write(txt,"(A,'-',i3.3)") TRIM(OutFileLabel),i_peak
        Else
            write(txt,"(A,'X')") TRIM(OutFileLabel)
        EndIf
        !write(txt,"(i3.3,'-'i2.2)") I_ant,i_peak
        i_eo=mod(Ant_IDs(I_ant,i_chunk),2)
        plot_offset=2*i_stat + 2 + i_eo
!        i_c = I_Ant-11*int(I_Ant/11)      ! color
        i_c = MODULO(i_c+1,11)      ! color
        i_type=MOD(Unique_StatID(i_stat),10)
        !
        lstyle=0
        k_peak=i_peak
        If(i_eo.eq.1) k_peak=i_peak_o
         Do i_stat=1, Nr_UniqueStat      ! Get station number from the Unique_StatID list
             If(Unique_StatID(i_stat).eq. Ant_Stations(I_ant,i_chunk)) Then
               If(Dropped(i_stat,k_peak) .gt. 0)  lstyle=4  !  ! Marks excluded station, see "ReImAtMax"
        !write(2,"(A, 40I5)") 'i_stat=',i_Stat, Ant_Stations(I_ant,i_chunk), i_eo, Dropped(i_stat,k_peak)
               exit
             EndIf
         Enddo
        !
         If(i_type.eq.0) Then
            Write(Unt,901) plot_offset,trim(DataFolder), TRIM(txt),I_Ant+1,lstyle,i_c
         ElseIf(i_type.eq.1) Then
            Write(Unt,903) plot_offset,trim(DataFolder), TRIM(txt),I_Ant+1,lstyle,i_c
         EndIf
901     Format('amove 4 ',F5.2,/'begin graph',/'  size 55 2',/'  vscale 1',&
            /'  hscale 1',/'   NoBox',&
            /'   xaxis min t_min max t_max',/'   x1axis off',&
            /'   yaxis  min -scl max scl',/'   y1axis off',&
            /'     data "',A,'LOFAR_Time',A,'L.dat" d1=c1,c',I0, &
            /'     d1 line lwidth 0.04 lstyle ',i1,' color MyCol',i0,'$',&
            /' end graph')
903     Format('amove 4 ',F5.2,/'begin graph',/'  size 55 2',/'  vscale 1',&
            /'  hscale 1',/'   NoBox',&
            /'   xaxis min 2*t_min max 2*t_max',/'   x1axis off',&
            /'   yaxis  min -scl max scl',/'   y1axis off',&
            /'     data "',A,'LOFAR_Time',A,'H.dat" d1=c1,c',I0, &
            /'     d1 line lwidth 0.04 lstyle ',i1,' color MyCol',i0,'$',&
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
902 Format('amove 4 ',F5.2,/'aline 59 ',F5.2,/'rmove 2 0 ',/'write "',A6,'"')
    Enddo
    !
    Close(unit=unt)
    !
    Call GLEplotControl(SpecialCmnd='gle -d pdf '//trim(file)//'.gle') !  Command to produce curtain plots
    Call GLEplotControl(SpecialCmnd='rm '//TRIM(file)//'.gle') !  Command to delete curtain files
    !
    return
End Subroutine GLEscript_CurtainPlot
! ===========================================================================
Subroutine GLEscript_AntCurtainPlot(unt, file, CurtainWidth, i_Peak, UsePeakNr, B_a, P_pulse, P_Bckg)
!  Plot the traces for the different antennas separated from eachother
   use constants, only : dp
    use DataConstants, only : Time_dim, DataFolder, OutFileLabel
    use Chunk_AntInfo, only : Ant_Stations, Ant_nr, Ant_IDs, Ant_pos, Ant_RawSourceDist
    use Chunk_AntInfo, only :  Nr_UniqueStat, Unique_StatID,Tot_UniqueAnt
    use DataConstants, only : Used_StationNr
    use StationMnemonics, only : Station_ID2Mnem
    use ThisSource, only : PeakPos, ChunkNr, ExclStatNr, Error_norm, CCorr_Err, CorrAntNrs   !(j_corr,i_Peak), Ant_IDs(CorrAntNrs(j_corr,i_eo,i_chunk),i_chunk)
   use ThisSource, only : PeakNrTotal, Dropped, PeakChiSQ, PeakRMS
    use GLEplots, only : GLEplotControl
    Implicit none
    Integer, intent(in) :: unt, CurtainWidth, i_Peak
    Logical, intent(in) :: UsePeakNr
    Character(Len=*), intent(in) :: file
    Real(dp), intent(in) :: B_a(*), P_pulse(*), P_Bckg(*)
    integer :: I_ant, i_stat, Stat_nr, Stations(1:Used_StationNr), MinLocAnt(0:Used_StationNr), i_time, i_dip, lstyle
    integer :: A(1), ch1, ch2, VLine, i_chunk, k, NrAntStat, MaxAntStat, NrOdd,V_offset, H_off, Layer, i_c
   Integer ::  i_src, i_Peak_o, i_eo, k_peak, i_type, j_corr0, j_corr1, j_corr, i_stUq(1:Nr_UniqueStat)
   Integer, parameter :: Yoff0=4
   Integer :: Ysize
   Integer, parameter :: Xside=5 ! X-width side pane
   Integer, parameter :: Xtime=55
   Integer, parameter :: Xcent=5
   Integer, parameter :: Xoffs1 =4   ! X-offset first side pane
   Integer, parameter :: Xofft1=Xoffs1 +Xside  ! X-offset first timetrace
   Integer, parameter :: Xoffc =Xofft1 +Xtime  ! X-offset central pane
   Integer, parameter :: Xofft2=Xoffc  +Xcent  ! X-offset second timetrace
   Integer, parameter :: Xoffs2=Xofft2+ Xtime   ! X-offset second side pane
   Integer, parameter :: Ysubsize=6
   real*8 :: Plot_scale, plot_offset,b, eoRatMax, SRNMax
   Logical :: AntPair
    Character(len=45) :: txt
    Character(len=6) :: Station_Mnem, color
   Integer, external :: SourceNr  ! source code in  LOFLI_InputHandling.f90
   !
   i_chunk=ChunkNr(i_peak)
   Plot_scale= 1.0 ! int(B) + 2.
   VLine=PeakPos(i_Peak)
   ch1=PeakPos(i_Peak)-CurtainWidth
   ch2=PeakPos(i_Peak)+CurtainWidth
   !vsize=Nr_UniqueStat*2 34
   !
   ! Check if there is another peak (with other i_eo) is associated with this source
   i_src=SourceNr(i_Peak)
   i_Peak_o=i_Peak
   If(i_src.gt.0) then
      Do k=i_Peak+1,PeakNrTotal
         If(i_src .eq. SourceNr(k) ) Then
            i_Peak_o=k
            Exit
         EndIf
      EndDo
   end if
   !
   ! Get sizes of different sections
   !
   i_stat=0
   Stations(:)=0
   MaxAntStat=0
   Layer=0
   i_dip=-1
   eoRatMax=0.
   SRNMax=10.  ! in order to have a reasonable lower limit
   MinLocAnt(0)=0
   Do I_ant=1,Ant_nr(i_chunk)
      k=MinLoc(Ant_Stations(:,i_chunk), mask=Ant_Stations(:,i_chunk).eq.Ant_Stations(I_ant,i_chunk), Dim=1)
      !write(2,*) I_ant,Ant_Stations(I_ant),A(1)
      If(k.eq.I_ant) then
         MinLocAnt(i_stat)=i_ant
         i_stat=i_stat+1
         ! i_stCh is the sequential station number for the stations in this chunk
         ! i_stUq corresponds to the sequential station number in Unique_StatID, i.e. all used stations for several chunks
         ! here i_stat corresponds to i_stCh
         Stations(i_stat)=Ant_Stations(I_ant,i_chunk)
         k=Tot_UniqueAnt(i_stat)    ! highest ranknr in Unique_SAI of unique_antenna for station Unique_StatID(i_stat)
         NrAntStat=k-MaxAntStat
         MaxAntStat=k
         NrOdd=0
         Do k=I_ant,MaxAntStat
            If(mod(Ant_IDs(k,i_chunk),2) .eq. 1)  NrOdd=NrOdd+1
         EndDo
         Do k=1,Nr_UniqueStat
            If( Unique_StatID(k) .ne. Stations(i_stat)) cycle
            i_stUq(i_stat)=k
            !write(2,*) '!GLEscript_AntCurtainPlot, stat_nr=', i_stat,k, Stations(i_stat)
            exit
 !GLEscript_AntCurtainPlot          26          22        1300        1300           1           0 RS210L
 !GLEscript_AntCurtainPlot          26          22        1300        1300           3           2 RS210L
         EndDo
         !Call Station_ID2Mnem(Stations(I_stat),Station_Mnem)
         !write(2,*) '!GLEscript_AntCurtainPlot, stat_nr=',Stat_nr,Station_Mnem, MaxAntStat, NrAntStat,NrOdd, &
         !   Ant_IDs(i_ant,i_chunk),Ant_IDs(MaxAntStat,i_chunk), layer
      EndIf
      SRNMax=MAX(SRNMax,P_pulse(i_ant)/P_Bckg(i_ant))
      !If(SRNMax.gt.1000.) write(2,*) '!SNRMax=',i_ant,SRNMax,P_pulse(i_ant),P_Bckg(i_ant)
      If(Ant_IDs(I_ant,i_chunk)/2 .ne. i_dip) Then
         Layer=Layer+1
         i_dip=Ant_IDs(I_ant,i_chunk)/2
      Else  ! other dipole at the same location as the previous, check polarization health
         i_dip=-1
         eoRatMax=MAX(eoRatMax,ABS(P_pulse(i_ant-1)/P_pulse(i_ant)-P_pulse(i_ant)/P_pulse(i_ant-1)), &
            ABS(B_a(i_ant-1)/B_a(i_ant)-B_a(i_ant)/B_a(i_ant-1)))
         !write(2,"(I3, 2F6.2, 2F6.1,' ;')",  ADVANCE='NO') Ant_IDs(I_ant,i_chunk)-1, B_a(i_ant-1)/B_a(i_ant), &
         !   P_pulse(i_ant-1)/P_pulse(i_ant), P_pulse(i_ant-1)/P_Bckg(i_ant-1), P_pulse(i_ant)/P_Bckg(i_ant)
      EndIf
   EndDo
   Ysize=Layer+1
   !
   !write(2,*) 'GLEscript_CurtainPlot:Nr_UniqueStat=', Nr_UniqueStat
   Open(UNIT=UNT,STATUS='unknown',ACTION='WRITE',FILE=trim(file)//'Ant.gle')
   Write(unt,"(A)") '! COMMAND:  gle -d pdf '//trim(file)//'Ant.gle'
   ! make sure:  export GLE_USRLIB=/Users/users/scholten/LOFLI/GLEsrc
   !Write(unt,"(A,I3,3(/A),2(/A,I0),(/A,F4.1))") 'size 117 ',Layer+9,'set font pstr fontlwidth 0.08 hei 1.2 just CC',&
   !Write(unt,"(A,I3,3(/A),2(/A,I0),(/A,F4.1))") 'size 150 ',Layer+9,'set font pstr fontlwidth 0.08 hei 1.2 just CC',&
   !  'include "DiscreteColor.gle"',&
   !  'set lwidth 0.2','t_min = ',ch1,'t_max = ',ch2,'scl = ',Plot_scale
   Write(unt,"('size ',I3,I4)") Xoffs2+Xside+3,Ysize+Yoff0+4
   Write(unt,"(3(/A))") 'set font pstr fontlwidth 0.08 hei 1.2 just CC',&
     'include "DiscreteColor.gle"',&
     'set lwidth 0.2'
   Write(unt,"(2(/A,I0),(/A,F4.1))") 't_min = ',ch1,'t_max = ',ch2,'scl = ',Plot_scale
   Write(unt,"(A,I0)") 'Xside= ',Xside ! X-width side pane
   Write(unt,"(A,I0)") 'Xtime= ',Xtime
   Write(unt,"(A,I0)") 'Xcent= ',Xcent
   Write(unt,"(A,I0)") 'Xoffs1 = ',Xoffs1   ! X-offset first side pane
   Write(unt,"(A,I0)") 'Xofft1= ',Xofft1 ! X-offset first timetrace
   Write(unt,"(A,I0)") 'Xoffc = ',Xoffc   ! X-offset central pane
   Write(unt,"(A,I0)") 'Xofft2= ',Xofft2  ! X-offset second timetrace
   Write(unt,"(A,I0)") 'Xoffs2= ',Xoffs2   ! X-offset second side pane
   Write(unt,"(A,I0)") 'Yoff0= ',Yoff0 ! Y-offset base
   Write(unt,"(A,I0)") 'Ysize= ',Ysize ! Y-heights main plots
   Write(unt,"(A,I0)") 'Ysubsize= ',Ysubsize ! Y-heights of sub-plots of the time traces
   !
   If(UsePeakNr) Then  !  "'//TRIM(OutFileLabel)//'"
      write(txt,"(A,i2,A,i3.3,A,2F5.1)") 'Chunk=',modulo(i_chunk,100),', peak=',i_peak,' even, RMS=',PeakRMS(i_Peak)
   Else
      write(txt,"(A)") TRIM(OutFileLabel)
   EndIf
   Write(unt,"(A,4(/A),(/A,I6,A,I6),5(/A))")  'amove Xoffs1 Yoff0','begin graph', & ! first side pane
     '   size Xside Ysize','   vscale 1','  hscale 1','   xaxis min 0 max ',INT(SRNMax)+1,' dticks ',Int(SRNMax/3)+1,&
     '   x2axis min 0 max 5 dticks 5. dsubticks 1. nofirst nolast','   yaxis  nticks 1','   ylabels off ', &
     '   xtitle "SNR"', ' end graph'
   Write(unt,"(A,(/A))")  'amove Xoffs1+Xside/5 Yoff0', 'set lwidth 0.05','rline 0 Ysize', 'set lwidth 0.2'
   Write(unt,"(A,3(/A))")  'amove Xoffs1 Yoff0+Ysize+2.5', 'marker cross 0.5 ', 'rmove 2 0', 'write SNR'
   Write(unt,"(A,3(/A))")  'amove Xoffs1 Yoff0+Ysize+1.5', 'marker dot 0.5 ', 'rmove 2 0', 'write Err'
   !
   Write(unt,"(A,2(/A),2(/A),8(/A))")  'amove Xofft1 Yoff0','begin graph', & ! first time trace pane
     '   size Xtime Ysize','   vscale 1','  hscale 1',&
     '   title  "Time Traces, '//TRIM(txt)//'"', &
     '   xtitle "time [5 ns] with shifts rounded off"', &
     '   xaxis min t_min max t_max ! nticks 10','   yaxis  min 0 max 2 dticks 5 dsubticks 1', &
     '   x2labels on','   ylabels off ',' end graph'
   Write(unt,"(A,i5.5,A,1(/A))")  'amove xg(',Vline,') Yoff0','rline 0 Ysize'
   Write(unt,"(A,i5,3(/A))") 'amove Xofft1+2 ',Yoff0+Ysize+3, 'write Ampl','amove Xofft1+Xtime-10 0','write Antenna'
   !
   Write(unt,"(A,4(/A),(/A,I6,A,I5,A,I5),2(/A))")  'amove Xoffc Yoff0','begin graph', & ! central pane for e/o ratio
     '   size Xcent Ysize','   vscale 1','  hscale 1', &
     '   xaxis min ',-(Int(eoRatMax)+1),' max ',Int(eoRatMax)+1,' dticks ',Int(eoRatMax/2)+1,'   ylabels off ', ' end graph'
   Write(unt,"(A,3(/A))")  'amove Xoffc+Xcent/2 Yoff0+Ysize+2', 'write "Even/Odd"'
   Write(unt,"(A,(/A))")  'amove Xoffc+Xcent/2 Yoff0', 'set lwidth 0.05','rline 0 Ysize', 'set lwidth 0.2'
   !
   If(UsePeakNr) write(txt,"(A,i2,A,i3.3,A,2F5.1)") 'Chunk=',modulo(i_chunk,100),', peak=',i_peak,' odd, RMS=',PeakRMS(i_Peak_o)
   Write(unt,"(A,2(/A),3(/A),8(/A))")  'amove Xofft2 Yoff0','begin graph', &   ! second time trace pane
     '   size Xtime Ysize','   vscale 1','  hscale 1',&
     '   title  "Time Traces, '//TRIM(txt)//'"', &
     '   xtitle "time [5 ns] with shifts rounded off"',&
     '   xaxis min t_min max t_max ! nticks 10','   yaxis  min 0 max 2 dticks 5 dsubticks 1', &
     '   x2labels on','   ylabels off ',' end graph'
   Write(unt,"(A,i5.5,A,1(/A))")  'amove xg(',Vline,') Yoff0','rline 0 Ysize'
   Write(unt,"(A,6(/A))" ) 'begin key','position tc','nobox','   offset 20. -1.5', &
      '   text "Excluded" lstyle 4 color black lwidth 0.1', &
      'end key','set just rb'
   Write(unt,"(A,4(/A),(/A,I6,A,I6),6(/A))")  'amove Xoffs2 Yoff0','begin graph', & ! second side pane
     '   size Xside Ysize','   vscale 1','  hscale 1','   xaxis min 0 max ',INT(SRNMax)+1,' dticks ',Int(SRNMax/3)+1,&
     '   x2axis min 0 max 5 dticks 5. dsubticks 1.','   yaxis  nticks 1','   ylabels off ', &
     '   xtitle "SNR"','   x2title "Err[ns]"', ' end graph'
   Write(unt,"(A,(/A))")  'amove Xoffs2+Xside/5 Yoff0', 'set lwidth 0.05','rline 0 Ysize', 'set lwidth 0.2'
   !
   ! Plot traces for the different antennas
   !
   i_stat=0
   V_offset=3 ! 2
   Layer=0
   j_corr0=0
   j_corr1=0
   Do I_ant=1,Ant_nr(i_chunk)
      !k=MinLoc(Ant_Stations(:,i_chunk), mask=Ant_Stations(:,i_chunk).eq.Ant_Stations(I_ant,i_chunk), Dim=1)
      !write(2,*) I_ant,Ant_Stations(I_ant),A(1)
      If(MinLocAnt(i_stat).eq.I_ant) then  ! Deal with the next station in the list
         i_stat=i_stat+1
         Call Station_ID2Mnem(Stations(I_stat),Station_Mnem)
         i_dip=-1 ! to make sure layer is increased
         Write(unt,902) MODULO(Layer+1,10)+1, Xoffs1, V_offset+Layer+1, Station_Mnem
      endif
902 Format('set color MyCol',i0,'$'/'amove ',I0,' ',I0,/'write "',A6,'"')
      !If(i_dip.eq. -1) write(2,"(/,'Pol ratio',A,':')",  ADVANCE='NO') Station_Mnem
      If(Ant_IDs(I_ant,i_chunk)/2 .ne. i_dip) Then
         Layer=Layer+1
         AntPair=.false.
      Else  ! other dipole at the same location as the previous, check polarization health
         AntPair=.true.
         !write(2,"(I3, 2F6.2, 2F6.1,' ;')",  ADVANCE='NO') Ant_IDs(I_ant,i_chunk)-1, B_a(i_ant-1)/B_a(i_ant), &
         !   P_pulse(i_ant-1)/P_pulse(i_ant), P_pulse(i_ant-1)/P_Bckg(i_ant-1), P_pulse(i_ant)/P_Bckg(i_ant)
      EndIf
      i_dip=Ant_IDs(I_ant,i_chunk)/2
      plot_offset=V_offset+Layer
      !write(2,"(A, 40I3)") 'i_stat=',i_Stat,Stat_nr, Ant_Stations(I_ant,i_chunk), Stations(1:Stat_nr)
      If(UsePeakNr) Then  !  "'//TRIM(OutFileLabel)//'"
         write(txt,"(A,'-',i3.3)") TRIM(OutFileLabel),i_peak
      Else
         write(txt,"(A,'X')") TRIM(OutFileLabel)
      EndIf
      !write(txt,"(i3.3,'-'i2.2)") I_ant,i_peak
      i_eo=mod(Ant_IDs(I_ant,i_chunk),2)
      !color='red'
      If(i_eo.eq.0) Then
         !color='blue'
         j_corr0=j_corr0+1
         j_corr=j_corr0
      Else
         j_corr1=j_corr1+1
         j_corr=j_corr1
      EndIf
      If(CorrAntNrs(j_corr,i_eo,i_chunk) .ne. i_ant) write(2,*) '!GLEscript_AntCurtainPlot, antenna', j_corr, i_eo, &
         i_ant, CorrAntNrs(j_corr,i_eo,i_chunk),Ant_IDs(i_ant,i_chunk),Ant_IDs(CorrAntNrs(j_corr,i_eo,i_chunk),i_chunk)
      !        i_c = I_Ant-11*int(I_Ant/11)      ! color
      i_type=MOD(Unique_StatID(i_stat),10)
      !
      lstyle=0
      k_peak=i_peak
      If(i_eo.eq.1) k_peak=i_peak_o
      If(Dropped(i_stUq(i_stat),k_peak) .gt. 0)  lstyle=4  !  ! Marks excluded station, see "ReImAtMax"
      ! i_stCh is the sequential station number for the stations in this chunk
      ! i_stUq corresponds to the sequential station number in Unique_StatID, i.e. all used stations for several chunks
      ! here i_stat corresponds to i_stCh
      !write(2,*) '!GLEscript_AntCurtainPlot', i_stat,i_stUq(i_stat), Unique_StatID(i_stUq(i_stat)),  &
      !   Stations(I_stat), k_peak,Dropped(i_stUq(i_stat),k_peak), Station_Mnem
 !GLEscript_AntCurtainPlot          26          22        1300        1300           1           0 RS210L
 !GLEscript_AntCurtainPlot          26          22        1300        1300           3           2 RS210L
      !
      i_c=MODULO(Layer,10)+1
      If(i_type.eq.0) Then
         Write(Unt,901) i_eo+1,plot_offset-Ysubsize/2., i_c, trim(DataFolder), TRIM(txt),I_Ant+1,lstyle,i_c, &
            Stations(I_stat)*100+ Ant_IDs(i_ant,i_chunk), B_a(i_ant)
      ElseIf(i_type.eq.1) Then
         Write(Unt,903) i_eo+1,plot_offset-Ysubsize/2., i_c, trim(DataFolder), TRIM(txt),I_Ant+1,lstyle, i_c, &
            Stations(I_stat)*100+ Ant_IDs(i_ant,i_chunk), B_a(i_ant)
         !Write(2,*) '!GLEscript_AntCurtainPlot, H:',H_off,plot_offset,trim(DataFolder), TRIM(txt),I_Ant+1, &
         !   lstyle,i_c, Station_Mnem, Ant_IDs(i_ant,i_chunk), B_a(i_ant)
      EndIf
901   Format('amove Xofft',I1,1x,F7.1,/'set color MyCol',i0,'$',/'begin graph',/'  size Xtime Ysubsize',/'  vscale 1',&
            /'  hscale 1',/'   NoBox',&
            /'   xaxis min t_min max t_max',/'   x1axis off',&
            /'   yaxis  min -scl max scl',/'   y1axis off',&
            /'     data "',A,'LOFAR_Time',A,'L.dat" d1=c1,c',I0, &
            /'     d1 line lwidth 0.14 lstyle ',i1,' color MyCol',i0,'$', &
            /' end graph' , &
            /' amove xg(t_max)-2 yg(0)'/' write ',I7 , &
            /' amove xg(t_min)+3 yg(0)'/' write ',g8.3)
!901   Format('amove Xofft',I1,1x,F5.2,/'begin graph',/'  size Xtime 4',/'  vscale 1',&
!            /'  hscale 1',/'   NoBox',&
!            /'   xaxis min t_min max t_max',/'   x1axis off',&
!            /'   yaxis  min -scl max scl',/'   y1axis off',&
!            /'     data "',A,'LOFAR_Time',A,'L.dat" d1=c1,c',I0, &
!            /'     d1 line lwidth 0.14 lstyle ',i1,' color MyCol',i0,'$',&
!            /' end graph')
903   Format('amove Xofft',I1,1x,F7.1,/'set color MyCol',i0,'$',/'begin graph',/'  size Xtime Ysubsize',/'  vscale 1',&
            /'  hscale 1',/'   NoBox',&
            /'   xaxis min 2*t_min max 2*t_max',/'   x1axis off',&
            /'   yaxis  min -scl max scl',/'   y1axis off',&
            /'     data "',A,'LOFAR_Time',A,'H.dat" d1=c1,c',I0, &
            /'     d1 line lwidth 0.14 lstyle ',i1,' color MyCol',i0,'$', &
            /' end graph' , &
            /' amove xg(2*t_max)-2 yg(0)'/' write ',I7 , &
            /' amove xg(2*t_min)+3 yg(0)'/' write ',g8.3)
      !
      !  plot SNR & error bar in side panels
      !B=(P_pulse(i_ant)/P_Bckg(i_ant))/SRNMax
      !If(B.gt.100. .or. B.lt.0.1) &
      !   write(2,*) '!Curtain, sidepanel,P', i_ant,P_pulse(i_ant), P_Bckg(i_ant), SRNMax, B
      !write(2,*) '!Curtain, sidepanel,C', j_corr,i_Peak, CCorr_Err(j_corr,i_Peak),Error_norm
      write(Unt, 911) i_eo+1, (P_pulse(i_ant)/P_Bckg(i_ant))/SRNMax, plot_offset, 'cross'
      If(CCorr_Err(j_corr,i_Peak)/Error_norm .lt. 5.) Then
         write(Unt, 911) i_eo+1, CCorr_Err(j_corr,i_Peak)/(5*Error_norm), plot_offset, 'dot'
         !write(2,*) '!Curtain, sidepanel-error', CCorr_Err(j_corr,i_Peak), Error_norm
      Else
         write(Unt, 911) i_eo+1, 0.99, plot_offset, 'circle'
      EndIf
911   Format('amove Xoffs',i0,'+Xside*',E9.3,1x,F7.1/'marker ',A,' 0.5 ')
      !
      !  plot even/odd strength ratio in central pane
      If(AntPair) Then
         !write(Unt, 910) Xoffc+Xcent*(1+(P_pulse(i_ant-1)/P_pulse(i_ant)-P_pulse(i_ant)/P_pulse(i_ant-1))/eoRatMax)/2, &
         !   plot_offset+3, 'dot'
         !write(Unt, 910) Xoffc+Xcent*(1+(B_a(i_ant-1)/B_a(i_ant)-B_a(i_ant)/B_a(i_ant-1))/eoRatMax)/2, &
         !   plot_offset+3, 'cross'
         write(Unt, 910) Xoffc+Xcent*(1+(P_pulse(i_ant-1)-P_pulse(i_ant))/(P_pulse(i_ant-1)+P_pulse(i_ant)))/2, &
            plot_offset, 'dot'
         write(Unt, 910) Xoffc+Xcent*(1+(B_a(i_ant-1)-B_a(i_ant))/(B_a(i_ant-1)+B_a(i_ant)))/2, &
            plot_offset, 'cross'
      EndIf
910   Format('amove ',F7.2,1x,F7.1,/'marker ',A,' 0.5 ')
   Enddo   !  I_ant=1,Ant_nr(i_chunk)
!    Do i_stat=1,Stat_nr
!        Call Station_ID2Mnem(Stations(I_stat),Station_Mnem)
!        plot_offset=2*i_stat  + 3
!        Write(unt,902) plot_offset,plot_offset,Station_Mnem
!902 Format('amove 4 ',F5.2,/'aline 59 ',F5.2,/'rmove 2 0 ',/'write "',A6,'"')
!    Enddo
    !
    Close(unit=unt)
    !
    Call GLEplotControl(SpecialCmnd='gle -d pdf '//trim(file)//'Ant.gle') !  Command to produce curtain plots
    Call GLEplotControl(SpecialCmnd='rm '//TRIM(file)//'Ant.gle') !  Command to delete curtain files
    !
    return
End Subroutine GLEscript_AntCurtainPlot
!====================================
Subroutine GLE_Corr()
    use DataConstants, only : Used_StationNr, DataFolder, OutFileLabel
    use DataConstants, only :  Ant_nrMax
    use Chunk_AntInfo
    use ThisSource, only : PeakNrTotal, CorrAntNrs, Nr_Corr, SafeHL, Peak_eo, ChunkNr, Peakpos, CCorr_Err, RefAntErr
    use ThisSource, only : CCorr_L, CCorr_H, CCorr_Val, SourcePos
    use FitParams, only : N_FitPar_max
    use constants, only : dp
    use StationMnemonics, only : Station_ID2Mnem
    use GLEplots, only : GLEplotControl
    Implicit none
    Integer :: J_Corr
    Integer ::  i_ant
    !
    Integer ::  i_Peak, i_eo, Station, i_c, i_chunk, lstyle, i_src, i_type
    Integer :: i, i_max, k, J_corr_st(0:Used_StationNr),Height,lr
    Logical :: nra,nrb
    character(len=6) :: Station_Mnem,Lab_a,Lab_b
    character(len=10) :: txt
    character(len=13) :: XCorrPlot
    Character(len=2) :: EveOdd(0:1)=(/'Ev','Od'/), tp(0:1)=(/'Th','Ph'/), ext
    Real(dp) :: plot_offset, lwidth
    Real(dp), parameter :: step_offset=1.0d0  ! 0.5d0
    Integer, external :: SourceNr  ! source code in  LOFLI_InputHandling.f90
    !
    Real(dp) :: CCorr(10,10,10)  ! used for plotting only
    !
    !Open(UNIT=10,STATUS='unknown',ACTION='WRITE',FILE='GLE-plots.sh')
    !If(polariz) i_tp2=1
    Do i_peak=1,PeakNrTotal
      i_eo=Peak_eo(i_peak)
      i_chunk=ChunkNr(i_peak)
      i_src=SourceNr(i_peak)
      If(i_src .le. 0 ) i_src=i_peak
      ext=EveOdd(i_eo)
      write(XCorrPlot,"(I3.3,'.gle')") i_peak
      !
      ! group all antennas of same station together, and store their numbers in 'J_corr_st'
      i=0
      J_corr_st(i)=1
      Station=Ant_Stations(CorrAntNrs(1,i_eo,i_chunk),i_chunk)
      Do j_corr=1,Nr_Corr(i_eo,i_chunk)  ! no need to distinguish LBA/HBA yet
         i_ant=CorrAntNrs(j_corr,i_eo,i_chunk)
         If(Ant_Stations(i_ant,i_chunk) .ne. station) then
             i=i+1
             J_corr_st(i)=j_corr ! Number of spectra for this station
             Station=Ant_Stations(i_ant,i_chunk)
         endif
      enddo
      i_max=i+1   ! total number of stations for which spectra that will be plotted
      J_corr_st(i_max)=Nr_Corr(i_eo,i_chunk)+1
      !
      Height=4 + (i_max+Nr_Corr(i_eo,i_chunk))*step_offset
      Open(UNIT=9,STATUS='unknown',ACTION='WRITE',FILE='XCP-'//TRIM(OutFileLabel)//trim(XCorrPlot))
      ! make sure:  export GLE_USRLIB=/Users/users/scholten/LOFLI/GLEsrc
      Write(9,"(A,i3, /A,f3.1,A, 2(/A),f3.1, /A, 3(/A,I0))") 'size 64 ',Height+8, &
         'set font pstr fontlwidth 0.08 hei 2.4*',step_offset,' just CC', &
         'include "DiscreteColor.gle"', 'set hei 1.4*',step_offset, &
         'set lwidth 0.1','t_min = ',-SafeHL(0),'t_max = ',SafeHL(0),'Height = ',Height
      !        'set lwidth 0.1','t_min = 13250 !12500 !0','t_max = 13350 !15000 !',Time_dim,'scl = ',Plot_scale
      write(9,"('amove 6 ',F5.1/,'write ',A)") Height+6.5,'"norm*100." '
      write(9,"('amove 55 ',F5.1/,'write ',A)") Height+6.5,'"Dipole #" '
      write(9,"('amove 15 1.5'/,'write ',A,I3,A)") '"Peak# ', i_peak,'"'
      write(9,"('amove 50 1.5'/,'write ',A,2(F6.1,','), F5.1,A)") '"@', SourcePos(1:3,i_peak)/1000.,' km"'
      Write(9,"(A,5(/A),I3,'=',i3,':',i5.5,'(',A,')',A,i3,A,     6(/A))") &
         'amove 4 4','begin graph','  size 55 Height','  vscale 1','  hscale 1',&
         '   title  "X-correlation traces, Peak# ',i_peak,i_chunk,Peakpos(i_peak),TRIM(ext),', src#',i_src,'"', &
         '   xtitle "time [5 ns]"', &! '   ytitle "Abs cross corr"',&
         '   xaxis min t_min max t_max ! nticks 10','   yaxis  min 0 max 2 dticks 5 dsubticks 1', &
         '   x2labels on',' end graph'
      Write(9,"(A,5(/A))" ) 'begin key','position tc','nobox','   offset 5. 0.0', &
         '   text "Excluded" lstyle 4 color black lwidth 0.1', &
         'end key'
      Write(9,"(A,4(/A),I2.2,A,/A)" ) 'begin key','position tc','nobox','   offset 14. 0.0', &
         '   text ">',SafeHL(0),' off" lstyle 2 color black lwidth 0.1', &
         'end key'
      Write(9,"(A,4(/A),I2.2,A,/A)" ) 'begin key','position tc','nobox','   offset 23. 0.0', &
         '   text ">',SafeHL(0)/4,' off" lstyle 9 color black lwidth 0.1', &
         'end key'
      !Write(9,"(/A)" ) 'begin key','position tr','nobox', & ! '   offset -0.2 -0.5', &
      !    '   text ">Safe/4" lstyle 9 color black lwidth 0.1', &
      !    'end key'
      Write(9,"('amove xg(',f4.1,') 4',/'aline xg(',f4.1,') Height+4')") RefAntErr(i_Peak), RefAntErr(i_Peak) !vertical line at zero
      !Write(9,"('amove xg(',f4.1,') 4',/'aline xg(',f4.1,') Height+4',/'set hei 0.7')") RefAntErr(i_Peak), RefAntErr(i_Peak) !vertical line at zero
      !
      plot_offset=4.+1.*step_offset
      nra=.true.
      Do i=1,i_max  ! Loop over stations
         i_type=Mod(Ant_Stations(CorrAntNrs(J_corr_st(i-1),i_eo,i_chunk),i_chunk),2)
         write(txt,"(i4.4,A,I3.3,A2)") Ant_Stations(CorrAntNrs(J_corr_st(i-1),i_eo,i_chunk),i_chunk),'_',i_peak,ext
         !write(2,*) i,txt,J_corr_st(i)-1-J_corr_st(i-1)
         Open(Unit=11,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//'LOFAR_Corr-'//txt//'.dat')
         write(11,*) '!',Ant_Stations(J_corr_st(i-1),i_chunk), &
             (Ant_IDs(CorrAntNrs(k,i_eo,i_chunk),i_chunk),k=J_corr_st(i-1), J_corr_st(i)-1)
         If(i_type.eq.0) Then
            Do k=-SafeHL(i_type),+SafeHL(i_type)
                write(11,*) k,(CCorr_L(k,J_corr_st(i-1):J_corr_st(i)-1,i_Peak))
            enddo
         Else
            Do k=-SafeHL(i_type),+SafeHL(i_type)
                write(11,*) k,(CCorr_H(k,J_corr_st(i-1):J_corr_st(i)-1,i_Peak))
            enddo
         EndIf
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
902 Format('set color black'/'amove 50 ',F6.2,/'aline 59 ',F6.2,/'rmove 0.5 0 ',/'write "',A6,'"')
         Write(9,903) plot_offset,Lab_b,plot_offset
903 Format('amove 0.5 ',F6.2,/'write "',A6,'"',/'amove 4 ',F6.2,/'rline 5 0')
         lwidth=0.1
         Do J_corr=J_corr_st(i-1),J_corr_st(i)-1
            !plot_offset=4.+(i+J_corr-1)*step_offset
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
            ! i_type+1 = scale time
            Write(9,901) plot_offset,i_type+1,i_type+1,trim(DataFolder), TRIM(txt),k,k+1,k,lwidth,lstyle, i_c, &
               i_c,   NINT(CCorr_Val(j_corr,i_Peak)*100),  Ant_IDs(CorrAntNrs(j_corr,i_eo,i_chunk),i_chunk)
901     Format('amove 4 ',F6.2,/'begin graph',/'  size 55 2',/'  vscale 1',&
            /'  hscale 1',/'   NoBox',&
            /'   xaxis min ',I1,'*t_min max ',I1,'*t_max',/'   x1axis off',&
            /'   yaxis  min 0 max 1',/'   y1axis off',&
            /'     data "',A,'LOFAR_Corr-',A,'.dat" d',I0,'=c1,c',I0,&
            /'     d',I0,' line lwidth ',f4.2,' lstyle ',i1,' color MyCol',i0,'$', &
            /' end graph' ,&
            /'set color MyCol',i0,'$',/ 'amove 7 yg(0)'/'write ',I3, &
            / 'amove 55 yg(0)'/'write ', I2)
            ! lwidth=0.05  ! Reset to normal width
            plot_offset=plot_offset+1*step_offset
         Enddo
         plot_offset=plot_offset+1*step_offset  ! creates an extra offset separating stations
      enddo   ! i=1,i_max
      !
      Close(unit=9) !
      Call GLEplotControl(SpecialCmnd='gle -d pdf XCP-'//TRIM(OutFileLabel)//trim(XCorrPlot)) !  Command to produce curtain plots
      !Call GLEplotControl(SpecialCmnd='gle -d jpg -r 300 XCP-'//TRIM(OutFileLabel)//trim(XCorrPlot)) ! .jgp is about 20* larger than .pdf and still poor quality
      Call GLEplotControl(SpecialCmnd='rm XCP-'//TRIM(OutFileLabel)//trim(XCorrPlot)) !  Command to delete curtain files
   enddo  !  i_peak=1,PeakNrTotal
   Return
End Subroutine GLE_Corr
!=================================
! ===========================================================================
Subroutine GLEscript_EFieldCurtain(unt, file, WWidth, i_chunk, FileA, Label, dChi_ap, dChi_at, Power_p, Power_t, Chi2pDF, VoxLoc)
!  To make curtain plots for the polarized fields for the antennas
   use constants, only : dp,pi
   !use DataConstants, only : DataFolder!, OutFileLabel
   !use Chunk_AntInfo, only : Ant_Stations, Ant_nr, Ant_IDs, Ant_pos, Ant_RawSourceDist, Nr_UniqueStat
   use Chunk_AntInfo, only : Ant_Stations, Nr_UniqueStat, Unique_StatID,  Ant_IDs, Ant_pos
   !use DataConstants, only : Station_nrMax
   Use Interferom_Pars, only : IntFer_ant,  Nr_IntferCh ! the latter gives # per chunk
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
   !Character(len=6) :: Station_Mnem
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
   !Write(2,*) 'GLEscript_EFieldCurtain:Nr_UniqueStat=',Nr_UniqueStat
   Open(UNIT=UNT,STATUS='unknown',ACTION='WRITE',FILE=trim(file)//'.gle')
   Write(unt,"(A)") '! COMMAND:  gle -d pdf '//trim(file)//'.gle'
    ! make sure:  export GLE_USRLIB=/Users/users/scholten/LOFLI/GLEsrc
   Write(unt,"(A,I2,3(/A),2(/A,I0))") 'size 63 ',Nr_UniqueStat+9,'set font pstr fontlwidth 0.08 hei 1.2 just CC',&
      'include "DiscreteColor.gle"',&
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
         Write(Unt,903) trim(FileA)//'PhiMod_'//TRIM(Label), dn, j_IntFer+1, dn, '4', MODULO(i_ant,11)
         counter=counter+1
         dChi_sp(i_stat)=dChi_sp(i_stat)*(counter-1)/counter + dChi_ap(j_IntFer)/counter  ! keep running mean
         StPowr_p(i_stat)=StPowr_p(i_stat)*(counter-1)/counter + Power_p(j_IntFer)/counter  ! keep running mean
      EndDo
      !write(2,*) 'dn:',dn,j_IntFer
      Write(unt,"(T6,'let d',i0,'=0.',/T6,'d',i0,' line lstyle 0 lwidth 0.01')") dn+1, dn+1
      Write(Unt,"(' end graph')")
      !Write(unt,"('amove ',F5.1,' yg(0)',/'rline ',F5.1' 0')") HOffSt_p, PlotW
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
         Write(Unt,903) trim(FileA)//'ThMod_'//TRIM(Label), dn, j_IntFer+1, dn, '4', MODULO(i_ant,11)
         counter=counter+1
         dChi_st(i_stat)=dChi_st(i_stat)*(counter-1)/counter + dChi_at(j_IntFer)/counter  ! keep running mean
         StPowr_t(i_stat)=StPowr_t(i_stat)*(counter-1)/counter + Power_t(j_IntFer)/counter  ! keep running mean
         Phi(i_stat)=180.+atan2( VoxLoc(2)-Ant_pos(2,i_ant,i_chunk) , VoxLoc(1)-Ant_pos(1,i_ant,i_chunk) ) *180./pi ! \phi=0 = south
      EndDo
      Write(unt,"(T6,'let d',i0,'=0.',/T6,'d',i0,' line lstyle 0 lwidth 0.01')") dn+1, dn+1
      Write(Unt,"(' end graph')")
      !Write(unt,"('amove ',F5.1,' yg(0)',/'rline ',F5.1' 0')") HOffSt_t, PlotW
   EndDo ! i_stat=1,Nr_UniqueStat
   !
901   Format('amove ',2F5.1,/'begin graph',/'  size 27 2',/'  vscale 1',&
         /'  hscale 1',/'   NoBox',&
         /'   xaxis min t_min max t_max',/'   x1axis off',&
         /'   yaxis ',/'   y1axis off')
!         /'   yaxis  min -scl max scl',/'   y1axis off')
902   Format('     data "',A,'.dat" d',i0,'=c1,c',I0,&
         /'     d',i0,' line lwidth 0.04 lstyle ',A,' color MyCol',i0,'$')   ! data
903   Format('     data "',A,'.dat" d',i0,'=c1,c',I0,&
         /'     d',i0,' line lwidth 0.1 lstyle ',A,' color MyCol',i0,'$')   ! model; thicker dots
    !
    Write(unt,"(A,2F5.1,/A,F6.2,A)") 'amove ',HOffSt_t,Nr_UniqueStat+7.5, 'write "Chi^2/dof=',Chi2pDF,'"'
    Write(unt,"('set lstyle 0',/'set lwidth 0.01',/'set hei 0.7')")
    !write(2,*) 'Chi2pDF',Chi2pDF
    !write(2,*) 'dChi_sp(i_stat):',dChi_sp(:)
    Do i_stat=1,Nr_UniqueStat  ! put labels
!   HOffSt_t=HOffSt_p+PlotW+Sep
        plot_offset=i_stat  + 3
        !Write(unt,"('amove ',2F5.1,/'rline ',F5.1' 0')") HOffSt_p, plot_offset,    PlotW  ! hline
        Write(unt,"('amove ',2F5.1,/A,A5,A)") HOffSt_p-1.5, plot_offset,'write "',Statn_ID2Mnem(Unique_StatID(i_stat)),'"'  ! text on left
        whiteness=1/(1+2*dChi_sp(i_stat)/Chi2pDF)
        Write(unt,904) HOffSt_p-4.0+PlotW, plot_offset+.5, 1.-whiteness, whiteness, dChi_sp(i_stat)            ! text in right side
        Write(unt,905) HOffSt_p-1.5+PlotW, plot_offset+.5,  sqrt(StPowr_p(i_stat))            ! text in right side
        !Write(unt,"('amove ',2F5.1,/'rline ',F5.1' 0')") HOffSt_t, plot_offset,    PlotW
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
End Subroutine GLEscript_EFieldCurtain
!====================================

!====================================
