!=================================
!==========================================
Subroutine EI_PrntFitPars(X)
    use FitParams
    !use FittingParameters
    use ThisSource, only : PeakPos, PeakNrTotal, ChunkNr, PeakChiSQ
    use Chunk_AntInfo, only : Unique_StatID, Start_time
    use StationMnemonics, only : Station_ID2Mnem
    use StationMnemonics, only : Statn_ID2Mnem
    use DataConstants, only : ChunkNr_dim
    use Constants, only : Sample
    Implicit none
    real ( kind = 8 ), intent(in) :: X(N_FitPar_max)
    integer ( kind = 4 ) :: i,j, i_Peak, i_chunk, i_fit
    Character(len=5) :: Station_Mnem
    Character(len=1) :: FitParam_Mnem(4)=(/'N','E','h','t'/)
    !
    If(N_FitStatTim .gt. 0) then
        Do i=1,N_FitStatTim
            Call Station_ID2Mnem(Unique_StatID(FitParam(i)),Station_Mnem)
            write(2,"(A,i2,'), ',A5,50F11.3)") 'Fit_TimeOffsets[samples](',i,Station_Mnem,(X(j),j= X_Offset(i-1),X_Offset(i)-1)
        enddo
    endif
    i_chunk=0
    If((N_FitPar-N_FitStatTim) .gt. 0) then
        Do i_Peak = 1,PeakNrTotal
            i_fit=N_FitStatTim+i_Peak
            If(i_chunk.ne.ChunkNr(i_Peak)) then
               i_chunk=ChunkNr(i_Peak)
               write(2,"(A,i2,A,I11,A,f10.3,A)") 'block=',i_chunk,', starting time=',Start_time(i_chunk),&
                  '[samples]=',Start_time(i_chunk)*Sample*1000.,'[ms]'
            Endif
            Write(2,"(i3,2x,i2,A,I6,A,3F8.3)", ADVANCE='NO') i_Peak,ChunkNr(i_Peak),&
                ', PeakPos',PeakPos(i_Peak),', source@', X( X_Offset(i_fit-1):X_Offset(i_fit)-1 )/1000.
            write(2,"(', chi^2/df=',F7.2)", ADVANCE='NO') PeakChiSQ(i_Peak)
            !If(Station_nrMax.gt.0 .and. (SUM(Dropped(:,i_Peak)).gt.0)) then
            !      write(2,"(' Excluded:')", ADVANCE='NO')
            !   Do i=1,Station_nrMax
            !      If(ExclStatNr(i,i_peak).eq.0) exit
            !      Write(2,"(1x,A5)", ADVANCE='NO') Statn_ID2Mnem(ExclStatNr(i,i_peak))
            !   Enddo
            !Endif
            write(2,*) ' '
        enddo  !  i_Peak = 1,PeakNrTotal
    Endif
End Subroutine EI_PrntFitPars
!==========================================
Subroutine EIPrntNewSources()
   use ThisSource, only : PeakNrTotal, ChunkNr, PeakChiSQ, PeakRMS, PeakPos, SourcePos
   use Chunk_AntInfo, only : Start_time
    !use StationMnemonics, only : Statn_ID2Mnem
   use Constants, only : Sample
   Implicit none
   !real ( kind = 8 ), intent(in) :: X(N_FitPar_max)
   integer ( kind = 4 ) :: i,j, i_Peak, i_chunk
   integer, external :: XIndx
   Character(len=5) :: Station_Mnem
   Character(len=1) :: FitParam_Mnem(4)=(/'N','E','h','t'/)
   !
   i_chunk=0
   Do i_Peak = 1,PeakNrTotal
   If(i_chunk.ne.ChunkNr(i_Peak)) then
      i_chunk=ChunkNr(i_Peak)
      write(2,"(A,i2,A,I11,A,f10.3,A)") 'block=',i_chunk,', starting time=',Start_time(i_chunk),&
         '[samples]=',Start_time(i_chunk)*Sample*1000.,'[ms]'
   Endif
   Enddo
   !
   Write(2,*) 'Nr,eo,Blk,PPos,(Northing,   Easting,   height, Del-t); RMS[ns], sqrt(chi^2/df), Excluded: '
   i_chunk=1
   j=0
   Do i_Peak = 1,PeakNrTotal
      Write(2,"(i3,2x,i2,I8)", ADVANCE='NO') i_Peak,ChunkNr(i_Peak),PeakPos(i_Peak)
      Write(2,"(3(F10.2,','))", ADVANCE='NO') SourcePos(:,i_Peak)
      write(2,"(';',F7.2)") PeakChiSQ(i_Peak)
   Enddo  !  i_Peak = 1,PeakNrTotal
End Subroutine EIPrntNewSources
!==========================================
Subroutine EIPrntCompactSource(i_Peak)
    use FitParams, only : station_nrMax
    !use FittingParameters
    use ThisSource, only : PeakPos, RefAntErr, PeakNrTotal, ChunkNr, PeakChiSQ, ExclStatNr, Dropped, SourcePos
    use Chunk_AntInfo, only : Unique_StatID
    use StationMnemonics, only : Statn_ID2Mnem
    Implicit none
    Integer, intent(in) :: i_Peak
    integer ( kind = 4 ) :: i
    Character(len=5) :: Station_Mnem
    !Character(len=1) :: FitParam_Mnem(4)=(/'N','E','h','t'/)
    !
   If(Station_nrMax.gt.0 .and. (SUM(Dropped(:,i_Peak)).gt.0)) then
      Do i=1,Station_nrMax
         If(ExclStatNr(i,i_peak).eq.0) exit
         Write(2,"(1x,A5)", ADVANCE='NO') Statn_ID2Mnem(ExclStatNr(i,i_peak))
      Enddo
      write(2,"(/,'exclude  ')", ADVANCE='NO')
      Do i=1,Station_nrMax
         If(Dropped(i,i_peak).eq.0) cycle
         Write(2,"(1x,A5)", ADVANCE='NO') Statn_ID2Mnem(Unique_StatID(i))
      Enddo
   Endif
   write(2,*) ' '
   Return
End Subroutine EIPrntCompactSource
!==========================================
Subroutine EIX2Source(X)
! Move fit results fron X to Sourcepos and FineOffset_
    use DataConstants, only : Station_nrMax
    use DataConstants, only : Polariz
    use ThisSource, only : SourcePos, RefAntErr, PeakNrTotal
    use Chunk_AntInfo, only : Unique_StatID, Nr_UniqueStat, Tot_UniqueAnt,  Unique_SAI
    use FitParams
    use constants, only : dp
    Implicit none
    real ( kind = 8 ), intent(in) :: X(N_FitPar_max)
    integer :: i_Peak
    integer, external :: XIndx
    integer :: i, k,j
    Real(dp) :: dt
    !
    !Call EI_PrntFitPars(X)
    !write(2,*) 'SourcePos:',SourcePos(:,1)
    !write(2,*) 'X_Offset:',N_FitStatTim+1 , N_FitPar,X_Offset(0:N_FitPar)
    !write(2,*) 'X:',X(1:10)
    Flush(unit=2)
    If(N_FitPar .gt. N_FitStatTim) then  ! sources are fitted also
        Do i=N_FitStatTim+1 , N_FitPar  != number of fitted sources=i_peak
         i_Peak=i-N_FitStatTim
         SourcePos(1:3,i_Peak) = X( X_Offset(i-1):X_Offset(i)-1 )
    !write(2,*) 'SourcePos(1,i_peak)', i_peak,SourcePos(:,i_peak)
        enddo
    endif
    !i_Peak=2
    !Write(*,*) 'SourcePos(2)=',SourcePos(1,i_Peak),XIndx(1,i_Peak),Dual
    !
    If(N_FitStatTim .gt. 0) then
        Do i=1, N_FitStatTim
            k = FitParam(i)   ! =Station number
            !write(2,*) 'X2Sourcei',i,k, X_Offset(i), Tot_UniqueAnt(k), Tot_UniqueAnt(k-1)
            !write(2,*) 'X2SourceX',X(X_Offset(i):X_Offset(i)+Tot_UniqueAnt(k)-Tot_UniqueAnt(k-1))
            !write(2,*) 'X2Sourceu',Unique_SAI(Tot_UniqueAnt(k-1)+1:Tot_UniqueAnt(k))
            If(Fit_AntOffset) then
                 Do j= 1,(Tot_UniqueAnt(k)-Tot_UniqueAnt(k-1))  ! there is same fit_timeoffset for even and odd
                    Fit_TimeOffsetAnt(Tot_UniqueAnt(k-1)+j)=X(X_Offset(i-1)+j-1)
                    !write(2,*) 'X2Source', i,k,j,X(X_Offset(i-1)+j), X_Offset(i-1), Tot_UniqueAnt(k-1)+j
                 Enddo
            Else
              Fit_TimeOffsetStat(k) = X(X_Offset(i-1))
            endif
            !write(2,*) 'X2SourceF', i,k,Fit_TimeOffsetAnt(Tot_UniqueAnt(k-1)+1:Tot_UniqueAnt(k))
            !write(2,*) 'X2SourceX', X_Offset(i-1),X_Offset(i)-1,X(X_Offset(i-1):X_Offset(i)-1)
        Enddo  ! i=1, N_FitStatTim
    endif
    !
    !Call EI_PrntFitPars(X)
    !write(2,*) 'x2s;Fit_TimeOffsetAnt',Fit_TimeOffsetAnt(150:180)
    !write(2,*) 'Fit_TimeOffsetStat=',Fit_TimeOffsetStat(1:10)
    return
End Subroutine EIX2Source
!-----------------------------------------------
Subroutine EI_PolarizPeak(i_Peak)
   use constants, only : dp
   use DataConstants, only : Ant_nrMax
   use Interferom_Pars, only :  Nr_IntFerMx, Nr_IntferCh ! the latter gives # per chunk
   use Interferom_Pars, only : W_ap, W_at, Cnu0, Cnu1, IntfNuDim, CTime_p, CTime_t, Noise_p, Noise_t, AntPeak_OffSt
   Use Interferom_Pars, only :  N_Smth, N_fit, Chi2pDF
   use ThisSource, only : ChunkNr, PeakPos, sourcepos, PeakChiSQ
   Implicit none
   Integer, intent(in) :: i_Peak
   Integer :: IntfBase, i_chunk, i_1, i_2, Outpt, j_IntFer, Windw
   Real(dp) :: ChiSq, DelChi(-N_fit:+N_fit,1:2*Nr_IntFerMx)
   Real(dp) :: PartChiSq(1:Nr_IntFerMx) !
   Integer :: PartChiSqInt(1:Nr_IntFerMx)
   Real(dp) :: FitDelay(1:Ant_nrMax), VoxLoc(1:3), del_1, del_2
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
   Call EI_PolSetUp(Nr_IntFerCh(i_chunk), IntfBase, i_chunk, SourcePos(1,i_Peak), &
      AntPeak_OffSt(1,i_Peak), Cnu0(0,1,i_Peak), Cnu1(0,1,i_Peak), W_ap(1,i_Peak), W_at(1,i_Peak))
   Windw=3*N_Smth
   Call EI_Weights(Nr_IntFerCh(i_chunk), IntfNuDim, i_chunk, SourcePos(1,i_Peak), &
      AntPeak_OffSt(1,i_Peak), W_ap(1,i_Peak), W_at(1,i_Peak))
   !write(2,*) 'EI_PolarizPeak:W_ap='
   write(Label,"('Pk ',i4.2)") i_Peak
   Call EI_PolGridDel(Nr_IntFerCh(i_chunk), FitDelay, IntfNuDim, i_chunk, SourcePos(1,i_Peak), &
      AntPeak_OffSt(1,i_Peak), W_ap(1,i_peak), W_at(1,i_peak), Cnu0(0,1,i_peak), Cnu1(0,1,i_peak), Outpt, DelChi, Label)
   PeakChiSQ(i_Peak)=Chi2pDF
   Call WriteDelChiPeak(i_chunk, DelChi,PartChiSq,PartChiSqInt)
   !j_IntFer=PartChiSqInt(Nr_IntFerCh(i_chunk))
   !write(Label,"(i2.2,i3.3)") i_Peak,j_IntFer
   !Call TimeTracePlot(j_IntFer, IntfBase, i_chunk, SourcePos(1,i_Peak), Windw, Label)
   !
   Return
   Del_1=2  ! in meter
   Del_2=2
   Outpt=1
   write(2,*) i_Peak, 'Del_1=',Del_1,', Del_2=',Del_2, ' [m]'
   Do i_1=-2,2
      Do i_2=-2,2
         VoxLoc(1)=SourcePos(1,i_Peak) + i_1*del_1
         VoxLoc(2)=SourcePos(2,i_Peak) + i_2*del_2
         VoxLoc(3)=SourcePos(3,i_Peak)
         write(Label,"('Gr ',I2,',',i2)") i_1,i_2
         Call EI_Weights(Nr_IntFerCh(i_chunk), IntfNuDim, i_chunk, VoxLoc(:),  &
            AntPeak_OffSt(1,i_Peak), W_ap(1,i_Peak), W_at(1,i_Peak))
         Call EI_PolGridDel(Nr_IntFerCh(i_chunk), FitDelay, IntfNuDim, i_chunk, VoxLoc(:), AntPeak_OffSt(1,i_Peak), &
            W_ap(1,i_peak), W_at(1,i_peak), Cnu0(0,1,i_peak), Cnu1(0,1,i_peak), Outpt, DelChi, Label)
      Enddo
   Enddo
   !
   Return
End Subroutine EI_PolarizPeak


!===================================
