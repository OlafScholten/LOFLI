    !include 'C:/OlafsUtil/Spline/MGMR3D_spline.f90'
    !Include '../../NumLib/Spline/MGMR3D_spline.f90'  ! ~/NumLib
    Include 'MGMR3D_spline.f90'  ! ../../NumLib
    !Include '/home/olaf/NumLib/Spline/MGMR3D_spline.f90'  ! should have been  ~/NumLib  which did not work
    ! Fitter included after the module, needed to properly define the function to optimize
    ! v17 : offset in Nr_TimeOffset reduced by 1 in XIndx
!=================================
!-----------------------------------------------
  Subroutine SetFitParamsLMA(X,first,FitPos)
    !use unque
    use FitParams, only : FitParam, X_Offset, Fit_TimeOffsetAnt, ParamScaleFac, Fit_TimeOffsetStat
    use FitParams, only : N_FitPar_max, N_FitPar, N_FitStatTim, Nr_TimeOffset, Fit_AntOffset
    use DataConstants, only : Station_nrMax, Ant_nrMax, RunMode, ChunkNr_dim
    use DataConstants, only : Polariz
    use ThisSource, only : Nr_Corr, CorrAntNrs, SourcePos, RefAntErr, PeakNrTotal
    use Chunk_AntInfo, only : Ant_Stations, Ant_nr, Unique_StatID, Nr_UniqueStat, Tot_UniqueAnt, Unique_SAI
    use StationMnemonics, only : Station_Mnem2ID
    implicit none
    Real(kind = 8), intent(out) :: X(N_FitPar_max)
    logical, intent(inout) :: First
    integer, intent(in) :: FitPos(4)
    integer :: i,j,k,j_corr, nxx, i_ant, i_Peak, i_eo, i_chunk
    Integer :: vec(Ant_nrMax)
    !Equivalence (Ant_Stations(1,1),vec(1))
    logical,dimension(Station_nrMax) :: mask
    integer, save ::  FP_s(0:N_FitPar_max), FP(0:4)
    Character(len=5) :: FP_MNem(0:N_FitPar_max)
    character*180 :: lname
    character*4 :: option
    Logical,save :: FitNoSources
    !integer,dimension(:),allocatable :: Station_IDs
    integer, external :: XIndx
    !
    If(first) then
        !Do i=1, Station_nrMax
        !    Read(*,*) FineOffset_StatID(i), FineOffset_Station(i)  ! units of samples
        !    ! After fitting the values of 'FineOffset_Station(i)' are replaced by the fitted values
        !    if(FineOffset_StatID(i) .le. 0) exit
        !Enddo
        !vec(:)=Ant_Stations(:,1)
        !
        FP_s=0
        !Call unique(vec,Station_IDs)
        !Unique_Nr=size(Station_IDs)
        !Unique_StatID(:)=Station_IDs(:)
        FP_MNem(:)=' '
        If(RunMode.eq.1) goto 1
        Call GetStationFitOption(FP_s, FitNoSources)
    endif
  1   Continue
    !
    N_FitStatTim=0
    Nr_TimeOffset=1
    if(FP_s(1).GT.0) then
     Do i=1,N_FitPar_max
        ! write(2,*) 'SetFitParamsLMA:',i,FP_s(i)
        if(FP_s(i) .le. 0) exit
        Do i_chunk=1, ChunkNr_dim
         Do i_eo=0,1
            If(any(FP_s(i)==Ant_Stations(CorrAntNrs(1:Nr_Corr(i_eo,i_chunk),i_eo,i_chunk),i_chunk)) ) goto 2  ! check if this station included in present fit
         Enddo
        Enddo
        exit
  2     Continue
        Do k=1, Nr_UniqueStat
          If(Unique_StatID(k).eq. FP_s(i)) then
            N_FitStatTim= N_FitStatTim + 1 ! position in 'FitParam'
            FitParam(N_FitStatTim)= k ! the timing of station "Unique_StatID(k)" will be fitted.
            X_Offset(N_FitStatTim)=Nr_TimeOffset   ! position in 'X' for this station, or the first antanna for this station
            If(Fit_AntOffset) then
               If(Polariz) Then
                 ! write(2,*) 'SetFitParam',k, Tot_UniqueAnt(k)-Tot_UniqueAnt(k-1)
                 Do j= 1,(Tot_UniqueAnt(k)-Tot_UniqueAnt(k-1))/2
                    X(Nr_TimeOffset)=(Fit_TimeOffsetAnt(Tot_UniqueAnt(k-1)+2*j-1)+ Fit_TimeOffsetAnt(Tot_UniqueAnt(k-1)+2*j))/2.
                    !write(2,*) 'timeoffseteven-odd', Tot_UniqueAnt(k-1)+2*j-1, Nr_TimeOffset,  &
                    !    Unique_SAI(Tot_UniqueAnt(k-1)+2*j-1), Unique_SAI(Tot_UniqueAnt(k-1)+2*j), X(Nr_TimeOffset)
                    Fit_TimeOffsetAnt(Tot_UniqueAnt(k-1)+2*j-1) = X(Nr_TimeOffset)
                    Fit_TimeOffsetAnt(Tot_UniqueAnt(k-1)+2*j) = X(Nr_TimeOffset) ! set equal to average
                    Nr_TimeOffset=Nr_TimeOffset+1
                 Enddo
               Else
                 Do j= 1,Tot_UniqueAnt(k)-Tot_UniqueAnt(k-1)
                    X(Nr_TimeOffset)=Fit_TimeOffsetAnt(Tot_UniqueAnt(k-1)+j)
                    Nr_TimeOffset=Nr_TimeOffset+1
                 Enddo
               EndIf
            Else
              X(Nr_TimeOffset)=Fit_TimeOffsetStat(k)
              Nr_TimeOffset=Nr_TimeOffset+1
            endif
            exit
          endif
        Enddo ! k=1, Nr_UniqueStat
     Enddo
    endif
    X_Offset(N_FitStatTim+1)=Nr_TimeOffset
    !
    N_FitPar=N_FitStatTim
    !write(*,*) 'FP',FP(:),N_FitPar
    !write(2,"(A,40('(',i2,i3,')'))") '(FitParam,X_Offset)=',(FitParam(i),X_Offset(i),i=1,N_FitStatTim+1)
    if(.not. FitNoSources .and. FitPos(1).gt.0) then
        Do i=1,4        ! aug19: changed from 3 to 4 where 4 is the error in the reference antenna due to noise for this peak
            if((FitPos(i) .le. 0) .or. (FitPos(i) .gt. 4) ) exit
            N_FitPar = N_FitPar +1
            FitParam(N_FitPar)=FitPos(i)
            !
            Do i_Peak=1,PeakNrTotal
                k=XIndx(i,i_Peak)
                  !write(2,*) 'from XIndx:', k,i,i_peak
                If(k .gt. N_FitPar_max) Then
                  write(2,*) 'Too many fit parameters, ',k, &
                     'exceeds max of', N_FitPar_max,' for i,i_peak=',i,i_Peak
                  stop 'SetFitParamsLMA; Too many fit parameters'
                Endif
                If(FitPos(i).eq.4) then
                    X( k ) = RefAntErr(i_Peak)
                    !write(2,*) 'SetFitParamsLMA:', i,i_peak,XIndx(i,i_Peak),RefAntErr(i_Peak)
                else
                    X( k ) = SourcePos(FitPos(i),i_Peak)
                EndIf
            enddo
            !
        enddo
    endif
    !
    First=.false.
    ParamScaleFac(:) = 1.
    !write(*,*) 'N_FitPar,N_FitStatTim', N_FitPar,N_FitStatTim
    !write(2,*) 'Fit_TimeOffsetAnt',Fit_TimeOffsetAnt(150:180)
    !write(2,"(A,10F11.2)") 'RefAntErr_SFP(i_Peak):',(RefAntErr(i_Peak),i_Peak=1,PeakNrTotal)
    !write(2,"(A,10F11.2)") 'X_SFP=',(X(i),i=1,10)
    !write(2,"(A,10F11.2)") 'X_SFP=',(X(i),i=11,19)
    !Call PrntFitPars(X)
    !
    return
  end subroutine SetFitParamsLMA
!==========================================
Subroutine PrntFitPars(X)
    use FitParams
    !use FittingParameters
    use ThisSource, only : PeakPos, PeakNrTotal, Peak_eo, ChunkNr, PeakChiSQ, PeakRMS, ExclStatNr, Dropped
    use Chunk_AntInfo, only : Unique_StatID, Start_time
    use StationMnemonics, only : Station_ID2Mnem
    use StationMnemonics, only : Statn_ID2Mnem
    !use DataConstants, only : ChunkNr_dim
    use Constants, only : Sample
    Implicit none
    real ( kind = 8 ), intent(in) :: X(N_FitPar_max)
    integer ( kind = 4 ) :: i,j, i_Peak, i_chunk
    integer, external :: XIndx
    Character(len=5) :: Station_Mnem
    Character(len=1) :: FitParam_Mnem(4)=(/'N','E','h','t'/)
    !
    If(N_FitStatTim .gt. 0) then
        Do i=1,N_FitStatTim
            Call Station_ID2Mnem(Unique_StatID(FitParam(i)),Station_Mnem)
            write(2,"(A,i2,'), ',A5,50F11.3)") 'Fit_TimeOffsets[samples](',i,Station_Mnem,(X(j),j= X_Offset(i),X_Offset(i+1)-1)
        enddo
    endif
    i_chunk=0
    If((N_FitPar-N_FitStatTim) .gt. 0) then
        Do i_Peak = 1,PeakNrTotal
            If(i_chunk.ne.ChunkNr(i_Peak)) then
               i_chunk=ChunkNr(i_Peak)
               write(2,"(A,i2,A,I11,A,f10.3,A)") 'block=',i_chunk,', starting time=',Start_time(i_chunk),&
                  '[samples]=',Start_time(i_chunk)*Sample*1000.,'[ms]'
            Endif
            Write(2,"(i3,2i2,A,I6,A)", ADVANCE='NO') i_Peak,Peak_eo(i_Peak),ChunkNr(i_Peak),&
                ', PeakPos',PeakPos(i_Peak),', source@'
            Do i=1,N_FitPar-N_FitStatTim
                Write(2,"(A1,'=',F9.2,',')", ADVANCE='NO') &
                    FitParam_Mnem(FitParam(N_FitStatTim+i)),X( XIndx(i,i_Peak) )
            Enddo
            write(2,"(' RMS[ns]=',F7.2,', sqrt(chi^2/df)=',F7.2)", ADVANCE='NO') PeakRMS(i_Peak),sqrt(PeakChiSQ(i_Peak))
            If(Station_nrMax.gt.0 .and. (SUM(Dropped(:,i_Peak)).gt.0)) then
               !If(WriteCalib) then
               !   write(2,"(/,'exclude  ')", ADVANCE='NO')
               !Else
                  write(2,"(' Excluded:')", ADVANCE='NO')
               !Endif
               Do i=1,Station_nrMax
                  If(ExclStatNr(i,i_peak).eq.0) exit
                  Write(2,"(1x,A5)", ADVANCE='NO') Statn_ID2Mnem(ExclStatNr(i,i_peak))
               Enddo
               !write(2,"(/,'exclude  ')", ADVANCE='NO')
               !Do i=1,Station_nrMax
               !   If(Dropped(i,i_peak).eq.0) cycle
               !   Write(2,"(1x,A5)", ADVANCE='NO') Statn_ID2Mnem(Unique_StatID(i))
               !Enddo
            Endif
            write(2,*) ' '
        enddo  !  i_Peak = 1,PeakNrTotal
    Endif
End Subroutine PrntFitPars
!==========================================
Subroutine PrntNewSources()
    !use FitParams, only : station_nrMax
    !use FittingParameters
    !use ThisSource, only : PeakPos, RefAntErr, PeakNrTotal, ChunkNr, PeakChiSQ, PeakRMS, ExclStatNr, Dropped, SourcePos
   use ThisSource, only : PeakNrTotal, ChunkNr, Peak_eo
   use Chunk_AntInfo, only : Start_time
    !use StationMnemonics, only : Statn_ID2Mnem
   use Constants, only : Sample
   use DataConstants, only : Polariz
   Implicit none
   !real ( kind = 8 ), intent(in) :: X(N_FitPar_max)
   integer ( kind = 4 ) :: i,j, i_Peak, i_chunk, i_eo
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
   If(Polariz) Then
      i_chunk=1
      j=0
      Do i_Peak = 1,PeakNrTotal
         If(i_chunk.eq.ChunkNr(i_Peak)) Then
            j=j+1
         Else
            i_eo=1
            Do i=i_Peak-j, i_Peak-1
               Call PrntCompactSource(i,i_eo)
            Enddo
            j=1
         EndIf
         i_chunk=ChunkNr(i_Peak)
         i_eo=0
         Call PrntCompactSource(i_Peak,i_eo)
      Enddo  !  i_Peak = 1,PeakNrTotal
      i_eo=1
      Do i=PeakNrTotal-j+1, PeakNrTotal
         Call PrntCompactSource(i,i_eo)
      Enddo
   Else
      Do i_Peak = 1,PeakNrTotal
         Call PrntCompactSource(i_Peak,Peak_eo(i_Peak))
      Enddo  !  i_Peak = 1,PeakNrTotal
   EndIf
End Subroutine PrntNewSources
!==========================================
Subroutine PrntCompactSource(i_Peak,i_eo)
   use constants, only : dp
   use FitParams, only : station_nrMax
   !use FittingParameters
   use ThisSource, only : PeakPos, RefAntErr, PeakNrTotal, ChunkNr, PeakChiSQ, PeakRMS, ExclStatNr, Dropped, SourcePos
   use Chunk_AntInfo, only : Unique_StatID, Ant_pos, Ant_RawSourceDist, RefAnt
   use StationMnemonics, only : Statn_ID2Mnem
   Implicit none
   Integer, intent(in) :: i_Peak,i_eo
   integer ( kind = 4 ) :: i,k, i_chunk
   Character(len=5) :: Station_Mnem
   Real(dp) :: RDist
   !Character(len=1) :: FitParam_Mnem(4)=(/'N','E','h','t'/)
   !
   ! translate peakposition in the reference antenna to one for a (virtual) antenna at the core
   i_chunk=ChunkNr(i_Peak)
   Call RelDist(SourcePos(:,i_Peak),Ant_pos(:,RefAnt(i_chunk,i_eo),i_chunk),RDist) !distance to ant - distance to core in samples
   k=PeakPos(i_Peak)-NINT(RDist-Ant_RawSourceDist(RefAnt(i_chunk,i_eo),i_chunk)) ! in samples due to signal travel distance to reference
   !
   Write(2,"('C',3i2,I8)", ADVANCE='NO') i_Peak,i_eo,ChunkNr(i_Peak),k
   Write(2,"(3(F10.2,','))", ADVANCE='NO') SourcePos(:,i_Peak)
   Write(2,"(F8.2)", ADVANCE='NO') RefAntErr(i_Peak)
   write(2,"(';',F7.2,',',F7.2)", ADVANCE='NO') PeakRMS(i_Peak),sqrt(PeakChiSQ(i_Peak))
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
End Subroutine PrntCompactSource
!==========================================
Subroutine X2Source(X)
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
    !Call PrntFitPars(X)
    !write(*,*) 'X2Source:',PeakNrTotal
    If(N_FitPar .gt. N_FitStatTim) then
        Do i=1 , N_FitPar-N_FitStatTim
            k=FitParam(N_FitStatTim+i) ! coordiate, x,y,z
            if(k.eq.4) then
                Do i_Peak=1,PeakNrTotal
                    RefAntErr(i_Peak) = X( XIndx(i,i_Peak) )
                enddo
            else
                Do i_Peak=1,PeakNrTotal
                    SourcePos(k,i_Peak) = X( XIndx(i,i_Peak) )
                enddo
            endif
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
               If(Polariz) Then ! should agree with convention for calculating the Jacobian for fitting
                 Do j= 1,(Tot_UniqueAnt(k)-Tot_UniqueAnt(k-1))/2  ! there is same fit_timeoffset for even and odd
                    ! dt=Fit_TimeOffsetAnt(Tot_UniqueAnt(k-1)+2*j-1)-Fit_TimeOffsetAnt(Tot_UniqueAnt(k-1)+2*j)
                    Fit_TimeOffsetAnt(Tot_UniqueAnt(k-1)+2*j-1)=X(X_Offset(i)+j-1)
                    Fit_TimeOffsetAnt(Tot_UniqueAnt(k-1)+2*j)=X(X_Offset(i)+j-1)
                    !write(2,*) 'X2Source', i,k,j,X(X_Offset(i)+j-1), X_Offset(i)
                 Enddo
               Else
                 Do j= 1,Tot_UniqueAnt(k)-Tot_UniqueAnt(k-1)
                    Fit_TimeOffsetAnt(Tot_UniqueAnt(k-1)+j)=X(X_Offset(i)+j-1)
                 Enddo
               EndIf
            Else
              Fit_TimeOffsetStat(k) = X(X_Offset(i))
            endif
            !write(2,*) 'FitOffset',i,k,Unique_StatID(k),Fit_TimeOffsetStat(k)
        Enddo  ! i=1, N_FitStatTim
    endif
    !
    !Call PrntFitPars(X)
    !write(2,*) 'x2s;Fit_TimeOffsetAnt',Fit_TimeOffsetAnt(150:180)
    !write(2,*) 'Fit_TimeOffsetStat=',Fit_TimeOffsetStat(1:10)
    return
End Subroutine X2Source
!===================================
Subroutine Find_unique_StatAnt()
!  Update list of unique stations and antennas.
    use FitParams, only : Fit_TimeOffsetStat, Fit_TimeOffsetAnt
    use DataConstants, only : Station_nrMax, Ant_nrMax, ChunkNr_dim
    use DataConstants, only : Polariz
    !use ThisSource, only : Fit_TimeOffsetStat, Unique_StatID
    use Chunk_AntInfo, only : Ant_Stations, Ant_nr, Ant_IDs
    use Chunk_AntInfo, only : Unique_StatID, Nr_UniqueStat, Unique_SAI, Nr_UniqueAnt, Tot_UniqueAnt
    use StationMnemonics, only : Statn_ID2Mnem
    use unque, only : unique
    Implicit none
    Integer :: vec(Ant_nrMax*(ChunkNr_dim+1))
    integer, save ::  i_chunk, k, StatID, i_stat, i_ant
    integer,dimension(:),allocatable :: Station_IDs  ! scratch array
    Logical, save :: First=.true.
    !
    !write(2,*) 'ChunkNr_dim=',ChunkNr_dim
    vec=0
    k=0
    If(Nr_UniqueStat.gt.0) then
      Vec(1:Nr_UniqueStat)=Unique_StatID(1:Nr_UniqueStat)
      k=Nr_UniqueStat
    Endif
    Do i_chunk=1, ChunkNr_dim
        !write(2,*) 'ant#=',Ant_nr(i_chunk),k
        vec(k+1:k+Ant_nr(i_chunk))=Ant_Stations(1:Ant_nr(i_chunk),i_chunk)
        k=k+Ant_nr(i_chunk)
    enddo
    Call unique(vec,Station_IDs)
    Nr_UniqueStat=size(Station_IDs)-1
    Write(2,*) 'Find_unique_StatAnt:Nr_UniqueStat=',Nr_UniqueStat
    Unique_StatID(1:Nr_UniqueStat)=Station_IDs(2:Nr_UniqueStat+1)  ! get rid of leading 0
    Deallocate(Station_IDs)
    !
    !Write(2,*) 'Nr_UniqueStat1=',Nr_UniqueStat, Unique_StatID(1:Nr_UniqueStat)
    !  Find Unique antenna SAI-numbers
    k=0
    If(Nr_UniqueAnt.gt.0) then
      Vec(1:Nr_UniqueAnt)=Unique_SAI(1:Nr_UniqueAnt)
      k=Nr_UniqueAnt
    Endif
    vec(k+1:Ant_nrMax*(ChunkNr_dim+1))=0
    Do i_chunk=1, ChunkNr_dim
      If(Polariz) Then  ! keep only even-odd pairs
         Do i_ant=1,Ant_nr(i_chunk)
            If((MOD(Ant_IDs(i_ant,i_chunk),2).eq.0) .and. ((Ant_IDs(i_ant,i_chunk)+1) .eq. Ant_IDs(i_ant+1,i_chunk))) then
              vec(k+1)= 1000*Ant_Stations(i_ant,i_chunk)+Ant_IDs(i_ant,i_chunk)
              vec(k+2)=vec(k+1)+1
              k=k+2
              !If(Ant_Stations(i_ant,i_chunk).eq.7) Write(2,*) 'Find_unique_StatAnt',k, vec(k-1),vec(k)
            EndIf
         Enddo
      Else
        vec(k+1:k+Ant_nr(i_chunk))= &  !  How can this work? should this rather be vec(k+1:k+Ant_nr(i_chunk)); it does work however!
            1000*Ant_Stations(1:Ant_nr(i_chunk),i_chunk) + Ant_IDs(1:Ant_nr(i_chunk),i_chunk)
        !write(2,*) 'Find_unique, CHECK!!!!', k, Ant_nr(i_chunk), Ant_IDs(1,i_chunk),'--',Ant_IDs(Ant_nr(i_chunk),i_chunk), &
        !       ' ; ',vec(k+1),'--',vec(Ant_nr(i_chunk)),'--', vec(k+Ant_nr(i_chunk))
        k=k+Ant_nr(i_chunk)
        !write(2,*) 'vec(k+1:k+Ant_nr(i_chunk)):',k,Ant_nr(i_chunk)
      Endif
    enddo
    Call unique(vec,Station_IDs)
    Nr_UniqueAnt=size(Station_IDs) -1
    if(Ant_nrMax.lt.Nr_UniqueAnt) then
      write(2,*) 'max number of unique antennas is exceeded',Nr_UniqueAnt,Ant_nrMax
      Nr_UniqueAnt=Ant_nrMax
    endif
    Unique_SAI(1:Nr_UniqueAnt)=Station_IDs(2:Nr_UniqueAnt+1)  ! get rid of leading 0
    Deallocate(Station_IDs)
    !
    !Write(2,*) 'Nr_UniqueStat1=',Nr_UniqueStat, Unique_SAI(1:5)
    !Write(2,*) 'Nr_UniqueStat2=',Nr_UniqueStat, Unique_SAI(56:68)
    !Write(2,*) 'Find_unique_StatAnt',First, ChunkNr_dim
    If(First ) then
      First=.false.
      Write(2,"(A,I3,2x,50(1x,I5))") 'Nr_UniqueStat=',Nr_UniqueStat, Unique_StatID(1:Nr_UniqueStat)
      write(2,"(A,I4,2x,50(1x,A5))") 'Nr_UniqueAnt=',Nr_UniqueAnt, (Statn_ID2Mnem(Unique_StatID(k)),k=1,Nr_UniqueStat)
    Endif
    !Do k=1,Nr_UniqueAnt/10 +1
    !write(2,*) k,Unique_SAI(10*(k-1)+1:10*k)
    !enddo
    !
    Tot_UniqueAnt(0)=0
    StatID=Unique_SAI(1)/1000
    i_stat=1
    Do k=1,Nr_UniqueAnt
        If(Unique_SAI(k)/1000 .ne. StatID) then
            StatID=Unique_SAI(k)/1000
            i_stat=i_stat+1
        endif
        Tot_UniqueAnt(i_stat)=k       ! highest ranknr of unique_antenna for station Unique_StatID(i_stat)
    Enddo
    !
    Fit_TimeOffsetStat(:)=0.
    Fit_TimeOffsetAnt(:)=0.
    !FitOffset_SAI(:)=Unique_SAI        ! SAI = Station+Antenna_ID = station_ID*1000 + Ant_nr
    !
    Return
    !
    Do i_stat=1,Nr_UniqueStat
        write(2,*) i_stat,Tot_UniqueAnt(i_stat),Unique_StatID(i_stat),Unique_SAI(Tot_UniqueAnt(i_stat-1)+1:Tot_UniqueAnt(i_stat))
    enddo
    !
    Return
End Subroutine Find_unique_StatAnt
!===================================
Integer Function XIndx(i,i_peak)
    use FitParams, only : FitParam, N_FitStatTim, Fit_PeakNrTotal, Nr_TimeOffset
    use DataConstants, only : PeakNr_dim, RunMode
    !use FittingParameters
    use ThisSource, only : PeakNrTotal, PeakPos, Dual, Unique_pos
    use unque, only : unique
    Implicit none
    logical, save :: First=.true.
    integer, intent(in) :: i, i_Peak ! i stands for the parameter label as in FitParam(N_FitStatTim+i)
    integer :: k, i_pk
    integer,dimension(:),allocatable :: vec_unique
    !integer, save :: Unique_pos(1:PeakNr_dim)
    !Integer :: i_chunk, i_eo, k
    If(first) then
       If(Dual) then
         Call unique(PeakPos(1:PeakNrTotal),vec_unique)
         Fit_PeakNrTotal=size(vec_unique)
         Unique_pos(1:Fit_PeakNrTotal)=vec_unique(1:Fit_PeakNrTotal)
         Deallocate(vec_unique)
         !write(2,*) 'Fit_PeakNrTotal=',Fit_PeakNrTotal
       else
         Fit_PeakNrTotal=PeakNrTotal
       endif
       First=.false.
    endif
    i_pk=i_Peak
    k=0
    If(Dual .and. (FitParam(N_FitStatTim+i).ne.4) ) then
        Do k=1,Fit_PeakNrTotal
            If(Unique_pos(k) .eq. PeakPos(i_Peak)) then
                i_pk=k
                exit
            endif
        enddo
        If(RunMode.eq.3) i_pk=1
        k=0
        If(i.gt.1) then
            If(any(FitParam(N_FitStatTim+1:N_FitStatTim+i-1)==4)) k=PeakNrTotal-Fit_PeakNrTotal
        endif
    endif
    !
    !i_eo=Peak_eo(i_peak)
    !i_chunk=ChunkNr(i_peak)
    !k=i_peak - TotPeakNr(i_eo,i_chunk) + PeakNr(i_eo,i_chunk)       ! number of peaks for this (i_eo,i_chunk)
    !write(*,*) 'XIndx', FitParam(N_FitStatTim+i),i,Fit_PeakNrTotal,k,i_Peak, i_Pk,Dual
    XIndx = Nr_TimeOffset-1 +(i-1)*Fit_PeakNrTotal+k + i_Pk
    !write(2,*) 'XIndx', Nr_TimeOffset,i , Fit_PeakNrTotal,k, i_Pk
    !
End Function XIndx
!++++++++++++++++++++++++++++++++++++++++++
