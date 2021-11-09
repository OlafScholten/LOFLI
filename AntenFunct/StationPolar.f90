Subroutine StatStokes(Station,DistMax)
!  Obtain average Stokes parameters for a whole station for all sources
    use DataConstants, only : StartTime_dim
    use Chunk_AntInfo, only : Ant_Stations, Ant_IDs, Ant_nr, Ant_pos, Nr_UniqueStat
    use ThisSource, only : Nr_Corr, Peakpos, PlotCCPhase, CCPhase_ave, Stat_pos, TotPeakNr
    use constants, only : dp
    Implicit none
    integer, intent(in) :: StatMax
    Real(dp), intent(in) :: DistMax
    real ( kind = 8 ) :: Dist
    integer :: i_ant, j_corr, i_eo, i_start, Station, i_Peak
    character(len=2) :: txt
    !
    Do i_start=1, StartTime_dim
        Do i_eo=0,1
            j_corr=0
            Do i_ant=1,Ant_nr(i_start)
                If(Ant_Stations(i_ant,i_start) .ne. Station) cycle
                if(mod(Ant_IDs(i_ant,i_start),2) .ne. i_eo) cycle       ! limit to odd antennas
                !
                Call GetCorrSingAnt( i_ant, J_Corr, i_eo, i_start)
            Enddo ! i_ant=1,Ant_nr(i_start)
            Nr_Corr(i_eo,i_start)=j_corr
            !
            Call AddTrace()
            !
        Enddo !     i_eo=0,1
    EndDo  !  i_start=1, StartTime_dim
End Subroutine StatStokes
!=====================================
Subroutine AddTrace()
    complex(dp), parameter :: tipi=2*ci*pi
   Call RFTransform_su(T2_dim)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   Do i_Peak=PeakNr1,TotPeakNr(i_eo,i_start)
   Cnu_s=0.
      Do j_corr=1,Nr_Corr(i_eo,i_start)
         RtMax = T_Offset(j_corr,i_Peak) + CCorr_max(j_corr,i_Peak)
         Sample_Offset = INT(RtMax) ! in units of sample size
         SubSample_Offset = RtMax - Sample_Offset ! in units of sample size
         !
         StLoc=PeakPos(i_Peak) + Sample_Offset - T2_dim/2 ! Check - sign 'Sample_Offset'
         RTTrace2(1:T2_dim)= Real(CTime_spectr(StLoc + 1:StLoc + T2_dim,i_ant,i_start))
         !Call RFTransform_CF(RTime1,Cnu1_s)
         Call RFTransform_CF(RTTrace2,Cnu2)
         !
         DO i=0,T2_dim/2
            Cnu_s(i)=Cnu_s(i) + Cnu2(i) * exp(-tipi*SubSample_Offset*i/T2_dim)/Nr_Corr(i_eo,i_start)
         ENDDO
         !
      Enddo ! j_corr=1,Nr_Corr(i_eo,i_start)
      Call RFTransform_CF2CT(Cnu_s,SumTrace(1,i_Peak) )
   Enddo  ! i_Peak
   Call DAssignFFT()                   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !
   Return
End Subroutine AddTrace
!=====================================
