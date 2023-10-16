Module ThisSource
    use Chunk_AntInfo, only : Ant_nrMax, Station_nrMax
    use constants, only : dp
    Integer, parameter :: Safety=10
    Real(dp) :: SourcePos(3)
    Real(dp) :: t_CCorr(-Safety:Safety), CCorr(-Safety:Safety,Ant_nrMax), CCorr_pp(-Safety:Safety,Ant_nrMax)
    Integer :: Sample_Offset(1:Ant_nrMax), CorrAntNrs(1:Ant_nrMax), Nr_Corr
    Real(dp) :: FineOffset_Station(1:Station_nrMax)
    Integer :: FineOffset_StatID(1:Station_nrMax)
End Module ThisSource
