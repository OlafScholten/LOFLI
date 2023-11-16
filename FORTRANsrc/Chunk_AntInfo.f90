Module Chunk_AntInfo
    use constants, only : dp
    Integer, parameter :: Time_dim=32768 ! 32768=2^15 ! 2048=2^11
    Integer, parameter :: Cnu_dim=Time_dim/2
    Integer, parameter :: Ant_nrMax=200, Station_nrMax=45
    Integer :: Ant_Stations(1:Ant_nrMax), Ant_IDs(1:Ant_nrMax)
    Integer :: Ant_nr, Ant_SamplOffSet(1:Ant_nrMax)
    Real(dp) :: Ant_pos(3,1:Ant_nrMax)
    Real(dp) :: Ant_RawSourceDist(1:Ant_nrMax) ! distances that have been taken into account already
    Complex(dp) :: CTime_spectr(1:Time_dim,1:Ant_nrMax)
    !Character(len=8) :: Ant_Mnems(1:Ant_nrMax)
End Module Chunk_AntInfo
