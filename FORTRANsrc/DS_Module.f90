!    Include 'ConstantsModules.f90'
!-----------------------------------
Module TrackData
   Use constants, only : dp
   Character(len=100) :: PreDefTrackFile
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
contains
End Module TrackData
!===============================================
