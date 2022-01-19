!    Include 'InterferometryOptSbRtns-v18d.f90'
    !
Subroutine EI_Run
    !+++++++++ What is used:
    !- Loop over antennas
    ! - make frequency spectra
    !- End antennaLoop
    !- Loop over all pixels
    !  - calculate phase shift Each antenna
    !   - Loop frequency
    !    - Sum(antenna) spectra
    ! - End pixelloop
    !+++++++++ Alternative:
    !- Loop over antennas
    ! - make frequency spectra
    ! - Loop over all pixels
    !    - calculate phase shift
    !    - Sum frequency: update spectra with this antenna
    ! - End pixelLoop
    !- End antennaLoop
    !++++++++++++++++++++++
   Use Interferom_Pars, only : polar, N_pix, d_loc, CenLocPol, CenLoc, PixLoc, Nr_IntFerMx, Nr_IntFerCh, IntFer_ant, IntfNuDim
   Use Interferom_Pars, only : AveInten, AveIntenE, AveIntenN, RimInten, MaxIntfInten, SumWindw !, MaxIntfIntenLoc
   Use Interferom_Pars, only : MaxSlcInten, MaxSmPow, IntPowSpec !
   Use Interferom_Pars, only : i_chunk, Nr_IntFerMx, Nr_IntFerCh
   Use Interferom_Pars, only : PixPowOpt
   use constants, only : dp, ci, pi !, Sample, Refrac, c_mps
   use StationMnemonics, only : Statn_ID2Mnem
   use Chunk_AntInfo, only : SignFlp_SAI, PolFlp_SAI, BadAnt_SAI
   !use GLEplots, only : GLEplotControl
   Use Interferom_Pars, only : Alloc_EInterfImag_Pars, DeAlloc_EInterfImag_Pars
   use mod_test
   Implicit none
   integer :: i_eo, i, Nr_IntFer !,j, i_s, i_loc(1), i_ant, j_corr, Station, j_IntFer, i_nu, nxx
   Real(dp) :: PixLocPol(1:3) !, SmPow
   Integer :: Date_T(8)
   integer :: i_N, i_E, i_h !, N_h,  N_hlow, N_Eup, N_Elow
   Complex(dp), allocatable :: CMCnu(:,:,:), CMTime_Pix(:,:)
   !
   CALL DATE_AND_TIME (Values=DATE_T)
   WRITE(2,"(1X,I2,':',I2.2,':',I2.2,'.',I3.3,A,i1)") (DATE_T(i),i=5,8), ' start Interferometry'
   WRITE(*,"(1X,I2,':',I2.2,':',I2.2,'.',I3.3,A)") (DATE_T(i),i=5,8), achar(27)//'[0m'
   !
   i_chunk=1
   Call EISelectAntennas(i_chunk)  ! select antennas for which there is an even and an odd one.
   Nr_IntFer=Nr_IntFerCh(i_chunk)
   Call PixBoundingBox()      !  allocate array space, sets IntfNuDim
   Call Alloc_EInterfImag_Pars
   Allocate( CMCnu(Nr_IntFer,0:IntfNuDim,1:3),  CMTime_pix(1:2*IntfNuDim,1:3) )
   Call EISetupSpec(Nr_IntFer, IntfNuDim, CMCnu)
   !Nr_IntFerMx=Nr_IntFerCh(i_chunk) since there is only a single chunk
   !
   !------------------------------------
   ! main loop over Voxel location
   MaxIntfInten=0.
   MaxSlcInten(:)=0.
   AveInten(:,:)=0.
   AveIntenE(:,:)=0.
   AveIntenN(:,:)=0.
   RimInten(:,:)=0.
   MaxSmPow(:)=0.
   Select Case(PixPowOpt)
      Case (1)
         write(2,*) 'PixPowOpt=1 : sum all polarizations weighted with alpha to compensate A^-1 == intensity of F vector'
      Case(2)
         write(2,*) 'PixPowOpt=2 : Sum all three polarizations, thus including longitudinal with the full weight'
      Case Default
         write(2,*) 'PixPowOpt=0=default : sum two transverse polarizations only'
   End Select
   Do i_h= N_pix(3,1), N_pix(3,2) ! or distance
      CALL DATE_AND_TIME (Values=DATE_T)
      WRITE(2,"(1X,I2,':',I2.2,':',I2.2,'.',I3.3,A,i3)") (DATE_T(i),i=5,8), ' height',i_h
      WRITE(*,"(A,i3,1X,I2,':',I2.2,':',I2.2,'.',I3.3,A)") ' height',i_h,(DATE_T(i),i=5,8), achar(27)//'[0m'
      Flush(unit=2)
      !
      Do i_N=N_pix(1,1),N_pix(1,2)   ! or Phi  ! Loop over pixels
      Do i_E=N_pix(2,1),N_pix(2,2)   ! or elevation angle
         If(polar) then
            PixLocPol(1)=CenLocPol(1)+i_N*d_loc(1) ! azimuth
            PixLocPol(2)=CenLocPol(2)+i_E*d_loc(2) ! elevation angle
            PixLocPol(3)=CenLocPol(3)+i_h*d_loc(3)
            If(PixLocPol(2).gt.pi/2 .or. PixLocPol(2).lt.0.) stop 'Interferometry; theta problem'
            Call Pol2Carth(PixLocPol,PixLoc)
         else
            PixLoc(1)=CenLoc(1)+i_N*d_loc(1)
            PixLoc(2)=CenLoc(2)+i_E*d_loc(2)
            PixLoc(3)=CenLoc(3)+i_h*d_loc(3)
         Endif
         If(PixLoc(3).lt.0.) stop 'Interferometry; height problem'
         !
         ! sum frequency spectra
         !If(i_N.eq.0 .and. i_E.eq.0) then
            !CALL DATE_AND_TIME (Values=DATE_T)
            !WRITE(2,"(1X,I2,':',I2.2,':',I2.2,'.',I3.3,A,i2)") (DATE_T(i),i=5,8), ' A'
            !WRITE(*,"(A,1X,I2,':',I2.2,':',I2.2,'.',I3.3,A)") ' A:',(DATE_T(i),i=5,8), achar(27)//'[0m'
         !endif
         !
         Call EIEngine(Nr_IntFer, IntfNuDim, CMCnu, CMTime_pix)         !InterfEngine
         !
         Call EIAnalyzePixelTTrace(i_N, i_E, i_h, SumWindw, IntfNuDim, CMTime_pix)
         !
      EndDo ! i_E
      EndDo ! i_N
   EndDo   ! Loop over distances or heights
   !
   DeAllocate( CMCnu, CMTime_pix )
   !----------------------------------------
   WRITE(2,"(1X,I2,':',I2.2,':',I2.2,'.',I3.3,A,i3)") (DATE_T(i),i=5,8), ' analyze'
   WRITE(*,"(A,1X,I2,':',I2.2,':',I2.2,'.',I3.3,A)") ' analyze',(DATE_T(i),i=5,8), achar(27)//'[0m'
   i_eo=0
   Call OutputIntfPowrTotal(IntFer_ant(1,1), i_eo)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   If(IntPowSpec) Call OutputIntfPowrMxPos(i_eo)
   Call OutputIntfSlices(i_eo)
   !
   CALL DATE_AND_TIME (Values=DATE_T)
   WRITE(2,"(1X,I2,':',I2.2,':',I2.2,'.',I3.3,A,i1)") (DATE_T(i),i=5,8), ' end Interferometry',i_eo
   Call DeAlloc_EInterfImag_Pars
   Return
    !
End Subroutine EI_Run
! =========================
! ==========================
Subroutine Inverse33(A,B, eig, PolBasis)
   ! B= A^-1, from https://en.wikipedia.org/wiki/Eigenvalue_algorithm
   ! derivation, see https://www.geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf
   use constants, only : dp, pi
   Implicit none
   Real(dp), intent(in) :: A(3,3)
   Real(dp), intent(out) :: B(3,3),eig(3), PolBasis(3,3)
   integer, parameter :: c(3)=(/ 2,3,1 /) , r(3)=(/ 3,1,2 /)
   Real(dp) :: Det
   integer :: i,j,k
   Real(dp) :: p,q,s, phi
   s = A(1,2)**2 + A(1,3)**2 + A(2,3)**2 !Eigen values real symmetric 3x3 matrix
   if (s .eq. 0) then !      % A is diagonal.
      eig(1) = A(1,1)
      eig(2) = A(2,2)
      eig(3) = A(3,3)
   else
      q = (A(1,1)+A(2,2)+A(3,3))/3      !         % trace(A) is the sum of all diagonal values
      p = sqrt(((A(1,1) - q)**2 + (A(2,2) - q)**2 + (A(3,3) - q)**2 + 2 * s)/6.)
      s=(A(1,1)-q)*(A(2,2)-q)*(A(3,3)-q) + 2.*A(1,2)*A(2,3)*A(1,3) - &
         (A(1,1)-q)*A(2,3)*A(2,3) - (A(2,2)-q)*A(1,3)*A(1,3) - (A(3,3)-q)*A(1,2)*A(1,2)
      s=s/(2.*p*p*p)
      if (s .le. -1) then
         phi = pi / 3
      elseif (s .ge. 1) then
         phi = 0
      else
         phi = acos(s) / 3
      endif
      !write(2,*) 'phi',phi,s
      !% the eigenvalues satisfy eig3 <= eig2 <= eig1
      eig(1) = q + 2 * p * cos(phi)
      eig(3) = q + 2 * p * cos(phi + (2*pi/3))
      eig(2) = 3 * q - eig(1) - eig(3)   !  % since trace(A) = eig1 + eig2 + eig3
   endif
   !Write(2,*) 'eigenvalues', eig(:)
   !
   Call EigvecS33(A,eig,PolBasis)
   !
   Det=0.
   B(:,:)=0.
   Do i=1,3 ! get the Adjoint of A (transpose of subdeterminant matrices)
   Do j=1,3
   B(j,i)=A(c(i),c(j))*A(r(i),r(j))-A(c(i),r(j))*A(r(i),c(j))
   Enddo
   Enddo
   Do i=1,3
   Det=Det + A(i,1)*B(i,1)
   Enddo
   B(:,:)=B(:,:)/Det
   ! check
   !write(2,*) 'check inverse:'
   !Do i=1,3
   !Do j=1,3
   !   s=0
   !   Do k=1,3
   !      s=s+B(i,k)*A(k,j)
   !      !write(2,*) k,s,B(i,k),A(k,j)
   !   enddo
   !   write(2,*) i,j,B(i,j), A(i,j),s
   !Enddo
   !Enddo
   !
   Return
End Subroutine Inverse33
Subroutine EigvecS33(A,E,ev)
!  construct the eigen vectors ev(:,k) given the eigenvalues E(k) of the real symmetric matrix A.
!  Explicitly orthogonnormalize to the eigenvector for the smallest eigenvalue to minimize num. errors,
!  which could be necessary since EigVal could be like (1,1,0), i.e. the first two almost degenerate.
!  eigenvalues are assumed to be ordered, smallest last.
   use constants, only : dp
   Implicit none
   Real(dp), intent(in) :: A(3,3),E(3)
   Real(dp), intent(out) :: ev(3,3)
   Real(dp) :: s
   integer :: k
  ! check
   Do k=1,3
      ev(1,k)=1.
      ev(2,k)=(E(k)*A(2,3)-A(1,1)*A(2,3)+A(1,2)*A(1,3)) / (A(1,2)*A(2,3)+(E(k)-A(2,2))*A(1,3))
      ev(3,k)=(E(k)*A(2,3)-A(1,1)*A(2,3)+A(1,2)*A(1,3)) / (A(1,3)*A(2,3)+(E(k)-A(3,3))*A(1,2))
      s=SUM(ev(:,k)*ev(:,k))
      ev(:,k)=ev(:,k)/sqrt(s)
   Enddo
   s=SUM(ev(:,1)*ev(:,3))
   ev(:,1)=(ev(:,1)-s*ev(:,3))/sqrt(1-s*s)  ! ev 1 is orthonormalized to 3
   ev(1,2)=ev(2,1)*ev(3,3)-ev(3,1)*ev(2,3)  !  2= 1 x 3
   ev(2,2)=ev(3,1)*ev(1,3)-ev(1,1)*ev(3,3)  !  2= 1 x 3
   ev(3,2)=ev(1,1)*ev(2,3)-ev(2,1)*ev(1,3)  !  2= 1 x 3
   !write(2,*) 'check eigenvexcors (=0)', (SUM(ev(:,k)*ev(:,k))-1, k=1,3)
   !
   Return
End Subroutine EigvecS33
!++++++++++++++++++++++++++
Subroutine EISelectAntennas(i_chunk)
   ! Select antennas within a range of 'AntennaRange' from the core
   ! And put the condition the both dipoles are at same location
   !--------------------------------------------
   use constants, only : dp !, pi, Sample, Refrac, c_mps
   use Chunk_AntInfo, only : Ant_Stations, Ant_IDs, Ant_nr, Ant_pos !, Ant_RawSourceDist
   use Chunk_AntInfo, only : Unique_StatID, Nr_UniqueStat, Unique_SAI, Nr_UniqueAnt, Tot_UniqueAnt
   Use Interferom_Pars, only : IntFer_ant, Nr_IntFerMx, Nr_IntferCh, Nmax_IntFer ! the latter gives # per chunk
   use FitParams, only : AntennaRange
   use DataConstants, only : Station_nrMax, Ant_nrMax, ChunkNr_dim
   use StationMnemonics, only : Statn_ID2Mnem
    use unque, only : unique
   Implicit none
   Integer, intent(IN) :: i_chunk
   real ( kind = 8 ) :: Dist
   integer :: i_ant, Station, j_IntFer, k ! , MLoc, Mloc_all
    Integer :: vec(Ant_nrMax*(ChunkNr_dim+1))
    integer,dimension(:),allocatable :: Station_IDs  ! scratch array
   !
   !Do i_chunk=1,ChunkNr_dim
   If (.not. allocated(IntFer_ant)) Allocate( IntFer_ant(Nmax_IntFer,1:ChunkNr_dim), Nr_IntFerCh(1:ChunkNr_dim))
      Station=0
      j_IntFer=0
      Nr_UniqueStat=0
      Do i_ant=1,Ant_nr(i_chunk)         ! Select the antennas for Interferometry
         If(Ant_Stations(i_ant,i_chunk) .ne. Station) then
           Dist=sqrt(sum((Ant_pos(:,i_ant,i_chunk)-Ant_pos(:,1,i_chunk))**2))
           Station=Ant_Stations(i_ant,i_chunk)
           !write(2,*) Statn_ID2Mnem(Station), Ant_pos(:,i_ant,i_chunk)/1000.
           ! write(2,*) '-->',j_IntFer+1, Statn_ID2Mnem(Station), Ant_IDs(i_ant,i_chunk)
         endif
         !write(2,*) 'i_ant, Station, Dist',i_ant, Station, Dist
         If(Dist .gt. AntennaRange*1000.) cycle  ! Limit to antennas near the superterp
         if(mod(Ant_IDs(i_ant,i_chunk),2) .ne. 0) cycle       ! First select the even,
         if(Ant_IDs(i_ant+1,i_chunk) .ne. (Ant_IDs(i_ant,i_chunk)+1) ) cycle       ! then check if the odd antenna is there too
         j_IntFer=j_IntFer+1
         IntFer_ant(j_IntFer,i_chunk)=i_ant  ! keep track of the even ones only
         !write(2,*) j_IntFer, Statn_ID2Mnem(Station), Ant_IDs(i_ant,i_chunk)
         If(j_IntFer .ge. Nmax_IntFer) then
            Write(2,*) 'max nr of interferometry-antennas reached, station=',Station, Ant_IDs(i_ant,i_chunk),', at distance ', Dist
            exit
         Endif
      EndDo
      Nr_IntferCh(i_chunk)=j_IntFer
      If(j_IntFer .gt. Nr_IntFerMx) Nr_IntFerMx=j_IntFer  ! =max number for all chunks
   !EndDo
   !
    vec(:)=0
    k=0
    If(Nr_UniqueStat.gt.0) then
      Vec(1:Nr_UniqueStat)=Unique_StatID(1:Nr_UniqueStat)
      k=Nr_UniqueStat
      !write(2,*) 'How is this possible?????'
    Endif
    !Do i_chunk=1, ChunkNr_dim
        !write(2,*) 'ant#=',Ant_nr(i_chunk),k
     vec(k+1:k+Nr_IntferCh(i_chunk))=Ant_Stations(IntFer_ant(1:Nr_IntferCh(i_chunk),i_chunk),i_chunk)
     k=k+Nr_IntferCh(i_chunk)
    !enddo
    Call unique(vec,Station_IDs)
    Nr_UniqueStat=size(Station_IDs)-1
    Unique_StatID(1:Nr_UniqueStat)=Station_IDs(2:Nr_UniqueStat+1)  ! get rid of leading 0
    Deallocate(Station_IDs)
    !
    !Write(2,*) 'Nr_UniqueStat1=',Nr_UniqueStat, Unique_StatID(1:Nr_UniqueStat)
    !  Find Unique antenna SAI-numbers
    k=0
    vec(:)=0
    If(Nr_UniqueAnt.gt.0) then
      Vec(1:Nr_UniqueAnt)=Unique_SAI(1:Nr_UniqueAnt)
      k=Nr_UniqueAnt
      !write(2,*) 'This should not be possible!?!??!??'
    Endif
    !Do i_chunk=1, ChunkNr_dim
   Do j_IntFer=1,Nr_IntferCh(i_chunk)
   i_ant=IntFer_ant(j_IntFer,i_chunk)
        vec(k+1)= 1000*Ant_Stations(i_ant,i_chunk)+Ant_IDs(i_ant,i_chunk)
        vec(k+2)=vec(k+1)+1
        k=k+2
   Enddo
    !enddo
    Call unique(vec,Station_IDs)
    Nr_UniqueAnt=size(Station_IDs) -1
    if(Ant_nrMax.lt.Nr_UniqueAnt) then
      write(2,*) 'max number of unique antennas is exceeded',Nr_UniqueAnt,Ant_nrMax
      Nr_UniqueAnt=Ant_nrMax
    endif
    Unique_SAI(1:Nr_UniqueAnt)=Station_IDs(2:Nr_UniqueAnt+1)  ! get rid of leading 0
    Deallocate(Station_IDs)
    !
   Return
End Subroutine EISelectAntennas
!------------------------
!------------------------
Subroutine EIPixBoundingBox()
   ! Check max time shift for neighboring pixels; should ideally differ by 1 sample for neighboring pixels
   ! Calculate approximate bounding box of pixel image
   ! Setup arrays
   !  Needs
   !     call SelectIntfAntennas  first
   ! ------------------------------------
   use constants, only : dp !, pi, Sample, Refrac, c_mps
   use DataConstants, only : Time_Dim
   use Chunk_AntInfo, only : Ant_pos, Ant_RawSourceDist  ! CTime_spectr, Ant_Stations, Ant_IDs, Ant_nr,
   Use Interferom_Pars, only : i_chunk, IntFer_ant, Nr_IntFerCh
   Use Interferom_Pars, only : CenLoc, CenLocPol, Polar, d_loc, Diff, N_pix
   Use Interferom_Pars, only : xMin, xMax, yMin, yMax, zMin, zMax
   Use Interferom_Pars, only : SumStrt, SumWindw, IntfDim, IntfNuDim, IntfLead, IntfBase
   Implicit none
   Real(dp) :: PixLoc(1:3), PixLocPol(1:3), t_shft, RDist
   Real(dp) :: t_n, t_p
   Integer :: i, i_ant, j_IntFer
   !
   PixLoc(:)=CenLoc(:)
   xMin=+99999 ; xMax=-99999; yMin=+99999; yMax=-99999; zMin=+99999; zMax=-99999 !; tMin, tMax,
   Do i=1,3
      If(polar) then
         PixLocPol(:)=CenLocPol(:)
         PixLocPol(i)=CenLocPol(i)+d_loc(i)
         Call Pol2Carth(PixLocPol,PixLoc)
      else
         PixLoc(:)=CenLoc(:)
         PixLoc(i)=CenLoc(i)+d_loc(i)
      Endif
      !write(2,*) 'PixLoc-carthesian',PixLoc(:)
      t_n=0.   ; t_p=0
      Do j_IntFer=1,Nr_IntFerCh(i_chunk)   ! Loop over selected antennas
         i_ant=IntFer_ant(j_IntFer,i_chunk)
         Call RelDist(PixLoc(1),Ant_pos(1,i_ant,i_chunk),RDist)
         t_shft=Rdist - Ant_RawSourceDist(i_ant,i_chunk)
         if(t_shft.lt.t_n) t_n=t_shft
         if(t_shft.gt.t_p) t_p=t_shft
         !write(2,*) i_ant,'t_shft',t_shft, Rdist, Ant_RawSourceDist(i_ant,i_chunk)
      EndDo ! j_IntFer
      write(2,"(A,i2,f7.2,f6.2,A,f9.5,A,3f7.2,A)") 'min & max time shift [samples] per pixel',i, t_n, t_p, &
         ', for d(i)=',d_loc(i),', d(N,E,h)=',Pixloc(:)-CenLoc(:),'[m]'
      diff(i) =t_p
      xMin=min(xmin,CenLoc(2)+N_pix(i,1)*abs(CenLoc(2)-PixLoc(2))); xMax=max(xMax,CenLoc(2)+N_pix(i,2)*abs(CenLoc(2)-PixLoc(2)))
      yMin=min(ymin,CenLoc(1)+N_pix(i,1)*abs(CenLoc(1)-PixLoc(1))); yMax=max(yMax,CenLoc(1)+N_pix(i,2)*abs(CenLoc(1)-PixLoc(1)))
      zMin=min(zmin,CenLoc(3)+N_pix(i,1)*abs(CenLoc(3)-PixLoc(3))); zMax=max(zMax,CenLoc(3)+N_pix(i,2)*abs(CenLoc(3)-PixLoc(3)))
      If(-t_n.gt.t_p) diff(i) =-t_n
      PixLoc(i)=CenLoc(i)  ! set back to center for Cartesian coordinates
   Enddo ! i
   If(zMin.lt.0.) zMin=0.
   !
   !-----------------------------------------------
   !  setup arrays
   IntfLead=10+(diff(1)*N_pix(1,2)+diff(2)*N_pix(2,2)+diff(3)*N_pix(3,2)) ! Estimate lead and trail buffer length,
   IntfNuDim=SumWindw/2 + IntfLead +1 ! Length of frequency trace used for interference
   i = shiftr( IntfNuDim, 9 )+1
   IntfNuDim=shiftl( i, 9 )  ! Make sure it approches some multiple of 2
   !write(2,"(A,o10,o10,i7)") 'i:',i,IntfNuDim,IntfNuDim
   IntfDim=2*IntfNuDim ! Length of time trace used for interference
   IntfLead=(IntfDim-SumWindw)/2
   IntfBase=SumStrt - IntfLead
   !write(*,"(A,o10,b20)") '(lead,dim):', IntfLead, IntfDim
   If((IntfBase+IntfDim .gt. Time_dim) .or. (IntfBase .lt. 1)) then
      write(2,*) 'dimensions for interferometry out of range: ', IntfBase, IntfBase+IntfDim, Time_dim
      write(2,*) 'input for interferometry: ', SumStrt, SumWindw
      stop 'Intferom dimension problem'
   Endif
   !
   !write(2,*) diff(:)
   !
   Return
End Subroutine EIPixBoundingBox
! =========================
!-----------------------------------------------
Subroutine EISetupSpec(Nr_IntFer, IntfNuDim, CMCnu)
   ! CMCnu is the frequency dependent Current Moment Contribution of each antenna pair
   ! Needs:
   !   ?
   !--------------------------------------------
   use constants, only : dp, pi, Sample, Refrac, c_mps
   use DataConstants, only : Time_Dim, Cnu_dim
   use Chunk_AntInfo, only : CTime_spectr, Ant_Stations, Ant_IDs, Ant_nr, Ant_pos !, Ant_RawSourceDist
   use Chunk_AntInfo, only : Start_time, TimeBase
   use Chunk_AntInfo, only : NormOdd, NormEven !Powr_eo,NAnt_eo
   use FFT, only : RFTransform_CF !, RFTransform_CF2CT
   use StationMnemonics, only : Statn_ID2Mnem
   Use Interferom_Pars, only : IntfBase, SumStrt, IntfLead, SumWindw  !, IntfPhaseCheck
   Use Interferom_Pars, only : i_chunk, IntFer_ant, Nr_IntFerCh, Nr_IntFerMx
   Use Interferom_Pars, only : IntFer_ant
   Use Interferom_Pars, only : CenLoc, t_shft, t_offsetPow, tMin, tMax
   Use Interferom_Pars, only : alpha, PolBasis !, PowerScale
   use AntFunCconst, only : Freq_min, Freq_max,Ji_p0,Ji_t0,Ji_p1,Ji_t1, Gain ! J_0p,J_0t,J_1p,J_1t,
   Implicit none
   Integer, intent(in) :: Nr_IntFer, IntfNuDim
   Complex(dp), intent(out) :: CMCnu(1:Nr_IntFer,0:IntfNuDim,1:3)
   integer :: i_ant, j_IntFer, MLoc, Mloc_all, i, j, i_freq, i_nu, IntfDim, inu1, inu2
   !Real(dp) :: angle, angle_all, HorAng, AziAng, D, RD,  AbsAmp, NRM, SmPow  ! elevation
   !Complex(dp) :: CNu_s(0:Cnu_dim)  ! CTime_s(1:Time_dim),
   !
   !Complex(dp) :: Ji_p0(Freq_min:Freq_max),Ji_t0(Freq_min:Freq_max),Ji_p1(Freq_min:Freq_max),Ji_t1(Freq_min:Freq_max)
   Real(dp) :: Cur2E(1:3,1:3), E2Cur(1:3,1:3)
   Real(dp) :: Ras(1:3), Vec_p(1:3), Vec_t(1:3), Aip(1:3), Ait(1:3)
   Real(dp) :: HorDist, Thet_r ,Phi_r, dfreq, D, W, dnu, nu
   Real(dp) :: AntSourceD(1:Nr_IntFer), Weight(1:Nr_IntFer), thet_d(1:Nr_IntFer), Phi_d(1:Nr_IntFer)
   Real(dp) :: Rtime(1:2*IntfNuDim), AntennaNorm=1.! e-4
   Complex(dp) :: Cnu0(0:IntfNuDim), Cnu1(0:IntfNuDim), Sp, St
   Real(dp) :: AntW(2), NRM
   Character(len=50) :: FMT_A
   !Logical :: GainFactor=.false.
   !Logical :: GainFactor=.true.
   !
   ! Get antenna function parameters
   IntfDim=2*IntfNuDim
   !Nr_IntFer=Nr_IntFerCh(i_chunk)
   !
   Cur2E(:,:)=0.
   write(2,*) 'Nr_IntFer=',Nr_IntFer
   Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
      i_ant=IntFer_ant(j_IntFer,i_chunk)  ! Refers to even antenna
      !write(2,*) i_ant, j_IntFer, Ant_pos(:,i_ant,i_chunk)
      Ras(1)=(CenLoc(1)-Ant_pos(1,i_ant,i_chunk))/1000.
      Ras(2)=(CenLoc(2)-Ant_pos(2,i_ant,i_chunk))/1000.
      Ras(3)=(CenLoc(3)-Ant_pos(3,i_ant,i_chunk))/1000.
      HorDist= Ras(1)*Ras(1) + Ras(2)*Ras(2)  ! =HYPOT(X,Y)
      AntSourceD(j_IntFer)=sqrt(HorDist + Ras(3)*Ras(3))
      HorDist=sqrt( HorDist ) ! =HYPOT(X,Y)
      Thet_r=atan(HorDist,Ras(3))  ! Zenith angle
      Phi_r=atan2( Ras(2) , Ras(1) )  ! \phi=0 = north
      Thet_d(j_IntFer)=Thet_r*180/pi
      Phi_d(j_IntFer)=Phi_r*180/pi
      Weight(j_IntFer)= AntSourceD(1)*AntSourceD(1)/Nr_IntFer  ! to normalize at unity for the reference antenna and A=(1,1,0)
      ! Calculate Cur2E matrix
      If((AntSourceD(j_IntFer)/AntSourceD(1)).gt.1. ) Then  ! The factor between E-field and signal is approx square of this
         Weight(j_IntFer)=Weight(j_IntFer) ! *exp( -(AntSourceD(j_IntFer)/AntSourceD(1))**2 )
      EndIf
      !Write(2,*) j_IntFer,Ant_IDs(i_ant,i_chunk), Statn_ID2Mnem(Ant_Stations(i_ant,i_chunk)), &
      !      AntSourceD(j_IntFer),Thet_d(j_IntFer),Phi_d(j_IntFer)
      D=AntSourceD(j_IntFer)*AntSourceD(j_IntFer)
      w=Weight(j_IntFer)/D ! to account for 1/R_{as}^2 in formula for A
      Do i=1,3
         Cur2E(i,i)=Cur2E(i,i)+w
         Do j=1,3
            Cur2E(i,j)=Cur2E(i,j)-w*Ras(i)*Ras(j)/D  ! to normalize to \hat{r}_{as}
         Enddo
      Enddo
   Enddo
   !
   !write(2,*) 'Trace A=', Cur2E(1,1)+Cur2E(2,2)+Cur2E(3,3)
   Call Inverse33(Cur2E,E2Cur, alpha, PolBasis)  ! PolBasis are normalized to unity
   FMT_A="(A,F7.3,A,'(',2(F7.3,','), F7.3 ,')')"
   write(2,FMT_A) 'a1:', alpha(1),', ev1:', PolBasis(:,1)
   write(2,FMT_A) 'a2:', alpha(2),', ev2:', PolBasis(:,2)
   FMT_A="(A,F7.3,A,'(',2(F7.3,','),F7.3,')',A,F7.3)"
   write(2,FMT_A) 'a3:', alpha(3),', ev3:', PolBasis(:,3), ', Trace A=',alpha(1)+alpha(2)+alpha(3)
   !
   !Call AntFun_Inv(thet_d(1),Phi_d(1)) ! sets ,Ji_p0,Ji_t0,Ji_p1,Ji_t1; Inverse Jones; gives _p & _t polarized fields when contracted with (0.1) voltages
   !w=SUM(ABS(Ji_p0(Freq_min:Freq_max))**2)+SUM(ABS(Ji_p1(Freq_min:Freq_max))**2)+ &
   !   SUM(ABS(Ji_t0(Freq_min:Freq_max))**2)+SUM(ABS(Ji_t1(Freq_min:Freq_max))**2)
   !PowerScale=1.e4/W  ! 100 to off-set the value at moderate height
   !write(2,*) 'weight (in power) inverse Jones=',PowerScale,thet_d(1),Phi_d(1),Freq_min,Freq_max
   !Flush(unit=2)
   !
   !PowerScale=1.
   !NormEven=sqrt(PowerScale*2.*Powr_eo(0)/(Powr_eo(0)+Powr_eo(1)))/100.  ! to undo the factor 100 (and more) that was introduced in antenna-read
   !NormOdd=sqrt(PowerScale*2.*Powr_eo(1)/(Powr_eo(0)+Powr_eo(1)))/100.  ! to undo the factor that was introduced in antenna-read
   write(2,*) 'Amplitude NormEven,odd',NormEven,NormOdd,', ratio=',NormEven/NormOdd
   Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
      Call AntFun_Inv(thet_d(j_IntFer),Phi_d(j_IntFer)) ! sets ,Ji_p0,Ji_t0,Ji_p1,Ji_t1; Inverse Jones; gives _p & _t polarized fields when contracted with (0.1) voltages
      Thet_r = Thet_d(j_IntFer)*pi/180  ! Zenith angle
      Phi_r  = Phi_d(j_IntFer)*pi/180   ! \phi=0 = north
      Vec_p(1)=sin(Phi_r)              ; Vec_p(2)=-cos(Phi_r)        ; Vec_p(3)=0.
      Vec_t(1)=-cos(Thet_r)*Vec_p(2)  ; Vec_t(2)=cos(Thet_r)*Vec_p(1) ; Vec_t(3)=-sin(Thet_r)
      !
      !Write(2,*) 'Vec_p(:)',Vec_p(:),Vec_t(:)
      Do i=1,3  ! current moment is calculated in basis given by PolBasis
         Aip(i)=SUM(PolBasis(:,i)*Vec_p(:))/alpha(i)
         Ait(i)=SUM(PolBasis(:,i)*Vec_t(:))/alpha(i)
      Enddo
      !
      i_ant=IntFer_ant(j_IntFer,i_chunk)  ! Refers to even antenna
      RTime(:)=REAL(CTime_spectr(IntfBase+1:IntfBase+IntfDim,i_ant,i_chunk))*NormEven
      Call RFTransform_CF(RTime,Cnu0(0))
      RTime(:)=REAL(CTime_spectr(IntfBase+1:IntfBase+IntfDim,i_ant+1,i_chunk))*NormOdd
      Call RFTransform_CF(RTime,Cnu1(0))
      !
      AntW(:)=0.
      w=AntennaNorm*Weight(j_IntFer)/AntSourceD(j_IntFer)
      dnu=100./IntfNuDim   ! [MHz] Jones matrix is stored on 1MHz grid
      Do i_nu=0,IntfNuDim
         nu=i_nu*dnu
         i_freq=Int(nu)
         dfreq=nu-i_freq
         If((i_freq .lt. Freq_min) .or. (i_freq .ge. Freq_max) ) then
            CMCnu(j_IntFer,i_nu,1:3)=0.
         Else
            Sp=((1.-dfreq)*Ji_p0(i_freq) + dfreq*Ji_p0(i_freq+1)) *Cnu0(i_nu) + &
               ((1.-dfreq)*Ji_p1(i_freq) + dfreq*Ji_p1(i_freq+1)) *Cnu1(i_nu)
            St=((1.-dfreq)*Ji_t0(i_freq) + dfreq*Ji_t0(i_freq+1)) *Cnu0(i_nu) + &
               ((1.-dfreq)*Ji_t1(i_freq) + dfreq*Ji_t1(i_freq+1)) *Cnu1(i_nu)
            CMCnu(j_IntFer,i_nu,:)=w*( Aip(:)* Sp + Ait(:)* St ) *Gain(i_freq)
         Endif
      Enddo
      !
   Enddo         ! loop antennas for Interferometry
   !
   !write(2,*) 'CMCnu(Nr_IntFer/2,1:2,IntfNuDim/2)', IntfNuDim/2,  IntFer_ant(Nr_IntFer/2),  CMCnu(Nr_IntFer/2,1,IntfNuDim/2) &
   !   , IntFer_ant(Nr_IntFer/2+1),CMCnu(Nr_IntFer/2+1,1,IntfNuDim/2)
   i_ant=IntFer_ant(1,1)
   t_shft=sqrt(SUM(CenLoc(:)*CenLoc(:)))*Refrac/c_mps ! in seconds due to signal travel distance
   write(2,*) 'Nr_IntFer:',Nr_IntFer,', ref. ant.= ', Ant_IDs(i_ant,i_chunk), Statn_ID2Mnem(Ant_Stations(i_ant,i_chunk)),&
      ', time difference with central pixel=',t_shft*1000.,'[ms]'
   !
   tMin=((Start_time(i_chunk)+SumStrt)*sample-t_shft)*1000.
   tMax=((Start_time(i_chunk)+SumStrt+SumWindw)*sample-t_shft)*1000.
   write(2,*) 'Source time window start@', tMin, '[ms], end@', tMax
   t_offsetPow=((Start_time(i_chunk)+SumStrt)*sample)*1000.-TimeBase
   Flush(unit=2)
   !
   Return
End Subroutine EISetupSpec
!-----------------------------------------------
! =========================
Subroutine EIEngine(Nr_IntFer, IntfNuDim, CMCnu, CMTime_pix)
   ! Phase shifts the Curr Mom Contr of the individual antennas to obtain Current Moment for a pixel in time
   !
   use constants, only : dp, ci, pi
   !use DataConstants, only : Time_Dim, Cnu_dim, DataFolder, OutFileLabel, Calibrations
   use Chunk_AntInfo, only : Ant_pos, Ant_RawSourceDist
   Use Interferom_Pars, only : i_chunk, IntFer_ant, IntfLead, PixLoc
   use AntFunCconst, only : Freq_min, Freq_max
   use FFT, only : RFTransform_CF2CT
   Implicit none
   Integer, intent(in) :: Nr_IntFer, IntfNuDim
   Complex(dp), intent(out) ::  CMTime_pix(1:2*IntfNuDim,1:3), CMCnu(Nr_IntFer,0:IntfNuDim,1:3)
   Integer :: i, i_ant, i_nu, j_IntFer,i_freq, inu1, inu2
   Complex(dp) :: Phase(Nr_IntFer), dPhase(Nr_IntFer), CMnu_pix(0:IntfNuDim,1:3)
   complex(dp), parameter :: ipi=ci*pi
   Real(dp) :: RDist, dt_AntPix, dnu
   !InterfEngine
   CMnu_pix(:,:)=0.
   dnu=100./IntfNuDim   ! [MHz] Jones matrix is stored on 1MHz grid
   inu1=Int(Freq_min/dnu)+1
   !inu2=Int(Freq_max/dnu-0.5)+1  ! set integration regime to frequencies within filter band
   inu2=Int(Freq_max/dnu)-1
   Do j_IntFer=1,Nr_IntFer   ! Loop over selected antennas
      i_ant=IntFer_ant(j_IntFer,i_chunk)
      Call RelDist(PixLoc(1),Ant_pos(1,i_ant,i_chunk),RDist)
      dt_AntPix=Rdist - Ant_RawSourceDist(i_ant,i_chunk)
      ! check dt_AntPix in range
      If(ABS(dt_AntPix).gt.IntfLead) then
         write(2,*) 'Warning, time-shift=',dt_AntPix,'[samples] out of range for pixel',PixLoc
         Return
      Endif
      dphase(j_IntFer)= exp(-ipi*dt_AntPix/IntfNuDim)
      Phase(j_IntFer)=exp(-ipi*dt_AntPix*inu1/IntfNuDim)
   EndDo ! j_IntFer
   !write(2,*) PixLoc(:),j_IntFer, dphase(j_IntFer-5:j_IntFer)
   !Phase(:)=1.d0
   ! Cnu produced in Subroutine SetupIntfSpec, should be modified for version B
   Do i_nu=inu1,inu2   ! Increment frequency spectrum with this antenna
      Do i=1,3
      CMnu_pix(i_nu,i)= SUM(phase(:)*CMCnu(:,i_nu,i))
      Enddo
      phase(:)=phase(:)*dphase(:)
   EndDo
   !
   Do i=1,3  ! convert to time
      Call RFTransform_CF2CT(CMnu_pix(0,i),CMTime_pix(1,i) )
   Enddo
   !
   Return
End Subroutine EIEngine
!-----------------------------------------------
Subroutine EIAnalyzePixelTTrace(i_N, i_E, i_h, SumWindw, IntfNuDim, CMTime_pix)
   ! Analyze the trace of this pixel
   !--------------------------------------------
   use constants, only : dp, pi !ci !, , Sample, Refrac, c_mps
   Use Interferom_Pars, only : i_chunk, Nr_IntFerMx, IntfLead, PixLoc, polar, N_pix
   Use Interferom_Pars, only : IntPowSpec, MaxSmPow, N_smth, smooth
   Use Interferom_Pars, only : MaxSmPowQ, MaxSmPowU, MaxSmPowV, MaxSmPowI3, MaxSmPowU1, MaxSmPowV1, MaxSmPowU2, MaxSmPowV2
   Use Interferom_Pars, only : AveInten, AveIntenE, AveIntenN, RimInten, MaxIntfInten, MaxIntfIntenLoc
   Use Interferom_Pars, only : SlcInten, NrSlices, SliceLen, MaxSlcInten, MaxSlcIntenLoc
   Use Interferom_Pars, only : NrPixSmPowTr, MaxSmPowGrd, PixSmPowTr, RefSmPowTr
   Use Interferom_Pars, only : SumStrt, CenLoc, StartTime_ms, t_shft, d_loc
   Use Interferom_Pars, only : alpha, PolBasis, PixPowOpt
   use FitParams, only : AntennaRange
   use DataConstants, only : DataFolder, OutFileLabel
   Implicit none
   Integer, intent(in) :: i_N, i_E, i_h, SumWindw, IntfNuDim
   Complex(dp), intent(in) :: CMTime_pix(1:2*IntfNuDim,1:3)
   integer :: i, j, i_s, i_eo, m, n  !_loc(1), i_ant, j_IntFer, MLoc, Mloc_all
   Real(dp) :: IntfInten, SmPow, PixelPower(1:SumWindw)
   Complex :: Stk(3,3)  !  3D Stokes
   !
   i_eo=0
   ! sum over 3 orientations of current moment
   !write(2,*) '1',CMTime_pix(IntfLead+1:IntfLead+5,1)
   !write(2,*) '2',CMTime_pix(IntfLead+1:IntfLead+5,2)
   !
   !SMpow=0.
   !PixPowOpt =0 : sum two transverse polarizations only
   !PixPowOpt =1 : sum all polarizations weighted with alpha to compensate A^-1 == intensity of F vector
   !PixPowOpt =2 : Sum all three polarizations, thus including longitudinal with the full weight
   Select Case(PixPowOpt)
      Case (1)
         Do j=1,SumWindw
         PixelPower(j)=SUM((ABS(CMTime_pix(IntfLead+j,:)))**2*alpha(:))
         Enddo
      Case(2)
         Do j=1,SumWindw
         PixelPower(j)=SUM((ABS(CMTime_pix(IntfLead+j,:)))**2)
         Enddo
      Case Default
         Do j=1,SumWindw
         PixelPower(j)=((ABS(CMTime_pix(IntfLead+j,1)))**2+(ABS(CMTime_pix(IntfLead+j,2)))**2) !/Nr_IntFer**2
         Enddo
   End Select
   !write(2,*) i_N, i_E, i_h, SMpow/SumWindw, PixLoc(:)
   If((i_N.eq.0) .and. (i_E.eq.0) .and. (i_h.eq.0)) Then
      i_s=0
      If(polar) i_s=1
      !write(2,*) 'SumStrt,SumWindw,',SumStrt,SumWindw,IntfLead
      OPEN(UNIT=29,STATUS='unknown',ACTION='WRITE',FILE=trim(DataFolder)//TRIM(OutFileLabel)//'_EISpec.csv')
         Write(29,"(f9.3,3(',',f8.3),',',f5.1,',',i2,',',F8.3,' !')") StartTime_ms, CenLoc/1000., AntennaRange, &
            i_s, t_shft*1000.
         write(29,"(g12.4,5(',',i5),',',i8,',',i5,6(',',f9.4),' !')") MaxIntfInten, N_pix(1,2), N_pix(2,1), N_pix(2,2), &
            N_pix(3,1), N_pix(3,2), SumStrt,SumWindw, d_loc(1:2)*180./pi ,d_loc(3), d_loc(1:3)
      Do i_s=1, SumWindw
         Do i=1,3
            Do j=1,3
               Stk(i,j)= CMTime_pix(IntfLead+i_s,i)*Conjg(CMTime_pix(IntfLead+i_s,j))
            Enddo
         Enddo
         SMpow=100./(Stk(1,1)+Stk(2,2))
         write(29,"(i6,',',g12.4,8(',',f7.1))") i_s, Real(Stk(1,1)+Stk(2,2)), &
            Real(Stk(1,1)-Stk(2,2))*SMpow, 2*Real(Stk(1,2))*SMpow, 2*Imag(Stk(1,2))*SMpow, Real(Stk(3,3))*SMpow ! &
!            , Real(Stk(3,3))*SMpow, &
!            Real(Stk(2,3))*SMpow, Real(Stk(3,1))*SMpow, &
!            Imag(Stk(2,3))*SMpow, Imag(Stk(3,1))*SMpow
      Enddo
      Close(Unit=29)
   EndIf
   !write(2,*) 'PixelPower',PixelPower(1:10)
   !
   !stop
   ! Fold intensity with smooth function and determine sequential positions at the peaks of the smoothing function
   If(IntPowSpec) Then
      Do i=1,NrPixSmPowTr ! N_smth+1,SumWindw-N_smth,N_smth
         i_s=1+i*N_smth
         SmPow=0.
         Do j=-N_smth,N_smth ! fold with smoothing function
            SmPow=SmPow+smooth(j)*PixelPower(i_s+j)
         Enddo
         !If( SmPow.gt. 0.6*RefSmPowTr(i) ) SmPow=0.
         PixSmPowTr(i,i_N, i_E, i_h)=SmPow
         If( SmPow.gt. MaxSmPow(i) ) Then
            MaxSmPow(i)=SmPow
            MaxSmPowGrd(1,i)=i_N
            MaxSmPowGrd(2,i)=i_E
            MaxSmPowGrd(3,i)=i_h
            ! Calculate and store stokes parameters for the maximum
            Stk(:,:)=0.
            Do j=-N_smth,N_smth ! fold with smoothing function
               !SmPow=SmPow+smooth(j)*PixelPower(i_s+j)
               Do m=1,3
                  Do n=1,3
                     Stk(m,n)= Stk(m,n) + smooth(j)*CMTime_pix(IntfLead+i_s+j,m)*Conjg(CMTime_pix(IntfLead+i_s+j,n))
                  Enddo
               Enddo
            Enddo
            MaxSmPowQ(i)=Real(Stk(1,1)-Stk(2,2))/SMpow
            MaxSmPowU(i)= 2*Real(Stk(1,2))/SMpow
            MaxSmPowV(i)=  2*Imag(Stk(1,2))/SMpow
            MaxSmPowI3(i)=  Real(Stk(3,3))/SMpow
            MaxSmPowU1(i)= 2*Real(Stk(1,3))/SMpow
            MaxSmPowV1(i)=  2*Imag(Stk(1,3))/SMpow
            MaxSmPowU2(i)= 2*Real(Stk(2,3))/SMpow
            MaxSmPowV2(i)=  2*Imag(Stk(2,3))/SMpow
      !write(2,*) 'MaxSmPow(i)',i, i_N, i_E, i_h, SmPow, N_smth, i_s
         EndIf
      Enddo
   EndIf
   ! Extract interesting numbers such as pixel intensity
   IntfInten=0.
   i=0
   SlcInten(:)=0.
   Do i_s=1,NrSlices
      Do j=1,SliceLen
         i=i+1
         SlcInten(i_s)=SlcInten(i_s) + PixelPower(i)
      Enddo
      SlcInten(i_s)=SlcInten(i_s)/SliceLen
      IntfInten = IntfInten + SlcInten(i_s)
      If(SlcInten(i_s) .gt. MaxSlcInten(i_s)) then  ! find max intensity location
         MaxSlcInten(i_s)=SlcInten(i_s)
         MaxSlcIntenLoc(i_s,:)=PixLoc(:)
      Endif
      SlcInten(i_s)=SlcInten(i_s)*SlcInten(i_s) ! /SliceLen
      AveInten(i_s,i_h)=AveInten(i_s,i_h) + SlcInten(i_s)
      AveIntenE(i_s,i_h)=AveIntenE(i_s,i_h) + i_E*SlcInten(i_s)
      AveIntenN(i_s,i_h)=AveIntenN(i_s,i_h) + i_N*SlcInten(i_s)
      If((i_E.eq.N_pix(2,1)) .or. (i_E.eq.N_pix(2,2)) .or. (i_N.eq.N_pix(1,1)) .or. (i_N.eq.N_pix(1,2))) Then ! rim of picture
         RimInten(i_s,i_h)=RimInten(i_s,i_h) + SlcInten(i_s)
      Endif
   Enddo
   IntfInten = IntfInten/NrSlices
!   If(IntfInten .gt. MaxIntfInten) then  ! find max intensity location
!      MaxIntfInten=IntfInten
!      MaxIntfIntenLoc(:)=PixLoc(:)
!   Endif
   !
   i=0  ! Store totals information
   PixSmPowTr(i,i_N, i_E, i_h)=IntfInten
   !write(2,*) i_N, i_E, i_h, IntfInten
   !stop
   If( IntfInten .gt. MaxSmPow(i) ) Then
      MaxSmPow(i)=IntfInten
      MaxSmPowGrd(1,i)=i_N
      MaxSmPowGrd(2,i)=i_E
      MaxSmPowGrd(3,i)=i_h
      MaxIntfInten=IntfInten
      MaxIntfIntenLoc(:)=PixLoc(:)
      !write(2,*) i_N, i_E, i_h, IntfInten, PixLoc(:)
   EndIf
   !
   Return
End Subroutine EIAnalyzePixelTTrace
!-----------------------------------------------
