!=================================
!=================================
!    Include 'MappingUtilitiesModules.f90'

Subroutine ITRF2LOFARConstruct()
    use constants, only : dp,pi,RotMat,center_CS002
    implicit none
!    real(dp), intent(out) :: RotMat(3,3) = reshape( (/ 0.8056597 ,  0.        ,  0.59237863 , &
!                              -0.05134149, -0.99623707,  0.06982658 , &
!                              -0.59014956,  0.08667007,  0.80262806 /), (/ 3,3/))
!'''
! 52.99006128232752 deg North and 6.350944049999953 deg East
!rotation_matrix = np.array([[-0., -0.8056597, 0.59237863],
!       [0.99623707, 0.05134149, 0.06982658],
!       [-0.08667007, 0.59014956, 0.80262806]])
!
!# Center of station CS002 in ITRF coordinates, substract this before rotating
!    real(dp) :: center_CS002(3)= (/ 3826577.5, 461021.3125, 5064893.0 /)
    real(dp) :: x,a,b,Rotmatphi(3,3),Rotmatth(3,3),LOFAR(3)
    real(dp) :: theta=(52.915141-90.d0)*pi/180., phi=6.8698134123520900*pi/180. ! corrects for omegaxomegaxr-->0.19878 deg
!    real(dp) :: theta=(52.729828630863530-90.d0)*pi/180., phi=6.8698134123520900*pi/180. ! CS002 in Z-direction

    integer :: i,j,k
!    RotMat = (/ 0.8056597 ,  0. , 0.59237863 , -0.05134149, -0.99623707,  0.06982658 ,  -0.59014956,  0.08667007,  0.80262806 /)
    Rotmatphi=0.
    RotMatphi(1,1)=cos(phi)  ; RotMatphi(1,2)=sin(phi)
    RotMatphi(2,1)=-sin(phi) ; RotMatphi(2,2)=cos(phi)
    RotMatphi(3,3)=1
    RotMatTh=0.
    RotMatTh(2,2)=1.
    RotMatTh(1,1)=-cos(theta)  ; RotMatTh(1,3)=-sin(theta)
    RotMatTh(3,1)=-sin(theta) ; RotMatTh(3,3)=cos(theta)
    Do i=1,3
        Do j=1,3
            x=0.
            Do k=1,3
                x=x+RotMatTh(i,k)*RotMatPhi(k,j)
            enddo
            RotMat(i,j)=x
        enddo
    enddo

    a=0.
    Do i=1,3
        a=a + center_CS002(i)**2
        x=0.
        Do j=1,3
            x=x + RotMat(i,j) * center_CS002(j)
        enddo
        LOFAR(i)=x
    enddo
!    write(2,*) 'cs002=',LOFAR
    b=center_CS002(1)**2 + center_CS002(2)**2
    a=sqrt(a)
    b=sqrt(b)
!    write(2,*) 'phi=',b,(center_CS002(2)/b),asin(center_CS002(2)/b) *180./pi
!    write(2,*) 'theta=',a,(center_CS002(3)/a),asin(center_CS002(3)/a) *180./pi,asin(b/a)*180./pi
!    write(2,*) 'RotMat=',RotMat
    return
End Subroutine ITRF2LOFARConstruct
!-------------------------------------
!-----------------------------------------
Subroutine ITRF2LOFAR(ITRF,LOFAR)
    use constants, only : dp,pi,RotMat,center_CS002
    implicit none
    real(dp), intent(in) :: ITRF(3)
    real(dp), intent(out) :: LOFAR(3)
!    real(dp) :: RotMat(3,3) = reshape( (/ 0.8056597 ,  0.        ,  0.59237863 , &
!                              -0.05134149, -0.99623707,  0.06982658 , &
!                              -0.59014956,  0.08667007,  0.80262806 /), (/ 3,3/))
!'''
! 52.99006128232752 deg North and 6.350944049999953 deg East
!rotation_matrix = np.array([[-0., -0.8056597, 0.59237863],
!       [0.99623707, 0.05134149, 0.06982658],
!       [-0.08667007, 0.59014956, 0.80262806]])
!
!# Center of station CS002 in ITRF coordinates, substract this before rotating
!    real(dp) :: center_CS002(3)= (/ 3826577.5, 461021.3125, 5064893.0 /)
    real(dp) :: x
    integer :: i,j
    Do i=1,3
        x=0.
        Do j=1,3
            x=x + RotMat(i,j) * (ITRF(j)-center_CS002(j))
        enddo
        LOFAR(i)=x  ! 1=North, 2=East, 3=vertical(plumbline)
    enddo
    return
End Subroutine ITRF2LOFAR
!================================================
Subroutine RelDist(Source,LFRAnt,RDist)
!   RDist is distance relative to center of CS002, in units of time-samples
!  RDist= dt_(core-source) - dt_(ant-source)=dt_{cs} - dt_{as}
! t_core=t_s + dt_cs ; t_ant=t_s + dt_as  ; t_s=t_ant-dt_as
! thus: t_core=t_ant + Rdist
    use constants, only : dp,sample,c_mps
    implicit none
    real(dp), intent(in) :: Source(3)
    real(dp), intent(in) :: LFRAnt(3)
    real(dp), intent(out) :: RDist
    Real(dp) :: D
   Real(dp) :: IndxRefrac
   Real(dp), external :: RefracIndex
    !write(2,*) 'source',source
    !write(2,*) 'LFRAnt',LFRAnt
    D=sum((Source(:)-LFRAnt(:))*(Source(:)-LFRAnt(:)))
    RDist = sqrt(D)
    D=sum(Source(:)*Source(:))  ! core=center station CS002 is at zero
    RDist = RDist - sqrt(D)         ! units of [meter]
    RDist = Rdist*RefracIndex( Source(3) )/(c_mps*sample) ! Convert to units of [samples]
    Return
End Subroutine RelDist
!================================================
Real(kind=8) Function tShift_ms(Source)  !
    use constants, only : dp,sample,c_mps
    implicit none
    real(dp), intent(in) :: Source(1:3)
    Real(dp), external :: RefracIndex
    !real(dp), intent(out) :: tShift
    tShift_ms = RefracIndex(Source(3))*sqrt(Source(1)*Source(1)+Source(2)*Source(2)+Source(3)*Source(3))*1000.d0/c_mps
    Return
End Function tShift_ms
!================================================    Real(dp), external ::
Real(kind=8) Function tShift_smpl(Source)  !
    use constants, only : dp,sample,c_mps
    implicit none
    real(dp), intent(in) :: Source(1:3)
    Real(dp), external :: RefracIndex
    !real(dp), intent(out) :: tShift
    tShift_smpl = RefracIndex(Source(3))*sqrt(Source(1)*Source(1)+Source(2)*Source(2)+Source(3)*Source(3))/(c_mps*sample)
    Return
End Function tShift_smpl
!================================================
Real(kind=8) Function RefracIndex(h)
!        Parameters of atmospheric density as used in Corsika
    use constants, only : dp,HeightCorrectIndxRef
    implicit none
    integer, parameter :: AtmHei_Dim=1000
    Real(dp), parameter :: AtmHei_step=10.d0 ! [m]
    real(dp), intent(in) :: h
    Real(dp), save :: Xi(0:AtmHei_Dim)
    Real(dp) :: hi
    Logical, save :: First=.true.
    integer :: i
    !
    If (First) Then ! Calculate atmosphere (from MGMR3D)
      First=.false.
      Call AveIndxRefr(AtmHei_dim,AtmHei_step, xi)
    EndIf
    If(HeightCorrectIndxRef) Then
       hi=h/AtmHei_step
       i=INT(hi)
       If(i.ge.AtmHei_Dim) then
         RefracIndex=xi(AtmHei_Dim)
       ElseIf(i.gt.0) then
         !RefracIndex=xi(i)*(hi-i) + xi(i+1)*(i+1.-hi)
         RefracIndex=xi(i)*(i+1.-hi) + xi(i+1)*(hi-i) ! Corrected 7 oct 2023
       Else
         RefracIndex=Xi(0)
       EndIf
    Else
       RefracIndex=Xi(0)
    EndIf
    !Write(2,*) 'Function RefracIndex:',h,RefracIndex
End Function RefracIndex
!=========================================
Subroutine AveIndxRefr(AtmHei_dim,AtmHei_step, xi)
   use constants, only : dp
   implicit none
   integer, intent(in) :: AtmHei_Dim
   Real(dp), intent(in) :: AtmHei_step ! [m]
   Real(dp), intent(out) :: Xi(0:AtmHei_Dim)
   !REAL(dp), parameter :: Refrac=1.000267394d0 ! Index of refraction at sea level (20^0, 100% hum) From https://emtoolbox.nist.gov/Wavelength/Documentation.asp
   REAL(dp), parameter :: Refrac= 1.0003     ! Index of refraction at sea level
   real(dp) :: height
   real(dp):: mass,b,c,RefOrho
   integer :: i
   !
   xi(0)=Refrac
   mass=0.
   RefOrho=(xi(0)-1.d0)*9941.8638d0/1222.6562d0
   do i = 1, AtmHei_dim
      height=i*AtmHei_step  ! distance to ground along shower axis
       if (height.ge.10d3) then
            b = 1305.5948; c = 6361.4304
       elseif (height.ge.4d3) then
            b = 1144.9069d0; c = 8781.5355d0
       else
            b = 1222.6562d0; c = 9941.8638d0
       endif
      mass = mass + (b/c) * exp(-height/c) * AtmHei_step ! assume density is about constant
      !write(2,*) height, a, PenDepth(0)-X_rh-RPenDepth(i)
      !
      ! calculate  averaged refractivity per meter
      xi(i) = 1.d0+ RefOrho* mass/height   ! mean index of refraction for rays from height to 0
    end do
    !write(2,*) 'Xi:',Xi(0),Xi(1),Xi(10),Xi(100),Xi(1000)
    Return
End Subroutine AveIndxRefr
!=========================================
Real(kind=8) Function SubRelDist(SrcPos,i_ant,i_chunk)
   use Chunk_AntInfo, only : Ant_pos, Ant_RawSourceDist
   use constants, only : dp
   Implicit none
   Real(dp), intent(in) ::  SrcPos(3)
   Integer, intent(in) :: i_ant,i_chunk
   Real(dp) :: RDist
      Call RelDist(SrcPos(1),Ant_pos(1,i_ant,i_chunk),RDist)  ! source position may have changed compared to previous
      SubRelDist=Rdist - Ant_RawSourceDist(i_ant,i_chunk)
   !write(2,*) 'SubRelDist:',Rdist, Ant_RawSourceDist(i_ant,i_chunk)
End Function SubRelDist
!---------------------------------
Pure Subroutine SetSmooth(N_Smth, ParabolicSmooth, Smooth)
! Subroutine SetSmooth(N_Smth, ParabolicSmooth, Smooth)
   ! smooth(j)=weight factors for power (amplitude^2) per sample, used in calculating voxel power as in:
   ! Stk(m,n)= Stk(m,n) + smooth(j)*CMTime_pix(IntfLead+i_s+j,m)*Conjg(CMTime_pix(IntfLead+i_s+j,n))
   use constants, only : dp
   Implicit none
   Integer, intent(in) :: N_Smth
   logical, intent(in) :: ParabolicSmooth
   Real(dp), intent(out) :: Smooth(-N_smth:N_smth)
   integer :: i
   !Allocate( Smooth(-N_smth:N_smth) )
   Smooth(:)=0.
   Smooth(0)=1.
   !write(2,*) 'SetSmooth:',N_smth, ParabolicSmooth
   Smooth(N_smth/2)=0.5  ! needed for block profile for even values of N_smth
   Do i=1,N_smth
      If(ParabolicSmooth) Then
         Smooth(i)=(1.-(Real(i)/N_smth)**2)**2.5  ! Parabola to power to make it =0.5 at N_smth/2
      Else
         If(i.lt.(N_smth+1)/2) Smooth(i)=1.      ! Block gives indistinguishable results from parabola; parabola has a tiny bit smaller volume around peak.
         Smooth(-i)=Smooth(i)
      EndIf
   EndDo
   !SmPow=SUM(Smooth(:))
   Smooth(:)=Smooth(:)/SUM(Smooth(:))
   Return
End Subroutine SetSmooth
!-----------------------------------
Function random_stdnormal() Result(x)
!  https://masuday.github.io/fortran_tutorial/random.html
! General interest: https://en.wikibooks.org/wiki/Fortran/Fortran_procedures_and_functions
   implicit none
   real :: x
   real,parameter :: pi=3.14159265
   real :: u1,u2
   ! call random_number(r) gives 0=< r < 1 i.e. including 0, excluding 1
   call random_number(u1)
   call random_number(u2)
   x = sqrt(-2*log(1-u1))*cos(2*pi*u2)
   Return
end Function random_stdnormal
!-----------------------------------
Subroutine random_stdnormal3D(x)
!  https://masuday.github.io/fortran_tutorial/random.html
! General interest: https://en.wikibooks.org/wiki/Fortran/Fortran_procedures_and_functions
   implicit none
   real, intent(out) :: x(1:3)
   real :: R,D, Epsln=epsilon(Epsln)
   real :: u1,u2,u3
   Real :: random_stdnormal
   ! call random_number(r) gives 0=< r < 1 i.e. including 0, excluding 1
   call random_number(u1)  ! may include zero
   call random_number(u2)
   call random_number(u3)
   R=sqrt((0.5-u1)**2+(0.5-u2)**2+(0.5-u3)**2+Epsln) ! to prevent zero
   D=random_stdnormal()! can be zero
   D=((abs(d))**(1/3.))  ! to have distances distributed like [d^2 x gaussian(d)]
   !D=sqrt(abs(d)) * Space_width  ! to have distances distributed like [d x gaussian(d)]
   x(1)= D*(0.5-u1)/R
   x(2)= D*(0.5-u2)/R
   x(3)= D*(0.5-u3)/R
   Return
end Subroutine random_stdnormal3D
