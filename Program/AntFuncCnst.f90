Module AntFunCconst
    !   From Module constants:
    !  integer, parameter:: dp=kind(0.d0)
    !  COMPLEX(dp), parameter :: CI=CMPLX(0._dp,1._dp)
    !  REAL(dp), parameter :: pi=2.d0*asin(1.d0)
   use constants, only : dp
   Integer, parameter :: Freq_max=90
   Integer, parameter :: Freq_min=20
   Integer, parameter :: j_dim=4 ! Number of odd Legendre polynominals that are used
   Integer, parameter :: m_dim=2 ! Number of odd Chebyshev polynominals that are used
   Real(dp), save :: Gain(Freq_min:Freq_max)
   Complex(dp), save :: Ant_p(Freq_min:Freq_max,1:j_dim,1:m_dim)
   Complex(dp), save :: Ant_t(Freq_min:Freq_max,1:j_dim,1:m_dim)
   Complex(dp) :: J_0p(Freq_min:Freq_max),J_0t(Freq_min:Freq_max),J_1p(Freq_min:Freq_max),J_1t(Freq_min:Freq_max)
   Complex(dp) :: Ji_p0(Freq_min:Freq_max),Ji_t0(Freq_min:Freq_max),Ji_p1(Freq_min:Freq_max),Ji_t1(Freq_min:Freq_max)
   Logical, save :: AntTst=.false.
End Module AntFunCconst
!
