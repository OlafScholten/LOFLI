subroutine spline (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function ispline can be used for interpolation
!======================================================================
implicit none
integer, intent(in) :: n
double precision, intent(in) :: x(n), y(n)
double precision, intent(out) :: b(n), c(n), d(n)
integer i, j, gap
double precision h

    gap = n-1
    ! check input
    if ( n < 2 ) return
    if ( n < 3 ) then   ! n=2
      b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.
      return
    end if
    !
    ! step 1: preparation
    !
    d(1) = x(2) - x(1)
    c(2) = (y(2) - y(1))/d(1)
    do i = 2, gap
      d(i) = x(i+1) - x(i)
      b(i) = 2.0*(d(i-1) + d(i))
      c(i+1) = (y(i+1) - y(i))/d(i)
      c(i) = c(i+1) - c(i)
    end do
    !
    ! step 2: end conditions
    !
    b(1) = -d(1)
    b(n) = -d(n-1)
    c(1) = 0.0
    c(n) = 0.0
    if(n /= 3) then
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
    end if
    !
    ! step 3: forward elimination
    !
    do i = 2, n
      h = d(i-1)/b(i-1)
      b(i) = b(i) - h*d(i-1)
      c(i) = c(i) - h*c(i-1)
    end do
    !
    ! step 4: back substitution
    !
    c(n) = c(n)/b(n)
    do j = 1, gap
      i = n-j
      c(i) = (c(i) - d(i)*c(i+1))/b(i)
    end do
    !
    ! step 5: compute spline coefficients
    !
    b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
    do i = 1, gap
      b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
      d(i) = (c(i+1) - c(i))/d(i)
      c(i) = 3.*c(i)
    end do
    c(n) = 3.0*c(n)
    d(n) = d(n-1)
end subroutine spline

  function ispline(u, x, y, b, c, d, n)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
implicit none
double precision ispline
integer, intent(in) :: n
double precision, intent(in) :: u, x(n), y(n), b(n), c(n), d(n)
integer i, j, k
double precision dx

    ! if u is ouside the x() interval take a boundary value (left or right)
    if(u <= x(1)) then
      ispline = y(1)
      return
    end if
    if(u >= x(n)) then
      ispline = y(n)
      return
    end if

    !*
    !  binary search for for i, such that x(i) <= u <= x(i+1)
    !*
    i = 1
    j = n+1
    do while (j > i+1)
      k = (i+j)/2
      if(u < x(k)) then
        j=k
        else
        i=k
       end if
    end do
    !*
    !  evaluate spline interpolation
    !*
    dx = u - x(i)
    ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline
! ============================================
subroutine spline_cubic_set( n, t, y, ypp )
!subroutine spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )

!*****************************************************************************80
!
!! SPLINE_CUBIC_SET computes the second derivatives of a piecewise cubic spline.
!
!  Discussion:
!
!    For data interpolation, the user must call SPLINE_CUBIC_SET to
!    determine the second derivative data, passing in the data to be
!    interpolated, and the desired boundary conditions.
!
!    The data to be interpolated, plus the SPLINE_CUBIC_SET output,
!    defines the spline.  The user may then call SPLINE_CUBIC_VAL to
!    evaluate the spline at any point.
!
!    The cubic spline is a piecewise cubic polynomial.  The intervals
!    are determined by the "knots" or abscissas of the data to be
!    interpolated.  The cubic spline has continous first and second
!    derivatives over the entire interval of interpolation.
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A(IVAL)
!             + B(IVAL) * ( T - T(IVAL) )
!             + C(IVAL) * ( T - T(IVAL) )^2
!             + D(IVAL) * ( T - T(IVAL) )^3
!
!    If we assume that we know the values Y(*) and YPP(*), which represent
!    the values and second derivatives of the spline at each knot, then
!    the coefficients can be computed as:
!
!      A(IVAL) = Y(IVAL)
!      B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C(IVAL) = YPP(IVAL) / 2
!      D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!    Since the first derivative of the spline is
!
!      SPL'(T) =     B(IVAL)
!              + 2 * C(IVAL) * ( T - T(IVAL) )
!              + 3 * D(IVAL) * ( T - T(IVAL) )^2,
!
!    the requirement that the first derivative be continuous at interior
!    knot I results in a total of N-2 equations, of the form:
!
!      B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1))
!      + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))^2 = B(IVAL)
!
!    or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
!
!      ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!      - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
!      + YPP(IVAL-1) * H(IVAL-1)
!      + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
!      =
!      ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
!
!    or
!
!      YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
!      + YPP(IVAL) * H(IVAL)
!      =
!      6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!
!    Boundary conditions must be applied at the first and last knots.
!    The resulting tridiagonal system can be solved for the YPP values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 June 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points; N must be
!    at least 2.
!
!    Input, real ( kind = 8 ) T(N), the points where data is specified.
!    The values should be distinct, and increasing.
!
!    Input, real ( kind = 8 ) Y(N), the data values to be interpolated.
!
!    Input, integer ( kind = 4 ) IBCBEG, the left boundary condition flag:
!    0: the spline should be a quadratic over the first interval;
!    1: the first derivative at the left endpoint should be YBCBEG;
!    2: the second derivative at the left endpoint should be YBCBEG;
!    3: Not-a-knot: the third derivative is continuous at T(2).
!
!    Input, real ( kind = 8 ) YBCBEG, the left boundary value, if needed.
!
!    Input, integer ( kind = 4 ) IBCEND, the right boundary condition flag:
!    0: the spline should be a quadratic over the last interval;
!    1: the first derivative at the right endpoint should be YBCEND;
!    2: the second derivative at the right endpoint should be YBCEND;
!    3: Not-a-knot: the third derivative is continuous at T(N-1).
!
!    Input, real ( kind = 8 ) YBCEND, the right boundary value, if needed.
!
!    Output, real ( kind = 8 ) YPP(N), the second derivatives of
!    the cubic spline.
!
  implicit none

  integer, intent(in) :: n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  real ( kind = 8 ) a3(n)
  real ( kind = 8 ) a4(n)
  real ( kind = 8 ) a5(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ),parameter :: ibcbeg=3
  integer ( kind = 4 ),parameter :: ibcend=3
  integer ( kind = 4 ) info
  real ( kind = 8 ), intent(in) :: t(n)
  real ( kind = 8 ), intent(in) :: y(n)
  real ( kind = 8 ),parameter :: ybcbeg=0
  real ( kind = 8 ),parameter :: ybcend=0
  real ( kind = 8 ), intent(out) :: ypp(n)
!
!  Check.
!
  if ( n <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The number of knots must be at least 2.'
    write ( *, '(a,i8)' ) '  The input value of N = ', n
    stop !1
  end if

  do i = 1, n - 1
    if ( t(i+1) <= t(i) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
      write ( *, '(a)' ) '  The knots must be strictly increasing, but'
      write ( *, '(a,i8,a,g14.6)' ) '  T(',  i,') = ', t(i)
      write ( *, '(a,i8,a,g14.6)' ) '  T(',i+1,') = ', t(i+1)
      stop !1
    end if
  end do
!
!  Zero out the matrix.
!
  a1(1:n) = 0.0D+00
  a2(1:n) = 0.0D+00
  a3(1:n) = 0.0D+00
  a4(1:n) = 0.0D+00
  a5(1:n) = 0.0D+00
!
!  Set the first equation.
!
  if ( ibcbeg == 0 ) then
    b(1) = 0.0D+00
    a3(1) =  1.0D+00
    a4(1) = -1.0D+00
  else if ( ibcbeg == 1 ) then
    b(1) = ( y(2) - y(1) ) / ( t(2) - t(1) ) - ybcbeg
    a3(1) = ( t(2) - t(1) ) / 3.0D+00
    a4(1) = ( t(2) - t(1) ) / 6.0D+00
  else if ( ibcbeg == 2 ) then
    b(1) = ybcbeg
    a3(1) = 1.0D+00
    a4(1) = 0.0D+00
  else if ( ibcbeg == 3 ) then
    b(1) = 0.0D+00
    a3(1) = - ( t(3) - t(2) )
    a4(1) =   ( t(3)        - t(1) )
    a5(1) = - (        t(2) - t(1) )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The boundary flag IBCBEG must be 0, 1, 2 or 3.'
    write ( *, '(a,i8)' ) '  The input value is IBCBEG = ', ibcbeg
    stop !1
  end if
!
!  Set the intermediate equations.
!
  do i = 2, n - 1
    b(i) = ( y(i+1) - y(i) ) / ( t(i+1) - t(i) ) &
         - ( y(i) - y(i-1) ) / ( t(i) - t(i-1) )
    a2(i) = ( t(i+1) - t(i)   ) / 6.0D+00
    a3(i) = ( t(i+1) - t(i-1) ) / 3.0D+00
    a4(i) = ( t(i)   - t(i-1) ) / 6.0D+00
  end do
!
!  Set the last equation.
!
  if ( ibcend == 0 ) then
    b(n) = 0.0D+00
    a2(n) = -1.0D+00
    a3(n) = 1.0D+00
  else if ( ibcend == 1 ) then
    b(n) = ybcend - ( y(n) - y(n-1) ) / ( t(n) - t(n-1) )
    a2(n) = ( t(n) - t(n-1) ) / 6.0D+00
    a3(n) = ( t(n) - t(n-1) ) / 3.0D+00
  else if ( ibcend == 2 ) then
    b(n) = ybcend
    a2(n) = 0.0D+00
    a3(n) = 1.0D+00
  else if ( ibcend == 3 ) then
    b(n) = 0.0D+00
    a1(n) = - ( t(n) - t(n-1) )
    a2(n) =   ( t(n)          - t(n-2) )
    a3(n) = - (        t(n-1) - t(n-2) )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The boundary flag IBCEND must be 0, 1, 2 or 3.'
    write ( *, '(a,i8)' ) '  The input value is IBCEND = ', ibcend
    stop !1
  end if
!
!  Special case:
!    N = 2, IBCBEG = IBCEND = 0.
!
  if ( n == 2 .and. ibcbeg == 0 .and. ibcend == 0 ) then

    ypp(1) = 0.0D+00
    ypp(2) = 0.0D+00
!
!  Solve the linear system.
!
  else

    call penta ( n, a1, a2, a3, a4, a5, b, ypp )

  end if

  return
end subroutine spline_cubic_set
subroutine penta ( n, a1, a2, a3, a4, a5, b, x )

!*****************************************************************************80
!
!! PENTA solves a pentadiagonal system of linear equations.
!
!  Discussion:
!
!    The matrix A is pentadiagonal.  It is entirely zero, except for
!    the main diagaonal, and the two immediate sub- and super-diagonals.
!
!    The entries of Row I are stored as:
!
!      A(I,I-2) -> A1(I)
!      A(I,I-1) -> A2(I)
!      A(I,I)   -> A3(I)
!      A(I,I+1) -> A4(I)
!      A(I,I-2) -> A5(I)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 June 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Cheney, Kincaid,
!    Numerical Mathematics and Computing,
!    1985, pages 233-236.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), A3(N), A4(N), A5(N), the nonzero
!    elements of the matrix.  Note that the data in A2, A3 and A4
!    is overwritten by this routine during the solution process.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
  implicit none

  integer, intent(in) :: n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  real ( kind = 8 ) a3(n)
  real ( kind = 8 ) a4(n)
  real ( kind = 8 ) a5(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmult

  do i = 2, n - 1
    xmult = a2(i) / a3(i-1)
    a3(i) = a3(i) - xmult * a4(i-1)
    a4(i) = a4(i) - xmult * a5(i-1)
    b(i) = b(i) - xmult * b(i-1)
    xmult = a1(i+1) / a3(i-1)
    a2(i+1) = a2(i+1) - xmult * a4(i-1)
    a3(i+1) = a3(i+1) - xmult * a5(i-1)
    b(i+1) = b(i+1) - xmult * b(i-1)
  end do

  xmult = a2(n) / a3(n-1)
  a3(n) = a3(n) - xmult * a4(n-1)
  x(n) = ( b(n) - xmult * b(n-1) ) / a3(n)
  x(n-1) = ( b(n-1) - a4(n-1) * x(n) ) / a3(n-1)
  do i = n - 2, 1, -1
    x(i) = ( b(i) - a4(i) * x(i+1) - a5(i) * x(i+2) ) / a3(i)
  end do

  return
end subroutine penta
subroutine spline_cubic_val ( n, t, y, ypp, tval, yval) !, ypval, yppval )

!*****************************************************************************80
!
!! SPLINE_CUBIC_VAL evaluates a piecewise cubic spline at a point.
!
!  Discussion:
!
!    SPLINE_CUBIC_SET must have already been called to define the
!    values of YPP.
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A
!             + B * ( T - T(IVAL) )
!             + C * ( T - T(IVAL) )^2
!             + D * ( T - T(IVAL) )^3
!
!    Here:
!      A = Y(IVAL)
!      B = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C = YPP(IVAL) / 2
!      D = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) T(N), the knot values.
!
!    Input, real ( kind = 8 ) Y(N), the data values at the knots.
!
!    Input, real ( kind = 8 ) YPP(N), the second derivatives of the
!    spline at the knots.
!
!    Input, real ( kind = 8 ) TVAL, a point, typically between T(1) and
!    T(N), at which the spline is to be evalulated.  If TVAL lies outside
!    this range, extrapolation is used.
!
!    Output, real ( kind = 8 ) YVAL, YPVAL, YPPVAL, the value of the spline, and
!    its first two derivatives at TVAL.
!
  implicit none

  integer, intent(in) :: n

  real ( kind = 8 ) dt
  real ( kind = 8 ) h
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real ( kind = 8 ), intent(in) :: t(n)
  real ( kind = 8 ), intent(in) :: tval
  real ( kind = 8 ), intent(in) :: y(n)
  real ( kind = 8 ), intent(in) :: ypp(n)
!  real ( kind = 8 ) yppval
!  real ( kind = 8 ) ypval
  real ( kind = 8 ), intent(out) :: yval
!
  !
    If(tval.lt.t(1) ) then
        yval=y(1)*exp(tval-t(1))
        return
    elseif(tval.gt. t(n)) then
        yval=y(n)*exp(t(n)-tval)
        return
    else
        !  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
        !  Values below T(1) or above T(N) use extrapolation.
        !
          call r8vec_bracket ( n, t, tval, left, right )
        !
        !  Evaluate the polynomial.
        !
          dt = tval - t(left)
          h = t(right) - t(left)
          !
          yval = y(left) &
               + dt * ( ( y(right) - y(left) ) / h &
                      - ( ypp(right) / 6.0D+00 + ypp(left) / 3.0D+00 ) * h &
               + dt * ( 0.5D+00 * ypp(left) &
               + dt * ( ( ypp(right) - ypp(left) ) / ( 6.0D+00 * h ) ) ) )

        !  ypval = ( y(right) - y(left) ) / h &
        !       - ( ypp(right) / 6.0D+00 + ypp(left) / 3.0D+00 ) * h &
        !       + dt * ( ypp(left) &
        !       + dt * ( 0.5D+00 * ( ypp(right) - ypp(left) ) / h ) )
        !
        !  yppval = ypp(left) + dt * ( ypp(right) - ypp(left) ) / h
    endif
  return
end subroutine spline_cubic_val

subroutine r8vec_bracket ( n, x, xval, left, right )

!*****************************************************************************80
!
!! R8VEC_BRACKET searches a sorted R8VEC for successive brackets of a value.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input, real ( kind = 8 ) X(N), an array sorted into ascending order.
!
!    Input, real ( kind = 8 ) XVAL, a value to be bracketed.
!
!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xval

  do i = 2, n - 1

    if ( xval < x(i) ) then
      left = i - 1
      right = i
      return
    end if

   end do

  left = n - 1
  right = n

  return
end subroutine r8vec_bracket
subroutine Spline_Average2_Shift ( n, t, y00, ypp00, y01, ypp01, y10, ypp10, y11, ypp11, w1b, w2b, tshft, m, tas, yas )
!
!          DEx =(1.-did)*((1.-CD_d)*Ex_nub(inu,id  ,CD_i  ) + CD_d*Ex_nub(inu,id  ,CD_i+1)) &
!                 +  did*((1.-CD_d)*Ex_nub(inu,id+1,CD_i  ) + CD_d*Ex_nub(inu,id+1,CD_i+1))
!   yas(tas(i))=(1-wb)*ya(tas(i)-tshft) + wb*yb(tas(i)-tshft) where tshft>0
!   given ya(t(i)) and yb(t(i))
!   where in general the grids [t(i)] and [tas(i)-tshft] are displaced a non-integer value
!       which necessitates the use of spline interpolation
!
  implicit none
  integer, intent(in) :: n,m
  real ( kind = 8 ), intent(in) :: t(n),y00(n),ypp00(n),y01(n),ypp01(n),y10(n),ypp10(n),y11(n),ypp11(n),w1b,w2b,tshft,tas(m)
  real ( kind = 8 ), intent(out) :: yas(m)
  real ( kind = 8 ) :: w1a,w2a,wy(n),wypp(n)
!
    w1a=1.d0-w1b
    w2a=1.d0-w2b
    wy(:)  =w1a*(w2a*y00(:)    + w2b*y01(:)  ) + w1b*(w2a*y10(:)    + w2b*y11(:)  )
    wypp(:)=w1a*(w2a*ypp00(:)  + w2b*ypp01(:)) + w1b*(w2a*ypp10(:)  + w2b*ypp11(:))
    call Spline_Shift ( n, t, wy, wypp, tshft, m, tas, yas )
    return
end subroutine Spline_Average2_Shift
subroutine Spline_Average_Shift ( n, t, ya, yppa, yb, yppb, wb, tshft, m, tas, yas )
!
!   Sample call:
!        call  Spline_Average_Shift(tTrace_dim_b, t_tb(1), &
!            Ex_tb(1,id,CD_i), Ex_spld(1,id,CD_i), Ex_tb(1,id,CD_i+1), Ex_spld(1,id,CD_i+1), &
!            CD_d, tshft, tTrace_dim_o, t_to(1), DEx )
!
!   yas(tas(i))=(1-wb)*ya(tas(i)-tshft) + wb*yb(tas(i)-tshft) where tshft>0
!   given ya(t(i)) and yb(t(i))
!   where in general the grids [t(i)] and [tas(i)-tshft] are displaced a non-integer value
!       which necessitates the use of spline interpolation
!
  implicit none
  integer, intent(in) :: n,m
  real ( kind = 8 ), intent(in) :: t(n),ya(n),yppa(n),yb(n),yppb(n),wb,tshft,tas(m)
  real ( kind = 8 ), intent(out) :: yas(m)
  real ( kind = 8 ) :: wa,wy(n),wypp(n)
!
    wa=1.d0-wb
    wy(:)=wa*ya(:)    + wb*yb(:)
    wypp(:)=wa*yppa(:)    + wb*yppb(:)
    call Spline_Shift ( n, t, wy, wypp, tshft, m, tas, yas )
    return
end subroutine Spline_Average_Shift
subroutine Spline_Shift ( n, t, y, ypp, tshft, m, ts, ys )
!
! To set-up spline polynominals :
!           call spline_cubic_set( tTrace_dim_b, t_tb(1), Ex_tb(1,idi,CD_i), Ex_spld(1,idi,CD_i) )
!
!   yas(ts(i))=y( ts(i)-tshft )     where tshft>0
!   given y(t(i))
!   where in general the grids [t(i)] and [ts(i)-tshft] are displaced a non-integer value
!       which necessitates the use of spline interpolation
!
  implicit none
  integer, intent(in) :: n,m
  real ( kind = 8 ), intent(in) :: t(n),y(n),ypp(n),tshft,ts(m)
  real ( kind = 8 ), intent(out) :: ys(m)
  real ( kind = 8 ) :: yl,yr,yppl,yppr
  real ( kind = 8 ) :: dt,tval
  real ( kind = 8 ) h
  integer ( kind = 4 ) right,i,ii
!
    right=2
    Do ii=1,m-1
        if (t(1).le.ts(ii)-tshft) exit
        ys(ii)=0.d0
    enddo
    Do i=ii,m
        tval=ts(i)-tshft
        Do while (t(right).le.tval)
            right=right+1
            if(right.gt.n) goto 9
        enddo
!
!  Evaluate the spline polynomial.
!
        dt = tval - t(right-1)
        h = t(right) - t(right-1)
        yl  =y(right-1)
        yr  =y(right)
        yppl=ypp(right-1)
        yppr=ypp(right)
        !
        ys(i) = yl &
           + dt * ( ( yr-yl ) / h - ( yppr / 6.0D+00 + yppl / 3.0D+00 ) * h &
                    + dt * ( 0.5D+00 * yppl &
                             + dt * ( ( yppr - yppl ) / ( 6.0D+00 * h ) ) ) )
    enddo
  return
9   continue
    Do ii=i,m
        ys(ii)=0.d0
    enddo
    return
end subroutine Spline_Shift
