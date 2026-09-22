module former_mat_util


  contains 

  subroutine matscl(a, s, b)
    ! Multiplies a scalar, s, to a 3-by-3 matrices a.
  
    use module_kind_variables, only: kind_double
    implicit none
  
    real(kind_double), intent(in)   :: a(3, 3), s
    real(kind_double), intent(out)  :: b(3, 3)
  
    integer  :: j, i
  
    do j = 1, 3
      do i = 1, 3
        b(i, j) = a(i, j)*s
      end do
    end do
  end subroutine matscl


  subroutine matmul(a, b, c)
    ! Multiplies 3-by-3 matrices a & b, and stores the result in c.
  
    USE module_kind_variables, ONLY: kind_double
    implicit none
  
    real(kind_double), intent(in)   :: a(3, 3), b(3, 3)
    real(kind_double), intent(out)  :: c(3, 3)
  
    integer  :: j, i, k
  
    do j = 1, 3
      do i = 1, 3
        c(i, j) = 0d0
        do k = 1, 3
          c(i, j) = c(i, j) + a(i, k)*b(k, j)
        end do
      end do
    end do
  end subroutine matmul
  
  
  subroutine mattrp(a, at)
    ! Transposes a 3-by-3 matrix, a, and stores the result in at.
  
    USE module_kind_variables, ONLY: kind_double
  
    implicit none
  
    real(kind_double), intent(in)   :: a(3, 3)
    real(kind_double), intent(out)  :: at(3, 3)
  
    integer  :: j, i
  
    do j = 1, 3
      do i = 1, 3
        at(i, j) = a(j, i)
      end do
    end do
  end subroutine mattrp
  
  subroutine matinv(a, ai)
    ! Invert a 3-by-3 matrix a, and store the result in ai
  
    USE module_kind_variables, ONLY: kind_double
    implicit none
  
    real(kind_double), intent(in)   :: a(3, 3)
    real(kind_double), intent(out)  :: ai(3, 3)
  
    integer  :: j, jm, jp, i, im, ip
    real(kind_double)   :: deta, detai
  
  
    ! Transpose cofactor matrix
    do j = 1, 3
      jm = mod(j + 1, 3) + 1
      jp = mod(j, 3) + 1
      do i = 1, 3
        im = mod(i + 1, 3) + 1
        ip = mod(i, 3) + 1
        ai(i, j) = a(jp, ip)*a(jm, im) - a(jp, im)*a(jm, ip)                      ! transposed cofactor matrix
      end do
    end do
  
    ! Determinant
    deta = 0.d0
    do i = 1, 3
      deta = deta + a(i, 1)*ai(1, i)
    end do
  
    detai = 1d0/deta
  
    ! Inverse matrix
    do j = 1, 3
      do i = 1, 3
        ai(i, j) = ai(i, j)*detai
      end do
    end do
  end subroutine matinv
  
  
  subroutine matcof(a, c)
    ! Given a 3-by-3 matrix a, calculates its cofactor matrix and stores the result in c.
  
    USE module_kind_variables, ONLY: kind_double
  
    implicit none
  
    real(kind_double), intent(in)   :: a(3, 3)
    real(kind_double), intent(out)  :: c(3, 3)
  
    integer  :: j, jm, jp, i, im, ip
  
    do j = 1, 3
      jm = mod(j + 1, 3) + 1
      jp = mod(j, 3) + 1
      do i = 1, 3
        im = mod(i + 1, 3) + 1
        ip = mod(i, 3) + 1
        c(i, j) = a(ip, jp)*a(im, jm) - a(im, jp)*a(ip, jm)
      end do
    end do
  end subroutine matcof
  
  
  real(kind(0.0d0)) function detmat(a)
    ! Returns the determinant of a 3-by-3 matrices a.
  
    USE module_kind_variables,  ONLY: kind_double
  
    implicit none
  
    real(kind_double), intent(in)   :: a(3, 3)
  
    detmat = a(1, 1)*(a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)) &
             - a(2, 1)*(a(1, 2)*a(3, 3) - a(1, 3)*a(3, 2)) &
             + a(3, 1)*(a(1, 2)*a(2, 3) - a(1, 3)*a(2, 2))
  
  end function detmat
  
  
  subroutine boxmat()
    ! ets up matrices associated with the MD box.
  
    USE module_kind_variables, ONLY: kind_double
    use ondm_gen_com_m
    implicit none
  
    real(kind_double), dimension(3, 3)    :: atit, sgm, ah
  
    ! Areal tensor, SGM
    call mattrp(at, ah)
    call matcof(ah, sgm)
    ! MD-box volumeb
    volu = ah(1, 1)*sgm(1, 1) + ah(2, 1)*sgm(2, 1) + ah(3, 1)*sgm(3, 1)
    call matscl(sgm, 1d0/volu, atit)
    ! Inverse MD-box tensor, HI
    call mattrp(atit, ati)
  
  end subroutine boxmat
  
  
  !ACEC! subroutine eigen3(a, d, v)
  !ACEC!   ! Diagonalizes a real symmetric 3x3 matrix.
  !ACEC!   !     a(3,3): Input matrix
  !ACEC!   !     d(3):   Return eigenvalues
  !ACEC!   !     v(3,3): Return eigenvectors--v(*,i)  =  i-th eigenvector
  !ACEC! 
  !ACEC!   USE module_kind_variables, ONLY: kind_double
  !ACEC! 
  !ACEC!   implicit none
  !ACEC! 
  !ACEC!   real(kind_double)   :: a(3, 3)
  !ACEC!   real(kind_double)   :: v(3, 3), d(3)
  !ACEC!   integer  :: j, i
  !ACEC! 
  !ACEC! 
  !ACEC!   ! Find eigenvalues
  !ACEC!   call jacobi(a, 3, 3, d, v, i)
  !ACEC!   ! Restore the original matrix
  !ACEC!   do i = 1, 3
  !ACEC!     do j = i + 1, 3
  !ACEC!       a(i, j) = a(j, i)
  !ACEC!     end do
  !ACEC!   end do
  !ACEC! end subroutine eigen3
  
  subroutine jacobi(a, n, np, d, v, nrot)
    ! Diagonalizes a real symmetric matrix by Jacobi transformation. From "Numerical Recipes".
  
    USE module_kind_variables, ONLY: kind_double
  
    implicit none
  
    integer  :: n, np, nrot
    real(kind_double)   :: a(np, np)
    real(kind_double)   :: v(np, np), d(np)
  
    integer  :: j, i, iq, ip, nmax
    parameter(nmax=500)
    real(kind_double)   :: sm, tresh, g, h, t, theta, c, s, tau
    real(kind_double)   :: b(nmax), z(nmax)
  
    do ip = 1, n
      do iq = 1, n
        v(ip, iq) = 0.
      end do
      v(ip, ip) = 1.
    end do
    do ip = 1, n
      b(ip) = a(ip, ip)
      d(ip) = b(ip)
      z(ip) = 0.
    end do
    nrot = 0
    do i = 1, 50
      sm = 0.
      do ip = 1, n - 1
        do iq = ip + 1, n
          sm = sm + abs(a(ip, iq))
        end do
      end do
      if (sm == 0.) return
      if (i < 4) then
        tresh = 0.2*sm/n**2
      else
        tresh = 0.
      end if
      do ip = 1, n - 1
        do iq = ip + 1, n
          g = 100.*abs(a(ip, iq))
          if ((i > 4) .and. (abs(d(ip)) + &
                             g == abs(d(ip))) .and. (abs(d(iq)) + g == abs(d(iq)))) then
            a(ip, iq) = 0.
          else if (abs(a(ip, iq)) > tresh) then
            h = d(iq) - d(ip)
            if (abs(h) + g == abs(h)) then
              t = a(ip, iq)/h
            else
              theta = 0.5*h/a(ip, iq)
              t = 1./(abs(theta) + sqrt(1.+theta**2))
              if (theta < 0.) t = -t
            end if
            c = 1./sqrt(1 + t**2)
            s = t*c
            tau = s/(1.+c)
            h = t*a(ip, iq)
            z(ip) = z(ip) - h
            z(iq) = z(iq) + h
            d(ip) = d(ip) - h
            d(iq) = d(iq) + h
            a(ip, iq) = 0.
            do j = 1, ip - 1
              g = a(j, ip)
              h = a(j, iq)
              a(j, ip) = g - s*(h + g*tau)
              a(j, iq) = h + s*(g - h*tau)
            end do
            do j = ip + 1, iq - 1
              g = a(ip, j)
              h = a(j, iq)
              a(ip, j) = g - s*(h + g*tau)
              a(j, iq) = h + s*(g - h*tau)
            end do
            do j = iq + 1, n
              g = a(ip, j)
              h = a(iq, j)
              a(ip, j) = g - s*(h + g*tau)
              a(iq, j) = h + s*(g - h*tau)
            end do
            do j = 1, n
              g = v(j, ip)
              h = v(j, iq)
              v(j, ip) = g - s*(h + g*tau)
              v(j, iq) = h + s*(g - h*tau)
            end do
            nrot = nrot + 1
          end if
        end do
      end do
      do ip = 1, n
        b(ip) = b(ip) + z(ip)
        d(ip) = b(ip)
        z(ip) = 0.
      end do
    end do
    write (6, *) 'too many iterations in jacobi'
  
  end subroutine jacobi
  
  subroutine myrnd(rnd, dseed)
    ! Random-number generator.
  
    USE module_kind_variables, ONLY: kind_double
  
    implicit none
  
    real(kind_double)   :: dseed
    real(kind_double)   :: rnd
  
    real(kind_double)   :: d2p31m, d2p31
    data d2p31m/2147483647d0/
    data d2p31/2147483648d0/
  
  
    dseed = dmod(16807d0*dseed, d2p31m)
    rnd = dseed/d2p31
  end subroutine myrnd
  
  
  subroutine cspline(n, x, y, b, c, d)
    ! the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
    ! for a cubic interpolating spline
    ! s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
    ! for  x(i) .le. x .le. x(i+1)
    !
    !  input
    !    n = the number of data points or knots (n.ge.2)
    !    x = the abscissas of the knots in strictly increasing order
    !    y = the ordinates of the knots
    !  output
    !    b, c, d  = arrays of spline coefficients as defined above.
    !
    !  using  p  to denote differentiation,
    !    y(i) = s(x(i))
    !    b(i) = sp(x(i))
    !    c(i) = spp(x(i))/2
    !    d(i) = sppp(x(i))/6  (derivative from the right)
    !
    !  the accompanying function subprogram seval can be used
    !  to evaluate the spline.
  
    USE module_kind_variables, ONLY: kind_double
    implicit none
  
    integer, intent(in)  :: n
    real(kind_double), intent(in)   :: x(n)
    real(kind_double), intent(in)   :: y(n)
    real(kind_double), intent(inout)      :: b(n)
    real(kind_double), intent(inout)      :: c(n)
    real(kind_double), intent(inout)      :: d(n)
  
    integer  :: nm1
  
    ! write (6,*)'spline'
    nm1 = n - 1
    if (n < 2) return
    if (n >= 3) then
  
      ! set up tridiagonal system
      ! b = diagonal, d = offdiagonal, c = right hand side.
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      d(2:nm1) = x(3:nm1 + 1) - x(2:nm1)
      b(2:nm1) = 2.*(d(:nm1 - 1) + d(2:nm1))
      c(3:nm1 + 1) = (y(3:nm1 + 1) - y(2:nm1))/d(2:nm1)
      c(2:nm1) = c(3:nm1 + 1) - c(2:nm1)
      ! end conditions.  third derivatives at  x(1)  and  x(n)
      ! obtained from divided differences
      b(1) = -d(1)
      b(n) = -d(n - 1)
      c(1) = 0.
      c(n) = 0.
      if (n /= 3) then
        c(1) = c(3)/(x(4) - x(2)) - c(2)/(x(3) - x(1))
        c(n) = c(n - 1)/(x(n) - x(n - 2)) - c(n - 2)/(x(n - 1) - x(n - 3))
        c(1) = c(1)*d(1)**2/(x(4) - x(1))
        c(n) = -c(n)*d(n - 1)**2/(x(n) - x(n - 3))
      end if
  
      ! forward elimination
      !$€$!  do i = 2, n
      !$€$!    t = d(i - 1)/b(i - 1)
      !$€$!    b(i) = b(i) - t*d(i - 1)
      !$€$!    c(i) = c(i) - t*c(i - 1)
      !$€$!  end do
      !$€$!  
      !$€$!  ! back substitution
      !$€$!  c(n) = c(n)/b(n)
      !$€$!  do ib = 1, nm1
      !$€$!    i = n - ib
      !$€$!    c(i) = (c(i) - d(i)*c(i + 1))/b(i)
      !$€$!  end do
  
      ! c(i) is now the sigma(i) of the text
      ! compute polynomial coefficients
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
      b(:nm1) = (y(2:nm1 + 1) - y(:nm1))/d(:nm1) - d(:nm1)*(c(2:nm1 + 1) + 2.*c(:nm1))
      d(:nm1) = (c(2:nm1 + 1) - c(:nm1))/d(:nm1)
      c(:nm1) = 3.*c(:nm1)
      c(n) = 3.*c(n)
      d(n) = d(n - 1)
      return
    end if
  
    b(1) = (y(2) - y(1))/(x(2) - x(1))
    c(1) = 0.
    d(1) = 0.
    b(2) = b(1)
    c(2) = 0.
    d(2) = 0.
  
  end subroutine cspline
  
  
  SUBROUTINE DCSSMO(H, N, TNODE, G, WGS, RHO, GSMO, B, C, D)                ! DCS 10
    ! COMPUTES THE DISCRETE NATURAL CUBIC SPLINE DEFINED ON
    ! THE INTERVAL (TNODE(1),TNODE(N)) WHICH
    ! SMOOTHS THROUGH THE DATA (TNODE(I),G(I)),I=1,2,...,N. N MUST BE 2 OR GREATER.
  ! THE NODES MUST SATISFY TNODE(I).LT.TNODE(I+1).
  ! THE SOLUTION S(T) FOR T IN THE INTERVAL (TNODE(I),TNODE(I+1)) IS GIVEN BY
  !     S(T) = GSMO(I)+B(I)*(T-TNODE(I))+C(I)*(T-TNODE(I))**2+D(I)*(T-TNODE(I))**3
  !
  !  INPUT  PARAMETERS(NONE OF THE INPUT PARAMETERS ARE CHANGED
  !         BY THIS SUBROUTINE)
  !
  !  H     - THE STEP SIZE USED FOR THE DISCRETE CUBIC SPLINE
  !  N     - NUMBER OF NODES (TNODE) AND DATA VALUES(G)
  !  TNODE - REAL ARRAY CONTAINING THE NODES (TNODE(I).LT.
  !          TNODE(I+1)).
  !  G     - REAL ARRAY CONTAINING THE DATA VALUES.
  !  WGS   - REAL ARRAY CONTAINING THE WEIGHTS WGS(I)
  !          CORRESPONDING TO THE DATA (TNODE(I),G(I)).
  !  RHO   - SIMPLE REAL VARIABLE CONTAINING THE POSITIVE
  !          PARAMETER FOR VARYING THE SMOOTHNESS OF THE FIT.
  !          IF RHO IS SMALL SMOOTHNESS IS EMPHASIZED.
  !          IF RHO IS LARGE DATA FITTING IS EMPHASIZED.
  !
  !  OUTPUT PARAMETERS
  !
  !  GSMO  - REAL ARRAY CONTAINING THE SMOOTHED VALUES OF
  !          THE DATA G(I),I=1,2,....,N.
  !  B     - REAL ARRAY CONTAINING THE COEFFICIENTS B(I) FOR
  !          THE TERMS (T-TNODE(I)).
  !  C     - REAL ARRAY CONTAINING THE COEFFICIENTS C(I) FOR
  !          THE TERMS (T-TNODE(I))**2.
  !  D     - REAL ARRAY CONTAINING THE COEFFICIENTS D(I) FOR
  !          THE TERMS (T-TNODE(I))**3.
  use module_kind_variables, only: kind_double
  implicit none
  integer  :: N, N1, N2, N3, I, I1, I2, J, K, K1, K2, K3
  real(kind_double)      :: H, H2, H3, R6, HI3, HI4, HI5, &
                          HK1, HK2, ETA3, ETA4, BETA3, BETA4, BETA5, &
                          EPS3, EPS4, H2DHI, P
  real(kind_double)     :: TNODE(N), G(N), WGS(N), GSMO(N), B(N), C(N), D(N)
  real(kind_double)      :: RHO
  IF (N .EQ. 2) GO TO 180
  N1 = N - 1
  N2 = N1 - 1
  N3 = N2 - 1
  !  THE RIGHT HAND SIDE OF THE LINEAR SYSTEM FOR THE
  !  C(I)'S WILL NOW BE CONSTRUCTED.
  DO I = 1, N             ! 10
    C(I) = G(I)
  end do                  ! 10
  DO I = 1, N1            ! 20
    C(I) = (C(I + 1) - C(I))/(TNODE(I + 1) - TNODE(I))
  end do                  ! 20
  DO I = 1, N2            ! 30
    C(I) = 3.0*(C(I + 1) - C(I))
  end do                  ! 30
  !  THE RIGHT HAND SIDE IS NOW IN ARRAY C.
  !  THE P.D. 5 BANDED SYMMETRIC MATRIX WILL NOW BE CONSTRUCTED.
  !  THE THREE NEEDED DIAGONALS WILL BE STORED IN ARRAYS GSMO,B,D.
  H2 = H*H
  H3 = H2*H
  R6 = 6.0*H3/RHO
  HI3 = TNODE(2) - TNODE(1)
  HI4 = TNODE(3) - TNODE(2)
  ETA3 = HI3 + HI3 + H2/HI3
  BETA3 = R6/(WGS(1)*HI3)
  BETA4 = R6/(WGS(2)*HI3*HI4)
  EPS3 = (BETA3 + BETA4*HI4)/HI3
  H2DHI = H2/HI4
  ETA4 = HI4 + HI4 + H2DHI
  IF (N .EQ. 3) GO TO 60
  HI5 = TNODE(4) - TNODE(3)
  BETA5 = R6/(WGS(3)*HI4*HI5)
  EPS4 = (BETA4*HI3 + BETA5*HI5)/HI4
  GSMO(1) = ETA3 + ETA4 + BETA4 + BETA4 + EPS3 + EPS4
  P = H2DHI + BETA4 + BETA5 + EPS4
  B(1) = HI4 - P
  IF (N .EQ. 4) GO TO 50
  DO I = 2, N3            ! 40
    HI3 = HI4
    HI4 = HI5
    HI5 = TNODE(I + 3) - TNODE(I + 2)
    ETA3 = ETA4
    H2DHI = H2/HI4
    ETA4 = HI4 + HI4 + H2DHI
    BETA3 = BETA4
    BETA4 = BETA5
    BETA5 = R6/(WGS(I + 2)*HI4*HI5)
    EPS3 = EPS4
    EPS4 = (BETA4*HI3 + BETA5*HI5)/HI4
    D(I - 1) = BETA4
    P = H2DHI + BETA4 + BETA5 + EPS4
    B(I) = HI4 - P
    GSMO(I) = ETA3 + ETA4 + BETA4 + BETA4 + EPS3 + EPS4
  end do                  ! 40
50 HI3 = HI4
  HI4 = HI5
  ETA3 = ETA4
  ETA4 = HI4 + HI4 + H2/HI4
  BETA4 = BETA5
  EPS3 = EPS4
60 BETA5 = R6/(WGS(N)*HI4)
  EPS4 = (BETA4*HI3 + BETA5)/HI4
  GSMO(N2) = ETA3 + ETA4 + BETA4 + BETA4 + EPS3 + EPS4
  !  THE P.D. 5 BANDED SYMMETRIC MATRIX IS COMPLETE.
  !  THE SYSTEM OF LINEAR EQUATION WILL NOW BE SOLVED FOR THE C(I)'S.
  IF (N .GT. 3) GO TO 70
  C(1) = C(1)/GSMO(1)
  GO TO 150
70 IF (N .GT. 4) GO TO 80
  C(1) = (C(1)*GSMO(2) - C(2)*B(1))/(GSMO(1)*GSMO(2) - B(1)**2)
  C(2) = (C(2) - C(1)*B(1))/GSMO(2)
  GO TO 150
  !  THIS SOLVE THE 5 BANDED SYSTEM WHEN K=N-2.GT.3.
80 K = N2
  K1 = K - 1
  K2 = K1 - 1
  K3 = K2 - 1
  !  THE 5 BANDED MATRIX WILL NOW BE FACTORED.
  B(1) = B(1)/GSMO(1)
  D(1) = D(1)/GSMO(1)
  P = GSMO(1)*B(1)
  GSMO(2) = GSMO(2) - P*B(1)
  B(2) = (B(2) - P*D(1))/GSMO(2)
  IF (K .EQ. 3) GO TO 110
  D(2) = D(2)/GSMO(2)
  IF (K .EQ. 4) GO TO 100
  DO I = 3, K2            ! 90
    I1 = I - 1
    I2 = I1 - 1
    P = GSMO(I1)*B(I1)
    GSMO(I) = GSMO(I) - GSMO(I2)*(D(I2)**2) - P*B(I1)
    B(I) = (B(I) - P*D(I1))/GSMO(I)
    D(I) = D(I)/GSMO(I)
  end do                  ! 90
100 P = GSMO(K2)*B(K2)
  GSMO(K1) = GSMO(K1) - GSMO(K3)*(D(K3)**2) - P*B(K2)
  B(K1) = (B(K1) - P*D(K2))/GSMO(K1)
110 GSMO(K) = GSMO(K) - GSMO(K2)*(D(K2)**2) - GSMO(K1)*(B(K1)**2)
  !  FACTORIZATION COMPLETE.
  !  CARRY OUT FORWARD  AND BACKWARD SUBSTITUTION.
  C(2) = C(2) - B(1)*C(1)
  DO I = 3, K             ! 120
    I1 = I - 1
    I2 = I - 2
    C(I) = C(I) - B(I1)*C(I1) - D(I2)*C(I2)
  end do                  ! 120
  DO I = 1, K             ! 130
    C(I) = C(I)/GSMO(I)
  end do                  ! 130
  C(K1) = C(K1) - B(K1)*C(K)
  DO I = 2, K1            ! 140
    J = K - I
    C(J) = C(J) - B(J)*C(J + 1) - D(J)*C(J + 2)
  end do                  ! 140
  !  THE 5 BANDED SYSTEM HAS BEEN SOLVED.THE SOLUTION IS IN
  !  ARRAY C. THE COEFFICIENTS GSMO, B, C, AND D WILL NOW BE SET UP.
150 C(N) = 0.0
  D(N) = 0.0
  C(N1) = C(N2)
  HK1 = TNODE(N) - TNODE(N1)
  D(N1) = -C(N1)/(3.0*HK1)
  GSMO(N) = G(N) + R6*D(N1)/WGS(N)
  IF (N .EQ. 3) GO TO 170
  DO I = 2, N2            ! 160
    K = N - I
    K1 = K + 1
    HK2 = HK1
    HK1 = TNODE(K1) - TNODE(K)
    C(K) = C(K - 1)
    D(K) = (C(K1) - C(K))/(3.0*HK1)
    GSMO(K1) = G(K1) - R6*(D(K1) - D(K))/WGS(K1)
    B(K1) = (GSMO(K1 + 1) - GSMO(K1))/HK2 - HK2*(C(K1) + C(K1) + C(K1 + 1))/3.0
  end do                  ! 160
170 C(1) = 0.0
  HK2 = HK1
  HK1 = TNODE(2) - TNODE(1)
  D(1) = (C(2) - C(1))/(3.0*HK1)
  GSMO(2) = G(2) - R6*(D(2) - D(1))/WGS(2)
  GSMO(1) = G(1) - R6*D(1)/WGS(1)
  B(2) = (GSMO(3) - GSMO(2))/HK2 - HK2*(C(2) + C(2) + C(3))/3.0
  B(1) = (GSMO(2) - GSMO(1))/HK1 - HK1*(C(1) + C(1) + C(2))/3.0
  ! THE DISCRETE CUBIC SMOOTHING SPLINE IS NOW COMPLETE.
  RETURN

  ! THE TRIVIAL CASE WHEN N=2 IS HANDLED HERE.
180 GSMO(1) = G(1)
  GSMO(2) = G(2)
  B(1) = (G(2) - G(1))/(TNODE(2) - TNODE(1))
  C(1) = 0.0
  D(1) = 0.0

END subroutine



SUBROUTINE DCSINT(IENT, H, N, TNODE, G, END1, ENDN, B, C, D)
  !  THIS SUBROUTINE COMPUTES THE DISCRETE CUBIC SPLINE
  !  DEFINED ON THE INTERVAL (TNODE(1),TNODE(N)),WHICH INTER-
  !  POLATES THE DATA (TNODE(I),G(I)),I=1,2,...,N. WE REQUIRE
  !  THAT TNODE(I).LT.TNODE(I+1). END1 AND ENDN CONTAIN THE
  !  VALUES OF THE END CONDITIONS BEING USED.
  !
  !  IF IENT=1,THE FIRST CENTRAL DIVIDED DIFFERENCE END
  !  CONDITIONS ARE BEING USED.
  !
  !  IF IENT=2,THE SECOND CENTRAL DIVIDED DIFFERENCE END
  !  CONDITIONS ARE BEING USED.
  !
  !  IF IENT=3,THE PERODIC END CONDITIONS ARE BEING USED.
  !  FOR THIS CASE THE CONTENTS OF G(N), END1,AND ENDN ARE
  !  IGNORED.
  !
  !  FOR ALL THREE END CONDITIONS N MUST BE GREATER THAN OR
  !  EQUAL TO 2.
  !
  !  THE DISCRETE CUBIC SPLINE IS REPRESENTED BY PIECEWISE
  !  CUBIC POLYNOMIALS. FOR T IN THE INTERNAL (TNODE(I),TNODE
  !  (I+1)) THE CUBIC SPLINE IS
  !
  !     S(T)=G(I)+B(I)*(T-TNODE(I))
  !          +C(I)*(T-TNODE(I))**2
  !          +D(I)*(T-TNODE(I))**3
  !
  !  INPUT PARAMETERS (NONE OF THESE PARAMETERS
  !                          ARE CHANGED BY THIS SUBROUTINE.)
  !
  !  IENT  - SPECIFIES END CONDITIONS WHICH ARE IN EFFECT.
  !  H     - THE STEP SIZE USED FOR THE DISCRETE CUBIC SPLINE.
  !  N     - NUMBER OF NODES (TNODE) AND DATA VALUES (G).
  !          (N.GE.2)
  !  TNODE - REAL ARRAY CONTAINING THE NODES (TNODE(I).LT.
  !          TNODE(I+1)).
  !  G     - REAL ARRAY CONTAINING THE INTERPOLATING DATA.
  !  END1  - END CONDITION VALUE AT TNODE(1).
  !  ENDN  - END CONDITION VALUE AT TNODE(N).
  !
  !  OUTPUT PARAMETERS
  !
  !  B     - REAL ARRAY CONTAINING COEFFICIENTS OF
  !          (T-TNODE(I)),I=1,2,....,N-1.
  !  C     - REAL ARRAY CONTAINING COEFFICIENTS OF
  !          (T-TNODE(I))**2,I=1,2,...,N-1.
  !  D     - REAL ARRAY CONTAINING COEFFICIENTS OF
  !          (T-TNODE(I))**3,I=1,2,...,N-1.
  !
  !  SPECIAL CASES ARE ACCOUNTED FOR HERE.

  USE module_kind_variables, ONLY: kind_double
  integer :: N 
  real(kind_double)   :: TNODE(N), G(N), B(N), C(N), D(N)
  real(kind_double)   :: H, H2, END1, ENDN, GN, HI, H2DHI, G2, &
                         ETA2, ETA1, DEN, T, C1, C2, G1, GAMMA, CT, &
                         BS1   
  integer :: IENT, N1, N2, I,  L, ITRANS, L1, LI 


  IF (N .EQ. 2 .AND. IENT .EQ. 3) GO TO 220
  H2 = H*H
  N1 = N - 1
  IF (N .EQ. 2 .AND. IENT .EQ. 2) GO TO 180
  N2 = N1 - 1
  !  THE SYMMETRIC TRIDIAGONAL(OR NEAR TRIDIAGONAL) LINEAR
  !  SYSTEM WILL NOW BE SET UP FOR THE APPROPRIATE END
  !  CONDITIONS
  HI = TNODE(2) - TNODE(1)
  H2DHI = H2/HI
  ETA2 = HI + HI + H2DHI
  GO TO(10, 40, 60), IENT
  !  IENT=1 - FIRST CENTRAL DIVIDED DIFFERENCE END CONDITONS
  !  SET UP.
10 B(1) = ETA2
  D(1) = HI - H2DHI
  G2 = (G(2) - G(1))/HI
  C(1) = 3.0*(G2 - END1)
  IF (N .EQ. 2) GO TO 30
  DO I = 2, N1            ! 20
    ETA1 = ETA2
    G1 = G2
    HI = TNODE(I + 1) - TNODE(I)
    H2DHI = H2/HI
    ETA2 = HI + HI + H2DHI
    B(I) = ETA1 + ETA2
    D(I) = HI - H2DHI
    G2 = (G(I + 1) - G(I))/HI
    C(I) = 3.0*(G2 - G1)
  end do                  ! 20CONTINUE
30 B(N) = ETA2
  C(N) = 3.0*(ENDN - G2)
  L = N
  !  SET UP FOR (1) FIRST CENTRAL DIVIDED DIFFERENCE END
  !  CONDITING COMPLETE.THE LINEAR EQUATIONS WILL BE NOW
  !  SOLVED.
  GO TO 110
  !  IENT=2 - SECOND CENTRAL DIVIDED DIFFERENCE END CONDITIONS
  !  SET UP.
40 GAMMA = HI - H2DHI
  G2 = (G(2) - G(1))/HI
  DO I = 1, N2            ! 50
    ETA1 = ETA2
    G1 = G2
    HI = TNODE(I + 2) - TNODE(I + 1)
    H2DHI = H2/HI
    ETA2 = HI + HI + H2DHI
    B(I) = ETA1 + ETA2
    D(I) = HI - H2DHI
    G2 = (G(I + 2) - G(I + 1))/HI
    C(I) = 3.0*(G2 - G1)
  end do                  ! 50
  C(1) = C(1) - GAMMA*END1/2.0
  HI = TNODE(N) - TNODE(N1)
  GAMMA = HI - H2/HI
  C(N2) = C(N2) - GAMMA*ENDN/2.0
  !  STEP UP FOR (2) SECOND CENTRAL DIVIDED DIFFERENCE
  !  END CONDITIONS COMPLETE. THE LINEAR EQUATIONS WILL NOW
  !  BE SOLVED.
  IF (N2 .EQ. 1) GO TO 100
  L = N2
  GO TO 110
  !  IENT=3 - PERIODIC END CONDITIONS SET UP.
60 B(1) = ETA2
  D(1) = HI - H2DHI
  DO I = 2, N1            ! 70
    ETA1 = ETA2
    HI = TNODE(I + 1) - TNODE(I)
    H2DHI = H2/HI
    ETA2 = HI + HI + H2DHI
    B(I) = ETA1 + ETA2
    D(I) = HI - H2DHI
    C(I) = 0.0
  end do                  ! 70
  CT = D(N1)
  C(1) = CT
  C(N1) = CT
  B(1) = B(1) + ETA2 - CT
  B(N1) = B(N1) - CT
  L = N1
  ITRANS = 1
  GO TO 120
80 G1 = (G(N1) - G(N2))/(TNODE(N1) - TNODE(N2))
  G2 = (G(1) - G(N1))/(TNODE(N) - TNODE(N1))
  DEN = (1.0 + C(1) + C(N1))
  CT = 3.0*(G2 - G1)
  BS1 = CT*C(N1)
  C(N1) = CT
  DO I = 1, N2            ! 90
    HI = TNODE(I + 1) - TNODE(I)
    G1 = G2
    G2 = (G(I + 1) - G(I))/HI
    CT = 3.0*(G2 - G1)
    BS1 = BS1 + CT*C(I)
    C(I) = CT
  end do                  ! 90
  BS1 = BS1/DEN
  C(1) = C(1) - BS1
  C(N1) = C(N1) - BS1
  ITRANS = 0
  GO TO 140
  !  THE SET UP AND MOST OF THE DETAILS FOR SOLVING THE LINEAR EQUQTION FOR (3)
  !  THE PERIODIC END CONDITIONS ARE COMPLETED.
  !  THE LINEAR EQUATION ARE SOLVED HERE.
100 C(2) = C(1)/B(1)
  GO TO 180
110 ITRANS = 0
120 L1 = L - 1
  DO I = 1, L1            ! 130
    T = D(I)/B(I)
    B(I + 1) = B(I + 1) - T*D(I)
    D(I) = T
  end do                  ! 130
  !  THE LINEAR EQUATION SOLVER IS ENTERED AT THIS POINT
  !  IF THE LU FACTORIZATION HAS ALREADY BEEN DONE.
140 DO I = 1, L1            ! 150
    C(I + 1) = C(I + 1) - D(I)*C(I)
  end do                  ! 150
  C(L) = C(L)/B(L)
  DO I = 1, L1            ! 160
    LI = L - I
    C(LI) = C(LI)/B(LI) - D(LI)*C(LI + 1)
  end do                  ! 160
  IF (ITRANS .GE. 1) GO TO 80
  !  THE LINEAR SYSTEM HAS BEEN SOLVE FOR THE C-VECTOR
  IF (IENT .EQ. 3) C(N) = C(1)
  IF (IENT .NE. 2) GO TO 190
  DO I = 1, N2            ! 170
    LI = N - I
    C(LI) = C(LI - 1)
  end do                  ! 170
180 C(1) = END1/2.0
  C(N) = ENDN/2.0
190 C1 = C(1)
  C2 = C(2)
  HI = TNODE(2) - TNODE(1)
  IF (N .EQ. 2) GO TO 210
  DO I = 1, N2            ! 200
    B(I) = (G(I + 1) - G(I))/HI - HI*(C1 + C1 + C2)/3.0
    D(I) = (C2 - C1)/(3.0*HI)
    HI = TNODE(I + 2) - TNODE(I + 1)
    C1 = C2
    C2 = C(I + 2)
  end do                  ! 200
210 GN = G(1)
  IF (IENT .NE. 3) GN = G(N)
  B(N1) = (GN - G(N1))/HI - HI*(C1 + C1 + C2)/3.0
  D(N1) = (C2 - C1)/(3.0*HI)
  ! THE INTERPOLATING DISCRETE CUBIC SPLINE HAS BEEN CONSTRUCTED.
  RETURN

  ! THE FOLLOWING HANDLES THE TRIVIAL PERIODIC CASE
  ! (IENT.EQ.3) WHEN N.EQ.2.
220 B(1) = 0.0
  C(1) = 0.0
  D(1) = 0.0

END subroutine

end module former_mat_util
