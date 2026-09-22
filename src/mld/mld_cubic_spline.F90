! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   MiLaDy (Machine Learning Dynamics) 
! HND X   
! HND X
! HND X   Milady was written and designed by Mihai-Cosmin Marinica, Alexandra M. Goryaeva 
! HND X   Copyright 2015-2025.
! HND X
! HND X   Portions of MiLady were written by Wesley UnnToc, Clovis Lapointe, Anruo Zhong, 
! HND X   Jacopo Baima, Anida Khizar, Christian van Wambeke 
! HND X
! HND X   MiLaDy is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   MiLaDy is distributed at https://github.com/ai-atoms/milady in the hope that it 
! HND X   will be useful for non-commercial academic research, but WITHOUT ANY WARRANTY;  
! HND X   without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
! HND X   PURPOSE. See the ASL for more details.
! HND X
! HND X   You should have received a copy of the ASL along with this program
! HND X   (e.g. in a LICENSE.md file); if not, you can write to the original licensors,
! HND X   Mihai-Cosmin Marinica and Alexandra M. Goryaeva. The ASL is also published at
! HND X   http://github.com/ai-atoms/milady/ASL
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   A.M. Goryaeva et al. Comput. Mater. Sci. 166, 200-209 (2019)
! HND X   A.M. Goryaeva et al. Phys. Rev. Mater. 5: 103803 (2021)
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "../MLD_MACROS.INC"

module module_spline_interpolation_mine
    use, intrinsic :: iso_fortran_env, dp=>real64
    implicit none

    ! Define a type for the cubic spline
    type :: CubicSpline
        integer :: n = 0                             ! Number of points
        real(dp)  ::  delta                          ! Spacing between points
        real(dp) :: r_cut_in, r_cut_out              ! Cutoff radii of splinning ... 
        real(dp), allocatable :: x(:), y(:), y_prime(:)  ! Points and function values
        real(dp), allocatable :: a(:), b(:), c(:), d(:)  ! Spline coefficients
    contains
        procedure :: init => initSpline
        procedure :: compute => computeSpline
        procedure :: evaluate => evaluateSpline
        procedure :: derivative => derivativeSpline
    end type CubicSpline

contains

subroutine initSpline(this, r_cut_in, r_cut_out, npoints, xgrid, yfunc, d_yfunc)
    class(CubicSpline), intent(inout) :: this
    real(dp), intent(in) :: r_cut_in, r_cut_out
    integer, intent(in) :: npoints
    real(dp), dimension(:), allocatable, intent(out) :: yfunc, d_yfunc, xgrid  
    real(dp) ::  x_val
    integer :: i

    if (npoints < 2) then
        error stop "CubicSpline: At least two points are required"
    end if

    ! Initialize the number of points and allocate arrays
    this%n = npoints
    if (allocated(this%x))   deallocate(this%x)   ; allocate(this%x(this%n))
    if (allocated(this%y))   deallocate(this%y)   ; allocate(this%y(this%n))
    if (allocated(this%y_prime)) deallocate(this%y_prime) ; allocate(this%y_prime(this%n))
    !allocate(this%x(this%n), this%y(this%n), this%y_prime(this%n))
    if (allocated(this%a))   deallocate(this%a)   ; allocate(this%a(this%n))
    if (allocated(this%b))   deallocate(this%b)   ; allocate(this%b(this%n))
    if (allocated(this%c))   deallocate(this%c)   ; allocate(this%c(this%n))
    if (allocated(this%d))   deallocate(this%d)   ; allocate(this%d(this%n-1))
    !allocate(this%a(this%n), this%b(this%n), this%c(this%n), this%d(this%n-1))
    this%r_cut_in = r_cut_in
    this%r_cut_out = r_cut_out
    this%delta = (this%r_cut_out - this%r_cut_in) / dble(npoints - 1)

    ! Fill x, y, and y_prime arrays
    do i = 1, this%n
        x_val = this%r_cut_in + dble(i-1)*this%delta
        this%x(i) = x_val
    end do

    if (allocated(xgrid))   deallocate(xgrid)   ; allocate(xgrid(this%n))
    if (allocated(yfunc))   deallocate(yfunc)   ; allocate(yfunc(this%n))
    if (allocated(d_yfunc)) deallocate(d_yfunc) ; allocate(d_yfunc(this%n))
    xgrid(:) = this%x(:) 

end subroutine initSpline


subroutine computeSpline(this, yfunc, d_yfunc)
  class(CubicSpline), intent(inout) :: this 
  real(dp), dimension(:), allocatable, intent(in) :: yfunc, d_yfunc 
  real(dp) :: fp_in, fp_out
  logical :: is_constant
  integer :: i 

  ! Check if yfunc is constant and d_yfunc is zero
  is_constant = .true.
  do i = 2, this%n
      if  (  &
           !abs((yfunc(i) - yfunc(1)).gt.1.d-25) .or.(abs(d_yfunc(i) - 0.0_dp).gt.1.d-20)) then
           abs(yfunc(i) - yfunc(1)) > 1.d-25 .or. abs(d_yfunc(i)) > 1.d-20) then
          is_constant = .false.
          exit
      end if
  end do
  if (is_constant) then
      this%y(:) = yfunc(:)
      this%y_prime(:) = d_yfunc(:)
      this%a(:) = yfunc(:)
      this%b(:) = 0.0_dp
      this%c(:) = 0.0_dp
      this%d(:) = 0.0_dp
      return
  end if
  

  this%y(:) = yfunc(:) 
  this%y_prime(:) = d_yfunc(:) 

  fp_in = this%y_prime(1)
  fp_out = this%y_prime(this%n)
  call computeClampedCoefficients(this%x, this%y, this%n, this%a, this%b, this%c, this%d, fp_in, fp_out)

end subroutine computeSpline 

subroutine computeClampedCoefficients(x, y, n, a, b, c, d, fp1, fpn)
    integer, intent(in) :: n
    real(dp), dimension(n), intent(in) :: x, y
    real(dp), dimension(n), intent(out) :: a, b, c
    real(dp), dimension(n-1), intent(out) :: d
    real(dp), intent(in) :: fp1, fpn

    !local variables ... 
    real(dp), dimension(n-1) :: h
    real(dp), dimension(n) :: alpha, l, mu, z
    integer :: i

    a = y
    do i = 1, n-1
        h(i) = x(i+1) - x(i)
    end do

    ! Adjusted for clamped splines
    alpha(1) = 3_dp*(a(2) - a(1))/h(1) - 3_dp*fp1
    alpha(n) = 3_dp*fpn - 3_dp*(a(n) - a(n-1))/h(n-1)

    do i = 2, n-1
        alpha(i) = 3_dp*(a(i+1) - a(i))/h(i) - 3_dp*(a(i) - a(i-1))/h(i-1)
    end do

    l(1) = 2.0_dp*h(1)
    !old mu(1) = 0.5_dp
    mu(1) = h(1)/l(1)
    z(1) = alpha(1) / l(1)

    do i = 2, n-1
        !TODOace_crade mais efficace. 
        if (abs(z(i-1)) < 1.0d-290) then
            z(i-1) = 0.0_dp
        end if
        l(i) = 2.0_dp*(x(i+1) - x(i-1)) - h(i-1)*mu(i-1)
        mu(i) = h(i) / l(i)
        z(i) = (alpha(i) - h(i-1)*z(i-1)) / l(i)
    end do

    l(n) = h(n-1)*(2.0_dp - mu(n-1))
    z(n) = (alpha(n) - h(n-1)*z(n-1)) / l(n)
    c(n) = z(n)
    b(n)= 0.0_dp

    ! Back substitution loop adjusted for clamped spline
    do i = n-1, 1, -1
        !TODOace_crade mais efficace. 
        if (abs(c(i+1)) < 1.0d-290) then
            c(i+1) = 0.0_dp
        end if
        c(i) = z(i) - mu(i)*c(i+1)
        b(i) = (a(i+1) - a(i))/h(i) - h(i)*(c(i+1) + 2.0_dp*c(i))/3.0_dp
        d(i) = (c(i+1) - c(i))/(3.0_dp*h(i))
    end do
end subroutine computeClampedCoefficients

real(dp) function evaluateSpline(this, rr) result(y_val)
    !use, intrinsic :: iso_fortran_env, dp=>real64 
    use, intrinsic :: ieee_arithmetic
    class(CubicSpline), intent(in) :: this
    real(dp), intent(in) :: rr
    integer :: i
    real(dp) :: h

    y_val = 0.0
    if (rr < this%r_cut_in .or. rr > this%r_cut_out) then
        write(*,*)  "!!!!!!!!!!!!!!!! evaluate rr is outside the range of x", rr, this%r_cut_in, this%r_cut_out
        return
    endif

    ! Calculate the interval index directly
    i = min(int(floor((rr - this%r_cut_in) / this%delta) + 1), this%n - 1)
    h = rr - this%x(i)
    ! Evaluate the cubic spline polynomial
    y_val = this%a(i) + h*(this%b(i) + h*(this%c(i) + this%d(i)*h))

    !TODOace_crade mais efficace.
    if (dabs(y_val) < epsilon(1.d0)*100.d0) then
        !write(*,*) "❌ Spline output overflow or invalid for rr =", rr
        y_val = 0.0_dp  ! Or use a sentinel value or last valid value
    end if

end function evaluateSpline

real(dp) function derivativeSpline(this, rr) result(dy_val)
    class(CubicSpline), intent(in) :: this
    real(dp), intent(in) :: rr
    integer :: i
    real(dp) :: h

    dy_val = 0.0
  
    if (rr < this%r_cut_in .or. rr > this%r_cut_out) then
        ! rr is outside the range of x
        write(*,*)  "!!!!!!!!!!!!!!!! derivative rr is outside the range of x", rr, this%r_cut_in, this%r_cut_out
        return
    endif

    ! Calculate the interval index directly
    i = min(int(floor((rr - this%r_cut_in) / this%delta) + 1), this%n - 1)
    h = rr - this%x(i)
    ! Evaluate the derivative of the cubic spline polynomial
    dy_val = this%b(i) + 2.0_dp*this%c(i)*h + 3.0_dp*this%d(i)*h**2
        !TODOace_crade mais efficace.
    if (dabs(dy_val) < epsilon(1.d0)*100.d0) then
        !write(*,*) "❌ Spline output overflow or invalid for rr =", rr
        dy_val = 0.0_dp  ! Or use a sentinel value or last valid value
    end if
end function derivativeSpline

end module module_spline_interpolation_mine
