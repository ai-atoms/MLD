! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   MiLaDy (Machine Learning Dynamics) 
! HND X   
! HND X
! HND X   Milady was written and designed by Mihai-Cosmin Marinica, Alexandra M. Goryaeva 
! HND X   Copyright 2015-2024.
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

module lbfgs_teta

  real(kind(1.d0))     :: eps, xtol
  !     EPS     is a positive DOUBLE PRECISION variable that must be set by
  !             the user, and determines the accuracy with which the solution
  !             is to be found. The subroutine terminates when
  !
  !                         ||G|| < EPS max(1,||X||),
  !
  !             where ||.|| denotes the Euclidean norm.
  !
  !     XTOL    is a  positive DOUBLE PRECISION variable that must be set by
  !             the user to an estimate of the machine precision (e.g.
  !             10**(-16) on a SUN station 3/60). The line search routine will
  !             terminate if the relative width of the interval of uncertainty
  !             is less than XTOL.

  logical  :: diagco

  !     DIAGCO  is a LOGICAL variable that must be set to .TRUE. if the
  !             user  wishes to provide the diagonal matrix Hk0 at each
  !             iteration. Otherwise it should be set to .FALSE., in which
  !             case  LBFGS will use a default value described below. If
  !             DIAGCO is set to .TRUE. the routine will return at each
  !             iteration of the algorithm with IFLAG=2, and the diagonal
  !             matrix Hk0  must be provided in the array DIAG.

  integer  :: iflag, icall, call_max_lbfgs, m_hess, iprint(2)
  real(kind(1.d0)), dimension(:), allocatable  :: diag_lbfgs, w_lbfgs

  !internal set-up of the LBFGS subroutine ...
  real(kind(1.d0))     :: GTOL, STPMIN, STPMAX
  integer  :: MP, LP

end module lbfgs_teta


! call driver_lbfgs(20)


module lbfgs_mod
  ! useful for test subroutines signatures
  implicit none
contains
  include 'lbfgs.f'       ! yeah for old fortran code
end module


subroutine driver_lbfgs(dim_n)
  use lbfgs_teta
  use lbfgs_mod

  implicit none

  real(kind(1.d0)), dimension(:), allocatable  :: x_coord, grad_functx
  real(kind(1.d0))     :: functx
  integer  :: dim_n

  ! The driver for LBFGS must always declare LB2 as EXTERNAL
  interface
    subroutine func_and_grad(y, y_grad, x, n)
      implicit none
      integer, intent(in)  :: n
      real(kind(1.d0)), dimension(n), intent(in)   :: x(n)
      real(kind(1.d0)), intent(out)    :: y
      real(kind(1.d0)), dimension(n), intent(out)  :: y_grad(n)
    end subroutine func_and_grad
  end interface


  if (allocated(x_coord)) deallocate (x_coord); allocate (x_coord(dim_n))
  if (allocated(grad_functx)) deallocate (grad_functx); allocate (grad_functx(dim_n))

  call set_up_the_lbfgs(dim_n)
  call put_parameters_into_x(x_coord, dim_n)

  do icall = 1, call_max_lbfgs
    call func_and_grad(functx, grad_functx, x_coord, dim_n)
    call LBFGS(dim_n, m_hess, x_coord, functx, grad_functx, diagco, diag_lbfgs, iprint, eps, xtol, w_lbfgs, iflag)
    if (iflag .le. 0) exit
  end do

  call put_x_into_parameters(x_coord, dim_n)

  if (icall .ge. call_max_lbfgs) then
    write (*, *) 'WARNING: the convergency is not achieved in LBFGS in call_max_lbfgs gradient claculations', call_max_lbfgs
    if (iflag < 0) then
      write (*, *) 'LBFGS failed with iflag and info', iflag
      stop
    end if
    if (iflag == 1) then
      write (*, *) 'You should increase the maximum number of iterations call_max_lbfgs', call_max_lbfgs
    end if
  end if

end subroutine driver_lbfgs

subroutine put_parameters_into_x(x, n)
  use module_nlinear, only: order_nlinear, alpha_nl
  use temporary_data_cov, only: dim_xdesc
  implicit none

  integer, intent(in)  :: n
  real(kind(1.d0)), dimension(n), intent(inout)      :: x(n)

  if (order_nlinear == 2) then
    x(1:dim_xdesc) = alpha_nl(1, 1:dim_xdesc)
    x(dim_xdesc + 1:n) = alpha_nl(2, 1:dim_xdesc)
  end if

end subroutine put_parameters_into_x


subroutine put_x_into_parameters(x, n)
  use module_nlinear, only: order_nlinear, alpha_nl
  use temporary_data_cov, only: dim_xdesc
  implicit none

  integer, intent(in)  :: n
  real(kind(1.d0)), dimension(n), intent(in)   :: x(n)

  if (order_nlinear == 2) then
    alpha_nl(1, 1:dim_xdesc) = x(1:dim_xdesc)
    alpha_nl(2, 1:dim_xdesc) = x(1 + dim_xdesc:2*dim_xdesc)
  end if

end subroutine put_x_into_parameters



subroutine func_and_grad(y, y_grad, x, n)
  use module_objective_nl, only: jobj_T, jobj_T_d
  use temporary_data_cov, only: dim_xdesc
  use module_nlinear, only: order_nlinear, alpha_nl
  implicit none

  integer, intent(in)  :: n
  real(kind(1.d0)), dimension(n), intent(in)   :: x(n)
  real(kind(1.d0)), intent(out)    :: y
  real(kind(1.d0)), dimension(n), intent(out)  :: y_grad(n)

  if (order_nlinear == 2) then
    alpha_nl(1, 1:dim_xdesc) = x(1:dim_xdesc)
    alpha_nl(2, 1:dim_xdesc) = x(1 + dim_xdesc:2*dim_xdesc)
  end if

  call evaluate_objective
  y = jobj_T
  y_grad = jobj_T_d

end subroutine func_and_grad



subroutine func_and_grad_old(y, y_grad, x, n)
  implicit none

  integer, intent(in)  :: n
  real(8), dimension(n), intent(in)      :: x(n)
  real(8), intent(out) :: y
  real(8), dimension(n), intent(out)     :: y_grad(n)
  real(8)  :: tmp
  integer  :: i

  tmp = 0.d0
  do i = 1, n
    tmp = tmp + (x(i) - dble(i)**2)**4
    y_grad(i) = 4.d0*(x(i) - dble(i)**2)**3
  end do
  y = tmp

end subroutine func_and_grad_old


subroutine set_up_the_lbfgs(N)
  use lbfgs_teta
  use module_lbfgs_input, only: lbfgs_xtol, lbfgs_eps, &
                                lbfgs_m_hess, lbfgs_max_steps, &
                                lbfgs_print, lbfgs_gtol

  implicit none

  integer, intent(in)  :: N
  integer  :: nwork

  if (allocated(diag_lbfgs)) deallocate (diag_lbfgs)
  allocate (diag_lbfgs(n))
  ! orig to 5.
  !m_hess=100
  m_hess = lbfgs_m_hess
  nwork = N*(2*m_hess + 1) + 2*m_hess
  if (allocated(w_lbfgs)) deallocate (w_lbfgs); allocate (w_lbfgs(nwork))
  eps = lbfgs_eps
  !eps = 1.0D-3
  xtol = lbfgs_xtol
  !XTOL= 1.0D-16
  call_max_lbfgs = lbfgs_max_steps
  !call_max_lbfgs=20000

  DIAGCO = .FALSE.
  !COSBFGS was 1.0D-16
  ICALL = 0
  IFLAG = 0

  !IPRINT(1)=100
  !IPRINT(2)=0
  iprint(1) = lbfgs_print(1)
  iprint(2) = lbfgs_print(2)
  !GTOL org was on 9.d-1.
  !GTOL 1.d-1 nice
  !GTOL=4.d-1
  gtol = lbfgs_gtol
  STPMIN = 1.0D-20
  STPMAX = 1.0D+20
  MP = 6
  LP = 6

end subroutine set_up_the_lbfgs



