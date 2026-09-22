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

module temporary_data_cov
  use module_kind_variables, only: kind_double
  implicit none

  integer  :: dim_xdesc
  integer  :: dim_xdesc_patch
  ! for composed descriptors.
  integer  :: dim_xdesc1, dim_xdesc2
  integer  :: dim_extra, dim_train, dim_valid

  ! dim_data - dimension of the database
  ! dim_data_train + dim_data_valid= dim_data
  integer  :: dim_data, dim_data_train, dim_data_test, dim_data_train_local
  integer  :: dim_data_constraints

  real(kind=kind(1.d0)), dimension(:), allocatable   :: yfunc_average

  integer, parameter   :: dim_yfunc = 1
  real(kind=kind(1.d0)), dimension(:), allocatable   :: yfunc               ! vector of dimension M, size of the database
  real(kind=kind(1.d0)), dimension(:), allocatable   :: yfunc_train, yfunc_valid                     ! vector of dimension M, size of the database
  real(kind=kind(1.d0)), dimension(:), allocatable   :: y_test
  real(kind=kind(1.d0)), dimension(:), allocatable   :: y_extra             ! vector of dimension M, size of the database

  !if dim_yfunc > 1 probably the following functions should be used:

  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: yfunc_nd      ! vector of dimension YxM, M = size of the database, Y the componenets of the  y function

  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: matfor
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: xdesc         ! vetor of dimension (dim_desc, M) desc size x size of the database
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: xdesc_valid, xdesc_test   ! vetor of dimension (dim_desc, P) desc size x size of the database
  real(kind=kind(1.d0)), dimension(:), allocatable   :: error_train, error_valid, error_test
  real(kind=kind(1.d0)), dimension(:), allocatable   :: mlocal, kbfunc, eigen_values
  real(kind=kind(1.d0)), dimension(:), allocatable   :: d_mlocal            ! the values of the derivatives respect with the length of kernel
  integer, dimension(:), allocatable     :: u_local, v_local

  integer  :: i_start_cov, i_final_cov, i_local_cov
  integer, dimension(:), allocatable     :: neigh_cov, list_cov

end module temporary_data_cov




module k_cross_validation
  implicit none

  logical  :: kcross
  integer  :: n_kcross
  integer, allocatable, dimension(:)     :: dim_kcell_valid, dim_kcell_train, imin_kcell, imax_kcell, rdm

contains

  subroutine shuffle(n, a)

    integer, intent(in)  :: n
    integer, dimension(n), intent(out)     :: a
    integer  :: i, randpos, temp
    real     :: r

    a = (/(i, i=1, n)/)
    do i = n, 2, -1
      call random_number(r)
      randpos = int(r*i) + 1
      temp = a(randpos)
      a(randpos) = a(i)
      a(i) = temp
    end do
  end subroutine shuffle

end module k_cross_validation
