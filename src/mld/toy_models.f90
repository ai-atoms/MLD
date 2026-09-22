! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   MiLaDy (Machine Learning Dynamics) 
! HND X   
! HND X
! HND X   Milady was written and designed by Mihai-Cosmin Marinica, Alexandra M. Goryaeva 
! HND X   Copyright 2015-2023.
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

module toy_models

  implicit none

  integer, parameter   :: double = 8
  real(kind=kind(1.d0)), parameter :: pi = 4.d0*datan(1.d0)

  interface generate_random_gaussian
    module procedure random_number_vector, &
      random_number_scalar
  end interface generate_random_gaussian

  interface generate_random_uniform
    module procedure random_vector_uniform, &
      random_scalar_uniform
  end interface generate_random_uniform

contains

  subroutine random_vector_uniform(gau, dim_gau)
    integer, intent(in)  :: dim_gau
    real(kind=kind(1.d0)), dimension(dim_gau), intent(out)   :: gau
    real(kind=kind(1.d0))      :: u1
    integer  :: i

    do i = 1, dim_gau
      call random_number(u1)
      gau(i) = u1
    end do
  end subroutine random_vector_uniform


  subroutine random_scalar_uniform(rvalue)
    real(kind=kind(1.d0)), intent(out)     :: rvalue

    call random_number(rvalue)
  end subroutine random_scalar_uniform


  subroutine random_number_vector(gau, dim_gau)    ! random gaussian noise between -1 and 1.
    integer, intent(in)  :: dim_gau
    real(kind=kind(1.d0)), dimension(dim_gau), intent(out)   :: gau
    real(kind=kind(1.d0))      :: u1, u2, b1
    integer  :: i           ! ,iatom,i

    gau(:) = 0
    do i = 1, dim_gau
      call random_number(u1)
      call random_number(u2)
      b1 = sqrt(-2.*log(u1))*cos(2.*pi*u2)
      gau(i) = b1
    end do
  end subroutine random_number_vector



  subroutine random_number_scalar(value)           ! random gaussian noise between -1 and 1
    real(kind=kind(1.d0)), intent(out)     :: value
    real(kind=kind(1.d0))      :: u1, u2, b1

    value = 0
    call random_number(u1)
    call random_number(u2)
    b1 = sqrt(-2.*log(u1))*cos(2.*pi*u2)
    value = b1
  end subroutine random_number_scalar


  subroutine toy_nD(x, y, L, dim_xdesc)
    integer, intent(in)  :: dim_xdesc
    real(kind=kind(1.d0)), intent(in)      :: x(dim_xdesc), L
    real(kind=kind(1.d0)), intent(out)     :: y
    real(kind=kind(1.d0))      :: temp
    integer  :: i

    if (dim_xdesc == 1) y = dsin(2.d0*pi*x(1)/L)
    if (dim_xdesc == 2) y = dsin(x(1)*2.d0*pi/L)*dcos(x(2)*2.d0*pi/L)
    if (dim_xdesc > 2) then
      temp = 0
      do i = 3, dim_xdesc
        temp = temp + dexp((x(i) - L/dble(i - 2))**2/(2.D0*L))
      end do
      y = dsin(x(1)*2.d0*pi/L)*dcos(x(2)*2.d0*pi/L) + temp + 10000.0
    end if
  end subroutine toy_nD

end module toy_models



subroutine set_toy_model()
  ! use module_kind_variables, ONLY:  kind_double
  use ml_in_ndm_module

  use temporary_data_cov, ONLY: yfunc, yfunc_nd, xdesc, dim_data, dim_xdesc, dim_yfunc, yfunc_average
  ! use extrapolation
  use k_cross_validation
  use set_limits

  implicit none

  dim_xdesc = nd_fingerprint
  if (allocated(yfunc)) deallocate (yfunc)
  allocate (yfunc(dim_data))
  if (allocated(xdesc)) deallocate (xdesc)
  allocate (xdesc(dim_xdesc, dim_data))

  ! call build_database_with_function(dim_data,dim_xdesc,yfunc, xdesc)
  if (dim_yfunc /= 1) then
    if (rangml == 0) write (6, *) 'ML: < ml > ERROR for the toy_model case dim_yfunc should be 1', dim_yfunc
    stop
  end if

  if (allocated(yfunc_average)) deallocate (yfunc_average)
  allocate (yfunc_average(dim_yfunc))
  if (allocated(yfunc_nd)) deallocate (yfunc_nd)
  allocate (yfunc_nd(dim_yfunc, dim_data))

  ! TODOtoy what is that
  yfunc_nd(dim_yfunc, 1:dim_data) = yfunc(1:dim_data)
  call recentrate_database(dim_data, dim_yfunc, yfunc_nd, yfunc_average)
  yfunc(:) = yfunc_nd(dim_yfunc, :)

  if (rangml == 0) write (6, *) 'set_toy_model: dim_data ', dim_data
end subroutine set_toy_model
