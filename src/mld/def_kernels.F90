
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

module def_kernels
  USE module_kind_variables, ONLY: kind_double
  implicit none

  ! real(kind_double), dimension(:,:), allocatable :: xdesc ! first  index is the dimension of the descriptor
  real(kind_double)    :: length_kse
  real(kind_double)    :: sigma_kse
  real(kind_double), parameter     :: l_kou = 1.d0

  integer, parameter   :: kernel_se = 1, &
                          kernel_ou = 2, &
                          kernel_mc = 3, &
                          kernel_so = 4
  integer  :: kappa_soap

contains

  real(kind_double) function kernel_square_exp(a_local, b_local, n_dimension_desc)
    implicit none
    integer  :: n_dimension_desc
    real(kind_double), dimension(n_dimension_desc)     :: a_local, b_local

    kernel_square_exp = sigma_kse*dexp(-real(SUM((a_local(:) - b_local(:))*(a_local(:) - b_local(:))))/(2.d0*length_kse))
  end function kernel_square_exp


  real(kind_double) function d_length_kernel_square_exp(a_local, b_local, n_dimension_desc)
    implicit none
    integer  :: n_dimension_desc
    real(kind_double), dimension(n_dimension_desc)     :: a_local, b_local
    real(kind_double)    :: temp

    temp = SUM((a_local(:) - b_local(:))*(a_local(:) - b_local(:)))/(2.d0*length_kse)

    d_length_kernel_square_exp = sigma_kse*dexp(-temp) - temp/length_kse
  end function d_length_kernel_square_exp



  real(kind_double) function kernel_ornstein_uhlenbeck(a_local, b_local, n_dimension_desc)
    implicit none
    integer  :: n_dimension_desc
    real(kind_double), dimension(n_dimension_desc)     :: a_local, b_local

    kernel_ornstein_uhlenbeck = dexp(-dsqrt((SUM((a_local(:) - b_local(:))*(a_local(:) - b_local(:)))))/dble(l_kou))
  end function kernel_ornstein_uhlenbeck


  real(kind_double) function kernel_matern_class(a_local, b_local, n_dimension_desc)
    implicit none
    integer  :: n_dimension_desc
    real(kind_double), dimension(n_dimension_desc)     :: a_local, b_local

    kernel_matern_class = dexp(-dsqrt((SUM((a_local(:) - b_local(:))*(a_local(:) - b_local(:)))))/dble(l_kou))
  end function kernel_matern_class


  real(kind_double) function kernel_soap(a_local, b_local, n_dimension_desc)
    implicit none
    integer  :: n_dimension_desc
    real(kind_double), dimension(n_dimension_desc)     :: a_local, b_local

    kernel_soap = dot_product(a_local(:), b_local(:))**kappa_soap
  end function kernel_soap

end module def_kernels
