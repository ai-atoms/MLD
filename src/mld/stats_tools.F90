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

#include "../MLD_MACROS.INC"

module math

  use mld_logger
  use ml_in_ndm_module, only: rangml
  use module_kind_variables, only: kind_double

contains

  subroutine sample_1D_gaussian_scalar_box_muller(x1, mu_1d, sigma_1d)
    ! return a scalar x1 that smaple a 1D Gaussian N(mu_1d, sigma_1d)
    implicit none
    real(kind_double), intent(in)      :: mu_1d, sigma_1d
    real(kind_double), intent(out)     :: x1
    real(kind_double)      :: u1, u2, low_limit, mag, x2
    real(kind_double), parameter :: two_pi = 6.283185307179586d0          ! numpy value
  
    low_limit = 2.d0*epsilon(1.d0)
    u1 = low_limit
    do while ((u1 <= low_limit))
      call random_number(u1)
    end do
    call random_number(u2)
    mag = sigma_1d*dsqrt(-2.d0*log(u1))
    x1 = mag*cos(two_pi*u2) + mu_1d
    x2 = mag*sin(two_pi*u2) + mu_1d
  end subroutine sample_1D_gaussian_scalar_box_muller



  subroutine sample_1D_gaussian_vector_box_muller(dim_x, x1, mu_1d, sigma_1d)
    ! return a vector x1(dim_x) that sample a 1D Gaussian N(mu_1d, sigma_1d)
    implicit none
    integer, intent(in)  :: dim_x
    real(kind_double), intent(in)      :: mu_1d, sigma_1d
    real(kind_double), intent(out)     :: x1(:)
    ! local variables
    real(kind_double)      :: u1(dim_x), u2(dim_x), low_limit, mag(dim_x), x2(dim_x)
    real(kind_double), parameter :: two_pi = 6.283185307179586d0          ! numpy value
  
    low_limit = 2.d0*epsilon(1.d0)
    u1= low_limit
    do while ((minval(u1) <= low_limit))
      call random_number(u1)
    end do
    call random_number(u2)
    mag(:) = sigma_1d*dsqrt(-2.d0*log(u1(:)))
    x1(:) = mag(:)*cos(two_pi*u2(:)) + mu_1d
    x2(:) = mag(:)*sin(two_pi*u2(:)) + mu_1d
  end subroutine sample_1D_gaussian_vector_box_muller


  subroutine sample_dD_gaussian_vector_box_muller(no_sample_x, dim_d, x_mat, mu, sigma)
  ! input: mu(dim_d), sigma(dim_d, dim_d): mean and covariance matrix, respectively of the Gaussian distribution
  ! input: no_sample_x - number of x samples put into the matrix x_mat(dim_d, no_sample_x)
  ! output: matrix x_mat(dim_d, no_sample_x) with no_sample_x samples of the Gaussian distribution
  integer, intent(in) :: no_sample_x, dim_d
  real(kind_double), dimension(:), intent(inout) :: mu(:)
  real(kind_double), dimension(:,:), intent(inout) ::  sigma(:, :)
  ! dimension of the matrix x_mat(dim_d, no_sample_x)
  real(kind_double), dimension(:,:), allocatable, intent(out) :: x_mat
  real(kind_double), parameter :: two_pi = 6.283185307179586d0          ! numpy value

  integer :: info, ii, jj
  real(kind_double), dimension(:,:), allocatable :: ll_mat, sigma_local
  real(kind_double), dimension(dim_d) :: zz_temp 
  real(kind_double) :: low_limit, u1, u2 

  _NAMECURRENT_("sample_dD_gaussian_vector_box_muller")
  _MLD_BEGIN_

  allocate(sigma_local(dim_d, dim_d))
  allocate(x_mat(dim_d, no_sample_x))
  sigma_local(:,:) = sigma(:,:)
  ! Cholesky decomposition
  call  dpotrf('L', dim_d, sigma_local, dim_d, info)
  if (info /= 0) then
    call log_error("ML error: the Cholesky decomposition failed in "//NAMECURRENT) 
    stop "STOP: Cholesky decomposition failed!"
  end if 

  ! Copy lower triangular part to L
  allocate(ll_mat(dim_d, dim_d))
  ll_mat(:,:) = 0.0d0
  do jj = 1, dim_d
    do ii = jj, dim_d
      ll_mat(ii,jj) = sigma_local(ii,jj)
    end do
  end do



  do jj =1, no_sample_x
    ! Sample from standard normal distribution using Box-Muller transform
    do ii = 1, dim_d
      low_limit = 2.d0*epsilon(1.d0)
      u1= low_limit
      do while (u1 <= low_limit)
        call random_number(u1)
      end do
      call random_number(u2)
      zz_temp(ii) = dsqrt(-2.0d0 * log(u1)) * cos(two_pi*u2)
    end do
    ! compute the sample of x with right mean and covariance

    call dtrtrs('Lower', 'Transpose', 'Non-unit', dim_d, 1, ll_mat, dim_d, zz_temp, dim_d, info)
    if (info /= 0) then
      call log_error("ML error: solving Fourier point failed  "//NAMECURRENT) 
  end if 
    !call dtrtrs('Lower', 'No transpose', 'Non-unit', D, 1, L, D, z, D, info)
    x_mat(:,jj) = zz_temp(:) + mu(:)

    !x_mat(:,jj) = matmul(ll_mat, zz_temp) + mu(:) 
    !call dgemv('N', dim_d, dim_d, 1.0d0, sigma_in, dim_d, zz_temp, 1, 1.d0, mu_temp, 1)
    ! x_mat(:,jj) = mu_temp(:)
  end do 

  call log_info("ML info: the random sampling of the "//vtoa(dim_d)//"-dimensional Gaussian with "//vtoa(no_sample_x)//" samples is completed")
  
  deallocate(sigma_local, ll_mat)
  
  _MLD_END_
  end subroutine sample_dD_gaussian_vector_box_muller

  subroutine sample_1D_laplace_scalar(mean, width, number)
    use mld_logger
    implicit none 
    real(kind_double), intent(in) ::  mean, width
    real(kind_double), intent(out)  ::  number
    real(kind_double) :: uu 
    _NAMECURRENT_("sample_1D_laplace_scalar")
  
    _MLD_BEGIN_
  
    if (width <= 0.0d0) then
      call log_warning ("ML error: the width of Laplace distribution can be only positive "//NAMECURRENT) 
   end if
   
   call random_number(uu)
   
   if (uu < 0.5d0) then
      number = mean + width*log(2.d0*uu)
   else
      number = mean - width*log(2.d0*(1.0d0-uu))
   end if
   _MLD_END_ 
  
  end subroutine sample_1D_laplace_scalar


  subroutine sample_1D_laplace_vector(mean, width, number)
    use mld_logger
    implicit none 
    real(kind_double), intent(in) ::  mean, width
    real(kind_double), intent(inout), dimension(:)  ::  number
    real(kind_double),dimension(:), allocatable :: uu 
    integer :: ii 
    _NAMECURRENT_("sample_1D_laplace_vector")
  
    _MLD_BEGIN_
  
    if (width <= 0.0d0) then
      call log_warning ("ML error: the width of Laplace distribution can be only positive "//NAMECURRENT) 
   end if
   
   allocate(uu(size(number,1)))
  
   call random_number(uu)
   
   do ii = 1, size(number,1)
     if (uu(ii) < 0.5d0) then
        number(ii) = mean + width*log(2.d0*uu(ii))
     else
        number(ii) = mean - width*log(2.d0*(1.0d0-uu(ii)))
     end if
   end do 
  
   deallocate(uu)
  
   _MLD_END_ 
  
  end subroutine sample_1D_laplace_vector


  subroutine sample_1D_cauchy_scalar(mean, width, number)
    use mld_logger
    use ml_in_ndm_module, only: one_pi
    implicit none 
    real(kind_double), intent(in) ::  mean, width
    real(kind_double), intent(out)  ::  number
    real(kind_double) :: uu 
    _NAMECURRENT_("sample_1D_cauchy_scalar")
  
    _MLD_BEGIN_
  
    if (width <= 0.0d0) then
      call log_warning ("ML error: the width of Cauchy distribution can be only positive "//NAMECURRENT) 
    end if
  
    call random_number(uu)
  
    number = mean + width * tan (one_pi*(uu - 0.5d0))
  
    _MLD_END_ 
  
  end subroutine sample_1D_cauchy_scalar

  subroutine sample_1D_cauchy_vector(mean, width, number)
    use mld_logger
    use ml_in_ndm_module, only: one_pi
    implicit none 
    real(kind_double), intent(in) ::  mean, width
    real(kind_double), dimension(:), intent(inout)  ::  number
    real(kind_double), dimension(:), allocatable  :: uu 
    _NAMECURRENT_("sample_1D_cauchy_vector")
  
    _MLD_BEGIN_
  
    if (width <= 0.0d0) then
      call log_warning ("ML error: the width of Cauchy distribution can be only positive "//NAMECURRENT) 
    end if
  
    allocate(uu(size(number,1)) )

    call random_number(uu)
  
    number(:) = mean + width * tan (one_pi*(uu(:) - 0.5d0))
  
    deallocate(uu)
  
    _MLD_END_ 

  end subroutine sample_1D_cauchy_vector
  



FUNCTION matdet(A) result(det)

    implicit none

    REAL(kind(0.d0)), dimension(:, :), intent(in)      :: A
    REAL(kind(0.d0))     :: det

    IF ((size(A, 1) .NE. 3) .AND. (size(A, 2) .NE. 3)) &
      STOP '< MatDet >: size of matrix to invert should be equal to 3'

    det = a(1, 1)*a(2, 2)*a(3, 3) + a(1, 2)*a(2, 3)*a(3, 1) &
          + a(1, 3)*a(2, 1)*a(3, 2) - a(1, 3)*a(2, 2)*a(3, 1) &
          - a(1, 1)*a(2, 3)*a(3, 2) - a(1, 2)*a(2, 1)*a(3, 3)

  END FUNCTION matdet

  function my_exp(x) result(y) 
    implicit none 
    real(kind_double), intent(in) :: x
    real(kind_double) :: y 
    real(kind_double) :: lne, lnh 
    
    lne = log(tiny(x))
    lnh = log(huge(x))
    if (x < 0.d0  ) then 

      if ( x < 10+lne ) then 
         y = 0.d0  
      else 
         y = exp(x)
      end if 
    end if 

    if (x >= 0.d0 ) then 

      if ( x > lnh ) then 
          y = huge(x)    
      else 
          y = exp(x)
      end if      
    end if 

    ! debug write(*,*) 'xxx', x,lne, y 


  end function my_exp 

  SUBROUTINE matinv_gen(A, B)

    implicit none

    REAL(kind(0.d0)), dimension(:, :), intent(in)      :: A
    REAL(kind(0.d0)), dimension(:, :), intent(out)     :: B

    !REAL(kind(0.d0)) :: matdet
    REAL(kind(0.d0))     :: invdet

    ! Variables used with Lapack subroutine
    INTEGER  :: i, n, INFO
    INTEGER, dimension(:), allocatable     :: IPIV
    REAL(kind(0.d0)), dimension(:, :), allocatable     :: tempA, invA
    if (size(A, 1) .NE. size(A, 2)) &
      STOP '< MatInv >: matrix to invert is not square'

    SELECT CASE (size(A, 1))
    CASE (1)
      B(1, 1) = 1.d0/A(1, 1)

    CASE (2)
      invdet = 1.d0/(A(1, 1)*A(2, 2) - A(1, 2)*A(2, 1))
      B(1, 1) = A(2, 2)*invdet
      B(2, 2) = A(1, 1)*invdet
      B(1, 2) = -A(2, 1)*invdet
      B(2, 1) = -A(1, 2)*invdet

    CASE (3)
      invdet = 1.d0/matdet(A)

      b(1, 1) = a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)
      b(2, 1) = a(2, 3)*a(3, 1) - a(2, 1)*a(3, 3)
      b(3, 1) = a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1)

      b(1, 2) = a(3, 2)*a(1, 3) - a(3, 3)*a(1, 2)
      b(2, 2) = a(3, 3)*a(1, 1) - a(3, 1)*a(1, 3)
      b(3, 2) = a(3, 1)*a(1, 2) - a(3, 2)*a(1, 1)

      b(1, 3) = a(1, 2)*a(2, 3) - a(1, 3)*a(2, 2)
      b(2, 3) = a(1, 3)*a(2, 1) - a(1, 1)*a(2, 3)
      b(3, 3) = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)

      b(1:3, 1:3) = b(1:3, 1:3)*invdet

    CASE DEFAULT
      ! STOP '< MatInv >: size of matrix to invert should be less than 3 (lapack not implemented)'
      ! Use Lapack subroutine DGESV
      n = size(A, 1)
      allocate (tempA(n, n))
      tempA(1:n, 1:n) = A(1:n, 1:n)                    ! Needed in order to not modify A
      allocate (IPIV(1:n))
      allocate (invA(n, n))
      invA(1:n, 1:n) = 0.d0
      do i = 1, n
        invA(i, i) = 1.d0
      end do
      call DGESV(n, n, tempA, n, IPIV, invA, n, INFO)
      if (INFO .NE. 0) then
        write (0, *) 'Result of DGESV subroutine: INFO=', INFO
        write (0, *) ' < 0: if INFO = -i, the i-th argument had an illegal value'
        write (0, *) ' > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization'
        write (0, *) 'has been completed, but the factor U is exactly'
        write (0, *) 'singular, so the solution could not be computed.'
        STOP '< MatInv >'
      end if
      B(1:n, 1:n) = invA(1:n, 1:n)
      deallocate (tempA, invA, IPIV)

    END SELECT

  END SUBROUTINE matinv_gen

  subroutine hermite(n, x, y, d_y)

    implicit none 
    integer, intent(in)  :: n
    real(kind(0.d0)), intent(in)     :: x
    real(kind(0.d0)), intent(out)    :: y, d_y

    select case (n)

    case (1)

      y = x
      d_y = 1.d0

    case (2)

      y = x**2 - 1.d0
      d_y = 2.d0*x

    case (3)

      y = x**3 - 3.d0*x
      d_y = 3.d0*x**2 - 3.d0

    case (4)

      y = x**4 - 6.d0*x**2 + 3.d0
      d_y = 4.d0*x**3 - 12.d0*x

    case default

      if (rangml == 0) write (6, *) 'NO IMPLEMENTATION FOR THIS ORDER OF HERMITE-POLYNOMICAL-CHAOS  CASE, n is  :', n
      stop 'in hermite. oder of hermite poly'

    end select

  end subroutine hermite

  subroutine dreal_matmul(A, nrow_A, ncol_A, B, nrow_B, ncol_B, C, nrow_C, ncol_C)
    ! perform C = A*B using dgemm ....
    ! On input:
    !           the matrix A  of size  nrow_A x ncol_A
    !           the matrix B  of size  nrow_B x ncol_B
    ! On output
    !           the matrix C  of size nrow_C, ncol_C
    ! Observation:
    !            ncol_A = nrow_B
    !            nrow_A = nrow_C
    !            ncol_B = nrow_C
    implicit none 
    integer, intent(in)  :: nrow_A, ncol_A, nrow_B, ncol_B, nrow_C, ncol_C
    real(kind_double), dimension(:, :), intent(in)   :: A
    real(kind_double), dimension(:, :), intent(in)   :: B
    real(kind_double), dimension(:, :), intent(out)  :: C

    if ((ncol_A /= nrow_B) .or. (ncol_B /= ncol_C .or. (nrow_A /= nrow_C))) then
      write (6, *) "ML: multiplication wrong. stop"
      write (6, *) "ML: nrow_A, ncol_A", nrow_A, ncol_A
      write (6, *) "ML: nrow_B, ncol_B", nrow_B, ncol_B
      write (6, *) "ML: nrow_C, ncol_C", nrow_C, ncol_C

      stop "dreal_matmul wrong dimensions"
    end if

    call dgemm('N', 'N', nrow_A, ncol_B, ncol_A, 1.d0, A, nrow_A, B, nrow_B, 0.d0, C, nrow_C)

  end subroutine dreal_matmul

  subroutine dreal_inverse(dim_mat, mat, mat_inv)

    use ml_in_ndm_module, only: rangml
    implicit none 
    integer, intent(in)  :: dim_mat
    real(kind_double), dimension(dim_mat, dim_mat), intent(in) :: mat
    real(kind_double), dimension(dim_mat, dim_mat), intent(out)      :: mat_inv
    ! driver inversion
    integer  :: dim_work, dim_ipiv, info
    real(kind_double), dimension(:), allocatable   :: work
    integer, dimension(:), allocatable     :: ipiv
    real(kind_double), dimension(dim_mat, dim_mat) :: id_matrix, c_matrix
    integer  :: i


    id_matrix(:, :) = 0.d0
    do i = 1, dim_mat
      id_matrix(i, i) = 1.d0
    end do

    mat_inv(:, :) = mat(:, :)
    dim_work = max(1, 2*dim_mat)
    dim_ipiv = max(1, dim_mat)
    info = 0
    if (allocated(work)) deallocate (work); allocate (work(dim_work))
    if (allocated(ipiv)) deallocate (ipiv); allocate (ipiv(dim_ipiv))
    call dsytrf('L', dim_mat, mat_inv, dim_mat, ipiv, work, dim_work, info)   ! factorization
    !debug if (rangml==0) write (6,'("ML: dsytrf info for factorization...:", i6)') info
    if (info /= 0) then
      if (rangml == 0) write (6, '("ML: WARNING problems factorization  dsytrf info is ...:", i6)') info
    end if
    !dim_work=max(1, 2*dim_mat)
    !if (allocated(work)) deallocate(work) ; allocate (work(dim_work))
    call dsytri('L', dim_mat, mat_inv, dim_mat, ipiv, work, info)             ! inversion
    !debug if (rangml==0) write (6,'("ML: dsytri info for inversion.......:", i6)') info
    if (info /= 0) then
      if (rangml == 0) write (6, '("ML: WARNING problems inversion dsytri info is ...:", i6)') info
    end if
    call dsymm('L', 'L', dim_mat, dim_mat, 1.d0, mat_inv, dim_mat, id_matrix, dim_mat, 0.d0, c_matrix, dim_mat)
    mat_inv = c_matrix

  end subroutine dreal_inverse

  subroutine diago_serial_real_double(matfor,  eigenvectors, eigenvalues)
    implicit none 
    real(kind_double), dimension(:, :), intent(in) :: matfor
    real(kind_double), dimension(:, :), intent(out)      :: eigenvectors
    real(kind_double), dimension(:), intent(out)      :: eigenvalues
    !diago driver
    real(kind=kind(1.d0)), dimension(:), allocatable   :: work, w
    character*1    :: uplo
    integer  :: info, lwork,nmat 

    nmat = size(matfor,1)
    eigenvectors(:, :) = matfor(:, :)
    uplo = 'u'
    info = 0
    allocate (work(1), w(nmat))
    lwork = -1
    call DSYEV('V', 'L', nmat, eigenvectors, nmat, w, work, lwork, info)
    lwork = int(work(1)+2) 
    deallocate (work)
    allocate (work(lwork))
    call DSYEV('V', 'L', nmat, eigenvectors, nmat, w, work, lwork, info)

    if (info /= 0) then
      write (*, *) 'WARNING there are problems in the diago with eigenvalues and eigenvectors'
    end if

    eigenvalues(:) = w(:)
  end subroutine diago_serial_real_double


  subroutine serial_pseudo_inverse (mat, rank_mat, mat_pinv)
    use ml_in_ndm_module, only: svd_rcond, rangml 
    use mld_logger
    implicit none 
    real(kind_double), dimension(:,:), intent(inout) :: mat
    real(kind_double), dimension(:,:), allocatable, intent(inout) :: mat_pinv
    integer :: rank_mat 
    ! declaration for SVD driver ... 
    integer :: Mline, Ncolm
    integer :: lwork, LDA, LDU, LDVT,  info, ii
    integer :: rsize, icount
    real(kind_double), dimension(:), allocatable   :: work, s_svd
    real(kind_double), dimension(:,:), allocatable   :: matUU, matVT, matUUk, matVTk, mat_copy
!
    real(kind_double), parameter :: one =1.d0, zero = 0.d0 
    real(kind_double), external :: dlange
!
    _NAMECURRENT_("serial_pseudo_inverse")

    _MLD_BEGIN_
    Mline = size(mat,1)
    Ncolm = size(mat,2)

    LDA = max(Mline, 1)
    LDU = Mline 
    LDVT = min(Mline, Ncolm)
    rsize = min(Mline, Ncolm)
      
    if (allocated(s_svd)) deallocate (s_svd); allocate (s_svd(rsize))
    if (allocated(matUU)) deallocate (matUU); allocate (matUU(Mline, rsize))
    if (allocated(matVT)) deallocate (matVT); allocate (matVT(rsize, Ncolm))
    if (allocated(mat_pinv)) deallocate (mat_pinv); allocate (mat_pinv(Ncolm, Mline))
    if (allocated(mat_copy)) deallocate (mat_copy); allocate (mat_copy(Mline, Ncolm))
    mat_copy(:,:) = mat(:,:)
    lwork = -1
    if (allocated(work)) deallocate (work); allocate (work(1))
    s_svd(:) = 0
    !SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
    call dgesvd('S', 'S', Mline, Ncolm, mat_copy, LDA, s_svd, matUU,  LDU, matVT, LDVT, work, lwork, info)
    if (info < 0 ) then
      call log_info("ML error: the work array was not allocated properly in SVD in"//NAMECURRENT) 
    end if 
    lwork = int(work(1)) + 2
    if (allocated(work)) deallocate (work); allocate (work(lwork))
    call dgesvd('S', 'S', Mline, Ncolm, mat_copy, LDA, s_svd, matUU,  LDU, matVT, LDVT, work, lwork, info)
    if (info < 0 ) then
      call log_info("ML error: the serial SVD decomposition failed in "//NAMECURRENT) 
    end if 
    deallocate(mat_copy) 

    if (svd_rcond < 0) then
      svd_rcond = 100.d0*epsilon(1.d0)
    end if

    icount = 0 
    do ii = 1, rsize
      if (s_svd(ii) >= svd_rcond) then
        icount = icount + 1 
      end if  
    end do   
    rank_mat = icount 

    if (rangml == 0) then
      write (6, '("ML: dgesvd SVD decomposition info Mline Ncolm lwork RCOND RANK ....:", i4, i8, i8, i7, d20.10,i6)') &
        info, Mline, Ncolm, lwork, svd_rcond, rank_mat
    end if

    allocate(matUUk(Mline, rank_mat))
    matUUk(:,:) = matUU(:,1:rank_mat)
    do ii = 1, rank_mat
      matUUk(:,ii) = matUU(:,ii) / s_svd(ii)
    end do 
    deallocate(matUU)

    allocate(matVTk(rank_mat, Ncolm))
    matVTk(:,:) = matVT(1:rank_mat, :)
    deallocate(matVT)
    
    !  dgemm (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
    ! C(mxn) = A(mxk)*B(kxn)
    ! A+(nxm) = op(VT)(nxr)* op(U)(rxM)
    call dgemm('T', 'T', Ncolm, Mline, rank_mat, one, matVTk, rank_mat, matUUk, Mline, zero, mat_pinv, Ncolm )
    !mat_pinv = matmul(transpose(matVTk), transpose(matUUk) )
    deallocate (s_svd, work, matUUk, matVTk)

    !!----------------------test pseudo inverse -------------------- 
    !allocate(mtmp(Ncolm, Ncolm))
    !!mtmp(nxn)=pinvA(nxm)*A(mxn)
    !call dgemm('N', 'N', Ncolm, Ncolm, Mline, one, mat_pinv, size(mat_pinv,1), mat, size(mat,1), zero, mtmp, size(mtmp,1))
    !!mtmp = matmul(mat_pinv, mat) 
    !allocate(mtest(Mline, Ncolm))
    !!mtest(mxn) = a(mxn)*mtmp(nxn)
    !call dgemm('N', 'N', Mline, Ncolm, Ncolm,one, mat, size(mat,1), mtmp, size(mtmp,1), zero, mtest, size(mtest,1))
    !!mtest = matmul(mat, mtmp) 
    !deallocate(mtmp)
    !allocate(mdiff(Mline, Ncolm))
    !mdiff(:,:) = mtest(:,:) - mat(:,:)
    !deallocate(mtest) 
    !if (allocated(work)) deallocate(work)
    !allocate(work(Mline))
    !norm = dlange('F', Mline, Ncolm, mdiff, size(mdiff,1), work )
    !deallocate(work, mdiff)
    !write(*,*) 'test difference pinv_mat - mat', norm 
    !norm = dlange('F', Mline, Ncolm, mat, size(mat,1), work )
    !write(*,*) 'initial norm of matrix', norm
    !!---------------------end test pseudo inverse ------------------ 

    _MLD_END_
  end subroutine serial_pseudo_inverse 

  subroutine serial_svd_filter (mat, rank_mat, mat_filter, determinant, determinant_filter)
    use ml_in_ndm_module, only: rangml 
    use mld_logger
    implicit none 
    real(kind_double), dimension(:,:), intent(inout) :: mat
    real(kind_double), dimension(:,:), allocatable, intent(inout) :: mat_filter
    real(kind_double), intent(out) :: determinant, determinant_filter
    integer :: rank_mat 
    ! declaration for SVD driver ... 
    integer :: Mline, Ncolm
    integer :: lwork, LDA, LDU, LDVT,  info, ii
    integer :: rsize, icount
    real(kind_double), dimension(:), allocatable   :: work, s_svd
    real(kind_double), dimension(:,:), allocatable   :: matUU, matVT, matUUk, matVTk, mat_copy
    
    real(kind_double), parameter :: one =1.d0, zero = 0.d0 
    real(kind_double), external :: dlange
!
    real(kind_double) :: filter_threshold
    integer :: fsvd 
    _NAMECURRENT_("serial_pseudo_inverse")

    _MLD_BEGIN_
    Mline = size(mat,1)
    Ncolm = size(mat,2)

    LDA = max(Mline, 1)
    LDU = Mline 
    LDVT = min(Mline, Ncolm)
    rsize = min(Mline, Ncolm)
      
    if (allocated(s_svd)) deallocate (s_svd); allocate (s_svd(rsize))
    if (allocated(matUU)) deallocate (matUU); allocate (matUU(Mline, rsize))
    if (allocated(matVT)) deallocate (matVT); allocate (matVT(rsize, Ncolm))
    if (allocated(mat_filter)) deallocate (mat_filter); allocate (mat_filter(Ncolm, Mline))
    if (allocated(mat_copy)) deallocate (mat_copy); allocate (mat_copy(Mline, Ncolm))
    mat_copy(:,:) = mat(:,:)
    lwork = -1
    if (allocated(work)) deallocate (work); allocate (work(1))
    s_svd(:) = 0
    !SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
    call dgesvd('S', 'S', Mline, Ncolm, mat_copy, LDA, s_svd, matUU,  LDU, matVT, LDVT, work, lwork, info)
    if (info < 0 ) then
      call log_info("ML error: the work array was not allocated properly in SVD in"//NAMECURRENT) 
    end if 
    lwork = int(work(1)) + 2
    if (allocated(work)) deallocate (work); allocate (work(lwork))
    call dgesvd('S', 'S', Mline, Ncolm, mat_copy, LDA, s_svd, matUU,  LDU, matVT, LDVT, work, lwork, info)
    if (info < 0 ) then
      call log_info("ML error: the serial SVD decomposition failed in "//NAMECURRENT) 
    end if 
    deallocate(mat_copy) 

    
    !!! filter_threshold = 1.d5*epsilon(1.d0)
    filter_threshold = 1.d-2
    

    icount = 0 
    
    determinant = 1.d0
    determinant_filter = 1.d0

   if (rangml==0) open(newunit=fsvd, file='svd_eigenvalues.dat', status='unknown') 
    do ii = 1, rsize
      if (rangml == 0) then  
        write (fsvd, '(es20.10)') s_svd(ii)
      end if   
      determinant = determinant * s_svd(ii)
      if (s_svd(ii) >= filter_threshold) then
        determinant_filter = determinant_filter * s_svd(ii)
        icount = icount + 1 
      end if  
    end do   
    if (rangml==0) close(fsvd)
    rank_mat = icount 

    if (rangml == 0) then
      write (6, '("ML: matrix info  det, det_filter, RCOND RANK ....:", i4, 3d20.10,i6)') &
        info, determinant, determinant_filter, filter_threshold, rank_mat
    end if

    allocate(matUUk(Mline, rank_mat))
    matUUk(:,:) = matUU(:,1:rank_mat)
    do ii = 1, rank_mat
      matUUk(:,ii) = matUU(:,ii) * s_svd(ii)
    end do 
    deallocate(matUU)

    allocate(matVTk(rank_mat, Ncolm))
    matVTk(:,:) = matVT(1:rank_mat, :)
    deallocate(matVT)
    
    !  dgemm (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
    ! C(mxn) = A(mxk)*B(kxn)
    ! A+(nxm) = op(VT)(nxr)* op(U)(rxM)
    call dgemm('T', 'T', Ncolm, Mline, rank_mat, one, matVTk, rank_mat, matUUk, Mline, zero, mat_filter, Ncolm )

    deallocate (s_svd, work, matUUk, matVTk)



    _MLD_END_
  end subroutine serial_svd_filter 



  subroutine serial_determinant_symmetric_general (mat, det)
    use, intrinsic :: iso_fortran_env, dp=>real64
    use ml_in_ndm_module, only: rangml
    use mld_logger
    implicit none 
    real(kind_double), dimension(:,:), intent(in) :: mat
    real(kind_double), intent(out) :: det
    ! driver inversion
    integer  :: dim_work, dim_ipiv, info
    real(kind_double), dimension(:), allocatable   :: work
    integer, dimension(:), allocatable     :: ipiv
    real(kind_double), dimension(:,:), allocatable   :: mat_copy
    integer  :: i
    _NAMECURRENT_("serial_determinant_symmetric")

    _MLD_BEGIN_
    if (size(mat,1) /= size(mat,2)) then
      call log_warning("ML error: the matrix is not square in "//NAMECURRENT) 
    end if 
    if (allocated(mat_copy)) deallocate (mat_copy); allocate (mat_copy(size(mat,1), size(mat,2)))
    mat_copy(:,:) = mat(:,:)
    dim_work = max(1, 2*size(mat,1))
    dim_ipiv = max(1, size(mat,1))
    info = 0
    if (allocated(work)) deallocate (work); allocate (work(dim_work))
    if (allocated(ipiv)) deallocate (ipiv); allocate (ipiv(dim_ipiv))
    call dsytrf('L', size(mat,1), mat_copy, size(mat,1), ipiv, work, dim_work, info)   ! factorization
    if (info /= 0) then
      if (rangml == 0) write (6, '("ML: WARNING problems factorization  dsytrf info is ...:", i6)') info
    end if
    det = 1.0_dp
    do  while (i <= size(mat,1)) 
      if (ipiv(i) > 0 ) then 
        det = det * mat_copy(i,i)
        i =i + 1 
      else if (ipiv(i) < 0 .and. i < size(mat,1)) then 
        det = det * mat_copy(i,i) * mat_copy(i+1,i+1) - mat_copy(i,i+1)**2
        i =i + 2 ! Skip next diagonal element
      else 
        ! handle error 
        info = -1 
        exit 
      end if        
    end do
    if (info < 0) then
      call log_critical("ML error: the determinant of the matrix dsytrf failed in  "//NAMECURRENT)
      stop "factorization  dsytrf failed"
    end if
    deallocate (mat_copy, work, ipiv)
    _MLD_END_
  end subroutine serial_determinant_symmetric_general

SUBROUTINE ComputeDeterminantQR(A, N, Det)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  REAL(8), DIMENSION(N,N), INTENT(IN) :: A
  REAL(8), INTENT(OUT) :: Det
  REAL(8), DIMENSION(N,N) :: QR
  REAL(8), DIMENSION(N) :: TAU
  real(kind_double), dimension(:), allocatable :: work
  INTEGER :: I, INFO, lwork 
  _NAMECURRENT_("ComputeDeterminantQR")
  ! Copy A to QR since DGEQRF overwrites its input
  QR = A
  call log_warning("ML: here we get the determinant using the QR decomposition. There is no way to have the sign. ")
  ! QR decomposition of the matrix
  IF (ALLOCATED(work)) DEALLOCATE(work)
  ALLOCATE(work(1))
  INFO = 0
  lwork = -1
  call dgeqrf(N, N, QR, N, TAU, work, lwork, info)
  lwork = int(work(1)) + 2
  DEALLOCATE(work)
  ALLOCATE(work(lwork))
  call dgeqrf(N, N, QR, N, TAU, work, lwork, info)

  IF (INFO /= 0) THEN
     call log_warning("ML error: the QR fact dgqrf failed in  "//NAMECURRENT)
     Det = 1.0D0
     RETURN
  END IF

  ! Compute the determinant as the product of the diagonal elements of R
  Det = 1.0D0
  DO I = 1, N
     Det = Det * QR(I, I)
  END DO

END SUBROUTINE ComputeDeterminantQR


  SUBROUTINE ComputeDeterminantLU(A, N, Det)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  REAL(8), DIMENSION(N,N), INTENT(IN) :: A
  REAL(8), INTENT(OUT) :: Det
  REAL(8), DIMENSION(N,N) :: LU
  INTEGER, DIMENSION(N) :: IPIV
  INTEGER :: I, INFO, SIGN

  ! Copy A to LU since DGETRF overwrites its input with the LU decomposition
  LU = A

  ! LU decomposition of the matrix, DGETRF computes an LU decomposition
  CALL DGETRF(N, N, LU, N, IPIV, INFO)

  IF (INFO /= 0) THEN
     ! Matrix is singular and its determinant is zero
     Det = 0.0D0
     RETURN
  END IF

  ! Compute the determinant as the product of the diagonal elements of LU
  Det = 1.0D0
  SIGN = 1
  DO I = 1, N
     Det = Det * LU(I, I)
     IF (IPIV(I) /= I) SIGN = -SIGN  ! Adjust for row interchanges
  END DO

  Det = SIGN * Det

END SUBROUTINE ComputeDeterminantLU


  subroutine serial_determinant_symmetric_positive (mat, det)
    use, intrinsic :: iso_fortran_env, dp=>real64
    use ml_in_ndm_module, only: rangml
    use mld_logger
    implicit none 
    real(kind_double), dimension(:,:), intent(in) :: mat
    real(kind_double), intent(out) :: det
    ! driver inversion
    integer  :: info
    real(kind_double), dimension(:,:), allocatable   :: mat_copy
    integer  :: i
    _NAMECURRENT_("serial_determinant_symmetric_positive")

    _MLD_BEGIN_
    if (size(mat,1) /= size(mat,2)) then
      call log_warning("ML error: the matrix is not square in "//NAMECURRENT) 
    end if 
    if (allocated(mat_copy)) deallocate (mat_copy); allocate (mat_copy(size(mat,1), size(mat,2)))
    mat_copy(:,:) = mat(:,:)
    info = 0
    call dpotrf('L', size(mat,1), mat_copy, size(mat,1), info)   ! factorization
    if (info /= 0) then
      if (rangml == 0) write (6, '("ML: WARNING problems factorization  dpotrf info is ...:", i6)') info
    end if
    det = 1.0_dp
    do  while (i <= size(mat,1)) 
      det = det * mat_copy(i,i)**2
      i = i + 1
    end do 
  end subroutine serial_determinant_symmetric_positive


  subroutine gen_hash_key_with_basis (nu_val, size_of_basis, ivec, hash_key) 
    ! generate the hash key of a vector for some basis: e. g. ivec = (1,2,3) and basis = 4
    ! then the hash key is 1 + 2*4 + 3*4^2 = 1 + 8 + 48 = 57
    integer, intent(in)  :: nu_val, size_of_basis, ivec(:)
    integer, intent(out) :: hash_key
    !local variable 
    integer :: ii 
    hash_key = 0
    do ii = 1, nu_val
      hash_key = hash_key + ivec(ii)* (size_of_basis)**(ii-1)  
    end do
  end subroutine gen_hash_key_with_basis

  subroutine get_max_hash_key_with_basis(nu_val, size_of_basis, nmax)
    ! generate the max hash key for a basis and a vector of size nu_val
    ! e. g. nu_val = 3 and basis = 4 then the max hash key is 4^3 + 4^2 + 4 = 64 + 16 + 4 = 84
    integer, intent(in)  :: nu_val, size_of_basis
    integer, intent(out) :: nmax
    !local variable 
    integer :: ii 
    nmax = 0 
    do ii = 1, nu_val
      nmax = nmax + size_of_basis*size_of_basis**(ii-1)  
    end do
  end subroutine get_max_hash_key_with_basis

end module math

module module_serial_linear_solver
  use module_kind_variables, only: kind_double
  use module_optimization, only: optimize_weights_db, optimize_weights_chem
  implicit none
  contains

  subroutine serial_lsystem_by_home(AA, yy, Mline, Ncolm, w_params)
    
    use mld_logger
    implicit none
    real(kind_double), intent(inout)  :: AA(:,:), yy(:,:)
    real(kind_double), intent(inout) ::  w_params(:,:)
    integer, intent(in)  :: Mline, Ncolm
    real(kind_double), dimension(:), allocatable   :: tau, work
    integer, dimension(:), allocatable :: ipiv 
    integer :: lwork, info
    character(len=80) :: cinfo 
    _NAMECURRENT_("serial_lsystem_by_home")
  
    _MLD_BEGIN_
  
    if  (Mline > Ncolm) then 
        allocate(tau(Ncolm))
      ! QR factorization of A
      lwork=-1
      if (allocated(work)) deallocate (work); allocate (work(1))
      call dgeqrf(Mline, Ncolm, AA, Mline, tau, work, lwork, info)
      lwork = int(work(1)) + 2
      if (allocated(work)) deallocate (work); allocate (work(lwork))
      call dgeqrf(Mline, Ncolm, AA, Mline, tau, work, lwork, info)
      if (info < 0) then 
        write(cinfo,'(i8)') info
        call log_warning("ML error: the QR fact dgeqrf failed in  "//NAMECURRENT//", info is "//trim((cinfo)))
      end if
      ! Multiply Q^T with y (results stored in y)
      lwork=-1
      if (allocated(work)) deallocate (work); allocate (work(1))
      call dormqr('L', 'T', Mline, 1, Ncolm, AA, Mline, tau, yy, Mline, work, lwork, info)
      lwork = int(work(1)) + 2
      if (allocated(work)) deallocate (work); allocate (work(lwork))
      call dormqr('L', 'T', Mline, 1, Ncolm, AA, Mline, tau, yy, Mline, work, lwork, info)
      if (info < 0) then 
        write(cinfo,'(i8)') info
        call log_warning("ML error: the multiplication of Q^T with y failed in  "//NAMECURRENT//", info is "//trim((cinfo)))
      end if

      ! Solve triangular system using the R factor of A and the transformed y
      call dtrtrs('U', 'N', 'N', Ncolm, 1, AA, Mline, yy, Mline, info)
      if (info < 0) then 
        write(cinfo,'(i8)') info
        call log_warning("ML error: the triangular system dtrtrs failed in  "//NAMECURRENT//", info is "//trim((cinfo)))
      end if
      w_params(:, 1) = yy(1:Ncolm, 1)
      deallocate(tau, work)
    end if 
  
    if (Mline < Ncolm) then 
      allocate(tau(Mline))
      ! LQ factorization of A
      lwork=-1
      if (allocated(work)) deallocate (work); allocate (work(1))
      call dgelqf(Mline, Ncolm, AA, Mline, tau, work, lwork, info)
      lwork = int(work(1)) + 2
      if (allocated(work)) deallocate (work); allocate (work(lwork))
      call dgelqf(Mline, Ncolm, AA, Mline, tau, work, lwork, info)
      if (info < 0) then 
        write(cinfo,'(i8)') info
        call log_warning("ML error: the LQ fact dgefq failed in  "//NAMECURRENT//", info is "//trim((cinfo)))
      end if 

      ! Solve triangular system using the L factor of A
      call dtrtrs('L', 'N', 'N', Mline, 1, AA, Mline, yy, Mline, info)
      if (info < 0) then 
        write(cinfo,'(i8)') info
        call log_warning("ML error: the triangular system dtrtrs failed in  "//NAMECURRENT//", info is "//trim((cinfo)))
      end if

      ! Multiply Q with y to get the final solution
      lwork=-1
      if (allocated(work)) deallocate (work); allocate (work(1))
      call dormlq('L', 'N', Ncolm, 1, Mline, AA, Mline, tau, yy, Mline, work, lwork, info)
      lwork = int(work(1)) + 2
      if (allocated(work)) deallocate (work); allocate (work(lwork))
      call dormlq('L', 'N', Ncolm, 1, Mline, AA, Mline, tau, yy, Mline, work, lwork, info)      
      if (info < 0) then 
        write(cinfo,'(i8)') info
        call log_warning("ML error: the multiplication of Q with y failed in  "//NAMECURRENT//", info is "//trim((cinfo)))
      end if

      deallocate(tau, work)

      w_params(:, 1) = yy(1:Ncolm, 1)

    end if 

    if (Mline == Ncolm) then 
      
      allocate(ipiv(Mline))
      
      ! LU decomposition of A
      call dgetrf(Mline, Ncolm, AA, Mline, ipiv, info)
      if (info /= 0) then 
        write(cinfo,'(i8)') info
        call log_warning("ML error: the LU fact dgetrf failed in  "//NAMECURRENT//", info is "//trim((cinfo)))
        call messages_xgetrf(info)
      end if

      ! Solve using the LU factors
      call dgetrs('N', Ncolm, 1, AA, Mline, ipiv, yy, Mline, info)
      if (info /= 0) then 
        write(cinfo,'(i8)') info
        call log_warning("ML error: the triangular system dgetrs failed in  "//NAMECURRENT//", info is "//trim((cinfo)))
        call messages_xgetrs(info)
      end if

      w_params(:, 1) = yy(1:Ncolm, 1)
      deallocate(ipiv)

    end if 

    !if (.not. (optimize_weights_db.or.optimize_weights_chem)) then 
    !if (rangml == 0) write (6, '("ML: after QR...............................")')
    !end if 

    _MLD_END_
  end subroutine serial_lsystem_by_home


  subroutine  messages_xgetrf(info)
    use mld_logger
    integer, intent(in) :: info
    _NAMECURRENT_("messages_xgetrf")
    _MLD_BEGIN_ 

    if (info < 0) then 
      call log_warning(" If info = "//vtoa(info)//" parameter "//vtoa(-info)//" had an illegal value.")
    end if 

    if (info > 0) then 
      call log_warning("If info = "//vtoa(info)//" uii is 0. The factorization has been completed, but U is exactly singular.")  
      call log_warning(" Division by 0 will occur if you use the factor U for solving a system of linear equations.")
    end if                     

    _MLD_END_
  end  subroutine  messages_xgetrf

  subroutine  messages_xgetrs(info)
    use mld_logger
    integer, intent(in) :: info
    _NAMECURRENT_("messages_xgetrf")
    _MLD_BEGIN_ 

    if (info < 0) then 
      call log_warning(" If info = "//vtoa(info)//" parameter "//vtoa(-info)//" had an illegal value.")
    end if 

    if (info > 0) then 
      call log_warning("Unknwon error for dgetrs")
    end if                     

    _MLD_END_
  end  subroutine  messages_xgetrs



  !$! Here a proof that the brute invesion is the worst method to solve a linear system:
  !$! subroutine serial_lsystem_by_home (phi, Cmat, yy, w_params)
  !$!   use math, only : dreal_inverse  
  !$!   real(kind_double), intent(in)    :: phi(:,:), Cmat(:,:), yy(:,:)
  !$!   real(kind_double), intent(inout) :: w_params(:,:)
  !$!   real(kind=kind(1.d0)), allocatable, dimension(:, :)      :: phi_inv, rhs
  !$!   
  !$!   if ( size(phi,1) /= size(Cmat,1) ) then
  !$!     write(6,'("ML error: size mismatch in serial_lsystem_by_home", 2i8)') size(phi,1), size(Cmat,1)
  !$!     stop "size mismatch in serial_lsystem_by_home" 
  !$!   end if 
  !$! 
  !$!   if (allocated(phi_inv)) deallocate (phi_inv); allocate (phi_inv(size(Cmat, 1), size(Cmat, 1)))
  !$!   call dreal_inverse(size(phi, 1), phi, phi_inv)
  !$!   if (allocated(rhs)) deallocate (rhs); allocate (rhs(size(Cmat, 1), size(yy, 2)))
  !$!   rhs = matmul(Cmat, yy)
  !$!   w_params(:, :) = matmul(phi_inv, rhs)
  !$! end subroutine serial_lsystem_by_home

  subroutine serial_lsystem_by_qr (AA, yy, Mline, Ncolm, w_params)
    use ml_in_ndm_module, only: rangml
    real(kind_double), intent(in)  :: AA(:,:), yy(:,:)
    real(kind_double), intent(inout) ::  w_params(:,:)
    integer, intent(in)  :: Mline, Ncolm
    integer :: lwork, LDA, LDB,  info
    real(kind_double), dimension(:), allocatable   :: work


    LDA = max(Mline, 1)
    LDB = max(Mline, max(Ncolm, 1))
    lwork = -1
    if (allocated(work)) deallocate (work); allocate (work(1))
    call dgels('N', Mline, Ncolm, size(yy, 2), AA, &
               LDA, yy, LDB, work, lwork, info)

    lwork = int(work(1)) + 2
    if (allocated(work)) deallocate (work); allocate (work(lwork))
    call dgels('N', Mline, Ncolm, size(yy, 2), AA, &
               LDA, yy, LDB, work, lwork, info)

    if (.not. (optimize_weights_db.or.optimize_weights_chem)) then 
    if (rangml == 0) write (6, '("ML: after QR...............................")')
    end if 

    if (info > 0) then
      write (*, *) 'The diagonal element ', INFO, ' of the triangular '
      write (*, *) 'factor of A is zero, so that A does not have full '
      write (*, *) 'rank; the least squares solution could not be '
      write (*, *) 'computed.'
      STOP
    end if
    w_params(:, 1) = yy(1:Ncolm, 1)
    deallocate (work)
  end subroutine serial_lsystem_by_qr

  subroutine serial_lsystem_by_svd (AA, yy, Mline, Ncolm, w_params, svd_rcond_local)
    use ml_in_ndm_module, only: rangml ! TODO cvw, w_params
    real(kind_double), intent(in)  :: AA(:,:), yy(:,:)
    real(kind_double), intent(inout) ::  w_params(:,:)
    integer, intent(inout)  :: Mline, Ncolm
    real(kind_double), intent(in) :: svd_rcond_local 
    integer :: lwork, liwork, LDA, LDB, RANK, info
    real(kind_double) :: RCOND
    real(kind_double), dimension(:), allocatable   :: work, s_svd
    integer, dimension(:), allocatable     :: iwork

    !svd_rcond_local = svd_rcond
    !svd_rcond_local = -1.d0 ! this is a very small value to make sure that we keep all the singular values. We can adjust it later if needed.
    RCOND = svd_rcond_local
    LDA = max(Mline, 1)
    LDB = max(Mline, max(Ncolm, 1))


    lwork = -1

    if (allocated(work)) deallocate (work); allocate (work(1))
    if (allocated(iwork)) deallocate (iwork); allocate (iwork(1))
    if (allocated(s_svd)) deallocate (s_svd); allocate (s_svd(min(Mline, Ncolm)))
    s_svd(:) = 0
    call dgelsd(Mline, Ncolm, size(yy, 2), AA, LDA, yy, LDB, s_svd, RCOND, RANK, work, lwork, iwork, info)
    lwork = int(work(1)) + 2
    liwork = int(iwork(1)) + 2
    if (allocated(work)) deallocate (work); allocate (work(lwork))
    if (allocated(iwork)) deallocate (iwork); allocate (iwork(liwork))

    RCOND = svd_rcond_local
    call dgelsd(Mline, Ncolm, size(yy, 2), AA, LDA, yy, LDB, s_svd, RCOND, RANK, work, lwork, iwork, info)
    !if (rangml==0) write (6,'("ML: SVD kernel inversion LDA, LDB, RANK", 3i8)') LDA, LDB, RANK
    if (.not. (optimize_weights_db.or.optimize_weights_chem)) then 
    if (rangml == 0) write (6, '("ML: after SVD...............................")')
    end if 
    !if (rangml==0) write (6,'("ML: SVD linear inversion LDA, LDB, RANK", 3i6)') LDA, LDB, RANK
    if (.not. (optimize_weights_db.or.optimize_weights_chem)) then 
    if (rangml == 0) then
      write (6, '("ML: dgelsd SVD inversion info Mline Ncolm lwork RCOND RANK ....:", i4, i8, i8, i7, d20.10,i6)') &
        info, Mline, Ncolm, lwork, RCOND, rank
    end if
    end if 

    w_params(:, 1) = yy(1:Ncolm, 1)
    deallocate (s_svd, work, iwork)
    !if (allocated(yy)) deallocate(yy)
    !if (allocated(AA)) deallocate (AA)

    if (info > 0) then
      if (rangml == 0) then
        write (6, '("ML: dgelsd inversion info lwork liwork RANK ....:", i4, 2i6,i5)') &
          info, lwork, liwork, rank
        write (*, *) 'The diagonal element ', INFO, ' of the triangular '
        write (*, *) 'factor of A is zero, so that A does not have full '
        write (*, *) 'rank; the least squares solution could not be '
        write (*, *) 'computed.'
        stop
      end if
    end if
  end subroutine serial_lsystem_by_svd

  subroutine serial_lsystem_by_ortho (AA, yy, Mline, Ncolm, w_params)
    use ml_in_ndm_module, only: svd_rcond, rangml
    real(kind_double), intent(in)  :: AA(:,:), yy(:,:)
    real(kind_double), intent(inout) ::  w_params(:,:)
    integer, intent(in)  :: Mline, Ncolm
    integer :: lwork, LDA, LDB, RANK, info
    real(kind_double) :: RCOND
    real(kind_double), dimension(:), allocatable   :: work
    integer, dimension(:), allocatable     :: JPVT

    RCOND = svd_rcond
    LDA = max(Mline, 1)
    LDB = max(Mline, max(Ncolm, 1))

    RCOND = svd_rcond
    lwork = -1
    if (allocated(work)) deallocate (work); allocate (work(1))
    if (allocated(JPVT)) deallocate (JPVT); allocate (JPVT(Ncolm))
    JPVT(:) = 0
    call dgelsy(Mline, Ncolm, size(yy, 2), AA, LDA, yy, LDB, JPVT, RCOND, RANK, work, lwork, info)
    lwork = int(work(1)) + 2
    if (allocated(work)) deallocate (work); allocate (work(lwork))
    RCOND = svd_rcond
    call dgelsy(Mline, Ncolm, size(yy, 2), AA, LDA, yy, LDB, JPVT, RCOND, RANK, work, lwork, info)

    w_params(:, 1) = yy(1:Ncolm, 1)
    deallocate (JPVT, work)
    if (.not. (optimize_weights_db.or.optimize_weights_chem)) then 
    if (rangml == 0) then
      write (6, '("ML: dgelsy inversion info Mline Ncolm lwork RCOND RANK ....:", i4, i8, i8, i7, d20.10,i6)') &
        info, Mline, Ncolm, lwork, RCOND, rank
    end if
    end if 

    if (info > 0) then
      write (*, *) 'The diagonal element ', INFO, ' of the triangular '
      write (*, *) 'factor of A is zero, so that A does not have full '
      write (*, *) 'rank; the least squares solution could not be '
      write (*, *) 'computed.'
      stop
    end if
  end subroutine serial_lsystem_by_ortho



  subroutine serial_lsystem_by_constraints (Cmat, Bmat,ymat, zmat, w_params)
    use ml_in_ndm_module, only: rangml
    ! implicit none
    real(kind_double), intent(in) :: Cmat(:,:), Bmat(:,:), w_params(:,:), ymat(:,:), zmat(:,:)
    integer ::  lwork, info
    real(kind_double), dimension(:), allocatable   :: work
    !call dgglse(m, n, p, A = Amat^T, lda, B=Bmat^T, ldb, c=ymat, d=zmat, x, work, lwork, info)
    !integer m: The number of rows of the matrix A (m ≥ 0).
    !integer n: The number of columns of the matrices A and B (n ≥ 0).
    !integer p:  The number of rows of the matrix B (0 ≤ p ≤ n ≤ m+p).
    !integer lda: max(1,m)
    !integer ldb: max(1,p)
    !lwork= max(1,m+n+p), if lwork=-1 this is calculated automatically
    lwork = max(1, size(Cmat, 2) + size(Cmat, 1) + size(Bmat, 2) + size(Bmat, 1))
    if (allocated(work)) deallocate (work); allocate (work(lwork))
    if (size(Bmat, 2) > size(Cmat, 1)) then
      if (rangml == 0) then
        write (6, *) 'The constraints problem is ill posed.'
        write (6, *) 'The number of constraints should be lower than the dimesion of the descriptor.'
        write (6, *) 'xdesc dimesion ', size(Cmat, 1), 'number of constraints', size(Bmat, 2)
      end if
    end if
    call dgglse(size(Cmat, 2), size(Cmat, 1), size(Bmat, 2), transpose(Cmat), max(1, size(Cmat, 2)), &
                transpose(Bmat), max(1, size(Bmat, 2)), ymat(:, 1), zmat(:, 1), w_params(:, 1), work, lwork, info)
    write (*, *) 'LSE + Constraints with info', info
  end subroutine serial_lsystem_by_constraints

end module module_serial_linear_solver




module module_sample_rand_ker
  use module_kind_variables, only: kind_double

  contains
  
  subroutine get_rff_sigma(rff_type,  mu, sigma, dim_desc, dim_KF, omegas, bcoeffs)
    !use mpi
    use mld_mpi, only: comm_mld, mld_critical_abort

    ! ---------> sample p(x) ~ exp(-0.5*xx^T sigma xx)
    ! rff_type : type of distribution 1 Gaussian ; 2 Laplace, 3 Cauchy 
    ! mu       : \mu of distribution 
    ! sigma    : defines a value for sigma matrix 
    ! omegas   : the random omega's: D x K (or F)  size
    ! bcoeffs  : vectors of b between 0:2 pi, of dimension K (or F) 
    use ml_in_ndm_module, only: rangml, two_pi
    use module_kernel, only: krff_gaussian, krff_cauchy, krff_laplace
    use math, only: sample_dD_gaussian_vector_box_muller
    use mld_logger
    implicit none 
    integer, intent(in) :: rff_type, dim_desc, dim_KF
    real(kind_double), dimension(:), intent(inout) :: mu
    real(kind_double), dimension(:,:), intent(inout) ::  sigma

    real(kind_double), dimension(:,:), allocatable, intent(inout) :: omegas 
    real(kind_double), dimension(:), allocatable, intent(inout) :: bcoeffs 
!
    _NAMECURRENT_("get_rff_sigma")

    _MLD_BEGIN_

    !$! if (dabs(mu) > 1.d-10) then 
    !$!   call log_warning("ML: the median of Random Feature distribution can be only 0. Now in "//NAMECURRENT//" is "// vtoa(mu))
    !$!   call log_warning("ML: it will be assigned to 0. ")
    !$! end if 
    if (rangml == 0) then
      select case (rff_type)
        case (krff_gaussian)
          call sample_dD_gaussian_vector_box_muller(dim_KF, dim_desc, omegas, mu, sigma)
        case (krff_laplace)
          call log_error("ML error: this choice of krff_type Laplace is not possible: " // vtoa(rff_type))
          call mld_critical_abort("ML error: random Fourier Feature in  " // NAMECURRENT)   
        case (krff_cauchy)
          call log_error("ML error: this choice of krff_type Cauchy is not possible: " // vtoa(rff_type))
          call mld_critical_abort("ML error: random Fourier Feature in  " // NAMECURRENT)   
        case default 
          call log_error("ML error: the choice of krff_type: " // vtoa(rff_type))
          call mld_critical_abort("ML error: there is no choice for this random Fourier Feature in " // NAMECURRENT)   
      end select
      omegas(:, :) = omegas(:,:)
    end if
    !TORC! call MPI_BCAST(omegas, size(omegas, 1)*size(omegas,2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, codeml)
    call comm_mld%bcast(0, omegas)

  ! the phase of random kernel:
  if (allocated(bcoeffs)) deallocate (bcoeffs); allocate (bcoeffs(dim_KF))
  if (rangml==0) then 
    call random_number(bcoeffs)
    bcoeffs(:) = bcoeffs(:)*two_pi
  end if 
  !TORC! call MPI_BCAST(bcoeffs, size(bcoeffs,1), MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, codeml) 
  call comm_mld%bcast(0, bcoeffs)

  _MLD_END_

end subroutine get_rff_sigma  

  subroutine get_rff_val_sigma(rff_type, mu, dim_desc, dim_KF, omegas, bcoeffs)

    !use mpi
    use mld_mpi, only: comm_mld, mld_critical_abort

    ! --------->
    ! rff_type : type of distribution 1 Gaussian ; 2 Laplace, 3 Cauchy 
    ! mu       : \mu of distribution 
    ! sigma    : defines a value for scalar sigma
    ! omegas   : the random omega's: D x K (or F)  size
    ! bcoeffs  : vectors of b between 0:2 pi, of dimension K (or F) 
    use ml_in_ndm_module, only: rangml, two_pi
    use module_kernel, only: krff_gaussian, krff_cauchy, krff_laplace
    use math, only: sample_1D_gaussian_vector_box_muller, sample_1D_laplace_vector, sample_1D_cauchy_vector
    use mld_logger
    implicit none 
    integer, intent(in) :: rff_type, dim_desc, dim_KF
    real(kind_double), intent(in) :: mu
    real(kind_double), dimension(:,:), allocatable, intent(inout) :: omegas 
    real(kind_double), dimension(:), allocatable, intent(inout) :: bcoeffs 
    real(kind_double), dimension(:), allocatable :: tmp_kernel 
    real(kind_double) :: tmp_mu, tmp_sigma 
    integer :: ik 
    _NAMECURRENT_("get_rff_val_sigma")

    _MLD_BEGIN_

    if (allocated(tmp_kernel)) deallocate (tmp_kernel); allocate (tmp_kernel(dim_desc))
    if (dabs(mu) > 1.d-10) then 
      call log_warning("ML: the median of Random Feature distribution can be only 0. Now in "//NAMECURRENT//" is "// vtoa(mu))
      call log_warning("ML: it will be assigned to 0. ")
    end if 
    tmp_mu = 0.d0   
    tmp_sigma = 1.d0
    do ik = 1, dim_KF
      if (rangml == 0) then
        select case (rff_type)
          case (krff_gaussian)
            call sample_1D_gaussian_vector_box_muller(dim_desc, tmp_kernel, tmp_mu, tmp_sigma)
          case (krff_laplace)
            call sample_1D_laplace_vector(tmp_mu, tmp_sigma, tmp_kernel)
          case (krff_cauchy)
            call sample_1D_cauchy_vector(tmp_mu, tmp_sigma, tmp_kernel)
          case default 
            call log_error("ML error: the choise of krff_type: " // vtoa(rff_type))
            call mld_critical_abort("ML error: there is no choice for this random Fourier Feature in " // NAMECURRENT)   
        end select   
      end if
      !TORC! call MPI_BCAST(tmp_kernel, size(tmp_kernel, 1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, codeml)
      call comm_mld%bcast(0, tmp_kernel)
      omegas(:, ik) = tmp_kernel(:)
    end do


  ! the phase of random kernel:
  if (allocated(bcoeffs)) deallocate (bcoeffs); allocate (bcoeffs(dim_KF))
  if (rangml==0) then 
    call random_number(bcoeffs)
    bcoeffs(:) = bcoeffs(:)*two_pi
  end if 
  !TORC! call MPI_BCAST(bcoeffs, size(bcoeffs,1), MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, codeml) 
  call comm_mld%bcast(0, bcoeffs)

  _MLD_END_

end subroutine get_rff_val_sigma


subroutine get_rff_diag_sigma(rff_type, mu, sigma, dim_desc, dim_KF, omegas, bcoeffs)
    !use mpi
    use mld_mpi, only: comm_mld, mld_critical_abort
    ! --------->
    ! rff_type : type of distribution 1 Gaussian ; 2 Laplace, 3 Cauchy 
    ! mu       : \mu of distribution 
    ! sigma    : defines the diagonal elements of the covariance in vector sigma
    ! omegas   : the random omega's: D x K (or F)  size
    ! bcoeffs  : vectors of b between 0:2 pi, of dimension K (or F) 
    use ml_in_ndm_module, only: rangml, two_pi
    use module_kernel, only: krff_gaussian, krff_cauchy, krff_laplace
    use math, only: sample_1D_gaussian_vector_box_muller, sample_1D_laplace_vector, sample_1D_cauchy_vector
    use mld_logger
    implicit none 
    integer, intent(in) :: rff_type, dim_desc, dim_KF
    real(kind_double), intent(in) :: mu
    real(kind_double), dimension(:), allocatable, intent(in) :: sigma
    real(kind_double), dimension(:,:), allocatable, intent(inout) :: omegas 
    real(kind_double), dimension(:), allocatable, intent(inout) :: bcoeffs 
    real(kind_double), dimension(:), allocatable :: tmp_kernel 
    real(kind_double) :: tmp_mu, tmp_sigma 
    integer :: ik, ii 
    _NAMECURRENT_("get_rff_diag_sigma")

    _MLD_BEGIN_

    if (allocated(tmp_kernel)) deallocate (tmp_kernel); allocate (tmp_kernel(dim_desc))
    if (dabs(mu) > 1.d-10) then 
      call log_warning("ML: the median of Random Feature distribution can be only 0. Now in "//NAMECURRENT//" is "// vtoa(mu))
      call log_warning("ML: it will be assigned to 0. ")
    end if 
    tmp_mu = 0.d0   
    tmp_sigma = 1.d0
    do ik = 1, dim_KF
      if (rangml == 0) then
        select case (rff_type)
          case (krff_gaussian)
            call sample_1D_gaussian_vector_box_muller(dim_desc, tmp_kernel, tmp_mu, tmp_sigma)
          case (krff_laplace)
            call sample_1D_laplace_vector(tmp_mu, tmp_sigma, tmp_kernel)
          case (krff_cauchy)
            call sample_1D_cauchy_vector(tmp_mu, tmp_sigma, tmp_kernel)
          case default 
            call log_error("ML error: the choise of krff_type: " // vtoa(rff_type))
            call mld_critical_abort("ML error: there is no choice for this random Fourier Feature in " // NAMECURRENT)   
        end select   
      end if
      !TORC! call MPI_BCAST(tmp_kernel, size(tmp_kernel, 1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, codeml)
      call comm_mld%bcast(0, tmp_kernel)
      do ii = 1, dim_desc
        if (dabs(sigma(ii)) > 1.d-30) then
          omegas(ii, ik) = tmp_kernel(ii)/(sigma(ii)*dsqrt(2.d0))
        else
          omegas(ii, ik) = 0.d0
        end if
      end do
    end do


  ! the phase of random kernel:

  if (allocated(bcoeffs)) deallocate (bcoeffs); allocate (bcoeffs(dim_KF))
  if (rangml==0) then 
    call random_number(bcoeffs)  
    bcoeffs(:) = bcoeffs(:)*two_pi
  end if 
  !TORC! call MPI_BCAST(bcoeffs, size(bcoeffs, 1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, codeml)
  call comm_mld%bcast(0, bcoeffs)
  
  _MLD_END_

end subroutine get_rff_diag_sigma

end module module_sample_rand_ker


subroutine correlation_coef(y1, y2, n, r)
  implicit none
  integer, intent(in)  :: n
  real(kind=kind(1.d0)), intent(in)      :: y1(n), y2(n)
  real(kind=kind(1.d0)), intent(out)     :: r

  real(kind=kind(1.d0))      :: y1_m, y2_m, sigma_y1, sigma_y2,tmp_cross
  if ((n==0).or.(n==1)) then
    r = 1.d0 
    return 
  else   
    y1_m = 1.d0/dble(n)*SUM(y1(:))
    y2_m = 1.d0/dble(n)*SUM(y2(:))
    sigma_y1 = dsqrt(SUM((y1(:) - y1_m)**2)/dble(n))
    sigma_y2 = dsqrt(SUM((y2(:) - y2_m)**2)/dble(n))
    tmp_cross = SUM((y1(:) - y1_m)*(y2(:) - y2_m))
    if ( (dabs(sigma_y2) <= 1.e-20).and.(dabs(tmp_cross) <= 1.e-20) ) then 
      r = 1.d0
    else    
      r = SUM((y1(:) - y1_m)*(y2(:) - y2_m))/dble(n)/(sigma_y1*sigma_y2)
    end if   
  end if 
end subroutine



subroutine determination_coef(y_measure, y_predicted, n, r)
  implicit none
  integer, intent(in)  :: n
  real(kind=kind(1.d0)), intent(in)      :: y_measure(n), y_predicted(n)    ! y1-> database, y2-> fitted function
  real(kind=kind(1.d0)), intent(out)     :: r

  real(kind=kind(1.d0))      :: ym_m, yp_m, sigma_ym, sigma_yp


  ym_m = 1.d0/dble(n)*SUM(y_measure(:))
  yp_m = 1.d0/dble(n)*SUM(y_predicted(:))
  sigma_yp = dsqrt(SUM((y_predicted(:) - yp_m)**2)/dble(n))
  sigma_ym = dsqrt(SUM((y_measure(:) - ym_m)**2)/dble(n))
  if (n > 1) then 
    r = SUM((y_predicted(:) - ym_m)**2)/SUM((y_measure(:) - ym_m)**2)
  else 
    r = 1.d0
  end if     
end subroutine determination_coef



subroutine rmse_mae(y_measure, y_predicted, n, rmse, mae)
  implicit none
  integer, intent(in)  :: n
  real(kind=kind(1.d0)), intent(in)      :: y_measure(n), y_predicted(n)    ! y1-> database, y2-> fitted function
  real(kind=kind(1.d0)), intent(out)     :: rmse, mae

  mae = 1.d0/dble(n)*SUM(dabs(y_measure(:) - y_predicted(:)))
  rmse = dsqrt(1.d0/dble(n)*SUM((y_measure(:) - y_predicted(:))**2))
end subroutine rmse_mae


subroutine recentrate_database(nd_local_data, dim_yfunc, yfunc, yfunc_average)
  implicit none
  integer, intent(in)  :: nd_local_data, dim_yfunc
  real(kind=kind(1.d0))      :: yfunc(dim_yfunc, nd_local_data), yfunc_average(dim_yfunc)
  integer  :: i

  do i = 1, dim_yfunc
    yfunc_average(i) = SUM(yfunc(i, :))/dble(nd_local_data)
  end do
  do i = 1, dim_yfunc
    yfunc(i, :) = yfunc(i, :) - yfunc_average(i)
  end do
end subroutine recentrate_database


subroutine print_message(iproc, text)
  use ml_in_ndm_module, only: rangml
  implicit none
  integer, intent(in)  :: iproc
  character(len=150), intent(in)   :: text

  if (rangml == iproc) write (6, *) text
end subroutine print_message



integer function nitems2(line)
  character, intent(in)      :: line*(*)
  integer i, n, toks

  i = 1;
  n = len_trim(line)
  toks = 0
  nitems2 = 0
  do while (i <= n)
    do while (line(i:i) == ' ')
      i = i + 1
      if (n < i) return
    end do
    toks = toks + 1
    nitems2 = toks
    do
      i = i + 1
      if (n < i) return
      if (line(i:i) == ' ') exit
    end do
  end do
end function nitems2



integer function nitems(line)
  ! number of space-separated items in a line
  implicit none
  character line*(*)
  logical back
  integer length
  integer  :: k

  back = .true.
  line = ' '//line
  length = len_trim(line)
  k = index(line(1:length), ' ', back)
  ! write (*,*) 'nnnnnniiit', k, len_trim(line)
  if (k == 0) then
    nitems = 0
    return
  end if

  nitems = 1
  do
    ! starting with the right most blank space,
    ! look for the next non-space character down
    ! indicating there is another item in the line
    do
      if (k <= 0) exit
      if (line(k:k) == ' ') then
        k = k - 1
        cycle
      else
        nitems = nitems + 1
        exit
      end if
    end do

    ! once a non-space character is found,
    ! skip all adjacent non-space character
    do
      if (k <= 0) exit
      if (line(k:k) /= ' ') then
        k = k - 1
        cycle
      end if
      exit
    end do
    if (k <= 0) exit
  end do
end function nitems



integer function itest_there_is_a_number(line)
  implicit none
  character      :: line*(*)
  character(1), dimension(10)      :: numbers
  character(80)  :: trimline
  integer  :: i

  numbers(1:10) = (/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/)
  trimline = trim(adjustl(line))
  itest_there_is_a_number = 0
  do i = 1, 10
    if (numbers(i) == trimline(1:1)) itest_there_is_a_number = 1
  end do
end function itest_there_is_a_number


function randperm(num)
! This function returns a random permutation of the integers from 1 to num.
  use data_type, only: IB, RP
  implicit none
  integer(kind=IB), intent(in)     :: num
  integer(kind=IB)     :: number, i, j, k
  integer(kind=IB), dimension(num) :: randperm
  real(kind=RP), dimension(num)    :: rand2
  intrinsic random_number
  call random_number(rand2)
  do i = 1, num
    number = 1
    do j = 1, num
      if (rand2(i) > rand2(j)) then
        number = number + 1
      end if
    end do
    do k = 1, i - 1
      if (rand2(i) <= rand2(k) .and. rand2(i) >= rand2(k)) then
        number = number + 1
      end if
    end do
    randperm(i) = number
  end do
end function randperm



module module_eigenvalue_small_matrix
  use iso_fortran_env,  dp => real64
  use mld_logger, only: log_critical, vtoa, mld_verbose, log_debug
  implicit none

  private 

    integer, parameter :: i_tb_svd = 1,  &    ! pure SVD driver 
                          i_tb_syevd=2,  &    ! DSYVED driver 
                          i_tb_syev=3,   &    ! DSYEV driver
                          i_tb_syevr = 4      ! DSYEVR driver


  type  eigen_symmetric 
    logical :: sinit 
    integer :: tb_svd 
    ! the symmetric matrix 
    integer :: dim_AA 
    real(dp), dimension(:,:), allocatable :: AA
    ! the eigenvalues  - SVD-like order 
    real(dp), dimension(:), allocatable  :: eigval
    ! columnwise eigenvectors
    real(dp), dimension(:,:), allocatable  :: eigvec
    integer :: rank

    ! compulsory for internal lapack ...    
    real(dp), dimension(:,:), allocatable :: svdU, svdVT, eigZZ
    real(dp),   dimension(:), allocatable :: svdS, svdWORK, eigWW, eigWWout
    integer,    dimension(:), allocatable :: svdIWORK, issupZ
    integer :: svdINFO, svdLWORK,  svdLIWORK, svdM, svdN 

    contains 
    procedure :: init
    procedure :: evaluate
    procedure :: destroy 
  end type eigen_symmetric 

  type(eigen_symmetric) :: eigensym 
  public :: eigen_symmetric, eigensym

  

  contains

  subroutine init (this, tb_svd, dim_hh)
    class(eigen_symmetric), intent(inout) :: this
    integer, intent(in) :: dim_hh
    integer, intent(in) :: tb_svd
    real(dp)  :: VL, VU, ABSTOL 
    integer :: IL, IU, MM 
    
    _NAMECURRENT_("eigensym%init")
    _MLD_BEGIN_

    if (allocated(this%svdU))  deallocate(this%svdU)  ; allocate(this%svdU(dim_hh, dim_hh))
    if (allocated(this%svdS))  deallocate(this%svdS)  ; allocate(this%svdS(dim_hh))
    if (allocated(this%svdVT)) deallocate(this%svdVT) ; allocate(this%svdVT(dim_hh, dim_hh ))
    if (allocated(this%eigWW)) deallocate(this%eigWW) ; allocate(this%eigWW(dim_hh))
    if (allocated(this%eigWWout)) deallocate(this%eigWWout) ; allocate(this%eigWWout(dim_hh))

    if (allocated(this%eigval)) deallocate(this%eigval) ; allocate(this%eigval(dim_hh))
    if (allocated(this%eigvec)) deallocate(this%eigvec) ; allocate(this%eigvec(dim_hh, dim_hh))
    if (tb_svd == i_tb_syevr) then
      if (allocated(this%eigZZ)) deallocate(this%eigZZ) ; allocate(this%eigZZ(dim_hh, dim_hh))
      if (allocated(this%issupZ)) deallocate(this%issupZ) ; allocate(this%issupZ(2*dim_hh))
    end if


    if (allocated(this%AA)) deallocate(this%AA) ; allocate(this%AA(dim_hh, dim_hh))
    this%tb_svd = tb_svd
    this%AA(:,:) = 0.d0 
    this%dim_AA = dim_hh
    this%svdM = dim_hh
    this%svdN = dim_hh 
    this%sinit = .true. 
    select case (this%tb_svd)

      case  (i_tb_svd)  
        this%svdLWORK = -1 
        if (allocated(this%svdWORK)) deallocate(this%svdWORK) ; allocate(this%svdWORK(1))
        call dgesvd('A', 'A', this%svdM, this%svdN, this%AA, this%svdM, this%svdS, this%svdU, this%svdM, &
                     this%svdVT, this%svdN, this%svdWORK, this%svdLWORK, this%svdINFO)
        if (allocated(this%svdWORK)) deallocate(this%svdWORK) ; allocate(this%svdWORK(this%svdLWORK))
      case (i_tb_syevd)
        this%svdLWORK = -1
        this%svdLIWORK = -1
        if (allocated(this%svdWORK))  deallocate(this%svdWORK)  ; allocate(this%svdWORK(1))
        if (allocated(this%svdIWORK)) deallocate(this%svdIWORK) ; allocate(this%svdIWORK(1))
        call dsyevd('V', 'U', this%svdN, this%AA, this%svdN, this%eigWW, this%svdWORK, this%svdLWORK, & 
                     this%svdIWORK, this%svdLIWORK, this%svdINFO)
        this%svdLWORK = int(this%svdWORK(1))
        this%svdLIWORK = int(this%svdIWORK(1))
        if (allocated(this%svdWORK))  deallocate(this%svdWORK)  ; allocate(this%svdWORK( this%svdLWORK))
        if (allocated(this%svdIWORK)) deallocate(this%svdIWORK) ; allocate(this%svdIWORK(this%svdLIWORK))
      case (i_tb_syev)
        this%svdLWORK = -1
        if (allocated(this%svdWORK))  deallocate(this%svdWORK)  ; allocate(this%svdWORK(1))
        call dsyev('V', 'U', this%svdN, this%AA, this%svdN, this%eigWW, this%svdWORK, this%svdLWORK,  this%svdINFO)
        this%svdLWORK = int(this%svdWORK(1))
        if (allocated(this%svdWORK))  deallocate(this%svdWORK)  ; allocate(this%svdWORK( this%svdLWORK))
      case (i_tb_syevr)
        this%svdLWORK = -1
        this%svdLIWORK = -1
        VL = 1.d0
        VU = 1.d0
        IL = 1
        IU = 1 
        MM = 1 
        ABSTOL = -1.d0 
        if (allocated(this%svdWORK))  deallocate(this%svdWORK)  ; allocate(this%svdWORK(1))
        if (allocated(this%svdIWORK)) deallocate(this%svdIWORK) ; allocate(this%svdIWORK(1))
        call dsyevr('V', 'A', 'U', this%svdN, this%AA, this%svdN, VL, VU, IL, IU, ABSTOL, MM, this%eigWW, this%eigZZ, this%svdN, this%issupZ, this%svdWORK, this%svdLWORK, &
                     this%svdIWORK, this%svdLIWORK, this%svdINFO)
        call log_debug('Eigenvalue info: '//vtoa(this%svdINFO))
        if (this%svdINFO /= 0) then 
          call log_critical('Eigen dsyevr I  failed in ...'//NAMECURRENT//"with the error code:  "//vtoa(this%svdINFO)//" !!!")
          stop 'eigenvalue info error'
        end if              
        this%svdLWORK = int(this%svdWORK(1))
        this%svdLIWORK = int(this%svdIWORK(1))
        if (allocated(this%svdWORK))  deallocate(this%svdWORK)  ; allocate(this%svdWORK( this%svdLWORK))
        if (allocated(this%svdIWORK)) deallocate(this%svdIWORK) ; allocate(this%svdIWORK(this%svdLIWORK))  
    end select

    _MLD_END_
  end subroutine init

  subroutine evaluate (this, hh)
    use, intrinsic :: ieee_arithmetic
    class(eigen_symmetric), intent(inout) :: this
    real(dp) , dimension(:,:), intent(in) :: hh 
    real(dp) :: rcond_local
    integer :: irang_svd, ii, ian, idx_reverse
    real(dp), dimension(:), allocatable :: abs_svdS
    integer, dimension(:), allocatable :: indx_abs 
    real(dp)  :: VL, VU, ABSTOL 
    integer :: IL, IU, MM 
    

    _NAMECURRENT_("eigensym%evaluate")
    _MLD_BEGIN_

    this%AA(:,:) = hh(:,:)
    select case (this%tb_svd) 
      case  (i_tb_svd)  
        call dgesvd('A', 'A', this%svdM, this%svdN, this%AA, this%svdM, this%svdS, this%svdU, this%svdM, &
                    this%svdVT, this%svdN, this%svdWORK, this%svdLWORK, this%svdINFO)
        if (this%svdINFO /= 0) then 
          call log_critical('SVD dgesvd failed in ...'//NAMECURRENT//"with the error code:  "//vtoa(this%svdINFO)//" !!!")
          stop 'svd info error'
        end if
      case (i_tb_syevd)
        call dsyevd('V', 'U', this%svdN, this%AA, this%svdN, this%eigWW, this%svdWORK, this%svdLWORK, &
                     this%svdIWORK, this%svdLIWORK, this%svdINFO)
        call log_debug('Eigenvalue info: '//vtoa(this%svdINFO))
        if (this%svdINFO /= 0) then 
          call log_critical('Eigen dsyevd failed in ...'//NAMECURRENT//"with the error code:  "//vtoa(this%svdINFO)//" !!!")
          stop 'eigenvalue info error'
        end if
      case (i_tb_syev)
        call dsyev('V', 'U', this%svdN, this%AA, this%svdN, this%eigWW, this%svdWORK, this%svdLWORK,  this%svdINFO)
        call log_debug('Eigenvalue info: '//vtoa(this%svdINFO))
        if (this%svdINFO /= 0) then 
          call log_critical('Eigen dsyev in ...'//NAMECURRENT//"with the error code:  "//vtoa(this%svdINFO)//" !!!")
          stop 'eigenvalue info error'
        end if
      case (i_tb_syevr)
        VL = 1.d0
        VU = 1.d0
        IL = 1
        IU = 1
        ABSTOL = -1.D0
        MM = 1 

        call dsyevr('V', 'A', 'U', this%svdN, this%AA, this%svdN, VL, VU, IL, IU, ABSTOL, MM, this%eigWW, this%eigZZ, this%svdN, this%issupZ, this%svdWORK, this%svdLWORK, &
                     this%svdIWORK, this%svdLIWORK, this%svdINFO)
   
        call log_debug('Eigenvalue info: '//vtoa(this%svdINFO))
        if (this%svdINFO /= 0) then 
          call log_critical('Eigen dsyevr failed in ...'//NAMECURRENT//"with the error code:  "//vtoa(this%svdINFO)//" !!!")
          stop 'eigenvalue info error'
        end if  
    end select

    if (this%tb_svd /= i_tb_svd) then
      if (allocated(abs_svdS)) deallocate(abs_svdS) ; allocate(abs_svdS(this%dim_AA))
      abs_svdS(:) = dabs(this%eigWW(:))
      if (allocated(indx_abs)) deallocate(indx_abs) ; allocate(indx_abs(this%dim_AA))
      call indexx(this%dim_AA, abs_svdS, indx_abs)
      if (this%tb_svd == i_tb_syevr) then
        this%AA(:,:) = this%eigZZ(:,:)  
      end if

      !write(*,*) '-------------------------'
      do ii = 1, this%dim_AA
         idx_reverse = indx_abs(this%dim_AA - ii + 1)
         this%eigWWout(ii) = this%eigWW(idx_reverse)
         !this%svdS(ii) = this%eigWW(idx_reverse)
         !this%svdS(ii) = abs_svdS(idx_reverse)
         !this%svdS(ii) = sign(1.d0, this%eigWW(idx_reverse)) * sqrt(abs_svdS(idx_reverse))
         this%svdS(ii) = abs_svdS(idx_reverse)**2 
         !this%svdS(ii) = this%eigWW(idx_reverse)
         !this%svdS(ii) = sign(1.d0, this%eigWW(idx_reverse)) * sqrt(abs_svdS(idx_reverse)**2 + 0.1d0)
         

         !this%svdU(:,ii) = sign(1.d0, this%eigWW(idx_reverse)) * this%AA(:, idx_reverse)     
         this%svdU(:,ii) = this%AA(:, idx_reverse)     
         !write(*,*) 'svdS', ii, this%svdS(ii) 
      end do 
    end if 

    rcond_local = dabs(this%svdS(1))*1.d-14 
    irang_svd = 0
    do ian = 1, this%dim_AA
      if (dabs(this%svdS(ian)) > rcond_local) then 
        irang_svd = irang_svd + 1
      else 
        this%svdS(ian) = 0.d0
      end if   
    end do
    this%eigval = this%svdS(:)
    this%eigvec = this%svdU(:,:)
    this%rank = irang_svd

    _MLD_END_

  end subroutine evaluate 

  subroutine destroy (this)
    class(eigen_symmetric), intent(inout) :: this


    if (allocated(this%svdU))  deallocate(this%svdU)  
    if (allocated(this%svdS))  deallocate(this%svdS)  
    if (allocated(this%svdVT)) deallocate(this%svdVT) 
    if (allocated(this%eigWW)) deallocate(this%eigWW) 
    if (allocated(this%eigWWout)) deallocate(this%eigWWout) 

    if (allocated(this%AA)) deallocate(this%AA) 
    if (allocated(this%svdIWORK)) deallocate(this%svdIWORK)
    if (allocated(this%svdWORK)) deallocate(this%svdWORK)
    if (allocated(this%issupZ)) deallocate(this%issupZ)
    if (allocated(this%eigZZ)) deallocate(this%eigZZ)
    
     
    this%dim_AA = -1
    this%svdM = -1
    this%svdN = -1 
 

    if (allocated(this%eigval)) deallocate(this%eigval)
    if (allocated(this%eigvec)) deallocate(this%eigvec)
    
  end subroutine destroy


end module module_eigenvalue_small_matrix

module module_svd_small_matrix
  use iso_fortran_env,  dp => real64
  use mld_logger, only: log_critical, vtoa, mld_verbose, log_debug
  implicit none

  private 

  type  eigen_svd 
    ! the symmetric matrix 
    integer :: dim_AA 
    real(dp), dimension(:,:), allocatable :: AA
    ! the eigenvalues  - SVD-like order 
    real(dp), dimension(:), allocatable  :: eigval
    ! columnwise eigenvectors
    real(dp), dimension(:,:), allocatable  :: eigvecU
    real(dp), dimension(:,:), allocatable  :: eigvecVT
    integer :: rank

    ! compulsory for internal lapack ...    
    real(dp), dimension(:,:), allocatable :: svdU, svdVT
    real(dp),   dimension(:), allocatable :: svdS, svdWORK, eigWW, vv_ref
    integer,    dimension(:), allocatable :: svdIWORK
    integer :: svdINFO, svdLWORK,  svdLIWORK, svdM, svdN 

    contains 
    procedure :: init
    procedure :: evaluate
    procedure :: destroy 
  end type eigen_svd 

  type(eigen_svd) :: eigensvd 
  public :: eigen_svd, eigensvd

  contains

  subroutine init (this, dim_hh)
    class(eigen_svd), intent(inout) :: this
    !logical, intent(in) :: tb_svd
    integer, intent(in) :: dim_hh
    
    _NAMECURRENT_("eigensvd%init")
    _MLD_BEGIN_

    if (allocated(this%svdU))  deallocate(this%svdU)  ; allocate(this%svdU(dim_hh, dim_hh))
    if (allocated(this%svdS))  deallocate(this%svdS)  ; allocate(this%svdS(dim_hh))

    if (allocated(this%eigval)) deallocate(this%eigval) ; allocate(this%eigval(dim_hh))
    if (allocated(this%eigvecU)) deallocate(this%eigvecU) ; allocate(this%eigvecU(dim_hh, dim_hh))
    if (allocated(this%eigvecVT)) deallocate(this%eigvecVT) ; allocate(this%eigvecVT(dim_hh, dim_hh))

    if (allocated(this%svdVT)) deallocate(this%svdVT) ; allocate(this%svdVT(dim_hh, dim_hh ))
    if (allocated(this%eigWW)) deallocate(this%eigWW) ; allocate(this%eigWW(dim_hh))
    if (allocated(this%vv_ref)) deallocate(this%vv_ref) ; allocate(this%vv_ref(dim_hh))
    call make_default_vref(dim_hh, this%vv_ref)

    if (allocated(this%AA)) deallocate(this%AA) ; allocate(this%AA(dim_hh, dim_hh))
    
    this%AA(:,:) = 0.d0 
    this%dim_AA = dim_hh
    this%svdM = dim_hh
    this%svdN = dim_hh 

    
    this%svdLWORK = -1 
    if (allocated(this%svdWORK)) deallocate(this%svdWORK) ; allocate(this%svdWORK(1))
    call dgesvd('A', 'A', this%svdM, this%svdN, this%AA, this%svdM, this%svdS, this%svdU, this%svdM, &
                 this%svdVT, this%svdN, this%svdWORK, this%svdLWORK, this%svdINFO)
    this%svdLWORK = int(this%svdWORK(1))
    if (allocated(this%svdWORK)) deallocate(this%svdWORK) ; allocate(this%svdWORK(this%svdLWORK))
    

    _MLD_END_
  end subroutine init

  subroutine evaluate (this, hh)
    class(eigen_svd), intent(inout) :: this
    real(dp) , dimension(:,:), intent(in) :: hh 
    real(dp) :: rcond_local
    integer :: irang_svd, ian 
    _NAMECURRENT_("eigensym%evaluate")
    _MLD_BEGIN_

    this%AA(:,:) = hh(:,:)
    call dgesvd('A', 'A', this%svdM, this%svdN, this%AA, this%svdM, this%svdS, this%svdU, this%svdM, &
                this%svdVT, this%svdN, this%svdWORK, this%svdLWORK, this%svdINFO)
    call log_debug('Eigenvalue info: '//vtoa(this%svdINFO))

    if (this%svdINFO /= 0) then 
      call log_critical('SVD failed in ...'//NAMECURRENT//"with the error code:  "//vtoa(this%svdINFO)//" !!!")
      stop 'eigenvalue info error'
    end if

    call svd_enforce_signs_by_ref(this%svdS, this%svdU, this%svdVT, this%vv_ref)

    rcond_local = dabs(this%svdS(1))*1.e-14_dp  
    irang_svd = 0
    do ian = 1, this%dim_AA
      if (dabs(this%svdS(ian)) > rcond_local) then 
        irang_svd = irang_svd + 1
      else 
        this%svdS(ian) = 0.d0
      end if   
    end do

    
    this%eigval = this%svdS(:)
    this%eigvecU = this%svdU(:,:)
    this%eigvecVT = this%svdVT(:,:)
    this%rank = irang_svd

    _MLD_END_

  end subroutine evaluate 

  subroutine destroy (this)
    class(eigen_svd), intent(inout) :: this


    if (allocated(this%svdU))  deallocate(this%svdU)  
    if (allocated(this%svdS))  deallocate(this%svdS)  
    if (allocated(this%svdVT)) deallocate(this%svdVT) 
    if (allocated(this%eigWW)) deallocate(this%eigWW) 

    if (allocated(this%AA)) deallocate(this%AA) 
    
     
    this%dim_AA = -1
    this%svdM = -1
    this%svdN = -1 
 
    if (allocated(this%svdWORK)) deallocate(this%svdWORK)

    if (allocated(this%eigval)) deallocate(this%eigval)
    if (allocated(this%eigvecU)) deallocate(this%eigvecU)
    if (allocated(this%eigvecVT)) deallocate(this%eigvecVT)
    
  end subroutine destroy

  subroutine make_default_vref(n, v_ref)
    integer, intent(in) :: n
    real(dp), intent(inout) :: v_ref(:)
    real(dp) :: nr

    v_ref(:) = 1.0_dp
    nr = sqrt(sum(v_ref*v_ref))
    if (nr > 0.0_dp) v_ref = v_ref / nr
  end subroutine make_default_vref


  subroutine svd_enforce_signs_by_ref(s, U, VT, v_ref,  tol)
    ! Enforce a deterministic sign on singular vectors using a fixed reference v_ref.
    ! Also handles degenerate subspaces by diagonalizing a reference operator.
    !
    ! Inputs/Outputs:
    !   s    : (min(m,n)) singular values
    !   U    : (m x m) left singular vectors (modified in place)
    !   VT   : (n x n) right singular vectors transposed (modified in place)
    !   v_ref: (n)     fixed reference direction in feature space (need not be unit; we normalize)
    !
    ! Optional:
    !   tol  : small tolerance to detect near-zero projection; default 1e-14
    !
    real(dp), intent(in) :: s(:)
    real(dp), intent(inout) :: U(:,:), VT(:,:)     ! U(m,m), VT(n,n)
    real(dp), intent(inout)    :: v_ref(:)            ! length n = size(VT,2)
    real(dp), intent(in), optional :: tol

    integer :: m, n, k, kstop, idx
    real(dp) :: nr, proj, mytol

    ! New variables for degeneracy handling
    integer :: i, j, d, info, lwork, p, q, q1, q2
    real(dp) :: tol_degen, sum_val
    real(dp), allocatable :: V_block(:,:), H(:,:), evals_H(:), evecs_H(:,:), work(:), U_block(:,:), VT_block(:,:)

    m = size(U, 1)
    n = size(VT,2)

    if (size(U,2) < min(m,n)) stop 'svd_enforce_signs_by_ref: U has too few columns'
    if (size(VT,1) < min(m,n)) stop 'svd_enforce_signs_by_ref: VT has too few rows'
    if (size(v_ref) /= n)      stop 'svd_enforce_signs_by_ref: v_ref has wrong length'

    kstop = min(min(m,n), size(VT,1))
    kstop = min(kstop, size(s))

    mytol = 1.0e-14_dp
    if (present(tol)) mytol = tol

    ! Handle degenerate subspaces
    if (size(s) > 0) then
       tol_degen = 1.0e-12_dp * s(1)
       i = 1
       do while (i <= kstop)
          j = i + 1
          do while (j <= kstop)
             if (abs(s(i) - s(j)) >= tol_degen) exit
             j = j + 1
          end do
          
          if (j > i + 1) then
             d = j - i
             allocate(V_block(n, d))
             allocate(H(d, d))
             allocate(evecs_H(d, d))
             allocate(evals_H(d))
             allocate(work(3*d))

             ! Extract V_block (columns are basis vectors)
             ! VT contains V^T. So rows of VT are basis vectors.
             do q = 1, d
                do p = 1, n
                   V_block(p, q) = VT(i+q-1, p)
                end do
             end do

             ! Construct H = V^T M V with M = diag(1..n)
             H = 0.0_dp
             do q1 = 1, d
                do q2 = 1, d
                   sum_val = 0.0_dp
                   do p = 1, n
                      sum_val = sum_val + V_block(p, q1) * real(p, dp) * V_block(p, q2)
                   end do
                   H(q1, q2) = sum_val
                end do
             end do

             ! Diagonalize H
             lwork = 3*d
             call dsyev('V', 'U', d, H, d, evals_H, work, lwork, info)
             evecs_H = H

             ! Update VT (rows i..j-1)
             allocate(VT_block(d, n))
             VT_block = 0.0_dp
             do p = 1, n
                do q1 = 1, d
                   sum_val = 0.0_dp
                   do q2 = 1, d
                      sum_val = sum_val + evecs_H(q2, q1) * VT(i+q2-1, p)
                   end do
                   VT_block(q1, p) = sum_val
                end do
             end do
             do q = 1, d
                VT(i+q-1, :) = VT_block(q, :)
             end do
             deallocate(VT_block)

             ! Update U (columns i..j-1)
             allocate(U_block(m, d))
             U_block = 0.0_dp
             do p = 1, m
                do q1 = 1, d
                   sum_val = 0.0_dp
                   do q2 = 1, d
                      sum_val = sum_val + U(p, i+q2-1) * evecs_H(q2, q1)
                   end do
                   U_block(p, q1) = sum_val
                end do
             end do
             do q = 1, d
                U(:, i+q-1) = U_block(:, q)
             end do
             deallocate(U_block)

             deallocate(V_block, H, evecs_H, evals_H, work)
          end if
          i = j
       end do
    end if

    do k = 1, kstop
       proj = dot_product(VT(k,1:n), v_ref)
       if (proj < 0.0_dp) then
          U(:,k)  = -U(:,k)
          VT(k,:) = -VT(k,:)
       else if (abs(proj) <= mytol) then
          ! tie-breaker: use the largest-magnitude entry rule
          idx = maxloc(abs(VT(k,1:n)), dim=1)
          if (VT(k, idx) < 0.0_dp) then
             U(:,k)  = -U(:,k)
             VT(k,:) = -VT(k,:)
          end if
       end if
    end do

    
  end subroutine svd_enforce_signs_by_ref

end module module_svd_small_matrix


module module_svd_small_gen_matrix
  use iso_fortran_env,  dp => real64
  use mld_logger, only: log_critical, vtoa, mld_verbose, log_debug
  implicit none

  private 

  type  eigen_svd 
    ! the general matrix 
    integer :: dim_M_AA, dim_N_AA, dim_S 
    real(dp), dimension(:,:), allocatable :: AA
    ! the eigenvalues  - SVD-like order 
    real(dp), dimension(:), allocatable  :: eigval
    ! columnwise eigenvectors
    real(dp), dimension(:,:), allocatable  :: eigvecU
    real(dp), dimension(:,:), allocatable  :: eigvecVT
    integer :: rank

    ! compulsory for internal lapack ...    
    real(dp), dimension(:,:), allocatable :: svdU, svdVT
    real(dp),   dimension(:), allocatable :: svdS, svdWORK, eigWW, vv_ref 
    integer,    dimension(:), allocatable :: svdIWORK
    integer :: svdINFO, svdLWORK,  svdLIWORK, svdM, svdN 

    contains 
    procedure :: init
    procedure :: evaluate
    procedure :: destroy 
  end type eigen_svd 

  type(eigen_svd) :: eigensvd 
  public :: eigen_svd, eigensvd

  contains

  subroutine init (this, dim_M_AA, dim_N_AA)
    class(eigen_svd), intent(inout) :: this
    !logical, intent(in) :: tb_svd
    integer, intent(in) :: dim_M_AA, dim_N_AA
    
    _NAMECURRENT_("eigensvd%init")
    _MLD_BEGIN_

    if (allocated(this%svdU))  deallocate(this%svdU)  ; allocate(this%svdU(dim_M_AA, dim_M_AA))
    this%dim_S = min(dim_M_AA, dim_N_AA)
    if (allocated(this%svdS))  deallocate(this%svdS)  ; allocate(this%svdS(this%dim_S))

    if (allocated(this%eigval)) deallocate(this%eigval) ; allocate(this%eigval(this%dim_S))
    if (allocated(this%eigvecU)) deallocate(this%eigvecU) ; allocate(this%eigvecU(dim_M_AA, dim_M_AA))
    if (allocated(this%eigvecVT)) deallocate(this%eigvecVT) ; allocate(this%eigvecVT(dim_N_AA, dim_N_AA))
    if (allocated(this%vv_ref)) deallocate(this%vv_ref) ; allocate(this%vv_ref(dim_N_AA))

    if (allocated(this%svdVT)) deallocate(this%svdVT) ; allocate(this%svdVT(dim_N_AA, dim_N_AA ))
    if (allocated(this%eigWW)) deallocate(this%eigWW) ; allocate(this%eigWW(this%dim_S))

    if (allocated(this%AA)) deallocate(this%AA) ; allocate(this%AA(dim_M_AA, dim_N_AA))
    
    this%AA(:,:) = 0.d0 
    this%dim_M_AA = dim_M_AA
    this%svdM = dim_M_AA
    this%svdN = dim_N_AA

    
    this%svdLWORK = -1 
    if (allocated(this%svdWORK)) deallocate(this%svdWORK) ; allocate(this%svdWORK(1))
    call dgesvd('A', 'A', this%svdM, this%svdN, this%AA, this%svdM, this%svdS, this%svdU, this%svdM, &
                 this%svdVT, this%svdN, this%svdWORK, this%svdLWORK, this%svdINFO)
    this%svdLWORK = int(this%svdWORK(1))
    if (allocated(this%AA)) deallocate(this%AA) ; allocate(this%AA(dim_M_AA, dim_N_AA))
    if (allocated(this%svdWORK)) deallocate(this%svdWORK) ; allocate(this%svdWORK(this%svdLWORK))
    ! write(*,*) "!!!!!!!!!!!!!!! SVD LWORK is ", this%svdLWORK

    call make_default_vref(this%dim_N_AA, this%vv_ref)
    _MLD_END_
  end subroutine init

  subroutine evaluate (this, hh)
    class(eigen_svd), intent(inout) :: this
    real(dp) , dimension(:,:), intent(in) :: hh 
    real(dp) :: rcond_local
    integer :: irang_svd, ian 

    _NAMECURRENT_("eigensym%evaluate")
    _MLD_BEGIN_

    this%AA(:,:) = hh(:,:)
    call dgesvd('A', 'A', this%svdM, this%svdN, this%AA, this%svdM, this%svdS, this%svdU, this%svdM, &
                this%svdVT, this%svdN, this%svdWORK, this%svdLWORK, this%svdINFO)

    call log_debug('Eigenvalue info: '//vtoa(this%svdINFO))

    if (this%svdINFO /= 0) then 
      call log_critical('SVD failed in ...'//NAMECURRENT//"with the error code:  "//vtoa(this%svdINFO)//" !!!")
      stop 'eigenvalue info error'
    end if

    !first call svd_enforce_signs(this%svdS, this%svdU, this%svdVT)


    rcond_local = dabs(this%svdS(1))*1.e-14_dp  
    call svd_enforce_signs_by_ref(this%svdS, this%svdU, this%svdVT, this%vv_ref, rcond_local)
    irang_svd = 0
    do ian = 1, size(this%svdS)
      if (dabs(this%svdS(ian)) > rcond_local) then 
        irang_svd = irang_svd + 1
      else 
        this%svdS(ian) = 0.d0
      end if   
    end do

    
    this%eigval = this%svdS(:)
    this%eigvecU = this%svdU(:,:)
    this%eigvecVT = this%svdVT(:,:)
    this%rank = irang_svd
    ! write(*,*) "11111111: eigvecVT", this%eigvecVT(1,:)
    ! write(*,*) "11111111  svdVT", this%svdVT(1,:)
    ! write(*,*) "22222222  svdVT", this%svdVT(2,:)

    _MLD_END_

  end subroutine evaluate 

  subroutine make_default_vref(n, v_ref)
    integer, intent(in) :: n
    real(dp), intent(inout) :: v_ref(:)
    real(dp) :: nr

    v_ref(:) = 1.0_dp
    nr = sqrt(sum(v_ref*v_ref))
    if (nr > 0.0_dp) v_ref = v_ref / nr
  end subroutine make_default_vref


  !subroutine svd_enforce_signs_by_ref(U, VT, v_ref, kmax, tol)
  subroutine svd_enforce_signs_by_ref(s, U, VT, v_ref,  tol)
    ! Enforce a deterministic sign on singular vectors using a fixed reference v_ref.
    ! Also handles degenerate subspaces by diagonalizing a reference operator.
    !
    ! Inputs/Outputs:
    !   s    : (min(m,n)) singular values
    !   U    : (m x m) left singular vectors (modified in place)
    !   VT   : (n x n) right singular vectors transposed (modified in place)
    !   v_ref: (n)     fixed reference direction in feature space (need not be unit; we normalize)
    !
    ! Optional:
    !   tol  : small tolerance to detect near-zero projection; default 1e-14
    !
    real(dp), intent(in) :: s(:)
    real(dp), intent(inout) :: U(:,:), VT(:,:)     ! U(m,m), VT(n,n)
    real(dp), intent(inout)    :: v_ref(:)            ! length n = size(VT,2)
    real(dp), intent(in), optional :: tol

    integer :: m, n, k, kstop, idx
    real(dp) :: nr, proj, mytol

    ! New variables for degeneracy handling
    integer :: i, j, d, info, lwork, p, q, q1, q2
    real(dp) :: tol_degen, sum_val
    real(dp), allocatable :: V_block(:,:), H(:,:), evals_H(:), evecs_H(:,:), work(:), U_block(:,:), VT_block(:,:)

    m = size(U, 1)
    n = size(VT,2)

    if (size(U,2) < min(m,n)) stop 'svd_enforce_signs_by_ref: U has too few columns'
    if (size(VT,1) < min(m,n)) stop 'svd_enforce_signs_by_ref: VT has too few rows'
    if (size(v_ref) /= n)      stop 'svd_enforce_signs_by_ref: v_ref has wrong length'

    kstop = min(min(m,n), size(VT,1))
    kstop = min(kstop, size(s))

    mytol = 1.0e-14_dp
    if (present(tol)) mytol = tol

    ! Handle degenerate subspaces
    if (size(s) > 0) then
       tol_degen = 1.0e-12_dp * s(1)
       i = 1
       do while (i <= kstop)
          j = i + 1
          do while (j <= kstop)
             if (abs(s(i) - s(j)) >= tol_degen) exit
             j = j + 1
          end do
          
          if (j > i + 1) then
             d = j - i
             allocate(V_block(n, d))
             allocate(H(d, d))
             allocate(evecs_H(d, d))
             allocate(evals_H(d))
             allocate(work(3*d))

             ! Extract V_block (columns are basis vectors)
             ! VT contains V^T. So rows of VT are basis vectors.
             do q = 1, d
                do p = 1, n
                   V_block(p, q) = VT(i+q-1, p)
                end do
             end do

             ! Construct H = V^T M V with M = diag(1..n)
             H = 0.0_dp
             do q1 = 1, d
                do q2 = 1, d
                   sum_val = 0.0_dp
                   do p = 1, n
                      sum_val = sum_val + V_block(p, q1) * real(p, dp) * V_block(p, q2)
                   end do
                   H(q1, q2) = sum_val
                end do
             end do

             ! Diagonalize H
             lwork = 3*d
             call dsyev('V', 'U', d, H, d, evals_H, work, lwork, info)
             evecs_H = H

             ! Update VT (rows i..j-1)
             allocate(VT_block(d, n))
             VT_block = 0.0_dp
             do p = 1, n
                do q1 = 1, d
                   sum_val = 0.0_dp
                   do q2 = 1, d
                      sum_val = sum_val + evecs_H(q2, q1) * VT(i+q2-1, p)
                   end do
                   VT_block(q1, p) = sum_val
                end do
             end do
             do q = 1, d
                VT(i+q-1, :) = VT_block(q, :)
             end do
             deallocate(VT_block)

             ! Update U (columns i..j-1)
             allocate(U_block(m, d))
             U_block = 0.0_dp
             do p = 1, m
                do q1 = 1, d
                   sum_val = 0.0_dp
                   do q2 = 1, d
                      sum_val = sum_val + U(p, i+q2-1) * evecs_H(q2, q1)
                   end do
                   U_block(p, q1) = sum_val
                end do
             end do
             do q = 1, d
                U(:, i+q-1) = U_block(:, q)
             end do
             deallocate(U_block)

             deallocate(V_block, H, evecs_H, evals_H, work)
          end if
          i = j
       end do
    end if

    do k = 1, kstop
       proj = dot_product(VT(k,1:n), v_ref)
       if (proj < 0.0_dp) then
          U(:,k)  = -U(:,k)
          VT(k,:) = -VT(k,:)
       else if (abs(proj) <= mytol) then
          ! tie-breaker: use the largest-magnitude entry rule
          idx = maxloc(abs(VT(k,1:n)), dim=1)
          if (VT(k, idx) < 0.0_dp) then
             U(:,k)  = -U(:,k)
             VT(k,:) = -VT(k,:)
          end if
       end if
    end do

    
  end subroutine svd_enforce_signs_by_ref


  subroutine svd_enforce_signs(s, U, VT)
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in)    :: s(:)
    real(dp), intent(inout) :: U(:,:), VT(:,:) ! U(m,m), VT(n,n)
    integer :: k, n, idx
    real(dp), allocatable :: row(:)
    n = size(VT,2)
    do k = 1, size(s)
      row = VT(k,1:n)                    ! v_k^T
      idx = maxloc(abs(row), dim=1)      ! position of |max| component
      if (row(idx) < 0.0_dp) then
        U(:,k)  = -U(:,k)
        VT(k,:) = -VT(k,:)
      end if
    end do
  end subroutine

  subroutine destroy (this)
    class(eigen_svd), intent(inout) :: this


    if (allocated(this%svdU))  deallocate(this%svdU)  
    if (allocated(this%svdS))  deallocate(this%svdS)  
    if (allocated(this%svdVT)) deallocate(this%svdVT) 
    if (allocated(this%eigWW)) deallocate(this%eigWW) 

    if (allocated(this%AA)) deallocate(this%AA) 
    
     
    this%dim_M_AA = -1
    this%dim_N_AA = -1
    this%dim_S = -1
    this%rank = -1
    this%svdM = -1
    this%svdN = -1 
 
    if (allocated(this%svdWORK)) deallocate(this%svdWORK)

    if (allocated(this%eigval)) deallocate(this%eigval)
    if (allocated(this%eigvecU)) deallocate(this%eigvecU)
    if (allocated(this%eigvecVT)) deallocate(this%eigvecVT)
    
  end subroutine destroy

end module module_svd_small_gen_matrix
