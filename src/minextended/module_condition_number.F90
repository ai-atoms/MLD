
#include "../MLD_MACROS.INC"

module module_condition_number
  use mld_logger, only: vtoa, log_info, mld_verbose, log_critical, log_debug 
  use mld_mpi, only: mld_mpi_abort 
  implicit none

contains

  function compute_condition_number_mat_rss(a) result(cond_num)
    use module_kind_variables, only: kind_double 
    real(kind_double), intent(in) :: A(:,:)
    real(kind_double) :: cond_num
    real(kind_double), allocatable :: work(:)
    integer, allocatable :: iwork(:)
    real(kind_double) :: anorm, rcond
    integer :: n,  info
    character :: uplo
    real(kind_double), external  :: dlansy
    real(kind_double), dimension(:,:), allocatable :: Acopy 


    _NAMECURRENT_('compute_condition_number_mat_rss')
    _MLD_BEGIN_ 
    n = size(A, dim=1)   ! Assuming a is square
    if (n /= size(A,2)) then 
      call log_critical('Matrix is not square in ' //NAMECURRENT)
      call mld_mpi_abort('stop non-square in '//NAMECURRENT)  
    end if 

    uplo = 'u'   ! Using upper triangular part of a for symmetry

    ! Allocate workspace
    allocate(work(3*n), iwork(n))

    ! Compute the infinity-norm of a
    anorm = dlansy('F', uplo, size(A, dim=1), A, size(A,dim=1), work)

    ! Perform Cholesky factorization
    allocate(Acopy, source=A)
    call dpotrf(uplo, size(Acopy, dim=1), Acopy, size(Acopy, dim=1), info)

    if (info == 0) then
      ! Compute the reciprocal of the condition number
      call dpocon(uplo, size(Acopy, dim=1), Acopy, size(Acopy, dim=1), anorm, rcond, work, iwork, info)
      cond_num = 1.0 / rcond

    else
      print *, "Failed to factorize a"
      cond_num = -1.0  ! Indicate an error
    endif

    ! Deallocate
    deallocate(work, iwork, Acopy)

    _MLD_END_

  end function compute_condition_number_mat_rss

end module module_condition_number

