

program MinimizerExample
  use MinimizerModule, only: MinimizerType
  use SteepestDescentModule, only: SteepestDescentType
  use ConjugateGradientModule, only : ConjugateGradientType
  use ConjugateGradientExtModule, only : ConjugateGradientExtType
  use LbgfsExtModule, only: LbfgsExtType
  use FunctionModule, only: FunctionType
  use FuncQuadraticModule, only: FuncQuadraticType
  use FuncCubicModule, only: FuncCubicType
  use module_lbfgs_nocedal, only : lbfgs_m_hess, lbfgs_xtol
  use module_condition_number, only : compute_condition_number_mat_rss
  use module_random_matrix, only : generate_random_vector, generate_spd_matrix
  
  implicit none
  
  integer, parameter :: ndim = 3
  real(kind=8), dimension(:,:), allocatable  :: A
  real(kind=8), dimension(:), allocatable  :: b, x, xini, xend, df
  real(kind=8) :: f_value, cond_number
  !$! class(FunctionType), allocatable :: func
  class(FuncQuadraticType), allocatable :: qfunc
  class(FuncCubicType), allocatable ::  cfunc
  !$! class(MinimizerType), allocatable :: minimizer
  class(MinimizerType), allocatable :: minimizer, cgExt
  
  ! Define stadard fix toy model for A and b for QuadraticFunction
  allocate(b(ndim), x(ndim), df(ndim), A(ndim, ndim) )
  A = reshape([2.0, 0.5, 0.0, 0.5, 3.0, 0.0, 0.0, 0.0, 1.0], shape(A))
  b = [1.0, -2.0, 0.5]
  ! Define an initial guess for x
  x = [1.0, 2.0, 3.0]


  ! Define random A and b
  !call generate_spd_matrix(ndim, A)
  !cond_number = compute_condition_number_mat_rss (A) 
  !write(*,'("The condition number of A matrix is ....:", es25.14)') cond_number 
  !if (allocated(b))  deallocate(b)  ; allocate(b(ndim))
  !call generate_random_vector(ndim, b) 
  !x = [1.0, 2.0, 3.0]


  allocate(xini, source=x)
  allocate(xend, source=x)
  
  ! Allocate and initialize quadratic function object
  allocate(FuncQuadraticType :: qfunc)
  qfunc%fdim = ndim 
  allocate(qfunc%A, source=A)
  allocate(qfunc%b, source=b)
  call qfunc%init(x)
  
  allocate(ConjugateGradientType :: minimizer)
  write(*,*) 'Conjugate Gradient  first minimization ...  '
  call minimizer%init(max_iter=10000, tol=1.d-6, xx=x)  
  xini = minimizer%xx_ini
  call minimizer%move(x, qfunc)  
  write(*,*) '|<---------------------------------->|'

  write(*,*) 'Conjugate Gradient  second minimization ...  '
  x = xini 
  call minimizer%init(max_iter=10, tol=1.d-6, xx=x)  
  call minimizer%move(x, qfunc)  
  xend  = minimizer%xx_end 
  deallocate(minimizer)

  write(*,*) 'Conjugate Gradient  third minimization ...  '

  allocate(ConjugateGradientExtType :: cgExt)
  x = xend 
  call cgExt%init(max_iter=10000, tol=1.d-6, xx=x)  
  call cgExt%move(x, qfunc)  
  write(*,*) '|<---------------------------------->|'


  write(*,*) 'LBFGS  first minimization ...  '
  allocate(LbfgsExtType :: minimizer)
  x = xend 
  lbfgs_m_hess = 3 
  lbfgs_xtol = 1.d-12
  call minimizer%init(max_iter=10000, tol=1.d-6, xx=x)  
  call minimizer%move(x, qfunc) 
  write(*,*) '|<---------------------------------->|'


! Allocate and initialize steepest descent minimizer
!ifort   -C    -I${MKLROOT}/include/intel64/lp64 -I"${MKLROOT}/include"  module_lbfgs_nocedal.F90 class_functions.F90 maquette.f90   ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl

  
end program MinimizerExample
