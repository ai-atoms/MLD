module MinimizerModule
  implicit none

  public :: MinimizerType
  
  type, abstract :: MinimizerType
    integer :: max_iter
    real(kind=8) :: tol
    real(kind=8), dimension(:), allocatable :: xx_ini    
    real(kind=8), dimension(:), allocatable :: xx_end
    contains 
    procedure(init_shape), deferred :: init 
    procedure(move_shape), deferred :: move      
  end type MinimizerType

  interface
  !
    subroutine init_shape(this,  max_iter, tol, xx)
      import :: MinimizerType
      class(MinimizerType), intent(inout)   :: this
      real(kind=8), dimension(:), intent(inout) :: xx
      integer, intent(in)  :: max_iter 
      real(kind=8), intent(in)  :: tol
    end subroutine init_shape
    !
    subroutine move_shape(this, xx, func)
      use FunctionModule, only : FunctionType
      import :: MinimizerType
      class(MinimizerType), intent(inout)  :: this
      class(FunctionType), intent(inout) :: func 
      real(kind=8), dimension(:), intent(inout) :: xx
    end subroutine move_shape
    !
  end interface 

end module MinimizerModule  

module SteepestDescentModule 
  use MinimizerModule, only: MinimizerType
  implicit none 
  private 
  public ::  SteepestDescentType

  type, extends(MinimizerType) :: SteepestDescentType
    real(kind=8) :: sd_step 
    contains 
    procedure :: init => init_steepest_descent
    procedure :: move => move_steepest_descent
  end type SteepestDescentType

  contains 

  subroutine init_steepest_descent (this, max_iter, tol, xx)
    class(SteepestDescentType), intent(inout) :: this 
    real(kind=8), dimension(:), intent(inout) :: xx
    integer, intent(in)  :: max_iter 
    real(kind=8), intent(in)  :: tol
    
    this%max_iter = max_iter 
    this%tol = tol 
    this%sd_step = 1.d-10  
    if (allocated(this%xx_ini)) deallocate(this%xx_ini)
    allocate(this%xx_ini, source=xx)

  end subroutine init_steepest_descent 

  subroutine move_steepest_descent(this, xx, func)
    use FunctionModule, only: FunctionType
    class(SteepestDescentType), intent(inout) :: this
    class(FunctionType), intent(inout) :: func
    real(kind=8), dimension(:), intent(inout) :: xx
    integer :: ii 
    real(kind=8), dimension(:), allocatable  :: df 
    real(kind=8) :: eta, f_value, norm2

    eta = this%sd_step  
    allocate(df(func%fdim))  
    ! Perform steepest descent optimization
    do ii = 1, this%max_iter
      call func%evaluate(xx, f_value)
      call func%gradient(xx, df)
      xx = xx - eta * df
      if (norm2(df) < this%tol) exit
      if (mod(ii, 1)==0) write(*,'("step value ...: ", i9, 2es25.8)') ii, norm2(df), f_value
    end do  

    if (allocated(this%xx_end)) deallocate(this%xx_end)
    allocate(this%xx_end, source=xx)

    deallocate(df)
  end subroutine move_steepest_descent
  
end module SteepestDescentModule


module ConjugateGradientModule
  use MinimizerModule, only: MinimizerType
  implicit none 
  private 
  public :: ConjugateGradientType 

  type, extends(MinimizerType) :: ConjugateGradientType
    logical  :: precondition  
    real(kind=8) :: alpha, beta 
    contains 
    procedure :: init => init_conjugate_gradient
    procedure :: move => move_conjugate_gradient
  end type ConjugateGradientType

  contains 
  
  subroutine init_conjugate_gradient (this, max_iter, tol, xx)
    class(ConjugateGradientType), intent(inout) :: this 
    real(kind=8), dimension(:), intent(inout) :: xx
    integer, intent(in)  :: max_iter 
    real(kind=8), intent(in)  :: tol
    
    this%max_iter = max_iter 
    this%tol = tol
    !CG part ...  
    this%precondition = .false. 
    this%alpha=0.005d0  
    this%beta=0.9d0   
    if (allocated(this%xx_ini)) deallocate(this%xx_ini)
    allocate(this%xx_ini, source=xx)

  end subroutine init_conjugate_gradient

  subroutine move_conjugate_gradient(this, xx, func)
    use FunctionModule, only: FunctionType
    implicit none 
    class(ConjugateGradientType), intent(inout) :: this
    class(FunctionType), intent(inout) :: func
    real(kind=8), dimension(:), intent(inout) :: xx


    integer :: ii, xdim, iline_search
    real(kind=8) :: alpha, beta, f_curr, f_new   
    real(kind=8), dimension(:), allocatable :: xx_curr, xx_new, xx_temp, & 
                                               grad_curr, grad_new, direction
    ! Initialize x_curr, grad_curr, and direction
    if ( size(xx) /= func%fdim ) then
       write(*,*) 'size of the function and size of x is incompatible', size(xx), func%fdim
       stop 'minimization in conjugate_gradient'  
    end if 
    xdim = size(xx,1)
    allocate(xx_curr(xdim), xx_new(xdim), xx_temp(xdim), &
             grad_curr(xdim), grad_new(xdim),           &
             direction(xdim) ) 
    xx_curr(:) = xx(:)
    !grad_curr = numerical_gradient(funk, x_curr)
    call  func%gradient(xx_curr, grad_curr)
    direction = -grad_curr
    ! Conjugate gradient optimization loop
    do ii = 1, this%max_iter
      ! Compute step size by line search ... 
      ! alpha = line_search(funk, x_curr, grad_curr, direction)
      alpha = this%alpha
      beta = this%beta 
      call func%evaluate(xx_curr, f_curr)
      do iline_search = 1, 10
         xx_temp(:) = xx_curr(:) + alpha * (xx_new(:) - direction(:))
         call func%evaluate(xx_temp, f_new)
         if (f_new < f_curr) exit
         alpha = alpha * beta
      end do
          
      ! Update x_new and grad_new
      xx_new(:) = xx_curr(:) + alpha * direction(:)
      !del grad_new = numerical_gradient(funk, x_new)
      call func%gradient(xx_new, grad_new)  
      call func%evaluate(xx_new, f_new) 
      ! Check convergence
      if (norm2(grad_new) < this%tol ) exit
      if (mod(ii, 10)==0) write(*,'("cg step value ...: ", i9, 2es25.8)') ii, norm2(grad_new), f_new
      ! Update beta and direction
      beta = dot_product(grad_new, grad_new - grad_curr) / dot_product(grad_curr, grad_curr)
      direction(:) = -grad_new(:) + beta * direction(:)
      ! Update x_curr and grad_curr
      xx_curr(:) = xx_new(:)
      grad_curr(:) = grad_new(:)
    end do
    
    write(*,'("cg step value ...: ", i9, 2es25.8)') ii, norm2(grad_new), f_new

    if (allocated(this%xx_end)) deallocate(this%xx_end)
    allocate(this%xx_end, source=xx_curr)

    ! Deallocate memory
    deallocate(xx_curr, xx_new, grad_curr, grad_new, direction)
  end subroutine move_conjugate_gradient

end module ConjugateGradientModule 


module AdamDescentModule 
  use MinimizerModule, only: MinimizerType
  use module_kind_variables, only: kind_double
  implicit none 
  private 
  public ::  AdamDescentType

  type, extends(MinimizerType) :: AdamDescentType
    real(kind_double) :: learning_rate
    real(kind_double) :: beta1, beta2, epsilon
    real(kind_double), dimension(:), allocatable :: m, v
    integer :: t

    contains 
    procedure :: init => init_adam_sdescent
    procedure :: move => move_adam_sdescent
  end type AdamDescentType

  contains 

  subroutine init_adam_sdescent (this, max_iter, tol, xx)
    use mld_logger, only: log_info
    class(AdamDescentType), intent(inout) :: this 
    real(kind=8), dimension(:), intent(inout) :: xx
    integer, intent(in)  :: max_iter 
    real(kind=8), intent(in)  :: tol
    
    this%max_iter = max_iter 
    this%tol = tol 
    this%learning_rate = 1.d-6
    this%beta1 = 0.9d0
    this%beta2 = 0.999d0
    this%epsilon = 1.d-10
    this%t = 0
    if (allocated(this%m)) deallocate(this%m)
    allocate(this%m, source=xx)
    if (allocated(this%v)) deallocate(this%v)
    allocate(this%v, source=xx)  
    if (allocated(this%xx_ini)) deallocate(this%xx_ini)
    allocate(this%xx_ini, source=xx)
    this%xx_ini = xx 
    this%m = 0.d0 
    this%v = 0.d0
    call log_info('Adam minimizer        step           deriv               value')

  end subroutine init_adam_sdescent 

  subroutine update_adam(this, xx, df)
    class(AdamDescentType), intent(inout) :: this
    real(kind_double), dimension(:), intent(inout) :: xx
    real(kind_double), dimension(:), intent(in) :: df
    real(kind_double), dimension(size(df)) :: m_hat, v_hat

    this%t = this%t + 1
    this%m = this%beta1 * this%m + (1.d0 - this%beta1) * df
    this%v = this%beta2 * this%v + (1.d0 - this%beta2) * df**2

    m_hat = this%m / (1.d0 - this%beta1**this%t)
    v_hat = this%v / (1.d0 - this%beta2**this%t)
    xx = xx - this%learning_rate * m_hat / (sqrt(v_hat) + this%epsilon)

  end subroutine update_adam

  subroutine move_adam_sdescent(this, xx, func)
    use FunctionModule, only: FunctionType
    use mld_logger, only: log_info, vtoa
    class(AdamDescentType), intent(inout) :: this
    class(FunctionType), intent(inout) :: func
    real(kind_double), dimension(:), intent(inout) :: xx
    integer :: ii 
    real(kind_double), dimension(:), allocatable  :: df 
    real(kind_double) :: f_value, norm2
    character(len=100) :: log_tmp
    
    allocate(df(func%fdim))  

    do ii = 1, this%max_iter
      call func%evaluate(xx, f_value)
      call func%gradient(xx, df)
      call update_adam(this, xx, df)
      if (norm2(df) < this%tol) exit
      if (mod(ii, 10)==0) then
        !write(*,'("step value ...: ", i9, 2es25.8)') ii, norm2(df), f_value
        write(log_tmp, '(i9, 2es20.8)') ii, norm2(df), f_value
        call log_info('                 '//trim(log_tmp)) 
      end if 
    end do  

    if (allocated(this%xx_end)) deallocate(this%xx_end)
    allocate(this%xx_end, source=xx)
    this%xx_end = xx

    deallocate(df)
  end subroutine move_adam_sdescent
  
end module AdamDescentModule
