module module_random_matrix
  implicit none

contains

  subroutine generate_random_matrix(ndim, mdim, mat)
    integer, intent(in) :: ndim, mdim
    real(kind=8), allocatable, intent(out) :: mat(:,:)
    integer :: i, j

    ! Allocate space for the matrix
    allocate(mat(ndim, mdim))

    ! Generate random matrix
    call random_seed()
    do i = 1, ndim
      do j = 1, mdim
        call random_number(mat(i, j))
      end do
    end do
  end subroutine generate_random_matrix

  subroutine generate_spd_matrix(ndim, SPD)
    integer, intent(in) :: ndim
    real(kind=8), allocatable, intent(out) :: SPD(:,:)
    real(kind=8) :: A(ndim, ndim)
    integer :: i, j

    ! Generate random matrix A
    call random_seed()
    do i = 1, ndim
      do j = 1, ndim
        call random_number(A(i, j))
      end do
    end do

    ! Allocate space for the SPD matrix
    allocate(SPD(ndim, ndim))

    ! Create a symmetric positive definite matrix: SPD = A' * A
    SPD = matmul(transpose(A), A)
  end subroutine generate_spd_matrix

  subroutine generate_random_vector(vdim, vec)
    integer, intent(in) :: vdim
    real(kind=8), allocatable, intent(out) :: vec(:)
    integer :: i

    ! Allocate space for the vector
    allocate(vec(vdim))

    ! Generate random vector
    call random_seed()
    do i = 1, vdim
      call random_number(vec(i))
    end do
  end subroutine generate_random_vector

end module module_random_matrix


module FunctionModule
  implicit none
  private
  public :: FunctionType
  
  type, abstract  :: FunctionType
    ! Abstract methods
    integer :: fdim 
    real(kind=8), dimension(:),  allocatable :: xinit 
    real(kind=8), dimension(:),  allocatable :: params
    contains 
    procedure(init_shape)    , deferred   :: init
    procedure(evaluate_shape), deferred   :: evaluate 
    procedure(gradient_shape), deferred   :: gradient
    procedure(hessian_shape) , deferred   :: hessian
  end type FunctionType

  interface
  !
    subroutine init_shape(this, xx)
      import :: FunctionType
      class(FunctionType), intent(inout)   :: this
      real(kind=8), dimension(:), intent(inout) :: xx
    end subroutine init_shape
    !
    subroutine evaluate_shape(this, xx, fvalue)
      import :: FunctionType
      class(FunctionType), intent(inout)  :: this
      real(kind=8), dimension(:), intent(in) :: xx
      real(kind=8), intent(out) :: fvalue
    end subroutine evaluate_shape
    !
    subroutine gradient_shape(this, xx, fgradient)
      import :: FunctionType
      class(FunctionType), intent(inout)  :: this
      real(kind=8), dimension(:), intent(in) :: xx
      real(kind=8), dimension(:), intent(out) :: fgradient
    end subroutine gradient_shape
    !
    subroutine hessian_shape(this, xx, fhessian)
      import :: FunctionType
      class(FunctionType), intent(inout)  :: this
      real(kind=8), dimension(:), intent(in) :: xx
      real(kind=8), dimension(:,:), intent(out) :: fhessian
    end subroutine hessian_shape
    !
  end interface

end module FunctionModule


module FuncQuadraticModule 
  use FunctionModule, only: FunctionType
  implicit none
  private 
  public :: FuncQuadraticType 

  type, extends(FunctionType) :: FuncQuadraticType
    real(kind=8), dimension(:,:), allocatable :: A
    real(kind=8), dimension(:), allocatable   :: b
  contains
    procedure :: init
    procedure :: evaluate
    procedure :: gradient 
    procedure :: hessian  
  end type  FuncQuadraticType
  
  contains 

  ! Quadratic test function .... 
  subroutine init(this, xx)
    class(FuncQuadraticType), intent(inout) :: this
    real(kind=8), dimension(:), intent(inout) :: xx 
    
    allocate(this%xinit, source=xx)
    this%xinit = xx
    this%fdim = size(xx,1) 
      if (.not.(allocated(this%A))) then 
        write(*,*) 'A matrix is not initialized for FuncQuadraticType'
        stop 'stop for A in FuncQuyadraticModule init ' 
    end if 

    if (.not.(allocated(this%b))) then 
        write(*,*) 'b vector is not initialized for FuncQuadraticType'
        stop 'stop for b in FuncQuadraticModule init' 
    end if 

  end subroutine init
  
  subroutine evaluate(this, xx, fvalue)
    class(FuncQuadraticType), intent(inout) :: this
    real(kind=8), dimension(:), intent(in) ::  xx
    real(kind=8), intent(out) :: fvalue
    integer :: i
  
    fvalue = 0.d0
    do i = 1, this%fdim
      fvalue = fvalue + 0.5d0 * xx(i) * (DOT_PRODUCT(this%A(:,i), xx(:)) - this%b(i))
    end do
  end subroutine evaluate
  ! 
  subroutine gradient(this, xx, fgradient)
    class(FuncQuadraticType), intent(inout)    :: this
    real(kind=8), dimension(:), intent(in)  :: xx
    real(kind=8), dimension(:), intent(out) :: fgradient
    integer :: i, j

    fgradient(:) = 0.d0
    do i = 1, this%fdim
      do j = 1, this%fdim
        fgradient(i) = fgradient(i) + this%A(i,j) * xx(j)
      end do
      fgradient(i) = fgradient(i) - this%b(i)
    end do
  end subroutine gradient
  ! 
  subroutine hessian(this, xx, fhessian)
    class(FuncQuadraticType), intent(inout) :: this
    real(kind=8), dimension(:), intent(in) :: xx
    real(kind=8), dimension(:,:), intent(out) :: fhessian
    integer :: i, j
    real(kind=8) :: fvalue

    fvalue = sum(xx) 
    
    do i = 1, this%fdim
      do j = 1, this%fdim
        fhessian(i,j) = this%A(i,j)
      end do
    end do
  end subroutine hessian

end module FuncQuadraticModule 


module FuncCubicModule 
  use FunctionModule, only: FunctionType
  implicit none
  private 
  public :: FuncCubicType 

  type, extends(FunctionType) :: FuncCubicType
    real(kind=8), dimension(:,:), allocatable :: A
    real(kind=8), dimension(:), allocatable :: b
    real(kind=8), dimension(:), allocatable :: c
  contains
    procedure :: init
    procedure :: evaluate
    procedure :: gradient
    procedure :: hessian
  end type FuncCubicType

  contains 

  !   init Cubic test function ...; 
  subroutine init(this, xx)
    class(FuncCubicType), intent(inout) :: this
    real(kind=8), dimension(:), intent(inout) :: xx

    allocate(this%xinit, source=xx)
    this%xinit = xx
    this%fdim = size(xx,1) 
    if (.not.(allocated(this%A))) then 
        write(*,*) 'A matrix is not initialized for FuncCubicType'
        stop 'stop for A in FuncCubicModule init ' 
    end if 

    if (.not.(allocated(this%b))) then 
        write(*,*) 'b vector is not initialized for FuncCubicType'
        stop 'stop for b in FuncCubicModule init' 
    end if 

    if (.not.(allocated(this%c))) then 
        write(*,*) 'c vector  is not initialized for FuncCubicType'
        stop 'stop for c in FuncCubicModule init' 
    end if 

  end subroutine init
  
  ! evaluate Cubic test function ...; 
  subroutine evaluate(this, xx, fvalue)
    class(FuncCubicType), intent(inout) :: this
    real(kind=8), dimension(:), intent(in) :: xx
    real(kind=8), intent(out) :: fvalue
    integer :: i
    
    fvalue = 0.d0
    do i = 1, this%fdim
      fvalue = fvalue + this%c(i) * xx(i)**3
    end do
  end subroutine evaluate
  
  subroutine gradient(this, xx, fgradient)
    class(FuncCubicType), intent(inout) :: this
    real(kind=8), dimension(:), intent(in) :: xx
    real(kind=8), dimension(:), intent(out) :: fgradient
    integer :: i
    
    fgradient(:) = 0.d0 
    do i = 1, this%fdim
      fgradient(i) = 3.0 * this%c(i) * xx(i)**2
    end do
  end subroutine gradient
  
  subroutine hessian(this, xx, fhessian)
    class(FuncCubicType), intent(inout) :: this
    real(kind=8), dimension(:), intent(in) :: xx
    real(kind=8), dimension(:,:), intent(out) :: fhessian
    integer :: i
    
    fhessian(:,:) = 0.d0 
    do i = 1, this%fdim
      fhessian(i,i) = 6.0 * this%c(i) * xx(i)
    end do
  end subroutine hessian

end module FuncCubicModule 


