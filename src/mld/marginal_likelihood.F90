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
#include "../MLD_MACROS.INC"

module module_MargLikeMAP0 
  use mld_logger, only: vtoa, log_info, mld_verbose, log_critical, log_debug 
  use FunctionModule, only: FunctionType
  use ml_in_ndm_module, only: one_pi 
  use module_ml_scalapack, only: scalapack_driver, dimc_sca_Amat, dimr_sca_Amat

  implicit none
  private 
  public :: FuncMargLikeMAP0 

  type, extends(FunctionType) :: FuncMargLikeMAP0
    real(kind=8), dimension(:), allocatable :: wparams 

    integer :: M   ! number of data in design matrix 
    integer :: P   ! size of the problem D+1 for LML, D2+D+1 for others etc  
    ! The covariance matrix for data and parameters. We consider only the diagonal terms.  
    real(kind=8), dimension(:), allocatable :: iSigmab_M
    real(kind=8), dimension(:), allocatable :: iSigmab_P
    ! The mean vector for data and parameters.
    real(kind=8), dimension(:), allocatable :: mub_P

    !for fast mpi computing: 
    integer :: start_mm, end_mm
  contains
    procedure :: init
    procedure :: evaluate
    procedure :: gradient 
    procedure :: hessian  
  end type  FuncMargLikeMAP0
  
  contains 

  ! Quadratic test function .... 
  subroutine init(this, xx)
    use ml_in_ndm_module, only: rangml
    use snap, only: w_params, Amat
    use mld_mpi, only: mld_mpi_abort 
    use set_limits, only: set_limit_for_configs
    implicit none 
    class(FuncMargLikeMAP0), intent(inout) :: this
    real(kind=8), dimension(:), intent(inout) :: xx 
    integer, dimension(:), allocatable ::  mm_start_on_proc, mm_end_on_proc, &
                                           mm_size_on_proc, mm_to_proc
    _NAMECURRENT_("FuncMargLikeMAP0.init")

    _MLD_BEGIN_

    if (allocated(this%xinit)) deallocate(this%xinit)
    allocate(this%xinit, source=xx)
    
    ! check if w_params are allocated 
    if (.not.(allocated(w_params))) then 
      call log_critical('w_params are not initialized in '//NAMECURRENT)
      call mld_mpi_abort('stop for wparams in '//NAMECURRENT) 
    end if 

    if (allocated(this%wparams)) deallocate(this%wparams)
    allocate(this%wparams, source=w_params(:,1))
    this%xinit = xx
    this%fdim = size(xx,1) 
    

    ! for the next version ... 
    !this%params = (wparams, b, A)
    !  this%fdim = dim(xx) + dim(A) + dim()

    !un! if (this%fdim /= size(this%wparams)) then 
    !un!   call log_critical('wparams and fdim have not the same dimension  '// vtoa(size(this%wparams)) //"  "// vtoa(this%fdim))
    !un!   call log_critical('MLD will stop in '//NAMECURRENT)
    !un!   call mld_mpi_abort('stop for size fdim/wparams in '//NAMECURRENT)  
    !un! end if 
    
    if (scalapack_driver) then 
      this%M = dimc_sca_Amat
      this%P = dimr_sca_Amat 
    else 
      this%P = size(Amat,1)
      this%M = size(Amat,2)
    end if 

    if (this%fdim /= this%P )  then 
      call log_critical('P and fdim arenot equal '// vtoa(this%P) //"  "// vtoa(this%fdim))
      call log_critical('MLD will stop in '//NAMECURRENT)
      call mld_mpi_abort('stop for size fdim/wparams in '//NAMECURRENT)  
    end if 


    if (allocated(this%iSigmab_M)) deallocate(this%iSigmab_M) ; allocate(this%iSigmab_M(this%M))
    if (allocated(this%iSigmab_P)) deallocate(this%iSigmab_P) ; allocate(this%iSigmab_P(this%P))
    if (allocated(this%mub_P)) deallocate(this%mub_P) ; allocate(this%mub_P(this%P))

    this%iSigmab_M(:) = 1.d0 ! inverse of matrix, i.e. 1/sigmaM**2 
    this%iSigmab_P(:) = 0.1d0 ! inverse of matrix, i.e. 1/sigmaP**2
    !this%mub_P(:) = w_params(:,1) 
    this%mub_P(:) = 0.d0  

    ! for fast mpi computing, set the limits with hacked subroutine for config. But thart can reused here: 
    call set_limit_for_configs(rangml, this%M, this%start_mm, this%end_mm, &
                                mm_start_on_proc, mm_end_on_proc, mm_size_on_proc, mm_to_proc)


    _MLD_END_
  end subroutine init
  
  subroutine evaluate(this, xx, fvalue)
    use snap, only: ymat, Amat  
    use module_kind_variables, only: kind_double
    class(FuncMargLikeMAP0), intent(inout) :: this
    real(kind_double), dimension(:), intent(in) ::  xx
    real(kind_double), intent(out) :: fvalue
    real(kind_double) :: t1, t2, t3, t4, t5, t6, rnorm 
    real(kind_double), dimension(:), allocatable :: y_copy, vec_tmp 
    real(kind_double), external :: dnrm2
    
    _NAMECURRENT_("FuncMargLikeMAP0.evaluate")
    ! in this model we find just the ML parameters i.e. we find the extremum of the 
    ! posterior probability: p(w | X , Y ) = p(Y | w, X) * p(w). 
    ! We are in MAP model (i.e. the normalizing marginal likelihood is taken 
    ! as a parameter. 
    _MLD_BEGIN_
    this%wparams = xx
    ! There are 6 terms to evaluate: 
    ! 
    t1 = 0.5d0 * dble(this%M)  * dlog (2.d0 * one_pi)
    t1 = 0.d0 

    t2 = 0.5d0 * sum(dlog(this%iSigmab_M(:)))   
    
    if (scalapack_driver) then
      t3 = 0.d0 !to complete E(\w)    
    else 
      allocate(vec_tmp, source=xx)
      vec_tmp = xx - this%mub_P
      rnorm = dnrm2(this%P, vec_tmp,1)
      t3 = 0.5d0 * rnorm**2 !to complete E(\w)    
    end if 
    
    t4 = 0.5d0 * dble(this%P)  * dlog (2.d0 * one_pi)
    t4 = 0.d0
    t5 = 0.5d0 * sum(dlog(this%iSigmab_P(:)))   
    
    if (scalapack_driver) then
      t6 = 0.d0 ! to complete e(\w)   
    else 
      if (allocated(y_copy)) deallocate(y_copy)
      allocate(y_copy, source=ymat(:,1))
      y_copy(:) = ymat(:,1)
      ! y_copy = y_copy - Amat * wparams 
      !          M x 1  - M x P * P x 1
      !TODO there are missing weights here ...
      call dgemv('T', this%P, this%M, -1.d0, Amat, this%P, xx, 1, 1.d0, y_copy, 1) 
      rnorm = dnrm2(this%M, y_copy, 1)
      t6 = 0.5d0 * rnorm**2
    end if 

    fvalue = t1 + t2 + t3 + t4 + t5 + t6
    !debug! if (rangml == 0 ) then
    !   write(*,'(a, 6es15.7)') 't1 t2 t3 t4 t5 t6  = ', t1, t2, t3, t4, t5, t6  
    !debug! end if 
    _MLD_END_ 
  end subroutine evaluate
  ! 
  subroutine gradient(this, xx, fgradient)
    use module_kind_variables, only: kind_double
    use snap, only: ymat, Amat
    !use mpi 
    class(FuncMargLikeMAP0), intent(inout)    :: this
    real(kind_double), dimension(:), intent(in)  :: xx
    real(kind_double), dimension(:), intent(out) :: fgradient
    integer :: mm, pp 
    real(kind_double) :: tmp_mm 
    real(kind_double), dimension(:), allocatable :: l1_gradient, l2_gradient, l3_gradient 
    _NAMECURRENT_("FuncMargLikeMAP0.gradient")
    _MLD_BEGIN_ 
    
    fgradient(1:this%fdim) = 0.d0
    allocate(l1_gradient, source=fgradient)
    allocate(l2_gradient, source=fgradient)
    allocate(l3_gradient, source=fgradient)

    ! there two terms each of them with two componet 
    ! t1 = t1a + t1b 
    ! t1 =  - ymat^T * \Sigma^{-1} * Amat^T + ww^T * Amat * \Sigma^{-1} * Amat^T 
    ! 1XP =  1 x M  * MxM * MxP             1xP * PxM  * MxM * MxP
    ! we will compute directly the transpose: 
    ! t1 =  - Amat * \Sigma^{-1} * ymat + Amat * \Sigma^{-1} * Amat^T * ww
    ! Px1 =   PxM  *   MxM       * Mx1     PxM * MxM * MxP * Px1
    if (scalapack_driver) then
      fgradient(:) = 0.d0 ! to complete  
    else 
      !t1a: \sum_m Amat_pm * (Sigma_M)_m * ymat_m  ! is Sigma_M diagonal 
      l1_gradient(:) = 0.d0 
      do mm = 1, this%M
        l1_gradient(:) = l1_gradient(:) - Amat(:,mm) * this%iSigmab_M(mm) * ymat(mm,1)
      end do
      fgradient(:) = fgradient(:) + l1_gradient(:)

      !fast! local_gradient(:)=0.d0 
      !fast! do mm = this%start_mm, this%end_mm
      !fast!   local_gradient(:) = local_gradient(:) - Amat(:,mm) * this%iSigmab_M(mm) * ymat(mm,1)
      !fast! end do
      !fast! call MPI_ALLREDUCE(MPI_IN_PLACE, local_gradient, this%P, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, codeml)
      !fast! fgradient(:) = fgradient(:) + local_gradient(:)


      !t1b: \sum_m \sum_p' Amat_pm * (Sigma_M)_m * Amat_p'm * wparams_p'  ! is Sigma_M diagonal

      l2_gradient(:)=0.d0
      do mm = 1, this%M
        tmp_mm  = 0.d0 
        do pp = 1, this%P
          tmp_mm = tmp_mm  +  Amat(pp,mm) * xx(pp) 
        end do  
        l2_gradient(:) = l2_gradient(:) + Amat(:,mm) * this%iSigmab_M(mm) * tmp_mm  
      end do 
      fgradient(:) = fgradient(:) + l2_gradient(:)
      
      !fast! local_gradient(:)=0.d0
      !fast! do mm = 1, this%start_mm, this%end_mm
      !fast!   tmp_mm  = 0.d0 
      !fast!   do pp = 1, this%P
      !fast!     tmp_mm = tmp_mm  +  Amat(pp,mm) * xx(pp) 
      !fast!   end do  
      !fast!   local_gradient(:) = local_gradient(:) + Amat(:,mm) * this%iSigmab_M(mm) * tmp_mm 
      !fast! end do 
      !fast! call MPI_ALLREDUCE(MPI_IN_PLACE, local_gradient, this%P, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, codeml)
      !fast! fgradient(:) = fgradient(:) + local_gradient(:)

      !t2:  (xx^T - mub) 
      do pp = 1, this%P 
          l3_gradient(pp) = (xx(pp) - this%mub_P(pp))*this%iSigmab_P(pp)
      end do 
      fgradient(:) = fgradient(:) + l3_gradient(:) 
    end if 

    !debug! write(*,'(a, 3es15.7 )') '....in gradient ....', norm2(fgradient), &
    !debug!      norm2(l1_gradient + l2_gradient), norm2(l3_gradient)
    
    _MLD_END_ 

  end subroutine gradient
  ! 
  subroutine hessian(this, xx, fhessian)
    use mld_mpi, only: mld_mpi_abort 
    class(FuncMargLikeMAP0), intent(inout) :: this
    real(kind=8), dimension(:), intent(in) :: xx
    real(kind=8), dimension(:,:), intent(out) :: fhessian
    integer :: pp 

    _NAMECURRENT_("FuncMargLikeMAP0.hessian")
    _MLD_BEGIN_
    call log_critical('Hessian is not implemented for FuncMargLikeMAP0')
    call mld_mpi_abort('stop for hessian in '//NAMECURRENT)
    pp = this%P 
    fhessian(:,:) = 0.d0*sum(xx) 
    _MLD_END_
  end subroutine hessian

end module module_MargLikeMAP0


module module_MargLikeAlphaBeta 
  use, intrinsic :: iso_fortran_env, dp=>real64
  use mld_logger, only: vtoa, log_info, mld_verbose, log_critical, log_debug 
  use FunctionModule, only: FunctionType
  use ml_in_ndm_module, only: one_pi 
  use module_ml_scalapack, only: scalapack_driver, dimc_sca_Amat, dimr_sca_Amat

  implicit none
  private 
  public :: FuncMargLikeAlphaBeta 

  type, extends(FunctionType) :: FuncMargLikeAlphaBeta
    real(kind=8), dimension(:), allocatable :: xx 

    integer :: M   ! number of data in design matrix 
    integer :: P   ! size of the problem D+1 for LML, D2+D+1 for others etc  
    integer :: dim_margl, dim_Pside, dim_Mside  ! number of MargLike parameters to be optimized
    ! The covariance matrix for data and parameters. We consider only the diagonal terms.  
    real(dp), dimension(:), allocatable :: iSigmab_M
    real(dp), dimension(:), allocatable :: iSigmab_P
    ! The mean vector for data and parameters.
    real(dp), dimension(:), allocatable :: mub_P
    real(dp), dimension(:), allocatable :: mP
    ! Smat and SmatI 
    real(dp), dimension(:,:), allocatable :: Smat, SmatI

    !for fast mpi computing: 
    integer :: start_mm, end_mm
  contains
    procedure :: init
    procedure :: evaluate
    procedure :: gradient 
    procedure :: hessian  
  end type  FuncMargLikeAlphaBeta
  
  contains 

  ! Quadratic test function .... 
  subroutine init(this, xx)
    use ml_in_ndm_module, only: rangml
    use snap, only:  Amat, fit_snap
    use mld_mpi, only: mld_mpi_abort 
    use set_limits, only: set_limit_for_configs
    implicit none 
    class(FuncMargLikeAlphaBeta), intent(inout) :: this
    real(kind=8), dimension(:), intent(inout) :: xx 
    integer, dimension(:), allocatable ::  mm_start_on_proc, mm_end_on_proc, &
                                           mm_size_on_proc, mm_to_proc
    _NAMECURRENT_("FuncMargLikeAlphaBeta.init")

    _MLD_BEGIN_


    if (scalapack_driver) then 
      this%M = dimc_sca_Amat
      this%P = dimr_sca_Amat
    else 
      this%P = size(Amat,1)
      this%M = size(Amat,2)
    end if 

    this%dim_Pside = 2*this%P
    this%dim_Mside = 3 
    this%dim_margl = this%dim_Pside + this%dim_Mside

    this%fdim = this%dim_margl 
    if (this%fdim /= size(xx)) then 
      call log_critical('xx and fdim have not the same dimension  '// vtoa(size(xx)) //"  "// vtoa(this%fdim))
      call mld_mpi_abort('stop for size fdim/xx in '//NAMECURRENT)  
    end if

    if (allocated(this%iSigmab_M)) deallocate(this%iSigmab_M) ; allocate(this%iSigmab_M(this%M))
    if (allocated(this%iSigmab_P)) deallocate(this%iSigmab_P) ; allocate(this%iSigmab_P(this%P))
    if (allocated(this%mub_P)) deallocate(this%mub_P) ; allocate(this%mub_P(this%P))
    if (allocated(this%mP)) deallocate(this%mP) ; allocate(this%mP(this%P))
    if (allocated(this%Smat)) deallocate(this%Smat) ; allocate(this%Smat(this%P, this%P))
    if (allocated(this%SmatI)) deallocate(this%SmatI) ; allocate(this%SmatI(this%P, this%P))

    this%iSigmab_M(:) = 1.d0 ! inverse of matrix, i.e. 1/sigmaM**2 
    this%iSigmab_P(:) = 0.1d0 ! inverse of matrix, i.e. 1/sigmaP**2
    this%mub_P(:) = 0.d0  
    
    if (allocated(this%xinit)) deallocate(this%xinit)
    allocate(this%xinit, source=xx) 
    allocate(this%xx, source=xx) 
    this%xinit = xx
    this%xx = xx 
    !$! xx(1:this%P) = this%iSigmab_P(:)
    !$! xx(this%P+1:2*this%P) = this%mub_P(:)
    !$! xx(2*this%P+1:3*this%P) = 0.01d0
    ! for fast mpi computing, set the limits with hacked subroutine for config. But thart can reused here: 
    call set_limit_for_configs(rangml, this%M, this%start_mm, this%end_mm, &
                                mm_start_on_proc, mm_end_on_proc, mm_size_on_proc, mm_to_proc)

    if (size(fit_snap) /= size(Amat,2)) then 
      call log_critical('fit_snap and Amat have not the same dimension  '// vtoa(size(fit_snap)) //"  "// vtoa(size(Amat,2)))
      call mld_mpi_abort('stop for size fit_snap/Amat in '//NAMECURRENT)  
    end if


    _MLD_END_
  end subroutine init

  subroutine deploy_sigma_M_P(this)
    use, intrinsic :: iso_fortran_env, dp=>real64
    use snap, only : fit_snap, Amat 
    class(FuncMargLikeAlphaBeta), intent(inout) :: this
    integer :: ii, pp 

    do ii = 1, size(Amat,2) 
      if (fit_snap(ii)%energy) this%iSigmab_M(ii) = this%xx(2*this%P+1)
      if (fit_snap(ii)%force)  this%iSigmab_M(ii) = this%xx(2*this%P+2)
      if (fit_snap(ii)%stress) this%iSigmab_M(ii) = this%xx(2*this%P+3)
    end do  


    do pp = 1, this%P 
      this%iSigmab_P(pp) = this%xx(pp)
      this%mub_P(pp) = this%xx(this%P+pp)
    end do

  end subroutine deploy_sigma_M_P


  
  subroutine evaluate(this, xx, fvalue)
    use, intrinsic :: iso_fortran_env, dp=>real64
    use snap, only: ymat, Amat  
    use module_kind_variables, only: kind_double
    use math, only: serial_determinant_symmetric_general, serial_determinant_symmetric_positive
    use mld_mpi, only: mld_mpi_abort
    class(FuncMargLikeAlphaBeta), intent(inout) :: this
    real(kind_double), dimension(:), intent(in) ::  xx
    real(kind_double), intent(out) :: fvalue
    real(dp), dimension(:,:), allocatable  :: Bmat, SmatI 
    real(dp), dimension(:), allocatable :: vP, mP
    integer :: mm, pp, info 
    real(dp) :: tmp, ddot 
    
    _NAMECURRENT_("FuncMargLikeAlphaBeta.evaluate")
    _MLD_BEGIN_

    this%xx = xx
    call deploy_sigma_M_P(this)

    fvalue=0.0_dp

    !T1: -y^T * Sigma_M * y
    !TODO_para ... 
    tmp = 0.0_dp 
    do mm = 1, this%M
      tmp = tmp - ymat(mm,1) * this%iSigmab_M(mm) * ymat(mm,1)
    end do
    fvalue = fvalue + tmp

    !T4: log det Sigma_M^-1 + log det Sigma_P^-1 - Mlog(2 pi) 
    tmp = 0.0_dp
    do mm = 1, this%M
      tmp = tmp + dlog(this%iSigmab_M(mm))
    end do

    do pp = 1, this%P
      tmp = tmp + dlog(this%iSigmab_P(pp))
    end do
    fvalue = fvalue + tmp - dble(this%M) * dlog(2.0_dp * one_pi)

    ! T3: - log det(S_P) = log det(S_P^-1)
    !  S_P^{-1} =  Amat^T * iSigma_M * Amat + iSigma_P
    !TODO_para 
    allocate(Bmat(this%P, this%M))
    Bmat(:,:) = 0.0_dp
    do mm = 1, this%M
      Bmat(:,mm) = Bmat(:,mm) + Amat(:,mm) * this%iSigmab_M(mm)
    end do
    ! S_P^{-1} =   Bmat * Amat + iSigma_P
    ! Cmat = Bmat * Amat
    allocate(SmatI(this%P, this%P))
    SmatI(:,:)=0.0_dp
    call dgemm('N', 'T', this%P, this%P, this%M, 1.0_dp, Bmat, this%P, Amat, this%P, 0.0_dp, SmatI, this%P)
    do pp = 1, this%P
      SmatI(pp,pp) = SmatI(pp,pp) + this%iSigmab_P(pp)
    end do
    call serial_determinant_symmetric_general(SmatI, tmp)
    fvalue = fvalue + tmp
    !T2:  m_p = S_P * v_P
    !T2:  v_P = Amat^T * iSigma_M * y + iSigma_P * mub_P
    !TODO_para 
    allocate(vP(this%P))
    allocate(mP(this%P))
    vP(:) = 0.0_dp
    do mm = 1, this%M
      vP(:) = vP(:) + Amat(:,mm) * this%iSigmab_M(mm) * ymat(mm,1)
    end do
    do pp = 1, this%P
      vP(pp) = vP(pp) + this%iSigmab_P(pp) * this%mub_P(pp)
    end do

    ! we will try to have  m_p=S_p v_p as SmatI^-1 * v_p by Cholesky decomposition: 
    call dpotrf('U', this%P, SmatI, this%P, info)
    if (info /= 0) then 
      call log_critical('Cholesky decomposition failed in '//NAMECURRENT)
      call mld_mpi_abort('I stop for Cholesky decomposition in '//NAMECURRENT)
    end if
    call dpotrs('U', this%P, 1, SmatI, this%P, mP, this%P, info)
    if (info /= 0) then 
      call log_critical('Solving linear by Cholesky decomposition failed in '//NAMECURRENT)
      call mld_mpi_abort('II stop for Cholesky decomposition in '//NAMECURRENT)
    end if
    this%mP(:) = mP(:)
    tmp = ddot(this%P, vP, 1, mP, 1)
    fvalue = fvalue + tmp

    _MLD_END_ 
  end subroutine evaluate
  ! 
  subroutine gradient(this, xx, fgradient)
    use iso_fortran_env, dp=>real64
    use snap, only: ymat, Amat, fit_snap
    use mld_mpi, only: mld_mpi_abort
    !use mpi 
    class(FuncMargLikeAlphaBeta), intent(inout)    :: this
    real(dp), dimension(:), intent(in)  :: xx
    real(dp), dimension(:), intent(out) :: fgradient
    integer :: mm, pp 
    real(dp) ::  tmp_e, tmp_f, tmp_s, tmp  
    integer :: lwork, info, ddot 
    real(dp), dimension(:), allocatable :: work, WW, tmpvec
    real(dp), dimension(:), allocatable :: l1_gradient
    _NAMECURRENT_("FuncMargLikeAlphaBeta.gradient")
    _MLD_BEGIN_ 
    
    fgradient(1:this%fdim) = 0.d0
    allocate(l1_gradient, source=fgradient)

    if (scalapack_driver) then
      fgradient(:) = 0.d0 ! to complete  
      call log_critical('gradient is not implemented for scalapack driver in '//NAMECURRENT)
      call mld_mpi_abort('stop for gradient in '//NAMECURRENT)
    else 
      fgradient(:) = 0.d0
      !T1: -y^T * Sigma_M * y
      !l1_gradient(:) = 0.d0
      l1_gradient(1:this%P) = 0.0_dp
      l1_gradient(this%P+1:2*this%P) = 0.0_dp
      tmp_e = 0.0_dp
      tmp_f = 0.0_dp
      tmp_s = 0.0_dp
      !TODO_para ...
      do mm = 1, this%M
        tmp = - ymat(mm,1) * this%iSigmab_M(mm) * ymat(mm,1)
        if (fit_snap(mm)%energy) tmp_e = tmp_e  + tmp 
        if (fit_snap(mm)%force)  tmp_f = tmp_f  + tmp 
        if (fit_snap(mm)%stress) tmp_s = tmp_s  + tmp 
      end do
      l1_gradient(2*this%P+1) = tmp_e 
      l1_gradient(2*this%P+2) = tmp_f 
      l1_gradient(2*this%P+3) = tmp_s 

      fgradient(:) = fgradient(:) + l1_gradient(:)

      !T2: 
      l1_gradient(1:this%P) = 2.0_dp * xx(this%P+1:2*this%P) * this%mP(1:this%P) - this%mP(1:this%P)**2 
      l1_gradient(this%P+1:2*this%P) = 2.0_dp*xx(1:this%P)*this%mP(1:this%P)
      tmp_e=0.0_dp
      tmp_f=0.0_dp
      tmp_s=0.0_dp
      !TODO_para ...
      do mm=1, this%M
        tmp = ddot(this%P, Amat(:,mm), 1, this%mP, 1)
        if (fit_snap(mm)%energy) tmp_e = tmp_e + 2.0_dp * ymat(mm,1)*tmp - tmp**2
        if (fit_snap(mm)%force)  tmp_f = tmp_f + 2.0_dp * ymat(mm,1)*tmp - tmp**2
        if (fit_snap(mm)%stress) tmp_s = tmp_s + 2.0_dp * ymat(mm,1)*tmp - tmp**2
      end do 
      l1_gradient(2*this%P+1) = tmp_e
      l1_gradient(2*this%P+2) = tmp_f
      l1_gradient(2*this%P+3) = tmp_s
      fgradient(:) = fgradient(:) + l1_gradient(:)

      !T3:
      l1_gradient(1:this%P) = 1.0_dp
      l1_gradient(this%P+1:2*this%P) = 0.0_dp
      
      tmp_e=0.0_dp
      tmp_f=0.0_dp
      tmp_s=0.0_dp
      !TODO_para ...
      do mm =1, this%M
        tmpvec = Amat(:,mm)
        tmp = ddot(this%P, tmpvec, 1, tmpvec, 1) 
        if (fit_snap(mm)%energy) tmp_e = tmp_e + tmp
        if (fit_snap(mm)%force)  tmp_f = tmp_f + tmp
        if (fit_snap(mm)%stress) tmp_s = tmp_s + tmp
      end do 
      lwork=-1
      allocate(WW(this%P))
      allocate(work(1))
      call dsyev('N', 'U', this%P, this%SmatI, this%P, WW, work, lwork, info)
      lwork = int(work(1)) + 1
      deallocate(work)
      allocate(work(lwork))
      call dsyev('N', 'U', this%P, this%SmatI, this%P, WW, work, lwork, info)
      if (info /= 0) then 
        call log_critical('Eigenvalues computation failed in '//NAMECURRENT)
        call mld_mpi_abort('stop for eigenvalues in '//NAMECURRENT)
      end if
      tmp = 1.0_dp 
      do pp = 1, this%P
        tmp = tmp / WW(pp)
      end do
      !T4:
      !T1: -y^T * Sigma_M^1 * y
      l1_gradient(:) = 0.d0
      l1_gradient(1:this%P) = 0.0_dp
      l1_gradient(this%P+1:2*this%P) = 0.0_dp
      l1_gradient(2*this%P+1:3) = 1.0_dp/xx(2*this%P+1:3)
      
      fgradient(:) = fgradient(:) + l1_gradient(:)

    end if 

    _MLD_END_ 

  end subroutine gradient
  ! 
  subroutine hessian(this, xx, fhessian)
    use mld_mpi, only: mld_mpi_abort 
    class(FuncMargLikeAlphaBeta), intent(inout) :: this
    real(kind=8), dimension(:), intent(in) :: xx
    real(kind=8), dimension(:,:), intent(out) :: fhessian
    integer :: pp 

    _NAMECURRENT_("FuncMargLikeAlphaBeta.hessian")
    _MLD_BEGIN_
    call log_critical('Hessian is not implemented for FuncMargLikeMAP0')
    call mld_mpi_abort('stop for hessian in '//NAMECURRENT)
    pp = this%P
    fhessian(:,:) = 0.d0*sum(xx)
    _MLD_END_
  end subroutine hessian

end module module_MargLikeAlphaBeta





subroutine main_train_marginal_likelihood
  use iso_fortran_env, dp=>real64
  !use module_kind_variables, only: kind_double
  use mld_logger, only: vtoa, log_info, mld_verbose, log_critical, log_debug 
  use MinimizerModule, only: MinimizerType
  use SteepestDescentModule, only: SteepestDescentType
  use ConjugateGradientModule, only : ConjugateGradientType
  use ConjugateGradientExtModule, only : ConjugateGradientExtType
  use AdamDescentModule, only: AdamDescentType
  use LbgfsExtModule, only: LbfgsExtType
  use FunctionModule, only: FunctionType
  use module_MargLikeMAP0, only: FuncMargLikeMAP0
  use mld_mpi, only: mld_mpi_abort
  use snap, only: w_params, Amat 
  use main_mld_mod, only : main_train_mld
  use module_condition_number, only: compute_condition_number_mat_rss
  use module_lbfgs_nocedal, only : lbfgs_m_hess, lbfgs_xtol
  use module_ml_scalapack, only: scalapack_driver
  use module_MargLikeAlphaBeta, only: FuncMargLikeAlphaBeta
  implicit none
  real(dp), dimension(:), allocatable :: xini, xend, xx 
  class(FunctionType), allocatable :: func, funcAB
  class(MinimizerType), allocatable :: minimizer
  real(dp), dimension(:,:), allocatable :: symmat 
  real(dp) :: cond 
  character(len=80) :: ctmp 
  logical :: lMAP0, lAlphaBeta


  _NAMECURRENT_("main_train_marginal_likelihood")
  _MLD_BEGIN_

  ! ----> if e want precondition .... 
  ! this will build Amat and  Amat/sca_Amat and others babiolles 
  call main_train_mld

  ! ----> if we want to compute the condition number of Phi^T Phi ....
  if (scalapack_driver) then 
    call log_critical('Estimation of condition number of PhixPhi^T ....: not implemented for scalapack')
    call mld_mpi_abort('stop for condition number in '//NAMECURRENT)
  else 
    allocate(symmat(size(Amat,1), size(Amat,1)))
    call dsyrk('U', 'N', size(Amat,1), size(Amat,2), 1.d0, Amat, size(Amat,1), 0.d0, symmat, size(Amat,1))
    cond = compute_condition_number_mat_rss(symmat)
    deallocate(symmat)
    write(ctmp, '(es15.7)')  cond
    call log_info('Estimation of condition number of PhixPhi^T ....: '//trim(ctmp))
  end if 
  ! ----> END condition number of Phi^T Phi ....

  lMAP0=.false. 
  if (lMAP0) then 
  allocate(xini, source=w_params(:,1))
  allocate(xend, source=w_params(:,1))
  if (allocated(xx)) deallocate(xx) ; allocate(xx(size(xini)))
  ! Allocate and initialize quadratic function object
  allocate(FuncMargLikeMAP0 :: func)
  func%fdim = size(w_params,1)
  xini(:) = w_params(:,1)
  xx  = xini 
  call func%init(xx)
  !allocate(ConjugateGradientType :: minimizer)
  !allocate(SteepestDescentType :: minimizer )
  allocate(AdamDescentType :: minimizer )
  call log_info('min: Adam first minimization ...  '//NAMECURRENT) 
  call minimizer%init(max_iter=100, tol=1.d-4, xx=xx) 
  call minimizer%move(xx, func)  
  xend = minimizer%xx_end
  deallocate(minimizer)
  write(*,*) '|<---------------------------------->|'
  xx = xend
  call log_info('min: LBFGS second minimization ...  '//NAMECURRENT) 
  allocate(LbfgsExtType :: minimizer )
  lbfgs_m_hess = 40 
  lbfgs_xtol = 1.d-12
  call minimizer%init(max_iter=10000, tol=1.d-4, xx=xx) 
  call minimizer%move(xx, func)  
  xend = minimizer%xx_end
  deallocate(minimizer)
  w_params(:,1) = xend(:)
  end if !lMAP0 

  lAlphaBeta=.true.
  if (lAlphaBeta) then
   call log_info('start AlphaBeta model ... ')
   allocate(FuncMargLikeAlphaBeta :: funcAB)
   ! xini should be defined outside
   allocate(xini(2*size(w_params,1)+3))
   if (allocated(xx)) deallocate(xx) ; allocate(xx(size(xini)))
   xini=0.0_dp
   xx = xini
   call funcAB%init(xx)
   call log_info('end AlphaBeta model ... ')

   stop 'to be continued'
  end if !lAlphaBeta 

  _MLD_END_
end subroutine main_train_marginal_likelihood

