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

!-------------------------------------------------------------------------------
! module_serial_als  –  Serial ALS-Ridge solver with LAPACK
!
! Implements the block-preconditioning + ALS-Ridge fitting algorithm from
! MiladyNoteTechnique5.pdf, sections 1.10.1 (fixed alpha) and 1.10.2
! (learning alpha).
!
! This module is SELF-CONTAINED: it does not use any MiLaDy global state.
! All data enters through subroutine arguments.
!
! Public interface:
!   serial_als_ridge_solve  –  main entry point
!   serial_nnls             –  non-negative least squares (Lawson–Hanson)
!-------------------------------------------------------------------------------
module module_serial_als
  use module_kind_variables, only: kind_double
  implicit none
  private

  public :: serial_als_ridge_solve
  public :: serial_nnls

contains

  !-----------------------------------------------------------------------------
  ! serial_als_ridge_solve
  !
  ! Solve  y ≈ sum_{b=1}^{nu_max}  alpha_b  A^(b)  x^(b)
  ! by alternating ridge regression (w-step) and alpha regression (alpha-step).
  !
  ! If alpha_method == 0 (fixed), alpha is computed once from Frobenius/SVD
  ! preconditioning and kept constant (section 1.10.1).
  ! If alpha_method == 1 (learning), alpha is updated at each ALS iteration
  ! (section 1.10.2).
  !
  ! Arguments
  ! ---------
  ! AA(D,M)             : full design matrix  (D = descriptor dim, M = data)
  ! D, M                : dimensions of AA
  ! yy(M)               : target vector
  ! nu_max              : number of blocks
  ! block_part(2*nu_max) : block partition, block b occupies columns
  !                        block_part(2*b-1) : block_part(2*b)  of AA
  ! ridge_k             : k > 0, global ridge parameter
  ! rho                 : regularisation for the alpha sub-problem
  ! nsteps              : S, maximum number of ALS iterations
  ! tol                 : convergence tolerance
  ! alpha_method        : 0 = fixed alpha, 1 = learning alpha
  ! precond_type        : 0 = Frobenius norm, 1 = SVD
  ! use_nnls_alpha      : .true. => enforce alpha >= 0 via NNLS in alpha-step
  ! nnls_mode           : 0 = Lawson-Hanson (cold start), 1 = BK warm-start NNLS
  !
  ! Output
  ! ------
  ! w_out(D)            : fitted parameter vector
  ! alpha_out(nu_max)   : block renormalisation coefficients
  ! info                : 0 on success, >0 on failure
  !-----------------------------------------------------------------------------
  subroutine serial_als_ridge_solve(AA, D, M, yy, nu_max, block_part,  &
                                    ridge_k, rho, nsteps, tol,          &
                                    alpha_method, precond_type,         &
                                    use_nnls_alpha, nnls_mode,          &
                                    w_out, alpha_out, info)
    use mld_logger, only: log_info
    use mld_string, only: vtoa
    use module_input_als_fit, only: als_precond_flat, als_precond_svd, &
                                   als_nnls_bk_warm
    use module_nnls_bk, only: nnls_solver
    implicit none
    ! ---- arguments ----
    integer,  intent(in)    :: D, M, nu_max, nsteps
    integer,  intent(in)    :: alpha_method, precond_type
    logical,  intent(in)    :: use_nnls_alpha
    integer,  intent(in)    :: nnls_mode
    real(kind_double), intent(in)    :: ridge_k, rho, tol
    real(kind_double), intent(in)    :: AA(D, M)
    real(kind_double), intent(in)    :: yy(M)
    integer,  intent(in)    :: block_part(2*nu_max)
    real(kind_double), intent(out)   :: w_out(D)
    real(kind_double), intent(out)   :: alpha_out(nu_max)
    integer,  intent(out)   :: info

    ! ---- local variables ----
    real(kind_double), allocatable :: Atilde(:,:)     ! scaled design matrix (D,M)
    real(kind_double), allocatable :: lambda(:)       ! per-column ridge vector (D)
    real(kind_double), allocatable :: frob_nu(:)      ! squared Frobenius norms per block
    real(kind_double), allocatable :: alpha(:)        ! current alpha
    real(kind_double), allocatable :: alpha_prev(:)   ! previous alpha
    real(kind_double), allocatable :: w(:)            ! current w (D)
    real(kind_double), allocatable :: w_prev(:)       ! previous w (D)
    real(kind_double), allocatable :: BtB(:,:)        ! B^T B for alpha-step (nu_max x nu_max)
    real(kind_double), allocatable :: Bty(:)          ! B^T y for alpha-step (nu_max)
    real(kind_double), allocatable :: bvec(:,:)       ! b^(nu) vectors (M, nu_max)
    real(kind_double), allocatable :: alpha_rhs(:,:)  ! for alpha solve (nu_max x 1)
    real(kind_double), allocatable :: alpha_phi(:,:)  ! BtB + rho*I (nu_max x nu_max)
    integer :: b, ib, ie, Db, s, i
    real(kind_double) :: norm_w, norm_alpha, diff_w, diff_alpha
    real(kind_double) :: svd_max_sing
    real(kind_double), external :: dnrm2
    integer :: lapack_info
    ! BK warm-start NNLS solver (used when nnls_mode == als_nnls_bk_warm)
    type(nnls_solver) :: bk_solver
    logical, allocatable :: bk_passive(:)
    real(kind_double) :: bk_rnorm

    info = 0

    call log_info("ALS-Ridge serial: D=" // vtoa(D) // " M=" // vtoa(M) &
                  // " nu_max=" // vtoa(nu_max) // " nsteps=" // vtoa(nsteps) &
                  // " alpha_method=" // vtoa(alpha_method) // " precond_type=" // vtoa(precond_type))

    ! ---- Step A.1: Compute squared Frobenius norm per block ----
    allocate(frob_nu(nu_max))
    do b = 1, nu_max
      ib = block_part(2*b - 1)
      ie = block_part(2*b)
      frob_nu(b) = frobenius_sq(AA, D, M, ib, ie)
    end do

    ! ---- Step A.1: Build per-column ridge vector Lambda ----
    allocate(lambda(D))
    lambda(:) = 0.0d0
    if (precond_type == als_precond_flat) then
      ! Flat mode: uniform lambda = ridge_k (caller passes lambda_krr)
      lambda(:) = ridge_k
      call log_info("ALS-Ridge serial: flat precond, lambda = " // vtoa(ridge_k))
    else
      ! Frobenius-scaled: lambda_i = ridge_k * frob_nu(b) / Db
      do b = 1, nu_max
        ib = block_part(2*b - 1)
        ie = block_part(2*b)
        Db = ie - ib + 1
        if (Db > 0) then
          do i = ib, ie
            lambda(i) = ridge_k * frob_nu(b) / dble(Db)
          end do
        end if
      end do
      ! ---- Print per-block Lambda values ----
      do b = 1, nu_max
        ib = block_part(2*b - 1)
        ie = block_part(2*b)
        Db = ie - ib + 1
        if (Db > 0) then
          call log_info("ALS-Ridge serial: block " // vtoa(b) // &
                        " lambda = " // vtoa(lambda(ib)) //      &
                        " (ridge_k=" // vtoa(ridge_k) //         &
                        " frob=" // vtoa(frob_nu(b)) //          &
                        " Db=" // vtoa(Db) // ")")
        end if
      end do
    end if

    ! ---- Step A.2: Initial alpha ----
    allocate(alpha(nu_max), alpha_prev(nu_max))
    if (precond_type == als_precond_flat) then
      ! Flat mode: no preconditioning, alpha = 1
      alpha(:) = 1.0d0
    else if (precond_type == als_precond_svd) then
      ! SVD-based: alpha_nu = 1 / max_singular_value(A^(nu))
      do b = 1, nu_max
        ib = block_part(2*b - 1)
        ie = block_part(2*b)
        Db = ie - ib + 1
        call max_singular_value(AA(:, :), D, M, ib, ie, svd_max_sing)
        if (svd_max_sing > 0.0d0) then
          alpha(b) = 1.0d0 / svd_max_sing
        else
          alpha(b) = 1.0d0
        end if
      end do
    else
      ! Frobenius-based: alpha_nu = 1 / sqrt(f_nu / D_nu)
      do b = 1, nu_max
        ib = block_part(2*b - 1)
        ie = block_part(2*b)
        Db = ie - ib + 1
        if (Db > 0 .and. frob_nu(b) > 0.0d0) then
          alpha(b) = 1.0d0 / dsqrt(frob_nu(b) / dble(Db))
        else
          alpha(b) = 1.0d0
        end if
      end do
    end if

    ! ---- Build scaled matrix Atilde = [alpha_1 A^(1), ..., alpha_nu A^(nu)] ----
    allocate(Atilde(D, M))
    call build_Atilde(AA, D, M, nu_max, block_part, alpha, Atilde)

    ! ---- Allocate work arrays for the w-step ----
    allocate(w(D), w_prev(D))
    w(:) = 0.0d0

    ! ---- Allocate work arrays for the alpha-step ----
    allocate(bvec(M, nu_max))
    allocate(BtB(nu_max, nu_max), Bty(nu_max))

    ! ========================================================================
    ! Main ALS loop (section 1.10.2, steps B.1 – B.3)
    ! For alpha_method == 0 (fixed alpha), the loop runs a single w-step.
    ! ========================================================================
    do s = 0, nsteps - 1

      w_prev(:) = w(:)
      alpha_prev(:) = alpha(:)

      ! ---- B.1  w-step: ridge regression on the scaled matrix ----
      !   w <- (Atilde^T Atilde + Lambda)^{-1} Atilde^T y
      call ridge_solve_serial(Atilde, D, M, yy, lambda, w, lapack_info)
      if (lapack_info /= 0) then
        info = 1
        return
      end if

      ! If fixed alpha, we are done after one pass
      if (alpha_method == 0) exit

      ! ---- B.2  alpha-step: small regression ----
      !   b^(nu) = A^(nu) w^(nu)
      !   B = [b^(1) ... b^(nu_max)]  (M x nu_max)
      !   alpha <- (B^T B + rho I)^{-1} B^T y
      do b = 1, nu_max
        ib = block_part(2*b - 1)
        ie = block_part(2*b)
        ! b^(nu) = A^(nu) * w^(nu)
        ! A^(nu) is AA(ib:ie, 1:M), so b = A(ib:ie,:)^T * w(ib:ie)
        ! but AA is (D,M) so A^(nu) in the paper is AA(ib:ie,:)
        ! b^(nu) in R^M:   b_j = sum_{d=ib}^{ie}  AA(d,j) * w(d)
        call dgemv_wrapper(AA, D, M, ib, ie, w, bvec(:, b))
      end do

      ! B^T B
      ! B is (M, nu_max)
      call dgemm('T', 'N', nu_max, nu_max, M, 1.0d0, bvec, M, bvec, M, 0.0d0, BtB, nu_max)

      ! B^T y
      call dgemv('T', M, nu_max, 1.0d0, bvec, M, yy, 1, 0.0d0, Bty, 1)

      ! add rho*I to B^T B
      do i = 1, nu_max
        BtB(i, i) = BtB(i, i) + rho
      end do

      ! Solve for alpha
      if (use_nnls_alpha) then
        if (nnls_mode == als_nnls_bk_warm) then
          ! BK warm-start NNLS on pre-formed Gram system:
          !   G = BtB (with rho*I already added), c = Bty
          ! Solves  min_{alpha >= 0}  G*alpha = c  via active-set.
          ! Uses solve_gram_warm to avoid squaring the normal equations.
          bk_solver%tol = 1.0d-12
          bk_solver%jitter0 = 1.0d-14
          bk_solver%jitter_max = 1.0d-4
          if (.not. allocated(bk_passive)) allocate(bk_passive(nu_max))
          call bk_solver%solve_gram_warm(BtB, Bty, alpha, lapack_info, passive=bk_passive)
        else
          ! Lawson-Hanson NNLS (cold start)
          allocate(alpha_phi(nu_max, nu_max))
          alpha_phi(:,:) = BtB(:,:)
          call serial_nnls(alpha_phi, nu_max, nu_max, Bty, alpha, lapack_info)
          deallocate(alpha_phi)
        end if
      else
        ! standard solve: (BtB) alpha = Bty
        allocate(alpha_phi(nu_max, nu_max), alpha_rhs(nu_max, 1))
        alpha_phi(:,:) = BtB(:,:)
        alpha_rhs(:, 1) = Bty(:)
        call solve_small_system(alpha_phi, nu_max, alpha_rhs, alpha, lapack_info)
        deallocate(alpha_phi, alpha_rhs)
      end if

      if (lapack_info /= 0) then
        call log_info("ALS-Ridge serial: alpha solve failed at step " // vtoa(s) // " info=" // vtoa(lapack_info))
        info = 2
        return
      end if

      ! ---- Log alpha at each step ----
      call log_info("ALS-Ridge serial: step " // vtoa(s) // " alpha = [" // vtoa(alpha(1:nu_max)) // "]")

      ! ---- Rebuild Atilde with updated alpha ----
      call build_Atilde(AA, D, M, nu_max, block_part, alpha, Atilde)

      ! ---- B.3  Stopping criterion (skipped when tol < 0) ----
      if (tol >= 0.0d0) then
        norm_w = dnrm2(D, w, 1)
        diff_w = 0.0d0
        do i = 1, D
          diff_w = diff_w + (w(i) - w_prev(i))**2
        end do
        diff_w = dsqrt(diff_w)

        norm_alpha = dnrm2(nu_max, alpha, 1)
        diff_alpha = 0.0d0
        do i = 1, nu_max
          diff_alpha = diff_alpha + (alpha(i) - alpha_prev(i))**2
        end do
        diff_alpha = dsqrt(diff_alpha)

        if (norm_w > 0.0d0 .and. norm_alpha > 0.0d0) then
          if ((diff_w / norm_w < tol) .and. (diff_alpha / norm_alpha < tol)) exit
        end if
      end if

    end do ! s = 0, nsteps-1

    ! ---- Output ----
    ! The ALS routine delivers both w^(nu) and alpha_nu.
    ! For extrapolation the effective parameters are the product
    !   w_tilde^(nu) = alpha_nu * w^(nu)      (sec. 1.10.1 of the PDF)
    ! because at prediction time the raw descriptors A^(nu) are used
    ! (not the scaled Atilde^(nu) = alpha_nu * A^(nu)).
    w_out(:) = 0.0d0
    do b = 1, nu_max
      ib = block_part(2*b - 1)
      ie = block_part(2*b)
      do i = ib, ie
        w_out(i) = alpha(b) * w(i)
      end do
    end do
    alpha_out(:) = alpha(:)

    ! ---- Cleanup ----
    if (allocated(bk_passive)) deallocate(bk_passive)
    deallocate(frob_nu, lambda, alpha, alpha_prev, Atilde)
    deallocate(w, w_prev, bvec, BtB, Bty)

  end subroutine serial_als_ridge_solve


  !-----------------------------------------------------------------------------
  ! ridge_solve_serial
  !
  ! Solve  (A^T A + diag(lambda))  w  =  A^T y
  ! using LAPACK dsysv (symmetric indefinite factorisation).
  !
  ! A is (D, M),  y is (M),  lambda is (D),  w is (D).
  !-----------------------------------------------------------------------------
  subroutine ridge_solve_serial(A, D, M, y, lambda, w, info)
    implicit none
    integer,  intent(in)   :: D, M
    real(kind_double), intent(in)   :: A(D, M), y(M), lambda(D)
    real(kind_double), intent(out)  :: w(D)
    integer,  intent(out)  :: info

    real(kind_double), allocatable :: phi(:,:), rhs(:,:)
    integer :: i
    integer :: lwork
    real(kind_double), allocatable :: work(:)
    integer, allocatable :: ipiv(:)

    allocate(phi(D, D), rhs(D, 1), ipiv(D))

    ! phi = A * A^T  (since A is D x M, this gives D x D)
    call dgemm('N', 'T', D, D, M, 1.0d0, A, D, A, D, 0.0d0, phi, D)

    ! add diagonal regularisation
    do i = 1, D
      phi(i, i) = phi(i, i) + lambda(i)
    end do

    ! rhs = A * y
    call dgemv('N', D, M, 1.0d0, A, D, y, 1, 0.0d0, rhs(:,1), 1)

    ! workspace query
    allocate(work(1))
    lwork = -1
    call dsysv('U', D, 1, phi, D, ipiv, rhs, D, work, lwork, info)
    lwork = int(work(1)) + 2
    deallocate(work)
    allocate(work(lwork))

    ! solve
    call dsysv('U', D, 1, phi, D, ipiv, rhs, D, work, lwork, info)

    w(:) = rhs(:, 1)

    deallocate(phi, rhs, ipiv, work)
  end subroutine ridge_solve_serial


  !-----------------------------------------------------------------------------
  ! serial_nnls  –  Non-Negative Least Squares (Lawson–Hanson algorithm)
  !
  ! Solve   min || A x - b ||_2   subject to  x >= 0
  !
  ! A is (m, n),  b is (m),  x is (n).
  ! A and b are MODIFIED on output.
  !
  ! Reference: Lawson & Hanson, "Solving Least Squares Problems", Ch. 23.
  !-----------------------------------------------------------------------------
  subroutine serial_nnls(A, m, n, b, x, info)
    implicit none
    integer,  intent(in)    :: m, n
    real(kind_double), intent(inout) :: A(m, n)
    real(kind_double), intent(inout) :: b(m)
    real(kind_double), intent(out)   :: x(n)
    integer,  intent(out)   :: info

    integer, parameter :: MAX_ITER_NNLS = 3000

    real(kind_double), allocatable :: w_nnls(:), zz(:), Acopy(:,:)
    logical, allocatable :: passive(:)
    integer :: i, j, t_idx, nP, iter_main, iter_inner
    real(kind_double) :: w_max, alpha_nnls, ratio
    integer :: lapack_info, lwork_nnls
    real(kind_double), allocatable :: work_qr(:)
    integer :: nrhs_nnls

    info = 0
    allocate(w_nnls(n), zz(m), passive(n))
    x(:) = 0.0d0
    passive(:) = .false.

    ! w = A^T (b - A x)   (gradient of 0.5 ||b - Ax||^2 w.r.t. x, negated)
    ! initially x=0 so w = A^T b
    call dgemv('T', m, n, 1.0d0, A, m, b, 1, 0.0d0, w_nnls, 1)

    do iter_main = 1, MAX_ITER_NNLS

      ! find max w among non-passive variables
      w_max = -huge(1.0d0)
      t_idx = 0
      do j = 1, n
        if (.not. passive(j)) then
          if (w_nnls(j) > w_max) then
            w_max = w_nnls(j)
            t_idx = j
          end if
        end if
      end do

      ! if all w <= 0 for non-passive set, or no non-passive left, we are done
      if (t_idx == 0 .or. w_max <= 0.0d0) exit

      ! move t_idx into passive set
      passive(t_idx) = .true.

      ! inner loop: solve LS on passive set, fix negative components
      do iter_inner = 1, MAX_ITER_NNLS

        ! count passive variables
        nP = 0
        do j = 1, n
          if (passive(j)) nP = nP + 1
        end do

        if (nP == 0) exit

        ! solve unconstrained LS on passive columns: min || A_P z_P - b ||
        ! We use a temporary copy and LAPACK dgels
        allocate(Acopy(m, nP))
        nP = 0
        do j = 1, n
          if (passive(j)) then
            nP = nP + 1
            Acopy(:, nP) = A(:, j)
          end if
        end do

        zz(1:m) = b(1:m)
        nrhs_nnls = 1

        ! workspace query for dgels
        allocate(work_qr(1))
        lwork_nnls = -1
        call dgels('N', m, nP, nrhs_nnls, Acopy, m, zz, m, work_qr, lwork_nnls, lapack_info)
        lwork_nnls = int(work_qr(1)) + 2
        deallocate(work_qr)
        allocate(work_qr(lwork_nnls))
        call dgels('N', m, nP, nrhs_nnls, Acopy, m, zz, m, work_qr, lwork_nnls, lapack_info)
        deallocate(work_qr, Acopy)

        if (lapack_info /= 0) then
          info = lapack_info
          deallocate(w_nnls, zz, passive)
          return
        end if

        ! check if all z_P >= 0
        ! put z into x for passive, 0 for non-passive
        nP = 0
        alpha_nnls = huge(1.0d0)
        do j = 1, n
          if (passive(j)) then
            nP = nP + 1
            if (zz(nP) <= 0.0d0) then
              ! find interpolation ratio
              if (x(j) - zz(nP) > 0.0d0) then
                ratio = x(j) / (x(j) - zz(nP))
                if (ratio < alpha_nnls) alpha_nnls = ratio
              end if
            end if
          end if
        end do

        if (alpha_nnls >= 1.0d0) then
          ! all z_P >= 0, accept solution
          nP = 0
          do j = 1, n
            if (passive(j)) then
              nP = nP + 1
              x(j) = zz(nP)
            end if
          end do
          exit  ! exit inner loop
        end if

        ! interpolate: x <- x + alpha * (z - x)
        nP = 0
        do j = 1, n
          if (passive(j)) then
            nP = nP + 1
            x(j) = x(j) + alpha_nnls * (zz(nP) - x(j))
            ! remove near-zero variables from passive set
            if (dabs(x(j)) < 1.0d-14) then
              passive(j) = .false.
              x(j) = 0.0d0
            end if
          end if
        end do

      end do ! iter_inner

      ! recompute w = A^T (b - A x)
      ! residual = b - A x
      zz(1:m) = b(1:m)
      call dgemv('N', m, n, -1.0d0, A, m, x, 1, 1.0d0, zz, 1)
      call dgemv('T', m, n, 1.0d0, A, m, zz, 1, 0.0d0, w_nnls, 1)

      ! zero out w for passive variables (not needed for convergence check)
      do j = 1, n
        if (passive(j)) w_nnls(j) = 0.0d0
      end do

    end do ! iter_main

    deallocate(w_nnls, zz, passive)
  end subroutine serial_nnls


  !-----------------------------------------------------------------------------
  ! build_Atilde  –  Atilde = [alpha_1 A^(1), ..., alpha_nu A^(nu)]
  !-----------------------------------------------------------------------------
  subroutine build_Atilde(AA, D, M, nu_max, block_part, alpha, Atilde)
    implicit none
    integer,  intent(in)  :: D, M, nu_max
    real(kind_double), intent(in)  :: AA(D, M)
    integer,  intent(in)  :: block_part(2*nu_max)
    real(kind_double), intent(in)  :: alpha(nu_max)
    real(kind_double), intent(out) :: Atilde(D, M)
    integer :: b, ib, ie, i

    Atilde(:,:) = 0.0d0
    do b = 1, nu_max
      ib = block_part(2*b - 1)
      ie = block_part(2*b)
      do i = ib, ie
        Atilde(i, :) = alpha(b) * AA(i, :)
      end do
    end do
  end subroutine build_Atilde


  !-----------------------------------------------------------------------------
  ! frobenius_sq  –  squared Frobenius norm of a sub-block of A
  ! block = A(ib:ie, 1:M)
  !-----------------------------------------------------------------------------
  function frobenius_sq(A, D, M, ib, ie) result(fsq)
    implicit none
    integer,  intent(in) :: D, M, ib, ie
    real(kind_double), intent(in) :: A(D, M)
    real(kind_double) :: fsq
    integer :: i, j

    fsq = 0.0d0
    do j = 1, M
      do i = ib, ie
        fsq = fsq + A(i, j)**2
      end do
    end do
  end function frobenius_sq


  !-----------------------------------------------------------------------------
  ! max_singular_value  –  compute largest singular value of sub-block A(ib:ie,:)
  ! Uses LAPACK dgesvd with 'N','N' (no vectors, only singular values).
  !-----------------------------------------------------------------------------
  subroutine max_singular_value(A, D, M, ib, ie, smax)
    implicit none
    integer,  intent(in) :: D, M, ib, ie
    real(kind_double), intent(in) :: A(D, M)
    real(kind_double), intent(out) :: smax

    integer :: Db, mn, lwork, lapack_info
    real(kind_double), allocatable :: Acopy(:,:), S(:), work(:)

    Db = ie - ib + 1
    mn = min(Db, M)
    allocate(Acopy(Db, M), S(mn))

    Acopy(:,:) = A(ib:ie, 1:M)

    ! workspace query
    allocate(work(1))
    lwork = -1
    call dgesvd('N', 'N', Db, M, Acopy, Db, S, Acopy, Db, Acopy, Db, work, lwork, lapack_info)
    lwork = int(work(1)) + 2
    deallocate(work)
    allocate(work(lwork))

    call dgesvd('N', 'N', Db, M, Acopy, Db, S, Acopy, Db, Acopy, Db, work, lwork, lapack_info)

    if (lapack_info == 0 .and. mn > 0) then
      smax = S(1)
    else
      smax = 1.0d0
    end if

    deallocate(Acopy, S, work)
  end subroutine max_singular_value


  !-----------------------------------------------------------------------------
  ! dgemv_wrapper  –  b = A(ib:ie, :)^T * w(ib:ie)
  ! result is b(1:M)
  !-----------------------------------------------------------------------------
  subroutine dgemv_wrapper(A, D, M, ib, ie, w, b_out)
    implicit none
    integer,  intent(in) :: D, M, ib, ie
    real(kind_double), intent(in) :: A(D, M), w(D)
    real(kind_double), intent(out) :: b_out(M)
    integer :: Db

    Db = ie - ib + 1
    ! b = A(ib:ie,:)^T * w(ib:ie)
    ! A(ib:ie,:) is Db x M, so its transpose is M x Db
    call dgemv('T', Db, M, 1.0d0, A(ib, 1), D, w(ib), 1, 0.0d0, b_out, 1)
  end subroutine dgemv_wrapper


  !-----------------------------------------------------------------------------
  ! solve_small_system  –  solve phi * x = rhs for x
  ! phi is (n, n), rhs is (n, 1), uses dsysv.
  !-----------------------------------------------------------------------------
  subroutine solve_small_system(phi, n, rhs, x, info)
    implicit none
    integer,  intent(in)    :: n
    real(kind_double), intent(inout) :: phi(n, n), rhs(n, 1)
    real(kind_double), intent(out)   :: x(n)
    integer,  intent(out)   :: info

    integer :: lwork
    integer, allocatable :: ipiv(:)
    real(kind_double), allocatable :: work(:)

    allocate(ipiv(n))

    ! workspace query
    allocate(work(1))
    lwork = -1
    call dsysv('U', n, 1, phi, n, ipiv, rhs, n, work, lwork, info)
    lwork = int(work(1)) + 2
    deallocate(work)
    allocate(work(lwork))

    call dsysv('U', n, 1, phi, n, ipiv, rhs, n, work, lwork, info)

    x(:) = rhs(:, 1)

    deallocate(ipiv, work)
  end subroutine solve_small_system

end module module_serial_als
