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
! module_para_als  –  Parallel ALS-Ridge solver with ScaLAPACK
!
! Implements the block-preconditioning + ALS-Ridge fitting algorithm from
! MiladyNoteTechnique5.pdf, sections 1.10.1 (fixed alpha) and 1.10.2
! (learning alpha).
!
! This module is SELF-CONTAINED: it does not use any MiLaDy global state.
! All data enters through subroutine arguments.
! The NNLS alpha sub-problem is solved serially (nu_max is tiny),
! reusing the serial_nnls routine from module_serial_als.
!
! ScaLAPACK patterns follow the conventions of mld_scalapack.F90
! (allocation_scalapack_matrix, pdgemm, pdgels, pdposv, pdelset, pdelget,
!  INDXL2G, numroc, descinit, blockset_ml, blacs_barrier).
!
! Public interface:
!   para_als_ridge_solve  –  main entry point (ScaLAPACK)
!-------------------------------------------------------------------------------
module module_para_als
  use module_kind_variables, only: kind_double
  implicit none
  private

  public :: para_als_ridge_solve

contains

  !-----------------------------------------------------------------------------
  ! para_als_ridge_solve
  !
  ! Parallel version of the ALS-Ridge solver.
  !
  ! The design matrix sca_AA is already distributed on the BLACS grid.
  ! All heavy linear algebra (ridge w-step) is done in ScaLAPACK.
  ! The alpha-step involves a tiny nu_max × nu_max system, solved serially
  ! after an MPI reduction of B^T B and B^T y.
  !
  ! Arguments
  ! ---------
  ! sca_AA(l_dimr, l_dimc) : local tile of the distributed design matrix (D × M)
  ! desc_AA(9)             : ScaLAPACK descriptor for sca_AA
  ! D, M                   : global dimensions of the design matrix
  ! nbr_AA, nbc_AA         : block sizes (row, col) used for sca_AA
  ! sca_yy(l_dimr_y, l_dimc_y) : local tile of the distributed target vector (M × 1)
  ! desc_yy(9)             : ScaLAPACK descriptor for sca_yy
  ! nu_max                 : number of blocks
  ! block_part(2*nu_max)   : block partition (column ranges in global indexing)
  ! ridge_k                : k > 0, global ridge parameter
  ! rho                    : regularisation for the alpha sub-problem
  ! nsteps                 : S, maximum number of ALS iterations
  ! tol                    : convergence tolerance
  ! alpha_method           : 0 = fixed alpha, 1 = learning alpha
  ! precond_type           : 0 = Frobenius norm, 1 = SVD
  ! use_nnls_alpha         : .true. => enforce alpha >= 0 via NNLS in alpha-step
  ! context                : BLACS context handle
  ! nprow, npcol           : BLACS process grid dimensions
  ! myrow, mycol           : local process coordinates in BLACS grid
  !
  ! Output
  ! ------
  ! w_out(D)               : fitted parameter vector (replicated on all procs),
  !                          renormalized: w_out^(nu) = alpha_nu * w^(nu)
  !                          (eq. in sec. 1.10.1). These are the effective
  !                          coefficients for extrapolation since the raw
  !                          (un-scaled) descriptors A^(nu) are used at
  !                          prediction time.
  ! alpha_out(nu_max)      : block renormalisation coefficients (replicated)
  ! info                   : 0 on success, >0 on failure
  !-----------------------------------------------------------------------------
  subroutine para_als_ridge_solve(sca_AA, desc_AA, D, M, nbr_AA, nbc_AA,  &
                                  sca_yy, desc_yy,                         &
                                  nu_max, block_part,                      &
                                  ridge_k, rho, nsteps, tol,               &
                                  alpha_method, precond_type,              &
                                  use_nnls_alpha, nnls_mode,               &
                                  context, nprow, npcol, myrow, mycol,     &
                                  w_out, alpha_out, info)
    use module_serial_als, only: serial_nnls
    use module_scalapack_tools, only: blockset_ml
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
    real(kind_double), intent(in)    :: sca_AA(:,:)
    integer,  intent(in)    :: desc_AA(9)
    integer,  intent(in)    :: nbr_AA, nbc_AA
    real(kind_double), intent(in)    :: sca_yy(:,:)
    integer,  intent(in)    :: desc_yy(9)
    integer,  intent(in)    :: block_part(2*nu_max)
    integer,  intent(in)    :: context, nprow, npcol, myrow, mycol
    real(kind_double), intent(out)   :: w_out(D)
    real(kind_double), intent(out)   :: alpha_out(nu_max)
    integer,  intent(out)   :: info

    ! ---- local variables ----
    ! ScaLAPACK distributed matrices for the ridge w-step
    real(kind_double), allocatable :: sca_Atilde(:,:)  ! distributed scaled design matrix (D × M)
    integer :: desc_Atilde(9)
    integer :: l_dimr_Atilde, l_dimc_Atilde
    integer :: nbr_Atilde, nbc_Atilde

    real(kind_double), allocatable :: sca_phi(:,:)     ! distributed AtA + Lambda (D × D)
    integer :: desc_phi(9)
    integer :: l_dimr_phi, l_dimc_phi
    integer :: nbr_phi, nbc_phi

    real(kind_double), allocatable :: sca_rhs(:,:)     ! distributed rhs = At*y (D × 1)
    integer :: desc_rhs(9)
    integer :: l_dimr_rhs, l_dimc_rhs
    integer :: nbr_rhs, nbc_rhs

    ! Block norms and preconditioning
    real(kind_double), allocatable :: frob_nu(:)       ! squared Frobenius norm per block
    real(kind_double), allocatable :: lambda(:)        ! per-row ridge vector (D)
    real(kind_double), allocatable :: alpha(:)         ! current alpha (nu_max)
    real(kind_double), allocatable :: alpha_prev(:)    ! previous alpha (nu_max)
    real(kind_double), allocatable :: w(:)             ! current w (D) — replicated
    real(kind_double), allocatable :: w_prev(:)        ! previous w (D)

    ! Alpha sub-problem (small, solved serially after reduction)
    real(kind_double), allocatable :: BtB(:,:)         ! B^T B (nu_max × nu_max) replicated
    real(kind_double), allocatable :: Bty(:)           ! B^T y (nu_max) replicated
    real(kind_double), allocatable :: BtB_local(:,:)   ! local partial sum
    real(kind_double), allocatable :: Bty_local(:)     ! local partial sum
    real(kind_double), allocatable :: alpha_phi(:,:)   ! BtB + rho*I for solve
    real(kind_double), allocatable :: alpha_rhs(:,:)   ! Bty for solve

    ! Temporary vectors for b^(nu) computation
    real(kind_double), allocatable :: bvec_local(:,:)  ! local partial b vectors (M, nu_max)
    real(kind_double), allocatable :: bvec(:,:)        ! global b vectors (M, nu_max)

    ! SVD workspace for preconditioning
    ! (SVD computation is in compute_svd_alpha helper)

    ! Loop and BLACS variables
    integer :: b, ib, ie, Db, s, i
    real(kind_double) :: norm_w, norm_alpha, diff_w, diff_alpha
    real(kind_double), external :: dnrm2
    integer :: sca_info
    integer :: numroc  ! ScaLAPACK external function
    ! BK warm-start NNLS solver (used when nnls_mode == als_nnls_bk_warm)
    type(nnls_solver) :: bk_solver
    logical, allocatable :: bk_passive(:)
    real(kind_double) :: bk_rnorm

    info = 0

    call log_info("ALS-Ridge ScaLAPACK: D=" // vtoa(D) // " M=" // vtoa(M) &
                  // " nu_max=" // vtoa(nu_max) // " nsteps=" // vtoa(nsteps) &
                  // " alpha_method=" // vtoa(alpha_method) // " precond_type=" // vtoa(precond_type))

    ! ========================================================================
    ! Step A.1: Compute squared Frobenius norm per block (distributed)
    ! ========================================================================
    allocate(frob_nu(nu_max))
    call compute_distributed_frobenius(sca_AA, desc_AA, D, M, nbr_AA, nbc_AA, &
                                       nu_max, block_part, context,           &
                                       myrow, mycol, nprow, npcol, frob_nu)

    ! ========================================================================
    ! Step A.1 continued: Build per-row ridge vector Lambda
    ! ========================================================================
    allocate(lambda(D))
    lambda(:) = 0.0d0
    if (precond_type == als_precond_flat) then
      ! Flat mode: uniform lambda = ridge_k (caller passes lambda_krr)
      lambda(:) = ridge_k
      call log_info("ALS-Ridge ScaLAPACK: flat precond, lambda = " // vtoa(ridge_k))
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
          call log_info("ALS-Ridge ScaLAPACK: block " // vtoa(b) // &
                        " lambda = " // vtoa(lambda(ib)) //         &
                        " (ridge_k=" // vtoa(ridge_k) //            &
                        " frob=" // vtoa(frob_nu(b)) //             &
                        " Db=" // vtoa(Db) // ")")
        end if
      end do
    end if

    ! ========================================================================
    ! Step A.2: Initial alpha from preconditioning
    ! ========================================================================
    allocate(alpha(nu_max), alpha_prev(nu_max))
    if (precond_type == als_precond_flat) then
      ! Flat mode: no preconditioning, alpha = 1
      alpha(:) = 1.0d0
    else if (precond_type == als_precond_svd) then
      ! SVD-based: alpha_nu = 1 / max_singular_value(A^(nu))
      ! For the parallel case, we use pdlange('F',...) for Frobenius
      ! and fall back to serial SVD on the gathered sub-block.
      ! Since SVD preconditioning is expensive, we gather each sub-block
      ! (Db x M, typically small Db) to rank 0 and compute SVD there,
      ! then broadcast.
      call compute_svd_alpha(sca_AA, desc_AA, D, M, nbr_AA, nbc_AA,  &
                             nu_max, block_part, context,             &
                             myrow, mycol, nprow, npcol, alpha)
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

    ! ========================================================================
    ! Build distributed scaled matrix Atilde = [alpha_1 A^(1), ..., alpha_nu A^(nu)]
    ! ========================================================================
    ! Allocate sca_Atilde with same layout as sca_AA
    call blockset_ml(nbr_Atilde, nbr_AA, D, nprow, npcol)
    call blockset_ml(nbc_Atilde, nbc_AA, M, nprow, npcol)
    l_dimr_Atilde = numroc(D, nbr_Atilde, myrow, 0, nprow)
    l_dimc_Atilde = numroc(M, nbc_Atilde, mycol, 0, npcol)
    allocate(sca_Atilde(max(1, l_dimr_Atilde), max(1, l_dimc_Atilde)))
    call descinit(desc_Atilde, D, M, nbr_Atilde, nbc_Atilde, 0, 0, context, &
                  max(1, l_dimr_Atilde), sca_info)

    call build_sca_Atilde(sca_AA, desc_AA, D, M, nbr_AA,              &
                          nu_max, block_part, alpha,                   &
                          sca_Atilde, desc_Atilde,                     &
                          l_dimr_Atilde, l_dimc_Atilde,                &
                          nbr_Atilde, nbc_Atilde,                      &
                          myrow, mycol, nprow, npcol)

    ! ========================================================================
    ! Allocate w (replicated on all procs)
    ! ========================================================================
    allocate(w(D), w_prev(D))
    w(:) = 0.0d0

    ! Allocate alpha-step work arrays
    allocate(bvec_local(M, nu_max), bvec(M, nu_max))
    allocate(BtB(nu_max, nu_max), Bty(nu_max))
    allocate(BtB_local(nu_max, nu_max), Bty_local(nu_max))

    ! ========================================================================
    ! Main ALS loop (section 1.10.2, steps B.1 – B.3)
    ! For alpha_method == 0 (fixed alpha), the loop runs a single w-step.
    ! ========================================================================
    do s = 0, nsteps - 1

      w_prev(:) = w(:)
      alpha_prev(:) = alpha(:)

      ! ---- B.1  w-step: distributed ridge regression on Atilde ----
      !   w <- (Atilde^T Atilde + Lambda)^{-1} Atilde^T y
      call ridge_solve_scalapack(sca_Atilde, desc_Atilde, D, M,          &
                                 nbr_Atilde, nbc_Atilde,                 &
                                 sca_yy, desc_yy,                        &
                                 lambda, w,                              &
                                 context, nprow, npcol, myrow, mycol,    &
                                 sca_info)
      if (sca_info /= 0) then
        info = 1
        return
      end if

      ! If fixed alpha, we are done after one pass
      if (alpha_method == 0) exit

      ! ---- B.2  alpha-step: small regression ----
      !   b^(nu) = A^(nu) w^(nu)   for each block nu
      !   B = [b^(1) ... b^(nu_max)]  (M x nu_max)
      !   alpha <- (B^T B + rho I)^{-1} B^T y
      !
      ! Each b^(nu)_j = sum_{d in block nu} AA(d,j) * w(d)
      ! We compute local partial sums using the distributed sca_AA,
      ! then MPI_ALLREDUCE to get full b^(nu).
      call compute_bvec_distributed(sca_AA, desc_AA, D, M, nbr_AA, nbc_AA, &
                                    nu_max, block_part, w,                  &
                                    context, nprow, npcol, myrow, mycol,    &
                                    bvec)

      ! B^T B  and  B^T y  (these are small nu_max × nu_max, done serially)
      ! B is (M, nu_max), y is (M)
      call dgemm('T', 'N', nu_max, nu_max, M, 1.0d0, bvec, M, bvec, M, 0.0d0, BtB, nu_max)

      ! Gather y to a replicated vector for the alpha sub-problem
      ! y is already implicitly available through bvec; we need B^T y
      ! We gather y from the distributed sca_yy
      call compute_Bty(bvec, M, nu_max, sca_yy, desc_yy, D, &
                       context, nprow, npcol, myrow, mycol, Bty)

      ! add rho*I to B^T B
      do i = 1, nu_max
        BtB(i, i) = BtB(i, i) + rho
      end do

      ! Solve for alpha (small system, serial on all procs)
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
          call bk_solver%solve_gram_warm(BtB, Bty, alpha, sca_info, passive=bk_passive)
        else
          ! Lawson-Hanson NNLS (cold start)
          allocate(alpha_phi(nu_max, nu_max))
          alpha_phi(:,:) = BtB(:,:)
          call serial_nnls(alpha_phi, nu_max, nu_max, Bty, alpha, sca_info)
          deallocate(alpha_phi)
        end if
      else
        ! standard solve: (BtB) alpha = Bty
        allocate(alpha_phi(nu_max, nu_max), alpha_rhs(nu_max, 1))
        alpha_phi(:,:) = BtB(:,:)
        alpha_rhs(:, 1) = Bty(:)
        call solve_small_system_local(alpha_phi, nu_max, alpha_rhs, alpha, sca_info)
        deallocate(alpha_phi, alpha_rhs)
      end if

      if (sca_info /= 0) then
        call log_info("ALS-Ridge ScaLAPACK: alpha solve failed at step " // vtoa(s) // " info=" // vtoa(sca_info))
        info = 2
        return
      end if

      ! ---- Log alpha at each step ----
      call log_info("ALS-Ridge ScaLAPACK: step " // vtoa(s) // " alpha = [" // vtoa(alpha(1:nu_max)) // "]")

      ! ---- Rebuild distributed Atilde with updated alpha ----
      call build_sca_Atilde(sca_AA, desc_AA, D, M, nbr_AA,              &
                            nu_max, block_part, alpha,                   &
                            sca_Atilde, desc_Atilde,                     &
                            l_dimr_Atilde, l_dimc_Atilde,                &
                            nbr_Atilde, nbc_Atilde,                      &
                            myrow, mycol, nprow, npcol)

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

    ! ========================================================================
    ! Output: the ALS routine delivers both w^(nu) and alpha_nu.
    ! For extrapolation the effective parameters are the product
    !   w_tilde^(nu) = alpha_nu * w^(nu)      (sec. 1.10.1 of the PDF)
    ! because at prediction time the raw descriptors A^(nu) are used
    ! (not the scaled Atilde^(nu) = alpha_nu * A^(nu)).
    ! ========================================================================
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
    if (allocated(sca_Atilde)) deallocate(sca_Atilde)
    if (allocated(sca_phi))    deallocate(sca_phi)
    if (allocated(sca_rhs))    deallocate(sca_rhs)
    deallocate(frob_nu, lambda, alpha, alpha_prev)
    deallocate(w, w_prev)
    deallocate(bvec_local, bvec, BtB, Bty, BtB_local, Bty_local)

  end subroutine para_als_ridge_solve


  !-----------------------------------------------------------------------------
  ! ridge_solve_scalapack
  !
  ! Solve  (Atilde^T Atilde + diag(lambda)) w = Atilde^T y
  ! using ScaLAPACK pdposv (Cholesky, symmetric positive definite).
  !
  ! Steps:
  !   1. Allocate distributed phi (D×D) and rhs (D×1)
  !   2. phi = Atilde^T * Atilde  via pdgemm
  !   3. Add diagonal lambda to phi using pdelset
  !   4. rhs = Atilde^T * y       via pdgemm (D×M * M×1 = D×1, with 'T')
  !   5. Solve phi * w = rhs     via pdposv
  !   6. Gather w to all procs   via pdelget
  !-----------------------------------------------------------------------------
  subroutine ridge_solve_scalapack(sca_Atilde, desc_Atilde, D, M,        &
                                   nbr_Atilde, nbc_Atilde,               &
                                   sca_yy, desc_yy,                      &
                                   lambda, w,                            &
                                   context, nprow, npcol, myrow, mycol,  &
                                   info)
    use module_scalapack_tools, only: blockset_ml, pdgemm_nn
    implicit none

    integer,  intent(in)    :: D, M
    real(kind_double), intent(in)    :: sca_Atilde(:,:)
    integer,  intent(in)    :: desc_Atilde(9)
    integer,  intent(in)    :: nbr_Atilde, nbc_Atilde
    real(kind_double), intent(in)    :: sca_yy(:,:)
    integer,  intent(in)    :: desc_yy(9)
    real(kind_double), intent(in)    :: lambda(D)
    real(kind_double), intent(out)   :: w(D)
    integer,  intent(in)    :: context, nprow, npcol, myrow, mycol
    integer,  intent(out)   :: info

    ! Local ScaLAPACK matrices
    real(kind_double), allocatable :: sca_phi(:,:)   ! D × D
    integer :: desc_phi(9)
    integer :: l_dimr_phi, l_dimc_phi
    integer :: nbr_phi, nbc_phi

    real(kind_double), allocatable :: sca_rhs(:,:)   ! D × 1
    integer :: desc_rhs(9)
    integer :: l_dimr_rhs, l_dimc_rhs
    integer :: nbr_rhs, nbc_rhs

    integer :: ii, sca_info
    integer :: numroc
    real(kind_double), parameter :: one = 1.0d0, zero = 0.0d0
    integer :: n_retry
    real(kind_double) :: jitter_diag
    real(kind_double), parameter :: jitter0 = 1.0d-14
    real(kind_double), parameter :: jitter_max = 1.0d-4

    info = 0

    ! ---- Allocate phi (D × D) ----
    ! Use square block size for pdposv
    call blockset_ml(nbr_phi, nbr_Atilde, D, nprow, npcol)
    nbc_phi = nbr_phi  ! pdposv requires square block decomposition
    l_dimr_phi = numroc(D, nbr_phi, myrow, 0, nprow)
    l_dimc_phi = numroc(D, nbc_phi, mycol, 0, npcol)
    allocate(sca_phi(max(1, l_dimr_phi), max(1, l_dimc_phi)))
    sca_phi(:,:) = zero
    call descinit(desc_phi, D, D, nbr_phi, nbc_phi, 0, 0, context, &
                  max(1, l_dimr_phi), sca_info)

    ! ---- Allocate rhs (D × 1) ----
    nbr_rhs = nbr_phi
    nbc_rhs = 1
    l_dimr_rhs = numroc(D, nbr_rhs, myrow, 0, nprow)
    l_dimc_rhs = numroc(1, nbc_rhs, mycol, 0, npcol)
    allocate(sca_rhs(max(1, l_dimr_rhs), max(1, l_dimc_rhs)))
    sca_rhs(:,:) = zero
    call descinit(desc_rhs, D, 1, nbr_rhs, nbc_rhs, 0, 0, context, &
                  max(1, l_dimr_rhs), sca_info)

    ! ---- Robust solve loop: rebuild and retry pdposv with diagonal jitter ----
    n_retry = 0
    jitter_diag = 0.0d0
    do
      ! ---- phi = Atilde^T * Atilde ----
      ! Atilde is D × M, so Atilde^T * Atilde is  (M^T × D^T) ... 
      ! Actually: Atilde is (D rows, M cols) stored distributed.
      ! phi = Atilde * Atilde^T  gives (D × D) — correct since
      !   the design matrix A is D × M (short-fat), so A A^T is D × D.
      ! BUT in the formulation: w <- (A^T A + Lambda)^{-1} A^T y
      !   where A^T is M × D (the data matrix in standard ML notation).
      ! In MiLaDy, Amat is stored as (D × M), and the normal equation is
      !   phi = Amat * W * Amat^T = (D × M)(M × M)(M × D) = (D × D)
      ! So phi = Atilde * Atilde^T (D × D).
      sca_phi(:,:) = zero
      call pdgemm('N', 'T', D, D, M, one,               &
                  sca_Atilde, 1, 1, desc_Atilde,         &
                  sca_Atilde, 1, 1, desc_Atilde,         &
                  zero, sca_phi, 1, 1, desc_phi)

      ! ---- Add diagonal regularisation ----
      ! phi(ii,ii) += lambda(ii): uses local indexing for efficiency
      call add_diagonal_to_sca_matrix(sca_phi, desc_phi, D, lambda, &
                                       nbr_phi, myrow, mycol, nprow, npcol)

      ! Additional safety jitter for near-singular SPD cases
      if (jitter_diag > 0.0d0) then
        call add_scalar_diagonal_to_sca_matrix(sca_phi, D, jitter_diag, &
                                               nbr_phi, myrow, mycol, nprow, npcol)
      end if

      ! ---- rhs = Atilde * y ----
      ! Atilde is (D × M), y is (M × 1), rhs is (D × 1)
      sca_rhs(:,:) = zero
      call pdgemm_nn(D, 1, M, one, sca_Atilde, desc_Atilde, &
                     sca_yy, desc_yy, &
                     zero, sca_rhs, desc_rhs)

      ! ---- Solve phi * w = rhs via pdposv (Cholesky) ----
      call blacs_barrier(context, 'A')
      call pdposv('L', D, 1, sca_phi, 1, 1, desc_phi, sca_rhs, 1, 1, desc_rhs, sca_info)
      if (sca_info == 0) exit

      if (myrow == 0 .and. mycol == 0) then
        write(6, '("ALS-Ridge ScaLAPACK: pdposv failed, info = ", i8, " (retry=", i3, ", jitter=", ES12.4E3, ")")') &
             sca_info, n_retry, jitter_diag
      end if

      if (sca_info > 0 .and. jitter_diag < jitter_max) then
        if (jitter_diag == 0.0d0) then
          jitter_diag = jitter0
        else
          jitter_diag = 10.0d0 * jitter_diag
        end if
        n_retry = n_retry + 1
        cycle
      end if

      info = sca_info
      deallocate(sca_phi, sca_rhs)
      return
    end do

    ! ---- Gather w from distributed rhs to replicated vector ----
    w(:) = 0.0d0
    do ii = 1, D
      call pdelget('A', ' ', w(ii), sca_rhs, ii, 1, desc_rhs)
    end do

    deallocate(sca_phi, sca_rhs)

  end subroutine ridge_solve_scalapack


  !-----------------------------------------------------------------------------
  ! add_scalar_diagonal_to_sca_matrix
  !
  ! Add a scalar value to the diagonal entries of a distributed ScaLAPACK matrix.
  !-----------------------------------------------------------------------------
  subroutine add_scalar_diagonal_to_sca_matrix(sca_A, N, diag_shift, &
                                               nb, myrow, mycol, nprow, npcol)
    implicit none
    integer, intent(in) :: N, nb, myrow, mycol, nprow, npcol
    real(kind_double), intent(inout) :: sca_A(:,:)
    real(kind_double), intent(in) :: diag_shift

    integer :: ii, li, lj, rsrc, csrc
    integer :: INDXG2L, INDXG2P

    if (diag_shift == 0.0d0) return

    do ii = 1, N
      rsrc = INDXG2P(ii, nb, myrow, 0, nprow)
      csrc = INDXG2P(ii, nb, mycol, 0, npcol)
      if (rsrc == myrow .and. csrc == mycol) then
        li = INDXG2L(ii, nb, myrow, 0, nprow)
        lj = INDXG2L(ii, nb, mycol, 0, npcol)
        sca_A(li, lj) = sca_A(li, lj) + diag_shift
      end if
    end do

  end subroutine add_scalar_diagonal_to_sca_matrix


  !-----------------------------------------------------------------------------
  ! add_diagonal_to_sca_matrix
  !
  ! Add lambda(i) to the diagonal element (i,i) of a ScaLAPACK matrix.
  ! Uses local indexing for efficiency: only modifies locally-owned diagonal
  ! elements without MPI communication.
  !-----------------------------------------------------------------------------
  subroutine add_diagonal_to_sca_matrix(sca_A, desc_A, N, lambda, &
                                         nb, myrow, mycol, nprow, npcol)
    implicit none
    integer,  intent(in)    :: N, nb, myrow, mycol, nprow, npcol
    real(kind_double), intent(inout) :: sca_A(:,:)
    integer,  intent(in)    :: desc_A(9)
    real(kind_double), intent(in)    :: lambda(N)

    integer :: ii, li, lj, rsrc, csrc
    integer :: INDXG2L, INDXG2P

    ! For each global diagonal element (ii, ii), check if it is
    ! locally owned by this process and if so, add lambda(ii).
    do ii = 1, N
      rsrc = INDXG2P(ii, nb, myrow, 0, nprow)
      csrc = INDXG2P(ii, nb, mycol, 0, npcol)
      if (rsrc == myrow .and. csrc == mycol) then
        li = INDXG2L(ii, nb, myrow, 0, nprow)
        lj = INDXG2L(ii, nb, mycol, 0, npcol)
        sca_A(li, lj) = sca_A(li, lj) + lambda(ii)
      end if
    end do

  end subroutine add_diagonal_to_sca_matrix


  !-----------------------------------------------------------------------------
  ! build_sca_Atilde
  !
  ! Build the distributed scaled design matrix:
  !   Atilde(i, :) = alpha(b) * AA(i, :)   for i in block b
  !
  ! This is done locally: each process scales its own tile rows.
  !-----------------------------------------------------------------------------
  subroutine build_sca_Atilde(sca_AA, desc_AA, D, M, nbr_AA,            &
                              nu_max, block_part, alpha,                 &
                              sca_Atilde, desc_Atilde,                   &
                              l_dimr_Atilde, l_dimc_Atilde,              &
                              nbr_Atilde, nbc_Atilde,                    &
                              myrow, mycol, nprow, npcol)
    implicit none
    integer,  intent(in) :: D, M, nu_max, nbr_AA
    integer,  intent(in) :: nbr_Atilde, nbc_Atilde
    integer,  intent(in) :: l_dimr_Atilde, l_dimc_Atilde
    real(kind_double), intent(in) :: sca_AA(:,:)
    integer,  intent(in) :: desc_AA(9)
    integer,  intent(in) :: block_part(2*nu_max)
    real(kind_double), intent(in) :: alpha(nu_max)
    real(kind_double), intent(out) :: sca_Atilde(:,:)
    integer,  intent(in) :: desc_Atilde(9)
    integer,  intent(in) :: myrow, mycol, nprow, npcol

    integer :: il, jl, gi, b, ib, ie
    integer :: INDXL2G
    real(kind_double) :: scale_factor

    ! Atilde has the same layout as AA (same D, M, same block factors)
    ! For each local row il, find the global row gi, determine which block
    ! it belongs to, and scale accordingly.
    do jl = 1, l_dimc_Atilde
      do il = 1, l_dimr_Atilde
        gi = INDXL2G(il, nbr_Atilde, myrow, 0, nprow)
        ! find which block this global row belongs to
        scale_factor = 0.0d0
        do b = 1, nu_max
          ib = block_part(2*b - 1)
          ie = block_part(2*b)
          if (gi >= ib .and. gi <= ie) then
            scale_factor = alpha(b)
            exit
          end if
        end do
        sca_Atilde(il, jl) = scale_factor * sca_AA(il, jl)
      end do
    end do

  end subroutine build_sca_Atilde


  !-----------------------------------------------------------------------------
  ! compute_distributed_frobenius
  !
  ! Compute squared Frobenius norm of each block of the distributed matrix.
  ! Each process computes partial sums over its local tile, then we
  ! MPI_ALLREDUCE (via DGSUM2D) to get the global result.
  !-----------------------------------------------------------------------------
  subroutine compute_distributed_frobenius(sca_AA, desc_AA, D, M,        &
                                           nbr_AA, nbc_AA,              &
                                           nu_max, block_part, context, &
                                           myrow, mycol, nprow, npcol,  &
                                           frob_nu)
    implicit none
    integer,  intent(in) :: D, M, nu_max, nbr_AA, nbc_AA
    real(kind_double), intent(in) :: sca_AA(:,:)
    integer,  intent(in) :: desc_AA(9)
    integer,  intent(in) :: block_part(2*nu_max)
    integer,  intent(in) :: context, myrow, mycol, nprow, npcol
    real(kind_double), intent(out) :: frob_nu(nu_max)

    real(kind_double), allocatable :: frob_local(:)
    integer :: il, jl, gi, b, ib, ie
    integer :: l_dimr, l_dimc
    integer :: numroc, INDXL2G

    l_dimr = numroc(D, nbr_AA, myrow, 0, nprow)
    l_dimc = numroc(M, nbc_AA, mycol, 0, npcol)

    allocate(frob_local(nu_max))
    frob_local(:) = 0.0d0

    do jl = 1, l_dimc
      do il = 1, l_dimr
        gi = INDXL2G(il, nbr_AA, myrow, 0, nprow)
        ! find block
        do b = 1, nu_max
          ib = block_part(2*b - 1)
          ie = block_part(2*b)
          if (gi >= ib .and. gi <= ie) then
            frob_local(b) = frob_local(b) + sca_AA(il, jl)**2
            exit
          end if
        end do
      end do
    end do

    ! Global sum via BLACS
    frob_nu(:) = frob_local(:)
    call dgsum2d(context, 'A', ' ', nu_max, 1, frob_nu, nu_max, -1, -1)

    deallocate(frob_local)

  end subroutine compute_distributed_frobenius


  !-----------------------------------------------------------------------------
  ! compute_svd_alpha
  !
  ! Compute SVD-based alpha: alpha_b = 1 / max_singular_value(A^(b))
  ! Since each sub-block has Db rows (small) and M columns, we gather
  ! each sub-block to process (0,0), compute SVD serially, and broadcast.
  !-----------------------------------------------------------------------------
  subroutine compute_svd_alpha(sca_AA, desc_AA, D, M, nbr_AA, nbc_AA, &
                               nu_max, block_part, context,            &
                               myrow, mycol, nprow, npcol, alpha)
    implicit none
    integer,  intent(in) :: D, M, nu_max, nbr_AA, nbc_AA
    real(kind_double), intent(in) :: sca_AA(:,:)
    integer,  intent(in) :: desc_AA(9)
    integer,  intent(in) :: block_part(2*nu_max)
    integer,  intent(in) :: context, myrow, mycol, nprow, npcol
    real(kind_double), intent(out) :: alpha(nu_max)

    real(kind_double), allocatable :: Ablock(:,:), S(:), work_svd(:)
    integer :: b, ib, ie, Db, mn, lwork, ii, jj, sca_info
    real(kind_double) :: val

    do b = 1, nu_max
      ib = block_part(2*b - 1)
      ie = block_part(2*b)
      Db = ie - ib + 1

      if (Db <= 0) then
        alpha(b) = 1.0d0
        cycle
      end if

      mn = min(Db, M)

      ! Gather sub-block A^(b) to rank (0,0) via pdelget (scope='A')
      ! Db is small (number of descriptor components per body order),
      ! so the element-by-element gather is acceptable.
      if (myrow == 0 .and. mycol == 0) then
        allocate(Ablock(Db, M), S(mn))
      end if

      ! All procs participate in pdelget with scope 'A'
      do jj = 1, M
        do ii = 1, Db
          call pdelget('A', ' ', val, sca_AA, ib + ii - 1, jj, desc_AA)
          if (myrow == 0 .and. mycol == 0) then
            Ablock(ii, jj) = val
          end if
        end do
      end do

      ! Compute SVD on rank (0,0), then broadcast alpha(b)
      alpha(b) = 1.0d0
      if (myrow == 0 .and. mycol == 0) then
        allocate(work_svd(1))
        lwork = -1
        call dgesvd('N', 'N', Db, M, Ablock, Db, S, Ablock, Db, Ablock, Db, &
                    work_svd, lwork, sca_info)
        lwork = int(work_svd(1)) + 2
        deallocate(work_svd)
        allocate(work_svd(lwork))
        call dgesvd('N', 'N', Db, M, Ablock, Db, S, Ablock, Db, Ablock, Db, &
                    work_svd, lwork, sca_info)
        deallocate(work_svd)

        if (sca_info == 0 .and. mn > 0 .and. S(1) > 0.0d0) then
          alpha(b) = 1.0d0 / S(1)
        end if
        deallocate(Ablock, S)
      end if
    end do

    ! Broadcast alpha from rank (0,0) to all processes via dgebs2d/dgebr2d
    if (myrow == 0 .and. mycol == 0) then
      call dgebs2d(context, 'A', ' ', nu_max, 1, alpha, nu_max)
    else
      call dgebr2d(context, 'A', ' ', nu_max, 1, alpha, nu_max, 0, 0)
    end if

  end subroutine compute_svd_alpha


  !-----------------------------------------------------------------------------
  ! compute_bvec_distributed
  !
  ! Compute b^(nu)_j = sum_{d in block nu} AA(d,j) * w(d)  for j=1..M
  ! using the distributed sca_AA. Each process computes its local partial
  ! sum, then DGSUM2D gives the global result (replicated on all procs).
  !
  ! w is replicated, so each proc knows all w(d).
  !-----------------------------------------------------------------------------
  subroutine compute_bvec_distributed(sca_AA, desc_AA, D, M,              &
                                      nbr_AA, nbc_AA,                    &
                                      nu_max, block_part, w,             &
                                      context, nprow, npcol, myrow, mycol, &
                                      bvec)
    implicit none
    integer,  intent(in) :: D, M, nu_max, nbr_AA, nbc_AA
    real(kind_double), intent(in) :: sca_AA(:,:)
    integer,  intent(in) :: desc_AA(9)
    integer,  intent(in) :: block_part(2*nu_max)
    real(kind_double), intent(in) :: w(D)
    integer,  intent(in) :: context, nprow, npcol, myrow, mycol
    real(kind_double), intent(out) :: bvec(M, nu_max)

    real(kind_double), allocatable :: bvec_local(:,:)
    integer :: il, jl, gi, gj, b, ib, ie
    integer :: l_dimr, l_dimc
    integer :: numroc, INDXL2G

    l_dimr = numroc(D, nbr_AA, myrow, 0, nprow)
    l_dimc = numroc(M, nbc_AA, mycol, 0, npcol)

    allocate(bvec_local(M, nu_max))
    bvec_local(:,:) = 0.0d0

    ! Accumulate local contributions: for each local element (il, jl)
    ! of sca_AA, the global index is (gi, gj).
    ! Contribution to b^(b)_{gj} += AA(gi, gj) * w(gi)  if gi in block b.
    do jl = 1, l_dimc
      gj = INDXL2G(jl, nbc_AA, mycol, 0, npcol)
      do il = 1, l_dimr
        gi = INDXL2G(il, nbr_AA, myrow, 0, nprow)
        ! find block for global row gi
        do b = 1, nu_max
          ib = block_part(2*b - 1)
          ie = block_part(2*b)
          if (gi >= ib .and. gi <= ie) then
            bvec_local(gj, b) = bvec_local(gj, b) + sca_AA(il, jl) * w(gi)
            exit
          end if
        end do
      end do
    end do

    ! Global sum
    bvec(:,:) = bvec_local(:,:)
    call dgsum2d(context, 'A', ' ', M, nu_max, bvec, M, -1, -1)

    deallocate(bvec_local)

  end subroutine compute_bvec_distributed


  !-----------------------------------------------------------------------------
  ! compute_Bty
  !
  ! Compute B^T y where B is the replicated (M, nu_max) matrix of b-vectors
  ! and y is the distributed target vector sca_yy.
  !
  ! First gather y to a replicated vector, then compute B^T y serially.
  !-----------------------------------------------------------------------------
  subroutine compute_Bty(bvec, M, nu_max, sca_yy, desc_yy, D_global, &
                         context, nprow, npcol, myrow, mycol, Bty)
    implicit none
    integer,  intent(in) :: M, nu_max, D_global
    real(kind_double), intent(in) :: bvec(M, nu_max)
    real(kind_double), intent(in) :: sca_yy(:,:)
    integer,  intent(in) :: desc_yy(9)
    integer,  intent(in) :: context, nprow, npcol, myrow, mycol
    real(kind_double), intent(out) :: Bty(nu_max)

    real(kind_double), allocatable :: y_full(:)
    integer :: ii

    ! sca_yy is distributed as (M × 1) — gather to replicated vector
    ! Note: in MiLaDy, the target vector y has dimension M (= dim_data_train)
    ! but it's stored in sca_yy which has global dimension = M.
    ! We use pdelget with scope 'A' to gather all elements (M is the
    ! number of observations, not too large for the alpha sub-problem).
    allocate(y_full(M))
    do ii = 1, M
      call pdelget('A', ' ', y_full(ii), sca_yy, ii, 1, desc_yy)
    end do

    ! B^T y  where B is (M, nu_max)
    call dgemv('T', M, nu_max, 1.0d0, bvec, M, y_full, 1, 0.0d0, Bty, 1)

    deallocate(y_full)

  end subroutine compute_Bty


  !-----------------------------------------------------------------------------
  ! solve_small_system_local
  !
  ! Solve phi * x = rhs for x, where phi is (n, n) and rhs is (n, 1).
  ! Uses LAPACK dsysv. Runs identically on all procs (replicated data).
  !-----------------------------------------------------------------------------
  subroutine solve_small_system_local(phi, n, rhs, x, info)
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
  end subroutine solve_small_system_local

end module module_para_als
