! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   MiLaDy (Machine Learning Dynamics) 
! HND X   
! HND X
! HND X   Milady was written and designed by Mihai-Cosmin Marinica, Alexandra M. Goryaeva 
! HND X   Copyright 2015-2026.
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

!> @brief Low-rank chemical compression of ACE radial pair channels.
!>
!> For ACE with HSVD radials (ace_radial_chem = 3), each (k, l) channel has
!> S×S pair radial functions g_k^{lμν}(r). When S (number of species) is large
!> (25–108), this becomes prohibitive. This module compresses each channel to
!> rank Q via symmetric ALS factorization:
!>
!>   g_k^{lμν}(r_i) ≈ Σ_q A_{μq}^{kl} A_{νq}^{kl} u_q^{kl}(r_i)
!>
!> Storage drops from O(Nr S²) to O(SQ + Nr Q) per channel.
!>
!> Reference: MiladyNoteTechnique5.pdf, section 1.10.3
!>
!> Key memory constraint: the full g_k^{μν}(r_i) tensor is never formed.
!> All normal equations are accumulated column-by-column (loop over ν).
module module_chem_compress
  use iso_fortran_env, only: dp => real64
  use mld_logger
  use mld_mpi, only: mld_rank, mld_size, mld_mpi_abort, comm_mld
  use module_spline_interpolation_mine, only: CubicSpline
  implicit none
  private

  public :: chem_compress_t, chem_compressor

  type :: chem_compress_t
     integer :: S = 0             !< number of chemical species
     integer :: Q = 0             !< compression rank
     integer :: Nr = 0            !< radial grid points (= ace_npoints_spline)
     integer :: kmax = 0          !< max k index
     integer :: lmax_p1 = 0       !< lmax + 1
     integer :: n_channels = 0    !< total channels = kmax * (lmax+1)
     integer :: niter = 50        !< ALS iterations
     real(dp) :: lambda = 1.0d-6  !< ridge regularization

     !> Mixing matrices A(S, Q, n_channels) — one per (k,l) channel
     real(dp), allocatable :: A(:,:,:)

     !> Prototype radial splines u_spline(Q, n_channels) — for evaluation at arbitrary r
     type(CubicSpline), allocatable :: u_spline(:,:)

     logical :: initialized = .false.

   contains
     procedure :: init       => chem_compress_init
     procedure :: compute    => chem_compress_compute
      procedure :: compress_channel => chem_compress_channel
     procedure :: allreduce_compressed => chem_compress_allreduce
     procedure :: eval_rad   => chem_compress_eval_rad
     procedure :: eval_drad  => chem_compress_eval_drad
     procedure :: destroy    => chem_compress_destroy
  end type chem_compress_t

  !> Module-level instance (like radialace in module_ace_radial)
  type(chem_compress_t) :: chem_compressor

contains

  !> Initialize the compressor
  subroutine chem_compress_init(this, S, Q, Nr, kmax, lmax_p1, niter, lambda, r_cut_in, r_cut_out)
    class(chem_compress_t), intent(inout) :: this
    integer, intent(in) :: S, Q, Nr, kmax, lmax_p1, niter
    real(dp), intent(in) :: lambda
    real(dp), intent(in), optional :: r_cut_in, r_cut_out

    integer :: ich, qq
    real(dp), allocatable :: xg_dummy(:), yfunc(:), d_yfunc(:)

    _NAMECURRENT_("chem_compress_init")
    _MLD_BEGIN_

    this%S = S
    this%Q = min(Q, S)   ! Q cannot exceed S
    this%Nr = Nr
    this%kmax = kmax
    this%lmax_p1 = lmax_p1
    this%n_channels = kmax * lmax_p1
    this%niter = niter
    this%lambda = lambda

    if (allocated(this%A))        deallocate(this%A)
    if (allocated(this%u_spline)) deallocate(this%u_spline)

    allocate(this%A(S, this%Q, this%n_channels), source=0.0_dp)
    allocate(this%u_spline(this%Q, this%n_channels))

    ! Initialize ALL prototype splines with zero coefficients so that the
    ! allreduce buffer is valid even for non-owned channels.
    if (present(r_cut_in) .and. present(r_cut_out)) then
      do ich = 1, this%n_channels
         do qq = 1, this%Q
            call this%u_spline(qq, ich)%init(r_cut_in, r_cut_out, Nr, xg_dummy, yfunc, d_yfunc)
            yfunc(:) = 0.0_dp
            d_yfunc(:) = 0.0_dp
            call this%u_spline(qq, ich)%compute(yfunc, d_yfunc)
            if (allocated(xg_dummy)) deallocate(xg_dummy)
            if (allocated(yfunc))    deallocate(yfunc)
            if (allocated(d_yfunc))  deallocate(d_yfunc)
         end do
      end do
    end if

    this%initialized = .true.

    call log_info("CHEM_COMPRESS: init S="//vtoa(S)//" Q="//vtoa(this%Q)// &
                  " Nr="//vtoa(Nr)//" n_channels="//vtoa(this%n_channels)// &
                  " niter="//vtoa(niter)//" lambda="//vtoa(lambda))
    call log_info("CHEM_COMPRESS: compression factor ~ "// &
                  vtoa(dble(S*S)/dble(this%Q))//"x per channel")

    _MLD_END_
  end subroutine chem_compress_init


  !> Run ALS compression for all (k,l) channels.
  !>
  !> Reads g_k^{μν}(r_i) from existing splines column-by-column (loop over ν),
  !> runs ALS, and stores A(S,Q,ich) + Q prototype CubicSplines per channel.
  !>
  !> After this call, the compressed representation can be evaluated at any r
  !> via eval_rad / eval_drad.
  subroutine chem_compress_compute(this, ic_rad, radial_sp, xgrid, r_cut_in, r_cut_out)
    implicit none
    class(chem_compress_t), intent(inout) :: this
    integer, intent(in) :: ic_rad(:,:,:,0:)  ! (dim_mu, dim_mu, 1:kmax, 0:lmax)
    type(CubicSpline), intent(in) :: radial_sp(:)
    real(dp), intent(in) :: xgrid(:)
    real(dp), intent(in) :: r_cut_in, r_cut_out

    ! Local variables
    integer :: kk, ll, ich, iter, mu, nu, qq, ip, ispline, itask, mm
    real(dp), allocatable :: A_loc(:,:)       ! (S, Q) current channel mixing matrix
    real(dp), allocatable :: u_loc(:,:)       ! (Nr, Q) current channel prototype radials on grid
    real(dp), allocatable :: du_loc(:,:)      ! (Nr, Q) derivatives on grid
    real(dp), allocatable :: f_col(:,:)       ! (Nr, S) g^{:,ν}(r) for one ν, all r_i
    real(dp), allocatable :: df_col(:,:)      ! (Nr, S) dg^{:,ν}/dr for one ν, all r_i
    real(dp), allocatable :: AtA(:,:)         ! (Q, Q) = A^T A
    real(dp), allocatable :: BtB(:,:)         ! (Q, Q) = (A^T A) ⊙ (A^T A) (Hadamard)
    real(dp), allocatable :: Btv(:,:)         ! (Q, Nr) = B^T vec(R_i) for all i
    real(dp), allocatable :: w(:)             ! (Q) scratch
    real(dp), allocatable :: CtC(:,:)         ! (Q, Q) for A-update
    real(dp), allocatable :: XC(:,:)          ! (S, Q) for A-update
    real(dp), allocatable :: col_contrib(:,:) ! (Nr, Q) scratch
    real(dp), allocatable :: yfunc(:), d_yfunc(:), xg_dummy(:)  ! for spline init/compute
    real(dp), allocatable :: xg_safe(:)       ! clamped copy of xgrid
    real(dp) :: col_norm
    integer :: info_lap
    integer, allocatable :: ipiv(:)

    ! For MPI broadcast of spline coefficients
    integer :: n_proto_splines, isp_global
    real(dp), allocatable :: spline_buf(:)
    real(dp), allocatable :: A_buf(:)

    _NAMECURRENT_("chem_compress_compute")
    _MLD_BEGIN_

    if (.not. this%initialized) then
       call log_critical("CHEM_COMPRESS: compute called before init")
       call mld_mpi_abort("CHEM_COMPRESS: compute called before init")
    end if

    allocate(A_loc(this%S, this%Q))
    allocate(u_loc(this%Nr, this%Q),  source=0.0_dp)
    allocate(du_loc(this%Nr, this%Q), source=0.0_dp)
    allocate(f_col(this%Nr, this%S))
    allocate(df_col(this%Nr, this%S))
    allocate(AtA(this%Q, this%Q))
    allocate(BtB(this%Q, this%Q))
    allocate(Btv(this%Q, this%Nr))
    allocate(w(this%Q))
    allocate(CtC(this%Q, this%Q))
    allocate(XC(this%S, this%Q))
    allocate(col_contrib(this%Nr, this%Q))
    allocate(ipiv(this%Q))

    ! Build a clamped copy of xgrid to avoid floating-point boundary issues.
    ! The accumulated grid x(Nr) = r_cut_in + (Nr-1)*delta may overshoot r_cut_out
    ! by an ULP, triggering "rr is outside the range of x" warnings from CubicSpline.
    allocate(xg_safe(this%Nr))
    xg_safe(:) = xgrid(1:this%Nr)
    xg_safe(1) = max(xg_safe(1), r_cut_in)
    xg_safe(this%Nr) = min(xg_safe(this%Nr), r_cut_out)

    ! Initialize all A to zero so allreduce works
    this%A(:,:,:) = 0.0_dp

    ! Initialize ALL prototype splines on ALL ranks (zero coefficients)
    ! so the allreduce buffer is valid even for non-owned channels.
    do ich = 1, this%n_channels
       do qq = 1, this%Q
          call this%u_spline(qq, ich)%init(r_cut_in, r_cut_out, this%Nr, xg_dummy, yfunc, d_yfunc)
          ! init allocates yfunc/d_yfunc/xg_dummy, fill with zeros
          yfunc(:) = 0.0_dp
          d_yfunc(:) = 0.0_dp
          call this%u_spline(qq, ich)%compute(yfunc, d_yfunc)
          if (allocated(xg_dummy)) deallocate(xg_dummy)
          if (allocated(yfunc))    deallocate(yfunc)
          if (allocated(d_yfunc))  deallocate(d_yfunc)
       end do
    end do

    ! Round-robin distribution of channels over MPI ranks
    itask = -1
    do ll = 0, this%lmax_p1 - 1
    do kk = 1, this%kmax
      itask = itask + 1
      ich = ll * this%kmax + kk  ! 1-based channel index

      ! Round-robin: skip channels not owned by this rank
      if (mod(itask, mld_size) /= mld_rank) cycle

      call log_info("CHEM_COMPRESS: channel kk="//vtoa(kk)//" ll="//vtoa(ll)// &
                    " (ich="//vtoa(ich)//") on rank "//vtoa(mld_rank))

      ! === Initialize A from Gram-matrix SVD ===
      call init_A_from_columns(this, ic_rad, radial_sp, xg_safe, kk, ll, A_loc)

      u_loc(:,:) = 0.0_dp
      du_loc(:,:) = 0.0_dp

      ! === ALS iterations ===
      do iter = 1, this%niter

        ! --- Step 1: Update u(r_i) for all grid points ---
        ! BtB = (A^T A) ⊙ (A^T A)  — does not depend on r_i
        call dgemm('T', 'N', this%Q, this%Q, this%S, 1.0_dp, &
                    A_loc, this%S, A_loc, this%S, 0.0_dp, AtA, this%Q)
        BtB(:,:) = AtA(:,:) * AtA(:,:)
        do qq = 1, this%Q
           BtB(qq, qq) = BtB(qq, qq) + this%lambda
        end do

        ! Accumulate Btv column-by-column over ν
        Btv(:,:) = 0.0_dp
        do nu = 1, this%S
           do mu = 1, this%S
              ispline = ic_rad(mu, nu, kk, ll)
              do ip = 1, this%Nr
                 f_col(ip, mu) = radial_sp(ispline)%evaluate(xg_safe(ip))
              end do
           end do
           do ip = 1, this%Nr
              call dgemv('T', this%S, this%Q, 1.0_dp, A_loc, this%S, &
                         f_col(ip, :), 1, 0.0_dp, w, 1)
              do qq = 1, this%Q
                 Btv(qq, ip) = Btv(qq, ip) + A_loc(nu, qq) * w(qq)
              end do
           end do
        end do

        ! Solve BtB * u_i = Btv_i for all grid points at once
        call dsytrf('U', this%Q, BtB, this%Q, ipiv, w, this%Q, info_lap)
        if (info_lap /= 0) call log_warning("CHEM_COMPRESS: dsytrf u-update failed, info="//vtoa(info_lap))
        call dsytrs('U', this%Q, this%Nr, BtB, this%Q, ipiv, Btv, this%Q, info_lap)
        do ip = 1, this%Nr
           u_loc(ip, :) = Btv(:, ip)
        end do

        ! --- Step 2: Update A ---
        ! CtC = (U^T U) ⊙ (A^T A)
        call dgemm('T', 'N', this%Q, this%Q, this%Nr, 1.0_dp, &
                    u_loc, this%Nr, u_loc, this%Nr, 0.0_dp, CtC, this%Q)
        call dgemm('T', 'N', this%Q, this%Q, this%S, 1.0_dp, &
                    A_loc, this%S, A_loc, this%S, 0.0_dp, AtA, this%Q)
        CtC(:,:) = CtC(:,:) * AtA(:,:)
        do qq = 1, this%Q
           CtC(qq, qq) = CtC(qq, qq) + this%lambda
        end do

        ! XC = X_{(2)} C  accumulated column-by-column
        XC(:,:) = 0.0_dp
        do nu = 1, this%S
           do mu = 1, this%S
              ispline = ic_rad(mu, nu, kk, ll)
              do ip = 1, this%Nr
                 f_col(ip, mu) = radial_sp(ispline)%evaluate(xg_safe(ip))
              end do
           end do
           do qq = 1, this%Q
              col_contrib(:, qq) = u_loc(:, qq) * A_loc(nu, qq)
           end do
           call dgemm('T', 'N', this%S, this%Q, this%Nr, 1.0_dp, &
                      f_col, this%Nr, col_contrib, this%Nr, 1.0_dp, XC, this%S)
        end do

        ! Solve CtC * a_mu = xc_mu for each species μ
        call dsytrf('U', this%Q, CtC, this%Q, ipiv, w, this%Q, info_lap)
        if (info_lap /= 0) call log_warning("CHEM_COMPRESS: dsytrf A-update failed, info="//vtoa(info_lap))
        do mu = 1, this%S
           w(1:this%Q) = XC(mu, 1:this%Q)
           call dsytrs('U', this%Q, 1, CtC, this%Q, ipiv, w, this%Q, info_lap)
           A_loc(mu, 1:this%Q) = w(1:this%Q)
        end do

        ! Normalize columns of A, rescale u
        do qq = 1, this%Q
           col_norm = sqrt(sum(A_loc(:, qq)**2))
           if (col_norm > 1.0d-14) then
              u_loc(:, qq) = u_loc(:, qq) * col_norm * col_norm
              A_loc(:, qq) = A_loc(:, qq) / col_norm
           end if
        end do

      end do  ! ALS iter

      ! Compute du_loc (derivative of prototype radials)
      du_loc(:,:) = 0.0_dp
      do nu = 1, this%S
         do mu = 1, this%S
            ispline = ic_rad(mu, nu, kk, ll)
            do ip = 1, this%Nr
               df_col(ip, mu) = radial_sp(ispline)%derivative(xg_safe(ip))
            end do
         end do
         do qq = 1, this%Q
            do ip = 1, this%Nr
               du_loc(ip, qq) = du_loc(ip, qq) + A_loc(nu, qq) * &
                    dot_product(A_loc(:, qq), df_col(ip, :))
            end do
         end do
      end do

      ! Store A
      this%A(:, :, ich) = A_loc(:,:)

      ! Build cubic splines for each prototype radial q
      do qq = 1, this%Q
         call this%u_spline(qq, ich)%init(r_cut_in, r_cut_out, this%Nr, xg_dummy, yfunc, d_yfunc)
         ! init allocates yfunc/d_yfunc/xg_dummy, fill from ALS output
         yfunc(:)   = u_loc(:, qq)
         d_yfunc(:) = du_loc(:, qq)
         call this%u_spline(qq, ich)%compute(yfunc, d_yfunc)
         if (allocated(xg_dummy)) deallocate(xg_dummy)
         if (allocated(yfunc))    deallocate(yfunc)
         if (allocated(d_yfunc))  deallocate(d_yfunc)
      end do

    end do  ! kk
    end do  ! ll

    ! --- MPI broadcast: pack all A and spline coefficients, allreduce ---
    ! A: reduce in-place
    call comm_mld%sum(this%A)

    ! Spline coefficients: pack a/b/c/d for all prototype splines, allreduce
    n_proto_splines = this%Q * this%n_channels
    isp_global = n_proto_splines * (4 * this%Nr - 1)
    allocate(spline_buf(isp_global), source=0.0_dp)

    ! Pack
    mm = 0
    do ich = 1, this%n_channels
       do qq = 1, this%Q
          spline_buf(mm+1                : mm+  this%Nr  ) = this%u_spline(qq,ich)%a(1:this%Nr)
          spline_buf(mm+  this%Nr+1      : mm+2*this%Nr  ) = this%u_spline(qq,ich)%b(1:this%Nr)
          spline_buf(mm+2*this%Nr+1      : mm+3*this%Nr  ) = this%u_spline(qq,ich)%c(1:this%Nr)
          spline_buf(mm+3*this%Nr+1      : mm+4*this%Nr-1) = this%u_spline(qq,ich)%d(1:this%Nr-1)
          mm = mm + (4 * this%Nr - 1)
       end do
    end do

    call comm_mld%sum(spline_buf)

    ! Unpack
    mm = 0
    do ich = 1, this%n_channels
       do qq = 1, this%Q
          this%u_spline(qq,ich)%a(1:this%Nr)   = spline_buf(mm+1                : mm+  this%Nr  )
          this%u_spline(qq,ich)%b(1:this%Nr)   = spline_buf(mm+  this%Nr+1      : mm+2*this%Nr  )
          this%u_spline(qq,ich)%c(1:this%Nr)   = spline_buf(mm+2*this%Nr+1      : mm+3*this%Nr  )
          this%u_spline(qq,ich)%d(1:this%Nr-1) = spline_buf(mm+3*this%Nr+1      : mm+4*this%Nr-1)
          mm = mm + (4 * this%Nr - 1)
       end do
    end do
    deallocate(spline_buf)

    call log_info("CHEM_COMPRESS: ALS compression complete for all "// &
                  vtoa(this%n_channels)//" channels, "//vtoa(n_proto_splines)//" prototype splines built")

    deallocate(A_loc, u_loc, du_loc, f_col, df_col)
    deallocate(AtA, BtB, Btv, w, CtC, XC, col_contrib, ipiv)
    deallocate(xg_safe)

    _MLD_END_
  end subroutine chem_compress_compute


  !> Initialize A for a channel (kk, ll) from a thin SVD of representative columns.
  !>
  !> Strategy: pick the first min(Q, S) species as representative columns ν,
  !> evaluate g^{1:S, ν}(r_1..r_Nr), stack into (Nr*S, Q_rep) matrix,
  !> compute thin SVD, reshape left singular vectors into A(S, Q).
  subroutine init_A_from_columns(this, ic_rad, radial_sp, xgrid, kk, ll, A_loc)
    use module_spline_interpolation_mine, only: CubicSpline
    use module_svd_small_gen_matrix, only: eigen_svd
    implicit none
    class(chem_compress_t), intent(in) :: this
    integer, intent(in) :: ic_rad(:,:,:,0:)  ! (dim_mu, dim_mu, 1:kmax, 0:lmax)
    type(CubicSpline), intent(in) :: radial_sp(:)
    real(dp), intent(in) :: xgrid(:)
    integer, intent(in) :: kk, ll
    real(dp), intent(out) :: A_loc(this%S, this%Q)

    real(dp), allocatable :: G_sample(:,:)  ! (Nr, S) for one ν
    real(dp), allocatable :: GtG(:,:)       ! (S, S) Gram matrix
    real(dp) :: fval, dfval
    integer :: mu, nu, ip, qq, ispline, n_sample
    type(eigen_svd) :: svd_init

    ! Build A^T A ≈ average of g^T g over representative ν columns
    ! This avoids storing the full (Nr*S, S) matrix
    allocate(G_sample(this%Nr, this%S))
    allocate(GtG(this%S, this%S), source=0.0_dp)

    n_sample = min(this%Q, this%S)

    do nu = 1, this%S
       ! Evaluate g^{1:S, ν}(r_i) for all grid points
       do mu = 1, this%S
          ispline = ic_rad(mu, nu, kk, ll)
          do ip = 1, this%Nr
             G_sample(ip, mu) = radial_sp(ispline)%evaluate(xgrid(ip))
          end do
       end do
       ! Accumulate GtG += G_sample^T * G_sample
       call dgemm('T', 'N', this%S, this%S, this%Nr, 1.0_dp, &
                  G_sample, this%Nr, G_sample, this%Nr, 1.0_dp, GtG, this%S)
    end do

    ! SVD of GtG to get the principal directions
    call svd_init%init(this%S, this%S)
    call svd_init%evaluate(GtG)

    ! A = first Q right singular vectors (columns of V)
    ! svd_init%eigvecVT is (S, S), rows are right singular vectors
    do qq = 1, this%Q
       A_loc(:, qq) = svd_init%eigvecVT(qq, :)
    end do

    ! Normalize columns
    do qq = 1, this%Q
       fval = sqrt(sum(A_loc(:, qq)**2))
       if (fval > 1.0d-14) A_loc(:, qq) = A_loc(:, qq) / fval
    end do

    call svd_init%destroy()
    deallocate(G_sample, GtG)

  end subroutine init_A_from_columns


  !> Compress a single (kk, ll) channel from grid values directly.
  !>
  !> This is the on-the-fly variant: the S×S radial values on the grid are
  !> provided as arrays g_grid(Nr, S, S) and dg_grid(Nr, S, S) instead of
  !> being read from CubicSpline objects.
  !>
  !> Runs ALS and stores the result into this%A(:,:,ich) and this%u_spline(:,ich).
  !> This must be called for EACH channel; no allreduce is done here.
  !> Call allreduce_compressed() once after all channels are processed.
  subroutine chem_compress_channel(this, kk, ll, g_grid, dg_grid, xgrid, r_cut_in, r_cut_out)
    implicit none
    class(chem_compress_t), intent(inout) :: this
    integer, intent(in) :: kk, ll
    real(dp), intent(in) :: g_grid(:,:,:)    ! (Nr, S, S) — g_k^{l,mu,nu}(r_i)
    real(dp), intent(in) :: dg_grid(:,:,:)   ! (Nr, S, S) — dg/dr
    real(dp), intent(in) :: xgrid(:)
    real(dp), intent(in) :: r_cut_in, r_cut_out

    ! Local variables
    integer :: ich, iter, mu, nu, qq, ip, mm
    real(dp), allocatable :: A_loc(:,:)       ! (S, Q)
    real(dp), allocatable :: u_loc(:,:)       ! (Nr, Q)
    real(dp), allocatable :: du_loc(:,:)      ! (Nr, Q)
    real(dp), allocatable :: f_col(:,:)       ! (Nr, S) g^{:,ν}(r) for one ν
    real(dp), allocatable :: AtA(:,:)         ! (Q, Q)
    real(dp), allocatable :: BtB(:,:)         ! (Q, Q)
    real(dp), allocatable :: Btv(:,:)         ! (Q, Nr)
    real(dp), allocatable :: w(:)             ! (Q)
    real(dp), allocatable :: CtC(:,:)         ! (Q, Q)
    real(dp), allocatable :: XC(:,:)          ! (S, Q)
    real(dp), allocatable :: col_contrib(:,:) ! (Nr, Q)
    real(dp), allocatable :: W_mat(:,:)       ! (Nr, Q) batched dgemm result
    real(dp), allocatable :: yfunc(:), d_yfunc(:), xg_dummy(:)
    real(dp), allocatable :: GtG(:,:)         ! (S, S) Gram matrix for A init
    real(dp) :: col_norm, fval
    integer :: info_lap
    integer, allocatable :: ipiv(:)

    ! --- A initialisation from Gram-matrix SVD (inlined) ---
    ! We use the same approach as init_A_from_columns but read from g_grid.
    block
      use module_svd_small_gen_matrix, only: eigen_svd
      type(eigen_svd) :: svd_init
      real(dp), allocatable :: G_sample(:,:)

      allocate(G_sample(this%Nr, this%S))
      allocate(GtG(this%S, this%S), source=0.0_dp)
      allocate(A_loc(this%S, this%Q))

      do nu = 1, this%S
         G_sample(:, :) = g_grid(:, :, nu)
         call dgemm('T', 'N', this%S, this%S, this%Nr, 1.0_dp, &
                    G_sample, this%Nr, G_sample, this%Nr, 1.0_dp, GtG, this%S)
      end do

      call svd_init%init(this%S, this%S)
      call svd_init%evaluate(GtG)

      do qq = 1, this%Q
         A_loc(:, qq) = svd_init%eigvecVT(qq, :)
      end do
      do qq = 1, this%Q
         fval = sqrt(sum(A_loc(:, qq)**2))
         if (fval > 1.0d-14) A_loc(:, qq) = A_loc(:, qq) / fval
      end do

      call svd_init%destroy()
      deallocate(G_sample, GtG)
    end block

    ! --- Allocate ALS work arrays ---
    ich = ll * this%kmax + kk

    allocate(u_loc(this%Nr, this%Q),  source=0.0_dp)
    allocate(du_loc(this%Nr, this%Q), source=0.0_dp)
    allocate(f_col(this%Nr, this%S))
    allocate(AtA(this%Q, this%Q))
    allocate(BtB(this%Q, this%Q))
    allocate(Btv(this%Q, this%Nr))
    allocate(w(this%Q))
    allocate(CtC(this%Q, this%Q))
    allocate(XC(this%S, this%Q))
    allocate(col_contrib(this%Nr, this%Q))
    allocate(W_mat(this%Nr, this%Q))
    allocate(ipiv(this%Q))

    ! === ALS iterations ===
    do iter = 1, this%niter

      ! --- Step 1: Update u(r_i) ---
      call dgemm('T', 'N', this%Q, this%Q, this%S, 1.0_dp, &
                  A_loc, this%S, A_loc, this%S, 0.0_dp, AtA, this%Q)
      BtB(:,:) = AtA(:,:) * AtA(:,:)
      do qq = 1, this%Q
         BtB(qq, qq) = BtB(qq, qq) + this%lambda
      end do

      Btv(:,:) = 0.0_dp
      do nu = 1, this%S
         f_col(:, :) = g_grid(:, :, nu)
         ! W_mat(Nr, Q) = f_col(Nr, S) @ A_loc(S, Q) — one dgemm replaces Nr dgemv calls
         call dgemm('N', 'N', this%Nr, this%Q, this%S, 1.0_dp, &
                    f_col, this%Nr, A_loc, this%S, 0.0_dp, W_mat, this%Nr)
         ! Btv(q, ip) += A_loc(nu, q) * W_mat(ip, q)
         do qq = 1, this%Q
            do ip = 1, this%Nr
               Btv(qq, ip) = Btv(qq, ip) + A_loc(nu, qq) * W_mat(ip, qq)
            end do
         end do
      end do

      call dsytrf('U', this%Q, BtB, this%Q, ipiv, w, this%Q, info_lap)
      if (info_lap /= 0) call log_warning("CHEM_COMPRESS: dsytrf u-update failed, info="//vtoa(info_lap))
      call dsytrs('U', this%Q, this%Nr, BtB, this%Q, ipiv, Btv, this%Q, info_lap)
      do ip = 1, this%Nr
         u_loc(ip, :) = Btv(:, ip)
      end do

      ! --- Step 2: Update A ---
      call dgemm('T', 'N', this%Q, this%Q, this%Nr, 1.0_dp, &
                  u_loc, this%Nr, u_loc, this%Nr, 0.0_dp, CtC, this%Q)
      call dgemm('T', 'N', this%Q, this%Q, this%S, 1.0_dp, &
                  A_loc, this%S, A_loc, this%S, 0.0_dp, AtA, this%Q)
      CtC(:,:) = CtC(:,:) * AtA(:,:)
      do qq = 1, this%Q
         CtC(qq, qq) = CtC(qq, qq) + this%lambda
      end do

      XC(:,:) = 0.0_dp
      do nu = 1, this%S
         f_col(:, :) = g_grid(:, :, nu)
         do qq = 1, this%Q
            col_contrib(:, qq) = u_loc(:, qq) * A_loc(nu, qq)
         end do
         call dgemm('T', 'N', this%S, this%Q, this%Nr, 1.0_dp, &
                    f_col, this%Nr, col_contrib, this%Nr, 1.0_dp, XC, this%S)
      end do

      call dsytrf('U', this%Q, CtC, this%Q, ipiv, w, this%Q, info_lap)
      if (info_lap /= 0) call log_warning("CHEM_COMPRESS: dsytrf A-update failed, info="//vtoa(info_lap))
      do mu = 1, this%S
         w(1:this%Q) = XC(mu, 1:this%Q)
         call dsytrs('U', this%Q, 1, CtC, this%Q, ipiv, w, this%Q, info_lap)
         A_loc(mu, 1:this%Q) = w(1:this%Q)
      end do

      ! Normalize columns of A, rescale u
      do qq = 1, this%Q
         col_norm = sqrt(sum(A_loc(:, qq)**2))
         if (col_norm > 1.0d-14) then
            u_loc(:, qq) = u_loc(:, qq) * col_norm * col_norm
            A_loc(:, qq) = A_loc(:, qq) / col_norm
         end if
      end do

    end do  ! ALS iter

    ! Compute du_loc from dg_grid: du(ip,q) = sum_nu A(nu,q) * sum_mu A(mu,q) * dg(ip,mu,nu)
    du_loc(:,:) = 0.0_dp
    do nu = 1, this%S
       f_col(:, :) = dg_grid(:, :, nu)
       ! W_mat(Nr, Q) = f_col(Nr, S) @ A_loc(S, Q)
       call dgemm('N', 'N', this%Nr, this%Q, this%S, 1.0_dp, &
                  f_col, this%Nr, A_loc, this%S, 0.0_dp, W_mat, this%Nr)
       do qq = 1, this%Q
          do ip = 1, this%Nr
             du_loc(ip, qq) = du_loc(ip, qq) + A_loc(nu, qq) * W_mat(ip, qq)
          end do
       end do
    end do

    ! Store A
    this%A(:, :, ich) = A_loc(:,:)

    ! Build cubic splines for each prototype radial q
    do qq = 1, this%Q
       call this%u_spline(qq, ich)%init(r_cut_in, r_cut_out, this%Nr, xg_dummy, yfunc, d_yfunc)
       yfunc(:)   = u_loc(:, qq)
       d_yfunc(:) = du_loc(:, qq)
       call this%u_spline(qq, ich)%compute(yfunc, d_yfunc)
       if (allocated(xg_dummy)) deallocate(xg_dummy)
       if (allocated(yfunc))    deallocate(yfunc)
       if (allocated(d_yfunc))  deallocate(d_yfunc)
    end do

    deallocate(A_loc, u_loc, du_loc, f_col)
    deallocate(AtA, BtB, Btv, w, CtC, XC, col_contrib, W_mat, ipiv)

  end subroutine chem_compress_channel


  !> MPI allreduce of all compressed data (A matrices + prototype spline coefficients).
  !> Call once after all channels have been compressed via compress_channel.
  subroutine chem_compress_allreduce(this)
    implicit none
    class(chem_compress_t), intent(inout) :: this

    integer :: n_proto_splines, isp_global, ich, qq, mm
    real(dp), allocatable :: spline_buf(:)

    _NAMECURRENT_("chem_compress_allreduce")
    _MLD_BEGIN_

    ! A: reduce in-place
    call comm_mld%sum(this%A)

    ! Spline coefficients: pack a/b/c/d, allreduce, unpack
    n_proto_splines = this%Q * this%n_channels
    isp_global = n_proto_splines * (4 * this%Nr - 1)
    allocate(spline_buf(isp_global), source=0.0_dp)

    mm = 0
    do ich = 1, this%n_channels
       do qq = 1, this%Q
          spline_buf(mm+1                : mm+  this%Nr  ) = this%u_spline(qq,ich)%a(1:this%Nr)
          spline_buf(mm+  this%Nr+1      : mm+2*this%Nr  ) = this%u_spline(qq,ich)%b(1:this%Nr)
          spline_buf(mm+2*this%Nr+1      : mm+3*this%Nr  ) = this%u_spline(qq,ich)%c(1:this%Nr)
          spline_buf(mm+3*this%Nr+1      : mm+4*this%Nr-1) = this%u_spline(qq,ich)%d(1:this%Nr-1)
          mm = mm + (4 * this%Nr - 1)
       end do
    end do

    call comm_mld%sum(spline_buf)

    mm = 0
    do ich = 1, this%n_channels
       do qq = 1, this%Q
          this%u_spline(qq,ich)%a(1:this%Nr)   = spline_buf(mm+1                : mm+  this%Nr  )
          this%u_spline(qq,ich)%b(1:this%Nr)   = spline_buf(mm+  this%Nr+1      : mm+2*this%Nr  )
          this%u_spline(qq,ich)%c(1:this%Nr)   = spline_buf(mm+2*this%Nr+1      : mm+3*this%Nr  )
          this%u_spline(qq,ich)%d(1:this%Nr-1) = spline_buf(mm+3*this%Nr+1      : mm+4*this%Nr-1)
          mm = mm + (4 * this%Nr - 1)
       end do
    end do
    deallocate(spline_buf)

    call log_info("CHEM_COMPRESS: allreduce complete for "// &
                  vtoa(this%n_channels)//" channels, "//vtoa(n_proto_splines)//" prototype splines")

    _MLD_END_
  end subroutine chem_compress_allreduce


  !> Evaluate the compressed radial function g_k^{lμν}(r) at arbitrary distance.
  !>
  !> g_k^{lμν}(r) ≈ Σ_q A(μ,q,ich) * A(ν,q,ich) * u_q(r)
  !>
  !> @param[in]  mu    Central atom species index (1..S)
  !> @param[in]  nu    Neighbor atom species index (1..S)
  !> @param[in]  kk    ACE k index (1..kmax)
  !> @param[in]  ll    ACE l index (0..lmax)
  !> @param[in]  rr    Interatomic distance
  !> @return     gval  g_k^{lμν}(r)
  function chem_compress_eval_rad(this, mu, nu, kk, ll, rr) result(gval)
    class(chem_compress_t), intent(in) :: this
    integer, intent(in) :: mu, nu, kk, ll
    real(dp), intent(in) :: rr
    real(dp) :: gval

    integer :: qq, ich
    real(dp) :: coeff

    ich = ll * this%kmax + kk

    gval = 0.0_dp
    do qq = 1, this%Q
       coeff = this%A(mu, qq, ich) * this%A(nu, qq, ich)
       gval = gval + coeff * this%u_spline(qq, ich)%evaluate(rr)
    end do

  end function chem_compress_eval_rad


  !> Evaluate derivative of the compressed radial function d/dr g_k^{lμν}(r).
  !>
  !> @param[in]  mu    Central atom species index (1..S)
  !> @param[in]  nu    Neighbor atom species index (1..S)
  !> @param[in]  kk    ACE k index (1..kmax)
  !> @param[in]  ll    ACE l index (0..lmax)
  !> @param[in]  rr    Interatomic distance
  !> @return     dgval d/dr g_k^{lμν}(r)
  function chem_compress_eval_drad(this, mu, nu, kk, ll, rr) result(dgval)
    class(chem_compress_t), intent(in) :: this
    integer, intent(in) :: mu, nu, kk, ll
    real(dp), intent(in) :: rr
    real(dp) :: dgval

    integer :: qq, ich
    real(dp) :: coeff

    ich = ll * this%kmax + kk

    dgval = 0.0_dp
    do qq = 1, this%Q
       coeff = this%A(mu, qq, ich) * this%A(nu, qq, ich)
       dgval = dgval + coeff * this%u_spline(qq, ich)%derivative(rr)
    end do

  end function chem_compress_eval_drad


  !> Deallocate all arrays
  subroutine chem_compress_destroy(this)
    class(chem_compress_t), intent(inout) :: this

    if (allocated(this%A))        deallocate(this%A)
    if (allocated(this%u_spline)) deallocate(this%u_spline)
    this%initialized = .false.

  end subroutine chem_compress_destroy

end module module_chem_compress
