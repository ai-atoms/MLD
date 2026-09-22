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
#include "../../MLD_MACROS.INC"

! ============================================================================
! TB radial infrastructure  (tight-binding inspired descriptors)
!
! Builds HSVD radial functions for each TB channel (ss, pp, sp, dd).
! Each channel gets its own CubicSpline arrays built via SVD of a
! potential matrix, following the same HSVD recipe as ACE desc 300.
!
! f radial functions = HSVD radial splines indexed by (kk, ll, mua, muj)
! g radial function  = Annex smooth cutoff window (tb_g_type=1)
!
! Self-contained: does NOT import private types from module_ace_radial.
! ============================================================================

module module_tb_radial
  use iso_fortran_env, only: dp => real64
  use module_kind_variables, only: kind_double
  use module_spline_interpolation_mine, only: CubicSpline
  use mld_logger
  implicit none

  private

  ! ---------------------------------------------------------------------------
  ! Annex g-function: smooth cutoff window 
  ! g(r) = S((r - r_in) / w_in) * [1 - S((r - r_out + w_out) / w_out)]
  ! where S(t) = 6t^5 - 15t^4 + 10t^3  for 0<t<1, 0 for t<=0, 1 for t>=1.
  ! ---------------------------------------------------------------------------
  type, public :: tb_g_annex_t
    real(dp) :: r_in, w_in, r_out, w_out
    type(CubicSpline) :: spl 
    logical :: initialised = .false.
  contains
    procedure :: init   => g_annex_init
    procedure :: build  => g_annex_build
    procedure :: eval   => g_annex_eval
    procedure :: deval  => g_annex_deval
  end type tb_g_annex_t

  ! ---------------------------------------------------------------------------
  ! TB radial container for one channel
  ! Stores n_f CubicSplines for f, and one g-function object.
  ! ---------------------------------------------------------------------------
  type, public :: tb_channel_radial_t
    integer :: ch_id               ! channel id (1=ss, 2=pp, 3=sp, 4=dd)
    integer :: kmax, nmax, lmax
    ! Central (a) cutoff family: range of the envelope g(r_ai)
    real(dp) :: r_cut_in, r_cut_out, r_cut_width_in, r_cut_width_out
    ! Internal (ij) cutoff family: range of the hopping splines f(r_ij)
    ! (r_cut^(ij) of docs/Tbind_desc/main.tex; defaults inherit the (a) family)
    real(dp) :: r_cut_in_ij, r_cut_out_ij, r_cut_width_in_ij, r_cut_width_out_ij
    real(dp) :: lambda
    integer :: npoints             
    integer :: n_f                 ! total number of f radial splines
    integer :: dim_mu              ! number of species
    integer :: g_type = 1          ! 1=Annex window, 2=g=f (same HSVD function)

    type(CubicSpline), dimension(:), allocatable :: f_spline

    type(tb_g_annex_t) :: g_func

    ! ic_f(mua, muj, kk, ll) -> spline index
    integer, dimension(:,:,:,:), allocatable :: ic_f

    logical :: initialised = .false.
  contains
    procedure :: init      => channel_init
    procedure :: build_f   => channel_build_f_radials
    procedure :: build_g   => channel_build_g_radials
    procedure :: f_eval    => channel_f_eval
    procedure :: f_deval   => channel_f_deval
    procedure :: g_eval    => channel_g_eval
    procedure :: g_deval   => channel_g_deval
    procedure :: g_eval_kl => channel_g_eval_kl
    procedure :: g_deval_kl => channel_g_deval_kl
    procedure :: destroy   => channel_destroy
  end type tb_channel_radial_t

  ! Module-level radial containers per channel
  type(tb_channel_radial_t), public :: tb_rad_ss, tb_rad_pp, tb_rad_sp, tb_rad_dd

contains

  ! ===== Smooth polynomial S(t) =============================================
  pure function smooth_S(t) result(s)
    real(dp), intent(in) :: t
    real(dp) :: s
    if (t <= 0.0_dp) then
      s = 0.0_dp
    else if (t >= 1.0_dp) then
      s = 1.0_dp
    else
      s = 6.0_dp*t**5 - 15.0_dp*t**4 + 10.0_dp*t**3
    end if
  end function smooth_S

  pure function smooth_dS(t) result(ds)
    real(dp), intent(in) :: t
    real(dp) :: ds
    if (t <= 0.0_dp .or. t >= 1.0_dp) then
      ds = 0.0_dp
    else
      ds = 30.0_dp*t**4 - 60.0_dp*t**3 + 30.0_dp*t**2
    end if
  end function smooth_dS

  ! ===== HSVD potential function (local copy, same as ACE) ===================
  ! V_nl(r) = Z_j*Z_a/1000 * (r/r0)^l * exp(-alpha*(r-r0)^2)
  subroutine tb_V_nl_and_derivative(r_ja, Z_j, Z_a, r_Z_j, r_Z_a, n, l, &
                                     alpha0, beta, radius, V_nl, dV_nl_dr)
    implicit none
    integer, intent(in) :: n, l
    real(dp), intent(in) :: r_ja, alpha0, beta, r_Z_j, r_Z_a, Z_j, Z_a, radius
    real(dp), intent(out) :: V_nl, dV_nl_dr
    real(dp) :: alpha, prefactor, exp_term, term1, term2, r0, r_safe, r0_safe

    r0 = radius * (r_Z_a + r_Z_j)
    r0_safe = max(r0, 1.0d-12)
    r_safe  = max(r_ja, 1.0d-12)

    alpha = alpha0 * (Z_j**beta + Z_a**beta) / dble(n + l + 1)
    prefactor = (Z_j * Z_a) * (r_safe / r0_safe)**l / 1000.0_dp

    term1 = dble(l) / r_safe
    exp_term = exp(-alpha * (r_safe - r0_safe) ** 2)
    term2 = -2.0_dp * alpha * (r_safe - r0_safe)

    V_nl = prefactor * exp_term
    dV_nl_dr = prefactor * exp_term * (term1 + term2)
  end subroutine tb_V_nl_and_derivative

  ! ===== Annex g function ====================================================
  subroutine g_annex_init(this, r_in, w_in, r_out, w_out)
    class(tb_g_annex_t), intent(inout) :: this
    real(dp), intent(in) :: r_in, w_in, r_out, w_out
    this%r_in  = r_in
    this%w_in  = w_in
    this%r_out = r_out
    this%w_out = w_out
    this%initialised = .true.
  end subroutine g_annex_init

  subroutine g_annex_build(this, npoints, r_grid_in, r_grid_out)
    class(tb_g_annex_t), intent(inout) :: this
    integer, intent(in) :: npoints
    real(dp), intent(in) :: r_grid_in, r_grid_out
    real(dp), dimension(:), allocatable :: xgrid, yfunc, d_yfunc
    real(dp) :: rr, t_in, t_out, s_in, s_out, ds_in, ds_out
    integer :: ip

    call this%spl%init(r_grid_in, r_grid_out, npoints, xgrid, yfunc, d_yfunc)
    do ip = 1, npoints
      rr = xgrid(ip)
      ! Inner cutoff ramp-up
      if (this%w_in > 1.0d-14) then
        t_in = (rr - this%r_in) / this%w_in
        s_in = smooth_S(t_in)
        ds_in = smooth_dS(t_in) / this%w_in
      else
        s_in = 1.0_dp
        ds_in = 0.0_dp
      end if
      ! Outer cutoff ramp-down
      if (this%w_out > 1.0d-14) then
        t_out = (rr - this%r_out + this%w_out) / this%w_out
        s_out = smooth_S(t_out)
        ds_out = smooth_dS(t_out) / this%w_out
      else
        s_out = 0.0_dp
        ds_out = 0.0_dp
      end if
      yfunc(ip)   = s_in * (1.0_dp - s_out)
      d_yfunc(ip) = ds_in * (1.0_dp - s_out) - s_in * ds_out
    end do
    call this%spl%compute(yfunc, d_yfunc)
    if (allocated(xgrid))   deallocate(xgrid)
    if (allocated(yfunc))   deallocate(yfunc)
    if (allocated(d_yfunc)) deallocate(d_yfunc)
  end subroutine g_annex_build

  function g_annex_eval(this, rr) result(val)
    class(tb_g_annex_t), intent(in) :: this
    real(dp), intent(in) :: rr
    real(dp) :: val
    val = this%spl%evaluate(rr)
  end function g_annex_eval

  function g_annex_deval(this, rr) result(dval)
    class(tb_g_annex_t), intent(in) :: this
    real(dp), intent(in) :: rr
    real(dp) :: dval
    dval = this%spl%derivative(rr)
  end function g_annex_deval

  ! ===== Channel radial init =================================================
  subroutine channel_init(this, ch_id, kmax, nmax, lmax, dim_mu, &
                           r_cut_in, r_cut_out, r_cut_width_in, r_cut_width_out, &
                           lambda, npoints, g_type, &
                           r_cut_in_ij, r_cut_out_ij, r_cut_width_in_ij, r_cut_width_out_ij)
    class(tb_channel_radial_t), intent(inout) :: this
    integer, intent(in) :: ch_id, kmax, nmax, lmax, dim_mu, npoints
    integer, intent(in), optional :: g_type
    real(dp), intent(in) :: r_cut_in, r_cut_out, r_cut_width_in, r_cut_width_out, lambda
    ! Optional internal (ij) cutoff family for the hopping splines f(r_ij);
    ! when absent the central (a) family is inherited (legacy behaviour).
    real(dp), intent(in), optional :: r_cut_in_ij, r_cut_out_ij
    real(dp), intent(in), optional :: r_cut_width_in_ij, r_cut_width_out_ij
    integer :: mua, muj, kk, ll, icount

    this%ch_id    = ch_id
    this%g_type   = 1
    if (present(g_type)) this%g_type = g_type
    this%kmax     = kmax
    this%nmax     = nmax
    this%lmax     = lmax
    this%dim_mu   = dim_mu
    this%r_cut_in = r_cut_in
    this%r_cut_out = r_cut_out
    this%r_cut_width_in = r_cut_width_in
    this%r_cut_width_out = r_cut_width_out
    ! Internal (ij) family: inherit the (a) family unless explicitly provided
    this%r_cut_in_ij        = r_cut_in
    this%r_cut_out_ij       = r_cut_out
    this%r_cut_width_in_ij  = r_cut_width_in
    this%r_cut_width_out_ij = r_cut_width_out
    if (present(r_cut_in_ij))        this%r_cut_in_ij        = r_cut_in_ij
    if (present(r_cut_out_ij))       this%r_cut_out_ij       = r_cut_out_ij
    if (present(r_cut_width_in_ij))  this%r_cut_width_in_ij  = r_cut_width_in_ij
    if (present(r_cut_width_out_ij)) this%r_cut_width_out_ij = r_cut_width_out_ij
    this%lambda   = lambda
    this%npoints  = npoints

    ! Build ic_f index mapping for canonical unordered species pairs (mua <= muj):
    ! f^{AB} is stored once and aliased to f^{BA} so the radial is symmetric
    ! under neighbour permutation (i.e. f(mua,muj,...) == f(muj,mua,...)).
    if (allocated(this%ic_f)) deallocate(this%ic_f)
    allocate(this%ic_f(dim_mu, dim_mu, 1:kmax, 0:lmax))
    icount = 0
    do mua = 1, dim_mu
      do muj = mua, dim_mu
        do ll = 0, lmax
          do kk = 1, kmax
            icount = icount + 1
            this%ic_f(mua, muj, kk, ll) = icount
            this%ic_f(muj, mua, kk, ll) = icount   ! mirror swapped pair
          end do
        end do
      end do
    end do
    this%n_f = icount

    ! Allocate f splines
    if (allocated(this%f_spline)) deallocate(this%f_spline)
    allocate(this%f_spline(this%n_f))

    this%initialised = .true.
    call log_info("TB channel "//vtoa(ch_id)//": n_f="//vtoa(this%n_f)// &
                  " kmax="//vtoa(kmax)//" nmax="//vtoa(nmax)//" lmax="//vtoa(lmax)// &
                  " dim_mu="//vtoa(dim_mu))
  end subroutine channel_init

  ! ===== Build f radial functions using HSVD machinery =======================
  subroutine channel_build_f_radials(this, type_f_radial)
    use RadialFunctions, only: RadialFunction
    use module_svd_small_gen_matrix, only: eigen_svd
    use mld_mpi, only: mld_rank, mld_size, comm_mld
    use module_chemical_species, only: fix_type_to_periodic, periodic_table_element
    implicit none

    class(tb_channel_radial_t), intent(inout) :: this
    integer, intent(in) :: type_f_radial  ! polynomial radial type (1=powTcheb, etc.)

    type(RadialFunction) :: PolyRad_local
    type(eigen_svd) :: svd_local
    real(dp), dimension(:), allocatable :: xgrid, yfunc, d_yfunc
    real(dp), dimension(:,:), allocatable :: mat_vv, d_mat_vv
    real(dp), dimension(:), allocatable :: snorm
    real(dp), dimension(:,:), allocatable :: vk_local
    real(dp), dimension(:), allocatable :: sigmak_local
    real(dp), dimension(:,:), allocatable :: delta_local
    real(dp) :: rr, ffnn, dffnn, dtmp1, dtmp2
    real(dp) :: alpha, beta, radius, Z_j, Z_a, r_Z_j, r_Z_a, V_nl, dV_nl_dr
    integer :: ii, ip, mua, muj, kk, ll, nn, mm, dim_vv, ivect
    integer :: itask, ispline, e_nn, e_params
    ! Alpha-beta parameters
    integer :: no_alpha, no_beta
    integer :: ialpha, ibeta, no_of_parameters
    real(dp), dimension(:), allocatable :: alpha_arr, beta_arr

    _NAMECURRENT_("channel_build_f_radials")
    _MLD_BEGIN_

    ! Initialise polynomial radial basis on the internal (ij) range: the f
    ! splines are hopping functions of r_ij (equal to the (a) range by default)
    call PolyRad_local%init(type_f_radial, this%r_cut_in_ij, this%r_cut_out_ij, &
                            this%r_cut_width_in_ij, this%r_cut_width_out_ij, this%lambda, this%nmax)

    ! Build alpha-beta parameter grid (same as ACE init_alpha_beta_params_lin)
    no_alpha = 5
    no_beta  = 5
    no_of_parameters = no_alpha * no_beta
    allocate(alpha_arr(no_alpha), beta_arr(no_beta))
    do ialpha = 1, no_alpha
      if (no_alpha > 1) then
        alpha_arr(ialpha) = 0.05_dp + (0.2_dp - 0.05_dp) * dble(ialpha-1) / dble(no_alpha - 1)
      else
        alpha_arr(ialpha) = 0.05_dp
      end if
    end do
    do ibeta = 1, no_beta
      if (no_beta > 1) then
        beta_arr(ibeta) = 0.05_dp + (0.4_dp - 0.05_dp) * dble(ibeta-1) / dble(no_beta - 1)
      else
        beta_arr(ibeta) = 0.05_dp
      end if
    end do

    dim_vv = this%dim_mu * this%nmax * no_of_parameters
    allocate(mat_vv(this%npoints, dim_vv), d_mat_vv(this%npoints, dim_vv))
    allocate(vk_local(this%kmax, dim_vv), sigmak_local(this%kmax))
    allocate(snorm(dim_vv))
    allocate(delta_local(this%dim_mu, this%dim_mu))
    delta_local(:,:) = 0.0_dp
    do mm = 1, this%dim_mu
      delta_local(mm, mm) = 1.0_dp
    end do

    ! Initialise all f splines with zero coefficients (on the (ij) range)
    do ii = 1, this%n_f
      call this%f_spline(ii)%init(this%r_cut_in_ij, this%r_cut_out_ij, this%npoints, xgrid, yfunc, d_yfunc)
      yfunc(:) = 0.0_dp
      d_yfunc(:) = 0.0_dp
      call this%f_spline(ii)%compute(yfunc, d_yfunc)
    end do

    call log_info("TB ch"//vtoa(this%ch_id)//": HSVD mat "//vtoa(this%npoints)//"x"//vtoa(dim_vv)// &
                  " on "//vtoa(mld_size)//" ranks")

    ! MPI-distributed SVD loop over canonical unordered pairs (mua <= muj) and ll
    itask = -1
    do mua = 1, this%dim_mu
    do muj = mua, this%dim_mu
    do ll = 0, this%lmax
      itask = itask + 1
      if (mod(itask, mld_size) /= mld_rank) cycle

      Z_j = periodic_table_element(fix_type_to_periodic(muj))%Z
      Z_a = periodic_table_element(fix_type_to_periodic(mua))%Z
      r_Z_j = periodic_table_element(fix_type_to_periodic(muj))%covalent_radius / 100.0_dp
      r_Z_a = periodic_table_element(fix_type_to_periodic(mua))%covalent_radius / 100.0_dp

      ! Build the potential matrix for this (mua, muj, ll)
      do ip = 1, this%npoints
        rr = xgrid(ip)
        call PolyRad_local%evaluate(rr)
        ivect = 0
        do e_nn = 1, this%nmax
          ffnn  = PolyRad_local%radial(e_nn)
          dffnn = PolyRad_local%d_radial(e_nn)
          do ialpha = 1, no_alpha
          do ibeta = 1, no_beta
            alpha = alpha_arr(ialpha)
            beta  = beta_arr(ibeta)
            radius = 1.0_dp
            call tb_V_nl_and_derivative(rr, Z_j, Z_a, r_Z_j, r_Z_a, e_nn, ll, &
                                        alpha, beta, radius, V_nl, dV_nl_dr)
            do nn = 1, this%dim_mu
              ivect = ivect + 1
              mat_vv(ip, ivect)   = delta_local(nn, muj) * V_nl * ffnn
              d_mat_vv(ip, ivect) = delta_local(nn, muj) * (dV_nl_dr * ffnn + V_nl * dffnn)
            end do
          end do
          end do
        end do
      end do

      ! Normalise columns
      do mm = 1, dim_vv
        snorm(mm) = sqrt(sum(mat_vv(:,mm)**2))
        if (snorm(mm) <= 1.0d-14) snorm(mm) = 1.0_dp
        mat_vv(:,mm)   = mat_vv(:,mm)   / snorm(mm)
        d_mat_vv(:,mm) = d_mat_vv(:,mm) / snorm(mm)
      end do

      ! SVD
      call svd_local%init(size(mat_vv, 1), size(mat_vv, 2))
      call svd_local%evaluate(mat_vv)
      vk_local(1:this%kmax, 1:dim_vv) = svd_local%eigvecVT(1:this%kmax, 1:dim_vv)
      sigmak_local(1:this%kmax) = svd_local%eigval(1:this%kmax)

      ! Build splines for each kk (on the (ij) range)
      do kk = 1, this%kmax
        ispline = this%ic_f(mua, muj, kk, ll)
        call this%f_spline(ispline)%init(this%r_cut_in_ij, this%r_cut_out_ij, this%npoints, xgrid, yfunc, d_yfunc)
        do ip = 1, this%npoints
          dtmp1 = 0.0_dp
          dtmp2 = 0.0_dp
          do mm = 1, dim_vv
            dtmp1 = dtmp1 + vk_local(kk, mm) * mat_vv(ip, mm)
            dtmp2 = dtmp2 + vk_local(kk, mm) * d_mat_vv(ip, mm)
          end do
          yfunc(ip)   = dtmp1 / sigmak_local(kk)
          d_yfunc(ip) = dtmp2 / sigmak_local(kk)
        end do
        call this%f_spline(ispline)%compute(yfunc, d_yfunc)
      end do

    end do ! ll
    end do ! muj
    end do ! mua

    ! Allreduce all spline coefficients across MPI ranks
    do ii = 1, this%n_f
      call comm_mld%sum(this%f_spline(ii)%a)
      call comm_mld%sum(this%f_spline(ii)%b)
      call comm_mld%sum(this%f_spline(ii)%c)
      call comm_mld%sum(this%f_spline(ii)%d)
    end do

    deallocate(mat_vv, d_mat_vv, vk_local, sigmak_local, snorm, delta_local, alpha_arr, beta_arr)
    if (allocated(xgrid))   deallocate(xgrid)
    if (allocated(yfunc))   deallocate(yfunc)
    if (allocated(d_yfunc)) deallocate(d_yfunc)

    call log_info("TB ch"//vtoa(this%ch_id)//": f radials built ("//vtoa(this%n_f)//" splines)")
    _MLD_END_
  end subroutine channel_build_f_radials

  ! ===== Build g radial functions (Annex smooth window) ======================
  subroutine channel_build_g_radials(this)
    class(tb_channel_radial_t), intent(inout) :: this
    call this%g_func%init(this%r_cut_in, this%r_cut_width_in, this%r_cut_out, this%r_cut_width_out)
    call this%g_func%build(this%npoints, this%r_cut_in, this%r_cut_out)
    call log_info("TB ch"//vtoa(this%ch_id)//": g Annex window built")
  end subroutine channel_build_g_radials

  ! ===== Evaluation helpers ==================================================
  function channel_f_eval(this, mua, muj, kk, ll, rr) result(val)
    class(tb_channel_radial_t), intent(in) :: this
    integer, intent(in) :: mua, muj, kk, ll
    real(dp), intent(in) :: rr
    real(dp) :: val
    integer :: isp
    ! f(r_ij) lives on the internal (ij) range: h_{ij}=0 beyond r_cut_out_ij
    if (rr < this%r_cut_in_ij .or. rr > this%r_cut_out_ij) then
      val = 0.0_dp
      return
    end if
    isp = this%ic_f(mua, muj, kk, ll)
    val = this%f_spline(isp)%evaluate(rr)
  end function channel_f_eval

  function channel_f_deval(this, mua, muj, kk, ll, rr) result(dval)
    class(tb_channel_radial_t), intent(in) :: this
    integer, intent(in) :: mua, muj, kk, ll
    real(dp), intent(in) :: rr
    real(dp) :: dval
    integer :: isp
    if (rr < this%r_cut_in_ij .or. rr > this%r_cut_out_ij) then
      dval = 0.0_dp
      return
    end if
    isp = this%ic_f(mua, muj, kk, ll)
    dval = this%f_spline(isp)%derivative(rr)
  end function channel_f_deval

  function channel_g_eval(this, rr) result(val)
    class(tb_channel_radial_t), intent(in) :: this
    real(dp), intent(in) :: rr
    real(dp) :: val
    if (rr < this%r_cut_in .or. rr > this%r_cut_out) then
      val = 0.0_dp
      return
    end if
    val = this%g_func%eval(rr)
  end function channel_g_eval

  function channel_g_deval(this, rr) result(dval)
    class(tb_channel_radial_t), intent(in) :: this
    real(dp), intent(in) :: rr
    real(dp) :: dval
    if (rr < this%r_cut_in .or. rr > this%r_cut_out) then
      dval = 0.0_dp
      return
    end if
    dval = this%g_func%deval(rr)
  end function channel_g_deval

  ! ===== (kk,ll)-aware g evaluation (dispatches on g_type) ==================
  ! tb_g_type=1: g = Annex window (no dependence on kk,ll)
  ! tb_g_type=2: g = f (same HSVD spline as the bond integral)
  function channel_g_eval_kl(this, mua, muj, kk, ll, rr) result(val)
    class(tb_channel_radial_t), intent(in) :: this
    integer, intent(in) :: mua, muj, kk, ll
    real(dp), intent(in) :: rr
    real(dp) :: val
    if (this%g_type == 2) then
      val = this%f_eval(mua, muj, kk, ll, rr)
    else
      val = this%g_eval(rr)
    end if
  end function channel_g_eval_kl

  function channel_g_deval_kl(this, mua, muj, kk, ll, rr) result(dval)
    class(tb_channel_radial_t), intent(in) :: this
    integer, intent(in) :: mua, muj, kk, ll
    real(dp), intent(in) :: rr
    real(dp) :: dval
    if (this%g_type == 2) then
      dval = this%f_deval(mua, muj, kk, ll, rr)
    else
      dval = this%g_deval(rr)
    end if
  end function channel_g_deval_kl

  ! ===== Cleanup =============================================================
  subroutine channel_destroy(this)
    class(tb_channel_radial_t), intent(inout) :: this
    if (allocated(this%f_spline)) deallocate(this%f_spline)
    if (allocated(this%ic_f)) deallocate(this%ic_f)
    this%initialised = .false.
  end subroutine channel_destroy

end module module_tb_radial


! ============================================================================
! Spectral filter library for the equiv-B (Model A) descriptors
! (docs/README_tbind_filter_implementation.md, docs/Tbind_desc/main.tex App.)
!
! Families: polynomial (legacy raw lambda^q, handled OUTSIDE this module for
! strict backward compatibility), chebyshev, gaussian, lorentzian, resolvent,
! fermi, fermi_bin, bspline.
!
! All non-polynomial families act on the normalized eigenvalue
!     x = (lambda - lam0_c) / dlam_c
! with per-channel constants lam0_c, dlam_c built from the user-provided
! reference interval [lmin_c, lmax_c] widened by the spectral padding:
!     lam0_c = (lmax+lmin)/2,   dlam_c = (1+pad)*(lmax-lmin)/2 .
! tbf_eval returns both f_k(lambda) and df_k/dlambda (the 1/dlam_c chain
! factor included), which feeds the common divided-difference force kernel.
! ============================================================================
module module_tb_filter
  use iso_fortran_env, only: dp => real64
  use mld_logger
  implicit none

  private

  integer, parameter, public :: TBF_POLY  = 0, TBF_CHEB = 1, TBF_GAUSS = 2, &
                                 TBF_LOR   = 3, TBF_RESO = 4, TBF_FERMI = 5, &
                                 TBF_FBIN  = 6, TBF_BSPL = 7

  integer, public :: tbf_type = TBF_POLY   ! resolved family
  integer, public :: tbf_nf   = 0          ! scalar filters per channel (0 for poly)

  integer,  public :: tbf_order   = 10     ! chebyshev Q (T_0..T_Q)
  integer  :: tbf_ncenter = 8              ! K centers / thresholds / basis functions
  integer  :: tbf_degree  = 3              ! bspline degree (only 3 supported)
  real(dp), public :: tbf_wpar = 0.25_dp   ! sigma/eta/tau in x units
  logical,  public :: tbf_first_moment = .true.  ! gaussian: add x*G_k
  logical,  public :: tbf_re = .true., tbf_im = .true.  ! resolvent components
  real(dp), public :: tbf_lam0(4) = 0.0_dp, tbf_dlam(4) = 1.0_dp
  real(dp), allocatable, public :: tbf_centers(:)  ! centers/thresholds in x units
  real(dp), allocatable :: tbf_knots(:)    ! bspline clamped-uniform knots

  public :: tbf_init, tbf_eval, tbf_type_name

contains

  ! ===== helpers =============================================================
  subroutine tbf_lower(s)
    character(len=*), intent(inout) :: s
    integer :: i, ic
    do i = 1, len_trim(s)
      ic = iachar(s(i:i))
      if (ic >= iachar('A') .and. ic <= iachar('Z')) s(i:i) = achar(ic + 32)
    end do
  end subroutine tbf_lower

  function tbf_type_name(it) result(nm)
    integer, intent(in) :: it
    character(len=12) :: nm
    select case (it)
    case (TBF_POLY);  nm = 'polynomial'
    case (TBF_CHEB);  nm = 'chebyshev'
    case (TBF_GAUSS); nm = 'gaussian'
    case (TBF_LOR);   nm = 'lorentzian'
    case (TBF_RESO);  nm = 'resolvent'
    case (TBF_FERMI); nm = 'fermi'
    case (TBF_FBIN);  nm = 'fermi_bin'
    case (TBF_BSPL);  nm = 'bspline'
    case default;     nm = 'unknown'
    end select
  end function tbf_type_name

  ! ===== initialisation ======================================================
  subroutine tbf_init(type_str, order, ncenter, width_abs, width_ratio, &
                      first_moment, use_re, use_im, degree, pad, &
                      lmin, lmax, ham_active)
    character(len=*), intent(in) :: type_str
    integer,  intent(in) :: order, ncenter, degree
    real(dp), intent(in) :: width_abs, width_ratio, pad
    logical,  intent(in) :: first_moment, use_re, use_im
    real(dp), dimension(4), intent(in) :: lmin, lmax
    logical,  dimension(4), intent(in) :: ham_active

    character(len=64) :: ts
    real(dp) :: h, ratio
    integer  :: K, j, nthr, ich

    ts = adjustl(type_str)
    call tbf_lower(ts)

    select case (trim(ts))
    case ('polynomial', '');  tbf_type = TBF_POLY
    case ('chebyshev');       tbf_type = TBF_CHEB
    case ('gaussian');        tbf_type = TBF_GAUSS
    case ('lorentzian');      tbf_type = TBF_LOR
    case ('resolvent');       tbf_type = TBF_RESO
    case ('fermi');           tbf_type = TBF_FERMI
    case ('fermi_bin');       tbf_type = TBF_FBIN
    case ('bspline');         tbf_type = TBF_BSPL
    case default
      call log_critical("TB filter: unknown tb_filter_type '"//trim(type_str)// &
                        "' (allowed: polynomial, chebyshev, gaussian, lorentzian,"// &
                        " resolvent, fermi, fermi_bin, bspline)")
      stop
    end select

    if (tbf_type == TBF_POLY) then
      tbf_nf = 0   ! legacy raw-lambda^q path, handled by the trace machinery
      return
    end if

    ! ---- spectral normalization (mandatory for non-polynomial filters) ----
    do ich = 1, 4
      if (.not. ham_active(ich)) then
        tbf_lam0(ich) = 0.0_dp
        tbf_dlam(ich) = 1.0_dp
        cycle
      end if
      if (lmax(ich) <= lmin(ich)) then
        call log_critical("TB filter: tb_filter_lmin_list/tb_filter_lmax_list not set "// &
                          "for active channel "//vtoa(ich)//" (required for tb_filter_type='"// &
                          trim(ts)//"'). Run once with equiv-B active and read the observed "// &
                          "ranges from tbind_lambda_range.dat, then set the keywords.")
        stop
      end if
      tbf_lam0(ich) = 0.5_dp*(lmax(ich) + lmin(ich))
      tbf_dlam(ich) = (1.0_dp + pad)*0.5_dp*(lmax(ich) - lmin(ich))
    end do

    ! ---- per-type defaults ----
    tbf_first_moment = first_moment
    tbf_re = use_re
    tbf_im = use_im
    tbf_degree = degree

    select case (tbf_type)
    case (TBF_CHEB)
      tbf_order = order
      if (tbf_order < 0) tbf_order = 10
      tbf_nf = tbf_order + 1                       ! T_0 .. T_Q
    case (TBF_GAUSS)
      K = ncenter; if (K < 0) K = 8
      ratio = width_ratio; if (ratio < 0.0_dp) ratio = 1.0_dp
      call tbf_center_grid(K, K, h)
      tbf_wpar = merge(width_abs, ratio*h, width_abs > 0.0_dp)
      tbf_nf = K * merge(2, 1, tbf_first_moment)   ! (G_k[, xG_k]) per center
    case (TBF_LOR)
      K = ncenter; if (K < 0) K = 12
      ratio = width_ratio; if (ratio < 0.0_dp) ratio = 0.75_dp
      call tbf_center_grid(K, K, h)
      tbf_wpar = merge(width_abs, ratio*h, width_abs > 0.0_dp)
      tbf_nf = K
    case (TBF_RESO)
      K = ncenter; if (K < 0) K = 8
      ratio = width_ratio; if (ratio < 0.0_dp) ratio = 0.75_dp
      call tbf_center_grid(K, K, h)
      tbf_wpar = merge(width_abs, ratio*h, width_abs > 0.0_dp)
      if (.not. (tbf_re .or. tbf_im)) then
        call log_critical("TB filter resolvent: both real and imaginary parts disabled")
        stop
      end if
      tbf_nf = K * (merge(1,0,tbf_im) + merge(1,0,tbf_re))   ! (L_k[, Q_k])
    case (TBF_FERMI)
      K = ncenter; if (K < 0) K = 16
      ratio = width_ratio; if (ratio < 0.0_dp) ratio = 0.6_dp
      call tbf_center_grid(K, K, h)
      tbf_wpar = merge(width_abs, ratio*h, width_abs > 0.0_dp)
      tbf_nf = K
    case (TBF_FBIN)
      K = ncenter; if (K < 0) K = 16
      ratio = width_ratio; if (ratio < 0.0_dp) ratio = 0.6_dp
      nthr = K + 1                                 ! K bins need K+1 thresholds
      call tbf_center_grid(nthr, nthr, h)
      tbf_wpar = merge(width_abs, ratio*h, width_abs > 0.0_dp)
      tbf_nf = K
    case (TBF_BSPL)
      K = ncenter; if (K < 0) K = 16
      if (tbf_degree /= 3) then
        call log_critical("TB filter bspline: only tb_filter_degree = 3 is supported")
        stop
      end if
      if (K < 4) then
        call log_critical("TB filter bspline: tb_filter_ncenter must be >= 4")
        stop
      end if
      ! clamped uniform knot vector on [-1,1]: t(1:4)=-1, t(K+1:K+4)=+1,
      ! K-4 uniform interior knots
      if (allocated(tbf_knots)) deallocate(tbf_knots)
      allocate(tbf_knots(K+4))
      tbf_knots(1:4)     = -1.0_dp
      tbf_knots(K+1:K+4) =  1.0_dp
      do j = 1, K-4
        tbf_knots(4+j) = -1.0_dp + 2.0_dp*dble(j)/dble(K-3)
      end do
      tbf_nf = K
    end select

    tbf_ncenter = 0
    if (allocated(tbf_centers)) tbf_ncenter = size(tbf_centers)

    call log_info("TB filter: type="//trim(tbf_type_name(tbf_type))// &
                  "  n_filters/channel="//vtoa(tbf_nf)// &
                  "  width(x units)="//vtoa(tbf_wpar))
    do ich = 1, 4
      if (ham_active(ich)) then
        call log_info("TB filter ch"//vtoa(ich)//": lambda0="//vtoa(tbf_lam0(ich))// &
                      "  dlam(padded)="//vtoa(tbf_dlam(ich)))
      end if
    end do

  contains

    subroutine tbf_center_grid(npts, nstore, hout)
      integer,  intent(in)  :: npts, nstore
      real(dp), intent(out) :: hout
      integer :: i
      if (allocated(tbf_centers)) deallocate(tbf_centers)
      allocate(tbf_centers(nstore))
      if (npts == 1) then
        tbf_centers(1) = 0.0_dp
        hout = 2.0_dp
      else
        do i = 1, npts
          tbf_centers(i) = -1.0_dp + 2.0_dp*dble(i-1)/dble(npts-1)
        end do
        hout = 2.0_dp/dble(npts-1)
      end if
    end subroutine tbf_center_grid

  end subroutine tbf_init

  ! ===== evaluation: fv(k) = f_k(lambda), fp(k) = df_k/dlambda ===============
  subroutine tbf_eval(ich, lam, fv, fp)
    integer,  intent(in)  :: ich
    real(dp), intent(in)  :: lam
    real(dp), dimension(:), intent(out) :: fv, fp

    real(dp) :: x, dxdl, t, tm, tp2, u, um, up2, e, g, d2, w, ff, arg
    integer  :: k, q, iq, K0

    x    = (lam - tbf_lam0(ich)) / tbf_dlam(ich)
    dxdl = 1.0_dp / tbf_dlam(ich)

    select case (tbf_type)

    case (TBF_CHEB)
      ! T_q by recurrence; dT_q/dx = q U_{q-1}
      tm = 1.0_dp                 ! T_0
      t  = x                      ! T_1
      um = 1.0_dp                 ! U_0
      u  = 2.0_dp*x               ! U_1
      fv(1) = 1.0_dp;  fp(1) = 0.0_dp
      if (tbf_order >= 1) then
        fv(2) = x;  fp(2) = 1.0_dp*um*dxdl
      end if
      do q = 2, tbf_order
        tp2 = 2.0_dp*x*t - tm
        fv(q+1) = tp2
        fp(q+1) = dble(q)*u*dxdl        ! q * U_{q-1}
        tm = t;  t = tp2
        up2 = 2.0_dp*x*u - um
        um = u;  u = up2
      end do

    case (TBF_GAUSS)
      iq = 0
      do k = 1, size(tbf_centers)
        e = tbf_centers(k)
        g = exp(-0.5_dp*((x - e)/tbf_wpar)**2)
        iq = iq + 1
        fv(iq) = g
        fp(iq) = -(x - e)/(tbf_wpar**2)*g*dxdl
        if (tbf_first_moment) then
          iq = iq + 1
          fv(iq) = x*g
          fp(iq) = g*(1.0_dp - x*(x - e)/tbf_wpar**2)*dxdl
        end if
      end do

    case (TBF_LOR)
      w = tbf_wpar
      do k = 1, size(tbf_centers)
        e  = tbf_centers(k)
        d2 = (x - e)**2 + w*w
        fv(k) = w/d2
        fp(k) = -2.0_dp*w*(x - e)/(d2*d2)*dxdl
      end do

    case (TBF_RESO)
      w = tbf_wpar
      iq = 0
      do k = 1, size(tbf_centers)
        e  = tbf_centers(k)
        d2 = (x - e)**2 + w*w
        if (tbf_im) then
          iq = iq + 1
          fv(iq) = w/d2
          fp(iq) = -2.0_dp*w*(x - e)/(d2*d2)*dxdl
        end if
        if (tbf_re) then
          iq = iq + 1
          fv(iq) = (e - x)/d2
          fp(iq) = ((e - x)**2 - w*w)/(d2*d2)*dxdl
        end if
      end do

    case (TBF_FERMI)
      w = tbf_wpar
      do k = 1, size(tbf_centers)
        arg = (x - tbf_centers(k))/w
        ff = tbf_sigmoid(arg)
        fv(k) = ff
        fp(k) = ff*(1.0_dp - ff)/w*dxdl
      end do

    case (TBF_FBIN)
      ! W_k = F(E_k) - F(E_{k+1}) over K+1 thresholds
      w = tbf_wpar
      K0 = size(tbf_centers) - 1
      t = tbf_sigmoid((x - tbf_centers(1))/w)          ! F at first threshold
      tm = t*(1.0_dp - t)/w                            ! dF/dx at first threshold
      do k = 1, K0
        u  = tbf_sigmoid((x - tbf_centers(k+1))/w)
        um = u*(1.0_dp - u)/w
        fv(k) = t - u
        fp(k) = (tm - um)*dxdl
        t = u;  tm = um
      end do

    case (TBF_BSPL)
      call tbf_bspline_all(x, fv, fp)
      fp(1:tbf_nf) = fp(1:tbf_nf)*dxdl

    case default
      fv = 0.0_dp
      fp = 0.0_dp
    end select
  end subroutine tbf_eval

  pure function tbf_sigmoid(t) result(s)
    real(dp), intent(in) :: t
    real(dp) :: s
    if (t > 100.0_dp) then
      s = 1.0_dp
    else if (t < -100.0_dp) then
      s = 0.0_dp
    else
      s = 1.0_dp/(1.0_dp + exp(-t))
    end if
  end function tbf_sigmoid

  ! ---- cubic B-splines: values and x-derivatives of all K basis functions --
  subroutine tbf_bspline_all(x_in, fv, fp)
    real(dp), intent(in) :: x_in
    real(dp), dimension(:), intent(out) :: fv, fp
    real(dp) :: x, b2k, b2k1, denom1, denom2
    integer :: k, K0

    K0 = size(tbf_knots) - 4     ! number of cubic basis functions
    fv(1:K0) = 0.0_dp
    fp(1:K0) = 0.0_dp
    if (x_in < tbf_knots(1) .or. x_in > tbf_knots(K0+4)) return
    ! close the last interval (x == +1 belongs to the last basis function)
    x = min(x_in, tbf_knots(K0+4) - 1.0e-12_dp)

    do k = 1, K0
      fv(k) = tbf_bspl_B(k, 3, x)
      ! dB_{k,3}/dx = 3 [ B_{k,2}/(t_{k+3}-t_k) - B_{k+1,2}/(t_{k+4}-t_{k+1}) ]
      b2k  = tbf_bspl_B(k,   2, x)
      b2k1 = tbf_bspl_B(k+1, 2, x)
      denom1 = tbf_knots(k+3) - tbf_knots(k)
      denom2 = tbf_knots(k+4) - tbf_knots(k+1)
      if (denom1 > 1.0e-14_dp) fp(k) = fp(k) + 3.0_dp*b2k/denom1
      if (denom2 > 1.0e-14_dp) fp(k) = fp(k) - 3.0_dp*b2k1/denom2
    end do
  end subroutine tbf_bspline_all

  pure recursive function tbf_bspl_B(k, p, x) result(v)
    integer,  intent(in) :: k, p
    real(dp), intent(in) :: x
    real(dp) :: v, d1, d2
    if (p == 0) then
      if (tbf_knots(k) <= x .and. x < tbf_knots(k+1)) then
        v = 1.0_dp
      else
        v = 0.0_dp
      end if
      return
    end if
    v = 0.0_dp
    d1 = tbf_knots(k+p) - tbf_knots(k)
    if (d1 > 1.0e-14_dp) v = v + (x - tbf_knots(k))/d1*tbf_bspl_B(k, p-1, x)
    d2 = tbf_knots(k+p+1) - tbf_knots(k+1)
    if (d2 > 1.0e-14_dp) v = v + (tbf_knots(k+p+1) - x)/d2*tbf_bspl_B(k+1, p-1, x)
  end function tbf_bspl_B

end module module_tb_filter


! ============================================================================
! tbind potential I/O: export the descriptor to tbind.xml
!
! tbind.xml is the self-contained descriptor file consumed by LAMMPS.  It
! carries every parameter needed to rebuild the descriptor plus the HSVD
! radial splines themselves.  The splines are exported (not rebuilt on the
! LAMMPS side) because they come from an SVD whose singular vectors are only
! defined up to a sign: recomputing them elsewhere could silently flip the
! sign of an f channel and change the Hamiltonian.
! ============================================================================
module module_tbind_potio
  use iso_fortran_env, only: dp => real64
  use module_tb_radial, only: tb_channel_radial_t, tb_rad_ss, tb_rad_pp, tb_rad_sp, tb_rad_dd
  use mld_logger
  implicit none

  private
  public :: write_tbind_xml

contains

  subroutine write_int_list(unit, tag, arr)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: tag
    integer, dimension(:), intent(in) :: arr
    integer :: i
    write(unit, '(A)', advance='no') '    <'//trim(tag)//'>'
    do i = 1, size(arr)
      write(unit, '(I0,1X)', advance='no') arr(i)
    end do
    write(unit, '(A)') '</'//trim(tag)//'>'
  end subroutine write_int_list

  subroutine write_real_list(unit, tag, arr)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: tag
    real(dp), dimension(:), intent(in) :: arr
    integer :: i
    write(unit, '(A)', advance='no') '    <'//trim(tag)//'>'
    do i = 1, size(arr)
      write(unit, '(ES25.15E3,1X)', advance='no') arr(i)
    end do
    write(unit, '(A)') '</'//trim(tag)//'>'
  end subroutine write_real_list

  ! One cubic spline: a(n), b(n), c(n), d(n-1)
  subroutine write_spline(unit, indent, n, delta, a, b, c, d)
    integer, intent(in) :: unit, n
    character(len=*), intent(in) :: indent
    real(dp), intent(in) :: delta
    real(dp), dimension(:), intent(in) :: a, b, c, d
    integer :: j

    write(unit, '(A,I0,A)') indent//'<n>', n, '</n>'
    write(unit, '(A,ES25.15E3,A)') indent//'<delta>', delta, '</delta>'
    write(unit, '(A)', advance='no') indent//'<a>'
    do j = 1, size(a); write(unit, '(ES25.15E3,1X)', advance='no') a(j); end do
    write(unit, '(A)') '</a>'
    write(unit, '(A)', advance='no') indent//'<b>'
    do j = 1, size(b); write(unit, '(ES25.15E3,1X)', advance='no') b(j); end do
    write(unit, '(A)') '</b>'
    write(unit, '(A)', advance='no') indent//'<c>'
    do j = 1, size(c); write(unit, '(ES25.15E3,1X)', advance='no') c(j); end do
    write(unit, '(A)') '</c>'
    write(unit, '(A)', advance='no') indent//'<d>'
    do j = 1, size(d); write(unit, '(ES25.15E3,1X)', advance='no') d(j); end do
    write(unit, '(A)') '</d>'
  end subroutine write_spline

  subroutine write_channel(unit, rad, g_type, write_ij)
    integer, intent(in) :: unit, g_type
    type(tb_channel_radial_t), intent(in) :: rad
    logical, intent(in), optional :: write_ij
    integer :: isp
    logical :: do_ij

    do_ij = .false.
    if (present(write_ij)) do_ij = write_ij

    write(unit, '(A,I0,A)') '    <tb_channel ch_id="', rad%ch_id, '">'
    write(unit, '(A,I0,A)') '      <n_f>', rad%n_f, '</n_f>'
    write(unit, '(A,I0,A)') '      <kmax>', rad%kmax, '</kmax>'
    write(unit, '(A,I0,A)') '      <lmax>', rad%lmax, '</lmax>'
    write(unit, '(A,I0,A)') '      <dim_mu>', rad%dim_mu, '</dim_mu>'
    write(unit, '(A,I0,A)') '      <npoints>', rad%npoints, '</npoints>'
    write(unit, '(A,ES25.15E3,A)') '      <r_cut_in>', rad%r_cut_in, '</r_cut_in>'
    write(unit, '(A,ES25.15E3,A)') '      <r_cut_out>', rad%r_cut_out, '</r_cut_out>'
    ! Internal (ij) hopping range, emitted only when it differs from the (a)
    ! family (keeps legacy tbind.xml byte-identical for the LAMMPS reader)
    if (do_ij) then
      write(unit, '(A,ES25.15E3,A)') '      <r_cut_in_ij>', rad%r_cut_in_ij, '</r_cut_in_ij>'
      write(unit, '(A,ES25.15E3,A)') '      <r_cut_out_ij>', rad%r_cut_out_ij, '</r_cut_out_ij>'
    end if

    write(unit, '(A)') '      <f_splines>'
    do isp = 1, rad%n_f
      write(unit, '(A,I0,A)') '        <spline idx="', isp, '">'
      call write_spline(unit, '          ', rad%f_spline(isp)%n, rad%f_spline(isp)%delta, &
                        rad%f_spline(isp)%a, rad%f_spline(isp)%b, &
                        rad%f_spline(isp)%c, rad%f_spline(isp)%d)
      write(unit, '(A)') '        </spline>'
    end do
    write(unit, '(A)') '      </f_splines>'

    ! g_type=2 means g=f, so no separate window spline exists.
    if (g_type == 1) then
      write(unit, '(A)') '      <g_spline>'
      call write_spline(unit, '        ', rad%g_func%spl%n, rad%g_func%spl%delta, &
                        rad%g_func%spl%a, rad%g_func%spl%b, &
                        rad%g_func%spl%c, rad%g_func%spl%d)
      write(unit, '(A)') '      </g_spline>'
    end if

    write(unit, '(A)') '    </tb_channel>'
  end subroutine write_channel

  subroutine write_tbind_xml(fname)
    use module_tbind
    use module_tb_filter, only: tbf_type, tbf_nf, tbf_wpar, tbf_lam0, tbf_dlam, &
                                tbf_centers, tbf_order, tbf_first_moment, tbf_re, tbf_im, &
                                tbf_type_name
    use module_chemical_species, only: fix_no_of_elements
    character(len=*), intent(in) :: fname
    integer :: unit
    logical :: ij_differs

    ! The internal (ij) cutoff family is exported only when it differs from
    ! the central (a) family, and the equiv-B tag only when the model is
    ! active: legacy runs keep a byte-identical tbind.xml (the current
    ! LAMMPS reader does not yet parse the new tags).
    ij_differs = .false.
    if (allocated(tb_rcut_in_ij)) then
      ij_differs = any(abs(tb_rcut_in_ij - tb_rcut_in) > 1.0e-14_dp) .or. &
                   any(abs(tb_rcut_out_ij - tb_rcut_out) > 1.0e-14_dp) .or. &
                   any(abs(tb_rcut_width_in_ij - tb_rcut_width_in) > 1.0e-14_dp) .or. &
                   any(abs(tb_rcut_width_out_ij - tb_rcut_width_out) > 1.0e-14_dp)
    end if

    open(newunit=unit, file=fname, status='replace', action='write')

    write(unit, '(A)') '<?xml version="1.0"?>'
    write(unit, '(A)') '<tbind>'

    write(unit, '(A)') '  <tbind_params>'
    write(unit, '(A,I0,A)') '    <fix_no_of_elements>', fix_no_of_elements, '</fix_no_of_elements>'
    write(unit, '(A,I0,A)') '    <tbind_dim>', tbind_dim, '</tbind_dim>'
    write(unit, '(A,I0,A)') '    <tbind_dim_per_species>', tbind_dim_per_species, '</tbind_dim_per_species>'
    write(unit, '(A,I0,A)') '    <tb_g_type>', tb_g_type, '</tb_g_type>'
    write(unit, '(A,I0,A)') '    <tb_svd>', tb_svd, '</tb_svd>'
    write(unit, '(A,I0,A)') '    <tb_model_lambda>', merge(1, 0, tb_model_lambda), '</tb_model_lambda>'
    write(unit, '(A,I0,A)') '    <tb_model_trace>', merge(1, 0, tb_model_trace), '</tb_model_trace>'
    if (tb_model_equivb) then
      write(unit, '(A,I0,A)') '    <tb_model_equivb>', 1, '</tb_model_equivb>'
    end if
    ! Spectral-filter configuration (only for non-polynomial filters, so the
    ! legacy tbind.xml stays byte-identical; not parsed by LAMMPS yet)
    if (tb_model_equivb .and. .not. tb_equivb_poly) then
      write(unit, '(A)') '    <tb_filter_type>'//trim(tbf_type_name(tbf_type))//'</tb_filter_type>'
      write(unit, '(A,I0,A)') '    <tb_filter_nf>', tbf_nf, '</tb_filter_nf>'
      write(unit, '(A,I0,A)') '    <tb_filter_order>', tbf_order, '</tb_filter_order>'
      write(unit, '(A,ES25.15E3,A)') '    <tb_filter_wpar>', tbf_wpar, '</tb_filter_wpar>'
      write(unit, '(A,I0,A)') '    <tb_filter_first_moment>', merge(1, 0, tbf_first_moment), '</tb_filter_first_moment>'
      write(unit, '(A,I0,A)') '    <tb_filter_use_real>', merge(1, 0, tbf_re), '</tb_filter_use_real>'
      write(unit, '(A,I0,A)') '    <tb_filter_use_imag>', merge(1, 0, tbf_im), '</tb_filter_use_imag>'
      call write_real_list(unit, 'tb_filter_lam0', tbf_lam0)
      call write_real_list(unit, 'tb_filter_dlam', tbf_dlam)
      if (allocated(tbf_centers)) call write_real_list(unit, 'tb_filter_centers', tbf_centers)
    end if
    write(unit, '(A,I0,A)') '    <tb_ham_ss>', merge(1, 0, tb_ham_ss), '</tb_ham_ss>'
    write(unit, '(A,I0,A)') '    <tb_ham_pp>', merge(1, 0, tb_ham_pp), '</tb_ham_pp>'
    write(unit, '(A,I0,A)') '    <tb_ham_sp>', merge(1, 0, tb_ham_sp), '</tb_ham_sp>'
    write(unit, '(A,I0,A)') '    <tb_ham_dd>', merge(1, 0, tb_ham_dd), '</tb_ham_dd>'

    ! Effective (post-quantisation) values as used by the fit
    call write_int_list(unit, 'tb_kmax', tb_kmax)
    call write_int_list(unit, 'tb_nmax', tb_nmax)
    call write_int_list(unit, 'tb_lmax', tb_lmax)
    call write_int_list(unit, 'tb_kmax_desc', tb_kmax_desc)
    call write_int_list(unit, 'tb_lmax_desc', tb_lmax_desc)
    call write_int_list(unit, 'tb_nn_max', tb_nn_max)
    call write_int_list(unit, 'tb_power_trace', tb_power_trace)
    call write_real_list(unit, 'tb_lambda', tb_lambda_arr)
    call write_real_list(unit, 'tb_rcut_in', tb_rcut_in)
    call write_real_list(unit, 'tb_rcut_out', tb_rcut_out)
    call write_real_list(unit, 'tb_rcut_width_in', tb_rcut_width_in)
    call write_real_list(unit, 'tb_rcut_width_out', tb_rcut_width_out)
    if (ij_differs) then
      call write_real_list(unit, 'tb_rcut_in_ij', tb_rcut_in_ij)
      call write_real_list(unit, 'tb_rcut_out_ij', tb_rcut_out_ij)
      call write_real_list(unit, 'tb_rcut_width_in_ij', tb_rcut_width_in_ij)
      call write_real_list(unit, 'tb_rcut_width_out_ij', tb_rcut_width_out_ij)
    end if
    write(unit, '(A)') '  </tbind_params>'

    write(unit, '(A)') '  <tb_channels>'
    if (tb_ham_ss) call write_channel(unit, tb_rad_ss, tb_g_type, ij_differs)
    if (tb_ham_pp) call write_channel(unit, tb_rad_pp, tb_g_type, ij_differs)
    if (tb_ham_sp) call write_channel(unit, tb_rad_sp, tb_g_type, ij_differs)
    if (tb_ham_dd) call write_channel(unit, tb_rad_dd, tb_g_type, ij_differs)
    write(unit, '(A)') '  </tb_channels>'

    write(unit, '(A)') '</tbind>'
    close(unit)

    call log_info("write_tbind_xml: tbind descriptor written to "//trim(fname))
  end subroutine write_tbind_xml

end module module_tbind_potio


! ============================================================================
! s-s Hamiltonian  (Section 2 of Tbind_desc.pdf)
!
! h_{a;i,j}^{ss} = f(r_ij) * (g(r_ai) + g(r_aj))^2    for i /= j
!                   (2*g(r_ak))^2                        for i = j (= k)
!
! Descriptor = sorted eigenvalues of h, truncated to N_lambda.
! Forces via Hellmann-Feynman theorem.
! ============================================================================
module module_ss_hamiltonian
  use iso_fortran_env, only: dp => real64
  use module_eigenvalue_small_matrix, only: eigen_symmetric
  implicit none

  private

  type, public :: ss_hamiltonian
    integer :: dim_hh       ! matrix dimension = max_neigh_local
    integer :: dim_nn       ! number of neighbours
    integer :: tb_nn_G      ! eigenvalue cutoff (descriptor dimension per f channel)
    real(dp), dimension(:,:), allocatable :: hh
    real(dp), dimension(:), allocatable :: desc_matrix
    real(dp), dimension(:,:,:), allocatable :: force_desc_matrix
    ! Trace-based descriptors: Tr(h^p) for p=1..power_trace
    real(dp), dimension(:), allocatable :: trace_desc       ! (power_trace)
    real(dp), dimension(:,:,:), allocatable :: force_trace_desc  ! (3, power_trace, dim_nn)
    ! Generic equiv-B spectral-filter descriptors (non-polynomial filters)
    real(dp), dimension(:), allocatable :: equivb_desc            ! (nf)
    real(dp), dimension(:,:,:), allocatable :: force_equivb_desc  ! (3, nf, dim_nn)
    real(dp), dimension(:,:,:), allocatable :: deriv_eigen
    ! Raw eigenvalue derivatives d(lambda_n)/dx (without the 2*lambda factor),
    ! filled when generic spectral filters need them (tb_equivb_rawderiv)
    real(dp), dimension(:,:,:), allocatable :: deriv_eigen_raw
    real(dp), dimension(:,:), allocatable :: u_matrix
    real(dp), dimension(:), allocatable :: eigen_WW
    type(eigen_symmetric) :: eigsolver
    logical :: sinit = .false.
  contains
    procedure :: init            => ss_init
    procedure :: evaluate_energy => ss_evaluate_energy
    procedure :: evaluate_force  => ss_evaluate_force
    procedure :: compute_trace_energy => ss_compute_trace_energy
    procedure :: compute_trace_force  => ss_compute_trace_force
    procedure :: compute_filter_energy => ss_compute_filter_energy
    procedure :: compute_filter_force  => ss_compute_filter_force
  end type ss_hamiltonian

  type(ss_hamiltonian), public :: tb_ss

contains

  subroutine ss_init(this, tb_svd, tb_nn_G, max_neigh_local)
    class(ss_hamiltonian), intent(inout) :: this
    integer, intent(in) :: tb_svd, tb_nn_G, max_neigh_local

    ! Skip full reinit if matrix dimension unchanged (avoid repeated LAPACK workspace query)
    if (this%sinit .and. this%dim_hh == max_neigh_local .and. this%tb_nn_G == tb_nn_G) return

    this%dim_hh  = max_neigh_local
    this%dim_nn  = max_neigh_local
    this%tb_nn_G = tb_nn_G

    if (allocated(this%hh)) deallocate(this%hh)
    allocate(this%hh(this%dim_hh, this%dim_hh))

    if (allocated(this%deriv_eigen)) deallocate(this%deriv_eigen)
    allocate(this%deriv_eigen(3, this%dim_hh, this%dim_nn))
    if (allocated(this%deriv_eigen_raw)) deallocate(this%deriv_eigen_raw)
    allocate(this%deriv_eigen_raw(3, this%dim_hh, this%dim_nn))
    if (allocated(this%u_matrix)) deallocate(this%u_matrix)
    allocate(this%u_matrix(this%dim_hh, this%dim_hh))
    if (allocated(this%eigen_WW)) deallocate(this%eigen_WW)
    allocate(this%eigen_WW(this%dim_hh))

    if (allocated(this%desc_matrix)) deallocate(this%desc_matrix)
    allocate(this%desc_matrix(this%tb_nn_G))
    if (allocated(this%force_desc_matrix)) deallocate(this%force_desc_matrix)
    allocate(this%force_desc_matrix(3, this%tb_nn_G, this%dim_nn))

    call this%eigsolver%init(tb_svd, this%dim_hh)
    this%sinit = .true.
  end subroutine ss_init

  ! ---------------------------------------------------------------------------
  ! Build h^{ss} and compute eigenvalues 
  ! ---------------------------------------------------------------------------
  subroutine ss_evaluate_energy(this, rad, mua_a, mua_neigh, kk, ll, r_central, r_matrix)
    use module_tb_radial, only: tb_channel_radial_t
    class(ss_hamiltonian), intent(inout) :: this
    type(tb_channel_radial_t), intent(in) :: rad
    integer, intent(in) :: mua_a
    integer, dimension(:), intent(in) :: mua_neigh
    integer, intent(in) :: kk, ll
    real(dp), dimension(:), intent(in) :: r_central     ! r_{a,i} for i=1..dim_nn
    real(dp), dimension(0:,0:), intent(in) :: r_matrix  ! r_{i,j} for i,j=0..dim_nn

    integer :: kan, jan
    real(dp) :: rk, rj, rr, gk, gj, ff, gkgj
    ! Precomputed g-values for all neighbours (avoids repeated spline lookups)
    real(dp), dimension(:), allocatable :: g_arr

    ! Precompute g(r_{a,k}) for all neighbours
    allocate(g_arr(this%dim_nn))
    do kan = 1, this%dim_nn
      g_arr(kan) = rad%g_eval_kl(mua_a, mua_neigh(kan), kk, ll, r_central(kan))
    end do

    ! Build h^{ss}
    this%hh(:,:) = 0.0_dp
    do kan = 1, this%dim_nn
      gk = g_arr(kan)
      ! Diagonal: h(k,k) = g_k^2
      this%hh(kan, kan) = gk**2

      do jan = kan+1, this%dim_nn
        gj = g_arr(jan)
        rr = r_matrix(kan, jan)
        ff = rad%f_eval(mua_neigh(kan), mua_neigh(jan), kk, ll, rr)
        gkgj = gk * gj
        this%hh(kan, jan) = ff * gkgj
        this%hh(jan, kan) = this%hh(kan, jan)
      end do
    end do

    deallocate(g_arr)

    ! Eigendecomposition
    call this%eigsolver%evaluate(this%hh)

    ! Extract eigenvalues (eigsolver sorts them)
    if (this%tb_nn_G <= this%dim_hh) then
      this%desc_matrix(1:this%tb_nn_G) = this%eigsolver%eigval(1:this%tb_nn_G)
    else
      this%desc_matrix(1:this%dim_hh) = this%eigsolver%eigval(1:this%dim_hh)
      this%desc_matrix(this%dim_hh+1:this%tb_nn_G) = 0.0_dp
    end if

    this%u_matrix(1:this%dim_hh, 1:this%dim_hh) = this%eigsolver%eigvec(1:this%dim_hh, 1:this%dim_hh)
    this%eigen_WW(1:this%dim_hh) = this%eigsolver%eigWWout(1:this%dim_hh)
  end subroutine ss_evaluate_energy

  ! ---------------------------------------------------------------------------
  ! Force: Hellmann-Feynman for s-s channel
  !
  ! d_{k,chi} lambda_n = 2 u_{n,k} * sum_j u_{n,j} * d_{k,chi} h_{k,j}
  !                     - u_{n,k}^2 * d_{k,chi} h_{k,k}
  !
  ! d_{k,chi} h_{k,j} (j/=k) = f'(r_kj)*hat_kj_chi*g_k*g_j
  !                             + f(r_kj)*g'_k*hat_ka_chi*g_j
  ! d_{k,chi} h_{k,k} = 2*g_k*g'_k*hat_ka_chi
  ! ---------------------------------------------------------------------------
  subroutine ss_evaluate_force(this, rad, mua_a, mua_neigh, kk, ll, &
                                r_central, r_matrix, x_matrix, tmp_dxp)
    use module_tb_radial, only: tb_channel_radial_t
    use module_tbind, only: tb_model_trace, tb_model_equivb, tb_equivb_rawderiv
    class(ss_hamiltonian), intent(inout) :: this
    type(tb_channel_radial_t), intent(in) :: rad
    integer, intent(in) :: mua_a
    integer, dimension(:), intent(in) :: mua_neigh
    integer, intent(in) :: kk, ll
    real(dp), dimension(:), intent(in) :: r_central
    real(dp), dimension(0:,0:), intent(in) :: r_matrix
    real(dp), dimension(:,0:,0:), intent(in) :: x_matrix
    real(dp), dimension(:,0:), intent(in) :: tmp_dxp

    real(dp), dimension(3) :: cosk, tmp3
    real(dp) :: rk, rj, rr, gk, gj, d_gk, ff, d_ff
    integer :: kan, jan, nn, nn_max_deriv
    ! Precomputed g and g' values for all neighbours
    real(dp), dimension(:), allocatable :: g_arr, dg_arr

    ! Precompute g(r_{a,k}) and g'(r_{a,k}) for all neighbours
    allocate(g_arr(this%dim_nn), dg_arr(this%dim_nn))
    do kan = 1, this%dim_nn
      g_arr(kan)  = rad%g_eval_kl(mua_a, mua_neigh(kan), kk, ll, r_central(kan))
      dg_arr(kan) = rad%g_deval_kl(mua_a, mua_neigh(kan), kk, ll, r_central(kan))
    end do

    this%deriv_eigen(:,:,:) = 0.0_dp

    ! When the trace (or equiv-B, which reuses the trace machinery on this
    ! single-sector channel) model is active, we need derivatives for ALL
    ! eigenvalues (not just the truncated tb_nn_G used for lambda descriptors)
    if (tb_model_trace .or. tb_model_equivb) then
      nn_max_deriv = this%dim_hh
    else
      nn_max_deriv = min(this%tb_nn_G, this%dim_hh)
    end if

    do kan = 1, this%dim_nn
      rk = r_central(kan)
      cosk(:) = tmp_dxp(:, kan) / rk
      gk   = g_arr(kan)
      d_gk = dg_arr(kan)

      do nn = 1, nn_max_deriv
        tmp3(:) = 0.0_dp

        ! Diagonal term: d_k h_{k,k} = 2*g_k*g'_k*hat_ka
        ! off-diagonal contribution: 2*u_nk * sum_{j/=k} u_nj * d_k h_{k,j}
        do jan = 1, this%dim_nn
          if (jan == kan) cycle
          gj   = g_arr(jan)
          rr   = r_matrix(kan, jan)
          ff   = rad%f_eval(mua_neigh(kan), mua_neigh(jan), kk, ll, rr)
          d_ff = rad%f_deval(mua_neigh(kan), mua_neigh(jan), kk, ll, rr)
          ! d_k h_{k,j} = f'*hat_kj*g_k*g_j + f*g'_k*hat_ka*g_j
          tmp3(:) = tmp3(:) + this%u_matrix(jan, nn) * ( &
            d_ff * gk * gj * x_matrix(:, kan, jan) + &
            ff * d_gk * gj * cosk(:) )
        end do

        ! Combine: 2*u_nk*tmp3 + u_nk^2 * d_k h_{k,k}
        ! Hellmann-Feynman: d lambda_n / d x_k = u_nk^2 * d_k h_{kk} + 2*u_nk * sum_{j/=k} u_nj * d_k h_{kj}
        ! Then multiply by 2*lambda_n for d(lambda^2)/d x_k
        tmp3(:) = 2.0_dp * this%u_matrix(kan, nn) * tmp3(:) + &
                  this%u_matrix(kan, nn)**2 * 2.0_dp * gk * d_gk * cosk(:)
        if (tb_equivb_rawderiv) this%deriv_eigen_raw(:, nn, kan) = tmp3(:)
        this%deriv_eigen(:, nn, kan) = tmp3(:) * (2.0_dp * this%eigen_WW(nn))
      end do ! nn
    end do ! kan

    ! Pack into force_desc_matrix
    this%force_desc_matrix(:,:,:) = 0.0_dp
    if (this%tb_nn_G <= this%dim_hh) then
      this%force_desc_matrix(1:3, 1:this%tb_nn_G, 1:this%dim_nn) = &
        this%deriv_eigen(1:3, 1:this%tb_nn_G, 1:this%dim_nn)
    else
      this%force_desc_matrix(1:3, 1:this%dim_hh, 1:this%dim_nn) = &
        this%deriv_eigen(1:3, 1:this%dim_hh, 1:this%dim_nn)
    end if
    deallocate(g_arr, dg_arr)
  end subroutine ss_evaluate_force

  ! ---------------------------------------------------------------------------
  ! Trace-based energy descriptors for ss channel
  !
  ! Tr(h^p) = sum_{i=1}^{M} lambda_i^p
  !
  ! where M = dim_hh (all eigenvalues, not truncated).
  ! This uses the eigenvalues already computed by evaluate_energy.
  ! ---------------------------------------------------------------------------
  subroutine ss_compute_trace_energy(this, power_trace)
    class(ss_hamiltonian), intent(inout) :: this
    integer, intent(in) :: power_trace

    integer :: pp, ii, M

    if (allocated(this%trace_desc)) deallocate(this%trace_desc)
    allocate(this%trace_desc(power_trace))

    M = this%dim_hh  ! use ALL eigenvalues for trace

    do pp = 1, power_trace
      this%trace_desc(pp) = 0.0_dp
      do ii = 1, M
        this%trace_desc(pp) = this%trace_desc(pp) + this%eigen_WW(ii)**pp
      end do
    end do
  end subroutine ss_compute_trace_energy

  ! ---------------------------------------------------------------------------
  ! Trace-based force descriptors for ss channel
  !
  ! d/dx_k Tr(h^p) = p * sum_{i=1}^{M} lambda_i^{p-1} * d lambda_i / dx_k
  !
  ! where d lambda_i / dx_k is the raw eigenvalue derivative (before the
  ! 2*lambda_n factor used for lambda^2). We recover it from deriv_eigen which
  ! stores d(lambda^2)/dx_k = 2*lambda_n * d lambda_n / dx_k.
  ! ---------------------------------------------------------------------------
  subroutine ss_compute_trace_force(this, power_trace)
    class(ss_hamiltonian), intent(inout) :: this
    integer, intent(in) :: power_trace

    integer :: pp, ii, kan, M
    real(dp) :: lam, dlam_factor

    if (allocated(this%force_trace_desc)) deallocate(this%force_trace_desc)
    allocate(this%force_trace_desc(3, power_trace, this%dim_nn))

    M = this%dim_hh  ! use ALL eigenvalues

    this%force_trace_desc(:,:,:) = 0.0_dp

    do kan = 1, this%dim_nn
      do pp = 1, power_trace
        do ii = 1, M
          lam = this%eigen_WW(ii)
          ! deriv_eigen stores d(lambda^2)/dx_k = 2*lambda * d lambda/dx_k
          ! so d lambda/dx_k = deriv_eigen / (2*lambda) when lambda /= 0
          ! d Tr(h^p)/dx_k = p * lambda^{p-1} * d lambda/dx_k
          !                = p * lambda^{p-1} * deriv_eigen / (2*lambda)
          !                = p * lambda^{p-2} * deriv_eigen / 2
          ! For p=1: = deriv_eigen / (2*lambda) = d lambda / dx_k
          ! Special case: if lambda = 0, contribution is zero (lambda^{p-1} = 0 for p >= 2,
          ! and for p=1 the deriv_eigen is also zero since it has the 2*lambda factor)
          if (abs(lam) > 1.0e-30_dp) then
            dlam_factor = real(pp, dp) * lam**(pp-2) * 0.5_dp
            this%force_trace_desc(:, pp, kan) = this%force_trace_desc(:, pp, kan) + &
              dlam_factor * this%deriv_eigen(:, ii, kan)
          end if
        end do
      end do
    end do
  end subroutine ss_compute_trace_force

  ! ---------------------------------------------------------------------------
  ! Generic equiv-B spectral-filter descriptors (non-polynomial families):
  !   D_k = snorm * sum_n f_k(lambda_n)      (single s sector on this channel)
  ! ---------------------------------------------------------------------------
  subroutine ss_compute_filter_energy(this, nf, snorm)
    use module_tb_filter, only: tbf_eval
    use module_tbind, only: TB_CH_SS
    class(ss_hamiltonian), intent(inout) :: this
    integer, intent(in) :: nf
    real(dp), intent(in) :: snorm
    integer :: nn
    real(dp), dimension(nf) :: fv, fp

    if (allocated(this%equivb_desc)) deallocate(this%equivb_desc)
    allocate(this%equivb_desc(nf))
    this%equivb_desc(:) = 0.0_dp
    do nn = 1, this%dim_hh
      call tbf_eval(TB_CH_SS, this%eigen_WW(nn), fv, fp)
      this%equivb_desc(1:nf) = this%equivb_desc(1:nf) + fv(1:nf)
    end do
    this%equivb_desc(:) = this%equivb_desc(:)*snorm
  end subroutine ss_compute_filter_energy

  ! d D_k / d x = snorm * sum_n f_k'(lambda_n) * d lambda_n / d x
  ! (uses the raw eigenvalue derivatives, smooth for any differentiable f_k)
  subroutine ss_compute_filter_force(this, nf, snorm)
    use module_tb_filter, only: tbf_eval
    use module_tbind, only: TB_CH_SS
    class(ss_hamiltonian), intent(inout) :: this
    integer, intent(in) :: nf
    real(dp), intent(in) :: snorm
    integer :: nn, kan, k, M
    real(dp), dimension(nf) :: fv
    real(dp), dimension(:,:), allocatable :: fp_tab

    M = this%dim_hh
    if (allocated(this%force_equivb_desc)) deallocate(this%force_equivb_desc)
    allocate(this%force_equivb_desc(3, nf, this%dim_nn))
    this%force_equivb_desc(:,:,:) = 0.0_dp

    allocate(fp_tab(nf, M))
    do nn = 1, M
      call tbf_eval(TB_CH_SS, this%eigen_WW(nn), fv, fp_tab(:, nn))
    end do
    do kan = 1, this%dim_nn
      do nn = 1, M
        do k = 1, nf
          this%force_equivb_desc(:, k, kan) = this%force_equivb_desc(:, k, kan) + &
            fp_tab(k, nn)*this%deriv_eigen_raw(:, nn, kan)
        end do
      end do
    end do
    this%force_equivb_desc(:,:,:) = this%force_equivb_desc(:,:,:)*snorm
    deallocate(fp_tab)
  end subroutine ss_compute_filter_force

end module module_ss_hamiltonian


! ============================================================================
! p-p Hamiltonian  (Section 3 of Tbind_desc.pdf)
!
! 3*N_v x 3*N_v matrix in 3x3 blocks.
! Diagonal block: h_{k,k}^{pp} = g(r_ak)^2 * I_3
! Off-diagonal (Slater-Koster):
!   h_{k alpha, j beta} = [ hat_a*hat_b*(V_sigma - V_pi) + delta_ab*V_pi ] * g_k*g_j
!
! V_sigma and V_pi are DIFFERENT f radial functions:
!   V_sigma(r) = f_eval(mua,mua, kk, ll_sigma, r)
!   V_pi(r)    = f_eval(mua,mua, kk, ll_pi, r)
! where ll_sigma = 2*ll_pair and ll_pi = 2*ll_pair + 1 within the lmax range.
! ============================================================================
module module_pp_hamiltonian
  use iso_fortran_env, only: dp => real64
  use module_eigenvalue_small_matrix, only: eigen_symmetric
  implicit none

  private

  type, public :: pp_hamiltonian
    integer :: dim_hh       ! 3 * max_neigh_local
    integer :: dim_nn       ! max_neigh_local
    integer :: tb_nn_G      ! eigenvalue cutoff
    real(dp), dimension(:,:), allocatable :: hh
    real(dp), dimension(:), allocatable :: desc_matrix
    real(dp), dimension(:,:,:), allocatable :: force_desc_matrix
    ! Trace-based descriptors: Tr(h^p) for p=1..power_trace
    real(dp), dimension(:), allocatable :: trace_desc       ! (power_trace)
    real(dp), dimension(:,:,:), allocatable :: force_trace_desc  ! (3, power_trace, dim_nn)
    ! Generic equiv-B spectral-filter descriptors (non-polynomial filters)
    real(dp), dimension(:), allocatable :: equivb_desc            ! (nf)
    real(dp), dimension(:,:,:), allocatable :: force_equivb_desc  ! (3, nf, dim_nn)
    real(dp), dimension(:,:,:), allocatable :: deriv_eigen
    real(dp), dimension(:,:,:), allocatable :: deriv_eigen_raw  ! raw d(lambda)/dx
    real(dp), dimension(:,:), allocatable :: u_matrix
    real(dp), dimension(:), allocatable :: eigen_WW
    type(eigen_symmetric) :: eigsolver
    logical :: sinit = .false.
  contains
    procedure :: init            => pp_init
    procedure :: evaluate_energy => pp_evaluate_energy
    procedure :: evaluate_force  => pp_evaluate_force
    procedure :: compute_trace_energy => pp_compute_trace_energy
    procedure :: compute_trace_force  => pp_compute_trace_force
    procedure :: compute_filter_energy => pp_compute_filter_energy
    procedure :: compute_filter_force  => pp_compute_filter_force
  end type pp_hamiltonian

  type(pp_hamiltonian), public :: tb_pp

contains

  subroutine pp_init(this, tb_svd, tb_nn_G, max_neigh_local)
    class(pp_hamiltonian), intent(inout) :: this
    integer, intent(in) :: tb_svd, tb_nn_G, max_neigh_local

    ! Skip full reinit if matrix dimension unchanged (avoid repeated LAPACK workspace query)
    if (this%sinit .and. this%dim_nn == max_neigh_local .and. this%tb_nn_G == tb_nn_G) return

    this%dim_hh  = 3 * max_neigh_local
    this%dim_nn  = max_neigh_local
    this%tb_nn_G = tb_nn_G

    if (allocated(this%hh)) deallocate(this%hh)
    allocate(this%hh(this%dim_hh, this%dim_hh))

    if (allocated(this%deriv_eigen)) deallocate(this%deriv_eigen)
    allocate(this%deriv_eigen(3, this%dim_hh, this%dim_nn))
    if (allocated(this%deriv_eigen_raw)) deallocate(this%deriv_eigen_raw)
    allocate(this%deriv_eigen_raw(3, this%dim_hh, this%dim_nn))
    if (allocated(this%u_matrix)) deallocate(this%u_matrix)
    allocate(this%u_matrix(this%dim_hh, this%dim_hh))
    if (allocated(this%eigen_WW)) deallocate(this%eigen_WW)
    allocate(this%eigen_WW(this%dim_hh))

    if (allocated(this%desc_matrix)) deallocate(this%desc_matrix)
    allocate(this%desc_matrix(this%tb_nn_G))
    if (allocated(this%force_desc_matrix)) deallocate(this%force_desc_matrix)
    allocate(this%force_desc_matrix(3, this%tb_nn_G, this%dim_nn))

    call this%eigsolver%init(tb_svd, this%dim_hh)
    this%sinit = .true.
  end subroutine pp_init

  ! ---------------------------------------------------------------------------
  ! Build h^{pp} and compute eigenvalues
  ! ---------------------------------------------------------------------------
  subroutine pp_evaluate_energy(this, rad, mua_a, mua_neigh, kk, ll_sigma, ll_pi, &
                                 r_central, r_matrix, x_matrix)
    use module_tb_radial, only: tb_channel_radial_t
    class(pp_hamiltonian), intent(inout) :: this
    type(tb_channel_radial_t), intent(in) :: rad
    integer, intent(in) :: mua_a
    integer, dimension(:), intent(in) :: mua_neigh
    integer, intent(in) :: kk, ll_sigma, ll_pi
    real(dp), dimension(:), intent(in) :: r_central
    real(dp), dimension(0:,0:), intent(in) :: r_matrix
    real(dp), dimension(:,0:,0:), intent(in) :: x_matrix

    integer :: kan, jan, ix, iy
    real(dp) :: rr, gk, gj, ffs, ffp, gkgj
    real(dp), dimension(3) :: dxp
    ! Precomputed g-values for all neighbours (avoids repeated spline lookups)
    real(dp), dimension(:), allocatable :: g_arr

    ! Precompute g(r_{a,k}) for all neighbours
    ! For g_type=2, g depends on (kk, ll_sigma); for g_type=1, g_eval_kl ignores kk,ll
    allocate(g_arr(this%dim_nn))
    do kan = 1, this%dim_nn
      g_arr(kan) = rad%g_eval_kl(mua_a, mua_neigh(kan), kk, ll_sigma, r_central(kan))
    end do

    this%hh(:,:) = 0.0_dp

    ! Diagonal blocks: g_k^2 * I_3
    do kan = 1, this%dim_nn
      gk = g_arr(kan)
      do ix = 1, 3
        this%hh(3*(kan-1)+ix, 3*(kan-1)+ix) = gk**2
      end do
    end do

    ! Off-diagonal blocks: Slater-Koster pp
    do kan = 1, this%dim_nn
      gk = g_arr(kan)

      do jan = kan+1, this%dim_nn
        gj = g_arr(jan)
        rr = r_matrix(kan, jan)
        dxp(:) = x_matrix(:, kan, jan)

        ffs = rad%f_eval(mua_neigh(kan), mua_neigh(jan), kk, ll_sigma, rr)
        ffp = rad%f_eval(mua_neigh(kan), mua_neigh(jan), kk, ll_pi, rr)
        gkgj = gk * gj

        do ix = 1, 3
          do iy = 1, 3
            if (ix == iy) then
              ! Diagonal: hat_a^2 * V_sigma + (1 - hat_a^2) * V_pi
              this%hh(3*(kan-1)+ix, 3*(jan-1)+iy) = &
                (dxp(ix)**2 * ffs + (1.0_dp - dxp(ix)**2) * ffp) * gkgj
            else
              ! Off-diag: hat_a * hat_b * (V_sigma - V_pi)
              this%hh(3*(kan-1)+ix, 3*(jan-1)+iy) = &
                dxp(ix) * dxp(iy) * (ffs - ffp) * gkgj
            end if
            this%hh(3*(jan-1)+iy, 3*(kan-1)+ix) = this%hh(3*(kan-1)+ix, 3*(jan-1)+iy)
          end do
        end do
      end do ! jan
    end do ! kan

    ! Eigendecomposition
    call this%eigsolver%evaluate(this%hh)

    if (this%tb_nn_G <= this%dim_hh) then
      this%desc_matrix(1:this%tb_nn_G) = this%eigsolver%eigval(1:this%tb_nn_G)
    else
      this%desc_matrix(1:this%dim_hh) = this%eigsolver%eigval(1:this%dim_hh)
      this%desc_matrix(this%dim_hh+1:this%tb_nn_G) = 0.0_dp
    end if

    this%u_matrix(1:this%dim_hh, 1:this%dim_hh) = this%eigsolver%eigvec(1:this%dim_hh, 1:this%dim_hh)
    this%eigen_WW(1:this%dim_hh) = this%eigsolver%eigWWout(1:this%dim_hh)

    deallocate(g_arr)
  end subroutine pp_evaluate_energy

  ! ---------------------------------------------------------------------------
  ! Force: Hellmann-Feynman for p-p channel
  !
  ! Similar structure to ss but with 3x3 block derivatives involving
  ! direction cosines and Slater-Koster decomposition.
  ! ---------------------------------------------------------------------------
  subroutine pp_evaluate_force(this, rad, mua_a, mua_neigh, kk, ll_sigma, ll_pi, &
                                r_central, r_matrix, x_matrix, tmp_dxp)
    use module_tb_radial, only: tb_channel_radial_t
    use module_tbind, only: tb_model_trace, tb_model_equivb, tb_equivb_rawderiv
    class(pp_hamiltonian), intent(inout) :: this
    type(tb_channel_radial_t), intent(in) :: rad
    integer, intent(in) :: mua_a
    integer, dimension(:), intent(in) :: mua_neigh
    integer, intent(in) :: kk, ll_sigma, ll_pi
    real(dp), dimension(:), intent(in) :: r_central
    real(dp), dimension(0:,0:), intent(in) :: r_matrix
    real(dp), dimension(:,0:,0:), intent(in) :: x_matrix
    real(dp), dimension(:,0:), intent(in) :: tmp_dxp

    ! d_chi h_{3(k-1)+ix, 3(j-1)+iy}  for this particular derivative atom k
    real(dp) :: rk, rr, gk, gj, d_gk, ffs, ffp, dffs, dffp, gkgj
    real(dp), dimension(3) :: dxp, cosk
    real(dp) :: sk_pp_diag, sk_pp_off
    real(dp), dimension(3) :: d_sk_val, d_hat_ix, d_hat_iy
    real(dp), dimension(3) :: tmp_sum, diag_part
    integer :: kan, jan, nn, ix, iy, nn_max_deriv
    ! Precomputed g-values and derivatives for all neighbours
    real(dp), dimension(:), allocatable :: g_arr, dg_arr

    ! Precompute g(r_{a,k}) and g'(r_{a,k}) for all neighbours
    ! For g_type=2, g depends on (kk, ll_sigma); for g_type=1, g_eval_kl ignores kk,ll
    allocate(g_arr(this%dim_nn), dg_arr(this%dim_nn))
    do kan = 1, this%dim_nn
      g_arr(kan)  = rad%g_eval_kl(mua_a, mua_neigh(kan), kk, ll_sigma, r_central(kan))
      dg_arr(kan) = rad%g_deval_kl(mua_a, mua_neigh(kan), kk, ll_sigma, r_central(kan))
    end do

    this%deriv_eigen(:,:,:) = 0.0_dp

    ! When trace (or equiv-B: scaled trace on this single-sector channel) is
    ! active, we need derivatives for ALL eigenvalues
    if (tb_model_trace .or. tb_model_equivb) then
      nn_max_deriv = this%dim_hh
    else
      nn_max_deriv = min(this%tb_nn_G, this%dim_hh)
    end if

    do kan = 1, this%dim_nn
      rk = r_central(kan)
      cosk(:) = tmp_dxp(:, kan) / rk
      gk   = g_arr(kan)
      d_gk = dg_arr(kan)

      do nn = 1, nn_max_deriv

        ! Accumulate:
        !   out_term = 2 * sum_{jrow} u_{jrow,nn} * d_{k,chi} h_{3(k-1)+ix, jrow} * u_{3(k-1)+ix, nn}
        !   diag_term = sum_{ix,iy} u_{3(k-1)+ix,nn} * u_{3(k-1)+iy,nn} * d_{k,chi} h_{3(k-1)+ix, 3(k-1)+iy}

        ! Part 1: diagonal block derivative
        ! d_k [g_k^2 * delta_ab] = 2*g*g'*hat_ka * delta_ab
        diag_part(:) = 0.0_dp
        do ix = 1, 3
          diag_part(:) = diag_part(:) + this%u_matrix(3*(kan-1)+ix, nn)**2 * &
                         2.0_dp * gk * d_gk * cosk(:)
        end do

        ! Part 2: off-diagonal blocks
        tmp_sum(:) = 0.0_dp
        do jan = 1, this%dim_nn
          if (jan == kan) cycle
          gj   = g_arr(jan)
          rr   = r_matrix(kan, jan)
          ffs  = rad%f_eval(mua_neigh(kan), mua_neigh(jan), kk, ll_sigma, rr)
          dffs = rad%f_deval(mua_neigh(kan), mua_neigh(jan), kk, ll_sigma, rr)
          ffp  = rad%f_eval(mua_neigh(kan), mua_neigh(jan), kk, ll_pi, rr)
          dffp = rad%f_deval(mua_neigh(kan), mua_neigh(jan), kk, ll_pi, rr)
          dxp(:) = x_matrix(:, kan, jan)
          gkgj = gk * gj

          do ix = 1, 3
            ! Compute d_chi hat_ix = (delta_{chi,ix} - hat_ix*hat_chi) / r_kj
            d_hat_ix(:) = -dxp(ix) * dxp(:) / rr
            d_hat_ix(ix) = d_hat_ix(ix) + 1.0_dp / rr

            do iy = 1, 3
              if (ix == iy) then
                sk_pp_diag = dxp(ix)**2 * ffs + (1.0_dp - dxp(ix)**2) * ffp
                ! d_chi [SK * g_k*g_j]
                ! = d_chi[SK] * g_k*g_j + SK * g'_k*hat_ka_chi * g_j
                d_sk_val(:) = (2.0_dp * dxp(ix) * d_hat_ix(:) * (ffs - ffp) + &
                               dxp(ix)**2 * dffs * dxp(:) + (1.0_dp - dxp(ix)**2) * dffp * dxp(:)) * gkgj + &
                               sk_pp_diag * d_gk * gj * cosk(:)
              else
                sk_pp_off = dxp(ix) * dxp(iy) * (ffs - ffp)
                ! Compute d_chi hat_iy = (delta_{chi,iy} - hat_iy*hat_chi) / r_kj
                d_hat_iy(:) = -dxp(iy) * dxp(:) / rr
                d_hat_iy(iy) = d_hat_iy(iy) + 1.0_dp / rr
                ! d_chi [SK * g_k*g_j]
                d_sk_val(:) = ( &
                  d_hat_ix(:) * dxp(iy) * (ffs - ffp) + &
                  dxp(ix) * d_hat_iy(:) * (ffs - ffp) + &
                  dxp(ix) * dxp(iy) * (dffs - dffp) * dxp(:) &
                  ) * gkgj + &
                  sk_pp_off * d_gk * gj * cosk(:)
              end if

              ! Accumulate: u_{3(k-1)+ix, nn} * d_chi h_{3(k-1)+ix, 3(j-1)+iy} * u_{3(j-1)+iy, nn}
              tmp_sum(:) = tmp_sum(:) + &
                this%u_matrix(3*(kan-1)+ix, nn) * d_sk_val(:) * this%u_matrix(3*(jan-1)+iy, nn)
            end do ! iy
          end do ! ix
        end do ! jan

        ! Also add diagonal block contribution to tmp_sum:
        ! d_k h_{kk block} = 2*g*g'*hat_ka * I_3
        do ix = 1, 3
          tmp_sum(:) = tmp_sum(:) + &
            this%u_matrix(3*(kan-1)+ix, nn) * 2.0_dp * gk * d_gk * cosk(:) * this%u_matrix(3*(kan-1)+ix, nn)
        end do

        ! Hellmann-Feynman: d lambda_n = 2*sum - diag  
        ! Actually for symmetric: d lambda_n = sum_i sum_j u_in d h_ij u_jn
        ! The full expression is: d lambda_n = sum_{all i,j} u_in * d_k h_ij * u_jn
        ! But only terms where row i has k-index contribute (since d_k h only affects rows/cols with k)
        ! For SS: d_k lambda_n = 2*u_kn * sum_j u_jn * d_k h_kj - u_kn^2 * d_k h_kk
        !       = sum_j u_jn * d_k h_kj * u_kn + sum_j u_kn * d_k h_jk * u_jn - u_kn^2 * d_k h_kk
        ! For PP, the k-index covers rows 3(k-1)+1..3(k-1)+3
        ! d lambda_n = 2 * sum_{ix} sum_{jrow} u_{3(k-1)+ix,n} * d_k h_{3(k-1)+ix, jrow} * u_{jrow,n}
        !            - sum_{ix} sum_{iy} u_{3(k-1)+ix,n} * d_k h_{3(k-1)+ix, 3(k-1)+iy} * u_{3(k-1)+iy,n}

        tmp_sum(:) = 2.0_dp * tmp_sum(:) - diag_part(:)
        if (tb_equivb_rawderiv) this%deriv_eigen_raw(:, nn, kan) = tmp_sum(:)
        this%deriv_eigen(:, nn, kan) = tmp_sum(:) * (2.0_dp * this%eigen_WW(nn))
      end do ! nn
    end do ! kan

    ! Pack
    this%force_desc_matrix(:,:,:) = 0.0_dp
    if (this%tb_nn_G <= this%dim_hh) then
      this%force_desc_matrix(1:3, 1:this%tb_nn_G, 1:this%dim_nn) = &
        this%deriv_eigen(1:3, 1:this%tb_nn_G, 1:this%dim_nn)
    else
      this%force_desc_matrix(1:3, 1:this%dim_hh, 1:this%dim_nn) = &
        this%deriv_eigen(1:3, 1:this%dim_hh, 1:this%dim_nn)
    end if

    deallocate(g_arr, dg_arr)
  end subroutine pp_evaluate_force

  ! ---------------------------------------------------------------------------
  ! Trace-based energy descriptors for pp channel
  !
  ! Tr(h^p) = sum_{i=1}^{M} lambda_i^p   where M = dim_hh = 3*dim_nn
  ! ---------------------------------------------------------------------------
  subroutine pp_compute_trace_energy(this, power_trace)
    class(pp_hamiltonian), intent(inout) :: this
    integer, intent(in) :: power_trace

    integer :: pp, ii, M

    if (allocated(this%trace_desc)) deallocate(this%trace_desc)
    allocate(this%trace_desc(power_trace))

    M = this%dim_hh

    do pp = 1, power_trace
      this%trace_desc(pp) = 0.0_dp
      do ii = 1, M
        this%trace_desc(pp) = this%trace_desc(pp) + this%eigen_WW(ii)**pp
      end do
    end do
  end subroutine pp_compute_trace_energy

  ! ---------------------------------------------------------------------------
  ! Trace-based force descriptors for pp channel
  !
  ! d/dx_k Tr(h^p) = p * sum_{i=1}^{M} lambda_i^{p-1} * d lambda_i / dx_k
  ! ---------------------------------------------------------------------------
  subroutine pp_compute_trace_force(this, power_trace)
    class(pp_hamiltonian), intent(inout) :: this
    integer, intent(in) :: power_trace

    integer :: pp, ii, kan, M
    real(dp) :: lam, dlam_factor

    if (allocated(this%force_trace_desc)) deallocate(this%force_trace_desc)
    allocate(this%force_trace_desc(3, power_trace, this%dim_nn))

    M = this%dim_hh

    this%force_trace_desc(:,:,:) = 0.0_dp

    do kan = 1, this%dim_nn
      do pp = 1, power_trace
        do ii = 1, M
          lam = this%eigen_WW(ii)
          ! deriv_eigen stores d(lambda^2)/dx_k = 2*lambda * d lambda/dx_k
          ! d Tr(h^p)/dx_k = p * lambda^{p-2} * deriv_eigen / 2
          if (abs(lam) > 1.0e-30_dp) then
            dlam_factor = real(pp, dp) * lam**(pp-2) * 0.5_dp
            this%force_trace_desc(:, pp, kan) = this%force_trace_desc(:, pp, kan) + &
              dlam_factor * this%deriv_eigen(:, ii, kan)
          end if
        end do
      end do
    end do
  end subroutine pp_compute_trace_force

  ! ---------------------------------------------------------------------------
  ! Generic equiv-B spectral-filter descriptors (single p sector; norm 1/sqrt3)
  ! ---------------------------------------------------------------------------
  subroutine pp_compute_filter_energy(this, nf, snorm)
    use module_tb_filter, only: tbf_eval
    use module_tbind, only: TB_CH_PP
    class(pp_hamiltonian), intent(inout) :: this
    integer, intent(in) :: nf
    real(dp), intent(in) :: snorm
    integer :: nn
    real(dp), dimension(nf) :: fv, fp

    if (allocated(this%equivb_desc)) deallocate(this%equivb_desc)
    allocate(this%equivb_desc(nf))
    this%equivb_desc(:) = 0.0_dp
    do nn = 1, this%dim_hh
      call tbf_eval(TB_CH_PP, this%eigen_WW(nn), fv, fp)
      this%equivb_desc(1:nf) = this%equivb_desc(1:nf) + fv(1:nf)
    end do
    this%equivb_desc(:) = this%equivb_desc(:)*snorm
  end subroutine pp_compute_filter_energy

  subroutine pp_compute_filter_force(this, nf, snorm)
    use module_tb_filter, only: tbf_eval
    use module_tbind, only: TB_CH_PP
    class(pp_hamiltonian), intent(inout) :: this
    integer, intent(in) :: nf
    real(dp), intent(in) :: snorm
    integer :: nn, kan, k, M
    real(dp), dimension(nf) :: fv
    real(dp), dimension(:,:), allocatable :: fp_tab

    M = this%dim_hh
    if (allocated(this%force_equivb_desc)) deallocate(this%force_equivb_desc)
    allocate(this%force_equivb_desc(3, nf, this%dim_nn))
    this%force_equivb_desc(:,:,:) = 0.0_dp

    allocate(fp_tab(nf, M))
    do nn = 1, M
      call tbf_eval(TB_CH_PP, this%eigen_WW(nn), fv, fp_tab(:, nn))
    end do
    do kan = 1, this%dim_nn
      do nn = 1, M
        do k = 1, nf
          this%force_equivb_desc(:, k, kan) = this%force_equivb_desc(:, k, kan) + &
            fp_tab(k, nn)*this%deriv_eigen_raw(:, nn, kan)
        end do
      end do
    end do
    this%force_equivb_desc(:,:,:) = this%force_equivb_desc(:,:,:)*snorm
    deallocate(fp_tab)
  end subroutine pp_compute_filter_force

end module module_pp_hamiltonian


! ============================================================================
! s-p coupling Hamiltonian (Section 4 of Tbind_desc.pdf, eq. 52)
!
! 4 N_v x 4 N_v real-symmetric matrix in 4x4 blocks; per-neighbour orbital
! basis is (s, p_x, p_y, p_z). For one descriptor evaluation, the sp channel
! consumes FOUR consecutive angular indices on the sp radial basis:
!
!   ll_ss   = 4m     -> f_ss   (s-s bond integral)
!   ll_sp   = 4m+1   -> f_sp   (s-p coupling)
!   ll_pps  = 4m+2   -> f_pps  (pp sigma)
!   ll_ppp  = 4m+3   -> f_ppp  (pp pi)
!
! On-site block  (i = j = k):   diag( 4 g_k^2, 4 g_k^2, 4 g_k^2, 4 g_k^2 )
! Off-site block (i = k, j /= k), with chi_kj(:) = (r_k - r_j)/|r_k - r_j|:
!
!   block(s,s)     =  A_ss
!   block(s,p_b)   =  chi_kj(b) * A_sp
!   block(p_a,s)   = -chi_kj(a) * A_sp
!   block(p_a,p_a) =  chi_kj(a)^2 * A_sig + (1 - chi_kj(a)^2) * A_pi
!   block(p_a,p_b) =  chi_kj(a) * chi_kj(b) * (A_sig - A_pi)         (a /= b)
!
! with A_X = f_X(r_kj) * g(r_ak) * g(r_aj). A single g function (the sp
! channel's g) is used for all four amplitudes, consistent with how the pp
! channel collapses g_p_sigma / g_p_pi into one g.
!
! Forces use the Hellmann-Feynman formula (eq. 54 of Tbind_desc.pdf), with
! the sparse-block decomposition Q_k = {4(k-1)+1 .. 4(k-1)+4}.
! ============================================================================
module module_sp_hamiltonian
  use iso_fortran_env, only: dp => real64
  use module_eigenvalue_small_matrix, only: eigen_symmetric
  implicit none

  private

  ! 1/sqrt(3): Model A normalisation of the p-sector trace, Tr[Pi_p h^q]/sqrt(2l+1)
  real(dp), parameter :: inv_sqrt3 = 0.57735026918962576_dp

  type, public :: sp_hamiltonian
    integer :: dim_hh       ! 4 * max_neigh_local
    integer :: dim_nn       ! max_neigh_local
    integer :: tb_nn_G      ! eigenvalue cutoff
    real(dp), dimension(:,:), allocatable :: hh
    real(dp), dimension(:), allocatable :: desc_matrix
    real(dp), dimension(:,:,:), allocatable :: force_desc_matrix
    real(dp), dimension(:), allocatable :: trace_desc
    real(dp), dimension(:,:,:), allocatable :: force_trace_desc
    ! equiv-B (Model A of docs/Tbind_desc/main.tex): orbital-sector traces
    ! Tr[Pi_l h^q]/sqrt(2l+1) with polynomial filters. The sp channel has two
    ! sectors (s, p); descriptor layout is sector-inner: (q-1)*2+1 = s, +2 = p.
    real(dp), dimension(:), allocatable :: equivb_desc            ! (2*power)
    real(dp), dimension(:,:,:), allocatable :: force_equivb_desc  ! (3, 2*power, dim_nn)
    real(dp), dimension(:,:,:), allocatable :: WW_equivb          ! (dim_hh, dim_hh, 2*power)
    real(dp), dimension(:,:,:), allocatable :: deriv_eigen
    real(dp), dimension(:,:), allocatable :: u_matrix
    real(dp), dimension(:), allocatable :: eigen_WW
    type(eigen_symmetric) :: eigsolver
    logical :: sinit = .false.
  contains
    procedure :: init                  => sp_init
    procedure :: evaluate_energy       => sp_evaluate_energy
    procedure :: evaluate_force        => sp_evaluate_force
    procedure :: compute_trace_energy  => sp_compute_trace_energy
    procedure :: compute_trace_force   => sp_compute_trace_force
    procedure :: compute_equivb_energy => sp_compute_equivb_energy
    procedure :: compute_filter_energy => sp_compute_filter_energy
  end type sp_hamiltonian

  type(sp_hamiltonian), public :: tb_sp

contains

  subroutine sp_init(this, tb_svd, tb_nn_G, max_neigh_local)
    class(sp_hamiltonian), intent(inout) :: this
    integer, intent(in) :: tb_svd, tb_nn_G, max_neigh_local

    if (this%sinit .and. this%dim_nn == max_neigh_local .and. this%tb_nn_G == tb_nn_G) return

    this%dim_hh  = 4 * max_neigh_local
    this%dim_nn  = max_neigh_local
    this%tb_nn_G = tb_nn_G

    if (allocated(this%hh)) deallocate(this%hh)
    allocate(this%hh(this%dim_hh, this%dim_hh))

    if (allocated(this%deriv_eigen)) deallocate(this%deriv_eigen)
    allocate(this%deriv_eigen(3, this%dim_hh, this%dim_nn))
    if (allocated(this%u_matrix)) deallocate(this%u_matrix)
    allocate(this%u_matrix(this%dim_hh, this%dim_hh))
    if (allocated(this%eigen_WW)) deallocate(this%eigen_WW)
    allocate(this%eigen_WW(this%dim_hh))

    if (allocated(this%desc_matrix)) deallocate(this%desc_matrix)
    allocate(this%desc_matrix(this%tb_nn_G))
    if (allocated(this%force_desc_matrix)) deallocate(this%force_desc_matrix)
    allocate(this%force_desc_matrix(3, this%tb_nn_G, this%dim_nn))

    call this%eigsolver%init(tb_svd, this%dim_hh)
    this%sinit = .true.
  end subroutine sp_init

  ! ---------------------------------------------------------------------------
  ! Build h^{sp} (4 N_v x 4 N_v) and diagonalise
  ! ---------------------------------------------------------------------------
  subroutine sp_evaluate_energy(this, rad, mua_a, mua_neigh, kk, ll_ss, ll_sp, ll_pps, ll_ppp, &
                                 r_central, r_matrix, x_matrix)
    use module_tb_radial, only: tb_channel_radial_t
    class(sp_hamiltonian), intent(inout) :: this
    type(tb_channel_radial_t), intent(in) :: rad
    integer, intent(in) :: mua_a
    integer, dimension(:), intent(in) :: mua_neigh
    integer, intent(in) :: kk, ll_ss, ll_sp, ll_pps, ll_ppp
    real(dp), dimension(:), intent(in) :: r_central
    real(dp), dimension(0:,0:), intent(in) :: r_matrix
    real(dp), dimension(:,0:,0:), intent(in) :: x_matrix

    integer :: kan, jan, a, b, ia, ib, off_k, off_j
    real(dp) :: rr, gk, gj, gkgj, A_ss, A_sp, A_sig, A_pi
    real(dp), dimension(3) :: chikj
    real(dp), dimension(:), allocatable :: g_arr

    ! Precompute g(r_ak) for each neighbour (uses ll_ss for g_type=2; ignored for g_type=1)
    allocate(g_arr(this%dim_nn))
    do kan = 1, this%dim_nn
      g_arr(kan) = rad%g_eval_kl(mua_a, mua_neigh(kan), kk, ll_ss, r_central(kan))
    end do

    this%hh(:,:) = 0.0_dp

    ! Diagonal (on-site) blocks: 4 * g_k^2 on every diagonal entry of the 4x4 block
    do kan = 1, this%dim_nn
      gk = g_arr(kan)
      off_k = 4*(kan-1)
      this%hh(off_k+1, off_k+1) = 4.0_dp * gk*gk
      this%hh(off_k+2, off_k+2) = 4.0_dp * gk*gk
      this%hh(off_k+3, off_k+3) = 4.0_dp * gk*gk
      this%hh(off_k+4, off_k+4) = 4.0_dp * gk*gk
    end do

    ! Off-diagonal blocks (upper triangle); symmetric counterpart filled by transpose
    do kan = 1, this%dim_nn
      gk = g_arr(kan)
      off_k = 4*(kan-1)
      do jan = kan+1, this%dim_nn
        gj = g_arr(jan)
        off_j = 4*(jan-1)
        rr = r_matrix(kan, jan)
        chikj(:) = x_matrix(:, kan, jan)
        gkgj = gk * gj

        A_ss  = rad%f_eval(mua_neigh(kan), mua_neigh(jan), kk, ll_ss,  rr) * gkgj
        A_sp  = rad%f_eval(mua_neigh(kan), mua_neigh(jan), kk, ll_sp,  rr) * gkgj
        A_sig = rad%f_eval(mua_neigh(kan), mua_neigh(jan), kk, ll_pps, rr) * gkgj
        A_pi  = rad%f_eval(mua_neigh(kan), mua_neigh(jan), kk, ll_ppp, rr) * gkgj

        ! (s, s)
        this%hh(off_k+1, off_j+1) = A_ss
        ! (s, p_b) and (p_a, s)
        do b = 1, 3
          this%hh(off_k+1,   off_j+1+b) =  chikj(b) * A_sp
          this%hh(off_k+1+b, off_j+1  ) = -chikj(b) * A_sp
        end do
        ! (p_a, p_b)
        do a = 1, 3
          do b = 1, 3
            if (a == b) then
              this%hh(off_k+1+a, off_j+1+b) = chikj(a)*chikj(a) * A_sig + &
                                              (1.0_dp - chikj(a)*chikj(a)) * A_pi
            else
              this%hh(off_k+1+a, off_j+1+b) = chikj(a) * chikj(b) * (A_sig - A_pi)
            end if
          end do
        end do

        ! Symmetric transpose block (j,k) = (k,j)^T
        do ia = 1, 4
          do ib = 1, 4
            this%hh(off_j+ib, off_k+ia) = this%hh(off_k+ia, off_j+ib)
          end do
        end do
      end do
    end do

    deallocate(g_arr)

    ! Eigendecomposition (sorted by descending |lambda| inside eigsolver)
    call this%eigsolver%evaluate(this%hh)

    if (this%tb_nn_G <= this%dim_hh) then
      this%desc_matrix(1:this%tb_nn_G) = this%eigsolver%eigval(1:this%tb_nn_G)
    else
      this%desc_matrix(1:this%dim_hh) = this%eigsolver%eigval(1:this%dim_hh)
      this%desc_matrix(this%dim_hh+1:this%tb_nn_G) = 0.0_dp
    end if

    this%u_matrix(1:this%dim_hh, 1:this%dim_hh) = this%eigsolver%eigvec(1:this%dim_hh, 1:this%dim_hh)
    this%eigen_WW(1:this%dim_hh) = this%eigsolver%eigWWout(1:this%dim_hh)
  end subroutine sp_evaluate_energy

  ! ---------------------------------------------------------------------------
  ! Hellmann-Feynman force for s-p coupling (eq. 54 of Tbind_desc.pdf)
  !
  ! For each neighbour k and cartesian chi:
  !   d lambda_n / d k_chi = 2 * sum_{p in Q_k} sum_{q=1..4 N_v} u_pn u_qn (d_kchi h)_pq
  !                          -     sum_{p in Q_k} sum_{q in Q_k}  u_pn u_qn (d_kchi h)_pq
  !
  ! The on-site block contributes only to the Q_k x Q_k double sum (diagonal in
  ! the 4x4 block); off-site blocks contribute to the first term via j /= k.
  ! Translational invariance for k = a (central atom) is enforced in the caller.
  ! ---------------------------------------------------------------------------
  subroutine sp_evaluate_force(this, rad, mua_a, mua_neigh, kk, ll_ss, ll_sp, ll_pps, ll_ppp, &
                                r_central, r_matrix, x_matrix, tmp_dxp)
    use module_tb_radial, only: tb_channel_radial_t
    use module_tbind, only: tb_model_trace, tb_model_equivb, tb_power_trace, TB_CH_SP, &
                            tb_nf_equivb, tb_equivb_poly
    class(sp_hamiltonian), intent(inout) :: this
    type(tb_channel_radial_t), intent(in) :: rad
    integer, intent(in) :: mua_a
    integer, dimension(:), intent(in) :: mua_neigh
    integer, intent(in) :: kk, ll_ss, ll_sp, ll_pps, ll_ppp
    real(dp), dimension(:), intent(in) :: r_central
    real(dp), dimension(0:,0:), intent(in) :: r_matrix
    real(dp), dimension(:,0:,0:), intent(in) :: x_matrix
    real(dp), dimension(:,0:), intent(in) :: tmp_dxp

    real(dp), dimension(:), allocatable :: g_arr, dg_arr
    ! dh_block(alpha, beta, chi, jan) = d_kchi h_{Q_k[alpha], Q_jan[beta]} for jan /= kan
    real(dp), dimension(:,:,:,:), allocatable :: dh_block
    real(dp), dimension(3) :: cosk, chikj, d_alpha, d_beta
    real(dp) :: rk, rr, gk, gj, dgk, inv_rkj
    real(dp) :: f_ss, df_ss, f_sp, df_sp, f_sg, df_sg, f_pi, df_pi
    real(dp) :: A_ss, A_sp, A_sig, A_pi
    real(dp), dimension(3) :: dA_ss, dA_sp, dA_sig, dA_pi
    real(dp) :: sum_u2, u_alpha, u_beta, wdiag
    real(dp), dimension(3) :: diag_part, tmp_sum
    integer :: kan, jan, nn, a, b, ialp, ibet, off_k, off_j, nn_max_deriv
    integer :: iq, power_equivb

    allocate(g_arr(this%dim_nn), dg_arr(this%dim_nn))
    do kan = 1, this%dim_nn
      g_arr(kan)  = rad%g_eval_kl(mua_a, mua_neigh(kan), kk, ll_ss, r_central(kan))
      dg_arr(kan) = rad%g_deval_kl(mua_a, mua_neigh(kan), kk, ll_ss, r_central(kan))
    end do

    this%deriv_eigen(:,:,:) = 0.0_dp

    if (tb_model_trace) then
      nn_max_deriv = this%dim_hh
    else
      nn_max_deriv = min(this%tb_nn_G, this%dim_hh)
    end if

    ! equiv-B: build the filtered W matrices from the eigen-decomposition
    ! (the W contraction below does NOT use deriv_eigen, so nn_max_deriv is
    ! unaffected) and allocate the force store.  The legacy raw-polynomial
    ! path uses the exact phi_q recursion; every other spectral-filter family
    ! goes through the generic divided-difference builder.
    power_equivb = 0
    if (tb_model_equivb) then
      power_equivb = tb_nf_equivb(TB_CH_SP)
      if (tb_equivb_poly) then
        call sp_build_equivb_W(this, power_equivb)
      else
        call sp_build_filter_W(this, power_equivb)
      end if
      if (allocated(this%force_equivb_desc)) deallocate(this%force_equivb_desc)
      allocate(this%force_equivb_desc(3, 2*power_equivb, this%dim_nn))
      this%force_equivb_desc(:,:,:) = 0.0_dp
    end if

    allocate(dh_block(0:3, 0:3, 3, this%dim_nn))

    do kan = 1, this%dim_nn
      rk = r_central(kan)
      cosk(:) = tmp_dxp(:, kan) / rk
      gk  = g_arr(kan)
      dgk = dg_arr(kan)
      off_k = 4*(kan-1)

      ! ---- Precompute the 4x4 block derivative d_kchi h_{k,j} for each j /= k ----
      dh_block(:,:,:,:) = 0.0_dp
      do jan = 1, this%dim_nn
        if (jan == kan) cycle
        gj = g_arr(jan)
        rr = r_matrix(kan, jan)
        chikj(:) = x_matrix(:, kan, jan)
        inv_rkj = 1.0_dp / rr

        f_ss  = rad%f_eval (mua_neigh(kan), mua_neigh(jan), kk, ll_ss,  rr)
        df_ss = rad%f_deval(mua_neigh(kan), mua_neigh(jan), kk, ll_ss,  rr)
        f_sp  = rad%f_eval (mua_neigh(kan), mua_neigh(jan), kk, ll_sp,  rr)
        df_sp = rad%f_deval(mua_neigh(kan), mua_neigh(jan), kk, ll_sp,  rr)
        f_sg  = rad%f_eval (mua_neigh(kan), mua_neigh(jan), kk, ll_pps, rr)
        df_sg = rad%f_deval(mua_neigh(kan), mua_neigh(jan), kk, ll_pps, rr)
        f_pi  = rad%f_eval (mua_neigh(kan), mua_neigh(jan), kk, ll_ppp, rr)
        df_pi = rad%f_deval(mua_neigh(kan), mua_neigh(jan), kk, ll_ppp, rr)

        A_ss  = f_ss * gk * gj
        A_sp  = f_sp * gk * gj
        A_sig = f_sg * gk * gj
        A_pi  = f_pi * gk * gj
        ! d_kchi A = (df * chikj(chi) * gk + f * dg_k * cosk(chi)) * gj
        dA_ss (:) = (df_ss * gk * chikj(:) + f_ss * dgk * cosk(:)) * gj
        dA_sp (:) = (df_sp * gk * chikj(:) + f_sp * dgk * cosk(:)) * gj
        dA_sig(:) = (df_sg * gk * chikj(:) + f_sg * dgk * cosk(:)) * gj
        dA_pi (:) = (df_pi * gk * chikj(:) + f_pi * dgk * cosk(:)) * gj

        ! (s,s) entry
        dh_block(0, 0, :, jan) = dA_ss(:)

        ! (s, p_b) and (p_b, s) entries
        do b = 1, 3
          ! d_kchi chikj(b) = (delta_{b,chi} - chikj(b)*chikj(chi)) / r_kj
          d_beta(:) = -chikj(b) * chikj(:) * inv_rkj
          d_beta(b) = d_beta(b) + inv_rkj
          dh_block(0, b, :, jan) =  d_beta(:) * A_sp + chikj(b) * dA_sp(:)
          dh_block(b, 0, :, jan) = -d_beta(:) * A_sp - chikj(b) * dA_sp(:)
        end do

        ! (p_a, p_b) entries
        do a = 1, 3
          d_alpha(:) = -chikj(a) * chikj(:) * inv_rkj
          d_alpha(a) = d_alpha(a) + inv_rkj
          do b = 1, 3
            if (a == b) then
              dh_block(a, b, :, jan) = 2.0_dp * chikj(a) * d_alpha(:) * (A_sig - A_pi) &
                                     + chikj(a)*chikj(a) * dA_sig(:) &
                                     + (1.0_dp - chikj(a)*chikj(a)) * dA_pi(:)
            else
              d_beta(:) = -chikj(b) * chikj(:) * inv_rkj
              d_beta(b) = d_beta(b) + inv_rkj
              dh_block(a, b, :, jan) = (d_alpha(:) * chikj(b) + chikj(a) * d_beta(:)) * (A_sig - A_pi) &
                                     + chikj(a) * chikj(b) * (dA_sig(:) - dA_pi(:))
            end if
          end do
        end do
      end do  ! jan (precompute dh)

      ! ---- Hellmann-Feynman contraction per eigenvalue ----
      do nn = 1, nn_max_deriv
        ! On-site (Q_k x Q_k) diagonal block: d_kchi (4 g_k^2) = 8 g_k dg_k cosk(chi)
        sum_u2 = 0.0_dp
        do ialp = 0, 3
          sum_u2 = sum_u2 + this%u_matrix(off_k+1+ialp, nn)**2
        end do
        diag_part(:) = sum_u2 * 8.0_dp * gk * dgk * cosk(:)

        ! tmp_sum accumulates "sum1" = sum_{p in Q_k, q=1..4N_v} u_pn u_qn (dh)_pq
        ! Start with the Q_k x Q_k contribution from the on-site block
        tmp_sum(:) = diag_part(:)

        ! Off-site blocks (j /= k)
        do jan = 1, this%dim_nn
          if (jan == kan) cycle
          off_j = 4*(jan-1)
          do ialp = 0, 3
            u_alpha = this%u_matrix(off_k+1+ialp, nn)
            do ibet = 0, 3
              u_beta = this%u_matrix(off_j+1+ibet, nn)
              tmp_sum(:) = tmp_sum(:) + u_alpha * dh_block(ialp, ibet, :, jan) * u_beta
            end do
          end do
        end do

        ! d lambda_n / d k_chi = 2 * sum1 - sum2; descriptor stores lambda^2 -> *2*lambda
        this%deriv_eigen(:, nn, kan) = (2.0_dp * tmp_sum(:) - diag_part(:)) * &
                                       (2.0_dp * this%eigen_WW(nn))
      end do

      ! ---- equiv-B (Model A) force contraction: d F = Tr[ W d_kchi h ] ----
      ! Same sparse-block structure as Hellmann-Feynman, with the filtered
      ! matrix W in place of u_n u_n^T (reuses the dh_block precompute):
      !   Tr[W dh] = 2 sum_{j/=k} sum_{ab} W(Qk_a, Qj_b) dh_block(a,b,chi,j)
      !            + sum_a W(Qk_a, Qk_a) * d_kchi(4 g_k^2)
      if (tb_model_equivb) then
        do iq = 1, 2*power_equivb
          tmp_sum(:) = 0.0_dp
          do jan = 1, this%dim_nn
            if (jan == kan) cycle
            off_j = 4*(jan-1)
            do ialp = 0, 3
              do ibet = 0, 3
                tmp_sum(:) = tmp_sum(:) + &
                  this%WW_equivb(off_k+1+ialp, off_j+1+ibet, iq) * dh_block(ialp, ibet, :, jan)
              end do
            end do
          end do
          tmp_sum(:) = 2.0_dp * tmp_sum(:)
          ! On-site block derivative: d_kchi(4 g_k^2) = 8 g_k dg_k cosk on the diagonal
          wdiag = 0.0_dp
          do ialp = 0, 3
            wdiag = wdiag + this%WW_equivb(off_k+1+ialp, off_k+1+ialp, iq)
          end do
          tmp_sum(:) = tmp_sum(:) + wdiag * 8.0_dp * gk * dgk * cosk(:)
          if (mod(iq, 2) == 0) then
            ! p sector: apply the 1/sqrt(3) Model A normalisation
            this%force_equivb_desc(:, iq, kan) = tmp_sum(:) * inv_sqrt3
          else
            ! s sector (norm 1)
            this%force_equivb_desc(:, iq, kan) = tmp_sum(:)
          end if
        end do
      end if
    end do  ! kan

    ! Pack into force_desc_matrix with truncation / zero-padding
    this%force_desc_matrix(:,:,:) = 0.0_dp
    if (this%tb_nn_G <= this%dim_hh) then
      this%force_desc_matrix(1:3, 1:this%tb_nn_G, 1:this%dim_nn) = &
        this%deriv_eigen(1:3, 1:this%tb_nn_G, 1:this%dim_nn)
    else
      this%force_desc_matrix(1:3, 1:this%dim_hh, 1:this%dim_nn) = &
        this%deriv_eigen(1:3, 1:this%dim_hh, 1:this%dim_nn)
    end if

    deallocate(g_arr, dg_arr, dh_block)
  end subroutine sp_evaluate_force

  ! ---------------------------------------------------------------------------
  ! Trace-based descriptors (same pattern as ss/pp; works on all eigenvalues)
  ! ---------------------------------------------------------------------------
  subroutine sp_compute_trace_energy(this, power_trace)
    class(sp_hamiltonian), intent(inout) :: this
    integer, intent(in) :: power_trace
    integer :: pp, ii, M

    if (allocated(this%trace_desc)) deallocate(this%trace_desc)
    allocate(this%trace_desc(power_trace))
    M = this%dim_hh
    do pp = 1, power_trace
      this%trace_desc(pp) = 0.0_dp
      do ii = 1, M
        this%trace_desc(pp) = this%trace_desc(pp) + this%eigen_WW(ii)**pp
      end do
    end do
  end subroutine sp_compute_trace_energy

  subroutine sp_compute_trace_force(this, power_trace)
    class(sp_hamiltonian), intent(inout) :: this
    integer, intent(in) :: power_trace
    integer :: pp, ii, kan, M
    real(dp) :: lam, dlam_factor

    if (allocated(this%force_trace_desc)) deallocate(this%force_trace_desc)
    allocate(this%force_trace_desc(3, power_trace, this%dim_nn))
    M = this%dim_hh
    this%force_trace_desc(:,:,:) = 0.0_dp
    do kan = 1, this%dim_nn
      do pp = 1, power_trace
        do ii = 1, M
          lam = this%eigen_WW(ii)
          if (abs(lam) > 1.0e-30_dp) then
            dlam_factor = real(pp, dp) * lam**(pp-2) * 0.5_dp
            this%force_trace_desc(:, pp, kan) = this%force_trace_desc(:, pp, kan) + &
              dlam_factor * this%deriv_eigen(:, ii, kan)
          end if
        end do
      end do
    end do
  end subroutine sp_compute_trace_force

  ! ---------------------------------------------------------------------------
  ! equiv-B (Model A) energy descriptors for the sp channel
  !
  ! Sector-projected polynomial filters (docs/Tbind_desc/main.tex, Sec. 6.4):
  !   F_s(q) = Tr[Pi_s h^q]           = sum_n lambda_n^q w_s(n)
  !   F_p(q) = Tr[Pi_p h^q]/sqrt(3)   = (Tr(h^q) - F_s(q))/sqrt(3)
  ! with the s-sector mode weight w_s(n) = sum_k u_{n}(4(k-1)+1)^2
  ! (= u_n^T Pi_s u_n, the S_{a,n} diagonal ss entry summed over sites).
  ! Uses the eigen-decomposition already computed by evaluate_energy.
  ! ---------------------------------------------------------------------------
  subroutine sp_compute_equivb_energy(this, power)
    class(sp_hamiltonian), intent(inout) :: this
    integer, intent(in) :: power

    integer :: qq, nn, kan, M
    real(dp) :: feat_s, m_q
    real(dp), dimension(:), allocatable :: w_s, lam_pow

    if (allocated(this%equivb_desc)) deallocate(this%equivb_desc)
    allocate(this%equivb_desc(2*power))
    M = this%dim_hh

    allocate(w_s(M), lam_pow(M))
    do nn = 1, M
      w_s(nn) = 0.0_dp
      do kan = 1, this%dim_nn
        w_s(nn) = w_s(nn) + this%u_matrix(4*(kan-1)+1, nn)**2
      end do
    end do

    lam_pow(:) = 1.0_dp
    do qq = 1, power
      lam_pow(1:M) = lam_pow(1:M) * this%eigen_WW(1:M)   ! lambda^qq
      feat_s = 0.0_dp
      m_q    = 0.0_dp
      do nn = 1, M
        feat_s = feat_s + lam_pow(nn) * w_s(nn)
        m_q    = m_q    + lam_pow(nn)
      end do
      this%equivb_desc(2*(qq-1)+1) = feat_s
      this%equivb_desc(2*(qq-1)+2) = (m_q - feat_s) * inv_sqrt3
    end do

    deallocate(w_s, lam_pow)
  end subroutine sp_compute_equivb_energy

  ! ---------------------------------------------------------------------------
  ! Build the equiv-B filtered matrices W for the sp channel.
  !
  ! The force of a sector-projected polynomial feature is the exact Frechet
  ! derivative (no eigenvector derivatives, smooth at degeneracies; see
  ! docs/Tbind_desc/main.tex App. C, Daleckii-Krein):
  !
  !   d Tr[Pi h^q] = Tr[ W_q dh ],
  !   W_q = U C_q U^T,  (C_q)_nm = phi_q(lambda_n, lambda_m) * (U^T Pi U)_nm,
  !   phi_q(x, y) = sum_{s=0}^{q-1} x^s y^{q-1-s}  ( = (x^q - y^q)/(x - y) ).
  !
  ! For the complement sector Pi_p = I - Pi_s:
  !   C_q^p = diag(q lambda_n^{q-1}) - C_q^s.
  ! Storage: WW_equivb(:,:,(q-1)*2+1) = s sector, (q-1)*2+2 = p sector (unscaled;
  ! the 1/sqrt(3) normalisation is applied at contraction time).
  ! ---------------------------------------------------------------------------
  subroutine sp_build_equivb_W(this, power)
    class(sp_hamiltonian), intent(inout) :: this
    integer, intent(in) :: power

    integer :: M, qq, nn, mm, kan
    real(dp), dimension(:,:), allocatable :: Us, Ptil, CC, CCp, tmpMM, phi
    real(dp), dimension(:), allocatable :: lam_pow

    M = this%dim_hh
    if (allocated(this%WW_equivb)) deallocate(this%WW_equivb)
    allocate(this%WW_equivb(M, M, 2*power))

    allocate(Us(this%dim_nn, M), Ptil(M, M), CC(M, M), CCp(M, M), tmpMM(M, M), phi(M, M))
    allocate(lam_pow(M))

    ! s-orbital rows of U: Us(kan, n) = u_matrix(4*(kan-1)+1, n)
    do kan = 1, this%dim_nn
      Us(kan, 1:M) = this%u_matrix(4*(kan-1)+1, 1:M)
    end do
    ! Ptil = U^T Pi_s U (Pi_s = diagonal 0/1 projector on the s rows)
    Ptil = matmul(transpose(Us), Us)

    ! phi_1(n,m) = 1; recursion phi_{q+1}(n,m) = lambda_n * phi_q(n,m) + lambda_m^q
    phi(:,:) = 1.0_dp
    lam_pow(:) = 1.0_dp          ! lambda^{q-1} for q = 1

    do qq = 1, power
      ! s sector: C = phi o Ptil (elementwise)
      CC(:,:) = phi(:,:) * Ptil(:,:)
      ! p sector: C_p = diag(phi(n,n)) - C, with phi(n,n) = q * lambda^{q-1}
      CCp(:,:) = -CC(:,:)
      do nn = 1, M
        CCp(nn, nn) = CCp(nn, nn) + dble(qq) * lam_pow(nn)
      end do

      tmpMM = matmul(CC, transpose(this%u_matrix))
      this%WW_equivb(:, :, 2*(qq-1)+1) = matmul(this%u_matrix, tmpMM)
      tmpMM = matmul(CCp, transpose(this%u_matrix))
      this%WW_equivb(:, :, 2*(qq-1)+2) = matmul(this%u_matrix, tmpMM)

      if (qq < power) then
        lam_pow(1:M) = lam_pow(1:M) * this%eigen_WW(1:M)   ! lambda^qq
        do mm = 1, M
          phi(1:M, mm) = this%eigen_WW(1:M) * phi(1:M, mm) + lam_pow(mm)
        end do
      end if
    end do

    deallocate(Us, Ptil, CC, CCp, tmpMM, phi, lam_pow)
  end subroutine sp_build_equivb_W

  ! ---------------------------------------------------------------------------
  ! Generic equiv-B energy for arbitrary spectral filters f_k:
  !   D_{k,s} = sum_n f_k(lambda_n) w_s(n),
  !   D_{k,p} = ( sum_n f_k(lambda_n) - D_{k,s} ) / sqrt(3),
  ! with the s-sector mode weight w_s(n) = sum_kan u(4(kan-1)+1, n)^2.
  ! Layout identical to the polynomial path: (k-1)*2+1 = s, +2 = p.
  ! ---------------------------------------------------------------------------
  subroutine sp_compute_filter_energy(this, nf)
    use module_tb_filter, only: tbf_eval
    use module_tbind, only: TB_CH_SP
    class(sp_hamiltonian), intent(inout) :: this
    integer, intent(in) :: nf

    integer :: nn, kan, k, M
    real(dp) :: ws
    real(dp), dimension(nf) :: fv, fp, acc_s, acc_t

    if (allocated(this%equivb_desc)) deallocate(this%equivb_desc)
    allocate(this%equivb_desc(2*nf))
    M = this%dim_hh

    acc_s(:) = 0.0_dp
    acc_t(:) = 0.0_dp
    do nn = 1, M
      ws = 0.0_dp
      do kan = 1, this%dim_nn
        ws = ws + this%u_matrix(4*(kan-1)+1, nn)**2
      end do
      call tbf_eval(TB_CH_SP, this%eigen_WW(nn), fv, fp)
      acc_s(1:nf) = acc_s(1:nf) + fv(1:nf)*ws
      acc_t(1:nf) = acc_t(1:nf) + fv(1:nf)
    end do
    do k = 1, nf
      this%equivb_desc(2*(k-1)+1) = acc_s(k)
      this%equivb_desc(2*(k-1)+2) = (acc_t(k) - acc_s(k))*inv_sqrt3
    end do
  end subroutine sp_compute_filter_energy

  ! ---------------------------------------------------------------------------
  ! Generic equiv-B filtered W matrices (arbitrary spectral filters):
  !   W_k = U C_k U^T,  (C_k)_nm = f_k[lambda_n, lambda_m] * Ptil_nm,
  ! with the divided difference
  !   f_k[x,y] = (f_k(x)-f_k(y))/(x-y)  for |x-y| > tol,
  !            = ( f_k'(x)+f_k'(y) )/2  otherwise (stable at degeneracies).
  ! Complement (p) sector: C_k^p = diag(f_k'(lambda_n)) - C_k^s.
  ! Storage identical to the polynomial builder: (k-1)*2+1 = s, +2 = p.
  ! ---------------------------------------------------------------------------
  subroutine sp_build_filter_W(this, nf)
    use module_tb_filter, only: tbf_eval
    use module_tbind, only: TB_CH_SP
    class(sp_hamiltonian), intent(inout) :: this
    integer, intent(in) :: nf

    integer :: M, k, nn, mm, kan
    real(dp) :: dlam, tol, dd
    real(dp), dimension(:,:), allocatable :: Us, Ptil, CC, CCp, tmpMM
    real(dp), dimension(:,:), allocatable :: fv_tab, fp_tab

    M = this%dim_hh
    if (allocated(this%WW_equivb)) deallocate(this%WW_equivb)
    allocate(this%WW_equivb(M, M, 2*nf))

    allocate(Us(this%dim_nn, M), Ptil(M, M), CC(M, M), CCp(M, M), tmpMM(M, M))
    allocate(fv_tab(nf, M), fp_tab(nf, M))

    ! s-orbital rows of U and Ptil = U^T Pi_s U (as in the polynomial builder)
    do kan = 1, this%dim_nn
      Us(kan, 1:M) = this%u_matrix(4*(kan-1)+1, 1:M)
    end do
    Ptil = matmul(transpose(Us), Us)

    ! filter values and derivatives at all eigenvalues
    do nn = 1, M
      call tbf_eval(TB_CH_SP, this%eigen_WW(nn), fv_tab(:, nn), fp_tab(:, nn))
    end do

    do k = 1, nf
      ! divided-difference matrix times Ptil (elementwise)
      do mm = 1, M
        do nn = 1, M
          dlam = this%eigen_WW(nn) - this%eigen_WW(mm)
          tol  = 1.0e-7_dp*(1.0_dp + abs(this%eigen_WW(nn)) + abs(this%eigen_WW(mm)))
          if (abs(dlam) > tol) then
            dd = (fv_tab(k, nn) - fv_tab(k, mm))/dlam
          else
            dd = 0.5_dp*(fp_tab(k, nn) + fp_tab(k, mm))
          end if
          CC(nn, mm) = dd*Ptil(nn, mm)
          CCp(nn, mm) = -CC(nn, mm)
        end do
      end do
      do nn = 1, M
        CCp(nn, nn) = CCp(nn, nn) + fp_tab(k, nn)
      end do

      tmpMM = matmul(CC, transpose(this%u_matrix))
      this%WW_equivb(:, :, 2*(k-1)+1) = matmul(this%u_matrix, tmpMM)
      tmpMM = matmul(CCp, transpose(this%u_matrix))
      this%WW_equivb(:, :, 2*(k-1)+2) = matmul(this%u_matrix, tmpMM)
    end do

    deallocate(Us, Ptil, CC, CCp, tmpMM, fv_tab, fp_tab)
  end subroutine sp_build_filter_W

end module module_sp_hamiltonian


! ============================================================================
! d-d coupling Hamiltonian (Slater-Koster, 5x5 blocks per neighbour pair)
!
! Basis ordering (1..5): dxy, dyz, dzx, dx2-y2, d3z2-r2.
! Diagonal (on-site) block: 4 g_d(r_ak)^2 * I_5  (eq. 70 of Tbind_desc.pdf).
! Off-site block h_{kj}^{(d-d)} = A_s T_s + A_p T_p + A_d T_d, with
!   A_mu = f_ddmu(r_kj) * g_d(r_ak) * g_d(r_aj),  mu in {sigma, pi, delta},
! and T_mu are 5x5 SK angular matrices in the direction cosines (l,m,n)=chikj.
! ============================================================================
module module_dd_hamiltonian
  use iso_fortran_env, only: dp => real64
  use module_eigenvalue_small_matrix, only: eigen_symmetric
  implicit none

  private

  type, public :: dd_hamiltonian
    integer :: dim_hh       ! 5 * max_neigh_local
    integer :: dim_nn       ! max_neigh_local
    integer :: tb_nn_G
    real(dp), dimension(:,:), allocatable :: hh
    real(dp), dimension(:), allocatable :: desc_matrix
    real(dp), dimension(:,:,:), allocatable :: force_desc_matrix
    real(dp), dimension(:), allocatable :: trace_desc
    real(dp), dimension(:,:,:), allocatable :: force_trace_desc
    ! Generic equiv-B spectral-filter descriptors (non-polynomial filters)
    real(dp), dimension(:), allocatable :: equivb_desc            ! (nf)
    real(dp), dimension(:,:,:), allocatable :: force_equivb_desc  ! (3, nf, dim_nn)
    real(dp), dimension(:,:,:), allocatable :: deriv_eigen
    real(dp), dimension(:,:,:), allocatable :: deriv_eigen_raw  ! raw d(lambda)/dx
    real(dp), dimension(:,:), allocatable :: u_matrix
    real(dp), dimension(:), allocatable :: eigen_WW
    type(eigen_symmetric) :: eigsolver
    logical :: sinit = .false.
  contains
    procedure :: init                 => dd_init
    procedure :: evaluate_energy      => dd_evaluate_energy
    procedure :: evaluate_force       => dd_evaluate_force
    procedure :: compute_trace_energy => dd_compute_trace_energy
    procedure :: compute_trace_force  => dd_compute_trace_force
    procedure :: compute_filter_energy => dd_compute_filter_energy
    procedure :: compute_filter_force  => dd_compute_filter_force
  end type dd_hamiltonian

  type(dd_hamiltonian), public :: tb_dd

contains

  subroutine dd_init(this, tb_svd, tb_nn_G, max_neigh_local)
    class(dd_hamiltonian), intent(inout) :: this
    integer, intent(in) :: tb_svd, tb_nn_G, max_neigh_local

    if (this%sinit .and. this%dim_nn == max_neigh_local .and. this%tb_nn_G == tb_nn_G) return

    this%dim_hh  = 5 * max_neigh_local
    this%dim_nn  = max_neigh_local
    this%tb_nn_G = tb_nn_G

    if (allocated(this%hh)) deallocate(this%hh)
    allocate(this%hh(this%dim_hh, this%dim_hh))

    if (allocated(this%deriv_eigen)) deallocate(this%deriv_eigen)
    allocate(this%deriv_eigen(3, this%dim_hh, this%dim_nn))
    if (allocated(this%deriv_eigen_raw)) deallocate(this%deriv_eigen_raw)
    allocate(this%deriv_eigen_raw(3, this%dim_hh, this%dim_nn))
    if (allocated(this%u_matrix)) deallocate(this%u_matrix)
    allocate(this%u_matrix(this%dim_hh, this%dim_hh))
    if (allocated(this%eigen_WW)) deallocate(this%eigen_WW)
    allocate(this%eigen_WW(this%dim_hh))

    if (allocated(this%desc_matrix)) deallocate(this%desc_matrix)
    allocate(this%desc_matrix(this%tb_nn_G))
    if (allocated(this%force_desc_matrix)) deallocate(this%force_desc_matrix)
    allocate(this%force_desc_matrix(3, this%tb_nn_G, this%dim_nn))

    call this%eigsolver%init(tb_svd, this%dim_hh)
    this%sinit = .true.
  end subroutine dd_init

  ! ---------------------------------------------------------------------------
  ! Fill the 5x5 angular matrices T_sigma, T_pi, T_delta (upper triangle).
  ! Lower triangle filled by symmetry (T is symmetric in the in-pair orbital
  ! exchange because the block h^{d-d}_{kj} is symmetric under (alpha<->beta)).
  ! Direction cosines (l,m,n) are NOT assumed normalised individually; the
  ! formulas use them directly.
  ! ---------------------------------------------------------------------------
  pure subroutine dd_T_matrices(l, m, n, Ts, Tp, Td)
    real(dp), intent(in)  :: l, m, n
    real(dp), dimension(5,5), intent(out) :: Ts, Tp, Td

    real(dp), parameter :: s3 = 1.7320508075688772_dp
    real(dp) :: ll, mm, nn, v, u, w, lm, mn, ln, lmn

    ll = l*l;  mm = m*m;  nn = n*n
    v  = ll - mm
    u  = nn - 0.5_dp*(ll + mm)
    w  = ll + mm - nn
    lm = l*m;  mn = m*n;  ln = l*n;  lmn = l*m*n

    ! ----- T_sigma -----
    Ts(1,1) = 3.0_dp*ll*mm
    Ts(1,2) = 3.0_dp*l*mm*n
    Ts(1,3) = 3.0_dp*ll*mn
    Ts(1,4) = 1.5_dp*lm*v
    Ts(1,5) = s3*lm*u
    Ts(2,2) = 3.0_dp*mm*nn
    Ts(2,3) = 3.0_dp*lm*nn
    Ts(2,4) = 1.5_dp*mn*v
    Ts(2,5) = s3*mn*u
    Ts(3,3) = 3.0_dp*ll*nn
    Ts(3,4) = 1.5_dp*ln*v
    Ts(3,5) = s3*ln*u
    Ts(4,4) = 0.75_dp*v*v
    Ts(4,5) = 0.5_dp*s3*v*u
    Ts(5,5) = u*u

    ! ----- T_pi -----
    Tp(1,1) = ll + mm - 4.0_dp*ll*mm
    Tp(1,2) = ln*(1.0_dp - 4.0_dp*mm)
    Tp(1,3) = mn*(1.0_dp - 4.0_dp*ll)
    Tp(1,4) = 2.0_dp*lm*(mm - ll)
    Tp(1,5) = -2.0_dp*s3*lm*nn
    Tp(2,2) = mm + nn - 4.0_dp*mm*nn
    Tp(2,3) = lm*(1.0_dp - 4.0_dp*nn)
    Tp(2,4) = -mn*(1.0_dp + 2.0_dp*v)
    Tp(2,5) = s3*mn*w
    Tp(3,3) = ll + nn - 4.0_dp*ll*nn
    Tp(3,4) = ln*(1.0_dp - 2.0_dp*v)
    Tp(3,5) = s3*ln*w
    Tp(4,4) = ll + mm - v*v
    Tp(4,5) = -s3*nn*v
    Tp(5,5) = 3.0_dp*nn*(ll + mm)

    ! ----- T_delta -----
    Td(1,1) = nn + ll*mm
    Td(1,2) = ln*(mm - 1.0_dp)
    Td(1,3) = mn*(ll - 1.0_dp)
    Td(1,4) = 0.5_dp*lm*v
    Td(1,5) = 0.5_dp*s3*lm*(1.0_dp + nn)
    Td(2,2) = ll + mm*nn
    Td(2,3) = lm*(nn - 1.0_dp)
    Td(2,4) = mn*(1.0_dp + 0.5_dp*v)
    Td(2,5) = -0.5_dp*s3*mn*(ll + mm)
    Td(3,3) = mm + ll*nn
    Td(3,4) = -ln*(1.0_dp - 0.5_dp*v)
    Td(3,5) = -0.5_dp*s3*ln*(ll + mm)
    Td(4,4) = nn + 0.25_dp*v*v
    Td(4,5) = 0.25_dp*s3*(1.0_dp + nn)*v
    Td(5,5) = 0.75_dp*(ll + mm)**2

    call dd_symmetrize(Ts)
    call dd_symmetrize(Tp)
    call dd_symmetrize(Td)
  end subroutine dd_T_matrices

  pure subroutine dd_symmetrize(M)
    real(dp), dimension(5,5), intent(inout) :: M
    integer :: i, j
    do i = 1, 5
      do j = i+1, 5
        M(j,i) = M(i,j)
      end do
    end do
  end subroutine dd_symmetrize

  ! ---------------------------------------------------------------------------
  ! Derivatives of T_mu w.r.t. l, m, n (upper triangle), mu in {s,p,d}.
  ! ---------------------------------------------------------------------------
  pure subroutine dd_dT_matrices(l, m, n, &
                                  dTs_l, dTs_m, dTs_n, &
                                  dTp_l, dTp_m, dTp_n, &
                                  dTd_l, dTd_m, dTd_n)
    real(dp), intent(in) :: l, m, n
    real(dp), dimension(5,5), intent(out) :: dTs_l, dTs_m, dTs_n
    real(dp), dimension(5,5), intent(out) :: dTp_l, dTp_m, dTp_n
    real(dp), dimension(5,5), intent(out) :: dTd_l, dTd_m, dTd_n

    real(dp), parameter :: s3 = 1.7320508075688772_dp
    real(dp) :: ll, mm, nn, v, u, w, lm, mn, ln, lmn

    ll = l*l;  mm = m*m;  nn = n*n
    v  = ll - mm
    u  = nn - 0.5_dp*(ll + mm)
    w  = ll + mm - nn
    lm = l*m;  mn = m*n;  ln = l*n;  lmn = l*m*n

    ! ------ d T_sigma / d{l,m,n} ------
    dTs_l(1,1) = 6.0_dp*l*mm;        dTs_m(1,1) = 6.0_dp*ll*m;          dTs_n(1,1) = 0.0_dp
    dTs_l(1,2) = 3.0_dp*mm*n;        dTs_m(1,2) = 6.0_dp*lmn;            dTs_n(1,2) = 3.0_dp*l*mm
    dTs_l(1,3) = 6.0_dp*lmn;          dTs_m(1,3) = 3.0_dp*ll*n;          dTs_n(1,3) = 3.0_dp*ll*m
    dTs_l(1,4) = 1.5_dp*m*(3.0_dp*ll - mm); dTs_m(1,4) = 1.5_dp*l*(ll - 3.0_dp*mm); dTs_n(1,4) = 0.0_dp
    dTs_l(1,5) = s3*m*(u - ll);      dTs_m(1,5) = s3*l*(u - mm);        dTs_n(1,5) = 2.0_dp*s3*lmn
    dTs_l(2,2) = 0.0_dp;             dTs_m(2,2) = 6.0_dp*m*nn;          dTs_n(2,2) = 6.0_dp*mm*n
    dTs_l(2,3) = 3.0_dp*m*nn;        dTs_m(2,3) = 3.0_dp*l*nn;          dTs_n(2,3) = 6.0_dp*lmn
    dTs_l(2,4) = 3.0_dp*lmn;          dTs_m(2,4) = 1.5_dp*n*(ll - 3.0_dp*mm); dTs_n(2,4) = 1.5_dp*m*v
    dTs_l(2,5) = -s3*lmn;             dTs_m(2,5) = s3*n*(u - mm);        dTs_n(2,5) = s3*m*(u + 2.0_dp*nn)
    dTs_l(3,3) = 6.0_dp*l*nn;        dTs_m(3,3) = 0.0_dp;               dTs_n(3,3) = 6.0_dp*ll*n
    dTs_l(3,4) = 1.5_dp*n*(3.0_dp*ll - mm); dTs_m(3,4) = -3.0_dp*lmn;     dTs_n(3,4) = 1.5_dp*l*v
    dTs_l(3,5) = s3*n*(u - ll);      dTs_m(3,5) = -s3*lmn;               dTs_n(3,5) = s3*l*(u + 2.0_dp*nn)
    dTs_l(4,4) = 3.0_dp*l*v;         dTs_m(4,4) = -3.0_dp*m*v;          dTs_n(4,4) = 0.0_dp
    dTs_l(4,5) = 0.5_dp*s3*l*(2.0_dp*u - v); dTs_m(4,5) = -0.5_dp*s3*m*(2.0_dp*u + v); dTs_n(4,5) = s3*n*v
    dTs_l(5,5) = -2.0_dp*l*u;        dTs_m(5,5) = -2.0_dp*m*u;          dTs_n(5,5) = 4.0_dp*n*u

    ! ------ d T_pi / d{l,m,n} ------
    dTp_l(1,1) = 2.0_dp*l*(1.0_dp - 4.0_dp*mm); dTp_m(1,1) = 2.0_dp*m*(1.0_dp - 4.0_dp*ll); dTp_n(1,1) = 0.0_dp
    dTp_l(1,2) = n*(1.0_dp - 4.0_dp*mm);  dTp_m(1,2) = -8.0_dp*lmn; dTp_n(1,2) = l*(1.0_dp - 4.0_dp*mm)
    dTp_l(1,3) = -8.0_dp*lmn;              dTp_m(1,3) = n*(1.0_dp - 4.0_dp*ll); dTp_n(1,3) = m*(1.0_dp - 4.0_dp*ll)
    dTp_l(1,4) = 2.0_dp*m*(mm - 3.0_dp*ll); dTp_m(1,4) = 2.0_dp*l*(3.0_dp*mm - ll); dTp_n(1,4) = 0.0_dp
    dTp_l(1,5) = -2.0_dp*s3*m*nn;          dTp_m(1,5) = -2.0_dp*s3*l*nn; dTp_n(1,5) = -4.0_dp*s3*lmn
    dTp_l(2,2) = 0.0_dp; dTp_m(2,2) = 2.0_dp*m*(1.0_dp - 4.0_dp*nn); dTp_n(2,2) = 2.0_dp*n*(1.0_dp - 4.0_dp*mm)
    dTp_l(2,3) = m*(1.0_dp - 4.0_dp*nn); dTp_m(2,3) = l*(1.0_dp - 4.0_dp*nn); dTp_n(2,3) = -8.0_dp*lmn
    dTp_l(2,4) = -4.0_dp*lmn;              dTp_m(2,4) = -n*(1.0_dp + 2.0_dp*ll - 6.0_dp*mm); dTp_n(2,4) = -m*(1.0_dp + 2.0_dp*v)
    dTp_l(2,5) = 2.0_dp*s3*lmn;            dTp_m(2,5) = s3*n*(w + 2.0_dp*mm); dTp_n(2,5) = s3*m*(w - 2.0_dp*nn)
    dTp_l(3,3) = 2.0_dp*l*(1.0_dp - 4.0_dp*nn); dTp_m(3,3) = 0.0_dp; dTp_n(3,3) = 2.0_dp*n*(1.0_dp - 4.0_dp*ll)
    dTp_l(3,4) = n*(1.0_dp - 6.0_dp*ll + 2.0_dp*mm); dTp_m(3,4) = 4.0_dp*lmn; dTp_n(3,4) = l*(1.0_dp - 2.0_dp*v)
    dTp_l(3,5) = s3*n*(w + 2.0_dp*ll); dTp_m(3,5) = 2.0_dp*s3*lmn; dTp_n(3,5) = s3*l*(w - 2.0_dp*nn)
    dTp_l(4,4) = 2.0_dp*l*(1.0_dp - 2.0_dp*v); dTp_m(4,4) = 2.0_dp*m*(1.0_dp + 2.0_dp*v); dTp_n(4,4) = 0.0_dp
    dTp_l(4,5) = -2.0_dp*s3*l*nn;          dTp_m(4,5) = 2.0_dp*s3*m*nn; dTp_n(4,5) = -2.0_dp*s3*n*v
    dTp_l(5,5) = 6.0_dp*l*nn;              dTp_m(5,5) = 6.0_dp*m*nn; dTp_n(5,5) = 6.0_dp*n*(ll + mm)

    ! ------ d T_delta / d{l,m,n} ------
    dTd_l(1,1) = 2.0_dp*l*mm;        dTd_m(1,1) = 2.0_dp*ll*m;        dTd_n(1,1) = 2.0_dp*n
    dTd_l(1,2) = n*(mm - 1.0_dp);    dTd_m(1,2) = 2.0_dp*lmn;          dTd_n(1,2) = l*(mm - 1.0_dp)
    dTd_l(1,3) = 2.0_dp*lmn;          dTd_m(1,3) = n*(ll - 1.0_dp);    dTd_n(1,3) = m*(ll - 1.0_dp)
    dTd_l(1,4) = 0.5_dp*m*(3.0_dp*ll - mm); dTd_m(1,4) = 0.5_dp*l*(ll - 3.0_dp*mm); dTd_n(1,4) = 0.0_dp
    dTd_l(1,5) = 0.5_dp*s3*m*(1.0_dp + nn); dTd_m(1,5) = 0.5_dp*s3*l*(1.0_dp + nn); dTd_n(1,5) = s3*lmn
    dTd_l(2,2) = 2.0_dp*l;            dTd_m(2,2) = 2.0_dp*m*nn;        dTd_n(2,2) = 2.0_dp*mm*n
    dTd_l(2,3) = m*(nn - 1.0_dp);    dTd_m(2,3) = l*(nn - 1.0_dp);    dTd_n(2,3) = 2.0_dp*lmn
    dTd_l(2,4) = lmn;                 dTd_m(2,4) = n*(1.0_dp + 0.5_dp*ll - 1.5_dp*mm); dTd_n(2,4) = m*(1.0_dp + 0.5_dp*v)
    dTd_l(2,5) = -s3*lmn;             dTd_m(2,5) = -0.5_dp*s3*n*(ll + 3.0_dp*mm); dTd_n(2,5) = -0.5_dp*s3*m*(ll + mm)
    dTd_l(3,3) = 2.0_dp*l*nn;        dTd_m(3,3) = 2.0_dp*m;           dTd_n(3,3) = 2.0_dp*ll*n
    dTd_l(3,4) = 0.5_dp*n*(3.0_dp*ll - mm - 2.0_dp); dTd_m(3,4) = -lmn; dTd_n(3,4) = 0.5_dp*l*(ll - mm - 2.0_dp)
    dTd_l(3,5) = -0.5_dp*s3*n*(3.0_dp*ll + mm); dTd_m(3,5) = -s3*lmn; dTd_n(3,5) = -0.5_dp*s3*l*(ll + mm)
    dTd_l(4,4) = l*v;                dTd_m(4,4) = -m*v;               dTd_n(4,4) = 2.0_dp*n
    dTd_l(4,5) = 0.5_dp*s3*l*(1.0_dp + nn); dTd_m(4,5) = -0.5_dp*s3*m*(1.0_dp + nn); dTd_n(4,5) = 0.5_dp*s3*n*v
    dTd_l(5,5) = 3.0_dp*l*(ll + mm); dTd_m(5,5) = 3.0_dp*m*(ll + mm); dTd_n(5,5) = 0.0_dp

    call dd_symmetrize(dTs_l); call dd_symmetrize(dTs_m); call dd_symmetrize(dTs_n)
    call dd_symmetrize(dTp_l); call dd_symmetrize(dTp_m); call dd_symmetrize(dTp_n)
    call dd_symmetrize(dTd_l); call dd_symmetrize(dTd_m); call dd_symmetrize(dTd_n)
  end subroutine dd_dT_matrices

  ! ---------------------------------------------------------------------------
  ! Build h^{dd} (5 N_v x 5 N_v) and diagonalise
  ! ---------------------------------------------------------------------------
  subroutine dd_evaluate_energy(this, rad, mua_a, mua_neigh, kk, ll_sg, ll_pi, ll_dl, &
                                 r_central, r_matrix, x_matrix)
    use module_tb_radial, only: tb_channel_radial_t
    class(dd_hamiltonian), intent(inout) :: this
    type(tb_channel_radial_t), intent(in) :: rad
    integer, intent(in) :: mua_a
    integer, dimension(:), intent(in) :: mua_neigh
    integer, intent(in) :: kk, ll_sg, ll_pi, ll_dl
    real(dp), dimension(:), intent(in) :: r_central
    real(dp), dimension(0:,0:), intent(in) :: r_matrix
    real(dp), dimension(:,0:,0:), intent(in) :: x_matrix

    integer :: kan, jan, a, b, off_k, off_j
    real(dp) :: rr, gk, gj, gkgj, A_s, A_p, A_d
    real(dp) :: l, m, n
    real(dp), dimension(5,5) :: Ts, Tp, Td
    real(dp), dimension(:), allocatable :: g_arr

    allocate(g_arr(this%dim_nn))
    do kan = 1, this%dim_nn
      g_arr(kan) = rad%g_eval_kl(mua_a, mua_neigh(kan), kk, ll_sg, r_central(kan))
    end do

    this%hh(:,:) = 0.0_dp

    ! Diagonal (on-site) blocks: 4 g_k^2 * I_5
    do kan = 1, this%dim_nn
      gk = g_arr(kan)
      off_k = 5*(kan-1)
      do a = 1, 5
        this%hh(off_k+a, off_k+a) = 4.0_dp * gk*gk
      end do
    end do

    do kan = 1, this%dim_nn
      gk = g_arr(kan)
      off_k = 5*(kan-1)
      do jan = kan+1, this%dim_nn
        gj = g_arr(jan)
        off_j = 5*(jan-1)
        rr = r_matrix(kan, jan)
        l = x_matrix(1, kan, jan)
        m = x_matrix(2, kan, jan)
        n = x_matrix(3, kan, jan)
        gkgj = gk * gj

        A_s = rad%f_eval(mua_neigh(kan), mua_neigh(jan), kk, ll_sg, rr) * gkgj
        A_p = rad%f_eval(mua_neigh(kan), mua_neigh(jan), kk, ll_pi, rr) * gkgj
        A_d = rad%f_eval(mua_neigh(kan), mua_neigh(jan), kk, ll_dl, rr) * gkgj

        call dd_T_matrices(l, m, n, Ts, Tp, Td)

        do a = 1, 5
          do b = 1, 5
            this%hh(off_k+a, off_j+b) = A_s * Ts(a,b) + A_p * Tp(a,b) + A_d * Td(a,b)
            this%hh(off_j+b, off_k+a) = this%hh(off_k+a, off_j+b)
          end do
        end do
      end do
    end do

    deallocate(g_arr)

    call this%eigsolver%evaluate(this%hh)

    if (this%tb_nn_G <= this%dim_hh) then
      this%desc_matrix(1:this%tb_nn_G) = this%eigsolver%eigval(1:this%tb_nn_G)
    else
      this%desc_matrix(1:this%dim_hh) = this%eigsolver%eigval(1:this%dim_hh)
      this%desc_matrix(this%dim_hh+1:this%tb_nn_G) = 0.0_dp
    end if

    this%u_matrix(1:this%dim_hh, 1:this%dim_hh) = this%eigsolver%eigvec(1:this%dim_hh, 1:this%dim_hh)
    this%eigen_WW(1:this%dim_hh) = this%eigsolver%eigWWout(1:this%dim_hh)
  end subroutine dd_evaluate_energy

  ! ---------------------------------------------------------------------------
  ! Hellmann-Feynman force for d-d channel (eqs. 76-83 of Tbind_desc.pdf)
  ! Same sparse-block layout as sp: each neighbour k carries a 5-orbital block.
  ! ---------------------------------------------------------------------------
  subroutine dd_evaluate_force(this, rad, mua_a, mua_neigh, kk, ll_sg, ll_pi, ll_dl, &
                                r_central, r_matrix, x_matrix, tmp_dxp)
    use module_tb_radial, only: tb_channel_radial_t
    use module_tbind, only: tb_model_trace, tb_model_equivb, tb_equivb_rawderiv
    class(dd_hamiltonian), intent(inout) :: this
    type(tb_channel_radial_t), intent(in) :: rad
    integer, intent(in) :: mua_a
    integer, dimension(:), intent(in) :: mua_neigh
    integer, intent(in) :: kk, ll_sg, ll_pi, ll_dl
    real(dp), dimension(:), intent(in) :: r_central
    real(dp), dimension(0:,0:), intent(in) :: r_matrix
    real(dp), dimension(:,0:,0:), intent(in) :: x_matrix
    real(dp), dimension(:,0:), intent(in) :: tmp_dxp

    real(dp), dimension(:), allocatable :: g_arr, dg_arr
    real(dp), dimension(:,:,:,:), allocatable :: dh_block  ! (5,5,3,dim_nn)
    real(dp), dimension(5,5) :: Ts, Tp, Td
    real(dp), dimension(5,5) :: dTs_l, dTs_m, dTs_n
    real(dp), dimension(5,5) :: dTp_l, dTp_m, dTp_n
    real(dp), dimension(5,5) :: dTd_l, dTd_m, dTd_n
    real(dp), dimension(3,3) :: dlmn  ! d{l,m,n}(a) / dk_chi
    real(dp), dimension(3) :: cosk, chikj
    real(dp) :: rk, rr, gk, gj, dgk, inv_rkj, l, m, n
    real(dp) :: f_s, df_s, f_p, df_p, f_d, df_d
    real(dp) :: A_s, A_p, A_d
    real(dp), dimension(3) :: dA_s, dA_p, dA_d
    real(dp) :: sum_u2, u_alpha, u_beta, dTl_ab, dTm_ab, dTn_ab, T_combined
    real(dp), dimension(3) :: diag_part, tmp_sum, dh_abx
    integer :: kan, jan, nn, a, b, ialp, ibet, off_k, off_j, nn_max_deriv, chi

    allocate(g_arr(this%dim_nn), dg_arr(this%dim_nn))
    do kan = 1, this%dim_nn
      g_arr(kan)  = rad%g_eval_kl(mua_a, mua_neigh(kan), kk, ll_sg, r_central(kan))
      dg_arr(kan) = rad%g_deval_kl(mua_a, mua_neigh(kan), kk, ll_sg, r_central(kan))
    end do

    this%deriv_eigen(:,:,:) = 0.0_dp

    ! Trace or equiv-B (scaled trace on this single-sector channel): all eigenvalues
    if (tb_model_trace .or. tb_model_equivb) then
      nn_max_deriv = this%dim_hh
    else
      nn_max_deriv = min(this%tb_nn_G, this%dim_hh)
    end if

    allocate(dh_block(5, 5, 3, this%dim_nn))

    do kan = 1, this%dim_nn
      rk = r_central(kan)
      cosk(:) = tmp_dxp(:, kan) / rk
      gk  = g_arr(kan)
      dgk = dg_arr(kan)
      off_k = 5*(kan-1)

      dh_block(:,:,:,:) = 0.0_dp

      do jan = 1, this%dim_nn
        if (jan == kan) cycle
        gj = g_arr(jan)
        rr = r_matrix(kan, jan)
        chikj(1:3) = x_matrix(1:3, kan, jan)
        l = chikj(1);  m = chikj(2);  n = chikj(3)
        inv_rkj = 1.0_dp / rr

        f_s  = rad%f_eval (mua_neigh(kan), mua_neigh(jan), kk, ll_sg, rr)
        df_s = rad%f_deval(mua_neigh(kan), mua_neigh(jan), kk, ll_sg, rr)
        f_p  = rad%f_eval (mua_neigh(kan), mua_neigh(jan), kk, ll_pi, rr)
        df_p = rad%f_deval(mua_neigh(kan), mua_neigh(jan), kk, ll_pi, rr)
        f_d  = rad%f_eval (mua_neigh(kan), mua_neigh(jan), kk, ll_dl, rr)
        df_d = rad%f_deval(mua_neigh(kan), mua_neigh(jan), kk, ll_dl, rr)

        A_s = f_s * gk * gj
        A_p = f_p * gk * gj
        A_d = f_d * gk * gj
        ! d_kchi A_mu = (df_mu * gk * chikj(chi) + f_mu * dg_k * cosk(chi)) * gj
        dA_s(:) = (df_s * gk * chikj(:) + f_s * dgk * cosk(:)) * gj
        dA_p(:) = (df_p * gk * chikj(:) + f_p * dgk * cosk(:)) * gj
        dA_d(:) = (df_d * gk * chikj(:) + f_d * dgk * cosk(:)) * gj

        ! d_kchi (l,m,n)_a = delta_{a,chi}/r_kj - chikj(a)*chikj(chi)/r_kj
        do a = 1, 3
          do chi = 1, 3
            dlmn(a, chi) = -chikj(a) * chikj(chi) * inv_rkj
          end do
          dlmn(a, a) = dlmn(a, a) + inv_rkj
        end do

        call dd_T_matrices (l, m, n, Ts, Tp, Td)
        call dd_dT_matrices(l, m, n, dTs_l, dTs_m, dTs_n, &
                                       dTp_l, dTp_m, dTp_n, &
                                       dTd_l, dTd_m, dTd_n)

        do b = 1, 5
          do a = 1, 5
            dTl_ab = A_s*dTs_l(a,b) + A_p*dTp_l(a,b) + A_d*dTd_l(a,b)
            dTm_ab = A_s*dTs_m(a,b) + A_p*dTp_m(a,b) + A_d*dTd_m(a,b)
            dTn_ab = A_s*dTs_n(a,b) + A_p*dTp_n(a,b) + A_d*dTd_n(a,b)
            ! Radial-amplitude contribution
            dh_abx(:) = dA_s(:)*Ts(a,b) + dA_p(:)*Tp(a,b) + dA_d(:)*Td(a,b)
            ! Angular contribution via chain rule on (l,m,n)
            dh_abx(:) = dh_abx(:) + dTl_ab*dlmn(1,:) + dTm_ab*dlmn(2,:) + dTn_ab*dlmn(3,:)
            dh_block(a, b, :, jan) = dh_abx(:)
          end do
        end do
      end do  ! jan

      ! ---- Hellmann-Feynman contraction per eigenvalue ----
      do nn = 1, nn_max_deriv
        ! On-site Q_k x Q_k block: d_kchi (4 g_k^2) = 8 g_k dg_k cosk(chi), times I_5
        sum_u2 = 0.0_dp
        do ialp = 1, 5
          sum_u2 = sum_u2 + this%u_matrix(off_k+ialp, nn)**2
        end do
        diag_part(:) = sum_u2 * 8.0_dp * gk * dgk * cosk(:)

        tmp_sum(:) = diag_part(:)

        do jan = 1, this%dim_nn
          if (jan == kan) cycle
          off_j = 5*(jan-1)
          do ialp = 1, 5
            u_alpha = this%u_matrix(off_k+ialp, nn)
            do ibet = 1, 5
              u_beta = this%u_matrix(off_j+ibet, nn)
              tmp_sum(:) = tmp_sum(:) + u_alpha * dh_block(ialp, ibet, :, jan) * u_beta
            end do
          end do
        end do

        ! d lambda_n / d k_chi = 2*sum1 - sum2; descriptor is lambda^2 -> *2*lambda
        tmp_sum(:) = 2.0_dp * tmp_sum(:) - diag_part(:)
        if (tb_equivb_rawderiv) this%deriv_eigen_raw(:, nn, kan) = tmp_sum(:)
        this%deriv_eigen(:, nn, kan) = tmp_sum(:) * (2.0_dp * this%eigen_WW(nn))
      end do
    end do  ! kan

    this%force_desc_matrix(:,:,:) = 0.0_dp
    if (this%tb_nn_G <= this%dim_hh) then
      this%force_desc_matrix(1:3, 1:this%tb_nn_G, 1:this%dim_nn) = &
        this%deriv_eigen(1:3, 1:this%tb_nn_G, 1:this%dim_nn)
    else
      this%force_desc_matrix(1:3, 1:this%dim_hh, 1:this%dim_nn) = &
        this%deriv_eigen(1:3, 1:this%dim_hh, 1:this%dim_nn)
    end if

    deallocate(g_arr, dg_arr, dh_block)
  end subroutine dd_evaluate_force

  subroutine dd_compute_trace_energy(this, power_trace)
    class(dd_hamiltonian), intent(inout) :: this
    integer, intent(in) :: power_trace
    integer :: pp, ii, M

    if (allocated(this%trace_desc)) deallocate(this%trace_desc)
    allocate(this%trace_desc(power_trace))
    M = this%dim_hh
    do pp = 1, power_trace
      this%trace_desc(pp) = 0.0_dp
      do ii = 1, M
        this%trace_desc(pp) = this%trace_desc(pp) + this%eigen_WW(ii)**pp
      end do
    end do
  end subroutine dd_compute_trace_energy

  subroutine dd_compute_trace_force(this, power_trace)
    class(dd_hamiltonian), intent(inout) :: this
    integer, intent(in) :: power_trace
    integer :: pp, ii, kan, M
    real(dp) :: lam, dlam_factor

    if (allocated(this%force_trace_desc)) deallocate(this%force_trace_desc)
    allocate(this%force_trace_desc(3, power_trace, this%dim_nn))
    M = this%dim_hh
    this%force_trace_desc(:,:,:) = 0.0_dp
    do kan = 1, this%dim_nn
      do pp = 1, power_trace
        do ii = 1, M
          lam = this%eigen_WW(ii)
          if (abs(lam) > 1.0e-30_dp) then
            dlam_factor = real(pp, dp) * lam**(pp-2) * 0.5_dp
            this%force_trace_desc(:, pp, kan) = this%force_trace_desc(:, pp, kan) + &
              dlam_factor * this%deriv_eigen(:, ii, kan)
          end if
        end do
      end do
    end do
  end subroutine dd_compute_trace_force

  ! ---------------------------------------------------------------------------
  ! Generic equiv-B spectral-filter descriptors (single d sector; norm 1/sqrt5)
  ! ---------------------------------------------------------------------------
  subroutine dd_compute_filter_energy(this, nf, snorm)
    use module_tb_filter, only: tbf_eval
    use module_tbind, only: TB_CH_DD
    class(dd_hamiltonian), intent(inout) :: this
    integer, intent(in) :: nf
    real(dp), intent(in) :: snorm
    integer :: nn
    real(dp), dimension(nf) :: fv, fp

    if (allocated(this%equivb_desc)) deallocate(this%equivb_desc)
    allocate(this%equivb_desc(nf))
    this%equivb_desc(:) = 0.0_dp
    do nn = 1, this%dim_hh
      call tbf_eval(TB_CH_DD, this%eigen_WW(nn), fv, fp)
      this%equivb_desc(1:nf) = this%equivb_desc(1:nf) + fv(1:nf)
    end do
    this%equivb_desc(:) = this%equivb_desc(:)*snorm
  end subroutine dd_compute_filter_energy

  subroutine dd_compute_filter_force(this, nf, snorm)
    use module_tb_filter, only: tbf_eval
    use module_tbind, only: TB_CH_DD
    class(dd_hamiltonian), intent(inout) :: this
    integer, intent(in) :: nf
    real(dp), intent(in) :: snorm
    integer :: nn, kan, k, M
    real(dp), dimension(nf) :: fv
    real(dp), dimension(:,:), allocatable :: fp_tab

    M = this%dim_hh
    if (allocated(this%force_equivb_desc)) deallocate(this%force_equivb_desc)
    allocate(this%force_equivb_desc(3, nf, this%dim_nn))
    this%force_equivb_desc(:,:,:) = 0.0_dp

    allocate(fp_tab(nf, M))
    do nn = 1, M
      call tbf_eval(TB_CH_DD, this%eigen_WW(nn), fv, fp_tab(:, nn))
    end do
    do kan = 1, this%dim_nn
      do nn = 1, M
        do k = 1, nf
          this%force_equivb_desc(:, k, kan) = this%force_equivb_desc(:, k, kan) + &
            fp_tab(k, nn)*this%deriv_eigen_raw(:, nn, kan)
        end do
      end do
    end do
    this%force_equivb_desc(:,:,:) = this%force_equivb_desc(:,:,:)*snorm
    deallocate(fp_tab)
  end subroutine dd_compute_filter_force

end module module_dd_hamiltonian


! ============================================================================
! Main compute_tbind module and init_tbind subroutine
! ============================================================================
module module_compute_tbind

contains

subroutine compute_tbind(i_start_at, i_final_at, d_n_neigh, d_kind_neigh, iconf)
  USE module_kind_variables, ONLY: double, kind_double
#ifdef MLD_NDM
  use gen_com_m, ONLY: A2cm, lperiod
  use gen_com_m_ml, ONLY: imm, bg, at
  use tab_imm_m_ml, ONLY: xp
#else
  use ondm_gen_com_m, ONLY: imm, A2cm, lperiod, bg, at
  use ondm_tab_imm_m, ONLY: iwmax2, xp
#endif
  use ml_in_ndm_module, ONLY: imm_neigh, desc_forces, rangml, debug
  use module_neigh_local, only: r_cut, preallocate_neigh_ja, build_local_neighbours_ja
  use derived_types, only: config_real, config_desc
  use time_check_general, only: debug_time, MY_MPI_WTIME
  use module_tbind, only: tb_ham_ss, tb_ham_pp, tb_ham_sp, tb_ham_dd, tb_nn_max, tb_svd, &
                           tnn_tbind, TB_CH_SS, TB_CH_PP, TB_CH_SP, TB_CH_DD, &
                           tb_kmax_desc, tb_lmax_desc, &
                           tb_model_lambda, tb_model_trace, tb_model_equivb, tb_power_trace, &
                           tb_nf_equivb, tb_equivb_poly, tb_lam_obs_min, tb_lam_obs_max, &
                           tbind_dim_per_species
  use mld_mpi, only: mld_rank
  use mld_logger, only: log_critical, log_info, vtoa, mld_verbose, log_debug
  use module_tb_radial, only: tb_rad_ss, tb_rad_pp, tb_rad_sp, tb_rad_dd, tb_channel_radial_t
  use module_ss_hamiltonian, only: tb_ss
  use module_pp_hamiltonian, only: tb_pp
  use module_sp_hamiltonian, only: tb_sp
  use module_dd_hamiltonian, only: tb_dd
#ifdef MLD_NDM
  use notperiod_mod
#else
  use ondm_transform_coord, only: ondm_notperiod
#endif

  implicit none

  integer, intent(in) :: i_start_at, i_final_at
  integer, dimension(imm), intent(out) :: d_n_neigh
  integer, dimension(imm, imm_neigh), intent(out) :: d_kind_neigh
  real(kind_double), dimension(:,:), allocatable :: r_matrix_local, tmp_dxp, tmp_xp
  real(kind_double), dimension(:,:,:), allocatable :: x_matrix_local
  real(kind_double), dimension(:), allocatable :: r_central
  integer, dimension(:), allocatable :: i_type, i_type_db, i_central
  integer, optional :: iconf

  logical :: small
  real(kind_double), dimension(:,:), allocatable :: xpnp
  real(kind_double), dimension(3) :: dxp_ji, ds
  real(kind_double) :: t00, t11

  integer :: ja, iw2, ja_atom, max_neigh_local
  integer :: ian, jan, ixx
  real(kind_double) :: r2_ji, r_ji
  real(kind_double), dimension(3) :: tmp3
  logical :: desc_forces_local

  integer :: icnt_energy, icnt_force, ibegin, iend
  integer :: kk, ll, ll_sigma, ll_pi
  integer :: ll_q_ss, ll_q_sp, ll_q_pps, ll_q_ppp
  integer :: ll_t_sg, ll_t_pi, ll_t_dl
  integer :: mua_a
  integer :: dim_ss, dim_pp, dim_sp, dim_dd, offset_pp
  integer :: dim_lambda_ss, dim_trace_ss, dim_lambda_pp, dim_trace_pp
  integer :: dim_lambda_sp, dim_trace_sp
  integer :: dim_lambda_dd, dim_trace_dd
  integer :: dim_equivb_ss, dim_equivb_pp, dim_equivb_sp, dim_equivb_dd
  integer :: pp, iu_range

  ! Model A sector-trace normalisations 1/sqrt(2l+1) for the single-sector
  ! channels (p: l=1, d: l=2); the sp channel normalises inside its module.
  real(kind_double), parameter :: inv_sqrt3 = 0.57735026918962576d0
  real(kind_double), parameter :: inv_sqrt5 = 0.44721359549995794d0

  _NAMECURRENT_("compute_tbind")
  _MLD_BEGIN_

  if ((i_start_at == 0) .and. (i_final_at == 0)) then
    d_n_neigh(:) = 0
    d_kind_neigh(:,:) = 0
    config_desc(iconf)%energy(:,:) = 0.d0
    return
  end if

  desc_forces_local = desc_forces .and. (config_real(iconf)%has_force .or. config_real(iconf)%has_stress)

  small = .false.
  if (present(iconf)) then
    small = config_real(iconf)%small
  end if

  allocate(xpnp(3, imm))
  if (lperiod) then
    xpnp(:,:) = xp(:,:)
  else
#ifdef MLD_NDM
    call notperiod(imm, xp, xpnp, at, bg, .false.)
#else
    call ondm_notperiod(xp, xpnp)
#endif
  end if

  d_n_neigh(:) = 0
  d_kind_neigh(:,:) = 0
  config_desc(iconf)%energy(:,:) = 0.d0
  if (desc_forces_local) config_desc(iconf)%force(:,:,:,:) = 0.d0

  call preallocate_neigh_ja(r_central, i_type, i_central, tmp_dxp, tmp_xp, imm_neigh)

  ! Compute ss descriptor dimension per f-channel: lambda + trace + equiv-B parts
  dim_lambda_ss = 0
  dim_trace_ss  = 0
  dim_equivb_ss = 0
  if (tb_ham_ss) then
    if (tb_model_lambda) dim_lambda_ss = tb_nn_max(TB_CH_SS)
    if (tb_model_trace)  dim_trace_ss  = tb_power_trace(TB_CH_SS)
    if (tb_model_equivb) dim_equivb_ss = tb_nf_equivb(TB_CH_SS)        ! 1 sector (s)
  end if
  dim_ss = dim_lambda_ss + dim_trace_ss + dim_equivb_ss
  ! pp descriptor dimension per f-channel pair: lambda + trace + equiv-B parts
  dim_lambda_pp = 0
  dim_trace_pp  = 0
  dim_equivb_pp = 0
  if (tb_ham_pp) then
    if (tb_model_lambda) dim_lambda_pp = tb_nn_max(TB_CH_PP)
    if (tb_model_trace)  dim_trace_pp  = tb_power_trace(TB_CH_PP)
    if (tb_model_equivb) dim_equivb_pp = tb_nf_equivb(TB_CH_PP)        ! 1 sector (p)
  end if
  dim_pp = dim_lambda_pp + dim_trace_pp + dim_equivb_pp
  ! sp descriptor dimension per f-channel (4-tuple of l values)
  dim_lambda_sp = 0
  dim_trace_sp  = 0
  dim_equivb_sp = 0
  if (tb_ham_sp) then
    if (tb_model_lambda) dim_lambda_sp = tb_nn_max(TB_CH_SP)
    if (tb_model_trace)  dim_trace_sp  = tb_power_trace(TB_CH_SP)
    if (tb_model_equivb) dim_equivb_sp = 2 * tb_nf_equivb(TB_CH_SP)    ! 2 sectors (s, p)
  end if
  dim_sp = dim_lambda_sp + dim_trace_sp + dim_equivb_sp
  ! dd descriptor dimension per f-channel (triplet of l values: sigma,pi,delta)
  dim_lambda_dd = 0
  dim_trace_dd  = 0
  dim_equivb_dd = 0
  if (tb_ham_dd) then
    if (tb_model_lambda) dim_lambda_dd = tb_nn_max(TB_CH_DD)
    if (tb_model_trace)  dim_trace_dd  = tb_power_trace(TB_CH_DD)
    if (tb_model_equivb) dim_equivb_dd = tb_nf_equivb(TB_CH_DD)        ! 1 sector (d)
  end if
  dim_dd = dim_lambda_dd + dim_trace_dd + dim_equivb_dd

#ifdef MLD_NDM
  iw2 = 0
#else
  if (i_start_at == 1) iw2 = 0
  if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
#endif

  ! Pre-init eigensolvers with a reasonable initial size to avoid reinit overhead.
  ! The init will be called again inside the atom loop only if max_neigh_local changes.
  ! 64 is just a hardcoded initial guess for the number of neighbours.
  ! Note: eigensolver tb_nn_G uses the lambda dimension (eigenvalue truncation).
  ! For trace-only mode (no lambda), we still diagonalise to get eigenvalues.
  if (tb_ham_ss) call tb_ss%init(tb_svd, max(dim_lambda_ss, 1), 64)
  if (tb_ham_pp) call tb_pp%init(tb_svd, max(dim_lambda_pp, 1), 64)
  if (tb_ham_sp) call tb_sp%init(tb_svd, max(dim_lambda_sp, 1), 64)
  if (tb_ham_dd) call tb_dd%init(tb_svd, max(dim_lambda_dd, 1), 64)

  ! Main loop over atoms
  do ja = i_start_at, i_final_at
    ja_atom = ja

    if (debug_time) t00 = MY_MPI_WTIME()

    call build_local_neighbours_ja(iconf, ja, imm, xpnp, r_cut, d_n_neigh, d_kind_neigh, &
                                    r_central, i_type, i_type_db, i_central, tmp_dxp, tmp_xp, max_neigh_local, iw2)

    ! Build pair distance and direction matrices
    if (allocated(r_matrix_local)) deallocate(r_matrix_local)
    allocate(r_matrix_local(0:max_neigh_local, 0:max_neigh_local))
    if (allocated(x_matrix_local)) deallocate(x_matrix_local)
    allocate(x_matrix_local(3, 0:max_neigh_local, 0:max_neigh_local))

    r_matrix_local(:,:) = 0.0d0
    x_matrix_local(:,:,:) = 0.0d0

    do ian = 0, max_neigh_local
      do jan = ian+1, max_neigh_local
        dxp_ji(:) = tmp_xp(1:3, ian) - tmp_xp(1:3, jan)
#ifndef MLD_NDM
        if (.not. small) then
          ds = matmul(dxp_ji, bg)
          where ((ds > 0.5d0) .or. (ds < -0.5d0))
            ds(1:3) = ds(1:3) - dble(nint(ds(1:3)))
          end where
          dxp_ji = matmul(at, ds) / A2cm
        end if
#endif
        r2_ji = sum(dxp_ji(1:3)**2)
        r_ji = dsqrt(r2_ji)
        r_matrix_local(ian, jan) = r_ji
        r_matrix_local(jan, ian) = r_ji
        tmp3(:) = dxp_ji(:) / r_ji
        x_matrix_local(:, ian, jan) = tmp3(:)
        x_matrix_local(:, jan, ian) = -tmp3(:)
      end do
    end do

    ! Central atom species from database numbering (consistent with radial spline indexing)
    mua_a = i_type_db(0)

    ! descriptor index counter -- offset by central-species block (per-element linear model):
    ! atoms of species mua_a write into columns [(mua_a-1)*B + 1 .. mua_a*B] of the descriptor,
    ! where B = tbind_dim_per_species. Other species' blocks remain zero (already initialised).
    icnt_energy = (mua_a - 1) * tbind_dim_per_species

    ! ==================== s-s channel ====================
    if (tb_ham_ss) then
      ! Init hamiltonian for this atom's neighbourhood size
      call tb_ss%init(tb_svd, max(dim_lambda_ss, 1), max_neigh_local)

      do ll = 0, tb_lmax_desc(TB_CH_SS)
        do kk = 1, tb_kmax_desc(TB_CH_SS)
          ! Energy: always compute eigenvalues (needed for both lambda and trace)
          call tb_ss%evaluate_energy(tb_rad_ss, mua_a, i_type_db(1:max_neigh_local), kk, ll, &
                                     r_central(1:max_neigh_local), r_matrix_local)

          ! Observed spectral range diagnostic (for the filter normalization)
          if (tb_model_equivb) then
            tb_lam_obs_min(TB_CH_SS) = min(tb_lam_obs_min(TB_CH_SS), minval(tb_ss%eigen_WW(1:tb_ss%dim_hh)))
            tb_lam_obs_max(TB_CH_SS) = max(tb_lam_obs_max(TB_CH_SS), maxval(tb_ss%eigen_WW(1:tb_ss%dim_hh)))
          end if

          ! Forces: always compute eigenvalue derivatives (needed for both lambda and trace)
          if (desc_forces_local) then
            call tb_ss%evaluate_force(tb_rad_ss, mua_a, i_type_db(1:max_neigh_local), kk, ll, &
                                      r_central(1:max_neigh_local), r_matrix_local, x_matrix_local, tmp_dxp)
          end if

          ! --- Lambda part: eigenvalue descriptors ---
          if (tb_model_lambda) then
            config_desc(iconf)%energy(icnt_energy+1:icnt_energy+dim_lambda_ss, ja) = &
              tb_ss%desc_matrix(1:dim_lambda_ss)
            if (desc_forces_local) then
              ibegin = icnt_energy + 1
              iend   = icnt_energy + dim_lambda_ss
              do ian = 1, max_neigh_local
                do ixx = 1, 3
                  config_desc(iconf)%force(ibegin:iend, ja, ian, ixx) = &
                    tb_ss%force_desc_matrix(ixx, 1:dim_lambda_ss, ian)
                end do
              end do
            end if
            icnt_energy = icnt_energy + dim_lambda_ss
          end if

          ! --- Trace part: Tr(h^p) descriptors ---
          if (tb_model_trace) then
            call tb_ss%compute_trace_energy(dim_trace_ss)
            config_desc(iconf)%energy(icnt_energy+1:icnt_energy+dim_trace_ss, ja) = &
              tb_ss%trace_desc(1:dim_trace_ss)
            if (desc_forces_local) then
              call tb_ss%compute_trace_force(dim_trace_ss)
              ibegin = icnt_energy + 1
              iend   = icnt_energy + dim_trace_ss
              do ian = 1, max_neigh_local
                do ixx = 1, 3
                  config_desc(iconf)%force(ibegin:iend, ja, ian, ixx) = &
                    tb_ss%force_trace_desc(ixx, 1:dim_trace_ss, ian)
                end do
              end do
            end if
            icnt_energy = icnt_energy + dim_trace_ss
          end if

          ! --- equiv-B part (Model A): sector-filtered features
          !     Tr[Pi_l f_k(h)]/sqrt(2l+1). The ss channel has a single s
          !     sector (norm 1). Polynomial filters (default) coincide with
          !     the trace machinery (Sec. 6.7 of main.tex); other spectral
          !     filters go through the generic filter routines. ---
          if (tb_model_equivb) then
            if (tb_equivb_poly) then
              if (.not. tb_model_trace .or. dim_equivb_ss /= dim_trace_ss) &
                call tb_ss%compute_trace_energy(dim_equivb_ss)
              config_desc(iconf)%energy(icnt_energy+1:icnt_energy+dim_equivb_ss, ja) = &
                tb_ss%trace_desc(1:dim_equivb_ss)
            else
              call tb_ss%compute_filter_energy(dim_equivb_ss, 1.0d0)
              config_desc(iconf)%energy(icnt_energy+1:icnt_energy+dim_equivb_ss, ja) = &
                tb_ss%equivb_desc(1:dim_equivb_ss)
            end if
            if (desc_forces_local) then
              ibegin = icnt_energy + 1
              iend   = icnt_energy + dim_equivb_ss
              if (tb_equivb_poly) then
                if (.not. tb_model_trace .or. dim_equivb_ss /= dim_trace_ss) &
                  call tb_ss%compute_trace_force(dim_equivb_ss)
                do ian = 1, max_neigh_local
                  do ixx = 1, 3
                    config_desc(iconf)%force(ibegin:iend, ja, ian, ixx) = &
                      tb_ss%force_trace_desc(ixx, 1:dim_equivb_ss, ian)
                  end do
                end do
              else
                call tb_ss%compute_filter_force(dim_equivb_ss, 1.0d0)
                do ian = 1, max_neigh_local
                  do ixx = 1, 3
                    config_desc(iconf)%force(ibegin:iend, ja, ian, ixx) = &
                      tb_ss%force_equivb_desc(ixx, 1:dim_equivb_ss, ian)
                  end do
                end do
              end if
            end if
            icnt_energy = icnt_energy + dim_equivb_ss
          end if

        end do ! kk
      end do ! ll
    end if

    ! ==================== p-p channel ====================
    if (tb_ham_pp) then
      call tb_pp%init(tb_svd, max(dim_lambda_pp, 1), max_neigh_local)

      do ll = 0, tb_lmax_desc(TB_CH_PP) - 1, 2
        ll_sigma = ll
        ll_pi    = ll + 1
        do kk = 1, tb_kmax_desc(TB_CH_PP)
          ! Energy: always compute eigenvalues (needed for both lambda and trace)
          call tb_pp%evaluate_energy(tb_rad_pp, mua_a, i_type_db(1:max_neigh_local), kk, ll_sigma, ll_pi, &
                                     r_central(1:max_neigh_local), r_matrix_local, x_matrix_local)

          if (tb_model_equivb) then
            tb_lam_obs_min(TB_CH_PP) = min(tb_lam_obs_min(TB_CH_PP), minval(tb_pp%eigen_WW(1:tb_pp%dim_hh)))
            tb_lam_obs_max(TB_CH_PP) = max(tb_lam_obs_max(TB_CH_PP), maxval(tb_pp%eigen_WW(1:tb_pp%dim_hh)))
          end if

          ! Forces: always compute eigenvalue derivatives (needed for both lambda and trace)
          if (desc_forces_local) then
            call tb_pp%evaluate_force(tb_rad_pp, mua_a, i_type_db(1:max_neigh_local), kk, ll_sigma, ll_pi, &
                                      r_central(1:max_neigh_local), r_matrix_local, x_matrix_local, tmp_dxp)
          end if

          ! --- Lambda part: eigenvalue descriptors ---
          if (tb_model_lambda) then
            config_desc(iconf)%energy(icnt_energy+1:icnt_energy+dim_lambda_pp, ja) = &
              tb_pp%desc_matrix(1:dim_lambda_pp)
            if (desc_forces_local) then
              ibegin = icnt_energy + 1
              iend   = icnt_energy + dim_lambda_pp
              do ian = 1, max_neigh_local
                do ixx = 1, 3
                  config_desc(iconf)%force(ibegin:iend, ja, ian, ixx) = &
                    tb_pp%force_desc_matrix(ixx, 1:dim_lambda_pp, ian)
                end do
              end do
            end if
            icnt_energy = icnt_energy + dim_lambda_pp
          end if

          ! --- Trace part: Tr(h^p) descriptors ---
          if (tb_model_trace) then
            call tb_pp%compute_trace_energy(dim_trace_pp)
            config_desc(iconf)%energy(icnt_energy+1:icnt_energy+dim_trace_pp, ja) = &
              tb_pp%trace_desc(1:dim_trace_pp)
            if (desc_forces_local) then
              call tb_pp%compute_trace_force(dim_trace_pp)
              ibegin = icnt_energy + 1
              iend   = icnt_energy + dim_trace_pp
              do ian = 1, max_neigh_local
                do ixx = 1, 3
                  config_desc(iconf)%force(ibegin:iend, ja, ian, ixx) = &
                    tb_pp%force_trace_desc(ixx, 1:dim_trace_pp, ian)
                end do
              end do
            end if
            icnt_energy = icnt_energy + dim_trace_pp
          end if

          ! --- equiv-B part (Model A): the pp channel has a single p sector
          !     (norm 1/sqrt(3)). Polynomial filters reuse the traces; other
          !     spectral filters use the generic filter routines. ---
          if (tb_model_equivb) then
            if (tb_equivb_poly) then
              if (.not. tb_model_trace .or. dim_equivb_pp /= dim_trace_pp) &
                call tb_pp%compute_trace_energy(dim_equivb_pp)
              config_desc(iconf)%energy(icnt_energy+1:icnt_energy+dim_equivb_pp, ja) = &
                tb_pp%trace_desc(1:dim_equivb_pp) * inv_sqrt3
            else
              call tb_pp%compute_filter_energy(dim_equivb_pp, inv_sqrt3)
              config_desc(iconf)%energy(icnt_energy+1:icnt_energy+dim_equivb_pp, ja) = &
                tb_pp%equivb_desc(1:dim_equivb_pp)
            end if
            if (desc_forces_local) then
              ibegin = icnt_energy + 1
              iend   = icnt_energy + dim_equivb_pp
              if (tb_equivb_poly) then
                if (.not. tb_model_trace .or. dim_equivb_pp /= dim_trace_pp) &
                  call tb_pp%compute_trace_force(dim_equivb_pp)
                do ian = 1, max_neigh_local
                  do ixx = 1, 3
                    config_desc(iconf)%force(ibegin:iend, ja, ian, ixx) = &
                      tb_pp%force_trace_desc(ixx, 1:dim_equivb_pp, ian) * inv_sqrt3
                  end do
                end do
              else
                call tb_pp%compute_filter_force(dim_equivb_pp, inv_sqrt3)
                do ian = 1, max_neigh_local
                  do ixx = 1, 3
                    config_desc(iconf)%force(ibegin:iend, ja, ian, ixx) = &
                      tb_pp%force_equivb_desc(ixx, 1:dim_equivb_pp, ian)
                  end do
                end do
              end if
            end if
            icnt_energy = icnt_energy + dim_equivb_pp
          end if

        end do ! kk
      end do ! ll
    end if

    ! ==================== s-p coupling channel ====================
    if (tb_ham_sp) then
      call tb_sp%init(tb_svd, max(dim_lambda_sp, 1), max_neigh_local)

      ! Iterate over quartets of consecutive ll: (ss, sp, ppsigma, pppi)
      do ll = 0, tb_lmax_desc(TB_CH_SP) - 3, 4
        ll_q_ss  = ll
        ll_q_sp  = ll + 1
        ll_q_pps = ll + 2
        ll_q_ppp = ll + 3
        do kk = 1, tb_kmax_desc(TB_CH_SP)
          ! Energy: build h^{sp} and diagonalise (eigenvalues needed for lambda and trace)
          call tb_sp%evaluate_energy(tb_rad_sp, mua_a, i_type_db(1:max_neigh_local), &
                                     kk, ll_q_ss, ll_q_sp, ll_q_pps, ll_q_ppp, &
                                     r_central(1:max_neigh_local), r_matrix_local, x_matrix_local)

          if (tb_model_equivb) then
            tb_lam_obs_min(TB_CH_SP) = min(tb_lam_obs_min(TB_CH_SP), minval(tb_sp%eigen_WW(1:tb_sp%dim_hh)))
            tb_lam_obs_max(TB_CH_SP) = max(tb_lam_obs_max(TB_CH_SP), maxval(tb_sp%eigen_WW(1:tb_sp%dim_hh)))
          end if

          if (desc_forces_local) then
            call tb_sp%evaluate_force(tb_rad_sp, mua_a, i_type_db(1:max_neigh_local), &
                                      kk, ll_q_ss, ll_q_sp, ll_q_pps, ll_q_ppp, &
                                      r_central(1:max_neigh_local), r_matrix_local, x_matrix_local, tmp_dxp)
          end if

          ! --- Lambda part ---
          if (tb_model_lambda) then
            config_desc(iconf)%energy(icnt_energy+1:icnt_energy+dim_lambda_sp, ja) = &
              tb_sp%desc_matrix(1:dim_lambda_sp)
            if (desc_forces_local) then
              ibegin = icnt_energy + 1
              iend   = icnt_energy + dim_lambda_sp
              do ian = 1, max_neigh_local
                do ixx = 1, 3
                  config_desc(iconf)%force(ibegin:iend, ja, ian, ixx) = &
                    tb_sp%force_desc_matrix(ixx, 1:dim_lambda_sp, ian)
                end do
              end do
            end if
            icnt_energy = icnt_energy + dim_lambda_sp
          end if

          ! --- Trace part ---
          if (tb_model_trace) then
            call tb_sp%compute_trace_energy(dim_trace_sp)
            config_desc(iconf)%energy(icnt_energy+1:icnt_energy+dim_trace_sp, ja) = &
              tb_sp%trace_desc(1:dim_trace_sp)
            if (desc_forces_local) then
              call tb_sp%compute_trace_force(dim_trace_sp)
              ibegin = icnt_energy + 1
              iend   = icnt_energy + dim_trace_sp
              do ian = 1, max_neigh_local
                do ixx = 1, 3
                  config_desc(iconf)%force(ibegin:iend, ja, ian, ixx) = &
                    tb_sp%force_trace_desc(ixx, 1:dim_trace_sp, ian)
                end do
              end do
            end if
            icnt_energy = icnt_energy + dim_trace_sp
          end if

          ! --- equiv-B part (Model A): two sectors on sp,
          !     [Tr_s(f_k(h)), Tr_p(f_k(h))/sqrt(3)] per filter k (sector-inner
          !     layout). Forces come from the filtered-W contraction computed
          !     inside sp_evaluate_force (no eigenvector derivatives). ---
          if (tb_model_equivb) then
            if (tb_equivb_poly) then
              call tb_sp%compute_equivb_energy(tb_nf_equivb(TB_CH_SP))
            else
              call tb_sp%compute_filter_energy(tb_nf_equivb(TB_CH_SP))
            end if
            config_desc(iconf)%energy(icnt_energy+1:icnt_energy+dim_equivb_sp, ja) = &
              tb_sp%equivb_desc(1:dim_equivb_sp)
            if (desc_forces_local) then
              ibegin = icnt_energy + 1
              iend   = icnt_energy + dim_equivb_sp
              do ian = 1, max_neigh_local
                do ixx = 1, 3
                  config_desc(iconf)%force(ibegin:iend, ja, ian, ixx) = &
                    tb_sp%force_equivb_desc(ixx, 1:dim_equivb_sp, ian)
                end do
              end do
            end if
            icnt_energy = icnt_energy + dim_equivb_sp
          end if

        end do ! kk
      end do ! ll
    end if

    ! ==================== d-d coupling channel ====================
    if (tb_ham_dd) then
      call tb_dd%init(tb_svd, max(dim_lambda_dd, 1), max_neigh_local)

      ! Iterate over triplets of consecutive ll: (ddsigma, ddpi, dddelta)
      do ll = 0, tb_lmax_desc(TB_CH_DD) - 2, 3
        ll_t_sg = ll
        ll_t_pi = ll + 1
        ll_t_dl = ll + 2
        do kk = 1, tb_kmax_desc(TB_CH_DD)
          call tb_dd%evaluate_energy(tb_rad_dd, mua_a, i_type_db(1:max_neigh_local), &
                                     kk, ll_t_sg, ll_t_pi, ll_t_dl, &
                                     r_central(1:max_neigh_local), r_matrix_local, x_matrix_local)

          if (tb_model_equivb) then
            tb_lam_obs_min(TB_CH_DD) = min(tb_lam_obs_min(TB_CH_DD), minval(tb_dd%eigen_WW(1:tb_dd%dim_hh)))
            tb_lam_obs_max(TB_CH_DD) = max(tb_lam_obs_max(TB_CH_DD), maxval(tb_dd%eigen_WW(1:tb_dd%dim_hh)))
          end if

          if (desc_forces_local) then
            call tb_dd%evaluate_force(tb_rad_dd, mua_a, i_type_db(1:max_neigh_local), &
                                      kk, ll_t_sg, ll_t_pi, ll_t_dl, &
                                      r_central(1:max_neigh_local), r_matrix_local, x_matrix_local, tmp_dxp)
          end if

          if (tb_model_lambda) then
            config_desc(iconf)%energy(icnt_energy+1:icnt_energy+dim_lambda_dd, ja) = &
              tb_dd%desc_matrix(1:dim_lambda_dd)
            if (desc_forces_local) then
              ibegin = icnt_energy + 1
              iend   = icnt_energy + dim_lambda_dd
              do ian = 1, max_neigh_local
                do ixx = 1, 3
                  config_desc(iconf)%force(ibegin:iend, ja, ian, ixx) = &
                    tb_dd%force_desc_matrix(ixx, 1:dim_lambda_dd, ian)
                end do
              end do
            end if
            icnt_energy = icnt_energy + dim_lambda_dd
          end if

          if (tb_model_trace) then
            call tb_dd%compute_trace_energy(dim_trace_dd)
            config_desc(iconf)%energy(icnt_energy+1:icnt_energy+dim_trace_dd, ja) = &
              tb_dd%trace_desc(1:dim_trace_dd)
            if (desc_forces_local) then
              call tb_dd%compute_trace_force(dim_trace_dd)
              ibegin = icnt_energy + 1
              iend   = icnt_energy + dim_trace_dd
              do ian = 1, max_neigh_local
                do ixx = 1, 3
                  config_desc(iconf)%force(ibegin:iend, ja, ian, ixx) = &
                    tb_dd%force_trace_desc(ixx, 1:dim_trace_dd, ian)
                end do
              end do
            end if
            icnt_energy = icnt_energy + dim_trace_dd
          end if

          ! --- equiv-B part (Model A): the dd channel has a single d sector
          !     (norm 1/sqrt(5)). Polynomial filters reuse the traces; other
          !     spectral filters use the generic filter routines. ---
          if (tb_model_equivb) then
            if (tb_equivb_poly) then
              if (.not. tb_model_trace .or. dim_equivb_dd /= dim_trace_dd) &
                call tb_dd%compute_trace_energy(dim_equivb_dd)
              config_desc(iconf)%energy(icnt_energy+1:icnt_energy+dim_equivb_dd, ja) = &
                tb_dd%trace_desc(1:dim_equivb_dd) * inv_sqrt5
            else
              call tb_dd%compute_filter_energy(dim_equivb_dd, inv_sqrt5)
              config_desc(iconf)%energy(icnt_energy+1:icnt_energy+dim_equivb_dd, ja) = &
                tb_dd%equivb_desc(1:dim_equivb_dd)
            end if
            if (desc_forces_local) then
              ibegin = icnt_energy + 1
              iend   = icnt_energy + dim_equivb_dd
              if (tb_equivb_poly) then
                if (.not. tb_model_trace .or. dim_equivb_dd /= dim_trace_dd) &
                  call tb_dd%compute_trace_force(dim_equivb_dd)
                do ian = 1, max_neigh_local
                  do ixx = 1, 3
                    config_desc(iconf)%force(ibegin:iend, ja, ian, ixx) = &
                      tb_dd%force_trace_desc(ixx, 1:dim_equivb_dd, ian) * inv_sqrt5
                  end do
                end do
              else
                call tb_dd%compute_filter_force(dim_equivb_dd, inv_sqrt5)
                do ian = 1, max_neigh_local
                  do ixx = 1, 3
                    config_desc(iconf)%force(ibegin:iend, ja, ian, ixx) = &
                      tb_dd%force_equivb_desc(ixx, 1:dim_equivb_dd, ian)
                  end do
                end do
              end if
            end if
            icnt_energy = icnt_energy + dim_equivb_dd
          end if

        end do ! kk
      end do ! ll
    end if

    ! Translational invariance: force on central atom = - sum of neighbour forces
    if (desc_forces_local) then
      do ian = 1, max_neigh_local
        config_desc(iconf)%force(:, ja, 0, :) = config_desc(iconf)%force(:, ja, 0, :) - &
                                                 config_desc(iconf)%force(:, ja, ian, :)
      end do
    end if

    if (debug_time) then
      t11 = MY_MPI_WTIME()
      tnn_tbind = tnn_tbind + t11 - t00
    end if

  end do ! ja main loop

  ! Observed spectral range diagnostic: overwritten after every configuration,
  ! so at the end of the run the file holds the running extrema seen by rank 0.
  ! Used to set tb_filter_lmin_list / tb_filter_lmax_list for the
  ! non-polynomial spectral filters (padding is applied on top by the code).
  if (tb_model_equivb .and. mld_rank == 0) then
    open(newunit=iu_range, file='tbind_lambda_range.dat', status='replace', action='write')
    write(iu_range, '(A)') '# tbind observed eigenvalue ranges (rank 0 running extrema)'
    write(iu_range, '(A)') '# channel(1=ss,2=pp,3=sp,4=dd)   lambda_min   lambda_max'
    if (tb_ham_ss) write(iu_range, '(I4,2ES20.10)') TB_CH_SS, tb_lam_obs_min(TB_CH_SS), tb_lam_obs_max(TB_CH_SS)
    if (tb_ham_pp) write(iu_range, '(I4,2ES20.10)') TB_CH_PP, tb_lam_obs_min(TB_CH_PP), tb_lam_obs_max(TB_CH_PP)
    if (tb_ham_sp) write(iu_range, '(I4,2ES20.10)') TB_CH_SP, tb_lam_obs_min(TB_CH_SP), tb_lam_obs_max(TB_CH_SP)
    if (tb_ham_dd) write(iu_range, '(I4,2ES20.10)') TB_CH_DD, tb_lam_obs_min(TB_CH_DD), tb_lam_obs_max(TB_CH_DD)
    close(iu_range)
  end if

  deallocate(xpnp)
  if (allocated(r_matrix_local)) deallocate(r_matrix_local)
  if (allocated(x_matrix_local)) deallocate(x_matrix_local)
  _MLD_END_
end subroutine compute_tbind

end module module_compute_tbind


! ============================================================================
! init_tbind: parse parameters, build radial functions, compute dimensions
! ============================================================================
subroutine init_tbind
  use iso_fortran_env, only: dp => real64
  use module_kind_variables, only: kind_double
  use module_tbind
  use module_ace_desc, only: ace_radial_poly
  use module_chemical_species, only: fix_no_of_elements
  use mld_logger, only: log_info, log_critical, vtoa, mld_verbose, log_debug
  use module_tb_radial, only: tb_rad_ss, tb_rad_pp, tb_rad_sp, tb_rad_dd
  use module_tb_filter, only: tbf_init, tbf_type, tbf_nf, tbf_type_name, TBF_POLY

  implicit none

  integer :: ich, dim_mu, nn_max_default
  integer :: npoints_spline
  integer, dimension(TB_N_CHANNELS) :: eff_nn_max
  integer :: n_f_channels_ss, n_f_channels_pp, n_f_channels_sp, n_f_channels_dd
  integer :: kmax_desc_ss, lmax_desc_ss, kmax_desc_pp, lmax_desc_pp
  integer :: kmax_desc_sp, lmax_desc_sp
  integer :: kmax_desc_dd, lmax_desc_dd
  integer :: dim_lambda_ss, dim_trace_ss, dim_equivb_ss, dim_per_f_ss
  integer :: dim_lambda_pp, dim_trace_pp, dim_equivb_pp, dim_per_f_pp
  integer :: dim_lambda_sp, dim_trace_sp, dim_equivb_sp, dim_per_f_sp
  integer :: dim_lambda_dd, dim_trace_dd, dim_equivb_dd, dim_per_f_dd

  _NAMECURRENT_("init_tbind")
  _MLD_BEGIN_

  dim_mu = fix_no_of_elements
  npoints_spline = 500  ! 500 points is sufficient for spline interpolation

  ! Parse string lists into arrays
  call parse_tb_lists()

  ! Resolve the unified tbind_model keyword into the tb_model_* logicals
  ! ('' keeps the legacy tb_model_lambda / tb_model_trace / tb_model_equivb inputs)
  call parse_tbind_model()

  ! Resolve the equiv-B spectral-filter family
  ! (docs/README_tbind_filter_implementation.md). The default 'polynomial'
  ! keeps the legacy raw-lambda^q machinery bit-identical; every other family
  ! is served by module_tb_filter and requires the per-channel reference
  ! spectral intervals tb_filter_lmin_list / tb_filter_lmax_list.
  tb_equivb_poly = .true.
  tb_equivb_rawderiv = .false.
  tb_nf_equivb(:) = 0
  if (tb_model_equivb) then
    call tbf_init(tb_filter_type, tb_filter_order, tb_filter_ncenter, &
                  tb_filter_width, tb_filter_width_ratio, &
                  tb_filter_use_first_moment, tb_filter_use_real_part, tb_filter_use_imag_part, &
                  tb_filter_degree, tb_filter_spectral_padding, &
                  tb_filter_lmin, tb_filter_lmax, &
                  (/ tb_ham_ss, tb_ham_pp, tb_ham_sp, tb_ham_dd /))
    tb_equivb_poly = (tbf_type == TBF_POLY)
    tb_equivb_rawderiv = .not. tb_equivb_poly
    if (tb_equivb_poly) then
      ! legacy: powers from tb_power_trace, optionally overridden by tb_filter_order
      do ich = 1, TB_N_CHANNELS
        if (tb_filter_order > 0) then
          tb_nf_equivb(ich) = tb_filter_order
        else
          tb_nf_equivb(ich) = tb_power_trace(ich)
        end if
      end do
    else
      tb_nf_equivb(:) = tbf_nf
    end if
  end if

  ! Set default nn_max where -1
  nn_max_default = 20
  do ich = 1, TB_N_CHANNELS
    if (tb_nn_max(ich) < 0) then
      eff_nn_max(ich) = nn_max_default
    else
      eff_nn_max(ich) = tb_nn_max(ich)
    end if
  end do
  ! Write back effective values
  tb_nn_max(:) = eff_nn_max(:)

  ! Legacy aliases
  tb_bss    = tb_ham_ss
  tb_bpp    = tb_ham_pp
  tb_nn_G   = eff_nn_max(TB_CH_SS)
  k_param   = 1
  tb_ngf    = tb_nrf

  ! Compute descriptor dimension
  tbind_dim = 0

  ! ======= ss channel =======
  if (tb_ham_ss) then
    call log_info("TB: building ss channel radials (HSVD)")
    call tb_rad_ss%init(TB_CH_SS, tb_kmax(TB_CH_SS), tb_nmax(TB_CH_SS), tb_lmax(TB_CH_SS), &
                         dim_mu, tb_rcut_in(TB_CH_SS), tb_rcut_out(TB_CH_SS), &
                         tb_rcut_width_in(TB_CH_SS), tb_rcut_width_out(TB_CH_SS), &
                         tb_lambda_arr(TB_CH_SS), npoints_spline, tb_g_type, &
                         tb_rcut_in_ij(TB_CH_SS), tb_rcut_out_ij(TB_CH_SS), &
                         tb_rcut_width_in_ij(TB_CH_SS), tb_rcut_width_out_ij(TB_CH_SS))
    call tb_rad_ss%build_f(ace_radial_poly)
    if (tb_g_type == 1) call tb_rad_ss%build_g()
    if (tb_g_type == 2) call log_info('TB ss: g_type=2 (g=f, no Annex window)')
    ! Descriptor grid: use tb_kmax_desc/tb_lmax_desc (may be a subset of radial basis)
    kmax_desc_ss = tb_kmax_desc(TB_CH_SS)
    lmax_desc_ss = tb_lmax_desc(TB_CH_SS)
    if (kmax_desc_ss > tb_kmax(TB_CH_SS)) then
      call log_critical("TB ss: kmax_grid ("//vtoa(kmax_desc_ss)//") > kmax_list ("//vtoa(tb_kmax(TB_CH_SS))//")")
      stop
    end if
    if (lmax_desc_ss > tb_lmax(TB_CH_SS)) then
      call log_critical("TB ss: lmax_grid ("//vtoa(lmax_desc_ss)//") > lmax_list ("//vtoa(tb_lmax(TB_CH_SS))//")")
      stop
    end if
    ! ss: dimension per f-channel = lambda part + trace part + equiv-B part
    n_f_channels_ss = kmax_desc_ss * (lmax_desc_ss + 1)
    dim_lambda_ss = 0
    dim_trace_ss  = 0
    dim_equivb_ss = 0
    if (tb_model_lambda) dim_lambda_ss = eff_nn_max(TB_CH_SS)
    if (tb_model_trace)  dim_trace_ss  = tb_power_trace(TB_CH_SS)
    if (tb_model_equivb) dim_equivb_ss = tb_nf_equivb(TB_CH_SS)     ! 1 sector (s)
    dim_per_f_ss = dim_lambda_ss + dim_trace_ss + dim_equivb_ss
    tbind_dim = tbind_dim + n_f_channels_ss * dim_per_f_ss
    call log_info("TB ss radials: kmax="//vtoa(tb_kmax(TB_CH_SS))//" lmax="//vtoa(tb_lmax(TB_CH_SS))// &
                  " n_f="//vtoa(tb_rad_ss%n_f))
    call log_info("TB ss descriptor grid: kmax_grid="//vtoa(kmax_desc_ss)//" lmax_grid="//vtoa(lmax_desc_ss)// &
                  " n_f_desc="//vtoa(n_f_channels_ss)//" nn_max="//vtoa(eff_nn_max(TB_CH_SS))// &
                  " dim_lambda="//vtoa(n_f_channels_ss * dim_lambda_ss)// &
                  " dim_trace="//vtoa(n_f_channels_ss * dim_trace_ss)// &
                  " dim_equivb="//vtoa(n_f_channels_ss * dim_equivb_ss)// &
                  " dim="//vtoa(n_f_channels_ss * dim_per_f_ss))
  end if

  ! ======= pp channel =======
  if (tb_ham_pp) then
    call log_info("TB: building pp channel radials (HSVD)")
    ! For pp, lmax must be odd: sigma=ll, pi=ll+1 => lmax/2 pairs
    if (mod(tb_lmax(TB_CH_PP), 2) /= 1) then
      call log_info("TB pp: WARNING - lmax should be odd (sigma/pi pairs). Setting lmax="// &
                    vtoa(tb_lmax(TB_CH_PP)+1))
      tb_lmax(TB_CH_PP) = tb_lmax(TB_CH_PP) + 1
    end if
    call tb_rad_pp%init(TB_CH_PP, tb_kmax(TB_CH_PP), tb_nmax(TB_CH_PP), tb_lmax(TB_CH_PP), &
                         dim_mu, tb_rcut_in(TB_CH_PP), tb_rcut_out(TB_CH_PP), &
                         tb_rcut_width_in(TB_CH_PP), tb_rcut_width_out(TB_CH_PP), &
                         tb_lambda_arr(TB_CH_PP), npoints_spline, tb_g_type, &
                         tb_rcut_in_ij(TB_CH_PP), tb_rcut_out_ij(TB_CH_PP), &
                         tb_rcut_width_in_ij(TB_CH_PP), tb_rcut_width_out_ij(TB_CH_PP))
    call tb_rad_pp%build_f(ace_radial_poly)
    if (tb_g_type == 1) call tb_rad_pp%build_g()
    if (tb_g_type == 2) call log_info('TB pp: g_type=2 (g=f, no Annex window)')
    ! Descriptor grid for pp
    kmax_desc_pp = tb_kmax_desc(TB_CH_PP)
    lmax_desc_pp = tb_lmax_desc(TB_CH_PP)
    ! Ensure pp descriptor lmax is odd too
    if (mod(lmax_desc_pp, 2) /= 1) then
      call log_info("TB pp: WARNING - lmax_grid should be odd (sigma/pi pairs). Setting lmax_grid="// &
                    vtoa(lmax_desc_pp+1))
      lmax_desc_pp = lmax_desc_pp + 1
      tb_lmax_desc(TB_CH_PP) = lmax_desc_pp
    end if
    if (kmax_desc_pp > tb_kmax(TB_CH_PP)) then
      call log_critical("TB pp: kmax_grid ("//vtoa(kmax_desc_pp)//") > kmax_list ("//vtoa(tb_kmax(TB_CH_PP))//")")
      stop
    end if
    if (lmax_desc_pp > tb_lmax(TB_CH_PP)) then
      call log_critical("TB pp: lmax_grid ("//vtoa(lmax_desc_pp)//") > lmax_list ("//vtoa(tb_lmax(TB_CH_PP))//")")
      stop
    end if
    ! pp: dimension per f-channel = lambda part + trace part + equiv-B part
    n_f_channels_pp = kmax_desc_pp * ((lmax_desc_pp + 1) / 2)
    dim_lambda_pp = 0
    dim_trace_pp  = 0
    dim_equivb_pp = 0
    if (tb_model_lambda) dim_lambda_pp = eff_nn_max(TB_CH_PP)
    if (tb_model_trace)  dim_trace_pp  = tb_power_trace(TB_CH_PP)
    if (tb_model_equivb) dim_equivb_pp = tb_nf_equivb(TB_CH_PP)     ! 1 sector (p)
    dim_per_f_pp = dim_lambda_pp + dim_trace_pp + dim_equivb_pp
    tbind_dim = tbind_dim + n_f_channels_pp * dim_per_f_pp
    call log_info("TB pp radials: kmax="//vtoa(tb_kmax(TB_CH_PP))//" lmax="//vtoa(tb_lmax(TB_CH_PP))// &
                  " n_f="//vtoa(tb_rad_pp%n_f))
    call log_info("TB pp descriptor grid: kmax_grid="//vtoa(kmax_desc_pp)//" lmax_grid="//vtoa(lmax_desc_pp)// &
                  " n_f_desc="//vtoa(n_f_channels_pp)//" nn_max="//vtoa(eff_nn_max(TB_CH_PP))// &
                  " dim_lambda="//vtoa(n_f_channels_pp * dim_lambda_pp)// &
                  " dim_trace="//vtoa(n_f_channels_pp * dim_trace_pp)// &
                  " dim_equivb="//vtoa(n_f_channels_pp * dim_equivb_pp)// &
                  " dim="//vtoa(n_f_channels_pp * dim_per_f_pp))
  end if

  ! ======= sp channel =======
  if (tb_ham_sp) then
    call log_info("TB: building sp channel radials (HSVD)")
    ! For sp, the descriptor consumes quartets (ll_ss, ll_sp, ll_pps, ll_ppp) so the
    ! radial-basis lmax must satisfy (lmax + 1) mod 4 == 0, i.e. lmax in {3, 7, 11, ...}.
    if (mod(tb_lmax(TB_CH_SP) + 1, 4) /= 0) then
      call log_info("TB sp: WARNING - lmax must be of the form 4m+3 (ss/sp/ppsigma/pppi quartets). Setting lmax="// &
                    vtoa(((tb_lmax(TB_CH_SP) + 1) / 4 + 1) * 4 - 1))
      tb_lmax(TB_CH_SP) = ((tb_lmax(TB_CH_SP) + 1) / 4 + 1) * 4 - 1
    end if
    call tb_rad_sp%init(TB_CH_SP, tb_kmax(TB_CH_SP), tb_nmax(TB_CH_SP), tb_lmax(TB_CH_SP), &
                         dim_mu, tb_rcut_in(TB_CH_SP), tb_rcut_out(TB_CH_SP), &
                         tb_rcut_width_in(TB_CH_SP), tb_rcut_width_out(TB_CH_SP), &
                         tb_lambda_arr(TB_CH_SP), npoints_spline, tb_g_type, &
                         tb_rcut_in_ij(TB_CH_SP), tb_rcut_out_ij(TB_CH_SP), &
                         tb_rcut_width_in_ij(TB_CH_SP), tb_rcut_width_out_ij(TB_CH_SP))
    call tb_rad_sp%build_f(ace_radial_poly)
    if (tb_g_type == 1) call tb_rad_sp%build_g()
    if (tb_g_type == 2) call log_info('TB sp: g_type=2 (g=f, no Annex window)')
    ! Descriptor grid for sp
    kmax_desc_sp = tb_kmax_desc(TB_CH_SP)
    lmax_desc_sp = tb_lmax_desc(TB_CH_SP)
    if (mod(lmax_desc_sp + 1, 4) /= 0) then
      call log_info("TB sp: WARNING - lmax_grid must be of the form 4m+3. Setting lmax_grid="// &
                    vtoa(((lmax_desc_sp + 1) / 4 + 1) * 4 - 1))
      lmax_desc_sp = ((lmax_desc_sp + 1) / 4 + 1) * 4 - 1
      tb_lmax_desc(TB_CH_SP) = lmax_desc_sp
    end if
    if (kmax_desc_sp > tb_kmax(TB_CH_SP)) then
      call log_critical("TB sp: kmax_grid ("//vtoa(kmax_desc_sp)//") > kmax_list ("//vtoa(tb_kmax(TB_CH_SP))//")")
      stop
    end if
    if (lmax_desc_sp > tb_lmax(TB_CH_SP)) then
      call log_critical("TB sp: lmax_grid ("//vtoa(lmax_desc_sp)//") > lmax_list ("//vtoa(tb_lmax(TB_CH_SP))//")")
      stop
    end if
    ! sp: dimension per f-channel = lambda + trace + equiv-B; n_f counts quartets
    n_f_channels_sp = kmax_desc_sp * ((lmax_desc_sp + 1) / 4)
    dim_lambda_sp = 0
    dim_trace_sp  = 0
    dim_equivb_sp = 0
    if (tb_model_lambda) dim_lambda_sp = eff_nn_max(TB_CH_SP)
    if (tb_model_trace)  dim_trace_sp  = tb_power_trace(TB_CH_SP)
    if (tb_model_equivb) dim_equivb_sp = 2 * tb_nf_equivb(TB_CH_SP)    ! 2 sectors (s, p)
    dim_per_f_sp = dim_lambda_sp + dim_trace_sp + dim_equivb_sp
    tbind_dim = tbind_dim + n_f_channels_sp * dim_per_f_sp
    call log_info("TB sp radials: kmax="//vtoa(tb_kmax(TB_CH_SP))//" lmax="//vtoa(tb_lmax(TB_CH_SP))// &
                  " n_f="//vtoa(tb_rad_sp%n_f))
    call log_info("TB sp descriptor grid: kmax_grid="//vtoa(kmax_desc_sp)//" lmax_grid="//vtoa(lmax_desc_sp)// &
                  " n_f_desc="//vtoa(n_f_channels_sp)//" nn_max="//vtoa(eff_nn_max(TB_CH_SP))// &
                  " dim_lambda="//vtoa(n_f_channels_sp * dim_lambda_sp)// &
                  " dim_trace="//vtoa(n_f_channels_sp * dim_trace_sp)// &
                  " dim_equivb="//vtoa(n_f_channels_sp * dim_equivb_sp)// &
                  " dim="//vtoa(n_f_channels_sp * dim_per_f_sp))
  end if

  ! ======= dd channel =======
  if (tb_ham_dd) then
    call log_info("TB: building dd channel radials (HSVD)")
    ! For dd, the descriptor consumes triplets (ll_sg, ll_pi, ll_dl) so the
    ! radial-basis lmax must satisfy (lmax + 1) mod 3 == 0, i.e. lmax in {2, 5, 8, ...}.
    if (mod(tb_lmax(TB_CH_DD) + 1, 3) /= 0) then
      call log_info("TB dd: WARNING - lmax must be of the form 3m+2 (sigma/pi/delta triplets). Setting lmax="// &
                    vtoa(((tb_lmax(TB_CH_DD) + 1) / 3 + 1) * 3 - 1))
      tb_lmax(TB_CH_DD) = ((tb_lmax(TB_CH_DD) + 1) / 3 + 1) * 3 - 1
    end if
    call tb_rad_dd%init(TB_CH_DD, tb_kmax(TB_CH_DD), tb_nmax(TB_CH_DD), tb_lmax(TB_CH_DD), &
                         dim_mu, tb_rcut_in(TB_CH_DD), tb_rcut_out(TB_CH_DD), &
                         tb_rcut_width_in(TB_CH_DD), tb_rcut_width_out(TB_CH_DD), &
                         tb_lambda_arr(TB_CH_DD), npoints_spline, tb_g_type, &
                         tb_rcut_in_ij(TB_CH_DD), tb_rcut_out_ij(TB_CH_DD), &
                         tb_rcut_width_in_ij(TB_CH_DD), tb_rcut_width_out_ij(TB_CH_DD))
    call tb_rad_dd%build_f(ace_radial_poly)
    if (tb_g_type == 1) call tb_rad_dd%build_g()
    if (tb_g_type == 2) call log_info('TB dd: g_type=2 (g=f, no Annex window)')
    kmax_desc_dd = tb_kmax_desc(TB_CH_DD)
    lmax_desc_dd = tb_lmax_desc(TB_CH_DD)
    if (mod(lmax_desc_dd + 1, 3) /= 0) then
      call log_info("TB dd: WARNING - lmax_grid must be of the form 3m+2. Setting lmax_grid="// &
                    vtoa(((lmax_desc_dd + 1) / 3 + 1) * 3 - 1))
      lmax_desc_dd = ((lmax_desc_dd + 1) / 3 + 1) * 3 - 1
      tb_lmax_desc(TB_CH_DD) = lmax_desc_dd
    end if
    if (kmax_desc_dd > tb_kmax(TB_CH_DD)) then
      call log_critical("TB dd: kmax_grid ("//vtoa(kmax_desc_dd)//") > kmax_list ("//vtoa(tb_kmax(TB_CH_DD))//")")
      stop
    end if
    if (lmax_desc_dd > tb_lmax(TB_CH_DD)) then
      call log_critical("TB dd: lmax_grid ("//vtoa(lmax_desc_dd)//") > lmax_list ("//vtoa(tb_lmax(TB_CH_DD))//")")
      stop
    end if
    n_f_channels_dd = kmax_desc_dd * ((lmax_desc_dd + 1) / 3)
    dim_lambda_dd = 0
    dim_trace_dd  = 0
    dim_equivb_dd = 0
    if (tb_model_lambda) dim_lambda_dd = eff_nn_max(TB_CH_DD)
    if (tb_model_trace)  dim_trace_dd  = tb_power_trace(TB_CH_DD)
    if (tb_model_equivb) dim_equivb_dd = tb_nf_equivb(TB_CH_DD)     ! 1 sector (d)
    dim_per_f_dd = dim_lambda_dd + dim_trace_dd + dim_equivb_dd
    tbind_dim = tbind_dim + n_f_channels_dd * dim_per_f_dd
    call log_info("TB dd radials: kmax="//vtoa(tb_kmax(TB_CH_DD))//" lmax="//vtoa(tb_lmax(TB_CH_DD))// &
                  " n_f="//vtoa(tb_rad_dd%n_f))
    call log_info("TB dd descriptor grid: kmax_grid="//vtoa(kmax_desc_dd)//" lmax_grid="//vtoa(lmax_desc_dd)// &
                  " n_f_desc="//vtoa(n_f_channels_dd)//" nn_max="//vtoa(eff_nn_max(TB_CH_DD))// &
                  " dim_lambda="//vtoa(n_f_channels_dd * dim_lambda_dd)// &
                  " dim_trace="//vtoa(n_f_channels_dd * dim_trace_dd)// &
                  " dim_equivb="//vtoa(n_f_channels_dd * dim_equivb_dd)// &
                  " dim="//vtoa(n_f_channels_dd * dim_per_f_dd))
  end if

  ! Validate model selection
  if (.not. tb_model_lambda .and. .not. tb_model_trace .and. .not. tb_model_equivb) then
    call log_critical("TB: no active model (tb_model_lambda / tb_model_trace / tb_model_equivb "// &
                      "all .false., tbind_model empty). No descriptors!")
    stop
  end if
  if (tb_model_lambda) call log_info("TB: model_lambda = .true. (eigenvalue descriptors)")
  if (tb_model_trace) then
    call log_info("TB: model_trace = .true. (trace descriptors)")
    if (tb_ham_ss) call log_info("TB ss: trace power = "//vtoa(tb_power_trace(TB_CH_SS)))
    if (tb_ham_pp) call log_info("TB pp: trace power = "//vtoa(tb_power_trace(TB_CH_PP)))
    if (tb_ham_sp) call log_info("TB sp: trace power = "//vtoa(tb_power_trace(TB_CH_SP)))
    if (tb_ham_dd) call log_info("TB dd: trace power = "//vtoa(tb_power_trace(TB_CH_DD)))
  end if
  if (tb_model_equivb) then
    call log_info("TB: model_equivb = .true. (equiv-B = Model A: sector features "// &
                  "Tr[Pi_l f_k(h)]/sqrt(2l+1), tb_filter_type='"//trim(tbf_type_name(tbf_type))//"')")
    if (tb_ham_ss) call log_info("TB ss: equiv-B n_filters = "//vtoa(tb_nf_equivb(TB_CH_SS))//" (1 sector)")
    if (tb_ham_pp) call log_info("TB pp: equiv-B n_filters = "//vtoa(tb_nf_equivb(TB_CH_PP))//" (1 sector)")
    if (tb_ham_sp) call log_info("TB sp: equiv-B n_filters = "//vtoa(tb_nf_equivb(TB_CH_SP))//" (2 sectors: s, p)")
    if (tb_ham_dd) call log_info("TB dd: equiv-B n_filters = "//vtoa(tb_nf_equivb(TB_CH_DD))//" (1 sector)")
    if (tb_model_trace .and. tb_equivb_poly) then
      call log_info("TB WARNING: trace and polynomial equiv-B are both active. On the "// &
                    "single-sector channels (ss, pp, dd) the equiv-B features are the trace "// &
                    "features up to the constant 1/sqrt(2l+1): the linear system may be "// &
                    "rank deficient.")
    end if
    if (.not. tb_ham_sp .and. tb_equivb_poly) then
      call log_info("TB note: without the sp channel, polynomial equiv-B carries the same "// &
                    "information as the trace model (single orbital sector per channel).")
    end if
  end if

  ! Per-central-species block: enlarge total dimension by fix_no_of_elements so that
  ! each central species gets its own linear-weight block (per-element model).
  tbind_dim_per_species = tbind_dim
  tbind_dim = fix_no_of_elements * tbind_dim_per_species

  call log_info("TB: per-species block dim = "//vtoa(tbind_dim_per_species)// &
                "  N_species = "//vtoa(fix_no_of_elements)// &
                "  total descriptor dim = "//vtoa(tbind_dim))

  _MLD_END_

contains

  subroutine parse_tb_lists()
    implicit none
    integer :: jch

    ! Allocate parsed arrays
    if (allocated(tb_kmax)) deallocate(tb_kmax)
    allocate(tb_kmax(TB_N_CHANNELS))
    if (allocated(tb_nmax)) deallocate(tb_nmax)
    allocate(tb_nmax(TB_N_CHANNELS))
    if (allocated(tb_lmax)) deallocate(tb_lmax)
    allocate(tb_lmax(TB_N_CHANNELS))
    if (allocated(tb_kmax_desc)) deallocate(tb_kmax_desc)
    allocate(tb_kmax_desc(TB_N_CHANNELS))
    if (allocated(tb_lmax_desc)) deallocate(tb_lmax_desc)
    allocate(tb_lmax_desc(TB_N_CHANNELS))
    if (allocated(tb_lambda_arr)) deallocate(tb_lambda_arr)
    allocate(tb_lambda_arr(TB_N_CHANNELS))
    if (allocated(tb_rcut_in)) deallocate(tb_rcut_in)
    allocate(tb_rcut_in(TB_N_CHANNELS))
    if (allocated(tb_rcut_out)) deallocate(tb_rcut_out)
    allocate(tb_rcut_out(TB_N_CHANNELS))
    if (allocated(tb_rcut_width_in)) deallocate(tb_rcut_width_in)
    allocate(tb_rcut_width_in(TB_N_CHANNELS))
    if (allocated(tb_rcut_width_out)) deallocate(tb_rcut_width_out)
    allocate(tb_rcut_width_out(TB_N_CHANNELS))
    if (allocated(tb_rcut_in_ij)) deallocate(tb_rcut_in_ij)
    allocate(tb_rcut_in_ij(TB_N_CHANNELS))
    if (allocated(tb_rcut_out_ij)) deallocate(tb_rcut_out_ij)
    allocate(tb_rcut_out_ij(TB_N_CHANNELS))
    if (allocated(tb_rcut_width_in_ij)) deallocate(tb_rcut_width_in_ij)
    allocate(tb_rcut_width_in_ij(TB_N_CHANNELS))
    if (allocated(tb_rcut_width_out_ij)) deallocate(tb_rcut_width_out_ij)
    allocate(tb_rcut_width_out_ij(TB_N_CHANNELS))

    ! Parse integer lists (radial basis construction)
    read(tb_kmax_list, *) tb_kmax(1:TB_N_CHANNELS)
    read(tb_nmax_list, *) tb_nmax(1:TB_N_CHANNELS)
    read(tb_lmax_list, *) tb_lmax(1:TB_N_CHANNELS)

    ! Parse descriptor grid lists
    read(tb_kmax_grid, *) tb_kmax_desc(1:TB_N_CHANNELS)
    read(tb_lmax_grid, *) tb_lmax_desc(1:TB_N_CHANNELS)
    ! Default: if -1, use the same as the radial basis
    do jch = 1, TB_N_CHANNELS
      if (tb_kmax_desc(jch) < 0) tb_kmax_desc(jch) = tb_kmax(jch)
      if (tb_lmax_desc(jch) < 0) tb_lmax_desc(jch) = tb_lmax(jch)
    end do

    ! Parse real lists
    read(tb_lambda_list, *) tb_lambda_arr(1:TB_N_CHANNELS)
    read(tb_rcut_in_list, *) tb_rcut_in(1:TB_N_CHANNELS)
    read(tb_rcut_out_list, *) tb_rcut_out(1:TB_N_CHANNELS)
    read(tb_rcut_width_in_list, *) tb_rcut_width_in(1:TB_N_CHANNELS)
    read(tb_rcut_width_out_list, *) tb_rcut_width_out(1:TB_N_CHANNELS)

    ! Parse the internal (ij) cutoff family for the hopping f(r_ij) splines
    ! (r_cut^(ij) of docs/Tbind_desc/main.tex). Negative entries (default -1)
    ! inherit the corresponding central (a) value -> legacy behaviour unchanged.
    read(tb_rcut_in_ij_list, *) tb_rcut_in_ij(1:TB_N_CHANNELS)
    read(tb_rcut_out_ij_list, *) tb_rcut_out_ij(1:TB_N_CHANNELS)
    read(tb_rcut_width_in_ij_list, *) tb_rcut_width_in_ij(1:TB_N_CHANNELS)
    read(tb_rcut_width_out_ij_list, *) tb_rcut_width_out_ij(1:TB_N_CHANNELS)
    do jch = 1, TB_N_CHANNELS
      if (tb_rcut_in_ij(jch)        < 0.0d0) tb_rcut_in_ij(jch)        = tb_rcut_in(jch)
      if (tb_rcut_out_ij(jch)       < 0.0d0) tb_rcut_out_ij(jch)       = tb_rcut_out(jch)
      if (tb_rcut_width_in_ij(jch)  < 0.0d0) tb_rcut_width_in_ij(jch)  = tb_rcut_width_in(jch)
      if (tb_rcut_width_out_ij(jch) < 0.0d0) tb_rcut_width_out_ij(jch) = tb_rcut_width_out(jch)
    end do
    if (any(abs(tb_rcut_in_ij - tb_rcut_in) > 1.0d-14) .or. &
        any(abs(tb_rcut_out_ij - tb_rcut_out) > 1.0d-14) .or. &
        any(abs(tb_rcut_width_in_ij - tb_rcut_width_in) > 1.0d-14) .or. &
        any(abs(tb_rcut_width_out_ij - tb_rcut_width_out) > 1.0d-14)) then
      do jch = 1, TB_N_CHANNELS
        call log_info("TB ch"//vtoa(jch)//" (ij) cutoffs: rcut_in_ij="//vtoa(tb_rcut_in_ij(jch))// &
                      " rcut_out_ij="//vtoa(tb_rcut_out_ij(jch))// &
                      " width_in_ij="//vtoa(tb_rcut_width_in_ij(jch))// &
                      " width_out_ij="//vtoa(tb_rcut_width_out_ij(jch)))
      end do
      if (tb_g_type == 2) then
        call log_info("TB WARNING: tb_g_type=2 (g=f) with distinct (ij) cutoffs: the central "// &
                      "envelope g(r_ai) follows the f splines and therefore the (ij) range.")
      end if
    end if

    ! Parse trace power list
    if (allocated(tb_power_trace)) deallocate(tb_power_trace)
    allocate(tb_power_trace(TB_N_CHANNELS))
    read(tb_power_trace_list, *) tb_power_trace(1:TB_N_CHANNELS)

    ! Parse the reference spectral intervals for the equiv-B filters
    if (allocated(tb_filter_lmin)) deallocate(tb_filter_lmin)
    allocate(tb_filter_lmin(TB_N_CHANNELS))
    if (allocated(tb_filter_lmax)) deallocate(tb_filter_lmax)
    allocate(tb_filter_lmax(TB_N_CHANNELS))
    read(tb_filter_lmin_list, *) tb_filter_lmin(1:TB_N_CHANNELS)
    read(tb_filter_lmax_list, *) tb_filter_lmax(1:TB_N_CHANNELS)

    call log_info("TB radial basis: kmax="//vtoa(tb_kmax(1))//" "//vtoa(tb_kmax(2))// &
                  " "//vtoa(tb_kmax(3))//" "//vtoa(tb_kmax(4)))
    call log_info("TB radial basis: nmax="//vtoa(tb_nmax(1))//" "//vtoa(tb_nmax(2))// &
                  " "//vtoa(tb_nmax(3))//" "//vtoa(tb_nmax(4)))
    call log_info("TB radial basis: lmax="//vtoa(tb_lmax(1))//" "//vtoa(tb_lmax(2))// &
                  " "//vtoa(tb_lmax(3))//" "//vtoa(tb_lmax(4)))
    call log_info("TB descriptor grid: kmax="//vtoa(tb_kmax_desc(1))//" "//vtoa(tb_kmax_desc(2))// &
                  " "//vtoa(tb_kmax_desc(3))//" "//vtoa(tb_kmax_desc(4)))
    call log_info("TB descriptor grid: lmax="//vtoa(tb_lmax_desc(1))//" "//vtoa(tb_lmax_desc(2))// &
                  " "//vtoa(tb_lmax_desc(3))//" "//vtoa(tb_lmax_desc(4)))
  end subroutine parse_tb_lists

  ! ---------------------------------------------------------------------------
  ! Resolve the unified tbind_model keyword into the tb_model_* logicals.
  ! Accepted (case-insensitive) tokens, combined with '+':
  !   'lambda'  -> tb_model_lambda   (eigenvalue descriptors)
  !   'trace'   -> tb_model_trace    (trace descriptors)
  !   'equiv-B' -> tb_model_equivb   (Model A sector traces; 'equivb' and
  !                                   'equiv_b' are accepted aliases)
  ! Examples: 'lambda', 'trace', 'lambda+trace', 'equiv-B', 'trace+equiv-B'.
  ! An empty tbind_model (default) keeps the legacy logical inputs untouched.
  ! ---------------------------------------------------------------------------
  subroutine parse_tbind_model()
    implicit none
    character(len=64) :: work, token
    integer :: ip, istart
    logical :: has_lambda, has_trace, has_equivb

    work = adjustl(tbind_model)
    if (len_trim(work) == 0) return   ! legacy logicals stay in charge

    call lower_string(work)
    has_lambda = .false.
    has_trace  = .false.
    has_equivb = .false.

    istart = 1
    do
      ip = index(work(istart:), '+')
      if (ip == 0) then
        token = adjustl(work(istart:))
      else if (ip == 1) then
        token = ''
      else
        token = adjustl(work(istart:istart+ip-2))
      end if
      select case (trim(token))
      case ('lambda')
        has_lambda = .true.
      case ('trace')
        has_trace = .true.
      case ('equiv-b', 'equivb', 'equiv_b')
        has_equivb = .true.
      case ('')
        ! empty token (repeated or trailing '+'): ignore
      case default
        call log_critical("TB: unknown tbind_model token '"//trim(token)// &
                          "' in tbind_model='"//trim(tbind_model)// &
                          "' (allowed: lambda, trace, equiv-B, '+'-combinations)")
        stop
      end select
      if (ip == 0) exit
      istart = istart + ip
    end do

    tb_model_lambda = has_lambda
    tb_model_trace  = has_trace
    tb_model_equivb = has_equivb
    call log_info("TB: tbind_model = '"//trim(tbind_model)//"' -> lambda="// &
                  merge('T', 'F', tb_model_lambda)//" trace="//merge('T', 'F', tb_model_trace)// &
                  " equiv-B="//merge('T', 'F', tb_model_equivb))
  end subroutine parse_tbind_model

  subroutine lower_string(s)
    implicit none
    character(len=*), intent(inout) :: s
    integer :: i, ic
    do i = 1, len_trim(s)
      ic = iachar(s(i:i))
      if (ic >= iachar('A') .and. ic <= iachar('Z')) s(i:i) = achar(ic + 32)
    end do
  end subroutine lower_string

end subroutine init_tbind
