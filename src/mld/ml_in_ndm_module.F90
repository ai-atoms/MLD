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

#include "../MLD_MACROS.INC"


module data_type
   implicit none
   integer(kind=4), parameter :: IB = 4, RP = 8
   integer, save  :: icall
   integer  :: nunit
end module data_type

module module_units
    use, intrinsic :: iso_fortran_env, dp=>real64
    implicit none 
    real(dp), parameter :: two_pi = 8.0_dp*atan(1.0_dp)
    real(dp), parameter :: one_pi = 4.0_dp*atan(1.0_dp)

    real(dp), parameter :: CHARGE = 1.602176634e-19_dp ! 1 e in C units
    real(dp), parameter :: eV_TO_J = CHARGE            ! 1 eV in J units
    real(dp), parameter :: J_TO_eV = 1.0_dp/eV_TO_J    ! 1 J in eV units

    real(dp), parameter :: Erg_TO_J = 1.0e-7_dp ! 1 erg in J units
    real(dp), parameter :: J_TO_Erg = 1.0e+7_dp ! 1 J in erg units

    real(dp), parameter :: Erg_TO_eV = Erg_TO_J*J_TO_eV ! 1 erg in eV units
    real(dp), parameter :: eV_TO_Erg = 1.0_dp/Erg_TO_eV ! 1 eV in erg 
    
    real(dp), parameter :: kB_eVonK =  8.617333262e-5_dp ! Boltzmann constant in  eV K-1
    real(dp), parameter :: kB_JonK  = kB_eVonK*CHARGE       ! Boltzmann constant in J K-1

    real(dp), parameter :: K_TO_eV = kB_eVonK              ! 1 K in eV units 
    real(dp), parameter :: eV_TO_K = 1.0_dp/K_TO_eV     ! 1 eV in K units 

    real(dp), parameter :: A_TO_CM = 1.0e-8_dp           ! 1 A in cm units
    real(dp), parameter :: CM_TO_A = 1.0e+8_dp           ! 1 cm in A units

    real(dp), parameter :: h_eVs = 4.135667696e-15_dp    ! Planck constant in eV s
    real(dp), parameter :: hbar_eVs = h_eVs/two_pi       ! Reduced Planck constant in eV s
    real(dp), parameter :: THz_TO_eV= h_eVs * 1.0e+12_dp !THz in eV units
    real(dp), parameter :: eV_TO_Thz=1.0_dp/THz_TO_eV    !eV in THz units

    real(dp), parameter :: AMUNITSKG = 1.66053906892e-27_dp ! atomic mass unit in kg
    real(dp), parameter :: AMUNITSGR = AMUNITSKG / 1.0e-3_dp ! atomic mass unit in g
    real(dp), parameter :: AMASS = AMUNITSKG / CHARGE*1.0e+10_dp ! atomic mass unit in eV A^{-2} fs^2   
    
end module module_units 

module module_end_ml
   implicit none
contains

   subroutine end_ml (text, scalapack_driver, context, rangml)

      !TORC! use mpi
      use mld_mpi, only: codeml, comm_mld  
      !TORC! use my_mpi_subroutines, only:my_barrier_mld
      ! use gen_mpi
      integer, intent(in):: rangml, context
      character(len=*), intent(in)     :: text
      logical, intent(in):: scalapack_driver



      if (rangml == 0) then
         write (6, *) '  ', trim(text), '  '
         write (6, '("ML:   ----------------c u later alligator----------------  ")')
      end if

      if (rangml == 0) write (6, *) '************  END DE MACHINE LEARNING ****************'
#if(PARA)
      !TORC! call my_barrier_mld(codeml)
      call comm_mld%barrier
#endif

      if (scalapack_driver) then
         call blacs_barrier(context, 'A')
         call blacs_gridexit( context )
      end if
#if(PARA)
      !TORC! call my_barrier_mld(codeml)
      call comm_mld%barrier
      call MPI_FINALIZE(codeml)
#endif
      stop
   end subroutine end_ml
end module module_end_ml

module module_scalapack_interfaces
  implicit none 
  contains 

    subroutine my_pdgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt)
      INTEGER, INTENT(IN) :: m, n
      INTEGER, INTENT(IN) :: ia, ja, ib, jb
      INTEGER, INTENT(IN) :: ictxt
      INTEGER, INTENT(IN) :: desca(*), descb(*)
      DOUBLE PRECISION, INTENT(IN)  :: a(*)
      DOUBLE PRECISION, INTENT(OUT) :: b(*)

      ! Call the actual ScaLAPACK routine (assuming it's available and linked)
      call pdgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt)
    end subroutine my_pdgemr2d
  

end module module_scalapack_interfaces



module module_ml_scalapack
   logical :: scalapack_driver, debug_scalapack
   ! the id of proc returned by scalapck, both are used in parallel
   integer  :: iam, iproc_sca
   ! the total numbers of procs used by scalapack
   integer  :: nprocs_ml_sca
   ! the total number of row and col of the ScaLapack grid.
   integer  :: nprow, npcol
   ! the local index for row and col in the ScaLapack grid
   integer  :: myrow, mycol
   ! size of the block for row and col.
   integer :: nbr_predefined, nbc_predefined
   ! sca_Amat matrix is Scalapack version of design matrix
   integer  :: nbr_Amat, nbc_Amat
   real(kind=8), dimension(:, :), allocatable   :: sca_Amat
   integer, dimension(:), allocatable     :: desc_sca_Amat
   integer  :: dimr_sca_Amat, dimc_sca_Amat
   integer  :: l_dimr_sca_Amat, l_dimc_sca_Amat

   !sca_AmatT ... the transpose of sca_Amat. Why? Parce que c'est comme ca! 
   integer  :: nbr_AmatT, nbc_AmatT
   real(kind=8), dimension(:, :), allocatable   :: sca_AmatT
   integer, dimension(:), allocatable     :: desc_sca_AmatT
   integer  :: dimr_sca_AmatT, dimc_sca_AmatT
   integer  :: l_dimr_sca_AmatT, l_dimc_sca_AmatT

   integer  :: nbr_Amat_big, nbc_Amat_big
   real(kind=8), dimension(:, :), allocatable   :: sca_Amat_big
   integer, dimension(:), allocatable     :: desc_sca_Amat_big
   integer  :: dimr_sca_Amat_big, dimc_sca_Amat_big
   integer  :: l_dimr_sca_Amat_big, l_dimc_sca_Amat_big

   integer  :: nbr_Lmat, nbc_Lmat
   real(kind=8), dimension(:, :), allocatable   :: sca_Lmat
   integer, dimension(:), allocatable     :: desc_sca_Lmat
   integer  :: dimr_sca_Lmat, dimc_sca_Lmat
   integer  :: l_dimr_sca_Lmat, l_dimc_sca_Lmat

   ! sca_Qmat matrix is Scalapack version of Qmat from Amat = Lmat + Qmat
   integer  :: nbr_Qmat, nbc_Qmat
   real(kind=8), dimension(:, :), allocatable   :: sca_Qmat
   integer, dimension(:), allocatable     :: desc_sca_Qmat
   integer  :: dimr_sca_Qmat, dimc_sca_Qmat
   integer  :: l_dimr_sca_Qmat, l_dimc_sca_Qmat


   ! sca_Cmat matrix is Scalapack version of design matrix * weigths
   integer  :: nbr_Cmat, nbc_Cmat
   real(kind=8), dimension(:, :), allocatable   :: sca_Cmat
   integer, dimension(:), allocatable     :: desc_sca_Cmat
   integer  :: dimr_sca_Cmat, dimc_sca_Cmat
   integer  :: l_dimr_sca_Cmat, l_dimc_sca_Cmat

   ! sca_phi matrix is Scalapack version of Amat * Amat^T = D x D
   integer  :: nbr_phi, nbc_phi
   real(kind=8), dimension(:, :), allocatable   :: sca_phi
   integer, dimension(:), allocatable     :: desc_sca_phi
   integer  :: dimr_sca_phi, dimc_sca_phi
   integer  :: l_dimr_sca_phi, l_dimc_sca_phi

   ! sca_phia matrix is Scalapack version of atomic design matrix M_a x D.
   integer  :: nbr_phia, nbc_phia
   real(kind=8), dimension(:, :), allocatable   :: sca_phia
   integer, dimension(:), allocatable     :: desc_sca_phia
   integer  :: dimr_sca_phia, dimc_sca_phia
   integer  :: l_dimr_sca_phia, l_dimc_sca_phia


   ! sca_ymat  is the Scalapack version of ymat M x 1
   integer  :: nbr_ymat, nbc_ymat ! the last one is 1 but for symmetry of implemnetation we define that  as variable.
   real(kind=8), dimension(:, :), allocatable   :: sca_ymat
   integer, dimension(:), allocatable     :: desc_sca_ymat
   integer  :: dimr_sca_ymat, dimc_sca_ymat    ! the last one is 1, for symmetry
   integer  :: l_dimr_sca_ymat, l_dimc_sca_ymat ! the last one is 1, for symmetry

   ! sca_ymat_copy  is the Scalapack version of ymat M x 1
   integer  :: nbr_ymat_copy, nbc_ymat_copy ! the last one is 1 but for symmetry of implemnetation we define that  as variable.
   real(kind=8), dimension(:, :), allocatable   :: sca_ymat_copy
   integer, dimension(:), allocatable     :: desc_sca_ymat_copy
   integer  :: dimr_sca_ymat_copy, dimc_sca_ymat_copy    ! the last one is 1, for symmetry
   integer  :: l_dimr_sca_ymat_copy, l_dimc_sca_ymat_copy ! the last one is 1, for symmetry


   ! sca_ymat_qr_svd  is the Scalapack version of ymat_qr_svd M x 1
   integer  :: nbr_ymat_qr_svd, nbc_ymat_qr_svd ! the last one is 1
   real(kind=8), dimension(:, :), allocatable   :: sca_ymat_qr_svd
   integer, dimension(:), allocatable     :: desc_sca_ymat_qr_svd
   integer  :: dimr_sca_ymat_qr_svd, dimc_sca_ymat_qr_svd     ! the last one is 1, for symmetry
   integer  :: l_dimr_sca_ymat_qr_svd, l_dimc_sca_ymat_qr_svd ! the last one is 1, for symmetry


   integer  :: nbr_w_params, nbc_w_params ! the last one is 1
   real(kind=8), dimension(:, :), allocatable   :: sca_w_params
   integer, dimension(:), allocatable     :: desc_sca_w_params
   integer  :: dimr_sca_w_params, dimc_sca_w_params     ! the last one is 1, for symmetry
   integer  :: l_dimr_sca_w_params, l_dimc_sca_w_params  ! the last one is 1, for symmetry

   integer  :: context
   integer, parameter   :: descriptor_len = 9

end module module_ml_scalapack

module module_cur
   use module_kind_variables, only: kind_double
   integer :: cur_kval, cur_cval , cur_rval
   real(kind_double) :: cur_eps
   integer, parameter :: cur_no_of_samples = 10
   real(kind_double), dimension(:), allocatable :: pcol, prow
   logical, dimension(:), allocatable :: selcol, selrow
   real(kind_double), dimension(:,:), allocatable :: pcol_sample, prow_sample
   logical, dimension(:,:), allocatable :: selcol_sample, selrow_sample
   integer, dimension(:), allocatable :: norow_sample, nocol_sample, cur_info_selection

   integer :: col_no_of_selections, row_no_of_selections
   integer :: rank_sca_phia, rank_sca_phia_T

   type type_info_matrix
      integer :: iconf
      integer :: ia
   end type

   type(type_info_matrix), dimension(:), allocatable :: info_mat, info_mat_ref, info_mat_full

end module module_cur

module module_lbfgs_input
   real(kind(1.d0))     :: lbfgs_xtol, lbfgs_eps, lbfgs_gtol
   integer  :: lbfgs_m_hess, lbfgs_max_steps
   integer, dimension(2)      :: lbfgs_print
end module module_lbfgs_input


module time_check_general
   use module_kind_variables, only: kind_double
#if(PARA)
   use mpi
   use mld_mpi 
#endif

   implicit none

   logical  :: debug_time
   real(kind_double)    :: time_fit_params, &
      time_full_desc, &
      time_calc_desc, &
      time_neigh_desc, &
      time_fill_desc, &
      time_fill_desc_02, &
      time_fill_desc_02y, &
      time_fill_desc_02A, &
      time_fill_desc_02obj, &
      time_read_db
   real(kind_double)    :: time_scatter_blacs, &
      time_scatter_gather, &
      time_scatter_bcast, &
      time_scatter_dgemm, &
      time_scatter_barrier
   real(kind_double)    :: test_time_total, &
      test_time_desc, &
      test_time_eval, &
      test_time_neigh
   real(kind_double), dimension(20) :: time
   real(kind_double), dimension(20) :: tot_time, train_tot_time, test_tot_time
   integer, dimension(3)      :: it_counter

contains

!para! #if(PARA)
   function MY_MPI_WTIME()
      real(kind_double)    :: MY_MPI_WTIME
      MY_MPI_WTIME = MPI_WTIME()
   end function
!para! #else
!para!    function MY_MPI_WTIME()
!para!       real(kind_double)    :: MY_MPI_WTIME, res
!para!       call cpu_time(res)
!para!       MY_MPI_WTIME = res      ! TODO time fortan std ?
!para!    end function
!para! #endif

end module time_check_general





module module_afs
   use module_kind_variables, only: kind_double
   integer  :: afs_dim, n_rbf, n_rbf_afs, n_cheb
   integer  :: afs_type
   integer, parameter   :: afs_type_bartok = 1, afs_type_homemade = 2

   real(kind_double), allocatable, dimension(:) :: coeff_rbf_afs
   real(kind_double), allocatable, dimension(:, :)    :: W_afs

   ! new version
   real(kind_double), dimension(:, :), allocatable    :: drbf_afs, rbf_afs
end module module_afs



module module_bispectrum_so4
   use module_kind_variables, ONLY: kind_double

   double complex, parameter  :: czero = cmplx(0.d0, 0.d0, kind=kind(1.d0))
   integer  :: bisso4_cg_dim, bisso4_cg_full_dim, bisso4_cg_A_dim, bisso4_cg_B_dim, bisso4_cg_C_dim
   double complex, allocatable, dimension(:, :, :)    :: &
      cmm, cmm2, cmm3, &
      ZAcmm, ZAcmm2, ZAcmm3, &
      ZBcmm, ZBcmm2, ZBcmm3, &
      ZCcmm, ZCcmm2, ZCcmm3

   double precision, allocatable, dimension(:)  :: cg_A, cg_B, cg_C
   integer, allocatable, dimension(:)     :: bisso4_l, bisso4_l1, bisso4_l2
   integer  :: bisso4_dim
   real(kind_double)    :: r0, fcut, dfcut, fcut_w, dfcut_w, fcut_w_3ch, dfcut_w_3ch
   real(kind_double)    :: inv_r0, inv_r0_input
   logical  :: lbso4_diag
   type mml4D
      double complex, dimension(:), allocatable    :: vecall, vecall_w, vecall_w_3ch
   end type

   type(mml4D), dimension(:, :, :), allocatable, target     :: class_mml
end module module_bispectrum_so4



module module_so3
   use module_kind_variables, only: kind_double, kind_double_complex

   integer, parameter   :: radial_bartok = 1, &
      radial_sgg = 2, &
      radial_ace_bessel = 3

   integer  :: pow_so3_dim
   integer  :: radial_pow_so3, ini_rbf_so3, end_rbf_so3
   integer  ::  n_rbf_so3

   complex(kind_double_complex), allocatable, dimension(:, :, :, :)     :: dclmn
   complex(kind_double_complex), allocatable, dimension(:, :)     :: clmn
   complex(kind_double_complex), allocatable, dimension(:, :, :, :)     :: dclmn_w
   complex(kind_double_complex), allocatable, dimension(:, :)     :: clmn_w
   complex(kind_double_complex), allocatable, dimension(:, :, :, :)     :: dclmn_w_3ch
   complex(kind_double_complex), allocatable, dimension(:, :)     :: clmn_w_3ch
   real(kind_double), dimension(:), allocatable :: coeff_rbf_so3
   real(kind_double), dimension(:, :), allocatable    :: W_pow_so3

   real(kind_double), dimension(:), allocatable     :: chebT, d_chebT                  ! T Pftouny first kind
   real(kind_double), dimension(:), allocatable     :: chebU  ! T Pftouny second kind
   real(kind_double), dimension(:), allocatable     :: scgg, d_scgg
   real(kind_double), allocatable, dimension(:) :: radial_ace, d_radial_ace
   real(kind_double) :: r_cut_so3_in, r_cut_so3_width_in

end module module_so3

module module_tbind 
   use module_kind_variables, only: kind_double
   implicit none

   integer :: tbind_dim       ! total descriptor dimension (= tbind_dim_per_species * fix_no_of_elements)
   integer :: tbind_dim_per_species   ! single-species (core) descriptor dim; per-central-species block size

   ! Channel activation flags
   logical :: tb_ham_ss = .false.   ! activate s-s channel
   logical :: tb_ham_pp = .false.   ! activate p-p channel
   logical :: tb_ham_sp = .false.   ! activate s-p channel
   logical :: tb_ham_dd = .false.   ! activate d-d channel

   ! Number of channels (ss=1, pp=2, sp=3, dd=4)
   integer, parameter :: TB_CH_SS = 1, TB_CH_PP = 2, TB_CH_SP = 3, TB_CH_DD = 4
   integer, parameter :: TB_N_CHANNELS = 4

   ! Eigenvalue cutoff per channel (default 20, -1 means use default)
   integer, dimension(TB_N_CHANNELS) :: tb_nn_max = (/ -1, -1, -1, -1 /)

   ! ACE HSVD radial parameters per channel
   ! These select the k, n, l basis size for the radial functions f_n and g  
   character(len=5000) :: tb_kmax_list  = "8 8 8 8"
   character(len=5000) :: tb_nmax_list  = "14 4 2 1"
   character(len=5000) :: tb_lmax_list  = "3 2 2 1"
   character(len=5000) :: tb_lambda_list = "3.0 3.0 3.0 3.0"

   ! Descriptor grid parameters per channel
   ! These control which (k,l) pairs are used for the descriptor (can be a subset of the radial basis)
   ! If not set (default "-1"), they default to tb_kmax_list / tb_lmax_list
   character(len=5000) :: tb_kmax_grid = "-1 -1 -1 -1"
   character(len=5000) :: tb_lmax_grid = "-1 -1 -1 -1"

   ! Cutoff parameters per channel
   character(len=5000) :: tb_rcut_in_list        = "1.2d0 1.2d0 1.2d0 1.2d0"
   character(len=5000) :: tb_rcut_out_list       = "5.d0 5.d0 5.d0 5.d0"
   character(len=5000) :: tb_rcut_width_in_list  = "0.4d0 0.4d0 0.4d0 0.4d0"
   character(len=5000) :: tb_rcut_width_out_list = "0.5d0 0.5d0 0.5d0 0.5d0"

   ! Internal (i-j) cutoff family: r_cut^(ij) of docs/Tbind_desc/main.tex (Sec. 2.1).
   ! These act on the hopping splines f(r_ij) only (build range + evaluation guard);
   ! the central envelope g(r_ai) keeps the (a) family above.
   ! Default -1 = inherit the corresponding central (a) value -> legacy behaviour
   ! is bit-identical when the _ij keywords are not set.
   character(len=5000) :: tb_rcut_in_ij_list        = "-1 -1 -1 -1"
   character(len=5000) :: tb_rcut_out_ij_list       = "-1 -1 -1 -1"
   character(len=5000) :: tb_rcut_width_in_ij_list  = "-1 -1 -1 -1"
   character(len=5000) :: tb_rcut_width_out_ij_list = "-1 -1 -1 -1"

   ! Parsed arrays (allocated in init_tbind)
   integer, dimension(:), allocatable :: tb_kmax, tb_nmax, tb_lmax
   integer, dimension(:), allocatable :: tb_kmax_desc, tb_lmax_desc  ! descriptor grid (k,l) bounds
   real(kind_double), dimension(:), allocatable :: tb_lambda_arr
   real(kind_double), dimension(:), allocatable :: tb_rcut_in, tb_rcut_out
   real(kind_double), dimension(:), allocatable :: tb_rcut_width_in, tb_rcut_width_out
   real(kind_double), dimension(:), allocatable :: tb_rcut_in_ij, tb_rcut_out_ij
   real(kind_double), dimension(:), allocatable :: tb_rcut_width_in_ij, tb_rcut_width_out_ij

   ! Radial function dimensions
   integer :: tb_nrf = 1         ! number of f radial functions (n_f in the doc)
   integer :: tb_nrg = 1         ! number of g radial functions (n_g in the doc)

   ! Model selection: which descriptor components to include
   logical :: tb_model_lambda = .true.   ! eigenvalue-based descriptors (lambda^2, Section 8 of Tbind_desc.pdf)
   logical :: tb_model_trace  = .false.  ! trace-based descriptors (Tr(h^p), Section 7 of Tbind_desc.pdf)
   logical :: tb_model_equivb = .false.  ! equiv-B descriptors = Model A of docs/Tbind_desc/main.tex:
                                         ! orbital-sector traces Tr[Pi_l B_a^(cq)]/sqrt(2l+1) with
                                         ! polynomial filters f_q(lambda)=lambda^q, q=1..tb_power_trace

   ! Unified model keyword. When non-empty it OVERRIDES the three logicals above.
   ! Accepted values (case-insensitive): 'lambda', 'trace', 'equiv-B', and any
   ! '+'-combination such as 'lambda+trace' or 'trace+equiv-B'.
   ! Default '' keeps the legacy tb_model_lambda / tb_model_trace behaviour.
   character(len=64) :: tbind_model = ''

   ! Power of trace: max power p for Tr(h^p) per channel (ss, pp, sp, dd).
   ! The same powers are reused as the polynomial-filter orders of equiv-B.
   character(len=5000) :: tb_power_trace_list = "4 4 4 4"
   integer, dimension(:), allocatable :: tb_power_trace  ! parsed array (allocated in init_tbind)

   ! ------------------------------------------------------------------
   ! Spectral filters for the equiv-B (Model A) descriptors
   ! (docs/README_tbind_filter_implementation.md). The default
   ! 'polynomial' reproduces the raw lambda^q equiv-B implementation
   ! bit-identically; the other families act on the normalized
   ! eigenvalue x = (lambda - lambda0_c)/Delta_c per channel c.
   ! ------------------------------------------------------------------
   ! Filter family: polynomial (default) | chebyshev | gaussian |
   !                lorentzian | resolvent | fermi | fermi_bin | bspline
   character(len=64) :: tb_filter_type = 'polynomial'
   ! Polynomial/Chebyshev order Q. -1 = per-type default:
   !   polynomial -> tb_power_trace per channel (legacy powers q=1..p)
   !   chebyshev  -> Q=10 (T_0..T_Q)
   integer :: tb_filter_order = -1
   ! Number of centers/thresholds/basis functions K. -1 = per-type default:
   !   gaussian 8, lorentzian 12, resolvent 8, fermi 16, fermi_bin 16, bspline 16
   integer :: tb_filter_ncenter = -1
   ! Absolute width (sigma/eta/tau) in normalized x units; -1 = use ratio
   real(kind_double) :: tb_filter_width = -1.0d0
   ! Width relative to the center spacing h. -1 = per-type default:
   !   gaussian 1.0, lorentzian 0.75, resolvent 0.75, fermi 0.6, fermi_bin 0.6
   real(kind_double) :: tb_filter_width_ratio = -1.0d0
   logical :: tb_filter_use_first_moment = .true.   ! gaussian: add x*G_k filters
   logical :: tb_filter_use_real_part = .true.      ! resolvent: dispersive Q_k
   logical :: tb_filter_use_imag_part = .true.      ! resolvent: Lorentzian L_k
   integer :: tb_filter_degree = 3                  ! bspline degree (3 only)
   ! Padding applied to [lmin,lmax] when building the normalization:
   ! Delta_c = (1+pad)*(lmax-lmin)/2
   real(kind_double) :: tb_filter_spectral_padding = 0.15d0
   ! Reference spectral interval per channel (ss, pp, sp, dd); required for
   ! every non-polynomial filter. Run once with equiv-B active to obtain the
   ! observed ranges in 'tbind_lambda_range.dat', then paste them here (the
   ! padding above is applied on top).
   character(len=5000) :: tb_filter_lmin_list = "0 0 0 0"
   character(len=5000) :: tb_filter_lmax_list = "0 0 0 0"
   real(kind_double), dimension(:), allocatable :: tb_filter_lmin, tb_filter_lmax

   ! Resolved equiv-B bookkeeping (set in init_tbind):
   integer, dimension(TB_N_CHANNELS) :: tb_nf_equivb = 0  ! filters per channel
   logical :: tb_equivb_poly = .true.       ! legacy raw-polynomial path active
   logical :: tb_equivb_rawderiv = .false.  ! channels must store raw d(lambda)/dx
   ! Observed spectral range diagnostic (updated when equiv-B is active)
   real(kind_double), dimension(TB_N_CHANNELS) :: tb_lam_obs_min =  1.0d30
   real(kind_double), dimension(TB_N_CHANNELS) :: tb_lam_obs_max = -1.0d30

   ! g function type: 1 = Annex smoothing window (default), 2 = ACE g_k radial
   integer :: tb_g_type = 1

   ! Eigensolver type: 1 - SVD, 2 - DSYEVD (default), 3 - DSYEV, 4 - DSYEVR
   integer :: tb_svd = 2

   ! Timing accumulator
   real(kind_double) :: tnn_tbind = 0.0_kind_double

   ! Legacy aliases (kept for backward compatibility in compute_tbind internals)
   integer :: tb_nn_G          ! effective eigenvalue cutoff (set per channel)
   logical :: tb_bss, tb_bpp   ! internal copies from tb_ham_ss/tb_ham_pp
   integer :: tb_type_ff       ! legacy: radial type for f (unused in new code)
   integer :: tb_type_gg       ! legacy: radial type for g (unused in new code)
   integer :: tb_ngf            ! legacy alias for tb_nrf
   integer :: k_param           ! legacy

end module module_tbind

module module_ftnbody
   use module_kind_variables, only: kind_double
   integer, parameter :: MAX_ORDER_FTNBODY = 5
   integer :: ftnbody_dim

   integer :: dim_rff(MAX_ORDER_FTNBODY), dim_qbody(MAX_ORDER_FTNBODY)
   integer :: length_order = 0

   ! --- geometry model (docs/MLT5, sec. "Permutation-invariant Gram-coordinate
   !     Fourier representation" and sec. "PIP-RFF family of ftnbody
   !     descriptors") ---
   ! ftnbody_model:
   !   'poly'  = compact permutation-invariant polynomial q (historical, default;
   !             'pip' is accepted as an alias: it is the PIP level of the family)
   !   'gramm' = ordered Gram coordinates (u_1..u_m, c_pq) with exact orbit
   !             averaging over the S_{n-1} neighbour permutations after the
   !             Fourier map
   !   'cpip'  = compact mixed invariant set (Level B): selected radial/angular/
   !             mixed moments, plus (chem mode 2) colored chemistry-geometry
   !             moments inside the RFF input; chem modes 0 and 2 only
   !   'spip'  = systematic orbit-polynomial set (Level C): all S_{n-1} monomial
   !             orbits up to spip_degree_nbody; chem modes 0,1,2,3
   integer, parameter :: FTNBODY_MODEL_POLY  = 0
   integer, parameter :: FTNBODY_MODEL_GRAMM = 1
   integer, parameter :: FTNBODY_MODEL_CPIP  = 2
   integer, parameter :: FTNBODY_MODEL_SPIP  = 3
   character(len=16) :: ftnbody_model = 'poly'      ! namelist keyword string
   integer :: ftnbody_model_id = FTNBODY_MODEL_POLY ! parsed value
   ! Gram-model permutation tables, built in init_ftnbody:
   !   gram_nperm(n)      = (n-1)!            number of neighbour permutations
   !   gram_perm(p,ip,n)  = pi(p)             slot permutation ip of order n
   !   gram_qmap(i,ip,n)  = index map so that (Pi_pi q)_i = q(gram_qmap(i,ip,n))
   integer :: gram_nperm(MAX_ORDER_FTNBODY) = 1
   integer :: gram_perm(MAX_ORDER_FTNBODY-1, 24, MAX_ORDER_FTNBODY) = 0
   integer :: gram_qmap(10, 24, MAX_ORDER_FTNBODY) = 0

   ! --- cPIP model (ftnbody_model='cpip'): versioned compact basis v1 ---
   ! Per order n (m=n-1 neighbours): CPIP_MG_N(n) geometric moments and, in
   ! chem mode 2 only, CPIP_MC_N(n) colored chemistry-geometry moments appended
   ! to the RFF input (docs/MLT5, sec. "Level B: cPIP").
   !   n=2: (u);                       colored: (C^u, C^u2)
   !   n=3: (u1+u2, u1*u2, c12);       colored: (C^u, C^u2, C^uc, C^c, C^uuc)
   !   n=4,5: (Su, Su2, Su3, Sc, Sc2, Sc3, P7, P8, P9) with
   !     P7 = sum_i u_i sum_{j/=i} c_ij,  P8 = sum_i u_i prod_{j<k, j,k/=i} c_jk,
   !     P9 = sum_{i<j} u_i u_j c_ij;    colored: (C^u, C^u2, C^uc, C^c, C^uuc)
   integer, parameter, dimension(2:5) :: CPIP_MG_N = (/ 1, 3, 9, 9 /)
   integer, parameter, dimension(2:5) :: CPIP_MC_N = (/ 2, 5, 5, 5 /)

   ! --- sPIP model (ftnbody_model='spip'): monomial orbit tables ---
   ! For each order n, all monomials x^e in the primitive coordinates
   ! x = (u_1..u_m, c_12..c_{m-1,m}) with 1 <= |e| <= spip_degree_n(n) are
   ! grouped into S_m permutation orbits (built once in init_ftnbody):
   !   P_alpha(x) = (1/m!) sum_pi M_alpha(Pi_pi x)
   ! One representative exponent vector per orbit (rep_expo) plus the collapsed
   ! distinct-member list (mono_expo / mono_coef, coef = multiplicity/m!) so
   ! that P_alpha = sum_k mono_coef(k) * x^mono_expo(:,k).
   type :: spip_orbit_tab
      integer :: norb  = 0     ! M_n = number of distinct orbits
      integer :: nmono = 0     ! total collapsed members over all orbits
      integer, allocatable :: orb_ptr(:)      ! (norb+1) CSR pointer into members
      integer, allocatable :: mono_expo(:, :) ! (dimx, nmono) member exponents
      real(kind_double), allocatable :: mono_coef(:) ! (nmono) multiplicity/m!
      integer, allocatable :: rep_expo(:, :)  ! (dimx, norb) orbit representative
   end type spip_orbit_tab
   type(spip_orbit_tab) :: spip_tab(MAX_ORDER_FTNBODY)
   ! namelist keyword string (1-based per-order layout, like dim_fourier_nbody)
   character(len=80) :: spip_degree_nbody = " 0 4 4 4 3 "  ! D_n max degree, spip
   integer :: spip_degree_n(MAX_ORDER_FTNBODY) = 0         ! parsed

   ! --- multispecies chemistry modes (see docs/perspective_ftnbody.md) ---
   ! ftnbody_chem_mode:
   !   0 = single-species / current compact behaviour (default)
   !   1 = exact multispecies chemical channels (reference/debug)
   !   2 = low-rank multispecies chemical embedding (recommended, scalable)
   !   3 = hashed multispecies chemical channels (approximate)
   integer :: ftnbody_chem_mode = 0
   ! namelist keyword strings (1-based per-order layout, like dim_fourier_nbody):
   character(len=80) :: ftnbody_chem_rank     = " 0 0 0 0 0 "  ! R_n string, mode 2
   character(len=80) :: ftnbody_hash_channels = " 0 0 0 0 0 "  ! H_n string, mode 3
   integer :: ftnbody_chem_rank_n(MAX_ORDER_FTNBODY)     = 0  ! R_n, mode 2 (parsed)
   integer :: ftnbody_hash_channels_n(MAX_ORDER_FTNBODY) = 0  ! H_n, mode 3 (parsed)
   ! number of chemical channels actually used for each order:
   !   mode 0 -> 1 ; mode 1 -> N_chem(n,S) ; mode 2 -> R_n ; mode 3 -> H_n
   ! and dim_desc_body(n) = dim_rff(n) * ftnbody_n_channels(n)
   integer :: ftnbody_n_channels(MAX_ORDER_FTNBODY)    = 1
   ! low-rank role embeddings (mode 2): centre (S, R_n, order) and neighbour (S, R_n, slot, order)
   real(kind_double), dimension(:, :, :),    allocatable :: ftnbody_chem_center
   real(kind_double), dimension(:, :, :, :), allocatable :: ftnbody_chem_neigh

   logical :: init_mode_ftnbody = .false.
   logical :: ltmp_mean_ftnbody = .false.
   logical :: find_best_length = .false.
   ! when .true. (MD / prediction from a saved potential) the descriptor random
   ! features are loaded from random.xml instead of being (re-)estimated/generated
   logical :: ftnbody_read_from_potential = .false.

   real(kind_double) :: length_rff(MAX_ORDER_FTNBODY) , delta_rff(MAX_ORDER_FTNBODY)
   real(kind_double), dimension(:), allocatable :: ftnbody_phase_random_2b, ftnbody_phase_random_3b
   real(kind_double), dimension(:), allocatable :: ftnbody_phase_random_4b, ftnbody_phase_random_5b

   real(kind_double), dimension(:, :), allocatable :: ftnbody_omega_2b, ftnbody_omega_3b
   real(kind_double), dimension(:, :), allocatable :: ftnbody_omega_4b, transpose_ftnbody_omega_4b, ftnbody_omega_5b

   real(kind_double) :: count_2b, count_3b, count_4b, count_5b
   real(kind_double) :: number_2b_ftnb, number_3b_ftnb, number_4b_ftnb, number_5b_ftnb
   real(kind_double), dimension(:, :), allocatable :: covar_matrix_2b, covar_matrix_3b, covar_matrix_4b, covar_matrix_5b
   real(kind_double), dimension(:), allocatable :: length_rff_2b, length_rff_3b, length_rff_4b, length_rff_5b
   real(kind_double), dimension(:), allocatable :: mean_2b, mean_3b, mean_4b, mean_5b
   
   real(kind_double) :: tnn_ftbd, t2b_ftbd, t3b_ftbd, t4b_ftbd, t5b_ftbd, &
      tnn_ftbd_train, t2b_ftbd_train, t3b_ftbd_train, t4b_ftbd_train, t5b_ftbd_train, &
      t4b_inner_init, t4b_inner_desc, t4b_inner_deriv, t4b_inner_last

  real(kind_double) :: r_cut_ft2b, r_cut_width_ft2b, r_cut_ft3b, r_cut_width_ft3b, &
                       r_cut_ft4b, r_cut_width_ft4b, r_cut_ft5b, r_cut_width_ft5b

   real(kind_double), dimension(:), allocatable :: tmp_real
   real(kind_double), dimension(:,:,:), allocatable :: tmp_dreal

contains

   ! position of the angular coordinate c_pq (1 <= p < q <= m) inside the
   ! Gram vector q_G = (u_1..u_m, c_12, c_13, ..., c_{m-1,m}) (lexicographic)
   pure integer function gram_pair_pos(m, p, q) result(ipos)
      integer, intent(in) :: m, p, q
      ipos = m + (p - 1)*(2*m - p)/2 + (q - p)
   end function gram_pair_pos

   ! dimension of the invariant coordinate vector presented to the RFF map for
   ! body order nord, under the active geometry model (and, for cpip, chemistry
   ! mode: mode 2 appends the colored moments to the RFF input).
   ! For spip the orbit tables must already be built (init_ftnbody_spip_orbits).
   integer function ftnbody_dimq(nord) result(dq)
      integer, intent(in) :: nord
      select case (ftnbody_model_id)
      case (FTNBODY_MODEL_CPIP)
         dq = CPIP_MG_N(nord)
         if (ftnbody_chem_mode == 2) dq = dq + CPIP_MC_N(nord)
      case (FTNBODY_MODEL_SPIP)
         dq = spip_tab(nord)%norb
      case default          ! poly / gramm compact coordinates: d_n = n(n-1)/2
         dq = nord*(nord - 1)/2
      end select
   end function ftnbody_dimq

end module module_ftnbody

module module_ace_desc
   USE module_kind_variables, ONLY: kind_double
   integer :: ace_numax 
   integer :: ace_dim
   integer :: ace_chem ! type of ace  0 first incomplet chemical ace, 1 Aruo's from Ralf version standard version, 2 TS version 
                       ! 1 - default 
   integer, parameter :: ACE_CHEM_INCOMPLETE = 0, ACE_CHEM_STANDARD = 1, ACE_CHEM_TS = 2   
   integer :: pos_ace_chem ! position of the chemical in the list of ACEs, 0 for incomplete, ace_dim for standard, 3 for TS.                 
   integer, parameter :: NU_LIMIT_MAX = 10
   logical, dimension(:), allocatable :: l_ace_order
   logical :: l_ace_set_rcut
   logical :: init_mode_ace_radial  ! precompute radial functions 

   character(len=5000) :: ace_nmax_list, ace_lmax_list, ace_kmax_list, ace_lambda_list, & 
      ace_rcut_in_list, ace_rcut_out_list, ace_rcut_width_in_list, ace_rcut_width_out_list

   integer, dimension(:), allocatable :: ace_lmax 
   integer, dimension(:), allocatable :: ace_nmax 
   integer, dimension(:), allocatable :: ace_kmax 
   real(kind_double), dimension(:), allocatable :: ace_rcut_in, ace_rcut_out, ace_rcut_width_in, ace_rcut_width_out, ace_lambda 
   integer :: acenmax, acelmax, acekmax ! nmax and lmax for all body order together.   
   integer :: ace_radial_poly 
   integer :: ace_radial_chem 
   integer :: ace_npoints_spline
   integer :: ace_chem_low_rank       ! 0 = disabled, 1 = enable low-rank chemical compression
   integer :: ace_chem_low_rank_q     ! compression rank Q (default 4)
   integer :: ace_chem_low_rank_niter ! number of ALS iterations (default 50)
   real(kind_double) :: ace_chem_low_rank_lambda ! ridge regularization (default 1.0d-6)
   integer :: ace_svd_randomized           ! 0 = exact dgesvd (default), 1 = randomized (Halko/Martinsson-Rokhlin) truncated SVD for the HSVD/BLOCK_HSVD radial (S^2 SVDs -> big speedup for many species)
   integer :: ace_svd_randomized_oversample ! oversampling p added to kmax when sketching (default 10)
   integer :: ace_svd_randomized_power_iter ! number of power iterations to sharpen the spectrum (default 2)
   integer :: zetaace_order 
   integer, parameter :: MAX_ZETAACE_ORDER = 6 
   character(len=100), dimension(MAX_ZETAACE_ORDER) :: dim_delta_zetaace 
   ! the dimension of the zetaace_order x ace_numax 
   real(kind_double), dimension(:,:), allocatable   :: delta_zetaace 

   ! the way in which the CG is generated: 1 for draft way 2 Dusson-Ortner method. Default 2. 
   integer :: ace_gencg
   integer, parameter :: GENCG_DRAFT = 1, GENCG_DUSORT = 2
   integer :: nkmax_order1, nkmin_order1
   real(kind_double), dimension(:,:), allocatable :: delta_ace 
   real(kind_double) :: time_amat, time_amat_00, time_amat_01, time_amat_02, time_amat_03, time_cgord , & 
                        time_ylm, time_rad, time_unique 

   type basis_ace
      integer  :: nmax           ! maximum order of radial 
      integer  :: kmax           ! maximum order of radial functions
      integer  :: lmax           ! maximum angular momentum   
      integer  :: mumax          ! maximum number of species    
      logical  :: active         ! if the basis is active or not 
      real(kind_double) :: r_cut_in, r_cut_out, r_cut_width_in, r_cut_width_out, lambda 
   end type basis_ace

   type(basis_ace), dimension(:), allocatable  :: base_params 

end module module_ace_desc 

module module_kernel_zetabody 
   USE module_kind_variables, ONLY: kind_double
   integer, parameter :: MAX_ZETABODY_ORDER = 6
   character(len=100), dimension(MAX_ZETABODY_ORDER) :: dim_grid_zetabody, dim_length_zetabody, dim_delta_zetabody
   integer :: zetabody_dim, zeta
   integer :: zetabody_order 
   real(kind_double) :: time_for_desc_kernel_nb, tnn_nb, tuu_nb, tkk_nb
   real(kind_double) :: r_cut_z2b, r_cut_width_z2b
   real(kind_double) :: r_cut_z3b, r_cut_width_z3b
   integer, dimension(:,:), allocatable    :: ispec2 
   integer, dimension(:,:,:), allocatable  :: ispec3 

   real(kind_double), dimension(:), allocatable   :: delta_zetabody2, delta_zetabody3 
   real(kind_double)  :: length_zetabody2
   real(kind_double), dimension(3)  :: length_zetabody3 

   real(kind_double), dimension(:), allocatable :: tmp_zeta_real
   real(kind_double), dimension(:,:,:), allocatable :: tmp_zeta_dreal
   real(kind_double)    :: tmp_zeta_norm 
  
end module module_kernel_zetabody 

module module_body_desc
   USE module_kind_variables, ONLY: kind_double

   integer  :: bond_dist_transform
   integer, dimension(6)      :: dim_desc_body
   logical, dimension(6)      :: l_body_order
   integer, parameter   :: bond_dist_pure = 1, &
      bond_dist_exp = 2, &
      bond_dist_inverse = 3
   real(kind_double)    :: bond_beta, bond_dist_ann
   integer, dimension(5)      :: body_D_max
   integer, dimension(1)      :: deg_p1, deg_s1, deg_p2, deg_s2, deg_s3
   integer, dimension(3)      :: deg_p3
   integer, dimension(6)      :: deg_p4
   integer, dimension(12)     :: deg_s4
   real(kind(1.d0))     :: time_body(2:5), time_invariants(2:5), time_polynomial(2:5)
   integer, dimension(:), allocatable     :: bb4
   integer, target, dimension(:, :), allocatable      :: kk4

   logical  :: desc_forces_bond

   !integer, dimension(:), allocatable :: i_type
   !for desc
   real(kind(1.d0)), dimension(:), allocatable  :: tmp_real
   !for desc forces
   !real(kind(1.d0)), dimension(:,:), allocatable :: d_ur_central, d_r_fcut,  d_r_central tmpcos_dxp, tmp_dxp
   !real(kind(1.d0)), dimension(:,:), allocatable ::  tmpcos_dxp, tmp_dxp
   real(kind(1.d0)), dimension(:, :, :), allocatable  :: tmp_force

   type basis_body
      integer  :: n           ! body order
      integer  :: Dmax
      integer  :: b           ! the b from the secondary s_{n,b}
      integer  :: Dcurrent
      integer  :: i_in_body
      integer, dimension(:), allocatable     :: ka
   end type basis_body

   type(basis_body), dimension(:), allocatable  :: basis, basis2, basis3, basis4

end module module_body_desc


module module_nlinear
   use module_kind_variables, only: kind_double
   integer  :: order_nlinear, dim_allocate_nl
   logical  :: lnlinear_precond

   ! \eplsilon_a   -> epsa (Natoms)
   real(kind_double), dimension(:), allocatable  :: epsa
   ! \eplasion_a^n -> epsa_order(n=2:5, Natoms)
   real(kind_double), dimension(:, :), allocatable     :: epsa_order
   ! w_params_order -> (n=2:5, dim_desc)
   real(kind_double), dimension(:, :), allocatable     :: w_params_nlinear

   ! \grad_{alpha, k} \eplsilon_a  -> depsa(3, inn, Natoms)
   real(kind_double), dimension(:, :, :), allocatable  :: dw_epsa

   type nl_params
      real(kind_double), dimension(:, :), allocatable     :: parameters
   end type nl_params

   type(nl_params), dimension(2:5)  :: nl_wparams

   logical, dimension(2:5) :: l_nl_order
   integer, dimension(2:5), parameter :: nl_order = (/2,3,4,5/)
   real(kind_double), dimension(:, :), allocatable     :: alpha_nl
   real(kind_double), dimension(:, :), allocatable     :: beta_nl
   integer :: nparams_nlinear
end module module_nlinear

module module_kernel_2b
   USE module_kind_variables, ONLY: kind_double
   logical :: activate_k2b
   integer :: dim_kernel_2b   ! final dimension of 2b kernel
   integer :: np_radial_2b    ! number of point for radial distance
   integer :: n_types_pair_2b  ! number of pairs by species
   integer :: sparse_points_2b_build
   integer, dimension(:,:), allocatable :: delta_type
   real(kind_double) :: sigma_2b, delta_2b
   real(kind_double), dimension(:), allocatable :: zr_2b, zr_fcut_2b
   real(kind_double) :: time_for_desc_kernel_2b, train_time_for_desc_kernel_2b, &
      test_time_for_desc_kernel_2b, &
      tnn_2b, tuu_2b, tkk_2b

   real(kind_double) :: r_cut_2b, r_cut_width_2b

   type z_kernel_2b_points
      integer :: type1_z, type2_z
      integer :: irad
      real(kind_double) :: rr
   end type z_kernel_2b_points

   type(z_kernel_2b_points), dimension(:), allocatable :: zpoints_2b 

end module module_kernel_2b

module module_covariance
  use module_kind_variables, ONLY: kind_double
  logical :: train_covariance_matrix 
  real(kind_double), dimension(:,:), allocatable :: scovar_mat
end module module_covariance

module module_kernel
   USE module_kind_variables, ONLY: kind_double
   integer  :: dim_kernel

   logical  :: write_kernel_matrix

   real(kind_double)    :: kernel_power

   integer  :: kernel_type, krff_type  ! kernel type, one of those below ...
   integer, parameter   :: kernel_se = 1, &        ! square-exp - Gaussian
      kernel_ou = 2, &        ! ornstein_uhlenbeck
      kernel_mc = 3, &        ! matern_class
      kernel_po = 4, &        ! polynomial
      kernel_po_scaled = 5, & ! polynomial scaled
      kernel_maha = 6, &      ! by MAHA
      kernel_random_maha = 66, &     ! by MAHA Fourier
      kernel_random = 7, &    ! random Freature Fourier square-exp Gaussian 
      kernel_random_po = 44   ! polynomial random

   integer, parameter  ::  krff_gaussian = 1, &
      krff_cauchy = 2,  &
      krff_laplace = 3

   integer  :: kernel_dump, np_kernel_ref, np_kernel_full, np_omega
   integer, parameter   :: kernel_dump_by_mahalanobis_norm = 1, &
      kernel_dump_by_mahalanobis = 2, &
      kernel_dump_by_cur = 3, &
      kernel_dump_by_cur_maha = 4

   ! parameters of the kernel ...
   real(kind_double)    :: time_for_desc_kernel
   real(kind_double)    :: train_time_for_desc_kernel
   real(kind_double)    :: test_time_for_desc_kernel
   real(kind_double)    :: length_kernel
   real(kind_double)    :: sigma_kernel
   real(kind_double), parameter     :: l_kou = 1.d0
   !where the selected vectors are stores. To be ScaLpack in near future ...
   !if of dimension dim_descriptor x dim_kernel = (dim_desc, dim_kernel)

   real(kind_double)    :: power_mcd, norm_random_maha 
   real(kind_double), dimension(:, :), allocatable    :: global_kernel, draft_kernel
   real(kind_double), dimension(:, :), allocatable    :: maha_kernel
   real(kind_double), dimension(:), allocatable :: maha_norm_kernel

   !tenporary, vectors ...
   real(kind_double), dimension(:), allocatable :: tmp_kernel, min_ker, max_ker, mean_ker, var_ker

   real(kind_double), dimension(:), allocatable :: kernel_phase_random
   type kernel_state
      integer  :: ia          ! The iath atom from iconf
      integer  :: iconf       ! iconf order in the database
      character(len=15)    :: filename                 ! thename of the file
   end type kernel_state

   type(kernel_state), dimension(:), allocatable      :: info_kernel, info_kernel_full, &
      info_kernel_ref, info_kernel_copy

   type basis_random
      ! desc dimension
      integer  :: dim_xdesc
      ! dim of omega sampling
      integer  :: dim_omega
      ! omega's matrix dim_xdesc x dim_omega
      real(kind_double), dimension(:, :), allocatable    :: omega
   end type basis_random
   type(basis_random), dimension(:), allocatable      :: basis_random_po

   !AD length kernel for length computation
   real(kind_double), dimension(:), allocatable :: length_kernel_k, sigma_kse2_k, soverl_k
   logical :: multiple_kernel_length = .false.
end module module_kernel

module module_zbl
   use module_kind_variables, only: kind_double
   logical :: zbl_potential
   integer :: zbl_type
   integer, parameter :: zbl_mode_default_k2b = 1  ! k2b bridge (default)
   integer, parameter :: zbl_mode_alone       = 2  ! additive E_ZBL + E_MB
   integer :: type_rac_zbl
   integer, parameter :: rac_zbl_exp=1, rac_zbl_poly=2
   real(kind_double), parameter ::   epsilon0 = 55.26349406d-4 !in e^2/(eV Ang)
   real(kind_double), parameter :: pz_phi = 0.46848d0
   real(kind_double), parameter :: p1_phi=0.32825d0, e1_phi=-2.54931d0, &
      p2_phi=0.09219d0, e2_phi=-0.29182d0, &
      p3_phi=0.58110d0, e3_phi=-0.59231d0

   real(kind_double), parameter :: zz_exp = 0.23d0
   real(kind_double) :: r1_zbl
   real(kind_double) :: r2_zbl

   real(kind_double) :: rr_k2b
   real(kind_double), dimension(:,:), allocatable :: params_k2b_to_zbl

   ! pair-specific ZBL boundaries: r1_zbl_pair(i,j), r2_zbl_pair(i,j)
   ! dimensions: (fix_no_of_elements, fix_no_of_elements)
   real(kind_double), dimension(:,:), allocatable :: r1_zbl_pair, r2_zbl_pair
   ! pair-specific rr_k2b (= r_cut_pair_in for bridge mode)
   real(kind_double), dimension(:,:), allocatable :: rr_k2b_pair

end module module_zbl

module module_variance
   use module_kind_variables, only: kind_double
   logical :: force_variance, energy_variance, stress_variance
   real(kind_double) :: sigmam_e
   real(kind_double) :: sigmam_f
   real(kind_double) :: sigmam_s
   real(kind_double), dimension(:, :), allocatable      :: VInv
end module module_variance

module module_optimization
   use module_kind_variables, ONLY: kind_double
   integer :: itopt
   real(kind_double)    :: lambda_krr_fake = 1.d-8
   logical  :: optimize_weights_db, optimize_weights_chem
   logical :: optimize_weights_L1, optimize_weights_L2, optimize_weights_Le
   integer  :: optimize_ga_population
   character(len=80)    :: class_no_optimize_weights
   character(len=2), dimension(:), allocatable  :: no_class_weights
   integer  :: max_iter_optimize_weights
   real(kind_double)    :: factor_force_error, factor_energy_error, factor_stress_error

   ! dimensions for function in genetical algo
   integer  :: dim_weights_function, dim_weights_function_full, &
      dim_chemical_function

   real(kind(0.d0)), dimension(:), allocatable  :: tmp_weights
   ! reduced dimension of the objective function arguments ...
   real(kind(0.d0)), dimension(:), allocatable  :: lower_weights, upper_weights

   real(kind(0.d0)), dimension(:), allocatable  :: tmp_chemical
   ! reduced dimension of the objective function arguments ...
   real(kind(0.d0)), dimension(:), allocatable  :: lower_chemical, upper_chemical


end module module_optimization


module module_db_poscar
   integer  :: iread_energy
   character(len=180)   :: ref_energy_per_element
   real(kind=kind(0.d0)), dimension(:), allocatable   :: fix_ref_energy_per_element

   integer, dimension(:), allocatable ::  iconf_start_on_proc, iconf_end_on_proc, &
      iconf_size_on_proc, iconf_to_proc
   integer :: i_start_conf, i_final_conf
   integer :: procs_per_file
end module module_db_poscar

module module_db_setup
   use module_kind_variables, only: kind_double
   integer  :: iconf_data, iconf_data_train, iconf_data_test
   integer, dimension(:), allocatable     :: db_train, db_test

   character(:), allocatable  :: db_file, db_path
   integer  :: selection_type
   integer, parameter   :: selection_type_first = 1, &
      selection_type_last = 2, &
      selection_type_random = 3, &
      selection_type_first_start = 4
   !md_iconf is the corresponding iconf for molecular dynamics and is 1
   integer  :: md_iconf

   ! Interatomic-distance-based config filtering, applied uniformly to both
   ! the POSCAR and the XYZ/extxyz database-reading paths.
   ! drop_short_dist > 0  : threshold (Angstrom); any config containing an
   !                        atom pair closer than this is dropped from the
   !                        database (before random selection, so no gaps).
   ! drop_short_dist <= 0 : disabled, no config is ever dropped (default).
   real(kind_double)  :: drop_short_dist

end module module_db_setup


module module_error_by_class
   use module_kind_variables, ONLY: kind_double
   integer :: number_of_distinct_classes
   integer, allocatable, dimension(:) :: class_if_ene, class_if_for, class_if_str
   real(kind_double), allocatable, dimension(:) :: rmse_by_class_ene, mae_by_class_ene, &
      corr_by_class_ene, detr_by_class_ene
   real(kind_double), allocatable, dimension(:) :: rmse_by_class_for, mae_by_class_for, &
      corr_by_class_for, detr_by_class_for
   real(kind_double), allocatable, dimension(:) :: rmse_by_class_str, mae_by_class_str, &
      corr_by_class_str, detr_by_class_str

end module module_error_by_class

module module_chemical_species
   use module_kind_variables, only: kind_double
   integer, parameter   :: size_periodic_table = 109
   type element_periodic_table
      integer  :: Z
      character(len=2)     :: symbol
      character(len=50)    :: name, elec_conf
      real(kind_double)    :: mass, frozenp, meltp, density, ionization, covalent_radius
   end type element_periodic_table
   type(element_periodic_table), dimension(size_periodic_table)   :: periodic_table_element

   ! input variable
   integer  :: fix_no_of_elements
   !input variable
   character(len=5000)    ::  chemical_elements, weight_per_element, weight_per_element_3ch

   !input variable for Arnaud
   logical :: img_weighted
   integer :: img_num_ch
   character(len=5000) :: ws1, ws2, ws3, ws4
   real(kind_double), dimension(:,:), allocatable :: fix_wspecies
   real(kind_double) :: tnn_rdist, tcc01_rdist, tcc02_rdist, tcc03_rdist, tii_rdist
   ! end input variable for Arnaud.

   integer, dimension(:), allocatable     :: fix_type_to_periodic
   character(len=2), dimension(:), allocatable  :: fix_ch_elements
   real(kind=kind(0.d0)), dimension(:), allocatable   :: fix_Z_elements, fix_mass_elements, fix_covalent_radius_elements

   integer, dimension(:,:), allocatable :: map_species_half, map_species_full
   integer ::  size_species_full, size_species_half
end module module_chemical_species


module module_input_als_fit
!-------------------------------------------------------------------------------
! module_input_als_fit
!
! Input parameters for the ALS-Ridge block preconditioning fitting method
! (mld_fit_type = 5).  See MiladyNoteTechnique5.pdf, sections 1.10.1 and 1.10.2.
!
! This module is intentionally separated from ml_in_ndm_module so that the
! ALS-specific variables do not clutter the main parameter module.
! The fit_als constant itself stays in ml_in_ndm_module alongside the other
! fit_* enum values.
!-------------------------------------------------------------------------------

   use module_kind_variables, only: kind_double
   implicit none

   integer  :: als_nsteps          ! S: number of ALS outer iterations
   real(kind_double) :: als_tol    ! convergence tolerance for w and alpha
   real(kind_double) :: als_ridge_k  ! k > 0: global ridge parameter
   real(kind_double) :: als_rho    ! rho: regularisation for the alpha sub-problem
   integer  :: als_alpha_method    ! 0 = fixed alpha (sec 1.10.1), 1 = learning alpha (sec 1.10.2)
   integer, parameter :: als_alpha_fixed = 0, als_alpha_learning = 1
   integer  :: als_precond_type    ! 0 = Frobenius norm, 1 = SVD-based, 2 = flat (lambda=lambda_krr)
   integer, parameter :: als_precond_frobenius = 0, als_precond_svd = 1, als_precond_flat = 2
   logical  :: als_nnls_alpha      ! enforce alpha >= 0 via NNLS in the alpha sub-problem
   integer  :: als_nnls_mode       ! 0 = Lawson-Hanson (default), 1 = BK warm-start NNLS
   integer, parameter :: als_nnls_lawson_hanson = 0, als_nnls_bk_warm = 1
   integer  :: als_nu_max          ! number of blocks (body orders)
   ! block partition vector of size 2*als_nu_max:
   !   als_block_partition(2*i-1) = row start of block i
   !   als_block_partition(2*i)   = row end   of block i
   !   i = 1, ..., als_nu_max
   ! (indices refer to rows of the design matrix A, i.e. descriptor components)
   integer, dimension(:), allocatable :: als_block_partition

   ! Output: learned block renormalisation coefficients (size als_nu_max)
   ! Filled after training; w_params already contains alpha*w (the effective
   ! coefficients for extrapolation).  als_alpha stores the raw alpha_nu
   ! separately for diagnostics, logging and restart.
   real(kind_double), dimension(:), allocatable :: als_alpha

contains

!-------------------------------------------------------------------------------
! set_als_block_partition
!
! Compute als_block_partition and als_nu_max from the descriptor structure.
!
! Three modes:
!   (a) User-supplied partition — als_block_partition was already filled
!       (first entry /= 0):  return immediately.
!
!   (b) Default mode — nblocks == 0 (no per-block dims supplied):
!       single block [1, D_total].  als_nu_max = 1.
!
!   (c) Descriptor-driven partition — nblocks > 0, block_dims(:) supplied:
!       each block i spans a contiguous row range of the design matrix.
!       als_nu_max = nblocks, partition = cumulative sum of block_dims.
!
! The caller is responsible for assembling `block_dims` from the
! descriptor modules (module_body_desc or module_ace_desc / module_base_cnlm).
!
! Arguments:
!   D_total     — total number of descriptor rows (dim_xdesc or dimr_sca_Amat)
!   nblocks     — number of descriptor blocks (0 = default single block)
!   block_dims  — integer array(nblocks) with the dimension of each block.
!                 SUM(block_dims) must equal D_total.
!                 Ignored when nblocks == 0.
!
! On return:
!   als_nu_max          is set (or confirmed)
!   als_block_partition is (re-)allocated and filled
!-------------------------------------------------------------------------------
subroutine set_als_block_partition(D_total, nblocks, block_dims)

   implicit none
   integer, intent(in) :: D_total
   integer, intent(in) :: nblocks
   integer, intent(in), optional :: block_dims(:)

   integer :: inu, ipos, dsum

   ! (a) User explicitly provided the partition — keep it
   if (allocated(als_block_partition)) then
     if (size(als_block_partition) >= 2) then
       if (als_block_partition(1) /= 0) return
     end if
   end if

   ! (b) Default mode: single block covering all rows
   if (nblocks <= 0 .or. (.not. present(block_dims))) then
     als_nu_max = 1
     if (allocated(als_block_partition)) deallocate(als_block_partition)
     allocate(als_block_partition(2))
     als_block_partition(1) = 1
     als_block_partition(2) = D_total
     return
   end if

   ! (c) Descriptor-driven partition
   ! Sanity: check that the block dimensions sum to D_total
   dsum = 0
   do inu = 1, nblocks
     dsum = dsum + block_dims(inu)
   end do
   if (dsum /= D_total) then
     ! Mismatch — fall back to single block and warn
     write(6, '("WARNING set_als_block_partition: SUM(block_dims)=", i10, &
              &" /= D_total=", i10, ". Falling back to single block.")') dsum, D_total
     als_nu_max = 1
     if (allocated(als_block_partition)) deallocate(als_block_partition)
     allocate(als_block_partition(2))
     als_block_partition(1) = 1
     als_block_partition(2) = D_total
     return
   end if

   als_nu_max = nblocks
   if (allocated(als_block_partition)) deallocate(als_block_partition)
   allocate(als_block_partition(2 * als_nu_max))

   ipos = 1
   do inu = 1, als_nu_max
     als_block_partition(2*inu - 1) = ipos
     als_block_partition(2*inu)     = ipos + block_dims(inu) - 1
     ipos = ipos + block_dims(inu)
   end do

end subroutine set_als_block_partition

end module module_input_als_fit


module ml_in_ndm_module

   USE module_kind_variables, ONLY: kind_double
#ifdef MLD_NDM   
   use gen_com_m, only: rangml
   use gen_com_m_ml, only: im, imm 
#else
   use ondm_gen_com_m, only: rangml, im, imm 
#endif

#if(PARA)
   use mpi
   !use mld_mpi 
#endif

   implicit none

   real(kind_double), parameter     :: two_pi = 8.d0*atan(1.d0)
   real(kind_double), parameter     :: one_pi = 4.d0*atan(1.d0)
   ! real(kind_double), parameter     ::  pi = 3.141592654D0
   real(kind_double), dimension(3), parameter   :: vec_3d_zero = (/0.d0, 0.d0, 0.d0/)
   real(kind_double), parameter     :: sqrt_two = sqrt(2.d0)
   double complex, parameter  :: c_zero = (0.d0, 0.d0)

   integer  :: n_elements
   integer, dimension(:), allocatable     :: z_elements

   integer  :: n_pca

   integer  :: i_start_at, i_final_at
   real(kind_double), dimension(:, :), allocatable    :: design_mat          ! (D,M) D size of the database D values of x,
   !  M the dimension of fingerprint ...
   real(kind_double)    :: eta_max_g2, eta_min_g2
   real(kind_double)    :: rs_max_g2, rs_min_g2

   integer  :: pow_so4_dim
   
   integer  ::  l_max

   integer  :: mtp_dim

   integer  :: rmat_dim, dmilady_dim
   real(kind_double)    :: power_line, power_coeff_renorm
   !dmilady real(kind_double), allocatable, dimension(:, :)    :: dmilady
   !dmilady real(kind_double), allocatable, dimension(:, :, :, :)    :: dmilady_deriv

   integer  :: soap_dim, nspecies_soap
   real(kind_double), allocatable, dimension(:) :: rb_soap
   real(kind_double), allocatable, dimension(:, :)    ::  W_soap, S_factor_matrix
   integer, allocatable, dimension(:, :)  :: ns_soap_index
   integer  :: n_soap
   logical  :: lsoap, lsoap_fcut_wes
   logical  :: lsoap_diag, lsoap_norm, lsoap_lnorm
   double precision     :: alpha_soap, atom_sigma_soap, r_cut_width_soap
   real(kind_double), allocatable, dimension(:, :)    :: k_soap
   real(kind_double), allocatable, dimension(:, :, :) :: r_soap, at_soap
   integer, allocatable, dimension(:, :)  :: iwmax2_soap, indi2_soap

   real     :: j_max
   integer  :: jj_max

   real(kind_double)    :: svd_rcond


   double precision, allocatable, dimension(:, :, :, :, :, :)     :: cg_vector

   real(kind_double)    :: r_factorial(0:167)
   ! new version for g2 and g3
   integer  :: n_g2_eta, n_g2_rs, g2_dim
   real(kind_double), dimension(:), allocatable :: g2_eta, g2_rs

   integer  :: n_g3_eta, n_g3_zeta, n_g3_lambda, g3_dim
   real(kind_double), dimension(:), allocatable :: g3_eta, g3_zeta, g3_lambda
   integer  :: behler_dim

   integer              :: imm_neigh = 300
   integer  :: nd_data     ! dimension of the database M, internal
   integer  :: nd_fingerprint                       ! dimension of the feature space D, internal
   integer  :: ml_type     ! input
   integer  :: iread_ml    ! input
   integer  :: isave_ml    ! input
   integer  :: dim_kxx
   integer  :: descriptor_type                      ! input
   integer, parameter   :: descriptor_g2 = 1, &
      descriptor_g3 = 2, &
      descriptor_behler = 3, &
      descriptor_afs = 4, &
      descriptor_pow_so3 = 6, &
      descriptor_pow_so3_3body = 603, &
      descriptor_pow_so4 = 8, &
      descriptor_g2_pow_so4 = 18, &
      descriptor_bispectrum_so4 = 9, &
      descriptor_g2_afs = 14, &
      descriptor_g2_bispectrum_so4 = 19, &
      descriptor_milady = 77, &
      descriptor_mtp = 100, &
      descriptor_body = 200, &
      descriptor_ftnbody = 202, &
      descriptor_tbind = 203, &
      descriptor_ace = 300, & 
      descriptor_zetabody = 204  ! this is a real 2b 3b description of the system

   integer, parameter   :: ml_type_krr = 1, &       ! kernel versio
      ml_type_krr_old = 11, & ! this is to keep compatibility with old stuff.
      ml_type_basis = 0, &
      ml_type_nlinear = 3, &  ! n-linear version
      ml_type_gap = 2, &
      ml_type_descriptors = -1, &
      ml_type_analysis = -2

   integer  ::  seed
   real(kind_double)    :: lambda_krr, lambda_krr_2, min_lambda_krr, max_lambda_krr
   integer  :: n_values_lambda_krr
   character(len=4)     :: regularization_name
   real(kind_double), dimension(:), allocatable :: vector_lambda_krr
   logical  :: toy_model, debug, weighted_auto, weighted, weighted_3ch
   logical  :: write_desc, write_desc_dump, read_desc_dump, write_design_matrix, write_test_design_matrix
   integer :: desc_file_format
   integer, parameter ::  eml_type=1, &  ! the basic mode forced is  write_desc_dump, read_desc_dump are true
      csv_type=2, &
      npz_type=3, &
      hdf_type=4 
   logical :: hdf_positions   
   logical  :: marginal_likelihood, desc_forces

   !database related
   ! iconf_data = iconf_data_train + iconf_data_test

   character(len=4)     :: char_desc


   integer :: mld_dmtype 
   integer, parameter :: MD_MLD_DMTYPE = 181,  & 
                         ML_MLD_DMTYPE = 18 

   integer  :: mld_order
   integer, parameter   :: mld_linear = 1, &
                           mld_linear_extended = 11, &
                           mld_quadratic = 2, &
                           mld_polyc = 3, &
                           mld_kernel = 7

   integer  :: mld_type_quadratic
   integer, parameter   :: mld_type_quadratic_qnml = 1, &
                           mld_type_quadratic_qml = 2, &
                           mld_type_quadratic_bilinear = 3, & 
                           mld_type_quadratic_zaxa = 4, &
                           mld_type_quadratic_ZX = 5
   ! how to fit ...
   integer  :: mld_fit_type
   integer, parameter   :: fit_home_hb = 0, &
                           fit_lapack_qr = 1, &
                           fit_lapack_qr_constraints = 2, &
                           fit_lapack_ortho = 3, &
                           fit_lapack_svd = 4, &
                           fit_als = 5, &
                           fit_online_svd = 6
   ! number of training configs whose descriptor columns are batched
   ! together (per subworld group) before one collective broadcast +
   ! local DGEMM update of the ScaLAPACK Gram accumulator (fit_online_svd,
   ! see module_online_fit.F90). Amortizes the fixed per-collective
   ! latency/synchronization cost over more data. 1 = original
   ! per-config behaviour. Larger values reduce the number of
   ! world-wide MPI_Bcast calls proportionally, at the cost of a
   ! larger temporary D x (batch worth of observations) buffer on
   ! every rank -- keep conservative on memory-constrained allocations.
   integer  :: online_fit_batch_size

   logical  :: strict_behler

   real(kind(0.d0))     :: sign_stress, sign_stress_big_box

   !fit_home_regularization = 10
   ! regularization of not ...
   integer  :: mld_regularization_type
   integer, parameter   :: mld_regularization_type_home = 1
   character(len=2)     :: snap_class_constraints
   logical  :: train_only, train_time = .true.

   integer  :: mtp_poly_min, mtp_poly_max, mtp_rad_order
   real(kind_double)    :: tmp_val_desc_max

   integer  :: polyc_n_poly, polyc_n_hermite

   ! weights optimization
   integer :: type_of_loss
   integer, parameter :: loss_init = 1 , & ! the initial loss function
      loss_fair = 2 , &
      loss_per = 3


   real(kind=kind(0.d0)), dimension(:), allocatable   :: fix_weighted_for_element, fix_weighted_for_element_3ch
   real(kind=kind(0.d0)), dimension(:), allocatable   :: fix_weighted_for_element_ini, fix_weighted_for_element_3ch_ini


   integer  :: fix_no_of_elements_invisible
   character(len=80)    :: chemical_elements_invisible
   CHARACTER(len=2), dimension(:), allocatable  :: fix_ch_elements_invisible
   real(kind=kind(0.d0)), dimension(:), allocatable   :: fix_weighted_for_element_invisble
   integer, dimension(:), allocatable     :: fix_type_to_periodic_invisible
   real(kind=kind(0.d0)), dimension(:), allocatable   :: fix_Z_elements_invisible, fix_mass_elements_invisible
   real(kind_double)    :: rvois_ndm_ml, renorm_cov, renorm_mass


   logical  :: linvisible, lmask, mask_file

   character(len=2), dimension(:), allocatable  :: classes_train_for_sigma, classes_full_for_sigma
   character(len=80)    :: classes_for_mcd
   character(len=80)    :: dim_fourier_nbody , length_fourier_nbody, delta_fourier_nbody

   !fixed number of neighbour
   logical :: Nfix
   integer :: discrete_fix_N_rcut, fix_Nmax_neigh
   real(kind_double) :: delta_fix_N_rcut

#if(PARA)
   integer  :: iproc, nb_elements, itempproc
#endif

contains




   subroutine rotate(a, x, y)
      double precision, dimension(3), intent(in)   :: a
      double precision, dimension(3), intent(in)   :: x
      double precision, dimension(3), intent(out)  :: y
      double precision, dimension(3, 3)      :: rotx, rotz1, rotz2

      call mat_rotx(a(2)*one_pi, rotx)
      call mat_rotz(a(1)*2*one_pi, rotz1)
      call mat_rotz(a(3)*2*one_pi, rotz2)
      y = matmul(rotz2, matmul(rotx, matmul(rotz1, x)))
   end subroutine rotate



   subroutine mat_rotx(angle, rotx)
      double precision, intent(in)     :: angle
      double precision, dimension(3, 3), intent(out)     :: rotx

      rotx(1, 1) = 1
      rotx(1, 2) = 0
      rotx(1, 3) = 0
      rotx(2, 1) = 0
      rotx(2, 2) = cos(angle)
      rotx(2, 3) = -sin(angle)
      rotx(3, 1) = 0
      rotx(3, 2) = sin(angle)
      rotx(3, 3) = cos(angle)
   end subroutine mat_rotx



   subroutine mat_rotz(angle, rotz)
      double precision, intent(in)     :: angle
      double precision, dimension(3, 3), intent(out)     :: rotz

      rotz(1, 1) = cos(angle)
      rotz(1, 2) = -sin(angle)
      rotz(1, 3) = 0
      rotz(2, 1) = sin(angle)
      rotz(2, 2) = cos(angle)
      rotz(2, 3) = 0
      rotz(3, 1) = 0
      rotz(3, 2) = 0
      rotz(3, 3) = 1
   end subroutine mat_rotz



   subroutine prepare_factorial()
      r_factorial(0:167) = (/ &
         1.d0, &
         1.d0, &
         2.d0, &
         6.d0, &
         24.d0, &
         120.d0, &
         720.d0, &
         5040.d0, &
         40320.d0, &
         362880.d0, &
         3628800.d0, &
         39916800.d0, &
         479001600.d0, &
         6227020800.d0, &
         87178291200.d0, &
         1307674368000.d0, &
         20922789888000.d0, &
         355687428096000.d0, &
         6.402373705728d+15, &
         1.21645100408832d+17, &
         2.43290200817664d+18, &
         5.10909421717094d+19, &
         1.12400072777761d+21, &
         2.58520167388850d+22, &
         6.20448401733239d+23, &
         1.55112100433310d+25, &
         4.03291461126606d+26, &
         1.08888694504184d+28, &
         3.04888344611714d+29, &
         8.84176199373970d+30, &
         2.65252859812191d+32, &
         8.22283865417792d+33, &
         2.63130836933694d+35, &
         8.68331761881189d+36, &
         2.95232799039604d+38, &
         1.03331479663861d+40, &
         3.71993326789901d+41, &
         1.37637530912263d+43, &
         5.23022617466601d+44, &
         2.03978820811974d+46, &
         8.15915283247898d+47, &
         3.34525266131638d+49, &
         1.40500611775288d+51, &
         6.04152630633738d+52, &
         2.65827157478845d+54, &
         1.19622220865480d+56, &
         5.50262215981209d+57, &
         2.58623241511168d+59, &
         1.24139155925361d+61, &
         6.08281864034268d+62, &
         3.04140932017134d+64, &
         1.55111875328738d+66, &
         8.06581751709439d+67, &
         4.27488328406003d+69, &
         2.30843697339241d+71, &
         1.26964033536583d+73, &
         7.10998587804863d+74, &
         4.05269195048772d+76, &
         2.35056133128288d+78, &
         1.38683118545690d+80, &
         8.32098711274139d+81, &
         5.07580213877225d+83, &
         3.14699732603879d+85, &
         1.98260831540444d+87, &
         1.26886932185884d+89, &
         8.24765059208247d+90, &
         5.44344939077443d+92, &
         3.64711109181887d+94, &
         2.48003554243683d+96, &
         1.71122452428141d+98, &
         1.19785716699699d+100, &
         8.50478588567862d+101, &
         6.12344583768861d+103, &
         4.47011546151268d+105, &
         3.30788544151939d+107, &
         2.48091408113954d+109, &
         1.88549470166605d+111, &
         1.45183092028286d+113, &
         1.13242811782063d+115, &
         8.94618213078297d+116, &
         7.15694570462638d+118, &
         5.79712602074737d+120, &
         4.75364333701284d+122, &
         3.94552396972066d+124, &
         3.31424013456535d+126, &
         2.81710411438055d+128, &
         2.42270953836727d+130, &
         2.10775729837953d+132, &
         1.85482642257398d+134, &
         1.65079551609085d+136, &
         1.48571596448176d+138, &
         1.35200152767840d+140, &
         1.24384140546413d+142, &
         1.15677250708164d+144, &
         1.08736615665674d+146, &
         1.03299784882391d+148, &
         9.91677934870949d+149, &
         9.61927596824821d+151, &
         9.42689044888324d+153, &
         9.33262154439441d+155, &
         9.33262154439441d+157, &
         9.42594775983835d+159, &
         9.61446671503512d+161, &
         9.90290071648618d+163, &
         1.02990167451456d+166, &
         1.08139675824029d+168, &
         1.14628056373471d+170, &
         1.22652020319614d+172, &
         1.32464181945183d+174, &
         1.44385958320249d+176, &
         1.58824554152274d+178, &
         1.76295255109024d+180, &
         1.97450685722107d+182, &
         2.23119274865981d+184, &
         2.54355973347219d+186, &
         2.92509369349301d+188, &
         3.39310868445190d+190, &
         3.96993716080872d+192, &
         4.68452584975429d+194, &
         5.57458576120760d+196, &
         6.68950291344912d+198, &
         8.09429852527344d+200, &
         9.87504420083360d+202, &
         1.21463043670253d+205, &
         1.50614174151114d+207, &
         1.88267717688893d+209, &
         2.37217324288005d+211, &
         3.01266001845766d+213, &
         3.85620482362580d+215, &
         4.97450422247729d+217, &
         6.46685548922047d+219, &
         8.47158069087882d+221, &
         1.11824865119600d+224, &
         1.48727070609069d+226, &
         1.99294274616152d+228, &
         2.69047270731805d+230, &
         3.65904288195255d+232, &
         5.01288874827499d+234, &
         6.91778647261949d+236, &
         9.61572319694109d+238, &
         1.34620124757175d+241, &
         1.89814375907617d+243, &
         2.69536413788816d+245, &
         3.85437071718007d+247, &
         5.55029383273930d+249, &
         8.04792605747199d+251, &
         1.17499720439091d+254, &
         1.72724589045464d+256, &
         2.55632391787286d+258, &
         3.80892263763057d+260, &
         5.71338395644585d+262, &
         8.62720977423323d+264, &
         1.31133588568345d+267, &
         2.00634390509568d+269, &
         3.08976961384735d+271, &
         4.78914290146339d+273, &
         7.47106292628289d+275, &
         1.17295687942641d+278, &
         1.85327186949373d+280, &
         2.94670227249504d+282, &
         4.71472363599206d+284, &
         7.59070505394721d+286, &
         1.22969421873945d+289, &
         2.00440157654530d+291, &
         3.28721858553429d+293, &
         5.42391066613159d+295, &
         9.00369170577843d+297, &
         1.50361651486500d+300/)
   end subroutine prepare_factorial


end module ml_in_ndm_module



subroutine periodic_table()
   use module_chemical_species, only: periodic_table_element
   implicit none

   periodic_table_element(1)%Z = 1
   periodic_table_element(1)%mass = 1.008
   periodic_table_element(1)%name = "   Hydrogen   "
   periodic_table_element(1)%symbol = "H "
   periodic_table_element(1)%frozenp = -259
   periodic_table_element(1)%meltp = -253
   periodic_table_element(1)%density = 0.09000
   periodic_table_element(1)%elec_conf = "                 1s1"
   periodic_table_element(1)%ionization = 13.600
   periodic_table_element(2)%Z = 2
   periodic_table_element(2)%mass = 4.003
   periodic_table_element(2)%name = "    Helium    "
   periodic_table_element(2)%symbol = "He"
   periodic_table_element(2)%frozenp = -272
   periodic_table_element(2)%meltp = -269
   periodic_table_element(2)%density = 0.18000
   periodic_table_element(2)%elec_conf = "                 1s2"
   periodic_table_element(2)%ionization = 24.590
   periodic_table_element(3)%Z = 3
   periodic_table_element(3)%mass = 6.941
   periodic_table_element(3)%name = "   Lithium    "
   periodic_table_element(3)%symbol = "Li"
   periodic_table_element(3)%frozenp = 180
   periodic_table_element(3)%meltp = 1347
   periodic_table_element(3)%density = 0.53000
   periodic_table_element(3)%elec_conf = "             [He]2s1"
   periodic_table_element(3)%ionization = 5.390
   periodic_table_element(4)%Z = 4
   periodic_table_element(4)%mass = 9.012
   periodic_table_element(4)%name = "  Beryllium   "
   periodic_table_element(4)%symbol = "Be"
   periodic_table_element(4)%frozenp = 1278
   periodic_table_element(4)%meltp = 2970
   periodic_table_element(4)%density = 1.85000
   periodic_table_element(4)%elec_conf = "             [He]2s2"
   periodic_table_element(4)%ionization = 9.320
   periodic_table_element(5)%Z = 5
   periodic_table_element(5)%mass = 10.811
   periodic_table_element(5)%name = "    Boron     "
   periodic_table_element(5)%symbol = "B "
   periodic_table_element(5)%frozenp = 2300
   periodic_table_element(5)%meltp = 2550
   periodic_table_element(5)%density = 2.34000
   periodic_table_element(5)%elec_conf = "          [He]2s22p1"
   periodic_table_element(5)%ionization = 8.300
   periodic_table_element(6)%Z = 6
   periodic_table_element(6)%mass = 12.011
   periodic_table_element(6)%name = "    Carbon    "
   periodic_table_element(6)%symbol = "C "
   periodic_table_element(6)%frozenp = 3500
   periodic_table_element(6)%meltp = 4827
   periodic_table_element(6)%density = 2.26000
   periodic_table_element(6)%elec_conf = "          [He]2s22p2"
   periodic_table_element(6)%ionization = 11.260
   periodic_table_element(7)%Z = 7
   periodic_table_element(7)%mass = 14.007
   periodic_table_element(7)%name = "   Nitrogen   "
   periodic_table_element(7)%symbol = "N "
   periodic_table_element(7)%frozenp = -210
   periodic_table_element(7)%meltp = -196
   periodic_table_element(7)%density = 1.25000
   periodic_table_element(7)%elec_conf = "          [He]2s22p3"
   periodic_table_element(7)%ionization = 14.530
   periodic_table_element(8)%Z = 8
   periodic_table_element(8)%mass = 15.999
   periodic_table_element(8)%name = "    Oxygen    "
   periodic_table_element(8)%symbol = "O "
   periodic_table_element(8)%frozenp = -218
   periodic_table_element(8)%meltp = -183
   periodic_table_element(8)%density = 1.43000
   periodic_table_element(8)%elec_conf = "          [He]2s22p4"
   periodic_table_element(8)%ionization = 13.620
   periodic_table_element(9)%Z = 9
   periodic_table_element(9)%mass = 18.998
   periodic_table_element(9)%name = "   Fluorine   "
   periodic_table_element(9)%symbol = "F "
   periodic_table_element(9)%frozenp = -220
   periodic_table_element(9)%meltp = -188
   periodic_table_element(9)%density = 1.70000
   periodic_table_element(9)%elec_conf = "          [He]2s22p5"
   periodic_table_element(9)%ionization = 17.420
   periodic_table_element(10)%Z = 10
   periodic_table_element(10)%mass = 20.180
   periodic_table_element(10)%name = "     Neon     "
   periodic_table_element(10)%symbol = "Ne"
   periodic_table_element(10)%frozenp = -249
   periodic_table_element(10)%meltp = -246
   periodic_table_element(10)%density = 0.90000
   periodic_table_element(10)%elec_conf = "          [He]2s22p6"
   periodic_table_element(10)%ionization = 21.560
   periodic_table_element(11)%Z = 11
   periodic_table_element(11)%mass = 22.990
   periodic_table_element(11)%name = "    Sodium    "
   periodic_table_element(11)%symbol = "Na"
   periodic_table_element(11)%frozenp = 98
   periodic_table_element(11)%meltp = 883
   periodic_table_element(11)%density = 0.97000
   periodic_table_element(11)%elec_conf = "             [Ne]3s1"
   periodic_table_element(11)%ionization = 5.140
   periodic_table_element(12)%Z = 12
   periodic_table_element(12)%mass = 24.305
   periodic_table_element(12)%name = "  Magnesium   "
   periodic_table_element(12)%symbol = "Mg"
   periodic_table_element(12)%frozenp = 639
   periodic_table_element(12)%meltp = 1090
   periodic_table_element(12)%density = 1.74000
   periodic_table_element(12)%elec_conf = "             [Ne]3s2"
   periodic_table_element(12)%ionization = 7.650
   periodic_table_element(13)%Z = 13
   periodic_table_element(13)%mass = 26.982
   periodic_table_element(13)%name = "   Aluminum   "
   periodic_table_element(13)%symbol = "Al"
   periodic_table_element(13)%frozenp = 660
   periodic_table_element(13)%meltp = 2467
   periodic_table_element(13)%density = 2.70000
   periodic_table_element(13)%elec_conf = "          [Ne]3s23p1"
   periodic_table_element(13)%ionization = 5.990
   periodic_table_element(14)%Z = 14
   periodic_table_element(14)%mass = 28.086
   periodic_table_element(14)%name = "   Silicon    "
   periodic_table_element(14)%symbol = "Si"
   periodic_table_element(14)%frozenp = 1410
   periodic_table_element(14)%meltp = 2355
   periodic_table_element(14)%density = 2.33000
   periodic_table_element(14)%elec_conf = "          [Ne]3s23p2"
   periodic_table_element(14)%ionization = 8.150
   periodic_table_element(15)%Z = 15
   periodic_table_element(15)%mass = 30.974
   periodic_table_element(15)%name = "  Phosphorus  "
   periodic_table_element(15)%symbol = "P "
   periodic_table_element(15)%frozenp = 44
   periodic_table_element(15)%meltp = 280
   periodic_table_element(15)%density = 1.82000
   periodic_table_element(15)%elec_conf = "          [Ne]3s23p3"
   periodic_table_element(15)%ionization = 10.490
   periodic_table_element(16)%Z = 16
   periodic_table_element(16)%mass = 32.065
   periodic_table_element(16)%name = "    Sulfur    "
   periodic_table_element(16)%symbol = "S "
   periodic_table_element(16)%frozenp = 113
   periodic_table_element(16)%meltp = 445
   periodic_table_element(16)%density = 2.07000
   periodic_table_element(16)%elec_conf = "          [Ne]3s23p4"
   periodic_table_element(16)%ionization = 10.360
   periodic_table_element(17)%Z = 17
   periodic_table_element(17)%mass = 35.453
   periodic_table_element(17)%name = "   Chlorine   "
   periodic_table_element(17)%symbol = "Cl"
   periodic_table_element(17)%frozenp = -101
   periodic_table_element(17)%meltp = -35
   periodic_table_element(17)%density = 3.21000
   periodic_table_element(17)%elec_conf = "          [Ne]3s23p5"
   periodic_table_element(17)%ionization = 12.970
   periodic_table_element(18)%Z = 18
   periodic_table_element(18)%mass = 39.948
   periodic_table_element(18)%name = "    Argon     "
   periodic_table_element(18)%symbol = "Ar"
   periodic_table_element(18)%frozenp = -189
   periodic_table_element(18)%meltp = -186
   periodic_table_element(18)%density = 1.78000
   periodic_table_element(18)%elec_conf = "          [Ne]3s23p6"
   periodic_table_element(18)%ionization = 15.760
   periodic_table_element(19)%Z = 19
   periodic_table_element(19)%mass = 39.098
   periodic_table_element(19)%name = "  Potassium   "
   periodic_table_element(19)%symbol = "K "
   periodic_table_element(19)%frozenp = 64
   periodic_table_element(19)%meltp = 774
   periodic_table_element(19)%density = 0.86000
   periodic_table_element(19)%elec_conf = "             [Ar]4s1"
   periodic_table_element(19)%ionization = 4.340
   periodic_table_element(20)%Z = 20
   periodic_table_element(20)%mass = 40.078
   periodic_table_element(20)%name = "   Calcium    "
   periodic_table_element(20)%symbol = "Ca"
   periodic_table_element(20)%frozenp = 839
   periodic_table_element(20)%meltp = 1484
   periodic_table_element(20)%density = 1.55000
   periodic_table_element(20)%elec_conf = "             [Ar]4s2"
   periodic_table_element(20)%ionization = 6.110
   periodic_table_element(21)%Z = 21
   periodic_table_element(21)%mass = 44.956
   periodic_table_element(21)%name = "   Scandium   "
   periodic_table_element(21)%symbol = "Sc"
   periodic_table_element(21)%frozenp = 1539
   periodic_table_element(21)%meltp = 2832
   periodic_table_element(21)%density = 2.99000
   periodic_table_element(21)%elec_conf = "          [Ar]3d14s2"
   periodic_table_element(21)%ionization = 6.560
   periodic_table_element(22)%Z = 22
   periodic_table_element(22)%mass = 47.867
   periodic_table_element(22)%name = "   Titanium   "
   periodic_table_element(22)%symbol = "Ti"
   periodic_table_element(22)%frozenp = 1660
   periodic_table_element(22)%meltp = 3287
   periodic_table_element(22)%density = 4.54000
   periodic_table_element(22)%elec_conf = "          [Ar]3d24s2"
   periodic_table_element(22)%ionization = 6.830
   periodic_table_element(23)%Z = 23
   periodic_table_element(23)%mass = 50.942
   periodic_table_element(23)%name = "   Vanadium   "
   periodic_table_element(23)%symbol = "V "
   periodic_table_element(23)%frozenp = 1890
   periodic_table_element(23)%meltp = 3380
   periodic_table_element(23)%density = 6.11000
   periodic_table_element(23)%elec_conf = "          [Ar]3d34s2"
   periodic_table_element(23)%ionization = 6.750
   periodic_table_element(24)%Z = 24
   periodic_table_element(24)%mass = 51.996
   periodic_table_element(24)%name = "   Chromium   "
   periodic_table_element(24)%symbol = "Cr"
   periodic_table_element(24)%frozenp = 1857
   periodic_table_element(24)%meltp = 2672
   periodic_table_element(24)%density = 7.19000
   periodic_table_element(24)%elec_conf = "          [Ar]3d54s1"
   periodic_table_element(24)%ionization = 6.770
   periodic_table_element(25)%Z = 25
   periodic_table_element(25)%mass = 54.938
   periodic_table_element(25)%name = "  Manganese   "
   periodic_table_element(25)%symbol = "Mn"
   periodic_table_element(25)%frozenp = 1245
   periodic_table_element(25)%meltp = 1962
   periodic_table_element(25)%density = 7.43000
   periodic_table_element(25)%elec_conf = "          [Ar]3d54s2"
   periodic_table_element(25)%ionization = 7.430
   periodic_table_element(26)%Z = 26
   periodic_table_element(26)%mass = 55.845
   periodic_table_element(26)%name = "     Iron     "
   periodic_table_element(26)%symbol = "Fe"
   periodic_table_element(26)%frozenp = 1535
   periodic_table_element(26)%meltp = 2750
   periodic_table_element(26)%density = 7.87000
   periodic_table_element(26)%elec_conf = "          [Ar]3d64s2"
   periodic_table_element(26)%ionization = 7.900
   periodic_table_element(27)%Z = 27
   periodic_table_element(27)%mass = 58.933
   periodic_table_element(27)%name = "    Cobalt    "
   periodic_table_element(27)%symbol = "Co"
   periodic_table_element(27)%frozenp = 1495
   periodic_table_element(27)%meltp = 2870
   periodic_table_element(27)%density = 8.90000
   periodic_table_element(27)%elec_conf = "          [Ar]3d74s2"
   periodic_table_element(27)%ionization = 7.880
   periodic_table_element(28)%Z = 28
   periodic_table_element(28)%mass = 58.693
   periodic_table_element(28)%name = "    Nickel    "
   periodic_table_element(28)%symbol = "Ni"
   periodic_table_element(28)%frozenp = 1453
   periodic_table_element(28)%meltp = 2732
   periodic_table_element(28)%density = 8.90000
   periodic_table_element(28)%elec_conf = "          [Ar]3d84s2"
   periodic_table_element(28)%ionization = 7.640
   periodic_table_element(29)%Z = 29
   periodic_table_element(29)%mass = 63.546
   periodic_table_element(29)%name = "    Copper    "
   periodic_table_element(29)%symbol = "Cu"
   periodic_table_element(29)%frozenp = 1083
   periodic_table_element(29)%meltp = 2567
   periodic_table_element(29)%density = 8.96000
   periodic_table_element(29)%elec_conf = "         [Ar]3d104s1"
   periodic_table_element(29)%ionization = 7.730
   periodic_table_element(30)%Z = 30
   periodic_table_element(30)%mass = 65.390
   periodic_table_element(30)%name = "     Zinc     "
   periodic_table_element(30)%symbol = "Zn"
   periodic_table_element(30)%frozenp = 420
   periodic_table_element(30)%meltp = 907
   periodic_table_element(30)%density = 7.13000
   periodic_table_element(30)%elec_conf = "         [Ar]3d104s2"
   periodic_table_element(30)%ionization = 9.390
   periodic_table_element(31)%Z = 31
   periodic_table_element(31)%mass = 69.723
   periodic_table_element(31)%name = "   Gallium    "
   periodic_table_element(31)%symbol = "Ga"
   periodic_table_element(31)%frozenp = 30
   periodic_table_element(31)%meltp = 2403
   periodic_table_element(31)%density = 5.91000
   periodic_table_element(31)%elec_conf = "      [Ar]3d104s24p1"
   periodic_table_element(31)%ionization = 6.000
   periodic_table_element(32)%Z = 32
   periodic_table_element(32)%mass = 72.640
   periodic_table_element(32)%name = "  Germanium   "
   periodic_table_element(32)%symbol = "Ge"
   periodic_table_element(32)%frozenp = 937
   periodic_table_element(32)%meltp = 2830
   periodic_table_element(32)%density = 5.32000
   periodic_table_element(32)%elec_conf = "      [Ar]3d104s24p2"
   periodic_table_element(32)%ionization = 7.900
   periodic_table_element(33)%Z = 33
   periodic_table_element(33)%mass = 74.922
   periodic_table_element(33)%name = "   Arsenic    "
   periodic_table_element(33)%symbol = "As"
   periodic_table_element(33)%frozenp = 81
   periodic_table_element(33)%meltp = 613
   periodic_table_element(33)%density = 5.72000
   periodic_table_element(33)%elec_conf = "      [Ar]3d104s24p3"
   periodic_table_element(33)%ionization = 9.790
   periodic_table_element(34)%Z = 34
   periodic_table_element(34)%mass = 78.960
   periodic_table_element(34)%name = "   Selenium   "
   periodic_table_element(34)%symbol = "Se"
   periodic_table_element(34)%frozenp = 217
   periodic_table_element(34)%meltp = 685
   periodic_table_element(34)%density = 4.79000
   periodic_table_element(34)%elec_conf = "      [Ar]3d104s24p4"
   periodic_table_element(34)%ionization = 9.750
   periodic_table_element(35)%Z = 35
   periodic_table_element(35)%mass = 79.904
   periodic_table_element(35)%name = "   Bromine    "
   periodic_table_element(35)%symbol = "Br"
   periodic_table_element(35)%frozenp = -7
   periodic_table_element(35)%meltp = 59
   periodic_table_element(35)%density = 3.12000
   periodic_table_element(35)%elec_conf = "      [Ar]3d104s24p5"
   periodic_table_element(35)%ionization = 11.810
   periodic_table_element(36)%Z = 36
   periodic_table_element(36)%mass = 83.800
   periodic_table_element(36)%name = "   Krypton    "
   periodic_table_element(36)%symbol = "Kr"
   periodic_table_element(36)%frozenp = -157
   periodic_table_element(36)%meltp = -153
   periodic_table_element(36)%density = 3.75000
   periodic_table_element(36)%elec_conf = "      [Ar]3d104s24p6"
   periodic_table_element(36)%ionization = 14.000
   periodic_table_element(37)%Z = 37
   periodic_table_element(37)%mass = 85.468
   periodic_table_element(37)%name = "   Rubidium   "
   periodic_table_element(37)%symbol = "Rb"
   periodic_table_element(37)%frozenp = 39
   periodic_table_element(37)%meltp = 688
   periodic_table_element(37)%density = 1.63000
   periodic_table_element(37)%elec_conf = "             [Kr]5s1"
   periodic_table_element(37)%ionization = 4.180
   periodic_table_element(38)%Z = 38
   periodic_table_element(38)%mass = 87.620
   periodic_table_element(38)%name = "  Strontium   "
   periodic_table_element(38)%symbol = "Sr"
   periodic_table_element(38)%frozenp = 769
   periodic_table_element(38)%meltp = 1384
   periodic_table_element(38)%density = 2.54000
   periodic_table_element(38)%elec_conf = "             [Kr]5s2"
   periodic_table_element(38)%ionization = 5.690
   periodic_table_element(39)%Z = 39
   periodic_table_element(39)%mass = 88.906
   periodic_table_element(39)%name = "   Yttrium    "
   periodic_table_element(39)%symbol = "Y "
   periodic_table_element(39)%frozenp = 1523
   periodic_table_element(39)%meltp = 3337
   periodic_table_element(39)%density = 4.47000
   periodic_table_element(39)%elec_conf = "          [Kr]4d15s2"
   periodic_table_element(39)%ionization = 6.220
   periodic_table_element(40)%Z = 40
   periodic_table_element(40)%mass = 91.224
   periodic_table_element(40)%name = "  Zirconium   "
   periodic_table_element(40)%symbol = "Zr"
   periodic_table_element(40)%frozenp = 1852
   periodic_table_element(40)%meltp = 4377
   periodic_table_element(40)%density = 6.51000
   periodic_table_element(40)%elec_conf = "          [Kr]4d25s2"
   periodic_table_element(40)%ionization = 6.630
   periodic_table_element(41)%Z = 41
   periodic_table_element(41)%mass = 92.906
   periodic_table_element(41)%name = "   Niobium    "
   periodic_table_element(41)%symbol = "Nb"
   periodic_table_element(41)%frozenp = 2468
   periodic_table_element(41)%meltp = 4927
   periodic_table_element(41)%density = 8.57000
   periodic_table_element(41)%elec_conf = "          [Kr]4d45s1"
   periodic_table_element(41)%ionization = 6.760
   periodic_table_element(42)%Z = 42
   periodic_table_element(42)%mass = 95.940
   periodic_table_element(42)%name = "  Molybdenum  "
   periodic_table_element(42)%symbol = "Mo"
   periodic_table_element(42)%frozenp = 2617
   periodic_table_element(42)%meltp = 4612
   periodic_table_element(42)%density = 10.22000
   periodic_table_element(42)%elec_conf = "          [Kr]4d55s1"
   periodic_table_element(42)%ionization = 7.090
   periodic_table_element(43)%Z = 43
   periodic_table_element(43)%mass = 98.000
   periodic_table_element(43)%name = "  Technetium  "
   periodic_table_element(43)%symbol = "Tc"
   periodic_table_element(43)%frozenp = 2200
   periodic_table_element(43)%meltp = 4877
   periodic_table_element(43)%density = 10.22000
   periodic_table_element(43)%elec_conf = "          [Kr]4d55s2"
   periodic_table_element(43)%ionization = 7.280
   periodic_table_element(44)%Z = 44
   periodic_table_element(44)%mass = 101.070
   periodic_table_element(44)%name = "  Ruthenium   "
   periodic_table_element(44)%symbol = "Ru"
   periodic_table_element(44)%frozenp = 2250
   periodic_table_element(44)%meltp = 3900
   periodic_table_element(44)%density = 12.37000
   periodic_table_element(44)%elec_conf = "          [Kr]4d75s1"
   periodic_table_element(44)%ionization = 7.360
   periodic_table_element(45)%Z = 45
   periodic_table_element(45)%mass = 102.906
   periodic_table_element(45)%name = "   Rhodium    "
   periodic_table_element(45)%symbol = "Rh"
   periodic_table_element(45)%frozenp = 1966
   periodic_table_element(45)%meltp = 3727
   periodic_table_element(45)%density = 12.41000
   periodic_table_element(45)%elec_conf = "          [Kr]4d85s1"
   periodic_table_element(45)%ionization = 7.460
   periodic_table_element(46)%Z = 46
   periodic_table_element(46)%mass = 106.420
   periodic_table_element(46)%name = "  Palladium   "
   periodic_table_element(46)%symbol = "Pd"
   periodic_table_element(46)%frozenp = 1552
   periodic_table_element(46)%meltp = 2927
   periodic_table_element(46)%density = 12.02000
   periodic_table_element(46)%elec_conf = "            [Kr]4d10"
   periodic_table_element(46)%ionization = 8.340
   periodic_table_element(47)%Z = 47
   periodic_table_element(47)%mass = 107.868
   periodic_table_element(47)%name = "    Silver    "
   periodic_table_element(47)%symbol = "Ag"
   periodic_table_element(47)%frozenp = 962
   periodic_table_element(47)%meltp = 2212
   periodic_table_element(47)%density = 10.50000
   periodic_table_element(47)%elec_conf = "         [Kr]4d105s1"
   periodic_table_element(47)%ionization = 7.580
   periodic_table_element(48)%Z = 48
   periodic_table_element(48)%mass = 112.411
   periodic_table_element(48)%name = "   Cadmium    "
   periodic_table_element(48)%symbol = "Cd"
   periodic_table_element(48)%frozenp = 321
   periodic_table_element(48)%meltp = 765
   periodic_table_element(48)%density = 8.65000
   periodic_table_element(48)%elec_conf = "         [Kr]4d105s2"
   periodic_table_element(48)%ionization = 8.990
   periodic_table_element(49)%Z = 49
   periodic_table_element(49)%mass = 114.818
   periodic_table_element(49)%name = "    Indium    "
   periodic_table_element(49)%symbol = "In"
   periodic_table_element(49)%frozenp = 157
   periodic_table_element(49)%meltp = 2000
   periodic_table_element(49)%density = 7.31000
   periodic_table_element(49)%elec_conf = "      [Kr]4d105s25p1"
   periodic_table_element(49)%ionization = 5.790
   periodic_table_element(50)%Z = 50
   periodic_table_element(50)%mass = 118.710
   periodic_table_element(50)%name = "     Tin      "
   periodic_table_element(50)%symbol = "Sn"
   periodic_table_element(50)%frozenp = 232
   periodic_table_element(50)%meltp = 2270
   periodic_table_element(50)%density = 7.31000
   periodic_table_element(50)%elec_conf = "      [Kr]4d105s25p2"
   periodic_table_element(50)%ionization = 7.340
   periodic_table_element(51)%Z = 51
   periodic_table_element(51)%mass = 121.760
   periodic_table_element(51)%name = "   Antimony   "
   periodic_table_element(51)%symbol = "Sb"
   periodic_table_element(51)%frozenp = 630
   periodic_table_element(51)%meltp = 1750
   periodic_table_element(51)%density = 6.68000
   periodic_table_element(51)%elec_conf = "      [Kr]4d105s25p3"
   periodic_table_element(51)%ionization = 8.610
   periodic_table_element(52)%Z = 52
   periodic_table_element(52)%mass = 127.600
   periodic_table_element(52)%name = "  Tellurium   "
   periodic_table_element(52)%symbol = "Te"
   periodic_table_element(52)%frozenp = 449
   periodic_table_element(52)%meltp = 990
   periodic_table_element(52)%density = 6.24000
   periodic_table_element(52)%elec_conf = "      [Kr]4d105s25p4"
   periodic_table_element(52)%ionization = 9.010
   periodic_table_element(53)%Z = 53
   periodic_table_element(53)%mass = 126.905
   periodic_table_element(53)%name = "    Iodine    "
   periodic_table_element(53)%symbol = "I "
   periodic_table_element(53)%frozenp = 114
   periodic_table_element(53)%meltp = 184
   periodic_table_element(53)%density = 4.93000
   periodic_table_element(53)%elec_conf = "      [Kr]4d105s25p5"
   periodic_table_element(53)%ionization = 10.450
   periodic_table_element(54)%Z = 54
   periodic_table_element(54)%mass = 131.293
   periodic_table_element(54)%name = "    Xenon     "
   periodic_table_element(54)%symbol = "Xe"
   periodic_table_element(54)%frozenp = -112
   periodic_table_element(54)%meltp = -108
   periodic_table_element(54)%density = 5.90000
   periodic_table_element(54)%elec_conf = "      [Kr]4d105s25p6"
   periodic_table_element(54)%ionization = 12.130
   periodic_table_element(55)%Z = 55
   periodic_table_element(55)%mass = 132.906
   periodic_table_element(55)%name = "    Cesium    "
   periodic_table_element(55)%symbol = "Cs"
   periodic_table_element(55)%frozenp = 29
   periodic_table_element(55)%meltp = 678
   periodic_table_element(55)%density = 1.87000
   periodic_table_element(55)%elec_conf = "             [Xe]6s1"
   periodic_table_element(55)%ionization = 3.890
   periodic_table_element(56)%Z = 56
   periodic_table_element(56)%mass = 137.327
   periodic_table_element(56)%name = "    Barium    "
   periodic_table_element(56)%symbol = "Ba"
   periodic_table_element(56)%frozenp = 725
   periodic_table_element(56)%meltp = 1140
   periodic_table_element(56)%density = 3.59000
   periodic_table_element(56)%elec_conf = "             [Xe]6s2"
   periodic_table_element(56)%ionization = 5.210
   periodic_table_element(57)%Z = 57
   periodic_table_element(57)%mass = 138.906
   periodic_table_element(57)%name = "  Lanthanum   "
   periodic_table_element(57)%symbol = "La"
   periodic_table_element(57)%frozenp = 920
   periodic_table_element(57)%meltp = 3469
   periodic_table_element(57)%density = 6.15000
   periodic_table_element(57)%elec_conf = "          [Xe]5d16s2"
   periodic_table_element(57)%ionization = 5.580
   periodic_table_element(58)%Z = 58
   periodic_table_element(58)%mass = 140.116
   periodic_table_element(58)%name = "    Cerium    "
   periodic_table_element(58)%symbol = "Ce"
   periodic_table_element(58)%frozenp = 795
   periodic_table_element(58)%meltp = 3257
   periodic_table_element(58)%density = 6.77000
   periodic_table_element(58)%elec_conf = "       [Xe]4f15d16s2"
   periodic_table_element(58)%ionization = 5.540
   periodic_table_element(59)%Z = 59
   periodic_table_element(59)%mass = 140.908
   periodic_table_element(59)%name = " Praseodymium "
   periodic_table_element(59)%symbol = "Pr"
   periodic_table_element(59)%frozenp = 935
   periodic_table_element(59)%meltp = 3127
   periodic_table_element(59)%density = 6.77000
   periodic_table_element(59)%elec_conf = "          [Xe]4f36s2"
   periodic_table_element(59)%ionization = 5.470
   periodic_table_element(60)%Z = 60
   periodic_table_element(60)%mass = 144.240
   periodic_table_element(60)%name = "  Neodymium   "
   periodic_table_element(60)%symbol = "Nd"
   periodic_table_element(60)%frozenp = 1010
   periodic_table_element(60)%meltp = 3127
   periodic_table_element(60)%density = 7.01000
   periodic_table_element(60)%elec_conf = "          [Xe]4f46s2"
   periodic_table_element(60)%ionization = 5.530
   periodic_table_element(61)%Z = 61
   periodic_table_element(61)%mass = 145.000
   periodic_table_element(61)%name = "  Promethium  "
   periodic_table_element(61)%symbol = "Pm"
   periodic_table_element(61)%frozenp = 1100
   periodic_table_element(61)%meltp = 3000
   periodic_table_element(61)%density = 7.30000
   periodic_table_element(61)%elec_conf = "          [Xe]4f56s2"
   periodic_table_element(61)%ionization = 5.580
   periodic_table_element(62)%Z = 62
   periodic_table_element(62)%mass = 150.360
   periodic_table_element(62)%name = "   Samarium   "
   periodic_table_element(62)%symbol = "Sm"
   periodic_table_element(62)%frozenp = 1072
   periodic_table_element(62)%meltp = 1900
   periodic_table_element(62)%density = 7.52000
   periodic_table_element(62)%elec_conf = "          [Xe]4f66s2"
   periodic_table_element(62)%ionization = 5.640
   periodic_table_element(63)%Z = 63
   periodic_table_element(63)%mass = 151.964
   periodic_table_element(63)%name = "   Europium   "
   periodic_table_element(63)%symbol = "Eu"
   periodic_table_element(63)%frozenp = 822
   periodic_table_element(63)%meltp = 1597
   periodic_table_element(63)%density = 5.24000
   periodic_table_element(63)%elec_conf = "          [Xe]4f76s2"
   periodic_table_element(63)%ionization = 5.670
   periodic_table_element(64)%Z = 64
   periodic_table_element(64)%mass = 157.250
   periodic_table_element(64)%name = "  Gadolinium  "
   periodic_table_element(64)%symbol = "Gd"
   periodic_table_element(64)%frozenp = 1311
   periodic_table_element(64)%meltp = 3233
   periodic_table_element(64)%density = 7.90000
   periodic_table_element(64)%elec_conf = "       [Xe]4f75d16s2"
   periodic_table_element(64)%ionization = 6.150
   periodic_table_element(65)%Z = 65
   periodic_table_element(65)%mass = 158.925
   periodic_table_element(65)%name = "   Terbium    "
   periodic_table_element(65)%symbol = "Tb"
   periodic_table_element(65)%frozenp = 1360
   periodic_table_element(65)%meltp = 3041
   periodic_table_element(65)%density = 8.23000
   periodic_table_element(65)%elec_conf = "          [Xe]4f96s2"
   periodic_table_element(65)%ionization = 5.860
   periodic_table_element(66)%Z = 66
   periodic_table_element(66)%mass = 162.500
   periodic_table_element(66)%name = "  Dysprosium  "
   periodic_table_element(66)%symbol = "Dy"
   periodic_table_element(66)%frozenp = 1412
   periodic_table_element(66)%meltp = 2562
   periodic_table_element(66)%density = 8.55000
   periodic_table_element(66)%elec_conf = "         [Xe]4f106s2"
   periodic_table_element(66)%ionization = 5.940
   periodic_table_element(67)%Z = 67
   periodic_table_element(67)%mass = 164.930
   periodic_table_element(67)%name = "   Holmium    "
   periodic_table_element(67)%symbol = "Ho"
   periodic_table_element(67)%frozenp = 1470
   periodic_table_element(67)%meltp = 2720
   periodic_table_element(67)%density = 8.80000
   periodic_table_element(67)%elec_conf = "         [Xe]4f116s2"
   periodic_table_element(67)%ionization = 6.020
   periodic_table_element(68)%Z = 68
   periodic_table_element(68)%mass = 167.259
   periodic_table_element(68)%name = "    Erbium    "
   periodic_table_element(68)%symbol = "Er"
   periodic_table_element(68)%frozenp = 1522
   periodic_table_element(68)%meltp = 2510
   periodic_table_element(68)%density = 9.07000
   periodic_table_element(68)%elec_conf = "         [Xe]4f126s2"
   periodic_table_element(68)%ionization = 6.110
   periodic_table_element(69)%Z = 69
   periodic_table_element(69)%mass = 168.934
   periodic_table_element(69)%name = "   Thulium    "
   periodic_table_element(69)%symbol = "Tm"
   periodic_table_element(69)%frozenp = 1545
   periodic_table_element(69)%meltp = 1727
   periodic_table_element(69)%density = 9.32000
   periodic_table_element(69)%elec_conf = "         [Xe]4f136s2"
   periodic_table_element(69)%ionization = 6.180
   periodic_table_element(70)%Z = 70
   periodic_table_element(70)%mass = 173.040
   periodic_table_element(70)%name = "  Ytterbium   "
   periodic_table_element(70)%symbol = "Yb"
   periodic_table_element(70)%frozenp = 824
   periodic_table_element(70)%meltp = 1466
   periodic_table_element(70)%density = 6.90000
   periodic_table_element(70)%elec_conf = "         [Xe]4f146s2"
   periodic_table_element(70)%ionization = 6.250
   periodic_table_element(71)%Z = 71
   periodic_table_element(71)%mass = 174.967
   periodic_table_element(71)%name = "   Lutetium   "
   periodic_table_element(71)%symbol = "Lu"
   periodic_table_element(71)%frozenp = 1656
   periodic_table_element(71)%meltp = 3315
   periodic_table_element(71)%density = 9.84000
   periodic_table_element(71)%elec_conf = "      [Xe]4f145d16s2"
   periodic_table_element(71)%ionization = 5.430
   periodic_table_element(72)%Z = 72
   periodic_table_element(72)%mass = 178.490
   periodic_table_element(72)%name = "   Hafnium    "
   periodic_table_element(72)%symbol = "Hf"
   periodic_table_element(72)%frozenp = 2150
   periodic_table_element(72)%meltp = 5400
   periodic_table_element(72)%density = 13.31000
   periodic_table_element(72)%elec_conf = "      [Xe]4f145d26s2"
   periodic_table_element(72)%ionization = 6.830
   periodic_table_element(73)%Z = 73
   periodic_table_element(73)%mass = 180.948
   periodic_table_element(73)%name = "   Tantalum   "
   periodic_table_element(73)%symbol = "Ta"
   periodic_table_element(73)%frozenp = 2996
   periodic_table_element(73)%meltp = 5425
   periodic_table_element(73)%density = 16.65000
   periodic_table_element(73)%elec_conf = "      [Xe]4f145d36s2"
   periodic_table_element(73)%ionization = 7.550
   periodic_table_element(74)%Z = 74
   periodic_table_element(74)%mass = 183.840
   periodic_table_element(74)%name = "   Tungsten   "
   periodic_table_element(74)%symbol = "W "
   periodic_table_element(74)%frozenp = 3410
   periodic_table_element(74)%meltp = 5660
   periodic_table_element(74)%density = 19.35000
   periodic_table_element(74)%elec_conf = "      [Xe]4f145d46s2"
   periodic_table_element(74)%ionization = 7.860
   periodic_table_element(75)%Z = 75
   periodic_table_element(75)%mass = 186.207
   periodic_table_element(75)%name = "   Rhenium    "
   periodic_table_element(75)%symbol = "Re"
   periodic_table_element(75)%frozenp = 3180
   periodic_table_element(75)%meltp = 5627
   periodic_table_element(75)%density = 21.04000
   periodic_table_element(75)%elec_conf = "      [Xe]4f145d56s2"
   periodic_table_element(75)%ionization = 7.830
   periodic_table_element(76)%Z = 76
   periodic_table_element(76)%mass = 190.230
   periodic_table_element(76)%name = "    Osmium    "
   periodic_table_element(76)%symbol = "Os"
   periodic_table_element(76)%frozenp = 3045
   periodic_table_element(76)%meltp = 5027
   periodic_table_element(76)%density = 22.60000
   periodic_table_element(76)%elec_conf = "      [Xe]4f145d66s2"
   periodic_table_element(76)%ionization = 8.440
   periodic_table_element(77)%Z = 77
   periodic_table_element(77)%mass = 192.217
   periodic_table_element(77)%name = "   Iridium    "
   periodic_table_element(77)%symbol = "Ir"
   periodic_table_element(77)%frozenp = 2410
   periodic_table_element(77)%meltp = 4527
   periodic_table_element(77)%density = 22.40000
   periodic_table_element(77)%elec_conf = "      [Xe]4f145d76s2"
   periodic_table_element(77)%ionization = 8.970
   periodic_table_element(78)%Z = 78
   periodic_table_element(78)%mass = 195.078
   periodic_table_element(78)%name = "   Platinum   "
   periodic_table_element(78)%symbol = "Pt"
   periodic_table_element(78)%frozenp = 1772
   periodic_table_element(78)%meltp = 3827
   periodic_table_element(78)%density = 21.45000
   periodic_table_element(78)%elec_conf = "      [Xe]4f145d96s1"
   periodic_table_element(78)%ionization = 8.960
   periodic_table_element(79)%Z = 79
   periodic_table_element(79)%mass = 196.967
   periodic_table_element(79)%name = "     Gold     "
   periodic_table_element(79)%symbol = "Au"
   periodic_table_element(79)%frozenp = 1064
   periodic_table_element(79)%meltp = 2807
   periodic_table_element(79)%density = 19.32000
   periodic_table_element(79)%elec_conf = "     [Xe]4f145d106s1"
   periodic_table_element(79)%ionization = 9.230
   periodic_table_element(80)%Z = 80
   periodic_table_element(80)%mass = 200.590
   periodic_table_element(80)%name = "   Mercury    "
   periodic_table_element(80)%symbol = "Hg"
   periodic_table_element(80)%frozenp = -39
   periodic_table_element(80)%meltp = 357
   periodic_table_element(80)%density = 13.55000
   periodic_table_element(80)%elec_conf = "     [Xe]4f145d106s2"
   periodic_table_element(80)%ionization = 10.440
   periodic_table_element(81)%Z = 81
   periodic_table_element(81)%mass = 204.383
   periodic_table_element(81)%name = "   Thallium   "
   periodic_table_element(81)%symbol = "Tl"
   periodic_table_element(81)%frozenp = 303
   periodic_table_element(81)%meltp = 1457
   periodic_table_element(81)%density = 11.85000
   periodic_table_element(81)%elec_conf = "  [Xe]4f145d106s26p1"
   periodic_table_element(81)%ionization = 6.110
   periodic_table_element(82)%Z = 82
   periodic_table_element(82)%mass = 207.200
   periodic_table_element(82)%name = "     Lead     "
   periodic_table_element(82)%symbol = "Pb"
   periodic_table_element(82)%frozenp = 327
   periodic_table_element(82)%meltp = 1740
   periodic_table_element(82)%density = 11.35000
   periodic_table_element(82)%elec_conf = "  [Xe]4f145d106s26p2"
   periodic_table_element(82)%ionization = 7.420
   periodic_table_element(83)%Z = 83
   periodic_table_element(83)%mass = 208.980
   periodic_table_element(83)%name = "   Bismuth    "
   periodic_table_element(83)%symbol = "Bi"
   periodic_table_element(83)%frozenp = 271
   periodic_table_element(83)%meltp = 1560
   periodic_table_element(83)%density = 9.75000
   periodic_table_element(83)%elec_conf = "  [Xe]4f145d106s26p3"
   periodic_table_element(83)%ionization = 7.290
   periodic_table_element(84)%Z = 84
   periodic_table_element(84)%mass = 209.000
   periodic_table_element(84)%name = "   Polonium   "
   periodic_table_element(84)%symbol = "Po"
   periodic_table_element(84)%frozenp = 254
   periodic_table_element(84)%meltp = 962
   periodic_table_element(84)%density = 9.30000
   periodic_table_element(84)%elec_conf = "  [Xe]4f145d106s26p4"
   periodic_table_element(84)%ionization = 8.420
   periodic_table_element(85)%Z = 85
   periodic_table_element(85)%mass = 210.000
   periodic_table_element(85)%name = "   Astatine   "
   periodic_table_element(85)%symbol = "At"
   periodic_table_element(85)%frozenp = 302
   periodic_table_element(85)%meltp = 337
   periodic_table_element(85)%density = 0.00000
   periodic_table_element(85)%elec_conf = "  [Xe]4f145d106s26p5"
   periodic_table_element(85)%ionization = 9.300
   periodic_table_element(86)%Z = 86
   periodic_table_element(86)%mass = 222.000
   periodic_table_element(86)%name = "    Radon     "
   periodic_table_element(86)%symbol = "Rn"
   periodic_table_element(86)%frozenp = -71
   periodic_table_element(86)%meltp = -62
   periodic_table_element(86)%density = 9.73000
   periodic_table_element(86)%elec_conf = "  [Xe]4f145d106s26p6"
   periodic_table_element(86)%ionization = 10.750
   periodic_table_element(87)%Z = 87
   periodic_table_element(87)%mass = 223.000
   periodic_table_element(87)%name = "   Francium   "
   periodic_table_element(87)%symbol = "Fr"
   periodic_table_element(87)%frozenp = 27
   periodic_table_element(87)%meltp = 677
   periodic_table_element(87)%density = 0.00000
   periodic_table_element(87)%elec_conf = "             [Rn]7s1"
   periodic_table_element(87)%ionization = 4.070
   periodic_table_element(88)%Z = 88
   periodic_table_element(88)%mass = 226.000
   periodic_table_element(88)%name = "    Radium    "
   periodic_table_element(88)%symbol = "Ra"
   periodic_table_element(88)%frozenp = 700
   periodic_table_element(88)%meltp = 1737
   periodic_table_element(88)%density = 5.50000
   periodic_table_element(88)%elec_conf = "             [Rn]7s2"
   periodic_table_element(88)%ionization = 5.280
   periodic_table_element(89)%Z = 89
   periodic_table_element(89)%mass = 227.000
   periodic_table_element(89)%name = "   Actinium   "
   periodic_table_element(89)%symbol = "Ac"
   periodic_table_element(89)%frozenp = 1050
   periodic_table_element(89)%meltp = 3200
   periodic_table_element(89)%density = 10.07000
   periodic_table_element(89)%elec_conf = "          [Rn]6d17s2"
   periodic_table_element(89)%ionization = 5.170
   periodic_table_element(90)%Z = 90
   periodic_table_element(90)%mass = 232.038
   periodic_table_element(90)%name = "   Thorium    "
   periodic_table_element(90)%symbol = "Th"
   periodic_table_element(90)%frozenp = 1750
   periodic_table_element(90)%meltp = 4790
   periodic_table_element(90)%density = 11.72000
   periodic_table_element(90)%elec_conf = "          [Rn]6d27s2"
   periodic_table_element(90)%ionization = 6.310
   periodic_table_element(91)%Z = 91
   periodic_table_element(91)%mass = 231.036
   periodic_table_element(91)%name = " Protactinium "
   periodic_table_element(91)%symbol = "Pa"
   periodic_table_element(91)%frozenp = 1568
   periodic_table_element(91)%meltp = 0
   periodic_table_element(91)%density = 15.40000
   periodic_table_element(91)%elec_conf = "       [Rn]5f26d17s2"
   periodic_table_element(91)%ionization = 5.890
   periodic_table_element(92)%Z = 92
   periodic_table_element(92)%mass = 238.029
   periodic_table_element(92)%name = "   Uranium    "
   periodic_table_element(92)%symbol = "U "
   periodic_table_element(92)%frozenp = 1132
   periodic_table_element(92)%meltp = 3818
   periodic_table_element(92)%density = 18.95000
   periodic_table_element(92)%elec_conf = "       [Rn]5f36d17s2"
   periodic_table_element(92)%ionization = 6.190
   periodic_table_element(93)%Z = 93
   periodic_table_element(93)%mass = 237.000
   periodic_table_element(93)%name = "  Neptunium   "
   periodic_table_element(93)%symbol = "Np"
   periodic_table_element(93)%frozenp = 640
   periodic_table_element(93)%meltp = 3902
   periodic_table_element(93)%density = 20.20000
   periodic_table_element(93)%elec_conf = "       [Rn]5f46d17s2"
   periodic_table_element(93)%ionization = 6.270
   periodic_table_element(94)%Z = 94
   periodic_table_element(94)%mass = 244.000
   periodic_table_element(94)%name = "  Plutonium   "
   periodic_table_element(94)%symbol = "Pu"
   periodic_table_element(94)%frozenp = 640
   periodic_table_element(94)%meltp = 3235
   periodic_table_element(94)%density = 19.84000
   periodic_table_element(94)%elec_conf = "          [Rn]5f67s2"
   periodic_table_element(94)%ionization = 6.030
   periodic_table_element(95)%Z = 95
   periodic_table_element(95)%mass = 243.000
   periodic_table_element(95)%name = "  Americium   "
   periodic_table_element(95)%symbol = "Am"
   periodic_table_element(95)%frozenp = 994
   periodic_table_element(95)%meltp = 2607
   periodic_table_element(95)%density = 13.67000
   periodic_table_element(95)%elec_conf = "          [Rn]5f77s2"
   periodic_table_element(95)%ionization = 5.970
   periodic_table_element(96)%Z = 96
   periodic_table_element(96)%mass = 247.000
   periodic_table_element(96)%name = "    Curium    "
   periodic_table_element(96)%symbol = "Cm"
   periodic_table_element(96)%frozenp = 1340
   periodic_table_element(96)%meltp = 0
   periodic_table_element(96)%density = 13.50000
   periodic_table_element(96)%elec_conf = "             [Rn]xxx"
   periodic_table_element(96)%ionization = 5.990
   periodic_table_element(97)%Z = 97
   periodic_table_element(97)%mass = 247.000
   periodic_table_element(97)%name = "  Berkelium   "
   periodic_table_element(97)%symbol = "Bk"
   periodic_table_element(97)%frozenp = 986
   periodic_table_element(97)%meltp = 0
   periodic_table_element(97)%density = 14.78000
   periodic_table_element(97)%elec_conf = "             [Rn]xxx"
   periodic_table_element(97)%ionization = 6.200
   periodic_table_element(98)%Z = 98
   periodic_table_element(98)%mass = 251.000
   periodic_table_element(98)%name = " Californium  "
   periodic_table_element(98)%symbol = "Cf"
   periodic_table_element(98)%frozenp = 900
   periodic_table_element(98)%meltp = 0
   periodic_table_element(98)%density = 15.10000
   periodic_table_element(98)%elec_conf = "             [Rn]xxx"
   periodic_table_element(98)%ionization = 6.280
   periodic_table_element(99)%Z = 99
   periodic_table_element(99)%mass = 252.000
   periodic_table_element(99)%name = " Einsteinium  "
   periodic_table_element(99)%symbol = "Es"
   periodic_table_element(99)%frozenp = 860
   periodic_table_element(99)%meltp = 0
   periodic_table_element(99)%density = 0.00000
   periodic_table_element(99)%elec_conf = "             [Rn]xxx"
   periodic_table_element(99)%ionization = 6.420
   periodic_table_element(100)%Z = 100
   periodic_table_element(100)%mass = 257.000
   periodic_table_element(100)%name = "   Fermium    "
   periodic_table_element(100)%symbol = "Fm"
   periodic_table_element(100)%frozenp = 1527
   periodic_table_element(100)%meltp = 0
   periodic_table_element(100)%density = 0.00000
   periodic_table_element(100)%elec_conf = "             [Rn]xxx"
   periodic_table_element(100)%ionization = 6.500
   periodic_table_element(101)%Z = 101
   periodic_table_element(101)%mass = 258.000
   periodic_table_element(101)%name = " Mendelevium  "
   periodic_table_element(101)%symbol = "Md"
   periodic_table_element(101)%frozenp = 0
   periodic_table_element(101)%meltp = 0
   periodic_table_element(101)%density = 0.00000
   periodic_table_element(101)%elec_conf = "             [Rn]xxx"
   periodic_table_element(101)%ionization = 6.580
   periodic_table_element(102)%Z = 102
   periodic_table_element(102)%mass = 259.000
   periodic_table_element(102)%name = "   Nobelium   "
   periodic_table_element(102)%symbol = "No"
   periodic_table_element(102)%frozenp = 827
   periodic_table_element(102)%meltp = 0
   periodic_table_element(102)%density = 0.00000
   periodic_table_element(102)%elec_conf = "             [Rn]xxx"
   periodic_table_element(102)%ionization = 6.650
   periodic_table_element(103)%Z = 103
   periodic_table_element(103)%mass = 262.000
   periodic_table_element(103)%name = "  Lawrencium  "
   periodic_table_element(103)%symbol = "Lr"
   periodic_table_element(103)%frozenp = 1627
   periodic_table_element(103)%meltp = 0
   periodic_table_element(103)%density = 0.00000
   periodic_table_element(103)%elec_conf = "             [Rn]xxx"
   periodic_table_element(103)%ionization = 4.900
   periodic_table_element(104)%Z = 104
   periodic_table_element(104)%mass = 261.000
   periodic_table_element(104)%name = " Rutherfordium "
   periodic_table_element(104)%symbol = "Rf"
   periodic_table_element(104)%frozenp = 0
   periodic_table_element(104)%meltp = 0
   periodic_table_element(104)%density = 0.00000
   periodic_table_element(104)%elec_conf = "             [Rn]xxx"
   periodic_table_element(104)%ionization = 0.000
   periodic_table_element(105)%Z = 105
   periodic_table_element(105)%mass = 262.000
   periodic_table_element(105)%name = "   Dubnium    "
   periodic_table_element(105)%symbol = "Db"
   periodic_table_element(105)%frozenp = 0
   periodic_table_element(105)%meltp = 0
   periodic_table_element(105)%density = 0.00000
   periodic_table_element(105)%elec_conf = "             [Rn]xxx"
   periodic_table_element(105)%ionization = 0.000
   periodic_table_element(106)%Z = 106
   periodic_table_element(106)%mass = 266.000
   periodic_table_element(106)%name = "  Seaborgium  "
   periodic_table_element(106)%symbol = "Sg"
   periodic_table_element(106)%frozenp = 0
   periodic_table_element(106)%meltp = 0
   periodic_table_element(106)%density = 0.00000
   periodic_table_element(106)%elec_conf = "             [Rn]xxx"
   periodic_table_element(106)%ionization = 0.000
   periodic_table_element(107)%Z = 107
   periodic_table_element(107)%mass = 264.000
   periodic_table_element(107)%name = "   Bohrium    "
   periodic_table_element(107)%symbol = "Bh"
   periodic_table_element(107)%frozenp = 0
   periodic_table_element(107)%meltp = 0
   periodic_table_element(107)%density = 0.00000
   periodic_table_element(107)%elec_conf = "             [Rn]xxx"
   periodic_table_element(107)%ionization = 0.000
   periodic_table_element(108)%Z = 108
   periodic_table_element(108)%mass = 277.000
   periodic_table_element(108)%name = "   Hassium    "
   periodic_table_element(108)%symbol = "Hs"
   periodic_table_element(108)%frozenp = 0
   periodic_table_element(108)%meltp = 0
   periodic_table_element(108)%density = 0.00000
   periodic_table_element(108)%elec_conf = "             [Rn]xxx"
   periodic_table_element(108)%ionization = 0.000
   periodic_table_element(109)%Z = 109
   periodic_table_element(109)%mass = 268.000
   periodic_table_element(109)%name = "  Meitnerium  "
   periodic_table_element(109)%symbol = "Mt"
   periodic_table_element(109)%frozenp = 0
   periodic_table_element(109)%meltp = 0
   periodic_table_element(109)%density = 0.00000
   periodic_table_element(109)%elec_conf = "             [Rn]xxx"
   periodic_table_element(109)%ionization = 0.000

   periodic_table_element(1)%covalent_radius = 79.000000000
   periodic_table_element(2)%covalent_radius = 28.000000000
   periodic_table_element(3)%covalent_radius = 155.000000000
   periodic_table_element(4)%covalent_radius = 112.000000000
   periodic_table_element(5)%covalent_radius = 98.000000000
   periodic_table_element(6)%covalent_radius = 91.000000000
   periodic_table_element(7)%covalent_radius = 92.000000000
   periodic_table_element(8)%covalent_radius = 66.000000000
   periodic_table_element(9)%covalent_radius = 57.000000000
   periodic_table_element(10)%covalent_radius = 58.000000000
   periodic_table_element(11)%covalent_radius = 190.000000000
   periodic_table_element(12)%covalent_radius = 160.000000000
   periodic_table_element(13)%covalent_radius = 143.000000000
   periodic_table_element(14)%covalent_radius = 132.000000000
   periodic_table_element(15)%covalent_radius = 128.000000000
   periodic_table_element(16)%covalent_radius = 127.000000000
   periodic_table_element(17)%covalent_radius = 102.000000000
   periodic_table_element(18)%covalent_radius = 106.000000000
   periodic_table_element(19)%covalent_radius = 235.000000000
   periodic_table_element(20)%covalent_radius = 197.000000000
   periodic_table_element(21)%covalent_radius = 162.000000000
   periodic_table_element(22)%covalent_radius = 147.000000000
   periodic_table_element(23)%covalent_radius = 134.000000000
   periodic_table_element(24)%covalent_radius = 130.000000000
   periodic_table_element(25)%covalent_radius = 135.000000000
   periodic_table_element(26)%covalent_radius = 126.000000000
   periodic_table_element(27)%covalent_radius = 125.000000000
   periodic_table_element(28)%covalent_radius = 124.000000000
   periodic_table_element(29)%covalent_radius = 128.000000000
   periodic_table_element(30)%covalent_radius = 138.000000000
   periodic_table_element(31)%covalent_radius = 141.000000000
   periodic_table_element(32)%covalent_radius = 137.000000000
   periodic_table_element(33)%covalent_radius = 139.000000000
   periodic_table_element(34)%covalent_radius = 140.000000000
   periodic_table_element(35)%covalent_radius = 120.000000000
   periodic_table_element(36)%covalent_radius = 116.000000000
   periodic_table_element(37)%covalent_radius = 248.000000000
   periodic_table_element(38)%covalent_radius = 215.000000000
   periodic_table_element(39)%covalent_radius = 178.000000000
   periodic_table_element(40)%covalent_radius = 160.000000000
   periodic_table_element(41)%covalent_radius = 146.000000000
   periodic_table_element(42)%covalent_radius = 139.000000000
   periodic_table_element(43)%covalent_radius = 136.000000000
   periodic_table_element(44)%covalent_radius = 134.000000000
   periodic_table_element(45)%covalent_radius = 134.000000000
   periodic_table_element(46)%covalent_radius = 137.000000000
   periodic_table_element(47)%covalent_radius = 144.000000000
   periodic_table_element(48)%covalent_radius = 154.000000000
   periodic_table_element(49)%covalent_radius = 166.000000000
   periodic_table_element(50)%covalent_radius = 162.000000000
   periodic_table_element(51)%covalent_radius = 159.000000000
   periodic_table_element(52)%covalent_radius = 160.000000000
   periodic_table_element(53)%covalent_radius = 139.000000000
   periodic_table_element(54)%covalent_radius = 140.000000000
   periodic_table_element(55)%covalent_radius = 244.000000000
   periodic_table_element(56)%covalent_radius = 222.000000000
   periodic_table_element(57)%covalent_radius = 187.000000000
   periodic_table_element(58)%covalent_radius = 181.000000000
   periodic_table_element(59)%covalent_radius = 182.000000000
   periodic_table_element(60)%covalent_radius = 182.000000000
   periodic_table_element(61)%covalent_radius = 199.000000000
   periodic_table_element(62)%covalent_radius = 181.000000000
   periodic_table_element(63)%covalent_radius = 199.000000000
   periodic_table_element(64)%covalent_radius = 179.000000000
   periodic_table_element(65)%covalent_radius = 180.000000000
   periodic_table_element(66)%covalent_radius = 180.000000000
   periodic_table_element(67)%covalent_radius = 179.000000000
   periodic_table_element(68)%covalent_radius = 178.000000000
   periodic_table_element(69)%covalent_radius = 177.000000000
   periodic_table_element(70)%covalent_radius = 194.000000000
   periodic_table_element(71)%covalent_radius = 175.000000000
   periodic_table_element(72)%covalent_radius = 167.000000000
   periodic_table_element(73)%covalent_radius = 149.000000000
   periodic_table_element(74)%covalent_radius = 141.000000000
   periodic_table_element(75)%covalent_radius = 137.000000000
   periodic_table_element(76)%covalent_radius = 135.000000000
   periodic_table_element(77)%covalent_radius = 136.000000000
   periodic_table_element(78)%covalent_radius = 139.000000000
   periodic_table_element(79)%covalent_radius = 146.000000000
   periodic_table_element(80)%covalent_radius = 157.000000000
   periodic_table_element(81)%covalent_radius = 171.000000000
   periodic_table_element(82)%covalent_radius = 175.000000000
   periodic_table_element(83)%covalent_radius = 170.000000000
   periodic_table_element(84)%covalent_radius = 176.000000000
   periodic_table_element(85)%covalent_radius = 150.000000000
   periodic_table_element(86)%covalent_radius = 150.000000000
   periodic_table_element(87)%covalent_radius = 260.000000000
   periodic_table_element(88)%covalent_radius = 221.000000000
   periodic_table_element(89)%covalent_radius = 188.000000000
   periodic_table_element(90)%covalent_radius = 180.000000000
   periodic_table_element(91)%covalent_radius = 161.000000000
   periodic_table_element(92)%covalent_radius = 138.000000000
   periodic_table_element(93)%covalent_radius = 130.000000000
   periodic_table_element(94)%covalent_radius = 151.000000000
   periodic_table_element(95)%covalent_radius = 173.000000000
   periodic_table_element(96)%covalent_radius = 169.000000000
   periodic_table_element(97)%covalent_radius = 168.000000000
   periodic_table_element(98)%covalent_radius = 168.000000000
   periodic_table_element(99)%covalent_radius = 165.000000000
   periodic_table_element(100)%covalent_radius = 167.000000000
   periodic_table_element(101)%covalent_radius = 173.000000000
   periodic_table_element(102)%covalent_radius = 176.000000000
   periodic_table_element(103)%covalent_radius = 161.000000000
   periodic_table_element(104)%covalent_radius = 157.000000000
   periodic_table_element(105)%covalent_radius = 149.000000000
   periodic_table_element(106)%covalent_radius = 143.000000000
   periodic_table_element(107)%covalent_radius = 141.000000000
   periodic_table_element(108)%covalent_radius = 134.000000000
   periodic_table_element(109)%covalent_radius = 129.000000000

   !$! do ii = 1, 109
   !$!   write(234,*) '"',periodic_table_element(ii)%symbol,'"'
   !$! end do

end subroutine periodic_table




subroutine repport_time(imessage, t1, t2, text)
   use ml_in_ndm_module, only: rangml

   implicit none

   integer, intent(in)  :: imessage
   real(kind(1.d0)), intent(in)     :: t1, t2
   character(len=*), intent(in)     :: text
   integer, dimension(4)      :: len_target = (/60, 50, 40, 30/)
   integer  :: leng, idiff
   character(len=80)    :: chtmp


   leng = len_trim(text)
   if (leng + 4 >= len_target(imessage)) then
      ! " .: " : this is the minimal.
      idiff = 4
   else
      idiff = len_target(imessage) - leng - 1
   end if

   if (rangml == 0) then
      !write (chtmp, *) '(a', leng, ',1x,', idiff, '(".")', ',":",f15.3)'
      write (chtmp,'(a,i0,a,i0,a)') '(a', leng, ',1x,', idiff, '(".")'//',":",f15.3)'
      write (*, FMT=chtmp) trim(text), t2 - t1
   end if

end subroutine repport_time


module mesh_grid
   use module_kind_variables, only: kind_double
   implicit none

contains

   subroutine  linear_grid(npoints, xp_min, xp_max, xp_grid)
      implicit none
      integer, intent(in) :: npoints
      real(kind_double), intent(in) :: xp_min, xp_max
      real(kind_double), dimension(:), allocatable, intent(out) :: xp_grid
      integer :: ii

      if (allocated(xp_grid)) deallocate(xp_grid)
      allocate(xp_grid(npoints))
      if (npoints==1) then
         xp_grid(1) = xp_min
      else
         do ii = 1, npoints
            xp_grid(ii) = xp_min + dble(ii -1)*(xp_max-xp_min)/dble(npoints-1)
         end do
      end if

   end subroutine linear_grid

   subroutine  log_grid(npoints, xp_min, xp_max, xp_grid)
      implicit none
      integer, intent(in) :: npoints
      real(kind_double), intent(in) :: xp_min, xp_max
      real(kind_double), dimension(:), allocatable, intent(out) :: xp_grid
      real(kind_double) :: step
      integer :: ii

      if (allocated(xp_grid)) deallocate(xp_grid)
      allocate(xp_grid(npoints))
      if (npoints==1) then
         allocate(xp_grid(1))
         xp_grid(1) = xp_min
      else
         step = (log(xp_max) - log(xp_min))/dble(npoints - 1)
         do ii = 1, npoints
            xp_grid(ii) = exp(log(xp_min) + dble(ii-1)*step)
         end do
      end if

   end subroutine log_grid

end module mesh_grid
