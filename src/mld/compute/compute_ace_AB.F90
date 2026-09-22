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

module module_ace_radial  
  use iso_fortran_env, only: dp => real64
  use mld_logger
  use mld_mpi, only: mld_mpi_abort, mld_rank
  use RadialFunctions, only: RadialFunction
  use module_spline_interpolation_mine, only: CubicSpline

  implicit none 
  private
  type(RadialFunction) :: PolyRad
  type(CubicSpline), dimension(:), allocatable :: radial_spline
  
  ! ACE chemical radial descriptor type constants
  integer, parameter :: ACE_CHEM_RADIAL_RALF = 1
  integer, parameter :: ACE_CHEM_RADIAL_BLOCK_HSVD = 2
  integer, parameter :: ACE_CHEM_RADIAL_HSVD = 3
  integer, parameter :: ACE_CHEM_RADIAL_RANDPROJ = 5
  integer, parameter :: ACE_CHEM_RADIAL_CHEMMAP_HSVD = 4
  
  type type_dico_rad
    integer :: mua, muj, nn, ll, kk 
    logical :: lzero 
  end type type_dico_rad

  type type_alpha_beta_params 
    real(dp) :: alpha, beta, radius 
  end type type_alpha_beta_params

  type ace_radial
    integer:: kmax, nmax, lmax, dim_mu, dim_rad 
    integer:: type_f_radial, type_chem_radial
    integer :: type_fcut_in, type_fcut_out 
    integer :: npoints 
    real(dp) :: r_cut_out, r_cut_width_out, r_cut_in, r_cut_width_in
    integer,  dimension(:,:,:,:), allocatable :: ic_rad 
    type(type_dico_rad), dimension(:), allocatable :: dico_radial
    !

    type(type_alpha_beta_params), dimension(:), allocatable :: alpha_beta_params

    ! Chemical map parameters for ACE_CHEM_RADIAL_CHEMMAP_HSVD
    integer :: d_e     ! embedding dimension (input to maps)
    integer :: d_p     ! number of chemical output channels (output from maps)
    real(dp) :: gamma_l ! angular coupling parameter
    real(dp), dimension(:,:), allocatable :: W0        ! isotropic chemical map (d_p x 2*d_e)
    real(dp), dimension(:,:,:), allocatable :: W_aniso ! anisotropic chemical map (d_p x 2*d_e x d_p for each principal direction)
    real(dp), dimension(:), allocatable :: embedding_vec ! chemical embedding vector (dim_mu)
    real(dp), dimension(:,:), allocatable :: V_matrix  ! V^T from SVD (d_e x d_prop)
    real(dp), dimension(:), allocatable :: u_mean      ! mean property vector (d_prop)

    contains 
    procedure :: init => init_ace_radial
    procedure :: build => build_ace_radial
  end type ace_radial  
  type(ace_radial) :: radialace
  public :: radialace, radial_spline, ACE_CHEM_RADIAL_HSVD, ACE_CHEM_RADIAL_RALF, ACE_CHEM_RADIAL_BLOCK_HSVD, &
            ACE_CHEM_RADIAL_RANDPROJ, ACE_CHEM_RADIAL_CHEMMAP_HSVD

  contains 

  subroutine init_ace_radial(this, type_chem_radial, type_f_radial, kmax, nmax, lmax, dim_mu, &
                             r_cut_in, r_cut_out, r_cut_width_in, r_cut_width_out, lambda, npoints)
    class(ace_radial), intent(inout) :: this
    integer, intent(in) :: nmax, kmax, lmax, dim_mu
    integer, intent(in) :: type_f_radial ! pure radial function, Paftounty, Bessel etc 
    integer, intent(in) :: type_chem_radial ! How to treat the radial functions for different species
    integer, intent(in) :: npoints ! number of point on spline representation. 
    real(dp), intent(in) :: r_cut_in, r_cut_out, r_cut_width_in, r_cut_width_out, lambda
    integer :: icount, mua, muj, nn, kk, ll  
    this%nmax  = nmax
    this%kmax = kmax 
    this%lmax  = lmax 
    this%dim_mu =  dim_mu 
    this%type_f_radial = type_f_radial
    this%type_chem_radial = type_chem_radial
    this%r_cut_in = r_cut_in
    this%r_cut_out = r_cut_out
    this%r_cut_width_in = r_cut_width_in
    this%r_cut_width_out = r_cut_width_out
    this%npoints = npoints

    !type_chem_radial = 1 typical ACE d_mu(muj) = delta(mu, muj) and no l dependence. 
    if (this%type_chem_radial == ACE_CHEM_RADIAL_RALF ) then 
      this%dim_rad = (nmax+1)*(lmax+1)*dim_mu**2
      !this%dim_rad = (nmax+1)*dim_mu
      if (allocated(this%ic_rad)) deallocate(this%ic_rad) ; allocate(this%ic_rad(dim_mu, dim_mu, 0:nmax, 0:lmax)) 
      this%ic_rad = -777
      if (allocated(this%dico_radial)) deallocate(this%dico_radial) ; allocate(this%dico_radial(this%dim_rad))
      icount = 0 
      do mua = 1, this%dim_mu
        do muj = 1, this%dim_mu 
          do nn = 0, this%nmax 
            do ll = 0, this%lmax
              icount = icount + 1
              this%dico_radial(icount)%mua = mua
              this%dico_radial(icount)%muj = muj
              this%dico_radial(icount)%nn = nn
              this%dico_radial(icount)%kk = -999
              this%dico_radial(icount)%ll = ll
              this%dico_radial(icount)%lzero = .false.
              this%ic_rad(mua, muj, nn, ll) = icount  
            end do 
          end do     
        end do   
      end do
      if (allocated(radial_spline)) deallocate(radial_spline) ; allocate(radial_spline(icount))
    else if (this%type_chem_radial == ACE_CHEM_RADIAL_HSVD .or. &
         this%type_chem_radial == ACE_CHEM_RADIAL_BLOCK_HSVD .or. &
         this%type_chem_radial == ACE_CHEM_RADIAL_RANDPROJ) then
      !stop 'init type_chem_radial == 2 not implemented yet!' 
      ! For HSVD, dim_mu is 1 (no mu in basis). For Block-HSVD, dim_mu is mumax.
      ! However, ic_rad needs to be allocated with (dim_mu, mumax, ...) to handle neighbor species.
      ! But init_ace_radial only receives dim_mu.
      ! If HSVD, dim_mu=1. We need mumax.
      ! Assuming dim_mu passed here is the basis dimension.
      ! We need to know mumax. But it's not passed.
      ! Wait, for HSVD, dim_mu=1. But we loop muj up to dim_mu?
      ! If dim_mu=1, we only loop muj=1. This is wrong if we have multiple species.
      ! So for HSVD, we must assume dim_mu passed to init is mumax?
      ! But in base_cnlm_init we set dim_mu=1.
      ! This is a conflict.
      
      ! Let's assume for now that for HSVD, we allocate ic_rad as (1, dim_mu_real, ...)
      ! But we don't know dim_mu_real (mumax).
      
      ! Actually, if dim_mu=1 for HSVD, then ic_rad is (1, 1, ...).
      ! This implies we only support 1 species for neighbors? No.
      
      ! Revert to previous logic: dim_mu in init_ace_radial is mumax.
      ! But in base_cnlm_init, we set dim_mu=1 for HSVD.
      ! This means base_cnlm uses dim_mu=1, but radialace uses dim_mu=mumax?
      ! They are different objects.
      
      ! If radialace uses dim_mu=mumax, then ic_rad is (mumax, mumax).
      ! But compute_uniqueA uses mu=1 (from base_cnlm).
      ! So it accesses ic_rad(1, type_ia_db, ...).
      ! This works if ic_rad is allocated (mumax, mumax) and we just use the first row.
      
      this%dim_rad = this%kmax*(this%lmax+1)*this%dim_mu**2
      if (kmax == 0) then
         write(*,*) 'ACE_CHEM_RADIAL_HSVD/BLOCK_HSVD: kmax = 0 not allowed'
         stop 'critical error in init_ace_radial' 
      end if 
      !TODOace_radbody it is a symmetric vetor in mua and muj ? 
      if (allocated(this%ic_rad)) deallocate(this%ic_rad) ; allocate(this%ic_rad(dim_mu, dim_mu, 1:kmax, 0:lmax)) 
      this%ic_rad = -777
      if (allocated(this%dico_radial)) deallocate(this%dico_radial) ; allocate(this%dico_radial(this%dim_rad))
      
      icount = 0 
      do mua = 1, this%dim_mu
        do muj = 1, this%dim_mu 
          do kk = 1, this%kmax
            do ll = 0, this%lmax
              icount = icount + 1
              this%dico_radial(icount)%mua = mua
              this%dico_radial(icount)%muj = muj
              this%dico_radial(icount)%kk = kk
              this%dico_radial(icount)%nn = -999
              this%dico_radial(icount)%ll = ll
              this%dico_radial(icount)%lzero = .false.
              this%ic_rad(mua, muj, kk, ll) = icount  
            end do 
          end do     
        end do   
      end do
      if (allocated(radial_spline)) deallocate(radial_spline) ; allocate(radial_spline(icount))
    else if (this%type_chem_radial == ACE_CHEM_RADIAL_CHEMMAP_HSVD) then
      ! Chemical map version: Similar to HSVD but uses SVD-based chemical embedding
      this%dim_rad = this%kmax*(this%lmax+1)*this%dim_mu**2
      if (kmax == 0) then
         write(*,*) 'ACE_CHEM_RADIAL_CHEMMAP_HSVD: kmax = 0 not allowed'
         stop 'critical error in init_ace_radial' 
      end if 
      
      if (allocated(this%ic_rad)) deallocate(this%ic_rad) ; allocate(this%ic_rad(dim_mu, dim_mu, 1:kmax, 0:lmax)) 
      this%ic_rad = -777
      if (allocated(this%dico_radial)) deallocate(this%dico_radial) ; allocate(this%dico_radial(this%dim_rad))
      
      icount = 0 
      do mua = 1, this%dim_mu
        do muj = 1, this%dim_mu 
          do kk = 1, this%kmax
            do ll = 0, this%lmax
              icount = icount + 1
              this%dico_radial(icount)%mua = mua
              this%dico_radial(icount)%muj = muj
              this%dico_radial(icount)%kk = kk
              this%dico_radial(icount)%nn = -999
              this%dico_radial(icount)%ll = ll
              this%dico_radial(icount)%lzero = .false.
              this%ic_rad(mua, muj, kk, ll) = icount  
            end do 
          end do     
        end do   
      end do
      if (allocated(radial_spline)) deallocate(radial_spline) ; allocate(radial_spline(icount))
      
      ! Initialize chemical embedding parameters (defaults, will be set from input)
      this%d_e = 3       ! embedding dimension
      this%d_p = 6       ! number of principal chemical directions  
      this%gamma_l = 0.5_dp  ! angular coupling parameter
    end if 
    !TODOace radbody.  
    call PolyRad%init(this%type_f_radial, r_cut_in, r_cut_out, r_cut_width_in, r_cut_width_out, lambda, this%nmax)    
  end subroutine init_ace_radial

  subroutine build_ace_radial(this)
      !local variables
    !use module_chemical_species, only: fix_type_to_periodic, periodic_table_element
    use module_svd_small_gen_matrix, only: eigen_svd
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_IN_PLACE
    use mld_mpi, only: mld_rank, mld_size, comm_mld
    use time_check_general, only: MY_MPI_WTIME
    use module_ace_desc, only: ace_chem_low_rank, ace_chem_low_rank_q, &
                               ace_chem_low_rank_niter, ace_chem_low_rank_lambda, &
                               ace_svd_randomized, ace_svd_randomized_oversample, ace_svd_randomized_power_iter
    use module_chem_compress, only: chem_compressor
    implicit none 
    class(ace_radial), intent(inout) :: this
    integer :: kk, ii, ip , nn, ll, mua, muj, mm, dim_vv, ikount, itest 
    real(dp), dimension(:), allocatable :: xgrid, yfunc, d_yfunc, ymat_raw, d_ymat_raw
    real(dp) :: ffnn, dffnn, rr, dtmp1, dtmp2
    real(dp), dimension(:), allocatable :: vv, sigmak_local, snorm 
    real(dp), dimension(:,:), allocatable :: mat_vv, delta_local, d_mat_vv, vk_local, mat_pot
    integer,  dimension(:,:), allocatable :: mu_pair
    type(eigen_svd) :: svd_chem
    integer :: no_svd_to_perform
    integer :: rand_rank_local ! rank estimate returned by randomized_svd_topk
    ! MPI parallelisation of the HSVD/BLOCK_HSVD triple loop over (mua, muj, ll)
    integer :: itask, ntask, n_splines, ispline_global
    real(dp), dimension(:), allocatable :: spline_buf_local
    ! On-the-fly compression: temporary grid arrays for one (kk, ll) channel
    real(dp), dimension(:,:,:,:), allocatable :: g_grid, dg_grid   ! (npoints, S, S, kbatch)
    real(dp), dimension(:), allocatable :: xg_safe               ! clamped grid
    integer :: ich_task
    ! Ultra-low-memory compression: batched-kk approach
    type(CubicSpline), dimension(:), allocatable :: tmp_splines  ! (1)
    integer :: owner_rank, n_per_spline
    integer :: kbatch, kk_start, kk_end, kk_batch_size, kb
    integer(8) :: mem_per_kk, mem_cap_bytes
    ! MPI_Reduce to channel owner (replaces MPI_Allreduce)
    integer :: ierr_reduce, nelem_reduce
    ! Timing instrumentation for compression phases
    real(dp) :: t_svd_start, t_svd_total, t_allred_total, t_als_total, t_phase_start


    _NAMECURRENT_("compute_ace_radial")
    _MLD_BEGIN_

    !DEBUG_RADIAL ! rank 0 writes diagnostic dumps used to cross-check against python radial builders
    !DEBUG_RADIAL if (mld_rank == 0) then
    !DEBUG_RADIAL   open(unit=77, file='ace_radial_metadata.csv', status='replace', action='write')
    !DEBUG_RADIAL   open(unit=78, file='ace_radial_profiles.csv', status='replace', action='write')
    !DEBUG_RADIAL end if

    !TODOace is different species different answer. Or more complicated radial functions.
    if (this%type_chem_radial == ACE_CHEM_RADIAL_RALF) then 
    do ii  = 1, size(this%dico_radial)
      nn = this%dico_radial(ii)%nn
      ll = this%dico_radial(ii)%ll
      mua = this%dico_radial(ii)%mua
      muj = this%dico_radial(ii)%muj
      !TODOace radbody 
      call radial_spline(ii)%init(this%r_cut_in, this%r_cut_out, this%npoints, xgrid, yfunc, d_yfunc)
      do ip = 1, this%npoints
        rr = xgrid(ip)
        !TODOace radbody  
        call PolyRad%evaluate(rr)
        !TODOace: nmax 0 or 1. 
        ffnn = PolyRad%radial(nn)
        dffnn = PolyRad%d_radial(nn )
        yfunc(ip) = ffnn
        d_yfunc(ip) = dffnn
      end do 
      ! write(*,*) 'max min yfunc d_yfunc', maxval(yfunc), minval(yfunc), maxval(d_yfunc), minval(d_yfunc)

      !DEBUG_RADIAL if (mld_rank == 0) then
      !DEBUG_RADIAL   open(unit=77, file='ace_radial_profiles.csv', status='unknown', position='append')
      !DEBUG_RADIAL   do ip = 1, this%npoints
      !DEBUG_RADIAL     write(77, *) xgrid(ip), yfunc(ip), d_yfunc(ip)
      !DEBUG_RADIAL   end do
      !DEBUG_RADIAL   close(77)
      !DEBUG_RADIAL   open(unit=78, file='ace_radial_metadata.csv', status='unknown', position='append')
      !DEBUG_RADIAL   write(78, *) ii, nn, ll, mua, muj, 0.0_dp
      !DEBUG_RADIAL   close(78)
      !DEBUG_RADIAL end if

      call radial_spline(ii)%compute(yfunc, d_yfunc)

      !DEBUG_RADIAL if (mld_rank == 0) then 
      !DEBUG_RADIAL   ! For RALF: write nn instead of kk. The 4th column is interpreted as index (kk or nn).
      !DEBUG_RADIAL   write(78,'(6I7)') this%npoints, ii, ll, nn, mua, muj
      !DEBUG_RADIAL   do ip = 1, this%npoints
      !DEBUG_RADIAL        ! For RALF: write d_yfunc in the last column instead of potential
      !DEBUG_RADIAL        write(78,'(I0,",",ES23.15,",",ES23.15,",",ES23.15)') ip, xgrid(ip), yfunc(ip), d_yfunc(ip)
      !DEBUG_RADIAL   end do
      !DEBUG_RADIAL end if 
    end do

    else if ((this%type_chem_radial == ACE_CHEM_RADIAL_HSVD).or. &
         (this%type_chem_radial == ACE_CHEM_RADIAL_BLOCK_HSVD).or. &
         (this%type_chem_radial == ACE_CHEM_RADIAL_RANDPROJ)) then

      !alpha beta params ... 
      !TODOace_radbody include lambda  
      call init_alpha_beta_params_lin(this)
      !call init_alpha_beta_params_log(this)
      dim_vv =  this%dim_mu * this%nmax * size(this%alpha_beta_params) !Here version mu=mu_j
      ! dim_vv =  this%dim_mu * (this%dim_mu + 1) / 2 * this%nmax * size(this%alpha_beta_params)
      
      if (allocated(vv)) deallocate(vv) ; allocate(vv(dim_vv))
      if (allocated(mat_vv)) deallocate(mat_vv) ; allocate(mat_vv(this%npoints, dim_vv))
      if (allocated(d_mat_vv)) deallocate(d_mat_vv) ; allocate(d_mat_vv(this%npoints, dim_vv))
      if (allocated(vk_local)) deallocate(vk_local) ; allocate(vk_local(this%kmax, dim_vv))
      if (allocated(sigmak_local)) deallocate(sigmak_local) ; allocate(sigmak_local(this%kmax))
      
      !HSVDdebug 
      if (allocated(mat_pot)) deallocate(mat_pot) ; allocate(mat_pot(this%npoints, size(this%dico_radial)))
      if (allocated(snorm)) deallocate(snorm) ; allocate(snorm(dim_vv))

      if (allocated(delta_local)) deallocate(delta_local) ; allocate(delta_local(this%dim_mu, this%dim_mu)) !Here version mu=mu_j
      ! if (allocated(delta_local)) deallocate(delta_local) ; allocate(delta_local(this%dim_mu * (this%dim_mu + 1) / 2, this%dim_mu * (this%dim_mu + 1) / 2))
      delta_local(:,:) = 0.0_dp
      do mm = 1, this%dim_mu !Here version mu=mu_j
         delta_local(mm, mm) = 1.0_dp
        !  delta_local(1:mm, mm) = 1.0_dp
      end do   

      ! if (allocated(mu_pair)) deallocate(mu_pair) ; 
      ! allocate(mu_pair(this%dim_mu * (this%dim_mu + 1) / 2, 2))
      ! ikount = 0
      ! do mua = 1, this%dim_mu
      !   do muj = 1, mua
      !     ikount = ikount + 1
      !     mu_pair(ikount, 1) = mua
      !     mu_pair(ikount, 2) = muj
      !   end do 
      ! end do
      ! do mm = 1, this%dim_mu * (this%dim_mu + 1) / 2
      !   delta_local(mm, mm) = 1.0_dp
      !   do nn = 1, mm-1
      !     if (mu_pair(nn, 1) == mu_pair(mm, 1) .or. mu_pair(nn, 1) == mu_pair(mm, 2) .or. &
      !         mu_pair(nn, 2) == mu_pair(mm, 1) .or. mu_pair(nn, 2) == mu_pair(mm, 2)) then
      !       delta_local(nn, mm) = 0.5_dp
      !     end if
      !   end do
      ! end do      
         
      ntask = this%dim_mu * this%dim_mu * (this%lmax + 1) ! total number of independent SVD tasks
      no_svd_to_perform = ntask
      call log_info("ML:  there will be "//vtoa(size(this%dico_radial))//" radial functions for the ACE_CHEM_RADIAL_HSVD/BLOCK_HSVD")
      call log_info("ML:  the acek tensor reduction is "//vtoa(this%kmax))
      call log_info("ML:  the matrix from which the basis is extracted has the dimension "//vtoa(size(mat_vv, 1))//" x "//vtoa(size(mat_vv, 2)))
      call log_info("ML:  SVD tasks (mua,muj,ll) = "//vtoa(ntask)//"  distributed over "//vtoa(mld_size)//" MPI ranks")

      ! ----------------------------------------------------------------
      ! MPI parallelisation: each rank processes tasks itask such that
      !   mod(itask, mld_size) == mld_rank
      !
      ! When ace_chem_low_rank==1 the code processes one ll value at a
      ! time to avoid allocating all S²×kmax×(lmax+1) splines at once.
      ! Only S²×kmax temporary splines are alive per ll iteration.
      !
      ! When ace_chem_low_rank/=1 (standard path), all splines are
      ! allocated, computed, and allreduced in a single pass (original
      ! monolithic approach).
      ! ----------------------------------------------------------------

      if (ace_chem_low_rank == 1 .and. this%type_chem_radial /= ACE_CHEM_RADIAL_RANDPROJ) then
        ! ============================================================
        ! ULTRA-LOW-MEMORY PATH with BATCHED kk:
        !
        ! Outer loop: ll → kk-batches → (mua,muj)
        ! For each owned (mua,muj): SVD once per batch, evaluate all
        ! kk in batch → fill g_grid(:,:,:,kb). One allreduce per batch.
        !
        ! kbatch = min(kmax, floor(1GB / (2×Nr×S²×8)))
        !   → SVDs reduced by factor kbatch vs single-kk approach
        !   → For S=4:  kbatch=kmax (all kk at once), ~1 SVD per (mua,muj,ll)
        !   → For S=108: kbatch≈5, ~4 SVDs per (mua,muj,ll) instead of 20
        !
        ! Peak memory: 2 × kbatch × Nr × S² doubles (capped at ~1 GB)
        ! ============================================================
        call log_info("ML: ACE on-the-fly chemical compression (ultra-low-memory), Q="//vtoa(ace_chem_low_rank_q))
        call chem_compressor%init(this%dim_mu, ace_chem_low_rank_q, this%npoints, &
                                  this%kmax, this%lmax+1, ace_chem_low_rank_niter, &
                                  real(ace_chem_low_rank_lambda, dp), &
                                  this%r_cut_in, this%r_cut_out)

        ! Allocate one temporary spline
        allocate(tmp_splines(1))

        ! One dummy init to populate xgrid / yfunc / d_yfunc arrays
        call tmp_splines(1)%init(this%r_cut_in, this%r_cut_out, this%npoints, xgrid, yfunc, d_yfunc)

        ! Build clamped grid
        allocate(xg_safe(this%npoints))
        xg_safe(:) = xgrid(1:this%npoints)
        xg_safe(1) = max(xg_safe(1), this%r_cut_in)
        xg_safe(this%npoints) = min(xg_safe(this%npoints), this%r_cut_out)

        ! Compute kbatch: how many kk channels fit in ~1 GB
        ! Each kk needs 2 arrays of Nr×S² doubles (g_grid + dg_grid)
        mem_cap_bytes = 1073741824_8   ! 1 GB
        mem_per_kk = 2_8 * int(this%npoints, 8) * int(this%dim_mu, 8) * int(this%dim_mu, 8) * 8_8
        kbatch = max(1, int(mem_cap_bytes / mem_per_kk))
        kbatch = min(kbatch, this%kmax)
        call log_info("ML: ACE compression kbatch="//vtoa(kbatch)//" (kmax="//vtoa(this%kmax)// &
                      ", mem/kk="//vtoa(int(mem_per_kk/1048576_8))//" MB)")

        ! Allocate batched grid arrays: (npoints, S, S, kbatch)
        allocate(g_grid(this%npoints, this%dim_mu, this%dim_mu, kbatch))
        allocate(dg_grid(this%npoints, this%dim_mu, this%dim_mu, kbatch))

        t_svd_total = 0.0_dp
        t_allred_total = 0.0_dp
        t_als_total = 0.0_dp

        do ll = 0, this%lmax

          ! Process kk in batches
          do kk_start = 1, this%kmax, kbatch
            kk_end = min(kk_start + kbatch - 1, this%kmax)
            kk_batch_size = kk_end - kk_start + 1

            g_grid(:,:,:,1:kk_batch_size) = 0.0_dp
            dg_grid(:,:,:,1:kk_batch_size) = 0.0_dp

            ! Each rank fills its owned (mua,muj) entries for all kk in batch
            t_svd_start = MY_MPI_WTIME()
            itask = -1
            do mua = 1, this%dim_mu
            do muj = 1, this%dim_mu
              itask = itask + 1
              if (mod(itask, mld_size) /= mld_rank) cycle

              ! SVD: once per (mua,muj,ll) per batch
              call tmp_splines(1)%init(this%r_cut_in, this%r_cut_out, this%npoints, xgrid, yfunc, d_yfunc)
              call init_mat_params(this, this%ic_rad(mua, muj, 1, ll), xgrid, delta_local, mat_vv, d_mat_vv)
              do mm = 1, dim_vv
                snorm(mm) = sqrt( sum(mat_vv(:,mm)**2 ) )
                if (snorm(mm) <= 1.0d-14) snorm(mm) = 1.0d0
                mat_vv(:,mm)   = mat_vv(:,mm)   / snorm(mm)
                d_mat_vv(:,mm) = d_mat_vv(:,mm) / snorm(mm)
              end do
              if (ace_svd_randomized == 1) then
                call randomized_svd_topk(mat_vv, this%kmax, ace_svd_randomized_oversample, &
                                          ace_svd_randomized_power_iter, mua, muj, ll, &
                                          sigmak_local, vk_local, rand_rank_local)
                if (kk_start == 1 .and. this%kmax > rand_rank_local) then
                  call log_warning(NAMECURRENT//'ML: ACE_CHEM_RADIAL_HSVD (randomized): kmax > rank')
                  call log_warning('kmax='//vtoa(this%kmax)//' rank='//vtoa(rand_rank_local)// &
                                   ' ll='//vtoa(ll)//' mua='//vtoa(mua)//' muj='//vtoa(muj))
                end if
              else
                call svd_chem%init(size(mat_vv, 1), size(mat_vv, 2))
                call svd_chem%evaluate(mat_vv)
                if (kk_start == 1 .and. this%kmax > svd_chem%rank) then
                  call log_warning(NAMECURRENT//'ML: ACE_CHEM_RADIAL_HSVD: kmax > rank')
                  call log_warning('kmax='//vtoa(this%kmax)//' rank='//vtoa(svd_chem%rank)// &
                                   ' ll='//vtoa(ll)//' mua='//vtoa(mua)//' muj='//vtoa(muj))
                end if
                vk_local(1:this%kmax, 1:dim_vv) = svd_chem%eigvecVT(1:this%kmax, 1:dim_vv)
                sigmak_local(1:this%kmax) = svd_chem%eigval(1:this%kmax)
              end if

              ! Evaluate all kk in this batch
              do kk = kk_start, kk_end
                kb = kk - kk_start + 1  ! batch-local index

                call tmp_splines(1)%init(this%r_cut_in, this%r_cut_out, this%npoints, xgrid, yfunc, d_yfunc)
                do ip = 1, size(mat_vv, 1)
                  dtmp1 = 0.0_dp
                  dtmp2 = 0.0_dp
                  do mm = 1, dim_vv
                    dtmp1 = dtmp1 + vk_local(kk, mm) * mat_vv(ip, mm)
                    dtmp2 = dtmp2 + vk_local(kk, mm) * d_mat_vv(ip, mm)
                  end do
                  yfunc(ip)   = dtmp1 / sigmak_local(kk)
                  d_yfunc(ip) = dtmp2 / sigmak_local(kk)
                end do
                call tmp_splines(1)%compute(yfunc, d_yfunc)

                do ip = 1, this%npoints
                  g_grid(ip, mua, muj, kb)  = tmp_splines(1)%evaluate(xg_safe(ip))
                  dg_grid(ip, mua, muj, kb) = tmp_splines(1)%derivative(xg_safe(ip))
                end do
              end do  ! kk in batch

            end do  ! muj
            end do  ! mua
            t_phase_start = MY_MPI_WTIME()
            t_svd_total = t_svd_total + (t_phase_start - t_svd_start)

            ! Reduce each kk plane to its compression owner (not allreduce)
            nelem_reduce = this%npoints * this%dim_mu * this%dim_mu
            do kk = kk_start, kk_end
              kb = kk - kk_start + 1
              ich_task = (kk-1)*(this%lmax+1) + ll
              owner_rank = mod(ich_task, mld_size)
              call mld_reduce_to_owner(g_grid(1,1,1,kb), nelem_reduce, owner_rank, &
                                       comm_mld%comm, mld_rank)
              call mld_reduce_to_owner(dg_grid(1,1,1,kb), nelem_reduce, owner_rank, &
                                       comm_mld%comm, mld_rank)
            end do
            t_svd_start = MY_MPI_WTIME()
            t_allred_total = t_allred_total + (t_svd_start - t_phase_start)

            ! Compress each channel in the batch (only owner has valid g_grid)
            do kk = kk_start, kk_end
              kb = kk - kk_start + 1
              ich_task = (kk-1)*(this%lmax+1) + ll
              if (mod(ich_task, mld_size) == mld_rank) then
                call log_info("CHEM_COMPRESS: channel kk="//vtoa(kk)//" ll="//vtoa(ll)// &
                              " on rank "//vtoa(mld_rank))
                call chem_compressor%compress_channel(kk, ll, g_grid(:,:,:,kb), dg_grid(:,:,:,kb), xg_safe, &
                                                      this%r_cut_in, this%r_cut_out)
              end if
            end do  ! kk compress
            t_phase_start = MY_MPI_WTIME()
            t_als_total = t_als_total + (t_phase_start - t_svd_start)

          end do  ! kk_start batch
          call log_info("ML: ACE ll="//vtoa(ll)//" compressed (ultra-low-memory)")
        end do  ! ll

        ! MPI allreduce of all compressed data
        call chem_compressor%allreduce_compressed()

        deallocate(g_grid, dg_grid, xg_safe, tmp_splines)
        call log_info("ML: ACE on-the-fly compression complete")
        call log_info("ML: ACE compression timing: SVD="//vtoa(t_svd_total)//" s, Allreduce="// &
                      vtoa(t_allred_total)//" s, ALS="//vtoa(t_als_total)//" s")

      else
        ! ============================================================
        ! STANDARD PATH (ace_chem_low_rank /= 1):
        ! all S²×kmax×(lmax+1) splines allocated and allreduced at once
        ! ============================================================

        ! --- Initialise ALL splines on ALL ranks ---
        do mua = 1, this%dim_mu
        do muj = 1, this%dim_mu
        do ll = 0, this%lmax
          do kk = 1, this%kmax
            ii = this%ic_rad(mua, muj, kk, ll)
            call radial_spline(ii)%init(this%r_cut_in, this%r_cut_out, this%npoints, xgrid, yfunc, d_yfunc)
            radial_spline(ii)%a(:) = 0.0_dp
            radial_spline(ii)%b(:) = 0.0_dp
            radial_spline(ii)%c(:) = 0.0_dp
            radial_spline(ii)%d(:) = 0.0_dp
          end do
        end do
        end do
        end do

        ! --- Distributed SVD loop ---
        itask = -1
        do mua = 1, this%dim_mu
        do muj = 1, this%dim_mu
        do ll = 0, this%lmax
          itask = itask + 1
          if (mod(itask, mld_size) /= mld_rank) cycle

          ikount = 0
          do kk = 1, this%kmax
            ikount = ikount + 1
            ii = this%ic_rad(mua, muj, kk, ll)
            if (ikount == 1) then
              call radial_spline(ii)%init(this%r_cut_in, this%r_cut_out, this%npoints, xgrid, yfunc, d_yfunc)
              call init_mat_params(this, ii, xgrid, delta_local, mat_vv, d_mat_vv)
              do mm = 1, dim_vv
                snorm(mm) = sqrt( sum(mat_vv(:,mm)**2 ) )
                if (snorm(mm) <= 1.0d-14) snorm(mm) = 1.0d0
                mat_vv(:,mm)   = mat_vv(:,mm)   / snorm(mm)
                d_mat_vv(:,mm) = d_mat_vv(:,mm) / snorm(mm)
              end do
              if (this%type_chem_radial == ACE_CHEM_RADIAL_RANDPROJ) then
                ! Random-projection contraction (kACE Sec. III.F, eqs 36+38):
                ! replace the deterministic SVD contraction vectors t^(k) by
                ! generic Gaussian directions W^k ~ N(0,1/kmax) over the joint
                ! radial-chemical index; NO singular-value rescaling (sigma=1).
                call fill_random_projection_weights(this%kmax, dim_vv, mua, muj, ll, vk_local)
                sigmak_local(1:this%kmax) = 1.0_dp
              else if (ace_svd_randomized == 1) then
                call randomized_svd_topk(mat_vv, this%kmax, ace_svd_randomized_oversample, &
                                          ace_svd_randomized_power_iter, mua, muj, ll, &
                                          sigmak_local, vk_local, rand_rank_local)
                if (this%kmax > rand_rank_local) then
                  call log_warning(NAMECURRENT//'ML: ACE_CHEM_RADIAL_HSVD/BLOCK_HSVD (randomized): kmax > rank of the matrix')
                  call log_warning('kmax = '//vtoa(this%kmax)//' rank = '//vtoa(rand_rank_local)//' ll = '//vtoa(ll)// &
                                   ' kk = '//vtoa(kk)//' mua = '//vtoa(mua)//' muj = '//vtoa(muj)//' ')
                end if
                ! mat_pot needs the left singular vectors, which randomized_svd_topk
                ! never forms (they are unused downstream); left unset in this branch.
              else
                call svd_chem%init(size(mat_vv, 1), size(mat_vv, 2))
                call svd_chem%evaluate(mat_vv)
                if (this%kmax > svd_chem%rank) then
                  call log_warning(NAMECURRENT//'ML: ACE_CHEM_RADIAL_HSVD/BLOCK_HSVD: kmax > rank of the matrix')
                  call log_warning('kmax = '//vtoa(this%kmax)//' rank = '//vtoa(svd_chem%rank)//' ll = '//vtoa(ll)// &
                                   ' kk = '//vtoa(kk)//' mua = '//vtoa(mua)//' muj = '//vtoa(muj)//' ')
                end if
                vk_local(1:this%kmax, 1:dim_vv) = svd_chem%eigvecVT(1:this%kmax, 1:dim_vv)
                sigmak_local(1:this%kmax) = svd_chem%eigval(1:this%kmax)
                mat_pot(:,ii) = svd_chem%eigvecU(:, kk) * sigmak_local(kk) / (sigmak_local(kk))
              end if
            else
              call radial_spline(ii)%init(this%r_cut_in, this%r_cut_out, this%npoints, xgrid, yfunc, d_yfunc)
              if (this%type_chem_radial /= ACE_CHEM_RADIAL_RANDPROJ .and. ace_svd_randomized /= 1) then
                mat_pot(:,ii) = svd_chem%eigvecU(:, kk) * sigmak_local(kk) / (sigmak_local(kk))
              end if
            end if
            do ip = 1, size(mat_vv, 1)
              dtmp1 = 0.0_dp
              dtmp2 = 0.0_dp
              do mm = 1, dim_vv
                dtmp1 = dtmp1 + vk_local(kk, mm) * mat_vv(ip, mm)
                dtmp2 = dtmp2 + vk_local(kk, mm) * d_mat_vv(ip, mm)
              end do
              yfunc(ip)   = dtmp1 / sigmak_local(kk)
              d_yfunc(ip) = dtmp2 / sigmak_local(kk)
            end do
            call radial_spline(ii)%compute(yfunc, d_yfunc)
          end do  ! kk
        end do  ! ll
        end do  ! muj
        end do  ! mua

        ! --- Pack, allreduce, unpack all splines ---
        n_splines = size(radial_spline)
        ispline_global = n_splines * (4 * this%npoints - 1)
        if (allocated(spline_buf_local)) deallocate(spline_buf_local)
        allocate(spline_buf_local(ispline_global), source=0.0_dp)

        do ii = 1, n_splines
          mm = (ii-1) * (4 * this%npoints - 1)
          spline_buf_local(mm+1                : mm+  this%npoints  ) = radial_spline(ii)%a(1:this%npoints)
          spline_buf_local(mm+  this%npoints+1 : mm+2*this%npoints  ) = radial_spline(ii)%b(1:this%npoints)
          spline_buf_local(mm+2*this%npoints+1 : mm+3*this%npoints  ) = radial_spline(ii)%c(1:this%npoints)
          spline_buf_local(mm+3*this%npoints+1 : mm+4*this%npoints-1) = radial_spline(ii)%d(1:this%npoints-1)
        end do

        call comm_mld%sum(spline_buf_local)

        do ii = 1, n_splines
          mm = (ii-1) * (4 * this%npoints - 1)
          radial_spline(ii)%a(1:this%npoints)   = spline_buf_local(mm+1                : mm+  this%npoints  )
          radial_spline(ii)%b(1:this%npoints)   = spline_buf_local(mm+  this%npoints+1 : mm+2*this%npoints  )
          radial_spline(ii)%c(1:this%npoints)   = spline_buf_local(mm+2*this%npoints+1 : mm+3*this%npoints  )
          radial_spline(ii)%d(1:this%npoints-1) = spline_buf_local(mm+3*this%npoints+1 : mm+4*this%npoints-1)
        end do

        deallocate(spline_buf_local)

      end if  ! ace_chem_low_rank == 1
    
    else if (this%type_chem_radial == ACE_CHEM_RADIAL_CHEMMAP_HSVD) then
      ! Chemical map version: Use SVD-based linear chemical embedding
      ! No alpha_beta_params needed - chemical embedding replaces V_nl potentials
      call init_chemical_embedding(this)
      
      ! dim_vv = d_p * nmax (chemical channels × radial channels)
      dim_vv = this%d_p * this%nmax
      
      if (allocated(vv)) deallocate(vv) ; allocate(vv(dim_vv))
      if (allocated(mat_vv)) deallocate(mat_vv) ; allocate(mat_vv(this%npoints, dim_vv))
      if (allocated(d_mat_vv)) deallocate(d_mat_vv) ; allocate(d_mat_vv(this%npoints, dim_vv))
      if (allocated(vk_local)) deallocate(vk_local) ; allocate(vk_local(this%kmax, dim_vv))
      if (allocated(sigmak_local)) deallocate(sigmak_local) ; allocate(sigmak_local(this%kmax))
      if (allocated(snorm)) deallocate(snorm) ; allocate(snorm(dim_vv))
      
      call log_info("ML:  ACE_CHEM_RADIAL_CHEMMAP_HSVD with "//vtoa(size(this%dico_radial))//" radial functions")
      call log_info("ML:  kmax tensor reduction = "//vtoa(this%kmax))
      call log_info("ML:  chemical embedding dimension d_e = "//vtoa(this%d_e))
      call log_info("ML:  chemical output channels d_p = "//vtoa(this%d_p))
      call log_info("ML:  basis matrix dimension (npoints x dim_vv) = "//vtoa(size(mat_vv, 1))//" x "//vtoa(size(mat_vv, 2)))
      call log_info("ML:  dim_vv = d_p × nmax = "//vtoa(this%d_p)//" × "//vtoa(this%nmax)//" = "//vtoa(dim_vv))
      
      ii = 0
      do mua = 1, this%dim_mu
      do muj = 1, this%dim_mu
      do ll = 0, this%lmax
        ikount = 0
        do kk = 1, this%kmax
          ikount = ikount + 1
          ii = this%ic_rad(mua, muj, kk, ll)
          
          if (ikount == 1) then
            call radial_spline(ii)%init(this%r_cut_in, this%r_cut_out, this%npoints, xgrid, yfunc, d_yfunc)
            call init_mat_params_chemmap(this, ii, xgrid, mat_vv, d_mat_vv)
            
            ! Skip normalization for CHEMMAP_HSVD - chemical embedding provides proper scaling
            ! do mm = 1, dim_vv
            !   snorm(mm) = sqrt(sum(mat_vv(:,mm)**2))
            !   if (snorm(mm) <= 1.0d-14) snorm(mm) = 1.0d0
            !   mat_vv(:,mm) = mat_vv(:,mm) / snorm(mm)
            !   d_mat_vv(:,mm) = d_mat_vv(:,mm) / snorm(mm)
            ! end do
            
            ! SVD
            call svd_chem%init(size(mat_vv, 1), size(mat_vv, 2))
            call svd_chem%evaluate(mat_vv)
            
            if (this%kmax > svd_chem%rank) then
              call log_warning(NAMECURRENT//'ML: CHEMMAP_HSVD: kmax > matrix rank')
              call log_warning('kmax='//vtoa(this%kmax)//' rank='//vtoa(svd_chem%rank)// &
                               ' ll='//vtoa(ll)//' mua='//vtoa(mua)//' muj='//vtoa(muj))
            end if
            
            vk_local(1:this%kmax, 1:dim_vv) = svd_chem%eigvecVT(1:this%kmax, 1:dim_vv)
            sigmak_local(1:this%kmax) = svd_chem%eigval(1:this%kmax)
          else
            call radial_spline(ii)%init(this%r_cut_in, this%r_cut_out, this%npoints, xgrid, yfunc, d_yfunc)
          end if
          
          ! Project onto basis (no sigma normalization - chemical embedding provides scaling)
          do ip = 1, size(mat_vv, 1)
            dtmp1 = 0.0_dp
            dtmp2 = 0.0_dp
            do mm = 1, dim_vv
              dtmp1 = dtmp1 + vk_local(kk, mm) * mat_vv(ip, mm)
              dtmp2 = dtmp2 + vk_local(kk, mm) * d_mat_vv(ip, mm)
            end do
            yfunc(ip) = dtmp1
            d_yfunc(ip) = dtmp2
          end do
          
          call radial_spline(ii)%compute(yfunc, d_yfunc)
        end do ! kk
      end do ! ll
      end do ! muj
      end do ! mua

    end if

    _MLD_END_
  end subroutine build_ace_radial

  ! ------------------------------------------------------------------
  ! Random-projection chemical contraction
  !   (kACE Sec. III.F ; Darby et al. PRL 131, 028001, eqs (36)+(38)).
  !
  ! Fills vk_local(1:kmax,1:dim_vv) with generic random directions
  !     W^k_{mu n l}  ~  N(0, 1/K),      K = kmax
  ! over the joint radial-chemical index mm = (mu, n, params). These
  ! random weights REPLACE the deterministic HSVD contraction vectors
  ! t^(k) (paper eq. 40); the caller therefore applies NO singular-value
  ! rescaling (sigma_k is set to 1).  The number of channels K = kmax
  ! is the single accuracy knob: the reconstruction error decays as
  ! 1/sqrt(K) (Johnson-Lindenstrauss).
  !
  ! The entries are generated deterministically from the channel labels
  ! (kk, mm, mua, muj, ll) through a splitmix64 hash + Box-Muller
  ! transform, so they are reproducible run-to-run and identical on
  ! every MPI rank (each channel is owned/allreduced by a single rank).
  ! ------------------------------------------------------------------
  subroutine fill_random_projection_weights(kmax, dim_vv, mua, muj, ll, vk_local)
    implicit none
    integer, intent(in) :: kmax, dim_vv, mua, muj, ll
    real(dp), dimension(:,:), intent(inout) :: vk_local
    integer :: kk, mm
    real(dp) :: u1, u2, inv_sqrt_k
    real(dp), parameter :: two_pi = 6.283185307179586_dp

    inv_sqrt_k = 1.0_dp / sqrt(real(max(kmax, 1), dp))
    do kk = 1, kmax
      do mm = 1, dim_vv
        ! two independent uniforms from distinct hash streams
        u1 = rp_uniform(rp_hash(kk, mm, mua, muj, 2*ll + 1))
        u2 = rp_uniform(rp_hash(kk, mm, mua, muj, 2*ll + 2))
        u1 = max(u1, 1.0e-15_dp)                              ! avoid log(0)
        ! Box-Muller: standard normal N(0,1), then scale to N(0, 1/K)
        vk_local(kk, mm) = sqrt(-2.0_dp * log(u1)) * cos(two_pi * u2) * inv_sqrt_k
      end do
    end do
  end subroutine fill_random_projection_weights

  ! splitmix64-style hash of five integer labels -> 64-bit integer
  pure integer(kind=8) function rp_hash(a, b, c, d, e) result(h)
    integer, intent(in) :: a, b, c, d, e
    integer(kind=8), parameter :: GOLDEN = -7046029254386353131_8  ! 0x9E3779B97F4A7C15
    h = int(a, 8) + GOLDEN
    h = rp_mix( ieor(h, int(b, 8) * GOLDEN) )
    h = rp_mix( ieor(h, int(c, 8) * GOLDEN) )
    h = rp_mix( ieor(h, int(d, 8) * GOLDEN) )
    h = rp_mix( ieor(h, int(e, 8) * GOLDEN) )
  end function rp_hash

  ! splitmix64 finalising mix
  pure integer(kind=8) function rp_mix(x) result(z)
    integer(kind=8), intent(in) :: x
    z = x
    z = ieor(z, ishft(z, -30)) * (-4658895280553007687_8)  ! *0xBF58476D1CE4E5B9
    z = ieor(z, ishft(z, -27)) * (-7723592293110705685_8)  ! *0x94D049BB133111EB
    z = ieor(z, ishft(z, -31))
  end function rp_mix

  ! map a 64-bit hash to a uniform real in (0,1) using its low 53 bits
  pure real(dp) function rp_uniform(h) result(u)
    integer(kind=8), intent(in) :: h
    integer(kind=8) :: m
    m = iand(h, 9007199254740991_8)                         ! 2^53 - 1
    u = (real(m, dp) + 0.5_dp) / 9007199254740992.0_dp      ! / 2^53
  end function rp_uniform

  ! ------------------------------------------------------------------
  ! Randomized (Halko/Martinsson-Rokhlin) truncated SVD for the HSVD /
  ! BLOCK_HSVD radial-basis matrices AA (M x N), M = npoints,
  ! N = dim_vv. build_ace_radial only ever needs the k_target=kmax
  ! dominant right singular vectors/values (the left factor U is
  ! discarded downstream), yet this SVD is invoked once per
  ! (mua,muj,ll) triple, i.e. O(S^2*(lmax+1)) times for S chemical
  ! species -- exact dgesvd('A','A') additionally wastes O(M^2 N)
  ! flops forming the full M x M left factor that is never used.
  !
  ! Here we sketch the column space of AA with a random Gaussian test
  ! matrix Omega (N x l), l = kmax + oversample, optionally sharpen it
  ! with power iterations, orthonormalize the sketch Q (M x l), and
  ! diagonalize only the small projected matrix B = Q^T AA (l x N).
  ! Cost drops from O(M N min(M,N)) to O(M N l) with l << min(M,N).
  !
  ! The Gaussian sketch is generated from a deterministic hash of
  ! (mua,muj,ll), reusing the splitmix64 + Box-Muller machinery of
  ! fill_random_projection_weights above, so results are reproducible
  ! run-to-run (each (mua,muj,ll) task is computed by exactly one MPI
  ! rank, so no cross-rank consistency is required).
  ! ------------------------------------------------------------------
  subroutine randomized_svd_topk(AA, k_target, oversample, n_power_iter, seed1, seed2, seed3, &
                                  sigma_out, VT_out, rank_out)
    implicit none
    real(dp), dimension(:,:), intent(in) :: AA
    integer, intent(in) :: k_target, oversample, n_power_iter, seed1, seed2, seed3
    real(dp), dimension(:), intent(out) :: sigma_out
    real(dp), dimension(:,:), intent(out) :: VT_out
    integer, intent(out) :: rank_out

    integer :: Mrow, Ncol, l, it, ii, info, lwork, ktop
    real(dp), dimension(:,:), allocatable :: Omega, Yb, Zb, Bmat, Uhat, VThat
    real(dp), dimension(:), allocatable :: Sval, work, v_ref
    real(dp) :: rcond_local, proj, nr

    Mrow = size(AA, 1)
    Ncol = size(AA, 2)
    l = max(1, min(k_target + max(oversample, 0), Mrow, Ncol))

    allocate(Omega(Ncol, l))
    call fill_gaussian_hash(Ncol, l, seed1, seed2, seed3, Omega)

    allocate(Yb(Mrow, l))
    call dgemm('N', 'N', Mrow, l, Ncol, 1.0_dp, AA, Mrow, Omega, Ncol, 0.0_dp, Yb, Mrow)
    deallocate(Omega)
    call qr_orthonormalize(Yb, Mrow, l)

    if (n_power_iter > 0) then
      allocate(Zb(Ncol, l))
      do it = 1, n_power_iter
        call dgemm('T', 'N', Ncol, l, Mrow, 1.0_dp, AA, Mrow, Yb, Mrow, 0.0_dp, Zb, Ncol)
        call qr_orthonormalize(Zb, Ncol, l)
        call dgemm('N', 'N', Mrow, l, Ncol, 1.0_dp, AA, Mrow, Zb, Ncol, 0.0_dp, Yb, Mrow)
        call qr_orthonormalize(Yb, Mrow, l)
      end do
      deallocate(Zb)
    end if

    ! small projected matrix B = Q^T AA  (l x Ncol)
    allocate(Bmat(l, Ncol))
    call dgemm('T', 'N', l, Ncol, Mrow, 1.0_dp, Yb, Mrow, AA, Mrow, 0.0_dp, Bmat, l)
    deallocate(Yb)

    ! small dense SVD of the projected matrix
    allocate(Uhat(l, l), VThat(l, Ncol), Sval(l))
    lwork = -1
    allocate(work(1))
    call dgesvd('S', 'S', l, Ncol, Bmat, l, Sval, Uhat, l, VThat, l, work, lwork, info)
    lwork = int(work(1)) + 2
    deallocate(work) ; allocate(work(lwork))
    call dgesvd('S', 'S', l, Ncol, Bmat, l, Sval, Uhat, l, VThat, l, work, lwork, info)
    deallocate(work, Uhat, Bmat)

    if (info /= 0) then
      call log_critical('randomized_svd_topk: dgesvd on the projected matrix failed with info = '//vtoa(info))
      stop 'randomized_svd_topk: dgesvd failed'
    end if

    ! deterministic sign convention: project each right singular vector
    ! onto the fixed all-ones reference direction (same convention as
    ! svd_enforce_signs_by_ref in module_svd_small_gen_matrix)
    allocate(v_ref(Ncol))
    v_ref(:) = 1.0_dp
    nr = sqrt(sum(v_ref*v_ref))
    if (nr > 0.0_dp) v_ref = v_ref / nr
    do ii = 1, l
      proj = dot_product(VThat(ii,1:Ncol), v_ref)
      if (proj < 0.0_dp) VThat(ii,:) = -VThat(ii,:)
    end do
    deallocate(v_ref)

    rcond_local = dabs(Sval(1)) * 1.e-14_dp
    rank_out = 0
    do ii = 1, l
      if (dabs(Sval(ii)) > rcond_local) rank_out = rank_out + 1
    end do

    ktop = min(k_target, l)
    sigma_out(1:ktop) = Sval(1:ktop)
    VT_out(1:ktop, :) = VThat(1:ktop, :)
    if (l < k_target) then
      sigma_out(l+1:k_target) = 0.0_dp
      VT_out(l+1:k_target, :) = 0.0_dp
    end if
    deallocate(Sval, VThat)
  end subroutine randomized_svd_topk

  ! Fill an (Nrow x l) matrix with N(0,1) entries generated from a
  ! deterministic hash of (seed1,seed2,seed3,row,col) -- reproducible
  ! run-to-run, no shared RNG state needed across MPI ranks/tasks.
  subroutine fill_gaussian_hash(Nrow, l, seed1, seed2, seed3, Omega)
    implicit none
    integer, intent(in) :: Nrow, l, seed1, seed2, seed3
    real(dp), dimension(:,:), intent(out) :: Omega
    integer :: irow, jcol, combined
    real(dp) :: u1, u2
    real(dp), parameter :: two_pi = 6.283185307179586_dp

    do jcol = 1, l
      do irow = 1, Nrow
        combined = (jcol-1) * Nrow + irow
        u1 = rp_uniform(rp_hash(seed1, seed2, seed3, combined, 1))
        u2 = rp_uniform(rp_hash(seed1, seed2, seed3, combined, 2))
        u1 = max(u1, 1.0e-15_dp)
        Omega(irow, jcol) = sqrt(-2.0_dp * log(u1)) * cos(two_pi * u2)
      end do
    end do
  end subroutine fill_gaussian_hash

  ! In-place QR orthonormalization: X (Mx x l), Mx >= l, becomes the Q
  ! factor of its QR decomposition (orthonormal columns spanning the
  ! same range as the input X).
  subroutine qr_orthonormalize(X, Mx, l)
    implicit none
    integer, intent(in) :: Mx, l
    real(dp), dimension(:,:), intent(inout) :: X
    real(dp), dimension(:), allocatable :: tau, work
    integer :: lwork, info

    allocate(tau(l))
    lwork = -1
    allocate(work(1))
    call dgeqrf(Mx, l, X, Mx, tau, work, lwork, info)
    lwork = int(work(1)) + 2
    deallocate(work) ; allocate(work(lwork))
    call dgeqrf(Mx, l, X, Mx, tau, work, lwork, info)
    deallocate(work)

    lwork = -1
    allocate(work(1))
    call dorgqr(Mx, l, l, X, Mx, tau, work, lwork, info)
    lwork = int(work(1)) + 2
    deallocate(work) ; allocate(work(lwork))
    call dorgqr(Mx, l, l, X, Mx, tau, work, lwork, info)
    deallocate(work, tau)
  end subroutine qr_orthonormalize

  pure real(dp) function atanh_safe(x) result(y)
    real(dp), intent(in) :: x
    real(dp), parameter  :: eps = 1.0e-12_dp
    real(dp) :: xx
    xx = max(-1.0_dp+eps, min(1.0_dp-eps, x))
    y  = 0.5_dp * log( (1.0_dp+xx) / (1.0_dp-xx) )
  end function atanh_safe

  subroutine init_mat_params(this, ii, xgrid, delta_local, mat_vv, d_mat_vv)
    use module_chemical_species, only: fix_type_to_periodic, periodic_table_element
    implicit none 
    class(ace_radial), intent(inout) :: this
    integer, intent(in) :: ii
    real(dp), dimension(:), intent(in) :: xgrid
    real(dp), dimension(:,:), intent(inout) :: mat_vv, d_mat_vv 
    real(dp), dimension(:,:), intent(in) :: delta_local 
    real(dp), dimension(:), allocatable :: vv, d_vv 
    integer :: ip, mua, muj, nn, ll, e_nn, e_mu, e_params, ivect, idx_mu_pair
    real(dp) :: alpha, beta, radius 
    real(dp) :: Z_j, Z_a, r_Z_j, r_Z_a, rr, ffnn, dffnn, V_nl, dV_nl_dr 

    ! ivect = 0
    ! do mua = 1, this%dim_mu
    !   do muj = 1, mua
    !     ivect = ivect + 1
    !     mu_pair(ivect, 1) = mua
    !     mu_pair(ivect, 2) = muj
    !   end do 
    ! end do

    allocate(vv(size(mat_vv, 2)))
    allocate(d_vv(size(mat_vv, 2)))
    mua = this%dico_radial(ii)%mua
    muj = this%dico_radial(ii)%muj
    idx_mu_pair = (max(mua, muj)-1) * max(mua, muj) / 2 + min(mua, muj)

    ll = this%dico_radial(ii)%ll
    Z_j = periodic_table_element(fix_type_to_periodic(this%dico_radial(ii)%muj))%Z
    Z_a = periodic_table_element(fix_type_to_periodic(this%dico_radial(ii)%mua))%Z
    r_Z_j = periodic_table_element(fix_type_to_periodic(this%dico_radial(ii)%muj))%covalent_radius / 100.0_dp ! in Angstrom
    r_Z_a = periodic_table_element(fix_type_to_periodic(this%dico_radial(ii)%mua))%covalent_radius / 100.0_dp ! in Angstrom
    do ip = 1, this%npoints
      rr = xgrid(ip)
      !TODOace radbody  
      call PolyRad%evaluate(rr)
      !TODOace: nmax 0 or 1. 
      ivect = 0 
      
      do e_nn = 1, this%nmax
        ffnn = PolyRad%radial(e_nn)
        dffnn = PolyRad%d_radial(e_nn)
        do e_params = 1, size(this%alpha_beta_params)
          alpha = this%alpha_beta_params(e_params)%alpha
          beta = this%alpha_beta_params(e_params)%beta
          radius = this%alpha_beta_params(e_params)%radius  ! radius not used in V_nl_and_derivative
          call V_nl_and_derivative(rr, Z_j, Z_a, r_Z_j, r_Z_a, e_nn, ll, alpha, beta, radius, V_nl, dV_nl_dr)
          !call V_nl_and_derivative_cut(rr, Z_j, Z_a, r_Z_j, r_Z_a, e_nn, ll, alpha, beta, radius, V_nl, dV_nl_dr)
          !TODhea!
          do e_mu = 1, this%dim_mu !Here version mu=mu_j
            ivect = ivect + 1 
            vv(ivect) = delta_local(e_mu, muj) * V_nl * ffnn  
            d_vv(ivect) = delta_local(e_mu, muj) * (dV_nl_dr * ffnn + V_nl * dffnn) 
            !vv(ivect) =  V_nl * ffnn  
            !d_vv(ivect) =  (dV_nl_dr * ffnn + V_nl * dffnn) 
          end do   !e_mu
          ! do e_mu = 1, this%dim_mu * (this%dim_mu + 1) / 2
          !   ivect = ivect + 1 
          !   vv(ivect) = delta_local(e_mu, idx_mu_pair) * V_nl * ffnn  
          !   d_vv(ivect) = delta_local(e_mu, idx_mu_pair) * (dV_nl_dr * ffnn + V_nl * dffnn) 
          !   ! vv(ivect) = delta_local(e_mu, idx_mu_pair) * ffnn  
          !   ! d_vv(ivect) = delta_local(e_mu, idx_mu_pair) * dffnn
          !   !vv(ivect) =  V_nl * ffnn  
          !   !d_vv(ivect) =  (dV_nl_dr * ffnn + V_nl * dffnn) 
          ! end do   !e_mu
          !TODhea!
        end do !e_params 
      end do  !e_nn  
      mat_vv(ip, :) = vv(:)
      d_mat_vv(ip, :) = d_vv(:)
    end do 
    deallocate(vv, d_vv)
    return 
  end subroutine init_mat_params

  subroutine init_alpha_beta_params_lin(this)
    class(ace_radial), intent(inout) :: this
    integer :: iparams, no_alpha, no_beta, no_rad, no_of_parameters
    real(dp) :: alpha_min, alpha_max, beta_min, beta_max, rad_min, rad_max
    real(dp), dimension(:), allocatable  :: alpha, beta, rad 
    integer :: ialpha, ibeta, irad 
    ! if (allocated(this%alpha_beta_params)) deallocate(this%alpha_beta_params)
    ! allocate(this%alpha_beta_params(no_of_parameters))
    !old no_alpha= 5
    !old alpha_min = 0.05 !A--1
    !old alpha_max = 0.4 !A -1
    !old no_beta = 5
    !old beta_min = 0.3
    !old beta_max = 1.0

    no_alpha= 5
    alpha_min = 0.05_dp !A--1
    alpha_max = 0.2_dp !A -1
    no_beta = 5
    beta_min = 0.05_dp
    beta_max = 0.4_dp
    ! no_alpha= 5
    ! alpha_min = 0.1 !A--1
    ! alpha_max = 2 !A -1
    ! no_beta = 1
    ! beta_min = 0.5
    ! beta_max = 0.5
    no_rad = 1 
    rad_min = 1.0_dp 
    rad_max = 1.0_dp
    no_of_parameters = no_alpha * no_beta * no_rad 

  if (allocated(alpha)) deallocate(alpha) ; allocate(alpha(no_alpha))
  do ialpha = 1, no_alpha
    if (no_alpha > 1) then
      alpha(ialpha) = alpha_min + (alpha_max - alpha_min) * (dble(ialpha) - 1.0_dp) / dble(no_alpha - 1)
    else
      alpha(ialpha) = alpha_min
    end if
  end do   
  if (allocated(beta)) deallocate(beta) ; allocate(beta(no_beta))
  do ibeta = 1, no_beta
      if (no_beta > 1) then
        beta(ibeta) = beta_min + (beta_max - beta_min) * (dble(ibeta) - 1.0_dp) / dble(no_beta - 1)
      else
        beta(ibeta) = beta_min
      end if
  end do

  if (allocated(rad)) deallocate(rad) ; allocate(rad(no_rad)) 
    do irad = 1, no_rad
      if (no_rad > 1) then
        rad(irad) = rad_min + (rad_max - rad_min) * (dble(irad) - 1.0_dp) / dble(no_rad - 1)
      else
        rad(irad) = rad_min
      end if
    end do

    if (allocated(this%alpha_beta_params)) deallocate(this%alpha_beta_params)
    allocate(this%alpha_beta_params(no_of_parameters))
    iparams = 0
    do ialpha = 1, no_alpha
       do ibeta = 1, no_beta
         do irad = 1, no_rad       
         iparams = iparams + 1 
         this%alpha_beta_params(iparams)%alpha = alpha(ialpha) 
         this%alpha_beta_params(iparams)%beta = beta(ibeta)
         this%alpha_beta_params(iparams)%radius = rad(irad) 
        end do 
      end do
    end do   
  end subroutine init_alpha_beta_params_lin

  !===============================================================
  ! Initialize chemical embedding for ACE_CHEM_RADIAL_CHEMMAP_HSVD
  ! Computes SVD-based linear chemical embedding with:
  ! - W0: isotropic chemical map (d_e x dim_mu)
  ! - W_aniso: anisotropic chemical map (d_e x dim_mu x d_p)
  ! Based on atomic numbers Z_mu
  !===============================================================
  ! Section 1.6: Abstract Construction of Isotropic and Anisotropic Chemical Maps
  ! Following the LaTeX specification exactly
  !===============================================================
  subroutine init_chemical_embedding(this)
    use module_chemical_species, only: fix_type_to_periodic, periodic_table_element
    use module_svd_small_gen_matrix, only: eigen_svd
    implicit none
    class(ace_radial), intent(inout) :: this
    integer :: mu, i, p, d_prop, N_s
    real(dp), dimension(:,:), allocatable :: U_matrix  ! N_s x d_prop property matrix
    real(dp), dimension(:,:), allocatable :: U_c       ! centered property matrix
    real(dp), dimension(:), allocatable :: u_mean      ! mean property vector
    real(dp), dimension(:,:), allocatable :: V_matrix  ! d_prop x d_prop from SVD
    type(eigen_svd) :: svd_chem
    
    ! Following LaTeX: N_s = number of species, d_prop = number of properties
    N_s = this%dim_mu
    d_prop = 6  ! Using 6 properties as in LaTeX example: r, chi, E_coh, Z, M, rho
    
    ! Build U matrix (N_s x d_prop) - each row is u(μ)^T
    if (allocated(U_matrix)) deallocate(U_matrix)
    allocate(U_matrix(N_s, d_prop))
    
    ! Fill property matrix - using available properties from periodic table
    ! Note: Using simple properties for now. TODO: Add more physical properties
    do mu = 1, N_s
      U_matrix(mu, 1) = periodic_table_element(fix_type_to_periodic(mu))%Z      ! atomic number
      U_matrix(mu, 2) = periodic_table_element(fix_type_to_periodic(mu))%Z**2   ! Z squared
      U_matrix(mu, 3) = periodic_table_element(fix_type_to_periodic(mu))%mass   ! atomic mass
      U_matrix(mu, 4) = periodic_table_element(fix_type_to_periodic(mu))%density ! g/cm^3
      U_matrix(mu, 5) = periodic_table_element(fix_type_to_periodic(mu))%ionization ! eV
      U_matrix(mu, 6) = log(real(periodic_table_element(fix_type_to_periodic(mu))%mass, dp)) ! log(mass)
    end do
    
    ! Compute column-wise mean ū
    if (allocated(u_mean)) deallocate(u_mean)
    allocate(u_mean(d_prop))
    do i = 1, d_prop
      u_mean(i) = sum(U_matrix(:, i)) / real(N_s, dp)
    end do
    
    ! Center the matrix: U_c = U - ū (subtract mean from each row)
    if (allocated(U_c)) deallocate(U_c)
    allocate(U_c(N_s, d_prop))
    do mu = 1, N_s
      U_c(mu, :) = U_matrix(mu, :) - u_mean(:)
    end do
    
    call log_info("ML: Property matrix U: "//vtoa(N_s)//" species x "//vtoa(d_prop)//" properties")
    
    ! Perform SVD: U_c = Q Σ V^T
    call svd_chem%init(N_s, d_prop)
    call svd_chem%evaluate(U_c)
    
    call log_info("ML: SVD  CHEM MAPS rank = "//vtoa(svd_chem%rank))
    if (svd_chem%rank >= 1) call log_info("ML: Singular values: "//vtoa(svd_chem%eigval(1:min(3,svd_chem%rank))))
    
    ! Extract V^T matrix (first d_e rows of V^T from SVD)
    ! In eigen_svd, eigvecVT contains V^T
    ! Store V_matrix as d_e x d_prop for computing e(μ) = V^T (u(μ) - ū)
    if (allocated(this%V_matrix)) deallocate(this%V_matrix)
    allocate(this%V_matrix(this%d_e, d_prop))
    this%V_matrix = 0.0_dp
    do i = 1, min(this%d_e, svd_chem%rank)
      this%V_matrix(i, :) = svd_chem%eigvecVT(i, 1:d_prop)
    end do
    
    ! Store mean vector for embedding computation
    if (allocated(this%u_mean)) deallocate(this%u_mean)
    allocate(this%u_mean(d_prop))
    this%u_mean = u_mean
    
    if (allocated(V_matrix)) deallocate(V_matrix)
    allocate(V_matrix(d_prop, d_prop))
    V_matrix = 0.0_dp
    do i = 1, min(svd_chem%rank, d_prop)
      V_matrix(i, :) = svd_chem%eigvecVT(i, 1:d_prop)
    end do
    
    ! Build W0 (isotropic map): W0 = [I_{d_p x d_e}, I_{d_p x d_e}]
    ! This gives c_iso(μ_a, μ_j) = e(μ_a) + e(μ_j)
    if (allocated(this%W0)) deallocate(this%W0)
    allocate(this%W0(this%d_p, 2*this%d_e))  ! d_p x 2*d_e for concatenated z vector
    this%W0 = 0.0_dp
    
    ! W0 = [I, I] block structure (d_p rows, each picking up e_i from both blocks)
    do i = 1, min(this%d_p, this%d_e)
      this%W0(i, i) = 1.0_dp           ! First block: I
      this%W0(i, this%d_e + i) = 1.0_dp  ! Second block: I
    end do
    
    ! Build W_aniso (anisotropic map): uses principal directions v_k^T
    ! W_aniso encodes chemical contrast e(μ_a) - e(μ_j)
    ! Structure: [v_s1^T, -v_s1^T; v_s2^T, -v_s2^T; ...; 0, 0]
    ! Each slice W_aniso(:,:,p) is d_p × 2d_e
    if (allocated(this%W_aniso)) deallocate(this%W_aniso)
    allocate(this%W_aniso(this%d_p, 2*this%d_e, this%d_p))
    this%W_aniso = 0.0_dp
    
    ! Fill anisotropic map with principal directions
    if (this%d_p > 0) then
      do p = 1, min(this%d_p, svd_chem%rank)
        do i = 1, min(this%d_p, this%d_e)
          ! v_p^T row in first block, -v_p^T in second block to capture contrast
          this%W_aniso(i, i, p) = V_matrix(p, i)           ! +v_p^T for e(μ_a)
          this%W_aniso(i, this%d_e + i, p) = -V_matrix(p, i) ! -v_p^T for e(μ_j)
        end do
      end do
    end if
    
    call log_info("ML: Chemical embedding initialized (LaTeX Section 1.6)")
    call log_info("ML:   N_s (species) = "//vtoa(N_s))
    call log_info("ML:   d_prop (properties) = "//vtoa(d_prop))
    call log_info("ML:   d_e (embedding dim) = "//vtoa(this%d_e))
    call log_info("ML:   d_p (principal dirs) = "//vtoa(this%d_p))
    call log_info("ML:   gamma_l (ang. coupling) = "//vtoa(this%gamma_l))
    call log_info("ML:   W0 isotropic: "//vtoa(this%d_p)//" x "//vtoa(2*this%d_e))
    call log_info("ML:   W_aniso: "//vtoa(this%d_p)//" x "//vtoa(2*this%d_e)//" x "//vtoa(this%d_p))
    
  end subroutine init_chemical_embedding

  !===============================================================
  ! Initialize matrix parameters using chemical embedding
  ! For ACE_CHEM_RADIAL_CHEMMAP_HSVD
  ! Implements R_{(n,p)}^{l μ_a μ_j}(r) = P_n(r) c_p^{(l)}(μ_a, μ_j)
  ! where c^{(l)}(μ_a, μ_j) = W^{(l)} z(μ_a, μ_j)
  ! with W^{(l)} = γ_l W0 + (1-γ_l) W_aniso
  !===============================================================
  subroutine init_mat_params_chemmap(this, ii, xgrid, mat_vv, d_mat_vv)
    use module_chemical_species, only: fix_type_to_periodic, periodic_table_element
    implicit none
    class(ace_radial), intent(inout) :: this
    integer, intent(in) :: ii
    real(dp), dimension(:), intent(in) :: xgrid
    real(dp), dimension(:,:), intent(inout) :: mat_vv, d_mat_vv
    real(dp), dimension(:), allocatable :: vv, d_vv
    real(dp), dimension(:), allocatable :: e_mua, e_muj  ! embeddings e(μ_a), e(μ_j)
    real(dp), dimension(:), allocatable :: z_pair        ! concatenated [e(μ_a); e(μ_j)]
    real(dp), dimension(:), allocatable :: c_l           ! c^{(l)}(μ_a, μ_j)
    real(dp), dimension(:,:), allocatable :: W_l         ! W^{(l)} for this l
    real(dp), dimension(:), allocatable :: u_mua, u_muj, u_mean ! property vectors
    integer :: ip, mua, muj, nn, ll, e_nn, e_p, ivect, i, d_prop
    real(dp) :: rr, ffnn, dffnn, gamma_l
    
    allocate(vv(size(mat_vv, 2)))
    allocate(d_vv(size(mat_vv, 2)))
    
    mua = this%dico_radial(ii)%mua
    muj = this%dico_radial(ii)%muj
    ll = this%dico_radial(ii)%ll
    
    ! Compute chemical embeddings e(μ_a) and e(μ_j)
    ! e(μ) = V^T (u(μ) - ū)
    d_prop = 6  ! must match init_chemical_embedding
    allocate(u_mua(d_prop), u_muj(d_prop), u_mean(d_prop))
    allocate(e_mua(this%d_e), e_muj(this%d_e))
    
    ! Get property vectors (same as init_chemical_embedding)
    u_mua(1) = periodic_table_element(fix_type_to_periodic(mua))%Z
    u_mua(2) = periodic_table_element(fix_type_to_periodic(mua))%Z**2
    u_mua(3) = periodic_table_element(fix_type_to_periodic(mua))%mass
    u_mua(4) = periodic_table_element(fix_type_to_periodic(mua))%density
    u_mua(5) = periodic_table_element(fix_type_to_periodic(mua))%ionization
    u_mua(6) = log(real(periodic_table_element(fix_type_to_periodic(mua))%mass, dp))
    
    u_muj(1) = periodic_table_element(fix_type_to_periodic(muj))%Z
    u_muj(2) = periodic_table_element(fix_type_to_periodic(muj))%Z**2
    u_muj(3) = periodic_table_element(fix_type_to_periodic(muj))%mass
    u_muj(4) = periodic_table_element(fix_type_to_periodic(muj))%density
    u_muj(5) = periodic_table_element(fix_type_to_periodic(muj))%ionization
    u_muj(6) = log(real(periodic_table_element(fix_type_to_periodic(muj))%mass, dp))
    
    ! Compute embeddings e(μ) = V^T (u(μ) - ū)
    ! e_mua = V_matrix * (u_mua - u_mean)
    e_mua = matmul(this%V_matrix, u_mua - this%u_mean)
    e_muj = matmul(this%V_matrix, u_muj - this%u_mean)
    
    ! Build concatenated pair vector z(μ_a, μ_j) = [e(μ_a); e(μ_j)]
    allocate(z_pair(2*this%d_e))
    z_pair(1:this%d_e) = e_mua
    z_pair(this%d_e+1:2*this%d_e) = e_muj
    
    ! Compute W^{(l)} = γ_l W0 + (1-γ_l) W_aniso for this l
    gamma_l = exp(this%gamma_l * real(ll, dp))  ! Example: exponential decay
    allocate(W_l(this%d_p, 2*this%d_e))
    W_l = 0.0_dp
    
    ! W^{(l)} = γ_l W0 + (1-γ_l) W_aniso
    ! For l=0: pure isotropic (e(μ_a) + e(μ_j))
    ! For l>0: mix of isotropic and anisotropic (includes contrast e(μ_a) - e(μ_j))
    W_l = gamma_l * this%W0
    if (ll > 0 .and. ll <= this%d_p) then
      ! Use the anisotropic map for this specific angular channel
      W_l = W_l + (1.0_dp - gamma_l) * this%W_aniso(:, :, ll)
    else if (ll > this%d_p) then
      ! For l > d_p, use modulo to wrap around
      W_l = W_l + (1.0_dp - gamma_l) * this%W_aniso(:, :, mod(ll-1, this%d_p) + 1)
    end if
    
    ! Compute c^{(l)}(μ_a, μ_j) = W^{(l)} z(μ_a, μ_j)
    allocate(c_l(this%d_p))
    c_l = matmul(W_l, z_pair)
    
    ! Build radial channels: R_{(n,p)}^{l μ_a μ_j}(r) = P_n(r) c_p^{(l)}(μ_a, μ_j)
    do ip = 1, this%npoints
      rr = xgrid(ip)
      call PolyRad%evaluate(rr)
      
      ivect = 0
      do e_nn = 1, this%nmax
        ffnn = PolyRad%radial(e_nn)
        dffnn = PolyRad%d_radial(e_nn)
        
        do e_p = 1, this%d_p
          ivect = ivect + 1
          ! R_{(n,p)}^{l μ_a μ_j}(r) = P_n(r) c_p^{(l)}(μ_a, μ_j)
          vv(ivect) = ffnn * c_l(e_p)
          d_vv(ivect) = dffnn * c_l(e_p)
        end do
      end do
      
      mat_vv(ip, :) = vv(:)
      d_mat_vv(ip, :) = d_vv(:)
    end do
    
    deallocate(vv, d_vv, e_mua, e_muj, z_pair, c_l, W_l, u_mua, u_muj, u_mean)
  end subroutine init_mat_params_chemmap

!===============================================================
! Log-spaced grid for (alpha, beta) parameter pairs
! - alpha is spaced logarithmically between [alpha_min, alpha_max]
! - beta is spaced logarithmically between [beta_min, beta_max]
! - Handles edge cases: no_alpha=1 or no_beta=1, swapped bounds,
! and non-positive minima (clamped to tiny).
! - Keep units consistent with your potential:
! * Exponential-linear form: alpha in 1/Angstrom (Å^-1)
! * Gaussian form: alpha in 1/Angstrom^2 (Å^-2)
!===============================================================
subroutine init_alpha_beta_params_log(this)
class(ace_radial), intent(inout) :: this
integer :: iparams, no_alpha, no_beta, no_of_parameters
real(dp) :: alpha_min, alpha_max, beta_min, beta_max
real(dp) :: alpha, beta, t
real(dp) :: log_amin, log_amax, log_bmin, log_bmax
integer :: ialpha, ibeta
real(dp), parameter :: tiny = 1.0e-12_dp


! --- User-set ranges (examples) ---
! For linear-exponential potential (alpha in Å^-1):
! no_alpha = 5; alpha_min = 0.05_dp; alpha_max = 0.20_dp
! For Gaussian potential (alpha in Å^-2), e.g.:
! no_alpha = 5; alpha_min = 0.20_dp; alpha_max = 2.00_dp


no_alpha = 5
alpha_min = 0.0001_dp
alpha_max = 0.20_dp


no_beta = 5
beta_min = 0.05_dp
beta_max = 0.40_dp


! --- Sanity / robustness ---
if (alpha_min <= 0.0_dp) alpha_min = tiny
if (beta_min <= 0.0_dp) beta_min = tiny
if (alpha_max < alpha_min) then
t = alpha_min; alpha_min = alpha_max; alpha_max = t
end if
if (beta_max < beta_min) then
t = beta_min; beta_min = beta_max; beta_max = t
end if


log_amin = log(alpha_min)
log_amax = log(alpha_max)
log_bmin = log(beta_min)
log_bmax = log(beta_max)


no_of_parameters = no_alpha * no_beta
if (allocated(this%alpha_beta_params)) deallocate(this%alpha_beta_params)
allocate(this%alpha_beta_params(no_of_parameters))


iparams = 0
do ialpha = 1, no_alpha
if (no_alpha > 1) then
t = real(ialpha-1, dp) / real(no_alpha-1, dp)
alpha = exp( log_amin + t * (log_amax - log_amin) )
else
alpha = alpha_min
end if


do ibeta = 1, no_beta
if (no_beta > 1) then
t = real(ibeta-1, dp) / real(no_beta-1, dp)
beta = exp( log_bmin + t * (log_bmax - log_bmin) )
else
beta = beta_min
end if


iparams = iparams + 1
this%alpha_beta_params(iparams)%alpha = alpha
this%alpha_beta_params(iparams)%beta = beta
end do
end do


end subroutine init_alpha_beta_params_log


  ! Here version mu=mu_j
  subroutine V_nl_and_derivative(r_ja, Z_j, Z_a, r_Z_j, r_Z_a, n, l, alpha0, beta, radius, V_nl, dV_nl_dr)
    implicit none
    integer, intent(in) ::  n, l
    double precision, intent(in) :: r_ja, alpha0, beta, r_Z_j, r_Z_a, Z_j, Z_a, radius 
    double precision, intent(out) :: V_nl, dV_nl_dr
    double precision :: alpha, prefactor, exp_term, term1, term2, r0, r_safe, r0_safe

    r0 = radius * (r_Z_a + r_Z_j)
    r0_safe = max(r0, 1.0e-12_dp)
    r_safe = max(r_ja, 1.0e-12_dp)

    alpha = alpha0 * (Z_j**beta + Z_a**beta) / dble(n + l + 1) !v0
    prefactor = (Z_j * Z_a) * (r_safe / r0_safe)**l / 1000.0_dp !v0

    term1 = l / r_safe
    exp_term = exp(-alpha * (r_safe - r0_safe) ** 2)
    term2 = - 2.0_dp * alpha * (r_safe - r0_safe)

    V_nl = prefactor * exp_term
    ! write(*,*) 'V_nl: ', V_nl
    ! derivative terms
    dV_nl_dr = prefactor * exp_term * (term1 + term2)

    ! V_nl = log(1.0_dp + exp(prefactor * exp_term)) + 1.e-6! softplus(prefactor * exp_term)
    ! ! derivative terms
    ! dV_nl_dr = (prefactor * exp_term * (term1 + term2)) / (1.0_dp + exp(-prefactor * exp_term)) ! derivative of softplus(prefactor * exp_term)

    ! write(*,*) prefactor, exp_term, alpha0, alpha, beta, r_ja, r0, term1, term2, V_nl, dV_nl_dr
    ! write(*,*) alpha, r_ja - r0, r_ja / r0, exp_term, exp(-alpha * (r_ja / r0) ** 2)
    ! write(*,*) exp_term, exp(-alpha * (r_ja / r0) ** 2), prefactor, V_nl, dV_nl_dr
    return
  end subroutine V_nl_and_derivative

  !===============================================================
  ! Decreasing power form with smooth cutoff P_n(r) and no r0
  ! V(r) = [(Z_j Z_a)/Z_scale] * (r_cut/r)^l * exp(-alpha_nl * r) * P(r)
  ! - alpha_nl = alpha0 * (Z_j^beta + Z_a^beta) / (n + l)
  ! - P(r) : quintic smoothstep on [0, r_cut], 0 beyond r_cut (C^2 smooth)
  ! - Units: r, r_cut in Å; alpha0 in Å^-1; radii if used elsewhere in Å
  ! - Numerical guard: eps prevents r=0 singularity
  !===============================================================
  subroutine V_nl_and_derivative_cut(r, Z_j, Z_a, r_Z_j, r_Z_a, n, l, alpha0, beta, radius,  V, dVdr)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  integer, intent(in) :: n, l
  real(dp), intent(in) :: r, alpha0, beta, radius, Z_j, Z_a, r_Z_j, r_Z_a
  real(dp), intent(out) :: V, dVdr
  real(dp) :: alpha_nl, prefactor, powfac, dpow_dr
  real(dp) :: exp_term, dexp_dr, fc, dfc_dr, x, eps
  real(dp) :: Zj_beta, Za_beta,  r_cut, Z_scale
  
  r_cut = 5.0_dp  ! in Angstrom
  Z_scale = 1000.0_dp  ! to keep magnitudes reasonable
  eps = 1.0e-12_dp
  
  
  ! --- alpha_nl (Å^-1) ---
  Zj_beta = real(Z_j, dp)**beta
  Za_beta = real(Z_a, dp)**beta
  alpha_nl = alpha0 * (Zj_beta + Za_beta) / real(n + l, dp)
  
  
  ! --- chemical prefactor scaled to keep magnitudes reasonable ---
  prefactor = ( real(Z_j,dp) * real(Z_a,dp) / Z_scale )
  
  
  ! --- decreasing power: (r_cut / r)^l ---
  powfac = ( r_cut / max(r, eps) )**l
  dpow_dr = - real(l,dp) * powfac / max(r, eps)
  
  
  ! --- monotone exponential tail exp(-alpha * r) ---
  exp_term = exp( -alpha_nl * r )
  dexp_dr = -alpha_nl * exp_term
  
  
  ! --- quintic cutoff P(r) on [0, r_cut] ---
  if (r < r_cut) then
  x = r / r_cut
  fc = 1.0_dp - 10.0_dp*x**3 + 15.0_dp*x**4 - 6.0_dp*x**5
  dfc_dr = ( -30.0_dp*x**2 + 60.0_dp*x**3 - 30.0_dp*x**4 ) / r_cut
  else
  fc = 0.0_dp
  dfc_dr = 0.0_dp
  end if
  
  
  ! --- value and derivative (product rule) ---
  V = prefactor * powfac * exp_term * fc
  dVdr = prefactor * ( dpow_dr * exp_term * fc &
  + powfac * dexp_dr * fc &
  + powfac * exp_term * dfc_dr )
  
  
  end subroutine V_nl_and_derivative_cut

  ! subroutine V_nl_and_derivative(r_ja, Z_j, Z_a, r_Z_j, r_Z_a, n, l, alpha0, beta, V_nl, dV_nl_dr)
  !   implicit none
  !   integer, intent(in) ::  n, l
  !   double precision, intent(in) :: r_ja, alpha0, beta, r_Z_j, r_Z_a, Z_j, Z_a
  !   double precision, intent(out) :: V_nl, dV_nl_dr
  !   double precision :: alpha, prefactor, exp_term, term1, term2, r0 

  !   r0 = r_Z_a + r_Z_j
  !   alpha = alpha0 / dble(n + l)
  !   prefactor = beta * (r_ja / r0)**l 
  !   term1 = l / r_ja

  !   ! exp_term = exp(-alpha * (r_ja - r0) ** 2)
  !   ! term2 = - 2.0_dp * alpha * (r_ja - r0)
  !   exp_term = exp(-alpha * (r_ja - r0) ** 2)
  !   term2 = - 2.0_dp * alpha * (r_ja - r0)

  !   V_nl = prefactor * exp_term
  !   ! derivative terms
  !   dV_nl_dr = prefactor * exp_term * (term1 + term2)

  !   ! write(*,*) prefactor, exp_term, alpha0, alpha, beta, r_ja, r0, term1, term2, V_nl, dV_nl_dr
  !   ! write(*,*) alpha, r_ja - r0, r_ja / r0, exp_term, exp(-alpha * (r_ja / r0) ** 2)
  !   ! write(*,*) exp_term, exp(-alpha * (r_ja / r0) ** 2), prefactor, V_nl, dV_nl_dr
  !   return
  ! end subroutine V_nl_and_derivative

  !> Helper: MPI_Reduce to owner_rank, in-place on root.
  !! Isolates MPI_IN_PLACE in a separate routine so GNU gfortran f08
  !! strict type checking does not cross-compare sendbuf types across branches.
  subroutine mld_reduce_to_owner(buf, nelem, owner_rank, icomm, my_rank)
    use mpi, only: MPI_Reduce, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
    implicit none
    real(8), intent(inout) :: buf(*)
    integer, intent(in)    :: nelem, owner_rank, my_rank, icomm
    integer :: ierr
    if (my_rank == owner_rank) then
      call MPI_Reduce(MPI_IN_PLACE, buf, nelem, MPI_DOUBLE_PRECISION, MPI_SUM, &
                      owner_rank, icomm, ierr)
    else
      call MPI_Reduce(buf, buf, nelem, MPI_DOUBLE_PRECISION, MPI_SUM, &
                      owner_rank, icomm, ierr)
    end if
  end subroutine mld_reduce_to_owner

end module module_ace_radial

module module_ace_big_radial
  use mld_logger
  use mld_mpi
  use mld_string, only: vtoa
  implicit none 

  contains 


  subroutine init_ace_radial_basis 
      use module_kind_variables, only: kind_double
      use module_ace_desc, only: init_mode_ace_radial
!$!
!$!      use module_body_desc, only: dim_desc_body
!$!      use module_ftnbody, only: dim_qbody, init_mode_ftnbody, length_order, & 
!$!                                covar_matrix_4b, count_4b,  &
!$!                                length_rff_4b,    &
!$!                                ftnbody_omega_4b,&
!$!                                ftnbody_phase_random_4b, & 
!$!                                covar_matrix_4b, number_4b_ftnb, mean_4b , & 
!$!                                ltmp_mean_ftnbody
!$!      use module_kernel, only:  krff_type
!$!      use module_sample_rand_ker, only: get_rff_diag_sigma
!$!      use my_mpi_subroutines, only: subworlds_allreduce_from_evrywhere_double, &
!$!                                    subworlds_allreduce_from_evrywhere_vect_double, & 
!$!                                    subworlds_allreduce_from_evrywhere_matrix_double
!$!      implicit none
!$!      integer :: iiii, ii 
!$!      real(kind_double) :: mu
!$!      real(kind_double), allocatable, dimension(:) :: total_mean 
!$!      real(kind_double), allocatable, dimension(:,:) :: total_cov 
!$!      character(len=80) :: chstring, string  
!$!  
      _NAMECURRENT_("init_ace_radials_basis")
      _MLD_BEGIN_
!$!  
!$!      ! First step computing the means ...........................
      init_mode_ace_radial = .true.
!$!      ltmp_mean_ftnbody = .true. 
      call log_info("ML: ... ACE begins, in "//NAMECURRENT//" estimation of radial function: screen the database " )
      call main_compute_descriptors()
!$!
!$!  
!$!  
!$!      !----------if (l_body_order(4)) then
!$!        number_4b_ftnb = 0.d0 
!$!        call subworlds_allreduce_from_evrywhere_double(count_4b, number_4b_ftnb)
!$!  
!$!        if (allocated(total_mean)) deallocate(total_mean) 
!$!        allocate(total_mean(size(mean_4b,1)))
!$!  
!$!        call subworlds_allreduce_from_evrywhere_vect_double(mean_4b, total_mean)
!$!  
!$!        mean_4b = total_mean / number_4b_ftnb 
!$!  
!$!        call log_info("ML: ... 4b number ftnbody  : "//vtoa(number_4b_ftnb))
!$!        call log_info("ML: ... 4b mean   ftnbody  : "//vtoa(mean_4b))
!$!  
!$!  
!$!        count_4b = 0 
!$!        covar_matrix_4b(:,:)=0.d0 
!$!      !-----------end if 
!$!  
!$!      !  step computing the covariance matrix  ...........................
!$!      call log_info("ML: ... ftnbody mean estimated , in "//NAMECURRENT//" estimation of length, step 2. the covariance " )
!$!      init_mode_ftnbody = .true.
!$!      ltmp_mean_ftnbody = .false. 
!$!      call main_compute_descriptors()    
!$!      init_mode_ftnbody = .false. 
!$!      ltmp_mean_ftnbody = .false. 
!$!  
!$!      call log_info("ML: ... ftnbody completed the covariance " )
!$!  
!$!  
!$!      !------------if (l_body_order(4)) then
!$!        if (allocated(total_cov)) deallocate(total_cov) 
!$!        allocate(total_cov(size(covar_matrix_4b,1), size(covar_matrix_4b,2)))
!$!  
!$!        call subworlds_allreduce_from_evrywhere_matrix_double(covar_matrix_4b, total_cov)
!$!  
!$!        covar_matrix_4b(:,:) = total_cov(:,:) / number_4b_ftnb
!$!  
!$!        !covar_matrix_4b(:,:) = covar_matrix_4b(:,:) / SQRT(count_4b)
!$!  
!$!        do iiii = 1, dim_qbody(4)
!$!          length_rff_4b(iiii) = SQRT(covar_matrix_4b(iiii,iiii))
!$!        end do
!$!        length_rff_4b(:) = length_rff_4b(:)/2
!$!  
!$!        call log_info("ML: ... 4b length   ftnbody  : "//vtoa(length_rff_4b))
!$!  
!$!        length_order = 4
!$!        call  get_rff_diag_sigma(krff_type, mu, length_rff_4b, dim_qbody(4), dim_desc_body(4), ftnbody_omega_4b, ftnbody_phase_random_4b)
!$!      !-------------end if
!$!  
!$!  
      _MLD_END_
  
    end subroutine init_ace_radial_basis
  
end module module_ace_big_radial 
  
  
module module_AB_basis 
  use iso_fortran_env, only: dp => real64
  use mld_logger 
  implicit none 
  integer, parameter :: dp_complex = dp
  ! matrix of ylm computed for atom ja (neigh_ja%max_neigh, ilm)
  complex(dp_complex), dimension(:,:), allocatable  :: matylm 
  ! the same thing but derivatives (3, neigh_ja%max_neigh, ilm)
  complex(dp_complex), dimension(:,:,:), allocatable  :: d_matylm
  ! map l and m into ilm.  
  integer, dimension(:,:), allocatable :: lm_vec
  real(dp), dimension(:,:), allocatable :: rsp, drsp
  
  !---> matAa unique A basis for a given neighbourhood (uniqueA)  
  complex(dp_complex), dimension(:), allocatable :: matAa
  !---> d_matAa unique A basis for a given neighbourhood (3, neigh_ja%max_neigh, uniqueA)
  !nnew complex(dp_complex), dimension(:,:,:), allocatable :: d_matAa  
  !nnew 
  !---> d_matAa2 unique A basis for a given neighbourhood (3*neigh_ja%max_neigh, uniqueA)
  complex(dp_complex), dimension(:,:), allocatable :: d_matAa2  

  contains 

  subroutine compute_rad_for_neigh_ja(nu)
  ! compute the radial functions for the neighbours of atom ja.
  ! the radial functions are stored in:
  !   rsp(ia_n, is) of dimension #neighbours_ja x #(radial_spline)
  !   drsp(ia_n, is) the derivatives of rsp(ia_n, is) with respect to r
    use module_neigh_local, only: neigh_ja
    use module_ace_radial, only: radial_spline
    integer, intent(in) :: nu
    integer :: ia_n, is 
    real(dp) :: rr
    
    if (allocated(rsp)) deallocate(rsp) ; allocate(rsp(neigh_ja%max_neigh_local, size(radial_spline)))
    if (allocated(drsp)) deallocate(drsp) ; allocate(drsp(neigh_ja%max_neigh_local, size(radial_spline)))

    do ia_n = 1, neigh_ja%max_neigh_local
      rr = neigh_ja%r_central(ia_n)
      do is = 1, size(radial_spline)
        rsp(ia_n, is) = radial_spline(is)%evaluate(rr)
        drsp(ia_n, is) = radial_spline(is)%derivative(rr)
      end do
    end do

  end subroutine compute_rad_for_neigh_ja

  subroutine compute_ylm_for_neigh_ja (nu) 
  ! compute the Ylm for the neighbours of atom ja.
  ! the Ylm are stored in:
  !   matylm(ia_n, ilm) of dimension #neighbours_ja x #(ll,mm) combinations 
  !   d_matylm(3, ia_n, ilm) the derivatives with x_ja, y-ja, z_ja 
  !                          of dimension 3 x #neighbours_ja x #(ll,mm) combinations
  ! ilm is the index of the Ylm in the matrix lm_vec(ll, mm)
  ! the Ylm are computed for all the neighbours of atom ja.

    use module_neigh_local, only: neigh_ja
    use angular_functions, only: spherical_harm, grad_spherical_harm
    use module_ace_desc, only: base_params
    use mld_mpi, only: mld_rank
    integer, intent(in) :: nu 
    integer :: lmax, icount, ic1, ic2,  mm, ll, ia_n  
    real(dp), dimension(3)  :: dxp 
    complex(dp_complex) :: ctmp 
    complex(dp_complex) :: d_ctmp(3) 

    lmax = base_params(nu)%lmax 
  
    if (allocated(lm_vec)) deallocate(lm_vec)
    allocate(lm_vec(0:lmax, -lmax:lmax))
    !is orverestimated: there are (lmax+1)**2 Ylms and the matrix has more  (lmax+1)*(2*lmax+1)
    icount = 0 
    do ll = 0, lmax 
       do mm = 0,0
         icount = icount + 1  
         lm_vec(ll, mm) = icount  
       end do 

       do mm = 1, ll
         icount = icount + 1
         lm_vec(ll, mm) = icount
         icount = icount + 1 
         lm_vec(ll, -mm) = icount
       end do 
    end do    

    if (allocated(matylm)) deallocate(matylm) 
    allocate(matylm(neigh_ja%max_neigh_local, icount))

    if (allocated(d_matylm)) deallocate(d_matylm) 
    allocate(d_matylm(3,neigh_ja%max_neigh_local, icount))

    icount = 0 
    do ll = 0, lmax
     do mm = 0,0
        icount = icount + 1  
        do ia_n = 1, neigh_ja%max_neigh_local
          dxp(:) = neigh_ja%tmp_dxp(:,ia_n)
          ctmp = spherical_harm(ll, mm, dxp)
          d_ctmp(1:3) = grad_spherical_harm(ll, mm, dxp(1:3))
          matylm(ia_n, icount) = ctmp
          d_matylm(1:3, ia_n, icount) = d_ctmp(1:3)
        end do  
      end do 

      do mm = 1, ll
        ic1 = icount + 1
        ic2 = icount + 2  
        do ia_n = 1, neigh_ja%max_neigh_local
        dxp(:) = neigh_ja%tmp_dxp(:,ia_n)
        ctmp = spherical_harm(ll, mm, dxp)
        d_ctmp(1:3) = grad_spherical_harm(ll, mm, dxp(1:3))
        matylm(ia_n, ic1) = ctmp
        d_matylm(1:3, ia_n, ic1) = d_ctmp(1:3) 
        matylm(ia_n,ic2) = (-1.0_dp)**mm * conjg(ctmp) 
        d_matylm(1:3, ia_n, ic2) = (-1.0_dp)**mm * conjg(d_ctmp(1:3))
        end do  
        icount = icount + 2 
      end do 
    end do   
    
  end subroutine compute_ylm_for_neigh_ja

  subroutine compute_B_basis_atom(nu, baseB, d_baseB, desc_forces) 
    !use mld_mpi, only: mld_rank
    use module_base_cnlm, only: type_base_cnlm, base_cnlm
    use module_neigh_local, only: neigh_ja
    use time_check_general, only: time, debug_time, MY_MPI_WTIME !, tot_time 
    use module_ace_desc, only: time_amat, time_amat_00, time_amat_01, time_amat_02, time_amat_03, time_cgord
    integer, intent(in) :: nu 
    logical, intent(in) :: desc_forces
    real(dp), dimension(:), allocatable, intent(inout) :: baseB
    real(dp), dimension(:,:,:), allocatable, intent(inout) :: d_baseB

    integer :: ib, im, ia_n, ix 
    !integer, dimension(nu) :: nbold, mbold, lbold, mubold
    integer :: idx_nb, idx_mub, idx_lb0, idx_Li, idx_cg, dimr_M0mat, dim_baseB, icount, dim_munl   
    !integer :: idx_nblb0 
    complex(dp_complex) :: Aa 
    ! complex(dp_complex), dimension(:,:), allocatable  :: d_Aa 
    complex(dp_complex), dimension(:), allocatable :: Am  
    !nnew complex(dp_complex), dimension(:,:,:), allocatable :: d_Am
    !nnew   
    complex(dp_complex), dimension(:,:), allocatable :: d_Am2  
    real(dp), dimension(:), allocatable  :: cgm
    !old_target integer, dimension(:,:), pointer :: ptr_M0mat
    !old_target integer, dimension(:), pointer :: ptr_vec 
    !old_target integer, target, dimension(:,:), allocatable  :: target_M0mat
    !old_target integer,  dimension(:), allocatable  :: target_vec 

    integer :: k, inu, index_offset, nn3, sb, iLb 
    complex(dp_complex), dimension(nu) :: Atmp, partial_Aa
    
    !nnew complex(dp_complex),  dimension(:,:,:), allocatable :: d_Atemp
    !nnew 
    complex(dp_complex),  dimension(:,:), allocatable :: d_Atemp2
    complex(dp_complex), dimension(:), allocatable  :: dxnn
    complex(dp_complex) :: right_partial
    complex(dp_complex) :: c_one, c_zero  
    real(dp) :: tt01, tt02, tt03 
    ! complex(dp_complex) ::     baseB_2
    _NAMECURRENT_("compute_B_basis_atom")
    _MLD_BEGIN_ 

    nn3 = 3*neigh_ja%max_neigh_local
    c_one = cmplx((1.0_dp, 0.0_dp), kind=dp_complex)
    c_zero = cmplx((0.0_dp, 0.0_dp), kind=dp_complex)
    !nnew allocate(d_Atemp(3, neigh_ja%max_neigh_local, nu))
    !nnew 
    allocate(d_Atemp2(3*neigh_ja%max_neigh_local, nu))
    allocate(dxnn(3*neigh_ja%max_neigh_local))
    !integer :: idx_nb       ! this map into dico_nb 
    !integer :: idx_mub      ! this map into dico_mu
    !integer :: idx_lb0      ! this map into dico_lb0
    !integer :: idx_munblb0  ! this map into dico_mu_nnll 
    !integer :: idx_Li       ! this map into or dico_lb0%cg(:,i)  or dico_mub_nb_lb0%rpicg(:,i)
    dim_baseB = size(base_cnlm(nu)%dicoB)
    dim_munl = size(base_cnlm(nu)%dico_mub_nb_lb0)
    allocate(baseB(dim_baseB))
    allocate(d_baseB(dim_baseB, 3, neigh_ja%max_neigh_local))
    icount = 0 
    ib = 0 
    do sb = 1, dim_munl
      !idx_nb  = base_cnlm(nu)%dicoB(ib)%idx_nb
      !nbold = base_cnlm(nu)%dico_nb(idx_nb)%nbold
      idx_nb = base_cnlm(nu)%dico_mub_nb_lb0(sb)%idx_nb

      !idx_mub = base_cnlm(nu)%dicoB(ib)%idx_mub 
      !mubold = base_cnlm(nu)%dico_mu(idx_mub)%mubold 
      idx_mub = base_cnlm(nu)%dico_mub_nb_lb0(sb)%idx_mub

      !idx_lb0 = base_cnlm(nu)%dicoB(ib)%idx_lb0 
      !lbold = base_cnlm(nu)%dico_lb0(idx_lb0)%lbold 
      idx_lb0 = base_cnlm(nu)%dico_mub_nb_lb0(sb)%idx_lb0

      idx_cg = base_cnlm(nu)%dico_mub_nb_lb0(sb)%idx_cg


      if (debug_time) then
            time(10) = MY_MPI_WTIME()
      end if 

      dimr_M0mat = size(base_cnlm(nu)%dico_lb0(idx_lb0)%M0mat, 1) 
      allocate(Am(dimr_M0mat))
      allocate(d_Am2(dimr_M0mat, nn3))

      if (debug_time) then
            tt01 = MY_MPI_WTIME()
            time_amat_00 = time_amat_00 + tt01 - time(10)
      end if

      do im = 1, dimr_M0mat 
        if (debug_time) then
            tt01 = MY_MPI_WTIME()
        end if 
        icount = icount + 1 
        index_offset = (icount-1)*nu
        Aa = c_one  
        do inu = 1, nu
           k = base_cnlm(nu)%map_large_to_uniqueA(index_offset + inu) 
           Atmp(inu) =  matAa(k) 
           if (desc_forces) d_Atemp2(:, inu) = d_matAa2(:,k)
           partial_Aa (inu) = Aa  
           Aa = Aa * Atmp(inu)
        end do

        if (debug_time) then
            tt02 = MY_MPI_WTIME()
            time_amat_01 = time_amat_01 + tt02 - tt01
        end if

        right_partial = c_one    
        do inu = nu, 1, -1
            partial_Aa(inu) = partial_Aa(inu) * right_partial
            right_partial = right_partial * Atmp(inu)
        end do
        Am(im) = Aa

        if (debug_time) then
            tt03 = MY_MPI_WTIME()
            time_amat_02 = time_amat_02 + tt03 - tt02
        end if 

        !$! dxnn(:) = c_zero
        !$! do inu = 1, nu
        !$!    dxnn(:) = dxnn(:) + partial_Aa(inu) * d_Atemp2(:, inu)
        !$! end do
        !$! d_Am2(im, :) = dxnn(:) 
        !! this is lapack version of above. 
        if (desc_forces) then
        call zgemv('N', nn3, nu, c_one, d_Atemp2, nn3, partial_Aa, 1, c_zero, dxnn, 1)
        d_Am2(im, :) = dxnn(:)
        end if 

        if (debug_time) then
            tt02 = MY_MPI_WTIME()
            time_amat_03 = time_amat_03 + tt02 - tt03
        end if 
        
      end do   

      if (debug_time) then
            time(11) = MY_MPI_WTIME()
            time_amat = time_amat + time(11) - time(10)
      end if 

      allocate(cgm(dimr_M0mat))
      do iLb = 1, size(base_cnlm(nu)%dicoCG(idx_cg)%cg,2) 
         ib = ib + 1
         idx_Li = base_cnlm(nu)%dicoB(ib)%idx_Li
         cgm = base_cnlm(nu)%dicoCG(idx_cg)%cg(:,idx_Li)   
         baseB(ib) = real(DOT_PRODUCT(Am(:), cgm(:)), kind=dp) 
        !  baseB_2 = DOT_PRODUCT(Am(:), cgm(:))
        !  write(*,*) "version basis. n=", base_cnlm(nu)%dico_nb(idx_nb)%nbold(:), ", l=", base_cnlm(nu)%dico_lb0(idx_lb0)%lbold(:), ", baseB(ib) real / imaginary ", baseB_2%re, baseB_2%im

         if (desc_forces) then
         do ia_n = 1, neigh_ja%max_neigh_local
           do ix = 1,3 
           d_baseB(ib, ix, ia_n) = real(DOT_PRODUCT(d_Am2(:, ix + 3*(ia_n-1)), cgm(:)), kind=dp) 
           end do 
         end do
         else 
           d_baseB(ib, 1:3, :) = 0.0_dp
         end if 
      end do     
      deallocate(Am)
      deallocate(d_Am2)
      deallocate(cgm)
      if (debug_time) then
            tt03  = MY_MPI_WTIME()  
            time_cgord = time_cgord + tt03 - time(11)
      end if
    end do 
    !debug if (mld_rank==0) write(*,'("aaa", 3f18.10)') sum(d_baseB(:, 1, :)**2), sum(d_baseB(:, 2, :)**2), sum(d_baseB(:, 3, :)**2)

    _MLD_END_ 
  end subroutine compute_B_basis_atom

  subroutine compute_B_basis_atom_adjoint(nu, baseB, d_baseB, desc_forces) 
     use mld_mpi, only: mld_mpi_abort ! , mld_rank
     use module_base_cnlm, only: type_base_cnlm, base_cnlm
     use module_neigh_local, only: neigh_ja
     use time_check_general, only: time, debug_time, MY_MPI_WTIME ! , tot_time
     use module_ace_desc, only: time_amat, time_amat_00, time_amat_01, time_amat_02, time_amat_03, time_cgord, ace_kmax
     integer, intent(in) :: nu 
     logical, intent(in) :: desc_forces
     real(dp), dimension(:), allocatable, intent(inout) :: baseB
     real(dp), dimension(:,:,:), allocatable, intent(inout) :: d_baseB
     integer :: ib, ibase, sb, im, ia_n, ix, kk, num_unique
     integer :: idx_nb, idx_mub, idx_lb0, idx_Li, idx_cg, iLb, dimr_M0mat, dim_baseB, dim_munl, icount  
     !integer :: idx_nblb0 
     complex(dp_complex) :: Aa 
     complex(dp_complex), dimension(:), allocatable :: Am  
    !  complex(dp_complex), dimension(:,:), allocatable :: d_Am   ! , d_Amn 
     real(dp), dimension(:), allocatable :: rAm  
    !  real(dp), dimension(:,:), allocatable :: d_rAm

     real(dp), dimension(:), allocatable  :: cgm
     integer :: k, inu, index_offset, dim_k,nn3 
     integer :: acekmax_nu
     complex(dp_complex), dimension(nu) :: Atmp, partial_Aa  
     complex(dp_complex),  dimension(:), allocatable :: d_Atemp
     real(dp), dimension(:), allocatable ::  omega
     !vers_02 real(dp), dimension(:), allocatable ::  ptmpI, ptmpR
     !vers_02 real(dp), allocatable :: pA_re(:,:), pA_im(:,:), pM_re(:,:), pM_im(:,:)
     !vers_02 integer :: pnrow

     ! real(dp) :: dtmpBk
     complex(dp_complex) :: cdtmpBk
     
    !  real(dp), dimension(:), allocatable  :: dxnn
     complex(dp_complex), dimension(:), allocatable  :: dxnn
    !  complex(dp_complex), dimension(:), allocatable  :: cdxnn
     real(dp), dimension(:), allocatable  :: rdxnn
     real(dp), dimension(:,:), allocatable :: d_Am2  
     complex(dp_complex),  dimension(:,:), allocatable :: d_Atemp2


     complex(dp_complex) :: right_partial
     complex(dp_complex) :: c_one, c_zero
     real(dp) :: r_one   
     real(dp) :: tt01, tt02, tt03 
     type :: Bbasis_workspace
      complex(dp_complex), allocatable :: tmp(:), zres(:), cgmC(:)
     end type
     type(Bbasis_workspace) :: ws
     integer, dimension(nu) :: inu_to_k 
     integer, parameter :: i_one=1 
    !  complex(dp), allocatable :: y(:)

    !  complex(dp_complex) :: baseB_2
    _NAMECURRENT_("compute_B_basis_atom_adjoint")
    _MLD_BEGIN_ 

     acekmax_nu = ace_kmax(nu)
     nn3 = 3*neigh_ja%max_neigh_local
     c_one = cmplx((1.0_dp, 0.0_dp), kind=dp_complex)
     c_zero = cmplx((0.0_dp, 0.0_dp), kind=dp_complex)
     r_one = 1.0_dp
     dim_k = size(matAa) 
    !  write(*,*) "debug acekmax_nu, dim_k, nu, num_unique ", acekmax_nu, dim_k, nu
     num_unique = int(dim_k/acekmax_nu)
     allocate(d_Atemp(dim_k))
     allocate(omega(dim_k))
     allocate(dxnn(nn3))
     allocate(rdxnn(nn3))
    !  allocate(cdxnn(nn3))
     allocate(d_Atemp2(3*neigh_ja%max_neigh_local, nu))

     !integer :: idx_nb       ! this map into dico_nb 
     !integer :: idx_mub      ! this map into dico_mu
     !integer :: idx_lb0      ! this map into dico_lb0
     !integer :: idx_munblb0  ! this map into dico_mu_nnll 
     !integer :: idx_Li       ! this map into dico_lb0%Lmat(i,:) or dico_mub_nb_lb0%rpicg(:,i)
     dim_baseB = size(base_cnlm(nu)%dicoB) * acekmax_nu
     dim_munl = size(base_cnlm(nu)%dico_mub_nb_lb0)

     if (allocated(baseB)) deallocate(baseB) 
     allocate(baseB(dim_baseB))
     if (allocated(d_baseB)) deallocate(d_baseB)
     allocate(d_baseB(dim_baseB, 3, neigh_ja%max_neigh_local))

     
      ibase = 0
      do kk = 1, acekmax_nu
        icount = 0
        ib = 0 
        do sb = 1, dim_munl
          ! I think that this is not used at all 
          ! idx_nb  = base_cnlm(nu)%dico_mub_nb_lb0(sb)%idx_nb
          ! I think that this also is not used at all
          ! idx_mub = base_cnlm(nu)%dico_mub_nb_lb0(sb)%idx_mub 
          ! those two are used pletly often.
          idx_lb0 = base_cnlm(nu)%dico_mub_nb_lb0(sb)%idx_lb0 
          idx_cg =  base_cnlm(nu)%dico_mub_nb_lb0(sb)%idx_cg
          
          if (debug_time) then
                time(10) = MY_MPI_WTIME()
          end if 
          
          dimr_M0mat = size(base_cnlm(nu)%dico_lb0(idx_lb0)%M0mat, 1) 
          allocate(Am(dimr_M0mat))
          ! allocate(d_Am(dimr_M0mat, dim_k))
          allocate(d_Am2(dimr_M0mat, nn3))

          if (debug_time) then
                tt01 = MY_MPI_WTIME()
                time_amat_00 = time_amat_00 + tt01 - time(10)
          end if

          do im = 1, dimr_M0mat 
            if (debug_time) then
                tt01 = MY_MPI_WTIME()
            end if
            icount  = icount + 1 
            index_offset = (icount-1)*nu
            do inu = 1, nu 
              inu_to_k (inu) = base_cnlm(nu)%map_large_to_uniqueA(index_offset + inu) + (kk-1)*num_unique
            end do   
            Aa = c_one  
            do inu = 1, nu
                !k = base_cnlm(nu)%map_large_to_uniqueA(index_offset + inu) + (kk-1)*num_unique
                k = inu_to_k(inu) 
                Atmp(inu) =  matAa(k) 
                if (desc_forces) d_Atemp2(:, inu) = d_matAa2(:,k)
                partial_Aa (inu) = Aa  
                Aa = Aa * Atmp(inu)
            end do
    
            if (debug_time) then
                tt02 = MY_MPI_WTIME()
                time_amat_01 = time_amat_01 + tt02 - tt01
            end if
      
            right_partial = c_one    
            do inu = nu, 1, -1
                partial_Aa(inu) = partial_Aa(inu) * right_partial
                right_partial = right_partial * Atmp(inu)
            end do
      
            Am(im) = Aa
            if (debug_time) then
                tt03 = MY_MPI_WTIME()
                time_amat_02 = time_amat_02 + tt03 - tt02
            end if 
    
            ! if (desc_forces) then
            !   !omega(1:dim_k) = 0.0_dp
            !   d_Atemp(1:dim_k) = c_zero
            !   do inu = 1, nu
            !     !k = base_cnlm(nu)%map_large_to_uniqueA(index_offset + inu) + (kk-1)*num_unique   
            !     k = inu_to_k(inu)     
            !     !omega(k) = omega(k) + 1.0_dp 
            !     d_Atemp(k) = d_Atemp(k) + partial_Aa(inu)
            !   end do 
            !   d_Am(im, 1:dim_k) = d_Atemp(1:dim_k)
            ! end if 



            dxnn(:) = c_zero
            rdxnn(:) = 0.d0

            if (desc_forces) then
              do inu = 1, nu
                k = inu_to_k(inu) 
                dxnn(:) = dxnn(:) + partial_Aa(inu) * d_matAa2(:, k)
                ! call zaxpy(nn3, partial_Aa(inu), d_matAa2(:, k), 1, dxnn, 1)
                ! dxnn(:) = dxnn(:) + cdxnn(:)
              end do
              d_Am2(im, :) = real(dxnn, kind=dp)
            end if
            !!! this is lapack version of above. 
            ! if (desc_forces) then
            !   call zgemv('N', nn3, nu, c_one, d_Atemp2, nn3, partial_Aa, 1, c_zero, dxnn, 1)
            !   d_Am2(im, :) = dxnn(:)
            ! end if         
    
            if (debug_time) then
              tt02 = MY_MPI_WTIME()
              time_amat_03 = time_amat_03 + tt02 - tt03
            end if 
    
          end do ! im in M0mat   
    
          if (debug_time) then
                time(11) = MY_MPI_WTIME()
                time_amat = time_amat + time(11) - time(10)
          end if 
    
          !allocate(cgm(dimr_M0mat))
          allocate(cgm(dimr_M0mat), source = 0.0_dp)   ! everything = 0.0
          allocate(rAm(dimr_M0mat))
          rAm = real(Am, kind=dp)
          !  allocate(d_rAm(dimr_M0mat, dim_k))
          !  d_rAm = real(d_Am, kind=dp)
          !allocate(d_Amn(dimr_M0mat, nn3))

          !vers_02 pnrow = size(d_matAa2, 1)
          !vers_02 ! -- allocate and copy real/imag parts once --
          !vers_02 allocate(pA_re(dimr_M0mat,dim_k));  pA_re = real(d_Am)
          !vers_02 allocate(pA_im(dimr_M0mat,dim_k));  pA_im = aimag(d_Am)
          !vers_02 allocate(pM_re(pnrow,dim_k));        pM_re = real(d_matAa2)
          !vers_02 allocate(pM_im(pnrow,dim_k));        pM_im = aimag(d_matAa2)
          !vers_02 allocate(ptmpR(dim_k), ptmpI(dim_k))



          allocate(ws%tmp(dim_k), ws%zres(nn3))
          allocate(ws%cgmC(dimr_M0mat))

    
          do iLb = 1, size(base_cnlm(nu)%dicoCG(idx_cg)%cg, 2)
            ibase = ibase + 1
            ib = ib + 1
            idx_Li = base_cnlm(nu)%dicoB(ib)%idx_Li
            cgm = base_cnlm(nu)%dicoCG(idx_cg)%cg(:,idx_Li)  
            !  baseB(ib) = real(DOT_PRODUCT(cgm(:), Am(:)), kind=dp) 
            baseB(ibase) = DOT_PRODUCT(rAm(:), cgm(:))

            ! if (desc_forces) then
            !   do ia_n = 1, neigh_ja%max_neigh_local
            !     do ix = 1,3 
            !     d_baseB(ibase, ix, ia_n) = real(DOT_PRODUCT(d_Am2(:, ix + 3*(ia_n-1)), cgm(:)), kind=dp) 
            !     end do 
            !   end do
            ! else 
            !   d_baseB(ibase, 1:3, :) = 0.0_dp
            ! end if 
            if (desc_forces) then
              call dgemv('T', dimr_M0mat, nn3, 1.0_dp, d_Am2, dimr_M0mat, cgm, 1, 0.0_dp, rdxnn, 1)
              do ia_n = 1, neigh_ja%max_neigh_local
                do ix = 1, 3
                    d_baseB(ibase, ix, ia_n) = rdxnn(ix + 3*(ia_n-1))
                end do
              end do
            else 
              d_baseB(ibase, 1:3, :) = 0.0_dp
            end if 


            ! if (desc_forces) then
            !    cdxnn(:) = c_zero 
            !    do k = 1, dim_k 
            !      !dtmpBk  =  real(DOT_PRODUCT(d_Am(:,k), cgm(:)), kind=dp) 
            !      cdtmpBk  =  DOT_PRODUCT(cgm(:), d_Am(:,k))
            !      !dxnn(:) = dxnn(:) + dtmpBk * omega(k) *real(d_matAa2(:,k), kind=dp) 
            !      cdxnn(:) = cdxnn(:) + cdtmpBk * d_matAa2(:,k)
            !    end do
            !    dxnn(:) = real(cdxnn(:), kind=dp)
               
            !TODhea!
              !good! cdxnn(:) = 0.d0 
              !good! do k = 1, dim_k 
              !good!   cdtmpBk  =  DOT_PRODUCT(cgm(:), d_Am(:,k))
              !good!   cdxnn(:) = cdxnn(:) + cdtmpBk%re * d_matAa2(:,k)%re - cdtmpBk%im * d_matAa2(:,k)%im
              !good! end do
              !good! dxnn(:) = cdxnn(:)
              !!!! call fast_dxnn(nn3, dimr_M0mat, dim_k, cgm, d_Am, d_matAa2, dxnn, ws%cgmC, ws%tmp, ws%zres)

              !whynotwork! ws%cgmC = cmplx(cgm, 0.0_dp, kind=dp)
              !whynotwork! call zgemv('T', dimr_M0mat, dim_k, c_one, d_Am, dimr_M0mat, ws%cgmC, i_one, c_zero, ws%tmp, i_one)
              !whynotwork! ! zres = d_matAa2 * tmp
              !whynotwork! call zgemv('N', nn3, dim_k, c_one, d_matAa2, nn3, ws%tmp, i_one, c_zero, ws%zres, i_one)
              !whynotwork! dxnn = real(ws%zres)


              !vers_02 dxnn = 0.0_dp                         ! clear
              !vers_02 ! 1) tmpR = A_re^T * cgm
              !vers_02 call dgemv('T', dimr_M0mat, dim_k, r_one, pA_re, dimr_M0mat, cgm, 1, 0.0_dp, ptmpR, 1)
              !vers_02 ! 2) tmpI = A_im^T * cgm
              !vers_02 call dgemv('T', dimr_M0mat, dim_k, r_one, pA_im, dimr_M0mat, cgm, 1, 0.0_dp, ptmpI, 1)
              !vers_02 ! 3a)  dxnn =  M_re * tmpR
              !vers_02 call dgemv('N', pnrow, dim_k, r_one,  pM_re, pnrow, ptmpR, 1, 0.0_dp, dxnn, 1)
              !vers_02 ! 3b)  dxnn += -M_im * tmpI
              !vers_02 call dgemv('N', pnrow, dim_k,-1.0_dp,  pM_im, pnrow, ptmpI, 1, 1.0_dp, dxnn, 1)
              !TODhea! 
              
            !   do ia_n = 1, neigh_ja%max_neigh_local
            !     do ix = 1,3 
            !       d_baseB(ibase, ix, ia_n) = dxnn(ix + 3*(ia_n-1))
            !     end do
            !   end do
            ! else 
            !   d_baseB(ibase, 1:3, :) = 0.0_dp
            ! end if   
          end do  ! iLb in dicoCG
    
          deallocate(Am)
          !deallocate(d_Amn)
          deallocate(rAm)
          ! deallocate(d_Am)
          deallocate(d_Am2)
          !  deallocate(d_rAm)
          deallocate(cgm) 
          !vers_02 deallocate(pA_im, pA_re, pM_im, pM_re)
          !vers_02 deallocate(ptmpI, ptmpR)

          deallocate(ws%tmp, ws%zres, ws%cgmC)
          
          if (debug_time) then
                tt03  = MY_MPI_WTIME()  
                time_cgord = time_cgord + tt03 - time(11)
          end if 
        end do ! sb in dico_munl
      end do   ! kk in acekmax_nu
     if (dim_baseB /= ibase) then 
        call log_critical("Error in dimension of compute_B_basis_atom_adjoint")
        call log_critical("dim_baseB  "//vtoa(dim_baseB)//"  icount "//vtoa(ibase))
        call mld_mpi_abort("---ask developpers!---")
     end if
    !debug if (mld_rank==0) write(*,'("bbb", 3f18.10)') sum(d_baseB(:, 1, :)**2), sum(d_baseB(:, 2, :)**2), sum(d_baseB(:, 3, :)**2)


     _MLD_END_ 
  end subroutine compute_B_basis_atom_adjoint

  subroutine fast_dxnn(nrow, dimr, dim_k, cgm, d_Am, d_matAa2_in, dxnn, &
                     work_cgmC, work_tmp, work_zres)
    use iso_fortran_env, only : dp => real64
    implicit none
    integer, intent(in) :: nrow, dimr, dim_k
    real(dp)   , intent(in) :: cgm(dimr)                 ! real
    complex(dp), intent(in) :: d_Am(dimr,dim_k)          ! complex
    complex(dp), intent(in) :: d_matAa2_in(nrow,dim_k)   ! complex
    real(dp)   , intent(out):: dxnn(nrow)                ! real
    ! – reusable work buffers (already allocated) –
    complex(dp), intent(inout) :: work_cgmC(dimr), work_tmp(dim_k), work_zres(nrow)
    complex(dp), parameter :: CZ=(0.0_dp,0.0_dp),  CO=(1.0_dp,0.0_dp)
    integer, parameter :: ONE = 1
    ! real → complex copy once
    work_cgmC = cmplx(cgm, 0.0_dp, kind=dp)
    ! tmp = transpose(d_Am) * cgm   [no conjugate]
    call zgemv('T', dimr, dim_k, CO, d_Am, dimr, work_cgmC, ONE, CZ, work_tmp, ONE)
    ! zres = d_matAa2 * tmp
    call zgemv('N', nrow, dim_k, CO, d_matAa2_in, nrow, work_tmp, ONE, CZ, work_zres, ONE)
    dxnn = real(work_zres)
  end subroutine fast_dxnn


  !unsued_not_erase! subroutine compute_A_basis_atom_old(nu, icount,  Aa, d_Aa)
  !unsued_not_erase!   use module_base_cnlm, only: base_cnlm
  !unsued_not_erase!   use module_neigh_local, only: neigh_ja
  !unsued_not_erase!   integer, intent(in) :: nu, icount 
  !unsued_not_erase!   !integer, dimension(:), intent(in)  :: mubold, nbold, lbold, mbold  
  !unsued_not_erase!   complex(dp_complex), intent(inout)  :: Aa 
  !unsued_not_erase!   complex(dp_complex), dimension(:,:), allocatable, intent(inout)  :: d_Aa
  !unsued_not_erase!   !local variables
  !unsued_not_erase!   integer :: inu,  ia_n, k   
  !unsued_not_erase!   complex(dp_complex), dimension(nu) :: Atmp 
  !unsued_not_erase!   complex(dp_complex), dimension(3) :: cdxp
  !unsued_not_erase!   complex(dp_complex),  dimension(:,:,:), allocatable :: d_Atemp
  !unsued_not_erase!   _NAMECURRENT_("compute_A_basis_atom_new")
  !unsued_not_erase!   _MLD_BEGIN_ 
  !unsued_not_erase!   allocate(d_Aa(3, neigh_ja%max_neigh_local))
  !unsued_not_erase!   allocate(d_Atemp(3, neigh_ja%max_neigh_local, nu))
  !unsued_not_erase!   !Aa = 1.0_dp
  !unsued_not_erase!   !d_Aa= 1.d0 
  !unsued_not_erase!   do inu = 1, nu
  !unsued_not_erase!      k = base_cnlm(nu)%map_large_to_uniqueA((icount-1)*nu + inu) 
  !unsued_not_erase!      Atmp(inu) =  matAa(k)  
  !unsued_not_erase!      d_Atemp(:,:, inu) = d_matAa(:,:,k)  
  !unsued_not_erase!   end do
  !unsued_not_erase!   ! ! 
  !unsued_not_erase!   Aa = 1.0_dp
  !unsued_not_erase!   do inu = 1, nu  
  !unsued_not_erase!     Aa = Aa * Atmp(inu)         
  !unsued_not_erase!   end do 
  !unsued_not_erase!   ! ! 
  !unsued_not_erase!   do ia_n = 1, neigh_ja%max_neigh_local
  !unsued_not_erase!     cdxp(:) = (0.0_dp, 0.0_dp)
  !unsued_not_erase!     do inu = 1, nu
  !unsued_not_erase!       cdxp(:) = cdxp(:) +  Aa / Atmp(inu) *   d_Atemp(:,ia_n, inu)
  !unsued_not_erase!     end do 
  !unsued_not_erase!     d_Aa(:,ia_n) = cd    complex(dp), intent(inout) :: work_cgmC(dimr), work_tmp(dim_k), work_zres(nrow)xp(:)
  !unsued_not_erase!   end do 
  !unsued_not_erase!   _MLD_END_ 
  !unsued_not_erase! end subroutine compute_A_basis_atom_old


  !unused_not_erase! subroutine compute_A_basis_atom_very_old(nu, mubold, nbold, lbold, mbold,  Aa, d_Aa)
  !unused_not_erase!   !use module_base_cnlm, only: type_base_cnlm
  !unused_not_erase!   use module_neigh_local, only: neigh_ja
  !unused_not_erase!   use module_ace_desc, only: delta_ace
  !unused_not_erase!   use module_ace_radial, only: radialace
  !unused_not_erase!   use angular_functions, only: spherical_harm, grad_spherical_harm
  !unused_not_erase!   integer, intent(in) :: nu
  !unused_not_erase!   integer, dimension(:), intent(in)  :: mubold, nbold, lbold, mbold  
  !unused_not_erase!   complex(dp_complex), intent(inout)  :: Aa 
  !unused_not_erase!   complex(dp_complex), dimension(:,:), allocatable, intent(inout)  :: d_Aa
  !unused_not_erase!   !local variables
  !unused_not_erase!   integer :: inu, ispline, type_ja_db, type_ia_db, ia_n 
  !unused_not_erase!   integer :: mu, nn, ll, mm, ilm 
  !unused_not_erase!   complex(dp_complex) :: Atmp_nu, sharm_ll_mm  
  !unused_not_erase!   complex(dp_complex), dimension(nu) :: Atmp 
  !unused_not_erase!   real(dp), dimension(3) :: dxp
  !unused_not_erase!   complex(dp_complex), dimension(3) :: cdxp
  !unused_not_erase!   complex(dp_complex), dimension(3) :: grad_sharm_ll_mm
  !unused_not_erase!   real(dp) :: rad, drad, rr 
  !unused_not_erase!   complex(dp_complex),  dimension(:,:,:), allocatable :: d_Atemp
  !unused_not_erase!   _NAMECURRENT_("compute_A_basis_atom_old")
  !unused_not_erase!   _MLD_BEGIN_ 
  !unused_not_erase!   allocate(d_Aa(3, neigh_ja%max_neigh_local))
  !unused_not_erase!   allocate(d_Atemp(3, neigh_ja%max_neigh_local, nu))
  !unused_not_erase!   type_ja_db = neigh_ja%i_type_db(0)
  !unused_not_erase!   Aa = 1.0_dp 
  !unused_not_erase!   do inu = 1, nu
  !unused_not_erase!     mu = mubold(inu)
  !unused_not_erase!     nn = nbold(inu) 
  !unused_not_erase!     ll = lbold(inu)
  !unused_not_erase!     mm = mbold(inu)
  !unused_not_erase!     ilm = lm_vec(ll, mm)
  !unused_not_erase!     !
  !unused_not_erase!     Atmp_nu=(0.0_dp, 0.0_dp) 
  !unused_not_erase!     do ia_n = 1, neigh_ja%max_neigh_local
  !unused_not_erase!       type_ia_db = neigh_ja%i_type_db(ia_n)
  !unused_not_erase!       rr = neigh_ja%r_central(ia_n)
  !unused_not_erase!       dxp(:) = neigh_ja%tmp_dxp(:, ia_n)
  !unused_not_erase!     
  !unused_not_erase!       ispline = radialace%ic_rad(mu, type_ia_db, nn, ll)
  !unused_not_erase!        
  !unused_not_erase!       !sharm_ll_mm = spherical_harm(ll, mm, dxp(:))
  !unused_not_erase!       sharm_ll_mm = matylm(ia_n, ilm)
  !unused_not_erase!       
  !unused_not_erase!       !rad = radial_spline(ispline)%evaluate(rr)
  !unused_not_erase!       rad = rsp(ia_n, ispline)
  !unused_not_erase!       
  !unused_not_erase!       Atmp_nu = Atmp_nu +  delta_ace(mu, type_ia_db) * rad *sharm_ll_mm
  !unused_not_erase!       
  !unused_not_erase!       !grad_sharm_ll_mm(:) = grad_spherical_harm(ll, mm, dxp(:))
  !unused_not_erase!       grad_sharm_ll_mm(1:3) = d_matylm(1:3, ia_n, ilm)
  !unused_not_erase!       !drad = radial_spline(ispline)%derivative(rr)
  !unused_not_erase!       drad = drsp(ia_n, ispline)
  !unused_not_erase!       
  !unused_not_erase!       d_Atemp(:,ia_n, inu) = delta_ace(mu, type_ia_db)*sharm_ll_mm * drad * dxp(:) / rr + delta_ace(mu, type_ia_db)*rad * grad_sharm_ll_mm(:)   
  !unused_not_erase!     end do    
  !unused_not_erase!     Atmp(inu) = Atmp_nu 
  !unused_not_erase!     Aa = Aa * Atmp_nu         
  !unused_not_erase!   end do 
  !unused_not_erase!   do ia_n = 1, neigh_ja%max_neigh_local
  !unused_not_erase!     cdxp(:) = (0.0_dp, 0.0_dp)
  !unused_not_erase!     do inu = 1, nu
  !unused_not_erase!       cdxp(:) = cdxp(:) +  Aa / Atmp(inu) *   d_Atemp(:,ia_n, inu)
  !unused_not_erase!     end do 
  !unused_not_erase!     d_Aa(:,ia_n) = cdxp(:)
  !unused_not_erase!   end do 
  !unused_not_erase!   _MLD_END_ 
  !unused_not_erase! end subroutine compute_A_basis_atom_very_old 


  subroutine compute_uniqueA_neigh_ja (nu)
    use module_neigh_local, only: neigh_ja
    use module_base_cnlm, only: base_cnlm
    use module_ace_desc, only: delta_ace, ace_chem_low_rank
    use module_ace_radial, only:  radialace, radial_spline, ACE_CHEM_RADIAL_RALF, ACE_CHEM_RADIAL_BLOCK_HSVD, &
                    ACE_CHEM_RADIAL_HSVD, ACE_CHEM_RADIAL_RANDPROJ
    use module_chem_compress, only: chem_compressor
    integer, intent(in) :: nu
    integer :: mu, nn, ll, mm, ia_n, ilm, ispline, type_ia_db, type_ja_db, nuniqueA, k, mu_temp   
    complex(dp_complex) :: Atmp_nu, sharm_ll_mm
    real(dp), dimension(3) :: dxp
    complex(dp_complex), dimension(3) :: grad_sharm_ll_mm
    real(dp) :: rad, drad, rr
    complex(dp_complex) :: c_zero 
    !nnew
    integer :: iin

    c_zero = cmplx((0.0_dp, 0.0_dp), kind=dp_complex)
    nuniqueA= size(base_cnlm(nu)%uniqueA_tuples,2)

    ! if (allocated(matAa)) then
    !   write (*,*) "ML: matAa", matAa
    !   deallocate(matAa) 
    ! end if
    ! allocate(matAA(nuniqueA))
    if (allocated(matAa)) deallocate(matAa) 
    allocate(matAA(nuniqueA))
    !nnew if (allocated(d_matAa)) deallocate(d_matAa) ; allocate(d_matAa(3, neigh_ja%max_neigh_local, nuniqueA))
    !nnew 
    ! write (*,*) "ML: 3*neigh_ja%max_neigh_local", 3*neigh_ja%max_neigh_local, ", nuniqueA=", nuniqueA
    if (allocated(d_matAa2)) deallocate(d_matAa2)
    allocate(d_matAa2(3*neigh_ja%max_neigh_local, nuniqueA))
    !TODhea! 
    type_ja_db = neigh_ja%i_type_db(0)
    !TODhea!


    do k = 1, nuniqueA
      ! For HSVD, keep the species index from the uniqueA tuple to preserve cross-species channels.
      mu = base_cnlm(nu)%uniqueA_tuples(1,k)
      nn = base_cnlm(nu)%uniqueA_tuples(2,k) !this is actually kk when ACE_CHEM_RADIAL_HSVD
      ll = base_cnlm(nu)%uniqueA_tuples(3,k)
      mm = base_cnlm(nu)%uniqueA_tuples(4,k)
  
      
      mu_temp = mu
      if (radialace%type_chem_radial == ACE_CHEM_RADIAL_HSVD .or. &
          radialace%type_chem_radial == ACE_CHEM_RADIAL_RANDPROJ) then
        mu_temp = type_ja_db
      end if


      ilm = lm_vec(ll, mm)  
      Atmp_nu = c_zero
      do ia_n = 1, neigh_ja%max_neigh_local
        type_ia_db = neigh_ja%i_type_db(ia_n)
        rr = neigh_ja%r_central(ia_n)
        dxp(:) = neigh_ja%tmp_dxp(:, ia_n)
        ispline = radialace%ic_rad(mu_temp, type_ia_db, nn, ll)
     
        !sharm_ll_mm = spherical_harm(ll, mm, dxp(:))
        sharm_ll_mm = matylm(ia_n, ilm)    
        !TODhea!
        if (ace_chem_low_rank == 1 .and. radialace%type_chem_radial == ACE_CHEM_RADIAL_HSVD) then
          rad = chem_compressor%eval_rad(mu_temp, type_ia_db, nn, ll, rr)
        else
          rad = radial_spline(ispline)%evaluate(rr)
        end if
        !rad = rsp(ia_n, ispline)
        !TODhea! 

        Atmp_nu = Atmp_nu +  delta_ace(mu_temp, type_ia_db) * rad *sharm_ll_mm
    
        !grad_sharm_ll_mm(:) = grad_spherical_harm(ll, mm, dxp(:))
        grad_sharm_ll_mm(1:3) = d_matylm(1:3, ia_n, ilm)
        !TODhea!
        if (ace_chem_low_rank == 1 .and. radialace%type_chem_radial == ACE_CHEM_RADIAL_HSVD) then
          drad = chem_compressor%eval_drad(mu_temp, type_ia_db, nn, ll, rr)
        else
          drad = radial_spline(ispline)%derivative(rr)
        end if
        !drad = drsp(ia_n, ispline)
        !TODhea!    
        iin = 3*(ia_n-1)
        d_matAa2(iin+1:iin+3, k) = delta_ace(mu_temp, type_ia_db)*(sharm_ll_mm * drad * dxp(1:3) / rr + rad * grad_sharm_ll_mm(1:3))
      end do    
      matAa(k) = Atmp_nu 
    end do 
  end subroutine compute_uniqueA_neigh_ja


  subroutine compute_AB_body_one_neigh_ja(dim_baseBone, baseB, d_baseB, desc_forces_local)
    use ml_in_ndm_module, only: one_pi
    use module_neigh_local, only: neigh_ja
    !use module_base_cnlm, only: base_cnlm
    use module_ace_desc, only: delta_ace, base_params, nkmax_order1, nkmin_order1, ace_chem_low_rank
    use module_ace_radial, only:  radialace, radial_spline, ACE_CHEM_RADIAL_RALF, ACE_CHEM_RADIAL_BLOCK_HSVD, &
                    ACE_CHEM_RADIAL_HSVD, ACE_CHEM_RADIAL_RANDPROJ
    use module_chem_compress, only: chem_compressor
    implicit none 
    real(dp), dimension(:), allocatable, intent(inout) :: baseB
    real(dp), dimension(:,:,:), allocatable, intent(inout) :: d_baseB
    logical, intent(in) :: desc_forces_local
    integer, intent(inout) :: dim_baseBone

    integer :: mu, mu_temp, nn, ll, mm, ia_n, ix, ispline, type_ia_db, type_ja_db, k, kk   
    real(dp) :: Atmp_nu, sharm_ll_mm, rad, drad, rr
    real(dp), dimension(3) :: dxp
    real(dp), dimension(:), allocatable :: mat_dist
    real(dp), dimension(:,:), allocatable :: d_mat_dist
    !nnew
    integer :: iin
    
    ! if (radialace%type_chem_radial == ACE_CHEM_RADIAL_RALF ) then
    !   dim_baseBone = base_params(1)%mumax * (base_params(1)%nmax + 1)
    ! else if (radialace%type_chem_radial == ACE_CHEM_RADIAL_HSVD) then
    !   dim_baseBone = base_params(1)%mumax * (base_params(1)%kmax)
    ! else
    !   call log_warning("ML: in ACE the type_chem_radial is not yet implemented.")
    !   stop 'stop here in compute_AB_body_one_neigh_ja'
    ! end if
    dim_baseBone = base_params(1)%mumax * (nkmax_order1 + 1 - nkmin_order1)
    if (allocated(baseB)) deallocate(baseB) ; allocate(baseB(dim_baseBone))
    if (allocated(d_baseB)) deallocate(d_baseB) ; allocate(d_baseB(dim_baseBone, 3, neigh_ja%max_neigh_local))

    if (allocated(mat_dist)) deallocate(mat_dist) ; allocate(mat_dist(dim_baseBone))
    if (allocated(d_mat_dist)) deallocate(d_mat_dist) ; allocate(d_mat_dist(3*neigh_ja%max_neigh_local, dim_baseBone))
    type_ja_db = neigh_ja%i_type_db(0)

    k = 0 
    do mu = 1, base_params(1)%mumax

      
      ! CORRECT for HSVD nu=1:
      ! For nu=1, base_params(1)%mumax equals the number of species (not 1!).
      ! The loop creates one descriptor per species: mu=1,2,3,4 for 4 species.
      ! mu_temp MUST equal mu to access the correct per-species radials ic_rad(mu,...).
      ! This is different from nu≥2 where uniqueA_tuples has mu=1 always.
      mu_temp = mu   ! CORRECT: use mu from loop (varies by species)
      !mu_temp = 1   ! WRONG: hardcode gives bad RMSE ~2.9 eV/Å
      
      do nn = nkmin_order1, nkmax_order1
        k = k + 1
        ll = 0
        mm = 0  
        Atmp_nu = 0.0_dp 
        do ia_n = 1, neigh_ja%max_neigh_local
          type_ia_db = neigh_ja%i_type_db(ia_n)
          rr = neigh_ja%r_central(ia_n)
          dxp(:) = neigh_ja%tmp_dxp(:, ia_n)
          ispline = radialace%ic_rad(mu_temp, type_ia_db, nn, ll)
      
          !sharm_ll_mm = spherical_harm(ll, mm, dxp(:))
          sharm_ll_mm = 1.0_dp/sqrt(4.0_dp*one_pi)
      
            if (ace_chem_low_rank == 1 .and. radialace%type_chem_radial == ACE_CHEM_RADIAL_HSVD) then
            rad = chem_compressor%eval_rad(mu_temp, type_ia_db, nn, ll, rr)
          else
            rad = radial_spline(ispline)%evaluate(rr)
          end if
          !rad = rsp(ia_n, ispline)
      
          Atmp_nu = Atmp_nu +  delta_ace(mu_temp, type_ia_db) * rad *sharm_ll_mm
      
          !grad_sharm_ll_mm(:) = grad_spherical_harm(ll, mm, dxp(:))
          !grad_sharm_ll_mm(1:3) = d_matylm(1:3, ia_n, ilm)
            if (ace_chem_low_rank == 1 .and. radialace%type_chem_radial == ACE_CHEM_RADIAL_HSVD) then
            drad = chem_compressor%eval_drad(mu_temp, type_ia_db, nn, ll, rr)
          else
            drad = radial_spline(ispline)%derivative(rr)
          end if
          !drad = drsp(ia_n, ispline)
      
          !nnew d_matAa(:,ia_n, k) = delta_ace(mu, type_ia_db)*(sharm_ll_mm * drad * dxp(:) / rr + rad * grad_sharm_ll_mm(:))   
          !nnew 
          iin = 3*(ia_n-1)
          d_mat_dist(iin+1:iin+3, k) = delta_ace(mu_temp, type_ia_db)*(sharm_ll_mm * drad * dxp(1:3) / rr)
        end do    
        mat_dist(k) = Atmp_nu 
      end do 
      ! else if (radialace%type_chem_radial == ACE_CHEM_RADIAL_HSVD) then
      !   do kk = 1, base_params(1)%kmax
      !     k = k + 1
      !     ll = 0
      !     mm = 0  
      !     Atmp_nu = 0.0_dp 
      !     do ia_n = 1, neigh_ja%max_neigh_local
      !       type_ia_db = neigh_ja%i_type_db(ia_n)
      !       rr = neigh_ja%r_central(ia_n)
      !       dxp(:) = neigh_ja%tmp_dxp(:, ia_n)
      !       ispline = radialace%ic_rad(mu, type_ia_db, kk, ll)
      !       sharm_ll_mm = 1.0_dp/sqrt(4.0_dp*one_pi)
      !       rad = rsp(ia_n, ispline)
      !       Atmp_nu = Atmp_nu +  delta_ace(mu, type_ia_db) * rad *sharm_ll_mm
      !       drad = drsp(ia_n, ispline)
      !       iin = 3*(ia_n-1)
      !       d_mat_dist(iin+1:iin+3, k) = delta_ace(mu, type_ia_db)*(sharm_ll_mm * drad * dxp(1:3) / rr)
      !     end do    
      !     mat_dist(k) = Atmp_nu 
      !   end do 
      ! end if
    end do 

    baseB(:) = mat_dist(:)
    do ia_n = 1, neigh_ja%max_neigh_local
      do ix = 1,3 
        d_baseB(:, ix, ia_n) = d_mat_dist(ix + 3*(ia_n-1), :)
      end do 
    end do 

    ! do iin = 1, size(baseB)
    ! !  write(*,*) " n=", base_cnlm(nu)%dico_nb(idx_nb)%nbold(:), ", l=", base_cnlm(nu)%dico_lb0(idx_lb0)%lbold(:), ", baseB(iin)%re, baseB(iin)%im
    !   write(*,*) "baseB(ib) real / imaginary ", baseB(iin)%re, baseB(iin)%im
    ! end do
  end subroutine compute_AB_body_one_neigh_ja


end module module_AB_basis
