!===================================================================!
!     MiLaDy - Machine Learning Dynamics                            !
!                                                                   !
!     module_online_fit.F90                                         !
!     Online (streaming) Gram-matrix accumulation for fit_type=6    !
!                                                                   !
!     Accumulates G = A W A^T  and  h = A W y  one configuration   !
!     at a time, then solves the D x D normal equations via SVD.    !
!                                                                   !
!     Two modes:                                                    !
!       Serial   (scalapack_driver=.false.): dense G(D,D) per rank  !
!       Parallel (scalapack_driver=.true.):  block-cyclic sca_G     !
!         distributed over the BLACS grid — each rank holds only    !
!         its local tile (~D^2/p elements).                         !
!                                                                   !
!     See docs/README_ONLINE_FIT.md for full design document.       !
!===================================================================!

#include "../MLD_MACROS.INC"

module module_online_fit
   use module_kind_variables, only: kind_double
   use mld_logger
   implicit none
   private

   type, public :: online_gram_accumulator_t
      ! --- Serial mode (dense) ---
      real(kind_double), allocatable :: G(:,:)      ! D x D Gram matrix (dense, serial only)
      real(kind_double), allocatable :: h(:,:)      ! D x 1 RHS vector  (dense, serial only)
      ! --- ScaLAPACK mode (block-cyclic distributed) ---
      real(kind_double), allocatable :: sca_G(:,:)  ! local tile of D x D Gram matrix
      integer, allocatable           :: desc_G(:)   ! ScaLAPACK descriptor for sca_G
      integer :: l_dr_G = 0                         ! local row count of sca_G tile
      integer :: l_dc_G = 0                         ! local col count of sca_G tile
      integer :: nb_G = 0                           ! block size for sca_G
      real(kind_double), allocatable :: sca_h(:,:)  ! local tile of D x 1 RHS vector
      integer, allocatable           :: desc_h(:)   ! ScaLAPACK descriptor for sca_h
      integer :: l_dr_h = 0                         ! local row count of sca_h tile
      ! --- ScaLAPACK mode: per-config buffers (Strategy C) ---
      real(kind_double), allocatable :: Ac_config(:,:) ! D x mc, current config columns
      real(kind_double), allocatable :: yc_config(:,:) ! mc x 1, current config RHS
      integer :: mc_config = 0                         ! observations in current config
      ! --- Common ---
      integer :: D = 0                              ! = dim_design_line
      integer :: M_local = 0                        ! local observation count (cumulative)
      logical :: use_scalapack = .false.             ! mode flag
   contains
      procedure :: init       => gram_init
      procedure :: accumulate => gram_accumulate_config
      procedure :: scatter_config_to_blacs => gram_scatter_config_to_blacs
      procedure :: reduce_and_solve => gram_reduce_and_solve
      procedure :: cleanup    => gram_cleanup
   end type

   public :: main_train_online_fit

contains

   !----------------------------------------------------------------
   ! Initialise the accumulator.
   ! Serial mode: allocate dense G(D,D), h(D,1).
   ! ScaLAPACK mode: allocate block-cyclic sca_G, sca_h on BLACS grid.
   !----------------------------------------------------------------
   subroutine gram_init(self, D)
      use module_ml_scalapack, only: scalapack_driver, &
         myrow, mycol, nprow, npcol, context, iam
      use module_scalapack_tools, only: blockset_ml
      class(online_gram_accumulator_t), intent(inout) :: self
      integer, intent(in) :: D
      integer :: nb, info, lld, numroc
      real(kind_double) :: mem_MB

      self%D = D
      self%M_local = 0
      self%use_scalapack = scalapack_driver

      if (self%use_scalapack) then
         ! --- ScaLAPACK mode: distributed sca_G (D x D) ---
         nb = -1
         call blockset_ml(nb, -1, D, nprow, npcol)
         self%nb_G = nb

         self%l_dr_G = numroc(D, nb, myrow, 0, nprow)
         self%l_dc_G = numroc(D, nb, mycol, 0, npcol)

         if (allocated(self%desc_G)) deallocate(self%desc_G)
         allocate(self%desc_G(9))
         lld = max(1, self%l_dr_G)
         call descinit(self%desc_G, D, D, nb, nb, 0, 0, context, lld, info)

         if (allocated(self%sca_G)) deallocate(self%sca_G)
         allocate(self%sca_G(max(1, self%l_dr_G), max(1, self%l_dc_G)))
         self%sca_G = 0.0_kind_double

         ! --- Distributed sca_h (D x 1) ---
         self%l_dr_h = numroc(D, nb, myrow, 0, nprow)

         if (allocated(self%desc_h)) deallocate(self%desc_h)
         allocate(self%desc_h(9))
         lld = max(1, self%l_dr_h)
         call descinit(self%desc_h, D, 1, nb, 1, 0, 0, context, lld, info)

         if (allocated(self%sca_h)) deallocate(self%sca_h)
         allocate(self%sca_h(max(1, self%l_dr_h), 1))
         self%sca_h = 0.0_kind_double

         mem_MB = 8.0d0 * dble(self%l_dr_G) * dble(self%l_dc_G) / (1024.0d0*1024.0d0)
         call log_info("online_fit: ScaLAPACK Gram accumulator, D = " // vtoa(D) &
                       // ", block = " // vtoa(nb) &
                       // ", local tile = " // vtoa(self%l_dr_G) // " x " // vtoa(self%l_dc_G) &
                       // ", tile memory = " // vtoa(int(mem_MB)) // " MB")
      else
         ! --- Serial mode: dense G(D,D), h(D,1) ---
         if (allocated(self%G)) deallocate(self%G)
         if (allocated(self%h)) deallocate(self%h)
         allocate(self%G(D, D))
         allocate(self%h(D, 1))
         self%G = 0.0_kind_double
         self%h = 0.0_kind_double

         call log_info("online_fit: Dense Gram accumulator, D = " // vtoa(D) &
                       // ", G memory = " // vtoa(int(8.0d0*D*D / (1024.0d0*1024.0d0))) // " MB")
      end if
   end subroutine gram_init


   !----------------------------------------------------------------
   ! Accumulate one configuration into the Gram matrix.
   ! For each config: pack energy/force/stress descriptor columns
   ! into a small local buffer Ac(D, mc), apply weights, then
   ! rank-update G += Ac * diag(wc) * Ac^T  and  h += Ac * diag(wc) * yc
   !
   ! Serial mode:    DSYRK into dense G, DGEMV into dense h
   ! ScaLAPACK mode: tile-local accumulation into sca_G and sca_h
   !   — only the elements owned by this rank are computed and stored.
   !----------------------------------------------------------------
   subroutine gram_accumulate_config(self, iconf)
      use ml_in_ndm_module, only: desc_forces, &
         mld_order, mld_linear, mld_linear_extended, &
         mld_quadratic, mld_polyc, mld_kernel
      use derived_types, only: config_desc, config_real
      use snap, only: dim_design_line, dim_xdesc_linear
      use temporary_data_cov, only: dim_xdesc_patch
      use module_mld_quadratic, only: dim_xdesc_quadratic
      use module_mld_polyc, only: dim_xdesc_polyc
      use module_mld_kernel, only: dim_xdesc_kernel
      use module_kernel_2b, only: activate_k2b, dim_kernel_2b
      use module_db_poscar, only: iread_energy

      class(online_gram_accumulator_t), intent(inout) :: self
      integer, intent(in) :: iconf

      real(kind_double), allocatable :: Ac(:,:)   ! D x mc
      real(kind_double), allocatable :: yc(:)     ! mc
      real(kind_double), allocatable :: wc(:)     ! mc
      real(kind_double), allocatable :: tmp(:,:)  ! for buffer reallocation
      real(kind_double), allocatable :: tmp1(:,:) ! for buffer reallocation
      integer :: mc, D, col, ik, ix, new_cap

      _NAMECURRENT_("gram_accumulate_config")

      D = self%D

      ! --- Count observations for this config ---
      mc = 0
      if (config_real(iconf)%has_energy) mc = mc + 1
      if ((config_real(iconf)%has_force) .and. desc_forces) mc = mc + 3 * config_real(iconf)%nat
      if ((config_real(iconf)%has_stress) .and. desc_forces) mc = mc + 6

      if (mc == 0) return

      allocate(Ac(D, mc), yc(mc), wc(mc))
      Ac = 0.0_kind_double
      yc = 0.0_kind_double
      wc = 0.0_kind_double

      col = 0

      ! --- Pack energy column ---
      if (config_real(iconf)%has_energy) then
         col = col + 1
         call pack_energy_column(iconf, Ac(:, col))
         yc(col) = config_real(iconf)%energy(iread_energy)
         wc(col) = config_real(iconf)%w_e_model
      end if

      ! --- Pack force columns (3 per atom) ---
      if ((config_real(iconf)%has_force) .and. desc_forces) then
         do ik = 1, config_real(iconf)%nat
            do ix = 1, 3
               col = col + 1
               call pack_force_column(iconf, ik, ix, Ac(:, col))
               yc(col) = config_real(iconf)%force(ix, ik)
               wc(col) = config_real(iconf)%w_f_model
            end do
         end do
      end if

      ! --- Pack stress columns (6 Voigt) ---
      if ((config_real(iconf)%has_stress) .and. desc_forces) then
         do ix = 1, 6
            col = col + 1
            call pack_stress_column(iconf, ix, Ac(:, col))
            yc(col) = config_real(iconf)%stress(ix)
            wc(col) = config_real(iconf)%w_s_model
         end do
      end if

      ! --- Apply sqrt(weight) to Ac columns and yc ---
      ! G += Ac * diag(wc) * Ac^T  =  (Ac * sqrt(wc)) * (Ac * sqrt(wc))^T
      do ix = 1, col
         Ac(:, ix) = Ac(:, ix) * dsqrt(wc(ix))
         yc(ix) = yc(ix) * dsqrt(wc(ix))
      end do

      if (self%use_scalapack) then
         ! ============================================================
         ! ScaLAPACK mode (Strategy C, batched): append this config's
         ! weighted Ac/yc columns onto the running batch buffer.
         ! scatter_config_to_blacs is called by the main loop only
         ! every online_fit_batch_size configs (or at the last step),
         ! amortizing the collective broadcast + barrier cost over
         ! several configs instead of paying it once per config.
         ! Peak memory is O(D * batch worth of observations), bounded
         ! by online_fit_batch_size -- not O(D * all local configs).
         ! ============================================================
         if (self%mc_config == 0) then
            if (allocated(self%Ac_config)) deallocate(self%Ac_config)
            if (allocated(self%yc_config)) deallocate(self%yc_config)
            allocate(self%Ac_config(D, col))
            allocate(self%yc_config(col, 1))
            self%Ac_config(:, 1:col) = Ac(:, 1:col)
            self%yc_config(1:col, 1) = yc(1:col)
         else
            allocate(tmp(D, self%mc_config + col))
            tmp(:, 1:self%mc_config) = self%Ac_config(:, 1:self%mc_config)
            tmp(:, self%mc_config + 1 : self%mc_config + col) = Ac(:, 1:col)
            call move_alloc(tmp, self%Ac_config)

            allocate(tmp1(self%mc_config + col, 1))
            tmp1(1:self%mc_config, 1) = self%yc_config(1:self%mc_config, 1)
            tmp1(self%mc_config + 1 : self%mc_config + col, 1) = yc(1:col)
            call move_alloc(tmp1, self%yc_config)
         end if
         self%mc_config = self%mc_config + col

      else
         ! ============================================================
         ! Serial mode: standard BLAS accumulation into dense G and h
         ! ============================================================
         ! Rank-update: G += Ac * Ac^T  (DSYRK, upper triangle)
         call dsyrk('U', 'N', D, col, 1.0_kind_double, Ac, D, &
                    1.0_kind_double, self%G, D)

         ! RHS update: h += Ac * yc  (DGEMV)
         call dgemv('N', D, col, 1.0_kind_double, Ac, D, yc, 1, &
                    1.0_kind_double, self%h(:,1), 1)
      end if

      self%M_local = self%M_local + col

      deallocate(Ac, yc, wc)

   end subroutine gram_accumulate_config


   !----------------------------------------------------------------
   ! Pack the energy descriptor column for config iconf into vec(D)
   ! Mirrors the logic of train_fill_Amat_with_energy in mld_energy.F90
   !----------------------------------------------------------------
   subroutine pack_energy_column(iconf, vec)
      use ml_in_ndm_module, only: mld_order, &
         mld_linear, mld_linear_extended, &
         mld_quadratic, mld_polyc, mld_kernel
      use derived_types, only: config_desc, config_real
      use snap, only: dim_design_line, dim_xdesc_linear
      use temporary_data_cov, only: dim_xdesc_patch
      use module_mld_quadratic, only: dim_xdesc_quadratic
      use module_mld_polyc, only: dim_xdesc_polyc
      use module_mld_kernel, only: dim_xdesc_kernel
      use module_kernel_2b, only: activate_k2b, dim_kernel_2b

      integer, intent(in) :: iconf
      real(kind_double), intent(out) :: vec(:)

      vec = 0.0_kind_double
      vec(1) = dble(config_real(iconf)%nat)

      if ((mld_order == mld_linear) .or. (mld_order == mld_linear_extended)) then
         vec(2 : dim_xdesc_linear) = config_desc(iconf)%pack_energy_linear(1:dim_xdesc_linear - 1)
      end if

      if (mld_order == mld_quadratic) then
         vec(2 : dim_xdesc_linear) = config_desc(iconf)%pack_energy_linear(1:dim_xdesc_linear - 1)
         vec(dim_xdesc_linear + dim_xdesc_patch + 1 : dim_design_line) = &
            config_desc(iconf)%pack_energy_quadratic(dim_xdesc_linear : dim_xdesc_quadratic - 1)
      end if

      if (mld_order == mld_polyc) then
         vec(2 : dim_xdesc_linear) = config_desc(iconf)%pack_energy_linear(1:dim_xdesc_linear - 1)
         vec(dim_xdesc_linear + dim_xdesc_patch + 1 : dim_design_line) = &
            config_desc(iconf)%pack_energy_polyc(dim_xdesc_linear : dim_xdesc_polyc - 1)
      end if

      if (mld_order == mld_kernel) then
         vec(2 : dim_xdesc_linear) = config_desc(iconf)%pack_energy_kernel(1 : dim_xdesc_linear - 1)
         vec(dim_xdesc_linear + dim_xdesc_patch + 1 : dim_design_line) = &
            config_desc(iconf)%pack_energy_kernel(dim_xdesc_linear : dim_xdesc_kernel - 1)
      end if

      if (activate_k2b) then
         vec(dim_xdesc_linear + 1: dim_xdesc_linear + dim_kernel_2b) = &
            config_desc(iconf)%pack_energy_k2b(1: dim_kernel_2b)
      end if

   end subroutine pack_energy_column


   !----------------------------------------------------------------
   ! Pack one force descriptor column for config iconf, atom ik, direction ix
   ! Mirrors the logic of train_fill_Amat_with_force in mld_forces.F90
   !----------------------------------------------------------------
   subroutine pack_force_column(iconf, ik, ix, vec)
      use ml_in_ndm_module, only: mld_order, &
         mld_linear, mld_linear_extended, &
         mld_quadratic, mld_polyc, mld_kernel
      use derived_types, only: config_desc, config_real
      use snap, only: dim_design_line, dim_xdesc_linear
      use temporary_data_cov, only: dim_xdesc_patch
      use module_mld_quadratic, only: dim_xdesc_quadratic
      use module_mld_polyc, only: dim_xdesc_polyc
      use module_mld_kernel, only: dim_xdesc_kernel
      use module_kernel_2b, only: activate_k2b, dim_kernel_2b

      integer, intent(in) :: iconf, ik, ix
      real(kind_double), intent(out) :: vec(:)

      vec = 0.0_kind_double
      vec(1) = 0.0_kind_double

      if ((mld_order == mld_linear) .or. (mld_order == mld_linear_extended)) then
         vec(2 : dim_xdesc_linear) = config_desc(iconf)%pack_force_linear(1:dim_xdesc_linear-1, ix, ik)
      end if

      if (mld_order == mld_quadratic) then
         vec(2 : dim_xdesc_linear) = config_desc(iconf)%pack_force_linear(1:dim_xdesc_linear-1, ix, ik)
         vec(dim_xdesc_linear + dim_xdesc_patch + 1 : dim_design_line) = &
            config_desc(iconf)%pack_force_quadratic(dim_xdesc_linear:dim_xdesc_quadratic - 1, ix, ik)
      end if

      if (mld_order == mld_polyc) then
         vec(2 : dim_xdesc_linear) = config_desc(iconf)%pack_force_linear(1:dim_xdesc_linear-1, ix, ik)
         vec(dim_xdesc_linear + dim_xdesc_patch + 1 : dim_design_line) = &
            config_desc(iconf)%pack_force_polyc(dim_xdesc_linear : dim_xdesc_polyc - 1, ix, ik)
      end if

      if (mld_order == mld_kernel) then
         vec(2 : dim_xdesc_linear) = config_desc(iconf)%pack_force_linear(1:dim_xdesc_linear-1, ix, ik)
         vec(dim_xdesc_linear + dim_xdesc_patch + 1 : dim_design_line) = &
            config_desc(iconf)%pack_force_kernel(dim_xdesc_linear : dim_xdesc_kernel - 1, ix, ik)
      end if

      if (activate_k2b) then
         vec(dim_xdesc_linear + 1: dim_xdesc_linear + dim_kernel_2b) = &
            config_desc(iconf)%pack_force_k2b(1: dim_kernel_2b, ix, ik)
      end if

   end subroutine pack_force_column


   !----------------------------------------------------------------
   ! Pack one stress descriptor column for config iconf, Voigt component ix
   ! Mirrors the logic of train_fill_Amat_with_stress in mld_stress.F90
   !----------------------------------------------------------------
   subroutine pack_stress_column(iconf, ix, vec)
      use ml_in_ndm_module, only: mld_order, &
         mld_linear, mld_linear_extended, &
         mld_quadratic, mld_polyc, mld_kernel
      use derived_types, only: config_desc, config_real
      use snap, only: dim_design_line, dim_xdesc_linear
      use temporary_data_cov, only: dim_xdesc_patch
      use module_mld_quadratic, only: dim_xdesc_quadratic
      use module_mld_polyc, only: dim_xdesc_polyc
      use module_mld_kernel, only: dim_xdesc_kernel
      use module_kernel_2b, only: activate_k2b, dim_kernel_2b

      integer, intent(in) :: iconf, ix
      real(kind_double), intent(out) :: vec(:)

      vec = 0.0_kind_double
      vec(1) = 0.0_kind_double

      if ((mld_order == mld_linear) .or. (mld_order == mld_linear_extended)) then
         vec(2:dim_xdesc_linear) = config_desc(iconf)%pack_stress_linear(1:dim_xdesc_linear-1, ix)
      end if

      if (mld_order == mld_quadratic) then
         vec(2:dim_xdesc_linear) = config_desc(iconf)%pack_stress_linear(1:dim_xdesc_linear - 1, ix)
         vec(dim_xdesc_linear + dim_xdesc_patch + 1:dim_design_line) = &
            config_desc(iconf)%pack_stress_quadratic(dim_xdesc_linear:dim_xdesc_quadratic - 1, ix)
      end if

      if (mld_order == mld_polyc) then
         vec(2:dim_xdesc_linear) = config_desc(iconf)%pack_stress_linear(1:dim_xdesc_linear - 1, ix)
         vec(dim_xdesc_linear + dim_xdesc_patch + 1:dim_design_line) = &
            config_desc(iconf)%pack_stress_polyc(dim_xdesc_linear : dim_xdesc_polyc - 1, ix)
      end if

      if (mld_order == mld_kernel) then
         vec(2:dim_xdesc_linear) = config_desc(iconf)%pack_stress_linear(1:dim_xdesc_linear - 1, ix)
         vec(dim_xdesc_linear + dim_xdesc_patch + 1:dim_design_line) = &
            config_desc(iconf)%pack_stress_kernel(dim_xdesc_linear : dim_xdesc_kernel - 1, ix)
      end if

      if (activate_k2b) then
         vec(dim_xdesc_linear + 1: dim_xdesc_linear + dim_kernel_2b) = &
            config_desc(iconf)%pack_stress_k2b(1: dim_kernel_2b, ix)
      end if

   end subroutine pack_stress_column


   !----------------------------------------------------------------
   ! Reduce G and h across MPI subworlds, add regularization,
   ! then solve G * w = h via SVD.
   !
   ! Serial mode:  MPI_Reduce dense G(D,D) → serial SVD on rank 0
   ! ScaLAPACK mode: MPI_Reduce on local tiles (~D^2/p) →
   !   symmetrise → add lambda → ScaLAPACK SVD on sca_G
   !----------------------------------------------------------------
   subroutine gram_reduce_and_solve(self, w_params_out)
      use ml_in_ndm_module, only: rangml, lambda_krr, svd_rcond
      use module_ml_scalapack, only: scalapack_driver, &
         sca_phi, desc_sca_phi, dimr_sca_phi, dimc_sca_phi, &
         l_dimr_sca_phi, l_dimc_sca_phi, nbr_phi, nbc_phi, &
         sca_ymat_qr_svd, desc_sca_ymat_qr_svd, &
         dimr_sca_ymat_qr_svd, dimc_sca_ymat_qr_svd, &
         l_dimr_sca_ymat_qr_svd, l_dimc_sca_ymat_qr_svd, &
         nbr_ymat_qr_svd, nbc_ymat_qr_svd, &
         myrow, mycol, nprow, npcol, iam, context, &
         nprocs_ml_sca, iproc_sca
      use module_scalapack_tools, only: allocation_scalapack_matrix, blockset_ml
      use module_fit_ScaMatrix, only: scalapack_lsystem_by_SVD
      use module_serial_linear_solver, only: serial_lsystem_by_svd
      use my_mpi_subroutines, only: subworlds_allreduce_matrix_double, &
                                    subworlds_allreduce_vect_double
      use mld_mpi, only: mpi_comm_mld, comm_mld, mld_rank
#if(PARA)
      use mpi
#endif

      class(online_gram_accumulator_t), intent(inout) :: self
      real(kind_double), intent(out) :: w_params_out(:,:)

      integer :: D, i, rank_sca, ierr
      integer :: il, jl, nb
      real(kind_double) :: svd_rcond_local
      real(kind_double), allocatable :: phi_serial(:,:), rhs_serial(:,:)
      integer :: INDXG2L, INDXG2P

      _NAMECURRENT_("gram_reduce_and_solve")

      _MLD_BEGIN_

      D = self%D

      call log_info("online_fit: Gram accumulation completed, M_local = " // vtoa(self%M_local))

      if (self%use_scalapack) then
         ! =============================================================
         ! ScaLAPACK mode: sca_G and sca_h already accumulated via
         ! per-config scatter + PDGEMM in gram_scatter_config_to_blacs.
         ! Here we only add regularization and solve by SVD.
         ! =============================================================
         nb = self%nb_G

         call blacs_barrier(context, 'A')

         ! --- Add Tikhonov regularization to diagonal ---
         if (lambda_krr < 0.0_kind_double) lambda_krr = 0.0_kind_double
         do i = 1, D
            if (INDXG2P(i, nb, myrow, 0, nprow) == myrow .and. &
                INDXG2P(i, nb, mycol, 0, npcol) == mycol) then
               il = INDXG2L(i, nb, myrow, 0, nprow)
               jl = INDXG2L(i, nb, mycol, 0, npcol)
               self%sca_G(il, jl) = self%sca_G(il, jl) + lambda_krr
            end if
         end do

         ! --- Solve via ScaLAPACK SVD ---
         call log_info("online_fit: solving D x D system via ScaLAPACK SVD, D = " // vtoa(D))

         ! Copy sca_G -> sca_phi (SVD is destructive)
         dimr_sca_phi = D
         dimc_sca_phi = D
         nbr_phi = nb
         nbc_phi = nb
         l_dimr_sca_phi = self%l_dr_G
         l_dimc_sca_phi = self%l_dc_G

         if (allocated(sca_phi)) deallocate(sca_phi)
         allocate(sca_phi(max(1, l_dimr_sca_phi), max(1, l_dimc_sca_phi)))
         sca_phi(:,:) = self%sca_G(:,:)

         if (allocated(desc_sca_phi)) deallocate(desc_sca_phi)
         allocate(desc_sca_phi(9))
         desc_sca_phi(:) = self%desc_G(:)

         ! Copy sca_h -> sca_ymat_qr_svd
         dimr_sca_ymat_qr_svd = D
         dimc_sca_ymat_qr_svd = 1
         nbr_ymat_qr_svd = nb
         nbc_ymat_qr_svd = 1
         l_dimr_sca_ymat_qr_svd = self%l_dr_h
         l_dimc_sca_ymat_qr_svd = 1

         if (allocated(sca_ymat_qr_svd)) deallocate(sca_ymat_qr_svd)
         allocate(sca_ymat_qr_svd(max(1, l_dimr_sca_ymat_qr_svd), 1))
         sca_ymat_qr_svd(:,1) = self%sca_h(:,1)

         if (allocated(desc_sca_ymat_qr_svd)) deallocate(desc_sca_ymat_qr_svd)
         allocate(desc_sca_ymat_qr_svd(9))
         desc_sca_ymat_qr_svd(:) = self%desc_h(:)

         call blacs_barrier(context, 'A')

         call scalapack_lsystem_by_SVD(dimr_sca_phi, dimc_sca_phi, sca_phi, desc_sca_phi, &
            nbr_phi, nbc_phi, &
            dimr_sca_ymat_qr_svd, dimc_sca_ymat_qr_svd, &
            sca_ymat_qr_svd, desc_sca_ymat_qr_svd, w_params_out, &
            svd_rcond, rank_sca, context)

         call blacs_barrier(context, 'A')

         call log_info("online_fit: ScaLAPACK SVD solve done, rank = " // vtoa(rank_sca))

         ! Free SVD work arrays
         if (allocated(sca_phi)) deallocate(sca_phi)
         if (allocated(desc_sca_phi)) deallocate(desc_sca_phi)
         if (allocated(sca_ymat_qr_svd)) deallocate(sca_ymat_qr_svd)
         if (allocated(desc_sca_ymat_qr_svd)) deallocate(desc_sca_ymat_qr_svd)

      else
         ! =============================================================
         ! Serial mode: dense G, MPI reduce, serial SVD
         ! =============================================================
         !
         ! NOTE: pack_force_descriptor already MPI_ALLREDUCEs force
         ! descriptors within the subworld, so all ranks in the subworld
         ! hold identical G and h after gram_accumulate_config.
         ! Therefore we must NOT MPI_Reduce within the subworld —
         ! that would multiply G and h by procs_per_file.
         ! We only need the cross-subworld allreduce (masters_world)
         ! to sum contributions from different config ranges, then
         ! broadcast the result back to all subworld members.
         ! =============================================================

         ! --- Symmetrise G (DSYRK only fills upper triangle) ---
         do i = 1, D
            self%G(i+1:D, i) = self%G(i, i+1:D)
         end do

#if(PARA)
         ! --- Cross-subworld allreduce, then broadcast ---
         call subworlds_allreduce_matrix_double('G_online', self%G)
         call subworlds_allreduce_vect_double(self%h(:,1))
#endif

         ! --- Step 3: Add Tikhonov regularization ---
         if (lambda_krr < 0.0_kind_double) lambda_krr = 0.0_kind_double
         do i = 1, D
            self%G(i,i) = self%G(i,i) + lambda_krr
         end do

         ! --- Step 4: Serial SVD solve on rank 0, then broadcast ---
         call log_info("online_fit: solving D x D system via serial SVD, D = " // vtoa(D))

         if (rangml == 0) then
            allocate(phi_serial(D, D), rhs_serial(D, 1))
            phi_serial(:,:) = self%G(:,:)
            rhs_serial(:,:) = self%h(:,:)
            svd_rcond_local = svd_rcond
            call serial_lsystem_by_svd(phi_serial, rhs_serial, D, D, w_params_out, svd_rcond_local)
            deallocate(phi_serial, rhs_serial)
         end if

#if(PARA)
         call comm_mld%barrier()
         call comm_mld%bcast(0, w_params_out(:,1))
#endif
      end if

      call log_info("online_fit: solve completed, w_params norm = " // &
         vtoa(dsqrt(dot_product(w_params_out(:,1), w_params_out(:,1)))))

      _MLD_END_
   end subroutine gram_reduce_and_solve


   !----------------------------------------------------------------
   ! Deallocate accumulator arrays
   !----------------------------------------------------------------
   subroutine gram_cleanup(self)
      class(online_gram_accumulator_t), intent(inout) :: self

      if (allocated(self%G)) deallocate(self%G)
      if (allocated(self%h)) deallocate(self%h)
      if (allocated(self%sca_G)) deallocate(self%sca_G)
      if (allocated(self%sca_h)) deallocate(self%sca_h)
      if (allocated(self%desc_G)) deallocate(self%desc_G)
      if (allocated(self%desc_h)) deallocate(self%desc_h)
      if (allocated(self%Ac_config)) deallocate(self%Ac_config)
      if (allocated(self%yc_config)) deallocate(self%yc_config)
      self%D = 0
      self%M_local = 0
      self%mc_config = 0
      self%l_dr_G = 0
      self%l_dc_G = 0
      self%nb_G = 0
      self%l_dr_h = 0
      self%use_scalapack = .false.
   end subroutine gram_cleanup


   !----------------------------------------------------------------
   ! Strategy C: scatter one config's Ac/yc from each subworld
   ! master to the BLACS grid and accumulate via PDGEMM.
   !
   ! All P ranks must call this (BLACS operations are collective).
   ! For each subworld with mc_config > 0, we:
   !   1. pdgemr2d Ac_config(D, mc) from subworld master → sca_Ac
   !   2. pdgemr2d yc_config(mc, 1) from subworld master → sca_yc
   !   3. PDGEMM: sca_G += sca_Ac * sca_Ac^T
   !   4. PDGEMM: sca_h += sca_Ac * sca_yc
   !   5. Deallocate sca_Ac, sca_yc
   !
   ! After this call, Ac_config/yc_config are freed and mc_config=0.
   !----------------------------------------------------------------
   subroutine gram_scatter_config_to_blacs(self)
      use module_ml_scalapack, only: scalapack_driver, &
         myrow, mycol, nprow, npcol, iam, context, &
         nprocs_ml_sca, iproc_sca
      use mld_subworld, only: subworld, subrank, masters_world, nb_subworlds, id_subworld
      use mld_mpi, only: mpi_comm_mld, comm_mld, mld_rank
      use module_db_poscar, only: procs_per_file
      use time_check_general, only: MY_MPI_WTIME, &
         time_scatter_gather, time_scatter_bcast, time_scatter_dgemm, time_scatter_barrier
#if(PARA)
      use mpi
#endif

      class(online_gram_accumulator_t), intent(inout) :: self

      integer :: D, nb, gg, n_groups, root_proc, info, ierr
      integer :: mc_g
      integer, allocatable :: mc_per_group(:)
      real(kind_double), parameter :: one = 1.0_kind_double, zero = 0.0_kind_double
      ! --- For broadcast + local DGEMM path ---
      real(kind_double), allocatable :: Ac_full(:,:), yc_full(:,:)
      real(kind_double), allocatable :: Ac_rows(:,:), Ac_cols(:,:)
      integer :: il, jl, gi, gj, INDXL2G
      real(kind_double) :: t_sc_a, t_sc_b

      _NAMECURRENT_("gram_scatter_config_to_blacs")

      D = self%D
      nb = self%nb_G
      n_groups = nprocs_ml_sca / procs_per_file

      t_sc_a = MY_MPI_WTIME()

      ! --- Gather mc from each subworld ---
      allocate(mc_per_group(0:n_groups-1))
      mc_per_group(:) = 0
      if (subrank == 0) mc_per_group(id_subworld) = self%mc_config
#if(PARA)
      if (MPI_COMM_NULL /= masters_world) then
         call MPI_Allreduce(MPI_IN_PLACE, mc_per_group, n_groups, &
            MPI_INTEGER, MPI_SUM, masters_world, ierr)
      end if
      call MPI_Bcast(mc_per_group, n_groups, MPI_INTEGER, 0, subworld, ierr)
#endif
      t_sc_b = MY_MPI_WTIME()
      time_scatter_gather = time_scatter_gather + (t_sc_b - t_sc_a)

      ! --- Loop over subworlds: scatter each config's Ac, PDGEMM ---
      do gg = 0, n_groups - 1
         root_proc = gg * procs_per_file
         mc_g = mc_per_group(gg)
         if (mc_g == 0) cycle

         ! ============================================================
         ! Always use broadcast + local DGEMM to avoid MKL PB_CVMnpq
         ! integer divide-by-zero bug.  MKL's PDGEMM crashes when
         ! block sizes are large relative to the K dimension (mc_g),
         ! which happens for per-config matrices.  Since mc_g is at
         ! most 3*nat+7 (one config), broadcasting D*mc_g doubles is
         ! cheap, and the local DGEMM on each rank's tile of sca_G
         ! is correct and efficient.
         ! ============================================================

         ! Broadcast Ac_config (D x mc_g) and yc_config (mc_g x 1) from root
         allocate(Ac_full(D, mc_g), yc_full(mc_g, 1))
         if (mld_rank == root_proc) then
            Ac_full(:,:) = self%Ac_config(:, 1:mc_g)
            yc_full(:,:) = self%yc_config(1:mc_g, :)
         else
            Ac_full = zero
            yc_full = zero
         end if
         t_sc_a = MY_MPI_WTIME()
#if(PARA)
         call MPI_Bcast(Ac_full, D * mc_g, MPI_DOUBLE_PRECISION, &
                        root_proc, mpi_comm_mld, ierr)
         call MPI_Bcast(yc_full, mc_g, MPI_DOUBLE_PRECISION, &
                        root_proc, mpi_comm_mld, ierr)
#endif
         t_sc_b = MY_MPI_WTIME()
         time_scatter_bcast = time_scatter_bcast + (t_sc_b - t_sc_a)

         ! Each rank extracts the rows of Ac it owns (for sca_G rows)
         ! and the columns it owns (for sca_G columns), then does local DGEMM.
         !
         ! sca_G(il, jl) += sum_k  Ac_full(gi, k) * Ac_full(gj, k)
         !   where gi = INDXL2G(il, nb, myrow, 0, nprow)
         !         gj = INDXL2G(jl, nb, mycol, 0, npcol)
         !
         ! This is equivalent to: local_tile += Ac_rows * Ac_cols^T
         !   Ac_rows(l_dr_G, mc_g): rows of Ac_full owned by myrow
         !   Ac_cols(l_dc_G, mc_g): rows of Ac_full owned by mycol (= cols of G)

         if (self%l_dr_G > 0 .and. mc_g > 0) then
            allocate(Ac_rows(self%l_dr_G, mc_g))
            do il = 1, self%l_dr_G
               gi = INDXL2G(il, nb, myrow, 0, nprow)
               Ac_rows(il, :) = Ac_full(gi, :)
            end do
         end if

         if (self%l_dc_G > 0 .and. mc_g > 0) then
            allocate(Ac_cols(self%l_dc_G, mc_g))
            do jl = 1, self%l_dc_G
               gj = INDXL2G(jl, nb, mycol, 0, npcol)
               Ac_cols(jl, :) = Ac_full(gj, :)
            end do
         end if

         ! sca_G(local) += Ac_rows * Ac_cols^T
         if (self%l_dr_G > 0 .and. self%l_dc_G > 0 .and. mc_g > 0) then
            call dgemm('N', 'T', self%l_dr_G, self%l_dc_G, mc_g, one, &
                       Ac_rows, self%l_dr_G, Ac_cols, self%l_dc_G, &
                       one, self%sca_G, size(self%sca_G, 1))
         end if

         ! sca_h(local) += Ac_rows * yc_full
         if (self%l_dr_h > 0 .and. mc_g > 0) then
            call dgemm('N', 'N', self%l_dr_h, 1, mc_g, one, &
                       Ac_rows, self%l_dr_G, yc_full, mc_g, &
                       one, self%sca_h, size(self%sca_h, 1))
         end if

         if (allocated(Ac_rows)) deallocate(Ac_rows)
         if (allocated(Ac_cols)) deallocate(Ac_cols)
         deallocate(Ac_full, yc_full)
         t_sc_a = MY_MPI_WTIME()
         time_scatter_dgemm = time_scatter_dgemm + (t_sc_a - t_sc_b)

         ! No barrier here: MPI_Bcast is blocking and, per the MPI standard,
         ! it is safe to reuse/deallocate the buffers as soon as a rank's
         ! own call returns. Ac_full/yc_full/Ac_rows/Ac_cols are freshly
         ! (re)allocated every gg iteration, so there is no aliasing across
         ! iterations that would require an explicit synchronization point.
         ! With n_groups in the thousands and one lockstep step per
         ! training config, this fired ~(n_groups * n_configs) times on a
         ! production run (e.g. ~179,000 times at n_groups=1024,
         ! n_configs~175) -- a likely major, previously unmeasured
         ! contributor to the ~53000s "TRAIN TIME" of a 102-element,
         ! 8192-rank run where the rest of the scatter (bcast+dgemm+gather,
         ! now separately timed above) was orders of magnitude cheaper at
         ! small scale. time_scatter_barrier is kept (now always 0) so its
         ! absence is visible in the timing report.
      end do

      deallocate(mc_per_group)

      ! --- Free per-config buffers ---
      if (allocated(self%Ac_config)) deallocate(self%Ac_config)
      if (allocated(self%yc_config)) deallocate(self%yc_config)
      self%mc_config = 0

   end subroutine gram_scatter_config_to_blacs


   !----------------------------------------------------------------
   ! main_train_online_fit: The entry point for mld_fit_type=6.
   ! Called from main_train_mld() in mld.F90 before Amat allocation.
   ! Orchestrates: init -> config loop -> solve -> predict errors
   !----------------------------------------------------------------
   subroutine main_train_online_fit()
      use ml_in_ndm_module, only: rangml, debug, mld_fit_type, &
         tmp_val_desc_max, regularization_name, train_time, desc_forces, online_fit_batch_size
      use module_db_setup, only: iconf_data
      use derived_types, only: config_real
      use module_optimization, only: optimize_weights_db, optimize_weights_chem
      use module_mld_quadratic, only: i_e_fit_snap, i_f_fit_snap, i_s_fit_snap
      use snap, only: i_fit_snap, i_constraints_snap, w_params, dim_design_line
      use module_kind_variables, only: kind_double
      use time_check_general, only: time_fit_params, time_full_desc, time_calc_desc, &
         time_neigh_desc, time_fill_desc, time_fill_desc_02, MY_MPI_WTIME, train_tot_time, tot_time, &
         debug_time, time_read_db, time_fill_desc_02A, time_fill_desc_02obj, time_fill_desc_02y, &
         time_scatter_blacs, time_scatter_gather, time_scatter_bcast, time_scatter_dgemm, time_scatter_barrier
      use module_kernel, only: time_for_desc_kernel, train_time_for_desc_kernel
      use module_kernel_2b, only: activate_k2b, train_time_for_desc_kernel_2b, time_for_desc_kernel_2b, &
         tnn_2b, tuu_2b, tkk_2b, r_cut_2b, r_cut_width_2b
      use module_continuity_k2b_zbl, only: get_continuity_k2b_zbl
      use module_zbl, only: zbl_potential, zbl_type, zbl_mode_default_k2b, r1_zbl, rr_k2b
      use module_neigh_local, only: type_fcut
      use module_write_parameters, only: write_zbl_k2b, write_snap_parameters
      use module_chemical_species, only: tnn_rdist, tcc01_rdist, tcc02_rdist, tcc03_rdist, tii_rdist
      use module_ftnbody, only: tnn_ftbd, t2b_ftbd, t3b_ftbd, t4b_ftbd, t5b_ftbd, &
         tnn_ftbd_train, t2b_ftbd_train, t3b_ftbd_train, t4b_ftbd_train, t5b_ftbd_train
      use module_optimization, only: itopt
      use module_db_poscar, only: i_start_conf, i_final_conf, procs_per_file
      use set_limits, only: set_limit_for_configs_with_MPI_grid
      use mld_subworld
      use my_mpi_subroutines, only: subworlds_get_max_double
      use module_ml_scalapack, only: scalapack_driver
      use module_online_predict, only: online_predict_train_errors
      use mld_energy_mod, only: pack_energy_descriptor
      use mld_force_mod, only: pack_force_descriptor
      use mld_stress_mod, only: pack_stress_descriptor
      use temporary_data_cov, only: dim_xdesc
      use main_mld_mod, only: set_mld_dimension
      use mld_mpi, only: mpi_comm_mld
#if(PARA)
      use mpi
#endif

      implicit none

      type(online_gram_accumulator_t) :: gram
      integer :: i
      real(kind_double) :: tmp_val, val_max
      real(kind_double) :: t0, t1, t2, t3, t4, t5, t6, t7
      logical :: post_desc
      integer :: n_train_local, max_steps, step, k, ierr_mpi, batch_size
      integer, allocatable :: train_list(:)

      _NAMECURRENT_("main_train_online_fit")

      _MLD_BEGIN_

      train_time = .true.

      call log_info("online_fit: === Starting online Gram-matrix fit (mld_fit_type = 6) ===")

      ! --- Phase 0: Prepare database and dimensions ---
      call prepare_database()
      call init_subworld(procs_per_file)
      call set_limit_for_configs_with_MPI_grid(iconf_data, i_start_conf, i_final_conf)
      call prepare_train_dimensions()

      ! We need set_mld_dimension to compute dim_design_line, but NOT allocate full Amat
      call set_mld_dimension()

      ! Allocate only w_params (small D x 1), NOT Amat
      if (allocated(w_params)) deallocate(w_params)
      allocate(w_params(dim_design_line, 1))
      w_params = 0.0_kind_double

      call log_info("online_fit: dim_design_line (D) = " // vtoa(dim_design_line))
      call log_info("online_fit: config range [" // vtoa(i_start_conf) // ", " // vtoa(i_final_conf) // "]")

      ! --- Phase 1: Gram matrix accumulation ---
      call gram%init(dim_design_line)

      i_fit_snap = 0
      i_e_fit_snap = 0
      i_f_fit_snap = 0
      i_s_fit_snap = 0
      i_constraints_snap = 0

      tmp_val = -1.0_kind_double

      time_fit_params = 0.0_kind_double
      time_full_desc = 0.0_kind_double
      time_scatter_blacs = 0.0_kind_double
      time_scatter_gather = 0.0_kind_double
      time_scatter_bcast = 0.0_kind_double
      time_scatter_dgemm = 0.0_kind_double
      time_scatter_barrier = 0.0_kind_double
      time_calc_desc = 0.0_kind_double
      time_neigh_desc = 0.0_kind_double
      time_fill_desc = 0.0_kind_double
      time_fill_desc_02 = 0.0_kind_double
      time_fill_desc_02y = 0.0_kind_double
      time_fill_desc_02obj = 0.0_kind_double
      time_fill_desc_02A = 0.0_kind_double

      if (itopt == 0) time_read_db = 0.0_kind_double

      t0 = MY_MPI_WTIME()

      if (gram%use_scalapack) then
         ! =============================================================
         ! ScaLAPACK mode: lockstep config loop.
         ! Each step: subworlds compute descriptors in parallel, then
         ! all P ranks collectively scatter + PDGEMM each subworld's
         ! per-config Ac into sca_G/sca_h.  Eliminates Ac_local.
         ! =============================================================

         ! Build local list of train config indices
         n_train_local = 0
         do i = i_start_conf, i_final_conf
            if (config_real(i)%train) n_train_local = n_train_local + 1
         end do
         allocate(train_list(max(1, n_train_local)))
         k = 0
         do i = i_start_conf, i_final_conf
            if (config_real(i)%train) then
               k = k + 1
               train_list(k) = i
            end if
         end do

         ! All subworlds agree on the number of lockstep iterations
         max_steps = n_train_local
#if(PARA)
         call MPI_Allreduce(MPI_IN_PLACE, max_steps, 1, MPI_INTEGER, &
            MPI_MAX, mpi_comm_mld, ierr_mpi)
#endif

         batch_size = max(1, online_fit_batch_size)

         call log_info("online_fit: ScaLAPACK lockstep loop, n_train_local = " // &
            vtoa(n_train_local) // ", max_steps = " // vtoa(max_steps) // &
            ", batch_size = " // vtoa(batch_size))

         do step = 1, max_steps
            if (step <= n_train_local) then
               i = train_list(step)

               t1 = MY_MPI_WTIME()
               call test_if_config_is_small(i)
               call calc_neighbours(i)
               t2 = MY_MPI_WTIME()
               time_neigh_desc = time_neigh_desc + (t2 - t1)

               post_desc = .true.
               call compute_descriptors(i, post_desc)
               t3 = MY_MPI_WTIME()
               time_calc_desc = time_calc_desc + (t3 - t2)

               if (config_real(i)%has_energy) call pack_energy_descriptor(i)
               if (config_real(i)%has_force .and. desc_forces) call pack_force_descriptor(i)
               if (config_real(i)%has_stress .and. desc_forces) call pack_stress_descriptor(i)
               ! Appends this config's columns onto the running batch buffer
               ! (gram%mc_config accumulates across up to batch_size steps);
               ! a subworld with no more configs this step simply adds
               ! nothing and keeps whatever it has already accumulated.
               call gram%accumulate(i)

               if (tmp_val_desc_max >= tmp_val) tmp_val = tmp_val_desc_max

               t4 = MY_MPI_WTIME()
               time_fill_desc = time_fill_desc + (t4 - t3)
            end if

            ! Collective: scatter the accumulated batch to BLACS, PDGEMM.
            ! Only every batch_size steps (or on the final step, to flush
            ! any partial batch) -- amortizes the per-call broadcast and
            ! synchronization cost over several configs instead of one.
            if (mod(step, batch_size) == 0 .or. step == max_steps) then
               t1 = MY_MPI_WTIME()
               call gram%scatter_config_to_blacs()
               t2 = MY_MPI_WTIME()
               time_scatter_blacs = time_scatter_blacs + (t2 - t1)
            end if

            if (step <= n_train_local) then
               call train_deallocate_desc(train_list(step))
            end if
         end do

         deallocate(train_list)

      else
         ! =============================================================
         ! Serial mode: independent async config loop with DSYRK
         ! =============================================================
         do i = i_start_conf, i_final_conf
            if (.not. (config_real(i)%train)) cycle

            t1 = MY_MPI_WTIME()
            call test_if_config_is_small(i)
            call calc_neighbours(i)
            t2 = MY_MPI_WTIME()
            time_neigh_desc = time_neigh_desc + (t2 - t1)

            post_desc = .true.
            call compute_descriptors(i, post_desc)
            t3 = MY_MPI_WTIME()
            time_calc_desc = time_calc_desc + (t3 - t2)

            if (config_real(i)%has_energy) call pack_energy_descriptor(i)
            if (config_real(i)%has_force .and. desc_forces) call pack_force_descriptor(i)
            if (config_real(i)%has_stress .and. desc_forces) call pack_stress_descriptor(i)
            call gram%accumulate(i)

            if (tmp_val_desc_max >= tmp_val) then
               tmp_val = tmp_val_desc_max
            end if

            t4 = MY_MPI_WTIME()
            time_fill_desc = time_fill_desc + (t4 - t3)

            call train_deallocate_desc(i)
         end do
      end if

      call subworlds_get_max_double(tmp_val, val_max)

      t5 = MY_MPI_WTIME()
      time_full_desc = time_full_desc + (t5 - t0)

      if ((.not. (optimize_weights_db .or. optimize_weights_chem)) .or. (itopt == 0)) then
         if (rangml == 0) then
            write(6, *) 'ML: the max value of the descriptor is   :', val_max
         end if
      end if

      call log_info("online_fit: descriptor accumulation done")

      ! --- Phase 2: Solve G * w = h ---
      t6 = MY_MPI_WTIME()

      call gram%reduce_and_solve(w_params)

      t7 = MY_MPI_WTIME()
      time_fit_params = time_fit_params + (t7 - t6)

      call gram%cleanup()

      ! --- Phase 3: Compute train errors by re-reading configs ---
      call log_info("online_fit: === Computing train errors (re-reading descriptors) ===")
      regularization_name = ""
      call online_predict_train_errors(w_params)

      ! Close the training subworld.
      ! Test prediction (main_test_mld) is called from ml.F90 and manages
      ! its own subworld lifecycle.
      call close_subworld()

      ! --- Bookkeeping (same as end of main_train_mld) ---
      if (activate_k2b) then
         if (zbl_potential .and. zbl_type == zbl_mode_default_k2b) then
            call get_continuity_k2b_zbl(r1_zbl, rr_k2b, type_fcut, r_cut_2b, r_cut_width_2b)
            call write_zbl_k2b
         end if
         train_time_for_desc_kernel_2b = time_for_desc_kernel_2b
         time_for_desc_kernel_2b = 0.0_kind_double
         tnn_2b = 0.0_kind_double
         tuu_2b = 0.0_kind_double
         tkk_2b = 0.0_kind_double
      end if

      tnn_rdist = 0.0_kind_double
      tcc01_rdist = 0.0_kind_double
      tcc02_rdist = 0.0_kind_double
      tcc03_rdist = 0.0_kind_double
      tii_rdist = 0.0_kind_double

      tnn_ftbd_train = tnn_ftbd
      tnn_ftbd = 0.0_kind_double
      t2b_ftbd_train = t2b_ftbd
      t2b_ftbd = 0.0_kind_double
      t3b_ftbd_train = t3b_ftbd
      t3b_ftbd = 0.0_kind_double
      t4b_ftbd_train = t4b_ftbd
      t4b_ftbd = 0.0_kind_double
      t5b_ftbd_train = t5b_ftbd
      t5b_ftbd = 0.0_kind_double

      if (debug_time) then
         train_tot_time(:) = tot_time(:)
         tot_time = 0.0_kind_double
      end if

      call write_snap_parameters

      call log_info("online_fit: === Online fit completed ===")

      _MLD_END_
      return
   end subroutine main_train_online_fit

end module module_online_fit
