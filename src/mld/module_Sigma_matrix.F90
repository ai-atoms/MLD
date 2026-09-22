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

module module_Sigma_matrix

   use module_kind_variables, only: kind_double
   logical  :: first_passage_mcd, analysis_local_energy
   integer  :: dim_col_design_mcd, dim_col_design_full
   ! The design matrix in the form .... dim_xdesc X dim_col_design_***
   !real(kind_double), dimension(:, :), allocatable      :: Phi_mcd, Phi_full
   real(kind_double), dimension(:, :), allocatable      :: eigenvectors_Sigma_sample_mcd, Sigma_sample_mcd, Sigma_sample_full, eigenvectors_Sigma_sample_full
   real(kind_double), dimension(:, :), allocatable      :: Sigma_sample_mcd_inv, Sigma_sample_full_inv
   real(kind_double), dimension(:), allocatable   :: mean_mcd, mean_full, eigenvalues_Sigma_sample_mcd, eigenvalues_Sigma_sample_full



   integer  :: dim_col_design, dim_line_design
   real(kind_double), dimension(:, :), allocatable      :: Phi_mat
   real(kind_double), dimension(:, :), allocatable      :: Sigma_sample, Sigma_sample_inv
   real(kind_double), dimension(:), allocatable   :: mean

   real(kind_double), dimension(:), allocatable   :: right_vector_with_energy


   integer :: method_sigma
   integer, parameter ::  method_sigma_by_Phi=1, method_sigma_by_sum=2


contains

   subroutine compute_Sigma_class_by_sum(classes_for_sigma, no_of_classes_for_sigma, dimr_Phi)
#if(PARA)
      !TORC! use mpi
      use mld_mpi, only: mld_mpi_abort, comm_mld
#endif
      use ml_in_ndm_module, only: rangml, i_start_at, i_final_at, lambda_krr
      use module_db_setup, only: iconf_data
      use derived_types, only: config_desc, config_real
      use set_limits, only: set_limit_for_atoms
      !use module_Sigma_matrix, only : mean, Sigma_sample, Sigma_sample_inv
      use math, only: dreal_inverse
      use mld_logger
      !TORC! use my_mpi_subroutines, only: my_barrier_mld

      implicit none

      integer, intent(in)  :: no_of_classes_for_sigma, dimr_Phi
      character(len=2), dimension(no_of_classes_for_sigma), intent(in)     :: classes_for_sigma
      integer  :: d1, d2, ia, iclass, ic, id, ndim

      _NAMECURRENT_("compute_Sigma_class_by_sum")

      _MLD_BEGIN_

      ! Sigma is DXD
      if (allocated(Sigma_sample)) deallocate (Sigma_sample)
      allocate (Sigma_sample(dimr_Phi, dimr_Phi))

      call log_info("ML: allocate Sigma, the  sample covariance matrix......:   " //&
         vtoa(size(Sigma_sample, 1))//vtoa(size(Sigma_sample, 2)))

      if (size(Sigma_sample, 1) < 2) then
         call mld_mpi_abort("ML: ......... no covariance matrix for one/zero sample. Increase number of samples.")
      end if


      Sigma_sample(:, :) = 0.d0
      do iclass = 1, size(classes_for_sigma)
         do ic = 1, iconf_data
            if ((classes_for_sigma(iclass) == config_real(ic)%class)) then
               if (.not. (allocated(config_desc(ic)%energy))) cycle
               call set_limit_for_atoms(rangml, config_real(ic)%nat, i_start_at, i_final_at)

               if ((i_start_at == 0) .and. (i_final_at == 0)) then
                  if (.not. (allocated(config_desc(ic)%energy))) return
               else
                  if (.not. (allocated(config_desc(ic)%energy))) return
                  do ia = i_start_at, i_final_at
                     do d1 = 1, dimr_Phi
                        do d2 = 1, dimr_Phi
                           Sigma_sample(d1, d2) = Sigma_sample(d1, d2) + (config_desc(ic)%energy(d1, ia) - mean(d1) ) * &
                              (config_desc(ic)%energy(d2, ia) - mean(d2) )
                        end do
                     end do
                  end do
               end if
            end if ! if class selection.

         end do
      end do

#if(PARA)
      ndim = size(Sigma_sample, 1) * size(Sigma_sample, 2)
      !TORC! call MPI_ALLREDUCE(MPI_IN_PLACE, Sigma_sample, ndim, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, codeml)
      call comm_mld%sum(Sigma_sample)
      !TORC! call my_barrier_mld(codeml)
      call comm_mld%barrier()
#endif

      Sigma_sample(:, :) = Sigma_sample(:, :)/dble(size(Sigma_sample, 1) - 1)

      lambda_krr = 1.d-9

      if (allocated(Sigma_sample_inv)) deallocate (Sigma_sample_inv); allocate (Sigma_sample_inv(size(Sigma_sample, 1), size(Sigma_sample, 2)))
      if (lambda_krr >= 0.d0) then
         do id = 1, size(Sigma_sample, 1)
            Sigma_sample(id, id) = Sigma_sample(id, id) + lambda_krr
         end do
      end if

      call dreal_inverse(size(Sigma_sample, 1), Sigma_sample, Sigma_sample_inv)


      _MLD_END_
   end subroutine compute_Sigma_class_by_sum



end module module_Sigma_matrix


module module_Phi
   use module_kind_variables, only: kind_double
   use ml_in_ndm_module, only: rangml
   use mld_logger

   implicit none
   integer :: dimr_Phi, dimc_Phi

contains

   subroutine get_dimension_design_Phi_by_class(classes_for_sigma, mean, dim_line_design, dim_col_design)
      !------------------------------------------------------
      ! THis subroutine fix by the selection given by class:
      !  -- mean (dimrPhi)
      !  -- row (line) dimension dim_line_design, dimr_Phi
      !  -- column dimension     dim_col_design,  dimc_Phi
      !--------------------------------------------------------
      use temporary_data_cov, only: dim_xdesc
      use module_db_setup, only: iconf_data
      use derived_types, only: config_desc, config_real

      implicit none

      integer, intent(out) :: dim_col_design, dim_line_design
      character(len=*), dimension(:), intent(in)     :: classes_for_sigma
      real(kind_double), dimension(:), allocatable, intent(out)  :: mean
      integer  :: ic, iclass, icount

      _NAMECURRENT_("get_dimension_design_Phi_by_class")

      _MLD_BEGIN_
      dimr_Phi = dim_xdesc
      dim_line_design = dimr_Phi
      if (allocated(mean)) deallocate (mean); allocate (mean(dimr_Phi))
      mean(:) = 0.d0

      icount = 0
      do iclass = 1, size(classes_for_sigma)
         do ic = 1, iconf_data
            if ((classes_for_sigma(iclass) == config_real(ic)%class)) then
               if (.not. (allocated(config_desc(ic)%energy))) cycle
               icount = icount + config_real(ic)%nat
               mean(:) = mean(:) + SUM(config_desc(ic)%energy(:, :), DIM=2)
            end if
         end do
      end do
      dim_col_design = icount
      dimc_Phi = dim_col_design
      mean(:) = mean(:)/dimc_Phi
      _MLD_END_
   end subroutine get_dimension_design_Phi_by_class

   subroutine prepare_design_Phi_by_class(classes_for_sigma, mean, Phi)
      use module_db_setup, only: iconf_data
      use derived_types, only: config_desc, config_real

      implicit none

      character(len=*), dimension(:), intent(in)     :: classes_for_sigma
      real(kind_double), dimension(:,:), intent(inout)     :: Phi
      real(kind_double), dimension(:), intent(inout)     :: mean
      integer  :: ia, ic, iclass, icount

      _NAMECURRENT_("prepare_design_Phi_by_class")

      _MLD_BEGIN_
      icount = 0
      do iclass = 1, size(classes_for_sigma)
         do ic = 1, iconf_data
            if ((classes_for_sigma(iclass) == config_real(ic)%class)) then
               if (allocated(config_desc(ic)%energy)) then
                  do ia = 1, config_real(ic)%nat
                     icount = icount + 1
                     Phi(:, icount) = config_desc(ic)%energy(:, ia) - mean(:)
                  end do
               end if
            end if ! class
         end do
      end do
      _MLD_END_
   end subroutine prepare_design_Phi_by_class

   subroutine compute_serial_Phi (classes_for_sigma, dim_line_design, dim_col_design, mean, Phi_mat)
      integer, intent(inout) :: dim_line_design, dim_col_design
      character(len=*), dimension(:), intent(inout) :: classes_for_sigma
      real(kind_double), allocatable, dimension(:,:), intent(inout) :: Phi_mat
      real(kind_double), allocatable, dimension(:), intent(inout) :: mean

      _NAMECURRENT_("compute_serial_Phi")

      _MLD_BEGIN_
      call get_dimension_design_Phi_by_class(classes_for_sigma, mean, dim_line_design, dim_col_design)
      ! Here the design matrix is defined as DxM
      call log_info("ML: before allocate Phi matrix size.....:  "//vtoa(dimr_Phi)//vtoa(dimc_Phi))

      if (allocated(Phi_mat)) deallocate (Phi_mat)
      allocate (Phi_mat(dimr_Phi, dimc_Phi))

      call prepare_design_Phi_by_class(classes_for_sigma, mean, Phi_mat)
      call log_info("ML: after allocate Phi dimensions given by size...:  "//vtoa(size(Phi_mat, 1))//vtoa(size(Phi_mat, 2)))

      _MLD_END_
   end subroutine compute_serial_Phi

   subroutine compute_Sigma_class_by_Phi(Phi_mat, Sigma_sample, Sigma_sample_inv)

      use ml_in_ndm_module, only: lambda_krr
      use math, only: dreal_inverse, dreal_matmul
      use mld_logger, only: log_info
      use mld_string, only: vtoa
      use mld_mpi, only: mld_mpi_abort
      implicit none

      !integer, intent(in)  :: no_of_classes_for_sigma
      !character(len=2), dimension(no_of_classes_for_sigma), intent(in)     :: classes_for_sigma
      real(kind_double), dimension(:,:), allocatable, intent(inout)  :: Sigma_sample, Sigma_sample_inv
      real(kind_double), dimension(:,:), allocatable, intent(inout)  :: Phi_mat
      real(kind_double), dimension(:,:), allocatable :: transpose_phi
      integer  :: id

      ! Sigma is DXD
      dimr_Phi = size(Phi_mat, 1)
      if (allocated(Sigma_sample)) deallocate (Sigma_sample)
      allocate (Sigma_sample(dimr_Phi, dimr_Phi))

      call log_info("ML: allocate Sigma, the  sample covariance matrix......:   " //&
         vtoa(size(Sigma_sample, 1))//vtoa(size(Sigma_sample, 2)))


      allocate(transpose_phi(size(Phi_mat, 2), size(Phi_mat, 1)))
      transpose_phi = transpose(Phi_mat)
      call dreal_matmul(Phi_mat, size(Phi_mat, 1), size(Phi_mat, 2), &
         transpose_phi, size(Phi_mat, 2), size(Phi_mat, 1), &
         Sigma_sample, size(Phi_mat, 1), size(Phi_mat, 1))
      deallocate(transpose_phi)

      if (size(Sigma_sample, 1) < 2) then
         call mld_mpi_abort("ML: ......... no covariance matrix for one/zero sample. Increase number of samples.")
      end if


      Sigma_sample(:, :) = Sigma_sample(:, :)/dble(size(Sigma_sample, 1) - 1)

      lambda_krr = 1.d-9

      if (allocated(Sigma_sample_inv)) deallocate (Sigma_sample_inv)
      allocate (Sigma_sample_inv(size(Sigma_sample, 1), size(Sigma_sample, 2)))

      if (lambda_krr >= 0.d0) then
         do id = 1, size(Sigma_sample, 1)
            Sigma_sample(id, id) = Sigma_sample(id, id) + lambda_krr
         end do
      end if

      call dreal_inverse(size(Sigma_sample, 1), Sigma_sample, Sigma_sample_inv)

   end subroutine compute_Sigma_class_by_Phi


end module module_Phi



module module_compute_cur
   use module_kind_variables, only: kind_double

   implicit none
contains
!
!
   subroutine compute_leverage_column_score(sca_AA, desc_sca_AA, &
      nbr_AA, nbc_AA, dimr_sca_AA, dimc_sca_AA, &
      l_dimr_sca_AA, l_dimc_sca_AA,  &
      cur_kval, cur_cval, cur_eps, &
      rank_sca_AA, pcol)
      use module_fit_ScaMatrix, only: scalapack_SVD_decomposition
      use module_scalapack_tools, only: allocation_scalapack_matrix
      use module_ml_scalapack, only: debug_scalapack, myrow, mycol, nprow, npcol, context, iam
      use ml_in_ndm_module, only : svd_rcond

      real(kind_double), dimension(:,:), allocatable, intent(in) :: sca_AA
      integer, dimension(:) , allocatable, intent(in) :: desc_sca_AA
      integer, intent(in) :: nbr_AA, nbc_AA, dimr_sca_AA, dimc_sca_AA
      integer, intent(in) :: l_dimr_sca_AA, l_dimc_sca_AA
      integer , intent(inout) :: cur_kval, cur_cval
      real(kind_double), intent(in) ::  cur_eps
      integer, intent(out)  :: rank_sca_AA
      real(kind_double), dimension(:), allocatable, intent(inout) :: pcol
      real(kind_double), external :: pdlange

      ! local variables ...
      ! sca_VT utilities ...
      real(kind_double), dimension(:,:), allocatable :: sca_VT
      integer, dimension(:) , allocatable :: desc_sca_VT
      integer  ::  nbr_VT, nbc_VT, dimr_sca_VT, dimc_sca_VT
      integer  :: l_dimr_sca_VT, l_dimc_sca_VT
      real(kind_double), dimension(:,:), allocatable :: sca_UU
      ! sca_UU utilities ...
      integer, dimension(:), allocatable :: desc_sca_UU
      integer :: nbr_UU, nbc_UU, dimr_sca_UU, dimc_sca_UU
      integer :: l_dimr_sca_UU, l_dimc_sca_UU


      ! sca_vec utilities ...
      integer :: nbr_vec, nbc_vec, dimr_sca_vec, dimc_sca_vec
      integer :: l_dimr_sca_vec, l_dimc_sca_vec
      integer, dimension(:) , allocatable:: desc_sca_vec
      real(kind_double), dimension(:,:), allocatable :: sca_vec
      integer :: ik, isel
      real(kind_double) :: norm2

      if (svd_rcond < 0) then
         svd_rcond = 100.d0*epsilon(1.d0)
      end if

      call scalapack_SVD_decomposition(sca_AA, desc_sca_AA, &
         nbr_AA, nbc_AA, dimr_sca_AA, dimc_sca_AA, &
         l_dimr_sca_AA, l_dimc_sca_AA,  &
         sca_UU, desc_sca_UU, &
         nbr_UU, nbc_UU, dimr_sca_UU, dimc_sca_UU, &
         l_dimr_sca_UU, l_dimc_sca_UU,  &
         sca_VT, desc_sca_VT, &
         nbr_VT, nbc_VT, dimr_sca_VT, dimc_sca_VT, &
         l_dimr_sca_VT, l_dimc_sca_VT,  &
         svd_rcond, rank_sca_AA)

      ! let's define CUR stuff ...
      ! k value - k order of SVD decomposition
      !cur_kval  = rank_sca_AA

      ! rank_sca_AA should be lower that min(Mline, Nline).
      if (cur_kval < 0 ) then
         cur_kval = rank_sca_AA
      end if
      if (cur_kval > rank_sca_AA) then
         if (iam == 0) then
            write (6, '("ML: CUR column selections cannot be larger than the rank  ....:", i4, i9, i9, i9, d20.10,i6)') &
               cur_kval, rank_sca_AA
         end if
         cur_kval = rank_sca_AA
      end if

      if (cur_cval < 0) then
         cur_cval = int(dble(cur_kval) * log(dble(cur_kval)) /cur_eps**2)
      end if

      if ((iam ==0).and.debug_scalapack) then
         write(6,'("ML: CUR values for c-selection, k-dimension, cur_eps (in matrix) ......:",2i8,d15.5)') &
            cur_cval, cur_kval, cur_eps
      end if


      nbr_vec =  min(nbr_AA, nbc_AA)
      nbc_vec =  1
      dimr_sca_vec = cur_kval
      dimc_sca_vec = 1
      call allocation_scalapack_matrix (" sca_vec ", sca_vec, desc_sca_vec, dimr_sca_vec, dimc_sca_vec, &
         l_dimr_sca_vec, l_dimc_sca_vec, nbr_vec, nbc_vec, &
         myrow, mycol, nprow, npcol, context, iam  )

      isel = 0
      ! Mcolm is equal to dimr_sca_VT
      if (allocated(pcol)) deallocate(pcol)
      allocate(pcol(dimc_sca_VT))
      ! VT is the transpose(V) so it is ........
      !                        ----- M_a -----
      !                       |               |
      !                       |               |
      !  VT=  transpose(V) =  r               |
      !                       |               |
      !                       |               |
      !                        ---------------
      do ik = 1, dimc_sca_VT
         !old! call pdcopy(cur_kval, sca_VT, 1, ik, desc_sca_VT, 1, sca_vec, 1, 1, desc_sca_vec, 1)
         !old! call pdnrm2(cur_kval, norm2, sca_vec, 1, 1, desc_sca_vec, 1)
         norm2 = pdlange('F', cur_kval, 1, sca_VT, 1, ik, desc_sca_VT, 1)
         pcol(ik) = norm2**2 / dble(cur_kval)
         !d!if (pcol(ik).le.1.d-15) then
         !d!   write(*,*) 'zero', iam, rangml, ik, pcol(ik)
         !d!end if
         !d!pcol(ik) = min(1.d0, norm2**2*dble(cur_cval) / dble(cur_kval))
      end do


      ! deallocate SVD utilities ...
      if (allocated(sca_UU))  deallocate (sca_UU)
      if (allocated(desc_sca_UU))  deallocate (desc_sca_UU)
      if (allocated(sca_VT))  deallocate (sca_VT)
      if (allocated(desc_sca_VT))  deallocate (desc_sca_VT)

   end  subroutine compute_leverage_column_score

   subroutine compute_sca_phia (classes_for_sigma, mean)
      use module_Phi, only : get_dimension_design_Phi_by_class
      use module_ml_scalapack, only: sca_phia, desc_sca_phia, dimr_sca_phia, dimc_sca_phia, &
         nbr_phia, nbc_phia, l_dimr_sca_phia, l_dimc_sca_phia, &
         myrow, mycol, nprow, npcol, context, iam, &
         nbr_predefined, nbc_predefined
      use module_scalapack_tools, only: allocation_scalapack_matrix
      use module_db_setup, only: iconf_data
      use derived_types, only: config_real, config_desc
      use module_cur, only:  info_mat
      use mld_logger

      implicit none


      character(len=2), dimension(:), allocatable, intent(in) :: classes_for_sigma
      real(kind_double), dimension(:), allocatable, intent(inout) :: mean
      real(kind_double), dimension(:), allocatable  :: ptmp
      integer :: iclass, ic, ik, ia, icount

      ! pay attention that column and row are inversed in phia with respect the Amat and Phi.
      call get_dimension_design_Phi_by_class(classes_for_sigma, mean, dimr_sca_phia, dimc_sca_phia)

      call log_info("ML:>>>-------------------------------------------------------------------------")
      call log_info("ML: before allocate phia matrix size.....:  "//vtoa(dimr_sca_phia)//vtoa(dimc_sca_phia))

      if (allocated(info_mat)) deallocate(info_mat)
      allocate(info_mat(dimc_sca_phia))



      if (nbr_predefined > 0) then
         nbr_phia = nbr_predefined
      else
         nbr_phia = -1
      end if

      if (nbc_predefined > 0) then
         nbc_phia = nbc_predefined
      else
         nbc_phia = -1
      end if

      call allocation_scalapack_matrix (" sca_phia ", sca_phia, desc_sca_phia, &
         dimr_sca_phia, dimc_sca_phia, &
         l_dimr_sca_phia, l_dimc_sca_phia, nbr_phia, nbc_phia, &
         myrow, mycol, nprow, npcol, context, iam)
      if (allocated(ptmp)) deallocate(ptmp)
      allocate(ptmp(dimr_sca_phia))

      icount = 0
      do iclass = 1, size(classes_for_sigma)
         do ic = 1, iconf_data
            if ((classes_for_sigma(iclass) == config_real(ic)%class)) then
               if (allocated(config_desc(ic)%energy)) then
                  do ia = 1, config_real(ic)%nat
                     icount = icount + 1
                     ptmp = config_desc(ic)%energy(:, ia) - mean(:)
                     info_mat(icount)%iconf = ic
                     info_mat(icount)%ia= ia
                     do ik = 1, dimr_sca_phia
                        call pdelset(sca_phia, ik, icount, desc_sca_phia, ptmp(ik))
                     end do
                  end do
               end if
            end if ! class
         end do
      end do

   end subroutine compute_sca_phia

   subroutine perform_cur_selection(pcol,  cur_cval,  selcol, dimc_sca_AA, no_of_selections)
#if (PARA)
      !use mpi
      use mld_mpi, only: comm_mld
#endif
      use ml_in_ndm_module, only: rangml
      use module_kind_variables, only: kind_double
      use module_ml_scalapack, only:  debug_scalapack
      use mld_logger

      implicit none
      real(kind_double), dimension(:), allocatable, intent(in)  :: pcol
      integer, intent(in) ::  dimc_sca_AA
      integer, intent(in) ::  cur_cval
      logical,  dimension(:), allocatable, intent(out)  :: selcol
      integer, intent(out) :: no_of_selections
      !integer, dimension(:), allocatable, intent(out) :: cur_info_selection
      integer :: isel, ik
      real(kind_double) :: uu, pp

      ! CUR selection ...
      if (allocated(selcol)) deallocate(selcol) ; allocate(selcol(dimc_sca_AA))
      selcol(:) = .false.

      if (rangml == 0) then
         isel = 0
         do ik = 1, dimc_sca_AA
            call random_number(uu)
            pp = min(1.d0, pcol(ik)*dble(cur_cval))
            if (pp >= uu) then
               selcol(ik) = .true.
               isel  = isel + 1
            end if
         end do
         no_of_selections = isel
      end if ! rangml ==0

      !TORC! call MPI_BCAST(no_of_selections, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, codeml)
      call comm_mld%bcast(0, no_of_selections)
      !TORC! call MPI_BCAST(selcol, size(selcol, 1), MPI_LOGICAL, 0, MPI_COMM_WORLD, codeml)
      call comm_mld%bcast(0, selcol)

      if (debug_scalapack) call log_info("ML: CUR number of selections .....:  "//vtoa(no_of_selections))


   end subroutine perform_cur_selection


   subroutine perform_cur_maha_selection(pcol,  cur_cval,  selcol, dimc_sca_AA, no_of_selections)
#if (PARA)
      !use mpi
      use mld_mpi, only: comm_mld
#endif
      use ml_in_ndm_module, only: rangml
      use module_kind_variables, only: kind_double
      use module_ml_scalapack, only:  debug_scalapack
      use mld_logger
      use derived_types, only: config_desc
      use module_cur, only:  info_mat
      use module_kernel, only: power_mcd

      implicit none
      real(kind_double), dimension(:), allocatable, intent(inout)  :: pcol
      integer, intent(in) ::  dimc_sca_AA
      integer, intent(in) ::  cur_cval
      logical,  dimension(:), allocatable, intent(out)  :: selcol
      integer, intent(out) :: no_of_selections
      !integer, dimension(:), allocatable, intent(out) :: cur_info_selection
      integer :: isel, ik, ic, ia
      real(kind_double) :: uu, pp, maha

      ! CUR selection ...
      if (allocated(selcol)) deallocate(selcol) ; allocate(selcol(dimc_sca_AA))
      selcol(:) = .false.
      if (allocated(pcol)) deallocate(pcol) ; allocate(pcol(dimc_sca_AA))
      if (rangml == 0) then
         isel = 0
         do ik = 1, dimc_sca_AA
            ic = info_mat(ik)%iconf
            ia = info_mat(ik)%ia
            maha=config_desc(ic)%stat_dist_maha(ia)**power_mcd
            pcol(ik) = exp(-maha**2/2.d0)
            !pcol(ik) = min(1.d0,exp(-maha**2/2.d0))
         end do

         pp = dble(cur_cval)/sum(pcol)
         pcol(:) = pcol(:)*pp

         do ik = 1, dimc_sca_AA
            call random_number(uu)
            if (pcol(ik) >= uu) then
               selcol(ik) = .true.
               isel  = isel + 1
            end if
         end do
         no_of_selections = isel
         pcol(:) = pcol(:)/dble(cur_cval)
      end if ! rangml ==0

      !TORC! call MPI_BCAST(pcol, size(pcol, 1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, codeml)
      call comm_mld%bcast(0, pcol)
      !TORC! call MPI_BCAST(no_of_selections, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, codeml)
      call comm_mld%bcast(0, no_of_selections)
      !TORC! call MPI_BCAST(selcol, size(selcol, 1), MPI_LOGICAL, 0, MPI_COMM_WORLD, codeml)
      call comm_mld%bcast(0, selcol)

      if (debug_scalapack) call log_info("ML: CUR number of selections .....:  "//vtoa(no_of_selections))


   end subroutine perform_cur_maha_selection

   subroutine  compute_full_CUR_matrix(pcol, selcol, col_no_of_selections, &
      prow, selrow, row_no_of_selections, &
      cur_ksel, cur_rval, cur_cval, cur_eps, &
      sca_AA, desc_sca_AA, &
      dimr_sca_AA, dimc_sca_AA, &
      l_dimr_sca_AA, l_dimc_sca_AA, &
      nbr_AA, nbc_AA, &
      sca_CC, desc_sca_CC, &
      dimr_sca_CC, dimc_sca_CC, &
      l_dimr_sca_CC, l_dimc_sca_CC, &
      nbr_CC, nbc_CC , &
      rank_sca_CC, &
      sca_UU, desc_sca_UU, &
      dimr_sca_UU, dimc_sca_UU, &
      l_dimr_sca_UU, l_dimc_sca_UU, &
      nbr_UU, nbc_UU, &
      sca_RR, desc_sca_RR, &
      dimr_sca_RR, dimc_sca_RR, &
      l_dimr_sca_RR, l_dimc_sca_RR, &
      nbr_RR, nbc_RR, &
      rank_sca_RR, norm_decomposition )

      use module_ml_scalapack, only:   myrow, mycol, nprow, npcol, context, iam
      use module_scalapack_tools, only: allocation_scalapack_matrix, pdgemm_nn
      use module_fit_ScaMatrix, only: scalapack_pseudo_inverse
      use mld_logger

      implicit none

      real(kind_double), dimension(:), allocatable, intent(in) :: pcol, prow
      logical, dimension(:), allocatable, intent(in) :: selcol, selrow
      integer, intent(in) :: row_no_of_selections, col_no_of_selections, &
         cur_ksel, cur_rval, cur_cval

      real(kind_double), intent(in) :: cur_eps
      real(kind_double), intent(out) :: norm_decomposition

      real(kind_double), dimension(:, :), allocatable, intent(in) :: sca_AA
      integer, dimension(:), allocatable, intent(in)  :: desc_sca_AA
      integer, intent(in)   :: nbr_AA, nbc_AA, dimr_sca_AA, dimc_sca_AA
      integer, intent(in)   :: l_dimr_sca_AA, l_dimc_sca_AA

      real(kind_double), dimension(:, :), allocatable, intent(inout) :: sca_CC
      integer, dimension(:), allocatable, intent(inout)  :: desc_sca_CC
      integer, intent(out)   :: nbr_CC, nbc_CC, dimr_sca_CC, dimc_sca_CC
      integer, intent(out)   :: l_dimr_sca_CC, l_dimc_sca_CC, rank_sca_CC

      real(kind_double), dimension(:, :), allocatable, intent(inout) :: sca_UU
      integer, dimension(:), allocatable, intent(inout)  :: desc_sca_UU
      integer, intent(out)   :: nbr_UU, nbc_UU, dimr_sca_UU, dimc_sca_UU
      integer, intent(out)   :: l_dimr_sca_UU, l_dimc_sca_UU

      real(kind_double), dimension(:, :), allocatable, intent(inout) :: sca_RR
      integer, dimension(:), allocatable, intent(inout)  :: desc_sca_RR
      integer, intent(out)   :: nbr_RR, nbc_RR, dimr_sca_RR, dimc_sca_RR
      integer, intent(out)   :: l_dimr_sca_RR, l_dimc_sca_RR, rank_sca_RR

      ! local definitions ...
      real(kind_double), dimension(:, :), allocatable :: sca_pinvCC
      integer, dimension(:), allocatable  :: desc_sca_pinvCC
      integer   :: nbr_pinvCC, nbc_pinvCC, dimr_sca_pinvCC, dimc_sca_pinvCC
      integer   :: l_dimr_sca_pinvCC, l_dimc_sca_pinvCC

      real(kind_double), dimension(:, :), allocatable :: sca_pinvRR
      integer, dimension(:), allocatable  :: desc_sca_pinvRR
      integer   :: nbr_pinvRR, nbc_pinvRR, dimr_sca_pinvRR, dimc_sca_pinvRR
      integer   :: l_dimr_sca_pinvRR, l_dimc_sca_pinvRR

      real(kind_double), dimension(:, :), allocatable :: sca_matt
      integer, dimension(:), allocatable  :: desc_sca_matt
      integer   :: nbr_matt, nbc_matt, dimr_sca_matt, dimc_sca_matt
      integer   :: l_dimr_sca_matt, l_dimc_sca_matt

      real(kind_double), dimension(:, :), allocatable :: sca_AAcur
      integer, dimension(:), allocatable  :: desc_sca_AAcur
      integer   :: nbr_AAcur, nbc_AAcur, dimr_sca_AAcur, dimc_sca_AAcur
      integer   :: l_dimr_sca_AAcur, l_dimc_sca_AAcur

      real(kind_double), external :: pdlange
      integer, external :: indxl2g
      integer :: ii, jj, jglobal, iglobal
      real(kind_double), parameter :: one=1.d0, zero=0.d0
      real(kind_double) :: norm2
      !TODOcur
      integer :: izozo 
      real(kind_double) :: rzozo
      logical, dimension(size(selrow)) :: iv1_zozo 
      logical, dimension(size(selcol)) :: iv2_zozo 


      _NAMECURRENT_("compute_full_cur_matrix")

      _MLD_BEGIN_

      rzozo = cur_eps
      izozo = cur_ksel  
      iv1_zozo = selrow 
      iv2_zozo = selcol 
      !izozo = selcol 
      !izozo = selrow 

      ! building sca_CC matrix a Mxc matrix
      nbr_CC =  nbr_AA
      nbc_CC =  nbc_AA
      dimr_sca_CC = dimr_sca_AA
      dimc_sca_CC = col_no_of_selections
      call allocation_scalapack_matrix (" sca_CC in ", sca_CC, desc_sca_CC, dimr_sca_CC, dimc_sca_CC, &
         l_dimr_sca_CC, l_dimc_sca_CC, nbr_CC, nbc_CC, &
         myrow, mycol, nprow, npcol, context, iam  )

      !sca_CC(1:Mline:1:col_no_of_selections) = sca_AA(1:Mline,1:col_no_of_selections)
      !sca_CC(1:Mline:ik) = sca_AA(:,ik)/sqrt(pc^2)
      do jj=1,l_dimc_sca_CC
         jglobal = indxl2g(jj, nbc_CC, mycol, 0, npcol)
         do ii = 1, l_dimr_sca_CC
            sca_CC(ii,jj) = sca_AA(ii,jj)/(dble(cur_cval)*dsqrt(pcol(jglobal)))
         end do
      end do

      !call log_info('02_debug ... before pseudo inverse')
      call scalapack_pseudo_inverse(sca_CC, desc_sca_CC, nbr_CC, nbc_CC, dimr_sca_CC, dimc_sca_CC, &
         l_dimr_sca_CC, l_dimc_sca_CC, rank_sca_CC, &
         sca_pinvCC, desc_sca_pinvCC, nbr_pinvCC, nbc_pinvCC ,  &
         dimr_sca_pinvCC, dimc_sca_pinvCC, l_dimr_sca_pinvCC, l_dimc_sca_pinvCC )
      !call log_info('02_debug ... after pseudo inverse')
      ! building  sca_RR matrix a rxn matrix.
      nbr_RR =  nbr_AA
      nbc_RR =  nbc_AA
      dimr_sca_RR = row_no_of_selections
      dimc_sca_RR = dimc_sca_AA
      call allocation_scalapack_matrix (" sca_RR in ", sca_RR, desc_sca_RR, dimr_sca_RR, dimc_sca_RR, &
         l_dimr_sca_RR, l_dimc_sca_RR, nbr_RR, nbc_RR, &
         myrow, mycol, nprow, npcol, context, iam  )

      do jj=1,l_dimc_sca_RR
         do ii = 1, l_dimr_sca_RR
            iglobal = indxl2g(ii, nbr_RR, myrow, 0, nprow)
            sca_RR(ii,jj) = sca_AA(ii,jj)/(dble(cur_rval)*dsqrt(prow(iglobal)))
         end do
      end do
      !call log_info('03_debug ... before pseudo inverse')
      call scalapack_pseudo_inverse(sca_RR, desc_sca_RR, nbr_RR, nbc_RR, dimr_sca_RR, dimc_sca_RR, &
         l_dimr_sca_RR, l_dimc_sca_RR, rank_sca_RR, &
         sca_pinvRR, desc_sca_pinvRR, nbr_pinvRR, nbc_pinvRR ,  &
         dimr_sca_pinvRR, dimc_sca_pinvRR, l_dimr_sca_pinvRR, l_dimc_sca_pinvRR )
      !call log_info('02_debug ... after pseudo inverse')
      ! compute UU ... U(cxr) = C+(cxm)A(mxn)R+(nxr)

      ! matt (mxr) =  A(mxn)R+(nxr)
      nbr_matt =  nbr_AA
      nbc_matt =  nbc_AA
      dimr_sca_matt = dimr_sca_AA
      dimc_sca_matt = row_no_of_selections
      call allocation_scalapack_matrix (" sca_matt in ", sca_matt, desc_sca_matt, dimr_sca_matt, dimc_sca_matt, &
         l_dimr_sca_matt, l_dimc_sca_matt, nbr_matt, nbc_matt, &
         myrow, mycol, nprow, npcol, context, iam  )
      ! mxn = mxk * kxn
      ! call pdgemm(transa, transb, m, n, k, alpha, a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc)
      call pdgemm_nn(dimr_sca_matt, dimc_sca_matt, dimc_sca_AA, one, sca_AA, desc_sca_AA, &
         sca_pinvRR, desc_sca_pinvRR, &
         zero, sca_matt, desc_sca_matt)
      !call log_info('04_debug ... after pdgemm')
      nbr_UU =  nbc_AA
      nbc_UU =  nbr_AA
      dimr_sca_UU = col_no_of_selections
      dimc_sca_UU = row_no_of_selections
      call allocation_scalapack_matrix (" sca_UU in ", sca_UU, desc_sca_UU, dimr_sca_UU, dimc_sca_UU, &
         l_dimr_sca_UU, l_dimc_sca_UU, nbr_UU, nbc_UU, &
         myrow, mycol, nprow, npcol, context, iam  )
      ! U(cxr) = C+(cxm) matt(mxr)
      call pdgemm_nn(dimr_sca_UU, dimc_sca_UU, dimr_sca_matt, one, sca_pinvCC, desc_sca_pinvCC, &
         sca_matt, desc_sca_matt, &
         zero, sca_UU, desc_sca_UU)
      if (allocated(sca_matt)) deallocate(sca_matt)
      if (allocated(desc_sca_matt)) deallocate(desc_sca_matt)

      ! compute the norm of the reconstruction A - (A_{k,c,r} = CUR)
      nbr_matt =  nbc_AA
      nbc_matt =  nbr_AA
      dimr_sca_matt = col_no_of_selections
      dimc_sca_matt = dimc_sca_AA
      call allocation_scalapack_matrix (" sca_matt in 2 ", sca_matt, desc_sca_matt, dimr_sca_matt, dimc_sca_matt, &
         l_dimr_sca_matt, l_dimc_sca_matt, nbr_matt, nbc_matt, &
         myrow, mycol, nprow, npcol, context, iam  )
      ! mxn = mxk * kxn
      ! call pdgemm(transa, transb, m, n, k, alpha, a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc)
      ! matt (cxn) =  U(cxr)R(rxn)
      call pdgemm_nn(dimr_sca_matt, dimc_sca_matt, dimr_sca_RR, one, sca_UU, desc_sca_UU, &
         sca_RR, desc_sca_RR, &
         zero, sca_matt, desc_sca_matt)
      !call log_info('05_debug ... after pdgemm')

      call blacs_barrier(context,'A')
      nbr_AAcur =  nbr_AA
      nbc_AAcur =  nbc_AA
      dimr_sca_AAcur = dimr_sca_AA
      dimc_sca_AAcur = dimc_sca_AA
      call allocation_scalapack_matrix (" sca_AAcur in ", sca_AAcur, desc_sca_AAcur, dimr_sca_AAcur, dimc_sca_AAcur, &
         l_dimr_sca_AAcur, l_dimc_sca_AAcur, nbr_AAcur, nbc_AAcur, &
         myrow, mycol, nprow, npcol, context, iam  )

      !AAcur (mxn)=C(mxc) matt(cxn)
      call pdgemm_nn(dimr_sca_AAcur, dimc_sca_AAcur, dimr_sca_matt, one, sca_CC, desc_sca_CC, &
         sca_matt, desc_sca_matt, &
         zero, sca_AAcur, desc_sca_AAcur)

      !call log_info('06_debug ... after pdgemm')


      call blacs_barrier(context,'A')
      if (allocated(sca_matt)) deallocate(sca_matt)
      if (allocated(desc_sca_matt)) deallocate(desc_sca_matt)

      nbr_matt =  nbr_AA
      nbc_matt =  nbc_AA
      dimr_sca_matt = dimr_sca_AA
      dimc_sca_matt = dimc_sca_AA
      call allocation_scalapack_matrix (" sca_matt in 3", sca_matt, desc_sca_matt, dimr_sca_matt, dimc_sca_matt, &
         l_dimr_sca_matt, l_dimc_sca_matt, nbr_matt, nbc_matt, &
         myrow, mycol, nprow, npcol, context, iam  )
      do jj =1, l_dimc_sca_AA
         do ii = 1, l_dimr_sca_AA
            sca_matt(ii,jj) = sca_AA(ii,jj) - sca_AAcur(ii,jj)
         end do
      end do
      !norm2 =  pdlange('F', dimr_sca_AA, dimc_sca_AA, sca_AA, 1, 1, desc_sca_AA, 1)
      !write(*,*) 'AA norm ........:  ', norm2
      !norm2 =  pdlange('F', dimr_sca_AAcur, dimc_sca_AAcur, sca_AAcur, 1, 1, desc_sca_AAcur, 1)
      !write(*,*) 'AAcur norm ........:  ', norm2

      norm2 =  pdlange('F', dimr_sca_matt, dimc_sca_matt, sca_matt, 1, 1, desc_sca_matt, 1)
      norm_decomposition=norm2
      call blacs_barrier(context,'A')

   end subroutine  compute_full_CUR_matrix



   !dg! subroutine compute_cur_decomposition (classes_for_sigma, mean)
   subroutine compute_cur_decomposition(sca_phia, desc_sca_phia, dimr_sca_phia, dimc_sca_phia, &
      l_dimr_sca_phia, l_dimc_sca_phia, nbr_phia, nbc_phia, &
      cur_kval, cur_cval, cur_rval, cur_eps, &
      pcol, selcol,col_no_of_selections,  prow, selrow, row_no_of_selections, &
      cur_info_selection)


      use ml_in_ndm_module, only: kind_double
      use module_ml_scalapack, only: myrow, mycol, nprow, npcol, context, iam
      use module_cur, only: rank_sca_phia,  rank_sca_phia_T, &
         cur_no_of_samples, pcol_sample, prow_sample, selcol_sample, selrow_sample, &
         norow_sample, nocol_sample
      use module_scalapack_tools, only: allocation_scalapack_matrix
      use module_fit_ScaMatrix, only: scalapack_pseudo_inverse
      use mld_logger
      use mld_string

      implicit none

      real(kind_double), dimension(:, :), allocatable, intent(in) :: sca_phia
      integer, dimension(:), allocatable, intent(in)  :: desc_sca_phia
      integer, intent(in)   :: nbr_phia, nbc_phia, dimr_sca_phia, dimc_sca_phia
      integer, intent(in)   :: l_dimr_sca_phia, l_dimc_sca_phia
      integer, intent(inout) :: cur_kval
      integer, intent(inout) ::  cur_cval, cur_rval
      real(kind_double), intent(inout) ::  cur_eps
      real(kind_double), dimension(:), allocatable, intent(out) :: pcol, prow
      logical, dimension(:), allocatable, intent(out):: selcol, selrow
      integer, dimension(:), allocatable, intent(out):: cur_info_selection
      integer, intent(out) :: col_no_of_selections, row_no_of_selections

      !local variables ...
      real(kind_double), dimension(:, :), allocatable :: sca_phia_T
      integer, dimension(:), allocatable  :: desc_sca_phia_T
      integer   :: nbr_phia_T, nbc_phia_T, dimr_sca_phia_T, dimc_sca_phia_T
      integer   :: l_dimr_sca_phia_T, l_dimc_sca_phia_T

      real(kind_double), dimension(:, :), allocatable :: sca_CC
      integer, dimension(:), allocatable  :: desc_sca_CC
      integer   :: nbr_CC, nbc_CC, dimr_sca_CC, dimc_sca_CC
      integer   :: l_dimr_sca_CC, l_dimc_sca_CC, rank_sca_CC

      real(kind_double), dimension(:, :), allocatable :: sca_RR
      integer, dimension(:), allocatable  :: desc_sca_RR
      integer   :: nbr_RR, nbc_RR, dimr_sca_RR, dimc_sca_RR
      integer   :: l_dimr_sca_RR, l_dimc_sca_RR, rank_sca_RR

      real(kind_double), dimension(:, :), allocatable :: sca_UU
      integer, dimension(:), allocatable  :: desc_sca_UU
      integer   :: nbr_UU, nbc_UU, dimr_sca_UU, dimc_sca_UU
      integer   :: l_dimr_sca_UU, l_dimc_sca_UU

      integer ::   ik, ii, icur, isel
      integer, external :: indxl2g
      real(kind_double), dimension(:), allocatable :: cur_score_sample
      real(kind_double) :: cur_score


      real(kind_double) :: one=1.d0, zero=0.d0

      _NAMECURRENT_("compute_cur_decomposition")

      _MLD_BEGIN_


      !dn! call  compute_sca_phia (classes_for_sigma, mean)

      ! >>>>>>>>>>>>>>>>>>>>>compute levarage score for column ------------------!
      call  compute_leverage_column_score(sca_phia, desc_sca_phia, &
         nbr_phia, nbc_phia, dimr_sca_phia, dimc_sca_phia, &
         l_dimr_sca_phia, l_dimc_sca_phia, &
         cur_kval, cur_rval, cur_eps, &
         rank_sca_phia, pcol)

      if (allocated(pcol_sample)) deallocate(pcol_sample)
      allocate(pcol_sample(cur_no_of_samples, dimc_sca_phia))
      if (allocated(selcol_sample)) deallocate(selcol_sample)
      allocate(selcol_sample(cur_no_of_samples, dimc_sca_phia))
      if (allocated(nocol_sample)) deallocate(nocol_sample)
      allocate(nocol_sample(cur_no_of_samples))

      if (iam ==0) then
         write(6,'("ML: CUR values for c-selection, k-dimension, cur_eps ......:",2i8,d15.5)') &
            cur_cval, cur_kval, cur_eps
      end if

      do ii = 1, cur_no_of_samples
         call  perform_cur_selection(pcol, cur_cval, selcol, dimc_sca_phia, col_no_of_selections)
         pcol_sample(ii, :) = pcol(:)
         selcol_sample(ii,:) = selcol(:)
         nocol_sample(ii) = col_no_of_selections
         !call log_info("ML: CUR trail no and no of col sel .....:  "//vtoa(ii)//vtoa(col_no_of_selections))
      end do
      ! <<<<<<<<<<<<<<<<<<<<<compute levarage score for column ------------------!

      ! >>>>>>>>>>>>>>>>>>>>>compute levarage score for rows  -------------------!
      nbr_phia_T =  nbc_phia
      nbc_phia_T =  nbr_phia
      dimr_sca_phia_T = dimc_sca_phia
      dimc_sca_phia_T = dimr_sca_phia
      call allocation_scalapack_matrix (" sca_phia_T ", sca_phia_T, desc_sca_phia_T, dimr_sca_phia_T, dimc_sca_phia_T, &
         l_dimr_sca_phia_T, l_dimc_sca_phia_T, nbr_phia_T, nbc_phia_T, &
         myrow, mycol, nprow, npcol, context, iam  )
      ! PvTRAN( M, N, ALPHA, A, IA, JA, DESCA, BETA, C, IC, JC, DESCC )
      ! sub( C ) = beta * sub( C ) + alpha * op( sub( A ) )
      call pdtran(dimr_sca_phia_T, dimc_sca_phia_T, one, sca_phia, 1, 1 , desc_sca_phia, zero, sca_phia_T, 1, 1, desc_sca_phia_T)
      call  compute_leverage_column_score(sca_phia_T, desc_sca_phia_T, &
         nbr_phia_T, nbc_phia_T, dimr_sca_phia_T, dimc_sca_phia_T, &
         l_dimr_sca_phia_T, l_dimc_sca_phia_T, &
         cur_kval, cur_rval, cur_eps,  &
         rank_sca_phia_T, prow)
      if (rank_sca_phia /= rank_sca_phia_T) then
         call log_warning("ML: warning, somehow strange situation rank(A) /= rank(A^T). Weird, isn't it ? in subroutine "//NAMECURRENT)
      end if

      if (allocated(prow_sample)) deallocate(prow_sample)
      allocate(prow_sample(cur_no_of_samples, dimc_sca_phia_T))
      if (allocated(selrow_sample)) deallocate(selrow_sample)
      allocate(selrow_sample(cur_no_of_samples, dimc_sca_phia_T))
      if (allocated(norow_sample)) deallocate(norow_sample)
      allocate(norow_sample(cur_no_of_samples))

      if (iam ==0) then
         write(6,'("ML: CUR values for r-selection, k-dimension, cur_eps ......:",2i8,d15.5)') &
            cur_rval, cur_kval, cur_eps
      end if

      do ii = 1, cur_no_of_samples
         call  perform_cur_selection(prow, cur_rval, selrow, dimc_sca_phia_T, row_no_of_selections)
         prow_sample(ii, :) = prow(:)
         selrow_sample(ii,:) = selrow(:)
         norow_sample(ii) = row_no_of_selections
         call log_info("ML: CUR no and no of row and col sel ...:  "//vtoa(ii)//vtoa(row_no_of_selections)//vtoa(nocol_sample(ii)))
         !call log_info("ML: CUR trial no and no of row sel .....:  "//vtoa(ii)//vtoa(row_no_of_selections))
      end do
      ! <<<<<<<<<<<<<<<<<<<<<compute levarage score for rows  -------------------!
      if (allocated(cur_score_sample)) deallocate(cur_score_sample)
      allocate(cur_score_sample(cur_no_of_samples))
      do ii = 1, cur_no_of_samples
         pcol(:)=pcol_sample(ii,:)
         selcol(:) = selcol_sample(ii,:)
         prow(:)=prow_sample(ii,:)
         selrow(:)=selrow_sample(ii,:)
         col_no_of_selections = nocol_sample(ii)
         row_no_of_selections = norow_sample(ii)
         call  compute_full_CUR_matrix(pcol, selcol, col_no_of_selections, &
            prow, selrow, row_no_of_selections, &
            cur_kval, cur_rval, cur_cval, cur_eps, &
            sca_phia, desc_sca_phia, &
            dimr_sca_phia, dimc_sca_phia, &
            l_dimr_sca_phia, l_dimc_sca_phia, &
            nbr_phia, nbc_phia, &
            sca_CC, desc_sca_CC, &
            dimr_sca_CC, dimc_sca_CC, &
            l_dimr_sca_CC, l_dimc_sca_CC, &
            nbr_CC, nbc_CC , &
            rank_sca_CC, &
            sca_UU, desc_sca_UU, &
            dimr_sca_UU, dimc_sca_UU, &
            l_dimr_sca_UU, l_dimc_sca_UU, &
            nbr_UU, nbc_UU, &
            sca_RR, desc_sca_RR, &
            dimr_sca_RR, dimc_sca_RR, &
            l_dimr_sca_RR, l_dimc_sca_RR, &
            nbr_RR, nbc_RR, &
            rank_sca_RR, cur_score )

         cur_score_sample(ii) = cur_score
         call log_info("ML: CUR trial no and Frobenius score ...: "//vtoa(ii)//"  "//vtoa(cur_score))
      end do

      icur=minloc(cur_score_sample, dim=1)
      call log_info("ML: ----------- the winner is ----------> "//vtoa(icur))

      ! kernel selections ... prepare cur_info_selection
      col_no_of_selections = nocol_sample(icur)
      selcol(:) = selcol_sample(icur,:)
      if (allocated(cur_info_selection)) deallocate(cur_info_selection)
      allocate(cur_info_selection(col_no_of_selections))
      isel = 0
      do ik = 1, dimc_sca_phia
      if (selcol(ik) .eqv. .false.) cycle
         isel = isel + 1
         cur_info_selection(isel) = ik
      end do

      _MLD_END_
   end subroutine compute_cur_decomposition

   subroutine compute_cur_maha_decomposition(sca_phia, desc_sca_phia, dimr_sca_phia, dimc_sca_phia, &
      l_dimr_sca_phia, l_dimc_sca_phia, nbr_phia, nbc_phia, &
      cur_kval, cur_cval, cur_rval, cur_eps, &
      pcol, selcol,col_no_of_selections,  prow, selrow, row_no_of_selections, &
      cur_info_selection)


      use ml_in_ndm_module, only: kind_double
      use module_ml_scalapack, only: myrow, mycol, nprow, npcol, context, iam
      use module_cur, only:   rank_sca_phia_T, &
         cur_no_of_samples, pcol_sample, prow_sample, selcol_sample, selrow_sample, &
         norow_sample, nocol_sample
      use module_scalapack_tools, only: allocation_scalapack_matrix
      use module_fit_ScaMatrix, only: scalapack_pseudo_inverse
      use mld_logger
      use mld_string

      implicit none

      real(kind_double), dimension(:, :), allocatable, intent(in) :: sca_phia
      integer, dimension(:), allocatable, intent(in)  :: desc_sca_phia
      integer, intent(in)   :: nbr_phia, nbc_phia, dimr_sca_phia, dimc_sca_phia
      integer, intent(in)   :: l_dimr_sca_phia, l_dimc_sca_phia
      integer, intent(inout) :: cur_kval
      integer, intent(inout) ::  cur_cval, cur_rval
      real(kind_double), intent(inout) ::  cur_eps
      real(kind_double), dimension(:), allocatable, intent(out) :: pcol, prow
      logical, dimension(:), allocatable, intent(out):: selcol, selrow
      integer, dimension(:), allocatable, intent(out):: cur_info_selection
      integer, intent(out) :: col_no_of_selections, row_no_of_selections

      !local variables ...
      real(kind_double), dimension(:, :), allocatable :: sca_phia_T
      integer, dimension(:), allocatable  :: desc_sca_phia_T
      integer   :: nbr_phia_T, nbc_phia_T, dimr_sca_phia_T, dimc_sca_phia_T
      integer   :: l_dimr_sca_phia_T, l_dimc_sca_phia_T

      real(kind_double), dimension(:, :), allocatable :: sca_CC
      integer, dimension(:), allocatable  :: desc_sca_CC
      integer   :: nbr_CC, nbc_CC, dimr_sca_CC, dimc_sca_CC
      integer   :: l_dimr_sca_CC, l_dimc_sca_CC, rank_sca_CC

      real(kind_double), dimension(:, :), allocatable :: sca_RR
      integer, dimension(:), allocatable  :: desc_sca_RR
      integer   :: nbr_RR, nbc_RR, dimr_sca_RR, dimc_sca_RR
      integer   :: l_dimr_sca_RR, l_dimc_sca_RR, rank_sca_RR

      real(kind_double), dimension(:, :), allocatable :: sca_UU
      integer, dimension(:), allocatable  :: desc_sca_UU
      integer   :: nbr_UU, nbc_UU, dimr_sca_UU, dimc_sca_UU
      integer   :: l_dimr_sca_UU, l_dimc_sca_UU

      integer ::   ik, ii, icur, isel
      integer, external :: indxl2g
      real(kind_double), dimension(:), allocatable :: cur_score_sample
      real(kind_double) :: cur_score


      real(kind_double) :: one=1.d0, zero=0.d0

      _NAMECURRENT_("compute_cur_maha_decomposition")

      _MLD_BEGIN_


      !dn! call  compute_sca_phia (classes_for_sigma, mean)

      ! >>>>>>>>>>>>>>>>>>>>>compute levarage score for column ------------------!
      !call  compute_leverage_column_score(sca_phia, desc_sca_phia, &
      !        nbr_phia, nbc_phia, dimr_sca_phia, dimc_sca_phia, &
      !        l_dimr_sca_phia, l_dimc_sca_phia, &
      !        cur_kval, cur_rval, cur_eps, &
      !        rank_sca_phia, pcol)

      if (allocated(pcol_sample)) deallocate(pcol_sample)
      allocate(pcol_sample(cur_no_of_samples, dimc_sca_phia))
      if (allocated(selcol_sample)) deallocate(selcol_sample)
      allocate(selcol_sample(cur_no_of_samples, dimc_sca_phia))
      if (allocated(nocol_sample)) deallocate(nocol_sample)
      allocate(nocol_sample(cur_no_of_samples))

      if (iam ==0) then
         write(6,'("ML: CUR values for c-selection, k-dimension, cur_eps ......:",2i8,d15.5)') &
            cur_cval, cur_kval, cur_eps
      end if

      do ii = 1, cur_no_of_samples
         call  perform_cur_maha_selection(pcol, cur_cval, selcol, dimc_sca_phia, col_no_of_selections)
         pcol_sample(ii, :) = pcol(:)
         selcol_sample(ii,:) = selcol(:)
         nocol_sample(ii) = col_no_of_selections
         !call log_info("ML: CUR trail no and no of col sel .....:  "//vtoa(ii)//vtoa(col_no_of_selections))
      end do
      ! <<<<<<<<<<<<<<<<<<<<<compute levarage score for column ------------------!

      ! >>>>>>>>>>>>>>>>>>>>>compute levarage score for rows  -------------------!
      nbr_phia_T =  nbc_phia
      nbc_phia_T =  nbr_phia
      dimr_sca_phia_T = dimc_sca_phia
      dimc_sca_phia_T = dimr_sca_phia
      call allocation_scalapack_matrix (" sca_phia_T ", sca_phia_T, desc_sca_phia_T, dimr_sca_phia_T, dimc_sca_phia_T, &
         l_dimr_sca_phia_T, l_dimc_sca_phia_T, nbr_phia_T, nbc_phia_T, &
         myrow, mycol, nprow, npcol, context, iam  )
      ! PvTRAN( M, N, ALPHA, A, IA, JA, DESCA, BETA, C, IC, JC, DESCC )
      ! sub( C ) = beta * sub( C ) + alpha * op( sub( A ) )
      call pdtran(dimr_sca_phia_T, dimc_sca_phia_T, one, sca_phia, 1, 1 , desc_sca_phia, zero, sca_phia_T, 1, 1, desc_sca_phia_T)
      call  compute_leverage_column_score(sca_phia_T, desc_sca_phia_T, &
         nbr_phia_T, nbc_phia_T, dimr_sca_phia_T, dimc_sca_phia_T, &
         l_dimr_sca_phia_T, l_dimc_sca_phia_T, &
         cur_kval, cur_rval, cur_eps,  &
         rank_sca_phia_T, prow)
      ! if (rank_sca_phia /= rank_sca_phia_T) then
      !   call log_warning("ML: warning, somehow strange situation rank(A) /= rank(A^T). &
      !                     Weird, isn't it ? in subroutine "//NAMECURRENT)
      ! end if

      if (allocated(prow_sample)) deallocate(prow_sample)
      allocate(prow_sample(cur_no_of_samples, dimc_sca_phia_T))
      if (allocated(selrow_sample)) deallocate(selrow_sample)
      allocate(selrow_sample(cur_no_of_samples, dimc_sca_phia_T))
      if (allocated(norow_sample)) deallocate(norow_sample)
      allocate(norow_sample(cur_no_of_samples))

      if (iam ==0) then
         write(6,'("ML: CUR values for r-selection, k-dimension, cur_eps ......:",2i8,d15.5)') &
            cur_rval, cur_kval, cur_eps
      end if

      do ii = 1, cur_no_of_samples
         call  perform_cur_selection(prow, cur_rval, selrow, dimc_sca_phia_T, row_no_of_selections)
         prow_sample(ii, :) = prow(:)
         selrow_sample(ii,:) = selrow(:)
         norow_sample(ii) = row_no_of_selections
         call log_info("ML: CUR no and no of row and col sel ...:  "//vtoa(ii)//vtoa(row_no_of_selections)//vtoa(nocol_sample(ii)))
         !call log_info("ML: CUR trial no and no of row sel .....:  "//vtoa(ii)//vtoa(row_no_of_selections))
      end do
      ! <<<<<<<<<<<<<<<<<<<<<compute levarage score for rows  -------------------!
      if (allocated(cur_score_sample)) deallocate(cur_score_sample)
      allocate(cur_score_sample(cur_no_of_samples))
      do ii = 1, cur_no_of_samples
         pcol(:)=pcol_sample(ii,:)
         selcol(:) = selcol_sample(ii,:)
         prow(:)=prow_sample(ii,:)
         selrow(:)=selrow_sample(ii,:)
         col_no_of_selections = nocol_sample(ii)
         row_no_of_selections = norow_sample(ii)
         call  compute_full_CUR_matrix(pcol, selcol, col_no_of_selections, &
            prow, selrow, row_no_of_selections, &
            cur_kval, cur_rval, cur_cval, cur_eps, &
            sca_phia, desc_sca_phia, &
            dimr_sca_phia, dimc_sca_phia, &
            l_dimr_sca_phia, l_dimc_sca_phia, &
            nbr_phia, nbc_phia, &
            sca_CC, desc_sca_CC, &
            dimr_sca_CC, dimc_sca_CC, &
            l_dimr_sca_CC, l_dimc_sca_CC, &
            nbr_CC, nbc_CC , &
            rank_sca_CC, &
            sca_UU, desc_sca_UU, &
            dimr_sca_UU, dimc_sca_UU, &
            l_dimr_sca_UU, l_dimc_sca_UU, &
            nbr_UU, nbc_UU, &
            sca_RR, desc_sca_RR, &
            dimr_sca_RR, dimc_sca_RR, &
            l_dimr_sca_RR, l_dimc_sca_RR, &
            nbr_RR, nbc_RR, &
            rank_sca_RR, cur_score )

      cur_score_sample(ii) = cur_score
      call log_info("ML: CUR trial no and Frobenius score ...: "//vtoa(ii)//"  "//vtoa(cur_score))
    end do 

      icur=minloc(cur_score_sample, dim=1)
      call log_info("ML: ----------- the winner is ----------> "//vtoa(icur))

      ! kernel selections ... prepare cur_info_selection
      col_no_of_selections = nocol_sample(icur)
      selcol(:) = selcol_sample(icur,:)
      if (allocated(cur_info_selection)) deallocate(cur_info_selection)
      allocate(cur_info_selection(col_no_of_selections))
      isel = 0
      do ik = 1, dimc_sca_phia
      if (selcol(ik) .eqv. .false.) cycle
         isel = isel + 1
         cur_info_selection(isel) = ik
      end do

      _MLD_END_
   end subroutine compute_cur_maha_decomposition

end module module_compute_cur

