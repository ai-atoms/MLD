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

module module_post_descriptors
  use module_kind_variables, only: kind_double 
  implicit none 
  contains 
  subroutine compute_post_descriptors(icount,  ml_type, zbl_potential, i_start_at, i_final_at, imm, imm_neigh,  &
                                      desc_forces_local, subworld)
   use mpi
   use mld_mpi

  use ml_in_ndm_module, only: ml_type_krr
  use time_check_general, only: MY_MPI_WTIME
  use module_kernel, only: dim_kernel, time_for_desc_kernel
  use derived_types, only: config_desc, config_real
  use module_compute_kernel, only: compute_kernel, allocate_kernel
  use module_compute_zbl, only: compute_zbl, allocate_zbl, energy_force_stress_zbl
  integer, intent(inout) :: i_start_at, i_final_at
  integer, intent(in) :: icount, ml_type, imm, imm_neigh , subworld 
  logical, intent(in) :: desc_forces_local, zbl_potential
  real(kind_double) :: t00, t11 
  integer :: dim_reduce 

   ! complete with k2b kernel if needed
   !$! if (activate_k2b) then
   !$!    t00 = MY_MPI_WTIME()
   !$!    call allocate_kernel_k2b(icount, i_start_at, i_final_at, dim_kernel_2b, imm, imm_neigh, desc_forces_local)
   !$!    call compute_kernel_nbody(i_start_at, i_final_at, icount)
   !$!    dim_reduce = imm*dim_kernel_2b
   !$!    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%energy_k2b, dim_reduce, MPI_DOUBLE_PRECISION, MPI_SUM, subworld, codeml)
   !$!    t11 = MY_MPI_WTIME()
   !$!    time_for_desc_kernel_2b = time_for_desc_kernel_2b + (t11 - t00)
   !$! end if

   ! complete with pure kernel if needed
   if (ml_type == ml_type_krr) then
      t00 = MY_MPI_WTIME()
      call allocate_kernel(icount, i_start_at, i_final_at, dim_kernel, imm, imm_neigh, desc_forces_local)
      call compute_kernel(i_start_at, i_final_at, icount)
      dim_reduce = imm*dim_kernel
      call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%energy_kernel, dim_reduce, MPI_DOUBLE_PRECISION, MPI_SUM, subworld, codeml)
      t11 = MY_MPI_WTIME()
      time_for_desc_kernel = time_for_desc_kernel + (t11 - t00)
   end if

   !complete with zbl if nedded 
   if (zbl_potential) then
      call allocate_zbl(icount, i_start_at, i_final_at, desc_forces_local)
      call compute_zbl(i_start_at, i_final_at, config_real(icount)%n_neigh_zbl, config_real(icount)%kind_neigh_zbl, icount)
      call energy_force_stress_zbl (icount)
   end if
  end subroutine  compute_post_descriptors

end module module_post_descriptors

subroutine compute_descriptors(icount, post_desc)

   use mpi
   use mld_mpi !TODOMPI  , only: codeml, rangml 
   use my_mpi_subroutines, only: my_barrier_subworld

#ifdef MLD_NDM
   use gen_com_m_ml, ONLY: im, imm
#else
   use ondm_gen_com_m, ONLY: im, imm
#endif
   use ml_in_ndm_module, ONLY: debug, ml_type, &
      descriptor_type, descriptor_g2, descriptor_g3, descriptor_behler, &
      descriptor_afs, descriptor_pow_so3, &
      descriptor_pow_so4, descriptor_bispectrum_so4, descriptor_g2_bispectrum_so4, descriptor_g2_afs, &
      descriptor_mtp, descriptor_body, descriptor_pow_so3_3body, descriptor_ace, &
      i_start_at, i_final_at, &
      write_desc, write_desc_dump, read_desc_dump, &
      imm_neigh, descriptor_g2_pow_so4, &
      desc_forces, descriptor_milady, descriptor_ftnbody, descriptor_zetabody, &
      descriptor_tbind
   use temporary_data_cov, ONLY: dim_xdesc, dim_xdesc1, dim_xdesc2
   use set_limits, only: i_start_on_proc, i_end_on_proc, i_size_on_proc, &
      set_limit_for_atoms_with_MPI_grid
   use compute_g2_mod, only: compute_g2
   use compute_g3_mod, only: compute_g3
   use compute_dump_desc_mod, only: compute_dump_desc
   use compute_bispectrum_so4_mod, only: compute_bispectrum_so4
   use compute_afs_mod, only: compute_afs
   use compute_pow_so3_mod, only: compute_pow_so3

   !use angular_functions
   use derived_types, only: config_desc, config_real
   use module_kind_variables, only: kind_double

   use time_check_general, only: MY_MPI_WTIME
   use module_zbl, only: zbl_potential
   use module_kernel_2b, only: dim_kernel_2b, activate_k2b, time_for_desc_kernel_2b
   use module_chemical_species, only: img_weighted
   use mld_subworld, only: subworld, subrank, subworld_size
   use module_post_descriptors, only: compute_post_descriptors
   use module_compute_pow_so4, only: compute_pow_so4
   use module_compute_body_bonds_and_angles, only: compute_body_bonds_and_angles
   use module_compute_ftnbody, only: compute_ftnbody
   use module_compute_tbind, only: compute_tbind
   use module_compute_ace, only: compute_ace
   use module_compute_r_matrix, only: compute_r_matrix
   use module_compute_mtp, only: compute_mtp
   use module_compute_pow_so3_3body, only: compute_pow_so3_3body
   use module_compute_zetabody, only: compute_zetabody
   use mld_logger
   use mld_string
   use module_first_neighbour_pass, only: compute_max_n_neigh

   implicit none

   integer, intent(in)  :: icount
   logical,  intent(in) :: post_desc
   integer  :: dim_reduce, ip
   integer  :: saved_imm_neigh
   logical  :: desc_forces_local
   real(kind_double)    :: t00, t11
   _NAMECURRENT_("compute_descriptors")
   _MLD_BEGIN_

   ! depending on rang the index of atoms is distributed on procs
   desc_forces_local = desc_forces .and. (config_real(icount)%has_force .or. config_real(icount)%has_stress)

   ! set limits for atoms
   call set_limit_for_atoms_with_MPI_grid(config_real(icount)%nat, i_start_at, i_final_at)
   config_real(icount)%at_start = i_start_at
   config_real(icount)%at_final = i_final_at
   if (debug) then
      write (6, '("ML: compute_descriptors in compute_descriptor rangml im i_start_at, i_final_at  nf ",i4, 4i9)') rangml, im, i_start_at, i_final_at, i_final_at - i_start_at + 1
      if (subrank == 0) then
         do ip = 0, subworld_size - 1
            write (6, '("ML:i_size_on_proc ", 4i8)') ip, i_start_on_proc(ip), i_end_on_proc(ip), i_size_on_proc(ip)
         end do
      end if
      !call MPI_BARRIER(subworld, codeml)
      call my_barrier_subworld(codeml)
   end if

   ! Per-config imm_neigh: use actual max neighbour count from calc_neighbours
   ! Save original imm_neigh (needed by calc_neighbours for next config)
   saved_imm_neigh = imm_neigh
   call compute_max_n_neigh(icount, i_start_at, i_final_at)
   imm_neigh = config_real(icount)%max_n_neigh

   if (allocated(config_desc(icount)%n_neigh)) deallocate (config_desc(icount)%n_neigh); allocate (config_desc(icount)%n_neigh(imm))
   if (allocated(config_desc(icount)%kind_neigh)) deallocate (config_desc(icount)%kind_neigh); allocate (config_desc(icount)%kind_neigh(imm, imm_neigh))
   ! if (allocated(config_desc(icount)%incell)) deallocate(config_desc(icount)%incell) ; allocate(config_desc(icount)%incell(imm,imm_neigh))
   if (zbl_potential) then
      if (allocated(config_real(icount)%n_neigh_zbl)) deallocate (config_real(icount)%n_neigh_zbl); allocate (config_real(icount)%n_neigh_zbl(imm))
      if (allocated(config_real(icount)%kind_neigh_zbl)) deallocate (config_real(icount)%kind_neigh_zbl); allocate (config_real(icount)%kind_neigh_zbl(imm, imm_neigh))
   end if

   if (read_desc_dump) then

      config_desc(icount)%dim_desc1 = dim_xdesc1
      config_desc(icount)%dim_desc2 = dim_xdesc2
      config_desc(icount)%dim_desc = dim_xdesc

      if (.not.(img_weighted)) then
         if (allocated(config_desc(icount)%energy)) deallocate (config_desc(icount)%energy); allocate (config_desc(icount)%energy(dim_xdesc, imm))
         if (desc_forces_local) then
            if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
               if (allocated(config_desc(icount)%force)) deallocate (config_desc(icount)%force); allocate (config_desc(icount)%force(dim_xdesc, i_start_at:i_final_at, 0:imm_neigh, 3))
            end if
         end if
      end if

      call compute_dump_desc(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)
      call read_energy_descriptors_dump(icount)
#if(PARA)
      !call MPI_BARRIER(subworld, codeml)
      call my_barrier_subworld(codeml)
#endif
      call read_force_descriptors_dump(icount)
#if(PARA)
      !call MPI_BARRIER(subworld, codeml)
      call my_barrier_subworld(codeml)
#endif


   else

      ! here are mixed descriptors D1 + D2
      if ((descriptor_type == descriptor_behler) .or. &
         (descriptor_type == descriptor_g2_bispectrum_so4) .or. &
         (descriptor_type == descriptor_g2_pow_so4) .or. &
         (descriptor_type == descriptor_g2_afs)) then

         config_desc(icount)%dim_desc1 = dim_xdesc1
         config_desc(icount)%dim_desc2 = dim_xdesc2
         config_desc(icount)%dim_desc = dim_xdesc

         select case (descriptor_type)

          case (descriptor_behler)
            !D1 -> G2    D2 -> G3
            call compute_mixed_behler(icount)

          case (descriptor_g2_bispectrum_so4)
            !D1 -> G2    D2 -> bSO4
            call compute_mixed_g2_bispectrum_so4(icount)

          case (descriptor_g2_afs)
            !D1 -> G2    D2 -> bSO4
            call compute_mixed_g2_afs(icount)

          case (descriptor_g2_pow_so4)
            !D1 -> G2    D2 -> pSO4
            call compute_mixed_g2_pow_so4(icount)

         end select

      else                    ! here pure descriptors
         config_desc(icount)%dim_desc = dim_xdesc


         if (.not.(img_weighted)) then
            if (allocated(config_desc(icount)%energy)) deallocate (config_desc(icount)%energy)
            allocate (config_desc(icount)%energy(dim_xdesc, imm))


            if (desc_forces_local) then
               if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
                  if (allocated(config_desc(icount)%force)) deallocate (config_desc(icount)%force)
                  allocate (config_desc(icount)%force(dim_xdesc, i_start_at:i_final_at, 0:imm_neigh, 3))
               end if
            end if
         end if


         select case (descriptor_type)
          case (descriptor_g2)
            call compute_g2(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)

          case (descriptor_g3)
            call compute_g3(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)

          case (descriptor_afs)
            call compute_afs(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)

          case (descriptor_pow_so3)
            call compute_pow_so3(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)

          case (descriptor_pow_so3_3body)
            call compute_pow_so3_3body(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)

          case (descriptor_pow_so4)
            call compute_pow_so4(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)

          case (descriptor_body)
            call compute_body_bonds_and_angles(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)

          case (descriptor_ftnbody)
            call compute_ftnbody(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)

          case (descriptor_zetabody)
            call compute_zetabody(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)
         
          case (descriptor_tbind)
            call compute_tbind(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)

          case (descriptor_ace)
            call compute_ace(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)

          case (descriptor_bispectrum_so4)
            call compute_bispectrum_so4(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)
            if (desc_forces_local) then
               if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
                  config_desc(icount)%force(1:dim_xdesc, i_start_at:i_final_at, 0:imm_neigh, 1:3) = -config_desc(icount)%force(1:dim_xdesc, i_start_at:i_final_at, 0:imm_neigh, 1:3)
               end if
            end if

          case (descriptor_milady)
            if (img_weighted) then
               if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
                  if (allocated(config_desc(icount)%channel)) deallocate (config_desc(icount)%channel)
                  allocate (config_desc(icount)%channel(dim_xdesc, i_start_at:i_final_at))
               end if
            end if
            call compute_r_matrix(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)


          case (descriptor_mtp)
            call compute_mtp(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)

          case default
            if (rangml == 0) write (6, *) 'ML: No implementation for descriptor_type', descriptor_type
            stop "fatal in File: descriptor.f90, Subroutine: compute_descriptor"
         end select


      end if                  ! end compute descriptor type
   end if                  ! desc from read or calculus

   ! complete with kernel if needed
   if (activate_k2b) then
      t00 = MY_MPI_WTIME()
      call allocate_kernel_k2b(icount, i_start_at, i_final_at, dim_kernel_2b, imm, imm_neigh, desc_forces_local)
      call compute_kernel_k2b(i_start_at, i_final_at, icount)
      dim_reduce = imm*dim_kernel_2b
      call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%energy_k2b, dim_reduce, MPI_DOUBLE_PRECISION, MPI_SUM, subworld, codeml)
      t11 = MY_MPI_WTIME()
      time_for_desc_kernel_2b = time_for_desc_kernel_2b + (t11 - t00)
   end if

   if (post_desc) call compute_post_descriptors(icount,  ml_type, zbl_potential, &
                                                     i_start_at, i_final_at, imm, imm_neigh,  &
                                                     desc_forces_local, subworld)

#if (PARA)
   !call MPI_BARRIER(subworld, codeml)
   call my_barrier_subworld(codeml)
   ! all the neighbours
   call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%n_neigh, imm, MPI_INTEGER, MPI_SUM, subworld, codeml)
   call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%kind_neigh, imm*imm_neigh, MPI_INTEGER, MPI_SUM, subworld, codeml)
   if (.not.(img_weighted)) then
      dim_reduce = imm*dim_xdesc
      call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(icount)%energy, dim_reduce, MPI_DOUBLE_PRECISION, MPI_SUM, subworld, codeml)
   end if
   !call MPI_BARRIER(subworld, codeml)
   call my_barrier_subworld(codeml)
#endif

   !$! if (zbl_potential) then
   !$!    call allocate_zbl(icount, i_start_at, i_final_at, desc_forces_local)
   !$!    call compute_zbl(i_start_at, i_final_at, config_real(icount)%n_neigh_zbl, config_real(icount)%kind_neigh_zbl, icount)
   !$!    call energy_force_stress_zbl (icount)
   !$! end if

   call val_renormalize_descriptors(icount)

   ! WRITING DESCRIPTORS
   ! TODO when writing descriptors the force are packed. So computed two times
   if (write_desc .or. write_desc_dump) then
      call write_descriptors(icount)
   end if

#if(PARA)
   if (debug) then
      if (subrank == 0) write (6, '("ML: descriptor was computed in compute_descriptors")')
   end if
   !call MPI_BARRIER(subworld, codeml)
   call my_barrier_subworld(codeml)
#endif

   ! Restore original imm_neigh so that calc_neighbours for the next config
   ! still allocates r_ij, u_ij etc. with the full budget (300).
   imm_neigh = saved_imm_neigh

_MLD_END_ 
end subroutine compute_descriptors


subroutine val_renormalize_descriptors(iconf)
   use mpi
   use mld_mpi
   use module_kind_variables, only: kind_double
   use derived_types, only: config_desc
   use ml_in_ndm_module, only: tmp_val_desc_max, rangml
   use module_chemical_species, only: img_weighted
   use set_limits, only: i_start_on_proc, i_end_on_proc
   use mld_subworld

   implicit none
   integer, intent(in) :: iconf
   integer :: iproc, ia_start, ia_end
   real(kind_double), dimension(:), allocatable :: valproc

   tmp_val_desc_max = -1.d15
   if (img_weighted) then
      if (allocated(valproc)) deallocate(valproc)
      allocate(valproc(subworld_size))
      valproc(:)=0.d0

      do iproc = 0, subworld_size-1
         if (subrank == iproc) then
            ia_start = i_start_on_proc(iproc)
            ia_end = i_end_on_proc(iproc)
            if ((ia_start==0).and.(ia_end==0)) then
               valproc(iproc+1) = 0.d0
            else
               valproc(iproc+1) = MAXVAL(dabs(config_desc(iconf)%channel(:, ia_start:ia_end)))
            end if
         end if
      end do

      call MPI_ALLREDUCE(MPI_IN_PLACE, valproc, subworld_size, MPI_DOUBLE_PRECISION, MPI_SUM, subworld, codeml)
      tmp_val_desc_max = MAXVAL(dabs(valproc))
   else

      tmp_val_desc_max = MAXVAL(dabs(config_desc(iconf)%energy(:, :)))

   end if
end subroutine val_renormalize_descriptors


subroutine init_descriptors()

   !use mpi

   use ml_in_ndm_module, ONLY: debug, &
      descriptor_type, descriptor_g2, descriptor_g3, descriptor_behler, &
      descriptor_afs, descriptor_g2_afs, descriptor_pow_so3, &
      descriptor_pow_so4, descriptor_bispectrum_so4, descriptor_g2_bispectrum_so4, &
      descriptor_mtp, descriptor_body, descriptor_pow_so3_3body, descriptor_ace, descriptor_zetabody, &
      j_max, jj_max,  one_pi, &
      g2_dim, g3_dim, char_desc, descriptor_g2_pow_so4, &
      pow_so4_dim, mtp_dim,  &
      l_max,  &
      weighted, weighted_3ch, &
      descriptor_milady, rmat_dim, dmilady_dim, &
      descriptor_ftnbody, descriptor_zetabody, & 
      descriptor_tbind
   use module_bispectrum_so4, only:  inv_r0, inv_r0_input, bisso4_dim
   use temporary_data_cov, only: dim_xdesc, dim_xdesc1, dim_xdesc2
   use module_body_desc, only: dim_desc_body
   use module_afs, only: afs_dim
   use module_so3, only: n_rbf_so3, pow_so3_dim, radial_pow_so3, radial_bartok, radial_sgg
   use module_ace_desc, only: ace_dim, ace_chem, ACE_CHEM_INCOMPLETE, ACE_CHEM_STANDARD, ACE_CHEM_TS
   use compute_g2_mod, only: init_g2
   use compute_g3_mod, only: init_g3
   use compute_bispectrum_so4_mod, only: gen_dimension_for_bispectrum_so4, test_parameters_bso4
   use compute_afs_mod, only: init_afs_rbf
   use compute_pow_so3_mod, only: init_pow_so3, init_pow_so3_rbf
   use module_chemical_species, only: img_num_ch, img_weighted, fix_no_of_elements, &
      tnn_rdist, tcc01_rdist, tcc02_rdist, tii_rdist
   use module_ftnbody, only : ftnbody_dim, tnn_ftbd, t2b_ftbd, t3b_ftbd, t4b_ftbd, &
      t4b_inner_desc, t4b_inner_deriv, t4b_inner_init, t4b_inner_last,  find_best_length
   use module_kernel_zetabody, only: zetabody_dim   
   use module_tbind, only: tnn_tbind, tbind_dim
   use mld_logger
   use mld_string
   use mld_mpi, only: rangml, comm_mld
   !TORC! use my_mpi_subroutines, only: my_barrier_mld

   implicit none

   _NAMECURRENT_("init_descriptor")

   _MLD_BEGIN_
   if (debug) then
      !TORC! call my_barrier_mld(codeml)
      call comm_mld%barrier()
      if (rangml == 0) write (6, '("ML: descriptor was initialized in init_descriptors...")')
   end if


   select case (descriptor_type)

    case (descriptor_g2)
      call init_g2
      !call gen_param_behler
      !dim_xdesc = g2_dim
      if (weighted) then
         dim_xdesc = 2*g2_dim
         if (weighted_3ch) dim_xdesc = 3*g2_dim
      else
         dim_xdesc = g2_dim
      end if
      char_desc = 'bhg2'

    case (descriptor_g3)
      call init_g3
      !call gen_param_behler
      !if (weighted) then
      !    dim_xdesc = 2*g3_dim
      !else
      dim_xdesc = g3_dim
      !end if
      char_desc = 'bhg3'

    case (descriptor_behler)
      call init_g2
      dim_xdesc1 = g2_dim
      call init_g3
      dim_xdesc2 = g3_dim
      dim_xdesc = g2_dim + g3_dim
      char_desc = 'bhlr'

    case (descriptor_milady)
      if (img_weighted) then
         dim_xdesc = rmat_dim*rmat_dim*img_num_ch
      else
         dim_xdesc = rmat_dim*rmat_dim
      end if
      ! if (weighted) dim_xdesc = dim_xdesc*3
      dmilady_dim = rmat_dim*rmat_dim
      tnn_rdist = 0.d0
      tii_rdist = 0.d0
      tcc01_rdist = 0.d0
      tcc02_rdist = 0.d0
      char_desc = 'rmat'


    case (descriptor_tbind)
      call init_tbind
      dim_xdesc = tbind_dim
      char_desc = 'tbnd'
      tnn_tbind = 0.d0  

    case (descriptor_ftnbody)
      call init_ftnbody
      tnn_ftbd = 0.d0
      t2b_ftbd = 0.d0
      t3b_ftbd = 0.d0
      t4b_ftbd = 0.d0
      t4b_inner_desc = 0.d0
      t4b_inner_deriv = 0.d0
      t4b_inner_init = 0.d0
      t4b_inner_last = 0.d0
      char_desc = 'ftbd'
      dim_xdesc = ftnbody_dim

      if (find_best_length) then
         call init_length_ftnbody()
      end if



    case (descriptor_pow_so3)
      pow_so3_dim = int((1 + l_max))*n_rbf_so3
      if (radial_pow_so3 == radial_bartok) then
         pow_so3_dim = int((1 + l_max))*n_rbf_so3
      end if
      if (radial_pow_so3 == radial_sgg) then
         pow_so3_dim = int((1 + l_max))*(n_rbf_so3 + 1)
      end if
      call init_pow_so3
      call init_pow_so3_rbf()
      call compute_cg_vector(dble(l_max), 1)
      if (weighted) then
         dim_xdesc = 2*pow_so3_dim
         if (weighted_3ch) dim_xdesc = 3*pow_so3_dim
      else
         dim_xdesc = pow_so3_dim
      end if
      char_desc = 'pso3'

    case (descriptor_pow_so3_3body)
      call init_pow_so3_3body
      call init_pow_so3_rbf()
      call compute_cg_vector(dble(l_max), 1)
      if (weighted) then
         dim_xdesc = 2*pow_so3_dim
         if (weighted_3ch) dim_xdesc = 3*pow_so3_dim
      else
         dim_xdesc = pow_so3_dim
      end if
      char_desc = 'psb3'

    case (descriptor_ace)
      call init_ace
      !ACE_CHEM
      if (ace_chem == ACE_CHEM_INCOMPLETE) then
         dim_xdesc = ace_dim
         ! write(*,*) ace_dim 
         ! stop 
      else if (ace_chem == ACE_CHEM_STANDARD) then
         dim_xdesc = ace_dim * fix_no_of_elements
      else if (ace_chem == ACE_CHEM_TS) then
         call log_critical('ML: ACE_CHEM_TS not implemented yet')
         stop 'ML: ACE_CHEM_TS not implemented yet'
      else
         dim_xdesc = ace_dim * fix_no_of_elements
      end if
      char_desc = 'lace'

    case (descriptor_pow_so4)
      jj_max = int(2*j_max)
      inv_r0 = one_pi*inv_r0_input
      call compute_cg_vector(dble(j_max), 2)
      pow_so4_dim = int(2*j_max) + 1
      !dim_xdesc   = int(2*j_max)+1
      if (weighted) then
         dim_xdesc = 2*pow_so4_dim
         if (weighted_3ch) dim_xdesc = 3*pow_so4_dim
      else
         dim_xdesc = pow_so4_dim
      end if
      call init_pow_so4
      char_desc = 'pso4'
      call test_parameters_bso4


    case (descriptor_body)
      call init_body
      dim_xdesc = SUM(dim_desc_body(:))
      char_desc = 'bdys'
      if (weighted) stop 'that descriptor not implemented with weigthed = true'


   case(descriptor_zetabody)
      call init_kernel_zetabody
      dim_xdesc = zetabody_dim
      char_desc = 'zbdy'
      
    case (descriptor_afs)
      call init_afs_rbf()
      if (weighted) then
         dim_xdesc = 2*afs_dim
         if (weighted_3ch) dim_xdesc = 3*afs_dim
      else
         dim_xdesc = afs_dim
      end if
      char_desc = 'afsr'

    case (descriptor_g2_afs)
      call init_g2
      dim_xdesc1 = g2_dim
      call init_afs_rbf()
      dim_xdesc2 = afs_dim
      dim_xdesc = g2_dim + afs_dim
      char_desc = 'g2af'

    case (descriptor_g2_pow_so4)
      call init_g2
      dim_xdesc1 = g2_dim
      jj_max = int(2*j_max)
      inv_r0 = one_pi*inv_r0_input
      call compute_cg_vector(dble(j_max), 2)
      dim_xdesc = g2_dim + int(2*j_max) + 1
      pow_so4_dim = int(2*j_max) + 1
      dim_xdesc2 = int(2*j_max) + 1
      call init_pow_so4
      char_desc = 'g2p4'

    case (descriptor_bispectrum_so4)

      jj_max = int(2*j_max)
      inv_r0 = one_pi*inv_r0_input
      call compute_cg_vector(dble(j_max), 2)
      call gen_dimension_for_bispectrum_so4()
      if (weighted) then
         dim_xdesc = 2*bisso4_dim
         if (weighted_3ch) dim_xdesc = 3*bisso4_dim
      else
         dim_xdesc = bisso4_dim
      end if
      char_desc = 'bso4'
      call test_parameters_bso4

    case (descriptor_mtp)
      call gen_dimension_for_mtp()
      if (weighted) then
         dim_xdesc = mtp_dim
      else
         dim_xdesc = mtp_dim
      end if
      dim_xdesc = mtp_dim
      char_desc = 'mtp3'

    case (descriptor_g2_bispectrum_so4)
      call init_g2
      dim_xdesc1 = g2_dim
      jj_max = int(2*j_max)
      inv_r0 = one_pi*inv_r0_input
      call compute_cg_vector(dble(j_max), 2)
      !snap call gen_dimension_for_bispectrum_so4_all()
      ! Gabor version for which are taken only (J J_1 J_1) componenets
      !old if (lbso4_diag) then
      !old   call gen_dimension_for_bispectrum_so4_diagonal()
      !old else
      !old   call gen_dimension_for_bispectrum_so4_all()
      !old end if
      call gen_dimension_for_bispectrum_so4()
      dim_xdesc2 = bisso4_dim
      !if (weighted) then
      !    dim_xdesc=2*bisso4_dim+2*g2_dim
      !else
      dim_xdesc = bisso4_dim + g2_dim
      !end if
      if (rangml == 0) write (6, '("ML: hybrid descritor in init desc  bso4 + g2:  ", 2i5)') bisso4_dim, g2_dim
      char_desc = 'g2b4'
      call test_parameters_bso4
    case default
      if (rangml == 0) write (6, *) 'No implementation for descriptor_type...', descriptor_type
      stop "fatal in File: descriptor.f90, Subroutine: init_descriptor"
   end select

   call log_info('descriptor '//char_desc//' has the dimension '//vtoa(dim_xdesc))

   if (debug) then
      !TORC! call my_barrier_mld(codeml)
      call comm_mld%barrier()
      if (rangml == 0) write (6, '("ML: descriptor was initialized in init_descriptors...")')
   end if

   !TORC! call my_barrier_mld(codeml)
   call comm_mld%barrier()

   _MLD_END_
end subroutine init_descriptors


subroutine init_write_descriptors()

   use ml_in_ndm_module, ONLY: write_desc, write_desc_dump, write_design_matrix, write_test_design_matrix, & 
                               desc_file_format, hdf_type, eml_type, csv_type, npz_type
   use module_db_setup, only: db_path
   use mld_logger
   use mld_mpi, only: mld_rank, comm_mld
   use mld_subworld, only: subrank
#if(MLD_HDF5)
   use mld_hdf5
#endif
   implicit none
#if(MLD_HDF5)       
   integer :: hdf_error
   character(len=80), parameter :: fname = 'all_descriptors.h5'
#endif   

   _NAMECURRENT_("init_write_descriptors")
   _MLD_BEGIN_

   if (write_desc .or. write_desc_dump) then
      call log_info('ML: descriptors will be written')
   else
      call log_info('ML: descriptors will NOT be written')
   end if

   if (write_design_matrix) then
      call log_info('ML: design matrix will be written')
   end if 
   if (write_test_design_matrix) then
      call log_info('ML: test design matrix will be written')
   end if
   if (.not.(write_desc .or. write_desc_dump)) then 
   if (desc_file_format == hdf_type) then
      call log_info('ML: descriptors will be written in HDF5 format')
#if(MLD_HDF5)
      ! Check if HDF5 is already initialized
      call h5open_f(hdf_error)
      if (hdf_error /= 0) then
          call log_error('ML: Failed to initialize HDF5')
          stop 'ML: HDF5 initialization failed'
      end if      
      call log_info('ML: HDF5 successfully initialized')
#else
      call log_critical('ML: HDF5 format requested but HDF5 support not compiled in')
      stop 'ML: HDF5 support not available'
#endif
   else if (desc_file_format == eml_type) then
      call log_info('ML: descriptors will be written in EML format')
   else if (desc_file_format == csv_type) then
      call log_info('ML: descriptors will be written in CSV format')
   else if(desc_file_format == npz_type) then
      call log_info('ML: descriptors will be written in NPZ format')
      call log_warning('ML: NPZ format not working for more than 1 processor, use only for debugging')      
   else
      call log_critical('ML: desc_file_format not recognized in init_write_descriptors')
      stop 'ML: desc_file_format not recognized in init_write_descriptors'
   end if
   end if ! not write_desc or write_desc_dump

#if(MLD_HDF5)
    if (mld_rank == 0 ) then 
   !if (subrank == 0 ) then 
     if (desc_file_format == hdf_type) then
        hFile%fname = 'desc'//trim(adjustl(db_path))//trim(fname)
        call hFile%hdf5_create_file(hFile%fname)        ! H5F_ACC_TRUNC_F inside
     end if
     hFile%file_open = 1
   end if
   !call comm_mld%barrier

#endif   

   _MLD_END_
end subroutine init_write_descriptors

#if(MLD_HDF5)
subroutine end_write_descriptors
   use mld_hdf5
   use mld_mpi, only: mld_rank
   use mld_subworld, only: subrank
   use ml_in_ndm_module, only: desc_file_format, hdf_type

   if (desc_file_format == hdf_type) then 
     if (mld_rank == 0) then 
     !if (subrank == 0) then 
       ! 5. close the file after the last configuration
       call hFile%hdf5_close_file()
       hFile%file_open = 0
     end if 
   else 
     call log_info('ML: end writing descriptors')  
   end if 
end subroutine end_write_descriptors
#else
subroutine end_write_descriptors
   use mld_logger
   use module_db_setup, only: db_path
   implicit none

   _NAMECURRENT_("end_write_descriptors")
   _MLD_BEGIN_

   call log_info('ML: end writing descriptors')

   _MLD_END_
end subroutine end_write_descriptors
#endif


subroutine build_database_with_function(nd_local_data, dim_local_xdesc, yfunc, xdesc)
   use ml_in_ndm_module, ONLY: toy_model, seed, rangml
   use toy_models

   implicit none

   integer, intent(in)  :: nd_local_data, dim_local_xdesc
   !real(kind=kind(1.d0)), dimension(:,:), intent(out) ::  xdesc(dim_local_xdesc,nd_local_data)
   real(kind=kind(1.d0)), dimension(:, :), intent(out)      :: xdesc
   real(kind=kind(1.d0)), intent(out)     :: yfunc(nd_local_data)
   !local
   real(kind=kind(1.d0))      :: rvalue, leng, y, x(dim_local_xdesc), y_err(nd_local_data)
   integer  :: i, j


   ! internal length for data ... to see the units effetct.
   leng = 3.0

   if (toy_model) then
      ! generate nd_local_data random numbers between (0,L=leng)
      call random_seed(seed)
      do i = 1, nd_local_data
         do j = 1, dim_local_xdesc
            call random_number(rvalue)
            xdesc(j, i) = rvalue*leng
         end do
      end do
      do i = 1, nd_local_data
         x(:) = real(xdesc(:, i))
         call toy_nD(x, y, leng, dim_local_xdesc)
         yfunc(i) = y
      end do
      call generate_random_gaussian(y_err, nd_local_data)
      yfunc(:) = yfunc(:) + 0.02d0*y_err(:)
   else
      if (rangml == 0) write (*, *) 'ML error: <build_database_with_function> the toy_model=.F. in this subroutine is not yet implemented'
      stop 'stop in build_database_with_function'
   end if

end subroutine build_database_with_function



subroutine train_deallocate_desc(iconf)

   use ml_in_ndm_module, only: ml_type_descriptors, ml_type_analysis, ml_type, &
      mld_order, mld_kernel, mld_regularization_type, mld_regularization_type_home
   use derived_types, only: config_desc
   use module_kernel_2b, only: activate_k2b
   use module_chemical_species, only: img_weighted
   implicit none
   integer, intent(in)  :: iconf


   if (img_weighted) then
      if (allocated(config_desc(iconf)%channel)) deallocate (config_desc(iconf)%channel)
   end if
   if (ml_type == ml_type_analysis) then
      if (allocated(config_desc(iconf)%force)) deallocate (config_desc(iconf)%force)
      if (allocated(config_desc(iconf)%n_neigh)) deallocate (config_desc(iconf)%n_neigh)
      if (allocated(config_desc(iconf)%kind_neigh)) deallocate (config_desc(iconf)%kind_neigh)
   else
      !if (mld_order == mld_kernel) then
      !  if (allocated(config_desc(iconf)%force_kernel)) deallocate (config_desc(iconf)%force_kernel)
      !end if
      !TODOAcc

      if (.not.(mld_regularization_type == mld_regularization_type_home)) then

         if (allocated(config_desc(iconf)%n_neigh_ghost)) deallocate (config_desc(iconf)%n_neigh_ghost)
         if (allocated(config_desc(iconf)%incell)) deallocate (config_desc(iconf)%incell)
         if (allocated(config_desc(iconf)%kind_neigh)) deallocate (config_desc(iconf)%kind_neigh)
         if (allocated(config_desc(iconf)%kind_neigh_ghost)) deallocate (config_desc(iconf)%kind_neigh_ghost)
         if (allocated(config_desc(iconf)%kind_neigh_proc)) deallocate (config_desc(iconf)%kind_neigh_proc)

         if (allocated(config_desc(iconf)%n_neigh)) deallocate (config_desc(iconf)%n_neigh)
         if (allocated(config_desc(iconf)%energy)) deallocate (config_desc(iconf)%energy)
         if (allocated(config_desc(iconf)%energy_kernel)) deallocate (config_desc(iconf)%energy_kernel)
         if (allocated(config_desc(iconf)%energy_k2b)) deallocate (config_desc(iconf)%energy_k2b)
         if (allocated(config_desc(iconf)%channel)) deallocate (config_desc(iconf)%channel)
         if (allocated(config_desc(iconf)%energy1)) deallocate (config_desc(iconf)%energy1)
         if (allocated(config_desc(iconf)%energy2)) deallocate (config_desc(iconf)%energy2)

         if (allocated(config_desc(iconf)%pack_energy_linear)) deallocate (config_desc(iconf)%pack_energy_linear)
         if (allocated(config_desc(iconf)%pack_energy_quadratic)) deallocate (config_desc(iconf)%pack_energy_quadratic)
         if (allocated(config_desc(iconf)%pack_energy_polyc)) deallocate (config_desc(iconf)%pack_energy_polyc)
         if (allocated(config_desc(iconf)%pack_energy_kernel)) deallocate (config_desc(iconf)%pack_energy_kernel)
         if (allocated(config_desc(iconf)%pack_energy_k2b)) deallocate (config_desc(iconf)%pack_energy_k2b)

         if (allocated(config_desc(iconf)%force)) deallocate (config_desc(iconf)%force)
         if (allocated(config_desc(iconf)%force_kernel)) deallocate (config_desc(iconf)%force_kernel)
         if (allocated(config_desc(iconf)%force_k2b)) deallocate (config_desc(iconf)%force_k2b)
         if (allocated(config_desc(iconf)%force1)) deallocate (config_desc(iconf)%force1)
         if (allocated(config_desc(iconf)%force2)) deallocate (config_desc(iconf)%force2)


         if (allocated(config_desc(iconf)%pack_force_linear)) deallocate (config_desc(iconf)%pack_force_linear)
         if (allocated(config_desc(iconf)%pack_force_quadratic)) deallocate (config_desc(iconf)%pack_force_quadratic)
         if (allocated(config_desc(iconf)%pack_force_polyc)) deallocate (config_desc(iconf)%pack_force_polyc)
         if (allocated(config_desc(iconf)%pack_force_kernel)) deallocate (config_desc(iconf)%pack_force_kernel)
         if (allocated(config_desc(iconf)%pack_force_k2b)) deallocate (config_desc(iconf)%pack_force_k2b)

         if (allocated(config_desc(iconf)%pack_stress_linear)) deallocate (config_desc(iconf)%pack_stress_linear)
         if (allocated(config_desc(iconf)%pack_stress_quadratic)) deallocate (config_desc(iconf)%pack_stress_quadratic)
         if (allocated(config_desc(iconf)%pack_stress_polyc)) deallocate (config_desc(iconf)%pack_stress_polyc)
         if (allocated(config_desc(iconf)%pack_stress_kernel)) deallocate (config_desc(iconf)%pack_stress_kernel)
         if (allocated(config_desc(iconf)%pack_stress_k2b)) deallocate (config_desc(iconf)%pack_stress_k2b)
      end if

   end if


   if ((ml_type == ml_type_descriptors) .or. (ml_type == ml_type_analysis)) then
      if (allocated(config_desc(iconf)%type_neigh)) deallocate (config_desc(iconf)%type_neigh)
      if (allocated(config_desc(iconf)%n_neigh_ghost)) deallocate (config_desc(iconf)%n_neigh_ghost)
      if (allocated(config_desc(iconf)%kind_neigh_ghost)) deallocate (config_desc(iconf)%kind_neigh_ghost)
      if (allocated(config_desc(iconf)%kind_neigh_proc)) deallocate (config_desc(iconf)%kind_neigh_proc)
      if (allocated(config_desc(iconf)%pack_force_linear)) deallocate (config_desc(iconf)%pack_force_linear)
      if (allocated(config_desc(iconf)%pack_energy_linear)) deallocate (config_desc(iconf)%pack_energy_linear)
      if (allocated(config_desc(iconf)%pack_stress_linear)) deallocate (config_desc(iconf)%pack_stress_linear)
      if (mld_order == mld_kernel) then
         if (allocated(config_desc(iconf)%force_kernel)) deallocate (config_desc(iconf)%force_kernel)
      end if
      if (activate_k2b) then
         if (allocated(config_desc(iconf)%force)) deallocate (config_desc(iconf)%force)
         if (allocated(config_desc(iconf)%force1)) deallocate (config_desc(iconf)%force1)
         if (allocated(config_desc(iconf)%force2)) deallocate (config_desc(iconf)%force2)


         if (allocated(config_desc(iconf)%energy_k2b)) deallocate (config_desc(iconf)%energy_k2b)
         if (allocated(config_desc(iconf)%pack_energy_k2b)) deallocate (config_desc(iconf)%pack_energy_k2b)
         if (allocated(config_desc(iconf)%force_k2b)) deallocate (config_desc(iconf)%force_k2b)
         if (allocated(config_desc(iconf)%pack_stress_k2b)) deallocate (config_desc(iconf)%pack_stress_k2b)
      end if
   end if
   if (allocated(config_desc(iconf)%amask)) deallocate(config_desc(iconf)%amask)

end subroutine train_deallocate_desc


subroutine train_deallocate_real(iconf)

   use ml_in_ndm_module, only: ml_type_descriptors, ml_type_analysis, ml_type, &
      mld_regularization_type, mld_regularization_type_home
   use derived_types, only: config_real
   implicit none
   integer, intent(in)  :: iconf



   if (ml_type == ml_type_analysis) then
      if (allocated(config_real(iconf)%pos_cart)) deallocate (config_real(iconf)%pos_cart)
      if (allocated(config_real(iconf)%pos_crst)) deallocate (config_real(iconf)%pos_crst)
      if (allocated(config_real(iconf)%force)) deallocate (config_real(iconf)%force)
      if (allocated(config_real(iconf)%atomic_spin)) deallocate (config_real(iconf)%atomic_spin)
   else
      !if (mld_order == mld_kernel) then
      !  if (allocated(config_desc(iconf)%force_kernel)) deallocate (config_desc(iconf)%force_kernel)
      !end if
      !TODOAcc

      if (.not.(mld_regularization_type == mld_regularization_type_home)) then
         if (allocated(config_real(iconf)%pos_cart)) deallocate (config_real(iconf)%pos_cart)
         if (allocated(config_real(iconf)%pos_crst)) deallocate (config_real(iconf)%pos_crst)
         if (allocated(config_real(iconf)%force)) deallocate (config_real(iconf)%force)
         if (allocated(config_real(iconf)%atomic_spin)) deallocate (config_real(iconf)%atomic_spin)
      end if

   end if


   if ((ml_type == ml_type_descriptors) .or. (ml_type == ml_type_analysis)) then
      if (allocated(config_real(iconf)%pos_cart)) deallocate (config_real(iconf)%pos_cart)
      if (allocated(config_real(iconf)%pos_crst)) deallocate (config_real(iconf)%pos_crst)
      if (allocated(config_real(iconf)%force)) deallocate (config_real(iconf)%force)
      if (allocated(config_real(iconf)%atomic_spin)) deallocate (config_real(iconf)%atomic_spin)
   end if


end subroutine train_deallocate_real



subroutine status_allocate_config_desc (iconf)
   use ml_in_ndm_module, only: rangml
   use derived_types, only: config_desc
   integer, intent(in) :: iconf

   write(6, '("|desc--------------------------->", i5)') iconf

   write(6, '("n_neigh_ghost ..................: ", i5, l2)') rangml, allocated(config_desc(iconf)%n_neigh_ghost)
   write(6, '("incell .........................: ", i5, l2)') rangml, allocated(config_desc(iconf)%incell)
   write(6, '("kind_neigh_ghost ...............: ", i5, l2)') rangml, allocated(config_desc(iconf)%kind_neigh_ghost)
   write(6, '("kind_neigh_proc ................: ", i5, l2)') rangml, allocated(config_desc(iconf)%kind_neigh_proc)
   write(6, '("stat_dist_mcd ..................: ", i5, l2)') rangml, allocated(config_desc(iconf)%stat_dist_mcd)
   write(6, '("stat_dist_maha .................: ", i5, l2)') rangml, allocated(config_desc(iconf)%stat_dist_maha)
   write(6, '("stat_norm ......................: ", i5, l2)') rangml, allocated(config_desc(iconf)%stat_norm)
   write(6, '("stat_energy ....................: ", i5, l2)') rangml, allocated(config_desc(iconf)%stat_energy)
   write(6, '("stat_norm_mean .................: ", i5, l2)') rangml, allocated(config_desc(iconf)%stat_norm_mean)


   write(6, '("n_neigh ........................: ", i5, l2)') rangml, allocated(config_desc(iconf)%n_neigh)
   write(6, '("energy .........................: ", i5, l2)') rangml, allocated(config_desc(iconf)%energy)
   write(6, '("energy_kernel ..................: ", i5, l2)') rangml, allocated(config_desc(iconf)%energy_kernel)
   write(6, '("energy_k2b .....................: ", i5, l2)') rangml, allocated(config_desc(iconf)%energy_k2b)
   write(6, '("channel ........................: ", i5, l2)') rangml, allocated(config_desc(iconf)%channel)
   write(6, '("energy1 ........................: ", i5, l2)') rangml, allocated(config_desc(iconf)%energy1)
   write(6, '("energy2 ........................: ", i5, l2)') rangml, allocated(config_desc(iconf)%energy2)

   write(6, '("force ..........................: ", i5, l2)') rangml, allocated(config_desc(iconf)%force)
   write(6, '("force_kernel ...................: ", i5, l2)') rangml, allocated(config_desc(iconf)%force_kernel)
   write(6, '("force_k2b ......................: ", i5, l2)') rangml, allocated(config_desc(iconf)%force_k2b)
   write(6, '("force1 .........................: ", i5, l2)') rangml, allocated(config_desc(iconf)%force1)
   write(6, '("force2 .........................: ", i5, l2)') rangml, allocated(config_desc(iconf)%force2)

   write(6, '("pack_energy_linear .............: ", i5, l2)') rangml, allocated(config_desc(iconf)%pack_energy_linear)
   write(6, '("pack_energy_quadratic ..........: ", i5, l2)') rangml, allocated(config_desc(iconf)%pack_energy_quadratic)
   write(6, '("pack_energy_polyc ..............: ", i5, l2)') rangml, allocated(config_desc(iconf)%pack_energy_polyc)
   write(6, '("pack_energy_kernel .............: ", i5, l2)') rangml, allocated(config_desc(iconf)%pack_energy_kernel)
   write(6, '("pack_energy_k2b ................: ", i5, l2)') rangml, allocated(config_desc(iconf)%pack_energy_k2b)


   write(6, '("pack_force_linear ..............: ", i5, l2)') rangml, allocated(config_desc(iconf)%pack_force_linear)
   write(6, '("pack_force_quadratic ...........: ", i5, l2)') rangml, allocated(config_desc(iconf)%pack_force_quadratic)
   write(6, '("pack_force_polyc ...............: ", i5, l2)') rangml, allocated(config_desc(iconf)%pack_force_polyc)
   write(6, '("pack_force_kernel ..............: ", i5, l2)') rangml, allocated(config_desc(iconf)%pack_force_kernel)
   write(6, '("pack_force_k2b .................: ", i5, l2)') rangml, allocated(config_desc(iconf)%pack_force_k2b)

   write(6, '("pack_stress_linear .............: ", i5, l2)') rangml, allocated(config_desc(iconf)%pack_stress_linear)
   write(6, '("pack_stress_quadratic ..........: ", i5, l2)') rangml, allocated(config_desc(iconf)%pack_stress_quadratic)
   write(6, '("pack_stress_kernel .............: ", i5, l2)') rangml, allocated(config_desc(iconf)%pack_stress_kernel)
   write(6, '("pack_stress_polyc ..............: ", i5, l2)') rangml, allocated(config_desc(iconf)%pack_stress_polyc)
   write(6, '("pack_stress_k2b ................: ", i5, l2)') rangml, allocated(config_desc(iconf)%pack_stress_k2b)

end subroutine status_allocate_config_desc

subroutine status_allocate_config_real (iconf)
   use ml_in_ndm_module, only: rangml
   use derived_types, only: config_real
   integer, intent(in) :: iconf

   write(6, '("|real-------------------------->", i5)') iconf

   write(6, '("n_neigh ........................: ", i5, l2)') rangml, allocated(config_real(iconf)%n_neigh)
   write(6, '("n_neigh_zbl.....................: ", i5, l2)') rangml, allocated(config_real(iconf)%n_neigh_zbl)
   write(6, '("type_neigh .....................: ", i5, l2)') rangml, allocated(config_real(iconf)%type_neigh)
   write(6, '("kind_neigh .....................: ", i5, l2)') rangml, allocated(config_real(iconf)%kind_neigh)
   write(6, '("kind_neigh_zbl .................: ", i5, l2)') rangml, allocated(config_real(iconf)%kind_neigh_zbl)
   write(6, '("incell..........................: ", i5, l2)') rangml, allocated(config_real(iconf)%incell)
   write(6, '("proc_atom.......................: ", i5, l2)') rangml, allocated(config_real(iconf)%proc_atom)
   write(6, '("r_ij ...........................: ", i5, l2)') rangml, allocated(config_real(iconf)%r_ij)
   write(6, '("u_ij ...........................: ", i5, l2)') rangml, allocated(config_real(iconf)%u_ij)
   write(6, '("u_per ..........................: ", i5, l2)') rangml, allocated(config_real(iconf)%u_per)
   write(6, '("u_at ...........................: ", i5, l2)') rangml, allocated(config_real(iconf)%u_at)
   write(6, '("itype ..........................: ", i5, l2)') rangml, allocated(config_real(iconf)%itype)
   write(6, '("itype_db .......................: ", i5, l2)') rangml, allocated(config_real(iconf)%itype_db)
   write(6, '("itype_to_global ................: ", i5, l2)') rangml, allocated(config_real(iconf)%itype_to_global)
   write(6, '("fix_type_poscar_to_periodic ....: ", i5, l2)') rangml, allocated(config_real(iconf)%fix_type_poscar_to_periodic)
   write(6, '("pos_cart........................: ", i5, l2)') rangml, allocated(config_real(iconf)%pos_cart)
   write(6, '("pos_crst........................: ", i5, l2)') rangml, allocated(config_real(iconf)%pos_crst)
   write(6, '("force...........................: ", i5, l2)') rangml, allocated(config_real(iconf)%force)
   write(6, '("atomic_spin.....................: ", i5, l2)') rangml, allocated(config_real(iconf)%atomic_spin)
   write(6, '("ref_energy_per_element .........: ", i5, l2)') rangml, allocated(config_real(iconf)%ref_energy_per_element)
   write(6, '("mass_per_type ..................: ", i5, l2)') rangml, allocated(config_real(iconf)%mass_per_type)
   write(6, '("Z_per_type .....................: ", i5, l2)') rangml, allocated(config_real(iconf)%Z_per_type)
   write(6, '("covalent_radius_per_type .......: ", i5, l2)') rangml, allocated(config_real(iconf)%covalent_radius_per_type)
   write(6, '("weight_per_type ................: ", i5, l2)') rangml, allocated(config_real(iconf)%weight_per_type)
   write(6, '("weight_per_type_3ch ............: ", i5, l2)') rangml, allocated(config_real(iconf)%weight_per_type_3ch)
   write(6, '("invisible_per_type .............: ", i5, l2)') rangml, allocated(config_real(iconf)%invisible_per_type)
   write(6, '("wspecies_per_type_ch ...........: ", i5, l2)') rangml, allocated(config_real(iconf)%wspecies_per_type_ch)
   write(6, '("elocal_zbl .....................: ", i5, l2)') rangml, allocated(config_real(iconf)%elocal_zbl)
   write(6, '("flocal_zbl .....................: ", i5, l2)') rangml, allocated(config_real(iconf)%flocal_zbl)
   write(6, '("slocal_zbl .....................: ", i5, l2)') rangml, allocated(config_real(iconf)%slocal_zbl)
   write(6, '("szbl ...........................: ", i5, l2)') rangml, allocated(config_real(iconf)%szbl)
   write(6, '("fzbl ...........................: ", i5, l2)') rangml, allocated(config_real(iconf)%fzbl)


end subroutine status_allocate_config_real
