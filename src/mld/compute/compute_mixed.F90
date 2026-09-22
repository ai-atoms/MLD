! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   MiLaDy (Machine Learning Dynamics)
! HND X
! HND X
! HND X   Milady was written and designed by Mihai-Cosmin Marinica, Alexandra M. Goryaeva
! HND X   Copyright 2015-2024.
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

subroutine compute_mixed_behler(icount)


   use mld_mpi, only: codeml 
#ifdef MLD_NDM
   use gen_com_m_ml, ONLY: imm
#else
   use ondm_gen_com_m, ONLY: imm
#endif
   use ml_in_ndm_module, only: desc_forces, imm_neigh, i_start_at, i_final_at
   use derived_types, only: config_desc, config_real
   use temporary_data_cov, only: dim_xdesc, dim_xdesc1, dim_xdesc2

   use compute_g2_mod, only: compute_g2
   use compute_g3_mod, only: compute_g3
   !use mld_subworld, only: subworld
   use my_mpi_subroutines, only: my_barrier_subworld

   implicit none

   integer, intent(in)  :: icount
   logical  :: desc_forces_local


   desc_forces_local = desc_forces .and. (config_real(icount)%has_force .or. config_real(icount)%has_stress)
   ! Computing the first kind of desc ...
   if (allocated(config_desc(icount)%energy)) deallocate (config_desc(icount)%energy); allocate (config_desc(icount)%energy(dim_xdesc1, imm))
   if (desc_forces_local) then
      if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
         if (allocated(config_desc(icount)%force)) deallocate (config_desc(icount)%force); allocate (config_desc(icount)%force(dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 3))
      end if
   end if

   dim_xdesc = config_desc(icount)%dim_desc1
   config_desc(icount)%dim_desc = dim_xdesc
   call compute_g2(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)
#if(PARA)
   !call MPI_BARRIER(subworld, codeml)
   call my_barrier_subworld(codeml)
#endif
   if (allocated(config_desc(icount)%energy1)) deallocate (config_desc(icount)%energy1); allocate (config_desc(icount)%energy1(dim_xdesc1, imm))
   config_desc(icount)%energy1(:, :) = 0.d0
   config_desc(icount)%energy1(1:dim_xdesc1, 1:imm) = config_desc(icount)%energy(1:dim_xdesc1, 1:imm)
   if (desc_forces_local) then
      if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
         if (allocated(config_desc(icount)%force1)) deallocate (config_desc(icount)%force1); allocate (config_desc(icount)%force1(dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 3))
         config_desc(icount)%force1(:, :, :, :) = 0.d0
         config_desc(icount)%force1(1:dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 1:3) = config_desc(icount)%force(1:dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 1:3)
      end if
   end if

   ! Computing the second kind of desc ...
   if (allocated(config_desc(icount)%energy)) deallocate (config_desc(icount)%energy); allocate (config_desc(icount)%energy(dim_xdesc2, imm))
   if (desc_forces_local) then
      if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
         if (allocated(config_desc(icount)%force)) deallocate (config_desc(icount)%force); allocate (config_desc(icount)%force(dim_xdesc2, i_start_at:i_final_at, 0:imm_neigh, 3))
      end if
   end if

   dim_xdesc = config_desc(icount)%dim_desc2
   config_desc(icount)%dim_desc = dim_xdesc
   call compute_g3(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)
#if(PARA)
   !call MPI_BARRIER(subworld, codeml)
   call my_barrier_subworld(codeml)
#endif

   if (allocated(config_desc(icount)%energy2)) deallocate (config_desc(icount)%energy2); allocate (config_desc(icount)%energy2(dim_xdesc2, imm))
   config_desc(icount)%energy2(:, :) = 0.d0
   config_desc(icount)%energy2(1:dim_xdesc2, 1:imm) = config_desc(icount)%energy(1:dim_xdesc2, 1:imm)
   if (desc_forces_local) then
      if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
         if (allocated(config_desc(icount)%force2)) deallocate (config_desc(icount)%force2); allocate (config_desc(icount)%force2(dim_xdesc2, i_start_at:i_final_at, 0:imm_neigh, 3))
         config_desc(icount)%force2(:, :, :, :) = 0.d0
         config_desc(icount)%force2(1:dim_xdesc2, i_start_at:i_final_at, 0:imm_neigh, 1:3) = config_desc(icount)%force(1:dim_xdesc2, i_start_at:i_final_at, 0:imm_neigh, 1:3)
      end if
   end if

   ! Fill the final %energy and %force.
   dim_xdesc = config_desc(icount)%dim_desc1 + config_desc(icount)%dim_desc2
   config_desc(icount)%dim_desc = dim_xdesc
   if (allocated(config_desc(icount)%energy)) deallocate (config_desc(icount)%energy); allocate (config_desc(icount)%energy(dim_xdesc, imm))
   config_desc(icount)%energy(:, :) = 0.d0
   config_desc(icount)%energy(1:dim_xdesc1, 1:imm) = config_desc(icount)%energy1(1:dim_xdesc1, 1:imm)
   config_desc(icount)%energy(1 + dim_xdesc1:dim_xdesc, 1:imm) = config_desc(icount)%energy2(1:dim_xdesc2, 1:imm)
   deallocate (config_desc(icount)%energy1)
   deallocate (config_desc(icount)%energy2)

   if (desc_forces_local) then
      if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
         if (allocated(config_desc(icount)%force)) deallocate (config_desc(icount)%force); allocate (config_desc(icount)%force(dim_xdesc, i_start_at:i_final_at, 0:imm_neigh, 3))
         config_desc(icount)%force(:, :, :, :) = 0.d0
         config_desc(icount)%force(1:dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 1:3) = config_desc(icount)%force1(1:dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 1:3)
         config_desc(icount)%force(1 + dim_xdesc1:dim_xdesc, i_start_at:i_final_at, 0:imm_neigh, 1:3) = config_desc(icount)%force2(1:dim_xdesc2, i_start_at:i_final_at, 0:imm_neigh, 1:3)
         deallocate (config_desc(icount)%force1)
         deallocate (config_desc(icount)%force2)
      end if
   end if

end subroutine compute_mixed_behler



subroutine compute_mixed_g2_bispectrum_so4(icount)


   use mld_mpi, only: codeml  
#ifdef MLD_NDM
   use gen_com_m_ml, ONLY: imm
#else
   use ondm_gen_com_m, ONLY: imm
#endif
   use ml_in_ndm_module, only: desc_forces, imm_neigh, i_start_at, i_final_at
   use derived_types, only: config_desc, config_real
   use temporary_data_cov, only: dim_xdesc, dim_xdesc1, dim_xdesc2
   !use mld_subworld, only: subworld
   use compute_g2_mod, only: compute_g2
   use compute_bispectrum_so4_mod, only: compute_bispectrum_so4
   use my_mpi_subroutines, only: my_barrier_subworld
   implicit none

   integer, intent(in)  :: icount
   logical  :: desc_forces_local

   desc_forces_local = desc_forces .and. (config_real(icount)%has_force .or. config_real(icount)%has_stress)

   ! Computing the first kind of desc
   if (allocated(config_desc(icount)%energy)) deallocate (config_desc(icount)%energy)
   allocate (config_desc(icount)%energy(dim_xdesc1, imm))
   if (desc_forces_local) then
      if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
         if (allocated(config_desc(icount)%force)) deallocate (config_desc(icount)%force)
         allocate (config_desc(icount)%force(dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 3))
      end if
   end if

   dim_xdesc = config_desc(icount)%dim_desc1
   config_desc(icount)%dim_desc = dim_xdesc
   call compute_g2(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)
#if(PARA)
   !call MPI_BARRIER(subworld, codeml)
   call my_barrier_subworld(codeml)
#endif

   if (allocated(config_desc(icount)%energy1)) deallocate (config_desc(icount)%energy1)
   allocate (config_desc(icount)%energy1(dim_xdesc1, imm))
   config_desc(icount)%energy1(:, :) = 0.d0
   config_desc(icount)%energy1(1:dim_xdesc1, 1:imm) = config_desc(icount)%energy(1:dim_xdesc1, 1:imm)
   if (desc_forces_local) then
      if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
         if (allocated(config_desc(icount)%force1)) deallocate (config_desc(icount)%force1)
         allocate (config_desc(icount)%force1(dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 3))
         config_desc(icount)%force1(:, :, :, :) = 0.d0
         config_desc(icount)%force1(1:dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 1:3) = config_desc(icount)%force(1:dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 1:3)
      end if
   end if

   ! Computing the second kind of desc ...
   if (allocated(config_desc(icount)%energy)) deallocate (config_desc(icount)%energy)
   allocate (config_desc(icount)%energy(dim_xdesc2, imm))
   if (desc_forces_local) then
      if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
         if (allocated(config_desc(icount)%force)) deallocate (config_desc(icount)%force)
         allocate (config_desc(icount)%force(dim_xdesc2, i_start_at:i_final_at, 0:imm_neigh, 3))
      end if
   end if

   dim_xdesc = config_desc(icount)%dim_desc2
   config_desc(icount)%dim_desc = dim_xdesc
   call compute_bispectrum_so4(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)
#if(PARA)
   !call MPI_BARRIER(subworld, codeml)
   call my_barrier_subworld(codeml)
#endif

   if (allocated(config_desc(icount)%energy2)) deallocate (config_desc(icount)%energy2)
   allocate (config_desc(icount)%energy2(dim_xdesc2, imm))
   config_desc(icount)%energy2(:, :) = 0.d0
   config_desc(icount)%energy2(1:dim_xdesc2, 1:imm) = config_desc(icount)%energy(1:dim_xdesc2, 1:imm)
   if (desc_forces_local) then
      if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
         if (allocated(config_desc(icount)%force2)) deallocate (config_desc(icount)%force2)
         allocate (config_desc(icount)%force2(dim_xdesc2, i_start_at:i_final_at, 0:imm_neigh, 3))
         config_desc(icount)%force2(:, :, :, :) = 0.d0
         config_desc(icount)%force2(1:dim_xdesc2, i_start_at:i_final_at, 0:imm_neigh, 1:3) = config_desc(icount)%force(1:dim_xdesc2, i_start_at:i_final_at, 0:imm_neigh, 1:3)
      end if
   end if

   ! Fill the final %energy and %force.
   dim_xdesc = config_desc(icount)%dim_desc1 + config_desc(icount)%dim_desc2
   config_desc(icount)%dim_desc = dim_xdesc
   if (allocated(config_desc(icount)%energy)) deallocate (config_desc(icount)%energy)
   allocate (config_desc(icount)%energy(dim_xdesc, imm))
   config_desc(icount)%energy(:, :) = 0.d0
   config_desc(icount)%energy(1:dim_xdesc1, 1:imm) = config_desc(icount)%energy1(1:dim_xdesc1, 1:imm)
   config_desc(icount)%energy(1 + dim_xdesc1:dim_xdesc, 1:imm) = config_desc(icount)%energy2(1:dim_xdesc2, 1:imm)
   deallocate (config_desc(icount)%energy1)
   deallocate (config_desc(icount)%energy2)

   if (desc_forces_local) then
      if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
         if (allocated(config_desc(icount)%force)) deallocate (config_desc(icount)%force)
         allocate (config_desc(icount)%force(dim_xdesc, i_start_at:i_final_at, 0:imm_neigh, 3))
         config_desc(icount)%force(:, :, :, :) = 0.d0
         config_desc(icount)%force(1:dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 1:3) = config_desc(icount)%force1(1:dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 1:3)
         config_desc(icount)%force(1 + dim_xdesc1:dim_xdesc, i_start_at:i_final_at, 0:imm_neigh, 1:3) = -config_desc(icount)%force2(1:dim_xdesc2, i_start_at:i_final_at, 0:imm_neigh, 1:3)
         deallocate (config_desc(icount)%force1)
         deallocate (config_desc(icount)%force2)
      end if
   end if

end subroutine compute_mixed_g2_bispectrum_so4



subroutine compute_mixed_g2_pow_so4(icount)
#if(PARA)
   !use mpi
   use mld_mpi, only: codeml  
#endif
#ifdef MLD_NDM
   use gen_com_m_ml, ONLY: imm
#else
   use ondm_gen_com_m, ONLY: imm
#endif
   use ml_in_ndm_module, only: desc_forces, imm_neigh, i_start_at, i_final_at
   use derived_types, only: config_desc, config_real
   use temporary_data_cov, only: dim_xdesc, dim_xdesc1, dim_xdesc2
   use compute_g2_mod, only: compute_g2
   !use mld_subworld, only: subworld
   use my_mpi_subroutines, only: my_barrier_subworld
   use module_compute_pow_so4, only: compute_pow_so4
   implicit none

   integer, intent(in)  :: icount
   logical  :: desc_forces_local

   desc_forces_local = desc_forces .and. (config_real(icount)%has_force .or. config_real(icount)%has_stress)
   ! Computing the first kind of desc ...
   if (allocated(config_desc(icount)%energy)) deallocate (config_desc(icount)%energy)
   allocate (config_desc(icount)%energy(dim_xdesc1, imm))
   if (desc_forces_local) then
      if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
         if (allocated(config_desc(icount)%force)) deallocate (config_desc(icount)%force)
         allocate (config_desc(icount)%force(dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 3))
      end if
   end if

   dim_xdesc = config_desc(icount)%dim_desc1
   config_desc(icount)%dim_desc = dim_xdesc
   call compute_g2(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)
#if(PARA)
   !call MPI_BARRIER(subworld, codeml)
   call my_barrier_subworld(codeml)
#endif

   if (allocated(config_desc(icount)%energy1)) deallocate (config_desc(icount)%energy1)
   allocate (config_desc(icount)%energy1(dim_xdesc1, imm))
   config_desc(icount)%energy1(:, :) = 0.d0
   config_desc(icount)%energy1(1:dim_xdesc1, 1:imm) = config_desc(icount)%energy(1:dim_xdesc1, 1:imm)
   if (desc_forces_local) then
      if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
         if (allocated(config_desc(icount)%force1)) deallocate (config_desc(icount)%force1)
         allocate (config_desc(icount)%force1(dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 3))
         config_desc(icount)%force1(:, :, :, :) = 0.d0
         config_desc(icount)%force1(1:dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 1:3) = config_desc(icount)%force(1:dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 1:3)
      end if
   end if

   ! Computing the second kind of desc ...
   if (allocated(config_desc(icount)%energy)) deallocate (config_desc(icount)%energy)
   allocate (config_desc(icount)%energy(dim_xdesc2, imm))
   if (desc_forces_local) then
      if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
         if (allocated(config_desc(icount)%force)) deallocate (config_desc(icount)%force)
         allocate (config_desc(icount)%force(dim_xdesc2, i_start_at:i_final_at, 0:imm_neigh, 3))
      end if
   end if

   dim_xdesc = config_desc(icount)%dim_desc2
   config_desc(icount)%dim_desc = dim_xdesc
   call compute_pow_so4(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)
#if(PARA)
   !call MPI_BARRIER(subworld, codeml)
   call my_barrier_subworld(codeml)
#endif

   if (allocated(config_desc(icount)%energy2)) deallocate (config_desc(icount)%energy2)
   allocate (config_desc(icount)%energy2(dim_xdesc2, imm))
   config_desc(icount)%energy2(:, :) = 0.d0
   config_desc(icount)%energy2(1:dim_xdesc2, 1:imm) = config_desc(icount)%energy(1:dim_xdesc2, 1:imm)
   if (desc_forces_local) then
      if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
         if (allocated(config_desc(icount)%force2)) deallocate (config_desc(icount)%force2)
         allocate (config_desc(icount)%force2(dim_xdesc2, i_start_at:i_final_at, 0:imm_neigh, 3))
         config_desc(icount)%force2(:, :, :, :) = 0.d0
         config_desc(icount)%force2(1:dim_xdesc2, i_start_at:i_final_at, 0:imm_neigh, 1:3) = config_desc(icount)%force(1:dim_xdesc2, i_start_at:i_final_at, 0:imm_neigh, 1:3)
      end if
   end if

   ! Fill the final %energy and %force.
   dim_xdesc = config_desc(icount)%dim_desc1 + config_desc(icount)%dim_desc2
   config_desc(icount)%dim_desc = dim_xdesc
   if (allocated(config_desc(icount)%energy)) deallocate (config_desc(icount)%energy)
   allocate (config_desc(icount)%energy(dim_xdesc, imm))
   config_desc(icount)%energy(:, :) = 0.d0
   config_desc(icount)%energy(1:dim_xdesc1, 1:imm) = config_desc(icount)%energy1(1:dim_xdesc1, 1:imm)
   config_desc(icount)%energy(1 + dim_xdesc1:dim_xdesc, 1:imm) = config_desc(icount)%energy2(1:dim_xdesc2, 1:imm)
   deallocate (config_desc(icount)%energy1)
   deallocate (config_desc(icount)%energy2)

   if (desc_forces_local) then
      if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
         if (allocated(config_desc(icount)%force)) deallocate (config_desc(icount)%force)
         allocate (config_desc(icount)%force(dim_xdesc, i_start_at:i_final_at, 0:imm_neigh, 3))
         config_desc(icount)%force(:, :, :, :) = 0.d0
         config_desc(icount)%force(1:dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 1:3) = config_desc(icount)%force1(1:dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 1:3)
         config_desc(icount)%force(1 + dim_xdesc1:dim_xdesc, i_start_at:i_final_at, 0:imm_neigh, 1:3) = config_desc(icount)%force2(1:dim_xdesc2, i_start_at:i_final_at, 0:imm_neigh, 1:3)
         deallocate (config_desc(icount)%force1)
         deallocate (config_desc(icount)%force2)
      end if
   end if

end subroutine compute_mixed_g2_pow_so4



subroutine compute_mixed_g2_afs(icount)
#if(PARA)
   !use mpi
   use mld_mpi, only: codeml 
#endif
#ifdef MLD_NDM
   use gen_com_m_ml, ONLY: imm
#else
   use ondm_gen_com_m, ONLY: imm
#endif
   use ml_in_ndm_module, only: desc_forces, imm_neigh, i_start_at, i_final_at
   use derived_types, only: config_desc, config_real
   use temporary_data_cov, only: dim_xdesc, dim_xdesc1, dim_xdesc2
   use compute_g2_mod, only: compute_g2
   use compute_afs_mod, only: compute_afs
   !use mld_subworld, only: subworld
   use my_mpi_subroutines, only: my_barrier_subworld

   implicit none


   integer, intent(in)  :: icount
   logical  :: desc_forces_local


   desc_forces_local = desc_forces .and. (config_real(icount)%has_force .or. config_real(icount)%has_stress)
   ! Computing the first kind of desc ...
   if (allocated(config_desc(icount)%energy)) deallocate (config_desc(icount)%energy)
   allocate (config_desc(icount)%energy(dim_xdesc1, imm))
   if (desc_forces_local) then
      if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
         if (allocated(config_desc(icount)%force)) deallocate (config_desc(icount)%force)
         allocate (config_desc(icount)%force(dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 3))
      end if
   end if

   dim_xdesc = config_desc(icount)%dim_desc1
   config_desc(icount)%dim_desc = dim_xdesc
   call compute_g2(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)
#if(PARA)
   !call MPI_BARRIER(subworld, codeml)
   call my_barrier_subworld(codeml)
#endif

   if (allocated(config_desc(icount)%energy1)) deallocate (config_desc(icount)%energy1)
   allocate (config_desc(icount)%energy1(dim_xdesc1, imm))
   config_desc(icount)%energy1(:, :) = 0.d0
   config_desc(icount)%energy1(1:dim_xdesc1, 1:imm) = config_desc(icount)%energy(1:dim_xdesc1, 1:imm)
   if (desc_forces_local) then
      if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
         if (allocated(config_desc(icount)%force1)) deallocate (config_desc(icount)%force1)
         allocate (config_desc(icount)%force1(dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 3))
         config_desc(icount)%force1(:, :, :, :) = 0.d0
         config_desc(icount)%force1(1:dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 1:3) = config_desc(icount)%force(1:dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 1:3)
      end if
   end if

   ! Computing the second kind of desc ...
   if (allocated(config_desc(icount)%energy)) deallocate (config_desc(icount)%energy)
   allocate (config_desc(icount)%energy(dim_xdesc2, imm))
   if (desc_forces_local) then
      if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
         if (allocated(config_desc(icount)%force)) deallocate (config_desc(icount)%force)
         allocate (config_desc(icount)%force(dim_xdesc2, i_start_at:i_final_at, 0:imm_neigh, 3))
      end if
   end if

   dim_xdesc = config_desc(icount)%dim_desc2
   config_desc(icount)%dim_desc = dim_xdesc
   call compute_afs(i_start_at, i_final_at, config_desc(icount)%n_neigh, config_desc(icount)%kind_neigh, icount)

#if(PARA)
   !call MPI_BARRIER(subworld, codeml)
   call my_barrier_subworld(codeml)
#endif

   if (allocated(config_desc(icount)%energy2)) deallocate (config_desc(icount)%energy2)
   allocate (config_desc(icount)%energy2(dim_xdesc2, imm))
   config_desc(icount)%energy2(:, :) = 0.d0
   config_desc(icount)%energy2(1:dim_xdesc2, 1:imm) = config_desc(icount)%energy(1:dim_xdesc2, 1:imm)
   if (desc_forces_local) then
      if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
         if (allocated(config_desc(icount)%force2)) deallocate (config_desc(icount)%force2)
         allocate (config_desc(icount)%force2(dim_xdesc2, i_start_at:i_final_at, 0:imm_neigh, 3))
         config_desc(icount)%force2(:, :, :, :) = 0.d0
         config_desc(icount)%force2(1:dim_xdesc2, i_start_at:i_final_at, 0:imm_neigh, 1:3) = config_desc(icount)%force(1:dim_xdesc2, i_start_at:i_final_at, 0:imm_neigh, 1:3)
      end if
   end if

   ! Fill the final %energy and %force.
   dim_xdesc = config_desc(icount)%dim_desc1 + config_desc(icount)%dim_desc2
   config_desc(icount)%dim_desc = dim_xdesc
   if (allocated(config_desc(icount)%energy)) deallocate (config_desc(icount)%energy)
   allocate (config_desc(icount)%energy(dim_xdesc, imm))
   config_desc(icount)%energy(:, :) = 0.d0
   config_desc(icount)%energy(1:dim_xdesc1, 1:imm) = config_desc(icount)%energy1(1:dim_xdesc1, 1:imm)
   config_desc(icount)%energy(1 + dim_xdesc1:dim_xdesc, 1:imm) = config_desc(icount)%energy2(1:dim_xdesc2, 1:imm)
   deallocate (config_desc(icount)%energy1)
   deallocate (config_desc(icount)%energy2)

   if (desc_forces_local) then
      if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
         if (allocated(config_desc(icount)%force)) deallocate (config_desc(icount)%force)
         allocate (config_desc(icount)%force(dim_xdesc, i_start_at:i_final_at, 0:imm_neigh, 3))
         config_desc(icount)%force(:, :, :, :) = 0.d0
         config_desc(icount)%force(1:dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 1:3) = config_desc(icount)%force1(1:dim_xdesc1, i_start_at:i_final_at, 0:imm_neigh, 1:3)
         config_desc(icount)%force(1 + dim_xdesc1:dim_xdesc, i_start_at:i_final_at, 0:imm_neigh, 1:3) = -config_desc(icount)%force2(1:dim_xdesc2, i_start_at:i_final_at, 0:imm_neigh, 1:3)
         deallocate (config_desc(icount)%force1)
         deallocate (config_desc(icount)%force2)
      end if
   end if

end subroutine compute_mixed_g2_afs

!add something here for mixed descriptors