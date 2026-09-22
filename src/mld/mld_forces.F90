! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   MiLaDy (Machine Learning Dynamics)
! HND X
! HND X
! HND X   Milady was written and designed by Mihai-Cosmin Marinica, Alexandra M. Goryaeva
! HND X   Copyright 2015-2023.
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

module  mld_force_mod

   use mld_logger

contains

   subroutine pack_force_descriptor(iconf)

#if(PARA)
      use mpi
      use mld_mpi
#endif
      use ml_in_ndm_module, only:  rangml, desc_forces, &
         mld_order, mld_linear, mld_linear_extended, &
         mld_quadratic, mld_polyc, mld_kernel,  mld_type_quadratic, mld_type_quadratic_ZX
      use temporary_data_cov, only: dim_xdesc
      use derived_types, only: config_desc, config_real
      !use ondm_gen_com_m, only: imm
      use module_mld_quadratic, only: dim_xdesc_quadratic
      use module_mld_polyc, only: dim_xdesc_polyc
      use module_mld_kernel, only: dim_xdesc_kernel
      use module_kernel_2b, only : activate_k2b, dim_kernel_2b
      use mld_logger
      use mld_subworld, only: subworld

      implicit none

      !external :: para_snap_pack_force_descriptor
      integer, intent(in)  :: iconf
      integer  :: dim_reduce, inat

      _NAMECURRENT_("pack_force_descriptor")



      _MLD_BEGIN_

      inat = config_real(iconf)%nat

      if (activate_k2b) then
         if (allocated(config_desc(iconf)%pack_force_k2b)) deallocate (config_desc(iconf)%pack_force_k2b)
         allocate (config_desc(iconf)%pack_force_k2b(dim_kernel_2b, 3, inat))
         config_desc(iconf)%pack_force_k2b(:, :, :) = 0.d0
      end if


      if ((mld_order == mld_linear) .or. (mld_order == mld_linear_extended)) then
         if (allocated(config_desc(iconf)%pack_force_linear)) deallocate (config_desc(iconf)%pack_force_linear); allocate (config_desc(iconf)%pack_force_linear(dim_xdesc, 3, inat))
         config_desc(iconf)%pack_force_linear(:, :, :) = 0.d0
      end if

      if (mld_order == mld_quadratic) then
         if (allocated(config_desc(iconf)%pack_force_linear)) deallocate (config_desc(iconf)%pack_force_linear); allocate (config_desc(iconf)%pack_force_linear(dim_xdesc, 3, inat))
         if (allocated(config_desc(iconf)%pack_force_quadratic)) deallocate (config_desc(iconf)%pack_force_quadratic); allocate (config_desc(iconf)%pack_force_quadratic(dim_xdesc_quadratic - 1, 3, inat))
         config_desc(iconf)%pack_force_linear(:, :, :) = 0.d0
         config_desc(iconf)%pack_force_quadratic(:, :, :) = 0.d0
      end if


      if (mld_order == mld_polyc) then
         if (allocated(config_desc(iconf)%pack_force_linear)) deallocate (config_desc(iconf)%pack_force_linear); allocate (config_desc(iconf)%pack_force_linear(dim_xdesc, 3, inat))
         if (allocated(config_desc(iconf)%pack_force_polyc)) deallocate (config_desc(iconf)%pack_force_polyc); allocate (config_desc(iconf)%pack_force_polyc(dim_xdesc_polyc - 1, 3, inat))
         config_desc(iconf)%pack_force_linear(:, :, :) = 0.d0
         config_desc(iconf)%pack_force_polyc(:, :, :) = 0.d0
      end if


      if (mld_order == mld_kernel) then
         if (allocated(config_desc(iconf)%pack_force_linear)) deallocate (config_desc(iconf)%pack_force_linear); allocate (config_desc(iconf)%pack_force_linear(dim_xdesc, 3, inat))
         if (allocated(config_desc(iconf)%pack_force_kernel)) deallocate (config_desc(iconf)%pack_force_kernel); allocate (config_desc(iconf)%pack_force_kernel(dim_xdesc_kernel - 1, 3, inat))
         config_desc(iconf)%pack_force_linear(:, :, :) = 0.d0
         config_desc(iconf)%pack_force_kernel(:, :, :) = 0.d0
      end if

      call para_mld_pack_force_descriptor(iconf)



#if(PARA)

      ! WARNING : config_desc(iconf)%force is local to each subworld!! No reduction done on it (unlike the energy)
      ! so reduction on pack_force_linear (which is based on force) is mandatory!
      dim_reduce = dim_xdesc*inat*3
      if (desc_forces) call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(iconf)%pack_force_linear, dim_reduce, MPI_DOUBLE_PRECISION, MPI_SUM, subworld, codeml)
      if (activate_k2b) then
         dim_reduce = dim_kernel_2b*inat*3
         if (desc_forces) call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(iconf)%pack_force_k2b, dim_reduce, MPI_DOUBLE_PRECISION, MPI_SUM, subworld, codeml)
      end if

      if (mld_type_quadratic == mld_type_quadratic_ZX) then
         ! this call has absolutely need from ALLREDUCE on pack_force_linear and pack_force_k2b
         call quadratic_ZX_pack_force_descriptor(iconf)
      end if

      if ((mld_order == mld_linear) .or. (mld_order == mld_linear_extended)) then
      end if

      if (mld_order == mld_quadratic) then
         dim_reduce = (dim_xdesc_quadratic - 1)*inat*3
         if (desc_forces) call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(iconf)%pack_force_quadratic, dim_reduce, MPI_DOUBLE_PRECISION, MPI_SUM, subworld, codeml)
      end if


      if (mld_order == mld_polyc) then
         dim_reduce = (dim_xdesc_polyc - 1)*inat*3
         if (desc_forces) call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(iconf)%pack_force_polyc, dim_reduce, MPI_DOUBLE_PRECISION, MPI_SUM, subworld, codeml)
      end if


      if (mld_order == mld_kernel) then
         dim_reduce = (dim_xdesc_kernel - 1)*inat*3
         if (desc_forces) call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(iconf)%pack_force_kernel, dim_reduce, MPI_DOUBLE_PRECISION, MPI_SUM, subworld, codeml)
      end if



#endif

      _MLD_END_

   end subroutine pack_force_descriptor




   subroutine para_mld_pack_force_descriptor(iconf)

      !--------------------------------------------------------------------!
      ! Pack the descritor for forces once the derivatives are computed.   !
      ! Input: The number of configuration: iconf                          !
      ! Input: The derivatives through the config_desc(iconf)%force vector !
      ! Output: the vector config_desc(iconf)%pack_force(dim_of_desc,3,imm)!
      !--------------------------------------------------------------------!
      use derived_types, only: config_real
      use ml_in_ndm_module, only:  desc_forces, &
         mld_order, mld_linear, mld_linear_extended, &
         mld_quadratic, mld_polyc, &
         mld_kernel, mld_type_quadratic, mld_type_quadratic_ZX
      use module_kernel_2b, only: activate_k2b
      use mld_logger

      implicit none

      integer, intent(in)  :: iconf

      _NAMECURRENT_("para_mld_pack_force_descriptor")



      _MLD_BEGIN_
      if (.not. desc_forces) return
      if (.not. (config_real(iconf)%has_force)) return
      ! if (allocated(config_desc(iconf)%pack_force)) deallocate(config_desc(iconf)%pack_force ); allocate(config_desc(iconf)%pack_force (dim_xdesc,3,imm))

      if ((mld_order == mld_linear) .or. (mld_order == mld_linear_extended)) then
         call linear_pack_force_descriptor(iconf)
      end if

      if (mld_order == mld_quadratic) then
         if (mld_type_quadratic == mld_type_quadratic_ZX) then
            ! this is the first call ... in order to agredate global
            call linear_pack_force_descriptor(iconf)
         else
            call quadratic_pack_force_descriptor(iconf)
         end if
      end if

      if (mld_order == mld_polyc) then
         call polyc_pack_force_descriptor(iconf)
      end if


      if (mld_order == mld_kernel) then
         ! call linear_pack_force_descriptor(iconf)
         call kernel_pack_force_descriptor(iconf)
      end if

      if (activate_k2b) then
         ! this is good for k2b but also for Z X global local coupling. 
         call pack_force_k2b(iconf)
      end if

      _MLD_END_
   end subroutine para_mld_pack_force_descriptor




   subroutine linear_pack_force_descriptor(iconf)

      !--------------------------------------------------------------------!
      ! Pack the descritor for forces once the derivatives are computed.   !
      ! Input: The number of configuration: iconf                          !
      ! Input: The derivatives through the config_desc(iconf)%force vector !
      ! Output: the vector config_desc(iconf)%pack_force(dim_of_desc,3,imm)!
      !--------------------------------------------------------------------!
      use temporary_data_cov, only: dim_xdesc
      use derived_types, only: config_desc, config_real
      use ml_in_ndm_module, only: i_start_at, i_final_at
      use mld_logger

      implicit none

      integer, intent(in)  :: iconf
      integer  :: ik, ia, inn1
      logical  :: small

      _NAMECURRENT_("linear_pack_force_descriptor")



      _MLD_BEGIN_
      if ((i_start_at == 0) .and. (i_final_at == 0)) return
      small = config_real(iconf)%small
      ! if (allocated(config_desc(iconf)%pack_force)) deallocate(config_desc(iconf)%pack_force ); allocate(config_desc(iconf)%pack_force (dim_xdesc,3,imm))
      do ik = i_start_at, i_final_at                   ! derivative index \partial / \partial ik,alpha
         config_desc(iconf)%pack_force_linear(1:dim_xdesc, 1:3, ik) = config_desc(iconf)%pack_force_linear(1:dim_xdesc, 1:3, ik) - config_desc(iconf)%force(1:dim_xdesc, ik, 0, 1:3)
         ! icount1=0
         do inn1 = 1, config_desc(iconf)%n_neigh(ik)      ! sum over all atoms
            ia = config_desc(iconf)%kind_neigh(ik, inn1)
            !if (small) then
            !  if (.not.(config_desc(iconf)%incell(ik,inn1))) cycle
            !end if
            !if (small) config_desc(iconf)%pack_force(1:dim_xdesc,1:3,ia) = config_desc(iconf)%pack_force(1:dim_xdesc,1:3,ia) + config_desc(iconf)%force(1:dim_xdesc,ik,inn1,1:3 )
            config_desc(iconf)%pack_force_linear(1:dim_xdesc, 1:3, ia) = config_desc(iconf)%pack_force_linear(1:dim_xdesc, 1:3, ia) - config_desc(iconf)%force(1:dim_xdesc, ik, inn1, 1:3)
         end do
      end do

      _MLD_END_
   end subroutine linear_pack_force_descriptor




   subroutine kernel_pack_force_descriptor(iconf)

      !--------------------------------------------------------------------!
      ! Pack the descritor for forces once the derivatives are computed.   !
      ! Input: The number of configuration: iconf                          !
      ! Input: The derivatives through the config_desc(iconf)%force vector !
      ! Output: the vector config_desc(iconf)%pack_force(dim_of_desc,3,imm)!
      !--------------------------------------------------------------------!
      use temporary_data_cov, only: dim_xdesc
      use derived_types, only: config_desc, config_real
      use ml_in_ndm_module, only: i_start_at, i_final_at
      use module_kernel, only: dim_kernel
      use mld_logger

      implicit none

      integer, intent(in)  :: iconf
      integer  :: ik, ia, inn1
      logical  :: small

      _NAMECURRENT_("kernel_pack_force_descriptor")



      _MLD_BEGIN_
      ! write (*,*) 'FILE............   ', config_real(iconf)%class, config_real(iconf)%filename
      if ((i_start_at == 0) .and. (i_final_at == 0)) return
      small = config_real(iconf)%small
      ! if (allocated(config_desc(iconf)%pack_force))     deallocate(config_desc(iconf)%pack_force )     ; allocate(config_desc(iconf)%pack_force (dim_xdesc,3,imm))
      do ik = i_start_at, i_final_at                   ! derivative index \partial / \partial ik,alpha
         config_desc(iconf)%pack_force_linear(1:dim_xdesc, 1:3, ik) = config_desc(iconf)%pack_force_linear(1:dim_xdesc, 1:3, ik) - config_desc(iconf)%force(1:dim_xdesc, ik, 0, 1:3)
         config_desc(iconf)%pack_force_kernel(dim_xdesc + 1:dim_xdesc + dim_kernel, 1:3, ik) = config_desc(iconf)%pack_force_kernel(dim_xdesc + 1:dim_xdesc + dim_kernel, 1:3, ik) - config_desc(iconf)%force_kernel(1:dim_kernel, ik, 0, 1:3)
         ! icount1=0
         ! write (*,*) 'ik...', ik
         do inn1 = 1, config_desc(iconf)%n_neigh(ik)      ! sum over all atoms ...
            ia = config_desc(iconf)%kind_neigh(ik, inn1)
            !write (*,*) 'ik  inn1  ia', ik, inn1, ia, allocated(config_desc(iconf)%force_kernel)
            !if (small) then
            !  if (.not.(config_desc(iconf)%incell(ik,inn1))) cycle
            !end if
            !if (small) config_desc(iconf)%pack_force(1:dim_xdesc,1:3,ia) = config_desc(iconf)%pack_force(1:dim_xdesc,1:3,ia) + config_desc(iconf)%force(1:dim_xdesc,ik,inn1,1:3 )
            config_desc(iconf)%pack_force_linear(1:dim_xdesc, 1:3, ia) = config_desc(iconf)%pack_force_linear(1:dim_xdesc, 1:3, ia) - config_desc(iconf)%force(1:dim_xdesc, ik, inn1, 1:3)
            config_desc(iconf)%pack_force_kernel(dim_xdesc + 1:dim_xdesc + dim_kernel, 1:3, ia) = config_desc(iconf)%pack_force_kernel(dim_xdesc + 1:dim_xdesc + dim_kernel, 1:3, ia) - config_desc(iconf)%force_kernel(1:dim_kernel, ik, inn1, 1:3)
         end do
      end do


      !compatibility reasons
      config_desc(iconf)%pack_force_kernel(1:dim_xdesc, 1:3, :) = config_desc(iconf)%pack_force_linear(1:dim_xdesc, 1:3, :)
      ! if (allocated(config_desc(iconf)%force_kernel)) deallocate(config_desc(iconf)%force_kernel)

      _MLD_END_
   end subroutine kernel_pack_force_descriptor

   subroutine pack_force_k2b(iconf)

      !--------------------------------------------------------------------!
      ! Pack the descritor for forces once the derivatives are computed.   !
      ! Input: The number of configuration: iconf                          !
      ! Input: The derivatives through the config_desc(iconf)%force vector !
      ! Output: the vector config_desc(iconf)%pack_force(dim_of_desc,3,imm)!
      !--------------------------------------------------------------------!
      use derived_types, only: config_desc, config_real
      use ml_in_ndm_module, only:  i_start_at, i_final_at
      use module_kernel_2b, only: dim_kernel_2b
      use mld_logger

      implicit none

      integer, intent(in)  :: iconf
      integer  :: ik, ia, inn1
      logical  :: small

      _NAMECURRENT_("pack_force_k2b")

      _MLD_BEGIN_
      if ((i_start_at == 0) .and. (i_final_at == 0)) return
      small = config_real(iconf)%small
      do ik = i_start_at, i_final_at                   ! derivative index \partial / \partial ik,alpha
         config_desc(iconf)%pack_force_k2b(1:dim_kernel_2b, 1:3, ik) = config_desc(iconf)%pack_force_k2b(1:dim_kernel_2b, 1:3, ik) - config_desc(iconf)%force_k2b(1:dim_kernel_2b, ik, 0, 1:3)
         do inn1 = 1, config_desc(iconf)%n_neigh(ik)      ! sum over all atoms ...
            ia = config_desc(iconf)%kind_neigh(ik, inn1)
            config_desc(iconf)%pack_force_k2b(1:dim_kernel_2b, 1:3, ia) = config_desc(iconf)%pack_force_k2b(1:dim_kernel_2b, 1:3, ia) - config_desc(iconf)%force_k2b(1:dim_kernel_2b, ik, inn1, 1:3)
         end do
      end do

      _MLD_END_
   end subroutine pack_force_k2b

   subroutine quadratic_ZX_pack_force_descriptor(iconf)
      !--------------------------------------------------------------------!
      ! Pack the descritor for forces once the derivatives are computed.   !
      ! Input: The number of configuration: iconf                          !
      ! Input: The derivatives through the config_desc(iconf)%force vector !
      ! Output: the vector config_desc(iconf)%pack_force(dim_of_desc,3,imm)!
      !--------------------------------------------------------------------!
      use module_kind_variables, only: kind_double
      use temporary_data_cov, only: dim_xdesc
      use derived_types, only: config_desc
      use ml_in_ndm_module, only: i_start_at, i_final_at
      use module_kernel_2b, only: dim_kernel_2b                            
      use mld_logger

      implicit none

      integer, intent(in)  :: iconf
      integer  :: ix, ik, id1, id2, indx
      real(kind_double), dimension(dim_xdesc, 1) :: Dtmp
      real(kind_double), dimension(dim_kernel_2b, 1) :: Ztmp
      real(kind_double) :: temp, Znorm, Dnorm, tmp1(3), tmp2(3), tmp_ik_energy(3), tmp_ik_k2b(3) 

      !_NAMECURRENT_("quadratic_ZX_pack_force_descriptor")

         Dtmp(1:dim_xdesc, 1) = config_desc(iconf)%pack_energy_linear(1:dim_xdesc)
         Ztmp(1:dim_kernel_2b, 1) = config_desc(iconf)%pack_energy_k2b(1:dim_kernel_2b)
         Dnorm = dsqrt ( DOT_PRODUCT (Dtmp(1:dim_xdesc, 1), Dtmp(1:dim_xdesc, 1)))
         Znorm = dsqrt ( DOT_PRODUCT (Ztmp(1:dim_kernel_2b,1),Ztmp(1:dim_kernel_2b,1)))
         temp = (1.d0/Dnorm + 1.d0/Znorm)

         do  ik = i_start_at, i_final_at
           do ix =1,3 
             tmp_ik_energy(ix) = - dot_product(config_desc(iconf)%pack_energy_linear(:), config_desc(iconf)%pack_force_linear(:, ix, ik))/Dnorm**3
             tmp_ik_k2b(ix) = - dot_product(config_desc(iconf)%pack_energy_k2b(:), config_desc(iconf)%pack_force_k2b(:, ix, ik))/Znorm**3
           end do 

           do id1 = 1, dim_xdesc
              do id2 = 1, dim_kernel_2b
  
                 indx = dim_kernel_2b*(id1 - 1) + id2 - 1 + dim_xdesc + 1
                 tmp1(1:3)  = temp* (config_desc(iconf)%pack_force_linear(id1, 1:3, ik)*config_desc(iconf)%pack_energy_k2b(id2) &
                                     + config_desc(iconf)%pack_force_k2b(id2, 1:3, ik)*config_desc(iconf)%pack_energy_linear(id1)) 
                 tmp2(1:3) = Dtmp(id1, 1) * Ztmp(id2, 1)  * ( tmp_ik_energy(1:3) + tmp_ik_k2b(1:3)) 
                  config_desc(iconf)%pack_force_quadratic(indx, 1:3, ik) = config_desc(iconf)%pack_force_quadratic(indx, 1:3, ik) & 
                                -  tmp1(1:3)  - tmp2(1:3) 
                   
               end do      
            end do
         end do   

   end subroutine quadratic_ZX_pack_force_descriptor


   subroutine quadratic_pack_force_descriptor(iconf)
      !--------------------------------------------------------------------!
      ! Pack the descritor for forces once the derivatives are computed.   !
      ! Input: The number of configuration: iconf                          !
      ! Input: The derivatives through the config_desc(iconf)%force vector !
      ! Output: the vector config_desc(iconf)%pack_force(dim_of_desc,3,imm)!
      !--------------------------------------------------------------------!

      use temporary_data_cov, only: dim_xdesc
      use derived_types, only: config_desc
      use ml_in_ndm_module, only: i_start_at, i_final_at, & 
                                  mld_type_quadratic, mld_type_quadratic_zaxa
      use module_kernel_2b, only: dim_kernel_2b                            
      use mld_logger

      implicit none

      integer, intent(in)  :: iconf
      integer  :: ik, ia, id1, id2, indx, inn1

      _NAMECURRENT_("quadratic_pack_force_descriptor")



      _MLD_BEGIN_
      if ((i_start_at == 0) .and. (i_final_at == 0)) return
      ! if (allocated(config_desc(iconf)%pack_force)) deallocate(config_desc(iconf)%pack_force ); allocate(config_desc(iconf)%pack_force (dim_xdesc,3,imm))

      do ik = i_start_at, i_final_at                   ! derivative index \partial / \partial ik,alpha

         do id1 = 1, dim_xdesc
            indx = id1
            config_desc(iconf)%pack_force_quadratic(indx, 1:3, ik) = config_desc(iconf)%pack_force_quadratic(indx, 1:3, ik) - config_desc(iconf)%force(indx, ik, 0, 1:3)
            ! icount1=0
            do inn1 = 1, config_desc(iconf)%n_neigh(ik)      ! sum over all atoms ...
               ia = config_desc(iconf)%kind_neigh(ik, inn1)
               config_desc(iconf)%pack_force_quadratic(indx, 1:3, ia) = config_desc(iconf)%pack_force_quadratic(indx, 1:3, ia) - config_desc(iconf)%force(indx, ik, inn1, 1:3)
            end do
         end do

      if (mld_type_quadratic == mld_type_quadratic_zaxa) then 
           do id1 = 1, dim_xdesc
              do id2 = 1, dim_kernel_2b
  
                 indx = dim_kernel_2b*(id1 - 1) + id2 - 1 + dim_xdesc + 1
                 config_desc(iconf)%pack_force_quadratic(indx, 1:3, ik) = config_desc(iconf)%pack_force_quadratic(indx, 1:3, ik) &
                    - config_desc(iconf)%force(id1, ik, 0, 1:3)*config_desc(iconf)%energy_k2b(id2, ik) &
                    - config_desc(iconf)%force_k2b(id2, ik, 0, 1:3)*config_desc(iconf)%energy(id1, ik)
                 ! icount1=0
                 do inn1 = 1, config_desc(iconf)%n_neigh(ik)      ! sum over all atoms ...
                    ia = config_desc(iconf)%kind_neigh(ik, inn1)
                    config_desc(iconf)%pack_force_quadratic(indx, 1:3, ia) = config_desc(iconf)%pack_force_quadratic(indx, 1:3, ia) &
                       - config_desc(iconf)%force(id1, ik, inn1, 1:3)*config_desc(iconf)%energy_k2b(id2, ik) &
                       - config_desc(iconf)%force_k2b(id2, ik, inn1, 1:3)*config_desc(iconf)%energy(id1, ik)
                    !-  config_desc(iconf)%force(id1,ik,inn1,1:3 )*config_desc(iconf)%energy(id2,ia) &
                    !-  config_desc(iconf)%force(id2,ik,inn1,1:3 )*config_desc(iconf)%energy(id1,ia)
                 end do
  
              end do
           end do
      else 
           do id1 = 1, dim_xdesc
              do id2 = 1, dim_xdesc
  
                 indx = dim_xdesc*(id1 - 1) + id2 - 1 + dim_xdesc + 1
                 config_desc(iconf)%pack_force_quadratic(indx, 1:3, ik) = config_desc(iconf)%pack_force_quadratic(indx, 1:3, ik) &
                    - config_desc(iconf)%force(id1, ik, 0, 1:3)*config_desc(iconf)%energy(id2, ik) &
                    - config_desc(iconf)%force(id2, ik, 0, 1:3)*config_desc(iconf)%energy(id1, ik)
                 ! icount1=0
                 do inn1 = 1, config_desc(iconf)%n_neigh(ik)      ! sum over all atoms ...
                    ia = config_desc(iconf)%kind_neigh(ik, inn1)
                    config_desc(iconf)%pack_force_quadratic(indx, 1:3, ia) = config_desc(iconf)%pack_force_quadratic(indx, 1:3, ia) &
                       - config_desc(iconf)%force(id1, ik, inn1, 1:3)*config_desc(iconf)%energy(id2, ik) &
                       - config_desc(iconf)%force(id2, ik, inn1, 1:3)*config_desc(iconf)%energy(id1, ik)
                    !-  config_desc(iconf)%force(id1,ik,inn1,1:3 )*config_desc(iconf)%energy(id2,ia) &
                    !-  config_desc(iconf)%force(id2,ik,inn1,1:3 )*config_desc(iconf)%energy(id1,ia)
                 end do
  
              end do
           end do
         end if
         end do

      !compatibility reasons ...HERE IS NOT CORRECT - I THINK. NOW IT IS.
      ! do ik=i_start_at, i_final_at
      config_desc(iconf)%pack_force_linear(1:dim_xdesc, 1:3, :) = config_desc(iconf)%pack_force_quadratic(1:dim_xdesc, 1:3, :)
      ! end do

      _MLD_END_

   end subroutine quadratic_pack_force_descriptor




   subroutine polyc_pack_force_descriptor(iconf)
      !--------------------------------------------------------------------!
      ! Pack the descritor for forces once the derivatives are computed.   !
      ! Input: The number of configuration: iconf                          !
      ! Input: The derivatives through the config_desc(iconf)%force vector !
      ! Output: the vector config_desc(iconf)%pack_force(dim_of_desc,3,imm)!
      !--------------------------------------------------------------------!

      use temporary_data_cov, only: dim_xdesc
      use derived_types, only: config_desc, config_real
      use ml_in_ndm_module, only: i_start_at, i_final_at, polyc_n_poly, polyc_n_hermite
      use math, only: hermite
      use mld_logger

      implicit none

      integer, intent(in)  :: iconf
      real(kind=kind(0.d0))      :: y_her, d_y_her
      integer  :: ik, ia, nord, id1, id2, id3, indx, indx_last, inn1
      real(kind=kind(0.d0)), dimension(:, :, :), allocatable   :: etmp, d_etmp

      _NAMECURRENT_("polyc_pack_force_descriptor")



      _MLD_BEGIN_
      if ((i_start_at == 0) .and. (i_final_at == 0)) return
      ! if (allocated(config_desc(iconf)%pack_force)) deallocate(config_desc(iconf)%pack_force ); allocate(config_desc(iconf)%pack_force (dim_xdesc,3,imm))
      if (allocated(etmp)) deallocate (etmp); allocate (etmp(polyc_n_hermite, dim_xdesc, config_real(iconf)%nat))
      if (allocated(d_etmp)) deallocate (d_etmp); allocate (d_etmp(polyc_n_hermite, dim_xdesc, config_real(iconf)%nat))

      ! for energy we need from all the procs, in Lammps is different because we have the neighbouring topology.
      do ia = 1, config_real(iconf)%nat
         do nord = 1, polyc_n_hermite
            do id1 = 1, dim_xdesc
               call hermite(nord, config_desc(iconf)%energy(id1, ia), y_her, d_y_her)
               etmp(nord, id1, ia) = y_her
               d_etmp(nord, id1, ia) = d_y_her
            end do
         end do
      end do



      do ik = i_start_at, i_final_at                   ! derivative index \partial / \partial ik,alpha

         !do nord = 1 , polyc_n_hermite
         !    do id1 = 1, dim_xdesc
         !        call hermite(nord,  config_desc(iconf)%energy(id1,ik) , y_her, d_y_her)
         !        Dtmp(nord, id1)=y_her
         !        d_Dtmp(nord, id1)=d_y_her
         !    end do
         !end do

         indx_last = 0
         if (polyc_n_poly >= 1) then
            indx = indx_last
            do nord = 1, polyc_n_hermite
               do id1 = 1, dim_xdesc
                  indx = indx + 1
                  config_desc(iconf)%pack_force_polyc(indx, 1:3, ik) = config_desc(iconf)%pack_force_polyc(indx, 1:3, ik) &
                     - config_desc(iconf)%force(id1, ik, 0, 1:3)*d_etmp(nord, id1, ik)
                  !icount1=0
                  do inn1 = 1, config_desc(iconf)%n_neigh(ik)      ! sum over all atoms ...
                     ia = config_desc(iconf)%kind_neigh(ik, inn1)
                     config_desc(iconf)%pack_force_polyc(indx, 1:3, ia) = config_desc(iconf)%pack_force_polyc(indx, 1:3, ia) &
                        - config_desc(iconf)%force(id1, ik, inn1, 1:3)*d_etmp(nord, id1, ia)
                  end do
               end do
            end do
            indx_last = indx
         end if


         if (polyc_n_poly >= 2) then
            indx = indx_last
            do nord = 1, polyc_n_hermite
               do id1 = 1, dim_xdesc
                  do id2 = 1, dim_xdesc
                     indx = indx + 1
                     config_desc(iconf)%pack_force_polyc(indx, 1:3, ik) = config_desc(iconf)%pack_force_polyc(indx, 1:3, ik) &
                        - config_desc(iconf)%force(id1, ik, 0, 1:3)*d_etmp(nord, id1, ik)*etmp(nord, id2, ik) &
                        - config_desc(iconf)%force(id2, ik, 0, 1:3)*d_etmp(nord, id2, ik)*etmp(nord, id1, ik)
                     !icount1=0
                     do inn1 = 1, config_desc(iconf)%n_neigh(ik)      ! sum over all atoms ...
                        ia = config_desc(iconf)%kind_neigh(ik, inn1)
                        config_desc(iconf)%pack_force_polyc(indx, 1:3, ia) = config_desc(iconf)%pack_force_polyc(indx, 1:3, ia) &
                        ! -  config_desc(iconf)%force(id1,ik,inn1,1:3 )*d_etmp(nord,id1,ik)*etmp(nord,id2,ia) &
                        ! -  config_desc(iconf)%force(id2,ik,inn1,1:3 )*d_etmp(nord,id2,ik)*etmp(nord,id1,ia)
                        !
                           - config_desc(iconf)%force(id1, ik, inn1, 1:3)*d_etmp(nord, id1, ik)*etmp(nord, id2, ik) &
                           - config_desc(iconf)%force(id2, ik, inn1, 1:3)*d_etmp(nord, id2, ik)*etmp(nord, id1, ik)


                     end do
                  end do
               end do
            end do
            indx_last = indx
         end if


         if (polyc_n_poly >= 3) then
            indx = indx_last
            do nord = 1, polyc_n_hermite
               do id1 = 1, dim_xdesc
                  do id2 = 1, dim_xdesc
                     do id3 = 1, dim_xdesc
                        indx = indx + 1
                        config_desc(iconf)%pack_force_polyc(indx, 1:3, ik) = config_desc(iconf)%pack_force_polyc(indx, 1:3, ik) &
                           - config_desc(iconf)%force(id1, ik, 0, 1:3)*d_etmp(nord, id1, ik)*etmp(nord, id2, ik)*etmp(nord, id3, ik) &
                           - config_desc(iconf)%force(id2, ik, 0, 1:3)*d_etmp(nord, id2, ik)*etmp(nord, id1, ik)*etmp(nord, id3, ik) &
                           - config_desc(iconf)%force(id3, ik, 0, 1:3)*d_etmp(nord, id3, ik)*etmp(nord, id1, ik)*etmp(nord, id2, ik)
                        !icount1=0
                        do inn1 = 1, config_desc(iconf)%n_neigh(ik)      ! sum over all atoms ...
                           ia = config_desc(iconf)%kind_neigh(ik, inn1)
                           config_desc(iconf)%pack_force_polyc(indx, 1:3, ia) = config_desc(iconf)%pack_force_polyc(indx, 1:3, ia) &
                              - config_desc(iconf)%force(id1, ik, inn1, 1:3)*d_etmp(nord, id1, ik)*etmp(nord, id2, ik)*etmp(nord, id3, ik) &
                              - config_desc(iconf)%force(id2, ik, inn1, 1:3)*d_etmp(nord, id2, ik)*etmp(nord, id1, ik)*etmp(nord, id3, ik) &
                              - config_desc(iconf)%force(id3, ik, inn1, 1:3)*d_etmp(nord, id3, ik)*etmp(nord, id1, ik)*etmp(nord, id2, ik)
                           !-  config_desc(iconf)%force(id1,ik,inn1,1:3 )*d_etmp(nord,id1,ik)*etmp(nord,id2,ia)*etmp(nord,id3,ia) &
                           !-  config_desc(iconf)%force(id2,ik,inn1,1:3 )*d_etmp(nord,id2,ik)*etmp(nord,id1,ia)*etmp(nord,id3,ia) &
                           !-  config_desc(iconf)%force(id3,ik,inn1,1:3 )*d_etmp(nord,id3,ik)*etmp(nord,id1,ia)*etmp(nord,id2,ia)
                        end do
                     end do
                  end do
               end do
            end do
            indx_last = indx
         end if

      end do                  ! ik loop ...

      !compatibility reasons ......HERE IS NOT CORRECT - I THINK. NOW IT IS.
      ! do ik=i_start_at, i_final_at
      config_desc(iconf)%pack_force_linear(1:dim_xdesc, 1:3, :) = config_desc(iconf)%pack_force_polyc(1:dim_xdesc, 1:3, :)
      ! end do

      _MLD_END_
   end subroutine polyc_pack_force_descriptor




   subroutine train_fill_Amat_with_force(iconf, pack_opt)

      use module_kind_variables, only: kind_double
      use ml_in_ndm_module, only: desc_forces, &
         mld_order, mld_linear, mld_linear_extended, &
         mld_quadratic, mld_polyc, mld_kernel
      use temporary_data_cov, only: yfunc_train, dim_xdesc, dim_xdesc_patch
      use derived_types, only: config_desc, config_real
      use snap, only: Amat, ymat, i_fit_snap, fit_snap, weights_snap, dim_design_line, dim_xdesc_linear, y_zbl_train
      use mld_energy_mod
      use module_mld_quadratic, only: dim_xdesc_quadratic
      use module_mld_polyc, only: dim_xdesc_polyc
      use module_mld_kernel, only: dim_xdesc_kernel
      use module_init_ScaMatrix, only: init_ScaMatrix_with_LocalMatrix
      use module_kernel_2b, only: activate_k2b, dim_kernel_2b
      use module_zbl, only: zbl_potential, zbl_type, zbl_mode_alone
      use mld_logger

      implicit none

      integer, intent(in)  :: iconf
      integer  :: ik, ix, i_fit_mld_ini, i_fit_mld_end
      logical, intent(in), optional    :: pack_opt
      logical  :: lpack_tmp
      real(kind_double), dimension(:,:), allocatable :: local_mat

      _NAMECURRENT_("train_fill_Amat_with_force")



      _MLD_BEGIN_
      if (.not. desc_forces) return

      lpack_tmp = .true.
      if (present(pack_opt)) lpack_tmp = pack_opt
      if (.not. (config_real(iconf)%has_force)) return

      i_fit_mld_ini = i_fit_snap
      i_fit_mld_end = i_fit_snap + 3*config_real(iconf)%nat

      fit_snap(i_fit_mld_ini + 1:i_fit_mld_end)%force = .true.
      fit_snap(i_fit_mld_ini + 1:i_fit_mld_end)%iconf = iconf
      fit_snap(i_fit_mld_ini + 1:i_fit_mld_end)%weight = config_real(iconf)%w_f_model
      weights_snap(i_fit_mld_ini + 1:i_fit_mld_end) = config_real(iconf)%w_f_model


      ! if has the energy the desc are already computed. This is for case where is not energy ... but we need desc of energy.
      if (.not. (config_real(iconf)%has_energy)) call pack_energy_descriptor(iconf)
      if (lpack_tmp) call pack_force_descriptor(iconf)

      do ik = 1, config_real(iconf)%nat

         yfunc_train(i_fit_snap + 1:i_fit_snap + 3) = config_real(iconf)%force(1:3, ik)

         do ix =1,3
            fit_snap(i_fit_snap + ix)%ix = ix
            fit_snap(i_fit_snap + ix)%iatom = ik
         end do

         i_fit_snap = i_fit_snap + 3

      end do

      if (allocated(local_mat)) deallocate(local_mat) ; allocate(local_mat(dim_design_line, 3))

      i_fit_snap = i_fit_mld_ini

      local_mat(1, 1:3) = 0.d0

      do ik = 1, config_real(iconf)%nat
         if ((mld_order == mld_linear) .or. (mld_order == mld_linear_extended)) then
            local_mat(2 : dim_xdesc_linear, 1:3) = config_desc(iconf)%pack_force_linear(1:dim_xdesc_linear-1, 1:3, ik)
         end if


         if (mld_order == mld_quadratic) then
            local_mat(2 : dim_xdesc_linear, 1:3) = config_desc(iconf)%pack_force_linear(1:dim_xdesc_linear-1, 1:3, ik)
            local_mat(dim_xdesc_linear + dim_xdesc_patch + 1  : dim_design_line, 1:3) = config_desc(iconf)%pack_force_quadratic(dim_xdesc_linear:dim_xdesc_quadratic - 1, 1:3, ik)
         end if

         if (mld_order == mld_polyc) then
            local_mat(2 : dim_xdesc_linear, 1:3) = config_desc(iconf)%pack_force_linear(1:dim_xdesc_linear-1, 1:3, ik)
            local_mat(dim_xdesc_linear + dim_xdesc_patch + 1 : dim_design_line, 1:3) = config_desc(iconf)%pack_force_polyc(dim_xdesc_linear : dim_xdesc_polyc - 1, 1:3, ik)
         end if

         if (mld_order == mld_kernel) then
            local_mat(2 : dim_xdesc_linear, 1:3) = config_desc(iconf)%pack_force_linear(1:dim_xdesc_linear-1, 1:3, ik)
            local_mat(dim_xdesc_linear + dim_xdesc_patch + 1 : dim_design_line, 1:3) = config_desc(iconf)%pack_force_kernel(dim_xdesc_linear : dim_xdesc_kernel - 1, 1:3, ik)
         end if

         if (activate_k2b) then
            local_mat(dim_xdesc_linear + 1: dim_xdesc_linear + dim_kernel_2b, 1:3) = config_desc(iconf)%pack_force_k2b(1: dim_kernel_2b, 1:3, ik)
         end if

         Amat(1, i_fit_snap + 1:i_fit_snap + 3) = local_mat(1, 1:3)
         Amat(2:dim_design_line, i_fit_snap + 1:i_fit_snap + 3) = local_mat(2 : dim_design_line, 1:3)

         if (zbl_potential .and. zbl_type == zbl_mode_alone) then
            ymat(i_fit_snap+1:i_fit_snap+3, 1) = yfunc_train(i_fit_snap+1:i_fit_snap+3) - config_real(iconf)%fzbl(1:3,ik)
            y_zbl_train(i_fit_snap+1:i_fit_snap+3) = config_real(iconf)%fzbl(1:3,ik)
         else
            ymat(i_fit_snap+1:i_fit_snap+3, 1) = yfunc_train(i_fit_snap+1:i_fit_snap+3)
         end if

         i_fit_snap = i_fit_snap + 3

      end do

      i_fit_snap = i_fit_mld_end

      _MLD_END_

   end subroutine train_fill_Amat_with_force


!$! ! OLD STUFFS .UNUSED
!$! 
!$!    subroutine serial_snap_pack_force_descriptor(iconf)
!$! 
!$!       !--------------------------------------------------------------------!
!$!       ! Pack the descritor for forces once the derivatives are computed.   !
!$!       ! Input: The number of configuration: iconf                          !
!$!       ! Input: The derivatives through the config_desc(iconf)%force vector !
!$!       ! Output: the vector config_desc(iconf)%pack_force(dim_of_desc,3,imm)!
!$!       !--------------------------------------------------------------------!
!$!       use ml_in_ndm_module, only: desc_forces
!$!       use ondm_gen_com_m, only: imm
!$!       use temporary_data_cov, only: dim_xdesc
!$!       use derived_types, only: config_desc, config_real
!$!       use mld_logger
!$! 
!$!       implicit none
!$! 
!$!       integer, intent(in)  :: iconf
!$!       integer  :: ik, ia, ib, inn1, inn2
!$!       real(kind(0.d0))     :: temp(dim_xdesc, 3)
!$! 
!$!       if (.not. desc_forces) return
!$!       if (.not. (config_real(iconf)%has_force)) return
!$!       if (allocated(config_desc(iconf)%pack_force_linear)) deallocate (config_desc(iconf)%pack_force_linear); allocate (config_desc(iconf)%pack_force_linear(dim_xdesc, 3, imm))
!$! 
!$!       do ik = 1, config_real(iconf)%nat                ! derivative index \partial / \partial ik,alpha
!$!          temp(:, :) = 0.d0
!$!          temp(1:dim_xdesc, 1:3) = temp(1:dim_xdesc, 1:3) + config_desc(iconf)%force(1:dim_xdesc, ik, 0, 1:3)
!$!          do inn1 = 1, config_desc(iconf)%n_neigh(ik)      ! sum over all atoms ...
!$!             ia = config_desc(iconf)%kind_neigh(ik, inn1)
!$!             !if (ia==ik) then ! for big boxes we are never here. Only for small ...
!$!             !          temp(1:dim_xdesc, 1:3) = temp(1:dim_xdesc, 1:3) + config_desc(iconf)%force(1:dim_xdesc,ik,0,1:3 )
!$!             !debug stop 'for big never ever should be here'
!$!             ! else
!$!             do inn2 = 1, config_desc(iconf)%n_neigh(ia)
!$!                ib = config_desc(iconf)%kind_neigh(ia, inn2)
!$!                if (ib /= ik) cycle
!$!                if (config_real(iconf)%small) then
!$!                   if (SUM((config_real(iconf)%u_ij(ik, inn1, :) + config_real(iconf)%u_ij(ia, inn2, :))**2) .lt. 1.d-6) then
!$!                      temp(1:dim_xdesc, 1:3) = temp(1:dim_xdesc, 1:3) + config_desc(iconf)%force(1:dim_xdesc, ia, inn2, 1:3)
!$!                   end if
!$!                else
!$!                   temp(1:dim_xdesc, 1:3) = temp(1:dim_xdesc, 1:3) + config_desc(iconf)%force(1:dim_xdesc, ia, inn2, 1:3)
!$!                end if
!$!                !icount1 = icount1 + 1
!$!             end do
!$!             !end if
!$! 
!$!          end do
!$!          config_desc(iconf)%pack_force_linear(1:dim_xdesc, 1:3, ik) = -temp(1:dim_xdesc, 1:3)
!$!       end do
!$! 
!$!    end subroutine serial_snap_pack_force_descriptor


end module  mld_force_mod
