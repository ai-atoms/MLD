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

#include "../MLD_MACROS.INC"

module mld_energy_mod

   use mld_logger

contains

   subroutine pack_energy_descriptor(iconf)

      use derived_types, only: config_real
      use ml_in_ndm_module, only:  mld_order, &
         mld_linear, mld_linear_extended, &
         mld_quadratic, mld_polyc, mld_kernel
      use module_kernel_2b, only: activate_k2b
      implicit none
      integer, intent(in)  :: iconf

      _NAMECURRENT_("pack_energy_descriptor")


      _MLD_BEGIN_

      if (activate_k2b) then
         call kernel_2b_pack_energy(iconf)
      end if

      if ((mld_order == mld_linear) .or. (mld_order == mld_linear_extended)) then
         if (.not. (config_real(iconf)%has_energy)) return
      end if


      if ((mld_order == mld_linear) .or. (mld_order == mld_linear_extended)) then
         call linear_pack_energy_descriptor(iconf)
      end if

      if (mld_order == mld_quadratic) then
           !$! call linear_pack_energy_descriptor(iconf)
           call quadratic_pack_energy_descriptor(iconf)   
      end if


      if (mld_order == mld_kernel) then
         call kernel_pack_energy_descriptor(iconf)
      end if

      if (mld_order == mld_polyc) then
         call polyc_pack_energy_descriptor(iconf)
      end if



      _MLD_END_

   end subroutine pack_energy_descriptor




   subroutine linear_pack_energy_descriptor(iconf)

      use temporary_data_cov, only: dim_xdesc
      use derived_types, only: config_desc, config_real
      !use mpi

      implicit none
      integer, intent(in)  :: iconf
      !$! integer :: id1, ik 

      _NAMECURRENT_("linear_pack_energy_descriptor")

      _MLD_BEGIN_

      if (allocated(config_desc(iconf)%pack_energy_linear)) deallocate (config_desc(iconf)%pack_energy_linear); allocate (config_desc(iconf)%pack_energy_linear(dim_xdesc))
      config_desc(iconf)%pack_energy_linear(1:dim_xdesc) = SUM(config_desc(iconf)%energy(1:dim_xdesc, 1:config_real(iconf)%nat), dim=2)

      !$! config_desc(iconf)%pack_energy_linear(1:dim_xdesc) = 0.d0 
      !$! do ik = 1, config_real(iconf)%nat
      !$!    do id1 = 1, dim_xdesc
      !$!    config_desc(iconf)%pack_energy_linear(id1) = config_desc(iconf)%pack_energy_linear(id1) +  config_desc(iconf)%energy(id1, ik)
      !$!    end do 
      !$! end do 
      ! No need for all_reduce among the subworlds as all the processes in the subworld already has an up-to-date config_desc(iconf)%energy
      ! (all_reduce done at the end of computing the descriptors)
      _MLD_END_
   end subroutine linear_pack_energy_descriptor




   subroutine kernel_pack_energy_descriptor(iconf)

      use temporary_data_cov, only: dim_xdesc
      use derived_types, only: config_desc, config_real
      use module_kernel, only: dim_kernel
      use module_mld_kernel, only: dim_xdesc_kernel
      implicit none
      integer, intent(in)  :: iconf

      _NAMECURRENT_("kernel_pack_energy_descriptor")


      _MLD_BEGIN_
      if (allocated(config_desc(iconf)%pack_energy_linear)) deallocate (config_desc(iconf)%pack_energy_linear); allocate (config_desc(iconf)%pack_energy_linear(dim_xdesc))
      if (allocated(config_desc(iconf)%pack_energy_kernel)) deallocate (config_desc(iconf)%pack_energy_kernel); allocate (config_desc(iconf)%pack_energy_kernel(dim_xdesc_kernel - 1))

      config_desc(iconf)%pack_energy_linear(1:dim_xdesc) = SUM(config_desc(iconf)%energy(1:dim_xdesc, 1:config_real(iconf)%nat), dim=2)
      config_desc(iconf)%pack_energy_kernel(dim_xdesc + 1:dim_xdesc + dim_kernel) = SUM(config_desc(iconf)%energy_kernel(1:dim_kernel, 1:config_real(iconf)%nat), dim=2)
      config_desc(iconf)%pack_energy_kernel(1:dim_xdesc) = config_desc(iconf)%pack_energy_linear(1:dim_xdesc)

      if (allocated(config_desc(iconf)%energy_kernel)) deallocate (config_desc(iconf)%energy_kernel)
      _MLD_END_
   end subroutine kernel_pack_energy_descriptor

   subroutine kernel_2b_pack_energy(iconf)

      use derived_types, only: config_desc, config_real
      use module_kernel_2b, only: dim_kernel_2b
      implicit none
      integer, intent(in)  :: iconf

      _NAMECURRENT_("kernel_2b_pack_energy")


      _MLD_BEGIN_
      if (allocated(config_desc(iconf)%pack_energy_k2b)) deallocate (config_desc(iconf)%pack_energy_k2b)
      allocate (config_desc(iconf)%pack_energy_k2b(dim_kernel_2b))

      config_desc(iconf)%pack_energy_k2b(1:dim_kernel_2b) = SUM(config_desc(iconf)%energy_k2b(1:dim_kernel_2b, 1:config_real(iconf)%nat), dim=2)

      !AD if (allocated(config_desc(iconf)%energy_k2b)) deallocate (config_desc(iconf)%energy_k2b)
      _MLD_END_
   end subroutine kernel_2b_pack_energy




   ! this can be done PARA_UPDATE
   subroutine quadratic_pack_energy_descriptor(iconf)

      ! The energy is packed in the vectors pack_energy and pack_energy_quadratic
      !       the pack_energy_quadratic with the following structure
      !        |          ... pack_energy_quadratic ...  |
      !        | ...  dim_xdesc  ... | ... dim_xdesc*dim_xdesc ... |
      ! For compatibility reason we keep also pack_energy
      !            pack_energy =  | ...  dim_xdesc  ... |
      use module_kind_variables, only: kind_double
      use ml_in_ndm_module, only: mld_type_quadratic, mld_type_quadratic_zaxa, mld_type_quadratic_ZX
      use temporary_data_cov, only: dim_xdesc
      use derived_types, only: config_desc, config_real
      use module_mld_quadratic, only: dim_xdesc_quadratic
      use module_kernel_2b, only: dim_kernel_2b
      implicit none
      integer, intent(in)  :: iconf
      real(kind=kind(0.d0)), dimension(dim_xdesc, 1)     :: Dtmp
      real(kind=kind(0.d0)), dimension(dim_kernel_2b, 1) :: Ztmp
      real(kind_double) :: Dnorm, Znorm, temp 
      integer  :: ik, id1, id2, indx, INCX

      _NAMECURRENT_("quadratic_pack_energy_descriptor")


      _MLD_BEGIN_
      INCX=1 
      if (allocated(config_desc(iconf)%pack_energy_linear)) deallocate (config_desc(iconf)%pack_energy_linear); allocate (config_desc(iconf)%pack_energy_linear(dim_xdesc))
      if (allocated(config_desc(iconf)%pack_energy_quadratic)) deallocate (config_desc(iconf)%pack_energy_quadratic); allocate (config_desc(iconf)%pack_energy_quadratic(dim_xdesc_quadratic - 1))

      !if (allocated(config_desc(iconf)%energy_quadratic ))     deallocate(config_desc(iconf)%energy_quadratic ); allocate(config_desc(iconf)%energy_quadratic(dim_xdesc_quadratic-dim_xdesc-1, config_real(iconf)%nat))
      !The quadratic pack_energy_quadratic  part ....
      !     | ... dim_xdesc ... | ... pack_energy_quadratic(dim_xdesc*dim_xdesc) ... |
      config_desc(iconf)%pack_energy_quadratic(:) = 0.d0

      ! The linear pack_energy part .... for compatibility reasons ...
      !     | ... pack_energy(dim_xdesc) ... |   
      config_desc(iconf)%pack_energy_linear(1:dim_xdesc) = SUM(config_desc(iconf)%energy(1:dim_xdesc, 1:config_real(iconf)%nat), dim=2)
      config_desc(iconf)%pack_energy_quadratic(1:dim_xdesc) = config_desc(iconf)%pack_energy_linear(1:dim_xdesc)

      if (mld_type_quadratic == mld_type_quadratic_ZX) then 
          Dtmp(1:dim_xdesc, 1) = config_desc(iconf)%pack_energy_linear(1:dim_xdesc)
          Ztmp(1:dim_kernel_2b, 1) = config_desc(iconf)%pack_energy_k2b(1:dim_kernel_2b)
          !$! Dnorm = dnrm2 (dim_xdesc,  Dtmp(1:dim_xdesc, 1), INCX)
          !$! Znorm = dnrm2 (dim_kernel_2b, Ztmp(1:dim_kernel_2b, 1), INCX)
          Dnorm = dsqrt ( DOT_PRODUCT (Dtmp(1:dim_xdesc, 1), Dtmp(1:dim_xdesc, 1)))
          Znorm = dsqrt ( DOT_PRODUCT (Ztmp(1:dim_kernel_2b,1),Ztmp(1:dim_kernel_2b,1)))
          temp = (1.d0/Dnorm + 1.d0/Znorm)
          
          config_desc(iconf)%pack_energy_quadratic(1:dim_xdesc) = config_desc(iconf)%pack_energy_linear(1:dim_xdesc)
          do id1 = 1, dim_xdesc
             do id2 = 1, dim_kernel_2b
                indx = dim_kernel_2b*(id1 - 1) + id2 - 1 + dim_xdesc + 1
                config_desc(iconf)%pack_energy_quadratic(indx) = config_desc(iconf)%pack_energy_quadratic(indx) + Dtmp(id1, 1)*Ztmp(id2, 1)*temp
             end do
          end do
         
      else if (mld_type_quadratic == mld_type_quadratic_zaxa ) then

        do ik = 1, config_real(iconf)%nat
          Dtmp(1:dim_xdesc, 1) = config_desc(iconf)%energy(1:dim_xdesc, ik)
          Ztmp(1:dim_kernel_2b, 1) = config_desc(iconf)%energy_k2b(1:dim_kernel_2b, ik)
          !$! do id1 = 1, dim_xdesc
          !$!    indx = id1
          !$!    config_desc(iconf)%pack_energy_quadratic(indx) = config_desc(iconf)%pack_energy_quadratic(indx) + Dtmp(id1, 1)
          !$! end do
          do id1 = 1, dim_xdesc
             do id2 = 1, dim_kernel_2b
                indx = dim_kernel_2b*(id1 - 1) + id2 - 1 + dim_xdesc + 1
                config_desc(iconf)%pack_energy_quadratic(indx) = config_desc(iconf)%pack_energy_quadratic(indx) + Dtmp(id1, 1)*Ztmp(id2, 1)
             end do
          end do
        end do  

      else 


        do ik = 1, config_real(iconf)%nat
           Dtmp(1:dim_xdesc, 1) = config_desc(iconf)%energy(1:dim_xdesc, ik)
           !$! do id1 = 1, dim_xdesc
           !$!    indx = id1
           !$!    config_desc(iconf)%pack_energy_quadratic(indx) = config_desc(iconf)%pack_energy_quadratic(indx) + Dtmp(id1, 1)
           !$! end do
           do id1 = 1, dim_xdesc
              do id2 = 1, dim_xdesc
                 indx = dim_xdesc*(id1 - 1) + id2 - 1 + dim_xdesc + 1
                 config_desc(iconf)%pack_energy_quadratic(indx) = config_desc(iconf)%pack_energy_quadratic(indx) + Dtmp(id1, 1)*Dtmp(id2, 1)
              end do
           end do
  
        end do
      end if 
      ! The linear pack_energy part .... for compatibility reasons ...
      !     | ... pack_energy(dim_xdesc) ... |
      !$!  config_desc(iconf)%pack_energy_linear(1:dim_xdesc) = config_desc(iconf)%pack_energy_quadratic(1:dim_xdesc)

      _MLD_END_
   end subroutine quadratic_pack_energy_descriptor




   ! this can be done PARA_UPDATE
   subroutine polyc_pack_energy_descriptor(iconf)

      ! The energy is packed in the vectors pack_energy and pack_energy_quadratic
      !       the pack_energy_quadratic with the following structure
      !        |          ... pack_energy_quadratic ...  |
      !        | ...  dim_xdesc  ... | ... dim_xdesc*dim_xdesc ... |
      ! For compatibility reason we keep also pack_energy
      !            pack_energy =  | ...  dim_xdesc  ... |

      use ml_in_ndm_module, only:  polyc_n_poly, polyc_n_hermite
      use temporary_data_cov, only: dim_xdesc
      use derived_types, only: config_desc, config_real
      use module_mld_polyc, only: dim_xdesc_polyc
      use math, only: hermite

      implicit none

      integer, intent(in)  :: iconf
      real(kind=kind(0.d0)), dimension(polyc_n_hermite, dim_xdesc, 1)      :: Dtmp
      real(kind=kind(0.d0))      :: y_her, d_y_her
      integer  :: ik, id1, id2, id3, indx, nord

      _NAMECURRENT_("polyc_pack_energy_descriptor")


      _MLD_BEGIN_
      if (allocated(config_desc(iconf)%pack_energy_linear)) deallocate (config_desc(iconf)%pack_energy_linear); allocate (config_desc(iconf)%pack_energy_linear(dim_xdesc))
      if (allocated(config_desc(iconf)%pack_energy_polyc)) deallocate (config_desc(iconf)%pack_energy_polyc); allocate (config_desc(iconf)%pack_energy_polyc(dim_xdesc_polyc - 1))

      !The linear pack_energy part ....
      !     | ... pack_energy(dim_xdesc) ... | ... dim_xdesc*dim_xdesc ... |
      !config_desc(iconf)%pack_energy(1:dim_xdesc)=SUM(config_desc(iconf)%energy(1:dim_xdesc,1:config_real(iconf)%nat), dim=2)

      !The  energy_polyc  part ....
      !     | ... dim_xdesc ... | ... pack_energy_polyc(dim_xdesc_polyc) ... |
      config_desc(iconf)%pack_energy_polyc(:) = 0.d0
      do ik = 1, config_real(iconf)%nat

         do nord = 1, polyc_n_hermite
            do id1 = 1, dim_xdesc
               call hermite(nord, config_desc(iconf)%energy(id1, ik), y_her, d_y_her)
               Dtmp(nord, id1, 1) = y_her
            end do
         end do


         indx = 0
         !........... polyc_n_poly = 1
         if (polyc_n_poly >= 1) then
            do nord = 1, polyc_n_hermite
               do id1 = 1, dim_xdesc
                  indx = indx + 1
                  config_desc(iconf)%pack_energy_polyc(indx) = config_desc(iconf)%pack_energy_polyc(indx) + Dtmp(nord, id1, 1)
               end do
            end do
         end if

         !........... polyc_n_poly = 2
         if (polyc_n_poly >= 2) then
            do nord = 1, polyc_n_hermite
               do id1 = 1, dim_xdesc
                  do id2 = 1, dim_xdesc
                     indx = indx + 1
                     config_desc(iconf)%pack_energy_polyc(indx) = config_desc(iconf)%pack_energy_polyc(indx) + Dtmp(nord, id1, 1)*Dtmp(nord, id2, 1)
                  end do
               end do
            end do
         end if


         !........... polyc_n_poly = 3
         if (polyc_n_poly >= 3) then
            do nord = 1, polyc_n_hermite
               do id1 = 1, dim_xdesc
                  do id2 = 1, dim_xdesc
                     do id3 = 1, dim_xdesc
                        indx = indx + 1
                        config_desc(iconf)%pack_energy_polyc(indx) = config_desc(iconf)%pack_energy_polyc(indx) + Dtmp(nord, id1, 1)*Dtmp(nord, id2, 1)*Dtmp(nord, id3, 1)
                     end do
                  end do
               end do
            end do
         end if

      end do                  ! end loop over ik the atom index.

      config_desc(iconf)%pack_energy_linear(1:dim_xdesc) = config_desc(iconf)%pack_energy_polyc(1:dim_xdesc)
      _MLD_END_
   end subroutine polyc_pack_energy_descriptor



   subroutine train_fill_Amat_with_energy(iconf, pack_opt)

      use module_kind_variables, only: kind_double
      use ml_in_ndm_module, only:  mld_order, &
         mld_linear, mld_linear_extended, &
         mld_quadratic, mld_polyc, mld_kernel
      use temporary_data_cov, only: yfunc_train, dim_xdesc, dim_xdesc_patch
      use derived_types, only: config_desc, config_real

      use snap, only: ymat, Amat, i_fit_snap, fit_snap, weights_snap, dim_xdesc_linear, &
                      dim_design_line, y_zbl_train
      use module_mld_quadratic, only: dim_xdesc_quadratic
      use module_mld_polyc, only: dim_xdesc_polyc
      use module_mld_kernel, only: dim_xdesc_kernel
      use module_db_poscar, only: iread_energy
      use module_init_ScaMatrix, only: init_ScaMatrix_with_LocalVector
      use module_kernel_2b, only: dim_kernel_2b, activate_k2b
      use module_zbl, only: zbl_potential, zbl_type, zbl_mode_alone
      implicit none

      integer, intent(in)  :: iconf
      integer  :: i_fit_snap_end, last_dim
      ! integer :: i
      ! real(kind(0.d0)) :: tmp(1:dim_xdesc)
      logical, intent(in), optional    :: pack_opt
      logical  :: lpack_tmp
      real(kind_double), dimension(:), allocatable :: local_vec

      _NAMECURRENT_("train_fill_Amat_with_energy")


      _MLD_BEGIN_
      lpack_tmp = .true.
      if (present(pack_opt)) lpack_tmp = pack_opt

      if (.not. (config_real(iconf)%has_energy)) return
      i_fit_snap = i_fit_snap + 1
      i_fit_snap_end = i_fit_snap
      if (lpack_tmp) call pack_energy_descriptor(iconf)
      ! set y and the characteristic of the fit for this line:
      yfunc_train(i_fit_snap) = config_real(iconf)%energy(iread_energy)
      fit_snap(i_fit_snap)%iconf = iconf
      fit_snap(i_fit_snap)%energy = .true.
      fit_snap(i_fit_snap)%weight = config_real(iconf)%w_e_model
      weights_snap(i_fit_snap) = config_real(iconf)%w_e_model

      last_dim = dim_design_line

      if ((mld_order == mld_linear) .or. (mld_order == mld_linear_extended)) then
         !Amat = | N | ... dim_xdesc       ... |    ...  dim_xdesc_patch ... |  -> i_fit_snap
         !last_dim = dim_xdesc_linear + dim_xdesc_patch
         if (allocated(local_vec)) deallocate(local_vec) ; allocate(local_vec(last_dim))
         local_vec(2 : dim_xdesc_linear) = config_desc(iconf)%pack_energy_linear(1:dim_xdesc_linear - 1)
      end if

      if (mld_order == mld_quadratic) then
         !Amat = | N | ... dim_xdesc ... | ...  dim_xdesc_patch ... | ... dim_xdesc*dim_xdesc ... | -> i_fit_snap
         !last_dim = dim_xdesc_quadratic + dim_xdesc_patch
         if (allocated(local_vec)) deallocate(local_vec) ; allocate(local_vec(last_dim))
         local_vec(2 : dim_xdesc_linear) = config_desc(iconf)%pack_energy_linear(1:dim_xdesc_linear - 1)
         local_vec(dim_xdesc_linear + dim_xdesc_patch + 1 : last_dim) = config_desc(iconf)%pack_energy_quadratic(dim_xdesc_linear : dim_xdesc_quadratic - 1)
      end if

      if (mld_order == mld_polyc) then
         !Amat = | N | ... dim_xdesc ... | ...  dim_xdesc_patch ... |  ...  dim_xdesc_polyc  ...  | -> i_fit_snap
         !last_dim = dim_xdesc_polyc + dim_xdesc_patch
         if (allocated(local_vec)) deallocate(local_vec) ; allocate(local_vec(last_dim))
         local_vec(2 : dim_xdesc_linear) = config_desc(iconf)%pack_energy_linear(1:dim_xdesc_linear - 1)
         local_vec(dim_xdesc_linear + dim_xdesc_patch + 1 : last_dim) = config_desc(iconf)%pack_energy_polyc(dim_xdesc_linear : dim_xdesc_polyc - 1)
      end if

      if (mld_order == mld_kernel) then
         !Amat = | N | ... dim_xdesc ... |     ...  dim_xdesc_patch ... | ... dim_xdesc_kernel ... | -> i_fit_snap
         !last_dim = dim_xdesc_kernel + dim_xdesc_patch
         if (allocated(local_vec)) deallocate(local_vec) ; allocate(local_vec(last_dim))
         local_vec(2 : dim_xdesc_linear) = config_desc(iconf)%pack_energy_kernel(1 : dim_xdesc_linear - 1)
         local_vec(dim_xdesc_linear + dim_xdesc_patch + 1 : last_dim) = config_desc(iconf)%pack_energy_kernel(dim_xdesc_linear : dim_xdesc_kernel - 1)
      end if

      if (activate_k2b) then
         local_vec(dim_xdesc_linear + 1: dim_xdesc_linear + dim_kernel_2b) = config_desc(iconf)%pack_energy_k2b(1: dim_kernel_2b)
      end if

      local_vec(1) = dble(config_real(iconf)%nat)

      Amat(1 : last_dim, i_fit_snap ) = local_vec (1 : last_dim)

      if (zbl_potential .and. zbl_type == zbl_mode_alone) then
         ymat(i_fit_snap,1) = yfunc_train(i_fit_snap) - config_real(iconf)%ezbl
         y_zbl_train(i_fit_snap) = config_real(iconf)%ezbl
      else
         ymat(i_fit_snap, 1) = yfunc_train(i_fit_snap)
      end if

      i_fit_snap = i_fit_snap_end

      _MLD_END_
   end subroutine train_fill_Amat_with_energy



   subroutine train_fill_Bmat_constraints(iconf)

      use temporary_data_cov, only: dim_xdesc
      use derived_types, only: config_desc, config_real
      use snap, only: Bmat, zmat, i_constraints_snap
      use module_db_poscar, only: iread_energy
      implicit none
      integer, intent(in)  :: iconf

      _NAMECURRENT_("train_fill_Bmat_constraints")


      _MLD_BEGIN_
      ! integer :: i
      ! real(kind(0.d0)) :: tmp(1:dim_xdesc)

      if (.not. (config_real(iconf)%has_energy)) return
      i_constraints_snap = i_constraints_snap + 1
      Bmat(1, i_constraints_snap) = dble(config_real(iconf)%nat)
      Bmat(2:dim_xdesc + 1, i_constraints_snap) = SUM(config_desc(iconf)%energy(1:dim_xdesc, 1:config_real(iconf)%nat), dim=2)
      zmat(i_constraints_snap, 1) = config_real(iconf)%energy(iread_energy)
      _MLD_END_
   end subroutine train_fill_Bmat_constraints


end module
