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

module mld_stress_mod
#if(PARA)
   use mpi
   use mld_mpi 
#endif

   use mld_logger

   implicit none

contains


   subroutine pack_stress_descriptor(iconf)
      use ml_in_ndm_module, only: descriptor_type, &
         descriptor_g3, &
         descriptor_g2, &
         descriptor_behler, &
         descriptor_afs, &
         descriptor_bispectrum_so4, &
         descriptor_g2_bispectrum_so4, &
         descriptor_g2_pow_so4, &
         descriptor_g2_afs, &
         descriptor_mtp, &
         descriptor_pow_so3, &
         descriptor_pow_so4, &
         descriptor_milady, descriptor_body, descriptor_zetabody, &
         descriptor_pow_so3_3body, &
         descriptor_ace, descriptor_ftnbody, &
         descriptor_tbind, &
         mld_order, mld_linear, mld_quadratic, &
         mld_polyc, mld_kernel
      use temporary_data_cov, only: dim_xdesc
      use derived_types, only: config_desc
      use module_mld_quadratic, only: dim_xdesc_quadratic
      use module_mld_polyc, only: dim_xdesc_polyc
      use module_mld_kernel, only: dim_xdesc_kernel
      use module_kernel_2b, only: activate_k2b, dim_kernel_2b
      use mld_subworld, only : subworld

      integer, intent(in)  :: iconf
      integer  :: dim_reduce

      _NAMECURRENT_("snap_pack_stress_descriptor")



      _MLD_BEGIN_
      if ((descriptor_type == descriptor_afs) .or. &
         (descriptor_type == descriptor_g3) .or. &
         (descriptor_type == descriptor_g2) .or. &
         (descriptor_type == descriptor_behler) .or. &
         (descriptor_type == descriptor_bispectrum_so4) .or. &
         (descriptor_type == descriptor_g2_bispectrum_so4) .or. &
         (descriptor_type == descriptor_g2_pow_so4) .or. &
         (descriptor_type == descriptor_g2_afs) .or. &
         (descriptor_type == descriptor_pow_so3) .or. &
         (descriptor_type == descriptor_pow_so3_3body) .or. &
         (descriptor_type == descriptor_pow_so4) .or. &
         (descriptor_type == descriptor_milady) .or. &
         (descriptor_type == descriptor_mtp) .or. &
         (descriptor_type == descriptor_ace) .or. &
         (descriptor_type == descriptor_ftnbody) .or. &
         (descriptor_type == descriptor_zetabody) .or. &
         (descriptor_type == descriptor_tbind) .or. &
         (descriptor_type == descriptor_body)) then

         if (mld_order >= mld_linear) then
            if (allocated(config_desc(iconf)%pack_stress_linear)) deallocate (config_desc(iconf)%pack_stress_linear); allocate (config_desc(iconf)%pack_stress_linear(dim_xdesc, 6))
            config_desc(iconf)%pack_stress_linear(:, :) = 0.d0
         end if

         if (mld_order == mld_quadratic) then
            if (allocated(config_desc(iconf)%pack_stress_linear)) deallocate (config_desc(iconf)%pack_stress_linear); allocate (config_desc(iconf)%pack_stress_linear(dim_xdesc, 6))
            if (allocated(config_desc(iconf)%pack_stress_quadratic)) deallocate (config_desc(iconf)%pack_stress_quadratic); allocate (config_desc(iconf)%pack_stress_quadratic(dim_xdesc_quadratic - 1, 6))
            config_desc(iconf)%pack_stress_linear(:, :) = 0.d0
            config_desc(iconf)%pack_stress_quadratic(:, :) = 0.d0
         end if


         if (mld_order == mld_polyc) then
            if (allocated(config_desc(iconf)%pack_stress_linear)) deallocate (config_desc(iconf)%pack_stress_linear); allocate (config_desc(iconf)%pack_stress_linear(dim_xdesc, 6))
            if (allocated(config_desc(iconf)%pack_stress_polyc)) deallocate (config_desc(iconf)%pack_stress_polyc); allocate (config_desc(iconf)%pack_stress_polyc(dim_xdesc_polyc - 1, 6))
            config_desc(iconf)%pack_stress_linear(:, :) = 0.d0
            config_desc(iconf)%pack_stress_polyc(:, :) = 0.d0
         end if

         if (mld_order == mld_kernel) then
            if (allocated(config_desc(iconf)%pack_stress_linear)) deallocate (config_desc(iconf)%pack_stress_linear); allocate (config_desc(iconf)%pack_stress_linear(dim_xdesc, 6))
            if (allocated(config_desc(iconf)%pack_stress_kernel)) deallocate (config_desc(iconf)%pack_stress_kernel); allocate (config_desc(iconf)%pack_stress_kernel(dim_xdesc_kernel - 1, 6))
            config_desc(iconf)%pack_stress_linear(:, :) = 0.d0
            config_desc(iconf)%pack_stress_kernel(:, :) = 0.d0
         end if

         if (activate_k2b) then
            if (allocated(config_desc(iconf)%pack_stress_k2b)) deallocate (config_desc(iconf)%pack_stress_k2b); allocate (config_desc(iconf)%pack_stress_k2b(dim_kernel_2b, 6))
            config_desc(iconf)%pack_stress_k2b(:, :) = 0.d0

         end if

         call para_pack_stress_descriptor(iconf)
         !... old version ... call para_snap_pack_stress_descriptor(iconf)

#if(PARA)
         !if (mld_order >= mld_linear) then
         dim_reduce = dim_xdesc*6
         call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(iconf)%pack_stress_linear, dim_reduce, MPI_DOUBLE_PRECISION, MPI_SUM, subworld, codeml)
         !end if

         if (mld_order == mld_quadratic) then
            dim_reduce = (dim_xdesc_quadratic - 1)*6
            call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(iconf)%pack_stress_quadratic, dim_reduce, MPI_DOUBLE_PRECISION, MPI_SUM, subworld, codeml)
         end if


         if (mld_order == mld_polyc) then
            dim_reduce = (dim_xdesc_polyc - 1)*6
            call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(iconf)%pack_stress_polyc, dim_reduce, MPI_DOUBLE_PRECISION, MPI_SUM, subworld, codeml)
         end if

         if (mld_order == mld_kernel) then
            dim_reduce = (dim_xdesc_kernel - 1)*6
            call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(iconf)%pack_stress_kernel, dim_reduce, MPI_DOUBLE_PRECISION, MPI_SUM, subworld, codeml)
         end if

         if (activate_k2b) then
            dim_reduce = dim_kernel_2b*6
            call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(iconf)%pack_stress_k2b, dim_reduce, MPI_DOUBLE_PRECISION, MPI_SUM, subworld, codeml)
         end if


#endif

      else
         write (*, *) 'message pron snap_pack_stress_descriptor    : serial_snap_q_pack_stress_descriptor NO LONGER  used. STOP'
         stop ' deasctivate serial snap_q_pack_stress'
         !call serial_snap_q_pack_stress_descriptor(iconf)
         !... old version ... call serial_snap_pack_stress_descriptor(iconf)
      end if


      _MLD_END_

   end subroutine pack_stress_descriptor



   subroutine train_fill_Amat_with_stress(iconf, pack_opt)
      use module_kind_variables, only: kind_double
      use ml_in_ndm_module, only: mld_order, &
         mld_linear, mld_linear_extended, &
         mld_quadratic, mld_polyc, mld_kernel, desc_forces
      use temporary_data_cov, only: yfunc_train, dim_xdesc, dim_xdesc_patch
      use derived_types, only: config_desc, config_real
      use snap, only: Amat, ymat, i_fit_snap, fit_snap, weights_snap, dim_design_line, dim_xdesc_linear, y_zbl_train
      use module_mld_quadratic, only: dim_xdesc_quadratic
      use module_mld_polyc, only: dim_xdesc_polyc
      use module_mld_kernel, only: dim_xdesc_kernel
      use module_init_ScaMatrix, only: init_ScaMatrix_with_LocalMatrix
      use module_kernel_2b, only: activate_k2b, dim_kernel_2b
      use module_zbl, only: zbl_potential, zbl_type, zbl_mode_alone

      integer, intent(in)  :: iconf
      logical, intent(in), optional    :: pack_opt
      integer  :: i_fit_mld_ini, i_fit_mld_end, ix
      logical  :: lpack_tmp
      real(kind_double), dimension(:,:), allocatable :: local_mat

      _NAMECURRENT_("train_fill_Amat_with_stress")



      _MLD_BEGIN_
      if (.not. desc_forces) return
      lpack_tmp = .true.
      if (present(pack_opt)) lpack_tmp = pack_opt
      if (.not. (config_real(iconf)%has_stress)) return
      i_fit_mld_ini = i_fit_snap
      i_fit_mld_end = i_fit_snap + 6
      if (lpack_tmp) call pack_stress_descriptor(iconf)

      yfunc_train(i_fit_mld_ini + 1:i_fit_mld_ini + 6) = config_real(iconf)%stress(1:6)
      weights_snap(i_fit_mld_ini + 1:i_fit_mld_ini + 6) = config_real(iconf)%w_s_model

      fit_snap(i_fit_mld_ini + 1:i_fit_mld_ini + 6)%iconf = iconf
      fit_snap(i_fit_mld_ini + 1:i_fit_mld_ini + 6)%stress = .true.
      fit_snap(i_fit_mld_ini + 1:i_fit_mld_ini + 6)%weight = config_real(iconf)%w_s_model
      do ix = 1, 6
         fit_snap(i_fit_mld_ini + ix)%ix = ix
      end do


      if (allocated(local_mat)) deallocate(local_mat) ; allocate(local_mat(dim_design_line, 6))

      if ((mld_order == mld_linear) .or. (mld_order == mld_linear_extended)) then
         local_mat(1, 1:6) = 0.d0
         local_mat(2:dim_xdesc_linear, 1:6) = config_desc(iconf)%pack_stress_linear(1:dim_xdesc_linear-1, 1:6)
      end if

      if (mld_order == mld_quadratic) then
         local_mat(1, 1:6) = 0.d0
         local_mat(2:dim_xdesc_linear, 1:6) = config_desc(iconf)%pack_stress_linear(1:dim_xdesc_linear - 1, 1:6)
         local_mat(dim_xdesc_linear + dim_xdesc_patch +1:dim_design_line, 1:6) = config_desc(iconf)%pack_stress_quadratic(dim_xdesc_linear:dim_xdesc_quadratic - 1, 1:6)
      end if

      if (mld_order == mld_polyc) then
         local_mat(1, 1:6) = 0.d0
         local_mat(2:dim_xdesc_linear, 1:6) = config_desc(iconf)%pack_stress_linear(1:dim_xdesc_linear - 1, 1:6)
         local_mat(dim_xdesc_linear + dim_xdesc_patch +1:dim_design_line, 1:6) = config_desc(iconf)%pack_stress_polyc(dim_xdesc_linear:dim_xdesc_polyc - 1, 1:6)
      end if


      if (mld_order == mld_kernel) then
         local_mat(1, 1:6) = 0.d0
         local_mat(2:dim_xdesc_linear, 1:6) = config_desc(iconf)%pack_stress_linear(1:dim_xdesc_linear - 1, 1:6)
         local_mat(dim_xdesc_linear + dim_xdesc_patch +1:dim_design_line, 1:6) = config_desc(iconf)%pack_stress_kernel(dim_xdesc_linear:dim_xdesc_kernel - 1, 1:6)
      end if

      if (activate_k2b) then
         local_mat(dim_xdesc_linear + 1: dim_xdesc_linear + dim_kernel_2b, 1:6) = config_desc(iconf)%pack_stress_k2b(1: dim_kernel_2b, 1:6)
      end if

      Amat(1:dim_design_line, i_fit_mld_ini + 1:i_fit_mld_ini + 6) = local_mat(1:dim_design_line,1:6)
      if (zbl_potential .and. zbl_type == zbl_mode_alone) then
         ymat(i_fit_mld_ini + 1:i_fit_mld_end, 1) = yfunc_train(i_fit_mld_ini + 1:i_fit_mld_end) - config_real(iconf)%szbl(1:6)
         y_zbl_train(i_fit_mld_ini + 1:i_fit_mld_end) = config_real(iconf)%szbl(1:6)
      else
         ymat(i_fit_mld_ini + 1:i_fit_mld_end, 1) = yfunc_train(i_fit_mld_ini + 1:i_fit_mld_end)
      end if

      i_fit_snap = i_fit_mld_end

      _MLD_END_
      return
   end subroutine train_fill_Amat_with_stress



   subroutine para_pack_stress_descriptor(iconf)
      ! This subroutine has strong interaction with neighbours subroutines
      ! The one proposed by Milady and the one proposed by NDM depending if the cell box is small or large.
      ! small=.true. ->  MiLaDy
      ! small=.false. -> NDM

#ifdef MLD_NDM
      use gen_com_m, only: lperiod, A2cm
      use gen_com_m_ml, only: imm, bg, at
      use tab_imm_m_ml, only: xp
#else
      use ondm_gen_com_m, only: imm, lperiod, indi2, A2cm, bg, at
      use ondm_tab_imm_m, only: xp, iwmax2
#endif

      !MiLaDy_interaction
      use ml_in_ndm_module, only: sign_stress, sign_stress_big_box, &
         i_start_at, i_final_at, desc_forces, &
         mld_order, mld_linear, mld_linear_extended, &
         mld_quadratic, mld_polyc, mld_kernel, &
         polyc_n_poly, polyc_n_hermite,  & 
         mld_type_quadratic, mld_type_quadratic_zaxa, mld_type_quadratic_ZX
      use module_kernel, only: dim_kernel
      use derived_types, only: config_desc, config_real
      use temporary_data_cov, only: dim_xdesc
      use module_mld_quadratic, only: dim_xdesc_quadratic
      use module_mld_polyc, only: dim_xdesc_polyc
      use module_mld_kernel, only: dim_xdesc_kernel
      use module_neigh_local, only: r_cut
      use module_kernel_2b, only: dim_kernel_2b, activate_k2b
      use math, only: hermite
#ifdef MLD_NDM
      use notperiod_mod
#else
      use ondm_transform_coord, only: ondm_notperiod
#endif
      integer, intent(in)  :: iconf
      integer  :: ia, ic
      real(kind(0.d0))     :: rr, uu(3), ds(3)
      real(kind(0.d0)), dimension(:, :), allocatable     :: force_ij, st_temp, &
         force_ij_quadratic, st_temp_quadratic, &
         force_ij_polyc, st_temp_polyc, &
         force_ij_kernel, st_temp_kernel, &
         force_ij_k2b,  st_temp_k2b
      integer  :: id1, id2, id3, indx_last, indx, ik, nord
      integer  :: iw, iw1, iw2, ia_n

      real(kind(0.d0)), dimension(:, :), allocatable     :: xpnp
      logical  :: small
      real(kind=kind(0.d0))      :: y_her, d_y_her
      real(kind=kind(0.d0)), dimension(:, :), allocatable      :: etmp, d_etmp

      _NAMECURRENT_("para_snap_pack_stress_descriptor")



      _MLD_BEGIN_
      if (.not. desc_forces) return
      if ((i_start_at == 0) .and. (i_final_at == 0)) return
      if (.not. (config_real(iconf)%has_stress)) return

      if (allocated(etmp)) deallocate (etmp); allocate (etmp(polyc_n_hermite, dim_xdesc))
      if (allocated(d_etmp)) deallocate (d_etmp); allocate (d_etmp(polyc_n_hermite, dim_xdesc))

      if ((mld_order == mld_linear) .or. (mld_order == mld_linear_extended)) then
         if (allocated(force_ij)) deallocate (force_ij); allocate (force_ij(dim_xdesc, 3))
         if (allocated(st_temp)) deallocate (st_temp); allocate (st_temp(dim_xdesc, 6))
      end if

      if (mld_order == mld_quadratic) then
         if (allocated(force_ij_quadratic)) deallocate (force_ij_quadratic); allocate (force_ij_quadratic(dim_xdesc_quadratic - 1, 3))
         if (allocated(st_temp_quadratic)) deallocate (st_temp_quadratic); allocate (st_temp_quadratic(dim_xdesc_quadratic - 1, 6))
      end if

      if (mld_order == mld_polyc) then
         if (allocated(force_ij_polyc)) deallocate (force_ij_polyc); allocate (force_ij_polyc(dim_xdesc_polyc - 1, 3))
         if (allocated(st_temp_polyc)) deallocate (st_temp_polyc); allocate (st_temp_polyc(dim_xdesc_polyc - 1, 6))
      end if


      if (mld_order == mld_kernel) then
         if (allocated(force_ij_kernel)) deallocate (force_ij_kernel); allocate (force_ij_kernel(dim_xdesc_kernel - 1, 3))
         if (allocated(st_temp_kernel)) deallocate (st_temp_kernel); allocate (st_temp_kernel(dim_xdesc_kernel - 1, 6))
      end if

      if (activate_k2b) then
         if (allocated(force_ij_k2b)) deallocate (force_ij_k2b); allocate (force_ij_k2b(dim_kernel_2b, 3))
         if (allocated(st_temp_k2b)) deallocate (st_temp_k2b); allocate (st_temp_k2b(dim_kernel_2b, 6))
      end if


      small = config_real(iconf)%small
#ifdef MLD_NDM
!CRC    small=.true.
#endif
      if (allocated(xpnp)) deallocate (xpnp); allocate (xpnp(3, imm))
      if (lperiod) then
         xpnp(:, :) = xp(:, :)
      else
#ifdef MLD_NDM
        call notperiod(imm,xp, xpnp,at,bg,.false.)
#else
         call ondm_notperiod(xp, xpnp)
#endif
      end if
      ! call cryst_to_cart (imm, xpnp, bg, -1)

      iw2 = 0

      ! init stress
      if ((mld_order == mld_linear) .or. (mld_order == mld_linear_extended)) st_temp(:, :) = 0.d0
      if (mld_order == mld_quadratic) st_temp_quadratic(:, :) = 0.d0
      if (mld_order == mld_polyc) st_temp_polyc(:, :) = 0.d0
      if (mld_order == mld_kernel) st_temp_kernel(:, :) = 0.d0
      if (activate_k2b) st_temp_k2b(:, :) = 0.d0
      ! init stress

      ! debug write (*,*) i_start_at, i_final_at, rangml
#ifdef MLD_NDM
    iw2=0
!!$    if (i_start_at == 1) iw2 = 0
!!$    if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
#else
      if (i_start_at == 1) iw2 = 0
      if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
#endif
      do ia = i_start_at, i_final_at                   ! descriptor index
         !begin small box or not 1/
#ifdef MLD_NDM
!!$      if (small) then
#else
         if (small) then
#endif
            iw1 = 1
            iw2 = config_real(iconf)%n_neigh(ia)
#ifdef MLD_NDM
!!$      else
!!$        iw1 = iw2 + 1
!!$        iw2 = iwmax2(ia)
#else
         else
            iw1 = iw2 + 1
            iw2 = iwmax2(ia)
#endif
            !write (*,*) 'debug neigh', iw1, iw2
#ifdef MLD_NDM
!!$      end if
#else
         end if
#endif
         ia_n = 0
         do iw = iw1, iw2
#ifdef MLD_NDM
!CRC        if (small) then
#else
            if (small) then
#endif
               ic = config_real(iconf)%kind_neigh(ia, iw)
#ifdef MLD_NDM
!!$        else
!!$          ic = indi2(iw)
!!$          if (ic == ia) cycle
!CRC        end if
!        write(6,*)'PCRCN1',small
#else
            else
               ic = indi2(iw)
               if (ic == ia) cycle
            end if
#endif

            if (small) then
               rr = config_real(iconf)%r_ij(ia, iw)
               uu(:) = config_real(iconf)%u_ij(ia, iw, :)

               ik = config_real(iconf)%kind_neigh(ia, iw)

               if ((mld_order == mld_linear) .or. (mld_order == mld_linear_extended)) then
                  force_ij(1:dim_xdesc, 1:3) = config_desc(iconf)%force(1:dim_xdesc, ia, iw, 1:3)
               end if


               if (mld_order == mld_kernel) then
                  force_ij_kernel(1:dim_xdesc, 1:3) = config_desc(iconf)%force(1:dim_xdesc, ia, iw, 1:3)
                  force_ij_kernel(dim_xdesc + 1:dim_xdesc_kernel - 1, 1:3) = config_desc(iconf)%force_kernel(1:dim_kernel, ia, iw, 1:3)
               end if                  ! mld_kernel

               if (activate_k2b) then
                  force_ij_k2b(1:dim_kernel_2b, 1:3) = config_desc(iconf)%force_k2b(1:dim_kernel_2b, ia, iw, 1:3)
               end if                  ! mld_kernel


               if (mld_order == mld_quadratic) then
                  if (mld_type_quadratic == mld_type_quadratic_zaxa) then 
                    force_ij_quadratic(1:dim_xdesc, 1:3) = config_desc(iconf)%force(1:dim_xdesc, ia, iw, 1:3)
                    indx_last = 0
                    do id1 = 1, dim_xdesc
                       do id2 = 1, dim_kernel_2b
                          indx = dim_kernel_2b*(id1 - 1) + id2 - 1 + dim_xdesc + 1
                          force_ij_quadratic(indx, 1:3) = (config_desc(iconf)%force(id1, ia, iw, 1:3)*config_desc(iconf)%energy_k2b(id2, ik) + &
                             config_desc(iconf)%force_k2b(id2, ia, iw, 1:3)*config_desc(iconf)%energy(id1, ik))
                       end do
                    end do

                  else if  (mld_type_quadratic == mld_type_quadratic_ZX) then
                       call log_critical("ZX quadratic not yet implemented for stress. Avoid stress fit") 
                       stop 
                  else  
                    force_ij_quadratic(1:dim_xdesc, 1:3) = config_desc(iconf)%force(1:dim_xdesc, ia, iw, 1:3)
                    indx_last = 0
                    do id1 = 1, dim_xdesc
                       do id2 = 1, dim_xdesc
                          indx = dim_xdesc*(id1 - 1) + id2 - 1 + dim_xdesc + 1
                          force_ij_quadratic(indx, 1:3) = (config_desc(iconf)%force(id1, ia, iw, 1:3)*config_desc(iconf)%energy(id2, ik) + &
                             config_desc(iconf)%force(id2, ia, iw, 1:3)*config_desc(iconf)%energy(id1, ik))
                       end do
                    end do
                  end if 
               end if                  ! mld_quadratic

               if (mld_order == mld_polyc) then


                  do nord = 1, polyc_n_hermite
                     do id1 = 1, dim_xdesc
                        call hermite(nord, config_desc(iconf)%energy(id1, ik), y_her, d_y_her)
                        etmp(nord, id1) = y_her
                        d_etmp(nord, id1) = d_y_her
                     end do
                  end do

                  indx_last = 0
                  if (polyc_n_poly >= 1) then
                     indx = indx_last
                     do nord = 1, polyc_n_hermite
                        do id1 = 1, dim_xdesc
                           indx = indx + 1
                           force_ij_polyc(indx, 1:3) = config_desc(iconf)%force(id1, ia, iw, 1:3)*d_etmp(nord, id1)
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
                              !force_ij_polyc(indx,1:3) = (config_desc(iconf)%force(id1,ia,iw,1:3)*config_desc(iconf)%energy(id2,ik) + &
                              !                            config_desc(iconf)%force(id2,ia,iw,1:3)*config_desc(iconf)%energy(id1,ik) )
                              force_ij_polyc(indx, 1:3) = (config_desc(iconf)%force(id1, ia, iw, 1:3)*d_etmp(nord, id1)*etmp(nord, id2) + &
                                 config_desc(iconf)%force(id2, ia, iw, 1:3)*d_etmp(nord, id2)*etmp(nord, id1))
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
                                 !force_ij_polyc(indx,1:3) = (config_desc(iconf)%force(id1,ia,iw,1:3)*config_desc(iconf)%energy(id2,ik)*config_desc(iconf)%energy(id3,ik) + &
                                 !                            config_desc(iconf)%force(id2,ia,iw,1:3)*config_desc(iconf)%energy(id1,ik)*config_desc(iconf)%energy(id3,ik) + &
                                 !                            config_desc(iconf)%force(id3,ia,iw,1:3)*config_desc(iconf)%energy(id1,ik)*config_desc(iconf)%energy(id2,ik) )
                                 force_ij_polyc(indx, 1:3) = (config_desc(iconf)%force(id1, ia, iw, 1:3)*d_etmp(nord, id1)*etmp(nord, id2)*etmp(nord, id3) + &
                                    config_desc(iconf)%force(id2, ia, iw, 1:3)*d_etmp(nord, id2)*etmp(nord, id1)*etmp(nord, id3) + &
                                    config_desc(iconf)%force(id3, ia, iw, 1:3)*d_etmp(nord, id3)*etmp(nord, id1)*etmp(nord, id2))
                              end do
                           end do
                        end do
                     end do
                     indx_last = indx
                  end if

               end if                  ! mld_polyc small boxes

            else
               uu(1:3) = xpnp(1:3, ia) - xpnp(1:3, ic)
               ds = MatMul(uu, bg)
               WHERE ((ds .GT. 0.5d0) .OR. (ds .LT. -0.5d0))
                  ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
               END WHERE
               uu(:) = MatMul(at(:, :), ds(:))/A2cm
               rr = dsqrt(Sum(uu(1:3)**2))

               if (rr >= r_cut) cycle
               !write (*,*) rangml, rr, r_cut
               ia_n = ia_n + 1

               ik = config_desc(iconf)%kind_neigh(ia, ia_n)

               if ((mld_order == mld_linear) .or. (mld_order == mld_linear_extended)) then
                  force_ij(1:dim_xdesc, 1:3) = sign_stress_big_box*config_desc(iconf)%force(1:dim_xdesc, ia, ia_n, 1:3)
               end if


               if (mld_order == mld_kernel) then
                  force_ij_kernel(1:dim_xdesc, 1:3) = config_desc(iconf)%force(1:dim_xdesc, ia, ia_n, 1:3)
                  force_ij_kernel(dim_xdesc + 1:dim_xdesc_kernel - 1, 1:3) = config_desc(iconf)%force_kernel(1:dim_kernel, ia, ia_n, 1:3)
               end if                  ! mld_kernel

               if (activate_k2b) then
                  force_ij_k2b(1:dim_kernel_2b, 1:3) = config_desc(iconf)%force_k2b(1:dim_kernel_2b, ia, ia_n, 1:3)
               end if                  ! mld_kernel

               if (mld_order == mld_quadratic) then
                  force_ij_quadratic(1:dim_xdesc, 1:3) = sign_stress_big_box*config_desc(iconf)%force(1:dim_xdesc, ia, ia_n, 1:3)
                  indx_last = 0
                  if (mld_type_quadratic == mld_type_quadratic_zaxa) then 
                    do id1 = 1, dim_xdesc
                       do id2 = 1, dim_kernel_2b
                          indx = dim_kernel_2b*(id1 - 1) + id2 - 1 + dim_xdesc + 1
                          force_ij_quadratic(indx, 1:3) = sign_stress_big_box*(config_desc(iconf)%force(id1, ia, ia_n, 1:3)*config_desc(iconf)%energy_k2b(id2, ik) + &
                             config_desc(iconf)%force_k2b(id2, ia, ia_n, 1:3)*config_desc(iconf)%energy(id1, ik))
                       end do
                    end do
                  else 
                    do id1 = 1, dim_xdesc
                       do id2 = 1, dim_xdesc
                          indx = dim_xdesc*(id1 - 1) + id2 - 1 + dim_xdesc + 1
                          force_ij_quadratic(indx, 1:3) = sign_stress_big_box*(config_desc(iconf)%force(id1, ia, ia_n, 1:3)*config_desc(iconf)%energy(id2, ik) + &
                             config_desc(iconf)%force(id2, ia, ia_n, 1:3)*config_desc(iconf)%energy(id1, ik))
                       end do
                    end do
                  end if 
               end if                  ! mld_quadratic


               if (mld_order == mld_polyc) then

                  do nord = 1, polyc_n_hermite
                     do id1 = 1, dim_xdesc
                        call hermite(nord, config_desc(iconf)%energy(id1, ik), y_her, d_y_her)
                        etmp(nord, id1) = y_her
                        d_etmp(nord, id1) = d_y_her
                     end do
                  end do

                  indx_last = 0
                  if (polyc_n_poly >= 1) then
                     indx = indx_last
                     do nord = 1, polyc_n_hermite
                        do id1 = 1, dim_xdesc
                           indx = indx + 1
                           force_ij_polyc(indx, 1:3) = config_desc(iconf)%force(id1, ia, ia_n, 1:3)*d_etmp(nord, id1)
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
                              !force_ij_polyc(indx,1:3) = (config_desc(iconf)%force(id1,ia,ia_n,1:3)*config_desc(iconf)%energy(id2,ik) + &
                              !                            config_desc(iconf)%force(id2,ia,ia_n,1:3)*config_desc(iconf)%energy(id1,ik) )
                              force_ij_polyc(indx, 1:3) = (config_desc(iconf)%force(id1, ia, ia_n, 1:3)*d_etmp(nord, id1)*etmp(nord, id2) + &
                                 config_desc(iconf)%force(id2, ia, ia_n, 1:3)*d_etmp(nord, id2)*etmp(nord, id1))
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
                                 !force_ij_polyc(indx,1:3) = (config_desc(iconf)%force(id1,ia,ia_n,1:3)*config_desc(iconf)%energy(id2,ik)*config_desc(iconf)%energy(id3,ik) + &
                                 !                            config_desc(iconf)%force(id2,ia,ia_n,1:3)*config_desc(iconf)%energy(id1,ik)*config_desc(iconf)%energy(id3,ik) + &
                                 !                            config_desc(iconf)%force(id3,ia,ia_n,1:3)*config_desc(iconf)%energy(id1,ik)*config_desc(iconf)%energy(id2,ik) )
                                 force_ij_polyc(indx, 1:3) = (config_desc(iconf)%force(id1, ia, ia_n, 1:3)*d_etmp(nord, id1)*etmp(nord, id2)*etmp(nord, id3) + &
                                    config_desc(iconf)%force(id2, ia, ia_n, 1:3)*d_etmp(nord, id2)*etmp(nord, id1)*etmp(nord, id3) + &
                                    config_desc(iconf)%force(id3, ia, ia_n, 1:3)*d_etmp(nord, id3)*etmp(nord, id1)*etmp(nord, id2))
                              end do
                           end do
                        end do
                     end do
                     indx_last = indx
                  end if

               end if                  ! mld_polyc big boxes

            end if
            !rr     = config_real(iconf)%r_ij(ia,ic)
            !uu(:)  = config_real(iconf)%u_ij(ia,ic,:)

            if ((mld_order == mld_linear) .or. (mld_order == mld_linear_extended)) then
               st_temp(1:dim_xdesc, 1) = st_temp(1:dim_xdesc, 1) + force_ij(1:dim_xdesc, 1)*uu(1)/rr
               st_temp(1:dim_xdesc, 2) = st_temp(1:dim_xdesc, 2) + force_ij(1:dim_xdesc, 2)*uu(2)/rr
               st_temp(1:dim_xdesc, 3) = st_temp(1:dim_xdesc, 3) + force_ij(1:dim_xdesc, 3)*uu(3)/rr
               st_temp(1:dim_xdesc, 4) = st_temp(1:dim_xdesc, 4) + force_ij(1:dim_xdesc, 2)*uu(3)/rr
               st_temp(1:dim_xdesc, 5) = st_temp(1:dim_xdesc, 5) + force_ij(1:dim_xdesc, 1)*uu(3)/rr
               st_temp(1:dim_xdesc, 6) = st_temp(1:dim_xdesc, 6) + force_ij(1:dim_xdesc, 1)*uu(2)/rr
            end if

            if (mld_order == mld_kernel) then
               st_temp_kernel(1:dim_xdesc_kernel - 1, 1) = st_temp_kernel(1:dim_xdesc_kernel - 1, 1) + force_ij_kernel(1:dim_xdesc_kernel - 1, 1)*uu(1)/rr
               st_temp_kernel(1:dim_xdesc_kernel - 1, 2) = st_temp_kernel(1:dim_xdesc_kernel - 1, 2) + force_ij_kernel(1:dim_xdesc_kernel - 1, 2)*uu(2)/rr
               st_temp_kernel(1:dim_xdesc_kernel - 1, 3) = st_temp_kernel(1:dim_xdesc_kernel - 1, 3) + force_ij_kernel(1:dim_xdesc_kernel - 1, 3)*uu(3)/rr
               st_temp_kernel(1:dim_xdesc_kernel - 1, 4) = st_temp_kernel(1:dim_xdesc_kernel - 1, 4) + force_ij_kernel(1:dim_xdesc_kernel - 1, 2)*uu(3)/rr
               st_temp_kernel(1:dim_xdesc_kernel - 1, 5) = st_temp_kernel(1:dim_xdesc_kernel - 1, 5) + force_ij_kernel(1:dim_xdesc_kernel - 1, 1)*uu(3)/rr
               st_temp_kernel(1:dim_xdesc_kernel - 1, 6) = st_temp_kernel(1:dim_xdesc_kernel - 1, 6) + force_ij_kernel(1:dim_xdesc_kernel - 1, 1)*uu(2)/rr
            end if

            if (activate_k2b) then
               st_temp_k2b(1:dim_kernel_2b, 1) = st_temp_k2b(1:dim_kernel_2b, 1) + force_ij_k2b(1:dim_kernel_2b, 1)*uu(1)/rr
               st_temp_k2b(1:dim_kernel_2b, 2) = st_temp_k2b(1:dim_kernel_2b, 2) + force_ij_k2b(1:dim_kernel_2b, 2)*uu(2)/rr
               st_temp_k2b(1:dim_kernel_2b, 3) = st_temp_k2b(1:dim_kernel_2b, 3) + force_ij_k2b(1:dim_kernel_2b, 3)*uu(3)/rr
               st_temp_k2b(1:dim_kernel_2b, 4) = st_temp_k2b(1:dim_kernel_2b, 4) + force_ij_k2b(1:dim_kernel_2b, 2)*uu(3)/rr
               st_temp_k2b(1:dim_kernel_2b, 5) = st_temp_k2b(1:dim_kernel_2b, 5) + force_ij_k2b(1:dim_kernel_2b, 1)*uu(3)/rr
               st_temp_k2b(1:dim_kernel_2b, 6) = st_temp_k2b(1:dim_kernel_2b, 6) + force_ij_k2b(1:dim_kernel_2b, 1)*uu(2)/rr
            end if


            if (mld_order == mld_quadratic) then
               st_temp_quadratic(1:dim_xdesc_quadratic - 1, 1) = st_temp_quadratic(1:dim_xdesc_quadratic - 1, 1) + force_ij_quadratic(1:dim_xdesc_quadratic - 1, 1)*uu(1)/rr
               st_temp_quadratic(1:dim_xdesc_quadratic - 1, 2) = st_temp_quadratic(1:dim_xdesc_quadratic - 1, 2) + force_ij_quadratic(1:dim_xdesc_quadratic - 1, 2)*uu(2)/rr
               st_temp_quadratic(1:dim_xdesc_quadratic - 1, 3) = st_temp_quadratic(1:dim_xdesc_quadratic - 1, 3) + force_ij_quadratic(1:dim_xdesc_quadratic - 1, 3)*uu(3)/rr
               st_temp_quadratic(1:dim_xdesc_quadratic - 1, 4) = st_temp_quadratic(1:dim_xdesc_quadratic - 1, 4) + force_ij_quadratic(1:dim_xdesc_quadratic - 1, 2)*uu(3)/rr
               st_temp_quadratic(1:dim_xdesc_quadratic - 1, 5) = st_temp_quadratic(1:dim_xdesc_quadratic - 1, 5) + force_ij_quadratic(1:dim_xdesc_quadratic - 1, 1)*uu(3)/rr
               st_temp_quadratic(1:dim_xdesc_quadratic - 1, 6) = st_temp_quadratic(1:dim_xdesc_quadratic - 1, 6) + force_ij_quadratic(1:dim_xdesc_quadratic - 1, 1)*uu(2)/rr
            end if



            if (mld_order == mld_polyc) then
               st_temp_polyc(1:dim_xdesc_polyc - 1, 1) = st_temp_polyc(1:dim_xdesc_polyc - 1, 1) + force_ij_polyc(1:dim_xdesc_polyc - 1, 1)*uu(1)/rr
               st_temp_polyc(1:dim_xdesc_polyc - 1, 2) = st_temp_polyc(1:dim_xdesc_polyc - 1, 2) + force_ij_polyc(1:dim_xdesc_polyc - 1, 2)*uu(2)/rr
               st_temp_polyc(1:dim_xdesc_polyc - 1, 3) = st_temp_polyc(1:dim_xdesc_polyc - 1, 3) + force_ij_polyc(1:dim_xdesc_polyc - 1, 3)*uu(3)/rr
               st_temp_polyc(1:dim_xdesc_polyc - 1, 4) = st_temp_polyc(1:dim_xdesc_polyc - 1, 4) + force_ij_polyc(1:dim_xdesc_polyc - 1, 2)*uu(3)/rr
               st_temp_polyc(1:dim_xdesc_polyc - 1, 5) = st_temp_polyc(1:dim_xdesc_polyc - 1, 5) + force_ij_polyc(1:dim_xdesc_polyc - 1, 1)*uu(3)/rr
               st_temp_polyc(1:dim_xdesc_polyc - 1, 6) = st_temp_polyc(1:dim_xdesc_polyc - 1, 6) + force_ij_polyc(1:dim_xdesc_polyc - 1, 1)*uu(2)/rr
            end if

         end do
      end do

      if ((mld_order == mld_linear) .or. (mld_order == mld_linear_extended)) then
         ! write (*,*) '!!!!!!!!!!!!!!!!', iconf, dim_xdesc, size(config_desc), size( config_desc(1)%pack_stress,1), size( config_desc(1)%pack_stress,2)
         config_desc(iconf)%pack_stress_linear(1:dim_xdesc, 1:6) = sign_stress*st_temp(1:dim_xdesc, 1:6)*3.d0/(config_real(iconf)%volume)
      end if


      if (mld_order == mld_kernel) then
         ! write (*,*) '!!!!!!!!!!!!!!!!', iconf, dim_xdesc_kernel, size(config_desc), size( config_desc(1)%pack_stress,1), size( config_desc(1)%pack_stress,2)
         config_desc(iconf)%pack_stress_linear(1:dim_xdesc, 1:6) = sign_stress*st_temp_kernel(1:dim_xdesc, 1:6)*3.d0/(config_real(iconf)%volume)
         config_desc(iconf)%pack_stress_kernel(1:dim_xdesc_kernel - 1, 1:6) = sign_stress*st_temp_kernel(1:dim_xdesc_kernel - 1, 1:6)*3.d0/(config_real(iconf)%volume)
      end if

      if (activate_k2b) then
         config_desc(iconf)%pack_stress_k2b(1:dim_kernel_2b, 1:6) = sign_stress*st_temp_k2b(1:dim_kernel_2b, 1:6)*3.d0/(config_real(iconf)%volume)
      end if


      if (mld_order == mld_quadratic) then
         ! write (*,*) '!!!!!!!!!!!!!!!!', iconf, dim_xdesc_quadratic, size(config_desc), size( config_desc(1)%pack_stress,1), size( config_desc(1)%pack_stress,2)
         config_desc(iconf)%pack_stress_linear(1:dim_xdesc, 1:6) = sign_stress*st_temp_quadratic(1:dim_xdesc, 1:6)*3.d0/(config_real(iconf)%volume)
         config_desc(iconf)%pack_stress_quadratic(1:dim_xdesc_quadratic - 1, 1:6) = sign_stress*st_temp_quadratic(1:dim_xdesc_quadratic - 1, 1:6)*3.d0/(config_real(iconf)%volume)
      end if


      if (mld_order == mld_polyc) then
         !write (*,*) '!!!!!!!!!!!!!!!!', iconf, dim_xdesc_polyc, size(config_desc), size( config_desc(1)%pack_stress_polyc,1), size( config_desc(1)%pack_stress_polyc,2)
         config_desc(iconf)%pack_stress_linear(1:dim_xdesc, 1:6) = sign_stress*st_temp_polyc(1:dim_xdesc, 1:6)*3.d0/(config_real(iconf)%volume)
         config_desc(iconf)%pack_stress_polyc(1:dim_xdesc_polyc - 1, 1:6) = sign_stress*st_temp_polyc(1:dim_xdesc_polyc - 1, 1:6)*3.d0/(config_real(iconf)%volume)
      end if

      _MLD_END_

   end subroutine para_pack_stress_descriptor

end module


! OLD STUFFS


! subroutine serial_snap_pack_stress_descriptor(iconf)

! This subroutine has strong interaction with neighbours subroutines
! The one proposed by Milady and the one proposed by NDM depending if the cell box is small or large.
! small=.true. ->  MiLaDy
! small=.false. -> NDM
!
!use ondm_gen_com_m, only : imm, lperiod, indi2, A2cm, bg, at
!use ondm_tab_imm_m, only : xp, iwmax2
 !!MiLaDy_interaction
!use ml_in_ndm_module, only: r_cut, sign_stress, sign_stress_big_box
!use derived_types, only: config_desc, config_real
!use temporary_data_cov, only: dim_xdesc
!implicit none
!integer, intent(in) :: iconf
!integer :: ia, ic
!real(kind(0.d0)) :: st_temp(1:dim_xdesc,6), rr, uu(3), ds(3), force_ij(dim_xdesc,3)
!integer :: iw,iw1,iw2, ia_n
!
!real(kind(0.d0)), dimension(:,:), allocatable :: xpnp
!logical :: small
!if (.not.(config_real(iconf)%has_stress)) return
!
!if (allocated(config_desc(iconf)%pack_stress_linear))     deallocate(config_desc(iconf)%pack_stress )  ; allocate(config_desc(iconf)%pack_stress_linear(dim_xdesc,6))
!small= config_real(iconf)%small
!allocate(xpnp(3,imm))
!if (lperiod) then
!  xpnp(:,:)=xp(:,:)
!else
!call ondm_notperiod(xp,xpnp)
!end if
 !!call cryst_to_cart (imm, xpnp, bg, -1)
!
! iw2=0
! st_temp(:,:)=0.d0
! do ia=1, config_real(iconf)%nat ! descriptor index
!     !begin small box or not 1/
!     if (small) then
!      iw1=1
!      iw2=config_real(iconf)%n_neigh(ia)
!     else
!      iw1=iw2+1
!      iw2=iwmax2(ia)
!      !write (*,*) 'debug neigh', iw1, iw2
!     end if
!     ia_n=0
!     do iw = iw1, iw2
!        if (small) then
!           ic = config_real(iconf)%kind_neigh(ia,iw)
!        else
!           ic=indi2(iw)
!           if (ic==ia) cycle
!        end if
!
!        if (small) then
!          rr = config_real(iconf)%r_ij(ia,iw)
!          uu (:) = config_real(iconf)%u_ij(ia,iw,:)
!          force_ij (1:dim_xdesc,1:3) = config_desc(iconf)%force(1:dim_xdesc,ia,iw,1:3)
!        else
!          uu(1:3) = xpnp(1:3,ia) - xpnp(1:3,ic)
!          ds=MatMul(uu,bg)
!          WHERE ( (ds.GT.0.5d0).OR.(ds.LT.-0.5d0) )
!            ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
!          END WHERE
!          uu(:) = MatMul(at(:,:),ds(:))/A2cm
!          rr = dsqrt( Sum( uu(1:3)**2 ) )
!          if (rr >= r_cut) cycle
!          ia_n=ia_n+1
!          force_ij(1:dim_xdesc,1:3) = sign_stress_big_box*config_desc(iconf)%force(1:dim_xdesc,ia,ia_n,1:3)
!        end if
!                !rr     = config_real(iconf)%r_ij(ia,ic)
!                !uu(:)  = config_real(iconf)%u_ij(ia,ic,:)
!                st_temp(1:dim_xdesc, 1) = st_temp(1:dim_xdesc, 1) + force_ij(1:dim_xdesc,1)*uu(1)/rr
!                st_temp(1:dim_xdesc, 2) = st_temp(1:dim_xdesc, 2) + force_ij(1:dim_xdesc,2)*uu(2)/rr
!                st_temp(1:dim_xdesc, 3) = st_temp(1:dim_xdesc, 3) + force_ij(1:dim_xdesc,3)*uu(3)/rr
!                st_temp(1:dim_xdesc, 4) = st_temp(1:dim_xdesc, 4) + force_ij(1:dim_xdesc,2)*uu(3)/rr
!                st_temp(1:dim_xdesc, 5) = st_temp(1:dim_xdesc, 5) + force_ij(1:dim_xdesc,1)*uu(3)/rr
!                st_temp(1:dim_xdesc, 6) = st_temp(1:dim_xdesc, 6) + force_ij(1:dim_xdesc,1)*uu(2)/rr
!     end do
!  end do
!
!config_desc(iconf)%pack_stress_linear(1:dim_xdesc,1:6)=sign_stress*st_temp(1:dim_xdesc,1:6)*3.d0/(config_real(iconf)%volume)
!
!end subroutine serial_snap_pack_stress_descriptor
!
!
!subroutine old_para_snap_pack_stress_descriptor(iconf)

! This subroutine has strong interaction with neighbours subroutines
! The one proposed by Milady and the one proposed by NDM depending if the cell box is small or large.
! small=.true. ->  MiLaDy
! small=.false. -> NDM
!
!use ondm_gen_com_m, only : imm, lperiod, indi2, A2cm, bg, at
!use ondm_tab_imm_m, only : xp, iwmax2
 !!MiLaDy_interaction
!use ml_in_ndm_module, only: r_cut, sign_stress, sign_stress_big_box, i_start_at, i_final_at
!use derived_types, only: config_desc, config_real
!use temporary_data_cov, only: dim_xdesc
!implicit none
!integer, intent(in) :: iconf
!integer :: ia, ic
!real(kind(0.d0)) :: st_temp(1:dim_xdesc,6), rr, uu(3), ds(3), force_ij(dim_xdesc,3)
!integer :: iw,iw1,iw2, ia_n
!
!real(kind(0.d0)), dimension(:,:), allocatable :: xpnp
!logical :: small
!
!if ((i_start_at==0).and.(i_final_at==0)) return
!if (.not.(config_real(iconf)%has_stress)) return
!
!small= config_real(iconf)%small
!allocate(xpnp(3,imm))
!if (lperiod) then
!  xpnp(:,:)=xp(:,:)
!else
!call ondm_notperiod(xp,xpnp)
!end if
 !!call cryst_to_cart (imm, xpnp, bg, -1)
!
!
! iw2=0
! st_temp(:,:)=0.d0
! !debug write (*,*) i_start_at, i_final_at, rangml
! if (i_start_at==1)  iw2=0
! if (i_start_at > 1) iw2=iwmax2(i_start_at-1)
! do ia=i_start_at, i_final_at ! descriptor index!
!     !begin small box or not 1/
!     if (small) then
!      iw1=1
!      iw2=config_real(iconf)%n_neigh(ia)
!     else
!      iw1=iw2+1
!      iw2=iwmax2(ia)
!      !write (*,*) 'debug neigh', iw1, iw2
!     end if
!     ia_n=0
!     do iw = iw1, iw2
!        if (small) then
!           ic = config_real(iconf)%kind_neigh(ia,iw)
!        else
!           ic=indi2(iw)
!           if (ic==ia) cycle
!        end if
!
!        if (small) then
!          rr = config_real(iconf)%r_ij(ia,iw)
!          uu (:) = config_real(iconf)%u_ij(ia,iw,:)
!          force_ij (1:dim_xdesc,1:3) = config_desc(iconf)%force(1:dim_xdesc,ia,iw,1:3)
!        else
!          uu(1:3) = xpnp(1:3,ia) - xpnp(1:3,ic)
!          ds=MatMul(uu,bg)
!          WHERE ( (ds.GT.0.5d0).OR.(ds.LT.-0.5d0) )
!            ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
!          END WHERE
!          uu(:) = MatMul(at(:,:),ds(:))/A2cm
!          rr = dsqrt( Sum( uu(1:3)**2 ) )
!
!          if (rr >= r_cut) cycle
!          !write (*,*) rangml, rr, r_cut
!          ia_n=ia_n+1
!          force_ij(1:dim_xdesc,1:3) = sign_stress_big_box*config_desc(iconf)%force(1:dim_xdesc,ia,ia_n,1:3)
!        end if
!                !rr     = config_real(iconf)%r_ij(ia,ic)
!                !uu(:)  = config_real(iconf)%u_ij(ia,ic,:)
!                st_temp(1:dim_xdesc, 1) = st_temp(1:dim_xdesc, 1) + force_ij(1:dim_xdesc,1)*uu(1)/rr
!                st_temp(1:dim_xdesc, 2) = st_temp(1:dim_xdesc, 2) + force_ij(1:dim_xdesc,2)*uu(2)/rr
!                st_temp(1:dim_xdesc, 3) = st_temp(1:dim_xdesc, 3) + force_ij(1:dim_xdesc,3)*uu(3)/rr
!                st_temp(1:dim_xdesc, 4) = st_temp(1:dim_xdesc, 4) + force_ij(1:dim_xdesc,2)*uu(3)/rr
!                st_temp(1:dim_xdesc, 5) = st_temp(1:dim_xdesc, 5) + force_ij(1:dim_xdesc,1)*uu(3)/rr
!                st_temp(1:dim_xdesc, 6) = st_temp(1:dim_xdesc, 6) + force_ij(1:dim_xdesc,1)*uu(2)/rr
!     end do
!  end do
!
!config_desc(iconf)%pack_stress_linear(1:dim_xdesc,1:6)=sign_stress*st_temp(1:dim_xdesc,1:6)*3.d0/(config_real(iconf)%volume)
!
!end subroutine old_para_snap_pack_stress_descriptor
