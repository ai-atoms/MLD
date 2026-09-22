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

#include "../../MLD_MACROS.INC"



module kernel_tools
   use mld_logger
   use module_kind_variables, ONLY: kind_double
   implicit none

contains

   subroutine inner_kernel_random_se (xx, zz, phase_ik, &
      vfk_trans1, vfk_trans2, vfk_trans3, &
      desc_forces_local, fact_random, &
      ktemp_norm, dktemp_imm, vtemp0)
      real(kind_double), dimension(:), intent(in) :: xx, zz
      real(kind_double),  intent(in) :: phase_ik
      real(kind_double), dimension(:,:) :: vfk_trans1, vfk_trans2, vfk_trans3
      logical, intent(in) :: desc_forces_local
      real(kind_double), intent(in) :: fact_random

      real(kind_double), intent(out) :: ktemp_norm
      real(kind_double), dimension(:,:), intent(out) :: dktemp_imm
      real(kind_double), dimension(:), intent(out) ::  vtemp0
      real(kind_double), dimension(:), allocatable :: vtemp1, vtemp2, vtemp3
      real(kind_double), dimension(3) :: dktemp_norm
      integer :: ja_neigh, dim_xdesc, inn
      real(kind_double) :: dtmp_jaik
      _NAMECURRENT_("inner_kernel_random_se")
      _MLD_BEGIN_

      ja_neigh = size(vfk_trans1, dim=1)
      dim_xdesc= size(vfk_trans1, dim=2)

      dtmp_jaik = dot_product(xx, zz) + phase_ik
      ktemp_norm = fact_random * cos(dtmp_jaik)

      if (desc_forces_local) then
         allocate(vtemp1(ja_neigh), vtemp2(ja_neigh), vtemp3(ja_neigh))
         vtemp0(:) = 0.d0
         call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans1, ja_neigh, zz, 1, 0.d0, vtemp1, 1)
         call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans2, ja_neigh, zz, 1, 0.d0, vtemp2, 1)
         call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans3, ja_neigh, zz, 1, 0.d0, vtemp3, 1)
         do inn = 1, ja_neigh
            dktemp_norm(1) = -vtemp1(inn)*sin(dtmp_jaik) * fact_random
            dktemp_norm(2) = -vtemp2(inn)*sin(dtmp_jaik) * fact_random
            dktemp_norm(3) = -vtemp3(inn)*sin(dtmp_jaik) * fact_random
            vtemp0(:) = vtemp0(:) - dktemp_norm(1:3)
            !config_desc(iconf)%force_kernel(ik, ja, inn, 1:3) = dktemp_norm(1:3)
            dktemp_imm(1:3, inn) = dktemp_norm(1:3)
         end do

         !config_desc(iconf)%force_kernel(ik, ja, 0, 1:3) = vtemp0(1:3)
      end if


      _MLD_END_
   end subroutine inner_kernel_random_se

   subroutine  inner_kernel_po( xx, zz, vfk_trans1, vfk_trans2, vfk_trans3, &
      vtemp_ja1,  vtemp_ja2,  vtemp_ja3, &
      desc_forces_local, sqrt_ker_xx, dker_xx, ker_xx, sigma_kse2, l2kse, kernel_power, &
      ktemp_norm, dktemp_inn,  vtemp0)
      ! This is a normalized polynomial kernel of the form: 
      ! f(x,y) =  (sigma**2  + <x,y>/(2 l^2))^p  
      ! k(x,y) = f(x,y) / sqrt(f(x,x) f(y,y))
      implicit none
      real(kind_double), dimension(:,:), intent(in) :: vfk_trans1, vfk_trans2, vfk_trans3
      real(kind_double), dimension(:), intent(in) :: xx, zz, vtemp_ja1, vtemp_ja2, vtemp_ja3
      real(kind_double), intent(out)  :: ktemp_norm
      real(kind_double), dimension(:), intent(out) ::  vtemp0
      real(kind_double), dimension(:,:), intent(out) :: dktemp_inn
      real(kind_double), intent(in)  :: sqrt_ker_xx, dker_xx, ker_xx, sigma_kse2, l2kse
      logical, intent(in) :: desc_forces_local
      real(kind_double), intent(in) :: kernel_power
      real(kind_double) :: dot_xz, dtmp, ker_zz, sqrt_ker_zz, ktemp,  tmp_f1, aaa, bbb
      integer :: ja_neigh, dim_xdesc, inn
      real(kind_double), dimension(:), allocatable  :: vtemp1, vtemp2, vtemp3
      real(kind_double), dimension(3) :: dktemp, dktemp_ja, dktemp_norm
      real(kind_double), external :: ddot

      _NAMECURRENT_("inner_kernel_po")
      _MLD_BEGIN_
      ja_neigh = size(vfk_trans1, dim=1)
      dim_xdesc= size(vfk_trans1, dim=2)
      dot_xz = dot_product(xx(:), zz(:))
      !dot_xz = ddot(dim_xdesc, xx(:), 1, zz(:), 1)
      dtmp = sigma_kse2 + dot_xz/l2kse
      ktemp = dtmp**kernel_power

      ker_zz = (sigma_kse2 + dot_product(zz, zz)/l2kse)**kernel_power
      sqrt_ker_zz = dsqrt(ker_zz)

      aaa = (sqrt_ker_zz*sqrt_ker_xx)
      bbb = 2.d0*ker_xx
      ktemp_norm = ktemp/aaa
      if (desc_forces_local) then
         allocate(vtemp1(ja_neigh), vtemp2(ja_neigh), vtemp3(ja_neigh))
         vtemp0(:) = 0.d0
         tmp_f1 = dble(kernel_power)*dtmp**(kernel_power - 1)/l2kse
         call dgemv('N', ja_neigh, dim_xdesc, tmp_f1, vfk_trans1, ja_neigh, zz, 1, 0.d0, vtemp1, 1)
         call dgemv('N', ja_neigh, dim_xdesc, tmp_f1, vfk_trans2, ja_neigh, zz, 1, 0.d0, vtemp2, 1)
         call dgemv('N', ja_neigh, dim_xdesc, tmp_f1, vfk_trans3, ja_neigh, zz, 1, 0.d0, vtemp3, 1)
         do inn = 1, ja_neigh
            dktemp(1) = vtemp1(inn) !*tmp_f1
            dktemp(2) = vtemp2(inn) !*tmp_f1
            dktemp(3) = vtemp3(inn) !*tmp_f1
            !dktemp = \nabla_b ktilde (x,x) = \partial(ktilde(x,x))/\partial(x)  \nabla_b x
            dktemp_ja(1) = vtemp_ja1(inn)*dker_xx
            dktemp_ja(2) = vtemp_ja2(inn)*dker_xx
            dktemp_ja(3) = vtemp_ja3(inn)*dker_xx
            dktemp_norm(1:3) = (dktemp(1:3) - ktemp*dktemp_ja(1:3)/bbb)/aaa
            dktemp_inn(1:3, inn) = dktemp_norm(1:3)
            vtemp0(:) = vtemp0(:) - dktemp_norm(1:3)
         end do
      end if                  ! desc_forces_local

      _MLD_END_
   end subroutine inner_kernel_po


   subroutine  inner_kernel_po_scaled( xx, zz, vfk_trans1, vfk_trans2, vfk_trans3, &
      vtemp_ja1,  vtemp_ja2,  vtemp_ja3, &
      desc_forces_local, dot_xx, sigma_kse2,  kernel_power, &
      ktemp_norm, dktemp_inn,  vtemp0)
      ! This is scaled version of the normalize polynomial kernel of the form 
      ! f(x,z) = (<x,z>/(|x| |z|))^p
      ! k(x,z) = sigma^2  f(x,z)
      implicit none
      real(kind_double), dimension(:,:), intent(in) :: vfk_trans1, vfk_trans2, vfk_trans3
      real(kind_double), dimension(:), intent(in) :: xx, zz, vtemp_ja1, vtemp_ja2, vtemp_ja3
      real(kind_double), intent(out)  :: ktemp_norm
      real(kind_double), dimension(:), intent(out) ::  vtemp0
      real(kind_double), dimension(:,:), intent(out) :: dktemp_inn
      real(kind_double), intent(in)  :: sigma_kse2, dot_xx 
      logical, intent(in) :: desc_forces_local
      real(kind_double), intent(in) :: kernel_power
      real(kind_double) :: dot_xz, dtmp, ktemp,  tmp_f1, norm_xx, norm_zz, norm_inv, norm_invb
      integer :: ja_neigh, dim_xdesc, inn
      real(kind_double), dimension(:), allocatable  :: vtemp1, vtemp2, vtemp3
      real(kind_double), dimension(3) :: dktemp, dktemp_ja, dktemp_norm

      _NAMECURRENT_("inner_kernel_po_scaled")
      _MLD_BEGIN_
      ja_neigh = size(vfk_trans1, dim=1)
      dim_xdesc= size(vfk_trans1, dim=2)
      norm_xx = sqrt(dot_xx) !sqrt(dot_product(xx(:), xx(:)))
      norm_zz = sqrt(dot_product(zz(:), zz(:)))
      dot_xz = dot_product(xx(:), zz(:))
      !dot_xz = ddot(dim_xdesc, xx(:), 1, zz(:), 1)
      norm_inv = 1.d0/(norm_xx*norm_zz)
      norm_invb = 1.d0/(norm_xx**3*norm_zz)
      dtmp = dot_xz*norm_inv 
      ktemp = dtmp**kernel_power
      ! kernel(x,z) is : 
      ktemp_norm = sigma_kse2 * ktemp
      if (desc_forces_local) then
         allocate(vtemp1(ja_neigh), vtemp2(ja_neigh), vtemp3(ja_neigh))
         vtemp0(:) = 0.d0
         tmp_f1 = sigma_kse2 * dble(kernel_power)*dtmp**(kernel_power - 1)*norm_inv 
         call dgemv('N', ja_neigh, dim_xdesc, tmp_f1, vfk_trans1, ja_neigh, zz, 1, 0.d0, vtemp1, 1)
         call dgemv('N', ja_neigh, dim_xdesc, tmp_f1, vfk_trans2, ja_neigh, zz, 1, 0.d0, vtemp2, 1)
         call dgemv('N', ja_neigh, dim_xdesc, tmp_f1, vfk_trans3, ja_neigh, zz, 1, 0.d0, vtemp3, 1)
         do inn = 1, ja_neigh
            dktemp(1) = vtemp1(inn) !*tmp_f1
            dktemp(2) = vtemp2(inn) !*tmp_f1
            dktemp(3) = vtemp3(inn) !*tmp_f1
            !dktemp = \nabla_b ktilde (x,x) = \partial(ktilde(x,x))/\partial(x)  \nabla_b x
            dktemp_ja(1) = vtemp_ja1(inn)*dot_xz * norm_invb
            dktemp_ja(2) = vtemp_ja2(inn)*dot_xz * norm_invb
            dktemp_ja(3) = vtemp_ja3(inn)*dot_xz * norm_invb
            dktemp_norm(1:3) = dktemp(1:3) - dktemp_ja(1:3)
            dktemp_inn(1:3, inn) = dktemp_norm(1:3)
            vtemp0(:) = vtemp0(:) - dktemp_norm(1:3)
         end do
      end if                  ! desc_forces_local

      _MLD_END_
   end subroutine inner_kernel_po_scaled





   subroutine inner_kernel_se (xx, vfk_trans1, vfk_trans2, vfk_trans3, &
      desc_forces_local, soverl, sigma_kse2, l2kse, &
      ktemp_norm, dktemp_inn,  vtemp0)

      implicit none
      real(kind_double), dimension(:), intent(in) :: xx
      real(kind_double), dimension(:,:), intent(in) :: vfk_trans1, vfk_trans2, vfk_trans3
      logical, intent(in) :: desc_forces_local
      real(kind_double), intent(in) :: soverl, sigma_kse2, l2kse
      real(kind_double), intent(out)  :: ktemp_norm
      real(kind_double), dimension(:), intent(out) ::  vtemp0
      real(kind_double), dimension(:,:), intent(out) :: dktemp_inn
      real(kind_double), dimension(:), allocatable :: vtemp1, vtemp2, vtemp3
      real(kind_double), dimension(3) :: dktemp, dktemp_norm
      integer :: ja_neigh, dim_xdesc, inn
      real(kind_double) :: dtmp, ktemp, tmp_f1

      _NAMECURRENT_("inner_kernel_se")
      _MLD_BEGIN_
      ja_neigh =  size(vfk_trans1, dim=1)
      dim_xdesc= size(vfk_trans1, dim=2)

      ! etemp(:) = energy_ja(:) - draft_kernel(:, ik)
      dtmp = dot_product(xx(:), xx(:))
      ktemp = dexp(-dtmp/l2kse)
      ktemp_norm = sigma_kse2*ktemp
      if (desc_forces_local) then
         allocate(vtemp1(ja_neigh), vtemp2(ja_neigh), vtemp3(ja_neigh))
         vtemp0(:) = 0.d0
         tmp_f1 = -soverl*ktemp
         call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans1, ja_neigh, xx, 1, 0.d0, vtemp1, 1)
         call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans2, ja_neigh, xx, 1, 0.d0, vtemp2, 1)
         call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans3, ja_neigh, xx, 1, 0.d0, vtemp3, 1)
         do inn = 1, ja_neigh
            dktemp(1) = vtemp1(inn)*tmp_f1
            dktemp(2) = vtemp2(inn)*tmp_f1
            dktemp(3) = vtemp3(inn)*tmp_f1
            dktemp_norm(1:3) = dktemp(1:3)
            !config_desc(iconf)%force_kernel(ik, ja, inn, 1:3) = dktemp_norm(1:3)
            dktemp_inn(1:3, inn) = dktemp_norm(1:3)
            vtemp0(:) = vtemp0(:) - dktemp_norm(1:3)
         end do
         !config_desc(iconf)%force_kernel(ik, ja, 0, 1:3) = vtemp0(1:3)
      end if

      _MLD_END_
   end subroutine inner_kernel_se

   subroutine inner_random_po  (xx, dim_kernel, omegaDF_ik, &
      vfk_trans1, vfk_trans2, vfk_trans3, &
      desc_forces_local, ktemp_norm, dktemp_inn, vtemp0)
      real(kind_double), dimension(:), intent(in) :: xx
      integer, intent(in) :: dim_kernel
      real(kind_double), dimension(:,:), intent(in) :: omegaDF_ik, vfk_trans1, vfk_trans2, vfk_trans3
      logical, intent(in) :: desc_forces_local
      real(kind_double), intent(out) :: ktemp_norm
      real(kind_double), dimension(:,:), intent(out) :: dktemp_inn
      real(kind_double), dimension(:), intent(out) :: vtemp0

      real(kind_double), dimension(:), allocatable :: vtemp1, vtemp2, vtemp3, vtemp4
      real(kind_double), dimension(3) :: dktemp_norm
      real(kind_double), dimension(:), allocatable :: ene_tmp, zz
      real(kind_double), dimension(:, :), allocatable :: force_tmp1, force_tmp2, force_tmp3
      integer :: ja_neigh, dim_xdesc, inn, dim_omega, ii, iid
      real(kind_double) :: dtmp, o_norm
      _NAMECURRENT_("inner_random_po")
      _MLD_BEGIN_

      dim_omega = size(omegaDF_ik, dim=2)
      ja_neigh = size(vfk_trans1, dim=1)
      dim_xdesc= size(xx)
      allocate(zz(dim_xdesc))
      if (allocated(ene_tmp)) deallocate (ene_tmp); allocate (ene_tmp(dim_omega))

      call dgemv('T', dim_xdesc, dim_omega, 1.d0, omegaDF_ik, dim_xdesc, xx, 1, 0.d0, ene_tmp, 1)
      dtmp = 1.d0
      do ii = 1, dim_omega
         !!! ene_tmp(ii) = dot_product(xx(:), omegaDF_ik(:, ii))
         dtmp = dtmp*ene_tmp(ii)
      end do
      o_norm= 1.d0/sqrt(dble(dim_kernel))
      ktemp_norm = dtmp*o_norm
      if (desc_forces_local) then
         allocate(vtemp1(ja_neigh), vtemp2(ja_neigh), vtemp3(ja_neigh), vtemp4(ja_neigh))
         vtemp0(1:3) = 0.d0
         vtemp1(1:ja_neigh) = 0.d0
         vtemp2(1:ja_neigh) = 0.d0
         vtemp3(1:ja_neigh) = 0.d0
         if (allocated(force_tmp1)) deallocate (force_tmp1); allocate (force_tmp1(ja_neigh, dim_omega))
         if (allocated(force_tmp2)) deallocate (force_tmp2); allocate (force_tmp2(ja_neigh, dim_omega))
         if (allocated(force_tmp3)) deallocate (force_tmp3); allocate (force_tmp3(ja_neigh, dim_omega))

         ! force_temp1 = vfktrans1 x omegaDF_ik
         ! C = A x B
         ! ja_neigh x dim_omega = (ja_neigh x dim_xdesc) x (dim_xdesc x dim_omega)
         call dgemm('N', 'N', ja_neigh, dim_omega, dim_xdesc, 1.d0, vfk_trans1, ja_neigh, omegaDF_ik, dim_xdesc, 0.d0, force_tmp1, ja_neigh)
         call dgemm('N', 'N', ja_neigh, dim_omega, dim_xdesc, 1.d0, vfk_trans2, ja_neigh, omegaDF_ik, dim_xdesc, 0.d0, force_tmp2, ja_neigh)
         call dgemm('N', 'N', ja_neigh, dim_omega, dim_xdesc, 1.d0, vfk_trans3, ja_neigh, omegaDF_ik, dim_xdesc, 0.d0, force_tmp3, ja_neigh)

         ! do ii = 1, dim_omega
         !   zz(:) = omegaDF_ik(:, ii)
         !   call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans1, ja_neigh, zz, 1, 0.d0, vtemp4, 1)
         !   force_tmp1(1:ja_neigh, ii) = vtemp4(1:ja_neigh)
         !   call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans2, ja_neigh, zz, 1, 0.d0, vtemp4, 1)
         !   force_tmp2(1:ja_neigh, ii) = vtemp4(1:ja_neigh)
         !   call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans3, ja_neigh, zz, 1, 0.d0, vtemp4, 1)
         !   force_tmp3(1:ja_neigh, ii) = vtemp4(1:ja_neigh)
         ! end do
         do iid = 1, dim_omega
            dtmp = 1.d0
            do ii = 1, dim_omega
               if (ii == iid) cycle
               dtmp = dtmp*ene_tmp(ii)
            end do
            vtemp1(1:ja_neigh) = vtemp1(1:ja_neigh) + dtmp*force_tmp1(1:ja_neigh, iid)
            vtemp2(1:ja_neigh) = vtemp2(1:ja_neigh) + dtmp*force_tmp2(1:ja_neigh, iid)
            vtemp3(1:ja_neigh) = vtemp3(1:ja_neigh) + dtmp*force_tmp3(1:ja_neigh, iid)
         end do
         do inn = 1, ja_neigh
            dktemp_norm(1) = vtemp1(inn)
            dktemp_norm(2) = vtemp2(inn)
            dktemp_norm(3) = vtemp3(inn)
            vtemp0(:) = vtemp0(:) - dktemp_norm(1:3)
            ! config_desc(iconf)%force_kernel(ik, ja, inn, 1:3) = dktemp_norm(1:3)/sqrt(dble(dim_kernel))
            dktemp_inn(1:3, inn) = dktemp_norm(1:3)*o_norm
         end do
         !config_desc(iconf)%force_kernel(ik, ja, 0, 1:3) = vtemp0(1:3)/sqrt(dble(dim_kernel))
         vtemp0 = vtemp0(1:3)*o_norm
      end if

      _MLD_END_
   end subroutine inner_random_po

   subroutine inner_kernel_maha  (xx, mat_sigma_inv,  &
      vfk_trans1, vfk_trans2, vfk_trans3, kernel_power, desc_forces_local, &
      ktemp_norm, dktemp_inn, vtemp0)
      real(kind_double), dimension(:), intent(in) :: xx
      real(kind_double), dimension(:,:), intent(in) :: mat_sigma_inv, vfk_trans1, vfk_trans2, vfk_trans3
      real(kind_double), intent(in) :: kernel_power
      logical, intent(in) :: desc_forces_local
      real(kind_double), intent(out) :: ktemp_norm
      real(kind_double), dimension(:,:), intent(out) :: dktemp_inn
      real(kind_double), dimension(:), intent(out) :: vtemp0

      real(kind_double), dimension(:), allocatable :: vtemp1, vtemp2, vtemp3, energy_ja_maha
      real(kind_double) :: tmp_f1, dtmp_ja_ik, tmp_f2, dktemp(3), dktemp_norm(3)
      integer :: ja_neigh, dim_xdesc, inn
      _NAMECURRENT_("inner_kernel_maha")
      _MLD_BEGIN_

      ja_neigh = size(vfk_trans1, dim=1)
      dim_xdesc= size(xx)
      allocate(energy_ja_maha(dim_xdesc))
      !energy_ja_maha(:) = matmul(mat_sigma_inv(:, :), xx(:))
      !equivalent multipliction of mat_sigma_inv and xx with dgemv
      call dgemv('N', dim_xdesc, dim_xdesc, 1.d0, mat_sigma_inv, dim_xdesc, xx, 1, 0.d0, energy_ja_maha, 1)
      tmp_f1 = dot_product(xx(:), energy_ja_maha(:))
      dtmp_ja_ik = dexp(-tmp_f1**kernel_power/0.05d0)
      ktemp_norm = dtmp_ja_ik
      if (desc_forces_local) then
         allocate(vtemp1(ja_neigh), vtemp2(ja_neigh), vtemp3(ja_neigh))
         vtemp0(:) = 0.d0
         call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans1, ja_neigh, energy_ja_maha, 1, 0.d0, vtemp1, 1)
         call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans2, ja_neigh, energy_ja_maha, 1, 0.d0, vtemp2, 1)
         call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans3, ja_neigh, energy_ja_maha, 1, 0.d0, vtemp3, 1)
         !tmp_f2 = 2.d0*dble(kernel_power)*tmp_f1**(kernel_power - 1)
         tmp_f2 = -2.d0*20.d0*dble(kernel_power)*tmp_f1**(kernel_power - 1)*dtmp_ja_ik
         do inn = 1, ja_neigh
            dktemp(1) = tmp_f2*vtemp1(inn)
            dktemp(2) = tmp_f2*vtemp2(inn)
            dktemp(3) = tmp_f2*vtemp3(inn)
            dktemp_norm(1:3) = dktemp(1:3)                   ! /(sqrt_ker_xx*k_ik2) - ktemp_norm * dktemp_ja(1:3) / (2.d0*k_ja32)
            !config_desc(iconf)%force_kernel(ik, ja, inn, 1:3) = dktemp_norm(1:3)
            dktemp_inn(1:3,inn) = dktemp_norm(1:3)
            vtemp0(:) = vtemp0(:) - dktemp_norm(1:3)
         end do
         !config_desc(iconf)%force_kernel(ik, ja, 0, 1:3) = vtemp0(1:3)
      end if


      _MLD_END_
   end subroutine inner_kernel_maha

   subroutine inner_kernel_random_maha (xx, zz, phase_ik, &
      vfk_trans1, vfk_trans2, vfk_trans3, &
      desc_forces_local, fact_random, &
      ktemp_norm, dktemp_imm, vtemp0)
      real(kind_double), dimension(:), intent(in) :: xx, zz
      real(kind_double),  intent(in) :: phase_ik
      real(kind_double), dimension(:,:) :: vfk_trans1, vfk_trans2, vfk_trans3
      logical, intent(in) :: desc_forces_local
      real(kind_double), intent(in) :: fact_random

      real(kind_double), intent(out) :: ktemp_norm
      real(kind_double), dimension(:,:), intent(out) :: dktemp_imm
      real(kind_double), dimension(:), intent(out) ::  vtemp0
      real(kind_double), dimension(:), allocatable :: vtemp1, vtemp2, vtemp3
      real(kind_double), dimension(3) :: dktemp_norm
      integer :: ja_neigh, dim_xdesc, inn
      real(kind_double) :: dtmp_jaik
      _NAMECURRENT_("inner_kernel_random_maha")
      _MLD_BEGIN_

      ja_neigh = size(vfk_trans1, dim=1)
      dim_xdesc= size(vfk_trans1, dim=2)

      dtmp_jaik = dot_product(xx, zz) + phase_ik
      ktemp_norm = fact_random * cos(dtmp_jaik)
      if (desc_forces_local) then
         allocate(vtemp1(ja_neigh), vtemp2(ja_neigh), vtemp3(ja_neigh))
         vtemp0(:) = 0.d0
         call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans1, ja_neigh, zz, 1, 0.d0, vtemp1, 1)
         call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans2, ja_neigh, zz, 1, 0.d0, vtemp2, 1)
         call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans3, ja_neigh, zz, 1, 0.d0, vtemp3, 1)
         do inn = 1, ja_neigh
            dktemp_norm(1) = -vtemp1(inn)*sin(dtmp_jaik) * fact_random
            dktemp_norm(2) = -vtemp2(inn)*sin(dtmp_jaik) * fact_random
            dktemp_norm(3) = -vtemp3(inn)*sin(dtmp_jaik) * fact_random
            vtemp0(:) = vtemp0(:) - dktemp_norm(1:3)
            !config_desc(iconf)%force_kernel(ik, ja, inn, 1:3) = dktemp_norm(1:3)
            dktemp_imm(1:3, inn) = dktemp_norm(1:3)
         end do
         !config_desc(iconf)%force_kernel(ik, ja, 0, 1:3) = vtemp0(1:3)
      end if


      _MLD_END_
   end subroutine inner_kernel_random_maha

end module kernel_tools


module module_compute_kernel

contains 

! subroutines designed for kernel implementation.
subroutine allocate_kernel(iconf, i_start_at, i_final_at, dim_kernel, imm, imm_neigh, desc_forces_local)
   use derived_types, only: config_desc
   use mld_logger
   implicit none
   integer, intent(in) :: iconf, imm, imm_neigh
   integer, intent(in) :: i_start_at, i_final_at, dim_kernel
   logical, intent(in) :: desc_forces_local

   _NAMECURRENT_("allocate_kernel")
   _MLD_BEGIN_

   if (allocated(config_desc(iconf)%energy_kernel)) deallocate (config_desc(iconf)%energy_kernel)
   allocate (config_desc(iconf)%energy_kernel(dim_kernel, imm))
   if (desc_forces_local) then
      if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
         if (allocated(config_desc(iconf)%force_kernel)) deallocate (config_desc(iconf)%force_kernel)
         allocate (config_desc(iconf)%force_kernel(1:dim_kernel, i_start_at:i_final_at, 0:imm_neigh, 1:3))
      end if
   end if

   _MLD_END_
end subroutine allocate_kernel

subroutine compute_kernel(i_start_at, i_final_at, iconf)
   ! Compute the kernel and it is stack as a descriptor with the dimension dim_kernel
   use module_kind_variables, ONLY: kind_double
   use ml_in_ndm_module, only: desc_forces
   use module_kernel, only: dim_kernel, global_kernel, draft_kernel, &
      length_kernel, sigma_kernel, &
      kernel_type, kernel_se, kernel_po, kernel_po_scaled, kernel_random, &
      kernel_random_maha, &
      kernel_maha, kernel_power, &
      kernel_phase_random, &
      kernel_random_po, &
      basis_random_po, norm_random_maha
   use derived_types, only: config_desc, config_real
   use module_Sigma_matrix, only: Sigma_sample_mcd_inv
#ifdef MLD_NDM
   use tab_imm_m_ml, ONLY: iwmax2
#else
   use ondm_tab_imm_m, ONLY: iwmax2
#endif
   use temporary_data_cov, only: dim_xdesc
   use kernel_tools, only: inner_kernel_se, inner_kernel_random_se, &
      inner_random_po, inner_kernel_po, inner_kernel_po_scaled, &
      inner_kernel_maha, inner_kernel_random_maha
   use mld_logger

   implicit none

   _NAMECURRENT_("compute_kernel")


   integer, intent(in)  :: i_start_at, i_final_at
   integer, optional    :: iconf
   integer  :: ja, iw2, ik, inn, ja_neigh
   real(kind_double), dimension(:), allocatable :: etemp, energy_ja, energy_ja_mean, &
      energy_ja_maha, kernel_ik
   real(kind_double)    :: l2kse, sigma_kse2, soverl, ktemp_norm, &
      dot_xx, ker_xx, dker_xx, sqrt_ker_xx, fact_random
   real(kind_double), dimension(3)  :: vtemp0
   real(kind_double), dimension(:, :), allocatable    :: vfk_trans1, vfk_trans2, vfk_trans3
   real(kind_double), dimension(:), allocatable :: vtemp1, vtemp2, vtemp3, vtemp4, vtemp_ja1, vtemp_ja2, vtemp_ja3
   real(kind_double), dimension(:, :), allocatable    :: dktemp_inn
   real(kind_double), external :: ddot

   logical  :: desc_forces_local


   _MLD_BEGIN_
   if ((i_start_at == 0) .and. (i_final_at == 0)) then
      config_desc(iconf)%energy_kernel(:, :) = 0.d0
      return
   end if
   desc_forces_local = desc_forces .and. (config_real(iconf)%has_force .or. config_real(iconf)%has_stress)
   config_desc(iconf)%energy_kernel(:, :) = 0.d0
   if (desc_forces_local) config_desc(iconf)%force_kernel(:, i_start_at:i_final_at, :, :) = 0.d0
#ifdef MLD_NDM
   iw2=0
!JPC iw2 inutilise ??
!!$  if (i_start_at == 1) iw2 = 0
!!$  if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
#else
   if (i_start_at == 1) iw2 = 0
   if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
#endif

   if (allocated(etemp)) deallocate (etemp); allocate (etemp(dim_xdesc))
   if (allocated(energy_ja)) deallocate (energy_ja); allocate (energy_ja(dim_xdesc))
   if (kernel_type == kernel_maha) then
      if (allocated(energy_ja_mean)) deallocate (energy_ja_mean); allocate (energy_ja_mean(dim_xdesc))
      if (allocated(energy_ja_maha)) deallocate (energy_ja_maha); allocate (energy_ja_maha(dim_xdesc))
   end if

   if (allocated(kernel_ik)) deallocate (kernel_ik); allocate (kernel_ik(dim_xdesc))
   l2kse = 2.d0*length_kernel**2
   sigma_kse2 = sigma_kernel**2
   soverl = sigma_kernel**2/length_kernel**2
   if ((kernel_type == kernel_random).or.(kernel_type==kernel_random_maha)) then
      if ((kernel_type == kernel_random)) fact_random = dsqrt( 2.d0 / dble(dim_kernel) ) * sigma_kernel
      if (kernel_type == kernel_random_maha) fact_random = dsqrt( 2.d0 / dble(dim_kernel * norm_random_maha) ) * sigma_kernel
   end if

   do ja = i_start_at, i_final_at
      !$! TODOkernel for mean kernels
      !$! if (.not. ((kernel_type == kernel_random) .or. (kernel_type == kernel_random_po))) then
      !$!   energy_ja(:) = (config_desc(iconf)%energy(:, ja) - min_ker(:))/(max_ker(:) - min_ker(:)) + min_ker(:)
      !$!   energy_ja(:) = (config_desc(iconf)%energy(:, ja) - mean_ker(:))/var_ker(:)
      !$! end if
      energy_ja(:) = config_desc(iconf)%energy(:, ja)

      if (kernel_type == kernel_po) then
         dot_xx = dot_product(energy_ja(:), energy_ja(:))
         ker_xx = (sigma_kse2 + dot_xx/l2kse)**kernel_power
         sqrt_ker_xx = dsqrt(ker_xx)
         dker_xx = 2.d0*dble(kernel_power)*(sigma_kse2 + dot_xx/l2kse)**(kernel_power - 1)/l2kse
      end if

      if (kernel_type == kernel_po_scaled) then
         dot_xx = dot_product(energy_ja(:), energy_ja(:)) 
      end if 

      !why_is_here? ja_neigh = 0
      ja_neigh = config_desc(iconf)%n_neigh(ja)
      if (desc_forces_local) then

         if (ja_neigh /= 0 ) then
            if (allocated(dktemp_inn)) deallocate (dktemp_inn); allocate (dktemp_inn(3, ja_neigh))
            if (allocated(vtemp1)) deallocate (vtemp1); allocate (vtemp1(ja_neigh))
            if (allocated(vtemp2)) deallocate (vtemp2); allocate (vtemp2(ja_neigh))
            if (allocated(vtemp3)) deallocate (vtemp3); allocate (vtemp3(ja_neigh))
            if (allocated(vtemp4)) deallocate (vtemp4); allocate (vtemp4(ja_neigh))
            !
            if (allocated(vfk_trans1)) deallocate (vfk_trans1); allocate (vfk_trans1(ja_neigh, dim_xdesc))
            if (allocated(vfk_trans2)) deallocate (vfk_trans2); allocate (vfk_trans2(ja_neigh, dim_xdesc))
            if (allocated(vfk_trans3)) deallocate (vfk_trans3); allocate (vfk_trans3(ja_neigh, dim_xdesc))
            do inn = 1, ja_neigh
               ! TODOkernel not forget There are those alternatives ...
               ! fk_trans(:) =(config_desc(iconf)%force(:,ja,inn,ix) - min_ker(:))/(max_ker(:)-min_ker(:)) + min_ker(:)
               ! fk_trans(:) =(config_desc(iconf)%force(:,ja,inn,ix) - mean_ker(:))/var_ker(:)
               vfk_trans1(inn, :) = config_desc(iconf)%force(:, ja, inn, 1)
               vfk_trans2(inn, :) = config_desc(iconf)%force(:, ja, inn, 2)
               vfk_trans3(inn, :) = config_desc(iconf)%force(:, ja, inn, 3)
            end do

            if ((kernel_type == kernel_po).or.( kernel_type == kernel_po_scaled)) then
               if (allocated(vtemp_ja1)) deallocate (vtemp_ja1); allocate (vtemp_ja1(ja_neigh))
               if (allocated(vtemp_ja2)) deallocate (vtemp_ja2); allocate (vtemp_ja2(ja_neigh))
               if (allocated(vtemp_ja3)) deallocate (vtemp_ja3); allocate (vtemp_ja3(ja_neigh))
               !write(*,*) 'ikkk', ja_neigh
               call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans1, ja_neigh, energy_ja, 1, 0.d0, vtemp_ja1, 1)
               call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans2, ja_neigh, energy_ja, 1, 0.d0, vtemp_ja2, 1)
               call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans3, ja_neigh, energy_ja, 1, 0.d0, vtemp_ja3, 1)
            end if

         end if ! ja_neigh /= 0

      end if

      do ik = 1, dim_kernel
         kernel_ik(:) = global_kernel(:, ik)
         if (ja_neigh == 0) then
            config_desc(iconf)%energy_kernel(ik, ja) = 0.d0
            if (desc_forces_local) config_desc(iconf)%force_kernel(ik, ja, 0, 1:3) = 0.d0
         else
            select case (kernel_type)


             case (kernel_random_maha)
               ! where I should put the length kernel and normalization factor ??????
               call inner_kernel_random_maha(energy_ja, kernel_ik, kernel_phase_random(ik), &
                  vfk_trans1, vfk_trans2, vfk_trans3, &
                  desc_forces_local, fact_random, &
                  ktemp_norm, dktemp_inn, vtemp0)

             case (kernel_random)

               !$! dtmp_jaik = dot_product(energy_ja(:), kernel_ik(:)) + kernel_phase_random(ik)
               !$! ktemp_norm = fact_random * cos(dtmp_jaik)
               !$!
               !$! if (desc_forces_local) then
               !$!   vtemp0(:) = 0.d0
               !$!   call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans1, ja_neigh, kernel_ik, 1, 0.d0, vtemp1, 1)
               !$!   call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans2, ja_neigh, kernel_ik, 1, 0.d0, vtemp2, 1)
               !$!   call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans3, ja_neigh, kernel_ik, 1, 0.d0, vtemp3, 1)
               !$!   do inn = 1, ja_neigh
               !$!     dktemp_norm(1) = -vtemp1(inn)*sin(dtmp_jaik) * fact_random
               !$!     dktemp_norm(2) = -vtemp2(inn)*sin(dtmp_jaik) * fact_random
               !$!     dktemp_norm(3) = -vtemp3(inn)*sin(dtmp_jaik) * fact_random
               !$!     vtemp0(:) = vtemp0(:) - dktemp_norm(1:3)
               !$!     config_desc(iconf)%force_kernel(ik, ja, inn, 1:3) = dktemp_norm(1:3)
               !$!   end do
               !$!
               !$!   config_desc(iconf)%force_kernel(ik, ja, 0, 1:3) = vtemp0(1:3)
               !$! end if


               call inner_kernel_random_se(energy_ja, kernel_ik, kernel_phase_random(ik), &
                  vfk_trans1, vfk_trans2, vfk_trans3, &
                  desc_forces_local, fact_random, &
                  ktemp_norm, dktemp_inn, vtemp0)
               !$! if (desc_forces_local) then
               !$!   config_desc(iconf)%force_kernel(ik, ja, 0, 1:3) = vtemp0(1:3)
               !$!   do inn=1, ja_neigh
               !$!     config_desc(iconf)%force_kernel(ik, ja, inn, 1:3) = dktemp_inn(1:3, inn)
               !$!   end do
               !$! end if

             case (kernel_random_po)

               !$! dim_omega = basis_random_po(ik)%dim_omega
               !$! if (allocated(ene_tmp)) deallocate (ene_tmp); allocate (ene_tmp(dim_omega))
               !$! dtmp = 1.d0
               !$! do ii = 1, dim_omega
               !$!   ene_tmp(ii) = dot_product(energy_ja(:), basis_random_po(ik)%omega(:, ii))
               !$!   dtmp = dtmp*ene_tmp(ii)
               !$! end do
               !$! ktemp_norm = dtmp/sqrt(dble(dim_kernel))
               !$! !dtmp_jaik = dot_product(energy_ja(:), kernel_ik(:)) + kernel_phase_random(ik)
               !$! if (desc_forces_local) then
               !$!   vtemp0(1:3) = 0.d0
               !$!   vtemp1(1:ja_neigh) = 0.d0
               !$!   vtemp2(1:ja_neigh) = 0.d0
               !$!   vtemp3(1:ja_neigh) = 0.d0
               !$!   if (allocated(force_tmp1)) deallocate (force_tmp1); allocate (force_tmp1(ja_neigh, dim_omega))
               !$!   if (allocated(force_tmp2)) deallocate (force_tmp2); allocate (force_tmp2(ja_neigh, dim_omega))
               !$!   if (allocated(force_tmp3)) deallocate (force_tmp3); allocate (force_tmp3(ja_neigh, dim_omega))
               !$!   do ii = 1, dim_omega
               !$!     kernel_ik(:) = basis_random_po(ik)%omega(:, ii)
               !$!     call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans1, ja_neigh, kernel_ik, 1, 0.d0, vtemp4, 1)
               !$!     force_tmp1(1:ja_neigh, ii) = vtemp4(1:ja_neigh)
               !$!     call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans2, ja_neigh, kernel_ik, 1, 0.d0, vtemp4, 1)
               !$!     force_tmp2(1:ja_neigh, ii) = vtemp4(1:ja_neigh)
               !$!     call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans3, ja_neigh, kernel_ik, 1, 0.d0, vtemp4, 1)
               !$!     force_tmp3(1:ja_neigh, ii) = vtemp4(1:ja_neigh)
               !$!   end do
               !$!
               !$!   do iid = 1, dim_omega
               !$!     dtmp = 1.d0
               !$!     do ii = 1, dim_omega
               !$!       if (ii == iid) cycle
               !$!       dtmp = dtmp*ene_tmp(ii)
               !$!     end do
               !$!     vtemp1(1:ja_neigh) = vtemp1(1:ja_neigh) + dtmp*force_tmp1(1:ja_neigh, iid)
               !$!     vtemp2(1:ja_neigh) = vtemp2(1:ja_neigh) + dtmp*force_tmp2(1:ja_neigh, iid)
               !$!     vtemp3(1:ja_neigh) = vtemp3(1:ja_neigh) + dtmp*force_tmp3(1:ja_neigh, iid)
               !$!   end do
               !$!
               !$!   do inn = 1, ja_neigh
               !$!     dktemp_norm(1) = vtemp1(inn)
               !$!     dktemp_norm(2) = vtemp2(inn)
               !$!     dktemp_norm(3) = vtemp3(inn)
               !$!     vtemp0(:) = vtemp0(:) - dktemp_norm(1:3)
               !$!     config_desc(iconf)%force_kernel(ik, ja, inn, 1:3) = dktemp_norm(1:3)/sqrt(dble(dim_kernel))
               !$!   end do
               !$!
               !$!   config_desc(iconf)%force_kernel(ik, ja, 0, 1:3) = vtemp0(1:3)/sqrt(dble(dim_kernel))
               !$! end if

               call inner_random_po(energy_ja, dim_kernel, basis_random_po(ik)%omega(:, :), &
                  vfk_trans1, vfk_trans2, vfk_trans3, desc_forces_local, &
                  ktemp_norm, dktemp_inn, vtemp0)
               !$! if (desc_forces_local) then
               !$!   config_desc(iconf)%force_kernel(ik, ja, 0, 1:3) = vtemp0(1:3)
               !$!   do inn=1, ja_neigh
               !$!     config_desc(iconf)%force_kernel(ik, ja, inn, 1:3) = dktemp_inn(1:3, inn)
               !$!   end do
               !$! end if




             case (kernel_po)
               !$! dtmp_jaik = dot_product(energy_ja(:), kernel_ik(:))
               !$! !dtmp_jaik = ddot(dim_xdesc, energy_ja(:), 1, kernel_ik(:), 1)
               !$! dtmp = sigma_kse2 + dtmp_jaik/l2kse
               !$! ktemp = dtmp**kernel_power
               !$! k_ik = (sigma_kse2 + dot_product(kernel_ik(:), kernel_ik(:))/l2kse)**kernel_power
               !$! k_ik2 = dsqrt(k_ik)
               !$! ktemp_norm = ktemp/(k_ik2*sqrt_ker_xx)
               !$! if (desc_forces_local) then
               !$!   vtemp0(:) = 0.d0
               !$!   tmp_f1 = dble(kernel_power)*dtmp**(kernel_power - 1)/l2kse
               !$!   call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans1, ja_neigh, kernel_ik, 1, 0.d0, vtemp1, 1)
               !$!   call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans2, ja_neigh, kernel_ik, 1, 0.d0, vtemp2, 1)
               !$!   call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans3, ja_neigh, kernel_ik, 1, 0.d0, vtemp3, 1)
               !$!   do inn = 1, ja_neigh
               !$!     dktemp(1) = tmp_f1*vtemp1(inn)
               !$!     dktemp(2) = tmp_f1*vtemp2(inn)
               !$!     dktemp(3) = tmp_f1*vtemp3(inn)
               !$!     !dktemp = \nabla_b ktilde (x,x) = \partial(ktilde(x,x))/\partial(x)  \nabla_b x
               !$!     dktemp_ja(1) = vtemp_ja1(inn)*dker_xx
               !$!     dktemp_ja(2) = vtemp_ja2(inn)*dker_xx
               !$!     dktemp_ja(3) = vtemp_ja3(inn)*dker_xx
               !$!     dktemp_norm(1:3) = (dktemp(1:3) - ktemp*dktemp_ja(1:3)/(2.d0*ker_xx))/(sqrt_ker_xx*k_ik2)
               !$!     config_desc(iconf)%force_kernel(ik, ja, inn, 1:3) = dktemp_norm(1:3)
               !$!     vtemp0(:) = vtemp0(:) - dktemp_norm(1:3)
               !$!   end do
               !$!   config_desc(iconf)%force_kernel(ik, ja, 0, 1:3) = vtemp0(1:3)
               !$! end if                  ! desc_forces_local

               !!! dktemp_inn(:,:) = 0.d0
               !!! vtemp0(:) = 0.d0
               call inner_kernel_po(energy_ja, kernel_ik, vfk_trans1, vfk_trans2, vfk_trans3, &
                  vtemp_ja1,  vtemp_ja2,  vtemp_ja3, desc_forces_local, sqrt_ker_xx, dker_xx, ker_xx, &
                  sigma_kse2, l2kse, kernel_power, &
                  ktemp_norm, dktemp_inn,   vtemp0)
               !$! if (desc_forces_local) then
               !$!   do inn = 1, ja_neigh
               !$!     config_desc(iconf)%force_kernel(ik, ja, inn, 1:3) = dktemp_inn(1:3, inn)
               !$!   end do
               !$!   config_desc(iconf)%force_kernel(ik, ja, 0, 1:3) = vtemp0(1:3)
               !$! end if

             case (kernel_po_scaled)
                 call inner_kernel_po_scaled(energy_ja, kernel_ik, vfk_trans1, vfk_trans2, vfk_trans3, &
                   vtemp_ja1, vtemp_ja1, vtemp_ja2, desc_forces_local, & 
                   dot_xx, sigma_kse2, kernel_power, ktemp_norm, dktemp_inn, vtemp0)  

             case (kernel_se)

               !$! !debug etemp(:) = (energy_ja(:) - mean_ker(:))/var_ker(:) - kernel_ik(:)
               !$! etemp(:) = energy_ja(:) - draft_kernel(:, ik)
               !$! dtmp = dot_product(etemp(:), etemp(:))
               !$! ktemp = dexp(-dtmp/l2kse)
               !$! ktemp_norm = sigma_kse2*ktemp
               !$!
               !$! if (desc_forces_local) then
               !$!   vtemp0(:) = 0.d0
               !$!   call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans1, ja_neigh, etemp, 1, 0.d0, vtemp1, 1)
               !$!   call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans2, ja_neigh, etemp, 1, 0.d0, vtemp2, 1)
               !$!   call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans3, ja_neigh, etemp, 1, 0.d0, vtemp3, 1)
               !$!   tmp_f1 = -soverl*ktemp
               !$!   do inn = 1, ja_neigh
               !$!     dktemp(1) = tmp_f1*vtemp1(inn)
               !$!     dktemp(2) = tmp_f1*vtemp2(inn)
               !$!     dktemp(3) = tmp_f1*vtemp3(inn)
               !$!     dktemp_norm(1:3) = dktemp(1:3)
               !$!     config_desc(iconf)%force_kernel(ik, ja, inn, 1:3) = dktemp_norm(1:3)
               !$!     vtemp0(:) = vtemp0(:) - dktemp_norm(1:3)
               !$!   end do
               !$!   config_desc(iconf)%force_kernel(ik, ja, 0, 1:3) = vtemp0(1:3)
               !$! end if

               !debug etemp(:) = (energy_ja(:) - mean_ker(:))/var_ker(:) - kernel_ik(:)
               etemp(:) = energy_ja(:) - draft_kernel(:, ik)
               call inner_kernel_se(etemp,   &
                  vfk_trans1, vfk_trans2, vfk_trans3, &
                  desc_forces_local, soverl, sigma_kse2, l2kse, &
                  ktemp_norm, dktemp_inn,   vtemp0)
               !$! if (desc_forces_local) then
               !$!    config_desc(iconf)%force_kernel(ik, ja, 0, 1:3) = vtemp0(1:3)
               !$!    do inn = 1, ja_neigh
               !$!      config_desc(iconf)%force_kernel(ik, ja, inn, 1:3) = dktemp_inn(1:3, inn)
               !$!    end do
               !$! end if


             case (kernel_maha)
               !$! energy_ja_mean(:) = energy_ja(:) - draft_kernel(:, ik)
               !$! energy_ja_maha(:) = matmul(Sigma_sample_mcd_inv(:, :), energy_ja_mean(:))
               !$! ! energy_ja_maha(:) = energy_ja_mean(:)
               !$! tmp_f1 = dot_product(energy_ja_mean(:), energy_ja_maha(:))
               !$! dtmp_ja_ik = tmp_f1**kernel_power
               !$! dtmp_ja_ik = dexp(-tmp_f1**kernel_power)
               !$! ! k_ik = maha_norm_kernel(ik) ** kernel_power
               !$! ! k_ik2=dsqrt(k_ik)
               !$! ktemp_norm = dtmp_ja_ik ! / (sqrt_ker_xx*k_ik2)
               !$! if (desc_forces_local) then
               !$!   vtemp0(:) = 0.d0
               !$!   call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans1, ja_neigh, energy_ja_maha, 1, 0.d0, vtemp1, 1)
               !$!   call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans2, ja_neigh, energy_ja_maha, 1, 0.d0, vtemp2, 1)
               !$!   call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans3, ja_neigh, energy_ja_maha, 1, 0.d0, vtemp3, 1)
               !$!   tmp_f2 = 2.d0*dble(kernel_power)*tmp_f1**(kernel_power - 1)
               !$!   tmp_f2 = -2.d0*dble(kernel_power)*tmp_f1**(kernel_power - 1)*dexp(-tmp_f1**kernel_power)
               !$!   do inn = 1, ja_neigh
               !$!     dktemp(1) = tmp_f2*vtemp1(inn)
               !$!     dktemp(2) = tmp_f2*vtemp2(inn)
               !$!     dktemp(3) = tmp_f2*vtemp3(inn)
               !$!     dktemp_norm(1:3) = dktemp(1:3)                   ! /(sqrt_ker_xx*k_ik2) - ktemp_norm * dktemp_ja(1:3) / (2.d0*k_ja32)
               !$!     config_desc(iconf)%force_kernel(ik, ja, inn, 1:3) = dktemp_norm(1:3)
               !$!     vtemp0(:) = vtemp0(:) - dktemp_norm(1:3)
               !$!   end do
               !$!   config_desc(iconf)%force_kernel(ik, ja, 0, 1:3) = vtemp0(1:3)
               !$! end if

               energy_ja_mean(:) = energy_ja(:) - draft_kernel(:, ik)
               !energy_ja_mean(:) = energy_ja(:) - kernel_ik(:)
               call inner_kernel_maha(energy_ja_mean, Sigma_sample_mcd_inv, &
                  vfk_trans1, vfk_trans2, vfk_trans3, kernel_power, desc_forces_local, &
                  ktemp_norm, dktemp_inn, vtemp0)



            end select

            if (desc_forces_local) then
               config_desc(iconf)%force_kernel(ik, ja, 0, 1:3) = vtemp0(1:3)
               do inn = 1, ja_neigh
                  config_desc(iconf)%force_kernel(ik, ja, inn, 1:3) = dktemp_inn(1:3, inn)
               end do
            end if

            config_desc(iconf)%energy_kernel(ik, ja) = ktemp_norm
         end if ! ja_neigh ....


      end do                  ! ik - kernel
      !$! if (rangml==0) write (667,'(10e20.7, i8)') config_desc(iconf)%energy_kernel(1:10,ja), rangml
      !$! if (rangml==3) write (668,'(10e20.7, i8)') config_desc(iconf)%energy_kernel(1:10,ja), rangml
   end do                  ! ja - atoms

   _MLD_END_
end subroutine

end module module_compute_kernel



subroutine init_sample_kernel()
   use module_kind_variables, only: kind_double
   use module_kernel, only: dim_kernel, np_kernel_full, global_kernel, kernel_phase_random, &
      sigma_kernel, length_kernel, &
      krff_type
   use temporary_data_cov, ONLY: dim_xdesc
   use module_sample_rand_ker, only: get_rff_val_sigma
   use mld_logger

   implicit none
   real(kind_double)    :: mu

   _NAMECURRENT_("init_sample_kernel")

   _MLD_BEGIN_
   dim_kernel = np_kernel_full
   if (allocated(global_kernel)) deallocate (global_kernel); allocate (global_kernel(dim_xdesc, dim_kernel))

   ! --------->
   ! krff_type : type of distribution
   ! mu : \mu of distribution
   ! the random omega's: D x K size
   ! length_kernel : defines sigma**2 = length_kernel**2
   ! global_kernel is the place of omega's
   ! kernel_phase_random is the place of b's
   mu = 0.d0
   call  get_rff_val_sigma(krff_type, mu, dim_xdesc, dim_kernel, global_kernel, kernel_phase_random)
   global_kernel = global_kernel/(length_kernel*dsqrt(2.d0))

   call log_info('init_sample_kernel: length_kernel  = '//vtoa(length_kernel))
   call log_info('init_sample_kernel: sigma_kernel = '//vtoa(sigma_kernel))

   _MLD_END_

end subroutine init_sample_kernel

subroutine init_sample_kernel_maha()
   use module_kind_variables, only: kind_double
   use module_kernel, only: dim_kernel, np_kernel_full, global_kernel, kernel_phase_random, &
      sigma_kernel, length_kernel, &
      krff_type, draft_kernel, norm_random_maha
   use temporary_data_cov, ONLY: dim_xdesc
   use module_sample_rand_ker, only: get_rff_sigma
   use module_covariance, only : train_covariance_matrix
   use mod_covariance_matrix, only: obj_cov
   use ml_in_ndm_module, only: desc_forces, two_pi
   use module_Sigma_matrix, only: Sigma_sample_mcd, Sigma_sample_mcd_inv
   use math, only: serial_pseudo_inverse, serial_determinant_symmetric_general, &
      serial_determinant_symmetric_positive, ComputeDeterminantLU, &
      ComputeDeterminantQR, serial_svd_filter
   use mld_logger


   implicit none
   real(kind_double), dimension(:), allocatable  :: mu
   real(kind_double), dimension(:,:), allocatable  :: sigma, sigma_filter
   integer :: rank_cov
   real(kind_double) :: determinant_sigma, determinant_sigma_filter

   logical :: desc_forces_local
   character(len=200)  :: chlog
   _NAMECURRENT_("init_sample_kernel_maha")

   _MLD_BEGIN_

   train_covariance_matrix = .true.
   desc_forces_local = desc_forces
   desc_forces = .false.
   call main_compute_descriptors
   desc_forces = desc_forces_local

   allocate(mu(dim_xdesc))
   allocate(sigma(dim_xdesc, dim_xdesc))
   if (dim_xdesc /= size(mu)) stop 'init_sample_kernel_maha 01 '
   !if (dim_xdesc /= size(obj_cov%matrix, 1)) stop 'init_sample_kernel_maha 02 '

   call obj_cov%get_mu(mu)

   call obj_cov%get_matrix(sigma)

   if (allocated(Sigma_sample_mcd)) deallocate(Sigma_sample_mcd) ;   allocate(Sigma_sample_mcd(dim_xdesc, dim_xdesc))
   if (allocated(Sigma_sample_mcd_inv)) deallocate(Sigma_sample_mcd_inv) ;  allocate(Sigma_sample_mcd_inv(dim_xdesc, dim_xdesc))
   Sigma_sample_mcd = sigma
   call serial_pseudo_inverse(Sigma_sample_mcd, rank_cov, Sigma_sample_mcd_inv)

   call log_info("ML: ... sample covariance has the rank ."//vtoa(rank_cov))


   dim_kernel = np_kernel_full
   if (allocated(global_kernel)) deallocate (global_kernel); allocate (global_kernel(dim_xdesc, dim_kernel))
   if (allocated(draft_kernel)) deallocate (draft_kernel); allocate (draft_kernel(dim_xdesc, dim_kernel))

   ! --------->
   ! krff_type : type of distribution
   ! mu : \mu of distribution
   ! the random omega's: D x K size
   ! length_kernel : defines sigma**2 = length_kernel**2
   ! global_kernel is the place of omega's
   ! kernel_phase_random is the place of b's
   !mu = 0.d0
   mu=0.d0
   !$! sigma(:,:) =  0.d0
   !$! do ii = 1, dim_xdesc
   !$!   sigma(ii, ii) = 1.d0
   !$! end do
   call serial_svd_filter(sigma, rank_cov, sigma_filter, determinant_sigma, determinant_sigma_filter)

   call log_info("ML: ... sample covariance determinant brut   ->  "//vtoa(determinant_sigma))
   call log_info("ML: ... sample covariance determinant filter ->  "//vtoa(determinant_sigma_filter))
   call  ComputeDeterminantQR(sigma, dim_xdesc, determinant_sigma)
   call log_info("ML: ... sample covariance determinant 1 "//vtoa(determinant_sigma))

   call ComputeDeterminantLU(sigma, dim_xdesc, determinant_sigma)
   call log_info("ML: ... sample covariance determinant 2 "//vtoa(determinant_sigma))
   !sigma = sigma / determinant_sigma
   !sigma_filter = sigma_filter / determinant_sigma_filter
   call get_rff_sigma(krff_type, mu, sigma, dim_xdesc, dim_kernel, global_kernel, kernel_phase_random)
   norm_random_maha = dsqrt(determinant_sigma)*(length_kernel/two_pi)**dble(dim_xdesc/2)
   write(chlog, '(es20.10)') norm_random_maha
   call log_info("ML: ... norm Fourier factor  "//trim(chlog))
   norm_random_maha=1.d0

   draft_kernel = global_kernel
   global_kernel = global_kernel/dsqrt(length_kernel)
   write(chlog, '(es20.10)') length_kernel
   call log_info('init_sample_kernel_maha: length_kernel  = '//trim(chlog))
   write(chlog, '(es20.10)') sigma_kernel
   call log_info('init_sample_kernel_maha: sigma_kernel = '//trim(chlog))

   _MLD_END_

end subroutine init_sample_kernel_maha

subroutine init_sample_kernel_po()
   !use mpi
   use mld_mpi, only: comm_mld
   use ml_in_ndm_module, only: rangml
   use module_kind_variables, only: kind_double
   use module_kernel, only: dim_kernel, np_kernel_full, np_omega, global_kernel, basis_random_po
   use temporary_data_cov, ONLY: dim_xdesc
   use mld_logger
   implicit none
   real(kind_double), dimension(:), allocatable :: tmp_kernel
   integer  :: ik, io, ii
   integer  :: dim_omega
   _NAMECURRENT_("init_sample_kernel_po")
   _MLD_BEGIN_
   dim_kernel = np_kernel_full
   dim_omega = np_omega

   if (allocated(basis_random_po)) deallocate (basis_random_po); allocate (basis_random_po(dim_kernel))
   do ik = 1, dim_kernel
      if (allocated(basis_random_po(ik)%omega)) deallocate (basis_random_po(ik)%omega)
      allocate (basis_random_po(ik)%omega(dim_xdesc, dim_omega))
      basis_random_po(ik)%dim_omega = dim_omega
      basis_random_po(ik)%dim_xdesc = dim_xdesc
   end do
   ! allocate for compatibility reason. But no need.
   if (allocated(global_kernel)) deallocate (global_kernel); allocate (global_kernel(dim_xdesc, dim_kernel))
   global_kernel = 0.d0
   if (allocated(tmp_kernel)) deallocate (tmp_kernel); allocate (tmp_kernel(dim_xdesc))

   !TODOkernel parameters ...
   !TODOkernel mu = 0.d0
   !TODOkernel sigma = sigma_kernel
   do ik = 1, dim_kernel
      do io = 1, dim_omega
         if (rangml == 0) then
            call random_number(tmp_kernel)
            do ii = 1, dim_xdesc
               if (tmp_kernel(ii) < 0.5d0) then
                  tmp_kernel(ii) = -1.d0
               else
                  tmp_kernel(ii) = 1.d0
               end if
            end do
         end if
         !TORC! call MPI_BCAST(tmp_kernel, size(tmp_kernel, 1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, codeml)
         call comm_mld%bcast(0, tmp_kernel)
         basis_random_po(ik)%omega(:, io) = tmp_kernel(:)
      end do
   end do

   _MLD_END_
end subroutine init_sample_kernel_po

subroutine init_kernel()
   !---------------------------------------------------------------------
   !   For kernel other than the ones based on Covariance matrix
   !   Guassian, Polynomials etc
   !---------------------------------------------------------------------
   !   draft_kernel    dim_xdesc, dim_kernel -> kernel as it is read
   !                                            from disk
   !   global_kernel   dim_xdesc, dim_kernel -> kernel from which is
   !                          substracted mean and devided covariance
   !---------------------------------------------------------------------
   !   For Covariance based kernel
   !---------------------------------------------------------------------
   !  maha_kernel      dim_xdesc, dim_kernel ->  Sigma-1 x (z - mean)
   !  maha_norm_points_kernel   dim_kernel   ->  Maha norm with respect Sigma
   !---------------------------------------------------------------------

   use module_kind_variables, only: kind_double
   use ml_in_ndm_module, only: rangml
   use temporary_data_cov, only: dim_xdesc
   use module_kernel, only: dim_kernel, draft_kernel, global_kernel, info_kernel, &
      draft_kernel, min_ker, max_ker, mean_ker, var_ker, &
      kernel_type, kernel_maha, time_for_desc_kernel, &
      maha_kernel, maha_norm_kernel
   use module_Sigma_matrix, only: Sigma_sample_mcd, Sigma_sample_mcd_inv
   use module_end_ml, only: end_ml
   use module_ml_scalapack, only: scalapack_driver, context
   use math, only: serial_pseudo_inverse
   use mld_logger
   use mld_string
   use mld_mpi, only: comm_mld

   implicit none

   _NAMECURRENT_("init_kernel")


   logical  :: ok
   integer   :: kunit
   integer  :: dim_tmp, dim_tmp2, iconf_temp, ia_temp, itemp, ii
   character(len=15)    :: filename_temp
   character(len=80)    :: CHFMT, text
   real(kind_double), dimension(dim_xdesc)      :: temp
   real(kind_double)    :: atemp, factn
   integer :: rank_kernel



   _MLD_BEGIN_
   time_for_desc_kernel = 0.d0
   if (rangml == 0) then
      inquire (file="kernel_matrix.dat", exist=ok)
      if (ok) then
         open (file='kernel_matrix.dat', newunit=kunit, action='read', status='unknown')
         write (CHFMT, *) '( i6,', dim_xdesc, 'e20.10, i6, i6, a15 )'
         read (kunit, *) dim_kernel, dim_tmp
      else
         text = "the file kernel_matrix.dat is not there"
      end if
   end if

   call comm_mld%bcast(0, ok)
   !WTF! if(ok .eqv. .true.) then
   !WTF!    call end_ml(text, scalapack_driver, context, rangml)
   !WTF! end if

   !TORC! call MPI_BCAST(dim_kernel, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, codeml)
   call comm_mld%bcast(0, dim_kernel)

   if (rangml == 0) then
      write (6, '("ML: kernel shape ................", 2i7)') dim_kernel, dim_tmp - 4
      if (dim_xdesc /= (dim_tmp - 4)) then
         if (rangml == 0) write (6, *) 'dim descriptor', dim_xdesc, "dim recorded kernel", dim_tmp
         text = " ERROR the kernel recorded and descriptor not the same  dimension"
         call end_ml(text, scalapack_driver, context, rangml)
      end if
   end if

   if (allocated(global_kernel)) deallocate (global_kernel); allocate (global_kernel(dim_xdesc, dim_kernel))
   if (allocated(draft_kernel)) deallocate (draft_kernel); allocate (draft_kernel(dim_xdesc, dim_kernel))
   if (allocated(info_kernel)) deallocate (info_kernel); allocate (info_kernel(dim_kernel))
   if (allocated(mean_ker)) deallocate (mean_ker); allocate (mean_ker(dim_xdesc))
   if (allocated(var_ker)) deallocate (var_ker); allocate (var_ker(dim_xdesc))
   if (allocated(min_ker)) deallocate (min_ker); allocate (min_ker(dim_xdesc))
   if (allocated(max_ker)) deallocate (max_ker); allocate (max_ker(dim_xdesc))

   do ii = 1, dim_kernel
      if (rangml == 0) then
         read (kunit, CHFMT) itemp, temp(1:dim_xdesc), iconf_temp, ia_temp, filename_temp
      end if

      !TORC! call MPI_BCAST(itemp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, codeml)
      call comm_mld%bcast(0, itemp)
      !TORC! call MPI_BCAST(iconf_temp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, codeml)
      call comm_mld%bcast(0, iconf_temp)
      !TORC! call MPI_BCAST(ia_temp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, codeml)
      call comm_mld%bcast(0, ia_temp)
      !TORC! call MPI_BCAST(filename_temp, len(filename_temp), MPI_CHARACTER, 0, MPI_COMM_WORLD, codeml)
      call comm_mld%bcast(0, filename_temp)
      !TORC! call MPI_BCAST(temp, size(temp, 1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, codeml)
      call comm_mld%bcast(0, temp)

      draft_kernel(1:dim_xdesc, itemp) = temp(1:dim_xdesc)
      info_kernel(itemp)%ia = ia_temp
      info_kernel(itemp)%iconf = iconf_temp
      info_kernel(itemp)%filename = trim(filename_temp)
   end do

   if (rangml == 0) then
      close (kunit, status="keep")
      call log_info("ML: ... kernel read.")
   end if

   ! mean of the kernel
   mean_ker(:) = sum(draft_kernel(:, :), dim=2)/size(draft_kernel, 2)
   do ii = 1, dim_xdesc
      atemp = minval(draft_kernel(ii, :))
      min_ker(ii) = atemp
      atemp = maxval(draft_kernel(ii, :))
      max_ker(ii) = atemp
   end do

   ! variance of the kernel
   var_ker(:) = 0.d0
   do ii = 1, dim_kernel
      var_ker(:) = var_ker(:) + (draft_kernel(:, ii) - mean_ker(:))**2
   end do
   var_ker(:) = dsqrt(var_ker(:)/dble(dim_kernel))

   !$! ! min max rescaling
   !$! do ii = 1, dim_kernel
   !$!   global_kernel(:, ii) = (draft_kernel(:, ii) - min_ker(:))/(max_ker(:) - min_ker(:)) + min_ker(:)
   !$! end do

   ! variance rescaling
   if ((kernel_type == kernel_maha)) then
      do ii = 1, dim_kernel
         global_kernel(:, ii) = draft_kernel(:, ii) - mean_ker(:)
      end do
   else
      do ii = 1, dim_kernel
         global_kernel(:, ii) = (draft_kernel(:, ii) - mean_ker(:))/var_ker(:)
      end do
   end if


   if (kernel_type == kernel_maha) then

      ! This will be replaced by a proper computation from the kernel matrix
      if (allocated(Sigma_sample_mcd_inv)) deallocate (Sigma_sample_mcd_inv); allocate (Sigma_sample_mcd_inv(dim_xdesc, dim_xdesc))
      if (allocated(Sigma_sample_mcd)) deallocate (Sigma_sample_mcd); allocate (Sigma_sample_mcd(dim_xdesc, dim_xdesc))
      if (rangml == 0) then
         open (file='inverse_Sigma_mcd_matrix.mat', newunit=kunit, action='read', status='unknown')
         read (kunit, *) dim_tmp, dim_tmp2
         if (dim_tmp /= dim_xdesc) then
            write (6, *) "The MCD inverse covariance matrix was generated using diffrent descriptor dim_tmp readed versus dim_xdesc", dim_tmp, dim_xdesc
            stop 'in init_kernel read error inverse covariance'
         end if

         write (CHFMT, *) '(', dim_xdesc, 'e25.15)'
         do ii = 1, size(Sigma_sample_mcd_inv, 2)
            read (kunit, CHFMT) Sigma_sample_mcd_inv(:, ii)
         end do
         close (kunit)

         open (file='Sigma_mcd_matrix.mat', newunit=kunit, action='read', status='unknown')
         read (kunit, '(2i9)') dim_tmp, dim_tmp2

         if (dim_tmp /= dim_xdesc) then
            write (6, *) "The MCD covariance matrix was generated using diffrent descriptor dim_tmp readed versus dim_xdesc", dim_tmp, dim_xdesc
            stop 'in init_kernel read error covariance'
         end if
         write (CHFMT, *) '(', dim_xdesc, 'e25.15)'
         do ii = 1, size(Sigma_sample_mcd, 2)
            read (kunit, CHFMT) Sigma_sample_mcd(:, ii)
         end do
         close (kunit)
      end if                  ! rangml==0

      !TORC! call MPI_BCAST(Sigma_sample_mcd, dim_xdesc**2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, codeml)
      call comm_mld%bcast(0, Sigma_sample_mcd)
      !TORC! call MPI_BCAST(Sigma_sample_mcd_inv, dim_xdesc**2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, codeml)
      call comm_mld%bcast(0, Sigma_sample_mcd_inv)

      if (allocated(maha_kernel)) deallocate (maha_kernel); allocate (maha_kernel(dim_xdesc, dim_kernel))
      if (allocated(maha_norm_kernel)) deallocate (maha_norm_kernel); allocate (maha_norm_kernel(dim_kernel))
      ! ( (DxD) matrix with (KxD)^T)^T
      maha_kernel(1:dim_xdesc, 1:dim_kernel) = matmul(Sigma_sample_mcd_inv(1:dim_xdesc, 1:dim_xdesc), global_kernel(1:dim_xdesc, 1:dim_kernel))
      factn= dsqrt(dble(dim_kernel))
      call dgemm('N', 'T', dim_xdesc, dim_xdesc, dim_kernel, factn, global_kernel, dim_xdesc, global_kernel, dim_xdesc, 0.d0, Sigma_sample_mcd, dim_xdesc)
      call serial_pseudo_inverse(Sigma_sample_mcd, rank_kernel, Sigma_sample_mcd_inv)
      call log_info("ML: ... sample covariance of the kernel maha has the rank ."//vtoa(rank_kernel))

      do ii = 1, dim_kernel
         maha_norm_kernel(ii) = dot_product(global_kernel(ii, 1:dim_xdesc), maha_kernel(1:dim_xdesc, ii))
      end do

   end if                  ! kernel_type == kernel_maha


   _MLD_END_
end subroutine init_kernel

!subroutine compute_length()

 !use derived_types, only: config_desc
 !use module_kernel, only: length_kernel_k, sigma_kse2_k, soverl_k, multiple_kernel_length

 !implicit none
 !integer :: length_desc
 !multiple_kernel_length = .true.
 !length_desc = SIZE(config_desc%energy, dim = 2)


 !if (allocated(length_kernel_k)) deallocate (length_kernel_k); allocate (length_kernel_k(length_desc))
 !if (allocated(sigma_kse2_k)) deallocate (sigma_kse2_k); allocate (sigma_kse2_k(length_desc))
 !if (allocated(sigma_kse2_k)) deallocate (sigma_kse2_k); allocate (sigma_kse2_k(length_desc))

 !if (length_kernel == 0) then
 !call compute_length()
 !end if
!end subroutine compute_length
