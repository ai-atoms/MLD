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


module compute_bispectrum_so4_mod
   use module_kind_variables, ONLY: kind_double

   implicit none

contains


   ! subroutine compute_bispectrum_so4(i_start_at,i_final_at,d_n_neigh, d_kind_neigh, local_bispectrum_so4_out,local_bispectrum_so4_deriv_out, iconf)
   subroutine compute_bispectrum_so4(i_start_at, i_final_at, d_n_neigh, d_kind_neigh, iconf)
     use mld_logger
#ifdef MLD_NDM
      use gen_com_m, ONLY: A2cm, lperiod
      use gen_com_m_ml, ONLY: imm, bg, at
      use tab_imm_m_ml, ONLY: iwmax2, xp
#else
      use ondm_gen_com_m, ONLY: imm, lperiod
      use ondm_tab_imm_m, ONLY: iwmax2, xp
#endif
      use angular_functions, only: spherical_4d
      use ml_in_ndm_module, ONLY: jj_max, imm_neigh, lmask, &
         weighted, weighted_3ch,  desc_forces,  linvisible, Nfix, fix_Nmax_neigh
      use derived_types, only: config_real, config_desc

      use module_bispectrum_so4, only: lbso4_diag, czero, ZAcmm, ZBcmm, ZCcmm, ZAcmm2, ZBcmm2, ZCcmm2, ZAcmm3, ZBcmm3, ZCcmm3, &
         cmm, cmm2, cmm3, &      ! dcmm, dcmm2, dcmm3, &
         bisso4_dim, &
         class_mml, &
         fcut, dfcut, fcut_w, dfcut_w, fcut_w_3ch, dfcut_w_3ch
      use module_neigh_local, only: type_fcut, r_cut, r_cut_width, r_cut_in, r_cut_width_in, r_cut_pair_in, fcut_rij, &
         r_central, i_central, i_type, i_type_db, tmp_dxp, tmp_xp, max_neigh_local, out_max_neigh_local, &
         preallocate_neigh_ja, build_local_neighbours_ja, build_local_neighbours_ja_Nfix, &
         reallocate_neigh_ja, build_neigh_ja_type, iw2 

      use time_check_general, only: time, tot_time, debug_time, MY_MPI_WTIME
      use mesh_grid, only : linear_grid
#ifdef MLD_NDM
  use notperiod_mod
#else
      use ondm_transform_coord, only: ondm_notperiod
#endif

      integer, intent(in)  :: i_start_at, i_final_at
      integer, dimension(imm), intent(out)   :: d_n_neigh
      integer, dimension(imm, imm_neigh), intent(out)    :: d_kind_neigh
      !double precision,dimension(bisso4_dim,imm),intent(out) :: config_desc(iconf)%energy
      !double precision,dimension(bisso4_dim,imm, 0:imm_neigh, 3),intent(out) :: config_desc(iconf)%force
      integer, optional    :: iconf

      integer  :: n_axpy
      integer, parameter   :: incx = 1, incy = 1
      logical  :: small

      double complex, dimension(-jj_max:jj_max, -jj_max:jj_max, 0:jj_max)  :: Umm
      double complex, dimension(-jj_max:jj_max, -jj_max:jj_max, 0:jj_max, 3)     :: dUmm
      real(kind_double), dimension(:, :), allocatable   :: xpnp
      real(kind_double), dimension(3) :: dxp_ji
      integer  :: j, iw

      integer  :: l, l1, l2, m1, m2
      integer  :: ia, ja, i_bi, i_bi2, i_bi3, ia_n
      real(kind_double)     :: r2_ji, r_ji, etmp, etmp2, etmp3, rcut_fix
      double complex :: sum_bi, sum_bi_w, sum_bi_w_3ch
      !double complex,dimension(3) :: dsum_bi, dsum_bi_w, dsum_bi_w_3ch
      double complex, allocatable, dimension(:)    :: dsum_bi_ia, dsum_bi_ia_w, dsum_bi_ia_w_3ch
      real(kind_double)    :: factor_ia, factor_ia_3ch, factor_ja, factor_ja_3ch, &
         ta, tb, tmpll
      integer  :: l1_min, l1_max, it_cg, max_neigh, ja_atom 
      real(kind_double) :: fcut_in, dfcut_in, fcut_out, dfcut_out
      real(kind_double) :: local_r_cut_in
      logical  :: desc_forces_local
      integer, dimension(:), allocatable :: i_central_Nfix
      _NAMECURRENT_("compute_bispectrum_so4")
      _MLD_BEGIN_


      max_neigh = imm_neigh
      desc_forces_local = desc_forces .and. (config_real(iconf)%has_force .or. config_real(iconf)%has_stress)

      if ((i_start_at == 0) .and. (i_final_at == 0)) then
         d_n_neigh(:) = 0
         d_kind_neigh(:, :) = 0
         config_desc(iconf)%energy(:, :) = 0.d0
         !  config_desc(iconf)%force(:,:,:,:)=0.d0
         return
      end if


#ifdef MLD_NDM
    small = .true.
#else
      small = .false.
#endif 
      if (present(iconf)) then
         small = config_real(iconf)%small
      end if
!JPC : XPNP INUTILE
      ALLOCATE (xpnp(3, imm))
      if (lperiod) then
         xpnp(:, :) = xp(:, :)
      else
#ifdef MLD_NDM
      call notperiod(imm,xp, xpnp,at,bg,.false.)
#else
         call ondm_notperiod(xp, xpnp)
#endif
      end if
      !call cryst_to_cart (imm, xpnp, bg, -1)

      d_n_neigh(:) = 0
      d_kind_neigh(:, :) = 0
      config_desc(iconf)%energy(:, :) = 0.d0
      if (desc_forces_local) then
         config_desc(iconf)%force(:, :, :, :) = 0.d0
         !if (allocated(dcmm)) deallocate(dcmm) ; allocate(dcmm(-jj_max:jj_max,-jj_max:jj_max, 0:jj_max, 0:imm_neigh,1:3))
         if (allocated(dsum_bi_ia)) deallocate (dsum_bi_ia); allocate (dsum_bi_ia(imm_neigh*3))
         if (weighted) then
            ! if (allocated(dcmm2)) deallocate(dcmm2) ; allocate(dcmm2(-jj_max:jj_max,-jj_max:jj_max,0:jj_max,0:imm_neigh, 1:3))
            if (allocated(dsum_bi_ia_w)) deallocate (dsum_bi_ia_w); allocate (dsum_bi_ia_w(imm_neigh*3))
         end if
         if (weighted_3ch) then
            ! if (allocated(dcmm3)) deallocate(dcmm3) ; allocate(dcmm3(-jj_max:jj_max, -jj_max:jj_max, 0:jj_max, 0:imm_neigh,1:3))
            if (allocated(dsum_bi_ia_w_3ch)) deallocate (dsum_bi_ia_w_3ch); allocate (dsum_bi_ia_w_3ch(imm_neigh*3))
         end if

      end if

      !allocate new object for neighbours...
      call preallocate_neigh_ja(r_central, i_type, i_central, tmp_dxp, tmp_xp, imm_neigh)

      !if (i_start_at == 1) iw2 = 0
      !if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
      if (i_start_at == 1) iw2 = 0
#ifdef MLD_NDM
!JPC IW2 est l'indice du preimer voisin ca vaut toujours 1 pour small
!  if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
      iw2=0  
!JPC 
#else
      if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
#endif

      if (debug_time) time(1) = MY_MPI_WTIME()
      do ja = i_start_at, i_final_at
         ja_atom = ja 

         if (lmask) then
            if (.not.config_desc(iconf)%amask(ja)) cycle
         end if

         ! if (rangml==1)  write (6,*) config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ja)), config_real(iconf)%Z_per_type(config_real(iconf)%itype(ja))
         if (linvisible .and. weighted) then
            if (config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ja))) cycle
         end if

         if (weighted) then
            factor_ja = config_real(iconf)%weight_per_type(config_real(iconf)%itype(ja))
            if (weighted_3ch) factor_ja_3ch = config_real(iconf)%weight_per_type_3ch(config_real(iconf)%itype(ja))
         else
            factor_ja = 1.d0
         end if

!---------------neighbours----------------------
         if (Nfix) then 
              call build_local_neighbours_ja(iconf, ja, imm, xpnp, r_cut, d_n_neigh, d_kind_neigh, &
                                  r_central, i_type, i_type_db, i_central, tmp_dxp, tmp_xp, &
                                  max_neigh_local, iw2)            
             call build_local_neighbours_ja_Nfix(r_central, fix_Nmax_neigh, max_neigh_local, out_max_neigh_local,  i_central_Nfix, rcut_fix)
             !if (max_neigh_local > out_max_neigh_local) max_neigh_local = out_max_neigh_local 
          else 
             call build_local_neighbours_ja(iconf, ja, imm, xpnp, r_cut, d_n_neigh, d_kind_neigh, &
                                  r_central, i_type, i_type_db, i_central, tmp_dxp, tmp_xp, &
                                  max_neigh_local, iw2)   
               
         end if
         !d write(*,*) 'mmmmmax_neigh_local', max_neigh_local, out_max_neigh_local,  rcut_fix
         if (max_neigh_local == 0) cycle
         if (max_neigh_local > max_neigh) then
           call log_critical("ML: Fatal error in "//NAMECURRENT//" concerning the the number of neighbours. ")
           call log_critical("The number of max_neig_local is bigger than max_neigh "// vtoa(max_neigh_local) //" "// vtoa(max_neigh))
           call log_critical("possible solutions: decrease r_cut or increase max_neigh") 
           stop 'decrease rcut for that descriptor'
         end if
         call reallocate_neigh_ja (max_neigh_local, ja_atom, i_central, i_type, i_type_db, r_central, tmp_dxp)
         call build_neigh_ja_type (max_neigh_local, ja_atom, i_central, i_type, i_type_db, r_central, tmp_dxp)
!---------------neighbours----------------------



         ! cmm(:,:,:) = 0.d0
         ! if (weighted)     cmm2(:,:,:)=0.d0
         ! if (weighted_3ch) cmm3(:,:,:)=0.d0

         ! INIT cmm, cmm2, cmm3
         ! cmm(:,:,:) = 0.d0
         do j = 0, jj_max; do m1 = -j, j, 2; do m2 = -j, j, 2
                  cmm(m2, m1, j) = 0.d0
               end do; end do; end do
         if (weighted) then
            ! cmm2(:,:,:)=0.d0
            do j = 0, jj_max; do m1 = -j, j, 2; do m2 = -j, j, 2
                     cmm2(m2, m1, j) = 0.d0
                  end do; end do; end do

            if (weighted_3ch) then
               ! cmm3(:,:,:)=0.d0
               do j = 0, jj_max; do m1 = -j, j, 2; do m2 = -j, j, 2
                        cmm3(m2, m1, j) = 0.d0
                     end do; end do; end do
            end if
         end if

         !    do j = 0,jj_max
         !      do m1=-j,j,2
         !        cmm(m1,m1,j) = 1.d0
         !        if (weighted)  then
         !          cmm2(m1,m1,j) = factor_ja
         !          if (weighted_3ch) cmm3(m1,m1,j) = factor_ja_3ch
         !        end if
         !      end do
         !    end do


         ! this way of init gives 1/3 of BSO4 execution time ... WTF !!!
         !  if (desc_forces_local) then
         !    dcmm(:,:,:,:,:) = 0.d0
         !    if (weighted) then
         !      dcmm2(:,:,:,:,:) = 0.d0
         !      if (weighted_3ch)  dcmm3(:,:,:,:,:) = 0.d0
         !    end if
         !  end if
         ! dcmm(:,:,:,:,:) = 0.d0
         ! THIS initializatio reduce to reduce to 1/9 from BSO4 execution time
         !DCOS do j=0,jj_max
         !DCOS do m1 = -j,j,2
         !DCOS do m2 = -j,j,2
         !DCOS do ia_n = 0,imm_neigh
         !DCOS do ix=1,3
         !DCOS dcmm(m1,m2,j,ia_n,ix) = 0.d0
         !DCOS end do
         !DCOS end do
         !DCOS end do
         !DCOS end do
         !DCOS end do
         !What is crazy is that I do not need any initialization
         !DCOS if (weighted) then
         !DCOS   dcmm2(:,:,:,:,:) = 0.d0
         !DCOS     if (weighted_3ch)  dcmm3(:,:,:,:,:) = 0.d0
         !DCOS  end if
         ia_n = 0



         if (debug_time) time(7) = MY_MPI_WTIME()
         do iw = 1, max_neigh_local
            ia = i_central(iw)
            !if (ja == ia) cycle
            
            if (weighted) then
               if (linvisible) then
                  if (config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ia))) cycle
               end if
               factor_ia = config_real(iconf)%weight_per_type(config_real(iconf)%itype(ia))
               if (weighted_3ch) factor_ia_3ch = config_real(iconf)%weight_per_type_3ch(config_real(iconf)%itype(ia))
            else
               factor_ia = 1.d0
            end if
!#ifdef MLD_NDM
!!!$        if (small) then
!#else
!            if (small) then
!#endif
!               r_ji = config_real(iconf)%r_ij(ja, iw)
!               r2_ji = r_ji**2
!               dxp_ji(:) = config_real(iconf)%u_ij(ja, iw, :)
!#ifdef MLD_NDM
!!!$        else
!!!$          dxp_ji(1:3) = xpnp(1:3, ia) - xpnp(1:3, ja)
!!!$          ds = MatMul(dxp_ji, bg)
!!!$          WHERE ((ds .GT. 0.5d0) .OR. (ds .LT. -0.5d0))
!!!$            ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
!!!$          END WHERE
!!!$          dxp_ji = MatMul(at, ds)/A2cm
!!!$          r2_ji = Sum(dxp_ji(1:3)**2)
!!!$          r_ji = dsqrt(r2_ji)
!!!$        end if
!#else
!            else
!               dxp_ji(1:3) = xpnp(1:3, ia) - xpnp(1:3, ja)
!               ds = MatMul(dxp_ji, bg)
!               WHERE ((ds .GT. 0.5d0) .OR. (ds .LT. -0.5d0))
!                  ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
!               END WHERE
!               dxp_ji = MatMul(at, ds)/A2cm
!               r2_ji = Sum(dxp_ji(1:3)**2)
!               r_ji = dsqrt(r2_ji)
!            end if
!#endif

            r_ji = r_central(iw) 
            dxp_ji(:) = tmp_dxp(:,iw)
            r2_ji = r_ji**2 

            if (Nfix) then 
              if (r_ji .gt. rcut_fix) cycle
            else 
              if (r_ji .gt. r_cut) cycle
            end if 
            ia_n = ia_n + 1

            ! fcut = 0.5d0*(cos(one_pi*r_ji/r_cut) + 1.d0)
            !dfcut = -0.5d0*one_pi/r_cut*sin(one_pi*r_ji/r_cut)

            !TODOrcut
            !call fcut_rij(type_fcut, r_ji, r_cut, r_cut_width, desc_forces_local, fcut_out, dfcut_out)
            if (Nfix) then 
               call fcut_rij(2, r_ji, rcut_fix, r_cut_width, desc_forces_local, fcut_out, dfcut_out)
            else 
               call fcut_rij(2, r_ji, r_cut, r_cut_width, desc_forces_local, fcut_out, dfcut_out)
            end if


            if (Nfix) then 
               ! For fixed r_cut from n-neighbour fixed ... no inner cut.  
               fcut_in = 0.d0
               dfcut_in = 0.d0
            else 
               ! Pair-specific inner cutoff
               local_r_cut_in = r_cut_pair_in(i_type_db(0), i_type_db(iw))
               if (local_r_cut_in > 0 ) then
                  !TODOrcut
                  !call fcut_rij(type_fcut, r_ji, r_cut_in + r_cut_width_in, r_cut_width_in, desc_forces_local, fcut_in, dfcut_in)
                  !if (Nfix) then
                  !   call fcut_rij(3, r_ji, rcut_fix + r_cut_width_in, r_cut_width_in, desc_forces_local, fcut_in, dfcut_in)
                  !else 
                     call fcut_rij(3, r_ji, local_r_cut_in + r_cut_width_in, r_cut_width_in, desc_forces_local, fcut_in, dfcut_in)
                  !end if 
               else
                  fcut_in = 0.d0
                  dfcut_in = 0.d0
               end if
            end if 
            fcut = (1.d0 - fcut_in) * fcut_out
            if (desc_forces_local) dfcut =-dfcut_in * fcut_out +  (1.d0 - fcut_in) * dfcut_out


            if (weighted) then
               fcut_w = fcut*factor_ia
               dfcut_w = dfcut*factor_ia
               if (weighted_3ch) then
                  fcut_w_3ch = fcut*factor_ia_3ch
                  dfcut_w_3ch = dfcut*factor_ia_3ch
               end if
            end if

            if (debug_time) time(5) = MY_MPI_WTIME()
            call spherical_4d(ia_n, dxp_ji, r_ji, Umm, dUmm)
            if (debug_time) then
               time(6) = MY_MPI_WTIME()
               tot_time(2) = tot_time(2) + time(6) - time(5)
            end if

            !cmm2(:,:,:)= cmm2(:,:,:)+Umm(:,:,:)*fcut_w
            do j = 0, jj_max; do m1 = -j, j, 2; do m2 = -j, j, 2
                     cmm(m2, m1, j) = cmm(m2, m1, j) + Umm(m2, m1, j)*fcut
                  end do; end do; end do
            if (weighted) then
               !cmm2(:,:,:)= cmm2(:,:,:)+Umm(:,:,:)*fcut_w
               do j = 0, jj_max; do m1 = -j, j, 2; do m2 = -j, j, 2
                        cmm2(m2, m1, j) = cmm2(m2, m1, j) + Umm(m2, m1, j)*fcut_w
                     end do; end do; end do
               if (weighted_3ch) then
                  !cmm3(:,:,:)= cmm3(:,:,:)+Umm(:,:,:)*fcut_w_3ch
                  do j = 0, jj_max; do m1 = -j, j, 2; do m2 = -j, j, 2
                           cmm3(m2, m1, j) = cmm3(m2, m1, j) + Umm(m2, m1, j)*fcut_w_3ch
                        end do; end do; end do
               end if
            end if
            !---------------!
            ! ONCE WAS HERE !
            !---------------!
            !        if (desc_forces_local) then
            !          do ix=1,3
            !              dcmm(:,:,:,ia_n,ix) =  (dfcut*dxp_ji(ix)*Umm(:,:,:)/r_ji + fcut*dUmm(:,:,:,ix))
            !            if (weighted) then
            !              dcmm2(:,:,:,ia_n,ix) =  (dfcut_w    *dxp_ji(ix)*Umm(:,:,:)/r_ji + fcut_w    *dUmm(:,:,:,ix))
            !              if (weighted_3ch) dcmm3(:,:,:,ia_n,ix) =  (dfcut_w_3ch*dxp_ji(ix)*Umm(:,:,:)/r_ji + fcut_w_3ch*dUmm(:,:,:,ix))
            !            end if
            !          end do
            !        end if
            !---------------!
            ! ONCE WAS HERE !
            !---------------!
            d_kind_neigh(ja, ia_n) = ia
         end do                  ! iw


         do j = 0, jj_max
            do m1 = -j, j, 2
               cmm(m1, m1, j) = cmm(m1, m1, j) + 1.d0
               if (weighted) then
                  cmm2(m1, m1, j) = cmm2(m1, m1, j) + 1.d0         ! factor_ja
                  if (weighted_3ch) cmm3(m1, m1, j) = cmm3(m1, m1, j) + 1.d0                ! factor_ja_3ch
               end if
            end do
         end do



         if (debug_time) then
            time(8) = MY_MPI_WTIME()
            tot_time(4) = tot_time(4) + time(8) - time(7)
         end if

         if (debug_time) ta = MY_MPI_WTIME()
         !  end do
         call pack_cmm_into_ZcmmABC(jj_max)
         !call pack_Zcmm_into_ZAcmm_ZBcmm_ZCcmm(jj_max)
         if (debug_time) then
            tb = MY_MPI_WTIME()
            tot_time(8) = tot_time(8) + tb - ta
         end if



         if (debug_time) time(3) = MY_MPI_WTIME()
         d_n_neigh(ja) = ia_n
         !DCOS    if (desc_forces_local) then
         !DCOS      dcmm(:,:,:,0,:) = - sum(dcmm(:,:,:,:,:), dim=4)
         !DCOS      if (weighted) then
         !DCOS        dcmm2(:,:,:,0,:) = - sum(dcmm2(:,:,:,:,:), dim=4)
         !DCOS        if (weighted_3ch) dcmm3(:,:,:,0,:) = - sum(dcmm3(:,:,:,:,:), dim=4)
         !DCOS      end if
         !DCOS    end if

         i_bi = 0
         n_axpy = 3*ia_n
         it_cg = 0
         do l1 = 0, jj_max

            if (lbso4_diag) then
               !GaborLike l2=l1
               l1_min = l1
               l1_max = l1
            else
               !Tlike do l2=0,l1
               l1_min = 0
               l1_max = l1
            end if

            do l2 = l1_min, l1_max
               do l = abs(l1 - l2), min(jj_max, l1 + l2)

                  if (mod(l1 + l2 + l, 2) == 1) cycle
                  if (.not. (lbso4_diag)) then
                     if (l < l1) cycle       ! this comes from Aidan Thompson and SNAP
                  end if
                  i_bi = i_bi + 1
                  etmp = 0.d0
                  if (desc_forces_local) dsum_bi_ia(:) = czero
                  if (weighted) then
                     etmp2 = 0.d0
                     if (desc_forces_local) dsum_bi_ia_w(:) = czero
                     if (weighted_3ch) then
                        etmp3 = 0.d0
                        if (desc_forces_local) dsum_bi_ia_w_3ch(:) = czero
                     end if
                  end if

                  if (debug_time) time(9) = MY_MPI_WTIME()

                  do m1 = -l, l, 2
                     do m2 = -l, l, 2
                        !energy
                        sum_bi = ZAcmm(m1, m2, i_bi)
                        etmp = etmp + real(conjg(cmm(m1, m2, l))*sum_bi, kind=kind(0.d0))
                        if (weighted) then
                           sum_bi_w = ZAcmm2(m1, m2, i_bi)
                           etmp2 = etmp2 + real(conjg(cmm2(m1, m2, l))*sum_bi_w, kind=kind(0.d0))
                           if (weighted_3ch) then
                              sum_bi_w_3ch = ZAcmm3(m1, m2, i_bi)
                              etmp3 = etmp3 + real(conjg(cmm3(m1, m2, l))*sum_bi_w_3ch, kind=kind(0.d0))
                           end if
                        end if



                        if (desc_forces_local) then
                           !read_vecall => class_mml(m1,m2,l)%vecall
                           !call zaxpy(n_axpy,-sum_bi, read_vecall, incx, dsum_bi_ia,incy)
                           call zaxpy(n_axpy, -sum_bi, class_mml(m1, m2, l)%vecall, incx, dsum_bi_ia, incy)
                           !nullify(read_vecall)
                           !do ia=1,ia_n
                           !do ix=1,3
                           !dsum_bi_ia_all(1:3*ia_n) = dsum_bi_ia_all(1:3*ia_n)   - read_vecall(1:3*ia_n)*sum_bi
                           !end do
                           !end do
                           !forces
                           if (weighted) then
                              !read_vecall => class_mml(m1,m2,l)%vecall_w
                              !call zaxpy(n_axpy,-sum_bi_w, read_vecall, incx, dsum_bi_ia_w,incy)
                              call zaxpy(n_axpy, -sum_bi_w, class_mml(m1, m2, l)%vecall_w, incx, dsum_bi_ia_w, incy)
                              !nullify(read_vecall)
                              !do ix=1,3
                              !do ia=1,ia_n
                              !   dsum_bi_ia_w(ia,ix) = dsum_bi_ia_w(ia,ix)   - conjg(dcmm2(m1,m2,l,ia,ix))*sum_bi_w
                              !end do
                              !end do
                              if (weighted_3ch) then
                                 !read_vecall => class_mml(m1,m2,l)%vecall_w_3ch
                                 !call zaxpy(n_axpy,-sum_bi_w_3ch, read_vecall, incx, dsum_bi_ia_w_3ch,incy)
                                 call zaxpy(n_axpy, -sum_bi_w_3ch, class_mml(m1, m2, l)%vecall_w_3ch, incx, dsum_bi_ia_w_3ch, incy)
                                 !nullify(read_vecall)
                                 !do ix=1,3
                                 !do ia=1,ia_n
                                 !  !TOC dsum_bi_ia_w_3ch(ia,ix) = dsum_bi_ia_w_3ch(ia,ix)   - conjg(dcmm3(m1,m2,l,ia,ix))*sum_bi_w_3ch
                                 !end do
                                 !end do
                              end if
                           end if
                        end if

                     end do                  ! m2
                  end do                  ! m1


                  tmpll = dble(1 + l)/dble(1 + l1)
                  do m1 = -l1, l1, 2
                     do m2 = -l1, l1, 2
                        sum_bi = ZBcmm(m1, m2, i_bi)*tmpll
                        if (weighted) then
                           sum_bi_w = ZBcmm2(m1, m2, i_bi)*tmpll
                           if (weighted_3ch) sum_bi_w_3ch = ZBcmm3(m1, m2, i_bi)*tmpll
                        end if

                        if (desc_forces_local) then
                           !read_vecall => class_mml(m1,m2,l1)%vecall
                           !call zaxpy(n_axpy,-sum_bi, read_vecall, incx, dsum_bi_ia,incy)
                           call zaxpy(n_axpy, -sum_bi, class_mml(m1, m2, l1)%vecall, incx, dsum_bi_ia, incy)
                           !nullify(read_vecall)
                           !do ia=1,ia_n
                           !do ix=1,3
                           !dsum_bi_ia_all(1:3*ia_n) = dsum_bi_ia_all(1:3*ia_n)   - read_vecall(1:3*ia_n)*sum_bi
                           !end do
                           !end do
                           if (weighted) then
                              !read_vecall => class_mml(m1,m2,l1)%vecall_w
                              !call zaxpy(n_axpy,-sum_bi_w, read_vecall, incx, dsum_bi_ia_w,incy)
                              call zaxpy(n_axpy, -sum_bi_w, class_mml(m1, m2, l1)%vecall_w, incx, dsum_bi_ia_w, incy)
                              !nullify(read_vecall)
                              !do ix=1,3
                              !do ia=1,ia_n
                              !  !TOC dsum_bi_ia_w(ia,ix) = dsum_bi_ia_w(ia,ix)   - conjg(dcmm2(m1,m2,l1,ia,ix))*sum_bi_w
                              !end do
                              !end do

                              if (weighted_3ch) then
                                 !read_vecall => class_mml(m1,m2,l1)%vecall_w_3ch
                                 !call zaxpy(n_axpy,-sum_bi_w_3ch, read_vecall, incx, dsum_bi_ia_w_3ch,incy)
                                 call zaxpy(n_axpy, -sum_bi_w_3ch, class_mml(m1, m2, l1)%vecall_w_3ch, incx, dsum_bi_ia_w_3ch, incy)
                                 !nullify(read_vecall)
                                 !do ix=1,3
                                 !do ia=1,ia_n
                                 !  !TOC dsum_bi_ia_w_3ch(ia,ix) = dsum_bi_ia_w_3ch(ia,ix)   - conjg(dcmm3(m1,m2,l1,ia,ix))*sum_bi_w_3ch
                                 !end do
                                 !end do
                              end if
                           end if
                        end if

                     end do                  ! m2
                  end do                  ! m1

                  tmpll = dble(1 + l)/dble(1 + l2)
                  do m1 = -l2, l2, 2
                     do m2 = -l2, l2, 2
                        !sum_bi = ZCcmm(m1,m2,i_bi)*dble(1+l)/dble(1+l2)
                        sum_bi = ZCcmm(m1, m2, i_bi)*tmpll
                        if (weighted) then
                           sum_bi_w = ZCcmm2(m1, m2, i_bi)*tmpll
                           if (weighted_3ch) sum_bi_w_3ch = ZCcmm3(m1, m2, i_bi)*tmpll
                        end if

                        if (desc_forces_local) then
                           !read_vecall => class_mml(m1,m2,l2)%vecall
                           !call zaxpy(n_axpy,-sum_bi, read_vecall, incx, dsum_bi_ia,incy)
                           call zaxpy(n_axpy, -sum_bi, class_mml(m1, m2, l2)%vecall, incx, dsum_bi_ia, incy)
                           !nullify(read_vecall)
                           !do ia=1,ia_n
                           !do ix=1,3
                           !dsum_bi_ia_all(1:3*ia_n) = dsum_bi_ia_all(1:3*ia_n)   - read_vecall(1:3*ia_n)*sum_bi
                           !end do
                           !end do
                           if (weighted) then
                              !read_vecall => class_mml(m1,m2,l2)%vecall_w
                              !call zaxpy(n_axpy,-sum_bi_w, read_vecall, incx, dsum_bi_ia_w,incy)
                              call zaxpy(n_axpy, -sum_bi_w, class_mml(m1, m2, l2)%vecall_w, incx, dsum_bi_ia_w, incy)
                              !nullify(read_vecall)
                              !do ix=1,3
                              !do ia=1,ia_n
                              !  !TOC dsum_bi_ia_w(ia,ix) = dsum_bi_ia_w(ia,ix)   - conjg(dcmm2(m1,m2,l2,ia,ix))*sum_bi_w
                              !end do
                              !end do
                              if (weighted_3ch) then
                                 !read_vecall => class_mml(m1,m2,l2)%vecall_w_3ch
                                 !call zaxpy(n_axpy,-sum_bi_w_3ch, read_vecall, incx, dsum_bi_ia_w_3ch,incy)
                                 call zaxpy(n_axpy, -sum_bi_w_3ch, class_mml(m1, m2, l2)%vecall_w_3ch, incx, dsum_bi_ia_w_3ch, incy)
                                 !nullify(read_vecall)
                                 !do ix=1,3
                                 !do ia=1,ia_n
                                 !  !TOC dsum_bi_ia_w_3ch(ia,ix) = dsum_bi_ia_w_3ch(ia,ix)   - conjg(dcmm3(m1,m2,l2,ia,ix))*sum_bi_w_3ch
                                 !end do
                                 !end do
                              end if
                           end if
                        end if

                     end do                  ! m2
                  end do                  ! m1

                  if (debug_time) then
                     time(10) = MY_MPI_WTIME()
                     tot_time(6) = tot_time(6) + time(10) - time(9)
                  end if
                  config_desc(iconf)%energy(i_bi, ja) = etmp
                  if (weighted) then
                     i_bi2 = i_bi + bisso4_dim
                     i_bi3 = i_bi2 + bisso4_dim
                     config_desc(iconf)%energy(i_bi2, ja) = etmp2
                     if (weighted_3ch) config_desc(iconf)%energy(i_bi3, ja) = etmp3
                  end if

                  if (desc_forces_local) then
                     do ia = 1, ia_n
                        config_desc(iconf)%force(i_bi, ja, ia, 1) = real(dsum_bi_ia(3*ia - 2), kind=kind(0.d0))
                        config_desc(iconf)%force(i_bi, ja, ia, 2) = real(dsum_bi_ia(3*ia - 1), kind=kind(0.d0))
                        config_desc(iconf)%force(i_bi, ja, ia, 3) = real(dsum_bi_ia(3*ia), kind=kind(0.d0))
                     end do
                     if (weighted) then
                        do ia = 1, ia_n
                           config_desc(iconf)%force(i_bi2, ja, ia, 1) = real(dsum_bi_ia_w(3*ia - 2), kind=kind(0.d0))
                           config_desc(iconf)%force(i_bi2, ja, ia, 2) = real(dsum_bi_ia_w(3*ia - 1), kind=kind(0.d0))
                           config_desc(iconf)%force(i_bi2, ja, ia, 3) = real(dsum_bi_ia_w(3*ia), kind=kind(0.d0))
                        end do
                        if (weighted_3ch) then
                           do ia = 1, ia_n
                              config_desc(iconf)%force(i_bi3, ja, ia, 1) = real(dsum_bi_ia_w_3ch(3*ia - 2), kind=kind(0.d0))
                              config_desc(iconf)%force(i_bi3, ja, ia, 2) = real(dsum_bi_ia_w_3ch(3*ia - 1), kind=kind(0.d0))
                              config_desc(iconf)%force(i_bi3, ja, ia, 3) = real(dsum_bi_ia_w_3ch(3*ia), kind=kind(0.d0))
                           end do
                        end if
                     end if
                  end if

               end do                  ! l
            end do                  ! l2
         end do                  ! l1

         do ia = 1, ia_n
            if (desc_forces_local) config_desc(iconf)%force(:, ja, 0, :) = config_desc(iconf)%force(:, ja, 0, :) - config_desc(iconf)%force(:, ja, ia, :)
            !if (desc_forces_local) config_desc(iconf)%force(:,ja, 0,:) =  -  sum(config_desc(iconf)%force(:,ja, 1:ia_n,:), dim=2)
         end do                  ! ia from ia_m
         if (debug_time) then
            time(4) = MY_MPI_WTIME()
            tot_time(5) = tot_time(5) + time(4) - time(3)
         end if


         ! NON NORMALIZED VERSION
         if (ia_n == 0) then
            config_desc(iconf)%energy(:, ja) = 0.d0
         else
            config_desc(iconf)%energy(:, ja) = config_desc(iconf)%energy(:, ja)       ! /dble(ia_n)
         end if

         if (weighted) then
            config_desc(iconf)%energy(bisso4_dim + 1:2*bisso4_dim, ja) = config_desc(iconf)%energy(bisso4_dim + 1:2*bisso4_dim, ja)     ! *factor_ja
            if (weighted_3ch) config_desc(iconf)%energy(2*bisso4_dim + 1:3*bisso4_dim, ja) = config_desc(iconf)%energy(2*bisso4_dim + 1:3*bisso4_dim, ja)        ! *factor_ja_3ch
            if (desc_forces_local) then
               config_desc(iconf)%force(bisso4_dim + 1:2*bisso4_dim, ja, :, :) = config_desc(iconf)%force(bisso4_dim + 1:2*bisso4_dim, ja, :, :)                    ! *factor_ja
               if (weighted_3ch) config_desc(iconf)%force(2*bisso4_dim + 1:3*bisso4_dim, ja, :, :) = config_desc(iconf)%force(2*bisso4_dim + 1:3*bisso4_dim, ja, :, :)                       ! *factor_ja_3ch
            end if
         end if
      end do                  ! ja
      deallocate (xpnp)

      if (debug_time) then
         time(2) = MY_MPI_WTIME()
         tot_time(1) = tot_time(1) + time(2) - time(1)
      end if


      !if (debug_time) then
      !  call repport_time(2, 0.d0, tot_time(1), "ML: full")
      !  call repport_time(4, 0.d0, tot_time(2), "ML: Ujmm ")
      !  call repport_time(4, 0.d0, tot_time(4), "ML: Ujmm + nn")
      !  call repport_time(4, 0.d0, tot_time(8), "ML: time pack")
      !  call repport_time(4, 0.d0, tot_time(5), "ML: CG last")
      !  call repport_time(4, 0.d0, tot_time(6), "ML: inner CG ")
      !end if
      _MLD_END_

   end subroutine compute_bispectrum_so4



   subroutine gen_dimension_for_bispectrum_so4()

      use module_bispectrum_so4, only: lbso4_diag, bisso4_cg_dim,  &
         bisso4_cg_A_dim, bisso4_cg_B_dim, bisso4_cg_C_dim, &
         cg_A, cg_B, cg_C, &
         cmm, cmm2, cmm3, &
         ZAcmm, ZAcmm2, ZAcmm3, &
         ZBcmm, ZBcmm2, ZBcmm3, &
         ZCcmm, ZCcmm2, ZCcmm3, &
         class_mml, &
         bisso4_l, bisso4_l1, bisso4_l2, bisso4_dim

      use ml_in_ndm_module, ONLY: imm_neigh, weighted, weighted_3ch, jj_max,  cg_vector
      
      integer  :: l1, l2, l, j, l1_min, l1_max, m1, m2, ma, mb, m1ma, m2mb
      integer  :: i_bi, it_cg, it_cgA, it_cgB, it_cgC, it_umm


      ! pre-compute dimension of BSO4
      i_bi = 0
      it_cg = 0
      do l1 = 0, jj_max
         !l2=l1 ! following Gabor  only the diagonal elements are important
         if (lbso4_diag) then
            !GaborLike l2=l1
            l1_min = l1
            l1_max = l1
         else
            !Tlike do l2=0,l1
            l1_min = 0
            l1_max = l1
         end if
         do l2 = l1_min, l1_max
            !do l2=0, l1  ! in the end we want only the componenets with l1 <= l2 <= l
            do l = abs(l1 - l2), min(jj_max, l1 + l2)
               if (mod(l1 + l2 + l, 2) == 1) cycle              ! assure invariance par reflexion
               if (.not. (lbso4_diag)) then
                  if (l < l1) cycle       ! this comes from Thompson
               end if
               i_bi = i_bi + 1
               do m1 = -l, l, 2
                  do m2 = -l, l, 2
                     do ma = max(-l1, m1 - l2), min(l1, m1 + l2), 2
                        do mb = max(-l1, m2 - l2), min(l1, m2 + l2), 2
                           it_cg = it_cg + 1
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      bisso4_cg_dim = it_cg
      bisso4_dim = i_bi
      if (allocated(bisso4_l)) deallocate (bisso4_l); allocate (bisso4_l(bisso4_dim))
      if (allocated(bisso4_l1)) deallocate (bisso4_l1); allocate (bisso4_l1(bisso4_dim))
      if (allocated(bisso4_l2)) deallocate (bisso4_l2); allocate (bisso4_l2(bisso4_dim))



      ! pre-compute cgA,cgB,cgC in order to compute
      ! i) the dimensions 2) computation
      ! used for packing ZAcmm ZBcmm and ZCcmm

      i_bi = 0
      it_cgA = 0
      it_cgB = 0
      it_cgC = 0
      it_umm = 0
      do l1 = 0, jj_max
         !l2=l1 ! following Gabor  only the diagonal elements are important
         if (lbso4_diag) then
            !GaborLike l2=l1
            l1_min = l1
            l1_max = l1
         else
            !Tlike do l2=0,l1
            l1_min = 0
            l1_max = l1
         end if
         do l2 = l1_min, l1_max
            !do l2=0, l1
            do l = abs(l1 - l2), min(jj_max, l1 + l2)
               if (mod(l1 + l2 + l, 2) == 1) cycle              ! assure invariance par reflexion
               if (.not. (lbso4_diag)) then
                  if (l < l1) cycle       ! in the end we want only the componenets with l1 <= l2 <= l
               end if
               i_bi = i_bi + 1
               bisso4_l(i_bi) = l
               bisso4_l1(i_bi) = l1
               bisso4_l2(i_bi) = l2

               do m1 = -l, l, 2
                  do m2 = -l, l, 2
                     it_umm = it_umm + 1
                     do ma = max(-l1, m1 - l2), min(l1, m1 + l2), 2
                        do mb = max(-l1, m2 - l2), min(l1, m2 + l2), 2
                           it_cgA = it_cgA + 1
                        end do                  ! ma
                     end do                  ! mb
                  end do                  ! m1
               end do                  ! m2


               do m1 = -l1, l1, 2
                  do m2 = -l1, l1, 2
                     it_umm = it_umm + 1
                     do ma = max(-l, m1 - l2), min(l, m1 + l2), 2
                        do mb = max(-l, m2 - l2), min(l, m2 + l2), 2
                           it_cgB = it_cgB + 1
                        end do                  ! ma
                     end do                  ! mb
                  end do                  ! m1
               end do                  ! m2


               do m1 = -l2, l2, 2
                  do m2 = -l2, l2, 2
                     it_umm = it_umm + 1
                     do ma = max(-l1, m1 - l), min(l1, m1 + l), 2
                        do mb = max(-l1, m2 - l), min(l1, m2 + l), 2
                           it_cgC = it_cgC + 1
                        end do                  ! ma
                     end do                  ! mb
                  end do                  ! m1
               end do                  ! m2


            end do
         end do
      end do


      bisso4_cg_A_dim = it_cgA
      bisso4_cg_B_dim = it_cgB
      bisso4_cg_C_dim = it_cgC

      if (allocated(cg_A)) deallocate (cg_A); allocate (cg_A(bisso4_cg_A_dim))
      if (allocated(cg_B)) deallocate (cg_B); allocate (cg_B(bisso4_cg_B_dim))
      if (allocated(cg_C)) deallocate (cg_C); allocate (cg_C(bisso4_cg_C_dim))
      cg_A(:) = 0.d0
      cg_B(:) = 0.d0
      cg_C(:) = 0.d0

      i_bi = 0
      it_cgA = 0
      it_cgB = 0
      it_cgC = 0
      do l1 = 0, jj_max
         !l2=l1 ! following Gabor  only the diagonal elements are important
         if (lbso4_diag) then
            !GaborLike l2=l1
            l1_min = l1
            l1_max = l1
         else
            !Tlike do l2=0,l1
            l1_min = 0
            l1_max = l1
         end if
         do l2 = l1_min, l1_max
            do l = abs(l1 - l2), min(jj_max, l1 + l2)
               if (mod(l1 + l2 + l, 2) == 1) cycle              ! assure invariance par reflexion
               if (.not. (lbso4_diag)) then
                  if (l < l1) cycle       ! in the end we want only the componenets with l1 <= l2 <= l
               end if
               i_bi = i_bi + 1

               do m1 = -l, l, 2
                  do m2 = -l, l, 2
                     do ma = max(-l1, m1 - l2), min(l1, m1 + l2), 2
                        do mb = max(-l1, m2 - l2), min(l1, m2 + l2), 2
                           m1ma = m1 - ma
                           m2mb = m2 - mb
                           it_cgA = it_cgA + 1
                           cg_A(it_cgA) = cg_vector(l1, ma, l2, m1ma, l, m1)*cg_vector(l1, mb, l2, m2mb, l, m2)
                        end do                  ! ma
                     end do                  ! mb
                  end do                  ! m1
               end do                  ! m2


               do m1 = -l1, l1, 2
                  do m2 = -l1, l1, 2
                     do ma = max(-l, m1 - l2), min(l, m1 + l2), 2
                        do mb = max(-l, m2 - l2), min(l, m2 + l2), 2
                           m1ma = m1 - ma
                           m2mb = m2 - mb
                           it_cgB = it_cgB + 1
                           cg_B(it_cgB) = cg_vector(l, ma, l2, m1ma, l1, m1)*cg_vector(l, mb, l2, m2mb, l1, m2)
                        end do                  ! ma
                     end do                  ! mb
                  end do                  ! m1
               end do                  ! m2


               do m1 = -l2, l2, 2
                  do m2 = -l2, l2, 2
                     do ma = max(-l1, m1 - l), min(l1, m1 + l), 2
                        do mb = max(-l1, m2 - l), min(l1, m2 + l), 2
                           m1ma = m1 - ma
                           m2mb = m2 - mb
                           it_cgC = it_cgC + 1
                           cg_C(it_cgC) = cg_vector(l1, ma, l, m1ma, l2, m1)*cg_vector(l1, mb, l, m2mb, l2, m2)
                        end do                  ! ma
                     end do                  ! mb
                  end do                  ! m1
               end do                  ! m2


            end do
         end do
      end do



      if (allocated(class_mml)) deallocate (class_mml)
      allocate (class_mml(-jj_max:jj_max, -jj_max:jj_max, 0:jj_max))
      do j = 0, jj_max; do m1 = -j, j, 2; do m2 = -j, j, 2
               if (allocated(class_mml(m2, m1, j)%vecall)) deallocate (class_mml(m2, m1, j)%vecall)
               allocate (class_mml(m2, m1, j)%vecall(1:3*imm_neigh))
            end do; end do; end do

      if (allocated(cmm)) deallocate (cmm); allocate (cmm(-jj_max:jj_max, -jj_max:jj_max, 0:jj_max))
      if (allocated(ZAcmm)) deallocate (ZAcmm); allocate (ZAcmm(-jj_max:jj_max, -jj_max:jj_max, bisso4_dim))
      if (allocated(ZBcmm)) deallocate (ZBcmm); allocate (ZBcmm(-jj_max:jj_max, -jj_max:jj_max, bisso4_dim))
      if (allocated(ZCcmm)) deallocate (ZCcmm); allocate (ZCcmm(-jj_max:jj_max, -jj_max:jj_max, bisso4_dim))

      if (weighted) then

         if (allocated(cmm2)) deallocate (cmm2); allocate (cmm2(-jj_max:jj_max, -jj_max:jj_max, 0:jj_max))
         do j = 0, jj_max; do m1 = -j, j, 2; do m2 = -j, j, 2
                  if (allocated(class_mml(m2, m1, j)%vecall_w)) deallocate (class_mml(m2, m1, j)%vecall_w); allocate (class_mml(m2, m1, j)%vecall_w(1:3*imm_neigh))
               end do; end do; end do
         if (allocated(ZAcmm2)) deallocate (ZAcmm2); allocate (ZAcmm2(-jj_max:jj_max, -jj_max:jj_max, bisso4_dim))
         if (allocated(ZBcmm2)) deallocate (ZBcmm2); allocate (ZBcmm2(-jj_max:jj_max, -jj_max:jj_max, bisso4_dim))
         if (allocated(ZCcmm2)) deallocate (ZCcmm2); allocate (ZCcmm2(-jj_max:jj_max, -jj_max:jj_max, bisso4_dim))

         if (weighted_3ch) then

            if (allocated(cmm3)) deallocate (cmm3); allocate (cmm3(-jj_max:jj_max, -jj_max:jj_max, 0:jj_max))
            do j = 0, jj_max; do m1 = -j, j, 2; do m2 = -j, j, 2
                     if (allocated(class_mml(m2, m1, j)%vecall_w_3ch)) deallocate (class_mml(m2, m1, j)%vecall_w_3ch); allocate (class_mml(m2, m1, j)%vecall_w_3ch(1:3*imm_neigh))
                  end do; end do; end do
            if (allocated(ZAcmm3)) deallocate (ZAcmm3); allocate (ZAcmm3(-jj_max:jj_max, -jj_max:jj_max, bisso4_dim))
            if (allocated(ZBcmm3)) deallocate (ZBcmm3); allocate (ZBcmm3(-jj_max:jj_max, -jj_max:jj_max, bisso4_dim))
            if (allocated(ZCcmm3)) deallocate (ZCcmm3); allocate (ZCcmm3(-jj_max:jj_max, -jj_max:jj_max, bisso4_dim))

         end if

      end if

   end subroutine gen_dimension_for_bispectrum_so4



   subroutine test_parameters_bso4()
      use ml_in_ndm_module, only: rangml, j_max
      use module_bispectrum_so4, only : inv_r0, inv_r0_input

      if (rangml == 0) then
         if (j_max < 0) then
            write (6, '("ML: bi_so4 parametrization j_max should be larger or at least equals 0")')
            stop 'read_ml_file bi_so4 j_max negative'
         end if
         if (.not. ((inv_r0_input > 0.d0) .and. (inv_r0_input < 1.d0))) then
            write (6, *) 'ML:', inv_r0, inv_r0_input
            write (6, '("ML: bi_so4 parametrization, inv_r0_input should be larger than 0 and lower than 1")')
            stop 'read_ml_file bi_so4 inv_r0_input beyound the limits'
         end if
      end if
   end subroutine test_parameters_bso4


   subroutine pack_cmm_into_ZcmmABC(jj_max)
      use ml_in_ndm_module, only:  weighted, weighted_3ch
      use module_bispectrum_so4, only: lbso4_diag, czero, cmm, cmm2, cmm3, &
         cg_a, cg_B, cg_C, &
         ZAcmm, ZBcmm, ZCcmm, &
         ZAcmm2, ZBcmm2, ZCcmm2, &
         ZAcmm3, ZBcmm3, ZCcmm3


      integer, intent(in)  :: jj_max

      integer  :: i_bi, l1, l2, l, l1_min, l1_max, m1ma, m2mb
      integer  :: m1, m2, ma, mb
      integer  :: itfA, itfB, itfC
      double complex :: sum_bi, sum_bi_w, sum_bi_w_3ch
      real(kind_double)   :: cg_local

      ZAcmm(:, :, :) = czero
      ZBcmm(:, :, :) = czero
      ZCcmm(:, :, :) = czero

      i_bi = 0
      itfA = 0
      itfB = 0
      itfC = 0
      do l1 = 0, jj_max

         if (lbso4_diag) then
            !GaborLike l2=l1
            l1_min = l1
            l1_max = l1
         else
            !Tlike do l2=0,l1
            l1_min = 0
            l1_max = l1
         end if

         do l2 = l1_min, l1_max
            do l = abs(l1 - l2), min(jj_max, l1 + l2)

               if (mod(l1 + l2 + l, 2) == 1) cycle
               if (.not. (lbso4_diag)) then
                  if (l < l1) cycle       ! this comes from Aidan Thompson and SNAP
               end if
               i_bi = i_bi + 1

               do m1 = -l, l, 2
                  do m2 = -l, l, 2
                     sum_bi = czero
                     sum_bi_w = czero
                     sum_bi_w_3ch = czero
                     do ma = max(-l1, m1 - l2), min(l1, m1 + l2), 2
                        do mb = max(-l1, m2 - l2), min(l1, m2 + l2), 2
                           itfA = itfA + 1
                           m1ma = m1 - ma
                           m2mb = m2 - mb
                           cg_local = cg_A(itfA)
                           sum_bi = sum_bi + cg_local*cmm(ma, mb, l1)*cmm(m1ma, m2mb, l2)
                           if (weighted) then
                              sum_bi_w = sum_bi_w + cg_local*cmm2(ma, mb, l1)*cmm2(m1ma, m2mb, l2)
                              if (weighted_3ch) sum_bi_w_3ch = sum_bi_w_3ch + cg_local*cmm3(ma, mb, l1)*cmm3(m1ma, m2mb, l2)
                           end if
                        end do                  ! mb
                     end do                  ! ma
                     ZAcmm(m1, m2, i_bi) = sum_bi
                     if (weighted) then
                        ZAcmm2(m1, m2, i_bi) = sum_bi_w
                        if (weighted_3ch) ZAcmm3(m1, m2, i_bi) = sum_bi_w_3ch
                     end if
                  end do                  ! m2
               end do                  ! m1

               do m1 = -l1, l1, 2
                  do m2 = -l1, l1, 2
                     sum_bi = czero
                     sum_bi_w = czero
                     sum_bi_w_3ch = czero
                     do ma = max(-l, m1 - l2), min(l, m1 + l2), 2
                        do mb = max(-l, m2 - l2), min(l, m2 + l2), 2
                           itfB = itfB + 1
                           m1ma = m1 - ma
                           m2mb = m2 - mb
                           cg_local = cg_B(itfB)
                           sum_bi = sum_bi + cg_local*cmm(ma, mb, l)*cmm(m1ma, m2mb, l2)
                           if (weighted) then
                              sum_bi_w = sum_bi_w + cg_local*cmm2(ma, mb, l)*cmm2(m1ma, m2mb, l2)
                              if (weighted_3ch) sum_bi_w_3ch = sum_bi_w_3ch + cg_local*cmm3(ma, mb, l)*cmm3(m1ma, m2mb, l2)
                           end if
                        end do                  ! mb
                     end do                  ! ma
                     ZBcmm(m1, m2, i_bi) = sum_bi
                     if (weighted) then
                        ZBcmm2(m1, m2, i_bi) = sum_bi_w
                        if (weighted_3ch) ZBcmm3(m1, m2, i_bi) = sum_bi_w_3ch
                     end if
                  end do                  ! m2
               end do                  ! m1

               do m1 = -l2, l2, 2
                  do m2 = -l2, l2, 2
                     sum_bi = czero
                     sum_bi_w = czero
                     sum_bi_w_3ch = czero
                     do ma = max(-l1, m1 - l), min(l1, m1 + l), 2
                        do mb = max(-l1, m2 - l), min(l1, m2 + l), 2
                           itfC = itfC + 1
                           m1ma = m1 - ma
                           m2mb = m2 - mb
                           cg_local = cg_C(itfC)
                           sum_bi = sum_bi + cg_local*cmm(ma, mb, l1)*cmm(m1ma, m2mb, l)
                           if (weighted) then
                              sum_bi_w = sum_bi_w + cg_local*cmm2(ma, mb, l1)*cmm2(m1ma, m2mb, l)
                              if (weighted_3ch) sum_bi_w_3ch = sum_bi_w_3ch + cg_local*cmm3(ma, mb, l1)*cmm3(m1ma, m2mb, l)
                           end if
                        end do                  ! mb
                     end do                  ! ma
                     ZCcmm(m1, m2, i_bi) = sum_bi
                     if (weighted) then
                        ZCcmm2(m1, m2, i_bi) = sum_bi_w
                        if (weighted_3ch) ZCcmm3(m1, m2, i_bi) = sum_bi_w_3ch
                     end if
                  end do                  ! m2
               end do                  ! m1


            end do                  ! l
         end do                  ! l2  Tlike enddo
      end do                  ! l1

   end subroutine pack_cmm_into_ZcmmABC

end module
