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


module module_neigh_local
   USE module_kind_variables, ONLY: kind_double

   real(kind_double) :: r_cut, r_cut_width
   real(kind_double) :: r_cut_in, r_cut_width_in
   ! pair-specific inner cutoff: r_cut_pair_in(i,j)
   ! dimensions: (fix_no_of_elements, fix_no_of_elements)
   real(kind_double), dimension(:,:), allocatable :: r_cut_pair_in
   ! ja atom
   integer  :: ja_atom, iw2
   ! the number of neighbours of ja_number
   integer  :: max_neigh_local
   integer  :: out_max_neigh_local 

   integer, dimension(:), allocatable     :: i_type, i_type_db, i_central
   ! r_central   - all the distance with respect central ja atom
   ! ur_central  - the same thing as previuos but with transformated distances
   ! r_fcut      - the fcut function

   real(kind_double), dimension(:), allocatable :: r_central, ur_central, r_fcut
   ! tmp_xp(3,:) = xp(3,ia_n)
   ! tmp_dxp(3,:) = xp(3,ia_n) - xp(3,ja)
   ! d_r_central(3,:)  - the derivatives
   ! d_ur_central(3,:) - the derivatives
   ! d_r_fcut(3,:) - the derivatives
   real(kind_double), dimension(:, :), allocatable    :: &
      d_r_central, d_ur_central, d_r_fcut, tmpcos_dxp, tmp_dxp, tmp_xp

   integer :: type_fcut

   type type_neigh_ja
     integer :: max_neigh_local
     integer :: ja_atom
     integer, dimension(:), allocatable               :: i_type, i_central, i_type_db 
     real(kind_double), dimension(:), allocatable     :: r_central    ! the distances with respect the central atom  
     real(kind_double), dimension(:,:), allocatable   :: tmp_dxp  ! the disance per x, y, z with respect the central atom 
   end type type_neigh_ja 

   type(type_neigh_ja) :: neigh_ja
contains

   subroutine build_neigh_ja_type(l_max_neigh_local, l_ja_atom, l_i_central, &
                            l_i_type, l_i_type_db, l_r_central, l_tmp_dxp)
     use module_kind_variables, only: kind_double
     implicit none
     integer, intent(in) :: l_max_neigh_local, l_ja_atom
     real(kind_double), dimension(3, 0:l_max_neigh_local), intent(in)  :: l_tmp_dxp
     real(kind_double), dimension(1:l_max_neigh_local), intent(in)  :: l_r_central
     integer,  dimension(0:l_max_neigh_local), intent(in)  :: l_i_central, l_i_type, l_i_type_db 


      neigh_ja%ja_atom = l_ja_atom
      neigh_ja%max_neigh_local = l_max_neigh_local
      
      if (allocated(neigh_ja%r_central)) deallocate(neigh_ja%r_central) ; allocate(neigh_ja%r_central(l_max_neigh_local))
      neigh_ja%r_central(:) = l_r_central(:)

      if (allocated(neigh_ja%tmp_dxp)) deallocate(neigh_ja%tmp_dxp) ; allocate(neigh_ja%tmp_dxp(3, 0:l_max_neigh_local))
      neigh_ja%tmp_dxp(:,:) = l_tmp_dxp(:,:)

      if (allocated(neigh_ja%i_central)) deallocate(neigh_ja%i_central) ; allocate(neigh_ja%i_central(0:l_max_neigh_local))
      neigh_ja%i_central(0:l_max_neigh_local) = l_i_central(0:l_max_neigh_local) 

      if (allocated(neigh_ja%i_type)) deallocate(neigh_ja%i_type) ; allocate(neigh_ja%i_type(0:l_max_neigh_local))
      neigh_ja%i_type(0:l_max_neigh_local) = l_i_type(0:l_max_neigh_local)

      if (allocated(neigh_ja%i_type_db)) deallocate(neigh_ja%i_type_db) ; allocate(neigh_ja%i_type_db(0:l_max_neigh_local))
      neigh_ja%i_type_db(0:l_max_neigh_local) = l_i_type_db(0:l_max_neigh_local)

   end subroutine build_neigh_ja_type   

   subroutine reallocate_neigh_ja (l_max_neigh_local, l_ja_atom, l_i_central, &
      l_i_type, l_i_type_db, l_r_central, l_tmp_dxp)
      use module_kind_variables, only: kind_double
      implicit none
      integer, intent(in) :: l_max_neigh_local, l_ja_atom
      real(kind_double), allocatable, dimension(:, :), intent(inout)  :: l_tmp_dxp
      real(kind_double), allocatable, dimension(:), intent(inout)  :: l_r_central
      integer, allocatable, dimension(:), intent(inout)  :: l_i_central, l_i_type, l_i_type_db 
  
      !local variables ...
      real(kind_double), allocatable, dimension(:) :: tmp_real
      real(kind_double), allocatable, dimension(:, :)    :: tmp_2real
      integer, allocatable, dimension(:)     :: tmp_integer

      !resize l_r_central
      if (allocated(tmp_real)) deallocate (tmp_real); allocate (tmp_real(l_max_neigh_local))
      tmp_real(1:l_max_neigh_local) = l_r_central(1:l_max_neigh_local)
      if (allocated(l_r_central)) deallocate (l_r_central); allocate (l_r_central(l_max_neigh_local))
      l_r_central(1:l_max_neigh_local) = tmp_real(1:l_max_neigh_local)
      deallocate (tmp_real)

      !resize dxp
      if (allocated(tmp_2real)) deallocate (tmp_2real); allocate (tmp_2real(3, 0:l_max_neigh_local))
      tmp_2real(1:3, 0:l_max_neigh_local) = l_tmp_dxp(1:3, 0:l_max_neigh_local)
      if (allocated(l_tmp_dxp)) deallocate (l_tmp_dxp); allocate (l_tmp_dxp(3, 0:l_max_neigh_local))
      l_tmp_dxp(1:3, 0:l_max_neigh_local) = tmp_2real(1:3, 0:l_max_neigh_local)
      deallocate (tmp_2real)

      !resize l_i_central
      if (allocated(tmp_integer)) deallocate (tmp_integer); allocate (tmp_integer(0:l_max_neigh_local))
      tmp_integer(1:l_max_neigh_local) = l_i_central(1:l_max_neigh_local)
      tmp_integer(0) = l_ja_atom
      if (allocated(l_i_central)) deallocate (l_i_central); allocate (l_i_central(0:l_max_neigh_local))
      l_i_central(0:l_max_neigh_local) = tmp_integer(0:l_max_neigh_local)
      deallocate (tmp_integer)

      !resize l_i_type
      if (allocated(tmp_integer)) deallocate (tmp_integer); allocate (tmp_integer(0:l_max_neigh_local))
      tmp_integer(1:l_max_neigh_local) = l_i_type(1:l_max_neigh_local)
      tmp_integer(0) = l_i_type(0)
      if (allocated(l_i_type)) deallocate (l_i_type); allocate (l_i_type(0:l_max_neigh_local))
      l_i_type(0:l_max_neigh_local) = tmp_integer(0:l_max_neigh_local)
      deallocate (tmp_integer)

      !resize l_i_type_db
      if (allocated(tmp_integer)) deallocate (tmp_integer); allocate (tmp_integer(0:l_max_neigh_local))
      tmp_integer(1:l_max_neigh_local) = l_i_type_db(1:l_max_neigh_local)
      tmp_integer(0) = l_i_type_db(0)
      if (allocated(l_i_type_db)) deallocate (l_i_type_db); allocate (l_i_type_db(0:l_max_neigh_local))
      l_i_type_db(0:l_max_neigh_local) = tmp_integer(0:l_max_neigh_local)
      deallocate (tmp_integer)


   end subroutine reallocate_neigh_ja

   subroutine preallocate_neigh_ja (l_r_central, l_i_type, l_i_central, l_tmp_dxp, l_tmp_xp, l_imm_neigh)
      use module_kind_variables, only: kind_double
      implicit none

      real(kind_double), dimension(:), allocatable, intent(inout) :: l_r_central
      integer, dimension(:), allocatable, intent(inout) :: l_i_type, l_i_central
      real(kind_double), dimension(:,:), allocatable, intent(inout) :: l_tmp_dxp, l_tmp_xp
      integer, intent(in)  :: l_imm_neigh

      if (allocated(l_i_central)) deallocate (l_i_central); allocate (l_i_central(0:l_imm_neigh))
      if (allocated(l_i_type)) deallocate (l_i_type); allocate (l_i_type(0:l_imm_neigh))

      if (allocated(l_r_central)) deallocate (l_r_central); allocate (l_r_central(l_imm_neigh))
      if (allocated(l_tmp_dxp)) deallocate (l_tmp_dxp); allocate (l_tmp_dxp(3, 0:l_imm_neigh))
      if (allocated(l_tmp_xp)) deallocate (l_tmp_xp); allocate (l_tmp_xp(3, 0:l_imm_neigh))

      l_i_central(:) = 0
      l_i_type(:) = 0
      l_r_central(:) = 0.d0
      l_tmp_dxp(:,:)= 0.d0
      l_tmp_xp = 0.d0

   end subroutine preallocate_neigh_ja


   subroutine build_local_neighbours_ja(  l_iconf, l_ja, l_imm, l_xpnp, r_cut_local, &
      l_d_n_neigh, l_d_kind_neigh, &
      l_r_central, l_i_type, l_i_type_db, l_i_central, l_tmp_dxp, l_tmp_xp, &
      l_max_neigh_local, l_iw2)
#if(MLD_NDM)
      use gen_com_m_ml, ONLY: bg, at!, indi2
      use tab_imm_m_ml, ONLY: realloc_all_tab_imm,dealloc_all_tab_imm,alloc_all_tab_imm
#else
      use ondm_gen_com_m, ONLY: A2cm, bg, at, indi2
      use ondm_tab_imm_m, ONLY: iwmax2
#endif 
      use module_kind_variables, only: kind_double
      use ml_in_ndm_module, only: imm_neigh
      use derived_types, only: config_real

      implicit none

      integer, intent(in)  :: l_iconf, l_ja, l_imm
      real(kind_double), dimension(3, l_imm), intent(in)   :: l_xpnp
      real(kind_double), intent(in) :: r_cut_local
      integer, dimension(l_imm), intent(out)   :: l_d_n_neigh
      integer, dimension(l_imm, imm_neigh), intent(out)    :: l_d_kind_neigh
      real(kind_double), dimension(:), allocatable, intent(inout) :: l_r_central
      integer, dimension(:), allocatable, intent(inout) :: l_i_type, l_i_type_db, l_i_central
      real(kind_double), dimension(:,:), allocatable, intent(inout) :: l_tmp_dxp, l_tmp_xp
      integer, intent(out) :: l_max_neigh_local
      integer, intent(inout) :: l_iw2

      integer  :: iw, iw1, ia_n, ia, max_neigh
      real(kind_double)    :: r_ji, r2_ji
      real(kind_double), dimension(3)  :: dxp_ji, ds
      logical  :: small


      max_neigh = imm_neigh
      small = config_real(l_iconf)%small
      ! begin small box or not 1/
      if (allocated(l_r_central)) deallocate (l_r_central); allocate (l_r_central(max_neigh))
      if (allocated(l_i_central)) deallocate (l_i_central); allocate (l_i_central(0:max_neigh))
      if (allocated(l_i_type)) deallocate (l_i_type); allocate (l_i_type(0:max_neigh))
      if (allocated(l_i_type_db)) deallocate (l_i_type_db); allocate (l_i_type_db(0:max_neigh))
      if (allocated(l_tmp_dxp)) deallocate (l_tmp_dxp); allocate (l_tmp_dxp(3, 0:max_neigh))
      if (allocated(l_tmp_xp)) deallocate (l_tmp_xp); allocate (l_tmp_xp(3, 0:max_neigh))

      l_tmp_dxp(:,0)=0.d0
#if(MLD_NDM)
!!$   if (small) then
#else
      if (small) then 
#endif
         iw1 = 1
         l_iw2 = config_real(l_iconf)%n_neigh(l_ja)
         l_tmp_xp(:, 0) = config_real(l_iconf)%pos_cart(:, l_ja)
#if(MLD_NDM)
!!$    else
!!$      iw1 = l_iw2 + 1
!!$      l_iw2 = iwmax2(l_ja)
!!$      l_tmp_xp(:, 0) = l_xpnp(:, l_ja)
!!$    end if
#else
      else
         iw1 = l_iw2 + 1
         l_iw2 = iwmax2(l_ja)
         l_tmp_xp(:, 0) = l_xpnp(:, l_ja)
      end if
#endif 
      ! end small box or not 1/
#if(MLD_NDM)
!!$   if (small) then
#else
      if (small) then 
#endif
         l_tmp_xp(:, 0) = config_real(l_iconf)%pos_cart(:, l_ja)

#if(MLD_NDM)
!!$    else  
!!$      l_tmp_xp(:, 0) = l_xpnp(:, l_ja)  !debug_neigh !Ang_or_not
!!$    end if 
#else
      else
         l_tmp_xp(:, 0) = l_xpnp(:, l_ja)  !debug_neigh !Ang_or_not
      end if
#endif 
      ia_n = 0
      l_i_central(0) = l_ja
      l_i_type(0) = config_real(l_iconf)%itype(l_ja)
      l_i_type_db(0) = config_real(l_iconf)%itype_db(l_ja)
      do iw = iw1, l_iw2
         ! begin small box or not 2/
#if(MLD_NDM)
!!$   if (small) then
#else
      if (small) then 
#endif
            ia = config_real(l_iconf)%kind_neigh(l_ja, iw)
#if(MLD_NDM)
!!$      else
!!$        ia = indi2(iw)
!!$        if (l_ja == ia) cycle
!!$      end if
#else
         else
            ia = indi2(iw)
            if (l_ja == ia) cycle
         end if
#endif 
         ! end small box or not 2/


#if(MLD_NDM)
!!$   if (small) then
#else
      if (small) then 
#endif
            r_ji = config_real(l_iconf)%r_ij(l_ja, iw)
            r2_ji = r_ji**2
            dxp_ji(:) = config_real(l_iconf)%u_ij(l_ja, iw, :)
#if(MLD_NDM)
!!$      else
!!$        dxp_ji(1:3) = l_xpnp(1:3, ia) - l_xpnp(1:3, l_ja)
!!$        ds = MatMul(dxp_ji, bg)
!!$        WHERE ((ds .GT. 0.5d0) .OR. (ds .LT. -0.5d0))
!!$          ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
!!$        END WHERE
!!$        dxp_ji = MatMul(at, ds)/A2cm
!!$        r2_ji = Sum(dxp_ji(1:3)**2)
!!$        r_ji = dsqrt(r2_ji)
!!$      end if
#else
         else
            dxp_ji(1:3) = l_xpnp(1:3, ia) - l_xpnp(1:3, l_ja)
            ds = MatMul(dxp_ji, bg)
            WHERE ((ds .GT. 0.5d0) .OR. (ds .LT. -0.5d0))
               ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
            END WHERE
            dxp_ji = MatMul(at, ds)/A2cm
            r2_ji = Sum(dxp_ji(1:3)**2)
            r_ji = dsqrt(r2_ji)
         end if
#endif
         if (r_ji >= r_cut_local) cycle
         ia_n = ia_n + 1
         l_r_central(ia_n) = r_ji
         l_i_central(ia_n) = ia
         l_i_type(ia_n) = config_real(l_iconf)%itype(ia)
         l_i_type_db(ia_n) = config_real(l_iconf)%itype_db(ia)
         l_tmp_dxp(:, ia_n) = dxp_ji(:)
#if(MLD_NDM)
!!$   if (small) then
#else
      if (small) then 
#endif
            l_tmp_xp(:, ia_n) = config_real(l_iconf)%u_at(l_ja, iw, :)
#if(MLD_NDM)
!!$      else
!!$        l_tmp_xp(:, ia_n) = l_xpnp(:, ia) !debug_neigh /A2cm
!!$      end if
#else
         else
            ! WTFtbind 
            l_tmp_xp(:, ia_n) = l_xpnp(:, ia)   !debug_neigh /A2cm
         end if
#endif

         l_d_kind_neigh(l_ja, ia_n) = ia
      end do                  ! ia_n
      if (ia_n > max_neigh) then
         write(6,*) 'build_local_neighbours_ja: ia_n=', ia_n, ' > max_neigh=', max_neigh, ' for atom ja=', l_ja
         write(6,*) '  r_cut_local=', r_cut_local, ' small=', small, ' config=', trim(config_real(l_iconf)%filename)
         write(6,*) '  Verlet range: iw1=', iw1, ' l_iw2=', l_iw2, ' span=', l_iw2 - iw1 + 1
         write(6,*) '  rmin=', minval(l_r_central(1:min(ia_n,max_neigh))), &
                    ' rmax=', maxval(l_r_central(1:min(ia_n,max_neigh)))
         stop 'build_local_neighbours_ja: ia_n > imm_neigh. Increase imm_neigh in ml_in_ndm_module.F90'
      end if
      ! is the no of neighbours for the atom l_ja
      l_d_n_neigh(l_ja) = ia_n
      l_max_neigh_local = ia_n

   end subroutine build_local_neighbours_ja

   !subroutine build_local_neighbours_ja_Nfix(  l_iconf, l_ja, l_imm, l_xpnp, r_cut_local, &
   !   l_d_n_neigh, l_d_kind_neigh, &
   !   l_r_central, l_i_type, l_i_type_db, l_i_central, l_tmp_dxp, l_tmp_xp, &
   !   fix_max_neigh_local, l_max_neigh_local, indx_N_fix, fix_N_rcut)

   subroutine build_local_neighbours_ja_Nfix(l_r_central, fix_max_neigh_local,  l_max_neigh_local, l_out_max_neigh_local, indx_N_fix, fix_N_rcut)
#if(MLD_NDM)
      use gen_com_m_ml, ONLY: bg, at!, indi2
      use tab_imm_m_ml, ONLY: realloc_all_tab_imm,dealloc_all_tab_imm,alloc_all_tab_imm
#else
      use ondm_gen_com_m, ONLY: A2cm, bg, at, indi2
      use ondm_tab_imm_m, ONLY: iwmax2
#endif 
      use module_kind_variables, only: kind_double
      use ml_in_ndm_module, only: imm_neigh, discrete_fix_N_rcut, delta_fix_N_rcut !I should add these variables ... delta_fix_N_rcut, discrete_fix_N_rcut 
      use derived_types, only: config_real
      implicit none

      integer, intent(in) :: fix_max_neigh_local
      integer, intent(out) ::  l_out_max_neigh_local
      real(kind_double), dimension(:), allocatable, intent(inout) :: l_r_central
      integer, dimension(:), allocatable :: indx_central

      integer :: ii, ir, jj
      real(kind_double) :: delta_r, tmp_tmp_loss, r_cut_width_Nfix, tmp_rcut, loss_rcut, float_fix_max_neigh_local, tmp_loss, tmp_fcut, tmp_dfcut
      integer, intent(in) :: l_max_neigh_local

      integer, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: indx_N_fix 
      real(kind_double), INTENT(OUT) :: fix_N_rcut
      integer :: ja_fix_max, itest, jj_order    
   
   !case for rcut
   if (l_max_neigh_local > imm_neigh) then 
      write(*,*) 'Stop the calculation l_max_meigh_local too high compared to imm_neigh'
      stop 'decrease imm_neigh ....' 
   end if 


   !write(*,*) 'deb01'  , l_max_neigh_local, fix_max_neigh_local, imm_neigh 
   if (l_max_neigh_local < fix_max_neigh_local) then 
      !debug write(*,*) 'Put the logical array here'
      if (allocated(indx_N_fix)) deallocate (indx_N_fix)
      allocate (indx_N_fix(l_max_neigh_local))  
   
      !call for sorting 
      if (allocated(indx_central)) deallocate (indx_central)
      allocate (indx_central(l_max_neigh_local))
      !if (allocated(l_r_central_copy)) deallocate (l_r_central_copy)
      !allocate (l_r_central_copy(l_max_neigh_local))

      call indexx(l_max_neigh_local, l_r_central(1:l_max_neigh_local), indx_central(1:l_max_neigh_local))

      do ii = 1, l_max_neigh_local 
         indx_N_fix(ii) = indx_central(ii)
      end do

      l_out_max_neigh_local = l_max_neigh_local
      fix_N_rcut = l_r_central(l_max_neigh_local)

   else 
      if (allocated(indx_N_fix)) deallocate (indx_N_fix)
      allocate (indx_N_fix(fix_max_neigh_local)) 

      !call for sorting 
      if (allocated(indx_central)) deallocate (indx_central)
      allocate (indx_central(l_max_neigh_local))
      call indexx(l_max_neigh_local, l_r_central(1:l_max_neigh_local), indx_central(1:l_max_neigh_local))

      do ii = 1, fix_max_neigh_local 
         indx_N_fix(ii) = indx_central(ii) 
      end do
      
      ja_fix_max = indx_N_fix(fix_max_neigh_local) 

      !here is the rcut fitting part...
      !fix number of iterations for the moment ...
      float_fix_max_neigh_local = dble(fix_max_neigh_local)
      
      loss_rcut = 1e4
      !$! delta_r = l_r_central(fix_max_neigh_local)/dble(discrete_fix_N_rcut)
      !$! fix_N_rcut = l_r_central(fix_max_neigh_local)
      !$! delta_r = l_r_central(ja_fix_max)/dble(discrete_fix_N_rcut)
      delta_r = 2.d0 * delta_fix_N_rcut / dble(discrete_fix_N_rcut)
      fix_N_rcut = l_r_central(ja_fix_max)

      r_cut_width_Nfix = r_cut_width

      do ir = 1, discrete_fix_N_rcut ! to be defined 
         !$! tmp_rcut = (1.d0 - delta_fix_N_rcut)*l_r_central(ja_fix_max) + ir*2.d0*delta_fix_N_rcut*delta_r
         tmp_rcut = (l_r_central(ja_fix_max) - delta_fix_N_rcut) + delta_r * dble(ir) 
         tmp_loss = 0.d0
         do jj = 1, l_max_neigh_local 
            jj_order = indx_central(jj)
            tmp_fcut = 0.d0
            call fcut_rij(3, l_r_central(jj_order), tmp_rcut + r_cut_width_Nfix, r_cut_width_Nfix, .false., tmp_fcut, tmp_dfcut)
            tmp_loss = tmp_loss + tmp_fcut      
         end do

         if (abs(tmp_loss-float_fix_max_neigh_local) <= loss_rcut) then 
            loss_rcut = abs(tmp_loss-float_fix_max_neigh_local)
         end if

         itest = 0 
         if (loss_rcut <= 1.d0  ) then 
            fix_N_rcut = tmp_rcut 
            tmp_tmp_loss = tmp_loss
            itest= 1 
         end if   

         !d write(*,'("aaair ", 2i5, 4f14.6)') ir, itest, loss_rcut, tmp_loss,  tmp_rcut, fix_N_rcut

         if (itest == 1) exit 
      end do
      if (itest == 0) tmp_tmp_loss = tmp_loss
      !what is that .... 
      if (loss_rcut > 0.6d0) then 
         fix_N_rcut = l_r_central(ja_fix_max)
      end if

      tmp_loss = 0.d0
      tmp_rcut = fix_N_rcut
      do jj = 1, l_max_neigh_local 
            jj_order = indx_central(jj)
            tmp_fcut = 0.d0
            call fcut_rij(3, l_r_central(jj_order), tmp_rcut + r_cut_width_Nfix, r_cut_width_Nfix, .false., tmp_fcut, tmp_dfcut)
            tmp_loss = tmp_loss + tmp_fcut 
            !d write(*,'("ccccc ", i5, 5F18.7)')  jj, tmp_loss, tmp_tmp_loss, tmp_rcut 
            if ((tmp_loss) >= tmp_tmp_loss) then
              !d write(*,'("dddd ", i5, 2F18.7)')  jj, tmp_loss, tmp_tmp_loss
              l_out_max_neigh_local = jj 
              exit   
            end if      
      end do
      !d write(*,'("oooo:  ",2i5, 2f20.10)') fix_max_neigh_local, l_out_max_neigh_local, tmp_rcut, fix_N_rcut
      !d write(*,'(">>>,"2i5, 4f20.10)') fix_max_neigh_local, l_out_max_neigh_local, tmp_rcut, fix_N_rcut, tmp_loss, loss_rcut

   !write(*,*) fix_N_rcut, l_r_central(fix_max_neigh_local), loss_rcut

   end if 

end subroutine build_local_neighbours_ja_Nfix

   subroutine build_cos_and_fcut_ja(l_type_fcut, l_desc_forces_local, &
      l_r_cut, l_r_cut_width, l_max_neigh_local, l_r_central, l_tmp_dxp, &
      l_tmpcos_dxp, l_r_fcut, l_d_r_fcut)
      use module_kind_variables, only : kind_double
      use ml_in_ndm_module, only: one_pi
      !use module_neigh_local, only: r_cut
      !use module_neigh_local, only: max_neigh_local, tmp_dxp, r_central, tmpcos_dxp, &
      !                              d_r_central, d_r_fcut, r_fcut
      implicit none
      integer, intent(in) :: l_type_fcut
      real(kind_double), intent(inout)  ::  l_r_cut, l_r_cut_width
      integer, intent(in)  :: l_max_neigh_local
      real(kind_double), dimension(:), allocatable, intent(in) :: l_r_central
      real(kind_double), dimension(:,:), allocatable, intent(in) :: l_tmp_dxp
      real(kind_double), dimension(:), allocatable, intent(out) :: l_r_fcut
      real(kind_double), dimension(:,:), allocatable, intent(out) :: l_tmpcos_dxp, l_d_r_fcut
      logical, intent(in)  :: l_desc_forces_local

      integer  :: ii
      real(kind_double), dimension(:), allocatable :: xtmp
      real(kind_double) :: tmprr

      ! cosine dx/dr .....
      if (l_desc_forces_local) then
         if (allocated(l_tmpcos_dxp)) deallocate (l_tmpcos_dxp); allocate (l_tmpcos_dxp(3, l_max_neigh_local))
         do ii = 1, 3
            l_tmpcos_dxp(ii, 1:l_max_neigh_local) = l_tmp_dxp(ii, 1:l_max_neigh_local)/l_r_central(1:l_max_neigh_local)
         end do
      end if

      !fcut and derivatives .................
      if (allocated(l_r_fcut)) deallocate (l_r_fcut); allocate (l_r_fcut(l_max_neigh_local))
      if (l_desc_forces_local) then
         if (allocated(l_d_r_fcut)) deallocate (l_d_r_fcut); allocate (l_d_r_fcut(3, l_max_neigh_local))
         !if (allocated(l_d_r_central)) deallocate (l_d_r_central); allocate (l_d_r_central(3, l_max_neigh_local))
      end if

      if (l_type_fcut == 1) then
         ! ((r_ij / r_cut)^2 - 1)^2

         if (allocated(xtmp)) deallocate (xtmp); allocate (xtmp(l_max_neigh_local))
         xtmp(1: l_max_neigh_local) = (l_r_central(1:l_max_neigh_local)/l_r_cut)**2 - 1
         l_r_fcut(1:l_max_neigh_local) = xtmp(1:l_max_neigh_local)**2
         if (l_desc_forces_local) then
            do ii = 1, 3
               l_d_r_fcut(ii, 1:l_max_neigh_local) = l_tmpcos_dxp(ii, 1:l_max_neigh_local)* &
                  4.d0* xtmp(1:l_max_neigh_local) *l_r_central(1:l_max_neigh_local)/l_r_cut**2
            end do
         end if

      else if (l_type_fcut==2) then

         do ii = 1, l_max_neigh_local
            tmprr = one_pi*l_r_central(ii)/l_r_cut
            l_r_fcut(ii) = 0.5d0*(cos(tmprr) + 1.d0)
            if (l_desc_forces_local)  l_d_r_fcut(1:3,ii) = -0.5d0*one_pi/l_r_cut*sin(tmprr)*l_tmpcos_dxp(1:3, ii)
         end do


      else if (l_type_fcut==3) then

         do ii = 1, l_max_neigh_local
            tmprr = one_pi*(l_r_central(ii)  - l_r_cut + l_r_cut_width)/l_r_cut_width
            if (l_r_central(ii) < (l_r_cut -l_r_cut_width)) then
               l_r_fcut(ii) =  1.d0
               if (l_desc_forces_local)  l_d_r_fcut(1:3,ii) = 0.d0
            else
               l_r_fcut(ii) = 0.5d0*(cos(tmprr) + 1.d0)
               if (l_desc_forces_local)  l_d_r_fcut(1:3,ii) = -0.5d0*one_pi/l_r_cut_width*sin(tmprr)*l_tmpcos_dxp(1:3, ii)
            end if
         end do

      else if (l_type_fcut==4) then

         if (allocated(xtmp)) deallocate (xtmp); allocate (xtmp(l_max_neigh_local))
         xtmp(1: l_max_neigh_local) =  l_r_central(1:l_max_neigh_local) - l_r_cut
         l_r_fcut(1:l_max_neigh_local) = xtmp(1:l_max_neigh_local)**2
         if (l_desc_forces_local) then
            do ii = 1,3
               l_d_r_fcut(ii, 1:l_max_neigh_local) = 2.d0*l_tmpcos_dxp(ii, 1:l_max_neigh_local)*xtmp(1:l_max_neigh_local)
            end do
         end if



      end if


   end subroutine build_cos_and_fcut_ja

   subroutine fcut_rij(l_type_fcut, l_rr, l_r_cut, l_rcut_width, l_desc_forces_local, l_fcut, l_dfcut)
      use module_kind_variables, only: kind_double
      use ml_in_ndm_module, only: one_pi

      implicit none
      real(kind_double), intent(in) :: l_rr, l_r_cut, l_rcut_width
      real(kind_double), intent(out) :: l_fcut, l_dfcut
      logical, intent(in) :: l_desc_forces_local
      integer, intent(in) :: l_type_fcut

      real(kind_double) :: xx

      if ( l_rr > l_r_cut) then
         l_fcut = 0.d0
         l_dfcut = 0.d0
      else

         select case(l_type_fcut)
            !
          case(1)
            xx = (l_rr/l_r_cut)**2 - 1.d0
            l_fcut = xx**2
            if (l_desc_forces_local) l_dfcut = 4.d0*xx/l_r_cut**2

          case(2)
            xx =  (l_rr/l_r_cut)
            l_fcut = 0.5d0*(cos (one_pi * xx ) + 1.d0)
            if (l_desc_forces_local) l_dfcut = - 0.5d0 * one_pi * sin(one_pi*xx)/l_r_cut
            !write(*,*) 'fcut', rangml, xx, sin(pi*xx), l_dfcut, l_desc_forces_local
          case (3)
            if (l_rr < (l_r_cut - l_rcut_width)) then
               l_fcut = 1.d0
               l_dfcut = 0.d0
            else
               xx = one_pi*(l_rr - l_r_cut + l_rcut_width)/l_rcut_width
               l_fcut = 0.5d0*(cos(xx) + 1.d0)
               if (l_desc_forces_local) l_dfcut = -0.5d0* sin(xx) * one_pi / l_rcut_width
            end if

          case(4)
            xx = l_rr - l_r_cut
            l_fcut = xx*2
            if (l_desc_forces_local) l_dfcut = 2.d0 * xx

         end select

      end if
   end subroutine fcut_rij

   subroutine fcut_rij_inout( type_fcut_in, type_fcut_out, l_rr, l_r_cut, l_rcut_width, l_rcut_in, l_rcut_width_in, l_desc_forces_local, l_fcut, l_dfcut)
     use module_kind_variables, only: kind_double

     implicit none
     real(kind_double), intent(in) :: l_rr, l_r_cut, l_rcut_width, l_rcut_in, l_rcut_width_in
     real(kind_double), intent(out) :: l_fcut, l_dfcut 
     integer, intent(in) :: type_fcut_in, type_fcut_out
     real(kind_double) :: fcut_out, fcut_in, dfcut_out, dfcut_in
     logical, intent(in) :: l_desc_forces_local 


     call fcut_rij(type_fcut_out, l_rr, l_r_cut, l_rcut_width, l_desc_forces_local, fcut_out, dfcut_out)

     if (l_rcut_in > 0 ) then 
        call fcut_rij(type_fcut_in, l_rr, l_rcut_in + l_rcut_width_in, l_rcut_width_in, l_desc_forces_local, fcut_in, dfcut_in)
     else 
        fcut_in = 0.d0 
        dfcut_in = 0.d0 
     end if
     l_fcut = (1.d0 - fcut_in) * fcut_out 
     if (l_desc_forces_local) l_dfcut =-dfcut_in * fcut_out +  (1.d0 - fcut_in) * dfcut_out

   end subroutine fcut_rij_inout


   subroutine fcut_rij_second(l_type_fcut, l_rr, l_r_cut, l_rcut_width, l_fcut, l_dfcut, l_ddfcut)
      use module_kind_variables, only: kind_double
      use ml_in_ndm_module, only: one_pi

      implicit none
      real(kind_double), intent(in) :: l_rr, l_r_cut, l_rcut_width
      real(kind_double), intent(out) :: l_fcut, l_dfcut, l_ddfcut
      integer, intent(in) :: l_type_fcut

      real(kind_double) :: xx, inv_cut

      if ( l_rr > l_r_cut) then
         l_fcut = 0.d0
         l_dfcut = 0.d0
         l_ddfcut = 0.d0
      else

         select case(l_type_fcut)
            !
          case(1)
            inv_cut = 1.d0/l_r_cut**2
            xx = (l_rr/l_r_cut)**2 - 1.d0
            l_fcut = xx**2
            l_dfcut = 4.d0*l_rr*xx*inv_cut
            l_ddfcut = 4.d0* inv_cut *(2.d0*l_rr**2*inv_cut + xx )

          case(2)
            xx =  (l_rr/l_r_cut)
            l_fcut = 0.5d0*(cos (one_pi * xx ) + 1.d0)
            l_dfcut = - 0.5d0 * one_pi * sin(one_pi*xx)/l_r_cut
            l_ddfcut = - 0.5d0 * one_pi**2 * cos (one_pi * xx ) / l_r_cut**2

          case (3)
            if (l_rr < (l_r_cut - l_rcut_width)) then
               l_fcut = 1.d0
               l_dfcut = 0.d0
               l_ddfcut = 0.d0
            else
               xx = one_pi*(l_rr - l_r_cut + l_rcut_width)/l_rcut_width
               l_fcut = 0.5d0*(cos(xx) + 1.d0)
               l_dfcut = -0.5d0* sin(xx) * one_pi / l_rcut_width
               l_ddfcut = - 0.5d0* one_pi**2 * cos(xx) / l_rcut_width**2
            end if

          case(4)
            xx = l_rr - l_r_cut
            l_fcut = xx*2
            l_dfcut = 2.d0 * xx
            l_ddfcut = 2.d0
         end select

      end if
   end subroutine fcut_rij_second


   subroutine compute_ur_transformed_distances(l_desc_forces_local, l_max_neigh_local, l_r_central, l_tmp_dxp, l_tmpcos_dxp, l_ur_central, d_l_ur_central)
      use module_kind_variables, only: kind_double
      use module_body_desc, only: bond_dist_pure, bond_dist_exp, bond_dist_inverse, bond_dist_transform, &
         bond_beta, bond_dist_ann
      implicit none

      integer, intent(in) :: l_max_neigh_local
      logical, intent(in) :: l_desc_forces_local
      real(kind_double), dimension(:), intent(in) :: l_r_central
      real(kind_double), dimension(:,:), intent(in) :: l_tmp_dxp
      real(kind_double), dimension(:,:), allocatable, intent(inout) :: l_tmpcos_dxp
      real(kind_double), dimension(:), allocatable, intent(inout) :: l_ur_central
      real(kind_double), dimension(:,:), allocatable, intent(inout)  :: d_l_ur_central

      integer  :: ii

      ! cosine dx/dr .....
      if (l_desc_forces_local) then
         if (allocated(l_tmpcos_dxp)) deallocate (l_tmpcos_dxp); allocate (l_tmpcos_dxp(3, l_max_neigh_local))
         do ii = 1, 3
            l_tmpcos_dxp(ii, 1:l_max_neigh_local) = l_tmp_dxp(ii, 1:l_max_neigh_local)/l_r_central(1:l_max_neigh_local)
         end do
      end if

      ! compute l_ur_central ... the transformation of l_r_central ...
      if (allocated(l_ur_central)) deallocate (l_ur_central); allocate (l_ur_central(l_max_neigh_local))
      ! transform the radial function ...
      if (l_desc_forces_local) then
         if (allocated(d_l_ur_central)) deallocate (d_l_ur_central); allocate (d_l_ur_central(3, l_max_neigh_local))
      end if
      select case (bond_dist_transform)
       case (bond_dist_pure)
         l_ur_central(1:l_max_neigh_local) = l_r_central(1:l_max_neigh_local)
         if (l_desc_forces_local) then
            do ii = 1, 3
               d_l_ur_central(ii, 1:l_max_neigh_local) = l_tmpcos_dxp(ii, 1:l_max_neigh_local)
            end do
         end if
       case (bond_dist_exp)
         l_ur_central(1:l_max_neigh_local) = dexp(-bond_beta*l_r_central(1:l_max_neigh_local))
         if (l_desc_forces_local) then
            do ii = 1, 3
               d_l_ur_central(ii, 1:l_max_neigh_local) = l_tmpcos_dxp(ii, 1:l_max_neigh_local)*(-bond_beta)*l_ur_central(1:l_max_neigh_local)
            end do
         end if
       case (bond_dist_inverse)
         l_ur_central(1:l_max_neigh_local) = bond_dist_ann**bond_beta/l_r_central(1:l_max_neigh_local)**bond_beta
         if (l_desc_forces_local) then
            do ii = 1, 3
               d_l_ur_central(ii, 1:l_max_neigh_local) = (-bond_beta)*l_tmpcos_dxp(ii, 1:l_max_neigh_local)*l_ur_central(1:l_max_neigh_local)/l_r_central(1:l_max_neigh_local)
            end do
         end if
      end select

   end subroutine compute_ur_transformed_distances

end module module_neigh_local



subroutine test_if_config_is_small(iconf)
   ! test is a configuration is small comapred to r_cut.
   ! Input:
   !    r_cut 
   !    config_real(iconf)%cell
   ! Output:
   !        config_real(iconf)%bg_cell - reciprocal cell 
   !         config_real(iconf)%small
   !                                 T - box is small
   !                                 F - box is large
   !       config_real(iconf)%nxCell   - integer number of cells in x direction                        
   !       config_real(iconf)%nyCell   - integer number of cells in y direction
   !       config_real(iconf)%nzCell   - integer number of cells in z direction

   use module_kind_variables, only: kind_double
   use ml_in_ndm_module, only: rangml, debug
   use module_neigh_local, only: r_cut
   use derived_types, only: config_real
   use mld_logger
#ifdef MLD_NDM
   use recips_mod, only: recips
#else
   use ondm_transform_coord, only: ondm_recips
#endif

   implicit none

   _NAMECURRENT_("test_config_is_small")

   integer, intent(in)  :: iconf
   real(kind_double)    :: bval(3)
   integer  :: i


   _MLD_BEGIN_
#ifdef MLD_NDM
  call recips(config_real(iconf)%cell(1, 1), config_real(iconf)%cell(1, 2), config_real(iconf)%cell(1, 3), &
              config_real(iconf)%bg_cell(1, 1), config_real(iconf)%bg_cell(1, 2), config_real(iconf)%bg_cell(1, 3))
#else
   call ondm_recips(config_real(iconf)%cell(1, 1), config_real(iconf)%cell(1, 2), config_real(iconf)%cell(1, 3), &
      config_real(iconf)%bg_cell(1, 1), config_real(iconf)%bg_cell(1, 2), config_real(iconf)%bg_cell(1, 3))
#endif

   do i = 1, 3
      bval(i) = 1.d0/sqrt(sum(config_real(iconf)%bg_cell(:, i)**2))
   end do
   config_real(iconf)%nxCell = int(2.d0*r_cut/bval(1)) + 1
   config_real(iconf)%nyCell = int(2.d0*r_cut/bval(2)) + 1
   config_real(iconf)%nzCell = int(2.d0*r_cut/bval(3)) + 1

   if (max(config_real(iconf)%nxCell, config_real(iconf)%nyCell, config_real(iconf)%nzCell) .gt. 1) config_real(iconf)%small = .true.

   if (debug) then
      if (rangml == 0) write (*, '("ML: small or big box in test_if_config_is_small ",  i5,  3i4, l3)') iconf, config_real(iconf)%nxCell, &
         config_real(iconf)%nyCell, &
         config_real(iconf)%nzCell, &
         config_real(iconf)%small
   end if

   _MLD_END_
end subroutine test_if_config_is_small




subroutine calc_neighbours(iconf)
   ! compute neighbours using milady and NDM style
#if(MLD_NDM) 
  use NDM_ML,only:neighbours_large
#endif 
   use ml_in_ndm_module, only: debug, rangml
   use derived_types, only: config_real
   use mld_subworld
   use mld_logger

   implicit none

   integer, intent(in)  :: iconf

   _NAMECURRENT_("calc_neighbours")



   _MLD_BEGIN_


   !<<<starting the neighbours calculations
#if(MLD_NDM)
  if (debug) then
    if (rangml == 0) write (6, '("ML: the start of neighbours calc of:", a)') config_real(iconf)%filename
  end if
  call neighbours_ndm_layer(config_real(iconf))
  if (config_real(iconf)%small) then
    call neighbours_small(config_real(iconf))
  else
    call neighbours_large(config_real(iconf))
  end if
#else 
   if (debug) then
      if (subrank == 0) write (6, '("ML: the start of neighbours calc of:", a)') config_real(iconf)%filename
   end if
   call neighbours_ndm_layer(iconf)
   if (config_real(iconf)%small) then
      call neighbours_small(iconf)
   else
      call neighbours_large(iconf)
   end if
#endif 
   !<<<ending the neighbours calculations

   if (debug) then
      if (rangml == 0) write (6, '("ML: the out of neighbours calc of:", a)') config_real(iconf)%filename
   end if

   _MLD_END_

end subroutine calc_neighbours

#if(MLD_NDM)
subroutine neighbours_ndm_layer(cn2m)
#else
subroutine neighbours_ndm_layer(iconf)
#endif

   ! set-up the ndm using the values imported from iconf
   use ml_in_ndm_module, only: debug, rangml, ML_MLD_DMTYPE
   use module_neigh_local, only: r_cut
   use mld_logger
#if(MLD_NDM)
  use gen_com_m_ml, only:  im_glob, im, imm, at, bg,rvois
  use gen_com_m, only: dmtype, A2cm, umass,pi!, im_glob, im, imm, at, bg, A2cm, rvois, umass
  ! CETTE LIGNE POSE PROBLEME :
  use boxconfig,only:initbox
  USE setcell,only:setnox,setcellconf
  use derived_types,only:system_state
  use tab_imm_m_ml,only:xp,fp,ityp
  use tab_imm_m_ml, ONLY: realloc_all_tab_imm,dealloc_all_tab_imm,alloc_all_tab_imm
  use gen_com_m_ml,only:at,bg
  use cellconfig,only:caltabtc
  use NDM_ML,only  : atndm, conf_real2atconfig, cellndm,boxndm
!$
  use var_pot, only: ntyp,cm
!  use tab_imm_m_ml,only:xpp,ielat,ax,ityp,xp,fp,realloc_all_tab_imm,dealloc_all_tab_imm,alloc_all_tab_imm
#else
   use ondm_gen_com_m, only: dmtype, im_glob, im, imm, at, bg, A2cm, rvois, umass
   use derived_types, only: config_real
   use ondm_var_pot, only: ntyp, cm, rumax, na
   use ondm_tab_imm_m
#endif

   implicit none

#if(MLD_NDM)
  type(system_state)::cn2m 
!  real(kind(0.d0)), dimension(size(xpp, 1), size(xpp, 2))  :: xpp_copy, ax_copy
!  integer, dimension(size(ielat))  :: ielat_copy
  integer, dimension(cn2m%ntypes)      :: na_copy
  integer::ipbc(3),nvt,nvperat
#else
   integer, intent(in)  :: iconf
   !convert xp, at from A to cm
   real(kind(0.d0)), dimension(size(xpp, 1), size(xpp, 2))  :: xpp_copy, ax_copy
   integer, dimension(size(ielat))  :: ielat_copy
   integer, dimension(config_real(iconf)%ntypes)      :: na_copy
#endif


   real(kind(0.d0))     :: rumax_copy

   _NAMECURRENT_("neighbours_ndm_layer")

   _MLD_BEGIN_
  
#if(MLD_NDM) 
  ipbc(:)=1
  imm = cn2m%nat
  im = cn2m%nat
  im_glob = im            ! Mise a jour de im_glob pour divid
  ntyp = cn2m%ntypes

!!$  if (dmtype /= 18) then
!!$    xpp_copy(:, :) = xpp(:, :)
!!$    ax_copy(:, :) = ax(:, :)
!!$    ielat_copy(:) = ielat(:)
!!$    rumax_copy = rumax
!!$    na_copy(:ntyp) = na(:ntyp)
!!$  end if
#else
   imm = config_real(iconf)%nat
   im = config_real(iconf)%nat
   im_glob = im            ! Mise a jour de im_glob pour divid
   ntyp = config_real(iconf)%ntypes

   if (dmtype /= ML_MLD_DMTYPE) then
      xpp_copy(:, :) = xpp(:, :)
      ax_copy(:, :) = ax(:, :)
      ielat_copy(:) = ielat(:)
      rumax_copy = rumax
      !$! if (ALLOCATED(na_copy)) write(*,*) 'na_copy A'
      !$! if (ASSOCIATED(na)) write(*,*) 'na   A'
      !$! na_copy(:ntyp) = na(:ntyp)
      !call dealloc_all_tab_imm
      !call alloc_all_tab_imm(im)
   end if
#endif


   if (dmtype == ML_MLD_DMTYPE) then
      ! why I do not know TODO
      call dealloc_all_tab_imm
      call alloc_all_tab_imm(im) 
   else
      !TODOmd call realloc_all_tab_imm(im)
      call dealloc_all_tab_imm
      call alloc_all_tab_imm(im) 
   end if

#if(MLD_NDM) 
!!$  if (dmtype /= 18) then
!!$    ax(:, :) = ax_copy(:, :)
!!$    xpp(:, :) = xpp_copy(:, :)
!!$    ielat(:) = ielat_copy(:)
!!$  end if

  !debug write (*,*) 'IN NEIGH10', im, imm, ntyp,xp(1,1), cn2m%pos_cart(1,1)
  ityp(:) = cn2m%itype(:)
  xp(:, :) = cn2m%pos_cart(:, :)
  fp(:, :) = cn2m%force(:, :)
  at(:, :) = cn2m%cell(:, :)
  bg(:, :) = cn2m%bg_cell(:, :)

  if (debug) then
    if (rangml == 0) write (6, *) 'The cell', cn2m%cell(:, :)
  end if
#else
   if (dmtype /= ML_MLD_DMTYPE) then
      ax(:, :) = ax_copy(:, :)
      xpp(:, :) = xpp_copy(:, :)
      ielat(:) = ielat_copy(:)
   end if

   !debug write (*,*) 'IN NEIGH10', im, imm, ntyp,xp(1,1), config_real(iconf)%pos_cart(1,1)
   ityp(:) = config_real(iconf)%itype(:)
   xp(:, :) = config_real(iconf)%pos_cart(:, :)
   fp(:, :) = config_real(iconf)%force(:, :)
   at(:, :) = config_real(iconf)%cell(:, :)
   bg(:, :) = config_real(iconf)%bg_cell(:, :)
   if (debug) then
      if (rangml == 0) write (6, *) 'The cell', config_real(iconf)%cell(:, :)
   end if
#endif 

  
   !C_DEBUG cutoff  for the neighbours list update.
   rvois = r_cut*A2cm
   call convert_A2cm(1, .true., .true.)
   call alloc_typ_ml()

#if(MLD_NDM)

  if (dmtype /= ML_MLD_DMTYPE) then
!    na(:ntyp) = na_copy(:ntyp)
 !   rumax = rumax_copy
  end if
  cm(:) = cn2m%mass_per_type(:)*umass


!atndm,boxndm and cellndm are variables of this whole module. Will not work for  multiple simultaneous instances of  cn2m.
  ntyp = cn2m%ntypes
  nvperat = 4*Pi*(r_cut + 1.0)**3/3
  nvt=cn2m%nat*nvperat
  if (cn2m%small) then
     call atndm%init(imin=cn2m%nat,ltabvois=.false.,lreallocate=.true.,im_glob=cn2m%nat)
  else
     call atndm%init(imin=cn2m%nat,ltabvois=.true.,rvois=r_cut,lreallocate=.true.,nvois=nvt,im_glob=cn2m%nat)
  end if
  !write(6,*)'PCRC1'
  call conf_real2atconfig(atndm,cn2m,'xi')
  !write(6,*)'PCRC1.1'
!  call atndm%print
  call initbox(boxndm,cn2m%cell,ipbc)
!  call boxndm%print
  call cellndm%dealloc
  !write(6,*)' PCRC RRRRCCCCUUUUTTTT',r_cut, atndm%rvois
  call setnox(boxndm,cellndm,r_cut,lverbose=.false.)
  call setcellconf(cellndm,atndm,boxndm,r_cut,lverbose=.false.)
  call caltabtC(cellndm,atndm,.false.,boxndm)
!  call cellndm%print

#else

   if (dmtype /= ML_MLD_DMTYPE) then
      na(:ntyp) = na_copy(:ntyp)
      rumax = rumax_copy
   end if
   cm(:) = config_real(iconf)%mass_per_type(:)*umass
   !debug write (*,*) 'IN NEIGH3', rumax, rvois, im, imm, ntyp,xp(1,1), config_real(iconf)%pos_cart(1,1)
#endif 
   _MLD_END_

end subroutine neighbours_ndm_layer



subroutine neighbours_large(iconf)

#ifdef MLD_NDM
   use gen_com_m_ml, only: natperc, nox, noy, noz
#else
   use ondm_gen_com_m, only: natperc, nox, noy, noz
#endif
   use ml_in_ndm_module, only: debug, rangml
   use mld_logger

   implicit none

   integer, intent(in)  :: iconf
   integer :: izozo 

   _NAMECURRENT_("neighbours_large")



   _MLD_BEGIN_
   izozo = iconf 
   ! write NAMECURRENT, iconf

   !impact! call Deallocatecel()
   call DeallocateVeryAll_varpot
   call DeallocateVeryAll_gencomm

   nox = -1; noy = -1; noz = -1
   natperc = -1            ! Force le calcul de natperc dans divid

   call ondm_divid(0)

   if (debug) then
      if (rangml == 0) write (6, '("ML:.............divid 0   :",3i6)') nox, noy, noz
   end if

   call ondm_divid(1)

   if (debug) then
      if (rangml == 0) write (6, '("ML:.............divid 1   :",3i6)') nox, noy, noz
   end if

   call ondm_DynamicalAllocationCell()                   ! Reallocation du pointeur last(natperc,:noxyz) pour caltabt
   call ondm_neigcel()
   call ondm_caltabt()
   call ondm_caltabi()

   _MLD_END_
end subroutine neighbours_large




#if(MLD_NDM) 
subroutine neighbours_small(cn2m)
#else
subroutine neighbours_small(iconf)
#endif
   use module_kind_variables, ONLY: kind_double
   use ml_in_ndm_module, only: imm_neigh
   use module_neigh_local, only: r_cut
#if(MLD_NDM)
  use gen_com_m_ml, only: im
  use derived_types, only: system_state
#else
   use ondm_gen_com_m, only: im
   use derived_types, only: config_real
#endif
   use mld_logger
   use mld_string, only: vtoa

   implicit none

#if(MLD_NDM)
!  integer, intent(in)  :: iconf
  type(system_state)  :: cn2m
  integer :: imloc 

#else
   integer, intent(in)  :: iconf
#endif
   integer  :: ia, ja, c, k, n1, n2, n3, ibox
   double precision     :: r2, r_cut2
   real(kind_double), dimension(:, :), allocatable    :: xpnp
   integer  :: nsize1, nsize2, nsize3
   real(kind_double), dimension(1:3)      :: utemp

   _NAMECURRENT_("neighbours_small")



   _MLD_BEGIN_

#if(MLD_NDM) 
  imloc=cn2m%nat
  if (allocated(cn2m%n_neigh)) deallocate (cn2m%n_neigh); allocate (cn2m%n_neigh(imloc))
  if (allocated(cn2m%r_ij)) deallocate (cn2m%r_ij); allocate (cn2m%r_ij(imloc, imm_neigh))
  if (allocated(cn2m%u_per)) deallocate (cn2m%u_per); allocate (cn2m%u_per(imloc, imm_neigh, 3))
  if (allocated(cn2m%u_at)) deallocate (cn2m%u_at); allocate (cn2m%u_at(imloc, imm_neigh, 3))
  if (allocated(cn2m%type_neigh)) deallocate (cn2m%type_neigh); allocate (cn2m%type_neigh(imloc, imm_neigh))
  if (allocated(cn2m%kind_neigh)) deallocate (cn2m%kind_neigh); allocate (cn2m%kind_neigh(imloc, imm_neigh))
  if (allocated(cn2m%incell)) deallocate (cn2m%incell); allocate (cn2m%incell(imloc, imm_neigh))

  if (allocated(cn2m%u_ij)) deallocate (cn2m%u_ij); allocate (cn2m%u_ij(imloc, imm_neigh, 3))

  ALLOCATE (xpnp(3, cn2m%nat))
  ! if (lperiod) then
  xpnp(:, :) = cn2m%pos_cart
  ! else
  !   call notperiod(xp,xpnp)
  ! end if
  ! call cryst_to_cart (imm, xpnp, cn2m%bg_cell, -1)


  nsize1 = int(cn2m%nxCell/2) + 1
  nsize2 = int(cn2m%nyCell/2) + 1
  nsize3 = int(cn2m%nzCell/2) + 1


#else
   if (allocated(config_real(iconf)%n_neigh)) deallocate (config_real(iconf)%n_neigh); allocate (config_real(iconf)%n_neigh(im))
   if (allocated(config_real(iconf)%r_ij)) deallocate (config_real(iconf)%r_ij); allocate (config_real(iconf)%r_ij(im, imm_neigh))
   if (allocated(config_real(iconf)%u_per)) deallocate (config_real(iconf)%u_per); allocate (config_real(iconf)%u_per(im, imm_neigh, 3))
   if (allocated(config_real(iconf)%u_at)) deallocate (config_real(iconf)%u_at); allocate (config_real(iconf)%u_at(im, imm_neigh, 3))
   if (allocated(config_real(iconf)%type_neigh)) deallocate (config_real(iconf)%type_neigh); allocate (config_real(iconf)%type_neigh(im, imm_neigh))
   if (allocated(config_real(iconf)%kind_neigh)) deallocate (config_real(iconf)%kind_neigh); allocate (config_real(iconf)%kind_neigh(im, imm_neigh))
   if (allocated(config_real(iconf)%incell)) deallocate (config_real(iconf)%incell); allocate (config_real(iconf)%incell(im, imm_neigh))
   if (allocated(config_real(iconf)%u_ij)) deallocate (config_real(iconf)%u_ij); allocate (config_real(iconf)%u_ij(im, imm_neigh, 3))

   ALLOCATE (xpnp(3, config_real(iconf)%nat))
   ! if (lperiod) then
   xpnp(:, :) = config_real(iconf)%pos_cart
   ! else
   !   call ondm_notperiod(xp,xpnp)
   ! end if
   ! call cryst_to_cart (imm, xpnp, config_real(iconf)%bg_cell, -1)


   nsize1 = int(config_real(iconf)%nxCell/2) + 1
   nsize2 = int(config_real(iconf)%nyCell/2) + 1
   nsize3 = int(config_real(iconf)%nzCell/2) + 1
#endif

   r_cut2 = r_cut**2

#if(MLD_NDM)
  do ia = 1, cn2m%nat
    cn2m%incell(ia, :) = .false.
#else
   do ia = 1, config_real(iconf)%nat
      config_real(iconf)%incell(ia, :) = .false.
#endif
      c = 0
      ibox = 0
#if(MLD_NDM)
    do ja = 1, cn2m%nat
      ! dxp_ji(:) = cn2m%pos_cart(:, ia) - cn2m%pos_cart(:, ja)
      !    ds(:) = MatMul(dxp_ji, cn2m%bg_cell)
#else
      do ja = 1, config_real(iconf)%nat
         ! dxp_ji(:) = config_real(iconf)%pos_cart(:, ia) - config_real(iconf)%pos_cart(:, ja)
         !    ds(:) = MatMul(dxp_ji, config_real(iconf)%bg_cell)
#endif
         !    WHERE ((ds .GT. 0.5d0) .OR. (ds .LT. -0.5d0))
         !      ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
         !    END WHERE
#if(MLD_NDM)
      !    dxp_ji(:) = MatMul(cn2m%cell, ds)
#else
         !    dxp_ji(:) = MatMul(config_real(iconf)%cell, ds)
#endif
         !    r2_ji = SUM(dxp_ji(:)**2)

         do n1 = -nsize1, nsize1
            do n2 = -nsize2, nsize2
               do n3 = -nsize3, nsize3
                  r2 = 0.d0
                  do k = 1, 3
#if(MLD_NDM)
              utemp(k) = xpnp(k, ia) - (xpnp(k, ja) + cn2m%cell(k, 1)*dble(n1) + cn2m%cell(k, 2)*dble(n2) + cn2m%cell(k, 3)*dble(n3))

#else
                     utemp(k) = xpnp(k, ia) - (xpnp(k, ja) + config_real(iconf)%cell(k, 1)*dble(n1) + config_real(iconf)%cell(k, 2)*dble(n2) + config_real(iconf)%cell(k, 3)*dble(n3))
#endif
                     r2 = r2 + utemp(k)**2
                  end do
                  if (r2 .lt. 1d-15) cycle
                  if (r2 .gt. r_cut2) cycle
                  c = c + 1
                  !          if (dabs(r2 - r2_ji) .lt. 1.d-15) then
                  !            if (n1 == 0) .or. (n2 == 0)
                  !            ibox = ibox + 1
                  !            debug write (*,'("nnn", 4i5)') ia, ja, ibox, c
                  !          end if
#if(MLD_NDM)
            cn2m%type_neigh(ia, c) = cn2m%itype(ja)
            cn2m%kind_neigh(ia, c) = ja
            cn2m%incell(ia, c) = .true.
            cn2m%r_ij(ia, c) = dsqrt(r2)
            cn2m%u_per(ia, c, :) = -utemp(:) - xpnp(:, ja) + xpnp(:, ia)
            cn2m%u_at(ia, c, :) = -utemp(:) + xpnp(:, ia)
            cn2m%u_ij(ia, c, :) = -utemp(:)    ! /cn2m%r_ij(ia,c)

#else
                  config_real(iconf)%type_neigh(ia, c) = config_real(iconf)%itype(ja)
                  config_real(iconf)%kind_neigh(ia, c) = ja
                  config_real(iconf)%incell(ia, c) = .true.
                  config_real(iconf)%r_ij(ia, c) = dsqrt(r2)
                  config_real(iconf)%u_per(ia, c, :) = -utemp(:) - xpnp(:, ja) + xpnp(:, ia)
                  config_real(iconf)%u_at(ia, c, :) = -utemp(:) + xpnp(:, ia)
                  config_real(iconf)%u_ij(ia, c, :) = -utemp(:)    ! /config_real(iconf)%r_ij(ia,c)
#endif
               end do
            end do
         end do
      end do
      if (c > imm_neigh) then
         call log_critical("neighbours_small: atom "//vtoa(ia)//" has "//vtoa(c)//" neighbours > imm_neigh="//vtoa(imm_neigh))
         call log_critical("  r_cut="//vtoa(sqrt(r_cut2)))
#if(MLD_NDM)
#else
         call log_critical("  config file: "//trim(config_real(iconf)%filename))
#endif
         stop 'neighbours_small: c > imm_neigh. Increase imm_neigh in ml_in_ndm_module.F90'
      end if
#if(MLD_NDM)
    cn2m%n_neigh(ia) = c
#else
      config_real(iconf)%n_neigh(ia) = c
#endif
   end do

   deallocate (xpnp)
   return

   _MLD_END_

end subroutine neighbours_small


subroutine calc_volume(a1, a2, a3, volume)
   implicit none
   real(kind(0.d0)), dimension(3), intent(in)   :: a1, a2, a3
   real(kind(0.d0)), intent(out)    :: volume
   integer  :: i, j, k, l, s, iperm
   !-----------------------------------------------
   volume = 0.0
   i = 1
   j = 2
   k = 3
   s = 1.D0
   do iperm = 1, 3
      volume = volume + s*a1(i)*a2(j)*a3(k)
      l = i
      i = j
      j = k
      k = l
   end do
   i = 2
   j = 1
   k = 3
   s = -s
   do while (s < 0.D0)
      do iperm = 1, 3
         volume = volume + s*a1(i)*a2(j)*a3(k)
         l = i
         i = j
         j = k
         k = l
      end do
      i = 2
      j = 1
      k = 3
      s = -s
   end do
   volume = dabs(volume)

end subroutine calc_volume



subroutine deallocate_real_config(ifile)
   use ml_in_ndm_module, only: ml_type, ml_type_descriptors
   use derived_types, only: config_real
   implicit none
   integer, intent(in)  :: ifile


   if (allocated(config_real(ifile)%itype)) deallocate (config_real(ifile)%itype)
   if (allocated(config_real(ifile)%pos_cart)) deallocate (config_real(ifile)%pos_cart)
   if (allocated(config_real(ifile)%pos_crst)) deallocate (config_real(ifile)%pos_crst)
   if (allocated(config_real(ifile)%force)) deallocate (config_real(ifile)%force)
   if (allocated(config_real(ifile)%atomic_spin)) deallocate (config_real(ifile)%atomic_spin)


   if (allocated(config_real(ifile)%mass_per_type)) deallocate (config_real(ifile)%mass_per_type)
   if (allocated(config_real(ifile)%weight_per_type)) deallocate (config_real(ifile)%weight_per_type)
   if (allocated(config_real(ifile)%Z_per_type)) deallocate (config_real(ifile)%Z_per_type)
   if (allocated(config_real(ifile)%fix_type_poscar_to_periodic)) deallocate (config_real(ifile)%fix_type_poscar_to_periodic)


   if (allocated(config_real(ifile)%proc_atom)) deallocate (config_real(ifile)%proc_atom)
   if (ml_type == ml_type_descriptors) then
      !if (allocated(config_real(ifile)%stress)) deallocate(config_real(ifile)%stress)
      if (allocated(config_real(ifile)%type_neigh)) deallocate (config_real(ifile)%type_neigh)
      if (allocated(config_real(ifile)%kind_neigh)) deallocate (config_real(ifile)%kind_neigh)
      if (allocated(config_real(ifile)%n_neigh)) deallocate (config_real(ifile)%n_neigh)
      if (allocated(config_real(ifile)%r_ij)) deallocate (config_real(ifile)%r_ij)
      if (allocated(config_real(ifile)%u_ij)) deallocate (config_real(ifile)%u_ij)
      if (allocated(config_real(ifile)%u_per)) deallocate (config_real(ifile)%u_per)
   end if

end subroutine


subroutine compute_cos_and_fcut_ja(desc_forces_local)
   use module_neigh_local, only: r_cut
   use module_neigh_local, only: max_neigh_local, tmp_dxp, r_central, tmpcos_dxp, &
      d_r_central, d_r_fcut, r_fcut
   implicit none
   logical, intent(in)  :: desc_forces_local
   integer  :: ii

   ! cosine dx/dr .....
   if (allocated(tmpcos_dxp)) deallocate (tmpcos_dxp); allocate (tmpcos_dxp(3, max_neigh_local))
   do ii = 1, 3
      tmpcos_dxp(ii, 1:max_neigh_local) = tmp_dxp(ii, 1:max_neigh_local)/r_central(1:max_neigh_local)
   end do

   !fcut and derivatives .................
   if (allocated(r_fcut)) deallocate (r_fcut); allocate (r_fcut(max_neigh_local))
   r_fcut(1:max_neigh_local) = ((r_central(1:max_neigh_local)/r_cut)**2 - 1)**2
   if (desc_forces_local) then
      if (allocated(d_r_fcut)) deallocate (d_r_fcut); allocate (d_r_fcut(3, max_neigh_local))
      if (allocated(d_r_central)) deallocate (d_r_central); allocate (d_r_central(3, max_neigh_local))
   end if
   if (desc_forces_local) then
      do ii = 1, 3
         d_r_fcut(ii, 1:max_neigh_local) = tmpcos_dxp(ii, 1:max_neigh_local)* &
            4.d0*((r_central(1:max_neigh_local)/r_cut)**2 - 1)*r_central(1:max_neigh_local)/r_cut**2
      end do
   end if

end subroutine compute_cos_and_fcut_ja


subroutine convert_A2cm(i, xa, xpos)
  ! i= 1 convert positions from A  -> cm
  ! i=-1 convert positions from cm -> A

#ifdef MLD_NDM
  use gen_com_m, only: A2cm
  use gen_com_m_ml, only: at, bg
  use tab_imm_m_ml, only: xp
#else
  use ondm_gen_com_m, only: A2cm, at, bg
  use ondm_tab_imm_m, only: xp
#endif
  use ml_in_ndm_module, ONLY: rangml

  implicit none

  integer, intent(in)  :: i
  logical, intent(in)  :: xa, xpos

  if (i == 1) then        ! A2cm
    if (xa) then
      at(:, :) = at(:, :)*A2cm
      bg(:, :) = bg(:, :)/A2cm
    end if
    if (xpos) xp(:, :) = xp(:, :)*A2cm
  elseif (i == -1) then   ! cm2A
    if (xa) then
      at(:, :) = at(:, :)/A2cm
      bg(:, :) = bg(:, :)*A2cm
    end if
    if (xpos) xp(:, :) = xp(:, :)/A2cm
  else
    if (rangml == 0) write (*, *) "Wrong input : input = 1 or -1"
    stop "fatal in convert_A2cm"
  end if
end subroutine convert_A2cm
