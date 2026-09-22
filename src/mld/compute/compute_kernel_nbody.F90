! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   MiLaDy (Machine Learning Dynamics) 
! HND X   
! HND X
! HND X   Milady was written and designed by Mihai-Cosmin Marinica, Alexandra M. Goryaeva 
! HND X   Copyright 2015-2025.
! HND X
! HND X   Portions of MiLady were written by Wesley UnnToc, Clovis Lapointe, Anruo Zhong, 
! HND X   Alexandre Dezaphie, Jacopo Baima, Anida Khizar, Christian van Wambeke 
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

module grid_and_permutations_nbody
  use iso_fortran_env,  dp => real64
  use mld_logger
  implicit none 
  private 


  type type_tuples_type
    integer, dimension(:), allocatable :: types_db   
  end type type_tuples_type 

  type type_zpoints
    integer :: idx_tuple  
    integer :: idx_rtheta 
  end type type_zpoints

  type type_zpoins_by_tuple
    ! the zpoints put into a matrix by unique tuple dim(zpoint-body) x no_of_zpoints_for_that_unique_tuple
    real(dp), allocatable, dimension(:,:)  :: zz
    integer, allocatable,     dimension(:) :: dim_rtheta_grid
    integer :: idx_min, idx_max          
  end type type_zpoins_by_tuple

  type type_grid_body
    integer :: nu, dim_qnu   
    logical :: automatic
    real(dp) :: r_cut_in, r_cut_out 
    
    !----> tools to handle the sparse points 
    ! all zpoints are stored in a object of size dim(zpoints) computed in init_zpoints 
    type(type_zpoints), dimension(:), allocatable :: zpoints
    ! zpoints_by_tuple is an object which handle the  matrix of zpoints for that tuple 
    ! zpoints_by_typle(zpoint(s)%idx_tuple)%zz(:,zpoint(s)%idx_rtheta))
    type(type_zpoins_by_tuple), dimension(:), allocatable :: zpoints_by_tuple
    !----> end of tools to handle the sparse points

    !----> tools to handle the chemical space by tuples. 
    type(type_tuples_type), dimension(:), allocatable :: tuples 
    integer :: no_unique_tuples
    ! full tuples: for body \nu and no_of_elements is a matrix id_of_tuples \times \nu. 
    ! Where id of tuples is order of (no_of_elements)^\nu   
    integer, dimension(:,:),  allocatable :: mu_space_ini
    ! unique_classes( \nu , id_of_unique_tuples) are the unique classes for the above tuples
    integer, dimension(:,:), allocatable  :: unique_classes  ! unique_classes(dim_permutation, id_tuples)
    ! map_into_unique_classes map the id_of_tuples into  id_of_unique_tuples 
    integer,   dimension(:),  allocatable :: map_into_unique_classes
    ! map_into_tuples  map from  id_of_unique_tuples into  a member of class as id_of_tuples
    integer,   dimension(:),  allocatable :: map_into_tuples
    integer,   dimension(:),  allocatable :: map_hash_into_class
    ! ----> end of tools to handle the chemical space by tuples.

    ! delta_class is a matrix of size unique_classes x unique_classes
    integer, dimension(:,:),  allocatable :: delta_class
    contains 
    procedure :: init_grid_nbody
    procedure :: init_tuples_nbody
    procedure :: init_zpoints

  end type type_grid_body

  type(type_grid_body), dimension(:), allocatable :: grid_nbody 
  public :: grid_nbody 

  contains 
    subroutine init_grid_nbody(this, nu, r_cut_in, r_cut) 
      use module_units, only: two_pi
      implicit none 
      class(type_grid_body), intent(inout) :: this

      integer, intent(in) :: nu  
      real(dp) :: r_cut_in, r_cut
      
      this%nu = nu
      if (this%nu == 2) then 
        this%dim_qnu = 1
      else 
        this%dim_qnu = 3*this%nu - 6
      end if 

      this%r_cut_in = r_cut_in
      this%r_cut_out = r_cut

    end subroutine init_grid_nbody

    subroutine init_zpoints(this, automatic, char_dim_grid)
      implicit none 
      class(type_grid_body), intent(inout) :: this
      ! to generate the grid automatically ... or not.  
      logical, intent(in) :: automatic
      character(len=*), intent(in) :: char_dim_grid
      integer :: no_tuples, is, itmp , jj, itest, ii  

      this%automatic = automatic
      no_tuples = size(this%unique_classes,2)

      if (this%automatic) then
        ! here we will define zpoints_by_tuple
        if (allocated(this%zpoints_by_tuple)) deallocate(this%zpoints_by_tuple)
        allocate(this%zpoints_by_tuple(no_tuples))
        do is = 1, no_tuples 
          call generate_zpoints_by_grid(this, char_dim_grid, this%zpoints_by_tuple(is)%dim_rtheta_grid, this%zpoints_by_tuple(is)%zz)
        end do   
      else 
        call generate_zpoints_by_selection(this) 
        call log_critical("grid generation by selection ... not implemented yet !!!!!!!!!!")
        stop 'error in the grid_generation by selection'
      end if

      ! some test  ..............
      do is = 1, no_tuples 
         itest = 1 
         do ii = 1, size(this%zpoints_by_tuple(is)%dim_rtheta_grid)
           itest = itest * this%zpoints_by_tuple(is)%dim_rtheta_grid(ii)
         end do  
         if (size(this%zpoints_by_tuple(is)%zz,2) /= itest )   then 
           call log_critical("The size of grid is not the same as the numbet of zpoints by tuple. Fatal ! ")
           call log_critical("The number of zpoints for the first tuple is "// vtoa(size(this%zpoints_by_tuple(is)%zz,2)))
           call log_critical("The size of grid is ........................ "// vtoa(itest))
           stop 'error in the grid generation init_zpoints'
         end if
      end do


      itmp = 0 
      do is = 1, no_tuples 
        this%zpoints_by_tuple(is)%idx_min = itmp + 1
        itmp = itmp + size(this%zpoints_by_tuple(is)%zz,2)
        this%zpoints_by_tuple(is)%idx_max = itmp
      end do


      if (allocated(this%zpoints)) deallocate(this%zpoints)
      allocate(this%zpoints(itmp))
      itmp = 0
      do is = 1, no_tuples 
        do jj = 1, size(this%zpoints_by_tuple(is)%zz,2)
          itmp = itmp + 1
          this%zpoints(itmp)%idx_tuple = is
          this%zpoints(itmp)%idx_rtheta = jj
        end do 
      end do

    end subroutine init_zpoints

    subroutine generate_zpoints_by_grid(this, char_dim_grid, dim_rgrid, radang) 
      use mld_string, only: count_integers_in_string
      use module_units, only: two_pi
      use module_body_desc, only: bond_dist_transform, bond_beta, bond_dist_ann 
      implicit none 
      class(type_grid_body), intent(inout) :: this
      real(dp), allocatable, dimension(:,:), intent(inout) :: radang
      integer, allocatable,  dimension(:) :: dim_rgrid  
      character(len=*), intent(in)   :: char_dim_grid  
      real(dp), dimension(:,:), allocatable :: limits
      integer :: itest, dim_all_rgrid   
      real(dp) :: x_min, x_max, r_min, r_max 
      ! to generate the zpoints for the grid ...
      
      
      if (allocated(dim_rgrid)) deallocate(dim_rgrid) ;  allocate(dim_rgrid(this%dim_qnu))
      itest = count_integers_in_string(char_dim_grid) 
      if (itest /= this%dim_qnu) then 
        call log_critical("The number of integers in the string dim_grid_nbody is not equal to the body dimension") 
        call log_critical(" body dimension is "// vtoa(this%dim_qnu))
        call log_critical(" the number of integers  is "// vtoa(itest))
        call log_critical(" dim_grid_nbody is  " // trim(char_dim_grid))
        stop 'error in the grid generation'
      end if
      
      read(char_dim_grid,*) dim_rgrid
      dim_all_rgrid = product(dim_rgrid(:))
      !$! write(*,*) 'debug --> dim_rgrid', dim_rgrid
      !$! write(*,*) 'debug --> dim_all_rgrid', dim_all_rgrid

      if (allocated(radang)) deallocate(radang)
      allocate(radang(this%dim_qnu, dim_all_rgrid))
      if (allocated(limits)) deallocate(limits)
      allocate(limits(2, this%dim_qnu))

      r_min = this%r_cut_in
      r_max = this%r_cut_out
      
      ! Apply coordinate transformation to limits
      if (bond_dist_transform == 1) then ! Linear
         x_min = r_min
         x_max = r_max
      else if (bond_dist_transform == 2) then ! Exp: exp(-beta*r) is decreasing
         x_min = exp(-bond_beta * r_max) 
         x_max = exp(-bond_beta * r_min)
      else if (bond_dist_transform == 3) then ! Inverse
         if (r_max > 1.d-12) x_min = (bond_dist_ann**bond_beta) / (r_max**bond_beta)
         if (r_min > 1.d-12) x_max = (bond_dist_ann**bond_beta) / (r_min**bond_beta)
      endif

      if (this%nu == 3) then 
        ! limits for x1 (q1 = x1 + x2)
        limits(1,1) = 2.0_dp * x_min
        limits(2,1) = 2.0_dp * x_max
        !limits for x2 (q2 = x1 * x2)
        limits(1,2) = x_min * x_min
        limits(2,2) = x_max * x_max
        !limits for x3 (theta)
        limits(1,3) = 0.0_dp
        limits(2,3) = two_pi
      else if (this%nu == 2) then 
        ! limits for x1
        limits(1,1) = x_min
        limits(2,1) = x_max
      else  
        call log_critical("Automatic grid generation ... not implemented yet for nu /= 3")
        stop 
      end if
     
      !d! write(*,*) 'debug --> limits', size(limits,1), size(limits,2)
      !d! write(*,*) 'debug --> limits', limits 
      !d! write(*,*) dim_rgrid

      call generate_linear_grid(dim_rgrid, limits, radang)
      !d! write(*,*) 'debug --> dim(radang) 1 and 2', size(radang,1), size(radang,2)
      !d! write(*,*) radang
      !d! stop 
 
    end subroutine generate_zpoints_by_grid 


    subroutine generate_zpoints_by_selection(this)
      class(type_grid_body), intent(inout) :: this  
      ! some complicated selection ...CUR, MAHA, RANDOM, CLUSTERING ... 
      ! TO_DO_ZETABODY 
    end subroutine generate_zpoints_by_selection

    subroutine generate_linear_grid(dim_rgrid, limits, radang)
      use mesh_grid, only: linear_grid
      implicit none 
      integer, dimension(:), intent(in) :: dim_rgrid
      real(dp), dimension(:,:), intent(out) :: radang
      real(dp), dimension(:,:), intent(in) :: limits
      !
      type type_one_dim_grid
        real(dp), dimension(:), allocatable :: grid
      end type type_one_dim_grid
      ! 
      type(type_one_dim_grid), dimension(:), allocatable :: one_dim_grid  
      integer :: ii, jj, kk, icount 
      
      allocate(one_dim_grid(size(dim_rgrid)))
      do ii = 1, size(dim_rgrid)
        allocate(one_dim_grid(ii)%grid(dim_rgrid(ii)))
        call linear_grid(dim_rgrid(ii), limits(1,ii), limits(2,ii), one_dim_grid(ii)%grid)
      end do
      
      ! generate the grid ...

      if (.not.((size(dim_rgrid) == 1) .or. (size(dim_rgrid) == 3))) then
         call log_critical("zeta-grid is not defined for this dimension")
         stop "zeta-grid error 1or3" 
      end if 
      if (size(dim_rgrid) == 1) then 
         do ii = 1, dim_rgrid(1)
           radang(1, ii) = one_dim_grid(1)%grid(ii)
         end do  
      end if 
      if (size(dim_rgrid) == 3) then 
        icount = 0 
        do ii = 1, dim_rgrid(1)
          do jj = 1, dim_rgrid(2)
            do kk = 1, dim_rgrid(3)
              icount = icount + 1
              !$! radang(1, (ii-1)*dim_rgrid(2)*dim_rgrid(3) + (jj-1)*dim_rgrid(3) + kk) = one_dim_grid(1)%grid(ii)
              !$! radang(2, (ii-1)*dim_rgrid(2)*dim_rgrid(3) + (jj-1)*dim_rgrid(3) + kk) = one_dim_grid(2)%grid(jj)
              !$! radang(3, (ii-1)*dim_rgrid(2)*dim_rgrid(3) + (jj-1)*dim_rgrid(3) + kk) = one_dim_grid(3)%grid(kk)
              radang(1, icount) = one_dim_grid(1)%grid(ii)
              radang(2, icount) = one_dim_grid(2)%grid(jj)
              radang(3, icount) = one_dim_grid(3)%grid(kk)
            end do 
          end do 
        end do 
      end if
      

    end subroutine generate_linear_grid


    subroutine init_tuples_nbody(this, no_of_elements)
      use module_mu_space_first, only: generate_mu_space_ini, find_permutations_for_mu !, base_mu 
      use module_l_space_L0, only: classify_tuples_for_l
      use math, only: get_max_hash_key_with_basis, gen_hash_key_with_basis
      class(type_grid_body), intent(inout) :: this
      integer, dimension(:,:), allocatable :: mu_space_ini 
      integer :: num_classes

      integer, dimension(:,:), allocatable :: unique_classes  ! unique_classes(dim_permutation, id_tuples)
      integer, dimension(:), allocatable :: map_into_unique_classes, map_into_tuples !  hash_values(num_perm)

      integer, intent(inout) :: no_of_elements
      integer :: dim_mu, ii 

      integer, dimension(:,:), allocatable :: tmp_unique_classes  ! unique_classes(dim_permutation, id_tuples)
      integer, dimension(:), allocatable :: tmp_map_into_unique_classes, tmp_map_into_tuples !  hash_values(num_perm)
      integer, dimension(:,:), allocatable :: tmp_mu_space_ini
      integer :: new_size_1, new_size_2, icount, jj, nn, nu_val, max_hash_key, size_of_basis, hash_key 
      integer, dimension(:) , allocatable :: ivec 

      ! build a mu space ini for nu-1 with no_of_elements
      call generate_mu_space_ini(this%nu-1, no_of_elements, mu_space_ini, dim_mu)
      call classify_tuples_for_l (mu_space_ini, num_classes,  map_into_unique_classes, map_into_tuples, unique_classes)

      !!$! from mu_space_ini we get only the unique classes ... 
      !$! write(*,*) size(mu_space_ini,1), size(mu_space_ini,2)
      !$! do ii = 1,  size(mu_space_ini,1)
      !$!      write(*,*) ii, " (  ", mu_space_ini(ii,1) , mu_space_ini(ii,2), "  )  "
      !$! end do 

      !$! do ii = 1, size(unique_classes,2)
      !$!     write(*,*) ii, "   ii   : ", unique_classes(:,ii)
      !$! end do   !$! 
      !$! do ii = 1, size(map_into_unique_classes)
      !$!     write(*,*) ii, "   map_to_cla   : ", map_into_unique_classes(ii)
      !$! end do   !$! 
      !$! do ii = 1, size(map_into_tuples)
      !$!     write(*,*) ii, "   map_into_tuples   : ", map_into_tuples(ii)
      !$! end do 
      ! prepare tuples ...
      new_size_2 = size(unique_classes,2) * no_of_elements
      new_size_1 = size(unique_classes,1) + 1 

      if (allocated(tmp_unique_classes)) deallocate(tmp_unique_classes)  
      allocate(tmp_unique_classes(new_size_1, new_size_2))

      if (allocated(tmp_map_into_tuples)) deallocate(tmp_map_into_tuples)
      allocate(tmp_map_into_tuples(new_size_2))

      icount = 0 
      do jj = 1, no_of_elements
        do ii = 1, size(unique_classes,2)
            !write(*,*) ii, "   ii   : ", unique_classes(:,ii)
            icount = icount +1 
            do nn = 1, this%nu-1
                tmp_unique_classes(nn+1, icount) = unique_classes(nn,ii) 
            end do 
            tmp_unique_classes(1, icount) = jj
            !tmp_map_into_tuples(icount) = (map_into_tuples(ii)-1)*no_of_elements + jj
            tmp_map_into_tuples(icount) = map_into_tuples(ii) + (jj-1)*size(mu_space_ini,1)
        end do 
      end do 

      new_size_2 = size(mu_space_ini,2) + 1
      new_size_1 = size(mu_space_ini,1) * no_of_elements
      if (allocated(tmp_mu_space_ini)) deallocate(tmp_mu_space_ini)
      allocate(tmp_mu_space_ini(new_size_1, new_size_2))

      if (allocated(tmp_map_into_unique_classes)) deallocate(tmp_map_into_unique_classes)
      allocate(tmp_map_into_unique_classes(new_size_1))
      icount = 0 
      do jj = 1, no_of_elements
        do ii = 1, size(mu_space_ini,1)
            icount = icount +1 
            do nn = 1, this%nu-1
                tmp_mu_space_ini(icount, nn+1) = mu_space_ini(ii,nn) 
            end do 
            tmp_mu_space_ini(icount,1) = jj
            !tmp_map_into_unique_classes(icount) = (map_into_unique_classes(ii)-1)*no_of_elements + jj
            tmp_map_into_unique_classes(icount) = (map_into_unique_classes(ii)) + (jj-1)*size(unique_classes,2)
        end do 
      end do

      deallocate(map_into_unique_classes, map_into_tuples, unique_classes, mu_space_ini)
      allocate(this%map_into_unique_classes(size(tmp_map_into_unique_classes)), &
               this%map_into_tuples(size(tmp_map_into_tuples)), &
               this%unique_classes(size(tmp_unique_classes,1), size(tmp_unique_classes,2)), &
               this%mu_space_ini(size(tmp_mu_space_ini,1), size(tmp_mu_space_ini,2)))
      this%map_into_unique_classes = tmp_map_into_unique_classes
      this%map_into_tuples = tmp_map_into_tuples
      this%unique_classes = tmp_unique_classes
      this%mu_space_ini = tmp_mu_space_ini
      deallocate(tmp_map_into_unique_classes, tmp_map_into_tuples, tmp_unique_classes, tmp_mu_space_ini)

      !$! if (size(this%mu_space_ini,2) == 3) then  
      !$! do ii = 1,  size(this%mu_space_ini,1)
      !$!      write(*,*) ii, " ( x ini  ", this%mu_space_ini(ii,1) , this%mu_space_ini(ii,2), this%mu_space_ini(ii,3), "  )  "
      !$! end do 
      !$! do ii = 1, size(this%unique_classes,2)
      !$!     write(*,*) ii, " x ii   : ", this%unique_classes(:,ii)
      !$! end do 
      !$! do ii = 1, size(this%map_into_unique_classes)
      !$!     write(*,*) ii, " x map_to_cla   : ", this%map_into_unique_classes(ii)
      !$! end do 
      !$! do ii = 1, size(this%map_into_tuples)
      !$!     write(*,*) ii, " x map_into_tuples   : ", this%map_into_tuples(ii)
      !$! end do 
      !$! ! prepare tuples ...
      !$! end if 

      
      size_of_basis = no_of_elements  
      !$! write(*,*) 'debug --> size_of_basis', size_of_basis
      nu_val = size(this%mu_space_ini,2)

      allocate(ivec(nu_val))

      !$! write(*,*) 'debug --> nu_val', nu_val
      call get_max_hash_key_with_basis(nu_val, size_of_basis, max_hash_key)
      !$! write(*,*) 'debug --> max_hash_key', size_of_basis, max_hash_key
      if (allocated(this%map_hash_into_class))  deallocate(this%map_hash_into_class)
      allocate(this%map_hash_into_class(max_hash_key))
      this%map_hash_into_class(:) = -1 
      do ii = 1, size(this%mu_space_ini,1)
        ivec= this%mu_space_ini(ii,:)
        !$! write(*,*) ii, 'debug --> ivec', ivec
        call gen_hash_key_with_basis( nu_val, size_of_basis, ivec, hash_key)
        !$! write(*,*) ii, 'debug --> hash_key',  hash_key
        this%map_hash_into_class(hash_key) = this%map_into_unique_classes(ii)
      end do 

      this%no_unique_tuples = size(this%unique_classes,2)
      

      ! produce delta_class 
      if (allocated(this%delta_class)) deallocate(this%delta_class)
      allocate(this%delta_class(size(this%unique_classes,2), size(this%unique_classes,2))) 
      do ii  = 1, size(this%unique_classes,2)
        do jj = 1, size(this%unique_classes,2)
          if (ii == jj ) then 
            this%delta_class(ii,jj) = 1
          else 
            this%delta_class(ii,jj) = 0
          end if 
        end do 
      end do 

    end subroutine init_tuples_nbody 
  
end module grid_and_permutations_nbody 



subroutine compute_kernel_k2b(i_start_at, i_final_at, iconf)
  use module_kind_variables, ONLY: kind_double
  use ml_in_ndm_module, only: desc_forces,  rangml, imm_neigh
  use derived_types, only: config_desc, config_real
#ifdef MLD_NDM
  use tab_imm_m_ml, ONLY: xp
  use gen_com_m, only: lperiod
  use gen_com_m_ml, only: imm,at,bg
#else
  use ondm_gen_com_m, only: lperiod, imm
  use ondm_tab_imm_m, ONLY: iwmax2, xp
#endif
  use module_neigh_local, only : max_neigh_local, i_type, i_type_db, & 
                                 tmp_dxp, tmpcos_dxp,  i_central, r_central, tmp_xp, type_fcut, &
                                 r_fcut, d_r_fcut, &
                                 build_local_neighbours_ja, build_cos_and_fcut_ja, reallocate_neigh_ja
                                  
  use module_kernel_2b, only: dim_kernel_2b, np_radial_2b,  delta_2b, sigma_2b, zr_2b, zr_fcut_2b, & 
                              delta_type, zpoints_2b, tnn_2b, tuu_2b, tkk_2b, r_cut_2b, r_cut_width_2b
  use time_check_general, only: debug_time, MY_MPI_WTIME
  use math, only: my_exp 
  use mld_logger
#ifdef MLD_NDM
  use notperiod_mod

#else
  use ondm_transform_coord, only: ondm_notperiod
#endif
  use mld_subworld
  implicit none
  integer, intent(in)  :: i_start_at, i_final_at, iconf

  integer, dimension(imm)   :: d_n_neigh
  integer, dimension(imm, imm_neigh)   :: d_kind_neigh
  real(kind_double), dimension(:, :), allocatable  :: xpnp
  
  integer :: ja,  ir, ik, iw2, ja_atom, max_neigh, ia_n 
  integer :: type_db_ja, type_db_ia, type_z1, type_z2
  integer :: type_fact 
  real(kind_double) :: r_ji, tmp_kk
  logical  :: desc_forces_local, small 
  real(kind_double), dimension(:), allocatable :: energy_k2b
  real(kind_double), dimension(:,:), allocatable :: tmp_energy, temp
  real(kind_double), dimension(:,:,:), allocatable :: tmp_force
  real(kind_double) :: ttt(1:3), t00, t11, t22, t33 
  real(kind_double), external  :: sin, cos
  
  _NAMECURRENT_("compute_kernel_k2b")

  max_neigh = imm_neigh
  if ((i_start_at == 0) .and. (i_final_at == 0)) then
    config_desc(iconf)%energy_k2b(:, :) = 0.d0
    d_n_neigh(:) = 0
    d_kind_neigh(:, :) = 0
    return
  end if
  desc_forces_local = desc_forces .and. (config_real(iconf)%has_force .or. config_real(iconf)%has_stress)
  !cos config_desc(iconf)%energy_k2b(:, :) = 0.d0
  !cos if (desc_forces_local) config_desc(iconf)%force_k2b(:, :, :, :) = 0.d0

  small = config_real(iconf)%small

  ALLOCATE (xpnp(3, imm))
  if (lperiod) then
    xpnp(:, :) = xp(:, :)
  else
    call ondm_notperiod(xp, xpnp)
  end if

  d_n_neigh(:) = 0
  d_kind_neigh(:, :) = 0
  config_desc(iconf)%energy_k2b(:, :) = 0.d0
  if (desc_forces_local) config_desc(iconf)%force_k2b(:, :, :, :) = 0.d0

  if (i_start_at == 1) iw2 = 0
  if (i_start_at > 1)  iw2 = iwmax2(i_start_at - 1)

   
  do ja = i_start_at, i_final_at

    ja_atom = ja
    type_db_ja = config_real(iconf)%itype_db(ja)


    if (debug_time) t00 = MY_MPI_WTIME()
    call build_local_neighbours_ja( iconf, ja, imm, xpnp, r_cut_2b, d_n_neigh, d_kind_neigh, &
                                    r_central, i_type, i_type_db, i_central, tmp_dxp, tmp_xp, &
                                    max_neigh_local, iw2)



    if (max_neigh_local == 0) cycle
    if (max_neigh_local > max_neigh) then
        write (6, *) 'ML: Fatal error in compute_kernel_k2b, rank=', rangml, ' ja=', ja
        write (6, *) '  max_neigh_local=', max_neigh_local, ' > max_neigh=', max_neigh
        write (6, *) '  r_cut_2b=', r_cut_2b, ' small=', small
        write (6, *) '  n_neigh(ja) from config_real=', config_real(iconf)%n_neigh(ja)
        write (6, *) '  config file: ', trim(config_real(iconf)%filename)
      stop 'compute_kernel_k2b: max_neigh_local > imm_neigh. Increase imm_neigh in ml_in_ndm_module.F90'
    end if
    call reallocate_neigh_ja (max_neigh_local, ja_atom, i_central, i_type, i_type_db, r_central, tmp_dxp)
    if (debug_time) then 
      t11 =MY_MPI_WTIME()
      tnn_2b= tnn_2b + t11 - t00    
    end if 
  
    call build_cos_and_fcut_ja(type_fcut, desc_forces_local, &
                 r_cut_2b, r_cut_width_2b, max_neigh_local, r_central, tmp_dxp, &
                 tmpcos_dxp, r_fcut,  d_r_fcut)

    if (allocated(tmp_energy)) deallocate(tmp_energy)  
    allocate(tmp_energy(np_radial_2b, max_neigh_local))

    if (desc_forces_local) then 
      if (allocated(tmp_force)) deallocate(tmp_force)  
      allocate(tmp_force(np_radial_2b, max_neigh_local,3))
    end if 
 
    do ia_n = 1, max_neigh_local 
      r_ji = r_central(ia_n)
      do ir = 1, np_radial_2b 
        tmp_kk = delta_2b**2*zr_fcut_2b(ir) * exp(-0.5d0 * (r_ji - zr_2b(ir) )**2 / sigma_2b**2)
        tmp_energy(ir, ia_n) =  tmp_kk*r_fcut(ia_n)
        if (desc_forces_local) then
          tmp_force(ir, ia_n, 1:3) = tmp_kk*(d_r_fcut(1:3,ia_n) - (r_ji-zr_2b(ir))/sigma_2b**2*tmpcos_dxp(1:3, ia_n) *r_fcut(ia_n))
        end if   
      end do
    end do  !ia_n
        

    if (debug_time) then 
      t22 =MY_MPI_WTIME()
      tuu_2b= tuu_2b + t22 - t11    
    end if 


    if (allocated(energy_k2b)) deallocate(energy_k2b) ; allocate(energy_k2b(dim_kernel_2b))
    energy_k2b(1:dim_kernel_2b) = 0.d0 
    if (desc_forces_local) then
      if (allocated(temp)) deallocate(temp) 
      allocate(temp(dim_kernel_2b,3))
      temp(:,:) = 0.d0   
    end if   

    do ia_n = 1, max_neigh_local
      type_db_ia =  config_real(iconf)%itype_db(i_central(ia_n))
      !$! if (abs(type_db_ia) > 1 ) then 
      !$!   write(*,*) 'rrr1', rangml, id_subworld, type_db_ia, iconf, ia_n, i_central(ia_n)
      !$!   do ii = 1, size(config_real(iconf)%itype_db)
      !$!     write(*,*) 'rrr2', rangml, id_subworld, ii, config_real(iconf)%itype_db(ii), i_central(ia_n), i_size_on_proc(rangml)
      !$!   end do
      !$! end if 
      do ik = 1, dim_kernel_2b
        type_z1 = zpoints_2b(ik)%type1_z
        type_z2 = zpoints_2b(ik)%type2_z
        type_fact = delta_type(type_db_ia, type_z1)*delta_type(type_db_ja, type_z2) + delta_type(type_db_ia, type_z2)*delta_type(type_db_ja, type_z1)
        if (type_fact >= 1.d-15) then 
          ir = zpoints_2b(ik)%irad 
          energy_k2b(ik) = energy_k2b(ik) + type_fact* tmp_energy(ir, ia_n)
          if (desc_forces_local) then 
            ttt(1:3) = type_fact*tmp_force(ir, ia_n, 1:3)
            temp(ik, 1:3) = temp(ik,1:3) - ttt(1:3)
            config_desc(iconf)%force_k2b(ik, ja, ia_n, 1:3) = ttt (1:3)
          end if 
        end if 
      end do   ! ik 
    end do    ! ia_n second ... 

!$!     do ir = 1, np_radial_2b
!$!       write(*,'("LOLOL ", 1i5, 3f20.10)') type_db_ja,  tz1, tz2, ttt0 
!$!   end do 



    config_desc(iconf)%energy_k2b(1:dim_kernel_2b, ja) = energy_k2b(1:dim_kernel_2b) 
    if (desc_forces_local) then 
      do ik = 1, dim_kernel_2b
        config_desc(iconf)%force_k2b(ik, ja, 0, 1:3) = temp(ik,1:3)
      end do      
    end if

    if (debug_time) then 
      t33 =MY_MPI_WTIME()
      tkk_2b= tkk_2b + t33 - t22    
    end if 
 
  end do !ja

  _MLD_END_
end subroutine compute_kernel_k2b





subroutine init_kernel_k2b
  use module_kind_variables, only: kind_double
  use module_chemical_species, only: fix_no_of_elements
  use module_kernel_2b, only: np_radial_2b, dim_kernel_2b,  zpoints_2b, & 
                              zr_fcut_2b, zr_2b, delta_type, time_for_desc_kernel_2b, & 
                              tnn_2b, tuu_2b, tkk_2b
  use mesh_grid, only: linear_grid
  use module_neigh_local, only : type_fcut, r_cut, r_cut_width, build_cos_and_fcut_ja 
  use mld_logger
  implicit none 
  integer :: ii, jj , icount
  logical, parameter :: desc_forces_local = .false.  
  real(kind_double), dimension(:,:), allocatable :: tmp_dxp, d_zr_fcut  
  real(kind_double), dimension(:,:), allocatable :: tmpcos_dxp
  integer :: ir

  _NAMECURRENT_("init_kernel_k2b")

  dim_kernel_2b = np_radial_2b*fix_no_of_elements*(fix_no_of_elements+1)/2

  if (allocated(delta_type)) deallocate(delta_type) ; allocate(delta_type(fix_no_of_elements, fix_no_of_elements))
  do ii = 1, fix_no_of_elements
    do jj = 1, fix_no_of_elements
      if (ii == jj ) then 
        delta_type(jj, ii) = 1  
      else   
        delta_type(jj, ii) = 0  
      end if
    end do 
  end do     

  if (allocated(zpoints_2b)) deallocate(zpoints_2b) ; allocate(zpoints_2b(dim_kernel_2b))
  icount = 0  
  do ii = 1, fix_no_of_elements
    do ir = 1, np_radial_2b
      icount = icount + 1 
      zpoints_2b(icount)%type1_z = ii 
      zpoints_2b(icount)%type2_z = ii
      zpoints_2b(icount)%irad = ir 
      zpoints_2b(icount)%rr = 0.d0 
    end do  
  end do   
  
  do ii = 1, fix_no_of_elements
    do jj = ii + 1, fix_no_of_elements
      do ir = 1, np_radial_2b
        icount = icount + 1 
        zpoints_2b(icount)%type1_z = ii 
        zpoints_2b(icount)%type2_z = jj
        zpoints_2b(icount)%irad = ir
        zpoints_2b(icount)%rr = 0.d0 
      end do    
    end do    
  end do


  call linear_grid(np_radial_2b, 0.2d0, r_cut, zr_2b) 
  !define number of species pairs ... 
  !debug! write(*,'("type_fcut ..........................", i6, 2f20.10)') type_fcut, r_cut, r_cut_width

  call  build_cos_and_fcut_ja(type_fcut, desc_forces_local, &
    r_cut, r_cut_width, np_radial_2b, zr_2b, tmp_dxp, &
    tmpcos_dxp, zr_fcut_2b, d_zr_fcut)

  time_for_desc_kernel_2b = 0.d0 
  tnn_2b = 0.d0 
  tuu_2b = 0.d0 
  tkk_2b = 0.d0 

  _MLD_END_
end subroutine init_kernel_k2b 


subroutine allocate_kernel_k2b(iconf, i_start_at, i_final_at, dim_kernel_2b, imm, imm_neigh, desc_forces_local)
  use derived_types, only: config_desc 
  use mld_logger
  implicit none 
  integer, intent(in) :: iconf, imm_neigh, imm 
  integer, intent(in) :: i_start_at, i_final_at, dim_kernel_2b
  logical, intent(in) :: desc_forces_local 

  _NAMECURRENT_("allocate_kernel_k2b")

  if (allocated(config_desc(iconf)%energy_k2b)) deallocate (config_desc(iconf)%energy_k2b) 
  allocate (config_desc(iconf)%energy_k2b(dim_kernel_2b, imm))
  if (desc_forces_local) then
    if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
      if (allocated(config_desc(iconf)%force_k2b)) deallocate (config_desc(iconf)%force_k2b) 
      allocate (config_desc(iconf)%force_k2b(dim_kernel_2b, i_start_at:i_final_at, 0:imm_neigh, 1:3))
    end if
  end if
  _MLD_END_
end subroutine allocate_kernel_k2b  

module module_form_kernel_nbody
  implicit none 
  contains 
  subroutine kernel_nbody_exp_list(desc_forces_local, r_ji, ia_n, np_radial_2b, delta_2b, sigma_2b, zr_2b, zr_fcut_2b, r_fcut, d_r_fcut, &
                              tmpcos_dxp, &
                              tmp_energy, tmp_force)
    use module_kind_variables, only: kind_double
    implicit none 
    logical :: desc_forces_local 
    integer, intent(in) :: np_radial_2b , ia_n
    real(kind_double), intent(in) :: delta_2b, sigma_2b, r_ji 
    real(kind_double), dimension(:), allocatable, intent(in) :: zr_2b, zr_fcut_2b, &
                                                   r_fcut
    real(kind_double), dimension(:,:), allocatable, intent(in) :: tmpcos_dxp , d_r_fcut                                                 
    real(kind_double), dimension(:,:), allocatable, intent(inout) :: tmp_energy  
    real(kind_double), dimension(:,:,:), allocatable, intent(inout) :: tmp_force  

    integer :: ir 
    real(kind_double) :: tmp_kk 

      do ir = 1, np_radial_2b 
        tmp_kk = delta_2b**2*zr_fcut_2b(ir) * exp(-0.5d0 * (r_ji - zr_2b(ir) )**2 / sigma_2b**2)
        tmp_energy(ir, ia_n) =  tmp_kk*r_fcut(ia_n)
        if (desc_forces_local) then
          tmp_force(ir, ia_n, 1:3) = tmp_kk*(d_r_fcut(1:3,ia_n) - (r_ji-zr_2b(ir))/sigma_2b**2*tmpcos_dxp(1:3, ia_n) *r_fcut(ia_n))
        end if   
      end do

  end subroutine kernel_nbody_exp_list

  subroutine kernel_nbody_exp_rr(desc_forces_local, r_ji, np_radial_2b, delta_2b, sigma_2b, zr_2b, zr_fcut_2b, r_fcut, d_r_fcut, &
    tmp_energy, tmp_force)
    use module_kind_variables, only: kind_double
    implicit none 
    logical :: desc_forces_local 
    integer, intent(in) :: np_radial_2b
    real(kind_double), intent(in) :: delta_2b, sigma_2b, r_ji, r_fcut, d_r_fcut
    real(kind_double), dimension(:), allocatable, intent(in) :: zr_2b, zr_fcut_2b
    real(kind_double),  dimension(:), allocatable, intent(inout) :: tmp_energy  
    real(kind_double),  dimension(:), allocatable, intent(inout) :: tmp_force  

    integer :: ir 
    real(kind_double) :: tmp_kk 

    do ir = 1, np_radial_2b 
      tmp_kk = delta_2b**2*zr_fcut_2b(ir) * exp(-0.5d0 * (r_ji - zr_2b(ir) )**2 / sigma_2b**2)
      tmp_energy(ir) =  tmp_kk*r_fcut
      if (desc_forces_local) then
        tmp_force(ir) = tmp_kk*(d_r_fcut - (r_ji-zr_2b(ir))/sigma_2b**2*r_fcut)
        end if   
    end do

end subroutine kernel_nbody_exp_rr

subroutine kernel_nbody_exp_rr_second(r_ji, np_radial_2b, delta_2b, sigma_2b, zr_2b, zr_fcut_2b, r_fcut, d_r_fcut, dd_r_fcut, & 
                                      tmp_energy, tmp_force, tmp_second )
  use module_kind_variables, only: kind_double
  implicit none 
  integer, intent(in) :: np_radial_2b
  real(kind_double), intent(in) :: delta_2b, sigma_2b, r_ji, r_fcut, d_r_fcut, dd_r_fcut 
  real(kind_double), dimension(:), allocatable, intent(in) :: zr_2b, zr_fcut_2b
  real(kind_double),  dimension(:), allocatable, intent(inout) :: tmp_energy  
  real(kind_double),  dimension(:), allocatable, intent(inout) :: tmp_force  
  real(kind_double),  dimension(:), allocatable, intent(inout) :: tmp_second   

  integer :: ir 
  real(kind_double) :: tmp_kk 

  do ir = 1, np_radial_2b 
    tmp_kk = delta_2b**2*zr_fcut_2b(ir) * exp(-0.5d0 * (r_ji - zr_2b(ir) )**2 / sigma_2b**2)
    tmp_energy(ir) =  tmp_kk*r_fcut
    tmp_force(ir) = tmp_kk*(d_r_fcut - (r_ji-zr_2b(ir))/sigma_2b**2*r_fcut)
    tmp_second(ir) = tmp_kk*(- r_fcut + (r_ji-zr_2b(ir))**2/sigma_2b**2*r_fcut  -2.d0*(r_ji-zr_2b(ir))*d_r_fcut + sigma_2b**2*dd_r_fcut) /sigma_2b**2
  end do

end subroutine kernel_nbody_exp_rr_second 


subroutine kernel_2b_end  (type_ja, type_ia, desc_forces_local, l_type_fcut, rr, l_r_cut, l_rcut_width, energy, force, second)
  ! estimate the kernel and its radial derivative at the distance rr 
    use module_kind_variables, only: kind_double
    use module_neigh_local, only: fcut_rij, fcut_rij_second
    !use module_form_kernel_nbody, only: kernel_nbody_exp_rr, kernel_nbody_exp_rr_second
    use module_kernel_2b,only: np_radial_2b, delta_2b, sigma_2b, zr_2b, zr_fcut_2b, dim_kernel_2b, &
       zpoints_2b, delta_type
    use module_zbl, only: type_rac_zbl, rac_zbl_exp, rac_zbl_poly
    use snap, only: w_params, dim_xdesc_linear 
    use mld_logger, only: log_debug, log_warning, vtoa, mld_verbose, log_critical
    !use derived_types, only: typ_species_half
    implicit none 
    logical, intent(in) :: desc_forces_local 
    real(kind_double), intent(in) :: rr, l_r_cut, l_rcut_width
    integer, intent(in) :: l_type_fcut , type_ja, type_ia
    real(kind_double), intent(out) :: energy, force, second   
    
    !internal variables ...
    real(kind_double) :: rr_fcut, d_rr_fcut, dd_rr_fcut,  type_fact
    real(kind_double), dimension(:), allocatable :: tmp_energy, tmp_force, tmp_second
    integer :: ik, ir, type_z1, type_z2, ipini
    
    _NAMECURRENT_("estimate_kernel_2b")

    allocate(tmp_energy(np_radial_2b))
    allocate(tmp_force(np_radial_2b))
    allocate(tmp_second(np_radial_2b))

    tmp_second = 0.d0 

    select case(type_rac_zbl)
      case(rac_zbl_exp) 
        call fcut_rij(l_type_fcut, rr, l_r_cut, l_rcut_width, desc_forces_local, rr_fcut, d_rr_fcut)
        call kernel_nbody_exp_rr(desc_forces_local,rr, np_radial_2b, delta_2b, sigma_2b, zr_2b, zr_fcut_2b, rr_fcut,  & 
                                 d_rr_fcut, tmp_energy, tmp_force)
      case(rac_zbl_poly)                           
        call fcut_rij_second(l_type_fcut, rr, l_r_cut, l_rcut_width, rr_fcut, d_rr_fcut, dd_rr_fcut)
        call kernel_nbody_exp_rr_second(rr, np_radial_2b, delta_2b, sigma_2b, zr_2b, zr_fcut_2b, rr_fcut, & 
                               d_rr_fcut, dd_rr_fcut, tmp_energy, tmp_force, tmp_second)    
      case default                                              
        call log_critical('No implementation for this zbl '//vtoa(type_rac_zbl)//' in subroutine '//NAMECURRENT )
    end select 
    energy = 0.d0 
    force = 0.d0 
    second = 0.d0 
    ipini = dim_xdesc_linear

    do ik = 1, dim_kernel_2b
      type_z1 = zpoints_2b(ik)%type1_z
      type_z2 = zpoints_2b(ik)%type2_z
      type_fact = delta_type(type_ia, type_z1)*delta_type(type_ja, type_z2) + delta_type(type_ia, type_z2)*delta_type(type_ja, type_z1)
      if (type_fact >= 1.d-15) then 
        ir = zpoints_2b(ik)%irad 
        energy = energy + type_fact* tmp_energy(ir)*w_params(ipini + ik,1)
        if (desc_forces_local) then 
          force = force +   type_fact*tmp_force(ir)*w_params(ipini + ik,1)
          if (type_rac_zbl == rac_zbl_poly) second = second + type_fact*tmp_second(ir)*w_params(ipini + ik,1)
        end if 
      end if 
    end do
    
    deallocate(tmp_energy)
    deallocate(tmp_force)
    deallocate(tmp_second)
                      
    _MLD_END_

  end subroutine kernel_2b_end 


end module module_form_kernel_nbody 


module module_continuity_k2b_zbl
  implicit none 
  contains 

  subroutine  parameters_two_ends_function_continuity (r1, r2, f1, f1d, f2, f2d, params)
    !
    !  Having the values of two functions f1 and f2 in r1 and r2 , respectively: 
    !  f1(r1) = f1, 
    !  f1'(r1)= f1d,
    !  f2(r2) = f2, 
    !  f2'(r2)= f2d,
    !  provides the 4 parameters params (a,b,c,d) of the function f(x) = exp(a + b*x + c*x**2 + d*x**3) than 
    !  ensure the continuity of the function and the derivatives between those two points ...  
    !
    use module_kind_variables, only: kind_double
    use module_serial_linear_solver, only: serial_lsystem_by_svd
    use mld_logger, only: log_debug, log_warning, vtoa, mld_verbose
    implicit none 
  
    real(kind_double), intent(in) :: r1, r2, f1, f1d, f2, f2d
    real(kind_double), dimension(4,1), intent(out) :: params 

    real(kind_double), dimension(4,4) :: sysmat
    real(kind_double), dimension(4,1) :: sysyy
    integer :: nn 
    real(kind_double) :: svd_rcond_local

    _NAMECURRENT_("parameters_two_ends_function_continuity")



    nn = 4
    
    sysmat(:,1) = (/   1.d0,    1.d0,       0.d0,       0.d0 /) 
    sysmat(:,2) = (/     r1,      r2,       1.d0,       1.d0 /) 
    sysmat(:,3) = (/  r1**2,   r2**2,    2.d0*r1,    2.d0*r2 /) 
    sysmat(:,4) = (/  r1**3,   r2**3, 3.d0*r1**2, 3.d0*r2**2 /) 

    if (f1 <=  0.d0 ) then 
       call log_warning('ML: in '//NAMECURRENT//' f(r1) in r1 is negative i.e. ZBL negative at r1 : '//'f1 = '//vtoa(f1)//' r1 = '//vtoa(r1))
       call log_warning('ML: probably the continuity of ZBL function in r1 cannot be ensured. Ask the master for possible solutions') 
    end if 

    if (f2 <=  0.d0 ) then 
      call log_warning('ML: in '//NAMECURRENT//' f(r2) in r2 is negative i.e. kernel_2b negative at r2 : '//'f2 = '//vtoa(f2)//' r2 = '//vtoa(r2))
      call log_warning('ML: probably the continuity with kernel_2b function in r2 cannot be ensured. Ask the master for possible solutions.') 
    end if 


    sysyy(1:4,1) = (/ log(f1), log(f2), f1d/f1, f2d/f2 /)

    svd_rcond_local = -1.d0
    call serial_lsystem_by_svd ( sysmat, sysyy, nn, nn, params, svd_rcond_local)

    _MLD_END_

  end subroutine  parameters_two_ends_function_continuity 


  subroutine  parameters_two_ends_function_continuity_second (r1, r2, f1, f1d, f1dd, f2, f2d, f2dd, params)
    !
    !  Having the values of two functions f1 and f2 in r1 and r2 , respectively: 
    !  f1(r1) = f1, 
    !  f1'(r1)= f1d,
    !  f1''(r1) = f1dd,
    !  f2(r2) = f2, 
    !  f2'(r2)= f2d,
    !  f2''(r2)= f2dd,
    !  provides the 4 parameters params (a,b,c,d,e,f) of the function f(x) = a + b*x + c*x^2 + d*x^3 + e*x^4 + f*x^5 than 
    !  ensure the continuity of the function and the derivatives between those two points ...  
    !
    use module_kind_variables, only: kind_double
    use module_serial_linear_solver, only: serial_lsystem_by_svd
    use mld_logger, only: log_debug, log_warning, vtoa, mld_verbose
    implicit none 
    real(kind_double), intent(in) :: r1, r2, f1, f1d, f1dd, f2, f2d, f2dd 
    real(kind_double), dimension(6,1), intent(out) :: params 

    real(kind_double), dimension(6,6) :: sysmat
    real(kind_double), dimension(6,1) :: sysyy
    integer :: nn 
    real(kind_double) :: svd_rcond_local

    _NAMECURRENT_("parameters_two_ends_function_continuity")

    nn = 6 
    !this is nor correct ... 
    
    sysmat(:,1) = (/   1.d0,    1.d0,       0.d0,       0.d0,        0.d0,        0.d0 /) 
    sysmat(:,2) = (/     r1,      r2,       1.d0,       1.d0,        0.d0,        0.d0 /) 
    sysmat(:,3) = (/  r1**2,   r2**2,    2.d0*r1,    2.d0*r2,        2.d0,        2.d0 /) 
    sysmat(:,4) = (/  r1**3,   r2**3, 3.d0*r1**2, 3.d0*r2**2,     6.d0*r1,     6.d0*r2 /) 
    sysmat(:,5) = (/  r1**4,   r2**4, 4.d0*r1**3, 4.d0*r2**3, 12.d0*r1**2, 12.d0*r2**2 /) 
    sysmat(:,6) = (/  r1**5,   r2**5, 5.d0*r1**4, 5.d0*r2**4, 20.d0*r1**3, 20.d0*r2**3 /) 

    if (f1 <=  0.d0 ) then 
       call log_warning('ML: in '//NAMECURRENT//' f(r1) in r1 is negative i.e. ZBL negative at r1 : '//'f1 = '//vtoa(f1)//' r1 = '//vtoa(r1))
       call log_warning('ML: probably the continuity of ZBL function in r1 cannot be ensured. Ask the master for possible solutions') 
    end if 

    if (f2 <=  0.d0 ) then 
      call log_warning('ML: in '//NAMECURRENT//' f(r2) in r2 is negative i.e. kernel_2b negative at r2 : '//'f2 = '//vtoa(f2)//' r2 = '//vtoa(r2))
      call log_warning('ML: probably the continuity with kernel_2b function in r2 cannot be ensured. Ask the master for possible solutions.') 
    end if 


    sysyy(1:6,1) = (/ f1, f2, f1d, f2d, f1dd, f2dd  /)

    svd_rcond_local = -1.d0
    call serial_lsystem_by_svd ( sysmat, sysyy, nn, nn, params, svd_rcond_local)

    _MLD_END_

  end subroutine  parameters_two_ends_function_continuity_second




  subroutine  get_continuity_k2b_zbl(r1, r2, l_type_fcut, l_r_cut, l_rcut_width)
    use module_kind_variables, only : kind_double
    use module_chemical_species, only: size_species_half
    use derived_types, only : typ_species_half
    use module_zbl, only : params_k2b_to_zbl, type_rac_zbl, rac_zbl_exp, rac_zbl_poly
    use module_potential_zbl, only: zbl_end
    use module_form_kernel_nbody, only: kernel_2b_end
    use mld_logger
    implicit none 
    real(kind_double), intent(in) :: r1, r2, l_r_cut, l_rcut_width
    integer, intent(in) :: l_type_fcut 
    integer :: is 
    real(kind_double) :: Z_ja, Z_ia, f1, f1d, f1dd, f2, f2d, f2dd  
    real(kind_double), dimension(:,:), allocatable :: params
    integer :: type_db_ia, type_db_ja
    logical :: desc_forces_local 

    _NAMECURRENT_("get_continuity_k2b_zbl")
    desc_forces_local = .true. 
    f2dd =0.d0 
    f1dd =0.d0 
    type_rac_zbl = 2 


    select case (type_rac_zbl)
      case(rac_zbl_exp) 
        if (allocated(params_k2b_to_zbl)) deallocate(params_k2b_to_zbl)
        allocate(params_k2b_to_zbl(4, size_species_half))      
        if (allocated(params)) deallocate(params) ; allocate(params(4,1))
      case(rac_zbl_poly) 
        if (allocated(params_k2b_to_zbl)) deallocate(params_k2b_to_zbl)
        allocate(params_k2b_to_zbl(6, size_species_half))
        if (allocated(params)) deallocate(params) ; allocate(params(6,1))
      case default 
        call log_warning('No implementation for this zbl '//vtoa(type_rac_zbl)//' in subroutine '//NAMECURRENT )
        call log_critical('No implementation for this zbl '//vtoa(type_rac_zbl)//' in subroutine '//NAMECURRENT )
    end select    


    do is = 1, size_species_half
      Z_ja = typ_species_half(is)%Z1
      Z_ia = typ_species_half(is)%Z2
      type_db_ia = typ_species_half(is)%type1
      type_db_ja = typ_species_half(is)%type2 
      call zbl_end(Z_ja, Z_ia, r1, f1, f1d, f1dd) 
      call kernel_2b_end (type_db_ja, type_db_ia, desc_forces_local, l_type_fcut, r2, l_r_cut, l_rcut_width, f2, f2d, f2dd)
      select case(type_rac_zbl)
        case(rac_zbl_exp) 
           call parameters_two_ends_function_continuity(r1, r2, f1, f1d, f2, f2d, params)
           params_k2b_to_zbl(1:4,is) = params(1:4,1)
        case(rac_zbl_poly) 
           call parameters_two_ends_function_continuity_second(r1, r2, f1, f1d, f1dd, f2, f2d, f2dd, params)
           params_k2b_to_zbl(1:6,is) = params(1:6,1)
        case default 
           call log_warning('No implementation for this zbl '//vtoa(type_rac_zbl)//' in subroutine '//NAMECURRENT )
           call log_critical('No implementation for this zbl '//vtoa(type_rac_zbl)//' in subroutine '//NAMECURRENT )
      end select    
    end do 
    
    _MLD_END_

  end subroutine  get_continuity_k2b_zbl  

end   module module_continuity_k2b_zbl 
