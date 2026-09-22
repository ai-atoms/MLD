
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

!in_dev! subroutine main_get_database
!in_dev!   use ml_in_ndm_module, only: iconf_data, rangml
!in_dev!   use module_db_poscar, only: i_start_conf, i_final_conf, iconf_start_on_proc, &
!in_dev!                               iconf_end_on_proc, iconf_size_on_proc, iconf_to_proc
!in_dev! 
!in_dev!   implicit none 
!in_dev!   
!in_dev!   call prepare_database()
!in_dev!   call set_limit_for_configs(rangml, iconf_data, i_start_conf, i_final_conf, &
!in_dev!                                      iconf_start_on_proc, iconf_end_on_proc, &
!in_dev!                                      iconf_size_on_proc, iconf_to_proc)
!in_dev! 
!in_dev! end subroutine main_get_database 

module  mod_covariance_matrix
    use module_kind_variables, only: kind_double
    use mld_logger
    implicit none
    private 
    public :: covariance_matrix, obj_cov, get_matrix, get_mu 

    type :: covariance_matrix
        private
        real(kind_double), dimension(:,:), allocatable :: matrix
        real(kind_double), dimension(:), allocatable :: mean, mean_local 
        real(kind_double), dimension(:,:), allocatable :: data_local
        integer :: dim_matrix, dim_data_local, dim_data 
        integer :: ilocal 
    contains
        procedure :: init => init_covariance_matrix
        procedure :: compute_mean => compute_local_mean
        procedure :: compute_cov => compute_gather_covariance
        !procedure :: analyze => analyze_covariance_matrix
        procedure :: get_matrix
        procedure :: get_mu
    end type covariance_matrix

    type(covariance_matrix) :: obj_cov 

contains

    subroutine get_matrix(this, matrix_out)
      class(covariance_matrix), intent(inout) :: this
      real(kind_double), dimension(:,:), intent(out) :: matrix_out
      matrix_out = this%matrix
    end subroutine get_matrix

    subroutine get_mu(this, mu_out)
      class(covariance_matrix), intent(inout) :: this
      real(kind_double), dimension(:), intent(out) :: mu_out
      mu_out = this%mean
    end subroutine get_mu


    subroutine init_covariance_matrix(this,  config_real, dim_matrix, i_start_conf, i_final_conf)
      use my_mpi_subroutines, only: subworlds_allreduce_int, subworlds_allreduce_from_evrywhere_int
      use set_limits, only: set_limit_for_atoms_with_MPI_grid
      use derived_types, only: system_state 
      implicit none 
      class(covariance_matrix), intent(inout) :: this
      integer, intent(in) ::  dim_matrix, i_start_conf, i_final_conf
      type(system_state), dimension(:), intent(in) :: config_real
      integer :: ic, dlocal 
      integer :: istart, ifinal
      _NAMECURRENT_("init_covariance_matrix")
      _MLD_BEGIN_

      ! Initialize data
      this%dim_matrix = dim_matrix
      allocate(this%matrix(dim_matrix, dim_matrix))
      allocate(this%mean(dim_matrix))
      this%matrix(:,:) = 0.d0  
      this%mean(:) = 0.d0
      this%dim_data_local = 0
      this%ilocal = 0

      do ic = i_start_conf, i_final_conf
        if (.not. (config_real(ic)%train)) cycle
        if (.not.(config_real(ic)%has_energy)) cycle
        call set_limit_for_atoms_with_MPI_grid(config_real(ic)%nat, istart, ifinal)
        if ((ifinal == istart).and.(ifinal==0)) cycle
        this%dim_data_local = this%dim_data_local + ifinal - istart + 1 
      end do

      dlocal = this%dim_data_local 
      call subworlds_allreduce_from_evrywhere_int(dlocal, this%dim_data)
      allocate(this%data_local(dim_matrix, this%dim_data_local))
      this%data_local(:,:) = 0.d0
      allocate(this%mean_local(dim_matrix))
      this%mean_local(:) = 0.d0

      _MLD_END_ 
    end subroutine init_covariance_matrix

    subroutine compute_local_mean(this,  iconf, config_real, config_desc)
      use my_mpi_subroutines, only: subworlds_allreduce_from_evrywhere_vect_double
      use derived_types, only: descriptor_cell, system_state
      class(covariance_matrix), intent(inout) :: this
      integer, intent(in):: iconf
      type(descriptor_cell), dimension(:), intent(in) :: config_desc
      type(system_state), dimension(:), intent(in) :: config_real
      integer :: ja
      real(kind_double), dimension(:), allocatable :: mean 
      _NAMECURRENT_("compute_local_mean")
      _MLD_BEGIN_

      allocate(mean(this%dim_matrix))
      mean(:) = 0.d0
      if (this%dim_matrix /= size(config_desc(iconf)%energy, 1)) then
        call log_critical('Error: dimension of the descriptor does not match the dimension of the covariance matrix')
        stop
      end if     
      
      do ja=config_real(iconf)%at_start, config_real(iconf)%at_final
         if (ja == 0) cycle
         this%ilocal = this%ilocal + 1
         this%mean_local(:) = this%mean_local(:) + config_desc(iconf)%energy(:, ja) 
         this%data_local(:, this%ilocal) = config_desc(iconf)%energy(:, ja)
      end do

      _MLD_END_
    end subroutine compute_local_mean

    subroutine compute_gather_covariance(this)
      use my_mpi_subroutines, only: subworlds_allreduce_from_evrywhere_vect_double, & 
                                    subworlds_allreduce_from_evrywhere_matrix_double
      class(covariance_matrix), intent(inout) :: this
      integer :: ii, id, jd, mm 
      real(kind_double), dimension(:,:), allocatable :: local_matrix
      _NAMECURRENT_("compute_gather_covariance")
      _MLD_BEGIN_

      if (this%ilocal /= this%dim_data_local) then
        call log_critical("Error: the number of local data points does not match") 
        call log_critical("Error: the expected number in "//NAMECURRENT//vtoa(this%ilocal)//"  "//vtoa(this%dim_data_local)) 
        stop "in compute_gather_covariance"
      end if

      call subworlds_allreduce_from_evrywhere_vect_double(this%mean_local, this%mean)
      this%mean(:) = this%mean(:) / dble(this%dim_data)
      do ii = 1, this%dim_data_local
        this%data_local(:, ii) = this%data_local(:, ii) - this%mean(:)
      end do 

      allocate(local_matrix(this%dim_matrix, this%dim_matrix))
      local_matrix(:,:) = 0.d0
      ! Compute the covariance matrix
      do id = 1, this%dim_matrix
        do jd = 1, this%dim_matrix
          do mm = 1, this%dim_data_local
            local_matrix(id, jd) = local_matrix(id, jd) + this%data_local(id, mm) * this%data_local(jd, mm)
          end do
          local_matrix(id, jd) = local_matrix(id, jd) / dble(this%dim_data - 1)
        end do
      end do

      call subworlds_allreduce_from_evrywhere_matrix_double(local_matrix, this%matrix)

      deallocate(local_matrix)
      deallocate(this%data_local)

      _MLD_END_
    end subroutine compute_gather_covariance


    !$! subroutine analyze_covariance_matrix(this)
    !$!     class(covariance_matrix), intent(inout) :: this
    !$!     ! Code for analyzing the covariance matrix (e.g., finding eigenvalues)
    !$! end subroutine analyze_covariance_matrix

end module mod_covariance_matrix


subroutine main_compute_descriptors()
  use mld_logger, only: vtoa, log_info, mld_verbose, log_critical, log_debug 
  use ml_in_ndm_module, only: tmp_val_desc_max, lmask 
  use module_db_setup, only: iconf_data                             
  use temporary_data_cov, only: dim_data, dim_xdesc, dim_data_train, dim_data_test, dim_data_constraints
  use derived_types, only: config_real, config_desc
  use snap, only: i_fit_snap, dim_ene_train_lml, dim_force_train_lml, dim_stress_train_lml, &
                  i_constraints_snap
  use time_check_general, only: time_full_desc, time_calc_desc, time_neigh_desc, MY_MPI_WTIME               
  use set_limits, only: set_limit_for_configs
  use module_db_poscar, only: i_start_conf, i_final_conf, procs_per_file 
  use set_limits, only: set_limit_for_configs_with_MPI_grid
  use mld_subworld
  use module_covariance, only: train_covariance_matrix
  use mod_covariance_matrix, only: covariance_matrix, obj_cov 
  use module_ftnbody, only: init_mode_ftnbody
  use module_extxyz, only: db_xyz, read_xyz_config
  use module_json, only: db_json, read_json_config
#if(MLD_HDF5)
  use mld_hdf5, only: db_hdf5
#endif

  implicit none
  integer  :: i
  real(kind(1.d0))     :: tmp_val, t0, t1, t2, t3, t5 
  character(len=100) :: ctmp
  logical :: post_desc

  _NAMECURRENT_("main_compute_descriptors")

  _MLD_BEGIN_


  call prepare_database()
  dim_data = 0
  dim_data_constraints = 0
  dim_data_test = 0
  dim_data_train = 0
  dim_ene_train_lml = 0
  dim_force_train_lml = 0
  dim_stress_train_lml = 0

  time_neigh_desc = 0.d0
  time_full_desc = 0.d0 
  t0 = MY_MPI_WTIME()

  !computing the descriptors for the training part ...
  tmp_val = -1.d0
  !d! call set_limit_for_configs(rangml, iconf_data, i_start_conf, i_final_conf, &
  !d!                            iconf_start_on_proc, iconf_end_on_proc, &
  !d!                            iconf_size_on_proc, iconf_to_proc)
  !set-up on which proc is read each conf ...
  call init_subworld(procs_per_file)
  call set_limit_for_configs_with_MPI_grid(iconf_data, i_start_conf, i_final_conf)
  call prepare_train_dimensions()
  !call train_allocate_mld()
  i_fit_snap = 0
  i_constraints_snap = 0

  if (train_covariance_matrix) call obj_cov%init(config_real, dim_xdesc, i_start_conf, i_final_conf) 

  do i = i_start_conf, i_final_conf
    if (.not. (config_real(i)%train)) cycle
    t1 = MY_MPI_WTIME()
    !call read_poscar_sasha(trim(config_real(i)%filename),i)
    if (db_xyz) then
       call read_xyz_config(i)
    else if (db_json) then
       call read_json_config(i)
    else
#if(MLD_HDF5)
    if (db_hdf5) then
       call read_hdf5_file(i)
    else
       call read_poscar_sasha(i)
    end if
#else
    call read_poscar_sasha(i)
#endif
    end if
  
    call fix_poscar_weights
    call fix_atoms_weights(i)  
    t2 = MY_MPI_WTIME()
    time_neigh_desc = time_neigh_desc + (t2 - t1)

    if (lmask) then 
      call read_mask(i)
    !else 
    !  if (allocated(config_desc(i)%amask)) deallocate(config_desc(i)%amask)
    !  allocate(config_desc(i)%amask(config_real(i)%nat))
    end if     
    !config_real(i)%has_force=.false.
    call test_if_config_is_small(i)
    !config_real(i)%small=.true.
    call calc_neighbours(i)
    t2 = MY_MPI_WTIME()
    time_neigh_desc = time_neigh_desc + (t2 - t1)
    post_desc=.false. 
    call compute_descriptors(i, post_desc)
    if (train_covariance_matrix) then
      call  obj_cov%compute_mean(i, config_real, config_desc)
    end if
    !not_imp call val_renormalize_descritors(i)
    if (tmp_val_desc_max >= tmp_val) then
      !debug if (rangml) write (6,*) tmp_val_desc_max
      tmp_val = tmp_val_desc_max
    end if
    t3 = MY_MPI_WTIME()
    time_calc_desc = time_calc_desc + (t3 - t2)
    call train_deallocate_desc(i)
    call deallocate_real_config(i)
  end do
  t5 = MY_MPI_WTIME()
  time_full_desc = time_full_desc + (t5 - t0)

  if (train_covariance_matrix) then
    call obj_cov%compute_cov()
  end if

  !write(tmp_val, '(f20.10)' ) tmp_val_desc_max
  write(ctmp, '(es20.10)')  tmp_val 
  call log_info('ML: the max value of the descriptor is   :'// trim(ctmp))
  if  (.not.(init_mode_ftnbody)) call close_subworld()
  _MLD_END_ 

end subroutine main_compute_descriptors
