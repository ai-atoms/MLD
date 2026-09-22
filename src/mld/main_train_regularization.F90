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

subroutine main_train_regularization()

   !use mpi
   use mld_mpi, only: comm_mld

   use module_kind_variables, only: kind_double
   use ml_in_ndm_module, only: lambda_krr, rangml, debug, vector_lambda_krr, regularization_name
   use snap                ! TODOcvw, only: w_params
   use main_mld_mod, only: main_train_optimize_weights
   use mld_logger
   use time_check_general, only: MY_MPI_WTIME, time_fit_params
   use module_ml_scalapack, only: scalapack_driver
   !TORC! use my_mpi_subroutines,only: my_barrier_mld


   implicit none

   integer  :: i_lambda, ntest, irows
   character(len=70)    :: tmp_name
   real(kind_double) :: t7, t6
   logical :: post_desc

   _NAMECURRENT_("main_train_regularization")



   _MLD_BEGIN_
   call prepare_train_dimensions
   call prepare_test_dimensions
   ! the desc are computed on all procs but they are pack one one proc.
   post_desc = .true.
   call compute_descriptors_all_database(post_desc)
   ! write (*,*) '!!!!', dim_data_train

   if (rangml == 0) then
      write (6, '(">>>------------All descriptors are computed -------------------------------")')
   end if
   call prepare_train_Amatrix
   if (rangml == 0) then
      write (6, '(">>>--------------------Amatrix is filled------------------------------------")')
   end if

   do i_lambda = 1, size(vector_lambda_krr)
      ntest = 10000 + i_lambda
      write (tmp_name, '(i5)') ntest
      regularization_name = tmp_name(2:5)

      lambda_krr = vector_lambda_krr(i_lambda)
      if (rangml == 0) then
         write (6, '(">>>---------------------------------------------------------------------------")')
         write (6, '("ML:   ", a4, "     RIDGE REGRESSION              lambda_krr   ", e20.8)') regularization_name, lambda_krr
         write (6, '("------------------------------------------------------------------------------")')
      end if
      t6 = MY_MPI_WTIME()
      if (scalapack_driver) then
         call train_snap_get_parameters
      else
         if (rangml == 0) then   ! rangml=>0 running only on the proc 0 for trainning ... train 0milady ...
            if (debug) then
               if (rangml == 0) write (6, *) 'ML: train only on one proc ......'
            end if
            call train_snap_get_parameters
            if (debug) then
               if (rangml == 0) write (6, *) 'ML: train out on one proc ......'
            end if
         end if
      end if

      if (debug) then
         write (6, *) 'ML: init broadcast of all parameters ......', rangml
      end if
#if(PARA)
      if (.not.scalapack_driver) then
         !TORC! call my_barrier_mld(codeml)
         call comm_mld%barrier
         !TORC! call MPI_BCAST(w_params(:, 1), size(w_params, 1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, codeml)
         call comm_mld%bcast(0, w_params(:, 1))
         !TODO Amat - here are obsolete lines ...
         do irows = 1, size(Amat, 2)
            !TORC! call MPI_BCAST(Amat(:, irows), size(Amat, 1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, codeml)
            call comm_mld%bcast(0, Amat(:, irows))
         end do
      end if
#endif
      !TORC! call my_barrier_mld(codeml)
      call comm_mld%barrier
      t7 = MY_MPI_WTIME()
      time_fit_params = time_fit_params + (t7 - t6)


      if (scalapack_driver) then
         if (rangml==0) regularization_name = ""
         call get_train_errors
      else
         if (rangml == 0) then
            regularization_name = ""
            ! PAY ATTENTION THAT THIS IS also THE OPTIMIZATION ON WEIGTHS ... i DO NOT KONW why ....
            call get_train_errors
         end if                  ! rangml=>0 running on the proc 0 for trainning ... train 0milady ...
      end if
      if (rangml == 0) write (6, '("------------------------------------------------------------------------------")')

      !call get_train_errors
      !if (rangml == 0) write (6, '("------------------------------------------------------------------------------")')
      call get_test_errors
      if (rangml == 0) write (6, '("<<<--------------------------------------------------------------------------|")')
      if (rangml == 0) write (6, '("   ")')
   end do


   _MLD_END_

end subroutine main_train_regularization



subroutine prepare_train_Amatrix()

   use ml_in_ndm_module, only:  mld_fit_type, &
      fit_lapack_qr_constraints,  snap_class_constraints, &
      ml_type_analysis, ml_type, i_start_at, i_final_at
   use module_db_setup, only: iconf_data
   use derived_types, only: config_real
   use snap, only: i_fit_snap
   use main_mld_mod
   use mld_energy_mod, only: train_fill_Amat_with_energy
   use module_mld_quadratic, only: i_e_fit_snap, i_f_fit_snap, i_s_fit_snap
   use mld_logger

   implicit none

   integer  :: i

   _NAMECURRENT_("prepare_train_Amatrix")



   _MLD_BEGIN_
   !testing the training
   !debug write (*,*) dim_data, dim_data_test,  dim_data_train
   call train_allocate_mld()
   i_fit_snap = 0
   i_e_fit_snap = 0
   i_f_fit_snap = 0
   i_s_fit_snap = 0
   !computing the descriptors for the training part ...
   do i = 1, iconf_data
      if (.not. (config_real(i)%train)) cycle
      if (.not. (ml_type == ml_type_analysis)) then
         i_start_at = config_real(i)%at_start
         i_final_at = config_real(i)%at_final
         call train_fill_Amat_with_energy(i, .true.)
         call train_fill_Amat_with_force(i, .true.)
         call train_fill_Amat_with_stress(i, .true.)
         ! if constraints ...
         if ((mld_fit_type == fit_lapack_qr_constraints) .and. config_real(i)%class == snap_class_constraints) then
            call train_fill_Bmat_constraints(i)
         end if
      end if
   end do

   _MLD_END_

end subroutine prepare_train_Amatrix



subroutine prepare_train_dimensions()
   use module_kind_variables, only: kind_double
   use ml_in_ndm_module, only: rangml, type_of_loss, loss_init, loss_fair, loss_per, ml_type, ml_type_descriptors
   use module_db_setup, only: iconf_data
   use temporary_data_cov, only: dim_data, dim_data_train, dim_data_constraints, dim_data_train_local
   use derived_types, only: config_real, db_model
   use snap, only: dim_ene_train_lml, dim_force_train_lml, dim_stress_train_lml
   use module_optimization, only: itopt, optimize_weights_db,  factor_energy_error, factor_force_error, factor_stress_error
   use main_mld_mod
   use module_kind_variables, ONLY: kind_double
   use time_check_general, only: MY_MPI_WTIME, time_read_db
   use module_db_poscar, only: i_start_conf, i_final_conf
   use my_mpi_subroutines, only : subworlds_allreduce_int, subworlds_allreduce_vect_int
   use mld_logger
   use mld_subworld, only: id_subworld, subrank
   use module_extxyz, only: db_xyz, read_xyz_config
   use module_json, only: db_json, read_json_config
#if(MLD_HDF5)
   use mld_hdf5, only: db_hdf5
#endif
   !TORC! use my_mpi_subroutines, only: my_barrier_mld

   implicit none
   integer  :: i
   real(kind_double) :: no_ene, no_for, no_str, no_sum, tt0, tt1
   integer, dimension(:), allocatable :: tmp_nat


   _NAMECURRENT_("prepare_train_dimensions")


   _MLD_BEGIN_

   dim_data = 0
   dim_data_constraints = 0

   dim_data_train = 0
   dim_ene_train_lml = 0
   dim_force_train_lml = 0
   dim_stress_train_lml= 0

   config_real(:)%nat = 0

   !write(chlog,"(5i8)") rangml, subrank, id_subworld, i_start_conf, i_final_conf
   !call log_info("ML grid info w-rank, sub-rank, master_rank, i_start_conf, ifinal_conf "//trim(chlog))
   write(*,'("ML grid info w-rank, sub-rank, id-sub, i_start_conf, ifinal_conf :", 5i6)') rangml, subrank, id_subworld, i_start_conf, i_final_conf

   !set-up the dimensions of train matrix Amat
   tt0 = MY_MPI_WTIME()

   do i = i_start_conf, i_final_conf
      if (.not. (config_real(i)%train)) cycle
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
      ! fix weights if they are ...
      call fix_poscar_weights
      call fix_atoms_weights(i)
      call get_fit_dimensions(i)

      !major clean-up for large boxes ... 
      if (ml_type == ml_type_descriptors) then
        if (allocated(config_real(i)%pos_cart)) deallocate (config_real(i)%pos_cart)
        if (allocated(config_real(i)%pos_crst)) deallocate (config_real(i)%pos_crst)
        if (allocated(config_real(i)%force)) deallocate (config_real(i)%force)
        if (allocated(config_real(i)%atomic_spin)) deallocate (config_real(i)%atomic_spin)
        if (allocated(config_real(i)%itype)) deallocate (config_real(i)%itype)
        if (allocated(config_real(i)%itype_db)) deallocate (config_real(i)%itype_db)
      end if
   end do

   ! dim_data_train_local: dim_data_train on each MPI subgroup
   dim_data_train_local = dim_data_train
   !TORC! call my_barrier_mld(codeml)
   call comm_mld%barrier
   call subworlds_allreduce_int(dim_ene_train_lml)
   call subworlds_allreduce_int(dim_force_train_lml)
   call subworlds_allreduce_int(dim_stress_train_lml)
   call subworlds_allreduce_int(dim_data)
   call subworlds_allreduce_int(dim_data_train)

   ! everyone will need this info during train error
   allocate(tmp_nat(iconf_data))
   tmp_nat(:) = config_real(:)%nat
   call subworlds_allreduce_vect_int(tmp_nat)
   config_real(:)%nat = tmp_nat(:)
   deallocate(tmp_nat)

   tt1 = MY_MPI_WTIME()
   time_read_db = time_read_db + tt1 - tt0

   call set_mld_dimension()

   if (itopt == 0 ) then
      if (rangml == 0) then
         write (6, '("ML: the number of datapoints in train, n_train, n_E, n_F, n_S:   ", 4i8)') dim_data_train, dim_ene_train_lml, dim_force_train_lml, dim_stress_train_lml
         !write (6,'("ML: the number of datapoints in test,  n_test                :   ",  i8)') dim_data-dim_data_train
      end if
   end if

   if (dim_ene_train_lml > 0 ) then
      no_ene = dim_ene_train_lml
   else
      no_ene = 1
   end if

   if (dim_force_train_lml > 0 ) then
      no_for = dim_force_train_lml
   else
      no_for = 1
   end if


   if (dim_stress_train_lml > 0 ) then
      no_str = dim_stress_train_lml
   else
      no_str = 1
   end if

   no_sum = no_ene + no_for + no_str
   no_ene = dble(no_ene) / dble(no_sum)
   no_for = dble(no_for) / dble(no_sum)
   no_str = dble(no_str) / dble(no_sum)

   do i = i_start_conf, i_final_conf
      if (.not. (config_real(i)%train)) cycle

      if (dim_ene_train_lml > 0 ) then
         if (type_of_loss == loss_init) then
            config_real(i)%w_e_model = config_real(i)%w_e*factor_energy_error
            if (optimize_weights_db) config_real(i)%w_e_end_model = config_real(i)%w_e_end*factor_energy_error
         else if (type_of_loss == loss_fair ) then
            config_real(i)%w_e_model = config_real(i)%w_e / no_ene*factor_energy_error
            if (optimize_weights_db) config_real(i)%w_e_end_model = config_real(i)%w_e_end / no_ene*factor_energy_error
         else if (type_of_loss == loss_per ) then
            config_real(i)%w_e_model = config_real(i)%w_e / dble(config_real(i)%nat)*factor_energy_error
            if (optimize_weights_db) config_real(i)%w_e_end_model = config_real(i)%w_e_end / dble(config_real(i)%nat)*factor_energy_error
         else
            stop 'unknown loss type type_of_loss e'
         end if
      else
         config_real(i)%w_e_model = config_real(i)%w_e
         if (optimize_weights_db) config_real(i)%w_e_end_model = config_real(i)%w_e_end
      end if

      if (dim_force_train_lml > 0 ) then
         if (type_of_loss == loss_init) then
            config_real(i)%w_f_model = config_real(i)%w_f*factor_force_error
            if (optimize_weights_db) config_real(i)%w_f_end_model = config_real(i)%w_f*factor_force_error
         else if (type_of_loss == loss_fair ) then
            config_real(i)%w_f_model = config_real(i)%w_f / no_for*factor_force_error
            if (optimize_weights_db) config_real(i)%w_f_end_model = config_real(i)%w_f_end / no_for*factor_force_error
         else if (type_of_loss == loss_per ) then
            config_real(i)%w_f_model = config_real(i)%w_f / dble(3 * config_real(i)%nat)*factor_force_error
            if (optimize_weights_db) config_real(i)%w_f_end_model = config_real(i)%w_f_end  / dble(3 * config_real(i)%nat)*factor_force_error
         else
            stop 'unknown loss type type_of_loss f'
         end if
      else
         config_real(i)%w_f_model = config_real(i)%w_f
         if (optimize_weights_db)  config_real(i)%w_f_end_model = config_real(i)%w_f_end
      end if

      if (dim_stress_train_lml > 0 ) then
         if (type_of_loss == loss_init) then
            config_real(i)%w_s_model = config_real(i)%w_s*factor_stress_error
            if (optimize_weights_db) config_real(i)%w_s_end_model = config_real(i)%w_s_end*factor_stress_error
         else if (type_of_loss == loss_fair ) then
            config_real(i)%w_s_model = config_real(i)%w_s / no_str*factor_stress_error
            if (optimize_weights_db) config_real(i)%w_s_end_model = config_real(i)%w_s_end / no_str*factor_stress_error
         else if (type_of_loss == loss_per ) then
            config_real(i)%w_s_model = config_real(i)%w_s / dble(3 * config_real(i)%nat)*factor_stress_error
            if (optimize_weights_db) config_real(i)%w_s_end_model = config_real(i)%w_s_end / dble(3 * config_real(i)%nat)*factor_stress_error
         else
            stop 'unknown loss type type_of_loss s'
         end if
      else
         config_real(i)%w_s_model = config_real(i)%w_s
         if (optimize_weights_db) config_real(i)%w_s_end_model = config_real(i)%w_s_end
      end if

      db_model(config_real(i)%db_line)%w_e_model = config_real(i)%w_e_model
      if (optimize_weights_db) db_model(config_real(i)%db_line)%w_e_end_model = config_real(i)%w_e_end_model

      db_model(config_real(i)%db_line)%w_f_model = config_real(i)%w_f_model
      if (optimize_weights_db) db_model(config_real(i)%db_line)%w_f_end_model = config_real(i)%w_f_end_model

      db_model(config_real(i)%db_line)%w_s_model = config_real(i)%w_s_model
      if (optimize_weights_db) db_model(config_real(i)%db_line)%w_s_end_model = config_real(i)%w_s_end_model
   end do

   _MLD_END_

end subroutine prepare_train_dimensions



subroutine prepare_test_dimensions()
   ! prepare dimension of the fit from configuration by directly reading the poscar configurations
   use module_kind_variables, only: kind_double
   use ml_in_ndm_module, only: rangml, debug
   use temporary_data_cov, only: dim_data_test
   use derived_types, only: config_real
   use snap, only: dim_ene_test_snap, dim_force_test_snap, dim_stress_test_snap
   use main_mld_mod
   use time_check_general, only: MY_MPI_WTIME, time_read_db
   use module_db_poscar, only: i_start_conf, i_final_conf
   use module_db_setup, only: iconf_data
   use my_mpi_subroutines, only : subworlds_allreduce_int, subworlds_allreduce_vect_int
   use mld_logger
   use module_extxyz, only: db_xyz, read_xyz_config
   use module_json, only: db_json, read_json_config
#if(MLD_HDF5)
   use mld_hdf5, only: db_hdf5
#endif

   implicit none

   integer  :: i
   real(kind_double) :: tt0, tt1
   integer, allocatable :: tmp_nat(:)

   _NAMECURRENT_("prepare_test_dimensions")


   _MLD_BEGIN_
   dim_data_test = 0
   dim_ene_test_snap = 0
   dim_force_test_snap = 0
   dim_stress_test_snap = 0

   tt0=MY_MPI_WTIME()
   do i = i_start_conf, i_final_conf
      if (config_real(i)%train) cycle
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
      call get_fit_dimensions(i)
      ! PAY ATTENTION OF THE FACT THAT AT THE CURRENT POINT dim_data_train WILL BE WRONG
   end do

   call subworlds_allreduce_int(dim_data_test)
   call subworlds_allreduce_int(dim_ene_test_snap)
   call subworlds_allreduce_int(dim_force_test_snap)
   call subworlds_allreduce_int(dim_stress_test_snap)

   ! Make the number of atoms (nat) of every TEST config available on
   ! every rank. write_test_output / test_error_by_class run on rank 0
   ! over the whole database, and the force output relies on nat to
   ! know the per-config block size. Each test config is read by a
   ! single subworld, so summing over subworlds reconstructs nat
   ! globally. Train configs already carry a global nat (set in
   ! prepare_train_dimensions) and are left untouched here.
   allocate(tmp_nat(iconf_data))
   tmp_nat(:) = 0
   do i = i_start_conf, i_final_conf
      if (.not. config_real(i)%train) tmp_nat(i) = config_real(i)%nat
   end do
   call subworlds_allreduce_vect_int(tmp_nat)
   do i = 1, iconf_data
      if (.not. config_real(i)%train) config_real(i)%nat = tmp_nat(i)
   end do
   deallocate(tmp_nat)

   tt1=MY_MPI_WTIME()
   time_read_db = time_read_db + tt1 - tt0

   if (rangml == 0) then
      write (6, '("ML: the number of datapoints in  test,  n_test, n_E, n_F, n_S:   ", 4i8)') dim_data_test, dim_ene_test_snap, dim_force_test_snap, dim_stress_test_snap
      !write (6,'("ML: the number of datapoints in test,  n_test                :   ",  i8)') dim_data-dim_data_train
   end if
   if (debug) then
      if (rangml == 0) then
         write (6, '("ML: Dim for test dim_ene_test_snap, dim_force_test_snap, dim_stress_test_snap", 3i10 )') dim_ene_test_snap, dim_force_test_snap, dim_stress_test_snap
      end if
   end if

   _MLD_END_

end subroutine prepare_test_dimensions



subroutine get_train_errors()
   use ml_in_ndm_module, only: rangml, rangml,  &
                              regularization_name
   use temporary_data_cov, only: dim_data_train
   use snap, only: dim_ene_train_lml, dim_force_train_lml, dim_stress_train_lml, &
      i_e_train_snap, i_f_train_snap, i_s_train_snap, &
      y_e_train_snap, y_f_train_snap, y_s_train_snap, &
      y_e_train_base, y_f_train_base, y_s_train_base, &
      y_e_p_a_train_base, y_e_p_a_train_snap
   use module_optimization, only: itopt, optimize_weights_db, optimize_weights_chem
   use main_mld_mod
   use mld_logger

   implicit none

   integer  :: i, einp, finp, sinp

   _NAMECURRENT_("get_train_errors")



   _MLD_BEGIN_
   !testing the training

   if (dim_ene_train_lml > 0) then
      if (allocated(y_e_train_snap)) deallocate (y_e_train_snap); allocate (y_e_train_snap(dim_ene_train_lml))
      if (allocated(y_e_train_base)) deallocate (y_e_train_base); allocate (y_e_train_base(dim_ene_train_lml))
      if (allocated(y_e_p_a_train_snap)) deallocate (y_e_p_a_train_snap); allocate (y_e_p_a_train_snap(dim_ene_train_lml))
      if (allocated(y_e_p_a_train_base)) deallocate (y_e_p_a_train_base); allocate (y_e_p_a_train_base(dim_ene_train_lml))
      if (rangml==0) open (file="train_energy"//trim(regularization_name)//".out", newunit=einp, action='write')
   end if


   if (dim_force_train_lml > 0) then
      if (allocated(y_f_train_snap)) deallocate (y_f_train_snap); allocate (y_f_train_snap(dim_force_train_lml))
      if (allocated(y_f_train_base)) deallocate (y_f_train_base); allocate (y_f_train_base(dim_force_train_lml))
      if (rangml==0) open (file="train_force"//trim(regularization_name)//".out", newunit=finp, action='write')
   end if


   if (dim_stress_train_lml > 0) then
      if (allocated(y_s_train_snap)) deallocate (y_s_train_snap); allocate (y_s_train_snap(dim_stress_train_lml))
      if (allocated(y_s_train_base)) deallocate (y_s_train_base); allocate (y_s_train_base(dim_stress_train_lml))
      if (rangml==0)  open (file="train_stress"//trim(regularization_name)//".out", newunit=sinp, action='write')
   end if

   i_e_train_snap = 0
   i_f_train_snap = 0
   i_s_train_snap = 0
   do i = 1, dim_data_train
      call train_snap_compute_energy_force_stress(i, einp, finp, sinp)
   end do

   close(einp)
   close(finp)
   close(sinp)

   if ((.not. (optimize_weights_db.or.optimize_weights_chem)).or.(itopt==0)) then
      if (rangml == 0) then
         write (6, '("ML:        ERROR TYPE                             R^2             D             RMSE          MAE    ")')
      end if
   end if
   call train_error_mld()

   !$! call train_error_by_class

   !$! if (debug) then 
   !$!   write(6,*) "End of get train errors", rangml 
   !$! end if 


   _MLD_END_

end subroutine get_train_errors



subroutine compute_descriptors_all_database(post_desc)

   use ml_in_ndm_module, only: rangml,  tmp_val_desc_max
   use main_mld_mod
   use mld_energy_mod, only: pack_energy_descriptor
   use mld_stress_mod, only: pack_stress_descriptor
   use mld_force_mod, only: pack_force_descriptor
   use module_db_poscar, only: i_start_conf, i_final_conf
   use mld_logger

   implicit none
   logical, intent(in) :: post_desc

   integer :: i
   real(kind(0.d0)) :: tmp_val

   _NAMECURRENT_("compute_descriptors_all_database")



   _MLD_BEGIN_
   tmp_val = -1.d0
   do i = i_start_conf, i_final_conf
      call test_if_config_is_small(i)
      call calc_neighbours(i)
      !
      call compute_descriptors(i, post_desc)
      call pack_energy_descriptor(i)
      call pack_force_descriptor(i)
      call pack_stress_descriptor(i)
      !not_imp call val_renormalize_descritors(i)
      !
      if (tmp_val_desc_max >= tmp_val) then
         tmp_val = tmp_val_desc_max
      end if

      call train_deallocate_desc(i)
   end do

   if (rangml == 0) then
      write (6, *) 'ML: the max value of the descriptor is   :', tmp_val
   end if

   _MLD_END_

end subroutine compute_descriptors_all_database




subroutine get_test_errors()

   use ml_in_ndm_module, only: regularization_name
   use module_db_setup, only: iconf_data
   use derived_types, only: config_real
   use snap, only: y_e_test_snap, y_f_test_snap, y_s_test_snap, &
      y_e_test_base, y_f_test_base, y_s_test_base, &
      y_e_p_a_test_snap, y_e_p_a_test_base, &
      ene_snap, fp_snap, stress_snap, dim_ene_test_snap, &
      dim_force_test_snap, dim_stress_test_snap, &
      fp_snap
   use main_mld_mod
   use module_db_poscar, only: iread_energy
   use module_md_mld, only: md_allocate_mld_desc
   use mld_logger

   implicit none

   integer  :: i, ix, einp, finp, sinp
   integer  :: i_e_test_snap, i_f_test_snap, i_s_test_snap

   _NAMECURRENT_("get_test_errors")



   _MLD_BEGIN_
   !testing the training ...
   einp = 51
   finp = 52
   sinp = 53
   if (dim_ene_test_snap > 0) then
      open (file="test_energy"//regularization_name//".out", unit=einp, action='write')
      if (allocated(y_e_test_snap)) deallocate (y_e_test_snap); allocate (y_e_test_snap(dim_ene_test_snap))
      if (allocated(y_e_test_base)) deallocate (y_e_test_base); allocate (y_e_test_base(dim_ene_test_snap))
      if (allocated(y_e_p_a_test_snap)) deallocate (y_e_p_a_test_snap); allocate (y_e_p_a_test_snap(dim_ene_test_snap))
      if (allocated(y_e_p_a_test_base)) deallocate (y_e_p_a_test_base); allocate (y_e_p_a_test_base(dim_ene_test_snap))
   end if

   if (dim_force_test_snap > 0) then
      open (file="test_force"//regularization_name//".out", unit=finp, action='write')
      if (allocated(y_f_test_snap)) deallocate (y_f_test_snap); allocate (y_f_test_snap(dim_force_test_snap))
      if (allocated(y_f_test_base)) deallocate (y_f_test_base); allocate (y_f_test_base(dim_force_test_snap))
   end if

   if (dim_stress_test_snap > 0) then
      open (file="test_stress"//regularization_name//".out", unit=sinp, action='write')
      if (allocated(y_s_test_snap)) deallocate (y_s_test_snap); allocate (y_s_test_snap(dim_stress_test_snap))
      if (allocated(y_s_test_base)) deallocate (y_s_test_base); allocate (y_s_test_base(dim_stress_test_snap))
   end if

   i_e_test_snap = 0
   i_f_test_snap = 0
   i_s_test_snap = 0
   do i = 1, iconf_data
      if (config_real(i)%train) cycle
      !call read_poscar_sasha(i)
      call test_if_config_is_small(i)
      call calc_neighbours(i)
      call md_allocate_mld_desc(fp_snap)
      if (config_real(i)%has_energy) then
         i_e_test_snap = i_e_test_snap + 1
         call md_mld_compute_energy(i, .false.)
         y_e_test_snap(i_e_test_snap) = ene_snap
         y_e_test_base(i_e_test_snap) = config_real(i)%energy(iread_energy)
         y_e_p_a_test_snap(i_e_test_snap) = ene_snap/dble(config_real(i)%nat)
         y_e_p_a_test_base(i_e_test_snap) = config_real(i)%energy(iread_energy)/dble(config_real(i)%nat)
         write (einp, '(2e20.10," ", (a), "  ", (a))') ene_snap, y_e_test_base(i_e_test_snap), config_real(i)%class, config_real(i)%filename
      end if
      if (config_real(i)%has_force) then
         call md_mld_compute_force(i, .false.)
         do ix = 1, 3
            y_f_test_snap(i_f_test_snap + (ix - 1)*config_real(i)%nat + 1:i_f_test_snap + ix*config_real(i)%nat) = fp_snap(ix, 1:config_real(i)%nat)
            y_f_test_base(i_f_test_snap + (ix - 1)*config_real(i)%nat + 1:i_f_test_snap + ix*config_real(i)%nat) = config_real(i)%force(ix, 1:config_real(i)%nat)
         end do
         do ix = 1, 3*config_real(i)%nat
            write (finp, '(2e20.10," ", (a), "  ", (a))') y_f_test_snap(i_f_test_snap + ix), y_f_test_base(i_f_test_snap + ix), config_real(i)%class, config_real(i)%filename
         end do
         i_f_test_snap = i_f_test_snap + 3*config_real(i)%nat
      end if

      if (config_real(i)%has_stress) then
         call md_mld_compute_stress(i, .false.)
         do ix = 1, 3
            y_s_test_snap(i_s_test_snap + 1:i_s_test_snap + 6) = stress_snap(1:6)
            y_s_test_base(i_s_test_snap + 1:i_s_test_snap + 6) = config_real(i)%stress(1:6)
         end do
         do ix = 1, 6
            write (sinp, '(2e20.10," ", (a),"  ", (a))') y_s_test_snap(i_s_test_snap + ix), y_s_test_base(i_s_test_snap + ix), config_real(i)%class, config_real(i)%filename
         end do
         i_s_test_snap = i_s_test_snap + 6
      end if

      !call train_deallocate_desc(i)
   end do

   close (einp)
   close (finp)
   close (sinp)
   !compute the test
   call test_error_snap

   _MLD_END_

end subroutine get_test_errors
