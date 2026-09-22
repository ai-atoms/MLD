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

module module_evaluate_parameters
   use module_kind_variables, only: kind_double
   implicit none

contains

   subroutine product_w_params_Amat(outval, idata)
      use module_ml_scalapack, only: scalapack_driver, sca_Amat, desc_sca_Amat,  dimr_sca_Amat, &
         sca_w_params, desc_sca_w_params


      use snap, only : Amat, w_params

      integer, intent(in) :: idata
      real(kind_double), intent(inout) :: outval
      ! ToScaRemove
      !real(kind_double) :: one=1.d0
      real(kind_double) :: outval2
      !end ToScaRemove


      ! scaTODO ToScaRemove - put the right output.
      if (scalapack_driver) then
         ! mx1 = mxn * nx1 y = a * x
         ! 1x1 = 1xn * nx1
         !call pddot(n, dot, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
         call pddot(dimr_sca_Amat, outval2, sca_Amat, 1, idata, desc_sca_Amat, 1, &
            sca_w_params, 1, 1, desc_sca_w_params, 1)
         !if (iam==0) write(6,*) 'pddot after '

         !$! ! ToScaRemove
         !$! outval = dot_product (w_params(:,1), Amat(:,idata))
         !$! ! end ToScaRemove
      else
         outval = dot_product (w_params(:,1), Amat(:,idata))

         !$! ! ToScaRemove
         !$! outval2 = outval
         !$! ! end ToScaRemove
      end if


      !$! ! ToScaRemove  NOT FORGET TO SEND OUT outval  and remove outval2
      !$! if (dabs(outval-outval2) > 1.d-13 ) then
      !$!   if (iam == 0) write(6,*) 'ML WARNING WARNING : serious problems in module_evaluate_parameters '
      !$! f (iam==0) write(6,'(2ES30.15E3)') outval-outval2, outval
      !$!  if
      if (scalapack_driver) then
         outval = outval2
      end if
      !$! !end ToScaRemove

   end subroutine product_w_params_Amat

end module module_evaluate_parameters



module main_mld_mod

   !/------------------------------------------------------------\
   !                                                             !
   !                           TRAIN PART                        !
   !                                                             !
   !\------------------------------------------------------------/

#if(PARA)
   !use mpi
   use mld_mpi
#endif
   !TODO to remove this s**t
   !use snap_interface, only: train_fill_Amat_with_force
   use snap, only: Amat, Bmat, zmat, ymat, &
      i_fit_snap, i_constraints_snap, w_params, &
      dim_ene_constraints, dim_force_constraints, dim_stress_constraints, &
      dim_ene_train_lml, dim_force_train_lml, dim_stress_train_lml, &
      dim_ene_test_snap, dim_force_test_snap, dim_stress_test_snap, &
      dim_xdesc_linear, &
      dim_design_line, dim_design_column, &
      fit_snap, weights_snap, y_zbl_train, &
      ene_snap, stress_snap, fp_snap, &
      y_e_test_base, y_e_test_snap, y_e_p_a_test_base, y_e_p_a_test_snap, &
      y_e_train_base, y_e_train_snap, y_e_p_a_train_base, y_e_p_a_train_snap, &
      y_s_test_base, y_s_test_snap, y_s_train_base, y_s_train_snap, &
      y_f_test_base, y_f_test_snap, y_f_train_base, y_f_train_snap, &
      train_mae_energy, train_rmse_energy, train_mae_force, train_rmse_force, &
      train_mae_stress, train_rmse_force, train_rmse_stress, &
      i_e_train_snap, i_f_train_snap, i_s_train_snap

   use mld_energy_mod
   use mld_stress_mod     ! TODOcvw, only: i_fit_snap,  i_constraints_snap,  dim_force_constraints, dim_stress_constraints, Amat, w_params
   use mld_force_mod, only: train_fill_Amat_with_force
   use mld_logger

   implicit none

contains

   subroutine main_train_gather_mld()
      ! The purpose of this subroutine is to gather some variables from the different MPI groups. 
      ! The object garthered are: fit_snap, ymat, yfunc_train, weights_snap,  Amat

      !use mpi
      
      use snap, only: i_fit_snap, fit_snap_type, Amat, ymat, weights_snap, dim_design_line, y_zbl_train
      use temporary_data_cov, only: yfunc_train, dim_xdesc
      use module_ml_scalapack, only: scalapack_driver, sca_Amat, desc_sca_Amat, &
                                     nprocs_ml_sca
      use module_init_ScaMatrix, only: init_distributed_ScaMatrix_with_LocalMatrix_from_group, &
                                       init_distributed_ScaMatrix_with_LocalMatrix_from_proc, & 
                                       init_distributed_ScaMatrix_with_LocalMatrix_mine, & 
                                       init_distributed_ScaMatrix_with_LocalMatrix_Anida
      use my_mpi_subroutines, only: subworlds_allreduce_matrix_double, subworlds_allreduce_vect_double, &
                                    subworlds_allreduce_vect_int, subworlds_allreduce_vect_logical
      use mld_subworld
      use time_check_general, only: MY_MPI_WTIME, time_fill_desc_02A, time_fill_desc_02obj, time_fill_desc_02y
      use module_db_poscar, only: procs_per_file
      

      implicit none
      integer  :: start, end, group
      integer  :: i
      integer, dimension(:), allocatable :: i_fit_snap_on_MPI_group, i_fit_snap_offsets
      integer, dimension(:), allocatable :: dimI_local, dimJ_local, Iini, Jini
      real(kind=kind(0.d0)), dimension(:,:), allocatable :: tmp_mat
      real(kind=kind(0.d0)), dimension(:), allocatable :: tmp_arr
      integer, dimension(:), allocatable :: tmp_int
      logical, dimension(:), allocatable :: tmp_log
      real(kind=kind(0.d0)), dimension(:), allocatable :: tmp_real
      real(kind=kind(0.d0)) :: timea, timeb 

      !TORC! call my_barrier_mld(codeml)
      call comm_mld%barrier()

      ! assembling all local objects into global object:
      ! i_fit_snap_on_MPI_group contains the number of columns filled in each local object
      ! i_fit_snap_offsets contains the correct i_fit_snap (i.e column offset) for each MPI group in the global object
      timea = MY_MPI_WTIME()

      allocate (i_fit_snap_on_MPI_group(0:nb_subworlds - 1))
      allocate (i_fit_snap_offsets(0:nb_subworlds - 1))
      i_fit_snap_on_MPI_group(:) = 0
      i_fit_snap_on_MPI_group(id_subworld) = i_fit_snap
      call subworlds_allreduce_vect_int(i_fit_snap_on_MPI_group)
      i_fit_snap_offsets(0) = 0
      do i=1, nb_subworlds - 1
         i_fit_snap_offsets(i) = SUM(i_fit_snap_on_MPI_group(0:i-1))
      end do
      start = i_fit_snap_offsets(id_subworld)+1
      end = start + i_fit_snap - 1


      timeb = MY_MPI_WTIME()
      time_fill_desc_02y = time_fill_desc_02y + (timeb -timea)

      if (scalapack_driver) then
         allocate(dimI_local(nprocs_ml_sca)); allocate(dimJ_local(nprocs_ml_sca))
         allocate(Iini(nprocs_ml_sca)); allocate(Jini(nprocs_ml_sca))
         ! i_fit_snap_on_MPI_group and i_fit_snap_offsets contain data regarding MPI groups, for scalapack we need info on each MPI process
         dimI_local = dim_design_line; Iini = 0
         do i=1,nprocs_ml_sca
            group = (i-1) / procs_per_file
            dimJ_local(i) = i_fit_snap_on_MPI_group(group)
            Jini(i) = i_fit_snap_offsets(group)
         end do
         !$! if (mld_rank==0) then
         !$!   write(*,'("mld_rank ", i6, " i_fit_snap_offset       :",3i7)')  mld_rank, i_fit_snap_offsets(0:nb_subworlds - 1)
         !$!   write(*,'("mld_rank ", i6, " i_fit_snap_on_MPI_group :", 3i7)') mld_rank, i_fit_snap_on_MPI_group(0:nb_subworlds-1)
         !$! end if 
         !$! 
         !$! write(*, '("mld_rank ", i6," group sca ", 2i6," Amat_local ", 2i7 )') mld_rank, iproc_sca, int((iproc_sca)/procs_per_file), size(Amat,2), dimJ_local(iproc_sca+1)
       
         
         call init_distributed_ScaMatrix_with_LocalMatrix_from_group(sca_Amat, desc_sca_Amat, Iini, Jini, Amat, dimI_local, dimJ_local)

         deallocate(dimI_local); deallocate(dimJ_local)
         deallocate(Iini); deallocate(Jini)
      else
         ! translating the local portion into the correct slot in the global matrix

         allocate(tmp_mat(dim_design_line,i_fit_snap))
         tmp_mat(:,:) = Amat(:,1:i_fit_snap); 
         Amat(:,:) = 0.d0; 
         Amat(:,start:end) = tmp_mat(:,:)
         deallocate(tmp_mat)
         call subworlds_allreduce_matrix_double('Amat', Amat)
      end if

      timea = MY_MPI_WTIME()
      time_fill_desc_02A = time_fill_desc_02A + (timea - timeb) 

      ! translating the local portion into the correct slot in the global vector
      allocate(tmp_arr(i_fit_snap))

      tmp_arr(:) = ymat(1:i_fit_snap,1); ymat(:,1) = 0.d0; ymat(start:end,1) = tmp_arr(:)
      tmp_arr(:) = yfunc_train(1:i_fit_snap); yfunc_train(:) = 0.d0; yfunc_train(start:end) = tmp_arr(:)
      tmp_arr(:) = weights_snap(1:i_fit_snap); weights_snap(:) = 0.d0; weights_snap(start:end) = tmp_arr(:)
      tmp_arr(:) = y_zbl_train(1:i_fit_snap); y_zbl_train(:) = 0.d0; y_zbl_train(start:end) = tmp_arr(:)

      deallocate(tmp_arr)
      !debug write(*,*) start, end ,rangml, id_subworld, ymat(32385 ,1 )

      call subworlds_allreduce_matrix_double('ymat', ymat)
      call subworlds_allreduce_vect_double(yfunc_train)
      call subworlds_allreduce_vect_double(weights_snap)
      call subworlds_allreduce_vect_double(y_zbl_train)

      allocate(tmp_int(dim_design_column))
      tmp_int = 0; tmp_int(start:end) = fit_snap(1:i_fit_snap)%iconf
      call subworlds_allreduce_vect_int(tmp_int)
      fit_snap(:)%iconf = tmp_int

      !debugMCM! tmp_int = 0; tmp_int(start:end) = fit_snap(1:i_fit_snap)%force
      !debugMCM! call subworlds_allreduce_vect_int(tmp_int)      
      !debugMCM! fit_snap(:)%force = tmp_int
      allocate(tmp_log(dim_design_column))
      tmp_log = .false.; tmp_log(start:end) = fit_snap(1:i_fit_snap)%force
      call subworlds_allreduce_vect_logical(tmp_log)
      fit_snap(:)%force = tmp_log

      !debugMCM! tmp_int = 0; tmp_int(start:end) = fit_snap(1:i_fit_snap)%energy
      !debugMCM! call subworlds_allreduce_vect_int(tmp_int)
      !debugMCM! fit_snap(:)%energy = tmp_int
      tmp_log = .false.; tmp_log(start:end) = fit_snap(1:i_fit_snap)%energy
      call subworlds_allreduce_vect_logical(tmp_log)
      fit_snap(:)%energy = tmp_log
   
      !debugMCM! tmp_int = 0; tmp_int(start:end) = fit_snap(1:i_fit_snap)%stress
      !debugMCM! call subworlds_allreduce_vect_int(tmp_int)
      !debugMCM! fit_snap(:)%stress = tmp_int
      tmp_log = .false.; tmp_log(start:end) = fit_snap(1:i_fit_snap)%stress
      call subworlds_allreduce_vect_logical(tmp_log)
      fit_snap(:)%stress = tmp_log


      tmp_int = 0; tmp_int(start:end) = fit_snap(1:i_fit_snap)%iatom
      call subworlds_allreduce_vect_int(tmp_int)
      fit_snap(:)%iatom = tmp_int
      tmp_int = 0; tmp_int(start:end) = fit_snap(1:i_fit_snap)%ix
      call subworlds_allreduce_vect_int(tmp_int)
      fit_snap(:)%ix = tmp_int
      deallocate(tmp_int)

      ! write(*,*) 'all', fit_snap(34853)%energy,fit_snap(34853)%force,fit_snap(34853)%stress  

      allocate(tmp_real(dim_design_column))
      tmp_real = 0.d0; tmp_real(start:end) = fit_snap(1:i_fit_snap)%weight
      call subworlds_allreduce_vect_double(tmp_real)
      fit_snap(:)%weight = tmp_real
      deallocate(tmp_real)

      ! updating i_fit_snap so that it matches the number of columns of global Amat
      i_fit_snap = SUM(i_fit_snap_on_MPI_group(:))

      deallocate(i_fit_snap_on_MPI_group); deallocate(i_fit_snap_offsets)

      timeb = MY_MPI_WTIME()
      time_fill_desc_02obj = time_fill_desc_02obj + (timeb - timea)

   end subroutine main_train_gather_mld

   subroutine main_train_mld()

      use ml_in_ndm_module, only: rangml, debug, mld_fit_type, snap_class_constraints, &
         fit_lapack_qr_constraints, tmp_val_desc_max, regularization_name,  &
         write_design_matrix, ml_type, ml_type_krr, train_time
      use module_db_setup, only: iconf_data   
      use temporary_data_cov, only: dim_data_constraints
      use derived_types, only: config_real
      use module_optimization, only: optimize_weights_db, optimize_weights_chem
      use mld_energy_mod, only: train_fill_Amat_with_energy
      use module_mld_quadratic, only: i_e_fit_snap, i_f_fit_snap, i_s_fit_snap

      ! TODOcvw use descriptors_interface, only: compute_descriptors

      use snap, only: i_fit_snap, i_constraints_snap

      use module_kind_variables, only: kind_double
      use time_check_general, only: time_fit_params, time_full_desc, time_calc_desc, &
         time_neigh_desc, time_fill_desc,  time_fill_desc_02, MY_MPI_WTIME, train_tot_time, tot_time, &
         debug_time, time_read_db, time_fill_desc_02A, time_fill_desc_02obj, time_fill_desc_02y
      use module_kernel, only: time_for_desc_kernel, train_time_for_desc_kernel
      use module_ml_scalapack, only: scalapack_driver
      use module_kernel_2b, only: activate_k2b, train_time_for_desc_kernel_2b, time_for_desc_kernel_2b, &
         tnn_2b, tuu_2b, tkk_2b, r_cut_2b, r_cut_width_2b
      use module_continuity_k2b_zbl, only: get_continuity_k2b_zbl
      use module_zbl, only: zbl_potential, zbl_type, zbl_mode_default_k2b, r1_zbl, rr_k2b
      use module_neigh_local, only: type_fcut
      use module_write_parameters, only: write_zbl_k2b, write_snap_parameters
      use module_pair_cutoffs, only: write_pair_cutoffs_xml
      use module_chemical_species, only: tnn_rdist, tcc01_rdist, tcc02_rdist, tcc03_rdist, tii_rdist
      use module_ftnbody, only: tnn_ftbd, t2b_ftbd, t3b_ftbd, t4b_ftbd, t5b_ftbd, tnn_ftbd_train, t2b_ftbd_train, t3b_ftbd_train, t4b_ftbd_train, t5b_ftbd_train
      use module_optimization, only: itopt
      use module_db_poscar, only: i_start_conf, i_final_conf, procs_per_file
      use set_limits, only: set_limit_for_configs_with_MPI_grid
      use mld_subworld
      use module_write_design_matrix, only: dump_design_matrix, dump_design_matrix_test
      use my_mpi_subroutines, only: subworlds_get_max_double
      integer  :: i
      real(kind(1.d0))     :: tmp_val, val_max 
      real(kind_double)    :: t0, t1, t2, t3, t4, t5, t6, t7
      logical :: post_desc

      _NAMECURRENT_("main_train_mld")
      _MLD_BEGIN_
      train_time=.true.

      !testing the training ...
      call prepare_database()

      !set-up on which proc is read each conf ...
      call init_subworld(procs_per_file)
      call set_limit_for_configs_with_MPI_grid(iconf_data, i_start_conf, i_final_conf)
      call prepare_train_dimensions()
      call train_allocate_mld()

      i_fit_snap = 0
      i_e_fit_snap = 0
      i_f_fit_snap = 0
      i_s_fit_snap = 0
      i_constraints_snap = 0

      !computing the descriptors for the training part ...
      tmp_val = -1.d0

      time_fit_params = 0.d0
      time_full_desc = 0.d0
      time_calc_desc = 0.d0
      time_neigh_desc = 0.d0
      time_fill_desc = 0.d0
      time_fill_desc_02 = 0.d0
      time_fill_desc_02y = 0.d0 
      time_fill_desc_02obj = 0.d0 
      time_fill_desc_02A = 0.d0 

      if (itopt == 0) time_read_db = 0.d0

      !$! if (train_covariance_matrix) then 
      !$!    call compute_train_covariace_matrix(i_start_conf, i_final_conf) 
      !$! end if 

      t0 = MY_MPI_WTIME()
      do i = i_start_conf, i_final_conf
         if (.not. (config_real(i)%train)) cycle
         !config_real(i)%has_force=.false.
         t1 = MY_MPI_WTIME()
         call test_if_config_is_small(i)
         !  config_real(i)%small=.true.
         call calc_neighbours(i)
         t2 = MY_MPI_WTIME()
         time_neigh_desc = time_neigh_desc + (t2 - t1)
         post_desc=.true. 
         call compute_descriptors(i, post_desc)
         t3 = MY_MPI_WTIME()
         time_calc_desc = time_calc_desc + (t3 - t2)
         call train_fill_Amat_with_energy(i)
         call train_fill_Amat_with_force(i)
         call train_fill_Amat_with_stress(i)
         if ((mld_fit_type == fit_lapack_qr_constraints) .and. config_real(i)%class == snap_class_constraints) then
            call train_fill_Bmat_constraints(i)
         end if
         if (tmp_val_desc_max >= tmp_val) then
            tmp_val = tmp_val_desc_max
         end if
         t4 = MY_MPI_WTIME()
         time_fill_desc = time_fill_desc + (t4 - t3)
         call train_deallocate_desc(i)
         !debug call status_allocate_config_desc(i)
         !herecos! call train_deallocate_real(i)
         !call status_allocate_config_real(i)
      end do

      !debug write(*,*) "MAX.....  ", mld_rank,  subworld, subrank, tmp_val_desc_max, tmp_val 
      call subworlds_get_max_double(tmp_val, val_max)

      t4 = MY_MPI_WTIME()
      call main_train_gather_mld()
      t5 = MY_MPI_WTIME()
      time_fill_desc_02 = time_fill_desc_02 + (t5 - t4)
      time_full_desc = time_full_desc + (t5 - t0)

      if (write_design_matrix) call dump_design_matrix

      ! if (ml_type == ml_type_descriptors) return 

      if ((.not. (optimize_weights_db.or.optimize_weights_chem)).or.(itopt==0)) then
         if (rangml == 0) then
            write (6, *) 'ML: the max value of the descriptor is   :', val_max 
         end if
      end if

      if (mld_fit_type == fit_lapack_qr_constraints) then
         if (i_constraints_snap /= dim_data_constraints) then
            write (6, *) 'FATAL: Serious problems in imposing constraints', i_constraints_snap, dim_data_constraints
            stop 'The dimension of constraints are ill defined'
         end if
         if ((dim_force_constraints > 0) .or. (dim_stress_constraints > 0)) then
            write (6, *) 'FATAL: the constraints of forces and stress is not, YET, implemented. If you want to continue '// &
               'put the corresponding flags from TRUE to FALSE'
            stop ' not force or stress constraints'
         end if
      end if

      t6 = MY_MPI_WTIME()
      ! TODO parallel fitting ToScaRemove this if ...
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

! TODO if sequential
#if(PARA)
      if (.not.scalapack_driver) then
         !TORC! call my_barrier_mld(codeml)
         call comm_mld%barrier()
         !TORC! call MPI_BCAST(w_params(:, 1), size(w_params, 1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, codeml)
         call comm_mld%bcast(0, w_params(:,1))
         !TORC! TODO Amat - here are obsolete lines ...
         !TORC!  do irows = 1, size(Amat, 2)
         !TORC!    call MPI_BCAST(Amat(:, irows), size(Amat, 1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, codeml)
         !TORC!  end do
      end if
#endif
      !TORC! call my_barrier_mld(codeml)
      call comm_mld%barrier()

      t7 = MY_MPI_WTIME()
      time_fit_params = time_fit_params + (t7 - t6)


      !TODO regularization ToScaRemove .. I do not understand why the following
      ! subroutine is executed on only one proc.
      if (scalapack_driver) then
         if (rangml==0) regularization_name = ""
         call get_train_errors
      else
         if (.not.optimize_weights_chem) then
            if (rangml == 0) then
               regularization_name = ""
               ! PAY ATTENTION THAT THIS IS also THE OPTIMIZATION ON WEIGTHS ... i DO NOT KONW why ....
               call get_train_errors
            end if                  ! rangml=>0 running on the proc 0 for trainning ... train 0milady ...
         else
            regularization_name = ""
            call get_train_errors
         end if
      end if


      if (ml_type == ml_type_krr) then
         train_time_for_desc_kernel = time_for_desc_kernel
         time_for_desc_kernel = 0.d0
      end if

      if (activate_k2b) then
         if (zbl_potential .and. zbl_type == zbl_mode_default_k2b) then
            call get_continuity_k2b_zbl(r1_zbl, rr_k2b, type_fcut, r_cut_2b, r_cut_width_2b)
            call write_zbl_k2b
            !zbl_debug! np=1000
            !zbl_debug! call linear_grid(np, 0.8d0,3.0d0,xp_grid)
            !zbl_debug! do ii = 1, np
            !zbl_debug! rr = xp_grid(ii)
            !zbl_debug! call full_potential_zbl(.true., 1, rr, ee, ff)
            !zbl_debug! write(324, '(3e25.10)') rr, ee, ff
            !zbl_debug! !if (rr <= rr_k2b) then
            !zbl_debug! call kernel_2b_end(1, 1, .true., type_fcut, rr, r_cut_2b, r_cut_width_2b, ee, ff, ss )
            !zbl_debug! write(325,'(4e25.10)') rr, ee, ff, ss
            !zbl_debug! !end if
            !zbl_debug! call potential_zbl(.true.,74.d0, 74.d0, rr, ee, ff)
            !zbl_debug! write(327, '(4e25.10)') rr, ee, ff, ss
            !zbl_debug! call potential_zbl_second(74.d0, 74.d0, rr, ee, ff, ss)
            !zbl_debug! write(326, '(4e25.10)') rr, ee, ff, ss
            !zbl_debug! !cont_test
            !zbl_debug! end do
         end if
         train_time_for_desc_kernel_2b = time_for_desc_kernel_2b
         time_for_desc_kernel_2b = 0.d0
         tnn_2b = 0.d0
         tuu_2b = 0.d0
         tkk_2b = 0.d0
      end if
      tnn_rdist = 0.d0
      tcc01_rdist = 0.d0
      tcc02_rdist = 0.d0
      tcc03_rdist = 0.d0
      tii_rdist = 0.d0

      ! Write pair-specific cutoff parameters to XML for LAMMPS
      call write_pair_cutoffs_xml()

      tnn_ftbd_train=tnn_ftbd
      tnn_ftbd=0.d0

      t2b_ftbd_train=t2b_ftbd
      t2b_ftbd=0.d0

      t3b_ftbd_train=t3b_ftbd
      t3b_ftbd=0.d0

      t4b_ftbd_train=t4b_ftbd
      t4b_ftbd=0.d0

      t5b_ftbd_train=t5b_ftbd
      t5b_ftbd=0.d0

      if (optimize_weights_chem) then
         call deallocate_fit
         call train_deallocate_mld
      end if

      if (debug_time) then
         train_tot_time(:) =  tot_time(:)
         tot_time = 0.d0
      end if

      call write_snap_parameters
      call write_pair_cutoffs_xml()
      call close_subworld()


      return 

      _MLD_END_

   end subroutine main_train_mld

   subroutine main_train_optimize_weights
      use ml_in_ndm_module, only: rangml, regularization_name
      use module_ml_scalapack, only: scalapack_driver
      use module_optimization, only: optimize_weights_db, optimize_weights_chem

      _NAMECURRENT_("main_train_optimize_weights")

      if (optimize_weights_db) then
         if (rangml == 0) then
            if (scalapack_driver) then
               write(6,*) 'optimize_weights not yet implemented with scalapack drivers'
               write(6,*) 'put scalapack_driver to .false. '
               stop
            end if
            regularization_name = ""
            if (rangml == 0) write (6, '("ML: the errors before optimization ...")')
            call get_train_errors

            !(r)evolutianary part ...
            call mld_optimize_weights

            if (rangml == 0) write (6, '("ML: the errors after  optimization ...")')
            !call get_test_errors
         end if ! rangml ==0
      end if

      if (optimize_weights_chem) then
         regularization_name = ""
         if (rangml == 0) write (6, '("ML: the errors before optimization ...")')
         !call get_train_errors
         call mld_optimize_weights
      end if
      _MLD_END_

   end subroutine main_train_optimize_weights





   subroutine get_fit_dimensions(iconf)
      ! increment the dimension of the fit depending on the configuration and if
      ! the energy forces or stress are included in the fit
      ! Input: iconf
      !        config_real(iconf)%
      !                          %has_stress has_force has_energy
      ! itest=1 -> train
      ! itest=2 -> test
      ! Ouptput:    incremented value of
      !             dim_data, dim_data_train, dim_data_test

      use ml_in_ndm_module, only: desc_forces
      use temporary_data_cov, only: dim_data_test, dim_data_train
      use derived_types, only: config_real
      !TODOcvw use snap, only: dim_ene_train_lml, dim_force_train_lml, dim_stress_train_lml, dim_ene_test_snap, dim_force_test_snap, dim_stress_test_snap

      integer, intent(in)  :: iconf

      _NAMECURRENT_("get_fit_dimensions")


      _MLD_BEGIN_

      if (config_real(iconf)%train) then
         if ((config_real(iconf)%has_force) .and. (desc_forces)) then
            dim_data_train = dim_data_train + 3*config_real(iconf)%nat
            dim_force_train_lml = dim_force_train_lml + 3*config_real(iconf)%nat
         end if

         if (config_real(iconf)%has_energy) then
            dim_data_train = dim_data_train + 1
            dim_ene_train_lml = dim_ene_train_lml + 1
         end if

         if ((config_real(iconf)%has_stress) .and. (desc_forces)) then
            dim_data_train = dim_data_train + 6
            dim_stress_train_lml = dim_stress_train_lml + 6
         end if
      end if

      if (.not. config_real(iconf)%train) then
         if ((config_real(iconf)%has_force) .and. (desc_forces)) then
            dim_data_test = dim_data_test + 3*config_real(iconf)%nat
            dim_force_test_snap = dim_force_test_snap + 3*config_real(iconf)%nat
         end if

         if (config_real(iconf)%has_energy) then
            dim_data_test = dim_data_test + 1
            dim_ene_test_snap = dim_ene_test_snap + 1
         end if

         if ((config_real(iconf)%has_stress) .and. (desc_forces)) then
            dim_data_test = dim_data_test + 6
            dim_stress_test_snap = dim_stress_test_snap + 6
         end if
      end if

      ! dim_data =  dim_data_test + dim_data_train
      _MLD_END_
   end subroutine get_fit_dimensions



   subroutine train_get_constraints_dimensions(iconf)
      ! increment the dimension of the constraints matrix Bmat and Bvec (from the constraints Bmax x X)
      ! the energy forces or stress are included in the fit
      ! Input: iconf
      !        config_real(iconf)%
      !                          %has_stress has_force has_energy
      ! Ouptput:    incremented value of
      !             dim_data_constraints

      use ml_in_ndm_module, only: rangml
      use temporary_data_cov, only: dim_data_constraints
      use derived_types, only: config_real
      !TODOcvw use snap, only: dim_ene_constraints, dim_force_constraints, dim_stress_constraints

      integer, intent(in)  :: iconf

      _NAMECURRENT_("train_get_constraints_dimensions")



      _MLD_BEGIN_
      if (config_real(iconf)%train) then
         if (config_real(iconf)%has_force) then
            dim_data_constraints = dim_data_constraints + 3*config_real(iconf)%nat
            dim_force_constraints = dim_force_constraints + 3*config_real(iconf)%nat
         end if
         if (config_real(iconf)%has_energy) then
            dim_data_constraints = dim_data_constraints + 1
            dim_ene_constraints = dim_ene_constraints + 1
         end if
         if (config_real(iconf)%has_stress) then
            dim_data_constraints = dim_data_constraints + 6
            dim_stress_constraints = dim_stress_constraints + 6
         end if
      end if

      if (dim_data_constraints == 0) then
         if (rangml == 0) then
            write (6, *) 'WARNING: despite the fact that mld_fit_type is with constraints  we do not found any valuable constraint to impose:'
            write (6, *) 'WARNING: Possibles sources:'
            write (6, *) 'WARNING:   - the classes for which the constraints are selected have no train configutations'
            write (6, *) 'WARNING:   - the configurations that have been chosen have F F F trainning flags '
         end if
      end if

      if (dim_force_constraints > 0) then
         if (rangml == 0) then
            write (6, *) 'WARNING: I not recommend at all to impose constraints on forces. Unless you really know what you do: goooood luck!'
         end if
      end if
      _MLD_END_
   end subroutine train_get_constraints_dimensions


   subroutine set_mld_dimension()
      ! set-up the dimension of lien and column of the design matrix.
      ! In milady the design matrix
      ! line=something_related_to_the_dimension_of_the_descritor x column=the_number_of_data

      use ml_in_ndm_module, only: mld_order, mld_kernel, mld_linear, mld_linear_extended, &
         mld_quadratic, mld_polyc, &
         polyc_n_poly, polyc_n_hermite, mld_type_quadratic, mld_type_quadratic_zaxa, mld_type_quadratic_ZX
      use temporary_data_cov, only: dim_data_train, dim_xdesc, dim_xdesc_patch

      use snap, only: dim_xdesc_full

      !TODOcvw use snap, only: dim_xdesc_linear, dim_design_line, dim_design_column
      use module_mld_quadratic, only: dim_xdesc_quadratic, dim_quadratic
      use module_mld_polyc, only: dim_xdesc_polyc, dim_polyc
      use module_mld_kernel, only: dim_xdesc_kernel
      use module_kernel_2b, only: dim_kernel_2b, activate_k2b
      use module_snap_nlinear, only: dim_xdesc_nlinear, dim_nlinear
      use module_kernel, only: dim_kernel
      use module_nlinear, only: order_nlinear


      _NAMECURRENT_("set_mld_dimension")



      _MLD_BEGIN_
      !TODOoffset
      dim_xdesc_linear = 1 + dim_xdesc
      dim_xdesc_patch = 0
      if (activate_k2b) dim_xdesc_patch = dim_xdesc_patch + dim_kernel_2b

      select case (mld_order)

       case (mld_linear, mld_linear_extended)
         dim_xdesc_nlinear = dim_xdesc_linear + order_nlinear*dim_xdesc
         dim_nlinear = order_nlinear*dim_xdesc
         dim_design_column = dim_data_train
         dim_design_line = dim_xdesc_linear + dim_xdesc_patch
         dim_xdesc_full = dim_xdesc_linear + dim_xdesc_patch

       case (mld_quadratic)
         dim_xdesc_quadratic = dim_xdesc_linear + dim_xdesc**2
         dim_quadratic = dim_xdesc**2
         if ((mld_type_quadratic == mld_type_quadratic_zaxa).or.(mld_type_quadratic == mld_type_quadratic_ZX))  then
            dim_xdesc_quadratic = dim_xdesc_linear + dim_xdesc*dim_kernel_2b
            dim_quadratic = dim_xdesc*dim_kernel_2b
         end if
         dim_design_column = dim_data_train
         dim_design_line = dim_xdesc_quadratic + dim_xdesc_patch
         dim_xdesc_full = dim_xdesc_quadratic + dim_xdesc_patch

       case (mld_polyc)
         if (polyc_n_poly == 1) dim_xdesc_polyc = 1 + polyc_n_hermite*dim_xdesc
         if (polyc_n_poly == 2) dim_xdesc_polyc = 1 + polyc_n_hermite*dim_xdesc + polyc_n_hermite*dim_xdesc**2
         if (polyc_n_poly == 3) dim_xdesc_polyc = 1 + polyc_n_hermite*dim_xdesc + polyc_n_hermite*dim_xdesc**2 &
            + polyc_n_hermite*dim_xdesc**3
         dim_polyc = dim_xdesc_polyc - dim_xdesc_linear
         dim_design_column = dim_data_train
         dim_design_line = dim_xdesc_polyc + dim_xdesc_patch
         dim_xdesc_full = dim_xdesc_polyc + dim_xdesc_patch
       case (mld_kernel)
         dim_xdesc_kernel = dim_xdesc_linear + dim_kernel
         dim_design_column = dim_data_train
         dim_design_line = dim_xdesc_kernel + dim_xdesc_patch
         dim_xdesc_full = dim_xdesc_kernel + dim_xdesc_patch
       case default
         call log_critical('ML: this potential - not implemented in the test part, mld_order is  '// vtoa(mld_order)//' in '//NAMECURRENT )
         stop 'FATAL in set_mld_order with mld_order'
      end select
      _MLD_END_
   end subroutine set_mld_dimension



   subroutine train_allocate_mld()
      ! Allocate the matrix needed for snap fitting
      ! used only for training

      use ml_in_ndm_module, only: mld_fit_type, fit_lapack_qr_constraints
      use temporary_data_cov, only: dim_data_train, dim_data_train_local, &
         yfunc_train, dim_xdesc, dim_data_constraints
      use module_ml_scalapack, only: scalapack_driver, &
         myrow, mycol, nprow, npcol, iam, context, &
         sca_Amat, desc_sca_Amat, nbr_Amat, nbc_Amat, &
         dimr_sca_Amat, dimc_sca_Amat, l_dimr_sca_Amat, l_dimc_sca_Amat, &
         sca_ymat, desc_sca_ymat, nbr_ymat, nbc_ymat, &
         dimr_sca_ymat, dimc_sca_ymat, l_dimr_sca_ymat, l_dimc_sca_ymat
      use module_scalapack_tools, only: allocation_scalapack_matrix
      use module_optimization, only: itopt
      !TODOcvw use snap, only: w_params, Amat, ymat, Bmat, zmat, fit_snap, weights_snap,dim_design_line, dim_design_column


      _NAMECURRENT_("train_allocate_mld")



      _MLD_BEGIN_
      if ((dim_xdesc == 0) .or. (dim_data_train == 0)) then
         call log_critical("The problem is not well posed dim_xdesc, dim_data_train"//vtoa(dim_xdesc)//" "//vtoa(dim_data_train))
         stop 'zeros in allocate_snap'
      end if

      call set_mld_dimension()

      ! In milady the design matrix  line=something_related_to_the_dimension_of_the_descritor x column=the_number_of_datza
      ! Fit related large matrix ... that should be distributed
      if (scalapack_driver) then
         !sca_Amat initialisation and allocation: generic short verion.
         nbr_Amat = -1
         nbc_Amat = -1
         dimr_sca_Amat = dim_design_line
         dimc_sca_Amat = dim_design_column
         call allocation_scalapack_matrix (" sca_Amat ", sca_Amat, desc_sca_Amat, dimr_sca_Amat, dimc_sca_Amat, &
            l_dimr_sca_Amat, l_dimc_sca_Amat, nbr_Amat, nbc_Amat, &
            myrow, mycol, nprow, npcol, context, iam  )

         ! ymat vector distribution:   generic version.
         !         |   ...   |
         !         |   ...   |
         !  ymat = |  M x 1  |  M is columns number of Amat
         !         |   ...   |
         !         |   ...   |
         nbr_ymat = nbc_Amat
         nbc_ymat = 1
         dimr_sca_ymat = dim_design_column
         dimc_sca_ymat = 1
         call allocation_scalapack_matrix (" sca_ymat ", sca_ymat, desc_sca_ymat, &
            dimr_sca_ymat, dimc_sca_ymat, &
            l_dimr_sca_ymat, l_dimc_sca_ymat, nbr_ymat, nbc_ymat, &
            myrow, mycol, nprow, npcol, context, iam  )

         ! Amat local to each subworld
         if (allocated(Amat)) deallocate (Amat); allocate (Amat(dim_design_line, dim_data_train_local))
      else
         if (allocated(Amat)) deallocate (Amat); allocate (Amat(dim_design_line, dim_design_column))
      end if

      ! Fit related small matrix than can be shared on all procs
      if (allocated(w_params)) deallocate (w_params); allocate (w_params(dim_design_line, 1))

      if (allocated(yfunc_train)) deallocate (yfunc_train); allocate (yfunc_train(dim_design_column))
      if (allocated(ymat)) deallocate (ymat); allocate (ymat(dim_design_column, 1))


      if (allocated(fit_snap)) deallocate (fit_snap); allocate (fit_snap(dim_design_column))
      if (allocated(weights_snap)) deallocate (weights_snap); allocate (weights_snap(dim_design_column))
      if (allocated(y_zbl_train)) deallocate (y_zbl_train); allocate (y_zbl_train(dim_design_column))
      y_zbl_train(:) = 0.d0

      ! Constraints related
      if (mld_fit_type == fit_lapack_qr_constraints) then
         if (allocated(Bmat)) deallocate (Bmat); allocate (Bmat(dim_design_line, dim_data_constraints))
         if (allocated(zmat)) deallocate (zmat); allocate (zmat(dim_data_constraints, 1))
      end if


      ! Infos utilities
      if (itopt == 0 ) then
         !cosw if (rangml == 0) write (6, '("ML: The number of fit  parameters are exactly: ", i9)') size(w_params, 1)
         call log_info("ML: The number of fit parameters are exactly:  "// vtoa(size(w_params, 1)))
      end if
      _MLD_END_
   end subroutine train_allocate_mld

   subroutine train_deallocate_mld()
      ! Allocate the matrix needed for snap fitting
      ! used only for training

      use ml_in_ndm_module, only: mld_fit_type, fit_lapack_qr_constraints
      use temporary_data_cov, only: yfunc_train
      use module_ml_scalapack, only: scalapack_driver, sca_Amat, desc_sca_Amat,  &
         sca_ymat, desc_sca_ymat
      use module_scalapack_tools, only: allocation_scalapack_matrix
      use derived_types, only: config_real_copy


      _NAMECURRENT_("train_deallocate_mld")

      _MLD_BEGIN_


      if (scalapack_driver) then
         if (allocated(sca_Amat)) deallocate(sca_Amat)
         if (allocated(desc_sca_Amat)) deallocate(desc_sca_Amat)

         if (allocated(sca_ymat)) deallocate(sca_ymat)
         if (allocated(desc_sca_ymat)) deallocate(desc_sca_ymat)
      end if

      if (allocated(Amat)) deallocate (Amat)

      if (allocated(yfunc_train)) deallocate (yfunc_train)
      if (allocated(ymat)) deallocate (ymat)


      if (allocated(fit_snap)) deallocate (fit_snap)
      if (allocated(weights_snap)) deallocate (weights_snap)
      if (allocated(y_zbl_train)) deallocate (y_zbl_train)

      ! Constraints related
      if (mld_fit_type == fit_lapack_qr_constraints) then
         if (allocated(Bmat)) deallocate (Bmat)
         if (allocated(zmat)) deallocate (zmat)
      end if

      !if (allocated(config_desc)) deallocate(config_desc)
      !if (allocated(config_real)) deallocate(config_real)
      if (allocated(config_real_copy)) deallocate(config_real_copy)

      _MLD_END_
   end subroutine train_deallocate_mld


   subroutine train_snap_compute_energy_force_stress(idata, einp, finp, sinp)

      use ml_in_ndm_module, only: rangml
      use temporary_data_cov, only: yfunc_train
      use derived_types, only: config_real
      use module_evaluate_parameters, only: product_w_params_Amat
      use module_zbl, only: zbl_potential, zbl_type, zbl_mode_alone

      integer, intent(in)  :: idata, einp, finp, sinp
      real(kind(0.d0))     :: force_snap, stress_snap, local_ref_energy, local_ref_zbl, &
         local_force_zbl, local_stress_zbl

      _NAMECURRENT_("train_snap_compute_energy_force_stress")



      _MLD_BEGIN_
      if (fit_snap(idata)%energy) then
         !ref_energy = config_real(fit_snap(idata)%iconf)%ref_energy
         i_e_train_snap = i_e_train_snap + 1
         !ene_snap = dot_product(w_params(:, 1), Amat(:, idata)) + ref_energy
         call product_w_params_Amat(ene_snap, idata)

         !DZBL
         local_ref_energy =  config_real(fit_snap(idata)%iconf)%ref_energy
         local_ref_zbl = 0.d0
         if (zbl_potential .and. zbl_type == zbl_mode_alone) then
            local_ref_zbl    = y_zbl_train(idata)
         end if
         ene_snap = ene_snap + local_ref_zbl + local_ref_energy
         if (rangml==0) write (einp, '(2es20.10," ", (a)," ",(a))') ene_snap, yfunc_train(idata) + local_ref_energy, config_real(fit_snap(idata)%iconf)%class, config_real(fit_snap(idata)%iconf)%filename
         !debug! if (rangml==0) write (777, '(a,5es20.10," ",(a)," ",(a))') 'E ', ene_snap, yfunc_train(idata) + local_ref_energy, &
         !debug!                                                   local_ref_zbl, &
         !debug!                                                   ene_snap - local_ref_zbl, &
         !debug!                                                   (yfunc_train(idata) + local_ref_energy) - local_ref_zbl, &
         !debug!                                                   config_real(fit_snap(idata)%iconf)%class, config_real(fit_snap(idata)%iconf)%filename

         y_e_train_snap(i_e_train_snap) = ene_snap
         y_e_train_base(i_e_train_snap) = yfunc_train(idata) + local_ref_energy
         y_e_p_a_train_snap(i_e_train_snap) = ene_snap/dble(config_real(fit_snap(idata)%iconf)%nat)
         y_e_p_a_train_base(i_e_train_snap) = (yfunc_train(idata) + local_ref_energy)/dble(config_real(fit_snap(idata)%iconf)%nat)
      end if

      if (fit_snap(idata)%force) then
         i_f_train_snap = i_f_train_snap + 1
         !force_snap = dot_product(w_params(:, 1), Amat(:, idata))
         call product_w_params_Amat(force_snap, idata)

         local_force_zbl = 0.d0
         if (zbl_potential .and. zbl_type == zbl_mode_alone) then
            local_force_zbl = y_zbl_train(idata)
         end if
         force_snap = force_snap + local_force_zbl
         if (rangml==0) write (finp, '(2es20.10," ",(a)," ",(a) )') force_snap, yfunc_train(idata), config_real(fit_snap(idata)%iconf)%class, config_real(fit_snap(idata)%iconf)%filename

         y_f_train_snap(i_f_train_snap) = force_snap
         y_f_train_base(i_f_train_snap) = yfunc_train(idata)
      end if

      if (fit_snap(idata)%stress) then
         i_s_train_snap = i_s_train_snap + 1
         !stress_snap = dot_product(w_params(:, 1), Amat(:, idata))
         call product_w_params_Amat(stress_snap, idata)

         local_stress_zbl = 0.d0
         if (zbl_potential .and. zbl_type == zbl_mode_alone) then
            local_stress_zbl = y_zbl_train(idata)
         end if
         stress_snap = stress_snap + local_stress_zbl
         if (rangml==0) write (sinp, '(2es20.10," ",(a), " ", (a))') stress_snap, yfunc_train(idata), config_real(fit_snap(idata)%iconf)%class, config_real(fit_snap(idata)%iconf)%filename
         !debug! if (rangml==0) write (777, '(a,2es20.10," ",(a)," ",(a))') 'S ', stress_snap, yfunc_train(idata), &
         !debug!                                                   config_real(fit_snap(idata)%iconf)%class, config_real(fit_snap(idata)%iconf)%filename
         y_s_train_snap(i_s_train_snap) = stress_snap
         y_s_train_base(i_s_train_snap) = yfunc_train(idata) + local_stress_zbl
      end if
      _MLD_END_
   end subroutine train_snap_compute_energy_force_stress



   subroutine train_error_mld()
      use module_optimization, only: itopt, optimize_weights_chem, optimize_weights_db
      use ml_in_ndm_module, only: rangml


      real(kind(0.d0))     :: corr, detr, rmse, mae

      _NAMECURRENT_("train_error_mld")



      _MLD_BEGIN_
      if (dim_ene_train_lml > 0) then
         call correlation_coef(y_e_train_base, y_e_train_snap, dim_ene_train_lml, corr)
         call determination_coef(y_e_train_base, y_e_train_snap, dim_ene_train_lml, detr)
         call rmse_mae(y_e_train_base, y_e_train_snap, dim_ene_train_lml, rmse, mae)
         if ((.not. (optimize_weights_db.or.optimize_weights_chem)).or.(itopt==0)) then
            if (rangml == 0) then
               write (6, '("ML: stat  Total_Energy            train:", 4f15.6)') corr, detr, rmse, mae
            end if
         end if

         call correlation_coef(y_e_p_a_train_base, y_e_p_a_train_snap, dim_ene_train_lml, corr)
         call determination_coef(y_e_p_a_train_base, y_e_p_a_train_snap, dim_ene_train_lml, detr)
         call rmse_mae(y_e_p_a_train_base, y_e_p_a_train_snap, dim_ene_train_lml, rmse, mae)
         if ((.not. (optimize_weights_db.or.optimize_weights_chem)).or.(itopt==0)) then
            if (rangml == 0) then
               write (6, '("ML: stat  Total_Energy_per_atom   train:", 4f15.6)') corr, detr, rmse, mae
            end if
         end if
         train_rmse_energy = rmse
         train_mae_energy = mae
      end if

      if (dim_force_train_lml > 0) then
         call correlation_coef(y_f_train_base, y_f_train_snap, dim_force_train_lml, corr)
         call determination_coef(y_f_train_base, y_f_train_snap, dim_force_train_lml, detr)
         call rmse_mae(y_f_train_base, y_f_train_snap, dim_force_train_lml, rmse, mae)
         if ((.not. (optimize_weights_db.or.optimize_weights_chem)).or.(itopt==0)) then
            if (rangml == 0) then
               write (6, '("ML: stat  Force_each_component    train:", 4f15.6)') corr, detr, rmse, mae
            end if
         end if
         train_rmse_force = rmse
         train_mae_force = mae
      end if

      if (dim_stress_train_lml > 0) then
         call correlation_coef(y_s_train_base, y_s_train_snap, dim_stress_train_lml, corr)
         call determination_coef(y_s_train_base, y_s_train_snap, dim_stress_train_lml, detr)
         call rmse_mae(y_s_train_base, y_s_train_snap, dim_stress_train_lml, rmse, mae)
         if ((.not. (optimize_weights_db.or.optimize_weights_chem)).or.(itopt==0)) then
            if (rangml == 0) then
               write (6, '("ML: stat  Stress_virial           train:", 4f15.6)') corr, detr, rmse, mae
            end if
         end if
         train_rmse_stress = rmse
         train_mae_stress = mae
      end if
      _MLD_END_
   end subroutine train_error_mld



   subroutine train_error_by_class()
      use module_kind_variables, only: kind_double
      use ml_in_ndm_module, only: rangml,  classes_full_for_sigma
      use module_error_by_class, only: class_if_ene, class_if_for, class_if_str, &
         rmse_by_class_ene, mae_by_class_ene, &
         corr_by_class_ene, detr_by_class_ene, &
         rmse_by_class_for, mae_by_class_for, &
         corr_by_class_for, detr_by_class_for, &
         rmse_by_class_str, mae_by_class_str, &
         corr_by_class_str, detr_by_class_str, &
         number_of_distinct_classes

      use temporary_data_cov, only: dim_data_train, yfunc_train
      use derived_types, only : config_real
      use module_evaluate_parameters, only: product_w_params_Amat
      use module_optimization, only: itopt, optimize_weights_db, optimize_weights_chem
      use module_zbl, only: zbl_potential, zbl_type, zbl_mode_alone

      implicit none

      real(kind_double)  :: corr, detr, rmse, mae, tmp_ref_ene, &
         tmp_ene, tmp_for, tmp_str, tmp_ref_zbl
      character(len=2)  :: ctmp, class_tmp
      integer  :: ii, i_ene, i_for, i_str, iconfig, iclass, idata
      real(kind_double), allocatable, dimension(:) :: y_train, y_base, y_p_a_base, y_p_a_train
      _NAMECURRENT_("train_error_by_class")


      _MLD_BEGIN_
      if (allocated(class_if_ene)) deallocate(class_if_ene) ; allocate(class_if_ene(number_of_distinct_classes))
      if (allocated(class_if_for)) deallocate(class_if_for) ; allocate(class_if_for(number_of_distinct_classes))
      if (allocated(class_if_str)) deallocate(class_if_str) ; allocate(class_if_str(number_of_distinct_classes))
      class_if_ene(:) = 0
      class_if_for(:) = 0
      class_if_str(:) = 0
      if (allocated(rmse_by_class_ene)) deallocate(rmse_by_class_ene) ; allocate(rmse_by_class_ene(number_of_distinct_classes))
      if (allocated(mae_by_class_ene)) deallocate(mae_by_class_ene) ; allocate(mae_by_class_ene(number_of_distinct_classes))
      if (allocated(corr_by_class_ene)) deallocate(corr_by_class_ene) ; allocate(corr_by_class_ene(number_of_distinct_classes))
      if (allocated(detr_by_class_ene)) deallocate(detr_by_class_ene) ; allocate(detr_by_class_ene(number_of_distinct_classes))
      if (allocated(rmse_by_class_for)) deallocate(rmse_by_class_for) ; allocate(rmse_by_class_for(number_of_distinct_classes))
      if (allocated(mae_by_class_for)) deallocate(mae_by_class_for) ; allocate(mae_by_class_for(number_of_distinct_classes))
      if (allocated(corr_by_class_for)) deallocate(corr_by_class_for) ; allocate(corr_by_class_for(number_of_distinct_classes))
      if (allocated(detr_by_class_for)) deallocate(detr_by_class_for) ; allocate(detr_by_class_for(number_of_distinct_classes))
      if (allocated(rmse_by_class_str)) deallocate(rmse_by_class_str) ; allocate(rmse_by_class_str(number_of_distinct_classes))
      if (allocated(mae_by_class_str)) deallocate(mae_by_class_str) ; allocate(mae_by_class_str(number_of_distinct_classes))
      if (allocated(corr_by_class_str)) deallocate(corr_by_class_str) ; allocate(corr_by_class_str(number_of_distinct_classes))
      if (allocated(detr_by_class_str)) deallocate(detr_by_class_str) ; allocate(detr_by_class_str(number_of_distinct_classes))
      if ((.not. (optimize_weights_db.or.optimize_weights_chem)).or.(itopt==0)) then
         if (rangml == 0) write (*, '("ML:|---------------------------- Errors by class -----------------------------|")')
      end if
      ! count if energy in this class ...
      do  iclass = 1, number_of_distinct_classes
         ctmp = classes_full_for_sigma(iclass)
         ! here we will deal with energy ...
         i_ene = 0
         i_for = 0
         i_str = 0
         do idata = 1, dim_data_train
            iconfig = fit_snap(idata)%iconf
            class_tmp = config_real(iconfig)%class
            if (ctmp /= class_tmp ) cycle
            if (fit_snap(idata)%energy) then
               i_ene = i_ene + 1
            end if
            if (fit_snap(idata)%force) then
               i_for = i_for + 1
            end if
            if (fit_snap(idata)%stress) then
               i_str = i_str + 1
            end if
         end do

         ! ... energy by class ....
         if (i_ene > 0 ) then
            class_if_ene(iclass) = 1
            if (allocated(y_train)) deallocate(y_train) ; allocate(y_train(i_ene))
            if (allocated(y_base)) deallocate(y_base) ; allocate(y_base(i_ene))
            if (allocated(y_p_a_train)) deallocate(y_p_a_train) ; allocate(y_p_a_train(i_ene))
            if (allocated(y_p_a_base)) deallocate(y_p_a_base) ; allocate(y_p_a_base(i_ene))
            ii = 0
            do idata = 1, dim_data_train
               iconfig = fit_snap(idata)%iconf
               class_tmp = config_real(iconfig)%class
               if (ctmp /= class_tmp ) cycle
               if (fit_snap(idata)%energy) then
                  tmp_ref_ene = config_real(iconfig)%ref_energy
                  tmp_ref_zbl = 0.d0
                  if (zbl_potential .and. zbl_type == zbl_mode_alone) tmp_ref_zbl = config_real(iconfig)%ezbl
                  ii = ii + 1
                  call product_w_params_Amat(tmp_ene, idata)
                  tmp_ene = tmp_ene + tmp_ref_ene
                  if (zbl_potential .and. zbl_type == zbl_mode_alone) then
                     tmp_ene = tmp_ene + tmp_ref_zbl
                  end if
                  !write (einp, '(2es20.10," ", (a)," ",(a))') ene_snap, yfunc_train(idata) + ref_energy, config_real(fit_snap(idata)%iconf)%class, config_real(fit_snap(idata)%iconf)%filename

                  y_train(ii) = tmp_ene
                  y_base(ii) = yfunc_train(idata) + tmp_ref_ene
                  y_p_a_train (ii) = tmp_ene/dble(config_real(iconfig)%nat)
                  y_p_a_base(ii) = (yfunc_train(idata) + tmp_ref_ene)/dble(config_real(iconfig)%nat)
               end if
            end do ! idata

            call correlation_coef(y_base, y_train, i_ene, corr)
            call determination_coef(y_base, y_train, i_ene, detr)
            call rmse_mae(y_base, y_train, i_ene, rmse, mae)
            rmse_by_class_ene(iclass) = rmse
            mae_by_class_ene(iclass) = mae
            corr_by_class_ene(iclass) = corr
            detr_by_class_ene(iclass) = detr
            if ((.not. (optimize_weights_db.or.optimize_weights_chem)).or.(itopt==0)) then
               if (rangml==0) write(6,'("ML: class Total_Energy             ",(a), "  :", 4f15.6)') ctmp, corr_by_class_ene(iclass),   detr_by_class_ene(iclass),  rmse_by_class_ene(iclass), mae_by_class_ene(iclass)
            end if
            call correlation_coef(y_p_a_base, y_p_a_train, i_ene, corr)
            call determination_coef(y_p_a_base, y_p_a_train, i_ene, detr)
            call rmse_mae(y_p_a_base, y_p_a_train, i_ene, rmse, mae)
            if ((.not. (optimize_weights_db.or.optimize_weights_chem)).or.(itopt==0)) then
               if (rangml==0) write(6,'("ML: class Total_Energy_per_atom    ",(a), "  :", 4f15.6)') ctmp, corr,   detr,  rmse, mae
            end if
         else
            if ((.not. (optimize_weights_db.or.optimize_weights_chem)).or.(itopt==0)) then
               if (rangml==0) write(6,'("ML: class Energy                   ",(a), "  :", 4(a))') ctmp, "       --       "
               if (rangml==0) write(6,'("ML: class Total_Energy_per_atom    ",(a), "  :", 4(a))') ctmp, "       --       "
            end if
         end if

         ! ... forces by class ....
         if (i_for > 0 ) then
            class_if_for(iclass) = 1
            if (allocated(y_train)) deallocate(y_train) ; allocate(y_train(i_for))
            if (allocated(y_base)) deallocate(y_base) ; allocate(y_base(i_for))
            ii = 0
            do idata = 1, dim_data_train
               iconfig = fit_snap(idata)%iconf
               class_tmp = config_real(iconfig)%class
               if (ctmp /= class_tmp ) cycle
               if (fit_snap(idata)%force) then
                  ii = ii + 1
                  !tmp_for = dot_product(w_params(:, 1), Amat(:, idata))
                  call product_w_params_Amat(tmp_for, idata)

                  y_train(ii) = tmp_for
                  y_base(ii) = yfunc_train(idata)
               end if
            end do ! idata

            call correlation_coef(y_base, y_train, i_for, corr)
            call determination_coef(y_base, y_train, i_for, detr)
            call rmse_mae(y_base, y_train, i_for, rmse, mae)
            rmse_by_class_for(iclass) = rmse
            mae_by_class_for(iclass) = mae
            corr_by_class_for(iclass) = corr
            detr_by_class_for(iclass) = detr
            if ((.not. (optimize_weights_db.or.optimize_weights_chem)).or.(itopt==0)) then
               if (rangml==0) write(6,'("ML: class Force_each_component     ",(a), "  :", 4f15.6)') ctmp, corr_by_class_for(iclass),   detr_by_class_for(iclass),  rmse_by_class_for(iclass), mae_by_class_for(iclass)
            end if
         else
            if ((.not. (optimize_weights_db.or.optimize_weights_chem)).or.(itopt==0)) then
               if (rangml==0) write(6,'("ML: class Force_each_component     ",(a), "  :", 4(a))') ctmp, "       --       "
            end if
         end if

         ! ... stress by class ....
         if (i_str > 0 ) then
            class_if_for(iclass) = 1
            if (allocated(y_train)) deallocate(y_train) ; allocate(y_train(i_str))
            if (allocated(y_base)) deallocate(y_base) ; allocate(y_base(i_str))
            ii = 0
            do idata = 1, dim_data_train
               iconfig = fit_snap(idata)%iconf
               class_tmp = config_real(iconfig)%class
               if (ctmp /= class_tmp ) cycle
               if (fit_snap(idata)%stress) then
                  ii = ii + 1
                  !tmp_str = dot_product(w_params(:, 1), Amat(:, idata))
                  call product_w_params_Amat(tmp_str, idata)

                  y_train(ii) = tmp_str
                  y_base(ii) = yfunc_train(idata)
               end if
            end do ! idata

            call correlation_coef(y_base, y_train, i_str, corr)
            call determination_coef(y_base, y_train, i_str, detr)
            call rmse_mae(y_base, y_train, i_str, rmse, mae)
            rmse_by_class_str(iclass) = rmse
            mae_by_class_str(iclass) = mae
            corr_by_class_str(iclass) = corr
            detr_by_class_str(iclass) = detr
            if ((.not. (optimize_weights_db.or.optimize_weights_chem)).or.(itopt==0)) then
               if (rangml==0) write(6,'("ML: class Stress_virial            ",(a), "  :", 4f15.6)') ctmp, corr_by_class_str(iclass),   detr_by_class_str(iclass),  rmse_by_class_str(iclass), mae_by_class_str(iclass)
            end if
         else
            if ((.not. (optimize_weights_db.or.optimize_weights_chem)).or.(itopt==0)) then
               if (rangml==0) write(6,'("ML: class Stress_virial            ",(a), "  :", 4(a))') ctmp, "       --       "
            end if
         end if
         if ((.not. (optimize_weights_db.or.optimize_weights_chem)).or.(itopt==0)) then
            if (rangml==0) write(6,'(" -- ")')
         end if
      end do ! end of classes

      _MLD_END_
   end subroutine train_error_by_class




   subroutine test_error_snap()
      use ml_in_ndm_module, only: rangml
      use snap, only: dim_ene_test_snap, dim_force_test_snap, dim_stress_test_snap, &
         y_e_test_snap, y_f_test_snap, y_s_test_snap, &
         y_e_test_base, y_f_test_base, y_s_test_base, &
         y_e_p_a_test_snap, y_e_p_a_test_base, &
         test_rmse_energy, test_rmse_force, test_rmse_stress, &
         test_mae_energy, test_mae_force, test_mae_stress
      use module_optimization, only: optimize_weights_db, optimize_weights_chem

      implicit none

      real(kind(0.d0))     :: corr, detr, rmse, mae

      _NAMECURRENT_("test_error_snap")



      _MLD_BEGIN_
      if (dim_ene_test_snap > 0) then
         call correlation_coef(y_e_test_base, y_e_test_snap, dim_ene_test_snap, corr)
         call determination_coef(y_e_test_base, y_e_test_snap, dim_ene_test_snap, detr)
         call rmse_mae(y_e_test_base, y_e_test_snap, dim_ene_test_snap, rmse, mae)
         if (.not. (optimize_weights_db.or.optimize_weights_chem)) then
            if (rangml == 0) then
               write (6, '("ML: stat  Total_Energy             test:", 4f15.6)') corr, detr, rmse, mae
            end if
         end if

         call correlation_coef(y_e_p_a_test_base, y_e_p_a_test_snap, dim_ene_test_snap, corr)
         call determination_coef(y_e_p_a_test_base, y_e_p_a_test_snap, dim_ene_test_snap, detr)
         call rmse_mae(y_e_p_a_test_base, y_e_p_a_test_snap, dim_ene_test_snap, rmse, mae)
         if (.not. (optimize_weights_db.or.optimize_weights_chem)) then
            if (rangml == 0) then
               write (6, '("ML: stat  Total_Energy_per_atom    test:", 4f15.6)') corr, detr, rmse, mae
            end if
         end if
         test_rmse_energy = rmse
         test_mae_energy = mae
      end if

      if (dim_force_test_snap > 0) then

         call correlation_coef(y_f_test_base, y_f_test_snap, dim_force_test_snap, corr)
         call determination_coef(y_f_test_base, y_f_test_snap, dim_force_test_snap, detr)
         call rmse_mae(y_f_test_base, y_f_test_snap, dim_force_test_snap, rmse, mae)
         if (.not. (optimize_weights_db.or.optimize_weights_chem)) then
            if (rangml == 0) then
               write (6, '("ML: stat  Force_each_component     test:", 4f15.6)') corr, detr, rmse, mae
            end if
         end if
         test_rmse_force = rmse
         test_mae_force = mae

      end if

      if (dim_stress_test_snap > 0) then

         call correlation_coef(y_s_test_base, y_s_test_snap, dim_stress_test_snap, corr)
         call determination_coef(y_s_test_base, y_s_test_snap, dim_stress_test_snap, detr)
         call rmse_mae(y_s_test_base, y_s_test_snap, dim_stress_test_snap, rmse, mae)
         if (.not. (optimize_weights_db.or.optimize_weights_chem)) then
            if (rangml == 0) then
               write (6, '("ML: stat  Stress_virial            test:", 4f15.6)') corr, detr, rmse, mae
            end if
         end if
         test_rmse_stress = rmse
         test_mae_stress = mae

      end if
      !debug write (*,*) dim_data, dim_data_test,  dim_data_train
      _MLD_END_
   end subroutine test_error_snap


   subroutine test_error_by_class()
      use module_kind_variables, only: kind_double
      use ml_in_ndm_module, only: rangml, classes_full_for_sigma, desc_forces
      use module_db_setup, only: iconf_data   
      use module_error_by_class, only: class_if_ene, class_if_for, class_if_str, &
         rmse_by_class_ene, mae_by_class_ene, &
         corr_by_class_ene, detr_by_class_ene, &
         rmse_by_class_for, mae_by_class_for, &
         corr_by_class_for, detr_by_class_for, &
         rmse_by_class_str, mae_by_class_str, &
         corr_by_class_str, detr_by_class_str, &
         number_of_distinct_classes

      use snap, only: y_e_test_snap, y_f_test_snap, y_s_test_snap, &
         y_e_test_base, y_f_test_base, y_s_test_base, &
         y_e_p_a_test_snap, y_e_p_a_test_base
      use derived_types, only : config_real
      use module_optimization, only: optimize_weights_db, optimize_weights_chem
      implicit none

      real(kind_double)  :: corr, detr, rmse, mae
      character(len=2)  :: ctmp, class_tmp
      integer  :: ii, i_ene, i_for, i_str, ivec_e, ivec_f, ivec_s, iconfig, iclass, nattmp, ix
      real(kind_double), allocatable, dimension(:) :: y_train, y_base, y_p_a_base, y_p_a_train
      _NAMECURRENT_("train_error_by_class")


      _MLD_BEGIN_
      if (allocated(class_if_ene)) deallocate(class_if_ene) ; allocate(class_if_ene(number_of_distinct_classes))
      if (allocated(class_if_for)) deallocate(class_if_for) ; allocate(class_if_for(number_of_distinct_classes))
      if (allocated(class_if_str)) deallocate(class_if_str) ; allocate(class_if_str(number_of_distinct_classes))
      class_if_ene(:) = 0
      class_if_for(:) = 0
      class_if_str(:) = 0
      if (allocated(rmse_by_class_ene)) deallocate(rmse_by_class_ene) ; allocate(rmse_by_class_ene(number_of_distinct_classes))
      if (allocated(mae_by_class_ene)) deallocate(mae_by_class_ene) ; allocate(mae_by_class_ene(number_of_distinct_classes))
      if (allocated(corr_by_class_ene)) deallocate(corr_by_class_ene) ; allocate(corr_by_class_ene(number_of_distinct_classes))
      if (allocated(detr_by_class_ene)) deallocate(detr_by_class_ene) ; allocate(detr_by_class_ene(number_of_distinct_classes))
      if (allocated(rmse_by_class_for)) deallocate(rmse_by_class_for) ; allocate(rmse_by_class_for(number_of_distinct_classes))
      if (allocated(mae_by_class_for)) deallocate(mae_by_class_for) ; allocate(mae_by_class_for(number_of_distinct_classes))
      if (allocated(corr_by_class_for)) deallocate(corr_by_class_for) ; allocate(corr_by_class_for(number_of_distinct_classes))
      if (allocated(detr_by_class_for)) deallocate(detr_by_class_for) ; allocate(detr_by_class_for(number_of_distinct_classes))
      if (allocated(rmse_by_class_str)) deallocate(rmse_by_class_str) ; allocate(rmse_by_class_str(number_of_distinct_classes))
      if (allocated(mae_by_class_str)) deallocate(mae_by_class_str) ; allocate(mae_by_class_str(number_of_distinct_classes))
      if (allocated(corr_by_class_str)) deallocate(corr_by_class_str) ; allocate(corr_by_class_str(number_of_distinct_classes))
      if (allocated(detr_by_class_str)) deallocate(detr_by_class_str) ; allocate(detr_by_class_str(number_of_distinct_classes))
      if (.not. (optimize_weights_db.or.optimize_weights_chem)) then
         if (rangml == 0) write (*, '("ML:|---------------------------- Errors by class -----------------------------|")')
      end if
      ! count if energy in this class ...
      do  iclass = 1, number_of_distinct_classes
         ctmp = classes_full_for_sigma(iclass)
         ! here we will deal with energy ...
         i_ene = 0
         i_for = 0
         i_str = 0
         ivec_e = 0
         !TODO_subworld_problem 
         do iconfig = 1, iconf_data
            if(config_real(iconfig)%train) cycle
            class_tmp = config_real(iconfig)%class
            if (ctmp /= class_tmp ) cycle
            if (config_real(iconfig)%has_energy) then
               i_ene = i_ene + 1
            end if
            if (config_real(iconfig)%has_force .and. desc_forces) then
               i_for = i_for + 3*config_real(iconfig)%nat
            end if
            if ((config_real(iconfig)%has_stress).and.desc_forces) then
               i_str = i_str + 6
            end if
         end do
         ! ... energy by class ....
         if (i_ene > 0 ) then
            class_if_ene(iclass) = 1
            if (allocated(y_train)) deallocate(y_train) ; allocate(y_train(i_ene))
            if (allocated(y_base)) deallocate(y_base) ; allocate(y_base(i_ene))
            if (allocated(y_p_a_train)) deallocate(y_p_a_train) ; allocate(y_p_a_train(i_ene))
            if (allocated(y_p_a_base)) deallocate(y_p_a_base) ; allocate(y_p_a_base(i_ene))
            ii = 0
            do iconfig = 1, iconf_data
               if(config_real(iconfig)%train) cycle
               if (config_real(iconfig)%has_energy) then
                  ivec_e = ivec_e + 1
               end if
               class_tmp = config_real(iconfig)%class
               if (ctmp /= class_tmp ) cycle
               if (config_real(iconfig)%has_energy) then
                  ii = ii + 1
                  y_train(ii) = y_e_test_snap(ivec_e)
                  y_base(ii) = y_e_test_base(ivec_e)
                  y_p_a_train (ii) = y_e_p_a_test_snap(ivec_e)
                  y_p_a_base(ii) = y_e_p_a_test_base(ivec_e)
               end if
            end do ! idata

            call correlation_coef(y_base, y_train, i_ene, corr)
            call determination_coef(y_base, y_train, i_ene, detr)
            call rmse_mae(y_base, y_train, i_ene, rmse, mae)
            rmse_by_class_ene(iclass) = rmse
            mae_by_class_ene(iclass) = mae
            corr_by_class_ene(iclass) = corr
            detr_by_class_ene(iclass) = detr

            if (rangml==0) write(6,'("ML: class Total_Energy             ",(a), "  :", 4f15.6)') ctmp, corr_by_class_ene(iclass),   detr_by_class_ene(iclass),  rmse_by_class_ene(iclass), mae_by_class_ene(iclass)
            call correlation_coef(y_p_a_base, y_p_a_train, i_ene, corr)
            call determination_coef(y_p_a_base, y_p_a_train, i_ene, detr)
            call rmse_mae(y_p_a_base, y_p_a_train, i_ene, rmse, mae)
            if (rangml==0) write(6,'("ML: class Total_Energy_per_atom    ",(a), "  :", 4f15.6)') ctmp, corr,   detr,  rmse, mae
         else
            if (rangml==0) write(6,'("ML: class Energy                   ",(a), "  :", 4(a))') ctmp, "       --       "
            if (rangml==0) write(6,'("ML: class Total_Energy_per_atom    ",(a), "  :", 4(a))') ctmp, "       --       "
         end if

         ! ... forces by class ....
         ivec_f = 0
         if (i_for > 0 ) then
            class_if_for(iclass) = 1
            if (allocated(y_train)) deallocate(y_train) ; allocate(y_train(i_for))
            if (allocated(y_base)) deallocate(y_base) ; allocate(y_base(i_for))
            ii = 0
            do iconfig = 1, iconf_data
               if (config_real(iconfig)%train) cycle
               class_tmp = config_real(iconfig)%class
               if (ctmp == class_tmp ) then
                  if ((config_real(iconfig)%has_force) .and. desc_forces) then
                     nattmp = config_real(iconfig)%nat
                     do ix = 1, 3
                        y_train(ii + (ix - 1)*nattmp + 1:ii + ix*nattmp)  = y_f_test_snap(ivec_f + (ix - 1)*nattmp + 1 : ivec_f + ix*nattmp )
                        y_base(ii + (ix - 1)*nattmp + 1:ii + ix*nattmp) = y_f_test_base(ivec_f + (ix - 1)*nattmp + 1 : ivec_f + ix*nattmp )
                     end do

                     ii = ii + 3 * config_real(iconfig)%nat
                  end if !class

                  if ((config_real(iconfig)%has_force) .and. desc_forces) then
                     ivec_f = ivec_f + 3 * config_real(iconfig)%nat
                  end if
               end if  ! error on forces

            end do ! idata

            call correlation_coef(y_base, y_train, i_for, corr)
            call determination_coef(y_base, y_train, i_for, detr)
            call rmse_mae(y_base, y_train, i_for, rmse, mae)
            rmse_by_class_for(iclass) = rmse
            mae_by_class_for(iclass) = mae
            corr_by_class_for(iclass) = corr
            detr_by_class_for(iclass) = detr

            if (rangml==0) write(6,'("ML: class Force_each_component     ",(a), "  :", 4f15.6)') ctmp, corr_by_class_for(iclass),   detr_by_class_for(iclass),  rmse_by_class_for(iclass), mae_by_class_for(iclass)
         else
            if (rangml==0) write(6,'("ML: class Force_each_component     ",(a), "  :", 4(a))') ctmp, "       --       "
         end if

         ! ... stress by class ....
         ivec_s = 0
         if (i_str > 0 ) then
            class_if_for(iclass) = 1
            if (allocated(y_train)) deallocate(y_train) ; allocate(y_train(i_str))
            if (allocated(y_base)) deallocate(y_base) ; allocate(y_base(i_str))
            ii = 0
            do iconfig = 1, iconf_data
               if (config_real(iconfig)%train) cycle
               class_tmp = config_real(iconfig)%class
               if (ctmp == class_tmp ) then
                  if ((config_real(iconfig)%has_stress) .and. desc_forces) then
                     y_train( ii + 1 : ii + 6 )  = y_s_test_snap(ivec_s + 1 : ivec_s + 6 )
                     y_base( ii +  1 : ii + 6 ) =    y_s_test_base(ivec_s + 1 : ivec_s + 6)

                     ii = ii + 6
                  end if
               end if !class

               if (config_real(iconfig)%has_stress .and. desc_forces) then
                  ivec_s = ivec_s + 6
               end if
            end do ! iconf stress

            call correlation_coef(y_base, y_train, i_str, corr)
            call determination_coef(y_base, y_train, i_str, detr)
            call rmse_mae(y_base, y_train, i_str, rmse, mae)
            rmse_by_class_str(iclass) = rmse
            mae_by_class_str(iclass) = mae
            corr_by_class_str(iclass) = corr
            detr_by_class_str(iclass) = detr

            if (rangml==0) write(6,'("ML: class Stress_virial            ",(a), "  :", 4f15.6)') ctmp, corr_by_class_str(iclass),   detr_by_class_str(iclass),  rmse_by_class_str(iclass), mae_by_class_str(iclass)
         else
            if (rangml==0) write(6,'("ML: class Stress_virial            ",(a), "  :", 4(a))') ctmp, "       --       "
         end if
         if (.not. (optimize_weights_db.or.optimize_weights_chem)) then
            if (rangml==0) write(6,'(" -- ")')
         end if
      end do ! end of classes

      _MLD_END_
   end subroutine test_error_by_class


   !/------------------------------------------------------------\
   !                                                             !
   !                           MD   PART                         !
   !                                                             !
   !\------------------------------------------------------------/



   subroutine md_mld_compute_energy(iconf, pack_opt)

      use module_kind_variables, only: kind_double
      use ml_in_ndm_module, only: mld_order, mld_linear, mld_quadratic, mld_polyc, mld_kernel
      use derived_types, only: config_desc, config_real
      use module_mld_quadratic, only: dim_xdesc_quadratic
      use module_mld_polyc, only: dim_xdesc_polyc
      use module_mld_kernel, only: dim_xdesc_kernel
      use module_kernel_2b, only: activate_k2b, dim_kernel_2b
      use module_zbl, only: zbl_potential, zbl_type, zbl_mode_alone
      use snap, only: dim_xdesc_full, xdesc_emd
      use mld_logger, only: log_critical, vtoa
      use mld_mpi, only: mld_mpi_abort

      integer, intent(in)  :: iconf
      logical, intent(in), optional    :: pack_opt
      logical  :: lpack_tmp
      integer  :: ii
      real(kind_double)    :: tmp_ene, local_ref_energy, local_ref_zbl
      ! real(kind_double), dimension(:), allocatable :: xdesc_emd
      integer :: dim_full_partial

      _NAMECURRENT_("md_mld_compute_energy")



      _MLD_BEGIN_
      lpack_tmp = .true.
      if (present(pack_opt)) lpack_tmp = pack_opt

      if (lpack_tmp) call pack_energy_descriptor(iconf)


      tmp_ene = 0.d0
      do ii = 1, config_real(iconf)%nat
         tmp_ene = tmp_ene + config_real(iconf)%ref_energy_per_element(config_real(iconf)%itype(ii))
      end do

      if (allocated(xdesc_emd)) deallocate(xdesc_emd) ;  allocate(xdesc_emd (size(w_params,1)))

      dim_full_partial = 0
      xdesc_emd(1) = config_real(iconf)%nat
      xdesc_emd(2:dim_xdesc_linear) = config_desc(iconf)%pack_energy_linear(1:dim_xdesc_linear-1)
      dim_full_partial = dim_xdesc_linear
      if (activate_k2b)  then
         dim_full_partial = dim_full_partial + dim_kernel_2b
         xdesc_emd(dim_xdesc_linear + 1 : dim_xdesc_linear + dim_kernel_2b) = config_desc(iconf)%pack_energy_k2b(1:dim_kernel_2b)
      end if


      select case (mld_order)

       case (mld_linear)
         continue

       case (mld_quadratic)

         xdesc_emd(dim_full_partial + 1 : dim_xdesc_full) = config_desc(iconf)%pack_energy_quadratic(dim_xdesc_linear:dim_xdesc_quadratic - 1)
         dim_full_partial = dim_full_partial + dim_xdesc_quadratic - dim_xdesc_linear
       case (mld_polyc)
         !call polyc_fill_energy_from_xdesc_train(iconf)
         xdesc_emd(dim_full_partial + 1 :  dim_xdesc_full) = config_desc(iconf)%pack_energy_polyc(dim_xdesc_linear :dim_xdesc_polyc - 1)
         dim_full_partial = dim_full_partial + dim_xdesc_polyc - dim_xdesc_linear

       case (mld_kernel)
         !call polyc_fill_energy_from_xdesc_train(iconf)
         xdesc_emd(dim_full_partial + 1 : dim_xdesc_full) = config_desc(iconf)%pack_energy_kernel(dim_xdesc_linear : dim_xdesc_kernel - 1)
         dim_full_partial = dim_full_partial + dim_xdesc_kernel - dim_xdesc_linear

       case default
         call log_critical('ML: this potential - not implemented in the energy test part, mld order is  '// vtoa(mld_order)//' in '//NAMECURRENT )
         call mld_mpi_abort('milady will stop in '//NAMECURRENT //' with FATAL from energy mld_order')
      end select

      if (dim_full_partial /= dim_xdesc_full) then
         call log_critical('the local dimension of descritors '//   vtoa(dim_full_partial)  //' is different from the expected one ' &
            //vtoa(dim_xdesc_full)//' in '//NAMECURRENT)
         call mld_mpi_abort('milady will stop in energy test '//NAMECURRENT)
      end if

      ene_snap = dot_product(w_params(:, 1), xdesc_emd(:))

      local_ref_energy =  tmp_ene
      local_ref_zbl = 0.d0
      if (zbl_potential .and. zbl_type == zbl_mode_alone) then
         local_ref_zbl    = config_real(iconf)%ezbl
      end if
      ene_snap = ene_snap + local_ref_energy + local_ref_zbl

      _MLD_END_
   end subroutine md_mld_compute_energy



   subroutine md_mld_compute_force(iconf, pack_opt)

      use derived_types, only: config_real, config_desc
      use ml_in_ndm_module, only: mld_order, &
         mld_quadratic, mld_linear, mld_polyc, &
         mld_kernel
      use module_mld_quadratic, only: dim_xdesc_quadratic
      use module_mld_polyc, only: dim_xdesc_polyc
      use module_mld_kernel, only: dim_xdesc_kernel
      use module_kernel_2b, only: activate_k2b, dim_kernel_2b
      use module_zbl, only: zbl_potential, zbl_type, zbl_mode_alone
      use snap, only: dim_xdesc_full, xdesc_fmd 
      use mld_force_mod, only: pack_force_descriptor
      use mld_logger, only: log_critical, vtoa
      use mld_mpi, only: mld_mpi_abort

      integer, intent(in)  :: iconf
      integer  :: ik, ix
      logical, intent(in), optional    :: pack_opt
      logical  :: lpack_tmp
      ! real(kind_double), dimension(:,:,:), allocatable :: xdesc_fmd
      integer :: dim_full_partial
      _NAMECURRENT_("md_mld_compute_force")



      _MLD_BEGIN_
      !if (.not.desc_forces) return
      lpack_tmp = .true.
      if (present(pack_opt)) lpack_tmp = pack_opt

      !write (*,*) 'sizes in md_mld_compute_force', size(xdesc_fmd,1), size(xdesc_fmd,2), size(xdesc_fmd,3)
      if (lpack_tmp) call pack_force_descriptor(iconf)
      if (.not. (config_real(iconf)%has_energy)) call pack_energy_descriptor(iconf)

      if (allocated(xdesc_fmd)) deallocate(xdesc_fmd) ;   allocate (xdesc_fmd(size(w_params,1),config_real(iconf)%nat,3))
      xdesc_fmd(:, :, :) = 0.d0

      dim_full_partial = 0

      ! linear part ....
      do ik = 1, config_real(iconf)%nat
         xdesc_fmd(1, ik, 1:3) = 0.d0
         xdesc_fmd(2:dim_xdesc_linear, ik, 1:3) = config_desc(iconf)%pack_force_linear(1:dim_xdesc_linear-1, 1:3, ik)
         dim_full_partial = dim_xdesc_linear

         ! patch part ...
         if (activate_k2b)  then
            dim_full_partial = dim_full_partial + dim_kernel_2b
            xdesc_fmd(dim_xdesc_linear + 1 : dim_xdesc_linear + dim_kernel_2b, ik, 1:3) = config_desc(iconf)%pack_force_k2b(1:dim_kernel_2b, 1:3, ik)
         end if
      end do

      ! model part ...
      select case (mld_order)

       case (mld_linear)
         continue

       case (mld_quadratic)
         do ik = 1, config_real(iconf)%nat
            xdesc_fmd(dim_full_partial + 1 : dim_xdesc_full, ik, 1:3) = config_desc(iconf)%pack_force_quadratic(dim_xdesc_linear:dim_xdesc_quadratic - 1, 1:3 , ik)
         end do
         dim_full_partial = dim_full_partial + dim_xdesc_quadratic - dim_xdesc_linear

       case (mld_polyc)
         do ik = 1, config_real(iconf)%nat
            xdesc_fmd(dim_full_partial + 1 : dim_xdesc_full, ik, 1:3) = config_desc(iconf)%pack_force_polyc(dim_xdesc_linear:dim_xdesc_polyc - 1, 1:3 , ik)
         end do
         dim_full_partial = dim_full_partial + dim_xdesc_polyc - dim_xdesc_linear

       case (mld_kernel)
         do ik = 1, config_real(iconf)%nat
            xdesc_fmd(dim_full_partial + 1 : dim_xdesc_full, ik, 1:3) = config_desc(iconf)%pack_force_kernel(dim_xdesc_linear:dim_xdesc_kernel - 1, 1:3 , ik)
         end do
         dim_full_partial = dim_full_partial + dim_xdesc_kernel - dim_xdesc_linear

       case default
         call log_critical('ML: this potential - not implemented in the force test part, mld order is  '// vtoa(mld_order)//' in '//NAMECURRENT )
         call mld_mpi_abort('milady will stop in '//NAMECURRENT //' with FATAL from force mld_order')
      end select

      if (dim_full_partial /= dim_xdesc_full) then
         call log_critical('the local dimension of descriptors '//   vtoa(dim_full_partial)  //' is different from the expected one ' &
            //vtoa(dim_xdesc_full)//' in '//NAMECURRENT)
         call mld_mpi_abort('milady will stop in test force'//NAMECURRENT)
      end if

      do ik = 1, config_real(iconf)%nat                ! force index with respect that derivative ik
         do ix = 1, 3
            fp_snap(ix, ik) = dot_product(w_params(:, 1), xdesc_fmd(:, ik, ix))
         end do
      end do
      if (zbl_potential .and. zbl_type == zbl_mode_alone) then
         do ik = 1, config_real(iconf)%nat
            fp_snap(1:3, ik) = fp_snap(1:3, ik) + config_real(iconf)%fzbl(1:3, ik)
         end do
      end if
      _MLD_END_
   end subroutine md_mld_compute_force




   subroutine md_mld_compute_stress(iconf, pack_opt)
      ! This subroutine has strong interaction with neighbours subroutines
      ! The one proposed by Milady and the one proposed by NDM depending if the cell box is small or large.
      ! small=.true. ->  MiLaDy
      ! small=.false. -> NDM
      ! NDM_interaction
      ! MiLaDy_interaction
      use ml_in_ndm_module, only: mld_order, mld_linear, &
         mld_quadratic, mld_polyc, mld_kernel
      use derived_types, only: config_desc, config_real

      use module_mld_quadratic, only: dim_xdesc_quadratic
      use module_mld_polyc, only: dim_xdesc_polyc
      use module_mld_kernel, only: dim_xdesc_kernel
      use module_kernel_2b, only: activate_k2b, dim_kernel_2b
      use module_zbl, only: zbl_potential, zbl_type, zbl_mode_alone
      use snap, only: dim_xdesc_full, xdesc_smd
      use mld_logger, only: log_critical, vtoa
      use mld_mpi, only: mld_mpi_abort

      integer, intent(in)  :: iconf
      integer  :: ix
      logical, intent(in), optional    :: pack_opt
      logical  :: lpack_tmp
      ! real(kind_double), dimension(:,:), allocatable :: xdesc_smd
      integer :: dim_full_partial

      _NAMECURRENT_("md_mld_compute_stress")



      _MLD_BEGIN_
      ! if (.not.desc_forces) return
      lpack_tmp = .true.
      if (present(pack_opt)) lpack_tmp = pack_opt

      if (lpack_tmp) call pack_stress_descriptor(iconf)
      ! write (*,*) "md_mld_compute_stress", iconf, config_real(iconf)%volume

      if (allocated(xdesc_smd))  deallocate(xdesc_smd) ; allocate(xdesc_smd(size(w_params,1), 6))

      xdesc_smd(1, 1:6) = 0.d0
      ! linear part ....
      xdesc_smd(2:dim_xdesc_linear, 1:6) = config_desc(iconf)%pack_stress_linear (1:dim_xdesc_linear-1, 1:6)
      dim_full_partial = dim_xdesc_linear

      ! patch part ...
      if (activate_k2b)  then
         dim_full_partial = dim_full_partial + dim_kernel_2b
         xdesc_smd(dim_xdesc_linear + 1 : dim_xdesc_linear + dim_kernel_2b, 1:6) = config_desc(iconf)%pack_stress_k2b(1:dim_kernel_2b, 1:6)
      end if

      select case (mld_order)

       case (mld_linear)
         continue
       case (mld_quadratic)
         xdesc_smd(dim_full_partial + 1 : dim_xdesc_full, 1:6) = config_desc(iconf)%pack_stress_quadratic(dim_xdesc_linear :dim_xdesc_quadratic - 1, 1:6)
         dim_full_partial = dim_full_partial + dim_xdesc_quadratic - dim_xdesc_linear

       case (mld_polyc)
         xdesc_smd(dim_full_partial + 1 : dim_xdesc_full, 1:6) = config_desc(iconf)%pack_stress_polyc(dim_xdesc_linear :dim_xdesc_polyc - 1, 1:6)
         dim_full_partial = dim_full_partial + dim_xdesc_polyc - dim_xdesc_linear

       case (mld_kernel)
         xdesc_smd(dim_full_partial + 1 : dim_xdesc_full, 1:6) = config_desc(iconf)%pack_stress_kernel(dim_xdesc_linear :dim_xdesc_kernel - 1, 1:6)
         dim_full_partial = dim_full_partial + dim_xdesc_kernel - dim_xdesc_linear


       case default
         call log_critical('ML: this potential - not implemented in the stress test part, mld order is  '// vtoa(mld_order)//' in '//NAMECURRENT )
         call mld_mpi_abort('milady will stop in '//NAMECURRENT //' with FATAL from test mld_order')
      end select


      if (dim_full_partial /= dim_xdesc_full) then
         call log_critical('the local dimension of descriptors '//   vtoa(dim_full_partial)  //' is different from the expected one ' &
            //vtoa(dim_xdesc_full)//' in '//NAMECURRENT)
         call mld_mpi_abort('milady will stop in test force'//NAMECURRENT)
      end if


      do ix = 1, 6
         stress_snap(ix) = dot_product(w_params(:, 1), xdesc_smd(:, ix))
      end do
      if (zbl_potential .and. zbl_type == zbl_mode_alone) then
         stress_snap(1:6) = stress_snap(1:6) + config_real(iconf)%szbl(1:6)
      end if

      _MLD_END_
   end subroutine md_mld_compute_stress

   !/------------------------------------------------------------\
   !                                                             !
   !                 General/restart/input  PART                 !
   !                                                             !
   !\------------------------------------------------------------/


   subroutine read_parameters_for_md
      ! read the parameters file of the snap potential

      use snap, only: dim_xdesc_full, w_params, dim_design_line
      use ml_in_ndm_module, only: rangml,  mld_order, mld_linear, &
         ! char_desc, 
         mld_quadratic, mld_polyc, mld_kernel
      use mld_mpi, only: mld_mpi_abort   
      use module_write_parameters, only: build_name_of_parameters_file   

      implicit none

      character(len=80)    :: file_params
      character(len=120)  :: mld_name, lammps_name 
      integer  :: dim_pot, i
      logical  :: ok
      integer :: unitnumber
      _NAMECURRENT_("read_parameters_for_md")
      _MLD_BEGIN_

      if (allocated(w_params)) deallocate (w_params); allocate (w_params(dim_design_line, 1))
      if (dim_design_line /= dim_xdesc_full) then

         call log_critical("ML: MD_mode reading parameters critical inconsistency dim_xdesc_full not dim_design_line ")
         call log_critical("ML: MD_mode dim_xdesc_full"//vtoa(dim_xdesc_full)//" dim_design_line "//vtoa(dim_design_line))
         call mld_mpi_abort("stop with dim_xdesc_full in"//NAMECURRENT)
      end if 

      !file_params = trim(adjustl(char_desc//'_snap1_params.pot'))
      call build_name_of_parameters_file(mld_name, lammps_name)
      inquire (file=trim(mld_name), exist=ok)
      !inquire (file=file_params, exist=ok)
      if (ok) then
         open (file=trim(mld_name), newunit=unitnumber, action='read')
         read (unitnumber, *) dim_pot

         if (mld_order == mld_linear) then
            if (dim_pot /= (dim_xdesc_full)) then
               if (rangml == 0) write (6, *) 'ML: Fatal the parameters file and the descriptors no have the same dimension'
               stop "read_parameters_for_md linear dimension wrong"
            end if
         end if

         if (mld_order == mld_quadratic) then
            if (dim_pot /= (dim_xdesc_full)) then
               if (rangml == 0) write (6, *) 'ML: Fatal the parameters file and the descriptors no have the same dimension'
               stop "read_parameters_for_md quadratic dimension wrong"
            end if
         end if

         if (mld_order == mld_polyc) then
            if (dim_pot /= (dim_xdesc_full)) then
               if (rangml == 0) write (6, *) 'ML: Fatal the parameters file and the descriptors no have the same dimension'
               stop "read_parameters_for_md polyc dimension wrong"
            end if
         end if

         if (mld_order == mld_kernel) then
            if (dim_pot /= (dim_xdesc_full)) then
               if (rangml == 0) write (6, *) 'ML: Fatal the parameters file and the descriptors no have the same dimension'
               stop "read_parameters_for_md polyc dimension wrong"
            end if
         end if

         !TODOmd ... read only one proc ... then brodcast 
         call log_warning("TODOmd reading w_params in parallel ... problem")  
         call log_warning("TODOmd ... reading on mld_rank 0 then broadcast")
         do i = 1, dim_pot
            read (unitnumber, *) w_params(i, 1)
         end do
         close (unitnumber)
      else
         if (rangml == 0) write (6, *) 'ML: Fatal Error: the potential file is not there. The expected name is: ', file_params
         stop "read_parameters_snap missing the parameters file"
      end if
      _MLD_END_
   end subroutine read_parameters_for_md 

end module main_mld_mod


module module_md_mld
   use mld_logger
   implicit none
contains


   subroutine md_allocate_mld_desc(fp_snap)
      ! MD drivers for Milady
      ! used olny for md
      use module_kind_variables, only: kind_double
      use snap, only: dim_xdesc_linear, dim_xdesc_full
      use ml_in_ndm_module, only: rangml,  mld_order, mld_linear, mld_quadratic, mld_polyc,  &
         polyc_n_poly, polyc_n_hermite, mld_kernel, mld_linear_extended, &
         mld_type_quadratic, mld_type_quadratic_zaxa, mld_type_quadratic_ZX
      use temporary_data_cov, only:  dim_xdesc, dim_xdesc_patch
      use module_mld_quadratic, only: dim_xdesc_quadratic, dim_quadratic
      use module_mld_polyc, only: dim_xdesc_polyc, dim_polyc
      use module_kernel, only: dim_kernel
      use module_mld_kernel, only: dim_xdesc_kernel
      use module_nlinear, only: order_nlinear
      use module_snap_nlinear, only: dim_xdesc_nlinear, dim_nlinear
      use module_kernel_2b, only: dim_kernel_2b, activate_k2b
#ifdef MLD_NDM      
      use gen_com_m_ml, only: im, imm
#else
      use ondm_gen_com_m, only: im, imm
#endif
      implicit none
      real(kind_double), dimension(:,:), allocatable, intent(inout) :: fp_snap

      _NAMECURRENT_("md_allocate_mld_desc")

      _MLD_BEGIN_
      if (dim_xdesc == 0) then
         !if (rangml == 0) write (6, '("ML: dim_xdesc Fatal Error")')
         call log_critical("ML: dim_xdesc Fatal Error in "//NAMECURRENT)
         stop "dim_xdesc is zero in md_allocate_snap"
      end if

      if (im == 0) then
         if (rangml == 0) write (6, '("ML: im Fatal Error")')
         stop "im is zero in md_allocate_snap"
      end if

      if (im /= imm) then
         if (rangml == 0) write (6, *) 'im should be equal to imm', im, imm
      end if
      if (allocated(fp_snap)) deallocate (fp_snap); allocate (fp_snap(3, imm))

      dim_xdesc_linear = 1 + dim_xdesc
      dim_xdesc_patch = 0
      if (activate_k2b) dim_xdesc_patch = dim_xdesc_patch + dim_kernel_2b

      select case (mld_order)
       case (mld_linear, mld_linear_extended)
         dim_xdesc_nlinear = dim_xdesc_linear + order_nlinear*dim_xdesc
         dim_nlinear = order_nlinear*dim_xdesc
         dim_xdesc_full = dim_xdesc_linear + dim_xdesc_patch

       case (mld_quadratic)
         dim_xdesc_quadratic = dim_xdesc_linear + dim_xdesc**2
         dim_quadratic = dim_xdesc**2
         if ((mld_type_quadratic == mld_type_quadratic_zaxa).or.(mld_type_quadratic == mld_type_quadratic_ZX))  then
            dim_xdesc_quadratic = dim_xdesc_linear + dim_xdesc*dim_kernel_2b
            dim_quadratic = dim_xdesc*dim_kernel_2b
         end if        
         dim_xdesc_full = dim_xdesc_quadratic + dim_xdesc_patch

       case (mld_polyc)
         if (polyc_n_poly == 1) dim_xdesc_polyc = 1 + polyc_n_hermite*dim_xdesc
         if (polyc_n_poly == 2) dim_xdesc_polyc = 1 + polyc_n_hermite*dim_xdesc + polyc_n_hermite*dim_xdesc**2
         if (polyc_n_poly == 3) dim_xdesc_polyc = 1 + polyc_n_hermite*dim_xdesc + polyc_n_hermite*dim_xdesc**2 + polyc_n_hermite*dim_xdesc**3
         dim_polyc = dim_xdesc_polyc - dim_xdesc_linear
         dim_xdesc_full = dim_xdesc_polyc + dim_xdesc_patch

       case (mld_kernel)
         dim_xdesc_kernel = dim_xdesc_linear + dim_kernel
         dim_xdesc_full = dim_xdesc_kernel + dim_xdesc_patch

       case default
         call log_critical('ML: this potential - not implemented in the test part, mld order is  '// vtoa(mld_order)//' in '//NAMECURRENT )
         stop 'FATAL in md_allocate_snap with mld_order'
      end select


      _MLD_END_
   end subroutine md_allocate_mld_desc
end module module_md_mld

module module_test_mld

#if(PARA)
   !use mpi
   use mld_mpi
#endif

   use mld_logger
   implicit none


   !/------------------------------------------------------------\
   !                                                             !
   !                           TEST PART                         !
   !                                                             !
   !\------------------------------------------------------------/

contains

   subroutine compute_variance_matrix()
      use module_kind_variables, only: kind_double
      use ml_in_ndm_module, only: rangml
      use module_db_setup, only:  iconf_data
      use derived_types, only: config_desc, config_real
      use module_kernel, only: dim_kernel, info_kernel
      use module_variance, only: sigmam_f,  VInv

      implicit none
      external :: train_deallocate_desc
      real(kind_double), dimension(:, :), allocatable      :: KKmat
      real(kind_double), dimension(:, :), allocatable         :: MKmat
      real(kind_double), dimension(:, :), allocatable      :: tempmat
      real(kind_double), dimension(:, :), allocatable      :: tempmat2
      integer(4)  :: ii, ic,  ia, M
      real(kind_double)  :: ss
      character ::  fname*80
      integer :: lwork_svd, info_svd
      real(kind_double), dimension(:,:), allocatable :: u_svd, vt_svd, u
      real(kind_double), dimension(:), allocatable :: sigma_svd, work_svd
      logical :: post_desc 

      _NAMECURRENT_("compute_variance_matrix")
      _MLD_BEGIN_
      M=0
      do ic = 1, iconf_data
         if (.not. (config_real(ic)%train)) cycle
         M = M+config_real(ic)%nat
      end do
      !TORC! call MPI_ALLREDUCE(MPI_IN_PLACE, sigmam_f, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, codeml)
      !TODO:  TOTEST
      call comm_mld%sum(sigmam_f)
      !write (6, *) 'sigmam_f', sigmam_f


      ii=0
      if (allocated(MKmat)) deallocate(MKmat); allocate (MKmat(M,dim_kernel))
      do ic = 1, iconf_data
         if (.not. (config_real(ic)%train)) cycle
         call test_if_config_is_small(ic)
         call calc_neighbours(ic)
         post_desc=.true. 
         call compute_descriptors(ic, post_desc)
         do ia = 1, config_real(ic)%nat
            ii=ii+1
            ! MKmat(ii,:) = config_desc(ic)%energy_kernel(1:dim_kernel, ia)
            MKmat(ii,:) = config_desc(ic)%energy_kernel(:, ia)
         end do
         call train_deallocate_desc(ic)
      end do


      !maxv=0.d0
      !minv=1.d10
      !do ic = 1, M
      !    do ia = 1, dim_kernel
      !        maxv=max(MKmat(ic,ia),maxv)
      !        minv=min(MKmat(ic,ia),minv)
      !    end do
      !end do
      !write (6, *) '!!! MKmat !!!', maxv, minv

      if (allocated(tempmat)) deallocate(tempmat); allocate (tempmat(dim_kernel, dim_kernel))
      tempmat(:,:)=0.d0

      call dgemm('T', 'N', dim_kernel, dim_kernel, M, 1.d0, MKmat, M, MKmat, M, 0.d0, tempmat, dim_kernel)

      !  maxv=0.d0
      !  minv=1.d10
      !  do ic = 1, dim_kernel
      !      do ia = 1, dim_kernel
      !          maxv=max(tempmat(ic,ia),maxv)
      !          minv=min(tempmat(ic,ia),minv)
      !      end do
      !  end do
      !  write (6, *) '!!! MK^TxMK !!!', maxv, minv
      deallocate(MKmat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !write (6, *) '!!! 1 !!!', tempmat(1,1), tempmat(dim_kernel,dim_kernel)

      if (allocated(KKmat)) deallocate (KKmat); allocate (KKmat(dim_kernel, dim_kernel))
      KKmat(:,:)=0.d0
      do ii = 1, dim_kernel
         ic = info_kernel(ii)%iconf
         ia = info_kernel(ii)%ia
         call test_if_config_is_small(ic)
         call calc_neighbours(ic)
         post_desc = .true. 
         call compute_descriptors(ic, post_desc)

         ! KKmat(ii,:) = config_desc(ic)%energy_kernel(1:dim_kernel, ia)
         KKmat(ii,:) = config_desc(ic)%energy_kernel(:, ia)
         call train_deallocate_desc(ic)
      end do

      !  maxv=0.d0
      !  minv=1.d10
      !  do ic = 1, dim_kernel
      !      do ia = 1, dim_kernel
      !          maxv=max(KKmat(ic,ia),maxv)
      !          minv=min(KKmat(ic,ia),minv)
      !      end do
      !  end do
      !  write (6, *) '!!! KKmat !!!', maxv, minv


      if (allocated(tempmat2)) deallocate(tempmat2); allocate (tempmat2(dim_kernel, dim_kernel))
      tempmat2(:,:)=0.d0
      tempmat2(:,:)=dble(tempmat(:,:)/dble(sigmam_f)+KKmat(:,:))
      !write (6, *) '!!! tempmat, sigmam_f, KKmat, tempmat2 !!!',tempmat(1,5),sigmam_f,KKmat(1,5),tempmat2(1,5)
      !write (6, *) '!!! tempmat, sigmam_f, KKmat, tempmat2 !!!',tempmat(100,200),sigmam_f,KKmat(100,200),tempmat2(100,200)


      if (rangml == 0) then
         fname='V.dat'
         open(unit=3060,file=fname,status='unknown')
      end if

      !!!!!!!!!!!!!!!!!!!!!! symmetrization !!!!!!!!!!!!!!!!!!!!!
      ! minv=1.d10
      ! maxv=0.d0
      ! do ii = 1, dim_kernel
      !    do ic = ii, dim_kernel
      !        if (tempmat2(ii, ic)/=tempmat2(ic, ii)) then
      !            tempmat2(ii, ic) = 0.5d0*(tempmat2(ii, ic)+tempmat2(ic, ii))
      !            tempmat2(ic, ii) = tempmat2(ii, ic)
      !        end if
      !        if (rangml == 0) then
      !            write(3060,'(E25.12)') tempmat2(ii, ic)
      !        end if
      !        maxv=max(maxv, tempmat2(ii, ic))
      !        minv=min(minv, tempmat2(ii, ic))
      !    end do
      ! end do

      ! write (6, *) '!!! After symmetrization !!!',tempmat2(1,5),tempmat2(100,200), minv, maxv
      if (rangml == 0) close(3060)

      ! ia=0
      ! do ii = 1, dim_kernel
      !    do ic = 1, dim_kernel
      !       if (tempmat2(ii,ic)/=tempmat2(ii,ic)) then
      !           write (6, *) 'NaN in tempmat2', ii, ic
      !       end if
      !       if (tempmat2(ii,ic)/=0.d0) then
      !           ia=1
      !       end if
      !    end do
      ! end do
      ! if (ia==0) then
      !      write (6, *) 'All 0 in tempmat2'
      ! end if

      if (allocated(tempmat)) deallocate(tempmat)
      if (allocated(KKmat)) deallocate(KKmat)

      if (allocated(VInv)) deallocate (VInv); allocate (VInv(dim_kernel, dim_kernel))
      VInv(:,:)=tempmat2(:,:)

      !if (allocated(tempmat)) deallocate(tempmat);allocate(tempmat(dim_kernel,dim_kernel))
      !tempmat(:,:)=tempmat2(:,:)

      allocate(u_svd(dim_kernel,dim_kernel))
      allocate(vt_svd(dim_kernel,dim_kernel))
      allocate(sigma_svd(dim_kernel))
      allocate(work_svd(5*dim_kernel))
      allocate(u(dim_kernel,dim_kernel))
      u_svd(:,:)=0.d0
      vt_svd(:,:)=0.d0
      sigma_svd(:)=0.d0
      work_svd(:)=0.d0
      u(:,:)=0.d0

      call  dgesvd('A','A', dim_kernel, dim_kernel, VInv, dim_kernel, &
         sigma_svd, u_svd, dim_kernel, vt_svd, dim_kernel, work_svd, -1, info_svd)
      if (info_svd .ne. 0) write(6,'("ML: error for dgesvd1 info: ", i8)') info_svd

      lwork_svd = int(work_svd(1)) + 10

      if (allocated(work_svd)) deallocate(work_svd) ; allocate(work_svd(lwork_svd))
      call  dgesvd('A','A', dim_kernel, dim_kernel, tempmat2, dim_kernel, &
         sigma_svd, u_svd, dim_kernel, vt_svd, dim_kernel, work_svd,lwork_svd, info_svd)
      if (info_svd .ne. 0) write(6,'("ML: error for dgesvd2 info: ", i8)') info_svd

      write (6, *) '!!! sigma_svd !!!', sigma_svd(1), sigma_svd(2), sigma_svd(3)
      !if (rangml == 0) then
      !  do ii=1,dim_kernel
      !    write (6, *) '!!! sigma_svd !!!', sigma_svd(ii)
      !  end do
      !end if

      do ii=1,dim_kernel
         !if (sigma_svd(ii)>1.d-14) then
         if (sigma_svd(ii)>1.d-9) then
            ss=1.d0/sigma_svd(ii)
         else
            !if ((sigma_svd(ii)<1.d-10) .and. (rangml == 0)) then
            !write (6, *) '!!! Pseudo-inverse error !!! ', ii, ss
            ss=0.d0!sigma_svd(ii)
         end if
         !call dscal(dim_kernel, ss, u_svd(:,ii),1)
         u(:,ii)=ss*u_svd(:,ii)
      end do

      VInv(:,:)=0.d0
      call dgemm( "T", "T", dim_kernel, dim_kernel, dim_kernel, 1.d0, vt_svd, dim_kernel, u, dim_kernel, 0.d0, VInv, dim_kernel)

      if (allocated(u_svd)) deallocate (u_svd)
      if (allocated(vt_svd)) deallocate (vt_svd)
      if (allocated(sigma_svd)) deallocate (sigma_svd)
      if (allocated(work_svd)) deallocate (work_svd)
      if (allocated(u)) deallocate (u)


      _MLD_END_
   end subroutine compute_variance_matrix


   subroutine md_force_variance(iconf, pack_opt)
      use module_kernel, only: dim_kernel
      use module_variance, only: VInv
      use module_kind_variables, only: kind_double
      use temporary_data_cov, only: dim_xdesc
      use derived_types, only: config_real, config_desc
      use mld_force_mod, only: pack_force_descriptor
      use snap, only: fp_var
      ! use module_mld_quadratic, only: dim_xdesc_quadratic
      ! use module_mld_polyc, only: dim_xdesc_polyc
      ! use module_mld_kernel, only: dim_xdesc_kernel
      implicit none
      integer, intent(in)  :: iconf
      integer  :: ik, ix, ii
      logical, intent(in), optional    :: pack_opt
      logical  :: lpack_tmp
      real(kind_double), dimension(:,:,:), allocatable :: xdesc_fmd
      real(kind_double), dimension(:), allocatable :: tempmat

      _NAMECURRENT_("md_force_variance")
      _MLD_BEGIN_

      lpack_tmp = .true.
      if (present(pack_opt)) lpack_tmp = pack_opt

      if (lpack_tmp) call pack_force_descriptor(iconf)
      if (allocated(xdesc_fmd)) deallocate(xdesc_fmd); allocate(xdesc_fmd(dim_kernel,config_real(iconf)%nat,3))
      if (allocated(tempmat)) deallocate(tempmat); allocate (tempmat(dim_kernel))
      tempmat(:) = 0.d0
      xdesc_fmd(:, :, :) = 0.d0


      do ik = 1, config_real(iconf)%nat
         ! xdesc_fmd(1, ik, 1:3) = 0.d0
         xdesc_fmd(1:dim_kernel, ik, 1:3) = config_desc(iconf)%pack_force_kernel(dim_xdesc + 1:dim_xdesc + dim_kernel, 1:3, ik)
      end do

      do ik = 1, config_real(iconf)%nat
         do ix = 1, 3
            tempmat(:)=0.d0
            do ii = 1, dim_kernel
               tempmat(ii) = dot_product(VInv(:,ii), xdesc_fmd(:, ik, ix))
            end do
            fp_var(ix, ik) = dot_product(tempmat(:), xdesc_fmd(:, ik, ix))
         end do
      end do

      _MLD_END_
   end subroutine md_force_variance

   subroutine main_test_gather_mld(i_e_test_snap,i_f_test_snap,i_s_test_snap)

      !use mpi
      use snap, only: y_e_test_snap, y_e_test_base, y_e_p_a_test_snap, y_e_p_a_test_base, &
         y_f_test_snap, var_f_test_snap, y_f_test_base, &
         y_s_test_snap, y_s_test_base, &
         dim_ene_test_snap, dim_force_test_snap, dim_stress_test_snap, &
         e_test_snap, s_test_snap, f_test_snap
      use mld_subworld
      use my_mpi_subroutines, only: subworlds_allreduce_vect_double, subworlds_allreduce_vect_int

      implicit none
      integer, intent(inout)  :: i_e_test_snap, i_f_test_snap, i_s_test_snap
      integer, dimension(:,:), allocatable :: i_test_on_proc, i_test_offsets
      real(kind=kind(0.d0)), dimension(:), allocatable :: tmp
      integer, dimension(:), allocatable :: tmp_int, i_tmp 
      integer :: i, group, start, end

      !TORC! call my_barrier_mld(codeml)
      call comm_mld%barrier

      ! assembling all local vectors into global vectors:
      ! i_*_test_on_proc contains the number of elements filled in each local vector
      ! i_*_test_offsets contains the correct offset for each proc in the global object where to put the local object
      allocate (i_test_on_proc(3,0:nb_subworlds - 1))
      allocate (i_test_offsets(3,0:nb_subworlds - 1))
      !
      allocate (i_tmp(0:nb_subworlds - 1))
      !
      i_test_on_proc(:,:) = 0
      i_test_on_proc(1,id_subworld) = i_e_test_snap
      i_test_on_proc(2,id_subworld) = i_f_test_snap
      i_test_on_proc(3,id_subworld) = i_s_test_snap
      i_tmp(:)=0
      i_tmp(id_subworld) =  i_e_test_snap
      !call subworlds_allreduce_vect_int(i_test_on_proc(1,:))
      call subworlds_allreduce_vect_int(i_tmp)
      i_test_on_proc(1,:)=i_tmp(:)

      i_tmp(:)=0
      i_tmp(id_subworld) =  i_f_test_snap
      !call subworlds_allreduce_vect_int(i_test_snap_on_proc(2,:))
      call subworlds_allreduce_vect_int(i_tmp)
      i_test_on_proc(2,:)=i_tmp(:)

      i_tmp(:)=0
      i_tmp(id_subworld) =  i_s_test_snap
      !call subworlds_allreduce_vect_int(i_test_on_proc(3,:))
      call subworlds_allreduce_vect_int(i_tmp) 
      i_test_on_proc(3,:)=i_tmp(:)

      i_test_offsets(:,0) = 0
      do group=1, nb_subworlds-1
         do i=1,3
            i_test_offsets(i,group) = SUM(i_test_on_proc(i,0:group-1))
         end do
      end do

      ! translating the local portion into the correct slot in the global object and gathering it
      allocate(tmp(i_e_test_snap)); allocate(tmp_int(i_e_test_snap))
      start = i_test_offsets(1,id_subworld)+1
      end   = i_test_offsets(1,id_subworld)+i_e_test_snap
      if (dim_ene_test_snap > 0 ) then
        tmp(:) = y_e_test_snap(1:i_e_test_snap); y_e_test_snap(1:i_e_test_snap) = 0.d0; y_e_test_snap(start:end) = tmp(:)
        tmp(:) = y_e_test_base(1:i_e_test_snap); y_e_test_base(1:i_e_test_snap) = 0.d0; y_e_test_base(start:end) = tmp(:)
        tmp(:) = y_e_p_a_test_snap(1:i_e_test_snap); y_e_p_a_test_snap(1:i_e_test_snap) = 0.d0; y_e_p_a_test_snap(start:end) = tmp(:)
        tmp(:) = y_e_p_a_test_base(1:i_e_test_snap); y_e_p_a_test_base(1:i_e_test_snap) = 0.d0; y_e_p_a_test_base(start:end) = tmp(:)
        tmp_int(:) = e_test_snap(1:i_e_test_snap); e_test_snap(1:i_e_test_snap) = 0; e_test_snap(start:end) = tmp_int(:)
        call subworlds_allreduce_vect_double(y_e_test_snap)
        call subworlds_allreduce_vect_double(y_e_test_base)
        call subworlds_allreduce_vect_double(y_e_p_a_test_snap)
        call subworlds_allreduce_vect_double(y_e_p_a_test_base)
        call subworlds_allreduce_vect_int(e_test_snap)
      end if 
      deallocate(tmp); deallocate(tmp_int)

      allocate(tmp(i_f_test_snap)); allocate(tmp_int(i_f_test_snap))
      start = i_test_offsets(2,id_subworld) + 1
      end   = i_test_offsets(2,id_subworld) + i_f_test_snap
      if (dim_force_test_snap > 0 ) then 
        tmp(:) = y_f_test_snap(1:i_f_test_snap); y_f_test_snap(1:i_f_test_snap) = 0.d0; y_f_test_snap(start:end) = tmp(:)
        tmp(:) = var_f_test_snap(1:i_f_test_snap); var_f_test_snap(1:i_f_test_snap) = 0.d0; var_f_test_snap(start:end) = tmp(:)
        tmp(:) = y_f_test_base(1:i_f_test_snap); y_f_test_base(1:i_f_test_snap) = 0.d0; y_f_test_base(start:end) = tmp(:)
        tmp_int(:) = f_test_snap(1:i_f_test_snap); f_test_snap(1:i_f_test_snap) = 0; f_test_snap(start:end) = tmp_int(:)
        call subworlds_allreduce_vect_double(y_f_test_snap)
        call subworlds_allreduce_vect_double(var_f_test_snap)
        call subworlds_allreduce_vect_double(y_f_test_base)
        call subworlds_allreduce_vect_int(f_test_snap)
      end if 
      deallocate(tmp); deallocate(tmp_int)

      allocate(tmp(i_s_test_snap)); allocate(tmp_int(i_s_test_snap))
      start = i_test_offsets(3,id_subworld) + 1
      end   = i_test_offsets(3,id_subworld) + i_s_test_snap
      if (dim_stress_test_snap > 0 ) then 
        tmp(:) = y_s_test_snap(1:i_s_test_snap); y_s_test_snap(1:i_s_test_snap) = 0.d0; y_s_test_snap(start:end) = tmp(:)
        tmp(:) = y_s_test_base(1:i_s_test_snap); y_s_test_base(1:i_s_test_snap) = 0.d0; y_s_test_base(start:end) = tmp(:)
        tmp_int(:) = s_test_snap(1:i_s_test_snap); s_test_snap(1:i_s_test_snap) = 0; s_test_snap(start:end) = tmp_int(:)
        call subworlds_allreduce_vect_double(y_s_test_snap)
        call subworlds_allreduce_vect_double(y_s_test_base)
        call subworlds_allreduce_vect_int(s_test_snap)
      end if 
      deallocate(tmp); deallocate(tmp_int)

      ! updating i_*_test_snap so that it matches the global dimension of test
      i_e_test_snap = SUM(i_test_on_proc(1,:))
      i_f_test_snap = SUM(i_test_on_proc(2,:))
      i_s_test_snap = SUM(i_test_on_proc(3,:))

      deallocate(i_test_on_proc); deallocate(i_test_offsets)


   end subroutine main_test_gather_mld

   subroutine write_test_output(i_e_test_snap, i_f_test_snap, i_s_test_snap)

      use snap, only: dim_ene_test_snap, dim_force_test_snap, dim_stress_test_snap, &
         y_f_test_base, y_f_test_snap, y_e_test_base, y_e_test_snap, var_f_test_snap,&
         y_s_test_base, y_s_test_snap, &
         e_test_snap, s_test_snap, f_test_snap
      use module_variance, only: force_variance
      use derived_types, only: config_real

      integer, intent(in)  :: i_e_test_snap, i_f_test_snap, i_s_test_snap
      integer  :: i, ix, einp, finp, sinp, idata

      if (dim_ene_test_snap > 0) then
         open (file="test_energy.out", newunit=einp, action='write')
         do i=1,i_e_test_snap
            idata= e_test_snap(i)
            write (einp, '(2e20.10," ", (a), "  ", (a))') y_e_test_snap(i), y_e_test_base(i), config_real(idata)%class, config_real(idata)%filename
         end do
         close (einp)
      end if

      if (dim_force_test_snap > 0) then
         open (file="test_force.out", newunit=finp, action='write')
         i = 0
         do while(i < i_f_test_snap)
            idata = f_test_snap(i+1)
            do ix = 1, 3*config_real(idata)%nat
               if(force_variance) then
                  write (finp, '(3e20.10," ", (a), "  ", (a))')y_f_test_snap(i + ix), y_f_test_base(i + ix), &
                     var_f_test_snap(i + ix),config_real(idata)%class, config_real(idata)%filename
               else
                  write (finp, '(2e20.10," ", (a), "  ", (a))') y_f_test_snap(i + ix), y_f_test_base(i + ix), &
                     config_real(idata)%class, config_real(idata)%filename
               end if
            end do
            if(config_real(idata)%nat == 0) then
               i = i + 1
            else
               i = i + 3*config_real(idata)%nat
            end if
         end do
         close (finp)
      end if
      if (dim_stress_test_snap > 0) then
         open (file="test_stress.out", newunit=sinp, action='write')
         i = 0
         do while(i < i_s_test_snap)
            idata = s_test_snap(i+1)
            do ix = 1, 6
               write (sinp, '(2e20.10," ", (a),"  ", (a))') y_s_test_snap(i + ix), y_s_test_base(i + ix), config_real(idata)%class, config_real(idata)%filename
            end do
            i = i + 6
         end do
         close (sinp)
      end if

   end subroutine


   subroutine main_test_mld()

      use ml_in_ndm_module, only: rangml,   &
         desc_forces, &
         ml_type, ml_type_krr, train_time, write_test_design_matrix
      use module_db_setup, only:  iconf_data   
      use derived_types, only: config_real
      use module_kind_variables, only: kind_double
      use time_check_general, only: test_time_total, &
         test_time_desc, &
         test_time_eval, &
         test_time_neigh, MY_MPI_WTIME, debug_time, tot_time, test_tot_time
      use module_kernel, only: time_for_desc_kernel, test_time_for_desc_kernel
      use module_db_poscar, only: iread_energy, i_start_conf, i_final_conf,  procs_per_file
      use snap, only: w_params, dim_ene_test_snap, dim_force_test_snap, dim_stress_test_snap,  &
         ene_snap, stress_snap, fp_snap, &
         y_f_test_base, y_f_test_snap, y_e_test_base, y_e_test_snap, var_f_test_snap,&
         y_s_test_base, y_s_test_snap, &
         y_e_p_a_test_base, y_e_p_a_test_snap, fp_var, &
         e_test_snap, s_test_snap, f_test_snap, &
         xdesc_emd, xdesc_fmd, xdesc_smd
      use main_mld_mod, only: md_mld_compute_energy, md_mld_compute_force, md_mld_compute_stress, &
         test_error_snap, test_error_by_class
      use module_md_mld, only: md_allocate_mld_desc
      use module_variance, only: force_variance
      use module_kernel_2b, only: test_time_for_desc_kernel_2b, time_for_desc_kernel_2b, activate_k2b
      use set_limits, only: set_limit_for_configs_with_MPI_grid
      use module_write_design_matrix, only: dump_design_matrix_test, init_dump_desgin_matrix_test
      use mld_subworld


      implicit none

      !external :: prepare_test_dimensions,train_deallocate_desc
      integer  :: i, ix, ik, iwmat, dim_wmat, dunit, iunit 
      real(kind_double), dimension(:,:), allocatable :: wmat
      integer, dimension(:), allocatable :: wtags
      real(kind_double), dimension(:), allocatable :: wyy
      integer  :: i_e_test_snap, i_f_test_snap, i_s_test_snap
      real(kind_double)    :: t0, t1, t2, t3, t4, t5, &
         local_ref_energy, local_ref_zbl
      logical :: post_desc

      _NAMECURRENT_("main_test_mld")



      _MLD_BEGIN_
      train_time=.false.

      call init_subworld(procs_per_file)
      call set_limit_for_configs_with_MPI_grid(iconf_data, i_start_conf, i_final_conf)
      call prepare_test_dimensions()

      !testing the training
      if (dim_ene_test_snap > 0) then
         if (allocated(y_e_test_snap)) deallocate (y_e_test_snap); allocate (y_e_test_snap(dim_ene_test_snap)); y_e_test_snap = 0.d0
         if (allocated(y_e_test_base)) deallocate (y_e_test_base); allocate (y_e_test_base(dim_ene_test_snap)); y_e_test_base = 0.d0
         if (allocated(y_e_p_a_test_snap)) deallocate (y_e_p_a_test_snap); allocate (y_e_p_a_test_snap(dim_ene_test_snap)); y_e_p_a_test_snap = 0.d0
         if (allocated(y_e_p_a_test_base)) deallocate (y_e_p_a_test_base); allocate (y_e_p_a_test_base(dim_ene_test_snap)); y_e_p_a_test_base = 0.d0
         if (allocated(e_test_snap)) deallocate (e_test_snap); allocate (e_test_snap(dim_ene_test_snap)); e_test_snap = 0
      end if

      if (dim_force_test_snap > 0) then
         if (allocated(y_f_test_snap)) deallocate (y_f_test_snap); allocate (y_f_test_snap(dim_force_test_snap)); y_f_test_snap = 0.d0
         if (allocated(var_f_test_snap)) deallocate (var_f_test_snap); allocate (var_f_test_snap(dim_force_test_snap)); var_f_test_snap = 0.d0
         if (allocated(y_f_test_base)) deallocate (y_f_test_base); allocate (y_f_test_base(dim_force_test_snap)); y_f_test_base = 0.d0
         if (allocated(f_test_snap)) deallocate (f_test_snap); allocate (f_test_snap(dim_force_test_snap)); f_test_snap = 0
      end if


      if (dim_stress_test_snap > 0) then
         if (allocated(y_s_test_snap)) deallocate (y_s_test_snap); allocate (y_s_test_snap(dim_stress_test_snap)); y_s_test_snap = 0.d0
         if (allocated(y_s_test_base)) deallocate (y_s_test_base); allocate (y_s_test_base(dim_stress_test_snap)); y_s_test_base = 0.d0
         if (allocated(s_test_snap)) deallocate (s_test_snap); allocate (s_test_snap(dim_stress_test_snap)); s_test_snap = 0
      end if

      i_e_test_snap = 0
      i_f_test_snap = 0
      i_s_test_snap = 0
      test_time_total = 0.d0
      test_time_eval = 0.d0
      test_time_desc = 0.d0
      test_time_neigh = 0.d0
      t0 = MY_MPI_WTIME()

      if (force_variance) call compute_variance_matrix
      if (write_test_design_matrix) then 
         call init_dump_desgin_matrix_test (dunit, iunit)
      end if 
      
      do i = i_start_conf, i_final_conf
         if (config_real(i)%train) cycle
         !call read_poscar_sasha(trim(config_real(i)%filename),i)
         if (write_test_design_matrix) then 
           dim_wmat = 0
           if (config_real(i)%has_energy) dim_wmat = dim_wmat + 1
           if ((config_real(i)%has_force) .and. (desc_forces)) dim_wmat = dim_wmat + 3*config_real(i)%nat
           if ((config_real(i)%has_stress) .and. (desc_forces)) dim_wmat = dim_wmat + 6
           if (dim_wmat > 0) then 
             if (allocated(wmat)) deallocate (wmat); allocate (wmat(size(w_params,1),  dim_wmat))
             if (allocated(wtags)) deallocate (wtags); allocate (wtags(dim_wmat))
             if (allocated(wyy)) deallocate (wyy); allocate (wyy(dim_wmat))
             wmat(:,:) = 0.d0 ; wtags(:) = 0 ; wyy(:) = 0.d0
           end if
           iwmat = 0 
         end if 
         t1 = MY_MPI_WTIME()
         call test_if_config_is_small(i)
         call calc_neighbours(i)
         t2 = MY_MPI_WTIME()
         test_time_neigh = test_time_neigh + (t2 - t1)
         !!!write (6,*) 'before md_allocate', i, rangml
         call md_allocate_mld_desc(fp_snap)
         !!!write (6,*) 'after md_allocate', i, rangml
         post_desc=.true.
         call compute_descriptors(i, post_desc)
         t3 = MY_MPI_WTIME()
         test_time_desc = test_time_desc + (t3 - t2)

         !debug if (rangml==0) write (6,*) 'after',
         !mpi_rangml =0
         !mpi_rangml
         !if (rangml==0) then
         if (config_real(i)%has_energy) then
            i_e_test_snap = i_e_test_snap + 1
            e_test_snap(i_e_test_snap) = i
            call md_mld_compute_energy(i, pack_opt=.true.)
            local_ref_energy = config_real(i)%ref_energy
            local_ref_zbl = config_real(i)%ezbl
            y_e_test_snap(i_e_test_snap) = ene_snap
            y_e_test_base(i_e_test_snap) = config_real(i)%energy(iread_energy) + local_ref_energy
            y_e_p_a_test_snap(i_e_test_snap) = ene_snap/dble(config_real(i)%nat)
            y_e_p_a_test_base(i_e_test_snap) = (config_real(i)%energy(iread_energy) + local_ref_energy)/dble(config_real(i)%nat)

            !write (778, '(5e20.10," ", (a), "  ", (a))') ene_snap, y_e_test_base(i_e_test_snap), &
            !                                             ene_snap- y_e_test_base(i_e_test_snap), &
            !                                             local_ref_zbl, &
            !                                             local_ref_energy, &
            !                                             config_real(i)%class, config_real(i)%filename
            if (write_test_design_matrix) then 
               iwmat = iwmat + 1
               wmat(:, iwmat) = xdesc_emd 
               wtags(iwmat) = 1
               wyy(iwmat) = config_real(i)%energy(iread_energy) + local_ref_energy
            end if 
         end if
         if ((config_real(i)%has_force) .and. (desc_forces)) then
            call md_mld_compute_force(i, pack_opt=.true.)
            if (force_variance) then
               if (allocated(fp_var)) deallocate (fp_var); allocate(fp_var(3,config_real(i)%nat))
               call md_force_variance(i, pack_opt=.true.)
               do ix = 1, 3
                  y_f_test_snap(i_f_test_snap + (ix - 1)*config_real(i)%nat + 1:i_f_test_snap + ix*config_real(i)%nat) = fp_snap(ix, 1:config_real(i)%nat)
                  var_f_test_snap(i_f_test_snap + (ix - 1)*config_real(i)%nat + 1:i_f_test_snap + ix*config_real(i)%nat) = fp_var(ix, 1:config_real(i)%nat)
                  y_f_test_base(i_f_test_snap + (ix - 1)*config_real(i)%nat + 1:i_f_test_snap + ix*config_real(i)%nat) = config_real(i)%force(ix, 1:config_real(i)%nat)
               end do
            else
               do ix = 1, 3
                  y_f_test_snap(i_f_test_snap + (ix - 1)*config_real(i)%nat + 1:i_f_test_snap + ix*config_real(i)%nat) = fp_snap(ix, 1:config_real(i)%nat)
                  y_f_test_base(i_f_test_snap + (ix - 1)*config_real(i)%nat + 1:i_f_test_snap + ix*config_real(i)%nat) = config_real(i)%force(ix,1:config_real(i)%nat)
               end do
            end if
            f_test_snap(i_f_test_snap+1:i_f_test_snap + 3*config_real(i)%nat) = i
            i_f_test_snap = i_f_test_snap + 3*config_real(i)%nat

            if (write_test_design_matrix) then 
              do ik = 1, config_real(i)%nat
                do ix = 1, 3
                   iwmat = iwmat + 1
                   wmat(:, iwmat) = xdesc_fmd(:, ik, ix)
                   wtags(iwmat) = 2
                   wyy(iwmat) = config_real(i)%force(ix, ik)
                end do
              end do  
            end if
         end if

         if ((config_real(i)%has_stress) .and. (desc_forces)) then
            call md_mld_compute_stress(i, pack_opt=.true.)
            y_s_test_snap(i_s_test_snap + 1:i_s_test_snap + 6) = stress_snap(1:6)
            y_s_test_base(i_s_test_snap + 1:i_s_test_snap + 6) = config_real(i)%stress(1:6)
            s_test_snap(i_s_test_snap+1:i_s_test_snap + 6) = i
            i_s_test_snap = i_s_test_snap + 6
            if (write_test_design_matrix) then
              do ik = 1, 6
                iwmat = iwmat + 1
                wmat(:, iwmat) = xdesc_smd(:,ik)
                wtags(iwmat) = 3
                wyy(iwmat) = config_real(i)%stress(ik)
              end do  
            end if
         end if
         !end if !mpi_rangml
         if (write_test_design_matrix) call dump_design_matrix_test  (i, dunit, iunit, dim_wmat, wmat, wtags, wyy)
         call train_deallocate_desc(i)
         !mpi_rangml end if ! mpi_rangml ==0
         !end if
         t4 = MY_MPI_WTIME()
         test_time_eval = test_time_eval + (t4 - t3)
         !if (scalapack_driver) call blacs_barrier(context,'A')  ! NO barrier here as every proc is managing his own loop
      end do

      call comm_mld%barrier()

      call main_test_gather_mld(i_e_test_snap,i_f_test_snap,i_s_test_snap)
      if (rangml == 0) then
         call write_test_output(i_e_test_snap, i_f_test_snap, i_s_test_snap)
      end if

      !compute the test mpi_rangml==0
      if (rangml == 0) then
         call test_error_snap
         call test_error_by_class
      end if
      !call test_error_snap
      t5 = MY_MPI_WTIME()
      test_time_total = test_time_total + (t5 - t0)
      if (ml_type == ml_type_krr) then
         test_time_for_desc_kernel = time_for_desc_kernel
      end if
      if (activate_k2b) then
         test_time_for_desc_kernel_2b = time_for_desc_kernel_2b
      end if

      if (debug_time) then
         test_tot_time = tot_time
      end if

      if(allocated(e_test_snap)) deallocate(e_test_snap);
      if(allocated(f_test_snap)) deallocate(f_test_snap);
      if(allocated(s_test_snap)) deallocate(s_test_snap);
      call close_subworld()

      _MLD_END_
   end subroutine main_test_mld

end module module_test_mld

subroutine deallocate_fit
   use snap, only: Amat, Bmat, zmat, &
      fp_snap, fp_var, xdesc_emd, xdesc_fmd, xdesc_smd
   use module_ml_scalapack, only : sca_Amat, desc_sca_Amat, &
      sca_Amat_big, desc_sca_Amat_big, &
      sca_Lmat, desc_sca_Lmat, &
      sca_Qmat, desc_sca_Qmat, &
      sca_Cmat, desc_sca_Cmat,&
      sca_phia, desc_sca_phia, &
      sca_ymat, desc_sca_ymat, &
      sca_ymat_copy, desc_sca_ymat_copy, &
      sca_ymat_qr_svd, desc_sca_ymat_qr_svd, &
      sca_w_params, desc_sca_w_params
      
   implicit none

   if (allocated(Amat))    deallocate(Amat)
   if (allocated(Bmat))    deallocate(Bmat)
   if (allocated(zmat))    deallocate(zmat)
   if (allocated(fp_snap)) deallocate(fp_snap)
   if (allocated(fp_var))  deallocate(fp_var)

   if (allocated(sca_Amat)) deallocate(sca_Amat)
   if (allocated(desc_sca_Amat)) deallocate(desc_sca_Amat)
   if (allocated(sca_Amat_big)) deallocate(sca_Amat_big)
   if (allocated(desc_sca_Amat_big)) deallocate(desc_sca_Amat_big)
   if (allocated(sca_Lmat)) deallocate(sca_Lmat)
   if (allocated(desc_sca_Lmat)) deallocate(desc_sca_Lmat)
   if (allocated(sca_Qmat)) deallocate(sca_Qmat)
   if (allocated(desc_sca_Qmat)) deallocate(desc_sca_Qmat)
   if (allocated(sca_Cmat)) deallocate(sca_Cmat)
   if (allocated(desc_sca_Cmat)) deallocate(desc_sca_Cmat)
   if (allocated(sca_phia)) deallocate(sca_phia)
   if (allocated(desc_sca_phia)) deallocate(desc_sca_phia)
   if (allocated(sca_ymat)) deallocate(sca_ymat)
   if (allocated(desc_sca_ymat)) deallocate(desc_sca_ymat)
   if (allocated(sca_ymat_copy)) deallocate(sca_ymat_copy)
   if (allocated(desc_sca_ymat_copy)) deallocate(desc_sca_ymat_copy)
   if (allocated(sca_ymat_qr_svd)) deallocate(sca_ymat_qr_svd)
   if (allocated(desc_sca_ymat_qr_svd)) deallocate(desc_sca_ymat_qr_svd)
   if (allocated(sca_w_params)) deallocate(sca_w_params)
   if (allocated(desc_sca_w_params)) deallocate(desc_sca_w_params)


   if (allocated(xdesc_emd)) deallocate(xdesc_emd)
   if (allocated(xdesc_fmd)) deallocate(xdesc_fmd)
   if (allocated(xdesc_smd)) deallocate(xdesc_smd)

end subroutine deallocate_fit
