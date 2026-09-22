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



subroutine ml()

  use mld_mpi
  use mld_logger
  use mld_string
  use mld_unit
  use ml_in_ndm_module, only: rangml, toy_model, ml_type, ml_type_krr, ml_type_basis, &
                              ml_type_descriptors,  &
                              mld_regularization_type, mld_regularization_type_home, train_only, &
                              marginal_likelihood,   &
                              ml_type_analysis,   &
                              ml_type_krr_old, &
                              prepare_factorial, write_desc, write_desc_dump, &
                              mld_fit_type, fit_online_svd
                              
  use module_optimization, only: optimize_weights_db, optimize_weights_chem
  use module_kernel, only: kernel_type, kernel_random, kernel_random_po, kernel_random_maha
  use module_ml_scalapack, only: scalapack_driver, context
  use module_kernel_2b, only: activate_k2b, dim_kernel_2b
  use main_mld_mod, only: main_train_mld, main_train_optimize_weights
  use module_end_ml, only: end_ml
  use module_test_mld, only: main_test_mld
  use module_online_fit, only: main_train_online_fit

  implicit none

  _NAMECURRENT_("ml")

  character(len=:), allocatable    :: tmp

  _MLD_BEGIN_

  tmp = '------------------------'//nwl
  call log_info(tmp//'--- MACHINE LEARNING ---'//nwl//tmp)
  call read_ml_file()

  call prepare_factorial()
  call periodic_table()
  if (scalapack_driver) then 
    call set_up_scalapack
  end if 
   
  if (toy_model) call set_toy_model

  ! init section ......... 
  call init_descriptors
  
  ! init the kernel if there is some ...
  dim_kernel_2b = 0
  if (activate_k2b) then 
    call init_kernel_k2b
    call log_info('ML: 2-body kernel is activated and the dimension is '//vtoa(dim_kernel_2b))
  end if 

  if (ml_type == ml_type_krr) then
    if (kernel_type == kernel_random) then
      call init_sample_kernel
    else if (kernel_type == kernel_random_po) then
      call init_sample_kernel_po
    else if (kernel_type == kernel_random_maha) then 
      call init_sample_kernel_maha  
    else
      call init_kernel
    end if
  end if

  call init_write_descriptors

  !end init section .... 


  if ((ml_type == ml_type_basis) .or. (ml_type == ml_type_krr)) then
    !  fit with regularization
    if (mld_regularization_type == mld_regularization_type_home) then
      call main_train_regularization()
    else
      if (mld_fit_type == fit_online_svd) then
        ! Online Gram-matrix fit: bypasses full Amat allocation
        call main_train_online_fit()
        if (.not. train_only) call main_test_mld
        call repport_time_train_test
      else if (marginal_likelihood) then
        call main_train_marginal_likelihood()
      else
        ! regular fit one shot ... the one used all the time in the last years
        call main_train_mld
      end if
      if (mld_fit_type /= fit_online_svd) then
        if ((optimize_weights_db).or.(optimize_weights_chem)) then 
          call main_train_optimize_weights
        end if    
        if (.not. train_only) call main_test_mld
        call repport_time_train_test
      end if
    end if
  end if
  if (ml_type == ml_type_descriptors) then 
    call main_compute_descriptors()
    call repport_time_train_test()
  end if   
  if (ml_type == ml_type_analysis)   then 
     call main_descriptors_analysis()
     ! call repport_time_train_test()

  end if 
!  <----------------------------------------

  if (ml_type == ml_type_krr_old) then
    !oldw! ! we read how many files are in the repository given by the variable path in order to provide  nd_data
    !oldw! if (build_subdata) then
    !oldw!   call buildsubdata(nd_data)
    !oldw! end if
    !oldw! if (lsoap .or. acd) then
    !oldw!   size_indi2 = 0
    !oldw!   max_ntyp = 0
    !oldw!   if (allocated(data_im)) deallocate (data_im)
    !oldw!   allocate (data_im(0:ns_data))
    !oldw!   data_im(0) = 0
    !oldw! end if
    !oldw! ! sparsification
    !oldw! if (allocated(reject)) deallocate (reject)
    !oldw! allocate (reject(ns_data))
    !oldw! reject(:) = .false.
    !oldw! if (sparsification) then
    !oldw!   call try_sparsification(n_count)
    !oldw! else
    !oldw!   n_count = ns_data
    !oldw! end if
    !oldw! ! computing descriptors
    !oldw! i_count = 0
    !oldw! call wesley_fill_desc(n_count)
    !oldw! if (lsoap) call wesley_fill_soap(n_count)
    !oldw! if (acd) call wesley_fill_acd(n_count)
    !oldw! ! up to here we fill xdesc
    !oldw! if (rescale) then
    !oldw!   do i = 1, dim_data
    !oldw!     xdesc(:, i) = (xdesc(:, i) - minval(xdesc))*(s_max_r - s_min_r)/(maxval(xdesc) - minval(xdesc)) + s_min_r + &
    !oldw!                   (0d0, 1d0)*((xdesc(:, i) - minval(xdesc))*(s_max_i - s_min_i)/(maxval(xdesc) - minval(xdesc)) + s_min_i)
    !oldw!   end do
    !oldw! end if
    !oldw! 
    !oldw! if (rangml == 0) write (*, *) "dim_xdesc", dim_xdesc, "dim_data", dim_data, "sub_data", n_count    ! ns_data

    !oldw! if (debug) then
    !oldw!   do i = 1, dim_data
    !oldw!     if (rangml == 0) write (91, '(2D20.8, 12D16.8)') yfunc(i), abs(xdesc(:, i))**2, xdesc(:, i)
    !oldw!   end do
    !oldw! end if

    !oldw! if (allocated(yfunc_average)) deallocate (yfunc_average); allocate (yfunc_average(dim_yfunc))
    !oldw! 
    !oldw! ! TODO this part is completely strage, it seems that the database is not recentred.
    !oldw! allocate (yfunc_nd(dim_yfunc, dim_data))
    !oldw! yfunc_nd(dim_yfunc, 1:dim_data) = yfunc(1:dim_data)
    !oldw! call recentrate_database(dim_data, dim_yfunc, yfunc_nd, yfunc_average)
    !oldw! yfunc(:) = yfunc_nd(dim_yfunc, :)
    !oldw! ! TODO end of the part completely strange. That now is recentred


    !oldw! ! entering kcross, dim_data ', dim_data
    !oldw! ! fix dim_data = dim_data_train + dim_data_test (nf * dim_data)
    !oldw! call allocate_xdesc(n_frac, rangml)

    !oldw! if (kcross) then
    !oldw!   call set_kcross(dim_data_train, rangml)
    !oldw!   if (search_hyp) then
    !oldw!     if (rangml == 0) then
    !oldw!       inquire (file="res", exist=exist)
    !oldw!       if (exist) then
    !oldw!         open (63, file="res", status="old", position="append", action="write")
    !oldw!       else
    !oldw!         open (63, file="res", status="new", action="write")
    !oldw!       end if
    !oldw!     end if
    !oldw!   end if                  ! search_hyp
    !oldw!   if ((rangml == 0) .and. (kernel_type == 1)) write (6, *) 'lambda_krr', lambda_krr, 'length_kse', length_kse
    !oldw!   if (rangml == 0) then
    !oldw!     inquire (file="corr.dat", exist=exist)
    !oldw!     if (exist) then
    !oldw!       open (62, file="corr.dat", status="unknown", position="rewind", action="write")
    !oldw!     else
    !oldw!       open (62, file="corr.dat", status="new", position="append", action="write")
    !oldw!     end if
    !oldw!   end if
    !oldw!   r_mean = 0d0
    !oldw!   do ik = 1, n_kcross
    !oldw!     ! the original data is splitted in train and validate ...
    !oldw!     ! dim_train, dim_valid are the dimension for the ik process of the n_kcross
    !oldw!     call fill_xdesc_yfunc_kcross(ik, rangml)
    !oldw!     call set_limit_for_cov(rangml, dim_train, i_final_cov, i_start_cov)
    !oldw!     call build_kernel_matrix(xdesc_train, dim_xdesc, dim_train)
    !oldw!     ! debug
    !oldw!     ! #if(PARA)debu
    !oldw!     !     if (allocated(eigen_values)) deallocate(eigen_values)
    !oldw!     !     allocate(eigen_values(dim_train))
    !oldw!     !     call diago_scalapack(dim_train,eigen_values)
    !oldw!     ! #endif
    !oldw!     call solve_kernel_matrix(rangml)
    !oldw!     if (allocated(y_extra)) deallocate (y_extra)
    !oldw!     allocate (y_extra(dim_valid))
    !oldw! 
    !oldw!     call krr_extrapolation(dim_valid, xdesc_valid, y_extra, error_valid)
    !oldw!     call correlation_coef(yfunc_valid, y_extra, dim_valid, r_coeff)
    !oldw!     if (rangml == 0) write (6, '(A,I2,1x,3(G16.8,1x))') 'ML:  kcross ....:', ik, lambda_krr, length_kse, r_coeff
    !oldw!     do i = 1, dim_valid
    !oldw!       if (rangml == 0) write (62, '(I6,1x,3(G16.8,1x))') i, yfunc_valid(i), y_extra(i), error_valid(i)
    !oldw!     end do
    !oldw!     r_mean = r_mean + r_coeff
    !oldw!     ! write (*,*) rangml,im, dim_train, i_start_at, i_final_at,i_final_cov,i_start_cov
    !oldw!   end do                  ! ik and nk_cross
    !oldw!   r_mean = r_mean/n_kcross
    !oldw!   if (rangml == 0) write (6, '(A,1x,G16.8)') 'r_mean :', r_mean
    !oldw!   if (search_hyp .and. (rangml == 0)) write (63, '(3(G16.8))') lambda_krr, length_kse, r_mean
    !oldw!   close (62)
    !oldw! end if                  ! kcross


    !oldw! if (search_hyp) close (63)

    !oldw! !$! if (marginal_likelihood) call min_marginal_likelihood()
    !oldw! if ((.not. marginal_likelihood) .and. (.not. kcross)) then
    !oldw!   dim_valid = dim_data_train*0.2
    !oldw!   dim_train = dim_data_train - dim_valid
    !oldw!   call fill_xdesc_yfunc_oneshot(rangml)
    !oldw!   call set_limit_for_cov(rangml, dim_train, i_final_cov, i_start_cov)
    !oldw!   call set_unlimit_for_atoms(rangml, im, i_start_at, i_final_at)
    !oldw!   ! call set_neighbors_size_cov_ml (i_start_cov,i_final_cov,i_local_cov, dim_data)
    !oldw!   call build_kernel_matrix(xdesc_train, dim_xdesc, dim_train)
    !oldw!   call solve_kernel_matrix(rangml)
    !oldw!   if (allocated(y_extra)) deallocate (y_extra); allocate (y_extra(dim_valid))
    !oldw!   call krr_extrapolation(dim_valid, xdesc_valid, y_extra, error_valid)
    !oldw!   call correlation_coef(yfunc_valid, y_extra, dim_valid, r_coeff)
    !oldw!   if (rangml == 0) write (6, '(A,3(G16.8,1x))') 'ML:  oneshot ....:', lambda_krr, length_kse, r_coeff
    !oldw!   do i = 1, dim_valid
    !oldw!     if (rangml == 0) write (62, '(I6,1x,3(G16.8,1x))') i, yfunc_valid(i), y_extra(i), error_valid(i)
    !oldw!   end do
    !oldw!   if (krr_error) then
    !oldw!     allocate (error_valid(dim_data))
    !oldw!     call krr_extrapolation(dim_data, xdesc, yfunc_valid, error_valid)
    !oldw!   end if
    !oldw!   if (rangml == 0) write (*, *) 'Correlation coefficient for (both github) validation', r_coeff
    !oldw! end if

    !oldw! if (n_frac .gt. 0d0) then
    !oldw!   if (allocated(y_test)) deallocate (y_test); allocate (y_test(dim_data_test))
    !oldw!   if (allocated(error_test)) deallocate (error_test); allocate (error_test(dim_data_test))
    !oldw!   y_test(:) = yfunc(dim_data_train + 1:dim_data)
    !oldw!   call krr_extrapolation(dim_data_test, xdesc_test, y_test, error_test)
    !oldw!   call correlation_coef(yfunc(dim_data_train + 1:dim_data), y_test, dim_data_test, r_coeff)
    !oldw!   if (rangml == 0) write (6, *) 'Correlation coefficient for prediction :', r_coeff
    !oldw!   do i = 1, dim_data_test
    !oldw!     if (rangml == 0) write (63, '(I4,3(G16.8))') i, yfunc(dim_data_train + i), y_test(i), error_test(i)
    !oldw!   end do
    !oldw! end if
    !oldw! call log_info('Do not forget young Jedi'//nwl//'the database recentred was'//nwl// &
    !oldw!               'and also the results corrected should (add this value, you should)'//nwl// &
    !oldw!               'yfunc_average = '//vtoa(yfunc_average(1:dim_yfunc)))
    call log_warning('THIS BRANCH WAS COMPLETELY REMOVED. PLEASE TRY OTHER OPTIONS')
  end if                  ! here get out from ml_type_krr

  if (write_desc .or. write_desc_dump) then
     call end_write_descriptors  
  end if
  call end_ml("----Full end----", scalapack_driver, context, rangml)
  _MLD_END_
end subroutine ml


subroutine repport_time_train_test
  use ml_in_ndm_module, only: rangml, ml_type, ml_type_krr, ml_type_descriptors, train_time, train_only
  use time_check_general, only: time_fit_params, time_full_desc, time_calc_desc, &
                                time_neigh_desc, time_fill_desc, time_fill_desc_02,&
                                test_time_total, &
                                test_time_desc, &
                                test_time_eval, &
                                test_time_neigh, debug_time, test_tot_time, train_tot_time, &
                                time_read_db, &
                                time_fill_desc_02A, time_fill_desc_02obj, time_fill_desc_02y, &
                                time_scatter_blacs, time_scatter_gather, time_scatter_bcast, &
                                time_scatter_dgemm, time_scatter_barrier
  use module_kernel, only: train_time_for_desc_kernel, &
                           test_time_for_desc_kernel
  use module_kernel_2b, only: activate_k2b, train_time_for_desc_kernel_2b, &
                              test_time_for_desc_kernel_2b                    
  implicit none

  if (rangml == 0) write (6, '("ML:------------------REPORTING TIME---------------------------")')
  call repport_time(2, 0.d0, time_read_db, "ML: Full DB read time ")
  
  if (ml_type == ml_type_descriptors) then 
    call repport_time(2, 0.d0, time_full_desc, "ML: Full TIME + writting")
    call repport_time(3, 0.d0, time_neigh_desc, "ML: time for searching neighbours ")
    call repport_time(3, 0.d0, time_calc_desc, "ML: time for computing desc")
    if (debug_time)  then
      train_time = .true.  
      call repport_debug_time_descriptor(train_tot_time)
    end if   
  else  
    call repport_time(2, 0.d0, time_full_desc, "ML: TRAIN TIME")
    call repport_time(3, 0.d0, time_neigh_desc, "ML: time for searching neighbours ")
    call repport_time(3, 0.d0, time_calc_desc, "ML: time for computing desc")
    if (debug_time)  then
      train_time=.true.  
      call repport_debug_time_descriptor(train_tot_time)
    end if   
    if (ml_type == ml_type_krr) then
      call repport_time(4, 0.d0, time_calc_desc - train_time_for_desc_kernel, "ML: pure desc")
      call repport_time(4, 0.d0, train_time_for_desc_kernel, "ML: kernel desc")
    end if
    if (activate_k2b)  call repport_time(3, 0.d0, train_time_for_desc_kernel_2b, "ML: kernel 2-body")
  
    if (debug_time) call repport_debug_time_kernel()
    call repport_time(3, 0.d0, time_fill_desc, "ML: time for fillig the local design matrix")
    call repport_time(3, 0.d0, time_fill_desc_02, "ML: time for filling the global design matrix")
    call repport_time(4, 0.d0, time_fill_desc_02y, "ML: filling  global y")
    call repport_time(4, 0.d0, time_fill_desc_02A, "ML: filling  global Amat")
    call repport_time(4, 0.d0, time_fill_desc_02obj, "ML: filling  global obj")

    if (time_scatter_blacs > 0.d0) then
      call repport_time(3, 0.d0, time_scatter_blacs, "ML: time scattering config to BLACS (online fit)")
      call repport_time(4, 0.d0, time_scatter_gather, "ML: scatter: mc_per_group gather")
      call repport_time(4, 0.d0, time_scatter_bcast, "ML: scatter: MPI_Bcast(Ac_full)")
      call repport_time(4, 0.d0, time_scatter_dgemm, "ML: scatter: local extract + dgemm")
      call repport_time(4, 0.d0, time_scatter_barrier, "ML: scatter: barrier")
    end if

    call repport_time(3, 0.d0, time_fit_params, "ML: time for fitting params")
    if(.not.(train_only)) then 
    call repport_time(2, 0.d0, test_time_total, "ML: TEST TIME")
    call repport_time(3, 0.d0, test_time_neigh, "ML: time for searching neighbours ")
    call repport_time(3, 0.d0, test_time_desc, "ML: time for computing desc")
    end if 
    if (debug_time) then 
      train_time = .false.   
      if(.not.(train_only)) call repport_debug_time_descriptor(test_tot_time)
    end if 
      
    if (ml_type == ml_type_krr) then
      if(.not.(train_only)) call repport_time(4, 0.d0, test_time_desc - test_time_for_desc_kernel, "ML: pure desc")
      if(.not.(train_only)) call repport_time(4, 0.d0, test_time_for_desc_kernel, "ML: kernel desc")
    end if
    if(.not.(train_only)) then 
    if (activate_k2b)  call repport_time(3, 0.d0, test_time_for_desc_kernel_2b, "ML: kernel 2-body")
    end if 
    if (debug_time) call repport_debug_time_kernel()
  
    if(.not.(train_only)) call repport_time(3, 0.d0, test_time_eval, "ML: time for evaluation ")
  end if 
end subroutine repport_time_train_test



subroutine repport_debug_time_descriptor(tot_time)
  use module_kind_variables, only: kind_double
  use ml_in_ndm_module, only: descriptor_type, descriptor_bispectrum_so4, descriptor_milady, & 
                              descriptor_ftnbody, train_time
  use module_chemical_species, only: tcc01_rdist, tcc02_rdist, tcc03_rdist, tnn_rdist, tii_rdist
  use module_ftnbody, only: tnn_ftbd, t2b_ftbd, t3b_ftbd, t4b_ftbd, t5b_ftbd, &
                            tnn_ftbd_train, t2b_ftbd_train, t3b_ftbd_train, t4b_ftbd_train, t5b_ftbd_train, & 
                            t4b_inner_init, t4b_inner_desc, t4b_inner_deriv, t4b_inner_last
  implicit none 
  real(kind_double) :: timet
  real(kind_double), dimension(20), intent(in) :: tot_time 
  real(kind_double) :: l_tnn_ftbd, l_t2b_ftbd, l_t3b_ftbd, l_t4b_ftbd, l_t5b_ftbd

  if (descriptor_type == descriptor_bispectrum_so4) then 
    !call repport_time(2, 0.d0, tot_time(1), "ML: full")
    call repport_time(4, 0.d0, tot_time(2), "ML: Ujmm ")
    call repport_time(4, 0.d0, tot_time(4), "ML: Ujmm + nn")
    call repport_time(4, 0.d0, tot_time(8), "ML: time pack")
    call repport_time(4, 0.d0, tot_time(5), "ML: CG last")
    call repport_time(4, 0.d0, tot_time(6), "ML: inner CG ")
  end if   
  if (descriptor_type == descriptor_milady) then 
    timet= tcc01_rdist + tcc02_rdist + tcc03_rdist + tnn_rdist + tii_rdist
    call repport_time(4, 0.d0, timet, "ML: Total computing  ")
    call repport_time(4, 0.d0, tnn_rdist, "ML: nn search")
    call repport_time(4, 0.d0, tcc01_rdist, "ML: calc ja")
    call repport_time(4, 0.d0, tcc02_rdist, "ML: calc loop ia")
    call repport_time(4, 0.d0, tii_rdist, "ML: ordering ")
    call repport_time(4, 0.d0, tcc03_rdist, "ML: final ")

  end if   

  if (descriptor_type == descriptor_ftnbody) then 
    if (train_time) then 
      l_tnn_ftbd=tnn_ftbd_train
      l_t2b_ftbd=t2b_ftbd_train
      l_t3b_ftbd=t3b_ftbd_train
      l_t4b_ftbd=t4b_ftbd_train
      l_t5b_ftbd=t5b_ftbd_train
    else
      l_tnn_ftbd=tnn_ftbd
      l_t2b_ftbd=t2b_ftbd
      l_t3b_ftbd=t3b_ftbd
      l_t4b_ftbd=t4b_ftbd
      l_t5b_ftbd=t5b_ftbd
    end if   
    
    call repport_time(4, 0.d0, l_tnn_ftbd, "ML: nn search")
    call repport_time(4, 0.d0, l_t2b_ftbd, "ML: 1b and 2b time")
    call repport_time(4, 0.d0, l_t3b_ftbd, "ML: 3b time")
    call repport_time(4, 0.d0, l_t4b_ftbd, "ML: 4b time")
    call repport_time(4, 0.d0, l_t5b_ftbd, "ML: 5b time")  
    if (train_time) then
      call repport_time(3, 0.d0, t4b_inner_init, "ML: 4b init") 
      call repport_time(3, 0.d0, t4b_inner_desc, "ML: 4b desc") 
      call repport_time(3, 0.d0, t4b_inner_deriv, "ML: 4b deriv") 
      call repport_time(3, 0.d0, t4b_inner_last, "ML: 4b last") 
    end if
  end if   


end subroutine repport_debug_time_descriptor


subroutine repport_debug_time_kernel()
  use module_kernel_2b, only: activate_k2b, tnn_2b , tuu_2b, tkk_2b
  implicit none 
  if (activate_k2b) then    
    call repport_time(4, 0.d0, tnn_2b, "ML: searching nn")
    call repport_time(4, 0.d0, tuu_2b, "ML: utilities for nn")
    call repport_time(4, 0.d0, tkk_2b, "ML: kernel evaluation")
  end if   

end subroutine repport_debug_time_kernel 
