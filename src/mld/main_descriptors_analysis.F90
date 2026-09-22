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

subroutine main_descriptors_analysis


  !use mpi
  use mld_mpi, only: mld_mpi_abort, mld_size
  use ml_in_ndm_module, only: rangml, n_pca, classes_full_for_sigma, &
                              classes_train_for_sigma, desc_forces
  use module_db_setup, only: iconf_data, db_path
  use main_mld_mod, only: read_parameters_for_md 
  !use cholesky
  use temporary_data_cov, only: dim_xdesc
  use module_Sigma_matrix, only: first_passage_mcd, right_vector_with_energy, & 
                                 Phi_mat, Sigma_sample_mcd, Sigma_sample_full, &
                                 Sigma_sample_mcd_inv, Sigma_sample_full_inv, & 
                                 Sigma_sample, Sigma_sample_inv,  & 
                                 mean, mean_full, mean_mcd, & 
                                 eigenvalues_Sigma_sample_mcd, eigenvectors_Sigma_sample_mcd, &
                                 eigenvalues_Sigma_sample_full, eigenvectors_Sigma_sample_full, & 
                                 dim_line_design, dim_col_design, & 
                                 compute_Sigma_class_by_sum, &
                                 analysis_local_energy, & 
                                 method_sigma, method_sigma_by_Phi, method_sigma_by_sum
  use module_Phi, only: dimr_Phi, dimc_Phi, compute_serial_Phi, compute_Sigma_class_by_Phi, &
                        get_dimension_design_Phi_by_class  
  use module_kernel, only: write_kernel_matrix, kernel_dump, kernel_dump_by_mahalanobis_norm, &
                           kernel_dump_by_mahalanobis,  kernel_dump_by_cur, kernel_dump_by_cur_maha, &
                           info_kernel_full, info_kernel_ref, info_kernel, &
                           np_kernel_full, np_kernel_ref
  use module_write_parameters, only: dump_kernel_matrix
  use mld_unit
  use module_create_descdb, only: create_descdb
  use module_ml_scalapack, only: sca_phia, desc_sca_phia, dimr_sca_phia, dimc_sca_phia, & 
                                  nbr_phia, nbc_phia, l_dimr_sca_phia, l_dimc_sca_phia, &
                                  context, scalapack_driver
  use module_end_ml, only: end_ml
  use module_compute_cur, only: compute_cur_decomposition, compute_cur_maha_decomposition, compute_sca_phia 
  use module_cur, only :  info_mat, pcol, selcol, prow, selrow, &
                          col_no_of_selections, row_no_of_selections, &
                          cur_info_selection, &
                          cur_kval, cur_cval, cur_rval, cur_eps

  use math, only: diago_serial_real_double, serial_pseudo_inverse
  use mld_logger, only: log_info, log_critical
  use time_check_general, only: time, tot_time, MY_MPI_WTIME
  use module_db_poscar, only: i_start_conf, i_final_conf, procs_per_file
  use set_limits, only: set_limit_for_configs_with_MPI_grid
  use mld_subworld
  implicit none

  integer  :: ii, itmp, iref, nounit
  character(len=80)    :: namefile, text 
  logical  :: post_desc 
  _NAMECURRENT_("main_descriptor_analysis")

  _MLD_BEGIN_

  if (procs_per_file /= mld_size) then
    call log_warning("ML: procs_per_file is not equal to mld_size. This is not good for analysis mode. YET!")
    call log_warning("ML: procs_per_file has the value of "//trim(vtoa(procs_per_file))//" and mld_size has the value of "//trim(vtoa(mld_size)))
    call log_warning("ML: procs_per_file is set to mld_size.")
    procs_per_file = mld_size
  end if 
            
  analysis_local_energy=.false. 
  !creation of descDB directory. I do not understand why this is needed ... 
  time(1) = MY_MPI_WTIME()

  call prepare_database
  call create_descdb("desc", db_path)
  desc_forces= .false. 
  first_passage_mcd = .false.
     !set-up on which proc is read each conf ...

  call init_subworld(procs_per_file)
  call set_limit_for_configs_with_MPI_grid(iconf_data, i_start_conf, i_final_conf)
  call prepare_train_dimensions
  call prepare_test_dimensions
 
  call log_info("ML: in descriptors_analysis mode the descriptors force calculations are switch to .false.")  
  call log_info("ML: descriptors_analysis dim_xdesc = "//trim(vtoa(dim_xdesc)))
  ! in this subroutine the desc are computed on all procs but they are pack one one proc.
  time(2) = MY_MPI_WTIME()
  tot_time(2) = tot_time(2) + time(2) - time(1)  
  post_desc=.false.
  call compute_descriptors_all_database(post_desc)
  time(3) = MY_MPI_WTIME()
  tot_time(3) = tot_time(3) + time(3) - time(2)  

  call log_info ("ML:>>>------------All descriptors are computed -------------------------------")


  !TODO - TO remove this shit or add an option to compute this shit. 
  !TODO - desactivate by default.  
  if (analysis_local_energy) call read_parameters_for_md 
  
  if (analysis_local_energy) then 
    ! compute  the vector phi(D_ia)^T x E_ia the local energy for each atom 
    ! and all packer into right_vector_with_energy/compute_right_vector_with_energy
    if (allocated(right_vector_with_energy)) deallocate (right_vector_with_energy); 
    allocate (right_vector_with_energy(dim_xdesc))
    call compute_right_vector_with_energy(right_vector_with_energy, size(right_vector_with_energy, 1))
  end if 
  ! classes_train_for_sigma(:) =(/ '07', '08', '12'/)
  call log_info("ML: descriptors_analysis dim_line_design = "//trim(vtoa(dim_line_design))// &
               ", dim_col_design = "//trim(vtoa(dim_col_design))// &
               " => Phi_mat ≈ "//trim(vtoa(dim_line_design*dim_col_design*8/1024/1024))//" MB")

  call log_info("ML:>>>------------ Compute Sigma / stat dist - train  -------------------------")

  select case (method_sigma)
    case (method_sigma_by_Phi)
      ! THis is not scalapack ... not good :(
      call compute_serial_Phi(classes_train_for_sigma, dim_line_design, dim_col_design, mean, Phi_mat)
      call compute_Sigma_class_by_Phi(Phi_mat, Sigma_sample, Sigma_sample_inv)
    case (method_sigma_by_sum)
      call get_dimension_design_Phi_by_class(classes_train_for_sigma, mean, dimr_Phi, dimc_Phi)
      call compute_Sigma_class_by_sum(classes_train_for_sigma, size(classes_train_for_sigma), dimr_Phi)
    case default
      call log_critical("ML first passage will stop: method_sigma can be 1 or 2. Fix this in ml file.") 
      call mld_mpi_abort("stop in"//NAMECURRENT) 
  end select    
  
  !call  scalapack_pseudo_inverse (sca_AA, desc_sca_AA, & 
  !  nbr_AA, nbc_AA, dimr_sca_AA, dimc_sca_AA, &
  !  l_dimr_sca_AA, l_dimc_sca_AA, &
  !  rank_sca_AA, &
  !  sca_pinvAA, desc_sca_pinvAA, & 
  !  nbr_pinvAA, nbc_pinvAA, dimr_sca_pinvAA, dimc_sca_pinvAA, &
  !  l_dimr_sca_pinvAA, l_dimc_sca_pinvAA)
  call log_info("ML: dimr_Phi = "//trim(vtoa(dimr_Phi))// &
               " => covariance ≈ "//trim(vtoa(dimr_Phi*dimr_Phi*8/1024/1024))//" MB")
  if (allocated(Sigma_sample_mcd)) deallocate (Sigma_sample_mcd) 
  allocate (Sigma_sample_mcd(dimr_Phi, dimr_Phi))
  
  if (allocated(Sigma_sample_mcd_inv)) deallocate (Sigma_sample_mcd_inv) 
  allocate (Sigma_sample_mcd_inv(dimr_Phi, dimr_Phi))
  
  if (allocated(mean_mcd)) deallocate (mean_mcd); allocate (mean_mcd(dimr_Phi))
 
  Sigma_sample_mcd(:, :) = Sigma_sample(:, :)
  Sigma_sample_mcd_inv(:, :) = Sigma_sample_inv(:, :)
  mean_mcd(:) = mean(:)
  deallocate (Sigma_sample, Sigma_sample_inv, mean)
  if (method_sigma == method_sigma_by_Phi) deallocate(Phi_mat)

  first_passage_mcd = .true.
  call compute_mcd_dist(Sigma_sample_mcd_inv, size(Sigma_sample_mcd_inv, 1), mean_mcd)

  namefile = "train_sigma.mcd"
  open (newunit=nounit, file=trim(namefile), status='unknown', form='formatted')
  call write_mcd_distances(nounit)

  namefile = "train_sigma_selected.mcd"
  open (newunit=nounit, file=trim(namefile), status='unknown', form='formatted')
  call write_mcd_distances_for_some_classes(classes_train_for_sigma, size(classes_train_for_sigma), nounit)

  time(4) = MY_MPI_WTIME()
  tot_time(4) = tot_time(4) + time(4) - time(3)  


  call log_info("ML:>>>---------- Compute Sigma / stat dist - test/full  -----------------------")
  select case(method_sigma)
    case(method_sigma_by_Phi)
      call compute_serial_Phi(classes_full_for_sigma, dim_line_design, dim_col_design, mean, Phi_mat)
      call compute_Sigma_class_by_Phi( Phi_mat, Sigma_sample, Sigma_sample_inv)
    case(method_sigma_by_sum)
      call get_dimension_design_Phi_by_class(classes_full_for_sigma, mean, dimr_Phi, dimc_Phi)
      call compute_Sigma_class_by_sum(classes_full_for_sigma, size(classes_full_for_sigma), dimr_Phi)
    case default
      call log_critical("ML second passage will stop: methpd_sigma can be 1 or 2. Fix this in ml file.")
      call mld_mpi_abort("stop in "//NAMECURRENT) 
  end select
  
  if (allocated(Sigma_sample_full)) deallocate (Sigma_sample_full); allocate (Sigma_sample_full(dimr_Phi, dimr_Phi))
  if (allocated(Sigma_sample_full_inv)) deallocate (Sigma_sample_full_inv); allocate (Sigma_sample_full_inv(dimr_Phi, dimr_Phi))
  if (allocated(mean_full)) deallocate (mean_full); allocate (mean_full(dimr_Phi))
  Sigma_sample_full(:, :) = Sigma_sample(:, :)
  Sigma_sample_full_inv(:, :) = Sigma_sample_inv(:, :)
  mean_full(:) = mean(:)
  deallocate (Sigma_sample, Sigma_sample_inv, mean)
  if (method_sigma == method_sigma_by_Phi) deallocate(Phi_mat)

  first_passage_mcd = .false.
  call log_info("ML:>>>---------- Compute full Sigma MCD dist  -----------------------")
  call compute_mcd_dist(Sigma_sample_full_inv, size(Sigma_sample_full_inv, 1), mean_full)
  namefile = "full_sigma.mcd"
  open (newunit=nounit, file=trim(namefile), status='unknown', form='formatted')
  call write_mcd_distances(nounit)

  namefile = "full_sigma_selected.mcd"
  open (newunit=nounit, file=trim(namefile), status='unknown', form='formatted')
  call write_mcd_distances_for_some_classes(classes_train_for_sigma, size(classes_train_for_sigma), nounit)

  time(5) = MY_MPI_WTIME()
  tot_time(5) = tot_time(5) + time(5) - time(4)  


  ! SPECTRUM Sigma matrices
  if (allocated(eigenvalues_Sigma_sample_mcd)) deallocate (eigenvalues_Sigma_sample_mcd)
  allocate (eigenvalues_Sigma_sample_mcd(size(Sigma_sample_mcd, 1)))
  if (allocated(eigenvectors_Sigma_sample_mcd)) deallocate (eigenvectors_Sigma_sample_mcd) 
  allocate (eigenvectors_Sigma_sample_mcd(size(Sigma_sample_mcd, 1), size(Sigma_sample_mcd, 2)))

  call log_info("ML:>>>-------------- Diago  Sigma_sample ----------------------------")
  call diago_serial_real_double(Sigma_sample_mcd,  eigenvectors_Sigma_sample_mcd, eigenvalues_Sigma_sample_mcd)

  if (allocated(eigenvalues_Sigma_sample_full)) deallocate (eigenvalues_Sigma_sample_full); allocate (eigenvalues_Sigma_sample_full(size(Sigma_sample_full, 1)))
  if (allocated(eigenvectors_Sigma_sample_full)) deallocate (eigenvectors_Sigma_sample_full); allocate (eigenvectors_Sigma_sample_full(size(Sigma_sample_full, 1), size(Sigma_sample_full, 2)))
  call log_info("ML:>>>------------ Diago  Sigma_sample_full --------------------------")
  call diago_serial_real_double(Sigma_sample_full,  eigenvectors_Sigma_sample_full, eigenvalues_Sigma_sample_full)

  if (rangml == 0) then
    namefile = "train_sigma.eigenvalues"
    open (newunit=nounit, file=trim(namefile), status='unknown', form='formatted')
    do ii = 1, size(Sigma_sample_mcd, 1)
      write (nounit, '(e20.10)') eigenvalues_Sigma_sample_mcd(ii)
    end do
    close (nounit, status='keep')
  end if

  if (rangml == 0) then
    namefile = "full_sigma.eigenvalues"
    open (newunit=nounit, file=trim(namefile), status='unknown', form='formatted')
    do ii = 1, size(Sigma_sample_full, 1)
      write (nounit, '(e20.10)') eigenvalues_Sigma_sample_full(ii)
    end do
    close (nounit, status='keep')
  end if
  ! End SPECTRUM Sigma matrices

  time(6) = MY_MPI_WTIME()
  tot_time(6) = tot_time(6) + time(6) - time(5)  

  ! PCA Standard
  if (rangml ==0) then 
    namefile = "full_sigma.pca"
  call log_info("ML:>>>------------ Get PCA full  --------------------------")
    open (newunit=nounit, file=trim(namefile), status='unknown', form='formatted')
    call get_pca_projection(n_pca, mean_full, size(mean_full, 1), eigenvectors_Sigma_sample_full, &
                          size(eigenvalues_Sigma_sample_full), nounit)
    close (nounit, status='keep')
  end if
  
  if (rangml ==0) then 
    namefile = "train_sigma.pca"
    open (newunit=nounit, file=trim(namefile), status='unknown', form='formatted')
  call log_info("ML:>>>------------ Get PCA sample  -------------------------")
    call get_pca_projection(n_pca, mean_mcd, size(mean_mcd, 1), eigenvectors_Sigma_sample_mcd, &
                          size(eigenvalues_Sigma_sample_mcd), nounit)
    close (nounit, status='keep')
  end if 
  ! END PCA Standard
  if (analysis_local_energy) then 
    ! Write_funny distances
    if (rangml ==0) then 
      namefile = "train_sigma.stat"
      open (newunit=nounit, file=trim(namefile), status='unknown', form='formatted')
      call new_distances(eigenvalues_Sigma_sample_mcd, eigenvectors_Sigma_sample_mcd, mean_mcd, &
                       size(eigenvalues_Sigma_sample_mcd, 1), nounit)
      close (nounit, status='keep')
    end if 
  
    if (rangml==0) then 
      namefile = "full_sigma.stat"
      open (newunit=nounit, file=trim(namefile), status='unknown', form='formatted')
      call new_distances(eigenvalues_Sigma_sample_full, eigenvectors_Sigma_sample_full, mean_full, &
                       size(eigenvalues_Sigma_sample_full, 1), nounit)
      close (nounit, status='keep')
    end if
  end if 

  time(7) = MY_MPI_WTIME()
  tot_time(7) = tot_time(7) + time(7) - time(6)  

  ! Write_funny kernels ... 
  select case (kernel_dump)
    case (kernel_dump_by_mahalanobis_norm)
      call log_info("ML: The kernel is selected from MAHA distances normalized")
      call get_index_kernel_from_dist_normalized()
    case (kernel_dump_by_mahalanobis)
      call log_info("ML: The kernel is selected from draft MAHA distances")
      call get_index_kernel_from_dist()
    case (kernel_dump_by_cur)
      !call get_index_kernel_from_dist_normalized()
      call log_info("ML: The kernel is selected from CUR analysis")
      if (.not.(scalapack_driver)) then 
        call log_warning("ML: scalapack_driver is .false. for CUR decomposition. Will be set to .true.")
        scalapack_driver = .true. 
       
        call set_up_scalapack               
      end if 

      call log_info("ML: The Phia matrix for inner classes")
      call compute_sca_phia(classes_train_for_sigma, mean)
      cur_cval = np_kernel_ref
      call log_info("ML: The first CUR for inner kernel")
      call compute_cur_decomposition(sca_phia, desc_sca_phia, dimr_sca_phia, dimc_sca_phia, &
                                     l_dimr_sca_phia, l_dimc_sca_phia, nbr_phia, nbc_phia, &
                                     cur_kval, cur_cval, cur_rval, cur_eps, &
                                     pcol, selcol, col_no_of_selections,  prow, selrow, row_no_of_selections, &
                                     cur_info_selection)
      call log_info("ML: After the first CUR ...")
      
      if (allocated(info_kernel_ref)) deallocate(info_kernel_ref)
      allocate(info_kernel_ref(col_no_of_selections))
      do ii=1, col_no_of_selections
        info_kernel_ref(ii)%iconf = info_mat(cur_info_selection(ii))%iconf
        info_kernel_ref(ii)%ia =    info_mat(cur_info_selection(ii))%ia
      end do 
     
      call log_info("ML: The Phia matrix for full classes")
      call compute_sca_phia(classes_full_for_sigma, mean)
      cur_cval = np_kernel_full
      call log_info("ML: The second CUR for the full kernel")
      call compute_cur_decomposition(sca_phia, desc_sca_phia, dimr_sca_phia, dimc_sca_phia, &
                                     l_dimr_sca_phia, l_dimc_sca_phia, nbr_phia, nbc_phia, &
                                     cur_kval, cur_cval, cur_rval, cur_eps, &
                                     pcol, selcol, col_no_of_selections,  prow, selrow, row_no_of_selections, &
                                     cur_info_selection)
      
      if (allocated(info_kernel_full)) deallocate(info_kernel_full)
      allocate(info_kernel_full(col_no_of_selections))
      do ii=1, col_no_of_selections
        info_kernel_full(ii)%iconf = info_mat(cur_info_selection(ii))%iconf
        info_kernel_full(ii)%ia =    info_mat(cur_info_selection(ii))%ia
      end do 
      !TODO kernel remove points in double counting !!!!!
      itmp = size(info_kernel_ref) + size(info_kernel_full)
      if (allocated(info_kernel)) deallocate(info_kernel)
      allocate(info_kernel(itmp))
      do ii=1, size(info_kernel_ref)
        info_kernel(ii)%iconf = info_kernel_ref(ii)%iconf
        info_kernel(ii)%ia =    info_kernel_ref(ii)%ia
      end do 
      
      iref = size(info_kernel_ref)
      do ii=size(info_kernel_ref)+1, size(info_kernel_ref) +size(info_kernel_full)
        info_kernel(ii)%iconf = info_kernel_full(ii-iref)%iconf
        info_kernel(ii)%ia =    info_kernel_full(ii-iref)%ia
      end do 
      deallocate(info_kernel_full, info_kernel_ref)

    case (kernel_dump_by_cur_maha)
      !call get_index_kernel_from_dist_normalized()
      call log_info("ML: The kernel is selected from CUR analysis using MAHA sampling")
      if (.not.(scalapack_driver)) then 
        call log_warning("ML: scalapack_driver is .false. for CUR decomposition. Will be set to .true.")
        scalapack_driver = .true. 
        call set_up_scalapack               
      end if 

      call compute_sca_phia(classes_train_for_sigma, mean)
      cur_cval = np_kernel_ref
      call compute_cur_maha_decomposition(sca_phia, desc_sca_phia, dimr_sca_phia, dimc_sca_phia, &
                                     l_dimr_sca_phia, l_dimc_sca_phia, nbr_phia, nbc_phia, &
                                     cur_kval, cur_cval, cur_rval, cur_eps, &
                                     pcol, selcol, col_no_of_selections,  prow, selrow, row_no_of_selections, &
                                     cur_info_selection)
      
      if (allocated(info_kernel_ref)) deallocate(info_kernel_ref)
      allocate(info_kernel_ref(col_no_of_selections))
      do ii=1, col_no_of_selections
        info_kernel_ref(ii)%iconf = info_mat(cur_info_selection(ii))%iconf
        info_kernel_ref(ii)%ia =    info_mat(cur_info_selection(ii))%ia
      end do 

      call compute_sca_phia(classes_full_for_sigma, mean)
      cur_cval = np_kernel_full
      call compute_cur_maha_decomposition(sca_phia, desc_sca_phia, dimr_sca_phia, dimc_sca_phia, &
                                     l_dimr_sca_phia, l_dimc_sca_phia, nbr_phia, nbc_phia, &
                                     cur_kval, cur_cval, cur_rval, cur_eps, &
                                     pcol, selcol, col_no_of_selections,  prow, selrow, row_no_of_selections, &
                                     cur_info_selection)
      
      if (allocated(info_kernel_full)) deallocate(info_kernel_full)
      allocate(info_kernel_full(col_no_of_selections))
      do ii=1, col_no_of_selections
        info_kernel_full(ii)%iconf = info_mat(cur_info_selection(ii))%iconf
        info_kernel_full(ii)%ia =    info_mat(cur_info_selection(ii))%ia
      end do 
      !TODO kernel remove points in double counting !!!!!
      itmp = size(info_kernel_ref) + size(info_kernel_full)
      if (allocated(info_kernel)) deallocate(info_kernel)
      allocate(info_kernel(itmp))
      do ii=1, size(info_kernel_ref)
        info_kernel(ii)%iconf = info_kernel_ref(ii)%iconf
        info_kernel(ii)%ia =    info_kernel_ref(ii)%ia
      end do 
      
      iref = size(info_kernel_ref)
      do ii=size(info_kernel_ref)+1, size(info_kernel_ref) +size(info_kernel_full)
        info_kernel(ii)%iconf = info_kernel_full(ii-iref)%iconf
        info_kernel(ii)%ia =    info_kernel_full(ii-iref)%ia
      end do 
      deallocate(info_kernel_full, info_kernel_ref)

    case default
      write (6, '("ML: this selection of kernel is not implemented, kernel_dump is set to...",i6)') kernel_dump
      text = " --- ERROR on kernel dump in "//NAMECURRENT//" ERROR ---"
      call end_ml(text, scalapack_driver, context, rangml)
  end select

  if (write_kernel_matrix) then
    call dump_kernel_matrix

    namefile = "kernel_full_sigma.pca"
    if (rangml==0) open (newunit=nounit, file=trim(namefile), status='unknown', form='formatted')
    call get_pca_projection_kernel(n_pca, mean_full, size(mean_full, 1), eigenvectors_Sigma_sample_full, &
                                    size(eigenvalues_Sigma_sample_full), nounit)
    if (rangml==0) close (nounit, status='keep')

    namefile = "kernel_train_sigma.pca"
    if (rangml==0) open (newunit=nounit, file=trim(namefile), status='unknown', form='formatted')
    call get_pca_projection_kernel(n_pca, mean_mcd, size(mean_mcd, 1), eigenvectors_Sigma_sample_mcd, &
                                   size(eigenvalues_Sigma_sample_mcd),nounit)
    if (rangml==0) close (nounit, status='keep')
  end if


  time(10) = MY_MPI_WTIME()
  tot_time(1) = tot_time(1) + time(10) - time(1)
  tot_time(8) = tot_time(8) + time(10) - time(7)  


  call repport_time(2, 0.d0, tot_time(1), "ML: analysis total ")
  call repport_time(4, 0.d0, tot_time(2), "ML: init precomputing ")
  call repport_time(4, 0.d0, tot_time(3), "ML: init compute desc ")
  call repport_time(4, 0.d0, tot_time(4), "ML: Phi Sigma  train  ")
  call repport_time(4, 0.d0, tot_time(5), "ML: Phi Sigma  full   ")
  call repport_time(4, 0.d0, tot_time(6), "ML: Sigma spectrum    ")
  call repport_time(4, 0.d0, tot_time(7), "ML: PCA and dist      ")
  call repport_time(4, 0.d0, tot_time(8), "ML: Kernel selection  ")
  !call repport_time(6, 0.d0, tot_time(10), "ML: spher 3D :    cmm")
  

  if (rangml == 0) write (6, '("-------------------")')
  if (rangml == 0) write (6, '(" xxx  ")')

 _MLD_END_
end subroutine main_descriptors_analysis



subroutine get_index_kernel_from_dist()

  use ml_in_ndm_module, only: classes_train_for_sigma, rangml
  use module_db_setup, only: iconf_data 
  use derived_types, only: config_desc, config_real
  use module_kernel, only: info_kernel, info_kernel_copy, np_kernel_ref, np_kernel_full, power_mcd
  use mld_logger

  implicit none

  integer  :: ic, ia, icnt, iclass, ii, indx
  real(kind(1.d0)), dimension(:), allocatable  :: dist_mcd_full, dist_maha_full
  integer, dimension(:), allocatable     :: index_mcd_full, index_maha_full, ic_vector_full, ia_vector_full

  real(kind(1.d0)), dimension(:), allocatable  :: dist_mcd_train, dist_maha_train
  integer, dimension(:), allocatable     :: index_mcd_train, index_maha_train, ic_vector_train, ia_vector_train
  integer  :: np_train, np_full, np_train_tmp, np_full_tmp, size_full, size_train
  real(kind(1.d0))     :: dint, ddmin
  logical, dimension(:), allocatable     :: lhigh, llow, lvec

  _NAMECURRENT_("get_index_kernel_from_dist")



  _MLD_BEGIN_
  ! begin pack distances for train classes
  np_train_tmp = np_kernel_ref
  np_full_tmp = np_kernel_full
  if (allocated(info_kernel)) deallocate (info_kernel)
  allocate (info_kernel(np_full_tmp + np_train_tmp))


  icnt = 0
  do iclass = 1, size(classes_train_for_sigma)
    do ic = 1, iconf_data
      if ((classes_train_for_sigma(iclass) == config_real(ic)%class)) then
        if (.not. (allocated(config_desc(ic)%energy))) cycle
        do ia = 1, config_real(ic)%nat
          icnt = icnt + 1
        end do
      end if
    end do
  end do

  size_train = icnt

  if (allocated(dist_mcd_train)) deallocate (dist_mcd_train)
  allocate (dist_mcd_train(size_train))
  if (allocated(dist_maha_train)) deallocate (dist_maha_train)
  allocate (dist_maha_train(size_train))
  if (allocated(ic_vector_train)) deallocate (ic_vector_train)
  allocate (ic_vector_train(size_train))
  if (allocated(ia_vector_train)) deallocate (ia_vector_train)
  allocate (ia_vector_train(size_train))

  icnt = 0
  do iclass = 1, size(classes_train_for_sigma)
    do ic = 1, iconf_data
      if ((classes_train_for_sigma(iclass) == config_real(ic)%class)) then
        if (.not. (allocated(config_desc(ic)%energy))) cycle
        do ia = 1, config_real(ic)%nat
          icnt = icnt + 1
          dist_mcd_train(icnt) = config_desc(ic)%stat_dist_mcd(ia)
          dist_maha_train(icnt) = config_desc(ic)%stat_dist_maha(ia)**power_mcd
          ic_vector_train(icnt) = ic
          ia_vector_train(icnt) = ia
        end do
      end if
    end do
  end do


  if (allocated(index_mcd_train)) deallocate (index_mcd_train); allocate (index_mcd_train(size_train))
  if (allocated(index_maha_train)) deallocate (index_maha_train); allocate (index_maha_train(size_train))

  call indexx(icnt, dist_mcd_train, index_mcd_train)
  call indexx(icnt, dist_maha_train, index_maha_train)

  if (allocated(lhigh)) deallocate (lhigh); allocate (lhigh(size_train))
  if (allocated(llow)) deallocate (llow); allocate (llow(size_train))
  if (allocated(lvec)) deallocate (lvec); allocate (lvec(size_train))

  ddmin = minval(dist_maha_train)
  dint = (maxval(dist_maha_train) - ddmin)/dble(np_train_tmp)
  ! dint = (maxval(dist_maha_train) - ddmin)/dble(np_train_tmp - 1)

  icnt = 0
  do ii = 0, np_train_tmp - 1
    lhigh(:) = (dist_maha_train >= (dble(ii + 1)*dint + ddmin))
    llow(:) = (dist_maha_train <= (dble(ii)*dint + ddmin))
    lvec(:) = (llow(:) .eqv. lhigh(:))
    ! write (677,*) ii, count(lvec .eqv. .true., dim=1),(dble(ii)*dint + ddmin), (dble(ii+1)*dint + ddmin)
    if (count(lvec .eqv. .true., dim=1) > 0) then
      icnt = icnt + 1
      indx = findloc(lvec, .true., dim=1)
      info_kernel(icnt)%iconf = ic_vector_train(indx)
      info_kernel(icnt)%ia = ia_vector_train(indx)
      ! write (677,*) indx, ic_vector_train(indx), ia_vector_train(indx)
    end if
  end do

  np_train = icnt

  if (rangml == 0) then
    write (6, '("ML: train number of points, min ,max, step..........................:", i6, 3D18.9)') &
      np_train_tmp, minval(dist_maha_train), maxval(dist_maha_train), dint
    write (6, '("ML: The number of MCD train points requested and finally accepted ...:", i6, i6)') np_train_tmp, np_train
  end if
  ! end pack distances for train classes

  ! begin pack distances for ALL classes
  icnt = 0
  do ic = 1, iconf_data
    if (.not. (allocated(config_desc(ic)%energy))) cycle
    do ia = 1, config_real(ic)%nat
      icnt = icnt + 1
    end do
  end do
  size_full = icnt
  if (allocated(dist_mcd_full)) deallocate (dist_mcd_full); allocate (dist_mcd_full(size_full))
  if (allocated(dist_maha_full)) deallocate (dist_maha_full); allocate (dist_maha_full(size_full))
  if (allocated(ic_vector_full)) deallocate (ic_vector_full); allocate (ic_vector_full(size_full))
  if (allocated(ia_vector_full)) deallocate (ia_vector_full); allocate (ia_vector_full(size_full))

  icnt = 0
  do ic = 1, iconf_data
    if (.not. (allocated(config_desc(ic)%energy))) cycle
    do ia = 1, config_real(ic)%nat
      icnt = icnt + 1
      dist_mcd_full(icnt) = config_desc(ic)%stat_dist_mcd(ia)
      dist_maha_full(icnt) = config_desc(ic)%stat_dist_maha(ia)**power_mcd
      ic_vector_full(icnt) = ic
      ia_vector_full(icnt) = ia
    end do
  end do


  if (allocated(index_mcd_full)) deallocate (index_mcd_full); allocate (index_mcd_full(size_full))
  if (allocated(index_maha_full)) deallocate (index_maha_full); allocate (index_maha_full(size_full))
  call indexx(icnt, dist_mcd_full, index_mcd_full)
  call indexx(icnt, dist_maha_full, index_maha_full)




  if (allocated(lhigh)) deallocate (lhigh); allocate (lhigh(size_full))
  if (allocated(llow)) deallocate (llow); allocate (llow(size_full))
  if (allocated(lvec)) deallocate (lvec); allocate (lvec(size_full))

  ! if (allocated(info_kernel)) deallocate(info_kernel) ; allocate(info_kernel(np_all_tmp+np_train_tmp))

  ddmin = minval(dist_maha_full)
  dint = (maxval(dist_maha_full) - ddmin)/dble(np_full_tmp)
  ! dint = (maxval(dist_maha_full) - ddmin)/dble(np_full_tmp - 1)

  icnt = 0
  do ii = 0, np_full_tmp - 1
    lhigh(:) = (dist_maha_full >= (dble(ii + 1)*dint + ddmin))
    llow(:) = (dist_maha_full <= (dble(ii)*dint + ddmin))
    lvec(:) = (llow(:) .eqv. lhigh(:))
    if (count(lvec .eqv. .true., dim=1) > 0) then
      icnt = icnt + 1
      indx = findloc(lvec, .true., dim=1)
      info_kernel(np_train + icnt)%iconf = ic_vector_full(indx)
      info_kernel(np_train + icnt)%ia = ia_vector_full(indx)
    end if
  end do

  np_full = icnt

  if (rangml == 0) then
    write (6, '("ML: full number of points, min ,max, step..........................:", i6, 3D18.9)') &
      np_full_tmp, minval(dist_maha_full), maxval(dist_maha_full), dint
    write (6, '("ML: The number of MCD full points requested and finally accepted ...:", i6, i6)') np_full_tmp, np_full
  end if

  if (allocated(info_kernel_copy)) deallocate (info_kernel_copy)
  allocate (info_kernel_copy(np_full_tmp + np_train_tmp))
  info_kernel_copy(:) = info_kernel(:)
  if (allocated(info_kernel)) deallocate (info_kernel) 
  allocate (info_kernel(np_full + np_train))
  info_kernel(1:np_train + np_full) = info_kernel_copy(1:np_train + np_full)
  deallocate (info_kernel_copy)

  _MLD_END_
end subroutine get_index_kernel_from_dist


subroutine get_index_kernel_from_dist_normalized()

  use ml_in_ndm_module, only: classes_train_for_sigma, rangml, seed
  use module_db_setup, only: iconf_data
  use derived_types, only: config_desc, config_real
  use module_kernel, only: info_kernel, info_kernel_copy, np_kernel_ref, np_kernel_full, power_mcd
  use mld_logger

  implicit none

  integer  :: ic, ia, icnt, iclass, ii, indx, id, ibb
  real(kind(1.d0)), dimension(:), allocatable  :: dist_mcd_full, dist_maha_full
  integer, dimension(:), allocatable     :: index_mcd_full, index_maha_full, ic_vector_full, ia_vector_full

  real(kind(1.d0)), dimension(:), allocatable  :: dist_mcd_train, dist_maha_train
  integer, dimension(:), allocatable     :: index_mcd_train, index_maha_train, ic_vector_train, ia_vector_train
  integer  :: np_train, np_full,  np_train_tmp, np_full_tmp, size_full, size_train, np_rest, np_dist, np_dist_selected
  real(kind(1.d0))     :: dint, ddmin
  logical, dimension(:), allocatable     :: lhigh, llow, lvec
  integer , dimension(:), allocatable :: anum, indexlist 

_NAMECURRENT_("get_index_kernel_from_dist_normalized")


  _MLD_BEGIN_
  ! begin pack distances for train classes
  np_train_tmp = np_kernel_ref
  np_full_tmp = np_kernel_full
  if (allocated(info_kernel)) deallocate (info_kernel)
  allocate (info_kernel(np_full_tmp + np_train_tmp))


  icnt = 0
  do iclass = 1, size(classes_train_for_sigma)
    do ic = 1, iconf_data
      if ((classes_train_for_sigma(iclass) == config_real(ic)%class)) then
        if (.not. (allocated(config_desc(ic)%energy))) cycle
        do ia = 1, config_real(ic)%nat
          icnt = icnt + 1
        end do
      end if
    end do
  end do

  size_train = icnt

  if (allocated(dist_mcd_train)) deallocate (dist_mcd_train)
  allocate (dist_mcd_train(size_train))
  if (allocated(dist_maha_train)) deallocate (dist_maha_train)
  allocate (dist_maha_train(size_train))
  if (allocated(ic_vector_train)) deallocate (ic_vector_train)
  allocate (ic_vector_train(size_train))
  if (allocated(ia_vector_train)) deallocate (ia_vector_train)
  allocate (ia_vector_train(size_train))

  icnt = 0
  do iclass = 1, size(classes_train_for_sigma)
    do ic = 1, iconf_data
      if ((classes_train_for_sigma(iclass) == config_real(ic)%class)) then
        if (.not. (allocated(config_desc(ic)%energy))) cycle
        do ia = 1, config_real(ic)%nat
          icnt = icnt + 1
          dist_mcd_train(icnt) = config_desc(ic)%stat_dist_mcd(ia)
          dist_maha_train(icnt) = config_desc(ic)%stat_dist_maha(ia)**power_mcd
          ic_vector_train(icnt) = ic
          ia_vector_train(icnt) = ia
        end do
      end if
    end do
  end do


  if (allocated(index_mcd_train)) deallocate (index_mcd_train); allocate (index_mcd_train(size_train))
  if (allocated(index_maha_train)) deallocate (index_maha_train); allocate (index_maha_train(size_train))

  call indexx(icnt, dist_mcd_train, index_mcd_train)
  call indexx(icnt, dist_maha_train, index_maha_train)

  if (allocated(lhigh)) deallocate (lhigh); allocate (lhigh(size_train))
  if (allocated(llow)) deallocate (llow); allocate (llow(size_train))
  if (allocated(lvec)) deallocate (lvec); allocate (lvec(size_train))

  ddmin = minval(dist_maha_train)
  dint = (maxval(dist_maha_train) - ddmin)/dble(np_train_tmp)

  icnt = 0
  do ii = 0, np_train_tmp - 1
    lhigh(:) = (dist_maha_train >= (dble(ii + 1)*dint + ddmin))
    llow(:) = (dist_maha_train <= (dble(ii)*dint + ddmin))
    lvec(:) = (llow(:) .eqv. lhigh(:))
    ! write (677,*) ii, count(lvec .eqv. .true., dim=1),(dble(ii)*dint + ddmin), (dble(ii+1)*dint + ddmin)
    if (count(lvec .eqv. .true., dim=1) > 0) then
      icnt = icnt + 1
    end if
  end do  

  np_train = icnt
  np_rest = np_train_tmp - np_train
  icnt = 0
  do ii = 0, np_train_tmp - 1
    lhigh(:) = (dist_maha_train >= (dble(ii + 1)*dint + ddmin))
    llow(:) = (dist_maha_train <= (dble(ii)*dint + ddmin))
    lvec(:) = (llow(:) .eqv. lhigh(:))
    ! write (677,*) ii, count(lvec .eqv. .true., dim=1),(dble(ii)*dint + ddmin), (dble(ii+1)*dint + ddmin)
    np_dist = count(lvec .eqv. .true., dim=1)
    if (np_dist > 0) then
      np_dist_selected = 1 + nint(float(np_rest) * float(np_dist) / float(size_train))
      if (allocated(anum)) deallocate(anum) ; allocate(anum(np_dist_selected))
      if (allocated(indexlist)) deallocate(indexlist) ; allocate(indexlist(np_dist))
      call rks2(np_dist, np_dist_selected, seed, anum)

      do id = 1, np_dist
        indx = findloc(lvec, .true., dim=1)
        indexlist(id) = indx
        lvec(indx)=.false.
      end do

      do id = 1, np_dist_selected
        icnt = icnt + 1
        ibb = anum(id)
        indx = indexlist(ibb)
        info_kernel(icnt)%iconf = ic_vector_train(indx)
        info_kernel(icnt)%ia = ia_vector_train(indx)
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! do id = 1, np_dist_selected
      !   icnt = icnt + 1
      !   tmp_log=.false.
      !   do ibb=1, np_dist_selected
      !     indx = findloc(lvec, .true., dim=1)
      !     lvec(indx)=.false.
      !     if (ibb == anum(id)) then 
      !       lvec(indx)=.true.
      !       tmp_log = .true.
      !       cycle 
      !     end if 
      !   end do 
      !   if (tmp_log) then 
      !     info_kernel(icnt)%iconf = ic_vector_train(indx)
      !     info_kernel(icnt)%ia = ia_vector_train(indx)
      !   end if 
      ! end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      ! indx = findloc(lvec, .true., dim=1)
      ! info_kernel(icnt)%iconf = ic_vector_train(indx)
      ! info_kernel(icnt)%ia = ia_vector_train(indx)
      ! write (677,*) indx, ic_vector_train(indx), ia_vector_train(indx)
    end if
  end do

  np_train = icnt

  if (rangml == 0) then
    write (6, '("ML: train number of points, min ,max, step..........................:", i6, 3D18.9)') &
      np_train_tmp, minval(dist_maha_train), maxval(dist_maha_train), dint
    write (6, '("ML: The number of MCD train points requested and finally accepted ...:", i6, i6)') np_train_tmp, np_train
  end if
  ! end pack distances for train classes

  ! begin pack distances for ALL classes
  icnt = 0
  do ic = 1, iconf_data
    if (.not. (allocated(config_desc(ic)%energy))) cycle
    do ia = 1, config_real(ic)%nat
      icnt = icnt + 1
    end do
  end do
  size_full = icnt
  if (allocated(dist_mcd_full)) deallocate (dist_mcd_full); allocate (dist_mcd_full(size_full))
  if (allocated(dist_maha_full)) deallocate (dist_maha_full); allocate (dist_maha_full(size_full))
  if (allocated(ic_vector_full)) deallocate (ic_vector_full); allocate (ic_vector_full(size_full))
  if (allocated(ia_vector_full)) deallocate (ia_vector_full); allocate (ia_vector_full(size_full))

  icnt = 0
  do ic = 1, iconf_data
    if (.not. (allocated(config_desc(ic)%energy))) cycle
    do ia = 1, config_real(ic)%nat
      icnt = icnt + 1
      dist_mcd_full(icnt) = config_desc(ic)%stat_dist_mcd(ia)
      dist_maha_full(icnt) = config_desc(ic)%stat_dist_maha(ia)**power_mcd
      ic_vector_full(icnt) = ic
      ia_vector_full(icnt) = ia
    end do
  end do


  if (allocated(index_mcd_full)) deallocate (index_mcd_full); allocate (index_mcd_full(size_full))
  if (allocated(index_maha_full)) deallocate (index_maha_full); allocate (index_maha_full(size_full))
  call indexx(icnt, dist_mcd_full, index_mcd_full)
  call indexx(icnt, dist_maha_full, index_maha_full)




  if (allocated(lhigh)) deallocate (lhigh); allocate (lhigh(size_full))
  if (allocated(llow)) deallocate (llow); allocate (llow(size_full))
  if (allocated(lvec)) deallocate (lvec); allocate (lvec(size_full))

  ! if (allocated(info_kernel)) deallocate(info_kernel) ; allocate(info_kernel(np_all_tmp+np_train_tmp))

  ddmin = minval(dist_maha_full)
  dint = (maxval(dist_maha_full) - ddmin)/dble(np_full_tmp)
  ! dint = (maxval(dist_maha_full) - ddmin)/dble(np_full_tmp - 1)

  icnt = 0
  do ii = 0, np_full_tmp - 1
    lhigh(:) = (dist_maha_full >= (dble(ii + 1)*dint + ddmin))
    llow(:) = (dist_maha_full <= (dble(ii)*dint + ddmin))
    lvec(:) = (llow(:) .eqv. lhigh(:))
    if (count(lvec .eqv. .true., dim=1) > 0) then
      icnt = icnt + 1
    end if
  end do

  np_full = icnt
  np_rest = np_full_tmp - np_full
  icnt = 0
  do ii = 0, np_full_tmp - 1
    lhigh(:) = (dist_maha_full >= (dble(ii + 1)*dint + ddmin))
    llow(:) = (dist_maha_full <= (dble(ii)*dint + ddmin))
    lvec(:) = (llow(:) .eqv. lhigh(:))
    np_dist = count(lvec .eqv. .true., dim=1)
    if (np_dist > 0) then
      np_dist_selected = 1 + nint(float(np_rest) * float(np_dist) / float(size_full))
      if (allocated(anum)) deallocate(anum) ; allocate(anum(np_dist_selected))
      if (allocated(indexlist)) deallocate(indexlist) ; allocate(indexlist(np_dist))
      call rks2(np_dist, np_dist_selected, seed, anum)

      do id = 1, np_dist
        indx = findloc(lvec, .true., dim=1)
        indexlist(id) = indx
        lvec(indx)=.false.
      end do

      do id = 1, np_dist_selected
        icnt = icnt + 1
        ibb = anum(id)
        indx = indexlist(ibb)
        info_kernel(np_train + icnt)%iconf = ic_vector_full(indx)
        info_kernel(np_train + icnt)%ia = ia_vector_full(indx)
      end do

      ! do id = 1, np_dist_selected
      !   icnt = icnt + 1
      !   indx = findloc(lvec, .true., dim=1)
      !   info_kernel(np_train + icnt)%iconf = ic_vector_full(indx)
      !   info_kernel(np_train + icnt)%ia = ia_vector_full(indx)
      !   lvec(indx)=.false.
      ! end do
    end if
  end do

  np_full = icnt

  if (rangml == 0) then
    write (6, '("ML: full number of points, min ,max, step..........................:", i6, 3D18.9)') &
      np_full_tmp, minval(dist_maha_full), maxval(dist_maha_full), dint
    write (6, '("ML: The number of MCD full points requested and finally accepted ...:", i6, i6)') np_full_tmp, np_full
  end if

  if (allocated(info_kernel_copy)) deallocate (info_kernel_copy); allocate (info_kernel_copy(np_full_tmp + np_train_tmp))
  info_kernel_copy(:) = info_kernel(:)
  if (allocated(info_kernel)) deallocate (info_kernel); allocate (info_kernel(np_full + np_train))
  info_kernel(1:np_train + np_full) = info_kernel_copy(1:np_train + np_full)
  deallocate (info_kernel_copy)
  _MLD_END_
end subroutine get_index_kernel_from_dist_normalized



subroutine new_distances(eigenvalues, eigenvectors, mean, size_eigen, inamefile)

  use ml_in_ndm_module, only: rangml
  use module_db_setup, only: iconf_data, db_path
  use derived_types, only: config_real, config_desc
  use module_db_poscar, only: iread_energy
  use mld_unit
  use mld_logger

  implicit none

  integer, intent(in)  :: size_eigen, inamefile
  real(kind=kind(1.d0)), dimension(size_eigen, size_eigen) :: eigenvectors
  real(kind=kind(1.d0)), dimension(size_eigen) :: eigenvalues, mean

  integer  :: ic, ia, icount, jj, ifenergy, iftarget 
  real(kind=kind(1.d0))      :: tmp0, tmp1, tmp2, tmp22, tmp3, tmp4, ds(size_eigen), prod

  _NAMECURRENT_("new_distances")


  _MLD_BEGIN_
  if (rangml ==0) then 
    write (inamefile, '("# l^0        l^-1         l^-2         l^-3        l^1        l^2      local_energy")')
  end if 
  icount = 0
  do ic = 1, iconf_data

    if (rangml ==0) then 
      open (newunit=ifenergy, file="desc"//trim(adjustl(db_path))//"/"//config_real(ic)%filename(4:16)//".energy", status="unknown", form="formatted")
      open (newunit=iftarget, file="desc"//trim(adjustl(db_path))//"/"//config_real(ic)%filename(4:16)//".target", status="unknown", form="formatted")
    end if

    do ia = 1, config_real(ic)%nat
      ds(:) = (config_desc(ic)%energy(:, ia) - mean(:))
      tmp0 = 0.d0
      tmp1 = 0.d0
      tmp2 = 0.d0
      tmp22 = 0.d0
      tmp3 = 0.d0
      tmp4 = 0.d0
      do jj = 1, size(eigenvalues, 1)
        prod = (dot_product(ds(:), eigenvectors(:, jj)))**2
        tmp0 = tmp0 + prod
        tmp1 = tmp1 + prod/eigenvalues(jj)
        tmp2 = tmp2 + prod/eigenvalues(jj)**2
        tmp22 = tmp22 + prod/eigenvalues(jj)**3
        tmp3 = tmp3 + prod*eigenvalues(jj)
        tmp4 = tmp3 + prod*eigenvalues(jj)**2
      end do
      if (rangml ==0 ) then  
        write (inamefile, '(7e20.10)') tmp0, tmp1, tmp2, tmp22, tmp3, tmp4, config_desc(ic)%stat_energy(ia)
        write (ifenergy, '(e20.10)') config_desc(ic)%stat_energy(ia)
        write (iftarget, '(i6, 2e20.10,"   ", a2,"   ", a3,"   ", a)') ia, config_desc(ic)%stat_energy(ia), config_real(ic)%energy(iread_energy)/config_real(ic)%nat, config_real(ic)%class, config_real(ic)%klm, config_real(ic)%cnumber
      end if 
    end do
    if (rangml ==0 ) then 
      close (unit=ifenergy, status="keep")
      close (unit=iftarget, status="keep")
    end if 
  end do
  _MLD_END_
end subroutine new_distances





subroutine get_pca_projection(n_pca, mean, size_mean, eigmat, size_eigmat, inamefile)

  use ml_in_ndm_module, only: rangml
  use module_db_setup, only: iconf_data
  use derived_types, only: config_desc, config_real
  use mld_logger

  implicit none

  integer, intent(in)  :: n_pca, size_mean, size_eigmat, inamefile
  real(kind=kind(1.d0)), dimension(size_mean), intent(in)  :: mean
  real(kind=kind(1.d0)), dimension(size_eigmat, size_eigmat), intent(in)     :: eigmat
  real(kind=kind(1.d0)), dimension(n_pca)      :: tmp_proj
  integer  :: ip, ic, ia
  character(len=100)   :: chfmt

  _NAMECURRENT_("get_pca_projection")


  _MLD_BEGIN_
  if (size_mean /= size_eigmat) then
    if (rangml == 0) write (6, *) 'The sizes of mean and the sample covariance matrix are different: ', size_mean, size_eigmat
    if (rangml == 0) write (6, *) "Fatal. stop in subroutine"//NAMECURRENT
    stop 'stop in '//NAMECURRENT
  end if

  write (chfmt, *) '(', int(n_pca), 'e20.10,1x,i7,1x,e20.10,1x,a,1x,a)'
  do ic = 1, iconf_data
    if (.not. (allocated(config_desc(ic)%energy))) cycle
    do ia = 1, config_real(ic)%nat
      do ip = 1, n_pca
        tmp_proj(ip) = DOT_PRODUCT((config_desc(ic)%energy(:, ia) - mean(:)), eigmat(:, size_eigmat - ip + 1))
      end do
      write (inamefile, chfmt) tmp_proj(:), ia, config_desc(ic)%stat_dist_maha(ia), config_real(ic)%class, config_real(ic)%filename
    end do
  end do
  _MLD_END_
end subroutine get_pca_projection



subroutine get_pca_projection_kernel(n_pca, mean, size_mean, eigmat, size_eigmat, inamefile)
  use ml_in_ndm_module, only: rangml
  use derived_types, only: config_desc, config_real
  use module_kernel, only: info_kernel
  use mld_logger

  implicit none

  integer, intent(in)  :: n_pca, size_mean, size_eigmat, inamefile
  real(kind=kind(1.d0)), dimension(size_mean), intent(in)  :: mean
  real(kind=kind(1.d0)), dimension(size_eigmat, size_eigmat), intent(in)     :: eigmat
  real(kind=kind(1.d0)), dimension(n_pca)      :: tmp_proj
  integer  :: ip, ic, ia, ii
  character(len=100)   :: chfmt

  _NAMECURRENT_("get_pca_projection_kernel")


  _MLD_BEGIN_
  if (size_mean /= size_eigmat) then
    if (rangml == 0) write (6, *) 'The sizes of mean and the sample covariance matrix are different: ', size_mean, size_eigmat
    if (rangml == 0) write (6, *) "Fatal. stop in subroutine"//NAMECURRENT
    stop 'stop in '//NAMECURRENT
  end if

  write (chfmt, *) '(', int(n_pca), 'e20.10,1x,i7,1x,e20.10,1x,a,1x,a)'
  do ii = 1, size(info_kernel)
    ic = info_kernel(ii)%iconf
    ia = info_kernel(ii)%ia
    ! do ic = 1, iconf_data
    !   do ia=1, config_real(ic)%nat
    do ip = 1, n_pca
      tmp_proj(ip) = DOT_PRODUCT((config_desc(ic)%energy(:, ia) - mean(:)), eigmat(:, size_eigmat - ip + 1))
    end do
    if (rangml == 0 ) then 
      write (inamefile, chfmt) tmp_proj(:), ia, config_desc(ic)%stat_dist_maha(ia), config_real(ic)%class, config_real(ic)%filename
    end if 
    !   end do
    ! end do
  end do
  _MLD_END_
end subroutine get_pca_projection_kernel

subroutine write_mcd_distances_for_some_classes(classes_for_sigma, no_of_classes_for_sigma, inamefile)
  use ml_in_ndm_module, only: rangml
  use module_db_setup, only: iconf_data
  use derived_types, only: config_desc, config_real
  use module_Sigma_matrix, only: first_passage_mcd
  use mld_logger

  implicit none

  integer, intent(in)  :: no_of_classes_for_sigma, inamefile
  character(len=2), dimension(no_of_classes_for_sigma), intent(in)     :: classes_for_sigma
  integer  :: ic, ia, iclass

  _NAMECURRENT_("write_mcd_distances_for_some_classes")


  _MLD_BEGIN_
  ! writing only on one procs
  if (rangml == 0) then
    do iclass = 1, size(classes_for_sigma)
      do ic = 1, iconf_data
        if ((classes_for_sigma(iclass) == config_real(ic)%class)) then
          if (.not. (allocated(config_desc(ic)%energy))) cycle
          do ia = 1, config_real(ic)%nat
            if (first_passage_mcd) then
              write (inamefile, '(4e20.10, i7," ",a," ",a, " ",i8)') &
                config_desc(ic)%stat_norm(ia), config_desc(ic)%stat_norm_mean(ia), config_desc(ic)%stat_dist_mcd(ia), config_desc(ic)%stat_energy(ia), ia, config_real(ic)%class, trim(config_real(ic)%filename), ic
            else
              write (inamefile, '(4e20.10, i7," ",a," ",a, " ",i8)') &
                config_desc(ic)%stat_norm(ia), config_desc(ic)%stat_norm_mean(ia), config_desc(ic)%stat_dist_maha(ia), config_desc(ic)%stat_energy(ia), ia, config_real(ic)%class, trim(config_real(ic)%filename), ic
            end if
          end do
        end if
      end do
    end do
    close (inamefile, status='keep')
  end if                  ! rangml
  _MLD_END_
end subroutine write_mcd_distances_for_some_classes



subroutine write_mcd_distances(inamefile)
  use ml_in_ndm_module, only: rangml
  use module_db_setup, only: iconf_data
  use derived_types, only: config_desc, config_real
  use module_Sigma_matrix, only: first_passage_mcd
  use mld_logger

  implicit none

  integer, intent(in)  :: inamefile
  integer  :: ic, ia

  _NAMECURRENT_("write_mcd_distances")

  _MLD_BEGIN_
  ! writing only on one procs
  if (rangml == 0) then
    do ic = 1, iconf_data
      if (.not. (allocated(config_desc(ic)%energy))) cycle
      do ia = 1, config_real(ic)%nat
        if (first_passage_mcd) then
          write (inamefile, '(4e20.10, i7," ",a," ",a, " ",i8)') &
            config_desc(ic)%stat_norm(ia), config_desc(ic)%stat_norm_mean(ia), config_desc(ic)%stat_dist_mcd(ia), config_desc(ic)%stat_energy(ia), ia, config_real(ic)%class, trim(config_real(ic)%filename), ic
        else
          write (inamefile, '(4e20.10, i7," ",a," ",a, " ",i8)') &
            config_desc(ic)%stat_norm(ia), config_desc(ic)%stat_norm_mean(ia), config_desc(ic)%stat_dist_maha(ia), config_desc(ic)%stat_energy(ia), ia, config_real(ic)%class, trim(config_real(ic)%filename), ic
        end if
      end do
    end do
    close (inamefile, status='keep')
  end if                  ! rangml
  _MLD_END_
end subroutine write_mcd_distances



subroutine compute_mcd_dist_for_some_classes(classes_for_sigma, no_of_classes_for_sigma, Sigma_inv, dim_01, mean)
  use module_db_setup, only: iconf_data
  use derived_types, only: config_desc, config_real
  use mld_logger

  implicit none

  integer, intent(in)  :: dim_01, no_of_classes_for_sigma
  character(len=2), dimension(no_of_classes_for_sigma), intent(in)     :: classes_for_sigma
  real(kind=kind(1.d0)), dimension(dim_01, dim_01), intent(in)   :: Sigma_inv
  real(kind=kind(1.d0)), dimension(dim_01), intent(in)     :: mean
  integer  :: iclass, ic

  _NAMECURRENT_("compute_mcd_dist_for_some_class")


  _MLD_BEGIN_
  do iclass = 1, size(classes_for_sigma)
    do ic = 1, iconf_data
      if ((classes_for_sigma(iclass) == config_real(ic)%class)) then
        if (.not. (allocated(config_desc(ic)%energy))) cycle
        call compute_mcd_dist_for_one_config(ic, Sigma_inv, dim_01, mean)
        call compute_energy_per_atom_for_one_config(ic, Sigma_inv, dim_01)
      end if
    end do
  end do
  _MLD_END_
end subroutine compute_mcd_dist_for_some_classes



subroutine compute_mcd_dist(Sigma_inv, dim_01, mean)
  use module_db_setup, only: iconf_data
  use derived_types, only: config_desc
  use mld_logger

  implicit none

  integer, intent(in)  :: dim_01
  real(kind=kind(1.d0)), dimension(dim_01, dim_01), intent(in)   :: Sigma_inv
  real(kind=kind(1.d0)), dimension(dim_01), intent(in)     :: mean
  integer  :: ic

  _NAMECURRENT_("compute_mcd_dist")



  _MLD_BEGIN_
  do ic = 1, iconf_data
    if (.not. (allocated(config_desc(ic)%energy))) cycle
    call compute_mcd_dist_for_one_config(ic, Sigma_inv, dim_01, mean)
    call compute_energy_per_atom_for_one_config(ic, Sigma_inv, dim_01)
  end do
  _MLD_END_
end subroutine compute_mcd_dist



subroutine compute_right_vector_with_energy(right_vector, dim_right_vector)
  use module_db_setup, only: iconf_data
  use derived_types, only: config_desc, config_real
  use module_db_poscar, only: iread_energy

  implicit none

  integer, intent(in)  :: dim_right_vector
  real(kind=kind(1.d0)), dimension(dim_right_vector) :: right_vector(dim_right_vector)
  integer  :: ic, ia, icount

  right_vector(:) = 0.d0
  icount = 0
  do ic = 1, iconf_data
    do ia = 1, config_real(ic)%nat
      icount = icount + 1
      right_vector(:) = right_vector(:) + config_real(ic)%energy(iread_energy)*config_desc(ic)%energy(:, ia)
    end do
  end do
  right_vector(:) = right_vector(:)/dble(icount)
end subroutine compute_right_vector_with_energy



subroutine compute_energy_per_atom_for_one_config(iconf, Sigma_inv, dim_01)

#if(PARA)
  !TORC! use mpi
  use mld_mpi, only: comm_mld  
#endif
  use ml_in_ndm_module, only: rangml, i_start_at, i_final_at, mld_order, mld_linear
  use derived_types, only: config_desc, config_real
  use set_limits, only: set_limit_for_atoms
  use module_Sigma_matrix, only: right_vector_with_energy, analysis_local_energy
  use snap, only: w_params
  use temporary_data_cov, only: dim_xdesc
  use mld_logger
  !TORC! use my_mpi_subroutines, only: my_barrier_mld

  implicit none

  integer, intent(in)  :: iconf, dim_01
  real(kind=kind(1.d0)), dimension(dim_01, dim_01), intent(in)   :: Sigma_inv
  real(kind=kind(1.d0))      :: tmp
  integer  :: ia, d2

  _NAMECURRENT_("compute_energy_per_atom_for_one_config")



  _MLD_BEGIN_

  if (allocated(config_desc(iconf)%stat_energy)) deallocate (config_desc(iconf)%stat_energy); allocate (config_desc(iconf)%stat_energy(config_real(iconf)%nat))
  config_desc(iconf)%stat_energy(:) = 0.d0
  call set_limit_for_atoms(rangml, config_real(iconf)%nat, i_start_at, i_final_at)

  if ((i_start_at == 0) .and. (i_final_at == 0)) then
    config_desc(iconf)%stat_energy(:) = 0.d0
  else
    if (analysis_local_energy) then 
      do ia = i_start_at, i_final_at
        tmp = 0.d0
        do d2 = 1, dim_01
          tmp = tmp + dot_product(config_desc(iconf)%energy(:, ia), Sigma_inv(:, d2))*right_vector_with_energy(d2)
        end do
        config_desc(iconf)%stat_energy(ia) = tmp
        if (mld_order == mld_linear) then
          config_desc(iconf)%stat_energy(ia) = dble(1.d0)*w_params(1, 1) + dot_product(config_desc(iconf)%energy(1:dim_xdesc, ia), w_params(2:dim_xdesc + 1, 1))
        else
          if (rangml == 0) write (6, '("ML: THe analysis not yet implemented for other types than linear, mld_order is ... ...:", i6)') mld_order
          stop 'compute_energy_per_atom_for_one_config with mld_order'
        end if
      end do  !ia 
    end if 
  end if

#if(PARA)
  !TORC! call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(iconf)%stat_energy, config_real(iconf)%nat, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, codeml)
  call comm_mld%sum(config_desc(iconf)%stat_energy)
  
  !TORC! call my_barrier_mld(codeml)
  call comm_mld%barrier()
#endif

  _MLD_END_
end subroutine compute_energy_per_atom_for_one_config



subroutine compute_mcd_dist_for_one_config(iconf, Sigma_inv, dim_01, mean)

#if(PARA)
  !TORC! use mpi
  use mld_mpi, only: comm_mld
#endif
  use ml_in_ndm_module, only: rangml, i_start_at, i_final_at
  use derived_types, only: config_desc, config_real
  use set_limits, only: set_limit_for_atoms
  use module_Sigma_matrix, only: first_passage_mcd
  use mld_logger
  !TORC! use my_mpi_subroutines, only: my_barrier_mld

  implicit none

  integer, intent(in)  :: iconf, dim_01
  real(kind=kind(1.d0)), dimension(dim_01, dim_01), intent(in)   :: Sigma_inv
  real(kind=kind(1.d0)), dimension(dim_01), intent(in)     :: mean
  real(kind=kind(1.d0))      :: tmp
  integer  :: ia, d2

  _NAMECURRENT_("compute_mcd_dist")

  _MLD_BEGIN_

  if (first_passage_mcd) then
    if (allocated(config_desc(iconf)%stat_dist_mcd)) deallocate (config_desc(iconf)%stat_dist_mcd); allocate (config_desc(iconf)%stat_dist_mcd(config_real(iconf)%nat))
  else
    if (allocated(config_desc(iconf)%stat_dist_maha)) deallocate (config_desc(iconf)%stat_dist_maha); allocate (config_desc(iconf)%stat_dist_maha(config_real(iconf)%nat))
  end if
  if (allocated(config_desc(iconf)%stat_norm)) deallocate (config_desc(iconf)%stat_norm); allocate (config_desc(iconf)%stat_norm(config_real(iconf)%nat))
  if (allocated(config_desc(iconf)%stat_norm_mean)) deallocate (config_desc(iconf)%stat_norm_mean); allocate (config_desc(iconf)%stat_norm_mean(config_real(iconf)%nat))

  config_desc(iconf)%stat_norm(:) = 0.d0
  config_desc(iconf)%stat_norm_mean(:) = 0.d0
  if (first_passage_mcd) then
    config_desc(iconf)%stat_dist_mcd(:) = 0.d0
  else
    config_desc(iconf)%stat_dist_maha(:) = 0.d0
  end if

  call set_limit_for_atoms(rangml, config_real(iconf)%nat, i_start_at, i_final_at)

  if ((i_start_at == 0) .and. (i_final_at == 0)) then
    config_desc(iconf)%stat_norm(:) = 0.d0
    config_desc(iconf)%stat_norm_mean(:) = 0.d0
    if (first_passage_mcd) then
      config_desc(iconf)%stat_dist_mcd(:) = 0.d0
    else
      config_desc(iconf)%stat_dist_maha(:) = 0.d0
    end if

  else
    ! do ia=1,config_real(iconf)%nat
    do ia = i_start_at, i_final_at
      config_desc(iconf)%stat_norm(ia) = dot_product(config_desc(iconf)%energy(:, ia), config_desc(iconf)%energy(:, ia))
      config_desc(iconf)%stat_norm_mean(ia) = dot_product((config_desc(iconf)%energy(:, ia) - mean(:)), (config_desc(iconf)%energy(:, ia) - mean(:)))
      tmp = 0.d0
      do d2 = 1, dim_01
        ! tmp = tmp + config_desc(iconf)%energy(d1,ia)*dot_product(config_desc(iconf)%energy(:,ia),Sigma_inv(d1,:))
        tmp = tmp + dot_product((config_desc(iconf)%energy(:, ia) - mean(:)), Sigma_inv(:, d2))*(config_desc(iconf)%energy(d2, ia) - mean(d2))
      end do
      if (first_passage_mcd) then
        config_desc(iconf)%stat_dist_mcd(ia) = tmp
      else
        config_desc(iconf)%stat_dist_maha(ia) = tmp
      end if
    end do
  end if

  ! TODO if sequential
#if(PARA)
  !TORC! call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(iconf)%stat_norm, config_real(iconf)%nat, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, codeml)
  call comm_mld%sum(config_desc(iconf)%stat_norm)
  !TORC! call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(iconf)%stat_norm_mean, config_real(iconf)%nat, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, codeml)
  call comm_mld%sum(config_desc(iconf)%stat_norm_mean)
  if (first_passage_mcd) then
    !TORC! call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(iconf)%stat_dist_mcd, config_real(iconf)%nat, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, codeml)
    call comm_mld%sum(config_desc(iconf)%stat_dist_mcd)
  else
    !TORC! call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(iconf)%stat_dist_maha, config_real(iconf)%nat, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, codeml)
    call comm_mld%sum(config_desc(iconf)%stat_dist_maha)
  end if

  !TORC! call my_barrier_mld(codeml)
  call comm_mld%barrier()
#endif

  _MLD_END_

end subroutine compute_mcd_dist_for_one_config


 
